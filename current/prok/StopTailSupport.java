package prok;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

import dna.Data;
import dna.GeneticCode;
import prot.DigestSuffix;
import prot.StopTailDetector;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Optional stop-tail-only evidence for CallGenes, backed by a pinned detector resource.
 * Stop-tail representatives are evidence for terminal stops, not functional labels.
 * @author Yelan
 */
public final class StopTailSupport implements AutoCloseable {

	public static StopTailSupport load(final String manifest,final StopTailDetector.Method method,
			final double cutoff,final float middleBonus,final float stopBonus,final String[] outputs) throws IOException{
		if(method==null || !Double.isFinite(cutoff) || !Float.isFinite(middleBonus+stopBonus) ||
				middleBonus<0 || stopBonus<0){
			throw new IllegalArgumentException("Stop-tail support requires a method, finite cutoff, "
				+"and finite nonnegative bonuses");
		}
		final Manifest m=Manifest.load(manifest);
		m.protect(outputs);
		if(method==StopTailDetector.Method.HBM && m.hbm==null){
			throw new IllegalArgumentException("tailmethod=hbm requires HBM paths in the stop-tail manifest");
		}
		final StopTailDetector.HbmInput hbm=method==StopTailDetector.Method.HBM ?
			new StopTailDetector.HbmInput(m.hbm, m.hbmPin, m.provenance, m.provenancePin, m.background, m.backgroundPin) : null;
		final StopTailDetector.PairedSettings paired=new StopTailDetector.PairedSettings(
			m.profiles, m.profilesPin, m.sidecar, m.sidecarPin, m.kmers, m.kmersPin,
			m.buckets, m.mode, m.ddlK, m.seed, m.indexed, m.namespace, m.ddlTop, m.coverTop);
		final StopTailDetector.Reference reference;
		try{reference=StopTailDetector.load(m.consensus, m.consensusPin, paired, hbm);}
		catch(Exception e){throw new IOException("Could not load stop-tail detector manifest: "+manifest, e);}
		m.verifyPin();
		return new StopTailSupport(reference, m.pin, m.consensusPin, method, cutoff, middleBonus, stopBonus);
	}

	StopTailSupport(final StopTailDetector.Reference reference_,final String manifestPin_,final String consensusPin_,
			final StopTailDetector.Method method_,final double cutoff_,final float middleBonus_,final float stopBonus_){
		reference=reference_; method=method_; cutoff=cutoff_; cutoffText=Double.toString(cutoff_);
		middleBonus=middleBonus_; stopBonus=stopBonus_;
		manifestPin=manifestPin_; consensusPin=consensusPin_;
		policy=new StopTailDetector.Policy(0.7, 70, 10, true);
		pool=new Pool(reference, policy, method, cutoff, Shared.threads());
	}

	public Worker newWorker(){return new Worker();}

	@Override
	public void close(){pool.close();}

	public final class Worker {
		/** Scores one oriented longest-ORF batch, preserving input order for later start enumeration. */
		public StopEvidence[] prepareStops(final ArrayList<Orf> longest, final byte[] orientedBases, final GeneticCode code){
			try{
				final GeneticCode effective=code==null ? GeneticCode.forTable(11) : code;
				final StopTailDetector.Query[] queries=new StopTailDetector.Query[longest.size()];
				final StopEvidence[] results=new StopEvidence[longest.size()];
				for(int i=0; i<longest.size(); i++){
					final Orf orf=longest.get(i);
					assert(orf.type==ProkObject.CDS && orf.scaflen==orientedBases.length)
						: "Stop-tail batches must preserve the contig coordinate frame";
					if(orf.stopTruncated){results[i]=new StopEvidence(null, "synthetic_stop", effective.id(), 0); continue;}
					final int first=orf.flipped()==0 ? orf.start : orientedBases.length-1-orf.stop;
					final int last=orf.flipped()==0 ? orf.stop : orientedBases.length-1-orf.start;
					final byte[] peptide=ProteinSupport.translateTail(orientedBases, first, last, effective, TAIL_CAP);
					final String reason=ProteinSupport.translationReason(peptide);
					if(reason!=null){results[i]=new StopEvidence(null, reason, effective.id(), peptide.length);}
					else if(peptide.length<MIN_TAIL){results[i]=new StopEvidence(null, "short_tail", effective.id(), peptide.length);}
					else{queries[i]=new StopTailDetector.Query("stop_"+i, peptide);}
				}
				final StopTailDetector.Result[] hits=pool.assign(queries);
				for(int i=0; i<results.length; i++){
					if(queries[i]!=null){
						results[i]=new StopEvidence(hits[i], null, effective.id(), queries[i].length(), method, cutoff);
						searched++; candidates+=hits[i].candidates;
						if(hits[i].noKmers){noKmers++;}
						if(hits[i].candidates==0){noCandidates++;}
						final StopTailDetector.Work consensus=hits[i].blosumWork, hbm=hits[i].hbmWork;
						consensusScored+=consensus.scored; consensusTraces+=consensus.traces;
						consensusGeometryFail+=consensus.geometryFail; consensusGeometryAccepted+=consensus.accepted;
						hbmScored+=hbm.scored; hbmTraces+=hbm.traces;
						hbmGeometryFail+=hbm.geometryFail; hbmGeometryAccepted+=hbm.accepted;
						if(results[i].supported()){assigned++;}
					}else{unsupported++;}
				}
				return results;
			}catch(RuntimeException | Error failure){
				shared.KillSwitch.kill("Stop-tail batch failed before candidate scoring: "+failure);
				throw failure;
			}
		}
		/** Attaches one immutable stop decision to one admitted candidate start at that stop. */
		public void annotate(final Orf orf, final StopEvidence stop){
			try{
				assert(orf.type==ProkObject.CDS && orf.stopTailEvidence==null && stop!=null)
					: "Each CDS candidate inherits its existing stop evidence exactly once";
				final float nativeScore=orf.orfScore;
				final float orfAddition=stop.supported() ? middleBonus : 0;
				final float stopAddition=stop.supported() ? stopBonus : 0;
				ProteinSupport.addBonuses(orf, stopAddition, orfAddition);
				orf.stopTailEvidence=new Evidence(StopTailSupport.this, stop, stopAddition, orfAddition, nativeScore);
				if(orf.proteinBonus>0){boosted++; if(orf.proteinOnlyCandidate){rescuedStarts++;}}
			}catch(RuntimeException | AssertionError failure){
				shared.KillSwitch.kill("Stop-tail evidence attachment failed for "
					+orf.scafName+":"+orf.start+"-"+orf.stop+": "+failure);
				throw failure;
			}
		}
		public long searched, assigned, unsupported, withoutStarts, boosted, rescuedStarts;
		public long candidates, noKmers, noCandidates;
		public long consensusScored, consensusTraces, consensusGeometryFail, consensusGeometryAccepted;
		public long hbmScored, hbmTraces, hbmGeometryFail, hbmGeometryAccepted;
	}

	public static final class StopEvidence {
		StopEvidence(StopTailDetector.Result result_, String reason_, int table_, int length_){
			this(result_, reason_, table_, length_, null, UNSUPPORTED_CUTOFF);
		}
		StopEvidence(StopTailDetector.Result result_, String reason_, int table_, int length_,
				StopTailDetector.Method method_, double cutoff_){
			if(result_!=null && (method_==null || !Double.isFinite(cutoff_))){
				throw new IllegalArgumentException("Searched stop-tail evidence requires a method and finite cutoff");
			}
			if(result_==null && !Double.isFinite(cutoff_)){
				throw new IllegalArgumentException("Unsupported stop-tail evidence must use a finite ignored cutoff");
			}
			result=result_; reason=reason_; table=table_; length=length_; method=method_; cutoff=cutoff_;
		}
		public boolean supported(){return result!=null && hit().passes(cutoff);}
		public String reason(){return reason;}
		public int length(){return length;}
		public int candidates(){return result==null ? -1 : result.candidates;}
		public boolean noKmers(){return result!=null && result.noKmers;}
		public String representative(){return supported() ? hit().representative : null;}
		StopTailDetector.Hit hit(){return result.hit(method);}
		final StopTailDetector.Result result;
		final String reason;
		final int table, length;
		final StopTailDetector.Method method;
		final double cutoff;
		/** Ignored for explicit unavailable states; kept finite so StopEvidence never stores a NaN cutoff. */
		private static final double UNSUPPORTED_CUTOFF=0;
	}

	static final class Evidence {
		Evidence(StopTailSupport reference_, StopEvidence stop_, float stopAddition_, float orfAddition_, float nativeScore_){
			reference=reference_; stop=stop_; stopAddition=stopAddition_; orfAddition=orfAddition_; nativeScore=nativeScore_;
		}
		void append(final ByteBuilder row){
			final boolean assigned=stop.supported();
			row.append(";bbgap_status=").append(stop.reason!=null ? "unsupported" : assigned ? "stop_supported" : "no_hit");
			row.append(";bbgap_mode=tail;bbgap_reference=").append(reference.manifestPin);
			row.append(";bbgap_evidence_scope=stop_tail;bbgap_method=").append(reference.method.name);
			row.append(";bbgap_cutoff=").append(reference.cutoffText);
			row.append(";bbgap_translation_table=").append(stop.table);
			row.append(";bbgap_tail_length=").append(stop.length);
			row.append(";bbgap_tail_cap=").append(TAIL_CAP);
			row.append(";bbgap_native_score=").append(nativeScore, 6);
			row.append(";bbgap_call_bonus=").append(stopAddition+orfAddition, 6);
			row.append(";bbgap_stop_bonus=").append(stopAddition, 6);
			row.append(";bbgap_orf_bonus=").append(orfAddition, 6);
			if(stop.reason!=null){row.append(";bbgap_reason=").append(stop.reason);}
			if(!assigned){return;}
			final StopTailDetector.Hit hit=stop.hit();
			row.append(";bbgap_score=").append(hit.score, 8);
			row.append(";bbgap_score64=").append(hit.score64);
			row.append(";bbgap_tail_rank=").append(hit.rank);
			row.append(";bbgap_tail_rep_id="); AnnotationGff.escaped(row, hit.representative);
			row.append(";bbgap_tail_columns=").append(hit.columns);
			row.append(";bbgap_tail_paired=").append(hit.paired);
			row.append(";bbgap_tail_terminal_paired=").append(hit.terminalPaired);
			row.append(";bbgap_tail_query_start=").append(hit.queryStart);
			row.append(";bbgap_tail_query_end=").append(hit.queryEnd);
			row.append(";bbgap_tail_model_start=").append(hit.referenceStart);
			row.append(";bbgap_tail_model_end=").append(hit.referenceEnd);
			row.append(";bbgap_tail_candidates=").append(stop.result.candidates);
		}
		boolean protectsFromPruning(){return false;}
		final StopTailSupport reference;
		final StopEvidence stop;
		final float stopAddition, orfAddition, nativeScore;
	}

	static boolean protectsFromPruning(final Orf orf){
		assert(orf!=null) : "Soft pruning operates on existing native candidates";
		return orf.stopTailEvidence!=null && orf.stopTailEvidence.protectsFromPruning();
	}

	private static final class Pool implements AutoCloseable {
		Pool(final StopTailDetector.Reference reference_,final StopTailDetector.Policy policy_,
				final StopTailDetector.Method method_,final double cutoff_,final int threads_){
			reference=reference_; policy=policy_; method=method_; cutoff=cutoff_; threads=Tools.max(1, threads_);
			scratch=ThreadLocal.withInitial(reference::makeScratch);
			final AtomicInteger names=new AtomicInteger();
			pool=threads==1 ? null : Executors.newFixedThreadPool(threads, task->{
				final Thread worker=new Thread(task, "CallGenes-StopTail-"+names.incrementAndGet());
				worker.setDaemon(true); return worker;
			});
		}
		StopTailDetector.Result[] assign(final StopTailDetector.Query[] queries){
			final StopTailDetector.Result[] results=new StopTailDetector.Result[queries.length];
			if(queries.length==0){return results;}
			final int chunks=1+(queries.length-1)/BATCH;
			final AtomicInteger next=new AtomicInteger();
			final Callable<Void> task=()->{
				final StopTailDetector.Scratch local=scratch.get();
				for(int chunk=next.getAndIncrement(); chunk<chunks; chunk=next.getAndIncrement()){
					final int first=chunk*BATCH, end=first+Math.min(BATCH, queries.length-first);
					for(int i=first; i<end; i++){
						if(queries[i]!=null){
							results[i]=reference.assign(queries[i], policy, 1, local);
							results[i].hit(method).passes(cutoff);
						}
					}
				}
				return null;
			};
			try{
				if(pool==null){task.call();}
				else{
					final ArrayList<Callable<Void>> tasks=new ArrayList<Callable<Void>>();
					for(int t=0; t<Math.min(threads, chunks); t++){tasks.add(task);}
					for(Future<Void> result : pool.invokeAll(tasks)){result.get();}
				}
			}catch(InterruptedException e){
				Thread.currentThread().interrupt();
				throw new IllegalStateException("Interrupted while waiting for stop-tail detector workers", e);
			}catch(ExecutionException e){throw new IllegalStateException("A stop-tail detector worker failed", e.getCause());}
			catch(RuntimeException | Error e){throw e;}
			catch(Exception e){throw new IllegalStateException("Stop-tail detector batch failed", e);}
			for(int i=0; i<queries.length; i++){
				assert((queries[i]==null)==(results[i]==null)) :
					"Every searchable stop tail must publish one result at its original index: "+i;
			}
			return results;
		}
		@Override
		public void close(){if(pool!=null){pool.shutdown();}}
		final StopTailDetector.Reference reference;
		final StopTailDetector.Policy policy;
		final StopTailDetector.Method method;
		final double cutoff;
		final int threads;
		final ExecutorService pool;
		final ThreadLocal<StopTailDetector.Scratch> scratch;
	}

	private static final class Manifest {
		static Manifest load(final String file) throws IOException{
			final String resolved=Data.findPath(file, false);
			if(resolved==null){throw new IllegalArgumentException("Missing stop-tail manifest: "+file);}
			final Path manifest=Paths.get(resolved).toRealPath();
			final String pinBefore=DigestSuffix.file(manifest.toString());
			final HashMap<String,String> map=new HashMap<String,String>();
			for(String line : Files.readAllLines(manifest)){
				if(line.isEmpty() || line.charAt(0)=='#'){continue;}
				final String[] split=line.split("\t", -1);
				if(split.length!=2 || map.put(split[0].toLowerCase(Locale.ROOT), split[1])!=null){
					throw new IllegalArgumentException("Malformed stop-tail manifest row: "+line);
				}
			}
			for(String key : map.keySet()){
				if(!KNOWN_KEYS.contains(key)){throw new IllegalArgumentException("Unknown stop-tail manifest key: "+key);}
			}
			final String pinAfter=DigestSuffix.file(manifest.toString());
			if(!pinBefore.equals(pinAfter)){
				throw new IllegalStateException("Stop-tail manifest changed during load: "+manifest);
			}
			final Manifest m=new Manifest();
			m.manifest=manifest; m.root=manifest.getParent(); m.pin=pinBefore;
			m.consensus=m.path(map, "consensus"); m.consensusPin=m.pin(map, "consensussha80");
			m.profiles=m.path(map, "profiles"); m.profilesPin=m.pin(map, "profilessha80");
			m.sidecar=m.path(map, "sidecar"); m.sidecarPin=m.pin(map, "sidecarsha80");
			m.kmers=m.path(map, "kmers"); m.kmersPin=m.pin(map, "kmerssha80");
			m.hbm=m.path(map, "hbm", false); m.hbmPin=m.pin(map, "hbmsha80", false);
			m.provenance=m.path(map, "provenance", false); m.provenancePin=m.pin(map, "provenancesha80", false);
			m.background=m.path(map, "background", false); m.backgroundPin=m.pin(map, "backgroundsha80", false);
			if((m.hbm==null)!=(m.hbmPin==null) || (m.hbm==null)!=(m.provenance==null) || (m.hbm==null)!=(m.provenancePin==null)
					|| (m.hbm==null)!=(m.background==null) || (m.hbm==null)!=(m.backgroundPin==null)){
				throw new IllegalArgumentException("Stop-tail HBM manifest entries must be all present or all absent");
			}
			m.profileMembers=m.profileMembers();
			m.buckets=integer(map, "buckets", 512); m.mode=map.containsKey("mode") ? map.get("mode") : "nonempty";
			m.ddlK=integer(map, "ddlk", 5); m.seed=longValue(map, "seed", 0);
			m.indexed=bool(map, "indexed", true);
			m.namespace=map.containsKey("namespace") ? map.get("namespace") : "tail";
			m.ddlTop=integer(map, "ddltop", 25); m.coverTop=integer(map, "covertop", 25);
			return m;
		}
		void protect(final String[] outputs) throws IOException{
			for(String name : outputs){
				if(!Tools.isOutputFileName(name)){continue;}
				final Path out=new File(name).getCanonicalFile().toPath();
				if(out.toAbsolutePath().normalize().startsWith(root.toAbsolutePath().normalize())){
					throw new IllegalArgumentException("Output would be created under stop-tail resources: "+out+" in "+root);
				}
				if(!Files.exists(out)){continue;}
				for(String input : new String[]{
						manifest.toString(), consensus, profiles, sidecar, kmers, hbm, provenance, background}){
					if(input!=null && Files.exists(Paths.get(input)) && Files.isSameFile(out, Paths.get(input))){
						throw new IllegalArgumentException("Output aliases a stop-tail reference member: "+out+" -> "+input);
					}
				}
				for(String input : profileMembers){
					if(Files.isSameFile(out, Paths.get(input))){
						throw new IllegalArgumentException("Output aliases a stop-tail DDL profile member: "+out+" -> "+input);
					}
				}
			}
		}
		private String[] profileMembers() throws IOException{
			final ArrayList<String> members=new ArrayList<String>();
			final Path profileList=Paths.get(profiles);
			final Path parent=profileList.getParent();
			for(String line : Files.readAllLines(profileList)){
				final String[] fields=line.split("\t", -1);
				if(fields.length!=3){throw new IllegalArgumentException("Malformed DDL profile list row: "+line);}
				Path profile=Paths.get(fields[0]);
				if(!profile.isAbsolute()){profile=parent.resolve(profile);}
				members.add(profile.normalize().toString());
				members.add(profile.getParent().resolve("summary.tsv").normalize().toString());
			}
			return members.toArray(new String[members.size()]);
		}
		void verifyPin(){
			final String current=DigestSuffix.file(manifest.toString());
			if(!pin.equals(current)){
				throw new IllegalStateException("Stop-tail manifest changed during detector load: "+manifest);
			}
		}
		private String path(HashMap<String,String> map,String key){return path(map, key, true);}
		private String path(HashMap<String,String> map,String key,boolean required){
			final String value=map.get(key);
			if(value==null){if(required){throw new IllegalArgumentException("Stop-tail manifest lacks "+key);} return null;}
			final Path p=Paths.get(value);
			return (p.isAbsolute() ? p : root.resolve(p)).normalize().toString();
		}
		private String pin(HashMap<String,String> map,String key){return pin(map, key, true);}
		private String pin(HashMap<String,String> map,String key,boolean required){
			final String value=map.get(key);
			if(value==null){if(required){throw new IllegalArgumentException("Stop-tail manifest lacks "+key);} return null;}
			return value;
		}
		private static int integer(HashMap<String,String> map,String key,int def){
			return map.containsKey(key) ? Integer.parseInt(map.get(key)) : def;
		}
		private static long longValue(HashMap<String,String> map,String key,long def){
			return map.containsKey(key) ? Long.parseLong(map.get(key)) : def;
		}
		private static boolean bool(HashMap<String,String> map,String key,boolean def){
			if(!map.containsKey(key)){return def;}
			final String value=map.get(key);
			if("true".equalsIgnoreCase(value) || "t".equalsIgnoreCase(value)){return true;}
			if("false".equalsIgnoreCase(value) || "f".equalsIgnoreCase(value)){return false;}
			throw new IllegalArgumentException("Stop-tail manifest key "+key+" must be true or false, not "+value);
		}
		private static final HashSet<String> KNOWN_KEYS=new HashSet<String>(Arrays.asList(
			"consensus", "consensussha80", "profiles", "profilessha80", "sidecar", "sidecarsha80",
			"kmers", "kmerssha80", "hbm", "hbmsha80", "provenance", "provenancesha80",
			"background", "backgroundsha80", "buckets", "mode", "ddlk", "seed", "indexed",
			"namespace", "ddltop", "covertop"));
		Path manifest, root;
		String pin, consensus, consensusPin, profiles, profilesPin, sidecar, sidecarPin, kmers, kmersPin;
		String hbm, hbmPin, provenance, provenancePin, background, backgroundPin;
		String mode, namespace;
		String[] profileMembers;
		int buckets, ddlK, ddlTop, coverTop;
		long seed;
		boolean indexed;
	}

	static StopTailDetector.Method methodForName(final String name){
		if(name!=null){for(StopTailDetector.Method method : StopTailDetector.Method.values()){
			if(method.name.equalsIgnoreCase(name)){return method;}
		}}
		throw new IllegalArgumentException("tailmethod must be aai, blosum or hbm");
	}

	private final StopTailDetector.Reference reference;
	private final StopTailDetector.Policy policy;
	private final StopTailDetector.Method method;
	private final Pool pool;
	private final String manifestPin, consensusPin, cutoffText;
	private final double cutoff;
	private final float middleBonus, stopBonus;
	private static final int TAIL_CAP=200, MIN_TAIL=100, BATCH=32;
}
