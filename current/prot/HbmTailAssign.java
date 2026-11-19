package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.PreParser;
import shared.Shared;
import template.ThreadWaiter;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Streams one query shard through frozen tail profiles with existing F4 shortlisting.
 * All shortlisted candidates are scored before the highest-score eligible family wins.
 * No full-protein core gate, production early acceptance or query trimming is introduced.
 * Batch mode can split independently pinned shards across local workers.
 * @author Keqing
 */
public final class HbmTailAssign {

	public static void main(String[] args){
		try{
			final PreParser pp=new PreParser(args, HbmTailAssign.class, false);
			try{runForTest(pp.args);}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	/** Throws on invalid input rather than publishing a success marker for a partial shard. */
	public static void runForTest(String[] args) throws Exception{
		Shared.setThreads(1); Shared.AMINO_IN=true;
		final HashMap<String,String> options=HmmComparisonData.options(args);
		final int requested=Integer.parseInt(options.getOrDefault("shortlist", "50"));
		require(requested>=0, "shortlist must be positive or zero for exhaustive scoring");
		final HbmTailProfileRefiner.TailPolicy policy=parsePolicy(options.getOrDefault("tailpolicy", "construction"));
		if(options.containsKey("shards")){runBatch(options, requested, policy);}
		else{runSingle(options, requested, policy);}
	}

	private static void runSingle(HashMap<String,String> options,int requested,HbmTailProfileRefiner.TailPolicy policy) throws Exception{
		final HashSet<String> allowed=allowed("out", "expected", "shortlist", "tailpolicy");
		allowed.add("in"); allowed.add("insha80");
		for(String key : options.keySet()){require(allowed.contains(key), "Unknown tail-assignment option: "+key);}
		final Shard shard=new Shard("single", Long.parseLong(HmmComparisonData.required(options, "expected")),
			Paths.get(HmmComparisonData.required(options, "in")), options.get("insha80"),
			Paths.get(HmmComparisonData.required(options, "out")));
		require(shard.expected>0, "Expected tail-shard query count must be positive");
		final long setupStart=System.nanoTime();
		verifyReference(options); HbmProfilePilot.requireHash(shard.input.toString(), shard.inputSha80);
		require(!Files.exists(shard.result), "Tail assignment output must be fresh: "+shard.result);
		final Library library=new Library(options, policy);
		final int limit=requested==0 ? library.ids.length : Math.min(requested, library.ids.length);
		stage("setup", "single", 1, shard.expected, -1, library.ids.length, limit,
			"initial_reference_query_verification_and_library_construction", setupStart);
		assignOne(shard, null, true, options, requested, limit, library);
	}

	private static void runBatch(HashMap<String,String> options,int requested,HbmTailProfileRefiner.TailPolicy policy) throws Exception{
		final HashSet<String> allowed=allowed("out", "shards", "shardssha80", "shortlist", "t", "threads", "tailpolicy");
		for(String key : options.keySet()){require(allowed.contains(key), "Unknown tail-assignment option: "+key);}
		final long setupStart=System.nanoTime();
		verifyReference(options);
		final Path shardPath=Paths.get(HmmComparisonData.required(options, "shards"));
		HbmProfilePilot.requireHash(shardPath.toString(), HmmComparisonData.required(options, "shardssha80"));
		final ArrayList<Shard> shards=loadShards(shardPath);
		long queries=0;
		for(Shard shard : shards){
			HbmProfilePilot.requireHash(shard.input.toString(), shard.inputSha80);
			queries=Math.addExact(queries, shard.expected);
		}
		final Path out=Paths.get(HmmComparisonData.required(options, "out"));
		require(!Files.exists(out), "Tail-assignment batch output must be fresh: "+out);
		checkBatchOutput(out, shards);
		final int workers=workers(options, shards.size());
		if(workers>1){checkParallelOutputs(shards);}
		final Library library=new Library(options, policy);
		final int limit=requested==0 ? library.ids.length : Math.min(requested, library.ids.length);
		stage("setup", "batch", shards.size(), queries, -1, library.ids.length, limit,
			"initial_reference_query_verification_and_library_construction", setupStart);
		Files.createDirectory(out);
		if(workers<2){
			for(Shard shard : shards){assignOne(shard, shard.genome+"|", false, options, requested, limit, library);}
		}else{assignParallel(shards, options, requested, limit, library, workers);}
		verifyReference(options);
		HbmProfilePilot.requireHash(shardPath.toString(), options.get("shardssha80"));
		for(Shard shard : shards){HbmProfilePilot.requireHash(shard.input.toString(), shard.inputSha80);}
		final ByteStreamWriter summary=HmmComparisonData.writer(out.resolve("summary.tsv").toString());
		summary.println("shards\tqueries"); summary.println(new ByteBuilder().append(shards.size()).tab().append(queries));
		HmmComparisonData.close(summary);
		final ByteStreamWriter pass=HmmComparisonData.writer(out.resolve("PASS").toString());
		pass.println("HBM_TAIL_ASSIGN_BATCH_PASS"); HmmComparisonData.close(pass);
	}

	private static void assignParallel(ArrayList<Shard> shards,HashMap<String,String> options,int requested,
			int limit,Library library,int workerCount) throws Exception{
		final AtomicInteger next=new AtomicInteger();
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>(workerCount);
		final ArrayList<ProcessThread> started=new ArrayList<ProcessThread>(workerCount);
		for(int i=0; i<workerCount; i++){
			workers.add(new ProcessThread(shards, options, requested, limit, library, next, failure));
		}
		try{
			for(ProcessThread pt : workers){
				pt.start();
				started.add(pt);
			}
		}catch(Throwable e){failure.compareAndSet(null, e);}
		finally{ThreadWaiter.waitForThreadsToFinish(started);}
		if(failure.get()!=null){throw new RuntimeException("Tail assignment worker failed", failure.get());}
	}

	private static void assignOne(Shard shard,String prefix,boolean verifySharedBeforePublish,
			HashMap<String,String> options,int requested,int limit,Library library) throws Exception{
		final ProteinSearcher.ShortlistScratch shortlist=new ProteinSearcher.ShortlistScratch(library.sidecar.dims, library.ids.length, limit);
		final HbmTailCompetitiveGate.Scratch selection=new HbmTailCompetitiveGate.Scratch(library.ids.length, limit);
		final long[] familyCounts=new long[library.ids.length];
		final HashSet<String> seen=prefix==null ? null : new HashSet<String>();
		final long searchStart=System.nanoTime();
		Files.createDirectory(shard.result);
		final ByteStreamWriter calls=HmmComparisonData.writer(shard.result.resolve("assignments.tsv").toString());
		calls.println(library.policy==HbmTailProfileRefiner.TailPolicy.FINAL_BIDIRECTIONAL ? FINAL_HEADER : HEADER);
		final ByteBuilder row=new ByteBuilder();
		long count=0, assigned=0, residues=0, scored=0, recorded=0, rejected=0, nonpositive=0;
		final Streamer input=StreamerFactory.makeStreamer(FileFormat.testInput(shard.input.toString(), FileFormat.FASTA, null, true, false), null, true, -1);
		input.start();
		try{
			for(ListNum<Read> list=input.nextList(); list!=null && list.size()>0; list=input.nextList()){
				for(Read read : list){
					require(read.mate==null && read.bases!=null && read.bases.length>0, "Tail assignment requires nonempty unpaired proteins");
					final ProteinSequence query=new ProteinSequence(HbmCompetitiveAssign.queryId(read.id), read.bases);
					require(prefix==null || (query.id.startsWith(prefix) && seen.add(query.id)), "Duplicate query or wrong genome namespace: "+query.id);
					require(query.length()==read.bases.length, "Tail assignment cannot strip a terminal stop marker: "+query.id);
					validateTail(query);
					ProteinSearcher.scoreFamiliesF4(query, library.sidecar, shortlist);
					HbmCompetitiveAssign.top(shortlist.f4, shortlist.topIdx);
					final int rank=library.gate.select(query.enc, shortlist.topIdx, limit, selection);
					row.clear().append(query.id).tab();
					if(rank>=0){
						assigned++; familyCounts[rank]++;
						final HbmPositionModel.TailResult alignment=selection.alignment;
						final byte[] ref=library.refs[rank];
						int q=alignment.queryStart, r=alignment.start, lastQ=-1, lastR=-1, pairs=0, matches=0, terminal=0;
						for(byte op : alignment.path){
							if(op=='m'){
								lastQ=q; lastR=r;
								pairs++; if(query.enc[q]==ref[r]){matches++;}
								if(q>=query.length()-20 && r>=ref.length-20){terminal++;}q++; r++;
							}else if(op=='I'){q++;}else if(op=='D'){r++;}
							else{throw new IllegalArgumentException("Unexpected accepted tail trace operation");}
						}
						require(q==query.length() && r==ref.length, "Accepted trace must consume both tail endpoints before publication");
						if(library.policy==HbmTailProfileRefiner.TailPolicy.FINAL_BIDIRECTIONAL){
							require(lastQ==query.length()-1 && lastR==ref.length-1, "Final tail assignments require terminal paired endpoints");
						}
						row.append("ASSIGNED").tab().append(rank).tab().append(library.ids[rank]).tab().append(selection.score64)
							.tab().append(query.length()).tab().append(ref.length).tab().append(alignment.queryStart).tab().append(alignment.start)
							.tab();
						if(library.policy==HbmTailProfileRefiner.TailPolicy.FINAL_BIDIRECTIONAL){row.append(lastQ).tab().append(lastR).tab();}
						row.append(pairs).tab().append(matches).tab().append(alignment.path.length).tab().append(terminal);
					}else{
						row.append("ORPHAN\t-1\t.\t0").tab().append(query.length()).append("\t0\t-1\t-1");
						if(library.policy==HbmTailProfileRefiner.TailPolicy.FINAL_BIDIRECTIONAL){row.append("\t-1\t-1");}
						row.append("\t0\t0\t0\t0");
					}
					calls.println(row.tab().append(selection.scored).tab().append(selection.recorded).tab().append(selection.rejected).tab().append(selection.nonpositive));
					count++; residues=Math.addExact(residues, query.length());
					scored=Math.addExact(scored, selection.scored); recorded=Math.addExact(recorded, selection.recorded);
					rejected=Math.addExact(rejected, selection.rejected); nonpositive=Math.addExact(nonpositive, selection.nonpositive);
					if(count%1000==0){System.err.println("HBM_TAIL_ASSIGN_PROGRESS queries="+count+" assigned="+assigned);}
				}
			}
		}finally{input.close(); HmmComparisonData.close(calls);}
		stage("search_rows", shard.genome, 1, count, assigned, library.ids.length, limit,
			"input_stream_to_assignments_close_before_postverification", searchStart);
		require(!input.errorState(), "Tail assignment input reader failed");
		require(count==shard.expected, "Tail query count differs: expected="+shard.expected+" observed="+count);
		require(scored==Math.multiplyExact(count, limit) && recorded==assigned+rejected && recorded<=scored-nonpositive,
			"Every shortlisted score and recorded trace must have one conserved outcome");
		if(verifySharedBeforePublish){verifyReference(options);}
		HbmProfilePilot.requireHash(shard.input.toString(), shard.inputSha80);
		final ByteStreamWriter counts=HmmComparisonData.writer(shard.result.resolve("family_counts.tsv").toString());
		counts.println("rank\tfamily_id\tassigned"); long sum=0;
		for(int i=0; i<familyCounts.length; i++){
			sum=Math.addExact(sum, familyCounts[i]); counts.println(row.clear().append(i).tab().append(library.ids[i]).tab().append(familyCounts[i]));
		}
		HmmComparisonData.close(counts); require(sum==assigned && assigned<=count, "Family totals must conserve assigned queries");
		final ByteStreamWriter summary=HmmComparisonData.writer(shard.result.resolve("summary.tsv").toString());
		summary.println("queries\tassigned\torphaned\tresidues\tfamilies\tshortlist\tscored\trecorded\trejected\tnonpositive");
		summary.println(row.clear().append(count).tab().append(assigned).tab().append(count-assigned).tab().append(residues)
			.tab().append(library.ids.length).tab().append(limit).tab().append(scored).tab().append(recorded).tab().append(rejected).tab().append(nonpositive));
		HmmComparisonData.close(summary);
		final ByteStreamWriter recipe=HmmComparisonData.writer(shard.result.resolve("recipe.tsv").toString());
		recipe.println("in_sha80\t"+shard.inputSha80);
		for(String key : REFERENCE){recipe.println(key+"_sha80\t"+options.get(key+"sha80"));}
		if(hasPackedReference(options)){for(String key : PACKED_REFERENCE){recipe.println(key+"_sha80\t"+options.get(key+"sha80"));}}
		recipe.println("shortlist_requested\t"+requested+"\nshortlist_effective\t"+limit+
			"\nranking\texisting_F4_then_score64_ascii_family\nacceptance\t"+acceptance(library.policy)+
			"\nprofile\tlogodds_beta0.01_clip-4_11_gap4");
		HmmComparisonData.close(recipe);
		final ByteStreamWriter pass=HmmComparisonData.writer(shard.result.resolve("PASS").toString());
		pass.println("HBM_TAIL_ASSIGN_PASS"); HmmComparisonData.close(pass);
	}

	/** Owns detached scores/references after verifying the complete singleton and index bindings. */
	static final class Library {
		Library(HashMap<String,String> options,HbmTailProfileRefiner.TailPolicy policy_) throws Exception{
			libraryLoadsForTest++;
			final List<String[]> models=HmmComparisonData.rows(options.get("models"));
			require(models.size()>1 && String.join("\t", models.get(0)).equals(HbmTailLibrary.MODELS_HEADER), "Unsupported pinned tail-model manifest");
			final List<ProteinSequence> sequences=HbmTailLibrary.readTails(options.get("consensus"));
			final int n=sequences.size(); require(n==models.size()-1, "Tail consensus and model roster sizes differ");
			ids=FamilyShortlistSidecarBuilder.loadConsensusFamilyOrder(options.get("consensus"));
			require(ids.length==n, "Tail index rank count differs from consensus records");
			refs=new byte[n][]; profiles=new HbmPositionModel[n];
			final double[] background=NewFamilyHbmBuilder.readBackground(options.get("background"));
			final HashMap<String,byte[]> sequenceMap=new HashMap<String,byte[]>();
			final HashSet<String> seen=new HashSet<String>();
			for(int i=0; i<n; i++){
				final String[] row=models.get(i+1); final ProteinSequence ref=sequences.get(i);
				require(row.length==9 && row[0].equals(Integer.toString(i)) && row[1].equals(ids[i]) && ids[i].equals(ref.id) && seen.add(ref.id)
					&& sequenceMap.put(ref.id, ref.enc)==null, "Tail model/index identity differs at rank "+i);
				final long count=Long.parseLong(row[3]), residues=Long.parseLong(row[4]);
				require(count>=3 && count<=Integer.MAX_VALUE && residues>=100L*count && residues<=200L*count, "Invalid retained tail membership dimensions");
				refs[i]=ref.enc;
			}
			if(hasPackedReference(options)){
				final Path hbm=Paths.get(options.get("hbm"));
				require(!HbmCompactTextReader.matches(hbm) && !HbmDenseTextBundle.matches(hbm) && !HbmA48Codec.isEncoded(hbm),
					"Packed tail assignment requires native unpruned model.mqhb input");
				final byte[][] provenance=HbmBundleLoader.loadSemanticProvenance(options.get("provenance"));
				verifyPackedProvenance(options, provenance);
				final HbmPositionModel[] loaded=HbmBundleLoader.load(hbm, Arrays.asList(ids), id->sequenceMap.get(id), provenance)
					.positionModels("logodds", .01, true, -4, background);
				require(loaded.length==profiles.length, "Packed tail model count differs from consensus");
				for(int i=0; i<loaded.length; i++){profiles[i]=loaded[i];}
			}else{
				for(int i=0; i<n; i++){
					final String[] row=models.get(i+1); final ProteinSequence ref=sequences.get(i);
					final Path dir=Paths.get(row[2]); require(dir.isAbsolute(), "Tail singleton directory must be absolute");
					for(int j=0; j<ARTIFACTS.length; j++){HbmProfilePilot.requireHash(dir.resolve(ARTIFACTS[j]).toString(), row[5+j]);}
					HbmProfilePilot.requireHash(dir.resolve("background.tsv").toString(), options.get("backgroundsha80"));
					final List<ProteinSequence> singleton=HbmTailLibrary.readTails(dir.resolve("consensus.faa").toString());
					require(singleton.size()==1 && singleton.get(0).id.equals(ref.id) && Arrays.equals(singleton.get(0).enc, ref.enc), "Tail singleton and ranked consensus differ");
					final HbmBundleLoader.Loaded loaded=HbmBundleLoader.load(dir.resolve("refined.mqhb"), Collections.singletonList(ref.id), id->ref.enc,
						HbmBundleLoader.loadSemanticProvenance(dir.resolve("provenance.tsv").toString()));
					profiles[i]=loaded.positionModels("logodds", .01, true, -4, background)[0];
					for(int j=0; j<ARTIFACTS.length; j++){HbmProfilePilot.requireHash(dir.resolve(ARTIFACTS[j]).toString(), row[5+j]);}
				}
			}
			sidecar=FamilyShortlistSidecar.load(options.get("sidecar"), options.get("consensus"), options.get("kmers"));
			require(sidecar.schemaVersion==1 && sidecar.coveringK==5 && Arrays.equals(ids, sidecar.repIds), "Tail assignment needs matching five-mer indexes without production thresholds");
			policy=policy_;
			gate=new HbmTailCompetitiveGate(ids, refs, profiles, policy);
		}
		final String[] ids;
		final byte[][] refs;
		final HbmPositionModel[] profiles;
		final FamilyShortlistSidecar sidecar;
		final HbmTailCompetitiveGate gate;
		final HbmTailProfileRefiner.TailPolicy policy;
	}

	private static ArrayList<Shard> loadShards(Path path) throws Exception{
		final List<String[]> rows=HmmComparisonData.rows(path.toString());
		require(rows.size()>1 && String.join("\t", rows.get(0)).equals(SHARD_HEADER), "Unexpected tail-assignment shard roster schema");
		final HashSet<String> genomes=new HashSet<String>();
		final HashSet<Path> results=new HashSet<Path>();
		final ArrayList<Shard> shards=new ArrayList<Shard>(rows.size()-1);
		for(int i=1; i<rows.size(); i++){
			final String[] row=rows.get(i);
			require(row.length==5 && row[0].matches("[A-Za-z0-9_.-]+") && genomes.add(row[0]), "Malformed or duplicate genome namespace");
			final long count=Long.parseLong(row[1]); require(count>0, "Empty declared tail shard");
			final Path input=Paths.get(row[2]), result=Paths.get(row[4]);
			require(input.isAbsolute() && result.isAbsolute() && results.add(result.toAbsolutePath().normalize()),
				"Shard paths must be absolute with unique result directories");
			require(!Files.exists(result), "Tail assignment output must be fresh: "+result);
			DigestSuffix.requireSuffix(row[3], "query input pin");
			shards.add(new Shard(row[0], count, input, row[3], result));
		}
		return shards;
	}
	private static void checkBatchOutput(Path out,ArrayList<Shard> shards){
		final Path root=out.toAbsolutePath().normalize();
		for(Shard shard : shards){
			require(!shard.result.toAbsolutePath().normalize().equals(root), "Batch and shard outputs must be distinct: "+root);
		}
	}
	private static int workers(HashMap<String,String> options,int shards){
		final String t=options.get("t"), threads=options.get("threads");
		require(t==null || threads==null || t.equals(threads), "t and threads must agree");
		final String value=t==null ? threads : t;
		final int workers=value==null ? 1 : Integer.parseInt(value);
		require(workers>0, "threads must be positive");
		return Math.min(workers, shards);
	}
	private static void checkParallelOutputs(ArrayList<Shard> shards){
		final HashSet<Path> paths=new HashSet<Path>(shards.size()*2);
		for(Shard shard : shards){paths.add(shard.result.toAbsolutePath().normalize());}
		for(Path path : paths){
			for(Path parent=path.getParent(); parent!=null; parent=parent.getParent()){
				require(!paths.contains(parent), "Overlapping shard output paths cannot run in parallel: "+parent+" and "+path);
			}
		}
	}
	private static HashSet<String> allowed(String... names){
		final HashSet<String> allowed=new HashSet<String>(Arrays.asList(names));
		for(String key : REFERENCE){allowed.add(key); allowed.add(key+"sha80");}
		for(String key : PACKED_REFERENCE){allowed.add(key); allowed.add(key+"sha80");}
		return allowed;
	}
	static void verifyReference(HashMap<String,String> options){
		for(String key : REFERENCE){HbmProfilePilot.requireHash(HmmComparisonData.required(options, key), HmmComparisonData.required(options, key+"sha80"));}
		if(hasPackedReference(options)){
			for(String key : PACKED_REFERENCE){HbmProfilePilot.requireHash(HmmComparisonData.required(options, key), HmmComparisonData.required(options, key+"sha80"));}
		}
	}
	private static boolean hasPackedReference(HashMap<String,String> options){
		boolean found=false;
		for(String key : PACKED_REFERENCE){
			found|=(options.containsKey(key) || options.containsKey(key+"sha80"));
		}
		if(found){
			for(String key : PACKED_REFERENCE){
				require(options.containsKey(key) && options.containsKey(key+"sha80"), "Packed HBM option set requires hbm/provenance/packrecipe with sha80 pins");
			}
		}
		return found;
	}
	private static void verifyPackedProvenance(HashMap<String,String> options,byte[][] provenance) throws Exception{
		require(DigestSuffix.fromDigest(provenance[11]).equals(options.get("consensussha80")), "Packed consensus provenance differs from consensussha80");
		require(DigestSuffix.fromDigest(provenance[12]).equals(options.get("modelssha80")), "Packed source-corpus provenance differs from modelssha80");
		require(DigestSuffix.fromDigest(provenance[13]).equals(options.get("modelssha80")), "Packed cluster-membership provenance differs from modelssha80");
		require(DigestSuffix.fromDigest(provenance[14]).equals(options.get("packrecipesha80")), "Packed member-policy provenance differs from packrecipesha80");
		require(DigestSuffix.fromDigest(provenance[15]).equals(options.get("modelssha80")), "Packed member-manifest provenance differs from modelssha80");
		verifyRecipePins(options);
	}
	private static void verifyRecipePins(HashMap<String,String> options){
		final String[] keys={"models_sha80", "consensus_sha80", "background_sha80"};
		final String[] pins={options.get("modelssha80"), options.get("consensussha80"), options.get("backgroundsha80")};
		final boolean[] found=new boolean[keys.length];
		final ByteFile bf=ByteFile.makeByteFile(options.get("packrecipe"), true);
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				final String s=new String(line, StandardCharsets.UTF_8);
				final int eq=s.indexOf('=');
				if(eq<0){continue;}
				final String key=s.substring(0, eq);
				for(int i=0; i<keys.length; i++){
					if(key.equals(keys[i])){
						require(!found[i], "Packed tail recipe duplicates input pin "+key);
						found[i]=true;
						require(s.substring(eq+1).equals(pins[i]), "Packed tail recipe input pin differs for "+key);
					}
				}
			}
		}finally{require(!bf.close(), "Packed tail recipe read failed");}
		for(int i=0; i<keys.length; i++){require(found[i], "Packed tail recipe input pin differs for "+keys[i]);}
	}
	private static void validateTail(ProteinSequence query){
		require(query.length()>=100 && query.length()<=200, "Tail query must contain100..200 residues: "+query.id);
		for(byte b : query.enc){require(b>=0 && b<20, "Tail query contains an ambiguous or stop residue: "+query.id);}
	}
	private static HbmTailProfileRefiner.TailPolicy parsePolicy(String value){
		if("construction".equalsIgnoreCase(value)){return HbmTailProfileRefiner.TailPolicy.LEGACY_CONSTRUCTION;}
		if("final_bidirectional".equalsIgnoreCase(value)){return HbmTailProfileRefiner.TailPolicy.FINAL_BIDIRECTIONAL;}
		throw new IllegalArgumentException("Unknown tail policy: "+value);
	}
	private static String acceptance(HbmTailProfileRefiner.TailPolicy policy){
		switch(policy){
			case LEGACY_CONSTRUCTION: return "HbmTailProfileRefiner_construction";
			case FINAL_BIDIRECTIONAL: return "HbmTailProfileRefiner_final_bidirectional";
			default: throw new AssertionError(policy);
		}
	}
	private static void stage(String stage,String genome,int shards,long queries,long assigned,int families,int shortlist,
			String boundary,long startNanos){
		System.err.println("HBM_TAIL_ASSIGN_STAGE stage="+stage+" genome="+genome+" shards="+shards+
			" queries="+queries+" assigned="+assigned+" families="+families+" shortlist="+shortlist+
			" elapsed_seconds="+String.format(Locale.ROOT, "%.6f", (System.nanoTime()-startNanos)/1000000000.0)+
			" boundary="+boundary);
	}
	public static int libraryLoadsForTest(){return libraryLoadsForTest;}
	private static void require(boolean ok,String message){if(!ok){throw new IllegalArgumentException(message);}}
	private static final class Shard {
		Shard(String genome_,long expected_,Path input_,String inputSha80_,Path result_){
			genome=genome_; expected=expected_; input=input_; inputSha80=inputSha80_; result=result_;
		}
		final String genome;
		final long expected;
		final Path input, result;
		final String inputSha80;
	}
	private static final class ProcessThread extends Thread {
		ProcessThread(ArrayList<Shard> shards_,HashMap<String,String> options_,int requested_,int limit_,Library library_,
				AtomicInteger next_,AtomicReference<Throwable> failure_){
			shards=shards_; options=options_; requested=requested_; limit=limit_; library=library_; next=next_; failure=failure_;
			setDaemon(true);
		}
		@Override
		public void run(){
			try{
				for(int idx=next.getAndIncrement(); idx<shards.size() && failure.get()==null; idx=next.getAndIncrement()){
					final Shard shard=shards.get(idx);
					assignOne(shard, shard.genome+"|", false, options, requested, limit, library);
				}
			}catch(Throwable e){failure.compareAndSet(null, e);}
		}
		final ArrayList<Shard> shards;
		final HashMap<String,String> options;
		final int requested, limit;
		final Library library;
		final AtomicInteger next;
		final AtomicReference<Throwable> failure;
	}
	private static int libraryLoadsForTest=0;
	private static final String[] REFERENCE={"models", "consensus", "sidecar", "kmers", "background"};
	private static final String[] PACKED_REFERENCE={"hbm", "provenance", "packrecipe"};
	private static final String[] ARTIFACTS={"refined.mqhb", "consensus.faa", "provenance.tsv", "retained.faa"};
	private static final String SHARD_HEADER="genome_id\tqueries\tinput\tinput_sha80\tresult";
	private static final String HEADER="query\tstatus\trank\tfamily_id\tscore64\tquery_length\tmodel_length\tquery_start\tmodel_start\tpaired\tmatches\tcolumns\t"+
		"terminal_paired\tscored\trecorded\trejected\tnonpositive";
	private static final String FINAL_HEADER="query\tstatus\trank\tfamily_id\tscore64\tquery_length\tmodel_length\tquery_start\tmodel_start\t"+
		"query_paired_end\tmodel_paired_end\tpaired\tmatches\tcolumns\tterminal_paired\tscored\trecorded\trejected\tnonpositive";
	private HbmTailAssign(){}
}
