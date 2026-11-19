package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

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

/** Records configured scorers for every shortlisted stop-tail candidate, without acceptance gates.
 * Nonpositive outcomes remain explicit; a bounded shard streams through one worker.
 * Calibration and final panels are separate pinned inputs, never inferred from scores.
 * @author Keqing
 */
public final class HbmTailEvidence {

	public static void main(String[] args){
		try{
			final PreParser pp=new PreParser(args, HbmTailEvidence.class, false);
			try{runForTest(pp.args);}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	/** Publishes PASS only after count conservation and a second check of all input pins. */
	public static void runForTest(String[] args) throws Exception{
		Shared.setThreads(1); Shared.AMINO_IN=true;
		final HashMap<String,String> options=HmmComparisonData.options(args);
		final int mode=parseMode(options.getOrDefault("mode", "both"));
		final String modeName=modeName(mode);
		require(!controlRequested(options) || controlSupplied(options), "control and controlsha80 must be supplied together");
		require(mode!=MODE_BOTH || controlSupplied(options), "Paired evidence mode requires control and controlsha80");
		final boolean verifyControl=mode==MODE_BOTH || controlSupplied(options);
		final int requested=Integer.parseInt(options.getOrDefault("shortlist", "50"));
		require(requested>=0, "Evidence requires a nonnegative shortlist");
		final boolean batch=options.containsKey("shards") || options.containsKey("shardssha80");
		final boolean single=options.containsKey("in") || options.containsKey("insha80") || options.containsKey("expected");
		require(batch!=single, "Evidence requires either single in/insha80/expected or batch shards/shardssha80");
		if(batch){runBatch(options, mode, modeName, verifyControl, requested);}
		else{runSingle(options, mode, modeName, verifyControl, requested);}
	}

	private static void runSingle(HashMap<String,String> options,int mode,String modeName,boolean verifyControl,int requested) throws Exception{
		final HashSet<String> allowed=new HashSet<String>(Arrays.asList("out", "expected", "shortlist", "mode"));
		for(String key : allowedPins(mode)){allowed.add(key); allowed.add(key+"sha80");}
		for(String key : options.keySet()){require(allowed.contains(key), "Unknown evidence option: "+key);}
		final long setupStart=System.nanoTime();
		final long expected=Long.parseLong(HmmComparisonData.required(options, "expected"));
		require(expected>0, "Evidence requires a positive query count");
		final Shard shard=new Shard("single", expected, Paths.get(HmmComparisonData.required(options, "in")),
			HmmComparisonData.required(options, "insha80"), Paths.get(HmmComparisonData.required(options, "out")));
		verifyInputs(options, mode, true);
		require(!Files.exists(shard.result), "Evidence output must be fresh: "+shard.result);
		final Library library=new Library(options, mode);
		final int limit=requested==0 ? library.ids.length : Math.min(requested, library.ids.length);
		stage("setup", modeName, expected, -1, -1, -1, library.ids.length, limit, library.methodCount, setupStart);
		final long equalCells=verifyControl ? library.verifyControl(options, modeName) : 0;
		scoreOne(shard, true, options, mode, modeName, requested, limit, library, equalCells);
	}

	private static void runBatch(HashMap<String,String> options,int mode,String modeName,boolean verifyControl,int requested) throws Exception{
		final HashSet<String> allowed=new HashSet<String>(Arrays.asList("out", "shards", "shardssha80", "shortlist", "mode", "t", "threads"));
		for(String key : allowedPins(mode)){
			if(!key.equals("in")){allowed.add(key); allowed.add(key+"sha80");}
		}
		for(String key : options.keySet()){require(allowed.contains(key), "Unknown evidence option: "+key);}
		final long setupStart=System.nanoTime();
		verifyInputs(options, mode, false);
		final Path shardPath=Paths.get(HmmComparisonData.required(options, "shards"));
		HbmProfilePilot.requireHash(shardPath.toString(), HmmComparisonData.required(options, "shardssha80"));
		final ArrayList<Shard> shards=loadShards(shardPath);
		long queries=0;
		for(Shard shard : shards){
			HbmProfilePilot.requireHash(shard.input.toString(), shard.inputSha80);
			queries=Math.addExact(queries, shard.expected);
		}
		final Path out=Paths.get(HmmComparisonData.required(options, "out"));
		require(!Files.exists(out), "Evidence batch output must be fresh: "+out);
		checkBatchOutput(out, shards);
		checkParallelOutputs(shards);
		final int workers=workers(options, shards.size());
		final Library library=new Library(options, mode);
		final int limit=requested==0 ? library.ids.length : Math.min(requested, library.ids.length);
		stage("setup", modeName, queries, -1, -1, -1, library.ids.length, limit, library.methodCount, setupStart);
		final long equalCells=verifyControl ? library.verifyControl(options, modeName) : 0;
		Files.createDirectory(out);
		final long batchStart=System.nanoTime();
		if(workers<2){
			for(Shard shard : shards){scoreOne(shard, false, options, mode, modeName, requested, limit, library, equalCells);}
		}else{scoreParallel(shards, options, mode, modeName, requested, limit, library, equalCells, workers);}
		final long elapsed=System.nanoTime()-batchStart;
		verifyInputs(options, mode, false);
		if(verifyControl){verifyControlInput(options);}
		HbmProfilePilot.requireHash(shardPath.toString(), options.get("shardssha80"));
		for(Shard shard : shards){HbmProfilePilot.requireHash(shard.input.toString(), shard.inputSha80);}
		final long rows=Math.multiplyExact(Math.multiplyExact(queries, limit), library.methodCount);
		batchStage(modeName, shards.size(), queries, rows, workers, elapsed);
		final ByteStreamWriter summary=HmmComparisonData.writer(out.resolve("summary.tsv").toString());
		summary.println("shards\tqueries\tfamilies\tshortlist\tmethods\trows\tworkers\telapsed_seconds");
		summary.println(new ByteBuilder().append(shards.size()).tab().append(queries).tab().append(library.ids.length).tab()
			.append(limit).tab().append(library.methodCount).tab().append(rows).tab().append(workers).tab()
			.append(String.format(Locale.ROOT, "%.6f", elapsed/1000000000.0)));
		HmmComparisonData.close(summary);
		final ByteStreamWriter recipe=HmmComparisonData.writer(out.resolve("recipe.tsv").toString());
		recipe.println("shards_sha80\t"+options.get("shardssha80"));
		for(String key : verifiedPins(options, mode)){if(!key.equals("in")){recipe.println(key+"_sha80\t"+options.get(key+"sha80"));}}
		recipe.println("mode\t"+modeName+"\nthreads_requested\t"+options.getOrDefault("t", options.getOrDefault("threads", "1"))+
			"\nthreads_effective\t"+workers+"\nshortlist_requested\t"+requested+"\nshortlist_effective\t"+limit+
			"\nprofile\t"+profileName(mode)+"\nscore_scale\t64\nterminal_window\t20\nterminal_gap_coordinate\tbefore-operation cursors; exhausted coordinate remains in terminal window\nacceptance\tnone; preserve all candidate outcomes\ncontrol\t"+
			(verifyControl ? "exact native/pruned positional emission equality; not graph-scorer equivalence" : "not_run"));
		HmmComparisonData.close(recipe);
		final ByteStreamWriter pass=HmmComparisonData.writer(out.resolve("PASS").toString());
		pass.println("HBM_TAIL_EVIDENCE_BATCH_PASS"); HmmComparisonData.close(pass);
	}

	private static void scoreParallel(ArrayList<Shard> shards,HashMap<String,String> options,int mode,String modeName,int requested,
			int limit,Library library,long equalCells,int workerCount) throws Exception{
		final AtomicInteger next=new AtomicInteger();
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>(workerCount);
		final ArrayList<ProcessThread> started=new ArrayList<ProcessThread>(workerCount);
		for(int i=0; i<workerCount; i++){
			workers.add(new ProcessThread(shards, options, mode, modeName, requested, limit, library, equalCells, next, failure));
		}
		try{
			for(ProcessThread pt : workers){
				pt.start();
				started.add(pt);
			}
		}catch(Throwable e){failure.compareAndSet(null, e);}
		finally{ThreadWaiter.waitForThreadsToFinish(started);}
		if(failure.get()!=null){throw new RuntimeException("Evidence worker failed", failure.get());}
	}

	private static void scoreOne(Shard shard,boolean verifySharedBeforePublish,HashMap<String,String> options,int mode,
			String modeName,int requested,int limit,Library library,long equalCells) throws Exception{
		final HashSet<String> seen=new HashSet<String>();
		final ProteinSearcher.ShortlistScratch shortlist=new ProteinSearcher.ShortlistScratch(library.sidecar.dims, library.ids.length, limit);
		Files.createDirectory(shard.result);
		final ByteStreamWriter writer=HmmComparisonData.writer(shard.result.resolve("candidates.tsv").toString());
		writer.println(HEADER);
		final Streamer input=StreamerFactory.makeStreamer(FileFormat.testInput(shard.input.toString(), FileFormat.FASTA, null, true, false), null, true, -1);
		final ByteBuilder row=new ByteBuilder();
		long queries=0, residues=0, traces=0, rows=0;
		final long searchStart=System.nanoTime();
		input.start();
		try{
			for(ListNum<Read> list=input.nextList(); list!=null && list.size()>0; list=input.nextList()){
				for(Read read : list){
					require(read.mate==null && read.bases!=null, "Evidence requires unpaired amino-acid queries");
					final ProteinSequence query=new ProteinSequence(HbmCompetitiveAssign.queryId(read.id), read.bases);
					require(query.length()==read.bases.length, "Evidence cannot strip terminal stop markers");
					validateTail(query); require(seen.add(query.id), "Duplicate query ID: "+query.id);
					ProteinSearcher.scoreFamiliesF4(query, library.sidecar, shortlist);
					HbmCompetitiveAssign.top(shortlist.f4, shortlist.topIdx);
					for(int slot=0; slot<limit; slot++){
						final int rank=shortlist.topIdx[slot];
						for(int method=0; method<2; method++){
							if(library.profiles[method]==null){continue;}
							final HbmPositionModel model=library.profiles[method][rank];
							final byte[] ref=library.refs[rank];
							final HbmPositionModel.TailResult alignment=model.alignTail(query.enc, 200, true);
							row.clear().append(query.id).tab().append(METHODS[method]).tab().append(rank).tab().append(library.ids[rank])
								.tab().append(slot).tab().append(query.length()).tab().append(ref.length).tab();
							if(alignment==null){row.append("NONPOSITIVE\t0\t-1\t-1\t0\t0\t0\t0\t0\t0\t0\t0\t0\t.");}
							else{
								final Metrics m=measure(query.enc, ref, model, alignment); traces++;
								row.append("TRACE").tab().append(alignment.score).tab().append(alignment.queryStart).tab().append(alignment.start)
									.tab().append(m.paired).tab().append(m.matches).tab().append(alignment.path.length).tab().append(m.insertions)
									.tab().append(m.deletions).tab().append(m.terminalPaired).tab().append(m.terminalMatches)
									.tab().append(m.terminalScore64).tab().append(m.trailingGap).tab().append(alignment.path);
							}
							writer.println(row); rows++;
						}
					}
					queries++; residues=Math.addExact(residues, query.length());
					if(queries%100==0){System.err.println("HBM_TAIL_EVIDENCE_PROGRESS queries="+queries+" rows="+rows);}
				}
			}
		}finally{input.close(); HmmComparisonData.close(writer);}
		require(!input.errorState(), "Evidence query reader failed");
		require(queries==shard.expected, "Evidence query count differs: expected="+shard.expected+" actual="+queries);
		stage("search_rows", modeName, verifySharedBeforePublish ? null : shard.id, queries, residues, rows, traces,
			library.ids.length, limit, library.methodCount, searchStart);
		require(rows==Math.multiplyExact(Math.multiplyExact(queries, limit), library.methodCount) && traces<=rows,
			"All requested scorers must retain every candidate outcome");
		if(verifySharedBeforePublish){verifyInputs(options, mode, false);}
		if(verifySharedBeforePublish && controlSupplied(options)){verifyControlInput(options);}
		HbmProfilePilot.requireHash(shard.input.toString(), shard.inputSha80);
		final ByteStreamWriter summary=HmmComparisonData.writer(shard.result.resolve("summary.tsv").toString());
		summary.println("queries\tresidues\tfamilies\tshortlist\trows\ttraces\tnonpositive\tequal_emission_cells");
		summary.println(row.clear().append(queries).tab().append(residues).tab().append(library.ids.length).tab().append(limit)
			.tab().append(rows).tab().append(traces).tab().append(rows-traces).tab().append(equalCells));
		HmmComparisonData.close(summary);
		final ByteStreamWriter recipe=HmmComparisonData.writer(shard.result.resolve("recipe.tsv").toString());
		recipe.println("in_sha80\t"+shard.inputSha80);
		for(String key : verifiedPins(options, mode)){if(!key.equals("in")){recipe.println(key+"_sha80\t"+options.get(key+"sha80"));}}
		if(mode!=MODE_BOTH){recipe.println("mode\t"+modeName);}
		recipe.println("shortlist_requested\t"+requested+"\nshortlist_effective\t"+limit+"\nprofile\t"+profileName(mode)+"\nscore_scale\t64\nterminal_window\t20\nterminal_gap_coordinate\tbefore-operation cursors; exhausted coordinate remains in terminal window\nacceptance\tnone; preserve all candidate outcomes\ncontrol\t"+
			(equalCells>0 ? "exact native/pruned positional emission equality; not graph-scorer equivalence" : "not_run"));
		HmmComparisonData.close(recipe);
		final ByteStreamWriter pass=HmmComparisonData.writer(shard.result.resolve("PASS").toString());
		pass.println("HBM_TAIL_EVIDENCE_PASS"); HmmComparisonData.close(pass);
	}

	/** Reconstructs score and geometry from the complete suffix trace. Terminal score
	 * includes an operation when either BEFORE-operation cursor is in its final20
	 * window. An exhausted cursor remains in that window, so terminal gaps count.
	 * Clipped prefixes never contribute. Saved terminal score divides by64*20 later.
	 */
	public static Metrics measure(byte[] query,byte[] ref,HbmPositionModel model,HbmPositionModel.TailResult alignment){
		require(query!=null && ref!=null && model!=null && alignment!=null && alignment.path!=null && alignment.path.length>0,
			"Trace metrics require both sequences and a recorded nonempty alignment");
		require(model.length()==ref.length && alignment.queryEnd==query.length-1 && alignment.end==ref.length-1,
			"Trace geometry must consume the supplied query/model endpoints");
		int q=alignment.queryStart, r=alignment.start, score=0;
		require(q>=0 && q<query.length && r>=0 && r<ref.length, "Suffix start lies outside a sequence");
		final Metrics result=new Metrics();
		for(byte op : alignment.path){
			final boolean terminal=q>=query.length-20 || r>=ref.length-20;
			int value;
			if(op=='m'){
				require(q<query.length && r<ref.length, "Paired operation overruns a sequence");
				value=model.scoreAt(query[q], r); result.paired++;
				if(query[q]==ref[r]){result.matches++;}
				if(q>=query.length-20 && r>=ref.length-20){
					result.terminalPaired++; if(query[q]==ref[r]){result.terminalMatches++;}
				}
				q++; r++; result.trailingGap=0;
			}else if(op=='I'){
				require(q<query.length, "Insertion overruns query"); q++; result.insertions++; result.trailingGap++; value=-HbmPositionModel.GAP;
			}else if(op=='D'){
				require(r<ref.length, "Deletion overruns model"); r++; result.deletions++; result.trailingGap++; value=-HbmPositionModel.GAP;
			}else{throw new IllegalArgumentException("Unknown suffix trace operation: "+op);}
			score=Math.addExact(score, value);
			if(terminal){result.terminalScore64=Math.addExact(result.terminalScore64, value);}
		}
		require(q==query.length && r==ref.length && score==alignment.score, "Reconstructed trace must match DP score and both consumed endpoints");
		assert(result.paired+result.insertions+result.deletions==alignment.path.length) : "Every path operation must have one geometry category";
		return result;
	}

	/** Counts equal cells, failing at the first disagreement rather than hiding pruning effects. */
	public static long requireEqualMatrices(HbmPositionModel[] candidate,HbmPositionModel[] control){
		require(candidate.length==control.length && candidate.length>0, "Candidate/control family counts differ or are empty");
		long count=0;
		for(int i=0; i<candidate.length; i++){
			require(candidate[i].length()==control[i].length(), "Candidate/control lengths differ at rank "+i);
			for(int aa=0; aa<22; aa++){for(int p=0; p<candidate[i].length(); p++){
				require(candidate[i].scoreAt((byte)aa, p)==control[i].scoreAt((byte)aa, p), "Candidate/control emission differs at rank="+i+" aa="+aa+" column="+p); count++;
			}}
		}
		return count;
	}

	/** Owns detached immutable matrices after loading the verified compact and native models. */
	private static final class Library {
		Library(HashMap<String,String> options,int mode) throws Exception{
			libraryLoadsForTest++;
			final List<ProteinSequence> sequences=HbmTailLibrary.readTails(options.get("consensus"));
			ids=FamilyShortlistSidecarBuilder.loadConsensusFamilyOrder(options.get("consensus"));
			require(ids.length>0 && ids.length==sequences.size(), "Consensus roster must be nonempty and exact");
			refs=new byte[ids.length][]; byId=new HashMap<String,byte[]>();
			for(int i=0; i<ids.length; i++){
				final ProteinSequence ref=sequences.get(i); validateTail(ref);
				require(ids[i].equals(ref.id) && byId.put(ref.id, ref.enc)==null, "Consensus identity/order differs at rank "+i); refs[i]=ref.enc;
			}
			profiles=new HbmPositionModel[2][];
			HbmBundleLoader.Loaded candidate=null;
			byte[][] provenanceLoaded=null;
			double[] backgroundLoaded=null;
			if((mode&MODE_HBM)!=0){
				provenanceLoaded=HbmBundleLoader.loadSemanticProvenance(options.get("provenance"));
				backgroundLoaded=NewFamilyHbmBuilder.readBackground(options.get("background"));
				candidate=HbmBundleLoader.load(Paths.get(options.get("hbm")), Arrays.asList(ids), byId::get, provenanceLoaded);
				profiles[METHOD_HBM]=candidate.positionModels("logodds", .01, true, -4, backgroundLoaded);
			}
			if((mode&MODE_BLOSUM)!=0){
				profiles[METHOD_BLOSUM]=(candidate==null ? consensusProfiles(refs) :
					candidate.positionModels("consensus", .01, true, -4, backgroundLoaded));
			}
			sidecar=FamilyShortlistSidecar.load(options.get("sidecar"), options.get("consensus"), options.get("kmers"));
			require(sidecar.schemaVersion==1 && sidecar.coveringK==5 && Arrays.equals(ids, sidecar.repIds), "Evidence requires matching five-mer F4 family order");
			provenance=provenanceLoaded; background=backgroundLoaded;
			methodCount=Integer.bitCount(mode);
		}
		long verifyControl(HashMap<String,String> options,String modeName) throws Exception{
			final long start=System.nanoTime();
			verifyControlInput(options);
			final HbmBundleLoader.Loaded control=HbmBundleLoader.load(Paths.get(options.get("control")), Arrays.asList(ids), byId::get, provenance);
			final long equal=requireEqualMatrices(profiles[METHOD_HBM], control.positionModels("logodds", .01, true, -4, background));
			stage("control", modeName, -1, -1, -1, -1, ids.length, -1, 1, start);
			return equal;
		}
		final String[] ids;
		final byte[][] refs;
		final HashMap<String,byte[]> byId;
		final byte[][] provenance;
		final double[] background;
		final HbmPositionModel[][] profiles;
		final FamilyShortlistSidecar sidecar;
		final int methodCount;
	}
	public static final class Metrics {
		public int paired, matches, insertions, deletions, terminalPaired, terminalMatches, terminalScore64, trailingGap;
	}
	private static HbmPositionModel[] consensusProfiles(byte[][] refs){
		final HbmPositionModel[] profiles=new HbmPositionModel[refs.length];
		for(int i=0; i<refs.length; i++){profiles[i]=HbmPositionModel.fromConsensus(refs[i]);}
		return profiles;
	}
	private static ArrayList<Shard> loadShards(Path path) throws Exception{
		final List<String[]> rows=HmmComparisonData.rows(path.toString());
		require(rows.size()>1 && String.join("\t", rows.get(0)).equals(SHARD_HEADER), "Unexpected evidence shard roster schema");
		final HashSet<String> ids=new HashSet<String>();
		final HashSet<Path> results=new HashSet<Path>();
		final ArrayList<Shard> shards=new ArrayList<Shard>(rows.size()-1);
		for(int i=1; i<rows.size(); i++){
			final String[] row=rows.get(i);
			require(row.length==5 && row[0].matches("[A-Za-z0-9_.-]+") && ids.add(row[0]), "Malformed or duplicate shard ID");
			final long count=Long.parseLong(row[1]); require(count>0, "Empty declared evidence shard");
			final Path input=Paths.get(row[2]), result=Paths.get(row[4]);
			require(input.isAbsolute() && result.isAbsolute() && results.add(result.toAbsolutePath().normalize()),
				"Shard paths must be absolute with unique result directories");
			require(!Files.exists(result), "Evidence output must be fresh: "+result);
			DigestSuffix.requireSuffix(row[3], "query input pin");
			shards.add(new Shard(row[0], count, input, row[3], result));
		}
		return shards;
	}
	private static void checkBatchOutput(Path out,ArrayList<Shard> shards){
		final Path root=out.toAbsolutePath().normalize();
		for(Shard shard : shards){
			final Path result=shard.result.toAbsolutePath().normalize();
			require(!result.equals(root), "Batch and shard outputs must be distinct: "+root);
			require(!isParent(root, result) && !isParent(result, root), "Batch and shard outputs must not overlap: "+root+" and "+result);
		}
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
	private static boolean isParent(Path parent,Path child){
		for(Path p=child.getParent(); p!=null; p=p.getParent()){if(p.equals(parent)){return true;}}
		return false;
	}
	private static int workers(HashMap<String,String> options,int shards){
		final String t=options.get("t"), threads=options.get("threads");
		require(t==null || threads==null || t.equals(threads), "t and threads must agree");
		final String value=t==null ? threads : t;
		final int workers=value==null ? 1 : Integer.parseInt(value);
		require(workers>0, "threads must be positive");
		return Math.min(workers, shards);
	}
	private static void verifyInputs(HashMap<String,String> options,int mode,boolean includeIn){
		for(String key : requiredPins(mode)){
			if(includeIn || !key.equals("in")){
				HbmProfilePilot.requireHash(HmmComparisonData.required(options, key), HmmComparisonData.required(options, key+"sha80"));
			}
		}
	}
	private static void verifyControlInput(HashMap<String,String> options){
		HbmProfilePilot.requireHash(HmmComparisonData.required(options, "control"), HmmComparisonData.required(options, "controlsha80"));
	}
	private static String[] requiredPins(int mode){
		return mode==MODE_BLOSUM ? CONSENSUS_PINNED : mode==MODE_HBM ? HBM_PINNED : BOTH_PINNED;
	}
	private static String[] allowedPins(int mode){return mode==MODE_BLOSUM ? CONSENSUS_PINNED : HBM_CONTROL_PINNED;}
	private static String[] verifiedPins(HashMap<String,String> options,int mode){
		return ((mode&MODE_HBM)!=0 && controlSupplied(options)) ? HBM_CONTROL_PINNED : requiredPins(mode);
	}
	private static boolean controlRequested(HashMap<String,String> options){
		return options.containsKey("control") || options.containsKey("controlsha80");
	}
	private static boolean controlSupplied(HashMap<String,String> options){
		return options.containsKey("control") && options.containsKey("controlsha80");
	}
	private static int parseMode(String value){
		if("both".equalsIgnoreCase(value)){return MODE_BOTH;}
		if("blosum".equalsIgnoreCase(value)){return MODE_BLOSUM;}
		if("hbm".equalsIgnoreCase(value)){return MODE_HBM;}
		throw new IllegalArgumentException("Unknown evidence mode: "+value);
	}
	private static String modeName(int mode){return mode==MODE_BOTH ? "both" : (mode==MODE_BLOSUM ? "blosum" : "hbm");}
	private static String profileName(int mode){return mode==MODE_BLOSUM ? "blosum62_consensus" : "logodds_beta0.01_clip-4_11_gap4";}
	private static void stage(String stage,String mode,long queries,long residues,long rows,long traces,
			int families,int shortlist,int methods,long startNanos){
		stage(stage, mode, null, queries, residues, rows, traces, families, shortlist, methods, startNanos);
	}
	private static void stage(String stage,String mode,String shardId,long queries,long residues,long rows,long traces,
			int families,int shortlist,int methods,long startNanos){
		System.err.println("HBM_TAIL_EVIDENCE_STAGE stage="+stage+" mode="+mode+
			(shardId==null ? "" : " shard_id="+shardId)+" queries="+queries+
			" residues="+residues+" rows="+rows+" traces="+traces+" families="+families+
			" shortlist="+shortlist+" methods="+methods+" elapsed_seconds="+
			String.format(Locale.ROOT, "%.6f", (System.nanoTime()-startNanos)/1000000000.0));
	}
	private static void batchStage(String mode,int shards,long queries,long rows,int workers,long elapsedNanos){
		System.err.println("HBM_TAIL_EVIDENCE_BATCH stage=search_shards mode="+mode+" shards="+shards+
			" queries="+queries+" rows="+rows+" workers="+workers+" elapsed_seconds="+
			String.format(Locale.ROOT, "%.6f", elapsedNanos/1000000000.0));
	}
	public static int libraryLoadsForTest(){return libraryLoadsForTest;}
	private static void validateTail(ProteinSequence query){
		require(query.length()>=100 && query.length()<=200, "Evidence tail must contain100..200 residues: "+query.id);
		for(byte b : query.enc){require(b>=0 && b<20, "Evidence tail contains an ambiguous or stop residue: "+query.id);}
	}
	private static void require(boolean ok,String message){if(!ok){throw new IllegalArgumentException(message);}}
	private static final class Shard {
		Shard(String id_,long expected_,Path input_,String inputSha80_,Path result_){
			id=id_; expected=expected_; input=input_; inputSha80=inputSha80_; result=result_;
		}
		final String id;
		final long expected;
		final Path input, result;
		final String inputSha80;
	}
	private static final class ProcessThread extends Thread {
		ProcessThread(ArrayList<Shard> shards_,HashMap<String,String> options_,int mode_,String modeName_,int requested_,int limit_,
				Library library_,long equalCells_,AtomicInteger next_,AtomicReference<Throwable> failure_){
			shards=shards_; options=options_; mode=mode_; modeName=modeName_; requested=requested_; limit=limit_; library=library_;
			equalCells=equalCells_; next=next_; failure=failure_; setDaemon(true);
		}
		@Override
		public void run(){
			try{
				for(int idx=next.getAndIncrement(); idx<shards.size() && failure.get()==null; idx=next.getAndIncrement()){
					scoreOne(shards.get(idx), false, options, mode, modeName, requested, limit, library, equalCells);
				}
			}catch(Throwable e){failure.compareAndSet(null, e);}
		}
		final ArrayList<Shard> shards;
		final HashMap<String,String> options;
		final int mode, requested, limit;
		final String modeName;
		final Library library;
		final long equalCells;
		final AtomicInteger next;
		final AtomicReference<Throwable> failure;
	}
	private static int libraryLoadsForTest=0;
	private static final int METHOD_BLOSUM=0, METHOD_HBM=1, MODE_BLOSUM=1, MODE_HBM=2, MODE_BOTH=MODE_BLOSUM|MODE_HBM;
	private static final String[] CONSENSUS_PINNED={"in", "consensus", "sidecar", "kmers"};
	private static final String[] HBM_PINNED={"in", "hbm", "consensus", "provenance", "background", "sidecar", "kmers"};
	private static final String[] HBM_CONTROL_PINNED={"in", "hbm", "control", "consensus", "provenance", "background", "sidecar", "kmers"};
	private static final String[] BOTH_PINNED=HBM_PINNED;
	private static final String[] METHODS={"blosum", "hbm"};
	private static final String SHARD_HEADER="shard_id\tqueries\tinput\tinput_sha80\tresult";
	private static final String HEADER="query\tmethod\trank\tfamily_id\tshortlist_rank\tquery_length\tmodel_length\tstatus\tscore64\tquery_start\tmodel_start\tpaired\tmatches\tcolumns\tinsertions\tdeletions\tterminal_paired\tterminal_matches\tterminal_score64\ttrailing_gap\tpath";
	private HbmTailEvidence(){}
}
