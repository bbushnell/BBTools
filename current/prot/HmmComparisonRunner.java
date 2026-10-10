package prot;

import java.util.HashMap;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import aligner.SingleStateAlignerFlat2Amino;
import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Measures unchanged production assignment or existing pairwise kernels on a fixed panel. @author Keqing */
public final class HmmComparisonRunner {

	public static void main(String[] args){
		try{run(args);}catch(Throwable failure){
			failure.printStackTrace(System.err);
			// Output failures must terminate writer threads as well as the calling thread.
			System.exit(1);
		}
	}
	private static void run(String[] args) throws Exception{
		if(args.length==1 && args[0].equals("selftest")){selftest(); return;}
		final HashMap<String,String> opt=HmmComparisonData.options(args);
		final String mode=HmmComparisonData.required(opt,"mode");
		final boolean matched=parse.Parse.parseBoolean(opt.getOrDefault("matched","false"));
		if(!mode.equals("hbm") && !mode.equals("ordinary") && !mode.equals("blosum")){
			throw new IllegalArgumentException("Unknown alignment mode: "+mode);
		}
		final long start=System.nanoTime();
		final List<ProteinSequence> queries=ProteinSearch.readFasta(HmmComparisonData.required(opt,"in"));
		final String resource=HmmComparisonData.required(opt,"resources");
		if(!matched && (opt.containsKey("refs") || opt.containsKey("hbm") || opt.containsKey("provenance"))){throw new IllegalArgumentException("Artifact overrides are restricted to the frozen matched assay");}
		final String reps=opt.getOrDefault("refs",resource+"/consensus_reps_round2.faa.gz"), sets=resource+"/sets.tsv.gz";
		final String side=resource+"/magqc_sidecar_v1.tsv.gz";
		final ProteinSearcher.AssignmentBinding binding=mode.equals("hbm") && !matched ? ProteinSearcher.AssignmentBinding.loadSchema7(
			resource+"/family_thresholds_schema7_rare01_text.tsv.gz","8f941b1587695bd233b7",
			resource+"/roster_v4.tsv.gz",reps,resource+"/family_roles_v1.tsv.gz",
			resource+"/family_core_coordinates.tsv.gz",sets,side,
			resource+"/magqc_hbm_v1.rare01.hbmt.gz",resource+"/PROVENANCE_MANIFEST.tsv.gz") : null;
		final FamilyShortlistSidecar shortlist=binding==null && !matched ? FamilyShortlistSidecar.load(side,reps,sets) : null;
		final List<ProteinSequence> targets=binding==null ? ProteinSearch.readFasta(reps) : null;
		final HashMap<String,int[]> frozenPairs=new HashMap<String,int[]>();
		if(matched){
			final HashMap<String,Integer> rank=new HashMap<String,Integer>();
			for(int i=0; i<targets.size(); i++){
				if(rank.put(targets.get(i).id,i)!=null){throw new IllegalArgumentException("Duplicate matched reference: "+targets.get(i).id);}
			}
			for(String[] row : HmmComparisonData.rows(HmmComparisonData.required(opt,"pairs"))){
				if(row[0].equals("query")){continue;}
				if(row.length!=10){throw new IllegalArgumentException("Frozen shortlist row needs ten fields: "+row[0]);}
				final String[] ids=row[9].split(",",-1); final int[] indexes=new int[ids.length];
				if(ids.length!=50){throw new IllegalArgumentException("Expected full50 shortlist for "+row[0]);}
				final java.util.HashSet<String> unique=new java.util.HashSet<String>();
				for(int j=0; j<ids.length; j++){
					final Integer r=rank.get(ids[j]);
					if(r==null || !unique.add(ids[j])){throw new IllegalArgumentException("Unknown/duplicate candidate: "+ids[j]);}
					indexes[j]=r;
				}
				if(frozenPairs.put(row[0],indexes)!=null){throw new IllegalArgumentException("Duplicate shortlist query: "+row[0]);}
			}
			if(frozenPairs.size()!=queries.size()){throw new IllegalArgumentException("Frozen shortlist/query count differs");}
			for(ProteinSequence q : queries){if(!frozenPairs.containsKey(q.id)){throw new IllegalArgumentException("Missing frozen query: "+q.id);}}
		}
		final HbmBundleLoader.Loaded profile;
		if(matched && mode.equals("hbm")){
			final java.util.ArrayList<String> roster=new java.util.ArrayList<String>();
			final HashMap<String,byte[]> sequenceById=new HashMap<String,byte[]>();
			for(ProteinSequence target : targets){roster.add(target.id); sequenceById.put(target.id,target.enc);}
			profile=HbmBundleLoader.load(java.nio.file.Paths.get(opt.getOrDefault("hbm",resource+"/magqc_hbm_v1.rare01.hbmt.gz")),roster,
				id->sequenceById.get(id),HbmBundleLoader.loadSemanticProvenance(opt.getOrDefault("provenance",resource+"/PROVENANCE_MANIFEST.tsv.gz")));
		}else{profile=null;}
		if(targets!=null && shortlist!=null){
			if(targets.size()!=shortlist.nFamilies){throw new IllegalArgumentException("Consensus and shortlist counts differ");}
			for(int i=0; i<targets.size(); i++){
				if(!targets.get(i).id.equals(shortlist.repIds[i])){throw new IllegalArgumentException("Consensus order differs at "+i);}
			}
		}
		final byte[][] rawTargets=targets==null ? null : new byte[targets.size()][];
		if(rawTargets!=null){for(int i=0; i<targets.size(); i++){rawTargets[i]=HmmComparisonData.ascii(targets.get(i).enc);}}
		final int threads=Integer.parseInt(opt.getOrDefault("t","1"));
		if(threads<1){throw new IllegalArgumentException("t must be positive");}
		final Result[] results=new Result[queries.size()];
		final AtomicInteger next=new AtomicInteger();
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final Thread[] workers=new Thread[threads];
		final long loaded=System.nanoTime();
		for(int w=0; w<threads; w++){
			workers[w]=new Thread(()->{
				try{
					final ProteinSearcher searcher=new ProteinSearcher(); searcher.aligner="d55";
					final ProteinSearcher.ShortlistScratch scratch=matched ? null : (binding!=null ? binding.newScratch() :
						new ProteinSearcher.ShortlistScratch(shortlist.dims,shortlist.nFamilies,50));
					final SingleStateAlignerFlat2Amino ordinary=new SingleStateAlignerFlat2Amino();
					final ProteinSearcher.D55Metrics metrics=new ProteinSearcher.D55Metrics();
					for(int i=next.getAndIncrement(); i<queries.size() && failure.get()==null; i=next.getAndIncrement()){
						final ProteinSequence q=queries.get(i);
						if(binding!=null){
							final ProteinSearcher.FamilyAssignment a=searcher.assignFamily(binding,q,scratch,
								ProteinSearcher.AssignPolicy.BOUNDED_LOOKAHEAD,4);
							if(a.isAssigned()){finite(a.R,a.identity,a.overlap,a.hbmRaw,a.hbmPathRelative);}
							results[i]=a.isAssigned() ? new Result(a.repId,a.R,a.identity,a.overlap,
								a.rawScore,a.hbmRaw,a.hbmPathRelative,a.reason.name(),"NA") :
								new Result("NA",Double.NaN,Double.NaN,Double.NaN,Double.NaN,Double.NaN,Double.NaN,a.reason.name(),"NA");
							results[i].alignments=scratch.lastAlignedCandidates;
							results[i].tracebacks=scratch.lastTracebackAlignments;
						}else{
							if(!matched){
							ProteinSearcher.scoreFamiliesF4(q,shortlist,scratch);
							// Same unmasked F4 top50 for both pairwise arms; production HBM keeps its own gates.
							select(scratch.f4,scratch.topIdx);
							}
							final ByteBuilder candidates=new ByteBuilder();
							final byte[] rawQ=HmmComparisonData.ascii(q.enc);
							Result best=null;
							int alignments=0;
							for(int rank : matched ? frozenPairs.get(q.id) : scratch.topIdx){
								final ProteinSequence t=targets.get(rank);
								if(candidates.length()>0){candidates.append(',');} candidates.append(t.id);
								final double ratio=Math.min(q.length(),t.length())/(double)Math.max(q.length(),t.length());
								if(!matched && ratio<0.5){continue;}
								alignments++;
								final Result r;
								if(mode.equals("ordinary")){
									final double identity=ordinaryIdentity(ordinary,rawQ,rawTargets[rank]);
									finite(identity);
									r=new Result(t.id,identity,identity,ratio,Double.NaN,Double.NaN,Double.NaN,"CANDIDATE","");
								}else if(profile!=null){
									final AAAlignment a=GlocalAminoLinear.align(q.enc,t.enc,true);
									// In this experiment the user-selected shortlist is the only prefilter; production gates are deliberately absent.
									final float[] h=profile.score(rank,q.enc,a,true);
									if(h==null || h.length!=2){throw new IllegalStateException("Missing full-shortlist HBM score");}
									finite(h[0],h[1]);
									ProteinSearcher.d55Metrics(q,t,a,metrics);
									r=new Result(t.id,h[1],a.pident()/100.0,metrics.overlap,a.rawScore,h[0],h[1],"CANDIDATE","");
								}else{
									final AAAlignment a=GlocalAminoLinear.align(q.enc,t.enc);
									ProteinSearcher.d55Metrics(q,t,a,metrics);
									finite(metrics.R,metrics.overlap,a.pident());
									if(!matched && metrics.overlap<0.5){continue;}
									r=new Result(t.id,metrics.R,a.pident()/100.0,metrics.overlap,a.rawScore,Double.NaN,Double.NaN,"CANDIDATE","");
								}
								if(best==null || r.score>best.score || (r.score==best.score && r.family.compareTo(best.family)<0)){best=r;}
							}
							if(best==null){best=new Result("NA",Double.NEGATIVE_INFINITY,0,0,Double.NaN,Double.NaN,Double.NaN,"NO_CANDIDATE","");}
							best.candidates=candidates.toString(); results[i]=best;
							best.alignments=alignments;
						}
						if(i%10000==0){System.err.println(mode+" processed ordinal "+i+" / "+queries.size());}
					}
				}catch(Throwable ex){failure.compareAndSet(null,ex);}
			},"hmm-comparison-"+w);
			workers[w].start();
		}
		for(Thread worker : workers){worker.join();}
		if(failure.get()!=null){throw new RuntimeException("Comparison worker failed; no results published",failure.get());}
		final long scored=System.nanoTime();
		final ByteStreamWriter out=HmmComparisonData.writer(HmmComparisonData.required(opt,"out"));
		out.println("query\tfamily\tscore\tidentity\t"+(mode.equals("ordinary") ? "length_ratio" : "aligned_overlap")+
			"\traw\thbm_raw\thbm_relative\treason\tcandidates");
		for(int i=0; i<queries.size(); i++){
			final Result r=results[i];
			if(r==null){throw new IllegalStateException("Missing result: "+i);}
			final ByteBuilder row=new ByteBuilder().append(queries.get(i).id).tab().append(r.family);
			for(double value : new double[]{r.score,r.identity,r.overlap,r.raw,r.hbmRaw,r.hbmRelative}){
				row.tab(); if(Double.isFinite(value)){row.append(value,8);}else{row.append("NA");}
			}
			out.println(row.tab().append(r.reason).tab().append(r.candidates));
		}
		HmmComparisonData.close(out);
		final ByteStreamWriter timing=HmmComparisonData.writer(HmmComparisonData.required(opt,"timing"));
		timing.println("mode\tqueries\tthreads\tload_s\tcompute_s\ttotal_s");
		timing.println(new ByteBuilder().append(mode).tab().append(queries.size()).tab().append(threads).tab()
			.append((loaded-start)/1e9,6).tab().append((scored-loaded)/1e9,6).tab().append((System.nanoTime()-start)/1e9,6));
		HmmComparisonData.close(timing);
		if(opt.containsKey("work")){
			long alignments=0,tracebacks=0;
			for(Result r : results){alignments+=r.alignments; tracebacks+=r.tracebacks;}
			final ByteStreamWriter work=HmmComparisonData.writer(opt.get("work"));
			// Production's lastTracebackAlignments counts ADDITIONAL path-recorded DP calls, not the initial alignment's walk.
			work.println("mode\tqueries\tcandidate_alignments\tadditional_path_alignments\ttotal_dp_alignments");
			work.println(new ByteBuilder().append(mode).tab().append(queries.size()).tab().append(alignments).tab().append(tracebacks).tab().append(alignments+tracebacks));
			HmmComparisonData.close(work);
		}
	}

	/** Stable descending top-k; lower family rank wins exact F4 ties, matching production ordering. */
	private static void select(float[] scores,int[] top){
		assert(scores.length>=top.length && top.length>0) : "F4 shortlist cannot exceed its family array";
		java.util.Arrays.fill(top,-1);
		for(int i=0; i<scores.length; i++){
			for(int k=0; k<top.length; k++){
				if(top[k]<0 || scores[i]>scores[top[k]]){
					System.arraycopy(top,k,top,k+1,top.length-k-1); top[k]=i; break;
				}
			}
		}
	}

	/** Finite oracles for shortlist ties, residue round-trip and exact-match kernel contracts. */
	private static void selftest(){
		final int[] top=new int[3]; select(new float[]{1,4,4,-1,2},top);
		if(!java.util.Arrays.equals(top,new int[]{1,2,4})){throw new AssertionError("F4 tie/order oracle failed");}
		final String amino="ARNDCQEGHILKMFPSTWYV";
		final ProteinSequence q=new ProteinSequence("fixture",amino);
		if(!amino.equals(new String(HmmComparisonData.ascii(q.enc),java.nio.charset.StandardCharsets.US_ASCII))){
			throw new AssertionError("Amino-acid encoding round-trip failed");
		}
		final float identity=new SingleStateAlignerFlat2Amino().align(HmmComparisonData.ascii(q.enc),HmmComparisonData.ascii(q.enc));
		if(identity!=1){throw new AssertionError("Ordinary exact-match identity must equal one: "+identity);}
		final byte[] shortSeq="ARNDCQEGHILKMFPSTWYV".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
		final byte[] longSeq=("WWWWWWWWWWWWWWWWWWWW"+amino).getBytes(java.nio.charset.StandardCharsets.US_ASCII);
		if(ordinaryIdentity(new SingleStateAlignerFlat2Amino(),longSeq,shortSeq)!=1){
			throw new AssertionError("Unequal-length suffix match was truncated by implicit query/reference swap");
		}
		final ProteinSearcher.D55Metrics metrics=new ProteinSearcher.D55Metrics();
		ProteinSearcher.d55Metrics(q,q,GlocalAminoLinear.align(q.enc,q.enc),metrics);
		if(metrics.R!=1 || metrics.overlap!=1){throw new AssertionError("BLOSUM exact-match R/coverage must equal one");}
		System.err.println("HMM_COMPARISON_SELFTEST_PASS");
	}
	private static double ordinaryIdentity(SingleStateAlignerFlat2Amino aligner,byte[] q,byte[] t){
		//TODO: Probable engine bug - SSA2Amino.align swaps a longer query after computing to from the old reference.
		// Passing the shorter sequence first preserves its intended full reference bounds without altering the engine.
		return q.length<=t.length ? aligner.align(q,t) : aligner.align(t,q);
	}
	private static void finite(double... values){
		for(double value : values){if(!Double.isFinite(value)){throw new ArithmeticException("Nonfinite measured alignment score: "+value);}}
	}

	private static final class Result {
		Result(String family_,double score_,double identity_,double overlap_,double raw_,double hbmRaw_,
				double hbmRelative_,String reason_,String candidates_){
			family=family_; score=score_; identity=identity_; overlap=overlap_; raw=raw_;
			hbmRaw=hbmRaw_; hbmRelative=hbmRelative_; reason=reason_; candidates=candidates_;
		}
		final String family,reason;
		final double score,identity,overlap,raw,hbmRaw,hbmRelative;
		String candidates;
		int alignments,tracebacks;
	}
}
