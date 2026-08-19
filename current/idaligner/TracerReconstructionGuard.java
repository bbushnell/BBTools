package idaligner;

import java.util.*;

/**
 * GUARDS: the idaligner.Tracer traceback-reconstruction bug fixed Aug 19 2026 (the
 * blind score-derivative traceback overload -- disabled, not deleted, in Tracer.java
 * -- was provably non-injective under uniform +-1 scoring: a diagonal substitution
 * and an up insertion could produce the identical packed cell value whenever their
 * predecessor cells' position/deletion-count fields coincided, making the emitted
 * matchString silently wrong -- e.g. claiming a real 'C'/'A' mismatch was a match).
 * All seven idaligner aligners with a traceback path (ScrabbleAligner, BandedAligner,
 * DriftingAligner, GlocalAligner, ScrabbleAligner2, WobbleAligner, GlocalPlusAligner6)
 * were switched to Tracer's sequence-aware overload. This harness re-proves, on every
 * push, that the reconstruction invariant holds for all seven: for a synthetic
 * (query, ref) pair with KNOWN true bytes, walking stats.matchString against the real
 * query/ref must be self-consistent -- 'm' only where query[qi]==ref[r], 'S' only
 * where they differ, total query bases consumed == query.length, and the final
 * reference position consumed == stats.rStop. This is the ONLY invariant that
 * actually catches the bug (identity/score/rStart/rStop alone did not -- see coverage
 * note at the bottom).
 *
 * DESIGN: currently public-API-only (align...Static(byte[],byte[],AlignmentStats),
 * and stats.matchString/rStart/rStop/identity/score/matches/subs/ins/dels/ns, all
 * public) even though it now lives in the idaligner package itself (relocated here
 * from a standalone bbtools-dev harness, Aug 19 2026, per Brian's ruling: bbtools-dev
 * is for wrong-LANGUAGE artifacts, not a general test repo -- a Java test belongs in
 * BBTools). Kept public-API-only rather than reworked to use package-private access,
 * since nothing here currently needs it; being in-package now means a future check
 * COULD reach into Tracer's package-private internals (getTraceScore, header/trace
 * format) if one is ever needed, without another relocation.
 *
 * TODO (Noire, Aug 19 2026 review): test data here is pure uppercase ACGT only -- no
 * N, no lowercase/soft-masked bases, no other IUPAC codes. This guard does NOT exercise
 * the GlocalPlusAligner6 wire-time gate flagged in #18 (its fill encodes to long[] via
 * Factory.encodeLong while the sequence-aware traceback reads the raw byte[] directly --
 * possible desync on lowercase/masked input). When #18 closes that gate, add lowercase/
 * soft-masked/N-bearing cases HERE and this same harness will guard both concerns.
 *
 * STRICTNESS NOTE: the r-1==rStop check assumes matchString never ends on a trailing
 * 'I' past rStop (i.e. the last op is ref-consuming). Empirically true for all 7
 * aligners across all 1946 cases (0 failures) -- documented here as an assumption, not
 * a proven invariant, in case a future aligner's fill legitimately violates it.
 *
 * USAGE: java -cp <BBTOOLS_CURRENT_DIR> idaligner.TracerReconstructionGuard
 * Exits 1 with all failures printed to stderr if ANY case fails for ANY aligner.
 * Exits 0 and prints "ALL PASS" only if every case passes for every aligner.
 *
 * @author Neptune
 * @date August 19, 2026
 */
public class TracerReconstructionGuard {

	interface AlignFn { float run(byte[] query, byte[] ref, AlignmentStats stats); }

	static Map<String, AlignFn> ALIGNERS = new LinkedHashMap<>();
	static {
		ALIGNERS.put("ScrabbleAligner", (q, r, s) -> ScrabbleAligner.alignAndTraceStatic(q, r, s));
		ALIGNERS.put("BandedAligner", (q, r, s) -> BandedAligner.alignAndTraceStatic(q, r, s));
		ALIGNERS.put("DriftingAligner", (q, r, s) -> DriftingAligner.alignAndTraceStatic(q, r, s));
		ALIGNERS.put("GlocalAligner", (q, r, s) -> GlocalAligner.alignAndTraceStatic(q, r, s));
		ALIGNERS.put("ScrabbleAligner2", (q, r, s) -> ScrabbleAligner2.alignAndTraceStatic(q, r, s));
		ALIGNERS.put("WobbleAligner", (q, r, s) -> WobbleAligner.alignAndTraceStatic(q, r, s));
		ALIGNERS.put("GlocalPlusAligner6", (q, r, s) -> GlocalPlusAligner6.alignAndTraceStatic(q, r, s));
	}

	static byte[] randomSeq(int len, long seed){
		byte[] bases={'A','C','G','T'};
		byte[] out=new byte[len];
		long state=seed;
		for(int i=0;i<len;i++){
			state=state*6364136223846793005L+1442695040888963407L;
			int idx=(int)((state>>>33)&3);
			out[i]=bases[idx];
		}
		return out;
	}
	static byte flip(byte b, long seed){
		byte[] opts={'A','C','G','T'};
		byte pick; long s=seed;
		do{s=s*6364136223846793005L+1442695040888963407L; pick=opts[(int)((s>>>33)&3)];}while(pick==b);
		return pick;
	}

	/** The load-bearing invariant. Returns null if consistent, else a description of
	 * the first problem found. Also cross-checks the structural op-count identities
	 * (matches+subs+ins+ns==qLen, matches+subs+dels+ns==refAlnLength) against the
	 * stats fields setFromMatchString already computed, folding in what would
	 * otherwise be a separate boundary-consistency check. */
	static String checkReconstruction(AlignmentStats stats, byte[] query, byte[] ref){
		if(stats.matchString==null){return "matchString is null";}
		int qi=0, r=stats.rStart;
		for(byte op : stats.matchString){
			switch(op){
				case 'm':
					if(r<0||r>=ref.length||qi<0||qi>=query.length){return "'m' out of bounds at qi="+qi+" r="+r;}
					if(query[qi]!=ref[r]){return "'m' claimed but query["+qi+"]="+(char)query[qi]+" != ref["+r+"]="+(char)ref[r];}
					qi++; r++; break;
				case 'S':
					if(r<0||r>=ref.length||qi<0||qi>=query.length){return "'S' out of bounds at qi="+qi+" r="+r;}
					if(query[qi]==ref[r]){return "'S' claimed but query["+qi+"]==ref["+r+"]=="+(char)query[qi]+" (should be 'm')";}
					qi++; r++; break;
				case 'N': qi++; r++; break;
				case 'D': if(r<0||r>=ref.length){return "'D' out of bounds at r="+r;} r++; break;
				case 'I': if(qi<0||qi>=query.length){return "'I' out of bounds at qi="+qi;} qi++; break;
				default: return "unknown op '"+(char)op+"'";
			}
		}
		if(qi!=query.length){return "consumed "+qi+" query bases, expected "+query.length;}
		if(r-1!=stats.rStop){return "final ref position "+(r-1)+" != stats.rStop="+stats.rStop;}

		// Fold in the structural op-count / boundary identities (Noire's request to fold
		// in the rStart/rStop+counts invariant): these must ALWAYS hold for any valid
		// walk, and cross-check setFromMatchString's tally against the walk above.
		int refAlnLength=stats.rStop-stats.rStart+1;
		if(stats.matches+stats.subs+stats.ins+stats.ns!=query.length){
			return "structural identity broken: matches+subs+ins+ns="+(stats.matches+stats.subs+stats.ins+stats.ns)+" != qLen="+query.length;
		}
		if(stats.matches+stats.subs+stats.dels+stats.ns!=refAlnLength){
			return "structural identity broken: matches+subs+dels+ns="+(stats.matches+stats.subs+stats.dels+stats.ns)+" != refAlnLength="+refAlnLength;
		}
		return null;
	}

	static int failures=0;
	static int cases=0;
	// Per-aligner failure detail, for the coverage summary at the end.
	static Map<String, Integer> failuresByAligner=new LinkedHashMap<>();

	static void runCase(String aligner, AlignFn fn, String caseName, byte[] query, byte[] ref){
		cases++;
		AlignmentStats stats=new AlignmentStats(true);
		stats.doTrace=true;
		fn.run(query, ref, stats);
		String problem=checkReconstruction(stats, query, ref);
		if(problem!=null){
			failures++;
			failuresByAligner.merge(aligner, 1, Integer::sum);
			System.err.println("FAIL  aligner="+aligner+"  case="+caseName+"  -- "+problem+
				"  matchString="+(stats.matchString==null?"null":new String(stats.matchString))+
				"  rStart="+stats.rStart+" rStop="+stats.rStop);
		}
	}

	public static void main(String[] args){
		int len=76; // tRNA consensus-model scale, matches the diagnostics that found the bug
		double[] subRates={0.0,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50};
		int trialsPerRate=30;

		for(Map.Entry<String, AlignFn> e : ALIGNERS.entrySet()){
			String aligner=e.getKey();
			AlignFn fn=e.getValue();

			// Substitution sweep
			for(double subRate : subRates){
				for(int t=0; t<trialsPerRate; t++){
					long seed=100000L+t*97+(long)(subRate*1000);
					byte[] ref=randomSeq(len, seed);
					byte[] query=ref.clone();
					int nSubs=(int)(subRate*len);
					long mutSeed=seed+555;
					Random posRnd=new Random(seed*31+7);
					Set<Integer> used=new HashSet<>();
					for(int k=0;k<nSubs;k++){
						int p;
						do{p=posRnd.nextInt(len);}while(used.contains(p));
						used.add(p);
						query[p]=flip(query[p], mutSeed+k*13);
					}
					runCase(aligner, fn, "subRate="+subRate+" trial="+t, query, ref);
				}
			}

			// Indel cases (the exact shapes that first exposed the bug)
			byte[] ref3=randomSeq(41, 3001);
			byte[] shortDelQuery=new byte[ref3.length-3];
			System.arraycopy(ref3, 0, shortDelQuery, 0, 10);
			System.arraycopy(ref3, 13, shortDelQuery, 10, ref3.length-13);
			runCase(aligner, fn, "short deletion (3bp)", shortDelQuery, ref3);

			byte[] insBases=randomSeq(3, 30011);
			byte[] insQuery=new byte[ref3.length+3];
			System.arraycopy(ref3, 0, insQuery, 0, 10);
			System.arraycopy(insBases, 0, insQuery, 10, 3);
			System.arraycopy(ref3, 10, insQuery, 13, ref3.length-10);
			runCase(aligner, fn, "short insertion (3bp)", insQuery, ref3);

			byte[] longRef=randomSeq(300, 4001);
			byte[] longDelQuery=new byte[longRef.length-150];
			System.arraycopy(longRef, 0, longDelQuery, 0, 75);
			System.arraycopy(longRef, 225, longDelQuery, 75, longRef.length-225);
			runCase(aligner, fn, "long deletion (150bp of 300bp)", longDelQuery, longRef);

			for(int seed=5001; seed<5006; seed++){
				byte[] a=randomSeq(60, seed);
				byte[] b=randomSeq(80, seed+1);
				runCase(aligner, fn, "unrelated random pair seed="+seed, a, b);
			}
		}

		System.out.println();
		System.out.println("=== COVERAGE: "+ALIGNERS.size()+" aligners guarded: "+ALIGNERS.keySet()+" ===");
		System.out.println("=== "+cases+" cases run, "+failures+" failures ===");
		if(failures>0){
			for(Map.Entry<String,Integer> e : failuresByAligner.entrySet()){
				System.out.println("  "+e.getKey()+": "+e.getValue()+" failures");
			}
			System.out.println("REGRESSION DETECTED -- the Tracer reconstruction invariant is broken again.");
			System.exit(1);
		}else{
			System.out.println("ALL PASS");
		}
	}
}
