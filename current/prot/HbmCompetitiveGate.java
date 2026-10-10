package prot;

import java.util.Arrays;
import java.util.HashSet;

/**
 * Experimental highest-score eligible selection on a caller-supplied shortlist.
 * Positional scores determine ordering; whole-path identity and paired coverage
 * inside the family core determine eligibility. No production score cutoff is reused.
 * Instances are immutable; each calling thread owns its Scratch and query bytes.
 * @author Keqing
 */
public final class HbmCompetitiveGate {

	/** Binds detached references and cores to immutable positional profiles. */
	public HbmCompetitiveGate(String[] ids, byte[][] refs, HbmPositionModel[] profiles, int[] starts, int[] ends){
		require(ids!=null && refs!=null && profiles!=null && starts!=null && ends!=null, "Null competitive library component");
		final int n=ids.length;
		require(n>0 && refs.length==n && profiles.length==n && starts.length==n && ends.length==n, "Competitive library dimensions differ");
		repIds=ids.clone(); models=profiles.clone(); coreStart=starts.clone(); coreEnd=ends.clone(); references=new byte[n][];
		final HashSet<String> seen=new HashSet<String>();
		for(int i=0; i<n; i++){
			require(ids[i]!=null && !ids[i].isEmpty() && seen.add(ids[i]), "Blank or duplicate competitive family ID at "+i);
			for(int j=0; j<ids[i].length(); j++){require(ids[i].charAt(j)>32 && ids[i].charAt(j)<127, "Competitive IDs must be printable ASCII tokens for deterministic ties");}
			require(refs[i]!=null && profiles[i]!=null && starts[i]>=0 && ends[i]>=starts[i] && ends[i]<refs[i].length, "Invalid competitive core at "+i);
			references[i]=refs[i].clone(); models[i].requireConsensus(references[i]);
		}
	}

	/**
	 * Returns a family index, or -1 when no shortlisted family passes both gates.
	 * Every length-feasible candidate is scored. Recorded paths are evaluated in
	 * score/ASCII-ID order until the first eligible candidate; this is exactly
	 * the highest-score eligible result, not acceptance of an ineligible winner.
	 * Result fields and counters in work are reset on every call. The shortlist
	 * must contain distinct valid indexes; its construction is the caller's job.
	 */
	public int select(byte[] query, int[] shortlist, int count, Scratch work){
		require(query!=null && query.length>0 && shortlist!=null && work!=null, "Competitive selection requires query, shortlist and scratch");
		require(count>=0 && count<=shortlist.length && count<=work.order.length && work.seen.length==models.length,
			"Competitive shortlist/scratch dimensions differ");
		work.reset();
		int feasible=0;
		for(int k=0; k<count; k++){
			final int rank=shortlist[k];
			if(rank<0 || rank>=models.length || work.seen[rank]==work.generation){throw new IllegalArgumentException("Duplicate or invalid competitive family index: "+rank);}
			work.seen[rank]=work.generation;
			final int coreLength=coreEnd[rank]-coreStart[rank]+1;
			if(!lengthCanPass(query.length, coreLength)){work.lengthRejected++; continue;}
			final HbmPositionModel.Result result=models[rank].align(query, false); work.scored++;
			int at=feasible;
			while(at>0 && better(result.score, repIds[rank], work.scores[at-1], repIds[work.order[at-1]])){
				work.scores[at]=work.scores[at-1]; work.order[at]=work.order[at-1]; at--;
			}
			work.scores[at]=result.score; work.order[at]=rank; feasible++;
		}
		for(int k=0; k<feasible; k++){
			final int rank=work.order[k];
			final HbmPositionModel.Result result=models[rank].align(query, true); work.recorded++;
			require(result.score==work.scores[k], "Recording a positional path changed its optimum score");
			metrics(query, references[rank], result, coreStart[rank], coreEnd[rank], work.metrics);
			if(passes(work.metrics)){
				work.winner=rank; work.score64=result.score; work.start=result.start; work.end=result.end;
				return rank;
			}
		}
		work.metrics.clear();
		return -1;
	}

	/** Paired columns cannot exceed either query length or core length, even with gaps. */
	static boolean lengthCanPass(int queryLength, int coreLength){
		require(queryLength>0 && coreLength>0, "Positive sequence lengths are required for coverage");
		return (float)(Math.min(queryLength, coreLength)/(double)Math.max(queryLength, coreLength))>=MIN_COVERAGE;
	}

	/**
	 * Whole-path identity includes gap columns and excludes X identities, matching
	 * AAAlignment.pident. Coverage follows ProteinSearcher.constructionCoreMetrics:
	 * paired m columns inside the core, divided by each length or their maximum.
	 * Both metrics are rounded to float32 before the inclusive construction gates.
	 */
	static void metrics(byte[] query, byte[] ref, HbmPositionModel.Result alignment, int first, int last, Metrics out){
		require(query!=null && ref!=null && query.length>0 && alignment!=null && alignment.path!=null && alignment.path.length>0 && out!=null,
			"Competitive metrics require nonempty sequences and a recorded full-query path");
		require(first>=0 && last>=first && last<ref.length && alignment.start>=0 && alignment.end>=alignment.start && alignment.end<ref.length,
			"Competitive path/core coordinates exceed the reference");
		int q=0, r=alignment.start, identities=0, paired=0;
		for(byte op : alignment.path){
			if(op=='m'){
				require(q<query.length && r<ref.length, "Competitive paired column exceeds sequence bounds");
				if(query[q]==ref[r] && query[q]>=0 && query[q]<20){identities++;}
				if(r>=first && r<=last){paired++;}
				q++; r++;
			}else if(op=='I'){q++;}
			else if(op=='D'){r++;}
			else{throw new IllegalArgumentException("Unknown competitive path operation: "+(char)op);}
			require(q<=query.length && r<=ref.length, "Competitive path exceeds sequence bounds");
		}
		final int coreLength=last-first+1;
		require(q==query.length && r==alignment.end+1 && paired<=query.length && paired<=coreLength,
			"Competitive path must consume the full query and declared reference span exactly");
		out.identity=(float)(100.0*identities/alignment.path.length);
		out.coverageQ=(float)(paired/(double)query.length); out.coverageT=(float)(paired/(double)coreLength);
		out.coverage=(float)(paired/(double)Math.max(query.length, coreLength)); out.paired=paired;
	}

	static boolean passes(Metrics m){return m.identity>=MIN_IDENTITY && m.coverage>=MIN_COVERAGE;}
	static boolean better(int score, String id, int oldScore, String oldId){return score>oldScore || (score==oldScore && id.compareTo(oldId)<0);}
	private static void require(boolean ok, String message){if(!ok){throw new IllegalArgumentException(message);}}

	/** Reusable metrics; valid as winner evidence only when select returns a nonnegative index. */
	public static final class Metrics{
		public float identity, coverageQ, coverageT, coverage;
		public int paired;
		void clear(){identity=coverageQ=coverageT=coverage=Float.NaN; paired=0;}
	}

	/** One per worker; select overwrites prior results, arrays and counts. No per-candidate sorting objects. */
	public static final class Scratch{
		public Scratch(int families, int shortlistCapacity){
			require(families>0 && shortlistCapacity>0 && shortlistCapacity<=families, "Invalid competitive scratch capacity");
			seen=new int[families]; order=new int[shortlistCapacity]; scores=new int[shortlistCapacity];
		}
		private void reset(){
			if(generation==Integer.MAX_VALUE){Arrays.fill(seen, 0); generation=0;} generation++;
			winner=-1; score64=Integer.MIN_VALUE; start=end=-1; scored=recorded=lengthRejected=0; metrics.clear();
		}
		public int winner=-1, score64=Integer.MIN_VALUE, start=-1, end=-1, scored, recorded, lengthRejected;
		public final Metrics metrics=new Metrics();
		private final int[] seen, order, scores;
		private int generation;
	}

	public static final float MIN_IDENTITY=40.612846f, MIN_COVERAGE=0.8f;
	private final String[] repIds;
	private final byte[][] references;
	private final HbmPositionModel[] models;
	private final int[] coreStart, coreEnd;
}
