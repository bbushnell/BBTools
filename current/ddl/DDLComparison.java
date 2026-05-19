package ddl;

import cardinality.DynamicDemiLog;
import simd.Vector;

/**
 * Mutable result record for a DDL pairwise comparison.
 * Designed for reuse: populate via {@link #set}, compare to heap bottom,
 * copy via {@link #setFrom} only when keeping.  All fields are instance-based.
 *
 * Stores references to the query and ref {@link DDLRecord} rather than
 * duplicating their fields — avoids the antipattern of maintaining the
 * same data in two places.
 *
 * @author Noire, Brian Bushnell
 * @date May 17, 2026
 */
public class DDLComparison implements Comparable<DDLComparison> {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public DDLComparison(){}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Populates this comparison from raw DDL comparison results.
	 * @param cmp {lower, equal, higher, bothEmpty} from compareToDetailed or Vector.compareDDL
	 * @param qRec Query record (may be null for anonymous queries)
	 * @param rRec Reference record
	 * @param k K-mer length */
	public void set(int[] cmp, DDLRecord qRec, DDLRecord rRec, int k){
		lower=cmp[0];
		equal=cmp[1];
		higher=cmp[2];
		bothEmpty=cmp.length>3 ? cmp[3] : 0;

		ani=DynamicDemiLog.ani(lower, equal, higher, k);
		wkid=DynamicDemiLog.wkid(lower, equal, higher);
		containmentAB=DynamicDemiLog.containmentAB(lower, equal, higher);
		containmentBA=DynamicDemiLog.containmentBA(lower, equal, higher);
		completenessAB=DynamicDemiLog.completeness(lower, equal, higher);
		completenessBA=DynamicDemiLog.completenessBA(lower, equal, higher);

		queryRecord=qRec;
		refRecord=rRec;
		cachedScore=calcScore();
	}

	/** Populates from index match counts without a full comparison.
	 * WKID, containment, and completeness will be -1 (unavailable).
	 * Use this for turbo/index-only mode. */
	public void setFromIndex(int matchCount, DDLRecord qRec, DDLRecord rRec,
			int queryFilled, int refFilled, int k){
		equal=matchCount;
		lower=-1;
		higher=-1;
		bothEmpty=-1;
		int minDiv=Math.min(queryFilled, refFilled);
		if(minDiv>0 && matchCount>0){
			double c=Math.min(1.0, (double)matchCount/minDiv);
			ani=(float)Math.exp(Math.log(c)/k);
		}else{
			ani=0;
		}
		wkid=-1;
		containmentAB=-1;
		containmentBA=-1;
		completenessAB=-1;
		completenessBA=-1;
		queryRecord=qRec;
		refRecord=rRec;
		cachedScore=calcScore();
	}

	/** Performs a full DDL comparison and populates this record.
	 * @return this, for chaining */
	public DDLComparison compare(DDLRecord qRec, DDLRecord rRec, int k){
		int[] cmp=Vector.compareDDL(qRec.ddl.maxArray(), rRec.ddl.maxArray());
		set(cmp, qRec, rRec, k);
		return this;
	}

	/** Deep copy from another DDLComparison. */
	public void setFrom(DDLComparison other){
		lower=other.lower;
		equal=other.equal;
		higher=other.higher;
		bothEmpty=other.bothEmpty;

		ani=other.ani;
		wkid=other.wkid;
		containmentAB=other.containmentAB;
		containmentBA=other.containmentBA;
		completenessAB=other.completenessAB;
		completenessBA=other.completenessBA;

		queryRecord=other.queryRecord;
		refRecord=other.refRecord;
		ssuIdentity=other.ssuIdentity;
		cachedScore=other.cachedScore;
	}

	/** Reset all fields to default (uncompared) state. */
	public void clear(){
		lower=0; equal=0; higher=0; bothEmpty=0;
		ani=-1; wkid=-1;
		containmentAB=-1; containmentBA=-1;
		completenessAB=-1; completenessBA=-1;
		ssuIdentity=-1;
		queryRecord=null; refRecord=null;
		cachedScore=0;
	}

	/** Composite score weighting similarity by evidence strength.
	 * Uses WKID (raw similarity) scaled by sqrt(matches), so a reference
	 * with many matches at slightly lower similarity outranks one with
	 * few matches at nominally higher similarity.
	 * Cached in {@link #cachedScore} — computed once in set/setFrom, compared as a raw float. */
	private float calcScore(){
		float sim=wkid>=0 ? wkid : (ani>=0 ? ani : 0);
		return sim*(float)Math.sqrt(Math.max(0, equal));
	}

	/** Higher cached score is better; better comparisons sort first.
	 * Single float comparison — one instruction, no branches. */
	@Override
	public int compareTo(DDLComparison o){
		return cachedScore>o.cachedScore ? -1 : cachedScore<o.cachedScore ? 1 : 0;
	}

	/** Returns true if this comparison has real WKID/completeness data
	 * (not just index-approximate values). */
	public boolean hasFullMetrics(){
		return wkid>=0;
	}

	@Override
	public String toString(){
		return "DDLComparison(ani="+ani+", wkid="+wkid+", equal="+equal+
			", ref="+(refRecord!=null ? refRecord.name : "null")+
			", tid="+(refRecord!=null ? refRecord.taxID : -1)+")";
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/* Raw comparison buckets */
	public int lower;
	public int equal;
	public int higher;
	public int bothEmpty;

	/* Computed metrics */
	public float cachedScore;
	public float ani=-1;
	public float wkid=-1;
	public float containmentAB=-1;
	public float containmentBA=-1;
	public float completenessAB=-1;
	public float completenessBA=-1;

	/* SSU alignment */
	public float ssuIdentity=-1;

	/* Source records — access metadata through these, don't duplicate */
	public DDLRecord queryRecord;
	public DDLRecord refRecord;
}
