package assemble;

import structures.ByteBuilder;

/**
 * Mutable per-read detection, correction and rollback counters used by the
 * assembly and Bloom correctors. Callers update fields directly and clear this
 * object before processing another read. It owns no read or sequence data.
 * Method totals sum algorithm counters rather than deduplicating base positions;
 * this object is not synchronized and belongs to one correction worker.
 *
 * @author Brian Bushnell
 * @documentation Eru
 * @date Oct 1, 2016
 */
public class ErrorTracker{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a tracker with zero counts and no rollback request. */
	public ErrorTracker(){}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Resets detections, corrections, marking and rollback diagnostics for another read. */
	public void clear(){
		clearDetected();
		clearCorrected();

		rollback=false;
		rollbackSubstitute=false;
		rollbackTrigger=0;
		rollbackCorrectedSnapshot=0;
		suspected=0;
		marked=0;
	}

	/** Clears every detection counter, including dormant brute-force and substitute counts. */
	public void clearDetected(){
		detectedPincer=0;
		detectedTail=0;
		detectedBrute=0;
		detectedReassemble=0;
		detectedSubstitute=0;
	}

	/** Clears correction counts without changing detections or saved rollback diagnostics. */
	public void clearCorrected(){
		correctedPincer=0;
		correctedTail=0;
		correctedBrute=0;
		correctedReassembleInner=0;
		correctedReassembleOuter=0;
		correctedSubstitute=0;
	}

	/** Returns the sum of all correction counters, including the brute-force counter. */
	public int corrected(){
		return correctedPincer+correctedTail+correctedBrute+correctedReassembleInner+correctedReassembleOuter+correctedSubstitute;
	}

	/**
	 * Sums pincer, tail, reassembly and substitute detections. Unlike
	 * {@link #corrected()}, this excludes the dormant brute-force counter.
	 * @return Sum of the four included algorithm counters, not unique base positions
	 */
	public int detected(){
		return detectedPincer+detectedTail+detectedReassemble+detectedSubstitute;
	}

	/** Returns the sum of inner and outer reassembly corrections. */
	public int correctedReassemble(){
		return correctedReassembleInner+correctedReassembleOuter;
	}

	/**
	 * Formats suspected, active detection/correction and marked counters as labeled
	 * lines. Brute-force counters and rollback diagnostics are omitted.
	 * @return Newly allocated summary with no trailing newline
	 */
	@Override
	public String toString(){
		ByteBuilder sb=new ByteBuilder();
		sb.append("suspected         \t").append(suspected).nl();
		sb.append("detectedPincer    \t").append(detectedPincer).nl();
		sb.append("detectedTail      \t").append(detectedTail).nl();
//		sb.append("detectedBrute     \t").append(detectedBrute).nl();
		sb.append("detectedReassemble\t").append(detectedReassemble).nl();
		sb.append("detectedSubstitute\t").append(detectedSubstitute).nl();
		sb.append("correctedPincer   \t").append(correctedPincer).nl();
		sb.append("correctedTail     \t").append(correctedTail).nl();
//		sb.append("correctedBrute    \t").append(correctedBrute).nl();
		sb.append("correctedReassembleInner\t").append(correctedReassembleInner).nl();
		sb.append("correctedReassembleOuter\t").append(correctedReassembleOuter).nl();
		sb.append("correctedSubstitute\t").append(correctedSubstitute).nl();
		sb.append("marked            \t").append(marked);
		return sb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Caller-assigned estimate of possible errors before correction. */
	public int suspected;

	/** Detections reported by pincer correction. */
	public int detectedPincer;
	/** Detections reported by tail correction. */
	public int detectedTail;
	/** Dormant brute-force detections; omitted from detected() and toString(). */
	public int detectedBrute;
	/** Detections reported by reassembly correction. */
	public int detectedReassemble;
	/** Substitute-correction position/pass observations with a sampled alternate k-mer
	 * count above the corresponding nonnegative current count. The alternate need not
	 * clear the correction threshold or be applied; later passes may recount a site.
	 * Uses the windows in {@link Tadpole1#errorCorrectSubstitute} and its Tadpole2
	 * counterpart. These public counters do not enforce a detection/correction inequality.
	 * Historical substitute bookkeeping: Item 3c. */
	public int detectedSubstitute;

	/** Corrections reported by pincer correction. */
	public int correctedPincer;
	/** Corrections reported by tail correction. */
	public int correctedTail;
	/** Dormant brute-force corrections; included in corrected() but omitted from toString(). */
	public int correctedBrute;
	/** Reassembly corrections assigned to inner positions by the caller. */
	public int correctedReassembleInner;
	/** Reassembly corrections assigned to outer positions by the caller. */
	public int correctedReassembleOuter;
	/** Corrections reported by substitute correction. */
	public int correctedSubstitute;

	/** Number of bases marked by the caller's marking pass. */
	public int marked;

	/** Caller-controlled indication that this read's corrections must be or were rolled back. */
	public boolean rollback=false;
	/** Marks a Tadpole whole-read rollback that discards substitute corrections.
	 * Preserved by clearCorrected() for accounting, then reset by clear() (Item 3c). */
	public boolean rollbackSubstitute=false;
	/** Classified rollback trigger: 1=over-correction heuristic, 2=count contradiction.
	 * Zero means no classified trigger, including an already requested rollback.
	 * Diagnostic only: set by Tadpole1/2 correction and read by Tadpole's worker
	 * accounting; clearCorrected() preserves it (historical Item 3c). */
	public int rollbackTrigger=0;
	/** Snapshot of {@link #corrected()} before a whole-read rollback clears the counts.
	 * Supplies the rolled-back-read correction histogram and survives clearCorrected();
	 * clear() resets it (historical Item 3c). */
	public int rollbackCorrectedSnapshot=0;

}
