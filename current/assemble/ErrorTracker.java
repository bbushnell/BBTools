package assemble;

import structures.ByteBuilder;

/**
 * Tracks and manages sequencing errors detected during assembly processing.
 * Maintains statistics on error detection, correction attempts, and rollback operations
 * for quality assessment and debugging.
 *
 * @author Brian Bushnell
 * @documentation Eru
 * @date Oct 1, 2016
 */
public class ErrorTracker {
	
	public ErrorTracker(){
		
	}
	
	/** Resets all error tracking counters and flags to their initial state.
	 * Clears both detection and correction statistics along with control flags. */
	public void clear(){
		clearDetected();
		clearCorrected();
		
		rollback=false;
		rollbackMutate=false;
		rollbackTrigger=0;
		rollbackCorrectedSnapshot=0;
		suspected=0;
		marked=0;
	}
	
	/** Resets all error detection counters to zero.
	 * Clears pincer, tail, brute-force, and reassembly detection counts. */
	public void clearDetected(){
		detectedPincer=0;
		detectedTail=0;
		detectedBrute=0;
		detectedReassemble=0;
		detectedMutate=0;
	}

	/** Resets all error correction counters to zero.
	 * Clears pincer, tail, brute-force, and both types of reassembly correction counts. */
	public void clearCorrected() {
		correctedPincer=0;
		correctedTail=0;
		correctedBrute=0;
		correctedReassembleInner=0;
		correctedReassembleOuter=0;
		correctedMutate=0;
	}

	/**
	 * Returns the total number of errors corrected across all methods.
	 * Sums correction counts from pincer, tail, brute-force, both reassembly types, and mutate.
	 * @return Total count of corrected errors
	 */
	public int corrected(){
		return correctedPincer+correctedTail+correctedBrute+correctedReassembleInner+correctedReassembleOuter+correctedMutate; //Sum all correction counts
	}
	
	//TODO: POSSIBLE BUG - uses correctedTail instead of detectedTail
	/**
	 * Returns the total number of errors detected across methods.
	 * Sums detection counts from pincer, tail, reassembly, and mutate methods.
	 * Note: Does not include brute-force detection count.
	 * @return Total count of detected errors
	 */
	public int detected(){
		return detectedPincer+detectedTail+detectedReassemble+detectedMutate;
	}
	
	/**
	 * Returns the total number of errors corrected via reassembly methods.
	 * Sums both inner and outer reassembly correction counts.
	 * @return Total reassembly correction count
	 */
	public int correctedReassemble(){
		return correctedReassembleInner+correctedReassembleOuter; //Sum both reassembly types
	}
	
	/**
	 * Creates a formatted string representation of all tracking statistics.
	 * Displays suspected, detected, and corrected error counts in tabular format.
	 * @return Tab-separated statistics summary
	 */
	@Override
	public String toString(){
		ByteBuilder sb=new ByteBuilder();
		sb.append("suspected         \t").append(suspected).nl();
		sb.append("detectedPincer    \t").append(detectedPincer).nl();
		sb.append("detectedTail      \t").append(detectedTail).nl();
//		sb.append("detectedBrute     \t").append(detectedBrute).nl();
		sb.append("detectedReassemble\t").append(detectedReassemble).nl();
		sb.append("detectedMutate    \t").append(detectedMutate).nl();
		sb.append("correctedPincer   \t").append(correctedPincer).nl();
		sb.append("correctedTail     \t").append(correctedTail).nl();
//		sb.append("correctedBrute    \t").append(correctedBrute).nl();
		sb.append("correctedReassembleInner\t").append(correctedReassembleInner).nl();
		sb.append("correctedReassembleOuter\t").append(correctedReassembleOuter).nl();
		sb.append("correctedMutate   \t").append(correctedMutate).nl();
		sb.append("marked            \t").append(marked);
		return sb.toString();
	}

	public int suspected;
	
	public int detectedPincer;
	public int detectedTail;
	public int detectedBrute;
	public int detectedReassemble;
	/** Item 3c: positions where some alternate base yields a higher-confirmed depth than the
	 * current base (evidence of an error), independent of whether the position was actually
	 * corrected -- so detected&gt;=corrected always holds. Same window basis as the confirmation
	 * counts in {@link assemble.Tadpole1#errorCorrectMutate}. */
	public int detectedMutate;
	
	public int correctedPincer;
	public int correctedTail;
	public int correctedBrute;
	public int correctedReassembleInner;
	public int correctedReassembleOuter;
	public int correctedMutate;

	public int marked;

	public boolean rollback=false;
	/** Item 3c: true if this read's rollback (if any) discarded at least one mutate correction --
	 * lets "both" mode (reassemble+mutate) distinguish mutate's rollbacks from reassemble's. */
	public boolean rollbackMutate=false;
	/** Item 3c: which rollback trigger fired for this read (0=no rollback yet/none,
	 * 1=over-correction heuristic, 2=count-contradiction check). Diagnostic only -- set inside
	 * {@link assemble.Tadpole1#errorCorrect}, read back out by the per-thread accounting in
	 * {@link assemble.Tadpole}. */
	public int rollbackTrigger=0;
	/** Item 3c: snapshot of {@link #corrected()} taken immediately before a rollback clears it --
	 * feeds the rolled-back-reads corrected-count histogram. */
	public int rollbackCorrectedSnapshot=0;

}
