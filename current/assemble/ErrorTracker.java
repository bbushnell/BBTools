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
	
	/**
	 * Creates a new error tracking instance with all counters initialized to zero.
	 */
	public ErrorTracker(){
		
	}
	
	/** Resets all error tracking counters and flags to their initial state.
	 * Clears both detection and correction statistics along with control flags. */
	public void clear(){
		clearDetected();
		clearCorrected();
		
		rollback=false;
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
	}

	/** Resets all error correction counters to zero.
	 * Clears pincer, tail, brute-force, and both types of reassembly correction counts. */
	public void clearCorrected() {
		correctedPincer=0;
		correctedTail=0;
		correctedBrute=0;
		correctedReassembleInner=0;
		correctedReassembleOuter=0;
	}
	
	/**
	 * Returns the total number of errors corrected across all methods.
	 * Sums correction counts from pincer, tail, brute-force, and both reassembly types.
	 * @return Total count of corrected errors
	 */
	public int corrected(){
		return correctedPincer+correctedTail+correctedBrute+correctedReassembleInner+correctedReassembleOuter; //Sum all correction counts
	}
	
	//TODO: POSSIBLE BUG - uses correctedTail instead of detectedTail
	/**
	 * Returns the total number of errors detected across methods.
	 * Sums detection counts from pincer, tail, and reassembly methods.
	 * Note: Does not include brute-force detection count.
	 * @return Total count of detected errors
	 */
	public int detected(){
		return detectedPincer+detectedTail+detectedReassemble;
	}
	
	/**
	 * Returns the total number of errors corrected via reassembly methods.
	 * Sums both inner and outer reassembly correction counts.
	 * @return Total reassembly correction count
	 */
	public int correctedReassemble(){
		return correctedReassembleInner+correctedReassembleOuter; //Sum both reassembly types
	}
	
	@Override
	public String toString(){
		ByteBuilder sb=new ByteBuilder();
		sb.append("suspected         \t").append(suspected).nl();
		sb.append("detectedPincer    \t").append(detectedPincer).nl();
		sb.append("detectedTail      \t").append(detectedTail).nl();
//		sb.append("detectedBrute     \t").append(detectedBrute).nl();
		sb.append("detectedReassemble\t").append(detectedReassemble).nl();
		sb.append("correctedPincer   \t").append(correctedPincer).nl();
		sb.append("correctedTail     \t").append(correctedTail).nl();
//		sb.append("correctedBrute    \t").append(correctedBrute).nl();
		sb.append("correctedReassembleInner\t").append(correctedReassembleInner).nl();
		sb.append("correctedReassembleOuter\t").append(correctedReassembleOuter).nl();
		sb.append("marked            \t").append(marked);
		return sb.toString();
	}

	/** Number of errors suspected but not yet confirmed */
	public int suspected;
	
	/** Number of errors detected using pincer method */
	public int detectedPincer;
	/** Number of errors detected using tail method */
	public int detectedTail;
	/** Number of errors detected using brute-force method */
	public int detectedBrute;
	/** Number of errors detected requiring reassembly */
	public int detectedReassemble;
	
	/** Number of errors corrected using pincer method */
	public int correctedPincer;
	/** Number of errors corrected using tail method */
	public int correctedTail;
	/** Number of errors corrected using brute-force method */
	public int correctedBrute;
	/** Number of errors corrected via inner reassembly method */
	public int correctedReassembleInner;
	/** Number of errors corrected via outer reassembly method */
	public int correctedReassembleOuter;
	
	/** Number of sequences marked for special processing */
	public int marked;
	
	/** Flag indicating whether rollback operations should be performed */
	public boolean rollback=false;
	
}
