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
	
	public void clear(){
		clearDetected();
		clearCorrected();
		
		rollback=false;
		suspected=0;
		marked=0;
	}
	
	public void clearDetected(){
		detectedPincer=0;
		detectedTail=0;
		detectedBrute=0;
		detectedReassemble=0;
	}

	public void clearCorrected() {
		correctedPincer=0;
		correctedTail=0;
		correctedBrute=0;
		correctedReassembleInner=0;
		correctedReassembleOuter=0;
	}
	
	public int corrected(){
		return correctedPincer+correctedTail+correctedBrute+correctedReassembleInner+correctedReassembleOuter; //Sum all correction counts
	}
	
	//TODO: POSSIBLE BUG - uses correctedTail instead of detectedTail
	public int detected(){
		return detectedPincer+detectedTail+detectedReassemble;
	}
	
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

	public int suspected;
	
	public int detectedPincer;
	public int detectedTail;
	public int detectedBrute;
	public int detectedReassemble;
	
	public int correctedPincer;
	public int correctedTail;
	public int correctedBrute;
	public int correctedReassembleInner;
	public int correctedReassembleOuter;
	
	public int marked;
	
	public boolean rollback=false;
	
}
