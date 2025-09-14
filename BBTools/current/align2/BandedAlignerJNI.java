package align2;

import dna.AminoAcid;
import shared.KillSwitch;
import shared.Shared;

/**
 * @author Jonathan Rood
 * @date Jul 18, 2014
 *
 */
public class BandedAlignerJNI extends BandedAligner{

	static {
		Shared.loadJNI();
	}

	private native int alignForwardJNI(byte[] query, byte[] ref, int qstart, int rstart, int maxEdits, boolean exact, int maxWidth, byte[] baseToNumber, int[] returnVals);

	private native int alignForwardRCJNI(byte[] query, byte[] ref, int qstart, int rstart, int maxEdits, boolean exact, int maxWidth, byte[] baseToNumber, byte[] baseToComplementExtended, int[] returnVals);

	private native int alignReverseJNI(byte[] query, byte[] ref, int qstart, int rstart, int maxEdits, boolean exact, int maxWidth, byte[] baseToNumber, int[] returnVals);

	private native int alignReverseRCJNI(byte[] query, byte[] ref, int qstart, int rstart, int maxEdits, boolean exact, int maxWidth, byte[] baseToNumber, byte[] baseToComplementExtended, int[] returnVals);
	
	public static void main(String[] args){
		byte[] query=args[0].getBytes();
		byte[] ref=args[1].getBytes();
		int qstart=-1;
		int rstart=-1;
		int maxedits=big;
		int width=5;
		if(args.length>2){qstart=Integer.parseInt(args[2]);}
		if(args.length>3){rstart=Integer.parseInt(args[3]);}
		if(args.length>4){maxedits=Integer.parseInt(args[4]);}
		if(args.length>4){width=Integer.parseInt(args[5]);}
		
		BandedAlignerJNI ba=new BandedAlignerJNI(width);
		
		int edits;
		
		edits=ba.alignForward(query, ref, (qstart==-1 ? 0 : qstart), (rstart==-1 ? 0 : rstart), maxedits, true);
		System.out.println("Forward:    \tedits="+edits+", lastRow="+ba.lastRow+", score="+ba.score());
		System.out.println("***********************\n");
//
//		edits=ba.alignForwardRC(query, ref, (qstart==-1 ? query.length-1 : qstart), (rstart==-1 ? 0 : rstart), maxedits, true);
//		System.out.println("ForwardRC:  \tedits="+edits+", lastRow="+ba.lastRow+", score="+ba.score());
//		System.out.println("***********************\n");
		
		edits=ba.alignReverse(query, ref, (qstart==-1 ? query.length-1 : qstart), (rstart==-1 ? ref.length-1 : rstart), maxedits, true);
		System.out.println("Reverse:    \tedits="+edits+", lastRow="+ba.lastRow+", score="+ba.score());
		System.out.println("***********************\n");
		
//		edits=ba.alignReverseRC(query, ref, (qstart==-1 ? 0 : qstart), (rstart==-1 ? ref.length-1 : rstart), maxedits, true);
//		System.out.println("ReverseRC:  \tedits="+edits+", lastRow="+ba.lastRow+", score="+ba.score());
//		System.out.println("***********************\n");
	}
	
	public BandedAlignerJNI(int width_){
		super(width_);
		assert(big>maxWidth/2);
	}
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	@Override
	public int alignForward(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact){
		int[] returnVals = KillSwitch.allocInt1D(5);
		returnVals[0] = lastQueryLoc;
		returnVals[1] = lastRefLoc;
		returnVals[2] = lastRow;
		returnVals[3] = lastEdits;
		returnVals[4] = lastOffset;
		int edits = alignForwardJNI(query,ref,qstart,rstart,maxEdits,exact,maxWidth,AminoAcid.baseToNumber,returnVals);
		lastQueryLoc = returnVals[0];
		lastRefLoc = returnVals[1];
		lastRow = returnVals[2];
		lastEdits = returnVals[3];
		lastOffset = returnVals[4];
		return edits;
	}
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	@Override
	public int alignForwardRC(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact){
		int[] returnVals = KillSwitch.allocInt1D(5);
		returnVals[0] = lastQueryLoc;
		returnVals[1] = lastRefLoc;
		returnVals[2] = lastRow;
		returnVals[3] = lastEdits;
		returnVals[4] = lastOffset;
		int edits = alignForwardRCJNI(query,ref,qstart,rstart,maxEdits,exact,maxWidth,AminoAcid.baseToNumber,AminoAcid.baseToComplementExtended,returnVals);
		lastQueryLoc = returnVals[0];
		lastRefLoc = returnVals[1];
		lastRow = returnVals[2];
		lastEdits = returnVals[3];
		lastOffset = returnVals[4];
		return edits;
	}
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	@Override
	public int alignReverse(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact){
		int[] returnVals = KillSwitch.allocInt1D(5);
		returnVals[0] = lastQueryLoc;
		returnVals[1] = lastRefLoc;
		returnVals[2] = lastRow;
		returnVals[3] = lastEdits;
		returnVals[4] = lastOffset;
		int edits = alignReverseJNI(query,ref,qstart,rstart,maxEdits,exact,maxWidth,AminoAcid.baseToNumber,returnVals);
		lastQueryLoc = returnVals[0];
		lastRefLoc = returnVals[1];
		lastRow = returnVals[2];
		lastEdits = returnVals[3];
		lastOffset = returnVals[4];
		return edits;
	}
	
	/**
	 * @param query
	 * @param ref
	 * @param qstart
	 * @param rstart
	 * @return Edit distance
	 */
	@Override
	public int alignReverseRC(final byte[] query, final byte[] ref, final int qstart, final int rstart, final int maxEdits, final boolean exact){
		int[] returnVals = KillSwitch.allocInt1D(5);
		returnVals[0] = lastQueryLoc;
		returnVals[1] = lastRefLoc;
		returnVals[2] = lastRow;
		returnVals[3] = lastEdits;
		returnVals[4] = lastOffset;
		int edits = alignReverseRCJNI(query,ref,qstart,rstart,maxEdits,exact,maxWidth,AminoAcid.baseToNumber,AminoAcid.baseToComplementExtended,returnVals);
		lastQueryLoc = returnVals[0];
		lastRefLoc = returnVals[1];
		lastRow = returnVals[2];
		lastEdits = returnVals[3];
		lastOffset = returnVals[4];
		return edits;
	}
}
