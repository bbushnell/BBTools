package align2;

import java.util.Arrays;

/**
 * Standalone Needleman-Wunsch global alignment demonstration.
 * Scores exact byte matches +1, mismatches -1, and each gap base -1. Scoring
 * is fixed, not configurable. Both the whole query and the inclusive reference
 * window are consumed; ties prefer diagonal, then left, then up.
 *
 * The mutable matrix belongs to one caller. traceback returns the aligned query
 * (reference-consuming deletions are '-'), not a CIGAR or a two-sequence alignment.
 * A 2026-06-20 review reported no production callers; this review does not claim
 * a new package-wide caller audit or a BBMapS mapping-performance effect.
 *
 * @author Brian Bushnell, Collei
 * @date 2013
 */
public class NeedlemanWunsch{

	/** Demonstrates a full-reference alignment and prints matrices plus aligned query. */
	public static void main(String[] args){
		if(args.length!=2){throw new IllegalArgumentException("Expected query and reference sequence arguments");}
		final byte[] read=args[0].getBytes(), ref=args[1].getBytes();
		final NeedlemanWunsch nw=new NeedlemanWunsch(read.length, ref.length);
		nw.fill(read, ref, 0, ref.length-1);
		for(int row=0; row<nw.scores.length; row++){
			System.err.println(Arrays.toString(nw.scores[row]));
			System.err.println(Arrays.toString(nw.pointers[row]));
			System.err.println();
		}
		System.err.println(new String(nw.traceback(read, ref, 0, ref.length-1)));
	}

	/** Allocates capacity for maxRows_ query bases and maxColumns_ reference bases.
	 * The extra boundary row/column is allocated internally. Zero lengths are valid. */
	public NeedlemanWunsch(int maxRows_, int maxColumns_){
		assert(maxRows_>=0 && maxRows_<Integer.MAX_VALUE && maxColumns_>=0 && maxColumns_<Integer.MAX_VALUE) :
				"Matrix lengths must permit one boundary cell: "+maxRows_+" x "+maxColumns_;
		maxRows=maxRows_;
		maxColumns=maxColumns_;
		scores=new int[maxRows+1][maxColumns+1];
		pointers=new byte[maxRows+1][maxColumns+1];
		for(int i=0; i<maxColumns+1; i++){
			scores[0][i]=0-i;
			pointers[0][i]=LEFT;
		}
		for(int i=0; i<maxRows+1; i++){
			scores[i][0]=0-i;
			pointers[i][0]=UP;
		}
	}

	/** Fills scores/pointers for the whole query and an inclusive reference window.
	 * Empty windows use end=start-1. The active rectangle is overwritten on reuse;
	 * cells outside it are stale. Call traceback before another fill on this instance. */
	public void fill(byte[] read, byte[] ref, int refStartLoc, int refEndLoc){
		assert(read!=null && ref!=null) : "Global fill requires query and reference arrays";
		assert(refStartLoc>=0 && refStartLoc<=ref.length && refEndLoc>=refStartLoc-1 && refEndLoc<ref.length) :
				"Inclusive reference window is outside the array: "+refStartLoc+".."+refEndLoc+" / "+ref.length;
		rows=read.length;
		columns=refEndLoc-refStartLoc+1;
		assert(rows<=maxRows && columns<=maxColumns) :
				"Fill exceeds allocated matrix capacity: "+rows+" x "+columns+" > "+maxRows+" x "+maxColumns;
		for(int row=0; row<rows; row++){
			for(int col=0; col<columns; col++){
				final int match=(read[row]==ref[refStartLoc+col] ? 1 : -1);
				final int diag=match+scores[row][col];
				final int left=scores[row+1][col]-1;
				final int up=scores[row][col+1]-1;
				if(diag>=left && diag>=up){
					scores[row+1][col+1]=diag;
					pointers[row+1][col+1]=DIAG;
				}else if(left>=up){
					scores[row+1][col+1]=left;
					pointers[row+1][col+1]=LEFT;
				}else{
					scores[row+1][col+1]=up;
					pointers[row+1][col+1]=UP;
				}
			}
		}
	}

	/** Returns a newly allocated aligned query for the preceding fill.
	 * Arguments must describe that fill. Up moves emit query bases aligned to a
	 * reference gap; left moves emit '-'. Output may exceed either input length. */
	public byte[] traceback(byte[] read, byte[] ref, int refStartLoc, int refEndLoc){
		assert(read!=null && ref!=null) : "Traceback requires the preceding fill's input arrays";
		assert(refStartLoc>=0 && refStartLoc<=ref.length && refEndLoc>=refStartLoc-1 && refEndLoc<ref.length) :
				"Traceback reference window is outside the array: "+refStartLoc+".."+refEndLoc;
		int row=read.length, col=refEndLoc-refStartLoc+1;
		assert(row==rows && col==columns) : "Traceback dimensions differ from the last fill: "+row+" x "+col+" versus "+rows+" x "+columns;
		final byte[] out=new byte[row+col];
		int outPos=out.length;
		while(row>0 || col>0){
			final byte ptr=pointers[row][col];
			if(ptr==DIAG){
				out[--outPos]=read[--row];
				col--;
			}else if(ptr==LEFT){
				out[--outPos]='-';
				col--;
			}else{
				assert(ptr==UP) : "Traceback pointer is not DIAG/LEFT/UP: "+ptr+" at "+row+", "+col;
				out[--outPos]=read[--row];
			}
		}
		return outPos==0 ? out : Arrays.copyOfRange(out, outPos, out.length);
	}

	public final int maxRows;
	public final int maxColumns;
	private final int[][] scores;
	private final byte[][] pointers;
	private int rows, columns;

	/** Left consumes reference only; diagonal consumes both; up consumes query only. */
	public static final byte LEFT=0, DIAG=1, UP=2;
}
