package align2;

import java.util.Arrays;

import shared.Tools;
/**
 * Implementation of the Needleman-Wunsch global sequence alignment algorithm.
 * Uses dynamic programming to find optimal global alignment between two sequences
 * with configurable match/mismatch scoring and gap penalties.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class NeedlemanWunsch {
	
	
	/**
	 * Test program entry point that demonstrates alignment of two input sequences.
	 * Takes two command-line arguments as sequences and displays the scoring matrix
	 * and final alignment result.
	 * @param args Two sequence strings to align
	 */
	public static void main(String[] args){
		byte[] read=args[0].getBytes();
		byte[] ref=args[1].getBytes();
		NeedlemanWunsch nw=new NeedlemanWunsch(read.length, ref.length);
		nw.fill(read, ref, 0, ref.length-1);
		
		for(int row=0; row<nw.scores.length; row++){
			System.err.println(Arrays.toString(nw.scores[row]));
			System.err.println(Arrays.toString(nw.pointers[row]));
			System.err.println();
		}
		
		byte[] out=nw.traceback(read, ref,  0, ref.length-1);
		
		
		
		System.err.println(new String(out));
	}
	
	
	/**
	 * Constructs a Needleman-Wunsch aligner with specified matrix dimensions.
	 * Initializes scoring and pointer matrices with gap penalties along borders.
	 * Uses linear gap penalty of -1 per gap.
	 *
	 * @param maxRows_ Maximum number of rows (query sequence length + 1)
	 * @param maxColumns_ Maximum number of columns (reference sequence length + 1)
	 */
	public NeedlemanWunsch(int maxRows_, int maxColumns_){
		maxRows=maxRows_;
		maxColumns=maxColumns_;
		scores=new int[maxRows+1][maxColumns+1];
		pointers=new byte[maxRows+1][maxColumns+1];
		for(int i=0; i<maxColumns+1; i++){
			scores[0][i]=0-i;
		}
		for(int i=0; i<maxRows+1; i++){
			scores[i][0]=0-i;
		}
	}
	
//	public void initialize(int rows_, int columns_){
//		rows=rows_;
//		columns=columns_;
//		assert(rows<=maxRows);
//		assert(columns<=maxColumns);
//	}
	
	/**
	 * Fills the dynamic programming matrix for sequence alignment.
	 * Computes optimal scores using match/mismatch scoring (+1/-1) and gap penalty (-1).
	 * Updates both scoring matrix and pointer matrix for traceback.
	 *
	 * @param read Query sequence to align
	 * @param ref Reference sequence
	 * @param refStartLoc Start position in reference sequence
	 * @param refEndLoc End position in reference sequence (inclusive)
	 */
	public void fill(byte[] read, byte[] ref, int refStartLoc, int refEndLoc){
		rows=read.length;
		columns=refEndLoc-refStartLoc+1;
		System.err.println("rows = "+rows+", columns="+columns);
		
		for(int row=0; row<rows; row++){
			for(int col=0; col<columns; col++){
				System.err.println("row = "+row+", col="+col);
				int match=(read[row]==ref[refStartLoc+col] ? 1 : -1);
				int diag=match+scores[row][col];
				int left=scores[row+1][col]-1;
				int up=scores[row][col+1]-1;
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
	
	/**
	 * Performs traceback through pointer matrix to construct optimal alignment.
	 * Traces from bottom-right corner back to origin following optimal path.
	 * Returns aligned query sequence with gaps represented as '-' characters.
	 *
	 * @param read Original query sequence
	 * @param ref Original reference sequence
	 * @param refStartLoc Start position in reference sequence
	 * @param refEndLoc End position in reference sequence (inclusive)
	 * @return Aligned query sequence with gaps
	 */
	public byte[] traceback(byte[] read, byte[] ref, int refStartLoc, int refEndLoc){
		int row=read.length;
		int col=ref.length;
		
		byte[] out=new byte[Tools.max(row, col)];
		int outPos=out.length-1;
		
		while(row>0 || col>0){
			byte ptr=pointers[row][col];
			if(ptr==DIAG){
				out[outPos]=read[row-1];
				row--;
				col--;
				outPos--;
			}else if(ptr==LEFT){
				out[outPos]='-';
				col--;
				outPos--;
			}else{
				assert(ptr==UP);
//				out[outPos]='-';
				row--;
			}
		}
		return out;
	}
	
	/** Maximum number of rows in the alignment matrix */
	public final int maxRows;
	/** Maximum number of columns in the alignment matrix */
	public final int maxColumns;
	/** Dynamic programming matrix storing alignment scores */
	private final int[][] scores;
	/** Matrix storing traceback pointers for alignment reconstruction */
	private final byte[][] pointers;
	
	/** Pointer constant indicating gap in reference sequence (move up in matrix) */
	/** Pointer constant indicating match/mismatch (diagonal move in matrix) */
	/** Pointer constant indicating gap in query sequence (move left in matrix) */
	public static final byte LEFT=0, DIAG=1, UP=2;
	
	/** Current number of rows being used in the alignment matrix */
	private int rows;
	/** Current number of columns being used in the alignment matrix */
	private int columns;
	
}
