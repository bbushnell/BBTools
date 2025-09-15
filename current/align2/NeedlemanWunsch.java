package align2;

import java.util.Arrays;

import shared.Tools;
public class NeedlemanWunsch {
	
	
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
	
	public final int maxRows;
	public final int maxColumns;
	private final int[][] scores;
	private final byte[][] pointers;
	
	public static final byte LEFT=0, DIAG=1, UP=2;
	
	private int rows;
	private int columns;
	
}
