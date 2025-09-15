package aligner;

import dna.AminoAcid;
import shared.KillSwitch;
import shared.Tools;

/**
 * Based on SSAFlat2, but reduced to a 1D array. */
public final class SingleStateAlignerFlat2_1D implements Aligner, IDAligner {
	
	
	public SingleStateAlignerFlat2_1D(){}

	@Override
	public final String name() {return "SSAF2_1D";}
	@Override
	public long loops() {return -1;}
	@Override
	public void setLoops(long x) {};//Not supported
	@Override
	public final float align(byte[] a, byte[] b) {return align(a, b, null, 0, b.length-1, -9999);}
	@Override
	public final float align(byte[] a, byte[] b, int[] pos) {return align(a, b, pos, 0, b.length-1, -9999);}
	@Override
	public final float align(byte[] a, byte[] b, int[] posVector, int minScore) {return align(a, b, null, 0, b.length-1, minScore);}
	@Override
	public final float align(byte[] a, byte[] b, int[] pos, int from, int to) {return align(a, b, pos, from, to, -9999);}
	public float align(byte[] q, byte[] r, int[] pos, int from, int to, int minScore) {
		if(q.length>r.length && pos==null) {byte[] s=q; q=r; r=s;}
		assert(q.length<=r.length);
		int[] max=fillUnlimited(q, r, from, to, minScore);
		if(max==null){return 0;}

		final int rows=max[0];
		final int maxCol=max[1];
		final int maxState=max[2];
		final float id=tracebackIdentity(q, r, from, to, rows, maxCol, maxState, null);
		if(pos!=null) {
			final int[] score=score(q, r, from, to, rows, maxCol, maxState);
			final int rstart=Tools.max(score[1], from);
			final int rstop=Tools.min(score[2], to);
			pos[0]=rstart;
			pos[1]=rstop;
		}
		return id;
	}
	
	private void prefillTopRow(){
		final int qlen=rows;
		for(int i=0; i<=columns; i++){
			int x=columns-i+1;
			int qbases=qlen-x;
			
			//Minimal points to prefer a leftmost alignment
			set(0, i, qbases<=0 ? 0 : -qbases);
			
			//Forces consumption of query, but does not allow for insertions...
//			header[i]=qbases<=0 ? 0 : calcDelScoreOffset(qbases);
		}
	}
	
	private void prefillLeftColumnStartingAt(int i){
		set(0, 0, MODE_MATCH);
		i=Tools.max(1, i);
		for(int score=MODE_INS+(POINTS_INS*i); i<=rows; i++){//Fill column 0 with insertions
			score+=POINTS_INS;
			set(i, 0, score);
		}
	}
	
	private final int index(int row, int col){
		return col+row*rowMult;
	}
	
	private final int get(int row, int col){
		return matrix[col+row*rowMult];
	}
	
	private final void set(int row, int col, int value){
		matrix[col+row*rowMult]=value;
	}
	
	private void initialize(int rows_, int columns_){
		final int rows0=rows;
		final int columns0=columns;
		rows=rows_;
		columns=columns_;
		rowMult=columns+1;
		int size=(rows+1)*rowMult;
		if(matrix==null || size>matrix.length){
			matrix=KillSwitch.allocInt1D((rows+5)*(columns+10));
		}
		prefillTopRow();
		if(rows>rows0 || columns!=columns0){
			prefillLeftColumnStartingAt(1);
		}
	}
	
	/** return new int[] {rows, maxCol, maxState, maxScore, maxStart};
	 * Will not fill areas that cannot match minScore */
	@Override
	public final int[] fillLimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore){
		return fillUnlimited(read, ref, refStartLoc, refEndLoc, minScore);
	}
	
	@Override
	public final int[] fillUnlimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc){
		return fillUnlimited(read, ref, refStartLoc, refEndLoc, -999999);
	}
	
	/** return new int[] {rows, maxCol, maxState, maxScore, maxStart};
	 * Min score is optional */
	@Override
	public final int[] fillUnlimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore){
		initialize(read.length, refEndLoc-refStartLoc+1);
		
		//temporary, for finding a bug
		if((rows+1)*rowMult>matrix.length){
			throw new RuntimeException("rows="+rows+", cols="+columns+", len="+matrix.length+"\n"+new String(read)+"\n");
		}
		
		assert(refStartLoc>=0) : "Check that values are in-bounds before calling this function: "+refStartLoc;
		assert(refEndLoc<ref.length) : "Check that values are in-bounds before calling this function: "+refEndLoc+", "+ref.length;
		
		final int refOffset=refStartLoc-1;
		for(int row=1; row<=rows; row++){

			final byte qBase=read[row-1];
			for(int col=1; col<=columns; col++){
				
				final byte rBase=ref[refOffset+col];
				
				final boolean match=(qBase==rBase);
				final boolean defined=(qBase!='N' && rBase!='N');

				final int scoreFromDiag=get(row-1, col-1);
				final int scoreFromDel=get(row, col-1);
				final int scoreFromIns=get(row-1, col);
				
				final int diagScoreM=POINTS_MATCH;
				final int diagScoreS=POINTS_SUB;
				final int delScore=scoreFromDel+POINTS_DEL;
				final int insScore=scoreFromIns+POINTS_INS;
				
//				final int diagScore=scoreFromDiag+(defined ? (match ? diagScoreM : diagScoreS) : POINTS_NOREF);
				int diagScore=(match ? diagScoreM : diagScoreS);
				diagScore=scoreFromDiag+(defined ? diagScore : POINTS_NOREF);
				
				int score=diagScore>=delScore ? diagScore : delScore;
				score=score>=insScore ? score : insScore;
				
				set(row, col, score);
			}
			//iterationsUnlimited+=columns;
		}
		

		int maxCol=-1;
		int maxState=-1;
		int maxStart=-1;
		int maxScore=Integer.MIN_VALUE;
		
		for(int col=1; col<=columns; col++){
			int x=get(rows, col);
			if(x>maxScore){
				maxScore=x;
				maxCol=col;

//				assert(rows-1<read.length) : (rows-1)+", "+read.length;
//				assert(refOffset+col<ref.length) : refOffset+", "+col+", "+ref.length;
				maxState=getState(rows, col, read[rows-1], ref[refOffset+col]);
				maxStart=x;
			}
		}

//		System.err.println("Returning "+rows+", "+maxCol+", "+maxState+", "+maxScore+"; minScore="+minScore);
		return maxScore<minScore ? null : new int[] {rows, maxCol, maxState, maxScore, maxStart};
	}
	
	int getState(int row, int col, byte q, byte r){//zxvzxcv TODO: Fix - needs to find max
		final boolean match=(q==r);
		final boolean defined=(q!='N' && r!='N');
		
		final int scoreFromDiag=get(row-1, col-1);
		final int scoreFromDel=get(row, col-1);
		final int scoreFromIns=get(row-1, col);
//		final int score=packed[row][col];
		
		final int diagScoreM=POINTS_MATCH;
		final int diagScoreS=POINTS_SUB;
		final int delScore=scoreFromDel+POINTS_DEL;
		final int insScore=scoreFromIns+POINTS_INS;
		
		final int diagScore=scoreFromDiag+(defined ? (match ? diagScoreM : diagScoreS) : POINTS_NOREF);
		
//		int score2=diagScore>=delScore ? diagScore : delScore;
//		score2=score>=insScore ? score : insScore;
		
//		assert(score==score2) : score+", "+score2;
		
		if(diagScore>=delScore && diagScore>=insScore){
			return defined ? match ? MODE_MATCH : MODE_SUB : MODE_N;
		}else if(delScore>=insScore){
			return MODE_DEL;
		}
		return MODE_INS;
	}
	
	
	/** Generates the match string */
	@Override
	public final byte[] traceback(byte[] query, byte[] ref, int refStartLoc, int refEndLoc, int row, int col, int state){
//		assert(false);
		assert(refStartLoc<=refEndLoc) : refStartLoc+", "+refEndLoc;
		assert(row==rows);
		
		byte[] out=new byte[row+col-1]; //TODO if an out of bound crash occurs, try removing the "-1".
		int outPos=0;

//		assert(state==(packed[row][col]&MODEMASK));
		
		while(row>0 && col>0){
			byte q=query[row-1];
			byte r=ref[refStartLoc+col-1];
			boolean defined=(AminoAcid.isFullyDefined(q) && AminoAcid.isFullyDefined(r));
			state=getState(row, col, q, r);
			if(state==MODE_MATCH){
				col--;
				row--;
				out[outPos]=defined ? (byte)'m' : (byte)'N';
			}else if(state==MODE_SUB){
				col--;
				row--;
				out[outPos]=defined ? (byte)'S' : (byte)'N';
			}else if(state==MODE_N){
				col--;
				row--;
				out[outPos]='N';
			}else if(state==MODE_DEL){
				col--;
				out[outPos]='D';
			}else if(state==MODE_INS){
				row--;
//				out[outPos]='I';
				if(col>=0 && col<columns){
					out[outPos]='I';
				}else{
					out[outPos]='C';
					col--;
				}
			}else{
				assert(false) : state;
			}
			outPos++;
		}
		
		assert(row==0 || col==0);
		if(col!=row){//Not sure what this is doing
			while(row>0){
				out[outPos]='C';
				outPos++;
				row--;
				col--;
			}
			if(col>0){
				//do nothing
			}
		}
		
		//Shrink and reverse the string
		byte[] out2=new byte[outPos];
		for(int i=0; i<outPos; i++){
			out2[i]=out[outPos-i-1];
		}
		out=null;
		
		return out2;
	}
	
	@Override
	/** Generates identity;
	 * fills 'extra' with {match, sub, del, ins, N, clip} if present */
	public float tracebackIdentity(byte[] query, byte[] ref, int refStartLoc, int refEndLoc, int row, int col, int state, int[] extra){

//		assert(false);
		assert(refStartLoc<=refEndLoc) : refStartLoc+", "+refEndLoc;
		assert(row==rows);

//		assert(state==(packed[row][col]&MODEMASK));
		int match=0, sub=0, del=0, ins=0, noref=0, clip=0;
		
		while(row>0 && col>0){
			byte q=query[row-1];
			byte r=ref[refStartLoc+col-1];
			boolean defined=(AminoAcid.isFullyDefined(q) && AminoAcid.isFullyDefined(r));
			state=getState(row, col, q, r);
			if(state==MODE_MATCH){
				col--;
				row--;
				match+=(defined ? 1 : 0);
				noref+=(defined ? 0 : 1);
			}else if(state==MODE_SUB){
				col--;
				row--;
				sub+=(defined ? 1 : 0);
				noref+=(defined ? 0 : 1);
			}else if(state==MODE_N){
				col--;
				row--;
				noref++;
			}else if(state==MODE_DEL){
				col--;
				del++;
			}else if(state==MODE_INS){
				row--;
				boolean edge=(col<=1 || col>=columns);
				ins+=(edge ? 0 : 1);
				clip+=(edge ? 1 : 0);
			}else{
				assert(false) : state;
			}
		}
		
		assert(row==0 || col==0);
		if(col!=row){//Not sure what this is doing
			while(row>0){
				clip++;
				row--;
				col--;
			}
			if(col>0){
				//do nothing
			}
		}
		
		if(extra!=null){
			assert(extra.length==5);
			extra[0]=match;
			extra[1]=sub;
			extra[2]=del;
			extra[3]=ins;
			extra[4]=noref;
			extra[5]=clip;
		}
		
		float len=match+sub+ins+del+noref*0.1f;
		float id=match/Tools.max(1.0f, len);
		return id;
	}
	
	/** @return {score, bestRefStart, bestRefStop} */
	@Override
	public final int[] score(final byte[] read, final byte[] ref, final int refStartLoc, final int refEndLoc,
			final int maxRow, final int maxCol, final int maxState/*, final int maxScore, final int maxStart*/){
		
		int row=maxRow;
		int col=maxCol;
		int state=maxState;
		
		assert(maxRow>=0 && maxRow<matrix.length) :
			maxState+", "+maxRow+", "+maxCol+"\n"+new String(read)+"\n"+toString(ref, refStartLoc, refEndLoc); //TODO: Change to translate col and row
		
		int score=get(maxRow, maxCol); //Or zero, if it is to be recalculated
		
		if(row<rows){
			int difR=rows-row;
			int difC=columns-col;
			
			while(difR>difC){
				score+=POINTS_NOREF;
				difR--;
			}
			
			row+=difR;
			col+=difR;
			
		}
		
		assert(refStartLoc<=refEndLoc);
		assert(row==rows);

		
		final int bestRefStop=refStartLoc+col-1;
		
		while(row>0 && col>0){
			final byte q=read[row-1];
			final byte r=ref[refStartLoc+col-1];
//			final boolean defined=(AminoAcid.isFullyDefined(q) && AminoAcid.isFullyDefined(r));
			state=getState(row, col, q, r);
			if(state==MODE_MATCH){
				col--;
				row--;
			}else if(state==MODE_SUB){
				col--;
				row--;
			}else if(state==MODE_N){
				col--;
				row--;
			}else if(state==MODE_DEL){
				col--;
			}else if(state==MODE_INS){
				row--;
			}else{
				assert(false) : state;
			}
		}
//		assert(false) : row+", "+col;
		if(row>col){
			col-=row;
		}
		
		final int bestRefStart=refStartLoc+col;
		
//		System.err.println("t2\t"+score+", "+maxScore+", "+maxStart+", "+bestRefStart);
		int[] rvec;
		if(bestRefStart<refStartLoc || bestRefStop>refEndLoc){ //Suggest extra padding in cases of overflow
			int padLeft=Tools.max(0, refStartLoc-bestRefStart);
			int padRight=Tools.max(0, bestRefStop-refEndLoc);
			rvec=new int[] {score, bestRefStart, bestRefStop, padLeft, padRight};
		}else{
			rvec=new int[] {score, bestRefStart, bestRefStop};
		}
		return rvec;
	}
	
	
	/** Will not fill areas that cannot match minScore.
	 * @return {score, bestRefStart, bestRefStop}  */
	@Override
	public final int[] fillAndScoreLimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore){
		int a=Tools.max(0, refStartLoc);
		int b=Tools.min(ref.length-1, refEndLoc);
		assert(b>=a);
		int rlen=b-a+1;
		long size=(read.length+1L)*(rlen+1L);
		assert(b>=a);
		
//		if(b-a>=maxColumns){
//			System.err.println("Warning: Max alignment columns exceeded; restricting range. "+(b-a+1)+" > "+maxColumns);
//			assert(false) : refStartLoc+", "+refEndLoc;
//			b=Tools.min(ref.length-1, a+maxColumns-1);
//		}
		int[] max=fillLimited(read, ref, a, b, minScore);
//		return max==null ? null : new int[] {max[3], 0, max[1]};
		
		int[] score=(max==null ? null : score(read, ref, a, b, max[0], max[1], max[2]/*, max[3], max[4]*/));
		
		return score;
	}
	
	public static final String toString(byte[] ref, int startLoc, int stopLoc){
		StringBuilder sb=new StringBuilder(stopLoc-startLoc+1);
		for(int i=startLoc; i<=stopLoc; i++){sb.append((char)ref[i]);}
		return sb.toString();
	}
	
//	public int maxScoreByIdentity(int len, float identity){
//		assert(identity>=0 && identity<=1);
//		return (int)(len*(identity*POINTS_MATCH+(1-identity)*POINTS_SUB));
//	}
	
	@Override
	public int minScoreByIdentity(int len, float identity){
		assert(identity>=0 && identity<=1);
		
		int a=(int)(len*(identity*POINTS_MATCH+(1-identity)*POINTS_SUB));
		int b=(int)(len*(identity*POINTS_MATCH+(1-identity)*POINTS_INS));
		int c=(int)(len*(1*POINTS_MATCH+((1/(Tools.max(identity, 0.000001f)))-1)*POINTS_DEL));
		return Tools.min(a, b, c);
	}
	
	private static int calcDelScore(int len){
		if(len<=0){return 0;}
		int score=POINTS_DEL*len;
		return score;
	}
	
	@Override
	public int rows(){return rows;}
	@Override
	public int columns(){return columns;}

	private int[] matrix;
	
	public static final int MAX_SCORE=Integer.MAX_VALUE-2000;
	public static final int MIN_SCORE=0-MAX_SCORE; //Keeps it 1 point above "BAD".

	//For some reason changing MODE_DEL from 1 to 0 breaks everything
	private static final byte MODE_DEL=1;
	private static final byte MODE_INS=2;
	private static final byte MODE_SUB=3;
	private static final byte MODE_MATCH=4;
	private static final byte MODE_N=5;
	
	public static final int POINTS_NOREF=-15;
	public static final int POINTS_MATCH=100;
	public static final int POINTS_SUB=-50;
	public static final int POINTS_INS=-121;
	public static final int POINTS_DEL=-111;
	
//	public static final int BAD=MIN_SCORE-1;
	
	private int rows;
	private int columns;
	/** columns+1 **/
	private int rowMult;

	public boolean verbose=false;
	public boolean verbose2=false;
	
}
