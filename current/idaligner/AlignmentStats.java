package idaligner;

import structures.ByteBuilder;

/**
 * Results container for IDAligner.
 * 
 * @author Brian Bushnell
 * @contributor Collei
 * @date December 13, 2025
 */
public class AlignmentStats {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	public AlignmentStats() {}
	public AlignmentStats(boolean trace) {doTrace=trace;}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	// Reset method for reuse
	public void clear() {
		identity=0;
		matches=0;
		subs=0;
		ins=0;
		dels=0;
		ns=0;
		rStart=0;
		rStop=0;
		qLen=0;
		rLen=0;
		matchString=null;
	}

	// Helper to fill from the old int[] pos style
	public void setFromPos(int[] pos, float id) {
		if(pos!=null){
			rStart=pos[0];
			rStop=pos[1];
			score=(pos.length>2 ? pos[2] : -1);
			dels=(pos.length>3 ? pos[3] : -1);
//			matches=subs=ins=ns=-1;
		}
		identity=id;
	}

	// Helper to fill from the old int[] pos style
	public float setAndSolve(int[] pos, int qLen_, int rLen_) {
		rStart=pos[0];
		rStop=pos[1];
		score=pos[2];
		dels=pos[3];
		qLen=qLen_;
		rLen=rLen_;
		return solve();
	}
	
	public float setFromMatchString(byte[] matchString_) {
		assert(matchString_!=null);
		matchString=matchString_;
		matches=subs=ins=dels=ns=0;
		for(byte b : matchString) {
			switch(b) {
				case('m') : matches++; break;
				case('S') : subs++; break;
				case('I') : ins++; break;
				case('D') : dels++; break;
				case('N') : ns++; break;
				default : assert(false) : "Unknown symbol "+(char)b;
			}
		}
		score=matches-subs-ins-dels;
		int matchQLen=matches+subs+ins+ns;
		int matchRLen=matches+subs+dels+ns;
		return identity=(matches*1f)/(matches+subs+ins+dels+ns);
	}
	
	public void printOps() {
		System.err.println(this);
	}
	
	public String toString() {
		ByteBuilder bb=new ByteBuilder();
		bb.append("rStart=").append(rStart).nl();
		bb.append("rStop=").append(rStop).nl();
		bb.append("qLen=").append(qLen).nl();
		bb.append("rLen=").append(rLen).nl();
		bb.append("matches=").append(matches).nl();
		bb.append("refAlnLength=").append(rStop-rStart+1).nl();
		bb.append("rawScore=").append(score).nl();
		bb.append("deletions=").append(dels).nl();
		bb.append("matches=").append(matches).nl();
		bb.append("substitutions=").append(subs).nl();
		bb.append("insertions=").append(ins).nl();
		bb.append("identity=").append(identity, 6).nl();
		return bb.toString();
	}
	
	public float solve() {
		// Solve the system of equations:
		// 1. M+S+I=qLen
		// 2. M+S+D=refAlnLength
		// 3. Score=M-S-I-D
		
		// Calculate operation counts
		int refAlnLength=rStop-rStart+1;
		float insF=Math.max(0, qLen+dels-refAlnLength);
		float matchesF=((score+qLen+dels)/2f);
		float subsF=Math.max(0, qLen-matchesF-insF);
		identity=matches/(matchesF+subsF+insF+dels);
		ins=(int)insF;
		matches=(int)matchesF;
		subs=(int)subsF;
		return identity;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	public byte[] matchString; // Long Match string
	public float identity;
	public int matches, subs, ins, dels, ns;
	public int score; // Optional raw score
	public int rStart, rStop;
	public int qLen, rLen;
	public boolean doTrace=false;
	
}