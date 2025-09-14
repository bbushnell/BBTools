package aligner;

import prok.GeneCaller;
import stream.Read;

public class Alignment implements Comparable<Alignment>{
	
	public Alignment(Read r_){
		r=r_;
	}
	
	@Override
	public int compareTo(Alignment o) {
		return id>o.id ? 1 : id<o.id ? -1 : r.length()>o.r.length() ? 1 : r.length()<o.r.length() ? -1 : 0;
	}
	
	public float align(byte[] ref){
		id=align(r, ref);
		match=r.match;
		start=r.start;
		stop=r.stop;
		return id;
	}
	
	public static final float align(Read r, byte[] ref){
		SingleStateAlignerFlat2 ssa=GeneCaller.getSSA();
		final int a=0, b=ref.length-1;
		int[] max=ssa.fillUnlimited(r.bases, ref, a, b, 0);
		if(max==null){return 0;}
		
		final int rows=max[0];
		final int maxCol=max[1];
		final int maxState=max[2];
		
		//returns {score, bestRefStart, bestRefStop} 
		//padded: {score, bestRefStart, bestRefStop, padLeft, padRight};
		int[] score=ssa.score(r.bases, ref, a, b, rows, maxCol, maxState);
		int rstart=score[1];
		int rstop=score[2];
		r.start=rstart;
		r.stop=rstop;
		
		byte[] match=ssa.traceback(r.bases, ref, a, b, rows, maxCol, maxState);
		float id=Read.identity(match);
		return id;
	}
	
	public final Read r;
	public float id=-1;
	public byte[] match;
	public int start;
	public int stop;
	
}
