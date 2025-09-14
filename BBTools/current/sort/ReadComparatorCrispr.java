package sort;

import dna.Data;
import ml.CellNet;
import ml.CellNetParser;
import ml.ScoreSequence;
import shared.Parse;
import stream.Read;
import structures.SeqCountM;

/**
 * @author Brian Bushnell
 * @date Oct 27, 2014
 *
 */


public class ReadComparatorCrispr extends ReadComparator{
	
	private ReadComparatorCrispr(){}
	
	@Override
	public int compare(Read r1, Read r2) {
		return ascending*compare(r1, r2, true);
	}
	
	private SeqCountM getSCM(Read r) {
		if(r.obj!=null) {return (SeqCountM)r.obj;}
		String id=r.id;
		int x=id.indexOf("count=");
		int count=(x>=0 ? Parse.parseInt(id, x+6) : 0);
		int y=id.indexOf("score=");
		float score=(y>=0 ? Parse.parseFloat(id, y+6) : -1);
		if((y<0 || score==-1) && net!=null) {//TODO
			synchronized(ReadComparatorCrispr.class) {
				score=ScoreSequence.score(r.bases, vec, 0, net);
			}
		}
		SeqCountM scm=new SeqCountM(r.bases);
		scm.count=count;
		scm.score=score;
		r.obj=scm;
		return scm;
	}
	
	public int compare(Read r1, Read r2, boolean compareMates) {
		SeqCountM s1=getSCM(r1);
		SeqCountM s2=getSCM(r2);
		int x=s1.compareTo(s2);
		if(x!=0) {return x;}
//		if(r1.numericID!=r2.numericID){return r1.numericID>r2.numericID ? 1 : -1;}
		return r1.id.compareTo(r2.id);
	}

	@Override
	public void setAscending(boolean asc) {
		ascending=(asc ? 1 : -1);
	}
	
	public static synchronized void loadNet() {
		if(net!=null) {return;}
		setNet(CellNetParser.load(netFile));
	}
	
	public static synchronized void setNet(CellNet net_) {
		if(net_==null) {
			net=null;
			vec=null;
			return;
		}
		net=net_.copy(false);
		vec=new float[net.numInputs()];
	}
	
	private static String netFile=Data.findPath("?crispr.bbnet.gz", false);
	
	public static final ReadComparatorCrispr comparator=new ReadComparatorCrispr();
	private static CellNet net=null;
	private static float[] vec=null;
	
	int ascending=1;
}
