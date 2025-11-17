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
	
	/**
	 * Private constructor prevents instantiation; use the static comparator instance
	 */
	private ReadComparatorCrispr(){}
	
	@Override
	public int compare(Read r1, Read r2) {
		return ascending*compare(r1, r2, true);
	}
	
	/**
	 * Extracts or creates a SeqCountM object containing count and score metadata for a read.
	 * Parses count from the read ID using "count=" prefix and score using "score=" prefix.
	 * If score is missing or invalid and a neural network is loaded, computes score using
	 * the network in a thread-safe manner.
	 *
	 * @param r The read to extract metadata from
	 * @return SeqCountM object with count and score information
	 */
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
	
	/**
	 * Compares two reads based on their SeqCountM metadata.
	 * First compares using SeqCountM.compareTo(), then falls back to string ID comparison
	 * for tie-breaking.
	 *
	 * @param r1 First read to compare
	 * @param r2 Second read to compare
	 * @param compareMates Whether to include mate comparison (currently unused)
	 * @return Negative if r1 < r2, positive if r1 > r2, zero if equal
	 */
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
	
	/** Loads the neural network from the default CRISPR network file if not already loaded.
	 * Uses thread-safe initialization to prevent multiple loading attempts. */
	public static synchronized void loadNet() {
		if(net!=null) {return;}
		setNet(CellNetParser.load(netFile));
	}
	
	/**
	 * Sets the neural network for sequence scoring with thread-safe initialization.
	 * Creates a copy of the provided network and initializes the input vector array.
	 * Setting null clears both the network and vector.
	 * @param net_ The neural network to use, or null to clear
	 */
	public static synchronized void setNet(CellNet net_) {
		if(net_==null) {
			net=null;
			vec=null;
			return;
		}
		net=net_.copy(false);
		vec=new float[net.numInputs()];
	}
	
	/** Path to the default CRISPR neural network file (crispr.bbnet.gz) */
	private static String netFile=Data.findPath("?crispr.bbnet.gz", false);
	
	/** Singleton instance of the CRISPR read comparator */
	public static final ReadComparatorCrispr comparator=new ReadComparatorCrispr();
	/**
	 * Neural network used for sequence scoring when score metadata is unavailable
	 */
	private static CellNet net=null;
	/** Input vector for neural network scoring operations */
	private static float[] vec=null;
	
	/**
	 * Sorting direction multiplier: 1 for ascending order, -1 for descending order
	 */
	int ascending=1;
}
