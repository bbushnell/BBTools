package sort;

import dna.Data;
import ml.CellNet;
import ml.CellNetParser;
import ml.ScoreSequence;
import parse.Parse;
import stream.Read;
import structures.SeqCountM;



/**
 * Comparator for sorting and ranking DNA reads with optional neural network scoring.
 * Compares DNA reads based on embedded count and score metadata, with thread-safe
 * neural network model loading for sequence scoring when score metadata is missing.
 *
 * @author Brian Bushnell
 * @date Oct 27, 2014
 */
public class ReadComparatorCrispr extends ReadComparator{

	/** Private constructor; use the ascending/descending singletons. */
	private ReadComparatorCrispr(int mult_){mult=mult_;}

	@Override
	public int compare(Read r1, Read r2) {
		return mult*compare(r1, r2, true);
	}

	/** Returns the read's {@link SeqCountM} score record, building it on first use: parses count= and
	 * score= from the header, falling back to neural-network scoring when score is absent and a net is
	 * loaded. The result is cached on the read's obj field. */
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
	 * Compares two reads using SeqCountM metadata (count/score) and falls back to ID tie-breakers.
	 * @param r1 First read
	 * @param r2 Second read
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
	public ReadComparator getComparator(boolean asc){return asc ? ascending : descending;}

	@Override
	public final int ascendingMult() {return mult;}
	@Override
	public final String name() {return "Crispr";}

	/**
	 * Loads the default CRISPR neural network if not already loaded; thread-safe.
	 */
	public static synchronized void loadNet() {
		if(net!=null) {return;}
		setNet(CellNetParser.load(netFile));
	}

	/** Sets or clears the neural network used for sequence scoring; initializes input vector.
	 * @param net_ Network to use, or null to clear */
	public static synchronized void setNet(CellNet net_) {
		if(net_==null) {
			net=null;
			vec=null;
			return;
		}
		net=net_.copy(false);
		vec=new float[net.numInputs()];
	}

	/** Sort direction multiplier: +1 for ascending, -1 for descending (immutable). */
	private final int mult;

	/** Resource path to the default CRISPR scoring neural network. */
	private static String netFile=Data.findPath("?crispr.bbnet.gz", false);

	/** Ascending-order singleton. */
	public static final ReadComparatorCrispr ascending=new ReadComparatorCrispr(1);
	/** Descending-order singleton. */
	public static final ReadComparatorCrispr descending=new ReadComparatorCrispr(-1);
	/** Default (ascending) singleton; alias retained for selection/identity call sites. */
	public static final ReadComparatorCrispr comparator=ascending;
	/** Loaded neural network used to score sequences when a header score is missing; null if unloaded. */
	private static CellNet net=null;
	/** Reusable input-vector buffer sized to the network's input count. */
	private static float[] vec=null;
}
