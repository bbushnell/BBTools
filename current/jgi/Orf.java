package jgi;

import java.util.Arrays;

import shared.Tools;
import structures.CoverageArray;

/**
 * Represents an Open Reading Frame (ORF) with coverage analysis capabilities.
 * Calculates various coverage metrics including base coverage, depth statistics,
 * and provides comparison functionality for sorting ORFs.
 *
 * @author Brian Bushnell
 * @date May 13, 2013
 */
public class Orf implements Comparable<Orf>{
	
	/**
	 * Constructs an ORF with name, coordinates, and strand.
	 * @param name_ ORF identifier
	 * @param start_ Start position (0-based)
	 * @param stop_ Stop position (0-based, inclusive)
	 * @param strand_ Strand orientation
	 */
	public Orf(String name_, int start_, int stop_, byte strand_){
		name=name_;
		start=start_;
		stop=stop_;
		strand=strand_;
		assert(stop>start || (start==0 && stop==0));
	}
	
	@Override
	public String toString(){
		return name+"\t"+start+"\t"+stop+"\t"+strand;
	}

	/** Returns ORF length in bases (stop - start + 1).
	 * @return Length of the ORF */
	public int length(){return stop-start+1;}

	/** Calculates average coverage depth across the ORF (0 if length is 0).
	 * @return Average base depth */
	public double avgCoverage(){
		int len=length();
		return len<=0 ? 0 : baseDepth/(double)len;
	}
	
	/** Fraction of bases with coverage > 1 across the ORF (0 if length is 0).
	 * @return Coverage fraction in [0,1] */
	public double fractionCovered(){
		int len=length();
		return len<=0 ? 0 : baseCoverage/(double)len;
	}
	
	/**
	 * Reads coverage from a CoverageArray, populating depth stats (min/max/median/stdev) and returning per-base coverage.
	 * Counts bases with coverage > 1 toward coverage/depth.
	 * @param ca Coverage array source
	 * @return Per-base coverage array for the ORF, or null if invalid
	 */
	public int[] readCoverageArray(CoverageArray ca){
		
		final int len=length();
		if(len<1 || ca==null){return null;}
		final int[] array=new int[len];
		
		baseCoverage=0;
		baseDepth=0;
		minDepth=Integer.MAX_VALUE;
		maxDepth=0;
		medianDepth=0;
		stdevDepth=0;
		
		for(int i=start, j=0; i<=stop; i++, j++){
			int cov=ca.get(i);
			array[j]=cov;
			if(cov>1){
				baseCoverage++;
				baseDepth+=cov;
				minDepth=Tools.min(minDepth, cov);
				maxDepth=Tools.max(maxDepth, cov);
			}
		}
		if(baseDepth>0){
			Arrays.sort(array);
			medianDepth=array[array.length/2];
			stdevDepth=Tools.standardDeviation(array);
		}
		return array;
	}
	
	/**
	 * Compares ORFs for sorting by name, then by start position, stop position, and strand.
	 * Uses reverse ordering for positions (higher positions first).
	 * @param o The other ORF to compare to
	 * @return Negative, zero, or positive for less than, equal to, or greater than
	 */
	@Override
	public int compareTo(Orf o) {
		int x=name.compareTo(o.name);
		if(x!=0){return x;}
		x=o.start-start;
		if(x!=0){return x;}
		x=o.stop-stop;
		if(x!=0){return x;}
		return o.strand-strand;
	}
	
	@Override
	public boolean equals(Object o){return equals((Orf)o);}
	/**
	 * Tests equality by delegating to compareTo (all comparison fields equal).
	 * @param o Other ORF
	 * @return true if equal
	 */
	public boolean equals(Orf o){return compareTo(o)==0;}
	
	@Override
	public int hashCode(){return Integer.rotateLeft(name.hashCode(),16)^(start<<8)^(stop)^strand;}
	
	public String name;
	public int start;
	public int stop;
	public byte strand;

	public long baseCoverage;
	public long readDepth=0;
	public long baseDepth=0;
	public long minDepth=0;
	public long maxDepth=0;
	public long medianDepth=0;
	public double stdevDepth=0;
			
	
}
