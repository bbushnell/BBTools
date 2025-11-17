package jgi;

import java.util.Arrays;

import shared.Tools;
import structures.CoverageArray;

/**
 * This class is designed to help calculate coverage of ORFs
 * @author Brian Bushnell
 * @date May 13, 2013
 *
 */
public class Orf implements Comparable<Orf>{
	
	/**
	 * Constructs an ORF with position and strand information.
	 *
	 * @param name_ Name identifier for this ORF
	 * @param start_ Start position (0-based)
	 * @param stop_ Stop position (0-based, inclusive)
	 * @param strand_ Strand orientation as byte value
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

	/** Gets the length of the ORF in bases */
	public int length(){return stop-start+1;}

	/** Calculates average coverage depth across the ORF.
	 * @return Average base depth, or 0 if ORF has no length */
	public double avgCoverage(){
		int len=length();
		return len<=0 ? 0 : baseDepth/(double)len;
	}
	
	/** Calculates fraction of ORF bases with non-zero coverage.
	 * @return Fraction of bases covered (0.0 to 1.0), or 0 if ORF has no length */
	public double fractionCovered(){
		int len=length();
		return len<=0 ? 0 : baseCoverage/(double)len;
	}
	
	/**
	 * Reads coverage data from a CoverageArray and calculates statistics.
	 * Populates all coverage metrics including min/max/median depth and standard deviation.
	 * Only counts bases with coverage > 1 towards coverage statistics.
	 *
	 * @param ca CoverageArray containing per-base coverage data
	 * @return Array of coverage values for each base in the ORF, or null if invalid
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
	 * Tests equality with another ORF using compareTo method.
	 * @param o ORF to compare with
	 * @return true if ORFs are equal (compareTo returns 0)
	 */
	public boolean equals(Orf o){return compareTo(o)==0;}
	
	@Override
	public int hashCode(){return Integer.rotateLeft(name.hashCode(),16)^(start<<8)^(stop)^strand;}
	
	/** Name of ORF (not necessarily the name of its scaffold) */
	public String name;
	/** Start position of the ORF (0-based) */
	public int start;
	/** Stop position of the ORF (0-based, inclusive) */
	public int stop;
	/** Strand orientation of the ORF */
	public byte strand;

	/** Number of bases with nonzero coverage */
	public long baseCoverage;
	/** Number of reads mapped to this orf */
	public long readDepth=0;
	/** Number of bases mapped to this orf */
	public long baseDepth=0;
	/** Lowest base depth */
	public long minDepth=0;
	/** Highest base depth */
	public long maxDepth=0;
	/** Median base depth */
	public long medianDepth=0;
	/** Standard deviation of depth */
	public double stdevDepth=0;
			
	
}
