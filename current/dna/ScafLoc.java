package dna;

/**
 * Represents a scaffold location with a name identifier and genomic coordinates.
 * Simple data structure for storing scaffold-based genomic positions.
 * @author Brian Bushnell
 * @date Sep 24, 2013
 */
public class ScafLoc {
	
	public ScafLoc(String name_, int chrom_, int loc_){
		name=name_;
		chrom=chrom_;
		loc=loc_;
	}

	public String name;
	public int chrom;
	public int loc;
	
}
