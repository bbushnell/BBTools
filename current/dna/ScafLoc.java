package dna;

/**
 * @author Brian Bushnell
 * @date Sep 24, 2013
 *
 */
public class ScafLoc {
	
	/**
	 * Constructs a scaffold location with the specified coordinates.
	 * @param name_ Name identifier for the scaffold location
	 * @param chrom_ Chromosome or scaffold identifier
	 * @param loc_ Position within the chromosome/scaffold
	 */
	public ScafLoc(String name_, int chrom_, int loc_){
		name=name_;
		chrom=chrom_;
		loc=loc_;
	}

	/** Name identifier for the scaffold location */
	public String name;
	/** Chromosome or scaffold identifier */
	public int chrom;
	/** Position within the chromosome/scaffold */
	public int loc;
	
}
