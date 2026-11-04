package dna;

/**
 * Mutable scaffold-name entry locating its start within a concatenated chromosome.
 * Data.scafNameTable creates these from scaffoldNames/scaffoldLocs and shares the
 * entries with callers. This class performs no validation or synchronization.
 * @author Brian Bushnell
 * @date Sep 24, 2013
 */
public class ScafLoc{

	/** Stores the name, chromosome identifier and zero-based scaffold start as supplied. */
	public ScafLoc(String name_, int chrom_, int loc_){
		name=name_;
		chrom=chrom_;
		loc=loc_;
	}

	/** Scaffold identifier used as the lookup-table key when the entry is created. */
	public String name;
	/** Identifier of the chromosome containing this scaffold. */
	public int chrom;
	/** Zero-based start within the chromosome, including preceding padding/scaffolds. */
	public int loc;

}
