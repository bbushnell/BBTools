package sketch;

import java.util.HashMap;

/**
 * Represents a genomic comparison record containing taxonomic identifiers,
 * sequence metrics, and similarity scores for sketch-based genomic analysis.
 * Supports cloning and ANI-based comparison ordering for result processing.
 *
 * @author Brian Bushnell
 * @date 2013
 */
class Record implements Cloneable, Comparable<Record> {
	
	/**
	 * Constructs a Record from parsed result line data.
	 * Copies all taxonomic, sequence, and similarity metrics from the parser.
	 * @param parser ResultLineParser containing parsed comparison data
	 */
	Record(ResultLineParser parser){
		qTaxID=parser.qTaxID;
		rTaxID=parser.rTaxID;
		qBases=parser.qBases;
		rBases=parser.rBases;
		qSize=parser.qSize;
		rSize=parser.rSize;
		ani=parser.ani;
		ssu=parser.ssu;
		taxLevelExtended=parser.taxLevelExtended;
		text=parser.text;
	}
	
	/** Creates a deep copy of this Record using the Cloneable interface.
	 * @return Cloned Record instance, or null if cloning fails */
	Record copy(){
		try {
			return (Record) this.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	
	/**
	 * Processes and calculates SSU (Small Subunit) similarity if not already set.
	 * Retrieves 16S rRNA sequences for both query and reference taxa from SSUMap,
	 * then performs sequence alignment to determine SSU similarity score.
	 * Only processes when ssu <= 0 and both sequences are available in the map.
	 */
	public void processSSU(){
		if(verbose){System.err.println("Record.processSSU(): qTaxID="+qTaxID+", rTaxID="+rTaxID+", ssu="+(float)ssu);}
		final HashMap<Integer, byte[]> map=SSUMap.r16SMap;
		if(ssu<=0 && map!=null){
			byte[] qssu=map.get(qTaxID);
			byte[] rssu=map.get(rTaxID);
			if(qssu!=null && rssu!=null){
				ssu=SketchObject.align(qssu, rssu);
				if(verbose){System.err.println("Aligned; ssu="+(float)ssu);}
			}else{
				if(verbose){System.err.println("Missing: "+qTaxID+"="+qssu+", "+rTaxID+"="+rssu);}
			}
		}else{
			if(verbose){System.err.println("Skipping: ssu="+(float)ssu+", map="+(map==null ? "null" : map.size()));}
		}
	}
	
	@Override
	public int compareTo(Record o) {
		return ani>o.ani ? -1 : ani<o.ani ? 1 : 0;
	}
	
	/** Gets the SSU similarity score, returning 0 if no valid score exists.
	 * @return SSU similarity score, or 0.0 if ssu <= 0 */
	public double ssu(){return ssu<=0 ? 0 : ssu;}
	
	/** Query sequence taxonomic identifier */
	final int qTaxID;
	/** Reference sequence taxonomic identifier */
	final int rTaxID;
	/** Total number of bases in query sequence */
	final long qBases;
	/** Total number of bases in reference sequence */
	final long rBases;
	/** Query genome/sequence size estimate */
	final long qSize;
	/** Reference genome/sequence size estimate */
	final long rSize;
	/** Average Nucleotide Identity score between query and reference */
	final double ani;
	/** SSU (Small Subunit) similarity score; -1 indicates unprocessed */
	private double ssu=-1;
	/** Extended taxonomic level classification for this record */
	int taxLevelExtended;
	/** NCBI taxonomic correctness flag: 0=false, 1=true */
	int correctNCBI;//0 is false, 1 is true
	/** SSU taxonomic correctness flag: 0=false, 1=true */
	int correctSSU;//0 is false, 1 is true
	/** Flag indicating whether SSU sequence data is missing */
	boolean missingSSU;
	/** Raw text representation of the original result line */
	byte[] text;
	//All records are correctNCBI there is no later record with a lower taxLevel.
	//Records are correctSSU if:
	//Their SSUID is >= that of all subsequent records, and the next record (if present) has a SSUID
	
	/** Global flag to enable verbose debugging output for Record operations */
	static boolean verbose;
}