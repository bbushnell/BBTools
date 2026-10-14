package fileIO;

import dna.Data;
import dna.Gene;
import shared.Shared;

/**
 * One ungapped target-to-query segment produced by {@link ChainBlock#toLines()}.
 * Target endpoints are inclusive. Query endpoints follow target traversal order:
 * ascending for PLUS, descending for the reverse branch. Coordinates are stored
 * and translated literally; this class does not convert external coordinate formats.
 * Fields are mutable and caller-owned, with no synchronization or constructor validation.
 * Search accepts overlapping intervals sorted by target start within one chromosome;
 * it prefers the earliest target start, then the earliest destination at that start.
 *
 * @author Brian Bushnell, Shinobu
 */
public class ChainLine implements Comparable<ChainLine>{

	/*--------------------------------------------------------------*/
	/*----------------            Startup           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Prints translations through Data.ROOT_CHAIN's hg18ToHg19.over.chain file.
	 * Arguments are positional; coordinates are passed through without an offset.
	 * @param args Numeric chromosome token accepted by Gene.toChromosome, then positions
	 */
	public static void main(String[] args){

		int chrom=Gene.toChromosome(args[0]);

		ChainLine[][] lines=ChainBlock.loadChainLines(Data.ROOT_CHAIN+"hg18ToHg19.over.chain");

		for(int i=1; i<args.length; i++){
			int loc=Integer.parseInt(args[i]);
			int[] result=translate(loc, lines[chrom]);
			System.out.print(chrom+"\t+\t"+loc+"\t->\t");
			System.out.println(result==null ? "null" : result[0]+"\t"+Shared.strandCodes[result[1]]+"\t"+result[2]);
		}
		
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Copies all supplied values without checking or normalizing endpoints or strands.
	 *
	 * @param chromT Target chromosome number
	 * @param strandT Encoded target strand, normally Shared.PLUS; not used in translation
	 * @param startT Inclusive lower target endpoint
	 * @param stopT Inclusive upper target endpoint
	 * @param chromQ Query chromosome number
	 * @param strandQ Encoded query strand; Shared.PLUS selects ascending translation
	 * @param startQ Query coordinate corresponding to startT
	 * @param stopQ Query coordinate corresponding to stopT
	 */
	public ChainLine(int chromT, byte strandT, int startT, int stopT, int chromQ, byte strandQ, int startQ, int stopQ){
		tChrom=chromT;
		tStrand=strandT;
		tStart=startT;
		tStop=stopT;
		
		qChrom=chromQ;
		qStrand=strandQ;
		qStart=startQ;
		qStop=stopQ;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the eight stored fields separated by tabs, rendering strand codes as symbols. */
	@Override
	public String toString(){
		return tChrom+"\t"+Shared.strandCodes[tStrand]+"\t"+tStart+"\t"+tStop+"\t"+
		qChrom+"\t"+Shared.strandCodes[qStrand]+"\t"+qStart+"\t"+qStop;
	}
	
	
	/**
	 * Finds a segment containing loc, treating null and empty arrays as unmapped.
	 * Intervals must be sorted by ascending target start within one chromosome.
	 * Among containing segments, prefers the lowest tStart, then the lowest qStart
	 * (destination at that source start, also on reverse strands). Exact ties retain
	 * array order. This is not score selection or minimum translated position across
	 * mixed strands. Entries and keys must remain stable during the search.
	 * The historical method name is retained; overlapping intervals require a scan
	 * of the relevant prefix, with linear worst-case cost and no extra index storage.
	 * @param loc Target coordinate to search for
	 * @param array Searchable segment array, or null for an uncovered chromosome
	 * @return Preferred containing index, or -1 when no segment is found
	 */
	public static int binarySearch(int loc, ChainLine[] array){
		//FIXED [fileIO/ChainLine#001]: null array means the chromosome has zero chain blocks (Data.getChainLines /
		//ChainBlock.loadChainLines legitimately leaves lines[chrom]==null for an uncovered chromosome). Return -1
		//("not found") to match the contract and both callers' existing index<0/result==null "untranslatable" handling
		//(Translator.translate:132, Translator2 via translate()); without this, array.length below NPEs on a valid
		//liftover query for an uncovered chromosome — defeating the tools' own graceful no-mapping design.
		if(array==null){return -1;}
		return binarySearch(loc, array, 0, array.length-1);
	}
	
	
	/**
	 * Searches an inclusive index range with the same interval assumptions as
	 * {@link #binarySearch(int, ChainLine[])}. An empty range returns -1 before
	 * accessing the array; a nonempty range requires a nonnull array and valid bounds.
	 *
	 * @param loc Target coordinate to search for
	 * @param array Array of nonnull segments satisfying the search assumptions
	 * @param first Inclusive first index
	 * @param last Inclusive last index; less than first denotes an empty range
	 * @return Preferred containing index, or -1 when no segment is found
	 */
	public static int binarySearch(int loc, ChainLine[] array, int first, int last){
		//FIO-052: a later short interval cannot rule out an earlier enclosing interval.
		int best=-1;
		for(int i=first; i<=last; i++){
			final ChainLine line=array[i];
			if(line.tStart>loc || (best>=0 && line.tStart>array[best].tStart)){break;}
			if(line.contains(loc) && (best<0 || line.qStart<array[best].qStart)){best=i;}
		}
		return best;
	}
	
	/**
	 * Searches the supplied segments and translates through the selected segment.
	 * Uses the interval assumptions of {@link #binarySearch(int, ChainLine[])};
	 * it does not select among competing mappings by score or chain ID.
	 * @param loc Target coordinate, with no coordinate-format conversion
	 * @param array Searchable segments, or null for an uncovered chromosome
	 * @return Fresh {query chromosome, encoded query strand, position}, or null if unmapped
	 */
	public static int[] translate(int loc, ChainLine[] array){
		int index=binarySearch(loc, array);
		if(index<0){return null;}
		ChainLine cl=array[index];
		return cl.translate(loc);
	}
		
	/**
	 * Translates an inclusive target coordinate, rejecting out-of-range locations
	 * and query chromosomes outside 1..25. PLUS adds loc-tStart to qStart; every
	 * other strand code subtracts that offset and asserts qStart is at least qStop.
	 * Target strand is not applied. No query-bound, span-length or overflow check is made.
	 * @param loc Target coordinate to translate
	 * @return Fresh {qChrom, qStrand, position}, or null for a rejected location/chromosome
	 */
	public int[] translate(int loc){
		if(loc<tStart || loc>tStop){return null;}
//		assert(loc>=tStart && loc<=tStop);
		if(qChrom<1 || qChrom>25){return null;}
		if(qStrand==Shared.PLUS){
			return new int[]{qChrom, qStrand, qStart+loc-tStart};
		}else{
			assert(qStart>=qStop) : this;
			return new int[]{qChrom, qStrand, qStart-(loc-tStart)};
		}
	}
	
	
	/**
	 * Tests whether inclusive range [a,b] lies within the target endpoints.
	 * Asserts b is at least a; does not check chromosome, strand or query endpoints.
	 * @param a Inclusive lower coordinate
	 * @param b Inclusive upper coordinate
	 * @return true when both endpoint comparisons succeed
	 */
	public boolean contains(int a, int b){
		assert(b>=a);
		return a>=tStart && b<=tStop;
	}
	
	
	/**
	 * Tests an inclusive target coordinate without checking chromosome or strand.
	 * @param a Coordinate to test
	 * @return true if coordinate is within tStart to tStop range, false otherwise
	 */
	public boolean contains(int a){
		return a>=tStart && a<=tStop;
	}


	/**
	 * Orders by target chromosome, start, then stop using integer subtraction.
	 * Equal target chromosomes must have equal target strands when assertions are enabled.
	 * Query fields are ignored; equal ordering keys do not imply identical mappings.
	 * @param other Segment to compare
	 * @return First nonzero target-key difference, or zero when all ordering keys match
	 */
	@Override
	public int compareTo(ChainLine other){
		int temp;
		
		temp=tChrom-other.tChrom;
		if(temp!=0){return temp;}
		
		assert(tStrand==other.tStrand);
		
		temp=tStart-other.tStart;
		if(temp!=0){return temp;}

		temp=tStop-other.tStop;
		return temp;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Target chromosome number; used as the first ordering key. */
	public int tChrom;
	/** Encoded target strand; compared by assertion, not applied by translate. */
	public byte tStrand;
	/** Inclusive lower target endpoint. */
	public int tStart;
	/** Inclusive upper target endpoint. */
	public int tStop;
	
	/** Query chromosome number; translate accepts only 1..25. */
	public int qChrom;
	/** Encoded query strand; Shared.PLUS selects addition, other codes select subtraction. */
	public byte qStrand;
	/** Query coordinate at tStart, used as the translation origin. */
	public int qStart;
	/** Query coordinate at tStop; checked only by the reverse-branch ordering assertion. */
	public int qStop;
	
}
