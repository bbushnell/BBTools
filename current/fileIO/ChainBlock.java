package fileIO;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import dna.Gene;
import shared.Shared;

/**
 * Stores a chain header and its ungapped alignment chunks for legacy build conversion.
 * Header coordinates are parsed literally; {@link #toLines()} constructs inclusive
 * segments using the arithmetic documented there. The parser does not validate the
 * complete external format or reconcile header spans with chunk lengths and gaps.
 * Fields and chunk arrays are mutable, without synchronization. Chromosome grouping
 * uses Gene.chromCodes indices, while name parsing follows Gene's numeric-only policy.
 *
 * @author Brian Bushnell
 */
public class ChainBlock implements Comparable<ChainBlock>{

	/*--------------------------------------------------------------*/
	/*----------------            Startup           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Prints the loaded segments for chromosome indices 1 through 22.
	 * Prints one blank separator per index, including uncovered chromosomes.
	 * @param args Positional input filename at index zero
	 */
	public static void main(String[] args){
		ChainLine[][] lines=loadChainLines(args[0]);
		for(int i=1; i<=22; i++){
			//FIO-055: the loader represents an uncovered chromosome with a null segment array.
			final ChainLine[] array=lines[i];
			if(array!=null){
				for(ChainLine line : array){
					System.out.println(line);
				}
			}
			System.out.println();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses one tokenized header followed by chunk rows into newly allocated arrays.
	 * Asserts the header marker and three-field/nonfinal versus one-field/final row
	 * shapes. Header field count, positive lengths and span consistency are not checked.
	 * A header alone creates an empty chunk array. Input lists/row arrays are not retained;
	 * immutable name strings are retained directly. Unsupported names may fail conversion.
	 * @param list Nonempty list: 13-field header, then {size, target gap, query gap}
	 * rows ending with a {size} row
	 */
	public ChainBlock(List<String[]> list){
		
		String[] head=list.get(0);
		assert("chain".equals(head[0]));
		
		score=Long.parseLong(head[1]);
		
		tName=head[2];
		tChrom=toChromosome(head[2]);
		tSize=Integer.parseInt(head[3]);
		tStrand=Gene.toStrand(head[4]);
		tStart=Integer.parseInt(head[5]);
		tStop=Integer.parseInt(head[6]);

		qName=head[7];
		qChrom=toChromosome(head[7]);
		qSize=Integer.parseInt(head[8]);
		qStrand=Gene.toStrand(head[9]);
		qStart=Integer.parseInt(head[10]);
		qStop=Integer.parseInt(head[11]);
		
		chainID=Integer.parseInt(head[12]);
		
		chunks=new int[list.size()-1][];
		for(int i=1; i<list.size(); i++){
			String[] line=list.get(i);
			assert((i==list.size()-1) == (line.length==1));
			assert((i!=list.size()-1) == (line.length==3));
			chunks[i-1]=new int[line.length];
			for(int j=0; j<line.length; j++){
				chunks[i-1][j]=Integer.parseInt(line[j]);
			}
		}
		
	}
	
	/**
	 * Uses Gene's numeric chromosome parser, then tries the legacy U fallback on Exception.
	 * The fallback is also unsupported by that parser; assertion failures are not caught.
	 * @param s Chromosome token from the header
	 * @return Parsed numeric chromosome ID when conversion succeeds
	 */
	private static int toChromosome(String s){
		int result;
		try{
			result=Gene.toChromosome(s);
		}catch(Exception e){
			//TODO: Probable bug - FIO-051: Gene.toChromosome now accepts numeric names only, so
			//this legacy unknown-name fallback also throws. Reconcile name policy before repairing.
			result=Gene.toChromosome("U");
		}
		return result;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Expands chunks into fresh inclusive segments in input chunk order.
	 * Target traversal starts at tStart and increases by each size and following gap.
	 * PLUS query traversal starts at qStart and increases; every other strand starts
	 * at qStop-1 and decreases by each size and query gap. No qSize-based coordinate
	 * conversion is performed. tStrand is copied to the result without changing traversal.
	 * Does not check header endpoints against the final segment or validate lengths,
	 * gaps, array shapes, coordinate overflow or external reverse-strand conventions.
	 * @return New array with one new ChainLine per chunk, empty when there are no chunks
	 */
	public ChainLine[] toLines(){
		ChainLine[] out=new ChainLine[chunks.length];
		
		if(qStrand==Shared.PLUS){

			int tloc=tStart, qloc=qStart;
			for(int i=0; i<chunks.length; i++){
				int[] chunk=chunks[i];
				int tloc2=tloc+chunk[0]-1, qloc2=qloc+chunk[0]-1;
				out[i]=new ChainLine(tChrom, tStrand, tloc, tloc2, qChrom, qStrand, qloc, qloc2);
				if(chunk.length>1){
					tloc=tloc2+chunk[1]+1;
					qloc=qloc2+chunk[2]+1;
				}
			}
		}else{

			int tloc=tStart, qloc=qStop-1;
			for(int i=0; i<chunks.length; i++){
				int[] chunk=chunks[i];
				int tloc2=tloc+chunk[0]-1, qloc2=qloc-chunk[0]+1;
				out[i]=new ChainLine(tChrom, tStrand, tloc, tloc2, qChrom, qStrand, qloc, qloc2);
				if(chunk.length>1){
					tloc=tloc2+chunk[1]+1;
					qloc=qloc2-chunk[2]-1;
				}
			}
		}
		
		return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads blocks, groups by target chromosome and sorts each nonempty segment array.
	 * Empty chromosomes remain null. Merges all blocks without score/chain-ID selection
	 * or overlap validation. ChainLine search resolves overlaps by earliest source start,
	 * then earliest destination at that start; sorting here retains all alternatives.
	 * @param fname Input name interpreted by TextFile with subprocesses disabled
	 * @return New chromosome-indexed array, sized to Gene.chromCodes.length
	 */
	public static ChainLine[][] loadChainLines(String fname){
		ArrayList<ChainBlock> list=loadChainBlocks(fname);
		ChainBlock[][] blocks=splitChain(list);
		ChainLine[][] out=new ChainLine[blocks.length][];
		ArrayList<ChainLine> temp=new ArrayList<ChainLine>();
		for(int chrom=0; chrom<blocks.length; chrom++){
			temp.clear();
			ChainBlock[] cblocks=blocks[chrom];
			if(cblocks.length>0){
				for(ChainBlock block : cblocks){
					ChainLine[] blines=block.toLines();
					for(ChainLine line : blines){
						temp.add(line);
					}
				}
			}
			if(temp.size()>0){
				out[chrom]=temp.toArray(new ChainLine[temp.size()]);
				Arrays.sort(out[chrom]);
			}
		}
		return out;
	}
	
	
	/**
	 * Reads all TextFile logical lines, trims/splits whitespace and constructs blocks.
	 * Every one-token row terminates the current block, which is then parsed; blank
	 * handling comes from TextFile and there is no separate comment-line filter.
	 * Sorts complete blocks by their natural order. Rejects a nonempty unfinished
	 * tail rather than returning only the preceding complete blocks.
	 * Attempts reader closure after collection, including on read failure, before parsing.
	 * Reported close errors reject the result; a primary read failure remains primary.
	 * @param fname Input name interpreted by TextFile with subprocesses disabled
	 * @return New mutable list of parsed, sorted complete blocks
	 * @throws IllegalArgumentException If rows remain after the last block terminator
	 */
	public static ArrayList<ChainBlock> loadChainBlocks(String fname){
		//FIO-056: the static helper owns cleanup and preserves an existing collection failure.
		String[] lines=TextFile.toStringLines(fname);
		String[][] text=TextFile.doublesplitWhitespace(lines, true);
		
		ArrayList<ChainBlock> out=new ArrayList<ChainBlock>();
		ArrayList<String[]> current=new ArrayList<String[]>(40);
		for(int i=0; i<text.length; i++){
			String[] line=text[i];
			current.add(line);
			if(line.length==1){
				out.add(new ChainBlock(current));
				current.clear();
			}
		}
		//FIO-054: only a one-field size row finishes a block; leftover rows mean incomplete input.
		if(!current.isEmpty()){
			throw new IllegalArgumentException("Incomplete chain block at end of "+fname+
				"; buffered rows: "+current.size());
		}
		Shared.sort(out);
		return out;
	}
	
	
	/**
	 * Groups existing block references by tChrom, preserving input order within a group.
	 * Allocates every group, including empty arrays; does not copy blocks or sort them.
	 * @param list Nonnull blocks whose target IDs index Gene.chromCodes (currently 0..26)
	 * @return New outer/group arrays with the same mutable block objects
	 */
	public static ChainBlock[][] splitChain(ArrayList<ChainBlock> list){
		int[] size=new int[Gene.chromCodes.length];
		
		for(ChainBlock cb : list){size[cb.tChrom]++;}
		
		ChainBlock[][] out=new ChainBlock[size.length][];
		for(int i=0; i<out.length; i++){out[i]=new ChainBlock[size[i]];}
		
		Arrays.fill(size, 0);
		for(ChainBlock cb : list){
			out[cb.tChrom][size[cb.tChrom]]=cb;
			size[cb.tChrom]++;
		}
		
		return out;
	}
	

	/**
	 * Orders by target chromosome, target name, target start, then target stop.
	 * Integer keys use subtraction. Equal chromosome/name keys assert equal strands;
	 * score, query fields and chain ID do not break ties.
	 * @param other Block with a nonnull target name
	 * @return First nonzero key comparison, or zero when all ordering keys match
	 */
	@Override
	public int compareTo(ChainBlock other){
		int temp;
		
		temp=tChrom-other.tChrom;
		if(temp!=0){return temp;}
		
		temp=tName.compareTo(other.tName);
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

	/** Header score; not used to select or order mappings. */
	public long score;
	/** Literal target name; retained even when multiple names parse to the same ID. */
	public String tName;
	/** Parsed numeric target ID; grouping requires an index into Gene.chromCodes. */
	public int tChrom;
	/** Header target sequence size; retained without span validation. */
	public int tSize;
	/** Encoded target strand; copied to segments without changing target traversal. */
	public byte tStrand;
	/** Header target start; origin of ascending chunk traversal. */
	public int tStart;
	/** Header target end; stored literally and not reconciled with chunks. */
	public int tStop;

	/** Literal query name. */
	public String qName;
	/** Parsed numeric query ID; ChainLine.translate separately accepts only 1..25. */
	public int qChrom;
	/** Header query sequence size; not used by toLines. */
	public int qSize;
	/** Encoded query strand; PLUS increases query coordinates, other codes decrease them. */
	public byte qStrand;
	/** Header query start; used as the PLUS traversal origin. */
	public int qStart;
	/** Header query end; qStop-1 is the non-PLUS traversal origin. */
	public int qStop;
	
	/** Parsed header chain ID; not used for ordering or mapping selection. */
	public int chainID;
	
	/** Mutable {size, target gap, query gap} rows followed by a final {size} row. */
	public int[][] chunks;
	
	//chain 3303 chr1 247249719 + 13192499 13192587 chr1 249250621 - 236203315 236203403 109
	
//    *   score -- chain score
//    * tName -- chromosome (reference sequence)
//    * tSize -- chromosome size (reference sequence)
//    * tStrand -- strand (reference sequence)
//    * tStart -- alignment start position (reference sequence)
//    * tEnd -- alignment end position (reference sequence)
//    * qName -- chromosome (query sequence)
//    * qSize -- chromosome size (query sequence)
//    * qStrand -- strand (query sequence)
//    * qStart -- alignment start position (query sequence)
//    * qEnd -- alignment end position (query sequence)
//    * id -- chain ID
	
}
