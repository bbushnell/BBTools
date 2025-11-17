package fileIO;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import dna.Gene;
import shared.Shared;

/** For loading UCSC .chain files that convert one genome build to another. */
public class ChainBlock implements Comparable<ChainBlock>{
	
	
	/** Test method that loads and displays chain lines for chromosomes 1-22.
	 * @param args Command-line arguments, expects chain file path as first argument */
	public static void main(String args[]){
		ChainLine[][] lines=loadChainLines(args[0]);
		for(int i=1; i<=22; i++){
			for(ChainLine line : lines[i]){
				System.out.println(line);
			}
			System.out.println();
		}
	}


	/**
	 * Constructs a ChainBlock from parsed chain file lines.
	 * Parses the header line containing chain metadata and subsequent alignment chunks.
	 * The first element of the list must be the chain header line starting with "chain".
	 * @param list List of String arrays representing parsed chain file lines
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
	 * Converts chromosome name string to numeric chromosome code.
	 * Falls back to "U" (unknown) chromosome code if conversion fails.
	 * @param s Chromosome name string
	 * @return Numeric chromosome code
	 */
	private static int toChromosome(String s){
		int result;
		try{
			result=Gene.toChromosome(s);
		}catch(Exception e){
			result=Gene.toChromosome("U");
		}
		return result;
	}
	
	/**
	 * Converts this ChainBlock into an array of ChainLine objects.
	 * Processes alignment chunks based on query strand orientation.
	 * For plus strand queries, coordinates advance normally.
	 * For minus strand queries, query coordinates are processed in reverse.
	 *
	 * @return Array of ChainLine objects representing alignment segments
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
	
	
	/**
	 * Loads chain file and converts all blocks to organized ChainLine arrays.
	 * Loads chain blocks, organizes them by chromosome, and converts to sorted ChainLines.
	 * @param fname Chain file path
	 * @return Two-dimensional array of ChainLines organized by chromosome
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
	 * Loads and parses a chain file into ChainBlock objects.
	 * Reads the entire file, splits lines into tokens, and groups lines into chain blocks.
	 * Each block ends with a single-element line.
	 *
	 * @param fname Chain file path to load
	 * @return ArrayList of parsed ChainBlock objects sorted by position
	 */
	public static ArrayList<ChainBlock> loadChainBlocks(String fname){
		TextFile tf=new TextFile(fname, false);
		String[] lines=tf.toStringLines();
		tf.close();
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
		Shared.sort(out);
		return out;
	}
	
	
	/**
	 * Organizes ChainBlocks into arrays by target chromosome.
	 * Creates a two-dimensional array where each chromosome has its own ChainBlock array.
	 * @param list List of ChainBlocks to organize
	 * @return Two-dimensional array of ChainBlocks organized by chromosome
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
	

	@Override
	public int compareTo(ChainBlock other) {
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
	
	
	/** Chain alignment score */
	public long score;
	/** Target sequence name (chromosome) */
	public String tName;
	/** Target chromosome numeric code */
	public int tChrom;
	/** Target sequence size */
	public int tSize;
	/** Target sequence strand orientation */
	public byte tStrand;
	/** Target sequence alignment start position */
	public int tStart;
	/** Target sequence alignment end position */
	public int tStop;

	/** Query sequence name (chromosome) */
	public String qName;
	/** Query chromosome numeric code */
	public int qChrom;
	/** Query sequence size */
	public int qSize;
	/** Query sequence strand orientation */
	public byte qStrand;
	/** Query sequence alignment start position */
	public int qStart;
	/** Query sequence alignment end position */
	public int qStop;
	
	/** Unique identifier for this chain */
	public int chainID;
	
	/** Alignment chunk data arrays containing match lengths and gap sizes */
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
