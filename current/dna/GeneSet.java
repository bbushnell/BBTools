package dna;
import java.util.ArrayList;

import shared.Shared;
import shared.Tools;


/**
 * Represents a collection of gene transcripts sharing the same gene ID.
 * Aggregates transcript information including chromosomal location, strand,
 * and classification flags for pseudogenes and untranslated regions.
 *
 * @author Brian Bushnell
 * @date October 27, 2010
 */
public class GeneSet implements Comparable<GeneSet>{
	
	/** Program entry point that initializes the gene ID table.
	 * @param args Command-line arguments (not used) */
	public static void main(String[] args){
		Data.getGeneIDTable();
	}
	
	/**
	 * Constructs a GeneSet from a gene name and list of transcript variants.
	 * Aggregates transcript properties and determines genomic boundaries.
	 * All genes must share the same gene ID and chromosome.
	 *
	 * @param n The gene name/identifier
	 * @param g ArrayList of Gene objects representing transcript variants
	 */
	public GeneSet(String n, ArrayList<Gene> g){
		name=n;
		id=g.get(0).id;
		genes=g;
		chrom=g.get(0).chromosome;
		transcripts=genes.size();
		assert(transcripts>0);

		byte st=-1;
		
		boolean pse=true, unt=true;
		
		for(int i=0; i<transcripts; i++){
			Gene gene=g.get(i);
			
			assert(gene.id==id) : g;
			
			pse=(pse&&gene.pseudo);
			unt=(unt&&gene.untranslated);
			minStart=min((int)gene.txStart, minStart);
			maxEnd=max((int)gene.txStop, maxEnd);
			//				assert(st==-1 || st==gene.strand) : g;
			if(st==-1){st=gene.strand;}
			else if(st!=gene.strand){st=(byte) Tools.find("?", Shared.strandCodes);}
		}
		
		pseudo=pse;
		untranslated=unt;
		
		for(Gene gene : g){
			assert(pseudo==gene.pseudo || (!pseudo && gene.pseudo)) : g;
			assert(untranslated==gene.untranslated || (!untranslated && gene.untranslated)) : g;
//			assert(untranslated==gene.untranslated) : g;
		}
		
		strand=st;
	}
	
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append(name);
		while(sb.length()<10){sb.append(' ');}
		sb.append('\t');
		sb.append(padFront(transcripts+"",2)+" transcript"+(transcripts==1 ? " " : "s"));

		sb.append("\tchr"+chrom+" ("+minStart+" - "+maxEnd+"), '"+Shared.strandCodes[strand]+"'");

		return sb.toString();
	}

	/**
	 * Pads a string with leading spaces to reach specified width.
	 * @param num The string to pad
	 * @param width Target string width
	 * @return String padded with leading spaces
	 */
	private static final String padFront(String num, int width){
		String r=num;
		while(r.length()<width){r=" "+r;}
		return r;
	}

	/**
	 * Pads a string with trailing spaces to reach specified width.
	 * @param num The string to pad
	 * @param width Target string width
	 * @return String padded with trailing spaces
	 */
	private static final String padBack(String num, int width){
		String r=num;
		while(r.length()<width){r=r+" ";}
		return r;
	}

	/** Gene name or identifier */
	public final String name;
	/** Numeric gene ID shared by all transcripts in this set */
	public final int id;
	/** Chromosome number where this gene is located */
	public final byte chrom;
	/** Strand orientation of the gene (+ or - or ? for mixed) */
	public final byte strand;
	/** List of all Gene objects (transcripts) in this gene set */
	public final ArrayList<Gene> genes;
	/** Number of transcript variants in this gene set */
	public final int transcripts;
	
	/** True if all transcripts are untranslated */
	public final boolean untranslated;
	/** True if all transcripts are psuedogenes */
	public final boolean pseudo;

	/** Leftmost (minimum) start position across all transcripts */
	public int minStart=Integer.MAX_VALUE;
	/** Rightmost (maximum) end position across all transcripts */
	public int maxEnd=0;
	

	/**
	 * Tests if a genomic position intersects this gene's boundaries.
	 * @param point Genomic coordinate to test
	 * @return true if point falls within gene boundaries
	 */
	public boolean intersects(int point){
		return point>=minStart && point<=maxEnd;
	}
	/**
	 * Tests if a genomic interval overlaps this gene's boundaries.
	 * @param point1 Start coordinate of interval
	 * @param point2 End coordinate of interval
	 * @return true if interval overlaps gene boundaries
	 */
	public boolean intersects(int point1, int point2){
		return point2>=minStart && point1<=maxEnd;
	}


	@Override
	public int compareTo(GeneSet other) {
		if(chrom!=other.chrom){
			return chrom>other.chrom ? 1 : -1;
		}
		int x=minStart<other.minStart ? -1 : minStart>other.minStart ? 1 : 0;
		if(x!=0){return x;}
		return x=name.compareTo(other.name);
	}
	
	@Override
	public boolean equals(Object other){
		return equals((GeneSet)other);
	}
	
	/**
	 * Tests equality with another GeneSet using compareTo result.
	 * @param other GeneSet to compare
	 * @return true if gene sets are equal
	 */
	public boolean equals(GeneSet other){
		return compareTo(other)==0;
	}
	
	@Override
	public int hashCode(){
		return Integer.rotateLeft(name.hashCode(), 5)^chrom;
	}
	
	/**
	 * Returns the minimum of two integers.
	 * @param x First integer
	 * @param y Second integer
	 * @return The smaller value
	 */
	private static final int min(int x, int y){return x<y ? x : y;}
	/**
	 * Returns the maximum of two integers.
	 * @param x First integer
	 * @param y Second integer
	 * @return The larger value
	 */
	private static final int max(int x, int y){return x>y ? x : y;}
	/**
	 * Returns the minimum of two bytes.
	 * @param x First byte
	 * @param y Second byte
	 * @return The smaller value
	 */
	private static final byte min(byte x, byte y){return x<y ? x : y;}
	/**
	 * Returns the maximum of two bytes.
	 * @param x First byte
	 * @param y Second byte
	 * @return The larger value
	 */
	private static final byte max(byte x, byte y){return x>y ? x : y;}
	/**
	 * Returns the minimum of two longs.
	 * @param x First long
	 * @param y Second long
	 * @return The smaller value
	 */
	private static final long min(long x, long y){return x<y ? x : y;}
	/**
	 * Returns the maximum of two longs.
	 * @param x First long
	 * @param y Second long
	 * @return The larger value
	 */
	private static final long max(long x, long y){return x>y ? x : y;}
	/**
	 * Returns the minimum of two floats.
	 * @param x First float
	 * @param y Second float
	 * @return The smaller value
	 */
	private static final float min(float x, float y){return x<y ? x : y;}
	/**
	 * Returns the maximum of two floats.
	 * @param x First float
	 * @param y Second float
	 * @return The larger value
	 */
	private static final float max(float x, float y){return x>y ? x : y;}

}