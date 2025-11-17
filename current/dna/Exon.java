package dna;
import java.io.Serializable;
import java.util.HashMap;

import shared.Tools;


/**
 * Represents a genomic exon with coordinate information and region type flags.
 * Contains start and end positions, chromosome and strand information, and
 * flags indicating whether the exon contains UTR or CDS regions.
 * Supports merging, intersection testing, and distance calculations.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class Exon implements Comparable<Exon>, Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1890833345682913235L;


	/** Creates an uninitialized exon with default values.
	 * Start and end coordinates are set to -1, flags are false, and chromosome/strand are -1. */
	public Exon(){
		a=-1;
		b=-1;
		utr=false;
		cds=false;
		chromosome=-1;
		strand=-1;
	}
	
//	public Exon(String startPoint, String endPoint, String chrom){
//		this(startPoint, endPoint, chrom, "?");
//	}
//
//	public Exon(int startPoint, int endPoint, String chrom){
//		this(startPoint, endPoint, chrom, "?");
//	}
//
//	public Exon(int startPoint, int endPoint, byte chrom){
//		this(startPoint, endPoint, chrom, (byte)2);
//	}
	
	/**
	 * Creates an exon from string representations of coordinates and identifiers.
	 *
	 * @param startPoint Start coordinate as string
	 * @param endPoint End coordinate as string
	 * @param chrom Chromosome identifier as string
	 * @param strnd Strand identifier as string ("+", "-", or "?")
	 * @param utr_ Whether this exon contains UTR regions
	 * @param cds_ Whether this exon contains CDS regions
	 */
	public Exon(String startPoint, String endPoint, String chrom, String strnd, boolean utr_, boolean cds_){
		this(Integer.parseInt(startPoint), Integer.parseInt(endPoint), toChromosome(chrom), toStrand(strnd), utr_, cds_);
	}
	
	/**
	 * Creates an exon with integer coordinates and string identifiers.
	 *
	 * @param startPoint Start coordinate
	 * @param endPoint End coordinate
	 * @param chrom Chromosome identifier as string
	 * @param strnd Strand identifier as string ("+", "-", or "?")
	 * @param utr_ Whether this exon contains UTR regions
	 * @param cds_ Whether this exon contains CDS regions
	 */
	public Exon(int startPoint, int endPoint, String chrom, String strnd, boolean utr_, boolean cds_){
		this(startPoint, endPoint, toChromosome(chrom), toStrand(strnd), utr_, cds_);
	}
	
	/**
	 * Creates an exon with integer coordinates and byte identifiers.
	 *
	 * @param startPoint Start coordinate
	 * @param endPoint End coordinate
	 * @param chrom Chromosome identifier as byte
	 * @param strnd Strand identifier as byte (0=+, 1=-, 2=?)
	 * @param utr_ Whether this exon contains UTR regions
	 * @param cds_ Whether this exon contains CDS regions
	 */
	public Exon(int startPoint, int endPoint, byte chrom, byte strnd, boolean utr_, boolean cds_){
		a=startPoint;
		b=endPoint;
		chromosome=chrom;
		strand=strnd;
		utr=utr_;
		cds=cds_;
	}
	
	
	
	/**
	 * Merges two overlapping exons into a single exon spanning both ranges.
	 * The resulting exon covers the minimum start to maximum end coordinates.
	 * UTR and CDS flags are combined with logical OR.
	 *
	 * @param exon1 First exon to merge
	 * @param exon2 Second exon to merge
	 * @return New merged exon covering both input ranges
	 */
	public static Exon merge(Exon exon1, Exon exon2){
		assert(canMerge(exon1, exon2));
		return new Exon(min(exon1.a, exon2.a), max(exon1.b, exon2.b), exon1.chromosome, exon1.strand, exon1.cds||exon2.cds, exon1.utr||exon2.utr);
	}
	
	/**
	 * Tests whether two exons can be merged.
	 * Exons can be merged if they are on the same chromosome and overlap.
	 *
	 * @param exon1 First exon to test
	 * @param exon2 Second exon to test
	 * @return true if exons can be merged, false otherwise
	 */
	public static boolean canMerge(Exon exon1, Exon exon2){
		if(exon1.chromosome!=exon2.chromosome){return false;}
		return overlap(exon1.a, exon1.b, exon2.a, exon2.b);
	}
	
	
	/**
	 * Tests whether a point falls within this exon's coordinate range.
	 * @param point Coordinate to test
	 * @return true if point is within exon bounds (inclusive)
	 */
	public boolean intersects(int point){return point>=a && point<=b;}
	//Slow
	/**
	 * Tests whether a coordinate range overlaps with this exon.
	 * @param a2 Start coordinate of range to test
	 * @param b2 End coordinate of range to test
	 * @return true if ranges overlap, false otherwise
	 */
	public boolean intersects(int a2, int b2){
		assert(a2<=b2);
		return overlap(a, b, a2, b2);
	}
	
	/**
	 * Tests whether a coordinate range crosses the boundaries of this exon.
	 * Returns true if the range starts before and ends within, or starts within and ends after.
	 *
	 * @param a2 Start coordinate of range to test
	 * @param b2 End coordinate of range to test
	 * @return true if range crosses exon boundaries
	 */
	public boolean crosses(int a2, int b2){return (a2<a && b2>=a) || (a2<=b && b2>b);}
	/**
	 * Tests whether this exon completely contains the specified coordinate range.
	 * @param a2 Start coordinate of range to test
	 * @param b2 End coordinate of range to test
	 * @return true if this exon fully contains the specified range
	 */
	public boolean contains(int a2, int b2){return (a2>=a && b2<=b);}
	
	/**
	 * Tests whether a coordinate range intersects with this exon including nearby regions.
	 * Expands the test range by Data.NEAR on both sides before testing intersection.
	 *
	 * @param a Start coordinate of range to test
	 * @param b End coordinate of range to test
	 * @return true if range intersects exon or nearby regions
	 */
	public boolean intersectsNearby(int a, int b){
		return intersects(a-Data.NEAR, b+Data.NEAR);
	}
	
	/**
	 * Tests whether two coordinate ranges overlap.
	 *
	 * @param a1 Start of first range
	 * @param b1 End of first range
	 * @param a2 Start of second range
	 * @param b2 End of second range
	 * @return true if ranges overlap, false otherwise
	 */
	private static boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1;
	}
	
	/**
	 * Calculates the minimum distance from a coordinate range to either splice site of this exon.
	 * Tests distance to both start and end coordinates of the exon.
	 *
	 * @param x Start coordinate of range
	 * @param y End coordinate of range
	 * @return Minimum distance to either splice site
	 */
	public int distToSpliceSite(int x, int y){
		int distA=distToPoint(x, y, a);
		int distB=distToPoint(x, y, b);
		return min(distA, distB);
	}
	
	/**
	 * Calculates the minimum distance from a coordinate range to a specific point.
	 * Returns 0 if the point falls within the range.
	 *
	 * @param x Start coordinate of range
	 * @param y End coordinate of range
	 * @param point Target point coordinate
	 * @return Minimum distance from range to point, or 0 if point is within range
	 */
	public static int distToPoint(int x, int y, int point){
		assert(x<=y);
		if(y<=point){return point-y;}
		if(x>=point){return x-point;}
		return 0;
	}
	
	/**
	 * Converts string strand representation to byte encoding.
	 * @param s Strand string: "+" for forward, "-" for reverse, "?" for unknown
	 * @return Byte encoding: 0 for "+", 1 for "-", 2 for "?"
	 */
	public static byte toStrand(String s){
		byte r=2;
		if("-".equals(s)){
			r=1;
		}else if("+".equals(s)){
			r=0;
		}else{
			assert("?".equals(s));
		}
		return r;
	}
	
	/**
	 * Extracts chromosome number from string identifier.
	 * Skips non-digit characters and parses the numeric portion.
	 * @param s Chromosome identifier string (e.g., "chr1", "chromosome2")
	 * @return Chromosome number as byte
	 */
	public static byte toChromosome(String s){
		int i=0;
//		System.out.println(s);
		while(!Tools.isDigit(s.charAt(i))){i++;}
		return Byte.parseByte(s.substring(i));
	}
	
	/** Calculates the length of this exon in base pairs.
	 * @return Length as (end - start + 1) */
	public int length(){
		int r=(int)(b-a+1);
		assert(r>0);
		return r;
	}
	
	@Override
	public String toString(){
//		return "(chr"+chromosome+","+(strand==0 ? "+" : "-")+","+a+"~"+b+")";
		return "(chr"+chromosome+", "+a+" - "+b+", len "+length()+")";
	}
	
	@Override
	public int compareTo(Exon other){
		if(chromosome<other.chromosome){return -1;}
		if(chromosome>other.chromosome){return 1;}
		
		if(a<other.a){return -1;}
		if(a>other.a){return 1;}

		if(b<other.a){return -1;}
		if(b>other.a){return 1;}

		if(strand<other.strand){return -1;}
		if(strand>other.strand){return 1;}

		if(utr && !other.utr){return -1;}
		if(!utr && other.utr){return 1;}
		
		if(cds && !other.cds){return -1;}
		if(!cds && other.cds){return 1;}
		
		return 0;
	}
	
	@Override
	public boolean equals(Object other){
		return equals((Exon)other);
	}
	
	/**
	 * Tests equality with another exon.
	 * Exons are equal if all coordinates, chromosome, strand, and flags match.
	 * @param other Exon to compare against
	 * @return true if exons have identical properties
	 */
	public boolean equals(Exon other){
		return a==other.a && b==other.b && chromosome==other.chromosome && strand==other.strand && utr==other.utr && cds==other.cds;
	}
	
	@Override
	public int hashCode(){
		int xor=a^(Integer.rotateLeft(b, 16));
		xor^=Integer.rotateRight(chromosome, 6);
		return xor;
	}
	
	
	/**
	 * Returns the smaller of two integers.
	 * @param x First integer
	 * @param y Second integer
	 * @return Minimum value
	 */
	private static final int min(int x, int y){return x<y ? x : y;}
	/**
	 * Returns the larger of two integers.
	 * @param x First integer
	 * @param y Second integer
	 * @return Maximum value
	 */
	private static final int max(int x, int y){return x>y ? x : y;}
	
	/** Start coordinate of the exon (inclusive) */
	public final int a;
	/** End coordinate of the exon (inclusive) */
	public final int b;
	/** Whether this exon contains untranslated region (UTR) sequence */
	public final boolean utr;
	/** Whether this exon contains coding sequence (CDS) regions */
	public final boolean cds;
	/** Chromosome identifier as byte value */
	public final byte chromosome;
	/**
	 * Strand identifier: 0 for forward (+), 1 for reverse (-), 2 for unknown (?)
	 */
	public final byte strand;
	
	/** Hash table for exon lookup and deduplication */
	public static final HashMap<Exon,Exon> table=new HashMap<Exon,Exon>(65536);
}
