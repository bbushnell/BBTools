package dna;

import java.io.Serializable;
import java.util.HashMap;

import shared.Tools;
/**
 * An immutable genomic interval [a,b] (inclusive) on a given chromosome and strand,
 * tagged as UTR and/or CDS. Provides interval algebra (intersect/overlap/crosses/
 * contains) and distance-to-splice-site queries. Coordinates are inclusive, so
 * length()=b-a+1. Gene's exon parsers explicitly intern instances in {@link #table};
 * constructors do not intern. Instances are immutable, but the shared intern table
 * is mutable and requires external synchronization for concurrent modification.
 * @author Brian Bushnell
 */
public class Exon implements Comparable<Exon>, Serializable{

	/** Serialization identity retained for stored gene annotations. */
	private static final long serialVersionUID=1890833345682913235L;

	/** Creates a placeholder with both coordinates, chromosome and strand set to -1.
	 * Its equal endpoints give length one; it is not an empty-interval representation. */
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

	/** Parses inclusive integer endpoints, a numeric chromosome and a strand symbol.
	 * Coordinates are stored without a one-based/zero-based conversion. */
	public Exon(String startPoint, String endPoint, String chrom, String strnd, boolean utr_, boolean cds_){
		this(Integer.parseInt(startPoint), Integer.parseInt(endPoint), toChromosome(chrom), toStrand(strnd), utr_, cds_);
	}

	/** Stores inclusive endpoints and parses the chromosome and strand strings. */
	public Exon(int startPoint, int endPoint, String chrom, String strnd, boolean utr_, boolean cds_){
		this(startPoint, endPoint, toChromosome(chrom), toStrand(strnd), utr_, cds_);
	}

	/** Stores all fields as supplied; callers must provide ordered inclusive endpoints.
	 * UTR and CDS flags may both be true when the exon crosses a coding boundary. */
	public Exon(int startPoint, int endPoint, byte chrom, byte strnd, boolean utr_, boolean cds_){
		a=startPoint;
		b=endPoint;
		chromosome=chrom;
		strand=strnd;
		utr=utr_;
		cds=cds_;
	}

	/** Returns the union of overlapping same-chromosome exons, OR-ing their UTR/CDS flags.
	 * Requires {@link #canMerge(Exon, Exon)}; retains exon1's strand even if strands differ.
	 * Neither input nor the intern table is modified. */
	public static Exon merge(Exon exon1, Exon exon2){
		assert(canMerge(exon1, exon2));
		return new Exon(min(exon1.a, exon2.a), max(exon1.b, exon2.b), exon1.chromosome, exon1.strand, exon1.utr||exon2.utr, exon1.cds||exon2.cds);
	}

	/** Tests chromosome equality and inclusive overlap; ignores strand and annotation flags.
	 * Adjacent, nonoverlapping intervals are not mergeable. Both inputs must be nonnull. */
	public static boolean canMerge(Exon exon1, Exon exon2){
		if(exon1.chromosome!=exon2.chromosome){return false;}
		return overlap(exon1.a, exon1.b, exon2.a, exon2.b);
	}

	/** Whether point is between the inclusive endpoints. */
	public boolean intersects(int point){return point>=a && point<=b;}
	//Slow
	/** Whether the ordered inclusive query interval shares at least one coordinate. */
	public boolean intersects(int a2, int b2){
		assert(a2<=b2);
		return overlap(a, b, a2, b2);
	}

	/** Whether an ordered query spans an endpoint and extends strictly outside that side.
	 * A query equal to this exon does not cross; a strictly enclosing query does. */
	public boolean crosses(int a2, int b2){return (a2<a && b2>=a) || (a2<=b && b2>b);}
	/** Whether the ordered inclusive query is entirely within this exon. */
	public boolean contains(int a2, int b2){return (a2>=a && b2<=b);}

	/** Tests overlap after extending the query by Data.NEAR on both sides.
	 * The expanded endpoints must fit in int arithmetic. */
	public boolean intersectsNearby(int a, int b){
		return intersects(a-Data.NEAR, b+Data.NEAR);
	}

	/** Tests inclusive overlap of two ordered intervals. */
	private static boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1;
	}

	/** Minimum coordinate distance from the ordered query to either exon endpoint.
	 * Zero means the query contains an endpoint; strand does not affect the result. */
	public int distToSpliceSite(int x, int y){
		int distA=distToPoint(x, y, a);
		int distB=distToPoint(x, y, b);
		return min(distA, distB);
	}

	/** Coordinate distance from an ordered inclusive interval to point, or zero inside it.
	 * Endpoint differences must fit in int arithmetic. */
	public static int distToPoint(int x, int y, int point){
		assert(x<=y);
		if(y<=point){return point-y;}
		if(x>=point){return x-point;}
		return 0;
	}

	/** Maps '+' to 0, '-' to 1 and '?' to 2.
	 * Other values fail an assertion when enabled, otherwise they also return 2. */
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

	/** Parses a numeric chromosome id from a string, skipping any leading non-digit prefix (e.g. "chr7" -> 7).
	 * Requires at least one digit; a purely non-numeric name ("chrX"/"chrY"/"chrM") has no numeric id and is
	 * unsupported (BBTools is chromosome-number-centric; the human-era X/Y/M name map was retired). */
	public static byte toChromosome(String s){
		int i=0;
//		System.out.println(s);
		while(!Tools.isDigit(s.charAt(i))){i++;}
		return Byte.parseByte(s.substring(i));
	}

	/** Inclusive length b-a+1; requires a positive result fitting in int arithmetic. */
	public int length(){
		int r=(int)(b-a+1);
		assert(r>0);
		return r;
	}

	/** Returns string representation of this exon showing chromosome, coordinates, and length.
	 * @return String in format "(chrN, start - end, len length)" */
	@Override
	public String toString(){
//		return "(chr"+chromosome+","+(strand==0 ? "+" : "-")+","+a+"~"+b+")";
		return "(chr"+chromosome+", "+a+" - "+b+", len "+length()+")";
	}

	/** Orders by chromosome, start, end, strand, UTR and CDS; true flags sort first. */
	@Override
	public int compareTo(Exon other){
		if(chromosome<other.chromosome){return -1;}
		if(chromosome>other.chromosome){return 1;}

		if(a<other.a){return -1;}
		if(a>other.a){return 1;}

		if(b<other.b){return -1;}
		if(b>other.b){return 1;}

		if(strand<other.strand){return -1;}
		if(strand>other.strand){return 1;}

		if(utr && !other.utr){return -1;}
		if(!utr && other.utr){return 1;}

		if(cds && !other.cds){return -1;}
		if(!cds && other.cds){return 1;}

		return 0;
	}

	/** Delegates after a direct cast; unrelated types throw ClassCastException.
	 * The current typed implementation also throws for null. */
	@Override
	public boolean equals(Object other){
		return equals((Exon)other);
	}

	/** Compares all six fields; currently requires a nonnull argument. */
	public boolean equals(Exon other){
		return a==other.a && b==other.b && chromosome==other.chromosome && strand==other.strand && utr==other.utr && cds==other.cds;
	}

	/** Hashes coordinates and chromosome; omitted flags/strand can cause allowed collisions. */
	@Override
	public int hashCode(){
		int xor=a^(Integer.rotateLeft(b, 16));
		xor^=Integer.rotateRight(chromosome, 6);
		return xor;
	}

	/** Smaller of two coordinates. */
	private static final int min(int x, int y){return x<y ? x : y;}
	/** Larger of two coordinates. */
	private static final int max(int x, int y){return x>y ? x : y;}

	/** Inclusive start coordinate. */
	public final int a;
	/** Inclusive end coordinate. */
	public final int b;
	/** Whether any portion is annotated as untranslated. */
	public final boolean utr;
	/** Whether any portion overlaps the annotated coding interval. */
	public final boolean cds;
	/** Legacy signed-byte chromosome identifier. */
	public final byte chromosome;
	/** Strand code: 0 plus, 1 minus, 2 unknown; -1 for the default placeholder. */
	public final byte strand;

	/** Shared mutable intern pool used explicitly by Gene parsers; not thread-safe.
	 * Entries remain retained until callers remove them or clear the table. */
	public static final HashMap<Exon, Exon> table=new HashMap<Exon, Exon>(65536);
}
