package dna;
import java.io.Serializable;
import java.util.HashMap;

import shared.Tools;


public class Exon implements Comparable<Exon>, Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1890833345682913235L;


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
	
	public Exon(String startPoint, String endPoint, String chrom, String strnd, boolean utr_, boolean cds_){
		this(Integer.parseInt(startPoint), Integer.parseInt(endPoint), toChromosome(chrom), toStrand(strnd), utr_, cds_);
	}
	
	public Exon(int startPoint, int endPoint, String chrom, String strnd, boolean utr_, boolean cds_){
		this(startPoint, endPoint, toChromosome(chrom), toStrand(strnd), utr_, cds_);
	}
	
	public Exon(int startPoint, int endPoint, byte chrom, byte strnd, boolean utr_, boolean cds_){
		a=startPoint;
		b=endPoint;
		chromosome=chrom;
		strand=strnd;
		utr=utr_;
		cds=cds_;
	}
	
	
	
	public static Exon merge(Exon exon1, Exon exon2){
		assert(canMerge(exon1, exon2));
		return new Exon(min(exon1.a, exon2.a), max(exon1.b, exon2.b), exon1.chromosome, exon1.strand, exon1.cds||exon2.cds, exon1.utr||exon2.utr);
	}
	
	public static boolean canMerge(Exon exon1, Exon exon2){
		if(exon1.chromosome!=exon2.chromosome){return false;}
		return overlap(exon1.a, exon1.b, exon2.a, exon2.b);
	}
	
	
	public boolean intersects(int point){return point>=a && point<=b;}
	//Slow
	public boolean intersects(int a2, int b2){
		assert(a2<=b2);
		return overlap(a, b, a2, b2);
	}
	
	public boolean crosses(int a2, int b2){return (a2<a && b2>=a) || (a2<=b && b2>b);}
	public boolean contains(int a2, int b2){return (a2>=a && b2<=b);}
	
	public boolean intersectsNearby(int a, int b){
		return intersects(a-Data.NEAR, b+Data.NEAR);
	}
	
	private static boolean overlap(int a1, int b1, int a2, int b2){
		assert(a1<=b1 && a2<=b2) : a1+", "+b1+", "+a2+", "+b2;
		return a2<=b1 && b2>=a1;
	}
	
	public int distToSpliceSite(int x, int y){
		int distA=distToPoint(x, y, a);
		int distB=distToPoint(x, y, b);
		return min(distA, distB);
	}
	
	public static int distToPoint(int x, int y, int point){
		assert(x<=y);
		if(y<=point){return point-y;}
		if(x>=point){return x-point;}
		return 0;
	}
	
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
	
	public static byte toChromosome(String s){
		int i=0;
//		System.out.println(s);
		while(!Tools.isDigit(s.charAt(i))){i++;}
		return Byte.parseByte(s.substring(i));
	}
	
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
	
	public boolean equals(Exon other){
		return a==other.a && b==other.b && chromosome==other.chromosome && strand==other.strand && utr==other.utr && cds==other.cds;
	}
	
	@Override
	public int hashCode(){
		int xor=a^(Integer.rotateLeft(b, 16));
		xor^=Integer.rotateRight(chromosome, 6);
		return xor;
	}
	
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	
	public final int a;
	public final int b;
	public final boolean utr;
	public final boolean cds;
	public final byte chromosome;
	public final byte strand;
	
	public static final HashMap<Exon,Exon> table=new HashMap<Exon,Exon>(65536);
}
