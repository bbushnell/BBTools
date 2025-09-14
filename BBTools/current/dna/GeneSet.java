package dna;
import java.util.ArrayList;

import shared.Shared;
import shared.Tools;


public class GeneSet implements Comparable<GeneSet>{
	
	public static void main(String[] args){
		Data.getGeneIDTable();
	}
	
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

	private static final String padFront(String num, int width){
		String r=num;
		while(r.length()<width){r=" "+r;}
		return r;
	}

	private static final String padBack(String num, int width){
		String r=num;
		while(r.length()<width){r=r+" ";}
		return r;
	}

	public final String name;
	public final int id;
	public final byte chrom;
	public final byte strand;
	public final ArrayList<Gene> genes;
	public final int transcripts;
	
	/** True if all transcripts are untranslated */
	public final boolean untranslated;
	/** True if all transcripts are psuedogenes */
	public final boolean pseudo;

	public int minStart=Integer.MAX_VALUE;
	public int maxEnd=0;
	

	public boolean intersects(int point){
		return point>=minStart && point<=maxEnd;
	}
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
	
	public boolean equals(GeneSet other){
		return compareTo(other)==0;
	}
	
	@Override
	public int hashCode(){
		return Integer.rotateLeft(name.hashCode(), 5)^chrom;
	}
	
	private static final int min(int x, int y){return x<y ? x : y;}
	private static final int max(int x, int y){return x>y ? x : y;}
	private static final byte min(byte x, byte y){return x<y ? x : y;}
	private static final byte max(byte x, byte y){return x>y ? x : y;}
	private static final long min(long x, long y){return x<y ? x : y;}
	private static final long max(long x, long y){return x>y ? x : y;}
	private static final float min(float x, float y){return x<y ? x : y;}
	private static final float max(float x, float y){return x>y ? x : y;}

}