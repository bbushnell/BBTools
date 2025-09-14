package prok;

import java.util.ArrayList;
import java.util.HashMap;

import dna.AminoAcid;
import gff.GffLine;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * ORF means Open Reading Frame.
 * It starts at the first base of a start codon and ends at the last base of a stop codon.
 * The length is divisible by 3. 
 * @author Brian Bushnell
 * @date Sep 20, 2018
 *
 */
public class Orf extends PFeature {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** 
	 * Bases and coordinates are assumed to be the correct strand.
	 * Minus-strand ORFs can be flipped at the end of the constructor.
	 * @param scafName_
	 * @param start_
	 * @param stop_
	 * @param strand_
	 * @param frame_
	 * @param bases
	 * @param flip
	 */
	public Orf(String scafName_, int start_, int stop_, int strand_, int frame_, byte[] bases, boolean flip, int type_) {
		super(scafName_, start_, stop_, strand_, bases.length);
		frame=frame_;
		startCodon=getCodon(start, bases);
		stopCodon=getCodon(stop-2, bases);
		type=type_;
		
		if(flip && strand==Shared.MINUS){flip();}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Init Helpers         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Grab the codon starting at from.
	 * Assumes bases are in the correct strand
	 * @param from
	 * @param bases
	 * @return
	 */
	private static int getCodon(int from, byte[] bases){
		int codon=0;
		for(int i=0; i<3; i++){
//			assert(i+from<bases.length) : i+", "+from+", "+bases.length;
			byte b=bases[i+from];
			int x=AminoAcid.baseToNumber[b];
			codon=(codon<<2)|x;
		}
		return codon;
	}

	public float calcOrfScore(){
		return calcOrfScore(0);
	}

	/**
	 * The score of an ORF alone is a factor of the length, start score, stop score, and kmer score.
	 * The score of an ORF in the context of an overlapping gene also includes a penalty for the overlap length.
	 * @param overlap
	 * @return Calculated score
	 */
	public float calcOrfScore(int overlap){
		double a=Math.sqrt(Tools.max(f1, e1+startScore));
//		double b=Math.sqrt(f2/*Tools.max(f2, e2+stopScore)*/);//This is better, ignoring stopscore completely
		double b=Math.sqrt(Tools.max(f2, e2+0.35f*stopScore));
		double c=Tools.max(f3, e3+averageKmerScore());
		assert(a!=Double.NaN);
		assert(b!=Double.NaN);
		assert(c!=Double.NaN);
		c=4*Math.pow(c, 2.2);
		double d=(0.1*a*b*c*(Math.pow(length()-overlap, 2.5)-(overlap<1 ? 0 : Math.pow(overlap+50, 2))));//TODO: Adjust these constants
		if(d>0){d=Math.sqrt(d);}
		assert(d!=Double.NaN);
		return (float)d;
	}
	
	public float averageKmerScore(){
		return kmerScore/(length()-GeneModel.kInnerCDS-2); //This slightly affects score if kInnerCDS is changed
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean isValidPrev(Orf prev, int maxOverlap){
		if(prev.stop>=stop || prev.stop>=start+maxOverlap || prev.start>=start){return false;}
		if(prev.frame==frame && prev.strand==strand && prev.stop>=start){return false;}
		return true;
	}

	public float pathScore() {return Tools.max(pathScorePlus, pathScoreMinus);}
	public float pathScore(int prevStrand) {return prevStrand==0 ? pathScorePlus : pathScoreMinus;}

	public Orf prev(){return pathScorePlus>=pathScoreMinus ? prevPlus : prevMinus;}
	public Orf prev(int prevStrand){return prevStrand==0 ? prevPlus : prevMinus;}

	public int pathLength(int prevStrand){return prevStrand==0 ? pathLengthPlus : pathLengthMinus;}
	public int pathLength(){return pathScorePlus>=pathScoreMinus ? pathLengthPlus : pathLengthMinus;}
	
	/*--------------------------------------------------------------*/
	/*----------------           ToString           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * @param orfs A list of called features
	 * @param types Types of features to retain, e.g. "CDS,rRNA,tRNA"
	 * @return GffLines subdivided by type
	 */
	public static ArrayList<GffLine>[] toGffLinesByType(ArrayList<Orf> orfs, String types){
		String[] typeArray=types.split(",");
		@SuppressWarnings("unchecked")
		ArrayList<GffLine>[] lists=new ArrayList[typeArray.length];
		HashMap<String, ArrayList<GffLine>> map=new HashMap<String, ArrayList<GffLine>>(3*typeArray.length);
		for(int i=0; i<typeArray.length; i++){
			String type=typeArray[i];
			lists[i]=new ArrayList<GffLine>();
			map.put(type,  lists[i]);
		}
		for(Orf orf : orfs){
			String type=ProkObject.typeStrings2[orf.type];
			ArrayList<GffLine> glist=map.get(type);
			if(glist!=null) {
				glist.add(new GffLine(orf));
			}
		}
		return lists;
	}
	
	public static ArrayList<GffLine> toGffLines(ArrayList<Orf> orfs){
		if(orfs==null) {return null;}
		ArrayList<GffLine> list=new ArrayList<GffLine>(orfs.size());
		for(Orf orf : orfs){list.add(new GffLine(orf));}
		return list;
	}
	
	public String toStringFlipped(){
		if(strand==flipped()){
			return toString();
		}
		flip();
		String s=toString();
		flip();
		return s;
	}
	
	@Override
	public String toString(){
		return toGff();
	}
	
	public String toGff(){
		ByteBuilder bb=new ByteBuilder();
		appendGff(bb);
		return bb.toString();
	}
	
	public ByteBuilder appendGff(ByteBuilder bb){
		if(scafName==null){
			bb.append('.').tab();
		}else{
			for(int i=0, max=scafName.length(); i<max; i++){
				char c=scafName.charAt(i);
				if(c==' ' || c=='\t'){break;}
				bb.append(c);
			}
			bb.tab();
		}
		bb.append("BBTools").append('\t');
		bb.append(typeStrings2[type]).append('\t');
		bb.append(start+1).append('\t');
		bb.append(stop+1).append('\t');
		
		if(orfScore<0){bb.append('.').append('\t');}
		else{bb.append(orfScore, 2).append('\t');}
		
		bb.append(strand<0 ? '.' : Shared.strandCodes2[strand]).append('\t');
		
		bb.append('0').append('\t');

		//bb.append('.');
		bb.append(typeStrings[type]).append(',');
		if(type==0){
			bb.append("fr").append(frame).append(',');
		}
//		bb.append(startCodon).append(',');
//		bb.append(stopCodon).append(',');
		bb.append("startScr:").append(startScore, 3).append(',');
		bb.append("stopScr:").append(stopScore, 3).append(',');
		bb.append("innerScr:").append(averageKmerScore(), 3).append(',');
		bb.append("len:").append(length());
		if(type==0){
			bb.append(',');
			bb.append("start:").append(AminoAcid.codonToString(startCodon)).append(',');
			bb.append("stop:").append(AminoAcid.codonToString(stopCodon));
		}
		return bb;
	}
	
	public boolean isSSU(){return type==r16S || type==r18S;}
	public boolean is5S(){return type==r5S;}
	public boolean is16S(){return type==r16S;}
	public boolean is18S(){return type==r18S;}
	public boolean is23S(){return type==r23S;}
	public boolean isCDS(){return type==CDS;}
	public boolean isRRNA(){return type==r18S || type==r16S || type==r5S || type==r23S;}
	public boolean isTRNA(){return type==tRNA;}
	
	/*--------------------------------------------------------------*/
	/*----------------          Overrides           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public float score() {
		return orfScore;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final int frame;

	//These are not needed but nice for printing
	public final int startCodon;
	public final int stopCodon;
	
	public float startScore;
	public float stopScore;
	public float kmerScore;
	
	public float orfScore;

	//Path scores are for pathfinding phase
	
	public float pathScorePlus;
	public int pathLengthPlus=1;
	public Orf prevPlus;
	
	public float pathScoreMinus;
	public int pathLengthMinus=1;
	public Orf prevMinus;
	
	public final int type;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/

	/* for kinnercds=6 */ 
//	static float e1=0.1f; 
//	static float e2=-0.04f; 
//	static float e3=0.01f;//Decreasing this decreases TP but increases SNR
//	
//	static float f1=0.08f; 
//	static float f2=0.06f; 
//	static float f3=0.09f;

	/* for kinnercds=7 */ 
	static float e1=0.35f; 
	static float e2=-0.1f; 
	static float e3=-0.01f;//Decreasing this decreases TP but increases SNR
	
	static float f1=0.08f; 
	static float f2=0.02f; 
	static float f3=0.09f;
	
}
