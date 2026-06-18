package prok;

import java.util.ArrayList;
import java.util.HashMap;

import dna.AminoAcid;
import gff.GffLine;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * A scored open reading frame (gene candidate) on a scaffold: coordinates, reading frame, packed
 * start/stop codons, the component scores (start/stop/kmer) and the combined orfScore, plus the
 * per-strand path-finding fields (pathScore/prev/pathLength) used to chain ORFs into a gene model.
 * Extends PFeature for coordinate/strand handling.
 *
 * @author Brian Bushnell
 */
public class Orf extends PFeature {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Constructs an ORF (Open Reading Frame) with genomic coordinates and sequence data.
	 * Extracts start and stop codons from the provided bases and optionally flips
	 * coordinates for minus-strand features.
	 *
	 * @param scafName_ Scaffold/contig name containing this ORF
	 * @param start_ 0-based start position (first base of start codon)
	 * @param stop_ 0-based stop position (last base of stop codon)
	 * @param strand_ Strand orientation (plus or minus)
	 * @param frame_ Reading frame (0, 1, or 2)
	 * @param bases Nucleotide sequence containing the ORF
	 * @param flip Whether to flip coordinates if on minus strand
	 * @param type_ Feature type (CDS, tRNA, rRNA, etc.)
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
	 * Extracts a 3-base codon from the sequence starting at the specified position.
	 * Encodes the codon as a packed integer using base-to-number conversion.
	 *
	 * @param from Starting position in the sequence
	 * @param bases Nucleotide sequence array
	 * @return Packed integer representation of the codon
	 */
	private static int getCodon(int from, byte[] bases){
		int codon=0;
		for(int i=0; i<3; i++){
//			assert(i+from<bases.length) : i+", "+from+", "+bases.length;
			byte b=bases[i+from];//bounds guard disabled; safe for valid coords (ctor passes from=start or stop-2, both in-bounds)
			int x=AminoAcid.baseToNumber[b];
			codon=(codon<<2)|x;
		}
		return codon;
	}

	public float calcOrfScore(){
		return calcOrfScore(0);
	}

	/** Combines start/stop/inner-kmer sub-scores and ORF length into the heuristic orfScore.
	 * @param overlap Bases overlapping the previous ORF (subtracted as a length penalty).
	 * @return The combined ORF score. */
	public float calcOrfScore(int overlap){
		//Heuristic: geometric-ish blend of start-codon (a), stop-codon (b) and inner-kmer (c) sub-scores,
		//each floored by f1/f2/f3, then scaled by length^2.5 minus an overlap penalty.
		double a=Math.sqrt(Tools.max(f1, e1+startScore));
//		double b=Math.sqrt(f2/*Tools.max(f2, e2+stopScore)*/);//This is better, ignoring stopscore completely
		double b=Math.sqrt(Tools.max(f2, e2+0.35f*stopScore));
		double c=Tools.max(f3, e3+averageKmerScore());
		//TODO: Possible bug [prok/Orf#002] - assert(x!=Double.NaN) is ALWAYS true (NaN != anything, even NaN),
		//so these four asserts never catch a NaN. Correct form: assert(!Double.isNaN(x)). As written, a non-finite
		//c from averageKmerScore() (Orf#001) propagates silently into orfScore. Fix couples to Orf#001 and changes
		//silent->crash -> flag to Brian before applying.
		assert(a!=Double.NaN);
		assert(b!=Double.NaN);
		assert(c!=Double.NaN);
		c=4*Math.pow(c, 2.2);
		double d=(0.1*a*b*c*(Math.pow(length()-overlap, 2.5)-(overlap<1 ? 0 : Math.pow(overlap+50, 2))));//TODO: Adjust these constants
		if(d>0){d=Math.sqrt(d);}
		assert(d!=Double.NaN);
		return (float)d;
	}
	
	/** Mean inner-kmer score: total kmerScore divided by the number of inner kmers spanned by the CDS.
	 * @return kmerScore normalized by (length() - kInnerCDS - 2). */
	public float averageKmerScore(){
		//TODO: Possible bug [prok/Orf#001] - divisor (length()-kInnerCDS-2) is <=0 for ORFs of length<=kInnerCDS+2
		//(kInnerCDS defaults to 6 -> length<=8): ==0 yields +/-Infinity, <0 yields a wrong-signed score. Feeds
		//calcOrfScore (caught by Orf#002's dead NaN guard) and the GFF innerScr field. Reachability pends the
		//GeneCaller min-ORF-length filter. NOTE: active e1..f3 tuning block is labeled "kinnercds=7" but the
		//GeneModel.kInnerCDS default is 6 (GeneModel.java:1018) - a label/constant mismatch to confirm at GeneModel.
		return kmerScore/(length()-GeneModel.kInnerCDS-2); //This slightly affects score if kInnerCDS is changed
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Whether prev is a valid earlier ORF that this one may chain after, during path-finding.
	 * @param prev Candidate predecessor ORF. @param maxOverlap Max bases prev may overlap into this ORF.
	 * @return True iff prev legally precedes this ORF. */
	public boolean isValidPrev(Orf prev, int maxOverlap){
		//valid only if prev strictly precedes: ends before this ends, overlaps by <=maxOverlap, and starts earlier
		if(prev.stop>=stop || prev.stop>=start+maxOverlap || prev.start>=start){return false;}
		//and reject a same-frame, same-strand predecessor whose stop reaches into this ORF's start
		if(prev.frame==frame && prev.strand==strand && prev.stop>=start){return false;}
		return true;
	}

	//path-finding accessors below: the prevStrand arg uses 0=plus, anything else=minus (matches Shared.PLUS=0)
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
	 * Converts ORFs to GFF format lines grouped by feature type.
	 * Creates separate lists for each requested type (CDS, rRNA, tRNA, etc.).
	 *
	 * @param orfs List of ORFs to convert
	 * @param types Comma-separated string of feature types to include
	 * @return Array of GFF line lists, one per requested type
	 */
	public static ArrayList<GffLine>[] toGffLinesByType(ArrayList<Orf> orfs, String types){
		String[] typeArray=types.split(",");
		@SuppressWarnings("unchecked")
		ArrayList<GffLine>[] lists=new ArrayList[typeArray.length];
		HashMap<String, ArrayList<GffLine>> map=new HashMap<String, ArrayList<GffLine>>(3*typeArray.length);
		for(int i=0; i<typeArray.length; i++){
			String type=typeArray[i];
			lists[i]=new ArrayList<GffLine>();
			//TODO: Possible bug [prok/Orf#003] - a duplicate token in `types` makes this map.put overwrite the
			//earlier mapping, so matching ORFs route only to the LAST list for that type and the earlier
			//duplicate's list is returned empty (silent data loss, no warning). Latent LOW: callers pass fixed,
			//non-duplicate type sets.
			map.put(type,  lists[i]);
		}
		for(Orf orf : orfs){
			String type=ProkObject.typeStrings2[orf.type];
			ArrayList<GffLine> glist=map.get(type);
			if(glist!=null) {//null when this ORF's type wasn't among the requested `types` -> silently skip it
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
