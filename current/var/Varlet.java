package var;

import java.util.ArrayList;

import align2.QualityTools;
import dna.Gene;
import fileIO.TextFile;
import shared.Shared;
import shared.Tools;

/**
 * Detailed representation of a genetic variation with comprehensive read and quality metadata.
 * Stores and processes variant information from sequencing reads, capturing strand, quality,
 * and mapping details. Packs quality metrics into single integer using bitwise operations.
 * @author Brian Bushnell
 */
public class Varlet extends var.Variation {
	
	private static final long serialVersionUID = -606340580378991068L;

	public Varlet(int chrom_, byte strand_, int start_, int stop_, int matchStart_, int matchStop_, byte vType, String rf, String ca,
			 int varQuality_, int readQuality_, int mapScore_, int errors_, float expectedErrors_, int paired_, long readID_,
			 int readLen_,
			 int readStart_, int readStop_, int readCopies_, int headDist_, int tailDist_, int endDist_, int pairnum){
		super(chrom_, start_, stop_, vType, rf, ca);
		strand=strand_;
		
		setQvector(varQuality_, readQuality_, varQuality_, readQuality_);
		
		mapScore=mapScore_;
		errors=errors_;
		expectedErrors=expectedErrors_;
		paired=paired_;

		matchStart=matchStart_;
		matchStop=matchStop_;
		
		readID=readID_;
		readLen=readLen_;
		
		readStart=readStart_;
		readStop=readStop_;
		
		numReads=Tools.min(readCopies_, Short.MAX_VALUE);
		
		
		headDist=headDist_;
		tailDist=tailDist_;
		endDist=endDist_;
		
		if(pairnum==0){
			if(strand==Shared.PLUS){numPlusReads1=1;}
			else{numMinusReads1=1;}
		}else{
			if(strand==Shared.PLUS){numPlusReads2=1;}
			else{numMinusReads2=1;}
		}
		
		assert(pairnum==0 || pairnum==1) : pairnum+"\n"+this;
//		assert(readID_<Integer.MAX_VALUE) : readID_+"\n"+this;
		assert(readCopies_>=1) : readCopies_+"\n"+this;
//		assert(readCopies_<Short.MAX_VALUE) : readCopies_+"\n"+this;
		
		assert(endDist<=tailDist) : this;
		assert(endDist<=headDist) : this;
		
		assert(readStart<readStop) : this;
	}
	
	@Override
	public String toString(){return toText().toString();}
	
	/** Returns the header string for text output format.
	 * @return Tab-delimited header string for Varlet text representation */
	public static String header(){return textHeader().toString();}
	
	/**
	 * Creates the tab-delimited header for text output format.
	 * Contains column names for all Varlet fields in output order.
	 * @return StringBuilder containing the complete header line with tabs
	 */
	public static CharSequence textHeader(){
		StringBuilder sb=new StringBuilder(64);
		
		sb.append("#");
		sb.append("chrom").append('\t');
		sb.append("strand").append('\t');
		sb.append("readStart").append('\t');
		sb.append("readStop").append('\t');
		sb.append("varStart").append('\t');
		sb.append("varStop").append('\t');
		
		sb.append("type").append('\t');
		sb.append("mapScore").append('\t');
		sb.append("errors").append('\t');
		sb.append("expectedErrors").append('\t');
		sb.append("readID").append('\t');
		sb.append("readLen").append('\t');
		sb.append("headDist").append('\t');
		sb.append("tailDist").append('\t');
		sb.append("endDist").append('\t');

		sb.append("avgVarQuality").append('\t');
		sb.append("maxVarQuality").append('\t');
		sb.append("avgReadQuality").append('\t');
		sb.append("maxReadQuality").append('\t');
		
		sb.append("numReads").append('\t');
		sb.append("numSemiUnique").append('\t');
		sb.append("numUniqueReads").append('\t');
		sb.append("paired").append('\t');
		sb.append("plusReads1").append('\t');
		sb.append("minusReads1").append('\t');
		sb.append("plusReads2").append('\t');
		sb.append("minusReads2").append('\t');
		
		sb.append("ref").append('\t');
		sb.append("call");
		return sb;
	}
	
	/**
	 * Converts this Varlet to tab-delimited text representation.
	 * Includes all variant properties, quality metrics, and strand counts.
	 * @return StringBuilder with tab-separated values matching textHeader format
	 */
	public final StringBuilder toText(){
		StringBuilder sb=new StringBuilder(64);
		
		sb.append(chromosome).append('\t');
		sb.append(Shared.strandCodes[strand]).append('\t');
		sb.append(readStart).append('\t');
		sb.append(readStop).append('\t');
		sb.append(beginLoc).append('\t');
		sb.append(endLoc).append('\t');
		sb.append(Variation.varTypeMap[varType]).append('\t');
		
		sb.append(mapScore).append('\t');
		sb.append(errors).append('\t');
		sb.append(Tools.format("%.1f", expectedErrors)).append('\t');
		sb.append(readID).append('\t');
		sb.append(readLen).append('\t');
		sb.append(headDist).append('\t');
		sb.append(tailDist).append('\t');
		sb.append(endDist).append('\t');

		sb.append(avgVarQuality()).append('\t');
		sb.append(maxVarQuality()).append('\t');
		sb.append(avgReadQuality()).append('\t');
		sb.append(maxReadQuality()).append('\t');
		
		sb.append(numReads).append('\t');
		sb.append(numSemiUniqueReads).append('\t');
		sb.append(numUniqueReads).append('\t');
		sb.append(paired).append('\t');
		sb.append(numPlusReads1).append('\t');
		sb.append(numMinusReads1).append('\t');
		sb.append(numPlusReads2).append('\t');
		sb.append(numMinusReads2).append('\t');
		
		sb.append(ref==null || ref.length()==0 ? "." : ref).append('\t');
		sb.append(call==null || call.length()==0 ? "." : call);
		
//		if(coverageAtLoc>0){sb.append("\t"+coverageAtLoc);}
		return sb;
	}
	
	/**
	 * Reads Varlet objects from a text file.
	 * Skips comment lines starting with '#' and parses each data line.
	 * @param fname Path to the input text file
	 * @return ArrayList of parsed Varlet objects
	 */
	public static final ArrayList<Varlet> fromTextFile(String fname){
		TextFile tf=new TextFile(fname, false);
		ArrayList<Varlet> list=new ArrayList<Varlet>(2000);
		
		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(s.charAt(0)!='#'){
				Varlet v=Varlet.fromText(s);
				list.add(v);
			}
		}
		tf.close();
		list.trimToSize();
		return list;
	}
	
	/**
	 * Parses a single tab-delimited line into a Varlet object.
	 * Expects fields in the order defined by textHeader method.
	 * @param line Tab-delimited string containing variant data
	 * @return Parsed Varlet object
	 */
	public static final Varlet fromText(String line){
		String[] split=line.split("\t");
		
		int chrom=Byte.parseByte(split[0]);
		byte strand=Gene.toStrand(split[1]);
		int readStart=Integer.parseInt(split[2]);
		int readStop=Integer.parseInt(split[3]);
		int start=Integer.parseInt(split[4]);
		int stop=Integer.parseInt(split[5]);
		byte varType=Variation.varTypeMap2.get(split[6]);
		
		int mapScore=Integer.parseInt(split[7]);
		int errors=Integer.parseInt(split[8]);
		float expectedErrors=Float.parseFloat(split[9]);
		
		long readID=Integer.parseInt(split[10]);
		int readLen=Integer.parseInt(split[11]);
		int headDist=Integer.parseInt(split[13]);
		int tailDist=Integer.parseInt(split[14]);
		int endDist=Integer.parseInt(split[15]);

		int avgVarQuality=Integer.parseInt(split[16]);
		int maxVarQuality=Integer.parseInt(split[17]);
		int avgReadQuality=Integer.parseInt(split[18]);
		int maxReadQuality=Integer.parseInt(split[19]);
		int numReads=Integer.parseInt(split[20]);
		int numSemiUniqueReads=Integer.parseInt(split[21]);
		int numUniqueReads=Integer.parseInt(split[22]);
		int paired=Integer.parseInt(split[23]);
		int numPlusReads1=Integer.parseInt(split[24]);
		int numMinusReads1=Integer.parseInt(split[25]);
		int numPlusReads2=Integer.parseInt(split[26]);
		int numMinusReads2=Integer.parseInt(split[27]);
		
		String ref=split[28];
		String call=split[29];
		if(ref.length()==1 && ref.charAt(0)=='.'){ref=null;}
		if(call.length()==1 && call.charAt(0)=='.'){call=null;}
		
		
		
		Varlet v=new Varlet(chrom, strand, start, stop, -1, -1, varType, ref, call, avgVarQuality, avgReadQuality,
				mapScore, errors, expectedErrors, paired, readID, readLen, readStart, readStop, numReads,
				headDist, tailDist, endDist, 1);
		
		v.setQvector(avgVarQuality, avgReadQuality, maxVarQuality, maxReadQuality);
		v.numPlusReads1=numPlusReads1;
		v.numMinusReads1=numMinusReads1;
		v.numPlusReads2=numPlusReads2;
		v.numMinusReads2=numMinusReads2;
		v.numSemiUniqueReads=numSemiUniqueReads;
		v.numUniqueReads=numUniqueReads;
		
		return v;
	}
	

	@Override
	public boolean equals(Variation other){
//		assert(other.getClass()!=Varlet.class);
		return super.compareTo(other)==0;
	}
	
	//DO NOT enable this!  Varlets should use equality based on Variation data only.
//	public boolean equals(Varlet other){
//		return compareTo(other)==0;
//	}

	@Override
	public int compareTo(Variation other) {
//		if(other.getClass()==Varlet.class){} //not needed in practice
		return(compareTo((Varlet)other));
	}
	
	/**
	 * Compares this Varlet to another for sorting.
	 * Primary sort by chromosome, position, variant type, then call sequence.
	 * Secondary sort by read coordinates, strand, and quality.
	 *
	 * @param other The Varlet to compare against
	 * @return Negative, zero, or positive integer for less than, equal, or greater
	 */
	public int compareTo(Varlet other) {
		
//		int a=compareTo2(other);
//		int b=other.compareTo2(this);
//		assert(a==-b) : "\n"+a+", "+b+"\n"+Varlet.header()+"\n"+this+"\n"+other+"\n";
	
		if(chromosome!=other.chromosome){return chromosome-other.chromosome;}
		if(beginLoc!=other.beginLoc){return other.beginLoc>beginLoc ? -1 : 1;}
		if(endLoc!=other.endLoc){return other.endLoc>endLoc ? -1 : 1;}
		if(varType!=other.varType){return varType-other.varType;}
		if(varType==REF || varType==NOCALL){return 0;}

		if(call==null && other.call!=null){return -1;}
		if(call!=null && other.call==null){return 1;}
		if(call!=null && other.call!=null){
			int x=call.compareTo(other.call);
			if(x!=0){return x;}
		}

		if(readStart!=other.readStart){return readStart-other.readStart;}
		if(readStop!=other.readStop){return readStop-other.readStop;}
		if(strand!=other.strand){return strand-other.strand;}
		if(maxVarQuality()!=other.maxVarQuality()){return other.maxVarQuality()<maxVarQuality() ? -1 : 1;}
		
		return 0;
	}
	
	/** TODO: Add expected errors, tailDist, endDist */
	public int score(){
		int score=1000/(errors+1);
		score+=(int)(500/(expectedErrors+1));
		score+=Tools.max(0, (1000-(int)(16000*QualityTools.PROB_ERROR[maxReadQuality()])));
		score+=Tools.max(0, (1000-(int)(16000*QualityTools.PROB_ERROR[maxVarQuality()])));
		score+=10*Tools.min(35, maxVarQuality());
		score+=Tools.max(0, (200-(int)(8000*QualityTools.PROB_ERROR[avgVarQuality()])));
		score+=Tools.max(0, (200-(int)(8000*QualityTools.PROB_ERROR[avgReadQuality()])));
		score+=(1000-2000/(paired+2));
		score+=(500-1000/(numSemiUniqueReads+2));
		score+=(500-1000/(numUniqueReads+2));
		score+=(200-400/(numReads+2));
		score+=(50*Tools.min(20, tailDist));
		score+=(50*Tools.min(10, endDist));
		
		int lenFactor=Tools.min(readLen, 100);
		score+=(1000*lenFactor)/(lenFactor+100);
		
		score+=Tools.min(1000, (10*mapScore)/readLen); //TODO: This is temporary, until Read correctly supports mapLen in toText()
		score+=(1000-1000/(1+minStrandReads()));
		return score;
	}
	
	
	/**
	 * Packed quality vector containing avgVarQuality, avgReadQuality, maxVarQuality, maxReadQuality
	 */
	private int qvector;
	
	/** Gets average variant quality from packed qvector */
	public int avgVarQuality(){return qvector&0xFF;}
	/** Gets average read quality from packed qvector */
	public int avgReadQuality(){return (qvector>>8)&0xFF;};
	/** Gets maximum variant quality from packed qvector */
	public int maxVarQuality(){return (qvector>>16)&0xFF;};
	/** Gets maximum read quality from packed qvector */
	public int maxReadQuality(){return (qvector>>24)&0xFF;};
	
	/** Sets average variant quality in packed qvector.
	 * @param value Average variant quality score (0-255) */
	public void setAvgVarQuality(int value){
		qvector=((qvector&0xFFFFFF00)|(value&0xFF));
	}
	/** Sets average read quality in packed qvector.
	 * @param value Average read quality score (0-255) */
	public void setAvgReadQuality(int value){
		qvector=((qvector&0xFFFF00FF)|((value&0xFF)<<8));
	}
	/** Sets maximum variant quality in packed qvector.
	 * @param value Maximum variant quality score (0-255) */
	public void setMaxVarQuality(int value){
		qvector=((qvector&0xFF00FFFF)|((value&0xFF)<<16));
	}
	/** Sets maximum read quality in packed qvector.
	 * @param value Maximum read quality score (0-255) */
	public void setMaxReadQuality(int value){
		qvector=((qvector&0x00FFFFFF)|((value&0xFF)<<24));
	}
	/**
	 * Sets all quality values in packed qvector using bitwise operations.
	 *
	 * @param avq Average variant quality
	 * @param arq Average read quality
	 * @param mvq Maximum variant quality
	 * @param mrq Maximum read quality
	 */
	public void setQvector(int avq, int arq, int mvq, int mrq){
		qvector=mrq&0xFF;
		qvector=(qvector<<8)|(mvq&0xFF);
		qvector=(qvector<<8)|(arq&0xFF);
		qvector=(qvector<<8)|(avq&0xFF);
	}
	
	
	
	/** Mapping score of the read containing this variant */
	public int mapScore;
	/** Number of errors in the read alignment */
	public int errors;
	
	/** Expected number of errors based on quality scores */
	public float expectedErrors;

	/** Start position of the match within the alignment */
	public int matchStart;
	/** Stop position of the match within the alignment */
	public int matchStop;
	
	/** Start position of the read alignment on the reference */
	public int readStart;
	/** Stop position of the read alignment on the reference */
	public int readStop;

	/** Distance from variant to 5' end of read */
	public int headDist;
	/** Distance from variant to 3' end of read */
	public int tailDist;
	/** Distance from variant to nearest read end */
	public int endDist;
	
	/** Strand orientation of the read (+ or -) */
	public byte strand;
	/** Number of paired reads supporting this variant */
	public int paired;

	/** Unique identifier of the read containing this variant */
	public long readID;
	
	/** Length of read when used for calling vars; ie, after being trimmed. */
	public int readLen;
	
	/** Total number of reads supporting this variant */
	public int numReads;
	/** Number of semi-unique reads supporting this variant */
	public int numSemiUniqueReads=1;
	/** Number of unique reads supporting this variant */
	public int numUniqueReads=1;
	
	/** Varlets from read 1 mapped to plus strand */
	public int numPlusReads1=0;
	
	/** Varlets from read 1 mapped to minus strand */
	public int numMinusReads1=0;
	
	/** Varlets from read 2 mapped to plus strand */
	public int numPlusReads2=0;
	
	/** Varlets from read 2 mapped to minus strand */
	public int numMinusReads2=0;
	
	/** Number of reads1 and reads2 mapped to the plus strand */
	public int numPlusMappedReads(){
		return numPlusReads1+numPlusReads2;
	}
	
	/** Number of reads1 and reads2 from which the original molecule (i.e., read 1) mapped to the plus strand */
	public int numPlusOriginReads(){
		return numPlusReads1+numMinusReads2;
	}
	
	/** Number of reads1 and reads2 mapped to the minus strand */
	public int numMinusMappedReads(){
		return numMinusReads1+numMinusReads2;
	}
	
	/** Number of reads1 and reads2 from which the original molecule (i.e., read 1) mapped to the minus strand */
	public int numMinusOriginReads(){
		return numMinusReads1+numPlusReads2;
	}
	
	/** Gets minimum number of reads between plus and minus mapped strands.
	 * @return Minimum of numPlusMappedReads and numMinusMappedReads */
	public int minStrandReads(){return Tools.min(numPlusMappedReads(), numMinusMappedReads());}
	
//	public byte numStrands(){return (byte)((numPlusReads>0 ? 1 : 0)+(numMinusReads>0 ? 1 : 0));}
	/** Gets minimum read count among all four strand-pair combinations.
	 * @return Minimum of numPlusReads1, numMinusReads1, numPlusReads2, numMinusReads2 */
	public int minStrandReads4(){return Tools.min(numPlusReads1, numMinusReads1, numPlusReads2, numMinusReads2);}
	/**
	 * Gets second lowest number among four strand-pair combinations.
	 * Complex algorithm to find second minimum value.
	 * @return Second lowest read count among the four strand-pair values
	 */
	public int minStrandReads3(){//return second lowest number

		final int a, b, c, d;
		if(numPlusReads1<=numMinusReads1){a=numPlusReads1; b=numMinusReads1;}
		else{b=numPlusReads1; a=numMinusReads1;}
		if(numPlusReads2<=numMinusReads2){c=numPlusReads2; d=numMinusReads2;}
		else{d=numPlusReads2; c=numMinusReads2;}
		
		return Tools.min(b, d, (a>=c ? a : c));
		
	}
	/**
	 * Counts how many of the four strand-pair combinations have reads.
	 * @return Number of non-zero values among numPlusReads1, numMinusReads1,
	 * numPlusReads2, numMinusReads2
	 */
	public int strandReadCount(){
		return (numPlusReads1>0 ? 1 : 0)+(numMinusReads1>0 ? 1 : 0)+(numPlusReads2>0 ? 1 : 0)+(numMinusReads2>0 ? 1 : 0);
	}
	
	/** Determines which read pair this variant came from.
	 * @return 0 if from read1 pair, 1 if from read2 pair */
	public int pairNum(){
		return (numPlusReads1+numMinusReads1)>0 ? 0 : 1;
	}
	
}
