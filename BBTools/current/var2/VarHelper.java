package var2;

import java.util.ArrayList;

import align2.QualityTools;
import dna.AminoAcid;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.SamLine;

/**
 * Utility class providing helper methods for variant processing and output formatting.
 * Contains static methods for generating file headers, calculating scores, analyzing homopolymers,
 * and processing junction variants. Extracted from Var class for better code organization.
 * 
 * @author Brian Bushnell
 * @author Isla Winglet  
 * @date June 2025
 */
public class VarHelper {

	/*--------------------------------------------------------------*/
	/*----------------      Header Generation       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Generates comprehensive header for VAR format output files.
	 * Includes metadata about sequencing run, processing parameters, and column definitions.
	 * 
	 * @param properPairRate Fraction of reads that mapped as proper pairs
	 * @param totalQualityAvg Average base quality across all processed bases
	 * @param mapqAvg Average mapping quality across all reads
	 * @param rarity Minimum variant frequency threshold
	 * @param minAlleleFraction Minimum allele fraction for variant calling
	 * @param ploidy Expected ploidy level for organism
	 * @param reads Total number of reads processed
	 * @param pairs Number of paired reads in sequencing
	 * @param properPairs Number of properly paired reads
	 * @param bases Total number of bases processed
	 * @param ref Reference genome file path
	 * @return Complete VAR format header string
	 */
	public static String toVarHeader(double properPairRate, double totalQualityAvg, double mapqAvg, double rarity, double minAlleleFraction, int ploidy, 
			long reads, long pairs, long properPairs, long bases, String ref){
		StringBuilder sb=new StringBuilder();
		
		final double readLengthAvg=bases/Tools.max(1.0, reads);
		sb.append("#fileformat\tVar_"+Var.varFormat+"\n");
		sb.append("#BBMapVersion\t"+Shared.BBTOOLS_VERSION_STRING+"\n");
		sb.append("#ploidy\t"+ploidy+"\n");
		sb.append(Tools.format("#rarity\t%.5f\n", rarity));
		sb.append(Tools.format("#minAlleleFraction\t%.4f\n", minAlleleFraction));
		sb.append("#reads\t"+reads+"\n");
		sb.append("#pairedReads\t"+pairs+"\n");
		sb.append("#properlyPairedReads\t"+properPairs+"\n");
		sb.append(Tools.format("#readLengthAvg\t%.2f\n", readLengthAvg));
		sb.append(Tools.format("#properPairRate\t%.4f\n", properPairRate));
		sb.append(Tools.format("#totalQualityAvg\t%.4f\n", totalQualityAvg));
		sb.append(Tools.format("#mapqAvg\t%.2f\n", mapqAvg));
		if(ref!=null){sb.append("#reference\t"+ref+"\n");}
		
		// Column headers for variant data
		sb.append("#scaf\tstart\tstop\ttype\tcall\tr1p\tr1m\tr2p\tr2m\tpaired\tlengthSum");
		sb.append("\tmapq\tmapqMax\tbaseq\tbaseqMax\tedist\tedistMax\tid\tidMax");
		sb.append("\tcov\tminusCov");
		sb.append("\tnearbyVarCount");
		sb.append("\tflagged");
		sb.append("\tcontigEndDist");
		sb.append("\tphredScore");
		
		// Extended columns if detailed output is enabled
		if(Var.extendedText){
			sb.append("\treadCount\talleleFraction\trevisedAF\tstrandRatio\tbaseqAvg\tmapqAvg\tedistAvg\tidentityAvg");
			sb.append("\tedistScore\tidentityScore\tqualityScore\tpairedScore\tbiasScore\tcoverageScore\thomopolymerScore\tscore");
		}
		return sb.toString();
	}

	/**
	 * Generates basic column header for simplified VAR output.
	 * @return Basic VAR format header without metadata
	 */
	public static String toBasicHeader(){
		StringBuilder sb=new StringBuilder();
		sb.append("#scaf\tstart\tstop\ttype\tcall\tr1p\tr1m\tr2p\tr2m\tpaired\tlengthSum");
		sb.append("\tmapq\tmapqMax\tbaseq\tbaseqMax\tedist\tedistMax\tid\tidMax");
		sb.append("\tcov\tminusCov");
		sb.append("\tnearVarCount");
		sb.append("\tcontigEndDist");
		sb.append("\tphredScore");
		return sb.toString();
	}

	/**
	 * Generates comprehensive VCF format header with metadata and field definitions.
	 * Includes all INFO and FORMAT field descriptions required by VCF specification.
	 * 
	 * @param properPairRate Fraction of properly paired reads
	 * @param totalQualityAvg Average base quality score
	 * @param mapqAvg Average mapping quality score
	 * @param rarity Minimum variant frequency threshold
	 * @param minAlleleFraction Minimum allele fraction for calling
	 * @param ploidy Expected organism ploidy
	 * @param reads Total reads processed
	 * @param pairs Total paired reads
	 * @param properPairs Total properly paired reads
	 * @param bases Total bases processed
	 * @param ref Reference genome file path
	 * @param map Scaffold mapping information
	 * @param sampleName Sample identifier for output
	 * @param trimWhitespace Whether to trim scaffold names
	 * @return Complete VCF format header
	 */
	public static String toVcfHeader(double properPairRate, double totalQualityAvg, double mapqAvg,
			double rarity, double minAlleleFraction, int ploidy, long reads, long pairs,
			long properPairs, long bases, String ref, ScafMap map, String sampleName,
			boolean trimWhitespace) {
		StringBuilder sb=new StringBuilder();
		
		// VCF format and version information
		sb.append("##fileformat=VCFv4.2\n");
		sb.append("##BBMapVersion="+Shared.BBTOOLS_VERSION_STRING+"\n");
		sb.append("##ploidy="+ploidy+"\n");
		sb.append(Tools.format("##rarity=%.5f\n", rarity));
		sb.append(Tools.format("##minallelefraction=%.5f\n", minAlleleFraction));
		sb.append("##reads="+reads+"\n");
		sb.append("##pairedReads="+pairs+"\n");
		sb.append("##properlyPairedReads="+properPairs+"\n");
		sb.append(Tools.format("##readLengthAvg=%.3f\n", (bases/Tools.max(reads, 1.0))));
		sb.append(Tools.format("##properPairRate=%.5f\n", properPairRate));
		sb.append(Tools.format("##totalQualityAvg=%.3f\n", totalQualityAvg));
		sb.append(Tools.format("##mapqAvg=%.3f\n", mapqAvg));
		if(ref!=null){sb.append("##reference="+ref+"\n");}
		
		// Contig definitions for each scaffold
		for(Scaffold scaf : map.list){
			String name=scaf.name;
			if(trimWhitespace){name=Tools.trimWhitespace(name);}
			sb.append("##contig=<ID="+name+",length="+scaf.length+">\n");
		}
		
		// Filter and field definitions
		{
			sb.append("##FILTER=<ID=FAIL,Description=\"Fail\">\n");
			sb.append("##FILTER=<ID=PASS,Description=\"Pass\">\n");
			
			// INFO field definitions
			sb.append("##INFO=<ID=SN,Number=1,Type=Integer,Description=\"Scaffold Number\">\n");
			sb.append("##INFO=<ID=STA,Number=1,Type=Integer,Description=\"Start\">\n");
			sb.append("##INFO=<ID=STO,Number=1,Type=Integer,Description=\"Stop\">\n");
			sb.append("##INFO=<ID=TYP,Number=1,Type=String,Description=\"Type\">\n");
			
			sb.append("##INFO=<ID=R1P,Number=1,Type=Integer,Description=\"Read1 Plus Count\">\n");
			sb.append("##INFO=<ID=R1M,Number=1,Type=Integer,Description=\"Read1 Minus Count\">\n");
			sb.append("##INFO=<ID=R2P,Number=1,Type=Integer,Description=\"Read2 Plus Count\">\n");
			sb.append("##INFO=<ID=R2M,Number=1,Type=Integer,Description=\"Read2 Minus Count\">\n");
			
			sb.append("##INFO=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n");
			sb.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
			sb.append("##INFO=<ID=MCOV,Number=1,Type=Integer,Description=\"Minus Coverage\">\n");
			sb.append("##INFO=<ID=PPC,Number=1,Type=Integer,Description=\"Paired Count\">\n");

			sb.append("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Fraction\">\n");
			sb.append("##INFO=<ID=RAF,Number=1,Type=Float,Description=\"Revised Allele Fraction\">\n");
			sb.append("##INFO=<ID=LS,Number=1,Type=Integer,Description=\"Length Sum\">\n");

			sb.append("##INFO=<ID=MQS,Number=1,Type=Integer,Description=\"MAPQ Sum\">\n");
			sb.append("##INFO=<ID=MQM,Number=1,Type=Integer,Description=\"MAPQ Max\">\n");
			sb.append("##INFO=<ID=BQS,Number=1,Type=Integer,Description=\"Base Quality Sum\">\n");
			sb.append("##INFO=<ID=BQM,Number=1,Type=Integer,Description=\"Base Quality Max\">\n");
			
			sb.append("##INFO=<ID=EDS,Number=1,Type=Integer,Description=\"End Distance Sum\">\n");
			sb.append("##INFO=<ID=EDM,Number=1,Type=Integer,Description=\"End Distance Max\">\n");
			sb.append("##INFO=<ID=IDS,Number=1,Type=Integer,Description=\"Identity Sum\">\n");
			sb.append("##INFO=<ID=IDM,Number=1,Type=Integer,Description=\"Identity Max\">\n");

			sb.append("##INFO=<ID=NVC,Number=1,Type=Integer,Description=\"Nearby Variation Count\">\n");
			sb.append("##INFO=<ID=FLG,Number=1,Type=Integer,Description=\"Flagged\">\n");
			sb.append("##INFO=<ID=CED,Number=1,Type=Integer,Description=\"Contig End Distance\">\n");
			sb.append("##INFO=<ID=HMP,Number=1,Type=Integer,Description=\"Homopolymer Count\">\n");
			sb.append("##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand Bias\">\n");
			sb.append("##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Ref+, Ref-, Alt+, Alt-\">\n");
			
			// FORMAT field definitions
			sb.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
			sb.append("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
			sb.append("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n");
			sb.append("##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele Fraction\">\n");
			sb.append("##FORMAT=<ID=RAF,Number=1,Type=Float,Description=\"Revised Allele Fraction\">\n");
			sb.append("##FORMAT=<ID=NVC,Number=1,Type=Integer,Description=\"Nearby Variation Count\">\n");
			sb.append("##FORMAT=<ID=FLG,Number=1,Type=Integer,Description=\"Flagged\">\n");
			sb.append("##FORMAT=<ID=SB,Number=1,Type=Float,Description=\"Strand Bias\">\n");
			
			sb.append("##FORMAT=<ID=SC,Number=1,Type=Float,Description=\"Score\">\n");
			sb.append("##FORMAT=<ID=PF,Number=1,Type=String,Description=\"Pass Filter\">\n");
		}
		
		// Column headers
		sb.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
		if(sampleName!=null){sb.append('\t').append(sampleName);}
		return sb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------       Score Calculation      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Converts variant quality score to Phred scale.
	 * Applies scaling factors to convert internal scoring to standard Phred format.
	 * 
	 * @param score Internal variant quality score (0.0 to 1.0)
	 * @return Phred-scaled quality score
	 */
	public static double toPhredScore(double score){
		if(score==0){return 0;}
		score=score*0.998; // Apply scaling factor
		return 2.5*QualityTools.probErrorToPhredDouble(1-score);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Homopolymer Analysis      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Counts homopolymer length around a substitution position.
	 * Examines up to 4 bases in each direction for runs of identical bases.
	 * 
	 * @param bases Reference sequence
	 * @param pos Position of substitution
	 * @param base Base being substituted
	 * @return Length of homopolymer run (0-8)
	 */
	public static int homopolymerCountSub(final byte[] bases, final int pos, final byte base){
		if(pos<0 || pos>=bases.length){return 0;}
		if(!AminoAcid.isFullyDefined(base)){return 0;}
		
		// Count matching bases to the left
		int count1=0;
		for(int i=pos-1, lim=Tools.max(0, pos-4); i>=lim; i--){
			if(bases[i]==base){count1++;}
			else{break;}
		}
		
		// Count matching bases to the right
		int count2=0;
		for(int i=pos+1, lim=Tools.min(bases.length, pos+5); i<lim; i++){
			if(bases[i]==base){count2++;}
			else{break;}
		}
		assert(count1+count2<=8) : count1+", "+count2;
		
		// Add 1 if there are matches on both sides (indicates continuous run)
		return count1+count2+(count1>0 && count2>0 ? 1 : 0);
	}

	/**
	 * Counts homopolymer length extending leftward from position.
	 * @param bases Sequence data
	 * @param pos Starting position
	 * @param base Base to count
	 * @return Count of consecutive matching bases leftward
	 */
	public static int homopolymerCountLeft(final byte[] bases, final int pos, final byte base){
		if(pos<0 || bases[pos]!=base){return 0;}
		if(!AminoAcid.isFullyDefined(base)){return 0;}
		
		int count=0;
		for(int i=pos, lim=Tools.max(0, pos-3); i>=lim; i--){
			if(bases[i]==base){count++;}
			else{break;}
		}
		return count;
	}

	/**
	 * Counts homopolymer length extending rightward from position.
	 * @param bases Sequence data
	 * @param pos Starting position
	 * @param base Base to count
	 * @return Count of consecutive matching bases rightward
	 */
	public static int homopolymerCountRight(final byte[] bases, final int pos, final byte base){
		if(pos<0 || bases[pos]!=base){return 0;}
		if(!AminoAcid.isFullyDefined(base)){return 0;}
		
		int count=0;
		for(int i=pos, lim=Tools.min(bases.length, pos+4); i<lim; i++){
			if(bases[i]==base){count++;}
			else{break;}
		}
		return count;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Clipping Analysis       ----------------*/
	/*--------------------------------------------------------------*/

	/** 
	 * Counts left-clipped bases from match string (handles both short and long format).
	 * @param match Match string from read alignment
	 * @return Number of left-clipped bases
	 */
	private static int countLeftClip(byte[] match){
		byte mode=match[0];
		if(mode!='C'){return 0;} // No clipping if first character isn't 'C'
		
		int current=0;
		boolean hasDigit=false;
		
		// Parse clipping specification (e.g., "C12" or "CCCCCCCCCCCC")
		for(int mpos=0; mpos<match.length; mpos++){
			byte c=match[mpos];
			
			if(mode==c){
				current++; // Count 'C' characters
			}else if(Tools.isDigit(c)){
				current=(hasDigit ? current : 0)*10+c-'0'; // Parse number
				hasDigit=true;
			}else{
				break; // End of clipping specification
			}
		}
		return current;
	}

	/** 
	 * Counts right-clipped bases from match string (handles both short and long format).
	 * @param match Match string from read alignment
	 * @return Number of right-clipped bases
	 */
	private static int countRightClip(byte[] match){
		int mpos=match.length-1;
		boolean hasDigit=false;
		
		// Look for digits at end (short format)
		while(mpos>=0 && Tools.isDigit(match[mpos])){
			hasDigit=true;
			mpos--;
		}
		if(!hasDigit){return countRightClipLong(match);} // Use long format parser
		
		byte mode=match[mpos];
		if(mode!='C'){return 0;}
		int current=0;
		
		// Parse right clipping specification
		for(; mpos<match.length; mpos++){
			byte c=match[mpos];
			
			if(mode==c){
				current++;
			}else if(Tools.isDigit(c)){
				current=(hasDigit ? current : 0)*10+c-'0';
				hasDigit=true;
			}else{
				break;
			}
		}
		return current;
	}
	
	/** 
	 * Counts left clipping in long match format (e.g., "CCCCCCMMMM").
	 * @param longmatch Long format match string
	 * @return Number of left-clipped bases
	 */
	private static int countLeftClipLong(byte[] longmatch){
		for(int i=0; i<longmatch.length; i++){
			if(longmatch[i]!='C'){return i;}
		}
		return longmatch.length;
	}
	
	/** 
	 * Counts right clipping in long match format.
	 * @param longmatch Long format match string
	 * @return Number of right-clipped bases
	 */
	private static int countRightClipLong(byte[] longmatch){
		for(int i=0, j=longmatch.length-1; i<longmatch.length; i++, j--){
			if(longmatch[j]!='C'){return i;}
		}
		return longmatch.length;
	}

	/*--------------------------------------------------------------*/
	/*----------------     Junction Processing       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Identifies junction variants from clipped alignments.
	 * Detects left and right junction points where reads are clipped, indicating
	 * potential structural variations or assembly breaks.
	 * 
	 * @param r Read with alignment information
	 * @param sl SAM line with mapping details
	 * @param scafnum Scaffold number for variant coordinates
	 * @param containsVars Whether read contains other variants
	 * @param minClip Minimum clipping length to consider junction
	 * @return List of junction variants, or null if none found
	 */
	static ArrayList<Var> toJunctions(Read r, SamLine sl, final int scafnum,
			boolean containsVars, final int minClip){
		final byte[] match0=r.match, match;
		final byte[] bases=r.bases;
		final int start, stop;
		
		// Count clipping in original alignment
		int leftClip=countLeftClip(match0), rightClip=countRightClip(match0);
		
		if(leftClip==0 && rightClip==0){
			// Try soft-clipping if no hard clipping found
			int[] rvec=KillSwitch.allocInt1D(2);
			byte[] smatch=SoftClipper.softClipMatch(match0, minClip, false, r.start, r.stop, rvec);
			if(smatch==null){
				return null;
			}else{
				start=rvec[0];
				stop=rvec[1];
				match=smatch;
				leftClip=countLeftClip(match);
				rightClip=countRightClip(match);
			}
		}else{
			if(leftClip<minClip && rightClip<minClip){return null;}
			start=r.start;
			stop=r.stop;
			match=match0;
		}
		
		ArrayList<Var> list=new ArrayList<Var>();
		
		// Create left junction variant if sufficient clipping
		if(leftClip>=minClip){
			int bpos=leftClip-1; // Base position in read
			byte jcall=bases[bpos]; // Junction call base
			int jpos=start+leftClip; // Genomic position
			Var v=new Var(scafnum, jpos, jpos+1, jcall, Var.LJUNCT);
			v.add(r, bpos, bpos+1);
			list.add(v);
		}
		
		// Create right junction variant if sufficient clipping
		if(rightClip>=minClip){
			int bpos=bases.length-rightClip; // Base position in read
			byte jcall=bases[bpos]; // Junction call base
			int jpos=stop-rightClip+1; // Genomic position
			Var v=new Var(scafnum, jpos, jpos+1, jcall, Var.RJUNCT);
			v.add(r, bpos, bpos+1);
			list.add(v);
		}
		
		return list;
	}
}