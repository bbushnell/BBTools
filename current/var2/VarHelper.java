package var2;

import java.util.ArrayList;

import align2.QualityTools;
import dna.AminoAcid;
import ml.CellNet;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.SamLine;

/**
 * Shared variant helpers for header text, score conversion, bounded homopolymer
 * statistics and terminal-clipping junction calls. Extracted from Var.
 * Header methods return text without writing files; the junction helper creates
 * new Var observations without replacing the input read's match or coordinates.
 * @author Brian Bushnell
 * @author Isla
 * @contributor Sucrose
 * @date June 2025
 */
public class VarHelper{

	/*--------------------------------------------------------------*/
	/*----------------      Header Generation       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Formats caller-supplied statistics and VAR column names.
	 * Computes readLengthAvg as bases/max(1, reads), floors the displayed mincount
	 * at one, and includes the reference line only for nonnull ref. The final column
	 * line has no trailing newline; extended columns depend on Var.extendedText.
	 * @param properPairRate Fraction of reads that mapped as proper pairs
	 * @param totalQualityAvg Average base quality across all processed bases
	 * @param mapqAvg Average mapping quality across all reads
	 * @param rarity Minimum variant frequency threshold
	 * @param minAlleleFraction Minimum allele fraction for variant calling
	 * @param minAlleleDepth Allele-count threshold, displayed as at least one
	 * @param ploidy Expected ploidy level for organism
	 * @param reads Total number of reads processed
	 * @param pairs Number of paired reads in sequencing
	 * @param properPairs Number of properly paired reads
	 * @param bases Total number of bases processed
	 * @param ref Reference genome file path
	 * @return New header string, with internal newlines but no final newline
	 */
	public static String toVarHeader(double properPairRate, double totalQualityAvg, double mapqAvg, double rarity, double minAlleleFraction, int minAlleleDepth, int ploidy,
			long reads, long pairs, long properPairs, long bases, String ref){
		StringBuilder sb=new StringBuilder();

		final double readLengthAvg=bases/Tools.max(1.0, reads);
		sb.append("#fileformat\tVar_"+Var.varFormat+"\n");
		sb.append("#BBMapVersion\t"+Shared.BBTOOLS_VERSION_STRING+"\n");
		sb.append("#ploidy\t"+ploidy+"\n");
		sb.append(Tools.format("#rarity\t%.5f\n", rarity));
		sb.append(Tools.format("#minAlleleFraction\t%.4f\n", minAlleleFraction));
		sb.append("#mincount\t"+Tools.max(1, minAlleleDepth)+"\n");
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

	/** Returns the compact diagnostic column header, without a newline.
	 * Unlike toVarHeader, it omits flagged and extended columns and uses nearVarCount. */
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
	 * Formats a VCFv4.2-labelled BBTools header with its fixed INFO/FORMAT definitions,
	 * caller statistics and contigs in map order. Optionally includes reference and
	 * network metadata. Scaffold names are truncated at first whitespace when asked;
	 * other metadata is appended directly, without general escaping or validation.
	 * The column line always includes FORMAT, appends a sample only when nonnull,
	 * and has no final newline; the writer adds that newline.
	 * @param properPairRate Fraction of properly paired reads
	 * @param totalQualityAvg Average base quality score
	 * @param mapqAvg Average mapping quality score
	 * @param rarity Minimum variant frequency threshold
	 * @param minAlleleFraction Minimum allele fraction for calling
	 * @param minAlleleDepth Allele-count threshold, displayed as at least one
	 * @param ploidy Expected organism ploidy
	 * @param reads Total reads processed
	 * @param pairs Total paired reads
	 * @param properPairs Total properly paired reads
	 * @param bases Total bases processed
	 * @param ref Reference genome file path
	 * @param map Nonnull scaffold registry, read without copying or mutation
	 * @param sampleName Sample column name, or null to omit it
	 * @param trimWhitespace Whether to truncate names at their first whitespace
	 * @param net Optional network supplying its filename and cutoff
	 * @return New header string with no final newline
	 */
	public static String toVcfHeader(double properPairRate, double totalQualityAvg, double mapqAvg,
			double rarity, double minAlleleFraction, int minAlleleDepth, int ploidy, long reads, long pairs,
			long properPairs, long bases, String ref, ScafMap map, String sampleName,
			boolean trimWhitespace, CellNet net){
		StringBuilder sb=new StringBuilder();

		// VCF format and version information
		sb.append("##fileformat=VCFv4.2\n");
		sb.append("##BBMapVersion="+Shared.BBTOOLS_VERSION_STRING+"\n");
		sb.append("##ploidy="+ploidy+"\n");
		sb.append(Tools.format("##rarity=%.5f\n", rarity));
		sb.append(Tools.format("##minallelefraction=%.5f\n", minAlleleFraction));
		sb.append("##mincount="+Tools.max(1, minAlleleDepth)+"\n");
		sb.append("##reads="+reads+"\n");
		sb.append("##pairedReads="+pairs+"\n");
		sb.append("##properlyPairedReads="+properPairs+"\n");
		sb.append(Tools.format("##readLengthAvg=%.3f\n", (bases/Tools.max(reads, 1.0))));
		sb.append(Tools.format("##properPairRate=%.5f\n", properPairRate));
		sb.append(Tools.format("##totalQualityAvg=%.3f\n", totalQualityAvg));
		sb.append(Tools.format("##mapqAvg=%.3f\n", mapqAvg));
		if(ref!=null){sb.append("##reference="+ref+"\n");}
		// Neural network provenance (only when a net was used for scoring)
		if(net!=null){
			if(net.fname!=null){sb.append("##NNnetwork="+net.fname+"\n");}
			sb.append(Tools.format("##NNcutoff=%.6f\n", net.cutoff));
		}

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
			sb.append("##INFO=<ID=SCR,Number=1,Type=Float,Description=\"Composite (non-NN) variant score, phred-scaled\">\n");
			sb.append("##INFO=<ID=NNS,Number=1,Type=Float,Description=\"Raw neural-network score (0=worst, cutoff-centered, may exceed 1); present only when a network was used\">\n");
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
	 * Applies the output score transform: 2.5 times the Phred conversion of
	 * error value 1-0.998*score. Zero returns zero immediately; other boundary
	 * handling is delegated to QualityTools.probErrorToPhredDouble.
	 * @param score Internal variant score, normally between zero and one
	 * @return Scaled score, not an unadjusted Phred transform of 1-score
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
	 * Counts consecutive exact-byte matches to base in up to four positions on each
	 * side of pos, then adds one if both sides are nonempty. Does not compare the
	 * reference base at pos; this is a bounded flanking statistic, not a full run length.
	 * @param bases Nonnull reference sequence
	 * @param pos Zero-based position; out-of-range positions return zero
	 * @param base Allele base to match; not-fully-defined bases return zero
	 * @return Flanking count plus the two-sided bonus, from zero through nine
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

	/** Counts up to four exact matches to base, starting at pos and moving left.
	 * Includes pos; returns zero if pos is outside bases, its byte differs, or base
	 * is not fully defined. Requires a nonnull reference array. */
	public static int homopolymerCountLeft(final byte[] bases, final int pos, final byte base){
		if(pos<0 || pos>=bases.length || bases[pos]!=base){return 0;}
		if(!AminoAcid.isFullyDefined(base)){return 0;}

		int count=0;
		for(int i=pos, lim=Tools.max(0, pos-3); i>=lim; i--){
			if(bases[i]==base){count++;}
			else{break;}
		}
		return count;
	}

	/** Counts up to four exact matches to base, starting at pos and moving right.
	 * Includes pos; returns zero if pos is outside bases, its byte differs, or base
	 * is not fully defined. Requires a nonnull reference array. */
	public static int homopolymerCountRight(final byte[] bases, final int pos, final byte base){
		//pos can equal bases.length for an insertion whose right flank is the contig end
		//(caller passes stop+1); guard it like homopolymerCountSub before indexing bases[pos].
		if(pos<0 || pos>=bases.length || bases[pos]!=base){return 0;}
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

	/** Counts leading C operations in a nonempty expanded or valid count-coded match.
	 * A first symbol other than C returns zero; this is not an alphabet validator. */
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

	/** Counts trailing C bytes in an expanded match. The retained count-coded branch
	 * can overcount because its digit flag is already set before numeric accumulation;
	 * the current Var.toVars caller expands matches before reaching this helper. */
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
		//TODO: Probable bug - hasDigit is already true, so the C count prefixes the digits (C8 becomes 18).
		// Current junction caller supplies expanded matches; validate before reviving this branch (VAR2-012).
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

	/** Counts leading C bytes in an expanded match, returning zero for an empty array.
	 * Retained helper; the active left counter also handles count-coded matches. */
	private static int countLeftClipLong(byte[] longmatch){
		for(int i=0; i<longmatch.length; i++){
			if(longmatch[i]!='C'){return i;}
		}
		return longmatch.length;
	}

	/** Counts trailing C bytes in an expanded match, returning zero for an empty array. */
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
	 * Creates up to two junction observations from terminal clips meeting minClip.
	 * With neither end clipped, calls SoftClipper with mutation disabled. Otherwise
	 * uses the original clipping and coordinates. Does not replace r.match, r.start
	 * or r.stop; new variants accumulate evidence from r.
	 * Left calls use the last clipped read base at start+leftClip; right calls use
	 * the first right-clipped base at stop-rightClip+1. Both span one reference base.
	 * @param r Reference-oriented read with nonempty expanded match and bases;
	 * r.samline is needed when a variant accumulates evidence
	 * @param sl Currently unused; evidence reads r.samline instead
	 * @param scafnum Scaffold number assigned to the new variants
	 * @param containsVars Currently unused
	 * @param minClip Positive minimum terminal clip length in read bases
	 * @return New list of zero to two variants; null when pre-existing clips are both
	 * below threshold, or if the clipper were to return null
	 */
	static ArrayList<Var> toJunctions(Read r, SamLine sl, final int scafnum,
			boolean containsVars, final int minClip){
		final byte[] match0=r.match, match;
		final byte[] bases=r.bases;
		final int start, stop;

		// Count clipping in original alignment
		int leftClip=countLeftClip(match0), rightClip=countRightClip(match0);

		if(leftClip==0 && rightClip==0){
			// Try soft-clipping only when neither end already has C operations.
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
