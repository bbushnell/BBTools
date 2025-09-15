package var2;

import align2.MSA;
import dna.AminoAcid;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.SamLine;
import stream.SiteScore;

/**
 * Realigns reads using Multiple Sequence Alignment (MSA) to improve alignments
 * from non-affine-gap aligners. Particularly useful for reads with indels longer
 * than 1bp that were poorly aligned by simpler alignment algorithms.
 * 
 * Performs glocal alignment with padding around the original alignment region
 * and only retains realignments that improve the alignment score.
 * 
 * @author Brian Bushnell
 * @contributor Isla Winglet
 */
public class Realigner {
	
	/**
	 * Creates a Realigner with default parameters.
	 */
	public Realigner(){
		this(defaultMaxrows, defaultColumns, defaultPadding, defaultMsaType);
	}
	
	/**
	 * Creates a Realigner with specified MSA parameters.
	 * 
	 * @param maxrows_ Maximum number of rows for MSA matrix
	 * @param columns_ Maximum number of columns for MSA matrix  
	 * @param padding_ Padding bases to add around alignment region
	 * @param msaType_ MSA algorithm type string
	 */
	public Realigner(int maxrows_, int columns_, int padding_, String msaType_){
		maxrows=maxrows_;
		columns=columns_;
		padding=padding_;
		msaType=msaType_;
		msa=MSA.makeMSA(maxrows, columns+2, msaType);
	}
	
	/**
	 * Attempts to realign a read using scaffold lookup from the static map.
	 * 
	 * @param r Read to realign
	 * @param sl Corresponding SAM line
	 * @param unclip Whether to attempt unclipping of terminal indels
	 * @return true if realignment was successful and improved the alignment
	 */
	public boolean realign(Read r, SamLine sl, final boolean unclip){
		if(!r.mapped() || sl.supplementary() || !sl.primary()){return false;}
		Scaffold scaf=map.getScaffold(sl.rnameS());
		assert(scaf!=null) : sl.rnameS();
		return realign(r, sl, scaf, unclip);
	}
	
	/**
	 * Attempts to realign a read using the provided scaffold.
	 * 
	 * @param r Read to realign
	 * @param sl Corresponding SAM line
	 * @param scaf Scaffold containing reference sequence
	 * @param unclip Whether to attempt unclipping of terminal indels
	 * @return true if realignment was successful and improved the alignment
	 */
	public boolean realign(final Read r, final SamLine sl, final Scaffold scaf, final boolean unclip){
		return realign(r, sl, scaf.bases, unclip);
	}
	
	/**
	 * Core realignment method. Evaluates whether a read needs realignment based on
	 * alignment quality indicators, then performs MSA-based realignment if beneficial.
	 * 
	 * Only realigns reads with sufficient numbers of mismatches, clips, or indels.
	 * Retains realignments only if they improve the alignment score.
	 * 
	 * @param r Read to realign
	 * @param sl Corresponding SAM line (will be modified if realignment succeeds)
	 * @param ref Reference sequence bytes
	 * @param unclip Whether to attempt unclipping of terminal indels
	 * @return true if realignment was successful and improved the alignment
	 */
	public boolean realign(final Read r, final SamLine sl, final byte[] ref, final boolean unclip){
		if(!r.mapped() || sl.supplementary() || !sl.primary()){return false;}
		
		// Analyze alignment quality to determine if realignment is worthwhile
		int[] mSCNID=r.countMatchSymbols();
		int sumBad=mSCNID[1]+mSCNID[4]+mSCNID[5], sumIndel=mSCNID[4]+mSCNID[5];
		if(mSCNID[2]>0){//Has clips - continue with realignment
		}else if(sumBad>3){//Has many mismatches/indels - continue
		}else if(sumIndel>1 || (sumIndel>0 && mSCNID[1]>1)){//Complex indel pattern - continue
		}else{return false;}

		// Skip realignment for high-quality alignments with minimal issues
		if(mSCNID[1]<3 && mSCNID[2]==0 && mSCNID[4]<2 && mSCNID[5]<2 && sumBad<3 && sumIndel<2){return false;}
		
		// Resize MSA matrix if needed for this read length
		if(r.length()+2>msa.maxColumns+2*padding){
			msa=MSA.makeMSA(msa.maxRows, r.length()+2+r.length()/4+2*padding, msaType);
		}
		if(r.length()+2>msa.maxRows){
			msa=MSA.makeMSA(r.length()+2+r.length()/4+2*padding, msa.maxColumns, msaType);
		}

		// Calculate clipping information
		final int leadingClip=SamLine.countLeadingClip(sl.cigar, true, false);
		final int trailingClip=SamLine.countTrailingClip(sl.cigar, true, false);
		final int totalClip=leadingClip+trailingClip;
		final boolean clipped=totalClip>0;
		
		// Define padded alignment region
		final int start=sl.start(true, false);
		final int stop=sl.stop(start, true, false);
		final int paddedStart=start-padding, paddedStop=stop+padding;
		final int len0=paddedStop-paddedStart+1;
		final int a=0, b=len0-1;
		if(len0>=columns){return false;}
		
		realignmentsAttempted++;
		
		SiteScore ss=null;
		final byte[] rbases=makeRbases(ref, start, stop, padding);
		byte[] qbases=r.bases;
		
		// Reverse complement query bases if on minus strand
		if(sl.strand()==Shared.MINUS){
			qbases=AminoAcid.reverseComplementBases(qbases);
		}
		
		assert(!r.shortmatch()); //Otherwise convert it or change the score function.
		final int score0=msa.score(r.match);
		final int minScore=clipped ? 1 : Tools.max(1, score0-1);
		
		// Perform limited MSA alignment
		final int[] max;
		max=msa.fillLimited(qbases, rbases, a, b, minScore, null);
		if(max==null){return false;}

		final int[] score=msa.score(qbases, rbases, a, b, max[0], max[1], max[2], false);
		if(score==null){return false;}
		realignmentsSucceeded++;
		
		// Only keep realignments that improve the score
		if(score[0]<=minScore || (score[0]<=score0 && !clipped)){return false;}
		realignmentsImproved++;
		
		// Create new SiteScore from MSA results
		ss=new SiteScore(1, r.strand(), score[1], score[2], 1, score[0]);
		ss.setSlowScore(ss.quickScore);
		ss.score=ss.quickScore;
		ss.match=msa.traceback(qbases, rbases, a, b, score[3], score[4], score[5], false);
		assert(ss.match!=null);
		
		SiteScore oldSS2=ss.clone(); //Debugging backup
		oldSS2.match=ss.match.clone();
		
		// Correct coordinates for padded reference
		ss.start=ss.start+paddedStart;
		ss.stop=ss.stop+paddedStart;
		
		// Handle clipping and unclipping of terminal indels
		{
			int clipped2=ss.clipTipIndels(ref.length);
			if(unclip && clipped2>0 && ref!=null && oldSS2.match[oldSS2.match.length-1]=='Y'){
				ss.unclip(qbases, ref);
			}
		}
		
		// Reject alignments with problematic match symbols
		if(ss.matchContainsXY()){return false;}
		if(ss.matchContainsAB()){return false;}
		
		realignmentsRetained++;
		
		// Update read coordinates and match string
		r.start=ss.start;
		r.stop=ss.stop;
		r.match=ss.match;
		
		// Update SAM line with new alignment
		sl.pos=Tools.max(0, r.start)+1;
		sl.cigar=SamLine.toCigar14(r.match, r.start, r.stop, ref.length, qbases);
		sl.optional=null;
		
		// Convert to short match format for efficiency
		r.match=SamLine.cigarToShortMatch_old(sl.cigar, true);
		r.setShortMatch(true);
		r.toLongMatchString(true);
		
		return true;
	}
	
	/**
	 * Creates a padded reference sequence segment for MSA alignment.
	 * Extends the reference region with padding bases, using 'N' for out-of-bounds positions.
	 * 
	 * @param bases Full reference sequence
	 * @param start Original alignment start position
	 * @param stop Original alignment stop position  
	 * @param padding Number of bases to pad on each side
	 * @return Padded reference segment
	 */
	private static byte[] makeRbases(byte[] bases, int start, int stop, int padding) {
		final int start2=start-padding, stop2=stop+padding;
		byte[] out=new byte[stop2-start2+1];
		for(int opos=0, bpos=start2; opos<out.length; opos++, bpos++){
			byte b=(bpos<0 || bpos>=bases.length ? (byte)'N' : bases[bpos]);
			out[opos]=b;
		}
		return out;
	}
	
	/**
	 * Returns the MSA object used for alignment.
	 * 
	 * @return MSA instance
	 */
	public MSA msa(){return msa;}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Number of realignment attempts */
	long realignmentsAttempted=0;
	/** Number of successful MSA alignments */
	long realignmentsSucceeded=0;
	/** Number of realignments actually retained */
	long realignmentsRetained=0;  
	/** Number of realignments that improved scores */
	long realignmentsImproved=0;
	
	/** Maximum rows for MSA matrix */
	private int maxrows=602;
	/** Maximum columns for MSA matrix */
	private int columns=2000;
	/** Padding bases around alignment region */
	private int padding=100;
	/** MSA algorithm type */
	private String msaType;
	/** MSA instance for performing alignments */
	private MSA msa;
	
	/*--------------------------------------------------------------*/
	/*----------------       Static Fields          ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Default maximum rows for MSA matrix */
	public static int defaultMaxrows=603;
	/** Default maximum columns for MSA matrix */
	public static int defaultColumns=2000;
	/** Default padding around alignment region */
	public static int defaultPadding=200;
	/** Default MSA algorithm type */
	public static String defaultMsaType="MultiStateAligner11ts";
	/** Scaffold map for reference sequence lookup (initialized by CallVariants) */
	public static ScafMap map=null;
}