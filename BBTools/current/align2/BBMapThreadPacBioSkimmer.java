package align2;

import java.util.ArrayList;
import java.util.Arrays;

import bloom.BloomFilter;
import dna.AminoAcid;
import dna.Data;
import jgi.CoveragePileup;
import shared.Shared;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.SamLine;
import stream.SiteScore;

/**
 * Based on MapTestThread11f
 * Designed to skim and retain all sites above a threshold.
 * @author Brian Bushnell
 * @date Jul 10, 2012
 *
 */
public final class BBMapThreadPacBioSkimmer extends AbstractMapThread{

	static final int ALIGN_COLUMNS=BBIndexPacBioSkimmer.ALIGN_COLUMNS;
	static final int ALIGN_ROWS=4020;
	
	

	/** Don't trim for local alignments unless at least this many bases will be clipped */
	private final int LOCAL_ALIGN_TIP_LENGTH=1;
	/** Range is 0-1; a lower number makes trimming more aggressive */
	private final float LOCAL_ALIGN_MATCH_POINT_RATIO=0.75f;
	
	/** Ratio of the points for a match of a single base needed to declare unambiguous */
	public final float CLEARZONE_RATIOP=1.5f;
	public final float CLEARZONE_RATIO1=2.2f;
	public final float CLEARZONE_RATIO1b=2.8f;
	public final float CLEARZONE_RATIO1c=4.8f;
	public final float CLEARZONE_RATIO3=8f;
	/** Max allowed number of sites within 1 edit (excluding primary site) */
	public final int CLEARZONE_LIMIT1e=4;
	//public final int CLEARZONE1e;
	public final int CLEARZONEP;
	public final int CLEARZONE1;
	public final int CLEARZONE1b;
	public final int CLEARZONE1c;
	//public final int CLEARZONE1e;
	public final int CLEARZONE3;
	public final float INV_CLEARZONE3;
	public final float CLEARZONE1b_CUTOFF=0.92f;
	public final float CLEARZONE1c_CUTOFF=0.82f;
	
	public final BBIndexPacBioSkimmer index;
	
	
	private static int MIN_TRIM_SITES_TO_RETAIN_SINGLE=2;
	private static int MIN_TRIM_SITES_TO_RETAIN_PAIRED=1;
	
	/** TODO - perhaps I can rewrite cz3 to penalize reads that map similarly to more than the expected number of places */
	public static final boolean USE_CLEARZONE3=false;

	private static int EXPECTED_SITES=1;
	public static void setExpectedSites(int x){
		EXPECTED_SITES=x;
		MIN_TRIM_SITES_TO_RETAIN_SINGLE=Tools.max((int)(EXPECTED_SITES*4)+1, MIN_TRIM_SITES_TO_RETAIN_SINGLE);
		MIN_TRIM_SITES_TO_RETAIN_PAIRED=Tools.max((int)(EXPECTED_SITES*4)+1, MIN_TRIM_SITES_TO_RETAIN_PAIRED);
		MAX_TRIM_SITES_TO_RETAIN=EXPECTED_SITES*40+80;
	}
	
	@Override
	public final int ALIGN_COLUMNS(){return ALIGN_COLUMNS;}
	@Override
	public final int ALIGN_ROWS(){return ALIGN_ROWS;}
	@Override
	public final int maxReadLength(){return ALIGN_ROWS-1;}
	@Override
	final AbstractIndex index(){return index;}
	@Override
	final int CLEARZONE1(){return CLEARZONE1;}
	
	public BBMapThreadPacBioSkimmer(ConcurrentReadInputStream cris_, int keylen_,
			CoveragePileup pileup_, boolean SMITH_WATERMAN_, int THRESH_, int minChrom_,
			int maxChrom_, float keyDensity_, float maxKeyDensity_, float minKeyDensity_, int maxDesiredKeys_,
			boolean REMOVE_DUPLICATE_BEST_ALIGNMENTS_, boolean SAVE_AMBIGUOUS_XY_,
			float MINIMUM_ALIGNMENT_SCORE_RATIO_, boolean TRIM_LIST_, boolean MAKE_MATCH_STRING_, boolean QUICK_MATCH_STRINGS_,
			ConcurrentReadOutputStream outStream_, ConcurrentReadOutputStream outStreamMapped_, ConcurrentReadOutputStream outStreamUnmapped_, ConcurrentReadOutputStream outStreamBlack_,
			int SLOW_ALIGN_PADDING_, int SLOW_RESCUE_PADDING_, boolean DONT_OUTPUT_UNMAPPED_READS_, boolean DONT_OUTPUT_BLACKLISTED_READS_,
			int MAX_SITESCORES_TO_PRINT_, boolean PRINT_SECONDARY_ALIGNMENTS_,
			boolean REQUIRE_CORRECT_STRANDS_PAIRS_, boolean SAME_STRAND_PAIRS_, boolean KILL_BAD_PAIRS_, boolean RCOMP_MATE_,
			boolean PERFECTMODE_, boolean SEMIPERFECTMODE_, boolean FORBID_SELF_MAPPING_, int TIP_DELETION_SEARCH_RANGE_,
			boolean AMBIGUOUS_RANDOM_, boolean AMBIGUOUS_ALL_, int KFILTER_, float IDFILTER_, boolean TRIM_LEFT_, boolean TRIM_RIGHT_, boolean UNTRIM_, float TRIM_QUAL_, int TRIM_MIN_LEN_,
			boolean LOCAL_ALIGN_, boolean RESCUE_, boolean STRICT_MAX_INDEL_, String MSA_TYPE_, BloomFilter bloomFilter_){
		
		super(cris_,
				outStream_, outStreamMapped_, outStreamUnmapped_, outStreamBlack_,
				pileup_, SMITH_WATERMAN_, LOCAL_ALIGN_, REMOVE_DUPLICATE_BEST_ALIGNMENTS_,
				AMBIGUOUS_RANDOM_, AMBIGUOUS_ALL_, TRIM_LEFT_, TRIM_RIGHT_, UNTRIM_, TRIM_QUAL_, TRIM_MIN_LEN_, THRESH_,
				minChrom_, maxChrom_, KFILTER_, IDFILTER_, KILL_BAD_PAIRS_, SAVE_AMBIGUOUS_XY_,
				REQUIRE_CORRECT_STRANDS_PAIRS_,
				SAME_STRAND_PAIRS_, RESCUE_, STRICT_MAX_INDEL_, SLOW_ALIGN_PADDING_, SLOW_RESCUE_PADDING_,
				MSA_TYPE_, keylen_, PERFECTMODE_, SEMIPERFECTMODE_, FORBID_SELF_MAPPING_, RCOMP_MATE_,
				MAKE_MATCH_STRING_, DONT_OUTPUT_UNMAPPED_READS_, DONT_OUTPUT_BLACKLISTED_READS_, PRINT_SECONDARY_ALIGNMENTS_,
				QUICK_MATCH_STRINGS_, MAX_SITESCORES_TO_PRINT_, MINIMUM_ALIGNMENT_SCORE_RATIO_,
				keyDensity_, maxKeyDensity_, minKeyDensity_, maxDesiredKeys_,
				BBIndexPacBioSkimmer.MIN_APPROX_HITS_TO_KEEP, BBIndexPacBioSkimmer.USE_EXTENDED_SCORE,
				BBIndexPacBioSkimmer.BASE_HIT_SCORE, BBIndexPacBioSkimmer.USE_AFFINE_SCORE, BBIndexPacBioSkimmer.MAX_INDEL, TRIM_LIST_, TIP_DELETION_SEARCH_RANGE_, bloomFilter_);
		
		assert(SLOW_ALIGN_PADDING>=0);
		assert(!(RCOMP_MATE/* || FORBID_SELF_MAPPING*/)) : "RCOMP_MATE: TODO";
		
		if(SLOW_ALIGN || MAKE_MATCH_STRING){
//			msa=MSA.makeMSA(ALIGN_ROWS, ALIGN_COLUMNS, MSA_TYPE);
//			POINTS_MATCH=msa.POINTS_MATCH();
//			POINTS_MATCH2=msa.POINTS_MATCH2();
			CLEARZONE1=(int)(CLEARZONE_RATIO1*POINTS_MATCH2);
			CLEARZONE1b=(int)(CLEARZONE_RATIO1b*POINTS_MATCH2);
			CLEARZONE1c=(int)(CLEARZONE_RATIO1c*POINTS_MATCH2);
			CLEARZONEP=(int)(CLEARZONE_RATIOP*POINTS_MATCH2);
			CLEARZONE3=PENALIZE_AMBIG ? (int)(CLEARZONE_RATIO3*POINTS_MATCH2) : 0;
//			CLEARZONE1e=(int)(2*POINTS_MATCH2-POINTS_MATCH-msa.POINTS_SUB())+1;
		}else{
//			POINTS_MATCH=70;
//			POINTS_MATCH2=100;
//			msa=null;
			CLEARZONE1=0;
			CLEARZONE1b=0;
			CLEARZONE1c=0;
			CLEARZONEP=0;
			CLEARZONE3=0;
//			CLEARZONE1e=0;
		}
		INV_CLEARZONE3=(CLEARZONE3==0 ? 0 : 1f/CLEARZONE3);
		
		index=new BBIndexPacBioSkimmer(KEYLEN, minChrom, maxChrom, KFILTER, msa);
	}
	
	
	@Override
	public int trimList(ArrayList<SiteScore> list, boolean retainPaired, int maxScore, boolean specialCasePerfect, int minSitesToRetain, int maxSitesToRetain){
		if(list==null || list.size()==0){return -99999;}
		if(list.size()==1){return list.get(0).score;}

		boolean b=(list.size()>=minSitesToRetain);
		
		final int highestScore;
		if(USE_AFFINE_SCORE){
			
			Tools.trimSiteList(list, .10f, retainPaired, true, EXPECTED_SITES<3 ? EXPECTED_SITES : EXPECTED_SITES*2, maxSitesToRetain);
			
			highestScore=Tools.trimSiteList(list, .25f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
			if(highestScore==maxScore && specialCasePerfect){
				Tools.trimSiteList(list, .94f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
				if(list.size()>8){Tools.trimSiteList(list, .99f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
				return highestScore;
			}
			
			final int mstr2=(minSitesToRetain<=1 ? 1 : minSitesToRetain+1);
			final int mstr3=(minSitesToRetain+1)/2;
			
			assert(!b || list.size()>=EXPECTED_SITES);
			
			int N=(EXPECTED_SITES+3)/4;

			if(list.size()>2*N){
				Tools.trimSiteListByMax(list, (int)(.12f*maxScore), retainPaired, true, mstr3, maxSitesToRetain);
			}
			
			if(list.size()>3*N){
				Tools.trimSiteList(list, .28f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
				Tools.trimSiteListByMax(list, (int)(.14f*maxScore), retainPaired, true, mstr3, maxSitesToRetain);
			}
			
			if(list.size()>4*N){
				Tools.trimSiteList(list, .32f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
				Tools.trimSiteListByMax(list, (int)(.16f*maxScore), retainPaired, true, mstr3, maxSitesToRetain);
			}
			
			if(list.size()>5*N){
				Tools.trimSiteList(list, .36f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
				Tools.trimSiteListByMax(list, (int)(.17f*maxScore), retainPaired, true, mstr3, maxSitesToRetain);
			}
			
			if(list.size()>6*N){
				Tools.trimSiteList(list, .38f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
				Tools.trimSiteListByMax(list, (int)(.18f*maxScore), retainPaired, true, mstr3, maxSitesToRetain);
			}
			
			if(list.size()>8*N){
				Tools.trimSiteList(list, .40f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
				Tools.trimSiteListByMax(list, (int)(.20f*maxScore), retainPaired, true, mstr3, maxSitesToRetain);
			}
//			//			System.out.print(", "+list.size());
			if(list.size()>12*N){
				Tools.trimSiteList(list, .45f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
				Tools.trimSiteListByMax(list, (int)(.22f*maxScore), retainPaired, true, mstr3, maxSitesToRetain);
			}
//			//			System.out.print(", "+list.size());
			if(list.size()>16*N){
				Tools.trimSiteList(list, .50f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
				Tools.trimSiteListByMax(list, (int)(.24f*maxScore), retainPaired, true, mstr3, maxSitesToRetain);
			}
//			//			System.out.print(", "+list.size());
			if(list.size()>20*N){
				Tools.trimSiteList(list, .55f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
				Tools.trimSiteListByMax(list, (int)(.26f*maxScore), retainPaired, true, mstr3, maxSitesToRetain);
			}
			
			if(list.size()>24*N){
				Tools.trimSiteList(list, .60f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
				Tools.trimSiteListByMax(list, (int)(.28f*maxScore), retainPaired, true, mstr3, maxSitesToRetain);
			}
//			////			System.out.print(", "+list.size());
			if(list.size()>32*N){
				Tools.trimSiteList(list, .65f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
				Tools.trimSiteListByMax(list, (int)(.30f*maxScore), retainPaired, true, mstr3, maxSitesToRetain);
			}
//			////			System.out.print(", "+list.size());
			if(list.size()>40*N){Tools.trimSiteList(list, .70f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
//			//			System.out.print(", "+list.size());
			if(list.size()>48*N){Tools.trimSiteList(list, .75f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>64*N){Tools.trimSiteList(list, .80f, retainPaired, true, mstr2, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>80*N){Tools.trimSiteList(list, .85f, retainPaired, true, mstr2, maxSitesToRetain);}
			if(list.size()>96*N){Tools.trimSiteList(list, .90f, retainPaired, true, mstr2, maxSitesToRetain);}
			
			assert(!b || list.size()>=EXPECTED_SITES);
			

		}else if(BBIndexPacBioSkimmer.USE_EXTENDED_SCORE){
			highestScore=Tools.trimSiteList(list, .75f, retainPaired, true, minSitesToRetain, maxSitesToRetain);
			
			if(list.size()>8){Tools.trimSiteList(list, .8f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>16){Tools.trimSiteList(list, .85f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>24){Tools.trimSiteList(list, .90f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>36){Tools.trimSiteList(list, .92f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			////			System.out.print(", "+list.size());
			if(list.size()>40){Tools.trimSiteList(list, .94f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			////			System.out.print(", "+list.size());
			if(list.size()>48){Tools.trimSiteList(list, .96f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			////			System.out.print(", "+list.size());
			if(list.size()>56){Tools.trimSiteList(list, .97f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>64){Tools.trimSiteList(list, .98f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>80){Tools.trimSiteList(list, .99f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			

		}else{
			//			System.out.print("\n\nSize:\t"+list.size());


			highestScore=Tools.trimSiteList(list, .6f, retainPaired, true, minSitesToRetain, maxSitesToRetain);

			if(list.size()>12){Tools.trimSiteList(list, .65f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>16){Tools.trimSiteList(list, .7f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>24){Tools.trimSiteList(list, .74f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>28){Tools.trimSiteList(list, .8f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>32){Tools.trimSiteList(list, .85f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
			if(list.size()>48){Tools.trimSiteList(list, .90f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			////			System.out.print(", "+list.size());
			//			if(list.size()>40){Tools.trimSiteList(list, .95f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			////			System.out.print(", "+list.size());
			//			if(list.size()>48){Tools.trimSiteList(list, .96f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			////			System.out.print(", "+list.size());
			//			if(list.size()>56){Tools.trimSiteList(list, .97f, retainPaired, true, minSitesToRetain, maxSitesToRetain);}
			//			System.out.print(", "+list.size());
		}
		
		return highestScore;
	}
	
	
	@Override
	public void scoreSlow(final ArrayList<SiteScore> list, final byte[] basesP, final byte[] basesM,
			final int maxSwScore, final int maxImperfectSwScore){
		
		final int minMsaLimit;
		if(PAIRED){
			minMsaLimit=-CLEARZONE1e+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE*maxSwScore);
		}else{
			minMsaLimit=-CLEARZONE1e+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO*maxSwScore);
		}
		assert(Read.CHECKSITES(list, basesP, basesM, -1));
		
		int minMatch=Tools.max(-300, minMsaLimit-CLEARZONE3); //Score must exceed this to generate quick match string
		if(verbose){
			System.err.println("Slow-scoring.  maxSwScore="+maxSwScore+", maxImperfectSwScore="+maxImperfectSwScore+", minMsaLimit="+minMsaLimit+", minMatch="+minMatch);
		}
		for(int i=0; i<list.size(); i++){
			final SiteScore ss=list.get(i);
			assert(ss.lengthsAgree());
			final byte[] bases=(ss.strand==Shared.PLUS ? basesP : basesM);
			
			if(SEMIPERFECTMODE){
				assert(ss.stop-ss.start==bases.length-1);
				assert(ss.semiperfect);
			}
			
			if(verbose){System.err.println("\nSlow-scoring "+ss);}
			if(ss.stop-ss.start!=bases.length-1){
				assert(ss.stop-ss.start>bases.length-1) : bases.length+", "+ss.toText();
				assert(!ss.semiperfect) : "\n"+bases.length+", "+ss.toText()+", "+ss.perfect+", "+ss.semiperfect+", "+maxSwScore+"\n"+new String(basesP)+"\n";
				ss.setSlowScore(0);
				ss.semiperfect=false;
				ss.perfect=false;
			}
			
			final int swscoreNoIndel=ss.slowScore;
			int[] swscoreArray=null;
			
			boolean clipped=true, setLimits=false;
			if(swscoreNoIndel<maxImperfectSwScore && !ss.semiperfect){
				if(verbose && ss.stop-ss.start>4000){
					System.err.println(ss.toText());
					System.err.println(list.size());
					System.err.println();
				}
				
				int expectedLen=GapTools.calcGrefLen(ss);
				if(verbose){System.err.println("expectedLen="+expectedLen);}
				if(expectedLen>=EXPECTED_LEN_LIMIT){
					 //TODO: Alternately, I could kill the site.
					ss.setStop(ss.start+Tools.min(basesP.length+40, EXPECTED_LEN_LIMIT));
					if(verbose){System.err.println("expectedLen="+expectedLen+"; ss="+ss);}
				}
				
				int pad=SLOW_ALIGN_PADDING;
				final int minscore=Tools.max(swscoreNoIndel, minMsaLimit);
				final int minscore2=Tools.max(swscoreNoIndel-MSA.MIN_SCORE_ADJUST, minMsaLimit);
				if(verbose){System.err.println("Sent to msa with start="+ss.start+", stop="+ss.stop+", pad="+pad+", limit="+minscore+", gaps="+GapTools.toString(ss.gaps));}
				swscoreArray=msa.fillAndScoreLimited(bases, ss, pad, minscore);
				if(verbose){System.err.println("Received "+Arrays.toString(swscoreArray));}
				
				if(swscoreArray!=null && swscoreArray.length>6 && (swscoreArray[3]+swscoreArray[4]+expectedLen<EXPECTED_LEN_LIMIT)){
					int[] oldArray=swscoreArray.clone();
					assert(swscoreArray.length==8);
					int extraPadLeft=swscoreArray[6];
					int extraPadRight=swscoreArray[7];
					
					if(verbose){
						System.err.println("msa returned "+Arrays.toString(swscoreArray)+", re-running.");
						System.err.println("Added extra padding: "+ss.toText()+", "+Arrays.toString(oldArray));
					}
					
					ss.setLimits(ss.start-extraPadLeft, ss.stop+extraPadRight);
					pad=SLOW_ALIGN_PADDING+EXTRA_PADDING;
					if(verbose){System.err.println("Sent to msa with start="+ss.start+", stop="+ss.stop+", pad="+pad+", limit="+minscore+", gaps="+GapTools.toString(ss.gaps));}
					swscoreArray=msa.fillAndScoreLimited(bases, ss, pad, minscore);
					
					if(verbose){System.err.println("Result of extra padding: "+ss.toText()+", "+Arrays.toString(swscoreArray));}
					if(swscoreArray==null || swscoreArray[0]<oldArray[0]){
						if(verbose){
							System.err.println("Result was inferior.");
						}
						swscoreArray=oldArray;
					}
				}
				assert(ss.lengthsAgree());
				if(verbose){
					System.err.println(QUICK_MATCH_STRINGS+", "+(swscoreArray==null ? "null" : (swscoreArray.length+", "+swscoreArray[0]+" >=? "+minscore)));
					System.err.println("start="+ss.start+", stop="+ss.stop+", len="+ss.mappedLength());
				}
				if(QUICK_MATCH_STRINGS && swscoreArray!=null && swscoreArray.length==6 && swscoreArray[0]>=minscore2 && (PRINT_SECONDARY_ALIGNMENTS || (USE_SS_MATCH_FOR_PRIMARY && swscoreArray[0]>minMatch))){
					if(verbose){System.err.println("Generating match string.");}
					assert(swscoreArray.length==6) : swscoreArray.length;
					assert(swscoreArray[0]>=minscore2) : "\n"+Arrays.toString(swscoreArray)+"\n"+minscore+"\n"+minMatch;
					ss.match=msa.traceback(bases, Data.getChromosome(ss.chrom).array, ss.start-pad, ss.stop+pad, swscoreArray[3], swscoreArray[4], swscoreArray[5], ss.gaps!=null);
					if(ss.match!=null){
						assert(ss.pairedScore<1 || (ss.slowScore<=0 && ss.pairedScore>ss.quickScore ) || ss.pairedScore>ss.slowScore); //123
						ss.setLimits(swscoreArray[1], swscoreArray[2]);
						setLimits=true;
						assert(ss.lengthsAgree());
						clipped=ss.fixXY(bases, true, msa);
						assert(ss.pairedScore<1 || (ss.slowScore<=0 && ss.pairedScore>ss.quickScore ) || ss.pairedScore>ss.slowScore); //123
						clipped=ss.clipTipIndels(bases, basesM, 4, 10, msa) || clipped;
						assert(ss.pairedScore<1 || (ss.slowScore<=0 && ss.pairedScore>ss.quickScore ) || ss.pairedScore>ss.slowScore); //123
						assert(ss.lengthsAgree());
					}
				}else{
					ss.match=null;
				}
			}
			if(swscoreArray!=null && !setLimits){
				if(verbose){System.err.println("msa returned "+Arrays.toString(swscoreArray));}
				ss.setSlowScore(swscoreArray[0]);
				ss.setLimits(swscoreArray[1], swscoreArray[2]);
				assert(ss.lengthsAgree());
			}else{
				assert(swscoreNoIndel<=maxSwScore) : swscoreNoIndel+", "+maxImperfectSwScore+", "+maxSwScore+", "+new String(basesP);
				assert(clipped || swscoreNoIndel==-1 || msa.scoreNoIndels(bases, ss.chrom, ss.start)==swscoreNoIndel) :
					setLimits+", "+clipped+", "+(swscoreArray==null)+", "+
					swscoreNoIndel+" != "+msa.scoreNoIndels(bases, ss.chrom, ss.start)+"\n"+
					ss.toText()+"\n"+(ss.stop-ss.start)+", "+bases.length; //Slow
			}
			assert(ss.lengthsAgree());
			ss.setScore(ss.slowScore);
			minMatch=Tools.max(minMatch, ss.slowScore);
//			minMsaLimit=Tools.max(minMsaLimit, ss.slowScore-CLEARZONE3);
			assert(ss.slowScore<=maxSwScore);
			assert(!(ss.perfect && ss.slowScore<maxSwScore));
			ss.perfect=(ss.slowScore==maxSwScore);
			if(ss.perfect){ss.semiperfect=true;}
			else if(!ss.semiperfect){ss.setPerfect(bases);}
			
			if(verbose){System.err.println(" -> "+ss);}
		}
		
	}
	
	
	/** {group of correct hit (or -1), size of correct group, number of groups,
	 * number of elements, correctScore, maxScore, size of top group, num correct, firstElementCorrect,
	 * firstElementCorrectLoose, firstGroupCorrectLoose} */
	
	@Override
	public void calcStatistics1(final Read r, final int maxSwScore, final int maxPossibleQuickScore){
		if(OUTPUT_PAIRED_ONLY && r.mate!=null && !r.paired() && (r.mapped() || r.mate.mapped())){r.clearPairMapping();}
		final Read r2=r.mate;
		
		if(!r.mapped() && (r2==null || !r2.mapped())){
			bothUnmapped++;
			bothUnmappedBases+=r.length();
			if(r2!=null){
				bothUnmapped++;
				bothUnmappedBases+=r2.length();
			}
		}

		/* {number of correct (loose) sites, number of incorrect (loose) sites, number incorrect sites before last correct site,
		 * number of sites, correctScore, maxScore, firstElementCorrect, firstElementCorrectLoose, position of first correct element (or -1),
		 * sizeOfTopGroup, numTopCorrect} */
		int[] correctness=calcCorrectnessSkimmer(r, THRESH);
		
		int numCorrect=correctness[0];
		int numIncorrect=correctness[1];
		int numIncorrectPrior=correctness[2];
		int numSites=correctness[3];
		int correctScore=correctness[4];
		int topScore=correctness[5];
		boolean firstElementCorrect=(correctness[6]==1);
		boolean firstElementCorrectLoose=(correctness[7]==1);
		int positionOfFirstCorrect=correctness[8];
		int sizeOfTopGroup=correctness[9];
		int numTopCorrect=correctness[10];

		final int len1=r.length();
		final int len2=(r2==null ? 0 : r.length());
		
		assert(numSites==numCorrect+numIncorrect) : numSites+", "+numCorrect+", "+numIncorrect+", "+r.numSites();
		assert(numSites==r.numSites());

		totalNumCorrect1+=numCorrect;
		totalNumIncorrect1+=numIncorrect;
		totalNumIncorrectPrior1+=numIncorrectPrior;
		if(numCorrect>=EXPECTED_SITES){
			totalNumCapturedAllCorrect1++;
//			assert(numCorrect==EXPECTED_CORRECT_SITES) : numCorrect +", "+EXPECTED_CORRECT_SITES+", "+r.list;
			assert(r.sites.size()>=EXPECTED_SITES) : numCorrect;
			if(numTopCorrect==numCorrect){
				totalNumCapturedAllCorrectTop1++;
				if(numCorrect==numSites){
					totalNumCapturedAllCorrectOnly1++;
				}
			}
		}
		
		assert(numSites>0 == r.mapped());
		if(numSites>0){
			
			if(r.match!=null){
				int[] errors=r.countErrors(SamLine.INTRON_LIMIT);
				matchCountM1+=errors[0];
				matchCountS1+=errors[1];
				matchCountD1+=errors[2];
				matchCountI1+=errors[3];
				matchCountN1+=errors[4];
				
				readCountS1+=(errors[1]>0 ? 1 : 0);
				readCountD1+=(errors[2]>0 ? 1 : 0);
				readCountI1+=(errors[3]>0 ? 1 : 0);
				readCountN1+=(errors[4]>0 ? 1 : 0);
				readCountSplice1+=(errors[5]>0 ? 1 : 0);
				readCountE1+=((errors[1]>0 || errors[2]>0 || errors[3]>0)? 1 : 0);
			}
			
			
			mappedRetained1++;
			if(r.rescued()){
				if(r.strand()==Shared.PLUS){
					rescuedP1++;
				}else{
					rescuedM1++;
				}
			}
			if(r.paired()){
				numMated++;
				int inner;
				int outer;
				if(r.start<=r.mate.start){
					inner=r.mate.start-r.stop;
					outer=r.mate.stop-r.start;
				}else{
					inner=r.start-r.mate.stop;
					outer=r.stop-r.mate.start;
				}

				inner=Tools.min(MAX_PAIR_DIST, inner);
				inner=Tools.max(MIN_PAIR_DIST, inner);
				innerLengthSum+=inner;
				outerLengthSum+=outer;
				insertSizeSum+=(inner+r.length()+r.mateLength());
			}else if(r.mate!=null && r.mate.mapped()/*&& r.list!=null && r.list.size()>0*/){
				badPairs++;
			}
			
			if(r.perfect() || (maxSwScore>0 && r.topSite().slowScore==maxSwScore)){
				perfectMatch1++;
			}else if(SLOW_ALIGN){
				assert(r.topSite().slowScore<maxSwScore) : maxSwScore+"\t"+r.topSite().toText();
			}
			
			int foundSemi=0;
			for(SiteScore ss : r.sites){
				if(ss.perfect){
					perfectHitCount1++;
					assert(ss.semiperfect);
				}
				if(ss.semiperfect){
					semiPerfectHitCount1++;
					foundSemi=1;
				}
			}
			semiperfectMatch1+=foundSemi;
			if(foundSemi>0){semiperfectMatchBases1+=len1;}
			
			if(firstElementCorrect){
				if(r.strand()==Shared.PLUS){firstSiteCorrectP1++;}
				else{firstSiteCorrectM1++;}
				if(r.paired()){firstSiteCorrectPaired1++;}
				else{firstSiteCorrectSolo1++;}
				if(r.rescued()){firstSiteCorrectRescued1++;}
			}else{
				firstSiteIncorrect1++;
//				System.out.println("********");
//				System.out.println(r.toText(false));
//				System.out.println(r.mate.toText(false));
			}
			
			if(firstElementCorrectLoose){
				firstSiteCorrectLoose1++;
			}else{
				firstSiteIncorrectLoose1++;
			}

			siteSum1+=numSites;
			topSiteSum1+=sizeOfTopGroup;

			if(topScore==maxPossibleQuickScore){perfectHit1++;}
			if(sizeOfTopGroup==1){uniqueHit1++;}

			if(numCorrect>0){
				
				if(r.strand()==Shared.PLUS){truePositiveP1++;}
				else{truePositiveM1++;}
				totalCorrectSites1+=numCorrect;

				if(positionOfFirstCorrect==0){
					correctUniqueHit1++;
				}else{
					correctLowHit1++;
				}

			}else{

				falsePositive1++;
//				System.out.println("********");
//				System.out.println(r.toText(false));
//				System.out.println(r.mate.toText(false));
			}
		}else if(maxPossibleQuickScore==-1){
			lowQualityReadsDiscarded1++;
			lowQualityBasesDiscarded1+=len1;
			r.setDiscarded(true);
		}else{
			noHit1++;
		}
	}
	
	
	@Override
	public void calcStatistics2(final Read r, final int maxSwScore, final int maxPossibleQuickScore){

		int[] correctness=calcCorrectnessSkimmer(r, THRESH);
		
		int numCorrect=correctness[0];
		int numIncorrect=correctness[1];
		int numIncorrectPrior=correctness[2];
		int numSites=correctness[3];
		int correctScore=correctness[4];
		int topScore=correctness[5];
		boolean firstElementCorrect=(correctness[6]==1);
		boolean firstElementCorrectLoose=(correctness[7]==1);
		int positionOfFirstCorrect=correctness[8];
		int sizeOfTopGroup=correctness[9];
		int numTopCorrect=correctness[10];

		final int len=r.length();
		
		totalNumCorrect2+=numCorrect;
		totalNumIncorrect2+=numIncorrect;
		totalNumIncorrectPrior2+=numIncorrectPrior;
		if(numCorrect>=EXPECTED_SITES){
			totalNumCapturedAllCorrect2++;
			if(numTopCorrect==numCorrect){
				totalNumCapturedAllCorrectTop2++;
				if(numCorrect==numSites){
					totalNumCapturedAllCorrectOnly2++;
				}
			}
		}
		
		if(numSites>0){
			
			if(r.match!=null){
				int[] errors=r.countErrors(SamLine.INTRON_LIMIT);
				matchCountM2+=errors[0];
				matchCountS2+=errors[1];
				matchCountD2+=errors[2];
				matchCountI2+=errors[3];
				matchCountN2+=errors[4];
				
				readCountS2+=(errors[1]>0 ? 1 : 0);
				readCountD2+=(errors[2]>0 ? 1 : 0);
				readCountI2+=(errors[3]>0 ? 1 : 0);
				readCountN2+=(errors[4]>0 ? 1 : 0);
				readCountSplice2+=(errors[5]>0 ? 1 : 0);
				readCountE2+=((errors[1]>0 || errors[2]>0 || errors[3]>0)? 1 : 0);
			}
			
			mappedRetained2++;
			if(r.rescued()){
				if(r.strand()==Shared.PLUS){
					rescuedP2++;
				}else{
					rescuedM2++;
				}
			}
			
			if(r.perfect() || (maxSwScore>0 && r.topSite().slowScore==maxSwScore)){
				perfectMatch2++;
			}else if(SLOW_ALIGN){
				assert(r.topSite().slowScore<maxSwScore) : maxSwScore+"\t"+r.topSite().toText();
			}
			
			int foundSemi=0;
			for(SiteScore ss : r.sites){
				if(ss.perfect){
					perfectHitCount2++;
					assert(ss.semiperfect);
				}
				if(ss.semiperfect){
					semiPerfectHitCount2++;
					foundSemi=1;
				}
			}
			semiperfectMatch2+=foundSemi;
			if(foundSemi>0){semiperfectMatchBases2+=len;}

			if(firstElementCorrect){
				if(r.strand()==Shared.PLUS){firstSiteCorrectP2++;}
				else{firstSiteCorrectM2++;}
				if(r.paired()){firstSiteCorrectPaired2++;}
				else{firstSiteCorrectSolo2++;}
				if(r.rescued()){firstSiteCorrectRescued2++;}
			}else{
				firstSiteIncorrect2++;
//				System.out.println("********");
//				System.out.println(r.toText(false));
//				System.out.println(r.mate.toText(false));
			}
			
			if(firstElementCorrectLoose){
				firstSiteCorrectLoose2++;
			}else{
				firstSiteIncorrectLoose2++;
			}

			siteSum2+=numSites;
			topSiteSum2+=sizeOfTopGroup;

			if(topScore==maxPossibleQuickScore){perfectHit2++;}
			if(sizeOfTopGroup==1){uniqueHit2++;}

			if(numCorrect>0){

				if(r.strand()==Shared.PLUS){truePositiveP2++;}
				else{truePositiveM2++;}
				totalCorrectSites2+=numCorrect;

				if(positionOfFirstCorrect==0){
					correctUniqueHit2++;
				}else{
					correctLowHit2++;
				}

			}else{

				falsePositive2++;
				//				System.out.println("********");
				//				System.out.println(r.toText(false));
				//				System.out.println(r.mate.toText(false));
			}
		}else if(maxPossibleQuickScore==-1){
			lowQualityReadsDiscarded2++;
		}else{
			noHit2++;
		}
	}
	
	@Override
	public void processRead(final Read r, final byte[] basesM){
		if(idmodulo>1 && r.numericID%idmodulo!=1){return;}
		final byte[] basesP=r.bases;
		
//		System.err.print(" rd#"+r.numericID+" ");
//		if(r.numericID==25967){
//			verbose=true;
//			msa.verbose=true;
//			GapTools.verbose=true;
//			index.verbose=true;
//			tcr.verbose=true;
//		}
		
		if(verbose){System.err.println("\nProcessing "+r);}
		readsUsed1++;
		
		final int maxPossibleQuickScore=quickMap(r, basesM);
		if(verbose){System.err.println("\nQuick Map: \t"+r.sites);}
		
		if(maxPossibleQuickScore<0){
			r.sites=null;
			lowQualityReadsDiscarded1++;
			lowQualityBasesDiscarded1+=basesP.length;
			r.setDiscarded(true);
			return;
		}
		initialSiteSum1+=r.numSites();
		if(verbose){System.err.println("\ninitialSiteSum1: "+initialSiteSum1);}
		
		int maxSwScore=0;
		int maxImperfectSwScore=0;

		if(SLOW_ALIGN || USE_AFFINE_SCORE){
			maxSwScore=msa.maxQuality(r.length());
			maxImperfectSwScore=msa.maxImperfectScore(r.length());
		}
		
		if(TRIM_LIST && r.numSites()>1){
			if(MIN_TRIM_SITES_TO_RETAIN_SINGLE>1){Shared.sort(r.sites);}
			int highestQuickScore=trimList(r.sites, false, maxSwScore, true, MIN_TRIM_SITES_TO_RETAIN_SINGLE, MAX_TRIM_SITES_TO_RETAIN);
		}
		postTrimSiteSum1+=r.numSites();
		if(verbose){System.err.println("\nAfter trim: \t"+r.sites);}
		
		assert(Read.CHECKSITES(r, basesM));
		
		
		if(SLOW_ALIGN && r.numSites()>0){
			Tools.subsumeOverlappingSites(r.sites, true, false);
			int numNearPerfectScores=scoreNoIndels(r, basesP, basesM, maxSwScore, maxImperfectSwScore);

			Shared.sort(r.sites); //Puts higher scores first to better trigger the early exit based on perfect scores
			
//			int numPerfectScores=0;
//			if(numNearPerfectScores>0){
//				for(SiteScore ss : r.list){
//					if(ss.perfect){numPerfectScores++;}
//					else{break;}
//				}
//			}
			
			if(verbose){
				System.err.println("\nAfter scoreNoIndels: \t"+r.sites);
			}
			
			if(numNearPerfectScores<1){
				if(FIND_TIP_DELETIONS){findTipDeletions(r, basesP, basesM, maxSwScore, maxImperfectSwScore);}
			}
			
			if(verbose){
				System.err.println("\nAfter findTipDeletions: \t"+r.sites);
			}
			
			//TODO: This causes problems with perfect matches that are mapped to areas longer than the read length
			//***Above note should be resolved now, but needs to be verified.
			
			if(numNearPerfectScores<r.sites.size()){
				scoreSlow(r.sites, basesP, basesM, maxSwScore, maxImperfectSwScore);
			}
			
			if(STRICT_MAX_INDEL){
				int removed=removeLongIndels(r.sites, index.MAX_INDEL);
				if(r.numSites()==0){r.clearMapping();}
			}
			
			if(verbose){System.err.println("\nAfter scoreSlow: \t"+r.sites);}
			assert(Read.CHECKSITES(r, basesM, false));
		}


		if(r.numSites()>0){
			mapped1++;
			try {
				Tools.mergeDuplicateSites(r.sites, true, true);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				throw new RuntimeException("\n\n"+r.toText(false)+"\n\n");
			}
			try {
				Tools.removeOverlappingSites(r.sites, true);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				throw new RuntimeException("\n\n"+r.toText(false)+"\n\n");
			}
			Shared.sort(r.sites);
			assert(Read.CHECKSITES(r, basesM));
		}
		if(verbose){System.err.println("\nAfter merge: \t"+r.sites);}
		
		if(r.numSites()>1){
			SiteScore ss1=r.topSite();
			SiteScore ss2=r.sites.get(1);
			//Ensure no duplicates
			assert(ss1.chrom!=ss2.chrom || ss1.strand!=ss2.strand || ss1.start!=ss2.start || ss1.stop!=ss2.stop) : r.toText(false);
		}
		assert(Read.CHECKSITES(r, basesM));
		
		if(r.numSites()>1){
			assert(r.topSite().score==r.topSite().slowScore);
		}
		
		if(SLOW_ALIGN || USE_AFFINE_SCORE){r.setPerfectFlag(maxSwScore);}
		
		if(r.numSites()>1){
			final int clearzone=r.perfect() ? CLEARZONEP :
				r.topSite().score>=(int)(maxSwScore*CLEARZONE1b_CUTOFF) ? CLEARZONE1 :
					(r.topSite().score>=(int)(maxSwScore*CLEARZONE1c_CUTOFF) ? CLEARZONE1b : CLEARZONE1c);
			final int numBestSites1=Tools.countTopScores(r.sites, clearzone);
			if(numBestSites1>1){
				//Ambiguous alignment
				assert(r.sites.size()>1);
				boolean b=processAmbiguous(r.sites, true, AMBIGUOUS_TOSS, clearzone, false);
				r.setAmbiguous(b);
			}
		}
		
		if(SLOW_ALIGN || USE_AFFINE_SCORE){
			Tools.removeLowQualitySitesUnpaired2(r.sites, maxSwScore, MINIMUM_ALIGNMENT_SCORE_RATIO, EXPECTED_SITES);
			if(r.numSites()>1){
				SiteScore a=r.topSite();
				SiteScore b=r.sites.get(1);
				assert(a.score>=b.score);
				assert(a.score>=(int)(maxSwScore*MINIMUM_ALIGNMENT_SCORE_RATIO)) : a;
				assert(b.score>=(int)(maxSwScore*MINIMUM_ALIGNMENT_SCORE_RATIO)) : "\n"+a+"\t"+b+"\n"+(maxSwScore*MINIMUM_ALIGNMENT_SCORE_RATIO);
			}
		}
		if(verbose){System.err.println("\nAfter removal: \t"+r.sites);}
		
		if(r.numSites()==0){r.sites=null;r.mapScore=0;}
		r.setFromTopSite(AMBIGUOUS_RANDOM, true, MAX_PAIR_DIST);
		assert(Read.CHECKSITES(r, basesM));
		
		assert(r.gaps==null || r.gaps[0]==r.start && r.gaps[r.gaps.length-1]==r.stop);
		assert(r.sites==null || r.mapScore>0) : r.sites+", "+r.mapScore+"\n"+r;
		
		if(r.numSites()>1){
			assert(r.topSite().score==r.topSite().slowScore) : "\n"+r.toText(false)+"\n";
			assert(r.topSite().score==r.mapScore) : "\n"+r.toText(false)+"\n";
		}
		
		if(verbose){System.err.println("C: "+r);}
		
		//***$
		if(MAKE_MATCH_STRING && r.numSites()>0){
			if(USE_SS_MATCH_FOR_PRIMARY && r.topSite().match!=null){
				r.match=r.topSite().match;
			}else{
				if(r.sites.size()>1){
					assert(r.topSite().score>=r.sites.get(1).score) : "\n"+r.topSite().toText()+"\t<\t"+r.sites.get(1).toText()+"\n"+r.toText(false)+"\n";
				}
				int mapScore=r.mapScore;

				assert(r.mate!=null || r.numSites()==0 || r.topSite().score==r.mapScore) : "\n"+r.toText(false)+"\n";

				if(verbose){System.err.println("D: "+r);}
				
				{
					boolean firstIter=true;
					do{//
						if(!firstIter){
							Shared.sort(r.sites);
							r.setFromTopSite(AMBIGUOUS_RANDOM, true, MAX_PAIR_DIST);
						}
						genMatchString(r, basesP, basesM, maxImperfectSwScore, maxSwScore, true, true);
						assert(r.mate!=null || r.numSites()==0 || r.topSite().score==r.mapScore) : "\n"+r.toText(false)+"\n";
//						TODO: Fix this; it should never happen.
//						if(mapScore>r.mapScore){
//							System.err.println("genMatchString reduced mapping score: "+mapScore+" -> "+r.mapScore+" in read "+r.numericID);
//						}
						if(STRICT_MAX_INDEL && hasLongIndel(r.match, index.MAX_INDEL)){
							SiteScore ss=r.topSite();
							ss.score=r.mapScore=Tools.min(ss.score, -9999);
							ss.setSlowPairedScore(ss.score, ss.score);
						}
						r.topSite().score=r.topSite().slowScore;
						firstIter=false;
					}while(r.sites.size()>1 && r.topSite().score<r.sites.get(1).score);
				}
				
				if(r.numSites()>1){
					assert(r.topSite().score==r.topSite().slowScore) : "\n"+r.toText(false)+"\n";
					assert(r.topSite().score==r.mapScore) : "\n"+r.toText(false)+"\n";
					assert(Read.CHECKSITES(r, basesM));
				}

				if(verbose){System.err.println("E: "+r);}
			}
		}
		if(verbose){System.err.println("\nAfter match: \t"+r.sites);}
		
		if(r.numSites()>1){
			assert(r.topSite().score==r.topSite().slowScore) : "\n"+r.toText(false)+"\n";
			assert(r.topSite().score==r.mapScore) : "\n"+r.toText(false)+"\n";
			removeDuplicateBestSites(r);
		}
		if(r.numSites()>0){r.topSite().match=r.match;}
		
		
		
		if(r.sites!=null && r.mapScore<=0){
			if(!STRICT_MAX_INDEL && !Shared.anomaly){
				System.err.println("Note: Read "+r.id+" failed cigar string generation and will be marked as unmapped.\t"+(r.match==null)+"\t"+r.mapScore+"\t"+r.topSite()+"\t"+new String(r.bases));
				if(MSA.bandwidth>0 || MSA.bandwidthRatio>0 || MSA.flatMode){Shared.anomaly=true;}
			}
			r.mapScore=0;
			r.setMapped(false);
			r.sites=null;
		}
		
		
		
		//This block is to prevent an assertion from firing.  Generally caused by alignment being lost during match generation.
		//TODO: Fix cause.
		if(r.mapScore>0 && r.sites==null){
			if(!Shared.anomaly){System.err.println("Anomaly: mapScore>0 and list==null.\n"+r+"\n");}
			Shared.anomaly=true;
			r.clearMapping();
		}else if(r.mapScore<=0 && r.sites!=null){
			if(BANDWIDTH<1){
				if(!Shared.anomaly){System.err.println("Anomaly1: mapScore<=0 and list!=null.\n"+r+"\n");}
				Shared.anomaly=true;
			}
			r.clearMapping();
		}
		assert(r.sites==null || r.mapScore>0) :
			"\nmapScore = "+r.mapScore+"\nread = "+r.toText(false)+"\nscore thresh = "+(-100+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO*maxSwScore))+"\n"+
			"msa unlimited return = "+Arrays.toString(msa.fillAndScoreLimited(r.strand()==Shared.PLUS ? r.bases :
			AminoAcid.reverseComplementBases(r.bases), r.topSite(), Tools.max(SLOW_ALIGN_PADDING, 10), 0))+"\n"+
			"msa limited return = "+Arrays.toString(msa.fillAndScoreLimited(r.strand()==Shared.PLUS ? r.bases :
			AminoAcid.reverseComplementBases(r.bases), r.topSite(), Tools.max(SLOW_ALIGN_PADDING, 10), (-100+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO*maxSwScore))))+"\n\n"+
			"msa vert limit: "+msa.showVertLimit()+"\n\nmsa horz limit: "+msa.showHorizLimit()+"\n\n";
		
//		assert(r.list==null || r.mapScore>0) : r.mapScore+"\n"+r.list==null ? "null" : r.list.toString();
		
		if(USE_CLEARZONE3 && (CLEARZONE3>CLEARZONE1 || CLEARZONE3>CLEARZONEP)){
			boolean changed=applyClearzone3(r, CLEARZONE3, INV_CLEARZONE3);
			if(changed){
				int minScore=(int)(maxSwScore*MINIMUM_ALIGNMENT_SCORE_RATIO);
				if(r.mapScore<minScore){
					r.setAmbiguous(true);
					
					if(AMBIGUOUS_TOSS){
						r.sites.clear();
						r.clearSite();
						r.setMapped(false);
					}
				}
			}
		}
		
		if(r.ambiguous() && AMBIGUOUS_TOSS){r.sites=null; r.clearSite(); r.setMapped(false);}
		
		if(r.mapped() && r.numSites()>1 && PRINT_SECONDARY_ALIGNMENTS){
			ensureMatchStringsOnSiteScores(r, basesM, maxImperfectSwScore, maxSwScore);
			assert(Read.CHECKSITES(r, basesM));
		}
		
		assert(checkTopSite(r));
		if(r.mapped() && (LOCAL_ALIGN || r.containsXYC())){
			msa.toLocalAlignment(r, r.topSite(), basesM, r.containsXYC() ? 1 : LOCAL_ALIGN_TIP_LENGTH, LOCAL_ALIGN_MATCH_POINT_RATIO);
			assert(Read.CHECKSITES(r, basesM));
		}
		
		if(r.numSites()==0 || (!r.ambiguous() && r.mapScore<maxSwScore*MINIMUM_ALIGNMENT_SCORE_RATIO)){
			r.clearMapping();
		}
		postFilterRead(r, basesM, maxImperfectSwScore, maxSwScore);
		if(MAKE_MATCH_STRING){ensureMatchStringOnPrimary(r, basesM, maxImperfectSwScore, maxSwScore);}
		
		if(PENALIZE_AMBIG){
			int penalty=calcTipScorePenalty(r, maxSwScore, 7);
			applyScorePenalty(r, penalty);
		}
		
		if(CALC_STATISTICS){
			calcStatistics1(r, maxSwScore, maxPossibleQuickScore);
		}
	}
	
	
	/** Returns number of perfect pairs */
	@Override
	public int pairSiteScoresInitial(Read r, Read r2, boolean trim){
		
		if(r.numSites()<1 || r2.numSites()<1){return 0;}
		
		SiteScore.PCOMP.sort(r.sites);
		SiteScore.PCOMP.sort(r2.sites);
		
		for(SiteScore ss : r.sites){ss.setPairedScore(0);}
		for(SiteScore ss : r2.sites){ss.setPairedScore(0);}
		
//		ArrayList<SiteScorePair> pairs=new ArrayList<SiteScorePair>(Tools.min(8, Tools.min(r.list.size(), r2.list.size())));

		int maxPairedScore1=-1;
		int maxPairedScore2=-1;
		
		
//		for(SiteScore ss : r.list){
//			System.out.println(ss.toText());
//		}
		
//		int i=0, j=0;
		final int ilimit=r.sites.size()-1;
		final int jlimit=r2.sites.size()-1;
		final int maxReadLen=Tools.max(r.length(), r2.length());
		
//		final int outerDistLimit=MIN_PAIR_DIST+r.length()+r2.length();
		final int outerDistLimit=(Tools.max(r.length(), r2.length())*(OUTER_DIST_MULT))/OUTER_DIST_DIV;//-(SLOW_ALIGN ? 100 : 0);
		final int innerDistLimit=MAX_PAIR_DIST;//+(FIND_TIP_DELETIONS ? TIP_DELETION_SEARCH_RANGE : 0);
		final int expectedFragLength=AVERAGE_PAIR_DIST+r.length()+r2.length();
		
		int numPerfectPairs=0;
		
		for(int i=0, j=0; i<=ilimit && j<=jlimit; i++){
			SiteScore ss1=r.sites.get(i);
			SiteScore ss2=r2.sites.get(j);
			
			while(j<jlimit && (ss2.chrom<ss1.chrom || (ss2.chrom==ss1.chrom && ss1.start-ss2.stop>innerDistLimit))){
				j++;
				ss2=r2.sites.get(j);
			}

			for(int k=j; k<=jlimit; k++){
				ss2=r2.sites.get(k);

				if(ss2.chrom>ss1.chrom){break;}
				if(ss2.start-ss1.stop>innerDistLimit){break;}

//				int dist=0;
//
//				if(ss1.start<=ss2.start){
//					dist=ss2.start-ss1.stop;
//				}else if(ss1.start>ss2.start){
//					dist=ss1.start-ss2.stop;
//				}
				
				
//				int innerdist=0;
//				int outerdist=0;
//
//				if(ss1.start<=ss2.start){
//					innerdist=ss2.start-ss1.stop;
//					outerdist=ss2.stop-ss1.start;
//				}else if(ss1.start>ss2.start){
//					innerdist=ss1.start-ss2.stop;
//					outerdist=ss1.stop-ss2.start;
//				}
				
				final int innerdist, outerdist;
				//assert(!SAME_STRAND_PAIRS) : "TODO";
				
				if(REQUIRE_CORRECT_STRANDS_PAIRS){
					if(ss1.strand!=ss2.strand){
						if(ss1.strand==Shared.PLUS){
							innerdist=ss2.start-ss1.stop;
							outerdist=ss2.stop-ss1.start;
						}else{
							innerdist=ss1.start-ss2.stop;
							outerdist=ss1.stop-ss2.start;
						}
					}else{
						if(ss1.start<=ss2.start){
							innerdist=ss2.start-ss1.stop;
							outerdist=ss2.stop-ss1.start;
						}else{
							innerdist=ss1.start-ss2.stop;
							outerdist=ss1.stop-ss2.start;
						}
					}
				}else{
					if(ss1.start<=ss2.start){
						innerdist=ss2.start-ss1.stop;
						outerdist=ss2.stop-ss1.start;
					}else{
						innerdist=ss1.start-ss2.stop;
						outerdist=ss1.stop-ss2.start;
					}
				}
				
				assert(outerdist>=innerdist);

				if(outerdist>=outerDistLimit && innerdist<=innerDistLimit){
					
					boolean strandOK=((ss1.strand==ss2.strand)==SAME_STRAND_PAIRS);

					if(strandOK || !REQUIRE_CORRECT_STRANDS_PAIRS){
						
						boolean paired1=false, paired2=false;
						
						int deviation=absdif(AVERAGE_PAIR_DIST, innerdist);

						final int pairedScore1;
						final int pairedScore2;
						if(strandOK){
//							pairedScore1=ss1.score+ss2.score/2;
//							pairedScore2=ss2.score+ss1.score/2;
							
							pairedScore1=ss1.score+1+Tools.max(1, ss2.score/2-(((deviation)*ss2.score)/(32*expectedFragLength+100)));
							pairedScore2=ss2.score+1+Tools.max(1, ss1.score/2-(((deviation)*ss1.score)/(32*expectedFragLength+100)));
						}else{//e.g. a junction
							pairedScore1=ss1.score+Tools.max(0, ss2.score/16);
							pairedScore2=ss2.score+Tools.max(0, ss1.score/16);
						}

						if(pairedScore1>ss1.pairedScore){
							paired1=true;
							ss1.setPairedScore(Tools.max(ss1.pairedScore, pairedScore1));
							maxPairedScore1=Tools.max(ss1.score, maxPairedScore1);
							//						System.out.println("Paired "+ss1.toText()+" with "+ss2.toText());
						}else{
							//						System.out.println(ss1.toText()+" already paired.");
						}
						if(pairedScore2>ss2.pairedScore){
							paired2=true;
							ss2.setPairedScore(Tools.max(ss2.pairedScore, pairedScore2));
							maxPairedScore2=Tools.max(ss2.score, maxPairedScore2);
						}
						
						if(paired1 && paired2 && outerdist>=maxReadLen && deviation<=expectedFragLength && ss1.perfect && ss2.perfect){
							numPerfectPairs++; //Lower bound.  Some perfect pairs may be the same.
						}
						
//						ss1.pairedScore=Tools.max(ss1.pairedScore, pairedScore1);
//						ss2.pairedScore=Tools.max(ss2.pairedScore, pairedScore2);
//						maxPairedScore1=Tools.max(ss1.score, maxPairedScore1);
//						maxPairedScore2=Tools.max(ss2.score, maxPairedScore2);
					}
				}
			}
			
		}
		
		
		
		for(SiteScore ss : r.sites){
			if(ss.pairedScore>ss.score){ss.score=ss.pairedScore;}
			else{assert(ss.pairedScore==0);}
//			ss.score=ss.pairedScore=Tools.max(ss.pairedScore, ss.score);
		}
		for(SiteScore ss : r2.sites){
			if(ss.pairedScore>ss.score){ss.score=ss.pairedScore;}
			else{assert(ss.pairedScore==0);}
//			ss.score=ss.pairedScore=Tools.max(ss.pairedScore, ss.score);
		}
		
		if(trim){
			if(numPerfectPairs>MIN_TRIM_SITES_TO_RETAIN_PAIRED){
//				System.out.print(".");
				Shared.sort(r.sites);
				Shared.sort(r2.sites);
				Tools.trimSitesBelowCutoff(r.sites, (int)(maxPairedScore1*.1f), true, true, MIN_TRIM_SITES_TO_RETAIN_PAIRED, MAX_TRIM_SITES_TO_RETAIN);
				Tools.trimSitesBelowCutoff(r2.sites, (int)(maxPairedScore2*.1f), true, true, MIN_TRIM_SITES_TO_RETAIN_PAIRED, MAX_TRIM_SITES_TO_RETAIN);
			}else{
				if(r.sites.size()>MIN_TRIM_SITES_TO_RETAIN_PAIRED){
					Shared.sort(r.sites);
					Tools.trimSitesBelowCutoff(r.sites, (int)(maxPairedScore1*.1f), true, true, MIN_TRIM_SITES_TO_RETAIN_PAIRED, MAX_TRIM_SITES_TO_RETAIN);
				}
				if(r2.sites.size()>MIN_TRIM_SITES_TO_RETAIN_PAIRED){
					Shared.sort(r2.sites);
					Tools.trimSitesBelowCutoff(r2.sites, (int)(maxPairedScore2*.1f), true, true, MIN_TRIM_SITES_TO_RETAIN_PAIRED, MAX_TRIM_SITES_TO_RETAIN);
				}
			}
		}
		
//		if(pairs.isEmpty()){return null;}
//
//		ArrayList<SiteScore> temp=new ArrayList<SiteScore>(Tools.max(r.list.size(), r2.list.size()));
//
//		for(SiteScore ss : r.list){
//			if(ss.score>maxPairedScore1){temp.add(ss);}
//		}
//		for(SiteScorePair ssp : pairs){
//			temp.add(ssp.a);
//		}
//		r.list.clear();
//		r.list.addAll(temp);
//
//		for(SiteScore ss : r2.list){
//			if(ss.score>maxPairedScore2){temp.add(ss);}
//		}
//		for(SiteScorePair ssp : pairs){
//			temp.add(ssp.b);
//		}
//		r2.list.clear();
//		r2.list.addAll(temp);
//
//		return pairs;
		
		return numPerfectPairs;
	}
	
	
	@Override
	public void processReadPair(final Read r, final byte[] basesM1, final byte[] basesM2){
		if(idmodulo>1 && r.numericID%idmodulo!=1){return;}
		final Read r2=r.mate;
		assert(r2!=null);
		final byte[] basesP1=r.bases, basesP2=r2.bases;
		final int len1=(basesP1==null ? 0 : basesP1.length), len2=(basesP2==null ? 0 : basesP2.length);
		
		readsUsed1++;
		readsUsed2++;

		final int maxPossibleQuickScore1=quickMap(r, basesM1);
		final int maxPossibleQuickScore2=quickMap(r2, basesM2);
		
		if(verbose){
			System.err.println("\nAfter quick map:\nRead1:\t"+r+"\nRead2:\t"+r.mate);
		}
		
		if(maxPossibleQuickScore1<0 && maxPossibleQuickScore2<0){
			r.sites=null;
			r2.sites=null;
			lowQualityReadsDiscarded1++;
			lowQualityBasesDiscarded1+=len1;
			r.setDiscarded(true);
			lowQualityReadsDiscarded2++;
			lowQualityBasesDiscarded2+=len2;
			r2.setDiscarded(true);
			return;
		}
		
		//Not really needed due to subsumption
		Tools.mergeDuplicateSites(r.sites, true, false);
		Tools.mergeDuplicateSites(r2.sites, true, false);
		
		initialSiteSum1+=r.numSites();
		initialSiteSum2+=r2.numSites();
		
		//TODO: Fix this.  This is a workaround for an assertion error counting the number of reads used.
		//Discards need to be tracked separately for each end.
//		if(maxPossibleQuickScore2<0){lowQualityReadsDiscarded--;}
		
		final int maxSwScore1=msa.maxQuality(len1);
		final int maxImperfectSwScore1=msa.maxImperfectScore(len1);
		final int maxSwScore2=msa.maxQuality(len2);
		final int maxImperfectSwScore2=msa.maxImperfectScore(len2);
		
		//TODO: POSSIBLY block pairing across N blocks.
		pairSiteScoresInitial(r, r2, TRIM_LIST);
		if(verbose){System.err.println("\nAfter initial pair:\nRead1:\t"+r+"\nRead2:\t"+r2);}
		
		if(TRIM_LIST){

			if(MIN_TRIM_SITES_TO_RETAIN_PAIRED>1){
				if(r.numSites()>MIN_TRIM_SITES_TO_RETAIN_PAIRED){
					Shared.sort(r.sites);
					trimList(r.sites, true, maxSwScore1, false, MIN_TRIM_SITES_TO_RETAIN_PAIRED, MAX_TRIM_SITES_TO_RETAIN);
				}
				if(r2.numSites()>MIN_TRIM_SITES_TO_RETAIN_PAIRED){
					Shared.sort(r2.sites);
					trimList(r2.sites, true, maxSwScore2, false, MIN_TRIM_SITES_TO_RETAIN_PAIRED, MAX_TRIM_SITES_TO_RETAIN);
				}
			}
		}
		postTrimSiteSum1+=r.numSites();
		postTrimSiteSum2+=r2.numSites();
		
		{//Reset score to non-paired score
			if(r.sites!=null){
				for(SiteScore ss : r.sites){
					assert(ss.slowScore<=ss.quickScore);
					ss.score=ss.quickScore;
				}
			}
			if(r2.sites!=null){
				for(SiteScore ss : r2.sites){
					assert(ss.slowScore<=ss.quickScore);
					ss.score=ss.quickScore;
				}
			}
		}
		
		if(verbose){System.err.println("\nAfter trim:\nRead1:\t"+r.sites+"\nRead2:\t"+r2.sites);}
		
//		assert(Read.CHECKSITES(r, basesM1) && Read.CHECKSITES(r2, basesM2));
		
		if(SLOW_ALIGN){
			
			if(r.numSites()>0){
				Tools.subsumeOverlappingSites(r.sites, true, false);
				int numNearPerfectScores1=scoreNoIndels(r, basesP1, basesM1, maxSwScore1, maxImperfectSwScore1);
				Shared.sort(r.sites); //Puts higher scores first to better trigger the early exit based on perfect scores
				
				if(numNearPerfectScores1<MIN_TRIM_SITES_TO_RETAIN_PAIRED){
					if(FIND_TIP_DELETIONS){findTipDeletions(r, basesP1, basesM1, maxSwScore1, maxImperfectSwScore1);}
				}
				
				//TODO:
				//Note scoreSlow can be skipped under this circumstance:
				//When rescue is disabled, numNearPerfectScores>0, and there are no paired sites.
				scoreSlow(r.sites, basesP1, basesM1, maxSwScore1, maxImperfectSwScore1);
				if(STRICT_MAX_INDEL){
					int removed=removeLongIndels(r.sites, index.MAX_INDEL);
					if(r.numSites()==0){r.clearMapping();}
				}
				Tools.mergeDuplicateSites(r.sites, true, true);
			}
			
			if(r2.numSites()>0){
				Tools.subsumeOverlappingSites(r2.sites, true, false);
				int numNearPerfectScores2=scoreNoIndels(r2, basesP2, basesM2, maxSwScore2, maxImperfectSwScore2);
				Shared.sort(r2.sites); //Puts higher scores first to better trigger the early exit based on perfect scores
				
				if(numNearPerfectScores2<MIN_TRIM_SITES_TO_RETAIN_PAIRED){
					if(FIND_TIP_DELETIONS){findTipDeletions(r2, basesP2, basesM2, maxSwScore2, maxImperfectSwScore2);}
				}
				
				scoreSlow(r2.sites, basesP2, basesM2, maxSwScore2, maxImperfectSwScore2);
				if(STRICT_MAX_INDEL){
					int removed=removeLongIndels(r2.sites, index.MAX_INDEL);
					if(r2.numSites()<1){r2.clearMapping();}
				}
				Tools.mergeDuplicateSites(r2.sites, true, true);
			}
			
			
			if(verbose){System.err.println("\nAfter slow align:\nRead1:\t"+r+"\nRead2:\t"+r2);}
			assert(Read.CHECKSITES(r, basesM1, false) && Read.CHECKSITES(r2, basesM2, false));
			
			if(DO_RESCUE){
				int unpaired1=0;
				int unpaired2=0;
				if(r.sites!=null){
					for(SiteScore ss : r.sites){
						assert(ss.pairedScore<1 || ss.pairedScore>ss.quickScore || ss.pairedScore>ss.slowScore) :
							"\n"+ss.toText()+"\n"+r.toText(false)+"\n";
						if(ss.pairedScore==0){unpaired1++;}
					}
				}
				if(r2.sites!=null){
					for(SiteScore ss : r2.sites){
						assert(ss.pairedScore<1 || ss.pairedScore>ss.quickScore || ss.pairedScore>ss.slowScore) :
							"\n"+ss.toText()+"\n"+r2.toText(false)+"\n";
						if(ss.pairedScore==0){unpaired2++;}
					}
				}
				
				if(unpaired1>0 && r.numSites()>0){
					Shared.sort(r.sites);
					Tools.removeLowQualitySitesPaired(r.sites, maxSwScore1, MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE, MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE);
					rescue(r, r2, basesP2, basesM2, Tools.min(MAX_PAIR_DIST, 2*AVERAGE_PAIR_DIST+100));
					Tools.mergeDuplicateSites(r2.sites, true, true);
				}
				
				if(unpaired2>0 && r2.numSites()>0){
					Shared.sort(r2.sites);
					Tools.removeLowQualitySitesPaired(r2.sites, maxSwScore2, MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE, MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE);
					rescue(r2, r, basesP1, basesM1, Tools.min(MAX_PAIR_DIST, 2*AVERAGE_PAIR_DIST+100));
					Tools.mergeDuplicateSites(r.sites, true, true);
				}

				postRescueSiteSum1+=r.numSites();
				postRescueSiteSum2+=r2.numSites();
				
//				if(r.list!=null){Shared.sort(r.list);}
//				if(r2.list!=null){Shared.sort(r2.list);}
//
//				Tools.removeLowQualitySites(r.list, maxSwScore1, MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE, MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE);
//				Tools.removeLowQualitySites(r2.list, maxSwScore2, MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE, MINIMUM_ALIGNMENT_SCORE_RATIO_PRE_RESCUE);
				
				if(verbose){System.err.println("\nAfter rescue:\nRead1:\t"+r+"\nRead2:\t"+r2);}
				assert(Read.CHECKSITES(r, basesM1, false) && Read.CHECKSITES(r2, basesM2, false));
			}
		}else{
			Tools.mergeDuplicateSites(r.sites, true, false);
			Tools.mergeDuplicateSites(r2.sites, true, false);
			if(verbose){System.err.println("\nAfter merge:\nRead1:\t"+r+"\nRead2:\t"+r2);}
			assert(Read.CHECKSITES(r, basesM1, false) && Read.CHECKSITES(r2, basesM2, false));
		}
		
		try {
			Tools.removeOverlappingSites(r.sites, true);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			throw new RuntimeException("\n\n"+r.toText(false)+"\n\n");
		}
		try {
			Tools.removeOverlappingSites(r2.sites, true);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			throw new RuntimeException("\n\n"+r2.toText(false)+"\n\n");
		}
		
		if(r.numSites()>1){Shared.sort(r.sites);}
		if(r2.numSites()>1){Shared.sort(r2.sites);}
		
		if(false){//This block is optional, but increases correctness by a tiny bit. (or maybe not!)
			if(SLOW_ALIGN || USE_AFFINE_SCORE){
				Tools.removeLowQualitySitesPaired(r.sites, maxSwScore1, MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED, MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED);
				Tools.removeLowQualitySitesPaired(r2.sites, maxSwScore2, MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED, MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED);
			}

			pairSiteScoresFinal(r, r2, false, false, MAX_PAIR_DIST, AVERAGE_PAIR_DIST, SAME_STRAND_PAIRS, REQUIRE_CORRECT_STRANDS_PAIRS, MAX_TRIM_SITES_TO_RETAIN);
			
			if(r.numSites()>1){Shared.sort(r.sites);}
			if(r2.numSites()>1){Shared.sort(r2.sites);}
		}
		
		if(SLOW_ALIGN || USE_AFFINE_SCORE){
//			assert(false) : "Change to removeLowQualitySitesPaired2";
			//TODO Verify correctness
			Tools.removeLowQualitySitesPaired2(r.sites, maxSwScore1, MINIMUM_ALIGNMENT_SCORE_RATIO, MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED, EXPECTED_SITES);
			Tools.removeLowQualitySitesPaired2(r2.sites, maxSwScore2, MINIMUM_ALIGNMENT_SCORE_RATIO, MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED, EXPECTED_SITES);
		}
		
		//TODO
//		assert(false) : "Only verified up to this point.";
		pairSiteScoresFinal(r, r2, false, true, MAX_PAIR_DIST, AVERAGE_PAIR_DIST, SAME_STRAND_PAIRS, REQUIRE_CORRECT_STRANDS_PAIRS, MAX_TRIM_SITES_TO_RETAIN);
		
		if(r.numSites()>0){
			mapped1++;
			Shared.sort(r.sites);
		}
		if(r2.numSites()>0){
			mapped2++;
			Shared.sort(r2.sites);
		}
		assert(Read.CHECKSITES(r, basesM1) && Read.CHECKSITES(r2, basesM2));
		
		if(SLOW_ALIGN || USE_AFFINE_SCORE){
			r.setPerfectFlag(maxSwScore1);
			r2.setPerfectFlag(maxSwScore2);
//			assert(Read.CHECKSITES(r, basesM1) && Read.CHECKSITES(r2, basesM2));
		}
		

		if(r.numSites()>1){
			final int clearzone=r.perfect() ? CLEARZONEP :
				r.topSite().score>=(int)(maxSwScore1*CLEARZONE1b_CUTOFF) ? CLEARZONE1 :
					(r.topSite().score>=(int)(maxSwScore1*CLEARZONE1c_CUTOFF) ? CLEARZONE1b : CLEARZONE1c);
			int numBestSites1=Tools.countTopScores(r.sites, clearzone);
			if(numBestSites1>1){
				//Ambiguous alignment
				assert(r.sites.size()>1);
				
				boolean b=processAmbiguous(r.sites, true, AMBIGUOUS_TOSS, clearzone, false);
				r.setAmbiguous(b);
			}
//			assert(Read.CHECKSITES(r, basesM1));
		}

		if(r2.numSites()>1){
			final int clearzone=r2.perfect() ? CLEARZONEP :
				r2.topSite().score>=(int)(maxSwScore2*CLEARZONE1b_CUTOFF) ? CLEARZONE1 :
					(r2.topSite().score>=(int)(maxSwScore2*CLEARZONE1c_CUTOFF) ? CLEARZONE1b : CLEARZONE1c);
			int numBestSites2=Tools.countTopScores(r2.sites, clearzone);
			if(numBestSites2>1){
				//Ambiguous alignment
				assert(r2.sites.size()>1);
				
				boolean b=processAmbiguous(r2.sites, false, AMBIGUOUS_TOSS, clearzone, false);
				r2.setAmbiguous(b);
			}
//			assert(Read.CHECKSITES(r2, basesM2));
		}
		if(verbose){System.err.println("\nAfter ambiguous removal:\nRead1:\t"+r+"\nRead2:\t"+r2);}
		
		if(r.numSites()>0 && r2.numSites()>0){
			SiteScore ss1=r.topSite();
			SiteScore ss2=r2.topSite();
			if(canPair(ss1, ss2, len1, len2, REQUIRE_CORRECT_STRANDS_PAIRS, SAME_STRAND_PAIRS, MAX_PAIR_DIST)){
				assert(SLOW_ALIGN ? ss1.pairedScore>ss1.slowScore : ss1.pairedScore>ss1.quickScore) :
					"\n"+ss1.toText()+"\n"+ss2.toText()+"\n"+r.toText(false)+"\n"+r2.toText(false)+"\n\n"+
						r.mapped()+", "+r.paired()+", "+r.strand()+", "+r.ambiguous()+"\n\n"+r2.mapped()+", "+r2.paired()+", "+r2.strand()+", "+r2.ambiguous()+"\n\n";
				assert(SLOW_ALIGN ? ss2.pairedScore>ss2.slowScore : ss2.pairedScore>ss2.quickScore) :
					"\n"+ss1.toText()+"\n"+ss2.toText()+"\n"+r.toText(false)+"\n"+r2.toText(false)+"\n\n";
				r.setPaired(true);
				r.mate.setPaired(true);
			}
		}

		if(r.numSites()==0){r.sites=null;r.mapScore=0;}
		if(r2.numSites()==0){r2.sites=null;r2.mapScore=0;}
		
		r.setFromTopSite(AMBIGUOUS_RANDOM, true, MAX_PAIR_DIST);
		r2.setFromTopSite(AMBIGUOUS_RANDOM, true, MAX_PAIR_DIST);
		if(KILL_BAD_PAIRS){
			if(r.isBadPair(REQUIRE_CORRECT_STRANDS_PAIRS, SAME_STRAND_PAIRS, MAX_PAIR_DIST)){
				int x=r.mapScore/len1;
				int y=r2.mapScore/len2;
				if(x>=y){
					r2.clearAnswers(false);
				}else{
					r.clearAnswers(false);
				}
			}
		}
		if(verbose){System.err.println("\nAfter bad pair removal:\nRead1:\t"+r+"\nRead2:\t"+r2);}
		
		assert(r.sites==null || r.mapScore>0) : r.mapScore+"\n"+r.toText(false)+"\n\n"+r2.toText(false)+"\n";
		assert(r2.sites==null || r2.mapScore>0) : r2.mapScore+"\n"+r.toText(false)+"\n\n"+r2.toText(false)+"\n";
		if(MAKE_MATCH_STRING){
			if(r.numSites()>0){
				if(USE_SS_MATCH_FOR_PRIMARY && r.topSite().match!=null){
					r.match=r.topSite().match;
				}else{
					genMatchString(r, basesP1, basesM1, maxImperfectSwScore1, maxSwScore1, false, false);
					
					if(STRICT_MAX_INDEL && r.mapped()){
						if(hasLongIndel(r.match, index.MAX_INDEL)){
							r.clearMapping();
							r2.setPaired(false);
						}
					}
				}
//				assert(Read.CHECKSITES(r, basesM1));
			}
			if(r2.numSites()>0){
				if(USE_SS_MATCH_FOR_PRIMARY && r2.topSite().match!=null){
					r2.match=r2.topSite().match;
				}else{
					genMatchString(r2, basesP2, basesM2, maxImperfectSwScore2, maxSwScore2, false, false);
					
					if(STRICT_MAX_INDEL && r2.mapped()){
						if(hasLongIndel(r2.match, index.MAX_INDEL)){
							r2.clearMapping();
							r.setPaired(false);
						}
					}
				}
//				assert(Read.CHECKSITES(r2, basesM2));
			}
		}
		
		assert(checkTopSite(r)); // TODO remove this
		if(verbose){
			System.err.println("\nFinal:\nRead1:\t"+r+"\nRead2:\t"+r2);
			if(r.match!=null && r.shortmatch()){r.toLongMatchString(false);}
			if(r2.match!=null && r2.shortmatch()){r2.toLongMatchString(false);}
		}
		
		//Block to prevent assertion from firing.  Generally caused by alignment being lost during match generation.  TODO: Fix cause.
		if(r.mapScore>0 && r.sites==null){
			if(!Shared.anomaly){System.err.println("Anomaly: mapScore>0 and list==null.\n"+r+"\n");}
			Shared.anomaly=true;
			r.clearMapping();
			r2.setPaired(false);
		}else if(r.mapScore<=0 && r.sites!=null){
			if(!STRICT_MAX_INDEL && !Shared.anomaly){System.err.println("Anomaly2: mapScore<=0 and list!=null.\n"+r+"\n");}
			Shared.anomaly=true;
			r.clearMapping();
			r2.setPaired(false);
		}
		assert(checkTopSite(r)); // TODO remove this
		//Block to prevent assertion from firing.  Generally caused by alignment being lost during match generation.  TODO: Fix cause.
		if(r2.mapScore>0 && r2.sites==null){
			if(!Shared.anomaly){System.err.println("Anomaly: mapScore>0 and list==null.\n"+r+"\n");}
			Shared.anomaly=true;
			r2.clearMapping();
			r.setPaired(false);
		}else if(r2.mapScore<=0 && r2.sites!=null){
			if(!STRICT_MAX_INDEL && !Shared.anomaly){System.err.println("Anomaly3: mapScore<=0 and list!=null.\n"+r+"\n");}
			Shared.anomaly=true;
			r2.clearMapping();
			r.setPaired(false);
		}
		
		assert(r.sites==null || r.mapScore>0) :
			r.mapScore+"\t"+r.sites+"\n"+(-100+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED*maxSwScore1))+"\n"+
			Arrays.toString(msa.fillAndScoreLimited(r.strand()==Shared.PLUS ? r.bases :
			AminoAcid.reverseComplementBases(r.bases), r.topSite(), Tools.max(SLOW_ALIGN_PADDING, 80), 0))+"\n"+
			Arrays.toString(msa.fillAndScoreLimited(r.strand()==Shared.PLUS ? r.bases :
			AminoAcid.reverseComplementBases(r.bases), r.topSite(), Tools.max(SLOW_ALIGN_PADDING, 80), (-100+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED*maxSwScore1))))+"\n\n"+
			msa.showVertLimit()+"\n\n"+msa.showHorizLimit()+"\n\n"+r+"\n\n"+r2+"\n\n";
		assert(r2.sites==null || r2.mapScore>0) :
			r2.mapScore+"\t"+r2.sites+"\n"+(-100+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED*maxSwScore2))+"\n"+
			Arrays.toString(msa.fillAndScoreLimited(r2.strand()==Shared.PLUS ? r2.bases :
			AminoAcid.reverseComplementBases(r2.bases), r2.topSite(), Tools.max(SLOW_ALIGN_PADDING, 80), 0))+"\n"+
			Arrays.toString(msa.fillAndScoreLimited(r2.strand()==Shared.PLUS ? r2.bases :
			AminoAcid.reverseComplementBases(r2.bases), r2.topSite(), Tools.max(SLOW_ALIGN_PADDING, 80), (-100+(int)(MINIMUM_ALIGNMENT_SCORE_RATIO_PAIRED*maxSwScore2))))+"\n\n"+
			msa.showVertLimit()+"\n\n"+msa.showHorizLimit()+"\n\n"+r+"\n\n"+r2+"\n\n";
		
		assert(!r.mapped() || !MAKE_MATCH_STRING || r.match!=null) : "Note that sometimes, VERY RARELY, match string generation fails.";
		assert(checkTopSite(r)); // TODO remove this
		removeDuplicateBestSites(r);
		removeDuplicateBestSites(r2);
		
		if(DYNAMIC_INSERT_LENGTH && numMated>1000 && r.paired()){
			AVERAGE_PAIR_DIST=(int)(innerLengthSum*1f/numMated);
		}
		assert(checkTopSite(r)); // TODO remove this
		if(r.ambiguous() && AMBIGUOUS_TOSS){
			if(r.sites!=null){r.sites=null;}
			r.clearSite();
			r.setMapped(false);
			r.setPaired(false);
			r2.setPaired(false);
		}else if(r.mapped() && r.numSites()>1 && PRINT_SECONDARY_ALIGNMENTS){
			ensureMatchStringsOnSiteScores(r, basesM1, maxImperfectSwScore1, maxSwScore1);
			assert(Read.CHECKSITES(r, basesM1));
		}
		if(r2.ambiguous() && AMBIGUOUS_TOSS){
			if(r2.sites!=null){r2.sites=null;}
			r2.clearSite();
			r2.setMapped(false);
			r.setPaired(false);
			r2.setPaired(false);
		}else if(r2.mapped() && r2.numSites()>1 && PRINT_SECONDARY_ALIGNMENTS){
			ensureMatchStringsOnSiteScores(r2, basesM2, maxImperfectSwScore2, maxSwScore2);
			assert(Read.CHECKSITES(r2, basesM2));
		}
//		assert(Read.CHECKSITES(r, basesM1) && Read.CHECKSITES(r2, basesM2));

		assert(checkTopSite(r));
		if(r.mapped() && (LOCAL_ALIGN || r.containsXYC())){
			final SiteScore ss=r.topSite();
			ss.match=r.match;
			msa.toLocalAlignment(r, ss, basesM1, r.containsXYC() ? 1 : LOCAL_ALIGN_TIP_LENGTH, LOCAL_ALIGN_MATCH_POINT_RATIO);
//			System.err.println("\n\n*********\n\n"+r+"\n\n*********\n\n");
//			assert(Read.CHECKSITES(r, basesM1)); //TODO: This can fail; see bug#0001
		}
		
		assert(checkTopSite(r2));
		if(r2.mapped() && (LOCAL_ALIGN || r2.containsXYC())){
			final SiteScore ss=r2.topSite();
			ss.match=r2.match;
			msa.toLocalAlignment(r2, ss, basesM2, r2.containsXYC() ? 1 : LOCAL_ALIGN_TIP_LENGTH, LOCAL_ALIGN_MATCH_POINT_RATIO);
//			assert(Read.CHECKSITES(r2, basesM2)); //TODO: This can fail; see bug#0001
		}
		
		postFilterRead(r, basesM1, maxImperfectSwScore1, maxSwScore1);
		postFilterRead(r2, basesM2, maxImperfectSwScore2, maxSwScore2);
		if(MAKE_MATCH_STRING){
			ensureMatchStringOnPrimary(r, basesM1, maxImperfectSwScore1, maxSwScore1);
			ensureMatchStringOnPrimary(r2, basesM2, maxImperfectSwScore2, maxSwScore2);
		}
		
		if(CALC_STATISTICS){
			calcStatistics1(r, maxSwScore1, maxPossibleQuickScore1);
			calcStatistics2(r2, maxSwScore2, maxPossibleQuickScore2);
		}
	}
	
	/** {number of correct (loose) sites, number of incorrect (loose) sites, number incorrect sites before last correct site,
	 * number of sites, correctScore, maxScore, firstElementCorrect, firstElementCorrectLoose, position of first correct element (or -1),
	 * sizeOfTopGroup, numTopCorrect} */
	protected static int[] calcCorrectnessSkimmer(Read r, int thresh){
		//assume sorted.
		ArrayList<SiteScore> ssl=r.sites;
		
		if(ssl==null || ssl.isEmpty()){
			return new int[] {0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0};
		}

		SiteScore original=r.originalSite;
		assert((original==null) != (r.synthetic()));
		if(original==null){
			original=ssl.get(0);
		}
		
		int numCorrect=0;
		int numIncorrect=0;
		int numIncorrectPrior=0; //Prior to expected number of correct sites being found
		int numIncorrectPrior2=0; //Prior to last correct element actually found
		int numIncorrectPriorTemp=0;
		
		int correctScore=0;
		int maxScore=ssl.get(0).score;
		
		int positionOfFirstCorrect=-1;
		int firstElementCorrect=0;
		int firstElementCorrectLoose=0;
		
		int sizeOfTopGroup=0;
		
		int numTopCorrect=0;
		
		for(int i=0; i<ssl.size(); i++){
			SiteScore ss=ssl.get(i);
			if(ss.score==ssl.get(0).score){sizeOfTopGroup++;}
			

//			boolean b=isCorrectHit(ss, original.chrom, original.strand, original.start, 1, thresh);
			boolean b=isCorrectHit(ss, ss.chrom, original.strand, original.start, original.stop, thresh);
			boolean b2=isCorrectHitLoose(ss, ss.chrom, original.strand, original.start, original.stop, thresh+120);
			if(b){
				if(i==0){firstElementCorrect=1;}
			}
			if(b2){
				if(i==0){firstElementCorrectLoose=1;}
				numCorrect++;
				numIncorrectPrior2+=numIncorrectPriorTemp;
				numIncorrectPriorTemp=0;
				
				if(numIncorrect==0){numTopCorrect++;}
				if(positionOfFirstCorrect<0){
					positionOfFirstCorrect=i;
					correctScore=ss.score;
				}
			}else{
				numIncorrect++;
				numIncorrectPriorTemp++;
				if(numCorrect<EXPECTED_SITES){
					numIncorrectPrior++;
				}
			}
		}
		
		assert(sizeOfTopGroup>0 && sizeOfTopGroup<=ssl.size());
		return new int[] {numCorrect, numIncorrect, numIncorrectPrior2,
				ssl.size(), correctScore, maxScore, firstElementCorrect, firstElementCorrectLoose, positionOfFirstCorrect, sizeOfTopGroup, numTopCorrect};
	}
	
	
}
