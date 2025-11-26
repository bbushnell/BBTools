package align2;

import java.util.Arrays;

import dna.ChromosomeArray;
import dna.Data;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.SiteScore;

/**
 * @author Brian Bushnell
 * @date Jun 20, 2013
 *
 */
public abstract class MSA {
	
	/**
	 * Creates appropriate MSA implementation based on specified class name.
	 * Factory method supporting 7 specialized alignment variants with different performance characteristics.
	 * Automatically selects JNI acceleration when available and requested.
	 *
	 * @param maxRows_ Maximum number of alignment matrix rows
	 * @param maxColumns_ Maximum number of alignment matrix columns
	 * @param classname Name of MSA implementation class to instantiate
	 * @return MSA instance optimized for specified requirements
	 */
	public static final MSA makeMSA(int maxRows_, int maxColumns_, String classname){
		flatMode=false;
		if("MultiStateAligner9ts".equalsIgnoreCase(classname)){
			return new MultiStateAligner9ts(maxRows_, maxColumns_);
		}else if("MultiStateAligner10ts".equalsIgnoreCase(classname)){
			return new MultiStateAligner10ts(maxRows_, maxColumns_);
		}else if("MultiStateAligner11ts".equalsIgnoreCase(classname)){
			if(Shared.USE_JNI){
				return new MultiStateAligner11tsJNI(maxRows_, maxColumns_);
			}else{
				return new MultiStateAligner11ts(maxRows_, maxColumns_);
			}
		}else if("MultiStateAligner11tsJNI".equalsIgnoreCase(classname)){
			return new MultiStateAligner11tsJNI(maxRows_, maxColumns_);
		}else if("MultiStateAligner9PacBio".equalsIgnoreCase(classname)){
			return new MultiStateAligner9PacBio(maxRows_, maxColumns_);
		}else if("MultiStateAligner9Flat".equalsIgnoreCase(classname)){
			return new MultiStateAligner9Flat(maxRows_, maxColumns_);
		}else if("MultiStateAligner9XFlat".equalsIgnoreCase(classname)){
			flatMode=true;
			return new MultiStateAligner9XFlat(maxRows_, maxColumns_);
		}else{
			assert(false) : "Unhandled MSA type: "+classname;
			return new MultiStateAligner11ts(maxRows_, maxColumns_);
		}
	}
	
	/**
	 * Constructs MSA with specified matrix dimensions.
	 * @param maxRows_ Maximum number of alignment matrix rows
	 * @param maxColumns_ Maximum number of alignment matrix columns
	 */
	public MSA(int maxRows_, int maxColumns_){
		maxRows=maxRows_;
		maxColumns=maxColumns_;
	}
	
	/** return new int[] {rows, maxC, maxS, max};
	 * Will not fill areas that cannot match minScore */
	public abstract int[] fillLimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore, int[] gaps);
	
	
	/** return new int[] {rows, maxC, maxS, max};
	 * Will not fill areas that cannot match minScore */
	public abstract int[] fillUnlimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int[] gaps);
	
	/**
	 * Fills alignment matrix using quality scores for enhanced alignment precision.
	 * Quality-aware alignment incorporating base quality information into scoring.
	 *
	 * @param read Query sequence to align
	 * @param ref Reference sequence for alignment
	 * @param baseScores Quality scores for read bases
	 * @param refStartLoc Starting position in reference sequence
	 * @param refEndLoc Ending position in reference sequence
	 * @return Array containing [rows, maxC, maxS, max] with optimal alignment coordinates
	 */
	@Deprecated
	/** return new int[] {rows, maxC, maxS, max}; */
	public abstract int[] fillQ(byte[] read, byte[] ref, byte[] baseScores, int refStartLoc, int refEndLoc);

	
	/** @return {score, bestRefStart, bestRefStop} */
	/** Generates the match string */
	public abstract byte[] traceback(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int row, int col, int state, boolean gapped);
	
	
	/** Generates the match string */
	public abstract byte[] traceback2(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int row, int col, int state);
	
	/** @return {score, bestRefStart, bestRefStop} */
	public abstract int[] score(final byte[] read, final byte[] ref, final int refStartLoc, final int refEndLoc,
			final int maxRow, final int maxCol, final int maxState, boolean gapped);
	
	/** @return {score, bestRefStart, bestRefStop}, or {score, bestRefStart, bestRefStop, padLeft, padRight} if more padding is needed */
	public abstract int[] score2(final byte[] read, final byte[] ref, final int refStartLoc, final int refEndLoc,
			final int maxRow, final int maxCol, final int maxState);
	
	
	/** Will not fill areas that cannot match minScore.
	 * @return {score, bestRefStart, bestRefStop}  */
	public final int[] fillAndScoreLimited(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, int minScore, int[] gaps){
		int a=Tools.max(0, refStartLoc);
		int b=Tools.min(ref.length-1, refEndLoc);
		assert(b>=a);
		
		int[] score;
		
		if(verbose && b-a<500){
			System.err.println(new String(read));
			System.err.println(new String(ref, a, b-a));
		}
		
		if(gaps==null){
			if(verbose){
				System.err.println("no gaps");
			}
			if(b-a>=maxColumns){
				System.err.println("Warning: Max alignment columns exceeded; restricting range. "+(b-a+1)+" > "+maxColumns);
				assert(false) : refStartLoc+", "+refEndLoc;
				b=Tools.min(ref.length-1, a+maxColumns-1);
			}
			int[] max=fillLimited(read, ref, a, b, minScore, gaps);
			score=(max==null ? null : score(read, ref, a, b, max[0], max[1], max[2], false));
		}else{
			if(verbose){System.err.println("\ngaps: "+Arrays.toString(gaps)+"\n"+new String(read)+"\ncoords: "+refStartLoc+", "+refEndLoc);}
			int[] max=fillLimited(read, ref, a, b, minScore, gaps);
			if(verbose){System.err.println("max: "+Arrays.toString(max));}
//			score=(max==null ? null : score(read, grefbuffer, 0, greflimit, max[0], max[1], max[2], true));
			score=(max==null ? null : score(read, ref, a, b, max[0], max[1], max[2], true));
		}
		return score;
	}
	
	/**
	 * Convenience method for fillAndScoreLimited using SiteScore object.
	 * Extracts alignment parameters from SiteScore and delegates to main implementation.
	 *
	 * @param read Query sequence to align
	 * @param ss SiteScore containing reference coordinates and gap information
	 * @param thresh Score threshold for alignment computation
	 * @param minScore Minimum score threshold for matrix computation
	 * @return Array containing [score, bestRefStart, bestRefStop] or null if alignment fails
	 */
	public final int[] fillAndScoreLimited(byte[] read, SiteScore ss, int thresh, int minScore){
		return fillAndScoreLimited(read, ss.chrom, ss.start, ss.stop, thresh, minScore, ss.gaps);
	}
	
//	public final int[] translateScoreFromGappedCoordinate(int[] score)
	
	/**
	 * Chromosome-based alignment with automatic reference array retrieval.
	 * Uses Data.getChromosome to access reference sequence and applies threshold padding.
	 *
	 * @param read Query sequence to align
	 * @param chrom Chromosome number for reference lookup
	 * @param start Starting coordinate on chromosome
	 * @param stop Ending coordinate on chromosome
	 * @param thresh Threshold value for coordinate padding
	 * @param minScore Minimum score threshold for matrix computation
	 * @param gaps Gap array for gapped alignment, null for ungapped
	 * @return Array containing [score, bestRefStart, bestRefStop] or null if alignment fails
	 */
	public final int[] fillAndScoreLimited(byte[] read, int chrom, int start, int stop, int thresh, int minScore, int[] gaps){
		return fillAndScoreLimited(read, Data.getChromosome(chrom).array, start-thresh, stop+thresh, minScore, gaps);
	}
	
	/**
	 * Quality-aware alignment with fillQ integration.
	 * Deprecated method that performed quality-score-based alignment.
	 *
	 * @param read Query sequence to align
	 * @param ref Reference sequence for alignment
	 * @param refStartLoc Starting position in reference sequence
	 * @param refEndLoc Ending position in reference sequence
	 * @param baseScores Quality scores for read bases
	 * @return Always returns null in current implementation
	 */
	@Deprecated
	public final int[] fillAndScoreQ(byte[] read, byte[] ref, int refStartLoc, int refEndLoc, byte[] baseScores){
		int a=Tools.max(0, refStartLoc);
		int b=Tools.min(ref.length-1, refEndLoc);
		assert(b>=a);
		if(b-a>=maxColumns){
			System.err.println("Warning: Max alignment columns exceeded; restricting range. "+(b-a+1)+" > "+maxColumns);
			b=Tools.min(ref.length-1, a+maxColumns-1);
		}
		int[] max=fillQ(read, ref, baseScores, a, b);
//		int[] score=score(read, ref,  a, b, max[0], max[1], max[2]);
//		return score;
		return null;
	}
	
	/**
	 * SiteScore-based quality alignment convenience method.
	 * Deprecated wrapper around fillAndScoreQ with SiteScore parameter extraction.
	 *
	 * @param read Query sequence to align
	 * @param ss SiteScore containing reference coordinates
	 * @param thresh Score threshold for alignment computation
	 * @param baseScores Quality scores for read bases
	 * @return Result from fillAndScoreQ implementation
	 */
	@Deprecated
	public final int[] fillAndScoreQ(byte[] read, SiteScore ss, int thresh, byte[] baseScores){
		return fillAndScoreQ(read, ss.chrom, ss.start, ss.stop, thresh, baseScores);
	}
	
	/**
	 * Chromosome-based quality alignment with automatic reference retrieval.
	 * Deprecated method combining chromosome lookup with quality-aware alignment.
	 *
	 * @param read Query sequence to align
	 * @param chrom Chromosome number for reference lookup
	 * @param start Starting coordinate on chromosome
	 * @param stop Ending coordinate on chromosome
	 * @param thresh Threshold value for coordinate padding
	 * @param baseScores Quality scores for read bases
	 * @return Result from fillAndScoreQ implementation
	 */
	@Deprecated
	public final int[] fillAndScoreQ(byte[] read, int chrom, int start, int stop, int thresh, byte[] baseScores){
		return fillAndScoreQ(read, Data.getChromosome(chrom).array, start-thresh, stop+thresh, baseScores);
	}
	
	/**
	 * Calculates alignment score without allowing insertions or deletions.
	 * Simple alignment scoring using direct base-to-base comparison only.
	 *
	 * @param read Query sequence to score
	 * @param ss SiteScore containing reference coordinates
	 * @return Alignment score for ungapped alignment
	 */
	public final int scoreNoIndels(byte[] read, SiteScore ss){
		ChromosomeArray cha=Data.getChromosome(ss.chrom);
		return scoreNoIndels(read, cha.array, ss.start);
	}

	/**
	 * Chromosome-based ungapped alignment scoring.
	 * Retrieves reference sequence from chromosome and performs ungapped scoring.
	 *
	 * @param read Query sequence to score
	 * @param chrom Chromosome number for reference lookup
	 * @param refStart Starting position on reference sequence
	 * @return Alignment score for ungapped alignment
	 */
	public final int scoreNoIndels(byte[] read, final int chrom, final int refStart){
		ChromosomeArray cha=Data.getChromosome(chrom);
		return scoreNoIndels(read, cha.array, refStart);
	}
	
	/**
	 * Quality-aware ungapped alignment scoring using SiteScore.
	 * Incorporates base quality scores into ungapped alignment evaluation.
	 *
	 * @param read Query sequence to score
	 * @param ss SiteScore containing reference coordinates
	 * @param baseScores Quality scores for read bases
	 * @return Quality-weighted alignment score for ungapped alignment
	 */
	public final int scoreNoIndels(byte[] read, SiteScore ss, byte[] baseScores){
		ChromosomeArray cha=Data.getChromosome(ss.chrom);
		return scoreNoIndels(read, cha.array, baseScores, ss.start);
	}

	/**
	 * Chromosome-based quality-aware ungapped scoring.
	 * Combines chromosome lookup with quality-weighted ungapped alignment.
	 *
	 * @param read Query sequence to score
	 * @param chrom Chromosome number for reference lookup
	 * @param refStart Starting position on reference sequence
	 * @param baseScores Quality scores for read bases
	 * @return Quality-weighted alignment score for ungapped alignment
	 */
	public final int scoreNoIndels(byte[] read, final int chrom, final int refStart, byte[] baseScores){
		ChromosomeArray cha=Data.getChromosome(chrom);
		return scoreNoIndels(read, cha.array, baseScores, refStart);
	}

//	public final int scoreNoIndels(byte[] read, final int chrom, final int refStart){
	
	/** Calculates score based on an array from Index. */
	public abstract int calcAffineScore(int[] locArray, byte[] baseScores, byte[] bases);

	/** Calculates score based on an array from Index using a kfilter.  Slightly slower. */
	public abstract int calcAffineScore(int[] locArray, byte[] baseScores, byte[] bases, int minContig);

	/**
	 * Abstract method for ungapped alignment scoring.
	 * Implementation-specific ungapped alignment between read and reference.
	 *
	 * @param read Query sequence to align
	 * @param ref Reference sequence for alignment
	 * @param refStart Starting position in reference sequence
	 * @return Ungapped alignment score
	 */
	public abstract int scoreNoIndels(byte[] read, byte[] ref, final int refStart);
	/**
	 * SiteScore-enhanced ungapped alignment scoring.
	 * Default implementation throws RuntimeException, requiring subclass implementation.
	 *
	 * @param read Query sequence to align
	 * @param ref Reference sequence for alignment
	 * @param refStart Starting position in reference sequence
	 * @param ss SiteScore object for enhanced scoring context
	 * @return Ungapped alignment score
	 * @throws RuntimeException If method not implemented in subclass
	 */
	public int scoreNoIndels(byte[] read, byte[] ref, final int refStart, final SiteScore ss){
		throw new RuntimeException("Unimplemented method in class "+this.getClass());
	}

	/**
	 * Generates match string for ungapped alignment.
	 * Creates alignment representation without insertion or deletion operations.
	 *
	 * @param read Query sequence for alignment
	 * @param ref Reference sequence for alignment
	 * @param refStart Starting position in reference sequence
	 * @return Match string representing ungapped alignment operations
	 */
	public abstract byte[] genMatchNoIndels(byte[] read, byte[] ref, final int refStart);

	/**
	 * Quality-aware abstract ungapped alignment scoring.
	 * Implementation-specific quality-weighted ungapped alignment.
	 *
	 * @param read Query sequence to align
	 * @param ref Reference sequence for alignment
	 * @param baseScores Quality scores for read bases
	 * @param refStart Starting position in reference sequence
	 * @return Quality-weighted ungapped alignment score
	 */
	public abstract int scoreNoIndels(byte[] read, byte[] ref, byte[] baseScores, final int refStart);
	/**
	 * Quality-aware alignment with fillQ integration.
	 * Deprecated method that performed quality-score-based alignment.
	 *
	 * @param read Query sequence to align
	 * @param ref Reference sequence for alignment
	 * @param refStartLoc Starting position in reference sequence
	 * @param refEndLoc Ending position in reference sequence
	 * @param baseScores Quality scores for read bases
	 * @return Always returns null in current implementation
	 */
	public int scoreNoIndels(byte[] read, byte[] ref, byte[] baseScores, final int refStart, SiteScore ss){
		throw new RuntimeException("Unimplemented method in class "+this.getClass());
	}
	
	/**
	 * Quality-aware ungapped scoring with simultaneous match string generation.
	 * Combines scoring and match string creation for efficiency in quality-aware mode.
	 *
	 * @param read Query sequence to align
	 * @param ref Reference sequence for alignment
	 * @param baseScores Quality scores for read bases
	 * @param refStart Starting position in reference sequence
	 * @param matchReturn Array for returning generated match string
	 * @return Quality-weighted ungapped alignment score
	 */
	public abstract int scoreNoIndelsAndMakeMatchString(byte[] read, byte[] ref, byte[] baseScores, final int refStart, byte[][] matchReturn);
	
	/**
	 * Ungapped scoring with simultaneous match string generation.
	 * Efficient combined scoring and match string creation without quality information.
	 *
	 * @param read Query sequence to align
	 * @param ref Reference sequence for alignment
	 * @param refStart Starting position in reference sequence
	 * @param matchReturn Array for returning generated match string
	 * @return Ungapped alignment score
	 */
	public abstract int scoreNoIndelsAndMakeMatchString(byte[] read, byte[] ref, final int refStart, byte[][] matchReturn);
	
	/** Assumes match string is in long format */
	public final boolean toLocalAlignment(Read r, SiteScore ss, byte[] basesM, int minToClip, float matchPointsMult){
		final byte[] match=r.match, bases=(r.strand()==Shared.PLUS ? r.bases : basesM);
		if(match==null || match.length<1){return false;}
		
		assert(match==ss.match);
		assert(match==r.match);
		assert(r.start==ss.start);
		assert(r.stop==ss.stop);
		
		if(r.containsXY2()){
			if(verbose){System.err.println("\nInitial0:");}
			if(verbose){System.err.println("0: match="+new String(match));}
			if(verbose){System.err.println("0: r.start="+r.start+", r.stop="+r.stop+"; len="+bases.length+"; reflen="+(r.stop-r.start+1));}
			ss.fixXY(bases, false, this);
			r.start=ss.start;
			r.stop=ss.stop;
			if(verbose){System.err.println("\nAfter fixXY:");}
			if(verbose){System.err.println("0: match="+new String(match));}
			if(verbose){System.err.println("0: r.start="+r.start+", r.stop="+r.stop+"; len="+bases.length+"; reflen="+(r.stop-r.start+1));}
			assert(match==ss.match);
			assert(match==r.match);
			assert(r.start==ss.start);
			assert(r.stop==ss.stop);
			assert(ss.lengthsAgree()) : ss.mappedLength()+"!="+ss.matchLength()+"\n"+ss+"\n\n"+r+"\n";
		}
		assert(ss.lengthsAgree()) : ss.mappedLength()+"!="+ss.matchLength()+"\n"+ss+"\n\n"+r+"\n";
		
		int maxScore=-1;
		
		int startLocC=-1;
		int stopLocC=-1;
		int lastZeroC=0;
		
		int startLocM=-1;
		int stopLocM=-1;
		int lastZeroM=0;
		
		int startLocR=-1;
		int stopLocR=-1;
		int lastZeroR=0;
		
		byte mode=match[0], prevMode='0';
		int current=0, prevStreak=0;
		int cpos=0;
		int rpos=r.start;
		int score=0;

		if(verbose){System.err.println("\nInitial:");}
		if(verbose){System.err.println("A: r.start="+r.start+", r.stop="+r.stop+"; rpos="+rpos+"; len="+bases.length+"; reflen="+(r.stop-r.start+1));}
		if(verbose){System.err.println("A: match=\n"+new String(match));}
		if(verbose){System.err.println(new String(bases));}
		if(verbose){System.err.println(Data.getChromosome(r.chrom).getString(r.start, Tools.max(r.stop, r.start+bases.length-1)));}
		
		if(verbose){
			int calcscore=score(match);
			System.err.println("A: score="+r.mapScore+", ss.slowScore="+ss.slowScore+", calcscore="+calcscore);
//			assert(ss.slowScore<=calcscore); //May be lower due to ambig3.  I found a case where this line fails, possibly due to long deletions?
		}
		
		for(int mpos=0; mpos<match.length; mpos++){
			byte c=match[mpos];
			
			if(mode==c){
				current++;
			}else{
				if(mode=='m'){
					if(score<=0){
						score=0;
						lastZeroC=cpos;
						lastZeroM=mpos-current;
						lastZeroR=rpos;
					}
					int add=calcMatchScore(current);
					score+=(matchPointsMult*add);
//					if(prevMode=='N' || prevMode=='R'){score=score+POINTS_MATCH2()-POINTS_MATCH();} //Don't penalize first match after N
					cpos+=current;
					rpos+=current;
					if(score>maxScore){
						maxScore=score;
						startLocC=lastZeroC;
						startLocM=lastZeroM;
						startLocR=lastZeroR;
						stopLocC=cpos-1;
						stopLocM=mpos-1;
						stopLocR=rpos-1;
					}
				}else if(mode=='S'){
					score+=calcSubScore(current);
					if(prevMode=='N' || prevMode=='R'){score=score+POINTS_SUB2()-POINTS_SUB();} //Don't penalize first sub after N
					else if(prevMode=='m' && prevStreak<2){score=score+POINTS_SUBR()-POINTS_SUB();}
					cpos+=current;
					rpos+=current;
				}else if(mode=='D'){
					score+=calcDelScore(current, true);
					rpos+=current;
				}else if(mode=='I'){
					score+=calcInsScore(current);
					cpos+=current;
				}else if(mode=='C'){
					cpos+=current;
					rpos+=current;
				}else if(mode=='X' || mode=='Y'){
					score+=calcInsScore(current);//TODO: Consider changing XY to subs
					cpos+=current;
					rpos+=current;
				}else if(mode=='N'){
					score+=calcNocallScore(current);
					cpos+=current;
					rpos+=current;
				}else if(mode=='R'){
					score+=calcNorefScore(current);
					cpos+=current;
					rpos+=current;
				}else{
					assert(false) : "Unhandled symbol "+mode+"\n"+(char)mode+"\n"+new String(match)+"\n"+new String(bases);
				}
				if(verbose){System.err.println("mode "+(char)mode+"->"+(char)c+"; rpos="+rpos);}
				prevMode=mode;
				prevStreak=current;
				mode=c;
				current=1;
			}
		}
		if(current>0){
			assert(mode==match[match.length-1]);
			if(mode=='m'){
				if(score<=0){
					score=0;
					lastZeroC=cpos;
					lastZeroM=match.length-current;
					lastZeroR=rpos;
				}
				int add=calcMatchScore(current);
				score+=(matchPointsMult*add);
//				if(prevMode=='N' || prevMode=='R'){score=score+POINTS_MATCH2()-POINTS_MATCH();} //Don't penalize first match after N
				cpos+=current;
				rpos+=current;
				if(score>maxScore){
					maxScore=score;
					startLocC=lastZeroC;
					startLocM=lastZeroM;
					startLocR=lastZeroR;
					stopLocC=cpos-1;
					stopLocM=match.length-1;
					stopLocR=rpos-1;
				}
			}else if(mode=='S'){
				score+=calcSubScore(current);
				if(prevMode=='N' || prevMode=='R'){score=score+POINTS_SUB2()-POINTS_SUB();} //Don't penalize first sub after N
				else if(prevMode=='m' && prevStreak<2){score=score+POINTS_SUBR()-POINTS_SUB();}
				cpos+=current;
				rpos+=current;
			}else if(mode=='D'){
				score+=calcDelScore(current, true);
				rpos+=current;
			}else if(mode=='I'){
				score+=calcInsScore(current);
				cpos+=current;
			}else if(mode=='C'){
				cpos+=current;
				rpos+=current;
			}else if(mode=='X' || mode=='Y'){
				score+=calcInsScore(current);
				cpos+=current;
				rpos+=current;
			}else if(mode=='N'){
				score+=calcNocallScore(current);
				cpos+=current;
				rpos+=current;
			}else if(mode=='R'){
				score+=calcNorefScore(current);
				cpos+=current;
				rpos+=current;
			}else if(mode!=0){
				assert(false) : "Unhandled symbol "+mode+"\n"+(char)mode+"\n"+new String(match)+"\n"+new String(bases);
			}
			if(verbose){System.err.println("mode "+(char)mode+"->end; rpos="+rpos);}
		}
		
		if(startLocC<0 || stopLocC<0){
			//This can happen if there are zero matches.  Which would be rare, but I have seen it occur.
			r.clearMapping();
//			assert(false) : "Failed: "+startLocC+", "+stopLocC+"\n"+r+"\n"+r.mate+"\n"+r.toFastq()+"\n"+(r.mate==null ? "null" : r.mate.toFastq());
			return false;
		}
		
		
		if(verbose){System.err.println("A: r.start="+r.start+", r.stop="+r.stop+"; rpos="+rpos+"; len="+bases.length+"; reflen="+(r.stop-r.start+1));}
		
		assert(rpos==r.stop+1) : "\n\n\n"+rpos+"!="+(r.stop+1)+"\n"+r+"\n\n"+
			(r.topSite()==null ? "null" : r.topSite().mappedLength()+", "+r.topSite().matchLength()+", "+r.topSite().start+", "+r.topSite().stop+"\n"+r.topSite());
		
		if(verbose){System.err.println("B: rpos="+rpos+", startLocR="+startLocR+", stopLocR="+stopLocR);}
		
		int headTrimR=startLocC;
		int headTrimM=startLocM;
		int tailTrimR=bases.length-stopLocC-1;
		int tailTrimM=match.length-stopLocM-1;
		
		if(verbose){System.err.println("C: headTrimR="+headTrimR+", headTrimM="+headTrimM+", tailTrimR="+tailTrimR+", tailTrimM="+tailTrimM);}
		
		if(headTrimR<=minToClip && headTrimM<=minToClip){
			headTrimR=headTrimM=0;
		}
		if(tailTrimR<=minToClip && tailTrimM<=minToClip){
			tailTrimR=tailTrimM=0;
		}
		if(headTrimR==0 && headTrimM==0 && tailTrimR==0 && tailTrimM==0){
			return false;
		}
		//Do trimming
		final int headDelta=headTrimR-headTrimM;
		final int tailDelta=tailTrimR-tailTrimM;
		final byte[] match2;
		
		if(verbose){System.err.println("D: headTrimR="+headTrimR+", headTrimM="+headTrimM+", tailTrimR="+tailTrimR+", tailTrimM="+tailTrimM);}
		if(verbose){System.err.println("D: headDelta="+headDelta+", tailDelta="+tailDelta);}
		
		if(headDelta==0 && tailDelta==0){
			//Length-neutral trimming
			match2=match;
			for(int i=0; i<headTrimM; i++){match[i]='C';}
			for(int i=match.length-tailTrimM; i<match.length; i++){match[i]='C';}
		}else{
			final int newlen=match.length-headTrimM-tailTrimM+headTrimR+tailTrimR;
			match2=new byte[newlen];
			for(int i=0; i<headTrimR; i++){match2[i]='C';}
			for(int i=match2.length-tailTrimR; i<match2.length; i++){match2[i]='C';}
			for(int i=headTrimM, i2=headTrimR, lim=match2.length-tailTrimR; i2<lim; i++, i2++){
				match2[i2]=match[i];
			}
		}
		
		assert(ss==null || ((ss.start==r.start) && (ss.stop==r.stop) && (ss.strand==r.strand()) && (ss.chrom==r.chrom) && (ss.match==r.match))) :
			"\nr="+r+"\nr2="+r.mate+"\nss=\n"+ss+"\n"+(ss==null ? "ss is null" : ((ss.start==r.start)+", "+(ss.stop==r.stop)+", "+
			(ss.strand==r.strand())+", "+(ss.chrom==r.chrom)+", "+(ss.match==r.match)));
		
		if(headTrimR!=0){r.start=startLocR-headTrimR;}
		if(tailTrimR!=0){r.stop=stopLocR+tailTrimR;}
		r.match=match2;
		
		if(matchPointsMult!=1f){
			maxScore=score(match);
		}
		if(ss!=null){maxScore=Tools.max(maxScore, ss.slowScore);}
		r.mapScore=maxScore;

		if(verbose){System.err.println("E: r.start="+r.start+", r.stop="+r.stop);}
		
		if(ss!=null){
			assert(maxScore>=ss.slowScore) : maxScore+", "+ss.slowScore+"\n"+r.toFastq();
			ss.match=r.match;
			ss.setLimits(r.start, r.stop);
			int pairedScore=ss.pairedScore>0 ? Tools.max(ss.pairedScore+(maxScore-ss.slowScore), 0) : 0;
		}
		
		if(!ss.perfect && ss.isPerfect(bases)){
			ss.perfect=ss.semiperfect=true;
			r.setPerfect(true);
			Arrays.fill(r.match, (byte)'m');
			ss.setSlowScore(maxScore);
		}else if(!ss.semiperfect && ss.isSemiPerfect(bases)){
			ss.semiperfect=true;
			ChromosomeArray cha=Data.getChromosome(ss.chrom);
			r.match=ss.match=genMatchNoIndels(bases, cha.array, ss.start);
			return toLocalAlignment(r, ss, basesM, minToClip, matchPointsMult);
		}
		return true;
	}
	

	/** Works in short or long format. */
	public final int score(byte[] match){
		if(match==null || match.length<1){return 0;}
		
		byte mode=match[0], prevMode='0';
		int current=0, prevStreak=0;
		int score=0;
		boolean hasDigit=false;
		
		for(int mpos=0; mpos<match.length; mpos++){
			byte c=match[mpos];
			
			if(mode==c){
				current++;
			}else if(Tools.isDigit(c)){
				current=(hasDigit ? current : 0)*10+c-'0';
				hasDigit=true;
			}else{
				if(mode=='m'){
					score+=calcMatchScore(current);
//					if(prevMode=='N' || prevMode=='R'){score=score+POINTS_MATCH2()-POINTS_MATCH();} //Don't penalize first match after N
				}else if(mode=='S'){
					score+=calcSubScore(current);
					if(prevMode=='N' || prevMode=='R'){score=score+POINTS_SUB2()-POINTS_SUB();} //Don't penalize first sub after N
					else if(prevMode=='m' && prevStreak<2){score=score+POINTS_SUBR()-POINTS_SUB();}
				}else if(mode=='D'){
					score+=calcDelScore(current, true);
				}else if(mode=='I'){
					score+=calcInsScore(current);
				}else if(mode=='C'){
					//do nothing
				}else if(mode=='X' || mode=='Y'){
					score+=calcInsScore(current);
				}else if(mode=='N'){
					score+=calcNocallScore(current);
				}else if(mode=='R'){
					score+=calcNorefScore(current);
				}else{
					assert(false) : "Unhandled symbol "+mode+"\n"+(char)mode+"\n"+new String(match);
				}
				if(verbose){System.err.println("mode "+(char)mode+"->"+(char)c+"\tcurrent="+current+"\tscore="+score);}
				prevMode=mode;
				prevStreak=current;
				mode=c;
				current=1;
				hasDigit=false;
			}
		}
		if(current>0){
			assert(hasDigit || mode==match[match.length-1]);
			if(mode=='m'){
				score+=calcMatchScore(current);
//				if(prevMode=='N' || prevMode=='R'){score=score+POINTS_MATCH2()-POINTS_MATCH();} //Don't penalize first match after N
			}else if(mode=='S'){
				score+=calcSubScore(current);
				if(prevMode=='N' || prevMode=='R'){score=score+POINTS_SUB2()-POINTS_SUB();} //Don't penalize first sub after N
				else if(prevMode=='m' && prevStreak<2){score=score+POINTS_SUBR()-POINTS_SUB();}
			}else if(mode=='D'){
				score+=calcDelScore(current, true);
			}else if(mode=='I'){
				score+=calcInsScore(current);
			}else if(mode=='C'){
				//do nothing
			}else if(mode=='X' || mode=='Y'){
				score+=calcInsScore(current);
			}else if(mode=='N'){
				score+=calcNocallScore(current);
			}else if(mode=='R'){
				score+=calcNorefScore(current);
			}else if(mode!=0){
				assert(false) : "Unhandled symbol "+mode+"\n"+(char)mode+"\n"+new String(match);
			}
			if(verbose){System.err.println("mode "+(char)mode+"->end; score="+score);}
		}
		
		return score;
	}
	
	/**
	 * Calculates maximum possible quality score for given number of bases.
	 * @param numBases Number of bases in sequence
	 * @return Maximum achievable quality score
	 */
	public abstract int maxQuality(int numBases);
	
	/**
	 * Calculates maximum quality score based on actual base quality values.
	 * @param baseScores Array of base quality scores
	 * @return Maximum achievable quality score for these specific base scores
	 */
	public abstract int maxQuality(byte[] baseScores);
	
	/**
	 * Calculates maximum score for imperfect alignment with given number of bases.
	 * @param numBases Number of bases in sequence
	 * @return Maximum achievable score allowing for imperfect matches
	 */
	public abstract int maxImperfectScore(int numBases);
	
	/**
	 * Calculates maximum imperfect score based on actual base quality values.
	 * @param baseScores Array of base quality scores
	 * @return Maximum achievable imperfect score for these specific base scores
	 */
	public abstract int maxImperfectScore(byte[] baseScores);
	
	/**
	 * Formats integer array for display with fixed-width columns.
	 * Creates aligned string representation for debugging and visualization.
	 * @param a Integer array to format
	 * @return Formatted string with fixed-width columns
	 */
	public static final String toString(int[] a){
		
		int width=7;
		
		StringBuilder sb=new StringBuilder((a.length+1)*width+2);
		for(int num : a){
			String s=" "+num;
			int spaces=width-s.length();
			assert(spaces>=0) : width+", "+s.length()+", "+s+", "+num+", "+spaces;
			for(int i=0; i<spaces; i++){sb.append(' ');}
			sb.append(s);
		}
		
		return sb.toString();
	}
	
	/**
	 * Prints all modes of packed alignment matrix for debugging.
	 * Displays both time and score information for complete matrix visualization.
	 *
	 * @param packed 3D packed matrix containing alignment information
	 * @param readlen Length of read sequence
	 * @param reflen Length of reference sequence
	 * @param TIMEMASK Bit mask for extracting time information
	 * @param SCOREOFFSET Bit offset for extracting score information
	 */
	static void printMatrix(int[][][] packed, int readlen, int reflen, int TIMEMASK, int SCOREOFFSET){
		for(int mode=0; mode<packed.length; mode++){
			printMatrix(packed, readlen, reflen, TIMEMASK, SCOREOFFSET, mode);
		}
	}
	
	/**
	 * Prints specific mode of packed alignment matrix for debugging.
	 * Displays both time and score information for specified matrix mode.
	 *
	 * @param packed 3D packed matrix containing alignment information
	 * @param readlen Length of read sequence
	 * @param reflen Length of reference sequence
	 * @param TIMEMASK Bit mask for extracting time information
	 * @param SCOREOFFSET Bit offset for extracting score information
	 * @param mode Specific matrix mode to display
	 */
	static void printMatrix(int[][][] packed, int readlen, int reflen, int TIMEMASK, int SCOREOFFSET, int mode){
		final int ylim=Tools.min(readlen+1, packed[mode].length);
		final int xlim=Tools.min(reflen+1, packed[mode].length);
		for(int row=0; row<ylim; row++){
			System.out.println(toScorePacked(packed[mode][row], SCOREOFFSET, xlim));
		}
		System.out.println();
		for(int row=0; row<ylim; row++){
			System.out.println(toTimePacked(packed[mode][row], TIMEMASK, xlim));
		}
		System.out.println();
	}
	
	/**
	 * Extracts and formats time information from packed matrix row.
	 * Creates fixed-width display of time values for matrix visualization.
	 *
	 * @param a Packed matrix row containing time and score information
	 * @param TIMEMASK Bit mask for extracting time values
	 * @param lim Maximum number of elements to process
	 * @return Formatted string showing time values with fixed-width columns
	 */
	public static final String toTimePacked(int[] a, int TIMEMASK, int lim){
		int width=6;
		lim=Tools.min(lim, a.length);
		
		StringBuilder sb=new StringBuilder((a.length+1)*width+2);
		for(int j=0; j<lim; j++){
			int num=a[j]&TIMEMASK;
			String s=" "+num;
			int spaces=width-s.length();
			assert(spaces>=0) : width+", "+s.length()+", "+s+", "+num+", "+spaces;
			for(int i=0; i<spaces; i++){sb.append(' ');}
			sb.append(s);
		}
		
		return sb.toString();
	}
	
	/**
	 * Extracts and formats score information from packed matrix row.
	 * Creates fixed-width display of score values with overflow protection for matrix visualization.
	 *
	 * @param a Packed matrix row containing time and score information
	 * @param SCOREOFFSET Bit offset for extracting score values
	 * @param lim Maximum number of elements to process
	 * @return Formatted string showing score values with fixed-width columns
	 */
	public static final String toScorePacked(int[] a, int SCOREOFFSET, int lim){
		int width=6;
		lim=Tools.min(lim, a.length);
		
//		String minString=" -";
//		String maxString="  ";
//		while(minString.length()<width){minString+='9';}
//		while(maxString.length()<width){maxString+='9';}

		String minString=" -";
		String maxString=" +";
		while(minString.length()<width){minString=minString+' ';}
		while(maxString.length()<width){maxString=maxString+' ';}
		
		StringBuilder sb=new StringBuilder((a.length+1)*width+2);
		for(int j=0; j<lim; j++){
			int num=a[j]>>SCOREOFFSET;
			String s=" "+num;
			if(s.length()>width){s=num>0 ? maxString : minString;}
			int spaces=width-s.length();
			assert(spaces>=0) : width+", "+s.length()+", "+s+", "+num+", "+spaces;
			for(int i=0; i<spaces; i++){sb.append(' ');}
			sb.append(s);
		}
		
		return sb.toString();
	}
	
	/**
	 * Formats byte array for display with fixed-width columns.
	 * Creates aligned string representation for debugging byte sequences.
	 * @param a Byte array to format
	 * @return Formatted string with fixed-width columns
	 */
	public static final String toString(byte[] a){
		
		int width=6;
		
		StringBuilder sb=new StringBuilder((a.length+1)*width+2);
		for(int num : a){
			String s=" "+num;
			int spaces=width-s.length();
			assert(spaces>=0);
			for(int i=0; i<spaces; i++){sb.append(' ');}
			sb.append(s);
		}
		
		return sb.toString();
	}
	
	/**
	 * Extracts substring from byte array and converts to String.
	 * Creates string representation of specified byte array region.
	 *
	 * @param ref Source byte array
	 * @param startLoc Starting position (inclusive)
	 * @param stopLoc Ending position (inclusive)
	 * @return String representation of specified byte array region
	 */
	public static final String toString(byte[] ref, int startLoc, int stopLoc){
		StringBuilder sb=new StringBuilder(stopLoc-startLoc+1);
		for(int i=startLoc; i<=stopLoc; i++){sb.append((char)ref[i]);}
		return sb.toString();
	}
	
	/**
	 * Calculates score for consecutive matching bases.
	 * Uses differential scoring with higher points for first match and lower points for additional matches.
	 * @param len Number of consecutive matching bases
	 * @return Total score for match run of specified length
	 */
	public final int calcMatchScore(int len){
		assert(len>0) : len;
		return POINTS_MATCH()+(len-1)*POINTS_MATCH2();
	}
	
	/**
	 * Calculates penalty score for consecutive substitutions.
	 * Uses tiered penalty system with increasing penalties for longer substitution runs.
	 * @param len Number of consecutive substitutions
	 * @return Total penalty score for substitution run of specified length
	 */
	public final int calcSubScore(int len){
		assert(len>0) : len;
		final int lim3=LIMIT_FOR_COST_3();
		int score=POINTS_SUB();
		if(len>lim3){
			score+=(len-lim3)*POINTS_SUB3();
			len=lim3;
		}
		if(len>1){
			score+=(len-1)*POINTS_SUB2();
		}
		return score;
	}
	
	/**
	 * Calculates score for bases aligned to no-reference regions.
	 * @param len Number of bases in no-reference region
	 * @return Score for no-reference alignment
	 */
	public final int calcNorefScore(int len){return len*POINTS_NOREF();}
	
	/**
	 * Calculates score for bases aligned to no-call reference regions.
	 * @param len Number of bases in no-call region
	 * @return Score for no-call alignment
	 */
	public final int calcNocallScore(int len){return len*POINTS_NOCALL();}
	
	/**
	 * Calculates penalty score for deletion operations.
	 * @param len Length of deletion
	 * @param approximateGaps Whether to use approximate gap penalty calculation
	 * @return Penalty score for deletion of specified length
	 */
	public abstract int calcDelScore(int len, boolean approximateGaps);
	
	/**
	 * Calculates penalty score for insertion operations.
	 * @param len Length of insertion
	 * @return Penalty score for insertion of specified length
	 */
	public abstract int calcInsScore(int len);
	
	/**
	 * Converts minimum identity percentage to minimum ratio for specified MSA class.
	 * Delegates to appropriate MSA implementation for class-specific ratio calculation.
	 *
	 * @param minid Minimum identity percentage (0.0 to 1.0)
	 * @param classname Name of MSA implementation class
	 * @return Minimum ratio value corresponding to specified identity threshold
	 */
	public static final float minIdToMinRatio(double minid, String classname){
		if("MultiStateAligner9ts".equalsIgnoreCase(classname)){
			return MultiStateAligner9ts.minIdToMinRatio(minid);
		}else if("MultiStateAligner10ts".equalsIgnoreCase(classname)){
			return MultiStateAligner10ts.minIdToMinRatio(minid);
		}else if("MultiStateAligner11ts".equalsIgnoreCase(classname)){
			return MultiStateAligner11ts.minIdToMinRatio(minid);
		}else if("MultiStateAligner9PacBio".equalsIgnoreCase(classname)){
			return MultiStateAligner9PacBio.minIdToMinRatio(minid);
		}else if("MultiStateAligner9Flat".equalsIgnoreCase(classname)){
			return MultiStateAligner9Flat.minIdToMinRatio(minid);
		}else if("MultiStateAligner9XFlat".equalsIgnoreCase(classname)){
			return MultiStateAligner9XFlat.minIdToMinRatio(minid);
		}else{
			assert(false) : "Unhandled MSA type: "+classname;
			return MultiStateAligner11ts.minIdToMinRatio(minid);
		}
	}
	
	/** Gap buffer size constant from Shared configuration */
	static final int GAPBUFFER=Shared.GAPBUFFER;
	/** Secondary gap buffer size constant from Shared configuration */
	static final int GAPBUFFER2=Shared.GAPBUFFER2;
	/** Gap length constant from Shared configuration */
	static final int GAPLEN=Shared.GAPLEN;
	/** Minimum gap size constant from Shared configuration */
	static final int MINGAP=Shared.MINGAP;
	/** Gap cost penalty constant from Shared configuration */
	static final int GAPCOST=Shared.GAPCOST;
	/** Gap character constant from Shared configuration */
	static final byte GAPC=Shared.GAPC;
	
	/** Seemingly to clear out prior data from the gref.  Not sure what else it's used for. */
	static final int GREFLIMIT2_CUSHION=128; //Tools.max(GAPBUFFER2, GAPLEN);
	
	
	/**DO NOT MODIFY*/
	public abstract byte[] getGrefbuffer();

//	public final int[] vertLimit;
//	public final int[] horizLimit;

	/**
	 * Returns string representation of vertical alignment limits.
	 * Debugging method for displaying vertical matrix boundaries.
	 * @return CharSequence showing vertical alignment limits
	 */
	public abstract CharSequence showVertLimit();
	/**
	 * Returns string representation of horizontal alignment limits.
	 * Debugging method for displaying horizontal matrix boundaries.
	 * @return CharSequence showing horizontal alignment limits
	 */
	public abstract CharSequence showHorizLimit();

////	public static final int MODEBITS=2;
//	public static final int TIMEBITS=11;
//	public static final int SCOREBITS=32-TIMEBITS;
//	public static final int MAX_TIME=((1<<TIMEBITS)-1);
//	public static final int MAX_SCORE=((1<<(SCOREBITS-1))-1)-2000;
//	public static final int MIN_SCORE=0-MAX_SCORE; //Keeps it 1 point above "BAD".
//
////	public static final int MODEOFFSET=0; //Always zero.
////	public static final int TIMEOFFSET=0;
	/** Returns bit offset for score information in packed matrix format.
	 * @return Bit offset used to extract score values from packed integers */
	public abstract int SCOREOFFSET();
//
////	public static final int MODEMASK=~((-1)<<MODEBITS);
////	public static final int TIMEMASK=(~((-1)<<TIMEBITS))<<TIMEOFFSET;
//	public static final int TIMEMASK=~((-1)<<TIMEBITS);
//	public static final int SCOREMASK=(~((-1)<<SCOREBITS))<<SCOREOFFSET;
	
	/** Match/substitution mode constant for alignment state */
	static final byte MODE_MS=0;
	/** Deletion mode constant for alignment state */
	static final byte MODE_DEL=1;
	/** Insertion mode constant for alignment state */
	static final byte MODE_INS=2;
	/** Substitution mode constant for alignment state */
	static final byte MODE_SUB=3;
	
	//These are to allow constants to be overridden
	/** Returns scoring points for alignment to no-reference regions.
	 * @return Point value for no-reference alignment */
	public abstract int POINTS_NOREF();
	/** Returns scoring points for alignment to no-call reference bases.
	 * @return Point value for no-call alignment */
	public abstract int POINTS_NOCALL();
	/** Returns scoring points for first matching base in a run.
	 * @return Point value for initial match */
	public abstract int POINTS_MATCH();
	/** Returns scoring points for additional matching bases in a run.
	 * @return Point value for subsequent matches */
	public abstract int POINTS_MATCH2();
	/** Returns scoring points for compatible base alignments.
	 * @return Point value for compatible alignment */
	public abstract int POINTS_COMPATIBLE();
	/** Returns penalty points for first substitution in a run.
	 * @return Penalty value for initial substitution */
	public abstract int POINTS_SUB();
	/** Returns penalty points for reduced substitution scoring.
	 * @return Reduced penalty value for specific substitution contexts */
	public abstract int POINTS_SUBR();
	/** Returns penalty points for second substitution in a run.
	 * @return Penalty value for subsequent substitutions */
	public abstract int POINTS_SUB2();
	/** Returns penalty points for third and later substitutions in a run.
	 * @return Penalty value for extended substitution runs */
	public abstract int POINTS_SUB3();
	/** Returns penalty points for match-substitution transitions.
	 * @return Penalty value for match to substitution transitions */
	public abstract int POINTS_MATCHSUB();
	/** Returns penalty points for first insertion in a gap.
	 * @return Penalty value for gap opening (insertion) */
	public abstract int POINTS_INS();
	/** Returns penalty points for second insertion in a gap.
	 * @return Penalty value for gap extension (insertion) */
	public abstract int POINTS_INS2();
	/** Returns penalty points for third insertion in a gap.
	 * @return Penalty value for extended insertion gaps */
	public abstract int POINTS_INS3();
	/** Returns penalty points for fourth and later insertions in a gap.
	 * @return Penalty value for long insertion gaps */
	public abstract int POINTS_INS4();
	/** Returns penalty points for first deletion in a gap.
	 * @return Penalty value for gap opening (deletion) */
	public abstract int POINTS_DEL();
	/** Returns penalty points for second deletion in a gap.
	 * @return Penalty value for gap extension (deletion) */
	public abstract int POINTS_DEL2();
	/** Returns penalty points for third deletion in a gap.
	 * @return Penalty value for extended deletion gaps */
	public abstract int POINTS_DEL3();
	/** Returns penalty points for fourth deletion in a gap.
	 * @return Penalty value for long deletion gaps (fourth position) */
	public abstract int POINTS_DEL4();
	/** Returns penalty points for fifth and later deletions in a gap.
	 * @return Penalty value for very long deletion gaps */
	public abstract int POINTS_DEL5();
	/** Returns penalty points for deletions aligned to N bases in reference.
	 * @return Penalty value for deletions at ambiguous reference positions */
	public abstract int POINTS_DEL_REF_N();
	/** Returns penalty points for gap operations.
	 * @return General penalty value for gap alignment */
	public abstract int POINTS_GAP();

	/** Returns time slip value for matrix computation.
	 * @return Time slip parameter for alignment matrix calculations */
	public abstract int TIMESLIP();
	/** Returns bit mask for 5-bit operations.
	 * @return 5-bit mask value for bitwise operations */
	public abstract int MASK5();
	
	/** Returns barrier value for insertion state transitions.
	 * @return Barrier value for insertion alignment state */
	abstract int BARRIER_I1();
	/** Returns barrier value for deletion state transitions.
	 * @return Barrier value for deletion alignment state */
	abstract int BARRIER_D1();

	/** Returns length limit for applying third-tier penalty costs.
	 * @return Length threshold for tier 3 penalty application */
	public abstract int LIMIT_FOR_COST_3();
	/** Returns length limit for applying fourth-tier penalty costs.
	 * @return Length threshold for tier 4 penalty application */
	public abstract int LIMIT_FOR_COST_4();
	/** Returns length limit for applying fifth-tier penalty costs.
	 * @return Length threshold for tier 5 penalty application */
	public abstract int LIMIT_FOR_COST_5();
	
	/** Returns value representing invalid or bad alignment score.
	 * @return Sentinel value for invalid alignment results */
	public abstract int BAD();
	
	
	/** Maximum number of rows in alignment matrix */
	public final int maxRows;
	/** Maximum number of columns in alignment matrix */
	public final int maxColumns;

	/** Count of limited alignment iterations performed */
	public long iterationsLimited=0;
	/** Count of unlimited alignment iterations performed */
	public long iterationsUnlimited=0;

	/** Enable verbose debugging output for alignment operations */
	public boolean verbose=false;
	/**
	 * Enable additional verbose debugging output for detailed alignment analysis
	 */
	public boolean verbose2=false;

	/** Fixed bandwidth constraint for banded alignment algorithms */
	public static int bandwidth=0;
	/** Proportional bandwidth constraint as ratio of sequence length */
	public static float bandwidthRatio=0;
	/** Enable flat memory layout mode for cache-optimized alignment */
	public static boolean flatMode=false;
	
	/** Minimum score adjustment constant for alignment thresholding */
	public static final int MIN_SCORE_ADJUST=120;

}
