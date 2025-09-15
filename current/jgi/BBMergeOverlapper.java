package jgi;

import java.util.Arrays;

import assemble.Tadpole;
import assemble.Tadpole1;
import assemble.Tadpole2;
import dna.AminoAcid;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import structures.FloatList;
import ukmer.Kmer;

/**
 * @author Brian Bushnell
 * @date Apr 15, 2014
 *
 */
public final class BBMergeOverlapper {

	static{
		if(false && Shared.USE_JNI){//TODO: Disabled!
			Shared.loadJNI();
		}
	}
	
	private static native final int mateByOverlapJNI(byte[] a_bases, byte[] b_bases, byte[] a_quality, byte[] b_quality,
			float[] aprob, float[] bprob, int[] rvector, int minOverlap0, int minOverlap, int minInsert0, int margin,
			int maxMismatches0, int maxMismatches, int minq);

	private static native final int mateByOverlapRatioJNI_WithQualities(byte[] a_bases, byte[] b_bases, byte[] a_quality, byte[] b_quality,
			float[] aprob, float[] bprob, int[] rvector, int minOverlap0, int minOverlap, int minInsert0, int minInsert, float maxRatio,
			float margin, float offset);

	private static native final int mateByOverlapRatioJNI(byte[] a_bases, byte[] b_bases,
			int[] rvector, int minOverlap0, int minOverlap, int minInsert0, int minInsert, float maxRatio,
			float margin, float offset, float gIncr, float bIncr);
	
	protected static final int mateByOverlap(Read a, Read b, float[] aprob, float[] bprob, int[] rvector,
			int minOverlap0, final int minOverlap, final int minInsert0, int margin, final int maxMismatches0, final int maxMismatches, final int minq) {
		if(rvector==null){rvector=KillSwitch.allocInt1D(5);}
		final int x;
		if(false && Shared.USE_JNI){//TODO: Disabled!
			x=mateByOverlapJNI(a.bases, b.bases, a.quality, b.quality, aprob, bprob, rvector, minOverlap0, minOverlap, minInsert0, margin, maxMismatches0, maxMismatches, minq);
		}else if(!TAG_CUSTOM){
			x=mateByOverlapJava_unrolled(a, b, aprob, bprob, rvector, minOverlap0, minOverlap, minInsert0, margin, maxMismatches0, maxMismatches, minq);
		}else{
			x=mateByOverlapJava_unrolled_training(a, b, aprob, bprob, rvector, minOverlap0, minOverlap, minInsert0, margin, maxMismatches0+40, maxMismatches+40, minq);
		}
		return x;
	}
	
	public static final int mateByOverlapRatio(Read a, Read b, float[] aprob, float[] bprob, int[] rvector,
			int minOverlap0, int minOverlap, int minInsert0, int minInsert, final float maxRatio, final float minSecondRatio,
			final float margin, final float offset, final float gIncr, final float bIncr, boolean useQuality) {
		if(rvector==null){rvector=KillSwitch.allocInt1D(5);}
//		final boolean swapped;
//		if(a.length()>b.length()){
//			swapped=true;
//			a.swapBasesWithMate();
//			a.reverseComplement();
//			b.reverseComplement();
//		}else{
//			swapped=false;
//		}
//		assert(false) : useQuality+", "+a.quality;
		
		final int x;
		if(false && Shared.USE_JNI/* && !useQuality*/){
			if(useQuality && a.quality!=null && b.quality!=null){
				x=mateByOverlapRatioJNI_WithQualities(a.bases, b.bases, a.quality, b.quality, aprob, bprob, rvector, minOverlap0, minOverlap, minInsert0, minInsert, maxRatio, margin, offset);
			}else{
				x=mateByOverlapRatioJNI(a.bases, b.bases, rvector, minOverlap0, minOverlap, minInsert0, minInsert, maxRatio, margin, offset, gIncr, bIncr);
			}
		}else{
			if(useQuality && a.quality!=null && b.quality!=null){
				x=mateByOverlapRatioJava_WithQualities(a, b, aprob, bprob, rvector, minOverlap0, minOverlap, minInsert0, minInsert,
						maxRatio, minSecondRatio, margin, offset);
			}else{
				x=mateByOverlapRatioJava(a, b, rvector, minOverlap0, minOverlap, minInsert0, minInsert, maxRatio, minSecondRatio, margin, offset, gIncr, bIncr);
			}
		}
//		if(swapped){
//			a.swapBasesWithMate();
//			a.reverseComplement();
//			b.reverseComplement();
//		}
//		assert(x>-1) : x+", "+a.quality;
//		System.err.println("returning x="+x);
		return x;
	}
	
	protected static final int mateByOverlapRatioJava_WithQualities(Read a, Read b, float[] aprob, float[] bprob, int[] rvector,
			int minOverlap0, int minOverlap, int minInsert0, int minInsert, float maxRatio, final float minSecondRatio,
			final float margin, final float offset) {
		assert(rvector!=null);
		assert(margin>=1);
		minOverlap=Tools.max(4, minOverlap0, minOverlap);
		minOverlap0=Tools.mid(4, minOverlap0, minOverlap);
		if(rvector==null){rvector=KillSwitch.allocInt1D(5);}
		
		final byte[] abases=a.bases, bbases=b.bases, aqual=a.quality, bqual=b.quality;
		final int alen=abases.length, blen=bbases.length;
		final int minLength=Tools.min(alen, blen);
		
		assert(aqual!=null && bqual!=null);
		{//Set the correctness probabilities from the quality scores.
			for(int i=0; i<aqual.length; i++){aprob[i]=probCorrect3[aqual[i]];}
			for(int i=0; i<bqual.length; i++){bprob[i]=probCorrect3[bqual[i]];}
		}
		
		{
			final float x;
//			if(!BBMerge.TRAINING) {
				x=findBestRatio_WithQualities(a, b, aprob, bprob, minOverlap0, minOverlap, minInsert, maxRatio, offset);
//			}else{
//				x=findBestRatio_WithQualitiesTraining(a, b, aprob, bprob, minOverlap0, minOverlap, minInsert, maxRatio, offset);
//			}
			if(!TAG_CUSTOM && x>maxRatio){
				rvector[2]=minLength;
				rvector[4]=0;
				return -1;
			}
			maxRatio=Tools.min(maxRatio, x);
		}
		
//		final float altBadlimit=Tools.max(maxRatio, 0.07f)*2f*alen+2;
		final float margin2=(margin+offset)/minLength;

		int bestInsert=-1;
		int bestOverlap=-1;
		int bestBadInt=-1;
		float bestBad=minLength;
		float bestGood=0;
		float bestRatio=1;
		boolean ambig=false;

		float secondBestBad=0;
		float secondBestGood=0;
		float secondBestRatio=1;
		int secondBestInsert=0;
		int secondBestOverlap=0;
		int secondBestBadInt=-1;
		
		final float extraMult=(TAG_CUSTOM ? 80 : 1.2f);
		
		final int largestInsertToTest=(alen+blen-minOverlap0);
		final int smallestInsertToTest=minInsert0;
		for(int insert=largestInsertToTest; insert>=smallestInsertToTest; insert--){
			if(verbose){System.err.println("d\nTesting read "+a.numericID+", overlap "+insert+", insert "+(alen+blen-insert));}

			float good=0, bad=0;
			int badInt=0;
			
			final int istart=(insert<=blen ? 0 : insert-blen);
			final int jstart=(insert>=blen ? 0 : blen-insert);
			
			final int overlapLength=Tools.min(alen-istart, blen-jstart, insert);
//			final float badlimit=extraMult*Tools.max(minSecondRatio, Tools.min(altBadlimit, (Tools.min(bestRatio, maxRatio)))*margin*overlapLength)+1f;
//			final float badlimit=extraMult*Tools.max(minSecondRatio, Tools.min(bestRatio, maxRatio)*margin*overlapLength)+1f;
			final float badlimit=extraMult*(Tools.min(bestRatio, maxRatio)*margin*overlapLength)+1f+extraBadlimit;
//			final float badlimit=extraMult*Tools.max(minSecondRatio, Tools.min(bestRatio, maxRatio))*margin*overlapLength+1f;
			
			final int imax=istart+overlapLength;
			for(int i=istart, j=jstart; i<imax && bad<=badlimit; i++, j++){
				assert(i>=0 && i<alen && j>=0 && j<blen) : "\njstart="+jstart+", j="+j+
				", istart="+istart+", i="+i+" \n"+"insert="+insert+", overlap="+overlapLength+", a.length="+a.length()+
				", b.length="+b.length()+", bad="+bad+", badlimit="+badlimit+", good="+good;
				final byte ca=abases[i], cb=bbases[j];
				
				final float x=aprob[i]*bprob[j];

				if(ca==cb){good+=x;}
				else{
					bad+=x;
					badInt++;
				}
			}
			
//			if(verbose || true){
//				System.err.println("istart="+istart+", jstart="+jstart+", overlapLength="+overlapLength+", overlap="+overlap+", bestOverlap="+bestOverlap);
//				System.err.println("overlap="+overlap+", bad="+bad+", good="+good);
//				System.err.println("bestGood="+bestGood+", bestBad="+bestBad);
//				System.err.println();
//			}

			if(bad<=badlimit){
				if(!TAG_CUSTOM && bad==0 && good>minOverlap0 && good<minOverlap){
					rvector[2]=bestBadInt;//=(int)(bestBad+0.95f);
					rvector[4]=1;
					return -1;
				}

				float ratio=(bad+offset)/overlapLength;
//				System.err.println("*** ratio="+ratio+", bestRatio="+bestRatio);

				if(ratio<bestRatio*margin){

					ambig=(ratio*margin>=bestRatio || good<minOverlap);
					if(ratio<bestRatio){
						secondBestInsert=bestInsert;
						secondBestOverlap=bestOverlap;
						secondBestBad=bestBad;
						secondBestGood=bestGood;
						secondBestRatio=bestRatio;
						secondBestBadInt=bestBadInt;
						
						bestInsert=insert;
						bestOverlap=overlapLength;
						bestBad=bad;
						bestGood=good;
						bestRatio=ratio;
						bestBadInt=badInt;
					}
					else if(ratio<secondBestRatio){
						secondBestInsert=insert;
						secondBestOverlap=overlapLength;
						secondBestBad=bad;
						secondBestGood=good;
						secondBestRatio=ratio;
						secondBestBadInt=badInt;
					}
					
					if(!TAG_CUSTOM && ((ambig && bestRatio<margin2) || secondBestRatio<minSecondRatio)){
						rvector[2]=bestBadInt;//=(int)(bestBad+0.95f);
						rvector[4]=1;
						return -1;
					}
				}
			}
		}
		if(secondBestRatio<minSecondRatio){ambig=true;}
		if(MAKE_VECTOR) {
//			System.err.println("Read "+a.numericID+": "+a.obj+", "+b.obj);
			FloatList fl=b.fetchVector();
			assert(fl.size==3) : fl.size;
			fl.add((a.length()-100)*0.01f);
			fl.add((b.length()-100)*0.01f);
			fl.add(bestInsert*.004f);
			fl.add(bestOverlap/(bestOverlap+50f));
			fl.add((bestBad+1f)/(bestBad+5f));
			fl.add((bestGood+1f)/(bestGood+5f));
			fl.add(bestRatio);
			fl.add((bestBadInt+1f)/(bestBadInt+5f));

			fl.add(secondBestInsert*.004f);
			fl.add(secondBestOverlap/(secondBestOverlap+50f));
			fl.add((secondBestBad+1f)/(secondBestBad+5f));
			fl.add((secondBestGood+1f)/(secondBestGood+5f));
			fl.add(secondBestRatio);
			fl.add((secondBestBadInt)/(secondBestBadInt+5f));
			
			fl.add((secondBestRatio+1)/(bestRatio+1));
			fl.add(secondBestBad/(bestBad+8));
			fl.add(secondBestGood/(bestGood+8));
			fl.add((bestOverlap+1f)/(secondBestOverlap+bestOverlap+1f));
		}
		
		if(TAG_CUSTOM){
			final ByteBuilder sb=a.fetchBB();
			if(secondBestRatio<minSecondRatio){ambig=true;}
////			StringBuilder sb=new StringBuilder(a.id);
//
//			sb.append("_alen=").append(a.length()*0.01f, 2);
//			sb.append("_blen=").append(b.length()*0.01f, 2);
//			sb.append("_bi=").append(bestInsert*.01f, 2);
//			sb.append("_bo=").append(bestOverlap*0.1f, 1);
//			sb.append("_bb=").append(Tools.format("%.4f", bestBad*.1f));
//			sb.append("_bg=").append(Tools.format("%.4f", bestGood*.1f));
//			sb.append("_br=").append(Tools.format("%.4f", bestRatio));
//			sb.append("_bbi=").append(bestBadInt*.1f, 1);
//
//			sb.append("_sbi=").append(secondBestInsert*.01f, 2);
//			sb.append("_sbo=").append(secondBestOverlap*0.1f, 1);
//			sb.append("_sbb=").append(Tools.format("%.4f", secondBestBad*.1f));
//			sb.append("_sbg=").append(Tools.format("%.4f", secondBestGood*.1f));
//			sb.append("_sbr=").append(Tools.format("%.4f", secondBestRatio));
//			sb.append("_sbbi=").append(secondBestBadInt*.1f, 1);
//			sb.append("_rr=").append(secondBestRatio/bestRatio, 4);
//			sb.append("_badr=").append(secondBestBad/bestBad, 4);
//			sb.append("_goodr=").append(secondBestGood/bestGood, 4);
			
			float lenSum=a.length()+b.length();

			assert(sb.length()<40) : sb;
			sb.append((a.length()-100)*0.01f, 2).tab();
			sb.append((b.length()-100)*0.01f, 2).tab();
			sb.append(bestInsert*.004f, 3).tab();
			sb.append(bestOverlap/(bestOverlap+50f), 3).tab();
			sb.append((bestBad+1f)/(bestBad+5f), 5).tab();
			sb.append((bestGood+1f)/(bestGood+5f), 5).tab();
			sb.append(bestRatio, 4).tab();
			sb.append((bestBadInt+1f)/(bestBadInt+5f), 5).tab();

			sb.append(secondBestInsert*.004f, 3).tab();
			sb.append(secondBestOverlap/(secondBestOverlap+50f), 3).tab();
			sb.append((secondBestBad+1f)/(secondBestBad+5f), 5).tab();
			sb.append((secondBestGood+1f)/(secondBestGood+5f), 5).tab();
			sb.append(secondBestRatio, 4).tab();
			sb.append((secondBestBadInt)/(secondBestBadInt+5f), 5).tab();
			
			sb.append((secondBestRatio+1)/(bestRatio+1), 3).tab();
			sb.append(secondBestBad/(bestBad+8), 5).tab();
			sb.append(secondBestGood/(bestGood+8), 5).tab();
			sb.append((bestOverlap+1f)/(secondBestOverlap+bestOverlap+1f), 5).tab();
			
//			a.id=sb.toString();
			
			rvector[2]=bestBadInt;//=(int)(bestBad+0.95f);
			rvector[4]=(ambig ? 1 : 0);
			
			a.mapScore=bestInsert;//TODO: For testing
			return (bestInsert<0 ? -1 : bestInsert);
		}
		
		if(!ambig && bestRatio>maxRatio){bestInsert=-1;}
		
		rvector[2]=bestBadInt;//=(int)(bestBad+0.95f);
		rvector[4]=(ambig ? 1 : 0);

//		System.err.println("***C : "+bestOverlap+", "+ambig+", "+bestBad+", "+(bestOverlap<0 ? -1 : alen+blen-bestOverlap)+", "+
//				(bestOverlap<0 ? -1 : (bestOverlap<alen && alen>=blen) ? bestOverlap+alen-blen : bestOverlap)+", "+alen+", "+blen);
		return (bestInsert<0 ? -1 : bestInsert);
	}
	
	protected static final int mateByOverlapRatioJava(Read a, Read b, int[] rvector,
			int minOverlap0, int minOverlap, int minInsert0, int minInsert, float maxRatio, final float minSecondRatio,
			final float margin, final float offset, final float gIncr, final float bIncr) {
		assert(rvector!=null);
		assert(margin>=1);
		minOverlap=Tools.max(4, minOverlap0, minOverlap);
		minOverlap0=Tools.mid(4, minOverlap0, minOverlap);
		if(rvector==null){rvector=KillSwitch.allocInt1D(5);}
		
		final byte[] abases=a.bases, bbases=b.bases;
		final int alen=abases.length, blen=bbases.length;
		final int minLength=Tools.min(alen, blen);
		if(MAKE_VECTOR || TAG_CUSTOM){maxRatio=.7f;}
		{
			float x=findBestRatio(a, b, minOverlap0, minOverlap, minInsert, maxRatio, offset, gIncr, bIncr);
			if(verbose){
				System.err.println(x+", "+maxRatio+", "+Arrays.toString(rvector));
			}
			if(!TAG_CUSTOM && x>maxRatio){
				rvector[2]=minLength;
				rvector[4]=0;
				return -1;
			}
			maxRatio=Tools.min(maxRatio, x);
		}

//		final float altBadlimit=Tools.max(maxRatio, 0.07f)*2f*alen+1;
		final float margin2=(margin+offset)/minLength;
		final byte N='N';
		
		int bestInsert=-1;
		int bestOverlap=-1;
		int bestBadInt=-1;
		float bestBad=minLength;
		float bestGood=0;
		float bestRatio=1;
		boolean ambig=false;

		float secondBestBad=0;
		float secondBestGood=0;
		float secondBestRatio=1;
		int secondBestInsert=0;
		int secondBestOverlap=0;
		int secondBestBadInt=-1;
		
		final float extraMult=(TAG_CUSTOM || MAKE_VECTOR ? 4f : 1.2f);
		
		final int largestInsertToTest=(alen+blen-minOverlap0);
		final int smallestInsertToTest=minInsert0;
		for(int insert=largestInsertToTest; insert>=smallestInsertToTest; insert--){
//			if(verbose){System.err.println("Testing read "+a.numericID+", overlap "+insert+", insert "+(alen+blen-insert));}
			
			final int istart=(insert<=blen ? 0 : insert-blen);
			final int jstart=(insert>=blen ? 0 : blen-insert);
			final int overlapLength=Tools.min(alen-istart, blen-jstart, insert);
			
//			final float badlimit=Tools.min(altBadlimit, Tools.min(bestRatio, maxRatio)*margin*overlapLength);
			final float badlimit=extraMult*(Tools.min(bestRatio, maxRatio)*margin*overlapLength)+1f+extraBadlimit;
			float good=0, bad=0;
			int badInt=0;

			final int imax=istart+overlapLength;
			for(int i=istart, j=jstart; i<imax && bad<=badlimit; i++, j++){
				assert(i>=0 && i<alen && j>=0 && j<blen) : "\njstart="+jstart+", j="+j+
				", istart="+istart+", i="+i+" \n"+"insert="+insert+", overlap="+overlapLength+", a.length="+a.length()+
				", b.length="+b.length()+", bad="+bad+", badlimit="+badlimit+", good="+good;
				final byte ca=abases[i], cb=bbases[j];
				
				if(ca==cb){
					if(ca!=N){good+=gIncr;}
				}else{
					bad+=bIncr;
					badInt++;
				}
			}
			
//			if(verbose || true){
//				System.err.println("istart="+istart+", jstart="+jstart+", overlapLength="+overlapLength+", overlap="+overlap+", bestOverlap="+bestOverlap);
//				System.err.println("overlap="+overlap+", bad="+bad+", good="+good);
//				System.err.println("bestGood="+bestGood+", bestBad="+bestBad);
//				System.err.println();
//			}
			
			if(bad<=badlimit){
				if(!TAG_CUSTOM && bad==0 && good>minOverlap0 && good<minOverlap){
					rvector[2]=bestBadInt;//=(int)(bestBad+0.95f);
					rvector[4]=1;
					return -1;
				}

				float ratio=(bad+offset)/overlapLength;
//				System.err.println("*** ratio="+ratio+", bestRatio="+bestRatio);

				
				if(ratio<bestRatio*margin){

					ambig=(ratio*margin>=bestRatio || good<minOverlap);
					
					if(ratio<bestRatio){

						if(verbose){
							System.err.println("A: Set ambig="+ambig+": "+ratio+"*"+margin+"="+(ratio*margin)+">="+bestRatio+" || "+good+"<"+minOverlap);
						}
						secondBestInsert=bestInsert;
						secondBestOverlap=bestOverlap;
						secondBestBad=bestBad;
						secondBestRatio=bestRatio;
						secondBestBadInt=bestBadInt;
						
						bestInsert=insert;
						bestOverlap=overlapLength;
						bestBad=bad;
						bestRatio=ratio;
						bestBadInt=badInt;
					}
					else if(ratio<secondBestRatio){

						if(verbose){
							System.err.println("B: Set ambig="+ambig+": "+ratio+"*"+margin+"="+(ratio*margin)+">="+bestRatio+" || "+good+"<"+minOverlap);
						}
						secondBestInsert=insert;
						secondBestOverlap=overlapLength;
						secondBestBad=bad;
						secondBestRatio=ratio;
						secondBestBadInt=badInt;
					}
					
					if(!TAG_CUSTOM && ((ambig && bestRatio<margin2) || secondBestRatio<minSecondRatio)){
						rvector[2]=bestBadInt;//=(int)(bestBad+0.95f);
						rvector[4]=1;
						return -1;
					}
				}
			}
		}
		
		if(verbose){
			System.err.println("minSecondRatio="+minSecondRatio+", secondBestRatio="+secondBestRatio+", margin="+margin+", bestRatio="+bestRatio+", bestBad="+bestBad);
		}
		

		if(MAKE_VECTOR) {
//			System.err.println("Read "+a.numericID+": "+a.obj+", "+b.obj);
			FloatList fl=b.fetchVector();
			if(fl.size==3) {
				fl.add((a.length()-100)*0.01f); //alen
				fl.add((b.length()-100)*0.01f);
				fl.add(bestInsert*.004f);
				fl.add(bestOverlap/(bestOverlap+50f));
				fl.add((bestBad+1f)/(bestBad+5f));
				fl.add((bestGood+1f)/(bestGood+5f));
				fl.add(bestRatio);
				fl.add((bestBadInt+1f)/(bestBadInt+5f));

				fl.add(secondBestInsert*.004f);
				fl.add(secondBestOverlap/(secondBestOverlap+50f));
				fl.add((secondBestBad+1f)/(secondBestBad+5f));
				fl.add((secondBestGood+1f)/(secondBestGood+5f));
				fl.add(secondBestRatio);
				fl.add((secondBestBadInt)/(secondBestBadInt+5f));

				fl.add((secondBestRatio+1)/(bestRatio+1));
				fl.add(secondBestBad/(bestBad+8));
				fl.add(secondBestGood/(bestGood+8));
				fl.add((bestOverlap+1f)/(secondBestOverlap+bestOverlap+1f));
			}
		}
		
		if(TAG_CUSTOM){
			final ByteBuilder sb=(ByteBuilder) a.fetchBB();
			if(secondBestRatio<minSecondRatio){ambig=true;}

//			assert(sb.length()<40) : sb;
			if(sb.length()<40) {
				sb.append((a.length()-100)*0.01f, 2).tab();
				sb.append((b.length()-100)*0.01f, 2).tab();
				sb.append(bestInsert*.004f, 3).tab();
				sb.append(bestOverlap/(bestOverlap+50f), 3).tab();
				sb.append((bestBad+1f)/(bestBad+5f), 5).tab();
				sb.append((bestGood+1f)/(bestGood+5f), 5).tab();
				sb.append(bestRatio, 4).tab();
				sb.append((bestBadInt+1f)/(bestBadInt+5f), 5).tab();

				sb.append(secondBestInsert*.004f, 3).tab();
				sb.append(secondBestOverlap/(secondBestOverlap+50f), 3).tab();
				sb.append((secondBestBad+1f)/(secondBestBad+5f), 5).tab();
				sb.append((secondBestGood+1f)/(secondBestGood+5f), 5).tab();
				sb.append(secondBestRatio, 4).tab();
				sb.append((secondBestBadInt)/(secondBestBadInt+5f), 5).tab();
				
				sb.append((secondBestRatio+1)/(bestRatio+1), 3).tab();
				sb.append(secondBestBad/(bestBad+8), 5).tab();
				sb.append(secondBestGood/(bestGood+8), 5).tab();
				sb.append((bestOverlap+1f)/(secondBestOverlap+bestOverlap+1f), 5).tab();
			}
			
//			a.id=sb.toString();
			
			rvector[2]=bestBadInt;//=(int)(bestBad+0.95f);
			rvector[4]=(ambig ? 1 : 0);
			
			return (bestInsert<0 ? -1 : bestInsert);
		}
		
		if(!ambig && bestRatio>maxRatio){bestInsert=-1;}
		
		rvector[2]=bestBadInt;//=(int)(bestBad+0.95f);
		rvector[4]=(ambig ? 1 : 0);

//		System.err.println("***C : "+bestOverlap+", "+ambig+", "+bestBad+", "+(bestOverlap<0 ? -1 : alen+blen-bestOverlap)+", "+
//				(bestOverlap<0 ? -1 : (bestOverlap<alen && alen>=blen) ? bestOverlap+alen-blen : bestOverlap)+", "+alen+", "+blen);
		
		return (bestInsert<0 ? -1 : bestInsert);
	}
	
	protected static final float findBestRatio_WithQualities(Read a, Read b, final float[] aprob, final float[] bprob,
			final int minOverlap0, final int minOverlap, final int minInsert, final float maxRatio, final float offset) {
		final byte[] abases=a.bases, bbases=b.bases;
		final int alen=abases.length, blen=bbases.length;
		
		float bestRatio=maxRatio+0.0001f;
//		final float altBadlimit=Tools.max(maxRatio, 0.07f)*2f*alen+1;
		final float halfmax=maxRatio*0.5f;
		
		
		final int largestInsertToTest=(alen+blen-minOverlap); //TODO: test speed with minOverlap0
		final int smallestInsertToTest=minInsert;
		for(int insert=largestInsertToTest; insert>=smallestInsertToTest; insert--){
			if(verbose){System.err.println("a\nTesting read "+a.numericID+", overlap "+insert+", insert "+(alen+blen-insert));}
			
			final int istart=(insert<=blen ? 0 : insert-blen);
			final int jstart=(insert>=blen ? 0 : blen-insert);
			final int overlapLength=Tools.min(alen-istart, blen-jstart, insert);

//			final float badlimit=(Tools.min(altBadlimit, bestRatio*overlapLength));
			final float badlimit=bestRatio*overlapLength+extraBadlimit;
			float good=0, bad=0;
			
			final int imax=istart+overlapLength;
			for(int i=istart, j=jstart; i<imax && bad<=badlimit; i++, j++){
				assert(i>=0 && i<alen && j>=0 && j<blen) : "\njstart="+jstart+", j="+j+
				", istart="+istart+", i="+i+" \n"+"insert="+insert+", overlap="+overlapLength+", a.length="+a.length()+
				", b.length="+b.length()+", bad="+bad+", badlimit="+badlimit+", good="+good;
				final byte ca=abases[i], cb=bbases[j];
				
				final float x=aprob[i]*bprob[j];

				if(ca==cb){good+=x;}
				else{bad+=x;}
			}
			
			if(bad<=badlimit){
				if(bad==0 && good>minOverlap0 && good<minOverlap){
					return 100f;
				}
				
				float ratio=(bad+offset)/overlapLength;

				if(ratio<bestRatio){
					bestRatio=ratio;
					if(good>=minOverlap && ratio<halfmax){return bestRatio;}
				}
			}
		}
		
		return bestRatio;
	}
	
	protected static final float findBestRatio_WithQualitiesTraining(Read a, Read b, final float[] aprob, final float[] bprob,
			final int minOverlap0, final int minOverlap, final int minInsert, final float maxRatio, final float offset) {
		final byte[] abases=a.bases, bbases=b.bases;
		final int alen=abases.length, blen=bbases.length;
		
		float bestRatio=maxRatio+0.0001f;
		float bestRatio2=maxRatio+0.0001f;
		float bestRatio3=maxRatio+0.0001f;
//		final float altBadlimit=Tools.max(maxRatio, 0.07f)*2f*alen+1;
		final float halfmax=maxRatio*0.5f;
		
		int maxbad=99999;
		int maxbad2=99999;
		
		ByteBuilder bb=new ByteBuilder();//TODO, need expected length
		
		final int largestInsertToTest=(alen+blen-minOverlap); //TODO: test speed with minOverlap0
		final int smallestInsertToTest=minInsert;
		for(int insert=largestInsertToTest; insert>=smallestInsertToTest; insert--){
			if(verbose){System.err.println("a\nTesting read "+a.numericID+", overlap "+insert+", insert "+(alen+blen-insert));}
			
			final int istart=(insert<=blen ? 0 : insert-blen);
			final int jstart=(insert>=blen ? 0 : blen-insert);
			final int overlapLength=Tools.min(alen-istart, blen-jstart, insert);

//			final float badlimit=bestRatio*overlapLength*1.5f;//Made it a little bigger
			final float badlimit=bestRatio3*overlapLength+extraBadlimit;
			float good=0, bad=0;
			
			final int imax=istart+overlapLength;
			for(int i=istart, j=jstart; i<imax && bad<=badlimit; i++, j++){
				assert(i>=0 && i<alen && j>=0 && j<blen) : "\njstart="+jstart+", j="+j+
				", istart="+istart+", i="+i+" \n"+"insert="+insert+", overlap="+overlapLength+", a.length="+a.length()+
				", b.length="+b.length()+", bad="+bad+", badlimit="+badlimit+", good="+good;
				final byte ca=abases[i], cb=bbases[j];
				
				final float x=aprob[i]*bprob[j];

				if(ca==cb){good+=x;}
				else{bad+=x;}
			}
			
			float ratio=(bad+offset)/overlapLength;

			if(ratio<bestRatio){
				bestRatio3=bestRatio2;
				bestRatio2=bestRatio;
				bestRatio=ratio;
			}
		}
		
//		bb.append(alen).tab().append(blen).tab().append(maxRatio, 5).tab().append(maxbad).tab().append(bestRatio, 5);
//		bb.tab().append(bestRatio2, 5).tab().append(bestRatio3, 5)
		
		assert(false) : bb.length();
		return bestRatio;
	}
	
	protected static final float findBestRatio(Read a, Read b,
			final int minOverlap0, final int minOverlap, final int minInsert, final float maxRatio, final float offset, final float gIncr, final float bIncr) {
		final byte[] abases=a.bases, bbases=b.bases;
		final int alen=abases.length, blen=bbases.length;
		
		float bestRatio=maxRatio+0.0001f;
//		final float altBadlimit=Tools.max(maxRatio, 0.07f)*2f*alen+1;
		final float halfmax=maxRatio*0.5f;
		final byte N='N';
		
		
		final int largestInsertToTest=(alen+blen-minOverlap); //TODO: test speed with minOverlap0
		final int smallestInsertToTest=minInsert;
		for(int insert=largestInsertToTest; insert>=smallestInsertToTest; insert--){
//			if(verbose){System.err.println("Testing read "+a.numericID+", overlap "+insert+", insert "+(alen+blen-insert));}
			
			final int istart=(insert<=blen ? 0 : insert-blen);
			final int jstart=(insert>=blen ? 0 : blen-insert);
			final int overlapLength=Tools.min(alen-istart, blen-jstart, insert);

//			final float badlimit=(Tools.min(altBadlimit, bestRatio*overlapLength));
			final float badlimit=bestRatio*overlapLength+extraBadlimit;
			float good=0, bad=0;
			
			final int imax=istart+overlapLength;
			for(int i=istart, j=jstart; i<imax && bad<=badlimit; i++, j++){
				assert(i>=0 && i<alen && j>=0 && j<blen) : "\njstart="+jstart+", j="+j+
				", istart="+istart+", i="+i+" \n"+"insert="+insert+", overlap="+overlapLength+", a.length="+a.length()+
				", b.length="+b.length()+", bad="+bad+", badlimit="+badlimit+", good="+good;
				final byte ca=abases[i], cb=bbases[j];
				
				if(ca==cb){
					if(ca!=N){good+=gIncr;}
				}else{bad+=bIncr;}
			}

			if(bad<=badlimit){
				if(bad==0 && good>minOverlap0 && good<minOverlap){
					return 100f;
				}
				
				float ratio=(bad+offset)/overlapLength;

				if(ratio<bestRatio){
					bestRatio=ratio;
					if(good>=minOverlap && ratio<halfmax){return bestRatio;}
				}
			}
		}
		
		return bestRatio;
	}
	
	protected static final int mateByOverlapJava_unrolled(Read a, Read b, float[] aprob, float[] bprob, int[] rvector,
			int minOverlap0, final int minOverlap, final int minInsert0, int margin, final int maxMismatches0, final int maxMismatches, final int minq) {
		assert(rvector!=null);
		minOverlap0=Tools.min(Tools.max(1, minOverlap0), minOverlap);
		assert(maxMismatches<=maxMismatches0);
		margin=Tools.max(margin, 0);
		assert(maxMismatches>=margin);
		
		final byte[] abases=a.bases, bbases=b.bases;
		final byte[] aqual=a.quality, bqual=b.quality;
		final int alen=abases.length, blen=bbases.length;
		
		int bestOverlap=-1;
		int bestGood=-1;
		int bestBad=maxMismatches0;
		
		boolean ambig=false;
		final int maxOverlap=alen+blen-Tools.max(minOverlap, minInsert0);
//		assert(false) : minOverlap+", "+maxOverlap;
		
		if(aqual!=null && bqual!=null){
			for(int i=0; i<aqual.length; i++){aprob[i]=probCorrect3[aqual[i]];}
			for(int i=0; i<bqual.length; i++){bprob[i]=probCorrect3[bqual[i]];}
		}else{
			for(int i=0; i<alen; i++){aprob[i]=0.98f;}
			for(int i=0; i<blen; i++){bprob[i]=0.98f;}
		}
		
		final float minprob=probCorrect3[Tools.mid(1, minq, 41)];

		for(int overlap=Tools.max(minOverlap0, 0); overlap<maxOverlap; overlap++){
			if(verbose){System.err.println("c\nTesting read "+a.numericID+", overlap "+overlap+", insert "+(alen+blen-overlap));}

			int good=0, bad=0;

			int istart=(overlap<=alen ? 0 : overlap-alen);
			int jstart=(overlap<=alen ? alen-overlap : 0);

			{
				final int iters=Tools.min(overlap-istart, blen-istart, alen-jstart);
				final int imax=istart+iters;
				final int badlim=bestBad+margin+extraBadlimit;

				for(int i=istart, j=jstart; i<imax && bad<=badlim; i++, j++){
					assert(j>=0 && j<=alen && i>=0 && i<=blen) : "\njstart="+jstart+", j="+j+
					", istart="+istart+", i="+i+" \n"+"overlap="+overlap+", a.length="+alen+
					", b.length="+blen+", bad="+bad+", badlim="+badlim+", good="+good;
					final byte ca1=abases[j], cb1=bbases[i];
					final float pc=aprob[j]*bprob[i];
					
					if(pc<=minprob){//do nothing
					}else if(ca1==cb1){good++;}
					else{bad++;}
				}

				if(verbose){
					final int overlapLen=(imax-istart);
					System.err.println("overlapLen="+overlapLen+"; coordinates ("+jstart+"-"+(jstart+overlapLen)+"), ("+istart+"-"+imax+")");
					System.err.println(new String(abases, jstart, overlapLen));
					System.err.println(new String(bbases, istart, overlapLen));
				}

				if(verbose){System.err.println("overlap="+overlap+", bad="+bad+", good="+good+", badlim="+badlim+", bestOverlap="+
						bestOverlap+", bestGood="+bestGood+", bestBad="+bestBad+", ambig="+ambig+", mino="+minOverlap+", mino0="+minOverlap0+
						", margin="+margin+", maxMismatches="+maxMismatches);}
			}

			if(bad*2<good){
				if(verbose){System.err.print("a");}
				if(good>minOverlap){//Candidate
					if(verbose){System.err.print("b");}
					if(bad<=bestBad){

						if(verbose){System.err.print("c");}
						if(bad<bestBad || (bad==bestBad && good>bestGood)){//Current winner
							if(verbose){System.err.print("d");}
							if(bestBad-bad<margin){ambig=true;}
							bestOverlap=overlap;
							bestBad=bad;
							bestGood=good;
						}else if(bad==bestBad){
							if(verbose){System.err.print("e");}
							ambig=true;
						}

						if(ambig && bestBad<margin){
							if(verbose){System.err.print("f");}
							rvector[2]=bestBad;
							rvector[4]=(ambig ? 1 : 0);
							return -1;
						}
					}
				}else if(bad<margin){
					if(verbose){System.err.print("g");}
					ambig=true;
					rvector[2]=bestBad;
					rvector[4]=(ambig ? 1 : 0);
					return -1;
				}else{
					if(verbose){System.err.print("h");}
				}
			}
		}
		
		if(!ambig && bestBad>maxMismatches-margin){bestOverlap=-1;}
		
		rvector[2]=bestBad;
		rvector[4]=(ambig ? 1 : 0);
		
		if(verbose){System.err.println("bestOverlap="+
				bestOverlap+", bestGood="+bestGood+", bestBad="+bestBad+", ambig="+ambig+", mino="+minOverlap+", mino0="+minOverlap0+
				", margin="+margin+", maxMismatches="+maxMismatches);}
		
		return (bestOverlap<0 ? -1 : alen+blen-bestOverlap);
	}
	
	protected static final int mateByOverlapJava_unrolled_training(Read a, Read b, float[] aprob, float[] bprob, int[] rvector,
			int minOverlap0, final int minOverlap, final int minInsert0, int margin, final int maxMismatches0, final int maxMismatches, final int minq) {
		assert(rvector!=null);
		minOverlap0=Tools.min(Tools.max(1, minOverlap0), minOverlap);
		assert(maxMismatches<=maxMismatches0);
		margin=Tools.max(margin, 0);
		assert(maxMismatches>=margin);
		
		final byte[] abases=a.bases, bbases=b.bases;
		final byte[] aqual=a.quality, bqual=b.quality;
		final int alen=abases.length, blen=bbases.length;
		
		int bestOverlap=-1;
		int bestGood=-1;
		int bestBad=maxMismatches0;
		
		boolean ambig=false;
		final int maxOverlap=alen+blen-Tools.max(minOverlap, minInsert0);
//		assert(false) : minOverlap+", "+maxOverlap;
		assert(TAG_CUSTOM);
		
		if(aqual!=null && bqual!=null){
			for(int i=0; i<aqual.length; i++){aprob[i]=probCorrect3[aqual[i]];}
			for(int i=0; i<bqual.length; i++){bprob[i]=probCorrect3[bqual[i]];}
		}else{
			for(int i=0; i<alen; i++){aprob[i]=0.98f;}
			for(int i=0; i<blen; i++){bprob[i]=0.98f;}
		}
		
		final float minprob=probCorrect3[Tools.mid(1, minq, 41)];

		for(int overlap=Tools.max(minOverlap0, 0); overlap<maxOverlap; overlap++){
			if(verbose){System.err.println("c\nTesting read "+a.numericID+", overlap "+overlap+", insert "+(alen+blen-overlap));}

			int good=0, bad=0;

			int istart=(overlap<=alen ? 0 : overlap-alen);
			int jstart=(overlap<=alen ? alen-overlap : 0);

			{
				final int iters=Tools.min(overlap-istart, blen-istart, alen-jstart);
				final int imax=istart+iters;
				final int badlim=bestBad+margin+extraBadlimit;

				for(int i=istart, j=jstart; i<imax && bad<=badlim; i++, j++){
					assert(j>=0 && j<=alen && i>=0 && i<=blen) : "\njstart="+jstart+", j="+j+
					", istart="+istart+", i="+i+" \n"+"overlap="+overlap+", a.length="+alen+
					", b.length="+blen+", bad="+bad+", badlim="+badlim+", good="+good;
					final byte ca1=abases[j], cb1=bbases[i];
					final float pc=aprob[j]*bprob[i];
					
					if(pc<=minprob){//do nothing
					}else if(ca1==cb1){good++;}
					else{bad++;}
				}

				if(verbose){
					final int overlapLen=(imax-istart);
					System.err.println("overlapLen="+overlapLen+"; coordinates ("+jstart+"-"+(jstart+overlapLen)+"), ("+istart+"-"+imax+")");
					System.err.println(new String(abases, jstart, overlapLen));
					System.err.println(new String(bbases, istart, overlapLen));
				}

				if(verbose){System.err.println("overlap="+overlap+", bad="+bad+", good="+good+", badlim="+badlim+", bestOverlap="+
						bestOverlap+", bestGood="+bestGood+", bestBad="+bestBad+", ambig="+ambig+", mino="+minOverlap+", mino0="+minOverlap0+
						", margin="+margin+", maxMismatches="+maxMismatches);}
			}

			if(bad*0.0f<good){
				if(verbose){System.err.print("a");}
				if(good>minOverlap){//Candidate
					if(verbose){System.err.print("b");}
					if(bad<=bestBad){

						if(verbose){System.err.print("c");}
						if(bad<bestBad || (bad==bestBad && good>bestGood)){//Current winner
							if(verbose){System.err.print("d");}
							if(bestBad-bad<margin){ambig=true;}
							bestOverlap=overlap;
							bestBad=bad;
							bestGood=good;
						}else if(bad==bestBad){
							if(verbose){System.err.print("e");}
							ambig=true;
						}

						if(ambig && bestBad<margin){
							if(verbose){System.err.print("f");}
							rvector[2]=bestBad;
							rvector[4]=(ambig ? 1 : 0);
							return -1;
						}
					}
				}else if(bad<margin){
					if(verbose){System.err.print("g");}
					ambig=true;
					rvector[2]=bestBad;
					rvector[4]=(ambig ? 1 : 0);
					return -1;
				}else{
					if(verbose){System.err.print("h");}
				}
			}
		}
		
		if(!ambig && bestBad>maxMismatches-margin){bestOverlap=-1;}
		
		rvector[2]=bestBad;
		rvector[4]=(ambig ? 1 : 0);
		
		if(verbose){System.err.println("bestOverlap="+
				bestOverlap+", bestGood="+bestGood+", bestBad="+bestBad+", ambig="+ambig+", mino="+minOverlap+", mino0="+minOverlap0+
				", margin="+margin+", maxMismatches="+maxMismatches);}
		
		return (bestOverlap<0 ? -1 : alen+blen-bestOverlap);
	}
	
	
	/**
	 * TODO Use this
	 * @param a
	 * @param b
	 * @param overlap
	 * @return
	 */
	public static final float expectedMismatches(Read a, Read b, int overlap) {
		
		final byte[] abases=a.bases, bbases=b.bases, aqual=a.quality, bqual=b.quality;
		final int alen=abases.length, blen=bbases.length;
		final int istart=(overlap<=blen ? 0 : overlap-blen);
		final int jstart=(overlap<=alen ? alen-overlap : 0);
		
		float expected=0;
		float actual=0;
		
		if(aqual==null || bqual==null){return (overlap+0)/16;}

//		System.err.println(istart);
//		System.err.println(jstart);
//		System.err.println();
//
//		System.err.println(a.id);
//		System.err.println(overlap);
//		System.err.println(new String(a.bases));
//		System.err.println(new String(b.bases));
//		System.err.println();
//		for(int i=istart, j=jstart; i<overlap && i<alen && j<blen; i++, j++){
//			final byte ca=abases[i];
//			System.err.print((char)ca);
//		}
//		System.err.println();
//		for(int i=istart, j=jstart; i<overlap && i<alen && j<blen; i++, j++){
//			final byte cb=bbases[j];
//			System.err.print((char)cb);
//		}
//		System.err.println();
		
		for(int i=istart, j=jstart; i<overlap && i<alen && j<blen; i++, j++){
			final byte ca=abases[i], cb=bbases[j];
			final byte qa=aqual[i], qb=bqual[j];
			
			if(ca=='N' || cb=='N'){
				//do nothing
			}else{
				assert(AminoAcid.isFullyDefined(ca) && AminoAcid.isFullyDefined(cb)) :
					"A non-ACGTN base was detected.  Please rerun with the flag 'itn'.\n"+(char)ca+", "+(char)cb+"\n";
				float probC=probCorrect4[qa]*probCorrect4[qb];
				float probE=1-probC;
//				expected+=Tools.max(0.0005f, probE);
				expected+=probE;
				actual+=(ca==cb ? 0 : probC);
//				assert((probE==1) == (ca=='N' || cb=='N')) : ((char)ca)+", "+((char)cb)+", "+qa+", "+qb+", "+probC+", "+probE;
			}
		}
		
//		System.err.println("*expected:   \t"+expected);
//		System.err.println("*Actual:     \t"+actual);
//		System.err.println();
//
//		assert(a.id.equals("insert=142 /1") || a.id.equals("insert=263 /1")) : a.id;
		
		return expected;
	}
	
	/** Attempt at quantifying probability of an event like this.
	 * TODO: This returns an incorrect answer if reads are unequal lengths. */
	public static final float probability(Read a, Read b, int insert) {
		final byte[] abases=a.bases, bbases=b.bases, aqual=a.quality, bqual=b.quality;
		final int alen=abases.length, blen=bbases.length;
		final int istart=(insert<=blen ? 0 : insert-blen);
		final int jstart=(insert>=blen ? 0 : blen-insert);
		
		if(aqual==null || bqual==null){return 1;}
		
		float probActual=1;
		float probCommon=1;
		
//		float expected=0;
//		float actual=0;
//		int measuredOverlap=0;
		
//		assert(false) : "\n"+a.toFastq()+"\n"+b.toFastq()+"\n"+"istart="+istart+", jstart="+jstart+", insert="+insert+", alen="+alen+", blen="+blen;
		
		for(int i=istart, j=jstart; i<insert && i<alen && j<blen; i++, j++){
			final byte ca=abases[i], cb=bbases[j];
			final byte qa=aqual[i], qb=bqual[j];
			
			if(ca=='N' || cb=='N'){
				//do nothing
			}else{
				
//				System.err.println(((char)ca)+", "+((char)cb)+", "+i+", "+j);
				
				assert(AminoAcid.isFullyDefined(ca) && AminoAcid.isFullyDefined(cb)) :
					"A non-ACGTN base was detected.  Please rerun with the flag 'itn'.\n"+(char)ca+", "+(char)cb+"\n";
				float probC=probCorrect4[qa]*probCorrect4[qb];
				float probM=probC+(1-probC)*0.25f; //probability of matching
				float probE=1-probM;
				
				assert(probM>0) : qa+", "+qb+", "+probC+", "+probM+", "+probE;
				assert(probE>0) : qa+", "+qb+", "+probC+", "+probM+", "+probE;
				
				probCommon*=Tools.max(probM, probE);
				probActual*=(ca==cb ? probM : probE);
				
//				expected+=probE;
//				actual+=(ca==cb ? 0 : probM);
//				measuredOverlap++;
			}
		}
		
//		if(probActual>probCommon){
//			System.err.println("expected:   \t"+expected);
//			System.err.println("Actual:     \t"+actual);
//			System.err.println("probCommon: \t"+probCommon);
//			System.err.println("probActual: \t"+probActual);
//			System.err.println();
//			assert(false) : "\n"+a.toFastq()+"\n"+b.toFastq()+"\n";
//		}
		
		assert(probActual<=probCommon);

		return (float)Math.sqrt(probActual/probCommon); //sqrt is just so people don't need to type so many zeros.
	}
	
	protected static int minCoverage(final Read r, final Tadpole tadpole, final int k, int cutoff){
		if(k<32){
			return minCoverage(r, (Tadpole1)tadpole, k, cutoff);
		}else{
			return minCoverage(r, (Tadpole2)tadpole, k, cutoff);
		}
	}
	
	protected static int minCoverage(final Read r, final Tadpole1 tadpole, final int k, int cutoff){
		final byte[] bases=r.bases;
		if(bases==null || bases.length<k){return cutoff;}
		
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;
		int len=0;
		int min=cutoff;
		
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			final long x2=AminoAcid.baseToComplementNumber[b];

			//Update kmers
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

			//Handle Ns
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{len++;}

			if(len>=k){
				int cov=tadpole.getCount(kmer, rkmer);
				min=Tools.min(min, cov);
				if(min<cutoff){return min;}
			}
		}
		
		return min;
	}
	
	protected static int minCoverage(final Read r, final Tadpole2 tadpole, final int k, int cutoff){
		final byte[] bases=r.bases;
		if(bases==null || bases.length<k){return cutoff;}
		
		Kmer kmer=new Kmer(k);
		assert(kmer!=null);
		int min=cutoff;
		
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];

			//Update kmers
			kmer.addRight(b);

			if(kmer.len>=k){
				int cov=tadpole.getCount(kmer);
				min=Tools.min(min, cov);
				if(min<cutoff){return min;}
			}
		}
		
		return min;
	}
	
	protected static int calcMinOverlapByEntropy(byte[] bases, int k, short[] counts, int minscore){
		return Tools.max(calcMinOverlapByEntropyTail(bases, k, counts, minscore), calcMinOverlapByEntropyHead(bases, k, counts, minscore));
//		return calcMinOverlapByEntropyTail(bases, k, counts, minscore);
	}
	
	protected static int calcMinOverlapByEntropyTail(byte[] bases, int k, short[] counts, int minscore){
		final int bits=2*k;
		final int mask=~((-1)<<(bits));
		int kmer=0, len=0, ones=0, twos=0;
		
		if(counts==null){
			counts=localKmerCounts.get();
			if(counts==null){
				counts=new short[1<<(bits)];
				localKmerCounts.set(counts);
			}
		}
		
		Arrays.fill(counts, (short)0);
		
		for(int i=0, j=bases.length-1; i<bases.length; i++, j--){
			if(i<bases.length){
				final byte b=bases[j];
				if(!AminoAcid.isFullyDefined(b)){
					len=0;
					kmer=0;
				}else{
					len++;
					final int n=Dedupe.baseToNumber[b];
					kmer=((kmer<<2)|n)&mask;

					if(len>=k){
						counts[kmer]++;
						if(counts[kmer]==1){ones++;}
						else if(counts[kmer]==2){twos++;}
						if(ones*4+twos>=minscore){return i;}
					}
				}
			}
		}
		return bases.length+1;
	}
	
	protected static int calcMinOverlapByEntropyHead(byte[] bases, int k, short[] counts, int minscore){
		final int bits=2*k;
		final int mask=~((-1)<<(bits));
		int kmer=0, len=0, ones=0, twos=0;
		
		if(counts==null){
			counts=localKmerCounts.get();
			if(counts==null){
				counts=new short[1<<(bits)];
				localKmerCounts.set(counts);
			}
		}
		
		Arrays.fill(counts, (short)0);
		
		for(int i=0; i<bases.length; i++){
			if(i<bases.length){
				final byte b=bases[i];
				if(!AminoAcid.isFullyDefined(b)){
					len=0;
					kmer=0;
				}else{
					len++;
					final int n=Dedupe.baseToNumber[b];
					kmer=((kmer<<2)|n)&mask;

					if(len>=k){
						counts[kmer]++;
						if(counts[kmer]==1){ones++;}
						else if(counts[kmer]==2){twos++;}
						if(ones*4+twos>=minscore){return i;}
					}
				}
			}
		}
		return bases.length+1;
	}
	
	private static ThreadLocal<short[]> localKmerCounts=new ThreadLocal<short[]>();
	
	private static final int BAD_MULT=6;
	private static final int GOOD_MULT_1=8;
	private static final int GOOD_MULT_2=400;

	private static boolean TAG_CUSTOM=BBMerge.TAG_CUSTOM;
	static boolean MAKE_VECTOR=BBMerge.MAKE_VECTOR;

	private static int extraBadlimit=20;//(TAG_CUSTOM ? 40 : 20);//0 is fine when not using nets
	
	protected static final boolean verbose=false;
	
	private static final float[] probCorrect3=
		{0.000f, 0.251f, 0.369f, 0.499f, 0.602f, 0.684f, 0.749f, 0.800f, 0.842f, 0.874f,
		 0.900f, 0.921f, 0.937f, 0.950f, 0.960f, 0.968f, 0.975f, 0.980f, 0.984f, 0.987f,
		 0.990f, 0.992f, 0.994f, 0.995f, 0.996f, 0.997f, 0.997f, 0.998f, 0.998f, 0.999f,
		 0.999f, 0.999f, 0.999f, 0.999f, 1, 1, 1, 1, 1, 1,
		 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	
	private static final float[] probCorrect4=
		{0.0000f, 0.2501f, 0.3690f, 0.4988f, 0.6019f, 0.6838f, 0.7488f, 0.8005f, 0.8415f, 0.8741f,
		 0.9000f, 0.9206f, 0.9369f, 0.9499f, 0.9602f, 0.9684f, 0.9749f, 0.9800f, 0.9842f, 0.9874f,
		 0.9900f, 0.9921f, 0.9937f, 0.9950f, 0.9960f, 0.9968f, 0.9975f, 0.9980f, 0.9984f, 0.9987f,
		 0.9990f, 0.9992f, 0.9994f, 0.9995f, 0.9996f, 0.9997f, 0.9997f, 0.9998f, 0.9998f, 0.9999f,
		 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f,
		 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f, 0.9999f};
	
	private static final float[] probCorrect5=
		{0.20000f, 0.20567f, 0.36904f, 0.49881f, 0.60189f, 0.68377f, 0.74881f, 0.80047f, 0.84151f, 0.87411f,
		 0.90000f, 0.92057f, 0.93690f, 0.94988f, 0.96019f, 0.96838f, 0.97488f, 0.98005f, 0.98415f, 0.98741f,
		 0.99000f, 0.99206f, 0.99369f, 0.99499f, 0.99602f, 0.99684f, 0.99749f, 0.99800f, 0.99842f, 0.99874f,
		 0.99900f, 0.99921f, 0.99937f, 0.99950f, 0.99960f, 0.99968f, 0.99975f, 0.99980f, 0.99984f, 0.99987f,
		 0.99990f, 0.99992f, 0.99994f, 0.99995f, 0.99996f, 0.99997f, 0.99997f, 0.99998f, 0.99998f, 0.99999f,
		 0.99999f, 0.99999f, 0.99999f, 0.99999f, 0.99999f, 0.99999f, 0.99999f, 0.99999f, 0.99999f, 0.99999f};
	
}
