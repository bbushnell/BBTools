package bloom;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.BitSet;

import assemble.ErrorTracker;
import assemble.Rollback;
import dna.AminoAcid;
import kmer.AbstractKmerTable;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;
import ukmer.Kmer;

public class BloomFilterCorrector {
	
	public BloomFilterCorrector(BloomFilter filter_, int k_, int ksmall_) {
		filter=filter_;
		k=k_;
		ksmall=ksmall_;
		assert(ksmall<=k);
	}

	public int errorCorrect(Read r){
		initializeThreadLocals();
		int corrected=errorCorrect(r, localLeftCounts.get(), localRightCounts.get(), localLongList.get(),
				localIntList.get(), localIntList2.get(), localByteBuilder.get(), localByteBuilder2.get(), localTracker.get(), localBitSet.get());
		return corrected;
	}
	
	public int errorCorrect(Read r, final int[] leftCounts, final int[] rightCounts, LongList kmers, IntList counts, IntList counts2,
			final ByteBuilder bb, final ByteBuilder bb2, final ErrorTracker tracker, final BitSet bs, Kmer kmer, Kmer kmer2){
		return errorCorrect(r, leftCounts, rightCounts, kmers, counts, counts2, bb, bb2, tracker, bs);
	}
	
	public boolean hasErrorsFast(LongList kmers){
		if(kmers.size<1){return false;}
		int prev=-1;
		
		final int incr=Tools.mid(1, k/2, 9), mcc=minCountCorrect();
		for(int i=0; i<kmers.size; i+=incr){
			long kmer=kmers.get(i);
			if(kmer<0){
				return true;
			}
			long rkmer=rcomp(kmer);
			int count=getCount(kmer, rkmer);
			final int min=Tools.min(count, prev), max=Tools.max(count, prev);
			if(count<mcc || (i>0 && (isError(max+1, min-1)))){return true;}
			prev=count;
		}
		
		long kmer=kmers.get(kmers.size()-1);
		if(kmer<0){return true;}
		long rkmer=rcomp(kmer);
		int count=getCount(kmer, rkmer);
		final int min=Tools.min(count, prev), max=Tools.max(count, prev);
		return count<mcc || isError(max+1, min-1);
	}
	
	public int errorCorrect(Read r, final int[] leftCounts, final int[] rightCounts, LongList kmers, IntList counts, IntList counts2,
			final ByteBuilder bb, final ByteBuilder bb2, final ErrorTracker tracker, final BitSet bs){
		
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		tracker.clear();
		int valid=fillKmers(bases, kmers);
		if(valid<2){return 0;}
//		if(!r.containsUndefined() && !hasErrorsFast(kmers)){return 0;}//TODO: This does not seem useful
		

		fillCounts(bases, kmers, counts);
		final int possibleErrors=tracker.suspected=countErrors(counts, quals);
		if(possibleErrors<0){return 0;}
		final float expectedErrors=r.expectedErrors(true, r.length());
		final Rollback roll=ECC_ROLLBACK ? new Rollback(r, counts) : null;
		
		assert(counts.size>0) : counts.size;
		
		int correctedPincer=0;
		int correctedTail=0;
		int correctedBrute=0;
		int correctedReassemble=0;
		
		if(ECC_PINCER){
			correctedPincer+=errorCorrectPincer(bases, quals, leftCounts, rightCounts, kmers, counts, bb, tracker, errorExtensionPincer);
		}
		
		if(ECC_TAIL || ECC_ALL){
			int start=(ECC_ALL ? 0 : counts.size-k-1);
//			if(ECC_PINCER && tracker!=null && tracker.detected>correctedPincer){start=start-k;}
			correctedTail+=errorCorrectTail(bases, quals, leftCounts, rightCounts, kmers, counts, bb, tracker, start, errorExtensionTail);
			r.reverseComplement();
			valid=fillKmers(bases, kmers);
			counts.reverse();
			correctedTail+=errorCorrectTail(bases, quals, leftCounts, rightCounts, kmers, counts, bb, tracker, start, errorExtensionTail);
			r.reverseComplement();
			counts.reverse();
		}
		
		if(ECC_REASSEMBLE){
			if(verbose){System.err.println("Correcting "+possibleErrors+" errors.  Counts:\n"+counts);}
			if((correctedPincer<1 && correctedTail<1) || countErrors(counts, quals)>0){
				correctedReassemble=reassemble(bases, quals, rightCounts, counts, counts2, tracker, errorExtensionReassemble, bb, bb2, null, null, bs);
			}
			if(verbose){System.err.println("Corrected  "+correctedReassemble+" errors.  Counts:\n"+counts);}
		}
		assert(counts.size>0);
		
//		//123 For testing.
//		if(false && tracker.detected()>tracker.corrected()){
//			correctedBrute+=errorCorrectBruteForce(bases, quals, leftCounts, rightCounts, kmers, counts, bb, tracker, errorExtensionPincer);
//		}
		
		assert(correctedPincer+correctedTail+correctedReassemble+correctedBrute==tracker.corrected())
			: correctedPincer+", "+correctedTail+", "+correctedReassemble+", "+correctedBrute+", "+tracker;

		if(ECC_ROLLBACK && (tracker.corrected()>0 || tracker.rollback)){
			
			if(!tracker.rollback && quals!=null && tracker.corrected()>3){
				float mult=Tools.max(1, 0.5f*(0.5f+0.01f*r.length()));//1 for a 150bp read.
				if(countErrors(counts, quals)>0 && tracker.corrected()>mult+expectedErrors){tracker.rollback=true;}
				else if(tracker.corrected()>2.5f*mult+expectedErrors){tracker.rollback=true;}
			}
			
			IntList counts0=roll.counts0;
			for(int i=0; !tracker.rollback && i<counts.size; i++){
				int a=Tools.max(0, counts0.get(i));
				int b=Tools.max(0, counts.get(i));
				if(b<a-1 && !isSimilar(a, b)){
					if(verbose){System.err.println("Y: RID="+r.numericID+"; "+a+"->"+b+"\n"+counts0+"\n"+counts);}
					tracker.rollback=true;
				}
			}
			
			if(tracker.rollback){
				roll.rollback(r, counts);
				tracker.clearCorrected();
				return 0;
			}
		}
		
		if(MARK_BAD_BASES>0 && (!MARK_ERROR_READS_ONLY || countErrors(counts, quals)>0 ||
				r.expectedErrors(false, r.length())>3)){
			int marked=markBadBases(bases, quals, counts, bs, MARK_BAD_BASES, MARK_DELTA_ONLY, MARK_QUALITY);
			tracker.marked=marked;
		}
		
		return tracker.corrected();
	}
	
	/** Changes to N any base covered strictly by kmers with count below minCount */
	public final int markBadBases(final byte[] bases, final byte[] quals, final IntList counts, final BitSet bs,
			final int minCount, boolean deltaOnly, final byte markQuality){
		if(counts.size<1){return 0;}
		
		bs.clear();
		assert(counts.size==bases.length-k+1) : counts.size+", "+bases.length;
		
		for(int i=0; i<counts.size;){
			final int count=counts.get(i);
			if(count>=minCount){
				bs.set(i, i+k);
				i+=k;
			}else{
				i++;
			}
		}
		{//Last cycle
			final int i=counts.size-1;
			final int count=counts.get(i);
			if(count>=minCount){
				bs.set(i, i+k);
			}
		}
		
		final int card=bs.cardinality();
		final int toMark=bases.length-card;
		int marked=0;
		assert(card<=bases.length);
		
		int consecutiveBad=0;
		for(int i=0; i<bases.length; i++){
			if(bs.get(i)){
				consecutiveBad=0;
			}else{
				consecutiveBad++;
				boolean mark=((quals!=null && quals[i]>markQuality) || bases[i]!='N');
				if(mark && deltaOnly){
					mark=(consecutiveBad>=k) || bs.get(i+1) || (i>0 && bs.get(i-1));
				}
				if(mark){
					marked++;

					if(markQuality<1){
						bases[i]='N';
					}
					if(quals!=null){
						quals[i]=Tools.min(quals[i], (bases[i]=='N' ? 0 : markQuality));
					}
				}
				if(bases[i]=='N' || (quals!=null && quals[i]<=markQuality)){consecutiveBad=0;}
			}
		}
		
		return marked;
	}
	
	public void fillCounts(byte[] bases, LongList kmers, IntList counts){
		counts.clear();
		
		if(k==ksmall){
			fillCountsFromKmers(kmers, counts);
		}else{
			filter.fillCountsBig(bases, counts);
//			filter.fillCountsBig(kmers, counts);
		}
		
//		assert(counts.size==kmers.size) : counts.size+", "+kmers.size;
		if(smooth){
			smooth(kmers, counts, smoothWidth);
		}
	}
	
	private void fillCountsFromKmers(LongList kmers, IntList counts){
		for(int i=0; i<kmers.size; i++){
			long kmer=kmers.get(i);
			if(kmer>=0){
				long rkmer=rcomp(kmer);
				int count=getCount(kmer, rkmer);
				counts.add(count);
			}else{
				counts.add(0);
			}
		}
	}
	
	public void smooth(LongList kmerList, IntList countList, int width){
		final int size=countList.size;
		if(size<3){return;}
		final int[] counts=countList.array;
		final long[] kmers=kmerList.array;
		
//		for(int i=0; i<size; i++){
//			long kmer=kmers.get(i);
//			array[i]=Tools.min(array[i], Tools.max(maxLeftCount(kmer), maxRightCount(kmer)));
//		}
		
		counts[0]=Tools.min(counts[0], Tools.max(counts[1], maxLeftCount(kmers[0])));
		counts[size-1]=Tools.min(counts[size-1], Tools.max(counts[size-2], maxRightCount(kmers[size-1])));

		for(int i=1, max=size-1; i<max; i++){
			long ka=kmers[i-1], kc=kmers[i+1];
			if(ka>=0 && kc>=0){
				int a=counts[i-1], b=counts[i], c=counts[i+1];
				int maxCount=Tools.max(a, c);
				counts[i]=Tools.min(b, maxCount);
			}
		}
		
		//Smooth peaks 2-wide.
		if(width>=2){
			for(int i=1, max=size-2; i<max; i++){
				long ka=kmers[i-1], kd=kmers[i+2];
				if(ka>=0 && kd>=0){
					int a=counts[i-1], b=counts[i], c=counts[i+1], d=counts[i+2];
					int maxCount=Tools.max(a, d);
					counts[i]=Tools.min(b, maxCount);
					counts[i+1]=Tools.min(c, maxCount);
				}
			}
		}

		//Smooth peaks 3-wide.
		if(width>=3){
			boolean changed=false;
			for(int i=2, max=size-2; i<max; i++){
				long ka=kmers[i-2], kc=kmers[i+2];
				if(ka>=0 && kc>=0){
					int a=counts[i-2], b=counts[i], c=counts[i+2];
					int maxCount=Tools.max(a, c);
					if(maxCount<b){
						counts[i]=Tools.min(b, maxCount);
						changed=true;
					}
				}
			}
			if(changed){smooth(kmerList, countList, 2);}
		}
		
//		array[0]=Tools.min(array[0], array[1]);
//		array[size-1]=Tools.min(array[size-1], array[size-2]);
	}
	
	public final int reassemble(final byte[] bases, final byte[] quals, final int[] rightCounts, final IntList counts, final IntList counts2,
			final ErrorTracker tracker, final int errorExtension, final ByteBuilder bb, final ByteBuilder bb2, final Kmer kmer, final Kmer regenKmer, BitSet bs){
		if(bases.length<k+1+deadZone){return 0;}
		final ByteBuilder fromLeft=new ByteBuilder(bases.length);
		final ByteBuilder fromRight=new ByteBuilder(bases.length);
		
		int detected0=tracker.detectedReassemble;
		int corrected=reassemble_pass(bases, quals, fromLeft, fromRight, rightCounts, counts, counts2, tracker, errorExtension, kmer, regenKmer, bs);
		
		int correctedIncr=corrected;
		int detectedIncr=tracker.detectedReassemble-detected0;
		int uncorrected=detectedIncr-correctedIncr;
		
		for(int passes=1; passes<6 && correctedIncr>0 && uncorrected>0; passes++){//Without a pass limit this could, in rare cases, make an infinite loop
			tracker.detectedReassemble-=uncorrected;
			detected0=tracker.detectedReassemble;
			correctedIncr=reassemble_pass(bases, quals, fromLeft, fromRight, rightCounts, counts, counts2, tracker, errorExtension, kmer, regenKmer, bs);
			
			corrected+=correctedIncr;
			detectedIncr=tracker.detectedReassemble-detected0;
			uncorrected=detectedIncr-correctedIncr;
		}
		
		return corrected;
	}
	
	public final int reassemble_pass(final byte[] bases, final byte[] quals, final ByteBuilder fromLeft, final ByteBuilder fromRight,
			final int[] rightCounts, final IntList counts, final IntList counts2, final ErrorTracker tracker, final int errorExtension,
			final Kmer kmer, final Kmer kmer2, final BitSet bs){
		if(bases.length<k+1+deadZone){return 0;}

		fromLeft.clear();
		fromRight.clear();
		for(byte b : bases){
			fromLeft.append(b);
			fromRight.append(b);
		}
		
		assert(counts.size>0) : counts+", "+bases.length;
		
		counts2.clear();
		counts2.addAll(counts);
		reassemble_inner(fromLeft, quals, rightCounts, counts2, errorExtension, kmer, kmer2);
		
		fromRight.reverseComplementInPlace();
		counts2.clear();
		counts2.addAll(counts);
		counts2.reverse();
		
		reassemble_inner(fromRight, quals, rightCounts, counts2, errorExtension, kmer, kmer2);
		fromRight.reverseComplementInPlace();
		
//		System.err.println(new String(fromRight));
//		System.err.println(copy);
//		System.err.println();

		int correctedInner=0;
		int correctedOuter=0;
		int detectedInner=0;
		int detectedOuter=0;
		boolean rollback=false;
		
		for(int i=0; i<bases.length; i++){
			byte a=bases[i];
			byte b=fromLeft.get(i);
			byte c=fromRight.get(i);
			if(a!=b || a!=c){
				if(b==c){detectedInner++;}
				else{
					detectedOuter++;
					if(a!=b && a!=c){
						assert(b!=c);
						rollback=true;
					}
				}
			}
			if(b==a){fromLeft.set(i, (byte)0);}
			if(c==a){fromRight.set(i, (byte)0);}
		}
		
		final int detected=detectedInner+detectedOuter;
		tracker.detectedReassemble+=detected;
		if(rollback || detected==0){return 0;}
		bs.clear();
		
		int clearedLeft=clearWindow2(fromLeft, quals, windowLen, windowCount, windowQualSum/*, windowCountHQ, windowHQThresh*/);
		fromRight.reverseInPlace();
		Tools.reverseInPlace(quals);
		int clearedRight=clearWindow2(fromRight, quals, windowLen, windowCount, windowQualSum/*, windowCountHQ, windowHQThresh*/);
		fromRight.reverseInPlace();
		Tools.reverseInPlace(quals);
		
		for(int i=0; i<bases.length; i++){
			byte a=bases[i];
			byte b=fromLeft.get(i);
			byte c=fromRight.get(i);
			byte d=a;
			if(b==0 && c==0){
				//do nothing
			}else if(b==c){
				d=b;
			}else if(b==0){
				d=c;
			}else if(c==0){
				d=b;
			}else if(b!=c){
//				if(AminoAcid.isFullyDefined(a)){
//					quals[i]=(byte)Tools.max(2, q-3);
//				}
			}
			
			if(ECC_REQUIRE_BIDIRECTIONAL && b!=c && i>=k && i<bases.length-k){d=a;}//Clause to force pincer mode in the middle
			
			if(d!=a){
				byte q=(quals==null ? 30 : quals[i]);
				if(b==c){
					correctedInner++;
					q=(byte)Tools.mid(q+qIncreasePincer, qMinPincer, qMaxPincer);
				}else{
					correctedOuter++;
					q=(byte)Tools.mid(q+qIncreaseTail, qMinTail, qMaxTail);
				}
				if(!rollback){
					bs.set(i);
					bases[i]=d;
					if(quals!=null){quals[i]=q;}
				}
			}
		}
		
		if(rollback && correctedInner+correctedOuter>0){
			tracker.rollback=true;
			return 0;
		}
		
		{
			tracker.correctedReassembleInner+=correctedInner;
			tracker.correctedReassembleOuter+=correctedOuter;
		}
		int corrected=correctedOuter+correctedInner;
		
		if(corrected>0){
			regenerateCounts(bases, counts, kmer, bs);
			assert(counts.size>0);
		}
		
		return corrected;
	}
	
	public int regenerateCounts(byte[] bases, IntList counts, final Kmer dummy, BitSet changed){
		assert(!changed.isEmpty());
		final int firstBase=changed.nextSetBit(0), lastBase=changed.length()-1;
		final int ca=firstBase-k;
//		final int b=changed.nextSetBit(0);ca+kbig; //first base changed
		final int firstCount=Tools.max(firstBase-k+1, 0), lastCount=Tools.min(counts.size-1, lastBase); //count limit
//		System.err.println("ca="+ca+", b="+b+", lim="+lim);
//		System.err.println("Regen from count "+(ca+1)+"-"+lim);
		
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;
		int len=0;
		int valid=0;
		
//		System.err.println("ca="+ca+", b="+b+", lim="+lim+", "+counts);
		
		for(int i=Tools.max(0, firstBase-k+1), lim=Tools.min(lastBase+k-1, bases.length-1); i<=lim; i++){
			final byte base=bases[i];
			final long x=AminoAcid.baseToNumber[base];
			final long x2=AminoAcid.baseToComplementNumber[base];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{
				len++;
			}
			
			final int c=i-k+1;
			if(i>=firstBase){
				if(len>=k){
					valid++;
					int count=getCount(kmer, rkmer);
					counts.set(c, count);
				}else if(c>=0){
					counts.set(c, 0);
				}
			}
		}
		
		if(smooth){
			int[] array=counts.array;
			for(int i=1; i<counts.size-1; i++){
				int a=array[i-1], b=array[i], c=array[i+1];
				if(b>a && b>c){array[i]=Tools.max(a, c);}
			}
		}
		
		return valid;
	}
	
	private static int clearWindow2(final ByteBuilder bb, final byte[] quals, final int window,
			final int limit, final int qsumLimit/*, final int limitHQ, final byte hqThresh*/){
		final int len=bb.length;
		final byte[] array=bb.array;
		
		int cleared=0;
		int count=0, countHQ=0, qsum=0;
		for(int i=0, prev=-window; i<len; i++, prev++){
			byte b=array[i];
			
			if(b!=0 && (quals==null || quals[prev]>0)){
				count++;
				if(quals!=null){
					qsum+=quals[i];
//					if(quals!=null && quals[i]>=hqThresh){countHQ++;}
				}
				if(count>limit || qsum>qsumLimit /*|| countHQ>limitHQ*/){
					for(int j=Tools.max(0, i-window), lim=bb.length(); j<lim; j++){
						if(array[j]!=0){
							array[j]=0;
							cleared++;
						}
					}
					return cleared;
				}
			}
			if(prev>=0 && array[prev]>0 && (quals==null || quals[prev]>0)){
				count--;
				if(quals!=null){
					countHQ--;
//					 if(quals[prev]>=hqThresh){qsum-=quals[i];}
				}
			}
		}
		return cleared;
	}
	
	public int reassemble_inner(final ByteBuilder bb, final byte[] quals, final int[] rightCounts, final IntList counts,
			final int errorExtension, final Kmer kmer, final Kmer regenKmer){
		return reassemble_inner(bb, quals, rightCounts, counts, errorExtension);
	}
	
	public int reassemble_inner(final ByteBuilder bb, final byte[] quals, final int[] rightCounts, final IntList counts,
			final int errorExtension){
		final int length=bb.length();
		if(length<k+1+deadZone){return 0;}
		final byte[] bases=bb.array;
		
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;

		int detected=0;
		int corrected=0;
		int len=0;
		
		//a is the index of the rightmost base of the kmer
		//b=a+1 is the index of the next base
		//ca=a-k+1 is the index of the count of the kmer
		//cb=a-k+2 is the index of the count of the next kmer
		for(int a=0, lim=length-deadZone-1; a<lim; a++){
			
			//Generate the new kmer
			final byte aBase=bases[a];
			final long x=AminoAcid.baseToNumber[aBase];
			final long x2=AminoAcid.baseToComplementNumber[aBase];
			
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{


				//Now consider the next kmer
				kmer=((kmer<<2)|(long)x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				len++;

				if(verbose){
					System.err.println("len: "+len+" vs "+k+"; a="+a);
				}

				if(len>=k){

					final int b=a+1;
					final int ca=a-k+1;
					final int cb=ca+1;

					final int aCount=counts.get(ca);
					final int bCount=counts.get(cb);
					final byte qb=(quals==null ? 20 : quals[b]);

					if(verbose){
						System.err.println("ca="+ca+", cb="+cb+"; aCount="+aCount+", bCount="+bCount);
						System.err.println(isError(aCount, bCount, qb)+", "+isSimilar(aCount, ca-errorExtension, ca-1, counts)+
								", "+isError(aCount, ca+2, ca+k, counts));
					}

//					if(isError(aCount, bCount) && isSimilar(aCount, ca-errorExtension, ca-1, counts) && isError(aCount, ca+2, ca+k, counts)){
					if(isSubstitution(ca, errorExtension, qb, counts)){
						if(verbose){
							System.err.println("***Found error: "+aCount+", "+bCount);
						}
						//Assume a 1bp substitution; attempt to correct.

						int rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
						int rightMax=rightCounts[rightMaxPos];
						int rightSecondPos=Tools.secondHighestPosition(rightCounts);
						int rightSecond=rightCounts[rightSecondPos];

						byte base=bases[b];
						byte num=AminoAcid.baseToNumber[base];

						if(rightMax>=minCountExtend){
							detected++;
							if(num==rightMax){
								detected--;
								//							bases2[b]=base;
							}else if((isError(rightMax, rightSecond, qb) || !isJunction(rightMax, rightSecond)) && isSimilar(aCount, rightMax)){
								bases[b]=AminoAcid.numberToBase[rightMaxPos];
								corrected++;
								regenerateCounts(bases, counts, ca);
								if(verbose){System.err.println("Corrected error: "+num+"->"+rightMaxPos+". New counts:\n"+counts);}
							}

							//						else if(rightSecond>=minCountExtend && isJunction(rightMax, rightSecond) && isSimilar(aCount, rightSecond)
							//								&& !isSimilar(aCount, rightMax)){//This branch may not be very safe.
							//							bases2[b]=AminoAcid.numberToBase[rightSecondPos];
							//							corrected++;
							//							if(verbose){System.err.println("Corrected error.");}
							//						}
						}

					}else{
						if(verbose){
							System.err.println("Not an error: "+aCount+", "+bCount+
									";  "+isError(aCount, bCount, qb)+", "+isSimilar(aCount, a-errorExtension, a-1, counts)+", "+isError(aCount, a+2, a+k, counts));
						}
					}
				}
			}
		}
		
		return corrected;
	}
	
	protected final boolean isSubstitution(int ca, int errorExtension, byte qb, IntList counts){
		final int cb=ca+1;
		final int aCount=counts.get(ca);
		final int bCount=counts.get(cb);
		if(isError(aCount, bCount, qb) && isSimilar(aCount, ca-errorExtension, ca-1, counts) &&
				isError(aCount, ca+2, ca+k, counts)){
			final int cc=ca+k;
			final int cd=cc+1;
			if(cd<counts.size){
				final int cCount=counts.get(cc);
				final int dCount=counts.get(cd);
				if(isError(aCount, dCount) || isError(dCount, cCount, qb)){
					return true;
				}
			}else{return true;}
		}
		return false;
	}
	
	public final int countErrors(IntList counts, byte[] quals){
		int possibleErrors=0;
		for(int i=1; i<counts.size; i++){
			final int a=counts.get(i-1), b=counts.get(i);
			boolean error;
			if(quals!=null){
				error=isErrorBidirectional(a, b, quals[i-1], quals[i+k-1]);
			}else{
				error=isErrorBidirectional(a, b, (byte)20, (byte)20);
			}
			if(error){
				possibleErrors++;
				i+=k;
			}
		}
		return possibleErrors;
	}
	
	public int errorCorrectPincer(final byte[] bases, final byte[] quals, final int[] leftBuffer, final int[] rightBuffer,
			final LongList kmers, final IntList counts, final ByteBuilder bb, final ErrorTracker tracker, final int errorExtension){
		
		int detected=0;
		int corrected=0;
		
		//a is the index of the left kmer
		//b is a+1 (right-extension of left kmer)
		//c is d-1 (left-extension of right kmer)
		//d is the index of the right kmer
		//the base between the kmers is at a+k
		for(int a=0, d=k+1; d<counts.size; a++, d++){
			final int aCount=counts.get(a);
			final int bCount=counts.get(a+1);
			final int cCount=counts.get(d-1);
			final int dCount=counts.get(d);
			final byte qb=(quals==null ? 20 : quals[a+k]);
			if(isError(aCount, bCount, qb) && isError(dCount, cCount, qb) && isSimilar(aCount, dCount)){
				if(verbose){
					System.err.println("Found error: "+aCount+", "+bCount+", "+cCount+", "+dCount);
				}
				//Looks like a 1bp substitution; attempt to correct.
				detected++;
				int ret=correctSingleBasePincer(a, d, bases, quals, leftBuffer, rightBuffer, kmers, counts, bb, errorExtension);
				corrected+=ret;
				if(verbose){
					System.err.println("Corrected error.");
				}
			}else{
				if(verbose){
					System.err.println("Not an error: "+aCount+", "+bCount+", "+cCount+", "+dCount+
							";  "+isError(aCount, bCount, qb)+", "+isError(dCount, cCount, qb)+", "+isSimilar(aCount, dCount));
				}
			}
		}
		
//		if(detected==0 && counts.get(0)>2 && counts.get(counts.size-1)>2){
//			assert(!verbose);
//			verbose=true;
//			System.err.println("\n"+counts);
//			errorCorrectPincer(bases, quals, leftBuffer, rightBuffer, kmers, counts, bb, tracker);
//			assert(false);
//		}
		
		{
			tracker.detectedPincer+=detected;
			tracker.correctedPincer+=corrected;
		}
		
		return corrected;
	}

	public int errorCorrectTail(final byte[] bases, final byte[] quals, final int[] leftBuffer, final int[] rightBuffer,
			final LongList kmers, final IntList counts, final ByteBuilder bb, final ErrorTracker tracker, final int startPos, final int errorExtension){
		if(bases.length<k+2+errorExtension+deadZone){return 0;}
		int detected=0;
		int corrected=0;
		
		//a is the index of the left kmer
		//b is a+1
		//the base between the kmers is at a+k
		for(int a=Tools.max(startPos, errorExtension), lim=counts.size-deadZone-1; a<lim; a++){//errorExtension-1
			final int aCount=counts.get(a);
			final int bCount=counts.get(a+1);
			final byte qb=(quals==null ? 20 : quals[a+k]);
			if(isError(aCount, bCount, qb) && isSimilar(aCount, a-errorExtension, a-1, counts) && isError(aCount, a+2, a+k, counts)){
				if(verbose){
					System.err.println("Found error: "+aCount+", "+bCount);
				}
				//Assume like a 1bp substitution; attempt to correct.
				detected++;
				int ret=correctSingleBaseRight(a, bases, quals, leftBuffer, rightBuffer, kmers, counts, bb, errorExtension);
				corrected+=ret;
				if(verbose){
					System.err.println("Corrected error.");
				}
			}else{
				if(verbose){
					System.err.println("Not an error: "+aCount+", "+bCount+
							";  "+isError(aCount, bCount, qb)+", "+isSimilar(aCount, a-errorExtension, a-1, counts)+", "+isError(aCount, a+2, a+k, counts));
				}
			}
		}
		
//		if(detected==0 && counts.get(0)>2 && counts.get(counts.size-1)>2){
//			assert(!verbose);
//			verbose=true;
//			System.err.println("\n"+counts);
//			errorCorrectPincer(bases, quals, leftBuffer, rightBuffer, kmers, counts, bb, tracker);
//			assert(false);
//		}
		
		{
			tracker.detectedTail+=detected;
			tracker.correctedTail+=corrected;
		}
		
		return corrected;
	}
	
	private int correctSingleBasePincer(final int a, final int d, final byte[] bases, final byte[] quals, final int[] leftBuffer, final int[] rightBuffer,
			final LongList kmers, final IntList counts, final ByteBuilder bb, final int errorExtension){
		final byte leftReplacement, rightReplacement;
		final int loc=a+k;
		{
			bb.clear();
			final long kmer=kmers.get(a);
			final long rkmer=rcomp(kmer);
			int extension=extendToRight2(bb, null, rightBuffer, errorExtension, true, kmer, rkmer);
			if(extension<errorExtension){return 0;}
			for(int i=1; i<extension; i++){
				if(bb.get(i)!=bases[loc+i]){
					return 0;
				}
			}
			leftReplacement=bb.get(0);
		}
		{
			bb.clear();
			final long rkmer=kmers.get(d);
			final long kmer=rcomp(rkmer);
			int extension=extendToRight2(bb, null, rightBuffer, errorExtension, true, kmer, rkmer);
			if(extension<errorExtension){return 0;}
			bb.reverseComplementInPlace();
			for(int i=0; i<extension-1; i++){
				if(bb.get(i)!=bases[loc+i+1-extension]){
					return 0;
				}
			}
			rightReplacement=bb.get(extension-1);
		}
		if(leftReplacement!=rightReplacement){return 0;}
		if(bases[loc]==leftReplacement){return 0;}
		if(!isSimilar(a, leftReplacement, kmers, counts)){return 0;}
		
		bases[loc]=leftReplacement;
		assert(d==a+k+1);
		regenerateKmers(bases, kmers, counts, a);
		return 1;
	}
	
	private int correctSingleBaseRight(final int a, final byte[] bases, final byte[] quals, final int[] leftBuffer, final int[] rightBuffer,
			final LongList kmers, final IntList counts, final ByteBuilder bb, final int errorExtension0){
		final byte leftReplacement;
		final int loc=a+k;
		final int errorExtension=Tools.min(errorExtension0, bases.length-loc);
		{
			bb.clear();
			final long kmer=kmers.get(a);
			final long rkmer=rcomp(kmer);
			int extension=extendToRight2(bb, null, rightBuffer, errorExtension, true, kmer, rkmer);
			if(extension<errorExtension){return 0;}
			for(int i=1; i<extension; i++){
				if(bb.get(i)!=bases[loc+i]){
					return 0;
				}
			}
			leftReplacement=bb.get(0);
		}
		
		if(bases[loc]==leftReplacement){return 0;}
		if(!isSimilar(a, leftReplacement, kmers, counts)){return 0;}
		
		bases[loc]=leftReplacement;
		regenerateKmers(bases, kmers, counts, a);
		return 1;
	}
	
	private boolean isSimilar(int a, byte newBase, LongList kmers, IntList counts){
		final int shift=2*k;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=kmers.get(a);
				
		final long x=AminoAcid.baseToNumber[newBase];
		kmer=((kmer<<2)|x)&mask;
		long rkmer=rcomp(kmer);
		int count=getCount(kmer, rkmer);
		int aCount=counts.get(a);
		boolean similar=isSimilar(aCount, count);
		return similar;
	}
	
	protected final boolean isSimilar(final int a, int loc1, int loc2, final IntList counts){
		loc1=Tools.max(loc1, 0);
		loc2=Tools.min(loc2, counts.size-1);
		for(int i=loc1; i<=loc2; i++){
			if(!isSimilar(a, counts.get(i))){return false;}
		}
		return true;
	}
	
	protected final boolean isSimilar(final int a, final int b){
		int min=Tools.min(a, b);
		int max=Tools.max(a, b);
		int dif=max-min;
		assert(dif>=0);
		return (dif<pathSimilarityConstant || dif<max*pathSimilarityFraction);
	}
	
	protected final boolean isError(final int a, int loc1, int loc2, final IntList counts){
		loc1=Tools.max(loc1, 0);
		loc2=Tools.min(loc2, counts.size-1);
		for(int i=loc1; i<=loc2; i++){
			if(!isError(a, counts.get(i))){return false;}
		}
		return true;
	}
	
	protected final boolean isErrorBidirectional(final int a, final int b, final byte qa, final byte qb){
		return (a>=b ? isError(a, b, qb) : isError(b, a, qa));
	}
	
	protected final boolean isError(final int high, final int low){
		float em1;
		if(errorPath==1){
			em1=errorMult1;
		}else if(errorPath==2){
//		if(low<minCountCorrect && high>=minCountCorrect && high>2*low){return true;}
//		final float em1=4;//Tools.mid(errorMult1, 4, low-1); //4;//(low<3 ? 4 : Tools.min(errorMult1, 2*low));
			em1=Tools.mid(errorMult1, 4, low*1.6f-3);
//		final float em1=4;
		}else{throw new RuntimeException(""+errorPath);}
		
		return (low*em1<high || (low<=errorLowerConst && high>=Tools.max(minCountCorrect, low*errorMult2)));
	}
	
	protected final boolean isError(final int high, final int low, final byte q){
		float em1;
		if(errorPath==1){
			em1=errorMult1*(1+q*errorMultQFactor);
		}else if(errorPath==2){
			if(low<minCountCorrect && high>=minCountCorrect && q<20 && high>2*low){return true;}
			//		final float em1=errorMult1*(1+q*errorMultQFactor);
			em1=Tools.mid(errorMult1, 4, low*(q<=10 ? 1.6f : 2f)-3); //final float em1=4;//(low<3 ? 4 : Tools.min(errorMult1, 2*low));
			em1=em1*(1+q*errorMultQFactor);
		}else{throw new RuntimeException(""+errorPath);}
		
		return (low*em1<high || (low<=errorLowerConst && high>=Tools.max(minCountCorrect, low*errorMult2)));
	}
	
	public int extendToRight2(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int distance, boolean includeJunctionBase){
		if(verbose || verbose2){outstream.println("Entering extendToRight2 (no kmers).");}
		final int initialLength=bb.length();
		if(initialLength<k){return 0;}
		final int k2=k-1;
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0;
		long rkmer=0;
		
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts, to get the rightmost kmer */
		{
			int len=0;
			final byte[] bases=bb.array;
			for(int i=initialLength-k; i<initialLength; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{len++;}
				if(verbose){outstream.println("B: Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			}
			if(len<k){
				if(verbose || verbose2){outstream.println("Returning because len<k: "+len+"<"+k);}
				return 0;
			}
			else{assert(len==k);}
		}
		return extendToRight2(bb, leftCounts, rightCounts, distance, includeJunctionBase, kmer, rkmer);
	}
	
	/**
	 * Extend these bases to the right by at most 'distance'.
	 * Stops at right junctions only.
	 * Does not claim ownership.
	 */
	public int extendToRight2(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int distance, boolean includeJunctionBase,
			long kmer, long rkmer){
		if(verbose || verbose2){outstream.println("Entering extendToRight2 (with kmers).");}
		final int initialLength=bb.length();
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		
		/* Now the trailing kmer has been initialized. */
		
		long key=toValue(kmer, rkmer);
		int count=getCount(key);
		if(count<minCountCorrect){
			if(verbose || verbose2){outstream.println("Returning because count was too low: "+count+"<"+minCountCorrect);}
			return 0;
		}
		
		int leftMaxPos=0;
		int leftMax=minCountExtend;
		int leftSecondPos=1;
		int leftSecond=0;
		
		if(leftCounts!=null){
			leftMaxPos=fillLeftCounts(kmer, rkmer, leftCounts, mask, shift2);
			leftMax=leftCounts[leftMaxPos];
			leftSecondPos=Tools.secondHighestPosition(leftCounts);
			leftSecond=leftCounts[leftSecondPos];
		}
		
		int rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
		int rightMax=rightCounts[rightMaxPos];
		int rightSecondPos=Tools.secondHighestPosition(rightCounts);
		int rightSecond=rightCounts[rightSecondPos];
		
		if(verbose){
			outstream.println("kmer: "+toText(kmer)+", "+toText(rkmer));
			outstream.println("Counts: "+count+", "+Arrays.toString(rightCounts));
			outstream.println("rightMaxPos="+rightMaxPos);
			outstream.println("rightMax="+rightMax);
			outstream.println("rightSecondPos="+rightSecondPos);
			outstream.println("rightSecond="+rightSecond);
		}
		
		if(rightMax<minCountExtend){
			if(verbose || verbose2){outstream.println("Returning because rightMax was too low: "+rightMax+"<"+minCountExtend+"\n"+count+", "+Arrays.toString(rightCounts));}
			return 0;
		}
		if(isJunction(rightMax, rightSecond, leftMax, leftSecond)){
			if(verbose || verbose2){outstream.println("Returning because isJunction: "+rightMax+", "+rightSecond+"; "+leftMax+", "+leftSecond);}
			return 0;
		}
		
		final int maxLen=bb.length()+distance;
		
		while(bb.length()<maxLen){
			
			//Generate the new kmer
			final byte b=AminoAcid.numberToBase[rightMaxPos];
			final long x=rightMaxPos;
			final long x2=AminoAcid.numberToComplement[(int)x];
			
			final long evicted=(kmer>>>shift2); //Binary value that falls off the end.
			
			//Now consider the next kmer
			kmer=((kmer<<2)|(long)x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			
			key=toValue(kmer, rkmer);
			
			assert(getCount(key)==rightMax || rightMax==0);
			count=rightMax;
			
			assert(count>=minCountExtend) : count;
			
			if(leftCounts!=null){
				leftMaxPos=fillLeftCounts(kmer, rkmer, leftCounts, mask, shift2);
				leftMax=leftCounts[leftMaxPos];
				leftSecondPos=Tools.secondHighestPosition(leftCounts);
				leftSecond=leftCounts[leftSecondPos];
			}
			
			rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts, mask, shift2);
			rightMax=rightCounts[rightMaxPos];
			rightSecondPos=Tools.secondHighestPosition(rightCounts);
			rightSecond=rightCounts[rightSecondPos];
			
			if(verbose){
				outstream.println("kmer: "+toText(kmer)+", "+toText(rkmer));
				outstream.println("Counts: "+count+", "+Arrays.toString(rightCounts));
				outstream.println("rightMaxPos="+rightMaxPos);
				outstream.println("rightMax="+rightMax);
				outstream.println("rightSecondPos="+rightSecondPos);
				outstream.println("rightSecond="+rightSecond);
			}

			if(isJunction(rightMax, rightSecond, leftMax, leftSecond)){
				if(verbose){outstream.println("B: Breaking because isJunction("+rightMax+", "+rightSecond+", "+leftMax+", "+leftSecond+")");}
				if(includeJunctionBase && kmer>rkmer){
					bb.append(b);
					if(verbose){outstream.println("Added base "+(char)b);}
				}
				break;
			}
			
			if(leftCounts!=null && leftMaxPos!=evicted){
				if(verbose){outstream.println("B: Breaking because of hidden branch: leftMaxPos!=evicted ("+leftMaxPos+"!="+evicted+")" +
						"\nleftMaxPos="+leftMaxPos+", leftMax="+leftMax+", leftSecondPos="+leftSecondPos+", leftSecond="+leftSecond);}
				if(includeJunctionBase && kmer>rkmer){
					bb.append(b);
					if(verbose){outstream.println("Added base "+(char)b);}
				}
				break;
			}
			
			bb.append(b);
			if(verbose){outstream.println("Added base "+(char)b);}
			
			if(rightMax<minCountExtend){
				if(verbose || verbose2){outstream.println("C: Breaking because highest right was too low: "+rightMax+"<"+minCountExtend);}
				break;
			}
		}
		if(verbose || verbose2){System.err.println("Extended by "+(bb.length()-initialLength));}
		return bb.length()-initialLength;
	}
	
	public void regenerateKmers(byte[] bases, LongList kmers, IntList counts, final int a){
		final int loc=a+k;
		final int lim=Tools.min(counts.size, a+k+1);
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=kmers.get(a);
		long rkmer=rcomp(kmer);
		int len=k;
		
//		assert(false) : a+", "+loc+", "+lim;
		
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=loc, j=a+1; j<lim; i++, j++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			final long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{len++;}
			
			if(len>=k){
				assert(kmers.get(j)!=kmer);
				kmers.set(j, kmer);
				int count=getCount(kmer, rkmer);
				counts.set(j, count);
			}else{
				kmers.set(j, -1);
				counts.set(j, 0);
			}
		}
	}
	
	public int maxLeftCount(long kmer){
		long rkmer=rcomp(kmer);
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		rkmer=(rkmer<<2)&mask;
		kmer=(kmer>>>2);
		int max=-1, maxPos=0;
//		assert(false) : shift2+", "+k;
		for(int i=0; i<=3; i++){
			long rkmer2=rkmer|((long)AminoAcid.numberToComplement[i]);
			long kmer2=kmer|(((long)i)<<shift2);
			assert(kmer2==(kmer2&mask));
			assert(rkmer2==(rkmer2&mask));
			long key=toValue(rkmer2, kmer2);
			int count=getCount(key);
			count=Tools.max(count, 0);
			if(count>max){
				max=count;
				maxPos=i;
			}
		}
		return max;
	}
	
	public int maxRightCount(long kmer){
		long rkmer=rcomp(kmer);
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		if(verbose){outstream.println("fillRightCounts:   "+toText(kmer)+",   "+toText(rkmer));}
		kmer=(kmer<<2)&mask;
		rkmer=(rkmer>>>2);
		int max=-1, maxPos=0;
		
		for(int i=0; i<=3; i++){
			long kmer2=kmer|((long)i);
			long rkmer2=rkmer|(((long)AminoAcid.numberToComplement[i])<<shift2);
//			if(verbose){outstream.println("kmer:               "+toText(kmer2)+", "+toText(rkmer2));}
			assert(kmer2==(kmer2&mask));
			assert(rkmer2==(rkmer2&mask));
			assert(kmer2==rcomp(rkmer2));
			long key=toValue(kmer2, rkmer2);
			int count=getCount(key);
			count=Tools.max(count, 0);
			if(count>max){
				max=count;
				maxPos=i;
			}
		}
		return max;
	}
	
	public int fillLeftCounts(long kmer, long rkmer, int[] counts, long mask, int shift2){
		assert(kmer==rcomp(rkmer));
//		if(verbose){outstream.println("fillLeftCounts:    "+toText(kmer)+",   "+toText(rkmer));}
		rkmer=(rkmer<<2)&mask;
		kmer=(kmer>>>2);
		int max=-1, maxPos=0;
//		assert(false) : shift2+", "+k;
		for(int i=0; i<=3; i++){
			long rkmer2=rkmer|((long)AminoAcid.numberToComplement[i]);
			long kmer2=kmer|(((long)i)<<shift2);
//			if(verbose){outstream.println("kmer:             "+toText(kmer2)+",     "+toText(rkmer2));}
			assert(kmer2==(kmer2&mask));
			assert(rkmer2==(rkmer2&mask));
//			assert(kmer2==rcomp(rkmer2)) : "\n"+"kmer:      \t"+toText(rcomp(rkmer2))+", "+toText(rcomp(kmer2));
			long key=toValue(rkmer2, kmer2);
			int count=getCount(key);
			count=Tools.max(count, 0);
			counts[i]=count;
			if(count>max){
				max=count;
				maxPos=i;
			}
		}
		return maxPos;
	}
	
	public int fillRightCounts(long kmer, long rkmer, int[] counts, long mask, int shift2){
		assert(kmer==rcomp(rkmer));
		if(verbose){outstream.println("fillRightCounts:   "+toText(kmer)+",   "+toText(rkmer));}
		kmer=(kmer<<2)&mask;
		rkmer=(rkmer>>>2);
		int max=-1, maxPos=0;
		
		for(int i=0; i<=3; i++){
			long kmer2=kmer|((long)i);
			long rkmer2=rkmer|(((long)AminoAcid.numberToComplement[i])<<shift2);
//			if(verbose){outstream.println("kmer:               "+toText(kmer2)+", "+toText(rkmer2));}
			assert(kmer2==(kmer2&mask));
			assert(rkmer2==(rkmer2&mask));
			assert(kmer2==rcomp(rkmer2));
			long key=toValue(kmer2, rkmer2);
			int count=getCount(key);
			count=Tools.max(count, 0);
			counts[i]=count;
			if(count>max){
				max=count;
				maxPos=i;
			}
		}
		return maxPos;
	}
	
	protected final boolean isJunction(int rightMax, int rightSecond, int leftMax, int leftSecond){
		if(isJunction(rightMax, rightSecond)){return true;}
		return isJunction(leftMax, leftSecond);
	}
	
	protected final boolean isJunction(int max, int second){
		if(second<1 || second*branchMult1<max || (second<=branchLowerConst && max>=Tools.max(minCountExtend, second*branchMult2))){
			return false;
		}
		if(verbose){outstream.println("Breaking because second-highest was too high:\n" +
				"max="+max+", second="+second+", branchMult1="+branchMult1+"\n" +
				"branchLowerConst="+branchLowerConst+", minCountExtend="+minCountExtend+", branchMult2="+branchMult2+"\n" +
				(second*branchMult1<max)+", "+(second<=branchLowerConst)+", "+(max>=Tools.max(minCountExtend, second*branchMult2)));}
		return true;
	}
	
	/** Returns number of valid kmers */
	public int fillKmers(byte[] bases, LongList kmers){
		final int blen=bases.length;
		if(blen<k){return 0;}
		final int min=k-1;
		final int shift=2*k;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0;
		int len=0;
		int valid=0;

		kmers.clear();

		/* Loop through the bases, maintaining a forward kmer via bitshifts */
		for(int i=0; i<blen; i++){
			final byte b=bases[i];
			assert(b>=0) : Arrays.toString(bases);
			final long x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x<0){
				len=0;
				kmer=0;
			}else{len++;}
			if(i>=min){
				if(len>=k){
					kmers.add(kmer);
					valid++;
				}else{
					kmers.add(-1);
				}
			}
		}
		return valid;
	}
	
	/** Examines kmer counts around the merge borders to ensure the merge was not chimeric */ 
	public boolean mergeOK(Read merged, int len1, int len2, LongList kmers, final int width, final int thresh, final long highMult){
		final int len=merged.length();
		final int overlap=Tools.min(len, len1+len2-len);
		final byte[] bases=merged.bases;
		if(len<len1+width+1 || len<len2+width+1 || len<k+width){return true;}
//		corrector.fillCounts(bases, counts);
		int valid=fillKmers(bases, kmers);
		final int a=len-len2-1; //Base to left of first boundary
		final int b=a+1; //Base to right of first boundary 
		final int c=len1-1;
		final int d=c+1;

		final int ak=a-k+1; //kmer to left of first boundary
		final int bk=b; //kmer to right of first boundary
		final int ck=c-k+1; //kmer to left of second boundary
		final int dk=d; //kmer to right of second boundary
		
		//This is faster since fewer counts are looked up
		if(ak-width>=0 && ak+width<len){
			int min=getCount2(kmers.get(ak));
			for(int i=ak-width+1; i<ak; i++){min=Tools.min(min, getCount2(kmers.get(i)));}
			int min2=getCount2(kmers.get(ak+1));
			for(int i=ak+2; i<=ak+width; i++){min2=Tools.min(min2, getCount2(kmers.get(i)));}
			assert(min>=0 && min2>=0);
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}
		if(ck-width>=0 && ck+width<len){
			int min=getCount2(kmers.get(ck));
			for(int i=ck-width+1; i<ck; i++){min=Tools.min(min, getCount2(kmers.get(i)));}
			int min2=getCount2(kmers.get(ck+1));
			for(int i=ck+2; i<=ck+width; i++){min2=Tools.min(min2, getCount2(kmers.get(i)));}
			assert(min>=0 && min2>=0);
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}
		
		if(bk-width>=0 && bk+width<len){
			int min=getCount2(kmers.get(bk));
			for(int i=bk+1; i<bk+width+1; i++){min=Tools.min(min, getCount2(kmers.get(i)));}
			int min2=getCount2(kmers.get(bk-1));
			for(int i=bk-width; i<bk-1; i++){min2=Tools.min(min2, getCount2(kmers.get(i)));}
			assert(min>=0 && min2>=0);
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}
		if(dk-width>=0 && dk+width<len){
			int min=getCount2(kmers.get(dk));
			for(int i=dk+1; i<dk+width+1; i++){min=Tools.min(min, getCount2(kmers.get(i)));}
			int min2=getCount2(kmers.get(dk-1));
			for(int i=dk-width; i<dk-1; i++){min2=Tools.min(min2, getCount2(kmers.get(i)));}
			assert(min>=0 && min2>=0);
			if(min>=thresh && min2<=1){return false;}
			if(min2>0 && min>min2*highMult){return false;}
		}
		return true;
	}
	
	public int regenerateCounts(byte[] bases, IntList counts, final int ca){
		final int b=ca+k-1;
		final int lim=Tools.min(bases.length, b+k+1);
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;
		int len=0;
		int valid=0;
		
		final int clen=counts.size;
		
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts.
		 * i is an index in the base array, j is an index in the count array. */
		for(int i=b-k+1; i<lim; i++){
			final byte base=bases[i];
			final long x=AminoAcid.baseToNumber[base];
			final long x2=AminoAcid.baseToComplementNumber[base];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{
				len++;
			}
			
			if(i>=b){
				if(len>=k){
					valid++;
					int count=getCount(kmer, rkmer);
					counts.set(i-k+1, count);
				}else{
					counts.set(i-k+1, 0);
				}
			}
		}

		assert((counts.size==0 && bases.length<k) || counts.size==bases.length-k+1) : bases.length+", "+k+", "+counts.size;
		assert(clen==counts.size);
		
		return valid;
	}

	private final StringBuilder toText(long kmer){return AbstractKmerTable.toText(kmer, k);}
	private final long rcomp(long kmer){return AminoAcid.reverseComplementBinaryFast(kmer, k);}
	public final int getCount(long kmer, long rkmer){
		return (k==ksmall ? filter.getCount(kmer, rkmer) : filter.getCountBig(kmer));
	}
	public final int getCount(long key){
		return (k==ksmall ? filter.getCount(key) : filter.getCountBig(key));
	}
	public final int getCount2(long kmer){
		return kmer<0 ? 0 : (k==ksmall ? filter.getCount(toValue(kmer, rcomp(kmer))) : filter.getCountBig(kmer));
	}
	public final long toValue(long kmer, long rkmer){
		long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		return value;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       ThreadLocal Temps      ----------------*/
	/*--------------------------------------------------------------*/
	
	protected final void initializeThreadLocals(){
		if(localLeftCounts.get()!=null){return;}
		localLeftCounts.set(new int[4]);
		localRightCounts.set(new int[4]);
		localLongList.set(new LongList());
		localIntList.set(new IntList());
		localIntList2.set(new IntList());
		localByteBuilder.set(new ByteBuilder());
		localByteBuilder2.set(new ByteBuilder());
		localBitSet.set(new BitSet(300));
		localKmer.set(new Kmer(k));
		localKmer2.set(new Kmer(k));
		localTracker.set(new ErrorTracker());
	}
	
	protected ThreadLocal<int[]> localLeftCounts=new ThreadLocal<int[]>();
	protected ThreadLocal<int[]> localRightCounts=new ThreadLocal<int[]>();
	protected ThreadLocal<LongList> localLongList=new ThreadLocal<LongList>();
	protected ThreadLocal<IntList> localIntList=new ThreadLocal<IntList>();
	protected ThreadLocal<IntList> localIntList2=new ThreadLocal<IntList>();
	protected ThreadLocal<ByteBuilder> localByteBuilder=new ThreadLocal<ByteBuilder>();
	protected ThreadLocal<ByteBuilder> localByteBuilder2=new ThreadLocal<ByteBuilder>();
	protected ThreadLocal<BitSet> localBitSet=new ThreadLocal<BitSet>();
	private ThreadLocal<Kmer> localKmer=new ThreadLocal<Kmer>();
	private ThreadLocal<Kmer> localKmer2=new ThreadLocal<Kmer>();
	protected ThreadLocal<ErrorTracker> localTracker=new ThreadLocal<ErrorTracker>();
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	protected boolean ECC_PINCER=false;
	protected boolean ECC_TAIL=false;
	protected boolean ECC_ALL=false;
	protected boolean ECC_REASSEMBLE=true;
	protected boolean ECC_AGGRESSIVE=false;
	protected boolean ECC_CONSERVATIVE=false;
	protected boolean ECC_ROLLBACK=true;
	protected boolean ECC_REQUIRE_BIDIRECTIONAL=true;
	
	/** Mark bases as bad if they are completely covered by kmers with a count below this */
	protected int MARK_BAD_BASES=0;
	/** Only mark bad bases that are adjacent to good bases */
	protected boolean MARK_DELTA_ONLY=true;
	/** Only mark bad bases in reads that appear to have errors */
	protected boolean MARK_ERROR_READS_ONLY=true;
	/** Assign this quality score to marked bases */
	protected byte MARK_QUALITY=0;
	
	/*--------------------------------------------------------------*/
	
	BloomFilter filter;
	
	int k=31;
	int ksmall=31;
	final boolean rcomp=true;

	int minCountExtend=2;
	float branchMult1=20;
	float branchMult2=3;
	int branchLowerConst=3;
	
	int errorPath=1;
	float errorMult1=16;
	float errorMult2=2.6f;
	float errorMultQFactor=0.002f;
	int errorLowerConst=4;//3 seems fine
	int minCountCorrect=3;//5 is more conservative...
	int minCountCorrect(){return minCountCorrect;}
	int pathSimilarityConstant=3;
	float pathSimilarityFraction=0.45f;//0.3
	int errorExtensionReassemble=3;//default 2; higher is more conservative
	int errorExtensionPincer=3;//default 5; higher is more conservative
	int errorExtensionTail=8;//default 9; higher is more conservative
	int deadZone=0;
	int windowLen=12;
	int windowCount=6;
	int windowQualSum=80;
	
	byte qIncreasePincer=8;
	byte qMinPincer=24;
	byte qMaxPincer=32;
	
	byte qIncreaseTail=4;
	byte qMinTail=20;
	byte qMaxTail=28;

	boolean verbose=false;
	boolean verbose2=false;
	boolean smooth=true;
	int smoothWidth=3;
	
	PrintStream outstream=System.err;
	
}
