package cardinality;

import shared.Tools;

/**
 * Compare LL6 raw Mean vs ULL8 corrected Mean, both WITHOUT CF tables,
 * to isolate the effect of the history bit corrections on variance.
 */
public class ULL8Test {

	public static void main(String[] args){
		final int buckets=2048;
		final int trials=10000;
		final int[] cardinalities={100, 500, 1000, 2048, 4096, 10000, 20000, 50000, 100000};

		System.out.println("TrueCard\tLL6_RawMeanAbs\tULL8_CorrMeanAbs\tImprovement");

		for(int trueCard : cardinalities){
			double sumAbsLL6=0, sumAbsULL8=0;

			for(int trial=0; trial<trials; trial++){
				long seed=Tools.hash64shift((long)trial*1000000L+trueCard);

				// Track LL6 buckets manually (just byte values)
				byte[] ll6_buckets=new byte[buckets];
				byte[] ull8_buckets=new byte[buckets];

				for(int i=0; i<trueCard; i++){
					long val=Tools.hash64shift(seed+i);
					long key=Tools.hash64shift(val);
					int nlz=Long.numberOfLeadingZeros(key);
					int bucket=(int)(key&(buckets-1));
					int nlzStored=Math.min(nlz+1, 63);

					// LL6: just max
					int oldLL6=ll6_buckets[bucket]&0xFF;
					if(nlzStored>oldLL6){ll6_buckets[bucket]=(byte)nlzStored;}

					// ULL8: max + history bits
					int oldULL8=ull8_buckets[bucket]&0xFF;
					int oldNlzS=(oldULL8>0) ? (oldULL8>>2) : 0;
					if(nlzStored>oldNlzS){
						// Promotion: carry history bits
						int k=nlzStored-oldNlzS;
						int hist=0;
						if(oldULL8>0){
							if(k==1){hist=0x01; if((oldULL8&0x01)!=0) hist|=0x02;}
							else if(k==2){hist=0x02;}
						}
						ull8_buckets[bucket]=(byte)((nlzStored<<2)|hist);
					}else if(nlzStored<oldNlzS){
						int diff=oldNlzS-nlzStored;
						if(diff==1){ull8_buckets[bucket]=(byte)(oldULL8|0x01);}
						else if(diff==2){ull8_buckets[bucket]=(byte)(oldULL8|0x02);}
					}
				}

				// Compute raw Mean for both (same formula, no CF)
				double ll6mean=rawMean(ll6_buckets, buckets);
				double ull8mean=correctedMean(ull8_buckets, buckets);

				double errLL6=(ll6mean>0) ? Math.abs(ll6mean-trueCard)/(double)trueCard : 1.0;
				double errULL8=(ull8mean>0) ? Math.abs(ull8mean-trueCard)/(double)trueCard : 1.0;

				sumAbsLL6+=errLL6;
				sumAbsULL8+=errULL8;
			}

			double meanAbsLL6=sumAbsLL6/trials;
			double meanAbsULL8=sumAbsULL8/trials;
			double improvement=(meanAbsLL6-meanAbsULL8)/meanAbsLL6*100;

			System.out.printf("%d\t%.6f\t%.6f\t%+.2f%%\n",
				trueCard, meanAbsLL6, meanAbsULL8, improvement);
		}
	}

	/** Raw Mean: 2*(count+B)/B * count^2 / sum(2^(-nlz)) — no CF */
	static double rawMean(byte[] maxArray, int B){
		double sum=0;
		int count=0;
		for(int i=0; i<B; i++){
			int stored=maxArray[i]&0xFF;
			if(stored==0) continue;
			count++;
			int absNlz=stored-1;
			sum+=Math.pow(2.0, -absNlz);
		}
		if(count==0) return 0;
		return 2.0*((count+(double)B)/B)*((double)count*count)/sum;
	}

	/** Corrected Mean using ULL history bits with self-normalization */
	static double correctedMean(byte[] maxArray, int B){
		double corrSum=0, plainSum=0;
		int count=0;
		for(int i=0; i<B; i++){
			int stored=maxArray[i]&0xFF;
			if(stored==0) continue;
			count++;
			int absNlz=((stored>>2)&0x3F)-1;
			int state=stored&0x3;
			double base=Math.pow(2.0, -absNlz);
			plainSum+=base;
			corrSum+=base*STATE_MULT[state];
		}
		if(count==0) return 0;
		// Self-normalize: Mean from corrSum, scaled by normFactor
		// result = (2*(count+B)/B * count^2 / corrSum) * (corrSum/plainSum)
		//        = 2*(count+B)/B * count^2 / plainSum
		// That's the UNCORRECTED Mean! Self-normalization cancels out.
		//
		// Instead: the corrections should reduce VARIANCE not shift the MEAN.
		// The right approach: use corrSum directly with the Mean formula,
		// and let the fact that corrSum has less variance than plainSum
		// produce a better estimate. No normalization.
		return 2.0*((count+(double)B)/B)*((double)count*count)/corrSum;
	}

	// Unnormalized multipliers — let the CF tables handle bias
	static final double[] STATE_MULT={
		4.7745,  // state 00
		1.7108,  // state 01
		1.3056,  // state 10
		0.4713   // state 11
	};
}
