package cardinality;

import shared.Tools;
import java.lang.reflect.Field;

/**
 * Measure empirical P(0) — fraction of filled buckets with history bit = 0.
 * Reports at geometrically-spaced cardinalities for comparison with
 * MantissaCompare2 single-bucket simulation.
 */
public class HistoryP0Test {

	public static void main(String[] args) throws Exception {
		int buckets=2048, maxCard=2000000, numTrials=32;
		for(String arg : args){
			String[] ab=arg.split("=");
			if(ab[0].equals("buckets")){buckets=Integer.parseInt(ab[1]);}
			else if(ab[0].equals("maxcard")){maxCard=Integer.parseInt(ab[1]);}
			else if(ab[0].equals("trials")){numTrials=Integer.parseInt(ab[1]);}
		}

		// Geometrically-spaced sample points
		final int numPoints=60;
		final int[] samplePoints=new int[numPoints];
		for(int i=0; i<numPoints; i++){
			samplePoints[i]=(int)Math.round(Math.pow(10.0, 1.0+i*4.0/(numPoints-1)));
		}

		System.out.println("Card\tFilled\tHist0\tHist1\tP0\tAvgTier");

		Field regsF=TwinTailLogLog.class.getDeclaredField("regs");
		regsF.setAccessible(true);
		Field geF=TwinTailLogLog.class.getDeclaredField("globalExp");
		geF.setAccessible(true);

		// Average over multiple trials
		for(int sp=0; sp<numPoints; sp++){
			final int targetCard=samplePoints[sp];
			if(targetCard>maxCard){break;}
			long totalFilled=0, totalH0=0, totalH1=0;
			double totalTierSum=0;

			for(int trial=0; trial<numTrials; trial++){
				TwinTailLogLog ttll=new TwinTailLogLog(buckets, 31, trial*12345L+67890, 0);

				for(int card=1; card<=targetCard; card++){
					ttll.add(card+trial*1000000000L);
				}

				byte[] regs=(byte[])regsF.get(ttll);
				int globalExp=geF.getInt(ttll);

				// Count P(0) across filled buckets
				int filled=0, h0=0, h1=0;
				double tierSum=0;
				for(int i=0; i<buckets; i++){
					int b=regs[i]&0xFF;
					int localExp=(b>>>4)&0xF;
					if(b==0 && globalExp==0){continue;} // truly empty
					filled++;
					int absNlz=globalExp+localExp;
					tierSum+=absNlz;
					int ch=b&0xF;
					int histBit=((ch&1)|((ch>>>2)&1));
					if(histBit==0){h0++;}else{h1++;}
				}
				totalFilled+=filled;
				totalH0+=h0;
				totalH1+=h1;
				totalTierSum+=tierSum;
			}

			double p0=(totalFilled>0) ? (double)totalH0/totalFilled : 0;
			double avgTier=(totalFilled>0) ? totalTierSum/totalFilled : 0;
			System.out.printf("%d\t%d\t%d\t%d\t%.4f\t%.2f%n",
				targetCard, totalFilled/numTrials, totalH0/numTrials,
				totalH1/numTrials, p0, avgTier);
		}
	}
}
