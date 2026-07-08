package cardinality;

/**
 * Prices Brian's 1/3-bit: per family and size, the bump-counter value
 * histogram at end-of-stream plus the fraction pinned at saturation (c==3).
 * If LC families pile up at 3, a fifth counter state (the mixed-radix 1/3-bit
 * dividend) has a customer; if not, the bit belongs elsewhere.
 *
 * Usage: java cardinality.Fll53CtrHist [evalInst] [threads]
 *
 * @author Amber
 * @date July 2026
 */
public class Fll53CtrHist {

	public static void main(String[] args) throws Exception{
		final int evalInst=(args.length>0) ? Integer.parseInt(args[0]) : 64;
		final int threads=(args.length>1) ? Integer.parseInt(args[1])
			: Runtime.getRuntime().availableProcessors();
		CardinalityTracker.clampToAdded=false;

		System.out.println("#family\torgs\tc0%\tc1%\tc2%\tc3%(sat)\tfloorFull%");
		for(final int orgs : new int[]{192, 384, 768, 1536}){
			final int B=orgs*5;
			final ComplexityFarm.Spec[] evals={
				new ComplexityFarm.Spec(3, 244L*B, 4, 0),
				new ComplexityFarm.Spec(1, 98L*B, 1.0, 4),
				new ComplexityFarm.Spec(2, 1, 0.99, 196L*B),
				new ComplexityFarm.Spec(2, 100, 0.5, 390L*B),
				new ComplexityFarm.Spec(0, 146L*B, 1.0, 0),
				new ComplexityFarm.Spec(0, 98L*B, 2.0, 0),
			};
			final String[] names={"LC_minrand", "Zipf1.0", "onehot_k1",
				"onehot_k100", "unif_dup1", "unif_dup2"};
			for(int e=0; e<evals.length; e++){
				final ComplexityFarm.Spec sp=evals[e];
				final long[] hist=new long[4];
				final long[] full=new long[2];
				final Object lock=new Object();
				FllSkewTest.runPool(threads, evalInst, inst->{
					final Fll53[] sk=new Fll53[1];
					Fll53FarmW.runStream(sp, 909000L+inst*104729L, B,
						(f, distinct, adds)->{sk[0]=f;});
					final long[] h=new long[4];
					long fl=0;
					final Fll53 f=sk[0];
					for(int i=0; i<f.getNumOrgs(); i++){
						final int org=f.getOrg(i);
						h[org&Fll53.CTR_MASK]++;
						final int field=(org>>>Fll53.FIELD_SHIFT)&Fll53.FIELD_MASK;
						// floor bits nearly full = promotion imminent = antenna pressure
						if(Integer.bitCount(field&Fll53.LSB_MASK)>=4){fl++;}
					}
					synchronized(lock){
						for(int i=0; i<4; i++){hist[i]+=h[i];}
						full[0]+=fl; full[1]+=f.getNumOrgs();
					}
				});
				final double tot=hist[0]+hist[1]+hist[2]+hist[3];
				System.out.println(String.format("%s\t%d\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f",
					names[e], orgs, 100*hist[0]/tot, 100*hist[1]/tot,
					100*hist[2]/tot, 100*hist[3]/tot, 100.0*full[0]/full[1]));
			}
		}
	}
}
