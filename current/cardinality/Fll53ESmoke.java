package cardinality;

/** Smoke test for the Fll53E facade on referee-identical stream families.
 * Prints estimate vs truth at distinct-count doublings.  Not a benchmark. */
public class Fll53ESmoke {
	public static void main(String[] args){
		final String[] pats=(args.length>0) ? args[0].split(",")
			: new String[]{"fresh", "unifpool", "minrand"};
		final long maxTrue=(args.length>1) ? Long.parseLong(args[1]) : 2048L*1920;
		for(String pat : pats){
			final Fll53E e=new Fll53E(1024);
			final StreamSource src=StreamSource.make(pat, maxTrue, 12345L);
			final long maxDraws=src.maxDraws(maxTrue);
			System.out.println(pat+" stream, maxTrue="+maxTrue
				+", "+e.sizeBytes()+"B total:");
			long draws=0, nextCp=1024;
			double worstPastWarm=0;
			while(draws<maxDraws){
				e.add(src.next());
				draws++;
				if(src.distinct()>=nextCp){
					final long t=src.distinct();
					final long est=e.cardinality();
					final double err=100.0*(est-t)/t;
					System.out.printf("  adds=%d true=%d est=%d err=%+.2f%% cHat=%.3f%n",
						draws, t, est, err, e.complexity());
					if(t>=64L*1920){worstPastWarm=Math.max(worstPastWarm, Math.abs(err));}
					nextCp*=2;
				}
			}
			System.out.printf("  worst |err| past warmup: %.2f%%%n", worstPastWarm);
		}
	}
}
