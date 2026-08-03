package ddl;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import cardinality.DynamicDemiLog;
import shared.Timer;

/**
 * Histograms the frequency of DDL sketch bucket values across an entire
 * sketch set, binned by the top N bits of the 16-bit value.  Each bucket
 * stores a 16-bit value (0..65535); value 0 means an empty bucket and is
 * excluded.  Low values dominate the index when it holds many tiny/dense
 * sketches (like mitochondria), so this reveals how much of the index the
 * lowest values occupy and thus how much could be dropped to shrink it.
 *
 * Usage: DDLKeyHistogram <sketch.tsv.gz> [bits=6] [k=25] [t=32]
 *
 * @author Noire
 * @date August 3, 2026
 */
public class DDLKeyHistogram {

	public static void main(String[] args){
		String inFile=null;
		int k=25, threads=32, bits=6;
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("t") || a.equals("threads")){threads=Integer.parseInt(b);}
			else if(a.equals("bits")){bits=Integer.parseInt(b);}
			else if(a.equals("exponent") || a.equals("ebits")){DynamicDemiLog.setExponent(Integer.parseInt(b));}
			else if(inFile==null){inFile=arg;}
		}
		if(inFile==null){
			System.err.println("Usage: DDLKeyHistogram <sketch.tsv.gz> [bits=6] [k=25] [t=32]");
			System.exit(1);
		}
		assert(bits>=1 && bits<=16) : "bits must be in [1,16]: "+bits;

		Timer t=new Timer();
		System.err.println("Loading "+inFile+"...");
		final ArrayList<DDLRecord> records=DDLLoaderMT.loadFile(inFile, k, threads);
		System.err.println("Loaded "+records.size()+" records.");

		final int n=records.size();
		final int buckets=records.get(0).ddl.buckets;
		final int nBins=1<<bits;
		final int shift=16-bits;

		System.err.println("Histogramming "+n+" records ("+buckets+" buckets each) with "+threads+" threads...");

		final AtomicInteger nextQuery=new AtomicInteger(0);
		final long[][] threadHists=new long[threads][nBins];
		final long[] threadNonzero=new long[threads];

		Thread[] workers=new Thread[threads];
		for(int wi=0; wi<threads; wi++){
			final int wid=wi;
			workers[wi]=new Thread("KeyHist-"+wi){
				@Override
				public void run(){
					final long[] localHist=threadHists[wid];
					long localNonzero=0;
					while(true){
						final int qi=nextQuery.getAndIncrement();
						if(qi>=n){break;}
						final char[] max=records.get(qi).ddl.maxArray();
						for(int b=0; b<buckets; b++){
							final int v=max[b];
							if(v!=0){
								localHist[v>>>shift]++;
								localNonzero++;
							}
						}
					}
					threadNonzero[wid]=localNonzero;
				}
			};
			workers[wi].start();
		}
		for(Thread w : workers){
			try{w.join();}catch(InterruptedException e){}
		}

		final long[] hist=new long[nBins];
		long totalNonzero=0;
		for(int wi=0; wi<threads; wi++){
			for(int bin=0; bin<nBins; bin++){
				hist[bin]+=threadHists[wi][bin];
			}
			totalNonzero+=threadNonzero[wi];
		}

		t.stop();
		System.err.println("Done in "+t+".");

		final double denom=(totalNonzero>0 ? totalNonzero : 1);
		final StringBuilder sb=new StringBuilder();
		sb.append("Bin\tRangeLow\tRangeHigh\tCount\tFraction\tCumulative\n");
		long cumulative=0;
		for(int bin=0; bin<nBins; bin++){
			final int lo=bin<<shift;
			final int hi=((bin+1)<<shift)-1;
			cumulative+=hist[bin];
			sb.append(bin).append('\t').append(lo).append('\t').append(hi).append('\t')
				.append(hist[bin]).append('\t')
				.append(String.format("%.6f", hist[bin]/denom)).append('\t')
				.append(String.format("%.6f", cumulative/denom)).append('\n');
		}
		System.out.print(sb);

		System.err.println("Total records:       "+n);
		System.err.println("Buckets per record:  "+buckets);
		System.err.println("Total non-zero cells: "+totalNonzero);
	}
}
