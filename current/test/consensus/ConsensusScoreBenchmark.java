package test.consensus;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import consensus.BaseGraph;
import consensus.BaseGraphHelper;
import parse.Parse;
import stream.Read;

/** Matched in-memory scoring workload; output includes a consumed checksum. @author Shinobu */
public final class ConsensusScoreBenchmark{

	/** Sets reads and warmup explicitly; small runs calibrate size, not performance. */
	public static void main(final String[] args){
		long reads=1000000, warmup=100000;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException(arg);}
			if(pair[0].equals("reads")){reads=Parse.parseKMG(pair[1]);}
			else if(pair[0].equals("warmup")){warmup=Parse.parseKMG(pair[1]);}
			else{throw new IllegalArgumentException(arg);}
		}
		if(reads<1 || warmup<0){throw new IllegalArgumentException("Positive reads and nonnegative warmup required.");}
		final BaseGraph[] models=new BaseGraph[8];
		final Read[] queries=new Read[32];
		final byte[] alphabet="ACGT".getBytes(StandardCharsets.US_ASCII);
		for(int g=0; g<models.length; g++){
			final byte[] original=new byte[80];
			Arrays.fill(original, (byte)'A');
			final BaseGraph graph=new BaseGraph("model"+g, original, null, g, 0);
			for(int pos=0; pos<original.length; pos++){
				final int allele=(pos+g+(pos/7))&3;
				BaseGraphHelper.setCounts(graph, pos, allele==0 ? 9 : 1,
					allele==1 ? 9 : 1, allele==2 ? 9 : 1, allele==3 ? 9 : 1, 0);
			}
			BaseGraphHelper.initForScoring(graph);
			models[g]=graph;
		}
		for(int q=0; q<queries.length; q++){
			final byte[] bases=new byte[80], match=new byte[80];
			for(int pos=0; pos<bases.length; pos++){bases[pos]=alphabet[(pos+q+pos/(q+1))&3];}
			Arrays.fill(match, (byte)'m');
			final Read read=new Read(bases, null, "query"+q, q);
			read.match=match;
			read.start=0;
			read.stop=79;
			read.setMapped(true);
			read.setShortMatch(false);
			queries[q]=read;
		}
		sink=score(models, queries, warmup);
		final ThreadMXBean bean=ManagementFactory.getThreadMXBean();
		if(!bean.isCurrentThreadCpuTimeSupported()){throw new IllegalStateException("Current-thread CPU timing required.");}
		if(!bean.isThreadCpuTimeEnabled()){bean.setThreadCpuTimeEnabled(true);}
		final long cpuStart=bean.getCurrentThreadCpuTime(), start=System.nanoTime();
		final double sum=score(models, queries, reads);
		final long elapsed=System.nanoTime()-start, cpu=bean.getCurrentThreadCpuTime()-cpuStart;
		sink=sum;
		if(!Double.isFinite(sum) || cpu<0){throw new AssertionError("Finite score checksum and monotonic CPU clock required.");}
		System.out.println("reads\tbases\twall_ns\tcpu_ns\tchecksum");
		System.out.println(reads+"\t"+Math.multiplyExact(reads, 80L)+"\t"+elapsed+"\t"+cpu+"\t"+sum);
	}

	/** Rotate models and queries independently so the loop consumes different scores. */
	private static double score(final BaseGraph[] models, final Read[] queries, final long reads){
		assert(models.length==8 && queries.length==32) : "Bit masks below match the explicit model/query fixture sizes.";
		double sum=0;
		for(long i=0; i<reads; i++){
			sum+=models[(int)i&7].score(queries[(int)(i>>>3)&31], false, true);
		}
		return sum;
	}

	private static volatile double sink;
}
