package test.consensus;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;

import consensus.BaseGraph;
import consensus.BaseGraphHelper;
import consensus.ConsensusObject;
import parse.Parse;
import stream.Read;

/** Matched in-memory traversal workload, including complete output consumption.
 * Small defaults calibrate size; performance trials need at least ten CPU seconds.
 * Uses one thread and a fixed 1024-base graph with 9:1 allele mixtures, no indels
 * and no trimming. Measures traversal, allocation and checksum work together.
 * @author Brian Bushnell, Sucrose */
public final class ConsensusTraverseBenchmark{

	/** Accepts reads=traversal count and warmup=untimed traversal count. */
	public static void main(final String[] args){
		long reads=1000, warmup=100;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException(arg);}
			if(pair[0].equals("reads")){reads=Parse.parseKMG(pair[1]);}
			else if(pair[0].equals("warmup")){warmup=Parse.parseKMG(pair[1]);}
			else{throw new IllegalArgumentException(arg);}
		}
		if(reads<1 || warmup<0){throw new IllegalArgumentException("Positive reads and nonnegative warmup required.");}
		final long bases=Math.multiplyExact(reads, (long)LENGTH);
		final boolean oldQuality=BaseGraph.qualityByMS, oldTrimNs=ConsensusObject.trimNs;
		final float oldTrimDepth=ConsensusObject.trimDepthFraction, oldMaf=ConsensusObject.MAF_sub;
		try{
			BaseGraph.qualityByMS=true;
			ConsensusObject.trimNs=false;
			ConsensusObject.trimDepthFraction=0;
			ConsensusObject.MAF_sub=0.25f;
			final BaseGraph graph=makeGraph();
			verify(graph.traverse());
			final ThreadMXBean bean=ManagementFactory.getThreadMXBean();
			if(!bean.isCurrentThreadCpuTimeSupported()){throw new IllegalStateException("Current-thread CPU timing required.");}
			if(!bean.isThreadCpuTimeEnabled()){bean.setThreadCpuTimeEnabled(true);}
			sink=traverse(graph, warmup);
			final long cpuStart=bean.getCurrentThreadCpuTime(), start=System.nanoTime();
			final long sum=traverse(graph, reads);
			final long elapsed=System.nanoTime()-start, cpu=bean.getCurrentThreadCpuTime()-cpuStart;
			sink=sum;
			if(cpu<0){throw new AssertionError("Monotonic current-thread CPU clock required.");}
			System.out.println("reads\tbases\twall_ns\tcpu_ns\tchecksum");
			System.out.println(reads+"\t"+bases+"\t"+elapsed+"\t"+cpu+"\t"+sum);
		}finally{
			BaseGraph.qualityByMS=oldQuality;
			ConsensusObject.trimNs=oldTrimNs;
			ConsensusObject.trimDepthFraction=oldTrimDepth;
			ConsensusObject.MAF_sub=oldMaf;
		}
	}

	/** Replaces constructor pseudocounts with nine majority and one alternate observation. */
	private static BaseGraph makeGraph(){
		final byte[] original=new byte[LENGTH];
		for(int pos=0; pos<original.length; pos++){original[pos]=DNA[pos&3];}
		final BaseGraph graph=new BaseGraph("traversal", original, null, 0, 0);
		for(int pos=0; pos<LENGTH; pos++){
			final int major=(pos+pos/7)&3, minor=(major+1)&3;
			BaseGraphHelper.setCounts(graph, pos, count(0, major, minor), count(1, major, minor),
				count(2, major, minor), count(3, major, minor), 0);
		}
		return graph;
	}

	/** Returns the literal two-allele mixture, independent of reference identity. */
	private static int count(final int allele, final int major, final int minor){
		return allele==major ? 9 : allele==minor ? 1 : 0;
	}

	/** Checks the finite fixture before timing: (9-half of 1)/10 times39 rounds to33. */
	private static void verify(final Read read){
		if(read.bases.length!=LENGTH || read.quality==null || read.quality.length!=LENGTH){
			throw new AssertionError("Untrimmed fixture must emit1024 bases and numeric qualities.");
		}
		for(int pos=0; pos<LENGTH; pos++){
			if(read.bases[pos]!=DNA[(pos+pos/7)&3] || read.quality[pos]!=33){
				throw new AssertionError("Unexpected 9:1 consensus at position "+pos+": base="+read.bases[pos]+" quality="+read.quality[pos]);
			}
		}
	}

	/** Consumes every emitted base and quality in order; long checksum overflow is intentional. */
	private static long traverse(final BaseGraph graph, final long reads){
		assert(graph.ref.length==LENGTH) : "The fixed workload and reported base count require1024 reference positions.";
		long sum=1;
		for(long i=0; i<reads; i++){
			final Read read=graph.traverse();
			assert(read.bases.length==LENGTH && read.quality!=null && read.quality.length==LENGTH)
				: "No-indel/no-trim traversal must preserve every fixture position and quality.";
			for(int pos=0; pos<read.bases.length; pos++){
				sum=31*sum+(read.bases[pos]&255);
				sum=31*sum+(read.quality[pos]&255);
			}
		}
		return sum;
	}

	private static final int LENGTH=1024;
	private static final byte[] DNA={'A', 'C', 'G', 'T'};
	private static volatile long sink;
}
