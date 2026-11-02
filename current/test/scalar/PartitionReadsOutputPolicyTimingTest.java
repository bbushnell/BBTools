package test.scalar;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.net.URI;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.DigestInputStream;
import java.security.MessageDigest;
import java.util.Locale;
import java.util.zip.GZIPInputStream;

import com.sun.management.OperatingSystemMXBean;

import shared.Timer;

/** Deterministic PartitionReads raw/gzip fan-out timing and parity harness. */
public class PartitionReadsOutputPolicyTimingTest{

	/** Runs one variant, compression mode, fan-out, input size and repeat count. */
	public static void main(String[] args) throws Exception{
		final String variant=args[0], compression=args[1].toLowerCase(Locale.ROOT);
		final Path root=Paths.get(args[2]);
		final int ways=Integer.parseInt(args[3]);
		final int threads=Integer.parseInt(args[4]);
		final long reads=Long.parseLong(args[5]);
		final int length=Integer.parseInt(args[6]);
		final int rounds=(args.length>7 ? Integer.parseInt(args[7]) : 1);
		final String writerMode=(args.length>8 ? args[8].toLowerCase(Locale.ROOT) : "auto");
		require(ways>0 && threads>0 && reads>0 && length>0 && rounds>0,
				"ways, threads, reads, length, and rounds must be positive");
		require(compression.equals("raw") || compression.equals("gzip"), "compression must be raw or gzip");
		require(writerMode.equals("auto") || writerMode.equals("st") || writerMode.equals("zt"),
				"writer mode must be auto, st, or zt");
		requireEstimatedSize(reads, length, rounds);

		Files.createDirectories(root);
		final Path realRoot=root.toRealPath();
		final URI sourceUri=PartitionReadsOutputPolicyTimingTest.class
				.getProtectionDomain().getCodeSource().getLocation().toURI();
		final Path sourceRoot=Paths.get(sourceUri).toRealPath();
		final Path checkoutRoot=sourceRoot.getParent();
		require(checkoutRoot!=null, "could not find BBTools checkout root from "+sourceRoot);
		require(!realRoot.startsWith(checkoutRoot),
				"temporary timing directory must be outside "+checkoutRoot);

		final Path input=root.resolve("reads_"+reads+"x"+length+".fq");
		writeFastq(input, reads, length);
		final Count inputCount=new Count(reads, reads*length);
		final String extension=compression.equals("gzip") ? ".fq.gz" : ".fq";
		final OperatingSystemMXBean os=(OperatingSystemMXBean)ManagementFactory.getOperatingSystemMXBean();

		final ThreadSampler sampler=new ThreadSampler();
		Count outputCount=new Count(0, 0);
		long outputBytes=0, validationNanos=0, validationCpuNanos=0;
		final StringBuilder decodedSha80=new StringBuilder();
		long processNanos=0, processCpuNanos=0;
		final int sampledJavaThreadPeak;
		try{
			for(int i=0; i<rounds; i++){
				final Path roundDir=root.resolve(compression+"_ways"+ways);
				Files.createDirectories(roundDir);
				final String pattern=roundDir.resolve("part%"+extension).toString();
				final long cpu0=os.getProcessCpuTime();
				final long wall0=System.nanoTime();
				sampler.startSampling();
				try{
					final Object partition=construct(variant, new String[]{"in="+input,
							pattern, "mode=gc", "ways="+ways, "t="+threads,
							"multiwritermode="+writerMode, "pigz=f", "unpigz=f",
							"zl=1", "ow=t"});
					process(partition);
				}finally{
					processNanos+=System.nanoTime()-wall0;
					processCpuNanos+=os.getProcessCpuTime()-cpu0;
					sampler.stopSampling();
				}
				final long validateCpu0=os.getProcessCpuTime();
				final long validateWall0=System.nanoTime();
				final Count current=countOutputs(roundDir, extension, ways, reads, length, i, decodedSha80);
				validationNanos+=System.nanoTime()-validateWall0;
				validationCpuNanos+=os.getProcessCpuTime()-validateCpu0;
				outputCount=new Count(outputCount.reads+current.reads, outputCount.bases+current.bases);
				outputBytes+=sumOutputs(roundDir, extension, ways);
			}
		}finally{
			sampledJavaThreadPeak=sampler.finishAndPeak();
		}

		require(outputCount.reads==inputCount.reads*rounds, "output read count changed");
		require(outputCount.bases==inputCount.bases*rounds, "output base count changed");
		System.out.println("RESULT"
				+"\tvariant="+variant
				+"\tcompression="+compression
				+"\tways="+ways
				+"\tthreads="+threads
				+"\tmultiwritermode="+writerMode
				+"\treads="+reads
				+"\tlength="+length
				+"\trounds="+rounds
				+"\tinputReads="+inputCount.reads
				+"\tinputBases="+inputCount.bases
				+"\toutputReads="+outputCount.reads
				+"\toutputBases="+outputCount.bases
				+"\toutputBytes="+outputBytes
				+"\tdecodedSha80="+decodedSha80
				+"\tpigz=f"
				+"\tunpigz=f"
				+"\tziplevel=1"
				+"\treusedOutputPaths=true"
				+"\tsampledJavaThreadPeak="+sampledJavaThreadPeak
				+"\tprocessCpuSeconds="+seconds(processCpuNanos)
				+"\tprocessWallSeconds="+seconds(processNanos)
				+"\tvalidationCpuSeconds="+seconds(validationCpuNanos)
				+"\tvalidationWallSeconds="+seconds(validationNanos)
				+"\trepeatedConstructorOpenCloseIncluded="+(rounds>1));
	}

	/** Builds a deterministic FASTQ with a repeatable sweep over GC fractions. */
	private static void writeFastq(Path path, long reads, int length) throws Exception{
		try(BufferedWriter bw=Files.newBufferedWriter(path, StandardCharsets.UTF_8)){
			for(long i=0; i<reads; i++){
				bw.write("@read");
				bw.write(Long.toString(i));
				bw.newLine();
				final int gc=(int)(i%(length+1));
				for(int j=0; j<length; j++){
					bw.write(base(i, j, gc));
				}
				bw.newLine();
				bw.write('+');
				bw.newLine();
				for(int j=0; j<length; j++){bw.write('I');}
				bw.newLine();
			}
		}
	}

	/** Reflectively constructs the selected historical or current partitioner. */
	private static Object construct(String variant, String[] args) throws Exception{
		final Class<?> clazz=Class.forName(className(variant));
		final Constructor<?> c=clazz.getConstructor(String[].class);
		try{
			return c.newInstance((Object)args);
		}catch(InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof Error){throw (Error)cause;}
			if(cause instanceof RuntimeException){throw (RuntimeException)cause;}
			throw new RuntimeException(cause);
		}
	}

	/** Maps fixture variant numbers to production classes. */
	private static String className(String variant){
		if(variant.equals("1")){return "scalar.PartitionReads";}
		if(variant.equals("2")){return "scalar.PartitionReads2";}
		if(variant.equals("3")){return "scalar.PartitionReads3";}
		throw new IllegalArgumentException("Unknown PartitionReads variant: "+variant);
	}

	/** Invokes package-private processing without changing production visibility. */
	private static void process(Object target) throws Exception{
		final Method m=target.getClass().getDeclaredMethod("process", Timer.class);
		m.setAccessible(true);
		m.invoke(target, new Timer());
	}

	/** Counts and validates all expected output FASTQ shards for one round. */
	private static Count countOutputs(Path dir, String extension, int ways, long expectedReads,
			int length, int round, StringBuilder decodedSha80) throws Exception{
		Count sum=new Count(0, 0);
		require(expectedReads<=Integer.MAX_VALUE, "too many reads to track unique ids");
		final boolean[] seen=new boolean[(int)expectedReads];
		for(int i=0; i<ways; i++){
			final Path path=dir.resolve("part"+i+extension);
			final Count current=countFastq(path, seen, length, round, i, decodedSha80);
			sum=new Count(sum.reads+current.reads, sum.bases+current.bases);
		}
		for(int i=0; i<seen.length; i++){
			require(seen[i], "missing read"+i+" in "+dir);
		}
		return sum;
	}

	/** Totals compressed or raw bytes in all expected output shards for one round. */
	private static long sumOutputs(Path dir, String extension, int ways) throws Exception{
		long bytes=0;
		for(int i=0; i<ways; i++){
			bytes+=Files.size(dir.resolve("part"+i+extension));
		}
		return bytes;
	}

	/** Counts one raw or gzip FASTQ shard and validates every generated record. */
	private static Count countFastq(Path path, boolean[] seen, int length,
			int round, int part, StringBuilder decodedSha80) throws Exception{
		long reads=0, bases=0;
		final MessageDigest digest=MessageDigest.getInstance("SHA-256");
		try(BufferedReader br=new BufferedReader(new InputStreamReader(new DigestInputStream(open(path), digest),
				StandardCharsets.UTF_8))){
			String header;
			while((header=br.readLine())!=null){
				final String sequence=br.readLine(), plus=br.readLine(), quality=br.readLine();
				require(sequence!=null && plus!=null && quality!=null, "truncated FASTQ record in "+path);
				require(header.startsWith("@") && plus.startsWith("+"), "malformed FASTQ record in "+path);
				final int id=parseId(header, path);
				require(id<seen.length, "unexpected read"+id+" in "+path);
				require(!seen[id], "duplicate read"+id+" in "+path);
				require(sequence.length()==length, "sequence length changed for read"+id+" in "+path);
				require(quality.length()==length, "quality length changed for read"+id+" in "+path);
				final int gc=id%(length+1);
				for(int i=0; i<length; i++){
					if(sequence.charAt(i)!=base(id, i, gc)){fail("sequence changed for read"+id+" in "+path);}
					if(quality.charAt(i)!='I'){fail("quality changed for read"+id+" in "+path);}
				}
				seen[id]=true;
				reads++;
				bases+=sequence.length();
			}
		}
		if(decodedSha80.length()>0){decodedSha80.append(',');}
		decodedSha80.append("round").append(round).append("/part").append(part)
				.append('=').append(sha80(digest.digest()));
		return new Count(reads, bases);
	}

	/** Opens a raw or gzip shard based on its filename. */
	private static InputStream open(Path path) throws Exception{
		final InputStream in=Files.newInputStream(path);
		return path.toString().endsWith(".gz") ? new GZIPInputStream(in) : in;
	}

	/** Refuses peak fixture sizes likely to leave more than one GB of generated files. */
	private static void requireEstimatedSize(long reads, int length, int rounds){
		final long bytesPerRecord=2L*length+64;
		final long batches=2L;
		require(batches>0 && bytesPerRecord>0, "size estimate overflow");
		require(reads<900000000L/bytesPerRecord/batches,
				"estimated generated input/output exceeds the 1GB fixture budget");
		require(reads<=Integer.MAX_VALUE, "too many reads to track unique ids");
	}

	/** Calculates the deterministic base emitted for a read position. */
	private static char base(long read, int position, int gc){
		final boolean high=position<gc;
		return (char)((((read+position)&1)==0) ? (high ? 'G' : 'A') : (high ? 'C' : 'T'));
	}

	/** Parses @readN and bounds it to the generated input id range. */
	private static int parseId(String header, Path path){
		require(header.startsWith("@read"), "unexpected header in "+path+": "+header);
		int end=5;
		while(end<header.length() && Character.isDigit(header.charAt(end))){end++;}
		require(end>5, "header lacks numeric id in "+path+": "+header);
		require(end==header.length(), "header has trailing text in "+path+": "+header);
		final long id=Long.parseLong(header.substring(5, end));
		require(id<Integer.MAX_VALUE, "read id is too large in "+path+": "+header);
		return (int)id;
	}

	/** Returns the final 80 bits of a SHA-256 digest as hex. */
	private static String sha80(byte[] digest){
		final StringBuilder sb=new StringBuilder(20);
		for(int i=digest.length-10; i<digest.length; i++){
			sb.append(String.format(Locale.ROOT, "%02x", digest[i]&0xff));
		}
		return sb.toString();
	}

	/** Converts nanoseconds to a compact decimal second string. */
	private static String seconds(long nanos){
		return String.format(Locale.ROOT, "%.6f", nanos*1e-9);
	}

	/** These checks remain active under both -ea and -da. */
	private static void require(boolean valid, String message){
		if(!valid){fail(message);}
	}

	/** Throws a validation failure. */
	private static void fail(String message){
		throw new AssertionError(message);
	}

	/** Immutable read/base total. */
	private static class Count{
		Count(long reads_, long bases_){
			reads=reads_;
			bases=bases_;
		}

		final long reads;
		final long bases;
	}

	/** Samples live Java thread count while the partitioner is running. */
	private static class ThreadSampler extends Thread{
		ThreadSampler(){
			super("PartitionReadsOutputPolicyTimingTest thread sampler");
			setDaemon(true);
			start();
		}

		@Override
		public void run(){
			final ThreadMXBean threads=ManagementFactory.getThreadMXBean();
			while(running){
				if(sampling){
					sampledPeak=Math.max(sampledPeak, threads.getThreadCount());
				}
				try{Thread.sleep(1);}
				catch(InterruptedException e){/* continue until finish() */}
			}
		}

		void startSampling(){
			sampling=true;
		}

		void stopSampling(){
			sampling=false;
		}

		int finishAndPeak() throws InterruptedException{
			running=false;
			interrupt();
			join();
			return sampledPeak;
		}

		private volatile boolean sampling=false;
		private volatile boolean running=true;
		private int sampledPeak=0;
	}
}
