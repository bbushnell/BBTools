package test.consensus;

import java.lang.invoke.MethodHandle;
import java.lang.invoke.MethodHandles;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import consensus.Lilypad;
import parse.Parse;
import stream.SamLine;

/** Actual-method benchmark for Lilypad mixed-pair edge updates. @author Jean */
public final class LilypadHandleMixedPairBenchmark{

	/** Accepts reads=mixed-pair count and warmup=untimed mixed-pair count. */
	public static void main(final String[] args) throws Throwable{
		long reads=1000000, warmup=100000;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException(arg);}
			if(pair[0].equals("reads")){reads=Parse.parseKMG(pair[1]);}
			else if(pair[0].equals("warmup")){warmup=Parse.parseKMG(pair[1]);}
			else{throw new IllegalArgumentException(arg);}
		}
		if(reads<1 || warmup<0){throw new IllegalArgumentException("Positive reads and nonnegative warmup required.");}
		final boolean rnameAsBytes=SamLine.RNAME_AS_BYTES;
		try{
			final Runner warmupRunner=new Runner();
			final SamLine[] records=makeRecords();
			sink=warmupRunner.run(records, warmup).checksum;
			final ThreadMXBean bean=ManagementFactory.getThreadMXBean();
			if(!bean.isCurrentThreadCpuTimeSupported()){throw new IllegalStateException("Current-thread CPU timing required.");}
			if(!bean.isThreadCpuTimeEnabled()){bean.setThreadCpuTimeEnabled(true);}
			final Runner runner=new Runner();
			final long cpuStart=bean.getCurrentThreadCpuTime(), wallStart=System.nanoTime();
			final Result result=runner.run(records, reads);
			final long wall=System.nanoTime()-wallStart, cpu=bean.getCurrentThreadCpuTime()-cpuStart;
			sink=result.checksum;
			if(cpu<0){throw new AssertionError("Monotonic CPU clock required: "+cpu);}
			if(result.count!=reads || result.weight!=expectedWeight(reads)){
				throw new AssertionError("Expected count="+reads+" weight="+expectedWeight(reads)+
						", observed count="+result.count+" weight="+result.weight);
			}
			System.out.println("reads\tbases\twall_ns\tcpu_ns\tchecksum");
			System.out.println(reads+"\t"+(10*reads)+"\t"+wall+"\t"+cpu+"\t"+result.checksum);
		}finally{
			SamLine.RNAME_AS_BYTES=rnameAsBytes;
		}
	}

	/** Builds deterministic left/right mixed-pair records with ten bases per read. */
	private static SamLine[] makeRecords(){
		final SamLine[] records=new SamLine[64];
		for(int i=0; i<records.length; i++){
			final boolean left=(i&1)==1;
			final SamLine line=new SamLine();
			line.flag=(left ? 83 : 99);
			line.pos=3+(i&31);
			line.mapq=60-(i&7);
			line.setCigar("10M");
			line.pnext=5+(i&15);
			line.setRnameS("source");
			line.setRnextS("dest"+(i&15));
			line.setSeq("AAAAAAAAAA".getBytes(StandardCharsets.US_ASCII));
			records[i]=line;
		}
		return records;
	}

	/** Fresh actual Lilypad and Contig for each timed sample. */
	private static final class Runner{
		Runner() throws Throwable{
			lilypad=makeLilypad();
			set(lilypad, "maxPairDist", 3000);
			set(lilypad, "sameStrandPairs", false);
			final LinkedHashMap<?, ?> refMap=(LinkedHashMap<?, ?>)field(lilypad, "refMap");
			source=refMap.get("source");
		}

		Result run(final SamLine[] records, final long reads) throws Throwable{
			for(long i=0; i<reads; i++){HANDLE.invoke(source, records[(int)i&63]);}
			return result(source);
		}

		final Object lilypad;
		final Object source;
	}

	private static Object makeLilypad() throws Throwable{
		final Path dir=Files.createTempDirectory("lilypad-handle");
		final Path ref=dir.resolve("ref.fa");
		final Path sam=dir.resolve("empty.sam");
		final Path out=dir.resolve("out.fa");
		write(ref, refText());
		write(sam, "@HD\tVN:1.6\tSO:unsorted\n"+
				"@SQ\tSN:source\tLN:1000\n@SQ\tSN:dest0\tLN:1000\n");
		final Lilypad lilypad=new Lilypad(new String[]{"in="+sam, "ref="+ref, "out="+out, "ow=t", "append=f"});
		LOAD_REFERENCE.invoke(lilypad);
		return lilypad;
	}

	private static String refText(){
		final StringBuilder sb=new StringBuilder();
		sb.append(">source\n").append(repeat('A', 1000)).append('\n');
		for(int i=0; i<16; i++){
			sb.append(">dest").append(i).append('\n').append(repeat('C', 1000)).append('\n');
		}
		return sb.toString();
	}

	private static Result result(final Object contig){
		final Result left=result((LinkedHashMap<?, ?>)field(contig, "leftEdgeMap"));
		final Result right=result((LinkedHashMap<?, ?>)field(contig, "rightEdgeMap"));
		return new Result(left.checksum+31*right.checksum, left.count+right.count, left.weight+right.weight);
	}

	private static Result result(final LinkedHashMap<?, ?> map){
		long sum=17;
		long count=0;
		long weightSum=0;
		for(final Entry<?, ?> entry : map.entrySet()){
			final Object edge=entry.getValue();
			final long same=(Long)field(edge, "sameStrandCount");
			final long dif=(Long)field(edge, "difStrandCount");
			final long weight=(Long)field(edge, "weight");
			sum=31*sum+entry.getKey().hashCode();
			sum=31*sum+(Long)field(edge, "distanceSum");
			sum=31*sum+weight;
			sum=31*sum+same;
			sum=31*sum+dif;
			sum=31*sum+(Long)field(edge, "badCount");
			count+=same+dif;
			weightSum+=weight;
		}
		return new Result(sum, count, weightSum);
	}

	private static long expectedWeight(final long reads){
		long weight=(reads/64)*3616;
		for(long i=0, lim=reads&63; i<lim; i++){weight+=60-(i&7);}
		return weight;
	}

	private static final class Result{
		Result(final long checksum_, final long count_, final long weight_){
			checksum=checksum_;
			count=count_;
			weight=weight_;
		}
		final long checksum;
		final long count;
		final long weight;
	}

	private static Object field(final Object target, final String name){
		try{
			final Field field=target.getClass().getDeclaredField(name);
			field.setAccessible(true);
			return field.get(target);
		}catch(ReflectiveOperationException e){throw new AssertionError("Cannot read "+name, e);}
	}

	private static void set(final Object target, final String name, final Object value){
		try{
			final Field field=target.getClass().getDeclaredField(name);
			field.setAccessible(true);
			field.set(target, value);
		}catch(ReflectiveOperationException e){throw new AssertionError("Cannot set "+name, e);}
	}

	private static MethodHandle handle(final Class<?> type, final String name, final Class<?>... parameters){
		try{
			final Method method=type.getDeclaredMethod(name, parameters);
			method.setAccessible(true);
			return MethodHandles.lookup().unreflect(method);
		}catch(ReflectiveOperationException e){throw new ExceptionInInitializerError(e);}
	}

	private static Class<?> type(final String suffix){
		try{return Class.forName("consensus.Lilypad$"+suffix);}
		catch(ClassNotFoundException e){throw new ExceptionInInitializerError(e);}
	}

	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII));
	}

	private static String repeat(final char c, final int length){
		final char[] chars=new char[length];
		for(int i=0; i<length; i++){chars[i]=c;}
		return new String(chars);
	}

	private static final Class<?> CONTIG=type("Contig");
	private static final MethodHandle HANDLE=handle(CONTIG, "handleMixedPair", SamLine.class);
	private static final MethodHandle LOAD_REFERENCE=handle(Lilypad.class, "loadReferenceCustom");
	private static volatile long sink;
}
