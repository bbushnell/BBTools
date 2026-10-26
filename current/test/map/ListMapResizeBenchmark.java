package test.map;

import java.lang.management.ManagementFactory;
import java.util.ArrayList;
import java.util.Locale;
import map.LongArrayListHashMap;
import map.LongLongListHashMap;
import parse.Parse;
import structures.LongList;

/** Fixed-work growth benchmark with allocation counts for the two list-valued maps.
 * @author Shinobu
 */
public final class ListMapResizeBenchmark{

	/** Runs warmed single-threaded growth; allocation excludes warmup and reporting.
	 * @param args type=array|long, rounds=N, keys=N and values=N */
	public static void main(String[] args){
		String type="array";
		int rounds=512, keys=65536, values=8;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String a=arg.substring(0, eq).toLowerCase(Locale.ROOT), b=arg.substring(eq+1);
			if(a.equals("type")){type=b;}
			else if(a.equals("rounds")){rounds=Parse.parseIntKMG(b);}
			else if(a.equals("keys")){keys=Parse.parseIntKMG(b);}
			else if(a.equals("values")){values=Parse.parseIntKMG(b);}
			else{throw new IllegalArgumentException(arg);}
		}
		if(!type.equals("array") && !type.equals("long")){throw new IllegalArgumentException(type);}
		if(rounds<1 || keys<1 || values<1 || values>64){throw new IllegalArgumentException("Positive work sizes; values<=64");}
		run(type, 4, keys, values);
		final com.sun.management.OperatingSystemMXBean os=
			(com.sun.management.OperatingSystemMXBean)ManagementFactory.getOperatingSystemMXBean();
		final com.sun.management.ThreadMXBean threads=
			(com.sun.management.ThreadMXBean)ManagementFactory.getThreadMXBean();
		if(!threads.isThreadAllocatedMemorySupported()){throw new IllegalStateException("Allocation counter unavailable");}
		threads.setThreadAllocatedMemoryEnabled(true);
		final long tid=Thread.currentThread().getId(), allocated0=threads.getThreadAllocatedBytes(tid);
		final long cpu0=os.getProcessCpuTime(), start=System.nanoTime();
		final long sum=run(type, rounds, keys, values);
		final long elapsed=System.nanoTime()-start, cpu=os.getProcessCpuTime()-cpu0;
		final long allocated=threads.getThreadAllocatedBytes(tid)-allocated0;
		if(sum!=(long)rounds*keys*values*(values+1)/2){throw new AssertionError("Content checksum: "+sum);}
		System.out.printf(Locale.ROOT, "list-%s-growth\t%d\t%.6f\t%.6f\t%d%n",
			type, rounds, elapsed*1e-9, cpu*1e-9, sum);
		System.err.println("allocated_bytes="+allocated+" keys="+keys+" values="+values);
	}

	/** Dispatches outside the insertion loops and returns a content-size checksum. */
	private static long run(final String type, final int rounds, final int keys, final int values){
		assert(rounds>0 && keys>0 && values>0) : "Benchmark must execute positive fixed work";
		long sum=0;
		if(type.equals("array")){
			for(int r=0; r<rounds; r++){
				final LongArrayListHashMap<Long> map=new LongArrayListHashMap<Long>(4, 0.7f);
				for(int key=0; key<keys; key++){
					for(int value=0; value<values; value++){map.put(key, Long.valueOf(value));}
				}
				if(map.size()!=keys){throw new AssertionError("Missing array-list keys");}
				for(int key=0; key<keys; key++){
					final ArrayList<Long> list=map.get(key);
					if(list==null || list.size()!=values || list.get(0)!=0L || list.get(values-1)!=values-1L){
						throw new AssertionError("Array-list content oracle at "+key);
					}
					sum+=list.size();
					for(int j=0; j<values; j++){
						if(list.get(j)!=j){throw new AssertionError("Object-list value order");}
						sum+=list.get(j);
					}
				}
			}
		}else{
			for(int r=0; r<rounds; r++){
				final LongLongListHashMap map=new LongLongListHashMap(4, 0.7f);
				for(int key=0; key<keys; key++){
					for(int value=0; value<values; value++){map.put(key, (long)value);}
				}
				if(map.size()!=keys){throw new AssertionError("Missing long-list keys");}
				for(int key=0; key<keys; key++){
					final LongList list=map.get(key);
					if(list==null || list.size!=values || list.get(0)!=0 || list.get(values-1)!=values-1){
						throw new AssertionError("Long-list content oracle at "+key);
					}
					sum+=list.size;
					for(int j=0; j<values; j++){
						if(list.get(j)!=j){throw new AssertionError("Primitive-list value order");}
						sum+=list.get(j);
					}
				}
			}
		}
		return sum;
	}

	/** Static benchmark only. */
	private ListMapResizeBenchmark(){}
}
