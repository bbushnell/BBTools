package fun;

import java.io.BufferedReader;
import java.io.FileReader;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.TimeZone;

/**
 * Advanced file metadata processing utility with comprehensive statistical
 * analysis and time-based filtering. Processes pipe-delimited files containing
 * file metadata, extracting size and access time information to generate
 * detailed statistical reports including percentile distributions and
 * terabyte-scale summaries.
 *
 * @author Brian Bushnell
 */
public class Foo5 {

	public static void main(String[] args) throws Exception{
		final BufferedReader reader=new BufferedReader(new FileReader(args[0]));
		long sum=0;

		LongList list=new LongList(100000000);
		for(String line=reader.readLine(); line!=null; line=reader.readLine()) {
			final long v=processFast(line);
			if(v>=0) {
				sum+=lastSize;
				list.add(v);
			}
		}
		reader.close();
		
		final long tebi=1024L*1024L*1024L*1024L;
		final long tera=1000L*1000L*1000L*1000L;
		
		{
			final int[] idxArray=new int[] {10, 20, 30, 40, 50, 60, 70, 80, 90, 95};
			final int[] pairArray=idxArray.clone();
			for(int i=0; i<idxArray.length; i++) {pairArray[i]=(int)(idxArray[i]*.01*list.size);}
			long tsum=0;

			list.sort();
			for(int i=0, nextIdx=0, nextPair=pairArray[0]; i<list.size && nextIdx<idxArray.length; i++) {
				final long v=list.get(i);
				final long size=getSize(v), time=getTime(v);
				tsum+=size;
				if(i>=nextPair) {
					String s=idxArray[nextIdx]+" percent of files have not been accessed since: "+
							timeString(time*1000)+" ("+(tsum/tebi)+" tebibytes)";
					System.out.println(s);
					nextIdx++;
					if(nextIdx<idxArray.length) {nextPair=pairArray[nextIdx];}
				}
			}
		}
		
		for(int i=0; i<list.size; i++) {
			list.array[i]=getSize(list.array[i]);
		}
		list.sort();
		long mean=sum/list.size;
		long median=list.get((int)(list.size*0.5));
		System.out.println("total size: \t"+(sum/tera)+" TB \t("+sum+")"+"\t"+"("+((sum/tebi))+" tebibytes)");
		System.out.println("mean size:  \t"+mean+" bytes");
		System.out.println("P50 size:   \t"+median+" bytes");
		System.out.println("P80 size:   \t"+list.get((int)(list.size*0.8))+" bytes");
		System.out.println("P90 size:   \t"+list.get((int)(list.size*0.9))+" bytes");
		System.out.println("P95 size:   \t"+list.get((int)(list.size*0.95))+" bytes");
	}
	
	/**
	 * Converts Unix timestamp to formatted PST time string.
	 * @param time Unix timestamp in milliseconds
	 * @return Formatted date string in "yyyy-MM-dd HH:mm:ss" format (PST)
	 */
	static String timeString(long time){
		SimpleDateFormat sdf=new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		sdf.setTimeZone(TimeZone.getTimeZone("PST"));
		return sdf.format(new Date(time));
	}

	/**
	 * Fast parser for pipe-delimited file metadata lines. Extracts file size
	 * from 4th field and access time from 13th field, filtering for files
	 * marked as 'F' type in 7th field. Uses manual character parsing for
	 * optimal performance with large datasets.
	 * @param line Pipe-delimited metadata line to parse
	 * @return Combined time/size value, or -1 if file type is not 'F'
	 */
	static long processFast(String line) {
		final int delimiter='|';
		final int len=line.length();
		int a=0, b=0;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		assert(b>a) : "Missing term size: '"+new String(line)+"'";
		long size=parseLong(line, a, b);
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		assert(b>a) : "Missing term type: '"+new String(line)+"'";
		if(line.charAt(a)!='F') {return -1;}
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		b++;
		a=b;

		while(b<len && line.charAt(b)!=delimiter){b++;}
		assert(b>a) : "Missing term access time: '"+new String(line)+"'";
		long time=parseLong(line, a, b);
		b++;
		a=b;

		lastSize=size;
		return combine(time, size);
	}
	/** Most recently processed file size in bytes */
	static long lastSize=-1;
	
	/**
	 * Efficient integer parsing from substring without creating new String
	 * objects. Manually converts character digits to long value with support
	 * for negative numbers.
	 * @param array Source string containing numeric data
	 * @param a Start index (inclusive) of number to parse
	 * @param b End index (exclusive) of number to parse
	 * @return Parsed long value from specified substring range
	 */
	static long parseLong(String array, int a, int b){
		assert(b>a);
		long r=0;
		final byte z='0';
		long mult=1;
		if(array.charAt(a)=='-'){mult=-1; a++;}
		for(; a<b; a++){
			int x=(array.charAt(a)-z);
			assert(x<10 && x>=0) : x+" = "+array.charAt(a)+"\narray="+new String(array)+", start="+a+", stop="+b;
			r=(r*10)+x;
		}
		return r*mult;
	}
	
	/**
	 * Returns smaller of two long values.
	 * @param x First value to compare
	 * @param y Second value to compare
	 * @return Minimum of x and y
	 */
	static final long min(long x, long y){return x<y ? x : y;}
	
	/** Number of bits used for lower portion of combined time/size encoding */
	static final int LOWER_BITS=31;
	/** Number of bits used for compressed size mantissa */
	static final int MANTISSA_BITS=24;
	/** Number of bits used for size compression exponent */
	static final int EXP_BITS=LOWER_BITS-MANTISSA_BITS;
	/** Number of bits available for timestamp in combined encoding */
	static final int UPPER_BITS=64-MANTISSA_BITS;
	/** Bit mask for extracting lower 31 bits from combined value */
	static final long LOWER_MASK=~((-1L)<<LOWER_BITS);
	/** Bit mask for extracting 24-bit mantissa from compressed size */
	static final long MANTISSA_MASK=~((-1L)<<MANTISSA_BITS);
	/**
	 * Compresses file size using floating-point-like representation with
	 * 24-bit mantissa and 7-bit exponent. Values under 2^24 are stored
	 * uncompressed; larger values use mantissa/exponent encoding.
	 * @param raw Original file size in bytes
	 * @return Compressed representation fitting in 31 bits
	 */
	static final long compress(long raw) {
		if(raw<=MANTISSA_MASK){return raw;}
		int leading=Long.numberOfLeadingZeros(raw);
		int exp=UPPER_BITS-leading;
		assert(exp>=1);
		return (raw>>>exp)|(exp<<MANTISSA_BITS);
	}
	/**
	 * Decompresses size value from mantissa/exponent representation back to
	 * original byte count. Reverses the compression performed by compress().
	 * @param f Compressed size value from compress() method
	 * @return Decompressed file size in bytes (approximate for large values)
	 */
	static final long decompress(long f) {
		if(f<=MANTISSA_MASK){return f;}
		int exp=(int)(f>>>MANTISSA_BITS);
		assert(exp>=1);
		return (f&MANTISSA_MASK)<<exp;
	}
	/**
	 * Combines timestamp and compressed file size into single long value.
	 * Uses upper 33 bits for time and lower 31 bits for compressed size.
	 * @param time Unix timestamp in seconds
	 * @param size File size in bytes (will be compressed)
	 * @return Combined value encoding both time and size
	 */
	static final long combine(long time, long size) {
		return (time<<LOWER_BITS) | compress(size);
	}
	/**
	 * Extracts timestamp from combined time/size value.
	 * @param combined Value from combine() method
	 * @return Original timestamp in seconds
	 */
	static final long getTime(long combined) {
		return combined>>>LOWER_BITS;
	}
	/**
	 * Extracts and decompresses file size from combined time/size value.
	 * @param combined Value from combine() method
	 * @return Decompressed file size in bytes
	 */
	static final long getSize(long combined) {
		return decompress(combined&LOWER_MASK);
	}
	
	/**
	 * Resizable array implementation for long values with automatic capacity
	 * expansion and parallel sorting support. Optimized for large datasets
	 * with efficient memory management and overflow protection.
	 * @author Brian Bushnell
	 */
	static class LongList{
		
		/** Creates new LongList with specified initial capacity.
		 * @param initial Initial array capacity (must be positive) */
		public LongList(int initial){
			assert(initial>0) : initial;
			array=new long[initial];
		}
		
		/**
		 * Returns value at specified index.
		 * @param loc Array index to retrieve
		 * @return Long value at specified location
		 */
		public final long get(int loc){
			return array[loc];
		}
		
		/** Appends value to end of list, expanding capacity if necessary.
		 * @param x Long value to add to list */
		public final void add(long x){
			if(size>=array.length){
				resize(size*2L+1);
			}
			array[size]=x;
			size++;
		}
		
		private final void resize(final long size2){
			assert(size2>size) : size+", "+size2;
			final int size3=(int)min(MAX_ARRAY_LEN, size2);
			assert(size3>size) : "Overflow: "+size+", "+size2+" -> "+size3;
			array=Arrays.copyOf(array, size3);
		}
		
		/** Sorts occupied portion of array in ascending order using parallel
		 * sort algorithm for optimal performance on large datasets. */
		public void sort() {
			if(size>1){Arrays.parallelSort(array, 0, size);}
		}
		
		/** Internal storage array for long values */
		public long[] array;
		/** Highest occupied index plus 1, i.e., lowest unoccupied index */
		public int size=0;
		
	}
	
	/** Maximum safe array length to prevent JVM OutOfMemoryError */
	static final int MAX_ARRAY_LEN=Integer.MAX_VALUE-20;

}
