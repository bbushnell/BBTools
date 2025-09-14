package fun;

import java.io.BufferedReader;
import java.io.FileReader;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.TimeZone;

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
	
	static String timeString(long time){
		SimpleDateFormat sdf=new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		sdf.setTimeZone(TimeZone.getTimeZone("PST"));
		return sdf.format(new Date(time));
	}

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
	static long lastSize=-1;
	
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
	
	static final long min(long x, long y){return x<y ? x : y;}
	
	static final int LOWER_BITS=31;
	static final int MANTISSA_BITS=24;
	static final int EXP_BITS=LOWER_BITS-MANTISSA_BITS;
	static final int UPPER_BITS=64-MANTISSA_BITS;
	static final long LOWER_MASK=~((-1L)<<LOWER_BITS);
	static final long MANTISSA_MASK=~((-1L)<<MANTISSA_BITS);
	static final long compress(long raw) {
		if(raw<=MANTISSA_MASK){return raw;}
		int leading=Long.numberOfLeadingZeros(raw);
		int exp=UPPER_BITS-leading;
		assert(exp>=1);
		return (raw>>>exp)|(exp<<MANTISSA_BITS);
	}
	static final long decompress(long f) {
		if(f<=MANTISSA_MASK){return f;}
		int exp=(int)(f>>>MANTISSA_BITS);
		assert(exp>=1);
		return (f&MANTISSA_MASK)<<exp;
	}
	static final long combine(long time, long size) {
		return (time<<LOWER_BITS) | compress(size);
	}
	static final long getTime(long combined) {
		return combined>>>LOWER_BITS;
	}
	static final long getSize(long combined) {
		return decompress(combined&LOWER_MASK);
	}
	
	static class LongList{
		
		public LongList(int initial){
			assert(initial>0) : initial;
			array=new long[initial];
		}
		
		public final long get(int loc){
			return array[loc];
		}
		
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
		
		public void sort() {
			if(size>1){Arrays.parallelSort(array, 0, size);}
		}
		
		public long[] array;
		/** Highest occupied index plus 1, i.e., lowest unoccupied index */
		public int size=0;
		
	}
	
	static final int MAX_ARRAY_LEN=Integer.MAX_VALUE-20;

}
