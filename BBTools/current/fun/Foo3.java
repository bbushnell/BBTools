package fun;

import java.io.BufferedReader;
import java.io.FileReader;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.TimeZone;
import java.util.regex.Pattern;

public class Foo3 {

	public static void main(String[] args) throws Exception{
		final boolean slow=args.length<2 ? false : "slow".equalsIgnoreCase(args[1]);
		final BufferedReader reader=new BufferedReader(new FileReader(args[0]));
		long sum=0, lines=0, chars=0;

		LongList sizes=new LongList(100000000);
		ArrayList<Pair> list=new ArrayList<Pair>(10000000);
		for(String line=reader.readLine(); line!=null; line=reader.readLine()) {
			final Pair p=(slow ? processSlow(line) : processFast(line));
			if(p!=null) {
				sum+=p.size;
				lines++;
				chars+=line.length();
				list.add(p);
				sizes.add(p.size);
			}
		}
		
		reader.close();
		
		psort(list);
		final long tebi=1024L*1024L*1024L*1024L;
		final long tera=1000L*1000L*1000L*1000L;
		final int[] idxArray=new int[] {10, 20, 30, 40, 50, 60, 70, 80, 90, 95};
		final int[] pairArray=idxArray.clone();
		for(int i=0; i<idxArray.length; i++) {pairArray[i]=(int)(idxArray[i]*.01*list.size());}
		long tsum=0;
		for(int i=0, nextIdx=0, nextPair=pairArray[0]; i<list.size() && nextIdx<idxArray.length; i++) {
			Pair p=list.get(i);
			tsum+=p.size;
			if(i>=nextPair) {
				String s=idxArray[nextIdx]+" percent of files have not been accessed since: "+
						timeString(p.time*1000)+" ("+(tsum/tebi)+" tebibytes)";
				System.out.println(s);
				nextIdx++;
				if(nextIdx<idxArray.length) {nextPair=pairArray[nextIdx];}
			}
		}
		
		sizes.sort();
		long mean=sum/sizes.size;
		long median=sizes.get((int)(sizes.size*0.5));
		System.out.println("total size: \t"+(sum/tera)+" TB \t("+sum+")"+"\t"+"("+((sum/tebi))+" tebibytes)");
		System.out.println("mean size:  \t"+mean+" bytes");
		System.out.println("P50 size:   \t"+median+" bytes");
		System.out.println("P80 size:   \t"+sizes.get((int)(sizes.size*0.8))+" bytes");
		System.out.println("P90 size:   \t"+sizes.get((int)(sizes.size*0.9))+" bytes");
		System.out.println("P95 size:   \t"+sizes.get((int)(sizes.size*0.95))+" bytes");
	}
	
	static <T extends Comparable<? super T>> void psort(ArrayList<T> list) {
		@SuppressWarnings("unchecked")
		T[] array=list.toArray((T[])new Comparable[0]);
		list.clear();
		Arrays.parallelSort(array);
		for(T r : array){list.add(r);}
	}
	
	static String timeString(long time){
		SimpleDateFormat sdf=new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		sdf.setTimeZone(TimeZone.getTimeZone("PST"));
		return sdf.format(new Date(time));
	}
	
	/**Example line:<br>
	*
	120395044|1465795004|3191|11574|MEP_itags|1465795004|FAu|1|14297|14297
	|-rw-r-----|1496172205|1496172205|1048576|1516819694
	|%2Fdna%2F.snapshots%2F2023-05-17%2Fprojectdirs%2FMEP%2Fitags%2FiTagger3%2F
	data%2FDIVERSITY%2F1%2F1822.qiime%2Fbdiv_even1000%2Funweighted_unifrac_emperor
	_pcoa_plot%2Femperor_required_resources%2Fjs%2Fjquery.colorPicker.js
	*/
	static Pair processSlow(String line) {
		String[] split=pipePattern.split(line);
		if(split[6].charAt(0)!='F') {return null;}
		long size=Long.parseLong(split[3]);
		long time=Long.parseLong(split[11]);
		return new Pair(size, time);
	}

	static Pair processFast(String line) {
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
		if(line.charAt(a)!='F') {return null;}
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

		return new Pair(size, time);
	}
	
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
	
	public static final long min(long x, long y){return x<y ? x : y;}
	
	private static class Pair implements Comparable<Pair> {
		
		Pair(long size_, long time_){
			size=size_;
			time=time_;
		}
		
		@Override
		public int compareTo(Pair b) {
			return time>b.time ? 1 : time<b.time ? -1 : 0;
		}
		
		final long size, time;
	}
	
//	private static class TimeComparator implements Comparator<Pair> {
//		
//		private TimeComparator() {}
//		
//		@Override
//		public int compare(Pair a, Pair b) {
//			return a.time>b.time ? 1 : a.time<b.time ? -1 : 0;
//		}
//		
//		static final TimeComparator instance=new TimeComparator();
//		
//	}
	
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
	
	static final Pattern pipePattern=Pattern.compile("\\|");
	static final int MAX_ARRAY_LEN=Integer.MAX_VALUE-20;

}
