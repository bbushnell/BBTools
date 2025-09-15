package shared;

import java.io.File;
import java.util.Arrays;

import dna.Data;
import fileIO.ByteFile1;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import structures.LongList;

/**
 * @author Brian Bushnell
 * @date Oct 9, 2012
 *
 */
public class Primes {
	
	public static void main(String[] args){
		
		if(args.length==3){makePrimes(args);}
		else{
			

			System.out.println(primeAtLeast(100));
			System.out.println(primeAtLeast(1000));
			System.out.println(primeAtLeast(10000));
			System.out.println(primeAtLeast(100000));
			System.out.println(primeAtLeast(1000000));
			System.out.println(primeAtLeast(10000000));
			System.out.println(primeAtLeast(100000000));
			System.out.println(primeAtLeast(1000000000));
			System.out.println(primeAtLeast(10000000000L));
			System.out.println(primeAtLeast(100000000000L));
			System.out.println(primeAtLeast(1000000000000L));
			System.out.println(primeAtLeast(10000000000000L));
			System.out.println(primeAtLeast(100000000000000L));
			System.out.println(primeAtLeast(1000000000000000L));
			

			System.out.println(primeAtMost(100));
			System.out.println(primeAtMost(1000));
			System.out.println(primeAtMost(10000));
			System.out.println(primeAtMost(100000));
			System.out.println(primeAtMost(1000000));
			System.out.println(primeAtMost(10000000));
			System.out.println(primeAtMost(100000000));
			System.out.println(primeAtMost(1000000000));
			System.out.println(primeAtMost(10000000000L));
			System.out.println(primeAtMost(100000000000L));
			System.out.println(primeAtMost(1000000000000L));
			System.out.println(primeAtMost(10000000000000L));
			System.out.println(primeAtMost(100000000000000L));
			System.out.println(primeAtMost(1000000000000000L));
			
		}
		
	}
	
	
	public static void makePrimes(String[] args){
		
		
		String in=args[0];
		String out=args[1];
		double mult=Double.parseDouble(args[2]);
		assert(mult>=1);
		
		long next=1;
		
		if(!new File(in).exists()){throw new RuntimeException("File not found: "+in);}
		TextFile tf=new TextFile(in, true);
		TextStreamWriter tsw=new TextStreamWriter(out, true, false, false);
		tsw.start();
		
//		int cnt=0;
		
		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
//			cnt++;
//			if(cnt>10000){break;}
			long x=Long.parseLong(s.trim());
			
//			System.out.println("cnt="+cnt+", x="+x+", next="+next);
			
			if(x>=next){
				tsw.print(x+"\n");
				next=(long)(x*mult);
			}
		}
		tsw.poison();
		tf.close();
		
	}
	
	
	public static long primeAtLeast(long x){
		int loc=Arrays.binarySearch(primes, x);
		if(loc<0){
			loc=-(loc+1);
			assert(loc>=primes.length || primes[loc]>x) : x;
		}
		if(loc>=primes.length){//Out of bounds
			long d=(long)Math.pow(x, 0.4);
			long a=primeAtLeast(x/d);
			long b=primeAtLeast(x/a);
			long c=a*b;
			assert(c>=x && c<(x*9)/8) : x+", "+a+", "+b+", "+c+", "+d;
			return c;
		}
		while(primes[loc]<x){loc++;}
		return primes[loc];
		
//		for(long p : primes){
//			if(p>=x){return p;}
//		}
//		throw new RuntimeException("No primes big enough for "+x);
	}

	public static int primeAtMost(int x){
		return (int)primeAtMost((long)x);
	}
	
	public static long primeAtMost(final long x){
		final int loc0=Arrays.binarySearch(primes, x);
		int loc=loc0;
		if(loc<0){
			loc=-(loc+1);
			assert(loc>=primes.length || primes[loc]>x) : x;
		}
		assert(loc>=0) : loc+", "+x;
		if(loc>=primes.length){//Out of bounds
			final long d=(long)Math.pow(x, 0.4);
			final long a=primeAtMost(x/d);
			final long b=primeAtMost(x/a);
			final long c=a*b;
			assert(c<=x && c>(x*7)/8) : x+", "+a+", "+b+", "+c+", "+d;
			return c;
		}
		assert(loc>=0) : loc+", "+x;//-1, 0, -300302388, 2
		assert(x>=primes[0]) : loc0+", "+loc+", "+x+", "+primes[0];
		while(primes[loc]>x){loc--;}
		return primes[loc];
		
//		for(int i=primes.length-1; i>=0; i--){
//			if(primes[i]<=x){return primes[i];}
//		}
//		throw new RuntimeException("No primes small enough for "+x);
	}
	

//	/**
//	 * @return
//	 */
//	private static long[] fetchPrimes() {
//		String fname=Data.findPath("?primes.txt.gz");
//		
//		TextFile tf=new TextFile(fname, false);
//		String[] lines=tf.toStringLines();
//		long[] array=new long[lines.length];
//		for(int i=0; i<lines.length; i++){
//			array[i]=Long.parseLong(lines[i].trim());
//		}
//		return array;
//	}
	
	private static long[] fetchPrimes() {
		String fname=Data.findPath("?primes.txt.gz");
		
		ByteFile1 tf=new ByteFile1(fname, false);
		LongList list=new LongList();
		for(byte[] line=tf.nextLine(); line!=null; line=tf.nextLine()){
			list.add(Parse.parseLong(line));
		}
		return list.toArray();
	}
	
	public static final long[] primes=fetchPrimes();
	
}
