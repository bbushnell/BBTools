package shared;

import java.io.File;
import java.util.Arrays;

import dna.Data;
import fileIO.ByteFile1;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import parse.Parse;
import structures.LongList;

/**
 * Selects capacities from the sorted, sparse {@code ?primes.txt.gz} resource.
 * Table lookups select stored primes, not necessarily the nearest mathematical
 * prime. Beyond the table, the legacy fallback multiplies smaller selected values;
 * its result need not be prime. Callers must not mutate the public table.
 * Initialization reads the resource once; later lookups are read-only provided
 * callers respect that ownership constraint.
 *
 * @author Brian Bushnell
 * @date Oct 9, 2012
 */
public class Primes{

	/** Filters a prime file with three arguments, otherwise prints sample lookups.
	 * @param args Input path, output path and spacing multiplier, or any other count for the demo */
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

	/** Writes a progressively spaced subset of an ascending prime-number file.
	 * Selects the next value at least the previous selected value times the multiplier.
	 * Closes owned input and joins a successfully started writer before returning.
	 * Read, parse and output failures are propagated; cleanup failures are suppressed
	 * onto an earlier failure so its original cause remains visible.
	 * @param args Input path, overwriteable output path, and multiplier at least one */
	public static void makePrimes(String[] args){
		String in=args[0];
		String out=args[1];
		double mult=Double.parseDouble(args[2]);
		assert(mult>=1);

		long next=1;

		if(!new File(in).exists()){throw new RuntimeException("File not found: "+in);}
		TextFile tf=new TextFile(in, true);
		TextStreamWriter tsw=null;
		boolean writerStarted=false;
		Throwable failure=null;
		try{
			tsw=new TextStreamWriter(out, true, false, false);
			tsw.start();
			writerStarted=true;
			for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
				long x=Long.parseLong(s.trim());
				if(x>=next){
					tsw.print(x+"\n");
					next=(long)(x*mult);
				}
			}
		}catch(RuntimeException | Error e){
			failure=e;
			throw e;
		}finally{
			Throwable cleanupFailure=null;
			try{
				if(tf.close()){throw new RuntimeException("Error reading prime list from "+in);}
			}catch(RuntimeException | Error e){cleanupFailure=e;}
			try{
				// poisonAndWait requires a successfully started worker; otherwise it waits forever.
				if(writerStarted && tsw.poisonAndWait()){
					throw new RuntimeException("Error writing prime list to "+out);
				}
			}catch(RuntimeException | Error e){
				if(cleanupFailure==null){cleanupFailure=e;}
				else if(cleanupFailure!=e){cleanupFailure.addSuppressed(e);}
			}
			if(cleanupFailure!=null){
				if(failure!=null){
					if(failure!=cleanupFailure){failure.addSuppressed(cleanupFailure);}
				}else if(cleanupFailure instanceof RuntimeException){throw (RuntimeException)cleanupFailure;}
				else{throw (Error)cleanupFailure;}
			}
		}
	}

	/** Selects the first stored prime at least x, or uses the legacy product fallback.
	 * Values below the first stored prime select that first entry. The beyond-table
	 * fallback retains its first-factor choice and rounds the required second factor
	 * upward. It has no general primality guarantee and does not search alternative
	 * factor pairs when its selected product exceeds the long range.
	 * @param x Requested minimum capacity
	 * @return Selected table entry or fallback product
	 * @throws ArithmeticException If the selected fallback product exceeds Long.MAX_VALUE */
	public static long primeAtLeast(long x){
		int loc=Arrays.binarySearch(primes, x);
		if(loc<0){
			loc=-(loc+1);
			assert(loc>=primes.length || primes[loc]>x) : x;
		}
		if(loc>=primes.length){//Out of bounds
			final long d=(long)Math.pow(x, 0.4);
			final long a=primeAtLeast(x/d);
			// Ceiling division prevents a table hit at floor(x/a) from undershooting x.
			final long b=primeAtLeast(x/a+(x%a==0 ? 0 : 1));
			if(a>Long.MAX_VALUE/b){
				throw new ArithmeticException("Prime capacity product exceeds long range: request="+x+", factors="+a+", "+b);
			}
			final long c=a*b;
			// Compare the excess instead of forming x*9: floor(9*x/8)=x+floor(x/8).
			assert(c>=x && c-x<x/8) : "Sparse upper capacity must satisfy x<=result<floor(9*x/8): "+x+", "+a+", "+b+", "+c+", "+d;
			return c;
		}
		while(primes[loc]<x){loc++;}
		return primes[loc];

//		for(long p : primes){
//			if(p>=x){return p;}
//		}
//		throw new RuntimeException("No primes big enough for "+x);
	}

	/** Int-capacity convenience overload; x must be at least the first table entry. */
	public static int primeAtMost(int x){
		return (int)primeAtMost((long)x);
	}

	/** Selects the last stored prime at most x, or uses the legacy product fallback.
	 * Requests below the first table entry are invalid. The beyond-table fallback
	 * has no general primality guarantee; its factors are rounded downward.
	 * @param x Requested maximum capacity, at least the first stored prime
	 * @return Selected table entry or fallback product */
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
			// floor(7*x/8)=x-ceil(x/8), without an overflowing intermediate product.
			final long lower=x-(x/8+(x%8==0 ? 0 : 1));
			assert(c<=x && c>lower) : "Sparse lower capacity must satisfy floor(7*x/8)<result<=x: "+x+", "+a+", "+b+", "+c+", "+d;
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

	/** Reads the shared resource in file order and checks final input status.
	 * Always closes acquired input; a parse/read failure remains primary when cleanup
	 * also fails. Assumes a nonempty sorted list; does not sort or validate primality. */
	private static long[] fetchPrimes(){
		String fname=Data.findPath("?primes.txt.gz");

		ByteFile1 tf=new ByteFile1(fname, false);
		Throwable failure=null;
		try{
			LongList list=new LongList();
			for(byte[] line=tf.nextLine(); line!=null; line=tf.nextLine()){
				list.add(Parse.parseLong(line));
			}
			return list.toArray();
		}catch(RuntimeException | Error e){
			failure=e;
			throw e;
		}finally{
			try{
				if(tf.close()){throw new RuntimeException("Error reading prime resource "+fname);}
			}catch(RuntimeException | Error e){
				if(failure==null){throw e;}
				if(failure!=e){failure.addSuppressed(e);}
			}
		}
	}

	/** Shared sparse ascending table, initialized once from the resource; do not mutate. */
	public static final long[] primes=fetchPrimes();

}
