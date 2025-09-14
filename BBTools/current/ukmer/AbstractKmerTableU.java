package ukmer;

import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.Lock;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * @author Brian Bushnell
 * @date Oct 23, 2013
 *
 */
public abstract class AbstractKmerTableU {
	
	/*--------------------------------------------------------------*/
	/*----------------         Kmer methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns count */
	public abstract int increment(Kmer kmer);
	
	/** Returns number of entries created */
	public abstract int incrementAndReturnNumCreated(final Kmer kmer);

	public abstract int set(Kmer kmer, int value);
	
	public abstract int set(Kmer kmer, int[] vals);
	
	/** Returns number of kmers added */
	public abstract int setIfNotPresent(Kmer kmer, int value);

	/**
	 * Fetch the value associated with a kmer.
	 * @param kmer
	 * @return A value.  -1 means the kmer was not present.
	 */
	public abstract int getValue(Kmer kmer);
	
	/**
	 * Fetch the values associated with a kmer.
	 * @param kmer
	 * @param singleton A blank array of length 1.
	 * @return An array filled with values.  Values of -1 are invalid.
	 */
	public abstract int[] getValues(Kmer kmer, int[] singleton);

	public abstract boolean contains(Kmer kmer);
	
//	public abstract boolean contains(Kmer kmer, int v);
//
//	public abstract boolean contains(Kmer kmer, int[] vals);
//
//	public abstract Object get(Kmer kmer);
	
	public static final int compare(long[] key1, long[] key2){
		for(int i=0; i<key1.length; i++){
			long dif=key1[i]-key2[i];
			if(dif!=0){return (int)Tools.mid(-1, dif, 1);}
		}
		return 0;
	}
	
	public static final boolean equals(long[] key1, long[] key2){
		return compare(key1, key2)==0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Abstract Methods       ----------------*/
	/*--------------------------------------------------------------*/

	public abstract int getValue(long[] key, long xor);
	
//	/** Returns count */
//	public final int increment(long[] key){throw new RuntimeException();}
//
//	/** Returns number of entries created */
//	public final int incrementAndReturnNumCreated(final long[] key){throw new RuntimeException();}
//
//	public final int set(long[] key, int value){throw new RuntimeException();}
//
//	public final int set(long[] key, int[] vals){throw new RuntimeException();}
//
//	/** Returns number of kmers added */
//	public final int setIfNotPresent(long[] key, int value){throw new RuntimeException();}
//
//	/**
//	 * Fetch the value associated with a kmer.
//	 * @param kmer
//	 * @return A value.  -1 means the kmer was not present.
//	 */
//	final int getValue(long[] key){throw new RuntimeException();}
//
//	/**
//	 * Fetch the values associated with a kmer.
//	 * @param kmer
//	 * @param singleton A blank array of length 1.
//	 * @return An array filled with values.  Values of -1 are invalid.
//	 */
//	public final int[] getValues(long[] key, int[] singleton){throw new RuntimeException();}
//
//	public final boolean contains(long[] key){throw new RuntimeException();}
	
	public final boolean contains(Kmer kmer, int v){
		assert(TESTMODE);
		int[] set=getValues(kmer, new int[] {-1});
		if(set==null){return false;}
		for(int s : set){
			if(s==-1){break;}
			if(s==v){return true;}
		}
		return false;
	}
	
	public final boolean contains(Kmer kmer, int[] vals){
		assert(TESTMODE);
		int[] set=getValues(kmer, new int[] {-1});
		if(set==null){return false;}
		boolean success=true;
		for(int v : vals){
			if(v==-1){break;}
			success=false;
			for(int s : set){
				if(s==v){
					success=true;
					break;
				}
			}
			if(!success){break;}
		}
		return success;
	}

	public abstract void rebalance();

	public abstract long size();
	public abstract int arrayLength();
	public abstract boolean canRebalance();

	public abstract boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount);
	public abstract boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining);
	public abstract boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, final int maxcount, AtomicLong remaining);
	
	public abstract void fillHistogram(long[] ca, int max);
	public abstract void fillHistogram(SuperLongList sll);
	public abstract void countGC(long[] gcCounts, int max);
	
	public static final int gc(long kmer){
		int gc=0;
		while(kmer>0){
			long x=kmer&3;
			kmer>>>=2;
			if(x==1 || x==2){gc++;}
		}
		return gc;
	}
	
	Object get(Kmer kmer){return get(kmer.key());}
	abstract Object get(long[] key);
	abstract void resize();
	abstract boolean canResize();
	


	/**
	 * Removes entries with a value of zero or less.
	 * Rehashes the remainder.
	 * @return Number removed.
	 */
	abstract long regenerate(final int limit);

	final void lock(){getLock().lock();}
	final void unlock(){getLock().unlock();}
	final boolean tryLock(){return getLock().tryLock();}
	Lock getLock(){
		throw new RuntimeException("Unimplemented.");
	}
	
	/*--------------------------------------------------------------*/
	/*---------------       Allocation Methods      ----------------*/
	/*--------------------------------------------------------------*/
	
	final static AtomicIntegerArray allocAtomicInt(int len){
		return KillSwitch.allocAtomicInt(len);
	}
	
	final static long[] allocLong1D(int len){
		return KillSwitch.allocLong1D(len);
	}
	
	final static long[][] allocLong2D(int mult, int len){
		return KillSwitch.allocLong2D(mult, len);
	}
	
	final static int[] allocInt1D(int len){
		return KillSwitch.allocInt1D(len);
	}
	
	final static int[][] allocInt2D(int len){
		return KillSwitch.allocInt2D(len);
	}
	
	final static KmerNodeU[] allocKmerNodeArray(int len){
		KmerNodeU[] ret=null;
		try {
			ret=new KmerNodeU[len];
		} catch (OutOfMemoryError e) {
			synchronized(killMessage){
				e.printStackTrace();
				System.err.println(killMessage);
				KillSwitch.killSilent();
			}
		}
		return ret;
	}
	
	/*--------------------------------------------------------------*/
	/*---------------       Ownership Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Set the thread owning this kmer.  Return the new owner.
	 * Will only change the owner if newOwner is greater than current owner. */
	public abstract int setOwner(Kmer kmer, int newOwner);
	
	/** Reset owner to -1 if this is the current owner. */
	public abstract boolean clearOwner(Kmer kmer, int owner);
	
	/** Return the thread ID owning this kmer, or -1. */
	public abstract int getOwner(Kmer kmer);
	
	/** Create data structures needed for ownership representation */
	public abstract void initializeOwnership();
	
	/** Eliminate ownership data structures or set them to -1. */
	public abstract void clearOwnership();
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final StringBuilder toText(Kmer kmer){
		return toText(kmer.key(), kmer.k);
	}
	
	public static final StringBuilder toText(long[] array, int k){
		StringBuilder sb=new StringBuilder(k*array.length);
		for(int pos=0; pos<array.length; pos++){
			long kmer=array[pos];
			for(int i=k-1; i>=0; i--){
				int x=(int)((kmer>>(2*i))&3);
				sb.append((char)AminoAcid.numberToBase[x]);
			}
		}
		return sb;
	}

	static final StringBuilder toText(long[] array, int count, int k){
		StringBuilder sb=new StringBuilder(k+10);
		return toText(array, count, k, sb);
	}

	static final ByteBuilder toBytes(long[] array, int count, int k){
		ByteBuilder bb=new ByteBuilder(k+10);
		return toBytes(array, count, k, bb);
	}

	static final StringBuilder toText(long[] array, int[] values, int k){
		StringBuilder sb=new StringBuilder(k+10);
		return toText(array, values, k, sb);
	}

	static final ByteBuilder toBytes(long[] array, int[] values, int k){
		ByteBuilder bb=new ByteBuilder(k+10);
		return toBytes(array, values, k, bb);
	}
	
	static final StringBuilder toText(long[] array, int count, int k, StringBuilder sb){
		if(FASTA_DUMP){
			sb.append('>');
			sb.append(count);
			sb.append('\n');
			for(int i=0; i<array.length; i++){
				append(array[i], k, sb);
			}
		}else{
			for(int i=0; i<array.length; i++){
				append(array[i], k, sb);
			}
			sb.append('\t');
			sb.append(count);
		}
		return sb;
	}
	
	static final StringBuilder toText(long[] array, int[] values, int k, StringBuilder sb){
		if(FASTA_DUMP){
			sb.append('>');
			for(int i=0; i<values.length; i++){
				int x=values[i];
				if(x==-1){break;}
				if(i>0){sb.append(',');}
				sb.append(x);
			}
			sb.append('\n');
			for(int i=0; i<array.length; i++){
				append(array[i], k, sb);
			}
		}else{
			for(int i=0; i<array.length; i++){
				append(array[i], k, sb);
			}
			sb.append('\t');
			for(int i=0; i<values.length; i++){
				int x=values[i];
				if(x==-1){break;}
				if(i>0){sb.append(',');}
				sb.append(x);
			}
		}
		return sb;
	}
	
	private static final void append(long kmer, int k, StringBuilder sb){
		for(int i=k-1; i>=0; i--){
			int x=(int)((kmer>>(2*i))&3);
			sb.append((char)AminoAcid.numberToBase[x]);
		}
	}
	
	public static final ByteBuilder toBytes(long[] array, long count, int k, ByteBuilder sb){
		if(FASTA_DUMP){
			sb.append('>');
			sb.append(count);
			sb.append('\n');
			for(int i=0; i<array.length; i++){
				append(array[i], k, sb);
			}
		}else{
			for(int i=0; i<array.length; i++){
				append(array[i], k, sb);
			}
			sb.append('\t');
			sb.append(count);
		}
		return sb;
	}
	
	public static final ByteBuilder toBytes(long[] array, int[] values, int k, ByteBuilder sb){
		if(FASTA_DUMP){
			sb.append('>');
			for(int i=0; i<values.length; i++){
				int x=values[i];
				if(x==-1){break;}
				if(i>0){sb.append(',');}
				sb.append(x);
			}
			sb.append('\n');
			for(int i=0; i<array.length; i++){
				append(array[i], k, sb);
			}
		}else{
			for(int i=0; i<array.length; i++){
				append(array[i], k, sb);
			}
			sb.append('\t');
			for(int i=0; i<values.length; i++){
				int x=values[i];
				if(x==-1){break;}
				if(i>0){sb.append(',');}
				sb.append(x);
			}
		}
		return sb;
	}
	
	private static final void append(long kmer, int k, ByteBuilder sb){
		for(int i=k-1; i>=0; i--){
			int x=(int)((kmer>>(2*i))&3);
			sb.append((char)AminoAcid.numberToBase[x]);
		}
	}
	
	
//	static void appendKmerText(long kmer, int count, int k, StringBuilder sb){
//		sb.setLength(0);
//		toText(kmer, count, k, sb);
//		sb.append('\n');
//	}
	
	static void appendKmerText(long[] array, int count, int k, ByteBuilder bb){
		bb.setLength(0);
		toBytes(array, count, k, bb);
		bb.nl();
	}
	
	
	/** For buffered tables. */
	long flush(){
		throw new RuntimeException("Unsupported.");
	}
	
	/**
	 * This allocates the data structures in multiple threads.  Unfortunately, it does not lead to any speedup, at least for ARRAY type.
	 * @param ways
	 * @param tableType
	 * @param schedule
	 * @return Preallocated tables.
	 */
	public static final AbstractKmerTableU[] preallocate(int ways, int tableType, int[] schedule, int k, int kbig){

		final AbstractKmerTableU[] tables=new AbstractKmerTableU[ways];
		
		{
			final int t=Tools.max(1, Tools.min(Shared.threads(), 2, ways));
			final AllocThread[] allocators=new AllocThread[t];
			for(int i=0; i<t; i++){
				allocators[i]=new AllocThread(tableType, schedule, i, t, k, kbig, tables);
			}
			for(AllocThread at : allocators){at.start();}
			for(AllocThread at : allocators){
				while(at.getState()!=Thread.State.TERMINATED){
					try {
						at.join();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		
		synchronized(tables){
			for(int i=0; i<tables.length; i++){
				final AbstractKmerTableU akt=tables[i];
				if(akt==null){
					throw new RuntimeException("KmerTable allocation failed, probably due to lack of RAM: "+i+", "+tables.length);
				}
			}
		}
		
		return tables;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static class AllocThread extends Thread{
		
		AllocThread(int type_, int[] schedule_, int mod_, int div_,
				int k_, int kbig_, AbstractKmerTableU[] tables_){
			type=type_;
			schedule=schedule_;
			size=schedule[0];
			mod=mod_;
			div=div_;
			growable=schedule.length>1;
			tables=tables_;
			k=k_;
			kbig=kbig_;
		}
		
		@Override
		public void run(){
			for(int i=mod; i<tables.length; i+=div){
//				System.err.println("T"+i+" allocating "+i);
				final AbstractKmerTableU akt;
				if(type==FOREST1D){
					akt=new HashForestU(size, k, growable, false);
				}else if(type==ARRAY1D){
					akt=new HashArrayU1D(schedule, k, kbig);
//					akt=new HashArrayU1D(size, k, kbig, growable);
				}else if(type==NODE1D){
					throw new RuntimeException("Must use forest, table, or array data structure. Type="+type);
//					akt=new KmerNode2(-1, 0);
				}else if(type==FOREST2D){
					akt=new HashForestU(size, k, growable, true);
				}else if(type==ARRAY2D){
					akt=new HashArrayU2D(schedule, k, kbig);
//					akt=new HashArrayU2D(size, k, kbig, growable);
				}else if(type==NODE2D){
					throw new RuntimeException("Must use forest, table, or array data structure. Type="+type);
//					akt=new KmerNode(-1, 0);
				}else if(type==ARRAYH){
					akt=new HashArrayUHybrid(schedule, k, kbig);
//					akt=new HashArrayUHybrid(size, k, kbig, growable);
				}else{
					throw new RuntimeException("Must use forest, table, or array data structure. Type="+type);
				}
				synchronized(tables){
					tables[i]=akt;
				}
//				System.err.println("T"+i+" allocated "+i);
			}
		}
		
		private final int type;
		private final int[] schedule;
		private final int size;
		private final int mod;
		private final int div;
		private final int k;
		private final int kbig;
		private final boolean growable;
		final AbstractKmerTableU[] tables;
		
	}
	
	public WalkerU walk() {
		throw new RuntimeException("Unimplemented");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean FASTA_DUMP=true;
	public static boolean NUMERIC_DUMP=false;
	
	public static final boolean verbose=false; //slow
	public static final boolean TESTMODE=false; //slow
	
	public static final int UNKNOWN=0, ARRAY1D=1, FOREST1D=2, NODE1D=4, ARRAY2D=5, FOREST2D=6, NODE2D=8, ARRAYH=9;
	
	public static final int NOT_PRESENT=-1, HASH_COLLISION=-2;
	public static final int NO_OWNER=-1;
	
	private static final String killMessage=new String("\nThis program ran out of memory.  Try increasing the -Xmx flag and setting prealloc.");
	
}
