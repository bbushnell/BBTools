package kmer;

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
public abstract class AbstractKmerTable {
	
	/*--------------------------------------------------------------*/
	/*----------------       Abstract Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
//	/** Returns count */
//	public final int increment(long kmer){return increment(kmer, 1);}
	
	/** Returns count */
	public abstract int increment(final long kmer, final int incr);
	
//	/** Returns number of entries created */
//	public final int incrementAndReturnNumCreated(final long kmer){return incrementAndReturnNumCreated(kmer, 1);}
	
	/** Returns number of entries created.  Incr must be positive. */
	public abstract int incrementAndReturnNumCreated(final long kmer, final int incr);

	public abstract int set(long kmer, int value);

//	public abstract int set(long kmer, int[] vals);
	
	/** This is for IntList3 support with HashArrayHybridFast */
	public abstract int set(long kmer, int[] vals, int vlen);
	
	/** Returns number of kmers added */
	public abstract int setIfNotPresent(long kmer, int value);

	/**
	 * Fetch the value associated with a kmer.
	 * @param kmer
	 * @return A value.  -1 means the kmer was not present.
	 */
	public abstract int getValue(long kmer);
	
	/**
	 * Fetch the values associated with a kmer.
	 * @param kmer
	 * @param singleton A blank array of length 1.
	 * @return An array filled with values.  Values of -1 are invalid.
	 */
	public abstract int[] getValues(long kmer, int[] singleton);

	public abstract boolean contains(long kmer);
	
	public final boolean contains(long kmer, int v){
		assert(TESTMODE);
		int[] set=getValues(kmer, new int[] {-1});
		if(set==null){return false;}
		for(int s : set){
			if(s==-1){break;}
			if(s==v){return true;}
		}
		return false;
	}
	
	public final boolean contains(long kmer, int[] vals){
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
	public abstract boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, int maxcount, AtomicLong remaining);

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
	
	abstract Object get(long kmer);
	abstract void resize();
	abstract boolean canResize();
	


	/**
	 * Removes entries with a value of the limit or less.
	 * Rehashes the remainder.
	 * @return Number removed.
	 */
	abstract long regenerate(int limit);

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
	
	final static KmerNode[] allocKmerNodeArray(int len){
		KmerNode[] ret=null;
		try {
			ret=new KmerNode[len];
		} catch (OutOfMemoryError e) {
			synchronized(killMessage){
				e.printStackTrace();
				System.err.println(killMessage);
//				Shared.printMemory();
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
	public abstract int setOwner(long kmer, int newOwner);
	
	/** Reset owner to -1 if this is the current owner. */
	public abstract boolean clearOwner(long kmer, int owner);
	
	/** Return the thread ID owning this kmer, or -1. */
	public abstract int getOwner(long kmer);
	
	/** Create data structures needed for ownership representation */
	public abstract void initializeOwnership();
	
	/** Eliminate ownership data structures or set them to -1. */
	public abstract void clearOwnership();
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final StringBuilder toText(long kmer, int k){
		byte[] lookup=(Shared.AMINO_IN ? AminoAcid.numberToAcid : AminoAcid.numberToBase);
		int bits=(Shared.AMINO_IN ? 5 : 2);
		int mask=(Shared.AMINO_IN ? 31 : 3);
		StringBuilder sb=new StringBuilder(k);
		for(int i=k-1; i>=0; i--){
			int x=(int)((kmer>>(bits*i))&mask);
			sb.append((char)lookup[x]);
		}
		return sb;
	}

	static final StringBuilder toText(long kmer, int count, int k){
		StringBuilder sb=new StringBuilder(k+10);
		return toText(kmer, count, k, sb);
	}

	static final ByteBuilder toBytes(long kmer, int count, int k){
		ByteBuilder bb=new ByteBuilder(k+10);
		return toBytes(kmer, count, k, bb);
	}

	static final StringBuilder toText(long kmer, int[] values, int k){
		StringBuilder sb=new StringBuilder(k+10);
		return toText(kmer, values, k, sb);
	}

	static final ByteBuilder toBytes(long kmer, int[] values, int k){
		ByteBuilder bb=new ByteBuilder(k+10);
		return toBytes(kmer, values, k, bb);
	}
	
	static final StringBuilder toText(long kmer, int count, int k, StringBuilder sb){
		byte[] lookup=(Shared.AMINO_IN ? AminoAcid.numberToAcid : AminoAcid.numberToBase);
		int bits=(Shared.AMINO_IN ? 5 : 2);
		int mask=(Shared.AMINO_IN ? 31 : 3);
		if(FASTA_DUMP){
			sb.append('>');
			sb.append(count);
			sb.append('\n');
			for(int i=k-1; i>=0; i--){
				int x=(int)((kmer>>(bits*i))&mask);
				sb.append((char)lookup[x]);
			}
		}else{
			for(int i=k-1; i>=0; i--){
				int x=(int)((kmer>>(bits*i))&mask);
				sb.append((char)lookup[x]);
			}
			sb.append('\t');
			sb.append(count);
		}
		return sb;
	}
	
	static final StringBuilder toText(long kmer, int[] values, int k, StringBuilder sb){
		byte[] lookup=(Shared.AMINO_IN ? AminoAcid.numberToAcid : AminoAcid.numberToBase);
		int bits=(Shared.AMINO_IN ? 5 : 2);
		int mask=(Shared.AMINO_IN ? 31 : 3);
		if(FASTA_DUMP){
			sb.append('>');
			for(int i=0; i<values.length; i++){
				int x=values[i];
				if(x==-1){break;}
				if(i>0){sb.append(',');}
				sb.append(x);
			}
			sb.append('\n');
			for(int i=k-1; i>=0; i--){
				int x=(int)((kmer>>(bits*i))&mask);
				sb.append((char)lookup[x]);
			}
		}else{
			for(int i=k-1; i>=0; i--){
				int x=(int)((kmer>>(bits*i))&mask);
				sb.append((char)lookup[x]);
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
	
	public static final ByteBuilder toBytes(long kmer, long count, int k, ByteBuilder bb){
		byte[] lookup=(Shared.AMINO_IN ? AminoAcid.numberToAcid : AminoAcid.numberToBase);
		int bits=(Shared.AMINO_IN ? 5 : 2);
		int mask=(Shared.AMINO_IN ? 31 : 3);
		if(FASTA_DUMP){
			bb.append('>');
			bb.append(count);
			bb.nl();
			for(int i=k-1; i>=0; i--){
				int x=(int)((kmer>>(bits*i))&mask);
				bb.append(lookup[x]);
			}
//			assert(false) : kmer+"->\n"+bb+"\n"+AminoAcid.kmerToStringAA(kmer, k);
		}else if(NUMERIC_DUMP){
			bb.append(Long.toHexString(kmer));
			bb.tab();
			bb.append(count);
		}else{
			for(int i=k-1; i>=0; i--){
				int x=(int)((kmer>>(bits*i))&mask);
				bb.append(lookup[x]);
			}
			bb.tab();
			bb.append(count);
		}
		return bb;
	}
	
	public static final ByteBuilder toBytes(long kmer, int[] values, int k, ByteBuilder bb){
		byte[] lookup=(Shared.AMINO_IN ? AminoAcid.numberToAcid : AminoAcid.numberToBase);
		int bits=(Shared.AMINO_IN ? 5 : 2);
		int mask=(Shared.AMINO_IN ? 31 : 3);
		if(FASTA_DUMP){
			bb.append('>');
			for(int i=0; i<values.length; i++){
				int x=values[i];
				if(x==-1){break;}
				if(i>0){bb.append(',');}
				bb.append(x);
			}
			bb.nl();
			for(int i=k-1; i>=0; i--){
				int x=(int)((kmer>>(bits*i))&mask);
				bb.append(lookup[x]);
			}
		}else if(NUMERIC_DUMP){
			bb.append(Long.toHexString(kmer));
			bb.tab();
			for(int i=0; i<values.length; i++){
				int x=values[i];
				if(x==-1){break;}
				if(i>0){bb.append(',');}
				bb.append(x);
			}
		}else{
			for(int i=k-1; i>=0; i--){
				int x=(int)((kmer>>(bits*i))&mask);
				bb.append(lookup[x]);
			}
			bb.tab();
			for(int i=0; i<values.length; i++){
				int x=values[i];
				if(x==-1){break;}
				if(i>0){bb.append(',');}
				bb.append(x);
			}
		}
		return bb;
	}
	
//	static void appendKmerText(long kmer, int count, int k, StringBuilder sb){
//		sb.setLength(0);
//		toText(kmer, count, k, sb);
//		sb.append('\n');
//	}
	
	static void appendKmerText(long kmer, int count, int k, ByteBuilder bb){
		bb.setLength(0);
		toBytes(kmer, count, k, bb);
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
	 * @param mask
	 * @return The preallocated table
	 */
	public static final AbstractKmerTable[] preallocate(int ways, int tableType, int[] schedule, long mask){

		final AbstractKmerTable[] tables=new AbstractKmerTable[ways];
		
		{
			shared.Timer tm=new shared.Timer();
			final int t=Tools.max(1, Tools.min(Shared.threads(), 2, ways)); //More than 2 still improves allocation time, but only slightly; ~25% faster at t=4.
			final AllocThread[] allocators=new AllocThread[t];
			for(int i=0; i<t; i++){
				allocators[i]=new AllocThread(tableType, schedule, i, t, mask, tables);
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
			tm.stop();
			if(AbstractKmerTableSet.DISPLAY_PROGRESS){System.err.println(tm);}
		}
		
		synchronized(tables){
			for(int i=0; i<tables.length; i++){
				final AbstractKmerTable akt=tables[i];
				if(akt==null){
					throw new RuntimeException("KmerTable allocation failed, probably due to lack of RAM: "+i+", "+tables.length);
				}
			}
		}
		
		return tables;
	}
	
	public static final long makeMiddleMask(int k, boolean amino){
		final boolean odd=((k&1)==1);
		final int bitsPerSymbol=(amino ? 5 : 2);
		final int shift=bitsPerSymbol*((k-1)/2);
		final long middle=~((-1L)<<(odd ? bitsPerSymbol : 2*bitsPerSymbol));
		return ~(middle<<shift);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static class AllocThread extends Thread{
		
		AllocThread(int type_, int[] schedule_, int mod_, int div_,
				long mask_, AbstractKmerTable[] tables_){
			type=type_;
			schedule=schedule_;
			size=schedule[0];
			mod=mod_;
			div=div_;
			mask=mask_;
			growable=schedule.length>1;
			tables=tables_;
		}
		
		@Override
		public void run(){
			//Initialize tables
			long sum=0;

//			Shared.printMemory();}
			for(int i=mod; i<tables.length; i+=div){
//				System.err.println("T"+i+" allocating "+i);
				final AbstractKmerTable akt;
				if(type==FOREST1D){
					akt=new HashForest(size, growable, false);
				}else if(type==TABLE){
					akt=new KmerTable(size, growable);
				}else if(type==ARRAY1D){
					akt=new HashArray1D(schedule, mask);
//					long len=((HashArray1D)akt).arrayLength();
//					long mem=((HashArray1D)akt).calcMem();
//					sum+=mem;
//					akt=new HashArray1D(size, -1, mask, growable);//TODO: Set maxSize
				}else if(type==NODE1D){
					throw new RuntimeException("Must use forest, table, or array data structure. Type="+type);
//					akt=new KmerNode2(-1, 0);
				}else if(type==FOREST2D){
					akt=new HashForest(size, growable, true);
				}else if(type==TABLE2D){
					throw new RuntimeException("Must use forest, table, or array data structure. Type="+type);
				}else if(type==ARRAY2D){
					akt=new HashArray2D(schedule, mask);
//					akt=new HashArray2D(size, -1, mask, growable);//TODO: Set maxSize
				}else if(type==NODE2D){
					throw new RuntimeException("Must use forest, table, or array data structure. Type="+type);
//					akt=new KmerNode(-1, 0);
				}else if(type==ARRAYH){
					akt=new HashArrayHybrid(schedule, mask);
//					akt=new HashArrayHybrid(size, -1, mask, growable);//TODO: Set maxSize
				}else if(type==ARRAYHF){
					akt=new HashArrayHybridFast(schedule, mask);
//					akt=new HashArrayHybrid(size, -1, mask, growable);//TODO: Set maxSize
				}else{
					throw new RuntimeException("Must use forest, table, or array data structure. Type="+type);
				}
				synchronized(tables){
					tables[i]=akt;
				}
//				System.err.println("T"+i+" allocated "+i);
				
//				if(i%100==0){
//					System.err.println("Allocated "+(sum/1000000)+"MB");
//					Shared.printMemory();
//				}
			}
//			Shared.printMemory();
//			assert(false) : ("Allocated "+(sum/1000000)+"MB");
		}
		
		private final int type;
		private final int[] schedule;
		private final int size;
		private final int mod;
		private final int div;
		private final long mask;
		private final boolean growable;
		final AbstractKmerTable[] tables;
		
	}
	
	public Walker walk() {
		throw new RuntimeException("Unimplemented");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean FASTA_DUMP=true;
	public static boolean NUMERIC_DUMP=false;
	public static boolean TWO_PASS_RESIZE=false;
	
	public static final boolean verbose=false;
	public static final boolean TESTMODE=false; //123 SLOW!
	
	public static final int UNKNOWN=0, ARRAY1D=1, FOREST1D=2, TABLE=3, NODE1D=4, ARRAY2D=5, FOREST2D=6, TABLE2D=7, NODE2D=8, ARRAYH=9, ARRAYHF=10;
	
	public static final int NOT_PRESENT=-1, HASH_COLLISION=-2;
	public static final int NO_OWNER=-1;
	
	private static final String killMessage=new String("\nThis program ran out of memory.  Try increasing the -Xmx flag and setting prealloc.");
	
}
