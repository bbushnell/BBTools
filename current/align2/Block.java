package align2;

import java.io.File;
import java.io.Serializable;

import fileIO.LoadThread;
import fileIO.ReadWrite;
import shared.KillSwitch;

/**
 * Flat index hit storage with one start offset per key and a terminal offset.
 * Once built, key k occupies sites[starts[k]..starts[k+1]); marking its first
 * site -1 hides the list without changing offsets. Array contents remain mutable:
 * index builders fill them before readers use them. Final references do not
 * make a block immutable. Disk storage serializes the two arrays separately.
 * Access and in-place write compression require external ownership/synchronization.
 *
 * @author Brian Bushnell
 * @date Dec 23, 2012
 */
public class Block implements Serializable{

	/** Serialization version ID */
	private static final long serialVersionUID=-1638122096023589384L;

	/**
	 * Creates a Block with specified capacities.
	 * Requires numStarts to be a power of 2 of at least 2. The extra starts slot
	 * holds the terminal sites offset after index construction completes.
	 * @param numSites_ Number of hit position sites to allocate
	 * @param numStarts_ Number of start indices (must be power of 2)
	 */
	public Block(final int numSites_, final int numStarts_){
		numSites=numSites_;
		numStarts=numStarts_;
		sites=new int[numSites];
		starts=new int[numStarts+1];
		assert(Integer.bitCount(numStarts)==1 && Integer.bitCount(starts.length)==2) :
				"Index keyspace must be a power of two >=2 with a terminal start slot: keys="+numStarts+", slots="+starts.length;
	}

	/**
	 * Creates a Block with existing arrays.
	 * Retains the supplied arrays without copying; index builders continue filling them.
	 * @param sites_ Array of hit positions
	 * @param starts_ Start offsets, with length equal to a power-of-two keyspace plus 1
	 */
	public Block(final int[] sites_, final int[] starts_){
		sites=sites_;
		starts=starts_;
		numSites=sites.length;
		numStarts=starts.length-1;
		assert(Integer.bitCount(numStarts)==1 && Integer.bitCount(starts.length)==2) :
				"Index keyspace must be a power of two >=2 with a terminal start slot: keys="+numStarts+", slots="+starts.length;
	}

	/**
	 * Returns hit list for a given key.
	 * Creates a copy of the hit positions for legacy compatibility.
	 * @param key Index key for the hit list
	 * @return Array of hit positions, or null if empty
	 */
	public int[] getHitList(int key){
		int len=length(key);
		if(len==0){return null;}
		int start=starts[key];
		int[] r=KillSwitch.copyOfRange(sites, start, start+len);
		return r;
	}

	/**
	 * Returns hit list for a given range.
	 * Creates a copy of hit positions between start and stop indices.
	 *
	 * @param start Starting index in sites array
	 * @param stop Exclusive stopping index in sites array
	 * @return Array of hit positions, or null if empty
	 */
	public int[] getHitList(int start, int stop){
		int len=length(start, stop);
		if(len==0){return null;}
		assert(len>0) : len+", "+start+", "+stop;
		int[] r=KillSwitch.copyOfRange(sites, start, start+len);
		return r;
	}

	/**
	 * Returns multiple hit lists for given ranges.
	 * Allocates an outer array and a copied array for each nonempty retained list.
	 *
	 * @param start Array of starting indices
	 * @param stop Same-length array of exclusive stopping indices
	 * @return Array of hit lists corresponding to each start/stop pair
	 */
	public int[][] getHitLists(int[] start, int[] stop){
		int[][] r=new int[start.length][];
		for(int i=0; i<start.length; i++){r[i]=getHitList(start[i], stop[i]);}
		return r;
	}

	/**
	 * Returns the length of hit list for a given key.
	 * Returns 0 if the list is empty or marked as removed (first site == -1).
	 * @param key Index key for the hit list
	 * @return Number of hits in the list
	 */
	public int length(int key){
		int x=starts[key+1]-starts[key];
		if(x==0){return 0;}
		return sites[starts[key]]!=-1 ? x : 0; //Lists can be removed by making the first site -1.
	}

	/**
	 * Returns the length of hit list for a given range.
	 * Returns 0 if start equals stop or first site is marked as removed.
	 *
	 * @param start Starting index in sites array
	 * @param stop Exclusive stopping index in sites array
	 * @return Number of hits in the range
	 */
	public int length(int start, int stop){
		if(start==stop || sites[start]==-1){return 0;}
		return stop-start;
	}

	/**
	 * Serializes the Block to disk files.
	 * Writes sites array to fname and starts array to fname+"2.gz".
	 * The sites write is deferred to a worker thread. The enabled starts codec
	 * delta-encodes in place, serializes synchronously, then restores absolute
	 * offsets, including on failure. Do not access/mutate either array concurrently.
	 *
	 * @param fname Base filename for output
	 * @param overwrite Whether to overwrite existing files
	 * @return true after scheduling sites and writing starts, not a guarantee that
	 * background sites I/O finished; false on denied overwrite with assertions disabled
	 */
	public boolean write(final String fname, final boolean overwrite){
		final String fname2=fname+"2.gz";
		{
			File f=new File(fname);
			if(f.exists()){
				if(!overwrite){
					assert(false) : "Tried to overwrite file "+f.getAbsolutePath();
					return false;
				}
			}
			f=new File(fname2);
			if(f.exists()){
				if(!overwrite){
					assert(false) : "Tried to overwrite file "+f.getAbsolutePath();
					return false;
				}
			}
		}
		//TODO: Probable bug in dependency - ReadWrite.WriteObjectThread.run calls
		// addThread(-1) only after writeAsync succeeds (6e87b382, lines139-147).
		// A sites-write exception can strand the global writing counter; this
		// method does not receive that worker's failure. Failure propagation needs
		// a separate fileIO fix, not a success claim from this return value.
		ReadWrite.writeObjectInThread(sites, fname, allowSubprocess);
		if(!compress){
			ReadWrite.writeObjectInThread(starts, fname+"2.gz", allowSubprocess);
		}else{
			if(copyOnWrite){
				//align2/Block#001 FIXED: this DEAD branch (copyOnWrite is final=false) built the diff array x
				//but wrote `starts` (uncompressed) and left x[0]=0. Now matches compress()'s layout —
				//x[0]=absolute base, x[i]=diff — and writes x, so the read path would decompress it correctly.
				final int[] x;
				x=new int[starts.length];
				x[0]=starts[0];
				for(int i=1; i<x.length; i++){
					x[i]=starts[i]-starts[i-1];
				}
				ReadWrite.writeObjectInThread(x, fname+"2.gz", allowSubprocess);
			}else{
				//Race-free despite the in-place mutate-after-write: writeAsync() serializes `starts` fully
				//and synchronously (ObjectOutputStream.writeObject + close) before returning, so the
				//decompress() that follows cannot corrupt the written bytes. (Note: writeObjectInThread,
				//used for `sites` above, IS deferred — safe only because `sites` is never mutated here.)
				compress(starts);
				try{ReadWrite.writeAsync(starts, fname+"2.gz", allowSubprocess);}
				finally{decompress(starts);}
			}
		}
		return true;
	}

	/**
	 * Compresses array by converting absolute values to differences.
	 * Transforms each element to difference from previous element for space savings.
	 * @param x Array to compress in-place
	 */
	private static void compress(int[] x){
		for(int i=x.length-1; i>0; i--){
			x[i]=x[i]-x[i-1];
		}
	}

	/**
	 * Decompresses array by converting differences back to absolute values.
	 * Reverses the compression operation by computing cumulative sums.
	 * @param x Array to decompress in-place
	 */
	private static void decompress(int[] x){
		int sum=x[0];
		for(int i=1; i<x.length; i++){
			sum+=x[i];
			x[i]=sum;
		}
	}

	/**
	 * Deserializes a Block from disk files.
	 * Reads sites array from fname and starts array from fname+"2.gz".
	 * Decodes starts using the matching static codec, not a self-describing format.
	 * Wait for outstanding writes before reading; this method does not join writers.
	 *
	 * @param fname Base filename for input files
	 * @return Reconstructed Block object
	 */
	public static Block read(String fname){
		String fname2=fname+"2.gz";

		final int[] a, b;
		{
			//TODO: Probable dependency bug - LoadThread.run does not release its
			// active-reader count in finally if deserialization throws. This join
			// can finish with null output, while global reader waits may stay blocked.
			LoadThread<int[]> lta=LoadThread.load(fname, int[].class);
			b=ReadWrite.read(int[].class, fname2, false);
			lta.waitForThisToFinish();
			a=lta.output;
		}
//		{
//			LoadThread<int[]> lta=LoadThread.load(fname, int[].class);
//			LoadThread<int[]> ltb=LoadThread.load(fname2, int[].class);
//			lta.waitForThisToFinish();
//			ltb.waitForThisToFinish();
//			a=lta.output;
//			b=ltb.output;
//		}

//		int[] a=ReadWrite.read(int[].class, fname);
//		int[] b=ReadWrite.read(int[].class, fname2);

		assert(a!=null && b!=null) : a+", "+b;
		if(compress){
			int sum=b[0];
			for(int i=1; i<b.length; i++){
				sum+=b[i];
				b[i]=sum;
			}
		}
		Block r=new Block(a, b);
		assert(r.sites==a);
		assert(r.starts==b);
		return r;
	}

	public final int numSites;
	public final int numStarts;
	public final int[] sites;
	public final int[] starts;

	private static boolean allowSubprocess=false;
	private static final boolean compress=true;
	private static final boolean copyOnWrite=false;

}
