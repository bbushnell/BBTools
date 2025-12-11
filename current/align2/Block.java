package align2;

import java.io.File;
import java.io.Serializable;

import fileIO.LoadThread;
import fileIO.ReadWrite;
import shared.KillSwitch;

/**
 * Data structure for storing alignment hit positions and start indices.
 * Used by BBTools alignment framework to efficiently manage k-mer hit lists
 * with compressed storage and serialization support for index persistence.
 *
 * @author Brian Bushnell
 * @date Dec 23, 2012
 */
public class Block implements Serializable{
	
	/** Serialization version ID */
	private static final long serialVersionUID = -1638122096023589384L;
	
	/**
	 * Creates a Block with specified capacities.
	 * Requires numStarts to be a power of 2 for efficient indexing.
	 * @param numSites_ Number of hit position sites to allocate
	 * @param numStarts_ Number of start indices (must be power of 2)
	 */
	public Block(int numSites_, int numStarts_){
		numSites=numSites_;
		numStarts=numStarts_;
		sites=new int[numSites];
		starts=new int[numStarts+1];
		assert(Integer.bitCount(numStarts)==1 && Integer.bitCount(starts.length)==2) : numStarts;
	}
	
	/**
	 * Creates a Block with existing arrays.
	 * Takes ownership of the provided arrays without copying.
	 * @param sites_ Array of hit positions
	 * @param starts_ Array of start indices (length must be power of 2)
	 */
	public Block(int[] sites_, int[] starts_){
		sites=sites_;
		starts=starts_;
		numSites=sites.length;
		numStarts=starts.length-1;
		assert(Integer.bitCount(numStarts)==1 && Integer.bitCount(starts.length)==2) : numStarts;
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
	 * @param stop Stopping index in sites array
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
	 * Batch operation for retrieving multiple hit lists efficiently.
	 *
	 * @param start Array of starting indices
	 * @param stop Array of stopping indices
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
	 * @param stop Stopping index in sites array
	 * @return Number of hits in the range
	 */
	public int length(int start, int stop){
		if(start==stop || sites[start]==-1){return 0;}
		return stop-start;
	}
	
	/**
	 * Serializes the Block to disk files.
	 * Writes sites array to fname and starts array to fname+"2.gz".
	 * Optionally compresses starts array for space efficiency.
	 *
	 * @param fname Base filename for output
	 * @param overwrite Whether to overwrite existing files
	 * @return true if write succeeded, false otherwise
	 */
	public boolean write(String fname, boolean overwrite){
		String fname2=fname+"2.gz";
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
		ReadWrite.writeObjectInThread(sites, fname, allowSubprocess);
		if(!compress){
			ReadWrite.writeObjectInThread(starts, fname+"2.gz", allowSubprocess);
		}else{
			if(copyOnWrite){
				final int[] x;
				x=new int[starts.length];
				for(int i=1; i<x.length; i++){
					x[i]=starts[i]-starts[i-1];
				}
				ReadWrite.writeObjectInThread(starts, fname+"2.gz", allowSubprocess);
			}else{
				compress(starts);
				ReadWrite.writeAsync(starts, fname+"2.gz", allowSubprocess);
				decompress(starts);
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
	 * Automatically decompresses starts array if compression was used.
	 *
	 * @param fname Base filename for input files
	 * @return Reconstructed Block object
	 */
	public static Block read(String fname){
		String fname2=fname+"2.gz";
		
		final int[] a, b;
		{
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
