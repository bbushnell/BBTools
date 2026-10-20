package clump;

import java.io.Serializable;

import barcode.Barcode;
import hiseq.FlowcellCoordinate;
import shared.KillSwitch;
import shared.Tools;
import stream.Read;
import structures.IntList;

/**
 * Mutable pivot and optical metadata stored in a read's {@link Read#obj} field.
 * Positions identify the zero-based rightmost pivot base in the current read
 * orientation. Comparison and equality consult Clump's shared optical settings;
 * callers must keep those settings stable while sorting or comparing keys.
 * Key mutations and associated read mutations require exclusive ownership.
 *
 * @author Brian Bushnell
 * @date 2013
 */
class ReadKey implements Serializable, Comparable<ReadKey>{
	
//	public static ReadKey makeKeyIfNeeded(Read r){
//		if(r.obj==null){
//			return makeKey(r, true);
//		}
//		return (ReadKey)r.obj;
//	}
	
	/**
	 * Creates metadata for a read whose object slot is unused and whose orientation
	 * has not been swapped. Allocation failure terminates through KillSwitch.
	 *
	 * @param r The read to create a key for
	 * @param setObject Whether to store the key in r.obj
	 * @return New ReadKey for the read
	 */
	public static ReadKey makeKey(Read r, boolean setObject){
		assert(r.obj==null);
		try{
			ReadKey rk=new ReadKey(r);
			if(setObject){r.obj=rk;}
			return rk;
		}catch(OutOfMemoryError e){
			KillSwitch.memKill(e);
			throw new RuntimeException(e);
		}
	}
	
	/** Creates a ReadKey from a read with default k-mer parameters.
	 * @param r The read to create a key for */
	private ReadKey(Read r){
		this(r, 0, 0, true);
	}
	
	/**
	 * Creates a ReadKey with specific k-mer information and flowcell coordinates.
	 * Extracts lane, tile, x, y coordinates and UMI from read ID if needed.
	 *
	 * @param r The read to create a key for
	 * @param kmer_ K-mer signature value
	 * @param position_ Zero-based rightmost pivot-base position in the read
	 * @param plus_ True if k-mer is from forward strand
	 */
	private ReadKey(Read r, long kmer_, int position_, boolean plus_){
		kmer=kmer_;
		position=position_;
		clump=null;
		assert(!r.swapped());
		flipped=false;
		kmerMinusStrand=!plus_;

		//FIXED [clump/ReadKey#001]: was `clump.compareUMI` - accessing the STATIC Clump.compareUMI through the
		//instance field `clump`, which is null here (set at L67). Legal in Java (statics resolve at compile time,
		//no NPE) but a landmine: it reads as a null-deref bug, and would BECOME a real NPE if compareUMI were ever
		//made non-static. Changed to Clump.compareUMI (output-identical; matches the adjacent Clump.opticalOnly).
		if(Clump.opticalOnly || Clump.compareUMI){
			FlowcellCoordinate fc=FlowcellCoordinate.getFC();
			fc.setFrom(r.id);
			lane=fc.lane;
			tile=fc.tile;
			x=fc.x;
			y=fc.y;
			umi=fc.umi;
		}
//		expectedErrors=r.expectedErrorsIncludingMate(true);
	}
	
	/** Default constructor for subclasses */
	protected ReadKey(){}
	
	/**
	 * Updates pivot, position and strand without clearing other metadata.
	 * @param kmer_ New k-mer value
	 * @param position_ New zero-based rightmost pivot-base position
	 * @param minus_ True if k-mer is from minus strand
	 */
	public void set(long kmer_, int position_, boolean minus_){
		setKmer(kmer_);
		setPosition(position_);
//		setClump(null);
		kmerMinusStrand=minus_;
	}
	
	/**
	 * Sets the k-mer signature value.
	 * @param kmer_ New k-mer value
	 * @return The k-mer value that was set
	 */
	private long setKmer(long kmer_){
		return kmer=kmer_;
	}
	
	/**
	 * Sets the k-mer position in the read.
	 * @param position_ New position value
	 * @return The position value that was set
	 */
	private int setPosition(int position_){
		return position=position_;
	}
	
	/**
	 * Associates this ReadKey with a clump.
	 * @param clump_ The clump to associate with
	 * @return The clump that was set
	 */
	public Clump setClump(Clump clump_){
		return clump=clump_;
	}
	
	/**
	 * Sets the flipped state, ensuring it changes from current state.
	 * @param flipped_ New flipped state
	 * @return The flipped state that was set
	 */
	private boolean setFlipped(boolean flipped_){
		assert(flipped!=flipped_);
		return flipped=flipped_;
	}
	
	/** Clears pivot, position, clump and pivot-strand state. Retains flipped,
	 * optical coordinates, UMI, variants and the general-purpose flag. */
	public void clear(){
		setKmer(0);
		setPosition(0);
		setClump(null);
		kmerMinusStrand=false;
	}
	
	/**
	 * Flips the read to its reverse complement and updates position accordingly.
	 * Toggles the read's swapped bit and this key's flipped bit. For a full pivot,
	 * transforms its rightmost-base position; the pivot value and strand flag remain
	 * unchanged. Short reads retain their position.
	 * @param r The read to flip
	 * @param k K-mer length for position adjustment
	 */
	public void flip(Read r, int k){
		assert(r.swapped()==flipped);
		r.reverseComplement();
		r.setSwapped(!r.swapped());
		setFlipped(!flipped);
		if(r.length()>=k){setPosition(r.length()-position+k-2);}
		assert(r.swapped()==flipped);
	}
	
	/** Orders by descending pivot, plus before minus, then descending position.
	 * Optical mode also uses lane, optionally tile, and optionally early Y order.
	 * This ordering can distinguish keys that {@link #equals(ReadKey)} considers equal. */
	@Override
	public int compareTo(ReadKey b){
		if(kmer!=b.kmer){return kmer>b.kmer ? -1 : 1;} //Bigger kmers first...
		if(kmerMinusStrand!=b.kmerMinusStrand){return kmerMinusStrand ? 1 : -1;}
		if(position!=b.position){return position<b.position ? 1 : -1;}
		if(Clump.opticalOnly){
			if(lane!=b.lane){return lane-b.lane;}
			if(!spanTiles()){
				if(tile!=b.tile){return tile-b.tile;}
			}
			if(Clump.sortYEarly()){ //Improves speed slightly
				if(y!=b.y){return y-b.y;}
			}
		}
//		if(expectedErrors!=b.expectedErrors){
//			return expectedErrors>b.expectedErrors ? 1 : -1;//Higher quality first
//		}
		return 0;
	}
	
	/** Uses non-containment key equality. Null returns false; a foreign non-null
	 * type throws ClassCastException. Intended for homogeneous ReadKey collections. */
	@Override
	public boolean equals(Object b){
		return equals((ReadKey)b, false);
	}
	
	@Override
	//Consistent with equals (non-containment): hashes kmer+position+strand. Optical fields (lane/tile) are in
	//equals-under-opticalOnly but not here - fine, they only ADD inequality so equal->same-hash still holds.
	public int hashCode(){
		int x=(int)((kmer^position)&0xFFFFFFFFL);
		return kmerMinusStrand ? -x : x;
	}

	/** Tests geometric containment at equal pivot values, ignoring base mismatches
	 * and strand flags. Lengths and positions must describe the compared orientation. */
	public boolean physicallyContains(ReadKey b, int aLen, int bLen){
		if(bLen>aLen){return false;}
		if(kmer!=b.kmer){return false;}
		final int dif=position-b.position;
		int bStart=dif, bStop=dif+bLen;
		return bStart>=0 && bStop<=aLen;
	}
	
	/** True if this physically contains b as a prefix or suffix (ignoring mismatches).
	 * More restrictive than physicallyContains. */
	public boolean physicalAffix(ReadKey b, int aLen, int bLen){
		if(bLen>aLen){return false;}
		if(kmer!=b.kmer){return false;}
		final int dif=position-b.position;
		int bStart=dif, bStop=dif+bLen;
		return (bStart==0 || bStop==aLen) && bStart>=0 && bStop<=aLen;
	}
	
	/** Compares pivot and, outside containment mode, position and strand.
	 * Optical mode additionally requires matching lane and, unless tile spanning
	 * is enabled, tile. X/Y ordering, UMI and clump membership are not compared.
	 * @param b Other key; null is unequal
	 * @param containment Whether to ignore position and strand
	 * @return Whether the selected key fields match; not necessarily compareTo()==0
	 */
	public boolean equals(ReadKey b, boolean containment){
		if(b==null){return false;}
		if(kmer!=b.kmer){return false;}
		if(!containment && (kmerMinusStrand!=b.kmerMinusStrand || position!=b.position)){return false;}
		if(Clump.opticalOnly){
			if(lane!=b.lane){return false;}
			if(!spanTiles()){
				if(tile!=b.tile){return false;}
			}
		}
		return true;
	}
	/**
	 * Tests non-containment equality under the current optical settings.
	 * @param b ReadKey to compare with
	 * @return Whether pivot, position, strand and applicable optical fields match
	 */
	public boolean equals(ReadKey b){
		return equals(b, false);
	}
	
	@Override
	public String toString(){
		return position+","+(kmerMinusStrand ? ",t" : ",f")+","+kmer+"\t"+lane+","+tile+","+x+","+y;
	}
	
	/**
	 * Calculates physical distance between two ReadKeys on flowcell.
	 * Handles cross-tile distance calculation with configurable spanning modes.
	 * @param rkb ReadKey to calculate distance to
	 * @return Physical distance, or FlowcellCoordinate.big if incompatible
	 */
	public float distance(ReadKey rkb){
		if(lane!=rkb.lane){return FlowcellCoordinate.big;}
		
		long a=Tools.absdif(x, rkb.x), b=Tools.absdif(y, rkb.y);
		if(tile!=rkb.tile){
			if(spanTiles()){
				if(spanAdjacentOnly && Tools.absdif(tile, rkb.tile)>1){return FlowcellCoordinate.big;}
				if(spanTilesX && spanTilesY){
					return Tools.min(a, b);
				}else if(spanTilesX){
					return a;
				}else{
					return b;
				}
			}else{
				return FlowcellCoordinate.big;
			}
		}
		return (float)Math.sqrt(a*a+b*b);
	}
	
	/**
	 * Tests if another ReadKey is within specified distance.
	 * @param rkb ReadKey to test proximity to
	 * @param dist Exclusive distance threshold
	 * @return True if distance(rkb) is strictly below dist
	 */
	public boolean near(ReadKey rkb, float dist){
		return distance(rkb)<dist;
	}
	
	/**
	 * Tests proximity using minimum of X or Y coordinate differences.
	 * Requires equal lanes but ignores tile restrictions and accepts equality at
	 * the threshold. This differs from {@link #near(ReadKey, float)}.
	 *
	 * @param rkb ReadKey to test proximity to
	 * @param dist Maximum distance threshold
	 * @return True if minimum coordinate difference is within threshold
	 */
	public boolean nearXY(ReadKey rkb, float dist){
		if(lane!=rkb.lane){return false;}
		
		long a=Tools.absdif(x, rkb.x), b=Tools.absdif(y, rkb.y);
		return Tools.min(a, b)<=dist;
	}
	
	/**
	 * Tests if UMI sequences match within substitution tolerance.
	 * @param rkb ReadKey to compare UMI with
	 * @param maxSubs Maximum allowed substitutions
	 * @return True if UMIs match within substitution threshold
	 */
	public boolean umiMatches(ReadKey rkb, int maxSubs){
		if(umi==null || rkb.umi==null){return false;}
		int hdist=Barcode.hdist(umi, rkb.umi);
		return hdist<=maxSubs;
	}

	/** General purpose flag for marking ReadKeys during processing */
	public int flag;
	
	/** K-mer signature value for this read */
	public long kmer;
	/** Zero-based position of the rightmost pivot base in the current orientation. */
	public int position;
	/** True if the associated read has been reverse-complemented */
	public boolean flipped;
	/** True if the k-mer comes from the minus strand */
	public boolean kmerMinusStrand;
	/** Clump that this ReadKey belongs to */
	public Clump clump;
	/** List of variant positions for this read */
	public IntList vars;
//	public float expectedErrors;
	
	/** Flowcell coordinates, populated only when optical or UMI comparison needs them. */
	public int lane, tile, x, y;
	/** Unique molecular identifier sequence */
	public String umi;
	
	/** Historical approximate storage estimate used for clump memory planning. */
	public static final int overhead=overhead();
	/**
	 * Returns the historical fixed estimate, not a measured JVM object size.
	 * Does not include storage reachable through variants, UMI or clump references.
	 * @return Estimated bytes per key
	 */
	private static int overhead(){
		return 16+ //self
				1*(8)+ //kmer
				1*(4)+ //int fields
				2*(4)+ //booleans
				2*(8)+ //pointers
				4*(4); //flowcell coordinate
	}
	
	/** Tests if optical duplicate detection spans across tiles.
	 * @return True if spanning tiles in X or Y direction */
	public static boolean spanTiles(){return spanTilesX || spanTilesY;}
	/** Whether optical duplicate detection spans tiles in X direction */
	public static boolean spanTilesX=false;
	/** Whether optical duplicate detection spans tiles in Y direction */
	public static boolean spanTilesY=false;
	/** Whether to limit tile spanning to adjacent tiles only */
	public static boolean spanAdjacentOnly=false;
	
	/** Serialization compatibility identifier. */
	private static final long serialVersionUID=1L;
	
}
