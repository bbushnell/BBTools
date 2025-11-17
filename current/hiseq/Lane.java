package hiseq;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.concurrent.atomic.AtomicLongArray;

import fileIO.ByteStreamWriter;
import shared.Tools;

/**
 * Represents a sequencing lane in high-throughput sequencing datasets, managing
 * tiles, microtiles, and associated genomic metrics. Serves as a container for
 * organizing sequencing data at multiple granularities (lane, tile, microtile)
 * and provides methods for data aggregation, iteration, and high-depth genomic
 * kmer analysis.
 *
 * Key features include dynamic tile storage, thread-safe tile addition, genomic
 * metrics calculation, and parallel data structures using AtomicLongArray for
 * thread-safe depth, match, and count tracking.
 *
 * @author Brian Bushnell
 */
public class Lane implements Iterable<Tile> {
	
	/** Constructs a Lane with the specified lane number.
	 * @param lane_ The lane identifier number */
	public Lane(int lane_){
		lane=lane_;
	}
	
	/**
	 * Retrieves a MicroTile at the specified coordinates, creating it if needed.
	 * @param tile The tile index
	 * @param x The x coordinate within the tile
	 * @param y The y coordinate within the tile
	 * @return The MicroTile at the specified location
	 */
	public MicroTile getMicroTile(int tile, int x, int y){
		return getTile(tile).get(x, y, true);
	}
	
	/**
	 * Retrieves a MicroTile at the specified coordinates with optional creation.
	 * @param tile The tile index
	 * @param x The x coordinate within the tile
	 * @param y The y coordinate within the tile
	 * @param create Whether to create the MicroTile if it doesn't exist
	 * @return The MicroTile at the specified location, or null if create is false
	 * and the MicroTile doesn't exist
	 */
	public MicroTile getMicroTile(int tile, int x, int y, boolean create){
		return getTile(tile).get(x, y, create);
	}
	
	/**
	 * Retrieves the tile at the specified index, creating it if necessary.
	 * Expands the tiles array as needed to accommodate the index.
	 * @param index The tile index
	 * @return The Tile at the specified index
	 */
	public Tile getTile(int index){
		while(tiles.size()<=index){tiles.add(null);}
		Tile t=tiles.get(index);
		if(t==null){
			t=new Tile(lane, index);
			tiles.set(index, t);
		}
		return t;
	}

	/**
	 * Adds data from another lane to this lane in a thread-safe manner.
	 * Synchronizes on both source and destination tiles during addition.
	 * @param b The source lane to add data from
	 */
	public void add(Lane b) {
		for(Tile tb : b.tiles) {
			if(tb!=null) {
				synchronized(tb) {
					Tile ta=getTile(tb.tile);
					synchronized(ta) {
						ta.add(tb);
					}
				}
			}
		}
		addLists(b);
	}
	
	/**
	 * Adds atomic long array data from another lane to this lane's arrays.
	 * Merges depth sums, counts, match counts, and substitution counts.
	 * @param b The source lane containing data to add
	 */
	public void addLists(Lane b) {
		for(int i=0; i<longLists.length; i++) {
			for(int j=0; j<longLists[i].length; j++) {
				Tools.add(longLists[i][j], b.longLists[i][j]);
			}
		}
	}

	/**
	 * Prints lane data to a ByteStreamWriter with calculated high-depth genomic
	 * rate and error rate arrays.
	 * @param bsw The output stream writer
	 * @param k The kmer length for high-depth genomic calculation
	 * @param rerf Read error rate factors array
	 * @param berf Base error rate factors array
	 */
	public void print(ByteStreamWriter bsw, int k, double[] rerf, double[] berf) {
		double HG=calcHighDepthGenomic(k);
		for(Tile tile : tiles){
			if(tile!=null){
				bsw.print(tile.toText(k, HG, rerf, berf));
			}
		}
	}
	
	/**
	 * Calculates the rate of high-depth genomic kmers in the lane,
	 * using the observed rates of unique kmers and alignment errors.
	 * @param k Kmer length.
	 * @return HG Rate of high depth genomic kmers.
	 */
	public double calcHighDepthGenomic(int k) {
		long hits=0;
		long misses=0;
		long errors=0;
		long alignedBases=0;
		long mts=0;
		for(Tile t : this) {
			for(MicroTile mt : t) {
//				if(mt!=null) {
					hits+=mt.hits;
					misses+=mt.misses;
					errors+=mt.baseErrorCount;
					alignedBases+=mt.alignedBaseCount;
					mts++;
//				}
			}
		}
		if(mts<1 || alignedBases<1 || hits+misses<1) {return 0;}
		assert(alignedBases>0) : "No alignment data.";
		assert(hits+misses>0) : "No kmer data.";
		//PhiX error rate
		double E=errors/(double)alignedBases;
		//Per-base correctness chance
		double C=1-E;
		//Prob of kmer correctness
		double P=Math.pow(C, k);
		//Unique fraction
		double U=misses/(double)(hits+misses);
		//Non-unique fraction
		double NU=1-U;
		//High depth genomic kmer fraction
		double HG=Tools.mid(0.0001, 0.9999, NU/P); //Otherwise it can exceed 1 for high-depth things like PhiX.
		assert(HG>0 && HG<=1) : "\nE="+E+", P="+P+", U="+U+", NU="+NU+", HG="+HG
			+"\nmts="+mts+", hits="+hits+", misses="+misses+", errors="+errors+", bases="+alignedBases+", k="+k;
		return HG;
	}

	@Override
	public Iterator<Tile> iterator() {
		return new TileIterator();
	}
	
	private final class TileIterator implements Iterator<Tile> {

		/** Default constructor for TileIterator inner class. */
		TileIterator(){}
		
		@Override
		public boolean hasNext() {
			while(pos<tiles.size() && tiles.get(pos)==null) {pos++;}
			return pos<tiles.size();
		}

		@Override
		public Tile next() {
			Tile t=(hasNext() ? tiles.get(pos) : null);
			pos++;
			return t;
		}
		
		/** Current position in tile iteration */
		private int pos=0;
		
	}
	
	/** Dynamic storage for tiles in this lane */
	public ArrayList<Tile> tiles=new ArrayList<Tile>();

	/** Checks if this lane contains no tiles.
	 * @return true if the tiles list is empty, false otherwise */
	public boolean isEmpty() {
		return tiles.isEmpty();
	}

//	public LongList[] depthSums=new LongList[] {new LongList(151), new LongList(151)};
//	public LongList[] depthCounts=new LongList[] {new LongList(151), new LongList(151)};
//	public LongList[] matchCounts=new LongList[] {new LongList(151), new LongList(151)};
//	public LongList[] subCounts=new LongList[] {new LongList(151), new LongList(151)};

	/** Thread-safe arrays tracking depth sums across two dimensions */
	public AtomicLongArray[] depthSums=new AtomicLongArray[] 
			{new AtomicLongArray(500), new AtomicLongArray(500)};
	/** Thread-safe arrays tracking depth counts across two dimensions */
	public AtomicLongArray[] depthCounts=new AtomicLongArray[] 
			{new AtomicLongArray(500), new AtomicLongArray(500)};
	/** Thread-safe arrays tracking match counts across two dimensions */
	public AtomicLongArray[] matchCounts=new AtomicLongArray[] 
			{new AtomicLongArray(500), new AtomicLongArray(500)};
	/** Thread-safe arrays tracking substitution counts across two dimensions */
	public AtomicLongArray[] subCounts=new AtomicLongArray[] 
			{new AtomicLongArray(500), new AtomicLongArray(500)};
	
	/** Consolidated array of all atomic long arrays for bulk operations */
	public AtomicLongArray[][] longLists=new AtomicLongArray[][] {
		depthSums, depthCounts, matchCounts, subCounts};
	
	/** The lane identifier number */
	public int lane;
	
}
