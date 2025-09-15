package hiseq;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.concurrent.atomic.AtomicLongArray;

import fileIO.ByteStreamWriter;
import shared.Tools;
import structures.ByteBuilder;
import structures.LongList;

public class Lane implements Iterable<Tile> {
	
	public Lane(int lane_){
		lane=lane_;
	}
	
	public MicroTile getMicroTile(int tile, int x, int y){
		return getTile(tile).get(x, y, true);
	}
	
	public MicroTile getMicroTile(int tile, int x, int y, boolean create){
		return getTile(tile).get(x, y, create);
	}
	
	public Tile getTile(int index){
		while(tiles.size()<=index){tiles.add(null);}
		Tile t=tiles.get(index);
		if(t==null){
			t=new Tile(lane, index);
			tiles.set(index, t);
		}
		return t;
	}

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
	
	public void addLists(Lane b) {
		for(int i=0; i<longLists.length; i++) {
			for(int j=0; j<longLists[i].length; j++) {
				Tools.add(longLists[i][j], b.longLists[i][j]);
			}
		}
	}

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
		
		private int pos=0;
		
	}
	
	public ArrayList<Tile> tiles=new ArrayList<Tile>();

	public boolean isEmpty() {
		return tiles.isEmpty();
	}

//	public LongList[] depthSums=new LongList[] {new LongList(151), new LongList(151)};
//	public LongList[] depthCounts=new LongList[] {new LongList(151), new LongList(151)};
//	public LongList[] matchCounts=new LongList[] {new LongList(151), new LongList(151)};
//	public LongList[] subCounts=new LongList[] {new LongList(151), new LongList(151)};

	public AtomicLongArray[] depthSums=new AtomicLongArray[] 
			{new AtomicLongArray(500), new AtomicLongArray(500)};
	public AtomicLongArray[] depthCounts=new AtomicLongArray[] 
			{new AtomicLongArray(500), new AtomicLongArray(500)};
	public AtomicLongArray[] matchCounts=new AtomicLongArray[] 
			{new AtomicLongArray(500), new AtomicLongArray(500)};
	public AtomicLongArray[] subCounts=new AtomicLongArray[] 
			{new AtomicLongArray(500), new AtomicLongArray(500)};
	
	public AtomicLongArray[][] longLists=new AtomicLongArray[][] {
		depthSums, depthCounts, matchCounts, subCounts};
	
	public int lane;
	
}
