package hiseq;

import java.util.ArrayList;
import java.util.Iterator;

import structures.ByteBuilder;

/**
 * Two-dimensional grid-based container for MicroTile objects representing
 * spatial data in a sequencing lane and tile.
 * Manages a hierarchical 2D spatial data structure that allows dynamic creation
 * and retrieval of MicroTile instances organized in a grid-like arrangement.
 *
 * @author Brian Bushnell
 */
public class Tile implements Iterable<MicroTile> {
	
	/**
	 * Constructs a Tile with specified lane and tile identifiers.
	 * @param lane_ Lane number for this tile
	 * @param tile_ Tile number within the lane
	 */
	public Tile(int lane_, int tile_){
		lane=lane_;
		tile=tile_;
	}
	
	/**
	 * Retrieves or creates a MicroTile at the specified coordinates.
	 * Calculates grid indices by dividing coordinates by grid size and manages
	 * lazy initialization of the nested ArrayList structure.
	 * @param x X-coordinate position
	 * @param y Y-coordinate position
	 * @param create Whether to create the MicroTile if it doesn't exist
	 * @return MicroTile at the specified position, or null if not found and create=false
	 */
	public MicroTile get(int x, int y, boolean create){
		//When create=true this GROWS xlist/ylist (getIndex + ylist.add/set) with NO internal lock;
		//callers mutating a shared Tile concurrently must hold the owning FlowCell's monitor
		//(see FlowCell.getMicroTile(String) contract). x,y are non-negative pixel coords (format guarantees),
		//so xindex/yindex>=0 -> the /xSize bucketing and array indices are always valid.
		final int xindex=x/xSize, yindex=y/ySize;
		ArrayList<MicroTile> ylist=getIndex(xindex, create);
		if(ylist==null || (yindex>=ylist.size() && !create)) {return null;}
		while(yindex>=ylist.size()){ylist.add(null);}//ylist is padded with NULLs (unlike xlist's non-null lists)
		MicroTile mt=ylist.get(yindex);
		if(mt==null && create){
			//Inclusive cell bounds: cell xindex spans [xindex*xSize, (xindex+1)*xSize-1] with no gap/overlap
			//vs the neighbor, so the assert below (mt.contains(x,y)) always holds for the (x,y) that created it.
			mt=new MicroTile(lane, tile, xindex*xSize, (xindex+1)*xSize-1, yindex*ySize, (yindex+1)*ySize-1);
			ylist.set(yindex, mt);
		}
		assert(mt==null || mt.contains(x,  y)) : x+", "+y+", "+xindex+", "+yindex+", "+mt;
		return mt;
	}
	
	/**
	 * Retrieves or creates the y-list for a given x-index.
	 * Manages the nested ArrayList structure by expanding as needed.
	 * @param xindex X-index in the grid
	 * @param create Whether to create missing entries
	 * @return ArrayList of MicroTiles for the specified x-index, or null if not found
	 */
	private ArrayList<MicroTile> getIndex(int xindex, boolean create){
		if(!create && xindex>=xlist.size()) {return null;}
		while(xindex>=xlist.size()){xlist.add(new ArrayList<MicroTile>());}
		ArrayList<MicroTile> ylist=xlist.get(xindex);
		return ylist;
	}
	
	@Override
	public String toString(){
		return toText(31, 0, null, null).toString();
	}
	
	/**
	 * Converts all MicroTiles in this tile to text representation.
	 * Iterates through the nested structure and calls toText on each non-null MicroTile.
	 * @param k K-mer length parameter
	 * @param HG High-depth genomic parameter
	 * @param rerf Read error reference array
	 * @param berf Base error reference array
	 * @return ByteBuilder containing text representation of all MicroTiles
	 */
	public ByteBuilder toText(int k, double HG, double[] rerf, double[] berf){
		ByteBuilder bb=new ByteBuilder();
		for(ArrayList<MicroTile> ylist : xlist){
			if(ylist!=null){
				for(MicroTile mt : ylist){
					if(mt!=null){
						mt.toText(bb, k, HG, rerf, berf);
					}
				}
			}
		}
		return bb;
	}
	
	/**
	 * Adds all MicroTiles from another Tile to this one.
	 * Uses thread-safe synchronization on both source and destination MicroTiles
	 * to prevent race conditions during data merging.
	 * @param tb Source tile whose MicroTiles will be added to this tile
	 */
	public void add(Tile tb) {
		//Merge path. The mtb/mta locks below guard the per-MicroTile field accumulation ONLY, NOT the
		//get(...,true) structural growth of THIS tile's xlist/ylist -- that relies on the caller holding
		//the FlowCell monitor (FlowCell.add wraps this in synchronized(a=Lane); the fill path never calls
		//add()). 'x' (an xlist entry) is never null: getIndex only ever inserts non-null lists, so no
		//null-guard on x is needed here (unlike toText/mtList, whose if(ylist!=null) is defensive-only).
		//'mtb' CAN be null (ylist is null-padded in get()), hence the if(mtb!=null) guard is required.
		for(ArrayList<MicroTile> x : tb.xlist) {
			for(MicroTile mtb : x) {
//				System.err.println("Adding mt "+mtb.x1+" "+mtb.y1);
				if(mtb!=null) {
					synchronized(mtb) {
						MicroTile mta=get(mtb.x1, mtb.y1, true);
						synchronized(mta) {
							mta.add(mtb);
						}
					}
				}
			}
		}
	}

	@Override
	public Iterator<MicroTile> iterator() {return mtList().iterator();}
	
	/**
	 * Creates a flattened list of all non-null MicroTiles in this tile.
	 * Traverses the nested ArrayList structure and collects all MicroTiles
	 * into a single list for easier iteration.
	 * @return ArrayList containing all non-null MicroTiles
	 */
	public ArrayList<MicroTile> mtList() {
		ArrayList<MicroTile> list=new ArrayList<MicroTile>();
		for(ArrayList<MicroTile> ylist : xlist){
			if(ylist!=null){
				for(MicroTile mt : ylist){
					if(mt!=null){
						list.add(mt);
					}
				}
			}
		}
		return list;
	}
	
	/** Nested ArrayList structure storing MicroTiles in 2D grid arrangement */
	public ArrayList<ArrayList<MicroTile>> xlist=new ArrayList<ArrayList<MicroTile>>();
	
	/** Lane number identifier for this tile */
	public int lane;
	/** Tile number identifier within the lane */
	public int tile;
	/** Fixed width size for each MicroTile in the x-dimension */
	public static int xSize=500;
	/** Fixed height size for each MicroTile in the y-dimension */
	public static int ySize=500;
}
