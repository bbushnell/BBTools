package hiseq;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import shared.Tools;
import structures.FloatList;
import structures.Point;

/**
 * Represents an Illumina sequencing flow cell with comprehensive lane and
 * tile-level statistics management. Provides data aggregation, statistical
 * computation, and adaptive processing for genomic read analysis across
 * multiple lanes and tiles.
 *
 * @author Brian Bushnell
 */
public class FlowCell {
	
	/** Creates a FlowCell by loading data from a tile dump file.
	 * @param fname Path to the tile dump file to load */
	public FlowCell(String fname){
		TileDump.loadDump(fname, this);
	}
	
	/** Creates a FlowCell with the specified k-mer length.
	 * @param k_ K-mer length for sequence analysis */
	public FlowCell(int k_){k=k_;}
	
	/**
	 * Retrieves a MicroTile by parsing an Illumina read identifier.
	 * This method is NOT thread-safe due to shared IlluminaHeaderParser2.
	 * @param id Illumina read identifier string
	 * @return MicroTile corresponding to the parsed coordinates
	 */
	public MicroTile getMicroTile(String id) {//This method is NOT threadsafe
		ihp.parse(id);
		return getMicroTile(ihp);
	}
	
	/**
	 * Retrieves a MicroTile using a pre-parsed header parser.
	 * @param ihp Illumina header parser containing coordinate information
	 * @return MicroTile at the specified coordinates
	 */
	public MicroTile getMicroTile(IlluminaHeaderParser2 ihp){
		return getLane(ihp.lane()).getMicroTile(ihp.tile(), ihp.xPos(), ihp.yPos(), true);
	}
	
	/**
	 * Retrieves a MicroTile at the specified coordinates, creating if needed.
	 *
	 * @param lane Lane number
	 * @param tile Tile number
	 * @param x X coordinate
	 * @param y Y coordinate
	 * @return MicroTile at the specified location
	 */
	public MicroTile getMicroTile(int lane, int tile, int x, int y){
		return getLane(lane).getMicroTile(tile, x, y, true);
	}
	
	/**
	 * Retrieves a MicroTile at the specified coordinates.
	 *
	 * @param lane Lane number
	 * @param tile Tile number
	 * @param x X coordinate
	 * @param y Y coordinate
	 * @param create If true, creates MicroTile if it doesn't exist
	 * @return MicroTile at the specified location, or null if not found and create is false
	 */
	public MicroTile getMicroTile(int lane, int tile, int x, int y, boolean create){
		return getLane(lane).getMicroTile(tile, x, y, create);
	}
	
	/**
	 * Retrieves a Lane by number, creating new lanes as needed.
	 * Dynamically expands the lanes list to accommodate the requested lane.
	 * @param lane Lane number to retrieve
	 * @return Lane object at the specified index
	 */
	public Lane getLane(int lane){
		while(lanes.size()<=lane){lanes.add(new Lane(lanes.size()));}
		return lanes.get(lane);
	}

	/**
	 * Converts all MicroTiles in the FlowCell to a flat list.
	 * Traverses all lanes, tiles, and coordinate positions to collect
	 * non-null MicroTiles.
	 * @return ArrayList containing all MicroTiles in the FlowCell
	 */
	public ArrayList<MicroTile> toList() {
		ArrayList<MicroTile> mtList=new ArrayList<MicroTile>();
		for(Lane lane : lanes){
			if(lane!=null){
				for(Tile tile : lane.tiles){
					if(tile!=null){
						for(ArrayList<MicroTile> ylist : tile.xlist){
							if(ylist!=null){
								for(MicroTile mt : ylist){
									if(mt!=null){
										mtList.add(mt);
									}
								}
							}
						}
					}
				}
			}
		}
		return mtList;
	}
	
	/**
	 * Calculates comprehensive statistics for the entire FlowCell.
	 * Processes all MicroTiles to compute read counts, alignment rates,
	 * error rates, quality metrics, and linear regression formulas for
	 * unique percentage vs error rate relationships.
	 *
	 * @return ArrayList of all processed MicroTiles
	 */
	public ArrayList<MicroTile> calcStats(){
		ArrayList<MicroTile> mtList=toList();
		readsProcessed=basesProcessed=0;
		readsAligned=basesAligned=0;
		readErrors=baseErrors=0;
		for(MicroTile mt : mtList){
			mt.process();
			readsProcessed+=mt.readCount;
			basesProcessed+=mt.baseCount;
			readsAligned+=mt.alignedReadCount;
			basesAligned+=mt.alignedBaseCount;
			readErrors+=mt.readErrorCount;
			baseErrors+=mt.baseErrorCount;
		}
		final double mtDiv=Tools.max(1, mtList.size());
		avgReads=readsProcessed/mtDiv;
		avgAlignedReads=readsAligned/mtDiv;
		minCountToUse=(long)Tools.min(500, avgReads*0.25f);
		int toKeep=0;
		for(MicroTile mt : mtList){
			if(mt.readCount>=minCountToUse){toKeep++;}
		}
		
		FloatList avgQualityList=new FloatList(toKeep);
		FloatList avgUniqueList=new FloatList(toKeep);
		FloatList avgDepthList=new FloatList(toKeep);
		FloatList avgErrorFreeList=new FloatList(toKeep);
		FloatList avgPolyGList=new FloatList(toKeep);
		FloatList avgGList=new FloatList(toKeep);
		

		ArrayList<Point> readPoints=new ArrayList<Point>();
		ArrayList<Point> basePoints=new ArrayList<Point>();
		
		for(MicroTile mt : mtList){
			if(mt!=null && mt.readCount>=minCountToUse){
				double up=mt.uniquePercent();
				double rer=mt.readErrorRate();
				double ber=mt.baseErrorRate();
				avgQualityList.add((float)mt.averageReadQualityByProb());
				avgUniqueList.add((float)up);
				avgDepthList.add((float)mt.depth());
				avgErrorFreeList.add((float)mt.percentErrorFree());
				avgPolyGList.add((float)mt.polyGPercent());
				avgGList.add((float)mt.avgG());
				readPoints.add(new Point(up, rer));
				basePoints.add(new Point(up, Math.sqrt(ber)));
			}
		}
		
		if(readsAligned>1000) {
			Collections.sort(readPoints);
			Collections.sort(basePoints);
			uniqueToReadErrorRateFormula=Tools.linearRegression(readPoints, 0.001, 0.999);
			uniqueToBaseErrorRateFormula=Tools.linearRegression(basePoints, 0.001, 0.999);
//			System.err.println("Calculated "+Arrays.toString(uniqueToBaseErrorRateFormula)+" from "
//					+basePoints.size()+" points, "+basePoints.get(0)+"~"+basePoints.get(basePoints.size()-1));
		}else {
//			assert(false) : readsProcessed+", "+basesProcessed+", "+readsAligned+", "+readPoints.get(0);
			readPoints=basePoints=null;
			uniqueToReadErrorRateFormula=uniqueToBaseErrorRateFormula=null;
		}
		
		avgQuality=avgQualityList.mean();
		avgUnique=avgUniqueList.mean();
		avgDepth=avgDepthList.mean();
		avgErrorFree=avgErrorFreeList.mean();
		avgPolyG=avgPolyGList.mean();
		avgG=avgGList.mean();
		
		stdQuality=avgQualityList.stdev();
		stdUnique=avgUniqueList.stdev();
		stdDepth=avgDepthList.stdev();
		stdErrorFree=avgErrorFreeList.stdev();
		stdPolyG=avgPolyGList.stdev();
		stdG=avgGList.stdev();
		
		return mtList;
	}
	
	/**
	 * Adaptively widens tile dimensions to reach target read count.
	 * Iteratively increases tile sizes until average reads per tile
	 * meets or exceeds the target, or no improvement is achieved.
	 *
	 * @param targetReads Target average reads per tile
	 * @return FlowCell with widened tiles, or original if target already met
	 */
	public FlowCell widenToTargetReads(int targetReads){
		if(readsProcessed<1){
			System.err.println("Warning: Zero reads processed.");
			return this;
		}
		if(readsProcessed<targetReads){
			return this;
		}
		FlowCell fc=this;
		while(fc.avgReads<targetReads){
			int x2=Tile.xSize, y2=Tile.ySize;
			FlowCell fc2=fc.widen(true);
			fc2.calcStats();
			if(fc2.avgReads<=fc.avgReads){
//				unwiden();
				fc2=null;
				Tile.xSize=x2;
				Tile.ySize=y2;
				return fc;
			}
			fc=fc2;
		}
		return fc.setFrom(this);
	}
	
	/**
	 * Adaptively widens tile dimensions to reach target aligned read count.
	 * Similar to widenToTargetReads but focuses on aligned reads specifically.
	 * @param targetAlignedReads Target average aligned reads per tile
	 * @return FlowCell with widened tiles, or original if target already met
	 */
	public FlowCell widenToTargetAlignedReads(int targetAlignedReads){
		if(readsAligned<1){
			System.err.println("Warning: Zero aligned reads processed.");
			return this;
		}
		if(readsAligned<targetAlignedReads){
			System.err.println("Warning: Below target aligned reads ("+readsAligned+"<"+targetAlignedReads+")");
			return this;
		}
		FlowCell fc=this;
		while(fc.avgAlignedReads<targetAlignedReads){
			int x2=Tile.xSize, y2=Tile.ySize;
			FlowCell fc2=fc.widen(false);
			fc2.calcStats();
			if(fc2.avgAlignedReads<=fc.avgAlignedReads){
//				unwiden();
				fc2=null;
				Tile.xSize=x2;
				Tile.ySize=y2;
				return fc;
			}
			fc=fc2;
		}
		return fc.setFrom(this);
	}
	
//	/** Reverses tile widening by halving the larger dimension.
//	 * Reduces either xSize or ySize by half, choosing the larger dimension.
//	 * NOTE:  This is broken, causes stale microtiles. */
//	public void unwiden(){
//		if(Tile.xSize>Tile.ySize){Tile.ySize/=2;}
//		else{Tile.xSize/=2;}
//	}
	
	/**
	 * Widens tile dimensions by doubling the larger dimension.
	 * If x >= y, doubles y; otherwise doubles x.
	 * @param loud If true, prints widening information to stderr
	 * @return New FlowCell with wider tile dimensions
	 */
	public FlowCell widen(boolean loud){
		final int x=Tile.xSize, y=Tile.ySize;
		final int x2=x>=y ? x : 2*x;
		final int y2=x>=y ? 2*y : y;
		return widen(x2, y2, loud);
	}
	
	/**
	 * Widens tile dimensions to specified x and y sizes.
	 * Creates a new FlowCell with larger tiles and redistributes all
	 * existing MicroTile data to the new coordinate system.
	 *
	 * @param x2 New x dimension for tiles
	 * @param y2 New y dimension for tiles
	 * @param loud If true, prints widening information to stderr
	 * @return New FlowCell with specified tile dimensions
	 */
	public FlowCell widen(int x2, int y2, boolean loud){
//		assert(x2>Tile.xSize || y2>Tile.ySize);
		if(x2<=Tile.xSize && y2<=Tile.ySize) {return this;}
		Tile.xSize=Tools.max(x2, Tile.xSize);
		Tile.ySize=Tools.max(y2, Tile.ySize);
		if(loud) {System.err.println("Widening to "+Tile.xSize+"x"+Tile.ySize);}
		ArrayList<MicroTile> list=toList();
		FlowCell fc=new FlowCell(k).setFrom(this);
		for(MicroTile mt : list){
			MicroTile mt2=fc.getMicroTile(mt.lane, mt.tile, mt.x1, mt.y1);
			mt2.add(mt);
		}
		for(Lane lane : lanes) {
			if(lane!=null) {
				Lane lane2=fc.getLane(lane.lane);
				lane2.addLists(lane);
			}
		}
		return fc;
	}

	/**
	 * Applies spatial blur to MicroTile data by incorporating adjacent tiles.
	 * Each MicroTile is weighted 8x and combined with up, down, left, and right
	 * neighbors, then scaled by 0.125 to achieve ~50% above original counts.
	 */
	public void blur() {
		FlowCell fc=new FlowCell(k);
		fc.add(this);
		ArrayList<MicroTile> list=this.toList();
		for(MicroTile mt : list) {
			if(mt.readCount<1) {continue;}
			int tile=mt.tile, lane=mt.lane, x1=mt.x1, x2=mt.x2, y1=mt.y1, y2=mt.y2;
			int x0=x1-1, y0=y1-1, x3=x2+1, y3=y2+1;
//			MicroTile center=fc.getMicroTile(lane, tile, x1, y1, true);
			MicroTile up=(y0<0 ? null : fc.getMicroTile(lane, tile, x1, y0, false));
			MicroTile down=fc.getMicroTile(lane, tile, x1, y3, false);
			MicroTile left=(x0<0 ? null : fc.getMicroTile(lane, tile, x0, y1, false));
			MicroTile right=fc.getMicroTile(lane, tile, x3, y1, false);
//			mt.add(center);
//			mt.add(center);
//			mt.add(center);
//			mt.add(center);
//			mt.add(center);
//			mt.add(center);
//			mt.add(center);
//			if(up!=null) {mt.add(up);}
//			if(down!=null) {mt.add(down);}
//			if(left!=null) {mt.add(left);}
//			if(right!=null) {mt.add(right);}
			
			mt.multiplyBy(8);
			if(up!=null) {mt.add(up);}
			if(down!=null) {mt.add(down);}
			if(left!=null) {mt.add(left);}
			if(right!=null) {mt.add(right);}
			mt.multiplyBy(0.125); //Counts should end up around 50% above original
		}
	}
	
	/**
	 * Saves the FlowCell data to a tile dump file.
	 * @param fname Output filename for the dump
	 * @param overwrite If true, overwrites existing files
	 */
	public void dump(String fname, boolean overwrite) {
		TileDump.write(this, fname, overwrite);
	}
	
	/**
	 * Adds data from another FlowCell to this one.
	 * Thread-safe operation that synchronizes on both source and destination
	 * lanes during the merge operation.
	 * @param fcb FlowCell to merge into this one
	 */
	public void add(FlowCell fcb) {
		for(Lane b : fcb.lanes) {
			if(b!=null) {
				synchronized(b) {
					Lane a=getLane(b.lane);
					synchronized(a) {
						a.add(b);
					}
				}
			}
		}
	}
	
	/**
	 * Copies statistical data and metadata from another FlowCell.
	 * Transfers read counts, error counts, coordinate boundaries, and
	 * regression formulas from the source FlowCell.
	 *
	 * @param fc Source FlowCell to copy from
	 * @return This FlowCell with updated data
	 */
	FlowCell setFrom(FlowCell fc) {
		readsProcessed=fc.readsProcessed;
		basesProcessed=fc.basesProcessed;
		readsAligned=fc.readsAligned;
		basesAligned=fc.basesAligned;
		readErrors=fc.readErrors;
		baseErrors=fc.baseErrors;
		
		name=fc.name;
		xMin=fc.xMin;
		xMax=fc.xMax;
		yMin=fc.yMin;
		yMax=fc.yMax;
		tMin=fc.tMin;
		tMax=fc.tMax;
		k=fc.k;

		uniqueToReadErrorRateFormula=fc.uniqueToReadErrorRateFormula==null ? null : 
			Arrays.copyOf(fc.uniqueToReadErrorRateFormula, 2);
		uniqueToBaseErrorRateFormula=fc.uniqueToBaseErrorRateFormula==null ? null : 
			Arrays.copyOf(fc.uniqueToBaseErrorRateFormula, 2);
		return this;
	}

	/** Calculates the overall alignment rate for the FlowCell.
	 * @return Fraction of reads that aligned successfully */
	public double alignmentRate() {
		return readsAligned/(double)readsProcessed;
	}
	
	/** Calculates the overall base error rate for aligned reads.
	 * @return Base error rate as fraction of aligned bases */
	public double baseErrorRate() {
		return baseErrors/(double)(Tools.max(basesAligned, 1));
	}
	
	/** List of sequencing lanes in this FlowCell */
	public ArrayList<Lane> lanes=new ArrayList<Lane>();
	
	/** Total number of reads processed across all lanes */
	long readsProcessed;
	/** Total number of bases processed across all lanes */
	long basesProcessed;
	/** Total number of reads that aligned successfully */
	long readsAligned;
	/** Total number of bases in aligned reads */
	long basesAligned;
	/** Total number of read-level errors detected */
	long readErrors;
	/** Total number of base-level errors detected */
	long baseErrors;
	
	/** Name identifier for this FlowCell */
	String name;

	/** Minimum x coordinate boundary */
	int xMin=-1;
	/** Maximum x coordinate boundary */
	int xMax=-1;
	/** Minimum y coordinate boundary */
	int yMin=-1;
	/** Maximum y coordinate boundary */
	int yMax=-1;
	/** Minimum tile number boundary */
	int tMin=-1;
	/** Maximum tile number boundary */
	int tMax=-1;
	/** K-mer length used for sequence analysis */
	int k=31;

	/** Average number of reads per MicroTile */
	public double avgReads;
	/** Average number of aligned reads per MicroTile */
	public double avgAlignedReads;
	/** Minimum read count threshold for including MicroTiles in statistics */
	public double minCountToUse;
	
	/** Average quality score across all MicroTiles */
	public double avgQuality;
	/** Average unique read percentage across all MicroTiles */
	public double avgUnique;
	/** Average sequencing depth across all MicroTiles */
	public double avgDepth;
	/** Average percentage of error-free reads across all MicroTiles */
	public double avgErrorFree;
	/** Average poly-G percentage across all MicroTiles */
	public double avgPolyG;
	/** Average G-content across all MicroTiles */
	public double avgG;
	
	/** Standard deviation of quality scores across MicroTiles */
	public double stdQuality;
	/** Standard deviation of unique read percentages across MicroTiles */
	public double stdUnique;
	/** Standard deviation of sequencing depth across MicroTiles */
	public double stdDepth;
	/** Standard deviation of error-free percentages across MicroTiles */
	public double stdErrorFree;
	/** Standard deviation of poly-G percentages across MicroTiles */
	public double stdPolyG;
	/** Standard deviation of G-content across MicroTiles */
	public double stdG;
	
	/** Linear regression coefficients for unique percentage vs read error rate */
	public double[] uniqueToReadErrorRateFormula;
	/** Linear regression coefficients for unique percentage vs base error rate */
	public double[] uniqueToBaseErrorRateFormula;
	
	/** Header parser for extracting coordinates from read identifiers */
	private IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();
	
}
