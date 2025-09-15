package hiseq;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import shared.Tools;
import structures.FloatList;
import structures.Point;

public class FlowCell {
	
	public FlowCell(String fname){
		TileDump.loadDump(fname, this);
	}
	
	public FlowCell(int k_){k=k_;}
	
	public MicroTile getMicroTile(String id) {//This method is NOT threadsafe
		ihp.parse(id);
		return getMicroTile(ihp);
	}
	
	public MicroTile getMicroTile(IlluminaHeaderParser2 ihp){
		return getLane(ihp.lane()).getMicroTile(ihp.tile(), ihp.xPos(), ihp.yPos(), true);
	}
	
	public MicroTile getMicroTile(int lane, int tile, int x, int y){
		return getLane(lane).getMicroTile(tile, x, y, true);
	}
	
	public MicroTile getMicroTile(int lane, int tile, int x, int y, boolean create){
		return getLane(lane).getMicroTile(tile, x, y, create);
	}
	
	public Lane getLane(int lane){
		while(lanes.size()<=lane){lanes.add(new Lane(lanes.size()));}
		return lanes.get(lane);
	}

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
			FlowCell fc2=fc.widen(true);
			fc2.calcStats();
			if(fc2.avgReads<=fc.avgReads){
				unwiden();
				return fc;
			}
			fc=fc2;
		}
		return fc.setFrom(this);
	}
	
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
			FlowCell fc2=fc.widen(false);
			fc2.calcStats();
			if(fc2.avgAlignedReads<=fc.avgAlignedReads){
				unwiden();
				return fc;
			}
			fc=fc2;
		}
		return fc.setFrom(this);
	}
	
	public void unwiden(){
		if(Tile.xSize>Tile.ySize){Tile.ySize/=2;}
		else{Tile.xSize/=2;}
	}
	
	public FlowCell widen(boolean loud){
		final int x=Tile.xSize, y=Tile.ySize;
		final int x2=x>=y ? x : 2*x;
		final int y2=x>=y ? 2*y : y;
		return widen(x2, y2, loud);
	}
	
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
	
	public void dump(String fname, boolean overwrite) {
		TileDump.write(this, fname, overwrite);
	}
	
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

	public double alignmentRate() {
		return readsAligned/(double)readsProcessed;
	}
	
	public double baseErrorRate() {
		return baseErrors/(double)(Tools.max(basesAligned, 1));
	}
	
	public ArrayList<Lane> lanes=new ArrayList<Lane>();
	
	long readsProcessed;
	long basesProcessed;
	long readsAligned;
	long basesAligned;
	long readErrors;
	long baseErrors;
	
	String name;

	int xMin=-1;
	int xMax=-1;
	int yMin=-1;
	int yMax=-1;
	int tMin=-1;
	int tMax=-1;
	int k=31;

	public double avgReads;
	public double avgAlignedReads;
	public double minCountToUse;
	
	public double avgQuality;
	public double avgUnique;
	public double avgDepth;
	public double avgErrorFree;
	public double avgPolyG;
	public double avgG;
	
	public double stdQuality;
	public double stdUnique;
	public double stdDepth;
	public double stdErrorFree;
	public double stdPolyG;
	public double stdG;
	
	public double[] uniqueToReadErrorRateFormula;
	public double[] uniqueToBaseErrorRateFormula;
	
	private IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();
	
}
