package hiseq;

import java.util.ArrayList;
import java.util.Iterator;

import structures.ByteBuilder;

public class Tile implements Iterable<MicroTile> {
	
	public Tile(int lane_, int tile_){
		lane=lane_;
		tile=tile_;
	}
	
	public MicroTile get(int x, int y, boolean create){
		final int xindex=x/xSize, yindex=y/ySize;
		ArrayList<MicroTile> ylist=getIndex(xindex, create);
		if(ylist==null || (yindex>=ylist.size() && !create)) {return null;}
		while(yindex>=ylist.size()){ylist.add(null);}
		MicroTile mt=ylist.get(yindex);
		if(mt==null && create){
			mt=new MicroTile(lane, tile, xindex*xSize, (xindex+1)*xSize-1, yindex*ySize, (yindex+1)*ySize-1);
			ylist.set(yindex, mt);
		}
		assert(mt==null || mt.contains(x,  y)) : x+", "+y+", "+xindex+", "+yindex+", "+mt;
		return mt;
	}
	
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
	
	public void add(Tile tb) {
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
	
	public ArrayList<ArrayList<MicroTile>> xlist=new ArrayList<ArrayList<MicroTile>>();
	
	public int lane;
	public int tile;
	public static int xSize=500;
	public static int ySize=500;
}
