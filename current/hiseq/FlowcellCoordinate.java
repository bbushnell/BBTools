package hiseq;

/**
 * Parses and represents coordinate information for Illumina sequencing
 * flowcell data points. Extracts and manages lane, tile, x, and y
 * coordinates from Illumina sequencing read identifiers using an
 * IlluminaHeaderParser2. Provides thread-safe coordinate parsing and
 * comparison capabilities for sequencing data processing.
 *
 * @author Brian Bushnell
 */
public class FlowcellCoordinate implements Comparable<FlowcellCoordinate> {
	
	/** Creates a FlowcellCoordinate with default uninitialized values. */
	public FlowcellCoordinate() {}
	
	/** Creates a FlowcellCoordinate by parsing the specified identifier.
	 * @param id Illumina sequencing read identifier to parse */
	public FlowcellCoordinate(String id) {
		setFrom(id);
	}
	
//	public float distance(FlowcellCoordinate fc){ //Comment out due to being unused
//		assert(isSet());
//		assert(fc.isSet());
//
//		if(lane!=fc.lane){return big;}
//
//		long a=Tools.absdif(x, fc.x), b=Tools.absdif(y, fc.y);
//		if(tile!=fc.tile){
//			return spanTiles ? Tools.min(a, b) : big;
//		}
//		return (float)Math.sqrt(a*a+b*b);
//
//		//Hard to say...  could consider adjacent tiles?
////		if(tile!=fc.tile){
////			if(allowAdjacentTiles && Tools.absdif(tile, fc.tile)<2){return Tools.min(x-fc.x, y-fc.y);}
////			return big;
////		}
////
////		long a=x-fc.x, b=y-fc.y;
////		return (float)Math.sqrt(a*a+b*b);
//	}

	//2402:6:1101:6337:2237/1
	//MISEQ08:172:000000000-ABYD0:1:1101:18147:1925 1:N:0:TGGATATGCGCCAATT
	//HISEQ07:419:HBFNEADXX:1:1101:1238:2072
//	public void setFrom_old(String id){
//		final int lim=id.length();
//		
//		int i=0;
//		int current=0;
//		while(i<lim && id.charAt(i)!=' ' && id.charAt(i)!='/'){i++;}
//		if(i>=lim){i--;}
//		for(int semis=0; i>=0; i--){
//			if(id.charAt(i)==':'){
//				semis++;
//				if(semis==4){break;}
//			}
//		}
//		i++;
//		
//		assert(Tools.isDigit(id.charAt(i))) : id;
//		while(i<lim && Tools.isDigit(id.charAt(i))){
//			current=current*10+(id.charAt(i)-'0');
//			i++;
//		}
//		lane=current;
//		current=0;
//		i++;
//		
//		if(!Tools.isDigit(id.charAt(i))){//Hiseq 3000?
//			while(i<lim && id.charAt(i)!=':'){i++;}
//			i++;
//
//			assert(Tools.isDigit(id.charAt(i))) : id;
//			while(i<lim && Tools.isDigit(id.charAt(i))){
//				current=current*10+(id.charAt(i)-'0');
//				i++;
//			}
//			lane=current;
//			current=0;
//			i++;
//		}
//
//		assert(Tools.isDigit(id.charAt(i))) : id;
//		while(i<lim && Tools.isDigit(id.charAt(i))){
//			current=current*10+(id.charAt(i)-'0');
//			i++;
//		}
//		tile=current;
//		current=0;
//		i++;
//
//		assert(Tools.isDigit(id.charAt(i))) : id;
//		while(i<lim && Tools.isDigit(id.charAt(i))){
//			current=current*10+(id.charAt(i)-'0');
//			i++;
//		}
//		x=current;
//		current=0;
//		i++;
//
//		assert(Tools.isDigit(id.charAt(i))) : id;
//		while(i<lim && Tools.isDigit(id.charAt(i))){
//			current=current*10+(id.charAt(i)-'0');
//			i++;
//		}
//		y=current;
//		current=0;
//		i++;
//	}
	
	/**
	 * Parses flowcell coordinate information from an Illumina read identifier.
	 * Extracts lane, tile, x, y coordinates and optionally UMI from the header.
	 * @param id Illumina sequencing read identifier string
	 */
	public void setFrom(String id){
		ihp.parse(id);
		
		lane=ihp.lane();
		tile=ihp.tile();
		x=ihp.xPos();
		y=ihp.yPos();
		if(parseUMI) {umi=ihp.index3();}
	}
	
	/** Checks if all coordinate fields have valid values.
	 * @return true if lane, tile, x, and y are all non-negative */
	public boolean isSet(){
		return lane>=0 && tile>=0 && x>=0 && y>=0;
	}

	@Override
	public int compareTo(FlowcellCoordinate b) {
		if(lane!=b.lane){return lane-b.lane;}
		if(tile!=b.tile){return tile-b.tile;}
		if(y!=b.y){return y-b.y;}
		if(x!=b.x){return x-b.x;}
		return 0;
	}
	
	/** Lane number from the flowcell, initialized to -1 */
	public int lane=-1;
	/** Tile number within the lane, initialized to -1 */
	public int tile=-1;
	/** X coordinate on the tile, initialized to -1 */
	public int x=-1;
	/** Y coordinate on the tile, initialized to -1 */
	public int y=-1;
	/** Unique Molecular Identifier, parsed if parseUMI is enabled */
	public String umi;
	/** Parser for extracting coordinate information from Illumina headers */
	private final IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();
	
	/** Large distance value for coordinate calculations */
	public static final float big=10000000;
	/** Flag to enable UMI parsing from read identifiers */
	public static boolean parseUMI=false;
//	public static boolean spanTiles=false;

	
	/**
	 * Thread-safe retrieval of a FlowcellCoordinate instance.
	 * Uses ThreadLocal storage to provide one instance per thread.
	 * @return FlowcellCoordinate instance for the current thread
	 */
	public static FlowcellCoordinate getFC(){
		FlowcellCoordinate fc=localFC.get();
		if(fc==null){
			fc=new FlowcellCoordinate();
			localFC.set(fc);
		}
		return fc;
	}
	/** Thread-local storage for FlowcellCoordinate instances */
	private static final ThreadLocal<FlowcellCoordinate> localFC=new ThreadLocal<FlowcellCoordinate>();
	
}
