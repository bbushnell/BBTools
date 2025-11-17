package prok;

import java.util.ArrayList;

import fileIO.ByteFile;
import shared.Parse;
import shared.Tools;

/**
 * Parses gene model files containing pre-computed k-mer statistics for prokaryotic gene prediction.
 * Loads statistical models for CDS, tRNA, 16S, 23S, 5S, and 18S sequence types with frame-specific
 * k-mer frequency data used by prokaryotic gene annotation tools.
 * @author Brian Bushnell
 */
public class GeneModelParser {
	
	/**
	 * Constructs a parser for the specified gene model file.
	 * Loads all lines from the file and initializes an empty GeneModel for population.
	 * @param fname_ Path to the gene model file to parse
	 */
	GeneModelParser(String fname_){
		fname=fname_;
		lines=ByteFile.toLines(fname);
		gm=new GeneModel(false);
	}
	
	/** Checks if more lines are available for parsing.
	 * @return true if there are unprocessed lines remaining in the file */
	boolean hasMore(){
		return pos<lines.size();
	}
	
	/** Returns the next line from the file and advances the position cursor.
	 * @return The next line as a byte array, or null if at end of file */
	byte[] nextLine(){
		if(pos>=lines.size()){return null;}
		byte[] line=lines.get(pos);
		pos++;
		return line;
	}
	
	/** Path to the gene model file being parsed */
	final String fname;
	/** All lines from the gene model file stored as byte arrays */
	final ArrayList<byte[]> lines;
	/** The GeneModel object being populated during parsing */
	private final GeneModel gm;
	/** Current position in the lines array during parsing */
	int pos=0;
	
	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads a complete gene model from the specified file.
	 * Creates a parser instance and processes the entire file to build the model.
	 * @param fname Path to the gene model file
	 * @return A fully populated GeneModel with all statistics containers
	 */
	public static GeneModel loadModel(String fname) {
		GeneModelParser gmp=new GeneModelParser(fname);
		return gmp.parse();
	}
	
	/**
	 * Parses the entire gene model file and populates the GeneModel object.
	 * First processes header information, then parses 6 StatsContainer objects
	 * for different sequence types (CDS, tRNA, 16S, 23S, 5S, 18S).
	 * @return The fully populated GeneModel with all statistical data
	 */
	private GeneModel parse(){
		while(hasMore()){
			byte[] line=nextLine();
			boolean valid=parseHeader(line);
			if(!valid){
				pos--;
				break;
			}
		}//Done parsing headers
		
		ArrayList<StatsContainer> containers=new ArrayList<StatsContainer>();
		while(hasMore()){
			StatsContainer sc=parseContainer();
			if(sc!=null){
				containers.add(sc);
			}else{
				assert(false);
			}
		}
		
		assert(containers.size()==6) : containers.size();
		for(StatsContainer sc : containers){
			gm.allContainers[sc.type].setFrom(sc);
		}
//		gm.statsCDS.setFrom(containers.get(0);
//		gm.statstRNA=containers.get(1);
//		gm.stats16S=containers.get(2);
//		gm.stats23S=containers.get(3);
//		gm.stats5S=containers.get(4);
//		gm.stats18S=containers.get(5);

//		gm.stats16S.minIdentity=ProkObject.min16SIdentity;
//		gm.stats23S.minIdentity=ProkObject.min23SIdentity;
//		gm.stats18S.minIdentity=ProkObject.min18SIdentity;
		
		gm.setStatics();
		
		return gm;
	}
	
	/**
	 * Parses a single StatsContainer section from the gene model file.
	 * Reads container metadata (name, type, length statistics) followed by
	 * three FrameStats objects for inner, start, and stop regions.
	 * @return A populated StatsContainer with all frame statistics
	 */
	private StatsContainer parseContainer(){
		String name=null;
		int type=-1;
		long lengthCount=0;
		long lengthSum=0;
		for(byte[] line=nextLine(); line!=null; line=nextLine()){
			if(line[0]!='#'){
				pos--;
				break;
			}
			
			if(Tools.startsWith(line, "##")){
				//ignore
			}else if(Tools.startsWith(line, "#name")){
				name=parseString(line);
			}else if(Tools.startsWith(line, "#type")){
				type=parseInt(line);
			}else if(Tools.startsWith(line, "#count")){
				lengthCount=parseLong(line);
			}else if(Tools.startsWith(line, "#lengthSum")){
				lengthSum=parseLong(line);
			}else if(Tools.startsWith(line, "#contains")){
				break;
			}else{
				assert(false) : new String(line);
			}
		}
		
		ArrayList<FrameStats> list=new ArrayList<FrameStats>(3);
		for(int i=0; i<3; i++){
			FrameStats fs=parseStats();
			list.add(fs);
		}
		
		StatsContainer sc=new StatsContainer(type);
		assert(sc.name.equals(name)) : name+", "+sc.name+", "+type;
		sc.lengthCount=lengthCount;
		sc.lengthSum=lengthSum;
		
		sc.setInner(list.get(0));
		sc.setStart(list.get(1));
		sc.setStop(list.get(2));
		
		sc.calculate();
		assert(sc.inner!=null);
		return sc;
	}
	
	/**
	 * Parses a single FrameStats section containing k-mer frequency data.
	 * Reads metadata (name, k-mer length, frame count, offset) followed by
	 * frequency data for each frame in both forward and reverse orientations.
	 * @return A populated FrameStats object with k-mer frequency tables
	 */
	private FrameStats parseStats(){
		String name=null;
		int k=-1, frames=-1, offset=-1;
//		System.err.println("A");
		for(byte[] line=nextLine(); line!=null; line=nextLine()){
			if(line[0]!='#'){
				pos--;
//				System.err.println("B");
				assert(false) : new String(line);
				break;
			}
			
			if(Tools.startsWith(line, "##")){
				//ignore
			}else if(Tools.startsWith(line, "#name")){
				name=parseString(line);
			}else if(Tools.startsWith(line, "#k")){
				k=parseInt(line);
			}else if(Tools.startsWith(line, "#frames")){
				frames=parseInt(line);
			}else if(Tools.startsWith(line, "#offset")){
				offset=parseInt(line);
			}else if(Tools.startsWith(line, "#valid\tframe")){
//				assert(false);
//				System.err.println("C");
				break;
			}
//			System.err.println("D");
		}
//		assert(false);
//		System.err.println("E");
		
		FrameStats fs=new FrameStats(name, k, frames, offset);
		
		for(int i=0, max=2*fs.frames; i<max; i++){
			byte[] line=nextLine();
			fs.parseData(line);
		}
		return fs;
	}
	
	/**
	 * Extracts a string value from a tab-delimited line.
	 * Returns the portion after the first tab character.
	 * @param line The line to parse as byte array
	 * @return The string value after the tab separator
	 */
	private static String parseString(byte[] line){
		int idx=Tools.indexOf(line, '\t');
		String s=new String(line, idx+1, line.length-idx-1);
		return s;
	}
	/**
	 * Extracts an integer value from a tab-delimited line.
	 * Parses the portion after the first tab character as an integer.
	 * @param line The line to parse as byte array
	 * @return The integer value after the tab separator
	 */
	private static int parseInt(byte[] line){
		int idx=Tools.indexOf(line, '\t');
		return Parse.parseInt(line, idx+1, line.length);
	}
	/**
	 * Extracts a long value from a tab-delimited line.
	 * Parses the portion after the first tab character as a long.
	 * @param line The line to parse as byte array
	 * @return The long value after the tab separator
	 */
	private static long parseLong(byte[] line){
		int idx=Tools.indexOf(line, '\t');
		return Parse.parseLong(line, idx+1, line.length);
	}
	
//	public static void parseHeaderStatic(byte[] line){
//		
//		assert(line[0]=='#');
//		if(Tools.startsWith(line, "#k_inner")){
//			int x=(int)parseLong(line);
//			assert(x==innerKmerLength);
//			setInnerK(x);
//		}else if(Tools.startsWith(line, "#k_end")){
//			int x=(int)parseLong(line);
//			assert(x==endKmerLength);
//			setEndK(x);
//		}else if(Tools.startsWith(line, "#start_left_offset")){
//			int x=(int)parseLong(line);
//			assert(x==startLeftOffset);
//			setStartLeftOffset(x);
//		}else if(Tools.startsWith(line, "#start_right_offset")){
//			int x=(int)parseLong(line);
//			assert(x==startRightOffset);
//			setStartRightOffset(x);
//		}else if(Tools.startsWith(line, "#stop_left_offset")){
//			int x=(int)parseLong(line);
//			assert(x==stopLeftOffset);
//			setStopLeftOffset(x);
//		}else if(Tools.startsWith(line, "#stop_right_offset")){
//			int x=(int)parseLong(line);
//			assert(x==stopRightOffset);
//			setStopRightOffset(x);
//		}
//	}
	
	/**
	 * Parses header lines that contain file metadata and global statistics.
	 * Processes lines starting with '#' to extract file counts, tax IDs,
	 * scaffold counts, base counts, and gene counts into the GeneModel.
	 *
	 * @param line The header line to parse as byte array
	 * @return true if the line was a valid header, false if parsing should continue elsewhere
	 */
	public boolean parseHeader(byte[] line){
		if(line[0]!='#'){return false;}
		
		if(Tools.startsWith(line, "#BBMap")){
			//ignore
		}else if(Tools.startsWith(line, "##")){
			//ignore
		}else if(Tools.startsWith(line, "#files")){//Not necessary
			String[] split=new String(line).split("\t");
			try {
				gm.numFiles+=Integer.parseInt(split[1]);
			} catch (NumberFormatException e) {
				gm.numFiles+=split.length-1;//old style pgm
			}
//			for(String s : new String(line).split("\t")){
//				if(s.charAt(0)!='#'){
//					gm.fnames.add(s);
//				}
//			}
		}else if(Tools.startsWith(line, "#taxIDs")){//Can be made faster
			for(String s : new String(line).split("\t")){
				if(s.charAt(0)!='#'){
					gm.taxIds.add(Integer.parseInt(s));
				}
			}
		}else if(Tools.startsWith(line, "#scaffolds")){
			long x=parseLong(line);
			gm.readsProcessed=x;
		}else if(Tools.startsWith(line, "#bases")){
			long x=parseLong(line);
			gm.basesProcessed=x;
		}else if(Tools.startsWith(line, "#genes")){
			long x=parseLong(line);
			gm.genesProcessed=x;
		}else if(Tools.startsWith(line, "#GC")){
			//ignore
		}else if(Tools.startsWith(line, "#ACGTN")){
			String[] split=new String(line).split("\t");
			for(int i=0; i<gm.baseCounts.length; i++){
				gm.baseCounts[i]=Long.parseLong(split[i+1]);
			}
		}else{
			return false;
		}
		return true;
	}
	
}
