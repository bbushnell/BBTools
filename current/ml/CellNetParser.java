package ml;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import parse.LineParser2;
import parse.Parse;
import shared.Tools;
import structures.FloatList;
import structures.IntList;

/**
 * Parser for loading CellNet neural network configurations from text files.
 * Supports both dense and sparse weight representations with various metadata.
 * @author Brian Bushnell
 */
public class CellNetParser {
	
	/**
	 * Parses a CellNet configuration from the specified file.
	 * @param fname Path to the CellNet configuration file
	 * @param nullOnFailure If true, returns null when file doesn't exist instead of throwing exception
	 * @return Parsed CellNet object or null if file doesn't exist and nullOnFailure is true
	 */
	static CellNet parse(String fname, boolean nullOnFailure) {
		if(nullOnFailure && !new File(fname).exists()) {return null;}
		CellNetParser cnp=new CellNetParser(fname);
		return cnp.net;
	}

	/**
	 * Loads a CellNet from the specified file.
	 * @param fname Path to the CellNet configuration file
	 * @return Loaded CellNet object
	 */
	public static CellNet load(String fname) {return load(fname, true);}
	/**
	 * Loads a CellNet from the specified file with optional null handling.
	 * @param fname Path to the CellNet configuration file
	 * @param nullOnFailure If true, returns null when file doesn't exist instead of throwing exception
	 * @return Loaded CellNet object or null if file doesn't exist and nullOnFailure is true
	 */
	public static CellNet load(String fname, boolean nullOnFailure) {
		return parse(fname, nullOnFailure);
	}
	
	/** Constructs a parser for the specified CellNet file.
	 * @param fname_ Path to the CellNet configuration file */
	private CellNetParser(String fname_){
		this(ByteFile.toLines(fname_));
		fname=fname_;
	}
	
	/**
	 * Constructs a parser from pre-loaded file lines.
	 * Parses header information and creates the CellNet with appropriate edge parsing.
	 * @param lines_ List of file lines as byte arrays
	 */
	private CellNetParser(ArrayList<byte[]> lines_){
		lines=lines_;
		
		parseHeader();
		CellNet.DENSE=dense;
		
		net=new CellNet(dims, seed, density, density1, edgeBlockSize, commands);
		net.epochsTrained=epochs;
		net.samplesTrained=samples;
//		net.annealSeed=annealSeed;
		net.setCutoff(cutoff);
		assert(layers==net.layers);
		posFirstEdge=pos;
		if(dense) {
			parseEdgesDense();
		}else {
			parseEdgesSparse();
		}
		net.makeWeightMatrices();
	}
	
	/** Parses header lines containing network metadata and configuration.
	 * Processes version, layers, seed, density, dimensions, and other parameters. */
	public void parseHeader() {
		
		//TODO: This should really use LineParser instead.
		while(pos<lines.size()) {
			byte[] line=lines.get(pos);
			
			if(line.length<1){
				//ignore
			}else if(Tools.startsWith(line, "#")){//header
				if(Tools.startsWith(line, "##ctf") || Tools.startsWith(line, "#ctf")){
					cutoff=parseFloat(line);
				}else if(Tools.startsWith(line, "##")){
					//Comment; ignore
				}else if(Tools.startsWith(line, "#version")){
					version=parseInt(line);
				}else if(Tools.startsWith(line, "#layers")){
					layers=parseInt(line);
				}else if(Tools.startsWith(line, "#seed")){
					seed=parseLong(line);
				}else if(Tools.startsWith(line, "#annealseed")){
//					annealSeed=parseLong(line);
				}else if(Tools.startsWith(line, "#density1")){
					density1=parseFloat(line);
				}else if(Tools.startsWith(line, "#density")){
					density=parseFloat(line);
				}else if(Tools.startsWith(line, "#blocksize")){
					edgeBlockSize=parseInt(line);
//				}else if(Tools.startsWith(line, "#edgecount")){
//					edges=parseInt(line);
				}else if(Tools.startsWith(line, "#epochs")){
					epochs=parseLong(line);
				}else if(Tools.startsWith(line, "#samples")){
					samples=parseLong(line);
				}else if(Tools.startsWith(line, "#concise")){
					concise=true;
				}else if(Tools.startsWith(line, "#dense")){
					dense=true;
				}else if(Tools.startsWith(line, "#sparse")){
					dense=false;
				}else if(Tools.startsWith(line, "#dims")){
					dims=parseIntArray(line, delimiter, true);
					assert(layers==dims.length) : layers+", "+Arrays.toString(dims);
				}else if(Tools.startsWith(line, "#CL")){
					commands.add(new String(line));
				}else if(Tools.startsWith(line, "#edges")){
					if(line.length>7) {edges=parseInt(line);}
				}else if(Tools.startsWith(line, "#")){
					assert(false) : "\nUnexpected header line: '"+new String(line)+"'"
							+ "\nComments should start with ##\n";
				}
			}else{
				break; //A cell or edge
			}
			pos++;
		}
//		assert(false) : pos+", "+new String(lines.get(pos));
	}
	
	/** Parses dense weight representation where all connections are stored.
	 * Each cell line contains bias and all weights in sequential order. */
	private void parseEdgesDense() {
		assert(concise);
		pos=posFirstEdge;
		long numEdges=0;
		
		int numCells=(int) simd.Vector.sum(dims);
		LineParser2 lp=new LineParser2(delimiter);
		FloatList weights=new FloatList();
		
		while(pos<lines.size()) {
			byte[] line=lines.get(pos);
			
			if(line.length<1){
				//ignore
			}else if(Tools.startsWith(line, "##")){
				//ignore
			}else if(Tools.startsWith(line, 'C') || Tools.startsWith(line, 'W')){
				
				lp.set(line);
				lp.setBounds(0, 0);
				int cid=lp.parseInt();
//				assert(false) : cid;
				String s=lp.parseString();
				int type=Tools.find(s, Function.TYPES);
				assert(type>=0) : type+", "+s+"\n'"+new String(line)+"'";
				Cell c=net.list.get(cid);
				c.function=Function.getFunction(type);
				assert(c.function.type()==type);
				
				c.setBias(lp.parseFloat(), true);
				weights.clear();
				while(lp.hasMore()) {
					weights.add(lp.parseFloat());
				}
				c.weights=weights.toArray();
				c.deltas=new float[c.weights.length];
				assert(c.weights.length==c.id()-c.lpos-c.prevLayerStart) : new String(line)+"\n"+
						c.weights.length+", "+c.layer+", "+c.id()+", "+c.lpos+", "+c.prevLayerStart+", "+
						(c.id()-c.lpos-c.prevLayerStart);
			}else {
				assert(false) : new String(line);
			}
			pos++;
		}
		assert(CellNet.DENSE || net.check());
	}
	
	/** Parses sparse weight representation with explicit input indices.
	 * Uses separate lines for cell weights (C/W) and input connections (I/H). */
	private void parseEdgesSparse() {
		assert(concise);
		pos=posFirstEdge;
		long numEdges=0;
		
		int numCells=(int) simd.Vector.sum(dims);
		LineParser2 lp=new LineParser2(delimiter);
		FloatList weights=new FloatList();
		IntList inputs=new IntList();
		
		while(pos<lines.size()) {
			byte[] line=lines.get(pos);
			
			if(line.length<1){
				//ignore
			}else if(Tools.startsWith(line, "##")){
				//ignore
			}else if(Tools.startsWith(line, 'C') || Tools.startsWith(line, 'W')){
				
				lp.set(line);
				lp.setBounds(0, 0);
				int cid=lp.parseInt();
//				assert(false) : cid;
				String s=lp.parseString();
				int type=Tools.find(s, Function.TYPES);
				assert(type>=0) : type+", "+s+"\n'"+new String(line)+"'";
				Cell c=net.list.get(cid);
				assert(c.weights==null);
				c.function=Function.getFunction(type);
				assert(c.function.type()==type);
				
				c.setBias(lp.parseFloat(), true);
				weights.clear();
				while(lp.hasMore()) {
					weights.add(lp.parseFloat());
				}
				c.weights=weights.toArray();
				c.deltas=new float[c.weights.length];
				assert(c.inputs==null || c.inputs.length==c.weights.length) : 
					c.layer+", "+c.lpos+", "+c.inputs.length+", "+c.weights.length+"\n"+Arrays.toString(c.inputs);
			}else if(Tools.startsWith(line, 'I')){
				
				lp.set(line);
				lp.setBounds(0, 0);
				int cid=lp.parseInt();
				Cell c=net.list.get(cid);
				assert(c.inputs==null);
				inputs.clear();
				while(lp.hasMore()) {
					inputs.add(lp.parseInt());
				}
				c.inputs=inputs.toArray();
//				for(int i=0; i<c.inputs.length; i++) {c.inputs[i]-=c.prevLayerStart;}
				assert(c.weights==null || c.inputs.length==c.weights.length);
			}else if(Tools.startsWith(line, 'H')){
				
				lp.set(line);
				lp.setBounds(0, 0);
				int cid=lp.parseInt();
				Cell c=net.list.get(cid);
				assert(c.inputs==null);
				lp.advance();
				c.inputs=CellNet.fromHex(line, lp.a());
				assert(c.weights==null || c.inputs.length==c.weights.length);
			}else {
				assert(false) : new String(line);
			}
			pos++;
		}
		CellNet.makeOutputSets(net.net);
		assert(net.check());
	}
	
	/** Checks if more lines are available for parsing.
	 * @return true if more lines exist, false otherwise */
	boolean hasMore(){
		return pos<lines.size();
	}
	
	/** Returns the next line for parsing and advances the position.
	 * @return Next line as byte array or null if no more lines */
	byte[] nextLine(){
		if(pos>=lines.size()){return null;}
		byte[] line=lines.get(pos);
		pos++;
		return line;
	}
	
	/** Source filename for the parsed CellNet */
	String fname;
	/** Lines from the CellNet configuration file */
	final ArrayList<byte[]> lines;
	/** The parsed CellNet neural network */
	private final CellNet net;
	/** Random seed used for network initialization */
	long seed;
//	long annealSeed=-1;
	/** Network connection density (1.0 = fully connected) */
	float density=1f;
	/** First layer connection density */
	float density1=0f;
	/** Block size for edge processing */
	int edgeBlockSize=1;
	/** Total edge count (unused) */
	int edges=0;//unused
	/** Number of training epochs completed */
	long epochs=0;
	/** Number of training samples processed */
	long samples=0;
	/** File format version number */
	int version;
	/** Number of network layers */
	int layers=-1;
	/** Whether the format uses concise representation */
	boolean concise=false;
	/** Whether to use dense weight storage (true) or sparse (false) */
	boolean dense=true;
	/** Dimensions (neuron counts) for each layer */
	int[] dims;
	/** Current parsing position in the line list */
	int pos=0;
	/** Classification threshold cutoff value */
	float cutoff=0.5f;
	/** Position of the first edge/weight line in the file */
	final int posFirstEdge;
	/** Command lines used to generate this network */
	ArrayList<String> commands=new ArrayList<String>();
	
	/** Delimiter character used for parsing (space) */
	public static final byte delimiter=' ';
	
	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Extracts string value after the first delimiter in a line.
	 * @param line Line to parse
	 * @return String portion after delimiter
	 */
	private static String parseString(byte[] line){
		int idx=Tools.indexOf(line, delimiter);
		String s=new String(line, idx+1, line.length-idx-1);
		return s;
	}
	
	/**
	 * Extracts float value after the first delimiter in a line.
	 * @param line Line to parse
	 * @return Float value after delimiter
	 */
	private static float parseFloat(byte[] line){
		int idx=Tools.indexOf(line, delimiter);
		return Parse.parseFloat(line, idx+1, line.length);
	}
	
	/**
	 * Extracts integer value after the first delimiter in a line.
	 * @param line Line to parse
	 * @return Integer value after delimiter
	 */
	private static int parseInt(byte[] line){
		int idx=Tools.indexOf(line, delimiter);
		return Parse.parseInt(line, idx+1, line.length);
	}
	
	/**
	 * Extracts long value after the first delimiter in a line.
	 * @param line Line to parse
	 * @return Long value after delimiter
	 */
	private static long parseLong(byte[] line){
		int idx=Tools.indexOf(line, delimiter);
		return Parse.parseLong(line, idx+1, line.length);
	}
	
	/**
	 * Parses a delimited line into an array of integers.
	 *
	 * @param line Line containing delimited integer values
	 * @param delimiter Character separating values
	 * @param parseTitle If true, skips the first field as a title
	 * @return Array of parsed integer values
	 */
	public static int[] parseIntArray(final byte[] line, final byte delimiter, boolean parseTitle){
		int a=0, b=0;
		IntList list=new IntList(3);
		
		if(parseTitle) {
			while(b<line.length && line[b]!=delimiter){b++;}
			assert(b>a) : "Missing Title: "+new String(line);
			b++;
			a=b;
		}
		
		while(a<line.length) {
			while(b<line.length && line[b]!=delimiter){b++;}
			assert(b>a) : "Missing element "+list.size+": '"+new String(line)+"'";
			int x=Parse.parseInt(line, a, b);
//			assert(x>0) : new String(line);
			list.add(x);
			b++;
			a=b;
		}
		return list.toArray();
	}
	
	
}
