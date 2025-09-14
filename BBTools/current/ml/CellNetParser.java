package ml;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import shared.LineParser2;
import shared.Parse;
import shared.Tools;
import structures.FloatList;
import structures.IntList;

public class CellNetParser {
	
	static CellNet parse(String fname, boolean nullOnFailure) {
		if(nullOnFailure && !new File(fname).exists()) {return null;}
		CellNetParser cnp=new CellNetParser(fname);
		return cnp.net;
	}

	public static CellNet load(String fname) {return load(fname, true);}
	public static CellNet load(String fname, boolean nullOnFailure) {
		return parse(fname, nullOnFailure);
	}
	
	private CellNetParser(String fname_){
		this(ByteFile.toLines(fname_));
		fname=fname_;
	}
	
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
	
	private void parseEdgesDense() {
		assert(concise);
		pos=posFirstEdge;
		long numEdges=0;
		
		int numCells=(int) shared.Vector.sum(dims);
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
	
	private void parseEdgesSparse() {
		assert(concise);
		pos=posFirstEdge;
		long numEdges=0;
		
		int numCells=(int) shared.Vector.sum(dims);
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
	
	boolean hasMore(){
		return pos<lines.size();
	}
	
	byte[] nextLine(){
		if(pos>=lines.size()){return null;}
		byte[] line=lines.get(pos);
		pos++;
		return line;
	}
	
	String fname;
	final ArrayList<byte[]> lines;
	private final CellNet net;
	long seed;
//	long annealSeed=-1;
	float density=1f;
	float density1=0f;
	int edgeBlockSize=1;
	int edges=0;//unused
	long epochs=0;
	long samples=0;
	int version;
	int layers=-1;
	boolean concise=false;
	boolean dense=true;
	int[] dims;
	int pos=0;
	float cutoff=0.5f;
	final int posFirstEdge;
	ArrayList<String> commands=new ArrayList<String>();
	
	public static final byte delimiter=' ';
	
	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/
	
	private static String parseString(byte[] line){
		int idx=Tools.indexOf(line, delimiter);
		String s=new String(line, idx+1, line.length-idx-1);
		return s;
	}
	
	private static float parseFloat(byte[] line){
		int idx=Tools.indexOf(line, delimiter);
		return Parse.parseFloat(line, idx+1, line.length);
	}
	
	private static int parseInt(byte[] line){
		int idx=Tools.indexOf(line, delimiter);
		return Parse.parseInt(line, idx+1, line.length);
	}
	
	private static long parseLong(byte[] line){
		int idx=Tools.indexOf(line, delimiter);
		return Parse.parseLong(line, idx+1, line.length);
	}
	
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
