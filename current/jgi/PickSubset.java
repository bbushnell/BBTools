package jgi;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.LineParser1;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Select a maximally diverse subset from an all-to-all ANI comparison,
 * by removing a node with the highest pairwise ANI until the
 * desired number of nodes remain, or the desired ANI threshold
 * is reached.
 * Unlike RepresentativeSet, it does not use taxonomy or centroids.
 * @author Brian Bushnell
 * @date July 30, 2024
 *
 */
public class PickSubset {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		PickSubset x=new PickSubset(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public PickSubset(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, /*getClass()*/null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
		}

		validateParams();
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
		ffoutInvalid=FileFormat.testOutput(outInvalid, FileFormat.TXT, null, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Set any necessary Parser defaults here
		//parser.foo=bar;
//		parser.out1="stdout";
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("invalid") || a.equals("outu")){
				outInvalid=b;
			}else if(a.equals("nodes") || a.equals("files") || a.equals("genomes")){
				desiredNodes=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("ani") || a.equalsIgnoreCase("maxani")){
				desiredANI=Float.parseFloat(b);
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, out1)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
//		if(!ByteFile.FORCE_MODE_BF2){
//			ByteFile.FORCE_MODE_BF2=false;
//			ByteFile.FORCE_MODE_BF1=true;
//		}
	}
	
	/** Ensure parameter ranges are within bounds and required parameters are set */
	private boolean validateParams(){
//		assert(minfoo>0 && minfoo<=maxfoo) : minfoo+", "+maxfoo;
		assert(desiredANI>0 || desiredNodes>0) : "Desired ANI or genomes must be greater than 0.";
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create streams and process all data */
	void process(Timer t){
		
		ByteFile bf=ByteFile.makeByteFile(ffin1);
		ByteStreamWriter bsw=makeBSW(ffout1);
		ByteStreamWriter bswInvalid=makeBSW(ffoutInvalid);
		
//		assert(false) : "Header goes here.";
		if(bsw!=null){
//			assert(false) : "Header goes here.";
		}
		
		ArrayList<Edge> edges=readFile(bf);
		ArrayList<Node> chosen=pickNodes(edges, desiredNodes);
		if(bsw!=null) {
			for(Node n : chosen) {
				bsw.println(n.name);
			}
		}
		if(bswInvalid!=null) {
			for(Entry<String, Node> e : nodeMap.entrySet()) {
				Node n=e.getValue();
				if(!n.valid) {
					bswInvalid.println(n.name);
				}
			}
		}
		
		float maxANI=0;
		for(Node n : chosen) {
			for(Edge e : n.edges) {
				Node r=nodeMap.get(e.ref);
				if(r.valid) {
					maxANI=Tools.max(maxANI, e.ani);
				}
			}
		}
		
		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		if(bswInvalid!=null){errorState|=bswInvalid.poisonAndWait();}
		
		t.stop();
		
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
		outstream.println();
		outstream.println("Nodes In:          \t"+nodeMap.size());
		outstream.println("Nodes Out:         \t"+chosen.size());
		outstream.println("Maximum ANI:       \t"+maxANI);
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private ArrayList<Edge> readFile(ByteFile bf){
		byte[] line=bf.nextLine();
		ByteBuilder bb=new ByteBuilder();
		
		ArrayList<Edge> edges=new ArrayList<Edge>();
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				
				if(line[0]!='#'){
					Edge e=new Edge(line);
					Node query=nodeMap.get(e.query);
					if(query==null) {
						query=new Node(line);
						nodeMap.put(e.query, query);
					}
					if(!e.query.equals(e.ref)) {//if not a self edge
						query.edges.add(e);
						edges.add(e);
					}
				}
			}
			line=bf.nextLine();
		}
		return edges;
	}

	ArrayList<Node> pickNodes(ArrayList<Edge> list, int desiredNodes){
		int remaining=nodeMap.size();
		for(Node n : nodeMap.values()) {n.calcSum();}
		Collections.sort(list);
		Collections.reverse(list);
		for(Edge e : list) {
			if(remaining<=desiredNodes || e.ani<=desiredANI) {break;}
			if(e.valid) {
				Node q=nodeMap.get(e.query);
				Node r=nodeMap.get(e.ref);
				if(q.valid && r.valid) {
					if(q.compareTo(r)>=0) {
						q.setInvalid();
					}else {
						r.setInvalid();
					}
					remaining--;
				}
				e.valid=false;
			}
		}
		ArrayList<Node> chosen=new ArrayList<Node>(remaining);
		for(Node n : nodeMap.values()) {
			if(n.valid) {chosen.add(n);}
		}
		return chosen;
	}
	
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	private HashMap<String, Node> nodeMap=new HashMap<String, Node>();
	
	/*--------------------------------------------------------------*/
	
	private class Node implements Comparable<Node> {
		
		Node(byte[] line){
			lp.set(line);
			name=lp.parseString(0);
			size=lp.parseLong(3);
			bases=lp.parseLong(5);
		}
		
		float calcSum() {
			sum=0;
			for(Edge e : edges) {
				sum+=e.kid;
				sum+=e.ani*0.0001;//In case kid is absent
			}
			return sum;
		}
		
		void setInvalid() {
			assert(valid);
			valid=false;
			for(Edge e : edges) {e.valid=false;}
		}
		
		@Override
		public int compareTo(Node o) {
			if(sum!=o.sum) {return sum>o.sum ? 1 : -1;}
			return name.compareTo(o.name);
		}
		
		
		final String name;
		long size;
		long bases;
		float sum;
		boolean valid=true;
		ArrayList<Edge> edges=new ArrayList<Edge>();
	}
	
//	#Query	Ref	ANI	QSize	RefSize	QBases	RBases	QTaxID	RTaxID	KID	WKID	SSU
//	a.fa	b.fa	80.980	11178684	8473676	12055201	8514625	-1	-1	0.237	0.313	94.233
	
	private class Edge implements Comparable<Edge> {
		
		Edge(byte[] line){
			lp.set(line);
			query=lp.parseString(0);
			ref=lp.parseString(1);
			ani=lp.parseFloat(2);
			kid=(lp.terms()>9) ? lp.parseFloat(9) : 0;
			if(lp.terms()<11 || lp.termEquals('.', 11)) {
				ssu=0;
			}else {
				ssu=lp.parseFloat(11);
			}
		}
		
		@Override
		public int compareTo(Edge o) {
			if(ani!=o.ani) {return ani>o.ani ? 1 : -1;}
			if(kid!=o.kid) {return kid>o.kid ? 1 : -1;}
			int x=query.compareTo(o.query);
			if(x!=0) {return x;}
			return ref.compareTo(o.ref);
		}
		
		final String query;
		final String ref;
		float ani;
		float kid;
		float ssu;
		boolean valid=true;
	}
	
	private LineParser1 lp=new LineParser1('\t');
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;

	/** Primary output file path */
	private String out1=null;

	/** Junk output file path */
	private String outInvalid=null;
	
	private int desiredNodes=0;
	private float desiredANI=0;
	
	/*--------------------------------------------------------------*/
	
	private long edgesProcessed=0;
	private long nodesProcessed=0;
	private long edgesOut=0;
	private long nodesOut=0;
	
	private long linesProcessed=0;
	private long linesOut=0;
	private long bytesProcessed=0;
	private long bytesOut=0;
	
	private long maxLines=Long.MAX_VALUE;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Input File */
	private final FileFormat ffin1;
	/** Output File */
	private final FileFormat ffout1;
	/** Optional Output File for Junk */
	private final FileFormat ffoutInvalid;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	
}
