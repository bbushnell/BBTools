package ml;

import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.Vector;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import structures.LongList;
import tracker.PalindromeTracker;

/**
 * @author Brian Bushnell
 * @date Oct 6, 2014
 *
 */
public class ScoreSequence {

	/** Program entry point for sequence scoring pipeline.
	 * @param args Command-line arguments specifying input files, neural network, and scoring parameters */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		ScoreSequence x=new ScoreSequence(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs ScoreSequence instance by parsing command-line arguments.
	 * Initializes file formats, loads neural network, and configures scoring parameters.
	 * @param args Command-line arguments containing input/output paths and scoring options
	 */
	public ScoreSequence(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("net") || a.equals("nn")){
				netFile=b;
			}else if(a.equals("hist")){
				histFile=b;
			}else if(a.equals("width")){
				width=Integer.parseInt(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(a.equals("rcomp")){
				rcomp=Parse.parseBoolean(b);
			}else if(a.equals("parse")){
				parseHeader=Parse.parseBoolean(b);
			}else if(a.equals("cutoff")){
				cutoff=Float.parseFloat(b);
				filter=true;
			}else if(a.equals("highpass")){
				highpass=Parse.parseBoolean(b);
				filter=true;
			}else if(a.equals("lowpass")){
				highpass=!Parse.parseBoolean(b);
				filter=true;
			}else if(a.equals("filter")){
				filter=Parse.parseBoolean(b);
			}else if(a.equals("annotate") || a.equals("rename")){
				annotate=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				//				throw new RuntimeException("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			in1=parser.in1;
			out1=parser.out1;
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		ffnet=FileFormat.testOutput(netFile, FileFormat.BBNET, null, true, true, false, false);
		net=CellNetParser.load(netFile);
		assert(net!=null) : netFile;
		if(width<0) {width=(net.numInputs()-4)/4;}
		else {assert(width==(net.numInputs()-4)/4) : width+", "+net.numInputs();}
	}
	
	/**
	 * Main processing loop that scores sequences using the neural network.
	 * Reads input sequences, applies scoring model, optionally filters results,
	 * and generates output files and histograms.
	 * @param t Timer for tracking execution performance
	 */
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			cris.start();
		}

		int inputs=width*4+4;
		final ByteStreamWriter bsw=(ffout1==null ? null : new ByteStreamWriter(ffout1));
		if(bsw!=null) {
			bsw.start();
//			bsw.print("#dims\t"+inputs+"\t1\n");
		}
		
		final float[] vec=new float[width*4+4];
		long readsProcessed=0, basesProcessed=0;
		final ByteBuilder bb=new ByteBuilder();
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					readsProcessed+=r1.pairCount();
					basesProcessed+=r1.pairLength();

//					float result=(parseHeader ? Parse.parseFloat(r1.id, "result=", '\t') : 0);
					processRead(r1, bb, bsw, vec);
					processRead(r2, bb, bsw, vec);
				}

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		errorState=ReadWrite.closeStreams(cris) | errorState;
		if(verbose){outstream.println("Finished reading data.");}

		if(bsw!=null) {errorState=bsw.poisonAndWait() | errorState;}
		
		if(histFile!=null) {
			bb.clear();
			PalindromeTracker.append(bb, "#count\tneg\tpos", new LongList[] {mhist, phist}, 101);
			ReadWrite.writeStringInThread(bb, histFile, false);
		}
		
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
		outstream.println("Reads Out:          "+readsOut);
		assert(!errorState) : "An error was encountered.";
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Processes a single read by scoring it with the neural network.
	 * Converts sequence to vector representation, computes score, applies filtering,
	 * updates histograms, and optionally writes to output.
	 *
	 * @param r The read to process
	 * @param bb ByteBuilder for string operations
	 * @param bsw Output writer for filtered sequences
	 * @param vec Pre-allocated vector array for sequence encoding
	 * @return true if read passes filters and is written to output
	 */
	private boolean processRead(Read r, ByteBuilder bb, ByteStreamWriter bsw, float[] vec) {
		if(r==null) {return false;}
		float result=(parseHeader ? Parse.parseFloat(r.id, "result=", '\t') : 0);
		float score=score(r.bases, vec, k, net, rcomp);
		int iscore=Tools.mid(0, Math.round(score*100), 100);
		boolean pass=(!filter ? true : (score>=cutoff)==highpass);
		readsOut+=(pass ? 1 : 0);
		if(result<0.5f) {mhist.increment(iscore);}
		else {phist.increment(iscore);}
		if(bsw!=null && pass) {
			if(annotate) {r.id+=String.format("\tscore=%.4f", score);}
			bsw.print(r.toFasta(bb.clear()).nl());
		}
		return pass;
	}
	
	/**
	 * Scores a sequence using a neural network with optional reverse complement checking.
	 * Converts sequence to k-mer vector, feeds through network, and returns maximum score
	 * from forward and reverse orientations if reverse complement scoring is enabled.
	 *
	 * @param bases DNA sequence as byte array
	 * @param vec Pre-allocated vector for sequence encoding
	 * @param k K-mer size for sequence vectorization
	 * @param net Neural network classifier
	 * @param rcomp Whether to score reverse complement and take maximum
	 * @return Maximum neural network confidence score (0.0-1.0)
	 */
	public static float score(byte[] bases, float[] vec, int k, CellNet net, boolean rcomp) {
		SequenceToVector.fillVector(bases, vec, k);
		net.applyInput(vec);
		final float r;
		final float f=net.feedForward();
		if(rcomp) {
			Vector.reverseComplementInPlaceFast(bases);
			SequenceToVector.fillVector(bases, vec, k);
			net.applyInput(vec);
			r=net.feedForward();
			Vector.reverseComplementInPlaceFast(bases);
		}else {r=f;}
		return Tools.max(r, f);
	}
	
	/**
	 * Scores a sequence using a neural network in forward orientation only.
	 * Converts sequence to k-mer vector and feeds through network.
	 *
	 * @param bases DNA sequence as byte array
	 * @param vec Pre-allocated vector for sequence encoding
	 * @param k K-mer size for sequence vectorization
	 * @param net Neural network classifier
	 * @return Neural network confidence score (0.0-1.0)
	 */
	public static float score(byte[] bases, float[] vec, int k, CellNet net) {
		assert(vec!=null);
		SequenceToVector.fillVector(bases, vec, k);
		net.applyInput(vec);
		final float f=net.feedForward();
		return f;
	}
	
	/**
	 * Scores a subsequence using a neural network.
	 * Converts specified region of sequence to k-mer vector and feeds through network.
	 *
	 * @param bases DNA sequence as byte array
	 * @param vec Pre-allocated vector for sequence encoding
	 * @param k K-mer size for sequence vectorization
	 * @param net Neural network classifier
	 * @param from Starting position in sequence (inclusive)
	 * @param to Ending position in sequence (exclusive)
	 * @return Neural network confidence score (0.0-1.0)
	 */
	public static float score(byte[] bases, float[] vec, int k, CellNet net, int from, int to) {
		SequenceToVector.fillVector(bases, vec, k, from, to);
		net.applyInput(vec);
		final float f=net.feedForward();
		return f;
	}
	
	/*--------------------------------------------------------------*/
	
	/** Input FASTA/FASTQ file path */
	private String in1=null;
	/** Output file path for scored/filtered sequences */
	private String out1=null;
	/** Path to neural network model file */
	private String netFile=null;
	/** Path to histogram output file for score distributions */
	private String histFile=null;
	
	/** File format handler for input sequences */
	private final FileFormat ffin1;
	/** File format handler for output sequences */
	private final FileFormat ffout1;
	/** File format handler for neural network model */
	private final FileFormat ffnet;

	/** Histogram of scores for positive class sequences */
	private LongList phist=new LongList(101);
	/** Histogram of scores for negative class sequences */
	private LongList mhist=new LongList(101);
	
	/*--------------------------------------------------------------*/

	/** Number of reads written to output after filtering */
	private long readsOut=0;
	/** Maximum number of reads to process (-1 for unlimited) */
	private long maxReads=-1;
	/** Tracks whether any errors occurred during processing */
	private boolean errorState=false;

	/** K-mer size for sequence vectorization */
	private int k=0;
	/** Whether to score reverse complement and take maximum score */
	private boolean rcomp=false;
	/** Whether to parse expected result values from sequence headers */
	private boolean parseHeader=false;
	/** Width parameter for vector encoding (-1 for auto-detection from network) */
	private int width=-1;
	/** Neural network model for sequence scoring */
	private final CellNet net;

	/** Whether to filter sequences based on score cutoff */
	private boolean filter=false;
	/**
	 * Whether filtering retains high-scoring (true) or low-scoring (false) sequences
	 */
	private boolean highpass=true;
	/** Whether to add score annotations to output sequence headers */
	private boolean annotate=true;
	/** Score threshold for filtering sequences */
	private float cutoff=0.5f;
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and results */
	private java.io.PrintStream outstream=System.err;
	/** Whether to print verbose progress messages */
	public static boolean verbose=false;
	
}
