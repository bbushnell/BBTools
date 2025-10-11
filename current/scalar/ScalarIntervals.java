package scalar;

import java.io.File;
import java.util.ArrayList;

import clade.Clade;
import clade.Comparison;
import clade.SendClade;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sketch.DisplayParams;
import sketch.SendSketch;
import sketch.Sketch;
import sketch.SketchMakerMini;
import sketch.SketchObject;
import sketch.SketchTool;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;
import tax.TaxTree;
import tracker.KmerTracker;

/**
 * Calculates compositional scalar metrics from sequencing data.
 * Computes GC-independent metrics (HH, CAGA, strandedness, etc.) either globally
 * or using a sliding window to characterize within-genome variance.
 * Outputs mean and standard deviation for each metric.
 *
 * @author Brian Bushnell
 * @date Oct 2, 2025
 */
public class ScalarIntervals {

	/**
	 * Main entry point for the Scalars program.
	 * @param args Command-line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();

		//Create an instance of this class
		ScalarIntervals x=new ScalarIntervals(args);

		//Run the object
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructs a Scalars instance and parses command-line arguments.
	 * Supports windowed or global analysis of compositional metrics.
	 * @param args Command-line arguments
	 */
	public ScalarIntervals(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null/*getClass()*/, false);
			args=pp.args;
			outstream=pp.outstream;
		}

		Parser parser=new Parser();
		parser.out1="stdout.txt";
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("header") || a.equals("colheader") || a.equals("columnheader")){
				header=Parse.parseBoolean(b);
			}else if(a.equals("raw")){
				raw=Parse.parseBoolean(b);
			}else if(a.equals("window")){
				window=Parse.parseIntKMG(b);
			}else if(a.equals("interval")){
				interval=Parse.parseIntKMG(b);
			}else if(a.equals("shred")){
				interval=window=Parse.parseIntKMG(b);
			}else if(a.equals("break")){
				breakOnContig=Parse.parseBoolean(b);
			}else if(a.equals("sketch") | a.equals("bbsketch")){
				ScalarData.makeSketch=Parse.parseBoolean(b);
			}else if(a.equals("clade") || a.equals("quickclade")){
				ScalarData.makeClade=Parse.parseBoolean(b);
			}else if(a.equals("printname") || a.equals("printnames")){
				printName=Parse.parseBoolean(b);
			}else if(a.equals("parsetaxid") || a.equals("parsetax") || a.equals("parsetid")){
				ScalarData.parseTID=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
//				ReadWrite.verbose=verbose;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(new File(arg).exists()) {
				in.add(arg);
			}else{
				//				throw new RuntimeException("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
		}

		{//Process parser fields
			Parser.processQuality();

			maxReads=parser.maxReads;
			if(parser.in1!=null) {in.add(parser.in1);}
			out=parser.out1;
		}
		
		
		in=Tools.getFileOrFiles(in, true, false, false, false);
		ffout=FileFormat.testOutput(out, FileFormat.TXT, null, true, true, false, false);
	}

	/**
	 * Processes input reads and calculates compositional metrics.
	 * Either accumulates global dimer counts or builds histograms from sliding windows.
	 * @param t Timer for performance tracking
	 */
	void process(Timer t){

		readsProcessed=0;
		basesProcessed=0;
		if(verbose) {outstream.println("callingToIntervals");}
		
		if(verbose){outstream.println("Finished reading data; printing to "+out);}
		
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout);
		for(int i=0; i<in.size(); i++) {
			FileFormat ffin=FileFormat.testInput(in.get(i), FileFormat.FASTA, null, true, true);
			ScalarData data=toIntervals(ffin, window, interval, minlen, breakOnContig, maxReads);
			data.print(bsw, printName, header && i==0);
		}
		if(bsw!=null) {bsw.poison();}

		t.stop();
		if(printTime) {
			outstream.println();
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		}
		assert(!errorState) : "An error was encountered.";
	}
	
	public static ScalarData toIntervals(String fname, int window, int interval, int minlen, boolean breakOnContig, long maxReads) {
		if(verbose) {System.err.println("callingToIntervals(string)");}
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		return toIntervals(ff, window, interval, minlen, breakOnContig, maxReads);
	}
	
	public static ScalarData toIntervals(FileFormat ff, int window, int interval, int minlen, boolean breakOnContig, long maxReads) {
		if(verbose) {System.err.println("callingToIntervals(ff)");}
		int tid=(ScalarData.parseTID ? bin.BinObject.parseTaxID(ff.name()) : -1);
		final ConcurrentReadInputStream cris;
		cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ff, null);
		cris.start();
		ScalarData data=toIntervals(cris, window, interval, minlen, tid, breakOnContig, maxReads);
		boolean errorState=ReadWrite.closeStreams(cris);
		if(verbose){System.err.println("Finished reading data.");}
		if(errorState){System.err.println("Something went wrong reading "+ff.name());}
		return data;
	}
	
	public static ScalarData toIntervals(ConcurrentReadInputStream cris, int window, int interval, int minlen, 
		int tid, boolean breakOnContig, long maxReads) {
		if(ScalarData.makeSketch || ScalarData.makeClade) {
			return toIntervals_multi(cris, window, interval, minlen, tid, breakOnContig, maxReads);
		}
		if(verbose) {System.err.println("callingToIntervals(cris)");}
		ScalarData data=new ScalarData(printName, -1);
		final KmerTracker dimers=new KmerTracker(2, window);
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			while(ln!=null && reads!=null && reads.size()>0){
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx), r2=r1.mate;
					readsProcessed+=r1.pairCount();
					basesProcessed+=r1.pairLength();

					data.add(r1, dimers, null, interval, minlen, tid, breakOnContig);
					if(r2!=null) {data.add(r2, dimers, null, interval, minlen, tid, breakOnContig);}
				}

				cris.returnList(ln);
				if(verbose){System.err.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		if(verbose) {System.err.println("finished ToIntervals(cris)");}
		return data;
	}
	
	public static ScalarData toIntervals_multi(ConcurrentReadInputStream cris, 
		int window, int interval, int minlen, 
		int tid, boolean breakOnContig, long maxReads) {
		if(verbose) {System.err.println("callingToIntervals(cris)");}
		
		if(ScalarData.makeClade) {Clade.DELETE_COUNTS=false;}
		SketchMakerMini smm=null;
		if(ScalarData.makeSketch) {
			SketchObject.processSSU=false;
			SketchObject.postParse();
			DisplayParams params=new DisplayParams();
			params.format=DisplayParams.FORMAT_JSON;
			params.taxLevel=TaxTree.GENUS;
			params.maxRecords=2;
			SketchTool tool=new SketchTool(2000, params);
			smm=new SketchMakerMini(tool, SketchObject.ONE_SKETCH, params);
		}
		final KmerTracker dimers=new KmerTracker(2, window);
		ArrayList<ScalarData> sdlist=new ArrayList<ScalarData>();
//		System.err.println("Made KmerTracker "+2+", "+window);
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			while(ln!=null && reads!=null && reads.size()>0){
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					readsProcessed+=r1.pairCount();
					basesProcessed+=r1.pairLength();
					
					if(r1.length()<minlen) {continue;}
					ScalarData sd=new ScalarData(printName, -1);
					sd.add(r1, dimers, smm, interval, minlen, tid, breakOnContig);
					if(ScalarData.makeSketch) {
						if(r1.pairLength()>=ScalarData.minSketchSize) {
							sd.sketch=smm.toSketch(1);
						}else {smm.heap().clear(false);}
					}
					if(sd.clade!=null) {
						assert(!Clade.DELETE_COUNTS);
						sd.clade.finish();
					}
					sdlist.add(sd);
				}

				cris.returnList(ln);
				if(verbose){System.err.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		if(verbose) {System.err.println("finished ToIntervals(cris)");}
		
		ArrayList<Sketch> sketches=new ArrayList<Sketch>();
		ArrayList<Clade> clades=new ArrayList<Clade>();
		for(ScalarData sd : sdlist) {
			if(sd.sketch!=null) {sketches.add(sd.sketch);}
			if(sd.clade!=null) {clades.add(sd.clade);}
		}
		if(!clades.isEmpty()) {
			boolean success=SendClade.sendAndLabel(clades);
			for(ScalarData sd : sdlist) {
				if(sd.clade!=null && sd.clade.taxID>0) {
					sd.taxIDs.fill(sd.clade.taxID);
				}
			}
		}
		if(!sketches.isEmpty()) {
			boolean success=SendSketch.sendAndLabel(sketches, null);
			for(ScalarData sd : sdlist) {
				if(sd.sketch!=null && sd.sketch.taxID>0) {
					sd.taxIDs.fill(sd.sketch.taxID);
				}
			}
		}
		
		ScalarData data=sdlist.isEmpty() ? new ScalarData(false, -1) : sdlist.get(0);
		for(int i=1; i<sdlist.size(); i++) {
			ScalarData sd=sdlist.set(i, null);
			data.add(sd);
		}
		
		return data;
	}
	
//	public static void toIntervals(Read r, KmerTracker dimers, ScalarData sd, 
//			int interval, int minlen, int tid, boolean breakOnContig) {
//		if(verbose) {System.err.println("calling ToIntervals(read)");}
//		if(r==null || (breakOnContig && r.length()<minlen)) {return;}
//		if(breakOnContig) {dimers.clearAll();}
//		final byte[] bases=r.bases;
//		if(parseTID && tid<0) {tid=bin.BinObject.parseTaxID(r.name());}
//		if(dimers.window>0) {
//			for(byte b : bases) {
//				boolean newValid=dimers.addWindowed(b);
//				if(newValid && interval>0 && dimers.count()>=interval) {
//					toInterval(dimers, sd);
//					if(sd.taxID!=null) {sd.taxID.add(tid);}
//					if(sd.name!=null) {sd.name.add(r.name());}
//				}
//			}
//		}else {dimers.add(bases);}
//		if(verbose) {System.err.println("dimers.count()="+dimers.count()+", minlen="+minlen);}
//		if(dimers.count()>=minlen) {
//			toInterval(dimers, sd);
//			if(sd.taxID!=null) {sd.taxID.add(tid);}
//			if(sd.name!=null) {sd.name.add(r.name());}
//		}
//	}
//	
//	private static void toInterval(KmerTracker dimers, ScalarData data) {
//		if(verbose) {System.err.println("calling toInterval(dimers)");}
//		data.gc.add(dimers.GC());
//		data.hh.add(dimers.HH());
//		data.caga.add(dimers.CAGA());
//		dimers.resetCount();
//	}

	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/

	/** Input file path */
	private ArrayList<String> in=new ArrayList<String>();
	/** Output file path */
	private String out=null;

//	/** Input file format */
//	private final FileFormat ffin;
	/** Output file format */
	private final FileFormat ffout;
	/** Whether to print column headers */
	private boolean header=false;
	/** Whether to print row headers */
	private boolean rowheader=false;
	/** Whether to print timing information */
	private boolean printTime=false;
	private boolean raw=false;


	/** Window size for sliding window analysis (0 for global analysis) */
	private int window=0;
	private int interval=5000;
	private int minlen=500;
	private boolean breakOnContig=true;
	private static boolean printName=false;

	static long readsProcessed=0, basesProcessed=0;
	
	/*--------------------------------------------------------------*/

	/** Maximum number of reads to process (-1 for unlimited) */
	private long maxReads=-1;
	/** Whether an error occurred during processing */
	private boolean errorState=false;

	/*--------------------------------------------------------------*/

	/** Output stream for messages */
	private java.io.PrintStream outstream=System.err;
	/** Whether to print verbose progress messages */
	public static boolean verbose=false;

}
