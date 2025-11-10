package scalar;

import java.io.File;
import java.util.ArrayList;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
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
public class Scalars {

	/**
	 * Main entry point for the Scalars program.
	 * @param args Command-line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();

		//Create an instance of this class
		Scalars x=new Scalars(args);

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
	public Scalars(String[] args){

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
			}else if(a.equals("rowheader")){
				rowheader=Parse.parseBoolean(b);
			}else if(a.equals("raw")){
				raw=Parse.parseBoolean(b);
			}else if(a.equals("window")){
				window=Integer.parseInt(b);
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
		dimers=new KmerTracker(2, window);
	}

	/**
	 * Processes input reads and calculates compositional metrics.
	 * Either accumulates global dimer counts or builds histograms from sliding windows.
	 * @param t Timer for performance tracking
	 */
	void process(Timer t){

		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout);
		for(int i=0; i<in.size(); i++) {
			String fname=in.get(i);
			FileFormat ffin=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
			processOne(ffin);
			outputResults(bsw, header && i==0);
			Tools.fill(hist, 0);
		}
		errorState=bsw.poisonAndWait() | errorState;
		
		if(verbose){outstream.println("Finished reading data.");}

		t.stop();
		if(printTime) {
			outstream.println("Time:                         \t"+t);
			outstream.println("Reads Processed:    "+readsProcessed+" \t"+
				Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
		}
		assert(!errorState) : "An error was encountered.";
	}

	/**
	 * Processes input reads and calculates compositional metrics.
	 * Either accumulates global dimer counts or builds histograms from sliding windows.
	 * @param t Timer for performance tracking
	 */
	void processOne(FileFormat ffin){
		final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin, null);
		cris.start();

		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		dimers.clearAll();
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}

			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);
				readsProcessed+=r1.pairCount();
				basesProcessed+=r1.pairLength();

				if(window<1) {
					dimers.add(r1.bases);
					if(r1.mate!=null) {dimers.add(r1.mate.bases);}
				}else {
					addWindowed(r1.bases);
					if(r1.mate!=null) {addWindowed(r1.mate.bases);}
				}
			}

			cris.returnList(ln);
			if(verbose){outstream.println("Returned a list.");}
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
		errorState=ReadWrite.closeStreams(cris) | errorState;
	}

	/**
	 * Processes a sequence with sliding window analysis.
	 * For each valid window, calculates all compositional metrics and updates histograms.
	 * @param bases Sequence to process
	 */
	void addWindowed(byte[] bases) {
		for(byte b : bases) {
			boolean newValid=dimers.addWindowed(b);
			if(newValid) {
				hist[0][(int)(dimers.GC()*1024)]++;
				hist[1][(int)(dimers.strandedness()*1024)]++;
				hist[2][(int)(dimers.HH()*1024)]++;
				hist[3][(int)(dimers.PP()*1024)]++;
				hist[4][(int)(dimers.AAAT()*1024)]++;
				hist[5][(int)(dimers.CCCG()*1024)]++;
				hist[6][(int)(dimers.HMH()*1024)]++;
				hist[7][(int)(dimers.HHPP()*1024)]++;
				hist[8][(int)(dimers.ACTG()*1024)]++;
				hist[9][(int)(dimers.ACAG()*1024)]++;
				hist[10][(int)(dimers.CAGA()*1024)]++;
				hist[11][(int)(dimers.CCMCG()*1024)]++;
				hist[12][(int)(dimers.ATMTA()*1024)]++;
				hist[13][(int)(dimers.AT()*1024)]++;
			}
		}
	}

	/**
	 * Outputs calculated metrics to file.
	 * For global mode: outputs mean values for each metric.
	 * For windowed mode: outputs mean and standard deviation across all windows.
	 */
	private void outputResults(ByteStreamWriter bsw, boolean header){
		ByteBuilder bb=new ByteBuilder();
		final int decimals=ScalarData.decimals;
		if(rowheader && header) {bb.append("Header\t");}
		if(raw) {
			if(header){
//				bb.append("AA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n");
				bb.append("AA\tAC\tAG\tAT\tCA\tCC\tCG\tGA\tGC\tTA\n");
			}
			if(rowheader) {bb.append("kmers\t");}
			float mult=1f/Tools.sum(dimers.counts);
			for(int i=0; i<dimers.counts.length; i++) {
				int r=AminoAcid.reverseComplementBinary(i, 2);
				if(i<r) {
					bb.appendt((dimers.counts[i]+dimers.counts[r])*mult, decimals);
				}else if(i==r) {
					bb.appendt((dimers.counts[i])*mult, decimals);
				}
			}
			bb.set(bb.length()-1, '\n');
		}else{
			if(header){
				bb.append("GC\tSTR\tHH\tPP\tAAAT\tCCCG\tHMH\tHHPP\tACTG\tACAG\tCAGA\tCCMCG\tATMTA\tAT\n");
			}
			if(window<1) {
				if(rowheader) {bb.append("Mean\t");}
				bb.appendt(dimers.GC(), decimals);
				bb.appendt(dimers.strandedness(), decimals);
				bb.appendt(dimers.HH(), decimals);
				bb.appendt(dimers.PP(), decimals);
				bb.appendt(dimers.AAAT(), decimals);
				bb.appendt(dimers.CCCG(), decimals);
				bb.appendt(dimers.HMH(), decimals);
				bb.appendt(dimers.HHPP(), decimals);
				bb.appendt(dimers.ACTG(), decimals);
				bb.appendt(dimers.ACAG(), decimals);
				bb.appendt(dimers.CAGA(), decimals);
				bb.appendt(dimers.CCMCG(), decimals);
				bb.appendt(dimers.ATMTA(), decimals);
				bb.append(dimers.AT(), decimals);
				bb.nl();
			}else {
				if(rowheader) {bb.append("Mean\t");}
				for(int i=0; i<hist.length; i++) {
					bb.appendt(Tools.averageHistogram(hist[i])/1024, decimals);
				}
				bb.set(bb.length()-1, '\n');
				if(rowheader) {bb.append("STDev\t");}
				for(int i=0; i<hist.length; i++) {
					bb.appendt(Tools.standardDeviationHistogram(hist[i])/1024, decimals);
				}
				bb.set(bb.length()-1, '\n');
			}
		}
		bsw.print(bb);
	}

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

	/** Window size for sliding window analysis (0 for global analysis) */
	private int window=0;
	/** K-mer tracker for dimer counting */
	private final KmerTracker dimers;
	/** Whether to print column headers */
	private boolean header=false;
	/** Whether to print row headers */
	private boolean rowheader=false;
	/** Whether to print timing information */
	private boolean printTime=false;
	private boolean raw=false;
	/** Histograms for each metric in windowed mode (8 metrics, 1025 bins each) */
	private long[][] hist=new long[14][1025];
	private boolean breakOnContig=false;
	
	/*--------------------------------------------------------------*/

	/** Maximum number of reads to process (-1 for unlimited) */
	private long maxReads=-1;
	/** Whether an error occurred during processing */
	private boolean errorState=false;
	

	private long readsProcessed=0;
	private long basesProcessed=0;

	/*--------------------------------------------------------------*/

	/** Output stream for messages */
	private java.io.PrintStream outstream=System.err;
	/** Whether to print verbose progress messages */
	public static boolean verbose=false;

}
