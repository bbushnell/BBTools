package scalar;

import java.io.File;
import java.util.ArrayList;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
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
 * Computes GC content and dimer metrics, either per input file or in rolling
 * windows of valid dinucleotides. Global mode prints one metric row; windowed
 * mode prints histogram means and population standard deviations. Raw mode
 * reports reverse-complement-collapsed dimer fractions instead of those summaries
 * (the final retained window when windowed); no valid dimers produces NaN fractions.
 * An instance is single-use and unsynchronized; shared parser settings are global.
 *
 * @author Brian Bushnell
 * @date Oct 2, 2025
 */
public class Scalars{

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
		try{x.process(t);}
		finally{Shared.closeStream(x.outstream);}
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
			String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && (b.isEmpty() || b.equalsIgnoreCase("null"))){b=null;}

			if(a.equals("header") || a.equals("colheader") || a.equals("columnheader")){
				header=Parse.parseBoolean(b);
			}else if(a.equals("rowheader")){
				rowheader=Parse.parseBoolean(b);
			}else if(a.equals("raw")){
				raw=Parse.parseBoolean(b);
			}else if(a.equals("break")){
				breakOnContig=Parse.parseBoolean(b);
			}else if(a.equals("window")){
				window=Integer.parseInt(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(new File(arg).exists()){
				in.add(arg);
			}else{
				throw new IllegalArgumentException("Unknown parameter "+args[i]);
			}
		}

		{//Process parser fields
			if(parser.in2!=null){throw new IllegalArgumentException("Scalars does not support in2; use interleaved input or separate input files.");}
			Parser.processQuality();

			maxReads=parser.maxReads;
			if(parser.in1!=null){in.add(parser.in1);}
			out=parser.out1;
		}


		in=Tools.getFileOrFiles(in, true, false, false, false);
		ScalarData.checkFiles(in, out, parser.overwrite, parser.append);
		ffout=FileFormat.testOutput(out, FileFormat.TXT, null, true, parser.overwrite, parser.append, false);
		dimers=new KmerTracker(2, window);
	}

	/**
	 * Processes input reads and calculates compositional metrics.
	 * Either accumulates global dimer counts or builds histograms from sliding windows.
	 * @param t Timer for performance tracking
	 */
	void process(Timer t){

		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout);
		try{
			for(int i=0; i<in.size(); i++){
				String fname=in.get(i);
				FileFormat ffin=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
				processOne(ffin);
				outputResults(bsw, header && i==0);
				Tools.fill(hist, 0);
			}
		}finally{
			if(bsw!=null){errorState|=bsw.poisonAndWait();}
		}

		if(verbose){outstream.println("Finished reading data.");}

		t.stop();
		if(printTime){
			outstream.println("Time:                         \t"+t);
			outstream.println("Reads Processed:    "+readsProcessed+" \t"+
				Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
		}
		if(errorState){throw new RuntimeException("An error occurred reading or writing scalar data");}
	}

	/**
	 * Processes input reads and calculates compositional metrics.
	 * Either accumulates global dimer counts or builds histograms from sliding windows.
	 * Windowed state spans read boundaries unless break=t. Both mates are included.
	 * @param ffin Input descriptor; this method owns and closes its reader
	 */
	void processOne(FileFormat ffin){
		final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin, null);
		cris.start();

		try{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			dimers.clearAll();
			while(ln!=null && reads!=null && reads.size()>0){
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					readsProcessed+=r1.pairCount();
					basesProcessed+=r1.pairLength();
					if(window<1){
						dimers.add(r1.bases);
						if(r1.mate!=null){dimers.add(r1.mate.bases);}
					}else{
						addWindowed(r1.bases);
						if(r1.mate!=null){addWindowed(r1.mate.bases);}
					}
				}
				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());}
		}finally{
			errorState|=ReadWrite.closeStreams(cris);
		}
	}

	/**
	 * Processes a sequence with sliding window analysis.
	 * For each full valid-dimer window, calculates metrics and updates histograms.
	 * Ambiguity resets dimer construction but does not clear prior window entries.
	 * Undefined metric values retain the historical NaN-to-bin-zero behavior.
	 * @param bases Sequence to process
	 */
	void addWindowed(byte[] bases){
		if(breakOnContig){dimers.clearAll();}
		if(bases==null){return;}
		for(byte b : bases){
			boolean newValid=dimers.addWindowed(b);
			if(newValid){
				//NOTE [scalar/Scalars#001] RESOLVED — NOT a crash (verified against tracker/KmerTracker): hist is [14][1025] and
				//metrics bin as (int)(metric()*1024). I traced all 14 KmerTracker metrics — each is mathematically ∈[0,1] (numerator
				//⊆ or ±denominator), so the index is ∈[0,1024] and fits. The degenerate case (zero denominator) returns NaN, and
				//(int)(NaN*1024)==0 in Java → in-bounds bin 0. So no AIOOBE; the only artifact is that a degenerate window's NaN
				//metric lands in bin 0, slightly skewing that histogram. Downgraded from the earlier AIOOBE concern. LOW/cosmetic.
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
		if(bsw==null){return;}
		ByteBuilder bb=new ByteBuilder();
		final int decimals=ScalarData.decimals;
		if(rowheader && header){bb.append("Header\t");}
		if(raw){
			if(header){
//				bb.append("AA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n");
				bb.append("#AA\tAC\tAG\tAT\tCA\tCC\tCG\tGA\tGC\tTA\n");
			}
			if(rowheader){bb.append("kmers\t");}
			//NOTE [scalar/Scalars#002] LOW/dev: 1f/Tools.sum(dimers.counts) → Infinity when the input had no valid dimers (empty
			//file / all-N), so the raw counts print as NaN. Cosmetic (pattern a). Also the bb.set(bb.length()-1,'\n') calls (raw
			//and windowed) assume >=1 column was appended; safe in practice (2-mer always yields entries) but SIOOB if bb were empty.
			float mult=1f/Tools.sum(dimers.counts);
			for(int i=0; i<dimers.counts.length; i++){
				int r=AminoAcid.reverseComplementBinary(i, 2);
				if(i<r){
					bb.appendt((dimers.counts[i]+dimers.counts[r])*mult, decimals);
				}else if(i==r){
					bb.appendt((dimers.counts[i])*mult, decimals);
				}
			}
			bb.set(bb.length()-1, '\n');
		}else{
			if(header){
				bb.append("#GC\tSTR\tHH\tPP\tAAAT\tCCCG\tHMH\tHHPP\tACTG\tACAG\tCAGA\tCCMCG\tATMTA\tAT\n");
			}
			if(window<1){
				if(rowheader){bb.append("Mean\t");}
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
			}else{
				if(rowheader){bb.append("Mean\t");}
				for(int i=0; i<hist.length; i++){
					bb.appendt(Tools.averageHistogram(hist[i])/1024, decimals);
				}
				bb.set(bb.length()-1, '\n');
				if(rowheader){bb.append("STDev\t");}
				for(int i=0; i<hist.length; i++){
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
	private boolean header=true;
	/** Whether to print row headers */
	private boolean rowheader=false;
	/** Whether to print timing information */
	private boolean printTime=false;
	/** Emit ten reverse-complement-collapsed dimer fractions instead of metric summaries. */
	private boolean raw=false;
	/** Histograms for fourteen metrics; each contains bins 0 through 1024 inclusive. */
	private long[][] hist=new long[14][1025];
	/** Clear rolling-window state before each sequence when true; global mode is unaffected. */
	private boolean breakOnContig=false;

	/*--------------------------------------------------------------*/

	/** Maximum number of reads to process (-1 for unlimited) */
	private long maxReads=-1;
	/** Whether an error occurred during processing */
	private boolean errorState=false;


	/** Individual sequences processed, counting both mates. */
	private long readsProcessed=0;
	/** Bases processed, counting both mates and ambiguous bases. */
	private long basesProcessed=0;

	/*--------------------------------------------------------------*/

	/** Output stream for messages */
	private java.io.PrintStream outstream=System.err;
	/** Whether to print verbose progress messages */
	public static boolean verbose=false;

}
