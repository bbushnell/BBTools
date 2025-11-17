package driver;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Feb 21, 2017
 *
 */
public class PlotGC {
	
	/**
	 * Program entry point for GC content analysis.
	 * Creates PlotGC instance and executes the processing pipeline.
	 * @param args Command-line arguments for input/output files and parameters
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		PlotGC x=new PlotGC(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs PlotGC instance and parses command-line arguments.
	 * Initializes file formats, output streams, and processing parameters.
	 * Sets up input validation and error checking for file operations.
	 * @param args Command-line arguments containing file paths and options
	 */
	public PlotGC(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("interval")){
				interval=Integer.parseInt(b);
			}else if(a.equals("offset")){
				offset=Integer.parseInt(b);
			}else if(a.equals("printshortbins") || a.equals("psb")){
				printShortBins=Parse.parseBoolean(b);
			}else if(a.equals("out")){
				out1=b;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;
			
			extin=parser.extin;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TEXT, ".txt", true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, extin, true, true);
	}
	
	/**
	 * Main processing method that calculates GC content across sequences.
	 * Reads input sequences, applies sliding window analysis, and outputs results.
	 * Processes each sequence in intervals, counting ACGT bases and calculating GC percentage.
	 * @param t Timer for tracking execution performance
	 */
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null, null, null);
			cris.start();
			if(verbose){outstream.println("Started cris");}
		}

		final TextStreamWriter tsw;
		if(out1!=null){
			tsw=new TextStreamWriter(ffout1);
			tsw.start();
			tsw.println("name\tinterval\tstart\tstop\trunningStart\trunningStop\tgc");
		}else{tsw=null;}
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			long rStart=0, rStop=0;
			
			int[] acgt=KillSwitch.allocInt1D(4);
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					final int initialLength1=r1.length();
					
					readsProcessed++;
					basesProcessed+=initialLength1;
				}
				
				for(Read r : reads){
					if(r.bases!=null) {
						Arrays.fill(acgt, 0);
						byte[] bases=r.bases;
						int next=interval-1;
						int start=0, i=0;
						for(; i<bases.length; i++){
							int num=AminoAcid.baseToNumber[bases[i]];
							if(num>=0){acgt[num]++;}
							if(i>=next){
								int len=(i-start);
								rStop=rStart+len;
								String s=toGC(r.id, start, i, rStart, rStop, acgt);
								start=i+1;
								rStart=rStop+1;
								next=i+interval;
								Arrays.fill(acgt, 0);
								if(tsw!=null && s!=null){
									tsw.print(s);
								}
							}
						}
						if(printShortBins && i>start){
							int len=(i-start)-1;
							rStop=rStart+len;
							String s=toGC(r.id, start, i-1, rStart, rStop, acgt);
							start=i+1;
							rStart=rStop+1;
							next=i+interval;
							Arrays.fill(acgt, 0);
							if(tsw!=null && s!=null){
								tsw.print(s);
							}
						}
					}
				}
				
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadWrite.closeStreams(cris);
		
		if(tsw!=null){
			errorState|=tsw.poisonAndWait();
		}
		
		t.stop();
		
		if(showSpeed){
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		}
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * Converts base composition counts to formatted GC content string.
	 * Calculates GC percentage from ACGT counts and formats output line.
	 *
	 * @param name Sequence identifier
	 * @param start Start position within sequence (0-based)
	 * @param stop Stop position within sequence (0-based, inclusive)
	 * @param rstart Running start position across all sequences
	 * @param rstop Running stop position across all sequences
	 * @param acgt Array of base counts [A, C, G, T]
	 * @return Formatted tab-separated string with positions and GC percentage
	 */
	private String toGC(String name, int start, int stop, long rstart, long rstop, int[] acgt){
		int at=acgt[0]+acgt[3];
		int gc=acgt[1]+acgt[2];
		float sum=Tools.max(1, at+gc);
		float gcf=gc/sum;
		return Tools.format("%s\t%d\t%d\t%d\t%d\t%d\t%.3f\n", name, stop-start+1, start+offset, stop+offset, rstart+offset, rstop+offset, gcf);
	}
	
	/*--------------------------------------------------------------*/
	
	
	/*--------------------------------------------------------------*/
	
	/** Primary input file path for sequence data */
	private String in1=null;
	
	/** Output file path for GC content results */
	private String out1="stdout.txt";
	
	/** Input file extension override */
	private String extin=null;
	
	/*--------------------------------------------------------------*/

	/** Maximum number of reads to process (-1 for unlimited) */
	private long maxReads=-1;
	
	/** Window size in bases for GC content calculation */
	private int interval=1000;
	/** Position offset to add to output coordinates */
	private int offset=0;

	/** Whether to display processing speed statistics */
	private boolean showSpeed=false;
	/** Whether to output final partial bins shorter than full interval */
	private boolean printShortBins=true;
	
	/*--------------------------------------------------------------*/
	
	/** File format specification for primary input file */
	private final FileFormat ffin1;

	/** File format specification for output file */
	private final FileFormat ffout1;
	
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and errors */
	private PrintStream outstream=System.err;
	/** Global verbosity flag for detailed logging */
	public static boolean verbose=false;
	/** Flag indicating whether processing encountered errors */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files */
	private boolean append=false;
	
}
