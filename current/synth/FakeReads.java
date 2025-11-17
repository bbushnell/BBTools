package synth;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.Vector;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;
/**
 * @author Brian Bushnell
 * @date Sep 11, 2012
 *
 */
public class FakeReads {

	/** Program entry point.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		FakeReads x=new FakeReads(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs FakeReads instance and parses command-line arguments.
	 * Sets up input/output file formats and processing parameters.
	 * @param args Command-line arguments specifying input files, output files, and parameters
	 */
	public FakeReads(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(a.equals("passes")){
				assert(false) : "'passes' is disabled.";
//				passes=Integer.parseInt(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("addspacer") || a.equals("addspace") || a.equals("usespacer")){
				addSpacer=Parse.parseBoolean(b);
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("t") || a.equals("threads")){
				Shared.setThreads(b);
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				in1=b;
			}else if(a.equals("out") || a.equals("output") || a.equals("out1") || a.equals("output1")){
				out1=b;
			}else if(a.equals("out2") || a.equals("output2")){
				out2=b;
			}else if(a.equals("identifier") || a.equals("id")){
				identifier=b;
			}else if(a.equals("qfin") || a.equals("qfin1")){
				qfin1=b;
			}else if(a.equals("qfout") || a.equals("qfout1")){
				qfout1=b;
			}else if(a.equals("qfout2")){
				qfout2=b;
			}else if(a.equals("extin")){
				extin=b;
			}else if(a.equals("extout")){
				extout=b;
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.startsWith("minscaf") || a.startsWith("mincontig")){
				stream.FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
			}else if(a.equals("ml") || a.equals("minlen") || a.equals("minlength")){
				minReadLength=Integer.parseInt(b);
			}else if(a.equals("length") || a.equals("maxlen") || a.equals("length")){
				desiredLength=Integer.parseInt(b);
			}else if(a.equals("split")){
				SPLITMODE=Parse.parseBoolean(b);
			}else if(a.equals("overlap")){
				SPLITMODE=true;
				overlap=Integer.parseInt(b);
			}else if(in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				in1=arg;
				if(arg.indexOf('#')>-1 && !new File(arg).exists() && b!=null){//b!=null is implied
					in1=b.replace("#", "1");
				}
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
		}
		
		if(identifier==null){identifier="";}
		else{identifier=identifier+"_";}
		
		if(!addSpacer){spacer="";}
		
//		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
//			in2=in1.replace("#", "2");
//			in1=in1.replace("#", "1");
//		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
//		if(in2!=null){
//			if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
//			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
//		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
//			if(ReadWrite.isCompressed(in1)){ByteFile.FORCE_MODE_BF2=true;}
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
	}
	
	/**
	 * Main processing method that converts single-end reads to paired-end reads.
	 * Extracts sequences from both ends of input reads, reverse-complements the
	 * second read, and outputs as paired reads with proper mate relationships.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null, qfin1, null);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
		if(verbose){System.err.println("Input is "+(paired ? "paired" : "unpaired"));}

		ConcurrentReadOutputStream ros=null;
		if(out1!=null){
			final int buff=4;
			
			if(cris.paired() && out2==null && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(out1))) : "out1 and out2 have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, false);
			ros.start();
		}
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				ArrayList<Read> fake=new ArrayList<Read>(reads.size());

				for(int idx=0; idx<reads.size(); idx++){
					Read r=reads.get(idx);
					{
						readsProcessed++;
						basesProcessed+=r.length();
					}
					assert(r.mate==null);
					
					boolean remove=r.length()<minReadLength || (minReadLength+overlap)<2;
					
					if(remove){
						//Do nothing
					}else{
						int len=Tools.min(r.length(), desiredLength);
						if(SPLITMODE){len=Tools.min(r.length(), (r.length()+overlap+1)/2);}
						
						byte[] bases1=KillSwitch.copyOfRange(r.bases, 0, len);
						byte[] bases2=KillSwitch.copyOfRange(r.bases, r.length()-len, r.length());
						Vector.reverseComplementInPlaceFast(bases2);
						
						byte[] qual1=null;
						byte[] qual2=null;
						if(r.quality!=null){
							qual1=KillSwitch.copyOfRange(r.quality, 0, len);
							qual2=KillSwitch.copyOfRange(r.quality, r.quality.length-len, r.quality.length);
							Vector.reverseInPlace(qual2);
						}
						
//						public Read(byte[] s_, int chrom_, int start_, int stop_, String id_, byte[] quality_, long numericID_, int flags_){
						Read a=new Read(bases1, qual1, identifier+r.numericID+spacer+"/1", r.numericID, 0);
						Read b=new Read(bases2, qual2, identifier+r.numericID+spacer+"/2", r.numericID, Read.PAIRNUMMASK);
						a.mate=b;
						b.mate=a;
						fake.add(a);
					}
				}
				
				if(ros!=null){ros.add(fake, ln.id);}

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		if(errorState){
			throw new RuntimeException("FakeReads terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/** Tracks whether processing encountered any errors */
	public boolean errorState=false;
	
	/** Prefix added to read names in output */
	public String identifier=null;
	
	/** Primary input file path */
	private String in1=null;
	
	/** Whether to add spacer character between identifier and read number */
	private boolean addSpacer=true;
	/** Spacer character added between identifier and read number in output names */
	private String spacer=" ";
	
	/** Quality file input path for input1 */
	private String qfin1=null;

	/** Primary output file path for first reads in pairs */
	private String out1=null;
	/** Secondary output file path for second reads in pairs */
	private String out2=null;

	/** Quality file output path for first reads in pairs */
	private String qfout1=null;
	/** Quality file output path for second reads in pairs */
	private String qfout2=null;
	
	/** Input file extension override for format detection */
	private String extin=null;
	/** Output file extension override for format specification */
	private String extout=null;
	
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files instead of overwriting */
	private boolean append=false;
	
	/** Maximum number of reads to process (-1 for unlimited) */
	private long maxReads=-1;
	/** Minimum input read length required for processing */
	private int minReadLength=1;
	/** Target length for output paired reads */
	private int desiredLength=250;
	/** Overlap length between paired reads in split mode */
	private int overlap=50;
	/** Whether to use split mode for generating overlapping paired reads */
	private boolean SPLITMODE=false;
	
	/** File format specification for primary input */
	private final FileFormat ffin1;
	
	/** File format specification for primary output */
	private final FileFormat ffout1;
	/** File format specification for secondary output */
	private final FileFormat ffout2;
	
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and statistics */
	private PrintStream outstream=System.err;
	/** Enable verbose output for debugging and progress tracking */
	public static boolean verbose=false;
	
}

