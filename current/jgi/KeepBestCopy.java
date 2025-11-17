package jgi;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tax.GiToTaxid;
import tracker.ReadStats;

/**
 * Designed to keep the best copy of an SSU per organism.
 * @author Brian Bushnell
 * @date Oct 4, 2019
 *
 */
public class KeepBestCopy {
	
	/**
	 * Program entry point.
	 * Creates instance, processes input, and closes streams.
	 * @param args Command-line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();

		//Create an instance of this class
		KeepBestCopy x=new KeepBestCopy(args);

		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs KeepBestCopy instance with command-line arguments.
	 * Parses arguments, initializes file formats, validates input/output paths.
	 * Sets up shared configuration for threading and compression.
	 * @param args Command-line arguments including input/output files and options
	 */
	public KeepBestCopy(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		Shared.capBuffers(4); //Only for singlethreaded programs
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("maxlen")){
				maxLen=Integer.parseInt(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}
			
			else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
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
			qfin1=parser.qfin1;

			out1=parser.out1;
			qfout1=parser.qfout1;
			
			extin=parser.extin;
			extout=parser.extout;
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
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, extin, true, true);
	}
	
	/**
	 * Creates and starts a concurrent read input stream for processing.
	 * Configures stream with maximum read limit and input file format.
	 * @return Started ConcurrentReadInputStream ready for reading
	 */
	ConcurrentReadInputStream makeCris(){
		final ConcurrentReadInputStream cris;
		cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null, qfin1, null);
		cris.start();
		if(verbose){outstream.println("Started cris");}
		return cris;
	}
	
	/**
	 * Main processing method that reads input, processes sequences, and writes output.
	 * Maintains map of best copies per taxonomic ID and outputs final results.
	 * Tracks statistics for reads and bases processed versus output.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		final ConcurrentReadInputStream cris=makeCris();
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		long readsProcessed=0, readsOut=0;
		long basesProcessed=0, basesOut=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				final ArrayList<Read> listOut=new ArrayList<Read>(reads.size());
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					final int initialLength1=r1.length();
					final boolean keep=process(r1);
					
					readsProcessed++;
					basesProcessed+=initialLength1;
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

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=4;

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, qfout1, null, buff, null, false);
			ros.start();
		}else{ros=null;}
		
		{
			ArrayList<Read> list=new ArrayList<Read>(200);
			long ln=0;
			for(Entry<Integer, Read> e : map.entrySet()){
				Read r=e.getValue();
				list.add(r);
				readsOut++;
				basesOut+=r.length();
				if(list.size()>=200){
					if(ros!=null){ros.add(list, ln);}
					ln++;
					list=new ArrayList<Read>(200);
				}
			}
		}
		
		errorState|=ReadStats.writeAll();
		
		errorState|=ReadWrite.closeStream(ros);
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Processes a single read to determine if it should be kept.
	 * Extracts taxonomic ID from read header and compares against stored copy.
	 * Updates map if this is the first or best copy for the taxonomic ID.
	 *
	 * @param r The read to process
	 * @return true if the read was stored as the best copy, false otherwise
	 */
	private boolean process(Read r){
		int tid=GiToTaxid.parseTaxidNumber(r.id, '|');
		if(tid<0){return false;}
		Integer key=tid;
		Read old=map.get(key);
		if(old==null || isBetterThan(r, old)){
			map.put(key, r);
			return true;
		}
		return false;
	}
	
	/**
	 * Compares two reads to determine which is higher quality.
	 * Prioritizes shorter reads when both exceed maxLen threshold.
	 * Otherwise compares by defined base count (non-N bases) and total N count.
	 *
	 * @param r The candidate read to evaluate
	 * @param old The current stored read (may be null)
	 * @return true if r is better quality than old, false otherwise
	 */
	private boolean isBetterThan(Read r, Read old){
		if(old==null){return true;}
		int oldNs=r.countNocalls();
		int Ns=r.countUndefined();
		int oldDef=old.length()-oldNs;
		int def=r.length()-Ns;
		if(old.length()>maxLen && r.length()<old.length()){return true;}
		if(r.length()>maxLen && old.length()<r.length()){return false;}
		return def>oldDef || (def==oldDef && Ns<oldNs);
	}
	
	/*--------------------------------------------------------------*/
	
	/** Primary input file path */
	private String in1=null;
	
	/** Quality file path for input (if separate from sequence file) */
	private String qfin1=null;

	/** Primary output file path */
	private String out1=null;

	/** Quality file path for output (if separate from sequence file) */
	private String qfout1=null;
	
	/** Override input file extension for format detection */
	private String extin=null;
	/** Override output file extension for format detection */
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	/** Maximum sequence length threshold for quality comparison logic */
	int maxLen=1600;
	
	/** Maximum number of reads to process (-1 for unlimited) */
	private long maxReads=-1;
	
	/** Map storing best copy of each sequence by taxonomic ID */
	private LinkedHashMap<Integer, Read> map=new LinkedHashMap<Integer, Read>();
	
	/*--------------------------------------------------------------*/
	
	/** File format descriptor for primary input file */
	private final FileFormat ffin1;
	/** File format descriptor for primary output file */
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and logging */
	private PrintStream outstream=System.err;
	/** Enables verbose logging throughout the processing pipeline */
	public static boolean verbose=false;
	/** Tracks whether any errors occurred during processing */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files instead of overwriting */
	private boolean append=false;
	
}
