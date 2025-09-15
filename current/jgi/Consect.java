package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Generates a consensus from multiple error correction results.
 * 
 * @author Brian Bushnell
 * @date October 25, 2016
 *
 */
public class Consect {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		Consect x=new Consect(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Consect(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false; //Important; should be unset later
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("in") || a.equals("in")){
				assert(b!=null) : "Bad parameter: "+arg;
				in.clear();
				String[] split2=b.split(",");
				for(String s : split2){in.add(s);}
			}else if(a.equals("cq") || a.equals("changequality")){
				changeQuality=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(new File(arg).exists()){
				in.add(arg);
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
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

			out=parser.out1;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		//Ensure there is an input file
		if(in==null || in.size()<3){throw new RuntimeException("\nError - at least three input files are required; one original and two error-corrected.");}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in.toArray(new String[0]))){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Create output FileFormat objects
		ffout=FileFormat.testOutput(out, FileFormat.FASTQ, extout, true, overwrite, append, false);

		//Create input FileFormat objects
		ffin=new FileFormat[in.size()];
		for(int i=0; i<in.size(); i++){
			ffin[i]=FileFormat.testInput(in.get(i), FileFormat.FASTQ, extin, true, true);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Create a read input stream
		final ConcurrentReadInputStream[] crisa=new ConcurrentReadInputStream[ffin.length];
		for(int i=0; i<ffin.length; i++){
			crisa[i]=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin[i], null, null, null);
			crisa[i].start(); //Start the stream
			if(verbose){outstream.println("Started cris "+i);}
			assert(!crisa[i].paired());
		}
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros;
		if(ffout!=null){
			final int buff=4;
			
			ros=ConcurrentReadOutputStream.getStream(ffout, null, qfout, null, buff, null, false);
			ros.start(); //Start the stream
		}else{ros=null;}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read streams
		processInner(crisa, ros);
		
		if(verbose){outstream.println("Finished; closing streams.");}

		//Close the read streams
		for(int i=0; i<crisa.length; i++){
			errorState|=ReadWrite.closeStream(crisa[i]);
		}
		//Close the write streams
		errorState|=ReadWrite.closeStream(ros);
		
		//Report timing and results
		{
			t.stop();

			outstream.println("Errors Corrected:         \t"+corrections);
			outstream.println("Disagreements:            \t"+disagreements);
			outstream.println();

			outstream.println("Reads With Corrections:   \t"+readsWithCorrections);
			outstream.println("Reads With Disagreements: \t"+readsWithDisagreements);
			outstream.println();
			
			outstream.println("Reads Fully Corrected:    \t"+readsFullyCorrected);
			outstream.println("Reads Partly Corrected:   \t"+readsPartlyCorrected);
			outstream.println("Reads Not Corrected:      \t"+readsNotCorrected);
			outstream.println("Reads Error Free:         \t"+readsErrorFree);
			outstream.println();
			
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Iterate through the reads */
	void processInner(final ConcurrentReadInputStream[] crisa, final ConcurrentReadOutputStream ros){
		
		//Do anything necessary prior to processing
		
		{
			@SuppressWarnings("unchecked")
			ArrayList<Read>[] array=new ArrayList[crisa.length];
			long id=0;
			for(int i=0; i<crisa.length; i++){
				ListNum<Read> ln=crisa[i].nextList();
				if(ln!=null){
					if(verbose){outstream.println("Fetched "+(ln.list==null ? 0 : ln.list.size())+" reads.");}
					array[i]=ln.list;
					id=ln.id;
				}else{
					array[i]=null;
				}
			}
			if(verbose){outstream.println("Finished fetching block.");}
			
			//As long as there is a nonempty read list...
			while(array[0]!=null && array[0].size()>0){
				
				ArrayList<Read> readsOut=consensus(array);
				
				//Output reads to the output stream
				if(ros!=null){ros.add(readsOut, id);}
				
				//Notify the input stream that the list was used
				for(ConcurrentReadInputStream cris : crisa){
					cris.returnList(id, false);
					if(verbose){outstream.println("Returned a list.");}
				}
				
				for(int i=0; i<crisa.length; i++){
					ListNum<Read> ln=crisa[i].nextList();
					if(ln!=null){
						if(verbose){outstream.println("Fetched "+(ln.list==null ? 0 : ln.list.size())+" reads.");}
						array[i]=ln.list;
						id=ln.id;
					}else{
						array[i]=null;
					}
				}
				if(verbose){outstream.println("Finished fetching block.");}
			}
			
			//Notify the input stream that the final list was used
			for(ConcurrentReadInputStream cris : crisa){
				cris.returnList(id, true);
				if(verbose){outstream.println("Returned final list.");}
			}
		}
		
		//Do anything necessary after processing
		
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private ArrayList<Read> consensus(ArrayList<Read>[] array){
		Read[] row=new Read[array.length];
		for(int i=0, max=array[0].size(); i<max; i++){
			for(int j=0; j<array.length; j++){
				row[j]=array[j].get(i);
			}
			consensus(row);
		}
		return array[0];
	}
	
	/** The first read becomes the consensus of the other reads, where they all agree. */
	private int consensus(Read[] reads){
		Read original=reads[0];
		byte[] obases=original.bases;
		byte[] oquals=changeQuality ? null : original.quality;
		
		for(Read r : reads){assert(r.id.equals(original.id)) : "\n"+r.id+"\n"+original.id;}
		
		readsProcessed++;
		basesProcessed+=obases.length;
		
		int[] counts=KillSwitch.allocInt1D(4);
		int localCorrections=0;
		int localDisagreements=0;
		for(int pos=0; pos<obases.length; pos++){
			Arrays.fill(counts, 0);
			int sum=0, last=-1, tooShort=0;
			byte qmax=0;
			for(int col=1; col<reads.length; col++){
				Read r=reads[col];
				if(r.length()>pos){
					byte q=(r.quality==null ? 0 : r.quality[pos]);
					byte b=r.bases[pos];
					int x=AminoAcid.baseToNumber[b];
					if(x>=0){
						sum++;
						counts[x]++;
						last=x;
						qmax=Tools.max(q, qmax);
					}
				}else{tooShort++;}
			}
			if(sum>1){
				int maxIndex=Tools.maxIndex(counts);
				int max=counts[maxIndex];
				int x0=AminoAcid.baseToNumber[obases[pos]];
				
//				if(tooShort>0){
//					if(maxIndex!=x0){localDisagreements++;}
//				}else
				
				if(max==sum && tooShort==0){
					if(maxIndex!=x0){localCorrections++;}
					obases[pos]=AminoAcid.numberToBase[maxIndex];
					if(oquals!=null){
						oquals[pos]=qmax;
					}
				}else{
					localDisagreements++;
				}
			}
		}
		
		disagreements+=localDisagreements;
		corrections+=localCorrections;
		
		if(localDisagreements+localCorrections>0){
			
			if(localCorrections>0){
				readsWithCorrections++;
			}
			
			if(localDisagreements>0){
				readsWithDisagreements++;
				if(localCorrections>0){
					readsPartlyCorrected++;
				}else{
					readsNotCorrected++;
				}
			}else{
				readsFullyCorrected++;
			}
		}else{
			readsErrorFree++;
		}
		
		return localCorrections;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Input file paths */
	private ArrayList<String> in=new ArrayList<String>();

	/** Output file path */
	private String out=null;
	
	private String qfout=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	private boolean changeQuality=false;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	private long readsFullyCorrected=0;
	private long readsPartlyCorrected=0;
	private long readsNotCorrected=0;
	private long readsErrorFree=0;
	private long readsWithDisagreements=0;
	private long readsWithCorrections=0;
	
	private long disagreements=0;
	private long corrections=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat[] ffin;
	
	/** Primary output file */
	private final FileFormat ffout;
	
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
	/** This flag has no effect on singlethreaded programs */
	private final boolean ordered=false;
	
}
