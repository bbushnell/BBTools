package synth;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sort.Shuffle;
import sort.Shuffle.ShuffleThread;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;


/**
 * Generates artificial cross-contaminated data by mixing reads.
 * Takes input from multiple files, and writes output to the same number of files.
 * @author Brian Bushnell
 * @date Oct 27, 2014
 *
 */
public class CrossContaminate {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Program entry point for cross-contamination simulation.
	 * Creates CrossContaminate instance, processes files, and handles cleanup.
	 * @param args Command-line arguments for input/output files and parameters
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		CrossContaminate x=new CrossContaminate(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs CrossContaminate instance and parses command-line arguments.
	 * Configures input/output files, contamination parameters, and processing options.
	 * Validates file accessibility and parameter ranges before processing begins.
	 *
	 * @param args Command-line arguments containing file paths and parameters
	 * @throws RuntimeException if input files are missing or output files cannot be written
	 */
	public CrossContaminate(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean setInterleaved=false; //Whether it was explicitly set.
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		

		ArrayList<String> inTemp=new ArrayList<String>();
		ArrayList<String> outTemp=new ArrayList<String>();
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(parser.parseCommon(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("in")){
				assert(b!=null) : "Bad parameter: "+arg;
				String[] split2=b.split(",");
				for(String name : split2){
					inNames.add(name);
				}
			}else if(a.equals("out")){
				assert(b!=null) : "Bad parameter: "+arg;
				String[] split2=b.split(",");
				for(String name : split2){
					outNames.add(name);
				}
			}else if(a.equals("innamefile")){
				assert(b!=null) : "Bad parameter: "+arg;
				String[] split2=b.split(",");
				for(String name : split2){
					inTemp.add(name);
				}
			}else if(a.equals("outnamefile")){
				assert(b!=null) : "Bad parameter: "+arg;
				String[] split2=b.split(",");
				for(String name : split2){
					outTemp.add(name);
				}
			}else if(a.equals("shuffle")){
				shuffle=Parse.parseBoolean(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(a.equals("minsinks") || a.equals("ns")){
				minSinks=Integer.parseInt(b);
			}else if(a.equals("maxsinks") || a.equals("xs")){
				maxSinks=Integer.parseInt(b);
			}else if(a.equals("minprob") || a.equals("np")){
				minProb=Double.parseDouble(b);
			}else if(a.equals("maxprob") || a.equals("xp")){
				maxProb=Double.parseDouble(b);
			}else if(a.equals("showspeed")){
				showspeed=Parse.parseBoolean(b);
			}else if(a.equals("shufflethreads")){
				shufflethreads=Integer.parseInt(b);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=parser.overwrite;
			append=parser.append;

			setInterleaved=parser.setInterleaved;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		jgi.DecontaminateByNormalization.parseStringsFromFiles(inTemp);
		jgi.DecontaminateByNormalization.parseStringsFromFiles(outTemp);
		
		inNames.addAll(inTemp);
		outNames.addAll(outTemp);
		inTemp=outTemp=null;
		
		if(inNames.isEmpty() || inNames.size()!=outNames.size()){
			assert(false) : inNames+"\n"+outNames;
			throw new RuntimeException("Error - at least one input file is required, and # input files must equal # output files.");
		}
		
		assert(minSinks<=maxSinks);
		minSinks=Tools.max(0, minSinks);
		maxSinks=Tools.min(inNames.size()-1, maxSinks);
		assert(minSinks<=maxSinks) : minSinks+", "+maxSinks;
		
		assert(minProb<=maxProb);
		assert(minProb>=0 && maxProb<=1);

		minProbPow=Math.log(minProb);
		maxProbPow=Math.log(maxProb);
		
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(!Tools.testInputFiles(true, true, inNames.toArray(new String[0]))){
			outstream.println(outNames);
			throw new RuntimeException("Can't find some input files:\n"+inNames+"\n");
		}
		
		if(!Tools.testOutputFiles(overwrite, append, false, outNames.toArray(new String[0]))){
			outstream.println(outNames);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files.\n");
		}
		
		if(seed>0){randy.setSeed(seed);}
		
		vessels=makeVessels(outNames);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Main processing method that generates cross-contaminated data.
	 * Processes each input file sequentially, distributes reads to sink files,
	 * optionally shuffles output files, and reports timing statistics.
	 *
	 * @param t Timer for tracking execution time and generating performance reports
	 * @throws RuntimeException if processing fails or output becomes corrupted
	 */
	void process(Timer t){
		
		outstream.println("Processing data.");
		for(int i=0; i<inNames.size(); i++){
			try{
				processOneSource(i);
			}catch(Throwable e){
				System.err.println("Failed to open file "+inNames.get(i)+"\nException:"+e+"\n");
				errorState=true;
			}
		}
		
		for(Vessel v : vessels){
			errorState|=v.close();
		}
		
		if(shuffle){
			shuffle(shufflethreads);
		}
		
		t.stop();
		
		if(showspeed){
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		}
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	void shuffle(final int threads){
		outstream.println("Shuffling output in "+threads+" thread"+(threads==1 ? "." : "s."));
		Shuffle.showSpeed=Shuffle.printClass=false;
		Shuffle.setMaxThreads(threads);
		for(Vessel v : vessels){
			ShuffleThread st=new ShuffleThread(v.fname, null, v.fname, null, Shuffle.SHUFFLE, true);
			st.start();
		}
		Shuffle.waitForFinish();
	}
	
	/**
	 * Processes reads from a single source file and distributes them to sink files.
	 * Opens input stream, assigns contamination sinks, and routes each read
	 * probabilistically to appropriate output files based on contamination model.
	 * @param sourceNum Index of the source file in the input list to process
	 */
	void processOneSource(int sourceNum){
		String fname=inNames.get(sourceNum);
		
		FileFormat ffin=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, true);
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin, null);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		final boolean paired=cris.paired();
		if(verbose){
			if(!ffin.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		}
		
		ArrayList<Vessel> sinks=assignSinks(vessels, sourceNum);
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((r.mate!=null)==paired);
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					{
						readsProcessed++;
						basesProcessed+=initialLength1;
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=initialLength2;
					}
					
					addRead(r1, sinks);
				}

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadWrite.closeStream(cris);
	}
	
	/**
	 * Routes a single read to an appropriate sink file based on probability thresholds.
	 * Uses random sampling to select which vessel receives the read according to
	 * the pre-computed cumulative probability distribution.
	 *
	 * @param r The read to be written to a sink file
	 * @param list Ordered list of vessels with cumulative probability thresholds
	 */
	private void addRead(Read r, ArrayList<Vessel> list){
		double p=randy.nextDouble();
		for(Vessel v : list){
			if(p>=v.prob){
				v.bsw.println(r, true);
				r=null;
				break;
			}
		}
		assert(r==null) : p+"\n"+list;
	}
	
	/**
	 * Creates vessel objects for each output file.
	 * Each vessel manages a ByteStreamWriter for efficient output file writing.
	 * @param strings List of output file paths
	 * @return List of initialized vessels ready for writing
	 */
	private ArrayList<Vessel> makeVessels(ArrayList<String> strings){
		ArrayList<Vessel> list=new ArrayList<Vessel>(strings.size());
		for(String s : strings){
			Vessel v=new Vessel(s, true);
			list.add(v);
		}
		return list;
	}
	
	/**
	 * Assigns sink vessels for a source file with random contamination probabilities.
	 * Randomly selects subset of vessels as contamination targets, assigns each
	 * a probability weight, and creates cumulative probability distribution for
	 * efficient read routing during processing.
	 *
	 * @param list Complete list of available vessel containers
	 * @param sourceNum Index of source file to exclude from contamination targets
	 * @return Ordered list of vessels with cumulative probability thresholds
	 */
	private ArrayList<Vessel> assignSinks(ArrayList<Vessel> list, int sourceNum){
		int potential=list.size()-1;
		assert(potential>=minSinks && maxSinks<=potential) : potential+", "+minSinks+", "+maxSinks;
		int range=maxSinks-minSinks+1;
		
		int sinks=minSinks+(range>0 ? randy.nextInt(range) : 0);
		assert(sinks>=0);
		
		for(Vessel v : list){v.prob=0;}
		ArrayList<Vessel> sinklist=(ArrayList<Vessel>) list.clone();
		list=null;
		Vessel source=sinklist.remove(sourceNum);
		if(verbose || true){
			System.err.println("Source:   \t"+inNames.get(sourceNum));
			System.err.println("Sinks:    \t"+sinks);
		}
		
		while(sinklist.size()>sinks){
			int x=randy.nextInt(sinklist.size());
			sinklist.set(x, sinklist.get(sinklist.size()-1));
			sinklist.remove(sinklist.size()-1);
		}
//		if(verbose){System.err.println("Sinklist:\n"+sinklist);}
		
		{
			double probRange=maxProbPow-minProbPow;
			
			assert(probRange>=0) : minProb+", "+maxProb+", "+minProbPow+", "+maxProbPow+", "+probRange;
			
			double remaining=1.0;
			for(Vessel v : sinklist){
				double c=Math.pow(Math.E, minProbPow+randy.nextDouble()*probRange)*remaining;
				remaining-=c;
				v.prob=c;
			}
			source.prob=remaining;
			sinklist.add(source);
			if(verbose || true){System.err.println("Sinklist:\t"+sinklist+"\n");}
			double d=0;
			for(Vessel v : sinklist){
				d+=v.prob;
				v.prob=d;
			}
//			if(verbose){System.err.println("Sinklist:\t"+sinklist);}
			d=0;
			for(Vessel v : sinklist){
				double temp=v.prob;
				v.prob=d;
				d=temp;
			}
//			if(verbose){System.err.println("Sinklist:\t"+sinklist);}
		}
		Collections.reverse(sinklist);
		assert(sinklist.get(sinklist.size()-1).prob==0.0) : sinklist;
		
//		if(verbose){
//			System.err.println("Sinklist:\t"+sinklist);
//			System.err.println();
//		}
		if(verbose || true){System.err.println();}
		
		return sinklist;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Container for output file management with contamination probability.
	 * Encapsulates file format, stream writer, and probability threshold
	 * for efficient read routing during cross-contamination simulation.
	 */
	private class Vessel{
		
		/**
		 * Creates vessel for output file with configured stream writer.
		 * Initializes file format detection and starts ByteStreamWriter for writing.
		 * @param fname_ Output file path for this vessel
		 * @param allowSubprocess Whether to allow subprocess for file writing
		 */
		public Vessel(String fname_, boolean allowSubprocess){
			fname=fname_;
			ff=FileFormat.testOutput(fname, FileFormat.FASTQ, null, allowSubprocess, overwrite, append, false);
			bsw=new ByteStreamWriter(ff);
			bsw.start();
		}
		
		/** Closes the vessel's output stream and waits for completion.
		 * @return true if an error occurred during closing, false otherwise */
		public boolean close(){
			return bsw.poisonAndWait();
		}
		
		@Override
		public String toString(){
			return fname+", "+Tools.format("%.6f", prob);
		}
		
		/** Output file path for this vessel */
		final String fname;
		/** File format configuration for output writing */
		final FileFormat ff;
		/** Stream writer for efficient output file writing */
		final ByteStreamWriter bsw;
		
		/** Cumulative probability threshold for read assignment */
		double prob;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/*--------------------------------------------------------------*/

	/** List of input file paths to process */
	private ArrayList<String> inNames=new ArrayList<String>();
	/** List of output file paths for contaminated results */
	private ArrayList<String> outNames=new ArrayList<String>();
	
	/** Container objects managing output streams for each file */
	private ArrayList<Vessel> vessels;
	
	/*--------------------------------------------------------------*/

	/** Maximum number of reads to process per file (-1 for unlimited) */
	private long maxReads=-1;
	/** Random seed for reproducible contamination patterns (-1 for random) */
	private long seed=-1;
	
	/** Minimum number of contamination target files per source */
	private int minSinks=1;
	/** Maximum number of contamination target files per source */
	private int maxSinks=8;
	/** Minimum contamination probability for any sink file */
	private double minProb=0.000005;
	/** Maximum contamination probability for any sink file */
	private double maxProb=0.025;

	/** Natural logarithm of minimum probability for exponential distribution */
	private double minProbPow=Math.log(minProb);
	/** Natural logarithm of maximum probability for exponential distribution */
	private double maxProbPow=Math.log(maxProb);
	
//	private double root=3.0;
//
//	private double minProbRoot=Math.pow(minProb, 1/root);
//	private double maxProbRoot=Math.pow(maxProb, 1/root);
	
	/** Random number generator for contamination probability sampling */
	private final Random randy=new Random();
	
	/** Total number of reads processed across all files */
	long readsProcessed=0;
	/** Total number of bases processed across all reads */
	long basesProcessed=0;
	
	/** Number of threads to use for output file shuffling */
	private int shufflethreads=3;
	
	/** Whether to shuffle output files after contamination */
	private boolean shuffle=false;
	/** Whether to display processing speed statistics */
	private boolean showspeed=true;
	
	/*--------------------------------------------------------------*/
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and logging */
	private PrintStream outstream=System.err;
	/** Enable verbose output for debugging and detailed progress reporting */
	public static boolean verbose=false;
	/** Tracks whether any errors occurred during processing */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files instead of overwriting */
	private boolean append=false;
	
}
