package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;

import assemble.Tadpole;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import kmer.AbstractKmerTableSet;
import kmer.KmerTableSet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FastaReadInputStream;
import tracker.ReadStats;
import ukmer.KmerTableSetU;

/**
 * Builds a filtered k-mer set by repeatedly counting the current input,
 * requesting a bounded quota of retained k-mers to dump to the output file,
 * and filtering unmatched reads into the next pass's temporary input.
 * @author Brian Bushnell
 * @date Nov 22, 2013
 */
public class KmerFilterSetMaker{

	/**
	 * Creates and runs the k-mer filter workflow.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer(), t2=new Timer();
		t.start();
		t2.start();

		//Create a new filter instance.
		KmerFilterSetMaker x=new KmerFilterSetMaker(args);
		t2.stop();
//		outstream.println("Initialization Time:      \t"+t2);

		//Run the count and filter passes.
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}

	/**
	 * Parses local output, threshold, pass-count, length, and forwarded
	 * {@link KmerTableSet} arguments, then records the requested k-mer length.
	 * <p>The FASTA reader settings check is assertion-only; this constructor
	 * also changes global preprocessing and I/O flags without restoration.</p>
	 * @param args command line arguments
	 */
	public KmerFilterSetMaker(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		PreParser.silent=true;

		/* Set global defaults */
		ReadWrite.ZIPLEVEL=2;
		ReadWrite.USE_UNPIGZ=true;

		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}

		/* Initialize local variables with defaults */
		boolean useForest_=false, useTable_=false, useArray_=true;
		Parser parser=new Parser();

		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=Parse.splitOnFirst(arg, '='); // JGI-258: Preserve value text after the first '='.
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
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(parser.parseTrim(arg, a, b)){
				//do nothing
			}else if(a.equals("out") || a.equals("out1") || a.equals("outkmers") || a.equals("outk") || a.equals("dump")){
				kmerOutFile=b;
			}else if(a.equals("in") || a.equals("in1")){
				inFile=b;
			}else if(a.equals("initial") || a.equals("starting") || a.equals("initialset") || a.equals("startingset")){
				initialKmerFile=b;
			}else if(a.equals("temp") || a.equals("pattern")){
				outTemp=b;
				assert(outTemp.contains("#"));
			}else if(a.equals("mincount") || a.equals("min")){
				minCount=Integer.parseInt(b);
			}else if(a.equals("maxns")){
				maxNs=Integer.parseInt(b);
				AbstractKmerTableSet.maxNs=maxNs;
			}else if(a.equals("minlen")){
				minLen=Integer.parseInt(b);
				AbstractKmerTableSet.minLen=minLen;
			}else if(a.equals("passes") || a.equals("maxpasses")){
				maxPasses=Integer.parseInt(b);
			}else if(a.equals("minkmers") || a.equals("minkmersperpass") || a.equals("minkpp")){
				minKmersPerIteration=Integer.parseInt(b);
			}else if(a.equals("maxkmers") || a.equals("maxkmersperpass") || a.equals("maxkpp")){
				maxKmersPerIteration=Integer.parseInt(b);
				if(maxKmersPerIteration<1){
					maxKmersPerIteration=Integer.MAX_VALUE;
				}
			}
//			else if(a.equals("maxcounttodump") || a.equals("maxdump") || a.equals("maxcount")){
//				maxToDump=Integer.parseInt(b);
//			}

			else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=(b==null || b.equalsIgnoreCase("auto") ? Shared.threads() : Integer.parseInt(b));
			}else if(a.equals("verbose")){
				assert(false) : "Verbose flag is currently static final; must be recompiled to change.";
//				verbose=Parse.parseBoolean(b);
			}

			else if(KmerTableSet.isValidArgument(a)){
				tableArgs.add(arg);
			}

			else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		{//Process parser fields
			Parser.processQuality();
		}

		if(outTemp==null){
			outTemp=makeTempFile("rtemp_#", ReadWrite.getExtension(inFile));
		}

		/* Adjust I/O settings and filenames */

		assert(FastaReadInputStream.settingsOK());

		assert(kmerOutFile!=null) : "Kmer output file is required.";
		if(kmerOutFile!=null && !Tools.canWrite(kmerOutFile, overwrite)){throw new RuntimeException("Output file "+kmerOutFile+" already exists, and overwrite="+overwrite);}

		assert(THREADS>0);

		k=Tadpole.preparseK(args);
	}


	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/


	/**
	 * Runs all count/filter passes and throws if a k-mer dump error was latched.
	 * Filtering and summary output precede this final error check.
	 * @param t process timer stopped after the passes finish
	 */
	public void process(Timer t){

		/* Check for output file collisions; the return value is not consumed. */
		Tools.testOutputFiles(overwrite, append, false, kmerOutFile);

		/* Run the count, dump, and filter passes. */
		process2();

		/* Stop timer and calculate speed statistics */
		t.stop();

		/* Throw an exception if errors were detected */
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}


	/**
	 * Disables table statistics, runs all passes, then prints accumulated
	 * first-pass input counters and accumulated requested output quotas.
	 */
	public void process2(){

		/* Start phase timer */
		Timer t=new Timer();

		AbstractKmerTableSet.DISPLAY_STATS=false;

		/* Run count/filter passes. */
		runAllPasses(inFile, outTemp);

		t.stop();
		outstream.println();
		outstream.println("Input:                      \t"+readsIn+" reads \t\t"+basesIn+" bases.");
		outstream.println("Output:                     \t"+kmersOut+" kmers.");
		outstream.println("Passes:                     \t"+numPasses);

//		if(shave || rinse){
//			kmersRemoved=shave(shave, rinse, shaveDepth);
//		}
//		
//		outstream.println("\nFor K="+tables.kbig());
//		outstream.println("Unique Kmers:               \t"+tables.kmersLoaded);
//		if(shave || rinse){
//			outstream.println("After Shaving:              \t"+(tables.kmersLoaded-kmersRemoved));
//		}

		outstream.println();

		outstream.println("Time:                       \t"+t);
	}

	/**
	 * Copies a differently named initial k-mer set, then repeats primary
	 * counting and filtering until {@link #maxPasses} is reached or the selected
	 * cutoff is less than {@link #minCount}.  After each pass beyond the first,
	 * attempts to delete the input for that pass.  After normal loop completion,
	 * attempts to delete the last output, if any.  Deletion results are ignored
	 * and there is no finally cleanup.  Global table overwrite and append flags
	 * are changed without restoration.
	 * @param initialInputFile first primary input file
	 * @param tempPattern temporary read-file pattern containing {@code #}
	 * @return last selected cutoff, {@code -1}, or 20000000 if no pass ran
	 */
	int runAllPasses(String initialInputFile, String tempPattern){
		//Handle initial k-mer set
		if(initialKmerFile==null){
			AbstractKmerTableSet.overwrite=true;
			AbstractKmerTableSet.append=false;
		}else{
			AbstractKmerTableSet.overwrite=false;
			AbstractKmerTableSet.append=true;

			boolean same=initialKmerFile.equals(kmerOutFile) || new File(initialKmerFile).equals(new File(kmerOutFile));
			ReformatReads x=new ReformatReads(new String[]{
					"in="+initialKmerFile, "out="+(same ? null : kmerOutFile), "ow", "silent"});
			x.process(new Timer());
			initialSetSize=x.readsProcessed;
		}

		//Count and filter the primary input.
		int maxCount=20000000; //Initial array size
		String in=initialInputFile;
		String lastOut=null;
		for(int pass=0; pass<maxPasses && maxCount>=minCount; pass++){
			String out=tempPattern.replace("#", ""+pass);
			maxCount=runOnePass(in, out, kmerOutFile, maxCount, pass);
			if(pass>0){new File(in).delete();}
			in=lastOut=out;
			AbstractKmerTableSet.overwrite=false;
			AbstractKmerTableSet.append=true;
		}
		if(lastOut!=null){new File(lastOut).delete();}//This does not actually need to be created most of the time
		return maxCount;
	}

	/**
	 * Counts one input, selects a bounded quota from a histogram whose last bin
	 * may be capped or folded, requests dumping to the accumulated k-mer file,
	 * and filters unmatched reads.  Input counters are updated only on pass 0 and
	 * are not reset across repeated calls.
	 * @param inFile input reads for this pass
	 * @param outFile unmatched-read output for the following pass
	 * @param kmerFile accumulated k-mer output and BBDuk reference
	 * @param lastMaxSeen supplied final-bin cap for the histogram; initially
	 * 20000000 and the previous selected cutoff thereafter
	 * @param pass zero-based pass number
	 * @return selected minimum count to keep, or {@code -1} when the selected
	 * histogram frequency sum is less than 1
	 */
	int runOnePass(String inFile, String outFile, String kmerFile, int lastMaxSeen, int pass){
		Timer t=new Timer();

		@SuppressWarnings("unchecked")
		ArrayList<String> tableArgs2=(ArrayList<String>) tableArgs.clone();
		tableArgs2.add("k="+k);
		tableArgs2.add("in="+inFile);
//		tableArgs2.add("showstats=f");
//		tableArgs2.add("showprogress=f");

		AbstractKmerTableSet.DISPLAY_STATS=AbstractKmerTableSet.DISPLAY_PROGRESS=false;

		System.err.print("Pass "+pass+"  \t");

		//Count the input file.
		AbstractKmerTableSet tables;
		if(k<=31){//TODO: 123 add "false" to the clause to force KmerTableSetU usage.
			tables=new KmerTableSet(tableArgs2.toArray(new String[0]), 12);
		}else{
			tables=new KmerTableSetU(tableArgs2.toArray(new String[0]), 0);
		}
		tables.process(t);
		System.err.print(tables.readsIn+" reads \t"+tables.kmersIn+" kmers \t");
		if(pass==0){
			readsIn+=tables.readsIn;
			basesIn+=tables.basesIn;
			kmersIn+=tables.kmersIn;
		}
		numPasses++;

		//Summarize depth counts.
		long[] counts=tables.fillHistogram(lastMaxSeen);
		int max=0;
		for(int i=counts.length-1; i>=1; i--){
			if(counts[i]>0){
				max=i;
				break;
			}
		}
		System.err.print(max+" max depth \t");

		Selection selection=selectKmers(counts, max, minKmersPerIteration,
				maxKmersPerIteration, tables.readsIn);
		final long numGood=selection.numGood;
		final int minCountToKeep=selection.minCountToKeep;
		System.err.print(numGood+" high kmers \t");

		if(numGood<1){
			System.err.println(0+" retained");
			return -1;
		}

		final int numKept=selection.numKept;
//		final long outSize=initialSetSize+kmersOut+numKept;
//		System.err.println("\nmaxToKeep="+maxToKeep+", numKept="+numKept+", numGood="+numGood+
//				", maxKmersPerIteration="+maxKmersPerIteration+", tables.readsIn="+tables.readsIn);

		//Append the selected quota of retained k-mers to a file.
		if(numKept>0){
			// JGI-259: Preserve true-on-error dump status for the final process check.
			errorState|=tables.dumpKmersAsBytes_MT(kmerFile, minCountToKeep, Integer.MAX_VALUE, false, new AtomicLong(numKept));
//			tables.dumpKmersAsBytes(kmerFile, minCountToKeep, Integer.MAX_VALUE, false, new AtomicLong(numKept));
		}

		System.err.println(numKept+" retained");

		kmersOut+=numKept;

		//Filter to retain only unmatched reads.
		if(minCountToKeep>=1 && tables.readsIn>=1){
			ArrayList<String> bbdukArgs=new ArrayList<String>();
			bbdukArgs.add("in="+inFile);
			bbdukArgs.add("outu="+outFile);
			bbdukArgs.add("ref="+kmerFile); //Technically this can only be the latest kmers, in case this file gets big.
			bbdukArgs.add("k="+k);
			bbdukArgs.add("mm=f");
			bbdukArgs.add("rcomp="+tables.rcomp());
			bbdukArgs.add("silent=t");
			if(maxNs>=0 && maxNs<Integer.MAX_VALUE){bbdukArgs.add("maxns="+maxNs);}
			if(minLen>0){bbdukArgs.add("minlen="+minLen);}
			bbdukArgs.add("ordered");
			BBDuk.main(bbdukArgs.toArray(new String[0]));
		}

		return minCountToKeep;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Scans positive-depth bins downward to choose the lowest retained
	 * depth and a requested output quota.  The maximum can be zero if no positive
	 * bin exists, and depth 1 is skipped after any higher-depth k-mers are seen.
	 * @param counts count-frequency histogram
	 * @param max highest occupied histogram bin
	 * @param minKmersPerIteration minimum requested k-mers per pass
	 * @param maxKmersPerIteration maximum requested k-mers per pass
	 * @param readsIn reads counted by this pass, used to cap the request
	 * @return selected threshold and requested quota
	 */
	private static Selection selectKmers(final long[] counts, final int max,
			final int minKmersPerIteration, final int maxKmersPerIteration,
			final long readsIn){
		assert(counts!=null) : "KmerFilterSetMaker.runOnePass supplies a nonnull fillHistogram result.";
		assert(max>=0 && max<counts.length) : "KmerFilterSetMaker.runOnePass scans max inside the counts array; max="+
				max+", length="+counts.length;
		//Keep histogram frequencies wide until the bounded output quota is selected (JGI-212).
		long numGood=0;
		int minCountToKeep=-1;
		for(int i=max; i>=1 && numGood<minKmersPerIteration; i--){
			if(i==1 && numGood>0){break;}
			numGood+=counts[i];
			minCountToKeep=i;
		}
		final int maxToKeep=(int)Tools.min(maxKmersPerIteration, readsIn);
		final int numKept=(numGood<1 ? 0 : (int)Tools.min(maxToKeep, numGood));
		return new Selection(minCountToKeep, numGood, numKept);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Helper Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Immutable selected threshold and requested quota for one histogram pass. */
	private static final class Selection{

		Selection(final int minCountToKeep_, final long numGood_, final int numKept_){
			minCountToKeep=minCountToKeep_;
			numGood=numGood_;
			numKept=numKept_;
		}

		final int minCountToKeep;
		final long numGood;
		final int numKept;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

//	/** Hold kmers. */
//	private final AbstractKmerTableSet tables;

//	private boolean shave=false;
//	private boolean rinse=false;
//	private int shaveDepth=1;

	private ArrayList<String> tableArgs=new ArrayList<String>();

	private long basesIn=0;
	private long readsIn=0;

	private long kmersIn=0;
	private long kmersOut=0;
	private long numPasses;
	private long initialSetSize=0;

	private int maxPasses=3000;
	private int minCount=1;
	private int minKmersPerIteration=1;
	private int maxKmersPerIteration=2;
	private int maxNs=Integer.MAX_VALUE;
	private int minLen=1;

	private String kmerOutFile=null;

	private String inFile=null;
	private String initialKmerFile=null;
	private String outTemp=null;
	private String tempKmerFile=makeTempFile("ktemp",".fa");

	private boolean errorState=false;

	/**
	 * Generates a temporary path pattern string without creating a file or
	 * guaranteeing the name is unique.
	 * @param prefix file-name prefix
	 * @param ext extension with or without the leading dot
	 * @return generated file-name string
	 */
	private String makeTempFile(String prefix, String ext){
		if(!ext.startsWith(".")){ext="."+ext;}
//		try {
			return prefix+"_"+(System.nanoTime()&0xFFFFF)+"_"+((int)(10000000*(2+Math.random())))+ext;
//			return File.createTempFile("temp#", ".fq").getAbsolutePath();
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		return null;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/

	final int k;

	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/

	private static PrintStream outstream=System.err;
	public static boolean overwrite=true;
	public static boolean append=false;
	public static boolean DISPLAY_PROGRESS=false;
	public static final boolean verbose=false;
	public static int THREADS=Shared.threads();

}
