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
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FastaReadInputStream;
import tracker.ReadStats;
import ukmer.KmerTableSetU;

/**
 * @author Brian Bushnell
 * @date Nov 22, 2013
 *
 */
public class KmerFilterSetMaker {
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer(), t2=new Timer();
		t.start();
		t2.start();
		
		//Create a new CountKmersExact instance
		KmerFilterSetMaker x=new KmerFilterSetMaker(args);
		t2.stop();
//		outstream.println("Initialization Time:      \t"+t2);
		
		///And run it
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
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
	
	
	public void process(Timer t){
		
		/* Check for output file collisions */
		Tools.testOutputFiles(overwrite, append, false, kmerOutFile);
		
		/* Count kmers */
		process2();
		
		/* Stop timer and calculate speed statistics */
		t.stop();
		
		/* Throw an exception if errors were detected */
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	public void process2(){
		
		/* Start phase timer */
		Timer t=new Timer();
		
		AbstractKmerTableSet.DISPLAY_STATS=false;
		
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
	
	int runAllPasses(String initialInputFile, String tempPattern){
		//Handle initial input set
		if(initialKmerFile==null){
			AbstractKmerTableSet.overwrite=true;
			AbstractKmerTableSet.append=false;
		}else{
			AbstractKmerTableSet.overwrite=false;
			AbstractKmerTableSet.append=true;
			
			boolean same=initialKmerFile.equals(kmerOutFile) || new File(initialKmerFile).equals(new File(kmerOutFile));
			ReformatReads x=new ReformatReads(new String[] {
					"in="+initialKmerFile, "out="+(same ? null : kmerOutFile), "ow", "silent"});
			x.process(new Timer());
			initialSetSize=x.readsProcessed;
		}

		//Process primary input
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
		
		//Read input file
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
		
		//Summarize counts
		long[] counts=tables.fillHistogram(lastMaxSeen);
		int max=0;
		for(int i=counts.length-1; i>=1; i--){
			if(counts[i]>0){
				max=i;
				break;
			}
		}
		System.err.print(max+" max depth \t");
		
		//Determine minimum count to retain
		int numGood=0;
		int minCountToKeep=-1;
		for(int i=max; i>=1 && numGood<minKmersPerIteration; i--){
			if(i==1 && numGood>0){break;}
			numGood+=counts[i];
			minCountToKeep=i;
		}
		System.err.print(numGood+" high kmers \t");
		
		if(numGood<1){
			System.err.println(0+" retained");
			return -1;
		}
		
		final int maxToKeep=(int)Tools.min(maxKmersPerIteration, tables.readsIn);
		final int numKept=Tools.min(maxToKeep, numGood);
//		final long outSize=initialSetSize+kmersOut+numKept;
//		System.err.println("\nmaxToKeep="+maxToKeep+", numKept="+numKept+", numGood="+numGood+
//				", maxKmersPerIteration="+maxKmersPerIteration+", tables.readsIn="+tables.readsIn);

		//Append retained kmers to a file
		if(numKept>0){
			tables.dumpKmersAsBytes_MT(kmerFile, minCountToKeep, Integer.MAX_VALUE, false, new AtomicLong(numKept));
//			tables.dumpKmersAsBytes(kmerFile, minCountToKeep, Integer.MAX_VALUE, false, new AtomicLong(numKept));
		}
		
		System.err.println(numKept+" retained");
		
		kmersOut+=numKept;
		
		//Filter to retain only unmatched reads
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
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
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
	
	/** Kmer count output file */
	private String kmerOutFile=null;

	private String inFile=null;
	private String initialKmerFile=null;
	private String outTemp=null;
	private String tempKmerFile=makeTempFile("ktemp",".fa");
	
	private boolean errorState=false;
	
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
	
	/** Print messages to this stream */
	private static PrintStream outstream=System.err;
	/** Permission to overwrite existing files */
	public static boolean overwrite=true;
	/** Permission to append to existing files */
	public static boolean append=false;
	/** Display progress messages such as memory usage */
	public static boolean DISPLAY_PROGRESS=false;
	/** Verbose messages */
	public static final boolean verbose=false;
	/** Number of ProcessThreads */
	public static int THREADS=Shared.threads();
	
}
