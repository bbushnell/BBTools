package clump;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import bloom.KCountArray;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
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
import structures.Quantizer;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date June 20, 2014
 *
 */
public class KmerSplit {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		final boolean pigz=ReadWrite.USE_PIGZ, unpigz=ReadWrite.USE_UNPIGZ;
		final boolean oldFInt=FASTQ.FORCE_INTERLEAVED, oldTInt=FASTQ.TEST_INTERLEAVED;
		final int zl=ReadWrite.ZIPLEVEL;
		final float ztd=ReadWrite.ZIP_THREAD_MULT();
		final int mzt=ReadWrite.MAX_ZIP_THREADS();
		Timer t=new Timer();
		KmerSplit x=new KmerSplit(args);
		ReadWrite.ZIPLEVEL=Tools.min(ReadWrite.ZIPLEVEL, maxZipLevel);
		x.process(t);
		ReadWrite.USE_PIGZ=pigz;
		ReadWrite.USE_UNPIGZ=unpigz;
		ReadWrite.ZIPLEVEL=zl;
		ReadWrite.setZipThreadMult(ztd);
		ReadWrite.setZipThreads(mzt);
		FASTQ.FORCE_INTERLEAVED=oldFInt;
		FASTQ.TEST_INTERLEAVED=oldTInt;
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public KmerSplit(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=false;
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		boolean setInterleaved=false; //Whether it was explicitly set.
		Parser parser=new Parser();
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=KmerComparator.verbose=Parse.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
				assert(k>0 && k<32);
			}else if(a.equals("mincount") || a.equals("mincr")){
				minCount=Integer.parseInt(b);
			}else if(a.equals("groups") || a.equals("g") || a.equals("sets") || a.equals("ways")){
				groups=Integer.parseInt(b);
			}else if(a.equals("rename") || a.equals("addname")){
				//Do nothing
//				addName=Parse.parseBoolean(b);
			}else if(a.equals("shortname") || a.equals("shortnames")){
				if(b!=null && b.equals("shrink")){
					shrinkName=true;
				}else{
					shrinkName=false;
					shortName=Parse.parseBoolean(b);
				}
			}else if(a.equals("rcomp") || a.equals("reversecomplement")){
				//ignore rcomp=Parse.parseBoolean(b);
			}else if(a.equals("condense") || a.equals("consensus") || a.equals("concensus")){//Note the last one is intentionally misspelled
				//ignore
			}else if(a.equals("correct") || a.equals("ecc")){
				//ignore
			}else if(a.equals("passes")){
				int x=Integer.parseInt(b);
//				if(x>1){outstream.println("Warning: KmerSplit does not support multiple passes.");}
			}
			
			else if(a.equals("dedupe")){
				//ignore
			}else if(a.equals("entryfilter")){
				//ignore
			}else if(a.equals("markduplicates")){
				//ignore
			}else if(a.equals("markall")){
				//ignore
			}else if(a.equals("addcount") || a.equals("renamebycount")){
				//ignore
			}else if(a.equals("optical") || a.equals("opticalonly")){
				//ignore
			}else if(a.equals("dupesubs") || a.equals("duplicatesubs") || a.equals("dsubs") || a.equals("subs") || a.equals("s")){
				//ignore
			}else if(a.equals("dupedist") || a.equals("duplicatedistance") || a.equals("ddist") || a.equals("dist") || a.equals("opticaldist") || a.equals("distance")){
				//ignore
			}else if(a.equals("scanlimit") || a.equals("scan")){
				//ignore
			}else if(a.equals("removeallduplicates") || a.equals("allduplicates")){
				//ignore
			}else if(a.equals("allowns")){
				//ignore
			}else if(a.equals("containment") || a.equals("absorbcontainment") || a.equals("ac") || a.equals("contains")){
				//ignore
			}else if(a.equalsIgnoreCase("prefixOrSuffix") || a.equalsIgnoreCase("suffixOrPrefix") || a.equals("affix") || a.equals("pos")){
				//ignore
			}else if(a.equals("printduplicates")){
				//ignore
			}else if(a.equals("dupeidentity")){
				//ignore
			}else if(a.equals("dupesubrate") || a.equals("dsr") || a.equals("subrate")){
				//ignore
			}
			
			else if(a.equals("prefilter")){
				KmerReduce.prefilter=Parse.parseBoolean(b);
			}else if(a.equals("ecco")){
				ecco=Parse.parseBoolean(b);
			}else if(a.equals("seed")){
				KmerComparator.defaultSeed=Long.parseLong(b);
			}else if(a.equals("hashes")){
				KmerComparator.setHashes(Integer.parseInt(b));
			}else if(a.equals("border")){
				KmerComparator.defaultBorder=Integer.parseInt(b);
			}else if(a.equals("minprob")){
				KmerComparator.minProb=Float.parseFloat(b);
			}else if(a.equals("unpair")){
				unpair=Parse.parseBoolean(b);
			}else if(a.equals("repair")){
				//Do nothing
			}else if(a.equals("namesort") || a.equals("sort")){
				//Do nothing
			}else if(a.equals("fetchthreads")){
				//Do nothing
			}else if(a.equals("reorder") || a.equals("reorderclumps")){
				//reorder=Parse.parseBoolean(b);
			}else if(a.equals("reorderpaired") || a.equals("reorderclumpspaired")){
//				reorderpaired=Parse.parseBoolean(b);
			}
			
			
			else if(Clump.parseStatic(arg, a, b)){
				//Do nothing
			}
				
			else{
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

			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;

			out1=parser.out1;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		if(groups>2){ReadWrite.USE_PIGZ=false;}
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(!setInterleaved){
			assert(in1!=null) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(out1!=null){
			assert(out1.contains("%"));
			outArray=new String[groups];
			for(int i=0; i<groups; i++){
				outArray[i]=out1.replaceFirst("%", ""+i);
			}
			if(!Tools.testOutputFiles(overwrite, append, false, outArray)){
				outstream.println((out1==null)+", "+out1);
				throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
			}
			ffout=new FileFormat[groups];
			if(groups>1){ReadWrite.setZipThreadMult(Tools.min(0.5f, 2f/(groups+1)));}
			for(int i=0; i<groups; i++){
				ffout[i]=FileFormat.testOutput(outArray[i], FileFormat.FASTQ, extout, groups<10, overwrite, append, false);
			}
		}else{
			outArray=null;
			throw new RuntimeException("out is a required parameter.");
		}

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Count kmers */
	void preprocess(){
		if(minCount>1){
			table=ClumpTools.getTable(in1, in2, k, minCount);
		}
	}

	/** Create read streams and process all data */
	void process(Timer t){
		
		preprocess();
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			cris.start();
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		if(cris.paired() && (in1==null || !in1.contains(".sam") && !unpair)){
			outstream.println("Writing interleaved.");
		}

		final ConcurrentReadOutputStream ros[]=new ConcurrentReadOutputStream[groups];
		try {
			for(int i=0; i<groups; i++){
				final int buff=8;

				assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
				
				ros[i]=ConcurrentReadOutputStream.getStream(ffout[i], null, null, null, buff, null, false);
				ros[i].start();
			}
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
		
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read stream
		processInner(cris, ros);
		
		errorState|=ReadStats.writeAll();
		
		t.stop();
		
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		if(errorState){
			Clumpify.sharedErrorState=true;
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Collect and sort the reads */
	void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream[] ros){
		if(verbose){outstream.println("Making comparator.");}
		KmerComparator kc=new KmerComparator(k, false, false);
		if(verbose){outstream.println("Seed: "+kc.seed);}
		
		if(verbose){outstream.println("Splitting reads.");}
		splitReads(cris, ros, kc);
		lastMemProcessed=memProcessed;
		
		if(verbose){outstream.println("Done!");}
	}
	
	public void splitReads(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream[] ros, final KmerComparator kc){
		Timer t=new Timer();
		if(verbose){t.start("Making hash threads.");}
		final int threads=Shared.threads();
		ArrayList<HashThread> alht=new ArrayList<HashThread>(threads);
		for(int i=0; i<threads; i++){alht.add(new HashThread(i, cris, ros, kc));}
		
		if(verbose){outstream.println("Starting threads.");}
		for(HashThread ht : alht){ht.start();}
		
		
		if(verbose){outstream.println("Waiting for threads.");}
		/* Wait for threads to die */
		for(HashThread ht : alht){
			
			/* Wait for a thread to die */
			while(ht.getState()!=Thread.State.TERMINATED){
				try {
					ht.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			readsProcessed+=ht.readsProcessedT;
			basesProcessed+=ht.basesProcessedT;
			diskProcessed+=ht.diskProcessedT;
			memProcessed+=ht.memProcessedT;
		}
		
		if(verbose){outstream.println("Closing streams.");}
		errorState=ReadWrite.closeStreams(cris, ros)|errorState;
		if(verbose){t.stop("Split time: ");}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class HashThread extends Thread{

		HashThread(int id_, ConcurrentReadInputStream cris_, ConcurrentReadOutputStream[] ros_, KmerComparator kc_){
			id=id_;
			cris=cris_;
			ros=ros_;
			kc=kc_;
		}

		@Override
		public void run(){

			final boolean paired=cris.paired();
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			ArrayList<Read>[] array=new ArrayList[groups];
			for(int i=0; i<groups; i++){
				array[i]=new ArrayList<Read>(buffer);
			}
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				for(Read r : reads){
					if(!r.validated()){
						r.validate(true);
						if(r.mate!=null){r.mate.validate(true);}
					}
					readsProcessedT+=1+r.mateCount();
					basesProcessedT+=r.length()+r.mateLength();
					diskProcessedT+=r.countFastqBytes()+r.countMateFastqBytes();
					memProcessedT+=r.countBytes()+r.countMateBytes()+ReadKey.overhead;
					if(shrinkName){
						Clumpify.shrinkName(r);
						Clumpify.shrinkName(r.mate);
					}else if(shortName){
						Clumpify.shortName(r);
						Clumpify.shortName(r.mate);
					}
					
					if(quantizeQuality){
						Quantizer.quantize(r, r.mate);
					}
				}
				
				if(ecco){
					for(Read r : reads){
						if(r.mate!=null){BBMerge.findOverlapStrict(r, r.mate, true);}
					}
				}
				
				ArrayList<Read> hashList=reads;
				if(paired && unpair){
					hashList=new ArrayList<Read>(reads.size()*2);
					for(Read r1 : reads){
						Read r2=r1.mate;
						hashList.add(r1);
						hashList.add(r2);
						r1.mate=null;
						r2.mate=null;
					}
				}
				
				kc.hash(hashList, table, minCount, true);
				for(Read r : hashList){
					long kmer=((ReadKey)r.obj).kmer;
					long code=kc.hash(kmer);
					int code2=(int)(code%groups);
					assert(code2>=0 && code2<array.length) : code2+", "+groups+", "+array.length+", "+kmer+", "+r.obj+"\n"+r;
					array[code2].add(r);
					if(array[code2].size()>=buffer){
						ros[code2].add(array[code2], 0);
						array[code2]=new ArrayList<Read>(buffer);
					}
				}
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
			for(int i=0; i<groups; i++){
				if(!array[i].isEmpty()){
					ros[i].add(array[i], 0);
				}
			}
		}

		final int id;
		final ConcurrentReadInputStream cris;
		final ConcurrentReadOutputStream[] ros;
		final KmerComparator kc;
		static final int buffer=200;
		
		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
		protected long diskProcessedT=0;
		protected long memProcessedT=0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int k=31;
	int groups=16;
	int minCount=0;
	
	KCountArray table=null;
	
	/*--------------------------------------------------------------*/
	/*----------------          I/O Fields          ----------------*/
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String in2=null;

	private String out1=null;
	private String[] outArray=null;
	
	private String extin=null;
	private String extout=null;
	
	/*--------------------------------------------------------------*/
	
	protected long readsProcessed=0;
	protected long basesProcessed=0;
	protected long diskProcessed=0;
	protected long memProcessed=0;
	
	protected static long lastMemProcessed=0;
	
	private long maxReads=-1;
//	private boolean addName=false;
	boolean shortName=false;
	boolean shrinkName=false;
	boolean ecco=false;
	boolean unpair=false;
	
	static int maxZipLevel=2;

	static boolean quantizeQuality=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;
	
	private final FileFormat[] ffout;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
