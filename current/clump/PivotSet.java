package clump;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import bloom.KCountArray;
import bloom.ReadCounter;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Reduces reads to their feature kmer.
 * @author Brian Bushnell
 * @date August 19, 2016
 *
 */
public class PivotSet {

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		makeSet(args);
	}
	
	public static KCountArray makeSet(String[] args){
		final boolean pigz=ReadWrite.USE_PIGZ, unpigz=ReadWrite.USE_UNPIGZ;
		Timer t=new Timer();
		PivotSet x=new PivotSet(args);
		KCountArray kca=x.process(t, false);
		ReadWrite.USE_PIGZ=pigz;
		ReadWrite.USE_UNPIGZ=unpigz;
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
		
		return kca;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public PivotSet(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
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
			}else if(a.equals("ecco")){
				ecco=Parse.parseBoolean(b);
			}else if(a.equals("rename") || a.equals("addname")){
				//do nothing
			}else if(a.equals("rcomp") || a.equals("reversecomplement")){
				//do nothing
			}else if(a.equals("condense") || a.equals("consensus")){
				//do nothing
			}else if(a.equals("mincount") || a.equals("consensus")){
				minCount=Integer.parseInt(b);
			}else if(a.equals("correct") || a.equals("ecc")){
				//do nothing
			}else if(a.equals("groups") || a.equals("g") || a.equals("sets") || a.equals("ways")){
				//do nothing
			}else if(a.equals("seed")){
				KmerComparator.defaultSeed=Long.parseLong(b);
			}else if(a.equals("hashes")){
				KmerComparator.setHashes(Integer.parseInt(b));
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			in1=parser.in1;
			in2=parser.in2;
			
			extin=parser.extin;
		}
		
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
		
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private static long getCells(double fraction, int cbits){
		final long memory=Runtime.getRuntime().maxMemory();
		final long usable=(long)Tools.max(((memory-96000000)*.73), memory*0.45);
		final double filterMem=usable*fraction;
		return (long)((filterMem*8)/cbits);
	}
	
	/** Create read streams and process all data */
	public KCountArray process(Timer t, boolean amino){
		int cbits=2;
		while((1L<<cbits)<=minCount){cbits*=2;}
		int filterHashes=2;
		float fraction=0.1f;
		long cells=getCells(fraction, cbits);
		ReadCounter rc=new ReadCounter(k, true, ecco, false, amino);
		KCountArray kca=rc.makeKca(null, null, null, cbits, cells, filterHashes, 0, maxReads, 1, 1, 1, null, 0);
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			cris.start();
			if(verbose){outstream.println("Started cris");}
		}
		
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read stream
		processInner(cris, kca);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		errorState|=ReadStats.writeAll();
		errorState|=ReadWrite.closeStreams(cris);
		
		t.stop();

		outstream.println("Made filter:     \t"+kca.toShortString(filterHashes));
		outstream.println("Estimated pivots:      \t"+(long)kca.estimateUniqueKmers(filterHashes));
		outstream.println("Estimated pivots >1x:  \t"+(long)kca.estimateUniqueKmers(filterHashes, minCount));
		
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		if(errorState){
			Clumpify.sharedErrorState=true;
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
		return kca;
	}
	
	/** Manage threads */
	public static KCountArray makeKcaStatic(final ConcurrentReadInputStream cris, int k, int minCount, boolean amino){

		KmerComparator kc=new KmerComparator(k, false, false);
		int cbits=2;
		while((1L<<cbits)<=minCount){cbits*=2;}
		int filterHashes=2;
		float fraction=0.1f;
		long cells=getCells(fraction, cbits);
		ReadCounter rc=new ReadCounter(k, true, false, false, amino);
		KCountArray kca=rc.makeKca(null, null, null, cbits, cells, filterHashes, 0, -1, 1, 1, 1, null, 0);
		
		if(verbose){System.err.println("Making hash threads.");}
		final int threads=Shared.threads();
		ArrayList<HashThread> alht=new ArrayList<HashThread>(threads);
		for(int i=0; i<threads; i++){alht.add(new HashThread(cris, kc, kca, false));}
		
		if(verbose){System.err.println("Starting threads.");}
		for(HashThread ht : alht){ht.start();}
		
		if(verbose){System.err.println("Waiting for threads.");}
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
		}
		kca.shutdown();
		return kca;
	}
	
	/** Manage threads */
	public void processInner(final ConcurrentReadInputStream cris, KCountArray kca){
		if(verbose){outstream.println("Making comparator.");}
		KmerComparator kc=new KmerComparator(k, false, false);
		
		if(verbose){outstream.println("Making hash threads.");}
		final int threads=Shared.threads();
		ArrayList<HashThread> alht=new ArrayList<HashThread>(threads);
		for(int i=0; i<threads; i++){alht.add(new HashThread(cris, kc, kca, ecco));}
		
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
		}
		kca.shutdown();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static class HashThread extends Thread{

		HashThread(ConcurrentReadInputStream cris_, KmerComparator kc_, KCountArray kca_, boolean ecco_){
			cris=cris_;
			kc=kc_;
			kca=kca_;
			ecco=ecco_;
		}

		@Override
		public void run(){
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				for(Read r1 : reads){
					Read r2=r1.mate;
					readsProcessedT+=r1.pairCount();
					basesProcessedT+=r1.pairLength();
					if(ecco && r2!=null){
						if(r2!=null){BBMerge.findOverlapStrict(r1, r2, true);}
					}
					{
						final long kmer=kc.hash(r1, null, 0, false);
						if(kmer>=0){
							kca.increment(kmer);
						}
					}
					if(r2!=null){
						final long kmer=kc.hash(r2, null, 0, false);
						if(kmer>=0){
							kca.increment(kmer);
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
		
		final ConcurrentReadInputStream cris;
		final KmerComparator kc;
		final KCountArray kca;
		final boolean ecco;
		
		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int k=31;
	private int minCount=2;
	
	/*--------------------------------------------------------------*/
	/*----------------          I/O Fields          ----------------*/
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String in2=null;
	
	private String extin=null;
	
	/*--------------------------------------------------------------*/
	
	protected long readsProcessed=0;
	protected long basesProcessed=0;
	
	private long maxReads=-1;
	private boolean ecco=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	
}
