package cluster;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;

import dna.Data;
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

/**
 * @author Brian Bushnell
 * @date Feb 7, 2014
 *
 */
public class ReclusterByKmer {
	
	/**
	 * Program entry point for k-mer based read clustering.
	 * Initializes the clustering pipeline and processes input reads.
	 * @param args Command-line arguments specifying input files and clustering parameters
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		ReclusterByKmer x=new ReclusterByKmer(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs ReclusterByKmer instance with command-line configuration.
	 * Parses arguments to set up input/output files, k-mer parameters, and clustering modes.
	 * Initializes dual k-mer lengths (k1 and k2) and validates parameter ranges.
	 * @param args Command-line arguments for configuration
	 */
	public ReclusterByKmer(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean setInterleaved=false; //Whether it was explicitly set.
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		
		
		int k1_=12;
		int k2_=3;
		Parser parser=new Parser();
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(parser.parseFiles(arg, a, b)){
				//do nothing
			}else if(parser.parseCommon(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("build") || a.equals("genome")){
				Data.setGenome(Integer.parseInt(b));
			}else if(in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				in1=arg;
			}else{
				System.err.println("Unknown parameter "+args[i]);
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
			
			in1=parser.in1;
			in2=parser.in2;
//			qfin1=parser.qfin1;
//			qfin2=parser.qfin2;

			out1=parser.out1;
			out2=parser.out2;
//			qfout1=parser.qfout1;
//			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){System.err.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
//			if(ReadWrite.isCompressed(in1)){ByteFile.FORCE_MODE_BF2=true;}
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
		
		if(!setInterleaved){
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);

		/* Check for output file collisions */
		Tools.testOutputFiles(overwrite, append, false, out1, out2);

		k1=k1_;
		assert(k1>=-1 && k1<=15) : "k1 must lie between 1 and 15, inclusive (0 to disable)";
		k2=k2_;
		assert(k2>=-1 && k2<=6) : "k2 must lie between 1 and 6, inclusive (0 to disable)";

		arraylen1=(k1>0 ? ClusterTools.maxCanonicalKmer(k1)+1 : 0);
		arraylen2=(k2>0 ? ClusterTools.maxCanonicalKmer(k2)+1 : 0);
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Main processing method that reads input and writes output without clustering.
	 * Currently implements a pass-through mode that processes reads but does not
	 * perform actual clustering operations.
	 * @param t Timer for tracking execution time and performance metrics
	 */
	void process(Timer t){
		
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
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
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) : "out1 and out2 have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, null, null, buff, null, false);
			ros.start();
		}
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			System.err.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					{
						readsProcessed++;
						basesProcessed+=r1.length();
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=r2.length();
					}
					
					
					boolean remove=false;
					if(remove){reads.set(idx, null);}
				}
				
				if(ros!=null){ros.add(reads, ln.id);}

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
			throw new RuntimeException("ReformatReads terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Retrieves or creates a cluster with the specified ID.
	 * Thread-safe method that ensures cluster list capacity and initializes
	 * new clusters as needed with dual k-mer support.
	 *
	 * @param x The cluster ID to fetch or create
	 * @return The cluster instance for the specified ID
	 */
	Cluster fetchCluster(int x){
		if(x>clusterList.size()){
			synchronized(clusterList){
				clusterList.ensureCapacity(2*x);
				for(int i=clusterList.size(); i<x; i++){
					clusterList.add(new Cluster(i, k1, k2, arraylen1, arraylen2));
				}
				clusterList.notifyAll();
			}
		}
		Cluster c=clusterList.get(x);
		while(c==null){
			synchronized(clusterList){
				c=clusterList.get(x);
				assert(c!=null);
			}
		}
		return c;
	}
	
	/**
	 * Creates clusters; Generates kmer spectra for clusters
	 * @param t
	 */
	private void findKmerSpectra(Timer t){
		
		/* Create read input stream */
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2);
			cris.start(); //4567
			if(verbose){System.err.println("Started cris");}
		}
		
		/* Create ClusterThreads */
		ArrayList<ClusterThread> alct=new ArrayList<ClusterThread>(THREADS);
		
		for(int i=0; i<THREADS; i++){alct.add(new ClusterThread(i, CLUSTER_MODE_CREATE, -1, cris));}
		for(ClusterThread ct : alct){ct.start();}

		long readsIn=0, basesIn=0;
		
		/* Wait for threads to die, and gather statistics */
		for(ClusterThread ct : alct){
			while(ct.getState()!=Thread.State.TERMINATED){
				try {
					ct.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			readsIn+=ct.readsInT;
			basesIn+=ct.basesInT;
		}
		
		/* Shut down I/O streams; capture error status */
		errorState|=ReadWrite.closeStreams(cris);
	}
	
	/**
	 * Assign reads to clusters using additional kmer information.
	 * @param t
	 */
	private void recluster(Timer t){
		
		/* Create read input stream */
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2);
			cris.start(); //4567
			if(verbose){System.err.println("Started cris");}
		}
		
		/* Create ClusterThreads */
		ArrayList<ClusterThread> alct=new ArrayList<ClusterThread>(THREADS);
		for(int i=0; i<THREADS; i++){alct.add(new ClusterThread(i, CLUSTER_MODE_RECLUSTER, ambigMode, cris));}
		for(ClusterThread ct : alct){ct.start();}
		
		long readsIn=0, basesIn=0;
		
		/* Wait for threads to die, and gather statistics */
		for(ClusterThread ct : alct){
			while(ct.getState()!=Thread.State.TERMINATED){
				try {
					ct.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			readsIn+=ct.readsInT;
			basesIn+=ct.basesInT;
		}
		
		/* Shut down I/O streams; capture error status */
		errorState|=ReadWrite.closeStreams(cris);
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Worker thread for processing reads in clustering operations.
	 * Handles both cluster creation and reclustering modes with configurable
	 * ambiguity resolution strategies for optimal read assignment.
	 */
	private class ClusterThread extends Thread{
		
		/**
		 * Constructs a ClusterThread for read processing.
		 *
		 * @param id_ Thread identifier for this worker
		 * @param clusterMode_ Clustering operation mode (CREATE or RECLUSTER)
		 * @param ambigMode_ Strategy for handling ambiguous cluster assignments
		 * @param cris_ Input stream for reading sequences
		 */
		public ClusterThread(int id_, int clusterMode_, int ambigMode_, ConcurrentReadInputStream cris_){
			id=id_;
			ambigModeT=ambigMode_;
			clusterMode=clusterMode_;
			cris=cris_;
		}
		
		@Override
		public void run(){
			randy=(ambigModeT==AMBIG_MODE_RAND) ? Shared.threadLocalRandom() : null;
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//While there are more reads lists...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				//For each read (or pair) in the list...
				for(int i=0; i<reads.size(); i++){
					processRead(reads.get(i));
				}
				
				//Fetch a new read list
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln);
		}
		
		
		/**
		 * Processes a single read (and its mate if paired) for clustering.
		 * Creates ReadTag objects with metadata and delegates to appropriate
		 * clustering method based on the thread's operation mode.
		 * @param r1 The read to process (may have a mate in r1.mate)
		 */
		private void processRead(final Read r1){
			final Read r2=r1.mate;

			if(verbose){System.err.println("Considering read "+r1.id+" "+new String(r1.bases));}

			readsInT++;
			basesInT+=r1.length();
			if(r2!=null){
				readsInT++;
				basesInT+=r2.length();
			}

			final ReadTag rt1=new ReadTag(r1);
			final ReadTag rt2=(r2==null ? null : new ReadTag(r2));

			r1.obj=rt1;
			if(r2!=null){r2.obj=rt2;}

			if(clusterMode==CLUSTER_MODE_CREATE){
				addToCluster(r1, r2, rt1, rt2);
			}else if(clusterMode==CLUSTER_MODE_RECLUSTER){
				reCluster(r1, r2, rt1, rt2);
			}else{
				throw new RuntimeException("Unknown mode "+clusterMode);
			}
		}
		
		/**
		 * Adds read pair to appropriate clusters during initial cluster creation.
		 * If both reads belong to the same cluster, adds them once. If they belong
		 * to different clusters, adds them to both clusters.
		 *
		 * @param r1 First read of the pair
		 * @param r2 Second read of the pair (may be null for single-end)
		 * @param rt1 ReadTag metadata for first read
		 * @param rt2 ReadTag metadata for second read (may be null)
		 */
		private void addToCluster(Read r1, Read r2, ReadTag rt1, ReadTag rt2){
			final int cn1=rt1.cluster0;
			final int cn2=rt2==null ? cn1 : rt2.cluster0;
			if(cn1==cn2){
				Cluster c1=fetchCluster(cn1);
				c1.add(r1);
				c1.add(r2);
			}else{
				Cluster c1=fetchCluster(cn1);
				Cluster c2=fetchCluster(cn2);
				c1.add(r1);
				c1.add(r2);
				c2.add(r1);
				c2.add(r2);
			}
		}
		
		/**
		 * Reassigns read pair to optimal clusters based on k-mer similarity scores.
		 * Evaluates all existing clusters and selects best matches using configured
		 * ambiguity resolution strategy when reads have conflicting assignments.
		 *
		 * @param r1 First read of the pair
		 * @param r2 Second read of the pair (may be null for single-end)
		 * @param rt1 ReadTag metadata for first read
		 * @param rt2 ReadTag metadata for second read (may be null)
		 */
		private void reCluster(Read r1, Read r2, ReadTag rt1, ReadTag rt2){

			assert(false) : "TODO";

			Cluster bestCluster1=null;
			Cluster bestCluster2=null;

			float bestScore1=-999999999, bestScore1_2=-999999999;
			float bestScore2=-999999999, bestScore2_1=-999999999;

			for(Cluster c : clusterList){

				float score1=c.score(r1);
				float score2=c.score(r2);

				if(bestCluster1==null || score1>bestScore1){
					bestCluster1=c;
					bestScore1=score1;
					bestScore1_2=score2;
				}
				if(bestCluster2==null || score2>bestScore2){
					bestCluster2=c;
					bestScore2=score2;
					bestScore2_1=score1;
				}
			}

			if(r2==null){
				rt1.cluster1=bestCluster1.id;
			}else if(bestCluster1==bestCluster2){
				rt1.cluster1=rt2.cluster1=bestCluster1.id;
			}else{
				assert(r1!=null && r2!=null && bestCluster1!=bestCluster2);

				float a=bestScore1+bestScore1_2;
				float b=bestScore2+bestScore2_1;

				if(ambigModeT==AMBIG_MODE_BEST){
					if(a>=b){
						rt1.cluster1=rt2.cluster1=bestCluster1.id;
					}else{
						rt1.cluster1=rt2.cluster1=bestCluster2.id;
					}
				}else if(ambigModeT==AMBIG_MODE_BOTH){
					assert(false) : "TODO";
				}else if(ambigModeT==AMBIG_MODE_TOSS){
					rt1.cluster1=rt2.cluster1=-1;
				}else if(ambigModeT==AMBIG_MODE_RAND){
					if(a<0 || b<0){
						float c=0-(Tools.min(a, b))*1.5f;
						a=a+c;
						b=a+c;
					}
					float coin=randy.nextFloat()*(a+b);
					if(coin<=a){
						rt1.cluster1=rt2.cluster1=bestCluster1.id;
					}else{
						rt1.cluster1=rt2.cluster1=bestCluster2.id;
					}
				}
			}

		}

		/** Unique identifier for this worker thread */
		final int id;
		/** Clustering operation mode for this thread */
		final int clusterMode;
		/** Thread-local ambiguity resolution strategy */
		final int ambigModeT;
		/** Input stream for reading sequences in this thread */
		final ConcurrentReadInputStream cris;
		
		/** Random number generator for probabilistic cluster assignment */
		Random randy;
		
		/** Thread-local count of reads processed */
		long readsInT;
		/** Thread-local count of bases processed */
		long basesInT;
		
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/** Tracks whether any errors occurred during processing */
	public boolean errorState=false;
	
	/** Thread-safe list of all clusters created during processing */
	final ArrayList<Cluster> clusterList=new ArrayList<Cluster>(256);
	
	/** 'big' kmer */
	public final int k1;
	/** 'small' kmer */
	public final int k2;

	/**
	 * Array length for k1-based k-mer storage, calculated as maxCanonicalKmer(k1)+1
	 */
	public final int arraylen1;
	/**
	 * Array length for k2-based k-mer storage, calculated as maxCanonicalKmer(k2)+1
	 */
	public final int arraylen2;
	
	/** Primary input file path for reads */
	private String in1=null;
	/** Secondary input file path for paired reads */
	private String in2=null;

	/** Primary output file path for processed reads */
	private String out1=null;
	/** Secondary output file path for paired reads */
	private String out2=null;
	
	private String extin=null;
	private String extout=null;
	
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files instead of overwriting */
	private boolean append=false;
	
	/** Maximum number of reads to process (-1 for unlimited) */
	private long maxReads=-1;
	
	/** Strategy for handling ambiguous cluster assignments */
	private int ambigMode=AMBIG_MODE_RAND;
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;
	
	private final FileFormat ffout1;
	private final FileFormat ffout2;
	
	/** Output stream for status messages and logging */
	private PrintStream outstream=System.err;
	
	/** Number of worker threads for concurrent processing */
	private int THREADS=Shared.threads();
	
	/*--------------------------------------------------------------*/

	/** Enables verbose logging output for debugging */
	public static boolean verbose=false;

	/** Clustering mode constant for initial cluster creation */
	public static final int CLUSTER_MODE_CREATE=0;
	/** Clustering mode constant for reassigning reads to clusters */
	public static final int CLUSTER_MODE_RECLUSTER=1;
	/** Clustering mode constant for refining existing cluster assignments */
	public static final int CLUSTER_MODE_REFINE=2;
	
	/** Ambiguity resolution mode: assign to best-scoring cluster */
	public static final int AMBIG_MODE_BEST=0;
	/** Ambiguity resolution mode: assign to both competing clusters */
	public static final int AMBIG_MODE_BOTH=1;
	/** Ambiguity resolution mode: discard ambiguous reads */
	public static final int AMBIG_MODE_TOSS=2;
	/** Ambiguity resolution mode: randomly assign based on score weights */
	public static final int AMBIG_MODE_RAND=3;
	
}
