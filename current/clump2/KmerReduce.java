package clump2;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import assemble.AbstractRemoveThread;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
import kmer.KmerTableSet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Emits one representative k-base pivot sequence per input root read. Linked mates
 * can participate in overlap correction but are not emitted or counted separately.
 * Output reads have the constant name "1", the original numeric ID, and no quality
 * array; the requested output format controls their serialization.
 * Instances and shared parser/comparator settings require serialized configuration.
 * @author Brian Bushnell
 * @date Nov 10, 2015
 */
public class KmerReduce{

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parses and processes files, restoring pigz/unpigz flags on normal return.
	 * Other shared parser settings and exceptional exits are not restored. */
	public static void main(String[] args){
		final boolean pigz=ReadWrite.USE_PIGZ, unpigz=ReadWrite.USE_UNPIGZ;
		Timer t=new Timer();
		KmerReduce x=new KmerReduce(args);
		x.process(t);
		ReadWrite.USE_PIGZ=pigz;
		ReadWrite.USE_UNPIGZ=unpigz;
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Extracts valid k-mers from a read file by processing reads through KmerReduce.
	 * This unfinished entry point writes a temporary pivot file, then deliberately
	 * asserts before building the table. Do not use with assertions enabled.
	 * The assertion is retained as a historical unfinished-work guard.
	 * @param fname0 Input filename of reads
	 * @param k K-mer length
	 * @param cutoff Minimum count to retain in the final table
	 * @return Set of pivot k-mers meeting the count threshold
	 */
	public static KmerTableSet getValidKmersFromReads(final String fname0, int k, int cutoff){
		final String fname=fname0+"_"+(shared.Shared.random().nextLong()>>>1)+".fa.gz";
		assert(!new File(fname).exists());

		ArrayList<String> arglist=new ArrayList<String>();
		arglist.add("in="+fname0);
		arglist.add("out="+fname);
		arglist.add("k="+k);
		String[] args=arglist.toArray(new String[0]);

		main(args);

		//NOTE [clump/KmerReduce#002]: this method is DEAD (no callers anywhere in the tree) and INTENTIONALLY
		//disabled - the unconditional assert(false) below (dumping params) is a "not finished / do not use" guard,
		//like KmerSort3's no-arg main. Under -ea (always on) it always throws here, so getValidKmers below is never
		//reached via this path. Not a live bug; flagged as dead/incomplete so a future caller doesn't trip over it.
		assert(false) : fname+", "+k+", "+cutoff;
		KmerTableSet set=getValidKmers(fname, k, cutoff);
		File f=new File(fname);
		if(f.exists()){f.delete();}
		
		return set;
	}
	
	/**
	 * Builds a k-mer table from a file containing pivot k-mers.
	 * Applies optional prefiltering and count-based filtering to retain only k-mers meeting the cutoff.
	 * @param fname Input filename of pivot k-mers
	 * @param k K-mer length
	 * @param cutoff Minimum count to retain in the final table
	 * @return Set of pivot k-mers meeting the count threshold
	 */
	public static KmerTableSet getValidKmers(final String fname, int k, int cutoff){
		ArrayList<String> arglist=new ArrayList<String>();
		arglist.add("in="+fname);
		arglist.add("k="+k);
		if(cutoff>1 && prefilter){
			arglist.add("prefilter="+(cutoff-1));
		}
		
		String[] args=arglist.toArray(new String[0]);
		KmerTableSet set=new KmerTableSet(args, 12);
		
		Timer t=new Timer();
		
		set.process(t);
//		errorState|=set.errorState;
		assert(!set.errorState);
		t.stop();
		
		set.prefilterArray=null;
		AbstractRemoveThread.process(Shared.threads(), cutoff, Integer.MAX_VALUE, set, true);
		
		return set;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Parses command-line arguments and initializes processing parameters.
	 * Sets up file descriptors and mutates shared parser/compression settings.
	 * Streams are created and started only by process().
	 * @param args Command line arguments
	 */
	public KmerReduce(String[] args){
		
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
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
				assert(k>0 && k<32);
			}else if(a.equals("comparisons") || a.equals("c")){
				//do nothing
			}else if(a.equals("ecco")){
				ecco=Parse.parseBoolean(b);
			}else if(a.equals("rename") || a.equals("addname")){
				//do nothing
			}else if(a.equals("rcomp") || a.equals("reversecomplement")){
				//do nothing
			}else if(a.equals("condense") || a.equals("consensus")){
				//do nothing
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
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;
			in2=parser.in2;

			out1=parser.out1;
			
			extin=parser.extin;
			extout=parser.extout;
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

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		//Fixed [clump/KmerReduce#001]: reject either input as output before starting the reader.
		ClumpTools.validateFileNames(new FileFormat[] {ffin1, ffin2, ffout1});
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates read streams and processes all data.
	 * Sets up Streamer/Writer I/O, runs HashThread workers, and reports timing and error state.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		final Streamer cris;
		{
			cris=StreamerFactory.getReadInputStream(maxReads, true, ffin1, ffin2, null, null, -1);
			cris.start();
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
		
		final Writer ros;
		if(out1!=null){
			final int buff=Tools.max(4, Shared.threads());

			ros=WriterFactory.getStream(ffout1, null, buff, null, false, -1);
			ros.start();
		}else{ros=null;}
		
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read stream
		processInner(cris, ros);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		errorState|=ReadStats.writeAll();
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		t.stop();
		
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		if(errorState){
			Clumpify.sharedErrorState=true;
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * Manages multi-threaded processing using HashThread workers.
	 * Creates and joins workers, then aggregates root-read and root-base counts.
	 * Does not start or close the supplied streams; process() owns that lifecycle.
	 * @param cris Already-started input stream
	 * @param ros Already-started unordered output stream, or null to discard pivots
	 */
	public void processInner(final Streamer cris, final Writer ros){
		if(verbose){outstream.println("Making comparator.");}
		KmerComparator kc=new KmerComparator(k, false, false);
		
		if(verbose){outstream.println("Making hash threads.");}
		final int threads=Shared.threads();
		ArrayList<HashThread> alht=new ArrayList<HashThread>(threads);
		for(int i=0; i<threads; i++){alht.add(new HashThread(cris, ros, kc));}
		
		if(verbose){outstream.println("Starting threads.");}
		for(HashThread ht : alht){ht.start();}
		
		if(verbose){outstream.println("Waiting for threads.");}
		/* Wait for threads to die */
		for(HashThread ht : alht){
			
			/* Wait for a thread to die */
			while(ht.getState()!=Thread.State.TERMINATED){
				try{
					ht.join();
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
			readsProcessed+=ht.readsProcessedT;
			basesProcessed+=ht.basesProcessedT;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Consumes borrowed batches, emits new pivot reads, and returns every consumed list. */
	private class HashThread extends Thread{

		HashThread(Streamer cris_, Writer ros_, KmerComparator kc_){
			cris=cris_;
			ros=ros_;
			kc=kc_;
		}

		/** Processes root reads only; optional overlap correction can mutate both mates. */
		@Override
		public void run(){
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			while(ln!=null && reads!=null){//Streamer permits empty nonterminal batches.
				ArrayList<Read> out=new ArrayList<Read>(reads.size());
				for(Read r : reads){
					if(ecco && r.mate!=null){
						if(r.mate!=null){BBMerge.findOverlapStrict(r, r.mate, true);}
					}
					final long kmer=kc.hash(r, null, 0, false);
					readsProcessedT++;
					basesProcessedT+=r.length();
					if(kmer>=0){
						Read temp=new Read(toBytes(kmer), null, header, r.numericID);
						out.add(temp);
					}
				}
				if(ros!=null){ros.add(out, ln.id);}
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		final Streamer cris;
		final Writer ros;
		final KmerComparator kc;
		
		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
		
		private static final String header="1";
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Allocates and decodes exactly k bases from the low 2*k bits of a packed pivot. */
	public byte[] toBytes(final long kmer){
		byte[] dest=KillSwitch.allocByte1D(k);
		fill(kmer, dest, 0);
		return dest;
	}
	
	//Decodes a 2-bit-packed kmer back to ACGT bases, high-order 2 bits first (i=k-1 down to 0) so dest[0] is the
	//leftmost base - the inverse of the kmer-building shift in KmerComparator. Used to emit each read's pivot kmer
	//as a FASTA sequence (the whole point of KmerReduce: reduce a read to its representative kmer).
	/** Writes k A/C/G/T bytes at pos, most-significant selected base first.
	 * The caller must supply room for all k bytes; no terminator is written. */
	public void fill(final long kmer, final byte[] dest, int pos){
		for(int i=k-1; i>=0; i--, pos++){
			int x=(int)((kmer>>(2*i))&3);
			dest[pos]=AminoAcid.numberToBase[x];
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Pivot length, limited by the packed 2-bit representation to 1 through 31. */
	private int k=31;
	/** Whether getValidKmers uses a count prefilter when cutoff exceeds one. */
	static boolean prefilter=true;
	
	/*--------------------------------------------------------------*/
	/*----------------          I/O Fields          ----------------*/
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String in2=null;
	
	private String out1=null;
	
	private String extin=null;
	private String extout=null;
	
	/*--------------------------------------------------------------*/
	
	/** Processed root reads; linked mates are excluded. */
	protected long readsProcessed=0;
	/** Processed root-read bases; linked mates are excluded. */
	protected long basesProcessed=0;
	
	/** Input fragment limit, or a negative value for no configured limit. */
	private long maxReads=-1;
	/** Whether to perform strict overlap-based correction before choosing each pivot. */
	protected boolean ecco=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;

	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
