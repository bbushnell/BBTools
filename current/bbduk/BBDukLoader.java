package bbduk;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.atomic.AtomicLongArray;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * Index and loader for BBDuk
 * @author Brian Bushnell
 * @date November 19, 2025
 *
 */
public class BBDukLoader {
	
	/**
	 * Constructor.
	 * @param p Parser with command line arguments
	 */
	public BBDukLoader(BBDukParser p){
		
		silent=p.silent;
		json=p.json;
		ALLOW_LOCAL_ARRAYS=p.ALLOW_LOCAL_ARRAYS;
		tossJunk=p.tossJunk;
		dump=p.dump;
		useShortKmers=p.useShortKmers;
		hammingDistance=p.hammingDistance;
		editDistance=p.editDistance;
		maxSkip=p.maxSkip;
		minSkip=p.minSkip;
		rcomp=p.rcomp;
		k=p.k;
		k2=p.k2;
		mink=p.mink;
		ktrimLeft=p.ktrimLeft;
		ktrimRight=p.ktrimRight;
		ktrimN=p.ktrimN;
		ksplit=p.ksplit;
		
		useRefNames=p.useRefNames;
		amino=p.amino;
		bitsPerBase=p.bitsPerBase;
		shift2=p.shift2;
		mask=p.mask;
		kmask=p.kmask;
		symbolToNumber=p.symbolToNumber;
		symbolToNumber0=p.symbolToNumber0;
		symbolToComplementNumber0=p.symbolToComplementNumber0;
		trimByOverlap=p.trimByOverlap;
		
		outstream=BBDukParser.outstream;
		overwrite=BBDukParser.overwrite;
		DISPLAY_PROGRESS=BBDukParser.DISPLAY_PROGRESS;
		THREADS=BBDukParser.workers;
		REPLICATE_AMBIGUOUS=BBDukParser.REPLICATE_AMBIGUOUS;
		
		parser=p;
		index=(p.WAYS==BBDukIndexMod.WAYS ? new BBDukIndexMod(p) : 
			p.indexmask2 ? new BBDukIndexMask2(p) : new BBDukIndexMask(p));
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public synchronized BBDukIndex loadIndex(String in1_for_header){
		index.loadAssorted(in1_for_header);
		
		/* Start overall timer */
		Timer t=new Timer();
		
		loadKmers(t.time1);
		
		/* Stop timer and calculate speed statistics */
		t.stop();
		assert(loaded);
		assert(index.loaded());
		return index;
	}
	
	
	/**
	 * Core processing method that loads reference kmers and processes input reads.
	 * Fills kmer tables from reference sequences.
	 * @param startTime Start time in nanoseconds for timing calculations
	 */
	private synchronized void loadKmers(long startTime){
		assert(!loaded);
		assert(!index.loaded());
		
		/* Start phase timer */
		Timer t=new Timer();
		
		/* Fill tables with reference kmers */
		if((index.ref!=null && index.ref.length>0) || (index.literal!=null && index.literal.length>0)){
			final boolean oldTI=FASTQ.TEST_INTERLEAVED; //TODO: This needs to be changed to a non-static field, or somehow 'read mode' and 'ref mode' need to be distinguished.
			final boolean oldFI=FASTQ.FORCE_INTERLEAVED;
			final boolean oldPC=FASTQ.PARSE_CUSTOM;
			final boolean oldSplit=FastaReadInputStream.SPLIT_READS;
			final int oldML=FastaReadInputStream.MIN_READ_LEN;
			
			FASTQ.TEST_INTERLEAVED=false;
			FASTQ.FORCE_INTERLEAVED=false;
			FASTQ.PARSE_CUSTOM=false;
			FastaReadInputStream.SPLIT_READS=false;
			FastaReadInputStream.MIN_READ_LEN=1;
			
			
			storedKmers=spawnLoadThreads(index.ways());
			if(storedKmers<1 && index.altRefNames.size()>0){
				outstream.println("Ref had no kmers; using alt ref.");
				index.ref=index.altref;
				index.refNames=index.altRefNames;
				index.refScafCounts=new int[index.refNames.size()];
				index.scaffoldNames.clear();
				index.scaffoldLengths.clear();
				storedKmers=spawnLoadThreads(index.ways());
			}
			
			FASTQ.TEST_INTERLEAVED=oldTI;
			FASTQ.FORCE_INTERLEAVED=oldFI;
			FASTQ.PARSE_CUSTOM=oldPC;
			FastaReadInputStream.SPLIT_READS=oldSplit;
			FastaReadInputStream.MIN_READ_LEN=oldML;
			
			if(useRefNames){index.toRefNames();}
			t.stop();
		}
		

		index.refReads=refReads;
		index.refBases=refBases;
		index.refKmers=refKmers;
		index.storedKmers=storedKmers;

		loaded=true;
		index.setKmersLoaded();
		
		{
			long ram=freeMemory();
			parser.ALLOW_LOCAL_ARRAYS=ALLOW_LOCAL_ARRAYS=(index.scaffoldNames!=null && Tools.max(THREADS, 1)*3*8*index.scaffoldNames.size()<ram*5);
		}
		
		/* Dump kmers to text */
		if(dump!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(dump, overwrite, false, true);
			bsw.start();
			index.dump(bsw, 0, Integer.MAX_VALUE);
			bsw.poisonAndWait();
		}
		
		if(storedKmers<1 && (ktrimRight || ktrimLeft || ktrimN || ksplit)){
			outstream.println("******  WARNING! A KMER OPERATION WAS CHOSEN BUT NO KMERS WERE LOADED.  ******");
			if(index.ref==null && index.literal==null){
				outstream.println("******  YOU NEED TO SPECIFY A REFERENCE FILE OR LITERAL SEQUENCE.       ******\n");
			}else{
				outstream.println("******  PLEASE ENSURE K IS LESS THAN OR EQUAL TO REF SEQUENCE LENGTHS.  ******\n");
			}
			if(ktrimRight && trimByOverlap){
				parser.ktrimRight=ktrimRight=false;
			}else{
				assert(false) : "You can bypass this assertion with the -da flag.";
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	

	/**
	 * Fills tables with kmers from references, using multiple LoadThread.
	 * @return Number of kmers stored.
	 */
	private long spawnLoadThreads(int threads){
		Timer t=new Timer();
		if((index.ref==null || index.ref.length<1) && (index.literal==null || index.literal.length<1)){return 0;}
		long added=0;
		
		/* Create load threads */
		LoadThread[] loaders=new LoadThread[threads];
		for(int i=0; i<loaders.length; i++){
			loaders[i]=new LoadThread(i);
			loaders[i].start();
		}
		
		/* For each reference file... */
		int refNum=0;
		if(index.ref!=null){
			for(String refname : index.ref){

				/* Start an input stream */
				FileFormat ff=FileFormat.testInput(refname, FileFormat.FASTA, null, false, true);
				ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1L, false, ff, null, null, null, Shared.USE_MPI, true);
				cris.start(); //4567
				ListNum<Read> ln=cris.nextList();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);
				
				/* Iterate through read lists from the input stream */
				while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
					{
						/* Assign a unique ID number to each scaffold */
						ArrayList<Read> reads2=new ArrayList<Read>(reads);
						for(Read r1 : reads2){
							final Read r2=r1.mate;
							final Integer id=index.scaffoldNames.size();
							index.refScafCounts[refNum]++;
							index.scaffoldNames.add(r1.id==null ? id.toString() : r1.id);
							int len=r1.length();
							r1.obj=id;
							if(r2!=null){
								r2.obj=id;
								len+=r2.length();
							}
							index.scaffoldLengths.add(len);
						}
						
						if(REPLICATE_AMBIGUOUS){
							reads2=Tools.replicateAmbiguous(reads2, Tools.min(k, mink));
						}

						/* Send a pointer to the read list to each LoadThread */
						for(LoadThread lt : loaders){
							boolean b=true;
							while(b){
								try {
									lt.queue.put(reads2);
									b=false;
								} catch (InterruptedException e) {
									//TODO:  This will hang due to still-running threads.
									throw new RuntimeException(e);
								}
							}
						}
					}

					/* Dispose of the old list and fetch a new one */
					cris.returnList(ln);
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				/* Cleanup */
				cris.returnList(ln);
				errorState|=ReadWrite.closeStream(cris);
				refNum++;
			}
		}

		/* If there are literal sequences to use as references */
		if(index.literal!=null){
			ArrayList<Read> list=new ArrayList<Read>(index.literal.length);
			if(verbose){outstream.println("Adding literals "+Arrays.toString(index.literal));}

			/* Assign a unique ID number to each literal sequence */
			for(int i=0; i<index.literal.length; i++){
				final Integer id=index.scaffoldNames.size();
				final Read r=new Read(index.literal[i].getBytes(), null, id);
				index.refScafCounts[refNum]++;
				index.scaffoldNames.add(id.toString());
				index.scaffoldLengths.add(r.length());
				r.obj=id;
				list.add(r);
			}
			
			if(REPLICATE_AMBIGUOUS){
				list=Tools.replicateAmbiguous(list, Tools.min(k, mink));
			}

			/* Send a pointer to the read list to each LoadThread */
			for(LoadThread lt : loaders){
				boolean b=true;
				while(b){
					try {
						lt.queue.put(list);
						b=false;
					} catch (InterruptedException e) {
						//TODO:  This will hang due to still-running threads.
						throw new RuntimeException(e);
					}
				}
			}
		}
		
		/* Signal loaders to terminate */
		for(LoadThread lt : loaders){
			boolean b=true;
			while(b){
				try {
					lt.queue.put(POISON);
					b=false;
				} catch (InterruptedException e) {
					//TODO:  This will hang due to still-running threads.
					throw new RuntimeException(e);
				}
			}
		}
		
		/* Wait for loaders to die, and gather statistics */
		boolean success=true;
		for(LoadThread lt : loaders){
			while(lt.getState()!=Thread.State.TERMINATED){
				try {
					lt.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			added+=lt.addedT;
			refKmers+=lt.refKmersT;
			refBases+=lt.refBasesT;
			refReads+=lt.refReadsT;
			success&=lt.success;
		}
		if(!success){KillSwitch.kill("Failed loading ref kmers; aborting.");}
		
		//Correct statistics for number of threads, since each thread processes all reference data
		refKmers/=threads;
		refBases/=threads;
		refReads/=threads;
		
		index.scaffoldReadCounts=new AtomicLongArray(index.scaffoldNames.size());
		index.scaffoldBaseCounts=new AtomicLongArray(index.scaffoldNames.size());

		t.stop();
		if(DISPLAY_PROGRESS && !json){
			outstream.println("Added "+added+" kmers; time: \t"+t);
			Shared.printMemory();
			outstream.println();
		}
		
		if(verbose){
			ByteStreamWriter tsw=new ByteStreamWriter("stdout", false, false, false, FileFormat.TEXT);
			tsw.start();
			index.dump(tsw, 0, Integer.MAX_VALUE);
			tsw.poisonAndWait();
		}
		
		return added;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads kmers into a table.  Each thread handles all kmers X such that X%WAYS==tnum.
	 */
	private class LoadThread extends Thread{
		
		public LoadThread(final int tnum_){
			tnum=tnum_;
		}
		
		/**
		 * Get the next list of reads (or scaffolds) from the queue.
		 * @return List of reads
		 */
		private ArrayList<Read> fetch(){
			ArrayList<Read> list=null;
			while(list==null){
				try {
					list=queue.take();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			return list;
		}
		
		@Override
		public void run(){
			ArrayList<Read> reads=fetch();
			while(reads!=POISON){
				for(Read r1 : reads){
					assert(r1.pairnum()==0);
					final Read r2=r1.mate;

					final int rblen=(r1==null ? 0 : r1.length());
					final int rblen2=r1.mateLength();
					
					addedT+=addToMap(r1, rblen>20000000 ? k : rblen>5000000 ? 11 : rblen>500000 ? 2 : 0);
					if(r2!=null){
						addedT+=addToMap(r2, rblen2>20000000 ? k : rblen2>5000000 ? 11 : rblen2>500000 ? 2 : 0);
					}
				}
				reads=fetch();
			}
			
			index.rebalance(tnum);
			success=true;
		}

		/**
		 * Store the read's kmers in a table.
		 * @param r The current read to process
		 * @param skip Number of bases to skip between kmers
		 * @return Number of kmers stored
		 */
		private long addToMap(Read r, int skip){
			skip=Tools.max(minSkip, Tools.min(maxSkip, skip));
			final byte[] bases=r.bases;
			long kmer=0;
			long rkmer=0;
			long added=0;
			int len=0;
			
			if(bases!=null){
				refReadsT++;
				refBasesT+=bases.length;
			}
			if(bases==null || bases.length<k){return 0;}
			
			final int id=(Integer)r.obj;
			
			if(skip>1){ //Process while skipping some kmers
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					long x=symbolToNumber0[b];
					long x2=symbolToComplementNumber0[b];
					kmer=((kmer<<bitsPerBase)|x)&mask;
					rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
					if(isFullyDefined(b)){len++;}else{len=0; rkmer=0;}
					if(verbose){outstream.println("Scanning1 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=k){
						refKmersT++;
						if(len%skip==0){
							final long extraBase=(i>=bases.length-1 ? -1 : symbolToNumber[bases[i+1]]);
							added+=index.addToMap(kmer, rkmer, k, extraBase, id, kmask, hammingDistance, editDistance, tnum);
							if(useShortKmers){
								if(i==k2){added+=index.addToMapRightShift(kmer, rkmer, id, tnum);}
								if(i==bases.length-1){added+=index.addToMapLeftShift(kmer, rkmer, extraBase, id, tnum);}
							}
						}
					}
				}
			}else{ //Process all kmers
				for(int i=0; i<bases.length; i++){
					final byte b=bases[i];
					final long x=symbolToNumber0[b];
					final long x2=symbolToComplementNumber0[b];
					kmer=((kmer<<bitsPerBase)|x)&mask;
					rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
					if(isFullyDefined(b)){len++;}else{len=0; rkmer=0;}
					if(verbose){
						if(verbose){
							String fwd=new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k));
							String rev=AminoAcid.reverseComplementBases(fwd);
							String fwd2=kmerToString(kmer, Tools.min(len, k));
							outstream.println("fwd="+fwd+", fwd2="+fwd2+", rev="+rev+", kmer="+kmer+", rkmer="+rkmer);
							outstream.println("b="+(char)b+", x="+x+", x2="+x2+", bitsPerBase="+bitsPerBase+", shift2="+shift2);
							if(!amino){
								assert(AminoAcid.stringToKmer(fwd)==kmer) : fwd+", "+AminoAcid.stringToKmer(fwd)+", "+kmer+", "+len;
								if(len>=k){
									assert(rcomp(kmer, Tools.min(len, k))==rkmer);
									assert(rcomp(rkmer, Tools.min(len, k))==kmer);
									assert(AminoAcid.kmerToString(kmer, Tools.min(len, k)).equals(fwd));
									assert(AminoAcid.kmerToString(rkmer, Tools.min(len, k)).equals(rev)) : AminoAcid.kmerToString(rkmer, Tools.min(len, k))+" != "+rev+" (rkmer)";
								}
								assert(fwd.equalsIgnoreCase(fwd2)) : fwd+", "+fwd2; //may be unsafe
							}
							outstream.println("Scanning6 i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+", bases="+fwd+", rbases="+rev);
						}
					}
					if(len>=k){
						refKmersT++;
						final long extraBase=(i>=bases.length-1 ? -1 : symbolToNumber[bases[i+1]]);
						final long atm=index.addToMap(kmer, rkmer, k, extraBase, id, kmask, hammingDistance, editDistance, tnum);
						added+=atm;
						if(useShortKmers){
							if(i==k2){added+=index.addToMapRightShift(kmer, rkmer, id, tnum);}
							if(i==bases.length-1){added+=index.addToMapLeftShift(kmer, rkmer, extraBase, id, tnum);}
						}
					}
				}
			}
			return added;
		}
		
		/*--------------------------------------------------------------*/
		
		/** Number of kmers stored by this thread */
		public long addedT=0;
		/** Number of items encountered by this thread */
		public long refKmersT=0, refReadsT=0, refBasesT=0;
		/** Thread number; used to determine which kmers to store */
		public final int tnum;
		/** Buffer of input read lists */
		public final ArrayBlockingQueue<ArrayList<Read>> queue=new ArrayBlockingQueue<ArrayList<Read>>(32);
		
		/** Completed successfully */
		boolean success=false;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Current available memory */
	private static final long freeMemory(){
		Runtime rt=Runtime.getRuntime();
		return rt.freeMemory();
	}
	
	/**
	 * Computes reverse complement of kmer.
	 * @param kmer Input kmer
	 * @param len Kmer length
	 * @return Reverse complement kmer
	 */
	final long rcomp(long kmer, int len){
		return amino ? kmer : AminoAcid.reverseComplementBinaryFast(kmer, len);
	}
	
	/** For verbose / debugging output */
	final String kmerToString(long kmer, int k){
		return amino ? AminoAcid.kmerToStringAA(kmer, k) : AminoAcid.kmerToString(kmer, k);
	}
	
	/** Returns true if the symbol is not degenerate (e.g., 'N') for the alphabet in use. */
	final boolean isFullyDefined(byte symbol){
		return symbol>=0 && symbolToNumber[symbol]>=0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final BBDukIndex index;
	
	private boolean loaded=false;
	
	/** Has this class encountered errors while processing? */
	public boolean errorState;
	
	/** Number of reference reads processed for kmer loading */
	long refReads;
	/** Number of reference bases processed for kmer loading */
	long refBases;
	/** Number of reference kmers encountered during loading */
	long refKmers;
	/** Number of unique kmers actually stored in hash tables */
	long storedKmers;
	
	/** Set to false to force threads to share atomic counter arrays. */
	boolean ALLOW_LOCAL_ARRAYS;
	
	/*--------------------------------------------------------------*/
	/*----------------          Immutable           ----------------*/
	/*--------------------------------------------------------------*/

	private final BBDukParser parser;
	private final boolean silent;
	private final boolean json;
	
	final boolean tossJunk;
	
	/** Dump kmers here. */
	private final String dump;
	
	/** Attempt to match kmers shorter than normal k on read ends when doing kTrimming. */
	private final boolean useShortKmers;
	
	/** Store reference kmers with up to this many substitutions */
	private final int hammingDistance;
	/** Store reference kmers with up to this many edits (including indels) */
	private final int editDistance;
	/** Never skip more than this many consecutive kmers when hashing reference. */
	private final int maxSkip;
	/** Always skip at least this many consecutive kmers when hashing reference.
	 * 1 means every kmer is used, 2 means every other, etc. */
	private final int minSkip;
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Look for reverse-complements as well as forward kmers.  Default: true */
	private final boolean rcomp;
	
	/** Normal kmer length */
	private final int k;
	/** k-1; used in some expressions */
	private final int k2;
	/** Shortest kmer to use for trimming */
	private final int mink;
	/** Trim matching kmers and all bases to the left */
	private final boolean ktrimLeft;
	/** Trim matching kmers and all bases to the right */
	boolean ktrimRight;
	/** Don't trim, but replace matching kmers with a symbol (default N) */
	private final boolean ktrimN;
	/** Split into two reads around the kmer */
	private final boolean ksplit;
	/** Use names of reference files instead of scaffolds.
	 * Default: false. */
	private final boolean useRefNames;
	
	/*--------------------------------------------------------------*/
	/*-----------        Symbol-Specific Constants        ----------*/
	/*--------------------------------------------------------------*/

	/** True for amino acid data, false for nucleotide data */
	private final boolean amino;
	private final int bitsPerBase;
	private final int shift2;
	private final long mask;
	private final long kmask;
	
	/** Symbol code; -1 for undefined */
	private final byte[] symbolToNumber;
	/** Symbol code; 0 for undefined */
	private final byte[] symbolToNumber0;
	/** Complementary symbol code; 0 for undefined */
	private final byte[] symbolToComplementNumber0;
	
	/*--------------------------------------------------------------*/
	/*----------------         BBMerge Flags        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Trim implied adapters based on overlap, for reads with insert size shorter than read length */
	private final boolean trimByOverlap;
	
	/*--------------------------------------------------------------*/
	
	/** Print messages to this stream */
	private final PrintStream outstream;
	/** Permission to overwrite existing files */
	private final boolean overwrite;
	/** Display progress messages such as memory usage */
	private final boolean DISPLAY_PROGRESS;
	/** Number of ProcessThreads */
	private final int THREADS;
	/** Make unambiguous copies of ref sequences with ambiguous bases */
	private final boolean REPLICATE_AMBIGUOUS;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Verbose messages */
	public static final boolean verbose=false;
	/** Indicates end of input stream */
	private static final ArrayList<Read> POISON=new ArrayList<Read>(0);
	
}
