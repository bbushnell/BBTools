package align2;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/** Builds or loads flat k-mer blocks used by BBIndex.
 * Each block has four workers with disjoint key ranges. They count hits, share
 * one allocation, then fill their ranges without per-hit locking. A block is
 * published only after all workers complete successfully. Cached and generated
 * worker failures reach the caller after started workers have been drained.
 * Reference state and build/filter settings are global: callers must serialize
 * reconfiguration and keep chromosome arrays unchanged for the whole build.
 * Disk writes retain Block.write's asynchronous sites-file behavior; returning
 * from makeIndex is not a guarantee that all cache writes have completed.
 * @author Brian Bushnell, Collei
 * @date Dec 23, 2012
 *
 */
public class IndexMaker4{

	/**
	 * Creates k-mer indices for a range of chromosomes with multi-threaded processing.
	 * Initializes BlockMaker threads for each chromosome range and coordinates their execution.
	 * Successful blocks may remain in a supplied destination after a failed batch;
	 * a failed batch throws rather than returning that partially populated array.
	 *
	 * @param genome Genome build number for reference selection
	 * @param minChrom Minimum chromosome number to index
	 * @param maxChrom Maximum chromosome number to index
	 * @param k K-mer length for indexing
	 * @param CHROMBITS Number of bits allocated for chromosome encoding
	 * @param MAX_ALLOWED_CHROM_INDEX Maximum allowed chromosome index size
	 * @param CHROM_MASK_LOW Bitmask for low chromosome bits
	 * @param CHROM_MASK_HIGH Bitmask for high chromosome bits
	 * @param SITE_MASK Bitmask for genomic site extraction
	 * @param SHIFT_LENGTH Bit shift amount for encoding
	 * @param WRITE Whether to write indices to disk
	 * @param DISK_INVALID Whether to ignore existing disk indices
	 * @param index Pre-allocated index array or null
	 * @return Array of Block objects containing k-mer indices
	 */
	public static Block[] makeIndex(final int genome, int minChrom, int maxChrom, int k, int CHROMBITS,
			int MAX_ALLOWED_CHROM_INDEX, int CHROM_MASK_LOW, int CHROM_MASK_HIGH, int SITE_MASK, int SHIFT_LENGTH,
			boolean WRITE, boolean DISK_INVALID, Block[] index){
		Timer t=new Timer();

		MAX_CONCURRENT_BLOCKS=(Shared.LOW_MEMORY ? 1 : (Shared.WINDOWS ? (WRITE ? 1 : Tools.max(1, Shared.threads()/4)) : Tools.max(1, Shared.threads()/4)));

		minChrom=Tools.max(1, minChrom);
		if(genome>=0 && Data.GENOME_BUILD!=genome){
			Data.setGenome(genome);
			maxChrom=Tools.min(Data.numChroms, maxChrom);
		}

		assert(minChrom<=maxChrom) : "Index build requires a nonempty inclusive chromosome range: "+minChrom+".."+maxChrom;

		if(index==null){index=new Block[maxChrom+1];}

		ArrayList<BlockMaker> list=new ArrayList<BlockMaker>();

		Throwable failure=null;
		try{
			for(int i=1; i<=maxChrom;){
				if(i>=minChrom){
					final int a=minChrom(i, minChrom, CHROM_MASK_HIGH);
					final int b=maxChrom(i, minChrom, maxChrom, CHROM_MASK_LOW);
					assert(b>=i) : "Block upper bound must advance the chromosome loop: chromosome="+i+", end="+b;
					final BlockMaker idm=new BlockMaker(a, b, k, CHROMBITS, MAX_ALLOWED_CHROM_INDEX, CHROM_MASK_LOW, CHROM_MASK_HIGH, SITE_MASK, SHIFT_LENGTH, WRITE, DISK_INVALID, index);
					list.add(idm);
					incrementActiveBlocks(1);
					try{idm.start();}
					catch(RuntimeException|Error problem){incrementActiveBlocks(-1); throw problem;}
					i=b+1;
				}else{i++;}
			}
		}catch(RuntimeException|Error problem){failure=problem;}

		// Drain every started worker even when reservation or launch failed.
		boolean interrupted=false;
		for(BlockMaker cm : list){
			while(cm.isAlive()){
				try{cm.join();}
				catch(InterruptedException problem){
					interrupted=true;
					if(failure==null){failure=new RuntimeException("Interrupted while joining index block builder", problem);}
				}
			}
			if(cm.failure!=null){
				if(failure==null){failure=cm.failure;}
				else if(failure!=cm.failure){failure.addSuppressed(cm.failure);}
			}
		}
		if(interrupted){Thread.currentThread().interrupt();}
		if(failure!=null){rethrow(failure);}

		t.stop();

		return index;
	}

	/** Preserves the original unchecked failure after workers have been joined. */
	private static void rethrow(final Throwable problem){
		if(problem instanceof Error){throw (Error)problem;}
		if(problem instanceof RuntimeException){throw (RuntimeException)problem;}
		throw new RuntimeException(problem);
	}

	/**
	 * Creates a single Block index for a chromosome range.
	 * Note: Contains assertion failures suggesting this method is deprecated.
	 *
	 * @param minChrom Minimum chromosome number
	 * @param maxChrom Maximum chromosome number
	 * @param k K-mer length
	 * @param CHROMBITS Number of bits for chromosome encoding
	 * @param MAX_ALLOWED_CHROM_INDEX Maximum chromosome index size
	 * @param CHROM_MASK_LOW Low chromosome bitmask
	 * @param CHROM_MASK_HIGH High chromosome bitmask
	 * @param SITE_MASK Site extraction bitmask
	 * @param SHIFT_LENGTH Encoding shift amount
	 * @param WRITE Whether to write to disk
	 * @param DISK_INVALID Whether to ignore disk cache
	 * @param matrix Index matrix for storage
	 * @return Block containing k-mer index data
	 */
	public static Block makeBlock(int minChrom, int maxChrom, int k, int CHROMBITS, int MAX_ALLOWED_CHROM_INDEX,
			int CHROM_MASK_LOW, int CHROM_MASK_HIGH, int SITE_MASK, int SHIFT_LENGTH, boolean WRITE, boolean DISK_INVALID, Block[] matrix){
		//TODO [align2/IndexMaker4#002]: dead method - ZERO callers (grep finds only defs in IndexMaker4+IndexMaker5,
		//no invocations) and the two assert(false) guards make it unusable under -ea. Superseded single-block relic;
		//live path is makeIndex->BlockMaker. (Same dead method in IndexMaker5, the unfinished successor.)
		assert(false) : maxChrom+", "+MAX_ALLOWED_CHROM_INDEX;
		BlockMaker idm=new BlockMaker(minChrom, maxChrom, k, CHROMBITS, MAX_ALLOWED_CHROM_INDEX, CHROM_MASK_LOW, CHROM_MASK_HIGH, SITE_MASK, SHIFT_LENGTH, WRITE, DISK_INVALID, matrix);
		Block block=idm.makeArrays();

		assert(false) : maxChrom+", "+MAX_ALLOWED_CHROM_INDEX;

		if(verbose){
			for(int i=0; i<block.numStarts; i++){
				int[] array=block.getHitList(i);
				if(array==null){Data.sysout.println(i+": "+null);}
				else{Data.sysout.println(i+": "+Arrays.toString(array));}
			}
		}

		return block;
	}

	/**
	 * Thread-based worker for creating k-mer index blocks for chromosome ranges.
	 * Manages the complete lifecycle from reading cached indices to generating new ones.
	 * Coordinates multiple CountThread workers for parallel k-mer processing.
	 */
	private static class BlockMaker extends Thread{

		/**
		 * Constructs a BlockMaker for a specific chromosome range and configuration.
		 *
		 * @param minChrom_ Minimum chromosome to process
		 * @param maxChrom_ Maximum chromosome to process
		 * @param k K-mer length
		 * @param CHROMBITS_ Chromosome encoding bits
		 * @param MAX_ALLOWED_CHROM_INDEX_ Maximum chromosome index size
		 * @param CHROM_MASK_LOW_ Low chromosome bitmask
		 * @param CHROM_MASK_HIGH_ High chromosome bitmask
		 * @param SITE_MASK_ Site extraction bitmask
		 * @param SHIFT_LENGTH_ Encoding shift length
		 * @param WRITE_TO_DISK_ Whether to cache results on disk
		 * @param DISK_INVALID_ Whether to ignore existing disk cache
		 * @param matrix_ Index storage matrix
		 */
		public BlockMaker(int minChrom_, int maxChrom_, int k, int CHROMBITS_,
				int MAX_ALLOWED_CHROM_INDEX_, int CHROM_MASK_LOW_, int CHROM_MASK_HIGH_, int SITE_MASK_, int SHIFT_LENGTH_,
				boolean WRITE_TO_DISK_, boolean DISK_INVALID_, Block[] matrix_){

			KEYLEN=k;
			CHROMBITS=CHROMBITS_;
			KEYSPACE=1<<(2*KEYLEN);
			MAX_ALLOWED_CHROM_INDEX=MAX_ALLOWED_CHROM_INDEX_;
			WRITE_TO_DISK=WRITE_TO_DISK_;
			DISK_INVALID=DISK_INVALID_;

			CHROM_MASK_LOW=CHROM_MASK_LOW_;
			CHROM_MASK_HIGH=CHROM_MASK_HIGH_;
			SITE_MASK=SITE_MASK_;
			SHIFT_LENGTH=SHIFT_LENGTH_;

			minChrom=minChrom_;
			maxChrom=maxChrom_;
			matrix=matrix_;
		}

		@Override
		public void run(){
			try{makeArrays();}
			catch(RuntimeException|Error problem){failure=problem;}
			finally{incrementActiveBlocks(-1);}
		}

		/** Inspected by makeIndex after joining this worker. */
		private volatile Throwable failure;

		/** Wake both build barriers when any worker or its launcher fails. */
		private void failBuild(final Throwable problem, final int[] intercom){
			synchronized(intercom){
				if(failure==null){failure=problem;}
				intercom.notifyAll();
			}
		}

		/**
		 * Creates k-mer index arrays for the assigned chromosome range.
		 * First attempts to load from disk cache, then generates from scratch using CountThreads.
		 * Manages the two-phase process: counting k-mer occurrences, then filling index arrays.
		 * @return Block containing the completed k-mer index
		 */
		Block makeArrays(){

			if(!DISK_INVALID){
				String fname=fname(minChrom, maxChrom, KEYLEN, CHROMBITS);
				File f=new File(fname);

				if(f.exists() && new File(fname+"2.gz").exists()){
					Block x=Block.read(fname);
					if(matrix!=null){
						for(int i=baseChrom(minChrom); i<=maxChrom; i++){
							matrix[i]=x;
						}
					}
					return x;
				}else{
					synchronized(getClass()){
						Data.sysout.println("No index available; generating from reference genome: "+f.getAbsolutePath());
						if(WRITE_TO_DISK){
							String root=ReadWrite.parseRoot2(f.getAbsolutePath());
							File rf=new File(root);
							if(!rf.exists()){
								rf.mkdirs();
							}
						}
					}
				}
			}

			CountThread threads[]=new CountThread[4];
			int[] sizes=KillSwitch.allocInt1D(KEYSPACE+1);
			int[] intercom=KillSwitch.allocInt1D(4);
			Block[] indexHolder=new Block[1];

			try{
				for(int i=0; i<4; i++){
					threads[i]=new CountThread(i, sizes, intercom, indexHolder);
					threads[i].start();
				}
			}catch(RuntimeException|Error problem){failBuild(problem, intercom);}
			Data.sysout.println("Indexing threads started for block "+baseChrom(minChrom)+"-"+maxChrom);
			boolean interrupted=false;
			for(int i=0; i<threads.length; i++){
				if(threads[i]==null){continue;}
				while(threads[i].isAlive()){
					try{threads[i].join();}
					catch(InterruptedException problem){
						interrupted=true;
						failBuild(new RuntimeException("Interrupted while joining index count/fill workers", problem), intercom);
					}
				}
			}
			if(interrupted){Thread.currentThread().interrupt();}
			if(failure!=null){rethrow(failure);}
			Data.sysout.println("Indexing threads finished for block "+baseChrom(minChrom)+"-"+maxChrom);

			for(int i=sizes.length-2; i>=0; i--){
				sizes[i+1]=sizes[i];
			}
			sizes[0]=0;

			if(matrix!=null){
				for(int i=baseChrom(minChrom); i<=maxChrom; i++){
					matrix[i]=indexHolder[0];
				}
			}

			if(WRITE_TO_DISK){
				String fname=fname(minChrom, maxChrom, KEYLEN, CHROMBITS);
//				File f=new File(fname);
				indexHolder[0].write(fname, true);
			}

			return indexHolder[0];
		}

		/**
		 * Worker thread that processes k-mers for a specific base (A, C, G, T) range.
		 * Handles both counting phase (determining array sizes) and filling phase (populating indices).
			 * Uses the same periodic-key and optional modulo filters in both passes.
			 * Key ownership is the high two bits of the forward key, not chromosome.
		 */
		private class CountThread extends Thread{

			/**
			 * Constructs a CountThread for processing k-mers starting with a specific base.
			 *
			 * @param id_ Thread ID (0-3 for A, C, G, T)
			 * @param sizes_ Shared array for tracking k-mer counts
			 * @param intercom_ Inter-thread communication array
			 * @param indexHolder_ Shared holder for the final Block index
			 */
			public CountThread(int id_, int[] sizes_, int[] intercom_, Block[] indexHolder_){
				id=id_;
				sizes=sizes_;
				indexHolder=indexHolder_;
				intercom=intercom_;

				minIndex=(id<<(2*KEYLEN-2));
				maxIndex=(int)(((id+1L)<<(2*KEYLEN-2))-1);

				if(ALLOW_POLYMERS){
					banmask=-1; //The zero key (poly-A) is still excluded by the comparison.
				}else{
					//TODO [align2/IndexMaker4#001]: dead store - `b` (the base-id homopolymer key) is built here
					//but never read; banmask below uses the id-independent periodicity formula. Harmless vestige of
					//an earlier exact-homopolymer filter, superseded by the (key>>banshift)!=(key&banmask) clumpy filter.
					int b=0;
					for(int i=0; i<KEYLEN; i++){
						b<<=2;
						b=(b|id);
					}
					banmask=~((-1)<<((2*KEYLEN)-banshift));
				}
			}

			/** Thread identifier (0-3 for bases A, C, G, T) */
			private final int id;
			/** Shared array for counting k-mer occurrences */
			private final int[] sizes;
			/** Barrier monitor; slot1 counts finished counters,slot2 marks allocation. Slots0/3 are unused. */
			private final int[] intercom;
			/** Shared holder for the final Block index */
			private final Block[] indexHolder;
			/** Minimum k-mer index value handled by this thread */
			private final int minIndex;
			/** Maximum k-mer index value handled by this thread */
			private final int maxIndex;
			/** Bitmask for filtering homopolymer k-mers */
			private final int banmask;
			/** Bit shift amount for homopolymer filtering */
			private static final int banshift=4;

			@Override
			public void run(){
				try{build();}
				catch(InterruptedException problem){
					Thread.currentThread().interrupt();
					failBuild(new RuntimeException("Interrupted at index count/allocation barrier", problem), intercom);
				}catch(RuntimeException|Error problem){failBuild(problem, intercom);}
			}

			/** Disjoint key ranges share one allocation; neither barrier may outlive a failure. */
			private void build() throws InterruptedException{
				for(int i=minChrom; i<=maxChrom; i++){countSizes(i);}

				final Block b;
				synchronized(intercom){
					intercom[1]++;
					if(id==0){
						while(intercom[1]<4 && failure==null){
							intercom.wait();
						}
						if(failure!=null){rethrow(failure);}

						// Count in long so an oversized block fails before int wrap can
						// corrupt the allocation size or strand the fill workers.
						long sumL=0;
						for(int i=0; i<sizes.length; i++){
							int temp=sizes[i];
							sizes[i]=(int)sumL;
							sumL+=temp;
						}
						if(sumL>Shared.MAX_ARRAY_LEN){
							KillSwitch.kill("Block chr"+minChrom+"-"+maxChrom+" needs "+sumL+
								" k-mer slots, exceeding MAX_ARRAY_LEN="+Shared.MAX_ARRAY_LEN+
								". Use fewer chromosomes per block (lower chrombits) or split chromosomes shorter.");
						}
						final int sum=(int)sumL;

						if(USE_ALLOC_SYNC){
							synchronized(ALLOC_SYNC){//To allow contiguous memory allocation
								b=new Block(KillSwitch.allocInt1D(sum), sizes);
							}
						}else{
							b=new Block(KillSwitch.allocInt1D(sum), sizes);
						}
						indexHolder[0]=b;
						intercom[2]++;
						assert(intercom[2]==1) : "Only worker zero allocates the shared block; allocations="+intercom[2];
						intercom.notifyAll();
					}else{
						while(intercom[2]<1 && failure==null){
							if(intercom[1]>=4){intercom.notifyAll();}
							intercom.wait();
						}
						if(failure!=null){rethrow(failure);}
					}
				}

				for(int i=minChrom; i<=maxChrom; i++){fillArrays(i);}
			}

			/**
			 * First phase: counts occurrences of each k-mer in the specified chromosome.
			 * Scans chromosome sequence and increments counters for valid k-mers.
			 * Applies filtering for homopolymers and modulo constraints.
			 * @param chrom Chromosome number to process
			 */
			private void countSizes(final int chrom){

				ChromosomeArray ca=dna.Data.getChromosome(chrom);

				//			int baseChrom=baseChrom(chrom);

				if(ca.maxIndex>MAX_ALLOWED_CHROM_INDEX){
					throw new RuntimeException("Chrom "+chrom+": "+ca.maxIndex+" > "+MAX_ALLOWED_CHROM_INDEX);
				}

				final int skip=KEYLEN-1;
				assert(skip>0) : "\n*** The key length is too short.  For the flag set 'k=X', X should be between 8 and 15; it was set to "+KEYLEN+" ***\n";

				// Preserve the donor's exact window range: the old a<maxIndex-KEYLEN+1
				// loop examined through base maxIndex-1. Rolling avoids rebuilding all
				// KEYLEN bases at every valid start while retaining identical ownership.
				final byte[] array=ca.array;
				final byte[] baseToNumber=AminoAcid.baseToNumber;
				final int mask=KEYSPACE-1;
				int key=0, valid=0;
				for(int b=ca.minIndex; b<ca.maxIndex; b++){
					final int x=baseToNumber[array[b]];
					if(x<0){key=0; valid=0; continue;}
					key=((key<<2)|x)&mask;
					if(valid<KEYLEN){valid++;}
					if(valid<KEYLEN || key<minIndex || key>maxIndex){continue;}
					if((key>>banshift)!=(key&banmask) && (!USE_MODULO || key%MODULO==0 ||
							(AminoAcid.reverseComplementBinaryFast(key, KEYLEN))%MODULO==0)){
						//TODO: Probable extreme-input overflow - a single key with >Integer.MAX_VALUE
						// occurrences wraps before the long total-size check. Ordinary supported
						// blocks fit int storage; an oversized repetitive block needs an earlier failure.
						sizes[key]++;
					}
				}

			}

			/**
			 * Second phase: populates index arrays with genomic locations for each k-mer.
			 * Scans chromosome sequence and stores encoded (site, chromosome) pairs.
			 * @param chrom Chromosome number to process
			 */
			private void fillArrays(final int chrom){

				ChromosomeArray ca=dna.Data.getChromosome(chrom);

				int baseChrom=baseChrom(chrom);

				if(ca.maxIndex>MAX_ALLOWED_CHROM_INDEX){
					throw new RuntimeException("Chrom "+chrom+": "+ca.maxIndex+" > "+MAX_ALLOWED_CHROM_INDEX);
				}

				final int skip=KEYLEN-1;
				assert(skip>0) : "fillArrays requires KEYLEN>1 because site=b-(KEYLEN-1); KEYLEN="+KEYLEN;

				final int[] sites=indexHolder[0].sites;
				final byte[] array=ca.array;
				final byte[] baseToNumber=AminoAcid.baseToNumber;
				final int mask=KEYSPACE-1;
				int key=0, valid=0;
				for(int b=ca.minIndex; b<ca.maxIndex; b++){
					final int x=baseToNumber[array[b]];
					if(x<0){key=0; valid=0; continue;}
					key=((key<<2)|x)&mask;
					if(valid<KEYLEN){valid++;}
					if(valid<KEYLEN || key<minIndex || key>maxIndex){continue;}
					if((key>>banshift)!=(key&banmask) && (!USE_MODULO || key%MODULO==0 ||
							(AminoAcid.reverseComplementBinaryFast(key, KEYLEN))%MODULO==0)){
						final int a=b-skip;
						final int number=toNumber(a, chrom);
						assert(numberToChrom(number, baseChrom)==chrom) :
							"packed chromosome must round-trip for BBIndex block traversal; chrom="+chrom+", site="+a;
						assert(numberToSite(number)==a) :
							"packed site must round-trip before entering Block.sites; site="+a+", number="+number;
						final int loc=sizes[key];
						assert(sites[loc]==0) : "count/fill ownership mismatch at key="+key+", loc="+loc;
						sites[loc]=number;
						sizes[key]++;
					}
				}

			}

		}

		/** Encode a (location, chrom) pair to an index */
		public final int toNumber(int site, int chrom){
			int out=(chrom&CHROM_MASK_LOW);
			out=out<<SHIFT_LENGTH;
			out=(out|site);
			return out;
		}

		/** Decode an index to a location */
		public final int numberToSite(int number){
			return (number&SITE_MASK);
		}

		/** Decode an (index, baseChrom) pair to a chromosome */
		public final int numberToChrom(int number, int baseChrom){
			assert((baseChrom&CHROM_MASK_LOW)==0) : Integer.toHexString(number)+", baseChrom="+baseChrom;
			assert(baseChrom>=0) : Integer.toHexString(number)+", baseChrom="+baseChrom;

			int out=(number>>>SHIFT_LENGTH);

			out=out+(baseChrom&CHROM_MASK_HIGH);

			return out;
		}

		/**
		 * Calculates base chromosome for encoding operations.
		 * @param chrom Input chromosome number
		 * @return Base chromosome after applying high mask
		 */
		public final int baseChrom(int chrom){return Tools.max(0, chrom&CHROM_MASK_HIGH);}

		/** Length of k-mers for indexing */
		final int KEYLEN;
		/** Number of bits allocated for chromosome encoding */
		private final int CHROMBITS;
		/** Total number of possible k-mers (4^k) */
		private final int KEYSPACE;
		/** Maximum allowed chromosome index to prevent overflow */
		final int MAX_ALLOWED_CHROM_INDEX;
		/** Whether to cache index blocks to disk for reuse */
		public final boolean WRITE_TO_DISK;
		/** Whether to ignore existing disk caches and regenerate */
		public final boolean DISK_INVALID;

		/** Bitmask for extracting low chromosome bits */
		private final int CHROM_MASK_LOW;
		/** Bitmask for extracting high chromosome bits */
		private final int CHROM_MASK_HIGH;
		/** Bitmask for extracting genomic site positions */
		private final int SITE_MASK;
		/** Number of bits to shift for encoding operations */
		private final int SHIFT_LENGTH;

		/** Minimum chromosome number for this block */
		final int minChrom;
		/** Maximum chromosome number for this block */
		final int maxChrom;

		/** Storage matrix for completed index blocks */
		private final Block[] matrix;

	}

	/**
	 * Calculates minimum chromosome for a given range and mask.
	 *
	 * @param chrom Input chromosome number
	 * @param MINCHROM Minimum allowed chromosome
	 * @param CHROM_MASK_HIGH High chromosome bitmask
	 * @return Effective minimum chromosome for processing
	 */
	public static final int minChrom(int chrom, int MINCHROM, int CHROM_MASK_HIGH){return Tools.max(MINCHROM, chrom&CHROM_MASK_HIGH);}
	/**
	 * Calculates maximum chromosome for a given range and mask.
	 *
	 * @param chrom Input chromosome number
	 * @param MINCHROM Minimum allowed chromosome
	 * @param MAXCHROM Maximum allowed chromosome
	 * @param CHROM_MASK_LOW Low chromosome bitmask
	 * @return Effective maximum chromosome for processing
	 */
	public static final int maxChrom(int chrom, int MINCHROM, int MAXCHROM, int CHROM_MASK_LOW){return Tools.max(MINCHROM, Tools.min(MAXCHROM, chrom|CHROM_MASK_LOW));}

	/**
	 * Generates filename for index cache files using current genome build.
	 *
	 * @param minChrom Minimum chromosome number
	 * @param maxChrom Maximum chromosome number
	 * @param k K-mer length
	 * @param chrombits Chromosome encoding bits
	 * @return Filename for cached index file
	 */
	public static final String fname(int minChrom, int maxChrom, int k, int chrombits){
		return fname(minChrom, maxChrom, k, chrombits, Data.GENOME_BUILD);
	}

	/**
	 * Generates filename for index cache files with specific genome build.
	 *
	 * @param minChrom Minimum chromosome number
	 * @param maxChrom Maximum chromosome number
	 * @param k K-mer length
	 * @param chrombits Chromosome encoding bits
	 * @param build Genome build number
	 * @return Filename for cached index file
	 */
	public static final String fname(int minChrom, int maxChrom, int k, int chrombits, int build){
		String suffix="_index_k"+k+"_c"+chrombits+"_b"+build+".block";
		if(minChrom!=maxChrom){
			return Data.ROOT_INDEX+build+"/chr"+minChrom+"-"+maxChrom+suffix;
		}else{
			return Data.ROOT_INDEX+build+"/chr"+minChrom+suffix;
		}
	}

	/**
	 * Thread-safe method to increment the count of active index blocks.
	 * Enforces concurrency limits and provides thread coordination.
	 * @param i Increment amount (positive or negative)
	 */
	static void incrementActiveBlocks(int i){
		assert(i!=0) : "Slot accounting requires a nonzero acquire/release delta";
		synchronized(THREAD_SYNC){
			assert(ACTIVE_BLOCKS>=0) : "Every block release must match a reservation; active="+ACTIVE_BLOCKS;
			assert(ACTIVE_BLOCKS<=MAX_CONCURRENT_BLOCKS) : "Build concurrency must not be reconfigured during an active batch: active="+ACTIVE_BLOCKS+", limit="+MAX_CONCURRENT_BLOCKS;

			while(i>0 && ACTIVE_BLOCKS>0 && ACTIVE_BLOCKS>=MAX_CONCURRENT_BLOCKS){
				try{
					THREAD_SYNC.wait(10000);
				}catch(InterruptedException e){
					Thread.currentThread().interrupt();
					throw new RuntimeException("Interrupted while reserving an index block builder slot", e);
				}
			}
			ACTIVE_BLOCKS+=i;
			if(ACTIVE_BLOCKS<MAX_CONCURRENT_BLOCKS || i<0){THREAD_SYNC.notifyAll();}

			assert(ACTIVE_BLOCKS>=0) : "Block slot released without reservation: active="+ACTIVE_BLOCKS;
			assert(ACTIVE_BLOCKS<=MAX_CONCURRENT_BLOCKS) : "Reservation exceeded the block concurrency limit: active="+ACTIVE_BLOCKS+", limit="+MAX_CONCURRENT_BLOCKS;
		}
	}

	/** Whether to print detailed debugging output */
	public static boolean verbose=false;

	/** Whether to synchronize memory allocation for contiguous arrays */
	public static boolean USE_ALLOC_SYNC=false;
	/** Synchronization object for memory allocation coordination */
	static final String ALLOC_SYNC=new String("ALLOC_SYNC");
	/** Synchronization object for thread coordination */
	private static final String THREAD_SYNC=new String("THREAD_SYNC");

	/** Maximum number of index blocks that can be processed simultaneously */
	public static int MAX_CONCURRENT_BLOCKS=(Shared.LOW_MEMORY ? 1 : (Shared.WINDOWS ? 1 : Tools.max(1, Shared.threads()/4)));
	/** Current number of active index blocks being processed */
	private static int ACTIVE_BLOCKS=0;

	/** Disables periodic-key filtering except for the zero (poly-A) key. */
	public static boolean ALLOW_POLYMERS=false;
	/** Whether to apply modulo filtering to k-mers */
	public static boolean USE_MODULO=false;
	/** Modulo value for k-mer filtering when USE_MODULO is enabled */
	static final int MODULO=9;

}
