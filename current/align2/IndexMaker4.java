package align2;

import java.io.File;
import java.lang.Thread.State;
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


/**
 * @author Brian Bushnell
 * @date Dec 23, 2012
 *
 */
public class IndexMaker4 {
	
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
		
		assert(minChrom<=maxChrom);
		
		if(index==null){index=new Block[maxChrom+1];}
		
		ArrayList<BlockMaker> list=new ArrayList<BlockMaker>();
		
		for(int i=1; i<=maxChrom;){
			if(i>=minChrom){
				int a=minChrom(i, minChrom, CHROM_MASK_HIGH);
				int b=maxChrom(i, minChrom, maxChrom, CHROM_MASK_LOW);
				assert(b>=i);
				
				BlockMaker idm=new BlockMaker(a, b, k, CHROMBITS, MAX_ALLOWED_CHROM_INDEX, CHROM_MASK_LOW, CHROM_MASK_HIGH, SITE_MASK, SHIFT_LENGTH, WRITE, DISK_INVALID, index);
				list.add(idm);
				incrementActiveBlocks(1);
				idm.start();
				
				while(idm.getState()==State.NEW){}//wait
				
				i=b+1;
			}else{i++;}
		}
		
		for(BlockMaker cm : list){
			while(cm.getState()!=State.TERMINATED){
				try {
					cm.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		t.stop();
//		Data.sysout.println("Index gen time: \t"+t);
		
		return index;
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
//			assert(false) : maxChrom+", "+MAX_ALLOWED_CHROM_INDEX;
//			System.err.println(minChrom+"~"+maxChrom);
		}


		@Override
		public void run(){
			makeArrays();
			incrementActiveBlocks(-1);
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

			for(int i=0; i<4; i++){
				threads[i]=new CountThread(i, sizes, intercom, indexHolder);
				threads[i].start();
//				while(!threads[i].isAlive()){
//					//wait for these threads to start
//				}
			}
			Data.sysout.println("Indexing threads started for block "+baseChrom(minChrom)+"-"+maxChrom);
			for(int i=0; i<threads.length; i++){
				while(threads[i].getState()!=State.TERMINATED){
					try {
						threads[i].join();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
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
//				assert(!f.exists()) : "Tried to overwrite file "+f.getAbsolutePath();
				indexHolder[0].write(fname, true);
			}
			
			return indexHolder[0];
		}


		/**
		 * Worker thread that processes k-mers for a specific base (A, C, G, T) range.
		 * Handles both counting phase (determining array sizes) and filling phase (populating indices).
		 * Uses chromosome-specific filtering to avoid homopolymer sequences.
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
				idb=AminoAcid.numberToBase[id];
				sizes=sizes_;
				indexHolder=indexHolder_;
				intercom=intercom_;

				minIndex=(id<<(2*KEYLEN-2));
				maxIndex=(int)(((id+1L)<<(2*KEYLEN-2))-1);
				//Data.sysout.println("Thread "+id+" range is "+minIndex+", "+maxIndex);
				
				if(ALLOW_POLYMERS){
					banmask=-1; //poly-A still slips through
				}else{
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
			/** Byte value of the base this thread processes */
			private final int idb;
			/** Shared array for counting k-mer occurrences */
			private final int[] sizes;
			/** {sizeSum, #finishedCounting, #finishedAllocating, #finishedFilling} */
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

				//Data.sysout.println("Thread "+id+" counting sizes for ("+minChrom+", "+maxChrom+")");
				for(int i=minChrom; i<=maxChrom; i++){countSizes(i);}
				
				final Block b;
				synchronized(intercom){
					//Data.sysout.println("Thread "+id+" synced on intercom: "+Arrays.toString(intercom));
					intercom[1]++;
					if(id==0){
						while(intercom[1]<4){
							//Data.sysout.println("Thread "+id+" waiting on intercom: "+Arrays.toString(intercom));
							try {
								intercom.wait();
							} catch (InterruptedException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
						}
						
						int sum=0;
						for(int i=0; i<sizes.length; i++){
							int temp=sizes[i];
							sizes[i]=sum;
							sum+=temp;
						}
						
						if(USE_ALLOC_SYNC){
							synchronized(ALLOC_SYNC){//To allow contiguous memory allocation
								b=new Block(KillSwitch.allocInt1D(sum), sizes);
							}
						}else{
							b=new Block(KillSwitch.allocInt1D(sum), sizes);
						}
						indexHolder[0]=b;
						intercom[2]++;
						assert(intercom[2]==1);
						intercom.notifyAll();
					}else{
						while(intercom[2]<1){
							//Data.sysout.println("Thread "+id+" waiting on intercom: "+Arrays.toString(intercom));
							try {
								if(intercom[1]>=4){intercom.notify();}
								intercom.wait();
							} catch (InterruptedException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
						}
					}
				}

				//Data.sysout.println("Thread "+id+" filling arrays for ("+minChrom+", "+maxChrom+")");

				
				for(int i=minChrom; i<=maxChrom; i++){fillArrays(i);}
				//Data.sysout.println("Thread "+id+" finished.");
			}

			private void countSizes(final int chrom){

				//			System.err.println("Thread "+id+" using chr"+chrom+" for countSizes");
				ChromosomeArray ca=dna.Data.getChromosome(chrom);

				//			int baseChrom=baseChrom(chrom);

				if(ca.maxIndex>MAX_ALLOWED_CHROM_INDEX){
					throw new RuntimeException("Chrom "+chrom+": "+ca.maxIndex+" > "+MAX_ALLOWED_CHROM_INDEX);
				}

				final int max=ca.maxIndex-KEYLEN+1;
				final int skip=KEYLEN-1;
				assert(skip>0) : "\n*** The key length is too short.  For the flag set 'k=X', X should be between 8 and 15; it was set to "+KEYLEN+" ***\n";


				int start=ca.minIndex;
				while(start<max && ca.getNumber(start+skip)==-1){start+=skip;}
				while(start<max && ca.getNumber(start)==-1){start++;}

				//			Data.sysout.println("Entering hash loop.");

				// "a" is site start, "b" is site end
				final byte[] array=ca.array;
				for(int a=start, b=start+skip; a<max; a++, b++){
					if(array[a]==idb){
						int key=ca.getNumber(a, b);
//						if(key>=0 && (key>>banshift)!=(key&banmask) && (!USE_MODULO || key%MODULO==0)){
//							assert(key>=minIndex && key<=maxIndex) : "\n"+id+", "+ca.getNumber(a)+", "+(char)ca.get(a)+", "+key+", "+Integer.toHexString(key)+
//							", "+ca.getString(a, b)+"\n"+minIndex+", "+maxIndex+"\n";
//							sizes[key]++;
//						}
						if(key>=0 && (key>>banshift)!=(key&banmask) && (!USE_MODULO || key%MODULO==0 || (AminoAcid.reverseComplementBinaryFast(key, KEYLEN))%MODULO==0)){
							assert(key>=minIndex && key<=maxIndex) : "\n"+id+", "+ca.getNumber(a)+", "+(char)ca.get(a)+", "+key+", "+Integer.toHexString(key)+
							", "+ca.getString(a, b)+"\n"+minIndex+", "+maxIndex+"\n";
							sizes[key]++;
						}
					}
					//				Data.sysout.println("a="+a+", b="+b+", max="+max);
				}

				//			Data.sysout.println("Left hash loop.");

			}

			private void fillArrays(final int chrom){

				//			System.err.println("Thread "+id+" using chr"+chrom+" for fillArrays");
				ChromosomeArray ca=dna.Data.getChromosome(chrom);

				int baseChrom=baseChrom(chrom);

				if(ca.maxIndex>MAX_ALLOWED_CHROM_INDEX){
					throw new RuntimeException("Chrom "+chrom+": "+ca.maxIndex+" > "+MAX_ALLOWED_CHROM_INDEX);
				}

				final int max=ca.maxIndex-KEYLEN+1;
				final int skip=KEYLEN-1;
				assert(skip>0);


				int start=ca.minIndex;
				while(start<max && ca.getNumber(start+skip)==-1){start+=skip;}
				while(start<max && ca.getNumber(start)==-1){start++;}


//				//			Data.sysout.println("Entering hash loop.");
//				// "a" is site start, "b" is site end
//				int len=KEYLEN-1;
//				int keyB=ca.getNumber(start, start+skip-1);
//				final int mask=(KEYLEN==16 ? -1 : ~((-1)<<(2*KEYLEN)));
//				final byte[] array=ca.array;
//				final byte[] btn=AminoAcid.baseToNumber;
//				for(int a=start, b=start+skip; a<max; a++, b++){
//					int c=btn[array[b]];
//					if(c>=0){
//						keyB=((keyB<<2)|c);
//						len++;
//					}else{
//						len=0;
//					}
//					int key=keyB&mask;
//					if(len>=KEYLEN && /* array[a]==idb*/ key>=minIndex && key<=maxIndex){
////						int key=keyB&mask;
//						assert(key>=minIndex && key<=maxIndex);
//						int number=toNumber(a, chrom);
//						assert(numberToChrom(number, baseChrom)==chrom);
//						assert(numberToSite(number)==a);
//						index[key][sizes[key]]=number;
//						sizes[key]++;
//					}
//					//				Data.sysout.println("a="+a+", b="+b+", max="+max);
//				}


				//			Data.sysout.println("Entering hash loop.");
				// "a" is site start, "b" is site end
				
				int[] sites=indexHolder[0].sites;
				
				for(int a=start, b=start+skip; a<max; a++, b++){
					if(ca.array[a]==idb){
						int key=ca.getNumber(a, b);
						if(key>=0 && (key>>banshift)!=(key&banmask) && (!USE_MODULO || key%MODULO==0 || (AminoAcid.reverseComplementBinaryFast(key, KEYLEN))%MODULO==0)){
							assert(key>=minIndex && key<=maxIndex);
							int number=toNumber(a, chrom);
							assert(numberToChrom(number, baseChrom)==chrom);
							assert(numberToSite(number)==a);
							int loc=sizes[key];
							assert(sites[loc]==0);
							sites[loc]=number;
							sizes[key]++;
						}
					}
					//				Data.sysout.println("a="+a+", b="+b+", max="+max);
				}
				//			Data.sysout.println("Left hash loop.");

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
			//		assert(baseChrom<8) : Integer.toHexString(number)+", baseChrom="+baseChrom;

			int out=(number>>>SHIFT_LENGTH);

			out=out+(baseChrom&CHROM_MASK_HIGH);

			//		assert(out<8) : Integer.toHexString(number)+", baseChrom="+baseChrom;
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
		assert(i!=0);
		synchronized(THREAD_SYNC){
			assert(ACTIVE_BLOCKS>=0);
			assert(ACTIVE_BLOCKS<=MAX_CONCURRENT_BLOCKS);
			
			while(i>0 && ACTIVE_BLOCKS>0 && ACTIVE_BLOCKS>=MAX_CONCURRENT_BLOCKS){
				try {
					THREAD_SYNC.wait(10000);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			ACTIVE_BLOCKS+=i;
			if(ACTIVE_BLOCKS<MAX_CONCURRENT_BLOCKS || i<0){THREAD_SYNC.notifyAll();}
			
			assert(ACTIVE_BLOCKS>=0);
			assert(ACTIVE_BLOCKS<=MAX_CONCURRENT_BLOCKS);
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

	/** Whether to allow homopolymer k-mers in the index */
	public static boolean ALLOW_POLYMERS=false;
	/** Whether to apply modulo filtering to k-mers */
	public static boolean USE_MODULO=false;
	/** Modulo value for k-mer filtering when USE_MODULO is enabled */
	static final int MODULO=9;
	
}
