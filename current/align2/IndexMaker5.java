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
 * @date Jan 3, 2013
 *
 */
public class IndexMaker5 {
	
	
	/**
	 * Creates k-mer indices for specified chromosome range using parallel processing.
	 * Each chromosome range is processed by separate BlockMaker threads to generate
	 * hash-based lookup structures for sequence alignment.
	 *
	 * @param genome Genome build identifier for data loading
	 * @param minChrom Starting chromosome number (inclusive)
	 * @param maxChrom Ending chromosome number (inclusive)
	 * @param k K-mer length for indexing
	 * @param CHROMBITS Number of bits allocated for chromosome encoding
	 * @param MAX_ALLOWED_CHROM_INDEX Maximum allowed chromosome index value
	 * @param CHROM_MASK_LOW Low bits mask for chromosome encoding
	 * @param CHROM_MASK_HIGH High bits mask for chromosome encoding
	 * @param SITE_MASK Mask for site position encoding
	 * @param SHIFT_LENGTH Bit shift length for position encoding
	 * @param WRITE Whether to write index to disk
	 * @param DISK_INVALID Whether disk index is invalid and needs regeneration
	 * @param index Existing index array to populate (created if null)
	 * @return Array of Block objects containing the generated indices
	 */
	public static Block[] makeIndex(final int genome, int minChrom, int maxChrom, int k, int CHROMBITS,
			int MAX_ALLOWED_CHROM_INDEX, int CHROM_MASK_LOW, int CHROM_MASK_HIGH, int SITE_MASK, int SHIFT_LENGTH, boolean WRITE, boolean DISK_INVALID, Block[] index){
		Timer t=new Timer();
		
		MAX_CONCURRENT_BLOCKS=(Shared.WINDOWS ? 1 : Tools.max(1, Shared.threads()/4));
		
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
	 * Creates a single Block for the specified chromosome range.
	 * This method contains assertion failures and appears to be for debugging purposes.
	 *
	 * @param minChrom Starting chromosome number
	 * @param maxChrom Ending chromosome number
	 * @param k K-mer length for indexing
	 * @param CHROMBITS Number of bits for chromosome encoding
	 * @param MAX_ALLOWED_CHROM_INDEX Maximum chromosome index allowed
	 * @param CHROM_MASK_LOW Low chromosome mask
	 * @param CHROM_MASK_HIGH High chromosome mask
	 * @param SITE_MASK Site position mask
	 * @param SHIFT_LENGTH Position encoding shift length
	 * @param WRITE Whether to write to disk
	 * @param DISK_INVALID Whether disk version is invalid
	 * @param matrix Block matrix to update
	 * @return Generated Block object
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
	
	
	
	private static class BlockMaker extends Thread{

		/**
		 * Constructs a BlockMaker thread for processing a chromosome range.
		 * Initializes all bit masks and encoding parameters for k-mer indexing.
		 *
		 * @param minChrom_ Starting chromosome number
		 * @param maxChrom_ Ending chromosome number
		 * @param k K-mer length
		 * @param CHROMBITS_ Chromosome encoding bits
		 * @param MAX_ALLOWED_CHROM_INDEX_ Maximum chromosome index
		 * @param CHROM_MASK_LOW_ Low chromosome mask
		 * @param CHROM_MASK_HIGH_ High chromosome mask
		 * @param SITE_MASK_ Site position mask
		 * @param SHIFT_LENGTH_ Position shift length
		 * @param WRITE_TO_DISK_ Whether to write results to disk
		 * @param DISK_INVALID_ Whether disk version is invalid
		 * @param matrix_ Block matrix to populate
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
		 * First attempts to load existing index from disk, otherwise generates
		 * new index using parallel CountThread workers. Coordinates thread
		 * synchronization for size counting, array allocation, and filling phases.
		 *
		 * @return Generated Block containing the k-mer index arrays
		 */
		Block makeArrays(){
			
			{
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
			int[] sizes=new int[KEYSPACE+1];
			int[] intercom=KillSwitch.allocInt1D(4);
			Block[] indexHolder=new Block[1];

			for(int i=0; i<4; i++){
				threads[i]=new CountThread(i, sizes, intercom, indexHolder);
				threads[i].start();
//				while(!threads[i].isAlive()){
//					//wait for these threads to start
//				}
			}
			Data.sysout.println("Indexing threads started.");
			for(int i=0; i<threads.length; i++){
				if(threads[i].getState()!=State.TERMINATED){
					try {
						threads[i].join();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			Data.sysout.println("Threads finished.");
			
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


		private class CountThread extends Thread{

			/**
			 * Constructs a CountThread for parallel k-mer index generation.
			 * Each thread handles k-mers starting with a specific base (A,C,G,T).
			 * Calculates index range and polymeric k-mer banning parameters.
			 *
			 * @param id_ Thread identifier (0-3 for A,C,G,T respectively)
			 * @param sizes_ Shared array for k-mer counts per index
			 * @param intercom_ Inter-thread communication array for synchronization
			 * @param indexHolder_ Shared holder for the final Block object
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
					banned=-1;
					banmask=-1; //poly-A still slips through
				}else{
					int b=0;
					for(int i=0; i<KEYLEN; i++){
						b<<=2;
						b=(b|id);
					}
					banned=b;
					banmask=~((-1)<<((2*KEYLEN)-banshift));
				}
			}

			/** Thread identifier (0-3 for bases A,C,G,T respectively) */
			private final int id;
			/** Byte value of the base this thread processes */
			private final int idb;
			/** Shared array tracking k-mer counts for each index position */
			private final int[] sizes;
			/** {sizeSum, #finishedCounting, #finishedAllocating, #finishedFilling} */
			private final int[] intercom;
			/** Shared holder for the final generated Block object */
			private final Block[] indexHolder;
			/** Minimum k-mer index value this thread handles */
			private final int minIndex;
			/** Maximum k-mer index value this thread handles */
			private final int maxIndex;
			/** Bitmask pattern for polymeric k-mers to exclude from indexing */
			private final int banned;
			/** Bit mask used for polymeric k-mer detection and filtering */
			private final int banmask;
			/** Bit shift value for polymeric k-mer ban mask calculations */
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
								b=new Block(new int[sum], sizes);
							}
						}else{
							b=new Block(new int[sum], sizes);
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

			/**
			 * Counts k-mers for a specific chromosome to determine array sizes.
			 * Scans chromosome sequence data and increments size counters for valid
			 * k-mers that match this thread's assigned base and pass filtering criteria.
			 * @param chrom Chromosome number to process
			 */
			private void countSizes(final int chrom){

				//			System.err.println("Thread "+id+" using chr"+chrom+" for countSizes");
				ChromosomeArray ca=dna.Data.getChromosome(chrom);

				//			int baseChrom=baseChrom(chrom);

				if(ca.maxIndex>MAX_ALLOWED_CHROM_INDEX){
					throw new RuntimeException("Chrom "+chrom+": "+ca.maxIndex+" > "+MAX_ALLOWED_CHROM_INDEX);
				}

				final int max=ca.maxIndex-KEYLEN+1;
				final int skip=KEYLEN-1;
				assert(skip>0);


				int start=ca.minIndex;
				while(start<max && ca.getNumber(start+skip)==-1){start+=skip;}
				while(start<max && ca.getNumber(start)==-1){start++;}

				//			Data.sysout.println("Entering hash loop.");

				// "a" is site start, "b" is site end
				final byte[] array=ca.array;
				for(int a=start, b=start+skip; a<max; a++, b++){
					if(array[a]==idb){
						int key=ca.getNumber(a, b);
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

			/**
			 * Fills index arrays with genomic positions for a specific chromosome.
			 * Scans chromosome sequence and stores encoded (position, chromosome) pairs
			 * in the appropriate k-mer index locations based on sequence content.
			 * @param chrom Chromosome number to process
			 */
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
		 * Calculates base chromosome value using high chromosome mask.
		 * Used for chromosome encoding/decoding operations.
		 * @param chrom Chromosome number
		 * @return Base chromosome value after masking
		 */
		public final int baseChrom(int chrom){return Tools.max(0, chrom&CHROM_MASK_HIGH);}

		/** K-mer length for this block maker */
		final int KEYLEN;
		/** Number of bits allocated for chromosome encoding */
		private final int CHROMBITS;
		/** Total number of possible k-mer values (4^KEYLEN) */
		private final int KEYSPACE;
		/** Maximum chromosome index value allowed for processing */
		final int MAX_ALLOWED_CHROM_INDEX;
		/** Whether to write the generated index to disk */
		public final boolean WRITE_TO_DISK;
		/** Whether the existing disk index is invalid and needs regeneration */
		public final boolean DISK_INVALID;

		/** Bit mask for low chromosome bits in encoding */
		private final int CHROM_MASK_LOW;
		/** Bit mask for high chromosome bits in encoding */
		private final int CHROM_MASK_HIGH;
		/** Bit mask for genomic site position in encoding */
		private final int SITE_MASK;
		/** Number of bits to shift for position encoding */
		private final int SHIFT_LENGTH;

		/** Starting chromosome number for this block */
		final int minChrom;
		/** Ending chromosome number for this block */
		final int maxChrom;

		/** Array of Block objects to populate with generated indices */
		private final Block[] matrix;

	}

	/**
	 * Calculates minimum chromosome value using high chromosome mask.
	 *
	 * @param chrom Input chromosome number
	 * @param MINCHROM Minimum allowed chromosome value
	 * @param CHROM_MASK_HIGH High bits mask for chromosome
	 * @return Maximum of MINCHROM and masked chromosome value
	 */
	public static final int minChrom(int chrom, int MINCHROM, int CHROM_MASK_HIGH){return Tools.max(MINCHROM, chrom&CHROM_MASK_HIGH);}
	/**
	 * Calculates maximum chromosome value using low chromosome mask.
	 *
	 * @param chrom Input chromosome number
	 * @param MINCHROM Minimum allowed chromosome value
	 * @param MAXCHROM Maximum allowed chromosome value
	 * @param CHROM_MASK_LOW Low bits mask for chromosome
	 * @return Constrained chromosome value within allowed range
	 */
	public static final int maxChrom(int chrom, int MINCHROM, int MAXCHROM, int CHROM_MASK_LOW){return Tools.max(MINCHROM, Tools.min(MAXCHROM, chrom|CHROM_MASK_LOW));}
	
	/**
	 * Generates filename for index storage based on chromosome range and parameters.
	 *
	 * @param minChrom Starting chromosome number
	 * @param maxChrom Ending chromosome number
	 * @param k K-mer length
	 * @param chrombits Number of chromosome encoding bits
	 * @return Filename string for storing the index
	 */
	public static final String fname(int minChrom, int maxChrom, int k, int chrombits){
		String suffix="_index_k"+k+"_c"+chrombits+"_b"+Data.GENOME_BUILD+".blockB";
		if(minChrom!=maxChrom){
			return Data.ROOT_INDEX+Data.GENOME_BUILD+"/chr"+minChrom+"-"+maxChrom+suffix;
		}else{
			return Data.ROOT_INDEX+Data.GENOME_BUILD+"/chr"+minChrom+suffix;
		}
	}
	
	/**
	 * Thread-safe counter for active block processing threads.
	 * Implements blocking when maximum concurrent blocks is reached.
	 * @param i Increment value (positive to increase, negative to decrease)
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

	/** Whether to output verbose debugging information */
	public static boolean verbose=false;

	/** Whether to synchronize memory allocation for contiguous arrays */
	public static boolean USE_ALLOC_SYNC=false;
	/** Synchronization object for memory allocation operations */
	static final String ALLOC_SYNC=new String("ALLOC_SYNC");
	/** Synchronization object for thread coordination operations */
	private static final String THREAD_SYNC=new String("THREAD_SYNC");
	
	/** Maximum number of blocks that can be processed concurrently */
	public static int MAX_CONCURRENT_BLOCKS=(Shared.WINDOWS ? 1 : 2);
	/** Current number of actively processing blocks */
	private static int ACTIVE_BLOCKS=0;

	/** Whether to allow polymeric k-mers in the index */
	public static boolean ALLOW_POLYMERS=false;
	/** Whether to use modulo filtering for k-mer selection */
	public static boolean USE_MODULO=false;
	/** Modulo value used for k-mer filtering when enabled */
	private static final int MODULO=IndexMaker4.MODULO;
	
}
