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

			private final int id;
			private final int idb;
			private final int[] sizes;
			/** {sizeSum, #finishedCounting, #finishedAllocating, #finishedFilling} */
			private final int[] intercom;
			private final Block[] indexHolder;
			private final int minIndex;
			private final int maxIndex;
			private final int banned;
			private final int banmask;
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

		public final int baseChrom(int chrom){return Tools.max(0, chrom&CHROM_MASK_HIGH);}

		final int KEYLEN;
		private final int CHROMBITS;
		private final int KEYSPACE;
		final int MAX_ALLOWED_CHROM_INDEX;
		public final boolean WRITE_TO_DISK;
		public final boolean DISK_INVALID;

		private final int CHROM_MASK_LOW;
		private final int CHROM_MASK_HIGH;
		private final int SITE_MASK;
		private final int SHIFT_LENGTH;

		final int minChrom;
		final int maxChrom;

		private final Block[] matrix;

	}

	public static final int minChrom(int chrom, int MINCHROM, int CHROM_MASK_HIGH){return Tools.max(MINCHROM, chrom&CHROM_MASK_HIGH);}
	public static final int maxChrom(int chrom, int MINCHROM, int MAXCHROM, int CHROM_MASK_LOW){return Tools.max(MINCHROM, Tools.min(MAXCHROM, chrom|CHROM_MASK_LOW));}
	
	public static final String fname(int minChrom, int maxChrom, int k, int chrombits){
		String suffix="_index_k"+k+"_c"+chrombits+"_b"+Data.GENOME_BUILD+".blockB";
		if(minChrom!=maxChrom){
			return Data.ROOT_INDEX+Data.GENOME_BUILD+"/chr"+minChrom+"-"+maxChrom+suffix;
		}else{
			return Data.ROOT_INDEX+Data.GENOME_BUILD+"/chr"+minChrom+suffix;
		}
	}
	
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

	public static boolean verbose=false;

	public static boolean USE_ALLOC_SYNC=false;
	static final String ALLOC_SYNC=new String("ALLOC_SYNC");
	private static final String THREAD_SYNC=new String("THREAD_SYNC");
	
	public static int MAX_CONCURRENT_BLOCKS=(Shared.WINDOWS ? 1 : 2);
	private static int ACTIVE_BLOCKS=0;

	public static boolean ALLOW_POLYMERS=false;
	public static boolean USE_MODULO=false;
	private static final int MODULO=IndexMaker4.MODULO;
	
}
