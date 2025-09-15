package kmer;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ArrayBlockingQueue;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import jgi.BBMerge;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.IntList;
import structures.ListNum;

/**
 * @author Brian Bushnell
 * @date Mar 4, 2015
 *
 */
public class TableLoaderLockFree {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Timer t=new Timer();
		
		AbstractKmerTable[] tables=makeTables(AbstractKmerTable.ARRAY1D, 12, -1L, false, 1.0);
		
		int k=31;
		int mink=0;
		int speed=0;
		int hdist=0;
		int edist=0;
		boolean rcomp=true;
		boolean maskMiddle=false;
		
		//Create a new Loader instance
		TableLoaderLockFree loader=new TableLoaderLockFree(tables, k, mink, speed, hdist, edist, rcomp, maskMiddle);
		loader.setRefSkip(0);
		loader.hammingDistance2=0;
		loader.editDistance2=0;
		loader.storeMode(SET_IF_NOT_PRESENT);
		
		///And run it
		String[] refs=args;
		String[] literals=null;
		boolean keepNames=false;
		boolean useRefNames=false;
		long kmers=loader.processData(refs, literals, keepNames, useRefNames, false);
		t.stop();

		outstream.println("Time:     \t"+t);
		outstream.println("Return:   \t"+kmers);
		outstream.println("refKmers: \t"+loader.refKmers);
		outstream.println("refBases: \t"+loader.refBases);
		outstream.println("refReads: \t"+loader.refReads);
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	public TableLoaderLockFree(AbstractKmerTable[] tables_, int k_){
		this(tables_, k_, 0, 0, 0, 0, true, false);
	}
	
	public TableLoaderLockFree(AbstractKmerTable[] tables_, int k_, int mink_, int speed_, int hdist_, int edist_, boolean rcomp_, boolean maskMiddle_){
		tables=tables_;
		k=k_;
		k2=k-1;
		mink=mink_;
		rcomp=rcomp_;
		useShortKmers=(mink>0 && mink<k);
		speed=speed_;
		hammingDistance=hdist_;
		editDistance=edist_;
		middleMask=maskMiddle ? ~(3L<<(2*(k/2))) : -1L;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	
//	public static AbstractKmerTable[] makeTables(int tableType, int initialSize, long coreMask, boolean growable){
//		return AbstractKmerTable.preallocate(WAYS, tableType, initialSize, coreMask, growable);
//	}
	
//	public static AbstractKmerTable[] makeTables(int tableType, int[] schedule, long coreMask){
//		return AbstractKmerTable.preallocate(WAYS, tableType, schedule, coreMask);
//	}
	
	public static AbstractKmerTable[] makeTables(int tableType, int bytesPerKmer, long coreMask, 
			boolean prealloc, double memRatio){
		ScheduleMaker scheduleMaker=new ScheduleMaker(WAYS, bytesPerKmer, prealloc, memRatio);
		int[] schedule=scheduleMaker.makeSchedule();
		return AbstractKmerTable.preallocate(WAYS, tableType, schedule, coreMask);
	}
	
	ScheduleMaker scheduleMaker=new ScheduleMaker(WAYS, 12, false, 0.8);
	int[] schedule=scheduleMaker.makeSchedule();
	
	public long processData(String[] ref, String[] literal, boolean keepNames, boolean useRefNames, boolean ecc_){
		
		scaffoldNames=null;
		refNames=null;
		refScafCounts=null;
		scaffoldLengths=null;
		ecc=ecc_;
		
		if(keepNames){
			scaffoldNames=new ArrayList<String>();
			refNames=new ArrayList<String>();
			scaffoldLengths=new IntList();
			
			if(ref!=null){
				for(String s : ref){refNames.add(s);}
			}
			if(literal!=null){refNames.add("literal");}
			refScafCounts=new int[refNames.size()];
			
			if(useRefNames){toRefNames();}
		}
		
		return spawnLoadThreads(ref, literal);
	}
	
	public void setRefSkip(int x){setRefSkip(x, x);}
	
	public void setRefSkip(int min, int max){
		max=Tools.max(min, max);
		if(min==max){
			minRefSkip=maxRefSkip=min;
		}else{
			minRefSkip=min;
			maxRefSkip=max;
		}
		variableRefSkip=(minRefSkip!=maxRefSkip);
	}
	
	public void storeMode(final int x){
		assert(x==SET_IF_NOT_PRESENT || x==SET_ALWAYS || x==INCREMENT);
		storeMode=x;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	

	/**
	 * Fills tables with kmers from references, using multiple LoadThread.
	 * @return Number of kmers stored.
	 */
	private long spawnLoadThreads(String[] ref, String[] literal){
		Timer t=new Timer();
		if((ref==null || ref.length<1) && (literal==null || literal.length<1)){return 0;}
		long added=0;
		
		/* Create load threads */
		LoadThread[] loaders=new LoadThread[WAYS];
		for(int i=0; i<loaders.length; i++){
			loaders[i]=new LoadThread(i);
			loaders[i].start();
		}
		
		/* For each reference file... */
		int refNum=0;
		if(ref!=null){
			for(String refname : ref){

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
						if(scaffoldNames!=null){
							for(Read r1 : reads2){
								final Read r2=r1.mate;
								final Integer id=scaffoldNames.size();
								if(ecc && r1!=null && r2!=null){BBMerge.findOverlapStrict(r1, r2, true);}
								refScafCounts[refNum]++;
								scaffoldNames.add(r1.id==null ? id.toString() : r1.id);
								int len=r1.length();
								r1.obj=id;
								if(r2!=null){
									r2.obj=id;
									len+=r2.length();
								}
								scaffoldLengths.add(len);
							}
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
		if(literal!=null){
			ArrayList<Read> list=new ArrayList<Read>(literal.length);
			if(verbose){outstream.println("Adding literals "+Arrays.toString(literal));}

			/* Assign a unique ID number to each literal sequence */
			for(int i=0; i<literal.length; i++){
				if(scaffoldNames!=null){
					final Integer id=scaffoldNames.size();
					final Read r=new Read(literal[i].getBytes(), null, id);
					refScafCounts[refNum]++;
					scaffoldNames.add(id.toString());
					scaffoldLengths.add(r.length());
					r.obj=id;
				}else{
					final Read r=new Read(literal[i].getBytes(), null, i);
					list.add(r);
				}
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
		for(LoadThread lt : loaders){
			while(lt.getState()!=Thread.State.TERMINATED){
				try {
					lt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			added+=lt.addedT;
			refKmers+=lt.refKmersT;
			refBases+=lt.refBasesT;
			refReads+=lt.refReadsT;
		}
		//Correct statistics for number of threads, since each thread processes all reference data
		refKmers/=WAYS;
		refBases/=WAYS;
		refReads/=WAYS;

		t.stop();
		if(DISPLAY_PROGRESS){
			outstream.println("Added "+added+" kmers; time: \t"+t);
			Shared.printMemory();
			outstream.println();
		}
		
		if(verbose){
			TextStreamWriter tsw=new TextStreamWriter("stdout", false, false, false, FileFormat.TEXT);
			tsw.start();
			for(AbstractKmerTable table : tables){
				table.dumpKmersAsText(tsw, k, 1, Integer.MAX_VALUE);
			}
			tsw.poisonAndWait();
		}
		
		return added;
	}
	

	
	/**
	 * Fills the scaffold names array with reference names.
	 */
	public void toRefNames(){
		final int numRefs=refNames.size();
		for(int r=0, s=1; r<numRefs; r++){
			final int scafs=refScafCounts[r];
			final int lim=s+scafs;
			final String name=ReadWrite.stripToCore(refNames.get(r));
//			outstream.println("r="+r+", s="+s+", scafs="+scafs+", lim="+lim+", name="+name);
			while(s<lim){
//				outstream.println(r+", "+s+". Setting "+scaffoldNames.get(s)+" -> "+name);
				scaffoldNames.set(s, name);
				s++;
			}
		}
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
			map=tables[tnum];
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
					// TODO Auto-generated catch block
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
					
					addedT+=addToMap(r1, minRefSkip);
					if(r2!=null){addedT+=addToMap(r2, minRefSkip);}
				}
				reads=fetch();
			}
			
			if(map.canRebalance() && map.size()>2L*map.arrayLength()){
				map.rebalance();
			}
		}

		/**
		 * @param r The current read to process
		 * @param skip Number of bases to skip between kmers
		 * @return Number of kmers stored
		 */
		private long addToMap(Read r, int skip){
			if(variableRefSkip){
				int rblen=r.length();
				skip=(rblen>20000000 ? k : rblen>5000000 ? 11 : rblen>500000 ? 2 : 0);
				skip=Tools.mid(minRefSkip, maxRefSkip, skip);
			}
			final byte[] bases=r.bases;
			final int shift=2*k;
			final int shift2=shift-2;
			final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
			final long kmask=kMasks[k];
			long kmer=0;
			long rkmer=0;
			long added=0;
			int len=0;
			
			if(bases!=null){
				refReadsT++;
				refBasesT+=bases.length;
			}
			if(bases==null || bases.length<k){return 0;}
			
			final int id=(r.obj==null ? 1 : ((Integer)r.obj).intValue());
			
			if(skip>1){ //Process while skipping some kmers
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					long x=AminoAcid.baseToNumber[b];
					long x2=AminoAcid.baseToComplementNumber[b];
					kmer=((kmer<<2)|x)&mask;
					rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
					if(x<0){len=0; rkmer=0;}else{len++;}
					if(verbose){outstream.println("Scanning1 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=k){
						refKmersT++;
						if(len%skip==0){
							final long extraBase=(i>=bases.length-1 ? -1 : AminoAcid.baseToNumber[bases[i+1]]);
							added+=addToMap(kmer, rkmer, k, extraBase, id, kmask, hammingDistance, editDistance);
							if(useShortKmers){
								if(i==k2){added+=addToMapRightShift(kmer, rkmer, id);}
								if(i==bases.length-1){added+=addToMapLeftShift(kmer, rkmer, extraBase, id);}
							}
						}
					}
				}
			}else{ //Process all kmers
				for(int i=0; i<bases.length; i++){
					byte b=bases[i];
					long x=AminoAcid.baseToNumber[b];
					long x2=AminoAcid.baseToComplementNumber[b];
					kmer=((kmer<<2)|x)&mask;
					rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
					if(x<0){len=0; rkmer=0;}else{len++;}
					if(verbose){outstream.println("Scanning2 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
					if(len>=k){
						refKmersT++;
						final long extraBase=(i>=bases.length-1 ? -1 : AminoAcid.baseToNumber[bases[i+1]]);
						final long atm=addToMap(kmer, rkmer, k, extraBase, id, kmask, hammingDistance, editDistance);
						added+=atm;
//						assert(false) : atm+", "+map.contains(toValue(kmer, rkmer, kmask));
						if(useShortKmers){
							if(i==k2){added+=addToMapRightShift(kmer, rkmer, id);}
							if(i==bases.length-1){added+=addToMapLeftShift(kmer, rkmer, extraBase, id);}
						}
					}
				}
			}
			return added;
		}
		

		/**
		 * Adds short kmers on the left end of the read.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param extraBase Base added to end in case of deletions
		 * @param id Scaffold number
		 * @return Number of kmers stored
		 */
		private long addToMapLeftShift(long kmer, long rkmer, final long extraBase, final int id){
			if(verbose){outstream.println("addToMapLeftShift");}
			long added=0;
			for(int i=k-1; i>=mink; i--){
				kmer=kmer&rightMasks[i];
				rkmer=rkmer>>>2;
				long x=addToMap(kmer, rkmer, i, extraBase, id, kMasks[i], hammingDistance2, editDistance2);
				added+=x;
				if(verbose){
					if((toValue(kmer, rkmer, kMasks[i]))%WAYS==tnum){
						outstream.println("added="+x+"; i="+i+"; tnum="+tnum+"; Added left-shift kmer "+AminoAcid.kmerToString(kmer&~kMasks[i], i)+"; value="+(toValue(kmer, rkmer, kMasks[i]))+"; kmer="+kmer+"; rkmer="+rkmer+"; kmask="+kMasks[i]+"; rightMasks[i+1]="+rightMasks[i+1]);
						outstream.println("i="+i+"; tnum="+tnum+"; Looking for left-shift kmer "+AminoAcid.kmerToString(kmer&~kMasks[i], i));
						final long value=toValue(kmer, rkmer, kMasks[i]);
						if(map.contains(value)){outstream.println("Found "+value);}
					}
				}
			}
			return added;
		}
		

		/**
		 * Adds short kmers on the right end of the read.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param id Scaffold number
		 * @return Number of kmers stored
		 */
		private long addToMapRightShift(long kmer, long rkmer, final int id){
			if(verbose){outstream.println("addToMapRightShift");}
			long added=0;
			for(int i=k-1; i>=mink; i--){
				long extraBase=kmer&3L;
				kmer=kmer>>>2;
				rkmer=rkmer&rightMasks[i];
//				assert(Long.numberOfLeadingZeros(kmer)>=2*(32-i)) : Long.numberOfLeadingZeros(kmer)+", "+i+", "+kmer+", "+kMasks[i];
//				assert(Long.numberOfLeadingZeros(rkmer)>=2*(32-i)) : Long.numberOfLeadingZeros(rkmer)+", "+i+", "+rkmer+", "+kMasks[i];
				long x=addToMap(kmer, rkmer, i, extraBase, id, kMasks[i], hammingDistance2, editDistance2);
				added+=x;
				if(verbose){
					if((toValue(kmer, rkmer, kMasks[i]))%WAYS==tnum){
						outstream.println("added="+x+"; i="+i+"; tnum="+tnum+"; Added right-shift kmer "+AminoAcid.kmerToString(kmer&~kMasks[i], i)+"; value="+(toValue(kmer, rkmer, kMasks[i]))+"; kmer="+kmer+"; rkmer="+rkmer+"; kmask="+kMasks[i]+"; rightMasks[i+1]="+rightMasks[i+1]);
						outstream.println("i="+i+"; tnum="+tnum+"; Looking for right-shift kmer "+AminoAcid.kmerToString(kmer&~kMasks[i], i));
						final long value=toValue(kmer, rkmer, kMasks[i]);
						if(map.contains(value)){outstream.println("Found "+value);}
					}
				}
			}
			return added;
		}
		
		
		/**
		 * Adds this kmer to the table, including any mutations implied by editDistance or hammingDistance.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param len Kmer length
		 * @param extraBase Base added to end in case of deletions
		 * @param id Scaffold number
		 * @param kmask0
		 * @return Number of kmers stored
		 */
		private long addToMap(final long kmer, final long rkmer, final int len, final long extraBase, final int id, final long kmask0, final int hdist, final int edist){
			
			assert(kmask0==kMasks[len]) : kmask0+", "+len+", "+kMasks[len]+", "+Long.numberOfTrailingZeros(kmask0)+", "+Long.numberOfTrailingZeros(kMasks[len]);
			
			if(verbose){outstream.println("addToMap_A; len="+len+"; kMasks[len]="+kMasks[len]);}
			assert((kmer&kmask0)==0);
			final long added;
			if(hdist==0){
				final long key=toValue(kmer, rkmer, kmask0);
				if(failsSpeed(key)){return 0;}
				if(key%WAYS!=tnum){return 0;}
				if(verbose){outstream.println("addToMap_B: "+AminoAcid.kmerToString(kmer&~kMasks[len], len)+" = "+key);}
				if(storeMode==SET_IF_NOT_PRESENT){
					added=map.setIfNotPresent(key, id);
				}else if(storeMode==SET_ALWAYS){
					added=map.set(key, id);
				}else{
					assert(storeMode==INCREMENT);
					added=map.increment(key, 1);
				}
			}else if(edist>0){
//				long extraBase=(i>=bases.length-1 ? -1 : AminoAcid.baseToNumber[bases[i+1]]);
				added=mutate(kmer, rkmer, len, id, edist, extraBase);
			}else{
				added=mutate(kmer, rkmer, len, id, hdist, -1);
			}
			if(verbose){outstream.println("addToMap added "+added+" keys.");}
			return added;
		}
		
		/**
		 * Mutate and store this kmer through 'dist' recursions.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param id Scaffold number
		 * @param dist Number of mutations
		 * @param extraBase Base added to end in case of deletions
		 * @return Number of kmers stored
		 */
		private long mutate(final long kmer, final long rkmer, final int len, final int id, final int dist, final long extraBase){
			long added=0;
			
			final long key=toValue(kmer, rkmer, kMasks[len]);
			
			if(verbose){outstream.println("mutate_A; len="+len+"; kmer="+kmer+"; rkmer="+rkmer+"; kMasks[len]="+kMasks[len]);}
			if(key%WAYS==tnum){
				if(verbose){outstream.println("mutate_B: "+AminoAcid.kmerToString(kmer&~kMasks[len], len)+" = "+key);}
				int x;
				if(storeMode==SET_IF_NOT_PRESENT){
					x=map.setIfNotPresent(key, id);
				}else if(storeMode==SET_ALWAYS){
					x=map.set(key, id);
				}else{
					assert(storeMode==INCREMENT);
					x=map.increment(key, 1);
					x=(x==1 ? 1 : 0);
				}
				if(verbose){outstream.println("mutate_B added "+x+" keys.");}
				added+=x;
				assert(map.contains(key));
			}
			
			if(dist>0){
				final int dist2=dist-1;
				
				//Sub
				for(int j=0; j<4; j++){
					for(int i=0; i<len; i++){
						final long temp=(kmer&clearMasks[i])|setMasks[j][i];
						if(temp!=kmer){
							long rtemp=AminoAcid.reverseComplementBinaryFast(temp, len);
							added+=mutate(temp, rtemp, len, id, dist2, extraBase);
						}
					}
				}
				
				if(editDistance>0){
					//Del
					if(extraBase>=0 && extraBase<=3){
						for(int i=1; i<len; i++){
							final long temp=(kmer&leftMasks[i])|((kmer<<2)&rightMasks[i])|extraBase;
							if(temp!=kmer){
								long rtemp=AminoAcid.reverseComplementBinaryFast(temp, len);
								added+=mutate(temp, rtemp, len, id, dist2, -1);
							}
						}
					}

					//Ins
					final long eb2=kmer&3;
					for(int i=1; i<len; i++){
						final long temp0=(kmer&leftMasks[i])|((kmer&rightMasks[i])>>2);
						for(int j=0; j<4; j++){
							final long temp=temp0|setMasks[j][i-1];
							if(temp!=kmer){
								long rtemp=AminoAcid.reverseComplementBinaryFast(temp, len);
								added+=mutate(temp, rtemp, len, id, dist2, eb2);
							}
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
		
		/** Destination for storing kmers */
		private final AbstractKmerTable map;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Transforms a kmer into a canonical value stored in the table.  Expected to be inlined.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param lengthMask Bitmask with single '1' set to left of kmer
	 * @return Canonical value
	 */
	private final long toValue(long kmer, long rkmer, long lengthMask){
		assert(lengthMask==0 || (kmer<lengthMask && rkmer<lengthMask)) : lengthMask+", "+kmer+", "+rkmer;
		long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		return (value&middleMask)|lengthMask;
	}
	
	final boolean passesSpeed(long key){
		return speed<1 || ((key&Long.MAX_VALUE)%17)>=speed;
	}
	
	final boolean failsSpeed(long key){
		return speed>0 && ((key&Long.MAX_VALUE)%17)<speed;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Has this class encountered errors while processing? */
	public boolean errorState=false;
	
	/** How to associate values with kmers */
	private int storeMode=SET_IF_NOT_PRESENT;
	
	/** Hold kmers.  A kmer X such that X%WAYS=Y will be stored in keySets[Y] */
	public AbstractKmerTable[] tables;
	/** A scaffold's name is stored at scaffoldNames.get(id).
	 * scaffoldNames[0] is reserved, so the first id is 1. */
	public ArrayList<String> scaffoldNames;
	/** Names of reference files (refNames[0] is valid). */
	public ArrayList<String> refNames;
	/** Number of scaffolds per reference. */
	public int[] refScafCounts;
	/** scaffoldLengths[id] stores the length of that scaffold */
	public IntList scaffoldLengths=new IntList();
	
	/** Make the middle base in a kmer a wildcard to improve sensitivity */
	public final boolean maskMiddle=false;
	/** Correct errors via read overlap */
	public boolean ecc=false;
	
	/** Store reference kmers with up to this many substitutions */
	public final int hammingDistance;
	/** Store reference kmers with up to this many edits (including indels) */
	public final int editDistance;
	/** Store short reference kmers with up to this many substitutions */
	public int hammingDistance2=-1;
	/** Store short reference kmers with up to this many edits (including indels) */
	public int editDistance2=-1;
	/** Always skip at least this many consecutive kmers when hashing reference.
	 * 1 means every kmer is used, 2 means every other, etc. */
	private int minRefSkip=0;
	/** Never skip more than this many consecutive kmers when hashing reference. */
	private int maxRefSkip=0;
	
	private boolean variableRefSkip=false;
	
	/*--------------------------------------------------------------*/
	/*----------------          Statistics          ----------------*/
	/*--------------------------------------------------------------*/
	
	long refReads=0;
	long refBases=0;
	long refKmers=0;
	
	long storedKmers=0;
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Look for reverse-complements as well as forward kmers.  Default: true */
	private final boolean rcomp;
	/** AND bitmask with 0's at the middle base */
	private final long middleMask;
	
	/** Normal kmer length */
	private final int k;
	/** k-1; used in some expressions */
	private final int k2;
	/** Shortest kmer to use for trimming */
	private final int mink;
	/** Attempt to match kmers shorter than normal k on read ends when doing kTrimming. */
	private final boolean useShortKmers;
	
	/** Fraction of kmers to skip, 0 to 16 out of 17 */
	private final int speed;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Default initial size of data structures */
	private static final int initialSizeDefault=128000;
	
	/** Number of tables (and threads, during loading) */
	private static final int WAYS=7; //123
	/** Verbose messages */
	public static final boolean verbose=false; //123
	
	/** Print messages to this stream */
	private static PrintStream outstream=System.err;
	/** Display progress messages such as memory usage */
	public static boolean DISPLAY_PROGRESS=true;
	/** Indicates end of input stream */
	private static final ArrayList<Read> POISON=new ArrayList<Read>(0);
	/** Make unambiguous copies of ref sequences with ambiguous bases */
	public static boolean REPLICATE_AMBIGUOUS=false;
	
	/** x&clearMasks[i] will clear base i */
	private static final long[] clearMasks;
	/** x|setMasks[i][j] will set base i to j */
	private static final long[][] setMasks;
	/** x&leftMasks[i] will clear all bases to the right of i (exclusive) */
	private static final long[] leftMasks;
	/** x&rightMasks[i] will clear all bases to the left of i (inclusive) */
	private static final long[] rightMasks;
	/** x|kMasks[i] will set the bit to the left of the leftmost base */
	private static final long[] kMasks;
	
	public static final int SET_IF_NOT_PRESENT=1, SET_ALWAYS=2, INCREMENT=3;
	
	/*--------------------------------------------------------------*/
	/*----------------      Static Initializers     ----------------*/
	/*--------------------------------------------------------------*/
	
	static{
		clearMasks=new long[32];
		leftMasks=new long[32];
		rightMasks=new long[32];
		kMasks=new long[32];
		setMasks=new long[4][32];
		for(int i=0; i<32; i++){
			clearMasks[i]=~(3L<<(2*i));
		}
		for(int i=0; i<32; i++){
			leftMasks[i]=((-1L)<<(2*i));
		}
		for(int i=0; i<32; i++){
			rightMasks[i]=~((-1L)<<(2*i));
		}
		for(int i=0; i<32; i++){
			kMasks[i]=((1L)<<(2*i));
		}
		for(int i=0; i<32; i++){
			for(long j=0; j<4; j++){
				setMasks[(int)j][i]=(j<<(2*i));
			}
		}
	}
	
}
