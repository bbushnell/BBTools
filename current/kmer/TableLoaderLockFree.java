package kmer;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.TimeUnit;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import jgi.BBMerge;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Read;
import structures.IntList;
import structures.ListNum;

/**
 * Loads length-tagged short DNA k-mers into seven caller-supplied table shards.
 * Every worker scans every input batch but writes only keys assigned to its shard
 * by {@code key%WAYS}; this ownership avoids table locks during loading. It does
 * not make the supplied tables safe for concurrent external access. Callers keep
 * tables, configuration and queued reads stable until loading finishes.
 * Optional mutation expands reference keys; optional short keys come only from
 * complete terminal k-mer windows. Each call reports its own input statistics;
 * retaining scaffold metadata requires initially empty destination tables.
 * The mutation path retains keys without exact-key subsampling.
 *
 * @author Brian Bushnell
 * @date Mar 4, 2015
 */
public class TableLoaderLockFree{

	/**
	 * Loads the remaining command-line arguments as reference filenames using k=31.
	 * @param args PreParser options followed by reference paths
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

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Uses exact canonical keys without short-key generation or subsampling.
	 * @param tables_ Seven distinct table shards, borrowed and not cleared
	 * @param k_ Full k-mer length, 1 through 31
	 */
	public TableLoaderLockFree(AbstractKmerTable[] tables_, int k_){
		this(tables_, k_, 0, 0, 0, 0, true, false);
	}

	/**
	 * Borrows destination tables and records key-generation settings.
	 * Does not validate every setting or clear existing keys. Short-key mutation
	 * settings remain separately configurable through hammingDistance2/editDistance2.
	 *
	 * @param tables_ Seven distinct table shards, exclusively mutated by this loader
	 * @param k_ Full k-mer length, 1 through 31
	 * @param mink_ Enables terminal keys of lengths mink_ through k_-1 when {@code 0<mink_<k_}
	 * @param speed_ Exact-path remainder cutoff modulo 17; nonpositive disables sampling,
	 * 1..16 select subsets, and values at least 17 omit every exact-path key
	 * @param hdist_ Substitution recursion budget when edist_ is nonpositive
	 * @param edist_ Positive edit recursion budget overrides hdist_ and enables indels
	 * @param rcomp_ Select the larger forward/reverse-complement encoding when true
	 * @param maskMiddle_ Mask the middle base of full keys; short-key generation disables this
	 */
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
		//Keep the same short-end precedence and key encoding as TableReader/BBDukParser.
		maskMiddle=maskMiddle_ && !useShortKmers;
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

	/**
	 * Creates an array of k-mer tables with automatic size scheduling.
	 *
	 * @param tableType Type of table implementation to create
	 * @param bytesPerKmer Memory budget per k-mer
	 * @param coreMask Table-specific key-placement mask, not CPU affinity
	 * @param prealloc Whether to preallocate table memory
	 * @param memRatio Memory usage ratio for sizing
	 * @return Seven initialized table shards using the selected implementation
	 */
	public static AbstractKmerTable[] makeTables(int tableType, int bytesPerKmer, long coreMask,
			boolean prealloc, double memRatio){
		ScheduleMaker scheduleMaker=new ScheduleMaker(WAYS, bytesPerKmer, prealloc, memRatio);
		int[] schedule=scheduleMaker.makeSchedule();
		return AbstractKmerTable.preallocate(WAYS, tableType, schedule, coreMask);
	}

	/** Retained per-instance sizing helper; supplied destination tables are not allocated from it. */
	ScheduleMaker scheduleMaker=new ScheduleMaker(WAYS, 12, false, 0.8);
	/** Retained eager schedule, unused by this loader's loading path. */
	int[] schedule=scheduleMaker.makeSchedule();

	/**
	 * Replaces optional metadata, then adds reference keys to the existing tables.
	 * Processes reference files in order followed by literals. Input counters reset
	 * per call; existing keys and accumulated errorState are retained. Metadata IDs
	 * restart at one, so keeping names requires empty primary and victim storage.
	 * Callers may append without names or explicitly empty tables before retaining
	 * a fresh set of names. Producer/worker failures cancel loading and propagate after all
	 * started loader workers terminate; already inserted keys are not rolled back.
	 *
	 * @param ref Array of reference file paths to process
	 * @param literal Array of literal sequences to process as references
	 * @param keepNames Retain scaffold/reference metadata, reserving scaffold ID zero
	 * @param useRefNames Replace scaffold names with stripped reference names after loading
	 * @param ecc_ Request paired file-read overlap correction independently of name retention
	 * @return Number of distinct encoded keys newly stored during this invocation
	 * @throws IllegalStateException If keepNames is requested with nonempty destination tables
	 */
	public long processData(String[] ref, String[] literal, boolean keepNames, boolean useRefNames, boolean ecc_){

		if(keepNames){
			//KMR121: a new ID-to-name mapping must not reinterpret keys from an older load.
			for(int i=0; i<tables.length; i++){
				final AbstractKmerTable table=tables[i];
				if(table.size()!=0 || (table instanceof HashArray && ((HashArray)table).victims().size()!=0)){
					throw new IllegalStateException("Retaining scaffold names requires empty destination tables; shard="+i);
				}
			}
		}
		//KMR121: worker totals are divided by WAYS once for this invocation only.
		refReads=refBases=refKmers=0;
		scaffoldNames=null;
		refNames=null;
		refScafCounts=null;
		scaffoldLengths=null;
		ecc=ecc_;

		if(keepNames){
			scaffoldNames=new ArrayList<String>();
			refNames=new ArrayList<String>();
			scaffoldLengths=new IntList();
			//KMR115: real scaffold IDs start at one; lookup consumers reserve zero for no match.
			scaffoldNames.add("");
			scaffoldLengths.add(0);

			if(ref!=null){
				for(String s : ref){refNames.add(s);}
			}
			if(literal!=null){refNames.add("literal");}
			refScafCounts=new int[refNames.size()];
		}

		final long added=spawnLoadThreads(ref, literal);
		//KMR115: reference groups are known only after files and literals have been processed.
		if(keepNames && useRefNames){toRefNames();}
		return added;
	}

	/** Sets a fixed sampling interval measured within consecutive fully defined bases.
	 * @param x Values at most one inspect every full window; larger values require valid-run length divisible by x */
	public void setRefSkip(int x){setRefSkip(x, x);}

	/**
	 * Sets the bounds for a length-dependent sampling interval when they differ.
	 * The length-based choice (0, 2, 11 or k) is clamped to these bounds.
	 * @param min Lower interval bound; values at most one sample every full window
	 * @param max Upper interval bound, raised to min if smaller
	 */
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

	/** Sets the storage mode for k-mer insertion into tables.
	 * @param x Storage mode: SET_IF_NOT_PRESENT, SET_ALWAYS, or INCREMENT */
	public void storeMode(final int x){
		assert(x==SET_IF_NOT_PRESENT || x==SET_ALWAYS || x==INCREMENT);
		storeMode=x;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Broadcasts each shared read batch to seven shard-owning workers and joins them.
	 * Workers scan the same data; read/base/full-window totals are divided by seven.
	 * Caller processData resets counters first. Stream-close errors accumulate in errorState.
	 * Failures cancel pending batches and are rethrown after joining every started
	 * worker. Workers may finish their current batch after cancellation; no hard
	 * shutdown-latency bound or forced stop of a table operation is provided.
	 * A blocked input read/close still depends on the concrete Streamer to return.
	 * @return Number of newly stored keys, or zero for no inputs
	 */
	private long spawnLoadThreads(String[] ref, String[] literal){
		Timer t=new Timer();
		if((ref==null || ref.length<1) && (literal==null || literal.length<1)){return 0;}
		long added=0;

		final LoadThread[] loaders=new LoadThread[WAYS];
		final LoadState state=new LoadState(loaders);
		boolean interrupted=false;
		try{
			/* Construct every worker before publishing any through Thread.start(). */
			for(int i=0; i<loaders.length; i++){
				loaders[i]=new LoadThread(i, state);
			}
			for(LoadThread lt : loaders){state.checkFailure(); lt.start();}

			/* For each reference file... */
			int refNum=0;
			if(ref!=null){
				for(String refname : ref){
					state.checkFailure();

					/* Start an input stream */
					FileFormat ff=FileFormat.testInput(refname, FileFormat.FASTA, null, false, true);
					Streamer cris=StreamerFactory.getReadInputStream(-1, false, ff, null, -1);
					try{
						cris.start();
						ListNum<Read> ln=cris.nextList();
						ArrayList<Read> reads=(ln!=null ? ln.list : null);

						/* Iterate through read lists from the input stream */
						while(ln!=null && reads!=null){//ln!=null prevents a compiler potential null access warning
							state.checkFailure();
							{
								/* Correct requested pairs, then optionally assign scaffold metadata. */
								ArrayList<Read> reads2=new ArrayList<Read>(reads);
								if(ecc || scaffoldNames!=null){
									for(Read r1 : reads2){
										final Read r2=r1.mate;
										if(ecc && r2!=null){BBMerge.findOverlapStrict(r1, r2, true);}
										if(scaffoldNames!=null){
											final Integer id=scaffoldNames.size();
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
								}

								if(REPLICATE_AMBIGUOUS){
									reads2=Tools.replicateAmbiguous(reads2, Tools.min(k, mink));
								}

								/* Send a pointer to the read list to each LoadThread */
								state.broadcast(reads2);
							}
							ln=cris.nextList();
							reads=(ln!=null ? ln.list : null);
						}
					}catch(Throwable failure){
						state.fail(failure);
						throw failure;
					}finally{
						try{errorState|=ReadWrite.closeStream(cris);}
						catch(Throwable failure){state.fail(failure);}
					}
					state.checkFailure();
					refNum++;
				}
			}

			/* If there are literal sequences to use as references */
			if(literal!=null){
				state.checkFailure();
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
						//KMR114: retaining metadata must not omit the literal from the worker batch.
						list.add(r);
					}else{
						final Read r=new Read(literal[i].getBytes(), null, i);
						list.add(r);
					}
				}

				if(REPLICATE_AMBIGUOUS){
					list=Tools.replicateAmbiguous(list, Tools.min(k, mink));
				}

				/* Send a pointer to the read list to each LoadThread */
				state.broadcast(list);
			}

			/* Signal loaders to terminate */
			state.broadcast(POISON);
		}catch(Throwable failure){
			interrupted=(failure instanceof InterruptedException);
			state.fail(failure);
		}finally{
			if(Thread.interrupted()){
				interrupted=true;
				state.fail(new InterruptedException("Reference loading was interrupted."));
			}
			for(LoadThread lt : loaders){
				if(lt==null){continue;}
				while(lt.isAlive()){
					try{lt.join();}
					catch(InterruptedException failure){
						interrupted=true;
						state.fail(failure);
					}
				}
			}
			if(interrupted){Thread.currentThread().interrupt();}
			errorState|=state.failed();
		}
		state.checkFailure();

		/* Merge only successful workers after every join. */
		for(LoadThread lt : loaders){
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
	 * Replaces scaffold names with stripped reference filenames using refScafCounts.
	 * Requires populated metadata; the reserved scaffold slot zero remains unchanged.
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

	/** Owns cancellation for one invocation; no synchronization occurs per k-mer. */
	private static final class LoadState{

		/** Borrows the worker array, populated completely before workers start. */
		LoadState(final LoadThread[] loaders_){loaders=loaders_;}

		/** Publishes the first failure before waking workers blocked on their queues. */
		void fail(final Throwable cause){
			synchronized(this){
				if(failure!=null){return;}
				failure=cause;
			}
			for(LoadThread lt : loaders){
				if(lt!=null && lt!=Thread.currentThread() && lt.isAlive()){lt.interrupt();}
			}
		}

		/** Returns whether this invocation has failed, including producer cancellation. */
		boolean failed(){return failure!=null;}

		/** Preserves unchecked failure types; checked interruptions retain their cause. */
		void checkFailure(){
			final Throwable cause=failure;
			if(cause==null){return;}
			if(cause instanceof RuntimeException){throw (RuntimeException)cause;}
			if(cause instanceof Error){throw (Error)cause;}
			throw new RuntimeException("Reference loading failed.", cause);
		}

		/**
		 * Shares a batch with each worker, checking failure during backpressure.
		 * Timed offers prevent a failed worker's full queue from trapping the producer.
		 * The timeout is only a cancellation poll; a healthy slow worker is not rejected.
		 */
		void broadcast(final ArrayList<Read> reads) throws InterruptedException{
			for(LoadThread lt : loaders){
				checkFailure();
				while(!lt.queue.offer(reads, 100, TimeUnit.MILLISECONDS)){checkFailure();}
			}
			checkFailure();
		}

		/** All workers belonging to this invocation, including any not yet started. */
		private final LoadThread[] loaders;
		/** First observed failure; volatile publication makes batch-boundary checks safe. */
		private volatile Throwable failure;
	}

	/**
	 * Scans every queued read but mutates only table tnum for keys with key%WAYS==tnum.
	 * The producer owns read construction/correction before broadcasting; workers
	 * share the read objects and must not change them while scanning.
	 */
	private class LoadThread extends Thread{

		/** Constructs a LoadThread for a specific thread number.
		 * @param tnum_ Thread number determining which k-mers this thread processes
		 * @param state_ Shared cancellation and first failure for this loading invocation */
		public LoadThread(final int tnum_, final LoadState state_){
			tnum=tnum_;
			map=tables[tnum];
			state=state_;
		}

		/**
		 * Waits for the next batch unless a peer or producer has cancelled this invocation.
		 * @return Shared reads or the identity-tested POISON list
		 * @throws InterruptedException If interrupted while waiting
		 */
		private ArrayList<Read> fetch() throws InterruptedException{
			return state.failed() ? POISON : queue.take();
		}

		/** Consumes batches through POISON, then optionally rebalances its exclusively owned shard. */
		@Override
		public void run(){
			try{processBatches();}
			catch(Throwable failure){state.fail(failure);}
		}

		/** Processes accepted batches; cancellation may allow the current batch to finish. */
		private void processBatches() throws InterruptedException{
			ArrayList<Read> reads=fetch();
			while(reads!=POISON && !state.failed()){
				for(Read r1 : reads){
					assert(r1.pairnum()==0);
					final Read r2=r1.mate;

					addedT+=addToMap(r1, minRefSkip);
					if(r2!=null){addedT+=addToMap(r2, minRefSkip);}
				}
				reads=fetch();
			}

			if(!state.failed() && map.canRebalance() && map.size()>2L*map.arrayLength()){
				map.rebalance();
			}
		}

		/**
		 * Counts nonnull sequence reads/bases and fully defined full-length windows.
		 * Skipping controls insertion only, not those statistics. Short terminal keys
		 * are derived only from selected full windows; a read shorter than k yields none.
		 * @param r The current read to process
		 * @param skip Sampling interval, possibly replaced by the length-dependent setting
		 * @return Number of new keys stored in this shard while processing the read
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
					//Comprehension (verified, NOT a bug): on an N (x<0) only len and rkmer are reset, NOT kmer — the pervasive
					//BBTools idiom (see BBDuk2, QualityTools). Safe because the forward kmer is built by LEFT-shift, so any
					//N-garbage is shifted out of the 2k-bit mask window after k valid bases; and len>=k gates all use until
					//then. rkmer (right-shift, new base at the high end) is reset defensively so it's clean the moment len>=k.
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
		 * Adds shorter suffixes of the final full k-mer at the read's right end.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param extraBase Base added to end in case of deletions
		 * @param id Scaffold number
		 * @return New keys stored for lengths k-1 through mink
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
		 * Adds shorter prefixes of the first full k-mer at the read's left end.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param id Scaffold number
		 * @return New keys stored for lengths k-1 through mink
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
		 * Dispatches exact insertion or mutation expansion with invocation-local edit eligibility.
		 * hdist==0 uses the exact path only when edist is nonpositive; positive edist
		 * chooses edit recursion. Only exact insertion applies the speed filter.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param len Kmer length
		 * @param extraBase Base added to end in case of deletions
		 * @param id Scaffold number
		 * @param kmask0 Length marker, equal to kMasks[len]
		 * @param hdist Substitution budget when edist is nonpositive; zero selects exact insertion
		 * @param edist Positive budget selects substitution/indel recursion for this invocation
		 * @return New keys stored by exact insertion or mutation expansion
		 */
		private long addToMap(final long kmer, final long rkmer, final int len, final long extraBase, final int id, final long kmask0, final int hdist, final int edist){

			assert(kmask0==kMasks[len]) : kmask0+", "+len+", "+kMasks[len]+", "+Long.numberOfTrailingZeros(kmask0)+", "+Long.numberOfTrailingZeros(kMasks[len]);

			if(verbose){outstream.println("addToMap_A; len="+len+"; kMasks[len]="+kMasks[len]);}
			assert((kmer&kmask0)==0);
			final long added;
			if(hdist==0 && edist<=0){
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
					//KMR124 [kmer/TableLoaderLockFree#001]: report new keys, not their updated counts.
					added=map.incrementAndReturnNumCreated(key, 1);
				}
			}else if(edist>0){
//				long extraBase=(i>=bases.length-1 ? -1 : AminoAcid.baseToNumber[bases[i+1]]);
				added=mutate(kmer, rkmer, len, id, edist, extraBase, true);
			}else{
				added=mutate(kmer, rkmer, len, id, hdist, -1, false);
			}
			if(verbose){outstream.println("addToMap added "+added+" keys.");}
			return added;
		}

		/**
		 * Stores this key and recursively visits substitutions and eligible indels.
		 * May revisit a key through several mutation paths; storage mode determines
		 * whether repeated visits change its value. Indel eligibility belongs to this
		 * invocation and is preserved through recursion. No speed filtering occurs here.
		 * @param kmer Forward kmer
		 * @param rkmer Reverse kmer
		 * @param len Number of encoded bases
		 * @param id Scaffold number
		 * @param dist Remaining recursion budget; nonpositive values store only this key
		 * @param extraBase Base added to end in case of deletions
		 * @param allowIndels Whether this traversal may generate insertions and deletions
		 * @return New keys stored by this traversal; existing keys may have their counts updated
		 */
		private long mutate(final long kmer, final long rkmer, final int len, final int id, final int dist, final long extraBase, final boolean allowIndels){
			long added=0;

			//Mutation keys form a storage superset; TableReader still applies the shared query selector.
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
					//KMR124: count==1 does not prove creation when an existing key held zero.
					x=map.incrementAndReturnNumCreated(key, 1);
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
							added+=mutate(temp, rtemp, len, id, dist2, extraBase, allowIndels);
						}
					}
				}

				if(allowIndels){
					//Del
					if(extraBase>=0 && extraBase<=3){
						for(int i=1; i<len; i++){
							final long temp=(kmer&leftMasks[i])|((kmer<<2)&rightMasks[i])|extraBase;
							if(temp!=kmer){
								long rtemp=AminoAcid.reverseComplementBinaryFast(temp, len);
								added+=mutate(temp, rtemp, len, id, dist2, -1, allowIndels);
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
									added+=mutate(temp, rtemp, len, id, dist2, eb2, allowIndels);
							}
						}
					}
				}

			}

			return added;
		}

		/*--------------------------------------------------------------*/

		/** Number of new keys stored by this worker; existing-key count updates do not contribute. */
		public long addedT=0;
		/** Full valid windows, nonnull-sequence reads and bases scanned by this worker. */
		public long refKmersT=0, refReadsT=0, refBasesT=0;
		/** Thread number; used to determine which kmers to store */
		public final int tnum;
		/** Shared batches plus an identity-tested terminal marker; producer blocks when full. */
		public final ArrayBlockingQueue<ArrayList<Read>> queue=new ArrayBlockingQueue<ArrayList<Read>>(32);

		/** Destination for storing kmers */
		private final AbstractKmerTable map;
		/** Invocation-local cancellation and first failure, shared with the producer. */
		private final LoadState state;

	}

	/*--------------------------------------------------------------*/
	/*----------------         Key Helpers          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Adds the length marker to a forward or canonical encoded key.
	 * Canonicalization selects the larger orientation only when rcomp is enabled.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @param lengthMask Bitmask with single '1' set to left of kmer
	 * @return Selected orientation AND middleMask, then OR the length marker
	 */
	private final long toValue(long kmer, long rkmer, long lengthMask){
		assert(lengthMask==0 || (kmer<lengthMask && rkmer<lengthMask)) : lengthMask+", "+kmer+", "+rkmer;
		long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
		return (value&middleMask)|lengthMask;
	}

	/**
	 * Determines if k-mer passes speed filter for storage.
	 * @param key K-mer value to test
	 * @return true if k-mer should be stored based on speed setting
	 */
	final boolean passesSpeed(long key){
		return speed<1 || ((key&Long.MAX_VALUE)%17)>=speed;
	}

	/**
	 * Determines if k-mer fails speed filter and should be skipped.
	 * @param key K-mer value to test
	 * @return true if k-mer should be skipped based on speed setting
	 */
	final boolean failsSpeed(long key){
		return speed>0 && ((key&Long.MAX_VALUE)%17)<speed;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Accumulated input-close, producer and worker error flag; finalized after worker joins. */
	public boolean errorState=false;

	/** How to associate values with kmers */
	private int storeMode=SET_IF_NOT_PRESENT;

	/** Borrowed seven-shard destination; key X belongs to tables[X%WAYS]. */
	public AbstractKmerTable[] tables;
	/** Scaffold names indexed by assigned ID; slot zero is reserved and real IDs start at one. */
	public ArrayList<String> scaffoldNames;
	/** Names of reference files (refNames[0] is valid). */
	public ArrayList<String> refNames;
	/** Root-read counts per reference, with a final literal entry when literals are supplied. */
	public int[] refScafCounts;
	/** Bases per assigned root ID, including its mate when present; null when names are disabled. */
	public IntList scaffoldLengths=new IntList();

	/** Effective middle-base masking flag; short-key generation takes precedence. */
	public final boolean maskMiddle;
	/** Requested paired file-read overlap correction, independent of optional name retention. */
	public boolean ecc=false;

	/** Full-key substitution budget when the edit budget is nonpositive. */
	public final int hammingDistance;
	/** Positive full-key edit budget overrides substitutions-only traversal and permits indels. */
	public final int editDistance;
	/** Short-key substitution budget; negative default performs no further substitution recursion. */
	public int hammingDistance2=-1;
	/** Positive short-key edit budget overrides hammingDistance2 and permits short-key indels. */
	public int editDistance2=-1;
	/** Minimum sampling interval: values at most one use every full window, two every other valid-run position. */
	private int minRefSkip=0;
	/** Maximum interval for length-dependent sampling; fixed sampling uses minRefSkip. */
	private int maxRefSkip=0;

	/** Whether skip rate varies based on sequence length */
	private boolean variableRefSkip=false;

	/*--------------------------------------------------------------*/
	/*----------------          Statistics          ----------------*/
	/*--------------------------------------------------------------*/

	/** Current invocation's scanned reads with nonnull bases, including mates. */
	long refReads=0;
	/** Current invocation's scanned bases, before window selection. */
	long refBases=0;
	/** Current invocation's fully defined full-length windows before sampling/mutation; excludes short keys. */
	long refKmers=0;

	/** Retained unused statistic; loading returns its creation count without updating this field. */
	long storedKmers=0;

	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/

	/** Select the larger forward/reverse-complement encoding rather than storing both. */
	private final boolean rcomp;
	/** Clears the configured middle base; all ones when middle masking is disabled. */
	private final long middleMask;

	/** Normal kmer length */
	private final int k;
	/** k-1; used in some expressions */
	private final int k2;
	/** Shortest terminal key length when short-key generation is enabled. */
	private final int mink;
	/** Generate terminal keys shorter than k from selected complete end windows. */
	private final boolean useShortKmers;

	/** Exact-path cutoff on key remainder modulo 17; the mutation path retains keys without this filter. */
	private final int speed;

	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/

	/** Retained unused legacy initial-size constant. */
	private static final int initialSizeDefault=128000;

	/** Number of tables (and threads, during loading) */
	private static final int WAYS=7; //123
	/** Verbose messages */
	public static final boolean verbose=false; //123

	/** Print messages to this stream */
	private static PrintStream outstream=System.err;
	/** Display progress messages such as memory usage */
	public static boolean DISPLAY_PROGRESS=true;
	/** Terminal batch distinguished by identity, not by empty-list contents. */
	private static final ArrayList<Read> POISON=new ArrayList<Read>(0);
	/** Make unambiguous copies of ref sequences with ambiguous bases */
	public static boolean REPLICATE_AMBIGUOUS=false;

	/** {@code x&clearMasks[i]} clears base i. */
	private static final long[] clearMasks;
	/** {@code x|setMasks[base][position]} installs a base after its two bits have been cleared. */
	private static final long[][] setMasks;
	/** {@code x&leftMasks[i]} clears all bases to the right of i (exclusive). */
	private static final long[] leftMasks;
	/** {@code x&rightMasks[i]} clears all bases to the left of i (inclusive). */
	private static final long[] rightMasks;
	/** x|kMasks[i] will set the bit to the left of the leftmost base */
	private static final long[] kMasks;

	/** Insert only absent keys, delegate scalar set, or increment counts, respectively. */
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
