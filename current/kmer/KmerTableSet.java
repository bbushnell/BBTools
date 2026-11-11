package kmer;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicLong;

import assemble.Contig;
import bloom.KmerCountAbstract;
import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Primes;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.Streamer;
import stream.StreamerFactory;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;
import structures.LongList;
import tracker.ReadStats;
import ukmer.Kmer;

/**
 * Loads sharded short-key count tables for assembly and k-mer counting. DNA loading
 * uses two-bit keys; amino loading uses five-bit symbols with k capped at twelve.
 * Construction parses configuration and validates filenames but does not allocate
 * shards or load reads. Processing/allocation lifecycle is inherited from
 * AbstractKmerTableSet. Helpers for sequence ends, extension counts and ownership
 * are DNA-specific; several accessors require HashArray1D rather than arbitrary table types.
 * Configuration mutates shared parser/quality/output settings. Complete loading
 * before querying or claiming ownership, and exclude structural changes during
 * traversal. Worker error propagation remains incomplete (KMR011).
 * @author Brian Bushnell
 * @date Jun 22, 2015
 *
 */
public class KmerTableSet extends AbstractKmerTableSet{

	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer(), t2=new Timer();
		t.start();
		t2.start();
		KmerTableSet set=new KmerTableSet(args);
		t2.stop();
		outstream.println("Initialization Time:      \t"+t2);

		///And run it
		set.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}

	/**
	 * Creates the command-line loader using the default twelve-byte capacity estimate.
	 * @param args Command line arguments
	 */
	private KmerTableSet(String[] args){
		this(args, 12);//+5 if using ownership and building contigs
	}

	/**
	 * Parses options, computes budget estimates and validates input paths without loading.
	 * Expands a missing filename containing # into mate names; an existing literal
	 * filename is retained. Extra mate names become separate extra inputs.
	 * @param args Command line arguments, including shared parser options
	 * @param bytesPerKmer_ Positive bytes-per-key estimate used for budgeting, not a measured allocation
	 */
	public KmerTableSet(String[] args, final int bytesPerKmer_){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}//TODO - no easy way to close outstream

		/* Initialize local variables with defaults */
		Parser parser=new Parser();
		boolean prealloc_=false;
		int k_=31;
		int ways_=-1;
		int filterMax_=2;
		boolean ecco_=false, merge_=false;
		boolean rcomp_=true;
		double minProb_=defaultMinprob;
		int tableType_=AbstractKmerTable.ARRAY1D;
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			// Preserve '=' inside values; a missing or terminal delimiter still gives b=null.
			String[] split=Parse.splitOnFirst(arg, '=');
			String a=split[0].toLowerCase(Locale.ROOT);
			String b=split.length>1 ? split[1] : null;

			if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(parser.parseTrim(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("in1")){
				in1.clear();
				if(b!=null){
					String[] s=b.split(",");
					for(String ss : s){
						in1.add(ss);
					}
				}
			}else if(a.equals("in2")){
				in2.clear();
				if(b!=null){
					String[] s=b.split(",");
					for(String ss : s){
						in2.add(ss);
					}
				}
			}else if(a.equals("extra")){
				extra.clear();
				if(b!=null){
					String[] s=b.split(",");
					for(String ss : s){
						extra.add(ss);
					}
				}
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("initialsize")){
				initialSize=Parse.parseIntKMG(b);
			}else if(a.equals("showstats") || a.equals("stats")){
				showStats=Parse.parseBoolean(b);
			}else if(a.equals("ways")){
				ways_=Parse.parseIntKMG(b);
			}else if(a.equals("buflen") || a.equals("bufflen") || a.equals("bufferlength")){
				buflen=Parse.parseIntKMG(b);
			}else if(a.equals("k")){
				assert(b!=null) : "\nk needs an integer value from 1 to 31, such as k=27.  Default is 31.\n";
				k_=Parse.parseIntKMG(b);
			}else if(a.equals("threads") || a.equals("t")){
				Shared.setThreads(b);
			}else if(a.equals("showspeed") || a.equals("ss")){
				showSpeed=Parse.parseBoolean(b);
			}else if(a.equals("ecco")){
				ecco_=Parse.parseBoolean(b);
			}else if(a.equals("merge")){
				merge_=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
//				assert(false) : "Verbose flag is currently static final; must be recompiled to change.";
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("verbose2")){
//				assert(false) : "Verbose flag is currently static final; must be recompiled to change.";
				verbose2=Parse.parseBoolean(b);
			}else if(a.equals("minprob")){
				minProb_=Double.parseDouble(b);
			}else if(a.equals("minprobprefilter") || a.equals("mpp")){
				minProbPrefilter=Parse.parseBoolean(b);
			}else if(a.equals("minprobmain") || a.equals("mpm")){
				minProbMain=Parse.parseBoolean(b);
			}else if(a.equals("reads") || a.startsWith("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("prealloc") || a.equals("preallocate")){
				if(b==null || b.length()<1 || Character.isLetter(b.charAt(0))){
					prealloc_=Parse.parseBoolean(b);
				}else{
					preallocFraction=Tools.max(0, Double.parseDouble(b));
					prealloc_=(preallocFraction>0);
				}
			}else if(a.equals("prefilter")){
				if(b==null || b.length()<1 || !Tools.isDigit(b.charAt(0))){
					prefilter=Parse.parseBoolean(b);
				}else{
					filterMax_=Parse.parseIntKMG(b);
					prefilter=filterMax_>0;
				}
			}else if(a.equals("prefiltersize") || a.equals("prefilterfraction") || a.equals("pff")){
				prefilterFraction=Tools.max(0, Double.parseDouble(b));
				assert(prefilterFraction<=1) : "prefiltersize must be 0-1, a fraction of total memory.";
				prefilter=prefilterFraction>0;
			}else if(a.equals("prehashes") || a.equals("hashes")){
				prehashes=Parse.parseIntKMG(b);
			}else if(a.equals("prefilterpasses") || a.equals("prepasses")){
				assert(b!=null) : "Bad parameter: "+arg;
				if(b.equalsIgnoreCase("auto")){
					prepasses=-1;
				}else{
					prepasses=Parse.parseIntKMG(b);
				}
			}else if(a.equals("onepass")){
				onePass=Parse.parseBoolean(b);
			}else if(a.equals("passes")){
				int passes=Parse.parseIntKMG(b);
				onePass=(passes<2);
			}else if(a.equals("rcomp")){
				rcomp_=Parse.parseBoolean(b);
			}else if(a.equals("tabletype")){
				tableType_=Integer.parseInt(b);
			}

			else if(a.equalsIgnoreCase("filterMemoryOverride") || a.equalsIgnoreCase("filterMemory") ||
					a.equalsIgnoreCase("prefilterMemory") || a.equalsIgnoreCase("filtermem")){
				filterMemoryOverride=Parse.parseKMG(b);
			}

			else if(IGNORE_UNKNOWN_ARGS){
				//Do nothing
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		{//Process parser fields
			Parser.processQuality();

			qtrimLeft=parser.qtrimLeft;
			qtrimRight=parser.qtrimRight;
			trimq=parser.trimq;
			trimE=parser.trimE();

			minAvgQuality=parser.minAvgQuality;
			minAvgQualityBases=parser.minAvgQualityBases;
			amino=Shared.AMINO_IN;
			if(amino){k_=Tools.min(k_, 12);}
		}

		if(prepasses==0 || !prefilter){
			prepasses=0;
			prefilter=false;
		}

//		assert(false) : prepasses+", "+onePass+", "+prefilter;

		{
			long memory=Runtime.getRuntime().maxMemory();
			double xmsRatio=Shared.xmsRatio();
//			long tmemory=Runtime.getRuntime().totalMemory();
			usableMemory=(long)Tools.max(((memory-96000000)*(xmsRatio>0.97 ? 0.82 : 0.72)), memory*0.45);
			if(prepasses==0 || !prefilter){
				filterMemory0=filterMemory1=0;
			}else if(filterMemoryOverride>0){
				filterMemory0=filterMemory1=filterMemoryOverride;
			}else{
				double low=Tools.min(prefilterFraction, 1-prefilterFraction);
				double high=1-low;
				if(prepasses<0 || (prepasses&1)==1){//odd passes
					filterMemory0=(long)(usableMemory*low);
					filterMemory1=(long)(usableMemory*high);
				}else{//even passes
					filterMemory0=(long)(usableMemory*high);
					filterMemory1=(long)(usableMemory*low);
				}
			}
			tableMemory=(long)(usableMemory*.95-filterMemory0);
		}

		tableType=tableType_;
		prealloc=prealloc_;
		bytesPerKmer=bytesPerKmer_;
		if(ways_<1){
			long maxKmers=(2*tableMemory)/bytesPerKmer;
			long minWays=Tools.min(10000, maxKmers/Integer.MAX_VALUE);
			ways_=(int)Tools.max(31, (int)(Tools.min(96, Shared.threads())*2.5), minWays);
			ways_=(int)Primes.primeAtLeast(ways_);
			assert(ways_>0);
//			System.err.println("ways="+ways_);
		}

		/* Set final variables; post-process and validate argument combinations */

		onePass=onePass&prefilter;
		ways=ways_;
		filterMax=Tools.min(filterMax_, 0x7FFFFFFF);
		ecco=ecco_;
		merge=merge_;
		minProb=(float)minProb_;
		rcomp=rcomp_;
//		assert(false) : tableMemory+", "+bytesPerKmer+", "+prealloc+", "+preallocFraction;
		estimatedKmerCapacity=(long)((tableMemory*1.0/bytesPerKmer)*((prealloc ? preallocFraction : 0.81))*(HashArray.maxLoadFactorFinal*0.97));

//		System.err.println("tableMemory="+tableMemory+", bytesPerKmer="+bytesPerKmer+", estimatedKmerCapacity="+estimatedKmerCapacity);
		KmerCountAbstract.minProb=(minProbPrefilter ? minProb : 0);
		k=k_;
		k2=k-1;
		int bitsPerBase=(amino ? 5 : 2);
		int mask=(amino ? 31 : 3);
		coreMask=(AbstractKmerTableSet.MASK_CORE ? ~(((-1L)<<(bitsPerBase*(k_-1)))|mask) : -1L);

		middleMask=AbstractKmerTableSet.MASK_MIDDLE ? AbstractKmerTable.makeMiddleMask(k, amino) : -1L;
		cmMask=coreMask&middleMask;

		if(k<1 || k>31){throw new RuntimeException("\nk needs an integer value from 1 to 31, such as k=27.  Default is 31.\n");}

		requestedInitialSize=initialSize;
		if(initialSize<1){
			final long memOverWays=tableMemory/(bytesPerKmer*ways);
			final double mem2=(prealloc ? preallocFraction : 1)*tableMemory;
			initialSize=(prealloc || memOverWays<initialSizeDefault ? (int)Tools.min(2142000000, (long)(mem2/(bytesPerKmer*ways))) : initialSizeDefault);
			if(initialSize!=initialSizeDefault){
				System.err.println("Initial size set to "+initialSize);
			}
		}

		/* Adjust I/O settings and filenames */

		assert(FastaReadInputStream.settingsOK());

		if(in1.isEmpty()){
			//throw new RuntimeException("Error - at least one input file is required.");
		}

		for(int i=0; i<in1.size(); i++){
			String s=in1.get(i);
			if(s!=null && s.contains("#") && !new File(s).exists()){
				int pound=s.lastIndexOf('#');
				String a=s.substring(0, pound);
				String b=s.substring(pound+1);
				in1.set(i, a+1+b);
				in2.add(a+2+b);
			}
		}

		if(!extra.isEmpty()){
			ArrayList<String> temp=(ArrayList<String>) extra.clone();
			extra.clear();
			for(String s : temp){
				if(s!=null && s.contains("#") && !new File(s).exists()){
					int pound=s.lastIndexOf('#');
					String a=s.substring(0, pound);
					String b=s.substring(pound+1);
					extra.add(a+1+b);
					extra.add(a+2+b);
				}else{
					extra.add(s);
				}
			}
		}

		{
			boolean allowDuplicates=true;
			if(!Tools.testInputFiles(allowDuplicates, true, in1, in2, extra)){
				throw new RuntimeException("\nCan't read some input files.\n");
			}
		}

		if(DISPLAY_PROGRESS){
//			assert(false);
			outstream.println("Initial:");
			outstream.println("Ways="+ways+", initialSize="+initialSize+", prefilter="+(prefilter ? "t" : "f")+", prealloc="+(prealloc ? (""+preallocFraction) : "f"));
			Shared.printMemory();
			outstream.println();
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Drops the shard reference only; inherited allocation flags/counters are not reset. */
	@Override
	public void clear(){
		tables=null;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Builds a growth schedule and allocates shards; requires no current table array.
	 * Uses the parsed initial-size request as a schedule lower bound; nonpositive requests
	 * keep the scheduler default, and preallocation selects its maximum instead.
	 * Uses the core mask captured for this table set's alphabet; does not load reads. */
	@Override
	public void allocateTables(){
		assert(tables==null);
		tables=null;
		//KMR015: backend cell geometry uses the table set's captured DNA or amino core mask.
		final long coreMask=this.coreMask;

		//KMR024: forward the raw request, not the separate fallback computed for reporting.
		ScheduleMaker scheduleMaker=new ScheduleMaker(ways, bytesPerKmer, prealloc,
				(prealloc ? preallocFraction : 1.0), requestedInitialSize, (prefilter ? prepasses : 0), prefilterFraction, filterMemoryOverride);
		int[] schedule=scheduleMaker.makeSchedule();
		//assert(false) : preallocFraction+", "+ways+", "+Arrays.toString(schedule);
//		System.err.println("DEBUG "+preallocFraction+", "+ways+", "+Arrays.toString(schedule));
		tables=AbstractKmerTable.preallocate(ways, tableType, schedule, coreMask);

//		tables=AbstractKmerTable.preallocate(ways, tableType, initialSize, coreMask, (!prealloc || preallocFraction<1));
	}

	/**
	 * Loads one file or pair through Shared.threads() workers and joins every worker.
	 * Closes input and propagates recorded startup, worker or join failures before merging counters (KMR011).
	 * Restores caller interruption after cleanup. Partial table writes are not rolled back or retried.
	 * Other workers may finish their input; blocked backend reads, writes or close calls must return
	 * before cleanup can finish. This method does not force cancellation of those operations.
	 * @return Newly inserted keys reported by worker buffers; not the number of input windows
	 */
	@Override
	public long loadKmers(String fname1, String fname2){

		final int THREADS=Shared.threads();
		final ArrayList<LoadThread> alpt=new ArrayList<LoadThread>(THREADS);
		Streamer cris=null;
		Throwable failure=null;
		boolean interrupted=false;
		try{
			if(Thread.interrupted()){
				interrupted=true;
				throw new InterruptedException("Interrupted before loading k-mers");
			}
			/* Create read input stream */
			FileFormat ff1=FileFormat.testInput(fname1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(fname2, FileFormat.FASTQ, null, true, true);
			cris=StreamerFactory.getReadInputStream(maxReads, false, ff1, ff2, -1);
			cris.start();

//		/* Optionally skip the first reads, since initial reads may have lower quality */
//		if(skipreads>0){
//			long skipped=0;
//
//			ListNum<Read> ln=cris.nextList();
//			ArrayList<Read> reads=(ln!=null ? ln.list : null);
//
//			while(skipped<skipreads && reads!=null && reads.size()>0){
//				skipped+=reads.size();
//
//				cris.returnList(ln);
//				ln=cris.nextList();
//				reads=(ln!=null ? ln.list : null);
//			}
//			cris.returnList(ln);
//			if(reads==null || reads.isEmpty()){
//				ReadWrite.closeStreams(cris);
//				System.err.println("Skipped all of the reads.");
//				System.exit(0);
//			}
//		}

			/* Construct before starting; cleanup also handles a partially started list. */
			for(int i=0; i<THREADS; i++){alpt.add(new LoadThread(cris));}
			for(LoadThread pt : alpt){pt.start();}
		}catch(Throwable t){
			failure=t;
		}finally{
			/* NEW workers never started; terminated workers have published their failure. */
			for(LoadThread pt : alpt){
				while(pt.isAlive()){
					try{
						pt.join();
					}catch(InterruptedException e){
						interrupted=true;
						if(failure==null){failure=e;}
					}
				}
				if(failure==null){failure=pt.failure;}
			}
			if(Thread.interrupted()){
				interrupted=true;
				if(failure==null){failure=new InterruptedException("Interrupted while loading k-mers");}
			}
			try{
				if(cris!=null){errorState|=ReadWrite.closeStreams(cris);}
			}catch(Throwable t){
				if(failure==null){failure=t;}
			}finally{
				if(Thread.interrupted()){
					interrupted=true;
					if(failure==null){failure=new InterruptedException("Interrupted while loading k-mers");}
				}
				if(failure!=null){errorState=true;}
				if(interrupted){Thread.currentThread().interrupt();}
			}
		}
		if(failure!=null){
			if(failure instanceof Error){throw (Error)failure;}
			if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
			throw new RuntimeException("Failed to load k-mers", failure);
		}

		/* No recorded failure may reach the counter reduction. Stream error flags remain in errorState. */
		long added=0;
		for(LoadThread pt : alpt){
			added+=pt.added;

			readsIn+=pt.readsInT;
			basesIn+=pt.basesInT;
			lowqReads+=pt.lowqReadsT;
			lowqBases+=pt.lowqBasesT;
			readsTrimmed+=pt.readsTrimmedT;
			basesTrimmed+=pt.basesTrimmedT;
			kmersIn+=pt.kmersInT;
		}

		return added;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Consumes shared input batches, preprocesses mates and counts into a private HashBuffer.
	 * Flushes on normal completion and records any failure, including the final flush, for the coordinator.
	 */
	private class LoadThread extends Thread{

		/**
		 * Constructor
		 * @param cris_ Started shared input stream; coordinator owns closure after all workers terminate
		 */
		public LoadThread(Streamer cris_){
			cris=cris_;
			table=new HashBuffer(tables, buflen, k, false, true, KmerTableSet.this.amino);
			//TODO: handle MASK_MIDDLE in this thread.
		}

		/** Processes batches and final flush, recording failure for the coordinator after termination.
		 * Read counters include individual mates. */
		@Override
		public void run(){
			try{
				ListNum<Read> ln=cris.nextList();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);

				//While there are more reads lists...
				while(ln!=null && reads!=null){//ln!=null prevents a compiler potential null access warning

					//For each read (or pair) in the list...
					for(int i=0; i<reads.size(); i++){
						Read r1=reads.get(i);
						Read r2=r1.mate;

						if(!r1.validated()){r1.validate(true);}
						if(r2!=null && !r2.validated()){r2.validate(true);}

						if(verbose){System.err.println("Considering read "+r1.id+" "+new String(r1.bases));}

						readsInT++;
						basesInT+=r1.length();
						if(r2!=null){
							readsInT++;
							basesInT+=r2.length();
						}

						if(maxNs<Integer.MAX_VALUE){
							if(r1!=null && r1.countUndefined()>maxNs){r1.setDiscarded(true);}
							if(r2!=null && r2.countUndefined()>maxNs){r2.setDiscarded(true);}
						}

						//Determine whether to discard the reads based on average quality
						if(minAvgQuality>0){
							if(r1!=null && r1.quality!=null && r1.avgQuality(false, minAvgQualityBases)<minAvgQuality){r1.setDiscarded(true);}
							if(r2!=null && r2.quality!=null && r2.avgQuality(false, minAvgQualityBases)<minAvgQuality){r2.setDiscarded(true);}
						}

						if(r1!=null){
							if(qtrimLeft || qtrimRight){
								int x=TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq, trimE, 1, false);
								basesTrimmedT+=x;
								readsTrimmedT+=(x>0 ? 1 : 0);
							}
							if(r1.length()<k){r1.setDiscarded(true);}
						}
						if(r2!=null){
							if(qtrimLeft || qtrimRight){
								int x=TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq, trimE, 1, false);
								basesTrimmedT+=x;
								readsTrimmedT+=(x>0 ? 1 : 0);
							}
							if(r2.length()<k){r2.setDiscarded(true);}
						}

						if((ecco || merge) && r1!=null && r2!=null && !r1.discarded() && !r2.discarded()){
							if(merge){
								final int insert=BBMerge.findOverlapStrict(r1, r2, false);
								if(insert>0){
									r2.reverseComplementFast();
									r1=r1.joinRead(insert);
									r2=null;
								}
							}else if(ecco){
								BBMerge.findOverlapStrict(r1, r2, true);
							}
						}

						if(minLen>0){
							if(r1!=null && r1.length()<minLen){r1.setDiscarded(true);}
							if(r2!=null && r2.length()<minLen){r2.setDiscarded(true);}
						}

						if(r1!=null){
							if(r1.discarded()){
								lowqBasesT+=r1.length();
								lowqReadsT++;
							}else{
								long temp=addKmersToTable(r1);
								added+=temp;
								if(verbose){System.err.println("A: Added "+temp);}
							}
						}
						if(r2!=null){
							if(r2.discarded()){
								lowqBasesT+=r2.length();
								lowqReadsT++;
							}else{
								long temp=addKmersToTable(r2);
								added+=temp;
								if(verbose){System.err.println("B: Added "+temp);}
							}
						}
					}
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				long temp=table.flush();
				if(verbose){System.err.println("Flush: Added "+temp);}
				added+=temp;
			}catch(Throwable t){
				failure=t;
			}
		}

		/**
		 * Adds amino acid k-mers from a read to the hash tables.
		 * Handles protein sequences with 5-bit encoding per amino acid.
		 * @param r Read containing amino acid sequence
		 * @return New keys reported by the buffer; zero for null data, -1 for a sequence shorter than k
		 */
		private final int addKmersToTableAA(final Read r){
			if(r==null || r.bases==null){return 0;}
			final float minProb2=(minProbMain ? minProb : 0);
			final byte[] bases=r.bases;
			final byte[] quals=r.quality;
			final int shift=5*k;
			final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
			long kmer=0;
			int created=0;
			int len=0;

			if(bases==null || bases.length<k){return -1;}

			/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
			float prob=1;
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				final long x=AminoAcid.acidToNumber[b];

				//Update kmers
				kmer=((kmer<<5)|x)&mask;

				if(minProb2>0 && quals!=null){//Update probability
					prob=prob*PROB_CORRECT[quals[i]];
					if(len>=k){
						byte oldq=quals[i-k];
						prob=prob*PROB_CORRECT_INVERSE[oldq];
					}
				}

				//Handle Ns
				if(x<0){
					len=0;
					kmer=0;
					prob=1;
				}else{len++;}

				if(verbose){System.err.println("Scanning i="+i+", len="+len+", kmer="+kmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
				if(len>=k && prob>=minProb2){
					kmersInT++;
					final long key=kmer;
					if(!prefilter || prefilterArray.read(key)>filterMax2){
						int temp=table.incrementAndReturnNumCreated(key, 1);
						created+=temp;
						if(verbose){System.err.println("C: Added "+temp);}
					}
				}
			}

			return created;
		}

		/**
		 * Adds DNA k-mers from a read to the hash tables.
		 * Dispatches to appropriate method based on processing mode (amino/onePass/normal).
		 * @param r Read containing DNA sequence
		 * @return New keys reported by the buffer; zero for null data, -1 for a sequence shorter than k
		 */
		private final int addKmersToTable(final Read r){
			if(amino){return addKmersToTableAA(r);}
			if(onePass){return addKmersToTable_onePass(r);}
			if(r==null || r.bases==null){return 0;}
			final float minProb2=(minProbMain ? minProb : 0);
			final byte[] bases=r.bases;
			final byte[] quals=r.quality;
			final int shift=2*k;
			final int shift2=shift-2;
			final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
			long kmer=0;
			long rkmer=0;
			int created=0;
			int len=0;

			if(bases==null || bases.length<k){return -1;}

			/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
			float prob=1;
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];

				//Update kmers
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

				if(minProb2>0 && quals!=null){//Update probability
					prob=prob*PROB_CORRECT[quals[i]];
					if(len>=k){
						byte oldq=quals[i-k];
						prob=prob*PROB_CORRECT_INVERSE[oldq];
					}
				}

				//Handle Ns
				if(x<0){
					len=0;
					kmer=rkmer=0;
					prob=1;
				}else{len++;}

				if(verbose){System.err.println("Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
				if(len>=k && prob>=minProb2){
					kmersInT++;
					final long key=toValue(kmer, rkmer);
					if(!prefilter || prefilterArray.read(key)>filterMax2){
						int temp=table.incrementAndReturnNumCreated(key, 1);
						created+=temp;
						if(verbose){System.err.println("C: Added "+temp);}
					}
				}
			}

			return created;
		}

		/**
		 * Adds k-mers using single-pass prefiltering mode.
		 * Prefilter support and main-table counting obey independent probability flags.
		 * Counts main-quality-accepted input k-mers before rejecting insufficient depth.
		 * @param r Read to process
		 * @return Number of new main-table keys, or -1 for a read shorter than k
		 */
		private final int addKmersToTable_onePass(final Read r){
			assert(prefilter);
			if(r==null || r.bases==null){return 0;}
			final byte[] bases=r.bases;
			final byte[] quals=r.quality;
			//As in the normal loader, use one configured probability policy throughout this read.
			final boolean checkFilter=minProbPrefilter, checkMain=minProbMain;
			final boolean updateProb=minProb>0 && (checkFilter || checkMain) && quals!=null;
			final int shift=2*k;
			final int shift2=shift-2;
			final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
			long kmer=0;
			long rkmer=0;
			int created=0;
			int len=0;

			if(bases==null || bases.length<k){return -1;}

			/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
			float prob=1;
			for(int i=0; i<bases.length; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];

				//Update kmers
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

				if(updateProb){//Update probability
					prob=prob*PROB_CORRECT[quals[i]];
					if(len>=k){
						byte oldq=quals[i-k];
						prob=prob*PROB_CORRECT_INVERSE[oldq];
					}
				}

				//Handle Ns
				if(x<0){
					len=0;
					kmer=rkmer=0;
					prob=1;
				}else{len++;}

				if(verbose){System.err.println("Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
				if(len>=k){
					final boolean filterAccepted=(!checkFilter || prob>=minProb);
					final boolean mainAccepted=(!checkMain || prob>=minProb);
					if(!filterAccepted && !mainAccepted){continue;}
					final long key=toValue(kmer, rkmer);
					//Increment returns old support; a skipped increment must already exceed the threshold.
					final boolean supported=(filterAccepted ?
						prefilterArray.incrementAndReturnUnincremented(key, 1)>=filterMax2 :
						prefilterArray.read(key)>filterMax2);
					if(mainAccepted){
						kmersInT++;
						if(supported){
							int temp=table.incrementAndReturnNumCreated(key, 1);
							created+=temp;
							if(verbose){System.err.println("D: Added "+temp);}
						}
					}
				}
			}
			return created;
		}

		/*--------------------------------------------------------------*/

		/** Input read stream */
		private final Streamer cris;

		/** Hash buffer for distributing k-mers across tables */
		private final HashBuffer table;

		/** Written by this worker; read only after the coordinator observes its termination. */
		private Throwable failure;

		/** Newly created table keys reported by this worker, including its final buffer flush. */
		public long added=0;

		/** Number of reads processed by this thread */
		private long readsInT=0;
		/** Total bases processed by this thread */
		private long basesInT=0;
		/** Individual reads discarded for any configured quality/length/undefined-base reason. */
		private long lowqReadsT=0;
		/** Post-trimming bases in discarded individual reads. */
		private long lowqBasesT=0;
		/** Number of reads that were quality trimmed by this thread */
		private long readsTrimmedT=0;
		/** Total bases removed by quality trimming in this thread */
		private long basesTrimmedT=0;
		/** Valid main-quality-accepted windows before prefilter-depth rejection, including repeats. */
		private long kmersInT=0;

	}

	/*--------------------------------------------------------------*/
	/*----------------          Convenience         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Recomputes the next at-most-k forward DNA windows after an existing window at a.
	 * Uses kmers[a] as the seed, leaves that entry unchanged and updates both lists
	 * through min(counts.size-1,a+k). Undefined windows become key -1/count zero.
	 *
	 * @param bases Sequence bases
	 * @param kmers List to store generated k-mers
	 * @param counts List to store k-mer counts
	 * @param a Valid start index of the seed k-mer in already populated lists
	 */
	public void regenerateKmers(byte[] bases, LongList kmers, IntList counts, final int a){
		final int loc=a+k;
		final int lim=Tools.min(counts.size, a+k+1);
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=kmers.get(a);
		long rkmer=rcomp(kmer);
		int len=k;

//		assert(false) : a+", "+loc+", "+lim;

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=loc, j=a+1; j<lim; i++, j++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			final long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{len++;}

			if(len>=k){
				assert(kmers.get(j)!=kmer);
				kmers.set(j, kmer);
				int count=getCount(kmer, rkmer);
				counts.set(j, count);
			}else{
				kmers.set(j, -1);
				counts.set(j, 0);
			}
		}
	}

	/** Recomputes DNA counts overlapping the nonempty changed-base bitset.
	 * Keeps unaffected entries; invalid windows become zero. dummy is unused.
	 * @return Number of valid recomputed windows, not total count-list length */
	@Override
	public int regenerateCounts(byte[] bases, IntList counts, final Kmer dummy, BitSet changed){
		assert(!changed.isEmpty());
		final int firstBase=changed.nextSetBit(0), lastBase=changed.length()-1;
		final int ca=firstBase-k;
//		final int b=changed.nextSetBit(0);ca+kbig; //first base changed
		final int firstCount=Tools.max(firstBase-k+1, 0), lastCount=Tools.min(counts.size-1, lastBase); //count limit
//		System.err.println("ca="+ca+", b="+b+", lim="+lim);
//		System.err.println("Regen from count "+(ca+1)+"-"+lim);

		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;
		int len=0;
		int valid=0;

//		System.err.println("ca="+ca+", b="+b+", lim="+lim+", "+counts);

		for(int i=Tools.max(0, firstBase-k+1), lim=Tools.min(lastBase+k-1, bases.length-1); i<=lim; i++){
			final byte base=bases[i];
			final long x=AminoAcid.baseToNumber[base];
			final long x2=AminoAcid.baseToComplementNumber[base];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{
				len++;
			}

			final int c=i-k+1;
			if(i>=firstBase){
				if(len>=k){
					valid++;
					int count=getCount(kmer, rkmer);
					counts.set(c, count);
				}else if(c>=0){
					counts.set(c, 0);
				}
			}
		}

		return valid;
	}

	/** Clears/refills DNA counts through the last selected window (all windows if positions is null).
	 * Unselected or invalid windows become zero; absent keys are clamped to zero.
	 * @return Number of valid selected windows; dummy is unused */
	@Override
	public int fillSpecificCounts(byte[] bases, IntList counts, BitSet positions, final Kmer dummy){
		final int b=k-1;
		final int lim=(positions==null ? bases.length : Tools.min(bases.length, positions.length()+k-1));
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;
		int len=0;
		int valid=0;

		counts.clear();

		for(int i=0, j=-b; i<lim; i++, j++){
			final byte base=bases[i];
			final long x=AminoAcid.baseToNumber[base];
			final long x2=AminoAcid.baseToComplementNumber[base];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{
				len++;
			}

			if(j>=0){
				if(len>=k && (positions==null || positions.get(j))){
					valid++;
					int count=getCount(kmer, rkmer);
					assert(i-k+1==counts.size()) : "j="+j+", counts="+counts.size()+", b="+(b)+", (i-k+1)="+(i-k+1);
					assert(j==counts.size());
					counts.add(Tools.max(count, 0));
//					counts.set(i-k+1, count);
				}else{
					counts.add(0);
//					counts.set(i-k+1, 0);
				}
			}
		}

//		assert((counts.size==0 && bases.length<k) || counts.size==bases.length-k+1) : bases.length+", "+k+", "+counts.size;
		assert((counts.size==0 && bases.length<k) || counts.size==lim-k+1) : bases.length+", "+k+", "+counts.size;

		return valid;
	}

	/**
	 * Regenerates counts for k-mers in a specific region around position ca.
	 * Updates existing counts list without changing its size.
	 *
	 * @param bases Sequence bases
	 * @param counts List to update with new counts
	 * @param ca First k-mer start index to recompute; nonnegative with existing count storage
	 * @return Number of valid k-mers processed
	 */
	public int regenerateCounts(byte[] bases, IntList counts, final int ca){
		final int b=ca+k-1;
		final int lim=Tools.min(bases.length, b+k+1);
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;
		int len=0;
		int valid=0;

		final int clen=counts.size;

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts.
		 * i is an index in the base array, j is an index in the count array. */
		for(int i=b-k+1; i<lim; i++){
			final byte base=bases[i];
			final long x=AminoAcid.baseToNumber[base];
			final long x2=AminoAcid.baseToComplementNumber[base];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{
				len++;
			}

			if(i>=b){
				if(len>=k){
					valid++;
					int count=getCount(kmer, rkmer);
					counts.set(i-k+1, count);
				}else{
					counts.set(i-k+1, 0);
				}
			}
		}

		assert((counts.size==0 && bases.length<k) || counts.size==bases.length-k+1) : bases.length+", "+k+", "+counts.size;
		assert(clen==counts.size);

		return valid;
	}

	/** Fills forward DNA windows and -1 invalid markers, returning the valid-window count.
	 * Replaces the destination's logical contents; input shorter than k leaves it empty.
	 * Unused backing-array entries are retained. */
	public int fillKmers(byte[] bases, LongList kmers){
		final int blen=bases.length;
		//KMR010: even a short input replaces the previous logical window list.
		if(blen<k){kmers.clear(); return 0;}
		final int min=k-1;
		final int shift=2*k;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0;
		int len=0;
		int valid=0;

		kmers.clear();

		/* Loop through the bases, maintaining a forward kmer via bitshifts */
		for(int i=0; i<blen; i++){
			final byte b=bases[i];
			assert(b>=0) : Arrays.toString(bases);
			final long x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x<0){
				len=0;
				kmer=0;
			}else{len++;}
			if(i>=min){
				if(len>=k){
					kmers.add(kmer);
					valid++;
				}else{
					kmers.add(-1);
				}
			}
		}
		return valid;
	}

	/**
	 * Retrieves counts for each k-mer in the input list, considering reverse complements.
	 * Clears and refills counts; invalid negative keys yield zero. Absent nonnegative
	 * keys retain the selected table's absence sentinel instead of being clamped.
	 * @param kmers List of k-mers to query
	 * @param counts List to populate with k-mer counts
	 */
	public void fillCounts(LongList kmers, IntList counts){
		counts.clear();
		for(int i=0; i<kmers.size; i++){
			long kmer=kmers.get(i);
			if(kmer>=0){
				long rkmer=rcomp(kmer);
				int count=getCount(kmer, rkmer);
				counts.add(count);
			}else{
				counts.add(0);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Delegates threshold removal to each shard and sums removals; caller clears ownership first. */
	@Override
	public long regenerate(final int limit){
		long sum=0;
		for(AbstractKmerTable akt : tables){
			sum+=akt.regenerate(limit);
		}
		return sum;
	}

	/**
	 * Returns the specific table that should contain the given k-mer key.
	 * Uses hash distribution to determine the appropriate table index.
	 * @param key K-mer key value
	 * @return The assigned HashArray1D shard; the key need not be present
	 * @throws ClassCastException If a different table family was configured
	 */
	public HashArray1D getTableForKey(long key){
		return (HashArray1D) tables[kmerToWay(key)];
	}

	/** Returns the indexed live shard cast to HashArray1D; requires that configured family. */
	@Override
	public HashArray1D getTable(int tnum){
		return (HashArray1D) tables[tnum];
	}

	/** Builds the count histogram through HistogramMaker, including its worker-error handling. */
	@Override
	public long[] fillHistogram(int histMax){
		return HistogramMaker.fillHistogram(tables, histMax);
	}

	/** Initializes shard ownership through joined workers after loading; caller excludes mutation. */
	@Override
	public void initializeOwnership(){
		OwnershipThread.initialize(tables);
	}

	/** Clears shard ownership through joined workers; caller excludes concurrent claims. */
	@Override
	public void clearOwnership(){
		OwnershipThread.clear(tables);
	}

	/**
	 * Extracts the forward k-mer from the final k bases of a sequence buffer.
	 * @param bb Sequence buffer; only its logical length is examined
	 * @return Forward k-mer, or -1 if the sequence is too short or its suffix is undefined
	 */
	public long rightmostKmer(final ByteBuilder bb){
		return rightmostKmer(bb.array, bb.length());
	}

	/**
	 * Encodes the final k bases of the logical DNA sequence, without scanning earlier windows.
	 * Processes the last k bases to generate the terminal k-mer.
	 *
	 * @param bases Sequence bases
	 * @param blen Length of sequence to process
	 * @return Rightmost k-mer value, or -1 if no valid k-mer found
	 */
	public long rightmostKmer(final byte[] bases, final int blen){
		if(blen<k){return -1;}
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0;
		long rkmer=0;
		int len=0;

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts, to get the rightmost kmer */
		{
			for(int i=blen-k; i<blen; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{len++;}
				if(verbose){outstream.println("Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			}
		}

		if(len<k){return -1;}
		else{assert(len==k);}
		return kmer;
	}

	/**
	 * Encodes the first k bases of the logical DNA buffer, without scanning later windows.
	 * @param bb ByteBuilder containing sequence data
	 * @return Leftmost k-mer value, or -1 if no valid k-mer found
	 */
	public long leftmostKmer(final ByteBuilder bb){
		return leftmostKmer(bb.array, bb.length());
	}

	/**
	 * Encodes the first k bases of the logical DNA sequence, without scanning later windows.
	 * Processes the first k bases to generate the initial k-mer.
	 *
	 * @param bases Sequence bases
	 * @param blen Length of sequence to process
	 * @return Leftmost k-mer value, or -1 if no valid k-mer found
	 */
	public long leftmostKmer(final byte[] bases, final int blen){
		if(blen<k){return -1;}
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0;
		long rkmer=0;
		int len=0;

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts, to get the leftmost kmer */
		{
			for(int i=0; i<k; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{len++;}
				if(verbose){outstream.println("Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			}
		}

		if(len<k){return -1;}
		else{assert(len==k);}
		return kmer;
	}

	/**
	 * Attempts to claim ownership of all k-mers in sequence twice with different IDs.
	 * @param bb ByteBuilder containing sequence data
	 * @param id Base thread/owner ID
	 * @return true if both claim attempts succeeded
	 */
	public boolean doubleClaim(final ByteBuilder bb, final int id){
		return doubleClaim(bb.array, bb.length(), id);
	}

	/** Makes two best-effort claim passes, first with id then id+CLAIM_OFFSET.
	 * Failure leaves earlier claims intact; this is not an atomic whole-sequence claim. */
	public boolean doubleClaim(final byte[] bases, final int blength, final int id){
		boolean success=claim(bases, blength, id, true);
		if(verbose){outstream.println("success1="+success+", id="+id+", blength="+blength);}
		if(!success){return false;}
		success=claim(bases, blength, id+CLAIM_OFFSET, true);
		if(verbose){outstream.println("success2="+success+", id="+id+", blength="+blength);}
		return success;
	}

	/**
	 * Attempts to claim ownership of all k-mers in the ByteBuilder sequence.
	 *
	 * @param bb ByteBuilder containing sequence data
	 * @param id Thread/owner ID
	 * @param exitEarly Whether to stop on first failure
	 * @return true if all k-mers were successfully claimed
	 */
	public boolean claim(final ByteBuilder bb, final int id, final boolean exitEarly){
		return claim(bb.array, bb.length(), id, exitEarly);
	}

	/**
	 * Calculates average k-mer coverage for the requested sequence prefix.
	 * Absent k-mers contribute zero and remain in the valid-window denominator.
	 * Windows containing undefined bases are skipped; the input is unchanged.
	 *
	 * @param bases Sequence bases
	 * @param blength Length of sequence to analyze
	 * @return Average k-mer coverage, or 0 if no valid k-mers
	 */
	public float calcCoverage(final byte[] bases, final int blength){
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;
		int len=0;
		long sum=0, max=0;
		int kmers=0;

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=0; i<blength; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			final long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{len++;}
			if(len>=k){
				final int rawCount=getCount(kmer, rkmer);
				//Signed shift makes the mask zero for negative counts, all ones otherwise.
				final int count=rawCount&~(rawCount>>31);
				sum+=count;
				max=Tools.max(count, max);
				kmers++;
			}
		}
		return sum==0 ? 0 : sum/(float)kmers;
	}

	/**
	 * Calculates k-mer coverage statistics for a contig sequence.
	 * Updates coverage/maxCov/minCov unless the contig is shorter than k (early return).
	 * Absent k-mers contribute zero and remain in the valid-window denominator.
	 * Windows containing undefined bases are skipped. For contigs of at least k bases,
	 * no valid windows means all three statistics are zero. The base array is unchanged.
	 * @param contig Contig object to analyze and update
	 * @return Average coverage, or zero if no valid k-mers or the contig is shorter than k
	 */
	public float calcCoverage(final Contig contig){
		final byte[] bases=contig.bases;
		if(bases.length<k){return 0;}
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;
		int len=0;
		long sum=0;
		int max=0, min=Integer.MAX_VALUE;
		int kmers=0;

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			final long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{len++;}
			if(len>=k){
				//KMR022 [kmer/KmerTableSet#003]: absent keys have zero coverage, not the lookup sentinel -1.
				final int rawCount=getCount(kmer, rkmer);
				//Signed shift makes the mask zero for negative counts, all ones otherwise.
				final int count=rawCount&~(rawCount>>31);
				sum+=count;
				max=Tools.max(count, max);
				min=Tools.min(count, min);
				kmers++;
			}
		}
		contig.coverage=sum==0 ? 0 : sum/(float)kmers;
		contig.maxCov=max;
		contig.minCov=sum==0 ? 0 : min;
		return contig.coverage;
	}

	/**
	 * Attempts to claim ownership of all k-mers in the sequence for the given thread ID.
	 * With exitEarly=true, stops at the first conflicting owner without rollback.
	 * With exitEarly=false, visits every valid window and returns true even if a claim failed.
	 *
	 * @param bases Sequence bases
	 * @param blength Length of sequence
	 * @param id Thread/owner ID
	 * @param exitEarly Whether to stop processing on first failed claim
	 * @return true if all claims succeeded (or exitEarly is false)
	 */
	public boolean claim(final byte[] bases, final int blength, final int id, boolean exitEarly){
		if(verbose){outstream.println("Thread "+id+" claim start.");}
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;
		int len=0;
		boolean success=true;
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=0; i<blength && success; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			final long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{len++;}
			if(verbose){System.err.println("Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			if(len>=k){
				success=claim(kmer, rkmer, id/*, rid, i*/);
				success=(success || !exitEarly);
			}
		}
		return success;
	}

	/**
	 * Attempts to claim ownership of a specific k-mer/reverse-complement pair.
	 *
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @param id Thread/owner ID
	 * @return true if the key is absent or its observed owner equals id
	 */
	public boolean claim(final long kmer, final long rkmer, final int id/*, final long rid, final int pos*/){
		//TODO: rid and pos are just for debugging.
		final long key=toValue(kmer, rkmer);
		final int way=kmerToWay(key);
		final AbstractKmerTable table=tables[way];
		final int count=table.getValue(key);
		assert(count==-1 || count>0) : count;
//		if(verbose  /*|| true*/){outstream.println("Count="+count+".");}
		if(count<0){return true;}
		assert(count>0) : count;
		final int owner=table.setOwner(key, id);
		if(verbose){outstream.println("owner="+owner+".");}
//		assert(owner==id) : id+", "+owner+", "+rid+", "+pos;
		return owner==id;
	}

	/**
	 * Releases ownership of all k-mers in the ByteBuilder sequence.
	 * @param bb ByteBuilder containing sequence data
	 * @param id Thread/owner ID that currently owns the k-mers
	 */
	public void release(ByteBuilder bb, final int id){
		release(bb.array, bb.length(), id);
	}

	/**
	 * Releases ownership of all k-mers in the sequence for the given thread ID.
	 * @param bases Sequence bases
	 * @param blength Length of sequence
	 * @param id Thread/owner ID that currently owns the k-mers
	 */
	public void release(final byte[] bases, final int blength, final int id){
		if(verbose  /*|| true*/){outstream.println("*Thread "+id+" release start.");}
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;
		int len=0;
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=0; i<blength; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			final long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{len++;}
			if(verbose){System.err.println("Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			if(len>=k){
				release(kmer, rkmer, id);
			}
		}
	}

	/**
	 * Releases ownership of a specific k-mer/reverse-complement pair.
	 *
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @param id Thread/owner ID
	 * @return true if the count is below one or the matching owner was cleared
	 */
	public boolean release(final long kmer, final long rkmer, final int id){
		return release(toValue(kmer, rkmer), id);
	}

	/**
	 * Releases ownership of a k-mer key for the specified thread ID.
	 * @param key Canonical k-mer key
	 * @param id Thread/owner ID
	 * @return true if the count is below one or the matching owner was cleared
	 */
	public boolean release(final long key, final int id){
		final int way=kmerToWay(key);
		final AbstractKmerTable table=tables[way];
		final int count=table.getValue(key);
//		if(verbose  /*|| true*/){outstream.println("Count="+count+".");}
		if(count<1){return true;}
		return table.clearOwner(key, id);
	}

	/**
	 * Finds the maximum owner ID among all k-mers in the ByteBuilder sequence.
	 * @param bb ByteBuilder containing sequence data
	 * @param id Reference ID (stops searching if owner exceeds this)
	 * @return Highest owner observed before any early exit, or -1 if none
	 */
	public int findOwner(ByteBuilder bb, final int id){
		return findOwner(bb.array, bb.length(), id);
	}

	/**
	 * Finds the maximum owner ID among all k-mers in the sequence.
	 * Stops early if any owner ID exceeds the reference ID parameter.
	 *
	 * @param bases Sequence bases
	 * @param blength Length of sequence
	 * @param id Reference ID for early termination
	 * @return Maximum owner observed before termination, or -1 if none; not a whole-sequence maximum after early exit
	 */
	public int findOwner(final byte[] bases, final int blength, final int id){
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;
		int len=0;
		int maxOwner=-1;
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=0; i<blength; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			final long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{len++;}
			if(verbose){System.err.println("Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			if(len>=k){
				int owner=findOwner(kmer, rkmer);
				maxOwner=Tools.max(owner, maxOwner);
				if(maxOwner>id){break;}
			}
		}
		return maxOwner;
	}

	/**
	 * Finds the owner of a k-mer, considering its reverse complement.
	 * @param kmer K-mer to query
	 * @return Owner ID, or -1 if k-mer not found
	 */
	public int findOwner(final long kmer){
		return findOwner(kmer, rcomp(kmer));
	}

	/**
	 * Finds the owner of a specific k-mer/reverse-complement pair.
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @return Owner ID, or -1 if k-mer not found
	 */
	public int findOwner(final long kmer, final long rkmer){
		final long key=toValue(kmer, rkmer);
		final int way=kmerToWay(key);
		final AbstractKmerTable table=tables[way];
		final int count=table.getValue(key);
		if(count<0){return -1;}
		final int owner=table.getOwner(key);
		return owner;
	}

	/**
	 * Retrieves the count for a k-mer/reverse-complement pair.
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @return Stored count or the table's absence sentinel; negative results are not clamped
	 */
	public int getCount(long kmer, long rkmer){
		long key=toValue(kmer, rkmer);
		int way=kmerToWay(key);
		return tables[way].getValue(key);
	}

	/**
	 * Retrieves the count for a canonical k-mer key.
	 * @param key Canonical k-mer key
	 * @return Stored count or the table's absence sentinel; negative results are not clamped
	 */
	public int getCount(long key){
		int way=kmerToWay(key);
		return tables[way].getValue(key);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Fill Counts         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Fills counts array for k-mers extending to the right (3' direction).
	 * Uses fast implementation when core masking is enabled, otherwise uses safe version.
	 *
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @param counts Array to store counts for each possible extension (A,C,G,T)
	 * @param mask Bit mask for k-mer length
	 * @param shift2 Bit shift value for reverse complement operations
	 * @return Index of the first maximum in A,C,G,T order (zero when all counts are zero)
	 */
	public int fillRightCounts(long kmer, long rkmer, int[] counts, long mask, int shift2){
		if(FAST_FILL && MASK_CORE && k>2/*((k&1)==1)*/){
			return fillRightCounts_fast(kmer, rkmer, counts, mask, shift2);
		}else{
			return fillRightCounts_safe(kmer, rkmer, counts, mask, shift2);
		}
	}

	/**
	 * Fast implementation for filling right extension counts using core masking optimization.
	 * Reuses a primary probe start for four extensions; requires compatible HashArray
	 * shards and canonical DNA/core-mask configuration. Equal cores use the safe path.
	 *
	 * @param kmer0 Original forward k-mer
	 * @param rkmer0 Original reverse complement k-mer
	 * @param counts Array to store counts for each possible extension
	 * @param mask Bit mask for k-mer length
	 * @param shift2 Bit shift value for reverse complement operations
	 * @return Index of the first maximum in A,C,G,T order (zero when all counts are zero)
	 */
	public int fillRightCounts_fast(final long kmer0, final long rkmer0, int[] counts,
			long mask, int shift2){
//		assert((k&1)==1) : k;
		assert(MASK_CORE);
		final long kmer=(kmer0<<2)&mask;
		final long rkmer=(rkmer0>>>2);
		final long a=kmer&coreMask, b=rkmer&coreMask;

		int max=-1, maxPos=0;
		if(a==b){
			return fillRightCounts_safe(kmer0, rkmer0, counts, mask, shift2);
		}else if(a>b){
//			return fillRightCounts_safe(kmer0, rkmer0, counts, mask, shift2);
			final long core=a;
			final int way=kmerToWay(core);
//			final AbstractKmerTable table=tables[way];
			final HashArray table=(HashArray) tables[way];
			final int cell=table.kmerToCell(kmer);
			for(int i=0; i<=3; i++){
				final long key=kmer|i;
				final int count=Tools.max(0, table.getValue(key, cell));
//				final int count=Tools.max(0, table.getValue(key));
				counts[i]=count;
				if(count>max){
					max=count;
					maxPos=i;
				}
			}
		}else{
//			return fillRightCounts_safe(kmer0, rkmer0, counts, mask, shift2);
			//Use a higher increment for the left bits
			final long core=b;
			final int way=kmerToWay(core);
//			final AbstractKmerTable table=tables[way];
			final HashArray table=(HashArray) tables[way];
			final int cell=table.kmerToCell(rkmer);
			final long incr=1L<<(2*(k-1));
			long j=3*incr;
			for(int i=0; i<=3; i++, j-=incr){
				final long key=rkmer|j;
				final int count=Tools.max(0, table.getValue(key, cell));
//				final int count=Tools.max(0, table.getValue(key));
				counts[i]=count;
				if(count>max){
					max=count;
					maxPos=i;
				}
			}
		}
		return maxPos;
	}

	//TODO: Change this to take advantage of coreMask
	//Requires special handling of core palindromes.
	//Thus it would be easiest to just handle odd kmers, and K is normally 31 anyway.
	/**
	 * Safe implementation for filling right extension counts.
	 * Computes each DNA extension key separately and clamps absent counts to zero.
	 *
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @param counts Array to store counts for each possible extension
	 * @param mask Bit mask for k-mer length
	 * @param shift2 Bit shift value for reverse complement operations
	 * @return Index of the first maximum in A,C,G,T order (zero when all counts are zero)
	 */
	public int fillRightCounts_safe(long kmer, long rkmer, int[] counts, long mask, int shift2){
		assert(kmer==rcomp(rkmer));
		if(verbose){outstream.println("fillRightCounts:   "+toText(kmer)+",   "+toText(rkmer));}
		kmer=(kmer<<2)&mask;
		rkmer=(rkmer>>>2);
		int max=-1, maxPos=0;

		for(int i=0; i<=3; i++){
			long kmer2=kmer|((long)i);
			long rkmer2=rkmer|(((long)AminoAcid.numberToComplement[i])<<shift2);
			if(verbose){outstream.println("kmer:               "+toText(kmer2)+", "+toText(rkmer2));}
			assert(kmer2==(kmer2&mask));
			assert(rkmer2==(rkmer2&mask));
			assert(kmer2==rcomp(rkmer2));
			long key=toValue(kmer2, rkmer2);
			int way=kmerToWay(key);
			int count=tables[way].getValue(key);
			assert(count==NOT_PRESENT || count>=0);
			count=Tools.max(count, 0);
			counts[i]=count;
			if(count>max){
				max=count;
				maxPos=i;
			}
		}
		return maxPos;
	}

	/** Fills four DNA right-extension counts using only reverse-complement keys.
	 * For KmerCompressor; bypasses canonical selection and returns the first maximum index. */
	public int fillRightCountsRcompOnly(long kmer, long rkmer, int[] counts, long mask, int shift2){
		assert(kmer==rcomp(rkmer));
		if(verbose){outstream.println("fillRightCounts:   "+toText(kmer)+",   "+toText(rkmer));}
		kmer=(kmer<<2)&mask;
		rkmer=(rkmer>>>2);
		int max=-1, maxPos=0;

		for(int i=0; i<=3; i++){
			long kmer2=kmer|((long)i);
			long rkmer2=rkmer|(((long)AminoAcid.numberToComplement[i])<<shift2);
			if(verbose){outstream.println("kmer:               "+toText(kmer2)+", "+toText(rkmer2));}
			assert(kmer2==(kmer2&mask));
			assert(rkmer2==(rkmer2&mask));
			assert(kmer2==rcomp(rkmer2));
			long key=rkmer2;
			int way=kmerToWay(key);
			int count=tables[way].getValue(key);
			assert(count==NOT_PRESENT || count>=0);
			count=Tools.max(count, 0);
			counts[i]=count;
			if(count>max){
				max=count;
				maxPos=i;
			}
		}
		return maxPos;
	}

	/**
	 * Fills counts array for k-mers extending to the left (5' direction).
	 * Uses fast implementation when core masking is enabled, otherwise uses safe version.
	 *
	 * @param kmer Forward k-mer
	 * @param rkmer Reverse complement k-mer
	 * @param counts Array to store counts for each possible extension (A,C,G,T)
	 * @param mask Bit mask for k-mer length
	 * @param shift2 Bit shift value for reverse complement operations
	 * @return Index of the first maximum in A,C,G,T order (zero when all counts are zero)
	 */
	public int fillLeftCounts(long kmer, long rkmer, int[] counts, long mask, int shift2){
		if(FAST_FILL && MASK_CORE && k>2/*((k&1)==1)*/){
			return fillLeftCounts_fast(kmer, rkmer, counts, mask, shift2);
		}else{
			return fillLeftCounts_safe(kmer, rkmer, counts, mask, shift2);
		}
	}

	/**
	 * Fast implementation for filling left extension counts using core masking optimization.
	 * Reuses a primary probe start for four extensions; requires compatible HashArray
	 * shards and canonical DNA/core-mask configuration. Equal cores use the safe path.
	 *
	 * @param kmer0 Original forward k-mer
	 * @param rkmer0 Original reverse complement k-mer
	 * @param counts Array to store counts for each possible extension
	 * @param mask Bit mask for k-mer length
	 * @param shift2 Bit shift value for reverse complement operations
	 * @return Index of the first maximum in A,C,G,T order (zero when all counts are zero)
	 */
	public int fillLeftCounts_fast(final long kmer0, final long rkmer0, int[] counts,
			long mask, int shift2){
//		assert((k&1)==1) : k;
		assert(MASK_CORE);
		final long rkmer=(rkmer0<<2)&mask;
		final long kmer=(kmer0>>>2);
		final long a=kmer&coreMask, b=rkmer&coreMask;

		int max=-1, maxPos=0;
		if(a==b){
			return fillLeftCounts_safe(kmer0, rkmer0, counts, mask, shift2);
		}else if(a>b){
//			return fillLeftCounts_safe(kmer0, rkmer0, counts, mask, shift2);

			final long core=a;
			final int way=kmerToWay(core);
//			final AbstractKmerTable table=tables[way];
			final HashArray table=(HashArray) tables[way];
			final int cell=table.kmerToCell(kmer);
			final long incr=1L<<(2*(k-1));
			long j=0;//Must be declared as a long, outside of loop
			for(int i=0; i<=3; i++, j+=incr){
//				final long rkmer2=rkmer|((long)AminoAcid.numberToComplement[i]);
//				final long kmer2=kmer|(((long)i)<<shift2);
//				final long key2=toValue(rkmer2, kmer2);
//				int way2=kmerToWay(key2);
//				int count2=tables[way2].getValue(key2);
//				count2=Tools.max(count2, 0);

				final long key=kmer|j;
				final int count=Tools.max(0, table.getValue(key, cell));
//				int count=Tools.max(0, table.getValue(key));

//				assert(way==way2) : i+", "+way+", "+way2;
//				assert(key==key2) : i+", "+Long.toBinaryString(kmer)+", "+
//					Long.toBinaryString(key)+", "+Long.toBinaryString(key2)+
//					", "+Long.toBinaryString(j)+", "+Long.toBinaryString(incr);
//				assert(count==count2) : i+", "+count+", "+count2;

				counts[i]=count;
				if(count>max){
					max=count;
					maxPos=i;
				}
			}
			return maxPos;
		}else{
//			return fillLeftCounts_safe(kmer0, rkmer0, counts, mask, shift2);
			final long core=b;
			final int way=kmerToWay(core);
//			final AbstractKmerTable table=tables[way];
			final HashArray table=(HashArray) tables[way];
			final int cell=table.kmerToCell(rkmer);
			for(int i=0, j=3; i<=3; i++, j--){
				final long key=rkmer|j;
				final int count=Tools.max(0, table.getValue(key, cell));
//				final int count=Tools.max(0, table.getValue(key));

				counts[i]=count;
				if(count>max){
					max=count;
					maxPos=i;
				}
			}
			return maxPos;
		}
	}

	/**
	 * Computes left-extension counts separately and clamps absent counts to zero.
	 * Selects the forward extension when rcomp is false, or its canonical key otherwise.
	 *
	 * @param kmer0 Forward k-mer
	 * @param rkmer0 Reverse complement k-mer
	 * @param counts Array to store counts for each possible extension
	 * @param mask Bit mask for k-mer length
	 * @param shift2 Bit shift value for reverse complement operations
	 * @return Index of the first maximum in A,C,G,T order (zero when all counts are zero)
	 */
	public int fillLeftCounts_safe(final long kmer0, final long rkmer0, int[] counts, long mask, int shift2){
		assert(kmer0==rcomp(rkmer0));
		if(verbose){outstream.println("fillLeftCounts:    "+toText(kmer0)+",   "+toText(rkmer0));}
		long rkmer=(rkmer0<<2)&mask;
		long kmer=(kmer0>>>2);
		int max=-1, maxPos=0;
//		assert(false) : shift2+", "+k;
		for(int i=0; i<=3; i++){
			final long rkmer2=rkmer|((long)AminoAcid.numberToComplement[i]);
			final long kmer2=kmer|(((long)i)<<shift2);
			if(verbose){outstream.println("kmer:             "+toText(kmer2)+",     "+toText(rkmer2));}
			assert(kmer2==(kmer2&mask));
			assert(rkmer2==(rkmer2&mask));
			assert(kmer2==rcomp(rkmer2)) : "\n"+"kmer:      \t"+toText(rcomp(rkmer2))+", "+toText(rcomp(kmer2));
			//KMR023: forward first preserves strand-specific lookup; canonical selection remains symmetric.
			long key=toValue(kmer2, rkmer2);
			int way=kmerToWay(key);
			int count=tables[way].getValue(key);
			assert(count==NOT_PRESENT || count>=0);
			count=Tools.max(count, 0);
			counts[i]=count;
			if(count>max){
				max=count;
				maxPos=i;
			}
		}
		return maxPos;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Printing Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Emits formatted records sequentially and finishes the writer on normal traversal completion.
	 * @return Writer error state (true means error), or false for a null filename */
	@Override
	public boolean dumpKmersAsBytes(String fname, int mincount, int maxcount, boolean printTime, AtomicLong remaining){
		if(fname==null){return false;}
		Timer t=new Timer();

		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, append, true);
		bsw.start();
		for(AbstractKmerTable set : tables){
			set.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);
		}
		bsw.poisonAndWait();

		t.stop();
		if(printTime){outstream.println("Kmer Dump Time:             \t"+t);}
		return bsw.errorState;
	}

	/** Uses dump workers when at least three shards/threads are available and NUM_THREADS is not one.
	 * Otherwise delegates to the sequential path. Joins workers and finishes the writer.
	 * @return True for reported worker/writer error, false for success or a null filename */
	@Override
	public boolean dumpKmersAsBytes_MT(String fname, int mincount, int maxcount, boolean printTime, AtomicLong remaining){
		final int threads=Tools.min(Shared.threads(), tables.length);
		if(threads<3 || DumpThread.NUM_THREADS==1){return dumpKmersAsBytes(fname, mincount, maxcount, printTime, remaining);}

		if(fname==null){return false;}
		Timer t=new Timer();

		ByteStreamWriter bsw=new ByteStreamWriter(fname, overwrite, append, true);
		bsw.start();
		final boolean success=DumpThread.dump(k, mincount, maxcount, tables, bsw, remaining);
		//KMR118: finish the writer before combining independent worker and writer failure states.
		bsw.poisonAndWait();

		t.stop();
		if(printTime){outstream.println("Kmer Dump Time:             \t"+t);}
		return !success || bsw.errorState;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Recall Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Computes reverse complement of a k-mer using fast binary operations */
	private final long rcomp(long kmer){return AminoAcid.reverseComplementBinaryFast(kmer, k);}
	/** Converts k-mer to readable text representation */
	private final StringBuilder toText(long kmer){return AbstractKmerTable.toText(kmer, k);}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Selects a stored key using the configured orientation and masks. With rcomp off,
	 * returns the middle-masked forward argument. Otherwise compares combined-mask
	 * cores, breaking ties with the larger middle-masked full key.
	 * @param kmer Forward kmer
	 * @param rkmer Reverse kmer
	 * @return Canonical value
	 */
	public final long toValue(long kmer, long rkmer){
//		long value=(rcomp ? Tools.max(kmer, rkmer) : kmer);
//		return value;
		if(!rcomp){return kmer&middleMask;}
		final long a=kmer&cmMask, b=rkmer&cmMask;
		kmer&=middleMask;//Usually -1, so it does nothing.
		rkmer&=middleMask;
		return a>b ? kmer : a<b ? rkmer : Tools.max(kmer, rkmer);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the configured symbol count k. */
	@Override
	public int kbig(){return k;}
	/** Returns the planned byte budget for the even/odd prefilter pass. */
	@Override
	public long filterMemory(int pass){return ((pass&1)==0) ? filterMemory0 : filterMemory1;}
	/** Returns whether overlap correction was requested. */
	@Override
	public boolean ecco(){return ecco;}
	/** Returns whether left quality trimming was requested. */
	@Override
	public boolean qtrimLeft(){return qtrimLeft;}
	/** Returns whether right quality trimming was requested. */
	@Override
	public boolean qtrimRight(){return qtrimRight;}
	/** Returns the pre-trimming average-quality threshold. */
	@Override
	public float minAvgQuality(){return minAvgQuality;}
	/** Returns the constructor table-budget estimate in bytes. */
	@Override
	public long tableMemory(){return tableMemory;}
	/** Returns the constructor capacity estimate, not a measured live limit. */
	@Override
	public long estimatedKmerCapacity(){return estimatedKmerCapacity;}
	/** Returns the number of table shards, independently of worker count. */
	@Override
	public int ways(){return ways;}
	/** Returns whether key selection considers reverse complements. */
	@Override
	public boolean rcomp(){return rcomp;}

	/**
	 * Maps a k-mer to its assigned table partition using modular arithmetic.
	 * Uses core mask to ensure consistent distribution.
	 * @param kmer K-mer value to map
	 * @return Table index (0 to ways-1)
	 */
	public final int kmerToWay(final long kmer){
		final int way=(int)((kmer&cmMask)%ways);
		return way;
	}

	/** Live shards selected by {@code (key&cmMask)%ways}; null before allocation or after clear. */
	private AbstractKmerTable[] tables;

	/** Returns live shard storage; callers must not replace elements during active work. */
	public AbstractKmerTable[] tables(){return tables;}

	/** Manual override for prefilter memory allocation in bytes */
	public long filterMemoryOverride=0;

	/** Hash table implementation type (e.g., ARRAY1D) */
	public final int tableType; //AbstractKmerTable.ARRAY1D;

	/** Assumed bytes per key for capacity estimates; not measured heap use. */
	private final int bytesPerKmer;

	/** Parsed schedule lower bound, retained before initialSize is replaced with a reporting fallback. */
	private final int requestedInitialSize;

	/** Constructor heuristic for usable heap budget after reserves. */
	private final long usableMemory;
	/** Planned byte budget for even-numbered prefilter passes. */
	private final long filterMemory0;
	/** Planned byte budget for odd-numbered prefilter passes. */
	private final long filterMemory1;
	/** Constructor estimate of table byte budget, separate from actual allocation. */
	private final long tableMemory;
	/** Estimated maximum number of k-mers that can be stored */
	private final long estimatedKmerCapacity;

	/** Whether allocation starts from the requested preallocated fraction. */
	private final boolean prealloc;

	/** Number of table shards; loader thread count is selected separately. */
	public final int ways;

	/** Normal kmer length */
	public final int k;
	/** k-1; used in some expressions */
	public final int k2;

	/** Bit mask for core k-mer region excluding terminal bases */
	public final long coreMask;
	/** Bit mask for middle region of k-mers if enabled */
	public final long middleMask;
	/** Combined core and middle mask for k-mer processing */
	public final long cmMask;

	/** Look for reverse-complements as well as forward kmers.  Default: true */
	private final boolean rcomp;

	/** Quality-trim the left side */
	public final boolean qtrimLeft;
	/** Quality-trim the right side */
	public final boolean qtrimRight;
	/** Trim bases at this quality or below.  Default: 4 */
	public final float trimq;
	/** Error rate for trimming (derived from trimq) */
	private final float trimE;

	/** Minimum average quality checked before trimming; zero disables this filter. */
	public final float minAvgQuality;
	/** If positive, calculate average quality from the first X bases only.  Default: 0 */
	public final int minAvgQualityBases;

	/** Ignore kmers with probability of correctness less than this */
	public final float minProb;

	/** Correct via overlap */
	private final boolean ecco;

	/** Attempt to merge via overlap prior to counting kmers */
	private final boolean merge;

	/*--------------------------------------------------------------*/
	/*----------------            Walker            ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a Walker for iterating through all k-mers across all tables.
	 * @return A walker borrowing live shards; each table must support walking and structure must remain stable */
	public Walker walk(){
		return new WalkerKST();
	}

	/** Iterator for traversing all k-mers across all tables in the KmerTableSet.
	 * Advances through tables sequentially when current table is exhausted.
	 * Read key/value only after next() returns true; after exhaustion the delegate is null. */
	public class WalkerKST extends Walker{

		/** Starts at shard zero; requires a nonempty allocated table array. */
		WalkerKST(){
			w=tables[0].walk();
		}

		/**
		 * Fills this object with the next key and value.
		 * @return True if successful.
		 */
		public boolean next(){
			while(w!=null){
				if(w.next()){return true;}
				tnum++;
				w=(tnum<tables.length ? tables[tnum].walk() : null);
			}
			return false;
		}

		/** Returns the current k-mer from the active table walker */
		public long kmer(){return w.kmer();}
		/** Returns the current k-mer count from the active table walker */
		public int value(){return w.value();}

		/** Current table walker being used for iteration */
		private Walker w=null;

		/** current table number */
		private int tnum=0;
	}

}
