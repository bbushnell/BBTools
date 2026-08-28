package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.LongHashSet;
import map.LongIntMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FASTQ;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.LongList;
import structures.ListNum;

/**
 * Greedy set-cover kmer selection for tRNA (or any short-gene) covering sets.
 * Replaces the bash kmercountexact+bbduk descent loop with a single in-memory
 * pass: count all kmers across a pool of sequences, then iteratively select the
 * top-N by count, evict covered sequences, recount, and repeat until a coverage
 * target or kmer budget is reached.
 *
 * Two-tier ranking (Brian's design): each round takes the top 2*step kmers by
 * CURRENT count on the remaining pool, re-sorts them by ORIGINAL prevalence
 * (recorded before any eviction), and keeps the top step. This concentrates on
 * kmers shared between rare and common tRNAs rather than random sequence that
 * happens to be common among the rare stragglers.
 *
 * @author Neptune, Brian Bushnell
 */
public class CoveringSet {

	public static void main(String[] args){
		Timer t=new Timer();
		CoveringSet x=new CoveringSet(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public CoveringSet(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(a.equals("kdesign")){
				kDesign=Integer.parseInt(b);
			}else if(a.equals("step")){
				step=(int)Parse.parseKMG(b);
			}else if(a.equals("stepfraction") || a.equals("stepfrac") || a.equals("adaptivefraction")){
				stepFraction=Double.parseDouble(b);
			}else if(a.equals("minstepmult")){
				minStepMult=Double.parseDouble(b);
			}else if(a.equals("maxstepmult")){
				maxStepMult=Double.parseDouble(b);
			}else if(a.equals("step2boost")){
				step2Boost=Double.parseDouble(b);
			}else if(a.equals("maxkmers") || a.equals("budget")){
				maxKmers=(int)Parse.parseKMG(b);
			}else if(a.equals("mincov") || a.equals("mincoverage") || a.equals("target")){
				minCovFraction=Float.parseFloat(b);
			}else if(a.equals("extra") || a.equals("extra1")){
				extra=b;
			}else if(a.equals("copies")){
				copies=Integer.parseInt(b);
			}else if(a.equals("rcomp")){
				rcomp=Parse.parseBoolean(b);
			}else if(a.equals("partitions") || a.equals("numpartitions")){
				numPartitions=Integer.parseInt(b);
			}else if(a.equals("bufsize") || a.equals("partitionbuffer")){
				bufferSize=Integer.parseInt(b);
			}else if(parser.parse(arg, a, b)){
				//handled
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		in=parser.in1;
		out=parser.out1;
		if(in==null){throw new RuntimeException("Error - an input file is required.");}
		if(out==null){throw new RuntimeException("Error - an output file is required.");}

		if(kDesign<0){kDesign=k;}
		assert(kDesign>=k) : "kdesign must be >= k";
		assert(k>0 && k<=31) : "k must be 1..31";
		assert(kDesign<=31) : "kdesign must be <= 31";
		if(step<1){throw new IllegalArgumentException("step must be positive: "+step);}
		if(stepFraction<0 || stepFraction>1){
			throw new IllegalArgumentException("stepfraction must be in [0,1]: "+stepFraction);
		}
		if(minStepMult<=0 || maxStepMult<minStepMult || step2Boost<=0){
			throw new IllegalArgumentException("Require 0<minstepmult<=maxstepmult and step2boost>0: "+
					minStepMult+", "+maxStepMult+", "+step2Boost);
		}

		overwrite=parser.overwrite;
		ffin=FileFormat.testInput(in, FileFormat.FASTA, null, true, true);
		ffout=FileFormat.testOutput(out, FileFormat.FASTA, null, true, overwrite, false, false);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	void process(Timer t){
		ArrayList<byte[]> pool=loadPool();
		final int totalSeqs=pool.size();
		//numPartitions unset (<=0) means "scale with available threads" -- resolved
		//here rather than in field init since Shared.threads() reflects t= parsing
		//that happens after this object is constructed. Default is ~2x the thread
		//count (extra shards make it less likely two threads want to flush the same
		//partition at once, reducing lock contention) and forced ODD, since way()'s
		//kmer%numPartitions would otherwise alias against the structured low-bit
		//patterns of 2-bit-packed kmers -- an even modulus can route every kmer
		//ending in the same base to the same parity of partitions (Brian, 2026-08-28).
		if(numPartitions<=0){numPartitions=Tools.max(15, 2*Shared.threads())|1;}
		outstream.println("Pool: "+totalSeqs+" sequences, k="+k+
				(kDesign!=k ? " (design k="+kDesign+")" : "")+
				", step="+step+(stepFraction>0 ? ", stepfraction="+stepFraction+
				", stepbounds="+minStepMult+"-"+maxStepMult+", step2boost="+step2Boost : "")+
				", rcomp="+rcomp+", partitions="+numPartitions+", bufsize="+bufferSize);

		long totalBases=0;
		for(byte[] seq : pool){totalBases+=seq.length;}

		final long designMask=~((-1L)<<(2*kDesign));
		//totalBases is a safe (if loose) upper bound on distinct kDesign-mers -- used
		//only to pre-size the partition shards, never a correctness constraint.
		final PartitionedCounter originalCounts=countKmersMT(pool, null, kDesign, designMask, totalBases);
		outstream.println("Distinct "+kDesign+"-mers in pool: "+originalCounts.size());

		LongList selectedKmers=new LongList();
		LongHashSet selectedSet=new LongHashSet(maxKmers>0 ? maxKmers : 16384);
		boolean[] alive=new boolean[pool.size()];
		Arrays.fill(alive, true);
		int aliveCount=totalSeqs;
		int round=0;
		int currentStep=step;

		while(aliveCount>0){
			if(maxKmers>0 && selectedKmers.size>=maxKmers){break;}
			float covFrac=1f-(aliveCount/(float)totalSeqs);
			if(covFrac>=minCovFraction){break;}

			//Retain the historical near-target reduction when adaptive sizing is off.
			if(stepFraction<=0){
				final float gap=aliveCount/(float)totalSeqs;
				final float targetGap=1f-minCovFraction;
				if(gap<=2*targetGap && currentStep>Math.max(1, step/2)){
					currentStep=Math.max(1, step/2);
					outstream.println("  Step reduced to "+currentStep+" (approaching target, gap="+
							String.format("%.4f", gap)+", targetGap="+String.format("%.4f", targetGap)+")");
				}
			}
			final int roundStep=maxKmers>0 ? Tools.min(currentStep, maxKmers-selectedKmers.size) : currentStep;
			if(roundStep<=0){break;}

			//originalCounts.size() upper-bounds any later round's distinct count
			//(coverage only shrinks the alive set), so it's a safe, already-known
			//sizing hint -- no extra pool scan needed each round.
			PartitionedCounter currentCounts=countKmersMT(pool, alive, kDesign, designMask, originalCounts.size());
			if(currentCounts.isEmpty()){break;}

			final int candidateLimit=(int)Tools.min(Integer.MAX_VALUE, 2L*roundStep);
			long[] candidates=topNByCount(currentCounts, candidateLimit);
			long[] ranked=rankByOriginal(candidates, originalCounts, roundStep);

			int added=0;
			for(long kmer : ranked){
				if(selectedSet.contains(kmer)){continue;}
				selectedSet.add(kmer);
				selectedKmers.add(kmer);
				added++;
			}
			if(added==0){break;}

			int evicted=evictCoveredMT(pool, alive, selectedSet, kDesign, designMask);
			aliveCount-=evicted;
			round++;

			float cov=1f-(aliveCount/(float)totalSeqs);
			if(verbose || round%10==0){
				outstream.println("  Round "+round+": +"+added+" kmers (total "+
						selectedKmers.size+"), evicted "+evicted+
						", alive="+aliveCount+", coverage="+String.format("%.4f", cov));
			}
			if(stepFraction>0 && aliveCount>0){
				final int nextStep=adaptiveStep(added, evicted, aliveCount, round+1);
				if(verbose || nextStep!=currentStep){
					outstream.println("  Next step: "+nextStep+" (previous added="+added+
							", recovered="+evicted+", remaining="+aliveCount+")");
				}
				currentStep=nextStep;
			}
		}

		float finalCov=1f-(aliveCount/(float)totalSeqs);
		outstream.println("Selected "+selectedKmers.size+" "+kDesign+"-mers in "+round+
				" rounds, coverage="+String.format("%.6f", finalCov)+
				" ("+aliveCount+" uncovered of "+totalSeqs+")");

		ArrayList<Read> output;
		if(kDesign==k){
			output=kmersToReads(selectedKmers, k);
		}else{
			output=designToUseKmers(selectedKmers, kDesign, k);
		}
		outstream.println("Output: "+output.size()+" "+k+"-mer sequences");

		if(ffout!=null){
			//Writer, not ConcurrentReadOutputStream -- CRIS/CROS are being phased
			//out in favor of Streamer/Writer (Brian, 2026-08-26), same reasoning
			//as the loadFasta switch to Streamer above.
			Writer writer=WriterFactory.makeWriter(ffout);
			writer.start();
			writer.add(output, 0);
			errorState|=writer.poisonAndWait();
		}

		if(execPool!=null){execPool.shutdown();}

		t.stop();
		outstream.println("Time:   "+t);
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private ArrayList<byte[]> loadPool(){
		ArrayList<byte[]> pool=new ArrayList<>();
		loadFasta(in, pool);
		if(extra!=null){
			int base=pool.size();
			ArrayList<byte[]> extraSeqs=new ArrayList<>();
			loadFasta(extra, extraSeqs);
			for(int c=0; c<copies; c++){
				for(byte[] seq : extraSeqs){pool.add(seq);}
			}
			outstream.println("Added "+extraSeqs.size()+" extra sequences x"+copies+
					" copies (pool "+base+" -> "+pool.size()+")");
		}
		return pool;
	}

	/** Loads a FASTA file into memory via Streamer, not ConcurrentReadInputStream.
	 * CRIS spawns exactly one ReadThread for a single unpaired file, so the whole
	 * file's line-splitting/header-parsing/Read-construction runs on ONE core no
	 * matter how many cores are available -- confirmed pinning a single core at
	 * 99% for the entire load phase on a real 2.4Gbp input (Brian, 2026-08-26).
	 * StreamerFactory.makeStreamer(ff, pairnum, ordered, maxReads) with threads
	 * left at its default (-1, "auto") dispatches to the real parallel
	 * FastaStreamer when Shared.threads()>=8: one I/O thread splitting on '>'
	 * boundaries feeds an ordered queue that DEFAULT_THREADS=3 worker threads
	 * drain concurrently, each doing the actual Read construction/validation.
	 * (StreamerFactory.getReads(...) -- the tempting one-liner -- is NOT this:
	 * it hardcodes threads=1, which routes to FastaStreamerST, a single dedicated
	 * worker thread -- same one-core bottleneck as the old CRIS code, just a
	 * newer class name. Call makeStreamer directly to actually get the parallel
	 * path.) Below Shared.threads()<8 this transparently falls back to a
	 * single-threaded streamer, so this is never worse than the old behavior. */
	private void loadFasta(String fname, ArrayList<byte[]> dest){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		Streamer st=StreamerFactory.makeStreamer(ff, 0, true, -1);
		st.start();
		for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()){
			for(Read r : ln){
				if(r.bases!=null && r.length()>0){
					Tools.toUpperCase(r.bases);
					dest.add(r.bases);
				}
			}
		}
		st.close();
		//Reference-data load: a partial/corrupt read here silently produces a
		//wrong kmer pool with no error visible downstream, so crash loud instead
		//of continuing -- same idiom as StreamerFactory.getReads' own error path.
		if(st.errorState()){
			KillSwitch.kill("Error: a read error (corrupt or truncated input) occurred reading "
				+fname+"; aborting rather than building a covering set from partial/incorrect data.");
		}
	}

	/** Minimum pool size before bothering to fan out across threads -- below this,
	 * per-task submission overhead would exceed the work being parallelized. */
	private static final int MIN_PARALLEL_POOL=64;

	/** Worker pool for the per-round count/evict fan-out, sized by Shared.threads()
	 * and created lazily on first use. REUSED across every round of the greedy
	 * set-cover loop -- spawning fresh Threads per round cost more in thread-startup
	 * overhead than it saved once the algorithm's ~18+ rounds each paid that tax; a
	 * persistent pool amortizes it to (effectively) once per run, matching the
	 * pool()/Future pattern already used for the same reason in
	 * TrnaConsensusBuilder. */
	private ExecutorService pool(){
		if(execPool==null){execPool=Executors.newFixedThreadPool(Tools.max(1, Shared.threads()));}
		return execPool;
	}
	private ExecutorService execPool;

	/** Partitioned kmer counter (Brian's design, 2026-08-26). The previous
	 * multithreaded version gave each thread its own private LongIntMap, then
	 * folded all of them into one map at the end -- correct, but that final
	 * fold is inherently serial and, measured directly, ate roughly a third of
	 * the whole multithreaded runtime (worse the fewer times each kmer repeats,
	 * since repeats are what the fold-by-increment amortizes over). It also caps
	 * total distinct kmers at whatever one LongIntMap can hold.
	 * <p>
	 * This splits the kmer space itself across `numPartitions` independent
	 * LongIntMap shards, each guarded by its own lock. A thread never writes a
	 * shard directly -- it accumulates kmers into its own small per-partition
	 * buffer (a LocalBuffer, NOT shared with other threads) and only takes a
	 * shard's lock to flush a full batch. So lock acquisitions happen once per
	 * `bufferSize` kmers per partition, not once per kmer, keeping contention
	 * low even with many threads hammering the same run. Once every thread's
	 * range is processed and its buffers flushed, the counts are already in
	 * their final home -- there is no merge step at all. Splitting across
	 * shards also means total capacity is N separate LongIntMaps instead of
	 * one, so a corpus with more distinct kmers than a single table comfortably
	 * holds (the motivating case: RefSeq Plastid, 2.4Gbp, ~500M distinct
	 * kmers) is no longer a ceiling. */
	private final class PartitionedCounter {
		final LongIntMap[] shards;
		final Object[] locks;

		PartitionedCounter(long expectedDistinct){
			shards=new LongIntMap[numPartitions];
			locks=new Object[numPartitions];
			final int perShard=(int)Tools.max(256, expectedDistinct/numPartitions+1);
			for(int i=0; i<numPartitions; i++){
				shards[i]=new LongIntMap(perShard);
				locks[i]=new Object();
			}
		}

		private int way(long kmer){return (int)(kmer%numPartitions);}

		/** Per-thread accumulation buffer -- create one per worker via
		 * newLocalBuffer(), never share it across threads, and call flushAll()
		 * exactly once when that worker's assigned range is fully processed. */
		final class LocalBuffer {
			private final long[][] buf=new long[numPartitions][bufferSize];
			private final int[] size=new int[numPartitions];

			void add(long kmer){
				final int w=way(kmer);
				int s=size[w];
				buf[w][s]=kmer;
				s++;
				if(s>=bufferSize){flush(w, s); s=0;}
				size[w]=s;
			}

			private void flush(int w, int n){
				final long[] b=buf[w];
				final LongIntMap shard=shards[w];
				synchronized(locks[w]){
					for(int i=0; i<n; i++){shard.increment(b[i]);}
				}
			}

			void flushAll(){
				for(int w=0; w<numPartitions; w++){
					if(size[w]>0){flush(w, size[w]); size[w]=0;}
				}
			}
		}

		LocalBuffer newLocalBuffer(){return new LocalBuffer();}

		int get(long kmer){return shards[way(kmer)].get(kmer);}

		long size(){
			long sum=0;
			for(LongIntMap m : shards){sum+=m.size();}
			return sum;
		}

		boolean isEmpty(){return size()==0;}
	}

	/** Fills pool[from,to) (or only the alive ones there, if alive!=null) into a
	 * per-thread LocalBuffer. Rolls a forward kmer and (when rcomp is set) a
	 * reverse-complement kmer in lockstep and buffers only the canonical (max
	 * of the two) form -- same idiom as prok/RiboMaker.loadFilter -- so a motif
	 * and its reverse complement are counted as the SAME kmer instead of two
	 * independent ones. `buf` must be this thread's own LocalBuffer: safe to
	 * call from multiple threads concurrently as long as each call gets a
	 * distinct buffer and a disjoint [from,to) range of the same pool/alive. */
	private void fillRange(ArrayList<byte[]> pool, boolean[] alive, int from, int to,
			int kLen, long mask, PartitionedCounter.LocalBuffer buf){
		final int shift2=2*kLen-2;
		for(int i=from; i<to; i++){
			if(alive!=null && !alive[i]){continue;}
			byte[] bases=pool.get(i);
			long kmer=0, rkmer=0;
			int len=0;
			for(byte b : bases){
				int num=AminoAcid.baseToNumber[b];
				if(num>=0){
					kmer=((kmer<<2)|num)&mask;
					if(rcomp){
						long comp=AminoAcid.baseToComplementNumber[b];
						rkmer=((rkmer>>>2)|(comp<<shift2))&mask;
					}
					len++;
					if(len>=kLen){buf.add(rcomp ? Tools.max(kmer, rkmer) : kmer);}
				}else{len=0; kmer=0; rkmer=0;}
			}
		}
		buf.flushAll();
	}

	/** Multithreaded partitioned kmer counting -- see PartitionedCounter's own
	 * doc for why this replaced the previous private-map-per-thread-then-merge
	 * design. alive==null counts the whole pool (the one-time full-corpus
	 * count); a non-null alive[] counts only the currently-uncovered sequences
	 * (the per-round recount). expectedDistinct is only a SIZING HINT for the
	 * shards (never a correctness constraint) -- pass the best upper-bound
	 * estimate available; under- or over-estimating just costs a few extra
	 * resizes or a little unused capacity. Falls back to a single-threaded fill
	 * (still partitioned, so the return type and downstream code are uniform)
	 * when the pool is too small to be worth splitting across threads. */
	private PartitionedCounter countKmersMT(final ArrayList<byte[]> pool, final boolean[] alive,
			final int kLen, final long mask, long expectedDistinct){
		final PartitionedCounter counter=new PartitionedCounter(expectedDistinct);
		final int n=pool.size();
		final int threads=Tools.max(1, Shared.threads());
		if(threads<=1 || n<MIN_PARALLEL_POOL){
			fillRange(pool, alive, 0, n, kLen, mask, counter.newLocalBuffer());
			return counter;
		}
		final int chunk=(n+threads-1)/threads;
		final ArrayList<Future<?>> futures=new ArrayList<>(threads);
		for(int t=0; t<threads; t++){
			final int from=t*chunk, to=Tools.min(n, from+chunk);
			if(from>=to){continue;}
			final PartitionedCounter.LocalBuffer buf=counter.newLocalBuffer();
			futures.add(pool().submit(new Runnable(){
				@Override
				public void run(){fillRange(pool, alive, from, to, kLen, mask, buf);}
			}));
		}
		for(Future<?> f : futures){
			try{f.get();}catch(Exception e){throw new RuntimeException(e);}
		}
		return counter;
	}

	/** Returns the top N kmers by count without materializing or sorting the
	 * complete kmer space. Each shard is scanned independently in parallel into
	 * a bounded primitive heap, then the small shard results are merged. Numeric
	 * kmer order breaks count ties reproducibly, independent of shard count. */
	private long[] topNByCount(final PartitionedCounter counter, final int n){
		if(n<=0 || counter.isEmpty()){return new long[0];}
		final int limit=(int)Tools.min(n, counter.size());
		final int ways=counter.shards.length;
		final int threads=Tools.max(1, Shared.threads());
		final TopKHeap merged=new TopKHeap(limit);
		if(threads<=1 || ways<=1){
			for(int way=0; way<ways; way++){merged.add(topKFromShard(counter.shards[way], limit));}
		}else{
			final ArrayList<Future<TopKHeap>> futures=new ArrayList<>(ways);
			for(int way=0; way<ways; way++){
				final int w=way;
				futures.add(pool().submit(new Callable<TopKHeap>(){
					@Override
					public TopKHeap call(){return topKFromShard(counter.shards[w], limit);}
				}));
			}
			for(Future<TopKHeap> f : futures){
				try{merged.add(f.get());}catch(Exception e){throw new RuntimeException(e);}
			}
		}
		return merged.keysDescending();
	}

	/** Scans one immutable shard into a bounded top-K heap. */
	private static TopKHeap topKFromShard(LongIntMap map, int limit){
		final TopKHeap heap=new TopKHeap(Tools.min(limit, map.size()));
		final long[] keys=map.keys();
		final int[] counts=map.values();
		final long invalid=map.invalid();
		for(int cell=0; cell<keys.length; cell++){
			final long key=keys[cell];
			if(key!=invalid){heap.add(key, counts[cell], key);}
		}
		return heap;
	}

	/** From the candidate list, re-rank by original prevalence, keep top step. */
	private long[] rankByOriginal(long[] candidates, PartitionedCounter originalCounts, int keep){
		final TopKHeap heap=new TopKHeap(Tools.min(keep, candidates.length));
		for(int i=0; i<candidates.length; i++){
			heap.add(candidates[i], originalCounts.get(candidates[i]), candidates[i]);
		}
		return heap.keysDescending();
	}

	/** Projects the next batch size from the previous round's observed
	 * sequences-per-kmer recovery. Bounds and the optional round-2 boost are
	 * tunable because conservation changes the useful step scale. */
	private int adaptiveStep(int previousKmers, int previousRecovered, int remaining, int nextRound){
		double projected=previousKmers*(stepFraction*remaining)/Tools.max(1, previousRecovered);
		if(nextRound==2){projected*=step2Boost;}
		final long raw=projected>=Integer.MAX_VALUE ? Integer.MAX_VALUE : (long)Math.ceil(projected);
		final long min=scaledStep(minStepMult);
		final long max=Tools.max(min, scaledStep(maxStepMult));
		return (int)Tools.mid(min, max, raw);
	}

	private long scaledStep(double mult){
		final double x=Math.ceil(step*mult);
		return Tools.max(1L, x>=Integer.MAX_VALUE ? Integer.MAX_VALUE : (long)x);
	}

	/** Marks sequences as not-alive if they contain any selected (canonical) kmer,
	 * over pool[from,to). Only reads `selected` (a fixed snapshot for the whole
	 * round) and writes exclusively to its own index range of `alive`, so this is
	 * safe to run from multiple threads on disjoint ranges concurrently.
	 * Returns the number of newly evicted sequences in this range. */
	private int evictRange(ArrayList<byte[]> pool, boolean[] alive,
			LongHashSet selected, int from, int to, int kLen, long mask){
		int evicted=0;
		final int shift2=2*kLen-2;
		for(int i=from; i<to; i++){
			if(!alive[i]){continue;}
			byte[] bases=pool.get(i);
			long kmer=0, rkmer=0;
			int len=0;
			for(byte b : bases){
				int num=AminoAcid.baseToNumber[b];
				if(num>=0){
					kmer=((kmer<<2)|num)&mask;
					if(rcomp){
						long comp=AminoAcid.baseToComplementNumber[b];
						rkmer=((rkmer>>>2)|(comp<<shift2))&mask;
					}
					len++;
					if(len>=kLen && selected.contains(rcomp ? Tools.max(kmer, rkmer) : kmer)){
						alive[i]=false;
						evicted++;
						break;
					}
				}else{len=0; kmer=0; rkmer=0;}
			}
		}
		return evicted;
	}

	/** Multithreaded wrapper for evictRange -- same chunking/fallback/reused-pool
	 * strategy as countKmersMT. Threads write disjoint indices of `alive`, so no
	 * merge step is needed beyond summing the per-thread eviction counts. */
	private int evictCoveredMT(final ArrayList<byte[]> pool, final boolean[] alive,
			final LongHashSet selected, final int kLen, final long mask){
		final int n=pool.size();
		final int threads=Tools.max(1, Shared.threads());
		if(threads<=1 || n<MIN_PARALLEL_POOL){
			return evictRange(pool, alive, selected, 0, n, kLen, mask);
		}
		final int chunk=(n+threads-1)/threads;
		final ArrayList<Future<Integer>> futures=new ArrayList<>(threads);
		for(int t=0; t<threads; t++){
			final int from=t*chunk, to=Tools.min(n, from+chunk);
			if(from>=to){continue;}
			futures.add(pool().submit(new Callable<Integer>(){
				@Override
				public Integer call(){return evictRange(pool, alive, selected, from, to, kLen, mask);}
			}));
		}
		int total=0;
		for(Future<Integer> f : futures){
			try{total+=f.get();}catch(Exception e){throw new RuntimeException(e);}
		}
		return total;
	}

	/** Converts selected kmer longs into single-kmer FASTA reads. */
	private ArrayList<Read> kmersToReads(LongList kmers, int kLen){
		ArrayList<Read> list=new ArrayList<>(kmers.size);
		for(int i=0; i<kmers.size; i++){
			byte[] bases=longToBytes(kmers.get(i), kLen);
			list.add(new Read(bases, null, "kmer_"+i, i));
		}
		return list;
	}

	/** For k+1 design: each selected (k+1)-mer yields two overlapping k-mers.
	 * Deduplicate and output only the unique k-mers, preserving selection order. */
	private ArrayList<Read> designToUseKmers(LongList designKmers, int kD, int kU){
		final long useMask=~((-1L)<<(2*kU));
		LongHashSet seen=new LongHashSet(designKmers.size*2);
		LongList useKmers=new LongList();
		for(int i=0; i<designKmers.size; i++){
			long dkmer=designKmers.get(i);
			// A (kD)-mer contains (kD-kU+1) overlapping (kU)-mers
			for(int offset=0; offset<=kD-kU; offset++){
				long ukmer=(dkmer>>(2*offset))&useMask;
				if(!seen.contains(ukmer)){
					seen.add(ukmer);
					useKmers.add(ukmer);
				}
			}
		}
		outstream.println("Design "+kD+"-mers: "+designKmers.size+" -> use "+kU+"-mers: "+useKmers.size);
		ArrayList<Read> list=new ArrayList<>(useKmers.size);
		for(int i=0; i<useKmers.size; i++){
			byte[] bases=longToBytes(useKmers.get(i), kU);
			list.add(new Read(bases, null, "kmer_"+i, i));
		}
		return list;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Utility Methods       ----------------*/
	/*--------------------------------------------------------------*/

	private static byte[] longToBytes(long kmer, int kLen){
		byte[] bases=new byte[kLen];
		for(int i=kLen-1; i>=0; i--){
			bases[i]=AminoAcid.numberToBase[(int)(kmer&3)];
			kmer>>>=2;
		}
		return bases;
	}

	/** Primitive bounded min-heap. The root is the worst retained item: lower
	 * count, then later stable order. This avoids boxed entries in the hot rank. */
	private static final class TopKHeap {
		TopKHeap(int capacity_){
			capacity=capacity_;
			keys=new long[capacity];
			counts=new int[capacity];
			orders=new long[capacity];
		}

		void add(TopKHeap b){
			for(int i=0; i<b.size; i++){add(b.keys[i], b.counts[i], b.orders[i]);}
		}

		boolean add(long key, int count, long order){
			if(capacity<=0){return false;}
			if(size<capacity){
				keys[size]=key;
				counts[size]=count;
				orders[size]=order;
				siftUp(size);
				size++;
				return true;
			}
			if(!better(count, order, counts[0], orders[0])){return false;}
			keys[0]=key;
			counts[0]=count;
			orders[0]=order;
			siftDown(0);
			return true;
		}

		long[] keysDescending(){
			final long[] out=new long[size];
			for(int i=size-1; i>=0; i--){
				out[i]=keys[0];
				size--;
				if(size>0){
					keys[0]=keys[size];
					counts[0]=counts[size];
					orders[0]=orders[size];
					siftDown(0);
				}
			}
			return out;
		}

		private void siftUp(int child){
			final long key=keys[child], order=orders[child];
			final int count=counts[child];
			while(child>0){
				final int parent=(child-1)>>>1;
				if(!worse(count, order, counts[parent], orders[parent])){break;}
				keys[child]=keys[parent];
				counts[child]=counts[parent];
				orders[child]=orders[parent];
				child=parent;
			}
			keys[child]=key;
			counts[child]=count;
			orders[child]=order;
		}

		private void siftDown(int parent){
			final long key=keys[parent], order=orders[parent];
			final int count=counts[parent];
			for(int child=parent*2+1; child<size; child=parent*2+1){
				if(child+1<size && worse(counts[child+1], orders[child+1],
						counts[child], orders[child])){child++;}
				if(!worse(counts[child], orders[child], count, order)){break;}
				keys[parent]=keys[child];
				counts[parent]=counts[child];
				orders[parent]=orders[child];
				parent=child;
			}
			keys[parent]=key;
			counts[parent]=count;
			orders[parent]=order;
		}

		private static boolean better(int ac, long ao, int bc, long bo){
			return ac>bc || (ac==bc && ao<bo);
		}

		private static boolean worse(int ac, long ao, int bc, long bo){
			return ac<bc || (ac==bc && ao>bo);
		}

		private final int capacity;
		private final long[] keys;
		private final int[] counts;
		private final long[] orders;
		private int size=0;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in;
	private String out;
	private String extra;
	private final FileFormat ffin;
	private final FileFormat ffout;
	private boolean overwrite=true;
	private boolean errorState=false;

	private int k=17;
	private int kDesign=-1;
	private int step=500;
	/** Fraction of remaining sequences targeted per adaptive round; 0 keeps the
	 * historical fixed-step behavior. */
	private double stepFraction=0;
	private double minStepMult=0.1;
	private double maxStepMult=2;
	private double step2Boost=1;
	private int maxKmers=0;
	private float minCovFraction=0.999f;
	private int copies=10;
	private boolean rcomp=false;
	/** Number of PartitionedCounter shards. <=0 (default) resolves in process()
	 * to the nearest odd value >= max(15, 2*Shared.threads()) -- Brian's
	 * illustrative floor of 15, doubled against thread count so flush-time lock
	 * collisions between threads stay rare, and forced odd so kmer%numPartitions
	 * doesn't alias against 2-bit-packed kmers' structured low bits. */
	private int numPartitions=-1;
	/** Per-thread, per-partition batch size before a flush takes that
	 * partition's lock (Brian's illustrative 200). */
	private int bufferSize=200;

	private PrintStream outstream=System.err;
	private static boolean verbose=false;
}
