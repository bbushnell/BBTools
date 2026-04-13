package cardinality;

import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import parse.Parse;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import rand.FastRandomXoshiro;
import structures.ByteBuilder;
import structures.IntHashMap;
import template.Accumulator;
import template.ThreadWaiter;

/**
 * DLL4 Word-Based Simulator.
 *
 * Simulates DLL4 words to build per-tier state tables and/or validate
 * estimation accuracy against a loaded table.
 *
 * A "word" is a group of 4 DLL4 buckets (4 bits each = 16 bits total).
 * Each word's state is canonicalized: subtract the minimum register value
 * from all 4, sort ascending.  This always yields a form starting with 0,
 * giving C(18,3)=816 equivalence classes.  The absolute tier is
 * globalExp + wordMin, making canonical states comparable regardless of
 * how many words share the global floor.
 *
 * Recording modes control when observations are added to the state table:
 *   every - record state at every cardinality increment (original behavior).
 *           Biased toward long-lived states since they accumulate more samples.
 *   entry - record only when a word transitions to a new state.
 *           Each state visit contributes exactly one observation at the entry
 *           cardinality.  Produces lower stateAvg values than 'every'.
 *   both  - record at transitions: new state at trueCard (entry), AND
 *           old state at trueCard-1 (exit).  Averages entry and exit.
 *
 * Error tracking (enabled by table= flag):
 *   Loads a precomputed state table, simulates a full N-word DLL4,
 *   and at each state transition recomputes the running estimate
 *   (sum of per-word stateAvg lookups across all words).  Records
 *   absolute and signed error vs trueCard, both globally and per-tier.
 *
 * Run: java -ea cardinality.DLL4WordSimulator [iters=N] [threads=N]
 *      [maxTier=N] [words=N] [mode=every|entry|both] [table=file]
 *
 * @author Brian Bushnell, Chloe
 * @date April 10, 2026
 */
public class DLL4WordSimulator implements Accumulator<DLL4WordSimulator.SimThread> {

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	static final int BUCKETS_PER_WORD=4;// DLL4 packs 4 buckets per 16-bit word

	static final int MODE_EVERY=0;// Record at every cardinality step
	static final int MODE_ENTRY=1;// Record only on state transitions (new state)
	static final int MODE_BOTH=2; // Record on transitions (new + old state)

	static final int AVG_LINEAR=0;   // Arithmetic mean (original)
	static final int AVG_GEOMETRIC=1;// Geometric mean: exp(mean(log(card)))
	static final int AVG_HARMONIC=2; // Harmonic mean: n / sum(1/card)
	static final int AVG_BLEND=3;    // Weighted blend: (geoWeight*geo + harmWeight*harm) / (geoWeight+harmWeight)

	/*--------------------------------------------------------------*/
	/*----------------            Remap             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Maps each of the 65536 possible raw 16-bit register states to one of
	 * 816 canonical indices.  Canonicalization: subtract the minimum of the
	 * 4 nibbles from all 4 (so at least one is 0), then sort ascending.
	 * This makes {3,5,3,7} and {7,3,5,3} map to the same canonical form
	 * {0,0,2,4}, since order and absolute floor don't matter for estimation.
	 * Built once at class load; used by both simulator and evaluator.
	 */
	static final char[] REMAP=new char[65536];
	static final int NUM_CANONICAL;// 816 = C(18,3)
	/** Inverse: canonical index -> canonical packed key (for printing). */
	static final int[] CANONICAL_KEYS;

	static{
		IntHashMap map=new IntHashMap(1024);// canonKey -> compact index
		for(int state=0; state<65536; state++){
			int canonKey=canonize(state);
			assert (canonKey&0xF)==0 : "Canonical form must start with 0";
			if(!map.containsKey(canonKey)){
				map.put(canonKey, map.size());// Assign next compact index
			}
			REMAP[state]=(char)map.get(canonKey);// Store mapping
		}
		NUM_CANONICAL=map.size();
		assert(NUM_CANONICAL==816);// C(18,3)=816 equivalence classes
		CANONICAL_KEYS=new int[NUM_CANONICAL];
		int[] keys=map.keys();
		for(int k : keys){
			if(k==map.invalid()){continue;}// Skip IntHashMap sentinel values
			CANONICAL_KEYS[map.get(k)]=k;// Build inverse map for printing
		}
	}

	/** Canonicalizes a raw 16-bit register state: subtract min, sort ascending.
	 *  Returns a packed nibble key where the lowest nibble is always 0. */
	static final int canonize(int state){
		int r0=state&0xF;
		int r1=(state>>4)&0xF;
		int r2=(state>>8)&0xF;
		int r3=(state>>12)&0xF;
		int wMin=Math.min(Math.min(r0, r1), Math.min(r2, r3));
		// Subtract minimum so canonical form starts with 0
		int d0=r0-wMin, d1=r1-wMin, d2=r2-wMin, d3=r3-wMin;
		// Sort ascending via 5-comparator sorting network
		if(d0>d1){int t=d0; d0=d1; d1=t;}
		if(d2>d3){int t=d2; d2=d3; d3=t;}
		if(d0>d2){int t=d0; d0=d2; d2=t;}
		if(d1>d3){int t=d1; d1=d3; d3=t;}
		if(d1>d2){int t=d1; d1=d2; d2=t;}
		return d0|(d1<<4)|(d2<<8)|(d3<<12);// Pack sorted nibbles
	}

	/** Formats a packed canonical key as {a,b,c,d} for human-readable output. */
	static String keyToString(int key){
		return "{"+(key&0xF)+","+((key>>4)&0xF)+","
			+((key>>8)&0xF)+","+((key>>12)&0xF)+"}";
	}

	/** Thomas Wang's 64-bit hash finalizer.  Ensures uniform distribution
	 *  of NLZ values from sequential input keys. */
	static long hash64shift(long key){
		key=(~key)+(key<<21);
		key^=(key>>>24);
		key+=(key<<3)+(key<<8);
		key^=(key>>>14);
		key+=(key<<2)+(key<<4);
		key^=(key>>>28);
		key+=(key<<31);
		return key;
	}

	/** Packs 4 register values into a 16-bit raw key for REMAP lookup. */
	static int packWord(int[] wr){
		return wr[0]|(wr[1]<<4)|(wr[2]<<8)|(wr[3]<<12);
	}

	/** Returns the minimum of 4 register values (the word's local floor). */
	static int wordMin(int[] wr){
		return Math.min(Math.min(wr[0], wr[1]), Math.min(wr[2], wr[3]));
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	final int iters;   // Number of independent simulation runs
	final int threads; // Parallelism for simulation
	final int maxTier; // Maximum tier to record; tiers beyond this are sparse/truncated
	final int tiers;   // maxTier+2; array dimension for accumulation
	final int words;   // Words per simulated DLL4 (1=table building, 512=2048-bucket validation)
	final int mode;    // MODE_EVERY, MODE_ENTRY, or MODE_BOTH
	final int avgMode;// AVG_LINEAR, AVG_GEOMETRIC, or AVG_HARMONIC

	// Table-building accumulators [tier][canonicalIdx].
	// counts=number of observations; sums=sum of transformed values;
	// sumSq=sum of squared transformed values (for CV calculation).
	// sums2=second accumulator for blend mode (stores 1/card while sums stores log(card)).
	long[][] counts;
	double[][] sums;
	double[][] sumSq;
	double[][] sums2;

	// Blend weights for AVG_BLEND mode.
	double geoWeight=1, harmWeight=2;

	// Loaded state table for error tracking (null=table-building mode).
	// loadedTable[tier][canonicalIdx]=average per-word cardinality for that state.
	double[][] loadedTable;

	// Global error tracking accumulators.
	double logErrSum;     // sum of per-transition |est-true|/true
	double logSignedSum;  // sum of per-transition (est-true)/true
	double linErrNum;     // sum of |est-true| (for linear-weighted average)
	double linSignedNum;  // sum of (est-true)
	double linErrDen;     // sum of trueCard (denominator for linear averages)
	long transitionCount; // total transitions with error recorded

	// Per-tier error tracking: same metrics binned by transitioning word's tier.
	double[] tierLogErr;
	double[] tierLogSigned;
	double[] tierLinNum;
	double[] tierLinSigned;
	double[] tierLinDen;
	long[] tierTransitions;

	long belowFloorCount;// Hashes that fell below the early-exit mask or floor
	long totalAdds;      // Total hash attempts across all iterations
	long advanceCount;   // Number of globalExp advancement events

	boolean errorState=false;
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public DLL4WordSimulator(int iters, int threads, int maxTier, int words, int mode, int avgMode){
		this.iters=iters;
		this.threads=threads;
		this.maxTier=maxTier;
		this.tiers=maxTier+2;// Extra tier for overflow observations
		this.words=words;
		this.mode=mode;
		this.avgMode=avgMode;
		this.counts=new long[tiers][NUM_CANONICAL];
		this.sums=new double[tiers][NUM_CANONICAL];
		this.sumSq=new double[tiers][NUM_CANONICAL];
		this.sums2=(avgMode==AVG_BLEND) ? new double[tiers][NUM_CANONICAL] : null;
		this.tierLogErr=new double[tiers];
		this.tierLogSigned=new double[tiers];
		this.tierLinNum=new double[tiers];
		this.tierLinSigned=new double[tiers];
		this.tierLinDen=new double[tiers];
		this.tierTransitions=new long[tiers];
	}

	/*--------------------------------------------------------------*/
	/*----------------         Table Loading        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads a precomputed state table from a sparse TSV file.
	 * Format: #tier \t idx \t canonKey \t stateAvg \t count \t cv.
	 * Two-pass: first finds max tier to size the array, then fills values.
	 * Supports .gz files via ByteFile auto-detection.
	 * @param fname Path to the state table TSV (optionally gzipped).
	 */
	void loadTable(String fname){
		// First pass: determine array dimensions
		ByteFile bf=ByteFile.makeByteFile(fname, false);
		int maxLoadedTier=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length<1 || line[0]=='#'){continue;}
			String s=new String(line);
			String[] parts=s.split("\t");
			int tier=Integer.parseInt(parts[0]);
			if(tier>maxLoadedTier){maxLoadedTier=tier;}
		}
		bf.close();

		final int tableTiers=maxLoadedTier+1;
		double[][] table=new double[tableTiers][NUM_CANONICAL];

		// Second pass: populate the table
		bf=ByteFile.makeByteFile(fname, false);
		int entries=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length<1 || line[0]=='#'){continue;}
			String s=new String(line);
			String[] parts=s.split("\t");
			final int tier=Integer.parseInt(parts[0]);
			final int idx=Integer.parseInt(parts[1]);
			final double avg=Double.parseDouble(parts[3]);// stateAvg column
			if(idx<NUM_CANONICAL){
				table[tier][idx]=avg;
				entries++;
			}
		}
		bf.close();
		this.loadedTable=table;
		System.err.println("Loaded state table: "+entries+" entries, "+tableTiers+" tiers from "+fname);
	}

	/**
	 * Looks up the per-word cardinality estimate for a single word from the loaded table.
	 * Computes the absolute tier (globalExp+wordMin) and canonical index via REMAP,
	 * then returns the stateAvg.  Returns 0 for empty words at globalExp=0.
	 */
	double getWordEstimate(int[] reg, int globalExp){
		if(reg[0]==0 && reg[1]==0 && reg[2]==0 && reg[3]==0 && globalExp==0){return 0;}
		final int tier=globalExp+wordMin(reg);
		final int idx=REMAP[packWord(reg)];
		if(tier<loadedTable.length){return loadedTable[tier][idx];}
		return loadedTable[loadedTable.length-1][idx];// Beyond table: use last tier
	}

	/**
	 * Recomputes the full estimate from scratch by summing all per-word contributions.
	 * Called after globalExp advances (which changes all words simultaneously).
	 * @return Total cardinality estimate.
	 */
	double computeFullEstimate(int[][] regs, int globalExp, double[] wordContrib){
		double est=0;
		for(int w=0; w<words; w++){
			wordContrib[w]=getWordEstimate(regs[w], globalExp);
			est+=wordContrib[w];
		}
		return est;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/

	/** Spawns SimThreads and waits for completion.  Results are merged
	 *  via the Accumulator interface. */
	void simulate(){
		final int base=iters/threads;
		final int rem=iters%threads;
		ArrayList<SimThread> list=new ArrayList<>(threads);
		for(int t=0; t<threads; t++){
			list.add(new SimThread(this, t, base+(t<rem ? 1 : 0)));
		}
		boolean success=ThreadWaiter.startAndWait(list, this);
		errorState|=!success;
	}

	/** Merges per-thread accumulators into instance fields. */
	@Override
	public void accumulate(SimThread st){
		for(int tier=0; tier<tiers; tier++){
			for(int s=0; s<NUM_CANONICAL; s++){
				counts[tier][s]+=st.lCounts[tier][s];
				sums[tier][s]+=st.lSums[tier][s];
				sumSq[tier][s]+=st.lSumSq[tier][s];
				if(sums2!=null && st.lSums2!=null){sums2[tier][s]+=st.lSums2[tier][s];}
			}
		}
		belowFloorCount+=st.lBelow;
		totalAdds+=st.lTotal;
		advanceCount+=st.lAdvance;
		logErrSum+=st.lLogErr;
		logSignedSum+=st.lLogSigned;
		linErrNum+=st.lLinNum;
		linSignedNum+=st.lLinSigned;
		linErrDen+=st.lLinDen;
		transitionCount+=st.lTransitions;
		for(int tier=0; tier<tiers; tier++){
			tierLogErr[tier]+=st.lTierLogErr[tier];
			tierLogSigned[tier]+=st.lTierLogSigned[tier];
			tierLinNum[tier]+=st.lTierLinNum[tier];
			tierLinSigned[tier]+=st.lTierLinSigned[tier];
			tierLinDen[tier]+=st.lTierLinDen[tier];
			tierTransitions[tier]+=st.lTierTrans[tier];
		}
		errorState|=st.errorState;
	}

	@Override
	public ReadWriteLock rwlock(){return rwlock;}

	@Override
	public boolean success(){return !errorState;}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Class          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Worker thread that runs a subset of simulation iterations.
	 * Each thread independently simulates inserting up to 10M unique elements
	 * into a (words*4)-bucket DLL4, recording canonical state observations
	 * and optionally tracking estimation error against a loaded table.
	 * <p>
	 * Per-element state (wordIdx, oldTier, newTier, changed, etc.) is held
	 * as instance fields so that small methods can communicate without
	 * passing many parameters.
	 */
	static class SimThread extends Thread {

		/*--------------------------------------------------------------*/
		/*----------------        Thread Init           ----------------*/
		/*--------------------------------------------------------------*/

		SimThread(DLL4WordSimulator sim_, int tid_, int myIters_){
			sim=sim_;
			tid=tid_;
			myIters=myIters_;
			rng=new FastRandomXoshiro(tid+1);// Deterministic per-thread seed
			totalBuckets=BUCKETS_PER_WORD*sim.words;
			trackErrors=(sim.loadedTable!=null);
			singleWord=(sim.words==1);
			// Allocate per-thread accumulators
			lCounts=new long[sim.tiers][NUM_CANONICAL];
			lSums=new double[sim.tiers][NUM_CANONICAL];
			lSumSq=new double[sim.tiers][NUM_CANONICAL];
			lSums2=(sim.avgMode==AVG_BLEND) ? new double[sim.tiers][NUM_CANONICAL] : null;
			lTierLogErr=new double[sim.tiers];
			lTierLogSigned=new double[sim.tiers];
			lTierLinNum=new double[sim.tiers];
			lTierLinSigned=new double[sim.tiers];
			lTierLinDen=new double[sim.tiers];
			lTierTrans=new long[sim.tiers];
		}

		/*--------------------------------------------------------------*/
		/*----------------        Thread Body           ----------------*/
		/*--------------------------------------------------------------*/

		@Override
		public void run(){
			for(int iter=0; iter<myIters; iter++){
				simulateOneRun();
			}
			success=true;
		}

		/** Simulates one complete DLL4 lifecycle: insert up to 10M elements,
		 *  recording observations and errors at each step. */
		void simulateOneRun(){
			initRun();
			for(long trueCard=1; trueCard<=10_000_000L; trueCard++){
				processOneElement(trueCard);
				if(singleWord && newTier>sim.maxTier){break;}
			}
		}

		/** Resets per-iteration state for a fresh simulation run. */
		void initRun(){
			regs=new int[sim.words][4];
			globalExp=0;
			totalZeros=totalBuckets;
			eeMask=-1L;
			wordContrib=trackErrors ? new double[sim.words] : null;
			runningEst=0;
		}

		/** Processes one unique element: hash it, attempt a register update,
		 *  record observations, and optionally track estimation error. */
		void processOneElement(long trueCard){
			final long hash=hash64shift(rng.nextLong());
			lTotal++;
			changed=false;
			newTier=-1;

			if(!attemptUpdate(hash)){
				lBelow++;
			}
			recordObservation(trueCard);
			if(trackErrors && changed){
				trackError(trueCard);
			}
		}

		/*--------------------------------------------------------------*/
		/*----------------       Register Update        ----------------*/
		/*--------------------------------------------------------------*/

		/** Attempts to update a register from the given hash.
		 *  @return true if hash was above floor (not skipped), false if below. */
		boolean attemptUpdate(long hash){
			// Early-exit: unsigned compare filters hashes below floor threshold
			if(Long.compareUnsigned(hash, eeMask)>0){return false;}
			final int hashNlz=Long.numberOfLeadingZeros(hash);
			final int bucket=(int)((hash&0x7FFF_FFFFL)%totalBuckets);
			wordIdx=bucket/BUCKETS_PER_WORD;
			final int localBucket=bucket%BUCKETS_PER_WORD;
			final int relNlz=hashNlz-globalExp;// NLZ relative to current floor
			if(relNlz<0){return false;}

			final int newStored=Math.min(relNlz+1, 15);// DLL4 stores relNlz+1, capped at 15
			final int oldStored=regs[wordIdx][localBucket];
			if(newStored>oldStored){
				applyUpdate(localBucket, newStored, oldStored);
			}
			return true;// Above floor (even if no register change)
		}

		/** Applies a register update: snapshots old state, writes new value,
		 *  advances floor if needed, computes new state, and updates estimate. */
		void applyUpdate(int localBucket, int newStored, int oldStored){
			saveOldState();
			regs[wordIdx][localBucket]=newStored;

			final int prevGlobalExp=globalExp;
			if(oldStored==0){advanceFloor();}// Only check when a zero register becomes nonzero

			computeNewState();
			changed=(oldTier!=newTier || oldIdx!=newIdx);

			if(trackErrors && changed){
				updateEstimate(prevGlobalExp);
			}
		}

		/** Snapshots the canonical state of the affected word before modification. */
		void saveOldState(){
			final int[] wr=regs[wordIdx];
			oldIdx=REMAP[packWord(wr)];
			oldTier=globalExp+wordMin(wr);
		}

		/** Computes the canonical state of the affected word after modification.
		 *  Must be called after register write and floor advancement. */
		void computeNewState(){
			final int[] wr=regs[wordIdx];// Re-read; advancement may have decremented values
			newIdx=REMAP[packWord(wr)];
			newTier=globalExp+wordMin(wr);
		}

		/** Advances globalExp when all registers across all words are nonzero.
		 *  Decrements all registers by 1 per advance step, potentially cascading
		 *  if the decrement creates no new zeros. */
		void advanceFloor(){
			totalZeros--;
			while(totalZeros==0 && globalExp<64){
				globalExp++;
				eeMask>>>=1;// Tighten early-exit filter
				totalZeros=0;
				for(int w=0; w<sim.words; w++){
					for(int b=0; b<4; b++){
						regs[w][b]=Math.max(0, regs[w][b]-1);
						if(regs[w][b]==0){totalZeros++;}
					}
				}
				lAdvance++;
			}
		}

		/** Updates the running estimate after a word transitions.
		 *  If globalExp changed, all words shifted so we recompute from scratch.
		 *  Otherwise only the affected word's contribution is updated (O(1)). */
		void updateEstimate(int prevGlobalExp){
			if(globalExp!=prevGlobalExp){
				runningEst=sim.computeFullEstimate(regs, globalExp, wordContrib);
			}else{
				runningEst-=wordContrib[wordIdx];
				wordContrib[wordIdx]=sim.getWordEstimate(regs[wordIdx], globalExp);
				runningEst+=wordContrib[wordIdx];
			}
		}

		/*--------------------------------------------------------------*/
		/*----------------         Recording           ----------------*/
		/*--------------------------------------------------------------*/

		/** Records observations based on the current recording mode.
		 *  In 'every' mode, records at every cardinality step (single-word only).
		 *  In 'entry'/'both' modes, records only when the word's state changed. */
		void recordObservation(long trueCard){
			final double recordCard=(double)trueCard/sim.words;// Per-word effective cardinality
			if(sim.mode==MODE_EVERY){
				if(singleWord){recordCurrentState(regs[0], recordCard);}
			}else if(changed){
				recordEntry(recordCard);
				if(sim.mode==MODE_BOTH){recordExit(trueCard);}
			}
		}

		/** Transforms a cardinality value for the primary accumulator.
		 *  Linear: card.  Geometric/Blend: log(card).  Harmonic: 1/card. */
		double transform(double card){
			switch(sim.avgMode){
				case AVG_GEOMETRIC: case AVG_BLEND: return Math.log(card);
				case AVG_HARMONIC:  return 1.0/card;
				default:            return card;
			}
		}

		/** Accumulates an observation into the per-thread arrays.
		 *  In blend mode, also accumulates 1/card into lSums2. */
		void accumObs(int tier, int idx, double card){
			final double v=transform(card);
			lCounts[tier][idx]++;
			lSums[tier][idx]+=v;
			lSumSq[tier][idx]+=v*v;
			if(lSums2!=null){lSums2[tier][idx]+=1.0/card;}
		}

		/** Records the current state of a word (used in MODE_EVERY). */
		void recordCurrentState(int[] wr, double card){
			final int tier=globalExp+wordMin(wr);
			final int idx=REMAP[packWord(wr)];
			if(tier<sim.tiers){accumObs(tier, idx, card);}
		}

		/** Records the new state at entry cardinality (MODE_ENTRY and MODE_BOTH). */
		void recordEntry(double card){
			if(newTier<sim.tiers){accumObs(newTier, newIdx, card);}
		}

		/** Records the old state at exit cardinality (MODE_BOTH only).
		 *  Exit cardinality is (trueCard-1)/words: the last card the old state was valid. */
		void recordExit(long trueCard){
			final double exitCard=(double)(trueCard-1)/sim.words;
			if(exitCard>0 && oldTier<sim.tiers){accumObs(oldTier, oldIdx, exitCard);}
		}

		/*--------------------------------------------------------------*/
		/*----------------       Error Tracking        ----------------*/
		/*--------------------------------------------------------------*/

		/** Records estimation error at a transition point.
		 *  Computes absolute and signed relative error, accumulates globally
		 *  and per-tier (binned by the transitioning word's tier). */
		void trackError(long trueCard){
			if(trueCard<=10 || runningEst<=0){return;}// Skip very low cardinalities (noisy)
			final double diff=runningEst-trueCard;
			final double relErr=Math.abs(diff)/trueCard;
			final double signedRel=diff/trueCard;
			// Global accumulators
			lLogErr+=relErr;
			lLogSigned+=signedRel;
			lLinNum+=Math.abs(diff);
			lLinSigned+=diff;
			lLinDen+=trueCard;
			lTransitions++;
			// Per-tier accumulators
			if(newTier>=0 && newTier<sim.tiers){
				lTierLogErr[newTier]+=relErr;
				lTierLogSigned[newTier]+=signedRel;
				lTierLinNum[newTier]+=Math.abs(diff);
				lTierLinSigned[newTier]+=diff;
				lTierLinDen[newTier]+=trueCard;
				lTierTrans[newTier]++;
			}
		}

		/*--------------------------------------------------------------*/
		/*----------------        Thread Fields         ----------------*/
		/*--------------------------------------------------------------*/

		final DLL4WordSimulator sim;// Parent simulator (read-only access to config)
		final int tid;
		final int myIters;
		final FastRandomXoshiro rng;
		final int totalBuckets;    // BUCKETS_PER_WORD * words
		final boolean trackErrors; // True when a table is loaded
		final boolean singleWord;  // True when words==1

		// Per-iteration mutable state
		int[][] regs;       // Register values: regs[wordIdx][0..3]
		int globalExp;      // Shared floor exponent (minZeros equivalent)
		int totalZeros;     // Count of zero-valued registers across ALL words
		long eeMask;        // Early-exit mask: hashes above this are below floor
		double[] wordContrib;// Per-word contribution cache for incremental estimate updates
		double runningEst;  // Current total estimate (sum of all wordContrib)

		// Per-element mutable state (set by attemptUpdate, read by recording/tracking)
		boolean changed;    // True if the affected word's canonical state changed
		int wordIdx;        // Index of the word targeted by the current hash
		int oldTier, oldIdx;// Canonical state before the update
		int newTier, newIdx;// Canonical state after the update

		// Per-thread accumulators (merged into parent via accumulate())
		long[][] lCounts;
		double[][] lSums;
		double[][] lSumSq;
		double[][] lSums2;// Blend mode only: stores 1/card
		long lBelow, lTotal, lAdvance;
		double lLogErr, lLogSigned, lLinNum, lLinSigned, lLinDen;
		long lTransitions;
		double[] lTierLogErr, lTierLogSigned, lTierLinNum, lTierLinSigned, lTierLinDen;
		long[] lTierTrans;

		boolean success=false;
		boolean errorState=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Output            ----------------*/
	/*--------------------------------------------------------------*/

	/** Computes mean cardinality per (tier, canonicalIdx) from accumulated sums.
	 *  Linear: sum/n.  Geometric: exp(sum/n).  Harmonic: n/sum.
	 *  Blend: (geoWeight*geo + harmWeight*harm) / (geoWeight+harmWeight). */
	double[][] buildStateAvg(){
		double[][] avg=new double[tiers][NUM_CANONICAL];
		for(int tier=0; tier<tiers; tier++){
			for(int s=0; s<NUM_CANONICAL; s++){
				final long n=counts[tier][s];
				if(n>0){
					final double rawMean=sums[tier][s]/n;
					switch(avgMode){
						case AVG_GEOMETRIC: avg[tier][s]=Math.exp(rawMean); break;
						case AVG_HARMONIC:  avg[tier][s]=1.0/rawMean; break;
						case AVG_BLEND:{
							final double geo=Math.exp(rawMean);// sums has log(card)
							final double harm=n/sums2[tier][s];// sums2 has 1/card
							avg[tier][s]=(geoWeight*geo+harmWeight*harm)/(geoWeight+harmWeight);
							break;
						}
						default: avg[tier][s]=rawMean; break;
					}
				}
			}
		}
		return avg;
	}

	/**
	 * Computes coefficient of variation per (tier, canonicalIdx).
	 * For geometric mode: sqrt(var(log)) ≈ CV for log-normal distributions.
	 * For linear/harmonic: stddev/mean in transformed space.
	 * Defaults to 1.0 (uninformative) if too few observations.
	 */
	double[][] buildCV(){
		double[][] cv=new double[tiers][NUM_CANONICAL];
		for(int tier=0; tier<tiers; tier++){
			for(int s=0; s<NUM_CANONICAL; s++){
				cv[tier][s]=1.0;// Default: uninformative prior
				if(counts[tier][s]>1){
					final double mean=sums[tier][s]/counts[tier][s];
					final double variance=sumSq[tier][s]/counts[tier][s]-mean*mean;
					if(variance>0){
						if(avgMode==AVG_GEOMETRIC){
							cv[tier][s]=Math.sqrt(variance);// SD in log space ≈ CV
						}else if(mean!=0){
							cv[tier][s]=Math.sqrt(variance)/Math.abs(mean);
						}
					}else{
						cv[tier][s]=0.001;// All observations identical
					}
				}
			}
		}
		return cv;
	}

	/** Observation-count-weighted average CV for each tier. */
	double[] buildTierCV(double[][] cv){
		double[] tierCV=new double[tiers];
		for(int tier=0; tier<tiers; tier++){
			double wSum=0, wTot=0;
			for(int s=0; s<NUM_CANONICAL; s++){
				if(counts[tier][s]>0){
					wSum+=cv[tier][s]*counts[tier][s];
					wTot+=counts[tier][s];
				}
			}
			tierCV[tier]=wTot>0 ? wSum/wTot : 0;
		}
		return tierCV;
	}

	/** Prints summary statistics to stderr and sparse TSV state table to stdout. */
	void printResults(double[][] stateAvg, double[][] cv){
		System.err.println("=== DLL4 Word Simulator Results ===");
		System.err.println("iters="+iters+"  threads="+threads+"  maxTier="+maxTier
			+"  words="+words+"  mode="+modeName(mode));
		System.err.println("numCanonical="+NUM_CANONICAL);
		System.err.println("totalAdds="+totalAdds);
		System.err.printf("belowFloorRate=%.6f%n",
			totalAdds>0 ? (double)belowFloorCount/totalAdds : 0);
		System.err.printf("advanceCount  =%d%n", advanceCount);

		if(loadedTable!=null && transitionCount>0){
			System.err.println();
			System.err.println("=== Error Tracking ===");
			System.err.printf("transitions    =%d%n", transitionCount);
			System.err.printf("logAvgErr      =%.4f%% (mean |est-true|/true)%n", 100.0*logErrSum/transitionCount);
			System.err.printf("logSignedErr   =%+.4f%% (mean (est-true)/true)%n", 100.0*logSignedSum/transitionCount);
			System.err.printf("linAvgErr      =%.4f%% (sum|est-true| / sumTrue)%n", 100.0*linErrNum/linErrDen);
			System.err.printf("linSignedErr   =%+.4f%% (sum(est-true) / sumTrue)%n", 100.0*linSignedNum/linErrDen);
		}

		final double[] tierCV=buildTierCV(cv);
		final boolean hasErrors=(loadedTable!=null && transitionCount>0);
		System.err.println();
		if(hasErrors){
			System.err.printf("%-6s %-14s %-10s %-10s %-6s %-8s %-10s %-10s%n",
				"Tier", "AvgCard", "Growth", "Obs", "St", "CV", "LogErr%", "Signed%");
		}else{
			System.err.printf("%-6s %-14s %-10s %-10s %-6s %-8s%n",
				"Tier", "AvgCard", "Growth", "Obs", "St", "CV");
		}
		double prevAvg=0;
		for(int tier=0; tier<tiers; tier++){
			long obs=0;
			double totalSum=0;
			int statesUsed=0;
			for(int s=0; s<NUM_CANONICAL; s++){
				obs+=counts[tier][s];
				totalSum+=sums[tier][s];
				if(counts[tier][s]>0){statesUsed++;}
			}
			if(obs==0){continue;}
			final double avg=totalSum/obs;
			final double growth=(prevAvg>0) ? avg/prevAvg : Double.NaN;
			if(hasErrors && tierTransitions[tier]>0){
				final double tLogPct=100.0*tierLogErr[tier]/tierTransitions[tier];
				final double tSignedPct=100.0*tierLogSigned[tier]/tierTransitions[tier];
				System.err.printf("%-6d %-14.2f %-10.4f %-10d %-6d %-8.4f %-10.2f %+-10.2f%n",
					tier, avg, growth, obs, statesUsed, tierCV[tier], tLogPct, tSignedPct);
			}else{
				System.err.printf("%-6d %-14.2f %-10.4f %-10d %-6d %-8.4f%n",
					tier, avg, growth, obs, statesUsed, tierCV[tier]);
			}
			prevAvg=avg;
		}

		// Sparse TSV output to stdout: only populated (tier, idx) pairs
		System.err.println();
		ByteStreamWriter bsw=new ByteStreamWriter("stdout", true, false, false);
		bsw.start();
		bsw.println("#tier\tidx\tcanonKey\tstateAvg\tcount\tcv");
		ByteBuilder bb=new ByteBuilder(256);
		for(int tier=0; tier<tiers; tier++){
			for(int idx=0; idx<NUM_CANONICAL; idx++){
				if(counts[tier][idx]>0){
					bb.clear();
					bb.append(tier).append('\t').append(idx).append('\t');
					bb.append(CANONICAL_KEYS[idx]).append('\t');
					bb.append(String.valueOf(stateAvg[tier][idx])).append('\t');
					bb.append(counts[tier][idx]).append('\t');
					bb.append(String.valueOf(cv[tier][idx]));
					bsw.println(bb);
				}
			}
		}
		bsw.poisonAndWait();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns human-readable name for a recording mode constant. */
	static String modeName(int m){
		switch(m){
			case MODE_EVERY: return "every";
			case MODE_ENTRY: return "entry";
			case MODE_BOTH:  return "both";
			default: return "unknown";
		}
	}

	/** Parses a recording mode string from command-line args. */
	static int parseMode(String s){
		switch(s.toLowerCase()){
			case "every": case "all": return MODE_EVERY;
			case "entry": case "in":  return MODE_ENTRY;
			case "both":  case "io":  return MODE_BOTH;
			default:
				System.err.println("Unknown mode: "+s+", using 'every'");
				return MODE_EVERY;
		}
	}

	/** Parses an averaging mode string. */
	static int parseAvgMode(String s){
		switch(s.toLowerCase()){
			case "linear":    case "lin": case "arithmetic": return AVG_LINEAR;
			case "geometric": case "geo": case "log":        return AVG_GEOMETRIC;
			case "harmonic":  case "har": case "harm":       return AVG_HARMONIC;
			case "blend":     case "mix":                    return AVG_BLEND;
			default:
				System.err.println("Unknown avg mode: "+s+", using 'linear'");
				return AVG_LINEAR;
		}
	}

	/** Returns human-readable name for an averaging mode constant. */
	static String avgModeName(int m){
		switch(m){
			case AVG_LINEAR:    return "linear";
			case AVG_GEOMETRIC: return "geometric";
			case AVG_HARMONIC:  return "harmonic";
			case AVG_BLEND:     return "blend";
			default: return "unknown";
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		int iters=10000;
		int threads=8;
		int maxTier=12;
		int words=1;
		int mode=MODE_EVERY;
		String tableFile=null;
		int avgMode=AVG_LINEAR;
		double geoW=1, harmW=2;// Default blend weights: 1:2 geometric:harmonic

		for(String arg : args){
			String[] kv=arg.split("=", 2);
			if(kv.length!=2){throw new RuntimeException("Unknown parameter '"+arg+"'");}
			switch(kv[0]){
				case "iters":     iters=Parse.parseIntKMG(kv[1]);   break;
				case "threads":   threads=Integer.parseInt(kv[1]); break;
				case "maxTier":   maxTier=Integer.parseInt(kv[1]); break;
				case "words":     words=Integer.parseInt(kv[1]);   break;
				case "mode":      mode=parseMode(kv[1]);           break;
				case "table":     tableFile=kv[1];                 break;
				case "avg":       avgMode=parseAvgMode(kv[1]);     break;
				case "geoweight": geoW=Double.parseDouble(kv[1]);  break;
				case "harmweight":harmW=Double.parseDouble(kv[1]); break;
				default: throw new RuntimeException("Unknown parameter '"+arg+"'");
			}
		}

		System.err.println("DLL4WordSimulator: iters="+iters+" threads="+threads
			+" maxTier="+maxTier+" words="+words+" mode="+modeName(mode)
			+" avg="+avgModeName(avgMode)
			+(avgMode==AVG_BLEND ? " geo:harm="+geoW+":"+harmW : "")
			+" numCanonical="+NUM_CANONICAL
			+(tableFile!=null ? " table="+tableFile : ""));

		final long t0=System.currentTimeMillis();
		DLL4WordSimulator sim=new DLL4WordSimulator(iters, threads, maxTier, words, mode, avgMode);
		if(avgMode==AVG_BLEND){sim.geoWeight=geoW; sim.harmWeight=harmW;}
		if(tableFile!=null){sim.loadTable(tableFile);}
		sim.simulate();
		final long t1=System.currentTimeMillis();
		System.err.println("Simulation time: "+(t1-t0)+" ms");

		double[][] stateAvg=sim.buildStateAvg();
		double[][] cvArr=sim.buildCV();
		sim.printResults(stateAvg, cvArr);
	}

}
