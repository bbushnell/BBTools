package test.bloom;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Locale;
import java.util.Random;

import assemble.ErrorTracker;
import bloom.BloomFilter;
import bloom.BloomFilterCorrector;
import bloom.BloomFilterCorrector1;
import bloom.BloomFilterCorrector2;
import bloom.KCountArray7MTA;
import dna.AminoAcid;
import parse.Parse;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;
import ukmer.Kmer;

/** Standalone matched kernels for correction maintenance; not an end-to-end throughput test.
 * Workloads have literal expected results valid in both the pre-repair and repaired code.
 * Run each configuration in a fresh JVM. Timed kernels execute on the main thread;
 * table/input construction and reporting are outside the measurements.
 * Compare revisions within the same labeled configuration. Canonical and forward-only
 * query distributions differ and do not provide directly interchangeable speed figures.
 * @author Brian Bushnell, Sucrose */
public final class BloomCorrectorKernelBench{

	/** Selects one finite workload, warms it up, then measures main-thread CPU-bounded samples. */
	public static void main(final String[] args){
		String modeName=null, profile="supported";
		int k=31, seconds=10, warmup=3, samples=3;
		boolean canonical=true;
		for(String arg : args){
			final String[] pair=arg.split("=", 2);
			if(pair.length!=2){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String name=pair[0].toLowerCase(Locale.ROOT), value=pair[1];
			if(name.equals("mode")){modeName=value;}
			else if(name.equals("profile")){profile=value;}
			else if(name.equals("k")){k=Integer.parseInt(value);}
			else if(name.equals("seconds")){seconds=Integer.parseInt(value);}
			else if(name.equals("warmup")){warmup=Integer.parseInt(value);}
			else if(name.equals("samples")){samples=Integer.parseInt(value);}
			else if(name.equals("canonical")){canonical=Parse.parseBoolean(value);}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		final int mode=("mark".equals(modeName) ? MARK : "merge".equals(modeName) ? MERGE : "strand".equals(modeName) ? STRAND :
				"correction".equals(modeName) ? CORRECTION : "pass".equals(modeName) ? PASS :
				"passes".equals(modeName) ? QUERY_PASS : "counts".equals(modeName) ? BIG_COUNTS : -1);
		final boolean correctionWork=(mode==CORRECTION || mode==PASS);
		final boolean queryWork=(mode==QUERY_PASS || mode==BIG_COUNTS);
		check(mode>=0, "Use mode=mark|merge|strand|correction|pass|passes|counts.");
		check(seconds>=10 && warmup>=1 && samples>=1, "Samples require at least 10 main-thread CPU seconds, positive warmup and positive sample count.");
		check(mode==PASS ? k==31 || k==63 : mode==MERGE || mode==CORRECTION ? k==31 || k==62 || k==63 : k==31,
				"Pass accepts k31/63; merge/correction accept k31/62/63; other modes use k31.");
		check(mode==MARK ? profile.equals("supported") || profile.equals("low") || profile.equals("spaced") :
				mode==CORRECTION ? profile.equals("g3") || profile.equals("t2") : mode==PASS ? profile.equals("q30") :
				queryWork ? profile.equals("clean") || profile.equals("ambiguous") : profile.equals("supported"),
				"Mark profiles: supported|low|spaced; correction: g3|t2; pass: q30; passes/counts: clean|ambiguous; others: supported.");
		check(mode==STRAND || canonical, "Forward-only mode is a separate strand-query workload, not a merge/mark setting.");
		final ThreadMXBean clock=ManagementFactory.getThreadMXBean();
		check(clock.isCurrentThreadCpuTimeSupported(), "This benchmark requires an actual current-thread CPU clock.");
		if(!clock.isThreadCpuTimeEnabled()){clock.setThreadCpuTimeEnabled(true);}
		final int oldThreads=Shared.threads();
		final long oldCells=BloomFilter.OVERRIDE_CELLS;
		final boolean oldMask=Kmer.MASK_CORE;
		try{
			Shared.setThreads(1);
			BloomFilter.OVERRIDE_CELLS=CELLS;
			Kmer.MASK_CORE=false;
			final BloomCorrectorKernelBench bench=new BloomCorrectorKernelBench(mode, profile, k, canonical);
			check(bench.runBatch()==bench.expectedBatch, "Initial batch does not match the independently specified workload.");
			if(correctionWork){bench.validateCorrectionOutput();}
			if(queryWork){bench.validateQueryWork();}
			System.out.println("CONFIG mode="+modeName+" profile="+profile+" k="+k+
					" canonical="+(mode==MARK ? "na" : Boolean.toString(canonical))+" dual=false threads=1"+
					" pool="+(mode==STRAND ? bench.keys.length : correctionWork || queryWork ? 1 : READ_POOL)+
					" readLength="+(mode==MARK ? MARK_LENGTH : mode==MERGE ? MERGE_LENGTH : correctionWork ? CORRECTION_LENGTH : queryWork ? QUERY_LENGTH : 0)+
					" queryDistribution="+(mode!=STRAND ? "na" : canonical ? "both_orientations" : "already_canonical_subset")+
					" queryOrderSeed="+(mode==STRAND ? SEED+1 : 0)+
					" cells="+(mode==MARK ? 0 : CELLS)+" hashes="+(queryWork ? 1 : 3)+" seed="+SEED+" bits="+(queryWork ? 8 : 4)+
					(queryWork ? " ksmall="+(mode==BIG_COUNTS ? COUNT_SMALL_K : k)+" kbig=31 minConsecutive=3 threshold=1 loadedDepth=3" : "")+
					" unit="+(mode==STRAND ? "lookup" : "read")+" resetIncluded="+(mode==MARK || correctionWork || mode==BIG_COUNTS)+
					" correctionDepth="+bench.correctionDepth+" qualities="+(mode==PASS ? "uniform30_final28" : "na")+
					" resetSlots="+(mode==PASS ? "last_base_quality_count_and_tracker" : mode==CORRECTION ? "last_base_and_count" : mode==BIG_COUNTS ? "count_list_replaced_inside_kernel" : "na")+
					" batchUnits="+BATCH+" expectedBatchChecksum="+bench.expectedBatch+
					" warmupCpuSeconds="+warmup+" sampleCpuSeconds="+seconds+" samples="+samples+
					" assertions="+BloomCorrectorKernelBench.class.desiredAssertionStatus());
			bench.measure(clock, warmup, "warmup", 0);
			for(int sample=1; sample<=samples; sample++){bench.measure(clock, seconds, "sample", sample);}
		}finally{
			Shared.setThreads(oldThreads);
			BloomFilter.OVERRIDE_CELLS=oldCells;
			Kmer.MASK_CORE=oldMask;
		}
	}

	/** Allocates all state once; table population and oracle checks are outside measurement. */
	private BloomCorrectorKernelBench(final int mode_, final String profile_, final int k_, final boolean canonical){
		mode=mode_; profile=profile_; k=k_;
		final boolean correctionWork=(mode==CORRECTION || mode==PASS);
		final boolean queryWork=(mode==QUERY_PASS || mode==BIG_COUNTS);
		KCountArray7MTA.setSeed(SEED);
		filter=(mode==MARK ? null : queryWork ? new BloomFilter(mode==BIG_COUNTS ? COUNT_SMALL_K : k, k, 8, 1, 3, true, false, false, 0.01f) :
				new BloomFilter(null, null, null, k, k, 4, 3, 1, canonical, false, false, 0.01f, false));
		corrector=(queryWork ? null : (mode==MERGE || correctionWork) && k>31 ? new ReadyLongCorrector(filter, k) : new BloomFilterCorrector1(filter, k, k));
		originals=(mode==MARK ? new byte[READ_POOL][MARK_LENGTH] : null);
		bases=(mode==MARK ? new byte[MARK_LENGTH] : correctionWork ? new byte[CORRECTION_LENGTH] : queryWork ? new byte[QUERY_LENGTH] : null);
		qualities=(mode==MARK ? new byte[MARK_LENGTH] : mode==PASS ? new byte[CORRECTION_LENGTH] : null);
		counts=(mode==MARK ? new IntList(MARK_LENGTH-k+1) : correctionWork ? new IntList(CORRECTION_LENGTH-k+1) : mode==BIG_COUNTS ? new IntList(QUERY_LENGTH) : null);
		covered=(mode==MARK ? new BitSet(MARK_LENGTH) : mode==PASS ? new BitSet(CORRECTION_LENGTH) : null);
		reads=(mode==MERGE ? new Read[READ_POOL] : mode==QUERY_PASS ? new Read[1] : null);
		kmers=(mode==MERGE ? new LongList(MERGE_LENGTH) : mode==QUERY_PASS ? new LongList(QUERY_LENGTH) : null);
		keys=(mode==STRAND ? new long[canonical ? 2*QUERY_POOL : QUERY_POOL] : null);
		correctionBuilder=(correctionWork ? new ByteBuilder(bases) : null);
		rightCounts=(correctionWork ? new int[4] : null);
		correctionDepth=(mode==PASS ? 3 : mode==CORRECTION ? profile.equals("g3") ? 3 : 2 : 0);
		correctionBase=(byte)(mode==PASS || profile.equals("g3") ? 'G' : 'T');
		proposalLeft=(mode==PASS ? new ByteBuilder(CORRECTION_LENGTH) : null);
		proposalRight=(mode==PASS ? new ByteBuilder(CORRECTION_LENGTH) : null);
		proposalCounts=(mode==PASS ? new IntList(CORRECTION_LENGTH-k+1) : null);
		tracker=(mode==PASS ? new ErrorTracker() : null);
		queryInput=(queryWork ? new byte[QUERY_LENGTH] : null);
		queryExpectedCounts=(mode==BIG_COUNTS ? new int[QUERY_LENGTH-k+1] : null);
		queryCountBacking=(mode==BIG_COUNTS ? counts.array : null);
		expectedBatch=(mode==MARK ? prepareMarking() : mode==MERGE ? prepareMerges() :
				correctionWork ? prepareCorrection() : queryWork ? prepareQueryWork() : prepareQueries(canonical));
	}

	/** Prepares four equal-answer workloads for the direct streak and full-window count paths. */
	private long prepareQueryWork(){
		final int smallK=(mode==BIG_COUNTS ? COUNT_SMALL_K : k);
		final boolean ambiguous=profile.equals("ambiguous");
		Arrays.fill(bases, (byte)'C');
		if(mode==QUERY_PASS){
			Arrays.fill(bases, (byte)'A');
			Arrays.fill(bases, 0, k, (byte)'C');
			if(ambiguous){bases[k]='N'; Arrays.fill(bases, k+1, 2*k+1, (byte)'C');}
			else{Arrays.fill(bases, QUERY_LENGTH-k, QUERY_LENGTH, (byte)'C');}
			reads[0]=new Read(bases, null, "query-kernel", 0, false);
		}else if(ambiguous){bases[0]=bases[QUERY_LENGTH-1]='N';}
		System.arraycopy(bases, 0, queryInput, 0, bases.length);
		long canonicalC=0;
		for(int i=0; i<smallK; i++){canonicalC=(canonicalC<<2)|2;}//G...G is larger than C...C.
		filter.filter.increment(canonicalC, 3);
		check(filter.getCount(canonicalC)==3, "The one loaded canonical C/G key must have literal count three.");
		int hits=0;
		//Encode each window independently, outside clocks; reject collisions in every observed key.
		for(int start=0; start<=bases.length-smallK; start++){
			long forward=0, reverse=0;
			boolean defined=true, allC=true;
			for(int j=0; j<smallK; j++){
				final byte base=bases[start+j];
				final int code=(base=='A' ? 0 : base=='C' ? 1 : base=='G' ? 2 : base=='T' ? 3 : -1);
				if(code<0){defined=false; break;}
				allC&=base=='C';
				forward=(forward<<2)|code;
				reverse|=(long)(3-code)<<(2*j);
			}
			if(defined){
				check(filter.getCount(forward, reverse)==(allC ? 3 : 0), "Unexpected count/collision in an independently encoded query window.");
				if(allC){hits++;}
			}
		}
		if(mode==QUERY_PASS){
			check(hits==2, "Both pass workloads require exactly two isolated supported windows, below the required streak of three.");
			validateQueryWork();
			return BATCH;
		}
		for(int start=0; start<queryExpectedCounts.length; start++){
			boolean defined=true;
			for(int j=0; j<k; j++){if(bases[start+j]=='N'){defined=false;}}
			queryExpectedCounts[start]=(defined ? 3 : 0);
			if(defined){queryExpectedValid++;}
		}
		check(queryExpectedValid==(ambiguous ? 118 : 120), "Literal full-window geometry must match the independent base scan.");
		lastValid=filter.fillCountsBig(bases, counts);
		validateQueryWork();
		return (long)BATCH*(queryExpectedValid+queryExpectedCounts.length+queryExpectedCounts[0]+queryExpectedCounts[queryExpectedCounts.length-1]);
	}

	/** Checks complete inputs/results outside measurement; the buffered pass is an untimed control. */
	private void validateQueryWork(){
		check(Arrays.equals(bases, queryInput), "Query kernels must preserve every input base.");
		if(mode==QUERY_PASS){
			check(reads[0].bases==bases && reads[0].quality==null && reads[0].mate==null,
					"Direct query must retain the original unpaired read storage.");
			check(filter.passes(reads[0], 1) && filter.passes(reads[0], kmers, 1),
					"Both literal workloads must pass in direct and buffered query paths on either compared revision.");
		}else{
			check(lastValid==queryExpectedValid && counts.size==queryExpectedCounts.length && counts.array==queryCountBacking,
					"Full-window query must return the exact valid count/size and reuse its preallocated storage.");
			for(int i=0; i<queryExpectedCounts.length; i++){
				check(counts.array[i]==queryExpectedCounts[i], "Full-window count differs from its independent per-position oracle at "+i);
			}
		}
	}

	/** Validates equal-answer selector or full-pass workloads against a real C-only table. */
	private long prepareCorrection(){
		check(BloomFilterCorrector2.FORCE_MULT_SPLIT==null, "Correction timing requires natural long-window geometry.");
		Arrays.fill(bases, (byte)'C');
		if(mode==PASS){Arrays.fill(qualities, (byte)30);}
		final Kmer repeated=new Kmer(k);
		for(int i=0; i<k; i++){repeated.addRight((byte)'C');}
		final long key=(k<=31 ? filter.toKey(repeated.array1()[0], repeated.array2()[0]) : repeated.xor());
		filter.filter.increment(key, correctionDepth);
		final int windows=CORRECTION_LENGTH-k+1;
		final LongList stored=new LongList(windows);
		for(byte candidate : DNA){
			bases[CORRECTION_LENGTH-1]=candidate;
			check(corrector.fillKmers(bases, stored)==windows && stored.size()==windows,
					"Correction input must retain every valid window.");
			for(int i=0; i<windows; i++){
				final int expected=(i<windows-1 || candidate=='C' ? correctionDepth : 0);
				check(corrector.getCount2(stored.get(i))==expected, "A collision or encoding mismatch invalidated a correction candidate.");
			}
			if(mode==PASS){
				final byte[] reverse=bases.clone();
				AminoAcid.reverseComplementBasesInPlace(reverse);
				check(corrector.fillKmers(reverse, stored)==windows && stored.size()==windows, "Reverse proposals require complete dense windows.");
				for(int i=0; i<windows; i++){
					check(corrector.getCount2(stored.get(i))==(i>0 || candidate=='C' ? correctionDepth : 0),
							"A collision or encoding mismatch invalidated the reverse proposal table.");
				}
			}
		}
		for(int i=0; i<windows; i++){counts.add(correctionDepth);}
		resetCorrection();
		validateCorrectionInput();
		lastCorrections=(mode==PASS ? runPass() : corrector.reassemble_inner(correctionBuilder, null, rightCounts, counts, 0));
		validateCorrectionOutput();
		//Full input/output checks establish exactly which input slots need resetting.
		resetCorrection();
		validateCorrectionInput();
		return (long)BATCH*(1+'C'+correctionDepth+(mode==PASS ? 28 : 0));
	}

	/** Restores the altered input base/window and, for a full pass, its quality and tracker. */
	private void resetCorrection(){
		bases[CORRECTION_LENGTH-1]=correctionBase;
		counts.array[counts.size-1]=0;
		if(mode==PASS){qualities[CORRECTION_LENGTH-1]=30; tracker.clear();}
	}

	/** Checks every input base/count and exact backing capacity outside timing. */
	private void validateCorrectionInput(){
		check(correctionBuilder.array==bases && bases.length==CORRECTION_LENGTH && correctionBuilder.length()==CORRECTION_LENGTH,
				"Correction requires unchanged exact-capacity backing state.");
		check(counts.size==CORRECTION_LENGTH-k+1, "Correction count list must have one entry per window.");
		for(int i=0; i<bases.length; i++){
			check(bases[i]==(i==bases.length-1 ? correctionBase : 'C'), "Correction input base differs from its literal oracle.");
		}
		for(int i=0; i<counts.size; i++){
			check(counts.get(i)==(i==counts.size-1 ? 0 : correctionDepth), "Correction input count differs from its literal oracle.");
		}
		if(mode==PASS){
			for(byte quality : qualities){check(quality==30, "Pass input/reset must have uniform quality thirty.");}
			check(tracker.corrected()==0 && tracker.detected()==0 && !tracker.rollback, "Pass tracker reset must clear prior corrections.");
			check(proposalLeft.array.length==bases.length && proposalRight.array.length==bases.length,
					"Pass proposals must retain exact read-sized backing arrays.");
		}
	}

	/** Checks full output after preflight and each measured phase, outside its clocks. */
	private void validateCorrectionOutput(){
		check(lastCorrections==1, "Both compared versions must perform exactly one correction in this workload.");
		check(correctionBuilder.array==bases && bases.length==CORRECTION_LENGTH && correctionBuilder.length()==CORRECTION_LENGTH,
				"Correction must preserve exact backing state and logical length.");
		check(counts.size==CORRECTION_LENGTH-k+1, "Correction must preserve the dense count-list size.");
		for(byte base : bases){check(base=='C', "Corrected sequence must be entirely C.");}
		for(int i=0; i<counts.size; i++){check(counts.get(i)==correctionDepth, "Every corrected window must have the specified depth.");}
		if(mode==PASS){
			check(qualities.length==bases.length && proposalLeft.length()==bases.length && proposalRight.length()==bases.length,
					"Pass must preserve all sequence/quality lengths.");
			for(int i=0; i<bases.length; i++){
				check(qualities[i]==(i==bases.length-1 ? 28 : 30), "Only the terminal quality changes, to the one-sided cap28.");
				check(proposalLeft.get(i)==(i==bases.length-1 ? 'C' : 0) && proposalRight.get(i)==0,
						"Only the left scan supplies the terminal correction proposal.");
			}
			check(covered.cardinality()==1 && covered.get(bases.length-1), "Pass must mark only the corrected terminal base.");
			check(tracker.corrected()==1 && tracker.correctedReassembleOuter==1 && tracker.correctedReassembleInner==0 &&
					tracker.detectedReassemble==1 && !tracker.rollback, "Pass tracker must describe one retained one-sided correction.");
		}
	}

	/** Executes the complete directional pass with non-null qualities and reused caller scratch. */
	private int runPass(){
		return corrector.reassemble_pass(bases, qualities, proposalLeft, proposalRight, rightCounts, counts,
				proposalCounts, tracker, 0, null, null, covered);
	}

	/** Specifies support as all bases, no bases, or the literal intervals [0,31) and [62,93). */
	private long prepareMarking(){
		for(int i=0; i<MARK_LENGTH-k+1; i++){
			counts.add(profile.equals("supported") || (profile.equals("spaced") && (i==0 || i==62)) ? 5 : 0);
		}
		int marked=0;
		for(int i=0; i<MARK_LENGTH; i++){if(!supportedBase(i)){marked++;}}
		long poolChecksum=0;
		for(int r=0; r<READ_POOL; r++){
			for(int i=0; i<MARK_LENGTH; i++){originals[r][i]=DNA[(i+r)&3];}
			poolChecksum+=marked+(supportedBase(0) ? originals[r][0]+30 : 'N')+
					(supportedBase(MARK_LENGTH-1) ? originals[r][MARK_LENGTH-1]+30 : 'N');
			System.arraycopy(originals[r], 0, bases, 0, MARK_LENGTH);
			Arrays.fill(qualities, (byte)30);
			check(corrector.markBadBases(bases, qualities, counts, covered, 2, false, (byte)0)==marked,
					"Marking count must match the independent support intervals before timing.");
			for(int i=0; i<MARK_LENGTH; i++){
				final boolean supported=supportedBase(i);
				if(bases[i]!=(supported ? originals[r][i] : 'N') || qualities[i]!=(supported ? 30 : 0) || covered.get(i)!=supported){
					throw new AssertionError("Marking oracle failed before timing: profile="+profile+", read="+r+", base="+i);
				}
			}
		}
		return poolChecksum*(BATCH/READ_POOL);
	}

	/** Independent base-centric support specification, used only during setup. */
	private boolean supportedBase(final int pos){
		return profile.equals("supported") || (profile.equals("spaced") && (pos<31 || (pos>=62 && pos<93)));
	}

	/** Gives every distinct window in 32 deterministic reads count five; all tested ranges fit both versions. */
	private long prepareMerges(){
		final int windows=MERGE_LENGTH-k+1;
		check(MATE_LENGTH+WIDTH<windows && MERGE_LENGTH-MATE_LENGTH-k-WIDTH>=0,
				"All four merge windows must be in bounds in both baseline and repaired code.");
		final long[] loaded=new long[READ_POOL*windows];
		final Random random=new Random(SEED);
		final Kmer rolling=new Kmer(k);
		int size=0;
		for(int r=0; r<READ_POOL; r++){
			final byte[] sequence=new byte[MERGE_LENGTH];
			for(int i=0; i<sequence.length; i++){sequence[i]=DNA[random.nextInt(4)];}
			reads[r]=new Read(sequence, null, "kernel"+r, r, false);
			rolling.clear();
			for(byte base : sequence){
				rolling.addRight(base);
				if(rolling.len()>=k){
					loaded[size++]=(k<=31 ? filter.toKey(rolling.array1()[0], rolling.array2()[0]) : rolling.xor());
				}
			}
		}
		check(size==loaded.length, "All deterministic input windows must be valid DNA.");
		Arrays.sort(loaded);
		for(int i=0; i<loaded.length; i++){
			if(i==0 || loaded[i]!=loaded[i-1]){filter.filter.increment(loaded[i], 5);}
		}
		for(Read read : reads){
			check(corrector.fillKmers(read.bases, kmers)==windows && kmers.size()==windows,
					"Merge preparation must retain dense complete window positions.");
			for(int i=0; i<windows; i++){
				check(corrector.getCount2(kmers.get(i))==5, "A collision or encoding mismatch invalidated the uniform-count merge workload.");
			}
		}
		return BATCH;
	}

	/** Loads 64 distinct T...T keys; canonical mode queries both orientations, forward-only mode the high subset. */
	private long prepareQueries(final boolean canonical){
		final long interiorMask=(1L<<58)-1;
		for(int i=0; i<QUERY_POOL; i++){
			final long key=(3L<<60)|(((i*0x9E3779B97F4A7C15L)&interiorMask)<<2)|3;
			if(canonical){
				long reverse=0, remaining=key;
				for(int base=0; base<31; base++){reverse=(reverse<<2)|(3L^(remaining&3)); remaining>>>=2;}
				check(reverse<key, "T...T must be the larger strand; its reverse complement begins with A.");
				keys[2*i]=key;
				keys[2*i+1]=reverse;
			}else{keys[i]=key;}
			if((i&1)==0){filter.filter.increment(key, 7);}
		}
		for(int i=0; i<keys.length; i++){
			final int identity=(canonical ? i>>>1 : i), expected=((identity&1)==0 ? 7 : 0);
			final long stored=(canonical ? keys[i&~1] : keys[i]);
			check(filter.filter.read(stored)==expected && corrector.getCount2(keys[i])==expected,
					"A collision or strand mismatch invalidated the independently specified query pool.");
		}
		//Mix direction and hit/miss order deterministically, rather than alternate them.
		final Random order=new Random(SEED+1);
		for(int i=keys.length-1; i>0; i--){
			final int j=order.nextInt(i+1);
			final long saved=keys[i]; keys[i]=keys[j]; keys[j]=saved;
		}
		return (long)BATCH*7/2;
	}

	/** Runs only the selected kernel; all arrays, lists, reads and tables are reused. */
	private long runBatch(){
		long checksum=0;
		if(mode==MARK){
			for(int i=0; i<BATCH; i++){
				System.arraycopy(originals[i&(READ_POOL-1)], 0, bases, 0, MARK_LENGTH);
				Arrays.fill(qualities, (byte)30);
				checksum+=corrector.markBadBases(bases, qualities, counts, covered, 2, false, (byte)0);
				checksum+=bases[0]+bases[MARK_LENGTH-1]+qualities[0]+qualities[MARK_LENGTH-1];
			}
		}else if(mode==MERGE){
			for(int i=0; i<BATCH; i++){
				checksum+=corrector.mergeOK(reads[i&(READ_POOL-1)], MATE_LENGTH, MATE_LENGTH, kmers, WIDTH, 3, 80) ? 1 : 0;
			}
		}else if(mode==PASS){
			for(int i=0; i<BATCH; i++){
				resetCorrection();
				lastCorrections=runPass();
				if(lastCorrections!=1){throw new AssertionError("Timed pass must change exactly one base.");}
				checksum+=lastCorrections+bases[CORRECTION_LENGTH-1]+qualities[CORRECTION_LENGTH-1]+counts.array[counts.size-1];
			}
		}else if(mode==CORRECTION){
			for(int i=0; i<BATCH; i++){
				resetCorrection();
				lastCorrections=corrector.reassemble_inner(correctionBuilder, null, rightCounts, counts, 0);
				if(lastCorrections!=1){throw new AssertionError("Timed correction must change exactly one base.");}
				checksum+=lastCorrections+bases[CORRECTION_LENGTH-1]+counts.array[counts.size-1];
			}
		}else if(mode==QUERY_PASS){
			for(int i=0; i<BATCH; i++){checksum+=filter.passes(reads[0], 1) ? 1 : 0;}
		}else if(mode==BIG_COUNTS){
			for(int i=0; i<BATCH; i++){
				lastValid=filter.fillCountsBig(bases, counts);
				checksum+=lastValid+counts.size+counts.array[0]+counts.array[counts.size-1];
			}
		}else{
			for(int i=0; i<BATCH; i++){checksum+=corrector.getCount2(keys[i&(keys.length-1)]);}
		}
		return checksum;
	}

	/** Includes resets/checksum consumption; excludes setup, output and all subsequent phases. */
	private void measure(final ThreadMXBean clock, final int seconds, final String phase, final int sample){
		if(mode==QUERY_PASS || mode==BIG_COUNTS){validateQueryWork();}
		final long cpuStart=clock.getCurrentThreadCpuTime();
		check(cpuStart>=0, "The main-thread CPU clock must be enabled and available.");
		final long wallStart=System.nanoTime(), target=seconds*1000000000L;
		long cpuEnd, work=0, checksum=0;
		do{
			final long batchChecksum=runBatch();
			if(batchChecksum!=expectedBatch){
				throw new AssertionError("Kernel output changed during measurement: expected "+expectedBatch+", observed "+batchChecksum);
			}
			checksum+=batchChecksum;
			work+=BATCH;
			sink=checksum;
			cpuEnd=clock.getCurrentThreadCpuTime();
			check(cpuEnd>=cpuStart, "The current-thread CPU clock stopped or moved backward.");
		}while(cpuEnd-cpuStart<target);
		final long wall=System.nanoTime()-wallStart, cpu=cpuEnd-cpuStart;
		if(mode==CORRECTION || mode==PASS){validateCorrectionOutput();}
		if(mode==QUERY_PASS || mode==BIG_COUNTS){validateQueryWork();}
		System.out.printf(Locale.ROOT, "RESULT phase=%s sample=%d work=%d checksum=%d wall_seconds=%.6f thread_cpu_seconds=%.6f wall_ns_per_unit=%.3f cpu_ns_per_unit=%.3f%n",
				phase, sample, work, checksum, wall*1e-9, cpu*1e-9, wall/(double)work, cpu/(double)work);
	}

	/** Keeps setup and checksum failures active regardless of the JVM assertion mode. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Initializes normal scratch once without overriding any encoding or lookup method. */
	private static final class ReadyLongCorrector extends BloomFilterCorrector2{
		ReadyLongCorrector(final BloomFilter filter, final int k){super(filter, k, k); initializeThreadLocals();}
	}

	private final int mode, k;
	private final String profile;
	private final BloomFilter filter;
	private final BloomFilterCorrector corrector;
	private final byte[][] originals;
	private final byte[] bases, qualities;
	private final IntList counts;
	private final BitSet covered;
	private final Read[] reads;
	private final LongList kmers;
	private final long[] keys;
	private final long expectedBatch;
	private final ByteBuilder correctionBuilder;
	private final int[] rightCounts;
	private final int correctionDepth;
	private final byte correctionBase;
	private final ByteBuilder proposalLeft, proposalRight;
	private final IntList proposalCounts;
	private final ErrorTracker tracker;
	private int lastCorrections;
	/** Independent query input/count oracles and retained backing array; unused by older modes. */
	private final byte[] queryInput;
	private final int[] queryExpectedCounts, queryCountBacking;
	private int queryExpectedValid, lastValid;
	private static volatile long sink;
	private static final int MARK=0, MERGE=1, STRAND=2, CORRECTION=3, PASS=4, QUERY_PASS=5, BIG_COUNTS=6;
	private static final int READ_POOL=32, QUERY_POOL=64, BATCH=65536;
	private static final int MARK_LENGTH=150, MERGE_LENGTH=300, MATE_LENGTH=200, WIDTH=4;
	private static final int CORRECTION_LENGTH=150;
	private static final int QUERY_LENGTH=150, COUNT_SMALL_K=21;
	private static final long CELLS=2000003, SEED=20261105;
	private static final byte[] DNA={'A', 'C', 'G', 'T'};
}
