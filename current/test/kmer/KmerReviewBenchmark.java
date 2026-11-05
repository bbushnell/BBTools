package test.kmer;

import java.util.Locale;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.util.BitSet;

import assemble.Contig;
import kmer.AbstractKmerTable;
import kmer.AbstractKmerTableSet;
import kmer.HashArray1D;
import kmer.HashArray;
import kmer.HashArrayHybrid;
import kmer.HashArrayHybridFast;
import kmer.HashBuffer;
import kmer.HashForest;
import kmer.KmerNode;
import kmer.KmerTableSet;
import kmer.TableReader;
import kmer.Walker;
import stream.Read;
import structures.IntList;
import ukmer.AbstractKmerTableU;
import ukmer.HashArrayU;
import ukmer.HashArrayU1D;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/** Matched single-thread kernels for counts, canonical reversal and existing-value Node2D merges.
 * iterations=0 calibrates at least15seconds of main-thread CPU work; repeat the printed
 * fixed iteration count on both source revisions. This is not a whole-program benchmark.
 * node2d-merge measures steady four-value merges through HashForest, not insertion or resize.
 * hybrid-primary/victim and fast-primary/victim measure existing scalar count updates.
 * Only these four modes accept expect=baseline|fixed (default fixed); baseline victim
 * counts intentionally retain the known KMR021 error, with a different expected checksum.
 * long-fill-all repeatedly counts four unmasked long-key windows using fixed storage;
 * its checksum includes every returned valid count and all four output counts.
 * walker repeatedly drains complete public walkers over a stable fresh HashArray1D;
 * each walk allocates the Walker and its victim-node list.
 * reader-best-success measures the successful TableReader.findBestMatch path only.
 * reader-best-sampled-success uses the same work/oracles with speed1; all query
 * keys pass both the legacy and modulo17 selectors, preserving equal lookup work.
 * reader-markbits measures the unpadded successful TableReader.markBits path only,
 * preserving its per-call BitSet allocation.
 * coverage-present measures both all-present short-key calcCoverage overloads with
 * one reused byte array, one reused Contig and one fixed shard.
 * coverage-long-present does the same for the long-key ukmer overloads.
 * @author Shinobu, Sucrose, Jean */
public final class KmerReviewBenchmark {

	/** Selects one kernel and either a fixed work count or the initial calibration. */
	public static void main(String[] args){
		String mode="buffer";
		String expect=null;
		long iterations=0;
		for(String arg : args){
			final int equals=arg.indexOf('=');
			if(equals<1){throw new IllegalArgumentException("Expected mode= or iterations=: "+arg);}
			final String key=arg.substring(0, equals), value=arg.substring(equals+1);
			if(key.equals("mode")){mode=value;}
			else if(key.equals("iterations")){iterations=Long.parseLong(value);}
			else if(key.equals("expect")){expect=value;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(iterations<0 || iterations>1000000000000L){throw new IllegalArgumentException("Invalid work count: "+iterations);}
		if(expect!=null && (!isHybridMode(mode) || !(expect.equals("baseline") || expect.equals("fixed")))){
			throw new IllegalArgumentException("expect=baseline|fixed is supported only by the four hybrid count modes.");
		}
		if(!CPU_CLOCK.isCurrentThreadCpuTimeSupported()){throw new IllegalStateException("Current-thread CPU timing is unavailable.");}
		if(!CPU_CLOCK.isThreadCpuTimeEnabled()){CPU_CLOCK.setThreadCpuTimeEnabled(true);}
		final KmerReviewBenchmark bench=new KmerReviewBenchmark(mode, expect);
		bench.step(BLOCK);
		bench.flush();
		if(bench.mergeTable!=null){bench.verifyMerge();}
		if(bench.hybridTable!=null){bench.verifyHybrid(BLOCK);}
		if(bench.longFillSet!=null){bench.verifyLongFill(BLOCK);}
		long done=0;
		final long cpuStart=cpuTime();
		final long start=System.nanoTime();
		do{
			final int count=iterations==0 ? BLOCK : (int)Math.min(BLOCK, iterations-done);
			bench.step(count);
			done+=count;
		}while(iterations==0 ? cpuTime()-cpuStart<TARGET_NANOS : done<iterations);
		bench.flush();
		final long nanos=System.nanoTime()-start;
		bench.verify(done+BLOCK);
		blackhole=bench.checksum;
		if(bench.hybridTable!=null){
			System.err.println("hybrid_expect="+(bench.legacyHybrid ? "baseline" : "fixed")+
				", tier="+(bench.hybridVictim ? "victim" : "primary")+", created="+bench.hybridCreated);
		}
		System.out.printf(Locale.ROOT, "mode\titerations\tseconds\tchecksum\n%s\t%d\t%.6f\t%d%n",
				mode, done, nanos*1e-9, blackhole);
	}

	/** Keeps the known-wrong baseline oracle confined to the four explicitly named hybrid modes. */
	private static boolean isHybridMode(final String mode){
		return mode.equals("hybrid-primary") || mode.equals("hybrid-victim") ||
			mode.equals("fast-primary") || mode.equals("fast-victim");
	}

	/** Rejects unsupported/disabled timing rather than accepting a stalled wall-clock sample. */
	private static long cpuTime(){
		final long nanos=CPU_CLOCK.getCurrentThreadCpuTime();
		if(nanos<0){throw new IllegalStateException("Current-thread CPU timing did not return a measurement.");}
		return nanos;
	}

	/** Allocates a small stable live set; measured loops retain their arrays and keys. */
	private KmerReviewBenchmark(final String mode, final String expect){
		if(mode.equals("buffer") || mode.equals("sorted-buffer")){
			HashBuffer.SORT_BUFFERS=mode.equals("sorted-buffer");
			tables=new HashArray1D[31];
			for(int i=0; i<tables.length; i++){tables[i]=new HashArray1D(new int[] {5003}, -1L);}
			buffer=new HashBuffer(tables, 1024, 31, false, false);
			for(int i=0; i<keys.length; i++){keys[i]=(15485863L*i)&Long.MAX_VALUE;}
		}else if(mode.equals("victim")){
			tables=new HashArray1D[] {new HashArray1D(new int[] {100003}, Long.MAX_VALUE)};
			for(int i=0; i<keys.length; i++){keys[i]=100003L*i;}
		}else if(mode.equals("reverse-canonical")){
			rolling=new Kmer(95);
			for(int i=0; i<rolling.kbig; i++){rolling.addRightNumeric(i&3);}
		}else if(mode.equals("node2d-merge")){
			mergeTable=new HashForest(8191, false, true);
			mergeValues=new int[KEYS][];
			for(int i=0; i<KEYS; i++){
				keys[i]=i;
				final int[] values=new int[MERGE_VALUES+1];
				for(int j=0; j<MERGE_VALUES; j++){values[j]=i*MERGE_VALUES+j+1;}
				values[MERGE_VALUES]=-1;
				mergeValues[i]=values;
				if(mergeTable.set(keys[i], values, MERGE_VALUES)!=1){
					throw new AssertionError("Each initial forest key must be new: "+i);
				}
			}
		}else if(isHybridMode(mode)){
			hybridVictim=mode.endsWith("-victim");
			legacyHybrid="baseline".equals(expect);
			hybridTable=mode.startsWith("fast-") ? new HashArrayHybridFast(new int[] {100003}, -1L) :
				new HashArrayHybrid(new int[] {100003}, -1L);
			if(hybridVictim){
				for(int i=0; i<60; i++){
					if(hybridTable.incrementAndReturnNumCreated(100003L*i, 1)!=1){
						throw new AssertionError("Each primary filler must be new: "+i);
					}
				}
				hybridValueStorage=new int[KEYS][];
			}
			for(int i=0; i<KEYS; i++){
				keys[i]=hybridVictim ? 100003L*(60+i) : i;
				if(hybridTable.incrementAndReturnNumCreated(keys[i], 1)!=1){
					throw new AssertionError("Each timed hybrid key must be new: "+i);
				}
				if(hybridVictim){hybridValueStorage[i]=hybridTable.getValues(keys[i], null);}
			}
			hybridKeyStorage=hybridTable.array();
			hybridVictimStorage=hybridTable.victims().array();
			verifyHybrid(0);
		}else if(mode.equals("long-fill-all")){
			setupLongFill();
		}else if(mode.equals("walker")){
			setupWalker();
		}else if(mode.equals("reader-best-success")){
			setupReaderBest(0);
		}else if(mode.equals("reader-best-sampled-success")){
			setupReaderBest(1);
		}else if(mode.equals("reader-markbits")){
			setupReaderMarkBits();
		}else if(mode.equals("coverage-present")){
			setupCoveragePresent();
		}else if(mode.equals("coverage-long-present")){
			setupCoverageLongPresent();
		}else{throw new IllegalArgumentException("Unknown mode: "+mode);}
	}

	/** Installs the finite UKM110 unmasked fixture without invoking the main table allocator.
	 * This standalone mode fixes packing/orientation settings before constructing keys. */
	private void setupLongFill(){
		assert(longFillSet==null) : "The long-fill fixture must be constructed once, before warmup and timing.";
		Kmer.PACKED=false;
		Kmer.MASK_CORE=false;
		AbstractKmerTableSet.MASK_CORE=false;
		AbstractKmerTableSet.MASK_MIDDLE=false;
		longFillSet=new KmerTableSetU(new String[] {"in=/dev/null", "k=62", "ways=1", "threads=1",
			"prealloc=f", "prefilter=0", "onepass=f", "rcomp=f", "qin=33", "interleaved=f"}, 0);
		longFillTable=new HashArrayU1D(new int[] {101}, longFillSet.k, longFillSet.kbig);
		try{
			final Field field=KmerTableSetU.class.getDeclaredField("tables");
			field.setAccessible(true);
			field.set(longFillSet, new AbstractKmerTableU[] {longFillTable});
		}catch(ReflectiveOperationException e){
			throw new IllegalStateException("Unable to install the bounded long-fill fixture.", e);
		}
		longFillBases="ACGTTGCAACGAGTCCGATGATCGTACGTTAGCATCGATGCTAGCATGACCGTTAACGTTGCATT".
			getBytes(StandardCharsets.US_ASCII);
		if(longFillSet.kbig()!=62 || longFillSet.tables().length!=1 ||
				longFillTable.arrayLength()!=161 || longFillBases.length!=65){
			throw new AssertionError("Long-fill requires four k62 windows and one fixed161-slot table.");
		}
		longFillCounts=new IntList(4);
		longFillCountStorage=longFillCounts.array;
		longFillScratch=new Kmer(62);
		longFillForward=longFillScratch.array1();
		longFillReverse=longFillScratch.array2();
		for(int start=0; start<3; start++){
			longFillScratch.clear();
			for(int i=start; i<start+62; i++){longFillScratch.addRight(longFillBases[i]);}
			final int value=start==0 ? 7 : start==1 ? 11 : 13;
			if(longFillScratch.len()!=62 || longFillTable.set(longFillScratch, value)!=1){
				throw new AssertionError("Each long-fill fixture window must install a distinct complete key: "+start);
			}
		}
	}

	/** Runs a block divisible by the key period except for the final fixed-work tail. */
	private void step(final int count){
		assert(count>0 && count<=BLOCK) : "Measured blocks must advance the fixed work count within the configured bound.";
		if(buffer!=null){
			for(int i=0; i<count; i++){buffer.incrementAndReturnNumCreated(keys[i&(KEYS-1)], 1);}
		}else if(mergeTable!=null){
			int created=mergeCreated;
			for(int i=0; i<count; i++){
				final int index=i&(KEYS-1);
				created|=mergeTable.set(keys[index], mergeValues[index], MERGE_VALUES);
			}
			mergeCreated=created;
		}else if(hybridTable!=null){
			int created=hybridCreated;
			for(int i=0; i<count; i++){
				created|=hybridTable.incrementAndReturnNumCreated(keys[i&(KEYS-1)], 1);
			}
			hybridCreated=created;
		}else if(longFillSet!=null){
			long sum=checksum;
			int valid=0;
			for(int i=0; i<count; i++){
				valid=longFillSet.fillSpecificCounts(longFillBases, longFillCounts, null, longFillScratch);
				sum+=valid+longFillCounts.get(0)+2L*longFillCounts.get(1)+
					3L*longFillCounts.get(2)+4L*longFillCounts.get(3);
			}
			longFillValid=valid;
			checksum=sum;
		}else if(walkerTable!=null){
			stepWalker(count);
		}else if(readerBestReader!=null){
			stepReaderBest(count);
		}else if(readerMarkBitsReader!=null){
			stepReaderMarkBits(count);
		}else if(coverageSet!=null){
			stepCoveragePresent(count);
		}else if(coverageLongSet!=null){
			stepCoverageLongPresent(count);
		}else if(rolling!=null){
			long sum=checksum;
			for(int i=0; i<count; i++){
				rolling.addRightNumeric(i&3);
				final long hash=rolling.xor();
				rolling.rcomp();
				final long reverseHash=rolling.xor();
				if(reverseHash!=hash){throw new AssertionError("Canonical reversal changed the placement hash.");}
				rolling.rcomp();
				sum=(sum^hash)+1;
			}
			checksum=sum;
		}else{
			final HashArray1D table=tables[0];
			for(int i=0; i<count; i++){table.increment(keys[i&(KEYS-1)], 2);}
		}
	}

	/** Includes pending table updates in both the warmup boundary and timed work. */
	private void flush(){if(buffer!=null){buffer.flush();}}

	/** Checks every literal count and both physical storage tiers after the timed loop. */
	private void verify(final long total){
		if(longFillSet!=null){verifyLongFill(total); return;}
		if(walkerTable!=null){verifyWalker(total); return;}
		if(readerBestReader!=null){verifyReaderBest(total); return;}
		if(readerMarkBitsReader!=null){verifyReaderMarkBits(total); return;}
		if(coverageSet!=null){verifyCoveragePresent(total); return;}
		if(coverageLongSet!=null){verifyCoverageLongPresent(total); return;}
		if(hybridTable!=null){verifyHybrid(total); return;}
		if(mergeTable!=null){verifyMerge(); return;}
		if(rolling!=null){
			if(!rolling.verify(true)){throw new AssertionError("Rolling strands diverged.");}
			return;
		}
		long size=0;
		for(HashArray1D table : tables){size+=table.size()+table.victims().size();}
		if(size!=KEYS){throw new AssertionError("Lost or created keys: "+size);}
		if(buffer==null && tables[0].victims().size()!=KEYS-60){
			throw new AssertionError("Victim kernel did not exercise its intended storage tier.");
		}
		long sum=0;
		for(int i=0; i<KEYS; i++){
			final long expected=(total/KEYS+(i<total%KEYS ? 1 : 0))*(buffer==null ? 2 : 1);
			final int actual=buffer==null ? tables[0].getValue(keys[i]) : buffer.getValue(keys[i]);
			if(actual!=expected){throw new AssertionError("Count mismatch at key "+i+": "+actual+" versus "+expected);}
			sum+=(i+1L)*actual;
		}
		checksum=sum;
	}

	/** Checks the full vector, cumulative checksum and reused storage outside timing.
	 * total includes warmup so skipped or extra calls cannot silently pass the checksum. */
	private void verifyLongFill(final long total){
		assert(total>0) : "Long-fill verification follows at least the completed warmup block.";
		if(longFillValid!=4 || longFillCounts.size()!=4 || longFillCounts.get(0)!=7 ||
				longFillCounts.get(1)!=11 || longFillCounts.get(2)!=13 || longFillCounts.get(3)!=0 ||
				checksum!=72L*total){
			throw new AssertionError("Long-fill requires valid4, counts[7,11,13,0] and checksum72 per call: "+
				"valid="+longFillValid+", counts="+longFillCounts+", checksum="+checksum+", calls="+total);
		}
		if(longFillCounts.array!=longFillCountStorage || longFillScratch.array1()!=longFillForward ||
				longFillScratch.array2()!=longFillReverse || longFillTable.arrayLength()!=161 ||
				longFillSet.tables().length!=1 || longFillSet.tables()[0]!=longFillTable ||
				longFillTable.size()+longFillTable.victims().size()!=3){
			throw new AssertionError("Long-fill changed its preallocated storage or key count.");
		}
	}

	/** Installs four literal positive counts in a fresh table with zero-filled empty cells. */
	private void setupWalker(){
		assert(walkerTable==null) : "The walker fixture must be constructed once, before warmup and timing.";
		walkerTable=new HashArray1D(new int[] {5}, -1L);
		for(int i=0; i<WALKER_RECORDS; i++){
			final int actual=walkerTable.increment(WALKER_KEYS[i], WALKER_COUNTS[i]);
			if(actual!=WALKER_COUNTS[i]){
				throw new AssertionError("Walker literal count was not installed at key "+WALKER_KEYS[i]+
					": "+actual+" versus "+WALKER_COUNTS[i]);
			}
		}
		walkerKeyStorage=walkerTable.array();
		walkerValueStorage=walkerTable.values();
		verifyWalkerStorage();
	}

	/** Times complete primary-then-victim walks through the public cursor API. */
	private void stepWalker(final int count){
		long sum=checksum;
		for(int i=0; i<count; i++){
			long walkSum=0;
			int records=0;
			for(final Walker walker=walkerTable.walk(); walker.next(); ){
				walkSum+=walker.kmer()*WALKER_KEY_WEIGHT+walker.value()*WALKER_VALUE_WEIGHT;
				records++;
			}
			if(records!=WALKER_RECORDS || walkSum!=WALKER_CHECKSUM){
				throw new AssertionError("Walker emitted "+records+" records and checksum "+walkSum+
					" instead of "+WALKER_RECORDS+" records and checksum "+WALKER_CHECKSUM);
			}
			sum+=walkSum;
		}
		checksum=sum;
	}

	/** Checks the literal fresh table and accumulated full-walk checksum outside timing. */
	private void verifyWalker(final long total){
		assert(total>0) : "Walker verification follows at least the completed warmup block.";
		verifyWalkerStorage();
		if(checksum!=WALKER_CHECKSUM*total){
			throw new AssertionError("Walker checksum "+checksum+" does not equal "+
				WALKER_CHECKSUM+" per complete walk for "+total+" walks.");
		}
	}

	/** Keeps walker storage fixed and verifies every literal count after timed walking. */
	private void verifyWalkerStorage(){
		if(walkerTable.size()!=WALKER_RECORDS || walkerTable.victims().size()!=0 ||
				walkerTable.arrayLength()!=65 || walkerTable.array()!=walkerKeyStorage ||
				walkerTable.values()!=walkerValueStorage){
			throw new AssertionError("Walker table changed storage or key placement.");
		}
		for(int i=0; i<WALKER_RECORDS; i++){
			final int actual=walkerTable.getValue(WALKER_KEYS[i]);
			if(actual!=WALKER_COUNTS[i]){
				throw new AssertionError("Walker count mismatch at key "+WALKER_KEYS[i]+
					": "+actual+" versus "+WALKER_COUNTS[i]);
			}
		}
		int seen=0;
		for(final Walker walker=walkerTable.walk(); walker.next(); ){
			final long kmer=walker.kmer();
			final int value=walker.value();
			if(kmer<1 || kmer>WALKER_RECORDS){
				throw new AssertionError("Walker emitted unexpected key "+kmer+" with value "+value);
			}
			final int bit=1<<((int)kmer-1);
			if((seen&bit)!=0){throw new AssertionError("Walker emitted duplicate key "+kmer);}
			seen|=bit;
			if(value!=WALKER_COUNTS[(int)kmer-1]){
				throw new AssertionError("Walker emitted value "+value+" for key "+kmer+
					" instead of "+WALKER_COUNTS[(int)kmer-1]);
			}
		}
		if(seen!=(1<<WALKER_RECORDS)-1){throw new AssertionError("Walker emitted key mask "+seen);}
	}

	/** Installs the same AAA/ACG fixture for successful unsampled or speed1 lookup.
	 * @param speed Zero preserves the original mode; one exercises key filtering */
	private void setupReaderBest(final int speed){
		assert(readerBestReader==null) : "The reader-best fixture must be constructed once.";
		readerBestReader=new TableReader(3, 0, speed, 0, 0, false, false);
		readerBestReader.minHits=2;
		readerBestReader.countArray=new int[3];
		readerBestTables=new AbstractKmerTable[7];
		readerBestArrays=new HashArray1D[readerBestTables.length];
		readerBestKeyStorage=new long[readerBestTables.length][];
		readerBestValueStorage=new int[readerBestTables.length][];
		for(int i=0; i<readerBestTables.length; i++){
			final HashArray1D table=new HashArray1D(new int[] {17}, -1L);
			readerBestTables[i]=readerBestArrays[i]=table;
			readerBestKeyStorage[i]=table.array();
			readerBestValueStorage[i]=table.values();
		}
		if(readerBestTables[1].set(64, 1)!=1 || readerBestTables[0].set(70, 2)!=1){
			throw new AssertionError("Reader-best literal keys AAA and ACG must install as new IDs.");
		}
		readerBestRead=new Read("AAAACG".getBytes(StandardCharsets.US_ASCII), null, "AAAACG", 0);
		verifyReaderBestStorage();
	}

	/** Times only the successful findBestMatch path and requires immediate counter cleanup. */
	private void stepReaderBest(final int count){
		long sum=checksum;
		for(int i=0; i<count; i++){
			final int id=readerBestReader.findBestMatch(readerBestRead, readerBestTables);
			if(id!=1){throw new AssertionError("Reader-best expected ID 1, got "+id);}
			if(readerBestReader.countArray[1]!=0 || readerBestReader.countArray[2]!=0){
				throw new AssertionError("Reader-best left dirty counters after a successful lookup: "+
					readerBestReader.countArray[1]+", "+readerBestReader.countArray[2]);
			}
			sum+=id;
		}
		checksum=sum;
	}

	/** Checks successful path checksum, literal table state and retained fixed storage. */
	private void verifyReaderBest(final long total){
		assert(total>0) : "Reader-best verification follows at least the completed warmup block.";
		verifyReaderBestStorage();
		if(readerBestReader.minHits!=2 || readerBestReader.countArray.length!=3 ||
				readerBestReader.countArray[1]!=0 || readerBestReader.countArray[2]!=0 ||
				checksum!=total){
			throw new AssertionError("Reader-best expected clean counters and checksum "+total+
				", got counters "+readerBestReader.countArray[1]+", "+
				readerBestReader.countArray[2]+" and checksum "+checksum);
		}
	}

	/** Verifies literal AAA/ACG IDs and fixed victim-free shards outside timing. */
	private void verifyReaderBestStorage(){
		if(readerBestRead.length()!=6 || readerBestTables.length!=7 ||
				readerBestArrays.length!=readerBestTables.length ||
				readerBestReader.countArray.length!=3 ||
				readerBestArrays[1].getValue(64)!=1 || readerBestArrays[0].getValue(70)!=2){
			throw new AssertionError("Reader-best literal fixture changed.");
		}
		for(int i=0; i<readerBestArrays.length; i++){
			final HashArray1D table=readerBestArrays[i];
			final int expectedSize=(i<2 ? 1 : 0);
			if(table.arrayLength()!=77 || table.victims().size()!=0 || table.size()!=expectedSize ||
					table.array()!=readerBestKeyStorage[i] || table.values()!=readerBestValueStorage[i] ||
					readerBestTables[i]!=table){
				throw new AssertionError("Reader-best shard "+i+" changed size, storage or placement.");
			}
		}
	}

	/** Installs AAA in one of seven tiny shards for the unpadded markBits path. */
	private void setupReaderMarkBits(){
		assert(readerMarkBitsReader==null) : "The reader-markbits fixture must be constructed once.";
		readerMarkBitsReader=new TableReader(3, 0, 0, 0, 0, false, false);
		readerMarkBitsReader.trimPad=0;
		readerMarkBitsTables=new AbstractKmerTable[7];
		readerMarkBitsArrays=new HashArray1D[readerMarkBitsTables.length];
		readerMarkBitsKeyStorage=new long[readerMarkBitsTables.length][];
		readerMarkBitsValueStorage=new int[readerMarkBitsTables.length][];
		for(int i=0; i<readerMarkBitsTables.length; i++){
			final HashArray1D table=new HashArray1D(new int[] {17}, -1L);
			readerMarkBitsTables[i]=readerMarkBitsArrays[i]=table;
			readerMarkBitsKeyStorage[i]=table.array();
			readerMarkBitsValueStorage[i]=table.values();
		}
		if(readerMarkBitsTables[1].set(64, 1)!=1){
			throw new AssertionError("Reader-markbits literal key AAA must install as a new ID.");
		}
		readerMarkBitsRead=new Read("CCCAAACCC".getBytes(StandardCharsets.US_ASCII), null, "CCCAAACCC", 0);
		verifyReaderMarkBitsStorage();
		verifyReaderMarkBitsResult();
	}

	/** Times a successful markBits call while preserving its normal BitSet allocation. */
	private void stepReaderMarkBits(final int count){
		long sum=checksum;
		for(int i=0; i<count; i++){
			final BitSet bits=readerMarkBitsReader.markBits(readerMarkBitsRead, readerMarkBitsTables);
			final int cardinality=(bits==null ? 0 : bits.cardinality());
			if(cardinality!=3){throw new AssertionError("Reader-markbits expected cardinality 3, got "+cardinality);}
			sum+=cardinality;
		}
		checksum=sum;
	}

	/** Checks the accumulated cardinality and exact unpadded markBits oracle outside timing. */
	private void verifyReaderMarkBits(final long total){
		assert(total>0) : "Reader-markbits verification follows at least the completed warmup block.";
		verifyReaderMarkBitsStorage();
		verifyReaderMarkBitsResult();
		if(checksum!=3L*total){
			throw new AssertionError("Reader-markbits checksum "+checksum+" does not equal 3 per call for "+
				total+" calls.");
		}
	}

	/** Checks exact bits 3..5, one stable AAA key, seven tiny shards and an unchanged read. */
	private void verifyReaderMarkBitsStorage(){
		if(readerMarkBitsRead.length()!=9 ||
				!"CCCAAACCC".equals(new String(readerMarkBitsRead.bases, StandardCharsets.US_ASCII)) ||
				readerMarkBitsTables.length!=7 || readerMarkBitsArrays.length!=readerMarkBitsTables.length ||
				readerMarkBitsReader.trimPad!=0 || readerMarkBitsArrays[1].getValue(64)!=1){
			throw new AssertionError("Reader-markbits literal fixture changed.");
		}
		for(int i=0; i<readerMarkBitsArrays.length; i++){
			final HashArray1D table=readerMarkBitsArrays[i];
			final int expectedSize=(i==1 ? 1 : 0);
			if(table.arrayLength()!=77 || table.victims().size()!=0 || table.size()!=expectedSize ||
					table.array()!=readerMarkBitsKeyStorage[i] || table.values()!=readerMarkBitsValueStorage[i] ||
					readerMarkBitsTables[i]!=table){
				throw new AssertionError("Reader-markbits shard "+i+" changed size, storage or placement.");
			}
		}
	}

	/** Checks the untimed exact markBits oracle; the measured loop checks only cardinality. */
	private void verifyReaderMarkBitsResult(){
		final BitSet bits=readerMarkBitsReader.markBits(readerMarkBitsRead, readerMarkBitsTables);
		if(bits==null || bits.cardinality()!=3 || bits.length()!=6 || !bits.get(3) ||
				!bits.get(4) || !bits.get(5)){
			throw new AssertionError("Reader-markbits expected only bits 3..5, got "+bits);
		}
		if(!"CCCAAACCC".equals(new String(readerMarkBitsRead.bases, StandardCharsets.US_ASCII))){
			throw new AssertionError("Reader-markbits markBits mutated the reused read.");
		}
	}

	/** Installs a one-shard short-key fixture whose queried windows are all present. */
	private void setupCoveragePresent(){
		assert(coverageSet==null) : "The coverage-present fixture must be constructed once.";
		coverageSet=new KmerTableSet(new String[] {"in=/dev/null", "k=3", "ways=1",
			"threads=1", "prealloc=f", "prefilter=0", "onepass=f", "rcomp=t",
			"qin=33", "interleaved=f"}, 12);
		coverageTable=new HashArray1D(new int[] {257}, -1L);
		try{
			final Field field=KmerTableSet.class.getDeclaredField("tables");
			field.setAccessible(true);
			field.set(coverageSet, new AbstractKmerTable[] {coverageTable});
		}catch(ReflectiveOperationException e){
			throw new IllegalStateException("Unable to install the bounded coverage fixture.", e);
		}
		final byte[] pattern="AAAACGTTTACGAAA".getBytes(StandardCharsets.US_ASCII);
		coverageBases=new byte[COVERAGE_BASES];
		for(int i=0; i<coverageBases.length; i++){coverageBases[i]=pattern[i%pattern.length];}
		coverageContig=new Contig(coverageBases.clone(), "coverage-present", 1);
		coverageKeyStorage=coverageTable.array();
		coverageValueStorage=coverageTable.values();
		coverageWindows=coverageBases.length-coverageSet.kbig()+1;
		seedCoveragePresent(coverageBases, COVERAGE_COUNT);
		verifyCoverageStorage();
	}

	/** Seeds every valid window from bases into the reused one-way table. */
	private void seedCoveragePresent(final byte[] bases, final int count){
		final int k=coverageSet.kbig(), shift=2*k, shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;
		int len=0;
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			final long x=dna.AminoAcid.baseToNumber[b];
			final long x2=dna.AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{len++;}
			if(len>=k){
				final long key=coverageSet.toValue(kmer, rkmer);
				final int old=coverageTable.getValue(key);
				if(old<0){
					if(coverageTable.set(key, count)!=1){
						throw new AssertionError("Coverage seed key was not new at base "+i+": "+key);
					}
					coverageKeys++;
				}else if(old!=count){
					throw new AssertionError("Coverage seed changed count for key "+key+": "+old);
				}
			}
		}
	}

	/** Times the array and Contig calcCoverage overloads without per-iteration fixture allocation. */
	private void stepCoveragePresent(final int count){
		long sum=checksum;
		for(int i=0; i<count; i++){
			final float arrayCoverage=coverageSet.calcCoverage(coverageBases, coverageBases.length);
			final float contigCoverage=coverageSet.calcCoverage(coverageContig);
			if(arrayCoverage!=COVERAGE_COUNT || contigCoverage!=COVERAGE_COUNT ||
					coverageContig.coverage!=COVERAGE_COUNT || coverageContig.maxCov!=COVERAGE_COUNT ||
					coverageContig.minCov!=COVERAGE_COUNT){
				throw new AssertionError("Coverage-present expected all fields at "+COVERAGE_COUNT+
					", got array="+arrayCoverage+", contig="+contigCoverage+", coverage="+
					coverageContig.coverage+", max="+coverageContig.maxCov+", min="+coverageContig.minCov);
			}
			sum+=(long)arrayCoverage+(long)contigCoverage+(long)coverageContig.coverage+
				coverageContig.maxCov+coverageContig.minCov;
		}
		checksum=sum;
	}

	/** Checks accumulated outputs and the retained injected table after timing. */
	private void verifyCoveragePresent(final long total){
		assert(total>0) : "Coverage verification follows at least the completed warmup block.";
		verifyCoverageStorage();
		if(checksum!=COVERAGE_CHECKSUM*total){
			throw new AssertionError("Coverage checksum "+checksum+" does not equal "+
				COVERAGE_CHECKSUM+" per paired overload call for "+total+" calls.");
		}
	}

	/** Verifies stable one-way table storage, positive counts and all-present input. */
	private void verifyCoverageStorage(){
		if(coverageBases.length!=COVERAGE_BASES || coverageContig.bases.length!=coverageBases.length ||
				coverageWindows!=COVERAGE_BASES-2 || coverageSet.kbig()!=3 || coverageSet.ways()!=1 ||
				coverageSet.tables().length!=1 || coverageSet.tables()[0]!=coverageTable ||
				coverageTable.array()!=coverageKeyStorage || coverageTable.values()!=coverageValueStorage ||
				coverageTable.arrayLength()!=317 || coverageTable.victims().size()!=0 ||
				coverageTable.size()!=coverageKeys){
			throw new AssertionError("Coverage-present changed fixed table or input storage.");
		}
		for(int i=0; i<coverageBases.length; i++){
			if(coverageContig.bases[i]!=coverageBases[i]){
				throw new AssertionError("Coverage-present mutated the reused Contig bases at "+i);
			}
		}
		int valid=0;
		final int k=coverageSet.kbig(), shift=2*k, shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0, rkmer=0;
		int len=0;
		for(int i=0; i<coverageBases.length; i++){
			final byte b=coverageBases[i];
			final long x=dna.AminoAcid.baseToNumber[b];
			final long x2=dna.AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{len++;}
			if(len>=k){
				valid++;
				final int count=coverageSet.getCount(kmer, rkmer);
				if(count!=COVERAGE_COUNT){
					throw new AssertionError("Coverage-present missing key ending at base "+i+
						": "+count+" versus "+COVERAGE_COUNT);
				}
			}
		}
		if(valid!=coverageWindows){throw new AssertionError("Coverage-present valid windows changed: "+valid);}
	}

	/** Installs a one-key long-kmer fixture whose queried windows are all present. */
	private void setupCoverageLongPresent(){
		assert(coverageLongSet==null) : "The coverage-long fixture must be constructed once.";
		Kmer.PACKED=false;
		AbstractKmerTableSet.MASK_CORE=false;
		AbstractKmerTableSet.MASK_MIDDLE=false;
		coverageLongSet=new KmerTableSetU(new String[] {"in=/dev/null", "k=62", "ways=1",
			"threads=1", "prealloc=f", "prefilter=0", "onepass=f", "rcomp=f",
			"qin=33", "interleaved=f"}, 0);
		coverageLongTable=new HashArrayU1D(new int[] {257}, coverageLongSet.k, coverageLongSet.kbig);
		try{
			final Field field=KmerTableSetU.class.getDeclaredField("tables");
			field.setAccessible(true);
			field.set(coverageLongSet, new AbstractKmerTableU[] {coverageLongTable});
		}catch(ReflectiveOperationException e){
			throw new IllegalStateException("Unable to install the bounded long-coverage fixture.", e);
		}
		coverageLongBases=new byte[COVERAGE_LONG_BASES];
		for(int i=0; i<coverageLongBases.length; i++){coverageLongBases[i]='A';}
		coverageLongContig=new Contig(coverageLongBases.clone(), "coverage-long-present", 1);
		coverageLongScratch=new Kmer(coverageLongSet.kbig);
		coverageLongKeyStorage=coverageLongArrays(coverageLongTable);
		coverageLongValueStorage=coverageLongTable.values();
		final Kmer present=new Kmer(coverageLongSet.kbig);
		for(int i=0; i<coverageLongSet.kbig; i++){present.addRight('A');}
		if(coverageLongTable.set(present, COVERAGE_COUNT)!=1){
			throw new AssertionError("Coverage-long seed key must be new.");
		}
		coverageLongWindows=coverageLongBases.length-coverageLongSet.kbig+1;
		verifyCoverageLongStorage();
	}

	/** Times the long-key array and Contig calcCoverage overloads without allocation. */
	private void stepCoverageLongPresent(final int count){
		long sum=checksum;
		for(int i=0; i<count; i++){
			final float arrayCoverage=coverageLongSet.calcCoverage(coverageLongBases,
				coverageLongBases.length, coverageLongScratch);
			final float contigCoverage=coverageLongSet.calcCoverage(coverageLongContig,
				coverageLongScratch);
			if(arrayCoverage!=COVERAGE_COUNT || contigCoverage!=COVERAGE_COUNT ||
					coverageLongContig.coverage!=COVERAGE_COUNT ||
					coverageLongContig.maxCov!=COVERAGE_COUNT ||
					coverageLongContig.minCov!=COVERAGE_COUNT){
				throw new AssertionError("Coverage-long expected all fields at "+COVERAGE_COUNT+
					", got array="+arrayCoverage+", contig="+contigCoverage+", coverage="+
					coverageLongContig.coverage+", max="+coverageLongContig.maxCov+
					", min="+coverageLongContig.minCov);
			}
			sum+=(long)arrayCoverage+(long)contigCoverage+(long)coverageLongContig.coverage+
				coverageLongContig.maxCov+coverageLongContig.minCov;
		}
		checksum=sum;
	}

	/** Checks accumulated long-key outputs and the retained injected table after timing. */
	private void verifyCoverageLongPresent(final long total){
		assert(total>0) : "Coverage-long verification follows at least the completed warmup block.";
		verifyCoverageLongStorage();
		if(checksum!=COVERAGE_CHECKSUM*total){
			throw new AssertionError("Coverage-long checksum "+checksum+" does not equal "+
				COVERAGE_CHECKSUM+" per paired overload call for "+total+" calls.");
		}
	}

	/** Verifies stable long-key storage, the one positive key and all-present input. */
	private void verifyCoverageLongStorage(){
		if(coverageLongBases.length!=COVERAGE_LONG_BASES ||
				coverageLongContig.bases.length!=coverageLongBases.length ||
				coverageLongWindows!=COVERAGE_LONG_BASES-61 || coverageLongSet.kbig()!=62 ||
				coverageLongSet.ways()!=1 || coverageLongSet.tables().length!=1 ||
				coverageLongSet.tables()[0]!=coverageLongTable ||
				coverageLongArrays(coverageLongTable)!=coverageLongKeyStorage ||
				coverageLongTable.values()!=coverageLongValueStorage ||
				coverageLongTable.arrayLength()!=317 || coverageLongTable.victims().size()!=0 ||
				coverageLongTable.size()!=1){
			throw new AssertionError("Coverage-long changed fixed table or input storage.");
		}
		coverageLongScratch.clear();
		for(int i=0; i<coverageLongSet.kbig; i++){coverageLongScratch.addRight('A');}
		if(coverageLongScratch.len()!=coverageLongSet.kbig ||
				coverageLongSet.getCount(coverageLongScratch)!=COVERAGE_COUNT){
			throw new AssertionError("Coverage-long missing the stored all-A long key.");
		}
		for(int i=0; i<coverageLongBases.length; i++){
			if(coverageLongBases[i]!='A' || coverageLongContig.bases[i]!=coverageLongBases[i]){
				throw new AssertionError("Coverage-long mutated the reused input at "+i);
			}
		}
	}

	/** Returns borrowed long-key primary key storage for replacement checks. */
	private long[][] coverageLongArrays(final HashArrayU table){
		try{
			final Field field=HashArrayU.class.getDeclaredField("arrays");
			field.setAccessible(true);
			return (long[][])field.get(table);
		}catch(ReflectiveOperationException e){
			throw new IllegalStateException("Unable to inspect the bounded long-coverage fixture.", e);
		}
	}

	/** Checks literal placement, unchanged storage, count oracles and value tails outside timing. */
	private void verifyHybrid(final long total){
		assert(total>=0) : "Hybrid expected counts include only completed warmup and measured updates.";
		final long primary=hybridVictim ? 60 : KEYS, victims=hybridVictim ? KEYS : 0;
		if(hybridCreated!=0 || hybridTable.size()!=primary || hybridTable.victims().size()!=victims ||
				hybridTable.arrayLength()!=100063 || hybridTable.array()!=hybridKeyStorage ||
				hybridTable.victims().array()!=hybridVictimStorage){
			throw new AssertionError("Hybrid kernel changed its precreated storage or reported a new key.");
		}
		for(int i=0; i<hybridKeyStorage.length; i++){
			final long expected=i<primary ? (hybridVictim ? 100003L*i : i) : -1;
			if(hybridKeyStorage[i]!=expected){throw new AssertionError("Unexpected primary key at cell "+i);}
			if(hybridVictim && i<primary && hybridTable.getValue(expected)!=1){
				throw new AssertionError("An untimed primary filler changed at cell "+i);
			}
		}
		long sum=0;
		final int[] singleton=new int[1];
		for(int i=0; i<KEYS; i++){
			final long updates=total/KEYS+(i<total%KEYS ? 1 : 0);
			final long expected=legacyHybrid && hybridVictim ? 1 : 1+updates;
			if(expected>Integer.MAX_VALUE){throw new AssertionError("Hybrid work count exceeded the unsaturated oracle.");}
			final int actual=hybridTable.getValue(keys[i]);
			final int[] values=hybridTable.getValues(keys[i], singleton);
			if(actual!=expected || values==null || values.length<1 || values[0]!=expected){
				throw new AssertionError("Hybrid count mismatch at key "+i+": "+actual+" versus "+expected);
			}
			if(hybridVictim){
				final int tail=legacyHybrid && updates>0 ? 2 : -1;
				if(values!=hybridValueStorage[i] || values.length!=2 || values[1]!=tail ||
						hybridTable.victims().getValue(keys[i])!=expected){
					throw new AssertionError("Unexpected victim placement, allocation or tail at key "+i);
				}
			}else if(values!=singleton || hybridTable.victims().contains(keys[i])){
				throw new AssertionError("Primary count must remain direct scalar storage at key "+i);
			}
			sum+=(i+1L)*actual;
		}
		if(hybridTable.contains(hybridVictim ? 1 : KEYS)){
			throw new AssertionError("Hybrid kernel created a key outside the fixed input set.");
		}
		checksum=sum;
	}

	/** Checks the steady merge result outside timing, independently of borrowed input-array contents. */
	private void verifyMerge(){
		if(mergeTable.size()!=KEYS || mergeCreated!=0){
			throw new AssertionError("Existing-key merges must retain all keys and report no creation: "+
				mergeTable.size()+", "+mergeCreated);
		}
		long sum=0;
		for(int i=0; i<KEYS; i++){
			final int[] actual=mergeTable.getValues(i, null);
			if(actual==null || actual.length<MERGE_VALUES){throw new AssertionError("Missing values at key "+i);}
			for(int j=0; j<MERGE_VALUES; j++){
				final int expected=i*MERGE_VALUES+j+1;
				if(actual[j]!=expected){
					throw new AssertionError("Changed merge value at key "+i+", index "+j+": "+actual[j]+" versus "+expected);
				}
				sum+=(i+1L)*actual[j];
			}
			for(int j=MERGE_VALUES; j<actual.length; j++){
				if(actual[j]>=0){throw new AssertionError("Unexpected merge value at key "+i+", index "+j);}
			}
		}
		checksum=sum;
	}

	/** Deterministic key set; no allocation occurs in the update loops. */
	private final long[] keys=new long[KEYS];
	/** Backend tables, absent for the reversal kernel. */
	private HashArray1D[] tables;
	/** Thread-confined write buffer for the two buffered kernels. */
	private HashBuffer buffer;
	/** Reused rolling key for canonical reversal. */
	private Kmer rolling;
	/** Nonresizing multi-value forest used only by the steady merge kernel. */
	private HashForest mergeTable;
	/** Preallocated four-value lists with negative padding; forest nodes retain these arrays. */
	private int[][] mergeValues;
	/** Bitwise accumulation preserves any unexpected nonzero creation return for the untimed check. */
	private int mergeCreated;
	/** One nonresizing hybrid table, absent for all older modes. */
	private HashArray hybridTable;
	/** Selects an all-victim timed key set, with 60 untimed primary fillers. */
	private boolean hybridVictim;
	/** Enables only the explicitly requested known-wrong hybrid victim baseline oracle. */
	private boolean legacyHybrid;
	/** Retained key and bucket arrays detect unexpected replacement outside timed updates. */
	private long[] hybridKeyStorage;
	private KmerNode[] hybridVictimStorage;
	/** Initial two-slot victim arrays; warmup and measurement must retain these exact objects. */
	private int[][] hybridValueStorage;
	/** OR of creation returns; every timed and warmup update must report zero. */
	private int hybridCreated;
	/** Stable fresh scalar table walked completely once per measured iteration. */
	private HashArray1D walkerTable;
	/** Retained arrays detect unexpected resizing or regeneration in the walker fixture. */
	private long[] walkerKeyStorage;
	private int[] walkerValueStorage;
	/** Successful findBestMatch fixture with fixed shards and a reused query read. */
	private TableReader readerBestReader;
	private AbstractKmerTable[] readerBestTables;
	private HashArray1D[] readerBestArrays;
	private Read readerBestRead;
	/** Retained per-shard storage for the reader-best fixture. */
	private long[][] readerBestKeyStorage;
	private int[][] readerBestValueStorage;
	/** Successful markBits fixture with fixed shards and a reused query read. */
	private TableReader readerMarkBitsReader;
	private AbstractKmerTable[] readerMarkBitsTables;
	private HashArray1D[] readerMarkBitsArrays;
	private Read readerMarkBitsRead;
	/** Retained per-shard storage for the reader-markbits fixture. */
	private long[][] readerMarkBitsKeyStorage;
	private int[][] readerMarkBitsValueStorage;
	/** One injected short-key table set for all-present calcCoverage calls. */
	private KmerTableSet coverageSet;
	private HashArray1D coverageTable;
	/** Stable inputs for the two calcCoverage overloads. */
	private byte[] coverageBases;
	private Contig coverageContig;
	/** Retained table arrays detect unexpected replacement during coverage timing. */
	private long[] coverageKeyStorage;
	private int[] coverageValueStorage;
	/** Number of valid windows in the all-present coverage fixture. */
	private int coverageWindows;
	/** Number of unique canonical keys seeded in the fixed coverage table. */
	private int coverageKeys;
	/** One injected long-key table set for all-present calcCoverage calls. */
	private KmerTableSetU coverageLongSet;
	private HashArrayU1D coverageLongTable;
	/** Stable inputs and scratch for the two long-key calcCoverage overloads. */
	private byte[] coverageLongBases;
	private Contig coverageLongContig;
	private Kmer coverageLongScratch;
	/** Retained long-key table arrays detect unexpected replacement during timing. */
	private long[][] coverageLongKeyStorage;
	private int[] coverageLongValueStorage;
	/** Number of valid windows in the all-present long-key coverage fixture. */
	private int coverageLongWindows;
	/** Actual long-key table set with one injected fixed-size shard. */
	private KmerTableSetU longFillSet;
	private HashArrayU1D longFillTable;
	/** Immutable input sequence and reusable four-entry count output. */
	private byte[] longFillBases;
	private IntList longFillCounts;
	private int[] longFillCountStorage;
	/** One reused rolling key; retained arrays detect unexpected backing replacement. */
	private Kmer longFillScratch;
	private long[] longFillForward, longFillReverse;
	/** Last valid-window return, also included in every timed checksum contribution. */
	private int longFillValid;
	/** Retained reduction so results remain observable after timing. */
	private long checksum;
	private static volatile long blackhole;
	private static final long[] WALKER_KEYS={1, 2, 3, 4};
	private static final int[] WALKER_COUNTS={3, 5, 7, 11};
	private static final int WALKER_RECORDS=4;
	private static final long WALKER_CHECKSUM=248666L;
	private static final long WALKER_KEY_WEIGHT=1009L;
	private static final long WALKER_VALUE_WEIGHT=9176L;
	private static final int COVERAGE_COUNT=7;
	private static final int COVERAGE_BASES=128;
	private static final int COVERAGE_LONG_BASES=128;
	private static final long COVERAGE_CHECKSUM=COVERAGE_COUNT*5L;
	private static final ThreadMXBean CPU_CLOCK=ManagementFactory.getThreadMXBean();
	private static final int KEYS=4096;
	private static final int MERGE_VALUES=4;
	private static final int BLOCK=1048576;
	private static final long TARGET_NANOS=15000000000L;
}
