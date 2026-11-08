package test.bloom;

import java.io.InputStream;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.Properties;
import java.util.concurrent.TimeUnit;

import bloom.KCountArray;
import bloom.KCountArray7MTA;
import bloom.ReadCounter;
import shared.Shared;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;

/** Checks extra-only and finite primary/extra counting through real readers and the public makeKca API.
 * The unlimited cases preserve the original extra-only repair. Finite cases require
 * maxReads to limit only the primary input while every extra input is read completely
 * on the first, intermediate and final passes.
 * Each fresh child counts exact-k records, retaining raw failures and every observed table.
 * Observation delegates to real count methods without injecting counts.
 * Parent assertion mode is inherited by every child.
 * @author Brian Bushnell, Sucrose, Jean */
public final class ReadCounterExtraInputTest{

	/** Runs every case even when an earlier child's literal success expectation fails. */
	public static void main(final String[] args) throws Exception{
		if(args.length==2 && args[0].startsWith("child=") && args[1].startsWith("outdir=")){
			child(args[0].substring(6), Paths.get(args[1].substring(7)));
			return;
		}
		if(args.length!=1 || !args[0].startsWith("outdir=")){
			throw new IllegalArgumentException("Usage: outdir=<new directory>");
		}
		final Path root=Paths.get(args[0].substring(7)).toAbsolutePath();
		check(!Files.exists(root), "Fresh output prevents old results from satisfying a new capture.");
		Files.createDirectories(root);
		int failures=0;
		for(String name : CASES){failures+=runChild(root, name);}
		if(failures>0){throw new AssertionError("ReadCounter extra-input contracts failed "+failures+" of "+CASES.length+" cases.");}
		System.out.println("READ_COUNTER_EXTRA_INPUT_TEST_PASS cases="+CASES.length);
	}

	/** Bounds only the owned child; destruction is always a failed natural-exit check. */
	private static int runChild(final Path root, final String name) throws Exception{
		final Path dir=Files.createDirectory(root.resolve(name));
		final ArrayList<String> command=new ArrayList<String>();
		command.add(Paths.get(System.getProperty("java.home"), "bin", "java").toString());
		command.add(ReadCounterExtraInputTest.class.desiredAssertionStatus() ? "-ea" : "-da");
		command.add("-Xmx256m");
		command.add("-XX:+ExitOnOutOfMemoryError");
		try{
			Class.forName("jdk.incubator.vector.ByteVector", false, ReadCounterExtraInputTest.class.getClassLoader());
			command.add("--add-modules=jdk.incubator.vector");
		}catch(ClassNotFoundException | LinkageError unavailable){
			//Do not request an optional module that the parent cannot load.
		}
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(ReadCounterExtraInputTest.class.getName());
		command.add("child="+name);
		command.add("outdir="+dir);
		Files.write(dir.resolve("command.txt"), command, StandardCharsets.UTF_8);
		final Process process=new ProcessBuilder(command).redirectInput(Paths.get("/dev/null").toFile())
				.redirectOutput(dir.resolve("stdout.txt").toFile()).redirectError(dir.resolve("stderr.txt").toFile()).start();
		boolean completed=false;
		try{completed=process.waitFor(15, TimeUnit.SECONDS);}
		finally{
			if(process.isAlive()){
				process.destroy();
				if(!process.waitFor(2, TimeUnit.SECONDS)){
					process.destroyForcibly();
					check(process.waitFor(5, TimeUnit.SECONDS), "Fixture-owned child must terminate after forced destruction.");
				}
			}
		}
		final Properties result=new Properties();
		if(Files.exists(dir.resolve("result.properties"))){
			try(InputStream input=Files.newInputStream(dir.resolve("result.properties"))){result.load(input);}
		}
		final int exit=process.exitValue();
		final boolean pass=completed && exit==0 && "true".equals(result.getProperty("pass")) &&
				name.equals(result.getProperty("case")) &&
				Boolean.toString(ReadCounterExtraInputTest.class.desiredAssertionStatus()).equals(result.getProperty("assertions"));
		final String summary="case="+name+" completed="+completed+" exit="+exit+" result="+result+" pass="+pass;
		Files.write(dir.resolve("outcome.txt"), Arrays.asList(summary), StandardCharsets.UTF_8);
		System.out.println(summary);
		return pass ? 0 : 1;
	}

	/** Keeps inputs, failure causes, observed tables and restored settings even on the broken baseline. */
	private static void child(final String name, final Path dir) throws Exception{
		check(Arrays.asList(CASES).contains(name), "Unknown child case: "+name);
		final CaseSpec spec=caseSpec(name);
		final Path primaryPath=spec.primarySequences==null ? null : dir.resolve("primary.fq");
		final Path extraPath=spec.extraSequences==null ? null : dir.resolve("extra.fq");
		final byte[] primaryInput=writeFastq(primaryPath, spec.primarySequences);
		final byte[] extraInput=writeFastq(extraPath, spec.extraSequences);
		Files.write(dir.resolve("oracle.tsv"), spec.oracleRows(), StandardCharsets.UTF_8);
		final Properties result=new Properties();
		result.setProperty("case", name);
		result.setProperty("assertions", Boolean.toString(ReadCounterExtraInputTest.class.desiredAssertionStatus()));
		result.setProperty("passes", Integer.toString(spec.passes));
		result.setProperty("requestedLimit", Long.toString(spec.maxReads));
		result.setProperty("expectedTotal", Integer.toString(spec.expectedTotal));
		result.setProperty("inputRole", spec.inputRole());
		result.setProperty("primaryRecords", Integer.toString(spec.primarySequences==null ? 0 : spec.primarySequences.length));
		result.setProperty("extraRecords", Integer.toString(spec.extraSequences==null ? 0 : spec.extraSequences.length));
		result.setProperty("settings", "k=5 cells=4096 bits=8 hashes=1 seed_mask=0 threads=1 qin=33 minqual=0 minprob=0 thresholds=1,1 detectstepsize=1 canonical=true duplicates=true buffered=false sketch=false sampling=0 idmask=0 overlap=false prefilter=false dual_hash=false");
		final String primaryString=primaryPath==null ? null : primaryPath.toString();
		final String extraString=extraPath==null ? null : extraPath.toString();
		final ObservedCounter counter=new ObservedCounter(dir, result, spec, primaryString, extraString);
		KCountArray returned=null;
		Throwable failure=null;
		try{
			configure();
			returned=counter.makeKca(primaryString, null,
					extraString==null ? null : Arrays.asList(extraString), 8, 4096, 1, 0,
					spec.maxReads, spec.passes, 1, 1, null, 0);
			check(returned!=null && returned==counter.lastTable && counter.calls==spec.callRoles.length &&
					counter.completedPasses==spec.passes,
					"Every expected pass must count real input and return the final table; calls="+
					counter.calls+" completedPasses="+counter.completedPasses);
			check(ReadCounter.readsProcessed==spec.expectedTotal &&
					ReadCounter.keysCounted==spec.expectedTotal &&
					ReadCounter.increments==spec.expectedTotal,
					"Exact-k records require "+spec.expectedTotal+" reads/keys/increments; got "+
					ReadCounter.readsProcessed+"/"+ReadCounter.keysCounted+"/"+ReadCounter.increments);
			check(counter.allTablesMatch,
					"Each observed count call must leave exactly the expected literal depths in its table.");
		}catch(Throwable caught){failure=caught;}
		try{
			if(counter.lastTable!=null){
				final boolean exact=saveTable(counter.lastTable, dir.resolve("last-table.tsv"),
						result, "last", spec.finalSequences(), spec.absentSequences);
				if(returned!=null){check(exact, "Finalization must preserve the expected literal depth1 counts.");}
			}
		}catch(Throwable observed){failure=combine(failure, observed);}
		final boolean restored=ReadCounter.minQuality==7 && ReadCounter.maxReads==123 && !FastaReadInputStream.SPLIT_READS;
		result.setProperty("settingsRestored", Boolean.toString(restored));
		result.setProperty("restoredMinQuality", Byte.toString(ReadCounter.minQuality));
		result.setProperty("restoredMaxReads", Long.toString(ReadCounter.maxReads));
		result.setProperty("restoredFastaSplit", Boolean.toString(FastaReadInputStream.SPLIT_READS));
		try{
			check(restored, "makeKca must restore the caller's minQuality7, maxReads123 and FASTA split=false on every exit.");
		}catch(Throwable observed){failure=combine(failure, observed);}
		try{
			final boolean primaryPreserved=primaryPath==null || Arrays.equals(primaryInput, Files.readAllBytes(primaryPath));
			final boolean extraPreserved=extraPath==null || Arrays.equals(extraInput, Files.readAllBytes(extraPath));
			result.setProperty("primaryPreserved", Boolean.toString(primaryPreserved));
			result.setProperty("extraPreserved", Boolean.toString(extraPreserved));
			result.setProperty("inputsPreserved", Boolean.toString(primaryPreserved && extraPreserved));
			check(primaryPreserved && extraPreserved, "Counting must not change the saved literal FASTQ bytes.");
		}catch(Throwable observed){failure=combine(failure, observed);}
		result.setProperty("tableReturned", Boolean.toString(returned!=null));
		result.setProperty("countCalls", Integer.toString(counter.calls));
		result.setProperty("completedCalls", Integer.toString(counter.completedPasses));
		result.setProperty("readsProcessed", Long.toString(ReadCounter.readsProcessed));
		result.setProperty("keysCounted", Long.toString(ReadCounter.keysCounted));
		result.setProperty("increments", Long.toString(ReadCounter.increments));
		result.setProperty("pass", Boolean.toString(failure==null));
		result.setProperty("failure", failure==null ? "none" : failure.toString());
		try(OutputStream output=Files.newOutputStream(dir.resolve("result.properties"))){
			result.store(output, "Literal expectations unchanged across baseline and fixed captures");
		}
		if(failure!=null){throw new AssertionError("Extra-input contract failed: "+name, failure);}
		//No System.exit or reader cleanup here: leaked non-daemon workers must remain visible.
	}

	/** Writes one exact-k FASTQ record per sequence, or no file for an absent input role. */
	private static byte[] writeFastq(final Path path, final String[] sequences) throws Exception{
		if(sequences==null){return null;}
		final StringBuilder fastq=new StringBuilder();
		for(int i=0; i<sequences.length; i++){
			check(sequences[i].length()==K,
					"Each literal record must contain one exact k5 window: "+sequences[i]);
			canonical(sequences[i]);
			fastq.append('@').append("record").append(i).append('\n').append(sequences[i]).append("\n+\nIIIII\n");
		}
		final byte[] input=fastq.toString().getBytes(StandardCharsets.US_ASCII);
		Files.write(path, input);
		return input;
	}

	/** Records a sparse complete histogram (omitted bins are zero); never finalizes or changes the table. */
	private static boolean saveTable(final KCountArray table, final Path path, final Properties result,
			final String prefix, final String[] expectedSequences, final String[] absentSequences) throws Exception{
		final long[] expectedKeys=keys(expectedSequences);
		final long occupied=table.cellsUsed(1);
		final long[] histogram=table.transformToFrequency();
		final ArrayList<String> rows=new ArrayList<String>();
		boolean exact=occupied==expectedKeys.length && histogram.length>1 && histogram[1]==expectedKeys.length;
		rows.add("depth\tcells");
		for(int i=0; i<histogram.length; i++){if(histogram[i]!=0){rows.add(i+"\t"+histogram[i]);}}
		for(int i=0; i<expectedKeys.length; i++){
			final int depth=table.read(expectedKeys[i]);
			result.setProperty(prefix+".depth"+expectedKeys[i], Integer.toString(depth));
			exact&=depth==1;
		}
		for(int i=0; i<absentSequences.length; i++){
			final long key=canonical(absentSequences[i]);
			final int depth=table.read(key);
			result.setProperty(prefix+".absentDepth"+key, Integer.toString(depth));
			exact&=depth==0;
		}
		Files.write(path, rows, StandardCharsets.UTF_8);
		result.setProperty(prefix+".occupiedCells", Long.toString(occupied));
		result.setProperty(prefix+".allocatedCells", Long.toString(table.cells));
		result.setProperty(prefix+".histogramBins", Integer.toString(histogram.length));
		result.setProperty(prefix+".cellBits", Integer.toString(table.cellBits));
		result.setProperty(prefix+".literalCountsMatch", Boolean.toString(exact));
		return exact;
	}

	/** Observes actual count dispatch and completed tables; all input/counting work remains in super. */
	private static final class ObservedCounter extends ReadCounter{
		ObservedCounter(final Path dir_, final Properties result_, final CaseSpec spec_,
				final String primaryInput_, final String extraInput_){
			super(K, true, false, false, false, false);
			dir=dir_;
			result=result_;
			spec=spec_;
			primaryInput=primaryInput_;
			extraInput=extraInput_;
			detectStepsize=1;
		}

		@Override
		public KCountArray count(final String first, final String second, final KCountArray table) throws Exception{
			before(first, second, table, null, maxReads);
			final KCountArray counted=super.count(first, second, table);
			after(table, counted);
			return counted;
		}

		@Override
		public KCountArray count(final String first, final String second, final int bits,
				final KCountArray table, final KCountArray trusted, final long limit,
				final int threshold, final boolean conservative) throws Exception{
			before(first, second, table, trusted, limit);
			check(bits==8 && threshold==1, "Real multipass calls must preserve the fixture cell width and singleton threshold.");
			final KCountArray counted=super.count(first, second, bits, table, trusted, limit, threshold, conservative);
			after(table, counted);
			return counted;
		}

		/** Verifies call ordering and borrowed tables without pre-rejecting the old finite-extra limit bug. */
		private void before(final String first, final String second, final KCountArray table,
				final KCountArray trusted, final long limit){
			calls++;
			check(calls<=spec.callRoles.length,
					"makeKca must issue only the expected primary/extra count calls for "+spec.name+
					"; call="+calls+" expected="+spec.callRoles.length);
			final String role=spec.callRoles[calls-1];
			final String expectedInput=PRIMARY.equals(role) ? primaryInput : extraInput;
			final long expectedStaticLimit=spec.expectedStaticMaxReads[calls-1];
			result.setProperty("call"+calls+".role", role);
			result.setProperty("call"+calls+".primary", first==null ? "<null>" : first);
			result.setProperty("call"+calls+".limit", Long.toString(limit));
			result.setProperty("call"+calls+".expectedStaticLimit", Long.toString(expectedStaticLimit));
			check(expectedInput!=null && expectedInput.equals(first) && second==null &&
					minQuality==0 && maxReads==expectedStaticLimit,
					"Call "+calls+" must read the expected "+role+
					" input with no mate and the wrapper static settings.");
			if(PRIMARY.equals(role)){
				check(limit==spec.maxReads,
						"Primary calls must keep the explicit maxReads limit; limit="+limit+
						" expected="+spec.maxReads);
			}
			check(table!=null && table.cellBits==8 && table.maxValue==255,
					"Every count call needs an eight-bit destination table.");
			if(spec.newTable[calls-1]){
				check(table!=lastTable, "Call "+calls+" must start a fresh pass table.");
			}else{
				check(table==lastTable, "Call "+calls+" must continue the current pass table.");
			}
			check(trusted==completedPassTable,
					"Call "+calls+" may trust only the preceding completed pass table.");
			lastTable=table;
		}

		/** Records each cumulative call table and marks a pass complete only after its last input. */
		private void after(final KCountArray table, final KCountArray counted) throws Exception{
			final boolean exact=saveTable(table, dir.resolve("call"+calls+"-table.tsv"),
					result, "call"+calls, spec.expectedAfterCalls[calls-1], spec.absentSequences);
			check(counted==table,
					"Each actual count call must return the borrowed destination table.");
			allTablesMatch&=exact;
			if(spec.endOfPass[calls-1]){
				completedPasses++;
				completedPassTable=table;
			}
		}

		private final Path dir;
		private final Properties result;
		private final CaseSpec spec;
		private final String primaryInput;
		private final String extraInput;
		private int calls, completedPasses;
		private boolean allTablesMatch=true;
		private KCountArray lastTable, completedPassTable;
	}

	/** Holds one public makeKca vector and the cumulative table expected after each count call. */
	private static final class CaseSpec{
		CaseSpec(final String name_, final String[] primarySequences_, final String[] extraSequences_,
				final long maxReads_, final int passes_, final String[] callRoles_,
				final String[][] expectedAfterCalls_, final String[] absentSequences_,
				final long[] expectedStaticMaxReads_, final boolean[] newTable_,
				final boolean[] endOfPass_){
			name=name_;
			primarySequences=primarySequences_;
			extraSequences=extraSequences_;
			maxReads=maxReads_;
			passes=passes_;
			callRoles=callRoles_;
			expectedAfterCalls=expectedAfterCalls_;
			absentSequences=absentSequences_;
			expectedStaticMaxReads=expectedStaticMaxReads_;
			newTable=newTable_;
			endOfPass=endOfPass_;
			check(callRoles.length==expectedAfterCalls.length &&
					callRoles.length==expectedStaticMaxReads.length &&
					callRoles.length==newTable.length && callRoles.length==endOfPass.length,
					"Case vectors must describe every expected count call for "+name);
			expectedTotal=expectedTotal();
			for(String[] sequences : expectedAfterCalls){check(keysDistinct(concat(sequences, absentSequences)),
					"Literal oracles must have distinct canonical expected and absent keys for "+name);}
		}

		/** Returns the final cumulative exact-k sequences expected in the returned table. */
		String[] finalSequences(){return expectedAfterCalls[expectedAfterCalls.length-1];}

		/** Returns a compact role label for result.properties. */
		String inputRole(){
			return primarySequences==null ? "extra-only" : (extraSequences==null ? "primary" : "primary+extra");
		}

		/** Returns tabular final-depth oracle rows for baseline and fixed artifacts. */
		ArrayList<String> oracleRows(){
			final String[] sequences=finalSequences();
			final long[] expectedKeys=keys(sequences);
			final ArrayList<String> rows=new ArrayList<String>();
			rows.add("sequence\tcanonical_key\tfinal_depth");
			for(int i=0; i<sequences.length; i++){rows.add(sequences[i]+"\t"+expectedKeys[i]+"\t1");}
			return rows;
		}

		/** Sums per-call input cardinality from the cumulative expected tables. */
		private int expectedTotal(){
			int total=0, previous=0;
			for(int i=0; i<expectedAfterCalls.length; i++){
				if(newTable[i]){previous=0;}
				total+=expectedAfterCalls[i].length-previous;
				previous=expectedAfterCalls[i].length;
			}
			return total;
		}

		private final String name;
		private final String[] primarySequences;
		private final String[] extraSequences;
		private final long maxReads;
		private final int passes;
		private final int expectedTotal;
		private final String[] callRoles;
		private final String[][] expectedAfterCalls;
		private final String[] absentSequences;
		private final long[] expectedStaticMaxReads;
		private final boolean[] newTable;
		private final boolean[] endOfPass;
	}

	/** Sets independent child globals; caller sentinels test makeKca restoration without resetting counts. */
	private static void configure(){
		Locale.setDefault(Locale.ROOT);
		Shared.setThreads(1);
		FASTQ.ASCII_OFFSET=33;
		FASTQ.DETECT_QUALITY=FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		Read.VALIDATE_IN_CONSTRUCTOR=Read.CHANGE_QUALITY=true;
		Read.setMinCalledQuality(2);
		Read.setMaxCalledQuality(50);
		KCountArray7MTA.setSeed(0);
		ReadCounter.maxShortKmerLength=31;
		ReadCounter.MAX_COUNT_THREADS=1;
		ReadCounter.minProb=0;
		ReadCounter.KEEP_DUPLICATE_KMERS=ReadCounter.CANONICAL=true;
		ReadCounter.KMERS_PER_READ=ReadCounter.IDMASK=0;
		ReadCounter.BUFFERED=ReadCounter.SKETCH_MODE=ReadCounter.STORE_HASHED=ReadCounter.PREJOIN=false;
		ReadCounter.minQuality=7;
		ReadCounter.maxReads=123;
		FastaReadInputStream.SPLIT_READS=false;
		check(Read.capQuality((byte)40)==40, "Normal defined-base validation must preserve the literal Q40 qualities.");
		check(ReadCounter.keysCounted==0 && ReadCounter.increments==0 && ReadCounter.readsProcessed==0,
				"A fresh child must begin with zero totals; resetting would conceal state contamination.");
	}

	/** Builds the vector matching a child name; each vector keeps its own baseline artifacts. */
	private static CaseSpec caseSpec(final String name){
		if(name.equals("extra-one-pass")){
			return new CaseSpec(name, null, TWO_RECORDS, -1, 1,
					new String[]{EXTRA}, new String[][]{TWO_RECORDS}, EMPTY,
					new long[]{-1}, new boolean[]{true}, new boolean[]{true});
		}
		if(name.equals("extra-two-pass")){
			return new CaseSpec(name, null, TWO_RECORDS, -1, 2,
					new String[]{EXTRA, EXTRA}, new String[][]{TWO_RECORDS, TWO_RECORDS}, EMPTY,
					new long[]{-1, -1}, new boolean[]{true, true}, new boolean[]{true, true});
		}
		if(name.equals("primary-two-pass")){
			return new CaseSpec(name, TWO_RECORDS, null, -1, 2,
					new String[]{PRIMARY, PRIMARY}, new String[][]{TWO_RECORDS, TWO_RECORDS}, EMPTY,
					new long[]{-1, -1}, new boolean[]{true, true}, new boolean[]{true, true});
		}
		if(name.equals("primary-extra-one-pass-finite")){
			return new CaseSpec(name, LIMITED_PRIMARY, LIMITED_EXTRA, 1, 1,
					new String[]{PRIMARY, EXTRA},
					new String[][]{LIMITED_PRIMARY_FIRST, LIMITED_FINAL}, LIMITED_PRIMARY_ABSENT,
					new long[]{1, -1}, new boolean[]{true, false}, new boolean[]{false, true});
		}
		if(name.equals("extra-two-pass-finite")){
			return new CaseSpec(name, null, LIMITED_EXTRA, 1, 2,
					new String[]{EXTRA, EXTRA}, new String[][]{LIMITED_EXTRA, LIMITED_EXTRA}, EMPTY,
					new long[]{-1, -1}, new boolean[]{true, true}, new boolean[]{true, true});
		}
		if(name.equals("primary-extra-two-pass-finite")){
			return new CaseSpec(name, LIMITED_PRIMARY, LIMITED_EXTRA, 1, 2,
					new String[]{PRIMARY, EXTRA, PRIMARY, EXTRA},
					new String[][]{LIMITED_PRIMARY_FIRST, LIMITED_FINAL,
							LIMITED_PRIMARY_FIRST, LIMITED_FINAL}, LIMITED_PRIMARY_ABSENT,
					new long[]{1, -1, -1, -1}, new boolean[]{true, false, true, false},
					new boolean[]{false, true, false, true});
		}
		throw new IllegalArgumentException("Unknown child case: "+name);
	}

	/** Derives all canonical two-bit keys directly from literal DNA. */
	private static long[] keys(final String[] sequences){
		final long[] keys=new long[sequences.length];
		for(int i=0; i<sequences.length; i++){keys[i]=canonical(sequences[i]);}
		return keys;
	}

	/** Requires exact-depth oracles to use collision-free literal windows. */
	private static boolean keysDistinct(final String[] sequences){
		final long[] keys=keys(sequences);
		for(int i=0; i<keys.length; i++){
			for(int j=i+1; j<keys.length; j++){
				if(keys[i]==keys[j]){return false;}
			}
		}
		return true;
	}

	/** Joins two sequence arrays without exposing mutable vector construction to the cases. */
	private static String[] concat(final String[] first, final String[] second){
		final String[] joined=Arrays.copyOf(first, first.length+second.length);
		System.arraycopy(second, 0, joined, first.length, second.length);
		return joined;
	}

	/** Derives the two-bit forward/reverse keys directly from the literal DNA, without production hash helpers. */
	private static long canonical(final String sequence){
		long forward=0, reverse=0;
		for(int i=0; i<sequence.length(); i++){
			final int code="ACGT".indexOf(sequence.charAt(i));
			check(code>=0, "The fixed exact-k oracle contains only defined DNA.");
			forward=(forward<<2)|code;
			reverse|=(long)(3-code)<<(2*i);
		}
		return Math.max(forward, reverse);
	}

	/** Preserves the first failure and attaches later observation failures rather than replacing it. */
	private static Throwable combine(final Throwable original, final Throwable observed){
		if(original==null){return observed;}
		original.addSuppressed(observed);
		return original;
	}

	/** Keeps all oracle and setup checks active with assertions enabled or disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final int K=5;
	private static final String PRIMARY="primary", EXTRA="extra";
	private static final String[] TWO_RECORDS={"AAACC", "AACCC"};
	private static final String[] EMPTY={};
	private static final String[] LIMITED_PRIMARY={"AAACC", "AACCC", "CCGTA"};
	private static final String[] LIMITED_PRIMARY_FIRST={LIMITED_PRIMARY[0]};
	private static final String[] LIMITED_PRIMARY_ABSENT={LIMITED_PRIMARY[1], LIMITED_PRIMARY[2]};
	private static final String[] LIMITED_EXTRA={"CCCAA", "ACGTA", "TGCAT"};
	private static final String[] LIMITED_FINAL=concat(LIMITED_PRIMARY_FIRST, LIMITED_EXTRA);
	private static final String[] CASES={"extra-one-pass", "extra-two-pass", "primary-two-pass",
			"primary-extra-one-pass-finite", "extra-two-pass-finite", "primary-extra-two-pass-finite"};
}
