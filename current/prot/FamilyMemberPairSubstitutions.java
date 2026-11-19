package prot;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Locale;
import java.util.SplittableRandom;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.LongHashSet;
import parse.Parse;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import template.ThreadWaiter;

/** Samples distinct member pairs and retains their actual POSSUM40 alignment columns.
 * Emits both span-gated and ungated descriptive counts with Brian-selected LONG70/SHORT28 span gates as the primary population.
 * Every output is final study evidence; no dimer counting or matrix derivation occurs here.
 * @author Brian, Fischl
 */
public final class FamilyMemberPairSubstitutions{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private FamilyMemberPairSubstitutions(){}

	/** Parses the normal BBTools pre-parser contract and fails the whole process on worker errors. */
	public static void main(String[] args){
		try{
			final PreParser pp=new PreParser(args, FamilyMemberPairSubstitutions.class, false);
			try{run(pp.args);}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){
			failure.printStackTrace();
			System.exit(1);
		}
	}

	/** Loads one accepted bucket at a time; selected-pair memory never grows with the population. */
	public static void run(String[] args) throws Exception{
		final long started=System.nanoTime();
		final HashMap<String, String> options=parse(args);
		final String mode=options.containsKey("mode") ? options.get("mode") : "count";
		require(mode.equals("count") || mode.equals("pilot"), "Mode must be count or an explicit timing-only pilot");
		final boolean pilot=mode.equals("pilot");
		require(pilot==options.containsKey("pilotbuckets"), "Only pilot mode requires an explicit selected bucket list");
		final int[] buckets=selectedBuckets(integer(options, "buckets"), pilot ? required(options, "pilotbuckets") : null);
		final String view=required(options, "view");
		require(view.equals("long") || view.equals("short"), "View must identify the accepted LONG or SHORT population");
		require(required(options, "verification").equals("structural"), "This recipe explicitly inherits accepted input byte identity");
		final int threads=integer(options, "t"), span=view.equals("long") ? 70 : 28;
		final int expectedFamilies=integer(options, "expectedfamilies");
		final long expectedMembers=number(options, "expectedmembers"), seed=number(options, "seed");
		final long maxMembers=number(options, "maxbucketmembers"), maxBytes=number(options, "maxbucketbytes");
		require(threads>0 && threads<=64, "Counting requires 1..64 workers with independently owned scorer scratch");
		final MemberLengthRange range=MemberLengthRange.from(options);
		require(integer(options, "minmemberaa")== (view.equals("long") ? 100 : 40)
			&& integer(options, "maxmemberaa")== (view.equals("long") ? 100 : 99),
			"Current accepted views are LONG100 and SHORT40..99; do not silently change the population");
		final Path matrixPath=input(options, "matrixfile"), output=Paths.get(required(options, "out"));
		require(output.isAbsolute() && output.equals(output.normalize()) && !Files.exists(output)
			&& !Files.isSymbolicLink(output) && Files.isDirectory(output.getParent()), "Require a fresh canonical output directory");
		final String matrixPin=required(options, "matrixfilesha80");
		require(matrixPin.equals(ProteinSubstitutionMatrix.CUSTOM40_FILE_SHA80)
			&& matrixPin.equals(DigestSuffix.file(matrixPath.toString())), "Pair selection must use the accepted POSSUM40 bytes");
		final ProteinSubstitutionMatrix matrix=ProteinSubstitutionMatrix.load("custom40", matrixPath.toString());
		require(matrix.scale==4 && matrix.linearGap==16 && matrix.sha80.equals(ProteinSubstitutionMatrix.CUSTOM40_BINDING_SHA80),
			"The study fixes POSSUM40 scores, scale4 and gap16 before observing any pairs");

		Shared.AMINO_IN=true;
		Shared.setThreads(threads);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final Path bgzip=input(options, "bgzip");
		final String bgzipPin=required(options, "bgzipsha80");
		configureBgzf(bgzip, bgzipPin);
		final FamilyCalibration.MemberInputs inputs=new FamilyCalibration.MemberInputs(input(options, "families"),
			required(options, "familiessha80"), input(options, "chunks"), required(options, "chunkssha80"),
			input(options, "consensus"), required(options, "consensussha80"), expectedFamilies, expectedMembers,
			integer(options, "buckets"), integer(options, "batches"), integer(options, "maxreferenceaa"), range,
			FamilyCalibration.VerificationMode.STRUCTURAL);
		Files.createDirectory(output);
		final WorkerState[] scratch=new WorkerState[threads];
		for(int i=0; i<threads; i++){scratch[i]=new WorkerState(matrix, span, seed, view.equals("long") ? 0 : 1);}
		final Counts total=new Counts();
		final Counts replayed=new Counts();
		final ByteBuilder streams=new ByteBuilder("bucket\tstream\n");
		final ByteBuilder work=new ByteBuilder("bucket\tfamilies\tmembers\tpossible_pairs\tselected_pairs\tcap_hit_families\tempty_families\tsingleton_families\treturned\tfailed\tnominal_dp_area\tdecoded_member_record_bytes\tdecoded_stream_bytes\tcompressed_stream_bytes\n");
		final ByteBuilder timings=new ByteBuilder("phase\tbucket\twall_nanoseconds\tworker_cpu_nanoseconds\n");
		final long setupEnd=System.nanoTime();
		long selectedFamilies=0, selectedMembers=0, computeCpu=0;
		if(pilot){timing(timings, "setup", -1, setupEnd-started, -1);}
		for(int bucket:buckets){
			final long loadStart=System.nanoTime(), previousBytes=inputs.decodedMemberRecordBytes();
			final FamilyCalibration.MemberInput[] families=inputs.loadBucket(bucket, maxMembers, maxBytes);
			final long loadEnd=System.nanoTime();
			long members=0, possible=0, selected=0, capped=0, empty=0, singletons=0;
			for(FamilyCalibration.MemberInput family:families){
				final int n=family.members.length;
				final long pairs=pairCount(n);
				members=Math.addExact(members, n); possible=Math.addExact(possible, pairs);
				selected=Math.addExact(selected, Math.min(CAP, pairs));
				if(pairs>CAP){capped++;}
				if(n==0){empty++;}else if(n==1){singletons++;}
			}
			selectedFamilies=Math.addExact(selectedFamilies, families.length);
			selectedMembers=Math.addExact(selectedMembers, members);
			final ByteBuilder[] retained=new ByteBuilder[families.length];
			final AtomicInteger next=new AtomicInteger();
			final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>(threads);
			for(WorkerState state:scratch){workers.add(new ProcessThread(families, retained, next, state, pilot));}
			final long computeStart=System.nanoTime();
			ThreadWaiter.startAndWait(workers);
			final long computeEnd=System.nanoTime();
			long bucketCpu=0;
			for(ProcessThread worker:workers){
				if(worker.failure!=null){throw new IllegalStateException("Member-pair bucket failed; partial outputs are invalid", worker.failure);}
				require(worker.success, "Joined pair worker must explicitly report success");
				bucketCpu=addCpu(bucketCpu, worker.cpuNanos);
			}
			computeCpu=addCpu(computeCpu, bucketCpu);
			final String name=String.format(Locale.ROOT, "pairs_%03d.tsv.gz", bucket);
			final long writeStart=System.nanoTime();
			long decodedBytes=PAIR_HEADER.length()+1L;
			final ByteStreamWriter writer=writer(output.resolve(name));
			try{
				writer.println(PAIR_HEADER);
				for(int i=0; i<retained.length; i++){
					require(retained[i]!=null, "Each canonical family must finish before ordered bucket publication");
					decodedBytes=Math.addExact(decodedBytes, retained[i].length());
					writer.print(retained[i]);
					retained[i]=null;
				}
			}finally{require(!writer.poisonAndWait(), "Compressed pair-stream writer failed");}
			final long writeEnd=System.nanoTime();
			final Counts bucketCounts=readStream(output.resolve(name), matrix, span);
			final long replayEnd=System.nanoTime();
			require(bucketCounts.selected==selected, "Retained bucket must contain the selected pairs from its original families");
			replayed.add(bucketCounts);
			streams.append(bucket).tab().append(name).nl();
			if(pilot){
				work.append(bucket).tab().append(families.length).tab().append(members).tab().append(possible).tab()
					.append(selected).tab().append(capped).tab().append(empty).tab().append(singletons).tab()
					.append(selected-bucketCounts.failed).tab().append(bucketCounts.failed).tab().append(bucketCounts.nominalDpArea).tab()
					.append(inputs.decodedMemberRecordBytes()-previousBytes).tab().append(decodedBytes).tab()
					.append(Files.size(output.resolve(name))).nl();
				timing(timings, "load", bucket, loadEnd-loadStart, -1);
				timing(timings, "compute", bucket, computeEnd-computeStart, bucketCpu);
				timing(timings, "write_compress", bucket, writeEnd-writeStart, -1);
				timing(timings, "replay", bucket, replayEnd-writeEnd, -1);
			}
		}
		final long verifyStart=System.nanoTime();
		if(pilot){inputs.verifySelectedBuckets(buckets);}else{inputs.verifyComplete();}
		for(WorkerState state:scratch){total.add(state.counts);}
		total.verify();
		require(total.selected==replayed.selected && total.failed==replayed.failed && total.nominalDpArea==replayed.nominalDpArea
			&& Arrays.deepEquals(total.ordered, replayed.ordered)
			&& Arrays.deepEquals(total.alignments, replayed.alignments)
			&& Arrays.deepEquals(total.columns, replayed.columns)
			&& Arrays.deepEquals(total.eligible, replayed.eligible),
			"Retained compressed streams must reproduce every gated/ungated count before publication");
		require(total.families==selectedFamilies && total.members==selectedMembers,
			"Worker totals must conserve every member and family in the loaded buckets");
		require(pilot || total.families==expectedFamilies && total.members==expectedMembers,
			"Every accepted family and assigned member must be accounted for, including singleton families");
		require(matrixPin.equals(DigestSuffix.file(matrixPath.toString())), "POSSUM40 input changed during the counting pass");
		require(bgzipPin.equals(DigestSuffix.file(bgzip.toString())), "Bound bgzip executable changed during the counting pass");
		final ByteBuilder recipe=new ByteBuilder("key\tvalue\nview\t").append(view)
			.append("\nseed\t").append(seed).append("\npair_cap\t100\norientation\tlower_member_ordinal_is_query\n")
			.append("sampling\tuniform_distinct_unordered_pairs;Floyd;SplittableRandom;family_rank_view_mix64\n")
			.append("weighting\traw_eligible_columns_after_family_cap;no_inverse_sampling_weights\n")
			.append("alignment\tPOSSUM40;gap16;pin_right_trim_left;positive_returned_traces\n")
			.append("span\t").append(span).append("\nprimary_population\tspan_gated\n")
			.append("matrix_binding\t").append(matrix.binding()).append("\ninput_byte_identity\tinherited_structural\n")
			.append("compression\texternal_bgzip_BGZF_level9\ncompression_threads\t1\n")
			.append("bgzip\t").append(bgzip.toString()).append("\nbgzip_sha80\t").append(bgzipPin)
			.append("\nbgzip_version\t").append(Data.BGZIP_VERSION)
			.append("\ncompression_command\tbgzip -c -f -l 9 -@ 1\n")
			.append("bias\taccepted_family_selection;POSSUM40_alignment_selection;identity_band_selection;pair_cap\n")
			.append("dimer_counts\tnot_computed\n");
		if(pilot){
			recipe.append("mode\ttiming_pilot_only\npilot_buckets\t").append(required(options, "pilotbuckets"))
				.append("\noriginal_buckets\t").append(inputs.buckets).append("\noriginal_families\t").append(expectedFamilies)
				.append("\noriginal_members\t").append(expectedMembers).append("\nthreads\t").append(threads)
				.append("\nnominal_dp_area\tsum_selected_query_length_times_reference_length;not_visited_cells\n")
				.append("decoded_member_bytes\tdata_records_excluding_headers_and_line_terminators\n")
				.append("compute_cpu_baseline\t").append(threads!=1 ? "not_single_worker" : computeCpu<0 ? "unavailable"
					: computeCpu>=10000000000L ? "at_least_10_seconds" : "under_10_seconds")
				.append("\nscientific_derivation\tforbidden;full_population_outputs_required\n");
			timing(timings, "verification", -1, System.nanoTime()-verifyStart, -1);
			write(output.resolve("pilot.recipe.tsv"), recipe);
			write(output.resolve("pilot.work.tsv"), work);
			write(output.resolve("pilot.streams.tsv"), streams);
			// This total ends before timing-table/marker publication; GNU time covers the whole invocation.
			timing(timings, "total_before_timing_publication", -1, System.nanoTime()-started, computeCpu);
			write(output.resolve("pilot.timings.tsv"), timings);
			write(output.resolve("PILOT_PASS"), new ByteBuilder("FAMILY_MEMBER_PAIR_TIMING_PILOT_PASS science=false\n"));
			return;
		}
		writeCounts(output.resolve("ungated.ordered_counts.tsv"), total, 0);
		writeCounts(output.resolve("gated.ordered_counts.tsv"), total, 1);
		write(output.resolve("streams.tsv"), streams);
		write(output.resolve("summary.tsv"), total.summary(view));
		write(output.resolve("recipe.tsv"), recipe);
		write(output.resolve("PASS"), new ByteBuilder("FAMILY_MEMBER_PAIR_SUBSTITUTIONS_PASS\n"));
	}

	/*--------------------------------------------------------------*/
	/*----------------       Pair Selection         ----------------*/
	/*--------------------------------------------------------------*/

	/** Exact number of unordered non-self record pairs; identical sequences remain distinct records. */
	public static long pairCount(int members){
		require(members>=0, "A family cannot have a negative member count");
		return (long)members*(members-1)/2;
	}

	/** Force the supported external level9 branch rather than native BGZF's level5 clamp.
	 * ReadWrite invokes the PATH-selected executable; verify that it is the explicitly pinned file.
	 */
	public static void configureBgzf(Path binary, String pin) throws java.io.IOException{
		require(binary!=null && pin!=null && pin.matches("[0-9a-f]{20}"), "Compression requires an executable and recorded sha80");
		require(binary.isAbsolute() && Files.isRegularFile(binary) && Files.isExecutable(binary)
			&& binary.equals(binary.toRealPath()) && binary.equals(bgzipExecutable()),
			"PATH-selected bgzip must resolve to the declared canonical executable");
		require(pin.equals(DigestSuffix.file(binary.toString())), "Compression executable differs from its bound sha80");
		require(Data.BGZIP() && Data.BGZIP_VERSION_levelFlag && Data.BGZIP_VERSION_threadsFlag,
			"External bgzip must advertise both -l and -@; silent level/backend fallback is forbidden");
		ReadWrite.ZIPLEVEL=9;
		ReadWrite.ALLOW_ZIPLEVEL_CHANGE=false;
		ReadWrite.USE_BGZF=ReadWrite.USE_BGZIP=ReadWrite.FORCE_BGZIP=true;
		ReadWrite.ALLOW_NATIVE_BGZF=false;
		ReadWrite.setZipThreads(1);
	}

	/** Resolve the first executable that the existing PATH-based compression launcher will select. */
	public static Path bgzipExecutable() throws java.io.IOException{
		final String path=System.getenv("PATH");
		require(path!=null, "External compression requires an explicit executable search path");
		for(String entry:path.split(java.io.File.pathSeparator, -1)){
			final Path candidate=Paths.get(entry.isEmpty() ? "." : entry).resolve("bgzip");
			if(Files.isRegularFile(candidate) && Files.isExecutable(candidate)){return candidate.toRealPath();}
		}
		throw new IllegalArgumentException("No bgzip executable exists on PATH; native compression is not a level9 substitute");
	}

	/** Uniform bounded sampling without enumerating the quadratic pair population.
	 * Floyd's collision case inserts j, which cannot already occur in the preceding [0,j) set.
	 */
	public static long[] sample(int members, long seed){
		final long possible=pairCount(members);
		final int count=(int)Math.min(CAP, possible);
		final long[] selected=new long[count];
		if(possible<=CAP){
			for(int i=0; i<count; i++){selected[i]=i;}
		}else{
			final LongHashSet used=new LongHashSet(256);
			final SplittableRandom random=new SplittableRandom(seed);
			int at=0;
			for(long j=possible-count; j<possible; j++){
				final long candidate=random.nextLong(j+1);
				final long chosen=used.contains(candidate) ? j : candidate;
				require(used.add(chosen), "Floyd sampling must insert one new pair on every step");
				selected[at++]=chosen;
			}
			Arrays.sort(selected);
		}
		assert(selected.length==count) : "The per-family cap controls both work and retained-output memory";
		return selected;
	}

	/** Decodes lexicographic upper-triangle rank into (lower ordinal, higher ordinal). */
	public static void decodePair(int members, long rank, int[] pair){
		require(pair!=null && pair.length>=2 && rank>=0 && rank<pairCount(members), "Pair rank must lie in the family's non-self triangle");
		int low=0, high=members-1;
		while(low+1<high){
			final int middle=low+(high-low)/2;
			if(rowStart(members, middle)<=rank){low=middle;}else{high=middle;}
		}
		pair[0]=low;
		pair[1]=(int)(low+1+rank-rowStart(members, low));
		assert(pair[0]<pair[1] && pair[1]<members) : "Triangle inversion must preserve distinct ordered member coordinates";
	}

	/** Number of pairs preceding the row with lower member ordinal i. */
	private static long rowStart(int members, int i){return (long)i*(2L*members-i-1)/2;}

	/** Fixed SplitMix64 finalizer, used only to separate recorded per-view/per-family RNG streams. */
	public static long familySeed(long seed, int view, int rank){
		require((view==0 || view==1) && rank>=0, "Sampling identity requires an explicit view and canonical family rank");
		long z=seed+0x9e3779b97f4a7c15L*(2L*rank+view+1);
		z=(z^(z>>>30))*0xbf58476d1ce4e5b9L;
		z=(z^(z>>>27))*0x94d049bb133111ebL;
		return z^(z>>>31);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Trace Accounting       ----------------*/
	/*--------------------------------------------------------------*/

	/** Counts and serializes one selected pair. Both populations use the identical recorded path. */
	public static void countTrace(byte[] query, byte[] reference, HbmPositionModel.TailResult trace,
			ProteinSubstitutionMatrix matrix, int span, Counts counts, ByteBuilder row, TraceScratch scratch){
		require(query!=null && reference!=null && matrix!=null && counts!=null && row!=null && scratch!=null && span>0,
			"Trace accounting requires explicit sequences, matrix, span and independently owned outputs");
		counts.attempted++;
		counts.nominalDpArea=Math.addExact(counts.nominalDpArea, (long)query.length*reference.length);
		if(trace==null){
			counts.failed++;
			row.append("0\t-1\t-1\t-1\t-1\t0\t0\tfalse\t.\t.\n");
			return;
		}
		require(trace.path!=null && trace.path.length>0 && trace.score>0 && trace.queryStart>=0 && trace.start>=0
			&& trace.queryEnd==query.length-1 && trace.end==reference.length-1,
			"Positive right-pinned alignment must retain a complete path to both sequence ends");
		final boolean gated=trace.queryEnd-trace.queryStart+1>=span && trace.end-trace.start+1>=span;
		int q=trace.queryStart, r=trace.start, identities=0, score=0, eligible=0;
		final ByteBuilder queryText=scratch.query, referenceText=scratch.reference;
		queryText.clear(); referenceText.clear();
		final int[] cells=scratch.cells;
		Arrays.fill(cells, 0);
		for(byte op:trace.path){
			if(op=='m'){
				require(q<query.length && r<reference.length, "Paired trace column must stay inside both sequence coordinates");
				final byte a=query[q++], b=reference[r++];
				queryText.append(letter(a)); referenceText.append(letter(b));
				score+=matrix.score(a, b);
				if(Blosum62.isStandard(a) && Blosum62.isStandard(b)){
					cells[20*a+b]++; eligible++;
					if(a==b){identities++;}
				}
			}else if(op=='I'){
				require(q<query.length, "Insertion consumes one remaining query residue");
				queryText.append(letter(query[q++])); referenceText.append('-'); score-=matrix.linearGap;
			}else if(op=='D'){
				require(r<reference.length, "Deletion consumes one remaining reference residue");
				queryText.append('-'); referenceText.append(letter(reference[r++])); score-=matrix.linearGap;
			}else{throw new IllegalArgumentException("Unknown alignment operation: "+(char)op);}
		}
		require(q==query.length && r==reference.length && (long)score*HbmPositionModel.SCALE==trace.score,
			"Retained columns must reconstruct exact endpoints and the selected matrix's DP score");
		final int band=band(identities, trace.path.length);
		for(int gate=0; gate< (gated ? 2 : 1); gate++){
			counts.alignments[gate][band]++;
			counts.columns[gate][band]+=trace.path.length;
			counts.eligible[gate][band]+=eligible;
			for(int cell=0; cell<400; cell++){counts.ordered[gate][band][cell]+=cells[cell];}
		}
		row.append(score).tab().append(trace.queryStart).tab().append(trace.queryEnd).tab()
			.append(trace.start).tab().append(trace.end).tab().append(identities).tab().append(trace.path.length).tab()
			.append(gated ? "true" : "false").tab().append(queryText.array, queryText.length()).tab()
			.append(referenceText.array, referenceText.length()).nl();
	}

	/** Integer comparisons implement strict40/60 boundaries without floating-point edge ambiguity. */
	public static int band(int identities, int columns){
		require(columns>0 && identities>=0 && identities<=columns, "Identity counts must fit a nonempty whole alignment");
		return 100L*identities<40L*columns ? 0 : 100L*identities<60L*columns ? 1 : 2;
	}

	/** Retains standard residues and the canonical X ambiguity code; gaps are separate trace operations. */
	private static char letter(byte code){
		require(code>=0 && (code<20 || code==Blosum62.X_CODE), "Retained trace requires canonical or X encoded residues");
		return code<20 ? (char)AminoAcid.numberToAcid[code] : 'X';
	}

	/** Reconstructs alignment observations using only the final stream, never original member bytes.
	 * Prefix residues are inert X placeholders; the recorded endpoints exclude them from scoring.
	 */
	public static Counts readStream(Path path, ProteinSubstitutionMatrix matrix, int span){
		return readStream(path, matrix, span, null);
	}

	/** Synchronous access to a fully replay-validated positive trace.
	 * Arrays are reused by the next record: consumers must finish before returning and must not mutate or retain them.
	 * Null attempts remain in the returned census and are never delivered as positive observations.
	 */
	public interface TraceConsumer{
		void accept(byte[] query, byte[] reference, HbmPositionModel.TailResult trace, int band, boolean spanPass);
	}

	/** Replay each row once, then optionally expose its validated path to an adjacent-column consumer. */
	public static Counts readStream(Path path, ProteinSubstitutionMatrix matrix, int span, TraceConsumer consumer){
		require(path!=null && matrix!=null && span>0, "Stream replay requires an explicit matrix and span policy");
		final Counts counts=new Counts();
		final TraceScratch scratch=new TraceScratch();
		final LineParser1 parser=new LineParser1('\t');
		final ByteBuilder expected=new ByteBuilder();
		final byte[][] queries=new byte[101][], references=new byte[101][], paths=new byte[201][];
		final ByteFile reader=ByteFile.makeByteFile1(path.toString(), false);
		int previousFamily=-1;
		long previousPair=-1;
		try{
			final byte[] header=reader.nextLine();
			require(header!=null && Arrays.equals(header, PAIR_HEADER.getBytes(java.nio.charset.StandardCharsets.US_ASCII)),
				"Retained pair stream must have the exact versioned column schema");
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				parser.set(line);
				require(parser.terms()==16, "Pair stream requires16 fields, including explicit null-alignment placeholders");
				final int family=parser.parseInt(0), first=parser.parseInt(1), second=parser.parseInt(2);
				final long rank=parser.parseLong(3);
				final int qLength=parser.parseInt(4), rLength=parser.parseInt(5);
				require(family>=0 && first>=0 && second>first && rank>=0 && family>=previousFamily
					&& (family!=previousFamily || rank>previousPair), "Retained pair coordinates must be distinct and canonically ordered");
				require(qLength>0 && qLength<=100 && rLength>0 && rLength<=100, "Retained member lengths must fit current100-aa views");
				previousFamily=family; previousPair=rank;
				if(queries[qLength]==null){queries[qLength]=new byte[qLength];}
				if(references[rLength]==null){references[rLength]=new byte[rLength];}
				final byte[] query=queries[qLength], reference=references[rLength];
				Arrays.fill(query, Blosum62.X_CODE); Arrays.fill(reference, Blosum62.X_CODE);
				final int columns=parser.parseInt(12), score=parser.parseInt(6);
				require(columns>=0 && columns<=qLength+rLength, "A pair trace cannot consume more residues than its two inputs");
				HbmPositionModel.TailResult trace=null;
				if(columns>0){
					require(score>0 && parser.length(14)==columns && parser.length(15)==columns,
						"Both retained aligned strings must have one symbol per positive-scoring path column");
					if(paths[columns]==null){paths[columns]=new byte[columns];}
					final byte[] operations=paths[columns];
					int q=parser.parseInt(7), r=parser.parseInt(9);
					require(q>=0 && q<qLength && r>=0 && r<rLength, "Retained alignment starts must be inside both sequences");
					for(int i=0; i<columns; i++){
						final byte a=parser.parseByte(14, i), b=parser.parseByte(15, i);
						require(a!='-' || b!='-', "A trace column cannot contain two gaps");
						operations[i]=a=='-' ? (byte)'D' : b=='-' ? (byte)'I' : (byte)'m';
						if(a!='-'){
							require(q<qLength, "Retained query columns must fit their declared span");
							query[q++]=decodeLetter(a);
						}
						if(b!='-'){
							require(r<rLength, "Retained reference columns must fit their declared span");
							reference[r++]=decodeLetter(b);
						}
					}
					trace=new HbmPositionModel.TailResult(Math.multiplyExact(score, HbmPositionModel.SCALE),
						parser.parseInt(7), parser.parseInt(8), parser.parseInt(9), parser.parseInt(10), operations);
				}
				expected.clear().append(family).tab().append(first).tab().append(second).tab().append(rank).tab()
					.append(qLength).tab().append(rLength).tab();
				countTrace(query, reference, trace, matrix, span, counts, expected, scratch);
				counts.selected++;
				require(expected.length()==line.length+1, "Decoded pair must reproduce the full retained row, including state fields");
				for(int i=0; i<line.length; i++){
					require(expected.array[i]==line[i], "Retained score, endpoints, identity, span state or residues disagree with replay");
				}
				if(consumer!=null && trace!=null){
					consumer.accept(query, reference, trace, band(parser.parseInt(11), columns),
						trace.queryEnd-trace.queryStart+1>=span && trace.end-trace.start+1>=span);
				}
			}
		}finally{require(!reader.close(), "Compressed pair stream failed during read-back");}
		counts.verify();
		return counts;
	}

	/** Stream output contains only uppercase standard letters and X, never alternate ambiguity aliases. */
	private static byte decodeLetter(byte letter){
		if(letter=='X'){return Blosum62.X_CODE;}
		require(letter>='A' && letter<='Z', "Retained trace residue must be an uppercase amino-acid letter");
		final byte code=AminoAcid.acidToNumber[letter];
		require(code>=0 && code<20 && AminoAcid.numberToAcid[code]==letter, "Only standard residues or X are valid retained symbols");
		return code;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Workers             ----------------*/
	/*--------------------------------------------------------------*/

	/** Per-worker trace buffers; no new400-cell array or gapped-text builder per pair. */
	public static final class TraceScratch{
		final ByteBuilder query=new ByteBuilder(), reference=new ByteBuilder();
		final int[] cells=new int[400];
	}

	/** Reused across all buckets; query/reference encodings must have separate mutable storage. */
	private static final class WorkerState{
		WorkerState(ProteinSubstitutionMatrix matrix_, int span_, long seed_, int view_){
			matrix=matrix_; span=span_; seed=seed_; view=view_;
			scorer=new HbmPositionModel.ConsensusScorer(matrix);
		}
		/** Emits at most100 pairs for one family in increasing pair-rank order. */
		ByteBuilder process(FamilyCalibration.MemberInput family){
			assert(family!=null) : "Bucket ownership transfers every selected family to exactly one worker";
			final int n=family.members.length;
			final long possible=pairCount(n);
			counts.families++; counts.members=Math.addExact(counts.members, n);
			counts.possible=Math.addExact(counts.possible, possible);
			if(n==0){counts.empty++;}else if(n==1){counts.singletons++;}
			if(possible>CAP){counts.capped++;}
			final long[] selected=sample(n, familySeed(seed, view, family.rank));
			counts.selected+=selected.length;
			final ByteBuilder retained=new ByteBuilder();
			for(long rank:selected){
				decodePair(n, rank, pair);
				final byte[] query=queryBuffer.encode(family.members[pair[0]]).enc;
				final byte[] reference=referenceBuffer.encode(family.members[pair[1]]).enc;
				final HbmPositionModel.TailResult trace=scorer.resetConsensus(reference)
					.alignTail(query, Math.max(query.length, reference.length), true, true);
				retained.append(family.rank).tab().append(pair[0]).tab().append(pair[1]).tab().append(rank).tab()
					.append(query.length).tab().append(reference.length).tab();
				countTrace(query, reference, trace, matrix, span, counts, retained, traceScratch);
			}
			return retained;
		}
		final Counts counts=new Counts();
		final FamilyCalibration.Buffers queryBuffer=new FamilyCalibration.Buffers(), referenceBuffer=new FamilyCalibration.Buffers();
		final HbmPositionModel.ConsensusScorer scorer;
		final ProteinSubstitutionMatrix matrix;
		final TraceScratch traceScratch=new TraceScratch();
		final int[] pair=new int[2];
		final int span, view;
		final long seed;
	}

	/** Bounded bucket workers share only family indices; failures are rethrown after every join. */
	private static final class ProcessThread extends Thread{
		ProcessThread(FamilyCalibration.MemberInput[] families_, ByteBuilder[] retained_, AtomicInteger next_, WorkerState state_, boolean timed_){
			families=families_; retained=retained_; next=next_; state=state_; timed=timed_;
		}
		@Override public void run(){
			final ThreadMXBean bean=timed ? cpuBean() : null;
			final long cpuStart=FamilyCalibration.cpuTime(bean);
			try{
				for(int i=next.getAndIncrement(); i<families.length; i=next.getAndIncrement()){
					retained[i]=state.process(families[i]); families[i]=null;
				}
				success=true;
			}catch(Throwable caught){failure=caught; next.set(families.length);}
			finally{cpuNanos=FamilyCalibration.cpuElapsed(cpuStart, FamilyCalibration.cpuTime(bean));}
		}
		final FamilyCalibration.MemberInput[] families;
		final ByteBuilder[] retained;
		final AtomicInteger next;
		final WorkerState state;
		final boolean timed;
		long cpuNanos=-1;
		Throwable failure;
		boolean success;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Counts and Output      ----------------*/
	/*--------------------------------------------------------------*/

	/** Independent descriptive populations: ungated0 and span-gated1, with identical pair attempts. */
	public static final class Counts{
		/** Merge worker-local counters once, preserving overflow detection and population conservation. */
		public void add(Counts other){
			require(other!=null && other!=this, "Reduction requires a distinct worker's observations");
			families=Math.addExact(families, other.families); members=Math.addExact(members, other.members);
			possible=Math.addExact(possible, other.possible); selected=Math.addExact(selected, other.selected);
			attempted=Math.addExact(attempted, other.attempted); failed=Math.addExact(failed, other.failed);
			empty=Math.addExact(empty, other.empty); singletons=Math.addExact(singletons, other.singletons);
			capped=Math.addExact(capped, other.capped);
			nominalDpArea=Math.addExact(nominalDpArea, other.nominalDpArea);
			for(int g=0; g<2; g++){for(int b=0; b<3; b++){
				alignments[g][b]=Math.addExact(alignments[g][b], other.alignments[g][b]);
				columns[g][b]=Math.addExact(columns[g][b], other.columns[g][b]);
				eligible[g][b]=Math.addExact(eligible[g][b], other.eligible[g][b]);
				for(int c=0; c<400; c++){ordered[g][b][c]=Math.addExact(ordered[g][b][c], other.ordered[g][b][c]);}
			}}
		}
		/** Every attempt is either a failed alignment or exactly one ungated identity-band observation. */
		public void verify(){
			long aligned=0;
			for(int b=0; b<3; b++){
				aligned=Math.addExact(aligned, alignments[0][b]);
				for(int g=0; g<2; g++){
					long sum=0;
					for(long value:ordered[g][b]){require(value>=0, "Directed cell counts cannot be negative");sum=Math.addExact(sum, value);}
					require(sum==eligible[g][b] && sum<=columns[g][b], "Only standard paired columns enter the400-cell table");
				}
				require(alignments[1][b]<=alignments[0][b], "Span-gated observations are a subset of the same ungated traces");
			}
			require(attempted==selected && attempted==Math.addExact(failed, aligned), "Selected pair attempts must conserve explicit success/failure states");
		}
		/** One small population summary keeps failed attempts separate from identity bands. */
		public ByteBuilder summary(String view){
			verify();
			final ByteBuilder out=new ByteBuilder("view\tpopulation\tband\talignments\tcolumns\teligible_columns\n");
			for(int g=0; g<2; g++){for(int b=0; b<3; b++){
				out.append(view).tab().append(g==0 ? "ungated" : "span_gated").tab().append(BANDS[b]).tab()
					.append(alignments[g][b]).tab().append(columns[g][b]).tab().append(eligible[g][b]).nl();
			}}
			out.append("#families\t").append(families).append("\n#members\t").append(members)
				.append("\n#possible_pairs\t").append(possible).append("\n#selected_pairs\t").append(selected)
				.append("\n#failed_alignments\t").append(failed).append("\n#cap_hit_families\t").append(capped)
				.append("\n#empty_families\t").append(empty).append("\n#singleton_families\t").append(singletons).nl();
			return out;
		}
		public long families, members, possible, selected, attempted, failed, capped, empty, singletons, nominalDpArea;
		public final long[][] alignments=new long[2][3], columns=new long[2][3], eligible=new long[2][3];
		public final long[][][] ordered=new long[2][3][400];
	}

	/** Preserve the historical eight-field schema; the consensus endpoint is the second member here. */
	private static void writeCounts(Path path, Counts counts, int gate){
		final ByteStreamWriter writer=writer(path);
		final ByteBuilder row=new ByteBuilder();
		try{
			writer.println("band\tquery_aa\tconsensus_aa\tcount\tquery_monomer\tconsensus_monomer\talignments\taligned_pairs");
			for(int b=0; b<3; b++){
				final long[] query=new long[20], reference=new long[20];
				for(int a=0; a<20; a++){for(int c=0; c<20; c++){query[a]+=counts.ordered[gate][b][20*a+c];reference[c]+=counts.ordered[gate][b][20*a+c];}}
				for(int a=0; a<20; a++){for(int c=0; c<20; c++){
					row.clear().append(BANDS[b]).tab().append(letter((byte)a)).tab().append(letter((byte)c)).tab()
						.append(counts.ordered[gate][b][20*a+c]).tab().append(query[a]).tab().append(reference[c]).tab()
						.append(counts.alignments[gate][b]).tab().append(counts.eligible[gate][b]);
					writer.println(row);
				}}
			}
		}finally{require(!writer.poisonAndWait(), "Directed-count writer failed");}
	}

	/** Reject unknown/duplicate flags before touching inputs or outputs. */
	private static HashMap<String, String> parse(String[] args){
		final HashMap<String, String> options=new HashMap<String, String>();
		for(String arg:args){
			final String[] fields=Parse.splitOnFirst(arg, '=');
			require(fields.length==2 && !fields[1].isEmpty(), "Require nonempty flag=value arguments");
			final String key=fields[0].toLowerCase(Locale.ROOT);
			require(OPTIONS.contains(key) && options.put(key, fields[1])==null, "Unknown or repeated argument: "+key);
		}
		return options;
	}
	/** Require explicit values rather than infer scientific settings from defaults. */
	private static String required(HashMap<String, String> options, String key){
		final String value=options.get(key); require(value!=null, "Missing option: "+key); return value;
	}
	/** Parse nonnegative sizes and the recorded master seed. */
	private static long number(HashMap<String, String> options, String key){
		final long value=Long.parseLong(required(options, key)); require(value>=0, "Negative option: "+key); return value;
	}
	private static int integer(HashMap<String, String> options, String key){return Math.toIntExact(number(options, key));}
	/** Input aliases and relative paths must not bypass the launch binding. */
	private static Path input(HashMap<String, String> options, String key){
		final Path path=Paths.get(required(options, key));
		require(path.isAbsolute() && path.equals(path.normalize()) && Files.isRegularFile(path) && !Files.isSymbolicLink(path),
			"Require an explicit regular normalized input: "+key);
		return path;
	}
	/** Fresh-output writer with BBTools compression and failure reporting. */
	private static ByteStreamWriter writer(Path path){
		require(!Files.exists(path) && !Files.isSymbolicLink(path), "Refusing to replace retained study evidence: "+path);
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(), false, false, true); writer.start(); return writer;
	}
	private static void write(Path path, ByteBuilder data){
		final ByteStreamWriter writer=writer(path);
		try{writer.print(data);}finally{require(!writer.poisonAndWait(), "Study output failed: "+path);}
	}
	private static void require(boolean valid, String message){if(!valid){throw new IllegalArgumentException(message);}}

	/** Preserve original grid IDs; a null selection means the existing complete counting pass. */
	public static int[] selectedBuckets(int buckets, String selection){
		require(buckets>0, "Bucket selection requires the original positive grid size");
		if(selection==null){
			final int[] all=new int[buckets];
			for(int i=0; i<buckets; i++){all[i]=i;}
			return all;
		}
		final String[] tokens=selection.split(",", -1);
		require(tokens.length<=16, "A bounded timing pilot supports at most16 declared original buckets");
		final int[] selected=new int[tokens.length];
		int previous=-1;
		for(int i=0; i<tokens.length; i++){
			final int value=Integer.parseInt(tokens[i]);
			require(value>previous && value<buckets && tokens[i].equals(Integer.toString(value)),
				"Pilot bucket IDs must be canonical ascending distinct integers in the original grid");
			selected[i]=value; previous=value;
		}
		return selected;
	}

	/** Worker CPU is optional; never enable a VM-wide counter or substitute zero when unavailable. */
	private static ThreadMXBean cpuBean(){
		try{return ManagementFactory.getThreadMXBean();}
		catch(SecurityException | UnsupportedOperationException unavailable){return null;}
	}
	/** Preserve unavailable status while reducing independently measured worker intervals. */
	private static long addCpu(long left, long right){return left<0 || right<0 ? -1 : Math.addExact(left, right);}
	/** Emit integer elapsed intervals with explicit unavailable CPU; no NaN timing sentinel. */
	private static void timing(ByteBuilder rows, String phase, int bucket, long wall, long cpu){
		require(wall>=0 && cpu>=-1, "Pilot timings require nonnegative elapsed intervals or explicit unavailable CPU");
		rows.append(phase).tab().append(bucket).tab().append(wall).tab().append(cpu<0 ? "NA" : Long.toString(cpu)).nl();
	}

	public static final int CAP=100;
	public static final String[] BANDS={"<40", "40-60", ">=60"};
	public static final String PAIR_HEADER="family_rank\tmember1\tmember2\tpair_rank\tlength1\tlength2\tscore\tquery_start\tquery_end\treference_start\treference_end\tidentities\tcolumns\tspan_pass\tquery_aligned\treference_aligned";
	private static final ArrayList<String> OPTIONS=new ArrayList<String>(Arrays.asList("view", "families", "familiessha80",
		"chunks", "chunkssha80", "consensus", "consensussha80", "expectedfamilies", "expectedmembers", "buckets", "batches",
		"maxreferenceaa", "maxbucketmembers", "maxbucketbytes", "minmemberaa", "maxmemberaa", "matrixfile", "matrixfilesha80",
		"bgzip", "bgzipsha80", "verification", "seed", "t", "out", "mode", "pilotbuckets"));
}
