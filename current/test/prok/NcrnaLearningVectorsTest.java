package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;

import consensus.BaseGraph;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.LongHashSet;
import ml.DataLoader;
import ml.Matrix;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaFamily;
import prok.NcrnaLearningPopulation;
import prok.NcrnaLearningTables;
import prok.NcrnaLearningVectors;
import prok.NcrnaSupportIO;
import prok.NcrnaTextIO.Input;
import prok.NcrnaWorkers;
import shared.Shared;
import stream.FastaStreamerST;
import stream.FastaStreamerZT;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Native matrix-reader, row/label conservation, shard and compressed-format fixture.
 * Leaves small outputs under out= for a separate native regression-reader check.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaLearningVectorsTest{

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, null, false);
		try(AutoCloseable log=()->Shared.closeStream(pp.outstream)){
			run(pp.args);
		}
	}

	/** Main retains log ownership across native matrix and worker-failure controls. */
	private static void run(final String[] args) throws Exception{
		final Parser parser=new Parser();
		for(String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown fixture parameter: "+arg);}
		}
		if(parser.out1==null){throw new IllegalArgumentException("Fixture requires a fresh out= directory");}
		final Path root=Paths.get(parser.out1);
		Files.createDirectory(root);
		final String roster=root.resolve("roster.tsv").toString(), loci=root.resolve("loci.tsv").toString();
		final String graph=root.resolve("graph.tsv").toString(), exclude=root.resolve("exclude.txt").toString();
		final String fasta=root.resolve("contexts.fa").toString();
		NcrnaLearningPopulationTest.write(roster, NcrnaLearningPopulationTest.roster());
		NcrnaLearningPopulationTest.write(loci, NcrnaLearningPopulationTest.loci());
		NcrnaLearningPopulationTest.write(graph, NcrnaLearningPopulationTest.graph());
		NcrnaLearningPopulationTest.write(exclude, "gene0\ngene2\n");
		NcrnaLearningPopulationTest.write(fasta, contexts());
		final NcrnaLearningPopulation population=NcrnaLearningPopulation.load(NcrnaSupportIO.read(graph),
			new String[]{"common", "rare"}, roster, loci, exclude, 42);
		final byte[] consensus=new byte[100];
		Arrays.fill(consensus, (byte)'A');
		final NcrnaFamily family=new NcrnaFamily("fixture", new byte[][]{consensus},
			new BaseGraph[]{new BaseGraph("exact", consensus, null, 0, 0)}, new String[]{"exact"},
			new LongHashSet(), 17, 16, 20, 7, 1, false, 0, 0, 0, 0);
		final int buffer=Shared.bufferLen();
		final boolean trim=Shared.TRIM_READ_DESCRIPTION, upper=Read.TO_UPPER_CASE;
		final boolean force=ReadWrite.FORCE_BGZIP, bgzf=ReadWrite.USE_BGZF;
		Shared.setBufferLen(2);
		try{
			// Deliberately non-default settings make restoration falsifiable.
			Shared.TRIM_READ_DESCRIPTION=true;
			Read.TO_UPPER_CASE=false;
			ReadWrite.FORCE_BGZIP=false;
			ReadWrite.USE_BGZF=false;
			final NcrnaLearningTables tuning=NcrnaLearningTables.build(population, fasta, 3, 0, 2, 1, false, 3);
			final NcrnaLearningVectors train=NcrnaLearningVectors.write(population, fasta, family, tuning, 1, 1, 3, 3,
				root.resolve("train").toString());
			final NcrnaLearningTables all=NcrnaLearningTables.build(population, fasta, 3, 0, 2, 1, true, 2);
			final NcrnaLearningVectors refit=NcrnaLearningVectors.write(population, fasta, family, all, 1, 1, 3, 2,
				root.resolve("refit").toString());
			assert(Shared.TRIM_READ_DESCRIPTION && !Read.TO_UPPER_CASE && !ReadWrite.FORCE_BGZIP && !ReadWrite.USE_BGZF)
				: "Table/vector passes must restore reader and compressor settings after real transitions";
			check(population, train, false);
			check(population, refit, true);
			assert(refit.validationPath==null && Arrays.equals(refit.shardBases, new long[]{600, 600, 600}))
				: "Refit uses all eighteen eligible loci in three equal-length training shards, with no holdout claim";
			final NcrnaLearningVectors wide=NcrnaLearningVectors.write(population, fasta, family, all, 40, 40, 3, 1,
				root.resolve("wide").toString());
			checkWide(population, wide);
			checkStartup(fasta);
			final String expected="loaded "+train.rows[0]+" training samples, "+train.rows[3]+" external validation samples, dims=[10, 1]\n";
			NcrnaLearningPopulationTest.write(root.resolve("regression.expected").toString(), expected);
			System.out.println("Learning vectors passed: 108 narrow and 2376 clipped wide rows, complete locus accounting, balanced whole-locus shards, validation/refit separation, block handoff, native classifier loading and BGZF format");
		}finally{
			Shared.setBufferLen(buffer);
			Shared.TRIM_READ_DESCRIPTION=trim;
			Read.TO_UPPER_CASE=upper;
			ReadWrite.FORCE_BGZIP=force;
			ReadWrite.USE_BGZF=bgzf;
		}
	}

	/** Real native input is exhausted after first, later, and post-start failures.
	 * The same usage helper owns cohorts in the production tools; no process-wide
	 * worker factory or shared ThreadWaiter behavior is replaced.
	 */
	private static void checkStartup(final String fasta) throws Exception{
		final int threads=Shared.threads();
		final boolean simd=Shared.SIMD, low=Shared.LOW_MEMORY;
		try{
			Shared.SIMD=false;
			Shared.LOW_MEMORY=false;
			Shared.setThreads(3);
			checkStartup(fasta, false);
			Shared.setThreads(4);
			checkStartup(fasta, true);
		}finally{
			Shared.setThreads(threads);
			Shared.SIMD=simd;
			Shared.LOW_MEMORY=low;
		}
	}

	/** Scalar synchronous and background parsers exercise different terminal paths. */
	private static void checkStartup(final String fasta, final boolean background) throws Exception{
		for(int mode=0; mode<(background ? 5 : 4); mode++){
			final FileFormat format=FileFormat.testInput(fasta, FileFormat.FASTA, null, true, true);
			final Streamer input=mode==4 ? new FastaStreamerST(format, 0, -1){
				@Override
				public ListNum<Read> nextList(){
					if(first){
						first=false;
						Thread.currentThread().interrupt();
					}
					return super.nextList();
				}
				private boolean first=true;
			} : StreamerFactory.makeStreamer(format, null, false, -1, false, true, 1);
			assert(background ? input instanceof FastaStreamerST : input instanceof FastaStreamerZT)
				: "Explicit scalar policy and thread budgets must exercise the declared native input backend";
			final AtomicInteger consumed=new AtomicInteger();
			final Runnable consume=()->{
				for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
					consumed.addAndGet(batch.size());
				}
			};
			final Error expected=new AssertionError("Injected consumer startup failure");
			final boolean afterStart=mode==2;
			final ArrayList<Thread> threads=new ArrayList<Thread>();
			if(mode==1){threads.add(new Thread(consume));}
			threads.add(new Thread(consume){
				@Override
				public synchronized void start(){
					if(afterStart){super.start();}
					throw expected;
				}
			});
			final Thread neverStarted=new Thread(consume);
			threads.add(neverStarted);
			// FastaStreamerST records the deliberately injected interrupt as an error.
			// Keep native close semantics; the assertions below test cause and terminal state.
			try(AutoCloseable cleanup=()->{
				Thread.interrupted();
				input.close();
			}){
				input.start();
				if(mode==3){Thread.currentThread().interrupt();}
				Throwable caught=null;
				try(NcrnaWorkers cohort=new NcrnaWorkers()){
					cohort.startAndJoin(threads, input);
				}catch(RuntimeException | Error failure){caught=failure;}
				assert(caught==expected) : "Caller cleanup must retain the original startup failure";
				final boolean interrupted=Thread.interrupted();
				assert(interrupted==(mode>=3)) : "Observed interruption must survive cleanup; mode="+mode;
				for(final Thread worker:threads){
					assert(!worker.isAlive()) : "No started producer may outlive the scope owning its writers";
				}
				assert(neverStarted.getState()==Thread.State.NEW) : "Startup stops at the first failure";
				assert(!input.hasMore() && input.nextList()==null && input.readsProcessed()==21)
					: "Cleanup must reach the actual terminal after consuming all records, not an interrupted null";
				assert(consumed.get()==(mode==1 || mode==2 ? 21 : 0))
					: "A started cohort consumes every original record before cleanup; mode="+mode;
			}
		}
		System.out.println("Vector startup cleanup passed: background="+background+" first/later/post-start failures and interruptions retain cause, terminal input and joined producers");
	}

	private static void check(final NcrnaLearningPopulation population, final NcrnaLearningVectors vectors, final boolean refit){
		assert(vectors.records==21 && vectors.filtered==3 && vectors.unavailable==0 && vectors.noModel==0
			&& vectors.processed==18 && vectors.candidates==108 && vectors.outcomes[0]==108)
			: "Eighteen eligible loci emit six candidates each; excluded and partial loci remain in original accounting";
		final int[][] observed=new int[21][6];
		long total=0, positive=0;
		for(int shard=0; shard<vectors.rows.length; shard++){
			final String path=shard<3 ? vectors.trainPaths[shard] : vectors.validationPath;
			assert(ReadWrite.isBGZip(path)) : "Vector files must be BGZF, not a mislabeled uncompressed or ordinary gzip stream";
			checkText(path);
			final Matrix matrix=DataLoader.load(path, Integer.MAX_VALUE, false, 0, 0, false, 0)[0].matrix;
			assert(matrix.numInputs==10 && matrix.numOutputs==1
				&& matrix.inputs.length==vectors.rows[shard]) : "The actual classifier reader must accept every written row with the native dimensions";
			long labels=0;
			for(int row=0; row<matrix.inputs.length; row++){
				final float[] values=matrix.inputs[row];
				final int gene=Math.round(values[9]*100)-1;
				assert(gene>=0 && gene<21 && population.eligible(gene)) : "GC tags independently identify every emitted original locus";
				final boolean stop=values[4]==1;
				assert(stop || values[4]==0) : "Endpoint direction is binary";
				final int length=Math.round(values[8]*100);
				final int offset=stop ? length-100 : 100-length;
				assert(offset>=-1 && offset<=1) : "Only the declared three offsets per end may be emitted";
				final int label=(int)matrix.outputs[row][0];
				assert(matrix.outputs[row][0]==label && label==(offset==0 ? 1 : 0)) : "Only the original endpoint receives label one";
				assert(shard==vectors.shardByGene[gene] && (refit || (shard==3)==population.split().validationByGene[gene]))
					: "Every row stays with its original locus and declared train/validation membership";
				observed[gene][(stop ? 3 : 0)+offset+1]++;
				labels+=label;
			}
			assert(labels==vectors.positives[shard]) : "Native parsed labels must equal emitted positive counts";
			total+=matrix.inputs.length;
			positive+=labels;
		}
		for(int gene=0; gene<21; gene++){
			for(int count:observed[gene]){assert(count==(population.eligible(gene) ? 1 : 0)) : "Exactly one of each end/offset row per eligible locus";}
		}
		assert(total==108 && positive==36) : "Eighteen loci contribute six rows and two true endpoints each";
		final long[] lengths=vectors.shardBases.clone();
		Arrays.sort(lengths);
		assert(lengths[2]-lengths[0]<=100) : "Equal-length loci must balance training shard lengths within one locus";
	}

	/** One worker writes more than a full block per shard, exercising buffer transfer and reuse. */
	private static void checkWide(final NcrnaLearningPopulation population, final NcrnaLearningVectors vectors){
		assert(vectors.candidates==2916 && Arrays.equals(vectors.outcomes, new long[]{2376, 360, 180, 0}))
			: "Eighteen loci attempt 162 offsets each; twenty start and ten stop offsets exceed their actual flanks";
		final int[][] observed=new int[21][162];
		for(int shard=0; shard<3; shard++){
			checkText(vectors.trainPaths[shard]);
			final Matrix matrix=DataLoader.load(vectors.trainPaths[shard], Integer.MAX_VALUE, false, 0, 0, false, 0)[0].matrix;
			assert(matrix.inputs.length==792 && vectors.rows[shard]==792 && vectors.positives[shard]==12)
				: "Each equal shard has six whole loci, 132 emitted offsets per locus and two true labels each";
			for(int row=0; row<matrix.inputs.length; row++){
				final float[] values=matrix.inputs[row];
				final int gene=Math.round(values[9]*100)-1;
				final boolean stop=values[4]==1;
				final int length=Math.round(values[8]*100);
				final int offset=stop ? length-100 : 100-length;
				assert(gene>=0 && gene<21 && population.eligible(gene) && vectors.shardByGene[gene]==shard)
					: "Transferred blocks cannot move rows between loci or shards";
				assert(offset>=-40 && offset<=40 && matrix.outputs[row][0]==(offset==0 ? 1 : 0))
					: "Wide rows retain their declared offsets and exact-boundary labels";
				observed[gene][(stop ? 81 : 0)+offset+40]++;
			}
		}
		for(int gene=0; gene<21; gene++){
			for(int end=0; end<2; end++){
				for(int offset=-40; offset<=40; offset++){
					final boolean fits=end==0 ? offset>=-20 : offset<=30;
					assert(observed[gene][end*81+offset+40]==(population.eligible(gene) && fits ? 1 : 0))
						: "Every valid wide candidate occurs once; block flushing cannot duplicate or lose a row";
				}
			}
		}
	}

	private static void checkText(final String path){
		try(Input owner=new Input(path)){
			final ByteFile input=owner.reader;
			assert(Arrays.equals(input.nextLine(), "#dims\t10\t1".getBytes(StandardCharsets.US_ASCII))) : "Native dimensions header must precede all worker blocks";
			for(byte[] row=input.nextLine(); row!=null; row=input.nextLine()){
				int tabs=0;
				for(byte b:row){
					if(b=='\t'){tabs++;}
					else{assert(b=='-' || b=='.' || b>='0' && b<='9') : "Output numbers must use fixed decimals, never exponent notation";}
				}
				assert(tabs==10) : "Each row has ten input fields followed by one label";
			}
		}
	}

	private static String contexts(){
		final ByteBuilder b=new ByteBuilder();
		for(int gene=0; gene<21; gene++){
			b.append(">gene").append(gene).append(" source=genome").append(gene==20 ? 0 : gene)
				.append(" lflank=20 rflank=30 contig_gc=").append((gene+1)/100.0, 6).nl();
			for(int base=0; base<150; base++){b.append(base<20 ? 'G' : base<120 ? 'A' : 'T');}
			b.nl();
		}
		return b.toString();
	}
}
