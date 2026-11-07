package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.SequenceInput;
import prok.TrnaBoundaryFeatures.BoundaryType;
import prok.TrnaBoundaryFeatures.NinemerTable;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Screens endpoint-word geometry using the native population and table APIs.
 * Adapted from the Streamer and ByteFile templates: native sequence batches,
 * standard Parser arguments and owned byte writers. Each geometry builds tables
 * on the internal training partition, then streams all original contexts again.
 * Training uses the existing key-level leave-one-out view; validation never teaches.
 * Only one geometry's tables and primitive counters are retained at a time.
 * This measures count discrimination, not caller accuracy or network predictions.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaEndpointScreen{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, NcrnaEndpointScreen.class, false);
		try{
			new NcrnaEndpointScreen(pp.args).process(pp.outstream);
			timer.stop();
			pp.outstream.println("Time: "+timer);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Requires the complete original population; this screen has no exclusion option. */
	public NcrnaEndpointScreen(final String[] args){
		final Parser parser=new Parser();
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("graph")){graph=b;}
			else if(a.equals("original")){original=b;}
			else if(a.equals("roster")){roster=b;}
			else if(a.equals("loci")){loci=b;}
			else if(a.equals("splitseed")){seed=NcrnaNumbers.longValue(a, b);}
			else if(a.equals("radius")){radius=NcrnaNumbers.intValue(a, b);}
			else if(a.equals("ks")){
				final String[] terms=b.split(",");
				ks=new int[terms.length];
				for(int i=0; i<ks.length; i++){ks[i]=NcrnaNumbers.intValue(a, terms[i]);}
			}else if(a.equals("ratios")){
				final String[] terms=b.split(",");
				ratios=new double[terms.length];
				for(int i=0; i<ratios.length; i++){ratios[i]=Double.parseDouble(terms[i]);}
			}else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		contexts=parser.in1;
		check(contexts!=null && parser.out1!=null && graph!=null && original!=null && roster!=null && loci!=null,
			"Require in=contexts out=fresh_directory graph=original_graph original=consensus roster= and loci=");
		check(seed>=0 && radius>=1 && radius<Integer.MAX_VALUE && ks.length>0,
			"Require nonnegative splitseed, positive radius below integer maximum, and nonempty ks");
		check(parser.in2==null && parser.out2==null && !parser.append && parser.maxReads<0 && parser.samplerate==1,
			"Screen the complete unpaired original population; sampling or append would change denominators");
		Arrays.sort(ks);
		for(int i=0; i<ks.length; i++){
			check(ks[i]>=1 && ks[i]<=15 && (i==0 || ks[i]!=ks[i-1]), "Distinct ks must fit the native one-to-fifteen-base tables");
		}
		validateRatios(ratios);
		ratioLabels=new String[ratios.length];
		for(int i=0; i<ratios.length; i++){ratioLabels[i]=Double.toString(ratios[i]);}
		check(Tools.testInputFiles(false, true, contexts, graph, original, roster, loci), "Unreadable screen input");
		root=new File(parser.out1);
		check(!root.exists() && root.getAbsoluteFile().getParentFile().isDirectory(), "Output must be fresh beneath an existing directory");
		workers=Shared.threads();
		check(workers>0, "Native table counting needs at least one worker");
		ReadWrite.setZipThreads(workers);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Population          ----------------*/
	/*--------------------------------------------------------------*/

	public void process(final PrintStream log){
		population=NcrnaLearningPopulation.load(NcrnaSupportIO.read(graph),
			NcrnaSupport.loadLibrary(original).names, roster, loci, null, seed);
		check(root.mkdir(), "Could not create "+root);
		writePopulation();
		try(Output summary=new Output(file("summary.tsv")); Output curves=new Output(file("curves.tsv"))){
			summary.writer.println(SUMMARY_HEADER);
			curves.writer.println("k\tinside\toutside\tboundary\tpartition\tradius\tratio\teligible_loci\tdecisive_true\tdecisive_other\tabstain\ttrue_fraction\tother_fraction");
			for(final int k:ks){
				for(int outside=0; outside<=k; outside++){
					final int inside=k-outside;
					final NcrnaLearningTables tables=NcrnaLearningTables.build(population, contexts, inside, outside, inside, outside, false, workers);
					final Counts[] counts={new Counts(tables.start, false, radius, ratios), new Counts(tables.stop, false, radius, ratios),
						new Counts(tables.start, true, radius, ratios), new Counts(tables.stop, true, radius, ratios)};
					final long[][] accounting=score(tables, counts);
					for(int index=0; index<counts.length; index++){write(counts[index], accounting[index/2], tables.learning, summary, curves);}
					log.println("Screened k="+k+", inside="+inside+", outside="+outside+", original_loci="+population.ids.length);
				}
			}
		}
		log.println("NCRNA_ENDPOINT_SCREEN_COMPLETE original_loci="+population.ids.length+" no_training=true");
	}

	/** Persists the exact split and eligibility once, before any table is learned. */
	private void writePopulation(){
		try(Output output=new Output(file("population.tsv"))){
			output.writer.println("locus\tsource\tpartition\teligible");
			final ByteBuilder row=new ByteBuilder();
			for(int gene=0; gene<population.ids.length; gene++){
				output.writer.print(row.clear().append(population.ids[gene]).tab().append(population.sourceName(gene)).tab()
					.append(population.split().validationByGene[gene] ? "VALIDATION" : "TRAIN_LOO").tab()
					.append(population.eligible(gene) ? "true" : "false").nl());
			}
		}
	}

	/** Serial scoring borrows batches; native table counting above may use multiple workers. */
	private long[][] score(final NcrnaLearningTables tables, final Counts[] counts){
		final long[][] accounting=new long[2][6];
		final boolean[] seen=new boolean[population.ids.length];
		final NcrnaLearningContext context=new NcrnaLearningContext();
		final boolean trim=Shared.TRIM_READ_DESCRIPTION, upper=Read.TO_UPPER_CASE;
		Shared.TRIM_READ_DESCRIPTION=false;
		Read.TO_UPPER_CASE=true;
		try(SequenceInput owner=new SequenceInput(StreamerFactory.makeStreamer(
				FileFormat.testInput(contexts, FileFormat.FASTA, null, true, true), null, false, -1, false, true, 1))){
			final Streamer input=owner.stream;
			try{
				input.start();
				for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
					for(final Read read:batch){
						context.set(read, population);
						final int gene=context.gene, fold=population.split().validationByGene[gene] ? 1 : 0;
						check(!seen[gene] && context.status==tables.statusByGene[gene], "Duplicate context or changed native availability: "+population.ids[gene]);
						seen[gene]=true;
						accounting[fold][0]++;
						if(!population.eligible(gene)){accounting[fold][1]++;}
						else if(!context.hasFeatures()){accounting[fold][context.status+1]++;}
						else{
							accounting[fold][5]++;
							counts[fold*2].observe(read.bases, context.start, context.stop);
							counts[fold*2+1].observe(read.bases, context.start, context.stop);
						}
					}
				}
			}catch(Throwable failure){KillSwitch.exceptionKill(failure);}
		}finally{
			Shared.TRIM_READ_DESCRIPTION=trim;
			Read.TO_UPPER_CASE=upper;
		}
		for(final boolean found:seen){check(found, "Every original locus must remain in screening accounting");}
		for(final long[] a:accounting){check(a[0]==a[1]+a[2]+a[3]+a[4]+a[5], "Every original context has one eligibility/availability outcome");}
		return accounting;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Count Discrimination     ----------------*/
	/*--------------------------------------------------------------*/

	/** Primitive accumulator for one geometry/end/partition. Borrows a completed table.
	 * Each observe call represents one eligible feature-complete locus. Missing words
	 * and impossible candidate spans are excluded as competitors and counted explicitly.
	 * The curve denominator needs a valid true word and at least one valid other word.
	 */
	public static final class Counts{
		public Counts(final NcrnaEndpointTable table_, final boolean validation_, final int radius_, final double[] ratios_){
			check(table_!=null && radius_>=1 && radius_<Integer.MAX_VALUE, "Screening requires a completed table and positive finite radius");
			validateRatios(ratios_);
			table=table_;
			validation=validation_;
			radius=radius_;
			ratios=ratios_.clone();
			correct=new long[ratios.length];
			wrong=new long[ratios.length];
			loo=validation ? null : table.leaveOneOut();
		}

		/** Counts raw words at all valid candidate endpoints; +1 smoothing is applied only for ratios. */
		public void observe(final byte[] bases, final int start, final int stop){
			check(bases!=null && start>=0 && stop<bases.length && (long)stop-start+1>=16,
				"Screened contexts must satisfy native vector span geometry before candidate shifts");
			loci++;
			final boolean isStop=table.type==BoundaryType.STOP;
			if(loo!=null){loo.set(bases, isStop ? stop : start);}
			final NinemerTable counts=validation ? table : loo;
			int trueCount=-1, other=-1, trueStatus=VALID;
			for(int offset=-radius; offset<=radius; offset++){
				attempts++;
				final long left=(long)start+(isStop ? 0 : offset), right=(long)stop+(isStop ? offset : 0);
				if(left<0 || right>=bases.length || right-left+1<16){
					invalidSpan++;
					continue;
				}
				final int boundary=(int)(isStop ? right : left);
				final int position=TrnaBoundaryFeatures.kmerWindowStart(boundary, table.type, table.inside, table.outside);
				final int status=wordStatus(bases, position, table.k);
				if(offset==0){trueStatus=status;}
				if(status==MISSING){missingWords++;}
				else if(status==AMBIGUOUS){ambiguousWords++;}
				else{
					validWords++;
					final int raw=counts.rawCount(bases, position);
					assert(raw>=0) : "The native table/LOO contract supplies nonnegative counts for valid words: "+raw;
					if(offset==0){trueCount=raw;}
					else{other=Math.max(other, raw);}
				}
			}
			if(trueStatus==MISSING){
				trueMissing++;
				return;
			}
			if(trueStatus==AMBIGUOUS){
				trueAmbiguous++;
				return;
			}
			assert(trueCount>=0) : "Offset zero must survive the checked native context span";
			if(other<0){
				noNeighbors++;
				return;
			}
			eligible++;
			trueSum+=trueCount;
			otherSum+=other;
			if(trueCount>other){strictBest++;}
			else if(trueCount==other){tiedBest++;}
			else{notBest++;}
			final double center=(double)trueCount+1, neighbor=(double)other+1;
			for(int i=0; i<ratios.length; i++){
				if(center>=ratios[i]*neighbor){correct[i]++;}
				if(neighbor>=ratios[i]*center){wrong[i]++;}
			}
		}

		/** Checks independent locus and candidate conservation before publishing counters. */
		public void validate(){
			check(loci==eligible+trueMissing+trueAmbiguous+noNeighbors, "Every assessed locus needs exactly one screening disposition");
			check(attempts==loci*(2L*radius+1) && attempts==validWords+invalidSpan+missingWords+ambiguousWords,
				"Every attempted offset must be valid, outside the candidate span, missing a word, or ambiguous");
			check(eligible==strictBest+tiedBest+notBest, "Best-count outcomes must conserve the curve denominator");
			for(int i=0; i<ratios.length; i++){check(correct[i]+wrong[i]<=eligible, "Ratios above one make correct and wrong decisions mutually exclusive");}
		}

		public final NcrnaEndpointTable table;
		public final boolean validation;
		public final int radius;
		private final double[] ratios;
		private final NcrnaEndpointTable.LeaveOneOut loo;
		/** Curve counts aligned with the supplied ordered ratios; caller must not mutate. */
		public final long[] correct, wrong;
		public long loci, eligible, trueMissing, trueAmbiguous, noNeighbors;
		public long strictBest, tiedBest, notBest, trueSum, otherSum;
		public long attempts, validWords, invalidSpan, missingWords, ambiguousWords;
	}

	/** Availability only; actual packing/counts remain owned by the native table. */
	private static int wordStatus(final byte[] bases, final int position, final int k){
		if(position<0 || position>bases.length-k){return MISSING;}
		for(int i=position; i<position+k; i++){
			final int b=bases[i];
			if(b<0 || b>=AminoAcid.baseToNumber.length || AminoAcid.baseToNumber[b]<0){return AMBIGUOUS;}
		}
		return VALID;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final Counts c, final long[] a, final long learning, final Output summary, final Output curves){
		c.validate();
		check(c.loci==a[5], "Both endpoints must assess every feature-complete eligible locus");
		final ByteBuilder b=prefix(c);
		for(final long value:a){b.tab().append(value);}
		b.tab().append(learning).tab().append(c.table.counted()).tab().append(c.table.missing()).tab().append(c.table.ambiguous())
			.tab().append(c.eligible).tab().append(c.trueMissing).tab().append(c.trueAmbiguous).tab().append(c.noNeighbors)
			.tab().append(c.strictBest).tab().append(c.tiedBest).tab().append(c.notBest).tab().append(c.trueSum).tab().append(c.otherSum)
			.tab().append(c.attempts).tab().append(c.validWords).tab().append(c.invalidSpan).tab().append(c.missingWords).tab().append(c.ambiguousWords);
		summary.writer.print(b.nl());
		for(int ratio=0; ratio<ratios.length; ratio++){
			final ByteBuilder row=prefix(c).tab().append(ratioLabels[ratio]).tab().append(c.eligible)
				.tab().append(c.correct[ratio]).tab().append(c.wrong[ratio]).tab().append(c.eligible-c.correct[ratio]-c.wrong[ratio]);
			if(c.eligible==0){row.append("\tNA\tNA");}
			else{row.tab().append(c.correct[ratio]/(double)c.eligible, 8).tab().append(c.wrong[ratio]/(double)c.eligible, 8);}
			curves.writer.print(row.nl());
		}
	}

	private ByteBuilder prefix(final Counts c){
		return buffer.clear().append(c.table.k).tab().append(c.table.inside).tab().append(c.table.outside).tab()
			.append(c.table.type==BoundaryType.START ? "START" : "STOP").tab()
			.append(c.validation ? "VALIDATION" : "TRAIN_LOO").tab().append(radius);
	}

	private static void validateRatios(final double[] values){
		check(values!=null && values.length>0, "Supply at least one decisiveness ratio");
		for(int i=0; i<values.length; i++){
			check(Double.isFinite(values[i]) && values[i]>1 && (i==0 || values[i]>values[i-1]),
				"Ratios must be finite, greater than one and strictly increasing");
		}
	}

	private String file(final String name){return new File(root, name).getPath();}
	private static void check(final boolean condition, final String reason){if(!condition){throw new IllegalArgumentException(reason);}}

	private String contexts, graph, original, roster, loci;
	private long seed=-1;
	private int radius=10;
	private int[] ks={7, 9};
	private double[] ratios={1.5, 2, 3, 5, 10};
	/** One retained label per CLI value; report rounding must not relabel the measured threshold. */
	private final String[] ratioLabels;
	private final int workers;
	private final File root;
	private NcrnaLearningPopulation population;
	private final ByteBuilder buffer=new ByteBuilder();
	private static final int VALID=0, MISSING=1, AMBIGUOUS=2;
	private static final String SUMMARY_HEADER="k\tinside\toutside\tboundary\tpartition\tradius\toriginal_loci\tfiltered\tno_flank\tmalformed_span\tno_gc\tassessed_loci"
		+"\ttable_learning\ttable_counted\ttable_missing\ttable_ambiguous\teligible_loci\ttrue_missing\ttrue_ambiguous\tno_neighbors"
		+"\tstrict_best\ttied_best\tnot_best\ttrue_count_sum\tmax_other_count_sum\tattempted_offsets\tvalid_words\tinvalid_spans\tmissing_words\tambiguous_words";
}
