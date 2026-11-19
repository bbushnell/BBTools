package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import consensus.BaseGraph;
import consensus.BaseGraphHelper;
import map.LongHashSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaBoundaryJitter;
import prok.NcrnaSignedEndpointErrors;
import prok.NcrnaFamily;
import prok.NcrnaLearningPopulation;
import prok.NcrnaPositionalTables;
import prok.NcrnaPositionalBoundaryFeatures;
import prok.NcrnaPositionalVectors;
import prok.NcrnaFixedEndpointTrials;
import prok.RrnaPositionalKmerTable;
import prok.NcrnaSupportIO;
import prok.NcrnaTextIO.Input;
import shared.Shared;
import structures.ByteBuilder;

/** Literal vector dispositions, signed labels and parallel ordered-output oracle.
 * Only observed errors -1/+2 are sampled, so correction classes must be R+1/R-2.
 * One-worker and three-worker decompressed matrices must agree row-for-row.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaPositionalVectorsTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaPositionalVectorsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown positional-vector fixture option: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.maxReads<0 && parser.samplerate==1
				&& !parser.overwrite && !parser.append && new File(parser.out1).mkdir(), "Supply out=fresh_literal_fixture_directory");
			run(new File(parser.out1));
			pp.outstream.println("NCRNA_POSITIONAL_VECTORS_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final File root){
		final String roster=path(root, "roster.tsv"), loci=path(root, "loci.tsv"), graph=path(root, "graph.tsv");
		final String exclude=path(root, "exclude.txt"), fasta=path(root, "contexts.fa"), histogram=path(root, "errors.tsv");
		NcrnaLearningPopulationTest.write(roster, NcrnaLearningPopulationTest.roster());
		NcrnaLearningPopulationTest.write(loci, NcrnaLearningPopulationTest.loci());
		NcrnaLearningPopulationTest.write(graph, NcrnaLearningPopulationTest.graph());
		NcrnaLearningPopulationTest.write(exclude, "gene0\ngene2\n");
		NcrnaLearningPopulationTest.write(fasta, corpus());
		NcrnaLearningPopulationTest.write(histogram, "family\tend\tmeasure\toffset_bp\tloci\n"
			+"RF00174\t5prime\tsigned\t-1\t1\nRF00174\t5prime\tsigned\t2\t1\n"
			+"RF00174\t3prime\tsigned\t-1\t1\nRF00174\t3prime\tsigned\t2\t1\n");
		final NcrnaLearningPopulation population=NcrnaLearningPopulation.loadForRefit(NcrnaSupportIO.read(graph),
			new String[]{"common", "rare"}, roster, loci, exclude);
		final int[] radii={6, 10, 16};
		final NcrnaPositionalTables tables=NcrnaPositionalTables.build(population, fasta, "RF00174", radii, 1.0/262144);
		final NcrnaBoundaryJitter[] jitter=new NcrnaBoundaryJitter[3];
		for(int r=0; r<3; r++){jitter[r]=NcrnaBoundaryJitter.load(histogram, "RF00174", radii[r]);}
		final File first=new File(root, "one_worker"), parallel=new File(root, "three_workers");
		check(first.mkdir() && parallel.mkdir(), "Both finite output roots must be fresh");
		final NcrnaFamily family=family("cobalamin_rf00174");
		final NcrnaPositionalVectors one=NcrnaPositionalVectors.write(population, fasta, family, tables, jitter, 5, 74003, 74004, 1, first);
		final NcrnaPositionalVectors three=NcrnaPositionalVectors.write(population, fasta, family, tables, jitter, 5, 74003, 74004, 3, parallel);
		check(one.records==21 && one.filtered==3 && one.unavailable==3 && one.noModel==1 && one.processed==14,
			"All21originals remain:3filtered,3unavailable,1all-N/no-shortlist,14model-selected loci");
		check(three.records==one.records && three.processed==one.processed && three.noModel==one.noModel,
			"Worker count cannot change original-locus dispositions");
		for(int r=0; r<3; r++){
			for(int set=0; set<2; set++){
				for(int end=0; end<2; end++){
					final int file=NcrnaPositionalVectors.index(r, set, end);
					check(one.attempts[file]==70 && one.rows[file]==(end==0 ? 65 : 70)
						&& one.missingWord[file]==(end==0 ? 5 : 0) && one.invalidSpan[file]==0,
						"Fourteen processed loci give70trials/end/set/radius; the ambiguous start loses5rows only at that end");
					matrix(one.paths[file], radii[r], one.rows[file]);
					equalRows(one.paths[file], three.paths[file]);
				}
			}
		}
		equalRows(one.lociPath, three.lociPath);
		boolean noModel=false;
		try(Input input=new Input(one.lociPath)){
			final LineParser1 fields=new LineParser1('\t');
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				fields.set(line);
				if(fields.termEquals("gene4", 0)){
					check(fields.parseInt(2)==-1 && fields.termEquals("NO_SHORTLISTED_MODEL", 5), "The all-N original must remain in the locus ledger, not become a fake feature row");
					noModel=true;
				}
			}
		}
		check(noModel, "No-model loci cannot disappear from the original axis");
		final File wrong=new File(root, "wrong_member");
		check(wrong.mkdir(), "Wrong-member fixture root is fresh");
		boolean rejected=false;
		try{NcrnaPositionalVectors.write(population, fasta, family("cobalamin_rf01482"), tables, jitter, 5, 74003, 74004, 1, wrong);}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("registered member");}
		check(rejected && !new File(wrong, "vectors.stats.tsv").exists(), "Wrong family/table binding must fail before data output");
		final String ungroupedHistogram=path(root, "sam1_errors.tsv");
		NcrnaLearningPopulationTest.write(ungroupedHistogram, "family\tend\tmeasure\toffset_bp\tloci\n"
			+"sam1\t5prime\tsigned\t-1\t1\nsam1\t3prime\tsigned\t-1\t1\n");
		final NcrnaPositionalTables ungroupedTables=NcrnaPositionalTables.build(population, fasta, "sam1", new int[]{1}, 1.0/262144);
		final NcrnaBoundaryJitter[] ungroupedJitter={NcrnaBoundaryJitter.load(ungroupedHistogram, "sam1", 1)};
		final File ungroupedRoot=new File(root, "ungrouped");
		check(ungroupedRoot.mkdir(), "Ungrouped binding check needs fresh outputs");
		final NcrnaPositionalVectors ungrouped=NcrnaPositionalVectors.write(population, fasta, family("sam1"),
			ungroupedTables, ungroupedJitter, 5, 74003, 74004, 1, ungroupedRoot);
		check(ungrouped.records==21 && ungrouped.processed==14 && ungrouped.rows[0]==65 && ungrouped.rows[1]==70,
			"Ungrouped families use their stable registered ID, not the registry's empty group-member field");
		matrix(ungrouped.paths[0], 1, 65);
		matrix(ungrouped.paths[1], 1, 70);
		final NcrnaPositionalTables tripleTables=NcrnaPositionalTables.buildTriple(population, fasta, "RF00174",
			new int[]{10}, 7, 2, 8, 1, 1);
		final NcrnaBoundaryJitter[] tripleJitter={jitter[1]};
		final File tripleOne=new File(root, "triple_one"), tripleParallel=new File(root, "triple_parallel");
		check(tripleOne.mkdir() && tripleParallel.mkdir(), "Optional vector roots must be fresh");
		final NcrnaFamily tripleFamily=family("cobalamin_rf00174", true);
		final NcrnaPositionalVectors optional=NcrnaPositionalVectors.write(population, fasta, tripleFamily,
			tripleTables, tripleJitter, 5, 74003, 74004, 1, tripleOne);
		final NcrnaPositionalVectors optionalParallel=NcrnaPositionalVectors.write(population, fasta, tripleFamily,
			tripleTables, tripleJitter, 5, 74003, 74004, 3, tripleParallel);
		check(optional.records==21 && optional.filtered==3 && optional.unavailable==3 && optional.noModel==1
			&& optional.processed==14 && optionalParallel.processed==14,
			"Optional LOO and HBM features must preserve every original disposition and explicit unavailable/no-model locus");
		for(int set=0; set<2; set++){
			for(int end=0; end<2; end++){
				final int file=NcrnaPositionalVectors.index(0, set, end);
				check(optional.attempts[file]==70 && optional.rows[file]==(end==0 ? 65 : 70)
					&& optional.missingWord[file]==(end==0 ? 5 : 0) && optional.invalidSpan[file]==0,
					"Both optional training and same-locus monitoring retain the literal70trials and ambiguous-word limitations");
				matrix(optional.paths[file], 10, optional.rows[file], true);
				equalRows(optional.paths[file], optionalParallel.paths[file]);
			}
		}
		equalRows(optional.lociPath, optionalParallel.lociPath);
		refit(root, graph, roster, loci, fasta);
		interchange(root, graph, histogram);
	}

	/** Final refit includes formerly excluded complete loci and emits no outer matrices. */
	private static void refit(final File root, final String graph, final String roster, final String loci, final String fasta){
		final String histogram=path(root, "refit.errors.tsv");
		NcrnaLearningPopulationTest.write(histogram, "family\tend\tmeasure\toffset_bp\tloci\n"
			+"RF00174\t5prime\tsigned\t-11\t1\nRF00174\t5prime\tsigned\t0\t1\nRF00174\t5prime\tsigned\t11\t1\n"
			+"RF00174\t3prime\tsigned\t-11\t1\nRF00174\t3prime\tsigned\t0\t1\nRF00174\t3prime\tsigned\t11\t1\n");
		final NcrnaLearningPopulation population=NcrnaLearningPopulation.loadForRefit(NcrnaSupportIO.read(graph),
			new String[]{"common", "rare"}, roster, loci, null);
		final NcrnaSignedEndpointErrors errors=NcrnaSignedEndpointErrors.load(histogram, "RF00174");
		final NcrnaPositionalTables tables=NcrnaPositionalTables.buildFull(population, fasta, "RF00174",
			new int[]{10}, 7, 2, 8, 1, 1, errors);
		final NcrnaFamily family=family("cobalamin_rf00174", true);
		final File first=new File(root, "refit_one"), parallel=new File(root, "refit_three");
		check(first.mkdir() && parallel.mkdir(), "Refit fixture outputs must be fresh");
		final NcrnaPositionalVectors one=NcrnaPositionalVectors.writeRefit(population, fasta, family, tables, 5, 74003, 74004, 1, first);
		final NcrnaPositionalVectors three=NcrnaPositionalVectors.writeRefit(population, fasta, family, tables, 5, 74003, 74004, 3, parallel);
		check(one.records==21 && one.filtered==1 && one.processed==16 && one.paths.length==4,
			"Only the original partial locus is filtered; formerly withheld genes0/2 rejoin the sixteen usable refit loci");
		for(int file=0; file<4; file++){
			check(one.attempts[file]==80, "Sixteen refit loci each contribute five attempted draws to each end and set");
			equalRows(one.paths[file], three.paths[file]);
			try(Input matrix=new Input(one.paths[file])){
				check(Arrays.equals(matrix.reader.nextLine(), "#dims\t69\t22\t0".getBytes(StandardCharsets.US_ASCII)),
					"Product refit retains the full-error interchange input/output schema");
			}
		}
		check(!new File(first, "r10.evaluation.5prime.tsv").exists() && !new File(first, "r10.evaluation.3prime.tsv").exists(),
			"No evaluation data are invented when every complete eligible locus teaches");
		population.excluded[0]=true;
		final File rejectedRoot=new File(root, "refit_excluded");
		check(rejectedRoot.mkdir(), "Negative refit fixture root must be fresh");
		boolean rejected=false;
		try{NcrnaPositionalVectors.writeRefit(population, fasta, family, tables, 5, 74003, 74004, 1, rejectedRoot);}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("cannot withhold");}
		check(rejected && !new File(rejectedRoot, "r10.train.5prime.tsv").exists(), "Hidden refit exclusions fail before matrices are opened");
		population.excluded[0]=false;
	}

	/** A held-out-only C word must neither teach counts nor be subtracted from them.
	 * This would fail the LOO contribution precondition if evaluation used that view.
	 * Missing taxonomy is retained; row sidecars independently reconstruct features
	 * and labels and reject any cross-split locus or duplicate sample.
	 */
	private static void interchange(final File root, final String graph, final String histogram){
		final String roster=path(root, "exchange.roster.tsv"), loci=path(root, "exchange.loci.tsv");
		final String exclude=path(root, "exchange.exclude.txt"), fasta=path(root, "exchange.fa");
		NcrnaLearningPopulationTest.write(roster, NcrnaLearningPopulationTest.roster().replace("assembly2\t3\t", "assembly2\t.\t"));
		NcrnaLearningPopulationTest.write(loci, NcrnaLearningPopulationTest.loci().replace("assembly2\t3\t", "assembly2\t.\t"));
		NcrnaLearningPopulationTest.write(exclude, "gene0\ngene2\n");
		NcrnaLearningPopulationTest.write(fasta, corpus(true));
		final NcrnaLearningPopulation population=NcrnaLearningPopulation.loadForRefit(NcrnaSupportIO.read(graph),
			new String[]{"common", "rare"}, roster, loci, exclude);
		final NcrnaPositionalTables tables=NcrnaPositionalTables.buildTriple(population, fasta, "RF00174",
			new int[]{10}, 7, 2, 8, 1, 1);
		final int unique=RrnaPositionalKmerTable.encode("CGAAAAAAA".getBytes(StandardCharsets.US_ASCII), 0, 9);
		for(int site=0; site<21; site++){
			check(tables.five[0].counts[site][unique]==0, "The literal held-out-only word must have zero teaching votes at every offset");
		}
		final NcrnaBoundaryJitter[] jitter={NcrnaBoundaryJitter.load(histogram, "RF00174", 10)};
		final NcrnaFamily family=family("cobalamin_rf00174", true);
		final File oneRoot=new File(root, "exchange_one"), threeRoot=new File(root, "exchange_three");
		check(oneRoot.mkdir() && threeRoot.mkdir(), "Interchange fixture directories must be fresh");
		final NcrnaPositionalVectors one=NcrnaPositionalVectors.writeInterchange(population, fasta, family, tables,
			jitter, 5, 74003, 74004, 74005, 1, oneRoot);
		final NcrnaPositionalVectors three=NcrnaPositionalVectors.writeInterchange(population, fasta, family, tables,
			jitter, 5, 74003, 74004, 74005, 3, threeRoot);
		check(one.records==21 && one.filtered==3 && one.unavailable==3 && one.noModel==1 && one.processed==14
			&& one.evaluationLoci==2 && one.evaluationProcessed==2 && one.evaluationUnavailable==0 && one.evaluationNoModel==0,
			"Training dispositions stay unchanged while two complete excluded loci are evaluated, including missing taxonomy");
		for(int set=0; set<3; set++){
			for(int end=0; end<2; end++){
				final int file=one.matrixIndex(0, set, end);
				check(one.attempts[file]==(set==2 ? 10 : 70) && one.rows[file]==(set==2 ? 10 : end==0 ? 65 : 70)
					&& one.clippedWord[file]==0 && one.invalidSpan[file]==0,
					"The separate population axes retain literal sample counts and missing-word accounting");
				matrix(one.paths[file], 10, one.rows[file], true);
				alignedRows(one.paths[file], one.sidecarPaths[file], set, end, population, family, tables);
				equalRows(one.paths[file], three.paths[file]);
				equalRows(one.sidecarPaths[file], three.sidecarPaths[file]);
			}
		}
		equalRows(one.lociPath, three.lociPath);
		final String clippedHistogram=path(root, "exchange.clipped.errors.tsv");
		NcrnaLearningPopulationTest.write(clippedHistogram, "family\tend\tmeasure\toffset_bp\tloci\n"
			+"RF00174\t5prime\tsigned\t-2\t1\nRF00174\t5prime\tsigned\t0\t1\n"
			+"RF00174\t3prime\tsigned\t0\t2\n");
		final NcrnaBoundaryJitter[] clippedJitter={NcrnaBoundaryJitter.load(clippedHistogram, "RF00174", 10)};
		final File clippedRoot=new File(root, "exchange_clipped");
		check(clippedRoot.mkdir(), "Clipped-score-window fixture directory must be fresh");
		final NcrnaPositionalVectors clipped=NcrnaPositionalVectors.writeInterchange(population, fasta, family, tables,
			clippedJitter, 5, 74003, 74004, 74005, 1, clippedRoot);
		for(int set=0; set<3; set++){
			final int start=clipped.matrixIndex(0, set, 0), stop=clipped.matrixIndex(0, set, 1);
			check(clipped.clippedWord[start]>0 && clipped.rows[start]>0 && clipped.clippedWord[stop]==0
				&& clipped.rows[stop]==(set==2 ? 10 : 70), "Only the left score window can clip for the literal18-or20start/119stop candidates");
			try(Input input=new Input(clipped.sidecarPaths[start])){
				final LineParser1 fields=new LineParser1('\t');
				input.reader.nextLine();
				for(byte[] line; (line=input.reader.nextLine())!=null;){
					fields.set(line);
					check(fields.parseInt(7)==20, "Raw start18 clips the exclusive-left score window and must not acquire a fake feature row");
				}
			}
		}
		// An alias under another source name must not straddle the split.
		final String overlapRoster=path(root, "overlap.roster.tsv"), overlapLoci=path(root, "overlap.loci.tsv");
		NcrnaLearningPopulationTest.write(overlapRoster, NcrnaLearningPopulationTest.roster().replace("assembly19\t20\t", "assembly0\t1\t"));
		NcrnaLearningPopulationTest.write(overlapLoci, NcrnaLearningPopulationTest.loci().replace("assembly19\t20\t", "assembly0\t1\t"));
		final NcrnaLearningPopulation overlap=NcrnaLearningPopulation.loadForRefit(NcrnaSupportIO.read(graph),
			new String[]{"common", "rare"}, overlapRoster, overlapLoci, exclude);
		final NcrnaPositionalTables overlapTables=NcrnaPositionalTables.buildTriple(overlap, fasta, "RF00174", new int[]{10}, 7, 2, 8, 1, 1);
		final File rejectedRoot=new File(root, "exchange_overlap");
		check(rejectedRoot.mkdir(), "Overlap preflight fixture directory must be fresh");
		boolean rejected=false;
		try{NcrnaPositionalVectors.writeInterchange(overlap, fasta, family, overlapTables, jitter, 5, 74003, 74004, 74005, 1, rejectedRoot);}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("share assembly");}
		check(rejected && !new File(rejectedRoot, "r10.train.5prime.tsv").exists(), "Assembly aliases across splits must fail before vector output");
		fullErrors(root, population, fasta, family);
	}

	/** Full signed draws include usable distant truth, rather than forcing a 21-way answer. */
	private static void fullErrors(final File root, final NcrnaLearningPopulation population,
			final String fasta, final NcrnaFamily family){
		final String histogram=path(root, "full.errors.tsv");
		NcrnaLearningPopulationTest.write(histogram, "family\tend\tmeasure\toffset_bp\tloci\n"
			+"RF00174\t5prime\tsigned\t-11\t1\nRF00174\t5prime\tsigned\t0\t1\nRF00174\t5prime\tsigned\t11\t1\n"
			+"RF00174\t3prime\tsigned\t-11\t1\nRF00174\t3prime\tsigned\t0\t1\nRF00174\t3prime\tsigned\t11\t1\n");
		final NcrnaSignedEndpointErrors errors=NcrnaSignedEndpointErrors.load(histogram, "RF00174");
		final NcrnaPositionalTables tables=NcrnaPositionalTables.buildFull(population, fasta, "RF00174",
			new int[]{10}, 7, 2, 8, 1, 1, errors);
		final File oneRoot=new File(root, "full_one"), threeRoot=new File(root, "full_three");
		check(oneRoot.mkdir() && threeRoot.mkdir(), "Full-error fixture outputs must be fresh");
		final NcrnaPositionalVectors one=NcrnaPositionalVectors.writeInterchange(population, fasta, family,
			tables, null, 5, 74003, 74004, 74005, 1, oneRoot);
		final NcrnaPositionalVectors three=NcrnaPositionalVectors.writeInterchange(population, fasta, family,
			tables, null, 5, 74003, 74004, 74005, 3, threeRoot);
		long outside=0, zero=0;
		for(int set=0; set<3; set++){
			for(int end=0; end<2; end++){
				final int file=one.matrixIndex(0, set, end);
				equalRows(one.paths[file], three.paths[file]);
				equalRows(one.sidecarPaths[file], three.sidecarPaths[file]);
				check(one.attempts[file]==one.rows[file]+one.missingWord[file]+one.invalidSpan[file]
					&& one.invalidSpan[file]==one.crossedSpan[file]+one.outsideContext[file],
					"Full draws must remain emitted or explicitly unavailable, without redraw or hidden conditioning");
				final LineParser1 data=new LineParser1('\t'), metadata=new LineParser1('\t');
				long rows=0, fileOutside=0, fileZero=0;
				try(Input numbers=new Input(one.paths[file]); Input identities=new Input(one.sidecarPaths[file])){
					check(Arrays.equals(numbers.reader.nextLine(), "#dims\t69\t22\t0".getBytes(StandardCharsets.US_ASCII)),
						"Full-error matrices must explicitly declare the extra answer");
					identities.reader.nextLine();
					for(byte[] line; (line=numbers.reader.nextLine())!=null;){
						data.set(line);
						final byte[] identity=identities.reader.nextLine();
						check(identity!=null, "Every full-error numeric row needs its original identity");
						metadata.set(identity);
						check(data.terms()==91 && metadata.terms()==11, "The optional schema is exactly 69 inputs and 22 one-hot targets");
						final int gene=population.index(metadata.parseString(0));
						check(gene>=0 && !population.partial[gene] && population.excluded[gene]==(set==2),
							"Full-error rows must preserve the exact teaching/evaluation source split");
						final int truth=end==0 ? 20 : 119, raw=metadata.parseInt(end==0 ? 7 : 8);
						final int correction=truth-raw;
						check(correction==-11 || correction==0 || correction==11, "Only original measured signed errors may be drawn");
						final int expected=correction==0 ? 10 : 21;
						for(int feature=0; feature<69; feature++){
							check(Float.isFinite(data.parseFloat(feature)), "An out-of-range label needs usable finite features, not missing geometry");
						}
						for(int target=0; target<22; target++){
							check(data.parseInt(69+target)==(target==expected ? 1 : 0),
								"Reachable truth and either signed tail must have their independently reconstructed one-hot answer");
						}
						if(expected==21){
							outside++;
							fileOutside++;
						}else{
							zero++;
							fileZero++;
						}
						rows++;
					}
					check(identities.reader.nextLine()==null && rows==one.rows[file], "Numeric rows, identities and accounting must end together");
					check(fileOutside==one.outOfRangeWritten[file] && fileZero==one.zeroWritten[file]
						&& one.outOfRangeAttempts[file]>=fileOutside && one.zeroAttempts[file]>=fileZero,
						"Streaming class counts must match independently decoded labels and expose lost attempted mass");
				}
			}
		}
		check(outside>0 && zero>0 && one.clippedWord[0]>0,
			"The fixed toy draws must exercise usable tails, correct endpoints and unavailable clipped words separately");
		final String scores=path(root, "full.weighted_scores.tsv");
		tables.fiveWeighted[0][0].writeScores(scores);
		try(Input input=new Input(scores)){
			check(Arrays.equals(input.reader.nextLine(), "#NCRNA_SIGNED_BACKGROUND_SCORES_V1".getBytes(StandardCharsets.US_ASCII)),
				"Weighted scores must never claim the legacy uniform-background schema");
		}
		boolean rejected=false;
		try{RrnaPositionalKmerTable.loadScores(scores, "RF00174");}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "The legacy caller loader must not silently accept experimental weighted resources");
		fixedOuter(root, population, fasta, family, tables, one);
		invalidFullSpans(root, family);
	}

	/** Zero offsets reproduce original bytes; nonzero centers preserve original trials and explicit holds. */
	private static void fixedOuter(final File root, final NcrnaLearningPopulation population, final String fasta,
			final NcrnaFamily family, final NcrnaPositionalTables tables, final NcrnaPositionalVectors baseline){
		final int five=baseline.matrixIndex(0, 2, 0), three=baseline.matrixIndex(0, 2, 1);
		final String rows5=baseline.sidecarPaths[five], rows3=baseline.sidecarPaths[three];
		final long[] expected={baseline.rows[five], baseline.rows[three]};
		final String zero=path(root, "fixed.zero.offsets.tsv");
		NcrnaLearningPopulationTest.write(zero, offsetRows(0, 0));
		final NcrnaFixedEndpointTrials fixed=NcrnaFixedEndpointTrials.load(population, family, "RF00174",
			rows5, rows3, zero, 10, 5, expected);
		final File zeroRoot=new File(root, "fixed_zero");
		check(zeroRoot.mkdir(), "Zero-offset fixed outputs must be fresh");
		final NcrnaPositionalVectors zeroVectors=NcrnaPositionalVectors.writeInterchange(population, fasta, family,
			tables, null, 5, 74003, 74004, 74005, 3, zeroRoot, fixed);
		for(int file=0; file<baseline.paths.length; file++){
			equalRows(baseline.paths[file], zeroVectors.paths[file]);
			equalRows(baseline.sidecarPaths[file], zeroVectors.sidecarPaths[file]);
		}
		for(int end=0; end<2; end++){fixedLedger(zeroVectors, baseline, end, 0, 0, false);}
		final String unsupportedPath=path(root, "fixed.unsupported_zero.offsets.tsv");
		NcrnaLearningPopulationTest.write(unsupportedPath, offsetRows(0, 0).replace("a\t10\t10", "a\t9\t10"));
		final NcrnaFixedEndpointTrials unsupported=NcrnaFixedEndpointTrials.load(population, family, "RF00174",
			rows5, rows3, unsupportedPath, 10, 5, expected);
		check(!unsupported.supported[0] && unsupported.trainingSupport[0]==9 && unsupported.median5[0]==0 && unsupported.median3[0]==0,
			"Unsupported models are explicitly identifiable and unshifted, not silently assigned an estimated median");
		final String shifted=path(root, "fixed.shifted.offsets.tsv");
		NcrnaLearningPopulationTest.write(shifted, offsetRows(-2, 3));
		final NcrnaFixedEndpointTrials shiftedTrials=NcrnaFixedEndpointTrials.load(population, family, "RF00174",
			rows5, rows3, shifted, 10, 5, expected);
		final File shiftedRoot=new File(root, "fixed_shifted");
		check(shiftedRoot.mkdir(), "Shifted fixed outputs must be fresh");
		final NcrnaPositionalVectors shiftedVectors=NcrnaPositionalVectors.writeInterchange(population, fasta, family,
			tables, null, 5, 74003, 74004, 74005, 1, shiftedRoot, shiftedTrials);
		for(int file=0; file<4; file++){
			equalRows(baseline.paths[file], shiftedVectors.paths[file]);
			equalRows(baseline.sidecarPaths[file], shiftedVectors.sidecarPaths[file]);
		}
		for(int end=0; end<2; end++){
			fixedLedger(shiftedVectors, baseline, end, -2, 3, false);
			final int file=shiftedVectors.matrixIndex(0, 2, end);
			try(Input vectors=new Input(shiftedVectors.paths[file]); Input rows=new Input(shiftedVectors.sidecarPaths[file])){
				vectors.reader.nextLine();
				rows.reader.nextLine();
				final LineParser1 v=new LineParser1('\t'), p=new LineParser1('\t');
				for(byte[] line=vectors.reader.nextLine(); line!=null; line=vectors.reader.nextLine()){
					v.set(line);
					final byte[] metadata=rows.reader.nextLine();
					check(metadata!=null, "Every usable centered vector retains a row-aligned centered sidecar");
					p.set(metadata);
					final int correction=p.parseInt(end==0 ? 9 : 10)-p.parseInt(end==0 ? 7 : 8);
					check(v.terms()==91 && v.parseInt(69+(Math.abs(correction)>10 ? 21 : correction+10))==1,
						"The centered coordinate frame, not the original raw end, defines the new one-hot target");
				}
				check(rows.reader.nextLine()==null, "Centered sidecars and numeric matrices close together");
			}
		}
		for(final int median:new int[]{-1000, 1000}){
			final String badGeometry=path(root, "fixed.geometry"+median+".offsets.tsv");
			NcrnaLearningPopulationTest.write(badGeometry, offsetRows(median, 0));
			final NcrnaFixedEndpointTrials invalid=NcrnaFixedEndpointTrials.load(population, family, "RF00174",
				rows5, rows3, badGeometry, 10, 5, expected);
			final File invalidRoot=new File(root, "fixed_geometry"+median);
			check(invalidRoot.mkdir(), "Invalid-geometry fixed outputs must be fresh");
			final NcrnaPositionalVectors held=NcrnaPositionalVectors.writeInterchange(population, fasta, family,
				tables, null, 5, 74003, 74004, 74005, 3, invalidRoot, invalid);
			for(int end=0; end<2; end++){
				final int file=held.matrixIndex(0, 2, end);
				check(held.rows[file]==0 && held.attempts[file]==expected[end] && held.invalidSpan[file]==expected[end],
					"Every invalid shifted span is retained as a hold, even when no evaluation vector is usable");
				check((median<0 ? held.crossedSpan[file] : held.outsideContext[file])==expected[end],
					"Negative start median crosses the span; positive start median leaves the context");
				fixedLedger(held, baseline, end, median, 0, true);
			}
		}
		final String shortPath=path(root, "fixed.short.offsets.tsv");
		NcrnaLearningPopulationTest.write(shortPath, offsetRows(0, 110));
		final NcrnaFixedEndpointTrials shortTrials=NcrnaFixedEndpointTrials.load(population, family, "RF00174",
			rows5, rows3, shortPath, 10, 5, expected);
		final File shortRoot=new File(root, "fixed_short");
		check(shortRoot.mkdir(), "Below-minimum fixture needs fresh outputs");
		final NcrnaPositionalVectors tooShort=NcrnaPositionalVectors.writeInterchange(population, fasta, family,
			tables, null, 5, 74003, 74004, 74005, 1, shortRoot, shortTrials);
		check(tooShort.belowMinimumSpan[five]+tooShort.belowMinimumSpan[three]>0,
			"The literal110base stop correction must exercise a positive-length span below the loaded16base minimum");
		for(int end=0; end<2; end++){
			final int file=tooShort.matrixIndex(0, 2, end);
			check(tooShort.rows[file]==0 && tooShort.invalidSpan[file]==expected[end]
				&& tooShort.crossedSpan[file]+tooShort.belowMinimumSpan[file]==expected[end],
				"Crossed and positive-but-short centered spans both revert without dropping their original trials");
			fixedLedger(tooShort, baseline, end, 0, 110, true);
		}
		boolean rejected=false;
		try{NcrnaFixedEndpointTrials.load(population, family, "RF00174", rows5, rows3, zero, 10, 5,
			new long[]{expected[0]+1, expected[1]});}
		catch(IllegalArgumentException expectedFailure){rejected=true;}
		check(rejected, "A truncated frozen denominator cannot be accepted");
		final String badSupport=path(root, "fixed.unsupported.offsets.tsv");
		NcrnaLearningPopulationTest.write(badSupport, NcrnaFixedEndpointTrials.OFFSET_HEADER+"\nRF00174\ta\t9\t10\t1\t0\nRF00174\tg\t0\t10\t0\t0\n");
		rejected=false;
		try{NcrnaFixedEndpointTrials.load(population, family, "RF00174", rows5, rows3, badSupport, 10, 5, expected);}
		catch(IllegalArgumentException expectedFailure){rejected=true;}
		check(rejected, "Unsupported models cannot acquire an invented nonzero median");
	}

	/** Verify original metadata byte-for-byte and reconstruct both centered bounds independently. */
	private static void fixedLedger(final NcrnaPositionalVectors centered, final NcrnaPositionalVectors original,
			final int end, final int median5, final int median3, final boolean allHeld){
		final int file=original.matrixIndex(0, 2, end);
		final LineParser1 p=new LineParser1('\t'), baseline=new LineParser1('\t');
		long count=0;
		try(Input audit=new Input(centered.fixedAuditPaths[end]); Input before=new Input(original.sidecarPaths[file])){
			audit.reader.nextLine();
			before.reader.nextLine();
			for(byte[] line=audit.reader.nextLine(); line!=null; line=audit.reader.nextLine()){
				final byte[] raw=before.reader.nextLine();
				check(raw!=null && line.length>raw.length && line[raw.length]=='\t', "Every ledger row extends exactly one original trial");
				for(int i=0; i<raw.length; i++){check(line[i]==raw[i], "Original identities, raw spans and truth bytes never change");}
				p.set(line);
				baseline.set(raw);
				check(p.terms()==21 && p.parseLong(11)==++count && p.parseLong(12)==(long)baseline.parseInt(7)-median5
					&& p.parseLong(13)==(long)baseline.parseInt(8)-median3 && p.parseInt(14)==median5 && p.parseInt(15)==median3,
					"Centering subtracts the declared per-model medians in the same sense frame, retaining ordinal order");
				if(allHeld){check(!p.termEquals("VECTORIZED", 16), "Invalid original trials require explicit holds, never replacement vectors");}
				check(p.parseLong(17)==10 && p.parseInt(18)==10 && p.parseInt(19)==1 && p.parseInt(20)==16,
					"The ledger retains observed support, explicit minimum support and the loaded family's length guard");
			}
			check(before.reader.nextLine()==null && count==original.rows[file], "Every original trial remains, including all invalid holds");
		}
	}

	/** One supported model and one unsupported zero-offset model preserve exact resource order. */
	private static String offsetRows(final int five, final int three){
		return NcrnaFixedEndpointTrials.OFFSET_HEADER+"\nRF00174\ta\t10\t10\t"+five+"\t"+three+"\nRF00174\tg\t0\t10\t0\t0\n";
	}

	/** Longer loci provide measured-offset support; shorter contexts expose actual
	 * noncrossed out-of-context attempts, while independent ends also cross.
	 */
	private static void invalidFullSpans(final File root, final NcrnaFamily family){
		final String roster=path(root, "invalid.roster.tsv"), loci=path(root, "invalid.loci.tsv");
		final String graph=path(root, "invalid.graph.tsv"), exclude=path(root, "invalid.exclude.txt");
		final String fasta=path(root, "invalid.fa"), histogram=path(root, "invalid.errors.tsv");
		NcrnaLearningPopulationTest.write(roster, NcrnaLearningPopulationTest.roster());
		NcrnaLearningPopulationTest.write(graph, NcrnaLearningPopulationTest.graph());
		NcrnaLearningPopulationTest.write(exclude, "gene0\ngene2\ngene20\n");
		final ByteBuilder metadata=new ByteBuilder().append(NcrnaLearningPopulation.NATIVE_HEADER).nl();
		final ByteBuilder sequences=new ByteBuilder();
		for(int gene=0; gene<=20; gene++){
			final int source=gene==20 ? 0 : gene;
			final boolean shortContext=gene>=7 && gene<=11;
			final int left=shortContext ? 20 : 200, right=shortContext ? 30 : 200;
			metadata.append("gene").append(gene).append("\tgene").append(gene).append("\tgenome").append(source)
				.append("\tassembly").append(source).tab().append(source+1).append("\tbacteria\tP\tcontig\t+\t")
				.append(left+1).tab().append(left+100).append("\t100\t").append(left).tab().append(right).nl();
			final byte[] bases=new byte[left+100+right];
			Arrays.fill(bases, 0, left, (byte)'G');
			Arrays.fill(bases, left, left+100, (byte)'A');
			Arrays.fill(bases, left+100, bases.length, (byte)'T');
			sequences.append(">gene").append(gene).append(" source=genome").append(source)
				.append(" lflank=").append(left).append(" rflank=").append(right)
				.append(" contig_gc=").append(Double.toString(left/(double)bases.length)).nl().append(bases).nl();
		}
		NcrnaLearningPopulationTest.write(loci, metadata.toString());
		NcrnaLearningPopulationTest.write(fasta, sequences.toString());
		NcrnaLearningPopulationTest.write(histogram, "family\tend\tmeasure\toffset_bp\tloci\n"
			+"RF00174\t5prime\tsigned\t0\t3\nRF00174\t5prime\tsigned\t110\t3\n"
			+"RF00174\t3prime\tsigned\t-110\t2\nRF00174\t3prime\tsigned\t0\t2\nRF00174\t3prime\tsigned\t80\t2\n");
		final NcrnaLearningPopulation population=NcrnaLearningPopulation.loadForRefit(NcrnaSupportIO.read(graph),
			new String[]{"common", "rare"}, roster, loci, exclude);
		final NcrnaPositionalTables tables=NcrnaPositionalTables.buildFull(population, fasta, "RF00174",
			new int[]{10}, 7, 2, 8, 1, 1, NcrnaSignedEndpointErrors.load(histogram, "RF00174"));
		final File directory=new File(root, "full_invalid");
		check(directory.mkdir(), "Invalid-span full-mode fixture needs a fresh root");
		final NcrnaPositionalVectors vectors=NcrnaPositionalVectors.writeInterchange(population, fasta, family,
			tables, null, 5, 74003, 74004, 74005, 1, directory);
		long crossed=0, outside=0;
		for(int set=0; set<3; set++){
			for(int end=0; end<2; end++){
				final int file=vectors.matrixIndex(0, set, end);
				crossed+=vectors.crossedSpan[file];
				outside+=vectors.outsideContext[file];
				check(vectors.attempts[file]==5*(set==2 ? vectors.evaluationProcessed : vectors.processed)
					&& vectors.invalidSpan[file]==vectors.crossedSpan[file]+vectors.outsideContext[file]
					&& vectors.attempts[file]==vectors.rows[file]+vectors.missingWord[file]+vectors.invalidSpan[file]
					&& vectors.zeroAttempts[file]+vectors.outOfRangeAttempts[file]==vectors.attempts[file]
					&& vectors.zeroWritten[file]+vectors.outOfRangeWritten[file]==vectors.rows[file],
					"Crossed/outside attempts consume their original draw exactly once; they cannot be redrawn into rows");
				long decodedZero=0, decodedOutside=0;
				try(Input numbers=new Input(vectors.paths[file]); Input identities=new Input(vectors.sidecarPaths[file])){
					numbers.reader.nextLine();
					identities.reader.nextLine();
					final LineParser1 fields=new LineParser1('\t'), values=new LineParser1('\t');
					for(byte[] row; (row=numbers.reader.nextLine())!=null;){
						values.set(row);
						final byte[] identity=identities.reader.nextLine();
						check(identity!=null, "Every valid full-error row retains its paired coordinates");
						fields.set(identity);
						final String id=fields.parseString(0);
						final int gene=Integer.parseInt(id.substring(4));
						final boolean shortContext=gene>=7 && gene<=11;
						final int start=fields.parseInt(7), stop=fields.parseInt(8), length=shortContext ? 150 : 500;
						check(start>=0 && stop>=start && stop<length,
							"Neither a crossed nor an out-of-context span may be mislabeled as an out-of-range training row");
						final int correction=fields.parseInt(end==0 ? 9 : 10)-(end==0 ? start : stop);
						check(values.parseInt(69+(Math.abs(correction)>10 ? 21 : correction+10))==1,
							"Only valid geometry may carry its independently reconstructed reachable/out-of-range target");
						if(correction==0){decodedZero++;}
						else{
							check(Math.abs(correction)>10, "This literal fixture has only zero or distant measured errors");
							decodedOutside++;
						}
					}
					check(identities.reader.nextLine()==null, "Invalid spans cannot leave unmatched sidecar rows");
					check(decodedZero==vectors.zeroWritten[file] && decodedOutside==vectors.outOfRangeWritten[file],
						"Decoded valid coordinates/targets must independently recover the written class mixture despite actual geometry losses");
				}
			}
		}
		check(crossed>0 && outside>0, "This deterministic case must exercise both invalid-span classes with positive counts");
	}

	/** Independently joins each numeric row to the source/sample and true correction. */
	private static void alignedRows(final String matrix, final String sidecar, final int set, final int end,
			final NcrnaLearningPopulation population, final NcrnaFamily family, final NcrnaPositionalTables tables){
		final String[] sets={"train", "monitor", "evaluation"};
		final boolean[] seen=new boolean[21*5];
		final LineParser1 data=new LineParser1('\t'), metadata=new LineParser1('\t');
		final NcrnaPositionalBoundaryFeatures full=new NcrnaPositionalBoundaryFeatures(tables.fiveProbes[0], tables.threeProbes[0], false);
		try(Input numbers=new Input(matrix); Input identities=new Input(sidecar)){
			check(Arrays.equals(numbers.reader.nextLine(), "#dims\t69\t21\t0".getBytes(StandardCharsets.US_ASCII)), "Plain delivery must use the agreed unweighted native/GPU header");
			check(identities.reader.nextLine()!=null, "The row sidecar needs a header");
			for(byte[] line; (line=numbers.reader.nextLine())!=null;){
				final byte[] row=identities.reader.nextLine();
				check(row!=null && line.length>0 && line[0]!='#', "Each numeric data row needs exactly one identity row, with no inline metadata");
				data.set(line);
				metadata.set(row);
				check(data.terms()==90 && metadata.terms()==11, "Delivery rows must retain exactly90numbers and the separate eleven identity/coordinate fields");
				final int gene=population.index(metadata.parseString(0)), sample=metadata.parseInt(6);
				check(gene>=0 && sample>=0 && sample<5 && !seen[gene*5+sample], "A locus sample must occur at most once per end/set");
				seen[gene*5+sample]=true;
				check(metadata.termEquals(population.sourceName(gene), 1) && metadata.termEquals(population.assemblyName(gene), 2)
					&& metadata.termEquals(sets[set], 3) && !population.partial[gene]
					&& population.excluded[gene]==(set==2), "Source/assembly/split identities cannot move a locus across the teaching/evaluation boundary");
				final int start=metadata.parseInt(7), stop=metadata.parseInt(8);
				check(metadata.parseInt(9)==20 && metadata.parseInt(10)==119, "Sense truth coordinates must remain the independently known20..119");
				final int label=(end==0 ? 20-start : 119-stop)+10;
				check(data.parseInt(69+label)==1 && metadata.parseInt(4)==0 && metadata.termEquals("a", 5),
					"The one-hot class must invert this row's raw error and retain its selected model");
				if(set==2){
					final byte[] bases=contextBases(gene, true);
					check(full.fill(bases, start, stop, family.library[0], family.models[0], 1f, .5f, 0)==3,
						"Both held-out ends must have full-table features; the literal model has coverage1 at every position");
					for(int feature=0; feature<69; feature++){
						check(Math.abs(data.parseFloat(feature)-full.values[end][feature])<.000001,
							"Every held-out input must equal ordinary inference, with no self-subtraction");
					}
				}
			}
			check(identities.reader.nextLine()==null, "The identity sidecar cannot contain an unpaired trailing row");
		}
	}

	/** An observed error of-1/+2 means a correction of+1/-2, independent of radius. */
	private static void matrix(final String path, final int radius, final long expected){
		matrix(path, radius, expected, false);
	}
	private static void matrix(final String path, final int radius, final long expected, final boolean triple){
		final LineParser1 fields=new LineParser1('\t');
		final int sites=2*radius+1, extras=triple ? 3*sites : sites, inputs=extras+(triple ? 6 : 3);
		long count=0;
		try(Input input=new Input(path)){
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				if(line.length>0 && line[0]=='#'){continue;}
				fields.set(line);
				check(fields.terms()==inputs+sites, "Each positional matrix row must contain the declared feature/class widths");
				for(int feature=0; feature<inputs; feature++){check(Float.isFinite(fields.parseFloat(feature)), "Usable feature rows must be finite");}
				check(fields.parseFloat(extras)>.96f && fields.parseFloat(extras)<1.04f && fields.parseFloat(extras+1)==.5f,
					"Jittered100-base genes retain actual candidate/consensus length and declared whole-contig GC");
				int positive=-1, positives=0;
				for(int site=0; site<sites; site++){
					final int label=fields.parseInt(inputs+site);
					check(label==0 || label==1, "Targets must be explicit one-hot classes");
					if(label==1){
						positive=site;
						positives++;
					}
				}
				check(positives==1 && (positive==radius+1 || positive==radius-2), "True-end class must invert signed raw endpoint error, with exactly one positive");
				count++;
			}
		}
		check(count==expected, "All retained rows must be readable by the actual native numeric parser");
	}
	private static void equalRows(final String first, final String second){
		try(Input a=new Input(first); Input b=new Input(second)){
			for(byte[] left; (left=a.reader.nextLine())!=null;){
				check(Arrays.equals(left, b.reader.nextLine()), "Decompressed matrix/locus order and locus-keyed samples must be independent of worker scheduling");
			}
			check(b.reader.nextLine()==null, "Parallel output cannot add trailing rows");
		}
	}
	private static NcrnaFamily family(final String name){
		return family(name, false);
	}
	private static NcrnaFamily family(final String name, final boolean hbm){
		final byte[][] refs={new byte[100], new byte[100]};
		Arrays.fill(refs[0], (byte)'A');
		Arrays.fill(refs[1], (byte)'G');
		final LongHashSet seeds=new LongHashSet();
		seeds.add(0);
		final BaseGraph[] models=hbm ? new BaseGraph[2] : null;
		if(hbm){
			for(int model=0; model<2; model++){
				models[model]=new BaseGraph(model==0 ? "a" : "g", refs[model], null, model, 0);
				for(int position=0; position<100; position++){
					BaseGraphHelper.setCounts(models[model], position, model==0 ? 1 : 0, 0, model==1 ? 1 : 0, 0, 0);
				}
				BaseGraphHelper.initForScoring(models[model]);
			}
		}
		return new NcrnaFamily(name, refs, models, new String[]{"a n=1", "g n=1"}, seeds, 17, 16, 64, 6, 2, false, 1, 0, 0, 1, 0, 1, .6f, .6f);
	}
	private static String corpus(){
		return corpus(false);
	}
	private static String corpus(final boolean uniqueHeldoutWord){
		final ByteBuilder result=new ByteBuilder();
		for(int gene=0; gene<21; gene++){
			result.append(">gene").append(gene).append(" source=genome").append(gene==20 ? 0 : gene).append(" lflank=20");
			if(gene!=5){result.append(" rflank=30");}
			if(gene!=1){result.append(" contig_gc=0.5");}
			result.nl();
			result.append(contextBases(gene, uniqueHeldoutWord)).nl();
		}
		return result.toString();
	}
	/** Literal sequence source reused only to reconstruct a fixture's known bases. */
	private static byte[] contextBases(final int gene, final boolean uniqueHeldoutWord){
		final byte[] bases=new byte[gene==6 ? 10 : 150];
		if(gene==6){Arrays.fill(bases, (byte)'A');}
		else{
			Arrays.fill(bases, 0, 20, (byte)'G');
			Arrays.fill(bases, 20, 120, gene==4 ? (byte)'N' : (byte)'A');
			Arrays.fill(bases, 120, 150, (byte)'T');
			if(gene==3){bases[20]='N';}
			if(gene==0 && uniqueHeldoutWord){bases[18]='C';}
		}
		return bases;
	}
	private static String path(final File root, final String name){return new File(root, name).getPath();}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
}
