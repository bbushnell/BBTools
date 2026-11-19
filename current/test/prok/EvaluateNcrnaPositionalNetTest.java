package test.prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.EvaluateNcrnaPositionalNet;
import prok.EvaluateNcrnaPositionalNet.Result;
import prok.NcrnaTextIO;
import prok.ProjectNcrnaPositionalInputs;
import prok.RrnaPositionalKmerTable;
import shared.Shared;
import structures.ByteBuilder;

/** Literal class/coordinate oracles using actual constant and normalized dense nets.
 * These small rows test reading, inference and denominators, not biological benefit.
 * @author Brian Bushnell, Raiden
 */
public final class EvaluateNcrnaPositionalNetTest{

	/** Accepts only a fresh output directory for finite generated inputs. */
	public static void main(final String[] args){
		final PreParser pre=new PreParser(args, EvaluateNcrnaPositionalNetTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pre.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown fixture argument "+arg);
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && !parser.append
				&& !parser.overwrite && parser.maxReads<0 && parser.samplerate==1 && new File(parser.out1).mkdir(),
				"Require out=fresh_fixture_directory without input, append or filtering");
			run(new File(parser.out1));
			pre.outstream.println("NCRNA_POSITIONAL_NET_EVALUATION_TEST_PASS");
		}finally{Shared.closeStream(pre.outstream);}
	}

	/** One literal four-row oracle checks ties, signs and initially-correct denominator. */
	private static void run(final File directory){
		final CellNet tied=constant();
		tied.net[1][2].setBias(1, true);
		tied.net[1][3].setBias(1, true);
		final int[] truth={2, 3, 5, 10};
		compare(directory, truth);
		project(directory);
		exponents(directory);
		for(int end=0; end<2; end++){
			final String matrix=write(directory, "matrix"+end+".tsv", matrix(truth));
			final String rows=write(directory, "rows"+end+".tsv", sidecar(truth, end, "evaluation"));
			final Result result=EvaluateNcrnaPositionalNet.evaluate(tied, matrix, rows, end);
			check(result.trials==4 && result.exact==1 && result.within1==2 && result.within3==3,
				"First tied maximum is class2: literal errors0,-1,-3,-8 give exact1/near1=2/near3=3");
			check(result.errorSum==-12 && result.squaredError==74, "Prediction minus truth determines error sign and squared magnitude");
			check(result.alreadyCorrect==1 && result.movedAlreadyCorrect==1,
				"Only trueClass10 is initially correct; movement denominator is1, not all4trials");
			final LineParser1 summary=summary(result, end);
			check(summary.parseLong(8)==1 && summary.parseLong(9)==1 && summary.parseDouble(10)==100,
				"Movement of the one correct end is100percent regardless of numeric padding");
		}
		final int[] noneCorrect={2, 3};
		final String matrix=write(directory, "none.tsv", matrix(noneCorrect));
		final String rows=write(directory, "none.rows.tsv", sidecar(noneCorrect, 0, "evaluation"));
		final LineParser1 emptyDenominator=summary(EvaluateNcrnaPositionalNet.evaluate(tied, matrix, rows, 0), 0);
		check(emptyDenominator.parseLong(8)==0 && emptyDenominator.parseLong(9)==0 && emptyDenominator.termEquals("NA", 10),
			"No initially correct trials gives unavailable movement percentage");

		final CellNet normalized=constant();
		normalized.net[1][10].setBias(5, true);
		normalized.net[1][11].weights[0]=1;
		final float[] mean=new float[69], inverse=new float[69];
		Arrays.fill(inverse, 1);
		mean[0]=1;
		inverse[0]=2;
		normalized.setInputNormalization(mean, inverse);
		final String correctMatrix=write(directory, "correct.tsv", matrix(new int[]{10}));
		final String correctRows=write(directory, "correct.rows.tsv", sidecar(new int[]{10}, 0, "evaluation"));
		check(EvaluateNcrnaPositionalNet.evaluate(normalized, correctMatrix, correctRows, 0).exact==1
			&& normalized.getOutput(11)==4, "Raw3 must normalize once to4; twice would yield6 and choose wrong class11");
		final String exported=write(directory, "normalized.bbnet", normalized.toBytes().toString());
		final CellNet loaded=CellNetParser.load(exported, false);
		check(EvaluateNcrnaPositionalNet.evaluate(loaded, correctMatrix, correctRows, 0).exact==1
			&& loaded.getOutput(11)==4, "Actual exported network must retain normalization and finite output");
		fails(()->EvaluateNcrnaPositionalNet.evaluate(new CellNet(new int[]{68,21}, 58, 1, 0, 1, null), matrix, rows, 0));
		final CellNet invalid=constant();
		invalid.net[1][0].bias=Float.NaN;
		failsNonfinite(()->EvaluateNcrnaPositionalNet.evaluate(invalid, matrix, rows, 0));
		bad(directory, tied, "width", matrix(noneCorrect).replaceFirst("3\\t0", "3"), sidecar(noneCorrect, 0, "evaluation"));
		bad(directory, tied, "nonfinite", matrix(noneCorrect).replaceFirst("3\\t", "NaN\t"), sidecar(noneCorrect, 0, "evaluation"));
		bad(directory, tied, "header", matrix(noneCorrect).replace("69\t21\t0", "69\t21\t1"), sidecar(noneCorrect, 0, "evaluation"));
		bad(directory, tied, "training", matrix(noneCorrect), sidecar(noneCorrect, 0, "train"));
		bad(directory, tied, "labels", matrix(noneCorrect), sidecar(new int[]{3, 2}, 0, "evaluation"));
		bad(directory, tied, "short", matrix(noneCorrect), sidecar(new int[]{2}, 0, "evaluation"));
		bad(directory, tied, "long", matrix(new int[]{2}), sidecar(noneCorrect, 0, "evaluation"));
		bad(directory, tied, "onehot", matrix(new int[]{2}).replaceFirst("\\t1\\t", "\t0.5\t"), sidecar(new int[]{2}, 0, "evaluation"));
		bad(directory, tied, "no_target", matrix(new int[]{-1}), sidecar(new int[]{2}, 0, "evaluation"));
		fullErrors(directory);
	}

	/** Literal full-tail errors distinguish OOR answers from positional corrections. */
	private static void fullErrors(final File directory){
		final long[] corrections={0, 2, 3, -1, 20, -20};
		final CellNet moved=constant(69, 22);
		moved.net[1][12].setBias(1, true);
		for(int end=0; end<2; end++){
			final String matrix=write(directory, "full"+end+".tsv", fullMatrix(corrections));
			final String rows=write(directory, "full"+end+".rows.tsv", fullRows(corrections, end));
			final String raw=path(directory, "full"+end+".scores.tsv");
			final EvaluateNcrnaPositionalNet.FullComparison comparison;
			try(final NcrnaTextIO.Output output=new NcrnaTextIO.Output(raw)){
				comparison=EvaluateNcrnaPositionalNet.compareFull(moved, matrix, rows, end, output.writer);
			}
			final EvaluateNcrnaPositionalNet.FullResult baseline=comparison.doNothing, result=comparison.margin0;
			check(baseline.trials==6 && baseline.exact==1 && baseline.within1==2 && baseline.within3==4
				&& baseline.squaredBefore==814 && baseline.squaredAfter==814,
				"Literal corrections0,2,3,-1,20,-20 have full-tail squared error814, not an index21 surrogate");
			check(result.exact==1 && result.within1==2 && result.within3==4 && result.squaredAfter==822
				&& result.movedCorrect==1 && result.reachableWrong==3 && result.fixedReachableWrong==1
				&& result.outOfReach==2 && result.heldOutOfReach==0 && result.movedOutOfReach==2,
				"A +2 move gives errors2,0,-1,3,-18,22 and must retain distant harms in all-row RMS");
			check(comparison.margin01.squaredAfter==822 && result.reachableOorWins==0 && result.reachableOorVetoOnly==0,
				"A unit-winning reachable offset clears both declared margins without an OOR prediction");
			try(final NcrnaTextIO.Input scores=new NcrnaTextIO.Input(raw);
					final NcrnaTextIO.Input identities=new NcrnaTextIO.Input(rows)){
				scores.reader.nextLine();
				identities.reader.nextLine();
				final LineParser1 fields=new LineParser1('\t');
				for(int i=0; i<corrections.length; i++){
					final byte[] record=scores.reader.nextLine(), identity=identities.reader.nextLine();
					fields.set(record);
					check(fields.terms()==38 && fields.parseLong(12)==corrections[i]
						&& fields.parseInt(13)==12 && fields.parseInt(14)==12 && fields.parseInt(15)==12,
						"All22outputs accompany real signed correction and actual policy selections");
					check(Arrays.equals(identity, Arrays.copyOf(record, identity.length)), "Full raw scores retain exact source identities");
					for(int output=0; output<22; output++){
						check(fields.parseInt(16+output)==(output==12 ? 1065353216 : 0), "Every one of22FP32 outputs is lossless");
					}
				}
				check(new String(scores.reader.nextLine(), java.nio.charset.StandardCharsets.US_ASCII).equals("#COMPLETE\t6")
					&& scores.reader.nextLine()==null && identities.reader.nextLine()==null, "Full score footer requires exact complete pairing");
			}

			final CellNet tied=constant(69, 22);
			tied.net[1][12].setBias(1, true);
			tied.net[1][13].setBias(1, true);
			check(EvaluateNcrnaPositionalNet.compareFull(tied, matrix, rows, end, null).margin0.squaredAfter==814,
				"Unlike the preserved21path, a best-offset tie holds even at full-mode margin0");
			final CellNet oor=constant(69, 22);
			oor.net[1][12].setBias(1, true);
			oor.net[1][21].setBias(2, true);
			final EvaluateNcrnaPositionalNet.FullComparison wins=EvaluateNcrnaPositionalNet.compareFull(oor, matrix, rows, end, null);
			check(wins.margin0.squaredAfter==814 && wins.margin0.heldOutOfReach==2
				&& wins.margin0.reachableOorWins==4 && wins.margin0.reachableOorVetoOnly==0,
				"OOR wins hold the current end, expose false OOR on reachable truth and never move +11bases");
			oor.net[1][21].setBias(1, true);
			final EvaluateNcrnaPositionalNet.FullComparison equal=EvaluateNcrnaPositionalNet.compareFull(oor, matrix, rows, end, null);
			check(equal.margin0.squaredAfter==814 && equal.margin0.reachableOorWins==0 && equal.margin0.reachableOorVetoOnly==4,
				"OOR ties veto an otherwise valid offset but are separate from an OOR-head win");
			oor.net[1][12].setBias(.15f, true);
			oor.net[1][21].setBias(.1f, true);
			final EvaluateNcrnaPositionalNet.FullComparison margin=EvaluateNcrnaPositionalNet.compareFull(oor, matrix, rows, end, null);
			check(margin.margin0.squaredAfter==822 && margin.margin01.squaredAfter==814
				&& margin.margin0.reachableOorVetoOnly==0 && margin.margin01.reachableOorVetoOnly==4,
				"Positive margin independently gates against OOR, beyond current and runner-up");
			oor.net[1][12].setBias(.1f, true);
			oor.net[1][21].setBias(-1, true);
			final EvaluateNcrnaPositionalNet.FullComparison boundary=EvaluateNcrnaPositionalNet.compareFull(oor, matrix, rows, end, null);
			check(boundary.margin0.squaredAfter==822 && boundary.margin01.squaredAfter==814
				&& boundary.margin01.reachableOorWins==0 && boundary.margin01.reachableOorVetoOnly==0,
				"Equality at the declared margin is an ordinary confidence hold, not an OOR classification");
			final int endpoint=end;
			final String wrong=write(directory, "full_wrong"+end+".rows.tsv", fullRows(new long[]{0, 2, 3, -1, 2, -20}, end));
			fails(()->EvaluateNcrnaPositionalNet.compareFull(moved, matrix, wrong, endpoint, null));
			fails(()->EvaluateNcrnaPositionalNet.compareFull(constant(), matrix, rows, endpoint, null));
		}
		final CellNet normalized=constant(69, 22);
		normalized.net[1][10].setBias(5, true);
		normalized.net[1][11].weights[0]=1;
		final float[] mean=new float[69], inverse=new float[69];
		Arrays.fill(inverse, 1);
		mean[0]=1;
		inverse[0]=2;
		normalized.setInputNormalization(mean, inverse);
		final String net=write(directory, "full_normalized.bbnet", normalized.toBytes().toString());
		final String matrix=write(directory, "full_correct.tsv", fullMatrix(new long[]{0}));
		final String rows=write(directory, "full_correct.rows.tsv", fullRows(new long[]{0}, 0));
		check(EvaluateNcrnaPositionalNet.compareFull(CellNetParser.load(net, false), matrix, rows, 0, null).margin0.exact==1,
			"Actual22-output export/load must apply serialized raw-input normalization exactly once");
		final String output=path(directory, "full_cli.tsv");
		EvaluateNcrnaPositionalNet.main(new String[]{"in="+matrix, "rows="+rows, "net="+net, "end=5prime", "out="+output, "fullerrors=t"});
		check(new File(output).isFile() && !new File(output+".comparison.tsv").exists(),
			"Full mode emits one complete policy table, without a redundant comparison copy");
		final String wrongHeader=write(directory, "full_bad_header.tsv", fullMatrix(new long[]{0}).replace("69\t22\t0", "69\t21\t0"));
		fails(()->EvaluateNcrnaPositionalNet.compareFull(normalized, wrongHeader, rows, 0, null));
		final String training=write(directory, "full_training.rows.tsv", fullRows(new long[]{0}, 0).replace("evaluation", "train"));
		fails(()->EvaluateNcrnaPositionalNet.compareFull(normalized, matrix, training, 0, null));
		final String nonfinite=write(directory, "full_nonfinite.tsv", fullMatrix(new long[]{0}).replaceFirst("3\\t", "NaN\t"));
		fails(()->EvaluateNcrnaPositionalNet.compareFull(normalized, nonfinite, rows, 0, null));
		final String fractional=write(directory, "full_fractional.tsv", fullMatrix(new long[]{0}).replaceFirst("\\t1\\t", "\t0.5\t"));
		fails(()->EvaluateNcrnaPositionalNet.compareFull(normalized, fractional, rows, 0, null));
		final String distant=write(directory, "full_distant.tsv", fullMatrix(new long[]{20, -20}));
		final String distantRows=write(directory, "full_distant.rows.tsv", fullRows(new long[]{20, -20}, 0));
		final String[] record=EvaluateNcrnaPositionalNet.compareFull(normalized, distant, distantRows, 0, null)
			.summary(0).toString().split("\n")[1].split("\t");
		check(record[11].equals("NA") && record[14].equals("NA") && record[22].equals("NA") && record[24].equals("NA"),
			"Absent correct/reachable cohorts have unavailable percentages, not manufactured zero-percent performance");
	}

	/** Explicit22-unit targets preserve actual tail distances in the paired sidecar. */
	private static String fullMatrix(final long[] corrections){
		assert(corrections!=null && corrections.length>0) : "The literal matrix oracle needs a nonempty declared trial list";
		final ByteBuilder text=new ByteBuilder("#dims\t69\t22\t0\n");
		for(final long correction:corrections){
			text.append(3);
			for(int i=1; i<69; i++){text.tab().append(0);}
			final int target=Math.abs(correction)>10 ? 21 : (int)correction+10;
			for(int i=0; i<22; i++){text.tab().append(i==target ? 1 : 0);}
			text.nl();
		}
		return text.toString();
	}

	/** Both ends are valid inclusive spans; tails retain their actual signed distance. */
	private static String fullRows(final long[] corrections, final int end){
		assert(corrections!=null && corrections.length>0 && (end==0 || end==1)) : "Literal sidecars bind every matrix trial to one declared end";
		final ByteBuilder text=new ByteBuilder("locus_id\tsource\tassembly\tsplit\tmodel_index\tmodel\tsample\traw_start0\traw_stop0\ttrue_start0\ttrue_stop0\n");
		for(int i=0; i<corrections.length; i++){
			text.append("full_gene").append(i).append("\tsource\tassembly\tevaluation\t0\tmodel\t0\t")
				.append(end==0 ? 1000-corrections[i] : 1000).tab().append(end==1 ? 1999-corrections[i] : 1999)
				.append("\t1000\t1999\n");
		}
		return text.toString();
	}

	private static String path(final File directory, final String name){return new File(directory, name).getPath();}

	/** Compiled float literals independently bind the legal scientific-notation cases. */
	private static void exponents(final File directory){
		final String[] tokens={"3E1", "3e+1", "-2.3197112E-5"};
		final float[] expected={30, 30, -2.3197112e-5f};
		for(final int width:new int[]{24, 69}){
			final CellNet net=constant(width);
			for(int i=0; i<21; i++){net.net[1][i].setBias(i==10 ? 0 : -10, true);}
			net.net[1][10].weights[0]=1;
			final String rows=write(directory, "exponent"+width+".rows.tsv", sidecar(new int[]{10}, 0, "evaluation"));
			for(int i=0; i<tokens.length; i++){
				final String text=matrix(new int[]{10}, width).replaceFirst("3\\t", tokens[i]+"\t");
				final String input=write(directory, "exponent"+width+"_"+i+".tsv", text);
				check(EvaluateNcrnaPositionalNet.evaluate(net, input, rows, 0).exact==1
					&& Float.floatToRawIntBits(net.getOutput(10))==Float.floatToRawIntBits(expected[i]),
					"Exponent features must preserve the independently compiled FP32 value, including the real failing small negative token");
			}
			final String overflow=write(directory, "exponent_overflow"+width+".tsv",
				matrix(new int[]{10}, width).replaceFirst("3\\t", "1E99\t"));
			fails(()->EvaluateNcrnaPositionalNet.evaluate(net, overflow, rows, 0));
		}
	}

	/** A tied best score must hold at positive margin, with lossless saved scores. */
	private static void compare(final File directory, final int[] truth){
		for(final int width:new int[]{24, 69}){
			final CellNet net=constant(width);
			net.net[1][2].setBias(1, true);
			net.net[1][3].setBias(1, true);
			write(directory, "comparison"+width+".bbnet", net.toBytes().toString());
			for(int end=0; end<2; end++){
				final String stem="comparison"+width+"_"+end;
				final String matrix=write(directory, stem+".tsv", matrix(truth, width));
				final String rows=write(directory, stem+".rows.tsv", sidecar(truth, end, "evaluation"));
				final String raw=new File(directory, stem+".scores.tsv").getPath();
				final EvaluateNcrnaPositionalNet.Comparison result;
				try(final NcrnaTextIO.Output output=new NcrnaTextIO.Output(raw)){
					result=EvaluateNcrnaPositionalNet.compare(net, matrix, rows, end, output.writer);
				}
				check(result.ungated.exact==1 && result.ungated.squaredError==74,
					"Both widths must preserve first-maximum predictions and the original squared-error oracle");
				check(result.gated.exact==1 && result.gated.within1==1 && result.gated.within3==1
					&& result.gated.errorSum==20 && result.gated.squaredError==138 && result.gated.movedAlreadyCorrect==0,
					"Tied scores retain class10: literal errors8,7,5,0 have squared sum138 and move no correct end");
				check(result.doNothing.summary(end).toString().equals(result.gated.summary(end).toString())
					&& result.fixesKept==0 && result.fixesLost==1,
					"The tied net loses its one ungated correction; gating exactly equals doing nothing");
				try(final NcrnaTextIO.Input input=new NcrnaTextIO.Input(raw);
						final NcrnaTextIO.Input identities=new NcrnaTextIO.Input(rows)){
					final LineParser1 parser=new LineParser1('\t');
					check(input.reader.nextLine()!=null && identities.reader.nextLine()!=null, "Both saved streams have headers");
					for(int i=0; i<truth.length; i++){
						final byte[] line=input.reader.nextLine(), identity=identities.reader.nextLine();
						parser.set(line);
						check(parser.terms()==35 && parser.parseInt(11)==truth[i] && parser.parseInt(12)==2
							&& parser.parseInt(13)==10, "Raw record preserves truth and the two policy selections");
						check(Arrays.equals(identity, Arrays.copyOf(line, identity.length)) && line[identity.length]=='\t',
							"Each output starts with the exact original identity row, including raw/truth coordinates");
						for(int j=0; j<21; j++){
							check(parser.parseInt(14+j)==(j==2 || j==3 ? 1065353216 : 0),
								"Every output is stored losslessly: IEEE754 one is1065353216 and zero is0");
						}
					}
					check(new String(input.reader.nextLine(), java.nio.charset.StandardCharsets.US_ASCII).equals("#COMPLETE\t4")
						&& input.reader.nextLine()==null, "Only a complete four-trial raw stream has the terminal footer");
				}
			}
		}
	}

	/** Literal k-mer locations distinguish endpoint-base and word-start coordinates. */
	private static void project(final File directory){
		final byte[] bases=new byte[80];
		Arrays.fill(bases, (byte)'A');
		bases[11]=bases[42]='C';
		final String contexts=write(directory, "project.fa", ">g source=s lflank=20 rflank=20 contig_gc=0.4\n"
			+new String(bases, java.nio.charset.StandardCharsets.US_ASCII)+"\n");
		final String ref=write(directory, "project.ref.fa", ">m0\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n");
		final String roster=write(directory, "project.genomes.tsv",
			"index\tgenome\tassembly\ttaxid\tdomain\tphylum\tgff\tfasta\tloci\tgenome_bases\tcontigs\n"
			+"0\ts\ta\t.\t.\t.\tinput.gff\tinput.fa\t1\t80\t1\n");
		final ByteBuilder matrix=new ByteBuilder("#dims\t69\t21\t0\n");
		for(final int label:new int[]{10, 9}){
			matrix.append(0);
			for(int i=1; i<63; i++){matrix.append("\t0");}
			matrix.append("\t1.0000\t0.4000\t0.7500\t0\t0\t0");
			for(int i=0; i<21; i++){matrix.append(i==label ? "\t1.0000" : "\t0.0000");}
			matrix.nl();
		}
		final String input=write(directory, "project69.tsv", matrix.toString());
		for(int end=0; end<2; end++){
			final RrnaPositionalKmerTable.Table table=new RrnaPositionalKmerTable.Table("RF00174",
				end==0 ? "5prime" : "3prime", 9, end==0 ? 0 : -8, 1.0/262144, 10);
			table.add(bases, end==0 ? 20 : 59);
			table.finish();
			// Analytic lookup oracle: all-A=5; C first=7; C second=3.
			Arrays.fill(table.scores, 0);
			table.scores[0]=5;
			table.scores[65536]=7;
			table.scores[16384]=3;
			final String rowsText="locus_id\tsource\tassembly\tsplit\tmodel_index\tmodel\tsample"
				+"\traw_start0\traw_stop0\ttrue_start0\ttrue_stop0\n"
				+"g\ts\ta\tevaluation\t0\tm0\t0\t20\t59\t20\t59\n"
				+"g\ts\ta\tevaluation\t0\tm0\t1\t"+(end==0 ? "21\t59" : "20\t60")+"\t20\t59\n";
			final String rows=write(directory, "project"+end+".rows.tsv", rowsText);
			final String output=new File(directory, "project24_"+end+".tsv").getPath();
			check(ProjectNcrnaPositionalInputs.project(input, rows, contexts, roster, ref, table, end, 2, output)==2,
				"Projection conserves both original samples, permitting missing taxonomy");
			try(final NcrnaTextIO.Input projected=new NcrnaTextIO.Input(output);
					final NcrnaTextIO.Input original=new NcrnaTextIO.Input(input)){
				check(new String(projected.reader.nextLine(), java.nio.charset.StandardCharsets.US_ASCII).equals("#dims\t24\t21\t0"),
					"Projection declares exactly24features and the unchanged21classes");
				original.reader.nextLine();
				final LineParser1 values=new LineParser1('\t'), source=new LineParser1('\t');
				for(int row=0; row<2; row++){
					final byte[] line=projected.reader.nextLine(), before=original.reader.nextLine();
					values.set(line);
					source.set(before);
					check(values.terms()==45, "Each projected row has24features plus21targets");
					for(int i=0; i<21; i++){
						final int expected=row==0 && i==0 ? 3 : (i==(row==0 ? 1 : 0) ? 7 : 5);
						check(values.parseFloat(i)==expected,
							"Inside-K9 starts are10/11at5prime and41/42at3prime; jitter shifts the full score window once");
					}
					values.setBounds(21);
					source.setBounds(63);
					final int commonStart=values.a(), originalStart=source.a();
					source.setBounds(66);
					final int commonLength=source.a()-originalStart-1;
					check(Arrays.equals(Arrays.copyOfRange(line, commonStart, commonStart+commonLength),
						Arrays.copyOfRange(before, originalStart, originalStart+commonLength)),
						"Length/GC/identity tokens retain the producer's exact spelling");
					values.setBounds(24);
					source.setBounds(69);
					check(Arrays.equals(Arrays.copyOfRange(line, values.a(), line.length),
						Arrays.copyOfRange(before, source.a(), before.length)), "One-hot label tokens are copied verbatim");
				}
				check(projected.reader.nextLine()==null, "No duplicated or synthesized trial is emitted");
			}
			final int endpoint=end;
			final String wrong=write(directory, "project_wrong"+end+".rows.tsv", rowsText.replace("\tm0\t", "\twrong\t"));
			fails(()->ProjectNcrnaPositionalInputs.project(input, wrong, contexts, roster, ref, table, endpoint, 2,
				new File(directory, "wrong_project"+endpoint+".tsv").getPath()));
			fails(()->ProjectNcrnaPositionalInputs.project(input, rows, contexts, roster, ref, table, endpoint, 3,
				new File(directory, "missing_project"+endpoint+".tsv").getPath()));
		}
	}

	/** Parses the single small summary record for semantic numeric/NA checks. */
	private static LineParser1 summary(final Result result, final int end){
		final String row=result.summary(end).toString().split("\n")[1];
		return new LineParser1('\t').set(row.getBytes(java.nio.charset.StandardCharsets.US_ASCII));
	}

	/** A real69→21 dense linear network, with all incoming weights and biases zero. */
	private static CellNet constant(){
		return constant(69);
	}

	/** Uses the same explicit dense fixture format for either supported input width. */
	private static CellNet constant(final int width){
		return constant(width, 21);
	}

	/** Uses the same real dense-network parser for21- and22-answer finite cases. */
	private static CellNet constant(final int width, final int outputs){
		final ByteBuilder text=new ByteBuilder("##bbnet\n#version 1\n#concise\n#dense\n#density 1\n"
			+"#blocksize 1\n#seed 58\n#layers 2\n#dims ").append(width).space().append(outputs).append("\n#coding decimal\n");
		for(int output=width+1; output<=width+outputs; output++){
			text.append('C').append(output).append(" LINEAR 0");
			for(int input=0; input<width; input++){text.append(" 0");}
			text.nl();
		}
		final ArrayList<byte[]> lines=new ArrayList<byte[]>();
		for(final String line:text.toString().split("\n")){
			lines.add(line.getBytes(java.nio.charset.StandardCharsets.US_ASCII));
		}
		return CellNetParser.loadFromLines(lines);
	}

	/** Saves deliberately malformed inputs and requires a reader failure. */
	private static void bad(final File directory, final CellNet net, final String name, final String matrix, final String rows){
		final String data=write(directory, name+".tsv", matrix);
		final String metadata=write(directory, name+".rows.tsv", rows);
		fails(()->EvaluateNcrnaPositionalNet.evaluate(net, data, metadata, 0));
	}

	/** Finite raw vectors have first feature3 and otherwise zero features. */
	private static String matrix(final int[] truth){
		return matrix(truth, 69);
	}

	/** The schema width, features and targets agree independently of the net. */
	private static String matrix(final int[] truth, final int width){
		final ByteBuilder buffer=new ByteBuilder("#dims\t").append(width).append("\t21\t0\n");
		for(final int label:truth){
			buffer.append(3);
			for(int i=1; i<width; i++){buffer.tab().append(0);}
			for(int i=0; i<21; i++){buffer.tab().append(i==label ? 1 : 0);}
			buffer.nl();
		}
		return buffer.toString();
	}

	/** Explicit truth100..199, with the chosen raw end displaced by10-label. */
	private static String sidecar(final int[] truth, final int end, final String split){
		final ByteBuilder buffer=new ByteBuilder("locus_id\tsource\tassembly\tsplit\tmodel_index\tmodel\tsample"
			+"\traw_start0\traw_stop0\ttrue_start0\ttrue_stop0\n");
		for(int i=0; i<truth.length; i++){
			buffer.append("gene").append(i).append("\tsource\tassembly\t").append(split).append("\t0\tmodel\t0\t")
				.append(end==0 ? 110-truth[i] : 100).tab().append(end==1 ? 209-truth[i] : 199)
				.append("\t100\t199\n");
		}
		return buffer.toString();
	}

	/** Writes only small fresh fixture files with the native owner. */
	private static String write(final File root, final String name, final String text){
		final String path=new File(root, name).getPath();
		try(final NcrnaTextIO.Output output=new NcrnaTextIO.Output(path)){
			output.writer.print(text.getBytes(java.nio.charset.StandardCharsets.US_ASCII));
		}
		return path;
	}

	private static void fails(final Runnable action){
		boolean failed=false;
		try{action.run();}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Malformed schema, pairing, provenance or output must fail loudly");
	}

	/** Accepts only the native NaN guard or the scorer's later nonfinite rejection. */
	private static void failsNonfinite(final Runnable action){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		catch(AssertionError expected){
			final StackTraceElement[] trace=expected.getStackTrace();
			if(!"NaN, NaN, NaN".equals(expected.getMessage()) || trace.length==0
					|| !"ml.Cell".equals(trace[0].getClassName()) || !"setValue".equals(trace[0].getMethodName())){
				throw expected;
			}
			return;
		}
		throw new AssertionError("A nonfinite network output must be rejected");
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
}
