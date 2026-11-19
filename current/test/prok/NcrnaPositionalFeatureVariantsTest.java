package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Locale;

import ml.CellNet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaPositionalFeatureVariants;
import prok.EvaluateNcrnaPositionalNet;
import prok.EvaluateNcrnaPositionalNet.FullComparison;
import prok.NcrnaTextIO.Input;
import shared.Shared;
import structures.ByteBuilder;

/** Literal transforms plus the six-matrix exact-identity join and failure cases.
 * Uses the existing original-population fixture, without training/alignment.
 * @author Brian Bushnell, Ganyu
 */
public final class NcrnaPositionalFeatureVariantsTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaPositionalFeatureVariantsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, pair[0].toLowerCase(Locale.ROOT), pair.length>1 ? pair[1] : null), "Unknown fixture option");
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.maxReads<0
				&& parser.samplerate==1 && !parser.overwrite && !parser.append && new File(parser.out1).mkdir(), "Supply out=fresh_fixture_directory");
			run(new File(parser.out1));
			pp.outstream.println("NCRNA_POSITIONAL_FEATURE_VARIANTS_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final File root){
		check(NcrnaPositionalFeatureVariants.referenceLength(50, 50, 100)==0
			&& NcrnaPositionalFeatureVariants.referenceLength(75, 50, 100)==.5f
			&& NcrnaPositionalFeatureVariants.referenceLength(100, 50, 100)==1
			&& NcrnaPositionalFeatureVariants.referenceLength(50, 50, 50)==0,
			"Frozen reference endpoints, midpoint and equal-length range have literal scaling");
		check(NcrnaPositionalFeatureVariants.gc("ANNcGTN".getBytes(StandardCharsets.US_ASCII), 0, 6)==.5f,
			"GC counts recognized bases, not ambiguous symbols, and accepts lowercase");
		fails(()->NcrnaPositionalFeatureVariants.gc(new byte[]{'N'}, 0, 0));
		fails(()->NcrnaPositionalFeatureVariants.referenceLength(49, 50, 100));
		final String roster=path(root, "roster.tsv");
		final String loci=path(root, "loci.tsv");
		final String graph=path(root, "graph.tsv");
		final String exclude=path(root, "exclude.txt");
		final String contexts=path(root, "contexts.fa");
		final String original=path(root, "original.fa");
		final String ref=path(root, "ref.fa");
		write(roster, NcrnaLearningPopulationTest.roster());
		write(loci, NcrnaLearningPopulationTest.loci());
		write(graph, NcrnaLearningPopulationTest.graph());
		write(exclude, "gene2\n");
		final ByteBuilder contextText=new ByteBuilder();
		for(int gene=0; gene<=20; gene++){
			contextText.append(">gene").append(gene).append(" source=genome").append(gene==20 ? 0 : gene)
				.append(" lflank=20 rflank=30 contig_gc=0.5\n");
			for(int i=0; i<150; i++){contextText.append("ACGT".charAt(i&3));}
			contextText.nl();
		}
		write(contexts, contextText.toString());
		final ByteBuilder models=new ByteBuilder(">common\n");
		for(int i=0; i<50; i++){models.append('A');}
		models.append("\n>rare\n");
		for(int i=0; i<100; i++){models.append('G');}
		models.nl();
		write(ref, models.toString());
		write(original, models.toString());
		final File input=new File(root, "input");
		check(input.mkdir(), "Fresh input fixture");
		for(int set=0; set<3; set++){
			for(int end=0; end<2; end++){
				final String stem=stem(set, end);
				write(path(input, stem+".tsv"), matrix(set, end));
				write(path(input, stem+".rows.tsv"), NcrnaPositionalFeatureVariants.ROW_HEADER+'\n'+row(set));
			}
		}
		final String[] args={"in="+input, "out="+path(root, "good"), "contexts="+contexts,
			"ref="+ref, "original="+original, "graph="+graph, "roster="+roster, "loci="+loci,
			"exclude="+exclude, "samples=1", "expected=1,1,1,1,1,1"};
		new NcrnaPositionalFeatureVariants(args).process();
		for(int variant=0; variant<3; variant++){
			final String name=new String[]{"length", "gc", "identity"}[variant];
			final int width=new int[]{70, 71, 69}[variant];
			final String matrix=path(new File(root, "good/"+name), "r10.evaluation.5prime.tsv.gz");
			final String rows=path(input, "r10.evaluation.5prime.rows.tsv");
			final FullComparison result=EvaluateNcrnaPositionalNet.compareFull(new FixedNet(width), matrix, rows, 0, null, width);
			check(result.margin0.trials==1 && result.margin0.movedCorrect==1 && result.margin0.worse==1,
				"Each explicit feature width uses the same endpoint policy and signed-distance accounting");
			if(width!=69){fails(()->EvaluateNcrnaPositionalNet.compareFull(new FixedNet(width), matrix, rows, 0, null));}
		}
		replay(root);
		for(int set=0; set<3; set++){
			for(int end=0; end<2; end++){
				final String stem=stem(set, end);
				final LineParser1 before=new LineParser1('\t');
				before.set(matrix(set, end).split("\n")[1].getBytes(StandardCharsets.US_ASCII));
				for(int variant=0; variant<3; variant++){
					final String name=new String[]{"length", "gc", "identity"}[variant];
					final int width=new int[]{70, 71, 69}[variant];
					try(Input in=new Input(path(new File(root, "good/"+name), stem+".tsv.gz"))){
						check(new String(in.reader.nextLine(), StandardCharsets.US_ASCII).equals("#dims\t"+width+"\t22\t0"), "Projected width must match its schema");
						final LineParser1 after=new LineParser1('\t');
						after.set(in.reader.nextLine());
						check(after.terms()==width+22 && in.reader.nextLine()==null, "Exactly one preserved trial in fixture output");
						for(int i=0; i<63; i++){same(before, i, after, i);}
						for(int i=0; i<22; i++){same(before, 69+i, after, width+i);}
						if(variant==0){
							check(after.parseFloat(63)==0 && after.parseFloat(64)==(set==2 ? 1 : 1.5f), "Raw candidate length remains unclamped above one");
							for(int i=64; i<69; i++){same(before, i, after, i+1);}
						}else if(variant==1){
							same(before, 63, after, 63);
							same(before, 64, after, 64);
							check(after.parseFloat(65)==.5f && after.parseFloat(66)==0, "Span GC and reference GC are different independently reconstructed inputs");
							for(int i=65; i<69; i++){same(before, i, after, i+2);}
						}else{
							for(int i=63; i<69; i++){if(i!=65){same(before, i, after, i);}}
							check(after.parseFloat(65)==-.5f, "Identity below the caller floor is retained without clipping on jittered trials");
						}
					}
				}
			}
		}
		try(Input in=new Input(path(new File(root, "good"), "projection.tsv"))){
			in.reader.nextLine();
			int records=0;
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				final LineParser1 p=new LineParser1('\t');
				p.set(line);
				check(p.terms()==13 && p.parseLong(3)==1 && p.parseInt(6)==50 && p.parseInt(7)==100
					&& p.parseLong(8)==(p.termEquals("evaluation", 1) ? 0 : 1) && p.parseLong(9)==1, "Summary binds row counts, range and out-of-unit inputs");
				records++;
			}
			check(records==18, "Every variant/split/end has one summary row");
		}
		final String badRow=path(input, stem(0, 0)+".rows.tsv");
		for(int kind=0; kind<4; kind++){
			String changed=row(0);
			if(kind==0){changed=changed.replace("genome1", "genome2");}
			if(kind==1){changed=changed.replace("\tcommon\t", "\trare\t");}
			if(kind==2){changed=changed.replace("\t20\t119\n", "\t21\t119\n");}
			if(kind==3){changed+=row(0);}
			write(badRow, NcrnaPositionalFeatureVariants.ROW_HEADER+'\n'+changed);
			args[1]="out="+path(root, "bad"+kind);
			fails(()->new NcrnaPositionalFeatureVariants(args).process());
		}
		write(badRow, NcrnaPositionalFeatureVariants.ROW_HEADER+'\n'+row(0));
		final String badMatrix=path(input, stem(0, 0)+".tsv");
		write(badMatrix, matrix(0, 0).replace("2.5E-7", "2.5Eoops"));
		args[1]="out="+path(root, "bad_exponent");
		fails(()->new NcrnaPositionalFeatureVariants(args).process());
		write(badMatrix, matrix(2, 0));
		args[1]="out="+path(root, "bad_target");
		fails(()->new NcrnaPositionalFeatureVariants(args).process());
		final String duplicate=matrix(0, 0);
		write(badMatrix, duplicate+duplicate.substring(duplicate.indexOf('\n')+1));
		write(badRow, NcrnaPositionalFeatureVariants.ROW_HEADER+'\n'+row(0)+row(0));
		args[1]="out="+path(root, "duplicate_trial");
		fails(()->new NcrnaPositionalFeatureVariants(args).process());
	}

	/** Hand-derived six-trial counterexample: improvement ratio .5 while RMS worsens. */
	private static void replay(final File root){
		final ByteBuilder raw=new ByteBuilder("#NCRNA_FULL_ERROR_SCORES_V1\n#coding\tFP32_RAW_SIGNED_INT\n"
			+"#matrix\t/nonexistent/matrix\n#rows\t/nonexistent/rows\n#network\t/nonexistent/net\n#end\t5prime\n")
			.append(NcrnaPositionalFeatureVariants.ROW_HEADER)
			.append("\ttrue_class\ttrue_minus_raw_bp\tmargin0_class\tmargin01_class\targmax22_class");
		for(int i=0; i<22; i++){raw.append("\tscore_bits_").append(i);}
		raw.nl();
		final int[] correction={2, 0, 20, -20, 3, 5};
		final int[] chosen={12, 12, 20, 20, 10, 10};
		for(int trial=0; trial<6; trial++){
			final int truth=Math.abs(correction[trial])>10 ? 21 : correction[trial]+10;
			raw.append("gene").append(trial).append("\tsource\tassembly\tevaluation\t0\tcommon\t0\t")
				.append(100-correction[trial]).append("\t250\t100\t250\t").append(truth).tab().append(correction[trial])
				.tab().append(chosen[trial]).tab().append(chosen[trial]).tab().append(trial==5 ? 21 : chosen[trial]);
			for(int i=0; i<22; i++){
				final float score=trial==5 ? (i==21 ? 1 : i==12 ? .8f : 0) : (i==chosen[trial] ? 1 : 0);
				raw.tab().append(Float.floatToRawIntBits(score));
			}
			raw.nl();
		}
		raw.append("#COMPLETE\t6\n");
		final String saved=path(root, "saved.scores.tsv");
		write(saved, raw.toString());
		final FullComparison result=EvaluateNcrnaPositionalNet.replayFull(saved, 0);
		check(result.margin0.trials==6 && result.margin0.better==2 && result.margin0.worse==2 && result.margin0.tied==2
			&& result.margin0.improvedOutOfReach==1 && result.margin0.worsenedOutOfReach==1
			&& result.margin0.squaredBefore==838 && result.margin0.squaredAfter==1038
			&& result.margin0.exact==1 && result.doNothing.exact==1 && result.doNothing.tied==6,
			"Literal distances, not OOR labels, determine improvement/harm and full-tail RMS");
		final String summary=path(root, "replayed.tsv");
		EvaluateNcrnaPositionalNet.main(new String[]{"fullerrors=t", "replayfull=t", "extended=t",
			"in="+saved, "end=5prime", "out="+summary});
		try(Input input=new Input(summary)){
			final LineParser1 p=new LineParser1('\t');
			p.set(input.reader.nextLine());
			check(p.terms()==37 && p.termEquals("rms_change_bp", 27), "Extended output appends explicit paired-change fields");
			p.set(input.reader.nextLine());
			check(p.termEquals("NA", 31) && p.parseFloat(27)==0, "No changed distances means unavailable improvement ratio and zero RMS change");
			p.set(input.reader.nextLine());
			check(p.parseFloat(27)>0 && p.parseLong(28)==2 && p.parseLong(29)==2 && p.parseFloat(31)==.5f,
				"Equal improvement/worsening counts can still increase RMS; report both faithfully");
		}
		final String incomplete=path(root, "incomplete.scores.tsv");
		write(incomplete, raw.toString().replace("#COMPLETE\t6\n", ""));
		fails(()->EvaluateNcrnaPositionalNet.replayFull(incomplete, 0));
		final String extra=path(root, "extra.scores.tsv");
		write(extra, raw.toString()+"extra\n");
		fails(()->EvaluateNcrnaPositionalNet.replayFull(extra, 0));
		final String wrong=path(root, "wrong.scores.tsv");
		write(wrong, raw.toString().replaceFirst("\\t12\\t12\\t12\\t", "\t10\t12\t12\t"));
		fails(()->EvaluateNcrnaPositionalNet.replayFull(wrong, 0));
		fails(()->EvaluateNcrnaPositionalNet.replayFull(saved, 1));
	}

	private static final class FixedNet extends CellNet{
		FixedNet(final int width_){width=width_;scores[12]=1;}
		@Override
		public int numInputs(){return width;}
		@Override
		public int numOutputs(){return 22;}
		@Override
		public CellNet applyInput(final float[] input){
			NcrnaPositionalFeatureVariantsTest.check(input.length==width, "The explicit feature schema must reach the real network input");
			return this;
		}
		@Override
		public float feedForward(){return 1;}
		@Override
		public float getOutput(final int index){return scores[index];}
		private final int width;
		private final float[] scores=new float[22];
	}

	private static String matrix(final int set, final int end){
		final ByteBuilder b=new ByteBuilder("#dims\t69\t22\t0\n");
		for(int i=0; i<63; i++){
			if(i==0){b.append("-3.1e1\t");}
			else if(i==1){b.append("2.5E-7\t");}
			else{b.append(i-31).append(".0000\t");}
		}
		b.append(set==2 ? "2.0000" : "3.0000").append("\t5e-1\t2.5E-1\t0.1250\t0.5000\t0.8750\t");
		for(int i=0; i<22; i++){
			if(i>0){b.tab();}
			b.append(i==(set==2 ? 10 : 21) ? 1 : 0);
		}
		return b.nl().toString();
	}
	private static String row(final int set){
		final int gene=set==2 ? 2 : 1;
		return "gene"+gene+"\tgenome"+gene+"\tassembly"+gene+'\t'+SETS[set]+"\t0\tcommon\t0\t"
			+(set==2 ? "20\t119" : "0\t149")+"\t20\t119\n";
	}
	private static void same(final LineParser1 before, final int source, final LineParser1 after, final int target){
		check(before.parseString(source).equals(after.parseString(target)), "Unchanged feature/label bytes must be preserved exactly");
	}
	private static String stem(final int set, final int end){return "r10."+SETS[set]+'.'+(end==0 ? "5prime" : "3prime");}
	private static String path(final File root, final String file){return new File(root, file).getPath();}
	private static void write(final String path, final String text){NcrnaLearningPopulationTest.write(path, text);}
	private static void fails(final Runnable action){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Invalid retained row must fail loudly");
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private static final String[] SETS={"train", "monitor", "evaluation"};
}
