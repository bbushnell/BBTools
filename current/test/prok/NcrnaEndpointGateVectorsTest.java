package test.prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import ml.CellNet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaEndpointGateVectors;
import prok.NcrnaPositionalFeatureVariants;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Literal proposals and full population/score joins, without training a network.
 * @author Brian Bushnell, Ganyu
 */
public final class NcrnaEndpointGateVectorsTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaEndpointGateVectorsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, pair[0].toLowerCase(Locale.ROOT), pair.length>1 ? pair[1] : null), "Unknown gate fixture argument");
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Supply only out=fresh_directory and native runtime options");
			final File root=new File(parser.out1);
			check(root.mkdir(), "Fixture directory must be fresh");
			run(root);
			pp.outstream.println("NCRNA_ENDPOINT_GATE_VECTORS_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private static void run(final File root){
		final float[] scores=new float[22];
		check(NcrnaEndpointGateVectors.proposal(scores)==-1, "A complete reachable tie holds");
		scores[10]=1;
		check(NcrnaEndpointGateVectors.proposal(scores)==-1, "The unchanged endpoint winning bypasses the gate");
		scores[12]=2;
		scores[21]=3;
		check(NcrnaEndpointGateVectors.proposal(scores)==12, "An out-of-range head may be overridden by the learned gate");
		scores[11]=2;
		check(NcrnaEndpointGateVectors.proposal(scores)==-1, "A tied reachable proposal holds even with a higher non-positional head");
		scores[11]=Float.NaN;
		fails(()->NcrnaEndpointGateVectors.proposal(scores));
		scores[11]=0;
		final float[] input=new float[69];
		input[63]=2;
		input[64]=.4f;
		input[65]=.25f;
		input[66]=.1f;
		input[67]=.2f;
		input[68]=.3f;
		final float[] result=new float[32];
		NcrnaEndpointGateVectors.fill(scores, input, .6f, .25f, 98, 50, result);
		check(result[22]==2 && result[23]==1.96f && result[24]==.4f && result[25]==.6f && result[26]==.25f
			&& result[27]==.6f-.25f && result[28]==-.5f && result[29]==.1f && result[31]==.3f,
			"Original/proposed lengths and raw-span GC/identity/fuzz have distinct literal positions");
		for(int i=0; i<22; i++){check(result[i]==scores[i], "First-stage heads remain unmodified");}
		NcrnaEndpointGateVectors.fill(scores, input, .6f, .25f, -1, 50, result);
		check(result[23]==-.02f, "Crossed trial proposals remain explicit rather than silently clipped");

		final String roster=write(root, "roster.tsv", NcrnaLearningPopulationTest.roster());
		final String loci=write(root, "loci.tsv", NcrnaLearningPopulationTest.loci());
		final String graph=write(root, "graph.tsv", NcrnaLearningPopulationTest.graph());
		final String exclude=write(root, "exclude.txt", "gene2\n");
		final ByteBuilder model=new ByteBuilder(">common\n");
		for(int i=0; i<50; i++){model.append('A');}
		model.append("\n>rare\n");
		for(int i=0; i<100; i++){model.append('G');}
		model.nl();
		final String ref=write(root, "ref.fa", model.toString());
		final String original=write(root, "original.fa", model.toString());
		final ByteBuilder contexts=new ByteBuilder();
		final ByteBuilder mapping=new ByteBuilder("locus_id\tsource\tassembly\trole\tfold\tpartial\n");
		for(int gene=0; gene<=20; gene++){
			final int source=gene==20 ? 0 : gene;
			contexts.append(">gene").append(gene).append(" source=genome").append(source)
				.append(" lflank=20 rflank=30 contig_gc=0.5\n");
			for(int i=0; i<150; i++){contexts.append("ACGT".charAt(i&3));}
			contexts.nl();
			mapping.append("gene").append(gene).append("\tgenome").append(source).append("\tassembly").append(source)
				.append(gene==2 ? "\tOUTER_TEST\t-1\t" : "\tINTERNAL_FOLD\t0\t").append(gene==20 ? 1 : 0).nl();
		}
		final String contextPath=write(root, "contexts.fa", contexts.toString());
		final String map=write(root, "mapping.tsv", mapping.toString());
		final String net=write(root, "unused_net.txt", "No network is loaded for saved-score replay or a finite injected net.\n");
		final String heldRow="gene1\tgenome1\tassembly1\tFOLD_HELD\t0\tcommon\t0\t18\t119\t20\t119";
		final String outerRow=heldRow.replace("gene1\tgenome1\tassembly1\tFOLD_HELD", "gene2\tgenome2\tassembly2\tevaluation");
		final String numeric=matrix(102/50f);
		final String matrix=write(root, "matrix.tsv", numeric);
		final String rows=write(root, "rows.tsv", NcrnaPositionalFeatureVariants.ROW_HEADER+'\n'+heldRow+'\n');
		final String[] args={"in="+matrix, "rows="+rows, "net="+net, "contexts="+contextPath, "ref="+ref,
			"original="+original, "graph="+graph, "roster="+roster, "loci="+loci, "exclude="+exclude,
			"mapping="+map, "fold=0", "samples=1", "expected=1", "minlen=101", "end=5prime", "out="+path(root, "fold")};
		final FixedNet fixed=new FixedNet();
		new NcrnaEndpointGateVectors(args).process(fixed);
		check(fixed.calls==1, "Exactly one first-stage prediction per held row");
		try(Input vectors=new Input(path(root, "fold/vectors.tsv.gz"))){
			check(Arrays.equals(vectors.reader.nextLine(), "#dims\t32\t1\t0".getBytes()), "Gate width is32with one binary target");
			final LineParser1 p=new LineParser1('\t').set(vectors.reader.nextLine());
			check(p.terms()==33 && p.parseFloat(22)==102/50f && p.parseFloat(23)==2
				&& p.parseInt(32)==1 && vectors.reader.nextLine()==null,
				"A +2start correction is positive even when the101base assay minimum makes its100base proposal invalid");
		}
		try(Input summary=new Input(path(root, "fold/summary.tsv"))){
			summary.reader.nextLine();
			final LineParser1 p=new LineParser1('\t').set(summary.reader.nextLine());
			check(p.parseLong(2)==1 && p.parseLong(3)==1 && p.parseLong(4)==1 && p.parseLong(7)==1 && p.parseLong(8)==1,
				"Invalid improving proposals remain counted, with their original strict-improvement labels");
		}
		final String stopRow=heldRow.replace("\t18\t119\t20\t119", "\t20\t117\t20\t119");
		write(matrix, matrix(98/50f));
		write(rows, NcrnaPositionalFeatureVariants.ROW_HEADER+'\n'+stopRow+'\n');
		args[15]="end=3prime";
		args[16]="out="+path(root, "stop");
		new NcrnaEndpointGateVectors(args).process(new FixedNet());
		try(Input vectors=new Input(path(root, "stop/vectors.tsv.gz"))){
			vectors.reader.nextLine();
			final LineParser1 p=new LineParser1('\t').set(vectors.reader.nextLine());
			check(p.parseFloat(22)==98/50f && p.parseFloat(23)==2 && p.parseInt(32)==1,
				"A +2stop correction increases proposed length while a +2start correction decreases it");
		}
		write(matrix, numeric);
		args[15]="end=5prime";
		write(rows, NcrnaPositionalFeatureVariants.ROW_HEADER+'\n'+outerRow+'\n');
		final String raw=write(root, "scores.tsv", saved(matrix, rows, net, outerRow));
		final String[] replay=Arrays.copyOf(args, args.length+1);
		replay[11]="fold=-1";
		replay[16]="out="+path(root, "outer");
		replay[17]="scores="+raw;
		new NcrnaEndpointGateVectors(replay).process();
		check(!new File(path(root, "outer/first_stage.tsv.gz")).exists(), "Outer rows never repeat first-stage inference");
		write(raw, saved(matrix, rows, net, outerRow).replace("\t10\t10\t21\t", "\t10\t10\t12\t")
			.replace("\t"+Float.floatToRawIntBits(2f)+"\n", "\t"+Float.floatToRawIntBits(1f)+"\n"));
		replay[16]="out="+path(root, "outer_head_tie");
		new NcrnaEndpointGateVectors(replay).process();
		check(new File(path(root, "outer_head_tie/summary.tsv")).exists(),
			"A tied non-positional head retains the first reachable maximum in the saved22argmax field");
		args[16]="out="+path(root, "leak");
		fails(()->new NcrnaEndpointGateVectors(args).process(new FixedNet()));
		write(rows, NcrnaPositionalFeatureVariants.ROW_HEADER+'\n'+heldRow+'\n');
		args[13]="expected=2";
		args[16]="out="+path(root, "truncated");
		fails(()->new NcrnaEndpointGateVectors(args).process(new FixedNet()));
		check(!new File(path(root, "truncated/summary.tsv")).exists(), "A truncated paired input cannot publish completion");
		args[13]="expected=1";
		write(rows, NcrnaPositionalFeatureVariants.ROW_HEADER+'\n'+outerRow+'\n');
		write(raw, saved(matrix, rows, net, outerRow).replace("#COMPLETE\t1\n", ""));
		replay[16]="out="+path(root, "missing_terminal");
		fails(()->new NcrnaEndpointGateVectors(replay).process());
		write(raw, saved(matrix, rows, net, outerRow.replace("\tcommon\t", "\twrong\t")));
		replay[16]="out="+path(root, "wrong_join");
		fails(()->new NcrnaEndpointGateVectors(replay).process());
	}
	private static String matrix(final float ratio){
		final ByteBuilder text=new ByteBuilder("#dims\t69\t22\t0\n");
		for(int i=0; i<69; i++){
			if(i==0){text.append("2.5E-7");}
			else{text.append(i==63 ? ratio : i==64 ? .5f : i==65 ? .25f : 0, 18);}
			text.tab();
		}
		for(int i=0; i<22; i++){
			text.append(i==12 ? 1 : 0);
			if(i<21){text.tab();}
		}
		return text.nl().toString();
	}
	private static String saved(final String matrix, final String rows, final String net, final String identity){
		final ByteBuilder text=new ByteBuilder("#NCRNA_FULL_ERROR_SCORES_V1\n#coding\tFP32_RAW_SIGNED_INT\n")
			.append("#matrix\t").append(matrix).nl().append("#rows\t").append(rows).nl().append("#network\t").append(net)
			.append("\n#end\t5prime\n").append(NcrnaPositionalFeatureVariants.ROW_HEADER)
			.append("\ttrue_class\ttrue_minus_raw_bp\tmargin0_class\tmargin01_class\targmax22_class");
		for(int i=0; i<22; i++){text.append("\tscore_bits_").append(i);}
		text.nl().append(identity).append("\t12\t2\t10\t10\t21");
		for(int i=0; i<22; i++){text.tab().append(Float.floatToRawIntBits(i==12 ? 1 : i==21 ? 2 : 0));}
		return text.append("\n#COMPLETE\t1\n").toString();
	}
	private static final class FixedNet extends CellNet{
		@Override public int numInputs(){return 69;}
		@Override public int numOutputs(){return 22;}
		@Override
		public CellNet applyInput(final float[] input){
			NcrnaEndpointGateVectorsTest.check(input.length==69, "Frozen first-stage width");
			return this;
		}
		@Override public float feedForward(){calls++;return 1;}
		@Override public float getOutput(final int index){return index==12 ? 1 : index==21 ? 2 : 0;}
		int calls;
	}
	private static String path(final File root, final String name){return new File(root, name).getPath();}
	private static String write(final File root, final String name, final String text){
		final String file=path(root, name);
		write(file, text);
		return file;
	}
	private static void write(final String path, final String text){
		try(Output output=new Output(new ByteStreamWriter(path, true, false, false))){output.writer.print(text);}
	}
	private static void fails(final Runnable run){
		try{run.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Malformed or leaking gate input must fail loudly");
	}
	private static void check(final boolean valid, final String message){if(!valid){throw new AssertionError(message);}}
}
