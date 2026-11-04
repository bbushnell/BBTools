package test.prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;

import consensus.BaseGraph;
import consensus.BaseGraphHelper;
import dna.AminoAcid;
import fileIO.FileFormat;
import idaligner.QuantumAligner;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.ReduceNcrnaModels;
import prok.NcrnaSupport;
import prok.NcrnaSupportIO;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NucleotideHbmCodec;
import shared.KillSwitch;
import shared.Shared;
import stream.Read;
import stream.Writer;
import stream.WriterFactory;

/** Finite controls using the actual public consensus tool and Quantum engine.
 * An ambiguous parent supports two disagreeing alleles; a majority allele must
 * be reverted if it strands the minority. Identical members instead keep the
 * regenerated allele. A short group exercises the public tool's no-model result.
 * A supporting alternative receives a migrated member and is regenerated later.
 * Empty models and pre-existing unsupported genes remain explicit.
 * These synthetic controls are not resource construction or accuracy estimates.
 * Exercises the public fixed-library mode, including saved ORIGINAL/REDUCED support.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSequentialRegenerationTest{

	/** Runs five small controls under a fresh output directory. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaSequentialRegenerationTest.class, false);
		try{run(pp.args);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/** Uses native a/b options and the installed public consensus-builder locator. */
	private static void run(final String[] args) throws Exception{
		final Parser parser=new Parser();
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown fixture parameter: "+arg);}
		}
		check(parser.out1!=null, "Fixture requires out=");
		final File root=new File(parser.out1).getAbsoluteFile();
		check(root.mkdir(), "Fixture output directory must be fresh");
		final File recipe=new File(root, "consensus.config");
		try(Output output=new Output(recipe.getPath())){
			output.writer.print("passes=2\nminid=0.3\nendtrim=0\nstochastic=f\n");
		}
		test(root, recipe, "keep", 200, false);
		test(root, recipe, "revert", 200, true);
		test(root, recipe, "move", 200, true);
		test(root, recipe, "short", 40, false);
		finalMigration(root, recipe);
		graphControls(root);
		System.out.println("SEQUENTIAL_REGENERATION_TEST_PASS keep/revert/move/short/final-move/reduced-input");
	}

	/** Each original ambiguous model has exactly 75% ungapped expected identity to either allele. */
	private static void test(final File root, final File recipe,
			final String label, final int length, final boolean minority) throws Exception{
		final File directory=new File(root, label);
		check(directory.mkdir(), "Each control owns a fresh directory");
		final byte[] allele=new byte[length], parent=new byte[length], alternative=new byte[length], other=new byte[length];
		final byte[] alphabet={'A', 'C', 'G', 'T'};
		final Random random=new Random(71);
		for(int i=0; i<length; i++){
			final int base=random.nextInt(4);
			allele[i]=alphabet[base];
			parent[i]=(i&1)==0 ? (byte)'N' : allele[i];
			alternative[i]=(i&1)==0 ? alphabet[(base+1)&3] : allele[i];
			other[i]='T';
		}
		if(label.equals("move")){
			System.arraycopy(alternative, 0, other, 0, length);
			for(int i=0; i<length; i++){
				if(i%7==0 || i%7==3){other[i]=alphabet[(AminoAcid.baseToNumber[other[i]]+1)&3];}
			}
			final float parentIdentity=QuantumAligner.alignStatic(alternative, parent, null);
			final float alternateIdentity=QuantumAligner.alignStatic(alternative, other, null);
			check(alternateIdentity>=.70f && alternateIdentity<parentIdentity,
				"Migration control requires a supported alternative weaker than the original parent: "+alternateIdentity+" / "+parentIdentity);
		}
		final String fasta=new File(directory, "base.fa").getPath();
		final String hbm=new File(directory, "base.hbmc").getPath();
		final ArrayList<Read> references=new ArrayList<Read>();
		references.add(new Read(parent, null, "model0", 0));
		references.add(new Read(other, null, "model1", 1));
		writeLibrary(fasta, hbm, references);
		final ArrayList<Read> genes=new ArrayList<Read>();
		for(int i=0; i<3; i++){genes.add(new Read(allele.clone(), null, "member"+i, i));}
		if(minority){genes.add(new Read(alternative, null, "minority", genes.size()));}
		final byte[] unsupported=new byte[length];
		Arrays.fill(unsupported, (byte)'G');
		genes.add(new Read(unsupported, null, "unsupported", genes.size()));
		final String corpus=new File(directory, "genes.fa").getPath();
		writeFasta(corpus, genes);
		final File result=new File(directory, "result");
		ReduceNcrnaModels.main(new String[]{"in="+corpus, "ref="+fasta, "hbm="+hbm,
			"out="+result, "consensusconfig="+recipe, "reduce=f", "regenerate=t", "callerverify=f",
			"cutoff=0.70", "consensusheap=512m", "t=2"});
		checkSummary(new File(result, "regeneration/summary.tsv"), genes.size(), label);
		if(label.equals("short")){
			try(Input input=new Input(new File(result, "regeneration/models.tsv").getPath())){
				input.reader.nextLine();
				final LineParser1 parser=new LineParser1('\t');
				parser.set(input.reader.nextLine());
				check(parser.termEquals("model0", 0) && parser.termEquals("NO_CANDIDATE_REVERTED", 1),
					"The short member group must exercise the consensus builder's no-model branch");
			}
		}
		final NcrnaSupport.Library written=NcrnaSupport.loadLibrary(new File(result, "consensus.fa").getPath());
		check(Arrays.equals(written.names, new String[]{"model0", "model1"}), "Model names/order/count must remain fixed");
		if(!label.equals("move")){check(Arrays.equals(written.bases[1], other), "Empty model remains unchanged");}
		if(label.equals("keep") || label.equals("move")){
			check(Arrays.equals(written.bases[0], allele), "Identical supported members regenerate their exact allele");
			check(!Arrays.equals(parent, allele), "The accepted control must actually change the original consensus");
			if(label.equals("move")){check(Arrays.equals(written.bases[1], alternative), "A later model must regenerate from its newly arrived member");}
		}else{
			final NcrnaSupport.Library original=NcrnaSupport.loadLibrary(fasta);
			check(Arrays.deepEquals(original.bases, written.bases) && Arrays.equals(original.headers, written.headers),
				"Rollback must retain both complete FASTA records regardless of line wrapping");
			compare(hbm, new File(result, "models.hbmc").getPath());
		}
		System.out.println("SEQUENTIAL_CONTROL_PASS "+label);
	}

	/** A short second model later becomes a full consensus that attracts a non-member.
	 * The probe stays with model0 during both steps and moves only in the final scan.
	 * Supplying a native REDUCED graph also exercises the real runs' NONE/NA input path.
	 */
	private static void finalMigration(final File root, final File recipe) throws Exception{
		final File directory=new File(root, "final_move");
		check(directory.mkdir(), "Final-migration control needs its own fresh directory");
		final byte[] x=new byte[400];
		final byte[] alphabet={'A', 'C', 'G', 'T'};
		final Random random=new Random(71);
		for(int i=0; i<x.length; i++){x[i]=alphabet[random.nextInt(4)];}
		final byte[] z=x.clone(), y=x.clone();
		for(int i=100; i<=112; i+=3){z[i]=alphabet[(AminoAcid.baseToNumber[z[i]]+1)&3];}
		y[125]=alphabet[(AminoAcid.baseToNumber[y[125]]+1)&3];
		final ArrayList<Read> references=new ArrayList<Read>();
		references.add(new Read(x, null, "model0", 0));
		references.add(new Read(Arrays.copyOfRange(y, 100, 150), null, "model1", 1));
		final String fasta=new File(directory, "base.fa").getPath();
		final String hbm=new File(directory, "base.hbmc").getPath();
		writeLibrary(fasta, hbm, references);
		final ArrayList<Read> genes=new ArrayList<Read>();
		for(int i=0; i<3; i++){genes.add(new Read(z.clone(), null, "z"+i, genes.size()));}
		genes.add(new Read(x.clone(), null, "probe", genes.size()));
		for(int i=0; i<3; i++){genes.add(new Read(y.clone(), null, "y"+i, genes.size()));}
		final byte[] unsupported=new byte[x.length];
		Arrays.fill(unsupported, (byte)'G');
		genes.add(new Read(unsupported, null, "unsupported", genes.size()));
		final String corpus=new File(directory, "genes.fa").getPath();
		writeFasta(corpus, genes);
		final NcrnaSupport.Graph graph=NcrnaSupport.scan(FileFormat.testInput(corpus, FileFormat.FASTA,
			null, true, true), NcrnaSupport.loadLibrary(fasta), .70f, 2, 1).subset(new boolean[]{true, true});
		check(graph.reduced, "The supplied-input control must exercise REDUCED graph semantics");
		for(int row=0; row<graph.queryNames.length; row++){
			final String id=graph.queryNames[row];
			final int expected=id.equals("unsupported") ? -1 : id.startsWith("y") ? 1 : 0;
			check(graph.bestModels[row]==expected, "Controlled initial winner differs for "+id+": "+graph.bestModels[row]);
		}
		final String graphPath=new File(directory, "members.tsv").getPath();
		NcrnaSupportIO.write(graph, graphPath);
		final File result=new File(directory, "result");
		ReduceNcrnaModels.main(new String[]{"in="+corpus, "ref="+fasta, "hbm="+hbm,
			"graph="+graphPath, "out="+result, "consensusconfig="+recipe,
			"reduce=f", "regenerate=t", "callerverify=f", "cutoff=0.70", "consensusheap=512m", "t=2"});
		checkSummary(new File(result, "regeneration/summary.tsv"), genes.size(), "final_move");
		checkZeroInitialWork(result);
		final NcrnaSupport.Library written=NcrnaSupport.loadLibrary(new File(result, "consensus.fa").getPath());
		check(Arrays.equals(written.bases[0], z) && Arrays.equals(written.bases[1], y),
			"Each full regenerated consensus must equal its identical majority");
		boolean found=false;
		try(Input input=new Input(new File(result, "regeneration/series.assignments.tsv").getPath())){
			input.reader.nextLine();
			final LineParser1 parser=new LineParser1('\t');
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				if(parser.termEquals("probe", 0)){
					check(parser.termEquals("model0", 1) && parser.termEquals("model0", 2) && parser.termEquals("model1", 3),
						"The probe must migrate only during the final whole-corpus pass");
					found=true;
				}
			}
		}
		check(found, "Final assignments must include the probe");
		System.out.println("SEQUENTIAL_CONTROL_PASS final_move");
	}

	/** Header checks make the fixed native report columns an explicit test contract. */
	private static void checkSummary(final File file, final int genes, final String label){
		try(Input input=new Input(file.getPath())){
			final LineParser1 parser=new LineParser1('\t');
			parser.set(input.reader.nextLine());
			check(parser.termEquals("models", 1) && parser.termEquals("initial_supported", 12)
				&& parser.termEquals("initial_unsupported", 13) && parser.termEquals("candidate_supported", 14)
				&& parser.termEquals("kept", 2) && parser.termEquals("reverted", 3) && parser.termEquals("empty", 4)
				&& parser.termEquals("final_migrations", 6) && parser.termEquals("member_migrations", 5),
				"Read the declared native summary schema");
			parser.set(input.reader.nextLine());
			check(parser.parseInt(1)==2 && parser.parseInt(12)==genes-1
				&& parser.parseInt(13)==1 && parser.parseInt(14)==genes-1, "Retain every initially supported gene and separate the unsupported gene");
			final boolean move=label.equals("move"), finalMove=label.equals("final_move"), keep=label.equals("keep");
			final boolean bothKept=move || finalMove;
			check(parser.parseInt(2)==(bothKept ? 2 : keep ? 1 : 0) && parser.parseInt(3)==(bothKept || keep ? 0 : 1)
				&& parser.parseInt(4)==(bothKept ? 0 : 1), "Expected the declared keep/revert/empty decisions");
			check(parser.parseInt(5)==(move ? 1 : 0) && parser.parseInt(6)==(finalMove ? 1 : 0)
				&& parser.parseInt(7)==0 && parser.parseString(0).contains("not caller-verified"),
				"Preserve migration expectations and explicitly report zero fresh caller passes");
			check(input.reader.nextLine()==null, "One summary row per completed experiment");
		}
	}

	/** A loaded graph must produce the same subset while spending zero new initial alignments. */
	private static void graphControls(final File root){
		final File fixture=new File(root, "final_move");
		final String corpus=new File(fixture, "genes.fa").getPath();
		final NcrnaSupport.Library library=NcrnaSupport.loadLibrary(new File(fixture, "base.fa").getPath());
		final NcrnaSupport.Graph original=NcrnaSupport.scan(FileFormat.testInput(corpus, FileFormat.FASTA,
			null, true, true), library, .70f, 1, 1);
		final String saved=new File(fixture, "original.support.tsv").getPath();
		NcrnaSupportIO.write(original, saved);
		final File fresh=new File(root, "fresh_subset");
		reduce(fixture, null, .70f, fresh);
		for(String graph:new String[]{saved, new File(fixture, "members.tsv").getPath()}){
			final File loaded=new File(root, graph.equals(saved) ? "original_subset" : "reduced_subset");
			reduce(fixture, graph, .70f, loaded);
			compare(new File(fresh, "consensus.fa").getPath(), new File(loaded, "consensus.fa").getPath());
			compare(new File(fresh, "members.tsv").getPath(), new File(loaded, "members.tsv").getPath());
			checkZeroInitialWork(loaded);
		}
		fails(()->reduce(fixture, saved, .75f, new File(root, "bad_cutoff")), "cutoff");
		final String oldModel=original.modelNames[0];
		original.modelNames[0]="wrong_model";
		final String wrongModel=new File(fixture, "wrong_model.tsv").getPath();
		NcrnaSupportIO.write(original, wrongModel);
		original.modelNames[0]=oldModel;
		fails(()->reduce(fixture, wrongModel, .70f, new File(root, "bad_model")), "model identities");
		original.queryNames[0]="wrong_gene";
		final String wrongGene=new File(fixture, "wrong_gene.tsv").getPath();
		NcrnaSupportIO.write(original, wrongGene);
		fails(()->reduce(fixture, wrongGene, .70f, new File(root, "bad_gene")), "corpus gene");
		final String malformed=new File(fixture, "malformed.tsv").getPath();
		try(Output out=new Output(malformed)){out.writer.print("not a support graph\n");}
		fails(()->reduce(fixture, malformed, .70f, new File(root, "bad_format")), "graph");
		System.out.println("SAVED_SUPPORT_GRAPH_TEST_PASS original reduced parity zero-alignments wrong-model wrong-gene wrong-cutoff malformed");
	}

	private static void reduce(final File fixture, final String graph, final float cutoff, final File output){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+new File(fixture, "genes.fa"),
			"ref="+new File(fixture, "base.fa"), "cutoff="+cutoff, "target=0.8", "t=1", "out="+output));
		if(graph!=null){args.add("graph="+graph);}
		ReduceNcrnaModels.main(args.toArray(new String[args.size()]));
	}

	private static void checkZeroInitialWork(final File output){
		boolean zero=false, labeled=false;
		try(Input input=new Input(new File(output, "summary.txt").getPath())){
			final LineParser1 parser=new LineParser1('\t');
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				if(parser.termEquals("initial_pair_alignments", 0)){zero=parser.parseLong(1)==0;}
				if(parser.termEquals("initial_support", 0)){labeled=parser.parseString(1).contains("no fresh alignments performed");}
			}
		}
		check(zero && labeled, "Loaded graph work must be zero and distinguished from fresh alignments");
	}

	private static void fails(final Runnable action, final String expected){
		try{action.run();}catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(expected), "Graph rejection must explain the violated contract: "+failure);
			return;
		}
		throw new AssertionError("Invalid saved graph was accepted: "+expected);
	}

	/** Compares independent original and rolled-back native records, including every HBM count. */
	private static void compare(final String original, final String result){
		try(Input first=new Input(original); Input second=new Input(result)){
			for(byte[] a=first.reader.nextLine(); a!=null; a=first.reader.nextLine()){
				check(Arrays.equals(a, second.reader.nextLine()), "Rollback must retain the entire original paired resource: "+result);
			}
			check(second.reader.nextLine()==null, "Rollback must not append another model");
		}
	}

	/** Generates the literal synthetic input through the same public sequence I/O interface. */
	private static void writeFasta(final String path, final ArrayList<Read> reads){
		final Writer writer=WriterFactory.makeWriter(FileFormat.testOutput(path, FileFormat.FASTA, null,
			true, false, false, true), null, 1, null, false);
		try{
			writer.start();
			writer.add(reads, 0);
		}catch(Throwable failure){KillSwitch.exceptionKill(failure);}
		check(!writer.poisonAndWait(), "Synthetic FASTA output must complete");
	}

	/** Writes independent baseline sequences and their paired native HBM records. */
	private static void writeLibrary(final String fasta, final String hbm, final ArrayList<Read> references){
		writeFasta(fasta, references);
		final ArrayList<BaseGraph> models=new ArrayList<BaseGraph>();
		for(Read reference:references){
			final BaseGraph model=new BaseGraph(reference.id, reference.bases, null, reference.numericID, 0);
			BaseGraphHelper.initForScoring(model);
			models.add(model);
		}
		NucleotideHbmCodec.writeModels(models, hbm, false);
	}

	private static void check(final boolean condition, final String reason){
		if(!condition){throw new AssertionError(reason);}
	}
}
