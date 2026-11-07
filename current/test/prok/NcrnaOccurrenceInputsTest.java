package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.FileFormat;
import gff.CutGff;
import parse.Parse;
import prok.NcrnaLearningContext;
import prok.NcrnaLearningPopulation;
import prok.NcrnaOccurrenceInputs;
import prok.NcrnaSupport.Graph;
import prok.NcrnaSupportIO;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;

/** Literal extraction and learning-consumer fixture. No biological resource build.
 * The 80-base genome is A20 C39 G21; the two endpoint observations are independent
 * coordinate oracles. Inputs and outputs remain under the requested fresh root.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaOccurrenceInputsTest{

	public static void main(final String[] args){
		String out=null;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0], b=split.length>1 ? split[1] : null;
			if(a.equalsIgnoreCase("out")){out=b;}
			else{throw new IllegalArgumentException("Unknown test parameter: "+arg);}
		}
		check(out!=null && new File(out).mkdir(), "Fresh fixture directory required");
		new NcrnaOccurrenceInputsTest(out).run();
		System.err.println("NCRNA_OCCURRENCE_INPUTS_FINITE_PASS");
	}

	private NcrnaOccurrenceInputsTest(final String root_){root=root_;}

	private void run(){
		final String fasta=path("genome.fa"), roster=path("sources.tsv"), occurrences=path("occurrences.tsv");
		write(fasta, ">same\n"+repeat('A', 20)+repeat('C', 39)+repeat('G', 21)+"\n");
		write(roster, NcrnaLearningPopulation.ROSTER_HEADER+"\n0\tgA\tasmA\t1\tbacteria\tP\toriginal.gff\t"
			+fasta+"\t7\t80\t1\n1\tunused\tasmB\t2\tbacteria\tP\toriginal.gff\t"+fasta+"\t0\t80\t1\n");
		final String rows="occ1\tgA\tasmA\t1\tsame\t+\t1\t20\tRF00504\t11111111111111111111\t1\t80\n"
			+"occ2\tgA\tasmA\t1\tsame\t-\t60\t80\tRF00504\t11111111111111111111\t2\t80\n";
		write(occurrences, OCCURRENCE_HEADER+"\n"+rows);
		final String prepared=path("prepared");
		NcrnaOccurrenceInputs.main(new String[]{"in="+occurrences, "roster="+roster, "out="+prepared, "flank=5", "expected=2"});
		expectLine(prepared+"/summary.tsv", "2\t1\t2\t5\t51");
		expectLine(prepared+"/occurrence_loci.tsv", "ID=occ1\tocc1\tgA\tasmA\t1\tbacteria\tP\tsame\t+\t1\t20\t20\t80\t0\t5");
		expectLine(prepared+"/occurrence_loci.tsv", "ID=occ2\tocc2\tgA\tasmA\t1\tbacteria\tP\tsame\t-\t60\t80\t21\t80\t0\t5");
		final String contexts=path("contexts.fa"), cores=path("cores.fa");
		CutGff.main(new String[]{"config="+prepared+"/sources.config", "out="+contexts, "types=ncRNA", "flank=5",
			"gccontig=caller", "filename=f", "banpartial=f", "align=f", "minlen=1", "maxns=-1", "maxnfraction=-1", "t=1", "ow=f"});
		CutGff.main(new String[]{"in="+contexts, "gff="+prepared+"/core.gff", "out="+cores, "types=ncRNA", "flank=0",
			"gccontig=f", "filename=f", "banpartial=f", "align=f", "minlen=1", "maxns=-1", "maxnfraction=-1", "t=1", "ow=f"});
		final String graphPath=path("graph.tsv");
		write(graphPath, graph());
		final Graph graph=NcrnaSupportIO.read(graphPath);
		final String[] models={"model"};
		final NcrnaLearningPopulation population=NcrnaLearningPopulation.loadForRefit(graph, models,
			prepared+"/genomes.tsv", prepared+"/occurrence_loci.tsv", null);
		check(population.ids.length==2 && population.eligible(0) && population.eligible(1), "Every admitted occurrence remains present and eligible");
		check(Arrays.equals(population.length, new int[]{20, 21}) && Arrays.equals(population.leftFlank, new int[]{0, 0})
			&& Arrays.equals(population.rightFlank, new int[]{5, 5}), "Both-strand edge geometry must match literal lengths");
		final boolean trim=Shared.TRIM_READ_DESCRIPTION;
		Shared.TRIM_READ_DESCRIPTION=false;
		try{
			final ArrayList<Read> reads=ReadInputStream.toReads(contexts, FileFormat.FASTA, -1);
			check(reads.size()==2, "Both contexts must be extracted");
			final NcrnaLearningContext context=new NcrnaLearningContext();
			for(final Read read:reads){
				context.set(read, population);
				check(context.hasFeatures() && context.gc==0.75f && context.start==0, "Actual extracted headers must join and carry whole-contig caller GC");
				final String expected=context.gene==0 ? repeat('A', 20)+repeat('C', 5) : repeat('C', 21)+repeat('G', 5);
				check(Arrays.equals(read.bases, expected.getBytes(StandardCharsets.US_ASCII)), "Exact sense-oriented flanked bases");
			}
			final ArrayList<Read> coreReads=ReadInputStream.toReads(cores, FileFormat.FASTA, -1);
			check(coreReads.size()==2, "Both exact cores must survive the second native extraction");
			for(final Read read:coreReads){
				final boolean first=read.id.startsWith("ID=occ1 ");
				check(first || read.id.startsWith("ID=occ2 "), "Exact core ID must retain its admitted occurrence");
				check(Arrays.equals(read.bases, repeat(first ? 'A' : 'C', first ? 20 : 21).getBytes(StandardCharsets.US_ASCII)), "Literal exact core sequence");
			}
		}finally{Shared.TRIM_READ_DESCRIPTION=trim;}
		badMetadata(graph, models, prepared);
		badPreparation(roster, rows);
	}

	/** Mismatched names, lengths and reverse-strand flanks must fail, never be dropped. */
	private void badMetadata(final Graph graph, final String[] models, final String prepared){
		final String ledger=readText(prepared+"/occurrence_loci.tsv"), bad=path("bad-loci.tsv");
		write(bad, ledger.replace("ID=occ1\tocc1", "ID=occ1\tother"));
		fails(()->NcrnaLearningPopulation.loadForRefit(graph, models, prepared+"/genomes.tsv", bad, null));
		write(bad, ledger.replace("\t60\t80\t21\t80\t0\t5", "\t60\t80\t22\t80\t0\t5"));
		fails(()->NcrnaLearningPopulation.loadForRefit(graph, models, prepared+"/genomes.tsv", bad, null));
		write(bad, ledger.replace("\t60\t80\t21\t80\t0\t5", "\t60\t80\t21\t80\t5\t0"));
		fails(()->NcrnaLearningPopulation.loadForRefit(graph, models, prepared+"/genomes.tsv", bad, null));
		write(bad, ledger.replace("\tasmA\t", "\twrong\t"));
		fails(()->NcrnaLearningPopulation.loadForRefit(graph, models, prepared+"/genomes.tsv", bad, null));
	}

	private void badPreparation(final String roster, final String rows){
		final String bad=path("bad-occurrences.tsv"), output=path("bad-output");
		write(bad, OCCURRENCE_HEADER+"\n"+rows.replace("occ2\t", "occ1\t"));
		fails(()->NcrnaOccurrenceInputs.main(new String[]{"in="+bad, "roster="+roster, "out="+output, "expected=2", "flank=5"}));
		check(!new File(output).exists(), "Duplicate occurrence must fail before output");
		write(bad, OCCURRENCE_HEADER+"\n"+rows.replace("\tgA\t", "\tunknown\t"));
		fails(()->NcrnaOccurrenceInputs.main(new String[]{"in="+bad, "roster="+roster, "out="+output, "expected=2", "flank=5"}));
		check(!new File(output).exists(), "Unknown source must fail before output");
	}

	private static String graph(){
		final ByteBuilder b=new ByteBuilder().append("#ncrna_support_graph\t2\n#scope\tORIGINAL\n#aligner\tquantum\n#cutoff_a48\t")
			.appendFloatA48(0.7f).append("\n#models\t1\n#model\t0\tmodel\n#queries\t2\nquery_id\tbest_model_index\tbest_identity_a48\tsupporters\n");
		for(int i=1; i<=2; i++){b.append("ID=occ").append(i).append("\t0\t").appendFloatA48(0.8f).append("\t0:").appendFloatA48(0.8f).nl();}
		return b.toString();
	}
	private static String repeat(final char c, final int count){
		final char[] chars=new char[count];
		Arrays.fill(chars, c);
		return new String(chars);
	}
	private static void write(final String path, final String text){NcrnaLearningPopulationTest.write(path, text);}
	private static String readText(final String path){
		final ByteBuilder b=new ByteBuilder();
		try(Input input=new Input(path)){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){b.append(line).nl();}
		}
		return b.toString();
	}
	private static void expectLine(final String path, final String expected){
		check(readText(path).contains(expected+"\n"), "Expected literal row: "+expected);
	}
	private static void fails(final Runnable action){
		try{action.run();}
		catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Malformed occurrence metadata must fail");
	}
	private String path(final String name){return new File(root, name).getAbsolutePath();}
	private static void check(final boolean value, final String why){if(!value){throw new AssertionError(why);}}
	private final String root;
	private static final String OCCURRENCE_HEADER="occurrence_id\tgenome\tassembly\ttaxid\tseqid\tstrand\tstart1\tstop1\tmodel_id\tmodel_sha80\trepresentative_physical_id\tcontig_length";
}
