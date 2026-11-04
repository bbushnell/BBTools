package test.prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import consensus.BaseGraph;
import consensus.BaseGraphHelper;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaSequentialRegenerationExperiment;
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

/** Literal saved-call controls for the development paired-model restoration step.
 * Two lost references share one winner; restoration must occur once. A different
 * model remains regenerated, including its HBM counts. Zero-loss input is passed
 * through, while an already-original lost winner stops instead of repeating work.
 * These test saved-call transformation, not genomic calling or accuracy.
 * Preserved with the caller-free saved-GFF restoration mode. The normal public
 * reducer performs new caller passes and does not expose this standalone mode;
 * these controls remain useful when saved-call restoration is needed again.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCallerReversionTest{

	/** Owns diagnostic output while creating one finite native fixture. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaCallerReversionTest.class, false);
		try{run(pp.args);}
		finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final String[] args){
		final Parser parser=new Parser();
		for(String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown fixture argument: "+arg);}
		}
		check(parser.out1!=null, "Fixture requires a fresh out= directory");
		final File root=new File(parser.out1).getAbsoluteFile();
		check(root.mkdir(), "Fixture directory must be fresh");
		final File original=new File(root, "original"), current=new File(root, "current"), expected=new File(root, "expected");
		writeLibrary(original, new byte[]{'A', 'G'}, new boolean[]{false, false});
		writeLibrary(current, new byte[]{'C', 'T'}, new boolean[]{true, true});
		writeLibrary(expected, new byte[]{'A', 'T'}, new boolean[]{false, true});
		final ArrayList<Read> genes=new ArrayList<Read>();
		for(byte base:new byte[]{'A', 'G', 'T'}){
			final byte[] sequence=new byte[60];
			Arrays.fill(sequence, base);
			genes.add(new Read(sequence, null, "gene"+(char)base, genes.size()));
		}
		final String corpus=new File(root, "genes.fa").getPath();
		writeFasta(corpus, genes);
		final NcrnaSupport.Graph graph=NcrnaSupport.scan(FileFormat.testInput(corpus, FileFormat.FASTA,
			null, true, true), NcrnaSupport.loadLibrary(new File(original, "consensus.fa").getPath()), .70f, 2, 1)
			.subset(new boolean[]{true, true});
		final String graphPath=new File(root, "original.support.tsv").getPath();
		NcrnaSupportIO.write(graph, graphPath);
		write(new File(root, "truth.tsv"), "seqid\tfamily\tstart\tstop\tstrand\n"
			+"chr\tRF00162\t100\t159\t+\nchr\tRF00162\t300\t359\t+\nchr\tRF00162\t500\t559\t+\n");
		write(new File(root, "baseline.gff"), call(100, "model0")+call(300, "model1")+call(500, "model0"));
		write(new File(root, "missing.gff"), call(300, "model1"));
		write(new File(root, "unknown.gff"), call(100, "absent_model")+call(300, "model1")+call(500, "model0"));
		final File restored=new File(root, "restored");
		invoke(root, original, current, "baseline.gff", "missing.gff", restored);
		compare(new File(expected, "consensus.fa"), new File(restored, "consensus.fa"));
		compare(new File(expected, "models.hbmc"), new File(restored, "models.hbmc"));
		summary(restored, 1, 2, 0, 1, 1);
		final File passed=new File(root, "passed");
		invoke(root, original, restored, "baseline.gff", "baseline.gff", passed,
			"genes="+corpus, "graph="+graphPath, "expectedgenes=3", "cutoff=0.70", "t=2");
		compare(new File(restored, "consensus.fa"), new File(passed, "consensus.fa"));
		compare(new File(restored, "models.hbmc"), new File(passed, "models.hbmc"));
		summary(passed, 3, 0, 1, 0, 1);
		try(Input input=new Input(new File(passed, "support.summary.tsv").getPath())){
			final LineParser1 supportRow=new LineParser1('\t');
			supportRow.set(input.reader.nextLine());
			check(supportRow.termEquals("retention_pass", 7), "Support measurement has its own acceptance column");
			supportRow.set(input.reader.nextLine());
			check(supportRow.parseInt(0)==3 && supportRow.parseInt(1)==2 && supportRow.parseInt(2)==2
				&& supportRow.parseInt(3)==1 && supportRow.parseInt(4)==1 && supportRow.parseInt(5)==1 && supportRow.parseInt(7)==0,
				"Equal corpus support totals must not hide one lost original gene and one newly supported gene");
		}
		final File stalled=new File(root, "stalled");
		boolean stopped=false;
		try{
			invoke(root, original, restored, "baseline.gff", "missing.gff", stalled);
		}catch(IllegalStateException failure){
			check(failure.getMessage().contains("made no progress"), "The stop must explain the exhausted restoration set");
			stopped=true;
		}
		check(stopped && !new File(stalled, "consensus.fa").exists(), "No-progress must not export another allegedly useful proposal");
		summary(stalled, 1, 2, 1, 0, 1);
		boolean unknown=false;
		try{
			invoke(root, original, current, "unknown.gff", "missing.gff", new File(root, "unknown"));
		}catch(IllegalArgumentException failure){
			check(failure.getMessage().contains("Baseline winner is absent"), "Unknown winner must fail its identity join");
			unknown=true;
		}
		check(unknown, "A missing model identity cannot be guessed or silently skipped");
		System.out.println("CALLER_REVERSION_TEST_PASS shared-winner paired-counts passthrough no-progress unknown-model support-loss");
	}

	/** Exercises the actual shared development entry point, including its mode dispatch. */
	private static void invoke(final File root, final File original, final File current,
			final String baseline, final String calls, final File out, final String... extra){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(new String[]{"mode=protect",
			"in="+new File(current, "consensus.fa"), "hbm="+new File(current, "models.hbmc"),
			"ref="+new File(original, "consensus.fa"), "originalhbm="+new File(original, "models.hbmc"),
			"baseline="+new File(root, baseline), "calls="+new File(root, calls),
			"truth="+new File(root, "truth.tsv"), "family=SAM-I", "truthfamily=RF00162",
			"expectedmodels=2", "out="+out}));
		for(String arg:extra){args.add(arg);}
		try{
			NcrnaSequentialRegenerationExperiment.main(args.toArray(new String[args.size()]));
		}catch(RuntimeException failure){
			throw failure;
		}catch(Exception failure){
			throw new IllegalStateException("Unexpected checked failure in saved-call restoration", failure);
		}
	}

	/** Checks literal reference counts and deduplicated model-restoration counts. */
	private static void summary(final File root, final int tp, final int lost,
			final int before, final int added, final int after){
		try(Input input=new Input(new File(root, "summary.tsv").getPath())){
			final LineParser1 parser=new LineParser1('\t');
			parser.set(input.reader.nextLine());
			check(parser.termEquals("lost", 8) && parser.termEquals("newly_restored", 11), "Fixture must read the declared native summary schema");
			parser.set(input.reader.nextLine());
			check(parser.parseInt(0)==2 && parser.parseInt(1)==3 && parser.parseInt(2)==tp && parser.parseInt(8)==lost,
				"Saved-call grading must preserve the literal model/reference population");
			check(parser.parseInt(10)==before && parser.parseInt(11)==added && parser.parseInt(12)==after,
				"Two missing references sharing a winner require only one paired restoration");
			check(input.reader.nextLine()==null, "Exactly one summary row describes this step");
		}
	}

	/** Native paired inputs deliberately differ in both sequence and positional HBM counts. */
	private static void writeLibrary(final File root, final byte[] bases, final boolean[] regenerated){
		check(root.mkdir(), "Each literal library requires a fresh directory");
		final ArrayList<Read> reads=new ArrayList<Read>();
		final ArrayList<BaseGraph> graphs=new ArrayList<BaseGraph>();
		for(int model=0; model<2; model++){
			final byte[] sequence=new byte[60];
			Arrays.fill(sequence, bases[model]);
			final String name="model"+model+(regenerated[model] ? " n=5" : " n=4");
			reads.add(new Read(sequence, null, name, model));
			final BaseGraph graph=new BaseGraph(name, sequence, null, model, 0);
			for(int pos=0; pos<sequence.length; pos++){
				if(regenerated[model]){BaseGraphHelper.setCounts(graph, pos, 2, 3, 5, 7, 11);}
				else{BaseGraphHelper.setCounts(graph, pos, 11, 7, 5, 3, 2);}
			}
			BaseGraphHelper.initForScoring(graph);
			graphs.add(graph);
		}
		final String fasta=new File(root, "consensus.fa").getPath();
		writeFasta(fasta, reads);
		NucleotideHbmCodec.writeModels(graphs, new File(root, "models.hbmc").getPath(), false);
	}

	/** Writes a complete finite sequence list with the native sequence writer. */
	private static void writeFasta(final String fasta, final ArrayList<Read> reads){
		final Writer writer=WriterFactory.makeWriter(FileFormat.testOutput(fasta, FileFormat.FASTA,
			null, true, false, false, true), null, 1, null, false);
		try{
			writer.start();
			writer.add(reads, 0);
		}catch(Throwable failure){KillSwitch.exceptionKill(failure);}
		check(!writer.poisonAndWait(), "Native fixture FASTA must complete before it is consumed");
	}

	private static String call(final int start, final String model){
		return "chr\tCallGenes\tRNA\t"+start+"\t"+(start+59)+"\t.\t+\t.\tfamily:SAM-I;model:"+model+"\n";
	}

	private static void write(final File file, final String text){
		try(Output output=new Output(file.getPath())){output.writer.print(text);}
	}

	/** Compares all native record bytes, including every positional count and metadata field. */
	private static void compare(final File expected, final File actual){
		try(Input a=new Input(expected.getPath()); Input b=new Input(actual.getPath())){
			for(byte[] line=a.reader.nextLine(); line!=null; line=a.reader.nextLine()){
				check(Arrays.equals(line, b.reader.nextLine()), "Paired restoration differs from the independent literal library: "+actual);
			}
			check(b.reader.nextLine()==null, "Fixed model count must not grow during restoration");
		}
	}

	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
