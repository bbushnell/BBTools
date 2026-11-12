package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import prok.BuildNcrnaReference;
import prok.CompareNcrnaScenarios;
import prok.GradeNcrnaReference;
import prok.SelectNcrnaScenarios;
import prok.WriteNcrnaReference;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Literal selection and downstream matching oracle, including an empty query
 * input, aliases, disabled families and duplicate calls. No genome is opened.
 * @author Brian Bushnell, Raiden
 */
final class SelectNcrnaScenariosTest{
	SelectNcrnaScenariosTest(final File root_){root=root_;}

	void run(){
		check(root.mkdir(), "Selector fixture must be new");
		prepare();
		invoke("result");
		final String q='\t'+SelectNcrnaScenarios.QUALIFICATION;
		check(lines("result/summary.tsv")[1].equals("2\t1\t3\t6\t2\t4"+q),
			"Six global decisions yield two included decisions over both original inputs, including one empty query");
		check(Arrays.equals(lines("result/ga.tsv.gz"), new String[]{CompareNcrnaScenarios.HEADER, lowOne.trim(), empty.trim()}),
			"Selected GA rows remain byte-identical, including the empty input's missed reference");
		check(Arrays.equals(lines("result/margin.tsv.gz"), new String[]{CompareNcrnaScenarios.HEADER, highOne.trim(), empty.trim()}),
			"Both sides select the same keys; dispositions are not rewritten");
		final String[] counts=lines("result/selection_counts.tsv");
		check(find(counts, "FMN\tOUTSIDE_BASELINE_INPUT\t").equals("FMN\tOUTSIDE_BASELINE_INPUT\tscored\tscored\t1"+q), "Outside inputs stay counted");
		check(find(counts, "SAM_I\tDISABLED_FAMILY\t").equals("SAM_I\tDISABLED_FAMILY\tscored\tscored\t1"+q), "Disabled family's apparent miss is outside this baseline");
		check(find(counts, "RNaseP\tUNSUPPORTED_OR_UNDECLARED_FAMILY\t").contains("\tdont_care\tdont_care\t1\t"), "Explicit unsupported family is counted");
		check(find(counts, "tRNA\tUNSUPPORTED_OR_UNDECLARED_FAMILY\t").contains("\tcontext\tcontext\t1\t"), "Global tRNA context needs no invented caller mapping");
		check(find(lines("result/inputs.tsv"), "two\t").endsWith("\t0\t0\t0\t1"+q), "Zero calls do not erase a complete input or its selected reference");
		CompareNcrnaScenarios.main(new String[]{"in="+path("result/ga.tsv.gz"), "in2="+path("result/margin.tsv.gz"),
			"calls="+path("calls.tsv"), "familymap="+path("result/familymap.tsv"), "contigs="+path("contigs.tsv"),
			"out="+path("comparison"), "t=1", "pigz=f", "bgzip=f", "zipthreads=1"});
		final String actual=find(lines("comparison/summary.tsv"), "FMN\t");
		check(actual.equals("FMN\t2\t1\t1\t1\t2\t0\t0\t1\t0\t0\t0\t1\t1\t1\t2\t1\t0\t-1\t-1\t0\t2\t0.500000\t"+CompareNcrnaScenarios.QUALIFICATION),
			"One TP and one duplicate FP become two masked calls; the empty input retains its FN on both sides");
		bad("missing_gate", "Every canonical gate", "callerconfig="+write("missing.config", config.replace("fmn=t\n", "")));
		bad("duplicate_gate", "Missing/duplicate/conflicting", "callerconfig="+write("duplicate.config", config+"fmn=f\n"));
		bad("alias_gate", "Unsupported caller alias", "callerconfig="+write("alias.config", config+"generalncrna=f\n"));
		bad("only_gate", "Unsupported caller alias", "callerconfig="+write("only.config", config+"trnaonly=t\n"));
		bad("nested_config", "Nested or paired", "callerconfig="+write("nested.config", config+"config=other.config\n"));
		bad("euk_profile", "Every canonical gate", "callerconfig="+write("euk.config", config.replace("euk=f", "euk=t")));
		bad("empty_removed", "complete adapter/config population", "calls="+write("short.calls.tsv", calls.substring(0, calls.lastIndexOf("two\t"))));
		bad("both_short", "Shortened manifests", "calls="+path("short.calls.tsv"), "bindings="+write("short.bindings.tsv", bindings.substring(0, bindings.indexOf("\ntwo\t")+1)));
		bad("alias_mismatch", "aliases and ordered", "bindings="+write("bad.bindings.tsv", bindings.replace("original%25one", "different")));
		bad("asymmetric", "Asymmetric scenario", "in2="+write("short.margin.tsv", high.substring(0, high.lastIndexOf("outside\t"))));
		bad("different_key", "identical ordered", "in2="+write("key.margin.tsv", high.replace("same\tone", "changed\tone")));
		bad("different_family", "identical ordered", "in2="+write("family.margin.tsv", high.replace("\tSAM_I\t", "\tFMN\t")));
		bad("noncanonical_key", "canonical evidence TSV", "in="+write("escaped.ga.tsv", low.replace("same\t", "%73ame\t")),
			"in2="+write("escaped.margin.tsv", high.replace("same\t", "%73ame\t")));
		bad("duplicate_included", "Duplicate input-local decision", "in="+write("duplicate.ga.tsv", low+lowOne), "in2="+write("duplicate.margin.tsv", high+highOne));
		final String outside=row("outside", "three", "FMN", "scored");
		bad("duplicate_excluded", "Duplicate input-local decision", "in="+write("excluded.ga.tsv", low+outside), "in2="+write("excluded.margin.tsv", high+outside));
		badQueries("disabled_query", "model:SAM1_1,family:SAM-I", "Saved query contradicts a disabled family");
		badQueries("unknown_query", "model:UNKNOWN_1,family:FMN", "Unknown generic RNA model");
		badQueries("wrong_label", "model:FMN_1,family:SAM-I", "label conflicts with model");
		System.err.println("SELECT_NCRNA_SCENARIOS_TEST_PASS");
	}

	private void prepare(){
		write("one.gff", call("model:FMN_1,family:FMN")+call("model:FMN_duplicate,family:FMN")+call("model:rnasep_1,family:RNase P"));
		write("two.gff", "# empty RNA query remains a complete input\n");
		final String original=path("original%one.fa"), canonical=path("canonical%one.fa"), two=path("two.fa");
		calls=GradeNcrnaReference.MANIFEST_HEADER+"\none\t"+canonical+'\t'+path("one.gff")+"\ntwo\t"+two+'\t'+path("two.gff")+'\n';
		write("calls.tsv", calls);
		bindings=SelectNcrnaScenarios.BINDING_HEADER+'\n'+fields("one", "original", "assembly", original, canonical, "2", path("one.gff"), "GENERIC_RNA_ONLY_ORIGINAL_OTHER_FEATURES_RETAINED")
			+fields("two", "two", "assembly", two, two, "3", path("two.gff"), "GENERIC_RNA_ONLY_ORIGINAL_OTHER_FEATURES_RETAINED");
		write("bindings.tsv", bindings);
		config="in="+original+"\nin="+two+"\nncrna=t\nfmn=t\nsam1=f\neuk=f\nncrnafamily=fmn\n";
		write("caller.config", config);
		write("scope.tsv", SelectNcrnaScenarios.SCOPE_HEADER+"\nFMN\tncrna,fmn\ttrue\nSAM_I\tncrna,sam1\ttrue\nRNaseP\tncrna\tfalse\n");
		write("models.tsv", SelectNcrnaScenarios.MODEL_HEADER+"\nFMN\tFMN\tFMN\nSAM_I\tSAM1\tSAM-I\nRNaseP\trnasep\tRNase P\n");
		write("contigs.tsv", WriteNcrnaReference.CONTIG_HEADER+'\n'+fields("one", canonical, "chr", "100")+fields("two", two, "chr", "100"));
		lowOne=row("same", "one", "FMN", "scored");
		highOne=row("same", "one", "FMN", "dont_care");
		empty=row("same", "two", "FMN", "scored");
		final String remaining=empty+row("disabled", "one", "SAM_I", "scored")+row("unsupported", "one", "RNaseP", "dont_care")
			+row("context", "one", "tRNA", "context")+row("outside", "three", "FMN", "scored");
		low=CompareNcrnaScenarios.HEADER+'\n'+lowOne+remaining;
		high=CompareNcrnaScenarios.HEADER+'\n'+highOne+remaining;
		write("ga.tsv", low);
		write("margin.tsv", high);
	}

	/** Changes one query file and both matching bindings, preserving original inputs. */
	private void badQueries(final String name, final String attributes, final String expected){
		final String gff=write(name+".gff", call(attributes));
		bad(name, expected, "calls="+write(name+".calls.tsv", calls.replace(path("one.gff"), gff)),
			"bindings="+write(name+".bindings.tsv", bindings.replace(path("one.gff"), gff)));
	}
	private void invoke(final String name, final String... extra){
		final String[] base={"in="+path("ga.tsv"), "in2="+path("margin.tsv"), "calls="+path("calls.tsv"), "bindings="+path("bindings.tsv"),
			"scope="+path("scope.tsv"), "models="+path("models.tsv"), "callerconfig="+path("caller.config"), "out="+path(name),
			"t=1", "pigz=f", "bgzip=f", "zipthreads=1"};
		final String[] args=Arrays.copyOf(base, base.length+extra.length);
		System.arraycopy(extra, 0, args, base.length, extra.length);
		SelectNcrnaScenarios.main(args);
	}
	private void bad(final String name, final String expected, final String... extra){
		try{invoke(name, extra);}
		catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(expected), "Unexpected selector failure: "+failure);
			check(!new File(path(name+"/summary.tsv")).exists(), "Invalid subsets must not publish completion");
			return;
		}
		throw new AssertionError("Expected selector failure: "+expected);
	}
	private static String call(final String attributes){return "chr\tBBTools\tRNA\t10\t29\t.\t+\t.\t"+attributes+'\n';}
	private static String row(final String id, final String input, final String family, final String state){
		return fields(id, input, "chr", "10", "29", "+", family, "RESOLVED", state, "diagnostic_rule", "M", state.equals("scored") ? "." : "uncertain", "original:decision", "original:witness");
	}
	private static String fields(final String... fields){
		final ByteBuilder b=new ByteBuilder();
		for(int i=0; i<fields.length; i++){
			if(i>0){b.tab();}
			BuildNcrnaReference.field(b, fields[i]);
		}
		b.nl();
		return new String(b.array, 0, b.length(), StandardCharsets.UTF_8);
	}
	private String path(final String name){return new File(root, name).getAbsolutePath();}
	private String write(final String name, final String text){
		try(Output out=new Output(path(name))){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
		return path(name);
	}
	private String[] lines(final String name){
		final ArrayList<String> rows=new ArrayList<String>();
		try(Input in=new Input(path(name))){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){rows.add(new String(line, StandardCharsets.UTF_8));}
		}
		return rows.toArray(new String[rows.size()]);
	}
	private static String find(final String[] rows, final String prefix){
		for(final String row:rows){if(row.startsWith(prefix)){return row;}}
		throw new AssertionError("Missing expected selection row: "+prefix);
	}
	private static void check(final boolean valid, final String message){if(!valid){throw new AssertionError(message);}}
	private final File root;
	private String config, calls, bindings, low, high, lowOne, highOne, empty;
}
