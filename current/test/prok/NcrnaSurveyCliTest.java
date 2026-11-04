package test.prok;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.PreParser;
import prok.NcrnaFamilySurvey;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Real survey CLI/config/log checks on the existing three-genome fixture.
 * One case per JVM reclaims the baseline's deliberately observed selected-log leak
 * at exit. Scientific output parity and launcher JVM settings are checked separately.
 * @author Raiden */
public final class NcrnaSurveyCliTest{

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaSurveyCliTest.class, false);
		try(AutoCloseable log=()->Shared.closeStream(pp.outstream)){
			run(pp.args);
		}
	}
	private static void run(String[] args) throws Exception{
		String root=null, input=null, label=null, expectation=null;
		for(String arg:args){
			if(arg.startsWith("outdir=")){root=arg.substring(7);}
			else if(arg.startsWith("inputroot=")){input=arg.substring(10);}
			else if(arg.startsWith("case=")){label=arg.substring(5);}
			else if(arg.startsWith("expect=")){expectation=arg.substring(7);}
			else{throw new IllegalArgumentException("Unknown survey CLI fixture option: "+arg);}
		}
		check(root!=null && input!=null && Arrays.asList(CASES).contains(label)
			&& ("before".equals(expectation) || "after".equals(expectation)), "Declare fresh outdir, fixed inputroot, one case and before/after expectation");
		check(new File(root).mkdir(), "CLI case needs its own fresh output directory");
		final boolean after=expectation.equals("after"), borrowed=label.equals("borrowed"), fixture=label.equals("fixture_log");
		final boolean always=label.equals("explicit") || label.equals("default") || label.equals("config") || label.equals("config_override");
		final boolean newly=label.equals("kmg") || label.equals("false_outputs");
		final boolean success=always || after && newly;
		final String prefix=root+"/survey", log=new File(root, "selected.log").getCanonicalPath();
		final ArrayList<String> command=command(root, input, prefix, label);
		if(!borrowed){command.add("outstream="+log);}
		Shared.setThreads(3);ReadWrite.setZipThreads(1);
		final PrintStream prior=System.err;
		final ByteArrayOutputStream diagnostics=new ByteArrayOutputStream();
		Throwable problem=null;
		try(PrintStream capture=new PrintStream(diagnostics, true, "UTF-8")){
			System.setErr(capture);
			try{
				if(fixture){NcrnaFamilyExportTest.main(new String[]{"outstream="+log, "unknown=1"});}
				else{NcrnaFamilySurvey.main(command.toArray(new String[0]));}
			}catch(IllegalArgumentException failure){problem=failure;}
			System.err.println("SURVEY_CLI_BORROWED_STDERR_OPEN");
			check(!System.err.checkError(), "Tool cleanup must leave borrowed stderr usable");
		}finally{System.setErr(prior);}
		check((problem==null)==success, "Declared CLI outcome differs: "+label+" "+problem);
		if(after && (label.startsWith("nested") || label.equals("too_deep") || label.equals("inner_log"))){
			check(problem!=null && problem.getMessage().startsWith("Nested config= is unsupported"),
				"Only PreParser expands config files; nested directives must fail before input/output work");
		}
		if(label.equals("parse") || borrowed){check(problem!=null && problem.getMessage().equals("Unknown flag: unknown"), "Selected-log repair must preserve the original parser error");}
		if(label.equals("count")){check(problem!=null && problem.getMessage().equals("Corpus map count changed: 3 expected=4"), "Selected-log cleanup must preserve the original roster error");}
		if(label.equals("worker")){check(problem!=null && problem.getMessage().startsWith("GFF worker failed: "), "The actual worker failure must remain primary");}
		if(fixture){check(problem!=null && problem.getMessage().equals("Unknown fixture argument: unknown=1"), "Outer fixture cleanup preserves its parser error");}
		final int descriptors=descriptors(log), expectedDescriptors=!after && !success && !borrowed ? 1 : 0;
		check(descriptors==expectedDescriptors, "Selected log descriptors: "+label+" actual="+descriptors+", expected="+expectedDescriptors);
		final File selectedLog=new File(log);
		check(borrowed ? !selectedLog.exists() : selectedLog.isFile() && (fixture ? selectedLog.length()==0 : selectedLog.length()>0),
			"Export fixture parsing leaves an empty selected log; survey preprocessing emits diagnostics");
		final String selected=borrowed ? "" : text(log), captured=diagnostics.toString("UTF-8");
		if(success){
			final int threads=label.equals("default") ? after ? 3 : 8 : label.endsWith("override") ? 3 : 2;
			check(Shared.threads()==threads && ReadWrite.MAX_ZIP_THREADS()==1, "Survey worker budget changes independently of its preserved single decompressor: "+label);
			check(selected.contains("NCRNA_SURVEY_PASS")==after && captured.contains("NCRNA_SURVEY_PASS")==!after,
				"Success diagnostics must reach the selected log after CLI integration: "+label);
		}
		if(label.equals("worker")){
			check(selected.contains("Missing GFF:")==after && captured.contains("Missing GFF:")==!after,
				"Worker diagnostics must use the same selected log as the main invocation");
		}
		for(String suffix:SUFFIXES){check(new File(prefix+suffix).exists()==success, "Only a successful survey creates aggregate tables: "+label+" "+suffix);}
		check(!new File(root, "must_not_open.log").exists(), "Rejected nested config must not open its selected log");
		final ByteBuilder report=new ByteBuilder("case\texpectation\taccepted\tlog_descriptors\tthreads\tzip_threads\texception\n")
			.append(label).tab().append(expectation).tab().append(success).tab().append(descriptors).tab()
			.append(Shared.threads()).tab().append(ReadWrite.MAX_ZIP_THREADS()).tab().append(problem==null ? "none" : problem.getClass().getName()).nl();
		try(Output out=output(root+"/observations.tsv")){out.writer.print(report);}
		try(Output out=output(root+"/captured.stderr")){out.writer.print(captured);}
		System.out.println("NCRNA_SURVEY_CLI_PASS case="+label+" expectation="+expectation);
	}

	/** Native config expansion order is explicit; no generated Java or replacement parser. */
	private static ArrayList<String> command(String root, String input, String prefix, String label){
		assert(input!=null && prefix!=null):"Every CLI case shares the same fixed three-genome input and owns its output prefix";
		final ArrayList<String> base=new ArrayList<String>(Arrays.asList("map="+input+"/map.tsv", "root="+input,
			"expected="+(label.equals("count") ? 4 : 3), "outprefix="+prefix));
		if(!label.equals("default")){base.add("t="+(label.equals("kmg") ? "0.002k" : "2"));}
		if(label.equals("config") || label.equals("config_override") || label.startsWith("nested") || label.equals("too_deep") || label.equals("inner_log")){
			final String inner=root+"/inner.config", outer=root+"/outer.config", deep=root+"/deep.config";
			if(label.equals("inner_log")){base.add("outstream="+root+"/must_not_open.log");}
			write(inner, base);write(outer, Arrays.asList("config="+inner));write(deep, Arrays.asList("config="+outer));
			base.clear();base.add("config="+(label.equals("too_deep") ? deep : label.startsWith("nested") || label.equals("inner_log") ? outer : inner));
			if(label.endsWith("override")){base.add("t=3");}
		}else if(label.equals("false_outputs")){
			base.add("ow=f");base.add("app=f");
		}else if(label.equals("overwrite") || label.equals("hidden_overwrite")){
			base.add("ow=t");if(label.equals("hidden_overwrite")){base.add("ow=f");}
		}else if(label.equals("population")){
			base.add("reads=1");
		}else if(label.equals("parse") || label.equals("borrowed")){
			base.add("unknown=1");
		}else if(label.equals("worker")){
			check(new File(root, "absent_genomes").mkdir(), "Worker failure uses an existing but empty GFF root");
			base.add("root="+root+"/absent_genomes");base.add("loci=t");base.add("t=1");
		}
		return base;
	}
	private static void write(String path, Iterable<String> lines){
		try(Output out=output(path)){for(String line:lines){out.writer.println(line);}}
	}
	/** Keeps the fixture's native TEXT format and fresh-output policy. */
	private static Output output(final String path){
		return new Output(FileFormat.testOutput(path, FileFormat.TEXT, null, true, false, false, false));
	}
	/** Small retained diagnostic text only; scientific inputs remain byte-oriented. */
	private static String text(String path){
		final ByteBuilder result=new ByteBuilder();
		try(Input input=new Input(path)){for(byte[] row; (row=input.reader.nextLine())!=null;){result.append(row).nl();}}
		return result.toString();
	}
	private static int descriptors(String path) throws Exception{
		assert(path!=null):"Count only this case's selected log, never unrelated process descriptors";
		int count=0;
		try(DirectoryStream<Path> entries=Files.newDirectoryStream(new File("/proc/self/fd").toPath())){
			for(Path entry:entries){try{if(Files.readSymbolicLink(entry).toString().equals(path)){count++;}}
			catch(NoSuchFileException closedConcurrently){}}
		}return count;
	}
	private static final class Input implements AutoCloseable{
		Input(String path){reader=ByteFile.makeByteFile(FileFormat.testInput(path, FileFormat.TEXT, null, false, true));}
		@Override public void close(){check(!reader.close(), "Diagnostic reader must report its native I/O error");}
		final ByteFile reader;
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private static final String[] CASES={"explicit", "default", "kmg", "config", "config_override", "nested", "nested_override", "too_deep",
		"inner_log", "false_outputs", "overwrite", "hidden_overwrite", "population", "parse", "count", "worker", "borrowed", "fixture_log"};
	private static final String[] SUFFIXES={".aliases.tsv", ".excluded_regulatory.tsv", ".families.tsv", ".phyla.tsv", ".genomes.tsv", ".summary.tsv"};
}
