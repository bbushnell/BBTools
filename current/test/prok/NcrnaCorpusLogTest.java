package test.prok;

import java.io.File;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import fileIO.FileFormat;
import parse.PreParser;
import prok.NcrnaAnnotatedCorpus;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Measures native selected-log ownership after one real extractor invocation.
 * Run one case per JVM: process exit reclaims the deliberately observed baseline
 * log descriptor. No reflection, replacement writer or injected failure is used.
 * @author Raiden */
public final class NcrnaCorpusLogTest{

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaCorpusLogTest.class, false);
		try{run(pp.args);}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(String[] args) throws Exception{
		String root=null, roster=null, label=null, expect=null;
		for(String arg:args){
			if(arg.startsWith("outdir=")){root=arg.substring(7);}
			else if(arg.startsWith("roster=")){roster=arg.substring(7);}
			else if(arg.startsWith("case=")){label=arg.substring(5);}
			else if(arg.startsWith("expect=")){expect=arg.substring(7);}
			else{throw new IllegalArgumentException("Unknown selected-log fixture argument: "+arg);}
		}
		check(root!=null && roster!=null && ("before".equals(expect) || "after".equals(expect)),
			"Declare the fixed three-genome fifteen-locus roster, fresh outdir and before/after expectation");
		check("parse".equals(label) || "roster".equals(label) || "body".equals(label)
			|| "success".equals(label) || "borrowed".equals(label), "Select one declared lifecycle case per JVM");
		check(new File(root).mkdir(), "Each selected-log case must own a fresh output directory");
		final boolean success=label.equals("success"), borrowed=label.equals("borrowed");
		final String log=new File(root, "selected.log").getCanonicalPath(), prefix=root+"/corpus";
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("roster="+roster,
			"outprefix="+prefix, "family=RF00174", "expectedgenomes="+(label.equals("roster") ? 4 : 3),
			"expectedloci="+(label.equals("body") ? 14 : 15), "flank=2", "t=1", "gccontig=t"));
		if(!borrowed){command.add("outstream="+log);}
		if(label.equals("parse") || borrowed){command.add("unknown=1");}
		IllegalArgumentException problem=null;
		try{NcrnaAnnotatedCorpus.main(command.toArray(new String[0]));}
		catch(IllegalArgumentException failure){problem=failure;}
		if(success){check(problem==null, "The declared complete corpus must succeed: "+problem);}
		else{
			final String message=label.equals("roster") ? "Corpus membership must match the complete frozen survey, even in a pilot"
				: label.equals("body") ? "Extracted raw loci must reconcile with the survey: 15 expected=14" : "unknown";
			check(problem!=null && message.equals(problem.getMessage()),
				"Selected-log cleanup must preserve the exact original extractor failure: "+label+" "+problem);
		}
		final int descriptors=openDescriptors(log);
		final int expected=expect.equals("before") && !success && !borrowed ? 1 : 0;
		check(descriptors==expected, "Selected log descriptor count after "+label+": actual="+descriptors+", expected="+expected);
		check(borrowed ? !new File(log).exists() : new File(log).length()>0,
			"Only an explicitly selected log must be created and populated by native preprocessing");
		System.err.println("NCRNA_CORPUS_BORROWED_STDERR_OPEN");
		check(!System.err.checkError(), "Extractor cleanup must leave borrowed stderr usable after success or failure");
		for(String suffix:SUFFIXES){
			check(new File(prefix+suffix).exists()==(success || label.equals("body")),
				"Only the completed extraction and output-stage failure may create corpus outputs: "+label+" "+suffix);
		}
		final ByteBuilder report=new ByteBuilder("case\texpectation\topen_log_descriptors\toutcome\n")
			.append(label).tab().append(expect).tab().append(descriptors).tab().append(success ? "success" : problem.getMessage()).nl();
		try(Output output=new Output(FileFormat.testOutput(root+"/observations.tsv", FileFormat.TEXT, null, true, false, false, false))){
			output.writer.print(report);
		}
		System.out.println("NCRNA_CORPUS_LOG_PASS case="+label+" expectation="+expect+" descriptors="+descriptors);
	}

	/** Linux metadata only; exact canonical target excludes unrelated process descriptors. */
	private static int openDescriptors(String path) throws Exception{
		int count=0;
		try(DirectoryStream<Path> entries=Files.newDirectoryStream(new File("/proc/self/fd").toPath())){
			for(Path entry:entries){
				try{if(Files.readSymbolicLink(entry).toString().equals(path)){count++;}}
				catch(NoSuchFileException closedConcurrently){}
			}
		}
		return count;
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private static final String[] SUFFIXES={".all.fa.gz", ".full.fa.gz", ".contexts.fa.gz", ".loci.tsv.gz", ".genomes.tsv", ".summary.tsv"};
}
