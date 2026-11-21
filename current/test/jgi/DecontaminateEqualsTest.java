package test.jgi;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import driver.RenameAndMux;
import fileIO.ReadWrite;
import jgi.DemuxByName;
import jgi.DecontaminateByNormalization;
import jgi.FilterByCoverage;
import jgi.KmerNormalize;
import parse.Parse;
import shared.Shared;

/**
 * Finite first-equals parser baseline/fixed fixture for DecontaminateByNormalization children.
 * @author Jean
 */
public final class DecontaminateEqualsTest{

	/** Requires a new {@code outdir=...} and one {@code case=...}; {@code phase=baseline|fixed} defaults to baseline. */
	public static void main(final String[] args) throws Exception{
		String outdir=null, caseName=null, phase="baseline";
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			final String key=pair[0].toLowerCase(Locale.ROOT);
			if(key.equals("outdir")){outdir=pair[1];}
			else if(key.equals("case")){caseName=pair[1];}
			else if(key.equals("phase")){phase=pair[1];}
			else{throw new IllegalArgumentException("Unknown fixture option: "+arg);}
		}
		require(outdir!=null && !outdir.isEmpty(), "Set outdir=<new case directory>.");
		require(caseName!=null && !caseName.isEmpty(), "Set case=<name>.");
		require(phase.equals("baseline") || phase.equals("fixed"), "Set phase=baseline or phase=fixed.");
		require(outdir.indexOf(',')<0 && outdir.indexOf('=')<0 &&
				outdir.indexOf('\t')<0 && outdir.indexOf('\n')<0 && outdir.indexOf('\r')<0,
				"Fixture root must fit BBTools comma-delimited path arguments and one TSV field.");

		Shared.setThreads(1);
		ReadWrite.setZipThreads(1);
		ReadWrite.FORCE_PIGZ=false;
		ReadWrite.FORCE_BGZIP=false;
		ReadWrite.USE_BGZIP=false;
		ReadWrite.USE_BGZF=false;
		ReadWrite.USE_PIGZ=false;
		ReadWrite.USE_UNPIGZ=false;

		new DecontaminateEqualsTest(Paths.get(outdir).toAbsolutePath().normalize(),
			findCase(caseName), phase).run();
	}

	/** Runs the selected isolated probe and writes exact observations. */
	private void run() throws Exception{
		Files.createDirectory(root);
		report=root.resolve("observations.tsv");
		allowedEntries.add(root.toAbsolutePath().normalize());
		allowedEntries.add(report.toAbsolutePath().normalize());
		write(report, "case\tphase\tfield\texpected\tactual\tmatch\tfailureType\tfailureMessage\n");
		final LinkedHashMap<String, String> expected=new LinkedHashMap<String, String>();
		Throwable failure=null;
		try{
			if(spec.kind.equals("parent")){parent(expected);}
			else if(spec.kind.equals("rename")){rename(expected);}
			else if(spec.kind.equals("demux")){demux(expected);}
			else if(spec.kind.equals("filter")){filter(expected);}
			else if(spec.kind.equals("kmer")){kmer(expected);}
			else{throw new IllegalStateException(spec.kind);}
		}catch(final InvocationTargetException ite){
			failure=ite.getCause();
		}catch(final Throwable caught){
			failure=caught;
		}

		boolean all=true;
		if(spec.failureType==null){
			all=(failure==null);
		}else{
			all=spec.failureType.equals(type(failure)) && spec.failureMessage(root).equals(message(failure));
		}
		for(final Map.Entry<String, String> entry : expected.entrySet()){
			final String actual=observed.containsKey(entry.getKey()) ? observed.get(entry.getKey()) : MISSING;
			final boolean match=entry.getValue().equals(actual);
			all&=match;
			append(report, spec.name+'\t'+phase+'\t'+entry.getKey()+'\t'+entry.getValue()+'\t'+actual+
				'\t'+match+'\t'+type(failure)+'\t'+message(failure)+'\n');
		}
		final List<String> changedInputs=changedInputs();
		final boolean inputMatch=changedInputs.isEmpty();
		all&=inputMatch;
		append(report, spec.name+"\t"+phase+"\t<inputs>\t[]\t"+changedInputs+'\t'+inputMatch+
			'\t'+type(failure)+'\t'+message(failure)+'\n');
		final List<String> extraEntries=unexpectedEntries();
		final boolean extraMatch=extraEntries.isEmpty();
		all&=extraMatch;
		append(report, spec.name+"\t"+phase+"\t<extra-entries>\t[]\t"+extraEntries+'\t'+extraMatch+
			'\t'+type(failure)+'\t'+message(failure)+'\n');
		if(!all){throw new AssertionError("Unexpected first-equals observation; inspect "+report);}
		System.out.println("DECONTAMINATE_EQUALS_PASS phase="+phase+" case="+spec.name);
	}

	/** Exercises only the parent constructor and inspects its retained raw fields. */
	private void parent(final LinkedHashMap<String, String> expected) throws Exception{
		final ArrayList<String> args=new ArrayList<String>();
		args.add("ref="+plainRef());
		args.add("read="+plainRead());
		args.add("mapraw=f");
		if(spec.name.equals("parent_read_path")){
			final Path full=root.resolve("read=a=b=c.fq");
			args.set(1, "read="+full);
			expected.put("readNames", choose(s(root.resolve("read")), s(full)));
		}else if(spec.name.equals("parent_ref_path")){
			final Path full=root.resolve("ref=A=b=c.fa");
			args.set(0, "ref="+full);
			expected.put("refNames", choose(s(root.resolve("ref")), s(full)));
		}else if(spec.name.equals("parent_later_comma_member")){
			final Path later=root.resolve("later=A=b.fq");
			args.set(1, "read="+root.resolve("ordinary.fq")+","+later);
			args.set(0, "ref="+root.resolve("ordinary.fa")+","+root.resolve("later.fa"));
			expected.put("readNames", choose(s(root.resolve("ordinary.fq"))+"|"+s(root.resolve("later")),
				s(root.resolve("ordinary.fq"))+"|"+s(later)));
			expected.put("refNames", s(root.resolve("ordinary.fa"))+"|"+s(root.resolve("later.fa")));
		}else if(spec.name.equals("parent_list_path")){
			final Path readList=root.resolve("reads=list.txt");
			final Path refList=root.resolve("refs=list.txt");
			writeInput(readList, root.resolve("inside-read=a=b.fq")+"\n");
			writeInput(refList, root.resolve("inside-ref=a=b.fa")+"\n");
			args.clear();
			args.add("readnamefile="+readList);
			args.add("refnamefile="+refList);
			expected.put("readNames", choose(s(root.resolve("reads")), s(root.resolve("inside-read=a=b.fq"))));
			expected.put("refNames", choose(s(root.resolve("refs")), s(root.resolve("inside-ref=a=b.fa"))));
		}else if(spec.name.equals("parent_list_contents")){
			final Path readList=root.resolve("reads.txt");
			final Path refList=root.resolve("refs.txt");
			writeInput(readList, root.resolve("inside-read=a=b.fq")+"\n");
			writeInput(refList, root.resolve("inside-ref=A=b.fa")+"\n");
			args.clear();
			args.add("readnamefile="+readList);
			args.add("refnamefile="+refList);
			expected.put("readNames", s(root.resolve("inside-read=a=b.fq")));
			expected.put("refNames", s(root.resolve("inside-ref=A=b.fa")));
		}else if(spec.name.equals("parent_tempdir_path")){
			final Path full=root.resolve("tmp=A=b=c");
			args.add("tempdir="+full);
			expected.put("tempdir", choose(str(p(root.resolve("tmp"))+"/"), str(p(full)+"/")));
		}else if(spec.name.equals("parent_outdir_log_results")){
			final Path out=root.resolve("out=A=b");
			final Path log=root.resolve("log=A=b.txt");
			final Path results=root.resolve("results=A=b.txt");
			args.add("outdir="+out);
			args.add("log="+log);
			args.add("results="+results);
			expected.put("outdir", choose(str(p(root.resolve("out"))+"/"), str(p(out)+"/")));
			expected.put("logName", choose(s(root.resolve("log")), s(log)));
			expected.put("resultsName", choose(s(root.resolve("results")), s(results)));
		}else if(spec.name.startsWith("parent_results_")){
			if(spec.name.endsWith("bare")){args.add("results");}
			else if(spec.name.endsWith("empty")){args.add("results=");}
			else if(spec.name.endsWith("ordinary")){args.add("results=ordinary.txt");}
			else if(spec.name.endsWith("literal_null")){args.add("results=null");}
			expected.put("resultsName", spec.expectedValue);
		}else{throw new IllegalStateException(spec.name);}
		final DecontaminateByNormalization tool=new DecontaminateByNormalization(args.toArray(new String[0]));
		observe("readNames", ordered(READ_NAMES.get(tool)));
		observe("refNames", ordered(REF_NAMES.get(tool)));
		observe("tempdir", scalar(TEMPDIR.get(tool)));
		observe("outdir", scalar(OUTDIR.get(tool)));
		observe("logName", scalar(LOG_NAME.get(tool)));
		observe("resultsName", scalar(RESULTS_NAME.get(tool)));
	}

	/** Exercises RenameAndMux direct filename parsing without processing reads. */
	private void rename(final LinkedHashMap<String, String> expected) throws Exception{
		final Path input=fastq(root.resolve("rename-in.fq"), "ACGTACGT");
		final ArrayList<String> args=new ArrayList<String>();
		args.add("in="+input);
		args.add("out="+root.resolve("rename-out.fq"));
		args.add("ow=t");
		if(spec.name.equals("rename_in_path")){
			final Path full=root.resolve("rename-in=a=b.fq");
			fastq(root.resolve("rename-in"), "TGCATGCA");
			fastq(full, "ACGTACGT");
			args.set(0, "in="+full);
			expected.put("readPaths", choose(s(root.resolve("rename-in")), s(full)));
		}else if(spec.name.equals("rename_out_path")){
			final Path full=root.resolve("rename-out=A=b.fq");
			args.set(1, "out="+full);
			expected.put("out1", choose(s(root.resolve("rename-out")), s(full)));
		}else if(spec.name.startsWith("rename_out_")){
			if(spec.name.endsWith("omitted")){args.remove(1);}
			else if(spec.name.endsWith("bare")){args.set(1, "out");}
			else if(spec.name.endsWith("empty")){args.set(1, "out=");}
			else if(spec.name.endsWith("ordinary")){args.set(1, "out=ordinary.fq");}
			else if(spec.name.endsWith("literal_null")){args.set(1, "out=null");}
			if(spec.failureType==null){expected.put("out1", spec.expectedValue);}
		}else{throw new IllegalStateException(spec.name);}
		final RenameAndMux tool=new RenameAndMux(args.toArray(new String[0]));
		observe("readPaths", ordered(RENAME_READ_PATHS.get(tool)));
		observe("out1", scalar(RENAME_OUT1.get(tool)));
	}

	/** Exercises legacy DemuxByName direct input, pattern, and literal-name parsing. */
	private void demux(final LinkedHashMap<String, String> expected) throws Exception{
		final Path input=fastq(root.resolve("demux-in.fq"), "ACGTACGT");
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(
			"in="+input, "out="+root.resolve("%_demux.fq"), "names=alpha", "literalnames=t", "ow=t"));
		if(spec.name.equals("demux_in_path")){
			final Path full=root.resolve("demux-in=a=b.fq");
			fastq(root.resolve("demux-in"), "TGCATGCA");
			fastq(full, "ACGTACGT");
			args.set(0, "in="+full);
			expected.put("in1", choose(s(root.resolve("demux-in")), s(full)));
		}else if(spec.name.equals("demux_out_pattern")){
			final Path full=root.resolve("%_demux=A=b.fq");
			args.set(1, "out="+full);
			expected.put("out1", choose(s(root.resolve("%_demux")), s(full)));
		}else if(spec.name.equals("demux_literal_names")){
			args.set(2, "names=alpha,beta=A=b_");
			expected.put("names", choose(str("alpha")+"|"+str("beta"), str("alpha")+"|"+str("beta=A=b_")));
		}else if(spec.name.equals("demux_default_namefile")){
			final Path manifest=root.resolve("names.txt");
			writeInput(manifest, "alpha=A=b\n");
			args.set(2, "names="+manifest);
			args.remove(3);
			expected.put("literalNames", str("false"));
			expected.put("names", str("alpha=A=b"));
		}else if(spec.name.startsWith("demux_names_")){
			args.add("delimiter=whitespace");
			if(spec.name.endsWith("omitted")){args.remove(2);}
			else if(spec.name.endsWith("bare")){args.set(2, "names");}
			else if(spec.name.endsWith("empty")){args.set(2, "names=");}
			else if(spec.name.endsWith("ordinary")){args.set(2, "names=ordinary");}
			else if(spec.name.endsWith("literal_null")){args.set(2, "names=null");}
			expected.put("names", spec.expectedValue);
		}else{throw new IllegalStateException(spec.name);}
		final DemuxByName tool=new DemuxByName(args.toArray(new String[0]));
		observe("in1", scalar(DEMUX_IN1.get(tool)));
		observe("out1", scalar(DEMUX_OUT1.get(tool)));
		observe("names", sorted(DEMUX_NAMES.get(tool)));
		observe("literalNames", scalar(DEMUX_LITERAL_NAMES.get(tool)));
	}

	/** Exercises FilterByCoverage input, generated output/log/cov fields, and trim null controls. */
	private void filter(final LinkedHashMap<String, String> expected) throws Exception{
		final Path input=fasta(root.resolve("filter-in.fa"));
		final Path cov0=cov(root.resolve("cov0.txt"));
		final Path cov1=cov(root.resolve("cov1.txt"));
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(
			"in="+input, "cov0="+cov0, "cov1="+cov1, "out="+root.resolve("clean.fa"),
			"outd="+root.resolve("dirty.fa"), "log="+root.resolve("filter.log"), "ow=t"));
		if(spec.name.equals("filter_in_path")){
			final Path full=fasta(root.resolve("filter-in=A=b.fa"));
			fasta(root.resolve("filter-in"));
			args.set(0, "in="+full);
			expected.put("in1", choose(s(root.resolve("filter-in")), s(full)));
		}else if(spec.name.equals("filter_cov0_path")){
			final Path full=root.resolve("cov0=A=b.txt");
			cov(root.resolve("cov0"));
			cov(full);
			args.set(1, "cov0="+full);
			expected.put("covStatsBefore", choose(s(root.resolve("cov0")), s(full)));
		}else if(spec.name.equals("filter_cov1_path")){
			final Path full=root.resolve("cov1=A=b.txt");
			cov(root.resolve("cov1"));
			cov(full);
			args.set(2, "cov1="+full);
			expected.put("covStatsAfter", choose(s(root.resolve("cov1")), s(full)));
		}else if(spec.name.equals("filter_log_path")){
			final Path full=root.resolve("filter=A=b.log");
			args.set(5, "log="+full);
			expected.put("logfile", choose(s(root.resolve("filter")), s(full)));
		}else if(spec.name.equals("filter_out_path")){
			final Path full=root.resolve("clean=A=b.fa");
			args.set(3, "out="+full);
			expected.put("outclean", choose(s(root.resolve("clean")), s(full)));
		}else if(spec.name.equals("filter_outd_path")){
			final Path full=root.resolve("dirty=A=b.fa");
			args.set(4, "outd="+full);
			expected.put("outdirty", choose(s(root.resolve("dirty")), s(full)));
		}else if(spec.name.startsWith("filter_trim_")){
			if(spec.name.endsWith("bare")){args.add("trim");}
			else if(spec.name.endsWith("empty")){args.add("trim=");}
			else if(spec.name.endsWith("ordinary")){args.add("trim=7");}
			expected.put("trimEnds", spec.expectedValue);
		}else if(spec.name.equals("filter_log_literal_null")){
			args.set(5, "log=null");
			expected.put("logfile", spec.expectedValue);
		}else{throw new IllegalStateException(spec.name);}
		final FilterByCoverage tool=new FilterByCoverage(args.toArray(new String[0]));
		observe("in1", scalar(FILTER_IN1.get(tool)));
		observe("covStatsBefore", scalar(FILTER_COV0.get(tool)));
		observe("covStatsAfter", scalar(FILTER_COV1.get(tool)));
		observe("outclean", scalar(FILTER_OUTCLEAN.get(tool)));
		observe("outdirty", scalar(FILTER_OUTDIRTY.get(tool)));
		observe("logfile", scalar(FILTER_LOG.get(tool)));
		observe("trimEnds", scalar(FILTER_TRIM.get(tool)));
	}

	/** Runs KmerNormalize.main once on a bounded one-input, non-# workload. */
	private void kmer(final LinkedHashMap<String, String> expected) throws Exception{
		final Path input=fastq(root.resolve("kmer-in.fq"), "ACGTACGT");
		// k=4 stays above KCountArray's two-array minimum after KmerNormalize caps cells to 4^k.
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+input,
			"out="+root.resolve("kmer-out.fq"), "k=4", "cells=1024", "bits=2", "hashes=1",
			"passes=1", "t=1", "keepall=t", "ecc=f", "prefilter=f", "rename=f", "minprob=0",
			"ow=t", "qin=33", "pigz=f", "unpigz=f"));
		Path expectedOut=root.resolve("kmer-out.fq");
		String expectedSequence=str("ACGTACGT");
		if(spec.name.equals("kmer_in_path")){
			final Path full=fastq(root.resolve("kmer-in=A=b.fq"), "TTTTCCCC");
			fastq(root.resolve("kmer-in"), "ACGTACGT");
			args.set(0, "in="+full);
			expectedSequence=choose(str("ACGTACGT"), str("TTTTCCCC"));
		}else if(spec.name.equals("kmer_out_path")){
			final Path full=root.resolve("kmer-out=A=b.fq");
			args.set(1, "out="+full);
			expectedOut=phase.equals("baseline") ? root.resolve("kmer-out") : full;
		}else if(spec.name.startsWith("kmer_tmpdir_")){
			if(spec.name.endsWith("bare")){args.add("tmpdir");}
			else if(spec.name.endsWith("empty")){args.add("tmpdir=");}
			else if(spec.name.endsWith("ordinary")){args.add("tmpdir=ordinaryTmp");}
			else if(spec.name.endsWith("literal_null")){args.add("tmpdir=null");}
			expected.put("Shared.tmpdir", spec.name.endsWith("omitted") ?
				scalar(Shared.tmpdir()) : spec.expectedValue);
			expected.put("KmerNormalize.TMPDIR", scalar(KMER_TMPDIR.get(null)));
		}else{throw new IllegalStateException(spec.name);}
		allowOutput(expectedOut);
		KmerNormalize.main(args.toArray(new String[0]));
		observe("outputSequences", sorted(readSequences(expectedOut)));
		expected.put("outputSequences", expectedSequence);
		if(spec.name.startsWith("kmer_tmpdir_")){
			observe("Shared.tmpdir", scalar(Shared.tmpdir()));
			observe("KmerNormalize.TMPDIR", scalar(KMER_TMPDIR.get(null)));
		}
	}

	private String plainRead(){return p(root.resolve("plain-read.fq"));}
	private String plainRef(){return p(root.resolve("plain-ref.fa"));}

	private Path fastq(final Path path, final String sequence) throws Exception{
		writeInput(path, "@r\n"+sequence+"\n+\nIIIIIIII\n");
		return path;
	}

	private Path fasta(final Path path) throws Exception{
		writeInput(path, ">ref\nACGTACGTACGT\n");
		return path;
	}

	private Path cov(final Path path) throws Exception{
		writeInput(path, "#ID\tAvg_fold\tLength\tRef_GC\tCovered_percent\tCovered_bases\tPlus_reads\tMinus_reads\tRead_GC\n"+
			"ref\t9\t12\t0.5\t100\t12\t2\t2\t0.5\n");
		return path;
	}

	private List<String> readSequences(final Path path) throws Exception{
		if(!Files.exists(path)){return Collections.emptyList();}
		final ArrayList<String> list=new ArrayList<String>();
		InputStream is=Files.newInputStream(path);
		try{
			final BufferedReader reader=new BufferedReader(new InputStreamReader(is, StandardCharsets.US_ASCII));
			for(int line=0; ; line++){
				final String s=reader.readLine();
				if(s==null){break;}
				if((line&3)==1){list.add(s);}
			}
		}finally{is.close();}
		return list;
	}

	private List<String> changedInputs() throws Exception{
		final ArrayList<String> changed=new ArrayList<String>();
		for(final Map.Entry<Path, byte[]> entry : inputBytes.entrySet()){
			final Path path=entry.getKey();
			if(!Files.exists(path) || !Arrays.equals(entry.getValue(), Files.readAllBytes(path))){
				changed.add(rel(path));
			}
		}
		Collections.sort(changed);
		return changed;
	}

	private List<String> unexpectedEntries() throws Exception{
		final ArrayList<String> extra=new ArrayList<String>();
		collectUnexpected(root, extra);
		Collections.sort(extra);
		return extra;
	}

	private void collectUnexpected(final Path path, final ArrayList<String> extra) throws Exception{
		final Path abs=path.toAbsolutePath().normalize();
		if(!allowedEntries.contains(abs)){extra.add(rel(abs)+(Files.isDirectory(abs) ? "/" : ""));}
		if(Files.isDirectory(abs)){
			try(final DirectoryStream<Path> stream=Files.newDirectoryStream(abs)){
				for(final Path child : stream){collectUnexpected(child, extra);}
			}
		}
	}

	private void writeInput(final Path path, final String text) throws Exception{
		allowParents(path);
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		final Path abs=path.toAbsolutePath().normalize();
		allowedEntries.add(abs);
		inputBytes.put(abs, Files.readAllBytes(abs));
	}

	private void allowOutput(final Path path) throws Exception{
		allowParents(path);
		allowedEntries.add(path.toAbsolutePath().normalize());
	}

	private void allowParents(final Path path) throws Exception{
		final Path parent=path.toAbsolutePath().normalize().getParent();
		if(parent!=null){Files.createDirectories(parent); allowAncestry(parent);}
	}

	private void allowAncestry(final Path path){
		final Path abs=path.toAbsolutePath().normalize();
		allowedEntries.add(abs);
		if(!abs.equals(root) && abs.getParent()!=null && abs.startsWith(root)){allowAncestry(abs.getParent());}
	}

	private static CaseSpec findCase(final String name){
		for(final CaseSpec spec : CASES){if(spec.name.equals(name)){return spec;}}
		throw new IllegalArgumentException("Unknown case: "+name);
	}

	private void observe(final String field, final String value){observed.put(field, value);}

	private String choose(final String baseline, final String fixed){return phase.equals("baseline") ? baseline : fixed;}
	private static String p(final Path path){return path.toAbsolutePath().normalize().toString();}
	private static String s(final Path path){return str(p(path));}
	private static String str(final String s){return "<str:"+s+">";}
	private static String scalar(final Object o){return o==null ? NULL : str(String.valueOf(o));}

	private static String ordered(final Object o){
		if(o==null){return NULL;}
		final ArrayList<String> list=new ArrayList<String>();
		for(final Object x : (Iterable<?>)o){list.add(scalar(x));}
		return join(list, "|");
	}

	private static String sorted(final Object o){
		if(o==null){return NULL;}
		final ArrayList<String> list=new ArrayList<String>();
		for(final Object x : (Iterable<?>)o){list.add(scalar(x));}
		Collections.sort(list);
		return join(list, "|");
	}

	private String rel(final Path path){
		final Path abs=path.toAbsolutePath().normalize();
		return abs.startsWith(root) ? root.relativize(abs).toString() : abs.toString();
	}

	private static String type(final Throwable failure){return failure==null ? "none" : failure.getClass().getName();}

	private static String message(final Throwable failure){
		return failure==null ? "" : String.valueOf(failure.getMessage()).replace('\t', ' ').
			replace('\n', ' ').replace('\r', ' ');
	}

	private static String join(final List<String> list, final String delimiter){
		final StringBuilder sb=new StringBuilder();
		for(final String s : list){
			if(sb.length()>0){sb.append(delimiter);}
			sb.append(s);
		}
		return sb.toString();
	}

	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	private static void append(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);
	}

	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	private static Field field(final Class<?> clazz, final String name){
		try{
			final Field field=clazz.getDeclaredField(name);
			field.setAccessible(true);
			return field;
		}catch(final Exception e){
			throw new RuntimeException(e);
		}
	}

	private DecontaminateEqualsTest(final Path root_, final CaseSpec spec_, final String phase_){
		root=root_;
		spec=spec_;
		phase=phase_;
	}

	private final Path root;
	private final CaseSpec spec;
	private final String phase;
	private Path report;
	private final LinkedHashMap<String, String> observed=new LinkedHashMap<String, String>();
	private final HashSet<Path> allowedEntries=new HashSet<Path>();
	private final HashMap<Path, byte[]> inputBytes=new HashMap<Path, byte[]>();

	private static final String NULL="<java-null>";
	private static final String MISSING="<absent>";
	private static final String MISSING_OUT="java.lang.RuntimeException";
	private static final String MISSING_OUT_MESSAGE="Error - output destination is required.";

	private static final CaseSpec[] CASES={
		new CaseSpec("parent_read_path", "parent", null, null),
		new CaseSpec("parent_ref_path", "parent", null, null),
		new CaseSpec("parent_later_comma_member", "parent", null, null),
		new CaseSpec("parent_list_path", "parent", null, null),
		new CaseSpec("parent_list_contents", "parent", null, null),
		new CaseSpec("parent_tempdir_path", "parent", null, null),
		new CaseSpec("parent_outdir_log_results", "parent", null, null),
		new CaseSpec("parent_results_omitted", "parent", null, str("results.txt")),
		new CaseSpec("parent_results_bare", "parent", null, NULL),
		new CaseSpec("parent_results_empty", "parent", null, NULL),
		new CaseSpec("parent_results_ordinary", "parent", null, str("ordinary.txt")),
		new CaseSpec("parent_results_literal_null", "parent", null, NULL),
		new CaseSpec("rename_in_path", "rename", null, null),
		new CaseSpec("rename_out_path", "rename", null, null),
		new CaseSpec("rename_out_omitted", "rename", MISSING_OUT, MISSING_OUT_MESSAGE),
		new CaseSpec("rename_out_bare", "rename", MISSING_OUT, MISSING_OUT_MESSAGE),
		new CaseSpec("rename_out_empty", "rename", MISSING_OUT, MISSING_OUT_MESSAGE),
		new CaseSpec("rename_out_ordinary", "rename", null, str("ordinary.fq")),
		new CaseSpec("rename_out_literal_null", "rename", MISSING_OUT, MISSING_OUT_MESSAGE),
		new CaseSpec("demux_in_path", "demux", null, null),
		new CaseSpec("demux_out_pattern", "demux", null, null),
		new CaseSpec("demux_literal_names", "demux", null, null),
		new CaseSpec("demux_default_namefile", "demux", null, null),
		new CaseSpec("demux_names_omitted", "demux", null, ""),
		new CaseSpec("demux_names_bare", "demux", null, ""),
		new CaseSpec("demux_names_empty", "demux", null, ""),
		new CaseSpec("demux_names_ordinary", "demux", null, str("ordinary")),
		new CaseSpec("demux_names_literal_null", "demux", null, ""),
		new CaseSpec("filter_in_path", "filter", null, null),
		new CaseSpec("filter_cov0_path", "filter", null, null),
		new CaseSpec("filter_cov1_path", "filter", null, null),
		new CaseSpec("filter_log_path", "filter", null, null),
		new CaseSpec("filter_out_path", "filter", null, null),
		new CaseSpec("filter_outd_path", "filter", null, null),
		new CaseSpec("filter_trim_omitted", "filter", null, str("0")),
		new CaseSpec("filter_trim_bare", "filter", null, str("100")),
		new CaseSpec("filter_trim_empty", "filter", null, str("100")),
		new CaseSpec("filter_trim_ordinary", "filter", null, str("7")),
		new CaseSpec("filter_log_literal_null", "filter", null, NULL),
		new CaseSpec("kmer_in_path", "kmer", null, null),
		new CaseSpec("kmer_out_path", "kmer", null, null),
		new CaseSpec("kmer_tmpdir_omitted", "kmer", null, null),
		new CaseSpec("kmer_tmpdir_bare", "kmer", null, NULL),
		new CaseSpec("kmer_tmpdir_empty", "kmer", null, NULL),
		new CaseSpec("kmer_tmpdir_ordinary", "kmer", null, str("ordinaryTmp/")),
		new CaseSpec("kmer_tmpdir_literal_null", "kmer", null, NULL)
	};

	private static final Field READ_NAMES=field(DecontaminateByNormalization.class, "readNames");
	private static final Field REF_NAMES=field(DecontaminateByNormalization.class, "refNames");
	private static final Field TEMPDIR=field(DecontaminateByNormalization.class, "tempdir");
	private static final Field OUTDIR=field(DecontaminateByNormalization.class, "outdir");
	private static final Field LOG_NAME=field(DecontaminateByNormalization.class, "logName");
	private static final Field RESULTS_NAME=field(DecontaminateByNormalization.class, "resultsName");
	private static final Field RENAME_READ_PATHS=field(RenameAndMux.class, "readPaths");
	private static final Field RENAME_OUT1=field(RenameAndMux.class, "out1");
	private static final Field DEMUX_IN1=field(DemuxByName.class, "in1");
	private static final Field DEMUX_OUT1=field(DemuxByName.class, "out1");
	private static final Field DEMUX_NAMES=field(DemuxByName.class, "names");
	private static final Field DEMUX_LITERAL_NAMES=field(DemuxByName.class, "literalNames");
	private static final Field FILTER_IN1=field(FilterByCoverage.class, "in1");
	private static final Field FILTER_COV0=field(FilterByCoverage.class, "covStatsBefore");
	private static final Field FILTER_COV1=field(FilterByCoverage.class, "covStatsAfter");
	private static final Field FILTER_OUTCLEAN=field(FilterByCoverage.class, "outclean");
	private static final Field FILTER_OUTDIRTY=field(FilterByCoverage.class, "outdirty");
	private static final Field FILTER_LOG=field(FilterByCoverage.class, "logfile");
	private static final Field FILTER_TRIM=field(FilterByCoverage.class, "trimEnds");
	private static final Field KMER_TMPDIR=field(KmerNormalize.class, "TMPDIR");

	private static final class CaseSpec{
		CaseSpec(final String name_, final String kind_, final String failureType_, final String expectedValue_){
			name=name_;
			kind=kind_;
			failureType=failureType_;
			expectedValue=expectedValue_;
		}
		String failureMessage(final Path root){return expectedValue==null ? null : expectedValue.replace("%ROOT%", root.toString());}
		final String name;
		final String kind;
		final String failureType;
		final String expectedValue;
	}
}
