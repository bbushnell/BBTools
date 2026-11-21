package test.jgi;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import fileIO.ReadWrite;
import jgi.DemuxByName;
import jgi.DecontaminateByNormalization;
import parse.Parse;
import shared.Shared;

/**
 * Source-level DecontaminateByNormalization naming fixture for RenameAndMux to DemuxByName.
 * @author Jean
 */
public final class DecontaminateNamingTest{

	/** Requires a new {@code outdir=...} and {@code case=...}; {@code phase=baseline|fixed} defaults to baseline. */
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
		require(phase.equals("baseline") || phase.equals("fixed"),
			"Set phase=baseline or phase=fixed.");
		require(outdir.indexOf(',')<0 && outdir.indexOf('=')<0 &&
				outdir.indexOf('\t')<0 && outdir.indexOf('\n')<0 && outdir.indexOf('\r')<0,
				"Fixture paths must fit comma-delimited Decontaminate arguments and one TSV row.");

		Shared.setThreads(1);
		ReadWrite.setZipThreads(1);
		ReadWrite.FORCE_PIGZ=false;
		ReadWrite.FORCE_BGZIP=false;
		ReadWrite.USE_BGZIP=false;
		ReadWrite.USE_BGZF=false;
		ReadWrite.USE_PIGZ=false;
		ReadWrite.USE_UNPIGZ=false;

		new DecontaminateNamingTest(Paths.get(outdir).toAbsolutePath().normalize(),
			findCase(phase, caseName)).run();
	}

	/** Writes inputs, invokes real parent stages, and records memberships; legacy name-file coverage uses public DemuxByName.main. */
	private void run() throws Exception{
		Files.createDirectory(root);
		final Path temp=root.resolve("tmp");
		Files.createDirectory(temp);
		final Path report=root.resolve("observations.tsv");
		write(report, "case\toutput\texpected\tactual\tmatch\tfailurePhase\tfailureType\tfailureMessage\n");
		final ArrayList<String> reads=new ArrayList<String>();
		final ArrayList<String> refs=new ArrayList<String>();
		final LinkedHashMap<Path, byte[]> inputBytes=new LinkedHashMap<Path, byte[]>();
		for(final Map.Entry<String, String> entry : spec.cwdFiles.entrySet()){
			final Path file=Paths.get(entry.getKey());
			write(file, entry.getValue());
			inputBytes.put(file.toAbsolutePath().normalize(), Files.readAllBytes(file));
		}
		for(final Library library : spec.libraries){
			final Path read=root.resolve(library.relativeRead);
			final Path ref=root.resolve(library.refName);
			if(read.getParent()!=null){Files.createDirectories(read.getParent());}
			write(read, fastq(library.sequences));
			write(ref, ">"+library.name+"\nACGTACGTACGT\n");
			inputBytes.put(read, Files.readAllBytes(read));
			inputBytes.put(ref, Files.readAllBytes(ref));
			reads.add(read.toString());
			refs.add(ref.toString());
		}
		if(spec.publicDemuxNamefile){
			final Path manifest=root.resolve("names.txt");
			write(manifest, "a\nb\n");
			inputBytes.put(manifest, Files.readAllBytes(manifest));
		}

		String failurePhase="none";
		Throwable failure=null;
		try{
			failurePhase="constructor";
			final DecontaminateByNormalization tool=new DecontaminateByNormalization(
				args(reads, refs, temp));
			final ArrayList<String> parsedReadNames=readNames(tool);
			final Path merged=temp.resolve("merged.fq");

			failurePhase="renameAndMux_MT";
			RENAME.invoke(tool, parsedReadNames, merged.toString());

			failurePhase="demux";
			if(spec.publicDemuxNamefile){
				DemuxByName.main(new String[]{"in="+merged, "out="+temp.resolve("%_demuxed.fq.gz"),
					"names="+root.resolve("names.txt"), "t=1", "pigz=f", "unpigz=f"});
			}else{
				DEMUX.invoke(tool, merged.toString(), parsedReadNames);
			}
			failurePhase="none";
		}catch(final InvocationTargetException ite){
			failure=ite.getCause();
		}catch(final Throwable caught){
			failure=caught;
		}

		final boolean failureExpected=(spec.expectedFailurePhase!=null);
		final boolean failureMatch=(!failureExpected ? failure==null :
			spec.expectedFailurePhase.equals(failurePhase) &&
			spec.expectedFailureType.equals(type(failure)) &&
			spec.expectedFailureMessage(root).equals(message(failure)));
		boolean allMatch=failureMatch;
		for(final Map.Entry<String, List<String>> entry : spec.expected.entrySet()){
			final Path output=temp.resolve(entry.getKey()+"_demuxed.fq.gz");
			final List<String> expected=sorted(entry.getValue());
			final List<String> actual=readSequences(output);
			final boolean match=expected.equals(actual);
			allMatch&=match;
			append(report, spec.name+'\t'+entry.getKey()+'\t'+expected+'\t'+actual+'\t'+match+'\t'+
				failurePhase+'\t'+type(failure)+'\t'+message(failure)+'\n');
		}
		final List<String> extras=extraOutputs(temp, spec.expected, failureExpected);
		final boolean extraMatch=extras.isEmpty();
		allMatch&=extraMatch;
		append(report, spec.name+"\t<extra>\t[]\t"+extras+"\t"+extraMatch+'\t'+failurePhase+
			'\t'+type(failure)+'\t'+message(failure)+'\n');
		final List<String> changedInputs=changedInputs(inputBytes, root);
		final boolean inputMatch=changedInputs.isEmpty();
		allMatch&=inputMatch;
		append(report, spec.name+"\t<inputs>\t[]\t"+changedInputs+"\t"+inputMatch+'\t'+
			failurePhase+'\t'+type(failure)+'\t'+message(failure)+'\n');
		if(!allMatch){
			throw new AssertionError("Unexpected Decontaminate naming observation; inspect "+report);
		}
		System.out.println("DECONTAMINATE_NAMING_PASS case="+spec.name);
	}

	/** Looks up the one requested literal phase and case. */
	private static CaseSpec findCase(final String phase, final String name){
		for(final CaseSpec spec : CASES){
			if(spec.phase.equals(phase) && spec.name.equals(name)){return spec;}
		}
		throw new IllegalArgumentException("Unknown case: phase="+phase+", case="+name);
	}

	/** Fetches the constructor-parsed input list to exercise the real wrapper's argument path. */
	@SuppressWarnings("unchecked")
	private static ArrayList<String> readNames(final DecontaminateByNormalization tool) throws Exception{
		return (ArrayList<String>)READ_NAMES.get(tool);
	}

	/** Writes one four-line FASTQ record per supplied marker sequence. */
	private static String fastq(final String[] sequences){
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<sequences.length; i++){
			sb.append("@read").append(i).append('\n');
			sb.append(sequences[i]).append('\n');
			sb.append("+\n");
			for(int j=0; j<sequences[i].length(); j++){sb.append('I');}
			sb.append('\n');
		}
		return sb.toString();
	}

	/** Reads sequence lines from a possibly gzip-compressed FASTQ output. */
	private static List<String> readSequences(final Path path) throws Exception{
		if(!Files.exists(path)){return Collections.emptyList();}
		final ArrayList<String> list=new ArrayList<String>();
		InputStream is=Files.newInputStream(path);
		if(path.toString().endsWith(".gz")){is=new GZIPInputStream(is);}
		try{
			final BufferedReader reader=new BufferedReader(new InputStreamReader(is, StandardCharsets.US_ASCII));
			for(int line=0; ; line++){
				final String s=reader.readLine();
				if(s==null){break;}
				if((line&3)==1){list.add(s);}
			}
		}finally{is.close();}
		Collections.sort(list);
		return list;
	}

	/** Lists unexpected demux or, for constructor-failure cases, any temp outputs by filename. */
	private static List<String> extraOutputs(final Path temp, final Map<String, List<String>> expected,
			final boolean allTempFiles) throws Exception{
		final ArrayList<String> list=new ArrayList<String>();
		try(final DirectoryStream<Path> stream=Files.newDirectoryStream(temp, allTempFiles ? "*" : "*_demuxed.fq.gz")){
			for(final Path path : stream){
				if(Files.isDirectory(path)){continue;}
				final String name=path.getFileName().toString();
				if(allTempFiles){list.add(name);}
				else{
					final String core=name.substring(0, name.length()-"_demuxed.fq.gz".length());
					if(!expected.containsKey(core)){list.add(name);}
				}
			}
		}
		Collections.sort(list);
		return list;
	}

	/** Lists generated read/ref inputs whose exact bytes changed. */
	private static List<String> changedInputs(final Map<Path, byte[]> inputs, final Path root) throws Exception{
		final ArrayList<String> list=new ArrayList<String>();
		for(final Map.Entry<Path, byte[]> entry : inputs.entrySet()){
			final Path path=entry.getKey();
			if(!Files.exists(path) || !Arrays.equals(entry.getValue(), Files.readAllBytes(path))){
				list.add(path.startsWith(root) ? root.relativize(path).toString() : path.toString());
			}
		}
		return list;
	}

	/** Builds the real parent constructor arguments. */
	private String[] args(final ArrayList<String> reads, final ArrayList<String> refs, final Path temp){
		final ArrayList<String> list=new ArrayList<String>();
		list.add("read="+join(reads));
		list.add("ref="+join(refs));
		list.add("tempdir="+temp);
		list.add("t=1");
		list.add("pigz=f");
		list.add("unpigz=f");
		for(final String arg : spec.extraArgs){list.add(arg);}
		return list.toArray(new String[0]);
	}

	/** Writes a new ASCII file. */
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	/** Appends one ASCII report row. */
	private static void append(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);
	}

	/** Joins paths for legacy comma-list arguments. */
	private static String join(final ArrayList<String> list){
		final StringBuilder sb=new StringBuilder();
		for(final String s : list){
			if(sb.length()>0){sb.append(',');}
			sb.append(s);
		}
		return sb.toString();
	}

	/** Returns a sorted copy to avoid depending on stream flush order. */
	private static List<String> sorted(final List<String> list){
		final ArrayList<String> copy=new ArrayList<String>(list);
		Collections.sort(copy);
		return copy;
	}

	/** Formats an absent or present exception class for the report. */
	private static String type(final Throwable failure){
		return failure==null ? "none" : failure.getClass().getName();
	}

	/** Formats an absent or present exception message for the report. */
	private static String message(final Throwable failure){
		return failure==null ? "" : String.valueOf(failure.getMessage()).replace('\t', ' ').
			replace('\n', ' ').replace('\r', ' ');
	}

	/** Rejects fixture misuse. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	private DecontaminateNamingTest(final Path root_, final CaseSpec spec_){
		root=root_;
		spec=spec_;
	}

	private final Path root;
	private final CaseSpec spec;

	private static final String A0="ACGTACGA";
	private static final String A1="CCTTAACC";
	private static final String B0="GGTTGGTA";

	private static final CaseSpec[] CASES={
		baseline("control",
			new Library[]{new Library("a", "a.fq", A0, A1), new Library("b", "b.fq", B0)},
			map(new String[]{"a", "b"}, new String[][]{{A0, A1}, {B0}})),
		baseline("prefix",
			new Library[]{new Library("a", "a.fq", A0, A1), new Library("a_1", "a_1.fq", B0)},
			map(new String[]{"a", "a_1"}, new String[][]{{A0}, {A1, B0}})),
		baseline("duplicate",
			new Library[]{new Library("left", "left/a.fq", A0, A1), new Library("right", "right/a.fq", B0)},
			map(new String[]{"a"}, new String[][]{{A0, A1, B0}})),
		fixed("control",
			new Library[]{new Library("a", "a.fq", A0, A1), new Library("b", "b.fq", B0)},
			map(new String[]{"a", "b"}, new String[][]{{A0, A1}, {B0}})),
		fixed("prefix",
			new Library[]{new Library("a", "a.fq", A0, A1), new Library("a_1", "a_1.fq", B0)},
			map(new String[]{"a", "a_1"}, new String[][]{{A0, A1}, {B0}})),
		new CaseSpec("fixed", "duplicate",
			new Library[]{new Library("left", "left/a.fq", A0, A1), new Library("right", "right/a.fq", B0)},
			map(new String[0], new String[0][]), "constructor", "java.lang.IllegalArgumentException",
			"Duplicate read basename 'a': %ROOT%/left/a.fq and %ROOT%/right/a.fq. "+
			"Each read library needs a unique basename because output paths use it.",
			new String[0], new String[0][0]),
		fixed("trailing_underscore",
			new Library[]{new Library("a_", "a_.fq", A0, A1), new Library("b", "b.fq", B0)},
			map(new String[]{"a_", "b"}, new String[][]{{A0, A1}, {B0}})),
		new CaseSpec("fixed", "duplicate_opfn",
			new Library[]{new Library("left", "left/a.fq", A0, A1), new Library("right", "right/a.fq", B0)},
			map(new String[]{"a"}, new String[][]{{A0, A1}}), null, null, null,
			new String[]{"opfn=1"}, new String[0][0]),
		new CaseSpec("fixed", "ambient_cwd",
			new Library[]{new Library("a", "a.fq", A0, A1), new Library("b", "b.fq", B0)},
			map(new String[]{"a", "b"}, new String[][]{{A0, A1}, {B0}}), null, null, null,
			new String[0], new String[][]{{"a_", "different_literal\n"}}),
		new CaseSpec("fixed", "legacy_namefile",
			new Library[]{new Library("a", "a.fq", A0, A1), new Library("b", "b.fq", B0)},
			map(new String[]{"a", "b"}, new String[][]{{A0, A1}, {B0}}), null, null, null,
			new String[0], new String[0][0], true)
	};

	/** Builds a baseline-phase case. */
	private static CaseSpec baseline(final String name, final Library[] libraries,
			final Map<String, List<String>> expected){
		return new CaseSpec("baseline", name, libraries, expected, null, null, null,
			new String[0], new String[0][0]);
	}

	/** Builds a fixed-phase case that should complete successfully. */
	private static CaseSpec fixed(final String name, final Library[] libraries,
			final Map<String, List<String>> expected){
		return new CaseSpec("fixed", name, libraries, expected, null, null, null,
			new String[0], new String[0][0]);
	}

	/** Builds expected memberships in insertion order. */
	private static Map<String, List<String>> map(final String[] names, final String[][] sequences){
		final LinkedHashMap<String, List<String>> map=new LinkedHashMap<String, List<String>>();
		for(int i=0; i<names.length; i++){map.put(names[i], Arrays.asList(sequences[i]));}
		return map;
	}

	private static final Method RENAME;
	private static final Method DEMUX;
	private static final Field READ_NAMES;
	static{
		try{
			RENAME=DecontaminateByNormalization.class.getDeclaredMethod("renameAndMux_MT",
				ArrayList.class, String.class);
			DEMUX=DecontaminateByNormalization.class.getDeclaredMethod("demux",
				String.class, ArrayList.class);
			READ_NAMES=DecontaminateByNormalization.class.getDeclaredField("readNames");
			RENAME.setAccessible(true);
			DEMUX.setAccessible(true);
			READ_NAMES.setAccessible(true);
		}catch(final Exception e){
			throw new RuntimeException(e);
		}
	}

	private static final class CaseSpec{
		CaseSpec(final String phase_, final String name_, final Library[] libraries_,
				final Map<String, List<String>> expected_, final String expectedFailurePhase_,
				final String expectedFailureType_, final String expectedFailureMessage_,
				final String[] extraArgs_, final String[][] cwdFiles_){
			this(phase_, name_, libraries_, expected_, expectedFailurePhase_, expectedFailureType_,
				expectedFailureMessage_, extraArgs_, cwdFiles_, false);
		}
		CaseSpec(final String phase_, final String name_, final Library[] libraries_,
				final Map<String, List<String>> expected_, final String expectedFailurePhase_,
				final String expectedFailureType_, final String expectedFailureMessage_,
				final String[] extraArgs_, final String[][] cwdFiles_,
				final boolean publicDemuxNamefile_){
			phase=phase_;
			name=name_;
			libraries=libraries_;
			expected=expected_;
			expectedFailurePhase=expectedFailurePhase_;
			expectedFailureType=expectedFailureType_;
			expectedFailureMessage=expectedFailureMessage_;
			extraArgs=extraArgs_;
			cwdFiles=new LinkedHashMap<String, String>();
			for(final String[] pair : cwdFiles_){
				cwdFiles.put(pair[0], pair[1]);
			}
			publicDemuxNamefile=publicDemuxNamefile_;
		}
		String expectedFailureMessage(final Path root){
			return expectedFailureMessage==null ? null :
				expectedFailureMessage.replace("%ROOT%", root.toString());
		}
		final String phase;
		final String name;
		final Library[] libraries;
		final Map<String, List<String>> expected;
		final String expectedFailurePhase;
		final String expectedFailureType;
		final String expectedFailureMessage;
		final String[] extraArgs;
		final LinkedHashMap<String, String> cwdFiles;
		final boolean publicDemuxNamefile;
	}

	private static final class Library{
		Library(final String name_, final String relativeRead_, final String... sequences_){
			name=name_;
			relativeRead=relativeRead_;
			refName=name+".fa";
			sequences=sequences_;
		}
		final String name;
		final String relativeRead;
		final String refName;
		final String[] sequences;
	}
}
