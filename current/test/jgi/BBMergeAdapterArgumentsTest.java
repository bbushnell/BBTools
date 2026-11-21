package test.jgi;

import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import jgi.BBMerge;
import parse.Parse;
import stream.Read;
import stream.Streamer;
import stream.Writer;

/** Eleven fresh-child constructor cases for JGI087/JGI092.
 * Literal expected mate lists distinguish six mate-option repairs and a filename
 * containing '=' from four shared/null-clearing controls. Both the truncated prefix
 * and complete filename are valid FASTAs, so baseline silently loads the wrong list
 * rather than failing to open an input. Uppercase adapters are exactly21 bases.
 * Two cases also invoke the real worker adapter predicate, without starting a worker:
 * a wrong configured mate must not be rescued by the unconfigured mate. Literal raw
 * and aligned read/quality inputs permit full orientation-restoration checks.
 * This tests construction and that predicate, not merger accuracy, all file loading,
 * arbitrary equals-containing inputs, or the separate earlier mode-preparse pass.
 * @author Sucrose
 */
public final class BBMergeAdapterArgumentsTest{

	/** phase=baseline|fixed, root, new outdir; fixed reference is the successful
	 * baseline FIXTURE directory. Internal mode=child adds case=name. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && Arrays.asList("mode", "phase", "root", "outdir", "reference", "case").contains(pair[0]), "Expected documented flag=value: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate fixture option: "+pair[0]);
		}
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Use dev/bbmergeadapterargumentstest.sh with assertions enabled.");
		final String mode=options.containsKey("mode") ? options.get("mode") : "suite", phase=options.get("phase");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		final Path root=path(options.get("root")), out=path(options.get("outdir"));
		require(root.toString().indexOf('=')<0 && out.toString().indexOf('=')<0, "Keep fixture roots free of '='; only the adapter filename is the constructor suffix witness.");
		if(mode.equals("child")){
			require(!options.containsKey("reference") && Files.isDirectory(out), "Child uses an existing case directory and no reference.");
			child(root, out, new Spec(options.get("case")));
			return;
		}
		require(mode.equals("suite") && !options.containsKey("case"), "Use suite or internal child mode.");
		final Path reference=options.containsKey("reference") ? path(options.get("reference")) : null;
		require(phase.equals("baseline") ? reference==null : reference!=null, "Only fixed needs reference=<successful baseline fixture outdir>.");
		if(reference!=null){require(read(reference.resolve("summary.txt")).startsWith("BBMERGE_ADAPTER_ARGUMENTS_PASS phase=baseline cases=11 children=11 "), "Reference is not a successful eleven-child baseline fixture directory.");}
		final BBMergeAdapterArgumentsTest test=new BBMergeAdapterArgumentsTest(root, out, reference, phase);
		for(final String name : NAMES){test.run(new Spec(name));}
		test.check(test.children==11 && test.changed==7 && test.controls==4 && test.predicates==2, "Keep the eleven-case contract: seven changed configurations, four controls and two predicates.");
		test.check(test.sourceMatches==(phase.equals("baseline") ? 11 : 4), "Baseline configuration match count differs from literal oracles.");
		test.check(test.correctedMatches==(phase.equals("baseline") ? 4 : 11), "Corrected configuration match count differs from literal oracles.");
		final String summary="BBMERGE_ADAPTER_ARGUMENTS_"+(test.failures==0 ? "PASS" : "FAIL")+" phase="+phase+" cases=11 children="+test.children+
			" changed="+test.changed+" controls="+test.controls+" predicates="+test.predicates+" sourceMatches="+test.sourceMatches+
			" correctedMatches="+test.correctedMatches+" checks="+test.checks+" failures="+test.failures+"\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
		if(test.failures>0){throw new AssertionError(summary+"Inspect saved constructor arguments, fields and predicate evidence; never fit expected lists to actual loading.");}
	}

	private BBMergeAdapterArgumentsTest(final Path root_, final Path out_, final Path reference_, final String phase_) throws IOException{
		root=root_; out=out_; reference=reference_; phase=phase_;
		Files.createDirectory(out);
		write(out.resolve("cases.tsv"), "case\tchildExit\tbaselineConfigMatches\tcorrectedConfigMatches\tpredicate\tphaseMatches\n");
	}

	/** Preserves constructor inputs and separate wrong/correct field oracles. */
	private void run(final Spec spec) throws Exception{
		if(spec.baseline().equals(spec.corrected())){controls++;}else{changed++;}
		if(spec.probe){predicates++;}
		final Path work=Files.createDirectory(out.resolve(spec.name));
		write(work.resolve("A.fq"), "@adapter-case/1\n"+R1_BASES+"\n+\n56789:;<=>?@AB\n");
		write(work.resolve("B.fq"), "@adapter-case/2\n"+RAW_B+"\n+\nDEFGHIJKLMNOPQ\n");
		final ArrayList<String> inputs=new ArrayList<String>(Arrays.asList("A.fq", "B.fq"));
		if(spec.name.equals("file-equals")){
			write(work.resolve("prefix.fa"), ">prefix-A\n"+A+"\n");
			write(work.resolve("prefix.fa=full.fa"), ">full-B\n"+B+"\n");
			inputs.add("prefix.fa"); inputs.add("prefix.fa=full.fa");
		}
		write(work.resolve("expected-baseline.tsv"), spec.baseline());
		write(work.resolve("expected-fixed.tsv"), spec.corrected());
		final Map<String, byte[]> before=new LinkedHashMap<String, byte[]>();
		for(final String file : inputs){before.put(file, Files.readAllBytes(work.resolve(file)));}
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("bash", root.resolve("dev/bbmergeadapterargumentstest.sh").toString(),
			"mode=child", "phase="+phase, "root="+root, "outdir="+work, "case="+spec.name));
		write(work.resolve("child.command.tsv"), argumentsText(command));
		final ProcessBuilder builder=new ProcessBuilder(command).directory(root.toFile());
		builder.redirectOutput(work.resolve("child.stdout").toFile()); builder.redirectError(work.resolve("child.stderr").toFile());
		builder.environment().remove("_JAVA_OPTIONS"); builder.environment().remove("JDK_JAVA_OPTIONS");
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		final Process process=builder.start(); children++;
		final int code=process.waitFor(); write(work.resolve("child.exit"), code+"\n");
		check(code==0, spec.name+": constructor/loading/predicate child failed.");
		for(final String file : inputs){check(Arrays.equals(before.get(file), Files.readAllBytes(work.resolve(file))), spec.name+": constructor changed input file "+file);}
		if(code!=0){append(out.resolve("cases.tsv"), spec.name+"\t"+code+"\tfalse\tfalse\tNA\tfalse\n"); return;}
		final byte[] actual=Files.readAllBytes(work.resolve("actual.tsv"));
		final boolean source=Arrays.equals(actual, spec.baseline().getBytes(StandardCharsets.US_ASCII));
		final boolean corrected=Arrays.equals(actual, spec.corrected().getBytes(StandardCharsets.US_ASCII));
		if(source){sourceMatches++;} if(corrected){correctedMatches++;}
		boolean matches=phase.equals("baseline") ? source : corrected;
		check(matches, spec.name+": complete adapter-list contents/null state/verification flag differs from the literal phase oracle.");
		String predicate="NA";
		if(spec.probe){
			final String predicateText=read(work.resolve("predicate.tsv"));
			final String expected="wrongConfiguredStatus\t0\nmatchingUnconfiguredStatus\t2\nreturned\t"+phase.equals("baseline")+"\norientationRestored\ttrue\n";
			final boolean predicateMatches=predicateText.equals(expected);
			check(predicateMatches, spec.name+": actual worker predicate must change true to false with exact input restoration.");
			matches&=predicateMatches;
			check(read(work.resolve("orientation-before.tsv")).equals(expectedReadState()), spec.name+": predicate was not given the intended literal aligned mates and qualities.");
			check(read(work.resolve("orientation-before.tsv")).equals(read(work.resolve("orientation-after.tsv"))), spec.name+": worker did not restore bases/qualities.");
			predicate="MISSING";
			for(final String row : predicateText.split("\n")){if(row.startsWith("returned\t")){predicate=row.substring("returned\t".length()); break;}}
		}
		if(reference!=null){
			final Path old=reference.resolve(spec.name);
			for(final String file : inputs){same(work.resolve(file), old.resolve(file), spec.name+": input differs from baseline: "+file);}
			for(final String file : new String[]{"expected-baseline.tsv", "expected-fixed.tsv"}){same(work.resolve(file), old.resolve(file), spec.name+": configuration oracle changed: "+file);}
			if(!spec.probe && spec.baseline().equals(spec.corrected())){same(work.resolve("actual.tsv"), old.resolve("actual.tsv"), spec.name+": unchanged control configuration changed.");}
			if(spec.probe){same(work.resolve("orientation-before.tsv"), old.resolve("orientation-before.tsv"), spec.name+": predicate inputs changed between phases.");}
		}
		append(out.resolve("cases.tsv"), spec.name+"\t"+code+"\t"+source+"\t"+corrected+"\t"+predicate+"\t"+matches+"\n");
	}

	/** Reads actual constructor fields; only the tested constructor sets adapter state. */
	private static void child(final Path root, final Path work, final Spec spec) throws Exception{
		require(A.length()==21 && B.length()==21 && G.length()==21 && R1_BASES.length()==14 && RAW_B.length()==14 && ALIGNED_B.length()==14, "Keep adapter and predicate geometry literal and finite.");
		for(final String literal : new String[]{A, B, G}){require(!Files.exists(root.resolve(literal)), "Literal adapter must not accidentally name a file in the child working directory.");}
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in1="+work.resolve("A.fq"), "in2="+work.resolve("B.fq"), "t=1", "simd=f", "nn=f", "qtrim=f", "minlength=1", "ow=t"));
		args.addAll(Arrays.asList(spec.options(work)));
		write(work.resolve("constructor.command.tsv"), argumentsText(args));
		final BBMerge merge=new BBMerge(args.toArray(new String[0]));
		require(((Number)value(merge, "maxAdapterLength")).intValue()==21, "Adapter loading contract must retain the known maximum length21.");
		write(work.resolve("actual.tsv"), actualList(merge, "adapterList1", "mate1")+actualList(merge, "adapterList2", "mate2")+
			"verifyAdapters\t"+value(merge, "verifyAdapters")+"\n");
		if(spec.probe){predicate(work, merge, spec.name.equals("r1-only"));}
	}

	/** A bad configured tail cannot be rescued by a good unconfigured tail after repair. */
	private static void predicate(final Path work, final BBMerge merge, final boolean firstOnly) throws Exception{
		final Read a=new Read(R1_BASES.getBytes(StandardCharsets.US_ASCII), A_QUAL.clone(), "adapter-case/1", 0, false);
		final Read b=new Read(ALIGNED_B.getBytes(StandardCharsets.US_ASCII), ALIGNED_B_QUAL.clone(), "adapter-case/2", 0, false);
		a.mate=b; b.mate=a; b.setPairnum(1);
		final Read rawB=new Read(RAW_B.getBytes(StandardCharsets.US_ASCII), RAW_B_QUAL.clone(), "adapter-case/2", 0, false);
		final byte[] adapter=(firstOnly ? A : B).getBytes(StandardCharsets.US_ASCII);
		final int wrong=(Integer)invoke(BBMerge.class, null, "adapterIsValid", new Class<?>[]{Read.class, byte[].class, int.class, int.class, float.class}, firstOnly ? a : rawB, adapter, 8, 3, .6f);
		final int matching=(Integer)invoke(BBMerge.class, null, "adapterIsValid", new Class<?>[]{Read.class, byte[].class, int.class, int.class, float.class}, firstOnly ? rawB : a, adapter, 8, 3, .6f);
		require(wrong==0 && matching==2, "Literal six-base tails must prove bad configured and good unconfigured adapter matches before testing the worker.");
		final Class<?> workerClass=Class.forName("jgi.BBMerge$MateThread");
		final Constructor<?> constructor=workerClass.getDeclaredConstructor(BBMerge.class, Streamer.class, Writer.class, Writer.class, fileIO.ByteStreamWriter.class, boolean.class, boolean.class);
		constructor.setAccessible(true);
		final Object worker=constructor.newInstance(merge, null, null, null, null, true, false);
		require(((Thread)worker).getState()==Thread.State.NEW, "A predicate probe must not start a worker or any production stream.");
		final String before=readState(a, b); write(work.resolve("orientation-before.tsv"), before);
		final boolean returned=(Boolean)invoke(workerClass, worker, "verifyAdapters", new Class<?>[]{Read.class, Read.class, int.class}, a, b, 8);
		final String after=readState(a, b); write(work.resolve("orientation-after.tsv"), after);
		final boolean restored=before.equals(after);
		write(work.resolve("predicate.tsv"), "wrongConfiguredStatus\t"+wrong+"\nmatchingUnconfiguredStatus\t"+matching+"\nreturned\t"+returned+"\norientationRestored\t"+restored+"\n");
		require(restored, "Worker predicate must restore aligned mate2 and preserve every input base and quality.");
	}

	/** Emits every loaded sequence; expected configurations below use separate literals. */
	@SuppressWarnings("unchecked")
	private static String actualList(final BBMerge merge, final String field, final String label) throws Exception{
		final List<byte[]> list=(List<byte[]>)value(merge, field);
		if(list==null){return label+"\tNULL\n";}
		final StringBuilder text=new StringBuilder(label).append('\t').append(list.size()).append('\n');
		for(int i=0; i<list.size(); i++){text.append(label).append('.').append(i).append('\t').append(new String(list.get(i), StandardCharsets.US_ASCII)).append('\n');}
		return text.toString();
	}
	private static String expectedList(final String label, final String literal){return literal==null ? label+"\tNULL\n" : label+"\t1\n"+label+".0\t"+literal+"\n";}
	private static String expectedConfig(final String first, final String second){return expectedList("mate1", first)+expectedList("mate2", second)+"verifyAdapters\t"+(first!=null || second!=null)+"\n";}

	/** Fixed cases retain argument order, singular/plural aliases and null-clearing controls. */
	private static final class Spec{
		Spec(final String name_){
			name=name_;
			require(name!=null && Arrays.asList(NAMES).contains(name), "Unknown adapter constructor case: "+name);
			probe=name.equals("r1-only") || name.equals("r2-only");
			switch(name){
				case "r1-only": old1=A; old2=A; fixed1=A; fixed2=null; flags=new String[]{"adapter1="+A}; break;
				case "r2-only": old1=B; old2=B; fixed1=null; fixed2=B; flags=new String[]{"adapters2="+B}; break;
				case "shared-r1": old1=A; old2=A; fixed1=A; fixed2=G; flags=new String[]{"adapter="+G, "adapters1="+A}; break;
				case "shared-r2": old1=B; old2=B; fixed1=G; fixed2=B; flags=new String[]{"adapters="+G, "adapter2="+B}; break;
				case "r1-r2": old1=B; old2=B; fixed1=A; fixed2=B; flags=new String[]{"adapters1="+A, "adapter2="+B}; break;
				case "r2-r1": old1=A; old2=A; fixed1=A; fixed2=B; flags=new String[]{"adapters2="+B, "adapter1="+A}; break;
				case "shared": old1=G; old2=G; fixed1=G; fixed2=G; flags=new String[]{"adapters="+G}; break;
				case "generic-last": old1=G; old2=G; fixed1=G; fixed2=G; flags=new String[]{"adapter1="+A, "adapters2="+B, "adapter="+G}; break;
				case "file-equals": old1=A; old2=A; fixed1=B; fixed2=B; flags=null; break;
				case "clear-bare": old1=null; old2=null; fixed1=null; fixed2=null; flags=new String[]{"adapter="+G, "adapters"}; break;
				case "clear-empty": old1=null; old2=null; fixed1=null; fixed2=null; flags=new String[]{"adapters="+G, "adapters="}; break;
				default:throw new AssertionError("Unreachable case: "+name);
			}
		}
		String[] options(final Path work){return flags==null ? new String[]{"adapters="+work.resolve("prefix.fa=full.fa")} : flags.clone();}
		String baseline(){return expectedConfig(old1, old2);}
		String corrected(){return expectedConfig(fixed1, fixed2);}
		final String name, old1, old2, fixed1, fixed2;
		final String[] flags;
		final boolean probe;
	}

	private static String argumentsText(final List<String> args){final StringBuilder text=new StringBuilder(); for(int i=0; i<args.size(); i++){text.append(i).append('\t').append(args.get(i)).append('\n');} return text.toString();}
	private static String readState(final Read a, final Read b){return "mate\tbases\tqualities\nmate1\t"+new String(a.bases, StandardCharsets.US_ASCII)+"\t"+Arrays.toString(a.quality)+"\nmate2-aligned\t"+new String(b.bases, StandardCharsets.US_ASCII)+"\t"+Arrays.toString(b.quality)+"\n";}
	private static String expectedReadState(){return "mate\tbases\tqualities\nmate1\t"+R1_BASES+"\t"+Arrays.toString(A_QUAL)+"\nmate2-aligned\t"+ALIGNED_B+"\t"+Arrays.toString(ALIGNED_B_QUAL)+"\n";}
	private static Object value(final Object owner, final String name) throws Exception{final Field f=owner.getClass().getDeclaredField(name); f.setAccessible(true); return f.get(owner);}
	private static Object invoke(final Class<?> type, final Object owner, final String name, final Class<?>[] types, final Object... args) throws Exception{final Method m=type.getDeclaredMethod(name, types); m.setAccessible(true); return m.invoke(owner, args);}
	private void same(final Path a, final Path b, final String message) throws IOException{check(Arrays.equals(Files.readAllBytes(a), Files.readAllBytes(b)), message);}
	private void check(final boolean condition, final String message){checks++; if(!condition){failures++; System.err.println("CHECK_FAIL\t"+message);}}
	private static void require(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
	private static Path path(final String value){require(value!=null && !value.isEmpty(), "Missing required path."); return Paths.get(value).toAbsolutePath().normalize();}
	private static void write(final Path path, final String text) throws IOException{Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path path, final String text) throws IOException{Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);}
	private static String read(final Path path) throws IOException{return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);}

	private final Path root, out, reference;
	private final String phase;
	private int checks, failures, children, changed, controls, predicates, sourceMatches, correctedMatches;
	private static final String A="AAAAAAAAAAAAAAAAAAAAA", B="CCCCCCCCCCCCCCCCCCCCC", G="GGGGGGGGGGGGGGGGGGGGG";
	private static final String R1_BASES="ACGTTCGACCCCCC", RAW_B="TCGAACGTAAAAAA", ALIGNED_B="TTTTTTACGTTCGA";
	private static final byte[] A_QUAL={20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33};
	private static final byte[] RAW_B_QUAL={35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48};
	private static final byte[] ALIGNED_B_QUAL={48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35};
	private static final String[] NAMES={"r1-only", "r2-only", "shared-r1", "shared-r2", "r1-r2", "r2-r1", "shared", "generic-last", "file-equals", "clear-bare", "clear-empty"};
}
