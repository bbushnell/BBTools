package test.jgi;

import java.io.IOException;
import java.lang.reflect.Field;
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

import fileIO.FileFormat;
import jgi.BBDuk;
import kmer.AbstractKmerTable;
import structures.IntList;

/** Checks one repaired filename case per fresh JVM through the public constructor.
 * Explicit readable prefix shadows and full target inputs match the saved baseline.
 * Both complete mate names and descriptors are required, so fixing the split alone
 * cannot pass while equals-based in2 clearing remains. No intermediate repair,
 * process call, reference loading, owner-field
 * writes or global resets occur. Normal probing, validation and allocation remain.
 * Raw arguments precede PreParser's permitted null rewrite; the after array is separate.
 * Constructor state is not evidence of paired processing or output correctness.
 * @author Sucrose
 */
public final class BBDukEqualsPathTest{

	/** Use exactly scenario=NAME outdir=NEW with the development launcher. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final int split=arg.indexOf('=');
			require(split>0 && split<arg.length()-1, "Use nonempty flag=value arguments: "+arg);
			require(options.put(arg.substring(0, split), arg.substring(split+1))==null, "Duplicate option: "+arg);
		}
		require(options.size()==2 && options.keySet().containsAll(Arrays.asList("scenario", "outdir")), "Use exactly scenario and outdir.");
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbdukequalspathtest.sh with assertions enabled for normal constructor validation.");
		final int ordinal=Arrays.asList(NAMES).indexOf(options.get("scenario"));
		require(ordinal>=0, "Choose one of eight fixed scenarios: "+Arrays.toString(NAMES));
		new BBDukEqualsPathTest(path(options.get("outdir")), ordinal).run();
	}

	private BBDukEqualsPathTest(final Path out_, final int ordinal_) throws IOException{
		out=out_; ordinal=ordinal_; name=NAMES[ordinal];
		Files.createDirectory(out);
		final StringBuilder manifest=new StringBuilder("file\tbases\n");
		for(int i=0; i<INPUTS[ordinal].length; i++){
			final String file=INPUTS[ordinal][i]; write(out.resolve(file), inputText(file, i));
			manifest.append(file).append('\t').append(BASES[i]).append('\n');
		}
		write(out.resolve("fixtures.tsv"), manifest.toString());
		write(out.resolve("settings.tsv"), lines(SETTINGS));
		write(out.resolve("counter-schema.tsv"), lines(ZERO_COUNTERS));
		write(out.resolve("reflection.tsv"), "operation\tmember\nORDINARY_PUBLIC_CALL\tnew BBDuk(public arguments) only\nPROCESS_CALLS\t0\nREAD_ONLY\tBBDuk names, FileFormat names, settings, counters, tables and dummy scaffold\nFIELD_WRITES\tNONE\nGLOBAL_RESETS\tNONE; fresh JVM per case\n");
		assertAbsent();
	}

	private void run() throws Exception{
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(SETTINGS));
		args.add("literal=ACGT");
		args.add("in="+out.resolve(RAW_IN1[ordinal]));
		if(RAW_IN2[ordinal]!=null){args.add("in2="+(RAW_IN2[ordinal].contains(".fq") ? out.resolve(RAW_IN2[ordinal]) : RAW_IN2[ordinal]));}
		if(ordinal==4){args.add("out="+out.resolve("output=result.fq"));}
		final String[] raw=args.toArray(new String[args.size()]), passed=raw.clone(), expectedAfter=raw.clone();
		if(ordinal==6){expectedAfter[expectedAfter.length-1]="in2=";}
		write(out.resolve("constructor-args.txt"), lines(raw));
		final BBDuk owner=new BBDuk(passed);
		write(out.resolve("constructor-args-after.txt"), lines(passed));
		check(Arrays.equals(passed, expectedAfter), "Only PreParser's entire null sentinel is rewritten in the caller array; case="+name);
		boolean inputsUnchanged=true;
		for(int i=0; i<INPUTS[ordinal].length; i++){inputsUnchanged&=Arrays.equals(Files.readAllBytes(out.resolve(INPUTS[ordinal][i])), inputText(INPUTS[ordinal][i], i).getBytes(StandardCharsets.US_ASCII));}
		check(inputsUnchanged, "Normal constructor probing must leave every explicit prefix and target FASTQ byte unchanged; case="+name);
		final String[] observed={(String)field("in1").get(owner), (String)field("in2").get(owner), (String)field("out1").get(owner)};
		final String[] target={resolved(TARGET_IN1[ordinal]), resolved(TARGET_IN2[ordinal]), ordinal==4 ? resolved("output=result.fq") : null};
		write(out.resolve("paths.tsv"), pathTable(observed)); write(out.resolve("target-paths.tsv"), pathTable(target));
		check(Arrays.equals(observed, target), "Fixed construction must preserve each complete target path, including both mates; case="+name+" observed="+Arrays.toString(observed));
		final String[] descriptors={descriptor(owner, "ffin1"), descriptor(owner, "ffin2"), descriptor(owner, "ffout1")};
		write(out.resolve("descriptors.tsv"), pathTable(descriptors));
		check(Arrays.equals(descriptors, target), "Normal input and output descriptors must preserve BOTH complete target mate paths and output; case="+name);
		final boolean discrepancy=!Arrays.equals(observed, target);
		check(!discrepancy, "All eight fixed cases must match their literal full-path targets; case="+name);
		check(number(owner, "k")==4 && number(owner, "mink")==4 && !flag(owner, "useShortKmers") && !flag(owner, "maskMiddle") && flag(owner, "ktrimN") && !flag(owner, "ktrimLeft") && !flag(owner, "ktrimRight") && !flag(owner, "ksplit") && number(owner, "THREADS")==1 && number(owner, "maxReads")==1 && !flag(owner, "rcomp") && field("out2").get(owner)==null, "Use disabled short kmers, ordinary N masking, one worker and finite input without processing; case="+name);
		final StringBuilder parsed=new StringBuilder("field\tvalue\n");
		for(final String member : PARSED_FIELDS){parsed.append(member).append('\t').append(String.valueOf(field(member).get(owner))).append('\n');}
		write(out.resolve("parsed-settings.tsv"), parsed.toString());
		final StringBuilder counters=new StringBuilder("field\tvalue\n"); boolean allZero=true;
		for(final String member : ZERO_COUNTERS){final long value=number(owner, member); allZero&=value==0; counters.append(member).append('\t').append(value).append('\n');}
		write(out.resolve("counters.tsv"), counters.toString());
		check(allZero, "Public construction must not load references or process reads; case="+name);
		final AbstractKmerTable[] tables=(AbstractKmerTable[])field("keySets").get(owner); int allocated=0, empty=0;
		for(final AbstractKmerTable table : tables){if(table!=null){allocated++; if(table.size()==0){empty++;}}}
		final List<?> names=(List<?>)field("scaffoldNames").get(owner); final IntList lengths=(IntList)field("scaffoldLengths").get(owner);
		check(tables.length==7 && allocated==7 && empty==7 && names.size()==1 && "".equals(names.get(0)) && lengths.size==1 && lengths.get(0)==0 && field("ref").get(owner)==null && Arrays.equals((String[])field("literal").get(owner), new String[]{"ACGT"}), "Ordinary construction allocates seven empty tables and dummy scaffold0; no reference-loading lifecycle is called.");
		write(out.resolve("allocation.tsv"), "tables\tallocated\tempty\tscaffold_names\tscaffold_lengths\n"+tables.length+"\t"+allocated+"\t"+empty+"\t"+names.size()+"\t"+lengths.size+"\n");
		check(absent(), "Both unexpanded templates and both output names must remain absent; no biological output is opened; case="+name);
		write(out.resolve("absent-paths.tsv"), "file\tbefore\tafter\npair#\tABSENT\tABSENT\npair#=batch.fq\tABSENT\tABSENT\noutput\tABSENT\tABSENT\noutput=result.fq\tABSENT\tABSENT\n");
		require(checks==9, "Nine independent constructor contracts are checked; actual="+checks);
		write(out.resolve("case.tsv"), "scenario\toutcome\tpath_discrepancy\tchecks\n"+name+"\tCONSTRUCTED\t"+discrepancy+"\t"+checks+"\n");
		final String summary="BBDUK_EQUALS_PATH_FIXED_PASS version=1 scenario="+name+" constructorAttempts=1 constructed=1 pathDiscrepancies="+(discrepancy ? 1 : 0)+" processCalls=0 checks="+checks+" failures=0\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
	}

	private String inputText(final String file, final int index){return "@equals-path-"+name+"-"+file+"\n"+BASES[index]+"\n+\nIIIIIIII\n";}
	private String resolved(final String file){return file==null ? null : out.resolve(file).toString();}
	private boolean absent(){for(final String file : ABSENT){if(Files.exists(out.resolve(file))){return false;}} return true;}
	private void assertAbsent(){require(absent(), "Fixture directory must not contain unexpanded template or output names: "+out);}
	private static String pathTable(final String[] paths){return "field\tvalue\nin1\t"+paths[0]+"\nin2\t"+paths[1]+"\nout1\t"+paths[2]+"\n";}
	private static String descriptor(final BBDuk owner, final String member) throws Exception{final FileFormat value=(FileFormat)field(member).get(owner); return value==null ? null : value.name();}
	private static Field field(final String name) throws NoSuchFieldException{final Field result=BBDuk.class.getDeclaredField(name); result.setAccessible(true); return result;}
	private static long number(final BBDuk owner, final String name) throws Exception{return ((Number)field(name).get(owner)).longValue();}
	private static boolean flag(final BBDuk owner, final String name) throws Exception{return field(name).getBoolean(owner);}
	private void check(final boolean value, final String why){checks++; require(value, why);}
	private static void require(final boolean value, final String why){if(!value){throw new IllegalArgumentException(why);}}
	private static String lines(final String[] values){return String.join("\n", values)+"\n";}
	private static void write(final Path file, final String value) throws IOException{Files.write(file, value.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);}
	private static Path path(final String value){require(value!=null && !value.isEmpty() && value.indexOf('=')<0 && value.indexOf('\t')<0 && value.indexOf('\n')<0, "Use a nonempty base directory without equals, tabs or newlines; intentional equals belong only in fixture filenames."); return Paths.get(value).toAbsolutePath().normalize();}

	private final Path out;
	private final int ordinal;
	private final String name;
	private int checks;
	private static final String[] NAMES={"plain_single", "in1_equals", "in2_equals", "template_equals", "out_equals", "plain_pair", "null_in2", "empty_in2"};
	private static final String[][] INPUTS={{"plain1.fq"}, {"first", "first=one=two.fq"}, {"plain1.fq", "second", "second=two.fq"}, {"pair1", "pair2", "pair1=batch.fq", "pair2=batch.fq"}, {"plain1.fq"}, {"plain1.fq", "plain2.fq"}, {"plain1.fq"}, {"plain1.fq"}};
	private static final String[] BASES={"GGACGTGG", "CCACGTCC", "TTACGTTT", "AAACGTAA"};
	private static final String[] RAW_IN1={"plain1.fq", "first=one=two.fq", "plain1.fq", "pair#=batch.fq", "plain1.fq", "plain1.fq", "plain1.fq", "plain1.fq"};
	private static final String[] RAW_IN2={null, null, "second=two.fq", null, null, "plain2.fq", "null", ""};
	private static final String[] TARGET_IN1={"plain1.fq", "first=one=two.fq", "plain1.fq", "pair1=batch.fq", "plain1.fq", "plain1.fq", "plain1.fq", "plain1.fq"};
	private static final String[] TARGET_IN2={null, null, "second=two.fq", "pair2=batch.fq", null, "plain2.fq", null, null};
	private static final String[] ABSENT={"pair#", "pair#=batch.fq", "output", "output=result.fq"};
	private static final String[] SETTINGS={"k=4", "kmask=N", "mink=-1", "maskmiddle=f", "rcomp=f", "hdist=0", "edist=0", "qhdist=0", "qskip=1", "rskip=1", "speed=0", "rename=f", "minlen=1", "qtrim=f", "maxns=-1", "maq=0", "reads=1", "t=1", "simd=f", "qin=33", "array=t", "prealloc=f", "ow=f", "app=f"};
	private static final String[] ZERO_COUNTERS={"readsIn", "basesIn", "readsOut", "basesOut", "storedKmers", "refReads", "refBases", "refKmers", "readsKFiltered", "basesKFiltered", "readsKTrimmed", "basesKTrimmed", "readsQTrimmed", "basesQTrimmed", "readsFTrimmed", "basesFTrimmed", "readsQFiltered", "basesQFiltered", "readsNFiltered", "basesNFiltered", "readsEFiltered", "basesEFiltered", "readsPolyTrimmed", "basesPolyTrimmed", "readsTrimmedByOverlap", "basesTrimmedByOverlap", "readsTrimmedBySwift", "basesTrimmedBySwift", "badGcReads", "badGcBases", "badHeaderReads", "badHeaderBases"};
	private static final String[] PARSED_FIELDS={"k", "mink", "useShortKmers", "maskMiddle", "midMaskLen", "middleMask", "minlen", "minminlen", "minlen2", "ktrimN", "ktrimLeft", "ktrimRight", "ksplit", "THREADS", "maxReads", "rcomp", "out2"};
}
