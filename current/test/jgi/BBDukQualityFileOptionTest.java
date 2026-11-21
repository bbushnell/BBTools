package test.jgi;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
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

/** Checks one repaired quality-file option case per fresh public constructor.
 * All eight cases must construct successfully; every caught exception fails the
 * observer after its full trace is saved. Construction exposes all four
 * quality filenames and both sequence descriptors. Valid paired FASTA/QUAL fixtures,
 * normal file validation/probing and table allocation remain. No processing, stream
 * factory, reference loading, owner-field writes or artificial global resets occur.
 * Constructor observations do not
 * establish actual quality consumption, output bytes or I/O lifecycle correctness.
 * @author Sucrose
 */
public final class BBDukQualityFileOptionTest{

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
		require(assertions, "Use dev/bbdukqualityfileoptiontest.sh with assertions enabled for normal constructor validation.");
		final int ordinal=Arrays.asList(NAMES).indexOf(options.get("scenario"));
		require(ordinal>=0, "Choose one of eight fixed scenarios: "+Arrays.toString(NAMES));
		new BBDukQualityFileOptionTest(path(options.get("outdir")), ordinal).run();
	}

	private BBDukQualityFileOptionTest(final Path out_, final int ordinal_) throws IOException{
		out=out_; ordinal=ordinal_; name=NAMES[ordinal]; Files.createDirectory(out);
		for(int i=0; i<INPUTS.length; i++){write(out.resolve(INPUTS[i]), inputText(i));}
		write(out.resolve("fixtures.txt"), lines(INPUTS)); write(out.resolve("settings.tsv"), lines(SETTINGS));
		write(out.resolve("counter-schema.tsv"), lines(ZERO_COUNTERS));
		write(out.resolve("reflection.tsv"), "operation\tmember\nORDINARY_PUBLIC_CALL\tnew BBDuk(public arguments) only\nPROCESS_CALLS\t0\nFACTORY_CALLS\t0\nREAD_ONLY\tQuality filenames, sequence filenames and descriptors, settings, counters, tables and dummy scaffold\nFIELD_WRITES\tNONE\nGLOBAL_RESETS\tNONE; fresh JVM per case\n");
		require(absent(), "Fixture directory must not contain sequence or quality output paths: "+out);
	}

	private void run() throws Exception{
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(SETTINGS));
		args.add("literal=ACGT"); args.add("in="+out.resolve("read1.fa")); args.add("in2="+out.resolve("read2.fa"));
		args.add("out="+out.resolve("out1.fa")); args.add("out2="+out.resolve("out2.fa"));
		for(final String option : OPTIONS[ordinal]){final int split=option.indexOf('='); args.add(option.substring(0, split+1)+out.resolve(option.substring(split+1)));}
		final String[] raw=args.toArray(new String[args.size()]), passed=raw.clone();
		write(out.resolve("constructor-args.txt"), lines(raw));
		BBDuk owner=null; RuntimeException error=null;
		try{owner=new BBDuk(passed);}catch(final RuntimeException caught){error=caught;}
		write(out.resolve("constructor-args-after.txt"), lines(passed));
		final boolean rejected=error!=null;
		if(rejected){
			final StringWriter trace=new StringWriter(); final PrintWriter writer=new PrintWriter(trace);
			error.printStackTrace(writer); writer.flush(); write(out.resolve("parser-exception.txt"), trace.toString());
		}else{write(out.resolve("parser-exception.txt"), "NONE\n");}
		boolean unchanged=true;
		for(int i=0; i<INPUTS.length; i++){unchanged&=Arrays.equals(Files.readAllBytes(out.resolve(INPUTS[i])), inputText(i).getBytes(StandardCharsets.US_ASCII));}
		check(unchanged, "Normal construction must leave both FASTA and both matching QUAL inputs byte-exact; case="+name);
		check(absent(), "Construction must create no sequence output, quality output or literal template; case="+name);
		write(out.resolve("absent-paths.txt"), lines(ABSENT));
		check(Arrays.equals(raw, passed), "These nonnull quality arguments must remain complete in the caller array; case="+name);
		check(!rejected, "All eight fixed cases must construct, including explicit second quality output; case="+name+" actual="+error);
		final StringBuilder quality=new StringBuilder("field\tvalue\n"); boolean qualityOK=true;
		for(int i=0; i<QUALITY_FIELDS.length; i++){
			final String observed=(String)field(QUALITY_FIELDS[i]).get(owner), expected=QUALITY_PATHS[ordinal][i]==null ? null : out.resolve(QUALITY_PATHS[ordinal][i]).toString();
			qualityOK&=observed==null ? expected==null : observed.equals(expected); quality.append(QUALITY_FIELDS[i]).append('\t').append(observed).append('\n');
		}
		write(out.resolve("quality-paths.tsv"), quality.toString());
		check(qualityOK, "Every input and output quality field must keep its own option; qfin2 must not become qfout2, and # output expansion must retain both mates; case="+name);
		final StringBuilder paths=new StringBuilder("field\tvalue\n"); boolean pathsOK=true;
		for(int i=0; i<SEQUENCE_FIELDS.length; i++){
			final String observed=(String)field(SEQUENCE_FIELDS[i]).get(owner); final FileFormat descriptor=(FileFormat)field(FORMAT_FIELDS[i]).get(owner);
			pathsOK&=out.resolve(SEQUENCE_PATHS[i]).toString().equals(observed) && descriptor!=null && observed.equals(descriptor.name()) && descriptor.fasta();
			paths.append(SEQUENCE_FIELDS[i]).append('\t').append(observed).append('\n').append(FORMAT_FIELDS[i]).append('\t').append(descriptor==null ? null : descriptor.name()).append('\n');
		}
		write(out.resolve("sequence-paths.tsv"), paths.toString());
		check(pathsOK, "Both normal FASTA inputs and both absent sequence output descriptors must remain exact; case="+name);
		check(number(owner, "k")==4 && number(owner, "mink")==4 && !flag(owner, "useShortKmers") && !flag(owner, "maskMiddle") && flag(owner, "ktrimN") && !flag(owner, "ktrimLeft") && !flag(owner, "ktrimRight") && !flag(owner, "ksplit") && number(owner, "THREADS")==1 && number(owner, "maxReads")==1 && !flag(owner, "rcomp"), "Use ordinary no-short N masking and finite one-worker constructor configuration; case="+name);
		final StringBuilder parsed=new StringBuilder("field\tvalue\n");for(final String member : PARSED_FIELDS){parsed.append(member).append('\t').append(field(member).get(owner)).append('\n');} write(out.resolve("parsed-settings.tsv"), parsed.toString());
		final StringBuilder counters=new StringBuilder("field\tvalue\n"); boolean allZero=true;
		for(final String member : ZERO_COUNTERS){final long value=number(owner, member); allZero&=value==0; counters.append(member).append('\t').append(value).append('\n');} write(out.resolve("counters.tsv"), counters.toString());
		check(allZero, "Construction must not load references or process paired sequence/quality records; case="+name);
		final AbstractKmerTable[] tables=(AbstractKmerTable[])field("keySets").get(owner); int allocated=0, empty=0;
		for(final AbstractKmerTable table : tables){if(table!=null){allocated++; if(table.size()==0){empty++;}}}
		final List<?> names=(List<?>)field("scaffoldNames").get(owner); final IntList lengths=(IntList)field("scaffoldLengths").get(owner);
		check(tables.length==7 && allocated==7 && empty==7 && names.size()==1 && "".equals(names.get(0)) && lengths.size==1 && lengths.get(0)==0 && field("ref").get(owner)==null && Arrays.equals((String[])field("literal").get(owner), new String[]{"ACGT"}), "Ordinary constructor allocates seven empty tables and dummy scaffold0; no loading lifecycle is invoked.");
		write(out.resolve("allocation.tsv"), "tables\tallocated\tempty\tscaffold_names\tscaffold_lengths\n"+tables.length+"\t"+allocated+"\t"+empty+"\t"+names.size()+"\t"+lengths.size+"\n");
		require(checks==9, "Four common checks and five successful-state checks are required for every fixed case; actual="+checks);
		write(out.resolve("case.tsv"), "scenario\toutcome\tchecks\n"+name+"\tCONSTRUCTED\t"+checks+"\n");
		final String summary="BBDUK_QUALITY_FILE_FIXED_PASS version=1 scenario="+name+" constructorAttempts=1 constructed=1 unknownQfout2=0 processCalls=0 factoryCalls=0 checks="+checks+" failures=0\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
	}

	private String inputText(final int index){final int mate=(index&1)+1; return ">quality-option-"+name+"/"+mate+"\n"+(index==0 ? "GGACGTGG" : index==1 ? "CCACGTCC" : index==2 ? "40 39 38 37 36 35 34 33" : "32 31 30 29 28 27 26 25")+"\n";}
	private boolean absent(){for(final String file : ABSENT){if(Files.exists(out.resolve(file))){return false;}} return true;}
	private static Field field(final String name) throws NoSuchFieldException{final Field result=BBDuk.class.getDeclaredField(name); result.setAccessible(true); return result;}
	private static long number(final BBDuk owner, final String name) throws Exception{return ((Number)field(name).get(owner)).longValue();}
	private static boolean flag(final BBDuk owner, final String name) throws Exception{return field(name).getBoolean(owner);}
	private void check(final boolean value, final String why){checks++; require(value, why);}
	private static void require(final boolean value, final String why){if(!value){throw new IllegalArgumentException(why);}}
	private static String lines(final String[] values){return String.join("\n", values)+"\n";}
	private static void write(final Path file, final String value) throws IOException{Files.write(file, value.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);}
	private static Path path(final String value){require(value!=null && !value.isEmpty() && value.indexOf('=')<0 && value.indexOf('\t')<0 && value.indexOf('\n')<0, "Use a nonempty directory without equals, tabs or newlines for reproducible option fixtures."); return Paths.get(value).toAbsolutePath().normalize();}

	private final Path out;
	private final int ordinal;
	private final String name;
	private int checks;
	private static final String[] NAMES={"plain_pair", "qfin2_only", "qfin_pair", "qfout_alias", "qfout1_explicit", "qfout2_explicit", "quality_pair", "qfout_template"};
	private static final String[] INPUTS={"read1.fa", "read2.fa", "input1.qual", "input2.qual"};
	private static final String[] ABSENT={"out1.fa", "out2.fa", "output1.qual", "output2.qual", "quality#.qual", "quality1.qual", "quality2.qual"};
	private static final String[][] OPTIONS={{}, {"qfin2=input2.qual"}, {"qfin=input1.qual", "qfin2=input2.qual"}, {"qfout=output1.qual"}, {"qfout1=output1.qual"}, {"qfout2=output2.qual"}, {"qfin2=input2.qual", "qfout1=output1.qual", "qfout2=output2.qual"}, {"qfout=quality#.qual"}};
	private static final String[] QUALITY_FIELDS={"qfin1", "qfin2", "qfout1", "qfout2"};
	private static final String[][] QUALITY_PATHS={{null, null, null, null}, {null, "input2.qual", null, null}, {"input1.qual", "input2.qual", null, null}, {null, null, "output1.qual", null}, {null, null, "output1.qual", null}, {null, null, null, "output2.qual"}, {null, "input2.qual", "output1.qual", "output2.qual"}, {null, null, "quality1.qual", "quality2.qual"}};
	private static final String[] SEQUENCE_FIELDS={"in1", "in2", "out1", "out2"}, FORMAT_FIELDS={"ffin1", "ffin2", "ffout1", "ffout2"}, SEQUENCE_PATHS={"read1.fa", "read2.fa", "out1.fa", "out2.fa"};
	private static final String[] SETTINGS={"k=4", "kmask=N", "mink=-1", "maskmiddle=f", "rcomp=f", "hdist=0", "edist=0", "qhdist=0", "qskip=1", "rskip=1", "speed=0", "rename=f", "minlen=1", "qtrim=f", "maxns=-1", "maq=0", "reads=1", "t=1", "simd=f", "qin=33", "array=t", "prealloc=f", "ow=f", "app=f"};
	private static final String[] ZERO_COUNTERS={"readsIn", "basesIn", "readsOut", "basesOut", "storedKmers", "refReads", "refBases", "refKmers", "readsKFiltered", "basesKFiltered", "readsKTrimmed", "basesKTrimmed", "readsQTrimmed", "basesQTrimmed", "readsFTrimmed", "basesFTrimmed", "readsQFiltered", "basesQFiltered", "readsNFiltered", "basesNFiltered", "readsEFiltered", "basesEFiltered", "readsPolyTrimmed", "basesPolyTrimmed", "readsTrimmedByOverlap", "basesTrimmedByOverlap", "readsTrimmedBySwift", "basesTrimmedBySwift", "badGcReads", "badGcBases", "badHeaderReads", "badHeaderBases"};
	private static final String[] PARSED_FIELDS={"k", "mink", "useShortKmers", "maskMiddle", "midMaskLen", "middleMask", "minlen", "minminlen", "minlen2", "ktrimN", "ktrimLeft", "ktrimRight", "ksplit", "THREADS", "maxReads", "rcomp"};
}
