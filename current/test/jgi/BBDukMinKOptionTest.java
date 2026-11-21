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

import jgi.BBDuk;
import kmer.AbstractKmerTable;
import structures.IntList;

/** One fixed option case per fresh JVM, using the ordinary public constructor.
 * Nonpositive mink values disable automatic short-kmer selection; explicit usk
 * remains independent. Only the exact upper32 parser assertion is captured;
 * every other failure fails this observer. Normal input checks, format probing,
 * constructor static side effects and table allocation remain. No process(),
 * reference loading, read processing, state seeding or artificial resets occur.
 * Successful owners expose effective settings, empty tables and zero counters;
 * rejected construction has no usable owner state and is recorded as unavailable.
 * This checks parsing/configuration, not biological output or an I/O lifecycle.
 * @author Sucrose
 */
public final class BBDukMinKOptionTest{

	/** Use exactly scenario=NAME outdir=NEW; fixed only, in a fresh JVM. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final int split=arg.indexOf('=');
			require(split>0 && split<arg.length()-1, "Use nonempty flag=value arguments: "+arg);
			require(options.put(arg.substring(0, split), arg.substring(split+1))==null, "Duplicate option: "+arg);
		}
		require(options.size()==2 && options.keySet().containsAll(Arrays.asList("scenario", "outdir")), "Use exactly scenario and outdir.");
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbdukminkoptiontest.sh with assertions enabled; parser rejection is the observed contract.");
		final int ordinal=Arrays.asList(NAMES).indexOf(options.get("scenario"));
		require(ordinal>=0, "Choose one of the eight fixed scenarios: "+Arrays.toString(NAMES));
		new BBDukMinKOptionTest(path(options.get("outdir")), ordinal).run();
	}

	private BBDukMinKOptionTest(final Path out_, final int ordinal_) throws IOException{
		out=out_; ordinal=ordinal_; name=NAMES[ordinal]; input=out.resolve("input.fq");
		Files.createDirectory(out); write(input, inputText());
		write(out.resolve("settings.tsv"), String.join("\n", SETTINGS)+"\n");
		write(out.resolve("reflection.tsv"), "operation\tmember\nORDINARY_PUBLIC_CALL\tnew BBDuk(public arguments) only\nPROCESS_CALLS\t0\nREAD_ONLY\tBBDuk effective settings, counters, tables and dummy scaffold\nFIELD_WRITES\tNONE\nGLOBAL_RESETS\tNONE; fresh JVM per case\n");
		write(out.resolve("counter-schema.tsv"), String.join("\n", ZERO_COUNTERS)+"\n");
	}

	private void run() throws Exception{
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(SETTINGS));
		if(!OPTIONS[ordinal].isEmpty()){args.add(OPTIONS[ordinal]);}
		if(EXPLICIT_SHORT[ordinal]){args.add("usk=t");}
		args.add("literal=ACGT"); args.add("in="+input);
		write(out.resolve("constructor-args.txt"), String.join("\n", args)+"\n");
		BBDuk owner=null; AssertionError error=null;
		try{owner=new BBDuk(args.toArray(new String[args.size()]));}
		catch(final AssertionError caught){error=caught;}
		final boolean rejected=error!=null;
		if(rejected){
			final StringWriter trace=new StringWriter(); final PrintWriter writer=new PrintWriter(trace);
			error.printStackTrace(writer); writer.flush(); write(out.resolve("parser-assertion.txt"), trace.toString());
		}else{write(out.resolve("parser-assertion.txt"), "NONE\n");}
		check(Arrays.equals(Files.readAllBytes(input), inputText().getBytes(StandardCharsets.US_ASCII)), "Public construction must leave the finite input unchanged; case="+name);
		check(rejected==REJECTED[ordinal], "Fixed has seven successful constructions and the retained upper32 parser assertion; case="+name+" rejected="+rejected);
		String effective="NA", shortKmers="NA", middle="NA";
		if(rejected){
			check(error.getClass()==AssertionError.class && PARSER_MESSAGE.equals(error.getMessage()), "Capture only the known mink/kmin parser assertion, never another constructor failure; case="+name+" actual="+error);
			final StackTraceElement[] frames=error.getStackTrace();
			check(frames.length>0 && frames[0].getClassName().equals("jgi.BBDuk") && frames[0].getMethodName().equals("<init>") && "BBDuk.java".equals(frames[0].getFileName()) && frames[0].getLineNumber()==255, "The repaired range assertion remains at BBDuk.java255; case="+name);
			write(out.resolve("parsed-settings.tsv"), "UNAVAILABLE\n");
			write(out.resolve("counters.tsv"), "UNAVAILABLE\n");
			write(out.resolve("allocation.tsv"), "UNAVAILABLE\n");
		}else{
			final int expectedMink=(ordinal==4 ? 3 : 4); final boolean expectedShort=(ordinal==4 || ordinal==6 || ordinal==7);
			final int observedMink=number(owner, "mink").intValue(); final boolean observedShort=flag(owner, "useShortKmers"), observedMiddle=flag(owner, "maskMiddle");
			check(number(owner, "k").intValue()==4 && observedMink==expectedMink && observedShort==expectedShort && observedMiddle==!expectedShort && flag(owner, "ktrimN") && !flag(owner, "ktrimLeft") && !flag(owner, "ktrimRight") && !flag(owner, "ksplit"), "Constructor805/850 normalizes disabled values, preserves explicit usk=t and disables middle masking only for enabled short kmers; case="+name);
			final StringBuilder counters=new StringBuilder("field\tvalue\n"); boolean allZero=true;
			for(final String field : ZERO_COUNTERS){final long value=number(owner, field).longValue(); allZero&=value==0; counters.append(field).append('\t').append(value).append('\n');}
			write(out.resolve("counters.tsv"), counters.toString());
			check(allZero, "Public construction alone must not load references or process reads; case="+name);
			final AbstractKmerTable[] tables=(AbstractKmerTable[])field("keySets").get(owner);
			int allocated=0, empty=0;
			for(final AbstractKmerTable table : tables){if(table!=null){allocated++; if(table.size()==0){empty++;}}}
			final List<?> names=(List<?>)field("scaffoldNames").get(owner); final IntList lengths=(IntList)field("scaffoldLengths").get(owner);
			check(tables.length==7 && allocated==7 && empty==7 && names.size()==1 && "".equals(names.get(0)) && lengths.size==1 && lengths.get(0)==0 && field("ref").get(owner)==null && Arrays.equals((String[])field("literal").get(owner), new String[]{"ACGT"}), "Ordinary constructor preallocates seven empty tables and dummy scaffold0; reference loading belongs to process(), which is never called.");
			write(out.resolve("allocation.tsv"), "tables\tallocated\tempty\tscaffold_names\tscaffold_lengths\n"+tables.length+"\t"+allocated+"\t"+empty+"\t"+names.size()+"\t"+lengths.size+"\n");
			check(input.toString().equals(field("in1").get(owner)) && field("out1").get(owner)==null && field("out2").get(owner)==null && number(owner, "THREADS").intValue()==1 && number(owner, "maxReads").longValue()==1 && !flag(owner, "rcomp"), "Use one valid finite input, no biological output and one worker; normal constructor validation remains intact.");
			final StringBuilder parsed=new StringBuilder("field\tvalue\n");
			for(final String field : PARSED_FIELDS){parsed.append(field).append('\t').append(String.valueOf(field(field).get(owner))).append('\n');}
			write(out.resolve("parsed-settings.tsv"), parsed.toString());
			effective=String.valueOf(observedMink); shortKmers=String.valueOf(observedShort); middle=String.valueOf(observedMiddle);
		}
		require(checks==(rejected ? 4 : 6), "Fixed contract defines two common checks plus two rejection or four successful-state checks; actual="+checks);
		final String outcome=(rejected ? "PARSER_ASSERTION" : "CONSTRUCTED");
		write(out.resolve("case.tsv"), "scenario\toption\texplicit_usk\toutcome\teffective_mink\tuse_short_kmers\tmask_middle\tchecks\n"+name+"\t"+(OPTIONS[ordinal].isEmpty() ? "OMIT" : OPTIONS[ordinal])+"\t"+EXPLICIT_SHORT[ordinal]+"\t"+outcome+"\t"+effective+"\t"+shortKmers+"\t"+middle+"\t"+checks+"\n");
		final String summary="BBDUK_MINK_OPTION_FIXED_PASS version=1 scenario="+name+" constructorAttempts=1 constructed="+(rejected ? 0 : 1)+" parserAssertions="+(rejected ? 1 : 0)+" processCalls=0 checks="+checks+" failures=0\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
	}

	private String inputText(){return "@mink-option-"+name+"\nGGACGTGG\n+\nIIIIIIII\n";}
	private static Field field(final String name) throws NoSuchFieldException{final Field result=BBDuk.class.getDeclaredField(name); result.setAccessible(true); return result;}
	private static Number number(final BBDuk owner, final String name) throws Exception{return (Number)field(name).get(owner);}
	private static boolean flag(final BBDuk owner, final String name) throws Exception{return field(name).getBoolean(owner);}
	private void check(final boolean value, final String why){checks++; require(value, why);}
	private static void require(final boolean value, final String why){if(!value){throw new IllegalArgumentException(why);}}
	private static void write(final Path file, final String value) throws IOException{Files.write(file, value.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);}
	private static Path path(final String value){require(value!=null && !value.isEmpty() && value.indexOf('=')<0 && value.indexOf('\t')<0 && value.indexOf('\n')<0, "Use a nonempty path without equals, tabs or newlines for the legacy constructor."); return Paths.get(value).toAbsolutePath().normalize();}

	private final Path out, input;
	private final int ordinal;
	private final String name;
	private int checks;
	private static final String PARSER_MESSAGE="mink/kmin must be 1 to 31 when enabled; zero or negative values leave automatic short-kmer selection disabled; value=32";
	private static final String[] NAMES={"default", "negative", "zero", "alias_zero", "positive_three", "upper_32", "negative_explicit_short", "zero_explicit_short"};
	private static final String[] OPTIONS={"", "mink=-1", "mink=0", "kmin=0", "mink=3", "mink=32", "mink=-1", "mink=0"};
	private static final boolean[] EXPLICIT_SHORT={false, false, false, false, false, false, true, true};
	private static final boolean[] REJECTED={false, false, false, false, false, true, false, false};
	private static final String[] SETTINGS={"k=4", "kmask=N", "maskmiddle=t", "rcomp=f", "hdist=0", "edist=0", "qhdist=0", "qskip=1", "rskip=1", "speed=0", "rename=f", "minlen=1", "qtrim=f", "maxns=-1", "maq=0", "reads=1", "t=1", "simd=f", "qin=33", "array=t", "prealloc=f", "ow=f", "app=f"};
	private static final String[] ZERO_COUNTERS={"readsIn", "basesIn", "readsOut", "basesOut", "storedKmers", "refReads", "refBases", "refKmers", "readsKFiltered", "basesKFiltered", "readsKTrimmed", "basesKTrimmed", "readsQTrimmed", "basesQTrimmed", "readsFTrimmed", "basesFTrimmed", "readsQFiltered", "basesQFiltered", "readsNFiltered", "basesNFiltered", "readsEFiltered", "basesEFiltered", "readsPolyTrimmed", "basesPolyTrimmed", "readsTrimmedByOverlap", "basesTrimmedByOverlap", "readsTrimmedBySwift", "basesTrimmedBySwift", "badGcReads", "badGcBases", "badHeaderReads", "badHeaderBases"};
	private static final String[] PARSED_FIELDS={"k", "mink", "useShortKmers", "maskMiddle", "midMaskLen", "middleMask", "minlen", "minminlen", "minlen2", "ktrimN", "ktrimLeft", "ktrimRight", "ksplit", "THREADS", "maxReads", "rcomp", "in1", "out1", "out2"};
}
