package test.jgi;

import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;

import jgi.BBDuk;
import json.JsonObject;
import json.JsonParser;
import kmer.AbstractKmerTable;
import shared.Shared;
import structures.IntList;

/** One repaired JSON-key observation per fresh constructor and reporter invocation.
 * Public json=t construction initializes ordinary metadata. Exactly one reflective
 * call invokes the real private toJson(long), which naturally adds reporting members
 * to the owner's jsonStats object. This mutation is recorded, never seeded or reset.
 * Complete strict validation precedes typed native parsing. Keys retain punctuation
 * and trailing spaces; raw reports and all metadata values remain available. Time
 * varies, while commandLine contains the exact case path; no key normalization occurs.
 * No process, stream factory, counter seeding, owner-field write or global reset.
 * Zero-count serialization checks do not establish biological or lifecycle behavior.
 * @author Sucrose
 */
public final class BBDukJsonKeyTest{

	/** Use scenario=mask|left|right|filter and a new outdir via the dev launcher. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final int split=arg.indexOf('='); require(split>0 && split<arg.length()-1, "Use nonempty flag=value arguments: "+arg);
			require(options.put(arg.substring(0, split), arg.substring(split+1))==null, "Duplicate option: "+arg);
		}
		require(options.size()==2 && options.keySet().containsAll(Arrays.asList("scenario", "outdir")), "Use exactly scenario and outdir.");
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbdukjsonkeytest.sh with assertions enabled for normal constructor validation.");
		final int ordinal=Arrays.asList(NAMES).indexOf(options.get("scenario")); require(ordinal>=0, "Choose one of four baseline modes: "+Arrays.toString(NAMES));
		new BBDukJsonKeyTest(path(options.get("outdir")), ordinal).run();
	}

	private BBDukJsonKeyTest(final Path out_, final int ordinal_) throws IOException{
		out=out_; ordinal=ordinal_; name=NAMES[ordinal]; input=out.resolve("input.fq"); Files.createDirectory(out); write(input, inputText());
		write(out.resolve("settings.tsv"), lines(SETTINGS)); write(out.resolve("counter-schema.tsv"), lines(ZERO_COUNTERS));
		write(out.resolve("metadata-policy.tsv"), METADATA_POLICY);
		write(out.resolve("reflection.tsv"), "operation\tmember\nORDINARY_PUBLIC_CALL\tnew BBDuk(public arguments) only\nREPORT_CALL\tprivate BBDuk.toJson(long), Method.invoke exactly once\nNATURAL_MUTATION\tReporter adds members to the existing jsonStats object\nPROCESS_CALLS\t0\nFACTORY_CALLS\t0\nFIELD_WRITES\tNONE\nGLOBAL_RESETS\tNONE; fresh JVM per case\n");
	}

	private void run() throws Exception{
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(SETTINGS));
		if(ordinal<3){args.add(OPTIONS[ordinal]);} args.add("literal=ACGT"); args.add("in="+input);
		final String[] raw=args.toArray(new String[args.size()]), passed=raw.clone(); write(out.resolve("constructor-args.txt"), lines(raw));
		final BBDuk owner=new BBDuk(passed);
		check(number(owner, "k")==4 && number(owner, "mink")==4 && !flag(owner, "useShortKmers") && !flag(owner, "maskMiddle") && flag(owner, "ktrimN")==(ordinal==0) && flag(owner, "ktrimLeft")==(ordinal==1) && flag(owner, "ktrimRight")==(ordinal==2) && !flag(owner, "ksplit") && flag(owner, "json") && number(owner, "THREADS")==1 && number(owner, "maxReads")==1, "Each ordinary constructor must select exactly its requested mode and finite one-worker settings; case="+name);
		check(emptyState(owner), "Construction must leave all counters zero, seven tables empty and only dummy scaffold0; case="+name);
		final JsonObject stats=(JsonObject)field("jsonStats").get(owner);
		check(stats!=null && stats.omap!=null && stats.jmapSize()==0 && stats.omap.keySet().equals(new TreeSet<String>(Arrays.asList(METADATA))), "json=t PreParser initializes only the six ordinary scalar metadata members; case="+name);
		final Map<String, Object> before=new LinkedHashMap<String, Object>(stats.omap);
		final String beforeText=stats.toString(); write(out.resolve("metadata-before.json"), beforeText+"\n");
		final JsonObject beforeParsed=parseStrict(beforeText);
		check(validMetadata(beforeParsed, raw), "Native parsed metadata must retain exact command line, version, memory, assertions, Java version and JVM argument array types/values; case="+name);
		write(out.resolve("metadata-before.tsv"), typedTable(beforeParsed));
		final Method reporter=BBDuk.class.getDeclaredMethod("toJson", long.class); reporter.setAccessible(true);
		final long start=System.nanoTime();
		final String report=(String)reporter.invoke(owner, start); reporterInvocations++;
		write(out.resolve("report.json"), report+"\n");
		final JsonParser validator=new JsonParser(report);
		check(validator.validate() && !validator.errorState(), "Validate the complete returned JSON, including its final byte, before relying on permissive typed parsing; case="+name);
		final JsonParser parser=new JsonParser(report); final JsonObject parsed=parser.parseJsonObject();
		check(parsed!=null && !parser.errorState() && validReport(parsed, raw), "Native typed keys/values must exactly match the repaired report schema, including the matching bases suffix and existing reads punctuation; case="+name);
		check(!new JsonParser(report+" trailing").validate(), "The full-document gate must reject trailing text that parseJsonObject alone permits; case="+name);
		check(stats==field("jsonStats").get(owner) && stats.jmapSize()==0 && stats.omap.keySet().equals(expectedKeys()) && preservedMetadata(before, stats.omap), "The one real reporter call naturally adds eleven members to the same object and preserves all six original metadata values; case="+name);
		check(reporterInvocations==1 && emptyState(owner) && Arrays.equals(raw, passed) && Arrays.equals(Files.readAllBytes(input), inputText().getBytes(StandardCharsets.US_ASCII)) && input.toString().equals(field("in1").get(owner)) && field("in2").get(owner)==null && field("out1").get(owner)==null && field("out2").get(owner)==null, "Reporting must leave inputs, arguments, zero counters and ordinary allocation unchanged without process/factory calls; case="+name);
		write(out.resolve("constructor-args-after.txt"), lines(passed)); write(out.resolve("typed-keys.tsv"), typedTable(parsed));
		write(out.resolve("report-pins.tsv"), "file\tsha80\nmetadata-before.json\t"+sha80(beforeText+"\n")+"\nreport.json\t"+sha80(report+"\n")+"\n");
		final StringBuilder counters=new StringBuilder("field\tvalue\n"); for(final String member : ZERO_COUNTERS){counters.append(member).append('\t').append(number(owner, member)).append('\n');} write(out.resolve("counters.tsv"), counters.toString());
		write(out.resolve("allocation.tsv"), "tables\tallocated\tempty\tscaffold_names\tscaffold_lengths\n7\t7\t7\t1\t1\n");
		write(out.resolve("mutation.tsv"), "object_identity\tmetadata_before\tmembers_after\tmetadata_preserved\tadded_members\treporter_calls\nSAME\t6\t17\t6\t11\t"+reporterInvocations+"\n");
		write(out.resolve("validation.tsv"), "strict_complete\ttyped_parse_errors\ttrailing_rejected\ntrue\tfalse\ttrue\n");
		require(checks==9 && reporterInvocations==1, "Nine contracts and exactly one actual toJson invocation are required; checks="+checks+" calls="+reporterInvocations);
		write(out.resolve("case.tsv"), "scenario\tmode\treads_key_json\tintended_bases_key_json\tobserved_bases_key_json\tchecks\n"+name+"\t"+MODES[ordinal]+"\t"+jsonValue(READ_KEYS[ordinal])+"\t"+jsonValue(BASE_KEYS[ordinal])+"\t"+jsonValue(BASE_KEYS[ordinal])+"\t"+checks+"\n");
		final String summary="BBDUK_JSON_KEY_FIXED_PASS version=1 scenario="+name+" constructed=1 reporterCalls=1 wrongBasesKeys=0 processCalls=0 factoryCalls=0 checks="+checks+" failures=0\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
	}

	private boolean validReport(final JsonObject parsed, final String[] raw){
		if(parsed.omap==null || parsed.jmapSize()!=0 || !parsed.omap.keySet().equals(expectedKeys()) || !validMetadata(parsed, raw)){return false;}
		if(!MODES[ordinal].equals(parsed.omap.get("mode")) || !Long.valueOf(4).equals(parsed.omap.get("k"))){return false;}
		for(final String member : ZERO_MEMBERS){if(!Long.valueOf(0).equals(parsed.omap.get(member))){return false;}}
		if(!Long.valueOf(0).equals(parsed.omap.get(READ_KEYS[ordinal])) || !Long.valueOf(0).equals(parsed.omap.get(BASE_KEYS[ordinal]))){return false;}
		final Object elapsed=parsed.omap.get("time"); return elapsed instanceof Double && Double.isFinite((Double)elapsed) && (Double)elapsed>=0;
	}
	private Set<String> expectedKeys(){final Set<String> result=new TreeSet<String>(Arrays.asList(METADATA)); result.addAll(Arrays.asList("k", "mode", "time")); result.addAll(Arrays.asList(ZERO_MEMBERS)); result.add(READ_KEYS[ordinal]); result.add(BASE_KEYS[ordinal]); return result;}
	private static boolean validMetadata(final JsonObject parsed, final String[] raw){
		final Map<String, Object> map=parsed.omap; if(map==null){return false;}
		final Object memory=map.get("memory"), javaVersion=map.get("javaVersion"), jvm=map.get("JVM_args");
		return ("jgi.BBDuk "+Arrays.toString(raw)).equals(map.get("commandLine")) && Shared.BBTOOLS_VERSION_STRING.equals(map.get("version")) && Boolean.TRUE.equals(map.get("assertions")) && memory instanceof Long && ((Long)memory).longValue()==Shared.memTotal() && javaVersion instanceof Double && ((Double)javaVersion).doubleValue()==Shared.javaVersion && jvm instanceof Object[] && Arrays.equals((Object[])jvm, Shared.JVM_ARGS().toArray());
	}
	private static boolean preservedMetadata(final Map<String, Object> before, final Map<String, Object> after){for(final String key : METADATA){if(!after.containsKey(key) || !Objects.deepEquals(before.get(key), after.get(key))){return false;}} return before.size()==6;}
	private static JsonObject parseStrict(final String text){final JsonParser validator=new JsonParser(text); require(validator.validate() && !validator.errorState(), "Metadata must be complete standard JSON."); final JsonParser parser=new JsonParser(text); final JsonObject result=parser.parseJsonObject(); require(result!=null && !parser.errorState(), "Strictly validated metadata must parse without errors."); return result;}
	private static String typedTable(final JsonObject parsed){final StringBuilder result=new StringBuilder("key_json\ttype\tvalue_json\n"); for(final String key : new TreeSet<String>(parsed.omap.keySet())){final Object value=parsed.omap.get(key); result.append(jsonValue(key)).append('\t').append(value==null ? "null" : value.getClass().getSimpleName()).append('\t').append(jsonValue(value)).append('\n');} return result.toString();}
	private static String jsonValue(final Object value){final String object=new JsonObject("v", value).toText(null, 0, true).toString(); return object.substring(6, object.length()-1);}
	/** Binds retained UTF-8 report bytes to the strings strictly validated above. */
	private static String sha80(final String value) throws Exception{
		final byte[] digest=MessageDigest.getInstance("SHA-256").digest(value.getBytes(StandardCharsets.UTF_8)); final StringBuilder result=new StringBuilder(20);
		for(int i=digest.length-10; i<digest.length; i++){final int b=digest[i]&255; result.append("0123456789abcdef".charAt(b>>>4)).append("0123456789abcdef".charAt(b&15));} return result.toString();
	}
	private static boolean emptyState(final BBDuk owner) throws Exception{
		for(final String member : ZERO_COUNTERS){if(number(owner, member)!=0){return false;}}
		final AbstractKmerTable[] tables=(AbstractKmerTable[])field("keySets").get(owner); if(tables.length!=7){return false;}
		for(final AbstractKmerTable table : tables){if(table==null || table.size()!=0){return false;}}
		final List<?> names=(List<?>)field("scaffoldNames").get(owner); final IntList lengths=(IntList)field("scaffoldLengths").get(owner);
		return names.size()==1 && "".equals(names.get(0)) && lengths.size==1 && lengths.get(0)==0 && field("ref").get(owner)==null && Arrays.equals((String[])field("literal").get(owner), new String[]{"ACGT"});
	}
	private String inputText(){return "@json-key-"+name+"\nGGACGTGG\n+\nIIIIIIII\n";}
	private static Field field(final String name) throws NoSuchFieldException{final Field result=BBDuk.class.getDeclaredField(name); result.setAccessible(true); return result;}
	private static long number(final BBDuk owner, final String name) throws Exception{return ((Number)field(name).get(owner)).longValue();}
	private static boolean flag(final BBDuk owner, final String name) throws Exception{return field(name).getBoolean(owner);}
	private void check(final boolean value, final String why){checks++; require(value, why);}
	private static void require(final boolean value, final String why){if(!value){throw new IllegalArgumentException(why);}}
	private static String lines(final String[] values){return String.join("\n", values)+"\n";}
	private static void write(final Path file, final String value) throws IOException{Files.write(file, value.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);}
	private static Path path(final String value){require(value!=null && !value.isEmpty() && value.indexOf('=')<0 && value.indexOf('\t')<0 && value.indexOf('\n')<0, "Use a nonempty directory without equals, tabs or newlines for reproducible JSON fixtures."); return Paths.get(value).toAbsolutePath().normalize();}

	private final Path out, input;
	private final int ordinal;
	private final String name;
	private int checks, reporterInvocations;
	private static final String[] NAMES={"mask", "left", "right", "filter"}, OPTIONS={"kmask=N", "ktrim=l", "ktrim=r"}, MODES={"ktrimN", "ktrimLeft", "ktrimRight", "kFilter"};
	private static final String[] READ_KEYS={"readsKMasked: ", "readsKTrimmed:", "readsKTrimmed:", "readsKFiltered"}, BASE_KEYS={"basesKMasked: ", "basesKTrimmed:", "basesKTrimmed:", "basesKFiltered"};
	private static final String[] METADATA={"commandLine", "version", "memory", "assertions", "javaVersion", "JVM_args"};
	private static final String[] ZERO_MEMBERS={"readsIn", "basesIn", "readsRemoved", "basesRemoved", "readsOut", "basesOut"};
	private static final String[] SETTINGS={"k=4", "mink=-1", "maskmiddle=f", "rcomp=f", "hdist=0", "edist=0", "qhdist=0", "qskip=1", "rskip=1", "speed=0", "rename=f", "minlen=1", "qtrim=f", "maxns=-1", "maq=0", "reads=1", "t=1", "simd=f", "qin=33", "array=t", "prealloc=f", "ow=f", "app=f", "json=t"};
	private static final String[] ZERO_COUNTERS={"readsIn", "basesIn", "readsOut", "basesOut", "storedKmers", "refReads", "refBases", "refKmers", "readsKFiltered", "basesKFiltered", "readsKTrimmed", "basesKTrimmed", "readsQTrimmed", "basesQTrimmed", "readsFTrimmed", "basesFTrimmed", "readsQFiltered", "basesQFiltered", "readsNFiltered", "basesNFiltered", "readsEFiltered", "basesEFiltered", "readsPolyTrimmed", "basesPolyTrimmed", "readsTrimmedByOverlap", "basesTrimmedByOverlap", "readsTrimmedBySwift", "basesTrimmedBySwift", "badGcReads", "badGcBases", "badHeaderReads", "badHeaderBases"};
	private static final String METADATA_POLICY="key\thandling\ncommandLine\tExact constructor arguments retained; only case directory varies between runs\nversion\tExact Shared.BBTOOLS_VERSION_STRING retained\nmemory\tExact Shared.memTotal retained; no cross-run numeric normalization\nassertions\tBoolean true required and retained\njavaVersion\tExact Shared.javaVersion retained\nJVM_args\tExact Shared.JVM_ARGS collection retained as parsed Object[]\ntime\tFinite nonnegative Double retained; elapsed value is not compared across runs\n";
}
