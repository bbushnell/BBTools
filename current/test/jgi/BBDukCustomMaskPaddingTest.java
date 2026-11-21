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

import jgi.BBDuk;
import kmer.AbstractKmerTable;
import structures.IntList;

/** Four finite custom-symbol padding comparisons, one ordinary process() call per JVM.
 * Public arguments load literal ACGT without seeding tables, counters or globals.
 * Independent literals cover X internal/terminal padding, X zero padding and an
 * accepted positive-padding N control. X qualities are preserved, while N sets
 * qualities to zero. Complete FASTQ bytes, parsed settings, normal aggregated
 * counts and natural scaffold/table cleanup are checked without reparsing X.
 * Fully-covered, lowercase, short-kmer and negative-padding modes are excluded.
 * Baseline mode preserves the captured defects; fixed mode requires their repair.
 * @author Sucrose
 */
public final class BBDukCustomMaskPaddingTest{

	/** phase=baseline|fixed scenario=NAME outdir=NEW. Invoke once per fresh JVM. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final int split=arg.indexOf('=');
			require(split>0 && split<arg.length()-1, "Use nonempty flag=value arguments: "+arg);
			require(options.put(arg.substring(0, split), arg.substring(split+1))==null, "Duplicate option: "+arg);
		}
		require(options.size()==3 && options.keySet().containsAll(Arrays.asList("phase", "scenario", "outdir")), "Use exactly phase, scenario and outdir.");
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbdukcustommaskpaddingtest.sh with assertions enabled.");
		final String phase=options.get("phase");
		require("baseline".equals(phase) || "fixed".equals(phase), "Choose phase=baseline or phase=fixed.");
		final int ordinal=Arrays.asList(NAMES).indexOf(options.get("scenario"));
		require(ordinal>=0, "Choose one of the four named scenarios: "+Arrays.toString(NAMES));
		new BBDukCustomMaskPaddingTest(path(options.get("outdir")), phase, ordinal).run();
	}

	private BBDukCustomMaskPaddingTest(final Path out_, final String phase_, final int ordinal_) throws IOException{
		out=out_; phase=phase_; ordinal=ordinal_; name=NAMES[ordinal];
		Files.createDirectory(out); input=out.resolve("input.fq"); output=out.resolve("processed.fq");
		write(input, "@"+name+"\n"+INPUTS[ordinal]+"\n+\n"+qualityText(inputQuality())+"\n");
		write(out.resolve("settings.tsv"), String.join("\n", SETTINGS)+"\n");
		write(out.resolve("reflection.tsv"), "operation\tmember\nORDINARY_PUBLIC_CALL\tnew BBDuk(public arguments).process()\nREAD_ONLY\tBBDuk parsed settings, reference/scaffold setup and scalar statistics\nCOUNT_OBSERVATION\tbasesKTrimmed via kmask return to xsum to normal worker aggregation; NOT direct private return\nFIELD_WRITES\tNONE\n");
		write(out.resolve("observations.tsv"), "stage\treads_in\tbases_in\treads_out\tbases_out\tstored_kmers\tref_reads\tref_bases\tref_kmers\treads_ktrimmed\tbases_ktrimmed\tunrelated_counts\n");
		write(out.resolve("unrelated-counter-schema.tsv"), String.join("\n", ZERO_COUNTERS)+"\n");
		final StringBuilder oracle=new StringBuilder("scenario\tinput_bases\tinput_qualities\tsymbol\tpad\tbaseline_bases\tbaseline_qualities\tbaseline_caller_count\tintended_bases\tintended_qualities\tintended_caller_count\tcontrol\n");
		for(int i=0; i<NAMES.length; i++){oracle.append(NAMES[i]).append('\t').append(INPUTS[i]).append('\t').append(Arrays.toString(qualities(INPUTS[i].length()))).append('\t').append(SYMBOLS[i]).append('\t').append(PADS[i]).append('\t').append(BASELINE_BASES[i]).append('\t').append(Arrays.toString(BASELINE_QUALITY[i])).append('\t').append(BASELINE_COUNT[i]).append('\t').append(FIXED_BASES[i]).append('\t').append(Arrays.toString(FIXED_QUALITY[i])).append('\t').append(FIXED_COUNT[i]).append('\t').append(i>=2).append('\n');}
		write(out.resolve("oracles.tsv"), oracle.toString());
	}

	private void run() throws Exception{
		final boolean baseline="baseline".equals(phase);
		final String expectedBases=(baseline ? BASELINE_BASES[ordinal] : FIXED_BASES[ordinal]);
		final byte[] expectedQuality=(baseline ? BASELINE_QUALITY[ordinal] : FIXED_QUALITY[ordinal]);
		final int expectedCount=(baseline ? BASELINE_COUNT[ordinal] : FIXED_COUNT[ordinal]);
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(SETTINGS));
		args.add("kmask="+SYMBOLS[ordinal]); args.add("trimpad="+PADS[ordinal]); args.add("literal=ACGT"); args.add("in="+input); args.add("out="+output);
		write(out.resolve("constructor-args.txt"), String.join("\n", args)+"\n");
		final BBDuk owner=new BBDuk(args.toArray(new String[args.size()]));
		verifySettings(owner); observe(owner, "before", false);
		owner.process();
		observe(owner, "after", true);
		final long callerCount=number(owner, "basesKTrimmed").longValue();
		check(callerCount==expectedCount && number(owner, "readsKTrimmed").longValue()==1, "One read with only kmask active observes its returned count through the normal xsum/join aggregation; case="+name+" count="+callerCount);
		final List<String> lines=Files.readAllLines(output, StandardCharsets.US_ASCII);
		check(lines.size()==4 && lines.get(0).equals("@"+name) && lines.get(2).equals("+") && lines.get(1).length()==INPUTS[ordinal].length() && lines.get(3).length()==INPUTS[ordinal].length(), "Require one complete unpaired output FASTQ with unchanged ID and length; case="+name);
		final String actualBases=lines.get(1); final byte[] actualQuality=new byte[lines.get(3).length()];
		for(int i=0; i<actualQuality.length; i++){actualQuality[i]=(byte)(lines.get(3).charAt(i)-33);}
		check(actualBases.equals(expectedBases), "Output bases must equal the independent literal mask interval; case="+name+" actual="+actualBases);
		check(Arrays.equals(actualQuality, expectedQuality), "X qualities stay unchanged; only N replacements receive quality0; case="+name+" actual="+Arrays.toString(actualQuality));
		check(Arrays.equals(Files.readAllBytes(output), bytes("@"+name+"\n"+expectedBases+"\n+\n"+qualityText(expectedQuality)+"\n")) && Arrays.equals(Files.readAllBytes(input), bytes("@"+name+"\n"+INPUTS[ordinal]+"\n+\n"+qualityText(inputQuality())+"\n")), "Require exact output FASTQ bytes, final newline and unchanged input; do not reparse replacement symbols.");
		final boolean discrepancy=callerCount!=FIXED_COUNT[ordinal] || !actualBases.equals(FIXED_BASES[ordinal]) || !Arrays.equals(actualQuality, FIXED_QUALITY[ordinal]);
		check(discrepancy==(baseline && ordinal<2), baseline ? "Baseline has two X padding discrepancies; X zero-padding and accepted N controls retain their literals." : "Fixed mode must eliminate both X padding discrepancies and retain the two control literals.");
		int masked=0; for(int i=0; i<actualBases.length(); i++){if(actualBases.charAt(i)==SYMBOLS[ordinal]){masked++;}}
		write(out.resolve("case.tsv"), "scenario\tsymbol\tpad\tbases\tqualities\tcaller_count\tmasked_positions\tintended_bases\tintended_qualities\tintended_caller_count\tcontrol\tdiscrepancy\n"+name+"\t"+SYMBOLS[ordinal]+"\t"+PADS[ordinal]+"\t"+actualBases+"\t"+Arrays.toString(actualQuality)+"\t"+callerCount+"\t"+masked+"\t"+FIXED_BASES[ordinal]+"\t"+Arrays.toString(FIXED_QUALITY[ordinal])+"\t"+FIXED_COUNT[ordinal]+"\t"+(ordinal>=2)+"\t"+discrepancy+"\n");
		require(checks==20, "Source-derived contract is20 checks per independently launched case; actual="+checks);
		final String summary="BBDUK_CUSTOM_MASK_PADDING_"+(baseline ? "BASELINE" : "FIXED")+"_PASS version=1 phase="+phase+" scenario="+name+" owners=1 processCalls=1 callerCount="+callerCount+" maskedPositions="+masked+" controls="+(ordinal>=2 ? 1 : 0)+" discrepancies="+(discrepancy ? 1 : 0)+" checks="+checks+" failures=0\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
	}

	/** Read effective settings after public parsing, without seeding any owner state. */
	private void verifySettings(final BBDuk owner) throws Exception{
		check(flag(owner, "ktrimN") && !flag(owner, "ktrimLeft") && !flag(owner, "ktrimRight") && !flag(owner, "ksplit") && !flag(owner, "kmaskLowercase") && !flag(owner, "kmaskFullyCovered") && number(owner, "trimSymbol").intValue()==SYMBOLS[ordinal] && number(owner, "trimPad").intValue()==PADS[ordinal], "Use the exact public replacement character and padding; never reinterpret X as N.");
		check(number(owner, "k").intValue()==4 && number(owner, "mink").intValue()==4 && !flag(owner, "useShortKmers") && !flag(owner, "rcomp") && !flag(owner, "maskMiddle") && number(owner, "hammingDistance").intValue()==0 && number(owner, "editDistance").intValue()==0 && number(owner, "qHammingDistance").intValue()==0 && number(owner, "qSkip").intValue()==1 && number(owner, "minSkip").intValue()==1 && number(owner, "maxSkip").intValue()==1 && number(owner, "speed").intValue()==0, "Requested mink=-1 disables short kmers; BBDuk constructor normalizes effective mink to min(6,k)=4. Use exact forward k4 keys without middle-mask, speed or skipping.");
		check(number(owner, "THREADS").intValue()==1 && number(owner, "maxReads").longValue()==1 && number(owner, "minReadLength").intValue()==1 && number(owner, "maxReadLength").intValue()==Integer.MAX_VALUE && number(owner, "minLenFraction").doubleValue()==0 && !flag(owner, "rename") && !flag(owner, "skipR1") && !flag(owner, "skipR2") && number(owner, "restrictLeft").intValue()==0 && number(owner, "restrictRight").intValue()==0, "Process exactly one unpaired read without renaming, length rejection or restricted matching.");
		check(!flag(owner, "qtrimLeft") && !flag(owner, "qtrimRight") && !flag(owner, "trimClip") && number(owner, "minAvgQuality").doubleValue()==0 && number(owner, "minBaseQuality").intValue()==0 && number(owner, "maxNs").intValue()==-1 && number(owner, "maxNRate").doubleValue()==1 && number(owner, "minConsecutiveBases").intValue()<=0 && number(owner, "minBaseFrequency").doubleValue()==0, "Retain masked N/quality0 output; disable quality and undefined-base filtering/trimming.");
		check(number(owner, "forceTrimLeft").intValue()<=0 && number(owner, "forceTrimRight").intValue()<=0 && number(owner, "forceTrimRight2").intValue()<=0 && number(owner, "forceTrimModulo").intValue()<=0 && number(owner, "trimPolyA").intValue()<=0 && number(owner, "trimPolyGLeft").intValue()<=0 && number(owner, "trimPolyGRight").intValue()<=0 && number(owner, "trimPolyCLeft").intValue()<=0 && number(owner, "trimPolyCRight").intValue()<=0 && number(owner, "filterPolyG").intValue()<=0 && number(owner, "filterPolyC").intValue()<=0 && !flag(owner, "trimByOverlap") && !flag(owner, "swift") && !flag(owner, "ecc") && !flag(owner, "calcEntropy") && !flag(owner, "entropyMask") && !flag(owner, "entropyMark") && number(owner, "entropyTrim").intValue()==0 && !flag(owner, "filterGC"), "No force/polymer/overlap/Swift/ECC/entropy/GC operation may confound mask output or counters.");
		check(field("ref").get(owner)==null && Arrays.equals((String[])field("literal").get(owner), new String[]{"ACGT"}) && ((Object[])field("keySets").get(owner)).length==7 && ((List<?>)field("scaffoldNames").get(owner)).size()==1, "The ordinary constructor allocates seven tables and dummy scaffold0; normal process() must load the single literal as scaffold1.");
		final StringBuilder parsed=new StringBuilder("field\tvalue\n");
		for(final String name : PARSED_FIELDS){final Object value=field(name).get(owner); parsed.append(name).append('\t').append(value instanceof String[] ? Arrays.toString((String[])value) : String.valueOf(value)).append('\n');}
		write(out.resolve("parsed-settings.tsv"), parsed.toString());
	}

	private void observe(final BBDuk owner, final String stage, final boolean processed) throws Exception{
		final long readsIn=number(owner, "readsIn").longValue(), basesIn=number(owner, "basesIn").longValue(), readsOut=number(owner, "readsOut").longValue(), basesOut=number(owner, "basesOut").longValue();
		check(readsIn==(processed ? 1 : 0) && readsOut==(processed ? 1 : 0) && basesIn==(processed ? INPUTS[ordinal].length() : 0) && basesOut==(processed ? INPUTS[ordinal].length() : 0), "The sole normal input read must be emitted at unchanged length; stage="+stage);
		final long stored=number(owner, "storedKmers").longValue(), refReads=number(owner, "refReads").longValue(), refBases=number(owner, "refBases").longValue(), refKmers=number(owner, "refKmers").longValue();
		check(stored==(processed ? 1 : 0) && refReads==(processed ? 1 : 0) && refBases==(processed ? 4 : 0) && refKmers==(processed ? 1 : 0), "Normal literal loading, not seeded tables/counts, must produce exactly one k4 reference key; stage="+stage);
		final long[] unrelated=new long[ZERO_COUNTERS.length]; boolean zero=true;
		for(int i=0; i<unrelated.length; i++){unrelated[i]=number(owner, ZERO_COUNTERS[i]).longValue(); zero&=unrelated[i]==0;}
		check(zero && (processed || (number(owner, "readsKTrimmed").longValue()==0 && number(owner, "basesKTrimmed").longValue()==0)), "Unrelated trimming/filtering counters and initial mask statistics must stay zero; stage="+stage+" values="+Arrays.toString(unrelated));
		observeScaffolds(owner, stage, processed);
		append(out.resolve("observations.tsv"), stage+"\t"+readsIn+"\t"+basesIn+"\t"+readsOut+"\t"+basesOut+"\t"+stored+"\t"+refReads+"\t"+refBases+"\t"+refKmers+"\t"+number(owner, "readsKTrimmed")+"\t"+number(owner, "basesKTrimmed")+"\t"+Arrays.toString(unrelated)+"\n");
	}
	/** Observe normal table/scaffold release, never disable RELEASE_TABLES to retain state. */
	private void observeScaffolds(final BBDuk owner, final String stage, final boolean processed) throws Exception{
		final AbstractKmerTable[] tables=(AbstractKmerTable[])field("keySets").get(owner);
		boolean tableState=tables!=null && tables.length==7;
		if(tableState){for(final AbstractKmerTable table : tables){tableState&=(processed ? table==null : table!=null && table.size()==0);}}
		final List<?> names=(List<?>)field("scaffoldNames").get(owner);
		final IntList lengths=(IntList)field("scaffoldLengths").get(owner);
		final boolean scaffoldState=names!=null && (processed ? names.isEmpty() && lengths==null : names.size()==1 && "".equals(names.get(0)) && lengths!=null && lengths.size==1 && lengths.get(0)==0);
		check(flag(owner, "RELEASE_TABLES") && tableState && scaffoldState && Arrays.equals((int[])field("refScafCounts").get(owner), new int[]{processed ? 1 : 0}) && field("scaffoldReadCounts").get(owner)==null && field("scaffoldBaseCounts").get(owner)==null && !owner.errorState, "BBDuk process2 releases tables/scaffolds normally at1264/1276; constructor has dummy scaffold only, final refScafCounts retains the one literal; stage="+stage);
		write(out.resolve("scaffolds-"+stage+".tsv"), "stage\ttables\tnames\tlengths\tref_scaffolds\tatomic_counts\terror\n"+stage+"\t"+(processed ? "7_NULL\tEMPTY\tNULL\t[1]" : "7_EMPTY\tDUMMY\t[0]\t[0]")+"\tNULL\tfalse\n");
	}
	private byte[] inputQuality(){return qualities(INPUTS[ordinal].length());}
	private static byte[] qualities(final int length){assert(length>0 && length<=64) : "Literal qualities start at30 and Phred+33 permits at most93; length="+length; final byte[] result=new byte[length]; for(int i=0; i<length; i++){result[i]=(byte)(30+i);} return result;}
	private static String qualityText(final byte[] quality){final StringBuilder result=new StringBuilder(); for(final byte value : quality){assert(value>=0 && value<=93) : "Finite constructor FASTQ uses Phred+33; quality="+value; result.append((char)(value+33));} return result.toString();}
	private static Field field(final String name) throws NoSuchFieldException{final Field result=BBDuk.class.getDeclaredField(name); result.setAccessible(true); return result;}
	private static Number number(final BBDuk owner, final String name) throws Exception{return (Number)field(name).get(owner);}
	private static boolean flag(final BBDuk owner, final String name) throws Exception{return field(name).getBoolean(owner);}
	private void check(final boolean value, final String why){checks++; require(value, why);}
	private static void require(final boolean value, final String why){if(!value){throw new IllegalArgumentException(why);}}
	private static Path path(final String value){require(value!=null && !value.isEmpty() && value.indexOf('=')<0 && value.indexOf('\t')<0 && value.indexOf('\n')<0, "Use a nonempty path without equals,tabs,newlines for the legacy constructor."); return Paths.get(value).toAbsolutePath().normalize();}
	private static byte[] bytes(final String value){return value.getBytes(StandardCharsets.US_ASCII);}
	private static void write(final Path file, final String value) throws IOException{Files.write(file, bytes(value), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path file, final String value) throws IOException{Files.write(file, bytes(value), StandardOpenOption.APPEND);}

	private final Path out, input, output;
	private final String phase, name;
	private final int ordinal;
	private int checks;
	private static final String[] NAMES={"x_internal_pad1", "x_terminal_pad1", "x_internal_pad0", "n_internal_pad1"};
	private static final char[] SYMBOLS={'X', 'X', 'X', 'N'};
	private static final String[] INPUTS={"GGACGTGG", "GGACGT", "GGACGTGG", "GGACGTGG"};
	private static final int[] PADS={1, 1, 0, 1}, BASELINE_COUNT={4, 4, 4, 6}, FIXED_COUNT={6, 5, 4, 6};
	private static final String[] BASELINE_BASES={"GGAXXXXG", "GGAXXX", "GGXXXXGG", "GNNNNNNG"}, FIXED_BASES={"GXXXXXXG", "GXXXXX", "GGXXXXGG", "GNNNNNNG"};
	private static final byte[][] BASELINE_QUALITY={{30, 31, 32, 33, 34, 35, 36, 37}, {30, 31, 32, 33, 34, 35}, {30, 31, 32, 33, 34, 35, 36, 37}, {30, 0, 0, 0, 0, 0, 0, 37}};
	private static final byte[][] FIXED_QUALITY={{30, 31, 32, 33, 34, 35, 36, 37}, {30, 31, 32, 33, 34, 35}, {30, 31, 32, 33, 34, 35, 36, 37}, {30, 0, 0, 0, 0, 0, 0, 37}};
	private static final String[] SETTINGS={"k=4", "mink=-1", "maskfullycovered=f", "rcomp=f", "maskmiddle=f", "hdist=0", "edist=0", "qhdist=0", "qskip=1", "rskip=1", "speed=0", "rename=f", "minlen=1", "qtrim=f", "maxns=-1", "maq=0", "reads=1", "t=1", "simd=f", "qin=33", "qout=33", "array=t", "prealloc=f", "ow=f", "app=f"};
	private static final String[] ZERO_COUNTERS={"readsKFiltered", "basesKFiltered", "readsQTrimmed", "basesQTrimmed", "readsFTrimmed", "basesFTrimmed", "readsQFiltered", "basesQFiltered", "readsNFiltered", "basesNFiltered", "readsEFiltered", "basesEFiltered", "readsPolyTrimmed", "basesPolyTrimmed", "readsTrimmedByOverlap", "basesTrimmedByOverlap", "readsTrimmedBySwift", "basesTrimmedBySwift", "badGcReads", "badGcBases", "badHeaderReads", "badHeaderBases"};
	private static final String[] PARSED_FIELDS={"k", "mink", "ktrimN", "ktrimLeft", "ktrimRight", "ksplit", "trimSymbol", "trimPad", "kmaskLowercase", "kmaskFullyCovered", "useShortKmers", "rcomp", "maskMiddle", "hammingDistance", "editDistance", "qHammingDistance", "qSkip", "minSkip", "maxSkip", "speed", "THREADS", "maxReads", "minReadLength", "maxReadLength", "minLenFraction", "rename", "skipR1", "skipR2", "restrictLeft", "restrictRight", "qtrimLeft", "qtrimRight", "trimClip", "minAvgQuality", "minBaseQuality", "maxNs", "maxNRate", "minConsecutiveBases", "minBaseFrequency", "forceTrimLeft", "forceTrimRight", "forceTrimRight2", "forceTrimModulo", "trimPolyA", "trimPolyGLeft", "trimPolyGRight", "trimPolyCLeft", "trimPolyCRight", "filterPolyG", "filterPolyC", "trimByOverlap", "swift", "ecc", "calcEntropy", "entropyMask", "entropyMark", "entropyTrim", "filterGC", "literal"};
}
