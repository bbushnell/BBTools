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
import java.util.Map;

import jgi.BBDuk;
import kmer.AbstractKmerTable;
import kmer.HashArray;
import stream.Read;

/** Finite index-only check of legacy BBDuk's primary versus short edit modes.
 * Five cases compare literal complete sets over every length3/4 DNA word and all
 * seven shards. Two length4 witnesses check recursive insertion at budget2 and
 * its absence at budget1. Expected words are literal; no mutator or edit-distance
 * implementation generates the oracle. Each case constructs its own ordinary
 * owner and seven never-started loaders, invokes actual addToMap(Read,0) once per
 * loader, and then releases its owner references. Only normal constructor table
 * allocators may start and join. No queues, processing workers, reference streams,
 * reports, final/configuration/counter seeds or executable production changes.
 * @author Sucrose
 */
public final class BBDukShortEditModeTest{

	/** phase=baseline|fixed outdir=NEW; fixed requires reference=ORIGINAL_BASELINE_CASES. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final int split=arg.indexOf('=');
			require(split>0 && split<arg.length()-1, "Use documented nonempty flag=value arguments: "+arg);
			final String name=arg.substring(0, split), value=arg.substring(split+1);
			require(Arrays.asList("phase", "outdir", "reference").contains(name) && options.put(name, value)==null, "Unknown or duplicate flag: "+name);
		}
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbdukshorteditmodetest.sh with assertions enabled.");
		final String phase=options.get("phase");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		final Path reference=options.containsKey("reference") ? path(options.get("reference")) : null;
		require("baseline".equals(phase) ? reference==null : reference!=null, "Only fixed requires the original baseline fixture directory.");
		if(reference!=null){require(text(reference.resolve("summary.txt")).matches("BBDUK_SHORT_EDIT_MODE_PASS phase=baseline version=1 scenarios=7 owners=7 loaders=49 calls=49 completeCases=5 targetedCases=2 controls=4 discrepancies=3 checks=[0-9]+ failures=0\\n"), "Reference must be the successful original short-edit baseline fixture directory.");}
		final BBDukShortEditModeTest test=new BBDukShortEditModeTest(path(options.get("outdir")), reference, phase);
		try{test.exercise();}finally{test.saveDispatch();}
		test.finish();
	}

	private BBDukShortEditModeTest(final Path out_, final Path reference_, final String phase_) throws Exception{
		out=out_; reference=reference_; phase=phase_; baseline=phase.equals("baseline"); Files.createDirectory(out);
		final Class<?> type=Class.forName("jgi.BBDuk$LoadThread");
		constructor=type.getDeclaredConstructor(BBDuk.class, int.class); constructor.setAccessible(true);
		add=type.getDeclaredMethod("addToMap", Read.class, int.class); add.setAccessible(true);
		refReads=field(type, "refReadsT"); refBases=field(type, "refBasesT"); refKmers=field(type, "refKmersT"); added=field(type, "addedT");
		write(out.resolve("settings.tsv"), String.join("\n", SETTINGS)+"\n");
		write(out.resolve("reflection.tsv"), "operation\tmember\nREAD_ONLY\tBBDuk.keySets,storedKmers,k,mink,editDistance,editDistance2,hammingDistance,hammingDistance2\nCONSTRUCT_ONLY\tBBDuk.LoadThread\nINVOKE_SYNC\tLoadThread.addToMap(Read,int)\nOBSERVE\tLoadThread.refReadsT,refBasesT,refKmersT,addedT\nFIELD_WRITES\tNONE\n");
		write(out.resolve("cases.tsv"), "scenario\tkind\ttarget\tpresent\tintended_present\tobserved_keys\texpected_keys\tloaders\tpayloads\n");
		write(out.resolve("controls.tsv"), "scenario\ttarget\tpresent\tcomplete_keys\tloaders\tpayloads\n");
		final StringBuilder inputs=new StringBuilder("scenario\treference\tk\tmink\tedist\tedist2\thdist\thdist2\ttarget\tbaseline_present\tfixed_present\n");
		for(int i=0; i<SPECS.length; i++){
			final Spec spec=SPECS[i]; statuses[i]="SKIPPED";
			inputs.append(spec.name).append('\t').append(spec.bases).append('\t').append(spec.bases.length()).append('\t').append(spec.bases.length()-1).append('\t').append(spec.edist).append('\t').append(spec.edist2).append("\t0\t").append(spec.hdist2).append('\t').append(spec.target).append('\t').append(spec.before).append('\t').append(spec.after).append('\n');
		}
		write(out.resolve("inputs.tsv"), inputs.toString());
	}

	/** Each invocation retains its owner/loaders only until this case returns. */
	private void exercise() throws Exception{
		for(int i=0; i<SPECS.length; i++){statuses[i]="FAILED"; exercise(SPECS[i], i); statuses[i]="PASSED";}
		check(owners==7 && loaders==49 && calls==49 && controls==4 && discrepancies==(baseline ? 3 : 0), "Keep seven separate owners,49 NEW loaders/calls,four controls and exactly the phase's three mode discrepancies.");
		if(reference!=null){
			for(final String file : new String[]{"settings.tsv", "reflection.tsv", "inputs.tsv", "controls.tsv"}){same(file);}
			for(final Spec spec : SPECS){
				for(final String file : new String[]{"constructor-input.fq", "setup.tsv", "baseline-expected.txt", "fixed-expected.txt"}){same(spec.name+"/"+file);}
				if(spec.before==spec.after){same(spec.name+"/keys.tsv");}
			}
		}
	}

	private void exercise(final Spec spec, final int ordinal) throws Exception{
		final Path directory=out.resolve(spec.name); Files.createDirectory(directory);
		final String qualities=repeat('I', spec.bases.length());
		write(directory.resolve("constructor-input.fq"), "@constructor_input\n"+spec.bases+"\n+\n"+qualities+"\n");
		write(directory.resolve("setup.tsv"), "id\t"+spec.name+"\nnumeric_id\t"+ordinal+"\nscaffold_id\t1\nbases\t"+spec.bases+"\nqualities\t"+qualities+"\nk\t"+spec.bases.length()+"\nmink\t"+(spec.bases.length()-1)+"\nedist\t"+spec.edist+"\nedist2\t"+spec.edist2+"\nhdist\t0\nhdist2\t"+spec.hdist2+"\n");
		write(directory.resolve("baseline-expected.txt"), spec.complete ? spec.full+" "+spec.shortBefore+"\n" : spec.target+"\t"+spec.before+"\n");
		write(directory.resolve("fixed-expected.txt"), spec.complete ? spec.full+" "+spec.shortAfter+"\n" : spec.target+"\t"+spec.after+"\n");
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(SETTINGS));
		args.add("k="+spec.bases.length()); args.add("mink="+(spec.bases.length()-1)); args.add("edist="+spec.edist); args.add("edist2="+spec.edist2); args.add("hdist=0"); args.add("hdist2="+spec.hdist2);
		args.add("literal="+spec.bases); args.add("in="+directory.resolve("constructor-input.fq"));
		write(directory.resolve("constructor-args.txt"), String.join("\n", args)+"\n");
		final BBDuk tool=new BBDuk(args.toArray(new String[args.size()])); owners++;
		final AbstractKmerTable[] sets=(AbstractKmerTable[])field(BBDuk.class, "keySets").get(tool);
		check(sets.length==7 && field(BBDuk.class, "k").getInt(tool)==spec.bases.length() && field(BBDuk.class, "mink").getInt(tool)==spec.bases.length()-1, "Normal constructor must select seven shards and this case's k/mink.");
		check(field(BBDuk.class, "editDistance").getInt(tool)==spec.edist && field(BBDuk.class, "editDistance2").getInt(tool)==spec.edist2 && field(BBDuk.class, "hammingDistance").getInt(tool)==spec.edist && field(BBDuk.class, "hammingDistance2").getInt(tool)==Math.max(spec.hdist2, spec.edist2), "Normal parser promotes each Hamming budget independently; never seed configuration.");
		check(field(BBDuk.class, "storedKmers").getLong(tool)==0, "Direct-loader fixture starts before owner aggregation; no reference pipeline may run.");
		final Thread[] workers=new Thread[7];
		for(int shard=0; shard<7; shard++){
			check(sets[shard] instanceof HashArray, "Ordinary array=t creates HashArray tables with explicit primary/victim size APIs.");
			check(tableSize(sets[shard])==0, "Fresh constructor tables must be empty including victim forests.");
			workers[shard]=(Thread)constructor.newInstance(tool, Integer.valueOf(shard)); loaders++;
			worker(workers[shard], 0, 0, 0);
		}
		final byte[] expectedBases=bytes(spec.bases), expectedQuality=new byte[expectedBases.length]; Arrays.fill(expectedQuality, (byte)40);
		final Read read=new Read(expectedBases.clone(), expectedQuality.clone(), spec.name, ordinal); read.obj=Integer.valueOf(1);
		write(directory.resolve("loaders.tsv"), "shard\treturned_new_keys\tprimary_keys\tvictim_keys\tstate\tref_reads\tref_bases\tref_kmers\tadded_field\n");
		long returned=0, total=0;
		for(int shard=0; shard<7; shard++){
			final long result=((Long)add.invoke(workers[shard], read, Integer.valueOf(0))).longValue(); calls++; returned+=result;
			check(result>=0, "setIfNotPresent reports nonnegative new-key counts for a finite direct load.");
			worker(workers[shard], 1, expectedBases.length, 1);
			check(read.mate==null && spec.name.equals(read.id) && read.numericID==ordinal && Integer.valueOf(1).equals(read.obj) && Arrays.equals(read.bases, expectedBases) && Arrays.equals(read.quality, expectedQuality), "Loading must preserve the independent reference payload and assigned scaffold ID.");
			final HashArray table=(HashArray)sets[shard]; total+=tableSize(table);
			append(directory.resolve("loaders.tsv"), shard+"\t"+result+"\t"+table.size()+"\t"+table.victims().size()+"\t"+workers[shard].getState()+"\t"+refReads.getLong(workers[shard])+"\t"+refBases.getLong(workers[shard])+"\t"+refKmers.getLong(workers[shard])+"\t"+added.getLong(workers[shard])+"\n");
		}
		check(returned==total, "Actual new-key returns must equal primary plus victim keys; run() never accumulates addedT here.");
		check(field(BBDuk.class, "storedKmers").getLong(tool)==0, "Direct helper calls must not perform owner aggregation.");
		final String[] expected=spec.complete ? (spec.full+" "+(baseline ? spec.shortBefore : spec.shortAfter)).split(" ") : null;
		write(directory.resolve("keys.tsv"), "word\tkey\texpected_present\tshard_values\n");
		if(spec.complete){
			check(expected.length==(baseline ? spec.countBefore : spec.countAfter), "Literal complete-set cardinality must agree with the separately stated oracle.");
			for(int length=3; length<=4; length++){
				for(int value=0; value<(1<<(2*length)); value++){
					final String word=word(value, length); final boolean present=contains(expected, word);
					vector(directory, sets, word, present);
				}
			}
			check(total==expected.length, "Exhaustive length3/4 membership and complete primary/victim count must exclude any unlisted keys.");
		}else{
			vector(directory, sets, spec.target, baseline ? spec.before : spec.after);
			vector(directory, sets, spec.bases, true);// Primary full-length exact reference remains an unaffected control.
		}
		final boolean present=sets[(int)(key(spec.target)%7)].getValue(key(spec.target))==1;
		check(present==(baseline ? spec.before : spec.after), "Literal target membership must match this phase; target="+spec.target);
		if(present!=spec.after){discrepancies++;}
		append(out.resolve("cases.tsv"), spec.name+"\t"+(spec.complete ? "complete" : "targeted")+"\t"+spec.target+"\t"+present+"\t"+spec.after+"\t"+total+"\t"+(expected==null ? "NA" : expected.length)+"\t7NEW\tunchanged\n");
		if(spec.before==spec.after){controls++; append(out.resolve("controls.tsv"), spec.name+"\t"+spec.target+"\t"+present+"\t"+(expected==null ? "NA" : expected.length)+"\t7NEW\tunchanged\n");}
	}

	private void vector(final Path directory, final AbstractKmerTable[] sets, final String word, final boolean present) throws IOException{
		final long key=key(word); final int[] actual=new int[7], expected=new int[7]; Arrays.fill(expected, -1);
		if(present){expected[(int)(key%7)]=1;}
		for(int shard=0; shard<7; shard++){actual[shard]=sets[shard].getValue(key);}
		check(Arrays.equals(actual, expected), "Literal word/scaffold/shard mismatch: "+word+" key="+key+" actual="+Arrays.toString(actual)+" expected="+Arrays.toString(expected));
		append(directory.resolve("keys.tsv"), word+"\t"+key+"\t"+present+"\t"+Arrays.toString(actual)+"\n");
	}
	private void worker(final Thread worker, final long reads, final long bases, final long kmers) throws IllegalAccessException{
		check(worker.getState()==Thread.State.NEW && !worker.isAlive(), "Construct loaders only; never invoke start/run/fetch or queue paths.");
		check(refReads.getLong(worker)==reads && refBases.getLong(worker)==bases && refKmers.getLong(worker)==kmers && added.getLong(worker)==0, "Direct addToMap counts its one reference window but run() alone accumulates addedT.");
	}
	private static long tableSize(final AbstractKmerTable table){final HashArray array=(HashArray)table; return array.size()+array.victims().size();}
	/** Independent literal DNA packing, with a length bit; never call production mask/key helpers. */
	private static long key(final String word){long result=0; for(int i=0; i<word.length(); i++){final int base="ACGT".indexOf(word.charAt(i)); require(base>=0, "Oracle words must be literal DNA: "+word); result=(result<<2)|base;} return result|(1L<<(2*word.length()));}
	private static String word(int value, final int length){final char[] result=new char[length]; for(int i=length-1; i>=0; i--){result[i]="ACGT".charAt(value&3); value>>>=2;} return new String(result);}
	private static boolean contains(final String[] words, final String word){for(final String candidate : words){if(candidate.equals(word)){return true;}}return false;}
	private void same(final String relative) throws IOException{check(Arrays.equals(Files.readAllBytes(out.resolve(relative)), Files.readAllBytes(reference.resolve(relative))), "Frozen baseline literal/setup/control changed: "+relative);}
	private void saveDispatch() throws IOException{final StringBuilder text=new StringBuilder("scenario\tstatus\n"); for(int i=0; i<SPECS.length; i++){text.append(SPECS[i].name).append('\t').append(statuses[i]).append('\n');} write(out.resolve("dispatch.tsv"), text.toString());}
	private void finish() throws IOException{final String summary="BBDUK_SHORT_EDIT_MODE_PASS phase="+phase+" version=1 scenarios=7 owners="+owners+" loaders="+loaders+" calls="+calls+" completeCases=5 targetedCases=2 controls="+controls+" discrepancies="+discrepancies+" checks="+checks+" failures=0\n"; write(out.resolve("summary.txt"), summary); System.out.print(summary);}
	private void check(final boolean condition, final String why){checks++; require(condition, why);}
	private static Field field(final Class<?> type, final String name) throws NoSuchFieldException{final Field result=type.getDeclaredField(name); result.setAccessible(true); return result;}
	private static Path path(final String text){require(text!=null && !text.isEmpty() && text.indexOf('=')<0 && text.indexOf('\t')<0 && text.indexOf('\n')<0, "Use a nonempty path without equals, tabs or newlines for the legacy constructor."); return Paths.get(text).toAbsolutePath().normalize();}
	private static String repeat(final char value, final int length){final char[] values=new char[length]; Arrays.fill(values, value); return new String(values);}
	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static String text(final Path file) throws IOException{return new String(Files.readAllBytes(file), StandardCharsets.US_ASCII);}
	private static void write(final Path file, final String text) throws IOException{Files.write(file, bytes(text), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path file, final String text) throws IOException{Files.write(file, bytes(text), StandardOpenOption.APPEND);}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}

	private static final class Spec{
		Spec(final String name_, final String bases_, final int edist_, final int edist2_, final int hdist2_, final String full_, final String before_, final String after_, final int countBefore_, final int countAfter_, final boolean targetBefore_, final boolean targetAfter_){name=name_; bases=bases_; edist=edist_; edist2=edist2_; hdist2=hdist2_; full=full_; shortBefore=before_; shortAfter=after_; countBefore=countBefore_; countAfter=countAfter_; before=targetBefore_; after=targetAfter_; complete=full!=null; target=complete ? "AAC" : "ATAC";}
		final String name, bases, full, shortBefore, shortAfter, target;
		final int edist, edist2, hdist2, countBefore, countAfter;
		final boolean before, after, complete;
	}
	private final Path out, reference;
	private final String phase;
	private final boolean baseline;
	private final Constructor<?> constructor;
	private final Method add;
	private final Field refReads, refBases, refKmers, added;
	private final String[] statuses=new String[7];
	private int owners, loaders, calls, controls, discrepancies, checks;
	private static final String F1="ACGA";
	private static final String F19="AAGA ACAA ACCA ACGA ACGC ACGG ACGT ACTA AGGA ATGA CCGA GCGA TCGA AACG ACAG ACCG ACTG AGCG ATCG";
	private static final String S20="AAG ACA ACC ACG ACT AGA AGG ATG CAA CCA CCG CGA CGC CGG CGT CTA GCG GGA TCG TGA";
	private static final String S25=S20+" AAC AGC ATC CAG CTG";
	private static final Spec[] SPECS={
		new Spec("missed_short_indel", "ACGA", 0, 1, 0, F1, S20, S25, 21, 26, false, true),
		new Spec("unwanted_short_indel", "ACGA", 1, 0, 1, F19, S25, S20, 44, 39, true, false),
		new Spec("hamming_control", "ACGA", 0, 0, 1, F1, S20, S20, 21, 21, false, false),
		new Spec("matched_edit_control", "ACGA", 1, 1, 0, F19, S25, S25, 44, 44, true, true),
		new Spec("exact_control", "ACGA", 0, 0, 0, F1, "ACG CGA", "ACG CGA", 3, 3, false, false),
		new Spec("recursive_two_edits", "ACGTA", 0, 2, 0, null, null, null, -1, -1, false, true),
		new Spec("one_edit_absence_control", "ACGTA", 0, 1, 0, null, null, null, -1, -1, false, false)
	};
	private static final String[] SETTINGS={"ktrim=r", "maskmiddle=f", "rcomp=f", "speed=0", "rskip=1", "qskip=1", "qhdist=0", "qhdist2=0", "findbestmatch=f", "rename=f", "skipr1=f", "skipr2=f", "restrictleft=0", "restrictright=0", "minlen=1", "t=1", "simd=f", "qin=33", "qout=33", "array=t", "prealloc=f", "ow=f", "app=f"};
}
