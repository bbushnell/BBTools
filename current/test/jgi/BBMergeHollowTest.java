package test.jgi;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import assemble.KmerWindowDepth;
import assemble.Tadpole;
import assemble.Tadpole2;
import dna.AminoAcid;
import jgi.BBMerge;
import parse.PreParser;
import shared.Timer;
import stream.Read;
import stream.Streamer;
import stream.Writer;
import ukmer.AbstractKmerTableU;
import ukmer.Kmer;

/** Checks BBMerge's real constructor and paired extension path with controlled counts.
 * Reflection keeps test access out of the production API.
 * @author Brian Bushnell, Fischl */
public final class BBMergeHollowTest {

	/** Loads two small tables, then exercises both mates and repeated guarded attempts. */
	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, BBMergeHollowTest.class, false);
		String input=null;
		for(String arg : pp.args){
			final String[] split=arg.split("=", 2);
			if(split.length==2 && split[0].equalsIgnoreCase("in")){input=split[1];}
			else{throw new IllegalArgumentException("Expected in=<small fixture>, not "+arg);}
		}
		if(input==null){throw new IllegalArgumentException("A small input fixture is required.");}
		final BBMerge plain=make(input, "hollow=t"), guarded=make(input, "hollow=t", "hollowbranchstop=t");
		final Tadpole2 off=engine(plain), on=engine(guarded);
		check(bool(on, "hollow") && bool(on, "hollowBranchStop"), "Hollow/guard flags did not reach the embedded Tadpole.");
		check(!bool(off, "hollowBranchStop"), "Guard must remain off unless explicitly selected.");
		final Tadpole2 tuned=engine(make(input, "hollow=t", "prealloc=t", "hollowbranchratio=5", "hollowbranchsupport=7"));
		check(integer(tuned, "hollowBranchRatio")==5 && integer(tuned, "hollowBranchSupport")==7,
				"Threshold overrides did not reach the embedded Tadpole.");
		check(integer(tuned, "hashKmerMode")==1, "Hollow preallocation must retain resizable pair hashes.");
		check(!bool(engine(make(input, "hollow=f")), "hollow"), "Explicit solid mode unexpectedly enabled Hollow.");
		for(String[] bad : new String[][]{
			{"hollowbranchstop=t"}, {"hollow=t", "k=64"}, {"hollow=t", "hollowbranchstop=t", "k=66"},
			{"hollow=t", "extend2=0"}, {"hollow=t", "bridge=t"}, {"hollow=t", "ecct=t"},
			{"hollow=t", "eccb=t"}, {"hollow=t", "ecco=t"}, {"hollow=t", "kfilter=1"},
			{"hollow=t", "testmerge=t"}, {"hollow=t", "shave=t"}, {"hollow=t", "rinse=t"},
			{"hollow=t", "prefilter=1"}, {"hollow=t", "hashkmers=f"}, {"hollow=t", "hashkmers=fixed"},
			{"hollowbranchratio=0"}, {"hollowbranchsupport=0"}, {"hollow=t", "extend=20"}}){
			boolean rejected=false;
			try{make(input, bad);}catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "Unsupported combination accepted: "+Arrays.toString(bad));
		}
		off.loadKmers(new Timer());
		on.loadKmers(new Timer());
		final Class<?> workerClass=Class.forName("jgi.BBMerge$MateThread");
		final Constructor<?> constructor=workerClass.getDeclaredConstructor(BBMerge.class,
				Streamer.class, Writer.class, Writer.class, Writer.class, boolean.class, boolean.class);
		constructor.setAccessible(true);
		final Object offWorker=constructor.newInstance(plain, null, null, null, null, true, false);
		final Object onWorker=constructor.newInstance(guarded, null, null, null, null, true, false);
		final Method extend=workerClass.getDeclaredMethod("extendAndMerge", Read.class, Read.class,
				int.class, int.class, boolean.class, int[].class);
		extend.setAccessible(true);
		final Random random=new Random(42901);
		for(int span : new int[]{K+1, K-1}){
			for(int stop : new int[]{0, 4, 19}){
				for(int orientation=0; orientation<2; orientation++){
					final byte[] path=new byte[SEED+EXTENSION];
					for(int i=0; i<path.length; i++){path[i]=(byte)AminoAcid.numberToBase[random.nextInt(4)];}
					if(orientation==1){AminoAcid.reverseComplementBasesInPlace(path);}
					populate(off, path, span, stop);
					populate(on, path, span, stop);
					exercise(offWorker, extend, path, EXTENSION);
					exercise(onWorker, extend, path, stop);
				}
			}
		}
		System.out.println("BBMERGE_HOLLOW_TEST_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fixtures           ----------------*/
	/*--------------------------------------------------------------*/

	/** Uses normal CLI parsing; no production field is changed by the test. */
	private static BBMerge make(final String input, final String... options){
		assert(input!=null) : "BBMerge validates its actual input before the test accesses the embedded engine.";
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+input, "out=null", "t=1",
				"net=null", "k="+K, "extend2=20", "minprob=0", "prefilter=0", "prealloc=f"));
		args.addAll(Arrays.asList(options));
		return new BBMerge(args.toArray(new String[args.size()]));
	}

	/** Reads the actual constructor result, not a separately configured test engine. */
	private static Tadpole2 engine(final BBMerge merge) throws Exception{
		assert(merge!=null) : "A real BBMerge constructor must have completed before its engine is inspected.";
		final Field field=BBMerge.class.getDeclaredField("tadpole");
		field.setAccessible(true);
		return (Tadpole2)field.get(merge);
	}

	/** Reads a protected parser field without widening the production API. */
	private static Object value(final Tadpole engine, final String name) throws Exception{
		assert(engine!=null) : "Parser assertions require the same engine that BBMerge will call.";
		final Field field=Tadpole.class.getDeclaredField(name);
		field.setAccessible(true);
		return field.get(engine);
	}
	private static boolean bool(final Tadpole engine, final String name) throws Exception{return (Boolean)value(engine, name);}
	private static int integer(final Tadpole engine, final String name) throws Exception{return (Integer)value(engine, name);}

	/** Installs a unique nominal path plus one supported alternate-span ambiguity. */
	private static void populate(final Tadpole2 engine, final byte[] path, final int span, final int stop){
		assert(stop>=0 && stop<EXTENSION) : "The alternate branch must precede an otherwise available suffix base.";
		final AbstractKmerTableU[] tables=engine.tables.tables();
		final Kmer word=new Kmer(K, true), neighbor=new Kmer(K, true);
		for(int end=SEED; end<=path.length; end++){
			word.clear();
			for(int i=end-K; i<end; i++){word.addRight(path[i]);}
			put(tables, word, 30);
		}
		final int end=SEED+stop;
		KmerWindowDepth.setHollowQuery(word, path, end-span, span);
		put(tables, word, 10);
		final int chosen=AminoAcid.baseToNumber[path[end]];
		for(int base=0; base<4; base++){
			neighbor.setFrom(word); neighbor.addRightNumeric(base);
			put(tables, neighbor, base==chosen || base==((chosen+1)&3) ? 5 : 1);
		}
	}

	/** Exercises the real paired helper, mate orientation, qualities and per-call resets. */
	private static void exercise(final Object worker, final Method extend, final byte[] path, final int expected) throws Exception{
		assert(path.length==SEED+EXTENSION) : "Independent expected lengths must cover the complete requested extension.";
		final byte[] quality=new byte[SEED];
		Arrays.fill(quality, (byte)30);
		final Read first=new Read(Arrays.copyOf(path, SEED), quality.clone(), "first", 0);
		final Read second=new Read(Arrays.copyOf(path, SEED), quality.clone(), "second", 0);
		second.reverseComplementFast();
		final int[] counts=new int[2];
		extend.invoke(worker, first, second, EXTENSION, 1, false, counts);
		check(counts[0]==expected && counts[1]==expected, "Paired extension did not enforce the same stop on both mates.");
		check(Arrays.equals(first.bases, Arrays.copyOf(path, SEED+expected)), "First mate changed its original prefix or crossed the stop.");
		second.reverseComplementFast();
		check(Arrays.equals(second.bases, first.bases), "Second mate orientation or extension differs from the first.");
		check(first.quality.length==first.length() && second.quality.length==second.length()
				&& Arrays.equals(Arrays.copyOf(first.quality, SEED), quality)
				&& Arrays.equals(Arrays.copyOf(second.quality, SEED), quality), "Original qualities were changed or lengths diverged.");
		second.reverseComplementFast();
		if(expected<EXTENSION){
			extend.invoke(worker, first, second, 1, 1, false, counts);
			check(counts[0]==0 && counts[1]==0, "A repeated BBMerge call bypassed the guarded stop.");
		}
	}

	/** Stores fixture counts through normal pair-table insertion. */
	private static void put(final AbstractKmerTableU[] tables, final Kmer word, final int count){
		assert(count>0) : "A fixture key must not use the table's zero empty sentinel.";
		tables[word.mod(tables.length)].set(word, count);
	}

	/** Keeps comparisons fatal even when assertions are disabled accidentally. */
	private static void check(final boolean pass, final String message){
		if(!pass){throw new AssertionError(message);}
		checks++;
	}

	private static int checks;
	private static final int K=128, SEED=129, EXTENSION=20;
}
