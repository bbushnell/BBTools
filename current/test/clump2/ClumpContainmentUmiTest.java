package test.clump2;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import clump2.Clump;
import hiseq.FlowcellCoordinate;
import stream.Read;

/**
 * Finite removeDuplicates oracle for backward containment with UMI filtering.
 * @author Sucrose
 */
public final class ClumpContainmentUmiTest{

	/** Runs baseline or fixed expectations for backward-only UMI containment. */
	public static void main(String[] args) throws Exception{
		check(args.length==1, "Expected baseline or fixed");
		final boolean fixed=args[0].equals("fixed");
		check(fixed || args[0].equals("baseline"), "Unknown expectation: "+args[0]);

		final StaticState state=new StaticState();
		try{
			System.out.println("case\tcompareUmi\tmatchingUmi\tremoved\tafter");
			checkCase("different-umi", true, false, fixed ? 0 : 1, fixed ? 2 : 1);
			checkCase("matching-umi", true, true, 1, 1);
			checkCase("compareumi-off", false, false, 1, 1);
		}finally{
			state.restore();
		}
	}

	/** Exercises the real removeDuplicates entry point for one two-read clump. */
	private static void checkCase(final String name, final boolean compareUmi,
			final boolean matchingUmi, final int expectedRemoved,
			final int expectedSize) throws Exception{
		configure(compareUmi);
		final Clump clump=Clump.makeClump(KMER);
		final Read shortRead=read("ACGT", "AAAA", SHORT_POS, 1);
		final Read longRead=read("TACGTA", matchingUmi ? "AAAA" : "TTTT", LONG_POS, 2);
		check(umiMatches(shortRead, longRead)==matchingUmi,
				name+" constructed ReadKeys must reflect the requested UMI state");
		check(!Clump.equals(shortRead, longRead, 0, 0),
				name+" must not match in the forward short-vs-long orientation");
		check(Clump.equals(longRead, shortRead, 0, 0),
				name+" must match only when the backward scan compares long-vs-short");
		clump.add(shortRead);
		clump.add(longRead);

		final int removed=clump.removeDuplicates();
		check(removed==expectedRemoved,
				name+" expected "+expectedRemoved+" removed, got "+removed);
		check(clump.size()==expectedSize,
				name+" expected size "+expectedSize+", got "+clump.size());
		System.out.println(name+'\t'+compareUmi+'\t'+matchingUmi+'\t'+removed+'\t'+clump.size());
	}

	/** Sets duplicate-detection flags to the backward-containment case. */
	private static void configure(final boolean compareUmi) throws Exception{
		set(Clump.class, "allowNs", false);
		set(Clump.class, "markAll", false);
		set(Clump.class, "markOnly", false);
		set(Clump.class, "opticalOnly", false);
		set(Clump.class, "containment", true);
		set(Clump.class, "compareUMI", compareUmi);
		set(Clump.class, "umiSubs", 0);
		set(Clump.class, "affix", false);
		set(Clump.class, "printDuplicates", false);
		set(Clump.class, "renameByCount", false);
		set(Clump.class, "scanLimit", 5);
		set(Clump.class, "maxSubstitutions", 0);
		set(Clump.class, "dupeSubRate", 0f);
		set(Clump.class, "maxOpticalDistance", 40f);
		set(Clump.class, "forceSortXY", false);
		FlowcellCoordinate.parseUMI=compareUmi;
	}

	/** Creates a read with an attached package-private ReadKey. */
	private static Read read(final String bases, final String umi, final int pos,
			final long numericID) throws Exception{
		final byte[] baseArray=bases.getBytes(StandardCharsets.US_ASCII);
		final byte[] qualityArray=new byte[baseArray.length];
		Arrays.fill(qualityArray, (byte)40);
		final Read read=new Read(baseArray, qualityArray, id(umi), numericID);
		makeKey(read, pos);
		return read;
	}

	/** Builds the extended Illumina coordinate field parsed by index3(). */
	private static String id(final String umi){
		return "INST:1:FLOW:1:1101:100:100:"+umi+" 1:N:0:ACGT";
	}

	/** Creates and attaches the package-private ReadKey used by Clump. */
	private static void makeKey(final Read read, final int pos) throws Exception{
		final Class<?> readKey=Class.forName("clump2.ReadKey");
		final Method makeKey=readKey.getDeclaredMethod(
				"makeKey", Read.class, boolean.class);
		final Method set=readKey.getDeclaredMethod(
				"set", long.class, int.class, boolean.class);
		makeKey.setAccessible(true);
		set.setAccessible(true);
		set.invoke(makeKey.invoke(null, read, true), KMER, pos, false);
	}

	/** Tests the attached package-private ReadKey UMI state. */
	private static boolean umiMatches(final Read a, final Read b) throws Exception{
		final Class<?> readKey=Class.forName("clump2.ReadKey");
		final Method umiMatches=readKey.getDeclaredMethod(
				"umiMatches", readKey, int.class);
		umiMatches.setAccessible(true);
		return (Boolean)umiMatches.invoke(a.obj, b.obj, 0);
	}

	/** Writes a reflected static field. */
	private static void set(final Class<?> owner, final String name,
			final Object value) throws Exception{
		final Field field=owner.getDeclaredField(name);
		field.setAccessible(true);
		field.set(null, value);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final class StaticState{

		StaticState() throws Exception{
			values=new Object[FIELDS.length];
			for(int i=0; i<FIELDS.length; i++){values[i]=FIELDS[i].get(null);}
			parseUMI=FlowcellCoordinate.parseUMI;
		}

		void restore() throws Exception{
			for(int i=0; i<FIELDS.length; i++){FIELDS[i].set(null, values[i]);}
			FlowcellCoordinate.parseUMI=parseUMI;
		}

		private final Object[] values;
		private final boolean parseUMI;
	}

	/** Resolves and makes a Clump static field writable. */
	private static Field clumpField(final String name){
		try{
			final Field field=Clump.class.getDeclaredField(name);
			field.setAccessible(true);
			return field;
		}catch(ReflectiveOperationException e){
			throw new ExceptionInInitializerError(e);
		}
	}

	private static final long KMER=2L;
	private static final int SHORT_POS=2;
	private static final int LONG_POS=3;
	private static final Field[] FIELDS={
			clumpField("allowNs"),
			clumpField("markAll"),
			clumpField("markOnly"),
			clumpField("opticalOnly"),
			clumpField("containment"),
			clumpField("compareUMI"),
			clumpField("umiSubs"),
			clumpField("affix"),
			clumpField("printDuplicates"),
			clumpField("renameByCount"),
			clumpField("scanLimit"),
			clumpField("maxSubstitutions"),
			clumpField("dupeSubRate"),
			clumpField("maxOpticalDistance"),
			clumpField("forceSortXY")
	};

	private ClumpContainmentUmiTest(){}
}
