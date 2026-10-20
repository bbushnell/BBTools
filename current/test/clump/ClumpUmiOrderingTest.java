package test.clump;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import clump.Clump;
import clump.KmerComparator;
import hiseq.FlowcellCoordinate;
import stream.Read;

/**
 * Finite removeDuplicates oracle for UMI mismatches that precede later matches.
 * @author Sucrose
 */
public final class ClumpUmiOrderingTest{

	/** Runs baseline or fixed expectations for separated UMI matches. */
	public static void main(String[] args) throws Exception{
		check(args.length==1, "Expected baseline or fixed");
		final boolean fixed=args[0].equals("fixed");
		check(fixed || args[0].equals("baseline"), "Unknown expectation: "+args[0]);

		final StaticState state=new StaticState();
		try{
			KmerComparator.verbose=false;
			KmerComparator.compareSequence=true;
			System.out.println("case\tcompareUmi\tremoved\tretained\torder");
			checkOrdinarySeparated(fixed);
			checkOrdinaryAdjacent();
			checkOrdinaryCompareUmiOff();
			checkOrdinaryScanLimitBound();
			checkContainmentSeparated(fixed);
			checkContainmentAdjacent();
			checkContainmentCompareUmiOff();
		}finally{
			state.restore();
		}
	}

	/** Same-UMI ordinary duplicate separated by a different UMI after real sort. */
	private static void checkOrdinarySeparated(final boolean fixed) throws Exception{
		configure(false, true);
		final Clump clump=makeSorted(
				read("A_HIGH", "ACGT", UMI_A, 2, false, 40, 101),
				read("B_MID", "ACGT", UMI_B, 2, false, 30, 102),
				read("A_LOW", "ACGT", UMI_A, 2, false, 20, 103));
		checkSeparatedUmis(clump);
		checkCase("ordinary-separated", false, clump, "A_HIGH,B_MID,A_LOW",
				fixed ? 1 : 0, fixed ? "A_HIGH,B_MID" : "A_HIGH,B_MID,A_LOW");
	}

	/** Adjacent same-UMI ordinary duplicate, proving duplicate removal still works. */
	private static void checkOrdinaryAdjacent() throws Exception{
		configure(false, true);
		final Clump clump=makeSorted(
				read("A_HIGH", "ACGT", UMI_A, 2, false, 40, 101),
				read("A_MID", "ACGT", UMI_A, 2, false, 30, 102),
				read("B_LOW", "ACGT", UMI_B, 2, false, 20, 103));
		checkOrdinaryAdjacentUmis(clump);
		checkCase("ordinary-adjacent", false, clump, "A_HIGH,A_MID,B_LOW",
				1, "A_HIGH,B_LOW");
	}

	/** UMI filtering disabled: the separated ordinary clump collapses normally. */
	private static void checkOrdinaryCompareUmiOff() throws Exception{
		configure(false, true);
		final Clump clump=makeSorted(
				read("A_HIGH", "ACGT", UMI_A, 2, false, 40, 101),
				read("B_MID", "ACGT", UMI_B, 2, false, 30, 102),
				read("A_LOW", "ACGT", UMI_A, 2, false, 20, 103));
		checkSeparatedUmis(clump);
		checkCase("ordinary-compareumi-off", false, clump, "A_HIGH,B_MID,A_LOW",
				2, "A_HIGH");
	}

	/** UMI mismatches remain bounded by scanLimit when treated as unequals. */
	private static void checkOrdinaryScanLimitBound() throws Exception{
		configure(false, true);
		final Clump clump=makeSorted(
				read("A_HIGH", "ACGT", UMI_A, 2, false, 50, 101),
				read("C_NEXT", "ACGT", UMI_C, 2, false, 40, 102),
				read("B_MID", "ACGT", UMI_B, 2, false, 30, 103),
				read("A_LOW", "ACGT", UMI_A, 2, false, 20, 104));
		checkBoundedUmis(clump);
		checkCase("ordinary-scanlimit-bound", false, clump,
				"A_HIGH,C_NEXT,B_MID,A_LOW", 0, "A_HIGH,C_NEXT,B_MID,A_LOW");
	}

	/** Same-UMI backward containment separated by a different UMI after real sort. */
	private static void checkContainmentSeparated(final boolean fixed) throws Exception{
		configure(true, true);
		final Clump clump=makeSorted(
				read("A_LEFT", "ACGT", UMI_A, 2, false, 40, 201),
				read("B_RIGHT", "CGTA", UMI_B, 1, false, 30, 202),
				read("A_LONG", "TACGTA", UMI_A, 2, true, 20, 203));
		checkSeparatedUmis(clump);
		checkContainmentGeometry(clump, "containment-separated");
		checkCase("containment-separated", true, clump, "A_LEFT,B_RIGHT,A_LONG",
				fixed ? 1 : 0, fixed ? "B_RIGHT,A_LONG" : "A_LEFT,B_RIGHT,A_LONG");
	}

	/** Adjacent same-UMI backward containment, proving the fixed pass still works. */
	private static void checkContainmentAdjacent() throws Exception{
		configure(true, true);
		final Clump clump=makeSorted(
				read("B_LEFT", "ACGT", UMI_B, 2, false, 40, 201),
				read("A_RIGHT", "CGTA", UMI_A, 1, false, 30, 202),
				read("A_LONG", "TACGTA", UMI_A, 2, true, 20, 203));
		checkContainmentAdjacentUmis(clump);
		checkContainmentGeometry(clump, "containment-adjacent");
		checkCase("containment-adjacent", true, clump, "B_LEFT,A_RIGHT,A_LONG",
				1, "B_LEFT,A_LONG");
	}

	/** UMI filtering disabled: the backward container removes both contained reads. */
	private static void checkContainmentCompareUmiOff() throws Exception{
		configure(true, true);
		final Clump clump=makeSorted(
				read("A_LEFT", "ACGT", UMI_A, 2, false, 40, 201),
				read("B_RIGHT", "CGTA", UMI_B, 1, false, 30, 202),
				read("A_LONG", "TACGTA", UMI_A, 2, true, 20, 203));
		checkSeparatedUmis(clump);
		checkContainmentGeometry(clump, "containment-compareumi-off");
		checkCase("containment-compareumi-off", true, clump,
				"A_LEFT,B_RIGHT,A_LONG", 2, "A_LONG");
	}

	/** Exercises the public removeDuplicates entry point for a sorted clump. */
	private static void checkCase(final String name, final boolean containment,
			final Clump clump, final String expectedOrder,
			final int expectedRemoved, final String expectedAfter) throws Exception{
		configure(containment, !name.endsWith("off"));
		check(labels(clump).equals(expectedOrder),
				name+" expected comparator order "+expectedOrder+", got "+labels(clump));

		final int removed=clump.removeDuplicates();
		check(removed==expectedRemoved,
				name+" expected "+expectedRemoved+" removed, got "+removed);
		check(labels(clump).equals(expectedAfter),
				name+" expected retained reads "+expectedAfter+", got "+labels(clump));
		System.out.println(name+'\t'+!name.endsWith("off")+'\t'+removed+'\t'
				+labels(clump)+'\t'+expectedOrder);
	}

	/** Creates, sorts, and validates a Clump through the real KmerComparator. */
	private static Clump makeSorted(final Read... reads) throws Exception{
		final Clump clump=Clump.makeClump(KMER);
		for(Read read : reads){clump.add(read);}
		clump.sort(new KmerComparator(K, false, false));
		return clump;
	}

	/** Validates the separated A/B/A UMI pattern after comparator sorting. */
	private static void checkSeparatedUmis(final Clump clump) throws Exception{
		check(umiMatches(clump.get(0), clump.get(2)),
				"outer reads must share UMI after sorting");
		check(!umiMatches(clump.get(0), clump.get(1)),
				"intervening read must have a different UMI");
	}

	/** Validates the A/A/B UMI pattern after ordinary comparator sorting. */
	private static void checkOrdinaryAdjacentUmis(final Clump clump) throws Exception{
		check(umiMatches(clump.get(0), clump.get(1)),
				"adjacent reads must share UMI after sorting");
		check(!umiMatches(clump.get(0), clump.get(2)),
				"control read must have a different UMI");
	}

	/** Validates the B/A/A UMI pattern after containment comparator sorting. */
	private static void checkContainmentAdjacentUmis(final Clump clump) throws Exception{
		check(umiMatches(clump.get(1), clump.get(2)),
				"adjacent reads must share UMI after sorting");
		check(!umiMatches(clump.get(0), clump.get(2)),
				"control read must have a different UMI");
	}

	/** Validates the A/B/C/A UMI pattern for the scanLimit control. */
	private static void checkBoundedUmis(final Clump clump) throws Exception{
		check(umiMatches(clump.get(0), clump.get(3)),
				"outer reads must share UMI after sorting");
		check(!umiMatches(clump.get(0), clump.get(1)),
				"first intervening read must have a different UMI");
		check(!umiMatches(clump.get(0), clump.get(2)),
				"second intervening read must have a different UMI");
	}

	/** Verifies backward-containment geometry before duplicate removal mutates reads. */
	private static void checkContainmentGeometry(final Clump clump,
			final String name){
		final Read left=clump.get(0);
		final Read right=clump.get(1);
		final Read container=clump.get(2);
		check(Clump.contains(container, left, 0),
				name+" long read should contain left read");
		check(Clump.contains(container, right, 0),
				name+" long read should contain right read");
		check(!Clump.contains(left, container, 0),
				name+" left read should not contain long read");
		check(!Clump.contains(right, container, 0),
				name+" right read should not contain long read");
	}

	/** Sets duplicate-detection flags to the UMI-ordering test shape. */
	private static void configure(final boolean containment,
			final boolean compareUmi) throws Exception{
		set(Clump.class, "allowNs", false);
		set(Clump.class, "markAll", false);
		set(Clump.class, "markOnly", false);
		set(Clump.class, "opticalOnly", false);
		set(Clump.class, "containment", containment);
		set(Clump.class, "compareUMI", compareUmi);
		set(Clump.class, "umiSubs", 0);
		set(Clump.class, "affix", false);
		set(Clump.class, "printDuplicates", false);
		set(Clump.class, "renameByCount", false);
		set(Clump.class, "scanLimit", 1);
		set(Clump.class, "maxSubstitutions", 0);
		set(Clump.class, "dupeSubRate", 0.1f);
		set(Clump.class, "maxOpticalDistance", 40f);
		set(Clump.class, "forceSortXY", false);
		set(Clump.class, "k", K);
		FlowcellCoordinate.parseUMI=compareUmi;
	}

	/** Creates a read with an attached package-private ReadKey. */
	private static Read read(final String label, final String bases,
			final String umi, final int pos, final boolean minus,
			final int quality, final long numericID) throws Exception{
		final byte[] baseArray=bases.getBytes(StandardCharsets.US_ASCII);
		final byte[] qualityArray=new byte[baseArray.length];
		Arrays.fill(qualityArray, (byte)quality);
		final Read read=new Read(baseArray, qualityArray, id(umi, label), numericID);
		makeKey(read, pos, minus);
		return read;
	}

	/** Builds the extended Illumina coordinate field parsed by index3(). */
	private static String id(final String umi, final String label){
		return "INST:1:FLOW:1:1101:100:100:"+umi+" 1:N:0:"+label;
	}

	/** Creates and attaches the package-private ReadKey used by Clump. */
	private static void makeKey(final Read read, final int pos,
			final boolean minus) throws Exception{
		final Class<?> readKey=Class.forName("clump.ReadKey");
		final Method makeKey=readKey.getDeclaredMethod(
				"makeKey", Read.class, boolean.class);
		final Method set=readKey.getDeclaredMethod(
				"set", long.class, int.class, boolean.class);
		makeKey.setAccessible(true);
		set.setAccessible(true);
		final Object key=makeKey.invoke(null, read, true);
		set.invoke(key, KMER, pos, minus);
		checkPivot(read, pos, minus);
	}

	/** Tests that the synthetic key is attached to a valid canonical G pivot. */
	private static void checkPivot(final Read read, final int pos,
			final boolean minus){
		check(pos>=0 && pos<read.length(),
				read.id+" pivot position "+pos+" outside length "+read.length());
		final byte expected=(byte)(minus ? 'C' : 'G');
		check(read.bases[pos]==expected,
				read.id+" expected "+(char)expected+" at pivot position "+pos+
				", got "+(char)read.bases[pos]);
	}

	/** Tests the attached package-private ReadKey UMI state. */
	private static boolean umiMatches(final Read a, final Read b) throws Exception{
		final Class<?> readKey=Class.forName("clump.ReadKey");
		final Method umiMatches=readKey.getDeclaredMethod(
				"umiMatches", readKey, int.class);
		umiMatches.setAccessible(true);
		return (Boolean)umiMatches.invoke(a.obj, b.obj, 0);
	}

	/** Returns the sorted test labels in clump order. */
	private static String labels(final Clump clump){
		final StringBuilder sb=new StringBuilder();
		for(Read read : clump){
			if(sb.length()>0){sb.append(',');}
			sb.append(read.id.substring(read.id.lastIndexOf(':')+1));
		}
		return sb.toString();
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
			verbose=KmerComparator.verbose;
			compareSequence=KmerComparator.compareSequence;
		}

		void restore() throws Exception{
			for(int i=0; i<FIELDS.length; i++){FIELDS[i].set(null, values[i]);}
			FlowcellCoordinate.parseUMI=parseUMI;
			KmerComparator.verbose=verbose;
			KmerComparator.compareSequence=compareSequence;
		}

		private final Object[] values;
		private final boolean parseUMI;
		private final boolean verbose;
		private final boolean compareSequence;
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
	private static final int K=1;
	private static final String UMI_A="AAAA";
	private static final String UMI_B="TTTT";
	private static final String UMI_C="CCCC";
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
			clumpField("forceSortXY"),
			clumpField("k")
	};

	private ClumpUmiOrderingTest(){}
}
