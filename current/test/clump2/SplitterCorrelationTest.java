package test.clump2;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;

import clump2.Clump;
import stream.Read;
import structures.IntList;

/**
 * Finite oracle for Splitter correlated-partner shared-count reporting.
 * @author Sucrose, Shinobu
 */
public final class SplitterCorrelationTest{

	/** Runs baseline or fixed expectations for a non-final best partner. */
	public static void main(String[] args) throws Exception{
		check(args.length==1, "Expected baseline or fixed");
		final boolean fixed=args[0].equals("fixed");
		check(fixed || args[0].equals("baseline"), "Unknown expectation: "+args[0]);

		System.out.println("case\tbestVar2\treportedShared\texpectedShared");
		checkCase("non-final-best", fixed ? 3 : 1, VAR_A, VAR_B, VAR_C);
		checkCase("final-best", 3, VAR_A, VAR_B, VAR_NONE);
	}

	/** Checks a partner-count case with optional non-final singleton noise. */
	private static void checkCase(final String name, final int expectedShared,
			final int varA, final int varB, final int varC) throws Exception{
		final LinkedHashMap<Integer, ArrayList<Read>> map=makeMap(varA, varB, varC);
		final Clump clump=Clump.makeClump(3);
		addRead(clump, map, 0, varA, varB);
		addRead(clump, map, 1, varA, varB);
		if(varC>=0){
			addRead(clump, map, 2, varA, varB, varC);
			addRead(clump, map, 3, varC);
		}else{
			addRead(clump, map, 2, varA, varB);
		}
		//Reference-only reads make the supplied A/B/C variants non-consensus alleles.
		for(int id=4; id<8; id++){addRead(clump, map, id);}
		setMaxima(clump);

		final IntList collection=new IntList(8);
		final int[] vector=new int[5];
		final int bestVar2=examineVar(varA, map.get(varA), collection, vector, map);

		check(bestVar2==varB, "Expected var B as best partner, got "+bestVar2);
		check(vector[2]==varB, "Expected vector var2 B, got "+vector[2]);
		check(vector[3]==expectedShared, "Expected reported shared count "+
				expectedShared+", got "+vector[3]);

		System.out.println(name+'\t'+bestVar2+'\t'+vector[3]+'\t'+expectedShared);
	}

	/** Creates a variant-to-read map in the order the correlated scan would use. */
	private static LinkedHashMap<Integer, ArrayList<Read>> makeMap(final int... vars){
		final LinkedHashMap<Integer, ArrayList<Read>> map=
				new LinkedHashMap<Integer, ArrayList<Read>>();
		for(int var : vars){
			if(var>=0){map.put(var, new ArrayList<Read>());}
		}
		return map;
	}

	/** Adds a read containing the supplied variant list to every matching map entry. */
	private static void addRead(final Clump clump,
			final LinkedHashMap<Integer, ArrayList<Read>> map, final int id,
			final int... vars) throws Exception{
		final Read read=new Read(basesFor(vars).getBytes(StandardCharsets.US_ASCII),
				null, "splitter-"+id, id);
		makeKey(read);
		setVars(read, vars);
		clump.add(read);
		for(int var : vars){
			final ArrayList<Read> list=map.get(var);
			if(list!=null){list.add(read);}
		}
	}

	/** Creates and attaches the package-private ReadKey used by Splitter. */
	private static void makeKey(final Read read) throws Exception{
		final Class<?> readKey=Class.forName("clump2.ReadKey");
		final Method makeKey=readKey.getDeclaredMethod(
				"makeKey", Read.class, boolean.class);
		final Method set=readKey.getDeclaredMethod(
				"set", long.class, int.class, boolean.class);
		makeKey.setAccessible(true);
		set.setAccessible(true);
		set.invoke(makeKey.invoke(null, read, true), 3L, READ_POS, false);
	}

	/** Installs a literal variant list on the read's package-private key. */
	private static void setVars(final Read read, final int... vars) throws Exception{
		final IntList list=new IntList(vars.length);
		for(int var : vars){list.add(var);}
		final Field varsField=Class.forName("clump2.ReadKey").getDeclaredField("vars");
		varsField.setAccessible(true);
		varsField.set(read.obj(), list);
	}

	/** Sets the private clump-width cache used by Splitter's coordinate checks. */
	private static void setMaxima(final Clump clump) throws Exception{
		final Method setMaxima=Clump.class.getDeclaredMethod("setMaxima");
		setMaxima.setAccessible(true);
		setMaxima.invoke(clump);
	}

	/** Calls package-private Splitter.examineVar. */
	private static int examineVar(final int var, final ArrayList<Read> list,
			final IntList collection, final int[] vector,
			final LinkedHashMap<Integer, ArrayList<Read>> map) throws Exception{
		final Class<?> splitter=Class.forName("clump2.Splitter");
		final Method examineVar=splitter.getDeclaredMethod("examineVar",
				int.class, ArrayList.class, IntList.class, int[].class,
				LinkedHashMap.class);
		examineVar.setAccessible(true);
		return (Integer)examineVar.invoke(null, var, list, collection, vector, map);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Creates a read that contains the supplied variants at their clump columns. */
	private static String basesFor(final int... vars){
		final char[] bases=new char[READ_LENGTH];
		Arrays.fill(bases, 'T');
		for(int var : vars){
			bases[var>>2]=BASES[var&0x3];
		}
		return new String(bases);
	}

	private SplitterCorrelationTest(){}

	private static final int VAR_A=(10<<2)|0;
	private static final int VAR_B=(20<<2)|1;
	private static final int VAR_C=(30<<2)|2;
	private static final int VAR_NONE=-1;
	private static final int READ_LENGTH=40;
	private static final int READ_POS=0;
	private static final char[] BASES={'A', 'C', 'G', 'T'};
}
