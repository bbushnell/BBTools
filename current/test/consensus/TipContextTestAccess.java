package test.consensus;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import consensus.BaseNode;
import consensus.TipContextGuard;
import consensus.TipContextPanel;

/** Test-local access to existing package-private tip algorithms.
 * Keeps regression tests in test.consensus without widening production APIs.
 * Wrappers retain the actual production objects and arrays; no algorithm is
 * reimplemented here. Intended for deterministic tests, not benchmark timing.
 * @author Brian Bushnell, Sucrose */
final class TipContextTestAccess{

	/** Prevents instantiation of the access facade. */
	private TipContextTestAccess(){}

	/** Reads the real pseudocount-adjusted support calculation. */
	static int observations(BaseNode node, byte base){
		return (Integer)call(TipContextGuard.class, "observations", new Class<?>[]{BaseNode.class, byte.class}, node, base);
	}

	/** Searches the real candidate set, then refreshes the test's read-only field view. */
	static void scan(byte[] target, byte[] donor, int min, int max, CandidateSet set){
		call(TipContextGuard.class, "scan", new Class<?>[]{byte[].class, byte[].class, int.class, int.class, SET},
			target, donor, min, max, set.value);
		set.refresh();
	}

	/** Invokes the real context-origin check. */
	static boolean matchesContext(Candidate candidate, byte[] context, int offset, int length, boolean left){
		return (Boolean)call(TipContextGuard.class, "matchesContext",
			new Class<?>[]{CANDIDATE, byte[].class, int.class, int.class, boolean.class},
			candidate.value, context, offset, length, left);
	}

	/** Invokes the real evidence decision without copying candidate state. */
	static String decision(int original, int context, int pos, byte to, int length, CandidateSet candidates, boolean contextOK){
		return (String)call(TipContextGuard.class, "decision",
			new Class<?>[]{int.class, int.class, int.class, byte.class, int.class, SET, boolean.class},
			original, context, pos, to, length, candidates.value, contextOK);
	}

	/** Delegates orientation to the production reverse-complement helper. */
	static byte[] reverse(byte[] bases){
		return (byte[])call(TipContextGuard.class, "reverse", new Class<?>[]{byte[].class}, (Object)bases);
	}

	/** Builds a context with the production candidate and coordinate logic. */
	static byte[] context(byte[] target, Candidate donor, boolean left){
		return (byte[])call(TipContextPanel.class, "context", new Class<?>[]{byte[].class, CANDIDATE, boolean.class},
			target, donor.value, left);
	}

	/** Combines proposals with the production conservative merge. */
	static byte[] combine(byte[] original, byte[] left, byte[] right){
		return (byte[])call(TipContextPanel.class, "combine", new Class<?>[]{byte[].class, byte[].class, byte[].class},
			original, left, right);
	}

	/** Invokes a static helper and preserves its original runtime failure type. */
	private static Object call(Class<?> owner, String name, Class<?>[] parameters, Object... args){
		try{
			final Method method=owner.getDeclaredMethod(name, parameters);
			method.setAccessible(true);
			return method.invoke(null, args);
		}catch(InvocationTargetException e){throw failure(e.getCause());}
		catch(ReflectiveOperationException e){throw new AssertionError("Cannot access "+owner.getName()+"."+name, e);}
	}

	/** Constructs the actual production nested type rather than a parallel model. */
	private static Object create(Class<?> type, Class<?>[] parameters, Object... args){
		try{
			final Constructor<?> constructor=type.getDeclaredConstructor(parameters);
			constructor.setAccessible(true);
			return constructor.newInstance(args);
		}catch(InvocationTargetException e){throw failure(e.getCause());}
		catch(ReflectiveOperationException e){throw new AssertionError("Cannot construct "+type.getName(), e);}
	}

	/** Reads a test-visible field from its real production object. */
	private static Object field(Object value, String name){
		try{
			final Field field=value.getClass().getDeclaredField(name);
			field.setAccessible(true);
			return field.get(value);
		}catch(ReflectiveOperationException e){throw new AssertionError("Cannot read "+name, e);}
	}

	/** Returns the original runtime cause, including production assertions. */
	private static RuntimeException failure(Throwable cause){
		if(cause instanceof Error){throw (Error)cause;}
		if(cause instanceof RuntimeException){return (RuntimeException)cause;}
		throw new AssertionError("Unexpected checked failure from tip helper", cause);
	}

	/** Resolves only the named existing production nested classes. */
	private static Class<?> type(String suffix){
		try{return Class.forName("consensus.TipContextGuard$"+suffix);}
		catch(ClassNotFoundException e){throw new ExceptionInInitializerError(e);}
	}

	/** Read-only overlap view retaining the real candidate. */
	static final class Candidate{
		Candidate(byte[] bases, int overlap){this(create(CANDIDATE, new Class<?>[]{byte[].class, int.class}, bases, overlap));}
		private Candidate(Object value_){value=value_; overlap=(Integer)field(value, "overlap");}
		private final Object value;
		final int overlap;
	}

	/** Mutable production candidate set with refreshed observations for assertions. */
	static final class CandidateSet{
		CandidateSet(){value=create(SET, new Class<?>[0]); refresh();}
		private void refresh(){
			count=(Integer)field(value, "count");
			final Object candidate=field(value, "only");
			only=candidate==null ? null : new Candidate(candidate);
		}
		private final Object value;
		int count;
		Candidate only;
	}

	/** Shares the real depth array so missing/zero/veto cases exercise production storage. */
	static final class KmerEvidence{
		KmerEvidence(int length){
			value=create(KMERS, new Class<?>[]{int.class}, length);
			depths=(int[])field(value, "depths");
		}
		int support(int pos){
			try{
				final Method method=KMERS.getDeclaredMethod("support", int.class);
				method.setAccessible(true);
				return (Integer)method.invoke(value, pos);
			}catch(InvocationTargetException e){throw failure(e.getCause());}
			catch(ReflectiveOperationException e){throw new AssertionError("Cannot read real kmer support", e);}
		}
		private final Object value;
		final int[] depths;
	}

	private static final Class<?> CANDIDATE=type("Candidate");
	private static final Class<?> SET=type("CandidateSet");
	private static final Class<?> KMERS=type("KmerEvidence");
}
