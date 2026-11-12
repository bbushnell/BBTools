package test.align2;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Checks the fraction-exclusion setters used by legacy align2 indexes.
 * Fixed mode verifies that invalid fractions are rejected before the shared static
 * thresholds change. Baseline mode records the old contract: negative, NaN and
 * infinity assert under -ea, while a large finite fraction returns and stores
 * negative histogram indexes. No reference, index, mapper, or disk cache is built.
 * @author Brian Bushnell, Jean
 */
public final class FractionExclusionGuardTest{

	/** Accepts mode=fixed (default) or mode=baseline, plus launcher JVM flags. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, FractionExclusionGuardTest.class, false);
		try{
			boolean baseline=false;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				if(pair.length!=2 || !pair[0].equalsIgnoreCase("mode") ||
						(!pair[1].equalsIgnoreCase("baseline") && !pair[1].equalsIgnoreCase("fixed"))){
					throw new IllegalArgumentException("Expected mode=fixed or mode=baseline: "+arg);
				}
				baseline=pair[1].equalsIgnoreCase("baseline");
			}
			boolean assertions=false;
			assert(assertions=true);
			if(!assertions){throw new IllegalStateException("Fraction baseline observations require -ea, as supplied by the launcher.");}
			if(baseline){runBaseline();}else{runFixed();}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Static fixture entry point only. */
	private FractionExclusionGuardTest(){}

	/** Runs the post-fix contract over accepted and rejected fraction cases. */
	private static void runFixed() throws Exception{
		final ByteBuilder out=new ByteBuilder(4096);
		out.append("class\tcase\tpassed\n");
		int checks=0;
		for(Target target : TARGETS){
			final Snapshot original=target.read();
			try{
				checks+=expectAccepted(out, target, "typical_0.03", 0.03f, new Snapshot(0.03f, 895, 931, 958, 970, 916));
				checks+=expectAccepted(out, target, "zero", 0f, new Snapshot(0f, 1000, 1000, 1000, 1000, 1000));
				checks+=expectAccepted(out, target, "safe_two_sevenths", 2f/7f, new Snapshot(2f/7f, 0, 342, 599, 714, 199));
				checks+=expectRejected(out, target, "negative", -0.1f);
				checks+=expectRejected(out, target, "nan", Float.NaN);
				checks+=expectRejected(out, target, "positive_infinity", Float.POSITIVE_INFINITY);
				checks+=expectRejected(out, target, "negative_percentile", 0.5f);
			}finally{target.write(original);}
		}
		System.out.print(out);
		System.out.println("FRACTION_EXCLUSION_GUARD_FIXED checks="+checks+" failed=0");
	}

	/** Runs the original contract controls and the known large-fraction defect. */
	private static void runBaseline() throws Exception{
		final ByteBuilder out=new ByteBuilder(4096);
		out.append("class\tcase\tpassed\n");
		int checks=0;
		for(Target target : TARGETS){
			final Snapshot original=target.read();
			try{
				checks+=expectAccepted(out, target, "typical_0.03", 0.03f, new Snapshot(0.03f, 895, 931, 958, 970, 916));
				checks+=expectAccepted(out, target, "zero", 0f, new Snapshot(0f, 1000, 1000, 1000, 1000, 1000));
				checks+=expectAccepted(out, target, "safe_two_sevenths", 2f/7f, new Snapshot(2f/7f, 0, 342, 599, 714, 199));
				checks+=expectAssertion(out, target, "negative", -0.1f);
				checks+=expectAssertion(out, target, "nan", Float.NaN);
				checks+=expectAssertion(out, target, "positive_infinity", Float.POSITIVE_INFINITY);
				checks+=expectAccepted(out, target, "negative_percentile", 0.5f, new Snapshot(0.5f, -750, -149, 300, 500, -399));
			}finally{target.write(original);}
		}
		System.out.print(out);
		System.out.println("FRACTION_EXCLUSION_GUARD_BASELINE checks="+checks+" failed=0");
	}

	/** Verifies an accepted setter call and its exact shared-threshold snapshot. */
	private static int expectAccepted(final ByteBuilder out, final Target target, final String name,
			final float value, final Snapshot expected) throws ReflectiveOperationException{
		final Throwable thrown=target.call(value);
		if(thrown!=null){throw new AssertionError(target.name+" "+name+" threw "+thrown, thrown);}
		final Snapshot observed=target.read();
		if(!observed.equals(expected)){
			throw new AssertionError(target.name+" "+name+" expected "+expected+" but observed "+observed);
		}
		return pass(out, target, name);
	}

	/** Verifies fixed-code rejection leaves the shared thresholds unchanged. */
	private static int expectRejected(final ByteBuilder out, final Target target, final String name,
			final float value) throws ReflectiveOperationException{
		final Snapshot before=target.read();
		final Throwable thrown=target.call(value);
		if(!(thrown instanceof IllegalArgumentException)){
			throw new AssertionError(target.name+" "+name+" expected IllegalArgumentException but observed "+thrown);
		}
		final Snapshot after=target.read();
		if(!after.equals(before)){
			throw new AssertionError(target.name+" "+name+" changed shared state from "+before+" to "+after);
		}
		return pass(out, target, name);
	}

	/** Verifies baseline assertion rejection leaves shared thresholds unchanged. */
	private static int expectAssertion(final ByteBuilder out, final Target target, final String name,
			final float value) throws ReflectiveOperationException{
		final Snapshot before=target.read();
		final Throwable thrown=target.call(value);
		if(!(thrown instanceof AssertionError)){
			throw new AssertionError(target.name+" "+name+" expected AssertionError but observed "+thrown);
		}
		final Snapshot after=target.read();
		if(!after.equals(before)){
			throw new AssertionError(target.name+" "+name+" changed shared state from "+before+" to "+after);
		}
		return pass(out, target, name);
	}

	/** Records one successful check row. */
	private static int pass(final ByteBuilder out, final Target target, final String name){
		out.append(target.name).tab().append(name).tab().append("true\n");
		return 1;
	}

	/** Reflective handle for one legacy index class and its static thresholds. */
	private static final class Target{
		Target(final String name_) throws ReflectiveOperationException{
			name=name_;
			clazz=Class.forName("align2."+name);
			setter=clazz.getDeclaredMethod("setFractionToExclude", float.class);
			fraction=field("FRACTION_GENOME_TO_EXCLUDE");
			drop=field("MIN_INDEX_TO_DROP_LONG_HIT_LIST");
			average=field("MAX_AVERAGE_LIST_TO_SEARCH");
			average2=field("MAX_AVERAGE_LIST_TO_SEARCH2");
			single=field("MAX_SINGLE_LIST_TO_SEARCH");
			shortest=field("MAX_SHORTEST_LIST_TO_SEARCH");
		}

		/** Calls the public setter and unwraps its user-visible cause. */
		Throwable call(final float value) throws ReflectiveOperationException{
			try{
				setter.invoke(null, value);
				return null;
			}catch(InvocationTargetException ite){
				return ite.getCause();
			}
		}

		/** Reads the mutable static fields touched by the setter. */
		Snapshot read() throws ReflectiveOperationException{
			return new Snapshot(fraction.getFloat(null), drop.getInt(null), average.getInt(null),
					average2.getInt(null), single.getInt(null), shortest.getInt(null));
		}

		/** Restores the mutable static fields after one target class finishes. */
		void write(final Snapshot s) throws ReflectiveOperationException{
			fraction.setFloat(null, s.fraction);
			drop.setInt(null, s.drop);
			average.setInt(null, s.average);
			average2.setInt(null, s.average2);
			single.setInt(null, s.single);
			shortest.setInt(null, s.shortest);
		}

		private Field field(final String fieldName) throws ReflectiveOperationException{
			final Field field=clazz.getDeclaredField(fieldName);
			field.setAccessible(true);
			return field;
		}

		final String name;
		final Class<?> clazz;
		final Method setter;
		final Field fraction;
		final Field drop;
		final Field average;
		final Field average2;
		final Field single;
		final Field shortest;
	}

	/** Immutable copy of the setter-controlled fraction and thresholds. */
	private static final class Snapshot{
		Snapshot(final float fraction_, final int drop_, final int average_, final int average2_,
				final int single_, final int shortest_){
			fraction=fraction_;
			drop=drop_;
			average=average_;
			average2=average2_;
			single=single_;
			shortest=shortest_;
		}

		@Override
		public boolean equals(final Object obj){
			if(!(obj instanceof Snapshot)){return false;}
			final Snapshot s=(Snapshot)obj;
			return Float.floatToIntBits(fraction)==Float.floatToIntBits(s.fraction) && drop==s.drop &&
					average==s.average && average2==s.average2 && single==s.single && shortest==s.shortest;
		}

		@Override
		public int hashCode(){
			int hash=Float.floatToIntBits(fraction);
			hash=31*hash+drop;
			hash=31*hash+average;
			hash=31*hash+average2;
			hash=31*hash+single;
			hash=31*hash+shortest;
			return hash;
		}

		@Override
		public String toString(){
			return "Snapshot(f="+fraction+", drop="+drop+", average="+average+
					", average2="+average2+", single="+single+", shortest="+shortest+")";
		}

		final float fraction;
		final int drop;
		final int average;
		final int average2;
		final int single;
		final int shortest;
	}

	private static final Target[] TARGETS;
	static{
		try{
			TARGETS=new Target[] {
				new Target("BBIndex5"),
				new Target("BBIndexAcc"),
				new Target("BBIndexPacBio"),
				new Target("BBIndexPacBioSkimmer")
			};
		}catch(ReflectiveOperationException roe){
			throw new ExceptionInInitializerError(roe);
		}
	}
}
