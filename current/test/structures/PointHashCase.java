package test.structures;

import java.util.HashMap;
import java.util.HashSet;

import structures.Point;

/** Literal equals/hashCode contract oracle for finite Point coordinates.
 * @author Jean
 */
public final class PointHashCase{

	/** Runs baseline or fixed expectations for Point hashCode compatibility. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final boolean fixed="expect=fixed".equals(args[0]);
		check(fixed || "expect=baseline".equals(args[0]), "Unknown expectation: "+args[0]);
		final String expect=(fixed ? "fixed" : "baseline");

		final Pair finite=pairFor(1.5, -2.25, 1.5, -2.25);
		final Pair signedZero=pairFor(0.0, -0.0, -0.0, 0.0);

		checkObservation(expect, "finite_equal_object", "true",
				observeEquals(new Point(1.5, -2.25), new Point(1.5, -2.25)));
		checkObservation(expect, "finite_different_x", "false",
				observeEquals(new Point(1.5, -2.25), new Point(1.75, -2.25)));
		checkObservation(expect, "finite_different_y", "false",
				observeEquals(new Point(1.5, -2.25), new Point(1.5, -2.0)));
		checkObservation(expect, "signed_zero_equal", "true",
				observeEquals(signedZero.a, signedZero.b));
		checkObservation(expect, "object_subclass_equal", "false",
				observeEquals(new Point(1.5, -2.25), new ChildPoint(1.5, -2.25)));
		checkObservation(expect, "typed_subclass_equal", "true",
				observeTypedEquals(new Point(1.5, -2.25), new ChildPoint(1.5, -2.25)));

		checkObservation(expect, "finite_equal_hash", fixed ? "true" : "false",
				observeHashEqual(finite.a, finite.b));
		checkObservation(expect, "signed_zero_hash", fixed ? "true" : "false",
				observeHashEqual(signedZero.a, signedZero.b));

		checkObservation(expect, "hashset_finite_lookup", fixed ? "true" : "false",
				observeHashSetLookup(finite.a, finite.b));
		checkObservation(expect, "hashmap_finite_lookup", fixed ? "stored" : "null",
				observeHashMapLookup(finite.a, finite.b));
		checkObservation(expect, "hashset_signed_zero_lookup", fixed ? "true" : "false",
				observeHashSetLookup(signedZero.a, signedZero.b));
		checkObservation(expect, "hashmap_signed_zero_lookup", fixed ? "stored" : "null",
				observeHashMapLookup(signedZero.a, signedZero.b));
	}

	/** Calls Point.equals(Object) and returns its result or thrown type. */
	private static String observeEquals(final Point a, final Point b){
		try{
			return Boolean.toString(a.equals((Object)b));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}

	/** Calls Point.equals(Point) and returns its result or thrown type. */
	private static String observeTypedEquals(final Point a, final Point b){
		try{
			return Boolean.toString(a.equals(b));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}

	/** Compares Point hashCodes and returns the result or thrown type. */
	private static String observeHashEqual(final Point a, final Point b){
		try{
			return Boolean.toString(a.hashCode()==b.hashCode());
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}

	/** Inserts one Point into a HashSet and looks up an equal Point. */
	private static String observeHashSetLookup(final Point stored, final Point query){
		try{
			final HashSet<Point> set=new HashSet<Point>();
			set.add(stored);
			return Boolean.toString(set.contains(query));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}

	/** Inserts one Point into a HashMap and looks up an equal Point. */
	private static String observeHashMapLookup(final Point stored, final Point query){
		try{
			final HashMap<Point, String> map=new HashMap<Point, String>();
			map.put(stored, "stored");
			return String.valueOf(map.get(query));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}

	/**
	 * Returns equal points with distinct baseline identity hashes.
	 * Object identity hash collisions are possible, so baseline rows select a
	 * concrete non-colliding pair instead of assuming the first pair differs.
	 */
	private static Pair pairFor(final double aX, final double aY, final double bX,
			final double bY){
		for(int i=0; i<128; i++){
			final Pair pair=new Pair(new Point(aX, aY), new Point(bX, bY));
			check(pair.a.equals(pair.b), "Pair must be equal: "+pair.a+", "+pair.b);
			if(System.identityHashCode(pair.a)!=System.identityHashCode(pair.b)){return pair;}
		}
		throw new AssertionError("Inconclusive Point hashCode baseline: could not find "+
				"non-colliding identity hashes for equal Points in 128 attempts");
	}

	/** Checks and prints a value-or-exception outcome. */
	private static void checkObservation(final String expect, final String name,
			final String expected, final String observed){
		check(expected.equals(observed), name+" expected "+expected+", got "+observed);
		System.out.println("PASS\t"+expect+'\t'+name+'\t'+observed);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Pair of equal Points to compare through hash-based collections. */
	private static final class Pair{
		Pair(final Point a_, final Point b_){
			a=a_;
			b=b_;
		}

		final Point a;
		final Point b;
	}

	/** Same-coordinate subclass for Object.equals runtime-class checks. */
	private static final class ChildPoint extends Point{
		ChildPoint(final double x_, final double y_){
			super(x_, y_);
		}
	}

}
