package test.align2;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Checks one mapper constructor fraction option in an isolated JVM.
 * The outer recipe runs one process per mapper/case pair; this fixture does not
 * load an index, open input/output streams, or spawn subprocesses.
 * @author Brian Bushnell, Jean
 */
public final class MapperFractionOptionTest{

	/** Accepts explicit mode=baseline|fixed, mapper=..., and case=... arguments. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, MapperFractionOptionTest.class, false);
		try{
			String mode=null, mapper=null, name=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				if(pair.length!=2){throw new IllegalArgumentException("Expected key=value argument: "+arg);}
				if(pair[0].equalsIgnoreCase("mode")){mode=pair[1];}
				else if(pair[0].equalsIgnoreCase("mapper")){mapper=pair[1];}
				else if(pair[0].equalsIgnoreCase("case")){name=pair[1];}
				else{throw new IllegalArgumentException("Unexpected argument: "+arg);}
			}
			if(!"baseline".equalsIgnoreCase(mode) && !"fixed".equalsIgnoreCase(mode)){
				throw new IllegalArgumentException("Expected explicit mode=baseline or mode=fixed.");
			}
			boolean assertions=false;
			assert(assertions=true);
			if(!assertions){throw new IllegalStateException("Constructor fraction checks require -ea.");}
			run(mode.equalsIgnoreCase("baseline"), Target.forName(mapper), Case.forName(name));
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Static fixture entry point only. */
	private MapperFractionOptionTest(){}

	/** Runs one constructor case and prints the observed state transition. */
	private static void run(final boolean baseline, final Target target, final Case kase) throws Exception{
		final Snapshot original=target.read();
		final ByteBuilder out=new ByteBuilder(512);
		out.append("mode\tmapper\tcase\targ\texpected\tobserved\tbefore\tafter\tpassed\n");
		try{
			target.setSentinel();
			final Snapshot before=target.read();
			if(!before.equals(SENTINEL)){
				throw new AssertionError(target.mapperName+" "+kase.name+" could not set sentinel: "+before);
			}
			final Throwable thrown=target.construct(kase.arg);
			final Snapshot after=target.read();
			final String observed=outcome(thrown);
			final String expected=expected(baseline, kase);
			final Snapshot expectedAfter=expectedAfter(kase);
			if(!observed.equals(expected)){
				throw new AssertionError(target.mapperName+" "+kase.name+" expected "+expected+" but observed "+observed, thrown);
			}
			if(!after.equals(expectedAfter)){
				throw new AssertionError(target.mapperName+" "+kase.name+" expected "+expectedAfter+" but observed "+after);
			}
			out.append(baseline ? "baseline" : "fixed").tab().append(target.mapperName).tab()
					.append(kase.name).tab().append(kase.arg).tab().append(expected).tab()
					.append(observed).tab().append(before.toString()).tab().append(after.toString())
					.tab().append(1).nl();
			System.out.print(out);
			System.out.println("MAPPER_FRACTION_OPTION mapper="+target.mapperName+" case="+kase.name+
					" mode="+(baseline ? "baseline" : "fixed")+" passed=true");
		}finally{target.write(original);}
	}

	/** Returns the exact expected outcome token for one mode and case. */
	private static String expected(final boolean baseline, final Case kase){
		if(kase==Case.FGTE_NAN){return baseline ? "return" : "IllegalArgumentException:fgte must not be NaN";}
		if(kase==Case.EF_NAN){return "NumberFormatException:For input string: \"nan\"";}
		return "return";
	}

	/** Returns the exact expected setter-controlled state after one case. */
	private static Snapshot expectedAfter(final Case kase){
		return kase==Case.FGTE_ZERO ? ZERO : SENTINEL;
	}

	/** Converts a constructor result to a stable, human-readable outcome token. */
	private static String outcome(final Throwable thrown){
		if(thrown==null){return "return";}
		final String message=thrown.getMessage();
		return thrown.getClass().getSimpleName()+":"+(message==null ? "" : message);
	}

	/** One command-line constructor case. */
	private enum Case{
		FGTE_NAN("nan", "fgte=NaN"),
		FGTE_ZERO("zero", "fgte=0"),
		FGTE_NEGATIVE("negative", "fgte=-1"),
		EF_NAN("efnan", "ef=NaN");

		Case(final String name_, final String arg_){
			name=name_;
			arg=arg_;
		}

		static Case forName(final String name){
			for(final Case c : values()){
				if(c.name.equalsIgnoreCase(name)){return c;}
			}
			throw new IllegalArgumentException("Unknown case: "+name);
		}

		final String name;
		final String arg;
	}

	/** Reflective handle for one mapper constructor and its paired index fields. */
	private static final class Target{
		Target(final String mapperName_, final String indexName_) throws ReflectiveOperationException{
			mapperName=mapperName_;
			indexName=indexName_;
			constructor=Class.forName("align2."+mapperName).getConstructor(String[].class);
			index=Class.forName("align2."+indexName);
			setter=index.getDeclaredMethod("setFractionToExclude", float.class);
			fraction=field("FRACTION_GENOME_TO_EXCLUDE");
			drop=field("MIN_INDEX_TO_DROP_LONG_HIT_LIST");
			average=field("MAX_AVERAGE_LIST_TO_SEARCH");
			average2=field("MAX_AVERAGE_LIST_TO_SEARCH2");
			single=field("MAX_SINGLE_LIST_TO_SEARCH");
			shortest=field("MAX_SHORTEST_LIST_TO_SEARCH");
		}

		/** Creates the mapper and unwraps any user-visible constructor failure. */
		Throwable construct(final String arg) throws ReflectiveOperationException{
			try{
				constructor.newInstance((Object)new String[] {arg});
				return null;
			}catch(final InvocationTargetException ite){
				return ite.getCause();
			}
		}

		/** Sets the same valid sentinel for every mapper/case JVM. */
		void setSentinel() throws ReflectiveOperationException{
			setter.invoke(null, 0.03f);
		}

		/** Reads the mutable static fields touched by the fraction setter. */
		Snapshot read() throws ReflectiveOperationException{
			return new Snapshot(fraction.getFloat(null), drop.getInt(null), average.getInt(null),
					average2.getInt(null), single.getInt(null), shortest.getInt(null));
		}

		/** Restores the mutable static fields before this JVM exits. */
		void write(final Snapshot s) throws ReflectiveOperationException{
			fraction.setFloat(null, s.fraction);
			drop.setInt(null, s.drop);
			average.setInt(null, s.average);
			average2.setInt(null, s.average2);
			single.setInt(null, s.single);
			shortest.setInt(null, s.shortest);
		}

		/** Finds one setter-controlled static field. */
		private Field field(final String fieldName) throws ReflectiveOperationException{
			final Field field=index.getDeclaredField(fieldName);
			field.setAccessible(true);
			return field;
		}

		static Target forName(final String name) throws ReflectiveOperationException{
			if(name==null){throw new IllegalArgumentException("Missing mapper.");}
			if(name.equalsIgnoreCase("bbmap5")){return new Target("BBMap5", "BBIndex5");}
			if(name.equalsIgnoreCase("acc") || name.equalsIgnoreCase("bbmapacc")){
				return new Target("BBMapAcc", "BBIndexAcc");
			}
			if(name.equalsIgnoreCase("pacbio") || name.equalsIgnoreCase("bbmappacbio")){
				return new Target("BBMapPacBio", "BBIndexPacBio");
			}
			if(name.equalsIgnoreCase("skimmer") || name.equalsIgnoreCase("bbmappacbioskimmer")){
				return new Target("BBMapPacBioSkimmer", "BBIndexPacBioSkimmer");
			}
			throw new IllegalArgumentException("Unknown mapper: "+name);
		}

		final String mapperName;
		final String indexName;
		final Constructor<?> constructor;
		final Class<?> index;
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

	private static final Snapshot SENTINEL=new Snapshot(0.03f, 895, 931, 958, 970, 916);
	private static final Snapshot ZERO=new Snapshot(0f, 1000, 1000, 1000, 1000, 1000);
}
