package test.tracker;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;

import icecream.ZMW;
import stream.Read;
import tracker.EntropyTracker;

/** Invokes the actual two caller predicates without starting their I/O workers.
 * Reflection is limited to the callers' private worker constructor/method so
 * the production API need not be broadened for a test.
 * @author Brian Bushnell, Jean */
public final class EntropyCallerProbe {

	public static void main(String[] args) throws Exception{
		assert(args.length==2 && args[0].startsWith("in=") &&
			(args[1].equals("mode=observe") || args[1].equals("mode=verify"))) :
			"Use in=<existing small FASTA fixture> and mode=observe|verify.";
		System.out.println("caller\tcase\tthreshold\texpected\tactual\tpass");
		for(String name : new String[]{"icecream.IceCreamFinder", "icecream.ReformatPacBio"}){
			Class<?> outer=Class.forName(name);
			Object owner=outer.getConstructor(String[].class).newInstance((Object)new String[]{
				args[0], "out=null", "entropy=0.5", "entropyk=2", "entropywindow=4", "mmf=0.75", "t=1"});
			Class<?> worker=Class.forName(name+"$ProcessThread");
			Constructor<?>[] constructors=worker.getDeclaredConstructors();
			assert(constructors.length==1) : "Worker fixture assumes one constructor: "+name;
			Constructor<?> ctor=constructors[0];
			Class<?>[] types=ctor.getParameterTypes();
			assert(types[0]==outer && types[types.length-1]==int.class) : "Unexpected worker signature: "+name;
			Object[] parameters=new Object[types.length];
			parameters[0]=owner;
			parameters[types.length-1]=Integer.valueOf(0);
			ctor.setAccessible(true);
			Object instance=ctor.newInstance(parameters);
			Field field=worker.getDeclaredField("eTracker");
			field.setAccessible(true);
			EntropyTracker tracker=(EntropyTracker)field.get(instance);
			assert(tracker.k()==2 && tracker.windowBases()==4 && tracker.cutoff()==0.5f) :
				"Caller must use the declared JT006 settings, not unverified parser defaults: "+name;
			Method predicate=worker.getDeclaredMethod("flagLowEntropyReads", ZMW.class, float.class, int.class, float.class);
			predicate.setAccessible(true);
			check(name, "gap_at_six", instance, predicate, "AAAAAANNNNAAAAAA", 6, true);
			check(name, "gap_at_seven", instance, predicate, "AAAAAANNNNAAAAAA", 7, false);
			check(name, "gap_at_eight", instance, predicate, "AAAAAANNNNAAAAAA", 8, false);
			check(name, "gap_at_nine", instance, predicate, "AAAAAANNNNAAAAAA", 9, false);
			check(name, "gap_at_ten", instance, predicate, "AAAAAANNNNAAAAAA", 10, false);
			check(name, "continuous_control", instance, predicate, "AAAAAAAAAAAAAAAA", 8, true);
			check(name, "high_gap_control", instance, predicate, "AAAAAACGTCGAAAAAA", 8, false);
			check(name, "all_N_control", instance, predicate, "NNNNNN", 1, false);
		}
		System.err.println("caller_checks="+checks+" failures="+failures);
		if(args[1].equals("mode=verify") && failures>0){throw new AssertionError("Caller failures: "+failures);}
	}

	/** Thresholds6/7/8/9/10 surround the independently agreed six-base maximum. */
	private static void check(String caller, String name, Object worker, Method predicate,
			String sequence, int threshold, boolean discard) throws Exception{
		assert(threshold>0 && threshold<=sequence.length()) : "Use an absolute boundary unaffected by the fractional cap.";
		Read read=new Read(sequence.getBytes(StandardCharsets.US_ASCII), null, "movie/1/0_"+sequence.length(), 1);
		ZMW reads=new ZMW(1);
		reads.add(read);
		final int found=((Integer)predicate.invoke(worker, reads, 0.5f, threshold, 1.0f)).intValue();
		final String expected=(discard ? "1,true,true" : "0,false,false");
		final String actual=found+","+read.discarded()+","+read.junk();
		final boolean pass=expected.equals(actual);
		checks++;
		if(!pass){failures++;}
		System.out.println(caller+"\t"+name+"\t"+threshold+"\t"+expected+"\t"+actual+"\t"+pass);
	}

	private static int checks, failures;
}
