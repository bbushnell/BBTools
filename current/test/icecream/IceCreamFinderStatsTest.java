package test.icecream;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;

import icecream.IceCreamFinder;
import icecream.ZMWStreamer;
import json.JsonObject;
import shared.Timer;
import stream.ConcurrentReadOutputStream;

/** Literal fractional statistics and disabled-adapter construction oracles.
 * No worker or input stream is started; reflection preserves the private API.
 * @author Brian Bushnell, Shinobu
 */
public final class IceCreamFinderStatsTest{

	/** Runs with expect=baseline|fixed and outdir=new-directory.
	 * @param args Expected implementation and a fresh directory for the literal input
	 * @throws Exception if fixture setup or reflection fails */
	public static void main(String[] args) throws Exception{
		String expect=null;
		Path directory=null;
		for(String arg : args){
			if(arg.startsWith("expect=")){expect=arg.substring(7);}
			else if(arg.startsWith("outdir=")){directory=Paths.get(arg.substring(7));}
			else{throw new IllegalArgumentException("Unknown fixture argument: "+arg);}
		}
		if(directory==null || !("baseline".equals(expect) || "fixed".equals(expect))){
			throw new IllegalArgumentException("Use expect=baseline|fixed outdir=new-directory");
		}
		final boolean baseline=expect.equals("baseline");
		Files.createDirectory(directory);
		final Path input=directory.resolve("input.fq");
		Files.write(input, "@movie/1/0_4\nACGT\n+\nIIII\n".getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		final IceCreamFinder finder=new IceCreamFinder(new String[]{"in="+input, "threads=1",
				"parsecustom=t", "jni=f", "alignadapter=f", "alignrc=f", "trim=f"});
		set(finder, "readsProcessed", 9L); set(finder, "basesProcessed", 36L);
		set(finder, "readsOut", 4L); set(finder, "basesOut", 16L);
		set(finder, "ZMWs", 4L);
		set(finder, "truePositiveReads", 1L); set(finder, "trueNegativeReads", 1L);
		set(finder, "falsePositiveReads", 2L); set(finder, "falseNegativeReads", 1L);
		set(finder, "ambiguousReads", 4L);
		set(finder, "truePositiveZMWs", 1L); set(finder, "trueNegativeZMWs", 1L);
		set(finder, "falsePositiveZMWs", 1L); set(finder, "ambiguousZMWs", 1L);
		set(finder, "iceCreamRatio", 1.0); set(finder, "ratiosCounted", 1L);
		final Method textMethod=IceCreamFinder.class.getDeclaredMethod("toText", Timer.class);
		textMethod.setAccessible(true);
		final String text=textMethod.invoke(finder, new Timer()).toString();
		Files.write(directory.resolve("stats.txt"), text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		final Method jsonMethod=IceCreamFinder.class.getDeclaredMethod("toJson", Timer.class);
		jsonMethod.setAccessible(true);
		final JsonObject json=(JsonObject)jsonMethod.invoke(finder, new Timer());
		Files.write(directory.resolve("stats.json"), json.toString().getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		final ArrayList<Double> textRatios=new ArrayList<Double>();
		for(String line : text.split("\n")){
			if(line.startsWith("SNR:")){textRatios.add(Double.parseDouble(line.substring(4).trim()));}
		}
		if(textRatios.size()!=2){throw new AssertionError("Expected separate read and ZMW SNR lines");}
		System.out.println("case\texpected\tactual");
		// Literal decibel values for read ratio6/7 and well ratio3/2.
		check("text-reads", baseline ? Double.NEGATIVE_INFINITY : -0.67, textRatios.get(0));
		check("text-wells", baseline ? 4.77 : 1.76, textRatios.get(1));
		check("json-reads", baseline ? Double.NEGATIVE_INFINITY : -0.6694678963061322, json.getDouble("SNR_Reads"));
		check("json-wells", baseline ? 0.0 : 1.7609125905568124, json.getDouble("SNR_ZMWs"));
		checkDisabledAdapter(input, baseline);
	}

	/** Assigns existing counters without executing the sequence classifier. */
	private static void set(IceCreamFinder finder, String name, Object value) throws Exception{
		final Field field=IceCreamFinder.class.getDeclaredField(name);
		field.setAccessible(true);
		field.set(finder, value);
	}

	/** Compares finite numeric or explicitly observed baseline-infinity results. */
	private static void check(String name, double expected, Double actual){
		System.out.println(name+"\t"+expected+"\t"+actual);
		if(actual==null || !(Double.compare(expected, actual)==0 || Math.abs(expected-actual)<1e-12)){
			throw new AssertionError(name+": expected "+expected+", found "+actual);
		}
	}

	/** Disabling adapter matching must still allow a trimming worker to be built. */
	private static void checkDisabledAdapter(Path input, boolean baseline) throws Exception{
		final IceCreamFinder finder=new IceCreamFinder(new String[]{"in="+input,
				"threads=1", "adapter=null", "trim=t", "jni=f"});
		final Class<?> worker=Class.forName("icecream.IceCreamFinder$ProcessThread");
		final Constructor<?> constructor=worker.getDeclaredConstructor(IceCreamFinder.class,
				ZMWStreamer.class, ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class,
				ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class, int.class);
		constructor.setAccessible(true);
		boolean nullFailure=false;
		try{
			constructor.newInstance(finder, null, null, null, null, null, 0);
		}catch(InvocationTargetException e){
			if(!(e.getCause() instanceof NullPointerException)){throw e;}
			nullFailure=true;
		}
		System.out.println("disabled-adapter-null-failure\t"+baseline+"\t"+nullFailure);
		if(nullFailure!=baseline){throw new AssertionError("Unexpected disabled-adapter worker-construction outcome");}
	}

	/** Static-only fixture. */
	private IceCreamFinderStatsTest(){}
}
