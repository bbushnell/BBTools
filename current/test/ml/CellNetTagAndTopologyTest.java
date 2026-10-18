package test.ml;

import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Locale;

import ml.CellNet;
import ml.CellNetParser;
import ml.Function;

/** Regression controls for CellNet custom tags and topology density selection. */
public final class CellNetTagAndTopologyTest{

	/** Runs one baseline/fixed assertion case, or prints a deterministic topology fixture. */
	public static void main(String[] args) throws Exception{
		if(args.length<2 || args.length>3){throw new IllegalArgumentException("baseline|fixed case [fresh-directory]");}
		final boolean fixed=args[0].equals("fixed");
		check(fixed || args[0].equals("baseline"), "Unknown phase");
		final String name=args[1];
		if(name.equals("tags")){
			tags(fixed);
		}else if(name.equals("parser")){
			parser(fixed);
		}else if(name.equals("copy")){
			copy(fixed);
		}else if(name.equals("setfrom")){
			setFrom(fixed);
		}else if(name.equals("layerdensity")){
			layerDensity(fixed);
		}else if(name.equals("locale")){
			locale(fixed);
		}else if(name.equals("dense-topology")){
			topology(args, true);
			return;
		}else if(name.equals("sparse-topology")){
			topology(args, false);
			return;
		}else if(name.equals("rawmap")){
			rawMap();
		}else{throw new IllegalArgumentException(name);}
		System.out.println("CellNetTagAndTopologyTest PASS "+name+" "+args[0]);
	}

	/** Repeated differently cased writes collapse to one value and latest spelling after the fix. */
	private static void tags(final boolean fixed){
		final CellNet net=network(true);
		net.setTag("MiXeD", "old");
		net.setTag("MIXED", "new");
		check("new".equals(net.getTag("mixed")), "lowercase lookup kept stale value");
		check("new".equals(net.getTag("MiXeD"))==fixed, "wrong mixed-case lookup behavior");
		assertSpelling(net, fixed);
	}

	/** The parser feeds duplicate mixed-case custom headers through the same canonical path. */
	private static void parser(final boolean fixed){
		final CellNet net=parse("##bbnet\n#version 1\n#concise\n#dense\n#density 1\n#blocksize 1\n#seed 1\n"+
			"#layers 2\n#dims 2 1\n##MiXeD old\n##MIXED new\n#coding decimal\nC3 LINEAR 0 1 1\n");
		check("new".equals(net.getTag("mixed")), "parser did not keep the lowercase duplicate value");
		check("new".equals(net.getTag("MiXeD"))==fixed, "wrong parser mixed-case lookup behavior");
		assertSpelling(net, fixed);
		final CellNet interleaved=parse("##bbnet\n#version 1\n#concise\n#dense\n#density 1\n#blocksize 1\n#seed 1\n"+
			"#layers 2\n#dims 2 1\n##MiXeD old\n##MIXED middle\n##MiXeD new\n#coding decimal\nC3 LINEAR 0 1 1\n");
		check((fixed ? "new" : "middle").equals(interleaved.getTag("mixed")),
			"parser did not preserve the expected duplicate value");
		final String reserved=parse("##bbnet\n#version 1\n#concise\n#dense\n#density 1\n#blocksize 1\n"+
			"#seed 1\n#layers 2\n#dims 2 1\n##FPR custom\n#coding decimal\nC3 LINEAR 0 1 1\n")
			.toBytes().toString();
		check(reserved.contains(fixed ? "##FPR custom\n" : "##fpr custom\n"),
			"mixed-case reserved-looking custom tag serialized unexpectedly:\n"+reserved);
	}

	/** copy() preserves canonical values and output spelling without sharing tag maps. */
	private static void copy(final boolean fixed){
		final CellNet source=network(true);
		source.setTag("MiXeD", "old");
		source.setTag("MIXED", "new");
		final CellNet copy=source.copy(false);
		source.setTag("mixed", "source-only");
		check("new".equals(copy.getTag("mixed")), "copy did not preserve the copied lowercase value");
		check("new".equals(copy.getTag("MiXeD"))==fixed, "wrong copy mixed-case lookup behavior");
		assertSpelling(copy, fixed);
	}

	/** setFrom() preserves canonical values and output spelling without sharing tag maps. */
	private static void setFrom(final boolean fixed){
		final CellNet source=network(true);
		source.setTag("MiXeD", "old");
		source.setTag("MIXED", "new");
		final CellNet dest=source.copy(false);
		dest.setTag("mixed", "dest-only");
		dest.setFrom(source, false);
		source.setTag("mixed", "source-only");
		check("new".equals(dest.getTag("mixed")), "setFrom did not preserve the copied lowercase value");
		check("new".equals(dest.getTag("MiXeD"))==fixed, "wrong setFrom mixed-case lookup behavior");
		assertSpelling(dest, fixed);
		dest.setFrom(dest, false);
		check("new".equals(dest.getTag("mixed")), "self setFrom erased tags");
	}

	/** Raw public map replacement remains supported when callers write lowercase canonical keys. */
	private static void rawMap(){
		final CellNet net=network(true);
		net.tags=new LinkedHashMap<String,String>();
		net.tags.put("mixed", "raw");
		check("raw".equals(net.getTag("MIXED")), "lowercase direct map write was not queryable");
		check(count(net.toBytes().toString(), "##mixed raw")==1, "lowercase direct map write was not serialized");
	}

	/** Locale.ROOT canonicalization is independent from Turkish default-locale lowercasing. */
	private static void locale(final boolean fixed){
		final Locale previous=Locale.getDefault();
		try{
			Locale.setDefault(new Locale("tr", "TR"));
			final CellNet set=network(true);
			set.setTag("TITLE", "setter");
			check("setter".equals(set.getTag("title"))==fixed, "wrong setter behavior under Turkish locale");
			final CellNet parsed=parse("##bbnet\n#version 1\n#concise\n#dense\n#density 1\n#blocksize 1\n"+
				"#seed 1\n#layers 2\n#dims 2 1\n##TITLE parser\n#coding decimal\nC3 LINEAR 0 1 1\n");
			check("parser".equals(parsed.getTag("title"))==fixed, "wrong parser behavior under Turkish locale");
		}finally{
			Locale.setDefault(previous);
		}
	}

	/** The helper now reports full final-layer density, density1 first, and density elsewhere. */
	private static void layerDensity(final boolean fixed) throws Exception{
		final Method method=CellNet.class.getDeclaredMethod("layerDensity",
			float.class, float.class, int.class, int.class);
		method.setAccessible(true);
		check(callLayerDensity(method, 0.25f, 0.75f, 3, 4)==(fixed ? 1f : 0.25f),
			"wrong final-layer density");
		check(callLayerDensity(method, 0.25f, 0.75f, 1, 4)==0.75f, "density1 no longer controls first hidden layer");
		check(callLayerDensity(method, 0.25f, 0.75f, 2, 4)==0.25f, "middle layer no longer uses default density");
	}

	/** Prints generated bytes for baseline/fixed comparison under the same runtime. */
	private static void topology(final String[] args, final boolean dense) throws Exception{
		if(args.length!=3){throw new IllegalArgumentException("topology cases require a fresh-directory");}
		final float[][] rates=functionRates();
		final boolean[] flags=flags();
		try{
			useLinearFunctionRates();
			CellNet.codingA48Out=true;
			CellNet.CONCISE=true;
			CellNet.OUT_DENSE=false;
			CellNet.OUT_SPARSE=false;
			CellNet.OUT_HEX=false;
			CellNet.OUT_DELTA_A48=false;
			final CellNet net=new CellNet(new int[] {4, 5, 3, 2}, 6789L, 0.35f, 0.65f, 2,
				new ArrayList<String>(), dense);
			net.randomize();
			final Path root=Paths.get(args[2]);
			Files.createDirectory(root);
			final Path out=root.resolve(dense ? "dense-topology.bbnet" : "sparse-topology.bbnet");
			Files.write(out, net.toBytes().toString().getBytes(StandardCharsets.US_ASCII));
			System.out.println("CellNetTagAndTopologyTest PASS "+args[1]+" "+args[0]);
		}finally{
			CellNet.codingA48Out=flags[0];
			CellNet.CONCISE=flags[1];
			CellNet.OUT_DENSE=flags[2];
			CellNet.OUT_SPARSE=flags[3];
			CellNet.OUT_HEX=flags[4];
			CellNet.OUT_DELTA_A48=flags[5];
			restoreFunctionRates(rates);
		}
	}

	/** Builds a small randomized net with the copied-row state needed by copy/setFrom. */
	private static CellNet network(final boolean dense){
		final float[][] rates=functionRates();
		try{
			useLinearFunctionRates();
			final CellNet net=new CellNet(new int[] {2, 3, 1}, 123L, 0.5f, 0.75f, 1,
				new ArrayList<String>(), dense);
			net.randomize();
			return net;
		}finally{
			restoreFunctionRates(rates);
		}
	}

	/** Parses literal fixture lines through the production parser. */
	private static CellNet parse(final String text){
		final ArrayList<byte[]> lines=new ArrayList<byte[]>();
		for(String line : text.split("\n")){lines.add(line.getBytes(StandardCharsets.US_ASCII));}
		return CellNetParser.loadFromLines(lines);
	}

	/** Invokes CellNet.layerDensity through reflection. */
	private static float callLayerDensity(final Method method, final float density, final float density1,
			final int layer, final int layers) throws Exception{
		return ((Float)method.invoke(null, density, density1, layer, layers)).floatValue();
	}

	/** Ensures serialization follows the phase's expected custom tag spelling. */
	private static void assertSpelling(final CellNet net, final boolean fixed){
		final String serialized=net.toBytes().toString();
		check(count(serialized, fixed ? "##MIXED new" : "##mixed new")==1,
			"expected tag spelling/value not emitted once:\n"+serialized);
		check(count(serialized, fixed ? "##mixed new" : "##MIXED new")==0,
			"unexpected tag spelling emitted:\n"+serialized);
		check(count(serialized, "##MiXeD old")==0, "stale mixed-case tag emitted:\n"+serialized);
	}

	/** Counts exact serialized lines. */
	private static int count(final String text, final String line){
		int count=0;
		for(String current : text.split("\n")){
			if(current.equals(line)){count++;}
		}
		return count;
	}

	/** Captures all six global writer settings touched by topology serialization. */
	private static boolean[] flags(){
		return new boolean[] {CellNet.codingA48Out, CellNet.CONCISE, CellNet.OUT_DENSE,
			CellNet.OUT_SPARSE, CellNet.OUT_HEX, CellNet.OUT_DELTA_A48};
	}

	/** Captures the shared activation-selection weights touched by fixture randomization. */
	private static float[][] functionRates(){
		return new float[][] {Function.TYPE_RATES.clone(),
			Function.TYPE_RATES_CUM==null ? null : Function.TYPE_RATES_CUM.clone()};
	}

	/** Selects deterministic hidden activations through Function's supported normalization path. */
	private static void useLinearFunctionRates(){
		Arrays.fill(Function.TYPE_RATES, 0);
		Function.TYPE_RATES[Function.LINEAR]=1;
		Function.TYPE_RATES_CUM=null;
		Function.normalizeTypeRates();
	}

	/** Restores shared activation-selection weights after fixture randomization. */
	private static void restoreFunctionRates(final float[][] rates){
		System.arraycopy(rates[0], 0, Function.TYPE_RATES, 0, rates[0].length);
		Function.TYPE_RATES_CUM=rates[1];
	}

	/** Throws independently of assertion enablement when a test oracle fails. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
