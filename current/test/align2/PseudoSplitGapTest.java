package test.align2;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;

import align2.BBMapS;
import align2.HybridMaxIndelConfig;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.SiteScore;
import structures.ByteBuilder;

/**
 * Checks the finite pseudo-split helper boundary and the ordered retrymaxindel
 * parse sequence without loading an index or running mapper workers.
 * @author Brian Bushnell, Jean
 */
public final class PseudoSplitGapTest{

	/** Accepts explicit mode=baseline|fixed and case=helper|options arguments. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, PseudoSplitGapTest.class, false);
		try{
			String mode=null, name=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				if(pair.length!=2){throw new IllegalArgumentException("Expected key=value argument: "+arg);}
				if(pair[0].equalsIgnoreCase("mode")){mode=pair[1];}
				else if(pair[0].equalsIgnoreCase("case")){name=pair[1];}
				else{throw new IllegalArgumentException("Unexpected argument: "+arg);}
			}
			if(!"baseline".equalsIgnoreCase(mode) && !"fixed".equalsIgnoreCase(mode)){
				throw new IllegalArgumentException("Expected explicit mode=baseline or mode=fixed.");
			}
			boolean assertions=false;
			assert(assertions=true);
			if(!assertions){throw new IllegalStateException("Pseudo split checks require -ea.");}
			final boolean baseline=mode.equalsIgnoreCase("baseline");
			final Case kase=Case.forName(name);
			if(kase==Case.HELPER){runHelperSuite(baseline);}
			else if(kase==Case.OPTIONS){runOptionsCase(baseline);}
			else{throw new IllegalStateException("Unhandled case: "+kase);}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Static fixture entry point only. */
	private PseudoSplitGapTest(){}

	/** Runs the private helper against one finite table of synthetic site lists. */
	private static void runHelperSuite(final boolean baseline) throws Exception{
		final Method helper=Class.forName("align2.BBMapThread")
				.getDeclaredMethod("smallestPseudoSplitGap", Read.class, int.class);
		helper.setAccessible(true);
		final ByteBuilder out=new ByteBuilder(512);
		out.append("mode\tcase\tmax_gap\texpected\tobserved\tpassed\n");
		for(final HelperCase kase : HelperCase.values()){
			final int observed=((Integer)helper.invoke(null, kase.read(), kase.maxGap)).intValue();
			final int expected=kase.expected(baseline);
			if(observed!=expected){
				throw new AssertionError(kase.name+" expected "+expected+" but observed "+observed);
			}
			out.append(baseline ? "baseline" : "fixed").tab().append(kase.name).tab()
					.append(kase.maxGap).tab().append(expected).tab().append(observed).tab()
					.append(1).nl();
		}
		System.out.print(out);
		System.out.println("PSEUDO_SPLIT_GAP case=helper mode="+
				(baseline ? "baseline" : "fixed")+" checks=9 failed=0");
	}

	/** Verifies the ordered parser bypass and the literal immutable config boundary. */
	private static void runOptionsCase(final boolean baseline) throws Exception{
		final BBMapS mapper=new BBMapS(new String[] {
				"hybridmaxindel=t",
				"maxindel=50",
				"maxindel2=100",
				"retrymaxindel2=2147483647",
				"retrymaxindel=2147483647"
		});
		final int retry1=intField("retryMaxIndel1").getInt(mapper);
		final int retry2=intField("retryMaxIndel2").getInt(mapper);
		if(retry1!=Integer.MAX_VALUE || retry2!=Integer.MAX_VALUE){
			throw new AssertionError("Expected parsed retry MAX/MAX but observed "+retry1+"/"+retry2);
		}
		final HybridMaxIndelConfig parsedConfig=(HybridMaxIndelConfig)
				configField().get(mapper);
		if(parsedConfig==null || parsedConfig.lowPrimary!=50 ||
				parsedConfig.lowSum!=100 ||
				parsedConfig.highPrimary!=Integer.MAX_VALUE ||
				parsedConfig.highSum!=Integer.MAX_VALUE ||
				parsedConfig.retryMinMapq!=0){
			throw new AssertionError("Unexpected parsed config: "+parsedConfig);
		}
		final HybridMaxIndelConfig config=
				new HybridMaxIndelConfig(50, 100, retry1, retry2);
		if(config.lowPrimary!=50 || config.lowSum!=100 ||
				config.highPrimary!=Integer.MAX_VALUE || config.highSum!=Integer.MAX_VALUE ||
				config.retryMinMapq!=0){
			throw new AssertionError("Unexpected config: "+config+
					", retryMinMapq="+config.retryMinMapq);
		}
		final ByteBuilder out=new ByteBuilder(256);
		out.append("mode\tcase\tretry1\tretry2\tlow_primary\tlow_sum\thigh_primary\thigh_sum")
				.append("\tparsed_high_primary\tparsed_high_sum\tpassed\n")
				.append(baseline ? "baseline" : "fixed").tab().append("options").tab()
				.append(retry1).tab().append(retry2).tab()
				.append(config.lowPrimary).tab().append(config.lowSum).tab()
				.append(config.highPrimary).tab().append(config.highSum).tab()
				.append(parsedConfig.highPrimary).tab().append(parsedConfig.highSum).tab()
				.append(1).nl();
		System.out.print(out);
		System.out.println("PSEUDO_SPLIT_GAP case=options mode="+
				(baseline ? "baseline" : "fixed")+" checks=1 failed=0");
	}

	/** Returns one reflected package-private AbstractMapper option field. */
	private static Field intField(final String name) throws ReflectiveOperationException{
		final Field field=Class.forName("align2.AbstractMapper").getDeclaredField(name);
		field.setAccessible(true);
		return field;
	}

	/** Returns the reflected package-private AbstractMapper per-invocation config field. */
	private static Field configField() throws ReflectiveOperationException{
		final Field field=Class.forName("align2.AbstractMapper")
				.getDeclaredField("hybridMaxIndelConfig");
		field.setAccessible(true);
		return field;
	}

	/** Builds one synthetic site for the private helper. */
	private static SiteScore site(final int start, final int left, final int right){
		final SiteScore ss=new SiteScore(1, (byte)0, start, start, 1, 1);
		ss.setPseudoCoverage(left, right);
		return ss;
	}

	/** Builds a read with the supplied synthetic sites. */
	private static Read read(final SiteScore... sites){
		final Read r=new Read(new byte[] {(byte)'A'}, null, "synthetic", 1);
		if(sites!=null){
			r.sites=new ArrayList<SiteScore>(sites.length);
			for(final SiteScore ss : sites){r.sites.add(ss);}
		}
		return r;
	}

	/** One top-level fixture mode, each run in a separate JVM by the outer recipe. */
	private enum Case{
		HELPER("helper"),
		OPTIONS("options");

		Case(final String name_){name=name_;}

		static Case forName(final String name){
			for(final Case c : values()){
				if(c.name.equalsIgnoreCase(name)){return c;}
			}
			throw new IllegalArgumentException("Unknown case: "+name);
		}

		final String name;
	}

	/** One private-helper expectation. */
	private enum HelperCase{
		NULL_NO_HIT("null", 100, 0, 0, readNullSites()),
		ONE_SITE("one_site", 100, 0, 0,
				PseudoSplitGapTest.read(site(0, 100, 0))),
		GAP50_EXCLUDED("gap50", 50, 0, 0,
				PseudoSplitGapTest.read(site(0, 100, 0), site(50, 0, 100))),
		MIN51_OF_51_75("min51", 100, 51, 51,
				PseudoSplitGapTest.read(site(0, 100, 0), site(51, 0, 100),
				site(75, 0, 100))),
		ENDPOINT77_IGNORES78("endpoint77", 77, 77, 77,
				PseudoSplitGapTest.read(site(0, 100, 0), site(77, 0, 100),
				site(78, 0, 100))),
		MAX_MINUS_ONE("max_minus_one", Integer.MAX_VALUE-1, Integer.MAX_VALUE-1,
				Integer.MAX_VALUE-1, PseudoSplitGapTest.read(site(0, 100, 0),
				site(Integer.MAX_VALUE-1, 0, 100))),
		MAX_NO_VALID_COVERAGE("max_no_valid_coverage", Integer.MAX_VALUE,
				Integer.MIN_VALUE, 0,
				PseudoSplitGapTest.read(site(0, 94, 0), site(Integer.MAX_VALUE, 0, 94))),
		MAX_GAP51("max_gap51", Integer.MAX_VALUE, Integer.MIN_VALUE, 51,
				PseudoSplitGapTest.read(site(0, 100, 0), site(51, 0, 100))),
		MAX_ACTUAL("max_actual", Integer.MAX_VALUE, Integer.MIN_VALUE,
				Integer.MAX_VALUE,
				PseudoSplitGapTest.read(site(0, 100, 0), site(Integer.MAX_VALUE, 0, 100)));

		HelperCase(final String name_, final int maxGap_, final int baseline_,
				final int fixed_, final Read read_){
			name=name_;
			maxGap=maxGap_;
			baseline=baseline_;
			fixed=fixed_;
			read=read_;
		}

		Read read(){return read;}

		int expected(final boolean old){return old ? baseline : fixed;}

		final String name;
		final int maxGap;
		final int baseline;
		final int fixed;
		final Read read;
	}

	/** Builds a read whose candidate list is absent. */
	private static Read readNullSites(){
		final Read r=new Read(new byte[] {(byte)'A'}, null, "synthetic", 1);
		r.sites=null;
		return r;
	}
}
