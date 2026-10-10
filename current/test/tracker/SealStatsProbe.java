package test.tracker;

import tracker.SealStats;
import tracker.SealStats.SealStatsLine;

/** Literal JT009 observations using real empty/singleton/multi-entry statistics files.
 * Empty countNonprimary expectations are a proposal, not an approved contract.
 * @author Brian Bushnell, Jean */
public final class SealStatsProbe {

	public static void main(String[] args){
		assert(args.length==1 && args[0].startsWith("fixtures=")) : "JT009 requires the committed fixture directory.";
		final String root=args[0].substring("fixtures=".length());
		System.out.println("case\tproposed_or_control\tobserved\tagrees");
		final SealStats empty=new SealStats(root+"/tracker_seal_empty.tsv");
		show("empty_primary_control", "null", describe(empty.primary()));
		show("empty_nonmatching_control", "!sample,0,0", describe(empty.countNonmatching("sample")));
		show("empty_nonprimary_proposal", "!sample,0,0", nonprimary(empty, "sample"));
		show("empty_null_name_proposal", "!null,0,0", nonprimary(empty, null));

		final SealStats single=new SealStats(root+"/tracker_seal_singleton.tsv");
		show("singleton_primary_control", "alpha,5,50", describe(single.primary()));
		show("singleton_nonprimary_control", "!alpha,0,0", nonprimary(single, "sample"));
		show("singleton_missing_control", "!missing,5,50", describe(single.countNonmatching("missing")));

		final SealStats multi=new SealStats(root+"/tracker_seal_multi.tsv");
		show("multi_primary_control", "alpha,5,50", describe(multi.primary()));
		show("multi_nonprimary_control", "!alpha,5,90", nonprimary(multi, "sample"));
		show("multi_unused_name_control", "!alpha,5,90", nonprimary(multi, "beta"));
		show("multi_missing_control", "!missing,10,140", describe(multi.countNonmatching("missing")));
		show("multi_totals_control", "10,140,10,140",
			multi.totalReads+","+multi.totalBases+","+multi.matchedReads+","+multi.matchedBases);
		final SealStatsLine aggregate=multi.countNonprimary("sample");
		show("multi_repeated_control", "!alpha,5,90", describe(aggregate));
		aggregate.reads+=100;
		show("multi_sources_control", "alpha,5,50;beta,3,60;gamma,2,30",
			describe(multi.map.get("alpha"))+";"+describe(multi.map.get("beta"))+";"+describe(multi.map.get("gamma")));
	}

	/** Record the known empty-map failure while letting other failures propagate. */
	private static String nonprimary(SealStats stats, String name){
		try{return describe(stats.countNonprimary(name));}
		catch(NullPointerException e){return "NullPointerException";}
	}

	private static String describe(SealStatsLine line){
		return line==null ? "null" : line.name+","+line.reads+","+line.bases;
	}

	private static void show(String name, String expected, String observed){
		assert(name!=null && expected!=null && observed!=null) : "Each JT009 observation needs its literal named expectation.";
		System.out.println(name+"\t"+expected+"\t"+observed+"\t"+expected.equals(observed));
	}
}
