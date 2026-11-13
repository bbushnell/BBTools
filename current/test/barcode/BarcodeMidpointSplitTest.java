package test.barcode;

import java.util.Locale;

import barcode.Barcode;
import barcode.BarcodeStats;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Checks literal midpoint slices and the expected/observed BarcodeStats callers.
 * Baseline mode qualifies the original wrong slices and short-name exceptions;
 * fixed mode requires every input base to survive, with odd remainders on the right.
 * Positive-delimiter controls retain their explicit, possibly unequal boundaries.
 * @author Brian Bushnell, Sucrose
 */
public final class BarcodeMidpointSplitTest{

	/** Requires mode=baseline|fixed; case=all|api|stats|metadata selects independent finite checks. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, BarcodeMidpointSplitTest.class, false);
		try{
			String mode=null, testCase="all";
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(key.equals("case")){testCase=value;}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			check("baseline".equals(mode) || "fixed".equals(mode), "Set mode=baseline or mode=fixed.");
			check("all".equals(testCase) || "api".equals(testCase) || "stats".equals(testCase)
					|| "metadata".equals(testCase), "Set case=all, api, stats or metadata.");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Midpoint checks require -ea for production caller preconditions.");
			Shared.setThreads(1);
			final boolean baseline="baseline".equals(mode);
			if("all".equals(testCase) || "api".equals(testCase)){checkApi(baseline);}
			if("all".equals(testCase) || "stats".equals(testCase)){checkStats(baseline);}
			if("all".equals(testCase) || "metadata".equals(testCase)){checkMetadata(baseline);}
			System.out.println("BARCODE_MIDPOINT mode="+mode+" case="+testCase+" failed=0");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Expected strings are literals, independent of the production split arithmetic. */
	private static void checkApi(final boolean baseline){
		checkPair("empty", "", (char)0, "", "", null, null, baseline);
		checkPair("single", "A", (char)0, "", "A", "A", null, baseline);
		checkPair("short-even", "AC", (char)0, "A", "C", "AC", "", baseline);
		checkPair("short-odd", "ACG", (char)0, "A", "CG", "AC", "", baseline);
		checkPair("even", "AAAACCCC", (char)0, "AAAA", "CCCC", "AAAAC", "CCC", baseline);
		checkPair("odd", "AAAACCCCC", (char)0, "AAAA", "CCCCC", "AAAAC", "CCC", baseline);
		checkPair("delimiter-even", "AAAA+CCCC", '+', "AAAA", "CCCC", "AAAA", "CCCC", baseline);
		checkPair("delimiter-unequal", "AA+CCCC", '+', "AA", "CCCC", "AA", "CCCC", baseline);
		System.out.println("api_cases=8");
	}

	/** Exercises both direct strings and count-carrying child construction for each side. */
	private static void checkPair(final String label, final String name, final char delimiter,
			final String fixedLeft, final String fixedRight, final String oldLeft, final String oldRight,
			final boolean baseline){
		final Barcode source=new Barcode(name, 7);
		final String left=baseline ? oldLeft : fixedLeft, right=baseline ? oldRight : fixedRight;
		checkSlice(source, delimiter, 0, left, false, label);
		checkSlice(source, delimiter, 1, right, false, label);
		checkSlice(source, delimiter, 0, left, true, label);
		checkSlice(source, delimiter, 1, right, true, label);
		if(!baseline || delimiter>0){
			final String joined=delimiter>0 ? left+delimiter+right : left+right;
			check(name.equals(joined), label+": extracted components must preserve the original name.");
		}
		check(source.name.equals(name) && source.count()==7, label+": extraction mutated its source.");
		System.out.println(label+" left="+show(left)+" right="+show(right));
	}

	/** Null expectations mean a bounds failure specifically from the production extraction helper. */
	private static void checkSlice(final Barcode source, final char delimiter, final int pairnum,
			final String expected, final boolean child, final String label){
		String actual=null;
		StringIndexOutOfBoundsException failure=null;
		try{
			if(child){
				final Barcode result=source.getBarcodeForPairnum(pairnum, delimiter);
				actual=result.name;
				check(result!=source && result.count()==7, label+": child must copy the count into an independent object.");
			}else{actual=source.getStringForPairnum(pairnum, delimiter);}
		}catch(final StringIndexOutOfBoundsException e){failure=e;}
		check((failure!=null)==(expected==null), label+": pair="+pairnum+", child="+child+
				", expected="+show(expected)+", observed="+show(actual)+", failure="+failure);
		if(failure!=null){
			check(splitFrame(failure), label+": bounds failure must originate in Barcode.getStringForPairnum.");
		}else{check(expected.equals(actual), label+": pair="+pairnum+", expected="+show(expected)+", observed="+show(actual));}
	}

	/** Public Stats construction reaches the fallback for expected names and weighted observations. */
	private static void checkStats(final boolean baseline){
		checkStatsPair("stats-even", "AAAACCCC", (char)0, 4, 4, 4, 4,
				baseline ? "AAAAC" : "AAAA", baseline ? "CCC" : "CCCC");
		checkStatsPair("stats-odd", "AAAACCCCC", (char)0, 4, 5, 4, 5,
				baseline ? "AAAAC" : "AAAA", baseline ? "CCC" : "CCCCC");
		checkStatsPair("stats-delimiter", "AA+CCCC", '+', 2, 4, 2, 4, "AA", "CCCC");
		System.out.println("stats_cases=3");
	}

	/** Checks inferred lengths for uniform names, mixed component lengths and empty data. */
	private static void checkMetadata(final boolean baseline){
		checkStatsPair("metadata-unset", "AAAACCCCC", (char)0, 0, 0,
				baseline ? 0 : 4, baseline ? 0 : 5,
				baseline ? "AAAAC" : "AAAA", baseline ? "CCC" : "CCCCC");
		checkStatsPair("metadata-conflict", "AAAACCCC", (char)0, 3, 5,
				baseline ? 3 : 4, baseline ? 5 : 4,
				baseline ? "AAAAC" : "AAAA", baseline ? "CCC" : "CCCC");
		checkMixedMetadata(baseline, true);
		checkMixedMetadata(baseline, false);
		checkEmptyMetadata(baseline);
		System.out.println("metadata_cases=5");
	}

	/** Each collection can independently make an otherwise uniform component length unknown. */
	private static void checkMixedMetadata(final boolean baseline, final boolean observedMixed){
		final String label=observedMixed ? "metadata-observed-mixed" : "metadata-expected-mixed";
		final BarcodeStats parent=new BarcodeStats(0, 2, label);
		parent.length1=3;
		parent.length2=5;
		parent.addExpectedCode("AAAACCCC");
		parent.increment("AAAACCCC", 7);
		if(observedMixed){parent.increment("AAAATGGGGG", 11);}
		else{parent.addExpectedCode("AAAACCCCC");}
		final BarcodeStats left=parent.makeLeft(), right=parent.makeRight();
		final String firstLeft=baseline ? "AAAAC" : "AAAA", firstRight=baseline ? "CCC" : "CCCC";
		if(observedMixed){
			checkLayout(left, baseline ? 3 : 0,
					new String[]{firstLeft, baseline ? "AAAATG" : "AAAAT"}, new long[]{7, 11},
					new String[]{firstLeft}, label+" left");
			checkLayout(right, baseline ? 5 : 0,
					new String[]{firstRight, baseline ? "GGGG" : "GGGGG"}, new long[]{7, 11},
					new String[]{firstRight}, label+" right");
		}else{
			checkLayout(left, baseline ? 3 : 4, new String[]{firstLeft}, new long[]{7},
					new String[]{firstLeft}, label+" left");
			checkLayout(right, baseline ? 5 : 0, new String[]{firstRight}, new long[]{7},
					baseline ? new String[]{"CCC"} : new String[]{"CCCC", "CCCCC"}, label+" right");
		}
		check(parent.length1==3 && parent.length2==5 && parent.codeMap.get("AAAACCCC").count()==7,
				label+": parent metadata or first observation changed.");
		check(parent.codeMap.size()==(observedMixed ? 2 : 1) && parent.expectedCodeList.size()==(observedMixed ? 1 : 2),
				label+": parent membership changed.");
		if(observedMixed){check(parent.codeMap.get("AAAATGGGGG").count()==11, label+": parent second count changed.");}
		for(final Barcode expected : parent.expectedCodeList){check(expected.count()==0, label+": expected source count changed.");}
		System.out.println(label+" child_lengths="+left.length1+","+right.length1);
	}

	/** No extracted names means unknown lengths; parent metadata cannot supply a midpoint. */
	private static void checkEmptyMetadata(final boolean baseline){
		final BarcodeStats parent=new BarcodeStats(0, 2, "metadata-empty");
		parent.length1=3;
		parent.length2=5;
		final BarcodeStats left=parent.makeLeft(), right=parent.makeRight();
		checkLayout(left, baseline ? 3 : 0, new String[0], new long[0], new String[0], "metadata-empty left");
		checkLayout(right, baseline ? 5 : 0, new String[0], new long[0], new String[0], "metadata-empty right");
		check(parent.length1==3 && parent.length2==5 && parent.codeMap.isEmpty() && parent.expectedCodeList.isEmpty(),
				"Empty splitting must not mutate its parent.");
		System.out.println("metadata-empty child_lengths="+left.length1+","+right.length1);
	}

	/** Literal memberships and weights distinguish metadata derived from only one collection. */
	private static void checkLayout(final BarcodeStats child, final int length, final String[] observed,
			final long[] counts, final String[] expected, final String label){
		check(child.barcodesPerRead==1 && child.length1==length && child.length2==0, label+": incorrect child lengths.");
		check(child.codeMap.size()==observed.length && child.expectedCodeList.size()==expected.length
				&& child.expectedCodeMap.size()==expected.length, label+": incorrect component membership.");
		for(int i=0; i<observed.length; i++){
			final Barcode actual=child.codeMap.get(observed[i]);
			check(actual!=null && actual.count()==counts[i], label+": incorrect count for "+observed[i]);
		}
		for(final String name : expected){
			final Barcode actual=child.expectedCodeMap.get(name);
			check(actual!=null && actual.count()==0 && child.expectedCodeList.contains(actual),
					label+": incorrect expected component "+name);
		}
	}

	/** Checks literal keys, retained weights/length metadata and independent child objects. */
	private static void checkStatsPair(final String label, final String name, final char delimiter,
			final int parentLeftLength, final int parentRightLength,
			final int leftLength, final int rightLength, final String leftName, final String rightName){
		final BarcodeStats parent=new BarcodeStats(delimiter, 2, label);
		parent.length1=parentLeftLength;
		parent.length2=parentRightLength;
		parent.addExpectedCode(name);
		parent.increment(name, 7);
		final BarcodeStats left=parent.makeLeft(), right=parent.makeRight();
		checkStatsChild(left, parent, leftName, leftLength, label+" left");
		checkStatsChild(right, parent, rightName, rightLength, label+" right");
		check(parent.codeMap.size()==1 && parent.codeMap.get(name).count()==7,
				label+": splitting must retain the original observation.");
		check(parent.expectedCodeList.size()==1 && parent.expectedCodeList.get(0).name.equals(name)
				&& parent.expectedCodeList.get(0).count()==0, label+": splitting mutated the expected source.");
		check(parent.length1==parentLeftLength && parent.length2==parentRightLength,
				label+": splitting mutated the parent lengths.");
		System.out.println(label+" left="+show(leftName)+" right="+show(rightName)+" observations=7"+
				" parent_lengths="+parent.length1+","+parent.length2+" child_lengths="+left.length1+","+right.length1);
	}

	/** Expected and observed entries stay independent and retain their separate initial counts. */
	private static void checkStatsChild(final BarcodeStats child, final BarcodeStats parent,
			final String name, final int length, final String label){
		check(child!=parent && child.barcodesPerRead==1 && child.length1==length && child.length2==0,
				label+": child layout metadata changed.");
		check(child.codeMap.size()==1 && child.codeMap.containsKey(name)
				&& child.codeMap.get(name).count()==7, label+": incorrect observed name or count.");
		check(child.expectedCodeList.size()==1 && child.expectedCodeMap.size()==1
				&& child.expectedCodeMap.get(name)==child.expectedCodeList.get(0), label+": incorrect expected name.");
		check(child.expectedCodeList.get(0).count()==0 && child.expectedCodeList.get(0)!=child.codeMap.get(name)
				&& child.expectedCodeList.get(0)!=parent.expectedCodeList.get(0)
				&& child.codeMap.get(name)!=parent.codeMap.values().iterator().next(), label+": child objects must be independent.");
	}

	/** Rejects exceptions from fixture code or a different production method. */
	private static boolean splitFrame(final Throwable failure){
		for(final StackTraceElement frame : failure.getStackTrace()){
			if(frame.getClassName().equals("barcode.Barcode") && frame.getMethodName().equals("getStringForPairnum")){return true;}
		}
		return false;
	}

	/** Distinguishes empty components from expected exceptions in saved output. */
	private static String show(final String value){return value==null ? "<bounds>" : "'"+value+"'";}

	/** Keeps literal oracle failures visible even if downstream assertions are disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	private BarcodeMidpointSplitTest(){}
}
