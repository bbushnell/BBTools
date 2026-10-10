package assemble;

/** Exact parser and dispatcher regressions for multi-k Tadpole. */
public class TadpoleMultiUnitTest {

	public static void main(String[] args){
		int failures=0;
		failures+=run("orderedPhasesPreserveOccurrences", TadpoleMultiUnitTest::orderedPhasesPreserveOccurrences);
		failures+=run("orderedGraphDefaultsAndOverrides", TadpoleMultiUnitTest::orderedGraphDefaultsAndOverrides);
		failures+=run("orderedPhaseFilters", TadpoleMultiUnitTest::orderedPhaseFilters);
		failures+=run("shorthandExpandsFromLargest", TadpoleMultiUnitTest::shorthandExpandsFromLargest);
		failures+=run("explicitAssemblyReclassifiesShorthand", TadpoleMultiUnitTest::explicitAssemblyReclassifiesShorthand);
		failures+=run("explicitPhaseListsOverrideShorthand", TadpoleMultiUnitTest::explicitPhaseListsOverrideShorthand);
		failures+=run("fusionMismatchOption", TadpoleMultiUnitTest::fusionMismatchOption);
		failures+=run("fusionDeadEndOption", TadpoleMultiUnitTest::fusionDeadEndOption);
		failures+=run("fusionConflictOption", TadpoleMultiUnitTest::fusionConflictOption);
		failures+=run("fusionCoverageRatioOption", TadpoleMultiUnitTest::fusionCoverageRatioOption);
		failures+=run("finalGraphCoverageRatioOption", TadpoleMultiUnitTest::finalGraphCoverageRatioOption);
		failures+=run("singleKHasNoImplicitJoining", TadpoleMultiUnitTest::singleKHasNoImplicitJoining);
		failures+=run("noneDisablesPhases", TadpoleMultiUnitTest::noneDisablesPhases);
		failures+=run("explicitAutoUsesShorthand", TadpoleMultiUnitTest::explicitAutoUsesShorthand);
		failures+=run("duplicateKValuesAreDeduplicated", TadpoleMultiUnitTest::duplicateKValuesAreDeduplicated);
		failures+=run("defaultGraphKIsAssemblyK", TadpoleMultiUnitTest::defaultGraphKIsAssemblyK);
		failures+=run("explicitAutoGraphKIsAssemblyK", TadpoleMultiUnitTest::explicitAutoGraphKIsAssemblyK);
		failures+=run("intermediateGraphKIsAccepted", TadpoleMultiUnitTest::intermediateGraphKIsAccepted);
		failures+=run("invalidFuseKIsRejected", TadpoleMultiUnitTest::invalidFuseKIsRejected);
		failures+=run("unusedGraphKIsRejected", TadpoleMultiUnitTest::unusedGraphKIsRejected);
		failures+=run("lowDepthDiagnosticRequestsFinalGraph", TadpoleMultiUnitTest::lowDepthDiagnosticRequestsFinalGraph);
		failures+=run("lowDepthDiagnosticStages", TadpoleMultiUnitTest::lowDepthDiagnosticStages);
		failures+=run("graphClassificationOptions", TadpoleMultiUnitTest::graphClassificationOptions);
		failures+=run("logicalGraphClasses", TadpoleMultiUnitTest::logicalGraphClasses);
		failures+=run("sweepLengthDefaultsAndDisable", TadpoleMultiUnitTest::sweepLengthDefaultsAndDisable);
		failures+=run("omniAnchorPassesToFinalGraph", TadpoleMultiUnitTest::omniAnchorPassesToFinalGraph);
		failures+=run("explicitBlanketControl", TadpoleMultiUnitTest::explicitBlanketControl);
		failures+=run("hashBridgeSelection", TadpoleMultiUnitTest::hashBridgeSelection);
		failures+=run("dispatcherDetectsListsAndPhaseFlags", TadpoleMultiUnitTest::dispatcherDetectsListsAndPhaseFlags);
		System.out.println(failures==0 ? "ALL TESTS PASSED" : failures+" TEST(S) FAILED");
		if(failures>0){System.exit(1);}
	}

	private static int run(final String name, final Test test){
		try{
			test.run();
			System.out.println("PASS: "+name);
			return 0;
		}catch(Throwable t){
			System.err.println("FAIL: "+name);
			t.printStackTrace();
			return 1;
		}
	}

	/** A repeated K represents another pass, not a duplicate to discard. */
	private static void orderedPhasesPreserveOccurrences(){
		final TadpoleMulti.Config c=config("korder=input", "packed=t", "k=96,124,64,32,96");
		check(c.ordered && c.assembleK==96, "Input order must select the first assembly K");
		checkArray(c.phaseKs, 96, 124, 64, 32, 96);
		checkArray(c.fuseKs, 64, 32);
		check(c.graphK==96, "Repeated final96 must supply the graph");
		check(!c.bridgeInitial, "Default initial assembly must not repeat bridge discovery");
		final TadpoleMulti.Config moved=config("korder=input", "packed=t", "k=32,64,96,124", "assemblek=96");
		checkArray(moved.phaseKs, 96, 32, 64, 124);
		check(moved.graphK==124, "Moving the initial K must preserve remaining phase order");
	}

	/** Mixed-order lists expose the difference between assembly-K and last-K defaults. */
	private static void orderedGraphDefaultsAndOverrides(){
		final TadpoleMulti.Config c=config("korder=input", "packed=t", "k=96,124,32,64");
		checkArray(c.phaseKs, 96, 124, 32, 64);
		check(c.graphK==64, "Default graph must reuse final64, not initial96");
		check(config("korder=input", "k=96,124,32,64", "graphk=auto").graphK==64,
				"Explicit auto must also select final64");
		check(config("korder=input", "k=96,124,32,64", "graphk=96").graphK==96,
				"Explicit graphk override was lost");
		final TadpoleMulti.Config normal=new TadpoleMulti.Config(new String[]{"k=96,124,32,64", "out=test.fa"});
		check(normal.ordered && normal.graphK==64, "Ordered processing must be the actual default");
		check(config("korder=input", "k=32,64,96", "assemblek=96").graphK==96,
				"Explicit assembly relocation must not change the user's last listed graph K");
	}

	/** Explicit phase lists filter occurrences; they cannot silently add an unordered K. */
	private static void orderedPhaseFilters(){
		final TadpoleMulti.Config c=config("korder=input", "k=96,124,32,64,32", "fusek=64", "bridgek=124,32");
		checkArray(c.phaseKs, 96, 124, 32, 64, 32);
		checkArray(c.fuseKs, 64);
		checkArray(c.bridgeKs, 124, 32);
		expectFailure("absent from k", "korder=input", "k=96,64", "bridgek=32");
		final TadpoleMulti.Config noList=config("korder=input", "assemblek=96", "fusek=32,64", "bridgek=124,64");
		checkArray(noList.phaseKs, 96, 32, 64, 124);
		check(noList.graphK==124, "Without k, default graph uses the last constructed phase");
		expectFailure("requires a positive", "korder=input", "k=96,64,");
	}

	private static void shorthandExpandsFromLargest(){
		final TadpoleMulti.Config c=config("k=31,63,95,127");
		check(c.assembleK==127, "Wrong assembly k: "+c.assembleK);
		checkArray(c.fuseKs, 95, 63, 31);
		checkArray(c.bridgeKs, 95, 63, 31);
	}

	private static void explicitAssemblyReclassifiesShorthand(){
		final TadpoleMulti.Config c=config("k=63,127,63,31,95", "assemblek=95");
		check(c.assembleK==95, "Wrong assembly k: "+c.assembleK);
		checkArray(c.fuseKs, 63, 31);
		checkArray(c.bridgeKs, 127, 95, 63, 31);
	}

	private static void explicitPhaseListsOverrideShorthand(){
		final TadpoleMulti.Config c=config("k=128,96,64", "assemblek=96",
				"fusek=64", "bridgek=32,128,64,96");
		check(c.assembleK==96, "Wrong assembly k: "+c.assembleK);
		checkArray(c.fuseKs, 64);
		checkArray(c.bridgeKs, 128, 96, 64, 32);
	}

	/** Conflict detection needs a mismatch budget and must stay out of table-loader arguments. */
	private static void fusionConflictOption(){
		check(!config("k=95,63").fuseConflicts, "Conflict restriction changed the default");
		expectFailure("requires fusemaxmismatches", "k=95,63", "fuseconflicts=t");
		final TadpoleMulti.Config c=config("k=95,63", "fuseconflicts=t", "fusemaxmismatches=1");
		check(c.fuseConflicts, "Conflict restriction was not parsed");
		check(!config("k=95,63", "fuseconflicts=t", "fuseconflicts=f").fuseConflicts,
				"Explicit false did not restore legacy selection");
		for(String arg : new TadpoleMulti(c).makeArgs(63, false)){
			check(!arg.startsWith("fuseconflicts="), "Conflict policy leaked into a table loader");
		}
		check(TadpoleMulti.hasMultipleK(new String[]{"k=95", "fuseconflicts=t"}),
				"Conflict option was not routed to the multi-K controller");
	}

	/** Experimental depth ratios are bounded, opt-in, and consumed only by the multi-K controller. */
	private static void fusionCoverageRatioOption(){
		check(config("k=95,63").fuseCoverageRatio==0, "Coverage veto changed the default");
		final TadpoleMulti.Config c=config("k=95,63", "fusecoverageratio=1.75");
		check(c.fuseCoverageRatio==1.75f, "Coverage ratio was not parsed");
		check(config("k=95,63", "fusecoverageratio=1.75", "fusecoverageratio=0").fuseCoverageRatio==0,
				"Explicit zero did not disable the coverage veto");
		for(String value : new String[]{"-1", "0.99", "NaN", "Infinity"}){
			expectFailure("fusecoverageratio", "k=95,63", "fusecoverageratio="+value);
		}
		for(String arg : new TadpoleMulti(c).makeArgs(63, false)){
			check(!arg.startsWith("fusecoverageratio="), "Coverage ratio leaked into a table loader");
		}
		check(TadpoleMulti.hasMultipleK(new String[]{"k=95", "fusecoverageratio=1.75"}),
				"Coverage ratio was not routed to the multi-K controller");
	}

	/** Final-only coverage filtering must not leak into initial assembly or intermediate loaders. */
	private static void finalGraphCoverageRatioOption(){
		check(config("k=95,63").graphMergeCoverageRatio==0, "Final graph coverage filtering changed the default");
		final TadpoleMulti.Config c=config("korder=input", "k=96,124,64,32", "sweeplen=0",
				"graphmergecoverageratio=1.75");
		check(c.graphMergeCoverageRatio==1.75f && c.finalGraphNeeded() && c.graphK==32,
				"Final graph coverage policy did not request the last phase's graph");
		checkArray(c.phaseKs, 96, 124, 64, 32);
		check(c.fuseCoverageRatio==0, "Final-graph policy changed the independent fusion policy");
		for(boolean initial : new boolean[]{false, true}){
			for(String arg : new TadpoleMulti(c).makeArgs(initial ? 96 : 32, initial)){
				check(!arg.startsWith("graphmergecoverageratio="), "Final coverage policy leaked into a loader");
			}
		}
		final TadpoleMulti.Config disabled=config("k=95,63", "sweeplen=0",
				"graphmergecoverageratio=1.75", "graphmergecoverageratio=0");
		check(disabled.graphMergeCoverageRatio==0 && !disabled.finalGraphNeeded(),
				"Explicit zero did not restore the original graph requirements");
		for(String value : new String[]{"-1", "0.99", "NaN", "Infinity"}){
			expectFailure("graphmergecoverageratio", "k=95,63", "graphmergecoverageratio="+value);
		}
		check(TadpoleMulti.hasMultipleK(new String[]{"k=95", "graphmergecoverageratio=1.75"}),
				"Final-graph coverage option was not dispatched to the controller");
	}

	/** Historical eligibility is opt-in, parsed once, and never passed to a table loader. */
	private static void fusionDeadEndOption(){
		check(!config("k=95,63").fuseDeadEndsOnly, "Dead-end restriction changed the default");
		final TadpoleMulti.Config c=config("k=95,63", "fusedeadends=t");
		check(c.fuseDeadEndsOnly, "Dead-end restriction was not parsed");
		check(!config("k=95,63", "fusedeadends=t", "fusedeadends=f").fuseDeadEndsOnly,
				"Explicit false did not restore legacy fusion eligibility");
		for(String arg : new TadpoleMulti(c).makeArgs(63, false)){
			check(!arg.startsWith("fusedeadends="), "Fusion eligibility leaked into a table loader");
		}
		check(TadpoleMulti.hasMultipleK(new String[]{"k=95", "fusedeadends=t"}),
				"Dead-end option was not routed to the multi-K controller");
	}

	/** The guard is opt-in and must not leak into individual table loaders. */
	private static void fusionMismatchOption(){
		check(config("k=95,63").fuseMaxMismatches==-1, "Legacy fusion default changed");
		for(int allowance=0; allowance<=1; allowance++){
			final TadpoleMulti.Config c=config("k=95,63", "fusemaxmismatches="+allowance);
			check(c.fuseMaxMismatches==allowance, "Fusion mismatch budget not parsed");
			for(String arg : new TadpoleMulti(c).makeArgs(63, false)){
				check(!arg.startsWith("fusemaxmismatches="), "Fusion option leaked into a table loader");
			}
		}
		boolean rejected=false;
		try{config("k=95,63", "fusemaxmismatches=-2");}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Invalid negative fusion mismatch budget accepted");
		check(TadpoleMulti.hasMultipleK(new String[]{"k=95", "fusemaxmismatches=1"}),
				"Fusion option was not routed to the multi-K controller");
	}

	private static void singleKHasNoImplicitJoining(){
		final TadpoleMulti.Config c=config("k=96", "graphcover=t");
		check(c.assembleK==96, "Wrong assembly k: "+c.assembleK);
		checkArray(c.fuseKs);
		checkArray(c.bridgeKs);
	}

	private static void noneDisablesPhases(){
		final TadpoleMulti.Config c=config("k=128,96,64", "fusek=none", "bridgek=f");
		checkArray(c.fuseKs);
		checkArray(c.bridgeKs);
	}

	private static void explicitAutoUsesShorthand(){
		final TadpoleMulti.Config c=config("k=64,128,96", "assemblek=auto",
				"fusek=auto", "bridgek=auto");
		check(c.assembleK==128, "Wrong auto assembly k: "+c.assembleK);
		checkArray(c.fuseKs, 96, 64);
		checkArray(c.bridgeKs, 96, 64);
	}

	private static void duplicateKValuesAreDeduplicated(){
		final TadpoleMulti.Config c=config("k=31,63,31,95,63");
		check(c.assembleK==95, "Wrong assembly k: "+c.assembleK);
		checkArray(c.fuseKs, 63, 31);
		checkArray(c.bridgeKs, 63, 31);
	}

	private static void defaultGraphKIsAssemblyK(){
		final TadpoleMulti.Config c=config("k=31,63,95", "simpleomnitigs=t");
		check(c.simpleOmnitigs && !c.graphCover, "Wrong graph mode");
		check(c.graphK==95, "Default graphk was not assembly k: "+c.graphK);
	}

	private static void explicitAutoGraphKIsAssemblyK(){
		final TadpoleMulti.Config c=config("k=31,63,95", "simpleomnitigs=t", "graphk=auto");
		check(c.graphK==95, "Explicit auto graphk was not assembly k: "+c.graphK);
	}

	private static void intermediateGraphKIsAccepted(){
		final TadpoleMulti.Config c=config("k=31,95,127", "graphcover=t", "graphk=63");
		check(c.graphCover && !c.simpleOmnitigs, "Wrong graph-cover mode");
		check(c.graphK==63, "Intermediate graphk changed: "+c.graphK);
	}

	private static void invalidFuseKIsRejected(){
		expectFailure("fusek values must be shorter", "assemblek=95", "fusek=127,95,63");
	}

	private static void unusedGraphKIsRejected(){
		expectFailure("graphk requires", "k=31,63", "graphk=47", "sweeplen=0");
	}

	private static void lowDepthDiagnosticRequestsFinalGraph(){
		final TadpoleMulti.Config c=config("k=31,63,95", "ldcd=t", "graphk=63",
				"ldcmaxlen=700", "ldcmaxcov=4", "ldcfrac=0.15");
		check(c.finalGraphNeeded() && !c.graphOperations(), "Diagnostic did not request a final graph");
		check(c.graphK==63, "Diagnostic graph k changed: "+c.graphK);
		check(c.lowDepthContigMaxLen==700, "Wrong diagnostic maximum length");
		check(c.lowDepthContigMaxCov==4, "Wrong diagnostic absolute depth");
		check(c.lowDepthContigFraction==0.15f, "Wrong diagnostic relative depth");
		for(String s : c.common){
			check(!s.toLowerCase().startsWith("ldc") && !s.toLowerCase().startsWith("lowdepth"),
					"Final-only diagnostic argument leaked into an intermediate phase: "+s);
		}
	}

	private static void lowDepthDiagnosticStages(){
		final TadpoleMulti.Config normal=config("k=31,63,95", "ldcd=t");
		check(!normal.earlyLowDepthDiag() && normal.finalLowDepthDiag(), "Boolean diagnostic did not default to final");
		final TadpoleMulti.Config both=config("k=31,63,95", "ldcd=both");
		check(both.earlyLowDepthDiag() && both.finalLowDepthDiag(), "Both-stage diagnostic was not enabled");
		final TadpoleMulti.Config early=config("k=31,63,95", "ldcdstage=early", "sweeplen=0");
		check(early.earlyLowDepthDiag() && !early.finalLowDepthDiag(), "Early-only diagnostic was not enabled");
		check(!early.finalGraphNeeded(), "Early-only diagnostic unnecessarily requested a final graph");
		check(TadpoleMulti.hasMultipleK(new String[] {"k=95", "ldcd=both", "out=x"}),
				"Stage-valued diagnostic did not select the multi-phase dispatcher");
	}

	private static void graphClassificationOptions(){
		final TadpoleMulti.Config c=config("k=63,95", "graphclassify=t", "graphk=63",
				"gclmc=5", "gclf=0.18", "gcmf=0.45", "gchf=3", "emitsuspect=t", "emitbranchedconnected=t",
				"ecm=3", "evictsuspect=t", "evictmulticonnected=t", "evictdepth=low,semilow", "eca=4",
				"ldce=t", "ldcmaxcov=4", "ldcfrac=0.2", "ldctopology=notloop");
		check(c.finalGraphNeeded(), "Graph classification did not request a final graph");
		check(c.classifyGraphContigs && c.emitTerminal && c.emitBranchedTerminal && c.emitUnanchored && c.emitLoopback,
				"Graph-class emission aliases were not parsed");
		check(c.emitBranchedConnected, "Branched-connected emission was not parsed");
		check(c.evictTerminal && c.evictBranchedTerminal && c.evictUnanchored && c.evictLoopback,
				"Graph-class eviction aliases were not parsed");
		check(c.evictMultiConnected && c.evictGraphDepthMask==3, "Graph-class matrix eviction was not parsed");
		check(c.emitConnectedMax==3 && c.evictConnectedAbove==4, "Graph-class hop thresholds were not parsed");
		check(c.graphClassLowMaxCov==5 && c.graphClassLowFraction==0.18f &&
				c.graphClassMediumFraction==0.45f && c.graphClassHighFraction==3,
				"Graph-class depth boundaries were not parsed");
		check(c.lowDepthContigMaxCov==4 && c.lowDepthContigFraction==0.2f,
				"Graph-class depth thresholds were not parsed");
		check(c.lowDepthContigTopology==1, "Not-loop topology was not parsed");
		for(String s : c.common){
			final String lower=s.toLowerCase();
			check(!lower.startsWith("graphclass") && !lower.startsWith("classify")
					&& !lower.startsWith("emitsuspect") && !lower.startsWith("evictsuspect")
					&& !lower.startsWith("ecm") && !lower.startsWith("eca")
					&& !lower.startsWith("ldce") && !lower.startsWith("evictlowdepth"),
					"Final graph-class argument leaked into an intermediate phase: "+s);
		}
	}

	private static void logicalGraphClasses(){
		final TadpoleMulti.Config terminal=config("k=63,95", "graphclassify=t", "evictgraphclass=terminal");
		check(terminal.evictGraphTopologyMask==
				((1<<Contig.GRAPH_TERMINAL)|(1<<Contig.GRAPH_BRANCHED_TERMINAL)),
				"Logical terminal class was not expanded");
		final TadpoleMulti.Config connected=config("k=63,95", "graphclassify=t", "evictgraphclass=connected");
		check(connected.evictGraphTopologyMask==
				((1<<Contig.GRAPH_CONNECTED)|(1<<Contig.GRAPH_BRANCHED_CONNECTED)|
						(1<<Contig.GRAPH_MULTI_CONNECTED)),
				"Logical connected class was not expanded");
		final TadpoleMulti.Config selfLoop=config("k=63,95", "graphclassify=t", "evictgraphclass=self-loop");
		check(selfLoop.evictGraphTopologyMask==1<<Contig.GRAPH_SELF_LOOP,
				"Logical self-loop class was not parsed");
	}

	private static void sweepLengthDefaultsAndDisable(){
		final TadpoleMulti.Config normal=config("k=63,95");
		check(normal.sweepContigLen==500, "Default sweeplen was not 500");
		check(normal.finalGraphNeeded(), "Default sweep did not request a final graph");
		final TadpoleMulti.Config disabled=config("k=63,95", "sweeplen=0");
		check(disabled.sweepContigLen==0, "sweeplen=0 was not preserved");
		check(!disabled.graphClassificationRequested(), "sweeplen=0 requested graph classification without another graph option");
		final String[] longestArgs=new TadpoleMulti(disabled).makeArgs(disabled.assembleK, true);
		final String[] bridgeArgs=new TadpoleMulti(normal).makeArgs(normal.assembleK, false);
		check(contains(longestArgs, "sweeplen=0"), "Initial multi-k assembly enabled graph sweep");
		check(contains(bridgeArgs, "sweeplen=0"), "Intermediate bridge phase enabled graph sweep");
	}

	private static void omniAnchorPassesToFinalGraph(){
		final TadpoleMulti.Config c=config("k=63,95", "graphcover=t", "omnianchor=124");
		final String[] args=new TadpoleMulti(c).makeArgs(c.assembleK, true);
		check(contains(args, "omnianchor=124"), "Explicit omnianchor was not passed to Tadpole.");
	}

	private static void explicitBlanketControl(){
		final TadpoleMulti.Config c=config("k=63,95", "graphclassify=t", "retainshortcontigs=f");
		check(c.retainShortContigsSet && !c.retainShortContigs, "Explicit blanket control was not preserved");
		final String[] longestArgs=new TadpoleMulti(c).makeArgs(c.assembleK, true);
		int found=0;
		for(String s : longestArgs){if(s.equals("retainshortcontigs=false")){found++;}}
		check(found==1, "Explicit blanket control was overridden or duplicated");
		final TadpoleMulti.Config auto=config("k=63,95", "graphclassify=t");
		final String[] autoArgs=new TadpoleMulti(auto).makeArgs(auto.assembleK, true);
		found=0;
		for(String s : autoArgs){if(s.equals("retainshortcontigs=t")){found++;}}
		check(found==1, "Automatic graph retention was not enabled exactly once");
	}

	private static void hashBridgeSelection(){
		check(config("assemblek=95", "bridgek=127").useHashBridgeTables(),
				"Default long-k bridge did not select compact hashes");
		check(config("assemblek=95", "bridgek=127").bridgeHashMode()==Tadpole.HASH_PAIR,
				"Default long-k bridge did not select a resizable hash pair");
		check(config("assemblek=95", "bridgek=127", "prealloc=t").bridgeHashMode()==Tadpole.HASH_FIXED,
				"Preallocated long-k bridge did not select a fixed fingerprint");
		check(config("assemblek=95", "bridgek=127", "prealloc=t", "hashkmers=pair").bridgeHashMode()==Tadpole.HASH_PAIR,
				"Explicit hashkmers=pair did not override automatic fixed storage");
		check(config("assemblek=95", "bridgek=127", "hashkmers=fixed").bridgeHashMode()==Tadpole.HASH_FIXED,
				"Explicit hashkmers=fixed was not selected");
		check(!config("assemblek=95", "bridgek=127", "hashkmers=f").useHashBridgeTables(),
				"hashkmers=f did not disable compact hashes");
		check(config("assemblek=95", "bridgek=64", "prealloc=t").bridgeHashMode()==Tadpole.HASH_FIXED,
				"Preallocated k=64 bridge did not select a fixed fingerprint");
		check(!config("assemblek=95", "bridgek=64").useHashBridgeTables(),
				"Resizable k=64 bridge selected hash storage without a memory benefit");
		check(config("assemblek=95", "bridgek=127", "wash=t").useHashBridgeTables(),
				"wash=t incorrectly disabled tip-seeded compact washing");
		check(config("assemblek=95", "bridgek=127", "wash=f").useHashBridgeTables(),
				"wash=f incorrectly disabled compact hashes");
		check(!config("assemblek=95", "bridgek=127", "maxcr=10").useHashBridgeTables(),
				"Maximum-count pruning incorrectly selected a non-enumerable table");
		check(config("assemblek=95", "bridgek=127", "mincr=2").useHashBridgeTables(),
				"Minimum-count pruning incorrectly disabled compact hashes");
	}

	private static void dispatcherDetectsListsAndPhaseFlags(){
		check(TadpoleMulti.hasMultipleK(new String[] {"in=x", "K=31,63", "out=y"}),
				"K list was not detected");
		check(TadpoleMulti.hasMultipleK(new String[] {"in=x", "k=95", "bridgek=127,63", "out=y"}),
				"Phase-specific K flag was not detected");
		check(TadpoleMulti.hasMultipleK(new String[] {"in=x", "assemblek=95", "out=y"}),
				"assemblek was not detected");
		check(!TadpoleMulti.hasMultipleK(new String[] {"in=x", "k=63", "extra=a,b", "out=y"}),
				"Non-k comma list triggered multi-k dispatch");
	}

	private static TadpoleMulti.Config config(final String... args){
		final String[] withOut=new String[args.length+2];
		withOut[0]="korder=legacy";// Existing fixtures explicitly retain their historical parser contract.
		System.arraycopy(args, 0, withOut, 1, args.length);
		withOut[args.length+1]="out=test.fa";
		return new TadpoleMulti.Config(withOut);
	}

	private static boolean contains(final String[] array, final String value){
		for(String s : array){if(value.equals(s)){return true;}}
		return false;
	}

	private static void checkArray(final int[] observed, final int... expected){
		check(observed.length==expected.length, "Wrong k count: "+observed.length+" != "+expected.length);
		for(int i=0; i<expected.length; i++){
			check(observed[i]==expected[i], "Wrong k at "+i+": "+observed[i]+" != "+expected[i]);
		}
	}

	private static void expectFailure(final String text, final String... args){
		try{
			config(args);
			throw new AssertionError("Expected failure containing: "+text);
		}catch(RuntimeException e){
			check(e.getMessage()!=null && e.getMessage().contains(text), "Wrong failure: "+e.getMessage());
		}
	}

	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private interface Test {void run();}
}
