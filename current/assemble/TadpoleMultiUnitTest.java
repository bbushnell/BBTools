package assemble;

/** Exact parser and dispatcher regressions for multi-k Tadpole. */
public class TadpoleMultiUnitTest {

	public static void main(String[] args){
		int failures=0;
		failures+=run("shorthandExpandsFromLargest", TadpoleMultiUnitTest::shorthandExpandsFromLargest);
		failures+=run("explicitAssemblyReclassifiesShorthand", TadpoleMultiUnitTest::explicitAssemblyReclassifiesShorthand);
		failures+=run("explicitPhaseListsOverrideShorthand", TadpoleMultiUnitTest::explicitPhaseListsOverrideShorthand);
		failures+=run("singleKHasNoImplicitJoining", TadpoleMultiUnitTest::singleKHasNoImplicitJoining);
		failures+=run("noneDisablesPhases", TadpoleMultiUnitTest::noneDisablesPhases);
		failures+=run("explicitAutoUsesShorthand", TadpoleMultiUnitTest::explicitAutoUsesShorthand);
		failures+=run("duplicateKValuesAreDeduplicated", TadpoleMultiUnitTest::duplicateKValuesAreDeduplicated);
		failures+=run("defaultGraphKIsAssemblyK", TadpoleMultiUnitTest::defaultGraphKIsAssemblyK);
		failures+=run("explicitAutoGraphKIsAssemblyK", TadpoleMultiUnitTest::explicitAutoGraphKIsAssemblyK);
		failures+=run("intermediateGraphKIsAccepted", TadpoleMultiUnitTest::intermediateGraphKIsAccepted);
		failures+=run("invalidFuseKIsRejected", TadpoleMultiUnitTest::invalidFuseKIsRejected);
		failures+=run("unusedGraphKIsRejected", TadpoleMultiUnitTest::unusedGraphKIsRejected);
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
		expectFailure("graphk requires", "k=31,63", "graphk=47");
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
		final String[] withOut=new String[args.length+1];
		System.arraycopy(args, 0, withOut, 0, args.length);
		withOut[args.length]="out=test.fa";
		return new TadpoleMulti.Config(withOut);
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
