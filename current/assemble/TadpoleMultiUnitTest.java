package assemble;

/** Exact parser and dispatcher regressions for multi-k Tadpole. */
public class TadpoleMultiUnitTest {

	public static void main(String[] args){
		int failures=0;
		failures+=run("ascendingKIsSortedDescending", TadpoleMultiUnitTest::ascendingKIsSortedDescending);
		failures+=run("arbitraryKOrderIsSortedDescending", TadpoleMultiUnitTest::arbitraryKOrderIsSortedDescending);
		failures+=run("defaultGraphKIsShortest", TadpoleMultiUnitTest::defaultGraphKIsShortest);
		failures+=run("explicitAutoGraphKIsShortest", TadpoleMultiUnitTest::explicitAutoGraphKIsShortest);
		failures+=run("intermediateGraphKIsAccepted", TadpoleMultiUnitTest::intermediateGraphKIsAccepted);
		failures+=run("duplicateKIsRejected", TadpoleMultiUnitTest::duplicateKIsRejected);
		failures+=run("outOfRangeGraphKIsRejected", TadpoleMultiUnitTest::outOfRangeGraphKIsRejected);
		failures+=run("unusedGraphKIsRejected", TadpoleMultiUnitTest::unusedGraphKIsRejected);
		failures+=run("dispatcherDetectsOnlyKLists", TadpoleMultiUnitTest::dispatcherDetectsOnlyKLists);
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

	private static void ascendingKIsSortedDescending(){
		final TadpoleMulti.Config c=config("k=31,63,95,127");
		checkKmers(c, 127, 95, 63, 31);
	}

	private static void arbitraryKOrderIsSortedDescending(){
		final TadpoleMulti.Config c=config("k=63,127,31,95");
		checkKmers(c, 127, 95, 63, 31);
	}

	private static void defaultGraphKIsShortest(){
		final TadpoleMulti.Config c=config("k=31,63,95", "simpleomnitigs=t");
		check(c.simpleOmnitigs && !c.graphCover, "Wrong graph mode");
		check(c.graphK==31, "Default graphk was not shortest: "+c.graphK);
	}

	private static void explicitAutoGraphKIsShortest(){
		final TadpoleMulti.Config c=config("k=31,63,95", "simpleomnitigs=t", "graphk=auto");
		check(c.graphK==31, "Explicit auto graphk was not shortest: "+c.graphK);
	}

	private static void intermediateGraphKIsAccepted(){
		final TadpoleMulti.Config c=config("k=31,95,127", "graphcover=t", "graphk=63");
		check(c.graphCover && !c.simpleOmnitigs, "Wrong graph-cover mode");
		check(c.graphK==63, "Intermediate graphk changed: "+c.graphK);
	}

	private static void duplicateKIsRejected(){
		expectFailure("Duplicate kmer length", "k=31,63,31");
	}

	private static void outOfRangeGraphKIsRejected(){
		expectFailure("graphk must be between", "k=31,63,95", "omnitigs=t", "graphk=127");
	}

	private static void unusedGraphKIsRejected(){
		expectFailure("graphk requires", "k=31,63", "graphk=47");
	}

	private static void dispatcherDetectsOnlyKLists(){
		check(TadpoleMulti.hasMultipleK(new String[] {"in=x", "K=31,63", "out=y"}),
				"K list was not detected");
		check(!TadpoleMulti.hasMultipleK(new String[] {"in=x", "k=63", "extra=a,b", "out=y"}),
				"Non-k comma list triggered multi-k dispatch");
	}

	private static TadpoleMulti.Config config(final String... args){
		final String[] withOut=new String[args.length+1];
		System.arraycopy(args, 0, withOut, 0, args.length);
		withOut[args.length]="out=test.fa";
		return new TadpoleMulti.Config(withOut);
	}

	private static void checkKmers(final TadpoleMulti.Config c, final int... expected){
		check(c.kmers.length==expected.length, "Wrong k count");
		for(int i=0; i<expected.length; i++){
			check(c.kmers[i]==expected[i], "Wrong k at "+i+": "+c.kmers[i]+" != "+expected[i]);
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
