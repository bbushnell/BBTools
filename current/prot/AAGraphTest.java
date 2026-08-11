package prot;

/**
 * Synthetic tests for {@link AAGraph}: verifies (1) the consensus grows into the
 * padding when members share a terminal extension the pivot lacks, and (2) raising
 * the ANI (identity-inverted) ceiling pulls the consensus toward a low-identity
 * minority that count-majority voting would otherwise bury.
 *
 * <p>Run: {@code java -cp current --add-modules jdk.incubator.vector prot.AAGraphTest}</p>
 *
 * @author Eru
 */
public final class AAGraphTest {

	private static int passed=0, failed=0;

	public static void main(String[] args){
		testGrowthIntoPadding();
		testAniCeilingPullsDivergent();
		testLeadingOverhangSafe();
		testDepthTrimRemovesSmearedTail();
		System.err.println("\n"+passed+" passed, "+failed+" failed.");
		if(failed>0){System.exit(1);}
	}

	/** Members that all extend the pivot on the C-terminus should grow the consensus. */
	private static void testGrowthIntoPadding(){
		final byte[] pivot=enc("ACDEFGHIKL");
		final AAGraph g=new AAGraph(pivot, 6);
		final byte[] member=enc("ACDEFGHIKLMNPQ");//pivot + C-terminal MNPQ
		for(int i=0; i<3; i++){addMember(g, member);}

		final String cons=g.consensusString();
		check("growth: consensus contains the shared extension MNPQ", cons.contains("MNPQ"));
		check("growth: consensus grew beyond the pivot length", cons.length()>10);
		System.err.println("  pivot=ACDEFGHIKL  consensus="+cons);
	}

	/** A low-identity minority is buried by count, but surfaces once ANI weighting is on. */
	private static void testAniCeilingPullsDivergent(){
		//Pivot == the count-majority sequence; a 2-member minority diverges at cols 1-8.
		final String pivotStr="ACDEFGHIKL";
		final byte[] pivot=enc(pivotStr);
		final byte[] majority=enc(pivotStr);          //identical, 100% identity
		final byte[] divergent=enc("AWWWWWWWWL");      //A..L shared, W elsewhere -> ~20% identity

		//OFF: plain count majority -> consensus should track the majority (no W columns).
		final AAGraph off=new AAGraph(pivot, 4);
		for(int i=0; i<3; i++){addMember(off, majority);}
		for(int i=0; i<2; i++){addMember(off, divergent);}
		final String consOff=off.consensusString();
		final int wOff=countChar(consOff, 'W');

		//ON: identity-inverted weighting with a high ceiling -> minority is up-weighted.
		final AAGraph on=new AAGraph(pivot, 4);
		on.weightByIdentity=true;
		on.identityCeiling=100f;
		for(int i=0; i<3; i++){addMember(on, majority);}
		for(int i=0; i<2; i++){addMember(on, divergent);}
		final String consOn=on.consensusString();
		final int wOn=countChar(consOn, 'W');

		System.err.println("  ANI off: "+consOff+"  (W="+wOff+")");
		System.err.println("  ANI on : "+consOn+"  (W="+wOn+")");
		check("ANI: high ceiling pulls consensus toward the divergent minority", wOn>wOff);
		check("ANI: without weighting the majority wins", wOff==0);
	}

	/**
	 * A member whose N-terminal overhang exceeds the padding forces a leading insertion
	 * (no node to attach to). This must not crash (regression for the glocalTraceback qStart
	 * bug + the null-prevNode path); the consensus should still recover the shared core.
	 */
	private static void testLeadingOverhangSafe(){
		final byte[] pivot=enc("ACDEFGHIKLMNPQRST");
		final AAGraph g=new AAGraph(pivot, 1);//pad=1: overhang can't all fit as padding matches
		final byte[] member=enc("WWWWWWWWWWACDEFGHIKLMNPQRST");//10-residue N-terminal extension
		boolean crashed=false;
		try{
			for(int i=0; i<3; i++){addMember(g, member);}
		}catch(Throwable e){crashed=true; System.err.println("  (exception: "+e+")");}
		final String cons=(crashed ? "" : g.consensusString());
		check("leading overhang beyond pad does not crash", !crashed);
		check("leading overhang: consensus still recovers the shared core",
				cons.contains("ACDEFGHIKLMNPQRST"));
		System.err.println("  pad=1 overhang consensus="+cons);
	}

	/**
	 * A weakly-supported extension (few members reach into the padding) must not leave a
	 * "core, then padding X, then a couple residues" tail: depth-trimming strips the whole
	 * low-depth padding tail, leaving a clean core. Regression for the XXXX-tail bug.
	 */
	private static void testDepthTrimRemovesSmearedTail(){
		final byte[] pivot=enc("ACDEFGHIKLMNPQRST");
		final AAGraph g=new AAGraph(pivot, 6);
		g.trimDepthFraction=0.1f;
		for(int i=0; i<30; i++){addMember(g, pivot);}         //30 core members -> high-depth core
		addMember(g, enc("ACDEFGHIKLMNPQRSTWYV"));            //1 sparse low-depth extension
		final String cons=g.consensusString();
		check("depth-trim: sparse low-depth padding tail removed", cons.equals("ACDEFGHIKLMNPQRST"));
		check("depth-trim: no padding X survives in the consensus", cons.indexOf('X')<0);
		System.err.println("  depth-trim consensus="+cons);
	}

	/*--------------------------------------------------------------*/

	/** Encodes a member and adds it to the graph via a glocal alignment to the pivot. */
	private static void addMember(AAGraph g, byte[] memberEnc){
		final AAAlignment aln=AAAligner.alignGlocal(memberEnc, g.pivot, true);
		g.add(memberEnc, aln);
	}

	private static byte[] enc(String s){return Blosum62.encode(s.getBytes(), s);}

	private static int countChar(String s, char c){
		int n=0;
		for(int i=0; i<s.length(); i++){if(s.charAt(i)==c){n++;}}
		return n;
	}

	private static void check(String name, boolean ok){
		System.err.println((ok ? "PASS" : "FAIL")+": "+name);
		if(ok){passed++;}else{failed++;}
	}
}
