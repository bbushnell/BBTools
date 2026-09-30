package align2;

/**
 * Allocation-free, pure-Java score of a gapless m/S/N operation sequence.
 * Uses the unweighted MSA operation-score contract, including run-transition
 * adjustments; this is not the historical scoreNoIndels fast-score contract.
 * MultiStateAligner11ts uses it only when the experimental
 * bbmap3.nativeGaplessScore system property is enabled at aligner construction.
 * The name reflects that experiment, not a JNI dependency in this class.
 * No mapper accuracy or performance qualification is implied.
 * @author Collei
 */
public final class NativeGaplessScore{
	private NativeGaplessScore(){}

	/**
	 * Scores query against consecutive reference positions without allocating a trace.
	 * Literal uppercase N in either sequence, and positions outside reference,
	 * produce N operations. All other bytes use literal equality, including lowercase
	 * and non-ACGT symbols; this is not an IUPAC compatibility test. The operation
	 * classification follows MultiStateAligner11ts.genMatchNoIndels, with widened
	 * coordinate arithmetic so start+i cannot wrap into the reference.
	 * Null sequences or an empty query return zero. A scoring model is required
	 * for nonempty input; no model state or sequence content is modified.
	 */
	public static int score(MSA msa, byte[] query, byte[] reference, int start){
		if(query==null || reference==null || query.length==0){return 0;}
		byte mode=0, previous=0;
		int length=0, previousLength=0, total=0;
		for(int i=0; i<query.length; i++){
			long position=(long)start+i;
			byte r=position<0 || position>=reference.length ? (byte)'N' : reference[(int)position];
			byte q=query[i], op=q=='N' || r=='N' ? (byte)'N' : q==r ? (byte)'m' : (byte)'S';
			if(op==mode){length++;}
			else{
				if(length>0){total+=run(msa, mode, length, previous, previousLength);}
				previous=mode;
				previousLength=length;
				mode=op;
				length=1;
			}
		}
		return total+run(msa, mode, length, previous, previousLength);
	}

	/**
	 * Scores one nonempty run. As in MSA.score(byte[]), a substitution after N
	 * receives the continued-substitution cost, while one after a single match
	 * receives POINTS_SUBR. Leading substitutions use the ordinary initial cost.
	 * Only m/S/N are emitted by score, so insertion, deletion and R runs are absent.
	 */
	private static int run(MSA msa, byte mode, int length, byte previous, int previousLength){
		if(mode=='m'){return msa.calcMatchScore(length);}
		if(mode=='N'){return msa.calcNocallScore(length);}
		assert(mode=='S');
		int score=msa.calcSubScore(length);
		if(previous=='N'){score+=msa.POINTS_SUB2()-msa.POINTS_SUB();}
		else if(previous=='m' && previousLength==1){score+=msa.POINTS_SUBR()-msa.POINTS_SUB();}
		return score;
	}
}
