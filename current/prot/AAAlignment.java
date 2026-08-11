package prot;

/**
 * Result of a single local protein alignment (one HSP).
 *
 * <p>All coordinates are 0-based and inclusive, referring to positions in the
 * original query and target residue arrays. The BLAST-tab emitter converts them
 * to 1-based inclusive on output.</p>
 *
 * @author Eru
 */
public final class AAAlignment {

	/** Raw BLOSUM62 score of the alignment (sum of column scores). */
	public final int rawScore;
	/** First aligned query residue (0-based, inclusive). */
	public final int qStart;
	/** Last aligned query residue (0-based, inclusive). */
	public final int qStop;
	/** First aligned target residue (0-based, inclusive). */
	public final int tStart;
	/** Last aligned target residue (0-based, inclusive). */
	public final int tStop;
	/** Number of identical aligned columns (standard residues only; X never counts). */
	public final int identities;
	/** Number of aligned non-gap columns that are not identities. */
	public final int mismatches;
	/** Number of gap openings (contiguous gap runs) in either sequence. */
	public final int gapOpens;
	/** Alignment length: total columns, including gap columns. */
	public final int length;
	/**
	 * Per-column operations, left-to-right over the aligned region: 'm' = match/substitution
	 * (consumes a query and a reference residue), 'D' = deletion (reference residue vs a gap in
	 * the query), 'I' = insertion (query residue vs a gap in the reference). Null unless the
	 * alignment was requested with path recording ({@link AAAligner#align(byte[],byte[],boolean)});
	 * this is what an alignment/consensus graph needs to add a member.
	 */
	public final byte[] match;

	/**
	 * Constructs an immutable alignment result.
	 * @param rawScore Raw alignment score.
	 * @param qStart First aligned query residue (0-based).
	 * @param qStop Last aligned query residue (0-based).
	 * @param tStart First aligned target residue (0-based).
	 * @param tStop Last aligned target residue (0-based).
	 * @param identities Identical column count.
	 * @param mismatches Non-identical non-gap column count.
	 * @param gapOpens Gap-run count.
	 * @param length Total alignment columns.
	 * @param match Per-column m/D/I op string, or null if path recording was not requested.
	 */
	public AAAlignment(int rawScore, int qStart, int qStop, int tStart, int tStop,
			int identities, int mismatches, int gapOpens, int length, byte[] match){
		this.rawScore=rawScore;
		this.qStart=qStart;
		this.qStop=qStop;
		this.tStart=tStart;
		this.tStop=tStop;
		this.identities=identities;
		this.mismatches=mismatches;
		this.gapOpens=gapOpens;
		this.length=length;
		this.match=match;
	}

	/**
	 * Percent identity = 100 * identities / alignment length (gap columns included).
	 * @return Percent identity in [0,100].
	 */
	public final double pident(){return length==0 ? 0 : (100.0*identities)/length;}

	/**
	 * Bitscore from the raw score using gapped-BLOSUM62 Karlin-Altschul params.
	 * @return Normalized bitscore.
	 */
	public final double bitscore(){
		return (Blosum62.LAMBDA*rawScore-Blosum62.LN_K)/Blosum62.LN_2;
	}

	/**
	 * Karlin-Altschul E-value for a given effective search space. This uses the
	 * raw search space (query length times total database residues) WITHOUT the
	 * edge-effect length correction, so it is an approximation of a rigorous
	 * BLAST E-value (documented as approximate; see ProteinSearcher).
	 *
	 * @param searchSpace Effective search space (residue-pairs).
	 * @return Expected number of alignments this good or better by chance.
	 */
	public final double evalue(double searchSpace){
		return searchSpace*Blosum62.K*Math.exp(-Blosum62.LAMBDA*rawScore);
	}
}
