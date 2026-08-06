package prot;

/**
 * The copy-number distribution of one protein family across a set of genomes:
 * how many genomes carry it exactly 0, 1, 2, 3, or 4-or-more times.
 *
 * <p>This is the shape the freeze (§4c) requires — the <b>full distribution</b>,
 * not a mean, because the fraction of genomes carrying a family <b>exactly
 * once</b> is what identifies a single-copy marker and a mean copy count hides
 * it.</p>
 *
 * <p>The bin edges (0, 1, 2, 3, &ge;4) are the Phase-3 default; the freeze marks
 * the exact edges as a later decision, so they are collected here in one place.</p>
 *
 * @author Eru
 */
public final class CopyNumberDistribution {

	/** Number of bins: exactly 0, 1, 2, 3, and 4-or-more copies. */
	public static final int BINS=5;
	/** Index of the terminal "4 or more copies" bin. */
	public static final int PLUS_BIN=4;

	/** bins[c] = genomes carrying this family exactly c times (bins[4] = 4+). */
	public final int[] bins=new int[BINS];

	/**
	 * Records one genome that carries this family {@code copies} times.
	 * @param copies Number of copies in that genome (&ge;0; folded into the 4+ bin).
	 */
	public final void add(final int copies){
		if(copies<0){throw new RuntimeException("Negative copy count: "+copies);}
		bins[Math.min(copies, PLUS_BIN)]++;
	}

	/** Total genomes tallied (the prevalence/fraction denominator). @return Genome count. */
	public final int denominator(){
		int sum=0;
		for(final int b : bins){sum+=b;}
		return sum;
	}

	/** Genomes carrying the family at least once. @return Present-genome count. */
	public final int present(){return denominator()-bins[0];}

	/**
	 * Fraction of genomes carrying the family one or more times.
	 * @return Prevalence in [0,1], or 0 when no genomes were tallied.
	 */
	public final double prevalence(){
		final int d=denominator();
		return d==0 ? 0.0 : present()/(double)d;
	}

	/**
	 * Fraction of genomes carrying the family exactly once — the single-copy shape.
	 * @return Fraction in [0,1], or 0 when no genomes were tallied.
	 */
	public final double fractionExactlyOnce(){
		final int d=denominator();
		return d==0 ? 0.0 : bins[1]/(double)d;
	}
}
