package prok;

/** Worker-owned, leave-one-locus-out scores over a finalized positional table.
 * Counts are borrowed read-only. Each set removes one locus's valid word and
 * denominator at EVERY training offset, including repeated words. It caches the
 * union of score positions reachable by jitter within the table's radius, so
 * successive training and same-locus monitoring rows need only array copies.
 * This training view never participates in ordinary genome inference.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaPositionalLeaveOneOut{

	/** Borrows training counts; score-only inference tables cannot exclude a locus. */
	public NcrnaPositionalLeaveOneOut(final RrnaPositionalKmerTable.Table table_){
		this(table_, null);
	}

	/** Optional measured-background composition; correction geometry stays fixed. */
	public NcrnaPositionalLeaveOneOut(final RrnaPositionalKmerTable.Table table_,
			final NcrnaWeightedEndpointTable background){
		if(table_==null || !table_.finished || table_.counts==null || table_.records<2){
			throw new IllegalArgumentException("Positional leave-one-out requires finalized counts from at least two loci");
		}
		if(background!=null && background.geometry!=table_){
			throw new IllegalArgumentException("Optional weighted background must borrow this exact correction table");
		}
		table=table_;
		weighted=background==null ? null : background.leaveOneOut();
		words=new int[table.sites];
		denominators=new double[table.sites];
		positions=new float[Math.addExact(Math.multiplyExact(4, table.radius), 1)];
		missingPrefix=new int[Math.addExact(positions.length, 1)];
	}

	/** Binds a counted locus until the next set; bases are borrowed without mutation.
	 * Ambiguous or clipped words contributed neither a count nor a valid denominator.
	 * A remaining offset with no observations fails, like a literal table rebuild.
	 */
	public void set(final byte[] bases, final int trueEnd){
		ready=false;
		boundBases=null;
		if(bases==null || trueEnd<0 || trueEnd>=bases.length){
			throw new IllegalArgumentException("Leave-one-out needs a real sense-oriented true endpoint");
		}
		if(weighted!=null){weighted.set(bases, trueEnd);}
		else{
			final double prior=table.alpha*table.space;
			for(int site=0; site<table.sites; site++){
				final int word=wordAt(bases, (long)trueEnd+table.offset+site-table.radius);
				words[site]=word;
				if(word>=0 && table.counts[site][word]<1){
					throw new IllegalArgumentException("Excluded locus did not contribute its positional word at offset "+(site-table.radius));
				}
				final long remaining=table.valid[site]-(word>=0 ? 1 : 0);
				if(remaining<=0){
					throw new IllegalArgumentException("Leave-one-out leaves no observed bases at offset "+(site-table.radius));
				}
				denominators[site]=remaining+prior;
			}
		}
		endpoint=trueEnd;
		ready=true;
		try{
			missingPrefix[0]=0;
			for(int position=0; position<positions.length; position++){
				final int word=wordAt(bases, (long)trueEnd+table.offset+position-2L*table.radius);
				missingPrefix[position+1]=missingPrefix[position]+(word<0 ? 1 : 0);
				positions[position]=word<0 ? 0 : finiteScore(word);
			}
			boundBases=bases;
		}catch(RuntimeException | Error failure){
			ready=false;
			throw failure;
		}
	}

	/** Returns the same log ratio as rebuilding the table without the bound locus. */
	public double scoreKey(final int word){
		if(!ready || word<0 || word>=table.space){
			throw new IllegalArgumentException("Bind a counted locus before requesting a valid leave-one-out key");
		}
		if(weighted!=null){return weighted.scoreKey(word);}
		final int anchor=table.radius;
		final double frequency=(table.counts[anchor][word]-(words[anchor]==word ? 1 : 0)+table.alpha)/denominators[anchor];
		double background=0;
		for(int site=0; site<table.sites; site++){
			if(site==anchor){continue;}
			background+=(table.counts[site][word]-(words[site]==word ? 1 : 0)+table.alpha)/denominators[site];
		}
		final double value=Math.log(frequency/(background/(table.sites-1)));
		assert(Double.isFinite(value)) : "Positive smoothed frequencies must yield finite training inputs; word="+word;
		return value;
	}

	/** Copies the score window for a candidate within radius of the bound endpoint.
	 * Candidate coordinates are inclusive sense endpoints, without the table offset.
	 * Returns false without copying if any word is clipped or ambiguous. The
	 * prefix counts make availability constant time without sentinel score values.
	 */
	public boolean scoreWindow(final int candidateEnd, final float[] output, final int from){
		final long delta=(long)candidateEnd-endpoint;
		if(!ready || delta<-table.radius || delta>table.radius || output==null
				|| from<0 || from>output.length-table.sites){
			throw new IllegalArgumentException("Training jitter and output window must fit the bound positional radius");
		}
		final int first=(int)delta+table.radius;
		if(missingPrefix[first+table.sites]!=missingPrefix[first]){return false;}
		System.arraycopy(positions, first, output, from, table.sites);
		return true;
	}

	/** Scores an explicitly requested window at any signed distance from the locus.
	 * In-range candidates reuse the unchanged cached path. Distant candidates use
	 * the same whole-locus-excluded key scores; they do not enlarge the table's
	 * training offsets or correction radius. A missing word returns false; the
	 * incomplete output must not be used. Candidate
	 * span validity belongs to the vector writer, not this single-end word view.
	 * The bases borrowed by set must remain unchanged until the next binding;
	 * even a failed binding releases that reference and invalidates the old view.
	 */
	public boolean scoreDistantWindow(final int candidateEnd, final float[] output, final int from){
		if(!ready || boundBases==null || output==null || from<0 || from>output.length-table.sites){
			throw new IllegalArgumentException("Bind an unchanged counted locus and provide room for its complete score window");
		}
		final long delta=(long)candidateEnd-endpoint;
		if(delta>=-table.radius && delta<=table.radius){
			return scoreWindow(candidateEnd, output, from);
		}
		for(int site=0; site<table.sites; site++){
			final long position=(long)candidateEnd+table.offset+site-table.radius;
			final int word=wordAt(boundBases, position);
			if(word<0){return false;}
			output[from+site]=finiteScore(word);
		}
		return true;
	}

	/** A valid word's computed score cannot be treated as an unavailable observation. */
	private float finiteScore(final int word){
		final float value=(float)scoreKey(word);
		if(!Float.isFinite(value)){throw new IllegalArgumentException("A measured leave-one-out score must fit the finite neural input domain");}
		return value;
	}

	/** Missing and ambiguous positions share the scoring sentinel, never a count. */
	private int wordAt(final byte[] bases, final long position){
		return position<0 || position+table.k>bases.length ? -1
			: RrnaPositionalKmerTable.encode(bases, (int)position, table.k);
	}

	private final RrnaPositionalKmerTable.Table table;
	/** Null keeps the historical uniform-background arithmetic exactly as before. */
	private final NcrnaWeightedEndpointTable.View weighted;
	/** This locus's packed word at each training offset, or -1 when not counted. */
	private final int[] words;
	/** Per-offset remaining valid population plus the unchanged total prior. */
	private final double[] denominators;
	/** Scores spanning true endpoint minus twice the radius through plus twice it. */
	private final float[] positions;
	/** Missing-word prefix counts guard unused cache slots before any array copy. */
	private final int[] missingPrefix;
	private int endpoint;
	private boolean ready;
	/** Read-only borrowed locus, released on every binding attempt; never copied. */
	private byte[] boundBases;
}
