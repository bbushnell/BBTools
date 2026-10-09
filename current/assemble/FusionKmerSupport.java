package assemble;

import dna.AminoAcid;
import ukmer.Kmer;

/**
 * Checks a same-K path through the retained flanks and overlap of an exact fusion.
 * Both directions must select the proposed continuation without a significant
 * branch. This checks graph consistency, not support from a single spanning read.
 * Reused rolling keys avoid allocating a sequence or a record for each join.
 * @author Fischl
 */
final class FusionKmerSupport {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Retains immutable counts; callers use the checker only in serial fusion phases. */
	FusionKmerSupport(final Evidence counts_, final int k_, final int minDepth_){
		this(counts_, k_, minDepth_, 0);
	}

	/** A zero minimum preserves K-base flanks; a positive value can only extend them. */
	FusionKmerSupport(final Evidence counts_, final int k_, final int minDepth_, final int minFlank_){
		if(counts_==null || k_<2 || minDepth_<1 || minFlank_<0){
			throw new IllegalArgumentException("Fusion support needs counts, K>=2, depth>=1 and flank>=0.");
		}
		counts=counts_;
		key=new Kmer(k_);
		neighbor=new Kmer(k_);
		k=key.kbig;
		minDepth=minDepth_;
		flank=Math.max(k, minFlank_);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Checks the full overlap plus F retained bases on each flank, where F>=K.
	 * There are overlap+2F-K+1 words (65 for overlap=K=F=32). Insufficient flank
	 * context rejects rather than silently weakening the test. Each internal
	 * transition is checked forward and backward; outward directions beyond the
	 * bounded window are irrelevant. Discarded tails never enter the product.
	 */
	boolean supported(final Contig source, final boolean sourceReverse, final int sourceTrim,
			final Contig dest, final boolean destReverse, final int destTrim, final int overlap){
		assert(source!=null && dest!=null && sourceTrim>=0 && destTrim>=0 && overlap>0) :
				"Cross-k support must receive the positive exact-overlap geometry before a merge.";
		evaluations++;
		final int sourceEnd=source.length()-sourceTrim, destStart=destTrim+overlap;
		if(sourceEnd-overlap<flank || dest.length()-destStart<flank){
			rejected++;
			return false;
		}
		final int start=sourceEnd-overlap-flank, length=Math.toIntExact(overlap+2L*flank);
		key.clear();
		for(int i=0; i<length; i++){
			final byte base=productBase(source, sourceReverse, start, dest, destReverse, destStart, overlap+flank, i);
			if(AminoAcid.baseToNumber[base]<0){rejected++; return false;}
			key.addRight(base);
			if(i<k-1){continue;}
			words++;
			if(counts.count(key)<minDepth){
				rejected++;
				return false;
			}
			if(i>=k && !unbranched(productBase(source, sourceReverse, start, dest,
					destReverse, destStart, overlap+flank, i-k), false)){
				rejected++; return false;
			}
			if(i+1<length && !unbranched(productBase(source, sourceReverse, start, dest,
					destReverse, destStart, overlap+flank, i+1), true)){
				rejected++; return false;
			}
		}
		return true;
	}

	/** The intended step must dominate alternatives as well as pass the branch-ratio rule. */
	private boolean unbranched(final byte expectedBase, final boolean right){
		final int expected=AminoAcid.baseToNumber[expectedBase];
		if(expected<0){return false;}
		assert(key.len>=k) : "Neighbor queries require a complete word from the retained fusion window.";
		int intended=0, alternative=0;
		for(int b=0; b<4; b++){
			neighbor.setFrom(key);
			if(right){neighbor.addRightNumeric(b);}
			else{neighbor.addLeftNumeric(b);}
			final int depth=Math.max(0, counts.count(neighbor));
			if(b==expected){intended=depth;}
			else{alternative=Math.max(alternative, depth);}
		}
		return intended>=minDepth && intended>alternative && !counts.isJunction(intended, alternative);
	}

	/** Reads the virtual product without copying the contigs or including trimmed bases. */
	private static byte productBase(final Contig source, final boolean sr, final int start,
			final Contig dest, final boolean dr, final int destStart, final int split, final int pos){
		assert(pos>=0) : "Fusion window offsets must index retained sequence, not discarded tails.";
		return pos<split ? baseAt(source, sr, start+pos) : baseAt(dest, dr, destStart+pos-split);
	}

	/** Borrows one live table and its existing error-aware branch policy. */
	interface Evidence extends TadpoleGraph.Counts {
		/** True when the alternative is too strong to treat as a sequencing-error branch. */
		boolean isJunction(int max, int second);
	}

	/** Reads an oriented base without modifying either contig. */
	private static byte baseAt(final Contig c, final boolean reverse, final int pos){
		assert(pos>=0 && pos<c.length()) : "Fusion window must stay inside retained contig sequence: "+pos;
		return reverse ? AminoAcid.baseToComplementExtended[c.bases[c.length()-1-pos]] : c.bases[pos];
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final Evidence counts;
	private final Kmer key, neighbor;
	final int k, minDepth, flank;
	long evaluations, rejected, words;
}
