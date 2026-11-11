package assemble;

import dna.AminoAcid;
import ukmer.Kmer;

/**
 * Checks a same-K path through the retained flanks and overlap of an exact fusion.
 * Both directions must select the proposed continuation without a significant
 * branch. This checks graph consistency, not support from a single spanning read.
 * Reused rolling keys avoid allocating a sequence or a record for each join.
 * The caller keeps the borrowed count provider live and stable throughout each
 * serial phase. Scratch keys and counters make an instance unsafe for concurrent
 * calls; this class does not own or release the underlying table.
 * @author Fischl
 */
public final class FusionKmerSupport{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrows a stable count provider for serial fusion checks, retaining one
	 * effective K of context on each side. See the four-argument constructor. */
	public FusionKmerSupport(final Evidence counts_, final int k_, final int minDepth_){
		this(counts_, k_, minDepth_, 0);
	}

	/** Creates reusable oriented keys using the provider's hollow/explicit layout.
	 * A zero minimum preserves K-base flanks; a positive value can only extend them.
	 * @param counts_ Borrowed count provider and existing branch policy
	 * @param k_ Requested word length, at least two; effective K is the key's kbig
	 * @param minDepth_ Positive minimum depth for each word and intended continuation
	 * @param minFlank_ Nonnegative minimum retained bases per flank; effective
	 * context is max(effective K, minFlank_)
	 * @throws IllegalArgumentException If the provider is null or a range is invalid
	 */
	public FusionKmerSupport(final Evidence counts_, final int k_, final int minDepth_, final int minFlank_){
		if(counts_==null || k_<2 || minDepth_<1 || minFlank_<0){
			throw new IllegalArgumentException("Fusion support needs counts, K>=2, depth>=1 and flank>=0.");
		}
		counts=counts_;
		key=new Kmer(k_, counts.hollow());
		neighbor=new Kmer(k_, counts.hollow());
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
	 * The caller has already established exact-overlap geometry. This method neither
	 * compares the two overlap strings nor merges or modifies either contig.
	 * A normal rejection stops the scan immediately; counters include that scanned
	 * prefix, not all words that would have been visited on success.
	 * @param source Contig contributing the retained left flank and overlap
	 * @param sourceReverse Whether to read source in reverse-complement orientation
	 * @param sourceTrim Bases discarded from the oriented source's right end
	 * @param dest Contig contributing the retained right flank
	 * @param destReverse Whether to read destination in reverse-complement orientation
	 * @param destTrim Bases discarded from the oriented destination's left end
	 * @param overlap Positive exact-overlap length in bases
	 * @return True when the whole window has sufficient defined-base depth and
	 * every internal step dominates alternatives under the borrowed branch policy
	 */
	public boolean supported(final Contig source, final boolean sourceReverse, final int sourceTrim,
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

	/** Tests all four neighbors of the current complete key in one direction.
	 * The intended step must meet minimum depth, strictly dominate alternatives
	 * and pass the provider's branch rule. Absent/negative depths are treated as zero.
	 * Reuses neighbor scratch without changing key.
	 * @param expectedBase Intended continuation in the oriented virtual product
	 * @param right True for an outgoing step, false for an incoming step
	 */
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

	/** Reads a window-relative position from two oriented slices without copying.
	 * Positions below split use source[start+pos]; the rest use
	 * dest[destStart+pos-split]. The source slice already includes the overlap. */
	private static byte productBase(final Contig source, final boolean sr, final int start,
			final Contig dest, final boolean dr, final int destStart, final int split, final int pos){
		assert(pos>=0) : "Fusion window offsets must index retained sequence, not discarded tails.";
		return pos<split ? baseAt(source, sr, start+pos) : baseAt(dest, dr, destStart+pos-split);
	}

	/** Borrows one live table and its existing error-aware branch policy. */
	public interface Evidence extends TadpoleGraph.Counts{
		/** True when the alternative is too strong to treat as a sequencing-error branch.
		 * @param max Depth of the intended continuation
		 * @param second Largest competing depth; the caller requires max greater than second
		 */
		boolean isJunction(int max, int second);
	}

	/** Reads a zero-based position in the requested orientation, complementing when
	 * reversed. The contig and its cached metadata are unchanged. */
	private static byte baseAt(final Contig c, final boolean reverse, final int pos){
		assert(pos>=0 && pos<c.length()) : "Fusion window must stay inside retained contig sequence: "+pos;
		return reverse ? AminoAcid.baseToComplementExtended[c.bases[c.length()-1-pos]] : c.bases[pos];
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed counts and branch policy; lifetime and stability are caller-owned. */
	private final Evidence counts;
	/** Rolling virtual-product word and independent scratch for its four neighbors. */
	private final Kmer key, neighbor;
	/** Effective word length and minimum accepted word/continuation depth. */
	final int k, minDepth;
	/** Retained bases per flank, at least the effective word length. */
	public final int flank;
	/** Cumulative calls and normal false returns; exceptions are not normal rejections. */
	long evaluations, rejected;
	/** Complete words visited, including the word whose depth or transition rejects a call. */
	public long words;
}
