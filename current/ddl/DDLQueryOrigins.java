package ddl;

import java.util.Arrays;

/** Reused residue origins for the final winning value of each query register.
 * Equal winning values may occur repeatedly; forward and reverse scans use
 * their earliest and latest occurrences respectively. No state is allocated
 * by ordinary sketch construction.
 * One worker owns this mutable scratch. Reset before observing a fresh sketch;
 * all returned/public arrays are borrowed and overwritten by subsequent use.
 * @author Fischl
 */
public final class DDLQueryOrigins{

	/** Allocate fixed register scratch; position buckets grow only for longer queries. */
	public DDLQueryOrigins(int width){
		if(width<1){throw new IllegalArgumentException("Origin scratch requires positive width");}
		first=new int[width];
		last=new int[width];
		bins=new int[width];
		next=new int[width];
		tags=new int[width];
	}

	/** Resets origin/tag state for a fresh zeroed sketch of the same width.
	 * Ordering output beyond size and old quadrant boundaries are not cleared here.
	 */
	public void reset(int width){
		if(width!=first.length){throw new IllegalArgumentException("Origin/sketch widths differ");}
		Arrays.fill(first, -1);
		Arrays.fill(last, -1);
		Arrays.fill(tags, 0);
		size=0;
		ambiguous=0;
	}

	/** Records a candidate before replacing its register; lower scores leave origins unchanged.
	 * Ties retain both extreme zero-based residue starts and the union of position tags.
	 */
	public void observe(int bin, int score, int oldScore, int start, int tag){
		assert(bin>=0 && bin<first.length && score>0 && start>=0 && tag>=0 && tag<16) :
			"Origins describe positive sketch values at valid K5 starts";
		if(score>oldScore){first[bin]=last[bin]=start; tags[bin]=1<<tag;}
		else if(score==oldScore){
			assert(first[bin]>=0) : "An equal positive value must already have a recorded winning origin";
			first[bin]=Math.min(first[bin], start);
			last[bin]=Math.max(last[bin], start);
			tags[bin]|=1<<tag;
		}
	}

	/** Counting-sort occupied bins by residue coordinate, then ascending bin on ties.
	 * Returns the reusable bin array; only the first size entries are populated.
	 * Forward uses earliest origins ascending; reverse uses latest origins descending.
	 */
	public int[] order(char[] sketch, int sequenceLength, boolean reverse){
		if(sketch==null || sketch.length!=first.length || sequenceLength<0){
			throw new IllegalArgumentException("Origin order requires matching sketch and sequence length");
		}
		if(heads.length<sequenceLength){heads=new int[sequenceLength];}
		Arrays.fill(heads, 0, sequenceLength, -1);
		size=0;
		for(int bin=sketch.length-1; bin>=0; bin--){
			if(sketch[bin]==0){continue;}
			final int position=reverse ? last[bin] : first[bin];
			if(position<0 || position>=sequenceLength){throw new IllegalStateException("Winning key lacks a valid residue origin");}
			next[bin]=heads[position];
			heads[position]=bin;
		}
		for(int step=0; step<sequenceLength; step++){
			final int position=reverse ? sequenceLength-1-step : step;
			for(int bin=heads[position]; bin>=0; bin=next[bin]){bins[size++]=bin;}
		}
		assert(size<=sketch.length) : "Each occupied register enters exactly one coordinate bucket";
		return bins;
	}

	/** Four tagged groups, with one visit per selected key even if encoded values alias.
	 * Checkpoints occur after each quadrant, rather than after equal key counts.
	 * stopFirst selects the lowest present tag, otherwise the highest; bins within a
	 * quadrant are ascending. Updates size, ambiguous and the four keysAt boundaries.
	 * @return Borrowed bin array with size valid entries
	 */
	public int[] orderQuadrants(char[] sketch, boolean stopFirst){
		if(sketch==null || sketch.length!=tags.length){throw new IllegalArgumentException("Quadrant/sketch widths differ");}
		size=0;
		ambiguous=0;
		for(int group=0; group<4; group++){
			final int quadrant=stopFirst ? group : 3-group;
			for(int bin=0; bin<sketch.length; bin++){
				if(sketch[bin]==0){continue;}
				final int mask=tags[bin];
				if(mask<1 || mask>15){throw new IllegalArgumentException("Quadrant ordering requires four tagged partitions");}
				final int selected=stopFirst ? Integer.numberOfTrailingZeros(mask) : 31-Integer.numberOfLeadingZeros(mask);
				if(selected!=quadrant){continue;}
				bins[size++]=bin;
				if(Integer.bitCount(mask)>1){ambiguous++;}
			}
			keysAt[group]=size;
		}
		assert(size<=sketch.length) : "An encoded key is assigned to its first eligible quadrant exactly once";
		return bins;
	}

	/** Earliest/latest zero-based residue starts, output bin order and winning-tag masks. */
	public final int[] first, last, bins, tags;
	/** Cumulative occupied-key counts after each group in the last quadrant ordering. */
	public final int[] keysAt=new int[4];
	/** Valid output entries and multi-tag winners counted by the last quadrant ordering. */
	public int size, ambiguous;
	/** Per-bin links within the coordinate counting-sort buckets. */
	private final int[] next;
	/** First bin at each residue coordinate; grows to the longest sequence seen. */
	private int[] heads=new int[0];
}
