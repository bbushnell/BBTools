package assemble;

import dna.AminoAcid;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/** Rolling original-read support for a saved long-kmer context and its next choice.
 * This qualifies a separately established contradiction; it never creates a veto.
 * Counts are immutable, and neither the saved query nor its seed is added here.
 * @author Brian Bushnell, Fischl */
public final class ShortKmerPathSupport {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Wraps an already loaded, explicit Solid64 or Solid96 table. */
	public ShortKmerPathSupport(final KmerTableSetU tables, final int contextSpan){
		this(tableLookup(tables), contextSpan, tables.kbig);
	}

	/** Preserves the original Solid64 lookup geometry. */
	public ShortKmerPathSupport(final CountLookup lookup_, final int contextSpan_){
		this(lookup_, contextSpan_, 64);
	}

	/** An explicit immutable lookup permits independent rolling-state fixtures.
	 * The selected exact word must fit entirely within the long context. */
	public ShortKmerPathSupport(final CountLookup lookup_, final int contextSpan_, final int supportK_){
		if(lookup_==null || (supportK_!=64 && supportK_!=96) || contextSpan_<supportK_){
			throw new IllegalArgumentException("Short support requires a lookup, k=64 or96, and a context containing that complete word.");
		}
		lookup=lookup_;
		contextSpan=contextSpan_;
		supportK=supportK_;
		weak=new boolean[contextSpan-supportK+1];
		word=new Kmer(supportK, false);
		neighbor=new Kmer(supportK, false);
		assert(word.kbig==supportK && !word.hollow) : "The rolling filter requires the requested exact Solid words, not rounded or Hollow geometry.";
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Initializes every contained short word in the prefix ending at length (exclusive).
	 * One O(contextSpan) scan initializes the ring; subsequent appends are O(1).
	 * Reset is required for each new read, orientation or nonsequential edit. */
	public void reset(final byte[] bases, final int length){
		if(bases==null || length<contextSpan || length>bases.length){
			throw new IllegalArgumentException("Short-support reset needs a complete long context inside the saved prefix.");
		}
		ready=false;
		word.clear();
		lowWindows=0; cursor=0; lastChosen=-1;
		int windows=0;
		for(int i=length-contextSpan; i<length; i++){
			final int base=AminoAcid.baseToNumber[bases[i]];
			if(base<0){throw new IllegalArgumentException("Ambiguous short-support context at position "+i);}
			word.addRightNumeric(base);
			if(i>=length-contextSpan+supportK-1){
				weak[windows]=count(word)<=1;
				if(weak[windows]){lowWindows++;}
				windows++;
			}
		}
		assert(windows==weak.length) : "The ring must contain contextSpan-supportK+1 complete Solid windows; missing slots would lose weak evidence.";
		ready=true;
	}

	/** Advances exactly one accepted saved base, replacing the oldest contained word. */
	public void advance(final int base){
		requireReady(base);
		word.addRightNumeric(base);
		if(weak[cursor]){lowWindows--;}
		weak[cursor]=count(word)<=1;
		if(weak[cursor]){lowWindows++;}
		cursor++;
		if(cursor==weak.length){cursor=0;}
		lastChosen=-1;
		assert(lowWindows>=0 && lowWindows<=weak.length) : "Replacing one ring slot must preserve the number of weak windows within its capacity.";
	}

	/** Permits a pre-existing long-span cut only with weak context or weak next-key support.
	 * Repeated candidate queries do not advance the context or change its ring. */
	public boolean allowsCut(final int candidate){
		requireReady(candidate);
		neighbor.setFrom(word);
		neighbor.addRightNumeric(candidate);
		lastChosen=count(neighbor);
		return lowWindows>0 || lastChosen<=1;
	}

	/** Returns the number of depth0/1 words entirely within the current long context. */
	public int lowWindowCount(){return lowWindows;}
	/** Returns the most recently queried next-key count, or-1 after reset/advance. */
	public int chosenCount(){return lastChosen;}

	/** Rejects use before reset or an undefined appended/candidate base. */
	private void requireReady(final int base){
		if(!ready || base<0 || base>3){
			throw new IllegalArgumentException("Reset short support before querying/appending a defined ACGT base.");
		}
	}

	/** Keeps lookup failures distinct from a genuine zero count. */
	private int count(final Kmer key){
		assert(key.kbig==supportK && !key.hollow && key.len>=supportK) : "The lookup must receive a complete exact Solid word matching the configured table geometry.";
		final int value=lookup.count(key);
		if(value<0){throw new IllegalArgumentException("Short-support lookup must normalize absent keys to zero.");}
		return value;
	}

	/** Validates geometry once and normalizes only the table's absent-key sentinel. */
	private static CountLookup tableLookup(final KmerTableSetU tables){
		if(tables==null || tables.hollow || (tables.kbig!=64 && tables.kbig!=96) || tables.hashOnly()){
			throw new IllegalArgumentException("Short support needs an explicit Solid64 or Solid96 count table.");
		}
		return new CountLookup(){
			@Override
			public int count(final Kmer key){return Math.max(0, tables.getCount(key));}
		};
	}

	/*--------------------------------------------------------------*/
	/*----------------            Types             ----------------*/
	/*--------------------------------------------------------------*/

	/** Immutable count source; implementations must not retain or mutate scratch keys. */
	public interface CountLookup {
		int count(Kmer key);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final CountLookup lookup;
	private final int contextSpan, supportK;
	private final boolean[] weak;
	private final Kmer word, neighbor;
	private int lowWindows, cursor, lastChosen=-1;
	private boolean ready;
}
