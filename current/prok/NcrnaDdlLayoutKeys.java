package prok;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import kmer.AbstractKmerTable;
import map.IntHashSet;
import shared.Tools;
import structures.IntList;
import structures.LongList;

/** Reusable standard DDL objects for one mixed sketch or two length parts.
 * Every part is made by the existing factory and filled by native insertion.
 * Borrowed key arrays serve member voting only; scoring uses the DDL API.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlLayoutKeys implements NcrnaDdlLayoutVotes.KeyPrefix{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaDdlLayoutKeys(final int bins_, final boolean partitioned_, final boolean includeShort_,
			final NcrnaDdlKeyFormat format_){
		this(bins_, partitioned_ ? bins_/3 : 0, includeShort_, format_, 7, 10, true);
	}

	public NcrnaDdlLayoutKeys(final int bins_, final int shortBins_, final boolean includeShort_,
			final NcrnaDdlKeyFormat format_, final int shortK_, final int longK_, final boolean maskedLong_){
		if(bins_<1 || bins_>32768 || shortBins_<0 || shortBins_>=bins_ || format_==null
				|| shortK_<1 || shortK_>=longK_ || longK_>15){
			throw new IllegalArgumentException("Native DDL layouts require explicit word lengths and valid part widths");
		}
		if(shortBins_==0 ? !powerOfTwo(bins_) : !powerOfTwo(shortBins_) || !powerOfTwo(bins_-shortBins_)){
			throw new IllegalArgumentException("Every actual DDL part must have a power-of-two width, not a custom range");
		}
		if(DynamicDemiLog.exponentBits()!=format_.exponentBits || format_.mantissaBits!=16-format_.exponentBits){
			throw new IllegalArgumentException("Configure the native DDL exponent before creating any layout objects");
		}
		bins=bins_;
		shortBins=shortBins_;
		includeShort=includeShort_;
		format=format_;
		shortK=shortK_;
		longK=longK_;
		maskedLong=maskedLong_;
		if(shortBins==0){
			parts=new DynamicDemiLog[]{DynamicDemiLog.create(bins, longK, 12345, 0, false, false)};
			maxima=parts[0].maxArray();
		}else{
			parts=new DynamicDemiLog[]{
				DynamicDemiLog.create(shortBins, shortK, 12345, 0, false, false),
				DynamicDemiLog.create(bins-shortBins, longK, 12345, 0, false, false)};
			maxima=new char[bins];
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Construction         ----------------*/
	/*--------------------------------------------------------------*/

	public void build(final Raw raw){
		if(raw==null || raw.shortK!=shortK || raw.longK!=longK || raw.maskedLong!=maskedLong){
			throw new IllegalArgumentException("The shared word reader must match the layout's lengths and mask");
		}
		buildWords(raw.shortWords.array, raw.shortWords.size, raw.longWords.array, raw.longWords.size);
	}

	/** The prefixes contain un-hashed forward words, never preassigned bins. */
	public void buildWords(final int[] shortWords, final int shortCount, final int[] longWords, final int longCount){
		if(shortWords==null || longWords==null || shortCount<0 || shortCount>shortWords.length
				|| longCount<0 || longCount>longWords.length){
			throw new IllegalArgumentException("Native DDL construction consumes only the live raw-word prefixes");
		}
		for(final DynamicDemiLog part:parts){part.clear();}
		maximaReady=prefixReady=false;
		if(shortBins==0){
			final DynamicDemiLog ddl=parts[0];
			if(includeShort){
				for(int i=0; i<shortCount; i++){ddl.hashAndStore(shortWords[i]^SHORT_TAG);}
			}
			for(int i=0; i<longCount; i++){ddl.hashAndStore(longWords[i]^LONG_TAG);}
		}else{
			if(includeShort){
				for(int i=0; i<shortCount; i++){parts[0].hashAndStore(shortWords[i]);}
			}
			for(int i=0; i<longCount; i++){parts[1].hashAndStore(longWords[i]);}
		}
	}

	/** A concatenated borrowed view is used only while choosing reference votes.
	 * It is never wrapped as one non-power-of-two DDL or used by a score loop. */
	public char[] maximumKeys(){
		if(!maximaReady){
			if(parts.length==2){
				System.arraycopy(parts[0].maxArray(), 0, maxima, 0, shortBins);
				System.arraycopy(parts[1].maxArray(), 0, maxima, shortBins, bins-shortBins);
			}
			maximaReady=true;
		}
		return maxima;
	}

	public DynamicDemiLog part(final int index){return parts[index];}
	public int parts(){return parts.length;}

	/** Native consensus winners are virtual votes, one nonempty winner per bin. */
	private void preparePrefix(){
		if(prefixReady){return;}
		maximumKeys();
		prefix.clear();
		for(int bin=0; bin<bins; bin++){
			if(maxima[bin]!=0){prefix.add((bin<<16)|maxima[bin]);}
		}
		prefixReady=true;
	}

	public int buckets(){return bins;}
	public int keyCount(){preparePrefix(); return prefix.size;}
	public int packedKey(final int index){
		preparePrefix();
		assert(index>=0 && index<prefix.size) : "Consensus voting reads only the live native-winner prefix";
		return prefix.array[index];
	}
	public boolean contains(final int packed){
		final int bin=packed>>>16;
		final int key=packed&65535;
		return bin<bins && key!=0 && maximumKeys()[bin]==key;
	}
	public int bins(){return bins;}
	public int shortBins(){return shortBins;}
	public int longBins(){return shortBins==0 ? 0 : bins-shortBins;}
	public int shortK(){return shortK;}
	public int longK(){return longK;}
	public boolean maskedLong(){return maskedLong;}
	public NcrnaDdlKeyFormat format(){return format;}
	public boolean sameGeometry(final NcrnaDdlLayoutKeys other){
		return other!=null && bins==other.bins && shortBins==other.shortBins && includeShort==other.includeShort
			&& shortK==other.shortK && longK==other.longK && maskedLong==other.maskedLong
			&& format.exponentBits==other.format.exponentBits && format.mantissaBits==other.format.mantissaBits;
	}

	private static boolean powerOfTwo(final int value){return value>0 && (value&(value-1))==0;}

	/*--------------------------------------------------------------*/
	/*----------------          Raw Words           ----------------*/
	/*--------------------------------------------------------------*/

	/** Reused raw words/hashes for one declared short/long/mask combination.
	 * Ambiguous bases interrupt the full span even when the central base is masked.
	 * The native forward codes retain G3/T2. Distinct lengths are tagged before
	 * hashing; compression may subsequently make them indistinguishable.
	 */
	public static final class Raw implements CardinalityTracker.KmerConsumer{

		public Raw(){this(7, 10, true);}

		public Raw(final int shortK_, final int longK_, final boolean maskedLong_){
			if(shortK_<1 || shortK_>=longK_ || longK_>15){
				throw new IllegalArgumentException("Shared forward words require distinct short/long spans up to15bases");
			}
			shortK=shortK_;
			longK=longK_;
			maskedLong=maskedLong_;
			middleMask=AbstractKmerTable.makeMiddleMask(longK, false, maskedLong ? 1 : 0);
		}

		public void build(final byte[] bases){build(bases, false);}

		/** The extra untagged hash is confined to requested historical diagnostics. */
		public void build(final byte[] bases, final boolean legacyExponents){
			if(bases==null){throw new IllegalArgumentException("Shared raw keys require the actual sequence bytes");}
			seen.clear();
			shortWords.clear();
			longWords.clear();
			shortHashes.clear();
			longHashes.clear();
			shortExponent=longExponent=-1;
			legacyShortExponent=legacyLongExponent=-1;
			collectLegacy=legacyExponents;
			CardinalityTracker.forEachKmer(bases, shortK, longK, false, true, this);
			assert(shortWords.size<=Math.max(0, bases.length-shortK+1)
				&& longWords.size<=Math.max(0, bases.length-longK+1))
				: "Each complete forward span contributes at most one distinct word per declared length";
		}

		/** Internal primitive sink; one callback per complete unmasked input span. */
		@Override
		public void accept(final long word, final int k){
			assert(k==shortK || k==longK) : "Raw layout hashing must retain its declared length namespace";
			add((int)(k==shortK ? word : word&middleMask), k==shortK, collectLegacy);
		}

		private void add(final int word, final boolean shortWord, final boolean legacyExponents){
			assert(word>=0 && (word&1<<30)==0)
				: "At most fifteen bases leave bit30 free for the transient raw-word length namespace";
			if(!seen.add(word|(shortWord ? 1<<30 : 0))){return;}
			final long hash=Tools.hash64shift(word^(shortWord ? SHORT_TAG : LONG_TAG)^12345L);
			final int exponent=Long.numberOfLeadingZeros(hash);
			final int legacy=legacyExponents ? Long.numberOfLeadingZeros(Tools.hash64shift(word^12345L)) : -1;
			if(shortWord){
				shortWords.add(word);
				shortHashes.add(hash);
				shortExponent=Math.max(shortExponent, exponent);
				legacyShortExponent=Math.max(legacyShortExponent, legacy);
			}else{
				longWords.add(word);
				longHashes.add(hash);
				longExponent=Math.max(longExponent, exponent);
				legacyLongExponent=Math.max(legacyLongExponent, legacy);
			}
		}

		public IntList words(final boolean shortWord){return shortWord ? shortWords : longWords;}
		public LongList hashes(final boolean shortWord){return shortWord ? shortHashes : longHashes;}
		public int maxExponent(final boolean shortWord){return shortWord ? shortExponent : longExponent;}
		public int legacyMaxExponent(final boolean shortWord){return shortWord ? legacyShortExponent : legacyLongExponent;}
		public final int shortK, longK;
		public final boolean maskedLong;
		private final long middleMask;
		private boolean collectLegacy;
		private final IntHashSet seen=new IntHashSet(256);
		private final IntList shortWords=new IntList(), longWords=new IntList();
		private final LongList shortHashes=new LongList(), longHashes=new LongList();
		private int shortExponent=-1, longExponent=-1, legacyShortExponent=-1, legacyLongExponent=-1;
	}


	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int bins, shortBins, shortK, longK;
	private final boolean includeShort, maskedLong;
	private final NcrnaDdlKeyFormat format;
	private final DynamicDemiLog[] parts;
	private final char[] maxima;
	private final IntList prefix=new IntList();
	private boolean maximaReady, prefixReady;
	private static final long SHORT_TAG=1L<<62, LONG_TAG=1L<<63;
}
