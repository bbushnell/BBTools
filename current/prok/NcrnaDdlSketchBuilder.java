package prok;

import java.util.Arrays;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import kmer.AbstractKmerTable;
import map.LongHashSet;
import map.IntHashSet;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import shared.Tools;
import structures.IntList;

/** Native forward8/masked10 quarter-short sketches. Each worker reuses one
 * builder; returned DDL objects own only their native key arrays. The blacklist
 * is borrowed read-only and never learned from a clustering query.
 * @author Brian, Raiden, Ganyu
 */
public final class NcrnaDdlSketchBuilder implements CardinalityTracker.KmerConsumer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaDdlSketchBuilder(final int bins_, final long seed_, final LongHashSet blacklist_){
		this(bins_, seed_, blacklist_, false);
	}

	/** Fixed encoded keys do not need or change the caller's global DDL exponent. */
	public static NcrnaDdlSketchBuilder forKeyArrays(final int bins, final long seed, final LongHashSet blacklist){
		return new NcrnaDdlSketchBuilder(bins, seed, blacklist, true);
	}

	private NcrnaDdlSketchBuilder(final int bins_, final long seed_, final LongHashSet blacklist_, final boolean arraysOnly){
		if((bins_!=512 && bins_!=1024) || !arraysOnly && DynamicDemiLog.exponentBits()!=4 || blacklist_==null){
			throw new IllegalArgumentException("Quarter-short sketches need explicit512/1024bins, native4+12keys and an explicit padding-word set");
		}
		bins=bins_;
		seed=seed_;
		blacklist=blacklist_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Sketches           ----------------*/
	/*--------------------------------------------------------------*/

	public DynamicDemiLog newSketch(){
		if(DynamicDemiLog.exponentBits()!=4){throw new IllegalStateException("Native DDL objects require exponent4; use fixed key arrays without changing a caller's global exponent");}
		return DynamicDemiLog.create(bins, 10, seed, 0, false, false);
	}

	/** Replaces reusable per-worker maxima with explicit native4+12 keys. */
	public void fillKeys(final byte[] sequence, final char[] output){
		fillKeys(sequence, output, null);
	}

	/** Optional query-only admission before maximum reduction. Reference profiles
	 * and the standard native key geometry remain unchanged. Null is the ordinary
	 * query path; the immutable existing-index lookup is shared across workers.
	 */
	public void fillKeys(final byte[] sequence, final char[] output, final NcrnaDdlReferenceKeys allowed){
		fillKeys(sequence, output, allowed, null);
	}
	/** Optional query-work census; keep disabled during query timing. */
	public void fillKeys(final byte[] sequence, final char[] output, final NcrnaDdlReferenceKeys allowed, final QueryStats stats){
		if(sequence==null || output==null || output.length!=bins || active!=null || activePrefix!=null || activeKeys!=null){
			throw new IllegalArgumentException("A worker fills one complete matching key array at a time");
		}
		if(allowed!=null && (allowed.bins!=bins || allowed.seed!=seed)){
			throw new IllegalArgumentException("Query reference-key filter must share this builder's exact bin width and hash seed");
		}
		if(stats!=null){
			if(stats.bins!=bins){throw new IllegalArgumentException("Query-work counters must match this builder's bin width");}
			stats.clear();
		}
		Arrays.fill(output, (char)0);
		activeKeys=output;
		activeReferenceKeys=allowed;
		activeStats=stats;
		try{CardinalityTracker.forEachKmer(sequence, 8, 10, false, true, this);}
		finally{
			activeKeys=null;
			activeReferenceKeys=null;
			activeStats=null;
		}
	}

	/** Fills caller-owned native scratch; no per-word object or intermediate
	 * hash/key array. Ambiguity interrupts the entire original span before masking. */
	public void fill(final byte[] sequence, final DynamicDemiLog output){
		if(sequence==null || output==null || output.buckets!=bins || active!=null || activePrefix!=null || activeKeys!=null){
			throw new IllegalArgumentException("A worker fills one complete matching native DDL at a time");
		}
		output.clear();
		active=output;
		try{CardinalityTracker.forEachKmer(sequence, 8, 10, false, true, this);}
		finally{active=null;}
	}

	/** Retains every distinct eligible consensus bin/key before maximum reduction.
	 * Multiple words compressed to the same pair receive the virtual bonus once. */
	public void fillPrefix(final byte[] sequence, final Prefix output){
		if(sequence==null || output==null || output.bins!=bins || active!=null || activePrefix!=null || activeKeys!=null){
			throw new IllegalArgumentException("Consensus prefix enumeration requires matching worker-owned scratch");
		}
		output.clear();
		activePrefix=output;
		try{CardinalityTracker.forEachKmer(sequence, 8, 10, false, true, this);}
		finally{activePrefix=null;}
	}

	@Override
	public void accept(final long word, final int k){
		assert((active!=null || activePrefix!=null || activeKeys!=null) && (k==8 || k==10))
			: "The shared forward iterator must preserve the declared8/10word lengths before native routing";
		final boolean shortWord=k==8;
		final long tagged=(shortWord ? word : word&LONG_MASK)
			|BuildNcrnaDdlLayoutCounts.wordTag(k, shortWord, true);
		if(blacklist.contains(tagged)){return;}
		final long hash=Tools.hash64shift(tagged^seed);
		final int partition=BuildNcrnaDdlLayoutCounts.partitionFor(hash, bins, 2, shortWord);
		if(active!=null){active.storeHashed(tagged, hash, partition, 2);}
		else{
			//Exact high-selector replacement used by DynamicDemiLog.storeHashed;
			//keep every encoded key here instead of that method's bin maximum.
			final int bucket=((int)hash&((bins-1)>>>2))|(partition<<(Integer.numberOfTrailingZeros(bins)-2));
			final int key=DynamicDemiLog.encodeKey(hash, 4, 12);
			if(activePrefix!=null){activePrefix.add((bucket<<16)|key);}
			else{
				if(activeStats!=null){activeStats.observe(bucket, key, activeKeys[bucket]);}
				if(key>activeKeys[bucket]){
					final boolean present=activeReferenceKeys==null || activeReferenceKeys.contains(bucket, key);
					if(activeStats!=null && activeReferenceKeys!=null){
						activeStats.postingLookups++;
						if(!present){activeStats.rejectedCandidates++;}
					}
					if(present){activeKeys[bucket]=(char)key;}
				}
			}
		}
	}

	/** Events after the existing padding-word blacklist. Occupied visits include
	 * repeats of the incumbent key; different-key contests count challenger words,
	 * not unique key species. contestedBins counts each affected bin once.
	 */
	public static final class QueryStats{
		public QueryStats(final int bins_){
			if(bins_!=512 && bins_!=1024){throw new IllegalArgumentException("Query census requires the native bin width");}
			bins=bins_;
			contested=new boolean[bins];
		}
		private void clear(){
			words=0;
			occupiedVisits=0;
			differentKeyContests=0;
			contestedBins=0;
			postingLookups=0;
			rejectedCandidates=0;
			Arrays.fill(contested, false);
		}
		private void observe(final int bin, final int key, final int incumbent){
			words++;
			if(incumbent!=0){
				occupiedVisits++;
				if(key!=incumbent){
					differentKeyContests++;
					if(!contested[bin]){
						contested[bin]=true;
						contestedBins++;
					}
				}
			}
		}
		public long words, occupiedVisits, differentKeyContests, postingLookups, rejectedCandidates;
		public int contestedBins;
		public final int bins;
		private final boolean[] contested;
	}

	public static final class Prefix implements NcrnaDdlLayoutVotes.KeyPrefix{
		public Prefix(final int bins_){
			if(bins_!=512 && bins_!=1024){throw new IllegalArgumentException("Consensus prefix width must match the explicit native512/1024layout");}
			bins=bins_;
		}
		/** Historical bonus: exactly one native maximum per nonempty bin. */
		public void setMaxima(final char[] maxima){
			if(maxima==null || maxima.length!=bins){throw new IllegalArgumentException("Consensus maxima must match this prefix width");}
			clear();
			for(int bin=0; bin<bins; bin++){
				if(maxima[bin]!=0){add((bin<<16)|maxima[bin]);}
			}
		}
		private void clear(){
			keys.clear();
			packed.clear();
		}
		private void add(final int key){if(keys.add(key)){packed.add(key);}}
		@Override public int buckets(){return bins;}
		@Override public int keyCount(){return packed.size;}
		@Override public int packedKey(final int index){return packed.get(index);}
		@Override public boolean contains(final int key){return keys.contains(key);}
		private final int bins;
		private final IntHashSet keys=new IntHashSet(256);
		private final IntList packed=new IntList();
	}

	/*--------------------------------------------------------------*/
	/*----------------          Blacklist           ----------------*/
	/*--------------------------------------------------------------*/

	/** Reuses the existing exported padding-table format. Scope is selected
	 * explicitly, never inferred from the filename or treated as an all-gene fit.
	 * Only the decided256short/2048long prefixes are retained. */
	public static LongHashSet readBlacklist(final String path, final String scope){
		if(path==null || scope==null || scope.isEmpty()){
			throw new IllegalArgumentException("The padding word file and its literal fit_side scope are required");
		}
		final LongHashSet words=new LongHashSet(4096);
		final int[] counts=new int[2];
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, BLACKLIST_HEADER, path);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				if(p.terms()!=8){throw new IllegalArgumentException("Padding table requires eight fields: "+path);}
				if(!p.termEquals("TRUE_WINDOW_PADDING", 0) || !p.termEquals(scope, 7)){continue;}
				final boolean shortWord=p.termEquals("SHORT8", 1);
				if(!shortWord && !p.termEquals("MASKED10", 1)){
					throw new IllegalArgumentException("Selected padding scope must contain only SHORT8 and MASKED10 words");
				}
				final int part=shortWord ? 0 : 1, limit=shortWord ? 256 : 2048;
				final int rank=p.parseInt(2);
				if(rank>limit){continue;}
				final long raw=p.parseLong(3), tagged=p.parseLong(4);
				final int k=shortWord ? 8 : 10;
				if(rank!=counts[part]+1 || raw<0 || raw>=(1L<<(2*k))
						|| !shortWord && (raw&LONG_MASK)!=raw
						|| tagged!=(raw|BuildNcrnaDdlLayoutCounts.wordTag(k, shortWord, true))
						|| p.parseLong(5)<1 || !words.add(tagged)){
					throw new IllegalArgumentException("Padding prefixes need ordered ranks, valid native forward words/tags and distinct positive-frequency entries: "+path);
				}
				counts[part]++;
			}
		}
		if(counts[0]!=256 || counts[1]!=2048){
			throw new IllegalArgumentException("Selected padding scope lacks the decided256short/2048long prefix: "+counts[0]+"/"+counts[1]);
		}
		assert(words.size()==2304) : "Distinct length tags must conserve both complete padding prefixes";
		return words;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public final int bins;
	public final long seed;
	private final LongHashSet blacklist;
	private DynamicDemiLog active;
	private Prefix activePrefix;
	private char[] activeKeys;
	private NcrnaDdlReferenceKeys activeReferenceKeys;
	private QueryStats activeStats;
	private static final long LONG_MASK=AbstractKmerTable.makeMiddleMask(10, false, 1);
	public static final String BLACKLIST_HEADER="definition\tword_form\trank\traw_word\ttagged_word\tfrequency\tforward_sequence\tfit_side";
}
