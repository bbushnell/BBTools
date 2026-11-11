package assemble;

import java.util.Arrays;
import ukmer.Kmer;

/**
 * Enumerates three substitution, one deletion and four insertion counts at a supplied locus.
 * This is evidence only: no localization, acceptance policy or input editing. Depth arrays
 * use A/C/G/T indices 0..3. A queried absent key has depth 0; unavailable or unprobed
 * alternatives have depth -1 after a probe call, including the unchanged substitution allele.
 * Multiple-window mode takes the minimum for each alternative, retaining any unavailable
 * window as a veto. It does not average counts or omit an unsupported window.
 *
 * One mutable instance belongs to one worker. Read arrays and the count lookup are borrowed;
 * the lookup must not retain or alter the reusable Kmer sequence, and its layout/core-mask
 * settings must remain fixed as required by CountLookup. Input/counts must stay stable while
 * a locus's lazy results are in use. Only consume results after a successful probe call;
 * exceptions can leave partial arrays/counters. beginRead enables additional reuse across
 * loci on an immutable array and must be repeated after any edit, even in the same array.
 * @author Fischl */
public final class LocalEditKmerProbe{

	/** Creates a single-window probe using the supplied effective table K. */
	public LocalEditKmerProbe(final int k_, final HomopolymerIndelProposal.CountLookup lookup_){this(k_, lookup_, 1);}

	/** Allocates reusable keys/results and, for multiple windows, one single-window child.
	 * @param k_ Effective table K, at least 5; silently rounded K values are rejected
	 * @param lookup_ Nonnull immutable count lookup, following CountLookup's key-ownership contract
	 * @param windows_ Positive odd maximum number of adjacent windows, no greater than K
	 */
	public LocalEditKmerProbe(final int k_, final HomopolymerIndelProposal.CountLookup lookup_, final int windows_){
		if(k_<5 || lookup_==null){throw new IllegalArgumentException("Probe requires K>=5 and an immutable count lookup.");}
		if(windows_<1 || windows_>k_ || (windows_&1)==0){throw new IllegalArgumentException("Probe window count must be positive, odd and <=K.");}
		k=k_; lookup=lookup_; key=new Kmer(k);
		if(key.kbig!=k){throw new IllegalArgumentException("Probe K must equal effective table K.");}
		windows=windows_; single=windows==1 ? null : new LocalEditKmerProbe(k, lookup);
	}

	/** Starts an optional immutable-array cache scope; no copy is made.
	 * Call again after any read edit, even when array identity is unchanged. Without this
	 * scope, each new input call invalidates its original-window cache. This invalidates
	 * lazy-call state but does not clear published depth arrays or lifetime counters.
	 * @param read Nonnull array that will remain unchanged until the scope is replaced
	 */
	public void beginRead(final byte[] read){
		assert(read!=null) : "A cached original window belongs to a specific immutable read.";
		immutableRead=read; cachedBases=null; bases=null; substitutionsComplete=false; indelsEvaluated=false;
		if(single!=null){single.beginRead(read);}
	}

	/** Invalidates cross-call reuse unless this array is the explicitly scoped immutable read. */
	private void startCall(final byte[] read){
		assert(read!=null) : "Validated probe input must exist before establishing an immutable cache scope.";
		if(immutableRead!=read){
			cachedBases=null;
			// The child is private; the input is stable for this outer call.
			if(single!=null){single.beginRead(read);}
		}
	}

	/** Resets locus results and probes substitutions, leaving indels lazy.
	 * Requires an uppercase A/C/G/T at position; all other bytes are undefined.
	 * Undefined retained bases make a candidate window unavailable. The input/table must
	 * remain unchanged until the next substitutions or indelsAt call; an explicit beginRead
	 * scope imposes the longer immutable-array lifetime described there.
	 * @param bases_ Borrowed read bases
	 * @param windowStart_ Start of an in-bounds original k-base window; also the preferred center start in multi-window mode
	 * @param position_ Original-read base index inside that window
	 */
	public void substitutions(final byte[] bases_, final int windowStart_, final int position_){
		final int original=prepare(bases_, windowStart_, position_);
		startCall(bases_);
		if(windows>1){probeWindows(false); substitutionsComplete=true; return;}
		if(!loadWindow(false)){substitutionsComplete=true; return;}
		substitutionsAvailable=true;
		try{
			for(int b=0; b<4; b++){
				if(b==original){continue;}
				key.substituteBase(position-windowStart, b);
				substitutionDepth[b]=count();
			}
		}finally{key.substituteBase(position-windowStart, original);}
		assert(queries==3) : "A fully defined called base has exactly three substitution alternatives.";
		substitutionsComplete=true;
	}

	/** Starts a fresh locus and probes only deletion/insertion alternatives.
	 * No preceding substitution call is required; substitution results are reset to unavailable.
	 * Coordinates and input ownership are the same as substitutions. This class does not
	 * decide whether a caller should also compare substitution or competing-indel evidence.
	 */
	public void indelsAt(final byte[] bases_, final int windowStart_, final int position_){
		prepare(bases_, windowStart_, position_); startCall(bases_); substitutionsComplete=true; indels();
	}

	/** Validates an ordinary window and resets per-locus results before probing. */
	private int prepare(final byte[] bases_, final int windowStart_, final int position_){
		return prepare(bases_, windowStart_, position_, false);
	}

	/** Resets per-locus state; an insertion may permit one extra complete right-edge window.
	 * The requested original position must still contain uppercase A/C/G/T. Lifetime
	 * window-build/roll counters are not reset. An invalid request invalidates lazy-call state.
	 * @return Numeric original allele 0..3
	 */
	private int prepare(final byte[] bases_, final int windowStart_, final int position_, final boolean insertedWindow){
		substitutionsComplete=false; bases=null;
		if(bases_==null || windowStart_<0 || (long)windowStart_+k>bases_.length+(insertedWindow ? 1L : 0L) ||
				position_<windowStart_ || position_>=windowStart_+k || position_>=bases_.length){
			throw new IllegalArgumentException("Substitution probes require an in-bounds kmer and a position inside it.");
		}
		final int original=code(bases_[position_]);
		if(original<0){throw new IllegalArgumentException("Three-alternative probe requires a called A/C/G/T base.");}
		bases=bases_; windowStart=windowStart_; position=position_; queries=0; indelsEvaluated=false;
		Arrays.fill(substitutionDepth, -1); Arrays.fill(insertionDepth, -1); deletionDepth=-1;
		substitutionsAvailable=deletionAvailable=insertionsAvailable=false;
		return original;
	}

	/** Completes lazy indel results for the current successful locus; repeated calls do no work.
	 * Deletion drops bases[position] and draws in bases[windowStart+k].
	 * Insertion inserts before position and trims the original last window base.
	 * Each candidate is exactly k bases; no candidate sequence array is allocated here.
	 * Requires successful substitutions or indelsAt preparation. After an indel failure,
	 * a new locus call is required; repeating indels cannot retry partial work.
	 */
	public void indels(){
		if(!substitutionsComplete){throw new IllegalStateException("Successfully probe substitutions first; failed probes require a fresh input call.");}
		if(indelsEvaluated){return;}
		try{probeIndels(); indelsEvaluated=true;}
		catch(RuntimeException failure){substitutionsComplete=false; throw failure;}
		catch(Error failure){substitutionsComplete=false; throw failure;}
	}

	/** Dispatches one or several windows without repeating substitution lookups. */
	private void probeIndels(){
		if(windows>1){probeWindows(true); return;}
		probeDeletion(); probeInsertions();
		assert(queries<=8) : "One locus requires at most three substitutions, one deletion and four inserted-base lookups.";
	}

	/** Queries one deletion if the original window and its extra incoming base are defined. */
	private void probeDeletion(){
		if((long)windowStart+k<bases.length){
			final int incoming=code(bases[windowStart+k]);
			if(incoming<0 || !loadWindow(false)){return;}
			final long removed=key.deleteBase(position-windowStart, incoming);
			try{deletionAvailable=true; deletionDepth=count();}
			finally{
				final long dropped=key.insertBase(position-windowStart, removed);
				assert(dropped==incoming) : "Restoring the deleted original base must evict the appended incoming base.";
			}
		}
	}

	/** Queries all four inserted alleles before position and restores the original key in finally. */
	private void probeInsertions(){
		if(!loadWindow(true)){return;}
		final long dropped=key.insertBase(position-windowStart, 0);
		insertionsAvailable=true;
		try{
			insertionDepth[0]=count();
			for(int b=1; b<4; b++){
				key.substituteBase(position-windowStart, b); insertionDepth[b]=count();
			}
		}finally{key.deleteBase(position-windowStart, dropped);}
	}

	/** Restores/caches only an ORIGINAL window. Virtual edits always undo themselves.
	 * For insertion alone, an undefined/missing final base may be a dummy A: it
	 * is dropped before querying. No lookup may use that padded original window. */
	private boolean loadWindow(final boolean insertion){
		assert(bases!=null && windowStart>=0) : "prepare supplies a bounded original/virtual window before loading.";
		if(cachedBases==bases && cachedStart==windowStart){return insertion || cachedFull;}
		if(cachedBases==bases && cachedFull && Math.abs((long)windowStart-cachedStart)<k && (long)windowStart+k<=bases.length){
			while(cachedStart<windowStart){
				final int b=code(bases[cachedStart+k]);
				if(b<0){cachedBases=null; break;}
				key.addRightNumeric(b); cachedStart++; windowRolls++;
			}
			while(cachedBases!=null && cachedStart>windowStart){
				final int b=code(bases[cachedStart-1]);
				if(b<0){cachedBases=null; break;}
				key.addLeftNumeric(b); cachedStart--; windowRolls++;
			}
			if(cachedBases!=null){return true;}
		}
		cachedBases=null; key.clearFast(); windowBuilds++;
		for(int i=0; i<k-1; i++){
			final int b=code(bases[windowStart+i]);
			if(b<0){return false;}key.addRightNumeric(b);
		}
		final int last=windowStart+k-1;
		final int b=last<bases.length ? code(bases[last]) : -1;
		key.addRightNumeric(Math.max(0, b));
		cachedBases=bases; cachedStart=windowStart; cachedFull=b>=0;
		return insertion || cachedFull;
	}

	/** Prefers starts centered on windowStart, shifting the group to keep every
	 * window in the edited read and spanning the edit. Use fewer distinct starts
	 * only if fewer exist. Undefined retained bases still veto a candidate; never
	 * drop a low count from the minimum. N=1 retains its legacy query contract. */
	private void probeWindows(final boolean indels){
		assert(single!=null && windows>1) : "Only multi-window probes delegate to the reusable single-window enumerator.";
		if(indels){
			probeWindowFamily(1); probeWindowFamily(2);
			deletionAvailable=deletionDepth>=0; insertionsAvailable=insertionDepth[0]>=0;
		}else{
			probeWindowFamily(0);
			for(int b=0; b<4; b++){substitutionsAvailable|=substitutionDepth[b]>=0;}
		}
		assert(queries<=8*windows) : "Each adjacent window contributes at most the same 3+1+4 candidate lookups.";
	}

	/** Family 0=S, 1=D, 2=I. Start ranges refer to complete virtual edited windows.
	 * S/I must contain the changed/inserted base; D must straddle its junction. */
	private void probeWindowFamily(final int family){
		assert(family>=0 && family<=2 && single!=null) : "Only the three single-base edit families have these virtual coordinate ranges.";
		final int low=Math.max(0, position-k+1);
		final int high=Math.min(position-(family==1 ? 1 : 0), bases.length-k+(family==2 ? 1 : family==1 ? -1 : 0));
		if(high<low){return;}
		final int number=Math.min(windows, high-low+1);
		final int start=Math.max(low, Math.min(windowStart-windows/2, high-number+1));
		assert(start>=low && start+number-1<=high) : "Every selected start must be a complete window spanning this family's edit.";
		for(int i=0; i<number; i++){
			if(family==0){single.substitutions(bases, start+i, position);}else{
				single.prepare(bases, start+i, position, family==2);
				if(family==1){single.probeDeletion();}else{single.probeInsertions();}
			}
			queries+=single.queries;
			if(family==0){
				for(int b=0; b<4; b++){
					substitutionDepth[b]=i==0 ? single.substitutionDepth[b] : Math.min(substitutionDepth[b], single.substitutionDepth[b]);
				}
			}else if(family==1){
				deletionDepth=i==0 ? single.deletionDepth : Math.min(deletionDepth, single.deletionDepth);
			}else{
				for(int b=0; b<4; b++){
					insertionDepth[b]=i==0 ? single.insertionDepth[b] : Math.min(insertionDepth[b], single.insertionDepth[b]);
				}
			}
		}
	}

	/** Queries the reusable edited key, maps absent -1 to 0 and rejects depths below -1.
	 * queries counts returned lookup calls, including an invalid negative result, but
	 * excludes a lookup that throws before returning. The lookup must not retain the key.
	 */
	private int count(){
		assert(key.len()>=k) : "Only complete edited kmers may reach the count table.";
		final int depth=lookup.count(key); queries++;
		if(depth<-1){throw new IllegalStateException("CountLookup permits depth>=0 or absent=-1, got "+depth);}
		return Math.max(0, depth);
	}

	/** Returns 0..3 only for uppercase A/C/G/T, otherwise -1. */
	private static int code(final byte b){return b=='A' ? 0 : b=='C' ? 1 : b=='G' ? 2 : b=='T' ? 3 : -1;}

	/*--------------------------------------------------------------*/
	/*----------------       Results and Cache       ----------------*/
	/*--------------------------------------------------------------*/

	/** Owned result arrays exposed for reading, indexed A/C/G/T; reset to -1 by each prepared locus. */
	public final int[] substitutionDepth=new int[4], insertionDepth=new int[4];
	/** Latest deletion depth (-1 unavailable) and successful-locus lookup count; partial failures may undercount. */
	public int deletionDepth, queries;
	/** Candidate availability, independent of whether its observed depth is positive. */
	public boolean substitutionsAvailable, deletionAvailable, insertionsAvailable;
	/** Lazy indel completion and permission to consume/probe the current locus. */
	private boolean indelsEvaluated, substitutionsComplete;
	/** Borrowed bases for the current prepared locus. */
	private byte[] bases;
	/** Explicit immutable-array scope and identity of the window currently held in key. */
	private byte[] immutableRead, cachedBases;
	/** Original-read start represented by the cached key. */
	private int cachedStart;
	/** Whether the cached original window is fully called, rather than insertion-only padding. */
	private boolean cachedFull;
	/** Lifetime full-load attempts and successful one-base rolling steps; beginRead does not reset them. */
	private long windowBuilds, windowRolls;
	/** Returns cumulative full-load attempts, using the private child in multi-window mode. */
	public long windowBuilds(){return single==null ? windowBuilds : single.windowBuilds;}
	/** Returns cumulative one-base rolls, using the private child in multi-window mode. */
	public long windowRolls(){return single==null ? windowRolls : single.windowRolls;}
	/** Preferred original window start and original-read edited position. */
	private int windowStart, position;
	/** Effective, fixed table K. */
	private final int k;
	/** Requested odd maximum number of windows per edit family. */
	private final int windows;
	/** Private reusable single-window worker, null for a single-window probe. */
	private final LocalEditKmerProbe single;
	/** Owned mutable query word, restored after virtual edits; lookup borrows it temporarily. */
	private final Kmer key;
	/** Borrowed immutable count provider using this key's fixed table layout/canonicalization. */
	private final HomopolymerIndelProposal.CountLookup lookup;
}
