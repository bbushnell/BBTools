package assemble;

import java.util.Arrays;
import ukmer.Kmer;

/** Reusable, non-mutating 3+1+4 local-edit count probes at a supplied position.
 * This enumerates evidence only: no localization, acceptance, or read editing.
 * Caller may stop after substitutions; indel probes are lazy and cached.
 * @author Fischl */
final class LocalEditKmerProbe {
	LocalEditKmerProbe(final int k_,final HomopolymerIndelProposal.CountLookup lookup_){this(k_,lookup_,1);}
	LocalEditKmerProbe(final int k_,final HomopolymerIndelProposal.CountLookup lookup_,final int windows_){
		if(k_<5 || lookup_==null){throw new IllegalArgumentException("Probe requires K>=5 and an immutable count lookup.");}
		if(windows_<1 || windows_>k_ || (windows_&1)==0){throw new IllegalArgumentException("Probe window count must be positive, odd and <=K.");}
		k=k_;lookup=lookup_;key=new Kmer(k);
		if(key.kbig!=k){throw new IllegalArgumentException("Probe K must equal effective table K.");}
		windows=windows_;single=windows==1 ? null : new LocalEditKmerProbe(k,lookup);
	}
	/** Optional reuse across calls over an immutable read. Call again after any
	 * read edit, even if the byte-array identity stays the same. Other callers
	 * retain per-call invalidation and may change their input between calls. */
	void beginRead(final byte[] read){
		assert(read!=null) : "A cached original window belongs to a specific immutable read.";
		immutableRead=read;cachedBases=null;bases=null;substitutionsComplete=false;indelsEvaluated=false;
		if(single!=null){single.beginRead(read);}
	}
	private void startCall(final byte[] read){
		assert(read!=null) : "Validated probe input must exist before establishing an immutable cache scope.";
		if(immutableRead!=read){
			cachedBases=null;
			// The child is private; the input is stable for this outer call.
			if(single!=null){single.beginRead(read);}
		}
	}
	/** Original k-base window [windowStart,windowStart+k), position inside it.
	 * Requires a called A/C/G/T at position (three, not four substitutions).
	 * Each other undefined base makes a particular candidate window unavailable.
	 * The input/table must not change until the next call to substitutions. */
	void substitutions(final byte[] bases_,final int windowStart_,final int position_){
		final int original=prepare(bases_,windowStart_,position_);
		startCall(bases_);
		if(windows>1){probeWindows(false);substitutionsComplete=true;return;}
		if(!loadWindow(false)){substitutionsComplete=true;return;}
		substitutionsAvailable=true;
		try{
			for(int b=0;b<4;b++){
				if(b==original){continue;}
				key.substituteBase(position-windowStart,b);
				substitutionDepth[b]=count();
			}
		}finally{key.substituteBase(position-windowStart,original);}
		assert(queries==3) : "A fully defined called base has exactly three substitution alternatives.";
		substitutionsComplete=true;
	}
	/** After a caller's multi-position substitution stage found no repair, probe
	 * only the five indels here; do not repeat the three substitution lookups. */
	void indelsAt(final byte[] bases_,final int windowStart_,final int position_){
		prepare(bases_,windowStart_,position_);startCall(bases_);substitutionsComplete=true;indels();
	}
	private int prepare(final byte[] bases_,final int windowStart_,final int position_){
		return prepare(bases_,windowStart_,position_,false);
	}
	/** An inserted base can create one additional complete window at the right. */
	private int prepare(final byte[] bases_,final int windowStart_,final int position_,final boolean insertedWindow){
		substitutionsComplete=false;bases=null;
		if(bases_==null || windowStart_<0 || (long)windowStart_+k>bases_.length+(insertedWindow ? 1L : 0L) ||
				position_<windowStart_ || position_>=windowStart_+k || position_>=bases_.length){
			throw new IllegalArgumentException("Substitution probes require an in-bounds kmer and a position inside it.");
		}
		final int original=code(bases_[position_]);
		if(original<0){throw new IllegalArgumentException("Three-alternative probe requires a called A/C/G/T base.");}
		bases=bases_;windowStart=windowStart_;position=position_;queries=0;indelsEvaluated=false;
		Arrays.fill(substitutionDepth,-1);Arrays.fill(insertionDepth,-1);deletionDepth=-1;
		substitutionsAvailable=deletionAvailable=insertionsAvailable=false;
		return original;
	}
	/** Deletion drops bases[position] and draws in bases[windowStart+k].
	 * Insertion inserts before position and trims the original last window base.
	 * Each resulting candidate is exactly k bases. No candidate sequence allocation. */
	void indels(){
		if(!substitutionsComplete){throw new IllegalStateException("Successfully probe substitutions first; failed probes require a fresh input call.");}
		if(indelsEvaluated){return;}
		try{probeIndels();indelsEvaluated=true;}
		catch(RuntimeException failure){substitutionsComplete=false;throw failure;}
		catch(Error failure){substitutionsComplete=false;throw failure;}
	}
	private void probeIndels(){
		if(windows>1){probeWindows(true);return;}
		probeDeletion();probeInsertions();
		assert(queries<=8) : "One locus requires at most three substitutions, one deletion and four inserted-base lookups.";
	}
	private void probeDeletion(){
		if((long)windowStart+k<bases.length){
			final int incoming=code(bases[windowStart+k]);
			if(incoming<0 || !loadWindow(false)){return;}
			final long removed=key.deleteBase(position-windowStart,incoming);
			try{deletionAvailable=true;deletionDepth=count();}
			finally{
				final long dropped=key.insertBase(position-windowStart,removed);
				assert(dropped==incoming) : "Restoring the deleted original base must evict the appended incoming base.";
			}
		}
	}
	private void probeInsertions(){
		if(!loadWindow(true)){return;}
		final long dropped=key.insertBase(position-windowStart,0);
		insertionsAvailable=true;
		try{
			insertionDepth[0]=count();
			for(int b=1;b<4;b++){
				key.substituteBase(position-windowStart,b);insertionDepth[b]=count();
			}
		}finally{key.deleteBase(position-windowStart,dropped);}
	}
	/** Restore/cache only an ORIGINAL window. Virtual edits always undo themselves.
	 * For insertion alone, an undefined/missing final base may be a dummy A: it
	 * is dropped before querying. No lookup may use that padded original window. */
	private boolean loadWindow(final boolean insertion){
		assert(bases!=null && windowStart>=0) : "prepare supplies a bounded original/virtual window before loading.";
		if(cachedBases==bases && cachedStart==windowStart){return insertion || cachedFull;}
		if(cachedBases==bases && cachedFull && Math.abs((long)windowStart-cachedStart)<k && (long)windowStart+k<=bases.length){
			while(cachedStart<windowStart){
				final int b=code(bases[cachedStart+k]);
				if(b<0){cachedBases=null;break;}
				key.addRightNumeric(b);cachedStart++;windowRolls++;
			}
			while(cachedBases!=null && cachedStart>windowStart){
				final int b=code(bases[cachedStart-1]);
				if(b<0){cachedBases=null;break;}
				key.addLeftNumeric(b);cachedStart--;windowRolls++;
			}
			if(cachedBases!=null){return true;}
		}
		cachedBases=null;key.clearFast();windowBuilds++;
		for(int i=0;i<k-1;i++){
			final int b=code(bases[windowStart+i]);
			if(b<0){return false;}key.addRightNumeric(b);
		}
		final int last=windowStart+k-1;
		final int b=last<bases.length ? code(bases[last]) : -1;
		key.addRightNumeric(Math.max(0,b));
		cachedBases=bases;cachedStart=windowStart;cachedFull=b>=0;
		return insertion || cachedFull;
	}
	/** Prefer starts centered on windowStart, shifting the group to keep every
	 * window in the edited read and spanning the edit. Use fewer distinct starts
	 * only if fewer exist. Undefined retained bases still veto a candidate; never
	 * drop a low count from the minimum. N=1 retains its legacy query contract. */
	private void probeWindows(final boolean indels){
		assert(single!=null && windows>1) : "Only multi-window probes delegate to the reusable single-window enumerator.";
		if(indels){
			probeWindowFamily(1);probeWindowFamily(2);
			deletionAvailable=deletionDepth>=0;insertionsAvailable=insertionDepth[0]>=0;
		}else{
			probeWindowFamily(0);
			for(int b=0;b<4;b++){substitutionsAvailable|=substitutionDepth[b]>=0;}
		}
		assert(queries<=8*windows) : "Each adjacent window contributes at most the same 3+1+4 candidate lookups.";
	}
	/** Family0=S,1=D,2=I. Start ranges refer to complete virtual edited windows.
	 * S/I must contain the changed/inserted base; D must straddle its junction. */
	private void probeWindowFamily(final int family){
		assert(family>=0 && family<=2 && single!=null) : "Only the three single-base edit families have these virtual coordinate ranges.";
		final int low=Math.max(0,position-k+1);
		final int high=Math.min(position-(family==1 ? 1 : 0),bases.length-k+(family==2 ? 1 : family==1 ? -1 : 0));
		if(high<low){return;}
		final int number=Math.min(windows,high-low+1);
		final int start=Math.max(low,Math.min(windowStart-windows/2,high-number+1));
		assert(start>=low && start+number-1<=high) : "Every selected start must be a complete window spanning this family's edit.";
		for(int i=0;i<number;i++){
			if(family==0){single.substitutions(bases,start+i,position);}
			else{
				single.prepare(bases,start+i,position,family==2);
				if(family==1){single.probeDeletion();}else{single.probeInsertions();}
			}
			queries+=single.queries;
			if(family==0){for(int b=0;b<4;b++){substitutionDepth[b]=i==0 ? single.substitutionDepth[b] : Math.min(substitutionDepth[b],single.substitutionDepth[b]);}}
			else if(family==1){deletionDepth=i==0 ? single.deletionDepth : Math.min(deletionDepth,single.deletionDepth);}
			else{for(int b=0;b<4;b++){insertionDepth[b]=i==0 ? single.insertionDepth[b] : Math.min(insertionDepth[b],single.insertionDepth[b]);}}
		}
	}
	private int count(){
		assert(key.len()>=k) : "Only complete edited kmers may reach the count table.";
		final int depth=lookup.count(key);queries++;
		if(depth<-1){throw new IllegalStateException("CountLookup permits depth>=0 or absent=-1, got "+depth);}
		return Math.max(0,depth);
	}
	private static int code(final byte b){return b=='A' ? 0 : b=='C' ? 1 : b=='G' ? 2 : b=='T' ? 3 : -1;}
	/** -1 means unavailable/not probed; a queried absent kmer has depth zero. */
	final int[] substitutionDepth=new int[4],insertionDepth=new int[4];
	int deletionDepth,queries;
	boolean substitutionsAvailable,deletionAvailable,insertionsAvailable;
	private boolean indelsEvaluated,substitutionsComplete;
	private byte[] bases;
	private byte[] immutableRead,cachedBases;
	private int cachedStart;
	private boolean cachedFull;
	private long windowBuilds,windowRolls;
	long windowBuilds(){return single==null ? windowBuilds : single.windowBuilds;}
	long windowRolls(){return single==null ? windowRolls : single.windowRolls;}
	private int windowStart,position;
	private final int k;
	private final int windows;
	private final LocalEditKmerProbe single;
	private final Kmer key;
	private final HomopolymerIndelProposal.CountLookup lookup;
}
