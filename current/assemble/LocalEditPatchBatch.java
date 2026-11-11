package assemble;

import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import shared.Shared;

/** Worker-reused transaction for already validated original-coordinate patches.
 * Does NOT discover corrections or validate biological support. Replacements
 * include their final quality bytes; all local S/I/D operations count against
 * the read-wide budget, not merely the number of replacement intervals.
 * Borrows the target Read and its original arrays, but copies staged replacement
 * bytes so callers may reuse their scratch. Original arrays must remain immutable
 * until commit or clear; identity checks detect replacement, not in-place writes.
 * One worker owns the transaction and target throughout. Installation is not an
 * atomic publication mechanism for concurrent readers.
 * @author Fischl */
public final class LocalEditPatchBatch{

	/** Abandons an earlier transaction and clears result counters before validating
	 * an unpaired, alignment-free Read with bases and matching or null qualities.
	 * Borrows the input without copying it; failed validation leaves no active batch. */
	public void reset(final Read input){
		clear(); clearCounters(); HomopolymerIndelEdit.requireEditable(input);
		if(input.bases==null || (input.quality!=null && input.quality.length!=input.bases.length)){
			throw new IllegalArgumentException("Patch transaction requires intact bases and matching or null qualities.");
		}
		read=input; bases=input.bases; quality=input.quality;
	}
	/** Copies a replacement into worker-owned staging storage. Original intervals
	 * [from,to) must be disjoint and starts strictly increase; from==to inserts.
	 * Same-start operations must be combined into one patch by the caller.
	 * Validates count ranges and net length change, not whether the supplied operation
	 * counts describe the actual edit history. Substitutions and cancelling indels still
	 * consume budget. Replacement qualities preserve original nullness and must match
	 * the replacement-base length.
	 * Any append failure abandons the entire transaction, never a partial list.
	 * @param from Inclusive start in the original bases
	 * @param to Exclusive end in the original bases
	 * @param replacement Final bases for this interval; copied before return
	 * @param replacementQuality Final qualities, or null exactly when the original has none
	 * @param subs Number of local substitution operations
	 * @param ins Number of locally inserted bases
	 * @param dels Number of locally deleted bases
	 */
	public void add(final int from, final int to, final byte[] replacement, final byte[] replacementQuality,
		final int subs, final int ins, final int dels){
		boolean success=false;
		try{
			if(read==null){throw new IllegalStateException("Reset before staging a patch.");}
			if(from<0 || to<from || to>bases.length || replacement==null){throw new IllegalArgumentException("Patch interval/replacement outside the original read: "+from+".."+to);}
			if((quality==null)!=(replacementQuality==null) || (replacementQuality!=null && replacementQuality.length!=replacement.length)){
				throw new IllegalArgumentException("A patch must preserve quality nullness and supply one final quality per replacement base.");
			}
			final long operations=(long)subs+ins+dels, total=(long)stagedSubs+stagedIns+stagedDels+operations;
			if(subs<0 || ins<0 || dels<0 || operations<1 || total>Integer.MAX_VALUE){throw new IllegalArgumentException("Patch operation counts must be positive in total and fit the read-level counter.");}
			if((long)replacement.length-(to-from)!=(long)ins-dels){throw new IllegalArgumentException("Replacement length must agree with the actual insertion/deletion counts.");}
			if(spans.size>0 && (from<=spans.array[spans.size-4] || from<spans.array[spans.size-3])){
				throw new IllegalArgumentException("Patches overlap or share an original start; combine dependent edits before staging.");
			}
			if((long)replacementBases.length()+replacement.length>Shared.MAX_ARRAY_LEN){throw new IllegalArgumentException("Replacement staging exceeds the supported primitive array size.");}
			final int offset=replacementBases.length();
			replacementBases.append(replacement);
			if(quality!=null){replacementQualities.append(replacementQuality);}
			spans.add(from); spans.add(to); spans.add(offset); spans.add(replacement.length);
			stagedSubs+=subs; stagedIns+=ins; stagedDels+=dels;
			assert(spans.size%4==0 && (quality==null || replacementBases.length()==replacementQualities.length())) :
				"Each staged interval has four primitive fields and matching base/quality pool offsets.";
			success=true;
		}finally{if(!success){clear();}}
	}
	/** Returns committed operation count, -1 if over budget, or 0 for no patches.
	 * Requires a nonnegative budget. An active attempt clears result counters and
	 * closes the transaction on success, budget refusal or exception. An empty batch
	 * preserves array identities; otherwise installs complete new bases/qualities only
	 * after copying, without altering original buffers. Replaced input arrays are
	 * rejected before installation. Retained original qualities and supplied replacement
	 * qualities are copied verbatim; this class does not assign confidence scores.
	 * @param maxEdits Read-wide operation budget, not a limit on patch intervals
	 * @return Committed operation count, or -1 for whole-list budget refusal
	 */
	public int commit(final int maxEdits){
		if(read==null){throw new IllegalStateException("No active patch transaction.");}
		clearCounters();
		try{
			if(maxEdits<0){throw new IllegalArgumentException("Read-level edit limit must be nonnegative.");}
			HomopolymerIndelEdit.requireEditable(read);
			if(read.bases!=bases || read.quality!=quality){throw new IllegalStateException("Read arrays changed during patch staging.");}
			final int operations=stagedSubs+stagedIns+stagedDels;
			if(operations>maxEdits){return -1;}
			if(spans.size==0){return 0;}
			final int length=LocalEditBatch.resultLength(bases.length, (long)stagedIns-stagedDels);
			final byte[] output=new byte[length], outputQuality=quality==null ? null : new byte[length];
			int sourceEnd=bases.length, destinationEnd=length;
			long originalCopied=0, replacementCopied=0;
			for(int i=spans.size-4; i>=0; i-=4){
				final int from=spans.array[i], to=spans.array[i+1], offset=spans.array[i+2], len=spans.array[i+3];
				final int tailLength=sourceEnd-to; destinationEnd-=tailLength;
				assert(tailLength>=0 && destinationEnd>=len) : "Ordered disjoint spans and net indel accounting must leave room for each tail and replacement.";
				System.arraycopy(bases, to, output, destinationEnd, tailLength);
				if(quality!=null){System.arraycopy(quality, to, outputQuality, destinationEnd, tailLength);}
				destinationEnd-=len;
				System.arraycopy(replacementBases.array, offset, output, destinationEnd, len);
				if(quality!=null){System.arraycopy(replacementQualities.array, offset, outputQuality, destinationEnd, len);}
				originalCopied+=tailLength; replacementCopied+=len; sourceEnd=from;
			}
			destinationEnd-=sourceEnd;
			assert(destinationEnd==0) : "Backward interval splicing must consume exactly the final array length.";
			System.arraycopy(bases, 0, output, 0, sourceEnd);
			if(quality!=null){System.arraycopy(quality, 0, outputQuality, 0, sourceEnd);}
			originalCopied+=sourceEnd;
			assert(originalCopied+replacementCopied==length) : "Each final base must come from exactly one original span or staged replacement.";
			copiedOriginalBases=originalCopied; copiedReplacementBases=replacementCopied; materializedBases=length;
			substitutions=stagedSubs; insertions=stagedIns; deletions=stagedDels;
			read.bases=output; read.quality=outputQuality;
			return operations;
		}finally{clear();}
	}
	/** Drops pending spans and borrowed arrays while retaining reusable storage and
	 * the last result counters. Does not modify the borrowed Read or its arrays. */
	void clear(){
		spans.clear(); replacementBases.clear(); replacementQualities.clear();
		read=null; bases=null; quality=null; stagedSubs=stagedIns=stagedDels=0;
	}
	/** Returns pending patch intervals, each of which may contain several operations. */
	public int size(){return spans.size/4;}
	/** Returns the caller-reported operation total of the pending patches. */
	public int editCount(){return stagedSubs+stagedIns+stagedDels;}
	/** Clears result counters without changing the pending transaction. */
	private void clearCounters(){copiedOriginalBases=copiedReplacementBases=materializedBases=0; substitutions=insertions=deletions=0;}
	/** Base counts copied from originals/replacements and their final materialized total.
	 * Excludes quality copies; published on successful nonempty commit and preserved by clear(). */
	public long copiedOriginalBases, copiedReplacementBases, materializedBases;
	/** Local operation counts published by the last successful nonempty commit.
	 * Reset and active commit attempts clear them; clear() alone preserves them. */
	public int substitutions, insertions, deletions;
	/** Pending caller-reported operation counts; their sum is checked before each append. */
	private int stagedSubs, stagedIns, stagedDels;
	/** Borrowed transaction target, null when no batch is active. */
	private Read read;
	/** Original array identities retained for commit validation and copying. */
	private byte[] bases, quality;
	/** Ordered quadruples: original start, original end, staging offset, replacement length. */
	private final IntList spans=new IntList(32);
	/** Owned replacement pools, retained across transactions with matching offsets for qualities. */
	private final ByteBuilder replacementBases=new ByteBuilder(256), replacementQualities=new ByteBuilder(256);
}
