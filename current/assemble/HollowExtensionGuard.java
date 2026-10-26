package assemble;

import dna.AminoAcid;
import ukmer.HollowKmerLookup;
import ukmer.KmerTableSetU;

/** Worker-local, allocation-free alternate-span branch checks for read extension.
 * Reset for each oriented read/call, then query before appending each candidate.
 * Counts remain immutable; this guard never corrects bases or rewinds a seed.
 * @author Brian Bushnell, Fischl */
public final class HollowExtensionGuard {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Uses the same loaded Hollow counts as the nominal extension engine. */
	public HollowExtensionGuard(final KmerTableSetU source, final int ratio, final int floor){
		this(new HollowKmerLookup(source), source.kbig, ratio, floor);
	}

	/** Creates reusable state; an explicit lookup also permits small deterministic fixtures. */
	public HollowExtensionGuard(final HollowKmerLookup lookup_, final int k_, final int ratio_, final int floor_){
		if(lookup_==null || k_<67 || ratio_<1 || floor_<1){
			throw new IllegalArgumentException("Hollow extension guard requires loaded counts, K>=67, ratio>=1 and floor>=1.");
		}
		lookup=lookup_;
		k=k_;
		ratio=ratio_;
		floor=floor_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts at a prefix whose terminal K bases are defined, as required by rightmostKmer.
	 * The extra K+1 context is optional: a missing base or preceding N disables
	 * only that spacing for the first step. After one accepted base it is defined.
	 * This one-time scan is O(K); all subsequent sequential checks are O(1) in K. */
	public void reset(final byte[] bases, final int length){
		if(bases==null || length<k || length>bases.length){
			throw new IllegalArgumentException("Hollow guard reset needs a complete terminal K window.");
		}
		for(int i=length-k; i<length; i++){
			if(AminoAcid.baseToNumber[bases[i]]<0){
				throw new IllegalArgumentException("Hollow guard requires a defined terminal K window; ambiguous position="+i);
			}
		}
		clearWarning();
		end=length;
		leftMinus=arm(bases, length-k+1);
		right=arm(bases, length-32);
		plusValid=length>k && AminoAcid.baseToNumber[bases[length-k-1]]>=0;
		// If context is missing, the synthetic leading A is never queried; it
		// shifts out when the first accepted base makes the K+1 window available.
		leftPlus=arm(bases, length-k)>>>2;
		if(plusValid){leftPlus|=(long)AminoAcid.baseToNumber[bases[length-k-1]]<<62;}
		assert(end-k+33<end) : "K>=67 keeps both left-arm ingress bases inside the existing prefix, not the candidate.";
	}

	/** Tests one proposed base before the caller commits it to the sequence.
	 * Repeated queries at the same endpoint do not advance state. A following
	 * endpoint may be exactly one base longer; its actual appended base advances
	 * the anchors, even when the ByteBuilder has reallocated its backing array.
	 * Any new read, orientation or discontinuous endpoint requires reset(). */
	public boolean reject(final byte[] bases, final int length, final int candidate){
		return reject(bases, length, candidate, true);
	}

	/** Clears call-local history even when the extension engine returns before reset(). */
	public void clearWarning(){firstWarning=-1;}

	/** Remembers the earliest hard contradiction, but does not stop or alter the path.
	 * Once recorded, later count queries are unnecessary until the next call. */
	public void observe(final byte[] bases, final int length, final int candidate){
		if(firstWarning<0 && reject(bases, length, candidate, false)){firstWarning=length;}
	}

	/** A completed request is retained; only a later stall activates the checkpoint.
	 * A candidate that was inspected but never appended cannot itself be removed. */
	public int retainedLength(final int length, final int limit){
		assert(length>=0 && limit>=length) : "The extension wrapper bounds retained sequence by its requested limit.";
		return length<limit && firstWarning>=0 && firstWarning<length ? firstWarning : length;
	}

	/** Shares rolling geometry while keeping delayed hard conflicts distinct from ambiguity stops. */
	private boolean reject(final byte[] bases, final int length, final int candidate, final boolean ambiguity){
		assert(candidate>=0 && candidate<4) : "The nominal extension engine chooses an ACGT index, never an ambiguous candidate.";
		if(end<0 || bases==null || length>bases.length || (length!=end && length!=end+1)){
			throw new IllegalArgumentException("Reset Hollow guard before changing reads or skipping endpoints; previous="+end+", next="+length);
		}
		if(length==end+1){advance(bases);}
		if(plusValid && rejectSpan(leftPlus, bases[end-k+31], candidate, ambiguity)){return true;}
		return rejectSpan(leftMinus, bases[end-k+33], candidate, ambiguity);
	}

	/** Rolls from the previous prefix using only bases already accepted by the caller. */
	private void advance(final byte[] bases){
		final int accepted=AminoAcid.baseToNumber[bases[end]];
		assert(accepted>=0) : "Sequential extension may append only a defined base; reset is required after any other sequence edit.";
		leftPlus=(leftPlus<<2)|AminoAcid.baseToNumber[bases[end-k+31]];
		leftMinus=(leftMinus<<2)|AminoAcid.baseToNumber[bases[end-k+33]];
		right=(right<<2)|accepted;
		plusValid=true;
		end++;
	}

	/** Reproduces the diagnostic hard-veto or immediate-ambiguity rule for one spacing. */
	private boolean rejectSpan(final long left, final byte ingress, final int candidate, final boolean ambiguity){
		final int incoming=AminoAcid.baseToNumber[ingress];
		assert(incoming>=0) : "The left-arm ingress lies within the already validated terminal K window.";
		final int context=lookup.count(left, right);
		if(context<floor){return false;}
		final long shiftedLeft=(left<<2)|incoming, shiftedRight=right<<2;
		int chosen=0, other=0;
		for(int base=0; base<4; base++){
			final int count=lookup.count(shiftedLeft, shiftedRight|base);
			if(base==candidate){chosen=count;}
			else{other=Math.max(other, count);}
		}
		// Same predicates as KmerWindowDepth.rejectPathStep/ambiguousPathStep.
		return other>=floor && (other>=(long)ratio*chosen ||
				(ambiguity && chosen>=floor && chosen<(long)ratio*other && other<(long)ratio*chosen));
	}

	/** Packs one defined 32-base arm once at reset, never allocating a Kmer. */
	private static long arm(final byte[] bases, final int start){
		assert(start>=0 && start+32<=bases.length) : "Endpoint encoding requires a complete 32-base arm inside the prefix.";
		long value=0;
		for(int i=start; i<start+32; i++){value=(value<<2)|AminoAcid.baseToNumber[bases[i]];}
		return value;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final HollowKmerLookup lookup;
	private final int k, ratio, floor;
	private long leftPlus, leftMinus, right;
	private int end=-1;
	private int firstWarning=-1;
	private boolean plusValid;
}
