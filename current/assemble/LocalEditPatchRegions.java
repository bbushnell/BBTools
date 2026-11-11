package assemble;

import java.util.Arrays;
import stream.Read;
import structures.IntList;
import structures.LongList;

/** Conservative region discovery for the experimental native local-edit batch mode.
 * Tadpole uses this worker-owned selector with {@link LocalEditPatchCorrector}
 * when {@code fixindelsbatch=t}. That opt-in route requires dedicated
 * {@code fixindels=t} correction, dense scanning ({@code fixindelsstride=1}),
 * no pair lookahead or neural veto, and a patch bound of at least 2*K+5 bases.
 * Enumerates verified original-read proposals, then joins overlapping copy
 * intervals. It does not invent proposals for wide troughs the native kernel
 * cannot solve. All emitted coordinates refer to the original FORWARD read.
 * Singleton regions also use the local executor; this is not a hybrid fast path.
 * One worker owns discovery, geometry and returned-list scratch. Counts and input
 * arrays are borrowed and must stay stable; selection never installs read edits.
 * Calls are serial, and callers must consume or copy a result before reuse.
 * @author Fischl */
public final class LocalEditPatchRegions implements LocalEditCorrector.ProposalSink{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Allocates discovery and geometry scratch over borrowed immutable counts.
	 * Sparse discovery is permitted here, but cannot share preflight with the dense
	 * executor through correctSelected. The native batch route supplies stride one.
	 * @param k Exactly representable table word length, at least five
	 * @param lookup Borrowed nonnull count provider, stable during each call
	 * @param windows Positive odd probe-window count, no greater than K
	 * @param competition Whether indels compete with supported substitutions
	 * @param stride Positive original-profile sampling stride
	 * @param maxPatchBases Maximum emitted copy-window length, at least 2*K+5
	 */
	public LocalEditPatchRegions(final int k, final HomopolymerIndelProposal.CountLookup lookup,
			final int windows, final boolean competition, final int stride, final int maxPatchBases){
		this.k=k; this.lookup=lookup; this.windows=windows; this.competition=competition; this.stride=stride; this.maxPatchBases=maxPatchBases;
		regions=new Regions(k, maxPatchBases);
		discovery=new LocalEditCorrector(k, lookup, windows, competition, stride);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returned worker-owned list is valid only until the next selection. Input
	 * arrays stay untouched. LocalEditPatchCorrector.correctSelected owns a fresh
	 * selection and can reuse its preflight without accepting stale evidence.
	 * The positive limit is the original burden-guard threshold, not a maximum
	 * number of delivered proposals or emitted regions. Metadata/shape validation
	 * precedes scratch reset. Discovery/finish failure clears geometry; its query
	 * counts are retained by the finally, and partial guard state may remain.
	 * @param input Unpaired alignment-free Read with intact bases and matching or null qualities
	 * @param limit Positive whole-read edit budget used by the initial burden guard
	 * @return Borrowed quadruples [copyFrom, copyTo, coreFrom, coreTo), possibly empty
	 */
	public IntList select(final Read input, final int limit){
		if(limit<1){throw new IllegalArgumentException("Region discovery requires a positive whole-read edit budget.");}
		HomopolymerIndelEdit.requireEditable(input);
		if(input.bases==null || (input.quality!=null && input.quality.length!=input.length())){throw new IllegalArgumentException("Expected intact bases and matching or null qualities.");}
		regions.reset(input.length()); initialRejected=false; profileQueries=probeQueries=verificationQueries=0;
		try{
			discovery.collect(input, limit, this);
			initialRejected=discovery.callStatus==LocalEditCorrector.CallStatus.INITIAL_LIMIT;
			assert(!initialRejected || regions.order.size==0) : "LocalEditCorrector must finish its original burden guard before delivering any proposal.";
			return regions.finish();
		}catch(RuntimeException | Error failure){regions.reset(input.length()); throw failure;}
		finally{profileQueries=discovery.profileQueries; probeQueries=discovery.probeQueries; verificationQueries=discovery.verificationQueries;}
	}

	/** Consumes a verified proposal in the discovery engine's canonical orientation.
	 * Stores geometry only; the base identity and threshold are not retained. Reverse
	 * proposals are reflected into original-forward coordinates. Always returns true
	 * to continue discovery, including when the geometry is deferred at a read edge.
	 */
	@Override
	public boolean accept(final byte[] bases, final boolean reverse, final LocalSingleBaseEdit.Operation op,
			final int position, final byte base, final int contextFrom, final int contextTo, final int threshold){
		assert(bases.length==regions.length && threshold>=4) : "Region selection consumes verified proposals from one immutable original profile.";
		regions.add(op, position, contextFrom, contextTo, reverse); return true;
	}

	/** Tests exact setting equality, lookup identity and stride one before sharing
	 * selection's original guard with a dense executor. Equal count values alone
	 * cannot establish that two providers have the same lifetime or stable contents.
	 */
	boolean matchesDenseExecutor(final int k_, final HomopolymerIndelProposal.CountLookup lookup_, final int windows_, final boolean competition_, final int maxPatchBases_){
		// Lookup identity is deliberate: matching sizes do not prove the same
		// immutable count pool. Dense preflight parity also requires stride one.
		return k==k_ && lookup==lookup_ && windows==windows_ && competition==competition_ && stride==1 && maxPatchBases==maxPatchBases_;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Primitive geometry separated from sequence discovery for adversarial tests.
	 * Context endpoints and normalized positions need not arrive in sorted order.
	 * Reset opens one original-read collection; finish sorts/unions it once and
	 * closes appends until the next reset. Lists are owned reusable worker scratch.
	 */
	public static final class Regions{
		/** Configures pure geometry with K>=5 and space for two flanks and the preferred
		 * five-base core. This class itself does not construct words or check table layout.
		 */
		public Regions(final int k_, final int maxPatchBases_){
			if(k_<5 || maxPatchBases_<2L*k_+5){throw new IllegalArgumentException("Region scratch needs K context on both sides of a five-base core.");}
			k=k_; maxPatchBases=maxPatchBases_;
		}

		/** Opens a collection for a nonnegative original length, clearing lists/counters
		 * while retaining allocated capacities. Invalid length leaves earlier state intact.
		 */
		public void reset(final int length_){
			if(length_<0){throw new IllegalArgumentException("Original length must be nonnegative.");}
			length=length_; pending.clear(); order.clear(); output.clear(); edgeDeferred=sizeDeferred=merged=0; finished=false;
		}

		/** Adds a proposal window with K immutable flank bases and optional mutable slack.
		 * Position and dependency interval address the supplied canonical orientation;
		 * an insertion occupies a gap, while substitution/deletion consume one base.
		 * Edge-ineligible proposals increment edgeDeferred. Valid windows, including
		 * oversized ones, remain pending until finish unions connected copy intervals.
		 * The dependency interval must be nonempty and within the original length.
		 * @param op Nonnull proposed operation
		 * @param p Nonnegative base position below length, or insertion gap at most length
		 * @param contextFrom Inclusive dependency start
		 * @param contextTo Exclusive dependency end
		 * @param reverse Whether to reflect the resulting geometry into original-forward coordinates
		 */
		public void add(final LocalSingleBaseEdit.Operation op, final int p, final int contextFrom, final int contextTo, final boolean reverse){
			if(finished || op==null || p<0 || p>length || (op!=LocalSingleBaseEdit.Operation.INSERTION && p==length) ||
				contextFrom<0 || contextFrom>=contextTo || contextTo>length){throw new IllegalArgumentException("Expected a valid canonical proposal in an open original-coordinate collection.");}
			// Prefer two bases of mutable slack around the normalized event, but
			// trim only that OPTIONAL slack at true read boundaries. Keep K full
			// immutable bases on both sides and the complete event inside the core.
			// An insertion occupies a gap (including coreTo), not a consumed base.
			final boolean insertion=op==LocalSingleBaseEdit.Operation.INSERTION;
			long a=Math.max((long)k, (long)p-2), b=Math.min((long)length-k, (long)p+(insertion ? 2 : 3));
			if(a>b || p<a || (insertion ? p>b : p>=b)){edgeDeferred++; return;}
			long from=Math.min(a-k, contextFrom), to=Math.max(b+k, contextTo);
			if(from<0 || to>length){edgeDeferred++; return;}
			if(reverse){final long oldFrom=from, oldA=a; from=length-to; to=length-oldFrom; a=length-b; b=length-oldA;}
			final int index=pending.size;
			append(pending, (int)from, (int)to, (int)a, (int)b);
			order.add((from<<32)|(index&0xffffffffL));
		}

		/** Sorts windows by original start and unions overlapping or touching copies.
		 * Each component gets the smallest interval enclosing its cores, including any
		 * intervening gaps. Oversized copy components are omitted as a whole. Repeated
		 * calls return the same borrowed output without recounting. An exceptional
		 * first finish may leave partial closed state;
		 * select resets it before propagating that failure.
		 */
		public IntList finish(){
			if(finished){return output;}
			finished=true; Arrays.sort(order.array, 0, order.size);
			int from=-1, to=-1, a=-1, b=-1;
			for(int i=0; i<order.size; i++){
				final int n=(int)order.array[i], nextFrom=pending.array[n], nextTo=pending.array[n+1];
				if(from>=0 && nextFrom<=to){
					to=Math.max(to, nextTo); a=Math.min(a, pending.array[n+2]); b=Math.max(b, pending.array[n+3]); merged++;
				}else{
					if(from>=0){emit(from, to, a, b);}
					from=nextFrom; to=nextTo; a=pending.array[n+2]; b=pending.array[n+3];
				}
			}
			if(from>=0){emit(from, to, a, b);}
			assert(output.size%4==0) : "LocalEditPatchCorrector requires four original-coordinate fields per region.";
			return output;
		}

		/** Emits one disjoint component or counts one size deferral. Never splits a
		 * component whose windows may depend on each other's immutable evidence.
		 */
		private void emit(final int from, final int to, final int a, final int b){
			assert(from>=0 && to<=length && from<=a && a<=b && b<=to && a-from>=k && to-b>=k) : "Merged cores must retain the immutable flanks required by LocalEditPatchCorrector.";
			// Never split a connected component just to satisfy the size ceiling:
			// that could install two patches whose discovery evidence interferes.
			if(to-from>maxPatchBases){sizeDeferred++; return;}
			assert(output.size==0 || output.array[output.size-3]<from) : "Sorting and union must produce disjoint original copy intervals for one-splice materialization.";
			append(output, from, to, a, b);
		}

		/** Appends one copy/core quadruple; the caller owns validation and list lifetime. */
		private static void append(final IntList out, final int from, final int to, final int a, final int b){out.add(from); out.add(to); out.add(a); out.add(b);}

		/** Proposals deferred at edges, union components deferred for size, and absorbed
		 * windows. Merged counts include components later omitted for excessive size.
		 */
		public int edgeDeferred, sizeDeferred, merged;
		/** Original read length for both orientations and all emitted coordinates. */
		private int length;
		/** Whether finish has closed this collection to further proposals. */
		private boolean finished;
		/** Required immutable flank length and maximum emitted copy-window size. */
		private final int k, maxPatchBases;
		/** Pending and emitted original-forward copy/core quadruples; cleared on reset. */
		private final IntList pending=new IntList(), output=new IntList();
		/** Sort keys packing original copy start above the pending quadruple's offset. */
		private final LongList order=new LongList();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** True after discovery reports the original burden limit. Reset before discovery;
	 * earlier metadata/shape validation failures retain the previous value.
	 */
	public boolean initialRejected;
	/** Profile queries copied from discovery, including partial failed work. */
	public long profileQueries;
	/** Probe and verification queries copied from discovery, including partial failed work. */
	long probeQueries, verificationQueries;
	/** Worker-owned original-coordinate geometry and borrowed-result storage. */
	public final Regions regions;
	/** Worker-owned proposal engine; collect never installs edits into the input Read. */
	private final LocalEditCorrector discovery;
	/** Exact word length, probe-window count, profile stride and copy-window size ceiling. */
	private final int k, windows, stride, maxPatchBases;
	/** Whether supported substitutions still compete with indel proposals. */
	private final boolean competition;
	/** Borrowed immutable count provider; identity participates in dense-executor matching. */
	private final HomopolymerIndelProposal.CountLookup lookup;
}
