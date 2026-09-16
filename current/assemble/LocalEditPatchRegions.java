package assemble;

import java.util.Arrays;
import stream.Read;
import structures.IntList;
import structures.LongList;

/** Unrouted conservative region discovery for LocalEditPatchCorrector.
 * Enumerates verified original-read proposals, then joins overlapping copy
 * intervals. It does not invent proposals for wide troughs the native kernel
 * cannot solve. All emitted coordinates refer to the original FORWARD read.
 * Singleton regions also use the local executor; this is not a hybrid fast path.
 * @author Fischl */
final class LocalEditPatchRegions implements LocalEditCorrector.ProposalSink {

	LocalEditPatchRegions(final int k,final HomopolymerIndelProposal.CountLookup lookup,
		final int windows,final boolean competition,final int stride,final int maxPatchBases){
		this.k=k;this.lookup=lookup;this.windows=windows;this.competition=competition;this.stride=stride;this.maxPatchBases=maxPatchBases;
		regions=new Regions(k,maxPatchBases);
		discovery=new LocalEditCorrector(k,lookup,windows,competition,stride);
	}
	/** Returned worker-owned list is valid only until the next selection. Input
	 * arrays stay untouched. LocalEditPatchCorrector.correctSelected owns a fresh
	 * selection and can reuse its preflight without accepting stale evidence. */
	IntList select(final Read input,final int limit){
		if(limit<1){throw new IllegalArgumentException("Region discovery requires a positive whole-read edit budget.");}
		HomopolymerIndelEdit.requireEditable(input);
		if(input.bases==null || (input.quality!=null && input.quality.length!=input.length())){throw new IllegalArgumentException("Expected intact bases and matching or null qualities.");}
		regions.reset(input.length());initialRejected=false;profileQueries=probeQueries=verificationQueries=0;
		try{
			discovery.collect(input,limit,this);
			initialRejected=discovery.callStatus==LocalEditCorrector.CallStatus.INITIAL_LIMIT;
			assert(!initialRejected || regions.order.size==0) : "LocalEditCorrector must finish its original burden guard before delivering any proposal.";
			return regions.finish();
		}catch(RuntimeException | Error failure){regions.reset(input.length());throw failure;}
		finally{profileQueries=discovery.profileQueries;probeQueries=discovery.probeQueries;verificationQueries=discovery.verificationQueries;}
	}
	@Override public boolean accept(final byte[] bases,final boolean reverse,final LocalSingleBaseEdit.Operation op,
		final int position,final byte base,final int contextFrom,final int contextTo,final int threshold){
		assert(bases.length==regions.length && threshold>=4) : "Region selection consumes verified proposals from one immutable original profile.";
		regions.add(op,position,contextFrom,contextTo,reverse);return true;
	}
	boolean matchesDenseExecutor(final int k_,final HomopolymerIndelProposal.CountLookup lookup_,final int windows_,final boolean competition_,final int maxPatchBases_){
		// Lookup identity is deliberate: matching sizes do not prove the same
		// immutable count pool. Dense preflight parity also requires stride one.
		return k==k_ && lookup==lookup_ && windows==windows_ && competition==competition_ && stride==1 && maxPatchBases==maxPatchBases_;
	}
	/** Primitive geometry separated from sequence discovery for adversarial tests.
	 * Context endpoints and normalized positions need not arrive in sorted order. */
	static final class Regions {
		Regions(final int k_,final int maxPatchBases_){
			if(k_<5 || maxPatchBases_<2L*k_+5){throw new IllegalArgumentException("Region scratch needs K context on both sides of a five-base core.");}
			k=k_;maxPatchBases=maxPatchBases_;
		}
		void reset(final int length_){
			if(length_<0){throw new IllegalArgumentException("Original length must be nonnegative.");}
			length=length_;pending.clear();order.clear();output.clear();edgeDeferred=sizeDeferred=merged=0;finished=false;
		}
		void add(final LocalSingleBaseEdit.Operation op,final int p,final int contextFrom,final int contextTo,final boolean reverse){
			if(finished || op==null || p<0 || p>length || (op!=LocalSingleBaseEdit.Operation.INSERTION && p==length) ||
				contextFrom<0 || contextFrom>=contextTo || contextTo>length){throw new IllegalArgumentException("Expected a valid canonical proposal in an open original-coordinate collection.");}
			// Prefer two bases of mutable slack around the normalized event, but
			// trim only that OPTIONAL slack at true read boundaries. Keep K full
			// immutable bases on both sides and the complete event inside the core.
			// An insertion occupies a gap (including coreTo), not a consumed base.
			final boolean insertion=op==LocalSingleBaseEdit.Operation.INSERTION;
			long a=Math.max((long)k,(long)p-2),b=Math.min((long)length-k,(long)p+(insertion ? 2 : 3));
			if(a>b || p<a || (insertion ? p>b : p>=b)){edgeDeferred++;return;}
			long from=Math.min(a-k,contextFrom),to=Math.max(b+k,contextTo);
			if(from<0 || to>length){edgeDeferred++;return;}
			if(reverse){final long oldFrom=from,oldA=a;from=length-to;to=length-oldFrom;a=length-b;b=length-oldA;}
			final int index=pending.size;
			append(pending,(int)from,(int)to,(int)a,(int)b);
			order.add((from<<32)|(index&0xffffffffL));
		}
		IntList finish(){
			if(finished){return output;}
			finished=true;Arrays.sort(order.array,0,order.size);
			int from=-1,to=-1,a=-1,b=-1;
			for(int i=0;i<order.size;i++){
				final int n=(int)order.array[i],nextFrom=pending.array[n],nextTo=pending.array[n+1];
				if(from>=0 && nextFrom<=to){
					to=Math.max(to,nextTo);a=Math.min(a,pending.array[n+2]);b=Math.max(b,pending.array[n+3]);merged++;
				}else{
					if(from>=0){emit(from,to,a,b);}
					from=nextFrom;to=nextTo;a=pending.array[n+2];b=pending.array[n+3];
				}
			}
			if(from>=0){emit(from,to,a,b);}
			assert(output.size%4==0) : "LocalEditPatchCorrector requires four original-coordinate fields per region.";
			return output;
		}
		private void emit(final int from,final int to,final int a,final int b){
			assert(from>=0 && to<=length && from<=a && a<=b && b<=to && a-from>=k && to-b>=k) : "Merged cores must retain the immutable flanks required by LocalEditPatchCorrector.";
			// Never split a connected component just to satisfy the size ceiling:
			// that could install two patches whose discovery evidence interferes.
			if(to-from>maxPatchBases){sizeDeferred++;return;}
			assert(output.size==0 || output.array[output.size-3]<from) : "Sorting and union must produce disjoint original copy intervals for one-splice materialization.";
			append(output,from,to,a,b);
		}
		private static void append(final IntList out,final int from,final int to,final int a,final int b){out.add(from);out.add(to);out.add(a);out.add(b);}
		int edgeDeferred,sizeDeferred,merged;
		private int length;private boolean finished;
		private final int k,maxPatchBases;
		private final IntList pending=new IntList(),output=new IntList();
		private final LongList order=new LongList();
	}
	boolean initialRejected;
	long profileQueries,probeQueries,verificationQueries;
	final Regions regions;
	private final LocalEditCorrector discovery;
	private final int k,windows,stride,maxPatchBases;
	private final boolean competition;
	private final HomopolymerIndelProposal.CountLookup lookup;
}
