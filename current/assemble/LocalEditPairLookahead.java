package assemble;

import java.util.Arrays;
import stream.Read;
import structures.IntList;
import ukmer.Kmer;

/** Experimental evidence-only two-edit witness for a bracketed unresolved trough.
 * Caller supplies a canonical whole-read orientation and its original profile.
 * Neither input is mutated. A witness is not permission to install a policy:
 * this component is intentionally separate from the validated one-edit caller.
 * Worker-local; temporary full-read copies are acceptable only in this prototype.
 * Reuse serially with stable input arrays, table counts and kmer layout settings.
 * Consume witness fields only after propose returns true; exceptions may leave
 * partial state. The successful witness is a new array retained until the next call.
 * Public witness and counter fields are maintained by this worker; callers observe
 * them after a call and must not assign them or mutate the retained witness array.
 * @author Fischl */
public final class LocalEditPairLookahead{
	/** Retains a count source and allocates one-window probes and a default single-edit caller.
	 * K must be at least5 and exactly representable. Neither caller uses a neural gate.
	 * Borrowed lookup keys must not be retained or have their represented words changed. */
	public LocalEditPairLookahead(final int k_, final HomopolymerIndelProposal.CountLookup source){
		if(k_<5 || source==null){throw new IllegalArgumentException("Pair witness requires K>=5 and immutable exact counts.");}
		k=k_;
		lookup=new HomopolymerIndelProposal.CountLookup(){
			@Override public int count(final Kmer key){
				queries++; final int value=source.count(key);
				if(value<-1){throw new IllegalStateException("Invalid count in pair witness: "+value);}
				return value;
			}
		};
		probe=new LocalEditKmerProbe(k, lookup); second=new LocalEditCorrector(k, lookup); key=new Kmer(k);
	}
	/** Returns the first supported two-edit witness under the prioritized search.
	 * First-edit substitutions are searched before indels; a substitution witness
	 * returns without comparing indel-first outcomes. Distinct verified outcomes
	 * within the searched family make that trough ambiguous and skip it. The second
	 * edit uses the default substitution-first single-edit corrector, so uniqueness
	 * is scoped to this search, not all possible two-edit outcomes. Retains the
	 * first edit and complete witness; no recursion or unbounded path walk. */
	public boolean propose(final byte[] bases, final IntList counts){
		return propose(bases, counts, null);
	}
	/** Exclusions are increasing depth-start coordinates from this same original profile.
	 * Resets witness fields and counters before validating inputs. Profile size must be
	 * max(0,length-K+1), with depths at least -1. Edge, undefined and width>=2K regions
	 * are skipped. Excluded starts must be distinct and in range; null means none.
	 * A normal false return leaves no witness; no input bases or depths are changed. */
	public boolean propose(final byte[] bases, final IntList counts, final IntList excludedStarts){
		clear(); queries=0; regions=candidates=localSecondEdits=verifiedPairs=0;
		if(bases==null || counts==null || counts.size!=Math.max(0, bases.length-k+1)){
			throw new IllegalArgumentException("Pair witness needs the unchanged original kmer-start profile.");
		}
		for(int j=0; j<counts.size; j++){if(counts.get(j)<-1){throw new IllegalArgumentException("Invalid original depth: "+counts.get(j));}}
		if(excludedStarts!=null){for(int i=0; i<excludedStarts.size; i++){
			assert(excludedStarts.get(i)>=0 && excludedStarts.get(i)<counts.size && (i==0 || excludedStarts.get(i-1)<excludedStarts.get(i))) :
				"Exclusion merge requires distinct increasing starts in the unchanged canonical profile.";
		}}
		int excludedIndex=0;
		for(int cursor=0; cursor<counts.size;){
			if(counts.get(cursor)>=3){cursor++; continue;}
			final int a=cursor++; int originalMax=Math.max(0, counts.get(a));
			while(cursor<counts.size && counts.get(cursor)<3){originalMax=Math.max(originalMax, counts.get(cursor)); cursor++;}
			final int end=cursor, width=end-a;
			if(excludedStarts!=null){
				while(excludedIndex<excludedStarts.size && excludedStarts.get(excludedIndex)<a){excludedIndex++;}
				if(excludedIndex<excludedStarts.size && excludedStarts.get(excludedIndex)==a){continue;}
			}
			// Repeat-equivalent missing bases can make a two-error trough <=K.
			// Keep the finite upper bound; full two-edit evidence decides eligibility.
			if(a==0 || end==counts.size || width>=2L*k){continue;}
			final int from=a-1, to=end+k, p=a+k-1;
			boolean called=true;
			for(int i=from; i<to; i++){final byte b=bases[i]; if(b!='A' && b!='C' && b!='G' && b!='T'){called=false; break;}}
			if(!called){continue;}
			regions++; clear(); final int threshold=Math.max(4, 4*originalMax+1);
			probe.substitutions(bases, a, p);
			for(int b=0; b<4; b++){if(probe.substitutionDepth[b]>=threshold){consider(bases, LocalSingleBaseEdit.Operation.SUBSTITUTION, p, ALPHABET[b], from, to, threshold);}}
			if(ambiguous){clear(); continue;}
			if(witness!=null){return true;}
			probe.indelsAt(bases, a, p);
			if(probe.deletionDepth>=threshold){consider(bases, LocalSingleBaseEdit.Operation.DELETION, p, (byte)0, from, to, threshold);}
			for(int b=0; b<4; b++){if(probe.insertionDepth[b]>=threshold){consider(bases, LocalSingleBaseEdit.Operation.INSERTION, p, ALPHABET[b], from, to, threshold);}}
			// The first bad word can end just AFTER an extra repeated base.
			// Test the adjacent deletion uniformly; equivalent outcomes collapse below.
			probe.indelsAt(bases, a, p-1);
			if(probe.deletionDepth>=threshold){consider(bases, LocalSingleBaseEdit.Operation.DELETION, p-1, (byte)0, from, to, threshold);}
			if(!ambiguous && witness!=null){return true;}
			clear();
		}
		return false;
	}
	/** Applies the first edit to a shadow read, then asks the single-edit caller for one more.
	 * Rejects a second edit outside the local context, unsupported final windows and
	 * endpoints reachable in at most one edit. Distinct surviving sequences set ambiguity;
	 * equivalent sequences retain the first witness and its operation coordinates. */
	private void consider(final byte[] bases, final LocalSingleBaseEdit.Operation op, final int p, final byte b, final int from, final int to, final int threshold){
		candidates++;
		final Read shadow=new Read(bases, null, "pair-witness", 0, false);
		if(op==LocalSingleBaseEdit.Operation.DELETION){editor.apply(shadow, op, p);}else{editor.apply(shadow, op, p, b);}
		final int firstDelta=shadow.length()-bases.length;
		if(second.correctOne(shadow)!=1){return;}
		final int q=second.lastPosition;
		assert(q>=0 && (second.lastOperation==LocalSingleBaseEdit.Operation.INSERTION ? q<=bases.length+firstDelta : q<bases.length+firstDelta)) :
			"Second edit position is in the post-first-edit input, with insertion allowed at its final boundary.";
		if(q<from || q>=to+firstDelta){return;}
		localSecondEdits++;
		final int finalTo=to+shadow.length()-bases.length;
		final long bounds=verificationBounds(k, shadow.length(), op, p, second.lastOperation, q, from, finalTo);
		if(!verify(shadow.bases, (int)(bounds>>>32), (int)bounds, threshold)){return;}
		// Two operations are not necessarily two edits: MG2 read5030 reached a
		// one-deletion endpoint via SUB+DEL. Never commit that neutral detour.
		if(withinOneEdit(bases, shadow.bases)){return;}
		verifiedPairs++;
		if(witness==null){
			firstOperation=op; firstPosition=p; firstBase=b; witness=shadow.bases;
			secondOperation=second.lastOperation; secondPositionAfterFirst=q; contextStart=from; contextEnd=to; witnessContextEnd=finalTo;
		}else if(!Arrays.equals(witness, shadow.bases)){ambiguous=true;}
	}
	/** Union of context and both affected base intervals in final-read coordinates.
	 * p addresses the original read; q addresses the read after the first edit.
	 * length and the half-open context address the final read. Bounds enclose affected
	 * windows, not necessarily the smallest possible interval. Returns start in high32
	 * bits, exclusive end in low32; no per-candidate array. */
	public static long verificationBounds(final int k, final int length, final LocalSingleBaseEdit.Operation first, final int p,
		final LocalSingleBaseEdit.Operation second, final int q, final int contextFrom, final int contextTo){
		assert(k>0 && length>=k && contextFrom>=0 && contextTo<=length && contextTo-contextFrom>=k) :
			"Pair verification must retain the existing complete supported-flank context in the final read.";
		final int delta=second==LocalSingleBaseEdit.Operation.INSERTION ? 1 : second==LocalSingleBaseEdit.Operation.DELETION ? -1 : 0;
		//After the first edit, S/I affect starts p-K+1..p, D p-K+1..p-1.
		//Map the enclosing base interval through edit2; crossing windows created
		//by edit2 are separately covered by its own final affected interval.
		final int firstFrom=p-k+1, firstTo=p+k-(first==LocalSingleBaseEdit.Operation.DELETION ? 1 : 0);
		final int mappedFrom=firstFrom+(q<firstFrom ? delta : 0), mappedTo=firstTo+(q<firstTo ? delta : 0);
		final int secondFrom=q-k+1, secondTo=q+k-(second==LocalSingleBaseEdit.Operation.DELETION ? 1 : 0);
		final int from=Math.max(0, Math.min(contextFrom, Math.min(mappedFrom, secondFrom)));
		final int to=Math.min(length, Math.max(contextTo, Math.max(mappedTo, secondTo)));
		assert(from<=contextFrom && to>=contextTo) : "Expanded two-edit verification cannot omit any previously required context window.";
		return ((long)from<<32)|(to&0xffffffffL);
	}
	/** Exact unit-cost distance {@code <=1}, without allocations or an alignment matrix. */
	public static boolean withinOneEdit(final byte[] a, final byte[] b){
		assert(a!=null && b!=null) : "Witness minimality compares complete original and final sequences, not partial alignment spans.";
		if(a.length==b.length){
			int mismatches=0;
			for(int i=0; i<a.length; i++){if(a[i]!=b[i] && ++mismatches>1){return false;}}
			return true;
		}
		final byte[] longer=a.length>b.length ? a : b, shorter=a.length>b.length ? b : a;
		if(longer.length-shorter.length!=1){return false;}
		int i=0; while(i<shorter.length && longer[i]==shorter[i]){i++;}
		// The sole length-changing edit skips one longer-array base at the first
		// mismatch (or at the end); every remaining aligned base must agree.
		for(; i<shorter.length; i++){if(shorter[i]!=longer[i+1]){return false;}}
		return true;
	}
	/** Requires every complete window in the half-open final context to meet threshold.
	 * Undefined bytes reject the candidate; the reusable key and lookup counter change. */
	private boolean verify(final byte[] bases, final int from, final int to, final int threshold){
		assert(from>=0 && to<=bases.length && to-from>=k) :
			"Both edits are inside the original supported-flank context; mapped context must remain complete.";
		key.clearFast(); int measured=0;
		for(int i=from; i<to; i++){
			final byte b=bases[i]; if(b!='A' && b!='C' && b!='G' && b!='T'){return false;}
			key.addRight(bases[i]);
			if(key.len()>=k){measured++; if(lookup.count(key)<threshold){return false;}}
		}
		assert(measured==to-from-k+1) : "All kmers in the two-edit witness context must be verified, not just the boundary probe.";
		return true;
	}
	/** Clears retained witness/ambiguity fields without resetting accumulated call counters. */
	private void clear(){
		witness=null; firstOperation=secondOperation=null; firstPosition=secondPositionAfterFirst=contextStart=contextEnd=witnessContextEnd=-1; firstBase=0; ambiguous=false;
	}
	/** Complete two-edit result, retained by this worker; consume only after a true return. */
	public byte[] witness;
	/** Operations producing the retained witness, not a policy to install both at once. */
	public LocalSingleBaseEdit.Operation firstOperation, secondOperation;
	/** First edit and contextEnd are original coordinates; secondPositionAfterFirst
	 * addresses the temporary input to the second edit. witnessContextEnd includes both deltas. */
	public int firstPosition, secondPositionAfterFirst, contextStart, contextEnd, witnessContextEnd;
	/** First replacement/inserted allele in the supplied orientation; unused for deletion. */
	public byte firstBase;
	/** Count-source invocations this call, including an invocation that throws. */
	public long queries;
	/** Eligible troughs, first-edit attempts, local second edits and verified two-edit witnesses.
	 * These are search counters, not unique outcomes or biological corrections. */
	public int regions, candidates, localSecondEdits, verifiedPairs;
	/** Whether the current searched family produced different final sequences. */
	private boolean ambiguous;
	/** Exactly represented table span shared by the scratch helpers. */
	private final int k;
	/** Retained counting adapter over the borrowed source; also maintains queries. */
	private final HomopolymerIndelProposal.CountLookup lookup;
	/** One-window first-edit candidate probe. */
	private final LocalEditKmerProbe probe;
	/** Reused substitution-first caller, invoked without recursive pair lookahead. */
	private final LocalEditCorrector second;
	/** Transactional editor for the shadow read's first edit. */
	private final LocalSingleBaseEdit editor=new LocalSingleBaseEdit();
	/** Rolling word reused to verify the complete final context. */
	private final Kmer key;
	private static final byte[] ALPHABET={'A', 'C', 'G', 'T'};
}
