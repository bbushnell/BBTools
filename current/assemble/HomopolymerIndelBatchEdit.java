package assemble;

import shared.Shared;
import stream.Read;
import structures.IntList;

/** Transactional one-pass splice of approved disjoint-context original-coordinate edits.
 * This applies evidence supplied by a caller; it does not classify errors or count kmers.
 * @author Fischl
 */
final class HomopolymerIndelBatchEdit {
	HomopolymerIndelBatchEdit(){this(false);}
	HomopolymerIndelBatchEdit(final boolean singletons_){singletons=singletons_;}
	static final int WITHHELD=0, DIRECT=1, BASELINE=2, ISOLATED=3;

	/** Native acceptance order. DIRECT uses candidates; BASELINE/ISOLATED use out.
	 * The caller computed overlap while generating ordered original proposals.
	 * Rescue is never consulted after either old acceptance path succeeds. */
	int selectAccepted(final Read read, final IntList candidates, final boolean overlap,
			final int k, final int maxEdits, final boolean isolated, final IntList out){
		assert(candidates!=null && candidates.size%3==0 && maxEdits>0) :
			"Native batch selection requires complete original triples and a positive event-count budget.";
		if(candidates.size==0){return WITHHELD;}
		if(candidates.size/3<=maxEdits && !overlap){return DIRECT;}
		if(singletons && selectBaseline(candidates,3,k,maxEdits,out)){return BASELINE;}
		if(!isolated){return WITHHELD;}
		// A malformed proposal must not be hidden by discarding its overlap component.
		// Share apply's content/bounds/metadata validation, allowing context overlaps.
		validate(read,candidates,k,false);
		return HomopolymerIndelIsolation.select(candidates,k,maxEdits,out) ? ISOLATED : WITHHELD;
	}

	/** On combined-set rejection, copy the COMPLETE old non-singleton subset.
	 * Input stride3 is native triples; stride6 additionally carries dev evidence.
	 * No input mutation or arbitrary budget truncation. Empty/ineligible subsets
	 * leave output empty. Mechanical apply still validates complete Read contexts. */
	static boolean selectBaseline(final IntList candidates, final int stride, final int k, final int maxEdits, final IntList out){
		if(candidates==null || out==null || candidates==out || (stride!=3 && stride!=6) || candidates.size%stride!=0 || k<5 || maxEdits<1){
			throw new IllegalArgumentException("Baseline selection requires distinct lists, complete records, valid K and budget.");
		}
		out.clear(); long previousEnd=-1; boolean eligible=true;
		for(int i=0; i<candidates.size; i+=stride){
			final int start=candidates.get(i), end=candidates.get(i+1), delta=candidates.get(i+2);
			if(start<0 || end<=start || end-start>k-3 || (delta!=-1 && delta!=1) || (end-start==1 && delta!=1)){
				throw new IllegalArgumentException("Malformed approved candidate in baseline selection.");
			}
			if(end-start==1){continue;}
			final long from=(long)start-k-1, to=(long)end+k+1;
			if(from<0 || from<previousEnd){eligible=false;}
			previousEnd=to;
			out.add(start); out.add(end); out.add(delta);
		}
		if(!eligible || out.size==0 || out.size/3>maxEdits){out.clear(); return false;}
		assert(out.size%3==0 && out.size/3<=maxEdits) : "Fallback must be a complete budget-compliant original subset.";
		return true;
	}

	/** Validate before any allocation or mutation. Rescue may have overlapping
	 * contexts, but candidate maximal runs themselves must remain ordered/disjoint. */
	private long validate(final Read read, final IntList edits, final int k, final boolean disjoint){
		if(read==null || read.bases==null || edits==null || edits.size%3!=0 || k<5 ||
				(read.quality!=null && read.quality.length!=read.bases.length)){
			throw new IllegalArgumentException("Batch edit requires a read, matching qualities, valid K and complete triples.");
		}
		if(edits.size==0){return read.bases.length;}
		HomopolymerIndelEdit.requireEditable(read);
		final byte[] bases=read.bases;
		long length=bases.length, previousContextEnd=-1, previousRunEnd=-1;
		for(int i=0; i<edits.size; i+=3){
			final int start=edits.get(i), end=edits.get(i+1), delta=edits.get(i+2);
			final long from=(long)start-k-1, to=(long)end+k+1;
			if(start<0 || end>bases.length || (end-start<2 && !(singletons && end-start==1 && delta==1)) || end-start>k-3 || (delta!=-1 && delta!=1) ||
					from<0 || to>bases.length || start<previousRunEnd || (disjoint && from<previousContextEnd)){
				throw new IllegalArgumentException("Batch edits must have sorted, disjoint, complete original verification contexts: "+start+".."+end);
			}
			final byte b=bases[start];
			for(int p=(int)from; p<(int)to; p++){
				final byte x=bases[p];
				if((x!='A' && x!='C' && x!='G' && x!='T') || (p>=start && p<end && x!=b)){
					throw new IllegalArgumentException("Batch context must be defined and its run homogeneous.");
				}
			}
			if(bases[start-1]==b || bases[end]==b){throw new IllegalArgumentException("Batch edit run must be maximal.");}
			previousContextEnd=to; previousRunEnd=end; length+=delta;
		}
		assert(length>=0) : "Deleting at most one base per valid nonoverlapping run cannot exhaust the original sequence.";
		return length;
	}

	/** Triples (start,end,delta), ascending original coordinates; lists are never modified. */
	void apply(final Read read, final IntList edits, final int k){
		final long length=validate(read,edits,k,true);
		if(edits.size==0){return;}
		final byte[] bases=read.bases, qualities=read.quality;
		if(length<0 || length>Shared.MAX_ARRAY_LEN){throw new IllegalArgumentException("Batch output length exceeds supported array size.");}
		// Validate every candidate before allocating/installing either output array.
		final byte[] out=new byte[(int)length], outQ=qualities==null ? null : new byte[(int)length];
		int source=0, written=0;
		for(int i=0; i<edits.size; i+=3){
			final int start=edits.get(i), end=edits.get(i+1), delta=edits.get(i+2);
			final byte b=bases[start];
			final int pos=(b=='A' || b=='C') ? start : delta<0 ? end-1 : end;
			final int chunk=pos-source;
			assert(chunk>=0) : "Disjoint original contexts must produce monotone splice positions.";
			System.arraycopy(bases,source,out,written,chunk);
			if(outQ!=null){System.arraycopy(qualities,source,outQ,written,chunk);}
			written+=chunk;
			if(delta>0){out[written++]=b; source=pos;} // New quality slot remains raw0.
			else{source=pos+1;}
		}
		System.arraycopy(bases,source,out,written,bases.length-source);
		if(outQ!=null){System.arraycopy(qualities,source,outQ,written,bases.length-source);}
		assert(written+bases.length-source==out.length) : "All original suffix bytes and the declared edit deltas must explain final length.";
		read.bases=out; read.quality=outQ;
	}
	private final boolean singletons;
}
