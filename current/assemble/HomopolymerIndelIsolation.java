package assemble;

import structures.IntList;

/** Select size-one components of the original verification-context overlap graph.
 * Geometry only: the caller supplies approved proposals and must still validate
 * the actual read with HomopolymerIndelBatchEdit before applying any changes.
 * Native use is opt-in and only after existing acceptance paths reject.
 * @author Fischl
 */
final class HomopolymerIndelIsolation {

	/** Copy the complete isolated subset of sorted (start,end,delta) triples.
	 * False means empty or over budget; output is empty in either case. No input
	 * mutation, allocation, ranking, truncation, or net-delta budgeting. Malformed
	 * candidates throw, including discarded candidates. */
	static boolean select(final IntList candidates, final int k, final int maxEdits, final IntList out){
		if(candidates==null || out==null || candidates==out || candidates.size%3!=0 || k<5 || maxEdits<1){
			throw new IllegalArgumentException("Isolation needs distinct lists, complete triples, K>=5 and a positive edit budget.");
		}
		out.clear();
		long previousEnd=-1;
		for(int i=0; i<candidates.size; i+=3){
			final int start=candidates.get(i), end=candidates.get(i+1), delta=candidates.get(i+2);
			if((long)start-k-1<0 || end<=start || end-start>k-3 || start<previousEnd ||
					(delta!=-1 && delta!=1) || (end-start==1 && delta!=1)){
				throw new IllegalArgumentException("Isolation needs sorted disjoint run intervals, valid lengths/deltas and complete left flanks: "+start+".."+end);
			}
			previousEnd=end;
		}
		// All run intervals are ordered and disjoint; expanding each by the same
		// K+1 preserves both endpoint orders. Thus neighbors witness every overlap.
		for(int i=0; i<candidates.size; i+=3){
			final long from=(long)candidates.get(i)-k-1, to=(long)candidates.get(i+1)+k+1;
			final boolean left=i>0 && from<(long)candidates.get(i-2)+k+1;
			final boolean right=i+3<candidates.size && (long)candidates.get(i+3)-k-1<to;
			if(!left && !right){
				out.add(candidates.get(i)); out.add(candidates.get(i+1)); out.add(candidates.get(i+2));
			}
		}
		if(out.size==0 || out.size/3>maxEdits){out.clear(); return false;}
		assert(out.size%3==0 && out.size<=candidates.size && out.size/3<=maxEdits) :
			"Isolation must preserve whole input triples and obey the count-of-edits budget, never a net-delta budget.";
		return true;
	}
}
