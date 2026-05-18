package ddl;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.PriorityQueue;

/**
 * Fixed-size min-heap retaining the top-N best DDLComparison results.
 * Worst element sits at the root for O(1) eviction checks.
 * Copies on offer — never stores the caller's mutable working object.
 *
 * Modeled after clade.ComparisonHeap.
 *
 * @author Noire, Brian Bushnell
 * @date May 17, 2026
 */
public class DDLComparisonHeap {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public DDLComparisonHeap(int maxSize_){
		maxSize=maxSize_;
		heap=new PriorityQueue<>(Math.max(maxSize, 1), WORST_FIRST);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Offers a comparison to the heap.
	 * If the heap has room, copies and adds.
	 * If full, replaces the worst element only if the candidate is better.
	 * The caller's object is never retained — always copied via setFrom. */
	public void offer(DDLComparison candidate){
		if(heap.size()<maxSize){
			DDLComparison copy=new DDLComparison();
			copy.setFrom(candidate);
			heap.add(copy);
		}else{
			DDLComparison worst=heap.peek();
			if(worst!=null && candidate.compareTo(worst)<0){
				heap.poll();
				worst.setFrom(candidate);
				heap.add(worst);
			}
		}
	}

	/** Drains the heap into a sorted list (best first). */
	public ArrayList<DDLComparison> toList(){
		ArrayList<DDLComparison> list=new ArrayList<>(heap.size());
		while(!heap.isEmpty()){
			list.add(heap.poll());
		}
		for(int i=0, j=list.size()-1; i<j; i++, j--){
			DDLComparison tmp=list.get(i);
			list.set(i, list.get(j));
			list.set(j, tmp);
		}
		return list;
	}

	/** Returns the worst (lowest-scoring) element without removing it. */
	public DDLComparison worst(){return heap.peek();}

	public int size(){return heap.size();}

	public void clear(){heap.clear();}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final PriorityQueue<DDLComparison> heap;
	public final int maxSize;

	/** Comparator that puts the worst comparison at the root.
	 * Natural ordering is best-first (compareTo returns -1 for better),
	 * so we reverse it to put worst at root for eviction. */
	private static final Comparator<DDLComparison> WORST_FIRST=
		(a, b)->b.compareTo(a);
}
