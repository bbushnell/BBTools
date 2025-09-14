package bin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.PriorityQueue;

/**
 * Maintains a heap of the top N Comparisons encountered.
 * Useful for finding the best N matches rather than just the single best match.
 * 
 * @author Brian Bushnell
 * @contributor Isla (Highly-customized Claude instance)
 * @date April 19, 2024
 */
public class ComparisonHeap {
    
    /**
     * Comparator that reverses the natural ordering of Comparisons.
     * This keeps the worst comparison at the top of our heap.
     */
    private static class WorstFirstComparator implements Comparator<Comparison> {
        @Override
        public int compare(Comparison c1, Comparison c2) {
            // Invert the comparison result to get worst-first ordering
            return -c1.compareTo(c2);
        }
    }
    
    /**
     * Creates a new ComparisonHeap with the specified maximum size.
     * 
     * @param maxSize Maximum number of comparisons to keep
     */
    public ComparisonHeap(int maxSize) {
        this.maxSize = maxSize;
        // Create a heap ordered to keep the worst comparison at the top
        this.heap = new PriorityQueue<>(maxSize, new WorstFirstComparator());
    }
    
    /**
     * Offers a comparison to the heap. The comparison will be added if:
     * 1. The heap is not yet full, or
     * 2. The comparison is better than the worst one in the heap
     * 
     * @param comp The comparison to offer (will not be modified or stored directly)
     * @return true if the comparison was added, false otherwise
     */
    public boolean offer(Comparison comp) {
        if (heap.size() < maxSize) {
            // If the heap isn't full yet, create a new Comparison and add it
            Comparison newComp = new Comparison();
            newComp.setFrom(comp);
            heap.add(newComp);
            return true;
        } else {
            // The heap is full - compare with the worst element
            Comparison worst = heap.peek();
            // If comp is better than worst, comp.compareTo(worst) returns -1
            if (comp.compareTo(worst) < 0) {
                // Better than the worst - remove worst, update it with new values, and add back
                worst = heap.poll();
                worst.setFrom(comp);
                heap.add(worst);
                return true;
            }
        }
        return false;
    }
    
    /**
     * Returns a sorted list of the top comparisons (best first)
     */
    public ArrayList<Comparison> toList() {
    	ArrayList<Comparison> result = new ArrayList<Comparison>(heap);
        // Sort using natural ordering (compareTo) which puts best first
        Collections.sort(result);
        return result;
    }
    
    /**
     * Returns the number of comparisons currently in the heap
     */
    public int size() {
        return heap.size();
    }
    
    /**
     * Clears all comparisons from the heap
     */
    public void clear() {
        heap.clear();
    }
    
//    /**
//     * Returns the best comparison in the heap, or null if empty
//     */
//    public Comparison getBest() {
//        if (heap.isEmpty()) return null;
//        
//        // Find the best element in the heap by comparing all elements
//        Comparison best = null;
//        
//        for (Comparison comp : heap) {
//            if (best == null || comp.compareTo(best) < 0) {
//                best = comp;
//            }
//        }
//        
//        return best;
//    }
    
    public Comparison worst() {
    	return heap.peek();
    }
    
    private final PriorityQueue<Comparison> heap;
    private final int maxSize;
}