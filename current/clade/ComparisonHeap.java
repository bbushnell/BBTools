package clade;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.PriorityQueue;

/**
 * Maintains a heap of the top N Comparison objects (best-first retrieval).
 * Keeps worst element at the root to allow efficient replacement when better comparisons arrive.
 * Useful for retaining top matches rather than a single best match.
 * @author Brian Bushnell
 * @contributor Isla
 * @date April 19, 2024
 */
public class ComparisonHeap {
    
    private static class WorstFirstComparator implements Comparator<Comparison> {
        @Override
        public int compare(Comparison c1, Comparison c2) {
            // Invert the comparison result to get worst-first ordering
            return -c1.compareTo(c2);
        }
    }
    
    /** Creates a heap with a fixed maximum size; worst element sits at the top.
     * @param maxSize Maximum number of comparisons to retain */
    public ComparisonHeap(int maxSize) {
        this.maxSize = maxSize;
        // Create a heap ordered to keep the worst comparison at the top
        this.heap = new PriorityQueue<>(maxSize, new WorstFirstComparator());
    }
    
    /**
     * Offers a comparison; added if heap has room or if it is better than the current worst.
     * Copies the offered Comparison before storing.
     * @param comp Comparison candidate
     * @return true if added
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
     * Returns a sorted list of the top comparisons in best-first order.
     * The heap is converted to an ArrayList and sorted using natural ordering.
     * @return ArrayList of comparisons sorted from best to worst
     */
    public ArrayList<Comparison> toList() {
    	ArrayList<Comparison> result = new ArrayList<Comparison>(heap);
        // Sort using natural ordering (compareTo) which puts best first
        Collections.sort(result);
        return result;
    }
    
    /** Returns the number of comparisons currently stored.
     * @return Current heap size */
    public int size() {
        return heap.size();
    }
    
    /** Clears all stored comparisons from the heap. */
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
    
    /** Returns the worst comparison without removing it (root of the heap).
     * @return Worst comparison or null if empty */
    public Comparison worst() {
    	return heap.peek();
    }
    
    private final PriorityQueue<Comparison> heap;
    private final int maxSize;
}