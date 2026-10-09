package assemble;

/**
 * Runs the ordinary Tadpole pipeline with existing BubblePopper diagnostics.
 * Unlike verbose=t, this does not enable per-kmer traversal logging.
 * This development launcher changes logging only, not assembly decisions.
 * Existing name() logging can change cached GC/HH/CAGA header statistics.
 * @author Brian Bushnell, Fischl
 */
public class TadpoleMergeTrace {

	/** Passes the unchanged command line to Tadpole while logging graph merges. */
	public static void main(final String[] args){
		assert(args!=null) : "Tadpole dispatch requires a non-null argument array.";
		final boolean previous=BubblePopper.verbose;
		//TODO: BubblePopper.merge changes bases without invalidating Contig's scalar
		//caches; verbose name() calls expose stale GC/HH/CAGA in the final headers.
		BubblePopper.verbose=true;
		try{
			Tadpole.main(args);
		}finally{
			BubblePopper.verbose=previous;
		}
	}
}
