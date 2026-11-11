package assemble;

/**
 * Runs the ordinary Tadpole pipeline with existing BubblePopper diagnostics.
 * Unlike verbose=t, this does not enable per-kmer traversal logging.
 * This development launcher changes logging only, not assembly decisions.
 * Rebuilt contigs invalidate the scalar caches populated by header logging.
 * The logging flag is process-wide; run this launcher in its own JVM.
 * @author Brian Bushnell, Fischl
 */
public class TadpoleMergeTrace{

	/** Passes the unchanged command line to Tadpole while logging graph merges. */
	public static void main(final String[] args){
		assert(args!=null) : "Tadpole dispatch requires a non-null argument array.";
		final boolean previous=BubblePopper.verbose;
		BubblePopper.verbose=true;
		try{
			Tadpole.main(args);
		}finally{
			BubblePopper.verbose=previous;
		}
	}
}
