package assemble;

/**
 * Enables existing graph and overlap diagnostics for the experimental rescue driver.
 * No branch, fusion, or graph decision is changed. Legacy name() diagnostics can
 * affect cached header statistics, so compare sequence content rather than headers.
 * @author Fischl
 */
public final class MultiKBranchRescueTrace {

	/** Restores the process-wide diagnostic flag even when argument validation fails. */
	public static void main(final String[] args){
		assert(args!=null) : "MultiKBranchRescue requires a non-null command-line array.";
		final boolean previous=BubblePopper.verbose;
		BubblePopper.verbose=true;
		try{
			MultiKBranchRescue.main(args);
		}finally{
			BubblePopper.verbose=previous;
		}
	}
}
