package prok;

import java.util.ArrayList;
import stream.Read;

/** Optional whole-input iteration hook for the native CallGenes lifecycle.
 * One worker owns one complete input, its caller, and its mutable iteration state.
 * Implementations must execute the real caller for each policy and return the
 * designated baseline path for ordinary output. Native caller work counters then
 * describe all policy executions, not just that returned baseline.
 * @author Brian Bushnell, Raiden
 */
public interface CallGenesIteration{

	/** Called on the processing thread after native resources and its PGM load. */
	Worker create(String input, GeneCaller caller, GeneModel model);

	interface Worker{
		/** Borrows the read; every successful call must restore its base orientation. */
		ArrayList<Orf> call(Read read);
		/** Called once after the complete input, including a zero-read input. */
		void finish();
	}
}
