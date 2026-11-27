package prok;

/** Optional worker-owned decision after native endpoint processing and before
 * claims/path selection. Model ranking and walk stopping keep their old gates.
 * @author Brian Bushnell, Raiden
 */
public interface NcrnaCandidateAcceptance{

	/** Candidate admission only; zero keeps each native admission gate unchanged.
	 * Never substitute this for the walk's pass gate. */
	float floor();

	/** Receives sense-oriented coordinates/bases and the chosen global model.
	 * Implementations must preserve endpoints, model identity and path score. */
	boolean accept(Orf candidate, byte[] bases, int model, float identity);
}
