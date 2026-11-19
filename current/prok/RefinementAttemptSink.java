package prok;

/** Opt-in test instrumentation for rRNA consensus-alignment attempts.
 * Identity is -1f when no alignment was measured, or finite in [0,1]. */
public interface RefinementAttemptSink {
	void onAttempt(int consensusIndex, int totalConsensusCount, int preAttemptStart,
			int preAttemptStop, float preAttemptOrfScore, boolean accepted,
			int postAttemptStart, int postAttemptStop, int strand,
			String scaffold, String consensusLabel, int consensusLength, float identity,
			String reason, String candidateId);
}
