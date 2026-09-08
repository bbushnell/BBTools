package prok;

/** Opt-in test instrumentation for rRNA consensus-alignment attempts. */
interface RefinementAttemptSink {
	void onAttempt(int consensusIndex, int totalConsensusCount, int preAttemptStart,
			int preAttemptStop, float preAttemptOrfScore, boolean accepted,
			int postAttemptStart, int postAttemptStop, int strand,
			String scaffold, String consensusLabel, int consensusLength, float identity,
			String reason, String candidateId);
}
