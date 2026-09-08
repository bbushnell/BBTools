package prok;

import fileIO.ByteStreamWriter;

/** Opt-in TSV observer for rRNA consensus-alignment attempts. */
final class AttemptLogSink implements RefinementAttemptSink {
	AttemptLogSink(ByteStreamWriter bsw_){bsw=bsw_;}
	final ByteStreamWriter bsw;
	@Override
	public synchronized void onAttempt(int index, int total, int preStart, int preStop, float preScore, boolean accepted,
			int postStart, int postStop, int strand, String scaffold, String label, int consensusLength,
			float identity, String reason, String candidateId){
		final String safeScaffold=scaffold.replace('\t', ' ').replace('\n', ' ').replace('\r', ' ');
		final String safeLabel=label.replace('\t', ' ').replace('\n', ' ').replace('\r', ' ');
		bsw.println(candidateId+"\t"+index+"\t"+total+"\t"+preStart+"\t"+preStop+"\t"+postStart+"\t"+postStop
				+"\t"+strand+"\t"+safeScaffold+"\t"+safeLabel+"\t"+consensusLength+"\t"+preScore+"\t"+identity+"\t"+reason+"\t"+accepted);
	}
}
