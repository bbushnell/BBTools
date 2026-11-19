package prok;

import fileIO.ByteStreamWriter;

/** Opt-in TSV observer for rRNA consensus-alignment attempts. */
public final class AttemptLogSink implements RefinementAttemptSink {
	public AttemptLogSink(ByteStreamWriter bsw_){bsw=bsw_;}
	final ByteStreamWriter bsw;
	public static final String SCHEMA="#identity\tEMPTY_CELL_MEANS_UNMEASURED;OTHERWISE_FINITE_0_TO_1";
	public static final String COLUMNS="candidateId\tconsensusIndex\ttotalConsensusCount\tpreStart\tpreStop\tpostStart"
		+"\tpostStop\tstrand\tscaffold\tconsensusLabel\tconsensusLength\tpreOrfScore"
		+"\tidentity\treason\taccepted";
	public static final String HEADER=SCHEMA+"\n"+COLUMNS;
	@Override
	public synchronized void onAttempt(int index, int total, int preStart, int preStop, float preScore, boolean accepted,
			int postStart, int postStop, int strand, String scaffold, String label, int consensusLength,
			float identity, String reason, String candidateId){
		if(identity!=-1f && (!Float.isFinite(identity) || identity<0 || identity>1)){
			throw new IllegalArgumentException("Attempt identity must be -1 or a finite fraction: "+identity);
		}
		if(accepted && identity<0){
			throw new IllegalArgumentException("Accepted attempts require measured identity");
		}
		final String safeScaffold=scaffold.replace('\t', ' ').replace('\n', ' ').replace('\r', ' ');
		final String safeLabel=label.replace('\t', ' ').replace('\n', ' ').replace('\r', ' ');
		bsw.println(candidateId+"\t"+index+"\t"+total+"\t"+preStart+"\t"+preStop+"\t"+postStart+"\t"+postStop
				+"\t"+strand+"\t"+safeScaffold+"\t"+safeLabel+"\t"+consensusLength+"\t"+preScore+"\t"
				+(identity<0 ? "" : Float.toString(identity))+"\t"+reason+"\t"+accepted);
	}
}
