package prok;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Synchronized30-column generated-vote TSV. The driver owns header and writer lifecycle.
 * @author Keqing
 */
final class NcrnaVoteDiagLogSink implements NcrnaVoteDiagSink {
	NcrnaVoteDiagLogSink(ByteStreamWriter writer){
		if(writer==null){throw new IllegalArgumentException("Vote diagnostics require an output writer");}
		bsw=writer;
	}
	static final String HEADER="family\tcontig\tstrand\tscaflen\tpass\tgroup_id\tpolicy\tdisposition\treason"
		+"\tinput_start\tinput_end\traw_start\traw_end\trounded_start\trounded_end\tgenerated_start\tgenerated_end"
		+"\tslack\tfallback_pad\ttrained_voters\tmember_centers\tmember_keys\tmember_trained\tmember_parent_ids"
		+"\tleft_weight_sum\tright_weight_sum\tleft_pred_range\tright_pred_range\tleft_top10_sd\tright_top10_sd";

	@Override
	public synchronized void group(NcrnaVoteDiagnostics.Group g){
		assert(g!=null) : "A vote row must come from an executed proposal callback";
		bb.clear();safe(g.family);bb.tab();safe(g.contig);
		bb.tab().append(g.strand).tab().append(g.length).tab().append(g.pass).tab().append(g.id).tab();safe(g.policy);
		bb.tab();safe(g.disposition);bb.tab();safe(g.reason);
		integer(g.inputStart);integer(g.inputEnd);number(g.rawStart);number(g.rawEnd);
		integer(g.roundedStart);integer(g.roundedEnd);integer(g.start);integer(g.stop);
		bb.tab().append(g.slack).tab().append(g.pad).tab().append(g.voters);
		list(g.centers, false);bb.tab();
		if(g.keys.length==0){bb.append('.');}
		for(int i=0; i<g.keys.length; i++){
			if(i>0){bb.append(',');}
			if(g.keys[i]==NcrnaVoteDiagnostics.UNKNOWN){bb.append('.');}else{bb.append(Long.toHexString(g.keys[i]));}
		}
		list(g.trained, true);list(g.parents, true);
		number(g.leftWeight);number(g.rightWeight);number(g.leftRange);number(g.rightRange);number(g.leftSd);number(g.rightSd);
		bsw.print(bb.nl());
	}
	private void integer(long value){bb.tab();if(value==NcrnaVoteDiagnostics.UNKNOWN){bb.append('.');}else{bb.append(value);}}
	private void number(double value){bb.tab();if(Double.isNaN(value)){bb.append('.');}else{bb.appendSlow(value);}}
	private void list(int[] values, boolean missingNegative){
		bb.tab();if(values.length==0){bb.append('.');}
		for(int i=0; i<values.length; i++){
			if(i>0){bb.append(',');}
			if(missingNegative && values[i]<0){bb.append('.');}else{bb.append(values[i]);}
		}
	}
	private void safe(String text){
		if(text==null){bb.append('.');return;}
		for(int i=0; i<text.length(); i++){
			final char c=text.charAt(i);bb.append(c=='\t' || c=='\r' || c=='\n' ? ' ' : c);
		}
	}
	private final ByteStreamWriter bsw;
	private final ByteBuilder bb=new ByteBuilder();
}
