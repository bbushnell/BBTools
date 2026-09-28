package prok;

import fileIO.ByteStreamWriter;
import shared.KillSwitch;
import structures.ByteBuilder;

/** Synchronized21-column TSV sink; reuses one byte buffer per shared logger.
 * The driver writes HEADER once and owns writer startup/shutdown.
 * @author G11, Raiden
 */
final class NcrnaStageDiagLogSink implements NcrnaStageDiagSink {

	NcrnaStageDiagLogSink(ByteStreamWriter writer){
		if(writer==null){throw new IllegalArgumentException("Stage diagnostics require an output writer");}
		bsw=writer;
	}
	static final String HEADER="type\tfamily\tcontig\tstrand\tscaflen\tpass\twStart\twStop"
		+"\toutcome\tkhits\tshortlistSize\tmodelIndex\tmodelName\tbestId\tbestReId\tbestHbm"
		+"\torfStart\torfStop\ttrimSucceeded\tseedCount\tseedPositions";

	@Override
	public synchronized void seed(String family, String contig, int strand, int scaflen, int[] hits){
		assert(hits!=null) : KillSwitch.assertDie("Scavenger supplies a copied, possibly empty hit stream; null would conceal missing seed evidence");
		prefix("SEED",family,contig,strand,scaflen);
		for(int i=0; i<14; i++){bb.append("\t.");}
		bb.tab().append(hits.length).tab();
		if(hits.length==0){bb.append('.');}
		else{for(int i=0; i<hits.length; i++){if(i>0){bb.append(',');}bb.append(hits[i]);}}
		bsw.print(bb.nl());
	}

	@Override
	public synchronized void window(String family, String contig, int strand, int scaflen, int pass,
		int wStart, int wStop, String outcome, int khits, int shortlistSize,
		int modelIndex, String modelName, float bestId, float bestReId, float bestHbm,
		int orfStart, int orfStop, boolean trimSucceeded){
		assert(pass==1 || pass==2) : KillSwitch.assertDie("Stage rows must identify the actual primary or nearby-hit scavenger pass");
		prefix("WINDOW",family,contig,strand,scaflen);
		bb.tab().append(pass).tab().append(wStart).tab().append(wStop).tab();safe(outcome);
		bb.tab().append(khits).tab().append(shortlistSize).tab().append(modelIndex).tab();safe(modelName);
		//Keep float round-trip precision for inspecting decisions at a cutoff.
		bb.tab().appendSlow(bestId).tab().appendSlow(bestReId).tab().appendSlow(bestHbm).tab()
			.append(orfStart).tab().append(orfStop).tab().append(trimSucceeded).append("\t.\t.\n");
		bsw.print(bb);
	}

	private void prefix(String type,String family,String contig,int strand,int scaflen){
		assert(strand==0 || strand==1) : KillSwitch.assertDie("Diagnostic coordinates require an explicit strand to recover the genomic frame");
		bb.clear().append(type).tab();safe(family);bb.tab();safe(contig);
		bb.tab().append(strand).tab().append(scaflen);
	}
	private void safe(String text){
		if(text==null){bb.append('.');return;}
		for(int i=0; i<text.length(); i++){
			final char c=text.charAt(i);bb.append(c=='\t' || c=='\n' || c=='\r' ? ' ' : c);
		}
	}
	private final ByteStreamWriter bsw;
	private final ByteBuilder bb=new ByteBuilder();
}
