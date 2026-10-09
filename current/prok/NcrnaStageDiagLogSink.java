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
	@Override
	public synchronized void model(String family,String contig,int strand,int scaflen,int pass,
		int start,int stop,int modelIndex,String modelName,float identity){
		assert((pass==1 || pass==2) && modelIndex>=0 && Float.isFinite(identity)) : KillSwitch.assertDie(
			"Model evidence must identify an actual initial alignment and finite measured identity");
		prefix("MODEL",family,contig,strand,scaflen);
		bb.tab().append(pass).tab().append(start).tab().append(stop).append("\tALIGNMENT\t.\t.\t").append(modelIndex).tab();safe(modelName);
		// MODEL uses bestId for this model's actual identity; unused window fields
		// stay unavailable rather than inventing seed counts or acceptance decisions.
		bb.tab().appendSlow(identity).append("\t.\t.\t.\t.\t.\t.\t.\n");bsw.print(bb);
	}
	@Override
	public synchronized void path(String family,String contig,int strand,int scaflen,
		String modelName,int start,int stop,boolean retained){
		assert(start>=0 && stop>=start && stop<scaflen) : KillSwitch.assertDie(
			"DP diagnostic spans must fit the oriented contig to join verifier and emitted-call evidence");
		prefix("PATH",family,contig,strand,scaflen);
		bb.append("\t.\t.\t.\t").append(retained ? "DP_RETAINED" : "DP_INPUT");
		bb.append("\t.\t.\t.\t");safe(modelName);
		bb.append("\t.\t.\t.\t").append(start).tab().append(stop).append("\t.\t.\t.\n");
		bsw.print(bb);
	}
	private void safe(String text){
		if(text==null){bb.append('.');return;}
		for(int i=0; i<text.length(); i++){
			final char c=text.charAt(i);bb.append(c=='\t' || c=='\n' || c=='\r' ? ' ' : c);
		}
	}
	/** Opt-in supplemental TSV rows in captured caller stderr, leaving the stage
	 * file schema and all existing graders unchanged. Columns after the tag:
	 * family, contig, strand, pass, windowStart, windowStop, model, queryStart,
	 * queryStop, queryLength, refStart, refStop, refLength, leftClip, rightClip,
	 * identity, outcome. Local coordinates are inclusive and strand-oriented. */
	@Override
	public synchronized void modelClip(String family,String contig,int strand,int pass,int wStart,int wStop,
		int model,idaligner.EndClippedAligner.Result r,int minLen){
		assert(minLen>0) : KillSwitch.assertDie("Local alignment diagnostics need the real caller minimum to distinguish short winners");
		bb.clear().append("NCRNA_MODEL_CLIP\t");safe(family);bb.tab();safe(contig);
		bb.tab().append(strand).tab().append(pass).tab().append(wStart).tab().append(wStop).tab().append(model);
		if(r==null){bb.append("\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t0\tNO_POSITIVE_ALIGNMENT");}
		else{
			bb.tab().append(r.qStart).tab().append(r.qStop).tab().append(r.qLen).tab().append(r.rStart).tab().append(r.rStop).tab().append(r.rLen)
				.tab().append(r.qStart).tab().append(r.qLen-r.qStop-1).tab().appendSlow(r.identity).tab()
				.append(r.rStop-r.rStart+1<minLen ? "BEST_LOCAL_BELOW_MINLEN" : "SPAN_ELIGIBLE");
		}
		System.err.println(bb.toString());
	}
	/** Supplemental stderr columns: family,contig,strand,pass,wStart,wStop,model,
	 * worker,qLen,rLen,nativeScore,rStart,rStop,padLeft,padRight,X,Y,gapIdentity,
	 * flatIdentity,outcome,primaryNs,primaryCpuNs,allocationNs,workingMatrixBytes,
	 * maxQuery,maxRef,retainedMatrixBytes,quantumIdentity,quantumStart,quantumStop,
	 * quantumNs,quantumCpuNs. CPU -1 means unavailable. Matrix values are primitive
	 * payload, not heap peaks; max retained bytes per worker can be summed. */
	@Override
	public synchronized void modelPacBio(String family,String contig,int strand,int pass,int wStart,int wStop,
		int model,Euk18sPacBioAligner.Result r,float cutoff,int maxLen){
		assert(r!=null && r.comparedQuantum) : KillSwitch.assertDie("PacBio diagnostic rows require the paired actual Quantum attempt");
		bb.clear().append("NCRNA_MODEL_PACBIO\t");safe(family);bb.tab();safe(contig);
		bb.tab().append(strand).tab().append(pass).tab().append(wStart).tab().append(wStop).tab().append(model)
			.tab().append(r.worker).tab().append(r.qLen).tab().append(r.rLen).tab().append(r.score)
			.tab().append(r.rStart).tab().append(r.rStop).tab().append(r.padLeft).tab().append(r.padRight)
			.tab().append(r.leftOverhang).tab().append(r.rightOverhang).tab().appendSlow(r.identity)
			.tab().appendSlow(r.flatIdentity).tab().append(r.outcome(cutoff,maxLen))
			.tab().append(r.primaryNanos).tab().append(r.primaryCpuNanos).tab().append(r.allocationNanos)
			.tab().append(r.workingMatrixBytes).tab().append(r.maxQuery).tab().append(r.maxRef).tab().append(r.retainedMatrixBytes)
			.tab().appendSlow(r.quantumIdentity).tab().append(r.quantumStart).tab().append(r.quantumStop)
			.tab().append(r.quantumNanos).tab().append(r.quantumCpuNanos).nl();
		System.err.write(bb.array,0,bb.length);
	}
	private final ByteStreamWriter bsw;
	private final ByteBuilder bb=new ByteBuilder();
}
