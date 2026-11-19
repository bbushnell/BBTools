package prok;

import shared.KillSwitch;
import structures.ByteBuilder;

/** RF00174 experimental call carrying proposal provenance through final selection.
 * Only the enabled experiment creates this subtype; ordinary ORFs gain no fields.
 * Metadata changes neither acceptance, coordinates nor path score.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlOrf extends Orf{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaDdlOrf(final String name, final int start, final int stop,
			final int strand, final byte[] bases, final int outputType){
		super(name, start, stop, strand, 0, bases, false, outputType);
		assert(start>=0 && stop>=start && stop<bases.length)
			: KillSwitch.assertDie("Hybrid provenance belongs to a physical RNA candidate before strand flipping");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Attach only after the final primary/rescue winning model has been selected. */
	public void setProposal(final int count, final boolean fromSketch){
		assert(!attached && count>=-1 && count<=65535 && (count>=0 || fromSketch))
			: KillSwitch.assertDie("Attach one accepted-model count after rescue; -1 means unmeasured native support for sketch-only ordering");
		nativeCount=count;
		sketchFirst=fromSketch;
		attached=true;
		final NcrnaAlignmentWork work=NcrnaAlignmentWork.current();
		if(work!=null){work.proposal(count, fromSketch);}
	}

	/** Selected-format evidence retains the measured winner, never reconstructing identity from path score. */
	public void setSelectedProposal(final int rank, final float score, final float identity){
		setSelectedProposal(rank, score, identity, false);
	}

	/** Only explicit GC scoring may attach a negative adjusted winner score. */
	public void setGcSelectedProposal(final int rank, final float score, final float identity){
		setSelectedProposal(rank, score, identity, true);
	}

	private void setSelectedProposal(final int rank, final float score, final float identity, final boolean signed){
		assert(attached && !selected && rank>0 && Float.isFinite(score) && (signed || score>=0)
			&& Float.isFinite(identity) && identity>=0 && identity<=1)
			: KillSwitch.assertDie("Selected proposal evidence needs finite rank/identity/score; only the explicit GC-adjusted route permits negative scores");
		proposalRank=rank;
		proposalScore=score;
		acceptedIdentity=identity;
		selected=true;
		gcAdjusted=signed;
	}

	public float acceptedIdentity(){return acceptedIdentity;}
	public int proposalRank(){return proposalRank;}
	public float proposalScore(){return proposalScore;}

	@Override
	public ByteBuilder appendGff(final ByteBuilder b){
		assert(attached && "cobalamin_rf00174".equals(ncrnaFamily))
			: KillSwitch.assertDie("Only accepted experimental RF00174 calls have proposal provenance at final GFF output");
		super.appendGff(b);
		b.append(";ddl_stream=").append(nativeCount<0 ? "sketch_only" : sketchFirst ? "sketch" : "native")
			.append(";ddl_native_shared=");
		if(nativeCount<0){b.append("NA");}
		else{b.append(nativeCount);}
		if(selected){
			b.append(";ddl_rank=").append(proposalRank).append(";ddl_score=").append(proposalScore, 3)
				.append(";accepted_identity=").append(acceptedIdentity, 6);
			if(gcAdjusted){b.append(";ddl_gc_adjusted=true");}
		}
		return b;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private int nativeCount=-1;
	private boolean sketchFirst, attached;
	private boolean selected, gcAdjusted;
	private int proposalRank;
	private float proposalScore=Float.NaN, acceptedIdentity=Float.NaN;
}
