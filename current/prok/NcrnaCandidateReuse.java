package prok;

import java.util.ArrayList;

import shared.KillSwitch;

/** Contig-local pre-path candidates for a selected set of caller policies.
 * This stores candidates, never an accepted path. Native path selection mutates
 * predecessor/score fields and glycine reporting mutates endpoints; each reuse
 * restores discovery coordinates and a fresh path state. Lists are refilled
 * from retained references so merging new target candidates cannot grow the cache.
 * No candidate objects or primitive wrappers are allocated on restoration.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCandidateReuse{

	public NcrnaCandidateReuse(final int target_, final int families){
		this(new int[]{target_}, families);
	}

	/** Copies membership; changes to the caller's index array cannot stale the cache.
	 * The independent tRNA caller is never a selected generic-family target.
	 */
	public NcrnaCandidateReuse(final int[] targets_, final int families){
		if(families<1 || targets_==null || targets_.length<1){
			throw new IllegalArgumentException("Reuse requires at least one loaded target family");
		}
		targets=new boolean[families+1];
		for(final int family:targets_){
			if(family<0 || family>=families || targets[family]){
				throw new IllegalArgumentException("Reuse targets must be distinct loaded family indices: "+family);
			}
			targets[family]=true;
		}
		target=targets_.length==1 ? targets_[0] : -1;
		rna=new Batch[2*(families+1)];// Two strands per family, then the independent tRNA caller.
	}

	/** True for every varied member; the trailing independent tRNA slot is false. */
	public boolean isTarget(final int family){
		assert(family>=0 && family<targets.length)
			: KillSwitch.assertDie("Reuse membership uses loaded family indices, with independent tRNA last: "+family);
		return targets[family];
	}

	/** Null until the first policy supplies all six pre-path CDS frame lists. */
	public ArrayList<Orf>[] restoreCds(final GeneModel model){
		if(cds==null){return null;}
		assert(model==geneModel) : KillSwitch.assertDie("A contig's cached candidates belong to one unchanged GeneModel");
		for(int i=0; i<cds.length; i++){cdsLists[i]=cds[i].restore();}
		return cdsLists;
	}

	/** Capture before RNA merging, which mutates the ordinary frame lists. */
	public void captureCds(final ArrayList<Orf>[] lists, final GeneModel model){
		assert(cds==null && lists!=null && lists.length==6 && model!=null)
			: KillSwitch.assertDie("Capture exactly once, after CDS scoring and before RNA merging/path selection");
		geneModel=model;
		cdsLists=lists;
		cds=new Batch[lists.length];
		for(int i=0; i<lists.length; i++){cds[i]=new Batch(lists[i]);}
	}

	/** Slot is 2*family+strand; the final two slots belong to tRNA. */
	public ArrayList<Orf> restoreRna(final int slot){
		assert(slot>=0 && slot<rna.length) : KillSwitch.assertDie("RNA reuse must preserve the registry family and genomic strand");
		return rna[slot]==null ? null : rna[slot].restore();
	}

	/** Store forward-genomic raw bounds after the ordinary minus-strand flip. */
	public void captureRna(final int slot, final ArrayList<Orf> calls){
		assert(slot>=0 && slot<rna.length && rna[slot]==null && !isTarget(slot/2))
			: KillSwitch.assertDie("Only unchanged non-target candidates may be captured, once per family/strand");
		rna[slot]=new Batch(calls);
	}

	/** Original membership and discovery bounds survive list merging and reporting. */
	private static final class Batch{
		Batch(final ArrayList<Orf> candidates){
			list=candidates==null ? new ArrayList<Orf>(0) : candidates;
			orfs=list.toArray(new Orf[list.size()]);
			starts=new int[orfs.length];
			stops=new int[orfs.length];
			for(int i=0; i<orfs.length; i++){
				final Orf orf=orfs[i];
				assert(orf!=null && orf.prevPlus==null && orf.prevMinus==null)
					: KillSwitch.assertDie("Candidate capture must precede GeneCaller.findPath; accepted paths are not reusable inputs");
				starts[i]=orf.start;
				stops[i]=orf.stop;
			}
		}

		ArrayList<Orf> restore(){
			assert(starts.length==orfs.length && stops.length==orfs.length)
				: KillSwitch.assertDie("Every retained candidate needs both original discovery coordinates");
			list.clear();
			for(int i=0; i<orfs.length; i++){
				final Orf orf=orfs[i];
				orf.start=starts[i];
				orf.stop=stops[i];
				// findPath resets scores, but its root logic also requires null links.
				orf.prevPlus=null;
				orf.prevMinus=null;
				orf.pathScorePlus=0;
				orf.pathScoreMinus=0;
				orf.pathLengthPlus=1;
				orf.pathLengthMinus=1;
				list.add(orf);
			}
			return list;
		}

		final ArrayList<Orf> list;
		final Orf[] orfs;
		final int[] starts, stops;
	}

	/** Legacy single-target index, or -1 for a group. Use isTarget for membership. */
	public final int target;
	private final boolean[] targets;
	private final Batch[] rna;
	private Batch[] cds;
	private ArrayList<Orf>[] cdsLists;
	private GeneModel geneModel;
}
