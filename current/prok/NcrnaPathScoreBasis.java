package prok;

import java.util.ArrayList;
import java.util.Objects;

import map.ObjectIntMap;
import structures.FloatList;
import structures.IntList;

/** Worker-local score provenance keyed to native RNA candidate identity.
 * Primitive columns reuse their capacity between contigs; rescue overwrites the
 * same row. No CDS fields, genomic windows or extra formula evaluations occur.
 * @author Brian, Raiden
 */
public final class NcrnaPathScoreBasis{

	/** One binding per active scavenger, created outside the candidate hot path. */
	public NcrnaScavenger.PathScoreObserver bind(final NcrnaScavenger scorer){
		require(scorer!=null, "A score-basis binding requires the actual native scavenger");
		return (candidate, model, length, identity)->record(scorer, candidate, model, length, identity);
	}

	/** Drop native candidate references only after the preceding contig is consumed. */
	public void clear(){
		rows.clear();
		models.clear();
		lengths.clear();
		identities.clear();
		scores.clear();
		additives.clear();
		multipliers.clear();
		families.clear();
		names.clear();
	}

	private void record(final NcrnaScavenger scorer, final Orf candidate, final int model,
			final int length, final float identity){
		require(candidate!=null && model>=0 && model<scorer.modelCount(),
			"The assignment callback must retain the exact candidate and native model index");
		final String family=scorer.pathScoreFamily(model);
		if(!NcrnaFamilyRegistry.isCobalamin(family)){return;}
		final String name=scorer.pathScoreModel(model);
		final float a=scorer.pathScoreA(model);
		final float b=scorer.pathScoreB(model);
		require(candidate.type==ProkObject.RNA && name!=null && !name.isEmpty() && length>0
			&& Float.isFinite(identity) && Float.isFinite(a) && Float.isFinite(b) && Float.isFinite(candidate.orfScore),
			"A target score basis needs its actual model label, owner, finite coefficients and original aligned span/identity");
		int row=rows.get(candidate);
		if(row<0){
			row=models.size;
			rows.put(candidate, row);
			families.add(family);
			names.add(name);
		}else{
			families.set(row, family);
			names.set(row, name);
		}
		models.set(row, model);
		lengths.set(row, length);
		identities.set(row, identity);
		scores.set(row, candidate.orfScore);
		additives.set(row, a);
		multipliers.set(row, b);
	}

	/** Resolve the object while the native path seam owns it; publish scalars only. */
	void emit(final Orf candidate, final long id, final GeneCallerPathTrace.Observer observer){
		assert(candidate!=null && observer!=null) : "Path-basis emission runs inside an installed native observation seam";
		final int row=rows.get(candidate);
		if(row<0){
			require(!NcrnaFamilyRegistry.isCobalamin(candidate.ncrnaFamily),
				"Every returned target candidate must retain its own native assignment basis; coordinates cannot substitute for identity");
			return;
		}
		require(Objects.equals(candidate.ncrnaFamily, families.get(row))
			&& (candidate.trnaModel==null || candidate.trnaModel.equals(names.get(row)))
			&& Float.floatToRawIntBits(candidate.orfScore)==Float.floatToRawIntBits(scores.get(row)),
			"Final target owner/model/score must match the latest assignment, including a rescue update");
		observer.scoreBasis(id, families.get(row), names.get(row), models.get(row), lengths.get(row),
			identities.get(row), additives.get(row), multipliers.get(row), scores.get(row));
	}

	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private final ObjectIntMap<Orf> rows=new ObjectIntMap<Orf>(Orf.class);
	private final IntList models=new IntList(), lengths=new IntList();
	private final FloatList identities=new FloatList(), scores=new FloatList(),
		additives=new FloatList(), multipliers=new FloatList();
	private final ArrayList<String> families=new ArrayList<String>(), names=new ArrayList<String>();
}
