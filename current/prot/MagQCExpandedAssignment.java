package prot;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import map.IntHashMap;
import structures.ByteBuilder;

/** Experimental full8537 competition projected onto an unchanged4537-input model.
 * New-only and masked winners contribute no tracked copy; there is no second-choice fallback.
 * @author Brian Bushnell, Yoimiya */
public final class MagQCExpandedAssignment{

	/** All three explicit options are required together; absence preserves the legacy path. */
	public static boolean configured(HashMap<String,String> options){
		if(options==null){throw new IllegalArgumentException("Assignment options are required");}
		int found=0;for(String key:KEYS){if(options.containsKey(key)){found++;}}
		if(found!=0 && found!=KEYS.length){throw new IllegalArgumentException("Expanded assignment requires expandedmodels, expandedmodelsha80 and expandedidentitysha80 together");}
		if(found>0){
			DigestSuffix.requireSuffix(options.get(KEYS[1]),KEYS[1]);
			DigestSuffix.requireSuffix(options.get(KEYS[2]),KEYS[2]);
		}
		return found>0;
	}

	static MagQCExpandedAssignment load(HashMap<String,String> options) throws IOException{
		if(!configured(options)){throw new IllegalArgumentException("Expanded assignment was not explicitly configured");}
		final Schema7FamilyAcceptanceLoader.FeatureLayout layout=Schema7FamilyAcceptanceLoader.loadFeatureLayout(
			need(options,"profile"),need(options,"profilesha80"),need(options,"roster"),need(options,"rolemanifest"),
			need(options,"familylist"),need(options,"familylistsha80"));
		if(layout.ids.length!=4537){throw new IllegalArgumentException("This experiment preserves exactly4537 old feature ranks");}
		final ProteinFamilyAssigner library=ProteinFamilyAssigner.load(Paths.get(need(options,"expandedmodels")));
		if(library.familyCount()!=8537 || !library.modelPin().equals(need(options,"expandedmodelsha80")) ||
				!library.identityPin().equals(need(options,"expandedidentitysha80"))){
			throw new IllegalArgumentException("Expanded reference differs from the declared8537-family package");
		}
		final int[] ids=new int[library.familyCount()];final String[] reps=new String[ids.length];
		for(int i=0;i<ids.length;i++){ids[i]=library.familyId(i);reps[i]=library.representative(i);}
		return new MagQCExpandedAssignment(library,new Projection(layout.ids,layout.reps,layout.enabled,ids,reps),layout);
	}
	private MagQCExpandedAssignment(ProteinFamilyAssigner library_,Projection projection_,Schema7FamilyAcceptanceLoader.FeatureLayout layout_){
		assert(library_!=null && projection_!=null && layout_!=null) : "Expanded workers require authenticated reference and old feature layout";
		library=library_;projection=projection_;layout=layout_;
	}
	Worker newWorker(){return new Worker();}
	int familyCount(){return projection.familyCount();}
	void appendBindings(ByteBuilder out){
		assert(out!=null) : "Expanded assignment provenance belongs in the same completed panel report";
		out.append("#assignment_backend\texpanded_reference_projected\n#assignment_rule\t").append(ProteinFamilyAssigner.RULE)
			.append("\n#expandedmodelsha80\t").append(library.modelPin()).append("\n#expandedidentitysha80\t").append(library.identityPin())
			.append("\n#feature_roster_sha80\t").append(layout.rosterPin).append("\n#feature_roles_sha80\t").append(layout.rolesPin)
			.append("\n#assignment_projection\tpermanent_family_id_and_representative;new_or_masked_winner_has_no_feature;no_fallback\n");
	}
	final class Worker{
		private Worker(){aligner=library.newWorker();}
		int familyCount(){return projection.familyCount();}
		int assign(ProteinSequence protein){
			assert(protein!=null) : "Projection scores the exact caller translation, never a replacement query";
			final ProteinFamilyAssigner.Hit hit=aligner.assign(protein);
			return projection.project(hit.rank);
		}
		private final ProteinFamilyAssigner.Worker aligner;
	}

	/** Immutable mapping. Public for independent fixture oracles, not a synthetic assignment backend. */
	public static final class Projection{
		public Projection(int[] oldIds,String[] oldReps,boolean[] enabled,int[] newIds,String[] newReps){
			if(oldIds==null || oldReps==null || enabled==null || newIds==null || newReps==null || oldIds.length==0 ||
					oldIds.length!=oldReps.length || oldIds.length!=enabled.length || newIds.length!=newReps.length || newIds.length<oldIds.length){
				throw new IllegalArgumentException("Projection requires complete old/new identity arrays and old feature masks");
			}
			families=oldIds.length;ranks=new int[newIds.length];Arrays.fill(ranks,-1);
			final IntHashMap old=new IntHashMap(),seen=new IntHashMap();
			final HashSet<String> oldNames=new HashSet<String>(),newNames=new HashSet<String>();
			for(int i=0;i<oldIds.length;i++){
				if(oldIds[i]<0 || old.containsKey(oldIds[i]) || oldReps[i]==null || oldReps[i].isEmpty() || !oldNames.add(oldReps[i])){
					throw new IllegalArgumentException("Duplicate or invalid old family identity at rank "+i);
				}
				old.put(oldIds[i],i);
			}
			int covered=0;
			for(int i=0;i<newIds.length;i++){
				if(newIds[i]<0 || seen.containsKey(newIds[i]) || newReps[i]==null || newReps[i].isEmpty() || !newNames.add(newReps[i])){
					throw new IllegalArgumentException("Duplicate or invalid expanded family identity at rank "+i);
				}
				seen.put(newIds[i],i);final int rank=old.get(newIds[i]);
				if(rank>=0){
					if(!oldReps[rank].equals(newReps[i])){throw new IllegalArgumentException("Representative changed for permanent family "+newIds[i]);}
					covered++;if(enabled[rank]){ranks[i]=rank;}
				}else if(oldNames.contains(newReps[i])){throw new IllegalArgumentException("Old representative was rebound to a new permanent family ID");}
			}
			if(covered!=families){throw new IllegalArgumentException("Expanded reference is missing old family identities: "+covered+"/"+families);}
		}
		public int familyCount(){return families;}
		public int project(int winner){
			if(winner< -1 || winner>=ranks.length){throw new IllegalArgumentException("Expanded winner rank out of range: "+winner);}
			return winner<0 ? -1 : ranks[winner];
		}
		private final int families;private final int[] ranks;
	}
	private static String need(HashMap<String,String> options,String key){
		final String value=options.get(key);if(value==null || value.isEmpty()){throw new IllegalArgumentException("Missing "+key); }return value;
	}
	private static final String[] KEYS={"expandedmodels","expandedmodelsha80","expandedidentitysha80"};
	private final ProteinFamilyAssigner library;
	private final Projection projection;
	private final Schema7FamilyAcceptanceLoader.FeatureLayout layout;
}
