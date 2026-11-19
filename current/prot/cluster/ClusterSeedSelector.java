package prot.cluster;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Comparator;
import java.util.PriorityQueue;

import fileIO.ByteStreamWriter;
import prot.ProteinSequence;
import structures.ByteBuilder;

/** Global size selection across validated greedy codecs, without merging or rebuilding.
 * Codecs are visited sequentially; only topN candidates are retained across files.
 * Global IDs are original concatenated record ordinals; ties retain the earlier ordinal.
 * @author Brian Bushnell, Keqing
 */
public final class ClusterSeedSelector{

	/** Selects one original codec before its first ranked assignment. */
	public static RankedSelection selectRanked(final String codec,final int familyBudget,final double initialFraction,
			final Path members,final String membersPin,final ProteinSequence[] population,final int layer){
		return selectRanked(new String[]{codec},familyBudget,initialFraction,members,membersPin,population,layer);
	}

	/** Selects ranked seed models and binds them to one complete input population.
	 * Pinned membership hydration is restricted to one original codec; multiple
	 * serialized codecs must already carry validated immutable member descriptors.
	 * Selection retains explicit source entries, never reparses renamed models or
	 * constructs a source ordinal map from consensus state.
	 */
	public static RankedSelection selectRanked(final String[] codecs,final int familyBudget,final double initialFraction,
			final Path members,final String membersPin,final ProteinSequence[] population,final int layer){
		require(population!=null,"Ranked seed selection requires the complete original population");
		final boolean hydrated=members!=null || membersPin!=null;
		require(!hydrated || members!=null && membersPin!=null,
			"Seed hydration requires both seedmembers and seedmemberssha80");
		final Selection selected=select(codecs,familyBudget,initialFraction);
		if(!hydrated){return selectSerialized(selected,population);}
		require(selected.inputCodecs==1,"Pinned seedmembers cannot unambiguously hydrate multiple source codecs");
		final int[] ids=new int[selected.selected.length],counts=new int[ids.length];
		for(int i=0;i<ids.length;i++){
			final Seed seed=selected.selected[i];
			assert(seed.shard==0) : "One-codec hydration must not mix shard-local source-family namespaces";
			ids[i]=seed.entry.familyId;counts[i]=seed.entry.members;
		}
		final SeedMemberRankings rankings=SeedMemberRankings.read(members,membersPin,population,ids,counts,layer);
		final ClusterSketchCodec.Entry[] entries=new ClusterSketchCodec.Entry[ids.length];
		for(int i=0;i<entries.length;i++){
			entries[i]=ClusterSketchCodec.withMemberRanking(selected.models.families[i],rankings.family(ids[i]));
		}
		return new RankedSelection(ClusterSketchCodec.resident(selected.models.header,entries),rankings.inputs,
			selected.nextFamilyId);
	}

	/** Serialized entries already bind member descriptors; only the input cache is new. */
	private static RankedSelection selectSerialized(final Selection selected,final ProteinSequence[] population){
		for(ClusterSketchCodec.Entry entry:selected.models.families){
			require(entry.memberRanking!=null,"Ranked seed selection requires serialized member descriptors: "+entry.name);
		}
		return new RankedSelection(selected.models,new MemberRankInputs(population),selected.nextFamilyId);
	}

	/** Resident selected models and the exact original ordinal cache used to hydrate them. */
	public static final class RankedSelection{
		private RankedSelection(final ClusterSketchCodec.Loaded models_,final MemberRankInputs inputs_,final int nextId){
			assert(models_!=null && inputs_!=null && nextId>=models_.families.length) :
				"Ranked selection retains its checked population, selected roster and unreused family-ID cursor";
			models=models_;inputs=inputs_;nextFamilyId=nextId;
		}
		public final ClusterSketchCodec.Loaded models;
		public final MemberRankInputs inputs;
		public final int nextFamilyId;
	}

	/** Writes selected raw models, their source mapping and the unreused ID cursor. */
	public static void write(final String[] paths,final int familyBudget,final double initialFraction,final String output,
			final String roster,final String report,final boolean overwrite){
		require(output!=null && roster!=null && report!=null,"Selection needs out, roster and report");
		DdlClustStageB.PathRegistry.prepare().inputs(paths).outputs(output,roster,report);
		final Selection selection=select(paths,familyBudget,initialFraction);
		writeSelection(selection,output,roster,report,overwrite);
	}

	/** Returns the same globally ordered selection directly to a resident cake run. */
	static Selection select(final String[] paths,final int familyBudget,final double initialFraction){
		final int topN=requestedCount(familyBudget,initialFraction);
		require(paths!=null && paths.length>0,"Selection needs input codecs");
		final PriorityQueue<Seed> kept=new PriorityQueue<Seed>(Math.min(topN,1024),WORST_FIRST);
		ClusterSketchCodec.Header geometry=null;
		int sourceFamilies=0;
		long sourceMembers=0;
		for(int shard=0; shard<paths.length; shard++){
			final ClusterSketchCodec.Loaded loaded=ClusterSketchCodec.read(paths[shard]);
			if(geometry==null){geometry=loaded.header;}
			else{sameGeometry(geometry,loaded.header);}
			for(ClusterSketchCodec.Entry entry:loaded.families){
				require(entry.members>0,"Greedy selection requires positive source counts: "+entry.name);
				require(sourceFamilies<Integer.MAX_VALUE-1,"Global source IDs must leave room for the next-family cursor");
				final int ordinal=sourceFamilies++;
				sourceMembers=Math.addExact(sourceMembers,(long)entry.members);
				// Later equal-size entries cannot beat the earlier retained ordinal.
				if(kept.size()==topN && entry.members<=kept.peek().entry.members){continue;}
				if(kept.size()==topN){kept.poll();}
				kept.add(new Seed(ordinal,shard,entry));
			}
		}
		require(!kept.isEmpty(),"Global seed selection needs at least one nonempty source family");
		final Seed[] selected=kept.toArray(new Seed[kept.size()]);
		Arrays.sort(selected,BEST_FIRST);
		final ClusterSketchCodec.Header header=geometry.withFamilies(selected.length);
		final ClusterSketchCodec.Entry[] entries=new ClusterSketchCodec.Entry[selected.length];
		long selectedMembers=0;
		for(int i=0;i<selected.length;i++){
			final Seed seed=selected[i];
			entries[i]=ClusterSketchCodec.selectedEntry(seed.entry,seed.ordinal,"s"+seed.shard+"_f"+seed.entry.familyId);
			selectedMembers=Math.addExact(selectedMembers,(long)seed.entry.members);
		}
		require(selectedMembers<=sourceMembers,"Size selection cannot invent source membership");
		return new Selection(ClusterSketchCodec.resident(header,entries),selected,sourceFamilies,sourceMembers,
			selectedMembers,paths.length,familyBudget,initialFraction,topN);
	}

	private static void writeSelection(final Selection selection,final String output,final String roster,
			final String report,final boolean overwrite){
		final Seed[] selected=selection.selected;
		final ClusterSketchCodec.Header header=selection.models.header;
		ByteStreamWriter models=null,mapping=null,summary=null;
		try{
			models=open(output,overwrite);mapping=open(roster,overwrite);summary=open(report,overwrite);
			final ByteBuilder bytes=new ByteBuilder();
			ClusterSketchCodec.appendHeader(header,bytes);models.print(bytes);
			mapping.println("family_id\tsource_shard\tsource_family_id\tsource_name\tsource_members");
			for(Seed seed:selected){
				final String name="s"+seed.shard+"_f"+seed.entry.familyId;
				bytes.clear();ClusterSketchCodec.appendEntry(seed.entry,seed.ordinal,name,header,bytes);models.print(bytes);
				bytes.clear().append(seed.ordinal).tab().append(seed.shard).tab().append(seed.entry.familyId)
					.tab().append(seed.entry.name).tab().append(seed.entry.members);mapping.println(bytes);
			}
			summary.println("key\tvalue");
			summary.println("schema\tglobal_seed_selection_v1");
			summary.println("input_codecs\t"+selection.inputCodecs);
			summary.println("family_budget\t"+selection.budget);
			summary.println("initial_fraction\t"+selection.fraction);
			summary.println("source_families\t"+selection.nextFamilyId);
			summary.println("requested_families\t"+selection.requested);
			summary.println("selected_families\t"+selected.length);
			summary.println("source_members\t"+selection.sourceMembers);
			summary.println("selected_source_members\t"+selection.selectedMembers);
			summary.println("dissolved_source_members\t"+(selection.sourceMembers-selection.selectedMembers));
			summary.println("next_family_id\t"+selection.nextFamilyId);
		}finally{close(models,mapping,summary);}
	}

	static final class Selection{
		Selection(final ClusterSketchCodec.Loaded models_,final Seed[] selected_,final int next_,final long source_,
				final long kept_,final int inputs_,final int budget_,final double fraction_,final int requested_){
			assert(models_!=null && selected_.length==models_.families.length && next_>=selected_.length) :
				"Resident selection preserves the codec roster and unused global ID cursor";
			models=models_;selected=selected_;nextFamilyId=next_;sourceMembers=source_;selectedMembers=kept_;
			inputCodecs=inputs_;budget=budget_;fraction=fraction_;requested=requested_;
		}
		final ClusterSketchCodec.Loaded models;
		final int nextFamilyId,inputCodecs,budget,requested;
		final long sourceMembers,selectedMembers;
		final double fraction;
		private final Seed[] selected;
	}

	/** Decimal parameter arithmetic avoids a binary-rounding off-by-one before flooring. */
	public static int requestedCount(final int familyBudget,final double initialFraction){
		require(familyBudget>0 && Double.isFinite(initialFraction) && initialFraction>0,
			"Initial selection requires positive familybudget and finite positive initialfraction");
		final BigDecimal count=BigDecimal.valueOf(initialFraction).multiply(BigDecimal.valueOf(familyBudget))
			.setScale(0,RoundingMode.FLOOR);
		require(count.signum()>0 && count.compareTo(BigDecimal.valueOf(Integer.MAX_VALUE))<=0,
			"Initial roster size must be a positive bounded integer after flooring initialfraction times familybudget");
		return count.intValueExact();
	}

	/** All shards must describe the same effective model binding. */
	private static void sameGeometry(final ClusterSketchCodec.Header a,final ClusterSketchCodec.Header b){
		require(a.sameModelBinding(b),
			"Global selection cannot mix codec model bindings");
	}

	private static ByteStreamWriter open(final String path,final boolean overwrite){
		final ByteStreamWriter writer=new ByteStreamWriter(path,overwrite,false,false);writer.start();return writer;
	}

	/** Attempt every close even if an earlier writer failed. */
	private static void close(final ByteStreamWriter...writers){
		RuntimeException failure=null;
		for(ByteStreamWriter writer:writers){
			if(writer==null){continue;}
			try{if(writer.poisonAndWait()){throw new IllegalStateException("Global seed selection output failed");}}
			catch(RuntimeException e){if(failure==null){failure=e;}else{failure.addSuppressed(e);}}
		}
		if(failure!=null){throw failure;}
	}

	private static void require(final boolean condition,final String reason){
		if(!condition){throw new IllegalArgumentException(reason);}
	}

	private static final class Seed{
		Seed(final int ordinal_,final int shard_,final ClusterSketchCodec.Entry entry_){
			assert(ordinal_>=0 && shard_>=0 && entry_!=null) : "Selection provenance requires an existing source entry";
			ordinal=ordinal_;shard=shard_;entry=entry_;
		}
		final int ordinal,shard;
		final ClusterSketchCodec.Entry entry;
	}
	private static final Comparator<Seed> WORST_FIRST=new Comparator<Seed>(){
		@Override public int compare(final Seed a,final Seed b){
			final int count=Integer.compare(a.entry.members,b.entry.members);
			return count!=0 ? count : Integer.compare(b.ordinal,a.ordinal);
		}
	};
	private static final Comparator<Seed> BEST_FIRST=new Comparator<Seed>(){
		@Override public int compare(final Seed a,final Seed b){return WORST_FIRST.compare(b,a);}
	};
	private ClusterSeedSelector(){}
}
