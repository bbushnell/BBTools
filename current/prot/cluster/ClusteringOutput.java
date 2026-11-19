package prot.cluster;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Deterministic output in input order or final size-rank order, with stable family IDs.
 * @author Keqing
 */
public final class ClusteringOutput {

	private ClusteringOutput(){}

	public static void members(String path, boolean overwrite, ClusteringInput input, GreedyClusterer.Result result, int minLength){
		if(path==null){return;}
		final ByteStreamWriter writer=open(path, overwrite);
		try{
			writer.println("input_index\tinput_id\tfamily_id\tstatus");
			final ByteBuilder row=new ByteBuilder();
			for(int i=0; i<input.sequences.length; i++){
				final int family=result.familyFor(i);
				final String status=input.sequences[i].length()<minLength ? "filtered_length" : family<0 ? "unassigned" : "assigned";
				row.clear().append(i).tab().append(input.sequences[i].id).tab().append(family).tab().append(status);
				writer.println(row);
			}
		}finally{close(writer, path);}
	}

	public static void consensus(String path, boolean overwrite, GreedyClusterer.Result result){
		if(path==null){return;}
		final ByteStreamWriter writer=open(path, overwrite);
		try{
			final ByteBuilder row=new ByteBuilder();
			for(int rank=0; rank<result.families.length; rank++){
				final Family family=result.families[rank];
				row.clear().append('>').append(family.id).append(" members=").append(family.size())
					.append(" rank=").append(rank+1).nl().append(ClusteringInput.decode(family.consensusUnsafe())).nl();
				writer.print(row);
			}
		}finally{close(writer, path);}
	}

	/** One family file at a time bounds open handles even for a large final roster. */
	public static void clusters(String pattern, boolean overwrite, ClusteringInput input, GreedyClusterer.Result result){
		if(pattern==null){return;}
		if(pattern.indexOf('#')<0){throw new IllegalArgumentException("Cluster FASTA output requires # for its stable family ID");}
		for(Family family : result.families){
			final String path=pattern.replace("#", Integer.toString(family.familyId));
			final ByteStreamWriter writer=open(path, overwrite);
			try{
				final int[] members=family.memberIndexes(); java.util.Arrays.sort(members);
				final ByteBuilder row=new ByteBuilder();
				for(int member : members){row.clear(); input.appendFasta(member, row); writer.print(row);}
			}finally{close(writer, path);}
		}
	}

	public static void summary(String path, boolean overwrite, GreedyClusterer.Result result, long seed,
			int passes, double seconds, ClusteringMetrics metrics, ConsensusRefiner refiner, String[] configuration){
		if(path==null){return;}
		if(!Double.isFinite(seconds) || seconds<0){throw new IllegalArgumentException("Elapsed time must be finite and nonnegative");}
		final ByteStreamWriter writer=open(path, overwrite);
		try{
			writer.println("key\tvalue");
			line(writer, "input_records", result.inputs()); line(writer, "eligible_records", result.eligible);
			line(writer, "filtered_records", result.inputs()-result.eligible); line(writer, "assigned_records", result.assigned);
			line(writer, "unassigned_records", result.eligible-result.assigned); line(writer, "families", result.families.length);
			line(writer, "seed", seed); line(writer, "passes", passes); line(writer, "completed_passes", result.completedPasses);
			line(writer, "alignments", metrics.totalAlignments());
			line(writer, "sketch_updates", metrics.sketchUpdates()); line(writer, "consensus_updates", metrics.consensusUpdates());
			line(writer, "assignment_alignments", metrics.assignmentAlignments());
			line(writer, "construction_alignments", metrics.constructionAlignments());
			line(writer, "hit_filtered_candidates", metrics.hitFiltered());
			line(writer, "hydro_filtered_candidates", metrics.hydroFiltered());
			line(writer, "consensus_proposals", refiner.proposals()); line(writer, "consensus_proposals_rejected", refiner.rejected());
			line(writer, "consensus_proposals_unchanged", refiner.unchanged());
			line(writer, "consensus_rejected_member_checks", refiner.rejectedMemberChecks());
			writer.println(new ByteBuilder().append("elapsed_seconds\t").append(seconds, 6));
			for(String row : configuration){writer.println(row);}
		}finally{close(writer, path);}
	}

	private static void line(ByteStreamWriter writer, String key, long value){
		assert(key.indexOf('\t')<0) : "Summary keys form exactly one TSV field";
		writer.println(new ByteBuilder().append(key).tab().append(value));
	}

	private static ByteStreamWriter open(String path, boolean overwrite){
		final ByteStreamWriter writer=new ByteStreamWriter(path, overwrite, false, true);
		writer.start(); return writer;
	}

	private static void close(ByteStreamWriter writer, String path){
		if(writer.poisonAndWait()){throw new IllegalStateException("Failed clustering output: "+path);}
	}
}
