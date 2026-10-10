package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Reuses unchanged seed profiles and refines changed nonempty memberships through
 * the existing two-pass positional refiner. Equality is exact ID/encoded sequence.
 * The original background is passed explicitly and checked against each seed.
 * Empty families remain in the registry with a decision-required state; no model
 * is silently deleted or counted as a newly trained empty model.
 * @author Keqing
 */
public final class HbmReassignedProfiles {

	public static void main(String[] args){
		try{
			Shared.setThreads(1); Shared.AMINO_IN=true; Read.VALIDATE_IN_CONSTRUCTOR=false;
			final PreParser pp=new PreParser(args, HbmReassignedProfiles.class, false);
			try{run(HmmComparisonData.options(pp.args));}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	private static void run(HashMap<String,String> o) throws Exception{
		for(String key : o.keySet()){
			require(Arrays.asList("old", "oldsha80", "new", "newsha80", "models", "modelssha80", "background", "backgroundsha80", "runtime", "ranks", "out").contains(key), "Unknown reassigned-profile parameter");
		}
		for(String key : INPUTS){HbmProfilePilot.requireHash(HmmComparisonData.required(o, key), HmmComparisonData.required(o, key+"sha80"));}
		final List<String[]> oldRows=HmmComparisonData.rows(o.get("old")), newRows=HmmComparisonData.rows(o.get("new")), models=HmmComparisonData.rows(o.get("models"));
		require(oldRows.size()>1 && oldRows.size()==newRows.size() && oldRows.size()==models.size(), "Member/model roster sizes differ");
		require(String.join("\t", oldRows.get(0)).equals(OLD_HEADER) && String.join("\t", newRows.get(0)).equals(NEW_HEADER)
			&& String.join("\t", models.get(0)).equals(MODEL_HEADER), "Unsupported member/model manifests");
		final HashSet<String> ids=new HashSet<String>(), reps=new HashSet<String>();
		for(int i=1; i<oldRows.size(); i++){
			final String[] a=oldRows.get(i), b=newRows.get(i), h=models.get(i); final String rank=Integer.toString(i-1);
			require(a.length==7 && b.length==9 && h.length==6 && a[0].equals(rank) && b[0].equals(rank) && h[0].equals(rank), "Noncontiguous member/model rank");
			require(a[1].equals(b[1]) && a[2].equals(b[2]) && a[2].equals(h[1]) && a[6].equals(h[5]), "Seed/new family identity or seed-input binding differs");
			final int id=Integer.parseInt(a[1]); final long count=Long.parseLong(b[3]), raw=Long.parseLong(b[4]), encoded=Long.parseLong(b[5]);
			require(id>=0 && a[1].equals(Integer.toString(id)) && ids.add(a[1]) && reps.add(a[2]) && Long.parseLong(a[3])>0, "Invalid permanent identity or old membership count");
			require(count>=0 && raw>=encoded && encoded>=count && (count>0 || raw==0) && b[8].equals(count==0 ? "EMPTY" : "NONEMPTY"), "Invalid reassigned count/residue/empty state");
			require(Paths.get(a[5]).isAbsolute() && Paths.get(b[6]).isAbsolute() && Paths.get(h[2]).isAbsolute(), "Member/model paths must be absolute");
		}
		final double[] background=readBackground(o.get("background"));
		final int[] ranks=HbmProfileLibrary.readRanks(o.get("ranks"), oldRows.size()-1); Arrays.sort(ranks);
		final Path out=Paths.get(HmmComparisonData.required(o, "out")); require(!Files.exists(out), "Reassigned-profile output must be fresh"); Files.createDirectory(out);
		final ByteStreamWriter registry=HmmComparisonData.writer(out.resolve("families.tsv").toString());
		registry.println("active_index\tfamily_id\trep_id\told_members\tnew_members\tderived_raw_residues\tencoded_residues\tadded_ids\tremoved_ids\tchanged_sequences\tstate\tdirectory\tbundle_sha80\tconsensus_sha80\tnew_input_sha80");
		final ByteBuilder row=new ByteBuilder(); int rebuilt=0, reused=0, empty=0;
		try{
			for(int rank : ranks){
				final String[] a=oldRows.get(rank+1), b=newRows.get(rank+1), h=models.get(rank+1);
				final Path seed=Paths.get(h[2]);
				HbmProfilePilot.requireHash(seed.resolve("refined.mqhb").toString(), h[3]);
				HbmProfilePilot.requireHash(seed.resolve("consensus.faa").toString(), h[4]);
				HbmProfilePilot.requireHash(seed.resolve("background.tsv").toString(), o.get("backgroundsha80"));
				final ProteinSequence reference=readReference(seed.resolve("consensus.faa").toString());
				require(reference.id.equals(a[2]), "Seed consensus identity differs");
				final List<ProteinSequence> oldMembers=HbmMembershipCompare.readMembers(a[5], a[6], Long.parseLong(a[3]), Long.parseLong(a[4]), -1);
				final List<ProteinSequence> newMembers=HbmMembershipCompare.readMembers(b[6], b[7], Long.parseLong(b[3]), Long.parseLong(b[4]), Long.parseLong(b[5]));
				final long[] delta=HbmMembershipCompare.compare(oldMembers, newMembers);
				oldMembers.clear(); //Only the new encoded members are needed during refinement.
				String state; Path effective=seed;
				if(newMembers.isEmpty()){state="EMPTY_REQUIRES_DECISION"; empty++;}
				else if(delta[0]+delta[1]+delta[2]==0){state="REUSED"; reused++;}
				else{
					state="REBUILT"; rebuilt++; effective=out.resolve("rank_"+rank);
					refine(rank, a, b, h, reference, newMembers, background, o, effective);
				}
				HbmProfilePilot.requireHash(seed.resolve("refined.mqhb").toString(), h[3]); HbmProfilePilot.requireHash(seed.resolve("consensus.faa").toString(), h[4]);
				HbmProfilePilot.requireHash(seed.resolve("background.tsv").toString(), o.get("backgroundsha80"));
				HbmProfilePilot.requireHash(b[6], b[7]);
				registry.println(row.clear().append(rank).tab().append(a[1]).tab().append(a[2]).tab().append(a[3]).tab().append(b[3]).tab().append(b[4]).tab().append(b[5])
					.tab().append(delta[0]).tab().append(delta[1]).tab().append(delta[2]).tab().append(state).tab().append(effective.toString())
					.tab().append(DigestSuffix.file(effective.resolve("refined.mqhb").toString())).tab().append(DigestSuffix.file(effective.resolve("consensus.faa").toString())).tab().append(b[7]));
				System.err.println("HBM_REASSIGNED_PROFILE_PROGRESS rank="+rank+" state="+state);
			}
		}finally{HmmComparisonData.close(registry);}
		for(String key : INPUTS){HbmProfilePilot.requireHash(o.get(key), o.get(key+"sha80"));}
		final ByteStreamWriter summary=HmmComparisonData.writer(out.resolve("summary.tsv").toString());
		summary.println("families\trebuilt\treused\tempty_requires_decision"); summary.println(row.clear().append(ranks.length).tab().append(rebuilt).tab().append(reused).tab().append(empty));
		HmmComparisonData.close(summary);
		final ByteStreamWriter recipe=HmmComparisonData.writer(out.resolve("recipe.tsv").toString());
		for(String key : INPUTS){recipe.println(key+"_sha80\t"+o.get(key+"sha80"));}
		recipe.println("membership_equality\texact_id_and_encoded_sequence\nbackground_policy\tfrozen_original\nempty_policy\tdecision_required_parent_model_preserved");
		HmmComparisonData.close(recipe);
		final ByteStreamWriter pass=HmmComparisonData.writer(out.resolve("PASS").toString());
		pass.println("HBM_REASSIGNED_PROFILES_PASS final_composition_pending=true"); HmmComparisonData.close(pass);
	}

	/** Refines one changed family using an immutable singleton seed and sorted new members. */
	private static void refine(int rank, String[] old, String[] membersRow, String[] model, ProteinSequence reference, List<ProteinSequence> proteins,
			double[] background, HashMap<String,String> o, Path out) throws Exception{
		final Path seed=Paths.get(model[2]);
		final HbmBundleLoader.Loaded loaded=HbmBundleLoader.load(seed.resolve("refined.mqhb"), Collections.singletonList(reference.id), id->reference.enc,
			HbmBundleLoader.loadSemanticProvenance(seed.resolve("provenance.tsv").toString()));
		final byte[][] members=new byte[proteins.size()][]; int maximum=0;
		for(int i=0; i<members.length; i++){members[i]=proteins.get(i).enc; maximum=Math.max(maximum, members[i].length);}
		HbmProfileRefiner.Result result; int padding=20, retries=0;
		while(true){
			try{result=loaded.refineFamily(0, members, background, padding); break;}
			catch(HbmProfileRefiner.PaddingException small){
				if(padding>=maximum){throw small;}
				padding=(int)Math.min(maximum, Math.max(2L*padding, (long)padding+small.excluded)); retries++;
			}
		}
		require(result.members==members.length && result.residues==Long.parseLong(membersRow[5]), "Refined graph member/residue accounting differs");
		Files.createDirectory(out); Files.copy(Paths.get(o.get("background")), out.resolve("background.tsv"));
		final ByteStreamWriter fasta=HmmComparisonData.writer(out.resolve("consensus.faa").toString());
		HmmComparisonData.writeFasta(fasta, reference.id, result.graph.pivot); HmmComparisonData.close(fasta);
		final ByteStreamWriter roster=HmmComparisonData.writer(out.resolve("roster.tsv").toString());
		roster.println("rank\trep_id\tmembers"); roster.println(new ByteBuilder().append(0).tab().append(reference.id).tab().append(members.length)); HmmComparisonData.close(roster);
		final String runtime=HmmComparisonData.required(o, "runtime");
		final ByteStreamWriter recipe=HmmComparisonData.writer(out.resolve("recipe.txt").toString());
		recipe.println("experimental_reassigned_profile_v1\nbeta=0.01\nclipmin=-4\nclipmax=11\ngap=4\ntrimdepth=0.1\npasses=2\nfinal_padding=0\nmember_order=ascii_query_id\nfinal_terminal_policy=count_exclusions\nbackground_policy=frozen_original");
		recipe.println("source_rank="+rank+"\nfamily_id="+old[1]+"\nsource_rep_id="+reference.id+"\ngrowth_padding="+padding+"\nparent_hbm_sha80="+model[3]
			+"\nold_input_sha80="+old[6]+"\nnew_input_sha80="+membersRow[7]+"\nbackground_sha80="+o.get("backgroundsha80"));
		for(String name : new String[]{"HbmProfileRefiner", "HbmPositionModel", "HbmMembershipCompare"}){
			recipe.println(name+"_source_sha80="+DigestSuffix.file(runtime+"/current/prot/"+name+".java"));
			recipe.println(name+"_class_sha80="+DigestSuffix.file(runtime+"/current/prot/"+name+".class"));
		}
		HmmComparisonData.close(recipe);
		final byte[][] provenance=HbmProfilePilot.writeProvenance(out, runtime, membersRow[6], membersRow[6], o.get("new"), "HbmReassignedProfiles");
		HbmBundleBuilder.build(out.resolve("refined.mqhb"), Collections.singletonList(new HbmBundleBuilder.FamilyInput(reference.id, result.graph.pivot, result.graph)), provenance);
		final byte[] consensus=result.graph.pivot;
		final HbmBundleLoader.Loaded roundtrip=HbmBundleLoader.load(out.resolve("refined.mqhb"), Collections.singletonList(reference.id), id->consensus,
			HbmBundleLoader.loadSemanticProvenance(out.resolve("provenance.tsv").toString()));
		roundtrip.assertStructuralMatch(0, result.graph);
		HbmProfilePilot.requireHash(out.resolve("background.tsv").toString(), o.get("backgroundsha80"));
		final ByteStreamWriter summary=HmmComparisonData.writer(out.resolve("summary.tsv").toString());
		summary.println("active_index\tfamily_id\trep_id\tmembers\tderived_raw_residues\tencoded_residues\told_length\tnew_length\tpadding\tpadding_retries\tfinal_terminal_excluded\tbundle_sha80\tconsensus_sha80");
		summary.println(new ByteBuilder().append(rank).tab().append(old[1]).tab().append(reference.id).tab().append(members.length).tab().append(membersRow[4]).tab().append(result.residues)
			.tab().append(reference.length()).tab().append(consensus.length).tab().append(padding).tab().append(retries).tab().append(result.excludedTerminalResidues)
			.tab().append(DigestSuffix.file(out.resolve("refined.mqhb").toString())).tab().append(DigestSuffix.file(out.resolve("consensus.faa").toString())));
		HmmComparisonData.close(summary);
		final ByteStreamWriter pass=HmmComparisonData.writer(out.resolve("PASS").toString()); pass.println("HBM_REASSIGNED_FAMILY_PASS"); HmmComparisonData.close(pass);
	}

	/** Reads exactly one strict seed consensus; malformed records cannot be skipped. */
	private static ProteinSequence readReference(String file){
		final Streamer input=StreamerFactory.makeStreamer(FileFormat.testInput(file, FileFormat.FASTA, null, true, false), null, true, -1);
		ProteinSequence reference=null; input.start();
		try{
			for(ListNum<Read> list=input.nextList(); list!=null && list.size()>0; list=input.nextList()){
				for(Read read : list){
					require(reference==null && read.mate==null, "Expected exactly one unpaired seed consensus");
					reference=new ProteinSequence(HbmCompetitiveAssign.queryId(read.id), read.bases);
				}
			}
		}finally{input.close(); require(!input.errorState(), "Seed consensus reader failed");}
		require(reference!=null, "Empty seed consensus"); return reference;
	}
	private static double[] readBackground(String file){
		final List<String[]> rows=HmmComparisonData.rows(file);
		require(rows.size()==21 && String.join("\t", rows.get(0)).equals("residue_code\tprobability"), "Invalid background schema");
		final double[] out=new double[20]; double sum=0;
		for(int i=0; i<20; i++){
			final String[] r=rows.get(i+1); require(r.length==2 && r[0].equals(Integer.toString(i)), "Background residue order differs");
			out[i]=Double.parseDouble(r[1]); require(Double.isFinite(out[i]) && out[i]>0 && out[i]<=1, "Invalid background probability"); sum+=out[i];
		}
		require(Math.abs(sum-1)<1e-9, "Background probabilities must sum to one"); return out;
	}
	private static void require(boolean ok, String reason){if(!ok){throw new IllegalArgumentException(reason);}}
	private static final String[] INPUTS={"old", "new", "models", "background"};
	private static final String OLD_HEADER="active_index\tfamily_id\trep_id\tassigned_count\tsequence_bytes\tfasta_file\tfasta_sha80";
	private static final String NEW_HEADER="active_index\tfamily_id\trep_id\tassigned_count\tderived_raw_residues\tencoded_residues\tfasta_file\tfasta_sha80\tmembership_status";
	private static final String MODEL_HEADER="active_index\trep_id\tdirectory\tbundle_sha80\tconsensus_sha80\tinput_sha80";
}
