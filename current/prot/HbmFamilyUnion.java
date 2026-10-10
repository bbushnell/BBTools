package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import fileIO.ByteStreamWriter;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Composes verified old refinements and new seed profiles without changing their graphs.
 * The dense output order is old selections followed by new selections; permanent
 * family IDs and source indexes are retained separately. The common position-score
 * background is frozen, not recalculated from the enlarged model collection.
 * @author Brian Bushnell, Keqing
 */
public final class HbmFamilyUnion {

	public static void main(String[] args){
		try{
			Shared.setThreads(1); Shared.AMINO_IN=true;
			PreParser pp=new PreParser(args, HbmFamilyUnion.class, false);
			try{run(HmmComparisonData.options(pp.args));}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	private static void run(HashMap<String,String> o) throws Exception{
		for(String key : o.keySet()){
			require(Arrays.asList("resources", "oldmanifest", "oldmanifestsha80", "oldsource", "oldruntime", "oldfamilies",
				"newmanifest", "newmanifestsha80", "newsource", "newruntime", "newfamilies", "runtime", "out",
				"background", "backgroundsha80", "oldranks", "newranks").contains(key), "Unknown union parameter: "+key);
		}
		final Path out=Paths.get(HmmComparisonData.required(o, "out"));
		require(!Files.exists(out), "Union output must be a fresh directory: "+out);
		final HashMap<String,String> previous=new HashMap<String,String>();
		previous.put("resources", HmmComparisonData.required(o, "resources"));
		for(String key : new String[]{"manifest", "manifestsha80", "runtime"}){previous.put(key, HmmComparisonData.required(o, "old"+key));}
		final HbmProfilePilot.Session old=new HbmProfilePilot.Session(previous, true);
		final String manifest=HmmComparisonData.required(o, "newmanifest"), manifestPin=HmmComparisonData.required(o, "newmanifestsha80");
		HbmProfilePilot.requireHash(manifest, manifestPin);
		final List<String[]> rows=HmmComparisonData.rows(manifest);
		require(rows.size()>1 && String.join("\t", rows.get(0)).equals(MANIFEST_HEADER), "Unsupported new-family manifest");
		final HashSet<String> ids=new HashSet<String>(), reps=new HashSet<String>();
		for(int i=1; i<old.rows.size(); i++){checkIdentity(old.rows.get(i), i-1, ids, reps);}
		for(int i=1; i<rows.size(); i++){
			final String[] row=rows.get(i); checkIdentity(row, old.roster.size()+i-1, ids, reps);
			require(row[5].equals("rank_"+row[0]+".faa"), "Unexpected new-family input filename at "+row[0]);
		}
		final String background=HmmComparisonData.required(o, "background"), backgroundPin=HmmComparisonData.required(o, "backgroundsha80");
		HbmProfilePilot.requireHash(background, backgroundPin);
		final ByteBuilder bg=new ByteBuilder().append("residue_code\tprobability\n");
		for(int i=0; i<20; i++){bg.append(i).tab().append(Double.toString(old.background[i])).nl();}
		require(Arrays.equals(bg.toBytes(), Files.readAllBytes(Paths.get(background))), "Union background differs from the old frozen library background");
		final int[] oldRanks=HbmProfileLibrary.readRanks(o.get("oldranks"), old.roster.size());
		final int[] newRanks=HbmProfileLibrary.readRanks(o.get("newranks"), rows.size()-1);
		Arrays.sort(oldRanks); Arrays.sort(newRanks);
		final ArrayList<HbmBundleLoader.Loaded> models=new ArrayList<HbmBundleLoader.Loaded>();
		final ArrayList<String[]> entries=new ArrayList<String[]>();
		final ArrayList<Path> directories=new ArrayList<Path>();
		for(int rank : oldRanks){
			final Path dir=Paths.get(HmmComparisonData.required(o, "oldfamilies"), String.format(java.util.Locale.ROOT, "rank_%04d", rank));
			models.add(HbmProfileLibrary.validateFamily(dir, rank, old, HmmComparisonData.required(o, "oldsource")));
			entries.add(old.rows.get(rank+1)); directories.add(dir);
		}
		for(int rank : newRanks){
			final String[] row=rows.get(rank+1);
			final Path dir=Paths.get(HmmComparisonData.required(o, "newfamilies"), "rank_"+row[0]);
			models.add(validateNew(dir, row, o, manifest, manifestPin, backgroundPin));
			entries.add(row); directories.add(dir);
		}
		old.verify(); HbmProfilePilot.requireHash(manifest, manifestPin); HbmProfilePilot.requireHash(background, backgroundPin);
		final ArrayList<String> roster=new ArrayList<String>();
		final HashMap<String,byte[]> sequences=new HashMap<String,byte[]>();
		for(int i=0; i<entries.size(); i++){
			final String rep=entries.get(i)[2]; roster.add(rep);
			sequences.put(rep, ProteinSearch.readFasta(directories.get(i).resolve("consensus.faa").toString()).get(0).enc);
		}
		final HbmBundleLoader.Loaded combined=HbmBundleLoader.Loaded.combine(models, roster);
		Files.createDirectory(out);
		final ByteStreamWriter fasta=HmmComparisonData.writer(out.resolve("consensus.faa").toString());
		final ByteStreamWriter order=HmmComparisonData.writer(out.resolve("roster.tsv").toString());
		final ByteStreamWriter identity=HmmComparisonData.writer(out.resolve("family_identity.tsv").toString());
		final ByteStreamWriter artifacts=HmmComparisonData.writer(out.resolve("family_artifacts.tsv").toString());
		order.println("rank\trep_id\tmembers");
		identity.println("active_index\tfamily_id\tsource_active_index\trep_id\tmembers\traw_residues");
		artifacts.println("active_index\trep_id\tdirectory\tbundle_sha80\tconsensus_sha80\tinput_sha80");
		long members=0, residues=0;
		for(int i=0; i<entries.size(); i++){
			final String[] row=entries.get(i); final Path dir=directories.get(i);
			HmmComparisonData.writeFasta(fasta, row[2]+" active_index="+i+" family_id="+row[1], sequences.get(row[2]));
			order.println(new ByteBuilder().append(i).tab().append(row[2]).tab().append(row[3]));
			identity.println(new ByteBuilder().append(i).tab().append(row[1]).tab().append(row[0]).tab().append(row[2]).tab().append(row[3]).tab().append(row[4]));
			artifacts.println(new ByteBuilder().append(i).tab().append(row[2]).tab().append(dir.toString())
				.tab().append(DigestSuffix.file(dir.resolve("refined.mqhb").toString())).tab().append(DigestSuffix.file(dir.resolve("consensus.faa").toString())).tab().append(row[6]));
			members=Math.addExact(members, Long.parseLong(row[3])); residues=Math.addExact(residues, Long.parseLong(row[4]));
		}
		HmmComparisonData.close(fasta); HmmComparisonData.close(order); HmmComparisonData.close(identity); HmmComparisonData.close(artifacts);
		Files.copy(Paths.get(background), out.resolve("background.tsv"));
		final ByteStreamWriter recipe=HmmComparisonData.writer(out.resolve("recipe.txt").toString());
		recipe.println("experimental_hbm_family_union_v1\noperation=compose_unchanged_singletons\nbeta=0.01\nclipmin=-4\nclipmax=11\ngap=4\nbackground_policy=frozen_original_library");
		recipe.println("old_manifest_sha80="+old.manifestPin+"\nnew_manifest_sha80="+manifestPin+"\nbackground_sha80="+backgroundPin);
		recipe.println("family_identity_sha80="+DigestSuffix.file(out.resolve("family_identity.tsv").toString()));
		recipe.println("family_artifacts_sha80="+DigestSuffix.file(out.resolve("family_artifacts.tsv").toString()));
		HmmComparisonData.close(recipe);
		final String runtime=HmmComparisonData.required(o, "runtime");
		final byte[][] provenance=HbmProfilePilot.writeProvenance(out, runtime, out.resolve("family_artifacts.tsv").toString(),
			out.resolve("family_identity.tsv").toString(), manifest, "HbmFamilyUnion");
		combined.writeBundle(out.resolve("refined.mqhb"), provenance);
		final HbmBundleLoader.Loaded loaded=HbmBundleLoader.load(out.resolve("refined.mqhb"), roster, id->sequences.get(id),
			HbmBundleLoader.loadSemanticProvenance(out.resolve("provenance.tsv").toString()));
		loaded.assertStructuralEquivalent(combined);
		require(Arrays.equals(FamilyShortlistSidecarBuilder.loadConsensusFamilyOrder(out.resolve("consensus.faa").toString()), roster.toArray(new String[0])),
			"Sidecar consumer saw a different union order");
		final ByteStreamWriter summary=HmmComparisonData.writer(out.resolve("summary.tsv").toString());
		summary.println("old_families\tnew_families\tfamilies\tmembers\traw_residues\tbundle_sha80\tconsensus_sha80\tbackground_sha80");
		summary.println(new ByteBuilder().append(oldRanks.length).tab().append(newRanks.length).tab().append(roster.size()).tab().append(members).tab().append(residues)
			.tab().append(DigestSuffix.file(out.resolve("refined.mqhb").toString())).tab().append(DigestSuffix.file(out.resolve("consensus.faa").toString())).tab().append(backgroundPin));
		HmmComparisonData.close(summary);
		final ByteStreamWriter pass=HmmComparisonData.writer(out.resolve("PASS").toString()); pass.println("HBM_FAMILY_UNION_PASS"); HmmComparisonData.close(pass);
		System.err.println("HBM_FAMILY_UNION_PASS families="+roster.size()+" members="+members);
	}

	/** Rejects identity aliases because dense output indexes must retain unambiguous permanent IDs. */
	private static void checkIdentity(String[] row, int active, HashSet<String> ids, HashSet<String> reps){
		require(row.length==7 && row[0].equals(Integer.toString(active)), "Noncontiguous or malformed member manifest at "+active);
		final int id=Integer.parseInt(row[1]);
		require(id>=0 && row[1].equals(Integer.toString(id)) && ids.add(row[1]) && !row[2].isEmpty() && reps.add(row[2]),
			"Duplicate or invalid permanent family identity at "+active);
		require(Long.parseLong(row[3])>0 && Long.parseLong(row[4])>0, "Empty family membership at "+active);
		DigestSuffix.requireSuffix(row[6], "member input sha80");
	}

	/** Verifies the complete seed-builder output contract before accepting a new singleton. */
	private static HbmBundleLoader.Loaded validateNew(Path dir, String[] entry, HashMap<String,String> o,
		String manifest, String manifestPin, String backgroundPin) throws Exception{
		require(new String(Files.readAllBytes(dir.resolve("PASS")), StandardCharsets.UTF_8).trim().equals("NEW_FAMILY_HBM_PASS"), "Incomplete new singleton: "+dir);
		final List<String[]> rows=HmmComparisonData.rows(dir.resolve("summary.tsv").toString());
		require(rows.size()==2 && String.join("\t", rows.get(0)).equals(NEW_SUMMARY) && rows.get(1).length==12, "Unexpected new singleton summary: "+dir);
		final String[] row=rows.get(1);
		for(int i=0; i<5; i++){require(row[i].equals(entry[i]), "New singleton identity/count differs at field "+i+": "+dir);}
		final long encoded=Long.parseLong(row[5]);
		require(encoded>0 && encoded<=Long.parseLong(row[4]) && Long.parseLong(row[8])>=0 && Long.parseLong(row[8])<=encoded
			&& Long.parseLong(row[9])>=0 && Long.parseLong(row[9])<=encoded, "Invalid encoded/terminal counts: "+dir);
		HbmProfilePilot.requireHash(dir.resolve("refined.mqhb").toString(), row[10]);
		HbmProfilePilot.requireHash(dir.resolve("consensus.faa").toString(), row[11]);
		HbmProfilePilot.requireHash(dir.resolve("background.tsv").toString(), backgroundPin);
		final List<ProteinSequence> sequences=ProteinSearch.readFasta(dir.resolve("consensus.faa").toString());
		require(sequences.size()==1 && sequences.get(0).id.equals(entry[2]) && sequences.get(0).length()==Integer.parseInt(row[7]), "New singleton consensus differs: "+dir);
		final List<String[]> roster=HmmComparisonData.rows(dir.resolve("roster.tsv").toString());
		require(roster.size()==2 && roster.get(1).length==3 && roster.get(1)[0].equals("0") && roster.get(1)[1].equals(entry[2]) && roster.get(1)[2].equals(entry[3]), "New singleton roster differs: "+dir);
		final List<String> lines=Files.readAllLines(dir.resolve("recipe.txt"), StandardCharsets.UTF_8);
		require(!lines.isEmpty() && lines.get(0).equals("experimental_new_family_hbm_v1"), "Wrong new singleton recipe");
		final HashMap<String,String> recipe=new HashMap<String,String>();
		for(int i=1; i<lines.size(); i++){
			final int eq=lines.get(i).indexOf('=');
			require(eq>0 && recipe.put(lines.get(i).substring(0, eq), lines.get(i).substring(eq+1))==null, "Malformed or duplicate recipe field: "+dir);
		}
		for(String setting : new String[]{"seed=BLOSUM62_glocal_gap4_original_representative", "membership=all_no_cap", "beta=0.01", "clipmin=-4", "clipmax=11", "gap=4", "trimdepth=0.1", "profile_passes=2", "final_padding=0", "terminal_policy=count_exclusions",
			"source_active_index="+entry[0], "family_id="+entry[1], "rep_id="+entry[2], "input_sha80="+entry[6], "manifest_sha80="+manifestPin, "background_sha80="+backgroundPin}){
			final int eq=setting.indexOf('='); require(setting.substring(eq+1).equals(recipe.get(setting.substring(0, eq))), "New recipe setting differs: "+setting);
		}
		final String runtime=HmmComparisonData.required(o, "newruntime");
		for(String dependency : new String[]{"HbmPositionModel", "HbmProfileRefiner", "HbmProfilePilot"}){
			require(DigestSuffix.file(runtime+"/current/prot/"+dependency+".java").equals(recipe.get(dependency+"_source_sha80")), "New source dependency differs: "+dependency);
			require(DigestSuffix.file(runtime+"/current/prot/"+dependency+".class").equals(recipe.get(dependency+"_class_sha80")), "New class dependency differs: "+dependency);
		}
		final String input=Paths.get(HmmComparisonData.required(o, "newsource"), entry[5]).toString();
		HbmProfilePilot.requireHash(input, entry[6]);
		HbmProfilePilot.verifyProvenance(dir, runtime, input, input, manifest, "NewFamilyHbmBuilder");
		return HbmBundleLoader.load(dir.resolve("refined.mqhb"), Collections.singletonList(entry[2]), id->sequences.get(0).enc,
			HbmBundleLoader.loadSemanticProvenance(dir.resolve("provenance.tsv").toString()));
	}

	private static void require(boolean ok, String message){if(!ok){throw new IllegalArgumentException(message);}}
	private static final String MANIFEST_HEADER="active_index\tfamily_id\trep_id\tassigned_count\tsequence_bytes\tfasta_file\tfasta_sha80";
	private static final String NEW_SUMMARY="active_index\tfamily_id\trep_id\tmembers\traw_residue_bytes\tencoded_residues\tseed_length\tfinal_length\tseed_terminal_excluded\tfinal_terminal_excluded\tbundle_sha80\tconsensus_sha80";
}
