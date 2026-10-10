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
import fileIO.FileFormat;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Composes a complete reassigned-profile registry after checking every family
 * against independent membership and seed manifests. Unchanged and rebuilt
 * singletons are loaded natively; composition never edits their graphs.
 * Empty families stop composition explicitly until their disposition is decided.
 * @author Keqing
 */
public final class HbmReassignedUnion {

	public static void main(String[] args){
		try{
			Shared.setThreads(1); Shared.AMINO_IN=true; Read.VALIDATE_IN_CONSTRUCTOR=false;
			final PreParser pp=new PreParser(args, HbmReassignedUnion.class, false);
			try{run(HmmComparisonData.options(pp.args));}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	private static void run(HashMap<String,String> o) throws Exception{
		for(String key : o.keySet()){
			require(Arrays.asList("old", "oldsha80", "new", "newsha80", "models", "modelssha80", "shards", "shardssha80", "background", "backgroundsha80", "producer", "runtime", "out").contains(key), "Unknown reassigned-union parameter: "+key);
		}
		verifyInputs(o);
		final List<String[]> previous=rows(o.get("old"), OLD_HEADER), members=rows(o.get("new"), NEW_HEADER), seeds=rows(o.get("models"), MODEL_HEADER);
		final int count=members.size()-1;
		require(count>0 && previous.size()==members.size() && seeds.size()==members.size(), "Union member/seed roster counts differ");
		final HashSet<String> ids=new HashSet<String>(), reps=new HashSet<String>();
		for(int rank=0; rank<count; rank++){
			final String[] a=previous.get(rank+1), b=members.get(rank+1), h=seeds.get(rank+1); final String number=Integer.toString(rank);
			require(a.length==7 && b.length==9 && h.length==6 && a[0].equals(number) && b[0].equals(number) && h[0].equals(number), "Union ranks must be complete, unique and contiguous: "+rank);
			require(a[1].equals(b[1]) && a[2].equals(b[2]) && a[2].equals(h[1]) && a[6].equals(h[5]), "Union seed/member identity or input binding differs: "+rank);
			final int permanent=Integer.parseInt(b[1]); final long n=number(b[3]), raw=number(b[4]), encoded=number(b[5]);
			require(permanent>=0 && b[1].equals(Integer.toString(permanent)) && ids.add(b[1]) && !b[2].isEmpty() && reps.add(b[2]), "Duplicate or invalid permanent identity: "+rank);
			require(number(a[3])>0 && raw>=encoded && encoded>=n && (n>0 || raw==0) && b[8].equals(n==0 ? "EMPTY" : "NONEMPTY"), "Invalid family membership totals: "+rank);
			require(Paths.get(b[6]).isAbsolute() && Paths.get(h[2]).isAbsolute(), "Union input paths must be absolute");
			DigestSuffix.requireSuffix(b[7], "new input"); DigestSuffix.requireSuffix(h[3], "seed model"); DigestSuffix.requireSuffix(h[4], "seed consensus");
		}
		final String[][] registry=new String[count][];
		final List<String[]> shards=rows(o.get("shards"), SHARD_HEADER);
		final HashSet<Path> directories=new HashSet<Path>();
		for(int s=1; s<shards.size(); s++){
			final String[] shard=shards.get(s); require(shard.length==3, "Invalid profile shard manifest");
			final Path dir=Paths.get(shard[0]).normalize(); require(dir.isAbsolute() && directories.add(dir), "Duplicate or relative profile shard");
			require(text(dir.resolve("PASS")).equals("HBM_REASSIGNED_PROFILES_PASS final_composition_pending=true"), "Incomplete profile shard: "+dir);
			HbmProfilePilot.requireHash(dir.resolve("families.tsv").toString(), shard[1]);
			HbmProfilePilot.requireHash(dir.resolve("recipe.tsv").toString(), shard[2]);
			final HashMap<String,String> recipe=fields(dir.resolve("recipe.tsv"), '\t', false);
			for(String key : new String[]{"old", "new", "models", "background"}){
				require(o.get(key+"sha80").equals(recipe.get(key+"_sha80")), "Profile shard used different "+key+": "+dir);
			}
			require("exact_id_and_encoded_sequence".equals(recipe.get("membership_equality")) && "frozen_original".equals(recipe.get("background_policy"))
				&& "decision_required_parent_model_preserved".equals(recipe.get("empty_policy")), "Profile shard policy differs: "+dir);
			final List<String[]> familyRows=rows(dir.resolve("families.tsv").toString(), REGISTRY_HEADER);
			for(int i=1; i<familyRows.size(); i++){
				final String[] r=familyRows.get(i); require(r.length==15, "Invalid profile registry row"); final int rank=Integer.parseInt(r[0]);
				require(rank>=0 && rank<count && r[0].equals(Integer.toString(rank)) && registry[rank]==null, "Duplicate or invalid profile rank: "+r[0]);
				validateRegistry(r, previous.get(rank+1), members.get(rank+1), seeds.get(rank+1), dir);
				registry[rank]=r;
			}
		}
		int empty=0, rebuilt=0;
		for(int rank=0; rank<count; rank++){
			require(registry[rank]!=null, "Missing profile rank: "+rank);
			if(registry[rank][10].equals("EMPTY_REQUIRES_DECISION")){empty++;}
			else if(registry[rank][10].equals("REBUILT")){rebuilt++;}
		}
		require(empty==0, "Empty family policy is unresolved; composition would lose or invent membership. empty_families="+empty);
		final ArrayList<HbmBundleLoader.Loaded> loaded=new ArrayList<HbmBundleLoader.Loaded>();
		final ArrayList<String> roster=new ArrayList<String>(); final HashMap<String,byte[]> sequences=new HashMap<String,byte[]>();
		for(int rank=0; rank<count; rank++){
			final String[] r=registry[rank], member=members.get(rank+1), seed=seeds.get(rank+1); final Path dir=Paths.get(r[11]);
			HbmProfilePilot.requireHash(dir.resolve("refined.mqhb").toString(), r[12]); HbmProfilePilot.requireHash(dir.resolve("consensus.faa").toString(), r[13]);
			HbmProfilePilot.requireHash(dir.resolve("background.tsv").toString(), o.get("backgroundsha80"));
			final ProteinSequence ref=reference(dir.resolve("consensus.faa").toString()); require(ref.id.equals(r[2]), "Actual model consensus identity differs: "+rank);
			if(r[10].equals("REBUILT")){validateRebuild(dir, r, member, seed, o, ref.length());}
			else{HbmProfilePilot.requireHash(member[6], member[7]);}
			loaded.add(HbmBundleLoader.load(dir.resolve("refined.mqhb"), Collections.singletonList(r[2]), id->ref.enc,
				HbmBundleLoader.loadSemanticProvenance(dir.resolve("provenance.tsv").toString())));
			HbmProfilePilot.requireHash(dir.resolve("refined.mqhb").toString(), r[12]); HbmProfilePilot.requireHash(dir.resolve("consensus.faa").toString(), r[13]);
			roster.add(r[2]); sequences.put(r[2], ref.enc);
		}
		final HbmBundleLoader.Loaded combined=HbmBundleLoader.Loaded.combine(loaded, roster);
		final Path out=Paths.get(HmmComparisonData.required(o, "out")); require(!Files.exists(out), "Union output must be fresh"); Files.createDirectory(out);
		writeOutputs(out, members, registry, rebuilt, combined, roster, sequences, o);
		verifyInputs(o);
		for(int s=1; s<shards.size(); s++){
			HbmProfilePilot.requireHash(Paths.get(shards.get(s)[0], "families.tsv").toString(), shards.get(s)[1]);
			HbmProfilePilot.requireHash(Paths.get(shards.get(s)[0], "recipe.tsv").toString(), shards.get(s)[2]);
		}
		final ByteStreamWriter pass=HmmComparisonData.writer(out.resolve("PASS").toString()); pass.println("HBM_REASSIGNED_UNION_PASS"); HmmComparisonData.close(pass);
	}

	/** Cross-checks producer states/counts against independent old, new and seed manifests. */
	private static void validateRegistry(String[] r, String[] old, String[] member, String[] seed, Path shard){
		require(r[1].equals(member[1]) && r[2].equals(member[2]) && r[3].equals(old[3]) && r[4].equals(member[3]) && r[5].equals(member[4])
			&& r[6].equals(member[5]) && r[14].equals(member[7]), "Profile registry identity/count/input differs: "+r[0]);
		final long before=number(r[3]), after=number(r[4]), added=number(r[7]), removed=number(r[8]), changed=number(r[9]);
		require(added<=after && removed<=before && changed<=Math.min(before-removed, after-added) && before-removed==after-added, "Impossible profile membership delta: "+r[0]);
		final String expected=after==0 ? "EMPTY_REQUIRES_DECISION" : added==0 && removed==0 && changed==0 ? "REUSED" : "REBUILT";
		require(r[10].equals(expected), "Profile state contradicts membership delta: "+r[0]);
		if(expected.equals("REBUILT")){
			require(Paths.get(r[11]).normalize().equals(shard.resolve("rank_"+r[0])), "Rebuilt profile is outside its producing shard: "+r[0]);
		}else{
			require(r[11].equals(seed[2]) && r[12].equals(seed[3]) && r[13].equals(seed[4]), "Reused/empty seed pointer or hashes differ: "+r[0]);
		}
	}

	/** Revalidates the rebuilt singleton's actual provenance, fixed recipe and summary. */
	private static void validateRebuild(Path dir, String[] r, String[] member, String[] seed, HashMap<String,String> o, int length) throws Exception{
		require(text(dir.resolve("PASS")).equals("HBM_REASSIGNED_FAMILY_PASS"), "Incomplete rebuilt singleton: "+r[0]);
		final List<String[]> summary=rows(dir.resolve("summary.tsv").toString(), SUMMARY_HEADER);
		require(summary.size()==2 && summary.get(1).length==13, "Invalid rebuilt summary"); final String[] s=summary.get(1);
		for(int i=0; i<6; i++){require(s[i].equals(member[i]), "Rebuilt summary differs from membership at field "+i);}
		require(number(s[6])>0 && number(s[7])==length && number(s[8])>0 && number(s[9])>=0 && number(s[10])<=number(member[5])
			&& s[11].equals(r[12]) && s[12].equals(r[13]), "Rebuilt length/count/hash summary differs");
		final HashMap<String,String> recipe=fields(dir.resolve("recipe.txt"), '=', true);
		require("experimental_reassigned_profile_v1".equals(recipe.get("format")), "Unsupported rebuilt recipe");
		for(String setting : new String[]{"beta=0.01", "clipmin=-4", "clipmax=11", "gap=4", "trimdepth=0.1", "passes=2", "final_padding=0",
			"member_order=ascii_query_id", "final_terminal_policy=count_exclusions", "background_policy=frozen_original", "source_rank="+r[0], "family_id="+r[1],
			"source_rep_id="+r[2], "growth_padding="+s[8], "parent_hbm_sha80="+seed[3], "old_input_sha80="+seed[5], "new_input_sha80="+member[7], "background_sha80="+o.get("backgroundsha80")}){
			final int at=setting.indexOf('='); require(setting.substring(at+1).equals(recipe.get(setting.substring(0, at))), "Rebuilt recipe mismatch: "+setting);
		}
		final String producer=HmmComparisonData.required(o, "producer");
		for(String name : new String[]{"HbmProfileRefiner", "HbmPositionModel", "HbmMembershipCompare"}){
			HbmProfilePilot.requireHash(producer+"/current/prot/"+name+".java", recipe.get(name+"_source_sha80"));
			HbmProfilePilot.requireHash(producer+"/current/prot/"+name+".class", recipe.get(name+"_class_sha80"));
		}
		HbmProfilePilot.requireHash(member[6], member[7]);
		HbmProfilePilot.verifyProvenance(dir, producer, member[6], member[6], o.get("new"), "HbmReassignedProfiles");
	}

	/** Emits one coherent roster, consumer consensus and native bundle without graph edits. */
	private static void writeOutputs(Path out, List<String[]> members, String[][] registry, int rebuilt, HbmBundleLoader.Loaded combined,
			List<String> roster, HashMap<String,byte[]> sequences, HashMap<String,String> o) throws Exception{
		final ByteStreamWriter fasta=HmmComparisonData.writer(out.resolve("consensus.faa").toString()), order=HmmComparisonData.writer(out.resolve("roster.tsv").toString());
		final ByteStreamWriter identity=HmmComparisonData.writer(out.resolve("family_identity.tsv").toString()), artifacts=HmmComparisonData.writer(out.resolve("family_artifacts.tsv").toString());
		final ByteStreamWriter sources=HmmComparisonData.writer(out.resolve("members.tsv").toString());
		order.println("rank\trep_id\tmembers"); identity.println("active_index\tfamily_id\tsource_active_index\trep_id\tmembers\traw_residues"); artifacts.println(MODEL_HEADER); sources.println(OLD_HEADER);
		long total=0, raw=0, encoded=0; final ByteBuilder line=new ByteBuilder();
		for(int i=0; i<registry.length; i++){
			final String[] b=members.get(i+1), r=registry[i];
			HmmComparisonData.writeFasta(fasta, b[2]+" active_index="+i+" family_id="+b[1], sequences.get(b[2]));
			order.println(line.clear().append(i).tab().append(b[2]).tab().append(b[3]));
			identity.println(line.clear().append(i).tab().append(b[1]).tab().append(i).tab().append(b[2]).tab().append(b[3]).tab().append(b[4]));
			artifacts.println(line.clear().append(i).tab().append(b[2]).tab().append(r[11]).tab().append(r[12]).tab().append(r[13]).tab().append(b[7]));
			sources.println(line.clear().append(i).tab().append(b[1]).tab().append(b[2]).tab().append(b[3]).tab().append(b[4]).tab().append(b[6]).tab().append(b[7]));
			total=Math.addExact(total, number(b[3])); raw=Math.addExact(raw, number(b[4])); encoded=Math.addExact(encoded, number(b[5]));
		}
		HmmComparisonData.close(fasta); HmmComparisonData.close(order); HmmComparisonData.close(identity); HmmComparisonData.close(artifacts); HmmComparisonData.close(sources);
		Files.copy(Paths.get(o.get("background")), out.resolve("background.tsv")); HbmProfilePilot.requireHash(out.resolve("background.tsv").toString(), o.get("backgroundsha80"));
		final ByteStreamWriter recipe=HmmComparisonData.writer(out.resolve("recipe.txt").toString());
		recipe.println("experimental_reassigned_union_v1\noperation=compose_verified_singletons\nbackground_policy=frozen_original\nempty_policy=reject_unresolved\nbeta=0.01\nclipmin=-4\nclipmax=11\ngap=4");
		for(String key : INPUTS){recipe.println(key+"_sha80="+o.get(key+"sha80"));} HmmComparisonData.close(recipe);
		final byte[][] provenance=HbmProfilePilot.writeProvenance(out, HmmComparisonData.required(o, "runtime"), out.resolve("family_artifacts.tsv").toString(),
			out.resolve("family_identity.tsv").toString(), o.get("new"), "HbmReassignedUnion");
		combined.writeBundle(out.resolve("refined.mqhb"), provenance);
		final HbmBundleLoader.Loaded roundtrip=HbmBundleLoader.load(out.resolve("refined.mqhb"), roster, id->sequences.get(id), HbmBundleLoader.loadSemanticProvenance(out.resolve("provenance.tsv").toString()));
		roundtrip.assertStructuralEquivalent(combined);
		require(Arrays.equals(FamilyShortlistSidecarBuilder.loadConsensusFamilyOrder(out.resolve("consensus.faa").toString()), roster.toArray(new String[0])), "Union consumer family order differs");
		final ByteStreamWriter summary=HmmComparisonData.writer(out.resolve("summary.tsv").toString());
		summary.println("families\trebuilt\treused\tmembers\tderived_raw_residues\tencoded_residues\tbundle_sha80\tconsensus_sha80");
		summary.println(line.clear().append(registry.length).tab().append(rebuilt).tab().append(registry.length-rebuilt).tab().append(total).tab().append(raw).tab().append(encoded)
			.tab().append(DigestSuffix.file(out.resolve("refined.mqhb").toString())).tab().append(DigestSuffix.file(out.resolve("consensus.faa").toString()))); HmmComparisonData.close(summary);
	}

	private static ProteinSequence reference(String file){
		final Streamer input=StreamerFactory.makeStreamer(FileFormat.testInput(file, FileFormat.FASTA, null, true, false), null, true, -1);
		ProteinSequence result=null; input.start();
		try{
			for(ListNum<Read> list=input.nextList(); list!=null && list.size()>0; list=input.nextList()){
				for(Read read : list){require(result==null && read.mate==null, "Expected one unpaired consensus"); result=new ProteinSequence(HbmCompetitiveAssign.queryId(read.id), read.bases);}
			}
		}finally{input.close(); require(!input.errorState(), "Consensus input failed");}
		require(result!=null, "Empty consensus"); return result;
	}
	private static List<String[]> rows(String file, String header){
		final List<String[]> result=HmmComparisonData.rows(file); require(!result.isEmpty() && String.join("\t", result.get(0)).equals(header), "Unsupported table header: "+file); return result;
	}
	private static HashMap<String,String> fields(Path file, char separator, boolean format) throws Exception{
		final List<String> lines=Files.readAllLines(file, StandardCharsets.UTF_8); final HashMap<String,String> result=new HashMap<String,String>();
		require(!lines.isEmpty(), "Empty recipe: "+file); if(format){result.put("format", lines.get(0));}
		for(int i=format ? 1 : 0; i<lines.size(); i++){
			final String line=lines.get(i); final int at=line.indexOf(separator);
			require(at>0 && at<line.length()-1 && result.put(line.substring(0, at), line.substring(at+1))==null, "Invalid or duplicate recipe key: "+file);
		}return result;
	}
	private static String text(Path file) throws Exception{return new String(Files.readAllBytes(file), StandardCharsets.UTF_8).trim();}
	private static long number(String value){final long n=Long.parseLong(value); require(n>=0 && value.equals(Long.toString(n)), "Invalid count: "+value); return n;}
	private static void verifyInputs(HashMap<String,String> o){for(String key : INPUTS){HbmProfilePilot.requireHash(HmmComparisonData.required(o, key), HmmComparisonData.required(o, key+"sha80"));}}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final String[] INPUTS={"old", "new", "models", "shards", "background"};
	private static final String OLD_HEADER="active_index\tfamily_id\trep_id\tassigned_count\tsequence_bytes\tfasta_file\tfasta_sha80";
	private static final String NEW_HEADER="active_index\tfamily_id\trep_id\tassigned_count\tderived_raw_residues\tencoded_residues\tfasta_file\tfasta_sha80\tmembership_status";
	private static final String MODEL_HEADER="active_index\trep_id\tdirectory\tbundle_sha80\tconsensus_sha80\tinput_sha80";
	private static final String SHARD_HEADER="directory\tfamilies_sha80\trecipe_sha80";
	private static final String REGISTRY_HEADER="active_index\tfamily_id\trep_id\told_members\tnew_members\tderived_raw_residues\tencoded_residues\tadded_ids\tremoved_ids\tchanged_sequences\tstate\tdirectory\tbundle_sha80\tconsensus_sha80\tnew_input_sha80";
	private static final String SUMMARY_HEADER="active_index\tfamily_id\trep_id\tmembers\tderived_raw_residues\tencoded_residues\told_length\tnew_length\tpadding\tpadding_retries\tfinal_terminal_excluded\tbundle_sha80\tconsensus_sha80";
}
