package prot;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import structures.ByteBuilder;

/**
 * Builds the FROZEN candidate-A/B/C artifact-role manifest (plans/
 * CANDIDATE_A_ARTIFACT_ROLE_MANIFEST_DESIGN_v1.md, v3.1): one row per required BUILD -- a
 * fold-local rebuild for a selected focal family's own scoring role, or a fixed all-member decoy
 * artifact for a family used purely as a neighbor in someone else's challenge set -- computed as
 * pure set logic over four already-existing artifacts. Deliberately does NOT run the max-min
 * cover or any convergence loop (those are separately authorized, later steps); this tool only
 * decides WHICH artifacts must exist and WHAT training pool each one gets, before any ≤300
 * selection ever happens.
 * <p>
 * <b>Inputs:</b> {@code sample_manifest.tsv} ({@link SampleManifestGenerator}'s full 4432-family
 * output -- universe, selected focals, and every family's {@code cluster_path}/{@code
 * sidecar_path}/{@code group_partition_sha256}), {@code fold_manifest.tsv} ({@link
 * FoldManifestAggregator}'s output -- which groups are held out per fold, per focal), {@code
 * family_neighbors.tsv} ({@link FamilyNeighborManifest}'s output -- each focal's K=20 neighbor
 * families), and {@code family_members_v1.tsv} ({@link FamilyArtifactManifest}'s corpus-wide flat
 * membership table -- the all-member source for decoy rows).
 * <p>
 * <b>Build set:</b> let F be the selected focal families and N(x) focal x's neighbor list. Every
 * x&isin;F needs one FOCAL row per fold (excluding that fold's held-out groups, expanded from
 * x's own {@code cluster.tsv}). Every family y appearing in ANY N(x) needs exactly one DECOY row
 * (its full membership, from {@code family_members_v1.tsv}) -- including when y is ALSO a focal
 * family itself: being focal for its own scoring role never exempts a family from also needing
 * the plain all-member artifact used when it plays a decoy role elsewhere (the design doc's
 * "easy-to-miss case").
 * <p>
 * <b>Provenance binding (v3, Elly's item 1; UMP45: "ESSENTIAL... exactly what makes the role
 * manifest gatable"):</b> the output header records the schema version and a STREAMING SHA-256
 * (never {@link java.nio.file.Files#readAllBytes}, since {@code family_members_v1.tsv} covers the
 * full 23M-protein corpus and is multi-gigabyte) of all four upstream inputs. Every row's {@code
 * training_pool_sha256}/{@code heldout_groups_sha256}/{@code heldout_member_sha256} is the SHA-256
 * of a CANONICAL list construction (v3.1, Elly's item 6): duplicate IDs rejected outright,
 * non-ASCII IDs rejected outright (round 2, Elly -- {@code String.compareTo} only matches true
 * UTF-8 byte-lexical order for pure ASCII), lexical sort, exactly one LF after every ID including
 * the last. Every row's {@code group_partition_sha256} is carried from {@code sample_manifest.tsv}
 * AND independently re-verified against the actual on-disk {@code cluster_path}/{@code
 * sidecar_path} pair via {@link IdentityGroupOutputVerifier#verify} (v3.1, Elly's item 1/5;
 * cached per family, round 2, so a dual-role family is never re-verified twice) -- a
 * group-derived cover anchor (still unresolved, see the design doc §4.1) can shift if the
 * underlying grouping changes even while the flat member set stays byte-identical, so the
 * grouping itself must be bound, not just trusted from a copied string.
 * <p>
 * <b>Bounded-memory {@code family_members_v1.tsv} read (v3.1, Elly's item 5; round 2 fix):</b> the
 * reader never retains the whole 23M-ID corpus, and -- round 2's real fix -- never retains every
 * NEEDED family's full member list simultaneously either. It restricts to the decoy build set's
 * needed {@code rep_id}s (known up front from the build-set computation) and exploits the SAME
 * rep-contiguity convention {@link ConsensusRepBuilder#readClusters} already relies on for
 * mmseqs-style cluster TSVs: only one family's raw member list is ever open at a time, FINALIZED
 * into a tiny {@code {count,canonicalHash}} summary and DISCARDED the instant a different
 * {@code rep_id} appears, so peak memory is bounded by the single largest needed family, never the
 * sum across all needed families. Contiguity is VERIFIED, not assumed -- a needed {@code rep_id}
 * that reappears after being closed throws immediately rather than silently re-merging.
 * <p>
 * <b>Neighbor-list frozen contract (round 2, Elly's item 2):</b> {@code family_neighbors.tsv} is
 * parsed with the SAME rigor {@link NegativeQueryManifestGenerator#loadNeighbors} already applies
 * -- a single {@code #k\t20} header (duplicate header rejected), declared K required to equal the
 * FROZEN 20, and every focal required to have exactly K rows with neighbor_rank 1..K each
 * appearing exactly once, non-self, no duplicate neighbor family. Without this, a malformed or
 * truncated neighbor file (e.g. 19 rows) would silently change the build set with no error.
 * <p>
 * <b>Decoy member-count cross-check (round 2, Elly's item 3):</b> a decoy family's
 * {@code family_members_v1.tsv} row count is compared against its {@code sample_manifest.tsv}
 * {@code member_count} column and the run fails loud on any mismatch -- the design doc calls this
 * out explicitly (never assume the two independently-sourced membership descriptions agree).
 *
 * <p>Usage: {@code java -ea prot.ArtifactRoleManifestGenerator sample=sample_manifest.tsv
 *        fold=fold_manifest.tsv neighbors=family_neighbors.tsv members=family_members_v1.tsv
 *        out=artifact_role_manifest.tsv decoyusage=decoy_usage.tsv}
 * <br>Self-test: {@code java -ea prot.ArtifactRoleManifestGenerator selftest}
 *
 * @author Eru
 */
public class ArtifactRoleManifestGenerator {

	/** The frozen v5.2 fold-count cap (see NegativeQueryManifestGenerator's identical constant). */
	static final int MAX_FOLDS=10;
	/** The frozen candidate-independent neighbor-challenge K (see NegativeQueryManifestGenerator's
	 * identical constant) -- not a tunable, ALWAYS 20. */
	static final int FROZEN_K=20;

	public static void main(String[] args) throws Exception{
		if(args.length==1 && args[0].equalsIgnoreCase("selftest")){selftest(); return;}
		String sampleFile=null, foldFile=null, neighborsFile=null, membersFile=null,
			outFile=null, decoyUsageFile=null;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase();
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("sample")){sampleFile=b;}
			else if(a.equals("fold") || a.equals("foldmanifest")){foldFile=b;}
			else if(a.equals("neighbors")){neighborsFile=b;}
			else if(a.equals("members")){membersFile=b;}
			else if(a.equals("out")){outFile=b;}
			else if(a.equals("decoyusage")){decoyUsageFile=b;}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}
		if(sampleFile==null || foldFile==null || neighborsFile==null || membersFile==null
				|| outFile==null || decoyUsageFile==null){
			throw new RuntimeException("Required: sample= fold= neighbors= members= out= decoyusage=");
		}
		final Result r=process(sampleFile, foldFile, neighborsFile, membersFile, outFile, decoyUsageFile);
		System.err.println(r.summary);
	}

	static class UniverseFam {
		String repId, clusterPath, sidecarPath, groupPartitionSha256;
		long memberCount;
	}
	static class Sel { String repId; long groupCount; }
	static class FoldGroup { String groupId; long groupSize; int foldId; }
	static class RoleRow {
		String artifactId, familyId, role, clusterPath, sidecarPath, groupPartitionSha256,
			trainingPoolSource, trainingPoolSha256;
		int foldId;//-1 for decoy
		long trainingPoolMemberCount;
		long heldoutGroupCount=-1, heldoutMemberCount=-1;//-1 (rendered "-") for decoy
		String heldoutGroupsSha256, heldoutMemberSha256;//null (rendered "-") for decoy
	}
	static class Result { String summary; }
	/** Bounded per-family summary -- never the raw member list (round 2 fix, see class javadoc). */
	static class MemberSummary { long count; String sha256; }
	static class MemberStreamResult {
		final HashMap<String,MemberSummary> summaries=new HashMap<String,MemberSummary>();
		/** Largest size any SINGLE family's raw accumulator ever reached -- exposed so a selftest
		 * can directly prove boundedness (peak == the largest single needed family, never the sum
		 * across all needed families). */
		int maxOpenListSize=0;
	}

	static Result process(String sampleFile, String foldFile, String neighborsFile,
			String membersFile, String outFile, String decoyUsageFile) throws Exception{
		final HashMap<String,UniverseFam> universe=loadUniverse(sampleFile);
		final ArrayList<Sel> focals=loadSelectedFocals(sampleFile);
		final HashSet<String> focalIds=new HashSet<String>();
		for(Sel s : focals){focalIds.add(s.repId);}

		final HashMap<String,ArrayList<FoldGroup>> foldGroupsByFocal=loadFoldManifest(foldFile);
		requireExactSet(focalIds, foldGroupsByFocal.keySet(), foldFile);

		final HashMap<String,ArrayList<String>> neighborsByFocal=loadNeighbors(neighborsFile);
		requireExactSet(focalIds, neighborsByFocal.keySet(), neighborsFile);
		for(java.util.Map.Entry<String,ArrayList<String>> e : neighborsByFocal.entrySet()){
			for(String neighbor : e.getValue()){
				if(!universe.containsKey(neighbor)){
					throw new RuntimeException("Focal family '"+e.getKey()+"' names neighbor '"+neighbor
						+"' which does not exist in the universe ("+sampleFile+").");
				}
			}
		}

		// Decoy usage: decoyFamily -> sorted list of focal families that use it as a neighbor.
		// A family already focal for itself is NOT exempt from also needing a decoy row if it
		// appears in a DIFFERENT focal's neighbor list (the design doc's "easy-to-miss case").
		final TreeMap<String,ArrayList<String>> decoyUsage=new TreeMap<String,ArrayList<String>>();
		for(Sel focal : focals){
			for(String neighbor : neighborsByFocal.get(focal.repId)){
				decoyUsage.computeIfAbsent(neighbor, k -> new ArrayList<String>()).add(focal.repId);
			}
		}
		for(ArrayList<String> users : decoyUsage.values()){java.util.Collections.sort(users);}

		// Bounded-memory pass over family_members_v1.tsv for exactly the needed decoy families --
		// never retains a raw member list per family, only a tiny {count,hash} summary (round 2 fix).
		final MemberStreamResult memberStream=streamNeededFamilyMembers(membersFile, decoyUsage.keySet());

		final ArrayList<RoleRow> rows=new ArrayList<RoleRow>();
		// Cache so a dual-role family's group_partition_sha256 is verified exactly ONCE, not once
		// per role (round 2, Elly's item 5).
		final HashSet<String> partitionVerified=new HashSet<String>();

		for(Sel focal : focals){
			final UniverseFam fam=universe.get(focal.repId);
			verifyPartitionOnce(fam, sampleFile, partitionVerified);
			final ArrayList<FoldGroup> foldGroups=foldGroupsByFocal.get(focal.repId);
			if(foldGroups.size()!=focal.groupCount){
				throw new RuntimeException("Focal family '"+focal.repId+"': fold manifest has "
					+foldGroups.size()+" group rows, but the sample manifest declares group_count="
					+focal.groupCount+" -- the two files have diverged.");
			}
			int foldCount=0;
			for(FoldGroup g : foldGroups){foldCount=Math.max(foldCount, g.foldId+1);}
			final long expectedFoldCount=Math.min(focal.groupCount, MAX_FOLDS);
			if(foldCount!=expectedFoldCount){
				throw new RuntimeException("Focal family '"+focal.repId+"': observed foldCount="
					+foldCount+" but the frozen rule requires min(group_count="+focal.groupCount
					+", "+MAX_FOLDS+")="+expectedFoldCount+".");
			}
			final HashSet<Integer> distinctFoldIds=new HashSet<Integer>();
			for(FoldGroup g : foldGroups){distinctFoldIds.add(g.foldId);}
			for(int f=0; f<foldCount; f++){
				if(!distinctFoldIds.contains(f)){
					throw new RuntimeException("Focal family '"+focal.repId+"': fold IDs are not "
						+"contiguous from 0 -- missing fold "+f+" (observed folds: "+distinctFoldIds+").");
				}
			}

			final HashMap<Integer,ArrayList<String>> heldOutMembersByFold=expandGroups(
				fam.clusterPath, foldGroups);
			final ArrayList<String> allMembers=new ArrayList<String>();
			for(ArrayList<String> members : heldOutMembersByFold.values()){allMembers.addAll(members);}

			for(int f=0; f<foldCount; f++){
				final ArrayList<String> heldOut=heldOutMembersByFold.get(f);
				final ArrayList<String> heldOutGroupIds=new ArrayList<String>();
				for(FoldGroup g : foldGroups){if(g.foldId==f){heldOutGroupIds.add(g.groupId);}}
				final HashSet<String> heldOutSet=new HashSet<String>(heldOut);
				final ArrayList<String> training=new ArrayList<String>();
				for(String m : allMembers){if(!heldOutSet.contains(m)){training.add(m);}}

				final RoleRow row=new RoleRow();
				row.artifactId=focal.repId+"__focal_f"+f;
				row.familyId=focal.repId;
				row.role="focal_fold_local";
				row.foldId=f;
				row.clusterPath=fam.clusterPath;
				row.sidecarPath=fam.sidecarPath;
				row.groupPartitionSha256=fam.groupPartitionSha256;
				row.trainingPoolSource="cluster_path";
				row.trainingPoolMemberCount=training.size();
				row.trainingPoolSha256=canonicalListSha256(training);
				row.heldoutGroupCount=heldOutGroupIds.size();
				row.heldoutMemberCount=heldOut.size();
				row.heldoutGroupsSha256=canonicalListSha256(heldOutGroupIds);
				row.heldoutMemberSha256=canonicalListSha256(heldOut);
				rows.add(row);
			}
		}

		for(String decoyFamily : decoyUsage.keySet()){
			final UniverseFam fam=universe.get(decoyFamily);
			verifyPartitionOnce(fam, sampleFile, partitionVerified);
			final MemberSummary ms=memberStream.summaries.get(decoyFamily);
			// Round 2, Elly's item 3: the design doc explicitly requires this cross-check, never
			// assumed -- family_members_v1.tsv and sample_manifest.tsv are two independently
			// produced descriptions of "this family's members."
			if(ms.count!=fam.memberCount){
				throw new RuntimeException("Decoy family '"+decoyFamily+"': "+membersFile+" has "
					+ms.count+" member row(s) but "+sampleFile+" declares member_count="
					+fam.memberCount+" -- the two files have diverged.");
			}
			final RoleRow row=new RoleRow();
			row.artifactId=decoyFamily+"__decoy";
			row.familyId=decoyFamily;
			row.role="all_member_decoy";
			row.foldId=-1;
			row.clusterPath=fam.clusterPath;
			row.sidecarPath=fam.sidecarPath;
			row.groupPartitionSha256=fam.groupPartitionSha256;
			row.trainingPoolSource="family_members_v1";
			row.trainingPoolMemberCount=ms.count;
			row.trainingPoolSha256=ms.sha256;
			rows.add(row);
		}

		// Frozen sort order (this tool's own, no prior scope text to match): family ID lexical,
		// then fold_id numeric ascending -- decoy's foldId=-1 sorts before any real fold of the
		// SAME family, and a family with both a decoy row and focal rows lands in one contiguous
		// block (decoy first, then folds 0..foldCount-1).
		rows.sort((a, b) -> {
			final int c=a.familyId.compareTo(b.familyId);
			return (c!=0) ? c : Integer.compare(a.foldId, b.foldId);
		});

		final ByteBuilder out=new ByteBuilder(1<<20);
		out.append("#schema_version\t1\n");
		out.append("#sample_manifest_path\t").append(sampleFile).append("\t#sample_manifest_sha256\t")
			.append(streamingSha256Hex(sampleFile)).append('\n');
		out.append("#fold_manifest_path\t").append(foldFile).append("\t#fold_manifest_sha256\t")
			.append(streamingSha256Hex(foldFile)).append('\n');
		out.append("#family_neighbors_path\t").append(neighborsFile).append("\t#family_neighbors_sha256\t")
			.append(streamingSha256Hex(neighborsFile)).append('\n');
		out.append("#family_members_v1_path\t").append(membersFile).append("\t#family_members_v1_sha256\t")
			.append(streamingSha256Hex(membersFile)).append('\n');
		out.append("#artifact_id\tfamily_id\trole\tfold_id\tcluster_path\tsidecar_path\t"
			+"group_partition_sha256\ttraining_pool_source\ttraining_pool_member_count\t"
			+"training_pool_sha256\theldout_group_count\theldout_member_count\theldout_groups_sha256\t"
			+"heldout_member_sha256\n");
		for(RoleRow r : rows){
			out.append(r.artifactId).append('\t').append(r.familyId).append('\t').append(r.role)
				.append('\t').append(r.foldId).append('\t').append(r.clusterPath).append('\t')
				.append(r.sidecarPath).append('\t').append(r.groupPartitionSha256).append('\t')
				.append(r.trainingPoolSource).append('\t').append(r.trainingPoolMemberCount).append('\t')
				.append(r.trainingPoolSha256).append('\t')
				.append(r.heldoutGroupCount<0 ? "-" : String.valueOf(r.heldoutGroupCount)).append('\t')
				.append(r.heldoutMemberCount<0 ? "-" : String.valueOf(r.heldoutMemberCount)).append('\t')
				.append(r.heldoutGroupsSha256==null ? "-" : r.heldoutGroupsSha256).append('\t')
				.append(r.heldoutMemberSha256==null ? "-" : r.heldoutMemberSha256).append('\n');
		}
		writeFile(outFile, out.toString());
		writeSha256Sidecar(outFile);

		final ByteBuilder du=new ByteBuilder(1<<16);
		du.append("#decoy_family_id\tfocal_family_id\n");
		for(java.util.Map.Entry<String,ArrayList<String>> e : decoyUsage.entrySet()){
			for(String focalId : e.getValue()){du.append(e.getKey()).append('\t').append(focalId).append('\n');}
		}
		writeFile(decoyUsageFile, du.toString());
		writeSha256Sidecar(decoyUsageFile);

		final Result res=new Result();
		res.summary="ArtifactRoleManifestGenerator: "+focals.size()+" focal families, "+decoyUsage.size()
			+" distinct decoy families, "+rows.size()+" total build rows. Manifest written to "+outFile
			+"; decoy usage written to "+decoyUsageFile+".";
		return res;
	}

	/** Writes a {@code <path>.sha256} sidecar (Elly's strong recommendation), hashed from bytes
	 * READ BACK FROM DISK after the writer closes -- matching FoldManifestAggregator's and
	 * NegativeQueryManifestGenerator's established convention. Both output files here are small
	 * (bounded by the sample size, never corpus-wide), so a one-shot readAllBytes is safe -- unlike
	 * the multi-GB family_members_v1.tsv input, which uses streamingSha256Hex instead. */
	static void writeSha256Sidecar(String path) throws IOException{
		final byte[] bytes=java.nio.file.Files.readAllBytes(java.nio.file.Paths.get(path));
		writeFile(path+".sha256", sha256Hex(bytes)+"  "+path+"\n");
	}

	/** Verifies a family's claimed group_partition_sha256 exactly once per distinct family
	 * (round 2, Elly's item 5 -- a dual-role focal+decoy family must not be re-verified twice). */
	static void verifyPartitionOnce(UniverseFam fam, String sampleFile, HashSet<String> verified){
		if(!verified.add(fam.repId)){return;}
		verifyPartition(fam, sampleFile);
	}

	/** Independently re-verifies a family's claimed group_partition_sha256 against the actual
	 * on-disk cluster_path/sidecar_path pair (v3.1, Elly's item 1/5) -- never trusts the copied
	 * sample_manifest.tsv string blindly. */
	static void verifyPartition(UniverseFam fam, String sampleFile){
		if(fam.memberCount>Integer.MAX_VALUE){
			throw new RuntimeException("Family '"+fam.repId+"' has member_count="+fam.memberCount
				+" which exceeds Integer.MAX_VALUE -- IdentityGroupOutputVerifier.verify() cannot "
				+"accept it (round 2, Elly's item 5 -- guard the cast explicitly).");
		}
		final IdentityGroupOutputVerifier.Result vr=IdentityGroupOutputVerifier.verify(
			fam.clusterPath, fam.sidecarPath, (int)fam.memberCount);
		if(!vr.pass){
			throw new RuntimeException("Family '"+fam.repId+"': IdentityGroupOutputVerifier found "
				+vr.failures.size()+" failure(s) re-verifying its cluster.tsv/sidecar from primary "
				+"bytes: "+vr.failures.get(0));
		}
		if(!vr.rawHash.equals(fam.groupPartitionSha256)){
			throw new RuntimeException("Family '"+fam.repId+"': "+sampleFile+" claims "
				+"group_partition_sha256='"+fam.groupPartitionSha256+"' but the independently "
				+"recomputed hash of its cluster.tsv is '"+vr.rawHash+"' -- the sample manifest and "
				+"the actual grouping on disk have diverged.");
		}
	}

	/** Expands each fold's held-out group_id(s) into their real member lists via the focal
	 * family's own cluster.tsv, mirroring NegativeQueryManifestGenerator.expandHeldOutGroups: a
	 * member ID appearing twice (same or a different group) is rejected globally, and every
	 * cluster.tsv group must be covered by exactly one fold (both directions checked). */
	static HashMap<Integer,ArrayList<String>> expandGroups(String clusterPath, ArrayList<FoldGroup> foldGroups){
		final HashMap<String,Integer> foldByGroupId=new HashMap<String,Integer>();
		for(FoldGroup g : foldGroups){
			if(foldByGroupId.put(g.groupId, g.foldId)!=null){
				throw new RuntimeException("Duplicate group_id '"+g.groupId+"' in fold manifest for "
					+"this focal family.");
			}
		}
		final HashMap<Integer,ArrayList<String>> byFold=new HashMap<Integer,ArrayList<String>>();
		final HashMap<String,Integer> seenGroupMemberCount=new HashMap<String,Integer>();
		final HashSet<String> clusterGroupIds=new HashSet<String>();
		final HashSet<String> seenMemberIds=new HashSet<String>();
		final ByteFile bf=ByteFile.makeByteFile(clusterPath, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			lp.set(line);
			if(lp.terms()<2){throw new RuntimeException("Malformed cluster.tsv row: "+new String(line));}
			final String groupId=lp.parseString(0), memberId=lp.parseString(1);
			if(!seenMemberIds.add(memberId)){
				throw new RuntimeException("cluster.tsv "+clusterPath+" has member '"+memberId
					+"' appearing more than once (same or a different group) -- membership rows must "
					+"be unique.");
			}
			clusterGroupIds.add(groupId);
			final Integer foldId=foldByGroupId.get(groupId);
			if(foldId==null){continue;}
			byFold.computeIfAbsent(foldId, k -> new ArrayList<String>()).add(memberId);
			seenGroupMemberCount.merge(groupId, 1, Integer::sum);
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+clusterPath);}
		for(FoldGroup g : foldGroups){
			final Integer seen=seenGroupMemberCount.get(g.groupId);
			if(seen==null || seen!=g.groupSize){
				throw new RuntimeException("Group '"+g.groupId+"': fold manifest declares size "
					+g.groupSize+" but cluster.tsv "+clusterPath+" has "+(seen==null ? 0 : seen)
					+" member rows for it -- fold manifest and cluster.tsv have diverged.");
			}
		}
		final ArrayList<String> uncovered=new ArrayList<String>();
		for(String groupId : clusterGroupIds){
			if(!foldByGroupId.containsKey(groupId)){uncovered.add(groupId);}
		}
		if(!uncovered.isEmpty()){
			throw new RuntimeException("cluster.tsv "+clusterPath+" has "+uncovered.size()
				+" group(s) not covered by any fold in the fold manifest (e.g. '"+uncovered.get(0)
				+"') -- fold manifest and cluster.tsv have diverged.");
		}
		return byFold;
	}

	/** Bounded-memory pass over family_members_v1.tsv for exactly the needed rep_ids (v3.1, Elly's
	 * item 5; round 2's real fix): never retains the whole 23M-ID corpus, and never retains every
	 * needed family's full member list at once either -- ONLY one raw list is ever alive, finalized
	 * into a tiny {count,hash} summary and discarded the instant its family's block closes.
	 * Exploits real-file rep-contiguity (the same convention ConsensusRepBuilder.readClusters
	 * already relies on). Contiguity is VERIFIED, not assumed: a needed rep_id reappearing after
	 * being closed (a different rep_id intervened) throws immediately rather than silently
	 * re-merging. */
	static MemberStreamResult streamNeededFamilyMembers(String membersFile, java.util.Set<String> neededIds){
		final MemberStreamResult result=new MemberStreamResult();
		final HashSet<String> closed=new HashSet<String>();
		String currentRep=null;
		ArrayList<String> currentList=null;
		final ByteFile bf=ByteFile.makeByteFile(membersFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			lp.set(line);
			if(lp.terms()<2){throw new RuntimeException("Malformed members row: "+new String(line));}
			final String repId=lp.parseString(0), memberId=lp.parseString(1);
			if(repId.equals(currentRep)){
				if(currentList!=null){
					currentList.add(memberId);
					result.maxOpenListSize=Math.max(result.maxOpenListSize, currentList.size());
				}
				continue;
			}
			finalizeFamily(currentRep, currentList, neededIds, closed, result.summaries);
			currentRep=repId;
			if(neededIds.contains(repId)){
				if(closed.contains(repId)){
					throw new RuntimeException("family_members_v1.tsv "+membersFile+": needed family '"
						+repId+"' has rows that are NOT contiguous (it reappeared after being closed by "
						+"a different family's rows in between) -- membership rows for one family must "
						+"be grouped together.");
				}
				currentList=new ArrayList<String>();
				currentList.add(memberId);
				result.maxOpenListSize=Math.max(result.maxOpenListSize, currentList.size());
			}else{
				currentList=null;
			}
		}
		finalizeFamily(currentRep, currentList, neededIds, closed, result.summaries);
		if(bf.close()){throw new RuntimeException("I/O error reading "+membersFile);}
		for(String id : neededIds){
			if(!result.summaries.containsKey(id)){
				throw new RuntimeException("Needed decoy family '"+id+"' has no rows in "+membersFile+".");
			}
		}
		return result;
	}

	/** Finalizes ONE family's raw accumulator into a tiny {count,hash} summary and discards the
	 * list (the round-2 memory fix: the caller never retains more than one raw list at a time).
	 * No-op if rep is null, not needed, or already closed. */
	static void finalizeFamily(String rep, ArrayList<String> list, java.util.Set<String> neededIds,
			HashSet<String> closed, HashMap<String,MemberSummary> summaries){
		if(rep==null || !neededIds.contains(rep) || closed.contains(rep)){return;}
		final MemberSummary s=new MemberSummary();
		s.count=list.size();
		s.sha256=canonicalListSha256(list);
		summaries.put(rep, s);
		closed.add(rep);
	}

	/** Requires the two given ID sets be EXACTLY equal, both directions (matches
	 * NegativeQueryManifestGenerator.requireExactFocalSet's convention). */
	static void requireExactSet(HashSet<String> expected, java.util.Set<String> actual, String fileLabel){
		final ArrayList<String> missing=new ArrayList<String>(), extra=new ArrayList<String>();
		for(String id : expected){if(!actual.contains(id)){missing.add(id);}}
		for(String id : actual){if(!expected.contains(id)){extra.add(id);}}
		if(!missing.isEmpty() || !extra.isEmpty()){
			throw new RuntimeException(fileLabel+"'s focal-family set does not exactly match the "
				+"sample manifest's selected focals -- missing="+missing.size()+" ("
				+(missing.isEmpty() ? "" : missing.get(0))+"), extra="+extra.size()+" ("
				+(extra.isEmpty() ? "" : extra.get(0))+").");
		}
	}

	/** Canonical list-hash construction (v3.1, Elly's item 6, pinned exactly; round-2 correction):
	 * reject any duplicate ID (throw, never dedupe); reject any NON-ASCII ID (throw -- {@code
	 * String.compareTo} only matches true UTF-8 byte-lexical order for pure-ASCII strings, so a
	 * non-ASCII ID would silently break the "UTF-8 lexical order" claim rather than merely sort
	 * unexpectedly); sort by lexical order; join with exactly one LF after EVERY id including the
	 * last; SHA-256 over that exact byte sequence. Bounded per-family (never corpus-wide), so an
	 * in-memory build is safe here. */
	static String canonicalListSha256(ArrayList<String> ids){
		final ArrayList<String> sorted=new ArrayList<String>(ids);
		java.util.Collections.sort(sorted);
		final HashSet<String> seen=new HashSet<String>();
		final ByteBuilder bb=new ByteBuilder();
		for(String id : sorted){
			for(int i=0; i<id.length(); i++){
				if(id.charAt(i)>127){
					throw new RuntimeException("ID '"+id+"' contains a non-ASCII character -- "
						+"canonical list hashing requires pure-ASCII IDs (String.compareTo only "
						+"matches true UTF-8 byte-lexical order for ASCII).");
				}
			}
			if(!seen.add(id)){
				throw new RuntimeException("Duplicate ID '"+id+"' in a canonical list -- membership "
					+"must be unique.");
			}
			bb.append(id).append('\n');
		}
		return sha256Hex(bb.toBytes());
	}

	/** The FULL sample_manifest.tsv, every family, parsed by HEADER NAME. */
	static HashMap<String,UniverseFam> loadUniverse(String sampleFile){
		final HashMap<String,UniverseFam> map=new HashMap<String,UniverseFam>();
		final ByteFile bf=ByteFile.makeByteFile(sampleFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		HashMap<String,Integer> col=null;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='#'){
				if(col==null){
					final String h=new String(line, 1, line.length-1);
					col=new HashMap<String,Integer>();
					final String[] names=h.split("\t");
					for(int i=0; i<names.length; i++){col.put(names[i].trim().toLowerCase(), i);}
				}
				continue;
			}
			if(col==null){throw new RuntimeException("Data row before header in "+sampleFile);}
			lp.set(line);
			final UniverseFam f=new UniverseFam();
			f.repId=lp.parseString(idx(col, "rep_id", "family_id"));
			f.memberCount=lp.parseLong(idx(col, "member_count"));
			f.clusterPath=lp.parseString(idx(col, "cluster_path"));
			f.sidecarPath=lp.parseString(idx(col, "sidecar_path"));
			f.groupPartitionSha256=lp.parseString(idx(col, "group_partition_sha256"));
			if(map.put(f.repId, f)!=null){
				throw new RuntimeException("Duplicate rep_id '"+f.repId+"' in "+sampleFile);
			}
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+sampleFile);}
		return map;
	}

	/** sample_manifest.tsv's selected=1 rows only, taking rep_id and group_count. */
	static ArrayList<Sel> loadSelectedFocals(String sampleFile){
		final ArrayList<Sel> selected=new ArrayList<Sel>();
		final HashSet<String> seen=new HashSet<String>();
		final ByteFile bf=ByteFile.makeByteFile(sampleFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		HashMap<String,Integer> col=null;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='#'){
				if(col==null){
					final String h=new String(line, 1, line.length-1);
					col=new HashMap<String,Integer>();
					final String[] names=h.split("\t");
					for(int i=0; i<names.length; i++){col.put(names[i].trim().toLowerCase(), i);}
				}
				continue;
			}
			if(col==null){throw new RuntimeException("Data row before header in "+sampleFile);}
			lp.set(line);
			final int selectedFlag=lp.parseInt(idx(col, "selected"));
			if(selectedFlag!=0 && selectedFlag!=1){
				throw new RuntimeException("Malformed 'selected' value "+selectedFlag+" in "+sampleFile
					+" (must be exactly 0 or 1): "+new String(line));
			}
			if(selectedFlag!=1){continue;}
			final String repId=lp.parseString(idx(col, "rep_id", "family_id"));
			final long groupCount=lp.parseLong(idx(col, "group_count"));
			if(groupCount<2){
				throw new RuntimeException("Selected focal '"+repId+"' has group_count="+groupCount
					+" (<2) -- violates the frozen eligibility floor.");
			}
			if(!seen.add(repId)){
				throw new RuntimeException("Duplicate SELECTED rep_id '"+repId+"' in "+sampleFile+".");
			}
			final Sel s=new Sel();
			s.repId=repId; s.groupCount=groupCount;
			selected.add(s);
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+sampleFile);}
		return selected;
	}

	/** FoldManifestAggregator's output: #family_id\tgroup_id\tgroup_size\tfold_id\tgroup_hash. */
	static HashMap<String,ArrayList<FoldGroup>> loadFoldManifest(String foldManifestFile){
		final HashMap<String,ArrayList<FoldGroup>> byFocal=new HashMap<String,ArrayList<FoldGroup>>();
		final ByteFile bf=ByteFile.makeByteFile(foldManifestFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			lp.set(line);
			if(lp.terms()<5){throw new RuntimeException("Malformed fold-manifest row: "+new String(line));}
			final String familyId=lp.parseString(0);
			final FoldGroup g=new FoldGroup();
			g.groupId=lp.parseString(1);
			g.groupSize=lp.parseLong(2);
			g.foldId=lp.parseInt(3);
			final String groupHash=lp.parseString(4);
			if(g.groupSize<=0){throw new RuntimeException("Fold-manifest row for family '"+familyId
				+"' group '"+g.groupId+"' has nonpositive group_size "+g.groupSize+".");}
			if(g.foldId<0){throw new RuntimeException("Fold-manifest row for family '"+familyId
				+"' group '"+g.groupId+"' has negative fold_id "+g.foldId+".");}
			if(!groupHash.matches("[0-9a-f]{16}")){throw new RuntimeException("Fold-manifest row for "
				+"family '"+familyId+"' group '"+g.groupId+"' has a malformed group_hash '"+groupHash
				+"'.");}
			final String expectedHash=hex16(ReducedAlphabetSeedAssay.stableHash(g.groupId));
			if(!groupHash.equals(expectedHash)){
				throw new RuntimeException("Fold-manifest row for family '"+familyId+"' group '"
					+g.groupId+"' has group_hash '"+groupHash+"' but the real FNV-1a-64 of that "
					+"group_id is '"+expectedHash+"'.");
			}
			byFocal.computeIfAbsent(familyId, k -> new ArrayList<FoldGroup>()).add(g);
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+foldManifestFile);}
		return byFocal;
	}

	/** FamilyNeighborManifest's output: #focal_family\tneighbor_rank\tneighbor_family\t..., a
	 * single #k\t20 header. Round 2, Elly's item 2: enforces the frozen K=20 contract with the
	 * SAME rigor as NegativeQueryManifestGenerator.loadNeighbors -- a malformed/truncated file
	 * (e.g. 19 rows) previously changed the build set with no error at all. */
	static HashMap<String,ArrayList<String>> loadNeighbors(String neighborsFile){
		final HashMap<String,ArrayList<String>> byFocal=new HashMap<String,ArrayList<String>>();
		final HashMap<String,HashSet<Integer>> ranksSeenByFocal=new HashMap<String,HashSet<Integer>>();
		final HashMap<String,HashSet<String>> neighborsSeenByFocal=new HashMap<String,HashSet<String>>();
		int declaredK=-1;
		final ByteFile bf=ByteFile.makeByteFile(neighborsFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='#'){
				final String s=new String(line);
				if(s.startsWith("#k\t")){
					if(declaredK>=0){
						throw new RuntimeException(neighborsFile+" has more than one '#k\\t<K>' header "
							+"line -- ambiguous which is authoritative.");
					}
					declaredK=Integer.parseInt(s.substring(3).trim());
				}
				continue;
			}
			lp.set(line);
			if(lp.terms()<3){throw new RuntimeException("Malformed neighbor row: "+new String(line));}
			final String focal=lp.parseString(0);
			final int rank=lp.parseInt(1);
			final String neighbor=lp.parseString(2);
			if(neighbor.equals(focal)){
				throw new RuntimeException("Focal family '"+focal+"' lists itself as its own neighbor "
					+"(rank "+rank+") -- neighbors must be non-self.");
			}
			if(!ranksSeenByFocal.computeIfAbsent(focal, k -> new HashSet<Integer>()).add(rank)){
				throw new RuntimeException("Focal family '"+focal+"' has neighbor_rank "+rank
					+" recorded more than once.");
			}
			if(!neighborsSeenByFocal.computeIfAbsent(focal, k -> new HashSet<String>()).add(neighbor)){
				throw new RuntimeException("Focal family '"+focal+"' has neighbor family '"+neighbor
					+"' recorded more than once (duplicate neighbor).");
			}
			byFocal.computeIfAbsent(focal, k -> new ArrayList<String>()).add(neighbor);
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+neighborsFile);}
		if(declaredK<1){
			throw new RuntimeException(neighborsFile+" has no valid '#k\\t<K>' header line -- cannot "
				+"validate neighbor completeness without the declared K.");
		}
		if(declaredK!=FROZEN_K){
			throw new RuntimeException(neighborsFile+" declares K="+declaredK+", but this assay's K "
				+"is FROZEN at "+FROZEN_K+" -- refusing to proceed with a non-frozen neighbor count.");
		}
		for(java.util.Map.Entry<String,ArrayList<String>> e : byFocal.entrySet()){
			final String focal=e.getKey();
			if(e.getValue().size()!=declaredK){
				throw new RuntimeException("Focal family '"+focal+"' has "+e.getValue().size()
					+" neighbor rows, expected exactly the declared K="+declaredK+".");
			}
			final HashSet<Integer> ranks=ranksSeenByFocal.get(focal);
			for(int r=1; r<=declaredK; r++){
				if(!ranks.contains(r)){
					throw new RuntimeException("Focal family '"+focal+"' is missing neighbor_rank "+r
						+" (ranks must be exactly 1.."+declaredK+", each appearing once).");
				}
			}
		}
		return byFocal;
	}

	static int idx(HashMap<String,Integer> col, String... names){
		for(String n : names){final Integer i=col.get(n); if(i!=null){return i;}}
		throw new RuntimeException("Missing required column (any of): "+String.join(",", names));
	}

	static String hex16(long value){
		final char[] c=new char[16];
		for(int i=15, s=0; i>=0; i--, s+=4){final int d=(int)((value>>>s)&15); c[i]=(char)(d<10 ? '0'+d : 'a'+d-10);}
		return new String(c);
	}

	static String sha256Hex(byte[] data){
		try{
			final MessageDigest md=MessageDigest.getInstance("SHA-256");
			final byte[] digest=md.digest(data);
			final StringBuilder sb=new StringBuilder(digest.length*2);
			for(byte b : digest){sb.append(String.format("%02x", b));}
			return sb.toString();
		}catch(NoSuchAlgorithmException e){throw new RuntimeException(e);}
	}

	/** Streaming SHA-256 (v3, Elly + UMP45: "never Files.readAllBytes() for a multi-GB input" --
	 * family_members_v1.tsv covers the full 23M-protein corpus). Constant memory regardless of
	 * file size. */
	static String streamingSha256Hex(String path) throws IOException{
		final MessageDigest md;
		try{md=MessageDigest.getInstance("SHA-256");}
		catch(NoSuchAlgorithmException e){throw new RuntimeException(e);}
		try(BufferedInputStream in=new BufferedInputStream(new FileInputStream(path), 1<<20)){
			final byte[] buf=new byte[1<<20];
			int n;
			while((n=in.read(buf))>=0){md.update(buf, 0, n);}
		}
		final byte[] digest=md.digest();
		final StringBuilder sb=new StringBuilder(digest.length*2);
		for(byte b : digest){sb.append(String.format("%02x", b));}
		return sb.toString();
	}

	static void writeFile(String path, String content){
		final FileFormat ff=FileFormat.testOutput(path, FileFormat.TXT, null, false, true, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		bsw.print(new ByteBuilder().append(content));
		if(bsw.poisonAndWait()){throw new RuntimeException("I/O error writing "+path);}
	}

	static String readWhole(String path){
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		final StringBuilder sb=new StringBuilder();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){sb.append(new String(line)).append('\n');}
		bf.close();
		return sb.toString();
	}

	// ---------------- self-test ----------------
	// Synthetic-fixtures-only (production/canary data never touched). Real K=20 neighbor blocks
	// (the frozen contract, not a toy K=1) built from real IdentityGroupBuilder-produced
	// cluster.tsv/sidecar fixtures: focalA (2 groups -> 2 folds) and focalB (2 groups -> 2 folds)
	// are each other's... no -- focalB's neighbor list INCLUDES focalA (the dual-role case),
	// focalA's neighbor list includes decoyOnly (a plain single-role decoy), and both focals'
	// remaining 19 neighbor slots are filled by tiny single-member placeholder families.

	static void selftest() throws Exception{
		final String dir=System.getProperty("java.io.tmpdir")+"/armg_selftest_"+System.nanoTime();
		new java.io.File(dir).mkdirs();

		final ArrayList<String> allRepIds=new ArrayList<String>();
		final ArrayList<Boolean> allSelected=new ArrayList<Boolean>();
		final ArrayList<Long> allGroupCounts=new ArrayList<Long>();
		final ArrayList<Long> allMemberCounts=new ArrayList<Long>();
		final ByteBuilder membersTsv=new ByteBuilder();
		membersTsv.append("#rep_id\tmember_id\n");

		// focalA: 4 members, 2 groups {fa1,fa2} and {fa3,fa4} -> 2 folds.
		buildRealFamily(dir, "focalA", new String[]{"fa1", "fa2", "fa3", "fa4"},
			new String[][]{{"fa1", "fa2"}, {"fa3", "fa4"}});
		registerFamily(allRepIds, allSelected, allGroupCounts, allMemberCounts, membersTsv,
			"focalA", true, 2, new String[]{"fa1", "fa2", "fa3", "fa4"});

		// focalB: 4 members, 2 groups {fb1,fb2} and {fb3,fb4} -> 2 folds. focalB's neighbor list
		// includes focalA -- so focalA needs a decoy row IN ADDITION to its own 2 focal-fold rows.
		buildRealFamily(dir, "focalB", new String[]{"fb1", "fb2", "fb3", "fb4"},
			new String[][]{{"fb1", "fb2"}, {"fb3", "fb4"}});
		registerFamily(allRepIds, allSelected, allGroupCounts, allMemberCounts, membersTsv,
			"focalB", true, 2, new String[]{"fb1", "fb2", "fb3", "fb4"});

		// decoyOnly: never focal, only ever a neighbor -- a plain single-role decoy family.
		buildRealFamily(dir, "decoyOnly", new String[]{"d1", "d2"}, new String[][]{{"d1", "d2"}});
		registerFamily(allRepIds, allSelected, allGroupCounts, allMemberCounts, membersTsv,
			"decoyOnly", false, 1, new String[]{"d1", "d2"});

		// 19 tiny single-member placeholder families filling focalA's remaining neighbor slots,
		// and 19 more filling focalB's -- the frozen K=20 contract needs 20 real neighbors each.
		final String[] focalANeighbors=new String[FROZEN_K];
		focalANeighbors[0]="decoyOnly";
		for(int i=1; i<FROZEN_K; i++){
			final String id="nbrA"+(i-1);
			buildRealFamily(dir, id, new String[]{id+"_m"}, new String[][]{{id+"_m"}});
			registerFamily(allRepIds, allSelected, allGroupCounts, allMemberCounts, membersTsv,
				id, false, 1, new String[]{id+"_m"});
			focalANeighbors[i]=id;
		}
		final String[] focalBNeighbors=new String[FROZEN_K];
		focalBNeighbors[0]="focalA";
		for(int i=1; i<FROZEN_K; i++){
			final String id="nbrB"+(i-1);
			buildRealFamily(dir, id, new String[]{id+"_m"}, new String[][]{{id+"_m"}});
			registerFamily(allRepIds, allSelected, allGroupCounts, allMemberCounts, membersTsv,
				id, false, 1, new String[]{id+"_m"});
			focalBNeighbors[i]=id;
		}

		writeFile(dir+"/sample.tsv", sampleManifestRows(dir, allRepIds, allSelected, allGroupCounts,
			allMemberCounts));
		writeFile(dir+"/members.tsv", membersTsv.toString());

		writeFile(dir+"/fold_manifest.tsv", "#family_id\tgroup_id\tgroup_size\tfold_id\tgroup_hash\n"
			+"focalA\tfa1\t2\t0\t"+hex16(ReducedAlphabetSeedAssay.stableHash("fa1"))+"\n"
			+"focalA\tfa3\t2\t1\t"+hex16(ReducedAlphabetSeedAssay.stableHash("fa3"))+"\n"
			+"focalB\tfb1\t2\t0\t"+hex16(ReducedAlphabetSeedAssay.stableHash("fb1"))+"\n"
			+"focalB\tfb3\t2\t1\t"+hex16(ReducedAlphabetSeedAssay.stableHash("fb3"))+"\n");

		writeFile(dir+"/neighbors.tsv", "#k\t"+FROZEN_K+"\n#focal_family\tneighbor_rank\tneighbor_family\n"
			+neighborBlock("focalA", focalANeighbors) + neighborBlock("focalB", focalBNeighbors));

		final String out=dir+"/role.tsv", du=dir+"/decoy_usage.tsv";
		process(dir+"/sample.tsv", dir+"/fold_manifest.tsv", dir+"/neighbors.tsv", dir+"/members.tsv",
			out, du);
		final String outText=readWhole(out), duText=readWhole(du);

		// (1) focalA gets BOTH its 2 own focal-fold rows AND a decoy row (the dual-role case).
		int focalAFolds=0; boolean focalADecoy=false;
		for(String line : outText.split("\n")){
			if(line.startsWith("focalA__focal_f")){focalAFolds++;}
			if(line.startsWith("focalA__decoy\t")){focalADecoy=true;}
		}
		if(focalAFolds!=2){throw new RuntimeException("SELFTEST FAILED: expected 2 focalA fold rows, got "+focalAFolds);}
		if(!focalADecoy){throw new RuntimeException("SELFTEST FAILED: focalA (also a neighbor of focalB) "
			+"is missing its all-member decoy row -- the dual-role case is broken.");}
		System.err.println("  selftest[dual-role family gets BOTH focal-fold rows AND a decoy row]: PASS");

		// (2) decoyOnly gets exactly one decoy row and NO focal rows (never selected).
		int decoyOnlyRows=0; boolean decoyOnlyHasFocalRow=false;
		for(String line : outText.split("\n")){
			if(line.startsWith("decoyOnly__decoy\t")){decoyOnlyRows++;}
			if(line.startsWith("decoyOnly__focal_f")){decoyOnlyHasFocalRow=true;}
		}
		if(decoyOnlyRows!=1){throw new RuntimeException("SELFTEST FAILED: expected exactly 1 decoyOnly row, got "+decoyOnlyRows);}
		if(decoyOnlyHasFocalRow){throw new RuntimeException("SELFTEST FAILED: decoyOnly (never selected) has a focal row.");}
		System.err.println("  selftest[plain single-role decoy family]: PASS");

		// (3) focalB gets exactly 2 focal-fold rows and NO decoy row (never anyone's neighbor).
		int focalBFolds=0; boolean focalBDecoy=false;
		for(String line : outText.split("\n")){
			if(line.startsWith("focalB__focal_f")){focalBFolds++;}
			if(line.startsWith("focalB__decoy\t")){focalBDecoy=true;}
		}
		if(focalBFolds!=2){throw new RuntimeException("SELFTEST FAILED: expected 2 focalB fold rows, got "+focalBFolds);}
		if(focalBDecoy){throw new RuntimeException("SELFTEST FAILED: focalB (never anyone's neighbor) has a decoy row.");}
		System.err.println("  selftest[focal-only family, no decoy role needed]: PASS");

		// (4) decoy_usage.tsv records focalA used-by focalB, and decoyOnly used-by focalA.
		if(!duText.contains("focalA\tfocalB") || !duText.contains("decoyOnly\tfocalA")){
			throw new RuntimeException("SELFTEST FAILED: decoy_usage.tsv missing an expected "
				+"(decoy,focal) pair. Full text:\n"+duText);
		}
		System.err.println("  selftest[decoy_usage cross-reference]: PASS");

		// (5) heldout+training member sets partition focalA's 4 members exactly, per fold.
		for(String line : outText.split("\n")){
			if(!line.startsWith("focalA__focal_f")){continue;}
			final String[] cols=line.split("\t");
			final long trainingCount=Long.parseLong(cols[8]);
			final long heldoutCount=Long.parseLong(cols[11]);
			if(trainingCount+heldoutCount!=4){
				throw new RuntimeException("SELFTEST FAILED: focalA fold row's training("+trainingCount
					+")+heldout("+heldoutCount+") != 4 total members: "+line);
			}
		}
		System.err.println("  selftest[training+heldout partition every member exactly]: PASS");

		// (6) provenance header present with all 4 upstream hashes, and both output files' own
		// .sha256 sidecars exist and match an independent recomputation.
		if(!outText.contains("#schema_version\t1") || !outText.contains("#sample_manifest_sha256")
				|| !outText.contains("#fold_manifest_sha256") || !outText.contains("#family_neighbors_sha256")
				|| !outText.contains("#family_members_v1_sha256")){
			throw new RuntimeException("SELFTEST FAILED: provenance header incomplete.\n"+outText.substring(0, Math.min(400, outText.length())));
		}
		final String outSha=readWhole(out+".sha256").trim().split("\\s+")[0];
		final String duSha=readWhole(du+".sha256").trim().split("\\s+")[0];
		if(!outSha.equals(sha256Hex(java.nio.file.Files.readAllBytes(java.nio.file.Paths.get(out))))){
			throw new RuntimeException("SELFTEST FAILED: role manifest .sha256 sidecar does not match "
				+"an independent recomputation of the on-disk bytes.");
		}
		if(!duSha.equals(sha256Hex(java.nio.file.Files.readAllBytes(java.nio.file.Paths.get(du))))){
			throw new RuntimeException("SELFTEST FAILED: decoy_usage .sha256 sidecar does not match "
				+"an independent recomputation of the on-disk bytes.");
		}
		System.err.println("  selftest[provenance header + both output .sha256 sidecars]: PASS");

		// (7) group_partition_sha256 verification: corrupt a claimed hash, must throw.
		final String badSample=readWhole(dir+"/sample.tsv").replace(
			IdentityGroupOutputVerifier.hashFileRaw(dir+"/focalA.cluster.tsv"), new String(new char[64]).replace((char)0, '0'));
		writeFile(dir+"/sample_badhash.tsv", badSample);
		expectThrow(() -> process(dir+"/sample_badhash.tsv", dir+"/fold_manifest.tsv", dir+"/neighbors.tsv",
			dir+"/members.tsv", dir+"/out_bad1.tsv", dir+"/du_bad1.tsv"),
			"group_partition_sha256 mismatch rejected");

		// (8) non-contiguous repeated-family corruption in family_members_v1.tsv must be rejected.
		final String noncontigMembers=membersTsv.toString().replaceFirst(
			java.util.regex.Pattern.quote("focalA\tfa3\nfocalA\tfa4\n"),
			"decoyOnly\td1\nfocalA\tfa3\nfocalA\tfa4\n");//re-interleave decoyOnly mid-focalA block
		writeFile(dir+"/members_noncontig.tsv", noncontigMembers);
		expectThrow(() -> process(dir+"/sample.tsv", dir+"/fold_manifest.tsv", dir+"/neighbors.tsv",
			dir+"/members_noncontig.tsv", dir+"/out_bad2.tsv", dir+"/du_bad2.tsv"),
			"non-contiguous repeated-family in family_members_v1.tsv rejected");

		// (9) canonical list hash rejects a duplicate ID.
		expectThrow(() -> canonicalListSha256(new ArrayList<String>(java.util.Arrays.asList("a", "b", "a"))),
			"duplicate ID in canonical list rejected");

		// (10) canonical list hash rejects a non-ASCII ID (round 2, Elly's item 4).
		expectThrow(() -> canonicalListSha256(new ArrayList<String>(java.util.Arrays.asList("café"))),
			"non-ASCII ID in canonical list rejected");

		// (11) decoy member-count cross-check (round 2, Elly's item 3): family_members_v1.tsv
		// gains an EXTRA row for decoyOnly that sample_manifest.tsv's member_count doesn't know
		// about -- must be rejected, not silently trusted.
		final String extraMemberForDecoy=membersTsv.toString().replace(
			"decoyOnly\td2\n", "decoyOnly\td2\ndecoyOnly\td3\n");
		writeFile(dir+"/members_extra.tsv", extraMemberForDecoy);
		expectThrow(() -> process(dir+"/sample.tsv", dir+"/fold_manifest.tsv", dir+"/neighbors.tsv",
			dir+"/members_extra.tsv", dir+"/out_bad3.tsv", dir+"/du_bad3.tsv"),
			"decoy member-count cross-check rejected a family_members_v1.tsv/sample_manifest.tsv divergence");

		// (12) frozen K enforcement (round 2, Elly's item 2): a self-consistent-but-wrong K=19
		// neighbor file for focalA must be rejected outright.
		final ByteBuilder k19=new ByteBuilder();
		k19.append("#k\t19\n#focal_family\tneighbor_rank\tneighbor_family\n");
		for(int i=0; i<19; i++){k19.append("focalA\t").append(i+1).append('\t').append(focalANeighbors[i]).append('\n');}
		writeFile(dir+"/neighbors_k19.tsv", k19.toString());
		expectThrow(() -> process(dir+"/sample.tsv", dir+"/fold_manifest.tsv", dir+"/neighbors_k19.tsv",
			dir+"/members.tsv", dir+"/out_bad4.tsv", dir+"/du_bad4.tsv"),
			"declared K=19, not the frozen 20, rejected");

		// (13) duplicate neighbor_rank within a real K=20-declared block must be rejected.
		final ByteBuilder dupRank=new ByteBuilder();
		dupRank.append("#k\t").append(FROZEN_K).append("\n#focal_family\tneighbor_rank\tneighbor_family\n");
		for(int i=0; i<FROZEN_K; i++){
			final int rank=(i==FROZEN_K-1) ? FROZEN_K-1 : i+1;//last row reuses the previous rank
			dupRank.append("focalA\t").append(rank).append('\t').append(focalANeighbors[i]).append('\n');
		}
		writeFile(dir+"/neighbors_duprank.tsv", dupRank.toString());
		expectThrow(() -> process(dir+"/sample.tsv", dir+"/fold_manifest.tsv", dir+"/neighbors_duprank.tsv",
			dir+"/members.tsv", dir+"/out_bad5.tsv", dir+"/du_bad5.tsv"),
			"duplicate neighbor_rank in a real K=20 block rejected");

		// (14) bounded-memory proof (round 2, Elly's item 1): the peak SINGLE open accumulator
		// never exceeds the largest individual needed family (4, focalA) -- NOT the sum across all
		// 40 needed decoy families (2 + 19 + 4 + 19 = 44), which is what the old unbounded reader
		// would have implied by retaining every family's list simultaneously.
		final HashSet<String> allDecoyIds=new HashSet<String>();
		allDecoyIds.add("decoyOnly"); allDecoyIds.add("focalA");
		for(int i=0; i<19; i++){allDecoyIds.add("nbrA"+i); allDecoyIds.add("nbrB"+i);}
		final MemberStreamResult msr=streamNeededFamilyMembers(dir+"/members.tsv", allDecoyIds);
		if(msr.maxOpenListSize!=4){
			throw new RuntimeException("SELFTEST FAILED: peak single-family accumulator was "
				+msr.maxOpenListSize+", expected exactly 4 (focalA, the largest needed family) -- "
				+"bounded-memory reading is broken (a sum-of-all-needed-families bug would show 44).");
		}
		System.err.println("  selftest[bounded-memory family_members_v1.tsv reader]: PASS -- peak "
			+"single-family accumulator = 4 (focalA), never the 44-member sum across 40 needed "
			+"decoy families");

		// (15) determinism.
		final String out2=dir+"/role2.tsv", du2=dir+"/decoy_usage2.tsv";
		process(dir+"/sample.tsv", dir+"/fold_manifest.tsv", dir+"/neighbors.tsv", dir+"/members.tsv",
			out2, du2);
		if(!outText.equals(readWhole(out2))){
			throw new RuntimeException("SELFTEST FAILED: rerun did not reproduce byte-identical role manifest.");
		}
		System.err.println("  selftest[determinism]: PASS");

		System.err.println("SELFTEST PASS (dual-role family gets both focal AND decoy rows; plain "
			+"single-role decoy and focal-only families; decoy-usage cross-reference; exact training/"
			+"heldout member partition; provenance header + output sha256 sidecars; "
			+"group_partition_sha256 verification; non-contiguous family_members_v1.tsv corruption "
			+"rejected; duplicate-ID and non-ASCII-ID canonical-hash rejection; decoy member-count "
			+"cross-check; frozen-K=20 and duplicate-rank rejection; proven bounded-memory reading; "
			+"determinism -- all confirmed).");
	}

	interface ThrowingRunnable { void run() throws Exception; }
	static void expectThrow(ThrowingRunnable r, String label){
		boolean threw=false;
		try{ r.run(); }
		catch(Exception e){
			threw=true;
			System.err.println("    ["+label+"]: correctly caught -> "+e.getMessage());
		}
		if(!threw){throw new RuntimeException("SELFTEST FAILED: '"+label+"' should have thrown but did not.");}
	}

	static String neighborBlock(String focal, String[] neighbors){
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<neighbors.length; i++){
			sb.append(focal).append('\t').append(i+1).append('\t').append(neighbors[i]).append('\n');
		}
		return sb.toString();
	}

	/** Registers one selftest family into the shared accumulator lists used to build
	 * sample_manifest.tsv and family_members_v1.tsv fixtures together, so the two files can never
	 * silently disagree on membership the way a hand-duplicated pair of literals could. */
	static void registerFamily(ArrayList<String> repIds, ArrayList<Boolean> selected,
			ArrayList<Long> groupCounts, ArrayList<Long> memberCounts, ByteBuilder membersTsv,
			String repId, boolean isSelected, long groupCount, String[] members){
		repIds.add(repId); selected.add(isSelected); groupCounts.add(groupCount);
		memberCounts.add((long)members.length);
		for(String m : members){membersTsv.append(repId).append('\t').append(m).append('\n');}
	}

	/** Builds a real cluster.tsv+sidecar via IdentityGroupBuilder's actual writer for a family
	 * with the given intact groups (each inner array is one group's member IDs; a fully-connected
	 * edge clique within a group, no edges across groups, mirrors FoldManifestAggregator's own
	 * selftest fixture builder). */
	static void buildRealFamily(String dir, String repId, String[] allMembers, String[][] groups){
		final String ids=dir+"/"+repId+".ids";
		final String edges=dir+"/"+repId+"_edges.m8";
		final String cluster=dir+"/"+repId+".cluster.tsv";
		final StringBuilder idText=new StringBuilder();
		for(String m : allMembers){idText.append(m).append('\n');}
		writeFile(ids, idText.toString());
		final StringBuilder edgeText=new StringBuilder();
		for(String[] group : groups){
			for(int i=0; i<group.length; i++){
				for(int j=i+1; j<group.length; j++){
					edgeText.append(group[i]).append('\t').append(group[j])
						.append("\t100\t100\t1\t100\t100\t1\t100\t100\n");
				}
			}
		}
		writeFile(edges, edgeText.toString());
		IdentityGroupBuilder.main(new String[]{"ids="+ids, "edges="+edges, "out="+cluster,
			"minid=90", "minqcov=80", "mintcov=80"});
	}

	/** Builds a minimal real sample_manifest.tsv fixture (only the columns this tool reads),
	 * cluster_path/sidecar_path pointing at buildRealFamily's real output, group_partition_sha256
	 * the REAL raw hash of that cluster.tsv (so verifyPartition passes on the clean fixture), and
	 * member_count taken directly from registerFamily's own real member arrays (never re-derived
	 * from group_count, which is only valid when every group happens to have the same size). */
	static String sampleManifestRows(String dir, ArrayList<String> repIds, ArrayList<Boolean> selected,
			ArrayList<Long> groupCounts, ArrayList<Long> memberCounts){
		final StringBuilder sb=new StringBuilder();
		sb.append("#rep_id\tmember_count\tgroup_count\tselected\tcluster_path\tsidecar_path\t"
			+"group_partition_sha256\n");
		for(int i=0; i<repIds.size(); i++){
			final String cluster=dir+"/"+repIds.get(i)+".cluster.tsv";
			final String sidecar=cluster+".sidecar.tsv";
			sb.append(repIds.get(i)).append('\t').append(memberCounts.get(i)).append('\t')
				.append(groupCounts.get(i)).append('\t').append(selected.get(i) ? 1 : 0).append('\t')
				.append(cluster).append('\t').append(sidecar).append('\t')
				.append(IdentityGroupOutputVerifier.hashFileRaw(cluster)).append('\n');
		}
		return sb.toString();
	}
}
