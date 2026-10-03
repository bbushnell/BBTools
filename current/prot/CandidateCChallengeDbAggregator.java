package prot;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import structures.ByteBuilder;

/**
 * Candidate C Phase B2 (`CANDIDATE_C_BUILD_PIPELINE_DESIGN_v1.md` v4, sec 2/3b): the ONE non-array
 * job that collates Phase B's per-challenge-set {@code complete.sidecar} files into the single
 * published {@code candidate_c_challenge_db_manifest.tsv}. Mirrors {@link CandidateCArtifactAggregator}
 * (Phase A2) exactly, per the design's own text -- same independent-re-validation discipline, same
 * fail-closed required fields, same structured (never delimiter-string) uniqueness checks.
 * <p>
 * <b>The {@code complete.sidecar} contract (pinned here, since Phase B itself is a later delivery --
 * this is the contract it must honor, per sec 2 Phase B step 5):</b> a plain {@code key\tvalue} line
 * per field, written INSIDE the {@code <challenge_key>.hmmdb.d/} directory: {@code focal_family_id},
 * {@code fold_id}, {@code challenge_key}, {@code profile_count} (21 = 1 focal + {@value #FROZEN_K}
 * neighbors), {@code expected_artifact_id_set_sha256}, {@code expected_profile_order_sha256}
 * (BOTH upstream, RE-VALIDATED here against a fresh independent recomputation from the live
 * {@code artifact_role_manifest.tsv}/{@code family_neighbors.tsv} -- never assumed to still match just
 * because the sidecar exists), {@code concat_canonical_sha256}, {@code db_hmmdb_sha256},
 * {@code h3f_sha256}, {@code h3i_sha256}, {@code h3m_sha256}, {@code h3p_sha256} (ALL rehashed here
 * against the actual on-disk {@code db.hmmdb}/{@code db.hmmdb.h3{f,i,m,p}} files), {@code hmmpress_version}.
 * <p>
 * <b>Interpretive call, flagged for correction once Phase B actually exists (sec 3a/3b list
 * {@code concat_canonical_sha256} and {@code db_hmmdb_sha256} as two SEPARATE manifest columns):</b>
 * this phase treats them as the SAME underlying {@code db.hmmdb} file's hash recorded under two
 * provenance labels -- one being the build-time "what Phase B intended to publish" token, the other
 * the "what a downstream consumer reads" token -- and requires BOTH sidecar fields to independently
 * match the freshly-rehashed on-disk file. If Phase B's eventual design distinguishes them as two
 * physically different artifacts, this class's rehash step needs updating accordingly.
 * <p>
 * <b>Deterministic row order:</b> {@code focal_family_id} lexical, then {@code fold_id} numeric
 * ascending -- the SAME frozen key {@link CandidateCArtifactAggregator} and
 * {@link ArtifactRoleManifestGenerator} both use, for consistency across this project's manifests.
 * Uniqueness of that key is checked via a structured {@code HashMap<String,HashSet<Integer>>}, never
 * a delimiter-joined string (Elly's finding on A2 -- a hand-typed delimiter risks an invisible/control
 * byte slipping into source, as literally happened once in this same package).
 *
 * <p>
 * <b>Schema v4 / sec 1 staleness closure (`CANDIDATE_C_PHASE_AB_GATE_PLAN_v1.md` sec 6, Elly/UMP45
 * co-sealed):</b> everything above proves {@code db.hmmdb} is internally self-consistent with its OWN
 * sidecar and contains the correct 21 members in the correct order -- it does NOT prove those bytes are
 * still CURRENT, since Phase A/A2 can rebuild one or more {@code .hmm.canonical} files (e.g. after a
 * Phase-0 corpus edit) without changing any artifact ID, leaving a stale Phase-B directory that still
 * looks internally consistent. This class now REQUIRES {@code artifactmanifest=} (A2's current
 * {@code candidate_c_artifact_manifest.tsv}), resolves every challenge's 21 ids against it (validating
 * each row's {@code artifact_key} the same way {@link CandidateCChallengeDbBuilder} does: safe pattern,
 * recomputed-equal, live file non-empty, rehashed against the declared {@code hmm_canonical_sha256}),
 * rebuilds the current concat hash from those CURRENT bytes, and requires it equal the published
 * {@code db.hmmdb}'s actual hash. Output bumps to {@code #schema_version 4} and adds required header
 * provenance {@code #candidate_c_artifact_manifest_path}/{@code #candidate_c_artifact_manifest_sha256}.
 * The prior schema-3 seal is superseded for production by this contract change.
 *
 * <p>Usage: {@code java -ea prot.CandidateCChallengeDbAggregator rolemanifest=artifact_role_manifest.tsv
 *        neighbors=family_neighbors.tsv artifactmanifest=candidate_c_artifact_manifest.tsv
 *        sidecardir=<dir containing Phase-B's <challenge_key>.hmmdb.d/ subdirectories>
 *        out=candidate_c_challenge_db_manifest.tsv}
 * <br>Self-test: {@code java -ea prot.CandidateCChallengeDbAggregator selftest}
 *
 * @author Eru
 */
public final class CandidateCChallengeDbAggregator {

	/** The frozen candidate-independent neighbor-challenge K (see
	 * {@code ArtifactRoleManifestGenerator.FROZEN_K}/{@code NegativeQueryManifestGenerator}'s identical
	 * constant) -- not a tunable, ALWAYS 20. Independently redeclared per this project's convention of
	 * not sharing constants across tools via inheritance. */
	static final int FROZEN_K=20;
	/** Literal command template recorded once in the output header (sec 3b) -- documentation text,
	 * not a per-row substituted command. */
	static final String HMMPRESS_COMMAND_TEMPLATE="hmmpress db.hmmdb";

	public static void main(final String[] args) throws Exception{
		if(args.length==1 && args[0].equalsIgnoreCase("selftest")){selftest(); return;}
		String roleFile=null, neighborsFile=null, artifactManifestFile=null, sidecarDir=null, outFile=null;
		for(final String arg : args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase();
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("rolemanifest") || a.equals("role")){roleFile=b;}
			else if(a.equals("neighbors")){neighborsFile=b;}
			else if(a.equals("artifactmanifest")){artifactManifestFile=b;}
			else if(a.equals("sidecardir")){sidecarDir=b;}
			else if(a.equals("out")){outFile=b;}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}
		if(roleFile==null || neighborsFile==null || artifactManifestFile==null || sidecarDir==null
				|| outFile==null){
			throw new RuntimeException("Required: rolemanifest= neighbors= artifactmanifest= sidecardir= out=");
		}
		final Result r=process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, outFile);
		System.err.println(r.summary);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Data holders        ----------------*/
	/*--------------------------------------------------------------*/

	static class FocalRow { String artifactId, familyId; int foldId; }
	static class UniverseRow { String artifactId, familyId, role; int foldId; }
	static class RoleLoad {
		final ArrayList<FocalRow> focalRows=new ArrayList<FocalRow>();
		final HashMap<String,UniverseRow> byArtifactId=new HashMap<String,UniverseRow>();
	}
	static class RankedNeighbor { int rank; String family; }
	static class Sidecar {
		String focalFamilyId, challengeKey, expectedArtifactIdSetSha256, expectedProfileOrderSha256,
			concatCanonicalSha256, dbHmmdbSha256, h3fSha256, h3iSha256, h3mSha256, h3pSha256, hmmpressVersion;
		int foldId, profileCount;
	}
	static class ChallengeDbRow {
		String focalFamilyId, challengeKey, dbDir, expectedArtifactIdSetSha256, expectedProfileOrderSha256,
			concatCanonicalSha256, dbHmmdbSha256, h3fSha256, h3iSha256, h3mSha256, h3pSha256;
		int foldId, profileCount;
	}
	static class ArtifactManifestRow { String artifactId, artifactKey, hmmCanonicalPath, hmmCanonicalSha256; }
	static class Result { String summary; }

	/*--------------------------------------------------------------*/
	/*----------------           Pipeline           ----------------*/
	/*--------------------------------------------------------------*/

	static Result process(final String roleFile, final String neighborsFile, final String artifactManifestFile,
			final String sidecarDir, final String outFile) throws Exception{
		//Sec 1 step 5 (v4 schema amendment): the CURRENT candidate_c_artifact_manifest.tsv (A2's output),
		//loaded once up front exactly like roleFile/neighborsFile -- never per-challenge re-read, so every
		//challenge in one aggregation run is judged against the SAME snapshot of A2's output.
		final HashMap<String,ArtifactManifestRow> artifactManifest=loadArtifactManifest(artifactManifestFile);
		final RoleLoad roleLoad=loadRoleUniverse(roleFile);
		final ArrayList<FocalRow> focalRows=roleLoad.focalRows;

		//Structured (never delimiter-string) duplicate (family_id, fold_id) check, up front.
		final HashMap<String,HashSet<Integer>> foldsByFamily=new HashMap<String,HashSet<Integer>>();
		for(final FocalRow fr : focalRows){
			final HashSet<Integer> folds=foldsByFamily.computeIfAbsent(fr.familyId, k -> new HashSet<Integer>());
			if(!folds.add(fr.foldId)){
				throw new RuntimeException("Duplicate (family_id, fold_id) = ('"+fr.familyId+"', "+fr.foldId
					+") across more than one artifact_id in "+roleFile+".");
			}
		}

		final HashMap<String,ArrayList<String>> neighborsByFamily=loadNeighbors(neighborsFile);

		//Two-way set equality: distinct focal family_ids needing a challenge set vs family_neighbors.tsv's
		//own focal set -- never assume the two independently-produced files already agree.
		requireExactSet(foldsByFamily.keySet(), neighborsByFamily.keySet(),
			roleFile+"'s distinct focal_fold_local family_id set", neighborsFile);

		//Two-way set equality: expected challenge_key set vs the ACTUAL <key>.hmmdb.d directories found
		//on disk -- an extra directory (e.g. from a superseded manifest) is a hard failure, never
		//silently ignored, mirroring Phase A2's own convention.
		final HashSet<String> expectedChallengeKeys=new HashSet<String>();
		for(final FocalRow fr : focalRows){
			expectedChallengeKeys.add(CandidateCUtil.challengeKey(CandidateCUtil.artifactKey(fr.artifactId)));
		}
		final HashSet<String> actualChallengeKeys=listChallengeDbDirs(sidecarDir);
		requireExactSet(expectedChallengeKeys, actualChallengeKeys, "expected challenge_key set",
			sidecarDir+"/*.hmmdb.d");

		final ArrayList<ChallengeDbRow> results=new ArrayList<ChallengeDbRow>();
		String commonHmmpressVersion=null;
		for(final FocalRow fr : focalRows){
			final String artifactKey=CandidateCUtil.artifactKey(fr.artifactId);
			final String challengeKey=CandidateCUtil.challengeKey(artifactKey);
			final ArrayList<String> neighborFamilies=neighborsByFamily.get(fr.familyId);

			//Focal artifact_id first, then neighbor decoy artifact_ids in ascending neighbor_rank --
			//the pinned focal-then-neighbor-rank concat order (sec 3b/3c's construction style). Each
			//neighbor's decoy artifact is NEVER assumed to exist: it must be present in the full role
			//universe with role=="all_member_decoy" and fold_id==-1 -- Elly's finding. This loop also
			//builds the key->id translation map used later to independently re-derive db.hmmdb's
			//ACTUAL content, scoped to exactly this challenge set's expected artifacts.
			final ArrayList<String> orderedProfileIds=new ArrayList<String>(neighborFamilies.size()+1);
			final HashMap<String,String> keyToId=new HashMap<String,String>();
			orderedProfileIds.add(fr.artifactId);
			keyToId.put(artifactKey, fr.artifactId);
			for(final String neighborFamily : neighborFamilies){
				final String decoyId=neighborFamily+"__decoy";
				final UniverseRow decoyRow=roleLoad.byArtifactId.get(decoyId);
				if(decoyRow==null){
					throw new RuntimeException("Focal '"+fr.artifactId+"': neighbor family '"+neighborFamily
						+"'s required decoy artifact '"+decoyId+"' does not exist in "+roleFile+".");
				}
				if(!decoyRow.role.equals("all_member_decoy")){
					throw new RuntimeException("Focal '"+fr.artifactId+"': neighbor decoy artifact '"
						+decoyId+"' has role='"+decoyRow.role+"' in "+roleFile+", expected 'all_member_decoy'.");
				}
				if(!decoyRow.familyId.equals(neighborFamily)){
					throw new RuntimeException("Focal '"+fr.artifactId+"': neighbor decoy artifact '"
						+decoyId+"' has family_id='"+decoyRow.familyId+"' in "+roleFile+", expected '"
						+neighborFamily+"' -- a malformed row can otherwise claim this artifact_id while "
						+"actually representing a different family.");
				}
				if(decoyRow.foldId!=-1){
					throw new RuntimeException("Focal '"+fr.artifactId+"': neighbor decoy artifact '"
						+decoyId+"' has fold_id="+decoyRow.foldId+" in "+roleFile+", expected -1.");
				}
				orderedProfileIds.add(decoyId);
				keyToId.put(CandidateCUtil.artifactKey(decoyId), decoyId);
			}
			final String expectedOrderHash=CandidateCUtil.orderedIdListSha256(orderedProfileIds);
			final String expectedSetHash=ArtifactRoleManifestGenerator.canonicalListSha256(orderedProfileIds);

			final String dbDir=sidecarDir+"/"+challengeKey+".hmmdb.d";
			final Sidecar sc=loadSidecar(dbDir+"/complete.sidecar");

			if(!sc.challengeKey.equals(challengeKey)){
				throw new RuntimeException("Sidecar in '"+dbDir+"' declares challenge_key='"+sc.challengeKey
					+"', which does not match its own directory name.");
			}
			if(!sc.focalFamilyId.equals(fr.familyId)){
				throw new RuntimeException("Sidecar '"+challengeKey+"' declares focal_family_id='"
					+sc.focalFamilyId+"' but the role manifest maps this challenge_key to family '"
					+fr.familyId+"'.");
			}
			if(sc.foldId!=fr.foldId){
				throw new RuntimeException("Sidecar '"+challengeKey+"' declares fold_id="+sc.foldId
					+" but the role manifest declares fold_id="+fr.foldId+" for this artifact.");
			}
			if(sc.profileCount!=orderedProfileIds.size()){
				throw new RuntimeException("Sidecar '"+challengeKey+"' declares profile_count="
					+sc.profileCount+" but the expected challenge set has "+orderedProfileIds.size()
					+" profiles (1 focal + "+neighborFamilies.size()+" neighbors).");
			}

			//Independent re-validation, never assume Phase B's own resume check ran correctly or ran at
			//all for a stale sidecar: recompute the CURRENT expected set/order hashes from the live
			//artifact_role_manifest.tsv/family_neighbors.tsv and compare -- a mismatch means the
			//challenge composition changed (e.g. the neighbor list was regenerated) but Phase B was
			//never rerun to pick it up.
			if(!sc.expectedArtifactIdSetSha256.equals(expectedSetHash)){
				throw new RuntimeException("Challenge '"+challengeKey+"': sidecar expected_artifact_id_set_sha256='"
					+sc.expectedArtifactIdSetSha256+"' does not match the CURRENT independently-recomputed "
					+"value '"+expectedSetHash+"' -- the challenge set's membership changed, but Phase B "
					+"was never rerun to pick it up.");
			}
			if(!sc.expectedProfileOrderSha256.equals(expectedOrderHash)){
				throw new RuntimeException("Challenge '"+challengeKey+"': sidecar expected_profile_order_sha256='"
					+sc.expectedProfileOrderSha256+"' does not match the CURRENT independently-recomputed "
					+"value '"+expectedOrderHash+"' -- the challenge set's neighbor order changed, but "
					+"Phase B was never rerun to pick it up.");
			}

			//Non-empty requirement BEFORE hashing (Elly's finding, sec 2 Phase B step 4's own
			//requirement): an empty file whose hash happens to match a correspondingly "empty"
			//sidecar-recorded value would otherwise pass the byte-rehash check silently.
			final String dbHmmdbPath=dbDir+"/db.hmmdb";
			final String h3fPath=dbDir+"/db.hmmdb.h3f";
			final String h3iPath=dbDir+"/db.hmmdb.h3i";
			final String h3mPath=dbDir+"/db.hmmdb.h3m";
			final String h3pPath=dbDir+"/db.hmmdb.h3p";
			requireNonEmpty(dbHmmdbPath);
			requireNonEmpty(h3fPath);
			requireNonEmpty(h3iPath);
			requireNonEmpty(h3mPath);
			requireNonEmpty(h3pPath);

			//Independent re-hash of the actual on-disk pressed-database files -- a sidecar recording a
			//hash that no longer matches its own referenced file is rejected, never republished.
			final String actualDbHmmdbSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(dbHmmdbPath);
			requireRehashMatch(challengeKey, dbHmmdbPath, "concat_canonical_sha256", sc.concatCanonicalSha256,
				actualDbHmmdbSha256);
			requireRehashMatch(challengeKey, dbHmmdbPath, "db_hmmdb_sha256", sc.dbHmmdbSha256, actualDbHmmdbSha256);
			final String actualH3f=ArtifactRoleManifestGenerator.streamingSha256Hex(h3fPath);
			requireRehashMatch(challengeKey, h3fPath, "h3f_sha256", sc.h3fSha256, actualH3f);
			final String actualH3i=ArtifactRoleManifestGenerator.streamingSha256Hex(h3iPath);
			requireRehashMatch(challengeKey, h3iPath, "h3i_sha256", sc.h3iSha256, actualH3i);
			final String actualH3m=ArtifactRoleManifestGenerator.streamingSha256Hex(h3mPath);
			requireRehashMatch(challengeKey, h3mPath, "h3m_sha256", sc.h3mSha256, actualH3m);
			final String actualH3p=ArtifactRoleManifestGenerator.streamingSha256Hex(h3pPath);
			requireRehashMatch(challengeKey, h3pPath, "h3p_sha256", sc.h3pSha256, actualH3p);

			//Semantic content validation (Elly's finding): byte-rehashing alone proves db.hmmdb's bytes
			//match what the sidecar SAYS they are, never that those bytes are actually the CORRECT 21
			//profiles in the CORRECT order. Independently parse db.hmmdb's NAME lines in physical file
			//order and translate each back to its artifact_id via the per-challenge key->id map built
			//above (an unrecognized key is an immediate hard failure).
			final ArrayList<String> rawNames=parseHmmdbNamesInOrder(dbHmmdbPath);
			final ArrayList<String> actualIdsInFile=new ArrayList<String>(rawNames.size());
			for(final String key : rawNames){
				final String id=keyToId.get(key);
				if(id==null){
					throw new RuntimeException("Challenge '"+challengeKey+"': db.hmmdb at "+dbHmmdbPath
						+" contains a profile NAME '"+key+"' that is not one of this challenge set's "
						+orderedProfileIds.size()+" expected artifacts.");
				}
				actualIdsInFile.add(id);
			}
			if(actualIdsInFile.size()!=orderedProfileIds.size()){
				throw new RuntimeException("Challenge '"+challengeKey+"': db.hmmdb at "+dbHmmdbPath
					+" contains "+actualIdsInFile.size()+" profile(s) but the expected challenge set has "
					+orderedProfileIds.size()+" -- a duplicate NAME could mask a missing one even if the "
					+"membership SET happens to look complete.");
			}
			final HashSet<String> actualIdSet=new HashSet<String>(actualIdsInFile);
			final HashSet<String> expectedIdSet=new HashSet<String>(orderedProfileIds);
			requireExactSet(expectedIdSet, actualIdSet, "expected 21-artifact challenge set",
				"db.hmmdb's actual NAME-line profiles ("+dbHmmdbPath+")");
			final String actualOrderHash=CandidateCUtil.orderedIdListSha256(actualIdsInFile);
			if(!actualOrderHash.equals(expectedOrderHash)){
				throw new RuntimeException("Challenge '"+challengeKey+"': db.hmmdb's ACTUAL physical "
					+"NAME-line order at "+dbHmmdbPath+" does not match the expected focal-then-neighbor-"
					+"rank order (same membership, wrong order) -- hash '"+actualOrderHash
					+"' vs expected '"+expectedOrderHash+"'.");
			}

			if(commonHmmpressVersion==null){commonHmmpressVersion=sc.hmmpressVersion;}
			else if(!commonHmmpressVersion.equals(sc.hmmpressVersion)){
				throw new RuntimeException("Challenge '"+challengeKey+"': sidecar records hmmpress_version='"
					+sc.hmmpressVersion+"' but an earlier challenge set in this same build recorded '"
					+commonHmmpressVersion+"' -- every challenge set in one aggregation must share one "
					+"tool version.");
			}

			//Sec 1 (v4 schema amendment, the current-HMM-byte staleness closure): everything checked
			//above proves db.hmmdb is internally self-consistent with its OWN sidecar and contains the
			//correct 21 members in the correct order -- it does NOT prove those bytes are still CURRENT.
			//Phase A/A2 can rebuild one or more .hmm.canonical files (e.g. after a Phase-0 corpus edit)
			//without changing any artifact ID, leaving a stale Phase-B directory that still looks
			//internally consistent. Resolve every one of this challenge's 21 ids against the CURRENT
			//artifact manifest, rehash each live .hmm.canonical against its declared hash, then rebuild
			//the concat hash from those CURRENT bytes and require it to equal the published db.hmmdb's
			//actual hash -- catching the case where Phase A/A2 reran but Phase B never did.
			final ArrayList<String> hmmCanonicalPathsInOrder=new ArrayList<String>(orderedProfileIds.size());
			for(final String id : orderedProfileIds){
				final ArtifactManifestRow amr=artifactManifest.get(id);
				if(amr==null){
					throw new RuntimeException("Challenge '"+challengeKey+"': artifact '"+id+"' has no row "
						+"in "+artifactManifestFile+" -- Phase A2 has not produced it yet.");
				}
				if(!CandidateCUtil.ARTIFACT_KEY_PATTERN.matcher(amr.artifactKey).matches()){
					throw new RuntimeException("Challenge '"+challengeKey+"': "+artifactManifestFile
						+" declares artifact '"+id+"' has artifact_key='"+amr.artifactKey+"' which does not "
						+"match the required safe pattern "+CandidateCUtil.ARTIFACT_KEY_PATTERN+".");
				}
				final String recomputedKey=CandidateCUtil.artifactKey(id);
				if(!recomputedKey.equals(amr.artifactKey)){
					throw new RuntimeException("Challenge '"+challengeKey+"': "+artifactManifestFile
						+" declares artifact '"+id+"' has artifact_key='"+amr.artifactKey+"' but the "
						+"independently-recomputed FNV key is '"+recomputedKey+"' -- refusing to trust an "
						+"unverified copied key.");
				}
				if(!new java.io.File(amr.hmmCanonicalPath).exists() || new java.io.File(amr.hmmCanonicalPath).length()<=0){
					throw new RuntimeException("Challenge '"+challengeKey+"': the live .hmm.canonical at "
						+amr.hmmCanonicalPath+" declared by "+artifactManifestFile+" for artifact '"+id
						+"' is missing or empty.");
				}
				final String actualHmmCanonicalSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(
					amr.hmmCanonicalPath);
				if(!actualHmmCanonicalSha256.equals(amr.hmmCanonicalSha256)){
					throw new RuntimeException("Challenge '"+challengeKey+"': "+artifactManifestFile
						+" declares artifact '"+id+"' has hmm_canonical_sha256='"+amr.hmmCanonicalSha256
						+"' but the ACTUAL on-disk bytes at "+amr.hmmCanonicalPath+" hash to '"
						+actualHmmCanonicalSha256+"' -- stale or corrupt artifact-manifest row.");
				}
				hmmCanonicalPathsInOrder.add(amr.hmmCanonicalPath);
			}
			final String liveConcatSha256=streamingConcatSha256Hex(hmmCanonicalPathsInOrder);
			if(!liveConcatSha256.equals(actualDbHmmdbSha256)){
				throw new RuntimeException("Challenge '"+challengeKey+"': the CURRENT live concatenation of "
					+"this challenge's 21 .hmm.canonical files (per "+artifactManifestFile+") hashes to '"
					+liveConcatSha256+"', but the published db.hmmdb at "+dbHmmdbPath+" hashes to '"
					+actualDbHmmdbSha256+"' -- the artifact IDs/order are unchanged, but the underlying HMM "
					+"bytes changed (Phase A/A2 reran) and Phase B was never rerun to pick it up.");
			}

			final ChallengeDbRow row=new ChallengeDbRow();
			row.focalFamilyId=fr.familyId; row.foldId=fr.foldId; row.challengeKey=challengeKey;
			row.dbDir=dbDir; row.profileCount=orderedProfileIds.size();
			row.expectedArtifactIdSetSha256=expectedSetHash; row.expectedProfileOrderSha256=expectedOrderHash;
			row.concatCanonicalSha256=actualDbHmmdbSha256; row.dbHmmdbSha256=actualDbHmmdbSha256;
			row.h3fSha256=actualH3f; row.h3iSha256=actualH3i; row.h3mSha256=actualH3m; row.h3pSha256=actualH3p;
			results.add(row);
		}

		//Frozen sort order, matching CandidateCArtifactAggregator's/ArtifactRoleManifestGenerator's own
		//convention exactly: focal_family_id lexical, then fold_id numeric ascending. Uniqueness of this
		//key was already proven above (foldsByFamily), so the sort itself is safe to apply directly.
		results.sort((a, b) -> {
			final int c=a.focalFamilyId.compareTo(b.focalFamilyId);
			return (c!=0) ? c : Integer.compare(a.foldId, b.foldId);
		});

		writeManifest(outFile, results, roleFile, neighborsFile, artifactManifestFile,
			commonHmmpressVersion==null ? "-" : commonHmmpressVersion);

		final Result res=new Result();
		res.summary="CandidateCChallengeDbAggregator: "+results.size()+" challenge sets aggregated. "
			+"Manifest written to "+outFile+".";
		return res;
	}

	static void requireRehashMatch(final String challengeKey, final String path, final String fieldName,
			final String declared, final String actual){
		if(!declared.equals(actual)){
			throw new RuntimeException("Challenge '"+challengeKey+"': sidecar declares "+fieldName+"='"
				+declared+"' but the ACTUAL on-disk bytes at "+path+" hash to '"+actual+"' -- the file was "
				+"corrupted or modified after its sidecar was written.");
		}
	}

	/** Streams the exact raw-byte concatenation of the given files, in the given order, computing its
	 * SHA-256 WITHOUT ever holding more than one I/O buffer at a time (no reconstruction) --
	 * independently reimplemented here (matching {@link CandidateCChallengeDbBuilder}'s own copy, not
	 * calling it) per this project's re-derivation discipline for domain-logic staleness proofs. */
	static String streamingConcatSha256Hex(final ArrayList<String> pathsInOrder) throws Exception{
		final java.security.MessageDigest md=java.security.MessageDigest.getInstance("SHA-256");
		final byte[] buf=new byte[1<<20];
		for(final String path : pathsInOrder){
			try(java.io.BufferedInputStream in=new java.io.BufferedInputStream(
					new java.io.FileInputStream(path), 1<<20)){
				int n;
				while((n=in.read(buf))>=0){md.update(buf, 0, n);}
			}
		}
		return IdentityGroupOutputVerifier.toHex(md.digest());
	}

	/** Reads candidate_c_artifact_manifest.tsv (A2's output) by HEADER NAME, keyed by artifact_id --
	 * independently reimplemented here (matching {@link CandidateCChallengeDbBuilder}'s own copy, not
	 * calling it) per this project's re-derivation discipline. Rejects a duplicate artifact_id outright,
	 * never last-row-wins. */
	static HashMap<String,ArtifactManifestRow> loadArtifactManifest(final String file){
		final HashMap<String,ArtifactManifestRow> byArtifactId=new HashMap<String,ArtifactManifestRow>();
		final ByteFile bf=ByteFile.makeByteFile(file, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		HashMap<String,Integer> col=null;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='#'){
				final String h=new String(line, 1, line.length-1);
				final String[] names=h.split("\t");
				if(names.length>0 && names[0].trim().equalsIgnoreCase("artifact_id")){
					if(col!=null){throw new RuntimeException("Duplicate artifact_id header in "+file);}
					col=new HashMap<String,Integer>();
					for(int i=0; i<names.length; i++){col.put(names[i].trim().toLowerCase(), i);}
				}
				continue;
			}
			if(col==null){throw new RuntimeException("Data row before header in "+file);}
			lp.set(line);
			final ArtifactManifestRow r=new ArtifactManifestRow();
			r.artifactId=lp.parseString(col.get("artifact_id"));
			r.artifactKey=lp.parseString(col.get("artifact_key"));
			r.hmmCanonicalPath=lp.parseString(col.get("hmm_canonical_path"));
			r.hmmCanonicalSha256=lp.parseString(col.get("hmm_canonical_sha256"));
			if(byArtifactId.put(r.artifactId, r)!=null){
				throw new RuntimeException("Duplicate artifact_id '"+r.artifactId+"' in "+file+" -- "
					+"artifact_id must be unique.");
			}
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+file);}
		if(col==null){throw new RuntimeException("No artifact_id header found in "+file);}
		return byArtifactId;
	}

	/** Elly's finding: an empty file must never silently pass just because its hash happens to match a
	 * correspondingly "empty" sidecar-recorded value -- checked BEFORE hashing, matching sec 2 Phase B
	 * step 4's own "confirm all 4 .h3? files exist and are non-empty" requirement (extended here to
	 * db.hmmdb too, for the identical reason). */
	static void requireNonEmpty(final String path){
		final java.io.File f=new java.io.File(path);
		if(!f.exists()){throw new RuntimeException("Required file "+path+" does not exist.");}
		if(f.length()<=0){
			throw new RuntimeException("Required file "+path+" must be non-empty, but is "+f.length()
				+" bytes.");
		}
	}

	/** Independently parses a raw HMMER ASCII {@code .hmm}/{@code .hmmdb} file's {@code NAME} lines in
	 * the EXACT order they physically appear (never sorted) -- the semantic content validation this
	 * class performs depends on the ACTUAL file order, not any assumption about how it was built. */
	static ArrayList<String> parseHmmdbNamesInOrder(final String path){
		final ArrayList<String> names=new ArrayList<String>();
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			final String s=new String(line);
			if(s.startsWith("NAME  ")){
				final String name=s.substring(6).trim();
				if(name.isEmpty()){
					throw new RuntimeException("Malformed NAME line (empty name) in "+path+": "+s);
				}
				names.add(name);
			}
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+path);}
		return names;
	}

	/** Two-way set equality between two ID sets, matching
	 * {@code ArtifactRoleManifestGenerator.requireExactSet}'s established convention (independently
	 * re-implemented, not called -- a shared UTILITY shape, this class's own copy). */
	static void requireExactSet(final java.util.Set<String> expected, final java.util.Set<String> actual,
			final String expectedLabel, final String actualLabel){
		final ArrayList<String> missing=new ArrayList<String>(), extra=new ArrayList<String>();
		for(final String id : expected){if(!actual.contains(id)){missing.add(id);}}
		for(final String id : actual){if(!expected.contains(id)){extra.add(id);}}
		if(!missing.isEmpty() || !extra.isEmpty()){
			throw new RuntimeException(expectedLabel+" and "+actualLabel+" do not exactly match -- "
				+"missing="+missing.size()+" ("+(missing.isEmpty() ? "" : missing.get(0))+"), extra="
				+extra.size()+" ("+(extra.isEmpty() ? "" : extra.get(0))+").");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Loaders            ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads artifact_role_manifest.tsv by HEADER NAME, retaining the FULL row universe (every
	 * artifact_id/family_id/role/fold_id, whatever the role) -- Elly's finding: a neighbor family's
	 * decoy artifact must never be merely ASSUMED to exist just because its ID can be constructed; this
	 * loader is what lets that existence (and its role/fold_id) actually be proven. {@code focalRows}
	 * remains the {@code role=="focal_fold_local"} subset -- the pipeline's ultimate source of "which
	 * (family_id, fold_id) challenge sets must exist." */
	static RoleLoad loadRoleUniverse(final String file){
		final RoleLoad result=new RoleLoad();
		final ByteFile bf=ByteFile.makeByteFile(file, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		HashMap<String,Integer> col=null;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='#'){
				final String h=new String(line, 1, line.length-1);
				final String[] names=h.split("\t");
				if(names.length>0 && names[0].trim().equalsIgnoreCase("artifact_id")){
					if(col!=null){throw new RuntimeException("Duplicate artifact_id header in "+file);}
					col=new HashMap<String,Integer>();
					for(int i=0; i<names.length; i++){col.put(names[i].trim().toLowerCase(), i);}
				}
				continue;
			}
			if(col==null){throw new RuntimeException("Data row before header in "+file);}
			lp.set(line);
			final UniverseRow u=new UniverseRow();
			u.artifactId=lp.parseString(col.get("artifact_id"));
			u.familyId=lp.parseString(col.get("family_id"));
			u.role=lp.parseString(col.get("role"));
			u.foldId=lp.parseInt(col.get("fold_id"));
			if(result.byArtifactId.put(u.artifactId, u)!=null){
				throw new RuntimeException("Duplicate artifact_id '"+u.artifactId+"' in "+file);
			}
			if(u.role.equals("focal_fold_local")){
				final FocalRow r=new FocalRow();
				r.artifactId=u.artifactId; r.familyId=u.familyId; r.foldId=u.foldId;
				result.focalRows.add(r);
			}
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+file);}
		if(col==null){throw new RuntimeException("No artifact_id header found in "+file);}
		return result;
	}

	/** Reads family_neighbors.tsv (`FamilyNeighborManifest`'s output): a single {@code #k\t20} header,
	 * rows {@code focal_family\tneighbor_rank\tneighbor_family} in ANY file order -- independently
	 * reimplemented with the SAME rigor {@code ArtifactRoleManifestGenerator.loadNeighbors} already
	 * applies (declared K required to equal the FROZEN {@value #FROZEN_K}, every focal exactly K rows
	 * with neighbor_rank 1..K each appearing exactly once, non-self, no duplicate neighbor), PLUS this
	 * class's own requirement: the returned list is explicitly SORTED by neighbor_rank ascending (never
	 * assumed to already be in file order), since the order-sensitive
	 * {@code expected_profile_order_sha256} hash depends on it. */
	static HashMap<String,ArrayList<String>> loadNeighbors(final String file){
		final HashMap<String,ArrayList<RankedNeighbor>> raw=new HashMap<String,ArrayList<RankedNeighbor>>();
		final HashMap<String,HashSet<Integer>> ranksSeen=new HashMap<String,HashSet<Integer>>();
		final HashMap<String,HashSet<String>> neighborsSeen=new HashMap<String,HashSet<String>>();
		int declaredK=-1;
		final ByteFile bf=ByteFile.makeByteFile(file, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='#'){
				final String s=new String(line);
				if(s.startsWith("#k\t")){
					if(declaredK>=0){
						throw new RuntimeException(file+" has more than one '#k\\t<K>' header line -- "
							+"ambiguous which is authoritative.");
					}
					declaredK=Integer.parseInt(s.substring(3).trim());
				}
				continue;
			}
			lp.set(line);
			if(lp.terms()<3){throw new RuntimeException("Malformed neighbor row in "+file+": "+new String(line));}
			final String focal=lp.parseString(0);
			final int rank=lp.parseInt(1);
			final String neighbor=lp.parseString(2);
			if(neighbor.equals(focal)){
				throw new RuntimeException("Focal family '"+focal+"' lists itself as its own neighbor "
					+"(rank "+rank+") in "+file+" -- neighbors must be non-self.");
			}
			if(!ranksSeen.computeIfAbsent(focal, k -> new HashSet<Integer>()).add(rank)){
				throw new RuntimeException("Focal family '"+focal+"' has neighbor_rank "+rank
					+" recorded more than once in "+file+".");
			}
			if(!neighborsSeen.computeIfAbsent(focal, k -> new HashSet<String>()).add(neighbor)){
				throw new RuntimeException("Focal family '"+focal+"' has neighbor family '"+neighbor
					+"' recorded more than once (duplicate neighbor) in "+file+".");
			}
			final RankedNeighbor rn=new RankedNeighbor();
			rn.rank=rank; rn.family=neighbor;
			raw.computeIfAbsent(focal, k -> new ArrayList<RankedNeighbor>()).add(rn);
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+file);}
		if(declaredK<1){
			throw new RuntimeException(file+" has no valid '#k\\t<K>' header line -- cannot validate "
				+"neighbor completeness without the declared K.");
		}
		if(declaredK!=FROZEN_K){
			throw new RuntimeException(file+" declares K="+declaredK+", but this assay's K is FROZEN at "
				+FROZEN_K+" -- refusing to proceed with a non-frozen neighbor count.");
		}
		final HashMap<String,ArrayList<String>> ordered=new HashMap<String,ArrayList<String>>();
		for(final java.util.Map.Entry<String,ArrayList<RankedNeighbor>> e : raw.entrySet()){
			final String focal=e.getKey();
			final ArrayList<RankedNeighbor> list=e.getValue();
			if(list.size()!=declaredK){
				throw new RuntimeException("Focal family '"+focal+"' has "+list.size()+" neighbor rows in "
					+file+", expected exactly the declared K="+declaredK+".");
			}
			final HashSet<Integer> ranks=ranksSeen.get(focal);
			for(int r=1; r<=declaredK; r++){
				if(!ranks.contains(r)){
					throw new RuntimeException("Focal family '"+focal+"' is missing neighbor_rank "+r
						+" in "+file+" (ranks must be exactly 1.."+declaredK+", each appearing once).");
				}
			}
			list.sort((a, b) -> Integer.compare(a.rank, b.rank));
			final ArrayList<String> families=new ArrayList<String>(list.size());
			for(final RankedNeighbor rn : list){families.add(rn.family);}
			ordered.put(focal, families);
		}
		return ordered;
	}

	/** Lists every {@code *.hmmdb.d} directory's challenge_key (its dirname minus the extension) --
	 * the ACTUAL on-disk set, for two-way completeness checking. */
	static HashSet<String> listChallengeDbDirs(final String sidecarDir){
		final HashSet<String> keys=new HashSet<String>();
		final java.io.File[] files=new java.io.File(sidecarDir).listFiles();
		if(files!=null){
			for(final java.io.File f : files){
				if(!f.isDirectory()){continue;}
				final String name=f.getName();
				if(name.endsWith(".hmmdb.d")){keys.add(name.substring(0, name.length()-".hmmdb.d".length()));}
			}
		}
		return keys;
	}

	/** Reads one {@code complete.sidecar} (plain {@code key\tvalue} lines), rejecting a duplicate field
	 * outright (never last-value-wins -- Elly's finding on Phase A2, applied here from the start). */
	static Sidecar loadSidecar(final String path){
		final HashMap<String,String> fields=new HashMap<String,String>();
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			lp.set(line);
			if(lp.terms()<2){throw new RuntimeException("Malformed sidecar row in "+path+": "+new String(line));}
			final String key=lp.parseString(0);
			if(fields.put(key, lp.parseString(1))!=null){
				throw new RuntimeException("Sidecar "+path+" has duplicate field '"+key+"'.");
			}
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+path);}
		final Sidecar sc=new Sidecar();
		sc.focalFamilyId=requireField(fields, "focal_family_id", path);
		sc.foldId=Integer.parseInt(requireField(fields, "fold_id", path));
		sc.challengeKey=requireField(fields, "challenge_key", path);
		sc.profileCount=Integer.parseInt(requireField(fields, "profile_count", path));
		sc.expectedArtifactIdSetSha256=requireField(fields, "expected_artifact_id_set_sha256", path);
		sc.expectedProfileOrderSha256=requireField(fields, "expected_profile_order_sha256", path);
		sc.concatCanonicalSha256=requireField(fields, "concat_canonical_sha256", path);
		sc.dbHmmdbSha256=requireField(fields, "db_hmmdb_sha256", path);
		sc.h3fSha256=requireField(fields, "h3f_sha256", path);
		sc.h3iSha256=requireField(fields, "h3i_sha256", path);
		sc.h3mSha256=requireField(fields, "h3m_sha256", path);
		sc.h3pSha256=requireField(fields, "h3p_sha256", path);
		sc.hmmpressVersion=requireField(fields, "hmmpress_version", path);
		return sc;
	}

	private static String requireField(final HashMap<String,String> fields, final String key,
			final String path){
		final String v=fields.get(key);
		if(v==null || v.trim().isEmpty()){
			throw new RuntimeException("Sidecar "+path+" is missing required field '"+key+"'.");
		}
		return v;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	static void writeManifest(final String outFile, final ArrayList<ChallengeDbRow> results,
			final String roleFile, final String neighborsFile, final String artifactManifestFile,
			final String hmmpressVersion) throws Exception{
		final ByteBuilder out=new ByteBuilder(1<<20);
		out.append("#schema_version\t4\n");
		out.append("#artifact_role_manifest_path\t").append(roleFile)
			.append("\t#artifact_role_manifest_sha256\t")
			.append(ArtifactRoleManifestGenerator.streamingSha256Hex(roleFile)).append('\n');
		out.append("#family_neighbors_path\t").append(neighborsFile)
			.append("\t#family_neighbors_sha256\t")
			.append(ArtifactRoleManifestGenerator.streamingSha256Hex(neighborsFile)).append('\n');
		out.append("#candidate_c_artifact_manifest_path\t").append(artifactManifestFile)
			.append("\t#candidate_c_artifact_manifest_sha256\t")
			.append(ArtifactRoleManifestGenerator.streamingSha256Hex(artifactManifestFile)).append('\n');
		out.append("#hmmpress_version\t").append(hmmpressVersion).append("\t#hmmpress_command\t")
			.append(HMMPRESS_COMMAND_TEMPLATE).append('\n');
		out.append("#focal_family_id\tfold_id\tchallenge_key\tdb_dir\tprofile_count\t"
			+"expected_artifact_id_set_sha256\texpected_profile_order_sha256\tconcat_canonical_sha256\t"
			+"db_hmmdb_sha256\th3f_sha256\th3i_sha256\th3m_sha256\th3p_sha256\n");
		for(final ChallengeDbRow r : results){
			out.append(r.focalFamilyId).append('\t').append(r.foldId).append('\t').append(r.challengeKey)
				.append('\t').append(r.dbDir).append('\t').append(r.profileCount).append('\t')
				.append(r.expectedArtifactIdSetSha256).append('\t').append(r.expectedProfileOrderSha256)
				.append('\t').append(r.concatCanonicalSha256).append('\t').append(r.dbHmmdbSha256)
				.append('\t').append(r.h3fSha256).append('\t').append(r.h3iSha256).append('\t')
				.append(r.h3mSha256).append('\t').append(r.h3pSha256).append('\n');
		}
		writeFileAtomic(outFile, out.toString());
		writeSha256SidecarAtomic(outFile);
	}

	static void writeFileAtomic(final String path, final String content) throws Exception{
		final String tmp=path+".tmp";
		final FileFormat ff=FileFormat.testOutput(tmp, FileFormat.TXT, null, false, true, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		bsw.print(new ByteBuilder().append(content));
		if(bsw.poisonAndWait()){throw new RuntimeException("I/O error writing "+tmp);}
		Files.move(Paths.get(tmp), Paths.get(path), StandardCopyOption.REPLACE_EXISTING,
			StandardCopyOption.ATOMIC_MOVE);
	}

	static void writeSha256SidecarAtomic(final String path) throws Exception{
		final byte[] bytes=Files.readAllBytes(Paths.get(path));
		final String hash=ArtifactRoleManifestGenerator.sha256Hex(bytes);
		writeFileAtomic(path+".sha256", hash+"  "+path+"\n");
	}

	static String readWhole(final String path){
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		final StringBuilder sb=new StringBuilder();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){sb.append(new String(line)).append('\n');}
		bf.close();
		return sb.toString();
	}

	private static void writeFileForFixture(final String path, final String content){
		final FileFormat ff=FileFormat.testOutput(path, FileFormat.TXT, null, false, true, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		bsw.print(new ByteBuilder().append(content));
		if(bsw.poisonAndWait()){throw new RuntimeException("I/O error writing "+path);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Self-test            ----------------*/
	/*--------------------------------------------------------------*/

	private static String roleRow(final String artifactId, final String familyId, final String role,
			final int foldId){
		return artifactId+"\t"+familyId+"\t"+role+"\t"+foldId+"\n";
	}

	/** Builds FROZEN_K=20 distinct placeholder neighbor family names for a given prefix. */
	private static String[] frozenKNeighbors(final String prefix){
		final String[] names=new String[FROZEN_K];
		for(int i=0; i<FROZEN_K; i++){names[i]="nbr"+prefix+i;}
		return names;
	}

	/** Emits every neighbor family's {@code all_member_decoy}/{@code fold_id=-1} role row -- the proof
	 * of existence Elly's finding requires before a neighbor's decoy artifact_id may be used anywhere. */
	private static String decoyRowsForFamilies(final String[] families){
		final StringBuilder sb=new StringBuilder();
		for(final String f : families){sb.append(roleRow(f+"__decoy", f, "all_member_decoy", -1));}
		return sb.toString();
	}

	/** Emits one focal's full K=20 neighbor block. If {@code scramble} is true, rows are written in a
	 * REVERSED (not rank-ascending) file order -- proving {@link #loadNeighbors} sorts by rank itself
	 * rather than trusting file order. */
	private static String neighborBlock(final String focal, final String[] neighbors, final boolean scramble){
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<neighbors.length; i++){
			final int idx=scramble ? (neighbors.length-1-i) : i;
			sb.append(focal).append('\t').append(idx+1).append('\t').append(neighbors[idx]).append('\n');
		}
		return sb.toString();
	}

	/** Builds a synthetic but structurally-valid HMMER ASCII multi-profile concatenation: one minimal
	 * profile block per artifact_id, {@code NAME} lines carrying that id's real artifact_key, in the
	 * EXACT given order -- what {@link #parseHmmdbNamesInOrder} is meant to parse back out correctly. */
	private static String syntheticDbHmmdbContent(final ArrayList<String> orderedIds){
		final ArrayList<String> keys=new ArrayList<String>(orderedIds.size());
		for(final String id : orderedIds){keys.add(CandidateCUtil.artifactKey(id));}
		return syntheticDbHmmdbContentFromKeys(keys);
	}

	/** One minimal HMMER-ASCII profile block for a given NAME-line key -- the SAME literal bytes this
	 * class's synthetic db.hmmdb fixtures use per member, factored out so a per-artifact
	 * {@code .hmm.canonical} fixture file and the multi-profile {@code db.hmmdb} it concatenates into
	 * are built from the identical byte source (required for the sec 1 staleness-closure check to see a
	 * genuinely matching concatenation, not a coincidentally-equal one). */
	private static String singleHmmCanonicalBlockForKey(final String key){
		return "HMMER3/f [3.4 | Aug 2023]\nNAME  "+key+"\nLENG  10\n//\n";
	}

	/** Same as {@link #syntheticDbHmmdbContent} but takes raw NAME-line keys directly -- used by
	 * adversarial fixtures that need a key NOT derivable from any real artifact_id (a foreign/wrong
	 * member) or a deliberately reordered key sequence (a pure permutation). */
	private static String syntheticDbHmmdbContentFromKeys(final ArrayList<String> keys){
		final StringBuilder sb=new StringBuilder();
		for(final String key : keys){sb.append(singleHmmCanonicalBlockForKey(key));}
		return sb.toString();
	}

	/** Writes a real, hashable per-artifact {@code .hmm.canonical} fixture file for every id in
	 * {@code orderedIds} (content = {@link #singleHmmCanonicalBlockForKey} for that id's real
	 * artifact_key -- so their concatenation in order is byte-identical to
	 * {@link #syntheticDbHmmdbContent} over the same list) plus a matching
	 * {@code candidate_c_artifact_manifest.tsv} at {@code work+"/artifact_manifest.tsv"}, mirroring
	 * {@link CandidateCChallengeDbBuilder}'s own fixture schema exactly. Returns that manifest path. */
	private static String writeArtifactManifestFixture(final String work, final ArrayList<String> orderedIds)
			throws Exception{
		final String hmmDir=work+"/hmms";
		new java.io.File(hmmDir).mkdirs();
		final StringBuilder rows=new StringBuilder();
		for(final String id : orderedIds){
			final String key=CandidateCUtil.artifactKey(id);
			final String hmmPath=hmmDir+"/"+key+".hmm.canonical";
			writeFileForFixture(hmmPath, singleHmmCanonicalBlockForKey(key));
			final String hmmSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(hmmPath);
			rows.append(id).append('\t').append(key).append("\tfam\trole\t0\tpool\t1\tselhash\tfaahash\t")
				.append(hmmDir).append('/').append(key).append(".msa\tmsahash\t").append(hmmPath).append('\t')
				.append(hmmSha256).append('\n');
		}
		final String artifactManifestFile=work+"/artifact_manifest.tsv";
		writeFileForFixture(artifactManifestFile, "#artifact_id\tartifact_key\tfamily_id\trole\tfold_id\t"
			+"training_pool_sha256\tselected_member_count\tselected_member_sha256_ordered\tfaa_sha256\t"
			+"msa_path\tmsa_sha256\thmm_canonical_path\thmm_canonical_sha256\n"+rows);
		return artifactManifestFile;
	}

	/** A header-only artifact manifest (zero data rows) -- for fixtures whose target check fires strictly
	 * BEFORE this class's per-challenge artifact-manifest resolution ever runs (e.g. the two-way
	 * challenge-set-equality check, a malformed neighbors.tsv, a duplicate sidecar field): the file must
	 * still be loadable (this class loads it once, unconditionally, at the top of every run), but its
	 * content is never consulted for these tests. */
	private static String minimalArtifactManifestFixture(final String work) throws Exception{
		final String path=work+"/artifact_manifest.tsv";
		writeFileForFixture(path, "#artifact_id\tartifact_key\tfamily_id\trole\tfold_id\t"
			+"training_pool_sha256\tselected_member_count\tselected_member_sha256_ordered\tfaa_sha256\t"
			+"msa_path\tmsa_sha256\thmm_canonical_path\thmm_canonical_sha256\n");
		return path;
	}

	/** Writes a real {@code <challenge_key>.hmmdb.d/} directory: real db.hmmdb/.h3{f,i,m,p} bytes and a
	 * {@code complete.sidecar} whose recorded hashes genuinely match those files (never a hand-typed
	 * placeholder hash), so the aggregator's independent rehashing has real bytes to check against.
	 * {@code expectedSetHash}/{@code expectedOrderHash} are passed in explicitly so a test can supply a
	 * deliberately WRONG value to exercise the staleness checks. */
	private static void writeChallengeFixture(final String sidecarDir, final String challengeKey,
			final String focalFamilyId, final int foldId, final int profileCount,
			final String expectedSetHash, final String expectedOrderHash, final String dbHmmdbContent,
			final String h3fContent, final String h3iContent, final String h3mContent,
			final String h3pContent, final String hmmpressVersion) throws Exception{
		final String dbDir=sidecarDir+"/"+challengeKey+".hmmdb.d";
		new java.io.File(dbDir).mkdirs();
		final String dbHmmdbPath=dbDir+"/db.hmmdb";
		final String h3fPath=dbDir+"/db.hmmdb.h3f";
		final String h3iPath=dbDir+"/db.hmmdb.h3i";
		final String h3mPath=dbDir+"/db.hmmdb.h3m";
		final String h3pPath=dbDir+"/db.hmmdb.h3p";
		writeFileForFixture(dbHmmdbPath, dbHmmdbContent);
		writeFileForFixture(h3fPath, h3fContent);
		writeFileForFixture(h3iPath, h3iContent);
		writeFileForFixture(h3mPath, h3mContent);
		writeFileForFixture(h3pPath, h3pContent);
		final String dbHmmdbSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(dbHmmdbPath);
		final String h3fSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(h3fPath);
		final String h3iSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(h3iPath);
		final String h3mSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(h3mPath);
		final String h3pSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(h3pPath);
		writeFileForFixture(dbDir+"/complete.sidecar",
			"focal_family_id\t"+focalFamilyId+"\n"
			+"fold_id\t"+foldId+"\n"
			+"challenge_key\t"+challengeKey+"\n"
			+"profile_count\t"+profileCount+"\n"
			+"expected_artifact_id_set_sha256\t"+expectedSetHash+"\n"
			+"expected_profile_order_sha256\t"+expectedOrderHash+"\n"
			+"concat_canonical_sha256\t"+dbHmmdbSha256+"\n"
			+"db_hmmdb_sha256\t"+dbHmmdbSha256+"\n"
			+"h3f_sha256\t"+h3fSha256+"\n"
			+"h3i_sha256\t"+h3iSha256+"\n"
			+"h3m_sha256\t"+h3mSha256+"\n"
			+"h3p_sha256\t"+h3pSha256+"\n"
			+"hmmpress_version\t"+hmmpressVersion+"\n");
	}

	/** Builds a complete, VALID one-focal/one-fold fixture set (role.tsv WITH every neighbor's decoy
	 * row, neighbors.tsv, and a matching real challenge directory whose db.hmmdb carries real,
	 * correctly-ordered NAME lines) and returns the challenge_key -- shared setup for most tests. */
	private static String buildValidSingleChallengeFixture(final String work, final String artifactId,
			final String familyId, final int foldId, final boolean scrambleNeighbors) throws Exception{
		new java.io.File(work).mkdirs();
		final String[] neighbors=frozenKNeighbors(familyId);
		writeFileForFixture(work+"/role.tsv", "#artifact_id\tfamily_id\trole\tfold_id\n"
			+roleRow(artifactId, familyId, "focal_fold_local", foldId)
			+decoyRowsForFamilies(neighbors));
		writeFileForFixture(work+"/neighbors.tsv", "#k\t"+FROZEN_K+"\n"
			+neighborBlock(familyId, neighbors, scrambleNeighbors));
		final ArrayList<String> orderedProfileIds=new ArrayList<String>(FROZEN_K+1);
		orderedProfileIds.add(artifactId);
		for(final String n : neighbors){orderedProfileIds.add(n+"__decoy");}
		final String expectedSetHash=ArtifactRoleManifestGenerator.canonicalListSha256(orderedProfileIds);
		final String expectedOrderHash=CandidateCUtil.orderedIdListSha256(orderedProfileIds);
		final String artifactKey=CandidateCUtil.artifactKey(artifactId);
		final String challengeKey=CandidateCUtil.challengeKey(artifactKey);
		final String sidecarDir=work+"/sidecars";
		new java.io.File(sidecarDir).mkdirs();
		writeChallengeFixture(sidecarDir, challengeKey, familyId, foldId, orderedProfileIds.size(),
			expectedSetHash, expectedOrderHash, syntheticDbHmmdbContent(orderedProfileIds),
			"H3F_"+challengeKey, "H3I_"+challengeKey, "H3M_"+challengeKey, "H3P_"+challengeKey,
			"HMMER 3.4 (Aug 2023)");
		writeArtifactManifestFixture(work, orderedProfileIds);
		return challengeKey;
	}

	/** Builds a one-focal/one-fold fixture whose {@code db.hmmdb} carries EXACTLY the given
	 * {@code actualDbNames} keys (which may deliberately differ in membership/order from the CORRECT
	 * expected 21), while every sidecar field is otherwise self-consistent: {@code expected_artifact_
	 * id_set_sha256}/{@code expected_profile_order_sha256} reflect the TRUE upstream expectation (so
	 * the earlier staleness checks pass), and {@code concat_canonical_sha256}/{@code db_hmmdb_sha256}
	 * are the REAL hash of the (possibly wrong) db.hmmdb content actually written (so the byte-rehash
	 * checks pass too) -- isolating the semantic NAME-line content check as the ONLY thing that can
	 * still catch the problem. */
	private static void buildFixtureWithDbHmmdbNames(final String work, final String artifactId,
			final String familyId, final int foldId, final ArrayList<String> actualDbNames) throws Exception{
		new java.io.File(work).mkdirs();
		final String[] neighbors=frozenKNeighbors(familyId);
		writeFileForFixture(work+"/role.tsv", "#artifact_id\tfamily_id\trole\tfold_id\n"
			+roleRow(artifactId, familyId, "focal_fold_local", foldId)
			+decoyRowsForFamilies(neighbors));
		writeFileForFixture(work+"/neighbors.tsv", "#k\t"+FROZEN_K+"\n"+neighborBlock(familyId, neighbors, false));
		final ArrayList<String> correctOrderedIds=new ArrayList<String>(FROZEN_K+1);
		correctOrderedIds.add(artifactId);
		for(final String n : neighbors){correctOrderedIds.add(n+"__decoy");}
		final String expectedSetHash=ArtifactRoleManifestGenerator.canonicalListSha256(correctOrderedIds);
		final String expectedOrderHash=CandidateCUtil.orderedIdListSha256(correctOrderedIds);
		final String challengeKey=CandidateCUtil.challengeKey(CandidateCUtil.artifactKey(artifactId));
		final String dbDir=work+"/sidecars/"+challengeKey+".hmmdb.d";
		new java.io.File(dbDir).mkdirs();
		final String dbHmmdbPath=dbDir+"/db.hmmdb";
		writeFileForFixture(dbHmmdbPath, syntheticDbHmmdbContentFromKeys(actualDbNames));
		final String h3fPath=dbDir+"/db.hmmdb.h3f", h3iPath=dbDir+"/db.hmmdb.h3i",
			h3mPath=dbDir+"/db.hmmdb.h3m", h3pPath=dbDir+"/db.hmmdb.h3p";
		writeFileForFixture(h3fPath, "H3F"); writeFileForFixture(h3iPath, "H3I");
		writeFileForFixture(h3mPath, "H3M"); writeFileForFixture(h3pPath, "H3P");
		final String dbHash=ArtifactRoleManifestGenerator.streamingSha256Hex(dbHmmdbPath);
		writeFileForFixture(dbDir+"/complete.sidecar",
			"focal_family_id\t"+familyId+"\n"
			+"fold_id\t"+foldId+"\n"
			+"challenge_key\t"+challengeKey+"\n"
			+"profile_count\t"+correctOrderedIds.size()+"\n"
			+"expected_artifact_id_set_sha256\t"+expectedSetHash+"\n"
			+"expected_profile_order_sha256\t"+expectedOrderHash+"\n"
			+"concat_canonical_sha256\t"+dbHash+"\n"
			+"db_hmmdb_sha256\t"+dbHash+"\n"
			+"h3f_sha256\t"+ArtifactRoleManifestGenerator.streamingSha256Hex(h3fPath)+"\n"
			+"h3i_sha256\t"+ArtifactRoleManifestGenerator.streamingSha256Hex(h3iPath)+"\n"
			+"h3m_sha256\t"+ArtifactRoleManifestGenerator.streamingSha256Hex(h3mPath)+"\n"
			+"h3p_sha256\t"+ArtifactRoleManifestGenerator.streamingSha256Hex(h3pPath)+"\n"
			+"hmmpress_version\tHMMER 3.4 (Aug 2023)\n");
		//The artifact-manifest resolution/staleness closure runs AFTER the NAME-line semantic check this
		//fixture is designed to fail at, so it is never reached here -- built correctly anyway (matching
		//correctOrderedIds, not the deliberately-wrong actualDbNames) purely for fixture self-consistency.
		writeArtifactManifestFixture(work, correctOrderedIds);
	}

	static void selftest() throws Exception{
		final String dir=System.getProperty("java.io.tmpdir")+"/candc_b2_selftest_"+System.nanoTime();
		new java.io.File(dir).mkdirs();

		testBasicAggregation(dir);
		System.err.println("  selftest[basic aggregation: full K=20 challenge set, rehash, atomic "
			+"publish + sha sidecar]: PASS");

		testMissingChallengeDirCrashes(dir);
		System.err.println("  selftest[missing challenge_key directory crashes loud]: PASS");

		testExtraChallengeDirCrashes(dir);
		System.err.println("  selftest[extra/orphaned challenge_key directory crashes loud]: PASS");

		testStaleExpectedSetHashCrashes(dir);
		System.err.println("  selftest[stale expected_artifact_id_set_sha256 crashes loud]: PASS");

		testStaleExpectedOrderHashCrashes(dir);
		System.err.println("  selftest[stale expected_profile_order_sha256 crashes loud]: PASS");

		testCorruptedDbHmmdbOnDiskCrashes(dir);
		System.err.println("  selftest[on-disk db.hmmdb corruption vs sidecar hash crashes loud]: PASS");

		testCorruptedH3FileOnDiskCrashes(dir);
		System.err.println("  selftest[on-disk db.hmmdb.h3m corruption vs sidecar hash crashes loud]: PASS");

		testDuplicateSidecarFieldCrashes(dir);
		System.err.println("  selftest[duplicate field in a complete.sidecar crashes loud]: PASS");

		testNonFrozenKCrashes(dir);
		System.err.println("  selftest[family_neighbors.tsv declaring K!=20 crashes loud]: PASS");

		testDuplicateFamilyFoldCrashes(dir);
		System.err.println("  selftest[duplicate (family_id, fold_id) sort key crashes loud]: PASS");

		testScrambledNeighborOrderStillCorrect(dir);
		System.err.println("  selftest[neighbors.tsv written in reversed file order still sorts by "
			+"neighbor_rank for the order-sensitive hash]: PASS");

		testDeterministicRowOrderMultiFamily(dir);
		System.err.println("  selftest[deterministic row order across families/folds under "
			+"out-of-order input]: PASS");

		testConcatCanonicalFieldAloneCorruptedCrashes(dir);
		System.err.println("  selftest[concat_canonical_sha256 alone wrong (db_hmmdb_sha256 correct) "
			+"crashes loud]: PASS");

		testDbHmmdbFieldAloneCorruptedCrashes(dir);
		System.err.println("  selftest[db_hmmdb_sha256 alone wrong (concat_canonical_sha256 correct) "
			+"crashes loud]: PASS");

		testNeighborDecoyMissingCrashes(dir);
		System.err.println("  selftest[a neighbor's required decoy artifact absent from the role "
			+"manifest crashes loud]: PASS");

		testNeighborDecoyWrongRoleCrashes(dir);
		System.err.println("  selftest[a neighbor's decoy artifact with the wrong fold_id crashes loud]: PASS");

		testNeighborDecoyWrongFamilyCrashes(dir);
		System.err.println("  selftest[a neighbor's decoy artifact claiming the wrong family_id crashes "
			+"loud]: PASS");

		testWrongMemberInDbHmmdbCrashes(dir);
		System.err.println("  selftest[db.hmmdb containing a wrong/unexpected member (sidecar hashes "
			+"otherwise self-consistent) crashes loud]: PASS");

		testPurePermutationInDbHmmdbCrashes(dir);
		System.err.println("  selftest[db.hmmdb with the correct 21 members in the WRONG physical order "
			+"crashes loud]: PASS");

		testEmptyH3FileCrashes(dir);
		System.err.println("  selftest[an empty h3 file (hash self-consistent) crashes loud on the "
			+"explicit non-empty check]: PASS");

		testDuplicateArtifactManifestRowCrashes(dir);
		System.err.println("  selftest[sec 6: a duplicate artifact-manifest artifact_id crashes loud]: PASS");

		testArtifactManifestRowMissingForChallengeMemberCrashes(dir);
		System.err.println("  selftest[sec 6: a challenge member with no artifact-manifest row crashes "
			+"loud]: PASS");

		testArtifactManifestRowHashMismatchCrashes(dir);
		System.err.println("  selftest[sec 6: an artifact-manifest row's declared hash vs live HMM "
			+"mismatch crashes loud]: PASS");

		testStaleHmmBytesWithUnrebuiltPhaseBCrashes(dir);
		System.err.println("  selftest[sec 6 FLAGSHIP: current canonical-HMM byte change with identical "
			+"IDs/order and a stale (unrebuilt) Phase-B sidecar crashes loud]: PASS");

		System.err.println("CandidateCChallengeDbAggregator selftest: ALL PASS.");
	}

	private static void testBasicAggregation(final String dir) throws Exception{
		final String work=dir+"/basic";
		final String challengeKey=buildValidSingleChallengeFixture(work, "fZ__focal_f0", "fZ", 0, false);
		final String outFile=work+"/challenge_db_manifest.tsv";
		process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", work+"/sidecars", outFile);
		final String manifest=readWhole(outFile);
		if(!manifest.contains("fZ") || !manifest.contains(challengeKey) || !manifest.contains("21")){
			throw new RuntimeException("SELFTEST FAILED: manifest missing expected fields: "+manifest);
		}
		if(!manifest.contains("#schema_version\t4")){
			throw new RuntimeException("SELFTEST FAILED: manifest missing schema_version 4 header: "+manifest);
		}
		//Sec 6's own required check: the published header's path/hash must equal the EXACT current
		//artifact-manifest bytes -- rehashed independently here, never trusted from what this class
		//computed internally.
		final String artifactManifestFile=work+"/artifact_manifest.tsv";
		final String actualArtifactManifestSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(artifactManifestFile);
		if(!manifest.contains("#candidate_c_artifact_manifest_path\t"+artifactManifestFile)){
			throw new RuntimeException("SELFTEST FAILED: manifest missing the "
				+"#candidate_c_artifact_manifest_path header: "+manifest);
		}
		if(!manifest.contains("#candidate_c_artifact_manifest_sha256\t"+actualArtifactManifestSha256)){
			throw new RuntimeException("SELFTEST FAILED: published #candidate_c_artifact_manifest_sha256 "
				+"does not equal the independently-rehashed CURRENT artifact-manifest bytes ('"
				+actualArtifactManifestSha256+"'): "+manifest);
		}
		final byte[] published=Files.readAllBytes(Paths.get(outFile));
		final String expectedHash=ArtifactRoleManifestGenerator.sha256Hex(published);
		final String sidecarText=readWhole(outFile+".sha256");
		if(!sidecarText.contains(expectedHash)){
			throw new RuntimeException("SELFTEST FAILED: "+outFile+".sha256 does not match the actual "
				+"published manifest bytes: "+sidecarText);
		}
	}

	private static void testMissingChallengeDirCrashes(final String dir) throws Exception{
		final String work=dir+"/missing";
		new java.io.File(work).mkdirs();
		final String[] neighbors=frozenKNeighbors("fZ");
		writeFileForFixture(work+"/role.tsv", "#artifact_id\tfamily_id\trole\tfold_id\n"
			+roleRow("fZ__focal_f0", "fZ", "focal_fold_local", 0)+decoyRowsForFamilies(neighbors));
		writeFileForFixture(work+"/neighbors.tsv", "#k\t"+FROZEN_K+"\n"+neighborBlock("fZ", neighbors, false));
		minimalArtifactManifestFixture(work);//crash fires before this class ever consults it
		new java.io.File(work+"/sidecars").mkdirs();//no challenge dir written at all
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", work+"/sidecars", work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("do not exactly match")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a missing challenge "
					+"dir: "+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: a missing challenge dir did not crash.");}
	}

	private static void testExtraChallengeDirCrashes(final String dir) throws Exception{
		final String work=dir+"/extra";
		buildValidSingleChallengeFixture(work, "fZ__focal_f0", "fZ", 0, false);
		//An ORPHANED challenge dir from a superseded run.
		writeChallengeFixture(work+"/sidecars", "db_p_orphanorphanorphan", "stale", 0, 1,
			"stalesethash", "staleorderhash", "STALE", "STALE", "STALE", "STALE", "STALE",
			"HMMER 3.4 (Aug 2023)");
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", work+"/sidecars", work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("do not exactly match")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for an extra challenge "
					+"dir: "+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: an extra/orphaned challenge dir did not crash.");}
	}

	private static void testStaleExpectedSetHashCrashes(final String dir) throws Exception{
		final String work=dir+"/stalesethash";
		new java.io.File(work).mkdirs();
		final String[] neighbors=frozenKNeighbors("fZ");
		writeFileForFixture(work+"/role.tsv", "#artifact_id\tfamily_id\trole\tfold_id\n"
			+roleRow("fZ__focal_f0", "fZ", "focal_fold_local", 0)+decoyRowsForFamilies(neighbors));
		writeFileForFixture(work+"/neighbors.tsv", "#k\t"+FROZEN_K+"\n"+neighborBlock("fZ", neighbors, false));
		minimalArtifactManifestFixture(work);//crash fires before this class ever consults it
		final String artifactKey=CandidateCUtil.artifactKey("fZ__focal_f0");
		final String challengeKey=CandidateCUtil.challengeKey(artifactKey);
		final String sidecarDir=work+"/sidecars";
		//A WRONG expected_artifact_id_set_sha256 -- as if the neighbor list changed after this sidecar
		//was written but Phase B was never rerun.
		writeChallengeFixture(sidecarDir, challengeKey, "fZ", 0, FROZEN_K+1,
			"STALE_SET_HASH_FROM_BEFORE_NEIGHBOR_LIST_CHANGED", "irrelevant_order_hash_never_checked_first",
			"HMMDB", "H3F", "H3I", "H3M", "H3P", "HMMER 3.4 (Aug 2023)");
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", sidecarDir, work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("expected_artifact_id_set_sha256")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a stale expected-set "
					+"hash: "+e.getMessage());
			}
		}
		if(!crashed){
			throw new RuntimeException("SELFTEST FAILED: a stale expected_artifact_id_set_sha256 did not crash.");
		}
	}

	private static void testStaleExpectedOrderHashCrashes(final String dir) throws Exception{
		final String work=dir+"/staleorderhash";
		new java.io.File(work).mkdirs();
		final String[] neighbors=frozenKNeighbors("fZ");
		writeFileForFixture(work+"/role.tsv", "#artifact_id\tfamily_id\trole\tfold_id\n"
			+roleRow("fZ__focal_f0", "fZ", "focal_fold_local", 0)+decoyRowsForFamilies(neighbors));
		writeFileForFixture(work+"/neighbors.tsv", "#k\t"+FROZEN_K+"\n"+neighborBlock("fZ", neighbors, false));
		final ArrayList<String> orderedProfileIds=new ArrayList<String>(FROZEN_K+1);
		orderedProfileIds.add("fZ__focal_f0");
		for(final String n : neighbors){orderedProfileIds.add(n+"__decoy");}
		minimalArtifactManifestFixture(work);//crash fires before this class ever consults it
		//The SET hash is correct (same 21 members)...
		final String correctSetHash=ArtifactRoleManifestGenerator.canonicalListSha256(orderedProfileIds);
		final String artifactKey=CandidateCUtil.artifactKey("fZ__focal_f0");
		final String challengeKey=CandidateCUtil.challengeKey(artifactKey);
		final String sidecarDir=work+"/sidecars";
		//...but the ORDER hash is wrong -- e.g. the neighbor_rank assignment changed (a permutation),
		//which the set-only check would miss.
		writeChallengeFixture(sidecarDir, challengeKey, "fZ", 0, orderedProfileIds.size(),
			correctSetHash, "STALE_ORDER_HASH_FROM_A_DIFFERENT_PERMUTATION",
			"HMMDB", "H3F", "H3I", "H3M", "H3P", "HMMER 3.4 (Aug 2023)");
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", sidecarDir, work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("expected_profile_order_sha256")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a stale expected-order "
					+"hash: "+e.getMessage());
			}
		}
		if(!crashed){
			throw new RuntimeException("SELFTEST FAILED: a stale expected_profile_order_sha256 (set "
				+"correct, order wrong) did not crash.");
		}
	}

	private static void testCorruptedDbHmmdbOnDiskCrashes(final String dir) throws Exception{
		final String work=dir+"/corruptdb";
		final String challengeKey=buildValidSingleChallengeFixture(work, "fZ__focal_f0", "fZ", 0, false);
		//Corrupt db.hmmdb AFTER its sidecar was written -- the sidecar's recorded hash is now stale.
		writeFileForFixture(work+"/sidecars/"+challengeKey+".hmmdb.d/db.hmmdb", "CORRUPTED_AFTER_SIDECAR_WRITTEN");
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", work+"/sidecars", work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("corrupted or modified") || !e.getMessage().contains("db.hmmdb")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for corrupted on-disk "
					+"db.hmmdb: "+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: a corrupted on-disk db.hmmdb did not crash.");}
	}

	private static void testCorruptedH3FileOnDiskCrashes(final String dir) throws Exception{
		final String work=dir+"/corruptH3";
		final String challengeKey=buildValidSingleChallengeFixture(work, "fZ__focal_f0", "fZ", 0, false);
		writeFileForFixture(work+"/sidecars/"+challengeKey+".hmmdb.d/db.hmmdb.h3m", "CORRUPTED_AFTER_SIDECAR");
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", work+"/sidecars", work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("corrupted or modified") || !e.getMessage().contains("db.hmmdb.h3m")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for corrupted on-disk "
					+"db.hmmdb.h3m: "+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: a corrupted on-disk db.hmmdb.h3m did not crash.");}
	}

	private static void testDuplicateSidecarFieldCrashes(final String dir) throws Exception{
		final String work=dir+"/dupfield";
		new java.io.File(work).mkdirs();
		final String[] neighbors=frozenKNeighbors("fZ");
		writeFileForFixture(work+"/role.tsv", "#artifact_id\tfamily_id\trole\tfold_id\n"
			+roleRow("fZ__focal_f0", "fZ", "focal_fold_local", 0)+decoyRowsForFamilies(neighbors));
		writeFileForFixture(work+"/neighbors.tsv", "#k\t"+FROZEN_K+"\n"+neighborBlock("fZ", neighbors, false));
		final ArrayList<String> orderedProfileIds=new ArrayList<String>(FROZEN_K+1);
		orderedProfileIds.add("fZ__focal_f0");
		for(final String n : neighbors){orderedProfileIds.add(n+"__decoy");}
		minimalArtifactManifestFixture(work);//crash fires before this class ever consults it
		final String setHash=ArtifactRoleManifestGenerator.canonicalListSha256(orderedProfileIds);
		final String orderHash=CandidateCUtil.orderedIdListSha256(orderedProfileIds);
		final String artifactKey=CandidateCUtil.artifactKey("fZ__focal_f0");
		final String challengeKey=CandidateCUtil.challengeKey(artifactKey);
		final String dbDir=work+"/sidecars/"+challengeKey+".hmmdb.d";
		new java.io.File(dbDir).mkdirs();
		writeFileForFixture(dbDir+"/db.hmmdb", "HMMDB");
		writeFileForFixture(dbDir+"/db.hmmdb.h3f", "H3F");
		writeFileForFixture(dbDir+"/db.hmmdb.h3i", "H3I");
		writeFileForFixture(dbDir+"/db.hmmdb.h3m", "H3M");
		writeFileForFixture(dbDir+"/db.hmmdb.h3p", "H3P");
		final String dbHash=ArtifactRoleManifestGenerator.streamingSha256Hex(dbDir+"/db.hmmdb");
		//hmmpress_version recorded TWICE with different values.
		writeFileForFixture(dbDir+"/complete.sidecar",
			"focal_family_id\tfZ\n"
			+"fold_id\t0\n"
			+"challenge_key\t"+challengeKey+"\n"
			+"profile_count\t"+orderedProfileIds.size()+"\n"
			+"expected_artifact_id_set_sha256\t"+setHash+"\n"
			+"expected_profile_order_sha256\t"+orderHash+"\n"
			+"concat_canonical_sha256\t"+dbHash+"\n"
			+"db_hmmdb_sha256\t"+dbHash+"\n"
			+"h3f_sha256\t"+ArtifactRoleManifestGenerator.streamingSha256Hex(dbDir+"/db.hmmdb.h3f")+"\n"
			+"h3i_sha256\t"+ArtifactRoleManifestGenerator.streamingSha256Hex(dbDir+"/db.hmmdb.h3i")+"\n"
			+"h3m_sha256\t"+ArtifactRoleManifestGenerator.streamingSha256Hex(dbDir+"/db.hmmdb.h3m")+"\n"
			+"h3p_sha256\t"+ArtifactRoleManifestGenerator.streamingSha256Hex(dbDir+"/db.hmmdb.h3p")+"\n"
			+"hmmpress_version\tHMMER 3.4 (Aug 2023)\n"
			+"hmmpress_version\tSOME_OTHER_VERSION\n");
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", work+"/sidecars", work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("duplicate field")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a duplicate sidecar "
					+"field: "+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: a duplicate sidecar field did not crash.");}
	}

	private static void testNonFrozenKCrashes(final String dir) throws Exception{
		final String work=dir+"/nonfrozenk";
		new java.io.File(work).mkdirs();
		writeFileForFixture(work+"/role.tsv", "#artifact_id\tfamily_id\trole\tfold_id\n"
			+roleRow("fZ__focal_f0", "fZ", "focal_fold_local", 0));
		//Declares K=19, not the frozen 20.
		writeFileForFixture(work+"/neighbors.tsv", "#k\t19\n"
			+neighborBlock("fZ", java.util.Arrays.copyOf(frozenKNeighbors("fZ"), 19), false));
		minimalArtifactManifestFixture(work);//crash fires before this class ever consults it
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", work+"/sidecars", work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("FROZEN")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a non-frozen K: "
					+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: a non-frozen K did not crash.");}
	}

	private static void testDuplicateFamilyFoldCrashes(final String dir) throws Exception{
		final String work=dir+"/dupfamilyfold";
		new java.io.File(work).mkdirs();
		//Two DIFFERENT artifact_ids both claiming (family_id="fZ", fold_id=0) -- malformed input.
		writeFileForFixture(work+"/role.tsv", "#artifact_id\tfamily_id\trole\tfold_id\n"
			+roleRow("fZ__focal_f0", "fZ", "focal_fold_local", 0)
			+roleRow("fZ__focal_f0_dup", "fZ", "focal_fold_local", 0));
		final String[] neighbors=frozenKNeighbors("fZ");
		writeFileForFixture(work+"/neighbors.tsv", "#k\t"+FROZEN_K+"\n"+neighborBlock("fZ", neighbors, false));
		minimalArtifactManifestFixture(work);//crash fires before this class ever consults it
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", work+"/sidecars", work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("Duplicate (family_id, fold_id)")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a duplicate "
					+"(family_id, fold_id): "+e.getMessage());
			}
		}
		if(!crashed){
			throw new RuntimeException("SELFTEST FAILED: a duplicate (family_id, fold_id) did not crash.");
		}
	}

	/** Proves {@link #loadNeighbors} sorts by neighbor_rank itself: the fixture writes rows in
	 * REVERSED file order, yet the resulting expected_profile_order_sha256 must match what a
	 * rank-ascending construction produces -- a happy-path run (no crash) with the CORRECT order hash
	 * recorded in the sidecar. */
	private static void testScrambledNeighborOrderStillCorrect(final String dir) throws Exception{
		final String work=dir+"/scrambled";
		final String challengeKey=buildValidSingleChallengeFixture(work, "fZ__focal_f0", "fZ", 0, true);
		final String outFile=work+"/out.tsv";
		process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", work+"/sidecars", outFile);//must NOT throw
		final String manifest=readWhole(outFile);
		if(!manifest.contains(challengeKey)){
			throw new RuntimeException("SELFTEST FAILED: manifest missing the challenge_key: "+manifest);
		}
	}

	/** Multiple (family, fold) pairs across families, deliberately inserted out of order -- proves the
	 * published row order is the FROZEN (focal_family_id lexical, fold_id numeric ascending) key,
	 * matching {@link CandidateCArtifactAggregator}'s/{@link ArtifactRoleManifestGenerator}'s own
	 * convention, not insertion order. */
	private static void testDeterministicRowOrderMultiFamily(final String dir) throws Exception{
		final String work=dir+"/order";
		new java.io.File(work).mkdirs();
		//Deliberately out-of-order: fZ fold 1, then fA fold 0, then fZ fold 0.
		final String[] neighborsA=frozenKNeighbors("fA");
		final String[] neighborsZ=frozenKNeighbors("fZ");
		writeFileForFixture(work+"/role.tsv", "#artifact_id\tfamily_id\trole\tfold_id\n"
			+roleRow("fZ__focal_f1", "fZ", "focal_fold_local", 1)
			+roleRow("fA__focal_f0", "fA", "focal_fold_local", 0)
			+roleRow("fZ__focal_f0", "fZ", "focal_fold_local", 0)
			+decoyRowsForFamilies(neighborsA)+decoyRowsForFamilies(neighborsZ));
		writeFileForFixture(work+"/neighbors.tsv", "#k\t"+FROZEN_K+"\n"
			+neighborBlock("fA", neighborsA, false)+neighborBlock("fZ", neighborsZ, false));
		final String sidecarDir=work+"/sidecars";
		new java.io.File(sidecarDir).mkdirs();
		final String[][] rows={{"fZ__focal_f1", "fZ", "1"}, {"fA__focal_f0", "fA", "0"}, {"fZ__focal_f0", "fZ", "0"}};
		//The two fZ challenges (fold 0 and fold 1) share the SAME neighborsZ decoys -- dedupe before
		//writing one combined artifact manifest so the loader's own duplicate-artifact_id check doesn't
		//reject a decoy id that legitimately appears in more than one challenge's ordered list.
		final ArrayList<String> allIdsSeen=new ArrayList<String>();
		final HashSet<String> allIdsSeenSet=new HashSet<String>();
		for(final String[] row : rows){
			final String artifactId=row[0], familyId=row[1];
			final int foldId=Integer.parseInt(row[2]);
			final String[] neighbors=familyId.equals("fA") ? neighborsA : neighborsZ;
			final ArrayList<String> orderedProfileIds=new ArrayList<String>(FROZEN_K+1);
			orderedProfileIds.add(artifactId);
			for(final String n : neighbors){orderedProfileIds.add(n+"__decoy");}
			final String setHash=ArtifactRoleManifestGenerator.canonicalListSha256(orderedProfileIds);
			final String orderHash=CandidateCUtil.orderedIdListSha256(orderedProfileIds);
			final String challengeKey=CandidateCUtil.challengeKey(CandidateCUtil.artifactKey(artifactId));
			writeChallengeFixture(sidecarDir, challengeKey, familyId, foldId, orderedProfileIds.size(),
				setHash, orderHash, syntheticDbHmmdbContent(orderedProfileIds), "H3F_"+challengeKey,
				"H3I_"+challengeKey, "H3M_"+challengeKey, "H3P_"+challengeKey, "HMMER 3.4 (Aug 2023)");
			for(final String id : orderedProfileIds){
				if(allIdsSeenSet.add(id)){allIdsSeen.add(id);}
			}
		}
		writeArtifactManifestFixture(work, allIdsSeen);
		final String outFile=work+"/out.tsv";
		process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", sidecarDir, outFile);
		final String manifest=readWhole(outFile);
		final int idxA=manifest.indexOf("fA\t0\t");
		final int idxZ0=manifest.indexOf("fZ\t0\t");
		final int idxZ1=manifest.indexOf("fZ\t1\t");
		if(idxA<0 || idxZ0<0 || idxZ1<0){
			throw new RuntimeException("SELFTEST FAILED: manifest is missing an expected row: "+manifest);
		}
		if(!(idxA<idxZ0 && idxZ0<idxZ1)){
			throw new RuntimeException("SELFTEST FAILED: row order is not focal_family_id-lexical/"
				+"fold_id-ascending (expected fA:0, fZ:0, fZ:1 in that order): "+manifest);
		}
	}

	/** Elly's finding: {@code concat_canonical_sha256} and {@code db_hmmdb_sha256} are the SAME
	 * db.hmmdb file's hash but recorded under two provenance labels (pre-hmmpress vs. post-hmmpress
	 * rehash) -- both must be checked INDEPENDENTLY, not as a single OR'd condition. Here
	 * {@code db_hmmdb_sha256} is correct (matches the actual file) but {@code concat_canonical_sha256}
	 * alone is wrong -- must still crash on that specific field. */
	private static void testConcatCanonicalFieldAloneCorruptedCrashes(final String dir) throws Exception{
		final String work=dir+"/concatalone";
		buildFixtureWithOverride(work, "concat_canonical_sha256", "WRONG_CONCAT_HASH");
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", work+"/sidecars", work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("concat_canonical_sha256")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a lone corrupted "
					+"concat_canonical_sha256: "+e.getMessage());
			}
		}
		if(!crashed){
			throw new RuntimeException("SELFTEST FAILED: a lone corrupted concat_canonical_sha256 "
				+"(db_hmmdb_sha256 still correct) did not crash.");
		}
	}

	/** Mirror of the above: {@code concat_canonical_sha256} is correct, {@code db_hmmdb_sha256} alone
	 * is wrong -- must still crash on that specific field. */
	private static void testDbHmmdbFieldAloneCorruptedCrashes(final String dir) throws Exception{
		final String work=dir+"/dbhmmdbalone";
		buildFixtureWithOverride(work, "db_hmmdb_sha256", "WRONG_DB_HMMDB_HASH");
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", work+"/sidecars", work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("db_hmmdb_sha256")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a lone corrupted "
					+"db_hmmdb_sha256: "+e.getMessage());
			}
		}
		if(!crashed){
			throw new RuntimeException("SELFTEST FAILED: a lone corrupted db_hmmdb_sha256 "
				+"(concat_canonical_sha256 still correct) did not crash.");
		}
	}

	/** Shared setup for the two field-alone-corrupted tests: a valid one-focal/one-fold fixture whose
	 * {@code complete.sidecar} has ONE named field overridden to a wrong value, while every other field
	 * (including the OTHER db.hmmdb-hash field) stays correct against the actual on-disk bytes. */
	private static void buildFixtureWithOverride(final String work, final String overrideField,
			final String overrideValue) throws Exception{
		new java.io.File(work).mkdirs();
		final String[] neighbors=frozenKNeighbors("fZ");
		writeFileForFixture(work+"/role.tsv", "#artifact_id\tfamily_id\trole\tfold_id\n"
			+roleRow("fZ__focal_f0", "fZ", "focal_fold_local", 0)+decoyRowsForFamilies(neighbors));
		writeFileForFixture(work+"/neighbors.tsv", "#k\t"+FROZEN_K+"\n"+neighborBlock("fZ", neighbors, false));
		final ArrayList<String> orderedProfileIds=new ArrayList<String>(FROZEN_K+1);
		orderedProfileIds.add("fZ__focal_f0");
		for(final String n : neighbors){orderedProfileIds.add(n+"__decoy");}
		final String setHash=ArtifactRoleManifestGenerator.canonicalListSha256(orderedProfileIds);
		final String orderHash=CandidateCUtil.orderedIdListSha256(orderedProfileIds);
		final String challengeKey=CandidateCUtil.challengeKey(CandidateCUtil.artifactKey("fZ__focal_f0"));
		final String dbDir=work+"/sidecars/"+challengeKey+".hmmdb.d";
		new java.io.File(dbDir).mkdirs();
		final String dbHmmdbPath=dbDir+"/db.hmmdb";
		writeFileForFixture(dbHmmdbPath, syntheticDbHmmdbContent(orderedProfileIds));
		writeFileForFixture(dbDir+"/db.hmmdb.h3f", "H3F_CONTENT");
		writeFileForFixture(dbDir+"/db.hmmdb.h3i", "H3I_CONTENT");
		writeFileForFixture(dbDir+"/db.hmmdb.h3m", "H3M_CONTENT");
		writeFileForFixture(dbDir+"/db.hmmdb.h3p", "H3P_CONTENT");
		final String dbHash=ArtifactRoleManifestGenerator.streamingSha256Hex(dbHmmdbPath);
		final HashMap<String,String> fields=new HashMap<String,String>();
		fields.put("focal_family_id", "fZ");
		fields.put("fold_id", "0");
		fields.put("challenge_key", challengeKey);
		fields.put("profile_count", String.valueOf(orderedProfileIds.size()));
		fields.put("expected_artifact_id_set_sha256", setHash);
		fields.put("expected_profile_order_sha256", orderHash);
		fields.put("concat_canonical_sha256", dbHash);
		fields.put("db_hmmdb_sha256", dbHash);
		fields.put("h3f_sha256", ArtifactRoleManifestGenerator.streamingSha256Hex(dbDir+"/db.hmmdb.h3f"));
		fields.put("h3i_sha256", ArtifactRoleManifestGenerator.streamingSha256Hex(dbDir+"/db.hmmdb.h3i"));
		fields.put("h3m_sha256", ArtifactRoleManifestGenerator.streamingSha256Hex(dbDir+"/db.hmmdb.h3m"));
		fields.put("h3p_sha256", ArtifactRoleManifestGenerator.streamingSha256Hex(dbDir+"/db.hmmdb.h3p"));
		fields.put("hmmpress_version", "HMMER 3.4 (Aug 2023)");
		fields.put(overrideField, overrideValue);//exactly ONE field wrong, everything else correct.
		final StringBuilder sidecarText=new StringBuilder();
		for(final java.util.Map.Entry<String,String> e : fields.entrySet()){
			sidecarText.append(e.getKey()).append('\t').append(e.getValue()).append('\n');
		}
		writeFileForFixture(dbDir+"/complete.sidecar", sidecarText.toString());
		//Unreached before the requireRehashMatch crash this fixture targets, but built correctly (and
		//from the SAME bytes db.hmmdb was written from) for fixture self-consistency.
		writeArtifactManifestFixture(work, orderedProfileIds);
	}

	/** Elly's finding: a neighbor family's decoy artifact must NEVER be assumed to exist just because
	 * its ID can be constructed -- entirely absent from the role manifest must crash loud. */
	private static void testNeighborDecoyMissingCrashes(final String dir) throws Exception{
		final String work=dir+"/decoymissing";
		new java.io.File(work).mkdirs();
		final String[] neighbors=frozenKNeighbors("fZ");
		//role.tsv has ONLY the focal row -- no decoy rows at all.
		writeFileForFixture(work+"/role.tsv", "#artifact_id\tfamily_id\trole\tfold_id\n"
			+roleRow("fZ__focal_f0", "fZ", "focal_fold_local", 0));
		writeFileForFixture(work+"/neighbors.tsv", "#k\t"+FROZEN_K+"\n"+neighborBlock("fZ", neighbors, false));
		final String challengeKey=CandidateCUtil.challengeKey(CandidateCUtil.artifactKey("fZ__focal_f0"));
		//Only the DIRECTORY needs to exist for the upfront two-way set-equality check to pass -- the
		//crash happens before any sidecar is ever read.
		new java.io.File(work+"/sidecars/"+challengeKey+".hmmdb.d").mkdirs();
		minimalArtifactManifestFixture(work);//crash fires before this class ever consults it
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", work+"/sidecars", work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("does not exist in")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a missing neighbor "
					+"decoy: "+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: a missing neighbor decoy did not crash.");}
	}

	/** A neighbor's decoy artifact exists and has the correct role, but the WRONG fold_id (0 instead of
	 * the required -1) -- present but malformed must crash loud too, not just absent. */
	private static void testNeighborDecoyWrongRoleCrashes(final String dir) throws Exception{
		final String work=dir+"/decoywrongrole";
		new java.io.File(work).mkdirs();
		final String[] neighbors=frozenKNeighbors("fZ");
		final StringBuilder roleTsv=new StringBuilder("#artifact_id\tfamily_id\trole\tfold_id\n");
		roleTsv.append(roleRow("fZ__focal_f0", "fZ", "focal_fold_local", 0));
		for(int i=0; i<neighbors.length; i++){
			final int foldId=(i==0) ? 0 : -1;//the first neighbor's decoy has fold_id=0, not -1.
			roleTsv.append(roleRow(neighbors[i]+"__decoy", neighbors[i], "all_member_decoy", foldId));
		}
		writeFileForFixture(work+"/role.tsv", roleTsv.toString());
		writeFileForFixture(work+"/neighbors.tsv", "#k\t"+FROZEN_K+"\n"+neighborBlock("fZ", neighbors, false));
		final String challengeKey=CandidateCUtil.challengeKey(CandidateCUtil.artifactKey("fZ__focal_f0"));
		new java.io.File(work+"/sidecars/"+challengeKey+".hmmdb.d").mkdirs();
		minimalArtifactManifestFixture(work);//crash fires before this class ever consults it
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", work+"/sidecars", work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("fold_id=0") || !e.getMessage().contains("expected -1")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a malformed neighbor "
					+"decoy fold_id: "+e.getMessage());
			}
		}
		if(!crashed){
			throw new RuntimeException("SELFTEST FAILED: a malformed neighbor decoy fold_id did not crash.");
		}
	}

	/** Elly's finding: a decoy row can carry the correct artifact_id/role/fold_id while claiming the
	 * WRONG family_id -- e.g. artifact_id="nbr0__decoy" but family_id="someOtherFamily". Present but
	 * internally inconsistent must crash loud, the same as absent or wrong-role/fold. */
	private static void testNeighborDecoyWrongFamilyCrashes(final String dir) throws Exception{
		final String work=dir+"/decoywrongfamily";
		new java.io.File(work).mkdirs();
		final String[] neighbors=frozenKNeighbors("fZ");
		final StringBuilder roleTsv=new StringBuilder("#artifact_id\tfamily_id\trole\tfold_id\n");
		roleTsv.append(roleRow("fZ__focal_f0", "fZ", "focal_fold_local", 0));
		for(int i=0; i<neighbors.length; i++){
			//The first neighbor's decoy row claims artifact_id=<neighbors[0]>__decoy but declares a
			//DIFFERENT family_id than the neighbor it's supposed to represent.
			final String familyId=(i==0) ? "someOtherFamily" : neighbors[i];
			roleTsv.append(roleRow(neighbors[i]+"__decoy", familyId, "all_member_decoy", -1));
		}
		writeFileForFixture(work+"/role.tsv", roleTsv.toString());
		writeFileForFixture(work+"/neighbors.tsv", "#k\t"+FROZEN_K+"\n"+neighborBlock("fZ", neighbors, false));
		final String challengeKey=CandidateCUtil.challengeKey(CandidateCUtil.artifactKey("fZ__focal_f0"));
		new java.io.File(work+"/sidecars/"+challengeKey+".hmmdb.d").mkdirs();
		minimalArtifactManifestFixture(work);//crash fires before this class ever consults it
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", work+"/sidecars", work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("family_id='someOtherFamily'") || !e.getMessage().contains("expected '"+neighbors[0]+"'")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a wrong-family neighbor "
					+"decoy: "+e.getMessage());
			}
		}
		if(!crashed){
			throw new RuntimeException("SELFTEST FAILED: a wrong-family neighbor decoy did not crash.");
		}
	}

	/** Elly's finding, the core semantic gap: db.hmmdb contains a WRONG/unexpected member (one of the
	 * 21 NAME lines is a foreign key, not one of the real focal/neighbor artifacts), while every sidecar
	 * hash field is self-consistent with the ACTUAL (wrong) file -- only the independent NAME-line
	 * re-derivation can catch this, since byte-rehashing alone sees nothing wrong. */
	private static void testWrongMemberInDbHmmdbCrashes(final String dir) throws Exception{
		final String work=dir+"/wrongmember";
		final String[] neighbors=frozenKNeighbors("fZ");
		final ArrayList<String> correctOrderedIds=new ArrayList<String>(FROZEN_K+1);
		correctOrderedIds.add("fZ__focal_f0");
		for(final String n : neighbors){correctOrderedIds.add(n+"__decoy");}
		final ArrayList<String> actualDbNames=new ArrayList<String>(correctOrderedIds.size());
		for(final String id : correctOrderedIds){actualDbNames.add(CandidateCUtil.artifactKey(id));}
		actualDbNames.set(5, "p_deadbeefcafebabe");//substitute one member for a foreign/unexpected key.
		buildFixtureWithDbHmmdbNames(work, "fZ__focal_f0", "fZ", 0, actualDbNames);
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", work+"/sidecars", work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("not one of this challenge set's")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a wrong member in "
					+"db.hmmdb: "+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: a wrong member in db.hmmdb did not crash.");}
	}

	/** Elly's finding, the permutation half: db.hmmdb contains the CORRECT 21 members (same set) but in
	 * the WRONG physical order -- the set-based check alone would miss this; only the order-hash
	 * re-derivation catches it. */
	private static void testPurePermutationInDbHmmdbCrashes(final String dir) throws Exception{
		final String work=dir+"/permutation";
		final String[] neighbors=frozenKNeighbors("fZ");
		final ArrayList<String> correctOrderedIds=new ArrayList<String>(FROZEN_K+1);
		correctOrderedIds.add("fZ__focal_f0");
		for(final String n : neighbors){correctOrderedIds.add(n+"__decoy");}
		final ArrayList<String> actualDbNames=new ArrayList<String>(correctOrderedIds.size());
		for(final String id : correctOrderedIds){actualDbNames.add(CandidateCUtil.artifactKey(id));}
		Collections.reverse(actualDbNames);//same 21 members, WRONG order (focal ends up last, not first).
		buildFixtureWithDbHmmdbNames(work, "fZ__focal_f0", "fZ", 0, actualDbNames);
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", work+"/sidecars", work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("ACTUAL physical NAME-line")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a pure permutation in "
					+"db.hmmdb: "+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: a pure permutation in db.hmmdb did not crash.");}
	}

	/** Elly's finding: an h3 file that is EMPTY (0 bytes), with the sidecar updated to record the REAL
	 * hash of empty content (so the byte-rehash check alone would consider it "correct") -- only the
	 * explicit non-empty check can catch this. */
	private static void testEmptyH3FileCrashes(final String dir) throws Exception{
		final String work=dir+"/emptyh3";
		final String challengeKey=buildValidSingleChallengeFixture(work, "fZ__focal_f0", "fZ", 0, false);
		final String dbDir=work+"/sidecars/"+challengeKey+".hmmdb.d";
		final String h3mPath=dbDir+"/db.hmmdb.h3m";
		writeFileForFixture(h3mPath, "");
		final String emptyHash=ArtifactRoleManifestGenerator.streamingSha256Hex(h3mPath);
		final String sidecarPath=dbDir+"/complete.sidecar";
		final String oldSidecar=readWhole(sidecarPath);
		final StringBuilder newSidecar=new StringBuilder();
		for(final String line : oldSidecar.split("\n")){
			if(line.isEmpty()){continue;}
			if(line.startsWith("h3m_sha256\t")){newSidecar.append("h3m_sha256\t").append(emptyHash).append('\n');}
			else{newSidecar.append(line).append('\n');}
		}
		writeFileForFixture(sidecarPath, newSidecar.toString());
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", work+"/artifact_manifest.tsv", work+"/sidecars", work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("must be non-empty")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for an empty h3 file: "
					+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: an empty h3 file did not crash.");}
	}

	// Sec 6 negative, NOT exercised as an in-process fixture: a missing artifactmanifest= file crashes
	// loud via fileIO.ReadWrite.getRawInputStream's existing
	// shared.KillSwitch.exceptionKill(new RuntimeException("Can't find file "+fname)) path -- confirmed
	// by direct observation (a real run against a deleted manifest file prints that exact stack trace and
	// halts the JVM). This HALTS the process rather than throwing a normally-catchable exception up
	// through process(), so it cannot be asserted with a try/catch inside this same JVM the way every
	// other negative here is -- the identical reason KillSwitch exists (crash-loud-never-hang) also makes
	// it untestable via a plain fixture. Verified by inspection, not a fixture negative.

	/** Sec 6 negative: two artifact-manifest rows sharing one artifact_id must crash loud -- mirrors
	 * {@link CandidateCChallengeDbBuilder#loadArtifactManifest}'s identical convention. */
	private static void testDuplicateArtifactManifestRowCrashes(final String dir) throws Exception{
		final String work=dir+"/dupmanifestrow";
		buildValidSingleChallengeFixture(work, "fZ__focal_f0", "fZ", 0, false);
		final String artifactManifestFile=work+"/artifact_manifest.tsv";
		final String oldText=readWhole(artifactManifestFile);
		String firstDataRow=null;
		for(final String line : oldText.split("\n")){
			if(!line.isEmpty() && line.charAt(0)!='#'){firstDataRow=line; break;}
		}
		if(firstDataRow==null){
			throw new RuntimeException("SELFTEST SETUP FAILED: could not find a data row to duplicate.");
		}
		writeFileForFixture(artifactManifestFile, oldText+"\n"+firstDataRow+"\n");
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", artifactManifestFile, work+"/sidecars",
				work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("Duplicate artifact_id")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a duplicate "
					+"artifact-manifest row: "+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: a duplicate artifact-manifest row did not crash.");}
	}

	/** Sec 6 negative: a challenge member with NO row at all in the current artifact manifest (Phase A2
	 * has not produced it yet) must crash loud, distinct from a present-but-stale row. */
	private static void testArtifactManifestRowMissingForChallengeMemberCrashes(final String dir) throws Exception{
		final String work=dir+"/missingmemberrow";
		buildValidSingleChallengeFixture(work, "fZ__focal_f0", "fZ", 0, false);
		final String artifactManifestFile=work+"/artifact_manifest.tsv";
		final String oldText=readWhole(artifactManifestFile);
		final String targetId="nbrfZ5__decoy";
		final StringBuilder newText=new StringBuilder();
		for(final String line : oldText.split("\n")){
			if(line.startsWith(targetId+"\t")){continue;}
			if(!line.isEmpty()){newText.append(line).append('\n');}
		}
		writeFileForFixture(artifactManifestFile, newText.toString());
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", artifactManifestFile, work+"/sidecars",
				work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("has no row in") || !e.getMessage().contains("Phase A2 has not produced it yet")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a challenge member "
					+"missing from the artifact manifest: "+e.getMessage());
			}
		}
		if(!crashed){
			throw new RuntimeException("SELFTEST FAILED: a challenge member missing from the artifact "
				+"manifest did not crash.");
		}
	}

	/** Sec 6 negative: an artifact-manifest row's declared {@code hmm_canonical_sha256} not matching the
	 * ACTUAL on-disk bytes (a stale/corrupt row, the live file itself untouched) must crash loud --
	 * distinct from the sec-1 closure test below, where the row IS updated to match changed live bytes. */
	private static void testArtifactManifestRowHashMismatchCrashes(final String dir) throws Exception{
		final String work=dir+"/rowhashmismatch";
		buildValidSingleChallengeFixture(work, "fZ__focal_f0", "fZ", 0, false);
		final String artifactManifestFile=work+"/artifact_manifest.tsv";
		final String oldText=readWhole(artifactManifestFile);
		String targetLine=null;
		for(final String line : oldText.split("\n")){
			if(line.startsWith("nbrfZ3__decoy\t")){targetLine=line; break;}
		}
		if(targetLine==null){throw new RuntimeException("SELFTEST SETUP FAILED: could not find the target row.");}
		final String[] fields=targetLine.split("\t");
		final String realHash=fields[fields.length-1];
		final String corruptedLine=targetLine.substring(0, targetLine.length()-realHash.length())
			+"WRONG_HASH_ABCDEF";
		writeFileForFixture(artifactManifestFile, oldText.replace(targetLine, corruptedLine));
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", artifactManifestFile, work+"/sidecars",
				work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("stale or corrupt artifact-manifest row")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a stale/corrupt "
					+"artifact-manifest row: "+e.getMessage());
			}
		}
		if(!crashed){
			throw new RuntimeException("SELFTEST FAILED: a stale/corrupt artifact-manifest row did not crash.");
		}
	}

	/** The FLAGSHIP sec 1 closure test (sec 6's own required negative): current canonical-HMM byte
	 * change with IDENTICAL IDs/order and a STALE Phase-B sidecar. One neighbor's live {@code
	 * .hmm.canonical} bytes change AND the artifact-manifest row is updated to match (simulating a real
	 * Phase A/A2 rerun -- the row is internally consistent, not stale-vs-live), but {@code db.hmmdb}/its
	 * sidecar are NEVER rebuilt. Every check that ran under the old schema-3 contract (set/order hashes,
	 * sidecar-vs-actual-db.hmmdb rehash, NAME-line semantic proof) still PASSES here -- only the new
	 * live-concat-vs-published-db.hmmdb comparison can catch it. */
	private static void testStaleHmmBytesWithUnrebuiltPhaseBCrashes(final String dir) throws Exception{
		final String work=dir+"/staleunrebuilt";
		buildValidSingleChallengeFixture(work, "fZ__focal_f0", "fZ", 0, false);
		final String hmmDir=work+"/hmms";
		final String changedKey=CandidateCUtil.artifactKey("nbrfZ7__decoy");
		final String changedHmmPath=hmmDir+"/"+changedKey+".hmm.canonical";
		writeFileForFixture(changedHmmPath, "HMMER3/f [3.4 | Aug 2023]\nNAME  "+changedKey+"\nLENG  99\n//\n");
		final String newHash=ArtifactRoleManifestGenerator.streamingSha256Hex(changedHmmPath);
		final String artifactManifestFile=work+"/artifact_manifest.tsv";
		final String oldText=readWhole(artifactManifestFile);
		String targetLine=null;
		for(final String line : oldText.split("\n")){
			if(line.startsWith("nbrfZ7__decoy\t")){targetLine=line; break;}
		}
		if(targetLine==null){throw new RuntimeException("SELFTEST SETUP FAILED: could not find the target row.");}
		final String[] fields=targetLine.split("\t");
		final String oldHash=fields[fields.length-1];
		final String updatedLine=targetLine.substring(0, targetLine.length()-oldHash.length())+newHash;
		writeFileForFixture(artifactManifestFile, oldText.replace(targetLine, updatedLine));
		boolean crashed=false;
		try{
			process(work+"/role.tsv", work+"/neighbors.tsv", artifactManifestFile, work+"/sidecars",
				work+"/out.tsv");
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("Phase A/A2 reran")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for the sec-1 staleness "
					+"closure: "+e.getMessage());
			}
		}
		if(!crashed){
			throw new RuntimeException("SELFTEST FAILED: current .hmm.canonical bytes changed (IDs/order "
				+"unchanged, manifest row updated to match) but Phase B never rerun did NOT crash -- the "
				+"sec 1 staleness closure did not fire.");
		}
	}
}
