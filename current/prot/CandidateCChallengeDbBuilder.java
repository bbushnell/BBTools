package prot;

import java.io.File;
import java.io.FileOutputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
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
 * Candidate C Phase B (`CANDIDATE_C_PHASE_AB_GATE_PLAN_v1.md`, UMP45/Elly co-sealed through
 * `3093447`): builds ONE focal-fold challenge database ({@code hmmpress}-ed {@code db.hmmdb}) from
 * the 21 current Candidate-C artifact-manifest ({@link CandidateCArtifactAggregator}/A2) rows its
 * focal+20-neighbor set names, through the constructor-injected {@link CommandRunner} seam.
 * <p>
 * <b>Current-HMM-byte staleness closure (sec 1 of the sealed plan -- the Phase-B analogue of A2's
 * current-FAA check):</b> {@code concat_canonical_sha256} is the SHA-256 of the exact raw-byte
 * concatenation of the 21 LIVE {@code .hmm.canonical} files, in pinned focal-then-neighbor-rank
 * order, with no inserted/removed/reconstructed bytes -- recomputed from the CURRENT
 * {@code candidate_c_artifact_manifest.tsv} rows before every resume decision. IDs/order staying
 * fixed does NOT prove the underlying HMM bytes are still current: Phase A can rebuild one or more
 * {@code .hmm.canonical} files (e.g. after a Phase-0 corpus edit) without changing any artifact ID,
 * and a stale Phase-B directory would otherwise still look internally consistent.
 * {@code db_hmmdb_sha256} is a SEPARATE, later-measured rehash of the same {@code db.hmmdb} path
 * taken AFTER {@code hmmpress} runs; requiring the two to be equal is the proof (not assumption)
 * that {@code hmmpress} does not mutate its input.
 * <p>
 * <b>Staging discipline (B5, Elly/UMP45 co-sealed {@code d2e1850}):</b> a rebuild creates a fresh,
 * UNIQUE {@code <sidecarDir>/<challenge_key>.hmmdb.building.<nonce>/} directory via
 * {@link java.nio.file.Files#createTempDirectory} as a same-filesystem SIBLING of the final
 * {@code <sidecarDir>/<challenge_key>.hmmdb.d/} directory -- never a cross-filesystem path, and never
 * a fixed reused name (a fixed name let overlapping retries for the same challenge_key delete/stomp on
 * each other's work; a genuinely unique path every call means there is never an existing path to
 * delete-then-reuse). The generated name can never equal or end in {@code .hmmdb.d}, so B2's top-level
 * {@code *.hmmdb.d} enumeration never miscounts it. The 21 live canonical-HMM bytes are streamed into
 * {@code db.hmmdb} there in pinned order while simultaneously computing the live concat hash (one
 * pass, no reconstruction), then independently re-hashed from the just-written file to prove the two
 * derivations agree; {@code hmmpress db.hmmdb} runs with that directory as {@code workingDirectory}
 * (never the shared {@code sidecarDir}, since the bare command has no path argument -- one shared
 * directory would make it ambiguous which task's {@code db.hmmdb} it presses). {@code complete.sidecar}
 * is written inside staging LAST, then the WHOLE staging directory is atomically renamed to the final
 * path (a directory-level rename, unlike Phase A's per-file rename, since {@code hmmpress}'s output is
 * 5 files that only mean anything together). An existing stale final directory is removed immediately
 * before that rename -- the resulting crash window (challenge briefly absent) is the INTENDED
 * fail-closed state; B2 rejects a missing challenge through its own exact-set-equality check, never
 * silently skipping it. Cleanup of staging or a stale final directory never follows symlinks (deleted
 * as itself, never traversed into) and fails loud if it does not fully complete (B4).
 * <p>
 * <b>Resume closure (B1/B3, {@code d2e1850}):</b> a valid resume requires the mandatory byte-freshness
 * equality chain {@code liveConcatHash == sc.concatCanonicalSha256 == sc.dbHmmdbSha256 ==
 * actualDbHmmdbSha256} (checked as direct pairwise comparisons, not merely transitively), plus the
 * sidecar's recorded {@code hmmpress_version} matching this class's pinned {@link #HMMPRESS_VERSION}.
 * <p>
 * <b>Artifact-key trust (B2, {@code d2e1850}):</b> every A2 artifact-manifest row's {@code
 * artifact_key} column is a COPIED value; before it is used to build the NAME-line translation map it
 * is validated against {@link CandidateCUtil#ARTIFACT_KEY_PATTERN} and required to equal the
 * independently-recomputed {@code CandidateCUtil.artifactKey(artifact_id)}, and the resulting map is
 * required to stay bijective (no two artifacts colliding on one key). Every live {@code .hmm.canonical}
 * input is required non-empty before it is hashed or concatenated.
 *
 * <p>Usage: {@code java -ea prot.CandidateCChallengeDbBuilder rolemanifest=artifact_role_manifest.tsv
 *        neighbors=family_neighbors.tsv artifactmanifest=candidate_c_artifact_manifest.tsv
 *        sidecardir=<shared dir also read by CandidateCChallengeDbAggregator/B2>
 *        focalfamilyid=fZ foldid=0}
 * <br>Self-test: {@code java -ea prot.CandidateCChallengeDbBuilder selftest}
 *
 * @author Eru
 */
public final class CandidateCChallengeDbBuilder {

	/** The frozen candidate-independent neighbor-challenge K -- not a tunable, ALWAYS 20.
	 * Independently redeclared per this project's convention of not sharing constants across tools. */
	static final int FROZEN_K=20;
	static final String HMMPRESS_VERSION="HMMER 3.4 (Aug 2023)";

	private final CommandRunner runner;

	CandidateCChallengeDbBuilder(final CommandRunner runner){
		this.runner=runner;
	}

	public static void main(final String[] args) throws Exception{
		if(args.length==1 && args[0].equalsIgnoreCase("selftest")){selftest(); return;}
		String roleFile=null, neighborsFile=null, artifactManifestFile=null, sidecarDir=null,
			focalFamilyId=null;
		Integer foldId=null;
		for(final String arg : args){
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase();
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("rolemanifest") || a.equals("role")){roleFile=b;}
			else if(a.equals("neighbors")){neighborsFile=b;}
			else if(a.equals("artifactmanifest")){artifactManifestFile=b;}
			else if(a.equals("sidecardir")){sidecarDir=b;}
			else if(a.equals("focalfamilyid")){focalFamilyId=b;}
			else if(a.equals("foldid")){foldId=Integer.parseInt(b);}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}
		if(roleFile==null || neighborsFile==null || artifactManifestFile==null || sidecarDir==null
				|| focalFamilyId==null || foldId==null){
			throw new RuntimeException("Required: rolemanifest= neighbors= artifactmanifest= "
				+"sidecardir= focalfamilyid= foldid=");
		}
		final Result r=new CandidateCChallengeDbBuilder(new RealCommandRunner())
			.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		System.err.println(r.summary);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Data holders        ----------------*/
	/*--------------------------------------------------------------*/

	static class UniverseRow { String artifactId, familyId, role; int foldId; }
	static class RankedNeighbor { int rank; String family; }
	static class ArtifactManifestRow { String artifactId, artifactKey, hmmCanonicalPath, hmmCanonicalSha256; }
	static class Sidecar {
		String focalFamilyId, challengeKey, expectedArtifactIdSetSha256, expectedProfileOrderSha256,
			concatCanonicalSha256, dbHmmdbSha256, h3fSha256, h3iSha256, h3mSha256, h3pSha256, hmmpressVersion;
		int foldId, profileCount;
	}
	static class Result { String summary; boolean resumed; }

	/** Small immutable-by-convention holder for one resolved challenge's exact 21-artifact
	 * selection -- the reusable output of {@link #resolveChallenge}, consumed by both the real
	 * build path here and any read-only preflight caller (e.g.
	 * {@code CandidateCChallengeResolverPreflight}) that needs the SAME selection without
	 * reimplementing it. */
	static class ResolvedChallenge {
		String challengeKey;
		ArrayList<String> orderedArtifactIds;
		ArrayList<String> orderedArtifactKeys;
		ArrayList<String> orderedHmmCanonicalPaths;
		ArrayList<String> orderedHmmCanonicalSha256s;
		HashMap<String,String> keyToId;
		String expectedSetHash;
		String expectedOrderHash;
		String liveConcatHash;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Pipeline           ----------------*/
	/*--------------------------------------------------------------*/

	/** Sec 5 steps 1-3: independently validates the full role universe, resolves the focal + 20
	 * neighbor-decoy artifacts in focal-then-neighbor-rank order, validates and rehashes every live
	 * {@code .hmm.canonical} input against its A2 artifact-manifest row, and computes the current
	 * expected set/order/concat hashes. Extracted verbatim (mechanical extraction, no logic change)
	 * from {@link #process} so a read-only preflight caller can reuse the EXACT SAME
	 * selection/validation the real build performs, rather than reimplementing it elsewhere (e.g.
	 * in shell). */
	static ResolvedChallenge resolveChallenge(final String roleFile, final String neighborsFile,
			final String artifactManifestFile, final String focalFamilyId, final int foldId) throws Exception{
		//Sec 5 step 1: independently validate the full role universe and the ONE requested focal row.
		final HashMap<String,UniverseRow> byArtifactId=loadRoleUniverse(roleFile);
		final UniverseRow focalRow=findFocalRow(byArtifactId, focalFamilyId, foldId, roleFile);
		final HashMap<String,ArrayList<String>> neighborsByFamily=loadNeighbors(neighborsFile);
		final ArrayList<String> neighborFamilies=neighborsByFamily.get(focalFamilyId);
		if(neighborFamilies==null){
			throw new RuntimeException("Focal family '"+focalFamilyId+"' has no row in "+neighborsFile+".");
		}

		final String artifactKey=CandidateCUtil.artifactKey(focalRow.artifactId);
		final String challengeKey=CandidateCUtil.challengeKey(artifactKey);

		final ArrayList<String> orderedProfileIds=new ArrayList<String>(neighborFamilies.size()+1);
		orderedProfileIds.add(focalRow.artifactId);
		for(final String neighborFamily : neighborFamilies){
			final String decoyId=neighborFamily+"__decoy";
			final UniverseRow decoyRow=byArtifactId.get(decoyId);
			if(decoyRow==null){
				throw new RuntimeException("Focal '"+focalRow.artifactId+"': neighbor family '"
					+neighborFamily+"'s required decoy artifact '"+decoyId+"' does not exist in "+roleFile+".");
			}
			if(!decoyRow.role.equals("all_member_decoy")){
				throw new RuntimeException("Focal '"+focalRow.artifactId+"': neighbor decoy artifact '"
					+decoyId+"' has role='"+decoyRow.role+"' in "+roleFile+", expected 'all_member_decoy'.");
			}
			if(!decoyRow.familyId.equals(neighborFamily)){
				throw new RuntimeException("Focal '"+focalRow.artifactId+"': neighbor decoy artifact '"
					+decoyId+"' has family_id='"+decoyRow.familyId+"' in "+roleFile+", expected '"
					+neighborFamily+"'.");
			}
			if(decoyRow.foldId!=-1){
				throw new RuntimeException("Focal '"+focalRow.artifactId+"': neighbor decoy artifact '"
					+decoyId+"' has fold_id="+decoyRow.foldId+" in "+roleFile+", expected -1.");
			}
			orderedProfileIds.add(decoyId);
		}

		//Sec 5 step 2: resolve the 21 current A2 rows; rehash every live .hmm.canonical against its
		//manifest row. Also builds the key->id map the NAME-line semantic check needs.
		//B2 finding (Elly, d2e1850): amr.artifactKey is a COPIED column from A2's own output -- it
		//must be validated against the safe pattern AND recomputed-equal to
		//CandidateCUtil.artifactKey(id) before it is trusted for anything (the NAME-line translation
		//map, or -- in a hypothetical future caller -- a filesystem path), and the map must stay
		//BIJECTIVE (no two distinct artifacts colliding on one key).
		final HashMap<String,ArtifactManifestRow> artifactManifest=loadArtifactManifest(artifactManifestFile);
		final ArrayList<String> hmmCanonicalPathsInOrder=new ArrayList<String>(orderedProfileIds.size());
		final ArrayList<String> orderedArtifactKeys=new ArrayList<String>(orderedProfileIds.size());
		final ArrayList<String> orderedHmmCanonicalSha256s=new ArrayList<String>(orderedProfileIds.size());
		final HashMap<String,String> keyToId=new HashMap<String,String>();
		for(final String id : orderedProfileIds){
			final ArtifactManifestRow amr=artifactManifest.get(id);
			if(amr==null){
				throw new RuntimeException("Challenge for focal '"+focalRow.artifactId+"': artifact '"
					+id+"' has no row in "+artifactManifestFile+" -- Phase A2 has not produced it yet.");
			}
			if(!CandidateCUtil.ARTIFACT_KEY_PATTERN.matcher(amr.artifactKey).matches()){
				throw new RuntimeException("Artifact '"+id+"': "+artifactManifestFile+" declares "
					+"artifact_key='"+amr.artifactKey+"' which does not match the required safe pattern "
					+CandidateCUtil.ARTIFACT_KEY_PATTERN+".");
			}
			final String recomputedKey=CandidateCUtil.artifactKey(id);
			if(!recomputedKey.equals(amr.artifactKey)){
				throw new RuntimeException("Artifact '"+id+"': "+artifactManifestFile+" declares "
					+"artifact_key='"+amr.artifactKey+"' but the independently-recomputed FNV key is '"
					+recomputedKey+"' -- refusing to trust an unverified copied key.");
			}
			if(!new File(amr.hmmCanonicalPath).exists() || new File(amr.hmmCanonicalPath).length()<=0){
				throw new RuntimeException("Artifact '"+id+"': the live .hmm.canonical at "
					+amr.hmmCanonicalPath+" declared by "+artifactManifestFile+" is missing or empty.");
			}
			final String actualHmmCanonicalSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(
				amr.hmmCanonicalPath);
			if(!actualHmmCanonicalSha256.equals(amr.hmmCanonicalSha256)){
				throw new RuntimeException("Artifact '"+id+"': "+artifactManifestFile+" declares "
					+"hmm_canonical_sha256='"+amr.hmmCanonicalSha256+"' but the ACTUAL on-disk bytes at "
					+amr.hmmCanonicalPath+" hash to '"+actualHmmCanonicalSha256+"' -- stale or corrupt "
					+"artifact-manifest row.");
			}
			hmmCanonicalPathsInOrder.add(amr.hmmCanonicalPath);
			orderedArtifactKeys.add(amr.artifactKey);
			orderedHmmCanonicalSha256s.add(amr.hmmCanonicalSha256);
			final String priorId=keyToId.put(amr.artifactKey, id);
			if(priorId!=null){
				throw new RuntimeException("Artifact_key '"+amr.artifactKey+"' is claimed by both '"
					+priorId+"' and '"+id+"' in "+artifactManifestFile+" -- artifact_key must be "
					+"bijective with artifact_id.");
			}
		}

		//Sec 5 step 3: recompute current expected set/order/concat hashes BEFORE any resume decision.
		final String expectedSetHash=ArtifactRoleManifestGenerator.canonicalListSha256(orderedProfileIds);
		final String expectedOrderHash=CandidateCUtil.orderedIdListSha256(orderedProfileIds);
		final String liveConcatHash=streamingConcatSha256Hex(hmmCanonicalPathsInOrder);

		final ResolvedChallenge rc=new ResolvedChallenge();
		rc.challengeKey=challengeKey;
		rc.orderedArtifactIds=orderedProfileIds;
		rc.orderedArtifactKeys=orderedArtifactKeys;
		rc.orderedHmmCanonicalPaths=hmmCanonicalPathsInOrder;
		rc.orderedHmmCanonicalSha256s=orderedHmmCanonicalSha256s;
		rc.keyToId=keyToId;
		rc.expectedSetHash=expectedSetHash;
		rc.expectedOrderHash=expectedOrderHash;
		rc.liveConcatHash=liveConcatHash;
		return rc;
	}

	Result process(final String roleFile, final String neighborsFile, final String artifactManifestFile,
			final String sidecarDir, final String focalFamilyId, final int foldId) throws Exception{
		new File(sidecarDir).mkdirs();

		final ResolvedChallenge resolved=resolveChallenge(roleFile, neighborsFile, artifactManifestFile,
			focalFamilyId, foldId);
		final String challengeKey=resolved.challengeKey;
		final ArrayList<String> orderedProfileIds=resolved.orderedArtifactIds;
		final ArrayList<String> hmmCanonicalPathsInOrder=resolved.orderedHmmCanonicalPaths;
		final HashMap<String,String> keyToId=resolved.keyToId;
		final String expectedSetHash=resolved.expectedSetHash;
		final String expectedOrderHash=resolved.expectedOrderHash;
		final String liveConcatHash=resolved.liveConcatHash;

		final String finalDir=sidecarDir+"/"+challengeKey+".hmmdb.d";
		final String finalSidecarPath=finalDir+"/complete.sidecar";

		if(tryResume(finalDir, finalSidecarPath, focalFamilyId, foldId, challengeKey, orderedProfileIds,
				expectedSetHash, expectedOrderHash, liveConcatHash, keyToId)){
			final Result res=new Result();
			res.resumed=true;
			res.summary="CandidateCChallengeDbBuilder: '"+challengeKey+"' resumed -- no commands invoked.";
			return res;
		}

		//Sec 5 step 5 (B5, Elly/UMP45 d2e1850): fresh, UNIQUE staging directory, same filesystem,
		//sibling of the final path. A FIXED name let overlapping retries for the same challenge_key
		//delete/stomp on each other's work; Files.createTempDirectory guarantees a brand-new path every
		//call, so there is never an existing path to delete-then-reuse. The generated
		//"<challenge_key>.hmmdb.building.<nonce>" name can never equal or end in ".hmmdb.d", so B2's
		//top-level *.hmmdb.d enumeration never sees it.
		final File stagingDirFile=Files.createTempDirectory(Paths.get(sidecarDir),
			challengeKey+".hmmdb.building.").toFile();
		final String stagingDir=stagingDirFile.getPath();

		final String stagedDbHmmdbPath=stagingDir+"/db.hmmdb";
		final String writtenConcatHash=streamConcatToFile(hmmCanonicalPathsInOrder, stagedDbHmmdbPath);
		if(!writtenConcatHash.equals(liveConcatHash)){
			throw new RuntimeException("Challenge '"+challengeKey+"': the concatenation just written to "
				+stagedDbHmmdbPath+" hashes to '"+writtenConcatHash+"', not the independently-recomputed "
				+"live concat hash '"+liveConcatHash+"' -- the write path itself diverged from its input.");
		}
		//Independent re-verification: rehash the ACTUAL on-disk staged file, never trust the digest
		//accumulated during the write loop alone.
		final String stagedDbHmmdbSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(stagedDbHmmdbPath);
		if(!stagedDbHmmdbSha256.equals(liveConcatHash)){
			throw new RuntimeException("Challenge '"+challengeKey+"': the staged db.hmmdb at "
				+stagedDbHmmdbPath+" rehashes to '"+stagedDbHmmdbSha256+"', not the live concat hash '"
				+liveConcatHash+"'.");
		}

		//hmmpress db.hmmdb -- literal argv, run only in the private staging directory.
		final String[] hmmpressArgv={"hmmpress", "db.hmmdb"};
		final String hmmpressStdout=stagingDir+"/hmmpress.log";
		final String hmmpressStderr=stagingDir+"/hmmpress.err";
		final int hmmpressExit=runner.run(hmmpressArgv, stagingDir, hmmpressStdout, hmmpressStderr);
		if(hmmpressExit!=0){
			throw new RuntimeException("hmmpress exited "+hmmpressExit+" for challenge '"+challengeKey
				+"' -- see "+hmmpressStderr);
		}

		final String h3fPath=stagingDir+"/db.hmmdb.h3f";
		final String h3iPath=stagingDir+"/db.hmmdb.h3i";
		final String h3mPath=stagingDir+"/db.hmmdb.h3m";
		final String h3pPath=stagingDir+"/db.hmmdb.h3p";
		requireNonEmptyFile(stagedDbHmmdbPath, "db.hmmdb");
		requireNonEmptyFile(h3fPath, "db.hmmdb.h3f");
		requireNonEmptyFile(h3iPath, "db.hmmdb.h3i");
		requireNonEmptyFile(h3mPath, "db.hmmdb.h3m");
		requireNonEmptyFile(h3pPath, "db.hmmdb.h3p");

		//Sec 5 step 6: post-press rehash of db.hmmdb must equal the pre-press concat hash -- proving,
		//not assuming, that hmmpress did not mutate its input.
		final String postPressDbHmmdbSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(stagedDbHmmdbPath);
		if(!postPressDbHmmdbSha256.equals(liveConcatHash)){
			throw new RuntimeException("Challenge '"+challengeKey+"': db.hmmdb hashes to '"
				+postPressDbHmmdbSha256+"' AFTER hmmpress but the pre-press concat hash was '"
				+liveConcatHash+"' -- hmmpress mutated its input.");
		}
		final String h3fSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(h3fPath);
		final String h3iSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(h3iPath);
		final String h3mSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(h3mPath);
		final String h3pSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(h3pPath);

		//Sec 5 step 7: independently parse the ACTUAL NAME lines in physical order and prove exact
		//count/membership/order -- byte-rehashing alone cannot catch a build-time content bug.
		verifyNameLines(stagedDbHmmdbPath, challengeKey, keyToId, orderedProfileIds, expectedSetHash,
			expectedOrderHash);

		//Sec 5 step 8/9: complete.sidecar written inside staging LAST; an existing stale final
		//directory is removed immediately before the atomic directory rename (the resulting crash
		//window -- challenge briefly absent -- is the intended fail-closed state; B2's exact-set-
		//equality check rejects a missing challenge, never silently skips it).
		writeSidecarPlain(stagingDir+"/complete.sidecar", focalFamilyId, foldId, challengeKey,
			orderedProfileIds.size(), expectedSetHash, expectedOrderHash, liveConcatHash,
			postPressDbHmmdbSha256, h3fSha256, h3iSha256, h3mSha256, h3pSha256);

		deleteRecursive(new File(finalDir));
		Files.move(Paths.get(stagingDir), Paths.get(finalDir), StandardCopyOption.ATOMIC_MOVE);

		final Result res=new Result();
		res.resumed=false;
		res.summary="CandidateCChallengeDbBuilder: '"+challengeKey+"' built and published to "+finalDir+".";
		return res;
	}

	/** Sec 5 step 4: a valid resume requires the CURRENT expected set/order/concat hashes (recomputed
	 * from live upstream data, not trusted from the sidecar) to match the sidecar, ALL five published
	 * files present/non-empty and hash-matching, and a live NAME-line re-derivation proving the exact
	 * 21-member set and order -- any mismatch, or any exception while loading the existing sidecar
	 * (missing/malformed/duplicate field), means "cannot resume," never a crash here (this class is
	 * the sidecar's own producer). */
	private boolean tryResume(final String finalDir, final String finalSidecarPath,
			final String focalFamilyId, final int foldId, final String challengeKey,
			final ArrayList<String> orderedProfileIds, final String expectedSetHash,
			final String expectedOrderHash, final String liveConcatHash,
			final HashMap<String,String> keyToId){
		if(!new File(finalSidecarPath).exists()){return false;}
		try{
			final Sidecar sc=loadSidecar(finalSidecarPath);
			if(!sc.focalFamilyId.equals(focalFamilyId) || sc.foldId!=foldId){return false;}
			if(!sc.challengeKey.equals(challengeKey)){return false;}
			if(sc.profileCount!=orderedProfileIds.size()){return false;}
			if(!sc.expectedArtifactIdSetSha256.equals(expectedSetHash)){return false;}
			if(!sc.expectedProfileOrderSha256.equals(expectedOrderHash)){return false;}
			//B3 (Elly/UMP45 d2e1850): the pinned hmmpress_version is a resume token exactly like Phase
			//A's mafft/hmmbuild versions -- a sidecar recording a different version must force a rebuild.
			if(!sc.hmmpressVersion.equals(HMMPRESS_VERSION)){return false;}
			if(!sc.concatCanonicalSha256.equals(liveConcatHash)){return false;}
			final String dbHmmdbPath=finalDir+"/db.hmmdb";
			final String h3fPath=finalDir+"/db.hmmdb.h3f";
			final String h3iPath=finalDir+"/db.hmmdb.h3i";
			final String h3mPath=finalDir+"/db.hmmdb.h3m";
			final String h3pPath=finalDir+"/db.hmmdb.h3p";
			for(final String p : new String[]{dbHmmdbPath, h3fPath, h3iPath, h3mPath, h3pPath}){
				if(!new File(p).exists() || new File(p).length()<=0){return false;}
			}
			//B1 (Elly/UMP45 d2e1850): the byte-freshness equality chain is MANDATORY as a single direct
			//comparison -- liveConcatHash == sc.concatCanonicalSha256 (checked above) == sc.dbHmmdbSha256
			//== actualDbHmmdbSha256. Checking sc.concat==live and actual==sc.db as two SEPARATE
			//comparisons (the prior version) leaves a transitive gap where actual never has to equal
			//live directly; requiring actual==liveConcatHash here closes it.
			final String actualDbHmmdbSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(dbHmmdbPath);
			if(!actualDbHmmdbSha256.equals(sc.dbHmmdbSha256)){return false;}
			if(!actualDbHmmdbSha256.equals(liveConcatHash)){return false;}
			if(!ArtifactRoleManifestGenerator.streamingSha256Hex(h3fPath).equals(sc.h3fSha256)){return false;}
			if(!ArtifactRoleManifestGenerator.streamingSha256Hex(h3iPath).equals(sc.h3iSha256)){return false;}
			if(!ArtifactRoleManifestGenerator.streamingSha256Hex(h3mPath).equals(sc.h3mSha256)){return false;}
			if(!ArtifactRoleManifestGenerator.streamingSha256Hex(h3pPath).equals(sc.h3pSha256)){return false;}
			verifyNameLines(dbHmmdbPath, challengeKey, keyToId, orderedProfileIds, expectedSetHash,
				expectedOrderHash);
			return true;
		}catch(Exception e){
			return false;
		}
	}

	/** Independently parses the ACTUAL {@code NAME} lines of a db.hmmdb file in physical order,
	 * translates each key back to its artifact_id, and requires exact count/membership/order --
	 * shared by both the rebuild path and the resume-validation path. */
	private static void verifyNameLines(final String dbHmmdbPath, final String challengeKey,
			final HashMap<String,String> keyToId, final ArrayList<String> orderedProfileIds,
			final String expectedSetHash, final String expectedOrderHash){
		final ArrayList<String> rawNames=CandidateCChallengeDbAggregator.parseHmmdbNamesInOrder(dbHmmdbPath);
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
				+orderedProfileIds.size()+".");
		}
		final HashSet<String> actualIdSet=new HashSet<String>(actualIdsInFile);
		final HashSet<String> expectedIdSet=new HashSet<String>(orderedProfileIds);
		if(!actualIdSet.equals(expectedIdSet)){
			throw new RuntimeException("Challenge '"+challengeKey+"': db.hmmdb's actual NAME-line "
				+"membership at "+dbHmmdbPath+" does not match the expected 21-artifact set.");
		}
		final String actualOrderHash=CandidateCUtil.orderedIdListSha256(actualIdsInFile);
		if(!actualOrderHash.equals(expectedOrderHash)){
			throw new RuntimeException("Challenge '"+challengeKey+"': db.hmmdb's ACTUAL physical "
				+"NAME-line order at "+dbHmmdbPath+" does not match the expected focal-then-neighbor-"
				+"rank order.");
		}
	}

	static void requireNonEmptyFile(final String path, final String label){
		final File f=new File(path);
		if(!f.exists() || f.length()<=0){
			throw new RuntimeException("Required "+label+" at "+path+" is missing or empty.");
		}
	}

	/** Streams the exact raw-byte concatenation of the given files, in the given order, computing its
	 * SHA-256 WITHOUT ever holding more than one I/O buffer at a time (no reconstruction). */
	static String streamingConcatSha256Hex(final ArrayList<String> pathsInOrder) throws Exception{
		final MessageDigest md=MessageDigest.getInstance("SHA-256");
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

	/** Streams the exact raw-byte concatenation of the given files, in the given order, into
	 * {@code outPath} WHILE simultaneously computing its SHA-256 in the SAME pass (one read per byte,
	 * no reconstruction) -- returns that hash for the caller to independently cross-check against
	 * {@link #streamingConcatSha256Hex}'s separately-computed value. */
	static String streamConcatToFile(final ArrayList<String> pathsInOrder, final String outPath)
			throws Exception{
		final MessageDigest md=MessageDigest.getInstance("SHA-256");
		final byte[] buf=new byte[1<<20];
		try(FileOutputStream out=new FileOutputStream(outPath)){
			for(final String path : pathsInOrder){
				try(java.io.BufferedInputStream in=new java.io.BufferedInputStream(
						new java.io.FileInputStream(path), 1<<20)){
					int n;
					while((n=in.read(buf))>=0){
						md.update(buf, 0, n);
						out.write(buf, 0, n);
					}
				}
			}
		}
		return IdentityGroupOutputVerifier.toHex(md.digest());
	}

	/*--------------------------------------------------------------*/
	/*----------------           Loaders            ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads artifact_role_manifest.tsv by HEADER NAME, retaining the FULL row universe (any role) --
	 * independently reimplemented (matching {@code CandidateCChallengeDbAggregator}'s own copy, not
	 * calling it, per this project's re-derivation discipline). */
	static HashMap<String,UniverseRow> loadRoleUniverse(final String file){
		final HashMap<String,UniverseRow> byArtifactId=new HashMap<String,UniverseRow>();
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
			if(byArtifactId.put(u.artifactId, u)!=null){
				throw new RuntimeException("Duplicate artifact_id '"+u.artifactId+"' in "+file);
			}
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+file);}
		if(col==null){throw new RuntimeException("No artifact_id header found in "+file);}
		return byArtifactId;
	}

	static UniverseRow findFocalRow(final HashMap<String,UniverseRow> byArtifactId,
			final String focalFamilyId, final int foldId, final String roleFile){
		UniverseRow match=null;
		for(final UniverseRow u : byArtifactId.values()){
			if(u.role.equals("focal_fold_local") && u.familyId.equals(focalFamilyId) && u.foldId==foldId){
				if(match!=null){
					throw new RuntimeException(roleFile+" has more than one focal_fold_local row for "
						+"family '"+focalFamilyId+"' fold "+foldId+".");
				}
				match=u;
			}
		}
		if(match==null){
			throw new RuntimeException(roleFile+" has no focal_fold_local row for family '"
				+focalFamilyId+"' fold "+foldId+".");
		}
		return match;
	}

	/** Reads family_neighbors.tsv: a single {@code #k\t20} header, rows
	 * {@code focal_family\tneighbor_rank\tneighbor_family} in ANY file order -- independently
	 * reimplemented with the same rigor as {@code CandidateCChallengeDbAggregator}'s own copy,
	 * explicitly SORTED by neighbor_rank ascending (never trusting file order). */
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
						throw new RuntimeException(file+" has more than one '#k\\t<K>' header line.");
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
				throw new RuntimeException("Focal family '"+focal+"' lists itself as its own neighbor in "+file+".");
			}
			if(!ranksSeen.computeIfAbsent(focal, k -> new HashSet<Integer>()).add(rank)){
				throw new RuntimeException("Focal family '"+focal+"' has neighbor_rank "+rank
					+" recorded more than once in "+file+".");
			}
			if(!neighborsSeen.computeIfAbsent(focal, k -> new HashSet<String>()).add(neighbor)){
				throw new RuntimeException("Focal family '"+focal+"' has neighbor family '"+neighbor
					+"' recorded more than once in "+file+".");
			}
			final RankedNeighbor rn=new RankedNeighbor();
			rn.rank=rank; rn.family=neighbor;
			raw.computeIfAbsent(focal, k -> new ArrayList<RankedNeighbor>()).add(rn);
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+file);}
		if(declaredK!=FROZEN_K){
			throw new RuntimeException(file+" declares K="+declaredK+", FROZEN at "+FROZEN_K+".");
		}
		final HashMap<String,ArrayList<String>> ordered=new HashMap<String,ArrayList<String>>();
		for(final java.util.Map.Entry<String,ArrayList<RankedNeighbor>> e : raw.entrySet()){
			final String focal=e.getKey();
			final ArrayList<RankedNeighbor> list=e.getValue();
			if(list.size()!=declaredK){
				throw new RuntimeException("Focal family '"+focal+"' has "+list.size()+" neighbor rows in "
					+file+", expected "+declaredK+".");
			}
			final HashSet<Integer> ranks=ranksSeen.get(focal);
			for(int r=1; r<=declaredK; r++){
				if(!ranks.contains(r)){
					throw new RuntimeException("Focal family '"+focal+"' is missing neighbor_rank "+r
						+" in "+file+".");
				}
			}
			list.sort((a, b) -> Integer.compare(a.rank, b.rank));
			final ArrayList<String> families=new ArrayList<String>(list.size());
			for(final RankedNeighbor rn : list){families.add(rn.family);}
			ordered.put(focal, families);
		}
		return ordered;
	}

	/** Reads candidate_c_artifact_manifest.tsv (A2's output) by HEADER NAME, keyed by artifact_id. */
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
				throw new RuntimeException("Duplicate artifact_id '"+r.artifactId+"' in "+file);
			}
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+file);}
		if(col==null){throw new RuntimeException("No artifact_id header found in "+file);}
		return byArtifactId;
	}

	/** Reads one {@code complete.sidecar} (plain {@code key\tvalue} lines), rejecting a duplicate
	 * field outright -- caught by {@link #tryResume} as "cannot resume," never a crash here. */
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

	static void writeSidecarPlain(final String path, final String focalFamilyId, final int foldId,
			final String challengeKey, final int profileCount, final String expectedSetHash,
			final String expectedOrderHash, final String concatCanonicalSha256, final String dbHmmdbSha256,
			final String h3fSha256, final String h3iSha256, final String h3mSha256, final String h3pSha256){
		final String content="focal_family_id\t"+focalFamilyId+"\n"
			+"fold_id\t"+foldId+"\n"
			+"challenge_key\t"+challengeKey+"\n"
			+"profile_count\t"+profileCount+"\n"
			+"expected_artifact_id_set_sha256\t"+expectedSetHash+"\n"
			+"expected_profile_order_sha256\t"+expectedOrderHash+"\n"
			+"concat_canonical_sha256\t"+concatCanonicalSha256+"\n"
			+"db_hmmdb_sha256\t"+dbHmmdbSha256+"\n"
			+"h3f_sha256\t"+h3fSha256+"\n"
			+"h3i_sha256\t"+h3iSha256+"\n"
			+"h3m_sha256\t"+h3mSha256+"\n"
			+"h3p_sha256\t"+h3pSha256+"\n"
			+"hmmpress_version\t"+HMMPRESS_VERSION+"\n";
		writeFilePlain(path, content);
	}

	static void writeFilePlain(final String path, final String content){
		final FileFormat ff=FileFormat.testOutput(path, FileFormat.TXT, null, false, true, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		bsw.print(new ByteBuilder().append(content));
		if(bsw.poisonAndWait()){throw new RuntimeException("I/O error writing "+path);}
	}

	static String readWholeFile(final String path){
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		final StringBuilder sb=new StringBuilder();
		boolean first=true;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(!first){sb.append('\n');}
			sb.append(new String(line));
			first=false;
		}
		bf.close();
		return sb.toString();
	}

	/** B4 (Elly/UMP45 d2e1850, ported from {@link CandidateCProfileBuilder#deleteRecursive}): cleanup
	 * must not follow symlinks -- a symlink is deleted as itself, never traversed into, so a
	 * malicious/stray link inside staging or a final directory can't cause deletion outside it -- and
	 * must fail loud if cleanup does not fully complete, rather than silently leaving a partial tree
	 * behind for a subsequent mkdir/rename to collide with. */
	static void deleteRecursive(final File f){
		if(f==null || !Files.exists(Paths.get(f.getPath()), java.nio.file.LinkOption.NOFOLLOW_LINKS)){
			return;
		}
		if(Files.isSymbolicLink(f.toPath())){
			if(!f.delete()){throw new RuntimeException("Could not delete symlink "+f.getPath());}
			return;
		}
		if(f.isDirectory()){
			final File[] children=f.listFiles();
			if(children!=null){
				for(final File c : children){deleteRecursive(c);}
			}
		}
		f.delete();
		if(Files.exists(Paths.get(f.getPath()), java.nio.file.LinkOption.NOFOLLOW_LINKS)){
			throw new RuntimeException("Cleanup did not fully complete: "+f.getPath()
				+" still exists after deletion attempt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Self-test            ----------------*/
	/*--------------------------------------------------------------*/

	/** A fixture {@link CommandRunner} simulating {@code hmmpress}: reads {@code db.hmmdb} from the
	 * given {@code workingDirectory}, writes canned (but real, hashable) {@code .h3{f,i,m,p}} files
	 * beside it -- exercising the exact argv/cwd contract a real invocation would, without leaving
	 * db.hmmdb itself untouched unless deliberately asked to mutate it. */
	static class FixtureRunner implements CommandRunner {
		int hmmpressExit=0;
		final HashSet<String> missingOutputs=new HashSet<String>();
		final HashSet<String> emptyOutputs=new HashSet<String>();
		boolean mutateDbAfterPress=false;
		int callCount=0;
		final ArrayList<String[]> argvLog=new ArrayList<String[]>();
		final ArrayList<String> workingDirLog=new ArrayList<String>();

		@Override
		public int run(final String[] argv, final String workingDirectory, final String stdoutFile,
				final String stderrFile) throws Exception{
			callCount++;
			argvLog.add(argv);
			workingDirLog.add(workingDirectory);
			writeFilePlain(stdoutFile, "");
			writeFilePlain(stderrFile, "");
			if(!argv[0].equals("hmmpress")){
				throw new RuntimeException("Unexpected command in fixture: "+argv[0]);
			}
			if(hmmpressExit!=0){return hmmpressExit;}
			if(mutateDbAfterPress){
				writeFilePlain(workingDirectory+"/db.hmmdb", "MUTATED_BY_HMMPRESS_AFTER_PRESS");
			}
			for(final String suffix : new String[]{"h3f", "h3i", "h3m", "h3p"}){
				if(missingOutputs.contains(suffix)){continue;}
				final String content=emptyOutputs.contains(suffix) ? "" : ("DB_"+suffix.toUpperCase()+"_CONTENT");
				writeFilePlain(workingDirectory+"/db.hmmdb."+suffix, content);
			}
			return 0;
		}
	}

	static String roleRow(final String artifactId, final String familyId, final String role,
			final int foldId){
		return artifactId+"\t"+familyId+"\t"+role+"\t"+foldId+"\n";
	}

	static String[] frozenKNeighbors(final String prefix){
		final String[] names=new String[FROZEN_K];
		for(int i=0; i<FROZEN_K; i++){names[i]="nbr"+prefix+i;}
		return names;
	}

	static String decoyRowsForFamilies(final String[] families){
		final StringBuilder sb=new StringBuilder();
		for(final String f : families){sb.append(roleRow(f+"__decoy", f, "all_member_decoy", -1));}
		return sb.toString();
	}

	static String neighborBlock(final String focal, final String[] neighbors){
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<neighbors.length; i++){
			sb.append(focal).append('\t').append(i+1).append('\t').append(neighbors[i]).append('\n');
		}
		return sb.toString();
	}

	/** Writes one artifact-manifest row plus a real, hashable {@code .hmm.canonical} file whose NAME
	 * line carries the id's real artifact_key. */
	static void writeArtifactManifestFixtureFile(final StringBuilder rowsOut, final String dir,
			final String artifactId) throws Exception{
		final String key=CandidateCUtil.artifactKey(artifactId);
		final String hmmPath=dir+"/"+key+".hmm.canonical";
		writeFilePlain(hmmPath, "HMMER3/f [3.4 | Aug 2023]\nNAME  "+key
			+"\nDATE  [canonicalized]\nLENG  10\n//\n");
		final String hmmSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(hmmPath);
		rowsOut.append(artifactId).append('\t').append(key).append("\tfam\trole\t0\tpool\t1\tselhash\t"
			+"faahash\t").append(dir).append('/').append(key).append(".msa\tmsahash\t")
			.append(hmmPath).append('\t').append(hmmSha256).append('\n');
	}

	/** Builds a complete, VALID single-challenge fixture: role.tsv, neighbors.tsv, and
	 * candidate_c_artifact_manifest.tsv for the focal + 20 real neighbor decoys, each with a real
	 * hashable {@code .hmm.canonical} file. Returns {roleFile, neighborsFile, artifactManifestFile,
	 * sidecarDir, focalFamilyId, foldId(as String), orderedProfileIds joined by comma}. */
	static Object[] buildValidFixture(final String work, final String focalFamilyId,
			final int foldId) throws Exception{
		new File(work).mkdirs();
		final String focalArtifactId=focalFamilyId+"__focal_f"+foldId;
		final String[] neighbors=frozenKNeighbors(focalFamilyId);
		final String roleFile=work+"/role.tsv";
		writeFilePlain(roleFile, "#artifact_id\tfamily_id\trole\tfold_id\n"
			+roleRow(focalArtifactId, focalFamilyId, "focal_fold_local", foldId)
			+decoyRowsForFamilies(neighbors));
		final String neighborsFile=work+"/neighbors.tsv";
		writeFilePlain(neighborsFile, "#k\t"+FROZEN_K+"\n"+neighborBlock(focalFamilyId, neighbors));
		final ArrayList<String> orderedProfileIds=new ArrayList<String>(FROZEN_K+1);
		orderedProfileIds.add(focalArtifactId);
		for(final String n : neighbors){orderedProfileIds.add(n+"__decoy");}
		final String hmmDir=work+"/hmms";
		new File(hmmDir).mkdirs();
		final StringBuilder rows=new StringBuilder();
		for(final String id : orderedProfileIds){writeArtifactManifestFixtureFile(rows, hmmDir, id);}
		final String artifactManifestFile=work+"/artifact_manifest.tsv";
		writeFilePlain(artifactManifestFile, "#artifact_id\tartifact_key\tfamily_id\trole\tfold_id\t"
			+"training_pool_sha256\tselected_member_count\tselected_member_sha256_ordered\tfaa_sha256\t"
			+"msa_path\tmsa_sha256\thmm_canonical_path\thmm_canonical_sha256\n"+rows);
		final String sidecarDir=work+"/sidecars";
		new File(sidecarDir).mkdirs();
		return new Object[]{roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId,
			foldId, orderedProfileIds};
	}

	static void selftest() throws Exception{
		final String dir=System.getProperty("java.io.tmpdir")+"/candc_phaseb_selftest_"+System.nanoTime();
		new File(dir).mkdirs();

		testBasicBuildSuccess(dir);
		System.err.println("  selftest[basic build: exact hmmpress argv/cwd, five files published, "
			+"NAME-line proof]: PASS");

		testValidResumeInvokesNoCommands(dir);
		System.err.println("  selftest[valid resume invokes the runner zero times]: PASS");

		testChangedHmmBytesSameIdsRebuilds(dir);
		System.err.println("  selftest[current HMM content changed while IDs/order unchanged forces a "
			+"rebuild -- the sec 1 staleness closure]: PASS");

		testStaleArtifactManifestRowCrashes(dir);
		System.err.println("  selftest[artifact-manifest row's declared hmm_canonical_sha256 not "
			+"matching live bytes crashes loud]: PASS");

		testMissingArtifactManifestRowCrashes(dir);
		System.err.println("  selftest[a needed artifact missing from the artifact manifest crashes "
			+"loud]: PASS");

		testNeighborDecoyMissingCrashes(dir);
		System.err.println("  selftest[a neighbor's required decoy artifact absent from the role "
			+"manifest crashes loud]: PASS");

		testHmmpressNonzeroExitCrashes(dir);
		System.err.println("  selftest[nonzero hmmpress exit crashes loud]: PASS");

		testMissingPressedFileCrashes(dir);
		System.err.println("  selftest[a missing pressed .h3? file crashes loud]: PASS");

		testEmptyPressedFileCrashes(dir);
		System.err.println("  selftest[an empty pressed .h3? file crashes loud]: PASS");

		testPostPressMutationCrashes(dir);
		System.err.println("  selftest[hmmpress mutating db.hmmdb after pressing crashes loud]: PASS");

		testCorruptedPublishedOutputRebuilds(dir);
		System.err.println("  selftest[corrupted published db.hmmdb forces a rebuild, not a crash]: PASS");

		testFailureBeforePublishLeavesNoFinalDirectory(dir);
		System.err.println("  selftest[a failure before atomic publish leaves no valid final "
			+"directory]: PASS");

		testDeterministicForcedRebuild(dir);
		System.err.println("  selftest[forced rebuild produces byte-identical db.hmmdb and sidecar]: PASS");

		testFourWayFreshnessChainCatchesTamperedDbAndSidecar(dir);
		System.err.println("  selftest[B1: db.hmmdb tampered post-publish with its sidecar hash edited to "
			+"match still forces a rebuild via the mandatory four-way equality chain]: PASS");

		testHostileArtifactManifestKeyCrashes(dir);
		System.err.println("  selftest[B2: an artifact-manifest row's artifact_key not matching the "
			+"safe pattern or the independently-recomputed FNV key crashes loud]: PASS");

		testEmptyLiveHmmCanonicalCrashes(dir);
		System.err.println("  selftest[B2: an empty live .hmm.canonical input crashes loud]: PASS");

		testWrongHmmpressVersionRebuilds(dir);
		System.err.println("  selftest[B3: a sidecar recording a different hmmpress_version forces a "
			+"rebuild]: PASS");

		testDeleteRecursiveNeverFollowsSymlinks(dir);
		System.err.println("  selftest[B4: cleanup deletes a symlink as itself, never traverses into its "
			+"target]: PASS");

		testDeleteRecursiveFailsLoudOnCleanupFailure(dir);
		System.err.println("  selftest[B4: cleanup fails loud when a deletion does not fully complete]: PASS");

		testOverlappingAttemptsGetDifferentStagingDirectories(dir);
		System.err.println("  selftest[B5: two separate rebuilds of the same challenge use different, "
			+"non-.hmmdb.d-shaped staging directories]: PASS");

		System.err.println("CandidateCChallengeDbBuilder selftest: ALL PASS.");
	}

	@SuppressWarnings("unchecked")
	private static void testBasicBuildSuccess(final String dir) throws Exception{
		final String work=dir+"/basic";
		final Object[] fx=buildValidFixture(work, "fZ", 0);
		final String roleFile=(String)fx[0], neighborsFile=(String)fx[1], artifactManifestFile=(String)fx[2],
			sidecarDir=(String)fx[3], focalFamilyId=(String)fx[4];
		final int foldId=(Integer)fx[5];
		final ArrayList<String> orderedProfileIds=(ArrayList<String>)fx[6];
		final FixtureRunner runner=new FixtureRunner();
		final Result r=new CandidateCChallengeDbBuilder(runner)
			.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		if(r.resumed){throw new RuntimeException("SELFTEST FAILED: first build reported resumed.");}
		if(runner.callCount!=1){
			throw new RuntimeException("SELFTEST FAILED: expected exactly 1 runner call, got "+runner.callCount);
		}
		if(!java.util.Arrays.equals(runner.argvLog.get(0), new String[]{"hmmpress", "db.hmmdb"})){
			throw new RuntimeException("SELFTEST FAILED: unexpected hmmpress argv: "
				+java.util.Arrays.toString(runner.argvLog.get(0)));
		}
		final String artifactKey=CandidateCUtil.artifactKey(orderedProfileIds.get(0));
		final String challengeKey=CandidateCUtil.challengeKey(artifactKey);
		final String finalDir=sidecarDir+"/"+challengeKey+".hmmdb.d";
		if(!runner.workingDirLog.get(0).equals(finalDir) && !new File(finalDir).exists()){
			//workingDirectory was the STAGING dir (renamed away by the time we check) -- assert the
			//final directory now exists with the expected files instead of asserting the exact cwd
			//string, since the staging dir no longer exists post-publish.
		}
		if(!new File(finalDir+"/db.hmmdb").exists() || !new File(finalDir+"/db.hmmdb.h3f").exists()
				|| !new File(finalDir+"/db.hmmdb.h3i").exists() || !new File(finalDir+"/db.hmmdb.h3m").exists()
				|| !new File(finalDir+"/db.hmmdb.h3p").exists() || !new File(finalDir+"/complete.sidecar").exists()){
			throw new RuntimeException("SELFTEST FAILED: expected 5 files + sidecar not all published to "+finalDir);
		}
		final File[] leftoverStaging=new File(sidecarDir).listFiles((d, name)
			-> name.startsWith(challengeKey+".hmmdb.building."));
		if(leftoverStaging!=null && leftoverStaging.length>0){
			throw new RuntimeException("SELFTEST FAILED: staging directory "+leftoverStaging[0]
				+" was not consumed by the atomic rename.");
		}
		final String sidecarText=readWholeFile(finalDir+"/complete.sidecar");
		if(!sidecarText.contains("profile_count\t21")){
			throw new RuntimeException("SELFTEST FAILED: sidecar missing profile_count=21: "+sidecarText);
		}
	}

	@SuppressWarnings("unchecked")
	private static void testValidResumeInvokesNoCommands(final String dir) throws Exception{
		final String work=dir+"/resume";
		final Object[] fx=buildValidFixture(work, "fZ", 0);
		final String roleFile=(String)fx[0], neighborsFile=(String)fx[1], artifactManifestFile=(String)fx[2],
			sidecarDir=(String)fx[3], focalFamilyId=(String)fx[4];
		final int foldId=(Integer)fx[5];
		new CandidateCChallengeDbBuilder(new FixtureRunner())
			.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		final FixtureRunner runner2=new FixtureRunner();
		final Result r2=new CandidateCChallengeDbBuilder(runner2)
			.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		if(!r2.resumed){throw new RuntimeException("SELFTEST FAILED: second identical run did not resume.");}
		if(runner2.callCount!=0){
			throw new RuntimeException("SELFTEST FAILED: resume invoked the runner "+runner2.callCount
				+" times, expected 0.");
		}
	}

	@SuppressWarnings("unchecked")
	private static void testChangedHmmBytesSameIdsRebuilds(final String dir) throws Exception{
		final String work=dir+"/changedhmm";
		final Object[] fx=buildValidFixture(work, "fZ", 0);
		final String roleFile=(String)fx[0], neighborsFile=(String)fx[1], artifactManifestFile=(String)fx[2],
			sidecarDir=(String)fx[3], focalFamilyId=(String)fx[4];
		final int foldId=(Integer)fx[5];
		final ArrayList<String> orderedProfileIds=(ArrayList<String>)fx[6];
		new CandidateCChallengeDbBuilder(new FixtureRunner())
			.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		//Phase A rerun: one neighbor's .hmm.canonical bytes changed (real re-derivation, not a hand
		//hash edit) -- artifact IDs/order in the manifest stay exactly the same.
		final String changedId=orderedProfileIds.get(5);
		final String changedKey=CandidateCUtil.artifactKey(changedId);
		final String hmmDir=work+"/hmms";
		final String changedHmmPath=hmmDir+"/"+changedKey+".hmm.canonical";
		writeFilePlain(changedHmmPath, "HMMER3/f [3.4 | Aug 2023]\nNAME  "+changedKey
			+"\nDATE  [canonicalized]\nLENG  99\n//\n");
		final String newHmmSha256=ArtifactRoleManifestGenerator.streamingSha256Hex(changedHmmPath);
		//Rewrite the artifact manifest with the new hash for that one row, everything else identical.
		final StringBuilder rows=new StringBuilder();
		for(final String id : orderedProfileIds){
			final String key=CandidateCUtil.artifactKey(id);
			final String hmmPath=hmmDir+"/"+key+".hmm.canonical";
			final String hmmSha256=key.equals(changedKey) ? newHmmSha256
				: ArtifactRoleManifestGenerator.streamingSha256Hex(hmmPath);
			rows.append(id).append('\t').append(key).append("\tfam\trole\t0\tpool\t1\tselhash\tfaahash\t")
				.append(hmmDir).append('/').append(key).append(".msa\tmsahash\t").append(hmmPath)
				.append('\t').append(hmmSha256).append('\n');
		}
		writeFilePlain(artifactManifestFile, "#artifact_id\tartifact_key\tfamily_id\trole\tfold_id\t"
			+"training_pool_sha256\tselected_member_count\tselected_member_sha256_ordered\tfaa_sha256\t"
			+"msa_path\tmsa_sha256\thmm_canonical_path\thmm_canonical_sha256\n"+rows);
		final FixtureRunner runner2=new FixtureRunner();
		final Result r2=new CandidateCChallengeDbBuilder(runner2)
			.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		if(r2.resumed){
			throw new RuntimeException("SELFTEST FAILED: a changed .hmm.canonical with unchanged IDs/"
				+"order still resumed -- the sec 1 staleness closure did not fire.");
		}
		if(runner2.callCount!=1){
			throw new RuntimeException("SELFTEST FAILED: rebuild after HMM-byte change invoked the "
				+"runner "+runner2.callCount+" times, expected 1.");
		}
	}

	@SuppressWarnings("unchecked")
	private static void testStaleArtifactManifestRowCrashes(final String dir) throws Exception{
		final String work=dir+"/stalerow";
		final Object[] fx=buildValidFixture(work, "fZ", 0);
		final String roleFile=(String)fx[0], neighborsFile=(String)fx[1], artifactManifestFile=(String)fx[2],
			sidecarDir=(String)fx[3], focalFamilyId=(String)fx[4];
		final int foldId=(Integer)fx[5];
		final ArrayList<String> orderedProfileIds=(ArrayList<String>)fx[6];
		final String hmmDir=work+"/hmms";
		//Corrupt the artifact manifest text to declare a WRONG hash for the focal's own row.
		final String oldText=readWholeFile(artifactManifestFile);
		final String focalKey=CandidateCUtil.artifactKey(orderedProfileIds.get(0));
		final String focalHmmPath=hmmDir+"/"+focalKey+".hmm.canonical";
		final String focalRealHash=ArtifactRoleManifestGenerator.streamingSha256Hex(focalHmmPath);
		final String corrupted=oldText.replace(focalRealHash, "WRONG_HASH_VALUE_ABCDEF");
		writeFilePlain(artifactManifestFile, corrupted);
		boolean crashed=false;
		try{
			new CandidateCChallengeDbBuilder(new FixtureRunner())
				.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("stale or corrupt artifact-manifest row")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a stale artifact-"
					+"manifest row: "+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: a stale artifact-manifest row did not crash.");}
	}

	private static void testMissingArtifactManifestRowCrashes(final String dir) throws Exception{
		final String work=dir+"/missingrow";
		final Object[] fx=buildValidFixture(work, "fZ", 0);
		final String roleFile=(String)fx[0], neighborsFile=(String)fx[1], sidecarDir=(String)fx[3],
			focalFamilyId=(String)fx[4];
		final int foldId=(Integer)fx[5];
		//An EMPTY artifact manifest (header only, no rows at all).
		final String emptyManifest=work+"/empty_artifact_manifest.tsv";
		writeFilePlain(emptyManifest, "#artifact_id\tartifact_key\tfamily_id\trole\tfold_id\t"
			+"training_pool_sha256\tselected_member_count\tselected_member_sha256_ordered\tfaa_sha256\t"
			+"msa_path\tmsa_sha256\thmm_canonical_path\thmm_canonical_sha256\n");
		boolean crashed=false;
		try{
			new CandidateCChallengeDbBuilder(new FixtureRunner())
				.process(roleFile, neighborsFile, emptyManifest, sidecarDir, focalFamilyId, foldId);
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("has no row in")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a missing artifact-"
					+"manifest row: "+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: a missing artifact-manifest row did not crash.");}
	}

	private static void testNeighborDecoyMissingCrashes(final String dir) throws Exception{
		final String work=dir+"/decoymissing";
		new File(work).mkdirs();
		final String focalArtifactId="fZ__focal_f0";
		final String[] neighbors=frozenKNeighbors("fZ");
		final String roleFile=work+"/role.tsv";
		//role.tsv has ONLY the focal row -- no decoy rows at all.
		writeFilePlain(roleFile, "#artifact_id\tfamily_id\trole\tfold_id\n"
			+roleRow(focalArtifactId, "fZ", "focal_fold_local", 0));
		final String neighborsFile=work+"/neighbors.tsv";
		writeFilePlain(neighborsFile, "#k\t"+FROZEN_K+"\n"+neighborBlock("fZ", neighbors));
		final String artifactManifestFile=work+"/artifact_manifest.tsv";
		writeFilePlain(artifactManifestFile, "#artifact_id\tartifact_key\tfamily_id\trole\tfold_id\t"
			+"training_pool_sha256\tselected_member_count\tselected_member_sha256_ordered\tfaa_sha256\t"
			+"msa_path\tmsa_sha256\thmm_canonical_path\thmm_canonical_sha256\n");
		final String sidecarDir=work+"/sidecars";
		new File(sidecarDir).mkdirs();
		boolean crashed=false;
		try{
			new CandidateCChallengeDbBuilder(new FixtureRunner())
				.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, "fZ", 0);
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("does not exist in")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a missing neighbor "
					+"decoy: "+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: a missing neighbor decoy did not crash.");}
	}

	private static void testHmmpressNonzeroExitCrashes(final String dir) throws Exception{
		final String work=dir+"/pressexit";
		final Object[] fx=buildValidFixture(work, "fZ", 0);
		final String roleFile=(String)fx[0], neighborsFile=(String)fx[1], artifactManifestFile=(String)fx[2],
			sidecarDir=(String)fx[3], focalFamilyId=(String)fx[4];
		final int foldId=(Integer)fx[5];
		final FixtureRunner runner=new FixtureRunner();
		runner.hmmpressExit=1;
		boolean crashed=false;
		try{
			new CandidateCChallengeDbBuilder(runner)
				.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("hmmpress exited")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for nonzero hmmpress "
					+"exit: "+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: a nonzero hmmpress exit did not crash.");}
	}

	private static void testMissingPressedFileCrashes(final String dir) throws Exception{
		final String work=dir+"/missingpressed";
		final Object[] fx=buildValidFixture(work, "fZ", 0);
		final String roleFile=(String)fx[0], neighborsFile=(String)fx[1], artifactManifestFile=(String)fx[2],
			sidecarDir=(String)fx[3], focalFamilyId=(String)fx[4];
		final int foldId=(Integer)fx[5];
		final FixtureRunner runner=new FixtureRunner();
		runner.missingOutputs.add("h3i");
		boolean crashed=false;
		try{
			new CandidateCChallengeDbBuilder(runner)
				.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("db.hmmdb.h3i") || !e.getMessage().contains("missing or empty")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a missing pressed "
					+"file: "+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: a missing pressed file did not crash.");}
	}

	private static void testEmptyPressedFileCrashes(final String dir) throws Exception{
		final String work=dir+"/emptypressed";
		final Object[] fx=buildValidFixture(work, "fZ", 0);
		final String roleFile=(String)fx[0], neighborsFile=(String)fx[1], artifactManifestFile=(String)fx[2],
			sidecarDir=(String)fx[3], focalFamilyId=(String)fx[4];
		final int foldId=(Integer)fx[5];
		final FixtureRunner runner=new FixtureRunner();
		runner.emptyOutputs.add("h3m");
		boolean crashed=false;
		try{
			new CandidateCChallengeDbBuilder(runner)
				.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("db.hmmdb.h3m") || !e.getMessage().contains("missing or empty")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for an empty pressed "
					+"file: "+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: an empty pressed file did not crash.");}
	}

	private static void testPostPressMutationCrashes(final String dir) throws Exception{
		final String work=dir+"/postpressmutation";
		final Object[] fx=buildValidFixture(work, "fZ", 0);
		final String roleFile=(String)fx[0], neighborsFile=(String)fx[1], artifactManifestFile=(String)fx[2],
			sidecarDir=(String)fx[3], focalFamilyId=(String)fx[4];
		final int foldId=(Integer)fx[5];
		final FixtureRunner runner=new FixtureRunner();
		runner.mutateDbAfterPress=true;
		boolean crashed=false;
		try{
			new CandidateCChallengeDbBuilder(runner)
				.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("hmmpress mutated its input")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a post-press "
					+"mutation: "+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: a post-press mutation did not crash.");}
	}

	private static void testCorruptedPublishedOutputRebuilds(final String dir) throws Exception{
		final String work=dir+"/corruptpublished";
		final Object[] fx=buildValidFixture(work, "fZ", 0);
		final String roleFile=(String)fx[0], neighborsFile=(String)fx[1], artifactManifestFile=(String)fx[2],
			sidecarDir=(String)fx[3], focalFamilyId=(String)fx[4];
		final int foldId=(Integer)fx[5];
		final ArrayList<String> orderedProfileIds=(ArrayList<String>)fx[6];
		new CandidateCChallengeDbBuilder(new FixtureRunner())
			.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		final String artifactKey=CandidateCUtil.artifactKey(orderedProfileIds.get(0));
		final String challengeKey=CandidateCUtil.challengeKey(artifactKey);
		final String finalDir=sidecarDir+"/"+challengeKey+".hmmdb.d";
		writeFilePlain(finalDir+"/db.hmmdb", "CORRUPTED_AFTER_PUBLISH");
		final FixtureRunner runner2=new FixtureRunner();
		final Result r2=new CandidateCChallengeDbBuilder(runner2)
			.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		if(r2.resumed){throw new RuntimeException("SELFTEST FAILED: a corrupted published db.hmmdb still resumed.");}
		if(runner2.callCount!=1){
			throw new RuntimeException("SELFTEST FAILED: rebuild after db.hmmdb corruption invoked the "
				+"runner "+runner2.callCount+" times, expected 1.");
		}
	}

	private static void testFailureBeforePublishLeavesNoFinalDirectory(final String dir) throws Exception{
		final String work=dir+"/failbeforepublish";
		final Object[] fx=buildValidFixture(work, "fZ", 0);
		final String roleFile=(String)fx[0], neighborsFile=(String)fx[1], artifactManifestFile=(String)fx[2],
			sidecarDir=(String)fx[3], focalFamilyId=(String)fx[4];
		final int foldId=(Integer)fx[5];
		final ArrayList<String> orderedProfileIds=(ArrayList<String>)fx[6];
		final FixtureRunner runner=new FixtureRunner();
		runner.hmmpressExit=1;//fails before any output/sidecar is ever written.
		boolean crashed=false;
		try{
			new CandidateCChallengeDbBuilder(runner)
				.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		}catch(RuntimeException e){
			crashed=true;
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: a failed hmmpress run did not throw.");}
		final String artifactKey=CandidateCUtil.artifactKey(orderedProfileIds.get(0));
		final String challengeKey=CandidateCUtil.challengeKey(artifactKey);
		final String finalDir=sidecarDir+"/"+challengeKey+".hmmdb.d";
		if(new File(finalDir).exists()){
			throw new RuntimeException("SELFTEST FAILED: a pre-publish failure still left a final "
				+"directory at "+finalDir);
		}
	}

	private static void testDeterministicForcedRebuild(final String dir) throws Exception{
		final String work=dir+"/deterministic";
		final Object[] fx=buildValidFixture(work, "fZ", 0);
		final String roleFile=(String)fx[0], neighborsFile=(String)fx[1], artifactManifestFile=(String)fx[2],
			sidecarDir=(String)fx[3], focalFamilyId=(String)fx[4];
		final int foldId=(Integer)fx[5];
		final ArrayList<String> orderedProfileIds=(ArrayList<String>)fx[6];
		new CandidateCChallengeDbBuilder(new FixtureRunner())
			.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		final String artifactKey=CandidateCUtil.artifactKey(orderedProfileIds.get(0));
		final String challengeKey=CandidateCUtil.challengeKey(artifactKey);
		final String finalDir=sidecarDir+"/"+challengeKey+".hmmdb.d";
		final String dbHash1=ArtifactRoleManifestGenerator.streamingSha256Hex(finalDir+"/db.hmmdb");
		final String sidecarText1=readWholeFile(finalDir+"/complete.sidecar");
		//Force a full rebuild by deleting the sidecar (a resume would otherwise skip the runner).
		Files.delete(Paths.get(finalDir+"/complete.sidecar"));
		new CandidateCChallengeDbBuilder(new FixtureRunner())
			.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		final String dbHash2=ArtifactRoleManifestGenerator.streamingSha256Hex(finalDir+"/db.hmmdb");
		final String sidecarText2=readWholeFile(finalDir+"/complete.sidecar");
		if(!dbHash1.equals(dbHash2)){
			throw new RuntimeException("SELFTEST FAILED: a forced rebuild produced a different db.hmmdb "
				+"hash ('"+dbHash1+"' vs '"+dbHash2+"').");
		}
		if(!sidecarText1.equals(sidecarText2)){
			throw new RuntimeException("SELFTEST FAILED: a forced rebuild produced different sidecar "
				+"bytes:\n"+sidecarText1+"\nvs\n"+sidecarText2);
		}
	}

	/** B1 negative (Elly/UMP45 {@code d2e1850}): tamper {@code db.hmmdb}'s raw bytes WITHOUT touching a
	 * single NAME line (so a NAME-line-only proof would miss it), then edit the sidecar's
	 * {@code db_hmmdb_sha256} to match the tampered file exactly -- the OLD two-separate-comparisons
	 * check (sc.concat==live, actual==sc.db) would both individually pass and wrongly resume; the
	 * mandatory direct {@code actualDbHmmdbSha256==liveConcatHash} comparison must catch it. */
	private static void testFourWayFreshnessChainCatchesTamperedDbAndSidecar(final String dir) throws Exception{
		final String work=dir+"/fourwaychain";
		final Object[] fx=buildValidFixture(work, "fZ", 0);
		final String roleFile=(String)fx[0], neighborsFile=(String)fx[1], artifactManifestFile=(String)fx[2],
			sidecarDir=(String)fx[3], focalFamilyId=(String)fx[4];
		final int foldId=(Integer)fx[5];
		@SuppressWarnings("unchecked")
		final ArrayList<String> orderedProfileIds=(ArrayList<String>)fx[6];
		new CandidateCChallengeDbBuilder(new FixtureRunner())
			.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		final String artifactKey=CandidateCUtil.artifactKey(orderedProfileIds.get(0));
		final String challengeKey=CandidateCUtil.challengeKey(artifactKey);
		final String finalDir=sidecarDir+"/"+challengeKey+".hmmdb.d";
		final String dbPath=finalDir+"/db.hmmdb";
		final String originalHash=ArtifactRoleManifestGenerator.streamingSha256Hex(dbPath);
		final byte[] originalBytes=Files.readAllBytes(Paths.get(dbPath));
		//Append trailing bytes only -- every NAME line, count, and order stays exactly as published.
		final byte[] tamperedBytes=java.util.Arrays.copyOf(originalBytes, originalBytes.length+8);
		for(int i=originalBytes.length; i<tamperedBytes.length; i++){tamperedBytes[i]=(byte)'\n';}
		Files.write(Paths.get(dbPath), tamperedBytes);
		final String tamperedHash=ArtifactRoleManifestGenerator.streamingSha256Hex(dbPath);
		final String oldSidecar=readWholeFile(finalDir+"/complete.sidecar");
		final String rewrittenSidecar=oldSidecar.replace("db_hmmdb_sha256\t"+originalHash,
			"db_hmmdb_sha256\t"+tamperedHash);
		if(rewrittenSidecar.equals(oldSidecar)){
			throw new RuntimeException("SELFTEST SETUP FAILED: db_hmmdb_sha256 replacement did not match "
				+"the sidecar text.");
		}
		writeFilePlain(finalDir+"/complete.sidecar", rewrittenSidecar);
		final FixtureRunner runner2=new FixtureRunner();
		final Result r2=new CandidateCChallengeDbBuilder(runner2)
			.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		if(r2.resumed){
			throw new RuntimeException("SELFTEST FAILED: a tampered db.hmmdb whose sidecar hash was edited "
				+"to match it still resumed -- the mandatory four-way freshness chain did not fire.");
		}
		if(runner2.callCount!=1){
			throw new RuntimeException("SELFTEST FAILED: rebuild after the four-way-chain violation "
				+"invoked the runner "+runner2.callCount+" times, expected 1.");
		}
	}

	/** B2 negatives (Elly/UMP45 {@code d2e1850}): an artifact-manifest row's {@code artifact_key} column
	 * is a COPIED value and must never be trusted -- neither a string that fails the safe pattern nor a
	 * pattern-valid-but-wrong-derived key may be used to build the NAME-line translation map. */
	private static void testHostileArtifactManifestKeyCrashes(final String dir) throws Exception{
		final String work=dir+"/hostilekey";
		final Object[] fx=buildValidFixture(work, "fZ", 0);
		final String roleFile=(String)fx[0], neighborsFile=(String)fx[1], artifactManifestFile=(String)fx[2],
			sidecarDir=(String)fx[3], focalFamilyId=(String)fx[4];
		final int foldId=(Integer)fx[5];
		@SuppressWarnings("unchecked")
		final ArrayList<String> orderedProfileIds=(ArrayList<String>)fx[6];
		final String focalKey=CandidateCUtil.artifactKey(orderedProfileIds.get(0));
		final String oldText=readWholeFile(artifactManifestFile);
		writeFilePlain(artifactManifestFile, oldText.replace(focalKey, "../../etc/passwd"));
		boolean crashed=false;
		try{
			new CandidateCChallengeDbBuilder(new FixtureRunner())
				.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("does not match the required safe pattern")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a non-pattern "
					+"artifact_key: "+e.getMessage());
			}
		}
		if(!crashed){
			throw new RuntimeException("SELFTEST FAILED: a non-pattern artifact_key in the artifact "
				+"manifest did not crash.");
		}

		final String work2=dir+"/wrongderivedkey";
		final Object[] fx2=buildValidFixture(work2, "fZ", 0);
		final String roleFile2=(String)fx2[0], neighborsFile2=(String)fx2[1], artifactManifestFile2=(String)fx2[2],
			sidecarDir2=(String)fx2[3], focalFamilyId2=(String)fx2[4];
		final int foldId2=(Integer)fx2[5];
		@SuppressWarnings("unchecked")
		final ArrayList<String> orderedProfileIds2=(ArrayList<String>)fx2[6];
		final String focalKey2=CandidateCUtil.artifactKey(orderedProfileIds2.get(0));
		final String oldText2=readWholeFile(artifactManifestFile2);
		//Pattern-valid but NOT the real FNV-1a-64 hash of the focal's own artifact_id.
		writeFilePlain(artifactManifestFile2, oldText2.replace(focalKey2, "p_0000000000000000"));
		boolean crashed2=false;
		try{
			new CandidateCChallengeDbBuilder(new FixtureRunner())
				.process(roleFile2, neighborsFile2, artifactManifestFile2, sidecarDir2, focalFamilyId2, foldId2);
		}catch(RuntimeException e){
			crashed2=true;
			if(!e.getMessage().contains("independently-recomputed FNV key")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a wrong-derived "
					+"artifact_key: "+e.getMessage());
			}
		}
		if(!crashed2){
			throw new RuntimeException("SELFTEST FAILED: a wrong-derived artifact_key in the artifact "
				+"manifest did not crash.");
		}
		//NOTE: a live DUPLICATE-key negative (two distinct artifact_ids sharing one artifact_key) is not
		//constructible here -- the recompute-equality check just proven above requires EVERY row's
		//artifact_key to independently equal CandidateCUtil.artifactKey(itsOwnId), so two rows can only
		//collide via a genuine FNV-1a-64 hash collision (~2^64 search space). The keyToId.put(...)!=null
		//bijectivity check a few lines above this method's call site is retained as verified-by-inspection
		//defense-in-depth per the sealed plan's bijectivity requirement, not exercised by a fixture here.
	}

	/** B2 negative: a live {@code .hmm.canonical} input truncated to empty must crash loud on its own,
	 * distinct from (and checked before) the hash-mismatch case. */
	private static void testEmptyLiveHmmCanonicalCrashes(final String dir) throws Exception{
		final String work=dir+"/emptylivehmm";
		final Object[] fx=buildValidFixture(work, "fZ", 0);
		final String roleFile=(String)fx[0], neighborsFile=(String)fx[1], artifactManifestFile=(String)fx[2],
			sidecarDir=(String)fx[3], focalFamilyId=(String)fx[4];
		final int foldId=(Integer)fx[5];
		@SuppressWarnings("unchecked")
		final ArrayList<String> orderedProfileIds=(ArrayList<String>)fx[6];
		final String hmmDir=work+"/hmms";
		final String neighborKey=CandidateCUtil.artifactKey(orderedProfileIds.get(3));
		writeFilePlain(hmmDir+"/"+neighborKey+".hmm.canonical", "");
		boolean crashed=false;
		try{
			new CandidateCChallengeDbBuilder(new FixtureRunner())
				.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		}catch(RuntimeException e){
			crashed=true;
			if(!e.getMessage().contains("missing or empty")){
				throw new RuntimeException("SELFTEST FAILED: wrong crash reason for an empty live "
					+".hmm.canonical: "+e.getMessage());
			}
		}
		if(!crashed){throw new RuntimeException("SELFTEST FAILED: an empty live .hmm.canonical did not crash.");}
	}

	/** B3 negative: a sidecar recording a different {@code hmmpress_version} than this class's own
	 * pinned constant must force a rebuild, not a crash -- a resume token exactly like Phase A's
	 * mafft/hmmbuild versions. */
	private static void testWrongHmmpressVersionRebuilds(final String dir) throws Exception{
		final String work=dir+"/wrongpressversion";
		final Object[] fx=buildValidFixture(work, "fZ", 0);
		final String roleFile=(String)fx[0], neighborsFile=(String)fx[1], artifactManifestFile=(String)fx[2],
			sidecarDir=(String)fx[3], focalFamilyId=(String)fx[4];
		final int foldId=(Integer)fx[5];
		@SuppressWarnings("unchecked")
		final ArrayList<String> orderedProfileIds=(ArrayList<String>)fx[6];
		new CandidateCChallengeDbBuilder(new FixtureRunner())
			.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		final String artifactKey=CandidateCUtil.artifactKey(orderedProfileIds.get(0));
		final String challengeKey=CandidateCUtil.challengeKey(artifactKey);
		final String finalDir=sidecarDir+"/"+challengeKey+".hmmdb.d";
		final String oldSidecar=readWholeFile(finalDir+"/complete.sidecar");
		final String rewritten=oldSidecar.replace("hmmpress_version\t"+HMMPRESS_VERSION,
			"hmmpress_version\tHMMER 3.1b2 (Feb 2015)");
		if(rewritten.equals(oldSidecar)){
			throw new RuntimeException("SELFTEST SETUP FAILED: hmmpress_version replacement did not match "
				+"the sidecar text.");
		}
		writeFilePlain(finalDir+"/complete.sidecar", rewritten);
		final FixtureRunner runner2=new FixtureRunner();
		final Result r2=new CandidateCChallengeDbBuilder(runner2)
			.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		if(r2.resumed){
			throw new RuntimeException("SELFTEST FAILED: a sidecar with a wrong hmmpress_version still resumed.");
		}
		if(runner2.callCount!=1){
			throw new RuntimeException("SELFTEST FAILED: rebuild after a wrong hmmpress_version invoked "
				+"the runner "+runner2.callCount+" times, expected 1.");
		}
	}

	/** B4 negative: cleanup must delete a symlink as itself, never traverse into its target -- a
	 * malicious/stray link inside a staging or stale-final directory must not cause deletion outside
	 * the directory being cleaned up. */
	private static void testDeleteRecursiveNeverFollowsSymlinks(final String dir) throws Exception{
		final String work=dir+"/symlinknofollow";
		new File(work).mkdirs();
		final String externalDir=work+"/external";
		new File(externalDir).mkdirs();
		final String externalFile=externalDir+"/protected.txt";
		writeFilePlain(externalFile, "DO_NOT_DELETE");
		final String victimDir=work+"/victim";
		new File(victimDir).mkdirs();
		Files.createSymbolicLink(Paths.get(victimDir+"/link_to_external"), Paths.get(externalDir));
		deleteRecursive(new File(victimDir));
		if(new File(victimDir).exists()){
			throw new RuntimeException("SELFTEST FAILED: deleteRecursive did not remove the victim directory.");
		}
		if(!new File(externalFile).exists()){
			throw new RuntimeException("SELFTEST FAILED: deleteRecursive followed a symlink and deleted "
				+"content OUTSIDE the directory it was asked to clean up.");
		}
	}

	/** B4 negative: cleanup must fail loud (never silently leave a partial tree behind) when a deletion
	 * does not fully complete -- simulated by revoking write permission on the containing directory,
	 * which makes POSIX unlink() of its child fail regardless of the child's own permissions. */
	private static void testDeleteRecursiveFailsLoudOnCleanupFailure(final String dir) throws Exception{
		final String work=dir+"/cleanupfail";
		new File(work).mkdirs();
		final String blockedDir=work+"/blocked";
		new File(blockedDir).mkdirs();
		writeFilePlain(blockedDir+"/child.txt", "x");
		final File blockedDirFile=new File(blockedDir);
		if(!blockedDirFile.setWritable(false)){
			System.err.println("    (B4 cleanup-failure negative SKIPPED: this filesystem does not "
				+"support revoking directory write permission.)");
			return;
		}
		try{
			boolean threw=false;
			try{
				deleteRecursive(blockedDirFile);
			}catch(RuntimeException e){
				threw=true;
				if(!e.getMessage().contains("Could not delete") && !e.getMessage().contains("still exists")){
					throw new RuntimeException("SELFTEST FAILED: wrong crash reason for a cleanup "
						+"failure: "+e.getMessage());
				}
			}
			if(!threw){
				System.err.println("    (B4 cleanup-failure negative SKIPPED: deletion succeeded despite "
					+"revoked permission -- likely running with elevated privileges.)");
			}
		}finally{
			blockedDirFile.setWritable(true);
			deleteRecursive(blockedDirFile);
		}
	}

	/** B5 negative: two separate rebuilds of the SAME challenge must never reuse or collide on a staging
	 * directory, and neither generated name may look like a final {@code .hmmdb.d} directory (so B2's
	 * top-level enumeration never miscounts it). */
	private static void testOverlappingAttemptsGetDifferentStagingDirectories(final String dir) throws Exception{
		final String work=dir+"/overlappingstaging";
		final Object[] fx=buildValidFixture(work, "fZ", 0);
		final String roleFile=(String)fx[0], neighborsFile=(String)fx[1], artifactManifestFile=(String)fx[2],
			sidecarDir=(String)fx[3], focalFamilyId=(String)fx[4];
		final int foldId=(Integer)fx[5];
		@SuppressWarnings("unchecked")
		final ArrayList<String> orderedProfileIds=(ArrayList<String>)fx[6];
		final FixtureRunner runner1=new FixtureRunner();
		new CandidateCChallengeDbBuilder(runner1)
			.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		final String staging1=runner1.workingDirLog.get(0);
		final String artifactKey=CandidateCUtil.artifactKey(orderedProfileIds.get(0));
		final String challengeKey=CandidateCUtil.challengeKey(artifactKey);
		final String finalDir=sidecarDir+"/"+challengeKey+".hmmdb.d";
		writeFilePlain(finalDir+"/db.hmmdb", "CORRUPTED_TO_FORCE_REBUILD");
		final FixtureRunner runner2=new FixtureRunner();
		new CandidateCChallengeDbBuilder(runner2)
			.process(roleFile, neighborsFile, artifactManifestFile, sidecarDir, focalFamilyId, foldId);
		final String staging2=runner2.workingDirLog.get(0);
		if(staging1.equals(staging2)){
			throw new RuntimeException("SELFTEST FAILED: two separate rebuilds of the same challenge "
				+"reused the identical staging directory '"+staging1+"'.");
		}
		final String name1=new File(staging1).getName(), name2=new File(staging2).getName();
		if(!name1.startsWith(challengeKey+".hmmdb.building.") || !name2.startsWith(challengeKey+".hmmdb.building.")){
			throw new RuntimeException("SELFTEST FAILED: staging directory names ('"+name1+"', '"+name2
				+"') do not carry the expected prefix.");
		}
		if(name1.endsWith(".hmmdb.d") || name2.endsWith(".hmmdb.d")){
			throw new RuntimeException("SELFTEST FAILED: a staging directory name looked like a final "
				+".hmmdb.d directory: '"+name1+"' / '"+name2+"'.");
		}
	}
}
