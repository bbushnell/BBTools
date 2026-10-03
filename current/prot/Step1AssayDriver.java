package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import structures.ByteBuilder;

/**
 * Step-1 assay harness DRIVER for items 1 (exhaustive top-1/confusion) and 3a (positive ambiguity
 * margin/pass-fraction) only (plans/STEP1_ASSAY_HARNESS_DESIGN_v1.md #3/#4.1/#4.3/#7). Item 2
 * (shortlist recall) is explicitly OUT OF SCOPE here -- it is a named NOT_YET_COMPUTABLE gap
 * (blocked on Stage 1B seed-proxy machinery), not an ambiguity to resolve. Item 3b (negative-query
 * separation) is deferred too: it needs the same unresolved nearest-training-fold-member identity
 * stratification bins (design doc #5/#6) that 1/3a's own math does not depend on -- a real Brian
 * methodology decision, not something to invent here.
 * <p>
 * Reuses the already-sealed {@link Step1AssayMetrics} (classification/margin-support primitives)
 * and {@link FamilyArtifactScore} (margin/clearsMargin/flatScore) without reimplementing either.
 * Every I/O loader below is independently reimplemented per this project's established convention
 * -- it does not call {@link CandidateCChallengeDbBuilder}'s or {@link NegativeQueryManifestGenerator}'s
 * loaders, even though the shapes are the same (role-universe/neighbor resolution, group->member
 * expansion, needed-IDs-only FASTA extraction).
 * <p>
 * <b>Method-agnostic dispatch (design doc #2):</b> {@link ArtifactScorer#scoreChallenge} is
 * BATCH-oriented -- one call scores every held-out query in a focal/fold against every one of that
 * fold's 21 challenge-set artifacts at once, so a future HMM-based scorer can shell {@code hmmscan}
 * exactly ONCE per focal/fold batch rather than once per (query,artifact) pair. {@link
 * FlatConsensusScorer} is the one scorer implemented here (Candidate A, {@code fasta_consensus});
 * it loops in memory since it needs no external process.
 *
 * @author Eru
 */
public final class Step1AssayDriver {

	/** One held-out positive query: its member id and RAW (unencoded) ASCII protein sequence
	 * bytes -- method-agnostic; each {@link ArtifactScorer} does its own encoding internally
	 * (e.g. {@link FlatConsensusScorer} BLOSUM62-encodes via {@link ConsensusRepBuilder#encodeLenient},
	 * a future HMM scorer would write these bytes straight into a FASTA for {@code hmmscan}). */
	public static final class Query {
		public final String memberId;
		public final byte[] sequence;
		public Query(final String memberId, final byte[] sequence){
			this.memberId=memberId; this.sequence=sequence;
		}
	}

	/** Method-agnostic scoring dispatch (design doc #2). */
	public interface ArtifactScorer {
		/**
		 * Scores every query in {@code queries} against every artifact in {@code orderedArtifactIds}
		 * in ONE batched call -- {@code result[q][a]} is the score of {@code queries.get(q)} against
		 * {@code orderedArtifactIds.get(a)}, in the scorer's own native units (BSR points for A,
		 * bits for C -- never compared across methods). {@code focalFamilyId}/{@code foldId} are
		 * passed through so a real Candidate-C implementation can locate that focal/fold's pressed
		 * challenge database.
		 */
		double[][] scoreChallenge(ArrayList<Query> queries, ArrayList<String> orderedArtifactIds,
			String focalFamilyId, int foldId);

		/**
		 * Every file path this scorer instance will read to answer {@link #scoreChallenge} calls --
		 * runAssay hashes each one at run-start and re-hashes at run-end before publishing output,
		 * rejecting any drift (CP-11). Must be stable and COMPLETE before the first {@link
		 * #scoreChallenge} call (e.g. every {@code artifact_manifest.tsv} row's path for {@link
		 * FlatConsensusScorer}; a future Candidate-C scorer would return its pressed profile-database
		 * file(s) here too, per Elly's design requirement, 2026-09-02). Default returns an empty list
		 * (a synthetic/canned test scorer consumes no real files) -- a real scorer MUST override this.
		 */
		default ArrayList<String> consumedFilePaths(){return new ArrayList<String>();}
	}

	/** One held-out positive query's outcome. Retains the raw {@code sOwn}/{@code sMaxOther}
	 * (not merely the derived margin) so every §7 output row is independently re-derivable. */
	public static final class QueryResult {
		public final String memberId, focalFamilyId, cellId;
		public final int foldId;
		public final Step1AssayMetrics.Top1Bucket bucket;
		public final double sOwn, sMaxOther, marginNative;
		public final boolean passedMargin;
		QueryResult(String memberId, String focalFamilyId, int foldId, String cellId,
				Step1AssayMetrics.Top1Bucket bucket, double sOwn, double sMaxOther,
				double marginNative, boolean passedMargin){
			this.memberId=memberId; this.focalFamilyId=focalFamilyId; this.foldId=foldId; this.cellId=cellId;
			this.bucket=bucket; this.sOwn=sOwn; this.sMaxOther=sMaxOther;
			this.marginNative=marginNative; this.passedMargin=passedMargin;
		}
	}

	/** One family's §4.1/§4.3 family-first aggregate, pooling ALL its held-out queries ACROSS ALL
	 * its folds. */
	public static final class FamilySummary {
		public final String familyId, cellId;
		public final int totalHeldOut, correct, confused, tied;
		public final double selfConsistency, confusionMass, tiedRate, passFraction;
		FamilySummary(String familyId, String cellId, int totalHeldOut, int correct, int confused,
				int tied, double selfConsistency, double confusionMass, double tiedRate, double passFraction){
			this.familyId=familyId; this.cellId=cellId; this.totalHeldOut=totalHeldOut;
			this.correct=correct; this.confused=confused; this.tied=tied;
			this.selfConsistency=selfConsistency; this.confusionMass=confusionMass;
			this.tiedRate=tiedRate; this.passFraction=passFraction;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Pure scoring/classification (no I/O)     */
	/*--------------------------------------------------------------*/

	/**
	 * Scores one focal family's one fold: calls {@code scorer} ONCE (the batch seam), validates the
	 * returned matrix is exactly {@code heldOutQueries.size() x challengeSetArtifactIds.size()} and
	 * every entry finite (Elly's pin) BEFORE any metric call, then classifies/margins each row via
	 * the already-sealed {@link Step1AssayMetrics}/{@link FamilyArtifactScore} primitives. The focal
	 * artifact is always index 0 of {@code challengeSetArtifactIds} by construction (focal-then-
	 * neighbor-rank order, same convention as {@link CandidateCChallengeDbBuilder}).
	 */
	public static ArrayList<QueryResult> scoreFocalFold(final String focalFamilyId, final int foldId,
			final String cellId, final ArrayList<String> challengeSetArtifactIds,
			final ArrayList<Query> heldOutQueries, final ArtifactScorer scorer, final double marginMin){
		if(heldOutQueries.isEmpty()){
			throw new RuntimeException("Focal '"+focalFamilyId+"' fold "+foldId+" has zero held-out "
				+"queries -- a real fold/family must have at least one, this is a fixture/manifest bug.");
		}
		//CP-08: the challenge set's OWN size must be exactly 21 regardless of what the scorer returns --
		//a generic matrix-width check alone would accept a self-consistent-but-wrong-size list (Sayu/CP-08).
		if(challengeSetArtifactIds.size()!=21){
			throw new RuntimeException("Focal '"+focalFamilyId+"' fold "+foldId+": challenge set has "
				+challengeSetArtifactIds.size()+" artifacts, expected exactly 21 (1 focal + 20 neighbor decoys).");
		}
		final double[][] matrix=scorer.scoreChallenge(heldOutQueries, challengeSetArtifactIds,
			focalFamilyId, foldId);
		final int nQ=heldOutQueries.size(), nA=challengeSetArtifactIds.size();
		if(matrix.length!=nQ){
			throw new RuntimeException("Focal '"+focalFamilyId+"' fold "+foldId+": scorer returned "
				+matrix.length+" score rows, expected exactly "+nQ+" (one per held-out query).");
		}
		for(int i=0; i<nQ; i++){
			if(matrix[i].length!=nA){
				throw new RuntimeException("Focal '"+focalFamilyId+"' fold "+foldId+": score row "+i
					+" has "+matrix[i].length+" columns, expected exactly "+nA+" (one per challenge-set "
					+"artifact).");
			}
			for(int j=0; j<nA; j++){
				if(!Double.isFinite(matrix[i][j])){
					throw new RuntimeException("Focal '"+focalFamilyId+"' fold "+foldId+": score["+i
						+"]["+j+"]="+matrix[i][j]+" is not finite -- a missing/failed score must never "
						+"reach classification (see Step1AssayMetrics.requireFinite's identical contract).");
				}
			}
		}
		final int ownIndex=0;//focal-then-neighbor-rank order, focal is always rank 0.
		final ArrayList<QueryResult> results=new ArrayList<QueryResult>(nQ);
		for(int i=0; i<nQ; i++){
			final double[] scores=matrix[i];
			final Step1AssayMetrics.Top1Bucket bucket=Step1AssayMetrics.classifyTop1(scores, ownIndex);
			final double sOwn=scores[ownIndex];
			final double sMaxOther=Step1AssayMetrics.maxOtherScore(scores, ownIndex);
			final double margin=FamilyArtifactScore.margin(sOwn, sMaxOther);
			final boolean passed=FamilyArtifactScore.clearsMargin(sOwn, sMaxOther, marginMin);
			results.add(new QueryResult(heldOutQueries.get(i).memberId, focalFamilyId, foldId, cellId,
				bucket, sOwn, sMaxOther, margin, passed));
		}
		return results;
	}

	/**
	 * Pools ALL of one family's held-out queries ACROSS ALL its folds (design doc #4.1/#4.3
	 * family-first rule) into one summary. {@code self_consistency}/{@code confusion_mass}/
	 * {@code tied_rate} are independent fractions over the SAME {@code total_held_out} denominator
	 * (never each other's complement); the exhaustiveness invariant is a hard assert, not a hope.
	 */
	public static FamilySummary aggregateFamily(final String familyId, final String cellId,
			final ArrayList<QueryResult> allResultsAcrossThisFamilysFolds, final HashSet<Integer> expectedFoldIds){
		if(allResultsAcrossThisFamilysFolds.isEmpty()){
			throw new RuntimeException("Family '"+familyId+"' has zero held-out query results across "
				+"all its folds -- cannot aggregate.");
		}
		int correct=0, confused=0, tied=0, passed=0;
		final HashSet<String> seenMemberIds=new HashSet<String>();
		final HashSet<Integer> actualFoldIds=new HashSet<Integer>();
		for(final QueryResult r : allResultsAcrossThisFamilysFolds){
			if(!r.focalFamilyId.equals(familyId)){
				throw new RuntimeException("aggregateFamily('"+familyId+"') was given a QueryResult "
					+"for family '"+r.focalFamilyId+"' -- caller must pre-filter to one family's own results.");
			}
			//CP-04: every result must share the ONE supplied cellId -- a mixed-cell aggregate would
			//silently blend two different stratification cells into one summary.
			if(!r.cellId.equals(cellId)){
				throw new RuntimeException("aggregateFamily('"+familyId+"', cellId='"+cellId+"') was given a "
					+"QueryResult for member '"+r.memberId+"' with cellId='"+r.cellId+"' -- mixed-cell "
					+"aggregation is forbidden.");
			}
			//CP-04: a member id appearing more than once across this family's own held-out set (any
			//fold) signals fold_manifest/cluster.tsv corruption or a caller bug -- never silently
			//double-count it into the denominator.
			if(!seenMemberIds.add(r.memberId)){
				throw new RuntimeException("aggregateFamily('"+familyId+"'): member '"+r.memberId
					+"' appears more than once across this family's held-out results -- each held-out "
					+"member must be scored exactly once, in exactly one fold.");
			}
			actualFoldIds.add(r.foldId);
			switch(r.bucket){
				case CORRECT: correct++; break;
				case CONFUSED: confused++; break;
				case TIED: tied++; break;
				default: throw new RuntimeException("Unknown Top1Bucket: "+r.bucket);
			}
			if(r.passedMargin){passed++;}
		}
		//CP-04 (Elly's addition): the family must have covered EXACTLY its expected fold set from
		//fold_manifest.tsv -- neither a missing fold (silently under-scored) nor an unexpected extra
		//fold (a caller bug feeding in the wrong family/fold's results).
		if(!actualFoldIds.equals(expectedFoldIds)){
			throw new RuntimeException("aggregateFamily('"+familyId+"'): actual fold coverage "+actualFoldIds
				+" != expected fold coverage "+expectedFoldIds+" (from fold_manifest.tsv) -- every expected "
				+"fold must be scored exactly once, no more, no fewer.");
		}
		final int total=allResultsAcrossThisFamilysFolds.size();
		assert(correct+confused+tied==total) : "Family '"+familyId+"': correct("+correct+")+confused("
			+confused+")+tied("+tied+") != total_held_out("+total+") -- the three-bucket partition "
			+"must be exhaustive by construction (Step1AssayMetrics.classifyTop1's own contract).";
		if(correct+confused+tied!=total){
			throw new RuntimeException("Family '"+familyId+"': correct+confused+tied ("+(correct+confused+tied)
				+") != total_held_out ("+total+") -- exhaustiveness invariant violated.");
		}
		return new FamilySummary(familyId, cellId, total, correct, confused, tied,
			correct/(double)total, confused/(double)total, tied/(double)total, passed/(double)total);
	}

	/*--------------------------------------------------------------*/
	/*----------------      Top-level orchestration           -------*/
	/*--------------------------------------------------------------*/

	/**
	 * Items 1/3a end-to-end: iterates every {@code sample_manifest.tsv} {@code selected=1} family
	 * crossed with EXACTLY that family's own {@code fold_manifest.tsv} fold_id set (Elly's addition,
	 * 2026-09-01 -- proven by {@link #aggregateFamily}'s expected-fold-coverage assert), always
	 * binding the selected row's OWN {@code cluster_path} (never a separately-supplied one, closing
	 * CP-01 by construction -- there is no such parameter here), validates the resolved 21-item
	 * challenge set at the orchestration boundary (CP-08/CP-09), scores/aggregates per family, then
	 * writes both §7 output files with real provenance headers -- failing loud on any consumed-file
	 * hash drift between capture and publication (CP-11) before either file is written.
	 */
	public static void runAssay(final String sampleManifestFile, final String foldManifestFile,
			final String roleFile, final String neighborsFile, final String artifactManifestFile,
			final String memberSeqsFasta, final ArtifactScorer scorer, final String method,
			final double marginMin, final String scorerDescription, final String top1OutFile,
			final String marginOutFile){
		if(scorerDescription==null || scorerDescription.trim().isEmpty()){
			throw new RuntimeException("runAssay: scorerDescription must be non-empty -- provenance "
				+"requires a real scorer name/version/parameters description, never a placeholder.");
		}
		final LinkedHashMap<String,String> fixedInputs=new LinkedHashMap<String,String>();
		fixedInputs.put("sample_manifest", sampleManifestFile);
		fixedInputs.put("fold_manifest", foldManifestFile);
		fixedInputs.put("role_manifest", roleFile);
		fixedInputs.put("neighbors", neighborsFile);
		fixedInputs.put("artifact_manifest", artifactManifestFile);
		fixedInputs.put("member_seqs_fasta", memberSeqsFasta);
		final LinkedHashMap<String,String> startHashes=new LinkedHashMap<String,String>();
		for(final Map.Entry<String,String> e : fixedInputs.entrySet()){
			startHashes.put(e.getKey(), sha256File(e.getValue()));
		}

		//Elly's design requirement, 2026-09-02: provenance must tie to the actual bytes the SCORER
		//consumes, not merely to a manifest entry that might diverge from what's really read. The
		//scorer reports its own complete, stable file set BEFORE any scoring happens.
		final ArrayList<String> scorerPaths=scorer.consumedFilePaths();
		final LinkedHashMap<String,String> scorerFileHashesStart=new LinkedHashMap<String,String>();
		for(final String p : scorerPaths){
			validateSafePath(p, "scorer-consumed file");
			if(!scorerFileHashesStart.containsKey(p)){scorerFileHashesStart.put(p, sha256File(p));}
		}

		final ArrayList<FamilySelectionRow> selected=loadAllSelectedFamilies(sampleManifestFile);
		if(selected.isEmpty()){
			throw new RuntimeException("No selected=1 rows in "+sampleManifestFile+" -- nothing to assay.");
		}
		final LinkedHashMap<String,String> clusterHashesStart=new LinkedHashMap<String,String>();
		final ArrayList<FamilySummary> summaries=new ArrayList<FamilySummary>();
		final ArrayList<QueryResult> allResults=new ArrayList<QueryResult>();
		long expectedCells=0, observedCells=0;

		for(final FamilySelectionRow row : selected){
			//Elly's addition: validate the selected row's OWN cluster_path before hashing/expansion.
			validateSafePath(row.clusterTsvPath, "sample_manifest "+sampleManifestFile+" family '"+row.repId+"'");
			if(!clusterHashesStart.containsKey(row.clusterTsvPath)){
				clusterHashesStart.put(row.clusterTsvPath, sha256File(row.clusterTsvPath));
			}
			final HashSet<Integer> expectedFoldIds=loadExpectedFoldIds(foldManifestFile, row.repId);
			//Elly's finding, 2026-09-02: HashSet iteration order is NOT a language guarantee -- process
			//folds in a canonical SORTED order so output row order is deterministic regardless of
			//fold_manifest.tsv's file order or JVM hashing behavior. expectedFoldIds itself stays a
			//HashSet for aggregateFamily's order-independent equality check.
			final ArrayList<Integer> foldProcessingOrder=new ArrayList<Integer>(expectedFoldIds);
			java.util.Collections.sort(foldProcessingOrder);
			final ArrayList<QueryResult> familyResults=new ArrayList<QueryResult>();
			for(final int foldId : foldProcessingOrder){
				final ArrayList<String> challengeSet=loadChallengeSetIds(roleFile, neighborsFile, row.repId, foldId);
				final String independentFocalId=findFocalArtifactIdFromFile(roleFile, row.repId, foldId);
				validateChallengeSetBoundary(challengeSet, independentFocalId, row.repId, foldId);
				//CP-01: always the selected row's OWN cluster_path -- no other path is reachable here.
				final ArrayList<String> heldOutIds=loadHeldOutMemberIds(row.clusterTsvPath, foldManifestFile,
					row.repId, foldId);
				final HashMap<String,byte[]> seqs=extractSequences(memberSeqsFasta, new HashSet<String>(heldOutIds));
				final ArrayList<Query> queries=new ArrayList<Query>(heldOutIds.size());
				for(final String id : heldOutIds){queries.add(new Query(id, seqs.get(id)));}
				expectedCells+=(long)heldOutIds.size()*challengeSet.size();
				final ArrayList<QueryResult> foldResults=scoreFocalFold(row.repId, foldId, row.cellId, challengeSet,
					queries, scorer, marginMin);
				observedCells+=(long)foldResults.size()*challengeSet.size();
				familyResults.addAll(foldResults);
			}
			final FamilySummary summary=aggregateFamily(row.repId, row.cellId, familyResults, expectedFoldIds);
			summaries.add(summary);
			allResults.addAll(familyResults);
		}

		//CP-11: reject any drift in a consumed file's bytes between capture and publication.
		assertNoHashDrift(fixedInputs, startHashes, "Upstream input");
		assertNoPathHashDrift(clusterHashesStart, "cluster.tsv");
		assertNoPathHashDrift(scorerFileHashesStart, "scorer-consumed file");

		final ProvenanceHeader prov=new ProvenanceHeader(fixedInputs, startHashes, clusterHashesStart,
			scorerFileHashesStart, scorerDescription, expectedCells, observedCells);
		writeTop1ConfusionTsv(top1OutFile, summaries, method, prov);
		writeMargin3aTsv(marginOutFile, allResults, summaries, method, marginMin, prov);
	}

	/*--------------------------------------------------------------*/
	/*----------------              Loaders                ----------*/
	/*--------------------------------------------------------------*/

	static final class FamilySelectionRow { String repId, cellId, clusterTsvPath; }

	/** Parses sample_manifest.tsv by HEADER NAME (schema evolves; convention matches
	 * FoldManifestAggregator.loadSelected), requiring EXACTLY ONE row with rep_id==focalFamilyId AND
	 * selected==1 (Elly's pin) -- zero or more than one is a hard crash, never a silent pick. */
	static FamilySelectionRow findSelectedFamilyRow(final String sampleManifestFile, final String focalFamilyId){
		final ArrayList<FamilySelectionRow> matches=new ArrayList<FamilySelectionRow>();
		final ByteFile bf=ByteFile.makeByteFile(sampleManifestFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		HashMap<String,Integer> col=null;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='#'){
				if(col==null){
					final String h=new String(line, 1, line.length-1);
					final String[] names=h.split("\t");
					//The real header line is the ONE with more than one column and no digit-only
					//content -- #shortfall lines are also '#'-prefixed (SampleManifestGenerator's own
					//convention); skip any '#' line that doesn't look like the header itself.
					if(names.length>1 && names[0].trim().equalsIgnoreCase("rep_id")){
						col=new HashMap<String,Integer>();
						for(int i=0; i<names.length; i++){col.put(names[i].trim().toLowerCase(), i);}
					}
				}
				continue;
			}
			if(col==null){throw new RuntimeException("Data row before header in "+sampleManifestFile);}
			lp.set(line);
			final String repId=lp.parseString(idx(col, "rep_id"));
			if(!repId.equals(focalFamilyId)){continue;}
			final int selected=lp.parseInt(idx(col, "selected"));
			if(selected!=1){continue;}
			final FamilySelectionRow row=new FamilySelectionRow();
			row.repId=repId;
			row.cellId=lp.parseString(idx(col, "cell_id"));
			row.clusterTsvPath=lp.parseString(idx(col, "cluster_path"));
			matches.add(row);
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+sampleManifestFile);}
		if(matches.isEmpty()){
			throw new RuntimeException("No selected=1 row for rep_id='"+focalFamilyId+"' in "
				+sampleManifestFile+" -- exactly one is required.");
		}
		if(matches.size()>1){
			throw new RuntimeException(sampleManifestFile+" has "+matches.size()+" selected=1 rows for "
				+"rep_id='"+focalFamilyId+"' -- exactly one is required, the manifest is corrupt.");
		}
		return matches.get(0);
	}

	/** Every selected=1 row across the WHOLE sample_manifest.tsv (not just one family) -- runAssay's
	 * top-level iteration source. Rejects more than one selected=1 row for the same rep_id (manifest
	 * corruption), mirroring {@link #findSelectedFamilyRow}'s exact parsing convention. */
	static ArrayList<FamilySelectionRow> loadAllSelectedFamilies(final String sampleManifestFile){
		final ArrayList<FamilySelectionRow> result=new ArrayList<FamilySelectionRow>();
		final HashSet<String> seenRepIds=new HashSet<String>();
		final ByteFile bf=ByteFile.makeByteFile(sampleManifestFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		HashMap<String,Integer> col=null;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='#'){
				if(col==null){
					final String h=new String(line, 1, line.length-1);
					final String[] names=h.split("\t");
					if(names.length>1 && names[0].trim().equalsIgnoreCase("rep_id")){
						col=new HashMap<String,Integer>();
						for(int i=0; i<names.length; i++){col.put(names[i].trim().toLowerCase(), i);}
					}
				}
				continue;
			}
			if(col==null){throw new RuntimeException("Data row before header in "+sampleManifestFile);}
			lp.set(line);
			final int selected=lp.parseInt(idx(col, "selected"));
			if(selected!=1){continue;}
			final String repId=lp.parseString(idx(col, "rep_id"));
			if(!seenRepIds.add(repId)){
				throw new RuntimeException(sampleManifestFile+" has more than one selected=1 row for rep_id='"
					+repId+"' -- exactly one is required per family.");
			}
			final FamilySelectionRow row=new FamilySelectionRow();
			row.repId=repId;
			row.cellId=lp.parseString(idx(col, "cell_id"));
			row.clusterTsvPath=lp.parseString(idx(col, "cluster_path"));
			result.add(row);
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+sampleManifestFile);}
		if(col==null){throw new RuntimeException("No rep_id header found in "+sampleManifestFile);}
		return result;
	}

	static int idx(final HashMap<String,Integer> col, final String name){
		final Integer i=col.get(name);
		if(i==null){throw new RuntimeException("Missing required column: "+name);}
		return i;
	}

	/** cell_id for one focal family, bound directly from its ONE selected sample_manifest.tsv row
	 * (Elly's pin -- never deferred to a later join). */
	public static String loadCellId(final String sampleManifestFile, final String focalFamilyId){
		return findSelectedFamilyRow(sampleManifestFile, focalFamilyId).cellId;
	}

	static class UniverseRow { String artifactId, familyId, role; int foldId; }

	/** Focal + 20 neighbor-decoy artifact ids, focal-then-neighbor-rank order -- independently
	 * reimplemented from (not calling) {@link CandidateCChallengeDbBuilder#resolveChallenge}, same
	 * validation rigor: the literal fold-local focal artifact_id is SELECTED from the role-manifest
	 * universe (role=focal_fold_local, family_id, fold_id all matching), never constructed by naming
	 * convention (Elly's mandatory pin). */
	public static ArrayList<String> loadChallengeSetIds(final String roleFile, final String neighborsFile,
			final String focalFamilyId, final int foldId){
		final HashMap<String,UniverseRow> byArtifactId=loadRoleUniverse(roleFile);
		final UniverseRow focalRow=findFocalRow(byArtifactId, roleFile, focalFamilyId, foldId);
		final ArrayList<String> neighborFamilies=loadNeighborsForFocal(neighborsFile, focalFamilyId);
		final ArrayList<String> orderedIds=new ArrayList<String>(neighborFamilies.size()+1);
		orderedIds.add(focalRow.artifactId);
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
			orderedIds.add(decoyId);
		}
		return orderedIds;
	}

	/** Extracted from {@link #loadChallengeSetIds} (pure refactor, output-preserving): finds the ONE
	 * focal_fold_local row for (family,fold) in an already-loaded role universe. */
	static UniverseRow findFocalRow(final HashMap<String,UniverseRow> byArtifactId, final String roleFile,
			final String focalFamilyId, final int foldId){
		UniverseRow focalRow=null;
		for(final UniverseRow u : byArtifactId.values()){
			if(u.role.equals("focal_fold_local") && u.familyId.equals(focalFamilyId) && u.foldId==foldId){
				if(focalRow!=null){
					throw new RuntimeException(roleFile+" has more than one focal_fold_local row for "
						+"family '"+focalFamilyId+"' fold "+foldId+".");
				}
				focalRow=u;
			}
		}
		if(focalRow==null){
			throw new RuntimeException(roleFile+" has no focal_fold_local row for family '"
				+focalFamilyId+"' fold "+foldId+".");
		}
		return focalRow;
	}

	/** Independently re-derives the literal focal artifact_id for (family,fold) by RE-READING the
	 * role manifest from disk -- a genuine second lookup, not the cached value {@link
	 * #loadChallengeSetIds} already computed, so a caller (runAssay) can assert a previously-built
	 * challenge set's index 0 against it without trusting its own earlier call (CP-08/CP-09). */
	static String findFocalArtifactIdFromFile(final String roleFile, final String focalFamilyId, final int foldId){
		return findFocalRow(loadRoleUniverse(roleFile), roleFile, focalFamilyId, foldId).artifactId;
	}

	/** CP-08/CP-09 orchestration-boundary guard: the challenge set must be exactly 21 items and its
	 * index 0 must equal the independently-resolved literal focal artifact id -- never caller-trusted. */
	static void validateChallengeSetBoundary(final ArrayList<String> challengeSet, final String independentFocalId,
			final String familyId, final int foldId){
		if(challengeSet.size()!=21){
			throw new RuntimeException("Family '"+familyId+"' fold "+foldId+": challenge set has "
				+challengeSet.size()+" artifacts, expected exactly 21.");
		}
		if(!challengeSet.get(0).equals(independentFocalId)){
			throw new RuntimeException("Family '"+familyId+"' fold "+foldId+": challenge set index 0 is '"
				+challengeSet.get(0)+"', but the independently-resolved focal artifact is '"+independentFocalId
				+"' -- refusing to score against a mismatched focal binding.");
		}
	}

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

	/** family_neighbors.tsv: {@code #k\t20} header, rows {@code focal_family\tneighbor_rank\t
	 * neighbor_family}, sorted by neighbor_rank ascending -- returns only the ONE requested focal's
	 * ordered neighbor list. */
	static ArrayList<String> loadNeighborsForFocal(final String file, final String focalFamilyId){
		final TreeMap<Integer,String> byRank=new TreeMap<Integer,String>();
		final HashSet<String> neighborFamiliesSeen=new HashSet<String>();
		int declaredK=-1;
		final ByteFile bf=ByteFile.makeByteFile(file, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='#'){
				final String s=new String(line);
				if(s.startsWith("#k\t")){
					if(declaredK>=0){throw new RuntimeException(file+" has more than one '#k\\t<K>' header line.");}
					declaredK=Integer.parseInt(s.substring(3).trim());
				}
				continue;
			}
			lp.set(line);
			if(lp.terms()<3){throw new RuntimeException("Malformed neighbor row in "+file+": "+new String(line));}
			final String focal=lp.parseString(0);
			if(!focal.equals(focalFamilyId)){continue;}
			final int rank=lp.parseInt(1);
			final String neighbor=lp.parseString(2);
			if(neighbor.equals(focal)){
				throw new RuntimeException("Focal family '"+focal+"' lists itself as its own neighbor in "+file+".");
			}
			if(byRank.put(rank, neighbor)!=null){
				throw new RuntimeException("Focal family '"+focal+"' has neighbor_rank "+rank
					+" recorded more than once in "+file+".");
			}
			//CP-02: the SAME neighbor family at two different ranks would yield a duplicate decoy id in
			//the 21-item challenge set -- reject at parse time, mirroring CandidateCChallengeDbBuilder
			//.loadNeighbors's own neighborsSeen guard (that file, lines ~602-605).
			if(!neighborFamiliesSeen.add(neighbor)){
				throw new RuntimeException("Focal family '"+focal+"' has neighbor family '"+neighbor
					+"' recorded more than once (at different ranks) in "+file+".");
			}
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+file);}
		if(declaredK!=CandidateCChallengeDbBuilder.FROZEN_K){
			throw new RuntimeException(file+" declares K="+declaredK+", FROZEN at "
				+CandidateCChallengeDbBuilder.FROZEN_K+".");
		}
		if(byRank.size()!=declaredK){
			throw new RuntimeException("Focal family '"+focalFamilyId+"' has "+byRank.size()
				+" neighbor rows in "+file+", expected "+declaredK+".");
		}
		for(int r=1; r<=declaredK; r++){
			if(!byRank.containsKey(r)){
				throw new RuntimeException("Focal family '"+focalFamilyId+"' is missing neighbor_rank "
					+r+" in "+file+".");
			}
		}
		return new ArrayList<String>(byRank.values());
	}

	/** Every held-out TRUE member of ONE focal family's ONE fold: reads fold_manifest.tsv's rows for
	 * (focalFamilyId,foldId) to get the held-out group_ids (+declared sizes), then expands each via
	 * the family's own cluster.tsv (rep=group_id, member=member_id, no header) -- independently
	 * reimplemented from (not calling) {@link NegativeQueryManifestGenerator#expandHeldOutGroups},
	 * same declared-size-vs-actual-member-count validation, scoped to this one fold rather than
	 * every fold at once. */
	public static ArrayList<String> loadHeldOutMemberIds(final String clusterTsvPath,
			final String foldManifestPath, final String focalFamilyId, final int foldId){
		final HashMap<String,Long> targetGroupSizes=new HashMap<String,Long>();
		final ByteFile bfFold=ByteFile.makeByteFile(foldManifestPath, true);
		final LineParser1 lpFold=new LineParser1((byte)'\t');
		for(byte[] line=bfFold.nextLine(); line!=null; line=bfFold.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			lpFold.set(line);
			if(lpFold.terms()<4){throw new RuntimeException("Malformed fold_manifest row: "+new String(line));}
			final String familyId=lpFold.parseString(0);
			if(!familyId.equals(focalFamilyId)){continue;}
			final String groupId=lpFold.parseString(1);
			final long groupSize=lpFold.parseLong(2);
			final int rowFoldId=lpFold.parseInt(3);
			if(rowFoldId!=foldId){continue;}
			if(targetGroupSizes.put(groupId, groupSize)!=null){
				throw new RuntimeException("Duplicate group_id '"+groupId+"' for family '"+focalFamilyId
					+"' fold "+foldId+" in "+foldManifestPath+".");
			}
		}
		if(bfFold.close()){throw new RuntimeException("I/O error reading "+foldManifestPath);}
		if(targetGroupSizes.isEmpty()){
			throw new RuntimeException("No held-out groups for family '"+focalFamilyId+"' fold "+foldId
				+" in "+foldManifestPath+".");
		}

		final HashMap<String,ArrayList<String>> membersByGroup=new HashMap<String,ArrayList<String>>();
		final HashSet<String> seenMemberIds=new HashSet<String>();
		final ByteFile bfCluster=ByteFile.makeByteFile(clusterTsvPath, true);
		final LineParser1 lpCluster=new LineParser1((byte)'\t');
		for(byte[] line=bfCluster.nextLine(); line!=null; line=bfCluster.nextLine()){
			if(line.length==0){continue;}
			lpCluster.set(line);
			if(lpCluster.terms()<2){throw new RuntimeException("Malformed cluster.tsv row: "+new String(line));}
			final String groupId=lpCluster.parseString(0), memberId=lpCluster.parseString(1);
			if(!targetGroupSizes.containsKey(groupId)){continue;}
			if(!seenMemberIds.add(memberId)){
				throw new RuntimeException("cluster.tsv "+clusterTsvPath+" has member '"+memberId
					+"' appearing more than once among this fold's held-out groups.");
			}
			membersByGroup.computeIfAbsent(groupId, k -> new ArrayList<String>()).add(memberId);
		}
		if(bfCluster.close()){throw new RuntimeException("I/O error reading "+clusterTsvPath);}

		final ArrayList<String> result=new ArrayList<String>();
		for(final HashMap.Entry<String,Long> e : targetGroupSizes.entrySet()){
			final ArrayList<String> members=membersByGroup.get(e.getKey());
			final long actual=(members==null) ? 0 : members.size();
			if(actual!=e.getValue()){
				throw new RuntimeException("Group '"+e.getKey()+"': fold manifest declares size "
					+e.getValue()+" but cluster.tsv "+clusterTsvPath+" has "+actual+" member row(s) for "
					+"it -- fold manifest and cluster.tsv have diverged.");
			}
			result.addAll(members);
		}
		java.util.Collections.sort(result);
		return result;
	}

	/** Every distinct fold_id appearing in fold_manifest.tsv for one family -- both runAssay's fold
	 * iteration source and aggregateFamily's expected-coverage reference (Elly's addition, 2026-09-01:
	 * "prove exactly that set is scored once"). */
	static HashSet<Integer> loadExpectedFoldIds(final String foldManifestPath, final String familyId){
		final HashSet<Integer> foldIds=new HashSet<Integer>();
		final ByteFile bf=ByteFile.makeByteFile(foldManifestPath, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			lp.set(line);
			if(lp.terms()<4){throw new RuntimeException("Malformed fold_manifest row: "+new String(line));}
			if(!lp.parseString(0).equals(familyId)){continue;}
			foldIds.add(lp.parseInt(3));
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+foldManifestPath);}
		if(foldIds.isEmpty()){
			throw new RuntimeException("No fold_manifest.tsv rows for family '"+familyId+"' in "+foldManifestPath+".");
		}
		return foldIds;
	}

	/** ONE pass over a member-sequence FASTA, extracting RAW (unencoded) sequence bytes for needed
	 * ids only -- mirrors {@code NegativeQueryManifestGenerator.buildLengthIndex}'s needed-IDs-only
	 * streaming pattern, independently reimplemented here for full sequence bytes instead of just
	 * lengths. */
	public static HashMap<String,byte[]> extractSequences(final String fastaPath, final HashSet<String> neededIds){
		final HashMap<String,byte[]> seqs=new HashMap<String,byte[]>();
		final ByteFile bf=ByteFile.makeByteFile(fastaPath, true);
		String curId=null;
		ByteBuilder curSeq=null;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='>'){
				finalizeSeqRecord(curId, curSeq, neededIds, seqs, fastaPath);
				int end=line.length;
				for(int i=1; i<line.length; i++){if(line[i]==' ' || line[i]=='\t'){end=i; break;}}
				curId=new String(line, 1, end-1);
				curSeq=(neededIds.contains(curId)) ? new ByteBuilder() : null;
			}else if(curSeq!=null){
				curSeq.append(line);
			}
		}
		finalizeSeqRecord(curId, curSeq, neededIds, seqs, fastaPath);
		if(bf.close()){throw new RuntimeException("I/O error reading "+fastaPath);}
		//CP-03: fail closed at the corpus-binding boundary -- a partial map must never propagate a
		//null-sequence Query downstream. Name every missing id so the diagnostic is actionable.
		if(seqs.size()!=neededIds.size()){
			final ArrayList<String> missing=new ArrayList<String>();
			for(final String id : neededIds){
				if(!seqs.containsKey(id)){missing.add(id);}
			}
			java.util.Collections.sort(missing);
			throw new RuntimeException("extractSequences: "+missing.size()+" needed id(s) not found in "
				+fastaPath+": "+missing);
		}
		return seqs;
	}

	static void finalizeSeqRecord(final String curId, final ByteBuilder curSeq,
			final HashSet<String> neededIds, final HashMap<String,byte[]> seqs, final String fastaPath){
		if(curId==null || curSeq==null || !neededIds.contains(curId)){return;}
		if(seqs.containsKey(curId)){
			throw new RuntimeException("Member id '"+curId+"' appears more than once as a needed record "
				+"in "+fastaPath+" -- a duplicate here is a real corpus-integrity violation.");
		}
		seqs.put(curId, curSeq.toBytes());
	}

	/*--------------------------------------------------------------*/
	/*----------------      Candidate A scorer               --------*/
	/*--------------------------------------------------------------*/

	static final class ArtifactManifestRow { String artifactId, method, format, path; }

	/** Candidate A ({@code fasta_consensus}) scorer: binds each needed artifact_id to its
	 * {@code method}/{@code format}/{@code path} via {@code artifact_manifest.tsv} (design doc #2 --
	 * never an implicit/injected path map). Preloads only the manifest's id->row index at
	 * construction; lazily loads+BLOSUM62-encodes ({@link ConsensusRepBuilder#encodeLenient}, the
	 * project's existing encoder, not reinvented)+caches each artifact's sequence bytes the first
	 * time it is actually needed, then loops in memory calling the already-sealed {@link
	 * FamilyArtifactScore#flatScore} per (query,artifact) pair -- deterministic under input-order/
	 * thread perturbation since encoding and scoring are both pure functions of the artifact's own
	 * on-disk bytes. */
	public static final class FlatConsensusScorer implements ArtifactScorer {
		private final HashMap<String,ArtifactManifestRow> rowsById;
		private final HashMap<String,byte[]> encodedCache=new HashMap<String,byte[]>();

		public FlatConsensusScorer(final String artifactManifestFile){
			this.rowsById=loadArtifactManifest(artifactManifestFile);
		}

		private byte[] encodedConsensusFor(final String artifactId){
			byte[] enc=encodedCache.get(artifactId);
			if(enc!=null){return enc;}
			final ArtifactManifestRow row=rowsById.get(artifactId);
			if(row==null){
				throw new RuntimeException("FlatConsensusScorer: artifact_id '"+artifactId
					+"' has no row in the artifact manifest.");
			}
			if(!row.method.equals("A")){
				throw new RuntimeException("FlatConsensusScorer: artifact_id '"+artifactId
					+"' has method='"+row.method+"', expected 'A'.");
			}
			if(!row.format.equals("fasta_consensus")){
				throw new RuntimeException("FlatConsensusScorer: artifact_id '"+artifactId
					+"' has format='"+row.format+"', expected 'fasta_consensus'.");
			}
			final byte[] raw=readSingleSequenceFasta(row.path);
			enc=ConsensusRepBuilder.encodeLenient(raw);
			encodedCache.put(artifactId, enc);
			return enc;
		}

		/** Every path in the artifact manifest this scorer was constructed from -- stable and complete
		 * at construction time, sorted for deterministic provenance output (Elly's design requirement,
		 * 2026-09-02: provenance must tie to the actual bytes a scorer consumes, not merely to a
		 * manifest entry that might diverge from what's really read). */
		@Override
		public ArrayList<String> consumedFilePaths(){
			final ArrayList<String> paths=new ArrayList<String>();
			for(final ArtifactManifestRow r : rowsById.values()){paths.add(r.path);}
			java.util.Collections.sort(paths);
			return paths;
		}

		@Override
		public double[][] scoreChallenge(final ArrayList<Query> queries, final ArrayList<String> orderedArtifactIds,
				final String focalFamilyId, final int foldId){
			final byte[][] encodedArtifacts=new byte[orderedArtifactIds.size()][];
			for(int a=0; a<orderedArtifactIds.size(); a++){
				encodedArtifacts[a]=encodedConsensusFor(orderedArtifactIds.get(a));
			}
			final double[][] matrix=new double[queries.size()][orderedArtifactIds.size()];
			for(int q=0; q<queries.size(); q++){
				final byte[] encodedQuery=ConsensusRepBuilder.encodeLenient(queries.get(q).sequence);
				for(int a=0; a<orderedArtifactIds.size(); a++){
					matrix[q][a]=FamilyArtifactScore.flatScore(encodedQuery, encodedArtifacts[a]);
				}
			}
			return matrix;
		}
	}

	/** artifact_manifest.tsv (design doc #2): {@code #artifact_id\tmethod\tformat\tpath}. Validates
	 * the exact 4-column header, rejects a duplicate header, requires exactly 4 fields per row (never
	 * NF&lt;4 or NF&gt;4), rejects an empty required field, and validates each path's safety BEFORE any
	 * artifact is ever read (CP-05/CP-06/CP-07). */
	static HashMap<String,ArtifactManifestRow> loadArtifactManifest(final String file){
		final HashMap<String,ArtifactManifestRow> byId=new HashMap<String,ArtifactManifestRow>();
		final ByteFile bf=ByteFile.makeByteFile(file, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		boolean sawHeader=false;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='#'){
				final String h=new String(line, 1, line.length-1);
				final String[] names=h.split("\t", -1);
				if(names.length>0 && names[0].trim().equalsIgnoreCase("artifact_id")){
					if(sawHeader){throw new RuntimeException("Duplicate artifact_id header in "+file);}
					if(names.length!=4 || !names[1].trim().equalsIgnoreCase("method")
							|| !names[2].trim().equalsIgnoreCase("format") || !names[3].trim().equalsIgnoreCase("path")){
						throw new RuntimeException("Malformed artifact_manifest header in "+file
							+": expected exactly [artifact_id, method, format, path], got "
							+java.util.Arrays.toString(names));
					}
					sawHeader=true;
				}
				continue;
			}
			if(!sawHeader){throw new RuntimeException("Data row before a valid header in "+file);}
			lp.set(line);
			if(lp.terms()!=4){
				throw new RuntimeException("Malformed artifact_manifest row in "+file+" (expected exactly 4 "
					+"fields, got "+lp.terms()+"): "+new String(line));
			}
			final ArtifactManifestRow r=new ArtifactManifestRow();
			r.artifactId=lp.parseString(0);
			r.method=lp.parseString(1);
			r.format=lp.parseString(2);
			r.path=lp.parseString(3);
			if(r.artifactId.isEmpty() || r.method.isEmpty() || r.format.isEmpty() || r.path.isEmpty()){
				throw new RuntimeException("artifact_manifest row in "+file+" has an empty required field: "
					+new String(line));
			}
			validateSafePath(r.path, "artifact_manifest "+file+" artifact '"+r.artifactId+"'");
			if(byId.put(r.artifactId, r)!=null){
				throw new RuntimeException("Duplicate artifact_id '"+r.artifactId+"' in "+file);
			}
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+file);}
		if(!sawHeader){throw new RuntimeException("No artifact_id header found in "+file);}
		return byId;
	}

	/** Path-safety validation shared by artifact-manifest rows (CP-05/06/07) and the selected sample
	 * row's own cluster_path (Elly's addition, 2026-09-01): non-empty, absolute, not a symlink (never
	 * follow one), exists, and is a regular file (not a directory) -- checked BEFORE any read/hash. */
	static void validateSafePath(final String path, final String contextLabel){
		if(path.isEmpty()){
			throw new RuntimeException(contextLabel+" has an empty path.");
		}
		final Path p=Paths.get(path);
		if(!p.isAbsolute()){
			throw new RuntimeException(contextLabel+" path '"+path+"' is not absolute.");
		}
		if(Files.isSymbolicLink(p)){
			throw new RuntimeException(contextLabel+" path '"+path+"' is a symlink -- refusing to follow.");
		}
		if(!Files.exists(p)){
			throw new RuntimeException(contextLabel+" path '"+path+"' does not exist.");
		}
		//Sayu's catch, 2026-09-02: rejecting only directories still let a FIFO/device/socket through --
		//isRegularFile strictly subsumes "not a directory" and additionally rejects every special-file
		//type, so a blocking-read FIFO or a device node can never reach ByteFile.
		if(!Files.isRegularFile(p)){
			throw new RuntimeException(contextLabel+" path '"+path+"' is not a regular file (directory, "
				+"FIFO, device, or socket) -- refusing to read it.");
		}
	}

	/** Reads a single-sequence FASTA's raw ASCII residues (all sequence lines after the one header
	 * concatenated), requiring EXACTLY one record. */
	static byte[] readSingleSequenceFasta(final String path){
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		int records=0;
		ByteBuilder seq=new ByteBuilder();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='>'){records++;}
			else{seq.append(line);}
		}
		if(bf.close()){throw new RuntimeException("I/O error reading "+path);}
		if(records!=1){
			throw new RuntimeException(path+" has "+records+" FASTA record(s), expected exactly 1 "
				+"(single-sequence consensus per design doc #2).");
		}
		if(seq.length()==0){
			throw new RuntimeException(path+" has a header but zero residue bytes -- a consensus "
				+"sequence must not be empty.");
		}
		return seq.toBytes();
	}

	/*--------------------------------------------------------------*/
	/*----------------      Provenance / hashing              -------*/
	/*--------------------------------------------------------------*/

	/** Streaming SHA-256 of a file's bytes -- never loads the whole file into memory (some consumed
	 * files, e.g. the member-sequences FASTA, are corpus-scale). */
	static String sha256File(final String path){
		try{
			final MessageDigest md=MessageDigest.getInstance("SHA-256");
			try(java.io.InputStream is=Files.newInputStream(Paths.get(path))){
				final byte[] buf=new byte[1<<20];
				int n;
				while((n=is.read(buf))>0){md.update(buf, 0, n);}
			}
			final byte[] digest=md.digest();
			final StringBuilder sb=new StringBuilder(digest.length*2);
			for(final byte b : digest){
				final String hex=Integer.toHexString(b&0xff);
				if(hex.length()<2){sb.append('0');}
				sb.append(hex);
			}
			return sb.toString();
		}catch(final Exception e){
			throw new RuntimeException("Failed to SHA-256 "+path, e);
		}
	}

	/** CP-11: re-hashes every path in {@code paths} and rejects if any no longer matches the hash
	 * captured in {@code startHashes} (keyed the same way, e.g. role name -> path vs role name ->
	 * hash) -- the drift check that must run AFTER all reads/binding and BEFORE any output is
	 * published. */
	static void assertNoHashDrift(final Map<String,String> paths, final Map<String,String> startHashes,
			final String label){
		for(final Map.Entry<String,String> e : paths.entrySet()){
			final String nowHash=sha256File(e.getValue());
			final String was=startHashes.get(e.getKey());
			if(!nowHash.equals(was)){
				throw new RuntimeException(label+" '"+e.getKey()+"' ("+e.getValue()+") changed: hash was "
					+was+" at capture time, is "+nowHash+" now -- refusing to publish output against "
					+"drifted input.");
			}
		}
	}

	/** CP-11 variant for maps keyed BY PATH itself (e.g. {@code clusterHashesStart}, path->hash) --
	 * re-hashes each key (the path) and rejects if it no longer matches its captured value (the
	 * hash). */
	static void assertNoPathHashDrift(final Map<String,String> pathToStartHash, final String label){
		for(final Map.Entry<String,String> e : pathToStartHash.entrySet()){
			final String nowHash=sha256File(e.getKey());
			if(!nowHash.equals(e.getValue())){
				throw new RuntimeException(label+" '"+e.getKey()+"' changed: hash was "+e.getValue()
					+" at capture time, is "+nowHash+" now -- refusing to publish output against "
					+"drifted input.");
			}
		}
	}

	/** Everything the §7 output writers need to emit real provenance headers (design doc §7): every
	 * consumed upstream file's path+SHA-256 (fixed inputs keyed by role name; cluster.tsv files keyed
	 * by their own path since there is one per family), a caller-supplied scorer description, and the
	 * expected-vs-observed query x challenge-set cell counts computed during the run. */
	static final class ProvenanceHeader {
		final LinkedHashMap<String,String> fixedInputPaths, fixedInputHashes, clusterHashes, scorerFileHashes;
		final String scorerDescription;
		final long expectedCells, observedCells;
		ProvenanceHeader(final LinkedHashMap<String,String> fixedInputPaths, final LinkedHashMap<String,String> fixedInputHashes,
				final LinkedHashMap<String,String> clusterHashes, final LinkedHashMap<String,String> scorerFileHashes,
				final String scorerDescription, final long expectedCells, final long observedCells){
			this.fixedInputPaths=fixedInputPaths; this.fixedInputHashes=fixedInputHashes; this.clusterHashes=clusterHashes;
			this.scorerFileHashes=scorerFileHashes;
			this.scorerDescription=scorerDescription; this.expectedCells=expectedCells; this.observedCells=observedCells;
		}
	}

	/** Appends the full §7 provenance block (CP-10): one line per fixed upstream input, one line per
	 * distinct cluster.tsv consumed, the scorer description, and expected/observed cell counts --
	 * asserting the two counts agree (they always should: extractSequences now fails closed on any
	 * missing id, and scoreFocalFold asserts its own returned-matrix row count == query count, so a
	 * mismatch here means a real, otherwise-undetected wiring bug in the caller). */
	static void appendProvenanceHeader(final ByteBuilder bb, final ProvenanceHeader prov){
		for(final Map.Entry<String,String> e : prov.fixedInputPaths.entrySet()){
			bb.append("#input_").append(e.getKey()).append('\t').append(e.getValue()).append('\t')
				.append(prov.fixedInputHashes.get(e.getKey())).append('\n');
		}
		for(final Map.Entry<String,String> e : prov.clusterHashes.entrySet()){
			bb.append("#input_cluster_tsv\t").append(e.getKey()).append('\t').append(e.getValue()).append('\n');
		}
		for(final Map.Entry<String,String> e : prov.scorerFileHashes.entrySet()){
			bb.append("#input_scorer_file\t").append(e.getKey()).append('\t').append(e.getValue()).append('\n');
		}
		if(prov.scorerDescription==null || prov.scorerDescription.trim().isEmpty()){
			throw new RuntimeException("Provenance requires a non-empty scorer description.");
		}
		bb.append("#scorer\t").append(prov.scorerDescription).append('\n');
		bb.append("#expected_q_x_challenge_cells\t").append(prov.expectedCells).append('\n');
		bb.append("#observed_q_x_challenge_cells\t").append(prov.observedCells).append('\n');
		if(prov.expectedCells!=prov.observedCells){
			throw new RuntimeException("Provenance mismatch: expected "+prov.expectedCells
				+" query x challenge-set score cells but observed "+prov.observedCells+" -- refusing to "
				+"publish output whose own accounting doesn't close.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output writers             --------*/
	/*--------------------------------------------------------------*/

	/** §7 assay_top1_confusion.tsv: one row per family, with a real provenance header (CP-10). */
	public static void writeTop1ConfusionTsv(final String outFile, final ArrayList<FamilySummary> summaries,
			final String method, final ProvenanceHeader prov){
		final ByteBuilder bb=new ByteBuilder();
		bb.append("#schema_version\t1\n#method\t").append(method).append('\n');
		appendProvenanceHeader(bb, prov);
		bb.append("#cell_id\tfamily_id\tmethod\ttotal_held_out\ttop1_correct\ttop1_confused\ttop1_tied\t"
			+"self_consistency\tconfusion_mass\n");
		for(final FamilySummary s : summaries){
			bb.append(s.cellId).append('\t').append(s.familyId).append('\t').append(method).append('\t')
				.append(s.totalHeldOut).append('\t').append(s.correct).append('\t').append(s.confused)
				.append('\t').append(s.tied).append('\t').append(formatDecimal(s.selfConsistency))
				.append('\t').append(formatDecimal(s.confusionMass)).append('\n');
		}
		writeFile(outFile, bb.toString());
	}

	/** §7 assay_margin_3a.tsv: one row per positive query, plus family-first-then-macro pass-fraction
	 * summary lines (family lines, then per-cell macro = unweighted mean of that cell's family
	 * lines, then a secondary-only per-cell pooled line), with a real provenance header (CP-10). The
	 * pooled line counts {@code passedMargin} DIRECTLY from the raw per-cell-filtered results -- never
	 * {@code Math.round} of an already-rounded fraction (CP-12). */
	public static void writeMargin3aTsv(final String outFile, final ArrayList<QueryResult> allResults,
			final ArrayList<FamilySummary> summaries, final String method, final double marginMin,
			final ProvenanceHeader prov){
		final ByteBuilder bb=new ByteBuilder();
		bb.append("#schema_version\t1\n#method\t").append(method).append("\n#margin_min\t")
			.append(formatDecimal(marginMin)).append('\n');
		appendProvenanceHeader(bb, prov);
		bb.append("#cell_id\tfamily_id\tfold_id\tmember_id\tmethod\ts_own\ts_max_other\tmargin\tpassed_at_Mx\n");
		for(final QueryResult r : allResults){
			bb.append(r.cellId).append('\t').append(r.focalFamilyId).append('\t').append(r.foldId).append('\t')
				.append(r.memberId).append('\t').append(method).append('\t').append(formatDecimal(r.sOwn))
				.append('\t').append(formatDecimal(r.sMaxOther)).append('\t').append(formatDecimal(r.marginNative))
				.append('\t').append(r.passedMargin ? 1 : 0).append('\n');
		}
		for(final FamilySummary s : summaries){
			bb.append("#family_pass_fraction\t").append(s.cellId).append('\t').append(s.familyId).append('\t')
				.append(method).append('\t').append(formatDecimal(s.passFraction)).append('\t').append(s.totalHeldOut)
				.append('\n');
		}
		final HashMap<String,ArrayList<FamilySummary>> byCell=new HashMap<String,ArrayList<FamilySummary>>();
		for(final FamilySummary s : summaries){
			byCell.computeIfAbsent(s.cellId, k -> new ArrayList<FamilySummary>()).add(s);
		}
		final ArrayList<String> cellIdsSorted=new ArrayList<String>(byCell.keySet());
		java.util.Collections.sort(cellIdsSorted);
		for(final String cellId : cellIdsSorted){
			final ArrayList<FamilySummary> cellFamilies=byCell.get(cellId);
			double macroSum=0;
			for(final FamilySummary s : cellFamilies){macroSum+=s.passFraction;}
			final double macro=macroSum/cellFamilies.size();
			bb.append("#cell_pass_fraction_macro\t").append(cellId).append('\t').append(method).append('\t')
				.append(formatDecimal(macro)).append('\t').append(cellFamilies.size()).append('\n');
			//CP-12 fix: count directly from the raw per-cell-filtered results, never Math.round a fraction.
			int pooledPassed=0, pooledTotal=0;
			for(final QueryResult r : allResults){
				if(r.cellId.equals(cellId)){
					pooledTotal++;
					if(r.passedMargin){pooledPassed++;}
				}
			}
			bb.append("#cell_pass_fraction_pooled\t").append(cellId).append('\t').append(method).append('\t')
				.append(formatDecimal(pooledPassed/(double)pooledTotal)).append('\t').append(pooledTotal).append('\n');
		}
		writeFile(outFile, bb.toString());
	}

	static String formatDecimal(final double x){
		return String.format(java.util.Locale.ROOT, "%.12f", x);
	}

	static void writeFile(final String path, final String content){
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

	public static void main(final String[] args) throws Exception{
		if(args.length==1 && args[0].equalsIgnoreCase("selftest")){selftest(); return;}
		if(args.length>=1 && args[0].equalsIgnoreCase("runassay")){
			runAssayFromArgs(args);
			return;
		}
		throw new RuntimeException("Usage: java -ea prot.Step1AssayDriver selftest\n"
			+"   or: java -ea prot.Step1AssayDriver runassay samplemanifest=... foldmanifest=... "
			+"rolemanifest=... neighbors=... artifactmanifest=... memberseqs=... method=A marginmin=5.0 "
			+"scorerdesc='...' top1out=... marginout=...");
	}

	/** Item 1a's CLI entry point (Elly's/Sayu's finding 6): parses {@code flag=value} args and
	 * dispatches to {@link #runAssay}. Currently only {@code method=A} ({@link FlatConsensusScorer})
	 * is wired -- a future Candidate-C {@link ArtifactScorer} needs its own in-process caller since it
	 * must shell {@code hmmscan} per the design doc's own dispatch contract, not a generic CLI flag. */
	static void runAssayFromArgs(final String[] args){
		String sampleManifest=null, foldManifest=null, roleManifest=null, neighbors=null,
			artifactManifest=null, memberSeqs=null, method=null, scorerDesc=null, top1Out=null, marginOut=null;
		double marginMin=Double.NaN;
		for(int i=1; i<args.length; i++){
			final String arg=args[i];
			final int eq=arg.indexOf('=');
			if(eq<0){throw new RuntimeException("Bad argument (expected flag=value): "+arg);}
			final String flag=arg.substring(0, eq).toLowerCase();
			final String value=arg.substring(eq+1);
			switch(flag){
				case "samplemanifest": sampleManifest=value; break;
				case "foldmanifest": foldManifest=value; break;
				case "rolemanifest": roleManifest=value; break;
				case "neighbors": neighbors=value; break;
				case "artifactmanifest": artifactManifest=value; break;
				case "memberseqs": memberSeqs=value; break;
				case "method": method=value; break;
				case "marginmin": marginMin=Double.parseDouble(value); break;
				case "scorerdesc": scorerDesc=value; break;
				case "top1out": top1Out=value; break;
				case "marginout": marginOut=value; break;
				default: throw new RuntimeException("Unknown flag: "+flag);
			}
		}
		if(sampleManifest==null || foldManifest==null || roleManifest==null || neighbors==null
				|| artifactManifest==null || memberSeqs==null || method==null || Double.isNaN(marginMin)
				|| scorerDesc==null || top1Out==null || marginOut==null){
			throw new RuntimeException("runAssay requires: samplemanifest=,foldmanifest=,rolemanifest=,"
				+"neighbors=,artifactmanifest=,memberseqs=,method=,marginmin=,scorerdesc=,top1out=,marginout=");
		}
		if(!method.equals("A")){
			throw new RuntimeException("CLI runassay currently only supports method=A (FlatConsensusScorer); "
				+"a future Candidate-C ArtifactScorer must be wired in-process, not via this generic CLI.");
		}
		final ArtifactScorer scorer=new FlatConsensusScorer(artifactManifest);
		runAssay(sampleManifest, foldManifest, roleManifest, neighbors, artifactManifest, memberSeqs,
			scorer, method, marginMin, scorerDesc, top1Out, marginOut);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Self-test            ----------------*/
	/*--------------------------------------------------------------*/
	// Synthetic fixtures only -- neither candidate has real cluster artifacts yet. Reuses
	// CandidateCChallengeDbBuilder's own K=20 fixture helpers (frozenKNeighbors/roleRow/
	// decoyRowsForFamilies/neighborBlock, already widened to package-visible this session) rather
	// than re-deriving 20-neighbor boilerplate.

	/** Deterministic canned scorer for the end-to-end test: returns a caller-supplied fixed score
	 * row per query id, ignoring focalFamilyId/foldId/orderedArtifactIds (both folds of this
	 * fixture share the identical 21-artifact order, focal-then-neighbor, so one row per query id
	 * is unambiguous). */
	static final class CannedScorer implements ArtifactScorer {
		final HashMap<String,double[]> rowsByQueryId;
		CannedScorer(final HashMap<String,double[]> rowsByQueryId){this.rowsByQueryId=rowsByQueryId;}
		@Override
		public double[][] scoreChallenge(final ArrayList<Query> queries, final ArrayList<String> orderedArtifactIds,
				final String focalFamilyId, final int foldId){
			final double[][] matrix=new double[queries.size()][];
			for(int i=0; i<queries.size(); i++){
				final double[] row=rowsByQueryId.get(queries.get(i).memberId);
				if(row==null){throw new RuntimeException("CannedScorer: no canned row for query '"
					+queries.get(i).memberId+"'.");}
				matrix[i]=row;
			}
			return matrix;
		}
	}

	static void selftest() throws Exception{
		testLoaders();
		System.err.println("  selftest[loaders: findSelectedFamilyRow/loadCellId exactly-one-row "
			+"requirement, loadChallengeSetIds K=20 focal-then-neighbor order, loadHeldOutMemberIds "
			+"fold-scoped group expansion, extractSequences needed-IDs-only]: PASS");

		testFlatConsensusScorer();
		System.err.println("  selftest[FlatConsensusScorer: artifact_manifest method/format binding, "
			+"lazy encode+cache, scores cross-check against a direct FamilyArtifactScore.flatScore call]: PASS");

		testEndToEndTwoFoldFamily();
		System.err.println("  selftest[END-TO-END: 2 folds of one family, all 3 Top1Buckets exercised, "
			+"family-first pooling ACROSS folds, exhaustiveness invariant, output TSV row/summary "
			+"correctness]: PASS");

		testScoreFocalFoldDimensionAndFiniteGuards();
		System.err.println("  selftest[scoreFocalFold rejects a wrong-shaped or non-finite score "
			+"matrix before any metric call]: PASS");

		testRunAssayEndToEnd();
		System.err.println("  selftest[runAssay END-TO-END: full orchestration over sample/fold/role/"
			+"neighbors/artifact manifests, independently-recomputed SHA-256 provenance, expected==observed "
			+"Q x21 cell counts, direct-count pooled pass-fraction]: PASS");

		testCP01ClusterPathBoundToSelectedRow();
		testCP02DuplicateNeighborFamilyRejected();
		testCP03MissingSequenceIdRejected();
		testCP04AggregateFamilyHardening();
		testCP05to07ArtifactManifestHardening();
		testCP08And09ChallengeSetBoundary();
		testCP10ProvenanceMismatchRejected();
		testCP11HashDriftRejected();
		testCP12PooledPassFractionCountedDirectly();
		testFoldProcessingOrderDeterministic();
		System.err.println("  selftest[fold processing order is canonically ascending regardless of "
			+"fold_manifest.tsv's physical row order (Elly's finding, 2026-09-02)]: PASS");

		testArtifactFileProvenanceViaFlatConsensusScorer();
		testEmptyScorerDescriptionRejected();
		testZeroResidueFastaRejected();
		testNonRegularFileRejected();
		System.err.println("  selftest[Sayu's 2026-09-02 blockers: real FlatConsensusScorer artifact "
			+"files are individually SHA-256'd and drift-checked via the scorer's own consumedFilePaths() "
			+"contract, never merely trusted from a manifest entry; empty scorer description rejected; "
			+"zero-residue single-record FASTA rejected; validateSafePath rejects non-regular files "
			+"(FIFO/device/socket)]: PASS");
		System.err.println("  selftest[CP-01..CP-12 (records/STEP1ASSAYDRIVER_GATE_CP_MATRIX_v1.md): "
			+"cluster_path binding closure, duplicate neighbor family, missing sequence id, mixed-cell/"
			+"duplicate-member/fold-coverage aggregation, artifact-manifest header/row/path hardening, "
			+"21-item/focal-index-0 boundary, provenance count mismatch, hash drift, direct-count pooled "
			+"pass-fraction]: PASS");

		System.err.println("Step1AssayDriver SELFTEST PASS.");
	}

	static void testLoaders() throws Exception{
		final String dir=System.getProperty("java.io.tmpdir")+"/step1assaydriver_selftest_"+System.nanoTime();
		final String work=dir+"/loaders";
		new java.io.File(work).mkdirs();

		//-- findSelectedFamilyRow / loadCellId --
		final String sampleFile=work+"/sample.tsv";
		writeFile(sampleFile, "#rep_id\tselected\tcell_id\tcluster_path\n"
			+"fZ\t1\tQ1_F1\t"+work+"/fZ.cluster.tsv\n"
			+"fOther\t0\tQ2_F2\t"+work+"/fOther.cluster.tsv\n");
		final String cellId=loadCellId(sampleFile, "fZ");
		if(!cellId.equals("Q1_F1")){throw new RuntimeException("SELFTEST FAILED: loadCellId returned '"
			+cellId+"', expected 'Q1_F1'.");}
		expectThrows(()->loadCellId(sampleFile, "fOther"), "unselected family");
		expectThrows(()->loadCellId(sampleFile, "fNonexistent"), "family absent entirely");
		writeFile(work+"/sample_dup.tsv", "#rep_id\tselected\tcell_id\tcluster_path\n"
			+"fZ\t1\tQ1_F1\t"+work+"/fZ.cluster.tsv\n"
			+"fZ\t1\tQ2_F2\t"+work+"/fZ2.cluster.tsv\n");
		expectThrows(()->loadCellId(work+"/sample_dup.tsv", "fZ"), "two selected=1 rows for the same family");

		//-- loadChallengeSetIds (K=20, focal-then-neighbor order) --
		final String focalFamilyId="fZ";
		final String[] neighbors=CandidateCChallengeDbBuilder.frozenKNeighbors(focalFamilyId);
		final String roleFile=work+"/role.tsv";
		writeFile(roleFile, "#artifact_id\tfamily_id\trole\tfold_id\n"
			+CandidateCChallengeDbBuilder.roleRow(focalFamilyId+"__focal_f0", focalFamilyId, "focal_fold_local", 0)
			+CandidateCChallengeDbBuilder.decoyRowsForFamilies(neighbors));
		final String neighborsFile=work+"/neighbors.tsv";
		writeFile(neighborsFile, "#k\t20\n"+CandidateCChallengeDbBuilder.neighborBlock(focalFamilyId, neighbors));
		final ArrayList<String> challengeSet=loadChallengeSetIds(roleFile, neighborsFile, focalFamilyId, 0);
		if(challengeSet.size()!=21){throw new RuntimeException("SELFTEST FAILED: expected 21 challenge-set "
			+"ids, got "+challengeSet.size());}
		if(!challengeSet.get(0).equals("fZ__focal_f0")){throw new RuntimeException("SELFTEST FAILED: rank 0 "
			+"is '"+challengeSet.get(0)+"', expected the focal artifact.");}
		for(int i=0; i<20; i++){
			if(!challengeSet.get(i+1).equals(neighbors[i]+"__decoy")){
				throw new RuntimeException("SELFTEST FAILED: rank "+(i+1)+" is '"+challengeSet.get(i+1)
					+"', expected '"+neighbors[i]+"__decoy'.");
			}
		}
		expectThrows(()->loadChallengeSetIds(roleFile, neighborsFile, focalFamilyId, 1), "no focal row for fold 1");

		//-- loadHeldOutMemberIds --
		final String clusterFile=work+"/fZ.cluster.tsv";
		writeFile(clusterFile, "m1\tm1\n"+"m1\tm2\n"+"m3\tm3\n"+"m4\tm4\n"+"m4\tm5\n");
		final String foldManifestFile=work+"/fold_manifest.tsv";
		writeFile(foldManifestFile, "#family_id\tgroup_id\tgroup_size\tfold_id\tgroup_hash\n"
			+"fZ\tm1\t2\t0\t"+new String(new char[16]).replace((char)0, 'a')+"\n"
			+"fZ\tm3\t1\t0\t"+new String(new char[16]).replace((char)0, 'b')+"\n"
			+"fZ\tm4\t2\t1\t"+new String(new char[16]).replace((char)0, 'c')+"\n");
		final ArrayList<String> heldOutFold0=loadHeldOutMemberIds(clusterFile, foldManifestFile, "fZ", 0);
		final ArrayList<String> expectedFold0=new ArrayList<String>(java.util.Arrays.asList("m1","m2","m3"));
		java.util.Collections.sort(expectedFold0);
		if(!heldOutFold0.equals(expectedFold0)){
			throw new RuntimeException("SELFTEST FAILED: fold 0 held-out members "+heldOutFold0
				+" != expected "+expectedFold0);
		}
		final ArrayList<String> heldOutFold1=loadHeldOutMemberIds(clusterFile, foldManifestFile, "fZ", 1);
		if(!heldOutFold1.equals(new ArrayList<String>(java.util.Arrays.asList("m4","m5")))){
			throw new RuntimeException("SELFTEST FAILED: fold 1 held-out members "+heldOutFold1
				+" != expected [m4, m5]");
		}
		//declared-size-vs-actual divergence must crash loud.
		writeFile(work+"/fold_manifest_bad.tsv", "#family_id\tgroup_id\tgroup_size\tfold_id\tgroup_hash\n"
			+"fZ\tm1\t99\t0\t"+new String(new char[16]).replace((char)0, 'a')+"\n");
		expectThrows(()->loadHeldOutMemberIds(clusterFile, work+"/fold_manifest_bad.tsv", "fZ", 0),
			"declared group size does not match cluster.tsv");

		//-- extractSequences (needed-IDs-only) --
		final String fastaFile=work+"/seqs.fasta";
		writeFile(fastaFile, ">m1 some header text\nMKVL\nACDE\n>m2\nMKVLACDE\n>m3\nQQQQ\n>unrelated\nZZZZ\n");
		final HashSet<String> needed=new HashSet<String>(java.util.Arrays.asList("m1","m3"));
		final HashMap<String,byte[]> extracted=extractSequences(fastaFile, needed);
		if(extracted.size()!=2 || !extracted.containsKey("m1") || !extracted.containsKey("m3")){
			throw new RuntimeException("SELFTEST FAILED: extractSequences returned "+extracted.keySet()
				+", expected exactly {m1, m3}.");
		}
		if(!new String(extracted.get("m1")).equals("MKVLACDE")){
			throw new RuntimeException("SELFTEST FAILED: m1 sequence '"+new String(extracted.get("m1"))
				+"' != expected 'MKVLACDE' (multi-line concatenation).");
		}
	}

	static void testFlatConsensusScorer() throws Exception{
		final String dir=System.getProperty("java.io.tmpdir")+"/step1assaydriver_selftest_"+System.nanoTime();
		final String work=dir+"/flatconsensus";
		new java.io.File(work).mkdirs();
		final String consensusA=work+"/a1.fasta";
		writeFile(consensusA, ">a1\nMKVLACDEFGHIKLMNPQRSTVWY\n");
		final String consensusB=work+"/a2.fasta";
		writeFile(consensusB, ">a2\nQQQQQQQQQQQQQQQQQQQQQQQQ\n");
		final String artifactManifest=work+"/artifact_manifest.tsv";
		writeFile(artifactManifest, "#artifact_id\tmethod\tformat\tpath\n"
			+"art1\tA\tfasta_consensus\t"+consensusA+"\n"
			+"art2\tA\tfasta_consensus\t"+consensusB+"\n");
		final FlatConsensusScorer scorer=new FlatConsensusScorer(artifactManifest);
		final ArrayList<Query> queries=new ArrayList<Query>();
		queries.add(new Query("q1", "MKVLACDEFGHIKLMNPQRSTVWY".getBytes()));
		final ArrayList<String> ordered=new ArrayList<String>(java.util.Arrays.asList("art1", "art2"));
		final double[][] matrix=scorer.scoreChallenge(queries, ordered, "fZ", 0);
		final double directA=FamilyArtifactScore.flatScore(
			ConsensusRepBuilder.encodeLenient("MKVLACDEFGHIKLMNPQRSTVWY".getBytes()),
			ConsensusRepBuilder.encodeLenient("MKVLACDEFGHIKLMNPQRSTVWY".getBytes()));
		final double directB=FamilyArtifactScore.flatScore(
			ConsensusRepBuilder.encodeLenient("MKVLACDEFGHIKLMNPQRSTVWY".getBytes()),
			ConsensusRepBuilder.encodeLenient("QQQQQQQQQQQQQQQQQQQQQQQQ".getBytes()));
		if(matrix[0][0]!=directA || matrix[0][1]!=directB){
			throw new RuntimeException("SELFTEST FAILED: FlatConsensusScorer scores ["+matrix[0][0]+","
				+matrix[0][1]+"] != direct FamilyArtifactScore.flatScore calls ["+directA+","+directB+"].");
		}
		if(!(matrix[0][0]>matrix[0][1])){
			throw new RuntimeException("SELFTEST FAILED: an identical-sequence artifact should score "
				+"higher than an unrelated one -- got own="+matrix[0][0]+" other="+matrix[0][1]);
		}
		//Wrong method/format must crash loud.
		final String badManifest=work+"/artifact_manifest_bad.tsv";
		writeFile(badManifest, "#artifact_id\tmethod\tformat\tpath\nart1\tC\thmm_profile\t"+consensusA+"\n");
		final FlatConsensusScorer badScorer=new FlatConsensusScorer(badManifest);
		expectThrows(()->badScorer.scoreChallenge(queries, new ArrayList<String>(java.util.Arrays.asList("art1")),
			"fZ", 0), "wrong method/format for a Candidate-A scorer");
	}

	@SuppressWarnings("unchecked")
	static void testEndToEndTwoFoldFamily() throws Exception{
		final String dir=System.getProperty("java.io.tmpdir")+"/step1assaydriver_selftest_"+System.nanoTime();
		final String work=dir+"/e2e";
		new java.io.File(work).mkdirs();

		final String focalFamilyId="fZ";
		final String[] neighbors=CandidateCChallengeDbBuilder.frozenKNeighbors(focalFamilyId);
		final String roleFile=work+"/role.tsv";
		writeFile(roleFile, "#artifact_id\tfamily_id\trole\tfold_id\n"
			+CandidateCChallengeDbBuilder.roleRow(focalFamilyId+"__focal_f0", focalFamilyId, "focal_fold_local", 0)
			+CandidateCChallengeDbBuilder.roleRow(focalFamilyId+"__focal_f1", focalFamilyId, "focal_fold_local", 1)
			+CandidateCChallengeDbBuilder.decoyRowsForFamilies(neighbors));
		final String neighborsFile=work+"/neighbors.tsv";
		writeFile(neighborsFile, "#k\t20\n"+CandidateCChallengeDbBuilder.neighborBlock(focalFamilyId, neighbors));
		final String clusterFile=work+"/fZ.cluster.tsv";
		writeFile(clusterFile, "m1\tm1\nm1\tm2\nm3\tm3\nm4\tm4\nm4\tm5\n");
		final String foldManifestFile=work+"/fold_manifest.tsv";
		writeFile(foldManifestFile, "#family_id\tgroup_id\tgroup_size\tfold_id\tgroup_hash\n"
			+"fZ\tm1\t2\t0\t"+new String(new char[16]).replace((char)0, 'a')+"\n"+"fZ\tm3\t1\t0\t"+new String(new char[16]).replace((char)0, 'b')+"\n"
			+"fZ\tm4\t2\t1\t"+new String(new char[16]).replace((char)0, 'c')+"\n");
		final String sampleFile=work+"/sample.tsv";
		writeFile(sampleFile, "#rep_id\tselected\tcell_id\tcluster_path\nfZ\t1\tQ1_F1\t"+clusterFile+"\n");
		final String fastaFile=work+"/seqs.fasta";
		writeFile(fastaFile, ">m1\nAAAA\n>m2\nAAAA\n>m3\nCCCC\n>m4\nGGGG\n>m5\nGGGG\n");

		final String cellId=loadCellId(sampleFile, focalFamilyId);
		final ArrayList<QueryResult> allResults=new ArrayList<QueryResult>();
		final double marginMin=5.0;

		//-- Fold 0: m1=CORRECT (unique own max), m2=CONFUSED (unique other max), m3=TIED (own/other tie) --
		final ArrayList<String> challenge0=loadChallengeSetIds(roleFile, neighborsFile, focalFamilyId, 0);
		final ArrayList<String> heldOut0=loadHeldOutMemberIds(clusterFile, foldManifestFile, focalFamilyId, 0);
		final HashMap<String,byte[]> seqs0=extractSequences(fastaFile, new HashSet<String>(heldOut0));
		final ArrayList<Query> queries0=new ArrayList<Query>();
		for(final String id : heldOut0){queries0.add(new Query(id, seqs0.get(id)));}
		final HashMap<String,double[]> canned0=new HashMap<String,double[]>();
		canned0.put("m1", rowWithOwnHigh(21, 10.0, 2.0));//own(0)=10, all others=2 -> unique own max -> CORRECT, margin=8>=5 passes
		canned0.put("m2", rowWithOtherHigh(21, 1.0, 9.0, 5));//own(0)=1, index5=9, rest=2 -> unique OTHER max -> CONFUSED
		canned0.put("m3", rowWithTie(21, 7.0, 3));//own(0)=7 AND index3=7 tied for max -> TIED
		allResults.addAll(scoreFocalFold(focalFamilyId, 0, cellId, challenge0, queries0,
			new CannedScorer(canned0), marginMin));

		//-- Fold 1: m4, m5 both CORRECT --
		final ArrayList<String> challenge1=loadChallengeSetIds(roleFile, neighborsFile, focalFamilyId, 1);
		final ArrayList<String> heldOut1=loadHeldOutMemberIds(clusterFile, foldManifestFile, focalFamilyId, 1);
		final HashMap<String,byte[]> seqs1=extractSequences(fastaFile, new HashSet<String>(heldOut1));
		final ArrayList<Query> queries1=new ArrayList<Query>();
		for(final String id : heldOut1){queries1.add(new Query(id, seqs1.get(id)));}
		final HashMap<String,double[]> canned1=new HashMap<String,double[]>();
		canned1.put("m4", rowWithOwnHigh(21, 20.0, 1.0));
		canned1.put("m5", rowWithOwnHigh(21, 20.0, 1.0));
		allResults.addAll(scoreFocalFold(focalFamilyId, 1, cellId, challenge1, queries1,
			new CannedScorer(canned1), marginMin));

		if(allResults.size()!=5){throw new RuntimeException("SELFTEST FAILED: expected 5 total held-out "
			+"query results across both folds, got "+allResults.size());}

		final HashSet<Integer> expectedFoldIds=new HashSet<Integer>(java.util.Arrays.asList(0, 1));
		final FamilySummary summary=aggregateFamily(focalFamilyId, cellId, allResults, expectedFoldIds);
		if(summary.totalHeldOut!=5 || summary.correct!=3 || summary.confused!=1 || summary.tied!=1){
			throw new RuntimeException("SELFTEST FAILED: family summary total="+summary.totalHeldOut
				+" correct="+summary.correct+" confused="+summary.confused+" tied="+summary.tied
				+", expected total=5 correct=3 confused=1 tied=1.");
		}
		if(Math.abs(summary.selfConsistency-0.6)>1e-9 || Math.abs(summary.confusionMass-0.2)>1e-9
				|| Math.abs(summary.tiedRate-0.2)>1e-9){
			throw new RuntimeException("SELFTEST FAILED: self_consistency="+summary.selfConsistency
				+" confusion_mass="+summary.confusionMass+" tied_rate="+summary.tiedRate
				+", expected 0.6/0.2/0.2.");
		}
		//passFraction: m1(margin8>=5 pass), m2(margin=1-9=-8 fail), m3(margin=0, tie -> maxOther=7,
		//own=7, margin=0<5 fail), m4/m5(margin=19>=5 pass) -> 3 of 5 pass = 0.6
		if(Math.abs(summary.passFraction-0.6)>1e-9){
			throw new RuntimeException("SELFTEST FAILED: passFraction="+summary.passFraction+", expected 0.6.");
		}

		final ArrayList<FamilySummary> summaries=new ArrayList<FamilySummary>(java.util.Arrays.asList(summary));
		final LinkedHashMap<String,String> fixedPaths=new LinkedHashMap<String,String>();
		fixedPaths.put("sample_manifest", sampleFile);
		fixedPaths.put("fold_manifest", foldManifestFile);
		fixedPaths.put("role_manifest", roleFile);
		fixedPaths.put("neighbors", neighborsFile);
		fixedPaths.put("member_seqs_fasta", fastaFile);
		final LinkedHashMap<String,String> fixedHashes=new LinkedHashMap<String,String>();
		for(final Map.Entry<String,String> e : fixedPaths.entrySet()){fixedHashes.put(e.getKey(), sha256File(e.getValue()));}
		final LinkedHashMap<String,String> clusterHashes=new LinkedHashMap<String,String>();
		clusterHashes.put(clusterFile, sha256File(clusterFile));
		final ProvenanceHeader prov=new ProvenanceHeader(fixedPaths, fixedHashes, clusterHashes,
			new LinkedHashMap<String,String>(), "CannedScorer test fixture", 105, 105);
		final String top1Out=work+"/assay_top1_confusion.tsv";
		writeTop1ConfusionTsv(top1Out, summaries, "A", prov);
		final String top1Text=readWholeFile(top1Out);
		if(!top1Text.contains("Q1_F1\tfZ\tA\t5\t3\t1\t1\t")){
			throw new RuntimeException("SELFTEST FAILED: assay_top1_confusion.tsv missing the expected "
				+"family row -- got:\n"+top1Text);
		}
		if(!top1Text.contains("#input_sample_manifest\t"+sampleFile+"\t"+fixedHashes.get("sample_manifest"))){
			throw new RuntimeException("SELFTEST FAILED: assay_top1_confusion.tsv missing/wrong provenance "
				+"header for sample_manifest -- got:\n"+top1Text);
		}
		final String marginOut=work+"/assay_margin_3a.tsv";
		writeMargin3aTsv(marginOut, allResults, summaries, "A", marginMin, prov);
		final String marginText=readWholeFile(marginOut);
		int dataRows=0;
		for(final String line : marginText.split("\n")){
			if(!line.isEmpty() && !line.startsWith("#")){dataRows++;}
		}
		if(dataRows!=5){throw new RuntimeException("SELFTEST FAILED: assay_margin_3a.tsv has "+dataRows
			+" data rows, expected 5.");}
		if(!marginText.contains("#family_pass_fraction\tQ1_F1\tfZ\tA\t"+formatDecimal(0.6))){
			throw new RuntimeException("SELFTEST FAILED: missing/wrong #family_pass_fraction line -- got:\n"
				+marginText);
		}
		if(!marginText.contains("#cell_pass_fraction_macro\tQ1_F1\tA\t"+formatDecimal(0.6))){
			throw new RuntimeException("SELFTEST FAILED: missing/wrong #cell_pass_fraction_macro line "
				+"(single-family cell, macro must equal that family's own pass fraction) -- got:\n"+marginText);
		}
	}

	static void testScoreFocalFoldDimensionAndFiniteGuards() throws Exception{
		final ArrayList<String> challengeSet=new ArrayList<String>();
		for(int i=0; i<21; i++){challengeSet.add("art"+i);}
		final ArrayList<Query> queries=new ArrayList<Query>();
		queries.add(new Query("q1", "AAAA".getBytes()));
		final HashMap<String,double[]> wrongWidth=new HashMap<String,double[]>();
		wrongWidth.put("q1", new double[20]);//one short of 21
		expectThrows(()->scoreFocalFold("fZ", 0, "Q1_F1", challengeSet, queries,
			new CannedScorer(wrongWidth), 5.0), "score row width != challenge-set size");

		final HashMap<String,double[]> nonFinite=new HashMap<String,double[]>();
		final double[] row=new double[21];
		row[3]=Double.NaN;
		nonFinite.put("q1", row);
		expectThrows(()->scoreFocalFold("fZ", 0, "Q1_F1", challengeSet, queries,
			new CannedScorer(nonFinite), 5.0), "non-finite score entry");
	}

	static void testRunAssayEndToEnd() throws Exception{
		final String dir=System.getProperty("java.io.tmpdir")+"/step1assaydriver_selftest_"+System.nanoTime();
		final String work=dir+"/runassay_e2e";
		new java.io.File(work).mkdirs();

		final String focalFamilyId="fZ";
		final String[] neighbors=CandidateCChallengeDbBuilder.frozenKNeighbors(focalFamilyId);
		final String roleFile=work+"/role.tsv";
		writeFile(roleFile, "#artifact_id\tfamily_id\trole\tfold_id\n"
			+CandidateCChallengeDbBuilder.roleRow(focalFamilyId+"__focal_f0", focalFamilyId, "focal_fold_local", 0)
			+CandidateCChallengeDbBuilder.roleRow(focalFamilyId+"__focal_f1", focalFamilyId, "focal_fold_local", 1)
			+CandidateCChallengeDbBuilder.decoyRowsForFamilies(neighbors));
		final String neighborsFile=work+"/neighbors.tsv";
		writeFile(neighborsFile, "#k\t20\n"+CandidateCChallengeDbBuilder.neighborBlock(focalFamilyId, neighbors));
		final String clusterFile=work+"/fZ.cluster.tsv";
		writeFile(clusterFile, "m1\tm1\nm1\tm2\nm3\tm3\nm4\tm4\nm4\tm5\n");
		final String foldManifestFile=work+"/fold_manifest.tsv";
		writeFile(foldManifestFile, "#family_id\tgroup_id\tgroup_size\tfold_id\tgroup_hash\n"
			+"fZ\tm1\t2\t0\t"+new String(new char[16]).replace((char)0, 'a')+"\n"+"fZ\tm3\t1\t0\t"+new String(new char[16]).replace((char)0, 'b')+"\n"
			+"fZ\tm4\t2\t1\t"+new String(new char[16]).replace((char)0, 'c')+"\n");
		final String sampleFile=work+"/sample.tsv";
		writeFile(sampleFile, "#rep_id\tselected\tcell_id\tcluster_path\nfZ\t1\tQ1_F1\t"+clusterFile+"\n");
		final String fastaFile=work+"/seqs.fasta";
		writeFile(fastaFile, ">m1\nAAAA\n>m2\nAAAA\n>m3\nCCCC\n>m4\nGGGG\n>m5\nGGGG\n");
		final String artifactManifestFile=work+"/artifact_manifest.tsv";
		writeFile(artifactManifestFile, "#artifact_id\tmethod\tformat\tpath\n");

		final HashMap<String,double[]> canned=new HashMap<String,double[]>();
		canned.put("m1", rowWithOwnHigh(21, 10.0, 2.0));
		canned.put("m2", rowWithOtherHigh(21, 1.0, 9.0, 5));
		canned.put("m3", rowWithTie(21, 7.0, 3));
		canned.put("m4", rowWithOwnHigh(21, 20.0, 1.0));
		canned.put("m5", rowWithOwnHigh(21, 20.0, 1.0));
		final CannedScorer scorer=new CannedScorer(canned);

		final String top1Out=work+"/assay_top1_confusion.tsv";
		final String marginOut=work+"/assay_margin_3a.tsv";
		runAssay(sampleFile, foldManifestFile, roleFile, neighborsFile, artifactManifestFile, fastaFile,
			scorer, "A", 5.0, "CannedScorer test v1", top1Out, marginOut);

		final String top1Text=readWholeFile(top1Out);
		if(!top1Text.contains("Q1_F1\tfZ\tA\t5\t3\t1\t1\t")){
			throw new RuntimeException("SELFTEST FAILED: runAssay's top1 output missing the expected "
				+"family row -- got:\n"+top1Text);
		}
		for(final String[] pair : new String[][]{
				{"sample_manifest", sampleFile}, {"fold_manifest", foldManifestFile}, {"role_manifest", roleFile},
				{"neighbors", neighborsFile}, {"artifact_manifest", artifactManifestFile}, {"member_seqs_fasta", fastaFile}}){
			final String expectedHash=sha256File(pair[1]);
			if(!top1Text.contains("#input_"+pair[0]+"\t"+pair[1]+"\t"+expectedHash)){
				throw new RuntimeException("SELFTEST FAILED: runAssay's top1 output missing/wrong "
					+"independently-recomputed provenance for '"+pair[0]+"' -- expected hash "+expectedHash
					+", got:\n"+top1Text);
			}
		}
		if(!top1Text.contains("#input_cluster_tsv\t"+clusterFile+"\t"+sha256File(clusterFile))){
			throw new RuntimeException("SELFTEST FAILED: runAssay's top1 output missing/wrong cluster.tsv "
				+"provenance -- got:\n"+top1Text);
		}
		if(!top1Text.contains("#scorer\tCannedScorer test v1")){
			throw new RuntimeException("SELFTEST FAILED: runAssay's top1 output missing the scorer "
				+"description -- got:\n"+top1Text);
		}
		if(!top1Text.contains("#expected_q_x_challenge_cells\t105")
				|| !top1Text.contains("#observed_q_x_challenge_cells\t105")){
			throw new RuntimeException("SELFTEST FAILED: runAssay's top1 output missing/wrong Q x21 cell "
				+"counts (expected 105 = 5 held-out queries x 21 artifacts each) -- got:\n"+top1Text);
		}
		final String marginText=readWholeFile(marginOut);
		if(!marginText.contains("#cell_pass_fraction_pooled\tQ1_F1\tA\t"+formatDecimal(0.6)+"\t5")){
			throw new RuntimeException("SELFTEST FAILED: runAssay's margin output pooled line wrong -- "
				+"got:\n"+marginText);
		}
	}

	/** CP-01: proves the selected sample row's OWN cluster_path is what feeds held-out expansion --
	 * there is no API surface anywhere in this file that could substitute a different ("alien")
	 * cluster.tsv, so this is a structural-closure proof rather than a runtime-rejection test. */
	static void testCP01ClusterPathBoundToSelectedRow() throws Exception{
		final String dir=System.getProperty("java.io.tmpdir")+"/step1assaydriver_selftest_"+System.nanoTime();
		final String work=dir+"/cp01";
		new java.io.File(work).mkdirs();
		final String correctCluster=work+"/correct.cluster.tsv";
		writeFile(correctCluster, "m1\tm1\nm1\tm2\n");
		final String alienCluster=work+"/alien.cluster.tsv";
		writeFile(alienCluster, "m1\tXXXX\n");//deliberately wrong members -- must never be reachable
		final String foldManifestFile=work+"/fold_manifest.tsv";
		writeFile(foldManifestFile, "#family_id\tgroup_id\tgroup_size\tfold_id\tgroup_hash\nfZ\tm1\t2\t0\t"
			+new String(new char[16]).replace((char)0, 'a')+"\n");
		final String sampleFile=work+"/sample.tsv";
		writeFile(sampleFile, "#rep_id\tselected\tcell_id\tcluster_path\nfZ\t1\tQ1_F1\t"+correctCluster+"\n");
		final ArrayList<String> heldOut=loadHeldOutMemberIds(
			findSelectedFamilyRow(sampleFile, "fZ").clusterTsvPath, foldManifestFile, "fZ", 0);
		if(!heldOut.equals(new ArrayList<String>(java.util.Arrays.asList("m1", "m2")))){
			throw new RuntimeException("SELFTEST FAILED: expected held-out [m1, m2] from the DECLARED "
				+"cluster_path, got "+heldOut+" -- possible cluster_path binding leak (CP-01).");
		}
	}

	static void testCP02DuplicateNeighborFamilyRejected() throws Exception{
		final String dir=System.getProperty("java.io.tmpdir")+"/step1assaydriver_selftest_"+System.nanoTime();
		final String work=dir+"/cp02";
		new java.io.File(work).mkdirs();
		final StringBuilder distinct=new StringBuilder();
		for(int r=3; r<=20; r++){distinct.append("fZ\t").append(r).append("\tfN").append(r).append('\n');}
		final String neighborsFile=work+"/neighbors.tsv";
		writeFile(neighborsFile, "#k\t20\nfZ\t1\tfB\nfZ\t2\tfB\n"+distinct);
		expectThrows(()->loadNeighborsForFocal(neighborsFile, "fZ"),
			"duplicate neighbor family at two different ranks (CP-02)");
	}

	static void testCP03MissingSequenceIdRejected() throws Exception{
		final String dir=System.getProperty("java.io.tmpdir")+"/step1assaydriver_selftest_"+System.nanoTime();
		final String work=dir+"/cp03";
		new java.io.File(work).mkdirs();
		final String fastaFile=work+"/seqs.fasta";
		writeFile(fastaFile, ">m1\nAAAA\n");//m2 deliberately absent
		final HashSet<String> needed=new HashSet<String>(java.util.Arrays.asList("m1", "m2"));
		expectThrows(()->extractSequences(fastaFile, needed), "a needed id absent from the FASTA (CP-03)");
	}

	static void testCP04AggregateFamilyHardening() throws Exception{
		final QueryResult goodA=new QueryResult("m1", "fZ", 0, "Q1_F1", Step1AssayMetrics.Top1Bucket.CORRECT, 10, 2, 8, true);
		final QueryResult goodB=new QueryResult("m2", "fZ", 1, "Q1_F1", Step1AssayMetrics.Top1Bucket.CORRECT, 10, 2, 8, true);
		final HashSet<Integer> expectedFolds01=new HashSet<Integer>(java.util.Arrays.asList(0, 1));

		final FamilySummary ok=aggregateFamily("fZ", "Q1_F1",
			new ArrayList<QueryResult>(java.util.Arrays.asList(goodA, goodB)), expectedFolds01);
		if(ok.totalHeldOut!=2){
			throw new RuntimeException("SELFTEST FAILED: baseline aggregateFamily fixture did not aggregate as expected.");
		}

		final QueryResult mixedCell=new QueryResult("m3", "fZ", 0, "Q2_F2", Step1AssayMetrics.Top1Bucket.CORRECT, 10, 2, 8, true);
		expectThrows(()->aggregateFamily("fZ", "Q1_F1",
			new ArrayList<QueryResult>(java.util.Arrays.asList(goodA, mixedCell)), expectedFolds01), "mixed cellId (CP-04)");

		final QueryResult dupMember=new QueryResult("m1", "fZ", 1, "Q1_F1", Step1AssayMetrics.Top1Bucket.CORRECT, 10, 2, 8, true);
		expectThrows(()->aggregateFamily("fZ", "Q1_F1",
			new ArrayList<QueryResult>(java.util.Arrays.asList(goodA, dupMember)), expectedFolds01), "duplicate member id (CP-04)");

		expectThrows(()->aggregateFamily("fZ", "Q1_F1",
			new ArrayList<QueryResult>(java.util.Arrays.asList(goodA)), expectedFolds01),
			"missing expected fold coverage (CP-04)");

		final QueryResult extraFold=new QueryResult("m4", "fZ", 2, "Q1_F1", Step1AssayMetrics.Top1Bucket.CORRECT, 10, 2, 8, true);
		expectThrows(()->aggregateFamily("fZ", "Q1_F1",
			new ArrayList<QueryResult>(java.util.Arrays.asList(goodA, goodB, extraFold)), expectedFolds01),
			"unexpected extra fold coverage (CP-04)");
	}

	static void testCP05to07ArtifactManifestHardening() throws Exception{
		final String dir=System.getProperty("java.io.tmpdir")+"/step1assaydriver_selftest_"+System.nanoTime();
		final String work=dir+"/cp0507";
		new java.io.File(work).mkdirs();
		final String realFasta=work+"/real.fasta";
		writeFile(realFasta, ">a1\nMKVL\n");
		final ArrayList<Query> oneQuery=new ArrayList<Query>(java.util.Arrays.asList(new Query("q1", "AAAA".getBytes())));

		final String noHeader=work+"/no_header.tsv";
		writeFile(noHeader, "art1\tA\tfasta_consensus\t"+realFasta+"\n");
		expectThrows(()->loadArtifactManifest(noHeader), "data row before a valid header (CP-05)");

		final String wrongCol=work+"/wrong_col.tsv";
		writeFile(wrongCol, "#artifact_id\tmethod\tFORMAT_TYPO\tpath\nart1\tA\tfasta_consensus\t"+realFasta+"\n");
		expectThrows(()->loadArtifactManifest(wrongCol), "wrong header column name (CP-05)");

		final String dupHeader=work+"/dup_header.tsv";
		writeFile(dupHeader, "#artifact_id\tmethod\tformat\tpath\n#artifact_id\tmethod\tformat\tpath\n"
			+"art1\tA\tfasta_consensus\t"+realFasta+"\n");
		expectThrows(()->loadArtifactManifest(dupHeader), "duplicate header (CP-05)");

		final String shortRow=work+"/short_row.tsv";
		writeFile(shortRow, "#artifact_id\tmethod\tformat\tpath\nart1\tA\tfasta_consensus\n");
		expectThrows(()->loadArtifactManifest(shortRow), "row with too few fields (CP-06)");

		final String longRow=work+"/long_row.tsv";
		writeFile(longRow, "#artifact_id\tmethod\tformat\tpath\nart1\tA\tfasta_consensus\t"+realFasta+"\textra\n");
		expectThrows(()->loadArtifactManifest(longRow), "row with extra fields (CP-06)");

		final String dupId=work+"/dup_id.tsv";
		writeFile(dupId, "#artifact_id\tmethod\tformat\tpath\nart1\tA\tfasta_consensus\t"+realFasta+"\n"
			+"art1\tA\tfasta_consensus\t"+realFasta+"\n");
		expectThrows(()->loadArtifactManifest(dupId), "duplicate artifact_id (CP-06)");

		final String emptyField=work+"/empty_field.tsv";
		writeFile(emptyField, "#artifact_id\tmethod\tformat\tpath\nart1\tA\t\t"+realFasta+"\n");
		expectThrows(()->loadArtifactManifest(emptyField), "empty required field (CP-06)");

		final String relPath=work+"/rel_path.tsv";
		writeFile(relPath, "#artifact_id\tmethod\tformat\tpath\nart1\tA\tfasta_consensus\treal.fasta\n");
		expectThrows(()->loadArtifactManifest(relPath), "relative artifact path (CP-07)");

		final String dirPath=work+"/dir_path.tsv";
		writeFile(dirPath, "#artifact_id\tmethod\tformat\tpath\nart1\tA\tfasta_consensus\t"+work+"\n");
		expectThrows(()->loadArtifactManifest(dirPath), "directory artifact path (CP-07)");

		final String missingPath=work+"/missing_path.tsv";
		writeFile(missingPath, "#artifact_id\tmethod\tformat\tpath\nart1\tA\tfasta_consensus\t"
			+work+"/does_not_exist.fasta\n");
		expectThrows(()->loadArtifactManifest(missingPath), "nonexistent artifact path (CP-07)");

		try{
			final Path linkPath=Paths.get(work, "link.fasta");
			Files.createSymbolicLink(linkPath, Paths.get(realFasta));
			final String symlinkManifest=work+"/symlink_path.tsv";
			writeFile(symlinkManifest, "#artifact_id\tmethod\tformat\tpath\nart1\tA\tfasta_consensus\t"+linkPath+"\n");
			expectThrows(()->loadArtifactManifest(symlinkManifest), "symlinked artifact path (CP-07)");
		}catch(final java.io.IOException | UnsupportedOperationException e){
			System.err.println("  (skipped symlink sub-case of CP-07 -- unsupported in this environment: "
				+e.getClass().getSimpleName()+")");
		}

		final String twoRecordFasta=work+"/two_record.fasta";
		writeFile(twoRecordFasta, ">a\nAAAA\n>b\nCCCC\n");
		final String twoRecordManifest=work+"/two_record.tsv";
		writeFile(twoRecordManifest, "#artifact_id\tmethod\tformat\tpath\nart1\tA\tfasta_consensus\t"+twoRecordFasta+"\n");
		final FlatConsensusScorer twoRecordScorer=new FlatConsensusScorer(twoRecordManifest);
		expectThrows(()->twoRecordScorer.scoreChallenge(oneQuery, new ArrayList<String>(java.util.Arrays.asList("art1")), "fZ", 0),
			"two-record consensus FASTA (CP-07)");

		final String emptyFasta=work+"/empty.fasta";
		writeFile(emptyFasta, "");
		final String emptyFastaManifest=work+"/empty_fasta.tsv";
		writeFile(emptyFastaManifest, "#artifact_id\tmethod\tformat\tpath\nart1\tA\tfasta_consensus\t"+emptyFasta+"\n");
		final FlatConsensusScorer emptyFastaScorer=new FlatConsensusScorer(emptyFastaManifest);
		expectThrows(()->emptyFastaScorer.scoreChallenge(oneQuery, new ArrayList<String>(java.util.Arrays.asList("art1")), "fZ", 0),
			"empty consensus FASTA (CP-07)");
	}

	static void testCP08And09ChallengeSetBoundary() throws Exception{
		final ArrayList<String> good21=new ArrayList<String>();
		for(int i=0; i<21; i++){good21.add("art"+i);}

		final ArrayList<String> wrongSize=new ArrayList<String>(good21.subList(0, 20));
		expectThrows(()->validateChallengeSetBoundary(wrongSize, "art0", "fZ", 0),
			"challenge set with 20 items, not 21 (CP-08)");
		final ArrayList<String> tooMany=new ArrayList<String>(good21);
		tooMany.add("art21");
		expectThrows(()->validateChallengeSetBoundary(tooMany, "art0", "fZ", 0),
			"challenge set with 22 items, not 21 (CP-08)");

		final ArrayList<String> swapped=new ArrayList<String>(good21);
		swapped.set(0, "art1"); swapped.set(1, "art0");
		expectThrows(()->validateChallengeSetBoundary(swapped, "art0", "fZ", 0),
			"focal artifact not at index 0 (CP-09)");

		validateChallengeSetBoundary(good21, "art0", "fZ", 0);//must NOT throw -- genuinely correct.

		final ArrayList<Query> oneQ=new ArrayList<Query>(java.util.Arrays.asList(new Query("q1", "AAAA".getBytes())));
		final HashMap<String,double[]> matching20=new HashMap<String,double[]>();
		matching20.put("q1", new double[20]);
		expectThrows(()->scoreFocalFold("fZ", 0, "Q1_F1", wrongSize, oneQ, new CannedScorer(matching20), 5.0),
			"scoreFocalFold rejects a self-consistent-but-wrong-size challenge set (CP-08)");
	}

	static void testCP10ProvenanceMismatchRejected() throws Exception{
		final LinkedHashMap<String,String> emptyMap=new LinkedHashMap<String,String>();
		final ProvenanceHeader mismatched=new ProvenanceHeader(emptyMap, emptyMap, emptyMap, emptyMap, "test scorer", 10, 9);
		expectThrows(()->appendProvenanceHeader(new ByteBuilder(), mismatched),
			"expected/observed Q x21 cell count mismatch (CP-10)");
	}

	/** Sayu's finding, 2026-09-02: appendProvenanceHeader must also reject an empty scorer description. */
	static void testEmptyScorerDescriptionRejected() throws Exception{
		final LinkedHashMap<String,String> emptyMap=new LinkedHashMap<String,String>();
		final ProvenanceHeader blank=new ProvenanceHeader(emptyMap, emptyMap, emptyMap, emptyMap, "   ", 10, 10);
		expectThrows(()->appendProvenanceHeader(new ByteBuilder(), blank), "blank/whitespace-only scorer description");
		final String dir=System.getProperty("java.io.tmpdir")+"/step1assaydriver_selftest_"+System.nanoTime();
		new java.io.File(dir).mkdirs();
		final String sampleFile=dir+"/sample.tsv";
		writeFile(sampleFile, "#rep_id\tselected\tcell_id\tcluster_path\nfZ\t1\tQ1_F1\t"+dir+"/x.cluster.tsv\n");
		expectThrows(()->runAssay(sampleFile, sampleFile, sampleFile, sampleFile, sampleFile, sampleFile,
			new CannedScorer(new HashMap<String,double[]>()), "A", 5.0, "", "out1.tsv", "out2.tsv"),
			"runAssay itself rejects an empty scorerDescription before doing any work");
	}

	/** Sayu's finding, 2026-09-02: a single-record FASTA with a header but zero residue bytes must be
	 * rejected, not silently treated as a valid (empty) consensus. */
	static void testZeroResidueFastaRejected() throws Exception{
		final String dir=System.getProperty("java.io.tmpdir")+"/step1assaydriver_selftest_"+System.nanoTime();
		new java.io.File(dir).mkdirs();
		final String emptyResidueFasta=dir+"/empty_residue.fasta";
		writeFile(emptyResidueFasta, ">a1\n");
		expectThrows(()->readSingleSequenceFasta(emptyResidueFasta),
			"single FASTA record with a header but zero residue bytes");
	}

	/** Sayu's finding, 2026-09-02: validateSafePath must reject non-regular files (FIFOs), not just
	 * directories -- a blocking-read FIFO must never reach ByteFile. Skips gracefully if this
	 * environment can't create a FIFO (mkfifo unavailable). */
	static void testNonRegularFileRejected() throws Exception{
		final String dir=System.getProperty("java.io.tmpdir")+"/step1assaydriver_selftest_"+System.nanoTime();
		new java.io.File(dir).mkdirs();
		final String fifoPath=dir+"/a.fifo";
		try{
			final Process p=new ProcessBuilder("mkfifo", fifoPath).inheritIO().start();
			final int exit=p.waitFor();
			if(exit!=0 || !Files.exists(Paths.get(fifoPath))){
				System.err.println("  (skipped FIFO sub-case of CP-07 -- mkfifo unavailable in this environment)");
				return;
			}
		}catch(final Exception e){
			System.err.println("  (skipped FIFO sub-case of CP-07 -- "+e.getClass().getSimpleName()+")");
			return;
		}
		expectThrows(()->validateSafePath(fifoPath, "test"), "a FIFO (non-regular file) path (CP-07)");
	}

	static void testCP11HashDriftRejected() throws Exception{
		final String dir=System.getProperty("java.io.tmpdir")+"/step1assaydriver_selftest_"+System.nanoTime();
		final String work=dir+"/cp11";
		new java.io.File(work).mkdirs();
		final String file=work+"/drift.tsv";
		writeFile(file, "original content\n");
		final String startHash=sha256File(file);
		final LinkedHashMap<String,String> pathToHash=new LinkedHashMap<String,String>();
		pathToHash.put(file, startHash);
		assertNoPathHashDrift(pathToHash, "test file");//unchanged so far -- must NOT throw.
		writeFile(file, "mutated content\n");
		expectThrows(()->assertNoPathHashDrift(pathToHash, "test file"), "a file mutated after hash capture (CP-11)");

		final LinkedHashMap<String,String> paths=new LinkedHashMap<String,String>();
		paths.put("f", file);
		final LinkedHashMap<String,String> hashes=new LinkedHashMap<String,String>();
		hashes.put("f", startHash);//stale -- the file has since been mutated above.
		expectThrows(()->assertNoHashDrift(paths, hashes, "test file"), "role-keyed drift check also rejects (CP-11)");
	}

	/** Elly's finding, 2026-09-02: runAssay's fold-iteration order must be a canonical ascending sort,
	 * never a raw HashSet iteration -- proven here with fold_manifest.tsv rows PHYSICALLY UNORDERED
	 * (fold 1 declared before fold 0 in the file) and confirming the emitted output is still in
	 * ascending fold_id order, deterministically. */
	static void testFoldProcessingOrderDeterministic() throws Exception{
		final String dir=System.getProperty("java.io.tmpdir")+"/step1assaydriver_selftest_"+System.nanoTime();
		final String work=dir+"/fold_order";
		new java.io.File(work).mkdirs();

		final String focalFamilyId="fZ";
		final String[] neighbors=CandidateCChallengeDbBuilder.frozenKNeighbors(focalFamilyId);
		final String roleFile=work+"/role.tsv";
		writeFile(roleFile, "#artifact_id\tfamily_id\trole\tfold_id\n"
			+CandidateCChallengeDbBuilder.roleRow(focalFamilyId+"__focal_f0", focalFamilyId, "focal_fold_local", 0)
			+CandidateCChallengeDbBuilder.roleRow(focalFamilyId+"__focal_f1", focalFamilyId, "focal_fold_local", 1)
			+CandidateCChallengeDbBuilder.decoyRowsForFamilies(neighbors));
		final String neighborsFile=work+"/neighbors.tsv";
		writeFile(neighborsFile, "#k\t20\n"+CandidateCChallengeDbBuilder.neighborBlock(focalFamilyId, neighbors));
		final String clusterFile=work+"/fZ.cluster.tsv";
		writeFile(clusterFile, "m1\tm1\nm1\tm2\nm3\tm3\nm4\tm4\nm4\tm5\n");
		//fold_manifest rows deliberately PHYSICALLY UNORDERED: fold 1 (m4) declared BEFORE fold 0 (m1,m3).
		final String foldManifestFile=work+"/fold_manifest.tsv";
		writeFile(foldManifestFile, "#family_id\tgroup_id\tgroup_size\tfold_id\tgroup_hash\n"
			+"fZ\tm4\t2\t1\t"+new String(new char[16]).replace((char)0, 'c')+"\n"+"fZ\tm1\t2\t0\t"+new String(new char[16]).replace((char)0, 'a')+"\n"
			+"fZ\tm3\t1\t0\t"+new String(new char[16]).replace((char)0, 'b')+"\n");
		final String sampleFile=work+"/sample.tsv";
		writeFile(sampleFile, "#rep_id\tselected\tcell_id\tcluster_path\nfZ\t1\tQ1_F1\t"+clusterFile+"\n");
		final String fastaFile=work+"/seqs.fasta";
		writeFile(fastaFile, ">m1\nAAAA\n>m2\nAAAA\n>m3\nCCCC\n>m4\nGGGG\n>m5\nGGGG\n");
		final String artifactManifestFile=work+"/artifact_manifest.tsv";
		writeFile(artifactManifestFile, "#artifact_id\tmethod\tformat\tpath\n");

		final HashMap<String,double[]> canned=new HashMap<String,double[]>();
		for(final String id : new String[]{"m1", "m2", "m3", "m4", "m5"}){canned.put(id, rowWithOwnHigh(21, 10.0, 2.0));}
		final CannedScorer scorer=new CannedScorer(canned);

		final String top1Out=work+"/assay_top1_confusion.tsv";
		final String marginOut=work+"/assay_margin_3a.tsv";
		runAssay(sampleFile, foldManifestFile, roleFile, neighborsFile, artifactManifestFile, fastaFile,
			scorer, "A", 5.0, "order test", top1Out, marginOut);

		final String marginText=readWholeFile(marginOut);
		final ArrayList<Integer> seenFoldIds=new ArrayList<Integer>();
		for(final String line : marginText.split("\n")){
			if(line.isEmpty() || line.startsWith("#")){continue;}
			seenFoldIds.add(Integer.parseInt(line.split("\t")[2]));//column 2 = fold_id
		}
		final ArrayList<Integer> expectedAscending=new ArrayList<Integer>(seenFoldIds);
		java.util.Collections.sort(expectedAscending);
		if(!seenFoldIds.equals(expectedAscending)){
			throw new RuntimeException("SELFTEST FAILED: fold processing order is not canonically ascending "
				+"despite fold_manifest.tsv listing fold 1 before fold 0 -- got fold_id sequence "+seenFoldIds
				+", expected ascending "+expectedAscending+".");
		}
		if(!seenFoldIds.equals(new ArrayList<Integer>(java.util.Arrays.asList(0, 0, 0, 1, 1)))){
			throw new RuntimeException("SELFTEST FAILED: expected exactly [0,0,0,1,1] (3 fold-0 members "
				+"then 2 fold-1 members), got "+seenFoldIds+".");
		}
	}

	/** Sayu's finding, 2026-09-02: provenance must tie to the actual bytes a real scorer consumes, not
	 * merely to a manifest entry. Drives runAssay with a REAL FlatConsensusScorer over 21 real
	 * artifact FASTA files and confirms every one is independently-recomputed-hash-verified in the
	 * §7 provenance header (never checked via CannedScorer, whose consumedFilePaths() is empty by
	 * default). */
	static void testArtifactFileProvenanceViaFlatConsensusScorer() throws Exception{
		final String dir=System.getProperty("java.io.tmpdir")+"/step1assaydriver_selftest_"+System.nanoTime();
		final String work=dir+"/flatconsensus_provenance";
		new java.io.File(work).mkdirs();

		final String focalFamilyId="fZ";
		final String[] neighbors=CandidateCChallengeDbBuilder.frozenKNeighbors(focalFamilyId);
		final String roleFile=work+"/role.tsv";
		writeFile(roleFile, "#artifact_id\tfamily_id\trole\tfold_id\n"
			+CandidateCChallengeDbBuilder.roleRow(focalFamilyId+"__focal_f0", focalFamilyId, "focal_fold_local", 0)
			+CandidateCChallengeDbBuilder.decoyRowsForFamilies(neighbors));
		final String neighborsFile=work+"/neighbors.tsv";
		writeFile(neighborsFile, "#k\t20\n"+CandidateCChallengeDbBuilder.neighborBlock(focalFamilyId, neighbors));
		final String clusterFile=work+"/fZ.cluster.tsv";
		writeFile(clusterFile, "m1\tm1\n");
		final String foldManifestFile=work+"/fold_manifest.tsv";
		writeFile(foldManifestFile, "#family_id\tgroup_id\tgroup_size\tfold_id\tgroup_hash\nfZ\tm1\t1\t0\t"
			+new String(new char[16]).replace((char)0, 'a')+"\n");
		final String sampleFile=work+"/sample.tsv";
		writeFile(sampleFile, "#rep_id\tselected\tcell_id\tcluster_path\nfZ\t1\tQ1_F1\t"+clusterFile+"\n");
		final String fastaFile=work+"/seqs.fasta";
		writeFile(fastaFile, ">m1\nMKVLACDEFGHIKLMNPQRSTVWY\n");

		final ArrayList<String> artifactIds=new ArrayList<String>();
		artifactIds.add(focalFamilyId+"__focal_f0");
		for(final String n : neighbors){artifactIds.add(n+"__decoy");}
		final StringBuilder manifest=new StringBuilder("#artifact_id\tmethod\tformat\tpath\n");
		final ArrayList<String> artifactPaths=new ArrayList<String>();
		for(int i=0; i<artifactIds.size(); i++){
			final String p=work+"/artifact_"+i+".fasta";
			writeFile(p, ">a"+i+"\nACDEFGHIKLMNPQRSTVWY\n");
			artifactPaths.add(p);
			manifest.append(artifactIds.get(i)).append("\tA\tfasta_consensus\t").append(p).append('\n');
		}
		final String artifactManifestFile=work+"/artifact_manifest.tsv";
		writeFile(artifactManifestFile, manifest.toString());

		final FlatConsensusScorer scorer=new FlatConsensusScorer(artifactManifestFile);
		final String top1Out=work+"/top1.tsv";
		final String marginOut=work+"/margin.tsv";
		runAssay(sampleFile, foldManifestFile, roleFile, neighborsFile, artifactManifestFile, fastaFile,
			scorer, "A", 5.0, "FlatConsensusScorer real-file test", top1Out, marginOut);

		final String top1Text=readWholeFile(top1Out);
		for(final String p : artifactPaths){
			final String expectedHash=sha256File(p);
			if(!top1Text.contains("#input_scorer_file\t"+p+"\t"+expectedHash)){
				throw new RuntimeException("SELFTEST FAILED: missing/wrong scorer-consumed-file provenance "
					+"for '"+p+"' -- expected hash "+expectedHash+", got:\n"+top1Text);
			}
		}
	}

	static void testCP12PooledPassFractionCountedDirectly() throws Exception{
		final QueryResult pass=new QueryResult("m1", "fZ", 0, "Q1_F1", Step1AssayMetrics.Top1Bucket.CORRECT, 10, 2, 8, true);
		final QueryResult fail1=new QueryResult("m2", "fZ", 0, "Q1_F1", Step1AssayMetrics.Top1Bucket.CONFUSED, 1, 9, -8, false);
		final QueryResult fail2=new QueryResult("m3", "fZ", 0, "Q1_F1", Step1AssayMetrics.Top1Bucket.CONFUSED, 1, 9, -8, false);
		final ArrayList<QueryResult> raw=new ArrayList<QueryResult>(java.util.Arrays.asList(pass, fail1, fail2));
		//deliberately WRONG passFraction (0.5) vs the raw truth (1/3) -- reproduces Sayu's exact CP-12 fixture.
		final FamilySummary lying=new FamilySummary("fZ", "Q1_F1", 3, 1, 2, 0, 1.0/3, 2.0/3, 0.0, 0.5);
		final ArrayList<FamilySummary> summaries=new ArrayList<FamilySummary>(java.util.Arrays.asList(lying));
		final String dir=System.getProperty("java.io.tmpdir")+"/step1assaydriver_selftest_"+System.nanoTime();
		new java.io.File(dir).mkdirs();
		final LinkedHashMap<String,String> empty=new LinkedHashMap<String,String>();
		final ProvenanceHeader prov=new ProvenanceHeader(empty, empty, empty, empty, "test", 63, 63);
		final String out=dir+"/margin_cp12.tsv";
		writeMargin3aTsv(out, raw, summaries, "A", 5.0, prov);
		final String text=readWholeFile(out);
		if(!text.contains("#cell_pass_fraction_pooled\tQ1_F1\tA\t"+formatDecimal(1.0/3)+"\t3")){
			throw new RuntimeException("SELFTEST FAILED: CP-12 -- pooled pass-fraction line wrong (must be "
				+"1/3 from raw counts, never Math.round(0.5*3)=2) -- got:\n"+text);
		}
	}

	static double[] rowWithOwnHigh(final int n, final double ownVal, final double otherVal){
		final double[] row=new double[n];
		java.util.Arrays.fill(row, otherVal);
		row[0]=ownVal;
		return row;
	}

	static double[] rowWithOtherHigh(final int n, final double ownVal, final double otherHigh, final int otherIndex){
		final double[] row=new double[n];
		java.util.Arrays.fill(row, 2.0);
		row[0]=ownVal;
		row[otherIndex]=otherHigh;
		return row;
	}

	static double[] rowWithTie(final int n, final double tiedVal, final int otherTiedIndex){
		final double[] row=new double[n];
		java.util.Arrays.fill(row, 1.0);
		row[0]=tiedVal;
		row[otherTiedIndex]=tiedVal;
		return row;
	}

	interface ThrowingRunnable { void run() throws Exception; }

	static void expectThrows(final ThrowingRunnable r, final String label) throws Exception{
		try{
			r.run();
			throw new RuntimeException("SELFTEST FAILED: expected an exception for '"+label+"' but none was thrown.");
		}catch(RuntimeException expected){
			//correct
		}
	}
}
