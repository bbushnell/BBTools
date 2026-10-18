package test.prot;

import prot.AAGraph;
import prot.Blosum62;
import prot.FamilyShortlistScores;
import prot.FamilyShortlistSidecar;
import prot.HbmCompetitiveAssign;
import prot.HbmCompetitiveGate;
import prot.HbmPositionModel;
import prot.HbmTailGate;
import prot.ProteinFamilyAssigner;
import prot.ProteinSearcher;
import prot.ProteinSequence;
import prot.ProteinTailAssigner;

import java.nio.charset.StandardCharsets;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Objects;
import java.util.Random;

import dna.AminoAcid;
import shared.Shared;

/** Small real-profile tests of production tail preselection, scheduling and legacy parity. @author Sayu */
public final class ProteinTailSearchTest{
	public static void main(String[] args) throws Exception{
		policies(); lengthBounds(); scheduling(); preselection(); prefilter(); legacy();
		if(args.length>0){workerParity(args);}
		System.err.println("PROTEIN_TAIL_SEARCH_TEST_PASS policies=true capped_core=true seeds=true ties=true strict_reset=true late_winner=true counters=true reset=true legacy_trials=100");
	}
	private static void policies(){
		reject(()->new ProteinTailAssigner.SearchPolicy(null, 0, 0, 0));
		reject(()->fast(-1, 0, 4)); reject(()->fast(0, -0.1, 4)); reject(()->fast(0, 1.1, 4));
		reject(()->fast(0, Double.NaN, 4)); reject(()->fast(0, Double.POSITIVE_INFINITY, 4)); reject(()->fast(0, 0, -1));
		reject(()->new ProteinTailAssigner.SearchPolicy(ProteinTailAssigner.Mode.EXHAUSTIVE, 1, 0, 0));
		reject(()->new ProteinTailAssigner.SearchPolicy(ProteinTailAssigner.Mode.EXHAUSTIVE, 0, 0.1, 0));
		reject(()->new ProteinTailAssigner.SearchPolicy(ProteinTailAssigner.Mode.EXHAUSTIVE, 0, 0, 4));
		check(fast(0, 1, Integer.MAX_VALUE).lookahead==Integer.MAX_VALUE, "Large lookahead must not overflow the bounded stop");
	}
	private static void lengthBounds(){
		final Fixture f=new Fixture(new String[]{repeat('A', 100), repeat('A', 700), repeat('A', 700), repeat('A', 30)}, null);
		final int[] first={0, 0, 0, 0}, last={99, 699, 99, 29};
		final HbmCompetitiveGate gate=new HbmCompetitiveGate(f.ids[0], f.refs, f.models, first, last);
		final HbmTailGate.Policy p=HbmTailGate.EXPERIMENTAL;
		check(!gate.tailLengthCanPass(79, 0, p) && gate.tailLengthCanPass(80, 0, p), "Core coverage uses inclusive float32 .8 boundary");
		check(!gate.tailLengthCanPass(399, 1, p) && gate.tailLengthCanPass(400, 1, p), "Reference core must be capped to500, not full700");
		check(!gate.tailLengthCanPass(500, 2, p), "Core completely left of tail window is impossible");
		check(!gate.tailLengthCanPass(29, 3, p) && gate.tailLengthCanPass(30, 3, p), "Paired minimum is inclusive");
		check(gate.tailLengthCanPass(500, 3, p), "Long clipped query prefix must not invoke full-query mutual coverage");
	}
	private static void scheduling(){
		final Fixture f=new Fixture(new String[]{"AAAAC", "AAAAC", "CCCCC", "CCCCC", "CCCCC", "AAAAA", "CCCCC", "CCCCC", "CCCCC", "CCCCC", "CCCCC"},
			new String[]{"z", "a", "n2", "n3", "n4", "winner", "n6", "n7", "n8", "n9", "n10"});
		final byte[] q=enc("AAAAA");
		final ProteinTailAssigner.SearchScratch bounded=f.work(fast(0, 0, 4));
		final ProteinTailAssigner.Hit h=f.run(bounded, q);
		check(h.assigned() && h.representative.equals("a") && h.visited==5 && h.unvisited==6 && h.earlyStopped,
			"Equal-score ASCII winner must change without extending the lookahead boundary");
		check(h.accepted==2 && h.scored==5 && h.recorded==2 && h.rejected[HbmTailGate.NO_ALIGNMENT]==3,
			"Null candidates consume ranks and one DP fill, not a recorded path");
		check(h.runnerUpScore64==h.score64 && h.bestOtherScore64==h.score64, "Equal-score competitors were lost");
		final ProteinTailAssigner.Hit all=f.run(f.work(fast(0, 0, 50)), q);
		check(all.representative.equals("winner") && all.score64>h.score64 && all.visited==11 && !all.earlyStopped,
			"Bounded search must demonstrate a missed late higher supported score, not claim exhaustive parity");
		final ProteinTailAssigner.Hit zero=f.run(f.work(fast(0, 0, 0)), q);
		check(zero.visited==1 && zero.representative.equals("z") && zero.runnerUpScore64==-1 && zero.bestOtherScore64==-1,
			"Lookahead0 stops at first support and cannot claim an unvisited competitor");
		f.scores[0][5]=f.scores[0][0]-0.5f;
		final ProteinTailAssigner.Hit advance=f.run(bounded, q);
		check(advance.representative.equals("winner") && advance.visited==6, "Strict supported improvement at rank1 must reset to rank6 exclusive");
		final ProteinTailAssigner.Hit no=f.run(bounded, enc("WWWWW"));
		check(!no.assigned() && no.visited==11 && no.recorded==0 && no.accepted==0 && !no.earlyStopped
			&& Float.isNaN(no.identity) && no.runnerUpScore64==-1, "No support must scan all ranks and clear prior evidence");
		check(h.representative.equals("a") && h.rejected[HbmTailGate.NO_ALIGNMENT]==3 && h.identity==80f,
			"Later queries mutated detached hit metrics/rejections");
		final ProteinTailAssigner.Hit shortQuery=f.run(bounded, enc("AA"));
		check(shortQuery.lengthRejected==11 && shortQuery.scored==0 && shortQuery.shortlisted==0 && shortQuery.seedRejected==0,
			"Short query bypass must count all families once without reading stale features");
		check(f.run(bounded, q).representative.equals("winner"), "Scratch failed to recover after empty shortlist");
	}
	private static void preselection(){
		final String[] refs=new String[53]; Arrays.fill(refs, "AAAAA"); refs[0]="AA";
		final Fixture f=new Fixture(refs, null); final byte[] q=enc("AAAAA");
		Arrays.fill(f.seeds[0], 1); f.seeds[0][1]=0;
		ProteinTailAssigner.Hit h=f.run(f.work(fast(1, 0, 50)), q);
		check(h.lengthRejected==1 && h.seedRejected==1 && h.shortlisted==50 && h.visited==50 && h.accepted==50,
			"Preselection must precede top50 and refill from lower-ranked feasible families");
		check(h.rank==2 && h.scored==50 && h.recorded==50, "Stable retrieval ties and single-pass counts differ");
		Arrays.fill(f.scores[0], 1); h=f.run(f.work(fast(1, 0, 0)), q);
		check(h.rank==2, "Equal F4 must preserve canonical local rank order");
		Arrays.fill(f.seeds[0], 0); h=f.run(f.work(fast(1, 0, 4)), q);
		check(h.seedRejected==52 && h.lengthRejected==1 && h.scored==0 && !h.assigned(), "Disjoint mask counts or empty shortlist differ");
		f.qSizes[0]=0; h=f.run(f.work(fast(0, 0, 50)), q);
		check(h.shortlisted==50 && h.assigned(), "Seed0/density0 must not reject empty covering-word queries");
		h=f.run(f.work(fast(0, .5, 50)), q);
		check(h.seedRejected==52 && !h.assigned(), "Positive density cannot pass with no query words");
		f.qSizes[0]=4; f.seeds[0][3]=2; f.seeds[0][4]=1;
		h=f.run(f.work(fast(0, .5, 50)), q);
		check(h.shortlisted==1 && h.rank==3 && h.seedRejected==51, "Density boundary must include exactly2/4 and reject1/4");
		final Fixture a=new Fixture(new String[]{"AAAAA"}, new String[]{"z"});
		final Fixture b=new Fixture(new String[]{"AAAAA"}, new String[]{"a"});
		final ProteinTailAssigner.SearchScratch two=new ProteinTailAssigner.SearchScratch(POLICY, fast(0, 0, 0),
			new HbmCompetitiveGate[]{a.gate, b.gate}, new String[][]{a.ids[0], b.ids[0]}, new int[][]{{1}, {2}});
		h=two.assign(q, new float[][]{{1}, {1}}, new int[][]{{0}, {0}}, new int[]{0, 0});
		check(h.packIndex==0 && h.representative.equals("z"), "F4 ties use canonical pack order, not representative-ID order");
	}
	private static void legacy(){
		final String[] refs=new String[55];
		for(int i=0; i<refs.length; i++){refs[i]=i%3==0 ? "CCCCC" : i%3==1 ? "AAAAC" : "AAAAA";}
		final Fixture f=new Fixture(refs, null);
		//Strong positive profile with incompatible consensus: positive other != supported runner.
		final AAGraph misleading=new AAGraph(f.refs[0], 0);
		for(int i=0; i<100; i++){misleading.addTrace(enc("AAAAA"), 0, "mmmmm".getBytes(StandardCharsets.US_ASCII));}
		f.models[0]=profile(misleading);
		//deriveColumns uses allele frequencies, not observation depth. Pure-A
		//eligible columns would outscore the misleading 100/101 A columns.
		//One C trace plus the weight-1 A scaffold gives eligible columns p(A)=1/2.
		for(int i=2; i<refs.length; i+=3){
			final AAGraph mixed=new AAGraph(f.refs[i], 0);
			mixed.addTrace(enc("CCCCC"), 0, "mmmmm".getBytes(StandardCharsets.US_ASCII));
			f.models[i]=profile(mixed);
		}
		f.rebind();
		final ProteinTailAssigner.SearchScratch work=f.work(new ProteinTailAssigner.SearchPolicy(ProteinTailAssigner.Mode.EXHAUSTIVE, 0, 0, 0));
		final Random random=new Random(20261018);
		for(int trial=0; trial<100; trial++){
			for(int r=0; r<refs.length; r++){f.scores[0][r]=trial==0 ? refs.length-r : random.nextInt(10);}
			final byte[] q=enc(trial%3==0 ? "AAAAA" : trial%3==1 ? "WWWWW" : "AAAAC");
			final ProteinTailAssigner.Hit h=f.run(work, q), expected=oracle(f, q);
			check(h.rank==expected.rank && h.score64==expected.score64 && h.runnerUpScore64==expected.runnerUpScore64
				&& h.bestOtherScore64==expected.bestOtherScore64 && h.recorded==expected.recorded
				&& Arrays.equals(h.rejected, expected.rejected) && h.paired==expected.paired && h.pairedCore==expected.pairedCore
				&& Float.floatToIntBits(h.identity)==Float.floatToIntBits(expected.identity), "Original full50 differential failed trial="+trial);
			check(h.scored==50+h.recorded && h.visited==50 && h.unvisited==0 && h.lengthRejected==0 && h.seedRejected==0,
				"Legacy counters must include all50 matrix fills and each recorded rerun");
			if(trial==0){
				final HbmPositionModel.TailResult misleadingPath=f.gate.alignTail(q, 0, POLICY.maxLength, true);
				check(misleadingPath!=null && f.gate.evaluateTail(q, 0, misleadingPath, POLICY, new HbmTailGate.Metrics())==HbmTailGate.LOW_IDENTITY,
					"Misleading profile must yield positive but sequence-ineligible support");
				check(h.bestOtherScore64==misleadingPath.score && h.bestOtherScore64>h.score64,
					"Fixture must retain a higher positive ineligible other: other="+h.bestOtherScore64+", accepted="+h.score64+", misleading="+misleadingPath.score);
				final ProteinTailAssigner.Hit all=f.run(f.work(fast(0, 0, 50)), q);
				check(all.rank==h.rank && all.score64==h.score64 && all.bestOtherScore64==h.bestOtherScore64
					&& all.runnerUpScore64==h.runnerUpScore64 && all.scored==50, "FAST visited competitors must retain the rejected stronger profile too");
			}
		}
	}
	/** Cheap predicates must agree with the original post-F4 selection, including stale all-short features. */
	private static void prefilter(){
		final Fixture f=new Fixture(new String[]{"AA", "AAAAA", "CCCCC"}, null);
		for(ProteinTailAssigner.SearchPolicy policy : new ProteinTailAssigner.SearchPolicy[]{fast(0, 0, 4), fast(1, 0, 4), fast(0, .5, 4)}){
			final ProteinTailAssigner.SearchScratch work=f.work(policy);
			for(int qSize : new int[]{0, 2, 4}){
				f.qSizes[0]=qSize;
				for(int count=0; count<=qSize; count++){
					Arrays.fill(f.seeds[0], count);
					for(String text : new String[]{"AAAAA", "AA", "WWWWW", "AAAAA"}){
						final byte[] q=enc(text);
						final boolean present=work.hasCandidate(q, f.seeds, f.qSizes);
						final ProteinTailAssigner.Hit h=f.run(work, q);
						check(present==(h.shortlisted>0) && !h.compositionSkipped, "Pre-composition eligibility differs or scratch fixture claims actual bypass");
					}
				}
			}
			f.qSizes[0]=-1; Arrays.fill(f.seeds[0], -1);
			check(!work.hasCandidate(enc("AA"), f.seeds, f.qSizes), "All-short precheck read stale seed features");
		}
		final Fixture shortRefs=new Fixture(new String[]{"AA", "AA"}, null);
		final ProteinTailAssigner.SearchScratch shortWork=shortRefs.work(fast(0, 0, 4));
		check(!shortWork.hasCandidate(enc("AAAAA"), shortRefs.seeds, shortRefs.qSizes)
			&& shortRefs.run(shortWork, enc("AAAAA")).lengthRejected==2, "Long enough query must still reject all impossible reference cores");
		final ProteinTailAssigner.SearchScratch multi=new ProteinTailAssigner.SearchScratch(POLICY, fast(1, 0, 4),
			new HbmCompetitiveGate[]{shortRefs.gate, f.gate}, new String[][]{shortRefs.ids[0], f.ids[0]},
			new int[][]{shortRefs.families[0], f.families[0]});
		final int[][] seeds={{0, 0}, {0, 0, 1}}; final int[] sizes={0, 1};
		check(multi.hasCandidate(enc("AAAAA"), seeds, sizes), "An empty first pack must not hide a survivor in a later pack");
		seeds[1][2]=0;
		check(!multi.hasCandidate(enc("AAAAA"), seeds, sizes), "Multi-pack empty scan retained a prior survivor");
	}

	/** Optional cluster integration: models=path[,path...] uses only verified public package loaders. */
	private static void workerParity(String[] args) throws Exception{
		check(args.length==1 && args[0].startsWith("models="), "Expected only models=comma-separated package directories");
		Shared.setThreads(1);
		final String[] paths=args[0].substring(7).split(",", -1);
		final ProteinFamilyAssigner[] packs=new ProteinFamilyAssigner[paths.length];
		for(int p=0; p<packs.length; p++){packs[p]=ProteinFamilyAssigner.load(Paths.get(paths[p]));}
		int comparisons=workerParity(new ProteinFamilyAssigner[]{packs[0]});
		if(packs.length>1){comparisons+=workerParity(packs);}
		System.err.println("TAIL_COMPOSITION_BYPASS_PARITY_PASS comparisons="+comparisons+" packs="+packs.length+" full_hit_counters=true seed_scores=true reset=true");
	}
	private static int workerParity(ProteinFamilyAssigner[] packs){
		final HbmTailGate.Policy policy=HbmTailGate.EXPERIMENTAL;
		final ArrayList<byte[]> queries=new ArrayList<byte[]>();
		for(ProteinFamilyAssigner pack : packs){
			for(int rank : new int[]{0, pack.familyCount()-1}){
				final byte[] ref=pack.library().refs[rank];
				final int length=Math.min(policy.maxLength, ref.length);
				final byte[] raw=new byte[length];
				for(int i=0; i<length; i++){raw[i]=(byte)AminoAcid.numberToAcid[ref[ref.length-length+i]];}
				queries.add(raw); queries.add(repeat('X', policy.maxLength).getBytes(StandardCharsets.US_ASCII));
				queries.add("AA".getBytes(StandardCharsets.US_ASCII)); queries.add(raw);
			}
		}
		int comparisons=0, skipped=0, computed=0, supported=0;
		for(ProteinTailAssigner.SearchPolicy search : new ProteinTailAssigner.SearchPolicy[]{
				new ProteinTailAssigner.SearchPolicy(ProteinTailAssigner.Mode.EXHAUSTIVE, 0, 0, 0), fast(0, 0, 4), fast(1, 0, 4), fast(0, .5, 4)}){
			final ProteinTailAssigner.Worker worker=new ProteinTailAssigner(policy, search, packs).newWorker();
			final FullScoringOracle oracle=new FullScoringOracle(policy, search, packs);
			ProteinTailAssigner.Hit previous=null, previousExpected=null;
			for(int i=0; i<queries.size(); i++){
				final byte[] raw=queries.get(i); final ProteinSequence query=new ProteinSequence("q"+i, raw);
				final ProteinTailAssigner.Hit expected=oracle.assign(query), actual=worker.assign(query.id, raw);
				equal(expected, actual);
				check(actual.compositionSkipped==(search.mode==ProteinTailAssigner.Mode.FAST && expected.shortlisted==0),
					"Bypass flag must describe an actual zero-survivor FAST Worker return");
				if(previous!=null){equal(previousExpected, previous);}
				previous=actual; previousExpected=expected;
				if(actual.compositionSkipped){skipped++;}else{computed++;}
				if(actual.assigned()){supported++;} comparisons++;
			}
		}
		check(skipped>0 && computed>0 && supported>0, "Real-pack parity must exercise actual bypass, composition and support");
		return comparisons;
	}
	/** Original all-family F4/global-tail path, with independent scratch from the optimized Worker. */
	private static final class FullScoringOracle{
		FullScoringOracle(HbmTailGate.Policy policy, ProteinTailAssigner.SearchPolicy search, ProteinFamilyAssigner[] packs){
			final int n=packs.length;
			sidecars=new FamilyShortlistSidecar[n]; full=new ProteinSearcher.ShortlistScratch[n]; seeded=new ProteinSearcher.ShortlistScratch[n];
			hits=new int[n][]; scores=new float[n][]; composition=new float[n][]; medians=new double[n][]; spreads=new double[n][]; qSizes=new int[n];
			final HbmCompetitiveGate[] gates=new HbmCompetitiveGate[n]; final String[][] ids=new String[n][]; final int[][] families=new int[n][];
			for(int p=0; p<n; p++){
				final HbmCompetitiveAssign.Library lib=packs[p].library(); final FamilyShortlistSidecar sc=lib.sidecar;
				sidecars[p]=sc; full[p]=new ProteinSearcher.ShortlistScratch(sc.dims, sc.nFamilies, 50);
				seeded[p]=new ProteinSearcher.ShortlistScratch(sc.dims, sc.nFamilies, 50);
				hits[p]=full[p].sk; scores[p]=full[p].f4; composition[p]=full[p].compScores; medians[p]=sc.medLen; spreads[p]=sc.sdLnLen;
				gates[p]=lib.gate; ids[p]=lib.ids; families[p]=lib.familyIds;
			}
			work=new ProteinTailAssigner.SearchScratch(policy, search, gates, ids, families);
		}
		ProteinTailAssigner.Hit assign(ProteinSequence query){
			for(int p=0; p<sidecars.length; p++){
				final ProteinSearcher.ShortlistScratch s=seeded[p], original=full[p];
				ProteinSearcher.scoreFamiliesF4(query, sidecars[p], original);
				Arrays.fill(s.compScores, 12345f); Arrays.fill(s.f4, -12345f);
				ProteinSearcher.scoreFamiliesKmers(query, sidecars[p], s);
				check(s.qSize==original.qSize && Arrays.equals(s.sk, original.sk), "Seed-only counts/reset differ from original F4");
				for(int r=0; r<s.sk.length; r++){check(s.compScores[r]==12345f && s.f4[r]==-12345f, "Seed-only seam modified composition or F4");}
				ProteinSearcher.scoreFamiliesF4FromSeeds(query, sidecars[p], s);
				check(Arrays.equals(s.compScores, original.compScores) && Arrays.equals(s.f4, original.f4), "Seeded F4 changed all-family float scores");
				qSizes[p]=original.qSize;
			}
			FamilyShortlistScores.scoreTail(hits, composition, medians, spreads, query.length(), scores);
			return work.assign(query.enc, scores, hits, qSizes);
		}
		final FamilyShortlistSidecar[] sidecars;
		final ProteinSearcher.ShortlistScratch[] full, seeded;
		final int[][] hits; final int[] qSizes;
		final float[][] scores, composition; final double[][] medians, spreads;
		final ProteinTailAssigner.SearchScratch work;
	}
	private static void equal(ProteinTailAssigner.Hit a, ProteinTailAssigner.Hit b){
		check(a.packIndex==b.packIndex && a.rank==b.rank && a.familyId==b.familyId && Objects.equals(a.representative, b.representative)
			&& a.score64==b.score64 && a.runnerUpScore64==b.runnerUpScore64 && a.bestOtherScore64==b.bestOtherScore64
			&& a.paired==b.paired && a.pairedCore==b.pairedCore && a.terminalPaired==b.terminalPaired && a.terminalScore64==b.terminalScore64
			&& a.queryStart==b.queryStart && a.referenceStart==b.referenceStart && a.querySuffix==b.querySuffix && a.referenceSuffix==b.referenceSuffix
			&& Float.floatToIntBits(a.identity)==Float.floatToIntBits(b.identity) && Float.floatToIntBits(a.pairedFraction)==Float.floatToIntBits(b.pairedFraction)
			&& Float.floatToIntBits(a.coreCoverage)==Float.floatToIntBits(b.coreCoverage) && a.recorded==b.recorded && Arrays.equals(a.rejected, b.rejected)
			&& a.shortlisted==b.shortlisted && a.lengthRejected==b.lengthRejected && a.seedRejected==b.seedRejected
			&& a.visited==b.visited && a.scored==b.scored && a.accepted==b.accepted && a.unvisited==b.unvisited && a.earlyStopped==b.earlyStopped,
			"Full-F4 oracle differs from optimized Worker hit, metrics or old counters");
	}
	/** Independent record-all oracle, then score/ID selection; no SearchScratch selection logic. */
	private static ProteinTailAssigner.Hit oracle(Fixture f, byte[] q){
		final int[] p=new int[50], ranks=new int[50], rejected=new int[HbmTailGate.STATUS_COUNT];
		FamilyShortlistScores.top(f.scores, p, ranks);
		final HbmPositionModel.TailResult[] paths=new HbmPositionModel.TailResult[50];
		final boolean[] used=new boolean[50];
		for(int k=0; k<50; k++){paths[k]=f.gate.alignTail(q, ranks[k], POLICY.maxLength, true); if(paths[k]==null){rejected[HbmTailGate.NO_ALIGNMENT]++;}}
		final HbmTailGate.Metrics metric=new HbmTailGate.Metrics(), winner=new HbmTailGate.Metrics();
		int best=-1, runner=-1, recorded=0, bestScore=0;
		for(int step=0; step<50; step++){
			int next=-1;
			for(int k=0; k<50; k++){
				if(used[k] || paths[k]==null){continue;}
				if(next<0 || paths[k].score>paths[next].score || (paths[k].score==paths[next].score && f.ids[0][ranks[k]].compareTo(f.ids[0][ranks[next]])<0)){next=k;}
			}
			if(next<0){break;}
			used[next]=true; recorded++;
			final int status=f.gate.evaluateTail(q, ranks[next], paths[next], POLICY, metric);
			if(status!=HbmTailGate.SUPPORTED){rejected[status]++; continue;}
			if(best<0){best=next; bestScore=paths[next].score; winner.copyFrom(metric);}
			else{runner=paths[next].score; break;}
		}
		int other=-1;
		if(best>=0){for(int k=0; k<50; k++){if(k!=best && paths[k]!=null){other=Math.max(other, paths[k].score);}}}
		final int rank=best<0 ? -1 : ranks[best];
		return new ProteinTailAssigner.Hit(best<0 ? -1 : 0, rank, rank, rank<0 ? null : f.ids[0][rank], bestScore, runner, other,
			best<0 ? null : winner, recorded, rejected);
	}
	private static final class Fixture{
		Fixture(String[] sequences, String[] names){
			final int n=sequences.length; check(n>0, "Fixture needs a family");
			refs=new byte[n][]; models=new HbmPositionModel[n]; ids=new String[][]{new String[n]}; families=new int[][]{new int[n]};
			scores=new float[][]{new float[n]}; seeds=new int[][]{new int[n]}; qSizes=new int[]{1};
			for(int i=0; i<n; i++){
				refs[i]=enc(sequences[i]); models[i]=profile(new AAGraph(refs[i], 0));
				ids[0][i]=names==null ? "f"+(1000+i) : names[i]; families[0][i]=i; scores[0][i]=n-i;
			}
			rebind();
		}
		void rebind(){
			final int[] starts=new int[refs.length], ends=new int[refs.length];
			for(int i=0; i<refs.length; i++){ends[i]=refs[i].length-1;}
			gate=new HbmCompetitiveGate(ids[0], refs, models, starts, ends);
		}
		ProteinTailAssigner.SearchScratch work(ProteinTailAssigner.SearchPolicy search){
			return new ProteinTailAssigner.SearchScratch(POLICY, search, new HbmCompetitiveGate[]{gate}, ids, families);
		}
		ProteinTailAssigner.Hit run(ProteinTailAssigner.SearchScratch work, byte[] q){return work.assign(q, scores, seeds, qSizes);}
		final byte[][] refs; final HbmPositionModel[] models; final String[][] ids; final int[][] families, seeds;
		final float[][] scores; final int[] qSizes; HbmCompetitiveGate gate;
	}
	private static ProteinTailAssigner.SearchPolicy fast(int seeds, double density, int lookahead){return new ProteinTailAssigner.SearchPolicy(ProteinTailAssigner.Mode.FAST, seeds, density, lookahead);}
	private static HbmPositionModel profile(AAGraph graph){final double[] bg=new double[20]; Arrays.fill(bg, .05); return HbmPositionModel.deriveOne(graph, "logodds", .01, true, bg, -4);}
	private static byte[] enc(String s){return Blosum62.encode(s.getBytes(StandardCharsets.US_ASCII), s);}
	private static String repeat(char c, int n){final char[] chars=new char[n]; Arrays.fill(chars, c); return new String(chars);}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private static void reject(Runnable action){try{action.run();}catch(IllegalArgumentException expected){return;} throw new AssertionError("Invalid search policy was accepted");}
	private static final HbmTailGate.Policy POLICY=new HbmTailGate.Policy(10, 3, 1, 2, 1, 40, .6f);
}
