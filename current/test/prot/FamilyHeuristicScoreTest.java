package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import ddl.DDLFamilyHoldout;
import prot.DigestSuffix;
import prot.FamilyCalibrate;
import prot.FamilyHeuristicRanking;
import prot.FamilyHeuristicScore;
import prot.FamilyEval;
import prot.FamilyScoringPackage;
import prot.FamilyTailCascade;

/**
 * Checks caller-side fourth-root score binding, metadata pins and supported arms.
 * Builds a finite 61-family package, compares direct and CLI scoring, and retains
 * generated evidence. The synthetic sequences test software contracts, not accuracy.
 * Run once per JVM: the assertion count is shared and is not reset by main.
 * Documentation and formatting by Fischl; original authorship is unresolved.
 */
public final class FamilyHeuristicScoreTest{

	/** Runs score, query-length, roster and changed-input checks in a fresh directory.
	 * Prints the success marker only after all checks; generated files remain for inspection.
	 * @param args One nonexistent fixture directory whose parent already exists */
	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh fixture directory");
		final Path root=Paths.get(args[0]);
		check(!Files.exists(root), "Do not overwrite earlier score evidence");
		Files.createDirectory(root);
		final Path config=prepare(root), medians=root.resolve("ranked_medians.tsv"),
			hydro=root.resolve("family_hydro.tsv");
		final FamilyTailCascade.LoadedClassifier base=new FamilyScoringPackage(config).load(2);
		final FamilyHeuristicScore h0=new FamilyHeuristicScore(base, medians, pin(medians), null, null, 100, 149,
			"hitfourth");
		final double[] compactMedians=new double[base.familyCount()], compactHydro=new double[base.familyCount()];
		Arrays.fill(compactMedians, 200);
		Arrays.fill(compactHydro, 1.8);
		compactMedians[0]=100;
		final FamilyHeuristicScore compact=new FamilyHeuristicScore(base, compactMedians, compactHydro,
			100, 149, "hitfourth_h020_c150");
		check(compact.usesHydrophobicity() && compact.armName().equals("hitfourth_h020_c150"),
			"In-memory compact H/C arm binds selected medians and hydrophobicities");
		equal(compact.score(100, 0, 10, .5, 1.8), FamilyHeuristicRanking.score(16, 10, .5, 100, 100, 1));
		reject(()->new FamilyHeuristicScore(base, compactMedians, null, 100, 149, "hitfourth_h020"));
		reject(()->new FamilyHeuristicScore(base, compactMedians, compactHydro, 100, 149, "hitfourth"));
		check(!h0.usesHydrophobicity() && h0.armName().equals("hitfourth"), "H0 declares no query KD dependency");
		equal(h0.score(100, 0, 10, .5, -4.5), FamilyHeuristicRanking.score(6, 10, .5, 100, 100));
		equal(h0.score(100, 1, 10, .5, 4.5), FamilyHeuristicRanking.score(6, 10, .5, 100, 200));

		for(String arm:new String[]{"hitfourth_h005", "hitfourth_h010", "hitfourth_h020", "hitfourth_h040",
				"hitfourth_h080", "hitfourth_h160", "hitfourth_h320", "hitfourth_h020_c050",
				"hitfourth_h020_c150", "hitfourth_h020_c200"}){
			final FamilyHeuristicScore score=new FamilyHeuristicScore(base, medians, pin(medians), hydro, pin(hydro),
				100, 149, arm);
			check(score.usesHydrophobicity() && score.armName().equals(arm),
				"Selected H/cosine arm is exposed exactly: "+arm);
			equal(score.score(100, 0, 10, .5, 1.8), FamilyHeuristicRanking.score(armId(arm), 10, .5, 100, 100, 1));
			boolean missing=false;
			try{score.score(100, 0, 10, .5, Double.NaN);}catch(IllegalArgumentException expected){missing=true;}
			check(missing, "Enabled H rejects an undefined caller KD mean: "+arm);
		}

		final FamilyHeuristicScore shortScore=new FamilyHeuristicScore(base, medians, pin(medians), null, null, 40, 99,
			"hitfourth");
		shortScore.validateQueryLength(99);
		equal(shortScore.score(40, 0, 10, .5, 1.8), FamilyHeuristicRanking.score(6, 10, .5, 40, 100));
		reject(()->shortScore.validateQueryLength(100));
		reject(()->shortScore.score(100, 0, 10, .5, Double.NaN));
		reject(()->new FamilyHeuristicScore(base, medians, pin(medians), hydro, pin(hydro), 100, 149, "hitfourth"));
		reject(()->new FamilyHeuristicScore(base, medians, pin(medians), null, null, 100, 149, "hitfourth_h020"));
		reject(()->new FamilyHeuristicScore(base, medians, pin(medians), hydro, "00000000000000000000", 100, 149,
			"hitfourth_h020"));
		reject(()->new FamilyHeuristicScore(base, medians, pin(medians), null, null, 100, 149, "sumfourth"));
		reject(()->new FamilyHeuristicScore(base, medians, pin(medians), null, null, 0, 200, "hitfourth"));
		reject(()->new FamilyHeuristicScore(base, medians, pin(medians), null, null, 1, 201, "hitfourth"));
		final FamilyHeuristicScore fullRange=new FamilyHeuristicScore(base, medians, pin(medians), hydro, pin(hydro), 1,
			200, "hitfourth_h020_c150");
		for(int length:new int[]{1, 20, 40, 99, 100, 149, 150, 199, 200}){
			fullRange.validateQueryLength(length);
			equal(fullRange.score(length, 0, 10, .5, 1.8),
				10*Math.sqrt(Math.sqrt(Math.min(length, 100)/(double)Math.max(length, 100)))+.75+.2);
		}
		reject(()->fullRange.validateQueryLength(0));
		reject(()->fullRange.validateQueryLength(201));
		checkEvaluator(root, config, medians, hydro);

		final Path badMedian=root.resolve("bad_medians.tsv");
		write(badMedian, text(medians).replace("f1\t1\t200\t200\t200\t0\n", "f1\t2\t200\t200\t200\t0\n"));
		reject(()->new FamilyHeuristicScore(base, badMedian, pin(badMedian), null, null, 100, 149, "hitfourth"));
		final Path badHydro=root.resolve("bad_hydro.tsv");
		write(badHydro, text(hydro).replace("0\tf0\t1\t1.800000000000\n", "0\tf0\t1\t4.600000000000\n"));
		reject(()->new FamilyHeuristicScore(base, medians, pin(medians), badHydro, pin(badHydro), 100, 149,
			"hitfourth_h020"));
		final Path duplicateHydro=root.resolve("duplicate_hydro.tsv");
		write(duplicateHydro, text(hydro)+"61\tf0\t0\tNA\n");
		reject(()->new FamilyHeuristicScore(base, medians, pin(medians), duplicateHydro, pin(duplicateHydro), 100, 149,
			"hitfourth_h020"));
		final FamilyHeuristicScore pinned=new FamilyHeuristicScore(base, medians, pin(medians), null, null, 100, 149,
			"hitfourth");
		write(medians, text(medians)+"\n");
		reject(()->pinned.verify());
		System.out.println("FAMILY_HEURISTIC_SCORE_TEST_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------     Evaluator comparison     ----------------*/
	/*--------------------------------------------------------------*/

	/** The ordinary evaluator must bind the same explicit production policy and retain shorter capped tails. */
	private static void checkEvaluator(Path root, Path config, Path medians, Path hydro) throws Exception{
		final Path input=root.resolve("mixed_queries.faa"), out=root.resolve("heuristic_evaluation");
		final String[] queries={repeat('A', 20), repeat('A', 100), repeat('A', 200)+"*"};
		write(input, ">short\n"+queries[0]+"\n>middle\n"+queries[1]+"\n>full\n"+queries[2]+"\n");
		final FamilyTailCascade.LoadedClassifier loaded=new FamilyScoringPackage(config).loadWithTop(1, 100);
		final FamilyTailCascade.Worker direct=loaded.withHeuristicRanking(new FamilyHeuristicScore(
			loaded, medians, pin(medians), hydro, pin(hydro), 1, 200, "hitfourth_h020_c150")).newWorker();
		int accepted=0;
		for(String query:queries){
			if(direct.assign("oracle", query.getBytes(StandardCharsets.US_ASCII)).assigned()){accepted++;}
		}
		FamilyEval.main(new String[]{"in="+config, "pos="+input, "out="+out, "top=100", "heuristic=hitfourth_h020_c150",
			"lengthmedians="+medians, "lengthmedianssha80="+pin(medians),
			"familyhydro="+hydro, "familyhydrosha80="+pin(hydro),
			"minqueryaa=1", "maxqueryaa=200", "t=2"});
		final String[] row=Files.readAllLines(out.resolve("summary.tsv"), StandardCharsets.US_ASCII).get(1).split("\t");
		check(row[2].equals("3") && Integer.parseInt(row[3])==accepted,
			"Ordinary CLI must retain all lengths and match the direct pinned production heuristic");
		final String policy=text(out.resolve("evaluation_policy.tsv"));
		for(String entry:new String[]{"top\t100\n", "heuristic\thitfourth_h020_c150\n", "cosine_weight\t1.5\n",
			"hydrophobicity_weight\t0.2\n", "raw_hit_admission\tddl_ge_reject\n", "no_backfill\ttrue\n"}){
			check(policy.contains(entry), "Actual evaluator policy must record "+entry);
		}
	}

	/** Finds the exact public arm name or fails the fixture if it is absent. */
	private static int armId(String arm){
		for(int i=0; i<FamilyHeuristicRanking.ALL_ARMS; i++){
			if(FamilyHeuristicRanking.armName(i).equals(arm)){return i;}
		}
		throw new AssertionError(arm);
	}

	/*--------------------------------------------------------------*/
	/*----------------     Fixture preparation      ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates 61 one-member families: one 100-aa and sixty 200-aa alanine sequences.
	 * Writes literal length/hydrophobicity metadata, prepares one calibration round
	 * against a synthetic cysteine control, and returns the completed classifier config.
	 * Existing paths within root may be overwritten by the fixture writer. */
	private static Path prepare(Path root)throws Exception{
		final Path families=root.resolve("families.tsv"), members=root.resolve("members.tsv"),
			chunks=root.resolve("chunks.tsv");
		final Path consensus=root.resolve("consensus.faa"), bounds=root.resolve("bounds.tsv"),
			negative=root.resolve("control.faa");
		final Path panels=root.resolve("panels.tsv"), manifest=root.resolve("preparation.config");
		final StringBuilder fs=new StringBuilder(DDLFamilyHoldout.FAMILY_HEADER+"\n");
		final StringBuilder ms=new StringBuilder(DDLFamilyHoldout.MEMBER_HEADER+"\n");
		final StringBuilder cs=new StringBuilder(),
			ls=new StringBuilder("family\tmembers\tlength_p10\tlength_median\tlength_p90\tcount_under_200\n");
		final StringBuilder hs=new StringBuilder("family_rank\tfamily_id\tN\thydrophobicity_mean\n");
		for(int i=0; i<61; i++){
			final int length=i==0 ? 100 : 200;
			final String seq=repeat('A', length), name="f"+i;
			fs.append(i).append('\t').append(name).append("\t1\n");
			ms.append(i).append('\t').append(name).append("\tg|q").append(i).append("\tg\t").append(i).append('\t')
				.append(seq).append('\n');
			cs.append('>').append(name).append('\n').append(seq).append('\n');
			ls.append(name).append("\t1\t").append(length).append('\t').append(length).append('\t').append(length)
				.append('\t').append(i==0 ? 1 : 0).append('\n');
			hs.append(i).append('\t').append(name).append("\t1\t1.800000000000\n");
		}
		write(families, fs.toString());
		write(members, ms.toString());
		write(consensus, cs.toString());
		write(root.resolve("ranked_medians.tsv"), ls.toString());
		write(root.resolve("family_hydro.tsv"), hs.toString());
		write(chunks, DDLFamilyHoldout.CHUNK_HEADER+"\nb00\t0\t"+members+"\t"+pin(members)+"\t61\t12100\n");
		write(negative, ">fixture\n"+repeat('C', 100)+"\n");
		write(panels, "panel\tpath\tsha80\tqueries\nfixture\t"+negative+"\t"+pin(negative)+"\t1\n");
		final StringBuilder scalar=new StringBuilder("scalar\trows\tmin\tmax\tmin_family\tmax_family\n");
		for(String name:new String[]{"hydrophobicity", "net_charge", "aromaticity", "fraction_small", "pI_tail"}){
			scalar.append(name).append("\t79565\t-20\t20\tlo\thi\n");
		}
		write(bounds, scalar.toString());
		write(manifest, "schema=family_preparation_v1\nfamilies="+families+"\nfamiliessha80="+pin(families)+
			"\nchunks="+chunks+
			"\nchunkssha80="+pin(chunks)+"\nconsensus="+consensus+"\nconsensussha80="+pin(consensus)+
			"\nnegativepanels="+panels+"\nnegativepanelssha80="+pin(panels)+
			"\ncosinepanel=fixture\nexpectedfamilies=61\nexpectedmembers=61\n"+
			"buckets=1\nbatches=1\nmaxreferenceaa=1000\nmaxbucketmembers=64\nmaxbucketbytes=1048576\nmaxpreparedbytes=8388608\n"+
			"maxoutputbytes=8388608\nretention=uniform99\nuppertop=1\nminmemberaa=100\nmaxmemberaa=200\nminnegativeaa=100\n"+
			"maxnegativeaa=200\nk=5\nseed=0\nwidth=1024\nmatrix=blosum62\npositional=t\npartitions=4\noverlap=2\n"+
			"scalarbounds="+bounds+"\nscalarboundssha80="+pin(bounds)+
			"\nscalarscale=0.1\nblosumfloor=0\nscalarz=0\nmonomer=t\ntop=50\n"+
			"ddlreject=1\nddlaccept=1025\ncosine=ranking\nrankingfeatures=combined25\nscalargates=f\n"+
			"identityacceptfloor=30\n"+
			"blosumacceptfloor=0\nindexearlyreject=f\n");
		final Path out=root.resolve("package");
		new FamilyCalibrate(new String[]{"prepare="+manifest, "out="+out, "t=2", "rounds=1",
			"targetfpr=1", "perclass=f", "matrix=blosum62", "verification=sha80"}).process();
		check(Files.exists(out.resolve("PASS")), "The literal package must be prepared before heuristic scoring");
		return out.resolve("classifier.config");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Fixture helpers        ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns count copies of one character for a literal fixture sequence. */
	private static String repeat(char value, int count){
		final char[] text=new char[count];
		Arrays.fill(text, value);
		return new String(text);
	}

	/** Reads a small fixture artifact as ASCII. */
	private static String text(Path path)throws Exception{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}

	/** Creates or truncates a fixture file as ASCII; deliberate corruption cases reuse this writer. */
	private static void write(Path path, String text)throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII));
	}

	/** Computes the legacy SHA80 binding required by this fixture's metadata contract. */
	private static String pin(Path path){return DigestSuffix.file(path.toString());}
	/** Requires an absolute score error strictly below 1e-12. */
	private static void equal(double actual, double expected){
		check(Math.abs(actual-expected)<1e-12, "Unexpected heuristic score: "+actual+" expected="+expected);
	}

	/** Requires IllegalArgumentException; other exceptions propagate instead of counting as rejection. */
	private static void reject(ThrowingRunnable r)throws Exception{
		boolean failed=false;
		try{r.run();}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Invalid score fixture must fail");
	}

	/** Counts each attempted assertion and throws AssertionError on failure, independent of -ea. */
	private static void check(boolean valid, String why){
		checks++;
		if(!valid){throw new AssertionError(why);}
	}

	/** One invalid scoring operation, allowed to propagate checked I/O failures. */
	private interface ThrowingRunnable{void run() throws Exception;}

	/** Attempted assertions; retained for the single invocation in this JVM. */
	private static int checks;
}
