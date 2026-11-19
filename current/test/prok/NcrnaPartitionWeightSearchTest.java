package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Random;

import cardinality.DynamicDemiLog;
import ddl.DDLIndexSmall;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlPartitionKeys;
import prok.NcrnaPartitionWeightSearch;
import prok.NcrnaPartitionWeightSearch.Fit;
import prok.NcrnaPartitionWeightSearch.Metrics;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.OptimizeNcrnaPartitionWeights;
import prok.BuildNcrnaDdlProfiles;
import prok.RankDualNcrnaWindows;
import prok.RankCompactNcrnaWindows;
import shared.Shared;

/** Independent scalar partition counts, fit/judge isolation and native cache replay.
 * @author Ganyu
 */
public final class NcrnaPartitionWeightSearchTest{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaPartitionWeightSearchTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, pair[0].toLowerCase(java.util.Locale.ROOT), pair.length>1 ? pair[1] : null), "Unknown test parameter: "+arg);
			}
			check(parser.out1!=null && !new File(parser.out1).exists(), "Use out=fresh_fixture_directory");
			counts();
			plainCounts();
			largerBudget();
			ranking();
			topKMetrics();
			crossfit();
			check(new File(parser.out1).mkdirs(), "Create only the named fresh synthetic fixture root");
			SummarizeNcrnaLayoutFitsTest.run(parser.out1);
			final String fixture=fixture(parser.out1);
			integration(fixture);
			pp.outstream.println("NCRNA_PARTITION_WEIGHT_SEARCH_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private static void counts(){
		final byte[][] sequences={dna("AAAAAAACGTTGCAACTGAGTTTT"), dna("TTTTGGGCCCAAAGCTGACGT")};
		for(final int width:new int[]{512, 1024}){
			final NcrnaDdlPartitionKeys query=new NcrnaDdlPartitionKeys(7, width);
			final char[][] keys=new char[2][width];
			for(int m=0; m<2; m++){
				final char[] member=query.memberKeys(sequences[m]);
				System.arraycopy(member, 0, keys[m], 0, 3*width/4);
				final char[] consensus=query.consensusKeys(sequences[m]);
				System.arraycopy(consensus, 3*width/4, keys[m], 3*width/4, width/4);
			}
			final DDLIndexSmall index=new DDLIndexSmall(keys);
			final int[][] observed={new int[2], new int[2], new int[2]};
			for(final byte[] bases:new byte[][]{sequences[0], sequences[1], dna("NNNN"), dna(""), sequences[0]}){
				query.build(bases);
				query.countParts(index, observed);
				final int[][] expected=new int[3][2];
				for(int m=0; m<2; m++){
					for(int bucket=0; bucket<width; bucket++){
						for(int k=0; k<query.keyCount(); k++){
							final int packed=query.packedKey(k);
							if((packed>>>16)==bucket && (packed&65535)==keys[m][bucket]){
								final int part=bucket<width/4 ? 0 : bucket<3*width/4 ? 1 : 2;
								expected[part][m]++;
							}
						}
					}
				}
				for(int part=0; part<3; part++){
					for(int m=0; m<2; m++){check(observed[part][m]==expected[part][m], "Posting partition count differs from scalar raw bucket/key membership");}
				}
			}
			for(final int shortK:new int[]{0,8}){
				boolean rejected=false;
				try{new NcrnaDdlPartitionKeys(shortK, width).countParts(index, observed);}
				catch(IllegalArgumentException expected){rejected=true;}
				check(rejected, "Only the original short7 geometry may use its three-part count layout");
			}
			boolean rejected=false;
			try{query.countParts(index, new int[][]{observed[0], observed[0], observed[2]});}
			catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "Shared three-part output rows cannot alias");
			rejected=false;
			try{query.countParts(new DDLIndexSmall(new char[2][width+1]), observed);}
			catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "A three-part codec cannot be applied to a different index width");
		}
	}

	/** Literal equal keys match without a stored winner origin; bin boundaries
	 * alone define partitions. Includes signed postings and maximum ordinals. */
	private static void plainCounts(){
		for(final int width:new int[]{1, 85, 170, 171, 256, 342, 384, 512, 65536}){
			final char[][] profiles=new char[4][width];
			profiles[0][0]=65535;
			profiles[1][0]=65535;
			profiles[2][0]=32768;
			profiles[3][width-1]=1;
			final int[] query={65535, 32768, ((width-1)<<16)|1, 2, 0};
			final int[] totals=new int[4];
			final DDLIndexSmall index=new DDLIndexSmall(profiles);
			index.query(query, 4, totals);
			check(Arrays.equals(totals, new int[]{1, 1, 1, 1}),
				"Equal compressed keys match identically regardless of which word originally won");
			profiles[0][0]=0;
			index.query(query, 4, totals);
			check(totals[0]==1, "Producer scratch mutation cannot change copied postings");
			check(index.pointerSlots()==width && index.populatedBuckets()==(width==1 ? 1 : 2)
				&& index.emptyBuckets()==width-index.populatedBuckets()
				&& index.postingPayloadBytes()==16,
				"Payload consists only of key-derived postings and per-bin references");
			index.query(query, 0, totals);
			check(Arrays.equals(totals, new int[4]), "An empty live prefix clears previous scores");
			if(width>1){
				final int[] first=new int[4];
				final int[] second=new int[4];
				final int[] ends={width/2, width};
				index.query(query, 4, ends, new int[][]{first, second});
				check(Arrays.equals(first, new int[]{1, 1, 1, 0}) && Arrays.equals(second, new int[]{0, 0, 0, 1}),
					"Only the declared bin partition determines the score term");
				boolean rejected=false;
				try{index.query(query, 4, ends, new int[][]{first, first});}
				catch(IllegalArgumentException expected){rejected=true;}
				check(rejected, "Two partition outputs cannot alias");
			}
		}
		final char[][] large=new char[65536][1];
		large[65535][0]=65535;
		final DDLIndexSmall maximum=new DDLIndexSmall(large);
		final int[] scores=new int[65536];
		maximum.query(new int[]{65535}, 1, scores);
		check(scores[65535]==1 && scores[0]==0 && maximum.postingPayloadBytes()==4,
			"Both packed fields preserve unsigned65535, including posting int minus1");
		for(final int bad:new int[]{0, 1<<16|1}){
			boolean rejected=false;
			try{maximum.query(new int[]{bad}, 1, scores);}
			catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "Empty or out-of-layout keys must fail loudly");
		}
		System.err.println("NCRNA_DDL_PLAIN_KEY_INDEX_TEST_PASS");
	}

	private static void largerBudget(){
		final char[][] counts={{10, 1, 1, 10}, {10, 1, 1, 10}};
		final int[][] edges={{0}, {1}};
		final long[][][] masks={{{2}}, {{2}}};
		final NcrnaPartitionWeightSearch engine=new NcrnaPartitionWeightSearch(2, 2, counts, edges, masks,
			new int[]{0, 0}, new int[]{0, 1}, new int[]{0, 0}, 5000, 0);
		final Fit a=engine.fit(0);
		final Fit b=engine.fit(1);
		check(engine.randomDraws==5000 && engine.referenceTerm==0 && a.trials.size()>=5001,
			"The expanded budget must reach the actual scorer, fixing mixed-member scale rather than consensus");
		checkInitialPrefix(a, 2, 0, 5000);
		for(int i=0; i<=5000; i++){
			check(a.trials.get(i)==b.trials.get(i) && a.trials.get(i).weights[0]==1,
				"Both judging directions reuse all5000 initial draws with the declared fixed term");
		}
		check(a.winner.metrics[0].cappedSum==1 && a.winner.metrics[1].cappedSum==1,
			"A two-term mixed sketch can improve the literal held-out ranking");
		final NcrnaPartitionWeightSearch one=new NcrnaPartitionWeightSearch(2, 1,
			new char[][]{{1, 2}, {1, 2}}, edges, masks, new int[]{0, 0}, new int[]{0, 1}, new int[]{0, 0}, 5000, 0);
		final Fit fixed=one.fit(0);
		check(fixed.trials.size()==1 && fixed.winner.normalized()[0]==1,
			"A single pooled match count has no free coefficient and must not manufacture random dimensions");
		boolean rejected=false;
		try{new NcrnaPartitionWeightSearch(2, 2, counts, edges, masks,
			new int[]{0, 0}, new int[]{0, 1}, new int[]{0, 0}, 0, 0);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "A nonpositive search budget must fail before trial generation");
		System.err.println("NCRNA_DDL_EXPANDED_WEIGHT_BUDGET_TEST_PASS");
	}

	/** Pins the prior random stream and its insertion order independently of the
	 * membership structure; the fixed tested seed has no duplicate points. */
	private static void checkInitialPrefix(final Fit fit, final int terms, final int reference, final int draws){
		final double[] expected=new double[terms];
		Arrays.fill(expected, 1);
		check(Arrays.equals(fit.trials.get(0).weights, expected) && fit.trials.get(0).stage.equals("UNIT"),
			"Unit weights remain the first evaluated point after changing duplicate lookup");
		final Random random=new Random(20261113L);
		for(int i=1; i<=draws; i++){
			for(int term=0; term<terms; term++){
				if(term!=reference){expected[term]=StrictMath.exp(StrictMath.log(.1)+StrictMath.log(100)*random.nextDouble());}
			}
			check(Arrays.equals(fit.trials.get(i).weights, expected) && fit.trials.get(i).stage.equals("RANDOM"),
				"Exact seeded weights, fixed scale and draw order must survive faster duplicate lookup at draw "+i);
		}
	}

	private static void ranking(){
		check(NcrnaPartitionWeightSearch.firstRank(new double[]{0, 0, 0}, new long[]{4})==3,
			"Zero-score ties retain original ordinal, including acceptable tail models");
		check(NcrnaPartitionWeightSearch.firstRank(new double[]{3, 2, 5}, new long[]{3})==2,
			"Choose the highest-scoring acceptable model, then its rank among all models");
		check(NcrnaPartitionWeightSearch.firstRank(new double[]{3, 2, 5}, new long[]{0})==-1,
			"No acceptable model remains explicitly unavailable");
		final NcrnaPartitionWeightSearch search=new NcrnaPartitionWeightSearch(2,
			new char[][]{{10, 1, 0, 1, 4, 0}, {0, 1, 0, 0, 2, 0}}, new int[][]{{0}, {0}},
			new long[][][]{{{2}}, {{2}}}, new int[]{0, 0}, new int[]{0, 1}, new int[]{0, 4});
		final int[] ranks=search.ranks(1, 1);
		check(Arrays.equals(ranks, new int[]{1, -1}), "Parent rank is the minimum across windows; missing parents remain -1");
		check(search.evaluate(1, 1).metrics[2].cappedSum==385, "Missing rank costs cap384 in the fitting objective, not zero");
		final double[] many=new double[500];
		final long[] bits=new long[(many.length+63)/64];
		bits[499>>>6]=1L<<(499&63);
		check(NcrnaPartitionWeightSearch.firstRank(many, bits)==500, "Rank itself is not truncated before capped-objective reporting");
	}

	private static void topKMetrics(){
		final int[] expectedRanks={1, 10, 11, 50, 51, -1};
		final char[][] counts=new char[expectedRanks.length][64];
		final int[][] edges=new int[expectedRanks.length][1];
		final long[][][] masks=new long[expectedRanks.length][1][1];
		for(int i=0; i<expectedRanks.length; i++){
			Arrays.fill(counts[i], (char)1);
			edges[i][0]=i;
			if(expectedRanks[i]>0){masks[i][0][0]=1L<<(expectedRanks[i]-1);}
		}
		final NcrnaPartitionWeightSearch engine=new NcrnaPartitionWeightSearch(64, 1, counts, edges, masks,
			new int[6], new int[]{0, 1, 0, 1, 0, 1}, new int[]{0, 1, 0, 0, 1, 0}, 1, 0);
		final Metrics[] metrics=engine.evaluate(new double[]{1}).metrics;
		final Metrics all=metrics[2];
		check(all.parents==6 && all.unavailable==1 && all.availableRankSum==123 && all.cappedSum==507,
			"One-based ordinal ties preserve ranks1/10/11/50/51, while unavailable costs cap384 and has no available rank");
		check(all.within10==2 && all.within50==4,
			"Top-k reporting must include exact ranks10/50 and exclude11/51 and unavailable, independent of the larger cost cap");
		check(all.lowParents==4 && all.lowWithin10==1 && all.lowWithin50==3
			&& all.lowAvailableRankSum==62 && all.lowCappedSum==446,
			"Low-identity top-k counts use the same declared subset and never count its unavailable parent");
		check(metrics[0].within10+metrics[1].within10==all.within10
			&& metrics[0].within50+metrics[1].within50==all.within50
			&& metrics[0].lowWithin50+metrics[1].lowWithin50==all.lowWithin50,
			"The frozen opposite halves must conserve the complete judging population's top-k counts");
		System.err.println("NCRNA_DDL_TOP_K_METRICS_TEST_PASS");
	}

	private static NcrnaPartitionWeightSearch sample(final boolean reverseJudge){
		final char[] a={10, 1, 0, 1, 4, 0};
		final char[] b=reverseJudge ? new char[]{1, 4, 0, 10, 1, 0} : a.clone();
		return new NcrnaPartitionWeightSearch(2, new char[][]{a, b}, new int[][]{{0}, {1}},
			new long[][][]{{{2}}, {{2}}}, new int[]{0, 0}, new int[]{0, 1}, new int[]{0, 0});
	}

	private static void crossfit(){
		final NcrnaPartitionWeightSearch shared=sample(false);
		final Fit original=shared.fit(0);
		final Fit opposite=shared.fit(1);
		final Fit refit=shared.fit(-1);
		checkInitialPrefix(original, 3, 1, NcrnaPartitionWeightSearch.DRAWS);
		check(original.trials!=opposite.trials && original.trials!=refit.trials,
			"Fit-local refinement histories must remain separate even when initial scores are shared");
		for(int i=0; i<=NcrnaPartitionWeightSearch.DRAWS; i++){
			check(original.trials.get(i)==opposite.trials.get(i) && original.trials.get(i)==refit.trials.get(i),
				"Both halves and the final refit must reuse the unit and seeded random results without ranking them again");
		}
		final Fit changedJudge=sample(true).fit(0);
		check(original.winner.shortWeight==changedJudge.winner.shortWeight
			&& original.winner.consensusWeight==changedJudge.winner.consensusWeight
			&& original.trials.size()==changedJudge.trials.size(), "Judge-only changes cannot alter a fit winner or refinement budget");
		for(int i=0; i<original.trials.size(); i++){
			check(original.trials.get(i).shortWeight==changedJudge.trials.get(i).shortWeight
				&& original.trials.get(i).consensusWeight==changedJudge.trials.get(i).consensusWeight,
				"Every searched point must depend only on the seed and fitting half");
		}
		check(original.winner.metrics[0].cappedSum==1 && original.winner.metrics[1].cappedSum==1,
			"The literal fit can improve both independent halves from rank2 to rank1");
		check(changedJudge.winner.metrics[1].cappedSum==2 && changedJudge.trials.get(0).metrics[1].cappedSum==1,
			"A training improvement may hurt judgment and must remain visible");
		check(original.nearTrials>0 && original.range[0]<=original.winner.normalized()[0]
			&& original.range[1]>=original.winner.normalized()[0], "Evaluated normalized one-percent ranges must contain their own optimum");
		final NcrnaPartitionWeightSearch four=new NcrnaPartitionWeightSearch(2, 4,
			new char[][]{{10, 1, 8, 0, 1, 4, 0, 8}, {10, 1, 8, 0, 1, 4, 0, 8}},
			new int[][]{{0}, {1}}, new long[][][]{{{2}}, {{2}}},
			new int[]{0, 0}, new int[]{0, 1}, new int[]{0, 0});
		final Fit fit=four.fit(0);
		double normalizedSum=0;
		for(final double weight:fit.winner.normalized()){normalizedSum+=weight;}
		check(Math.abs(normalizedSum-1)<1e-12 && fit.winner.weights.length==4
			&& fit.winner.weights[1]==1 && fit.winner.metrics[0].cappedSum==1,
			"Four literal score terms are independently weighted; reported weights sum to one");
	}

	private static String fixture(final String root){
		final int exponent=DynamicDemiLog.exponentBits();
		try{
			DynamicDemiLog.setExponent(4);
			final RankSavedNcrnaWindowsTest.Fixture f=new RankSavedNcrnaWindowsTest.Fixture(root, "weight");
			RankDualNcrnaWindowsTest.coherent(f);
			RankDualNcrnaWindowsTest.baseline(f, "baseline.tsv.gz", false, false);
			// The original test has only one NOT_SEEN group. This new witness
			// explicitly makes beta NOT_SEEN before deriving its own frozen halves.
			relabel(f.dir+"/beta.labels.tsv", false);
			relabel(f.dir+"/baseline.tsv.gz", true);
			final RankDualNcrnaWindows.Config c=RankDualNcrnaWindowsTest.config(f, "one", 1);
			RankDualNcrnaWindows.run(c);
			c.identities=f.dir+"/one";
			c.output=f.dir+"/compact_partitions";
			RankCompactNcrnaWindows.runPartitions(c);
			final BuildNcrnaDdlProfiles.Config build=new BuildNcrnaDdlProfiles.Config();
			build.fasta=c.members;
			build.references=c.references;
			build.ledger=c.ledger;
			build.assignments=c.assignments;
			build.census=c.census;
			build.summary=c.summary;
			build.models=2;
			build.members=2;
			for(final int width:new int[]{512, 1024}){
				build.partitionBins=width;
				build.output=f.dir+"/compact_partitions_profiles"+width+".tsv";
				BuildNcrnaDdlProfiles.build(build);
			}
			BuildNcrnaDdlLayoutCountsTest.run(c, f.dir+"/compact_partitions_profiles512.tsv");
			return f.dir;
		}finally{DynamicDemiLog.setExponent(exponent);}
	}

	private static void relabel(final String path, final boolean betaOnly){
		final String saved=path+".original.gz";
		// Rename must preserve compression suffix for ByteFile format detection.
		final String archive=path.endsWith(".gz") ? saved : path+".original.tsv";
		check(new File(path).renameTo(new File(archive)), "Preserve the toy preimage before changing its declared group label");
		try(Input in=new Input(archive); Output out=new Output(path)){
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				final String text=new String(row, StandardCharsets.UTF_8);
				final String changed=!betaOnly || text.startsWith("beta\t") ? text.replace("\tSEEN", "\tNOT_SEEN") : text;
				out.writer.println(changed.getBytes(StandardCharsets.UTF_8));
			}
		}
	}

	private static void integration(final String root){
		final OptimizeNcrnaPartitionWeights.Config first=config(root, "/weight_search");
		OptimizeNcrnaPartitionWeights.run(first);
		final OptimizeNcrnaPartitionWeights.Config replay=config(root, "/weight_replay");
		replay.counts=first.output+".counts.tsv.gz";
		OptimizeNcrnaPartitionWeights.run(replay);
		for(final String suffix:new String[]{".trials.tsv.gz", ".fits.tsv", ".decision.tsv"}){
			try(Input a=new Input(first.output+suffix); Input b=new Input(replay.output+suffix)){
				for(byte[] row; (row=a.reader.nextLine())!=null;){check(Arrays.equals(row, b.reader.nextLine()), "Saved-count replay must reproduce every deterministic trial and fit row");}
				check(b.reader.nextLine()==null, "Replay cannot append rows");
			}
		}
		final String bad=root+"/bad_weight_counts.tsv";
		try(Input in=new Input(replay.counts); Output out=new Output(bad)){
			out.writer.println(in.reader.nextLine());
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){out.writer.println(row);}
		}
		final OptimizeNcrnaPartitionWeights.Config rejected=config(root, "/weight_rejected");
		rejected.counts=bad;
		boolean failed=false;
		try{OptimizeNcrnaPartitionWeights.run(rejected);}catch(IllegalArgumentException expected){failed=true;}
		check(failed && !new File(rejected.output+".trials.tsv.gz").exists(), "Missing counted window must fail before publishing search trials");
		NcrnaDdlLayoutCountsTest.run(root, first);
	}

	private static OptimizeNcrnaPartitionWeights.Config config(final String root, final String suffix){
		final OptimizeNcrnaPartitionWeights.Config c=new OptimizeNcrnaPartitionWeights.Config();
		c.left=root+"/compact_partitions_profiles512.tsv";
		c.right=root+"/compact_partitions_profiles1024.tsv";
		c.panels=root+"/panels.tsv";
		c.selection=root+"/selection.tsv";
		c.references=root+"/ref.fa";
		c.identities=root+"/one";
		c.ranks=root+"/compact_partitions.loci.tsv.gz";
		c.output=root+suffix;
		c.models=2;
		c.groups=3;
		c.randomDraws=200;
		return c;
	}
	private static byte[] dna(final String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static void check(final boolean condition, final String why){if(!condition){throw new AssertionError(why);}}
}
