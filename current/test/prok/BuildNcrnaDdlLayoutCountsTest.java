package test.prok;

import java.io.File;
import java.util.Arrays;
import java.util.ArrayList;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import kmer.AbstractKmerTable;
import map.LongHashSet;
import parse.LineParser1;

import prok.BuildNcrnaDdlLayoutCounts;
import prok.NcrnaDdlLayoutCounts;
import prok.NcrnaDdlSavedData;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.OptimizeNcrnaPartitionWeights;
import prok.RankDualNcrnaWindows;

/** Real producer/consumer interoperability on the original two-model witness.
 * No biological data, new identity calculation or alignment is introduced.
 * @author Brian Bushnell, Raiden
 */
final class BuildNcrnaDdlLayoutCountsTest{

	static void run(final RankDualNcrnaWindows.Config saved, final String control){
		partitionBits();
		gcStatistics();
		gcOrdering();
		gcReport(saved.output+"_gc_report");
		gcCostInference(saved.output+"_gc_cost");
		replicateUnits();
		replicateReport(saved.output+"_replicated_report");
		blacklistNamespace();
		final long c8=0x5555L, maskedC10=0x55555L&AbstractKmerTable.makeMiddleMask(10, false, 1);
		check(BuildNcrnaDdlLayoutCounts.prefixGc(new long[]{(1L<<32)|(0x7fffffffL-c8)}, 1, 8, false)==1
			&& BuildNcrnaDdlLayoutCounts.prefixGc(new long[]{(1L<<32)|(0x7fffffffL-maskedC10)}, 1, 10, true)==1,
			"Eight contributing short bases and nine contributing masked-ten bases must each give100percent GC for all-C words");
		mergeStudy(saved.output);
		final int[] padding=new int[2];
		BuildNcrnaDdlLayoutCounts.paddingBounds(100, 199, 0, 111, 160, padding);
		check(Arrays.equals(padding, new int[]{10,60}), "Forward flanks must exclude the one-based accepted span,including its endpoints");
		BuildNcrnaDdlLayoutCounts.paddingBounds(100, 199, 1, 111, 160, padding);
		check(Arrays.equals(padding, new int[]{40,90}), "Reverse-complemented saved windows must exchange and reverse the two physical flanks");
		BuildNcrnaDdlLayoutCounts.paddingBounds(100, 199, 1, 50, 300, padding);
		check(Arrays.equals(padding, new int[]{0,100}), "An accepted envelope covering the whole window leaves no background padding");
		final BuildNcrnaDdlLayoutCounts.Config c=config(saved, control, "_layout_counts");
		writeRecipes(c.recipes, null);
		final BuildNcrnaDdlLayoutCounts.Config gc=config(saved, control, "_gc_description");
		gc.recipes=c.recipes;
		gc.identities=saved.identities;
		gc.study="gcdescribe";
		BuildNcrnaDdlLayoutCounts.run(gc);
		check(rows(gc.output+".gc_models.tsv")==7 && rows(gc.output+".gc_family.tsv")==4,
			"The no-held-out description must retain both original model ordinals in all three sequence populations");
		final boolean vector=shared.Shared.SIMD;
		try{
			shared.Shared.SIMD=false;
			rebuiltFolds(saved, control);
		}finally{shared.Shared.SIMD=vector;}
		BuildNcrnaDdlLayoutCounts.main(new String[]{"recipes="+c.recipes,
			"panels="+c.panels, "selection="+c.selection, "ref="+c.references, "members="+c.members,
			"ledger="+c.ledger, "assignments="+c.assignments, "census="+c.census,
			"summary="+c.summary, "out="+c.output, "expectedmodels=2", "expectedgroups=3", "wordtag=leading", "t=1"});
		final NcrnaDdlLayoutCounts manifest=new NcrnaDdlLayoutCounts(c.output+".layouts.tsv", 2);
		check(manifest.schemaVersion==2 && manifest.layouts.length==ROWS.length && manifest.control(CONTROL)==1,
			"The producer must retain every declared layout and the explicitly named second-row legacy control");
		final NcrnaDdlSavedData data=NcrnaDdlSavedData.load(c.panels, c.selection, 2, 3);
		final int exponent=DynamicDemiLog.exponentBits();
		DynamicDemiLog.setExponent(4);
		try{
			legacyScores(c, control, manifest, data);
			leadingBitScores(c, manifest, data);
		}finally{DynamicDemiLog.setExponent(exponent);}
		long pooledChecksum=0;
		try(NcrnaDdlLayoutCounts.Reader reader=new NcrnaDdlLayoutCounts.Reader(c.output+".counts.tsv.gz", manifest, data)){
			int windows=-1;
			for(int layout=0; layout<ROWS.length; layout++){
				final NcrnaDdlLayoutCounts.Layout metadata=manifest.layouts[layout];
				check(metadata.terms==TERMS[layout] && metadata.referenceTerm==REFERENCES[layout]
					&& metadata.exponentBits==4 && metadata.mantissaBits==12 && metadata.overflowPolicy.equals("WRAP")
					&& metadata.memberKeyPayloadBytes==4L*metadata.memberBins && metadata.memberOriginPayloadBytes==0
					&& metadata.consensusKeyPayloadBytes==4L*metadata.consensusBins
					&& metadata.hashNamespace.equals(layout==1 ? "LEGACY_PARTITIONS" : "NATIVE_DDL")
					&& metadata.queryForm.equals("QUERY_DDL"),
					"Version-two metadata must retain declared plain-key terms, reference, full key width, wrap and query form");
				final char[][] counts=reader.next(true);
				if(windows<0){windows=counts.length;}
				check(counts.length==windows && windows>0, "Every layout must retain the identical original nonempty window roster");
				for(final char[] row:counts){
					check(row.length==8, "Both original models retain four interchange slots, with active terms first");
					for(int model=0; model<2; model++){
						for(int term=0; term<4; term++){
							check(term<metadata.terms ? row[4*model+term]<=LIMITS[layout][term] : row[4*model+term]==0,
								"Current counts must obey their explicit term geometry and leave every inactive slot zero");
							if(layout==3){pooledChecksum+=row[4*model+term];}
						}
					}
				}
			}
			reader.finish();
		}
		final BuildNcrnaDdlLayoutCounts.Config indexed=config(saved, control, "_layout_index");
		indexed.recipes=c.recipes;
		indexed.checkIndex="pooled256_long";
		indexed.indexSeconds=0;
		BuildNcrnaDdlLayoutCounts.run(indexed);
		final LineParser1 indexRow=new LineParser1('\t');
		try(Input input=new Input(indexed.output+".index.tsv")){
			check(Arrays.equals(input.reader.nextLine(), BuildNcrnaDdlLayoutCounts.INDEX_HEADER.getBytes(java.nio.charset.StandardCharsets.UTF_8)),
				"The index parity fixture must explicitly identify its measured population and timing contract");
			final byte[] row=input.reader.nextLine();
			check(row!=null, "The finite index fixture must publish its actual all-pair parity result");
			indexRow.set(row);
			check(indexRow.terms()==13 && indexRow.termEquals(indexed.checkIndex, 0) && indexRow.parseInt(1)==2
				&& indexRow.parseInt(2)==data.windowCount() && indexRow.parseLong(3)==2L*data.windowCount()
				&& indexRow.termEquals("PARITY_ONLY", 4) && indexRow.parseLong(5)==0
				&& indexRow.termEquals("NA", 6) && indexRow.termEquals("NA", 7)
				&& indexRow.parseLong(8)==pooledChecksum && indexRow.termEquals("FINITE_PARITY_NOT_TIMING", 12)
				&& input.reader.nextLine()==null,
				"Native object parity must preserve every saved score, including zeros, without labelling a finite fixture as timing");
		}
		check(!new File(indexed.output+".counts.tsv.gz").exists(), "Parity-only execution must not create an empty count cache");
		xorDefault(saved, control, c.recipes);
		final BuildNcrnaDdlLayoutCounts.Config seeded=config(saved, control, "_seed_noise");
		seeded.recipes=c.recipes;
		seeded.identities=saved.identities;
		seeded.study="seeds";
		final long previousSeed=CardinalityTracker.defaultSeed;
		BuildNcrnaDdlLayoutCounts.run(seeded);
		check(CardinalityTracker.defaultSeed==previousSeed, "The paired-seed study must restore process-wide seed state");
		try(Input input=new Input(seeded.output+".study.tsv")){
			check(input.reader.nextLine()!=null, "A paired-seed study requires a complete schema");
			final long[] seeds={12345,1,7,42,123,1337,20261114,8675309};
			int parents=-1;
			for(int row=0; row<16; row++){
				indexRow.set(input.reader.nextLine());
				if(row==0){parents=indexRow.parseInt(5);}
				check(indexRow.terms()==26 && indexRow.parseLong(25)==seeds[row%8]
					&& indexRow.parseInt(5)==parents
					&& indexRow.termEquals(row<8 ? "OLD_XOR_HIGH_BITS" : "LEADING_BIT", 0)
					&& indexRow.termEquals("NA", 22) && indexRow.termEquals("NA", 23),
					"Both tag forms must use each declared factory/default seed exactly once with absent-list GC preserved");
			}
			check(input.reader.nextLine()==null && rows(seeded.output+".seeds.summary.tsv")==13,
				"Six metrics per tag form summarize exactly sixteen paired observations");
			check(rows(seeded.output+".loci.tsv")==1+16*parents && rows(seeded.output+".locus_ranges.tsv")==1+2*parents
				&& rows(seeded.output+".tails.tsv")==3,
				"Every not-seen parent,including unavailable ones,retains all seed ranks and one exact tail decomposition per tag");
		}
		final BuildNcrnaDdlLayoutCounts.Config variance=config(saved, control, "_variance1024");
		variance.recipes=c.recipes;
		variance.identities=saved.identities;
		variance.study="variance1024";
		BuildNcrnaDdlLayoutCounts.run(variance);
		check(rows(variance.output+".study.tsv")==5 && rows(variance.output+".seeds.summary.tsv")==7,
			"The variance control is one XOR-only1024-bin case on the same four declared seeds");
		try(Input input=new Input(variance.output+".loci.tsv")){
			check(input.reader.nextLine()!=null, "Variance loci require explicit bins");
			indexRow.set(input.reader.nextLine());
			check(indexRow.termEquals("OLD_XOR_HIGH_BITS", 0) && indexRow.parseInt(1)==1024,
				"The explicit variance flag must reach the actual larger native DDL construction");
		}
		final BuildNcrnaDdlLayoutCounts.Config baseline=config(saved, control, "_leading_bit_baseline");
		final BuildNcrnaDdlLayoutCounts.Config partition=config(saved, control, "_partition_seeds");
		partition.recipes=c.recipes;
		partition.identities=saved.identities;
		partition.oldTags=true;
		partition.study="partitionseeds";
		BuildNcrnaDdlLayoutCounts.run(partition);
		check(rows(partition.output+".study.tsv")==17 && rows(partition.output+".seeds.summary.tsv")==13,
			"Both single-DDL selector partitions must retain all eight paired seeds and six ranking metrics");
		try(Input input=new Input(partition.output+".study.tsv")){
			check(input.reader.nextLine()!=null, "Partition seed rows need their explicit geometry and seed schema");
			for(int row=0; row<16; row++){
				indexRow.set(input.reader.nextLine());
				check(indexRow.terms()==26 && indexRow.termEquals(row<8 ? "SINGLE_DDL_256_256" : "SINGLE_DDL_128_384", 0)
					&& indexRow.termEquals("OLD_XOR_HIGH_BITS_RAW_PREHASH", 18),
					"Partitioning changes only native bin separation and preserves the mixed control's pre-hash length namespaces");
			}
		}
		baseline.recipes=c.recipes;
		baseline.identities=saved.identities;
		baseline.study="baseline";
		BuildNcrnaDdlLayoutCounts.run(baseline);
		try(Input input=new Input(baseline.output+".study.tsv")){
			check(Arrays.equals(input.reader.nextLine(), BuildNcrnaDdlLayoutCounts.STUDY_HEADER.getBytes(java.nio.charset.StandardCharsets.UTF_8)),
				"The fixed-recipe study must retain explicit populations,nonzero counts and experimental scope");
			indexRow.set(input.reader.nextLine());
			check(indexRow.terms()==25 && indexRow.termEquals("LEADING_BIT_BASELINE", 0)
				&& indexRow.parseInt(1)==0 && indexRow.parseInt(2)==0 && indexRow.parseInt(3)==2
				&& indexRow.parseInt(4)==data.windowCount() && indexRow.parseLong(14)>0
				&& Math.abs(indexRow.parseDouble(15)-indexRow.parseLong(14)/(double)data.windowCount())<0.00000000051
				&& indexRow.termEquals("NA", 22) && indexRow.termEquals("NA", 23)
				&& input.reader.nextLine()==null,
				"Baseline must score only the released512-bin recipe and count every saved window exactly once");
		}
		final OptimizeNcrnaPartitionWeights.Config fit=new OptimizeNcrnaPartitionWeights.Config();
		fit.panels=saved.panels;
		fit.selection=saved.selection;
		fit.references=saved.references;
		fit.identities=saved.identities;
		fit.ranks=saved.output+".loci.tsv.gz";
		fit.layouts=c.output+".layouts.tsv";
		fit.counts=c.output+".counts.tsv.gz";
		fit.control=CONTROL;
		fit.output=c.output+"_fit";
		fit.models=2;
		fit.groups=3;
		fit.randomDraws=7;
		OptimizeNcrnaPartitionWeights.run(fit);
		check(rows(fit.output+".decision.tsv")==ROWS.length+1,
			"The fitter must consume native producer counts and report every declared DDL-query layout without raw-rank validation");
		boolean rejected=false;
		try{BuildNcrnaDdlLayoutCounts.run(c);}
		catch(AssertionError expected){
			rejected=expected.getMessage()!=null && expected.getMessage().contains("exists and overwrite=false");
		}
		catch(IllegalArgumentException expected){
			rejected=expected.getMessage()!=null && expected.getMessage().contains("output fresh");
		}
		check(rejected, "Producer outputs cannot overwrite an existing successful count population");
		for(final String fault:new String[]{"duplicate", "missing_control", "bad_control", "bad_query_control", "bad_format", "bad_part"}){
			final BuildNcrnaDdlLayoutCounts.Config bad=config(saved, control, "_layout_reject_"+fault);
			writeRecipes(bad.recipes, fault);
			rejected=false;
			try{BuildNcrnaDdlLayoutCounts.run(bad);}catch(IllegalArgumentException expected){rejected=true;}
			check(rejected && !new File(bad.output+".layouts.tsv").exists()
				&& !new File(bad.output+".counts.tsv.gz").exists(),
				"Bad recipe must fail before publishing counts or metadata: "+fault);
		}
		System.err.println("NCRNA_DDL_LAYOUT_COUNT_PRODUCER_TEST_PASS");
	}

	private static void gcStatistics(){
		final double[] stats=NcrnaDdlSavedData.gcStatistics(1)[0];
		check(Double.isNaN(NcrnaDdlSavedData.gcSd(stats)), "An unobserved model must report missing SD rather than invented zero variance");
		NcrnaDdlSavedData.addGc(stats, 0);
		check(NcrnaDdlSavedData.gcSd(stats)==0, "A singleton has declared sample SD zero");
		NcrnaDdlSavedData.addGc(stats, 0.5);
		NcrnaDdlSavedData.addGc(stats, 1);
		NcrnaDdlSavedData.addGc(stats, Double.NaN);
		check(stats[0]==3 && stats[1]==0.5 && NcrnaDdlSavedData.gcSd(stats)==0.5
			&& stats[3]==0 && stats[4]==1 && stats[5]==1,
			"GC moments must preserve sample SD, extrema and a separate missing-observation count");
		check(shared.Tools.calcGC("NACGTN".getBytes(java.nio.charset.StandardCharsets.US_ASCII))==0.5f
			&& Float.isNaN(shared.Tools.calcGC("NNNN".getBytes(java.nio.charset.StandardCharsets.US_ASCII))),
			"The existing helper excludes ambiguous bases and makes all-ambiguous GC missing");
		System.err.println("NCRNA_DDL_GC_MOMENTS_TEST_PASS");
	}

	/** Vary selector/high hash bits while holding the routing source fixed.
	 * This checks the routing mask,not unconditional independence of the wrapped
	 * native exponent. Bit9 is a real selector at1024bins,so the old mask fails. */
	private static void partitionBits(){
		for(final int bins:new int[]{512,1024}){
			final int lowBits=Integer.numberOfTrailingZeros(bins);
			for(int i=0; i<1000; i++){
				final long hash=shared.Tools.hash64shift(i);
				final int selected=BuildNcrnaDdlLayoutCounts.partitionFor(hash, bins, 2, false);
				check(selected>=1 && selected<=3, "Long words must use precisely the three reserved quarters");
				for(int bit=0; bit<64; bit++){
					if(bit>=lowBits && bit<36){continue;}
					check(BuildNcrnaDdlLayoutCounts.partitionFor(hash^(1L<<bit), bins, 2, false)==selected,
						"Changing actual bucket-selector bits or hash bits36..63 must not affect the long-word quarter");
				}
				check(BuildNcrnaDdlLayoutCounts.partitionFor(hash, bins, 2, true)==0
					&& BuildNcrnaDdlLayoutCounts.partitionFor(hash, bins, 1, false)==1,
					"Short words own the first partition and long words own the second in the equal split");
				if(bins==512){
					check(selected==1+(int)((((hash>>>9)&0x7ffffffL)*3)>>>27),
						"The generalized routing must preserve every accepted512-bin partition result");
				}
			}
		}
	}

	private static void replicateUnits(){
		percentiles();
		excludedMemberPresence();
		final structures.IntHashMap held=new structures.IntHashMap();
		held.put(11, 1);
		final prok.NcrnaExactSequenceSet exact=new prok.NcrnaExactSequenceSet();
		exact.add("AACCGGTT".getBytes(java.nio.charset.StandardCharsets.US_ASCII));
		final byte[][] members={"TTTTAAAA".getBytes(java.nio.charset.StandardCharsets.US_ASCII),
			"aaccgguu".getBytes(java.nio.charset.StandardCharsets.US_ASCII),
			"AACCGGTA".getBytes(java.nio.charset.StandardCharsets.US_ASCII)};
		final boolean[] excluded=prok.NcrnaDdlFolds.exclude(new int[]{11,22,33}, members, held, exact);
		check(Arrays.equals(excluded, new boolean[]{true,true,false}),
			"A whole held-out assembly and its exact normalized duplicate in another assembly must be removed;one changed base must survive");
		check(BuildNcrnaDdlLayoutCounts.consensusBonus(0, 0, 0)==0
			&& BuildNcrnaDdlLayoutCounts.consensusBonus(0, 0.1, 0)==1
			&& BuildNcrnaDdlLayoutCounts.consensusBonus(0, 0.3, 7)==3,
			"Empty fit models distinguish fixed-zero from scaled minimum-one,and fractions use remaining voter counts with ceiling");
		final map.IntHashMap votes=new map.IntHashMap(512);
		final DynamicDemiLog empty=BuildNcrnaDdlLayoutCounts.chooseReference(votes, 512, 0, 9, true);
		for(final char key:empty.maxArray()){check(key==0, "No members plus fixed-zero must preserve a genuine empty native model sketch");}
		votes.increment((13<<16)|1234, BuildNcrnaDdlLayoutCounts.consensusBonus(0, 0.1, 0));
		final DynamicDemiLog consensus=BuildNcrnaDdlLayoutCounts.chooseReference(votes, 512, 0, 9, true);
		check(consensus.maxArray()[13]==1234, "The scaled minimum retains an actual nonzero consensus winner when no members remain");
		final double[] stats=prok.SummarizeNcrnaDdlReplicates.moments(new double[][]{{0,2},{10,12}});
		check(stats[0]==6 && Math.abs(stats[1]-Math.sqrt(104.0/3))<1e-12
			&& stats[2]==2 && stats[3]==50 && stats[4]==0 && Math.abs(stats[5]-Math.sqrt(26))<1e-12,
			"Crossed uncertainty must separate literal additive seed/split effects,not count the four cells as independent");
		final double[] fixed=prok.SummarizeNcrnaDdlReplicates.moments(new double[][]{{0,2}});
		check(fixed[0]==1 && fixed[2]==2 && fixed[3]==0 && fixed[5]==1,
			"The fixed-panel uncertainty has only hashseed variation and its ordinary sample standard error");
		final double[] interaction=prok.SummarizeNcrnaDdlReplicates.moments(new double[][]{{0,2},{2,0}});
		check(interaction[2]==0 && interaction[3]==0 && interaction[4]==4 && interaction[5]==1,
			"A pure seed-by-split interaction must retain its variance when estimated main components clip to zero");
		final double[] difference=prok.SummarizeNcrnaDdlReplicates.moments(new double[][]{{2,4},{12,14}});
		check(difference[0]==8 && Math.abs(difference[5]-Math.sqrt(26))<1e-12,
			"Nonzero paired differences retain their own crossed uncertainty rather than independent-method error propagation");
		System.err.println("NCRNA_DDL_EXCLUSION_AND_CROSSED_MOMENTS_TEST_PASS");
	}

	/** Uses the same presence counter called by the actual definition-0 path.
	 * A wins the unexcluded frequency tie;excluding its assembly and duplicates
	 * must change the selected prefix to the retained C word. */
	private static void excludedMemberPresence(){
		final byte[] a="AAAAAAAAAAAA".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
		final byte[] c="CCCCCCCCCCCC".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
		final byte[][][] models={{a,a,c},{a,c,c}};
		final structures.IntHashMap held=new structures.IntHashMap();
		held.put(11, 1);
		final prok.NcrnaExactSequenceSet exact=new prok.NcrnaExactSequenceSet();
		exact.add(a);
		final boolean[][] excluded={prok.NcrnaDdlFolds.exclude(new int[]{11,22,33}, models[0], held, exact),
			prok.NcrnaDdlFolds.exclude(new int[]{44,55,66}, models[1], held, exact)};
		check(Arrays.equals(excluded[0], new boolean[]{true,true,false}) && Arrays.equals(excluded[1], new boolean[]{true,false,false}),
			"The actual fold exclusion must remove held-out A and its cross-assembly exact duplicates while preserving both models' C members");
		for(final int k:new int[]{8,10,11}){
			final long[] before=BuildNcrnaDdlLayoutCounts.memberPresenceFrequencies(k, k!=8, models, null);
			final long[] after=BuildNcrnaDdlLayoutCounts.memberPresenceFrequencies(k, k!=8, models, excluded);
			final long expectedC=(((1L<<(2*k))-1)/3)&AbstractKmerTable.makeMiddleMask(k, false, k==8 ? 0 : 1);
			check(before.length==2 && (before[0]>>>32)==2 && (before[1]>>>32)==2
				&& 0x7fffffffL-(before[before.length-1]&0xffffffffL)==0,
				"Repeated words and duplicate members count once per model;the smaller raw A word wins the unexcluded tie");
			check(after.length==1 && (after[0]>>>32)==2 && 0x7fffffffL-(after[0]&0xffffffffL)==expectedC,
				"Excluded-member words must be absent from the actual positive-member prefix;the retained C word becomes the sole winner");
		}
		System.err.println("NCRNA_DDL_EXCLUDED_MEMBER_PRESENCE_PREFIX_TEST_PASS");
	}

	/** Runs the actual group join,assignment,exclusion and native reference builder.
	 * The tiny cache has two comparable assemblies,so it explicitly uses two folds
	 * and zero padding prefixes. Production reporting rejects that witness design. */
	private static void rebuiltFolds(final RankDualNcrnaWindows.Config saved, final String control){
		final BuildNcrnaDdlLayoutCounts.Config c=config(saved, control, "_rebuilt_folds");
		c.study="replicates";
		c.identities=saved.identities;
		c.oldTags=true;
		c.folds=2;
		c.longPrefix=0;
		Arrays.fill(c.shortPrefixes, 0);
		c.indexReference=true;
		c.sourceGroups=c.output+".source_groups.tsv";
		c.genomes=c.output+".genomes.tsv";
		try(Output groups=new Output(c.sourceGroups); Output genomes=new Output(c.genomes); Output recipes=new Output(c.recipes)){
			groups.writer.println(prok.NcrnaConstructionMembers.GROUP_HEADER);
			for(int i=0; i<3; i++){
				final String path=i==0 ? "/fixture/a.fa" : i==1 ? "/fixture/empty.fa" : "/fixture/b.fa";
				groups.writer.println("evaluation\t"+i+"\t"+path+"\t"+path+"\t.\t.\t.\t.\t.\t.\t"+(10+i)+"\t.\tfalse\tRESOLVED");
			}
			genomes.writer.println(prok.NcrnaConstructionMembers.GENOME_HEADER);
			genomes.writer.println("0\ttid_1_species\tfixture\t1\tBacteria\tfixture\t/fixture/a.gff\t/fixture/a.fa\t2\t1000\t1");
			recipes.writer.println(BuildNcrnaDdlLayoutCounts.RECIPE_HEADER);
			recipes.writer.println(ROWS[1]);
			for(int pair=0; pair<3; pair++){
				final int shortK=pair==2 ? 8 : 7, longK=pair==0 ? 9 : 10;
				recipes.writer.println("native_s"+shortK+"_l"+longK+"_mixed512_pooled_both_w100\t"+shortK+"\t"+longK
					+"\tMIDDLE_RIGHT\tMIXED\t512\t0\tPOOL_BOTH\t0\t100\t4\t12\tNATIVE_DDL\tfalse\tQUERY_DDL");
			}
		}
		BuildNcrnaDdlLayoutCounts.run(c);
		int modelZero=0, otherVoters=0;
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(c.output+".models.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(p.parseInt(1)==0){
					check(p.parseInt(2)==0, "The duplicate A gene in the other assembly must remove modelzero's voter in BOTH folds");
					modelZero++;
				}
				else{otherVoters+=p.parseInt(2);}
			}
		}
		check(modelZero==2 && otherVoters==1, "Whole-assembly removal and cross-assembly duplicate removal must both reach actual rebuilt models");
		int contributorFold=-1;
		try(Input in=new Input(c.output+".folds.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(p.parseInt(1)==10){contributorFold=p.parseInt(2);}
			}
		}
		check(contributorFold>=0, "The original donor assembly must have an actual saved fold assignment");
		boolean empty=false, scaled=false;
		try(Input in=new Input(c.output+".replicates.tsv")){
			in.reader.nextLine();
			int rows=0;
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==28, "Rebuilt rows retain their explicit original-assembly-contribution label");
				if(p.termEquals("ALL", 27) && p.termEquals("bonus_fixed_0", 3)){check(p.parseDouble(21)==0, "An excluded model with zero consensus bonus must retain no consensus keys"); empty=true;}
				if(p.termEquals("ALL", 27) && p.termEquals("bonus_tenths_1", 3)){check(p.parseDouble(21)>0, "Scaled minimum-one must actually insert native consensus keys into the empty-fit model"); scaled=true;}
				if(!p.termEquals("ALL", 27)){
					check(p.parseInt(12)==(p.parseInt(2)==contributorFold ? p.termEquals("CONTRIBUTOR", 27) ? 1 : 0
						: p.termEquals("NONCONTRIBUTOR", 27) ? 1 : 0),
						"Contribution uses the original donor group even after its voters were removed; empty subsets remain explicit");
					check(p.termEquals("NA", 21), "Profile-level retention is not a locus-subset measurement");
				}
				if(p.termEquals("CONSENSUS_INDEX", 3)){check(p.termEquals("NA", 21) && p.termEquals("NA", 22) && p.termEquals("NA", 23), "Absent index-reference retention must stay literal NA through native output");}
				rows++;
			}
			check(rows==300 && empty && scaled, "The two-fold native fixture emits all49methods and one reference in all three contribution populations");
		}
		final String roster=c.output+".reject.roster.tsv";
		final String rejectedTable=c.output+".reject.table.tsv";
		try(Input in=new Input(c.output+".replicates.tsv"); Output out=new Output(rejectedTable)){
			out.writer.println(in.reader.nextLine());
			boolean found=false;
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(!p.termEquals("s7_m9_mixed_512_padding", 3)){continue;}
				out.writer.println(row);
				found=true;
				break;
			}
			check(found, "The negative report fixture must use an actual emitted nonproduction padding row");
		}
		try(Output out=new Output(roster)){
			out.writer.println("assignment_seed\tddl_seed\ttable");
			out.writer.println(c.foldSeed+"\t"+c.seed+"\t"+rejectedTable);
		}
		boolean rejected=false;
		try{
			prok.SummarizeNcrnaDdlReplicates.main(new String[]{"in="+roster, "out="+c.output+".reject",
				"seeds="+c.seed+",7", "foldseeds="+c.foldSeed+",1", "t=1"});
		}catch(IllegalArgumentException expected){
			rejected=expected.getMessage().contains("Method identity must agree")
				&& expected.getMessage().contains("s7_m9_mixed_512_padding");
		}
		check(rejected && !new File(c.output+".reject.summary.tsv").exists(),
			"The finite two-fold/empty-prefix witness must never be accepted as the released five-fold/default-prefix study");
		System.err.println("NCRNA_DDL_REBUILT_FOLDS_END_TO_END_TEST_PASS");
		c.output+="_followup";
		c.study="replicatefollowup";
		BuildNcrnaDdlLayoutCounts.run(c);
		check(rows(c.output+".replicates.tsv")==121,
			"The actual follow-up runs only twenty settings across two toy folds and three contribution populations");
		try(Input in=new Input(c.output+".replicates.tsv")){
			in.reader.nextLine();
			boolean masked11=false, controlSeen=false;
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==30, "Follow-up rows retain exact rank histograms beside the unchanged aggregate schema");
				masked11|=p.parseInt(5)==11;
				controlSeen|=p.parseInt(10)==100 && p.parseInt(5)==10;
				final long[] histogram=new long[2], metrics=new long[4];
				p.setBounds(28);
				prok.SummarizeNcrnaDdlReplicates.addHistogram(row, p.a(), p.b(), histogram, metrics);
				check(metrics[0]==p.parseLong(12) && metrics[1]==p.parseLong(14),
					"The native scorer's saved distribution must conserve its actual locus count and rank sum");
			}
			check(masked11 && controlSeen, "Both the new masked11 iterator and explicitly repeated100controls must execute");
		}
		c.output+="_extended";
		c.study="replicatesizesextended";
		boolean missingBackground=false;
		try{BuildNcrnaDdlLayoutCounts.run(c);}
		catch(IllegalArgumentException expected){missingBackground=expected.getMessage().contains("Background blacklist cannot be learned from an empty");}
		check(missingBackground,
			"The22-cell fixed-prefix worker must reject this toy's absent padding background instead of silently substituting empty lists");
		c.output+="_index";
		c.study="replicateindex";
		BuildNcrnaDdlLayoutCounts.run(c);
		try(Input in=new Input(c.output+".index.tsv")){
			in.reader.nextLine();
			p.set(in.reader.nextLine());
			check(p.parseInt(1)==2 && p.parseInt(3)>=2 && p.parseInt(5)>0 && p.parseInt(5)<=2
				&& p.parseInt(7)==0 && in.reader.nextLine()==null,
				"The index-only path restores both toy loci with valid saved ranks and exactly assignment-invariant percentiles");
		}
		check(!new File(c.output+".replicates.tsv").exists(), "Index extraction must bypass the sketch-scoring output entirely");
		System.err.println("NCRNA_DDL_SAVED_INDEX_PERCENTILES_TEST_PASS");
		final byte[] padding=new byte[100000], alphabet={'A','C','G','T'};
		final java.util.Random random=new java.util.Random(7);
		for(int i=0; i<padding.length; i++){padding[i]=alphabet[random.nextInt(4)];}
		final RankSavedNcrnaWindowsTest.Fixture padded=new RankSavedNcrnaWindowsTest.Fixture(
			new File(c.output).getParent(), "padded_extension", padding);
		c.panels=padded.dir+"/panels.tsv";
		c.selection=padded.dir+"/selection.tsv";
		c.output=padded.dir+"/extension";
		c.study="replicatesizesextended";
		BuildNcrnaDdlLayoutCounts.run(c);
		check(rows(c.output+".replicates.tsv")==133,
			"The positive native fixture learns actual fit-side padding lists through16384short words and emits all22settings in both folds and all three populations");
		System.err.println("NCRNA_DDL_PADDED_EXTENSION_NATIVE_TEST_PASS");
		c.output=padded.dir+"/gc_ordering";
		c.study="replicategc";
		c.shortPrefixes[2]=256;
		c.longPrefix=2048;
		BuildNcrnaDdlLayoutCounts.run(c);
		check(rows(c.output+".replicates.tsv")==289 && rows(c.output+".gc_fit.tsv")==9,
			"GC production path runs48settings in both toy folds and three populations,with one fit-stat row per model/fold/source");
		gcIsolation(c, padding);
		try(Input in=new Input(c.output+".replicates.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==39 && p.parseInt(8)==256 && p.parseInt(9)==2048
					&& p.parseInt(10)==(p.parseInt(7)==512 ? 60 : 100),
					"Actual GC comparisons preserve the chosen padding sizes and each bin count's decided consensus bonus");
				final long[] histogram=new long[3], metrics=new long[4];
				p.setBounds(28);
				prok.SummarizeNcrnaDdlReplicates.addHistogram(row, p.a(), p.b(), histogram, metrics);
				check(metrics[0]==p.parseLong(12) && metrics[1]==p.parseLong(14),
					"GC native distributions include any all-filtered sentinel without dropping its locus");
			}
		}
		System.err.println("NCRNA_DDL_GC_NATIVE_ORDERING_TEST_PASS");
		final String gcTable=c.output+".replicates.tsv";
		c.gcFits=c.output+".workers.tsv";
		try(Output out=new Output(c.gcFits)){
			out.writer.println("assignment_seed\tddl_seed\ttable");
			out.writer.println(c.foldSeed+"\t12345\t"+gcTable);
		}
		c.output+="_losses";
		c.study="gclosses";
		BuildNcrnaDdlLayoutCounts.run(c);
		check(rows(c.output+".gc_curve.tsv")==101 && rows(c.output+".gc_choices.tsv")==21,
			"GC-only diagnostics retain both sources,25widths,two count policies,two fit folds and three descriptive choices");
		try(Input in=new Input(c.output+".gc_curve.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(p.parseInt(2)==30){
					for(int field=4; field<=8; field++){check(p.parseLong(field)==0,
						"Fewer than30own-source observations disables the gate,even for window fallback with many members");}
				}
			}
		}
		System.err.println("NCRNA_DDL_GC_LOSS_CURVES_NATIVE_TEST_PASS");
		c.output+="_all_data";
		c.study="gcalllosses";
		BuildNcrnaDdlLayoutCounts.run(c);
		check(rows(c.output+".gc_all_fit.tsv")==5 && rows(c.output+".gc_all_limits.tsv")==3
			&& rows(c.output+".gc_all_losses.tsv")==63,
			"All-data GC fits each model/source once and measures the complete SD/margin grid without ranking folds");
		boolean completeWindows=false, infiniteOnly=false, rangeLossless=false;
		try(Input in=new Input(c.output+".gc_all_fit.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				check(p.parseInt(0)==-1 && p.parseInt(4)>0, "All-data GC retains every original member instead of excluding judging assemblies or duplicates");
				if(p.parseInt(1)==0 && p.termEquals("window", 2)){
					check(p.parseInt(3)==10, "All-data window GC combines all ten real alpha/beta observations rather than fitting either half");
					completeWindows=true;
				}
			}
		}
		try(Input in=new Input(c.output+".gc_all_limits.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(p.parseLong(1)>0 && p.parseLong(1)==p.parseLong(4)){
					check(p.termEquals("NA", 5) && p.termEquals("INF", 6) && p.termEquals("INF", 7),
						"A real all-infinite sample has no finite maximum and no finite lossless SD width");
					infiniteOnly=true;
				}
			}
		}
		try(Input in=new Input(c.output+".gc_all_losses.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(p.termEquals("window", 0) && p.termEquals("MARGIN", 1) && p.parseDouble(2)==0){
					check(p.parseLong(5)==0 && p.parseLong(7)==0 && p.parseLong(9)==0,
						"All accepted all-data windows lie within their model's observed window GC range");
					rangeLossless=true;
				}
			}
		}
		check(completeWindows && infiniteOnly && rangeLossless, "The all-data native fixture must actually inspect complete statistics,infinite-only distances and lossless ranges");
		System.err.println("NCRNA_DDL_ALL_DATA_GC_LIMITS_NO_HOLDOUTS_TEST_PASS");
		c.gcAllFit=c.output+".gc_all_fit.tsv";
		c.gcAllLimits=c.output+".gc_all_limits.tsv";
		c.output+="_ranking_reject";
		c.study="replicategcall";
		boolean unavailable=false;
		try{BuildNcrnaDdlLayoutCounts.run(c);}
		catch(IllegalArgumentException expected){unavailable=expected.getMessage().contains("measured finite all-data lossless");}
		check(unavailable && !new File(c.output+".replicates.tsv").exists(),
			"All-data ranking must load the shared table and reject unavailable lossless limits before writing rows,never refit the current judging folds");
		System.err.println("NCRNA_DDL_SHARED_ALL_DATA_GC_LOAD_NO_REFIT_TEST_PASS");
		c.output+="_power";
		c.study="gcallpower";
		BuildNcrnaDdlLayoutCounts.run(c);
		check(rows(c.output+".gc_power.tsv")==73 && rows(c.output+".gc_power_histogram.tsv")==217
			&& rows(c.output+".gc_model_limits.tsv")==3 && rows(c.output+".gc_model_quantiles.tsv")==11,
			"Power runs the full36-condition roster on both explicit window populations without bins,seeds or GC refitting");
		boolean ownCount=false, noWindows=false, powerRange=false;
		try(Input in=new Input(c.output+".gc_model_limits.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(p.parseInt(0)==0){ownCount=p.parseLong(2)==10 && p.parseLong(3)==10 && p.termEquals("10_TO29", 7);}
				else{noWindows=p.parseLong(2)==0 && p.termEquals("NA", 5) && p.termEquals("NA", 6) && p.termEquals("NO_WINDOWS", 8);}
			}
		}
		try(Input in=new Input(c.output+".gc_power.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(p.termEquals("MARGIN", 1) && p.parseDouble(2)==0){
					check(p.parseLong(8)==0 && p.parseLong(10)==0 && p.parseLong(12)==0,
						"Observed window ranges preserve every accepted pair,window and locus even when other models are filtered");
					powerRange=true;
				}
			}
		}
		check(ownCount && noWindows && powerRange, "The native power fixture must inspect ten accepted observations and the separate no-window fallback model");
		int pairs=0;
		try(Input in=new Input(c.output+".gc_distances.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==24 && p.parseInt(0)==pairs && p.parseInt(11)==0 && p.parseLong(18)==10,
					"The primary ledger emits each of the ten accepted model-zero pairs once,in ordinal order,independently of parent links");
				pairs++;
			}
		}
		check(pairs==10, "Separate locus links must never multiply the primary accepted-pair census");
		c.output+="_additions";
		c.study="replicategcall";
		c.gcAdditions=true;
		BuildNcrnaDdlLayoutCounts.run(c);
		check(rows(c.output+".replicates.tsv")==73,
			"Only12new GC alternatives execute over the two actual toy folds and three contribution populations; controls are not rescored");
		try(Input in=new Input(c.output+".replicates.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				boolean named=false;
				for(int setting=0; setting<BuildNcrnaDdlLayoutCounts.GC_EXTRA_SETTINGS.length; setting++){
					named|=p.termEquals(BuildNcrnaDdlLayoutCounts.gcExtraMethod(p.parseInt(7), setting), 3);
				}
				check(named && p.terms()==39, "Every new native row belongs to the disjoint extension roster with complete failure histograms");
			}
		}
		System.err.println("NCRNA_DDL_GC_POWER_AND_NEW_ONLY_NATIVE_TEST_PASS");
		literalGcPower(saved, control, c.recipes);
	}

	/** A fixed query subset loads literal full-data statistics:40zero-GC and
	 * eight unit-GC accepted windows. Its first locus has G,A,G... windows;
	 * only its second window survives1SD. The SEEN beta locus is excluded. */
	private static void literalGcPower(final RankDualNcrnaWindows.Config saved, final String control, final String recipes){
		final RankSavedNcrnaWindowsTest.Fixture fixture=new RankSavedNcrnaWindowsTest.Fixture(
			new File(saved.output).getParent(), "gc_power_literal");
		final BuildNcrnaDdlLayoutCounts.Config c=config(saved, control, "_gc_power_literal");
		c.recipes=recipes;
		c.panels=fixture.dir+"/panels.tsv";
		c.selection=fixture.dir+"/selection.tsv";
		c.identities=saved.identities;
		c.output=fixture.dir+"/power";
		c.gcAllFit=fixture.dir+"/shared.gc_all_fit.tsv";
		c.gcAllLimits=fixture.dir+"/shared.gc_all_limits.tsv";
		c.study="gcallpower";
		c.fixedCheck=true;
		final double mean=1.0/6, sd=Math.sqrt(20.0/141), fallback=1+mean;
		final double maximum=(1-mean)/sd;
		double inclusive=maximum;
		while(!BuildNcrnaDdlLayoutCounts.gcWithin(1, new double[]{mean,sd,0,1}, inclusive, false)){inclusive=Math.nextUp(inclusive);}
		try(Output fit=new Output(c.gcAllFit); Output limits=new Output(c.gcAllLimits)){
			fit.writer.println(BuildNcrnaDdlLayoutCounts.GC_FIT_HEADER);
			fit.writer.println("-1\t0\twindow\t48\t1\t0\t0\t"+mean+'\t'+sd+"\t0\t1\tfalse\tfalse\t"+mean+"\t1");
			fit.writer.println("-1\t0\tmember\t48\t1\t0\t0\t0\t0\t0\t0\tfalse\tfalse\t"+mean+"\t1");
			fit.writer.println("-1\t1\twindow\t0\t1\t0\t0\t"+fallback+"\t0\t"+fallback+'\t'+fallback+"\ttrue\tfalse\t"+mean+"\t1");
			fit.writer.println("-1\t1\tmember\t0\t1\t0\t0\t1\t0\t1\t1\tfalse\tfalse\t"+mean+"\t1");
			limits.writer.println(BuildNcrnaDdlLayoutCounts.GC_ALL_LIMIT_HEADER);
			limits.writer.println("window\t48\t48\t0\t0\t"+maximum+'\t'+maximum+'\t'+inclusive+"\t0\tLITERAL_FULL_DATA_40_ZERO_8_ONE");
			limits.writer.println("member\t48\t40\t0\t8\t0\tINF\tINF\t1\tLITERAL_FULL_DATA_40_ZERO_8_ONE");
		}
		BuildNcrnaDdlLayoutCounts.run(c);
		final LineParser1 p=new LineParser1('\t');
		int populations=0, histogram=0, loci=0;
		try(Input in=new Input(c.output+".gc_power.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(!p.termEquals("SD", 1) || p.parseDouble(2)!=1){continue;}
				final boolean ranking=p.termEquals("FIXED_NOT_SEEN", 0);
				check(p.parseLong(4)==(ranking ? 7 : 11) && Math.abs(p.parseDouble(5)-(ranking ? 13.0/7 : 19.0/11))<1e-12
					&& p.parseDouble(6)==2 && p.parseLong(8)==(ranking ? 6 : 8) && p.parseLong(11)==1 && p.parseLong(12)==0,
					"Literal power distinguishes11saved/7ranking windows,nonzero removal means/median,and a retained fixed locus despite six rejected windows");
				check(p.parseLong(13)==(ranking ? 14 : 22) && p.parseLong(14)==(ranking ? 1 : 3)
					&& p.parseDouble(16)==0 && p.parseLong(18)==0,
					"Retained candidate totals are baseline minus removed;the rescued locus contributes no failed-search candidate cost");
				populations++;
			}
		}
		try(Input in=new Input(c.output+".gc_power_histogram.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(!p.termEquals("SD", 1) || p.parseDouble(2)!=1){continue;}
				final boolean ranking=p.termEquals("FIXED_NOT_SEEN", 0);
				final int removed=p.parseInt(3);
				check(p.parseLong(4)==(removed==0 ? 0 : removed==1 ? ranking ? 1 : 3 : ranking ? 6 : 8),
					"The actual0..2 removal histogram retains zero counts and exactly three/eight versus one/six windows");
				histogram++;
			}
		}
		try(Input in=new Input(c.output+".gc_survival.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(!p.termEquals("SD", 1) || p.parseDouble(2)!=1){continue;}
				check(p.termEquals("input_a", 5) && p.termEquals("p0", 6) && p.termEquals("true", 7),
					"The fixed locus survives through its second A window after the first G window fails;SEEN beta is not a fixed locus");
				loci++;
			}
		}
		check(populations==2 && histogram==6 && loci==1 && !new File(c.output+".gc_model_limits.tsv").exists(),
			"The literal witness must inspect both emitted populations,complete histograms and the original fixed locus without recomputing full-data limits");
		System.err.println("NCRNA_DDL_GC_POWER_LITERAL_HISTOGRAM_AND_SECOND_WINDOW_RESCUE_TEST_PASS");
	}

	/** Two loci and two candidates: a failed filtered locus costs one retained
	 * candidate,not the original failure sentinel of three. */
	private static void gcCostInference(final String prefix){
		final String summary=prefix+".summary.tsv", power=prefix+".power.tsv";
		try(Output out=new Output(summary)){
			out.writer.println("cohort\tmethod\tmetric\tassignments\thash_seeds\tmean\tstandard_deviation\thash_variance\tsplit_variance\tinteraction_variance\tmean_standard_error\tparents\tlow_parents\tassembly_contribution");
			for(final int bins:new int[]{512,1024}){
				for(int setting=0; setting<BuildNcrnaDdlLayoutCounts.GC_ALL_SETTINGS.length; setting++){
					final String method=BuildNcrnaDdlLayoutCounts.gcAllMethod(bins, setting);
					final boolean failed=BuildNcrnaDdlLayoutCounts.GC_ALL_SETTINGS[setting].equals("filter_z4");
					out.writer.println("REBUILT_FIVE_FOLD\t"+method+"\trank_sum\t2\t2\t"+(failed ? 4 : 3)+"\t0\t0\t0\t0\t0\t2\t1\tALL");
					out.writer.println("REBUILT_FIVE_FOLD\t"+method+"\tparents_all_acceptable_filtered\t2\t2\t"+(failed ? 1 : 0)+"\t0\t0\t0\t0\t0\t2\t1\tALL");
				}
			}
		}
		try(Output out=new Output(power)){
			out.writer.println(BuildNcrnaDdlLayoutCounts.GC_POWER_HEADER);
			out.writer.println("REBUILT_FIVE_FOLD\tNONE\t0\t2\t2\t0\t0\t2\t0\t2\t0\t2\t0\t4\t4\t2\t2\t0\t0\tLITERAL");
			out.writer.println("REBUILT_FIVE_FOLD\tSD\t4\t2\t2\t1\t1\t2\t1\t2\t1\t2\t1\t4\t2\t1\t1\t50\t1\tLITERAL");
			out.writer.println("REBUILT_FIVE_FOLD\tSD\t4.25\t2\t2\t1\t1\t2\t1\t2\t1\t2\t1\t4\t2\t1\t1\t50\t1\tLITERAL");
		}
		prok.SummarizeNcrnaDdlReplicates.writeGcCosts(summary, power, prefix);
		final LineParser1 p=new LineParser1('\t');
		int baseline=0, filtered=0, unmeasured=0;
		try(Input in=new Input(prefix+".gc_cost.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==20, "Every cost estimate carries its explicit candidate and locus denominators");
				if(p.termEquals("NONE", 2)){
					check(p.parseDouble(13)==3 && p.parseDouble(16)==0 && p.parseDouble(17)==0,
						"No filtering preserves the measured three-attempt baseline and saves zero");
					baseline++;
				}else if(p.parseDouble(3)==4){
					check(p.parseDouble(13)==2 && p.parseDouble(16)==1 && Math.abs(p.parseDouble(17)-100.0/3)<1e-8
						&& p.parseLong(18)==1 && p.parseDouble(6)==1 && p.parseDouble(7)==1,
						"One surviving rank plus one failed retained candidate costs two attempts,saving one of three;the sentinel is not an attempt count");
					filtered++;
				}else{
					check(p.termEquals("false", 12) && p.termEquals("NA", 13) && p.termEquals("NA", 16),
						"An unranked SD threshold retains measured candidate counts but no invented alignment estimate");
					unmeasured++;
				}
			}
		}
		check(baseline==2 && filtered==2 && unmeasured==2, "Both original bin widths must show baseline,filtered and unmeasured cases");
		boolean rejected=false;
		try{prok.SummarizeNcrnaDdlReplicates.writeGcCosts(summary, power, prefix);}
		catch(IllegalArgumentException | AssertionError expected){
			rejected=expected.getMessage().contains("overwrite=false") && expected.getMessage().contains(prefix+".gc_cost.tsv");
		}
		check(rejected, "Cost inference must refuse to overwrite a completed output or its source tables");
		System.err.println("NCRNA_DDL_GC_RETAINED_FAILURE_COST_LITERAL_TEST_PASS");
	}

	private static void gcOrdering(){
		check(BuildNcrnaDdlLayoutCounts.gcRemovedMedian(new long[]{1,0,0,1})==1.5
			&& BuildNcrnaDdlLayoutCounts.gcRemovedMedian(new long[]{0,1,2,0})==2,
			"Removed-model median includes zero and uses the midpoint for an even number of windows");
		final double[] combinedStats={0.5,0.125,0.375,0.625};
		check(BuildNcrnaDdlLayoutCounts.gcExtraKeep(0.625, combinedStats, 0)
			&& !BuildNcrnaDdlLayoutCounts.gcExtraKeep(Math.nextUp(0.625), combinedStats, 0)
			&& BuildNcrnaDdlLayoutCounts.gcExtraKeep(0.63, combinedStats, 3)
			&& BuildNcrnaDdlLayoutCounts.gcExtraScore(10, 2, 0)==9
			&& BuildNcrnaDdlLayoutCounts.gcExtraScore(10, 2, 1)==8
			&& BuildNcrnaDdlLayoutCounts.gcExtraScore(10, 2, 2)==6,
			"Combined ranking gates with the exact range/margin first and then subtracts the declared half/one/two-weight penalty");
		final double[] stats=NcrnaDdlSavedData.gcStatistics(1)[0];
		NcrnaDdlSavedData.addGc(stats, 0.5);
		final double mismatch=BuildNcrnaDdlLayoutCounts.gcDistance(0.25, stats[1], 0);
		check(mismatch==Double.POSITIVE_INFINITY && BuildNcrnaDdlLayoutCounts.gcDistance(0.5, stats[1], 0)==0,
			"Zero-SD mismatch is a distinct bottom tier; exact mean has zero penalty without epsilon");
		check(BuildNcrnaDdlLayoutCounts.gcScore(5, mismatch, 0)==5
			&& BuildNcrnaDdlLayoutCounts.gcScore(5, mismatch, 0.5)==Double.NEGATIVE_INFINITY,
			"w0 must bypass infinite distance rather than multiply zero by infinity");
		check(BuildNcrnaDdlLayoutCounts.gcRank(new double[]{-9,-9,Double.NEGATIVE_INFINITY}, new long[]{4}, new boolean[]{true,true,true})==3
			&& BuildNcrnaDdlLayoutCounts.gcRank(new double[]{0,0,0}, new long[]{2}, new boolean[]{true,false,true})==4
			&& BuildNcrnaDdlLayoutCounts.gcRank(new double[]{2,2,2}, new long[]{2}, new boolean[]{true,true,true})==2,
			"GC ranks retain signed scores,bottom-tier ordering,ordinal ties and an explicit all-filtered failure");
		final double[] effective={0.5,0,0.5,0.5}, missing={Double.NaN,Double.NaN,Double.NaN,Double.NaN};
		check(!BuildNcrnaDdlLayoutCounts.gcKeep(0.25, effective, 6)
			&& BuildNcrnaDdlLayoutCounts.gcKeep(0.5, effective, 8)
			&& BuildNcrnaDdlLayoutCounts.gcKeep(0.25, missing, 6)
			&& BuildNcrnaDdlLayoutCounts.gcKeep(Double.NaN, effective, 9),
			"Zero-SD filters are exact; absent fit observations or missing query GC cannot invent rejection evidence");
		final double[] spread={0.5,0.03125,0.4,0.6};
		for(final double value:new double[]{0.25,0.375,0.4375,0.5,0.5625,0.625,0.75}){
			for(int setting=6; setting<9; setting++){
				check(BuildNcrnaDdlLayoutCounts.gcAllKeep(value, spread, setting, 9)==BuildNcrnaDdlLayoutCounts.gcKeep(value, spread, setting),
					"All-data and historical2/3/4SD filters share the exact same inclusive arithmetic");
			}
		}
		check(BuildNcrnaDdlLayoutCounts.gcFirstWidth(0.25, effective, 29, 30)==0
			&& BuildNcrnaDdlLayoutCounts.gcFirstWidth(0.25, effective, 30, 30)==25
			&& BuildNcrnaDdlLayoutCounts.gcFirstWidth(0.5, effective, 30, 30)==0,
			"Thirty observations enable an exact zero-SD gate;29disable it and a mismatch has no finite lossless width");
		final double[] varied={0.5,0.03125,0.4,0.6};
		for(final double value:new double[]{0.25,0.375,0.4375,0.5,0.5625,0.625,0.75}){
			final int first=BuildNcrnaDdlLayoutCounts.gcFirstWidth(value, varied, 30, 0);
			for(int setting=6; setting<9; setting++){
				check((first<=(setting-5)*4)==BuildNcrnaDdlLayoutCounts.gcKeep(value, varied, setting),
					"First-passing-width search preserves the original inclusive z2/3/4 boundary decisions");
			}
		}
		System.err.println("NCRNA_DDL_GC_ZERO_SD_AND_FILTER_FAILURE_TEST_PASS");
	}

	/** Alpha has seven A-core and two G-core accepted windows; beta has one
	 * A-core window. Their different observed GC must stay on opposite fit sides. */
	private static void gcIsolation(final BuildNcrnaDdlLayoutCounts.Config c, final byte[] padding){
		final byte[] a=new byte[20+padding.length], g=new byte[a.length];
		Arrays.fill(a, 0, 20, (byte)'A');
		Arrays.fill(g, 0, 20, (byte)'G');
		System.arraycopy(padding, 0, a, 20, padding.length);
		System.arraycopy(padding, 0, g, 20, padding.length);
		final double beta=shared.Tools.calcGC(a), alpha=(7*beta+2*shared.Tools.calcGC(g))/9;
		check(Math.abs(alpha-beta)>0.000001, "The isolation witness requires distinguishable real fit/judge GC distributions");
		final LineParser1 p=new LineParser1('\t');
		int alphaFold=-1, betaFold=-1;
		try(Input in=new Input(c.output+".folds.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(p.parseInt(0)==0){alphaFold=p.parseInt(2);}
				if(p.parseInt(0)==2){betaFold=p.parseInt(2);}
			}
		}
		check(alphaFold>=0 && betaFold>=0 && alphaFold!=betaFold, "The two nonempty assemblies must occupy distinct judging folds");
		int checked=0;
		try(Input in=new Input(c.output+".gc_fit.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==15, "GC fit evidence retains actual model/fold counts and both fallback sources");
				if(p.parseInt(1)!=0 || !p.termEquals("window", 2)){continue;}
				final boolean judgingAlpha=p.parseInt(0)==alphaFold;
				check(p.parseInt(3)==(judgingAlpha ? 1 : 9)
					&& Math.abs(p.parseDouble(7)-(judgingAlpha ? beta : alpha))<0.00000000051,
					"Emitted model window GC must equal only the opposite assembly, not the combined or judging distribution");
				check(p.parseInt(4)==0, "The exact-duplicate held-out model must not regain member GC through the fallback arm");
				checked++;
			}
		}
		check(checked==2, "Both judging directions must execute the fit-GC isolation witness");
		System.err.println("NCRNA_DDL_GC_FIT_ISOLATION_TEST_PASS");
	}

	/** A known noncontributor fails every filter; no locus may disappear from p99. */
	private static void gcReport(final String prefix){
		final String roster=prefix+".roster.tsv";
		final structures.ByteBuilder line=new structures.ByteBuilder(8192);
		try(Output out=new Output(roster)){
			out.writer.println("assignment_seed\tddl_seed\ttable");
			for(int assignment=0; assignment<3; assignment++){
				for(int seed=0; seed<2; seed++){
					final long a=assignment==2 ? -1 : assignment+11, s=seed+101;
					final String path=prefix+"_"+assignment+"_"+seed+".tsv";
					out.writer.println(a+"\t"+s+"\t"+path);
					try(Output table=new Output(path)){
						table.writer.println(BuildNcrnaDdlLayoutCounts.GC_REPLICATE_HEADER);
						for(int fold=0; fold<(a<0 ? 2 : 5); fold++){
							for(final int width:new int[]{512,1024}){
								for(int setting=0; setting<24; setting++){
									for(int stratum=0; stratum<(a<0 ? 1 : 3); stratum++){
										final int count=stratum==0 ? 2 : 1, normal=seed+1;
										final int failures=setting%12>=6 && stratum!=1 ? 1 : 0;
										final int rankSum=(count-failures)*normal+failures*509, top=count-failures;
										line.clear().append(a).tab().append(s).tab().append(fold).tab().append(BuildNcrnaDdlLayoutCounts.gcMethod(width, setting))
											.append("\t8\t10\tquarter\t").append(width).append("\t256\t2048\t").append(width==512 ? 60 : 100).append("\t0")
											.tab().append(count).tab().append(count).tab().append(rankSum).tab().append(rankSum)
											.tab().append(top).tab().append(top).tab().append(top).tab().append(top)
											.append("\t0\tNA\tNA\tNA\t0\t0\tGC_FIXTURE\t").append(BuildNcrnaDdlLayoutCounts.CONTRIBUTION[stratum]);
										for(int low=0; low<2; low++){
											line.tab();
											for(int rank=1; rank<=509; rank++){
												if(rank>1){line.append(',');}
												line.append(rank==normal ? count-failures : rank==509 ? failures : 0);
											}
										}
										line.tab().append(stratum==0 ? 2 : 0).tab().append(stratum==0 && setting%12>=6 ? 1 : 0)
											.tab().append(stratum==0 ? 1 : 0).tab().append(stratum==0 && setting%12>=6 ? 1 : 0).append("\t0");
										for(final int value:new int[]{3,1,2,0}){line.tab().append(stratum==0 ? value : 0);}
										table.writer.print(line.nl());
									}
								}
							}
						}
					}
				}
			}
		}
		prok.SummarizeNcrnaDdlReplicates.main(new String[]{"gc=t", "in="+roster, "out="+prefix, "seeds=101,102", "foldseeds=11,12", "t=1"});
		final LineParser1 p=new LineParser1('\t');
		boolean found=false;
		try(Input in=new Input(prefix+".summary.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(p.termEquals("REBUILT_FIVE_FOLD", 0) && p.termEquals("gc_512_window_filter_z2", 1) && p.termEquals("ALL", 13)){
					if(p.termEquals("rank_p99", 2)){check(p.parseDouble(5)==509, "Failed acceptable models must remain in the pooled percentile tail"); found=true;}
					if(p.termEquals("parents_all_acceptable_filtered", 2)){check(p.parseDouble(5)==5, "All five failed noncontributor loci must remain explicit");}
				}
			}
		}
		check(found, "The GC reporter must actually emit the sentinel-aware pooled percentile");
		final String actual=prefix+"_0_0.tsv", control=prefix+".saved_control.tsv", corrupt=prefix+".wrong_control.tsv";
		try(Input in=new Input(actual); Output good=new Output(control); Output bad=new Output(corrupt)){
			in.reader.nextLine();
			good.writer.println(BuildNcrnaDdlLayoutCounts.FOLLOWUP_REPLICATE_HEADER);
			bad.writer.println(BuildNcrnaDdlLayoutCounts.FOLLOWUP_REPLICATE_HEADER);
			boolean changed=false;
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(!p.termEquals(BuildNcrnaDdlLayoutCounts.gcMethod(1024, 0), 3)){continue;}
				final structures.ByteBuilder correct=new structures.ByteBuilder(), wrong=new structures.ByteBuilder();
				for(int field=0; field<30; field++){
					String value=field==3 ? BuildNcrnaDdlLayoutCounts.sizeMethod(256, 2048) : p.parseString(field);
					if(field>=28){value=value.substring(0, value.length()-2);}
					if(field>0){correct.tab(); wrong.tab();}
					correct.append(value);
					wrong.append(field==14 && !changed ? Long.toString(p.parseLong(field)+1) : value);
				}
				changed=true;
				good.writer.print(correct.nl());
				bad.writer.print(wrong.nl());
			}
		}
		prok.SummarizeNcrnaDdlReplicates.checkGcControl(actual, control, prefix+".control");
		boolean rejected=false;
		try{prok.SummarizeNcrnaDdlReplicates.checkGcControl(actual, corrupt, prefix+".wrong");}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("changed a saved");}
		check(rejected, "An altered saved rank sum must fail exact GC-control parity rather than silently pass aggregate populations");
		System.err.println("NCRNA_DDL_GC_CROSSED_REPORT_TEST_PASS");
		gcAllReport(prefix+".all_data", roster);
	}

	/** Reuses synthetic saved rows to test64-method crossed reporting and control
	 * identity. The actual all-data scorer is checked by its retained native pilot. */
	private static void gcAllReport(final String prefix, final String previousRoster){
		final String roster=prefix+".roster.tsv";
		final LineParser1 p=new LineParser1('\t');
		final structures.ByteBuilder b=new structures.ByteBuilder(8192);
		final int[] previousSettings={0,1,2,3,4,5,6,7,8,8,8,0,9,10,10,11};
		String first=null, previousFirst=null;
		try(Input old=new Input(previousRoster); Output out=new Output(roster)){
			out.writer.println("assignment_seed\tddl_seed\ttable");
			old.reader.nextLine();
			for(byte[] binding; (binding=old.reader.nextLine())!=null;){
				p.set(binding);
				final long assignment=p.parseLong(0), seed=p.parseLong(1);
				final String source=p.parseString(2), target=source+".all_data.tsv";
				if(first==null){first=target; previousFirst=source;}
				out.writer.println(assignment+"\t"+seed+"\t"+target);
				try(Input in=new Input(source); Output table=new Output(target)){
					table.writer.println(BuildNcrnaDdlLayoutCounts.GC_REPLICATE_HEADER);
					in.reader.nextLine();
					for(byte[] row; (row=in.reader.nextLine())!=null;){
						p.set(row);
						final int width=p.parseInt(7);
						for(int arm=0; arm<2; arm++){
							for(int setting=0; setting<16; setting++){
								final int oldSetting=arm==0 && setting>=12 ? 0 : previousSettings[setting];
								if(!p.termEquals(BuildNcrnaDdlLayoutCounts.gcMethod(width, arm*12+oldSetting), 3)){continue;}
								p.setBounds(3);
								b.clear().append(row, 0, p.a()).append(BuildNcrnaDdlLayoutCounts.gcAllMethod(width, arm*16+setting))
									.append(row, p.b(), row.length-p.b());
								table.writer.print(b.nl());
							}
						}
					}
				}
			}
		}
		prok.SummarizeNcrnaDdlReplicates.main(new String[]{"gcall=t", "in="+roster, "out="+prefix, "seeds=101,102", "foldseeds=11,12", "t=1"});
		boolean seeded49=false;
		try(Input in=new Input(prefix+".summary.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(p.termEquals("REBUILT_FIVE_FOLD", 0) && p.termEquals(BuildNcrnaDdlLayoutCounts.gcAllMethod(1024, 17), 1)
					&& p.termEquals("rank_sum", 2) && p.termEquals("ALL", 13)){
					check(p.parseInt(4)==2 && p.parseDouble(5)==15 && p.parseDouble(6)>0,
						"GC method49 is a seeded alternative,not the legacy unseeded index reference;both seed observations must contribute");
					seeded49=true;
				}
			}
		}
		check(seeded49, "The64-method reporter fixture must inspect the actual second-bin member-weight0.5 estimate");
		prok.SummarizeNcrnaDdlReplicates.checkGcAllControl(first, previousFirst, prefix+".parity");
		gcExtraReport(prefix+".extension", roster, prefix+".summary.tsv");
		System.err.println("NCRNA_DDL_ALL_DATA_GC_CROSSED_REPORT_AND_CONTROL_PARITY_TEST_PASS");
	}

	/** Import saved observations once;extension tables contain no old setting. */
	private static void gcExtraReport(final String prefix, final String originalRoster, final String originalSummary){
		final String roster=prefix+".roster.tsv";
		final LineParser1 p=new LineParser1('\t');
		final structures.ByteBuilder b=new structures.ByteBuilder(8192);
		final map.ObjectSet<String> extra=new map.ObjectSet<String>(String.class);
		for(final int width:new int[]{512,1024}){
			for(int setting=0; setting<BuildNcrnaDdlLayoutCounts.GC_EXTRA_SETTINGS.length; setting++){
				extra.add(BuildNcrnaDdlLayoutCounts.gcExtraMethod(width, setting));
			}
		}
		try(Input old=new Input(originalRoster); Output out=new Output(roster)){
			out.writer.println(old.reader.nextLine());
			for(byte[] binding; (binding=old.reader.nextLine())!=null;){
				p.set(binding);
				final long assignment=p.parseLong(0), seed=p.parseLong(1);
				final String source=p.parseString(2), target=source+".extension.tsv";
				out.writer.println(assignment+"\t"+seed+"\t"+target);
				try(Input in=new Input(source); Output table=new Output(target)){
					table.writer.println(in.reader.nextLine());
					for(byte[] row; (row=in.reader.nextLine())!=null;){
						p.set(row);
						final int width=p.parseInt(7);
						if(!p.termEquals(BuildNcrnaDdlLayoutCounts.gcAllMethod(width, 0), 3)){continue;}
						p.setBounds(3);
						for(int setting=0; setting<BuildNcrnaDdlLayoutCounts.GC_EXTRA_SETTINGS.length; setting++){
							table.writer.print(b.clear().append(row, 0, p.a()).append(BuildNcrnaDdlLayoutCounts.gcExtraMethod(width, setting))
								.append(row, p.b(), row.length-p.b()).nl());
						}
					}
				}
			}
		}
		prok.SummarizeNcrnaDdlReplicates.main(new String[]{"gcall=t", "in="+originalRoster,
			"gcallfollowup="+roster, "out="+prefix, "seeds=101,102", "foldseeds=11,12", "t=1"});
		int imported=0, added=0;
		try(Input actual=new Input(prefix+".summary.tsv"); Input expected=new Input(originalSummary)){
			check(Arrays.equals(actual.reader.nextLine(), expected.reader.nextLine()), "An extension preserves the original report schema");
			for(byte[] row; (row=actual.reader.nextLine())!=null;){
				p.set(row);
				if(extra.contains(p.parseString(1))){
					added++;
					continue;
				}
				check(Arrays.equals(row, expected.reader.nextLine()),
					"Every imported estimate,including controls and SD4/6,must remain byte-identical with all saved seed/assignment observations");
				imported++;
			}
			check(expected.reader.nextLine()==null && imported>0 && added>0,
				"The extension must retain every original estimate and add actual new estimates");
		}
		for(int corruption=0; corruption<3; corruption++){
			final String bad=prefix+".bad"+corruption;
			try(Input in=new Input(roster); Input originals=new Input(originalRoster); Output out=new Output(bad+".roster.tsv")){
				out.writer.println(in.reader.nextLine());
				originals.reader.nextLine();
				final byte[] first=in.reader.nextLine(), original=originals.reader.nextLine();
				if(corruption==0){
					out.writer.println(first);
					out.writer.println(first);
				}
				else if(corruption==2){out.writer.println(original);}
				for(byte[] row; (row=in.reader.nextLine())!=null;){out.writer.println(row);}
			}
			boolean rejected=false;
			try{
				prok.SummarizeNcrnaDdlReplicates.main(new String[]{"gcall=t", "in="+originalRoster,
					"gcallfollowup="+bad+".roster.tsv", "out="+bad, "seeds=101,102", "foldseeds=11,12", "t=1"});
			}catch(IllegalArgumentException expected){rejected=true;}
			check(rejected && !new File(bad+".summary.tsv").exists(),
				"Duplicate or missing extension workers,and old controls relabeled as new work,must fail before reporting");
		}
		System.err.println("NCRNA_DDL_GC_EXTENSION_SAVED_ROWS_ONCE_TEST_PASS");
	}

	/** Unequal folds distinguish a pooled percentile from averaged percentiles. */
	private static void percentiles(){
		final double[][] best={{9,11},{9,11}}, near={{9,13},{9,13}}, outside={{11,13},{11,13}};
		final int[] choice=prok.SummarizeNcrnaDdlReplicates.chooseSize(new double[][][]{best,near,outside,near},
			new int[]{256,256,0,0}, new int[]{4096,256,256,512});
		check(Arrays.equals(choice, new int[]{0,3}),
			"Choose the smaller total list within one PAIRED SE;reject the yet smaller but precisely worse list,and prefer fewer short entries on equal totals");
		final long[] pooled=new long[100], metrics=new long[4];
		pooled[0]=99;
		pooled[9]=1;
		check(prok.SummarizeNcrnaDdlReplicates.percentile99(pooled)==1,
			"Nearest-rank99th percentile of100loci uses observation99,not an interpolated tail");
		pooled[99]=1;
		check(prok.SummarizeNcrnaDdlReplicates.percentile99(pooled)==10,
			"Pooling one high-rank locus with100prior loci yields rank10;averaging fold percentiles would incorrectly give50.5");
		final byte[] first="99,1,0".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
		final byte[] second="0,0,1".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
		final long[] bins=new long[3];
		prok.SummarizeNcrnaDdlReplicates.addHistogram(first, 0, first.length, bins, metrics);
		prok.SummarizeNcrnaDdlReplicates.addHistogram(second, 0, second.length, bins, metrics);
		check(Arrays.equals(bins, new long[]{99,1,1}) && Arrays.equals(metrics, new long[]{1,3,1,1}),
			"Byte-range pooling reuses bins while per-fold conservation scratch resets");
		for(final String bad:new String[]{"1,2", "1,2,3,4", "1,,2", "1,-1,2"}){
			boolean rejected=false;
			final byte[] bytes=bad.getBytes(java.nio.charset.StandardCharsets.US_ASCII);
			try{prok.SummarizeNcrnaDdlReplicates.addHistogram(bytes, 0, bytes.length, new long[3], metrics);}
			catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "Missing,extra,empty and negative rank bins cannot produce a percentile");
		}
	}

	/** Literal additive effects exercise whole-table fold pooling and uncertainty.
	 * Padding rows have a known nonzero paired effect;these are toy counts. */
	private static void replicateReport(final String prefix){
		final String roster=prefix+".roster.tsv";
		try(Output out=new Output(roster)){
			out.writer.println("assignment_seed\tddl_seed\ttable");
			for(int assignment=0; assignment<3; assignment++){
				for(int seed=0; seed<2; seed++){
					final long a=assignment==2 ? -1 : assignment+11, s=seed+101;
					final String path=prefix+"_"+assignment+"_"+seed+".tsv";
					out.writer.println(a+"\t"+s+"\t"+path);
					try(Output table=new Output(path)){
						table.writer.println(assignment==2 && seed==0 ? BuildNcrnaDdlLayoutCounts.LEGACY_REPLICATE_HEADER : BuildNcrnaDdlLayoutCounts.REPLICATE_HEADER);
						for(int fold=0; fold<(assignment==2 ? 2 : 5); fold++){
							for(int pair=0; pair<3; pair++){
								for(final int bins:new int[]{512,1024}){
									for(int layout=0; layout<3; layout++){
										for(int filtered=0; filtered<2; filtered++){
											final String name="s"+(pair==2 ? 8 : 7)+"_m"+(pair==0 ? 9 : 10)+"_"+new String[]{"mixed","half","quarter"}[layout]+"_"+bins+"_"+(filtered==0 ? "plain" : "padding");
											replicateRow(table, a, s, fold, name, pair, bins, layout, filtered, 100, 0, assignment, seed);
										}
									}
								}
							}
							for(final int bonus:new int[]{0,5,1000000}){replicateRow(table, a, s, fold, "bonus_fixed_"+bonus, 0, 512, 0, 0, bonus, 0, assignment, seed);}
							for(int tenth=1; tenth<=10; tenth++){replicateRow(table, a, s, fold, "bonus_tenths_"+tenth, 0, 512, 0, 0, 0, tenth/10.0, assignment, seed);}
							if(seed==0){replicateRow(table, a, -1, fold, "CONSENSUS_INDEX", 0, 512, 0, 0, 100, 0, assignment, 0);}
						}
					}
				}
			}
		}
		prok.SummarizeNcrnaDdlReplicates.main(new String[]{"in="+roster, "out="+prefix, "seeds=101,102", "foldseeds=11,12", "t=1"});
		try(Input in=new Input(prefix+".summary.tsv")){
			in.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			p.set(in.reader.nextLine());
			check(p.terms()==14 && p.termEquals("ALL", 13) && p.termEquals("rank_sum", 2) && p.parseInt(3)==2 && p.parseInt(4)==2
				&& p.parseDouble(5)==130 && p.parseDouble(7)==50 && p.parseDouble(8)==1250
				&& p.parseDouble(9)==0 && Math.abs(p.parseDouble(10)-Math.sqrt(650))<1e-8
				&& p.parseInt(11)==10 && p.parseInt(12)==5,
				"Whole report must pool five folds first and preserve the literal crossed effects and parent denominators");
			boolean contributor=false, noncontributor=false;
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(!p.termEquals("REBUILT_FIVE_FOLD", 0) || !p.termEquals("s7_m9_mixed_512_plain", 1) || !p.termEquals("rank_sum", 2)){continue;}
				if(p.termEquals("CONTRIBUTOR", 13)){
					check(p.parseDouble(5)==25 && p.parseDouble(10)==0 && p.parseInt(11)==5 && p.parseInt(12)==5,
						"The contributor subset pools its own five constant low-identity loci before uncertainty");
					contributor=true;
				}else if(p.termEquals("NONCONTRIBUTOR", 13)){
					check(p.parseDouble(5)==105 && Math.abs(p.parseDouble(10)-Math.sqrt(650))<1e-8
						&& p.parseInt(11)==5 && p.parseInt(12)==0,
						"The complementary subset conserves total rank mean and carries the literal crossed effects");
					noncontributor=true;
				}
			}
			check(contributor && noncontributor, "Whole-table reporting must contain both original contribution populations");
		}
		try(Input in=new Input(prefix+".paired.tsv")){
			in.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			p.set(in.reader.nextLine());
			check(p.parseDouble(6)==-10 && Math.abs(p.parseDouble(7)-Math.sqrt(12.5))<1e-8,
				"Nonzero paired method differences must preserve their own literal seed/split effects and crossed standard error");
		}
		followupReport(prefix, roster);
		rejectContributionReport(prefix, roster);
	}

	/** Reuses the accepted finite old roster and adds literal rank distributions. */
	private static void followupReport(final String prefix, final String originalRoster){
		final String roster=prefix+".followup.roster.tsv", report=prefix+".followup";
		try(Output out=new Output(roster)){
			out.writer.println("assignment_seed\tddl_seed\toriginal_table\ttable");
			for(int assignment=0; assignment<3; assignment++){
				for(int seed=0; seed<2; seed++){
					final long a=assignment==2 ? -1 : assignment+11, s=seed+101;
					final String original=prefix+"_"+assignment+"_"+seed+".tsv", path=report+"_"+assignment+"_"+seed+".tsv";
					out.writer.println(a+"\t"+s+"\t"+original+"\t"+path);
					try(Output table=new Output(path)){
						table.writer.println(BuildNcrnaDdlLayoutCounts.FOLLOWUP_REPLICATE_HEADER);
						for(int fold=0; fold<(a<0 ? 2 : 5); fold++){
							for(final int bins:new int[]{512,1024}){
								for(final int bonus:new int[]{60,80,100,120}){
									final String name="s8_m10_quarter_"+bins+"_padding"+(bonus==100 ? "" : "_bonus"+bonus);
									replicateRow(table, a, s, fold, name, 2, bins, 2, 1, bonus, 0, assignment, seed, true);
								}
								for(int layout=0; layout<3; layout++){
									for(int filter=0; filter<2; filter++){
										final String name="s8_m11_"+new String[]{"mixed","half","quarter"}[layout]+"_"+bins+"_"+(filter==0 ? "plain" : "padding");
										replicateRow(table, a, s, fold, name, 3, bins, layout, filter, 100, 0, assignment, seed, true);
									}
								}
							}
						}
					}
				}
			}
		}
		prok.SummarizeNcrnaDdlReplicates.main(new String[]{"in="+originalRoster, "followup="+roster,
			"out="+report, "seeds=101,102", "foldseeds=11,12", "t=1"});
		boolean all=false, low=false, control=false;
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(report+".summary.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(!p.termEquals("REBUILT_FIVE_FOLD", 0) || !p.termEquals("ALL", 13)){continue;}
				if(p.termEquals("s8_m11_mixed_512_plain", 1) && p.termEquals("rank_p99", 2)){
					check(p.parseDouble(5)==21 && Math.abs(p.parseDouble(6)-Math.sqrt(104.0/3))<1e-8,
						"Pooled per-cell tail ranks15/17/25/27 need their own mean and SD");
					all=true;
				}
				if(p.termEquals("s8_m11_mixed_512_plain", 1) && p.termEquals("low_rank_p99", 2)){
					check(p.parseDouble(5)==5 && p.parseDouble(6)==0, "Low-identity percentiles must use only their own loci");
					low=true;
				}
				if(p.termEquals("s8_m10_quarter_512_padding", 1) && p.termEquals("rank_p99", 2)){
					check(p.parseDouble(5)==23, "Repeated100controls supply only the new percentile while exact old metrics remain reused");
					control=true;
				}
			}
		}
		check(all && low && control, "The full follow-up report needs both tail populations and the repeated control");
		gridReport(prefix, originalRoster, roster);
		System.err.println("NCRNA_DDL_FOLLOWUP_AND_POOLED_PERCENTILE_TEST_PASS");
	}

	/** The72-cell production roster imports the winner once and preserves the
	 * no-blacklist cell. Toy means tie,so the literal one-SE choice is0/0. */
	private static void gridReport(final String prefix, final String originalRoster, final String followupRoster){
		final String grid=prefix+".grid.roster.tsv", definitions=prefix+".definitions.roster.tsv", extension=prefix+".extension.roster.tsv";
		for(int kind=0; kind<3; kind++){
			try(Output out=new Output(kind==0 ? grid : kind==1 ? definitions : extension)){
				out.writer.println("assignment_seed\tddl_seed\toriginal_table\ttable");
				for(int assignment=0; assignment<3; assignment++){
					for(int seed=0; seed<2; seed++){
						final long a=assignment==2 ? -1 : assignment+11, s=seed+101;
						final String path=prefix+".grid"+kind+"_"+assignment+"_"+seed+".tsv";
						out.writer.println(a+"\t"+s+"\t"+prefix+"_"+assignment+"_"+seed+".tsv\t"+path);
						try(Output table=new Output(path)){
							table.writer.println(BuildNcrnaDdlLayoutCounts.FOLLOWUP_REPLICATE_HEADER);
							for(int fold=0; fold<(a<0 ? 2 : 5); fold++){
								if(kind==0){
									for(final int shortSize:BuildNcrnaDdlLayoutCounts.GRID_SHORT){
										for(final int longSize:BuildNcrnaDdlLayoutCounts.GRID_LONG){
											if(shortSize==256 && longSize==4096){continue;}
											final String name=shortSize==0 && longSize==0 ? "s8_m10_quarter_1024_plain" : BuildNcrnaDdlLayoutCounts.sizeMethod(shortSize, longSize);
											replicateRow(table, a, s, fold, name, 2, 1024, 2, 1, 100, 0, assignment, seed, true, shortSize, longSize);
										}
									}
								}else if(kind==2){
									for(final int[] size:BuildNcrnaDdlLayoutCounts.SIZE_EXTENSION){
										replicateRow(table, a, s, fold, BuildNcrnaDdlLayoutCounts.sizeMethod(size[0], size[1]),
											2, 1024, 2, 1, 100, 0, assignment, seed, true, size[0], size[1]);
									}
								}else{
									for(final String definition:new String[]{"false","pooled","members"}){
										replicateRow(table, a, s, fold, "blacklist_"+definition+"_s0_l0", 2, 1024, 2, 1, 100, 0, assignment, seed, true, 0, 0);
									}
								}
							}
						}
					}
				}
			}
		}
		final String report=prefix+".grid_report";
		prok.SummarizeNcrnaDdlReplicates.main(new String[]{"in="+originalRoster, "followup="+followupRoster,
			"sizegrid="+grid, "definitions="+definitions, "shortsize=0", "longsize=0", "out="+report,
			"seeds=101,102", "foldseeds=11,12", "t=1"});
		check(rows(report+".sizes.tsv")==145, "The size table contains72cells in each population,including the single imported winner");
		int chosen=0;
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(report+".sizes.tsv")){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(!p.termEquals("true", 15)){continue;}
				check(p.parseInt(1)==0 && p.parseInt(2)==0, "Equal full means select the fewest entries,including the no-blacklist option");
				chosen++;
			}
		}
		check(chosen==2, "The full-population choice must be shown on both full and fixed tables");
		prok.SummarizeNcrnaDdlReplicates.main(new String[]{"in="+originalRoster, "followup="+followupRoster,
			"sizegrid="+grid, "sizeextension="+extension, "definitions="+definitions, "shortsize=0", "longsize=0", "out="+report+"_extended",
			"seeds=101,102", "foldseeds=11,12", "t=1"});
		check(BuildNcrnaDdlLayoutCounts.SIZE_EXTENSION.length==22 && rows(report+"_extended.sizes.tsv")==189,
			"The combined report preserves all72old cells and appends22distinct new cells in both populations");
		System.err.println("NCRNA_DDL_EXTENDED_SIZE_GRID_TEST_PASS");
		System.err.println("NCRNA_DDL_SIZE_GRID_AND_ONE_SE_SELECTION_TEST_PASS");
	}

	/** A changed subset must not publish totals; legacy rebuilt tables are barred. */
	private static void rejectContributionReport(final String prefix, final String roster){
		final String source=prefix+"_0_0.tsv", changed=prefix+".changed.tsv";
		final LineParser1 p=new LineParser1('\t');
		boolean damaged=false;
		try(Input in=new Input(source); Output out=new Output(changed)){
			out.writer.println(in.reader.nextLine());
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				if(!damaged && p.termEquals("CONTRIBUTOR", 27)){
					final structures.ByteBuilder b=new structures.ByteBuilder();
					for(int field=0; field<28; field++){
						if(field>0){b.tab();}
						if(field==14){b.append(p.parseLong(field)+1);}
						else{p.setBounds(field); b.append(row, p.a(), p.b()-p.a());}
					}
					out.writer.print(b.nl());
					damaged=true;
				}else{out.writer.println(row);}
			}
		}
		check(damaged, "The conservation fixture must actually alter a contributor rank sum");
		for(int fault=0; fault<2; fault++){
			final String badRoster=prefix+".bad_"+fault+".roster.tsv", badOutput=prefix+".bad_"+fault;
			try(Input in=new Input(roster); Output out=new Output(badRoster)){
				for(byte[] row; (row=in.reader.nextLine())!=null;){
					out.writer.println(new String(row, java.nio.charset.StandardCharsets.UTF_8).replace(source, fault==0 ? changed : prefix+".legacy.tsv"));
				}
			}
			if(fault==1){
				try(Output out=new Output(prefix+".legacy.tsv")){out.writer.println(BuildNcrnaDdlLayoutCounts.LEGACY_REPLICATE_HEADER);}
			}
			boolean rejected=false;
			try{
				prok.SummarizeNcrnaDdlReplicates.main(new String[]{"in="+badRoster, "out="+badOutput,
					"seeds=101,102", "foldseeds=11,12", "t=1"});
			}catch(IllegalArgumentException expected){
				rejected=expected.getMessage().contains(fault==0 ? "partition every pooled locus" : "Only accepted fixed-panel observations");
			}
			check(rejected && !new File(badOutput+".summary.tsv").exists(),
				"Wrong contribution totals and legacy rebuilt observations must fail before a report is published");
		}
	}

	private static void replicateRow(final Output out, final long assignment, final long seed, final int fold,
			final String method, final int pair, final int bins, final int layout, final int filtered,
			final int bonus, final double scale, final int a, final int s){
		replicateRow(out, assignment, seed, fold, method, pair, bins, layout, filtered, bonus, scale, a, s, false);
	}
	private static void replicateRow(final Output out, final long assignment, final long seed, final int fold,
			final String method, final int pair, final int bins, final int layout, final int filtered,
			final int bonus, final double scale, final int a, final int s, final boolean histogram){
		replicateRow(out, assignment, seed, fold, method, pair, bins, layout, filtered, bonus, scale, a, s, histogram, -1, -1);
	}
	private static void replicateRow(final Output out, final long assignment, final long seed, final int fold,
			final String method, final int pair, final int bins, final int layout, final int filtered,
			final int bonus, final double scale, final int a, final int s, final boolean histogram,
			final int shortOverride, final int longOverride){
		final int extra=method.endsWith("_padding") ? 1+(a==2 ? 0 : a)+s : 0;
		final int total=20+(seed<0 || a==2 ? 0 : 10*a)+(seed<0 ? 0 : 2*s)+extra;
		for(int stratum=0; stratum<(assignment<0 ? 1 : 3); stratum++){
			final structures.ByteBuilder b=new structures.ByteBuilder();
			b.append(assignment).tab().append(seed).tab().append(fold).tab().append(method)
				.tab().append(pair>=2 ? 8 : 7).tab().append(pair==0 ? 9 : pair==3 ? 11 : 10)
				.tab().append(new String[]{"mixed","half","quarter"}[layout]).tab().append(bins)
				.tab().append(shortOverride>=0 ? shortOverride : filtered==0 ? 0 : pair>=2 ? 256 : 64)
				.tab().append(longOverride>=0 ? longOverride : filtered==0 ? 0 : 4096)
				.tab().append(bonus).tab().append(scale, 9).tab().append(stratum==0 ? 2 : 1)
				.tab().append(stratum==2 ? 0 : 1).tab().append(stratum==0 ? total : stratum==1 ? 5 : total-5)
				.tab().append(stratum==2 ? 0 : 5).tab().append(stratum==2 ? 0 : 1).tab().append(stratum==0 ? 2 : 1)
				.tab().append(stratum==2 ? 0 : 1).tab().append(stratum==2 ? 0 : 1).tab().append(0)
				.append(seed<0 || stratum>0 ? "\tNA\tNA\tNA" : "\t0.5\t0.5\t0.5").append("\t0\t0\tSYNTHETIC_CROSSED_REPORT_FIXTURE");
			if(histogram || !(a==2 && s==0)){b.tab().append(BuildNcrnaDdlLayoutCounts.CONTRIBUTION[stratum]);}
			if(histogram){
				for(int low=0; low<2; low++){
					b.tab();
					for(int rank=1; rank<=508; rank++){
						if(rank>1){b.append(',');}
						b.append((stratum!=2 && rank==5 ? 1 : 0)+(low==0 && stratum!=1 && rank==total-5 ? 1 : 0));
					}
				}
			}
			b.nl();
			out.writer.print(b);
		}
	}

	/** The omitted public flag and explicit XOR flag must select the same actual
	 * scientific output;elapsed/CPU cost columns are deliberately not compared. */
	private static void xorDefault(final RankDualNcrnaWindows.Config saved, final String control, final String recipes){
		final String[] outputs=new String[2];
		for(int mode=0; mode<2; mode++){
			final BuildNcrnaDdlLayoutCounts.Config c=config(saved, control, "_xor_default_"+mode);
			outputs[mode]=c.output;
			final ArrayList<String> args=new ArrayList<String>(Arrays.asList("recipes="+recipes,
				"panels="+c.panels, "selection="+c.selection, "ref="+c.references, "members="+c.members,
				"ledger="+c.ledger, "assignments="+c.assignments, "census="+c.census, "summary="+c.summary,
				"out="+c.output, "identities="+saved.identities, "study=baseline", "expectedmodels=2", "expectedgroups=3", "t=1"));
			if(mode==1){args.add("wordtag=xor");}
			BuildNcrnaDdlLayoutCounts.main(args.toArray(new String[args.size()]));
		}
		final LineParser1 a=new LineParser1('\t'), b=new LineParser1('\t');
		try(Input first=new Input(outputs[0]+".study.tsv"); Input second=new Input(outputs[1]+".study.tsv")){
			check(Arrays.equals(first.reader.nextLine(), second.reader.nextLine()), "Both XOR entrypoints retain the same schema");
			a.set(first.reader.nextLine());
			b.set(second.reader.nextLine());
			for(int field=0; field<25; field++){
				if(field!=16 && field!=17){check(a.parseString(field).equals(b.parseString(field)), "Omitted versus explicit XOR changed scientific field="+field);}
			}
			check(a.termEquals("OLD_XOR_BASELINE", 0) && a.termEquals("OLD_XOR_HIGH_BITS_RAW_PREHASH", 18)
				&& first.reader.nextLine()==null && second.reader.nextLine()==null,
				"The public default really executes the restored XOR baseline exactly once");
		}
	}

	/** A literal listed word must be rejected in the actual insertion sink,while
	 * an unlisted word survives with its exact native key;includes Long.MIN_VALUE. */
	private static void blacklistNamespace(){
		final int exponent=DynamicDemiLog.exponentBits();
		DynamicDemiLog.setExponent(4);
		try{
			for(final boolean xor:new boolean[]{true,false}){
				for(final int k:new int[]{7,9,8,10}){
					final boolean shortWord=k==7 || k==8;
					final long literal=xor ? shortWord ? 0x4000000000000000L : 0x8000000000000000L
						: k==7 ? 0x4000L : k==8 ? 0x10000L : k==9 ? 0x40000L : 0x100000L;
					final long tag=BuildNcrnaDdlLayoutCounts.wordTag(k, shortWord, xor);
					check(tag==literal, "Each raw-word form must use its literal requested length namespace");
					final LongHashSet blacklist=new LongHashSet(4);
					blacklist.add(literal);
					final DynamicDemiLog observed=DynamicDemiLog.create(512, k, 12345, 0, false, false);
					final DynamicDemiLog expected=DynamicDemiLog.create(512, k, 12345, 0, false, false);
					check(!BuildNcrnaDdlLayoutCounts.storeWord(observed, tag, blacklist) && observed.filledBuckets()==0,
						"A listed zero raw word must be rejected before either XOR or leading-bit native hashing");
					check(BuildNcrnaDdlLayoutCounts.storeWord(observed, tag|1, blacklist), "A nonlisted raw word must survive the filter");
					expected.hashAndStore(literal|1);
					check(Arrays.equals(observed.maxArray(), expected.maxArray()) && observed.filledBuckets()>0,
						"The surviving word retains exact native hash insertion with no silent namespace bypass");
				}
			}
		}finally{DynamicDemiLog.setExponent(exponent);}
	}

	/** Literal complementary halves pin sums and a weighted mean with unequal
	 * window populations; an absent list remains NA rather than zero percent GC. */
	private static void mergeStudy(final String prefix){
		final String empty=prefix+"_study_empty.tsv";
		try(Output out=new Output(empty)){out.writer.println(BuildNcrnaDdlLayoutCounts.STUDY_HEADER);}
		boolean rejected=false;
		try{BuildNcrnaDdlLayoutCounts.mergeStudy(empty, prefix+"_empty_merge");}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !new File(prefix+"_empty_merge.merged.tsv").exists(),
			"A header-only study is not a measured result and must fail before publication");
		final String input=prefix+"_study_halves.tsv", output=prefix+"_study_merge";
		try(Output out=new Output(input)){
			out.writer.println(BuildNcrnaDdlLayoutCounts.STUDY_HEADER);
			out.writer.println("FALSE_SEEDED_WINDOWS\t0\t64\t2\t5\t3\t1\t12\t1\t2\t2\t7\t1\t1\t8\t1.6\t0\t0\tLEADING_BIT_RAW_PREHASH\tEXPLORATORY_SAVED_WINDOWS_NOT_CALLER\t0\t1\tNA\t0.5\t0.6");
			out.writer.println("FALSE_SEEDED_WINDOWS\t0\t64\t2\t6\t4\t1\t20\t2\t3\t2\t9\t1\t2\t12\t2\t0\t0\tLEADING_BIT_RAW_PREHASH\tEXPLORATORY_SAVED_WINDOWS_NOT_CALLER\t1\t0\tNA\t0.75\t0.6");
		}
		BuildNcrnaDdlLayoutCounts.mergeStudy(input, output);
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(output+".merged.tsv")){
			check(in.reader.nextLine()!=null, "A merged study needs an explicit output schema");
			p.set(in.reader.nextLine());
			check(p.terms()==22 && p.parseLong(4)==11 && p.parseLong(5)==7 && p.parseLong(6)==2
				&& p.parseLong(7)==32 && p.parseLong(8)==3 && p.parseLong(9)==5
				&& p.parseLong(10)==4 && p.parseLong(11)==16 && p.parseLong(12)==2 && p.parseLong(13)==3
				&& p.parseLong(14)==20 && Math.abs(p.parseDouble(15)-20.0/11)<0.00000000051
				&& p.termEquals("NA", 16) && p.termEquals("NA", 18)
				&& p.parseDouble(17)==0.5 && p.parseDouble(19)==0.75 && in.reader.nextLine()==null,
				"Disjoint halves sum rank/top-k counts,weight nonzero means by windows,and preserve each fitted list's GC");
		}
		crossSeedSummary(prefix, output+".merged.tsv");
	}

	/** Equal observed values can improve one seed and worsen another;comparison
	 * must use each matched anchor,not just the mean unfiltered result. */
	private static void crossSeedSummary(final String prefix, final String merged){
		final String roster=prefix+"_seed_roster.tsv", anchors=prefix+"_seed_anchors.tsv";
		try(Output out=new Output(roster)){
			out.writer.println("ddl_seed\tmerged_table");
			for(final int seed:new int[]{1,7}){
				final String dir=prefix+"_crossseed_"+seed;
				check(new File(dir).mkdir(), "Create only the named small-table fixture directory");
				try(Output seedFile=new Output(dir+"/ddl_seed.txt")){seedFile.writer.println(""+seed);}
				try(Input in=new Input(merged); Output copy=new Output(dir+"/table.tsv")){
					for(byte[] text; (text=in.reader.nextLine())!=null;){copy.writer.println(text);}
				}
				out.writer.println(seed+"\t"+new File(dir+"/table.tsv").getAbsolutePath());
			}
		}
		try(Output out=new Output(anchors)){
			out.writer.println(BuildNcrnaDdlLayoutCounts.STUDY_HEADER+"\tddl_seed");
			out.writer.println("OLD_XOR_HIGH_BITS\t0\t0\t2\t11\t7\t2\t40\t2\t5\t4\t20\t1\t2\t22\t2\t0\t0\tOLD_XOR_HIGH_BITS_RAW_PREHASH\tEXPLORATORY_SAVED_WINDOWS_NOT_CALLER\t-1\t-1\tNA\tNA\t0.6\t1");
			out.writer.println("OLD_XOR_HIGH_BITS\t0\t0\t2\t11\t7\t2\t24\t4\t5\t4\t12\t3\t4\t11\t1\t0\t0\tOLD_XOR_HIGH_BITS_RAW_PREHASH\tEXPLORATORY_SAVED_WINDOWS_NOT_CALLER\t-1\t-1\tNA\tNA\t0.6\t7");
		}
		BuildNcrnaDdlLayoutCounts.mergeSeeds(roster, anchors, prefix+"_crossseed");
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(prefix+"_crossseed.seed_summary.tsv")){
			check(in.reader.nextLine()!=null, "Cross-seed results need an explicit metric schema");
			p.set(in.reader.nextLine());
			check(p.terms()==15 && p.termEquals("available_rank_sum", 3) && p.parseInt(4)==2
				&& p.parseDouble(5)==32 && p.parseDouble(6)==32 && p.parseDouble(7)==32
				&& p.parseDouble(8)==32 && p.parseInt(9)==1 && p.parseInt(10)==0 && p.parseInt(11)==1
				&& p.parseDouble(12)==0 && p.parseDouble(13)==-8 && p.parseDouble(14)==8,
				"Matched-seed comparison must preserve one improvement and one loss even when both means tie");
		}
		check(rows(prefix+"_crossseed.seed_summary.tsv")==8, "One setting reports all seven metrics across the same seeds");
		final String comparisons=prefix+"_fit_anchor_roster.tsv";
		try(Output out=new Output(comparisons)){
			out.writer.println("ddl_seed\tcomparison_table");
			for(final int seed:new int[]{1,7}){
				final String path=new File(prefix+"_crossseed_"+seed+"/comparison.tsv").getAbsolutePath();
				fitAnchor(path, seed, "/literal/frozen_cohorts", 2);
				out.writer.println(seed+"\t"+path);
			}
		}
		BuildNcrnaDdlLayoutCounts.mergeSeeds(roster, comparisons, prefix+"_filter_only");
		check(rows(prefix+"_filter_only.seed_summary.tsv")==7,
			"Saved fit anchors expose six measured ranking metrics without inventing a nonzero-window baseline");
		try(Input in=new Input(prefix+"_filter_only.seed_summary.tsv")){
			check(in.reader.nextLine()!=null, "Saved fit anchors retain the paired report schema");
			p.set(in.reader.nextLine());
			check(p.parseDouble(5)==32 && p.parseDouble(8)==72 && p.parseDouble(12)==40
				&& p.parseDouble(13)==32 && p.parseDouble(14)==48 && p.parseInt(9)==2,
				"Filter-only deltas must use the saved8/10anchors64and80,not the old7/9anchors40and24");
		}
		final String bad=new File(prefix+"_crossseed_7/bad_comparison.tsv").getAbsolutePath();
		fitAnchor(bad, 7, "/different/cohorts", 2);
		try(Output out=new Output(prefix+"_bad_fit_roster.tsv")){
			out.writer.println("ddl_seed\tcomparison_table");
			out.writer.println("1\t"+new File(prefix+"_crossseed_1/comparison.tsv").getAbsolutePath());
			out.writer.println("7\t"+bad);
		}
		boolean rejected=false;
		try{BuildNcrnaDdlLayoutCounts.mergeSeeds(roster, prefix+"_bad_fit_roster.tsv", prefix+"_bad_fit_baseline");}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !new File(prefix+"_bad_fit_baseline.seed_summary.tsv").exists(),
			"A changed frozen cohort source cannot enter a filter-only paired report");
		try(Output out=new Output(prefix+"_swapped_fit_roster.tsv")){
			out.writer.println("ddl_seed\tcomparison_table");
			out.writer.println("1\t"+new File(prefix+"_crossseed_7/comparison.tsv").getAbsolutePath());
			out.writer.println("7\t"+new File(prefix+"_crossseed_1/comparison.tsv").getAbsolutePath());
		}
		rejected=false;
		try{BuildNcrnaDdlLayoutCounts.mergeSeeds(roster, prefix+"_swapped_fit_roster.tsv", prefix+"_swapped_fit_baseline");}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !new File(prefix+"_swapped_fit_baseline.seed_summary.tsv").exists(),
			"A saved fit baseline from another DDL seed must fail before paired output");
		final String wrongPopulation=new File(prefix+"_crossseed_7/wrong_population.tsv").getAbsolutePath();
		fitAnchor(wrongPopulation, 7, "/literal/frozen_cohorts", 3);
		try(Output out=new Output(prefix+"_wrong_population_roster.tsv")){
			out.writer.println("ddl_seed\tcomparison_table");
			out.writer.println("1\t"+new File(prefix+"_crossseed_1/comparison.tsv").getAbsolutePath());
			out.writer.println("7\t"+wrongPopulation);
		}
		rejected=false;
		try{BuildNcrnaDdlLayoutCounts.mergeSeeds(roster, prefix+"_wrong_population_roster.tsv", prefix+"_wrong_population_baseline");}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !new File(prefix+"_wrong_population_baseline.seed_summary.tsv").exists(),
			"A changed model population cannot be compared against the accepted filtered rows");
	}

	/** Two literal completed-comparison rows are deliberately unlike old tag-study anchors. */
	private static void fitAnchor(final String path, final int seed, final String cohorts, final int models){
		check(seed==1 || seed==7, "The independent baseline fixture has exactly two matched seeds");
		final String[] fields=new String[44];
		Arrays.fill(fields, "0");
		fields[0]="native_s8_l10_mixed512_pooled_both_w100";
		fields[1]="RF00174";
		fields[2]="8";
		fields[3]="10";
		fields[4]="MIDDLE_RIGHT";
		fields[5]="MIXED";
		fields[6]="512";
		fields[9]="POOL_BOTH";
		fields[11]="100";
		fields[12]="member_total";
		fields[13]="2";
		fields[14]="NATIVE_DDL";
		fields[15]="WRAP";
		fields[26]=seed==1 ? "20" : "24";
		fields[31]=""+models;
		fields[32]="7";
		fields[33]="4";
		fields[34]="5000";
		fields[35]="20261113";
		fields[36]=cohorts;
		fields[37]="QUERY_DDL";
		fields[38]=seed==1 ? "2" : "3";
		fields[39]=seed==1 ? "4" : "5";
		fields[40]=seed==1 ? "1" : "2";
		fields[41]=seed==1 ? "2" : "3";
		fields[42]="5";
		fields[43]=seed==1 ? "64" : "80";
		try(Output out=new Output(path)){
			out.writer.println(prok.SummarizeNcrnaLayoutFits.HEADER);
			out.writer.println(String.join("\t", fields));
		}
	}

	/** The previously built control supplies independent reference keys. Queries
	 * are native objects, never raw-key probes. Literal part boundaries prevent
	 * an accidental 384-bin wrapper, namespace tag, or wrong mask from passing. */
	private static void legacyScores(final BuildNcrnaDdlLayoutCounts.Config c, final String control,
			final NcrnaDdlLayoutCounts manifest, final NcrnaDdlSavedData data){
		final DynamicDemiLog[][] reference=new DynamicDemiLog[2][3];
		final LineParser1 p=new LineParser1('\t');
		final int[] starts={0, 128, 384}, ends={128, 384, 512};
		try(Input input=new Input(control)){
			check(input.reader.nextLine()!=null, "The independent control resource must retain its header");
			for(int model=0; model<2; model++){
				final byte[] row=input.reader.nextLine();
				check(row!=null, "The independent fixture has two original model rows");
				p.set(row);
				check(p.terms()==516 && p.parseInt(0)==model, "Control keys must retain the exact 512-bin legacy geometry");
				final char[] keys=new char[512];
				for(int bin=0; bin<512; bin++){keys[bin]=(char)p.parseInt(bin+4);}
				for(int part=0; part<3; part++){
					reference[model][part]=DynamicDemiLog.fromArray(Arrays.copyOfRange(keys, starts[part], ends[part]),
						model, null, part==0 ? 7 : 10, -1);
				}
			}
			check(input.reader.nextLine()==null, "No extra model may enter the independent reference");
		}
		final DynamicDemiLog[] query={DynamicDemiLog.create(128, 7, 12345, 0, false, false),
			DynamicDemiLog.create(256, 10, 12345, 0, false, false), DynamicDemiLog.create(128, 10, 12345, 0, false, false)};
		final long mask=AbstractKmerTable.makeMiddleMask(10, false, 1);
		final int[] scratch=new int[4];
		boolean matched=false, differed=false;
		try(NcrnaDdlLayoutCounts.Reader reader=new NcrnaDdlLayoutCounts.Reader(c.output+".counts.tsv.gz", manifest, data)){
			reader.next(false);
			final char[][] counts=reader.next(true);
			for(int window=0; window<data.windowCount(); window++){
				for(final DynamicDemiLog ddl:query){ddl.clear();}
				final byte[] bases=data.windowBases(window);
				CardinalityTracker.forEachKmer(bases, 7, 10, false, true, (word, k)->{
					final long value=k==7 ? word : word&mask;
					query[k==7 ? 0 : 1].hashAndStore(value);
					query[2].hashAndStore(value^(1L<<(k==7 ? 62 : 63)));
				});
				for(int model=0; model<2; model++){
					for(int part=0; part<3; part++){
						DynamicDemiLog.compareDetailed(query[part].maxArray(), reference[model][part].maxArray(), scratch);
						check(counts[window][4*model+part]==scratch[1],
							"Native control score differs from independently constructed query and original reference: window="+window+", model="+model+", part="+part);
						matched|=scratch[1]>0;
						differed|=scratch[1]<query[part].filledBuckets();
					}
				}
			}
			for(int layout=2; layout<ROWS.length; layout++){reader.next(false);}
			reader.finish();
		}
		check(matched && differed, "Control parity must exercise both matching and unequal nonempty native keys");
	}

	/** Independent literal A/G references and a byte-by-byte query encoder catch
	 * leading-zero loss,length-tag placement,mask placement and ambiguity leaks. */
	private static void leadingBitScores(final BuildNcrnaDdlLayoutCounts.Config c,
			final NcrnaDdlLayoutCounts manifest, final NcrnaDdlSavedData data){
		final DynamicDemiLog[] refs={DynamicDemiLog.create(512, 9, 12345, 0, false, false),
			DynamicDemiLog.create(512, 9, 12345, 0, false, false)};
		refs[0].hashAndStore(0x4000L);
		refs[0].hashAndStore(0x40000L);
		refs[1].hashAndStore(0x7fffL);
		refs[1].hashAndStore(0x7fcffL);
		final DynamicDemiLog query=DynamicDemiLog.create(512, 9, 12345, 0, false, false);
		final int[] compared=new int[4];
		try(NcrnaDdlLayoutCounts.Reader reader=new NcrnaDdlLayoutCounts.Reader(c.output+".counts.tsv.gz", manifest, data)){
			for(int i=0; i<ROWS.length-1; i++){reader.next(false);}
			final char[][] observed=reader.next(true);
			for(int w=0; w<data.windowCount(); w++){
				query.clear();
				final byte[] bases=data.windowBases(w);
				for(final int k:new int[]{7,9}){
					for(int end=k; end<=bases.length; end++){
						long word=0;
						boolean valid=true;
						for(int pos=end-k; pos<end; pos++){
							final int base=bases[pos]&~32;
							valid&=base=='A' || base=='C' || base=='T' || base=='G';
							word=(word<<2)|((base>>1)&3);
						}
						if(valid){query.hashAndStore((k==9 ? word&~0x300L : word)|(k==7 ? 0x4000L : 0x40000L));}
					}
				}
				for(int m=0; m<2; m++){
					DynamicDemiLog.compareDetailed(query.maxArray(), refs[m].maxArray(), compared);
					check(observed[w][4*m]==compared[1], "Leading-bit scores differ from literal native references at window="+w+", model="+m);
				}
			}
			reader.finish();
		}
	}

	private static BuildNcrnaDdlLayoutCounts.Config config(final RankDualNcrnaWindows.Config saved,
			final String control, final String suffix){
		assert(saved.models==2 && saved.groups==3) : "The producer witness reuses the original two-model/three-input saved fixture";
		final BuildNcrnaDdlLayoutCounts.Config c=new BuildNcrnaDdlLayoutCounts.Config();
		c.panels=saved.panels;
		c.selection=saved.selection;
		c.references=saved.references;
		c.members=saved.members;
		c.ledger=saved.ledger;
		c.assignments=saved.assignments;
		c.census=saved.census;
		c.summary=saved.summary;
		c.output=saved.output+suffix;
		c.recipes=c.output+".recipe.tsv";
		c.models=2;
		c.groups=3;
		c.oldTags=false;
		return c;
	}

	private static void writeRecipes(final String path, final String fault){
		try(Output out=new Output(path)){
			out.writer.println(BuildNcrnaDdlLayoutCounts.RECIPE_HEADER);
			for(int i=0; i<ROWS.length; i++){
				if(i==1 && "missing_control".equals(fault)){continue;}
				String row=ROWS[i];
				if(i==1 && "bad_control".equals(fault)){row=row.replace("384", "256");}
				if(i==1 && "bad_query_control".equals(fault)){row=row.replace("QUERY_DDL", "RAW_KEYS");}
				if(i==2 && "bad_part".equals(fault)){row=row.replace("\t192\t64\t", "\t192\t48\t");}
				if(i==0 && "bad_format".equals(fault)){row=row.replace("\t4\t12\t", "\t5\t10\t");}
				out.writer.println(row);
			}
			if("duplicate".equals(fault)){out.writer.println(ROWS[0]);}
		}
	}

	private static int rows(final String path){
		int rows=0;
		try(Input in=new Input(path)){while(in.reader.nextLine()!=null){rows++;}}
		return rows;
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private static final String CONTROL="legacy_member384_consensus128_both";
	private static final String[] ROWS={
		"mixed256_both256\t7\t10\tMIDDLE_RIGHT\tMIXED\t256\t0\tBOTH\t256\t0\t4\t12\tNATIVE_DDL\tfalse\tQUERY_DDL",
		CONTROL+"\t7\t10\tMIDDLE_RIGHT\tPARTITIONED\t384\t128\tBOTH\t128\t0\t4\t12\tLEGACY_PARTITIONS\ttrue\tQUERY_DDL",
		"part192_both128\t6\t9\tMIDDLE_RIGHT\tPARTITIONED\t192\t64\tBOTH\t128\t0\t4\t12\tNATIVE_DDL\tfalse\tQUERY_DDL",
		"pooled256_long\t7\t10\tMIDDLE_RIGHT\tMIXED\t256\t0\tPOOL_LONG\t0\t5\t4\t12\tNATIVE_DDL\tfalse\tQUERY_DDL",
		"member_only512\t8\t12\tMIDDLE_RIGHT\tMIXED\t512\t0\tNONE\t0\t0\t4\t12\tNATIVE_DDL\tfalse\tQUERY_DDL",
		"mixed128_both128\t7\t10\tMIDDLE_RIGHT\tMIXED\t128\t0\tBOTH\t128\t0\t4\t12\tNATIVE_DDL\tfalse\tQUERY_DDL",
		"part256_both512\t7\t10\tMIDDLE_RIGHT\tPARTITIONED\t256\t128\tBOTH\t512\t0\t4\t12\tNATIVE_DDL\tfalse\tQUERY_DDL",
		"mixed128_long512\t7\t10\tNONE\tMIXED\t128\t0\tLONG_ONLY\t512\t0\t4\t12\tNATIVE_DDL\tfalse\tQUERY_DDL",
		"native_s7_l9_mixed512_pooled_both_w100\t7\t9\tMIDDLE_RIGHT\tMIXED\t512\t0\tPOOL_BOTH\t0\t100\t4\t12\tNATIVE_DDL\tfalse\tQUERY_DDL"
	};
	private static final int[] TERMS={2,3,3,1,1,2,3,2,1};
	private static final int[] REFERENCES={0,1,1,0,0,0,1,0,0};
	private static final int[][] LIMITS={{256,256,0,0},{128,256,128,0},{64,128,128,0},{256,0,0,0},
		{512,0,0,0},{128,128,0,0},{128,128,512,0},{128,512,0,0},{512,0,0,0}};
}
