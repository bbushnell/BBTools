package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import cardinality.DynamicDemiLog;
import consensus.BaseGraph;
import ddl.DDLIndexSmall;
import kmer.AbstractKmerTable;
import idaligner.AlignmentStats;
import map.LongHashSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.BuildNcrnaDdlLayoutCounts;
import prok.BuildNcrnaDdlProfiles;
import prok.NcrnaClusterAssignment;
import prok.NcrnaClusterAssignment.Policy;
import prok.NcrnaDdlClusterer;
import prok.NcrnaDdlConsensus;
import prok.NcrnaDdlSketchBuilder;
import prok.NcrnaDdlLayoutVotes;
import prok.NcrnaDdlProfileResource;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.TrnaConsensusBuilder;
import shared.Shared;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

/** Finite policy boundaries, independent word enumeration, Quantum orientations,
 * real-member graph depth and full builder/loader/assignment smoke.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlClustererTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaDdlClustererTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, split[0], split.length>1 ? split[1] : null), "Unknown fixture argument: "+arg);
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append, "Supply out=fresh_native_fixture_directory");
			final Path root=Files.createDirectory(new File(parser.out1).toPath());
			DynamicDemiLog.setExponent(4);
			policies();
			kernels();
			onlineProfiles();
			consensus();
			clustering();
			builder(root);
			System.err.println("NCRNA_DDL_CLUSTERER_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void policies(){
		final byte[][] models=new byte[300][];
		for(int i=0; i<models.length; i++){models[i]=new byte[]{(byte)(i>>>8), (byte)i};}
		final int[] counts=new int[models.length];
		final float[] scores=new float[models.length];
		final NcrnaClusterAssignment.Alignment alignment=(query, reference)->scores[((reference[0]&255)<<8)|(reference[1]&255)];
		final byte[] query={'A'};
		scores[0]=0.8f;
		scores[60]=0.99f;
		NcrnaClusterAssignment match=new NcrnaClusterAssignment(Policy.TOP50, 100, 512);
		check(match.bestFamily(query, models, counts, 0.7f, alignment)==0 && match.alignments==50,
			"A passing top50 winner must deliberately outrank a better alignment outside the window");
		match=new NcrnaClusterAssignment(Policy.ALL, 100, 512);
		check(match.bestFamily(query, models, counts, 0.7f, alignment)==60 && match.alignments==300,
			"The all-consensus comparator must find the exhaustive best identity");
		scores[0]=0;
		match=new NcrnaClusterAssignment(Policy.TOP50, 100, 512);
		check(match.bestFamily(query, models, counts, 0.7f, alignment)==60 && match.alignments==300,
			"An all-top50-fail query with zero DDL hits must exhaust once and find the later passing model");
		Arrays.fill(scores, 0);
		scores[10]=0.8f;
		scores[80]=0.9f;
		scores[181]=0.99f;
		match=new NcrnaClusterAssignment(Policy.PATIENCE, 100, 512);
		check(match.bestFamily(query, models, counts, 0.7f, alignment)==80 && match.alignments==181,
			"Patience starts at first pass and resets only on improvement; the181st-index model is beyond100further trials");
		match=new NcrnaClusterAssignment(Policy.PATIENCE, 0, 512);
		check(match.bestFamily(query, models, counts, 0.7f, alignment)==10 && match.alignments==11,
			"Patience0 stops at the first passing model, never before any pass");
		Arrays.fill(scores, 0);
		match=new NcrnaClusterAssignment(Policy.PATIENCE, 100, 512);
		check(match.bestFamily(query, models, counts, 0.7f, alignment)<0 && match.alignments==300,
			"No passing identity must exhaust the entire roster regardless of patience");
		scores[299]=0.85f;
		match=new NcrnaClusterAssignment(Policy.PATIENCE, 100, 512);
		check(match.bestFamily(query, models, counts, 0.7f, alignment)==299 && match.alignments==300,
			"A first pass at the end of the roster remains assignable");
		Arrays.fill(scores, 0);
		counts[8]=3;
		scores[5]=scores[8]=0.8f;
		match=new NcrnaClusterAssignment(Policy.PATIENCE, 100, 512);
		check(match.bestFamily(query, models, counts, 0.7f, alignment)==5,
			"Equal passing identities use stable family IDs even when DDL order visits the larger ID first");
	}

	private static void kernels(){
		final byte[] sequence=bases("AAAAAAAAAAAAACGTTGCAGTCAGTACGATCGATGCTAGCTAGTNNACGTGCTAGCGATCGTACGATCGTAGCTAGCTAGTCGATGC");
		final LongHashSet blacklist=new LongHashSet(4);
		blacklist.add(BuildNcrnaDdlLayoutCounts.wordTag(8, true, true));
		for(final int bins:new int[]{512,1024}){
			for(final long seed:new long[]{1,12345}){
				final NcrnaDdlSketchBuilder builder=new NcrnaDdlSketchBuilder(bins, seed, blacklist);
				final DynamicDemiLog actual=builder.newSketch(), expected=builder.newSketch();
				builder.fill(sequence, actual);
				slow(sequence, expected, bins, seed, blacklist);
				check(Arrays.equals(actual.maxArray(), expected.maxArray()), "Quarter-routing keys differ from independent forward-window enumeration");
				final char[] fixedKeys=new char[bins];
				builder.fillKeys(sequence, fixedKeys);
				check(Arrays.equals(fixedKeys, expected.maxArray()), "Fixed native4+12 key arrays must equal actual DDL objects and independent enumeration");
				final int previous=DynamicDemiLog.exponentBits();
				try{
					DynamicDemiLog.setExponent(5);
					NcrnaDdlSketchBuilder.forKeyArrays(bins, seed, blacklist).fillKeys(sequence, fixedKeys);
					check(DynamicDemiLog.exponentBits()==5 && Arrays.equals(fixedKeys, expected.maxArray()),
						"Caller key arrays must retain the selected4+12 representation without changing the caller's global exponent5");
				}finally{DynamicDemiLog.setExponent(previous);}
				builder.fill(bases("NNNNNNNNNNNN"), actual);
				builder.fillKeys(bases("NNNNNNNNNNNN"), fixedKeys);
				for(final char value:actual.maxArray()){check(value==0, "Ambiguous spans and reused sketches must not retain old keys");}
				for(final char value:fixedKeys){check(value==0, "Reused fixed query keys must clear the preceding sequence");}
			}
		}
	}

	private static void slow(final byte[] sequence, final DynamicDemiLog output,
			final int bins, final long seed, final LongHashSet blacklist){
		for(final int k:new int[]{8,10}){
			for(int start=0; start+k<=sequence.length; start++){
				long word=0;
				boolean valid=true;
				for(int pos=start; pos<start+k; pos++){
					final int value=sequence[pos]=='A' ? 0 : sequence[pos]=='C' ? 1 : sequence[pos]=='T' ? 2 : sequence[pos]=='G' ? 3 : -1;
					if(value<0){valid=false;break;}
					word=(word<<2)|value;
				}
				if(!valid){continue;}
				if(k==10){word&=AbstractKmerTable.makeMiddleMask(10, false, 1);}
				final long tagged=word|BuildNcrnaDdlLayoutCounts.wordTag(k, k==8, true);
				if(blacklist.contains(tagged)){continue;}
				final long hash=Tools.hash64shift(tagged^seed);
				output.storeHashed(tagged, hash, BuildNcrnaDdlLayoutCounts.partitionFor(hash, bins, 2, k==8), 2);
			}
		}
	}

	private static void onlineProfiles(){
		final NcrnaDdlSketchBuilder builder=new NcrnaDdlSketchBuilder(512, 12345, new LongHashSet(4));
		final byte[] sequence=randomBases(500, 1701);
		final NcrnaDdlSketchBuilder.Prefix prefix=new NcrnaDdlSketchBuilder.Prefix(512);
		builder.fillPrefix(sequence, prefix);
		int bin=-1, low=0, high=0;
		for(int i=0; i<prefix.keyCount() && bin<0; i++){
			final int a=prefix.packedKey(i);
			for(int j=i+1; j<prefix.keyCount(); j++){
				final int b=prefix.packedKey(j);
				if(a>>>16==b>>>16 && (a&65535)!=(b&65535)){
					bin=a>>>16;
					low=Math.min(a&65535, b&65535);
					high=Math.max(a&65535, b&65535);
					break;
				}
			}
		}
		check(bin>=0, "The fixed consensus must exercise two distinct eligible keys in one bin");
		final char[] lower=new char[512], upper=new char[512];
		lower[bin]=(char)low;
		upper[bin]=(char)high;
		final NcrnaDdlLayoutVotes votes=new NcrnaDdlLayoutVotes(512);
		votes.add(lower);
		check(votes.finish(prefix, 100, true).keys[bin]==low,
			"A lower nonmaximum consensus key must receive its bonus and win by real-member support");
		votes.reset();
		final DynamicDemiLog consensus=builder.newSketch();
		builder.fill(sequence, consensus);
		final char[] selected=consensus.maxArray().clone();
		votes.addAndSelect(upper, prefix, 100, selected);
		check(Arrays.equals(selected, votes.finish(prefix, 100, true).keys), "First incremental member differs from full voting");
		votes.addAndSelect(lower, prefix, 100, selected);
		check(Arrays.equals(selected, votes.finish(prefix, 100, true).keys), "Second incremental member differs from full voting");
		votes.addAndSelect(lower, prefix, 100, selected);
		check(selected[bin]==low && Arrays.equals(selected, votes.finish(prefix, 100, true).keys),
			"The third member must update the current profile before a fourth query, not wait for a size threshold");
		final DynamicDemiLog old=DynamicDemiLog.fromArray(upper, 0, null, 10, -1);
		final DynamicDemiLog current=DynamicDemiLog.fromArray(lower, 0, null, 10, -1);
		check(current.filledBuckets()==1, "Native publication must recompute occupancy metadata");
		final DDLIndexSmall frozen=new DDLIndexSmall(new DynamicDemiLog[]{old});
		final DynamicDemiLog[] changed={current};
		final int[] counts=new int[1];
		frozen.queryReplacing(current, counts, changed, new int[]{0}, 1);
		check(counts[0]==1, "A changed row's current match must replace its frozen missing score");
		frozen.queryReplacing(old, counts, changed, new int[]{0,0}, 2);
		check(counts[0]==0, "Old matches must disappear; duplicate changed IDs must not add stale or repeated evidence");
	}

	private static void consensus(){
		final AlignmentStats direct=new AlignmentStats(true);
		direct.matchString=bases("mDIImm");
		direct.rStart=0;
		direct.rStop=3;
		final Read normalized=NcrnaDdlConsensus.toModelFrame(new Read(bases("ACGTA"), null, "direct", 0), direct, true);
		check(Arrays.equals(normalized.match, bases("mSImm")) && normalized.start==0 && normalized.stop==3,
			"Direct Quantum opposing indels must normalize without changing model/member consumption");
		final AlignmentStats swapped=new AlignmentStats(true);
		swapped.matchString=bases("mIDmmmm");
		swapped.rStart=1;
		swapped.rStop=6;
		final Read inverted=NcrnaDdlConsensus.toModelFrame(new Read(bases("AACGTAGG"), null, "swapped", 1), swapped, false);
		check(Arrays.equals(inverted.match, bases("mSmmmm")) && inverted.start==0 && inverted.stop==5 && inverted.length()==6,
			"Swapped Quantum opposing indels must invert into the actual model coordinate frame");
		final byte[] reference=randomBases(160, 11);
		final ArrayList<Read> identical=new ArrayList<Read>();
		for(int i=0; i<3; i++){identical.add(new Read(reference.clone(), null, "identical"+i, i));}
		final NcrnaDdlConsensus.Fit same=NcrnaDdlConsensus.refine(identical, reference, 0.9f, 2, false);
		final BaseGraph graph=same.graph("unit_votes");
		check(same.alignments==3 && same.identities.length==3, "Stable identical consensus needs exactly three actual placements and no reporting realignment");
		for(final float identity:same.identities){check(identity==1, "Reported member identity must be the retained final Quantum identity");}
		check(same.members()==3 && Arrays.equals(graph.original, same.bases), "Final consensus and HBM must share actual coordinates");
		for(int pos=0; pos<reference.length; pos++){
			check(graph.ref[pos].countSum==3 && graph.del[pos].countSum==0,
				"An external pivot must not become a fourth synthetic graph voter");
		}
		final ArrayList<Read> varied=new ArrayList<Read>();
		varied.add(new Read(reference.clone(), null, "reference", 0));
		varied.add(new Read(Arrays.copyOfRange(reference, 3, reference.length-2), null, "shorter", 1));
		final ByteBuilder longer=new ByteBuilder();
		longer.append("TT").append(reference).append("GG");
		varied.add(new Read(longer.toBytes(), null, "longer", 2));
		final ByteBuilder insertion=new ByteBuilder();
		insertion.append(Arrays.copyOfRange(reference, 0, 80)).append("ACGT").append(Arrays.copyOfRange(reference, 80, reference.length));
		varied.add(new Read(insertion.toBytes(), null, "insertion", 3));
		final NcrnaDdlConsensus.Fit fit=NcrnaDdlConsensus.refine(varied, reference, 0.8f, 2, false);
		check(fit.members()==4 && Arrays.equals(fit.bases, fit.graph("both_orientations").original),
			"Both Quantum operand orientations and indels must retain every member and the final model frame");
	}

	private static void clustering(){
		for(final Policy policy:Policy.values()){
			for(final int threads:new int[]{1,2}){
				final ArrayList<Read> reads=population();
				final NcrnaDdlClusterer clusterer=new NcrnaDdlClusterer(reads, policy, 100, 512, 12345,
					new LongHashSet(4), 0.95f, 0.95f, 2, 1, false, true, threads);
				final ArrayList<ArrayList<Read>> clusters=clusterer.cluster(null);
				check(clusterer.alignments==26 && clusterer.consensusAlignments==16,
					"Assignment and consensus work must count their distinct actual alignment calls for this fixed population");
				if(policy==Policy.ALL){
					check(clusterer.indexRebuilds==0 && clusterer.profileRebuilds==0 && clusterer.profileUpdates==0,
						"The brute-force baseline must not construct or charge unused DDL profiles/indexes");
				}
				check(clusters.size()==2 && clusters.get(0).size()==4 && clusters.get(1).size()==4,
					"A newly published family must immediately recruit exact later members, with stable frozen-round membership");
				for(int family=0; family<clusters.size(); family++){
					for(final Read read:clusters.get(family)){
						check(Arrays.equals(read.bases, clusterer.consensus(family)), "Frozen worker count changed the assigned family");
					}
				}
			}
		}
	}

	private static void builder(final Path root) throws Exception{
		final Path input=root.resolve("members.fa"), blacklist=root.resolve("padding.tsv");
		try(Output output=new Output(input.toString())){
			final ByteBuilder row=new ByteBuilder();
			for(final Read read:population()){
				row.clear().append('>').append(read.id).append("_tid_1\n").append(read.bases).nl();
				output.writer.print(row);
			}
		}
		try(Output output=new Output(blacklist.toString())){
			output.writer.println(NcrnaDdlSketchBuilder.BLACKLIST_HEADER);
			final ByteBuilder row=new ByteBuilder();
			for(final int k:new int[]{8,10}){
				final int limit=k==8 ? 256 : 2048;
				int rank=0;
				for(long raw=0; rank<limit; raw++){
					if(k==10 && (raw&AbstractKmerTable.makeMiddleMask(10, false, 1))!=raw){continue;}
					final long tagged=raw|BuildNcrnaDdlLayoutCounts.wordTag(k, k==8, true);
					row.clear().append("TRUE_WINDOW_PADDING\t").append(k==8 ? "SHORT8" : "MASKED10").tab()
						.append(++rank).tab().append(raw).tab().append(tagged).append("\t1\t").append(word(raw, k)).append("\tfixture\n");
					output.writer.print(row);
				}
			}
		}
		check(NcrnaDdlSketchBuilder.readBlacklist(blacklist.toString(), "fixture").size()==2304,
			"The existing padding-table reader must retain exactly the declared prefixes");
		for(final String policy:new String[]{"top50","patience","all"}){
			final String prefix=root.resolve(policy).toString();
			TrnaConsensusBuilder.main(new String[]{"in="+input, "out="+prefix+".fa", "outmodel="+prefix+".hbmc",
				"outassignments="+prefix+".assignments.tsv", "family=fixture", "clusterpolicy="+policy,
				"clusterpatience=100", "ddlbins=512", "paddingblacklist="+blacklist, "paddingfit=fixture",
				"clusterid=0.95", "recruitid=0.95", "mincluster=1", "mingroup=1", "reassignrounds=2", "passes=1", "endtrim=0", "t=2", "ow=f"});
			final byte[][] consensuses=TrnaConsensusBuilder.loadLibrary(prefix+".fa");
			final BaseGraph[] models=TrnaConsensusBuilder.loadModels(prefix+".hbmc");
			check(consensuses.length==2 && models.length==2, "The public builder must emit both actual families");
			for(int i=0; i<models.length; i++){
				check(Arrays.equals(consensuses[i], models[i].original), "Serialized HBM and consensus coordinates diverged");
			}
			int members=0;
			try(Input assignments=new Input(prefix+".assignments.tsv")){
				check(assignments.reader.nextLine()!=null, "Assignment output lacks its header");
				final LineParser1 p=new LineParser1('\t');
				for(byte[] row=assignments.reader.nextLine(); row!=null; row=assignments.reader.nextLine()){
					p.set(row);
					check(p.terms()==6 && p.termEquals("kept", 5), "Every finite fixture member must retain its explicit kept assignment");
					members++;
				}
			}
			check(members==8, "The end-to-end builder lost or duplicated a fixture member");
		}
		final Path generic=root.resolve("generic_members.fa");
		try(Output output=new Output(generic.toString())){
			final ByteBuilder row=new ByteBuilder();
			for(final Read read:population()){
				final int id=(int)read.numericID;
				output.writer.print(row.clear().append('>').append(read.id).append(" cmfrom=").append(id==0 || id==2 ? 45 : 1)
					.append(" cmto=").append(id==1 || id==2 ? 169 : 189).append(" cmlen=189").nl().append(read.bases).nl());
			}
		}
		for(final String policy:new String[]{"top50","patience"}){
			final String prefix=root.resolve("unfiltered_"+policy).toString();
			TrnaConsensusBuilder.main(new String[]{"in="+generic, "out="+prefix+".fa", "outmodel="+prefix+".hbmc",
				"outassignments="+prefix+".assignments.tsv", "outbuildstats="+prefix+".stats.tsv", "family=RF00174", "consensusprefix=cobalamin_RF00174_",
				"clusterpolicy="+policy, "ddlbins=1024", "paddingblacklist=none", "clusterid=0.95", "recruitid=0.95",
				"mincluster=1", "mingroup=1", "reassignrounds=2", "passes=1", "endtrim=0", "t=2", "ow=f"});
			final byte[][] consensuses=TrnaConsensusBuilder.loadLibrary(prefix+".fa");
			final BaseGraph[] models=TrnaConsensusBuilder.loadModels(prefix+".hbmc");
			check(consensuses.length==2 && models.length==2, "Explicit unfiltered1024-bin clustering must retain both nucleotide families");
			for(int i=0; i<models.length; i++){
				check(Arrays.equals(consensuses[i], models[i].original), "Unfiltered initial models must retain actual consensus/HBM coordinates");
			}
			check(statistic(prefix+".stats.tsv", "total_executed_alignments", "all")==50
				&& statistic(prefix+".stats.tsv", "assignment_alignments", "all")==26
				&& statistic(prefix+".stats.tsv", "clustering_consensus_alignments", "all")==16
				&& statistic(prefix+".stats.tsv", "final_consensus_alignments", "all")==8,
				"Builder reporting must count all actual assignment and trace placements once, including final HBM placements");
			check(statistic(prefix+".stats.tsv", "cluster_size", "4")==2
				&& statistic(prefix+".stats.tsv", "identity_nearest_thousandth", "1000")==8,
				"Final model-size and retained-identity histograms must conserve both clusters and all eight members");
			final String firstLabel="cobalamin_RF00174_unknown_c0 n=4";
			check(statistic(prefix+".stats.tsv", "cm_from_min", firstLabel)==1
				&& statistic(prefix+".stats.tsv", "cm_from_max", firstLabel)==45
				&& statistic(prefix+".stats.tsv", "cm_to_min", firstLabel)==169
				&& statistic(prefix+".stats.tsv", "cm_to_max", firstLabel)==189
				&& statistic(prefix+".stats.tsv", "cm_left44_only", firstLabel)==1
				&& statistic(prefix+".stats.tsv", "cm_right20_only", firstLabel)==1
				&& statistic(prefix+".stats.tsv", "cm_left44_and_right20", firstLabel)==1,
				"Per-model native-end spreads and joint end classes must preserve overlap instead of adding marginal populations");
			int genericMembers=0;
			try(Input assignments=new Input(prefix+".assignments.tsv")){
				assignments.reader.nextLine();
				final LineParser1 p=new LineParser1('\t');
				for(byte[] row=assignments.reader.nextLine(); row!=null; row=assignments.reader.nextLine()){
					p.set(row);
					check(p.terms()==6 && p.parseInt(1)==-1 && p.termEquals("kept", 5),
						"Generic nucleotide members must retain unknown taxonomy explicitly, without an invented source taxid");
					genericMembers++;
				}
			}
			check(genericMembers==8, "Generic member export must conserve every actual original record");
			selectedProfiles(prefix, generic, consensuses, blacklist);
		}
		orphanProfile(root, generic);
		for(final boolean suppliedNone:new boolean[]{false,true}){
			final Path rejectedOutput=root.resolve(suppliedNone ? "none_with_scope.fa" : "missing_blacklist.fa");
			final ArrayList<String> flags=new ArrayList<String>(Arrays.asList("in="+input, "out="+rejectedOutput,
				"clusterpolicy=top50", "ddlbins=1024", "t=1", "ow=f"));
			if(suppliedNone){
				flags.add("paddingblacklist=none");
				flags.add("paddingfit=fixture");
			}
			boolean rejected=false;
			try{new TrnaConsensusBuilder(flags.toArray(new String[0]));}
			catch(IllegalArgumentException expected){
				rejected=expected.getMessage().contains(suppliedNone ? "no learned paddingfit scope" : "padding word file");
			}
			check(rejected && !Files.exists(rejectedOutput),
				"Unfiltered clustering must be selected explicitly and cannot claim a learned fit scope or write a rejected result");
		}
		final Path singleton=root.resolve("singleton.fa");
		try(Output output=new Output(singleton.toString())){
			output.writer.print(new ByteBuilder().append(">singleton_tid_1\n").append(randomBases(160, 11)).nl());
		}
		for(final String policy:new String[]{"top50","patience","all"}){
			final String prefix=root.resolve("singleton_"+policy).toString();
			final ArrayList<String> flags=new ArrayList<String>(Arrays.asList("in="+singleton, "out="+prefix+".fa",
				"outmodel="+prefix+".hbmc", "clusterpolicy="+policy, "mincluster=1", "mingroup=1", "passes=1", "reassignrounds=1", "endtrim=0", "t=1", "ow=f"));
			if(!policy.equals("all")){
				flags.add("ddlbins=512");
				flags.add("paddingblacklist="+blacklist);
				flags.add("paddingfit=fixture");
			}
			TrnaConsensusBuilder.main(flags.toArray(new String[0]));
			final byte[][] sequences=TrnaConsensusBuilder.loadLibrary(prefix+".fa");
			final BaseGraph[] models=TrnaConsensusBuilder.loadModels(prefix+".hbmc");
			check(sequences.length==1 && models.length==1 && Arrays.equals(sequences[0], models[0].original),
				"Every selected new policy must build an actual singleton in its own model frame");
			for(int pos=0; pos<models[0].original.length; pos++){
				check(models[0].ref[pos].countSum+models[0].del[pos].countSum==1,
					"A singleton model must contain one real member observation");
			}
		}
		final Path malformed=root.resolve("bad_padding.tsv");
		Files.write(malformed, bases("wrong_header\n"));
		boolean rejected=false;
		try{NcrnaDdlSketchBuilder.readBlacklist(malformed.toString(), "fixture");}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Malformed padding metadata must fail loudly");
	}

	private static void selectedProfiles(final String prefix, final Path input, final byte[][] references, final Path blacklist) throws Exception{
		final String[] names=TrnaConsensusBuilder.lastLoadedNames.clone();
		for(final int bins:new int[]{512,1024}){
			for(final boolean filtered:new boolean[]{false,true}){
				final String output=prefix+"_selected_"+bins+"_"+filtered+".tsv";
				final ArrayList<String> flags=new ArrayList<String>(Arrays.asList("selected=t", "in="+input,
					"assignments="+prefix+".assignments.tsv", "ref="+prefix+".fa", "out="+output,
					"expectedinput=8", "expectedmembers=8", "expectedmodels=2", "partitionbins="+bins,
					"paddingblacklist="+(filtered ? blacklist : "none")));
				if(filtered){flags.add("paddingfit=fixture");}
				final LongHashSet words=filtered ? NcrnaDdlSketchBuilder.readBlacklist(blacklist.toString(), "fixture") : new LongHashSet(4);
				final int previous=DynamicDemiLog.exponentBits();
				try{
					DynamicDemiLog.setExponent(5);
					BuildNcrnaDdlProfiles.main(flags.toArray(new String[0]));
					final NcrnaDdlProfileResource resource=NcrnaDdlProfileResource.load(output, names, references);
					check(resource.selected() && resource.models()==2 && resource.buckets()==bins,
						"Selected export must bind the actual new library and its declared native layout");
					final NcrnaDdlSketchBuilder builder=NcrnaDdlSketchBuilder.forKeyArrays(bins, 12345, words);
					final char[][] expected=new char[2][bins];
					for(int model=0; model<2; model++){
						builder.fillKeys(references[model], expected[model]);
						check(resource.members(model)==4, "Virtual consensus bonuses must never become actual member counts");
						for(int bin=0; bin<bins; bin++){
							check(resource.key(model, bin)==expected[model][bin], "Identical real members plus all-key consensus support must preserve their exact native maxima");
						}
					}
					final NcrnaDdlProfileResource.Scorer scorer=resource.newScorer();
					final float[] scores=new float[2];
					final char[] query=new char[bins];
					final ByteBuilder padded=new ByteBuilder();
					padded.append("AAAAAAAAAAAA").append(references[1]);
					for(final byte[] sequence:new byte[][]{references[0], padded.toBytes(), bases("NNNNNNNNNNNN")}){
						builder.fillKeys(sequence, query);
						scorer.score(sequence, scores);
						for(int model=0; model<2; model++){
							int matches=0;
							for(int bin=0; bin<bins; bin++){if(query[bin]>0 && query[bin]==expected[model][bin]){matches++;}}
							check(scores[model]==matches, "Selected query scoring must count equal nonempty bin maxima, including the retained padding filter");
						}
					}
					check(DynamicDemiLog.exponentBits()==5, "Profile construction/load/query must restore or preserve the surrounding caller's global exponent");
				}finally{DynamicDemiLog.setExponent(previous);}
			}
		}
		final Path duplicate=new File(prefix+"_duplicate_members.fa").toPath();
		final byte[] original=Files.readAllBytes(input);
		try(Output out=new Output(duplicate.toString())){out.writer.print(new ByteBuilder().append(original).append(original));}
		final Path rejected=new File(prefix+"_duplicate_profiles.tsv").toPath();
		boolean failed=false;
		try{
			BuildNcrnaDdlProfiles.main(new String[]{"selected=t", "in="+duplicate, "assignments="+prefix+".assignments.tsv",
				"ref="+prefix+".fa", "out="+rejected, "expectedinput=8", "expectedmembers=8", "expectedmodels=2",
				"partitionbins=1024", "paddingblacklist=none"});
		}catch(IllegalArgumentException expected){failed=expected.getMessage().contains("exactly once");}
		check(failed && !Files.exists(rejected), "Duplicated FASTA genes must fail before a resource is published or receive repeated votes");
	}

	private static long statistic(final String path, final String kind, final String value){
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(path)){
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				if(p.termEquals(kind, 0) && p.termEquals(value, 1)){return p.parseLong(2);}
			}
		}
		throw new AssertionError("Missing build statistic: "+kind+" value="+value);
	}

	private static void orphanProfile(final Path root, final Path input) throws Exception{
		final Path seed=root.resolve("one_family_seed.fa");
		try(Output output=new Output(seed.toString())){
			output.writer.print(new ByteBuilder().append(">one_family\n").append(population().get(0).bases).nl());
		}
		final String prefix=root.resolve("seeded_orphans").toString();
		TrnaConsensusBuilder.main(new String[]{"in="+input, "seedin="+seed, "out="+prefix+".fa",
			"outassignments="+prefix+".assignments.tsv", "family=RF00174", "consensusprefix=cobalamin_RF00174_",
			"clusterpolicy=all", "clusterid=0.95", "recruitid=0.95", "mincluster=1", "mingroup=1",
			"reassignrounds=1", "passes=1", "endtrim=0", "t=1", "ow=f"});
		final byte[][] reference=TrnaConsensusBuilder.loadLibrary(prefix+".fa");
		final String[] names=TrnaConsensusBuilder.lastLoadedNames.clone();
		check(reference.length==1, "Fixed-seed construction must leave the other genuine family explicitly unmodeled");
		final String profiles=prefix+".profiles.tsv";
		BuildNcrnaDdlProfiles.main(new String[]{"selected=t", "in="+input, "ref="+prefix+".fa",
			"assignments="+prefix+".assignments.tsv", "out="+profiles, "expectedinput=8", "expectedmembers=4",
			"expectedmodels=1", "partitionbins=1024", "paddingblacklist=none"});
		final NcrnaDdlProfileResource resource=NcrnaDdlProfileResource.load(profiles, names, reference);
		check(resource.models()==1 && resource.members(0)==4, "All eight original rows must be checked while only four kept genes vote; orphans never acquire invented memberships");
		final byte[][] changed={reference[0].clone()};
		changed[0][0]=changed[0][0]=='A' ? (byte)'C' : (byte)'A';
		boolean rejected=false;
		try{NcrnaDdlProfileResource.load(profiles, names, changed);}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("content differs");}
		check(rejected, "A selected resource must bind actual consensus contents directly, without relying on filenames or hashes");
	}
	private static ArrayList<Read> population(){
		final ArrayList<Read> reads=new ArrayList<Read>();
		for(int family=0; family<2; family++){
			final byte[] sequence=randomBases(160, 11+family*6);
			for(int i=0; i<4; i++){reads.add(new Read(sequence.clone(), null, "sequence"+reads.size(), reads.size()));}
		}
		return reads;
	}
	private static byte[] randomBases(final int length, final long seed){
		final byte[] alphabet=bases("ACGT"), result=new byte[length];
		final Random random=new Random(seed);
		for(int i=0; i<length; i++){result[i]=alphabet[random.nextInt(4)];}
		return result;
	}
	private static String word(long raw, final int k){
		final char[] sequence=new char[k], alphabet={'A','C','T','G'};
		for(int i=k-1; i>=0; i--){sequence[i]=alphabet[(int)(raw&3)];raw>>>=2;}
		return new String(sequence);
	}

	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	private static byte[] bases(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
}
