package test.prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.ArrayList;

import parse.PreParser;
import parse.Parse;
import parse.Parser;
import prok.NcrnaReferenceCalls;
import prok.NcrnaFamily;
import prok.NcrnaShortlistTarget;
import prok.NcrnaScavenger;
import prok.GeneModel;
import prok.NcrnaCandidateReuse;
import prok.NcrnaJointPolicy;
import prok.NcrnaPatience;
import prok.NcrnaCpuProfile;
import prok.NcrnaReportedOrf;
import prok.NcrnaReferenceCalls.Call;
import prok.NcrnaReferenceCalls.Result;
import prok.NcrnaReferenceCalls.Rule;
import prok.NcrnaTextIO.Output;
import prok.Orf;
import prok.ProkObject;
import prok.TrnaKmerIndex;
import shared.Shared;
import structures.ByteBuilder;

/** Finite cutoff semantics and live/saved projection parity, with literal oracles.
 * Whole-caller equivalence remains a separate real-input check.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaPolicyIterationTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaPolicyIterationTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(java.util.Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown test parameter: "+arg);
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append, "Supply a fresh out=fixture.gff path with ow=f and append=f");
			shortlist();
			patience();
			targetSelection();
			candidateState();
			groupMembership();
			jointTable(parser.out1+".policies.tsv");
			projection(parser.out1);
		}finally{Shared.closeStream(pp.outstream);}
		System.out.println("NCRNA_POLICY_ITERATION_TEST_PASS");
	}

	private static void shortlist(){
		final byte[][] library={bytes("ACGT"), bytes("AA"), bytes("CC")};
		final TrnaKmerIndex index=new TrnaKmerIndex(library, 1, true, 0, 0, 0, 1);
		final byte[] query=bytes("ACGT");
		check(Arrays.equals(index.shortlist(query, 3), new int[]{0, 1, 2}), "Counts 4,2,2 rank by count then model ID");
		index.setAdaptiveFloor(3);
		check(Arrays.equals(index.shortlist(query, 3), new int[]{0}), "Changing the actual floor removes both count-two models");
		index.setAdaptiveFloor(2);
		check(Arrays.equals(index.shortlist(query, 3), new int[]{0, 1, 2}), "Equality to the adaptive floor is inclusive");
		check(Arrays.equals(index.shortlist(query, 2), new int[]{0, 1}), "Cardinality cap preserves the first of tied models");
		rejected(() -> index.setFixedMinHits(2), "An adaptive index cannot silently acquire a fixed cutoff");
		index.setAdaptiveFloor(5);
		check(Arrays.equals(index.shortlist(query, 3), new int[]{0}), "Multi-model keep-best survives above-best floor");
		final TrnaKmerIndex single=new TrnaKmerIndex(new byte[][]{bytes("ACGT")}, 1, true, 4, 0, 0, 1);
		check(single.shortlist(query, 1).length==1, "Single-model floor equality passes");
		single.setAdaptiveFloor(5);
		check(single.shortlist(query, 1).length==0, "Single-model prefilter has no keep-best exception");
		try{
			index.setAdaptiveFloor(Float.NaN);
			throw new AssertionError("NaN floor accepted");
		}
		catch(IllegalArgumentException expected){}
	}

	/** Literal attempt streams distinguish latching, raw identity and phase state. */
	private static void patience(){
		final NcrnaPatience state=new NcrnaPatience();
		state.reset(1);
		check(!state.observe(.9f, false, true) && !state.observe(.8f, false, true),
			"A stagnating identity below its model cutoff cannot stop before any pass");
		check(state.observe(.8f, true, true),
			"A lower-identity first pass under unequal cutoffs must retain pre-pass stagnation, not reset to the new winning model");
		state.reset(2);
		check(!state.observe(.8f, true, true) && !state.observe(.8f, true, true),
			"A tie consumes one attempt without resetting strict raw-identity improvement");
		check(state.observe(.8f, true, true), "Two equal identities must exhaust patience two rather than reset it");
		state.reset(2);
		check(!state.observe(.8f, true, true) && !state.observe(.8f, true, true), "Set one used attempt before the span rejection");
		check(!state.observe(.9f, false, false) && !state.observe(.85f, true, true),
			"Invalid spans consume nothing and cannot block a later eligible identity improvement");
		check(!state.observe(.7f, false, true) && state.observe(.6f, false, true),
			"After a pass, two eligible failures can stop without a current passing alignment");
		state.reset(0);
		check(!state.observe(.9f, true, false) && !state.observe(.5f, false, true)
			&& state.observe(.8f, true, true), "Zero stops only on the first span-eligible pass");
		state.reset(Integer.MAX_VALUE);
		check(!state.observe(1f, true, true) && !state.observe(.5f, false, true),
			"The whole-shortlist value cannot stop after a latched pass");
		state.reset(1);
		check(!state.observe(.5f, false, true) && !state.observe(.5f, false, true)
			&& !state.observe(.6f, true, true), "A new phase must forget the prior latch and best raw identity");
		rejected(() -> state.reset(-2), "Values below the explicit legacy sentinel must fail");
		for(final int limit:new int[]{-1, 0, 1, Integer.MAX_VALUE}){
			final long expected=limit==0 ? 1 : limit==1 ? 2 : 3;
			for(int route=0; route<3; route++){
				check(patienceWork(limit, false, false, route)==expected,
					"Actual Scrabble/Quantum/reused work must distinguish legacy/full, zero and one: patience="+limit+", route="+route);
			}
		}
		check(patienceWork(-1, true, false, 0)==1, "Legacy retains the ranked-fallback first-pass exit");
		check(patienceWork(0, true, true, 0)==3, "Exhaustive diagnostics bypass even zero patience");
	}

	/** All candidates fit one physical window; its first ranked model is exact. */
	private static long patienceWork(final int limit, final boolean ranked, final boolean oracle, final int route){
		final byte[][] library={bytes("ACGT"), bytes("AA"), bytes("CC")};
		final map.LongHashSet seeds=new map.LongHashSet();
		seeds.add(0L);
		final NcrnaScavenger caller=new NcrnaScavenger(library, null, null, seeds, 1, 1, 4,
			1, 3, false, 0, 0, 0, 1);
		caller.scavengePass2=false;
		caller.trimAlignmentExtent=false;
		caller.quantumThresh=route==0 ? 120 : 0;
		caller.reuseConsensusAlignment=route==2;
		caller.rankedModelFallback=ranked;
		caller.setShortlistPolicy(3, false, 0, 0, 0, 1, limit);
		caller.setMappingOracleExhaustive(oracle);
		final NcrnaCpuProfile profile=NcrnaCpuProfile.create(1);
		if(profile!=null){profile.begin(0, NcrnaCpuProfile.SCAVENGE);}
		try{caller.scavenge("patience_fixture", bytes("ACGT"), 0, new ArrayList<int[]>(), new int[]{0});}
		finally{if(profile!=null){profile.end();}}
		if(profile!=null){
			check(profile.alignmentWork(0).logical==caller.alignmentCount(),
				"Each fixture verification must be an actual logical alignment request");
		}
		return caller.alignmentCount();
	}

	/** A fixed cutoff changes the actual index without switching its mode. */
	private static void fixedShortlist(){
		final byte[][] library={bytes("ACGT"), bytes("AA"), bytes("CC")};
		final byte[] query=bytes("ACGT");
		final TrnaKmerIndex fixed=new TrnaKmerIndex(library, 1, false, 0, 0, 0, 1);
		fixed.setFixedMinHits(3);
		check(Arrays.equals(fixed.shortlist(query, 3), new int[]{0}), "Fixed cutoff 3 removes both count-two models");
		fixed.setFixedMinHits(2);
		check(Arrays.equals(fixed.shortlist(query, 3), new int[]{0, 1, 2}), "Fixed equality restores tied models in their original order");
		fixed.setFixedMinHits(5);
		check(Arrays.equals(fixed.shortlist(query, 3), new int[]{0}), "Fixed multi-model keep-best behavior is retained");
		rejected(() -> fixed.setAdaptiveFloor(2), "A fixed index cannot change adaptive floors");
		rejected(() -> fixed.setFixedMinHits(-1), "Negative fixed cutoff must fail");
		final TrnaKmerIndex single=new TrnaKmerIndex(new byte[][]{bytes("ACGT")}, 1, false, 0, 0, 0, 1);
		check(single.shortlist(bytes("NN"), 1).length==0, "No shared keys fail a positive single-model cutoff");
		single.setFixedMinHits(0);
		check(single.shortlist(bytes("NN"), 1).length==1, "Zero fixed cutoff is the actual direct-alignment control");
		single.setFixedMinHits(4);
		check(single.shortlist(query, 1).length==1, "Single-model fixed equality passes");
		single.setFixedMinHits(5);
		check(single.shortlist(query, 1).length==0, "Single-model fixed cutoff has no keep-best exception");
		final map.LongHashSet seeds=new map.LongHashSet();
		seeds.add(0L);// The provided A seed at position zero belongs to this sequence.
		final NcrnaScavenger caller=new NcrnaScavenger(library, null, null, seeds, 1, 1, 4,
			1, 3, false, 0, 0, 0, 1);
		caller.scavengePass2=false;
		caller.trimAlignmentExtent=false;
		caller.scavenge("fixture", query, 0, new ArrayList<int[]>(), new int[]{0});
		check(caller.alignmentCount()==3, "Fixed cutoff one must execute all three model verifications");
		caller.setFixedShortlistPolicy(2, 3);
		caller.scavenge("fixture", query, 0, new ArrayList<int[]>(), new int[]{0});
		check(caller.alignmentCount()==4, "Changing the worker cutoff must actually remove both count-two model alignments");
		check(!caller.adaptiveMinHits && caller.indexMinHitsDefault==3 && caller.indexTopN==2,
			"The worker's strict-cutoff metadata must follow its fixed index without changing mode");
		rejectedState(() -> caller.setShortlistPolicy(2, 3), "A fixed caller cannot consume an adaptive-floor setter");
		for(final String text:new String[]{"2.00000001", "2.5", "2147483648", "-1"}){
			final byte[] field=bytes(text);
			rejected(() -> NcrnaShortlistTarget.fixedHits(field, 0, field.length), "Fixed integer text must not round or overflow: "+text);
		}
		final byte[] large=bytes("16777217"), maximum=bytes("2147483647");
		check(NcrnaShortlistTarget.fixedHits(large, 0, large.length)==16777217,
			"An integer beyond float's exact range must retain every digit");
		check(NcrnaShortlistTarget.fixedHits(maximum, 0, maximum.length)==Integer.MAX_VALUE,
			"The largest valid fixed cutoff must parse without float conversion");
		caller.setFixedShortlistPolicy(2, 16777217);
		check(caller.indexMinHitsDefault==16777217, "The caller must retain the exact parsed integer cutoff");
	}

	/** Explicit member and view choices prevent a valid-looking wrong-family sweep. */
	private static void targetSelection(){
		fixedShortlist();
		final NcrnaFamily tpp=new NcrnaFamily("tpp", new byte[][]{bytes("ACGT")}, null,
			new String[]{"tpp_model"}, null, 17, 50, 232, 7, 361, true, 5, .1875f, 0, 1);
		final NcrnaFamily sam=new NcrnaFamily("sam1", new byte[][]{bytes("ACGT")}, null,
			new String[]{"SAM1_model"}, null, 17, 50, 378, 7, 238, false, 0, 0f, 0f, 1);
		sam.modelPatience=16;
		final Rule[] tppRules={new Rule("TPP", "tpp", "TPP", true), new Rule("SAM_I", "SAM1", "SAM-I", false)};
		final Rule[] samRules={new Rule("TPP", "tpp", "TPP", false), new Rule("SAM_I", "SAM1", "SAM-I", true)};
		final NcrnaShortlistTarget tppTarget=new NcrnaShortlistTarget("tpp", "target");
		final NcrnaShortlistTarget samTarget=new NcrnaShortlistTarget("sam1", "target");
		check(tppTarget.bind(Arrays.asList(tpp, sam), tppRules, "TPP")==0, "TPP keeps its original caller slot");
		check(samTarget.bind(Arrays.asList(tpp, sam), samRules, "SAM_I")==1, "SAM-I must select its own caller, not TPP");
		rejected(() -> samTarget.bind(Arrays.asList(tpp, sam), tppRules, "TPP"), "SAM-I with a TPP grading map must fail before calling");
		rejected(() -> samTarget.bind(Arrays.asList(tpp), samRules, "SAM_I"), "An absent target cannot silently select another caller");
		final Rule[] mislabeled={new Rule("TPP", "tpp", "TPP", false), new Rule("TPP", "SAM1", "SAM-I", true)};
		rejected(() -> samTarget.bind(Arrays.asList(tpp, sam), mislabeled, "TPP"),
			"Correct SAM-I model/display labels cannot authorize grading against TPP reference truth");
		tppTarget.validate(tpp, true, 361, 5, -1);
		samTarget.validate(sam, false, 238, 0, 1);
		rejected(() -> samTarget.validate(sam, true, 238, 0, 1), "Adaptive floor header is wrong for packaged fixed SAM-I");
		rejected(() -> samTarget.validate(sam, true, 361, 5, -1), "SAM-I must not inherit TPP baseline constants");
		final NcrnaFamily sixs=new NcrnaFamily("sixs_rf00013", new byte[][]{bytes("ACGT")}, null,
			new String[]{"6S_RF00013_model"}, null, 17, 60, 250, 7, 56, true, 1, 0, 0, 1);
		sixs.modelPatience=0;
		final NcrnaShortlistTarget sixsTarget=new NcrnaShortlistTarget("sixs_rf00013", "union");
		sixsTarget.validate(sixs, true, 56, 1, 1);
		rejected(() -> sixsTarget.validate(sixs, false, 56, 1, 1), "Adaptive 6S must not silently switch to a fixed policy table");
		for(final String aggregate:new String[]{"cobalamin", "6S", "SRP"}){
			rejected(() -> new NcrnaShortlistTarget(aggregate, "union"), "Aggregate names cannot choose a first member");
		}
		rejected(() -> new NcrnaShortlistTarget("unknown", "target"), "Unknown families cannot fall back to TPP");
		rejected(() -> new NcrnaShortlistTarget("sam1", "union"), "Standalone family is not a multi-member union");
		rejected(() -> new NcrnaShortlistTarget("cobalamin_rf00174", "target"), "Member-only queries must not be graded against union truth implicitly");
		final String[] members={"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"};
		final NcrnaFamily[] families=new NcrnaFamily[3];
		final Rule[] rules=new Rule[3];
		for(int i=0; i<members.length; i++){
			families[i]=new NcrnaFamily(members[i], new byte[][]{bytes("ACGT")}, null,
				new String[]{members[i]+"_model"}, null, 17, 109, 450);
			rules[i]=new Rule("cobalamin", members[i], members[i], true);
		}
		final NcrnaShortlistTarget member=new NcrnaShortlistTarget(members[1], "union");
		check(member.bind(Arrays.asList(families), rules, "cobalamin")==1,
			"A union-graded member sweep selects the named index while retaining every sibling's calls");
		final NcrnaShortlistTarget joint=new NcrnaShortlistTarget("cobalamin", "union", true);
		check(Arrays.equals(joint.bindTargets(Arrays.asList(families), rules, "cobalamin"), new int[]{0, 1, 2}),
			"Explicit joint mode must vary all three loaded cobalamin members");
		check(Arrays.equals(member.bindTargets(Arrays.asList(families), rules, "cobalamin"), new int[]{1}),
			"Legacy member mode must not silently expand to its whole group");
		rejected(() -> member.bind(Arrays.asList(families[0], families[1]), rules, "cobalamin"),
			"Missing sibling callers cannot masquerade as complete union grading");
	}

	private static void rejected(final Runnable action, final String reason){
		try{action.run();}
		catch(IllegalArgumentException expected){return;}
		throw new AssertionError(reason);
	}

	private static void rejectedState(final Runnable action, final String reason){
		try{action.run();}
		catch(IllegalStateException expected){return;}
		throw new AssertionError(reason);
	}

	private static void projection(final String path){
		final Rule[] rules={new Rule("TPP", "TPP", "TPP", true), new Rule("FMN", "FMN", "FMN", false)};
		final NcrnaReferenceCalls reader=new NcrnaReferenceCalls(rules);
		final Result live=reader.result("genome-A");
		final byte[] bases=new byte[100];
		Arrays.fill(bases, (byte)'A');
		final Orf[] orfs={new Orf("seq%20literal description", 9, 28, 1, 0, bases, false, ProkObject.RNA),
			new Orf("seq%20literal description", 9, 28, 1, 0, bases, false, ProkObject.RNA),
			new Orf("seq%20literal description", 40, 60, 0, 0, bases, false, ProkObject.RNA),
			new Orf("seq%20literal description", 0, 10, 0, 0, bases, false, ProkObject.CDS)};
		orfs[0].trnaModel="TPP_c1 n=7";
		orfs[0].ncrnaFamily="tpp";
		orfs[1].trnaModel="TPP_c1 n=7";
		orfs[1].ncrnaFamily="tpp";
		orfs[2].trnaModel="FMN_c1";
		orfs[2].ncrnaFamily="fmn";
		final ByteBuilder b=new ByteBuilder();
		try(Output out=new Output(path)){
			for(final Orf orf:orfs){
				reader.add(live, orf);
				b.clear();
				out.writer.print(orf.appendGff(b).nl());
			}
		}
		final Result saved=reader.read(path, "genome-A");
		check(live.features==4 && live.cds==1 && live.included==2 && live.excluded==1,
			"Two retained duplicates, one excluded RNA and one CDS conserve four features");
		check(Arrays.equals(live.ruleCounts, saved.ruleCounts) && saved.features==live.features,
			"Live and saved population accounting agrees");
		for(int i=0; i<live.calls.size(); i++){
			final Call a=live.calls.get(i), c=saved.calls.get(i);
			check(a.locus.seqid.equals(c.locus.seqid) && a.locus.family.equals(c.locus.family)
				&& a.locus.start==c.locus.start && a.locus.stop==c.locus.stop && a.locus.strand==c.locus.strand
				&& a.model.equals(c.model) && a.ruleIndex==c.ruleIndex && a.included==c.included,
				"Direct final-path projection must equal the native saved-GFF reader");
		}
		check(live.calls.get(0).locus.seqid.equals("seq%20literal") && live.calls.get(0).locus.start==10
			&& live.calls.get(0).locus.stop==29 && live.calls.get(0).locus.strand=='-',
			"Literal header token, one-based closed bounds and minus strand survive");
		final Orf unicode=new Orf("seqΩ", 9, 28, 0, 0, bases, false, ProkObject.RNA);
		unicode.trnaModel="TPP_c1";
		try{
			reader.add(reader.result("genome-A"), unicode);
			throw new AssertionError("Lossy live identifier accepted");
		}
		catch(IllegalArgumentException expected){check(expected.getMessage().contains("ASCII"), "Expected explicit byte-writer limitation");}
		try(Output out=new Output(path+".nonascii")){
			b.clear();
			out.writer.print(unicode.appendGff(b).nl());
		}
		try{
			reader.read(path+".nonascii", "genome-A");
			throw new AssertionError("Corrupt saved identifier accepted");
		}
		catch(IllegalArgumentException expected){check(expected.getMessage().contains("UTF-8"), "Expected corrupt UTF-8 diagnosis");}

	}

	/** Restore full candidate membership and raw state, not the previous winning path. */
	private static void candidateState(){
		final byte[] bases=new byte[100];
		Arrays.fill(bases, (byte)'A');
		final Orf cds=new Orf("fixture", 3, 47, 0, 0, bases, false, ProkObject.CDS);
		final Orf competing=new Orf("fixture", 30, 60, 0, 0, bases, false, ProkObject.RNA);
		final NcrnaReportedOrf glycine=new NcrnaReportedOrf("fixture", 20, 40, 1, bases);
		glycine.start=18;
		glycine.stop=43;
		glycine.captureRefinement(20, 40, 3, 100);
		glycine.flip();
		final int rawStart=glycine.start, rawStop=glycine.stop;
		final NcrnaCandidateReuse reuse=new NcrnaCandidateReuse(0, 2);
		final GeneModel pgm=new GeneModel(false);
		@SuppressWarnings("unchecked")
		final ArrayList<Orf>[] frames=new ArrayList[6];
		for(int i=0; i<frames.length; i++){frames[i]=new ArrayList<Orf>();}
		frames[0].add(cds);
		reuse.captureCds(frames, pgm);
		final ArrayList<Orf> rnas=new ArrayList<Orf>();
		rnas.add(glycine);
		reuse.captureRna(3, rnas);
		reuse.captureRna(2, null);
		frames[0].add(competing);// Ordinary RNA merge must not pollute the retained CDS population.
		cds.prevPlus=competing;
		cds.prevMinus=glycine;
		cds.pathScorePlus=91;
		cds.pathScoreMinus=87;
		cds.pathLengthPlus=7;
		cds.pathLengthMinus=5;
		glycine.prevPlus=cds;
		glycine.applyReportedBounds();
		check(glycine.start!=rawStart, "Fixture must exercise the deferred minus-strand reporting mutation");
		final ArrayList<Orf>[] restored=reuse.restoreCds(pgm);
		check(restored[0].size()==1 && restored[0].get(0)==cds && cds.prevPlus==null && cds.prevMinus==null
			&& cds.pathScorePlus==0 && cds.pathScoreMinus==0 && cds.pathLengthPlus==1 && cds.pathLengthMinus==1,
			"Each policy starts from the complete original CDS candidates with no previous winner chain");
		check(reuse.restoreRna(3).get(0)==glycine && glycine.start==rawStart && glycine.stop==rawStop && glycine.prevPlus==null,
			"Minus-strand RNA discovery bounds and path state survive reporting without another strand flip");
		check(reuse.restoreRna(2).isEmpty() && reuse.restoreRna(0)==null,
			"An empty unchanged family is reusable; target candidates remain absent from the cache");
		NcrnaReportedOrf.applyTo(rnas);
		check(glycine.start!=rawStart, "Final baseline output must recover its independently stored reporting bounds");
	}
	/** A group excludes all varied members while retaining unrelated RNA and tRNA. */
	private static void groupMembership(){
		final int[] members={0, 2};
		final NcrnaCandidateReuse reuse=new NcrnaCandidateReuse(members, 3);
		members[0]=1;
		check(reuse.isTarget(0) && reuse.isTarget(2) && !reuse.isTarget(1) && !reuse.isTarget(3),
			"Copied membership must exclude both varied members, retain the unrelated caller, and leave independent tRNA reusable");
		reuse.captureRna(2, null);
		reuse.captureRna(6, null);
		check(reuse.restoreRna(0)==null && reuse.restoreRna(4)==null
			&& reuse.restoreRna(2).isEmpty() && reuse.restoreRna(6).isEmpty(),
			"Selected members must be regenerated while empty unchanged-family and tRNA results remain reusable");
		rejected(() -> new NcrnaCandidateReuse(new int[]{0, 0}, 3), "Duplicate group members must fail");
		rejected(() -> new NcrnaCandidateReuse(new int[]{3}, 3), "Independent tRNA cannot be a generic-family target");
		rejected(() -> new NcrnaCandidateReuse(new int[0], 3), "An empty group cannot define a policy experiment");
	}

	/** Literal table distinguishes a fixed root from a later adaptive joint policy. */
	private static void jointTable(final String path){
		try(Output out=new Output(path)){
			out.writer.println(NcrnaJointPolicy.HEADER);
			out.writer.println("root\tsixs_rf00013\tfixed\t60\t0\t1\t0\t0\t-1");
			out.writer.println("root\tsixs_rf01685\tfixed\t60\t0\t1\t0\t0\t-1");
			out.writer.println("neighbor\tsixs_rf00013\tadaptive\t30\t4\t1\t0.25\t0.125\t0");
			out.writer.println("neighbor\tsixs_rf01685\tadaptive\t40\t3\t1\t0.5\t0.25\t2147483647");
		}
		final NcrnaJointPolicy[] policies=NcrnaJointPolicy.load(path);
		check(policies.length==2 && policies[0].size()==2 && policies[1].size()==2,
			"Contiguous member rows must become two complete policies, not four independent policies");
		final NcrnaJointPolicy.Member root=policies[0].member(0), neighbor=policies[1].member(1);
		check(!root.adaptive && root.topN==60 && root.floor==0 && root.minHits==1 && root.topFrac==0 && root.qFrac==0 && root.patience==-1,
			"The fixed input tuple remains unchanged when the later policy switches modes");
		check(neighbor.name.equals("sixs_rf01685") && neighbor.adaptive && neighbor.topN==40
			&& neighbor.floor==3 && neighbor.minHits==1 && neighbor.topFrac==.5f && neighbor.qFrac==.25f
			&& neighbor.patience==Integer.MAX_VALUE && policies[1].member(0).patience==0,
			"Each member retains all five active controls, including exact integer patience");
		rejected(() -> new NcrnaJointPolicy("duplicate", new NcrnaJointPolicy.Member[]{root, root}),
			"A duplicate member cannot count as complete group coverage");
		try(Output out=new Output(path+".missing")){
			out.writer.println(NcrnaJointPolicy.HEADER);
			out.writer.println("root\tsixs_rf00013\tfixed\t60\t0\t1\t0\t0\t-1");
			out.writer.println("root\tsixs_rf01685\tfixed\t60\t0\t1\t0\t0\t-1");
			out.writer.println("neighbor\tsixs_rf00013\tadaptive\t30\t4\t1\t0.25\t0.125\t0");
		}
		rejected(() -> NcrnaJointPolicy.load(path+".missing"), "An incomplete later member roster must fail before caller execution");
		try(Output out=new Output(path+".overflow")){
			out.writer.println(NcrnaJointPolicy.HEADER);
			out.writer.println("root\tsixs_rf00013\tfixed\t60\t0\t1\t0\t0\t2147483648");
		}
		rejected(() -> NcrnaJointPolicy.load(path+".overflow"), "Patience input must not wrap into the valid legacy -1 sentinel");
	}

	private static byte[] bytes(final String value){return value.getBytes(StandardCharsets.US_ASCII);}
	private static void check(final boolean condition, final String reason){if(!condition){throw new AssertionError(reason);}}
}
