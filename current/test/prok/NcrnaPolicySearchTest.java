package test.prok;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaJointPolicy;
import prok.NcrnaJointPolicy.Member;
import prok.NcrnaPolicySearch;
import prok.NcrnaPolicySearch.Axis;
import prok.NcrnaPolicySearch.Measurement;
import prok.OptimizeNcrnaShortlist;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NcrnaSearchCoverage.Counts;
import shared.Shared;
import structures.ByteBuilder;

/** Finite measured landscape, independent of the caller and search implementation.
 * The best first move changes after the first step; one-axis conclusions cannot
 * substitute for recomputing joint neighbors. Additional literal cases cover
 * free gains, ties, plateaus, unreachable targets and paired-reference limits.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaPolicySearchTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaPolicySearchTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown test parameter: "+arg);
			}
			interactions();
			selection();
			coverageSelection();
			gridValidation();
			check(parser.out1!=null, "Supply out=fresh_fixture_directory for driver restart tests");
			lifecycle(new File(parser.out1).getAbsoluteFile());
			pp.outstream.println("NCRNA_POLICY_SEARCH_FINITE_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void interactions(){
		// Rows are cap, columns are patience. Other active axes penalize this face.
		final long[][] tp={{90, 94, 95}, {96, 100, 99}, {97, 100, 100}};
		final long[][] cost={{100, 110, 125}, {130, 115, 90}, {90, 105, 110}};
		final int[][] expected={{0, 0}, {0, 1}, {1, 1}, {2, 1}};
		final NcrnaPolicySearch search=search();
		final boolean[][] visited=new boolean[3][3];
		int[] center=search.initial();
		int steps=0;
		for(; steps<10; steps++){
			check(steps<expected.length && center[0]==expected[steps][0] && center[4]==expected[steps][1],
				"Expected joint path "+Arrays.deepToString(expected)+", reached "+Arrays.toString(center)+" at step "+steps);
			check(!visited[center[0]][center[4]], "A strictly improving finite walk must not cycle");
			visited[center[0]][center[4]]=true;
			final ArrayList<int[]> points=new ArrayList<int[]>();
			final ArrayList<Measurement> values=new ArrayList<Measurement>();
			for(int a=0; a<search.dimensions(); a++){
				for(int d=-1; d<=1; d+=2){
					final int[] next=search.neighbor(center, a, d);
					if(next!=null && (next[1]+next[2]+next[3]>0 || !visited[next[0]][next[4]])){
						points.add(next);
						final int penalty=5*(next[1]+next[2]+next[3]);
						values.add(m(tp[next[0]][next[4]]-penalty, cost[next[0]][next[4]]+penalty));
					}
				}
			}
			if(steps==0){check(points.size()==5, "Every one of five nonconstant axes must be evaluated from the first center");}
			final int selected=NcrnaPolicySearch.choose(m(tp[center[0]][center[4]], cost[center[0]][center[4]]),
				values.toArray(new Measurement[values.size()]), 100, -1);
			if(selected<0){break;}
			center=points.get(selected);
		}
		check(steps==3 && tp[center[0]][center[4]]==100 && cost[center[0]][center[4]]==105,
			"Joint neighbor reevaluation must restore all 100 loci, then lower cost from 115 to 105");
		final Member selected=search.policy("chosen", center).member(0);
		check(selected.topN==3 && selected.floor==1 && selected.topFrac==0 && selected.qFrac==0
			&& selected.patience==1 && selected.adaptive,
			"The measured coordinates must materialize the exact native caller tuple");
		final int[] borrowed=search.initial();
		borrowed[0]=2;
		check(search.initial()[0]==0, "External coordinate mutation must not alter the saved starting center");
	}

	private static void selection(){
		check(NcrnaPolicySearch.choose(m(90, 100), new Measurement[]{m(94, 100), m(99, 110)}, 100, -1)==0,
			"A free recall gain takes precedence over a paid one");
		check(NcrnaPolicySearch.choose(m(90, 100), new Measurement[]{m(90, 80), m(95, 120)}, 100, -1)==0,
			"Loss-free savings are useful even before recall is restored");
		check(NcrnaPolicySearch.choose(m(90, 100), new Measurement[]{m(94, 110), m(94, 110)}, 100, -1)==0,
			"An exact tie uses stable neighbor order");
		check(NcrnaPolicySearch.choose(m(100, 100), new Measurement[]{m(99, 60), m(100, 95)}, 100, -1)==1,
			"Tightening cannot trade away the required recall for cheaper work");
		check(NcrnaPolicySearch.choose(m(100, 100), new Measurement[]{m(100, 100), m(99, 99)}, 100, -1)==-1,
			"An equal plateau or infeasible cheap neighbor stops locally rather than cycling");
		check(NcrnaPolicySearch.choose(m(90, 100), new Measurement[]{m(89, 80), m(90, 110)}, 100, -1)==-1,
			"An unreachable recall target must produce a finite local stop");
		final Measurement replacement=new Measurement(100, 0, 90, 3), preserved=new Measurement(99, 0, 110, 0);
		check(NcrnaPolicySearch.choose(m(98, 100), new Measurement[]{replacement, preserved}, 100, 0)==1,
			"Restored TP totals do not excuse explicit paired-reference losses under maxlost=0");
		check(NcrnaPolicySearch.choose(m(98, 100), new Measurement[]{replacement, preserved}, 100, -1)==0,
			"Without a paired-locus restriction, the same evidence must remain an explicit alternative");
	}

	private static void gridValidation(){
		rejected(() -> new Axis(0, 0, new double[]{1, 1.5}), "Caps cannot be fractional");
		rejected(() -> new Axis(0, 2, new double[]{.5, .50000000001}), "Distinct text that rounds to the same float is not two policies");
		rejected(() -> new Axis(0, 3, new double[]{Double.NaN}), "Nonfinite query fractions cannot enter a grid");
		rejected(() -> new Axis(0, 4, new double[]{-1, 0}), "Legacy patience is baseline-only, outside the optimization grid");
		rejected(() -> new Axis(0, 4, new double[]{1.5}), "Patience must retain an exact integer alignment count");
		rejected(() -> new Axis(0, 4, new double[]{2147483648d}), "Patience cannot overflow the caller's integer type");
		final Axis first=new Axis(0, 0, new double[]{1, 2, 3});
		rejected(() -> new NcrnaPolicySearch(start(), new Axis[]{first, first, first, first, first}),
			"Five repeated caps are not five independent controls");
		check(search().neighbor(new int[]{0, 0, 0, 0, 0}, 0, -1)==null, "Lower grid boundaries have no invented neighbor");
	}

	/** Literal coverage outcomes are independent of the selector implementation. */
	private static void coverageSelection(){
		final Measurement center=covered(100, 100, 0, 0, 0);
		final Measurement replacement=covered(100, 60, 1, 1, 0), retained=covered(100, 90, 0, 0, 0);
		check(NcrnaPolicySearch.choose(center, new Measurement[]{replacement, retained}, 100, -1, -1, -1)==0,
			"Disabled coverage must retain the original cheapest admissible choice");
		check(NcrnaPolicySearch.choose(center, new Measurement[]{replacement, retained}, 100, -1, 0, -1)==1,
			"A first detection in another input cannot refund the last-detection veto");
		check(NcrnaPolicySearch.choose(center, new Measurement[]{replacement}, 100, -1, 1, -1)==0,
			"A one-input allowance permits one explicit root-relative exception");
		check(NcrnaPolicySearch.choose(replacement, new Measurement[]{covered(100, 40, 2, 2, 0)}, 100, -1, 1, -1)==-1,
			"The allowance does not renew when a later center loses a second shipped input");
		check(NcrnaPolicySearch.choose(center, new Measurement[]{covered(101, 50, 0, 0, 1), retained}, 100, -1, 0, 0)==1,
			"A dominating union gain cannot bypass a member last-detection veto");
		check(NcrnaPolicySearch.choose(covered(98, 100, 1, 0, 0),
			new Measurement[]{covered(99, 110, 0, 0, 0)}, 100, -1, 0, -1)==0,
			"An infeasible starting center may take an admissible improving partial-recovery step below target");
		check(NcrnaPolicySearch.choose(center, new Measurement[]{covered(99, 20, 0, 0, 0)}, 100, -1, 0, -1)==-1,
			"Preserving coverage does not waive an already-attained total-TP floor");
		final Measurement locusLoss=new Measurement(100, 0, 60, 1, new Counts(0, 0, new long[0], new long[0]));
		check(NcrnaPolicySearch.choose(center, new Measurement[]{locusLoss}, 100, 0, 0, -1)==-1,
			"A same-input replacement still fails maxlost=0 even when last-loss coverage passes");
		rejected(() -> NcrnaPolicySearch.choose(m(100, 100), new Measurement[]{m(100, 60)}, 100, -1, 0, -1),
			"Missing coverage is unavailable, not zero");
		rejected(() -> NcrnaPolicySearch.choose(center, new Measurement[]{retained}, 100, -1, -2, -1),
			"A budget below the disabled sentinel cannot be accepted");
		final long[] borrowed={1};
		final Counts counts=new Counts(0, 0, borrowed, new long[]{0});
		borrowed[0]=0;
		check(counts.memberLastLoss(0)==1, "Published coverage owns its member vector rather than borrowing mutable client arrays");
	}

	private static Measurement covered(final long tp, final long cost, final long last, final long first, final long memberLast){
		return new Measurement(tp, 0, cost, 0, new Counts(last, first, new long[]{memberLast}, new long[]{0}));
	}

	private static NcrnaJointPolicy start(){
		return new NcrnaJointPolicy("start", new Member[]{new Member("tpp", true, 1, 1, 1, 0, 0, 0)});
	}

	/** Exercise the actual driver state machine with a labeled synthetic evaluator.
	 * Replaying a completed search must invoke the evaluator zero additional times.
	 */
	private static void lifecycle(final File directory){
		try{Files.createDirectory(directory.toPath());}
		catch(IOException e){throw new IllegalStateException("Driver fixtures need a fresh directory", e);}
		write(new File(directory, "caller.config"), "in=synthetic_unused.fa\n");
		write(new File(directory, "references.tsv"), "synthetic reference identity only\n");
		write(new File(directory, "familymap.tsv"), "synthetic projection identity only\n");
		// Literal fixture identities are not claimed to fingerprint biological data.
		write(new File(directory, "identity.tsv"), "source_commit\t0000000000000000000000000000000000000000\n"
			+"build_seal\t00000000000000000000\nresources_sha80\t00000000000000000000\n"
			+"inputs_sha80\t00000000000000000000\nreferences_sha80\t00000000000000000000\n");
		write(new File(directory, "initial.tsv"), NcrnaJointPolicy.HEADER+"\n"
			+"root\ttpp\tadaptive\t3\t1\t1\t0\t0.0001\t-1\nstart\ttpp\tadaptive\t1\t1\t1\t0\t0.0001\t0\n");
		write(new File(directory, "grid.tsv"), "member\tknob\tvalues\ntpp\ttopn\t1,2,3\n"
			+"tpp\tfloor\t1,2\ntpp\ttopfrac\t0,0.25\ntpp\tqfrac\t0,0.0001\n"
			+"tpp\tpatience\t0,1,2,4,8,16,32,2147483647\n");
		final String[] arguments={"caller="+new File(directory, "caller.config"),
			"reference="+new File(directory, "references.tsv"), "familymap="+new File(directory, "familymap.tsv"),
			"initial="+new File(directory, "initial.tsv"), "grid="+new File(directory, "grid.tsv"),
			"identity="+new File(directory, "identity.tsv"), "target=tpp", "grading=target", "targetrecall=1",
			"maxrounds=10", "out="+new File(directory, "search")};
		final SyntheticRunner runner=new SyntheticRunner();
		new OptimizeNcrnaShortlist(arguments, runner).process();
		check(runner.calls==4, "The literal four-center path must use four finite evaluator rounds");
		final String result=read(new File(directory, "search/result.tsv"));
		check(result.contains("synthetic_test\tLOCAL_POINT_MET_TARGET\t100\t0\t105\t0\t100\n"),
			"Driver must publish the independently expected local point and label it synthetic");
		final NcrnaJointPolicy[] chosen=NcrnaJointPolicy.load(new File(directory, "search/selected.tsv").getPath());
		final NcrnaJointPolicy[] firstRound=NcrnaJointPolicy.load(new File(directory, "search/round0/policies.tsv").getPath());
		final NcrnaJointPolicy[] grid=NcrnaJointPolicy.load(new File(directory, "search/grid_policies.tsv").getPath());
		boolean zero=false, small=false;
		for(final NcrnaJointPolicy policy:grid){
			zero|=policy.member(0).qFrac==0;
			small|=policy.member(0).qFrac==.0001f;
		}
		check(grid.length==17 && zero && small, "Actual grid file must retain all seventeen choices, including the complete eight-value patience grid");
		final int[] expectedCaps={1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
		final float[] expectedFloors={1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
		final int[] expectedPatience={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 4, 8, 16, 32, Integer.MAX_VALUE};
		for(int i=0; i<grid.length; i++){
			final Member member=grid[i].member(0);
			check(member.name.equals("tpp") && member.adaptive && member.topN==expectedCaps[i]
				&& member.floor==expectedFloors[i] && member.minHits==1 && member.topFrac==(i==6 ? .25f : 0)
				&& member.qFrac==(i==7 ? 0 : .0001f) && member.patience==expectedPatience[i],
				"Every reloaded grid choice must match the independent literal control tuple at row "+i);
		}
		check(chosen[1].member(0).qFrac==.0001f && firstRound[1].member(0).qFrac==.0001f,
			"Actual generated round and selected policy files must reload a small fraction without exponent parsing");
		check(chosen[1].member(0).topN==3 && chosen[1].member(0).floor==1 && chosen[1].member(0).patience==1,
			"Saved policy must reproduce the measured selected coordinates");
		new OptimizeNcrnaShortlist(with(arguments, "resume", "t"), runner).process();
		check(runner.calls==4 && read(new File(directory, "search/result.tsv")).equals(result),
			"Completed-round resume must preserve results without repeating evaluator work");
		new OptimizeNcrnaShortlist(with(with(arguments, "resume", "t"), "cost", "alignments"), runner).process();
		check(runner.calls==4 && read(new File(directory, "search/result.tsv")).equals(result),
			"Explicit default cost must resume the unchanged legacy contract and result");
		rejected(() -> new OptimizeNcrnaShortlist(with(with(arguments, "resume", "t"), "cost", "modelbases"), runner).process(),
			"Changing work units must reject old cached decisions before evaluator work");
		final String[] weighted=with(with(arguments, "out", new File(directory, "modelbases").getPath()), "cost", "modelbases");
		final SyntheticRunner weightedRunner=new SyntheticRunner(true);
		new OptimizeNcrnaShortlist(weighted, weightedRunner).process();
		check(read(new File(directory, "modelbases/result.tsv")).equals(
			"backend\tstatus\tTP\tFP\taligned_model_bases\tlost_root\ttarget_TP\n"
			+"synthetic_test\tLOCAL_POINT_MET_TARGET\t100\t0\t200\t0\t100\n"),
			"The independent model-length landscape selects cap two, unlike alignment-count cap three");
		final int weightedCalls=weightedRunner.calls;
		new OptimizeNcrnaShortlist(with(weighted, "resume", "t"), weightedRunner).process();
		check(weightedRunner.calls==weightedCalls, "Model-base resume must not rerun completed rounds");
		rejected(() -> new OptimizeNcrnaShortlist(with(with(weighted, "resume", "t"), "cost", "alignments"), weightedRunner).process(),
			"Model-base decisions cannot be resumed with alignment-count units");
		final String[] changed=with(with(arguments, "resume", "t"), "targetrecall", ".95");
		rejected(() -> new OptimizeNcrnaShortlist(changed, runner).process(), "A changed objective must invalidate resume");
		check(runner.calls==4, "Rejected resume must not execute a round");
		final File gridFile=new File(directory, "grid.tsv");
		final String originalGrid=read(gridFile);
		final File changedGrid=new File(directory, "changed_grid.tsv");
		write(changedGrid, originalGrid.replace("0,1,2,4,8,16,32,2147483647", "0,1,2,4,8,16,32,64,2147483647"));
		rejected(() -> new OptimizeNcrnaShortlist(with(with(arguments, "resume", "t"), "grid", changedGrid.getPath()), runner).process(),
			"Changing the patience grid must invalidate resume before evaluation");
		check(runner.calls==4, "Changed patience identity cannot repeat completed scientific work");
		final String[] limited=with(with(arguments, "out", new File(directory, "limited").getPath()), "maxrounds", "1");
		new OptimizeNcrnaShortlist(limited, runner).process();
		check(read(new File(directory, "limited/result.tsv")).contains("synthetic_test\tROUND_LIMIT\t94\t0\t110\t6\t100\n"),
			"At the round limit retain the best already measured neighbor, not the old center");
		final String[] interrupted=with(arguments, "out", new File(directory, "interrupted").getPath());
		final int[] calls={0};
		final OptimizeNcrnaShortlist.RoundRunner failure=(prefix, policies) -> {
			calls[0]++;
			throw new IllegalStateException("Deliberate interrupted evaluator fixture");
		};
		try{
			new OptimizeNcrnaShortlist(interrupted, failure).process();
			throw new AssertionError("Injected interruption did not reach the driver");
		}catch(IllegalStateException expected){check(calls[0]==1, "Failure must occur inside the owned round");}
		rejected(() -> new OptimizeNcrnaShortlist(with(interrupted, "resume", "t"), failure).process(),
			"An incomplete round must require inspection rather than an automatic repeat");
		check(calls[0]==1, "Incomplete-round rejection must preserve work without rerunning it");
	}

	private static final class SyntheticRunner implements OptimizeNcrnaShortlist.RoundRunner{
		SyntheticRunner(){this(false);}
		SyntheticRunner(final boolean modelBases_){modelBases=modelBases_;}
		@Override
		public void run(final File prefix, final ArrayList<NcrnaJointPolicy> policies){
			calls++;
			if(calls==1){check(policies.size()==7, "The first measured round needs the root, center and one neighbor on each of five active axes");}
			final long[][] tp={{90, 94, 95}, {96, 100, 99}, {97, 100, 100}};
			final long[][] cost={{100, 110, 125}, {130, 115, 90}, {90, 105, 110}};
			final ByteBuilder rows=new ByteBuilder();
			rows.append("policy\tTP\tFP\tFN\tscored_references\tlost_baseline_TP\tgained_TP\tkernel_logical_alignments");
			if(modelBases){rows.append("\taligned_model_bases");}
			rows.nl();
			for(final NcrnaJointPolicy policy:policies){
				final Member member=policy.member(0);
				final int a=member.topN-1, b=member.patience<0 ? 1 : Math.min(2, member.patience);
				final int penalty=5*((member.floor==1 ? 0 : 1)+(member.topFrac==0 ? 0 : 1)+(member.qFrac==.0001f ? 0 : 1));
				final long recovered=tp[a][b]-penalty;
				rows.append(policy.name).tab().append(recovered).append("\t0\t").append(100-recovered)
					.append("\t100\t").append(100-recovered).append("\t0\t").append(cost[a][b]+penalty);
				if(modelBases){rows.tab().append(100*(a+1)+(b==1 ? 0 : 100)+penalty);}
				rows.nl();
			}
			try(Output output=new Output(prefix.getPath()+".summary.tsv")){output.writer.print(rows);}
		}
		int calls;
		final boolean modelBases;
	}

	private static String[] with(final String[] arguments, final String key, final String value){
		for(int i=0; i<arguments.length; i++){
			if(arguments[i].startsWith(key+"=")){
				final String[] result=arguments.clone();
				result[i]=key+"="+value;
				return result;
			}
		}
		final String[] result=Arrays.copyOf(arguments, arguments.length+1);
		result[arguments.length]=key+"="+value;
		return result;
	}
	private static void write(final File path, final String text){
		try(Output output=new Output(path.getPath())){output.writer.print(new ByteBuilder(text.getBytes(StandardCharsets.UTF_8)));}
	}
	private static String read(final File path){
		final ByteBuilder text=new ByteBuilder();
		try(Input input=new Input(path.getPath())){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){text.append(line).nl();}
		}
		return new String(text.toBytes(), StandardCharsets.UTF_8);
	}
	private static NcrnaPolicySearch search(){
		return new NcrnaPolicySearch(start(), new Axis[]{new Axis(0, 0, new double[]{1, 2, 3}),
			new Axis(0, 1, new double[]{1, 2}), new Axis(0, 2, new double[]{0, .25}),
			new Axis(0, 3, new double[]{0, .0001}), new Axis(0, 4, new double[]{0, 1, 2})});
	}
	private static Measurement m(final long tp, final long cost){return new Measurement(tp, 0, cost, 0);}
	private static void rejected(final Runnable action, final String why){
		try{action.run();}
		catch(IllegalArgumentException expected){return;}
		throw new AssertionError(why);
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
