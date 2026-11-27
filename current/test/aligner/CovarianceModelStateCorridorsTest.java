package test.aligner;

import java.util.Arrays;
import aligner.CovarianceModel;
import aligner.CovarianceModelStateCorridors;
import aligner.CovarianceModelConsensusMap;
import aligner.CovarianceModelPlacementBands;
import aligner.CovarianceModelLengthBands;
import aligner.CovarianceModelBandedCyk;
import map.ObjectMap;
import idaligner.AlignmentStats;

/** Literal interval arithmetic independent of fitting, CM scores or path exports.
 * In particular, an interval outside the window must stay empty, and signed
 * offsets need not be ordered before adding their different raw anchors.
 * @author Brian Bushnell, Raiden
 */
public final class CovarianceModelStateCorridorsTest{

	public static void main(final String[] args){
		check(args.length==0, "Literal corridor arithmetic uses no input files or production models");
		final CovarianceModel model=new CovarianceModel("literal", "literal", 2, 10, 3, 3,
			new ObjectMap<String, String>(String.class, String.class), new float[]{.25f, .25f, .25f, .25f});
		model.type[0]=CovarianceModel.S;
		model.type[1]=CovarianceModel.MP;
		model.type[2]=CovarianceModel.E;
		final short[] offsets={1, -2, -1, 2, 0, 0, 0, 0, 0, 0, 0, 0};
		final boolean[] observed={true, false, true};
		final CovarianceModelStateCorridors profile=new CovarianceModelStateCorridors(model, offsets, observed);
		offsets[0]=99;
		observed[0]=false;
		check(profile.observed(0) && profile.offset(0, 0)==1 && !profile.observed(1),
			"The immutable resource must own its offsets and support flags after construction");
		final int[] bounds={2, 5, 5, 6};
		check(profile.project(0, 10, bounds, 32) && Arrays.equals(bounds, new int[]{3, 3, 4, 8}),
			"Separate low/high residuals turn the literal raw intervals into[3,3]and[4,8]");
		check(profile.contains(0, 10, bounds, 3, 4) && !profile.contains(0, 10, bounds, 2, 4),
			"Containment must use the actual projected endpoints");
		set(bounds, 2, 3, 5, 6);
		check(!profile.project(0, 10, bounds, 32) && !profile.contains(0, 10, bounds, 2, 5),
			"An inverted learned interval is empty rather than reordered or widened");
		set(bounds, 2, 2, 5, 5);
		check(profile.project(1, 10, bounds, 32) && Arrays.equals(bounds, new int[]{0, 10, 0, 10}),
			"An unseen state uses the complete uniform32fallback on this ten-base window");
		check(!profile.contains(1, 10, bounds, 8, 7), "The triangular DP domain forbids negative consumed lengths");
		set(bounds, 4, 6, 4, 6);
		check(profile.project(2, 10, bounds, 32) && profile.contains(2, 10, bounds, 5, 5)
			&& !profile.contains(2, 10, bounds, 4, 5), "END states consume zero residues even inside wide endpoint ranges");
		outside(model, 20);
		outside(model, -20);
		final CovarianceModelStateCorridors clipped=uniformOffsets(model, -2, 2, -2, 2);
		set(bounds, 0, 0, 9, 10);
		check(clipped.project(0, 10, bounds, 32) && Arrays.equals(bounds, new int[]{0, 2, 7, 10}),
			"A partially outside interval intersects the natural window edge");
		set(bounds, Integer.MAX_VALUE-1, Integer.MAX_VALUE-1, Integer.MAX_VALUE, Integer.MAX_VALUE);
		check(profile.project(1, Integer.MAX_VALUE, bounds, Integer.MAX_VALUE)
			&& Arrays.equals(bounds, new int[]{0, Integer.MAX_VALUE, 0, Integer.MAX_VALUE}),
			"Uniform expansion uses long arithmetic before intersecting integer coordinate limits");
		translatedAnchors(model);
		contextWindowSearches();
		searches();
		bifurcationDiagnostic();
		unionLayouts();
		mixedUnionPath();
		CovarianceModelSearchProfileTest.main(new String[0]);
		System.err.println("CM_STATE_CORRIDOR_INTERVAL_FIXTURE_PASS");
	}

	/** Hand arithmetic crosses both search edges; no fit or alignment is its oracle. */
	private static void translatedAnchors(final CovarianceModel model){
		final int[] bounds={-5, -5, 15, 15};
		final CovarianceModelStateCorridors inward=uniformOffsets(model, 8, 8, -8, -8);
		check(inward.projectUnclipped(0, 10, bounds, 32) && Arrays.equals(bounds, new int[]{3, 3, 7, 7}),
			"Raw anchors -5/15 plus residuals +8/-8 give3/7; preclipping to0/10 would wrongly give8/2");
		check(inward.contains(0, 10, bounds, 3, 7) && !inward.contains(0, 10, bounds, 8, 8),
			"Translated residuals must retain the literal target cell rather than a clipped-anchor substitute");
		set(bounds, -5, -5, 5, 5);
		check(!uniformOffsets(model, 2, 2, 0, 0).projectUnclipped(0, 10, bounds, 32),
			"Raw-5 plus2 remains outside; clamping the raw anchor first must not invent a start at2");
		set(bounds, 5, 5, 15, 15);
		check(!uniformOffsets(model, 0, 0, -2, -2).projectUnclipped(0, 10, bounds, 32),
			"Raw15 minus2 remains outside; clamping first must not invent an end at8");
		set(bounds, -5, 2, 8, 15);
		check(uniformOffsets(model, 1, 1, -1, -1).projectUnclipped(0, 10, bounds, 32)
			&& Arrays.equals(bounds, new int[]{0, 3, 7, 10}),
			"Partly outside translated intervals intersect the physical search domain only after their residuals");
		set(bounds, Integer.MIN_VALUE, Integer.MIN_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE);
		final CovarianceModelStateCorridors unseen=new CovarianceModelStateCorridors(model,
			new short[model.states()*4], new boolean[model.states()]);
		check(!unseen.projectUnclipped(0, 10, bounds, Integer.MAX_VALUE),
			"Even maximum uniform expansion leaves the minimum-integer start at-1; long arithmetic must preserve that empty intersection");
		set(bounds, -1, 0, 5, 6);
		boolean rejected=false;
		try{inward.project(0, 10, bounds, 32);}catch(final IllegalArgumentException expected){rejected=true;}
		check(rejected, "The original same-window API must continue rejecting outside raw anchors");
		set(bounds, 2, 1, 5, 6);
		rejected=false;
		try{inward.projectUnclipped(0, 10, bounds, 32);}catch(final IllegalArgumentException expected){rejected=true;}
		check(rejected, "A distinct-context API permits outside coordinates, never inverted raw intervals");
		System.err.println("CM_STATE_CORRIDOR_TRANSLATED_ANCHOR_FIXTURE_PASS");
	}

	/** A literal full-query trace is outside the smaller search window. The
	 * independently specified residuals recover AC; clipping first cannot do so.
	 */
	private static void contextWindowSearches(){
		final CovarianceModel model=CovarianceModelCykTest.pair();
		model.nodeType[0]=CovarianceModel.ROOT;
		model.nodeType[1]=CovarianceModel.MATP;
		model.nodeType[2]=CovarianceModel.END;
		model.consensusRight[0]=3;
		model.consensusLeft[2]=1;
		model.consensusRight[2]=2;
		Arrays.fill(model.emissionScore[1], Float.NEGATIVE_INFINITY);
		model.emissionScore[1][1]=9;
		final byte[] consensus={'A', 'C'}, rf={'x', 'x'};
		final CovarianceModelConsensusMap map=new CovarianceModelConsensusMap(model, consensus, consensus, rf);
		final AlignmentStats trace=new AlignmentStats(true);
		trace.qLen=2;
		trace.rLen=20;
		trace.rStart=15;
		trace.rStop=16;
		trace.setFromMatchString(new byte[]{'m', 'm'});
		final byte[] window={'T', 'T', 'A', 'C', 'G', 'G'};
		final CovarianceModelPlacementBands raw=map.composeWindowUnclipped(trace, 100, 100, 105, 1);
		final int[] bounds=new int[4];
		raw.stateBounds(0, bounds);
		check(Arrays.equals(bounds, new int[]{15, 15, 17, 17}),
			"The separate-context mapper retains raw15/17 on a six-base search domain before signed residuals");
		final CovarianceModelPlacementBands old=map.composeWindow(trace, 100, 100, 105, 1);
		old.stateBounds(0, bounds);
		check(Arrays.equals(bounds, new int[]{6, 6, 6, 6}),
			"The original same-window mapper preserves its existing edge-clipping contract");
		final CovarianceModelStateCorridors shifted=uniformOffsets(model, -13, -13, -13, -13);
		final CovarianceModelBandedCyk scorer=new CovarianceModelBandedCyk(model,
			new CovarianceModelLengthBands(model, 2), 10000, 10000);
		for(int repeat=0; repeat<3; repeat++){
			final CovarianceModelBandedCyk.EndpointResult result=scorer.searchOverlapping(window, 6, 2, 5, raw, shifted);
			check(result.score==9 && result.targetFrom==3 && result.targetTo==4,
				"Raw15/17 minus13 retains the independently specified AC at3..4 in actual learned DP layouts");
			check(!Float.isFinite(scorer.searchOverlapping(window, 6, 2, 5, old, shifted).score),
				"The old clipped raw coordinates with these residuals are a deliberate impossible counterexample");
			check(scorer.searchOverlapping(window, 6, 2, 5, raw, shifted, 1).score==9,
				"An empty outside glocal union component cannot discard a finite learned component");
			check(!Float.isFinite(scorer.searchOverlapping(window, 6, 2, 5, raw,
				uniformOffsets(model, 0, 0, 0, 0), 1).score),
				"Outside learned and glocal components together stay empty without fabricated edge cells");
			check(!Float.isFinite(scorer.searchOverlapping(window, 6, 2, 5, raw).score),
				"A small uniform radius cannot silently clamp an outside raw anchor onto the search edge");
		}
		final CovarianceModelPlacementBands translated=map.composeWindowUnclipped(trace, 100, 113, 118, 1);
		translated.stateBounds(0, bounds);
		check(Arrays.equals(bounds, new int[]{2, 2, 4, 4}),
			"Trace origin100 and search origin113 translate raw15/17 into2/4 exactly once");
		check(scorer.searchOverlapping(window, 6, 2, 5, translated, uniformOffsets(model, 0, 0, 0, 0)).score==9,
			"A translated in-window anchor needs zero residuals to recover the same literal AC grammar");
		check(scorer.searchOverlapping(window, 6, 2, 5, translated,
			uniformOffsets(model, 20, 20, 20, 20), 1).score==9,
			"A finite glocal component survives when the translated learned component is empty");
		System.err.println("CM_STATE_CORRIDOR_CONTEXT_WINDOW_SEARCH_FIXTURE_PASS");
	}

	/** Independent three-state grammar has only the AC pair, at bases3..4.
	 * Exercise actual layouts and recycled score/contact decks, including an
	 * empty END child under a nonempty root; empty-root tests alone miss that.
	 */
	private static void searches(){
		final CovarianceModel model=CovarianceModelCykTest.pair();
		model.nodeType[0]=CovarianceModel.ROOT;
		model.nodeType[1]=CovarianceModel.MATP;
		model.nodeType[2]=CovarianceModel.END;
		model.consensusRight[0]=3;
		model.consensusLeft[2]=1;
		model.consensusRight[2]=2;
		Arrays.fill(model.emissionScore[1], Float.NEGATIVE_INFINITY);
		model.emissionScore[1][1]=9;
		final byte[] window={'T', 'T', 'A', 'C', 'G', 'G'};
		final CovarianceModelPlacementBands placement=CovarianceModelPlacementBands.window(model, 6,
			new int[]{2, 3, 4}, new int[]{2, 3, 4}, 1);
		final CovarianceModelBandedCyk scorer=new CovarianceModelBandedCyk(model,
			new CovarianceModelLengthBands(model, 2), 10000, 10000);
		final CovarianceModelStateCorridors exact=uniformOffsets(model, 0, 0, 0, 0);
		final short[] endOutside=new short[12];
		Arrays.fill(endOutside, 8, 12, (short)20);
		final CovarianceModelStateCorridors emptyChild=new CovarianceModelStateCorridors(model, endOutside,
			new boolean[]{true, true, true});
		final CovarianceModelStateCorridors unseen=new CovarianceModelStateCorridors(model, endOutside,
			new boolean[]{true, true, false});
		for(int repeat=0; repeat<3; repeat++){
			final CovarianceModelBandedCyk.EndpointResult tight=scorer.searchOverlapping(window, 6, 2, 5, placement, exact);
			check(tight.score==9 && tight.targetFrom==3 && tight.targetTo==4 && tight.contacts!=0,
				"Literal learned zero-residual bounds retain AC at3..4 with interior contact");
			check(!Float.isFinite(scorer.searchOverlapping(window, 6, 2, 5, placement, emptyChild).score),
				"An empty END child must make the real grammar impossible, without a fabricated edge cell or stale score");
			check(scorer.searchOverlapping(window, 6, 2, 5, placement, unseen).score==9,
				"Unseen END states ignore stored offsets and use the uniform-radius fallback");
			check(!Float.isFinite(scorer.searchOverlapping(window, 6, 2, 5, placement,
				uniformOffsets(model, 20, 20, 20, 20)).score), "An empty learned root returns impossible instead of invoking exact CYK");
			final CovarianceModelBandedCyk.EndpointResult legacy=scorer.searchOverlapping(window, 6, 2, 5, placement);
			check(legacy.score==9 && legacy.targetFrom==3 && legacy.targetTo==4 && legacy.contacts==0,
				"The old uniform API must recover its original score, endpoints and contacts after learned/empty deck reuse");
			final CovarianceModelBandedCyk.DiagnosticResult traced=diagnosticParity(scorer, window, placement, exact);
			check(traced.states.size==3 && traced.states.get(0)==0 && traced.states.get(1)==1
				&& traced.states.get(2)==2 && traced.starts.get(0)==2 && traced.ends.get(0)==4,
				"Diagnostic traceback must follow the selected AC root in consumed-prefix coordinates, including its END cell");
			check(traced.emptyStates==0 && traced.emptyObservedStates==0,
				"The literal finite grammar allocates a nonempty layout for every observed state");
			final CovarianceModelBandedCyk.DiagnosticResult impossible=diagnosticParity(scorer, window, placement, emptyChild);
			check(impossible.states.size==0 && impossible.emptyStates==1 && impossible.emptyObservedStates==1,
				"An empty END layout makes the whole search impossible; no diagnostic path may be fabricated");
			diagnosticParity(scorer, window, placement, unseen);
			diagnosticParity(scorer, window, placement, uniformOffsets(model, 20, 20, 20, 20));
		}
		final byte[] tied={'T', 'T', 'A', 'C', 'A', 'C', 'G', 'G'};
		final CovarianceModelPlacementBands tiedPlacement=CovarianceModelPlacementBands.window(model, tied.length,
			new int[]{2, 3, 4}, new int[]{2, 3, 4}, 1);
		final CovarianceModelBandedCyk.DiagnosticResult tie=diagnosticParity(scorer, tied, tiedPlacement,
			uniformOffsets(model, -2, 4, -2, 4));
		check(tie.endpoint.score==9 && tie.endpoint.targetFrom==3 && tie.endpoint.targetTo==4,
			"Equal scoring AC hits must retain the ordinary root scan's earliest end rather than diagnostic traversal order");
		boolean rejected=false;
		try{scorer.searchOverlapping(window, 6, 2, 5, null, exact);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Learned search requires independently supplied anchors and cannot silently request exact fallback");
	}

	private static CovarianceModelBandedCyk.DiagnosticResult diagnosticParity(final CovarianceModelBandedCyk scorer,
			final byte[] window, final CovarianceModelPlacementBands placement, final CovarianceModelStateCorridors profile){
		final CovarianceModelBandedCyk.EndpointResult ordinary=scorer.searchOverlapping(window, window.length, 2, 5, placement, profile);
		final CovarianceModelBandedCyk.DiagnosticResult result=scorer.diagnoseOverlapping(window, window.length, 2, 5, placement, profile);
		final CovarianceModelBandedCyk.EndpointResult traced=result.endpoint;
		check(ordinary.score==traced.score && ordinary.targetFrom==traced.targetFrom && ordinary.targetTo==traced.targetTo
			&& ordinary.contacts==traced.contacts && ordinary.cells==traced.cells
			&& ordinary.unrestrictedScore==traced.unrestrictedScore && ordinary.unrestrictedFrom==traced.unrestrictedFrom
			&& ordinary.unrestrictedTo==traced.unrestrictedTo,
			"Trace-enabled diagnosis must preserve ordinary finite/impossible scores, root tie order, contact OR and allocated cells");
		int contacts=0;
		for(int step=0; step<result.states.size; step++){
			check(result.starts.get(step)>=0 && result.starts.get(step)<=result.ends.get(step)
				&& result.ends.get(step)<=window.length, "Exported diagnostic cells must be physical consumed-prefix intervals");
			contacts|=result.contactMasks.get(step);
		}
		check(contacts==traced.contacts, "The exported path's contact masks must explain the endpoint scorer's full-path contact OR");
		return result;
	}

	/** A B state needs both pinned child scores when its split is recovered after
	 * the fill. The independent grammar emits A in the left branch and T right.
	 */
	private static void bifurcationDiagnostic(){
		final CovarianceModel model=CovarianceModelCykTest.branch();
		model.nodeType[0]=CovarianceModel.ROOT;
		model.nodeType[1]=CovarianceModel.BIF;
		model.nodeType[2]=CovarianceModel.MATL;
		model.nodeType[3]=CovarianceModel.END;
		model.nodeType[4]=CovarianceModel.MATR;
		model.nodeType[5]=CovarianceModel.END;
		model.consensusRight[0]=3;
		model.consensusRight[1]=3;
		for(int node=2; node<6; node++){
			model.consensusLeft[node]=1;
			model.consensusRight[node]=2;
		}
		final byte[] window={'T', 'T', 'A', 'T', 'G', 'G'};
		final CovarianceModelPlacementBands placement=CovarianceModelPlacementBands.window(model, window.length,
			new int[]{2, 3, 4}, new int[]{2, 3, 4}, 1);
		final CovarianceModelBandedCyk scorer=new CovarianceModelBandedCyk(model,
			new CovarianceModelLengthBands(model, 2), 10000, 10000);
		final CovarianceModelBandedCyk.DiagnosticResult result=diagnosticParity(scorer, window, placement,
			uniformOffsets(model, 0, 0, 0, 0));
		check(result.endpoint.score==5 && result.endpoint.targetFrom==3 && result.endpoint.targetTo==4
			&& result.states.size==6 && result.states.get(1)==1 && result.states.get(2)==2
			&& result.states.get(3)==3 && result.states.get(4)==4 && result.states.get(5)==5,
			"Diagnostic B traceback must retain both A/T children and their END cells in left-first order");
	}

	private static void outside(final CovarianceModel model, final int shift){
		assert(shift==20 || shift==-20) : "Both external sides of the same ten-base window are independently enumerated";
		final CovarianceModelStateCorridors profile=uniformOffsets(model, shift, shift, shift, shift);
		final int[] bounds={2, 3, 4, 5};
		check(!profile.project(0, 10, bounds, 32), "A wholly outside learned interval must not become a point on the nearest edge");
		for(int start=0; start<=10; start++){
			for(int end=start; end<=10; end++){
				check(!profile.contains(0, 10, bounds, start, end), "An empty projection cannot admit any actual window cell");
			}
		}
	}

	/** Enumerate independent literal cell sets, including disjoint row segments,
	 * cross-corner exclusions, overlapping components and separated END rows.
	 */
	private static void unionLayouts(){
		final int[][] rectangles={{1, 2, 8, 10}, {5, 6, 8, 10}, {1, 2, 6, 7},
			{4, 5, 9, 10}, {0, 10, 0, 10}, {1, 2, 1, 2}, {5, 6, 5, 6}, {8, 9, 1, 2}};
		for(final int[] a:rectangles){
			for(final int[] b:rectangles){
				for(final int maximum:new int[]{0, 3, 10}){
					final CovarianceModelLengthBands.Layout union=CovarianceModelLengthBands.Layout.union(10, 0, maximum, a, b);
					int expected=0, iterated=0;
					for(int end=0; end<=10; end++){
						for(int size=0; size<=10; size++){
							final int start=end-size;
							final boolean kept=start>=0 && size<=maximum && (inside(start, end, a) || inside(start, end, b));
							check(union.index(end, size)==(kept ? expected : -1),
								"Packed union addresses must exactly enumerate the literal two-rectangle cell set without hull holes");
							if(kept){expected++;}
						}
						for(int size=union.lower(end); size<=union.upper(end); size=union.next(end, size)){
							check(union.index(end, size)==iterated++, "Row iteration must visit each packed union address once in score order");
						}
					}
					check(expected==union.cells && iterated==expected,
						"Union storage must allocate exactly the independently counted retained cells, including END row gaps");
				}
			}
		}
		final CovarianceModelLengthBands.Layout interior=CovarianceModelLengthBands.Layout.union(10, 0, 10,
			new int[]{1, 3, 5, 9}, new int[]{3, 5, 5, 9});
		check(interior.contactMask(3, 7)==0 && interior.contactMask(1, 7)==1 && interior.contactMask(5, 7)==2,
			"Covered component faces are internal; only the outer union faces retain their original direction bits");
		final CovarianceModelLengthBands.Layout corner=CovarianceModelLengthBands.Layout.union(10, 0, 10,
			new int[]{1, 3, 5, 9}, new int[]{1, 5, 5, 7});
		check(corner.index(8, 4)<0 && corner.contactMask(3, 7)==10 && corner.contactMask(3, 6)==0,
			"The missing upper-right cell exposes a concave start-high/end-high corner despite both axial neighbors being allowed");
		final CovarianceModelLengthBands.Layout natural=CovarianceModelLengthBands.Layout.union(10, 0, 10,
			new int[]{0, 10, 0, 10}, new int[]{0, 10, 0, 10});
		check(natural.contactMask(0, 10)==0 && natural.contactMask(5, 5)==0,
			"Natural window and triangular zero-length boundaries never acquire imposed-face contact bits");
		final int[] identical={1, 4, 5, 9};
		final CovarianceModelLengthBands.Layout rectangle=new CovarianceModelLengthBands.Layout(10, 0, 10, 1, 4, 5, 9);
		final CovarianceModelLengthBands.Layout same=CovarianceModelLengthBands.Layout.union(10, 0, 10, identical, identical);
		for(int end=5; end<=9; end++){
			for(int start=1; start<=4; start++){
				check(same.contactMask(start, end)==rectangle.contactMask(start, end),
					"Coincident union components must retain every legacy rectangle contact bit exactly");
			}
		}
	}
	private static boolean inside(final int start, final int end, final int[] bounds){
		assert(bounds.length==4) : "The independent fixture uses literal endpoint rectangles";
		return start>=bounds[0] && start<=bounds[1] && end>=bounds[2] && end<=bounds[3];
	}

	/** Neither complete region has a finite AC parse: its emitting root needs
	 * learned cells and its END child needs the glocal component. Only one union
	 * recurrence may combine those states; best-of-two whole searches cannot.
	 */
	private static void mixedUnionPath(){
		final CovarianceModel model=CovarianceModelCykTest.pair();
		model.nodeType[0]=CovarianceModel.ROOT;
		model.nodeType[1]=CovarianceModel.MATP;
		model.nodeType[2]=CovarianceModel.END;
		model.consensusRight[0]=3;
		model.consensusLeft[2]=1;
		model.consensusRight[2]=2;
		Arrays.fill(model.emissionScore[1], Float.NEGATIVE_INFINITY);
		model.emissionScore[1][1]=9;
		final byte[] window={'T', 'T', 'A', 'C', 'G', 'G', 'G', 'G', 'G', 'G'};
		final CovarianceModelPlacementBands placement=CovarianceModelPlacementBands.window(model, window.length,
			new int[]{4, 5, 8}, new int[]{4, 5, 8}, 2);
		final short[] offsets={-2, -2, -4, -4, -2, -2, -4, -4, 20, 20, 20, 20};
		final CovarianceModelStateCorridors learned=new CovarianceModelStateCorridors(model, offsets, new boolean[]{true, true, true});
		final CovarianceModelBandedCyk scorer=new CovarianceModelBandedCyk(model,
			new CovarianceModelLengthBands(model, 2), 10000, 10000);
		for(int repeat=0; repeat<3; repeat++){
			check(!Float.isFinite(scorer.searchOverlapping(window, 10, 3, 4, placement, learned).score)
				&& !Float.isFinite(scorer.searchOverlapping(window, 10, 3, 4, placement).score),
				"The two independent whole-region searches must both fail in this mixed-path counterexample");
			final CovarianceModelBandedCyk.EndpointResult joined=scorer.searchOverlapping(window, 10, 3, 4, placement, learned, 2);
			check(joined.score==9 && joined.targetFrom==3 && joined.targetTo==4 && joined.contacts!=0,
				"One exact union fill must retain the mixed learned-root/glocal-END path and its exposed imposed faces");
		}
		final short[] endGap={-2, -2, -4, -4, -2, -2, -4, -4, -5, -5, -5, -5};
		final CovarianceModelStateCorridors missingEnd=new CovarianceModelStateCorridors(model, endGap, new boolean[]{true, true, true});
		check(!Float.isFinite(scorer.searchOverlapping(window, 10, 3, 4, placement, missingEnd, 1).score),
			"The required END cell at3 lies between retained END0 and END4..6; zero initialization must never fabricate that gap cell");
		final byte[] rowGap={'G', 'G', 'G', 'G', 'G', 'G', 'A', 'C', 'G', 'G'};
		final CovarianceModelPlacementBands right=CovarianceModelPlacementBands.window(model, rowGap.length,
			new int[]{6, 7, 8}, new int[]{6, 7, 8}, 1);
		final short[] separated={-4, -4, 0, 0, -4, -4, 0, 0, 20, 20, 20, 20};
		final CovarianceModelStateCorridors gap=new CovarianceModelStateCorridors(model, separated, new boolean[]{true, true, true});
		final CovarianceModelBandedCyk.EndpointResult compact=scorer.searchOverlapping(rowGap, 10, 7, 8, right, gap, 1);
		check(compact.score==9 && compact.targetFrom==7 && compact.targetTo==8,
			"A real emitting/root row with retained lengths1..3 and6 must skip lengths4..5 and retain the correct packed score address");
	}
	private static CovarianceModelStateCorridors uniformOffsets(final CovarianceModel model, final int a, final int b, final int c, final int d){
		assert(Math.min(Math.min(a, b), Math.min(c, d))>=Short.MIN_VALUE
			&& Math.max(Math.max(a, b), Math.max(c, d))<=Short.MAX_VALUE) : "Literal offsets must fit the declared signed16-bit resource";
		final short[] values=new short[model.states()*4];
		for(int state=0; state<model.states(); state++){
			values[state*4]=(short)a;
			values[state*4+1]=(short)b;
			values[state*4+2]=(short)c;
			values[state*4+3]=(short)d;
		}
		final boolean[] observed=new boolean[model.states()];
		Arrays.fill(observed, true);
		return new CovarianceModelStateCorridors(model, values, observed);
	}
	private static void set(final int[] bounds, final int a, final int b, final int c, final int d){
		assert(bounds.length==4) : "Each independent projection uses one caller-owned four-value scratch buffer";
		bounds[0]=a;
		bounds[1]=b;
		bounds[2]=c;
		bounds[3]=d;
	}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
}
