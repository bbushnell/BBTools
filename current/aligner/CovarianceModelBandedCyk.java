package aligner;

import java.util.Arrays;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import structures.IntList;
import static aligner.CovarianceModel.*;

/** Full-input CYK within explicit bands; optional configured local begins/ends.
 * Global model scoring remains the default. Ordinary traceback choices
 * use one unsigned byte; B splits are recomputed from pinned child decks so a
 * split larger than255 is never truncated. Scores outside bands are impossible.
 * Banding may change the optimum; compare to the exact unbanded scorer.
 * Instances are thread-confined: endpoint scoring reuses mutable deck buffers.
 * @author Brian Bushnell, Raiden
 */
public final class CovarianceModelBandedCyk {

	public CovarianceModelBandedCyk(CovarianceModel model_, CovarianceModelLengthBands bands_, long maxScoreCells_, long maxTraceCells_){
		this(model_, bands_, maxScoreCells_, maxTraceCells_, false);
	}
	public CovarianceModelBandedCyk(CovarianceModel model_, CovarianceModelLengthBands bands_, long maxScoreCells_, long maxTraceCells_, boolean local){
		require(model_!=null && bands_!=null && bands_.model==model_, "Bands must belong to this exact immutable model, not another model with the same dimensions");
		require(maxScoreCells_>0 && maxTraceCells_>0, "Explicit positive score/trace limits protect band allocation");
		model=model_;bands=bands_;maxScoreCells=maxScoreCells_;maxTraceCells=maxTraceCells_;
		configuration=new CovarianceModelConfiguration(model, local);
		emissions=CovarianceModelAlphabet.expand(model);
		parents=new int[model.states()];pinned=new boolean[model.states()];
		int mostChildren=0;
		for(int v=0; v<model.states(); v++){
			final byte type=model.type[v];require(type>=D && type<=B, "Parsed model states are ordinary; virtual EL is configured separately");
			if(type==B){edge(v, model.childFirst[v]);edge(v, model.childCount[v]);pinned[model.childFirst[v]]=pinned[model.childCount[v]]=true;}
			else if(type==E){require(model.childCount[v]==0, "END states cannot have child edges");}
			else{require(model.childCount[v]>0 && model.childCount[v]<=(local ? 254 : 255) && model.transitionScore[v].length==model.childCount[v],
				"Byte choices reserve255 for unset and254 for local ends in local mode, state="+v);
				mostChildren=Math.max(mostChildren, model.childCount[v]);
				for(int c=0; c<model.childCount[v]; c++){edge(v, model.childFirst[v]+c);}}
		}
		childOffset=new int[mostChildren];
		childLow=new int[mostChildren];
		childHigh=new int[mostChildren];
	}

	public Result align(byte[] sequence){
		require(sequence!=null, "Missing sequence is not an empty input");
		final CovarianceModelLengthBands.Layout[] layout=new CovarianceModelLengthBands.Layout[model.states()];
		for(int v=0; v<layout.length; v++){layout[v]=bands.layout(v, sequence.length);}
		// Stored QDB tails were computed without local ends; their probability
		// bound does not describe the modified local distribution.
		return align(sequence, layout, configuration.local ? Double.NaN : bands.tailProbability);
	}
	/** An explicit candidate layout replaces QDB bounds; it must not silently
	 * inherit the model-tail probability, which says nothing about anchors. */
	Result align(byte[] sequence, CovarianceModelLengthBands.Layout[] layout, double tailProbability){
		final Fill fill=fill(sequence, layout, false);
		final int length=sequence.length;
		final int root=layout[0].index(length, length);
		final float score=root<0 ? Float.NEGATIVE_INFINITY : fill.storage.active[0][root];
		final CovarianceModelCyk.Result alignment=new CovarianceModelCyk.Result(score, fill.visited,
			length, 1, length, new IntList(), new IntList());
		if(Float.isFinite(score)){
			trace(alignment, fill.choices, layout, fill.storage.active, fill.bases, fill.localBegin);
		}
		return new Result(alignment, fill.traceCells, fill.storage.peak, tailProbability);
	}

	/** Best associated nonempty interval, without CM traceback. Null placement
	 * means an exact full-layout fallback, never an implicit QDB restriction.
	 * Equal scores retain earliest end, then shortest span, as in exact search. */
	public EndpointResult searchOverlapping(byte[] sequence, int maxWindow,
			int candidateFrom, int candidateTo, CovarianceModelPlacementBands placement){
		return searchOverlapping(sequence, maxWindow, candidateFrom, candidateTo, placement, null);
	}
	/** Optional learned bounds replace the uniform state intervals. The caller
	 * owns contact/impossibility fallback; this method never invokes exact CYK. */
	public EndpointResult searchOverlapping(byte[] sequence, int maxWindow,
			int candidateFrom, int candidateTo, CovarianceModelPlacementBands placement,
			CovarianceModelStateCorridors corridors){
		return searchOverlapping(sequence, maxWindow, candidateFrom, candidateTo, placement, corridors, -1);
	}
	/** Optional exact union with a second glocal region in one recurrence. This
	 * retains mixed paths and never admits cross-corner cells from a hull. */
	public EndpointResult searchOverlapping(byte[] sequence, int maxWindow,
			int candidateFrom, int candidateTo, CovarianceModelPlacementBands placement,
			CovarianceModelStateCorridors corridors, final int unionRadius){
		return searchOverlapping(sequence, maxWindow, candidateFrom, candidateTo, placement, corridors, unionRadius, null);
	}
	/** Optional worker-owned accounting leaves the ordinary recurrence unchanged.
	 * Cell geometry is counted after the timed search, without per-cell clocks. */
	public EndpointResult searchOverlapping(byte[] sequence, int maxWindow,
			int candidateFrom, int candidateTo, CovarianceModelPlacementBands placement,
			CovarianceModelStateCorridors corridors, final int unionRadius, final SearchProfile profile){
		require(sequence!=null && maxWindow>0 && candidateFrom<=candidateTo && !configuration.local,
			"Endpoint search needs a window, positive hit bound, candidate interval and global CM configuration");
		require(placement==null || placement.model==model && placement.length==sequence.length,
			"Placement must bind this exact CM and complete search window");
		require(corridors==null || placement!=null && corridors.belongsTo(model),
			"Learned bounds need independent placement anchors on this exact CM; absent anchors cannot request an exact fallback");
		require(unionRadius==-1 || corridors!=null && unionRadius>=0,
			"A union is explicit and needs learned bounds; legacy searches use no second component");
		final long begin=profile==null ? 0 : profile.now();
		final CovarianceModelLengthBands.Layout[] layout;
		if(placement!=null){layout=placement.layouts(corridors, unionRadius);}
		else{
			layout=new CovarianceModelLengthBands.Layout[model.states()];
			for(int v=0; v<layout.length; v++){
				layout[v]=new CovarianceModelLengthBands.Layout(sequence.length, 0, model.type[v]==E ? 0 : sequence.length);
			}
		}
		final long composed=profile==null ? 0 : profile.now();
		final Fill fill=fill(sequence, layout, true, profile==null ? null : profile.kernel);
		final long filled=profile==null ? 0 : profile.now();
		final EndpointResult result=selectEndpoints(fill, layout[0], maxWindow, candidateFrom, candidateTo);
		if(profile!=null){profile.record(model, layout, result, begin, composed, filled);}
		return result;
	}

	/** Explicit endpoint-search measurement, confined to one executing worker.
	 * representedCells counts initialized deck slots; updatedCells counts END
	 * initialization, bifurcation and ordinary recurrence cells, not transitions
	 * or bifurcation split iterations. Geometry counting is timed separately. */
	public static final class SearchProfile{
		public SearchProfile(){
			this(false);
		}
		/** Detailed fill accounting is opt-in; ordinary search clocks stay unchanged. */
		public SearchProfile(final boolean detailed){
			clock=ManagementFactory.getThreadMXBean();
			require(clock.isCurrentThreadCpuTimeSupported(), "Search profiling needs a current-thread CPU clock, never wall time");
			if(!clock.isThreadCpuTimeEnabled()){clock.setThreadCpuTimeEnabled(true);}
			kernel=detailed ? new KernelProfile(clock) : null;
		}
		private long now(){
			final long value=clock.getCurrentThreadCpuTime();
			require(value>=0, "Enabled search profiling requires a live current-thread CPU clock");
			return value;
		}
		private void record(final CovarianceModel model, final CovarianceModelLengthBands.Layout[] layouts,
				final EndpointResult result, final long begin, final long composed, final long filled){
			final long selected=now();
			assert(begin<=composed && composed<=filled && filled<=selected) :
				"Consecutive measurements on the same scoring thread must be monotonic";
			searches++;
			layoutCpu+=composed-begin;
			fillCpu+=filled-composed;
			selectionCpu+=selected-filled;
			lastResult=result;
			representedCells+=result.cells;
			long updates=0;
			for(int state=0; state<layouts.length; state++){
				final CovarianceModelLengthBands.Layout layout=layouts[state];
				long cells=layout.cells;
				final int delta=CovarianceModelCyk.delta(model.type[state]);
				// fill() skips d<delta in ordinary states. The triangular domain
				// already excludes d>j; END and bifurcation visit all their cells.
				if(model.type[state]!=E && model.type[state]!=B && delta>0){
					for(int end=layout.firstEnd; end<=layout.lastEnd && layout.width>0; end++){
						for(int length=0; length<delta; length++){
							if(layout.index(end, length)>=0){cells--;}
						}
					}
				}
				assert(cells>=0 && cells<=layout.cells) : "Updated cells are a subset of the actual allocated state domain";
				updates+=cells;
				if(kernel!=null){kernel.updatedCells[model.type[state]]+=cells;}
			}
			updatedCells+=updates;
			geometryCpu+=now()-selected;
		}
		private final ThreadMXBean clock;
		public final KernelProfile kernel;
		/** Last actual search; the refiner clears this at each proposal boundary. */
		public EndpointResult lastResult;
		public long searches, layoutCpu, fillCpu, selectionCpu, geometryCpu, representedCells, updatedCells;
	}

	/** State-boundary clocks only. Split trials count proposed k values, including
	 * those rejected by child geometry; they are not extra DP cells. END's second
	 * zero initialization belongs to END recurrence, after the initial clear. */
	public static final class KernelProfile{
		private KernelProfile(final ThreadMXBean clock_){clock=clock_;}
		private long now(){
			final long time=clock.getCurrentThreadCpuTime();
			require(time>=0, "Detailed fill clocks must run on the live scoring thread");
			return time;
		}
		private final ThreadMXBean clock;
		public long setupCpu, splitTrials;
		public final long[] states=new long[B+1], updatedCells=new long[B+1], clearedSlots=new long[B+1];
		public final long[] acquireCpu=new long[B+1], clearCpu=new long[B+1], recurrenceCpu=new long[B+1], releaseCpu=new long[B+1];
		public final long[] allocations=new long[B+1], allocatedCapacity=new long[B+1], reuses=new long[B+1], reusedCapacity=new long[B+1];
	}

	/** Uses the same root scan for endpoint scoring and diagnostic traceback.
	 * Trace fills have no propagated contact deck; their contacts come from the
	 * actual winning path below. No second fill or wider search is performed. */
	private EndpointResult selectEndpoints(final Fill fill, final CovarianceModelLengthBands.Layout root,
			final int maxWindow, final int candidateFrom, final int candidateTo){
		assert(fill!=null && root!=null && maxWindow>0 && candidateFrom<=candidateTo) :
			"Endpoint callers validate the window and candidate before selecting a finite root";
		float best=Float.NEGATIVE_INFINITY, unrestricted=Float.NEGATIVE_INFINITY;
		int from=0, to=0, unrestrictedFrom=0, unrestrictedTo=0, contacts=0;
		for(int j=Math.max(1, root.firstEnd); j<=root.lastEnd && root.width>0; j++){
			final int row=root.row(j), low=root.lower(j);
			final int last=Math.min(root.upper(j), maxWindow);
			for(int d=Math.max(1, low); d<=last; d=root.next(j, d)){
				final int index=root.union ? root.index(j, d) : row+d-low;
				if(index<0){continue;}
				final float value=fill.storage.active[0][index];
				if(!Float.isFinite(value)){continue;}
				if(value>unrestricted){
					unrestricted=value;
					unrestrictedFrom=j-d+1;
					unrestrictedTo=j;
				}
				final long overlap=Math.max(0L, (long)Math.min(j, candidateTo)-Math.max(j-d+1, candidateFrom)+1);
				if(value>best && 2*overlap>=Math.min((long)d, (long)candidateTo-candidateFrom+1)){
					best=value;
					from=j-d+1;
					to=j;
					contacts=fill.storage.activeEdges==null ? 0 : fill.storage.activeEdges[0][index]&255;
				}
			}
		}
		return new EndpointResult(best, from, to, contacts, fill.visited, fill.storage.peak,
			unrestricted, unrestrictedFrom, unrestrictedTo);
	}

	/** Traces the selected restricted optimum for contact diagnosis only. This
	 * retains ordinary choice bytes and bifurcation-child scores, so its CPU and
	 * memory are not endpoint-only performance measurements. Coordinates are
	 * consumed-prefix boundaries: each path cell covers[start,end). Callers can
	 * join these states to learned support, raw anchors and physical clipping.
	 * Impossible searches have empty path lists, even if some state layouts exist.
	 */
	public DiagnosticResult diagnoseOverlapping(final byte[] sequence, final int maxWindow,
			final int candidateFrom, final int candidateTo, final CovarianceModelPlacementBands placement,
			final CovarianceModelStateCorridors corridors){
		require(sequence!=null && maxWindow>0 && candidateFrom<=candidateTo && !configuration.local
			&& placement!=null && placement.model==model && placement.length==sequence.length
			&& corridors!=null && corridors.belongsTo(model),
			"Contact diagnosis requires the actual global CM, complete window, glocal anchors and learned resource");
		final CovarianceModelLengthBands.Layout[] layout=placement.layouts(corridors);
		final Fill fill=fill(sequence, layout, false);
		final EndpointResult best=selectEndpoints(fill, layout[0], maxWindow, candidateFrom, candidateTo);
		final IntList states=new IntList(), starts=new IntList(), ends=new IntList(), masks=new IntList();
		int contacts=0, empty=0, emptyObserved=0;
		for(int state=0; state<layout.length; state++){
			if(layout[state].cells==0){
				empty++;
				if(corridors.observed(state)){emptyObserved++;}
			}
		}
		if(Float.isFinite(best.score)){
			final IntList stack=new IntList();
			final int[] split=new int[1];
			push(stack, 0, best.targetTo, best.targetTo-best.targetFrom+1, -1, 0);
			while(stack.size>0){
				stack.pop(); // Side and parent are unnecessary for cell-only diagnostics.
				stack.pop();
				final int d=stack.pop(), j=stack.pop(), state=stack.pop();
				final int index=layout[state].index(j, d);
				require(index>=0, "A selected diagnostic path cannot leave its own allocated state layout: state="+state);
				final int mask=CovarianceModelPlacementBands.edgeMask(layout[state], j-d, j);
				states.add(state);
				starts.add(j-d);
				ends.add(j);
				masks.add(mask);
				contacts|=mask;
				final byte type=model.type[state];
				if(type==E){
					require(d==0, "A diagnostic END cell consumes zero residues, as in the scored recurrence");
				}else if(type==B){
					require(Float.isFinite(splitScore(state, j, d, layout, fill.storage.active, split)),
						"A diagnostic bifurcation must recover its finite scored split from pinned children");
					push(stack, model.childCount[state], j, split[0], -1, 0);
					push(stack, model.childFirst[state], j-split[0], d-split[0], -1, 0);
				}else{
					final int choice=fill.choices[state][index]&255;
					require(choice<model.childCount[state], "A finite diagnostic state must retain an ordinary child choice: state="+state);
					push(stack, model.childFirst[state]+choice, j-CovarianceModelCyk.rightDelta(type),
						d-CovarianceModelCyk.delta(type), -1, 0);
				}
			}
		}
		final EndpointResult endpoint=new EndpointResult(best.score, best.targetFrom, best.targetTo, contacts,
			best.cells, best.peakScoreCapacity, best.unrestrictedScore, best.unrestrictedFrom, best.unrestrictedTo);
		return new DiagnosticResult(endpoint, states, starts, ends, masks, empty, emptyObserved);
	}

	/** One recurrence serves full-input traceback and endpoint-only search.
	 * Endpoint decks recycle bifurcation children as soon as their last consumer
	 * finishes; only the traceback route pins them for later split recovery. */
	private Fill fill(byte[] sequence, CovarianceModelLengthBands.Layout[] layout, boolean endpointOnly){
		return fill(sequence, layout, endpointOnly, null);
	}
	private Fill fill(byte[] sequence, CovarianceModelLengthBands.Layout[] layout, boolean endpointOnly, final KernelProfile profile){
		long phase=profile==null ? 0 : profile.now();
		require(sequence!=null && layout!=null && layout.length==model.states(), "Candidate layouts must cover every state of the bound model");
		final byte[] bases=CovarianceModelAlphabet.encode(sequence);final int length=sequence.length, states=model.states();
		long visited=0, traceCells=0;
		for(int v=0; v<states; v++){
			require(layout[v]!=null && layout[v].length==length && (model.type[v]!=E || layout[v].max==0), "Layouts bind the input length and END states can only consume zero residues");visited+=layout[v].cells;
			if(!endpointOnly && model.type[v]!=E && model.type[v]!=B && !(configuration.local && v==0)){traceCells+=layout[v].cells;}}
		require(traceCells<=maxTraceCells, "Band traceback exceeds the declared byte-cell limit: "+traceCells+" > "+maxTraceCells);
		final byte[][] choices=endpointOnly ? null : new byte[states][];
		if(endpointOnly && endpointStorage==null){endpointStorage=new Decks(states, maxScoreCells, true);}
		final Decks storage=endpointOnly ? endpointStorage : new Decks(states, maxScoreCells, false);
		storage.reset();
		final int[] pending=parents.clone();
		final int[] split=endpointOnly ? new int[2] : null;
		float localBest=Float.NEGATIVE_INFINITY;int localBegin=-1;
		if(profile!=null){profile.setupCpu+=profile.now()-phase;}
		for(int v=states-1; v>=0; v--){
			final byte type=model.type[v];
			if(profile!=null){phase=profile.now();}
			final CovarianceModelLengthBands.Layout l=layout[v];
			final float[] s=storage.acquire(v, l.cells, profile, type);
			final byte[] edges=endpointOnly ? storage.activeEdges[v] : null;
			if(profile!=null){
				final long now=profile.now();
				profile.acquireCpu[type]+=now-phase;
				phase=now;
			}
			Arrays.fill(s, 0, l.cells, Float.NEGATIVE_INFINITY);
			if(profile!=null){
				final long now=profile.now();
				profile.clearCpu[type]+=now-phase;
				profile.clearedSlots[type]+=l.cells;
				phase=now;
			}
			if(configuration.local && v==0){
				final int root=l.index(length, length);if(root>=0){s[root]=localBest;}
			}else if(type==E){
				// Union rows have no retained cells in their gaps. Initialize only
				// real zero-length END addresses, never a rectangular surrogate.
				Arrays.fill(s, 0, l.cells, 0);
				if(endpointOnly){
					for(int j=l.firstEnd; j<=l.lastEnd && l.width>0; j++){
						final int index=l.index(j, 0);
						if(index>=0){edges[index]=(byte)CovarianceModelPlacementBands.edgeMask(l, j, j);}
					}
				}
			}
			else if(type==B){
				for(int j=l.firstEnd; j<=l.lastEnd && l.width>0; j++){final int row=l.row(j), low=l.lower(j);for(int d=low; d<=l.upper(j); d=l.next(j, d)){
					final int index=l.union ? l.index(j, d) : row+d-low;
					if(index<0){continue;}
					s[index]=splitScore(v, j, d, layout, storage.active, split, storage.activeEdges, profile);
					if(endpointOnly){edges[index]=(byte)(split[1]|CovarianceModelPlacementBands.edgeMask(l, j-d, j));}
				}}
			}else{
				final byte[] trace=endpointOnly ? null : new byte[l.cells];
				if(trace!=null){
					Arrays.fill(trace, (byte)255);
					choices[v]=trace;
				}
				final int delta=CovarianceModelCyk.delta(type), right=CovarianceModelCyk.rightDelta(type);
				final int first=model.childFirst[v], count=model.childCount[v];
				final float[] transitions=configuration.transitionScore[v];
				boolean directRows=true;
				for(int c=0; c<count; c++){
					if(layout[first+c].union){directRows=false;break;}
				}
				for(int j=Math.max(l.firstEnd, delta); j<=l.lastEnd && l.width>0; j++){
					final int row=l.row(j), low=l.lower(j), childEnd=j-right;
					if(directRows){prepareChildRows(first, count, childEnd, layout);}
					for(int d=Math.max(low, delta); d<=l.upper(j); d=l.next(j, d)){
						final int cell=l.union ? l.index(j, d) : row+d-low;
						if(cell<0){continue;}
						final int childSpan=d-delta;
						float best=configuration.local ? configuration.endScore[v]+configuration.elSelf*(d-delta) : Float.NEGATIVE_INFINITY;
						int choice=Float.isFinite(best) ? 254 : 255;
						int contact=0;
						for(int c=0; c<count; c++){
							final int child=first+c;
							final int index;
							if(directRows){
								if(childSpan<childLow[c] || childSpan>childHigh[c]){continue;}
								index=childOffset[c]+childSpan;
							}else{
								index=layout[child].index(childEnd, childSpan);
								if(index<0){continue;}
							}
							assert(storage.active[child]!=null):"A child's last forward consumer has not completed, state="+v+", child="+child;
							final float x=storage.active[child][index]+transitions[c];
							if(x>best){
								best=x;
								choice=c;
								if(endpointOnly){contact=storage.activeEdges[child][index]&255;}
							}
						}
						if(delta>0){best+=emission(v, bases, j-d+1, j);}
						s[cell]=best;
						if(trace!=null){trace[cell]=(byte)choice;}
						if(endpointOnly){edges[cell]=(byte)(contact|CovarianceModelPlacementBands.edgeMask(l, j-d, j));}
					}
				}
			}
			// Capture the full-input entry before recycling its score deck.
			// Descending states and strict > preserve Infernal's entry tie order.
			if(configuration.local && Float.isFinite(configuration.beginScore[v])){
				final int full=l.index(length, length);
				if(full>=0){final float value=s[full]+configuration.beginScore[v];if(value>localBest){localBest=value;localBegin=v;}}
			}
			if(profile!=null){
				final long now=profile.now();
				profile.recurrenceCpu[type]+=now-phase;
				phase=now;
			}
			if(type==B){
				consume(v, model.childFirst[v], pending, storage, !endpointOnly);
				consume(v, model.childCount[v], pending, storage, !endpointOnly);
			}else if(type!=E){
				for(int c=0; c<model.childCount[v]; c++){consume(v, model.childFirst[v]+c, pending, storage, !endpointOnly);}
			}
			if(v!=0 && pending[v]==0 && (endpointOnly || !pinned[v])){storage.release(v);}
			if(profile!=null){
				profile.releaseCpu[type]+=profile.now()-phase;
				profile.states[type]++;
			}
		}
		return new Fill(bases, choices, storage, visited, traceCells, localBegin);
	}

	/** Every cell of one parent row addresses the same child end. Cache only
	 * integer geometry, never score/contact arrays that Decks may release. Union
	 * children keep index(), since a row's missing middle is not one interval. */
	private void prepareChildRows(final int first, final int count, final int end,
			final CovarianceModelLengthBands.Layout[] layouts){
		assert(first>=0 && count<=childOffset.length && first+count<=layouts.length)
			: "Constructor-sized row scratch covers this ordinary state's exact child axis";
		for(int c=0; c<count; c++){
			final CovarianceModelLengthBands.Layout child=layouts[first+c];
			assert(!child.union) : "A union's row gap must use the exact index fallback";
			final int row=child.row(end);
			if(row<0){
				childLow[c]=1;
				childHigh[c]=0;
				childOffset[c]=0;
			}else{
				final int low=child.lower(end);
				childLow[c]=low;
				childHigh[c]=child.upper(end);
				childOffset[c]=row-low;
			}
		}
	}

	/** Increasing split lengths and strict > preserve the exact scorer's tie rule. */
	private float splitScore(int v, int j, int d, CovarianceModelLengthBands.Layout[] layouts, float[][] scores, int[] selected){
		return splitScore(v, j, d, layouts, scores, selected, null);
	}
	private float splitScore(int v, int j, int d, CovarianceModelLengthBands.Layout[] layouts,
			float[][] scores, int[] selected, byte[][] edges){
		return splitScore(v, j, d, layouts, scores, selected, edges, null);
	}
	private float splitScore(int v, int j, int d, CovarianceModelLengthBands.Layout[] layouts,
			float[][] scores, int[] selected, byte[][] edges, final KernelProfile profile){
		assert(model.type[v]==B):"Only bifurcations combine two child intervals";
		final int left=model.childFirst[v], right=model.childCount[v];
		final CovarianceModelLengthBands.Layout a=layouts[left], b=layouts[right];
		assert(scores[left]!=null && scores[right]!=null):"Both child decks must survive until their bifurcation is evaluated, state="+v;
		final int first=Math.max(Math.max(b.lower(j), d-a.max), j-a.lastEnd);
		final int last=Math.min(Math.min(Math.min(d, b.upper(j)), d-a.min), j-a.firstEnd);
		float best=Float.NEGATIVE_INFINITY;
		int choice=-1, contact=0;
		if(profile!=null){profile.splitTrials+=Math.max(0L, (long)last-Math.max(0, first)+1);}
		for(int k=Math.max(0, first); k<=last; k++){
			final int x=a.index(j-k, d-k), y=b.index(j, k);if(x<0 || y<0){continue;}
			final float value=scores[left][x]+scores[right][y];
			if(value>best){
				best=value;
				choice=k;
				if(edges!=null){contact=edges[left][x]|edges[right][y];}
			}
		}
		if(selected!=null){
			selected[0]=choice;
			if(edges!=null){selected[1]=contact;}
		}
		return best;
	}
	private void consume(int v, int child, int[] pending, Decks storage, boolean keepPinned){
		if(child==v){return;}
		require(pending[child]>0, "Each non-self child edge must be consumed exactly once: "+v+" -> "+child);
		if(--pending[child]==0 && (!keepPinned || !pinned[child])){storage.release(child);}
	}
	private void edge(int v, int child){
		require(child>=v && child<model.states(), "Descending evaluation requires forward children: "+v+" -> "+child);
		if(child==v){require(model.type[v]==IL || model.type[v]==IR, "Only consuming insertion states can self-loop");}
		else{parents[child]++;}
	}
	private void trace(CovarianceModelCyk.Result r, byte[][] choices, CovarianceModelLengthBands.Layout[] layouts, float[][] scores, byte[] bases, int localBegin){
		assert(Float.isFinite(r.score)):"Impossible roots cannot have a trace";
		final IntList stack=new IntList();push(stack, 0, bases.length-1, bases.length-1, -1, 0);
		final boolean[] matched=new boolean[model.clen+1];final int[] split=new int[1];
		while(stack.size>0){
			final int side=stack.pop(), parent=stack.pop(), d=stack.pop(), j=stack.pop(), v=stack.pop();
			final int i=j-d+1, n=r.states.size;final byte t=model.traceType(v);
			final int index=t==EL ? -1 : layouts[v].index(j, d);
			require(t==EL && configuration.local || index>=0, "Traceback must stay within allocated bands or the explicit local-end state, state="+v);
			int choice=-1;
			if(configuration.local && v==0){choice=CovarianceModelCyk.LOCAL_BEGIN;}
			else if(t==B){require(Float.isFinite(splitScore(v, j, d, layouts, scores, split)), "A selected B cell must retain a finite split");choice=split[0];}
			else if(t!=E && t!=EL){
				choice=choices[v][index]&255;
				if(configuration.local && choice==254){require(Float.isFinite(configuration.endScore[v]), "Local-end byte must name an eligible exit");choice=CovarianceModelCyk.LOCAL_END;}
				else{require(choice<model.childCount[v], "A finite ordinary trace requires a valid child offset");}
			}
			r.states.add(v);r.from.add(i);r.to.add(j);r.choice.add(choice);r.left.add(-1);r.right.add(-1);
			if(parent>=0){if(side==0){r.left.set(parent, n);}else{r.right.set(parent, n);}}
			if(t==EL){
				for(int pos=i; pos<=j; pos++){
					require(pos>0 && pos<=r.stateAtBase.length && r.stateAtBase[pos-1]<0, "EL emits only unassigned input residues");
					r.stateAtBase[pos-1]=v;r.modelPosition[pos-1]=-1;r.insertAnchor[pos-1]=-1;r.elResidues++;
				}continue;
			}
			if(choice==CovarianceModelCyk.LOCAL_BEGIN){
				require(localBegin>0 && localBegin<model.states() && Float.isFinite(configuration.beginScore[localBegin]), "Finite local root requires the retained full-input entry");
				r.localBeginState=localBegin;push(stack, localBegin, j, d, n, 0);continue;
			}
			if(t==E){require(d==0, "END states consume no residues");continue;}
			if(t==B){push(stack, model.childCount[v], j, choice, n, 1);push(stack, model.childFirst[v], j-choice, d-choice, n, 0);continue;}
			final int node=model.node[v];
			if(t==MP || t==ML || t==IL){assign(r, matched, i, v, t==IL ? 0 : model.consensusLeft[node], t==IL ? model.consensusLeft[node] : -1);}
			if(t==MP || t==MR || t==IR){assign(r, matched, j, v, t==IR ? 0 : model.consensusRight[node], t==IR ? model.consensusRight[node]-1 : -1);}
			push(stack, choice==CovarianceModelCyk.LOCAL_END ? model.states() : model.childFirst[v]+choice, j-CovarianceModelCyk.rightDelta(t), d-CovarianceModelCyk.delta(t), n, 0);
		}
		for(int v:r.stateAtBase){require(v>=0, "Full-input traceback must emit every input residue exactly once");}
		r.deleted=model.clen-r.matched;require(r.matched+r.inserted+r.elResidues==bases.length-1, "Match/insert/EL counts must conserve input length");
		final float[] values=new float[r.states.size];
		for(int n=values.length-1; n>=0; n--){
			final int v=r.states.get(n), left=r.left.get(n), right=r.right.get(n);final byte type=model.traceType(v);
			if(type==EL){values[n]=configuration.elSelf*(r.to.get(n)-r.from.get(n)+1);}
			else if(type==E){values[n]=0;}
			else if(type==B){require(left>n && right>n, "B trace must include both descendant subtrees");values[n]=values[left]+values[right];}
			else if(r.choice.get(n)==CovarianceModelCyk.LOCAL_BEGIN){require(v==0 && left>n, "Local root must retain its selected entry");values[n]=values[left]+configuration.beginScore[r.states.get(left)];}
			else{require(left>n, "An ordinary trace state must include its child");values[n]=values[left]+(r.choice.get(n)==CovarianceModelCyk.LOCAL_END ? configuration.endScore[v] : configuration.transitionScore[v][r.choice.get(n)]);
				if(CovarianceModelCyk.delta(type)>0){values[n]+=emission(v, bases, r.from.get(n), r.to.get(n));}}
		}
		require(values[0]==r.score, "The emitted trace must exactly rescore to the banded optimum");
	}
	private void assign(CovarianceModelCyk.Result r, boolean[] matched, int position, int state, int modelPosition, int anchor){
		require(position>0 && position<=r.stateAtBase.length && r.stateAtBase[position-1]<0, "Each input position must receive exactly one trace emission");
		r.stateAtBase[position-1]=state;r.modelPosition[position-1]=modelPosition;r.insertAnchor[position-1]=anchor;
		if(modelPosition==0){require(anchor>=0 && anchor<=model.clen, "Insertion anchor must be inside the model");r.inserted++;}
		else{require(modelPosition>0 && modelPosition<=model.clen && !matched[modelPosition], "Model match positions must be unique and bounded");matched[modelPosition]=true;r.matched++;
			if(r.modelFrom==0 || modelPosition<r.modelFrom){r.modelFrom=modelPosition;}r.modelTo=Math.max(r.modelTo, modelPosition);}
	}
	private float emission(int v, byte[] bases, int i, int j){
		assert(i>=1 && i<=j && j<bases.length):"Emitting states require a nonempty bounded interval";
		final byte type=model.type[v];return emissions[v][type==MP ? (bases[i]<<4)|bases[j] : (type==ML || type==IL ? bases[i] : bases[j])];
	}
	private static void push(IntList stack, int v, int j, int d, int parent, int side){
		assert(v>=0 && d>=0 && d<=j):"Trace stack intervals must be valid before coordinate assignment";
		stack.add(v);stack.add(j);stack.add(d);stack.add(parent);stack.add(side);
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}

	/** Selected restricted path and empty-layout counts; not a wider-search oracle. */
	public static final class DiagnosticResult{
		DiagnosticResult(final EndpointResult endpoint_, final IntList states_, final IntList starts_,
				final IntList ends_, final IntList masks_, final int empty_, final int emptyObserved_){
			endpoint=endpoint_;
			states=states_;
			starts=starts_;
			ends=ends_;
			contactMasks=masks_;
			emptyStates=empty_;
			emptyObservedStates=emptyObserved_;
			assert(states.size==starts.size && states.size==ends.size && states.size==contactMasks.size) :
				"Every diagnostic trace step retains state, consumed-prefix interval and contact mask together";
		}
		public final EndpointResult endpoint;
		public final IntList states, starts, ends, contactMasks;
		public final int emptyStates, emptyObservedStates;
	}

	/** Scalar endpoint result: no per-residue assignment or parsetree arrays. */
	public static final class EndpointResult {
		EndpointResult(float score_, int from, int to, int contacts_, long cells_, long capacity,
				float unrestricted, int unrestrictedFrom_, int unrestrictedTo_){
			score=score_;
			targetFrom=from;
			targetTo=to;
			contacts=contacts_;
			cells=cells_;
			peakScoreCapacity=capacity;
			unrestrictedScore=unrestricted;
			unrestrictedFrom=unrestrictedFrom_;
			unrestrictedTo=unrestrictedTo_;
			assert(Float.isFinite(score) ? from>0 && to>=from : from==0 && to==0):
				"Only finite endpoint scores can select a nonempty root interval";
		}
		public final float score, unrestrictedScore;
		public final int targetFrom, targetTo, contacts, unrestrictedFrom, unrestrictedTo;
		public final long cells, peakScoreCapacity;
	}
	private static final class Fill {
		Fill(byte[] bases_, byte[][] choices_, Decks storage_, long visited_, long traceCells_, int localBegin_){
			bases=bases_;
			choices=choices_;
			storage=storage_;
			visited=visited_;
			traceCells=traceCells_;
			localBegin=localBegin_;
		}
		final byte[] bases;
		final byte[][] choices;
		final Decks storage;
		final long visited, traceCells;
		final int localBegin;
	}

	/** Best-fit reuse of earlier band allocations. Drop unused buffers before
	 * growing; never count unreachable arrays as still resident or hide pinned ones. */
	private static final class Decks {
		Decks(int states, long limit_, boolean contacts){
			active=new float[states][];
			free=new float[states][];
			activeEdges=contacts ? new byte[states][] : null;
			freeEdges=contacts ? new byte[states][] : null;
			limit=limit_;
		}
		void reset(){
			for(int v=0; v<active.length; v++){
				if(active[v]!=null){release(v);}
			}
			peak=resident;
		}
		float[] acquire(int state, int cells, final KernelProfile profile, final int type){
			assert(cells>=0):"Band layouts validate array dimensions before allocation";
			int best=-1;for(int i=0; i<count; i++){if(free[i].length>=cells && (best<0 || free[i].length<free[best].length)){best=i;}}
			if(best>=0){
				active[state]=free[best];
				if(activeEdges!=null){activeEdges[state]=freeEdges[best];}
				if(profile!=null){
					profile.reuses[type]++;
					profile.reusedCapacity[type]+=active[state].length;
				}
				free[best]=free[--count];
				free[count]=null;
				if(freeEdges!=null){
					freeEdges[best]=freeEdges[count];
					freeEdges[count]=null;
				}
				return active[state];
			}
			for(int i=0; i<count; i++){
				resident-=free[i].length;
				free[i]=null;
				if(freeEdges!=null){freeEdges[i]=null;}
			}
			count=0;
			require(cells<=limit-resident, "Active/pinned score buffers exceed declared cell capacity: resident="+resident+", new="+cells+", limit="+limit);
			active[state]=new float[cells];
			if(profile!=null){
				profile.allocations[type]++;
				profile.allocatedCapacity[type]+=cells;
			}
			if(activeEdges!=null){activeEdges[state]=new byte[cells];}
			resident+=cells;
			peak=Math.max(peak, resident);
			return active[state];
		}
		void release(int state){
			require(active[state]!=null, "A score deck cannot be released twice: state="+state);
			free[count]=active[state];
			active[state]=null;
			if(activeEdges!=null){
				freeEdges[count]=activeEdges[state];
				activeEdges[state]=null;
			}
			count++;
		}
		final float[][] active, free;
		final byte[][] activeEdges, freeEdges;
		final long limit;
		int count;
		long resident, peak;
	}
	public static final class Result {
		Result(CovarianceModelCyk.Result alignment_, long traceCells_, long scoreCells_, double beta_){alignment=alignment_;traceCells=traceCells_;peakScoreCapacity=scoreCells_;tailProbability=beta_;}
		public final CovarianceModelCyk.Result alignment;
		public final long traceCells, peakScoreCapacity;
		public final double tailProbability;
	}
	private final CovarianceModel model;
	private final CovarianceModelConfiguration configuration;
	private final float[][] emissions;
	private final CovarianceModelLengthBands bands;
	private final int[] parents;private final boolean[] pinned;
	private final int[] childOffset, childLow, childHigh;
	private final long maxScoreCells, maxTraceCells;
	private Decks endpointStorage;
}
