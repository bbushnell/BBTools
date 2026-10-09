package aligner;

import java.util.Arrays;
import structures.IntList;
import static aligner.CovarianceModel.*;

/** Full-input CYK within explicit bands; optional configured local begins/ends.
 * Global model scoring remains the default. Ordinary traceback choices
 * use one unsigned byte; B splits are recomputed from pinned child decks so a
 * split larger than255 is never truncated. Scores outside bands are impossible.
 * Banding may change the optimum; compare to the exact unbanded scorer.
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
		for(int v=0; v<model.states(); v++){
			final byte type=model.type[v];require(type>=D && type<=B, "Parsed model states are ordinary; virtual EL is configured separately");
			if(type==B){edge(v, model.childFirst[v]);edge(v, model.childCount[v]);pinned[model.childFirst[v]]=pinned[model.childCount[v]]=true;}
			else if(type==E){require(model.childCount[v]==0, "END states cannot have child edges");}
			else{require(model.childCount[v]>0 && model.childCount[v]<=(local ? 254 : 255) && model.transitionScore[v].length==model.childCount[v],
				"Byte choices reserve255 for unset and254 for local ends in local mode, state="+v);
				for(int c=0; c<model.childCount[v]; c++){edge(v, model.childFirst[v]+c);}}
		}
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
		require(sequence!=null && layout!=null && layout.length==model.states(), "Candidate layouts must cover every state of the bound model");
		final byte[] bases=CovarianceModelAlphabet.encode(sequence);final int length=sequence.length, states=model.states();
		long visited=0, traceCells=0;
		for(int v=0; v<states; v++){
			require(layout[v]!=null && layout[v].length==length && (model.type[v]!=E || layout[v].max==0), "Layouts bind the input length and END states can only consume zero residues");visited+=layout[v].cells;
			if(model.type[v]!=E && model.type[v]!=B && !(configuration.local && v==0)){traceCells+=layout[v].cells;}}
		require(traceCells<=maxTraceCells, "Band traceback exceeds the declared byte-cell limit: "+traceCells+" > "+maxTraceCells);
		final byte[][] choices=new byte[states][];final Decks storage=new Decks(states, maxScoreCells);
		final int[] pending=parents.clone();
		float localBest=Float.NEGATIVE_INFINITY;int localBegin=-1;
		for(int v=states-1; v>=0; v--){
			final CovarianceModelLengthBands.Layout l=layout[v];final float[] s=storage.acquire(l.cells);storage.active[v]=s;
			Arrays.fill(s, 0, l.cells, Float.NEGATIVE_INFINITY);final byte type=model.type[v];
			if(configuration.local && v==0){
				final int root=l.index(length, length);if(root>=0){s[root]=localBest;}
			}else if(type==E){Arrays.fill(s, 0, l.cells, 0);}
			else if(type==B){
				for(int j=l.firstEnd; j<=l.lastEnd && l.width>0; j++){final int row=l.row(j), low=l.lower(j);for(int d=low; d<=l.upper(j); d++){
					s[row+d-low]=splitScore(v, j, d, layout, storage.active, null);
				}}
			}else{
				final byte[] trace=new byte[l.cells];Arrays.fill(trace, (byte)255);choices[v]=trace;
				final int delta=CovarianceModelCyk.delta(type), right=CovarianceModelCyk.rightDelta(type);
				for(int j=Math.max(l.firstEnd, delta); j<=l.lastEnd && l.width>0; j++){
					final int row=l.row(j), low=l.lower(j);for(int d=Math.max(low, delta); d<=l.upper(j); d++){
						float best=configuration.local ? configuration.endScore[v]+configuration.elSelf*(d-delta) : Float.NEGATIVE_INFINITY;
						int choice=Float.isFinite(best) ? 254 : 255;
						for(int c=0; c<model.childCount[v]; c++){
							final int child=model.childFirst[v]+c, index=layout[child].index(j-right, d-delta);
							if(index<0){continue;}
							assert(storage.active[child]!=null):"A child's last forward consumer has not completed, state="+v+", child="+child;
							final float x=storage.active[child][index]+configuration.transitionScore[v][c];
							if(x>best){best=x;choice=c;}
						}
						if(delta>0){best+=emission(v, bases, j-d+1, j);}
						final int index=row+d-low;s[index]=best;trace[index]=(byte)choice;
					}
				}
			}
			// Capture the full-input entry before recycling its score deck.
			// Descending states and strict > preserve Infernal's entry tie order.
			if(configuration.local && Float.isFinite(configuration.beginScore[v])){
				final int full=l.index(length, length);
				if(full>=0){final float value=s[full]+configuration.beginScore[v];if(value>localBest){localBest=value;localBegin=v;}}
			}
			if(type==B){consume(v, model.childFirst[v], pending, storage);consume(v, model.childCount[v], pending, storage);}
			else if(type!=E){for(int c=0; c<model.childCount[v]; c++){consume(v, model.childFirst[v]+c, pending, storage);}}
			if(v!=0 && pending[v]==0 && !pinned[v]){storage.release(v);}
		}
		final int root=layout[0].index(length, length);final float score=root<0 ? Float.NEGATIVE_INFINITY : storage.active[0][root];
		final CovarianceModelCyk.Result alignment=new CovarianceModelCyk.Result(score, visited, length, 1, length, new IntList(), new IntList());
		if(Float.isFinite(score)){trace(alignment, choices, layout, storage.active, bases, localBegin);}
		return new Result(alignment, traceCells, storage.peak, tailProbability);
	}

	/** Increasing split lengths and strict > preserve the exact scorer's tie rule. */
	private float splitScore(int v, int j, int d, CovarianceModelLengthBands.Layout[] layouts, float[][] scores, int[] selected){
		assert(model.type[v]==B):"Only bifurcations combine two child intervals";
		final int left=model.childFirst[v], right=model.childCount[v];
		final CovarianceModelLengthBands.Layout a=layouts[left], b=layouts[right];
		assert(scores[left]!=null && scores[right]!=null):"B children must remain pinned through traceback, state="+v;
		final int first=Math.max(Math.max(b.lower(j), d-a.max), j-a.lastEnd);
		final int last=Math.min(Math.min(Math.min(d, b.upper(j)), d-a.min), j-a.firstEnd);
		float best=Float.NEGATIVE_INFINITY;int choice=-1;
		for(int k=Math.max(0, first); k<=last; k++){
			final int x=a.index(j-k, d-k), y=b.index(j, k);if(x<0 || y<0){continue;}
			final float value=scores[left][x]+scores[right][y];if(value>best){best=value;choice=k;}
		}
		if(selected!=null){selected[0]=choice;}return best;
	}
	private void consume(int v, int child, int[] pending, Decks storage){
		if(child==v){return;}
		require(pending[child]>0, "Each non-self child edge must be consumed exactly once: "+v+" -> "+child);
		if(--pending[child]==0 && !pinned[child]){storage.release(child);}
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

	/** Best-fit reuse of earlier band allocations. Drop unused buffers before
	 * growing; never count unreachable arrays as still resident or hide pinned ones. */
	private static final class Decks {
		Decks(int states, long limit_){active=new float[states][];free=new float[states][];limit=limit_;}
		float[] acquire(int cells){
			assert(cells>=0):"Band layouts validate array dimensions before allocation";
			int best=-1;for(int i=0; i<count; i++){if(free[i].length>=cells && (best<0 || free[i].length<free[best].length)){best=i;}}
			if(best>=0){final float[] result=free[best];free[best]=free[--count];free[count]=null;return result;}
			for(int i=0; i<count; i++){resident-=free[i].length;free[i]=null;}count=0;
			require(cells<=limit-resident, "Active/pinned score buffers exceed declared cell capacity: resident="+resident+", new="+cells+", limit="+limit);
			final float[] result=new float[cells];resident+=cells;peak=Math.max(peak, resident);return result;
		}
		void release(int state){
			require(active[state]!=null, "A score deck cannot be released twice: state="+state);
			free[count++]=active[state];active[state]=null;
		}
		final float[][] active, free;final long limit;int count;long resident, peak;
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
	private final long maxScoreCells, maxTraceCells;
}
