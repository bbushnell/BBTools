package aligner;

import static aligner.CovarianceModel.*;

/** Candidate-specific endpoint corridors derived from a supplied global
 * sequence-to-consensus alignment. Oracle anchors are a best-case experiment,
 * not an independently available production filter. No QDB intersection is
 * imposed: the short-input QDB counterexample must remain representable.
 * @author Brian Bushnell, Raiden
 */
public final class CovarianceModelPlacementBands {
	/** Already composed CM-boundary intervals on the exact candidate coordinate
	 * axis. Copies protect a reusable scorer from later caller-buffer mutation. */
	public CovarianceModelPlacementBands(CovarianceModel model_, int length_, int[] low_, int[] high_, int radius_){
		this(model_, length_, low_, high_, radius_, false);
	}
	/** Free-terminal placement on a padded window; terminal anchors describe the
	 * predicted CM extent rather than forcing the whole window into the model. */
	public static CovarianceModelPlacementBands window(CovarianceModel model, int length,
			int[] low, int[] high, int radius){
		return new CovarianceModelPlacementBands(model, length, low, high, radius, true);
	}
	private CovarianceModelPlacementBands(CovarianceModel model_, int length_, int[] low_, int[] high_, int radius_, boolean freeEnds){
		this(model_, length_, low_, high_, radius_, freeEnds, false);
	}
	/** Only the explicit separate-context mapper supplies outside raw anchors. */
	static CovarianceModelPlacementBands windowUnclipped(final CovarianceModel model, final int length,
			final int[] low, final int[] high, final int radius){
		return new CovarianceModelPlacementBands(model, length, low, high, radius, true, true);
	}
	private CovarianceModelPlacementBands(final CovarianceModel model_, final int length_, final int[] low_,
			final int[] high_, final int radius_, final boolean freeEnds, final boolean unclipped_){
		require(model_!=null && length_>=0 && radius_>=0 && low_!=null && high_!=null
			&& low_.length==model_.clen+1 && high_.length==low_.length,
			"Composed anchors require one bounded interval for every CM boundary");
		model=model_;
		length=length_;
		radius=radius_;
		low=low_.clone();
		high=high_.clone();
		unclipped=unclipped_;
		require(!unclipped || freeEnds, "Outside context anchors are meaningful only for a free-terminal search window");
		for(int i=0; i<low.length; i++){
			require(low[i]<=high[i] && (unclipped || low[i]>=0 && high[i]<=length)
				&& (i==0 || low[i]>=low[i-1] && high[i]>=high[i-1]),
				"CM boundary intervals must stay ordered; outside coordinates require an explicit unclipped context, boundary="+i);
		}
		require(freeEnds || low[0]==0 && high[model.clen]==length,
			"Global CM placement must retain both ends of the complete refined candidate");
	}
	public CovarianceModelPlacementBands(CovarianceModel model_, byte[] sequence, byte[] aligned, byte[] rf, int radius_){
		require(model_!=null && sequence!=null && aligned!=null && rf!=null && aligned.length==rf.length,
			"Placement requires one original sequence and its same-width alignment/RF columns");
		require(radius_>=0, "Anchor radius is an explicit nonnegative experiment parameter");
		model=model_;length=sequence.length;radius=radius_;
		unclipped=false;
		low=new int[model.clen+1];high=new int[model.clen+1];int m=0, q=0;
		for(int col=0; col<rf.length; col++){
			final boolean match=letter(rf[col]);require(match || gap(rf[col]), "RF must describe a global consensus without local ends");
			if(match){require(++m<=model.clen, "RF contains more consensus columns than the bound model");}
			final byte b=aligned[col];
			if(CovarianceModelAlphabet.code(b)>=0){
				require(q<length && CovarianceModelAlphabet.code(b)==CovarianceModelAlphabet.code(sequence[q]), "Anchor alignment must reconstruct every original residue in order, position="+(q+1));q++;
			}else{require(gap(b), "Anchor alignment contains an unsupported residue or local end");}
			if(match){low[m]=q;}high[m]=q;
		}
		require(m==model.clen && q==length, "Global anchors must conserve both consensus length and the complete original input");
	}

	CovarianceModelLengthBands.Layout[] layouts(){
		return layouts(null);
	}
	CovarianceModelLengthBands.Layout[] layouts(final CovarianceModelStateCorridors corridors){
		return layouts(corridors, -1);
	}
	/** The optional second component is a true cell union, not an endpoint hull.
	 * Unseen states retain their existing uniform radius; supported learned states
	 * may additionally use the declared narrow glocal region in the same fill. */
	CovarianceModelLengthBands.Layout[] layouts(final CovarianceModelStateCorridors corridors, final int unionRadius){
		assert(low.length==model.clen+1):"Each model boundary needs a sequence placement";
		require(corridors==null || corridors.belongsTo(model), "Learned state bounds must bind this placement's exact CM");
		require(unionRadius==-1 || corridors!=null && unionRadius>=0 && unionRadius<=radius,
			"A narrow union component requires learned bounds and cannot exceed the unseen-state fallback radius");
		final CovarianceModelLengthBands.Layout[] out=new CovarianceModelLengthBands.Layout[model.states()];
		final int[] bounds=new int[4];
		final int[] narrow=unionRadius<0 ? null : new int[4];
		for(int v=0; v<out.length; v++){
			stateBounds(v, bounds);
			final byte st=model.type[v];
			final boolean joined=unionRadius>=0 && corridors.observed(v);
			boolean narrowPresent=joined;
			if(joined){
				narrow[0]=(int)Math.max(0L, (long)bounds[0]-unionRadius);
				narrow[1]=(int)Math.min(length, (long)bounds[1]+unionRadius);
				narrow[2]=(int)Math.max(0L, (long)bounds[2]-unionRadius);
				narrow[3]=(int)Math.min(length, (long)bounds[3]+unionRadius);
				narrowPresent=narrow[0]<=narrow[1] && narrow[2]<=narrow[3];
			}
			if(corridors!=null){
				final boolean present=unclipped ? corridors.projectUnclipped(v, length, bounds, radius)
					: corridors.project(v, length, bounds, radius);
				if(!present){
					out[v]=narrowPresent ? new CovarianceModelLengthBands.Layout(length, 0, st==E ? 0 : length,
						narrow[0], narrow[1], narrow[2], narrow[3]) : CovarianceModelLengthBands.Layout.empty(length);
					continue;
				}
			}else{
				bounds[0]=(int)Math.max(0L, (long)bounds[0]-radius);
				bounds[1]=(int)Math.min(length, (long)bounds[1]+radius);
				bounds[2]=(int)Math.max(0L, (long)bounds[2]-radius);
				bounds[3]=(int)Math.min(length, (long)bounds[3]+radius);
				if(bounds[0]>bounds[1] || bounds[2]>bounds[3]){
					out[v]=CovarianceModelLengthBands.Layout.empty(length);
					continue;
				}
			}
			out[v]=narrowPresent ? CovarianceModelLengthBands.Layout.union(length, 0, st==E ? 0 : length, bounds, narrow)
				: new CovarianceModelLengthBands.Layout(length, 0, st==E ? 0 : length, bounds[0], bounds[1], bounds[2], bounds[3]);
		}
		return out;
	}
	/** Unexpanded glocal-projected consumed-prefix boundaries for one CM state:
	 * start low/high, end low/high. A caller-owned four-int array is reused across
	 * trace steps; the radius is deliberately absent from these training anchors.
	 * An explicitly unclipped context placement can return outside coordinates;
	 * its layouts apply residuals before intersecting the search window.
	 */
	public void stateBounds(final int state, final int[] bounds){
		require(state>=0 && state<model.states() && bounds!=null && bounds.length>=4,
			"State-anchor export requires a real CM state and reusable four-value output");
		final int n=model.node[state];
		final byte nt=model.nodeType[n], st=model.type[state];
		final boolean insert=st==IL || st==IR;
		// CreateEmitMap stores emitted MAT positions and adjacent nonemitting
		// boundaries. Insert states lie inside the node's own emissions.
		final int a=model.consensusLeft[n]-(!insert && (nt==MATP || nt==MATL) ? 1 : 0);
		final int b=model.consensusRight[n]-(!insert && (nt==MATP || nt==MATR) ? 0 : 1);
		if(a<0 || a>b || b>model.clen){
			throw new IllegalArgumentException("CM node emission map must define an ordered subtree before projecting a state: state="+state+", boundaries="+a+","+b);
		}
		bounds[0]=low[a];
		bounds[1]=high[a];
		bounds[2]=low[b];
		bounds[3]=high[b];
	}
	/** Contacts with imposed start/end corridors, excluding natural input edges.
	 * A contact is a diagnostic, not proof that the exact optimum was excluded.
	 * Do not count the triangular d=0/d=j boundary or END's required zero length.
	 * Multiple edge directions at one trace step count once in steps. */
	static EdgeContacts edgeContacts(CovarianceModelCyk.Result trace, CovarianceModelLengthBands.Layout[] layout){
		return edgeContacts(trace, layout, false);
	}
	static EdgeContacts edgeContacts(CovarianceModelCyk.Result trace, CovarianceModelLengthBands.Layout[] layout, boolean local){
		require(trace!=null && layout!=null && trace.from.size==trace.states.size && trace.to.size==trace.states.size,
			"Band contacts require the actual traceback and its state layouts");
		require(!local || trace.choice.size==trace.states.size && trace.left.size==trace.states.size && trace.right.size==trace.states.size,
			"Local contact diagnostics require complete state/choice/child metadata to identify virtual EL");
		int steps=0, sl=0, sh=0, el=0, eh=0;
		for(int i=0; i<trace.states.size; i++){
			final int v=trace.states.get(i);
			if(local && v==layout.length){
				require(layout.length>0 && trace.choice.get(i)==-1 && trace.left.get(i)==-1 && trace.right.get(i)==-1
					&& trace.from.get(i)>=1 && trace.to.get(i)>=trace.from.get(i)-1 && trace.to.get(i)<=layout[0].length,
					"Virtual EL is a bounded terminal interval, including an empty interval");
				// EL has no allocated band. Its ordinary parent is still checked;
				// counting EL against an invented layout would fabricate contact.
				continue;
			}
			require(v>=0 && v<layout.length, "Trace state must index its own placement layout");
			final int mask=edgeMask(layout[v], trace.from.get(i)-1, trace.to.get(i));
			if(mask!=0){steps++;}if((mask&1)!=0){sl++;}if((mask&2)!=0){sh++;}if((mask&4)!=0){el++;}if((mask&8)!=0){eh++;}
		}
		return new EdgeContacts(steps, sl, sh, el, eh);
	}
	static int edgeMask(CovarianceModelLengthBands.Layout layout, int start, int end){
		require(layout!=null, "A contact mask needs the actual scored state layout");
		return layout.contactMask(start, end);
	}
	public static final class EdgeContacts {
		EdgeContacts(int steps_, int sl, int sh, int el, int eh){steps=steps_;startLow=sl;startHigh=sh;endLow=el;endHigh=eh;}
		public final int steps, startLow, startHigh, endLow, endHigh;
	}
	static final EdgeContacts NO_CONTACTS=new EdgeContacts(0, 0, 0, 0, 0);
	private static boolean letter(byte b){return b>='A' && b<='Z' || b>='a' && b<='z';}
	private static boolean gap(byte b){return b=='.' || b=='-';}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	final CovarianceModel model;
	public final int length, radius;
	final int[] low, high;
	private final boolean unclipped;
}
