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
		require(model_!=null && length_>=0 && radius_>=0 && low_!=null && high_!=null
			&& low_.length==model_.clen+1 && high_.length==low_.length,
			"Composed anchors require one bounded interval for every CM boundary");
		model=model_;length=length_;radius=radius_;low=low_.clone();high=high_.clone();
		for(int i=0; i<low.length; i++){
			require(low[i]>=0 && low[i]<=high[i] && high[i]<=length
				&& (i==0 || low[i]>=low[i-1] && high[i]>=high[i-1]),
				"CM boundary intervals must be ordered within the complete candidate, boundary="+i);
		}
		require(low[0]==0 && high[model.clen]==length,
			"Global CM placement must retain both ends of the complete refined candidate");
	}
	public CovarianceModelPlacementBands(CovarianceModel model_, byte[] sequence, byte[] aligned, byte[] rf, int radius_){
		require(model_!=null && sequence!=null && aligned!=null && rf!=null && aligned.length==rf.length,
			"Placement requires one original sequence and its same-width alignment/RF columns");
		require(radius_>=0, "Anchor radius is an explicit nonnegative experiment parameter");
		model=model_;length=sequence.length;radius=radius_;
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
		assert(low.length==model.clen+1):"Each model boundary needs a sequence placement";
		final CovarianceModelLengthBands.Layout[] out=new CovarianceModelLengthBands.Layout[model.states()];
		for(int v=0; v<out.length; v++){
			final int n=model.node[v];final byte nt=model.nodeType[n], st=model.type[v];
			final boolean insert=st==IL || st==IR;
			// CreateEmitMap/consensusMap stores emitted positions for MAT sides and
			// adjacent consensus boundaries for nonemitting sides. Insert states
			// follow the node's main state and therefore lie inside its emissions.
			final int a=model.consensusLeft[n]-(!insert && (nt==MATP || nt==MATL) ? 1 : 0);
			final int b=model.consensusRight[n]-(!insert && (nt==MATP || nt==MATR) ? 0 : 1);
			require(a>=0 && a<=b && b<=model.clen, "Node emission map must define an ordered consensus subtree, state="+v+", boundaries="+a+","+b);
			final int startLow=(int)Math.max(0L, (long)low[a]-radius), startHigh=(int)Math.min(length, (long)high[a]+radius);
			final int endLow=(int)Math.max(0L, (long)low[b]-radius), endHigh=(int)Math.min(length, (long)high[b]+radius);
			out[v]=new CovarianceModelLengthBands.Layout(length, 0, st==E ? 0 : length, startLow, startHigh, endLow, endHigh);
		}
		return out;
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
		require(layout!=null && start>=0 && start<=end && layout.index(end, end-start)>=0,
			"Contact coordinates must be consumed-prefix boundaries of a retained trace cell");
		return (layout.startLow>0 && start==layout.startLow ? 1 : 0)
			| (layout.startHigh<layout.length && start==layout.startHigh ? 2 : 0)
			| (layout.firstEnd>0 && end==layout.firstEnd ? 4 : 0)
			| (layout.lastEnd<layout.length && end==layout.lastEnd ? 8 : 0);
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
}
