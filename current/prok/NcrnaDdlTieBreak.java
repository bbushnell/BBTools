package prok;

import java.util.Arrays;
import structures.ByteBuilder;

/** Experimental secondary ordering only: integer hits, bounded penalty, ordinal.
 * Complete eligible orders and both single-term counterfactual orders are retained
 * before the cap. All scratch belongs to one worker; GC moments are read-only.
 * @author Brian, Ganyu
 */
public final class NcrnaDdlTieBreak{
	public NcrnaDdlTieBreak(final NcrnaWindowGcProfile profile_, final byte[][] consensus,
			final double gcWeight_, final double postingWeight_){
		require(profile_!=null && consensus!=null && consensus.length==profile_.models()
			&& Double.isFinite(gcWeight_) && Double.isFinite(postingWeight_) && gcWeight_>=0
			&& postingWeight_>=0 && gcWeight_+postingWeight_<1,
			"Tie penalties need a bound model axis and nonnegative total weight below one raw hit");
		profile=profile_;
		gcWeight=gcWeight_;
		postingWeight=postingWeight_;
		gcWeightText=Double.toString(gcWeight);
		postingWeightText=Double.toString(postingWeight);
		final int n=consensus.length;
		means=new double[n];
		fallback=new boolean[n];
		for(int m=0; m<n; m++){
			fallback[m]=!profile.known(m);
			means[m]=fallback[m] ? fraction(consensus[m]) : profile.mean(m);
		}
		gc=new double[n];
		posting=new double[n];
		penalty=new double[n];
		base=new int[n];
		gcOrder=new int[n];
		postingOrder=new int[n];
		order=new int[n];
		oldRank=new int[n];
		gcRank=new int[n];
		postingRank=new int[n];
		newRank=new int[n];
	}

	public int[] build(final byte[] window, final float[] hits_, final double[] sums,
			final int[] eligible, final int count){
		require(window!=null && hits_!=null && sums!=null && eligible!=null && hits_.length==means.length
			&& sums.length==means.length && count>=0 && count<=means.length && count<=eligible.length,
			"Tie ordering consumes complete raw scores and one distinct native-eligible prefix");
		hits=hits_;
		size=count;
		queryGc=fraction(window);
		Arrays.fill(oldRank, 0);
		Arrays.fill(gcRank, 0);
		Arrays.fill(postingRank, 0);
		Arrays.fill(newRank, 0);
		for(int i=0; i<count; i++){
			final int m=eligible[i];
			require(m>=0 && m<means.length && oldRank[m]==0, "Native eligibility must be distinct before tie sorting");
			oldRank[m]=-1;
			require(Float.isFinite(hits[m]) && hits[m]>=0 && hits[m]==(int)hits[m]
				&& Double.isFinite(sums[m]) && sums[m]>=0 && sums[m]<=hits[m],
				"Primary scores remain integer hits; exact-list specificity sums are bounded by those hits");
			gc[m]=Double.isFinite(queryGc) && Double.isFinite(means[m]) ? Math.min(1, Math.abs(queryGc-means[m])) : 0;
			posting[m]=hits[m]>0 ? sums[m]/hits[m] : 0;
			penalty[m]=gcWeight*gc[m]+postingWeight*posting[m];
			base[i]=gcOrder[i]=postingOrder[i]=order[i]=m;
		}
		sort(base, count, null);
		sort(gcOrder, count, gc);
		sort(postingOrder, count, posting);
		sort(order, count, penalty);
		for(int i=0; i<count; i++){
			oldRank[base[i]]=i+1;
			gcRank[gcOrder[i]]=i+1;
			postingRank[postingOrder[i]]=i+1;
			newRank[order[i]]=i+1;
		}
		return order;
	}

	/** Heap sort on primitive ordinals: no allocation and no floating hit subtraction. */
	private void sort(final int[] ids, final int count, final double[] secondary){
		assert(count>=0 && count<=ids.length) : "Live prefix "+count+" must fit ordinal capacity "+ids.length+" before heap indexing";
		for(int i=count/2-1; i>=0; i--){sift(ids, i, count, secondary);}
		for(int end=count-1; end>0; end--){
			final int value=ids[end];
			ids[end]=ids[0];
			ids[0]=value;
			sift(ids, 0, end, secondary);
		}
	}
	private void sift(final int[] ids, int root, final int count, final double[] secondary){
		assert(root>=0 && root<count && count<=ids.length) : "Heap root "+root+" must be inside [0,"+count+") with capacity "+ids.length+" before child indexing";
		while(2*root+1<count){
			int child=2*root+1;
			if(child+1<count && compare(ids[child], ids[child+1], secondary)<0){child++;}
			if(compare(ids[root], ids[child], secondary)>=0){return;}
			final int value=ids[root];
			ids[root]=ids[child];
			ids[child]=value;
			root=child;
		}
	}
	private int compare(final int a, final int b, final double[] secondary){
		final int primary=Integer.compare((int)hits[b], (int)hits[a]);
		if(primary!=0){return primary;}
		final int tie=secondary==null ? 0 : Double.compare(secondary[a], secondary[b]);
		return tie==0 ? Integer.compare(a, b) : tie;
	}
	public static double fraction(final byte[] bases){
		require(bases!=null, "GC requires real sequence bytes; an unknown fraction is separate from zero");
		int known=0;
		int cg=0;
		for(final byte b:bases){
			switch(b|32){
				case 'c': case 'g': cg++; known++; break;
				case 'a': case 't': known++; break;
				default: break;
			}
		}
		return known==0 ? Double.NaN : cg/(double)known;
	}

	/** One row per eligible model; ranks include tied groups beyond the caller cap.
	 * A profile's support/fallback flags and every term remain explicit. */
	public void append(final ByteBuilder out, final int input, final long read, final long window,
			final int winner, final boolean accepted){
		for(int i=0; i<size; i++){
			final int m=base[i];
			out.append(input).tab().append(read).tab().append(window).tab().append(m).tab().append((int)hits[m]).tab();
			if(Double.isFinite(queryGc)){out.append(queryGc, 12);}
			out.tab();
			if(Double.isFinite(means[m])){out.append(means[m], 12);}
			out.tab().append(profile.support(m)).tab().append(fallback[m]).tab()
				.append(gc[m], 12).tab().append(posting[m], 12).tab().append(gcWeightText).tab().append(postingWeightText).tab().append(penalty[m], 12).tab()
				.append(oldRank[m]).tab().append(gcRank[m]).tab().append(postingRank[m]).tab().append(newRank[m]).tab()
				.append(accepted).tab().append(accepted && m==winner).tab().append(Double.doubleToLongBits(penalty[m])).tab()
				.append(Double.doubleToLongBits(gc[m])).tab().append(Double.doubleToLongBits(posting[m])).nl();
		}
	}
	public int rank(final int model){return newRank[model];}
	public int originalRank(final int model){return oldRank[model];}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}
	public static final String HEADER="input_ordinal\tread_ordinal\twindow_ordinal\tmodel_ordinal\traw_hits\twindow_gc\tmodel_mean_gc\tknown_windows\tconsensus_fallback\tgc_difference\tposting_mean\tgc_weight\tposting_weight\tweighted_penalty\told_rank\tgc_only_rank\tposting_only_rank\tnew_rank\taccepted_window\twinning_model\tpenalty_bits\tgc_difference_bits\tposting_mean_bits";
	private final NcrnaWindowGcProfile profile;
	private final double gcWeight, postingWeight;
	private final String gcWeightText, postingWeightText;
	private final double[] means, gc, posting, penalty;
	private final boolean[] fallback;
	private final int[] base, gcOrder, postingOrder, order, oldRank, gcRank, postingRank, newRank;
	private float[] hits;
	private int size;
	private double queryGc;
}
