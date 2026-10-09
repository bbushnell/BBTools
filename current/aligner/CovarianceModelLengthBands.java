package aligner;

/** Model-only subsequence-length bounds stored in INFERNAL1/a CMs.
 * cm_file.c stores dmin2,dmin1,dmax1,dmax2. cm_qdband.c documents each
 * excluded length tail as <=beta; local-begin/no-local-end construction permits
 * global use. These are model-distribution bounds, not a guarantee that an
 * arbitrary sequence's best parse survives. No query anchors are used.
 * Sources: Infernal1.1.5 src/cm_file.c and src/cm_qdband.c, BandCalculationEngine.
 * @author Brian Bushnell, Raiden
 */
public final class CovarianceModelLengthBands {

	/** qdb=1 selects the tighter stored band; qdb=2 selects the wider stored band. */
	public CovarianceModelLengthBands(CovarianceModel model, int qdb){
		require(model!=null && (qdb==1 || qdb==2), "Choose one of the model's two declared QDB sets");
		final double beta1=beta(model, "QDBBETA1"), beta2=beta(model, "QDBBETA2");
		this.model=model;
		require(beta2<=beta1, "QDB set2 must have the smaller tail bound and wider intervals");
		set=qdb;tailProbability=qdb==1 ? beta1 : beta2;
		min=new int[model.states()];max=new int[model.states()];
		for(int v=0; v<min.length; v++){
			final int[] b=model.qdb[v];
			require(b!=null && b.length==4 && b[0]>=0 && b[0]<=b[1] && b[1]<=b[2] && b[2]<=b[3],
				"CM length bands must be ordered dmin2<=dmin1<=dmax1<=dmax2, state="+v);
			min[v]=b[qdb==1 ? 1 : 0];max[v]=b[qdb==1 ? 2 : 3];
			if(model.type[v]==CovarianceModel.E){require(min[v]==0 && max[v]==0, "END length distribution is concentrated at zero, state="+v);}
		}
	}
	private static double beta(CovarianceModel model, String key){
		final String value=model.header(key);require(value!=null, "Cannot report a tail probability without "+key);
		final double beta=Double.parseDouble(value);require(Double.isFinite(beta) && beta>0 && beta<.5, "Invalid model tail probability "+key+"="+value);return beta;
	}
	public Layout layout(int state, int length){
		require(state>=0 && state<min.length, "Band layout requires a declared CM state");return new Layout(length, min[state], max[state]);
	}

	/** Packed band cells in ascending end j, then subsequence length d. No cells
	 * exist for excluded lengths or d>j. Addresses are independent of array type. */
	static final class Layout {
		Layout(int length_, int min_, int max_){
			require(length_>=0 && min_>=0 && max_>=min_, "Band dimensions require nonnegative ordered lengths");
			length=length_;min=min_;max=Math.min(length_, max_);
			final long w=min>max ? 0 : (long)max-min+1;
			final long count=w*(w+1)/2+(long)(length-max)*w;
			require(count<=Integer.MAX_VALUE-8, "A band deck exceeds Java array indexing: "+count);
			width=(int)w;cells=(int)count;
			firstEnd=min>max ? 1 : min;lastEnd=min>max ? 0 : length;startLow=0;startHigh=length;offsets=null;
		}
		/** Candidate coordinates use consumed-prefix boundaries: start=j-d, end=j.
		 * Only rows inside the end corridor are represented; within each row the
		 * start corridor and optional length bounds are intersected. */
		Layout(int length_, int min_, int max_, int startLow_, int startHigh_, int endLow, int endHigh){
			require(length_>=0 && min_>=0 && max_>=min_, "Candidate layout needs ordered nonnegative lengths");
			require(startLow_>=0 && startHigh_>=startLow_ && startHigh_<=length_ && endLow>=0 && endHigh>=endLow && endHigh<=length_,
				"Candidate boundary corridors must be ordered and inside the original input");
			length=length_;min=min_;max=Math.min(max_, length);startLow=startLow_;startHigh=startHigh_;
			firstEnd=Math.max(endLow, min);lastEnd=endHigh;
			final long rows=Math.max(0L, (long)lastEnd-firstEnd+1);
			require(rows<Integer.MAX_VALUE-8, "Candidate row offsets must fit a Java array before allocation");
			offsets=new int[(int)rows+1];long count=0;int widest=0;
			for(int n=0; n<rows; n++){
				final int j=firstEnd+n;final long w=Math.max(0L, (long)upper(j)-lower(j)+1);count+=w;
				require(count<=Integer.MAX_VALUE-8, "Candidate deck exceeds Java array indexing: "+count);offsets[n+1]=(int)count;
				widest=Math.max(widest, (int)w);
			}
			width=widest;cells=(int)count;
		}
		int lower(int end){return Math.max(min, end-startHigh);}
		int upper(int end){return Math.min(Math.min(end, max), end-startLow);}
		/** First cell in an end-position row; -1 for an empty row. */
		int row(int end){
			if(offsets!=null){return end<firstEnd || end>lastEnd || lower(end)>upper(end) ? -1 : offsets[end-firstEnd];}
			if(width==0 || end<min || end>length){return -1;}
			final int n=Math.min(end, max)-min;
			return (int)((long)n*(n+1)/2+(long)Math.max(0, end-max)*width);
		}
		/** Missing band cells return -1; callers must use negative infinity. */
		int index(int end, int size){
			if(end<0 || end>length || size<lower(end) || size>upper(end)){return -1;}final int row=row(end);return row<0 ? -1 : row+size-lower(end);
		}
		final int length, min, max, width, cells;
		final int firstEnd, lastEnd, startLow, startHigh;
		private final int[] offsets;
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	public final int set;
	/** Declared per-tail probability from the stored CM, not a sensitivity estimate. */
	public final double tailProbability;
	final int[] min, max;
	final CovarianceModel model;
}
