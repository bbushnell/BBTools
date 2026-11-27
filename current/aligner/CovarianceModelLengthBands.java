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
	public static final class Layout {
		/** No retained cells. Keep an explicit empty layout rather than repairing
		 * an outside learned interval into a point on a physical window edge. */
		static Layout empty(final int length){return new Layout(length);}
		private Layout(final int length_){
			require(length_>=0, "An empty state still belongs to an actual nonnegative window length");
			length=length_;
			min=1;
			max=0;
			width=0;
			cells=0;
			firstEnd=1;
			lastEnd=0;
			startLow=0;
			startHigh=length;
			offsets=null;
			union=false;
			regionEndLow=firstEnd;
			regionEndHigh=lastEnd;
			otherStartLow=0;
			otherStartHigh=0;
			otherEndLow=0;
			otherEndHigh=0;
		}
		Layout(int length_, int min_, int max_){
			require(length_>=0 && min_>=0 && max_>=min_, "Band dimensions require nonnegative ordered lengths");
			length=length_;min=min_;max=Math.min(length_, max_);
			final long w=min>max ? 0 : (long)max-min+1;
			final long count=w*(w+1)/2+(long)(length-max)*w;
			require(count<=Integer.MAX_VALUE-8, "A band deck exceeds Java array indexing: "+count);
			width=(int)w;cells=(int)count;
			firstEnd=min>max ? 1 : min;lastEnd=min>max ? 0 : length;startLow=0;startHigh=length;offsets=null;
			union=false;
			regionEndLow=firstEnd;
			regionEndHigh=lastEnd;
			otherStartLow=0;
			otherStartHigh=0;
			otherEndLow=0;
			otherEndHigh=0;
		}
		/** Candidate coordinates use consumed-prefix boundaries: start=j-d, end=j.
		 * Only rows inside the end corridor are represented; within each row the
		 * start corridor and optional length bounds are intersected. */
		public Layout(int length_, int min_, int max_, int startLow_, int startHigh_, int endLow, int endHigh){
			require(length_>=0 && min_>=0 && max_>=min_, "Candidate layout needs ordered nonnegative lengths");
			require(startLow_>=0 && startHigh_>=startLow_ && startHigh_<=length_ && endLow>=0 && endHigh>=endLow && endHigh<=length_,
				"Candidate boundary corridors must be ordered and inside the original input");
			length=length_;min=min_;max=Math.min(max_, length);startLow=startLow_;startHigh=startHigh_;
			firstEnd=Math.max(endLow, min);lastEnd=endHigh;
			union=false;
			regionEndLow=endLow;
			regionEndHigh=endHigh;
			otherStartLow=0;
			otherStartHigh=0;
			otherEndLow=0;
			otherEndHigh=0;
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
		/** Exact union of two endpoint rectangles, intersected with the same
		 * triangular/length domain. Each row packs one or two intervals, with no
		 * allocation for cross-corner cells. One primitive array holds row offsets,
		 * low/high lengths and the optional gap; there are no per-cell objects.
		 * The supplied scratch bounds are copied into scalar fields, never retained.
		 */
		public static Layout union(final int length, final int min, final int max, final int[] a, final int[] b){
			require(a!=null && a.length>=4 && b!=null && b.length>=4,
				"An exact placement union needs both complete endpoint rectangles");
			return new Layout(length, min, max, a, b);
		}
		private Layout(final int length_, final int min_, final int max_, final int[] a, final int[] b){
			require(length_>=0 && min_>=0 && max_>=min_, "Union dimensions use the scorer's nonnegative length domain");
			for(int side=0; side<2; side++){
				final int[] bounds=side==0 ? a : b;
				require(bounds[0]>=0 && bounds[0]<=bounds[1] && bounds[1]<=length_
					&& bounds[2]>=0 && bounds[2]<=bounds[3] && bounds[3]<=length_,
					"Each union component must be an actual ordered rectangle inside the search window");
			}
			length=length_;
			min=min_;
			max=Math.min(length, max_);
			startLow=a[0];
			startHigh=a[1];
			regionEndLow=a[2];
			regionEndHigh=a[3];
			otherStartLow=b[0];
			otherStartHigh=b[1];
			otherEndLow=b[2];
			otherEndHigh=b[3];
			firstEnd=Math.max(min, Math.min(a[2], b[2]));
			lastEnd=Math.max(a[3], b[3]);
			union=true;
			final long rows=Math.max(0L, (long)lastEnd-firstEnd+1);
			require(rows*5+1<=Integer.MAX_VALUE-8, "Union row metadata must fit one primitive array before allocation");
			offsets=new int[(int)(rows*5+1)];
			long count=0;
			int widest=0;
			for(int n=0; n<rows; n++){
				final int j=firstEnd+n, at=n*5;
				final int lowA=j>=a[2] && j<=a[3] ? Math.max(min, j-a[1]) : 1;
				final int highA=j>=a[2] && j<=a[3] ? Math.min(Math.min(j, max), j-a[0]) : 0;
				final int lowB=j>=b[2] && j<=b[3] ? Math.max(min, j-b[1]) : 1;
				final int highB=j>=b[2] && j<=b[3] ? Math.min(Math.min(j, max), j-b[0]) : 0;
				final boolean hasA=lowA<=highA, hasB=lowB<=highB;
				final int low=hasA ? hasB ? Math.min(lowA, lowB) : lowA : hasB ? lowB : 1;
				final int high=hasA ? hasB ? Math.max(highA, highB) : highA : hasB ? highB : 0;
				final int gap=hasA && hasB ? (int)Math.max(0L, (long)Math.max(lowA, lowB)-Math.min(highA, highB)-1) : 0;
				final int gapLow=gap>0 ? Math.min(highA, highB)+1 : 0;
				final long kept=Math.max(0L, (long)high-low+1)-gap;
				offsets[at]=(int)count;
				offsets[at+1]=low;
				offsets[at+2]=high;
				offsets[at+3]=gapLow;
				offsets[at+4]=gap;
				count+=kept;
				require(count<=Integer.MAX_VALUE-8, "An exact union deck exceeds Java array indexing: "+count);
				widest=Math.max(widest, (int)kept);
			}
			offsets[(int)rows*5]=(int)count;
			width=widest;
			cells=(int)count;
		}
		public int lower(int end){
			return union ? end<firstEnd || end>lastEnd ? 1 : offsets[(end-firstEnd)*5+1] : Math.max(min, end-startHigh);
		}
		public int upper(int end){
			return union ? end<firstEnd || end>lastEnd ? 0 : offsets[(end-firstEnd)*5+2] : Math.min(Math.min(end, max), end-startLow);
		}
		/** Next retained length after this one; jump over the union's row gap. */
		public int next(final int end, final int size){
			assert(end>=firstEnd && end<=lastEnd && size<Integer.MAX_VALUE) : "Row iteration stays inside its allocated end range";
			final int next=size+1;
			if(!union){return next;}
			final int at=(end-firstEnd)*5, gapLow=offsets[at+3], gap=offsets[at+4];
			return next>=gapLow && next<gapLow+gap ? gapLow+gap : next;
		}
		/** First cell in an end-position row; -1 for an empty row. */
		public int row(int end){
			if(union){return end<firstEnd || end>lastEnd || lower(end)>upper(end) ? -1 : offsets[(end-firstEnd)*5];}
			if(offsets!=null){return end<firstEnd || end>lastEnd || lower(end)>upper(end) ? -1 : offsets[end-firstEnd];}
			if(width==0 || end<min || end>length){return -1;}
			final int n=Math.min(end, max)-min;
			return (int)((long)n*(n+1)/2+(long)Math.max(0, end-max)*width);
		}
		/** Missing band cells return -1; callers must use negative infinity. */
		public int index(int end, int size){
			if(union){
				if(end<firstEnd || end>lastEnd){return -1;}
				final int at=(end-firstEnd)*5, low=offsets[at+1], high=offsets[at+2];
				if(size<low || size>high){return -1;}
				final int gapLow=offsets[at+3], gap=offsets[at+4];
				if(size>=gapLow && size<gapLow+gap){return -1;}
				return offsets[at]+size-low-(size>=gapLow+gap ? gap : 0);
			}
			if(end<0 || end>length || size<lower(end) || size>upper(end)){return -1;}final int row=row(end);return row<0 ? -1 : row+size-lower(end);
		}
		/** Exposed imposed endpoint faces, excluding the natural window edges.
		 * Test the rectangle union before intersecting the triangular/END domain:
		 * those model constraints must not invent new endpoint contacts. A missing
		 * diagonal detects an exposed concave corner even if both axial neighbors
		 * are covered. Covered internal component faces are not union boundaries.
		 */
		public int contactMask(final int start, final int end){
			require(start>=0 && start<=end && index(end, end-start)>=0,
				"Contacts belong to an actual retained consumed-prefix cell");
			if(!union){return faces(start, end, startLow, startHigh, firstEnd, lastEnd);}
			int candidates=0;
			if(rectangle(start, end, startLow, startHigh, regionEndLow, regionEndHigh)){
				candidates|=faces(start, end, startLow, startHigh, regionEndLow, regionEndHigh);
			}
			if(rectangle(start, end, otherStartLow, otherStartHigh, otherEndLow, otherEndHigh)){
				candidates|=faces(start, end, otherStartLow, otherStartHigh, otherEndLow, otherEndHigh);
			}
			if(candidates==0){return 0;}
			int result=0;
			for(int ds=-1; ds<=1; ds++){
				for(int de=-1; de<=1; de++){
					final int directions=(ds<0 ? 1 : ds>0 ? 2 : 0)|(de<0 ? 4 : de>0 ? 8 : 0);
					if((directions&candidates)==0){continue;}
					final long s=(long)start+ds, e=(long)end+de;
					if(s<0 || s>length || e<0 || e>length){continue;}
					if(!rectangle((int)s, (int)e, startLow, startHigh, regionEndLow, regionEndHigh)
						&& !rectangle((int)s, (int)e, otherStartLow, otherStartHigh, otherEndLow, otherEndHigh)){
						result|=directions&candidates;
					}
				}
			}
			return result;
		}
		private int faces(final int start, final int end, final int sl, final int sh, final int el, final int eh){
			assert(sl>=0 && sh<=length && el>=0 && eh<=length) : "Placement faces were clipped to this exact window";
			return (sl>0 && start==sl ? 1 : 0)|(sh<length && start==sh ? 2 : 0)
				|(el>0 && end==el ? 4 : 0)|(eh<length && end==eh ? 8 : 0);
		}
		private static boolean rectangle(final int start, final int end, final int sl, final int sh, final int el, final int eh){
			return start>=sl && start<=sh && end>=el && end<=eh;
		}
		public final int length, min, max, width, cells;
		final int firstEnd, lastEnd, startLow, startHigh;
		private final int[] offsets;
		final boolean union;
		private final int regionEndLow, regionEndHigh, otherStartLow, otherStartHigh, otherEndLow, otherEndHigh;
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	public final int set;
	/** Declared per-tail probability from the stored CM, not a sensitivity estimate. */
	public final double tailProbability;
	final int[] min, max;
	final CovarianceModel model;
}
