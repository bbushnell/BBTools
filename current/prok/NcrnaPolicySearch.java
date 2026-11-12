package prok;

import java.util.ArrayList;

import map.ObjectIntMap;
import parse.LineParser1;
import prok.NcrnaJointPolicy.Member;
import prok.NcrnaTextIO.Input;

/** Finite coordinate search over complete adaptive ncRNA mapping policies.
 * Discovery and grading are supplied by the native caller, never approximated
 * here. A move improves measured recall below the requested target, or reduces
 * logical alignment cost while meeting that target. Dominating moves take
 * precedence; equal outcomes stop rather than wandering across a plateau.
 * This finds a local operating point, not a global optimum.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaPolicySearch{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaPolicySearch(final NcrnaJointPolicy start_, final Axis[] axes_){
		require(start_!=null && axes_!=null && axes_.length==KNOBS.length*start_.size(),
			"The search grid must declare five controls for every member");
		start=start_;
		axes=axes_.clone();
		initial=new int[axes.length];
		final boolean[] seen=new boolean[axes.length];
		for(int i=0; i<axes.length; i++){
			final Axis axis=axes[i];
			require(axis!=null && axis.member>=0 && axis.member<start.size(), "Grid member is outside the policy roster");
			final int slot=axis.member*KNOBS.length+axis.knob;
			require(!seen[slot], "A member/control pair appears twice in the grid");
			seen[slot]=true;
			final Member member=start.member(axis.member);
			require(member.adaptive, "Search starts must explicitly use adaptive mode; fixed/adaptive equivalence is a separate caller test");
			final double value=value(member, axis.knob);
			initial[i]=axis.index(value);
			require(initial[i]>=0, "Starting value is absent from the grid for "+member.name+" "+KNOBS[axis.knob]);
		}
	}

	/** Finite ordered choices; values are canonicalized to the actual caller types. */
	public static Axis[] loadGrid(final String path, final NcrnaJointPolicy start){
		final ObjectIntMap<String> members=new ObjectIntMap<String>(String.class);
		for(int i=0; i<start.size(); i++){members.put(start.member(i).name, i);}
		final ArrayList<Axis> axes=new ArrayList<Axis>();
		final LineParser1 parser=new LineParser1('\t'), numbers=new LineParser1(',');
		try(Input input=new Input(path)){
			final byte[] header=input.reader.nextLine();
			require(header!=null, "Missing grid header");
			parser.set(header);
			require(parser.terms()==3 && parser.termEquals("member", 0) && parser.termEquals("knob", 1)
				&& parser.termEquals("values", 2), "Grid header must be member, knob, values separated by tabs");
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				require(parser.terms()==3, "Each grid row needs member, knob and comma-separated values");
				final int member=members.get(parser.parseString(0));
				require(member>=0, "Grid names an unknown policy member");
				int knob=-1;
				for(int i=0; i<KNOBS.length; i++){if(parser.termEquals(KNOBS[i], 1)){knob=i;}}
				require(knob>=0, "Grid knob must be topn, floor, topfrac, qfrac or patience");
				// One retained small array per axis; no genome or sequence parsing here.
				numbers.set(parser.parseByteArray(2));
				final double[] values=new double[numbers.terms()];
				// Small startup grid only; preserve exponent/range semantics without
				// changing the native parser's global setting in parent or child JVMs.
				for(int i=0; i<values.length; i++){values[i]=Double.parseDouble(numbers.parseString(i));}
				axes.add(new Axis(member, knob, values));
			}
		}
		return axes.toArray(new Axis[axes.size()]);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Neighbors           ----------------*/
	/*--------------------------------------------------------------*/

	public int[] initial(){return initial.clone();}
	public int dimensions(){return axes.length;}

	/** Startup serialization cases, not scientific evaluations. Each axis value
	 * appears in a complete policy while other coordinates retain their start.
	 */
	public ArrayList<NcrnaJointPolicy> serializationCases(){
		final ArrayList<NcrnaJointPolicy> cases=new ArrayList<NcrnaJointPolicy>();
		for(int a=0; a<axes.length; a++){
			for(int v=0; v<axes[a].values.length; v++){
				final int[] point=initial();
				point[a]=v;
				cases.add(policy("grid"+a+"value"+v, point));
			}
		}
		return cases;
	}

	/** One adjacent value in one coordinate; null at that coordinate's boundary. */
	public int[] neighbor(final int[] center, final int axis, final int direction){
		validate(center);
		require(axis>=0 && axis<axes.length && (direction==-1 || direction==1), "A neighbor changes one declared axis by one grid step");
		final int next=center[axis]+direction;
		if(next<0 || next>=axes[axis].values.length){return null;}
		final int[] result=center.clone();
		result[axis]=next;
		return result;
	}

	/** Materializes the complete member vector; all unchanged values stay explicit. */
	public NcrnaJointPolicy policy(final String name, final int[] point){
		validate(point);
		final Member[] members=new Member[start.size()];
		for(int member=0; member<members.length; member++){
			final Member base=start.member(member);
			int cap=base.topN;
			int patience=base.patience;
			float floor=base.floor, topFrac=base.topFrac, qFrac=base.qFrac;
			for(int a=0; a<axes.length; a++){
				if(axes[a].member!=member){continue;}
				final double v=axes[a].values[point[a]];
				switch(axes[a].knob){
					case 0: cap=(int)v; break;
					case 1: floor=(float)v; break;
					case 2: topFrac=(float)v; break;
					case 3: qFrac=(float)v; break;
					case 4: patience=(int)v; break;
					default: throw new AssertionError("Unknown validated grid control");
				}
			}
			members[member]=new Member(base.name, true, cap, floor, base.minHits, topFrac, qFrac, patience);
		}
		return new NcrnaJointPolicy(name, members);
	}

	/** Small coordinate key, not a sequence key or a per-record allocation. */
	public String key(final int[] point){
		validate(point);
		final StringBuilder key=new StringBuilder();
		for(final int value:point){key.append(value).append('_');}
		return key.toString();
	}

	private void validate(final int[] point){
		require(point!=null && point.length==axes.length, "A search point needs one position per axis");
		for(int i=0; i<point.length; i++){
			require(point[i]>=0 && point[i]<axes[i].values.length, "Search position is outside its finite grid");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Selection           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the best admissible improving neighbor, or -1 on a local stop.
	 * maxLost applies to actual root-relative reference losses, independently of TP.
	 * A negative maxLost imposes no paired-locus constraint. Candidate order is the
	 * final deterministic tie breaker; callers generate axes in declared grid order.
	 */
	public static int choose(final Measurement center, final Measurement[] neighbors,
			final long targetTp, final long maxLost){
		return choose(center, neighbors, targetTp, maxLost, -1, -1);
	}

	/** Coverage budgets are independent root-relative admissibility checks. Gains
	 * elsewhere never refund a last-detection loss; disabled guards need no counts. */
	public static int choose(final Measurement center, final Measurement[] neighbors,
			final long targetTp, final long maxLost, final long maxLastLoss, final long maxMemberLastLoss){
		require(center!=null && neighbors!=null && targetTp>=0 && maxLost>=-1, "Selection needs a measured center and explicit objective");
		require(maxLastLoss>=-1 && maxMemberLastLoss>=-1, "Coverage budgets use -1 only to disable protection");
		int best=-1;
		for(int i=0; i<neighbors.length; i++){
			final Measurement next=neighbors[i];
			if(next==null || !constraintsMet(next, maxLost, maxLastLoss, maxMemberLastLoss)){continue;}
			if(!improves(center, next, targetTp)){continue;}
			if(best<0 || better(center, next, neighbors[best], targetTp)){best=i;}
		}
		return best;
	}

	/** Does not impose the TP floor on early recovery moves; improves() owns that rule. */
	public static boolean constraintsMet(final Measurement point, final long maxLost,
			final long maxLastLoss, final long maxMemberLastLoss){
		require(point!=null && maxLost>=-1 && maxLastLoss>=-1 && maxMemberLastLoss>=-1,
			"Admissibility needs a measured point and explicit loss budgets");
		if(maxLastLoss>=0 || maxMemberLastLoss>=0){
			require(point.coverage!=null, "Enabled coverage protection requires validated input counts");
		}
		if(maxLost>=0 && point.lost>maxLost){return false;}
		if(maxLastLoss>=0 && point.coverage.lastLoss>maxLastLoss){return false;}
		if(maxMemberLastLoss>=0){
			require(point.coverage.memberLastLoss.length>0, "Member protection requires every declared reference member");
			for(final long count:point.coverage.memberLastLoss){if(count>maxMemberLastLoss){return false;}}
		}
		return true;
	}

	private static boolean improves(final Measurement center, final Measurement next, final long target){
		if(dominates(next, center)){return true;}
		return center.tp<target ? next.tp>center.tp : next.tp>=target && next.cost<center.cost;
	}

	private static boolean better(final Measurement center, final Measurement a, final Measurement b, final long target){
		final boolean ad=dominates(a, center), bd=dominates(b, center);
		if(ad!=bd){return ad;}
		if(!ad){
			final double ar, br;
			if(center.tp<target){
				// Non-dominating recall gains necessarily have positive extra cost.
				assert(a.cost>center.cost && b.cost>center.cost) : "Zero-cost gains are selected as dominating moves before gain/cost division";
				ar=(a.tp-center.tp)/(double)(a.cost-center.cost);
				br=(b.tp-center.tp)/(double)(b.cost-center.cost);
			}else{
				assert(a.tp<center.tp && b.tp<center.tp) : "Loss-free savings are dominating moves, not loss-normalized ratios";
				ar=(center.cost-a.cost)/(double)(center.tp-a.tp);
				br=(center.cost-b.cost)/(double)(center.tp-b.tp);
			}
			if(ar!=br){return ar>br;}
		}
		if(a.tp!=b.tp){return a.tp>b.tp;}
		if(a.cost!=b.cost){return a.cost<b.cost;}
		if(a.lost!=b.lost){return a.lost<b.lost;}
		return a.fp<b.fp;
	}

	private static boolean dominates(final Measurement a, final Measurement b){
		return a.tp>=b.tp && a.cost<=b.cost && (a.tp>b.tp || a.cost<b.cost);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Types            ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Measurement{
		public Measurement(final long tp_, final long fp_, final long cost_, final long lost_){
			this(tp_, fp_, cost_, lost_, null);
		}
		public Measurement(final long tp_, final long fp_, final long cost_, final long lost_,
				final NcrnaSearchCoverage.Counts coverage_){
			require(tp_>=0 && fp_>=0 && cost_>=0 && lost_>=0, "Measured counts and root-relative losses must be nonnegative");
			tp=tp_;
			fp=fp_;
			cost=cost_;
			lost=lost_;
			coverage=coverage_;
		}
		public final long tp, fp, cost, lost;
		public final NcrnaSearchCoverage.Counts coverage;
	}

	public static final class Axis{
		public Axis(final int member_, final int knob_, final double[] values_){
			require(member_>=0 && knob_>=0 && knob_<KNOBS.length && values_!=null && values_.length>0, "An axis needs a member, known control and finite choices");
			member=member_;
			knob=knob_;
			values=values_.clone();
			for(int i=0; i<values.length; i++){
				final double value=values[i];
				require(Double.isFinite(value) && value>=0, "Grid values must be finite and nonnegative");
				if(knob==0 || knob==4){
					require(value>=(knob==0 ? 1 : 0) && value<=Integer.MAX_VALUE && value==(int)value,
						"Caps and patience must be exact bounded integers; legacy -1 is baseline-only");
				}else{
					values[i]=(float)value;
					require(Float.isFinite((float)value), "Grid cutoff overflows the caller's float type");
				}
				require(i==0 || values[i]>values[i-1], "Grid values must remain strictly increasing in the caller's actual numeric type");
			}
		}
		private int index(final double value){
			for(int i=0; i<values.length; i++){if(values[i]==value){return i;}}
			return -1;
		}
		public final int member, knob;
		private final double[] values;
	}

	private static double value(final Member member, final int knob){
		switch(knob){
			case 0: return member.topN;
			case 1: return member.floor;
			case 2: return member.topFrac;
			case 3: return member.qFrac;
			case 4: return member.patience;
			default: throw new IllegalArgumentException("Unknown grid control");
		}
	}
	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}
	private final NcrnaJointPolicy start;
	private final Axis[] axes;
	private final int[] initial;
	private static final String[] KNOBS={"topn", "floor", "topfrac", "qfrac", "patience"};
}
