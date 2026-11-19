package prot;

/** Explicit fragment-only gates on existing full-query, reference-local linear alignment.
 * No tail/member gates are inherited. Thresholds are supplied by a separately bound policy;
 * constructing this scorer makes no claim that those thresholds have been calibrated.
 * Both identity and full-query-normalized BLOSUM must meet their inclusive floors.
 * @author Qiqi
 */
public final class FragmentAlignment{

	public FragmentAlignment(final double identityFloor_, final double blosumFloor_){
		if(!Double.isFinite(identityFloor_) || identityFloor_<0 || identityFloor_>100 ||
				!Double.isFinite(blosumFloor_) || blosumFloor_<0){
			throw new IllegalArgumentException("Fragment identity and BLOSUM floors must be explicit finite nonnegative values; identity<=100");
		}
		identityFloor=identityFloor_; blosumFloor=blosumFloor_;
	}
	public Worker newWorker(){return new Worker();}

	/** Scratch belongs to one consuming thread. The supplied output retains rejected metrics
	 * as well as accepted metrics; caller-owned observations cannot be shared concurrently. */
	public final class Worker{
		public Reason measure(final ProteinSequence query, final ProteinSequence reference, final FamilyAlignmentMetric.Observation out){
			if(query==null || reference==null || out==null || query.length()<40 || query.length()>100 ||
					reference.length()<1 || reference.length()>100){
				throw new IllegalArgumentException("Fragment alignment requires40..100 query residues and1..100 retained reference residues");
			}
			for(byte aa:query.enc){
				if(!Blosum62.isStandard(aa)){throw new IllegalArgumentException("Fragment queries must be uninterrupted canonical residues");}
			}
			alignment.measure(query, reference, ProteinSearcher.selfScore(query), out);
			assert(out.measured() && out.queryStart==0 && out.queryEnd==query.length()-1 &&
				out.denominatorSelf==out.fullQuerySelf && out.pairedColumns<=out.columns) :
				"Fragment gates require all query residues, full-query self normalization and actual traceback pairs";
			if(out.pairedColumns<MINIMUM_PAIRED){return Reason.PAIRED_COLUMNS_SHORT;}
			if(out.identityPct()<identityFloor){return Reason.IDENTITY_LOW;}
			if(out.blosum()<blosumFloor){return Reason.BLOSUM_LOW;}
			return Reason.ACCEPT;
		}
		private final FamilyAlignmentMetric.Worker alignment=FamilyAlignmentMetric.LEGACY.newWorker();
	}
	public final double identityFloor, blosumFloor;
	public static final int MINIMUM_PAIRED=28;
	public enum Reason{ACCEPT, PAIRED_COLUMNS_SHORT, IDENTITY_LOW, BLOSUM_LOW}
}
