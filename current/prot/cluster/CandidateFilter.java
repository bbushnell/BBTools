package prot.cluster;

import dna.AminoAcid;
import prot.Blosum62;
import prot.ConsensusScalars;

/** Immutable cheap screens applied after ranking and before alignment.
 * Missing hydrophobicity leaves that screen unavailable; the hit/alignment gates
 * still apply. A Context belongs to one query worker and is reused.
 * @author Keqing
 */
public final class CandidateFilter {

	public CandidateFilter(int minimumHits,int greedyHits,double maximumHydroDelta,double greedyHydroDelta){
		this(minimumHits,greedyHits,maximumHydroDelta,greedyHydroDelta,false);
	}

	/** Optional reassignment flex uses the same occupancy rule; legacy callers keep a fixed floor. */
	public CandidateFilter(int minimumHits,int greedyHits,double maximumHydroDelta,double greedyHydroDelta,boolean passHitsFlex){
		if(minimumHits<0 || greedyHits<-1 || !Double.isFinite(maximumHydroDelta) ||
			!Double.isFinite(greedyHydroDelta) || maximumHydroDelta<0 || maximumHydroDelta>9 ||
			greedyHydroDelta<0 || greedyHydroDelta>9){
			throw new IllegalArgumentException("Nonnegative hit floors (greedy -1 means flex) and finite hydro caps in [0,9] required");
		}
		this.minimumHits=minimumHits; this.greedyHits=greedyHits;
		this.maximumHydroDelta=maximumHydroDelta; this.greedyHydroDelta=greedyHydroDelta;
		this.passHitsFlex=passHitsFlex;
	}

	/** Both released widths use the same absolute occupancy formula, without rescaling hits. */
	public void validateWidth(int width){
		if(width<1 || ((greedyHits<0 || passHitsFlex) && width!=512 && width!=1024)){
			throw new IllegalArgumentException("Flexible greedy or reassignment hit rules require width=512 or1024; choose fixed hit floors for other widths");
		}
	}

	public int requiredHits(boolean greedy,int occupiedBins){
		return requiredHits(greedy,occupiedBins,false);
	}

	/** earlyGreedy is captured once by the query, never read from live progress here. */
	public int requiredHits(boolean greedy,int occupiedBins,boolean earlyGreedy){
		if(earlyGreedy && (!greedy || greedyHits>=0)){
			throw new IllegalArgumentException("The early ramp applies only to flexible greedy screening");
		}
		if(!greedy && !passHitsFlex){return minimumHits;}
		if(greedy && greedyHits>=0){return greedyHits;}
		if(occupiedBins<0 || occupiedBins>1024){throw new IllegalArgumentException("Flex occupancy outside [0,1024]: "+occupiedBins);}
		return 2+(occupiedBins+200)/512-(earlyGreedy ? 1 : 0);
	}

	/** Hit failures take precedence so the two screening counts are disjoint. */
	public int rejection(boolean greedy,int hits,int occupiedBins,Context query,Family family){
		return rejection(greedy,hits,occupiedBins,query,family,false);
	}

	/** Uses one captured ramp state for both ranked candidates and competing seeds. */
	public int rejection(boolean greedy,int hits,int occupiedBins,Context query,Family family,boolean earlyGreedy){
		if(hits<0 || query==null || family==null){throw new IllegalArgumentException("Candidate screen needs nonnegative hits, query and family");}
		if(hits<requiredHits(greedy,occupiedBins,earlyGreedy)){return LOW_HITS;}
		final double reference=family.scalarsUnsafe()[ConsensusScalars.HYDRO];
		if(!Double.isFinite(reference) || !Double.isFinite(query.hydrophobicity)){
			throw new IllegalArgumentException("Hydrophobicity must be finite even when unavailable");
		}
		final double cap=greedy ? greedyHydroDelta : maximumHydroDelta;
		if(cap<9 && query.defined && (family.scalarMask()&(1<<ConsensusScalars.HYDRO))!=0 &&
			Math.abs(query.hydrophobicity-reference)>cap){
			return HYDRO;
		}
		return PASS;
	}

	/** Retains old API behavior without nonfinite disabling sentinels. */
	public static CandidateFilter disabled(){return DISABLED;}

	public static final class Context {
		/** Same residue order and native Kyte-Doolittle constants as ConsensusScalars,
		 * without calculating unrelated pI values for the cheap pre-alignment screen. */
		public void reset(byte[] encoded){
			if(encoded==null){throw new IllegalArgumentException("Missing encoded query");}
			double sum=0; int valid=0;
			for(int i=0; i<encoded.length; i++){
				final int aa=encoded[i];
				if(aa>=0 && aa<20){sum+=RESIDUE_HYDRO[aa]; valid++;}
				else if(aa!=Blosum62.X_CODE){throw new IllegalArgumentException("Invalid encoded query residue at "+i+": "+aa);}
			}
			defined=valid>0; hydrophobicity=defined ? sum/valid : 0;
		}
		public double hydrophobicity(){return hydrophobicity;}
		public boolean defined(){return defined;}
		private double hydrophobicity;
		private boolean defined;
	}

	/** Read the existing native scalar implementation once; do not duplicate its biochemical table. */
	private static double[] residueHydrophobicities(){
		final double[] values=new double[20],scratch=new double[ConsensusScalars.SCALAR_COUNT];
		final byte[] one=new byte[1];
		for(int aa=0; aa<values.length; aa++){
			one[0]=AminoAcid.numberToAcid[aa];
			final int mask=ConsensusScalars.computeScalarMask(one,1,scratch,0);
			if((mask&(1<<ConsensusScalars.HYDRO))==0 || !Double.isFinite(scratch[ConsensusScalars.HYDRO])){
				throw new IllegalStateException("Native hydrophobicity is unavailable for standard encoded residue "+aa);
			}
			values[aa]=scratch[ConsensusScalars.HYDRO];
		}
		return values;
	}

	public final int minimumHits,greedyHits;
	public final boolean passHitsFlex;
	public final double maximumHydroDelta,greedyHydroDelta;
	public static final int PASS=0,LOW_HITS=1,HYDRO=2;
	private static final double[] RESIDUE_HYDRO=residueHydrophobicities();
	private static final CandidateFilter DISABLED=new CandidateFilter(0,0,9,9);
}
