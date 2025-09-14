package bin;

import java.util.ArrayList;
import structures.IntHashSet;

import shared.Tools;

/**
 * Abstract base class for bin refinement strategies.
 * Accepts a potentially impure bin and returns null for no change,
 * or a list of refined bins if beneficial splits are found.
 * 
 * @author Brian Bushnell & UMP45
 * @date June 20, 2025
 */
abstract class AbstractRefiner extends BinObject {

	/**
	 * Attempts to refine/split the given bin.
	 * @param input Potentially impure cluster to analyze
	 * @return null if no refinement recommended, or ArrayList of 2+ bins if split beneficial
	 */
	abstract ArrayList<Bin> refine(Bin input);
	
	abstract ArrayList<IntHashSet> refineToIntSets(Bin input);

	/**
	 * Validates that a proposed split actually improves cluster quality.
	 * Override this for custom validation logic.
	 */
	protected boolean isSplitBeneficial(Bin original, ArrayList<Bin> splits) {
		if (splits == null || splits.size() < 2) return false; //Require at least 2 splits to be meaningful

		// Basic sanity checks
		long totalSize = 0;
		for (Bin bin : splits) {
			if (bin.numContigs() == 0) return false; //Reject empty splits
			totalSize += bin.size(); //Accumulate total size across all splits
		}

		// Conservation of mass
		if (totalSize != original.size()) return false;

		// Each split should be reasonably sized
		for (Bin bin : splits) {
			if (bin.size() < original.size() * 0.1f) return false; // No tiny fragments
		}

		return true;
	}

	/**
	 * Factory method to create a refiner instance using default type and parameters.
	 * @param oracle Truth data provider for supervised learning and validation
	 * @return Configured refiner instance using DEFAULT_TYPE with default parameters
	 */
	public static AbstractRefiner makeRefiner(Oracle oracle){
		return makeRefiner(oracle, DEFAULT_TYPE, null);
	}

	/**
	 * Factory method to create a refiner instance with specified type and default 
	 * parameters.
	 * @param oracle Truth data provider for supervised learning and validation
	 * @param type Refinement algorithm type (CRYSTAL, GRAPH, EVIDENCE, or ENSEMBLE)
	 * @return Configured refiner instance using default parameters for the specified type
	 */
	public static AbstractRefiner makeRefiner(Oracle oracle, int type){
		return makeRefiner(oracle, type, null);
	}
	
	/**
	 * Primary factory method for creating configured refiner instances with custom 
	 * parameters.
	 * 
	 * Instantiates the appropriate concrete refiner implementation based on the 
	 * specified type:
	 * - CRYSTAL: CrystalChamber using crystallization-based clustering
	 * - GRAPH: GraphRefiner using assembly graph community detection
	 * - EVIDENCE: EvidenceRefiner using DBSCAN clustering on coverage/composition 
	 *   features
	 * - ENSEMBLE: EnsembleRefiner combining multiple methods with consensus voting
	 * 
	 * @param oracle Truth data provider for supervised learning, validation, and ground 
	 *               truth comparison
	 * @param type Algorithm type constant (CRYSTAL, GRAPH, EVIDENCE, or ENSEMBLE)
	 * @param params Type-specific parameter object or null for defaults
	 * @return Fully configured refiner instance ready for bin processing
	 * @throws RuntimeException if type is not recognized
	 */
	public static AbstractRefiner makeRefiner(Oracle oracle, int type, RefinerParams params){
		if(type==CRYSTAL){return new CrystalChamber(oracle);}
		if(type==GRAPH){return new GraphRefiner(oracle, params!=null ? (GraphRefinerParams)params : new GraphRefinerParams());} //Possible bug: unsafe cast without type validation
		if(type==EVIDENCE){return new EvidenceRefiner(oracle, params!=null ? (EvidenceRefinerParams)params : new EvidenceRefinerParams());} //Possible bug: unsafe cast without type validation
		if(type==ENSEMBLE){return new EnsembleRefiner(oracle, params!=null ? (EnsembleRefinerParams)params : new EnsembleRefinerParams());} //Possible bug: unsafe cast without type validation
		throw new RuntimeException("Unknown refiner type: "+type);
	}
	
	public static int findType(String s) {
		int idx=Tools.find(s, types);
		assert(idx>=0) : "Can't find type "+s;
		return idx;
	}

	/** Algorithm type constants for supported refinement strategies */
	public static final int CRYSTAL=0, GRAPH=1, EVIDENCE=2, ENSEMBLE=3;
	
	/** Human-readable names corresponding to algorithm type constants */
	public static final String types[]={"CRYSTAL", "GRAPH", "EVIDENCE", "ENSEMBLE"};
	
	/** Default refinement algorithm type used when no specific type is requested */
	public static int DEFAULT_TYPE=CRYSTAL;
	
	/**
	 * Abstract base class for refiner-specific parameter configuration objects.
	 * 
	 * Provides a common interface for passing algorithm-specific parameters to refiner 
	 * instances
	 * while maintaining type safety through concrete subclasses. Each refiner 
	 * implementation
	 * defines its own parameter class extending this base to specify tunable algorithm 
	 * parameters
	 * such as clustering thresholds, iteration limits, and random seeds.
	 * 
	 * Supports deep copying to enable parameter mutation experiments without affecting 
	 * original configurations.
	 */
	public static abstract class RefinerParams{
		/** Creates a deep copy of this parameter configuration */
		public abstract RefinerParams copy();
	}
	
	/**
	 * Configuration parameters for graph-based bin refinement using community detection 
	 * algorithms.
	 * 
	 * Controls the behavior of GraphRefiner which analyzes assembly graph connectivity 
	 * patterns
	 * to identify natural breakpoints between different organisms that were incorrectly 
	 * binned together.
	 * Uses iterative community detection algorithms such as Louvain or Label Propagation 
	 * to partition
	 * contigs based on graph edge weights derived from paired-end links, shared k-mer 
	 * content, or other
	 * connectivity evidence.
	 */
	public static class GraphRefinerParams extends RefinerParams{
		/** Minimum edge weight threshold for including connections in community detection (0.0-1.0) */
		public float minEdgeWeight=0.3f;
		
		/** Maximum iterations for convergence in iterative community detection algorithms */
		public int maxIterations=50;
		
		/** Random seed for reproducible community detection results */
		public long seed=42;
		
		/** Default constructor using standard parameter values */
		public GraphRefinerParams(){}
		
		/**
		 * Constructor with custom parameter values.
		 * @param minEdgeWeight Minimum edge weight threshold for community detection
		 * @param maxIterations Maximum iterations for algorithm convergence
		 * @param seed Random seed for reproducible results
		 */
		public GraphRefinerParams(float minEdgeWeight, int maxIterations, long seed){
			this.minEdgeWeight=minEdgeWeight;
			this.maxIterations=maxIterations;
			this.seed=seed;
		}
		
		@Override
		public RefinerParams copy(){
			return new GraphRefinerParams(minEdgeWeight, maxIterations, seed);
		}
	}
	
	/**
	 * Configuration parameters for evidence-based bin refinement using DBSCAN density clustering.
	 * 
	 * Controls the behavior of EvidenceRefiner which applies DBSCAN clustering to 
	 * multi-dimensional
	 * feature vectors derived from contig properties such as coverage profiles, 
	 * tetranucleotide frequencies,
	 * GC content, and other compositional signatures. Identifies dense clusters of similar 
	 * contigs
	 * that likely originate from the same organism, revealing contamination or chimeric bins.
	 * 
	 * DBSCAN parameters determine the sensitivity and specificity of cluster detection in 
	 * high-dimensional
	 * feature space where euclidean distance represents biological similarity between 
	 * contigs.
	 */
	public static class EvidenceRefinerParams extends RefinerParams{
		/** DBSCAN epsilon parameter - maximum distance between contigs in same cluster */
		public float epsilon=0.4f;
		
		/** DBSCAN minPoints parameter - minimum contigs required to form dense cluster */
		public int minPoints=3;
		
		/** Minimum cluster size threshold for accepting a refined bin */
		public int minClusterSize=2;
		
		/** Random seed for reproducible clustering results */
		public long seed=123;
		
		/** Default constructor using standard DBSCAN parameter values */
		public EvidenceRefinerParams(){}
		
		/**
		 * Constructor with custom DBSCAN parameter values.
		 * @param epsilon Maximum euclidean distance between contigs in same cluster
		 * @param minPoints Minimum contigs required to form a dense cluster
		 * @param minClusterSize Minimum cluster size for accepting refined bins
		 * @param seed Random seed for reproducible results
		 */
		public EvidenceRefinerParams(float epsilon, int minPoints, int minClusterSize, long seed){
			this.epsilon=epsilon;
			this.minPoints=minPoints;
			this.minClusterSize=minClusterSize;
			this.seed=seed;
		}
		
		@Override
		public RefinerParams copy(){
			return new EvidenceRefinerParams(epsilon, minPoints, minClusterSize, seed);
		}
	}
	
	public static class EnsembleRefinerParams extends RefinerParams{
		/** Minimum consensus threshold for accepting ensemble splits (0.0-1.0) */
		public float consensusThreshold=0.6f;
		
		/** Minimum number of refinement methods that must agree on a split */
		public int minMethodsAgreeing=2;
		
		/** Random seed for reproducible ensemble results */
		public long seed=999;
		
		/** Parameter configuration for graph-based refinement component */
		public GraphRefinerParams graphParams;
		
		/** Parameter configuration for evidence-based refinement component */
		public EvidenceRefinerParams evidenceParams;
		
		/** Default constructor initializing component parameters with standard values */
		public EnsembleRefinerParams(){
			this.graphParams=new GraphRefinerParams();
			this.evidenceParams=new EvidenceRefinerParams();
		}
		
		/**
		 * Constructor with custom consensus parameters using default component parameters.
		 * @param consensusThreshold Minimum consensus threshold for accepting splits
		 * @param minMethodsAgreeing Minimum number of methods that must agree
		 * @param seed Random seed for reproducible results
		 */
		public EnsembleRefinerParams(float consensusThreshold, int minMethodsAgreeing, long seed){
			this();
			this.consensusThreshold=consensusThreshold;
			this.minMethodsAgreeing=minMethodsAgreeing;
			this.seed=seed;
		}
		
		@Override
		public RefinerParams copy(){
			EnsembleRefinerParams copy=new EnsembleRefinerParams(consensusThreshold, minMethodsAgreeing, seed);
			copy.graphParams=(GraphRefinerParams)graphParams.copy();
			copy.evidenceParams=(EvidenceRefinerParams)evidenceParams.copy();
			return copy;
		}
	}
	
}
