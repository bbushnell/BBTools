package bin;

import java.util.ArrayList;

import shared.Tools;
import structures.IntHashSet;

abstract class AbstractRefiner extends BinObject {

	/**
	 * Attempts to refine/split the given bin into multiple cleaner bins.
	 * Analyzes the bin for potential contamination or chimeric assemblies and
	 * returns refined bins if beneficial splits are identified.
	 *
	 * @param input Potentially impure bin to analyze for refinement opportunities
	 * @return null if no refinement recommended, or ArrayList of 2+ bins if split beneficial
	 */
	abstract ArrayList<Bin> refine(Bin input);
	
	abstract ArrayList<IntHashSet> refineToIntSets(Bin input);

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

	public static AbstractRefiner makeRefiner(Oracle oracle){
		return makeRefiner(oracle, DEFAULT_TYPE, null);
	}

	public static AbstractRefiner makeRefiner(Oracle oracle, int type){
		return makeRefiner(oracle, type, null);
	}
	
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

	public static final int CRYSTAL=0, GRAPH=1, EVIDENCE=2, ENSEMBLE=3;
	
	public static final String types[]={"CRYSTAL", "GRAPH", "EVIDENCE", "ENSEMBLE"};
	
	public static int DEFAULT_TYPE=CRYSTAL;
	
	public static abstract class RefinerParams{
		public abstract RefinerParams copy();
	}
	
	public static class GraphRefinerParams extends RefinerParams{
		public float minEdgeWeight=0.3f;
		
		public int maxIterations=50;
		
		public long seed=42;
		
		public GraphRefinerParams(){}
		
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
	
	public static class EvidenceRefinerParams extends RefinerParams{
		public float epsilon=0.4f;
		
		public int minPoints=3;
		
		public int minClusterSize=2;
		
		public long seed=123;
		
		public EvidenceRefinerParams(){}
		
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
		public float consensusThreshold=0.6f;
		
		public int minMethodsAgreeing=2;
		
		public long seed=999;
		
		public GraphRefinerParams graphParams;
		
		public EvidenceRefinerParams evidenceParams;
		
		public EnsembleRefinerParams(){
			this.graphParams=new GraphRefinerParams();
			this.evidenceParams=new EvidenceRefinerParams();
		}
		
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
