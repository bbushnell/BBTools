package prot;

import java.util.HashMap;
import java.util.List;

import fileIO.ByteFile;
import parse.LineParser1;

/**
 * Frozen protein-family search inputs for the FASTA-in MAG-QC seam.
 *
 * <p>The representatives remain the full FASTA database: a query may best hit
 * a representative outside the vector's top-N families, in which case the
 * shared CacheBuilder reduction deliberately excludes it.  {@link #repRank}
 * therefore contains only the family-list ranks below {@code topN}, exactly as
 * CacheBuilder does.</p>
 */
public final class FamilySearchResources {

	/** Frozen v1 family-assignment parameters shared by MagQCTool and ProteinSearch. */
	public static final int FROZEN_V1_K=5;
	public static final int FROZEN_V1_MIN_SEED_HITS=1;
	public static final double FROZEN_V1_EVALUE=0.001;
	public static final double FROZEN_V1_MIN_PIDENT=30;
	public static final double FROZEN_V1_MIN_COVERAGE=0.7;
	public static final int FROZEN_V1_MAX_TARGET_SEQS=25;
	public static final boolean FROZEN_V1_REDUCED_SEED=true;

	public final List<ProteinSequence> representatives;
	public final HashMap<String, Integer> repRank;

	private FamilySearchResources(final List<ProteinSequence> representatives,
			final HashMap<String, Integer> repRank){
		this.representatives=representatives;
		this.repRank=repRank;
	}

	/** Applies the single frozen v1 parameter set used by all production family searches. */
	public static void applyFrozenV1(final ProteinSearcher searcher){
		if(searcher==null){throw new RuntimeException("ProteinSearcher is required.");}
		searcher.k=FROZEN_V1_K;
		searcher.minSeedHits=FROZEN_V1_MIN_SEED_HITS;
		searcher.evalueCutoff=FROZEN_V1_EVALUE;
		searcher.minPident=FROZEN_V1_MIN_PIDENT;
		searcher.minCoverage=FROZEN_V1_MIN_COVERAGE;
		searcher.maxTargetSeqs=FROZEN_V1_MAX_TARGET_SEQS;
		searcher.reducedSeed=FROZEN_V1_REDUCED_SEED;
	}

	/** Loads representatives and their top-N vector-rank map from frozen files. */
	public static FamilySearchResources load(final String representativesFasta,
			final String familyList, final int topN){
		if(representativesFasta==null || familyList==null){throw new RuntimeException("Representative FASTA and family list are required.");}
		if(topN<1){throw new RuntimeException("topN must be positive: "+topN);}
		final List<ProteinSequence> representatives=ProteinSearch.readFasta(representativesFasta);
		final HashMap<String, Integer> repRank=new HashMap<String, Integer>(Math.max(16, topN*2));
		final LineParser1 lp=new LineParser1((byte)'\t');
		final ByteFile bf=ByteFile.makeByteFile(familyList, true);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			lp.set(line);
			if(lp.terms()<2){throw new RuntimeException("Malformed family-list row: "+new String(line));}
			final int rank=lp.parseInt(0);
			if(rank<0){throw new RuntimeException("Negative family rank: "+rank);}
			if(rank>=topN){continue;}
			final String rep=lp.parseString(1);
			if(repRank.put(rep, Integer.valueOf(rank))!=null){throw new RuntimeException("Duplicate family representative: "+rep);}
		}
		bf.close();
		if(repRank.isEmpty()){throw new RuntimeException("No top-"+topN+" family ranks loaded from "+familyList);}
		return new FamilySearchResources(representatives, repRank);
	}
}
