package prot;

import java.util.HashMap;
import java.util.Map;

/** Frozen v1.4 shared inputs computed from observations, never simulation truth.
 * Requiredness counts each found family once and counts its extra copies as zero.
 * The denominator is therefore the number of observed tracked protein copies.
 * @author Brian Bushnell, Yoimiya
 */
final class MagQCV14Features {

	/** Binds the population resource and explicit marker-population declarations once. */
	MagQCV14Features(final MagQCRequirednessTable table, final Map<String,int[]> markers,
			final double minimumPresence, final double minimumMulticopyCarriers){
		if(table==null || markers==null || !fraction(minimumPresence) || !fraction(minimumMulticopyCarriers)){
			throw new IllegalArgumentException("v1.4 inputs require a population table, marker bindings and thresholds in [0,1]");
		}
		families=table.families;
		global=presence(table.global); bacteria=presence(table.bacteria); archaea=presence(table.archaea);
		phylumPresence=new HashMap<String,double[]>();
		for(String phylum:table.phylumNames()){
			phylumPresence.put(phylum, presence(table.phylum(phylum)));
		}
		markerRanks=new HashMap<String,int[]>();
		for(Map.Entry<String,int[]> entry:markers.entrySet()){
			final int[] ranks=entry.getValue();
			if(entry.getKey()==null || entry.getKey().isEmpty() || ranks==null || ranks.length==0){
				throw new IllegalArgumentException("Marker declarations need a named phylum and family ranks");
			}
			final boolean[] seen=new boolean[families];
			for(int rank:ranks){
				if(rank<0 || rank>=families || seen[rank]){throw new IllegalArgumentException("Invalid or repeated marker family rank: "+rank);}
				seen[rank]=true;
			}
			markerRanks.put(entry.getKey(), ranks.clone());
		}
		multicopyInverseMean=new double[families];
		for(int rank=0; rank<families; rank++){
			final long carriers=table.global.carriers[rank];
			final long multiple=carriers-table.global.singletons[rank];
			if(carriers>0 && global[rank]>=minimumPresence && multiple/(double)carriers>=minimumMulticopyCarriers){
				multicopyInverseMean[rank]=1/table.global.carrierMean(rank);
			}
		}
	}

	/** Fills ten appended inputs; absent counts and unknown phyla have explicit finite fallbacks. */
	void compute(final MagQCPreparedBin bin, final double[] out){
		if(bin==null || bin.families.length!=families || out==null || out.length!=WIDTH){
			throw new IllegalArgumentException("v1.4 feature buffers must match the frozen family and input widths");
		}
		bin.validate();
		out[0]=Math.log(bin.contigCount())/LOG_TWO;
		out[1]=Math.log(bin.n50())/LOG_TWO;
		final boolean classified="classified".equals(bin.status);
		final double[] predicted=classified ? phylumPresence.get(bin.phylum) : null;
		final double[] phylum=predicted==null ? global : predicted;
		long copies=0;
		double universal=0, bacterial=0, archaeal=0, local=0, multipleSum=0;
		int multipleFound=0;
		for(int rank=0; rank<families; rank++){
			final int count=bin.families[rank];
			if(count<0){throw new IllegalArgumentException("Negative observed family count at "+rank);}
			if(count==0){continue;}
			copies+=count;
			universal+=global[rank]; bacterial+=bacteria[rank]; archaeal+=archaea[rank]; local+=phylum[rank];
			if(multicopyInverseMean[rank]>0){multipleSum+=count*multicopyInverseMean[rank]; multipleFound++;}
		}
		// CallGenes reports these disjoint feature classes. Unknown anticodons already belong to trna.
		final long calls=(long)bin.agg.cds+bin.agg.r16+bin.agg.r23+bin.agg.r5+bin.agg.rother+bin.agg.trna;
		out[2]=ratio(calls, copies);
		final int[] ranks=classified ? markerRanks.get(bin.phylum) : null;
		long markerCopies=copies;
		if(ranks!=null){markerCopies=0; for(int rank:ranks){markerCopies+=bin.families[rank];}}
		out[3]=ratio(calls, markerCopies);
		out[4]=ratio(universal, copies); out[5]=ratio(bacterial, copies);
		out[6]=ratio(archaeal, copies); out[7]=ratio(local, copies);
		out[8]=ratio(multipleSum, multipleFound);
		int anticodons=0;
		for(int count:bin.agg.anti){if(count>0){anticodons++;}}
		assert(bin.agg.anti.length==64) : "Anticodon fraction counts the 64 structural anticodons, excluding unknown";
		out[9]=anticodons/64.0;
		for(double value:out){
			if(!Double.isFinite(value) || value<0){throw new IllegalStateException("Nonfinite or negative v1.4 input for "+bin.id);}
		}
	}

	/** Precomputes presence fractions; zero-denominator populations are not valid resources. */
	private static double[] presence(final MagQCRequirednessTable.Population population){
		assert(population.genomes>0) : "Frozen population denominators must be positive to define presence";
		final double[] values=new double[population.carriers.length];
		for(int rank=0; rank<values.length; rank++){population.validate(rank); values[rank]=population.presence(rank);}
		return values;
	}

	/** A missing denominator carries no ratio evidence and is encoded as zero, never infinity. */
	private static double ratio(final double numerator, final long denominator){return denominator==0 ? 0 : numerator/denominator;}
	/** Validates tunable corpus thresholds before constructing a mask. */
	private static boolean fraction(final double value){return Double.isFinite(value) && value>=0 && value<=1;}

	private final int families;
	private final double[] global, bacteria, archaea, multicopyInverseMean;
	private final HashMap<String,double[]> phylumPresence;
	private final HashMap<String,int[]> markerRanks;
	static final int WIDTH=10;
	static final String COLUMNS="log2_contigs\tlog2_n50_bp\tcalled_to_tracked\tcalled_to_phylum_markers"
		+"\trequiredness_global\trequiredness_bacteria\trequiredness_archaea\trequiredness_phylum"
		+"\tmulticopy_ratio\tanticodon_fraction";
	private static final double LOG_TWO=Math.log(2);
}
