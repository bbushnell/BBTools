package align2;

import java.nio.charset.StandardCharsets;

import dna.Data;
import stream.CustomHeader;
import stream.Read;
import stream.SamLine;
import structures.ByteBuilder;

/** Serializes one mapped mate and its reciprocal pair evidence. */
public final class NeuralMapqPairedFeatureRow {

	private NeuralMapqPairedFeatureRow(){}

	public static ByteBuilder header(){
		final ByteBuilder bb=new ByteBuilder(2048);
		bb.append("#schema\t").append(NeuralMapqPairedFeatureSchema.SCHEMA_NAME).nl();
		bb.append("#oracle_tolerance\t").append(NeuralMapqPairedFeatureSchema.ORACLE_TOLERANCE).nl();
		bb.append("numeric_id\tanchor_pairnum\tread_name\terror\tanchor_correct\tmate_correct")
				.append("\tpair_both_correct\tpair_either_correct\tlegacy_mapq")
				.append("\ttruth_reference\ttruth_strand\ttruth_start\ttruth_stop")
				.append("\tmapped_reference\tmapped_strand\tmapped_start\tmapped_stop")
				.append("\tmate_truth_reference\tmate_truth_strand\tmate_truth_start\tmate_truth_stop")
				.append("\tmate_mapped_reference\tmate_mapped_strand\tmate_mapped_start\tmate_mapped_stop")
				.append("\ttruth_insert\tobserved_insert\taverage_pair_distance");
		for(final String name:NeuralMapqPairedFeatureSchema.names()){bb.tab().append(name);}
		return bb.nl();
	}

	public static void append(final Read anchor, final Read mate,
			final int averagePairDistance, final boolean requireCorrectStrands,
			final boolean sameStrandPairs, final float[] vector,
			final NeuralMapqPairedFeatureExtractor.Scratch scratch, final ByteBuilder bb){
		if(anchor==null || mate==null || bb==null){
			throw new IllegalArgumentException("Paired feature row requires anchor, mate, and output");
		}
		final CustomHeader truthA=new CustomHeader(anchor.id,anchor.pairnum());
		final CustomHeader truthM=new CustomHeader(mate.id,mate.pairnum());
		if(truthA.rname==null || truthM.rname==null){
			throw new IllegalArgumentException("mapqpairfeatures requires parseable paired synthetic headers: "+anchor.id);
		}
		final Location mappedA=location(anchor);
		final Location mappedM=location(mate);
		final boolean correctA=correct(mappedA,truthA);
		final boolean correctM=correct(mappedM,truthM);
		NeuralMapqPairedFeatureExtractor.fill(anchor,mate,averagePairDistance,
				requireCorrectStrands,sameStrandPairs,vector,scratch);
		final int observedInsert=NeuralMapqPairedFeatureExtractor.observedInsert(anchor,mate,
				requireCorrectStrands,sameStrandPairs);

		bb.append(anchor.numericID).tab().append(anchor.pairnum()).tab().append(anchor.id)
				.tab().append(correctA ? 0 : 1).tab().append(correctA ? 1 : 0)
				.tab().append(correctM ? 1 : 0).tab().append(correctA && correctM ? 1 : 0)
				.tab().append(correctA || correctM ? 1 : 0).tab().append(SamLine.toMapq(anchor,null));
		appendTruth(bb,truthA);appendLocation(bb,mappedA);
		appendTruth(bb,truthM);appendLocation(bb,mappedM);
		bb.tab().append(truthA.insert).tab().append(observedInsert).tab().append(averagePairDistance);
		for(final float value:vector){bb.tab().append(value,6,true);}
		bb.nl();
	}

	private static boolean correct(final Location mapped, final CustomHeader truth){
		return NeuralMapqOracle.isCorrectLoose(mapped.mapped,mapped.reference,mapped.strand,
				mapped.start,mapped.stop,truth.rname,(byte)truth.strand,truth.start,truth.stop,
				NeuralMapqPairedFeatureSchema.ORACLE_TOLERANCE);
	}

	private static Location location(final Read read){
		if(read==null || !read.mapped()){return Location.UNMAPPED;}
		if(read.chrom<1 || read.start<0 || read.stop<read.start ||
				!Data.isSingleScaffold(read.chrom,read.start,read.stop)){
			throw new IllegalArgumentException("Paired neural MAPQ site does not identify one scaffold: "+read.id);
		}
		final int scaffold=Data.scaffoldIndex(read.chrom,(read.start+read.stop)/2);
		final byte[] name=Data.scaffoldNames[read.chrom][scaffold];
		if(name==null){throw new IllegalArgumentException("Mapped scaffold name unavailable: "+read.id);}
		final int start=Data.scaffoldRelativeLoc(read.chrom,read.start,scaffold);
		return new Location(true,new String(name,StandardCharsets.UTF_8),read.strand(),
				start,start+read.stop-read.start);
	}

	private static void appendTruth(final ByteBuilder bb, final CustomHeader truth){
		bb.tab().append(truth.rname).tab().append(truth.strand).tab().append(truth.start).tab().append(truth.stop);
	}

	private static void appendLocation(final ByteBuilder bb, final Location mapped){
		bb.tab().append(mapped.reference).tab().append((int)mapped.strand)
				.tab().append(mapped.start).tab().append(mapped.stop);
	}

	private static final class Location {
		Location(final boolean mapped_,final String reference_,final byte strand_,final int start_,final int stop_){
			mapped=mapped_;reference=reference_;strand=strand_;start=start_;stop=stop_;
		}
		final boolean mapped;final String reference;final byte strand;final int start,stop;
		static final Location UNMAPPED=new Location(false,".",(byte)-1,-1,-1);
	}
}
