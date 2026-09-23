package align2;

import java.nio.charset.StandardCharsets;

import dna.Data;
import stream.CustomHeader;
import stream.Read;
import stream.SamLine;
import structures.ByteBuilder;

/** Serializes one synthetic mapped read and its single-end pilot features. */
public final class NeuralMapqFeatureRow {

	private NeuralMapqFeatureRow(){}

	public static ByteBuilder header(){
		final ByteBuilder bb=new ByteBuilder(1024);
		bb.append("#schema\t").append(NeuralMapqFeatureSchema.SCHEMA_NAME).nl();
		bb.append("#oracle_tolerance\t").append(NeuralMapqFeatureSchema.ORACLE_TOLERANCE).nl();
		bb.append("#max_mapq\t").append(NeuralMapqFeatureSchema.MAX_MAPQ).nl();
		bb.append("numeric_id\tread_name\terror\tcorrect_loose\tlegacy_mapq\ttruth_reference\ttruth_strand")
				.append("\ttruth_start\ttruth_stop\tmapped_reference\tmapped_strand")
				.append("\tmapped_start\tmapped_stop");
		for(final String name : NeuralMapqFeatureSchema.names()){bb.tab().append(name);}
		return bb.nl();
	}

	public static void append(final Read read, final float[] vector,
			final NeuralMapqFeatureExtractor.Scratch scratch, final ByteBuilder bb){
		if(read==null || bb==null){throw new IllegalArgumentException("Feature row requires read and output buffer");}
		final CustomHeader truth=new CustomHeader(read.id,0);
		if(truth.rname==null){throw new IllegalArgumentException("mapqfeatures requires a parseable synthetic header: "+read.id);}
		NeuralMapqFeatureExtractor.fill(read,vector,scratch);
		final int mappedScaffold=scaffoldIndex(read.chrom,read.start,read.stop);
		final byte[] mappedReference=Data.scaffoldNames[read.chrom][mappedScaffold];
		if(mappedReference==null){
			throw new IllegalArgumentException("Neural MAPQ mapped scaffold name is unavailable: "+read.id);
		}
		final int mappedStart=Data.scaffoldRelativeLoc(read.chrom,read.start,mappedScaffold);
		final int mappedStop=mappedStart+read.stop-read.start;
		final boolean correct=NeuralMapqOracle.isCorrectLoose(true,
				new String(mappedReference,StandardCharsets.UTF_8),read.strand(),
				mappedStart,mappedStop,truth.rname,(byte)truth.strand,truth.start,truth.stop,
				NeuralMapqFeatureSchema.ORACLE_TOLERANCE);

		bb.append(read.numericID).tab().append(read.id).tab().append(correct ? 0 : 1)
				.tab().append(correct ? 1 : 0).tab().append(SamLine.toMapq(read,null)).tab().append(truth.rname)
				.tab().append(truth.strand).tab().append(truth.start).tab().append(truth.stop)
				.tab().append(mappedReference).tab().append((int)read.strand())
				.tab().append(mappedStart).tab().append(mappedStop);
		for(final float value : vector){bb.tab().append(value,6,true);}
		bb.nl();
	}

	private static int scaffoldIndex(final int chrom, final int start, final int stop){
		if(chrom<1 || start<0 || stop<start || !Data.isSingleScaffold(chrom,start,stop)){
			throw new IllegalArgumentException("Neural MAPQ site does not identify one scaffold: chrom="+
					chrom+", start="+start+", stop="+stop);
		}
		return Data.scaffoldIndex(chrom,(start+stop)/2);
	}
}
