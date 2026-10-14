package assemble;

import java.util.ArrayList;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import parse.Parse;
import structures.ByteBuilder;

/**
 * Optional diagnostic feature and virtual-control output. Truth never enters
 * this writer. Control identities are a separate namespace and file schema.
 * @author Fischl
 */
final class SameKJoinExport implements AutoCloseable {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Paths must already be validated against every input and other output. */
	SameKJoinExport(final String features, final String controls, final int k_, final int context_,
			final SameKJoinFeatures.Moments moments_){
		if(features==null || k_<2 || context_<1 || moments_==null){throw new IllegalArgumentException("Feature export requires explicit schema and moments.");}
		k=k_; context=context_; moments=moments_;
		featureWriter=new ByteStreamWriter(features, false, false, false);
		controlWriter=controls==null ? null : new ByteStreamWriter(controls, false, false, false);
		featureWriter.start();
		if(controlWriter!=null){controlWriter.start(); controlWriter.print(CONTROL_HEADER+"\n");}
		row.append("#schema=").append(SameKJoinFeatures.VERSION).append(" unique=").append(moments.unique);
		row.append(" observations=").append(Double.toString(moments.observations));
		row.append(" squares=").append(Double.toString(moments.squares));
		featureWriter.print(row.nl());
		row.clear().append("kind\tevent");
		for(int i=0; i<SameKJoinFeatures.INPUTS; i++){row.tab().append('f').append(i);}
		featureWriter.print(row.nl());
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Uses original float coverage, not the rounded six-place raw diagnostic text. */
	void candidate(final long event, final Contig source, final boolean sr, final Contig dest,
			final boolean dr, final int[] sourceCounts, final int[] destCounts){
		if(source.length()<k+2 || dest.length()<k+2){
			throw new IllegalArgumentException("Feature export requires formed K+2 or longer unitigs; retain minextension=2.");
		}
		final int sourceInward=AminoAcid.baseToNumber[SameKJoinCollector.baseAt(source, sr, source.length()-k-1)];
		final int destInward=AminoAcid.baseToNumber[SameKJoinCollector.baseAt(dest, dr, k)];
		SameKJoinFeatures.fill(vector, k, sourceCounts, destCounts, sourceInward, destInward,
				source.coverage, dest.coverage, source.length()-k+1, dest.length()-k+1,
				moments.uniqueMean(), moments.observationMean());
		writeVector("candidate", event);
	}

	/** Samples and exports contexts without creating Contig objects or assuming positive truth. */
	SameKAcceptedControls.Result controls(final ArrayList<Contig> contigs, final TadpoleGraph.Counts counts,
			final int requested, final long seed){
		if(controlWriter==null){throw new IllegalStateException("Control output was not requested.");}
		final SameKAcceptedControls sampler=new SameKAcceptedControls(contigs, counts, k);
		final long[] event={0};
		final SameKAcceptedControls.Result result=sampler.collect(requested, seed, new SameKAcceptedControls.Sink(){
			@Override
			public void accept(final int id, final int cut, final Contig c, final double sourceMean,
					final double destMean, final int[] sc, final int[] dc, final int si, final int di){
				final long number=++event[0];
				SameKJoinFeatures.fill(vector, k, sc, dc, si, di, sourceMean, destMean,
						cut-k+1, c.length()-cut, moments.uniqueMean(), moments.observationMean());
				writeVector("control", number);
				row.clear().append(number).tab().append(k).tab().append(id).tab().append(cut);
				row.tab().append(cut).tab().append(c.length()-cut+k-1);
				row.tab().append(SameKJoinCollector.context(c, false, cut-k, cut));
				row.tab().append(SameKJoinCollector.context(c, false, cut-k+1, cut+1));
				row.tab().append(SameKJoinCollector.context(c, false, (int)Math.max(0L, (long)cut-k+1-context), cut));
				row.tab().append(SameKJoinCollector.context(c, false, cut-k+1, (int)Math.min(c.length(), (long)cut+context)));
				controlWriter.print(row.nl());
			}
		});
		assert(event[0]==result.sampled) : "Every sampled cut needs exactly one feature and context row.";
		row.clear().append("#controls population=").append(result.population).append(" sampled=").append(result.sampled);
		row.append(" seed=").append(seed).append(" requested=").append(requested);
		featureWriter.print(row.nl());
		return result;
	}

	/** Both populations use the same normalized floats and exact decimal roundtrip contract. */
	private void writeVector(final String kind, final long event){
		assert(event>0) : "Feature IDs join to consecutive population-specific context records.";
		row.clear().append(kind).tab().append(event);
		appendVector(row, vector);
		featureWriter.print(row.nl());
	}

	/** Flushes every started writer, including when a producer or another writer failed. */
	@Override
	public void close(){
		boolean failed=featureWriter.poisonAndWait();
		if(controlWriter!=null){failed|=controlWriter.poisonAndWait();}
		if(failed){throw new IllegalStateException("Same-K feature/control output failed.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Nine decimals preserve reachable normalized counts/lengths without exponent
	 * notation or temporary field Strings. Check the actual roundtrip, never silently
	 * round a future out-of-domain feature. Positive native counts/means are >=1.
	 */
	static void appendVector(final ByteBuilder out, final float[] values){
		if(values.length!=SameKJoinFeatures.INPUTS){throw new IllegalArgumentException("Wrong feature vector width.");}
		for(float value : values){
			if(!Float.isFinite(value) || value<0){throw new IllegalArgumentException("Invalid normalized feature.");}
			out.tab();
			final int start=out.length();
			out.append(value, 9);
			if(Parse.parseFloat(out.array, start, out.length())!=value){
				throw new IllegalArgumentException("Feature cannot roundtrip with schema precision: "+value);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final ByteStreamWriter featureWriter, controlWriter;
	private final int k, context;
	private final SameKJoinFeatures.Moments moments;
	private final ByteBuilder row=new ByteBuilder();
	private final float[] vector=new float[SameKJoinFeatures.INPUTS];
	static final String CONTROL_HEADER="event\tk\tunitig\tcut\tsource_len\tdest_len\tsource_word\tdest_word\tleft_context\tright_context";
}
