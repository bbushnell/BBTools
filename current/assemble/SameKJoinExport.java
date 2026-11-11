package assemble;

import java.util.ArrayList;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import parse.Parse;
import structures.ByteBuilder;

/**
 * Optional diagnostic feature and virtual-control output. Truth never enters
 * this writer. Control identities are a separate namespace and file schema.
 * One serial producer owns the reused row/vector scratch and closes the exporter.
 * Borrowed population moments must remain stable so header metadata and later
 * feature rows describe the same population. Output files are independent streams.
 * @author Fischl
 */
public final class SameKJoinExport implements AutoCloseable{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Paths must already be validated against every input and other output.
	 * Opens/starts the feature writer and optional control-context writer, then emits
	 * headers immediately. A failed control-writer constructor closes the unstarted
	 * feature writer. Startup/header failures attempt closure of both constructed
	 * writers, preserving the original exception with cleanup failures suppressed.
	 * Backend failures inside a writer's own constructor remain its responsibility. */
	SameKJoinExport(final String features, final String controls, final int k_, final int context_,
			final SameKJoinFeatures.Moments moments_){
		if(features==null || k_<2 || context_<1 || moments_==null){throw new IllegalArgumentException("Feature export requires explicit schema and moments.");}
		k=k_; context=context_; moments=moments_;
		featureWriter=new ByteStreamWriter(features, false, false, false);
		try{
			controlWriter=controls==null ? null : new ByteStreamWriter(controls, false, false, false);
		}catch(RuntimeException|Error failure){
			//Both workers are still NEW; ordinary poisonAndWait would wait for startup.
			try{
				if(featureWriter.closeBeforeStart()){
					failure.addSuppressed(new IllegalStateException("Feature output cleanup failed: "+features));
				}
			}catch(RuntimeException|Error cleanup){
				if(cleanup!=failure){failure.addSuppressed(cleanup);}
			}
			throw failure;
		}
		startAndWriteHeaders(featureWriter, controlWriter, moments, row);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Uses original float coverage, not the rounded six-place raw diagnostic text.
	 * Consumes borrowed oriented neighborhoods synchronously; does not retain them.
	 * The caller supplies the consecutive one-based candidate event identity. */
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

	/** Samples and exports contexts without creating Contig objects or assuming positive truth.
	 * Invoke once after candidate rows when producing input for SameKJoinModelGate:
	 * each invocation restarts control IDs at one and appends a terminal metadata row.
	 * This method does not close the streams; close remains the producer's responsibility. */
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

	/** Signals completion and waits for both writers, reporting accumulated I/O error flags. */
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
	 * Starts the constructor-owned, newly opened writers and emits their original
	 * headers using the exporter's scratch row. The control writer may be null.
	 * On unchecked failure, closes NEW writers directly and poisons/joins started
	 * workers, then rethrows the original failure. Partial files are not removed.
	 * Exposed within the package for deterministic lifecycle injection.
	 */
	static void startAndWriteHeaders(final ByteStreamWriter featureWriter, final ByteStreamWriter controlWriter,
			final SameKJoinFeatures.Moments moments, final ByteBuilder row){
		try{
			featureWriter.start();
			if(controlWriter!=null){controlWriter.start(); controlWriter.print(CONTROL_HEADER+"\n");}
			row.append("#schema=").append(SameKJoinFeatures.VERSION).append(" unique=").append(moments.unique);
			row.append(" observations=").append(Double.toString(moments.observations));
			row.append(" squares=").append(Double.toString(moments.squares));
			featureWriter.print(row.nl());
			row.clear().append("kind\tevent");
			for(int i=0; i<SameKJoinFeatures.INPUTS; i++){row.tab().append('f').append(i);}
			featureWriter.print(row.nl());
		}catch(RuntimeException|Error failure){
			try{closeAfterStartupFailure(featureWriter, failure);}
			finally{closeAfterStartupFailure(controlWriter, failure);}
			throw failure;
		}
	}

	/** Attempts one owned writer's cleanup, retaining its error on the original failure.
	 * Actual state matters: start() may throw after launching its worker. NEW writers
	 * must not receive poisonAndWait, which waits for a startup that never occurred. */
	private static void closeAfterStartupFailure(final ByteStreamWriter writer, final Throwable failure){
		if(writer==null){return;}
		try{
			final boolean failed=writer.getState()==Thread.State.NEW ? writer.closeBeforeStart() : writer.poisonAndWait();
			if(failed){failure.addSuppressed(new IllegalStateException("Exporter startup cleanup failed: "+writer.fname));}
		}catch(RuntimeException|Error cleanup){
			if(cleanup!=failure){failure.addSuppressed(cleanup);}
		}
	}

	/**
	 * Nine decimals preserve reachable normalized counts/lengths without exponent
	 * notation or temporary field Strings. Check the actual roundtrip, never silently
	 * round a future out-of-domain feature. Positive native counts/means are >=1.
	 * Appends one tab before each of exactly27 values to caller-owned output; leaves
	 * values unchanged. Rejection can leave a partially appended row in out.
	 */
	public static void appendVector(final ByteBuilder out, final float[] values){
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

	/** Owned feature stream and optional control-context stream; the latter may be null. */
	private final ByteStreamWriter featureWriter, controlWriter;
	/** Word length and extra retained context per side, in bases. */
	private final int k, context;
	/** Borrowed stable table population, used for both metadata and normalized features. */
	private final SameKJoinFeatures.Moments moments;
	/** Serial producer's reusable textual output scratch. */
	private final ByteBuilder row=new ByteBuilder();
	/** Reused normalized candidate/control feature vector. */
	private final float[] vector=new float[SameKJoinFeatures.INPUTS];
	/** Fixed control-context columns; candidate/control feature rows use a separate schema. */
	public static final String CONTROL_HEADER="event\tk\tunitig\tcut\tsource_len\tdest_len\tsource_word\tdest_word\tleft_context\tright_context";
}
