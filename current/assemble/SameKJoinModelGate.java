package assemble;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;
import structures.ByteBuilder;

/** Streams all candidate vectors into a truth-blind NN decision.
 * A serial caller scores candidates in source order, validates the remaining
 * controls with finish, then closes the native streams. Model/parser/vector/output
 * scratch belongs to this instance; it is not shared between concurrent callers.
 * Score rows are emitted during candidate validation, before final trailer checks.
 * The gate neither constructs graph edges nor establishes reference truth.
 * @author Fischl
 */
final class SameKJoinModelGate{

	/** Loads one 27-input/one-output model and validates metadata before starting the
	 * score trace. Paths must already be checked by the caller. Encoded K comes from k_;
	 * metadata supplies unique, observations and squares for the population means.
	 * Every candidate/control row must match that encoded K and those means.
	 * Failed setup closes acquired streams, retaining the original unchecked failure
	 * with cleanup failures suppressed. Partial output files are not removed. */
	SameKJoinModelGate(final String features, final String model, final String outscore, final int k_, final float cutoff_){
		if(k_<24 || k_>127 || !Float.isFinite(cutoff_) || cutoff_<0 || cutoff_>1){
			throw new IllegalArgumentException("This diagnostic model domain is K24–127; supply cutoff in[0,1].");
		}
		k=k_; cutoff=cutoff_;
		net=CellNetParser.load(model, false);
		if(net==null || net.numInputs()!=27 || net.numOutputs()!=1){throw new IllegalArgumentException("Expected a27-input binary model.");}
		net.simdFF=false;// Matches validation's layer setting, not a promise of scalar dot products.
		// Vector.fma independently follows simdfma; set it false for exact scalar trace parity.
		input=ByteFile.makeByteFile(features, true);
		ByteStreamWriter acquired=null;
		try{
			final byte[] first=input.nextLine();
			if(first==null){throw new IllegalArgumentException("Missing feature metadata.");}
			final String[] meta=new String(first, StandardCharsets.UTF_8).split(" ", -1);
			if(meta.length!=4 || !meta[0].equals("#schema="+SameKJoinFeatures.VERSION)){
				throw new IllegalArgumentException("Wrong same-K feature schema.");
			}
			final long unique=number(meta[1], "unique=");
			final double observations=decimal(meta[2], "observations="), squares=decimal(meta[3], "squares=");
			if(unique<0 || observations<unique || squares<observations || (unique==0 && (observations!=0 || squares!=0))){
				throw new IllegalArgumentException("Invalid exact depth moments.");
			}
			mean=SameKJoinFeatures.depth(unique==0 ? 0 : observations/unique);
			weighted=SameKJoinFeatures.depth(observations==0 ? 0 : squares/observations);
			encodedK=(float)(0.125*Math.log(k)*(1/Math.log(2)));
			final ByteBuilder header=new ByteBuilder("kind\tevent");
			for(int i=0; i<27; i++){header.tab().append('f').append(i);}
			if(!Arrays.equals(input.nextLine(), header.toBytes())){throw new IllegalArgumentException("Wrong feature columns.");}
			acquired=new ByteStreamWriter(outscore, false, false, false);
			acquired.start();
			acquired.print("event\tk\tscore_bits\tscore\taccepted\n");
		}catch(RuntimeException|Error failure){
			closeAfterConstructionFailure(input, acquired, failure);
			throw failure;
		}
		output=acquired;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Consumes and scores the next candidate, including those without resolved truth.
	 * Consumes the input row and queues its score-trace row before returning.
	 * @param event One-based event identity already checked by the merger
	 * @return Whether the finite probability meets the inclusive cutoff
	 */
	boolean accept(final int event){
		assert(event>0) : "Candidate event IDs are one-based and checked in source order by the merger.";
		readVector(input.nextLine(), "candidate", event);
		final float score=net.applyInput(vector).feedForward();
		if(!Float.isFinite(score) || score<0 || score>1){throw new IllegalStateException("Model output is not a finite probability: "+score);}
		final boolean accepted=score>=cutoff;
		output.print(buffer.clear().append(event).tab().append(k).tab().append(Float.floatToRawIntBits(score))
				.tab().appendSlow(score).tab().append(accepted ? 1 : 0).nl());
		return accepted;
	}

	/** Validates all remaining control rows and the single final trailer after the
	 * caller has scored every candidate. Controls are not evaluated by the model.
	 * Does not close either stream; the caller still invokes close on failure. */
	void finish(){
		int controls=0;
		boolean trailer=false;
		for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
			if(trailer){throw new IllegalArgumentException("Feature data follows the control trailer.");}
			if(line.length>0 && line[0]=='#'){
				final String[] terms=new String(line, StandardCharsets.UTF_8).split(" ", -1);
				if(terms.length!=5 || !terms[0].equals("#controls")){throw new IllegalArgumentException("Invalid control trailer.");}
				final long population=number(terms[1], "population="), sampled=number(terms[2], "sampled=");
				number(terms[3], "seed=");// Provenance only; it does not influence NN inference.
				final long requested=number(terms[4], "requested=");
				if(population<0 || requested<0 || sampled!=controls || sampled!=Math.min(population, requested)){
					throw new IllegalArgumentException("Control count disagrees with feature trailer.");
				}
				trailer=true;
			}else{readVector(line, "control", ++controls);}
		}
		if(!trailer){throw new IllegalArgumentException("Missing control trailer.");}
	}

	/** Reads exactly27 normalized values, checking population identity, K and global provenance. */
	private void readVector(final byte[] line, final String kind, final int event){
		if(line==null){throw new IllegalArgumentException("Missing "+kind+" feature row "+event);}
		fields.set(line);
		if(fields.terms()!=29 || !fields.termEquals(kind, 0) || fields.parseInt(1)!=event){
			throw new IllegalArgumentException("Feature population/ID/width mismatch at "+kind+" "+event);
		}
		for(int i=0; i<27; i++){
			vector[i]=fields.parseFloat(i+2);
			if(!Float.isFinite(vector[i]) || vector[i]<0){throw new IllegalArgumentException("Invalid normalized feature at "+event);}
		}
		if(vector[24]!=encodedK || vector[25]!=mean || vector[26]!=weighted){
			throw new IllegalArgumentException("Feature K/global moments disagree with source metadata.");
		}
	}

	/** Closes input and joins the score writer; does not consume or validate remaining rows.
	 * @return True if either native close reports an I/O error, causing the merger to refuse FASTA
	 */
	boolean close(){return input.close() | output.poisonAndWait();}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Attempts both acquired resources, preserving the primary failure and its cause.
	 * A writer can start before start() throws; use actual state, never return status.
	 * NEW writers need direct close because poisonAndWait waits for worker startup.
	 * Resources inside a writer constructor that never returns remain its responsibility. */
	private static void closeAfterConstructionFailure(final ByteFile input, final ByteStreamWriter writer, final Throwable failure){
		try{
			if(input.close()){failure.addSuppressed(new IllegalStateException("Model-gate feature input cleanup failed."));}
		}catch(RuntimeException|Error cleanup){
			if(cleanup!=failure){failure.addSuppressed(cleanup);}
		}finally{
			if(writer!=null){
				try{
					final boolean failed=writer.getState()==Thread.State.NEW ? writer.closeBeforeStart() : writer.poisonAndWait();
					if(failed){failure.addSuppressed(new IllegalStateException("Model-gate score output cleanup failed: "+writer.fname));}
				}catch(RuntimeException|Error cleanup){
					if(cleanup!=failure){failure.addSuppressed(cleanup);}
				}
			}
		}
	}

	/** Checks a named long field from the producer's cold metadata. */
	private static long number(final String text, final String prefix){
		if(!text.startsWith(prefix)){throw new IllegalArgumentException("Expected feature metadata "+prefix);}
		return Long.parseLong(text.substring(prefix.length()));
	}

	/** Exact table moments may be printed in exponent notation. */
	private static double decimal(final String text, final String prefix){
		if(!text.startsWith(prefix)){throw new IllegalArgumentException("Expected feature metadata "+prefix);}
		final double value=Double.parseDouble(text.substring(prefix.length()));
		if(!Double.isFinite(value) || value<0){throw new IllegalArgumentException("Invalid feature moment.");}
		return value;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Word length expected in this feature stream, in bases. */
	private final int k;
	/** Acceptance cutoff, K encoded from k_, and metadata-derived unique/observation-weighted means. */
	private final float cutoff, encodedK, mean, weighted;
	/** Owned inference network; layer SIMD is disabled but shared FMA remains independent. */
	private final CellNet net;
	/** Owned feature stream, positioned after its metadata and column header at construction. */
	private final ByteFile input;
	/** Owned asynchronous score-trace writer, closed after finish or a processing failure. */
	private final ByteStreamWriter output;
	/** Reused normalized vector for the current candidate or unscored control row. */
	private final float[] vector=new float[27];
	/** Reused parser over the current borrowed input line. */
	private final LineParser1 fields=new LineParser1('\t');
	/** Serial score-row formatting scratch. */
	private final ByteBuilder buffer=new ByteBuilder();
}
