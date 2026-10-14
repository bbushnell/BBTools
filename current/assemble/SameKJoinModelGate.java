package assemble;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;
import structures.ByteBuilder;

/** Streams all candidate vectors into a truth-blind NN decision. @author Fischl */
final class SameKJoinModelGate {

	/** Loads one model and validates feature metadata before starting the score trace. */
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
		}catch(RuntimeException failure){input.close(); throw failure;}
		output=new ByteStreamWriter(outscore, false, false, false);
		output.start();
		output.print("event\tk\tscore_bits\tscore\taccepted\n");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Scores every candidate, including those with no resolved reference label. */
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

	/** Validates unused control rows and the exact trailer; controls never create graph edges. */
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

	/** Closes both native streams; the merger refuses to emit FASTA if either fails. */
	boolean close(){return input.close() | output.poisonAndWait();}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

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

	private final int k;
	private final float cutoff, encodedK, mean, weighted;
	private final CellNet net;
	private final ByteFile input;
	private final ByteStreamWriter output;
	private final float[] vector=new float[27];
	private final LineParser1 fields=new LineParser1('\t');
	private final ByteBuilder buffer=new ByteBuilder();
}
