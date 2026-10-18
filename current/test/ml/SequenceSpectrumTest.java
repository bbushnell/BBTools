package test.ml;

import ml.CellNet;
import ml.NetFilter;
import ml.ScoreSequence;
import ml.SequenceToVector;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

/** Literal k=2 spectrum checks, separate from short-span entropy policy.
 * @author Brian Bushnell, Shinobu */
public final class SequenceSpectrumTest{

	/** Runs one named finite case without training or file I/O. */
	public static void main(String[] args){
		if(args.length!=1){throw new IllegalArgumentException("Expected one case name");}
		final String name=args[0];
		final float[] vec=new float[14];
		if(name.equals("spectrum_length")){
			SequenceToVector.fillVector(BASES, vec, 2);
			checkLength(vec);
		}else if(name.equals("spectrum_bins")){
			SequenceToVector.fillVector(BASES, vec, 2);
			checkBins(vec, false);
		}else if(name.equals("spectrum_span")){
			final byte[] padded=ascii("GGAAAAAAAAAACCCCCCCCCCTT");
			SequenceToVector.fillVector(padded, vec, 2, 2, 21);
			checkVector(vec);
		}else if(name.equals("score_sequence")){
			checkScore(ScoreSequence.score(BASES, vec, 2, new RecordingNet(), false));
		}else if(name.equals("netfilter_single")){
			checkScore(NetFilter.scoreSingle(BASES, vec, 2, new RecordingNet()));
		}else if(name.equals("direct_spectrum")){
			Arrays.fill(vec, 0, 4, 7);
			final int count=SequenceToVector.fillSpectrum(BASES, vec, 4, 0, BASES.length, 2);
			if(count!=19){throw new AssertionError("Expected 19 valid dimers, observed "+count);}
			checkPrefix(vec);
			checkBins(vec, false);
		}else if(name.equals("raw_control")){
			final float[] raw=new float[84];
			SequenceToVector.fillVector(BASES, raw, 0);
			close(0.8f, raw[0], "Raw length");
			close(0.5f, raw[1], "Raw GC");
			for(int i=0; i<20; i++){
				for(int j=0; j<4; j++){
					close(j==(i<10 ? 0 : 1) ? 1 : 0, raw[4+4*i+j], "Raw base "+i+" lane "+j);
				}
			}
		}else if(name.equals("empty_spectrum")){
			checkEmpty(ascii(""), vec);
		}else if(name.equals("short_spectrum")){
			checkEmpty(ascii("A"), vec);
		}else if(name.equals("undefined_spectrum")){
			checkEmpty(ascii("NNNNNNNNNNNNNNNNNNNN"), vec);
		}else if(name.equals("broken_spectrum")){
			checkEmpty(ascii("ANANANANANANANANANAN"), vec);
		}else{
			throw new IllegalArgumentException("Unknown case: "+name);
		}
		System.err.println("SequenceSpectrumTest PASS "+name);
	}

	/** Checks only spectrum bins for empty counts; metadata policy is not exercised. */
	private static void checkEmpty(final byte[] bases, final float[] vec){
		Arrays.fill(vec, 0, 4, 7);
		final int count=SequenceToVector.fillSpectrum(bases, vec, 4, 0, bases.length, 2);
		if(count!=0){throw new AssertionError("Expected zero valid dimers, observed "+count);}
		checkPrefix(vec);
		checkBins(vec, true);
	}

	/** The first four entries are not part of fillSpectrum's destination. */
	private static void checkPrefix(final float[] vec){
		for(int i=0; i<4; i++){close(7, vec[i], "Untouched prefix "+i);}
	}

	/** A10C10 has 19 dimers: nine AA, one AC, nine CC. */
	private static void checkBins(final float[] vec, final boolean empty){
		for(int i=0; i<10; i++){
			final int count=(i==0 || i==5 ? 9 : i==1 ? 1 : 0);
			final float expected=(empty ? 0 : (float)(count*2.5/19));
			close(expected, vec[4+i], "Spectrum bin "+i);
		}
	}

	/** Generator spectrum length is valid-count*0.25/10 for k=2. */
	private static void checkLength(final float[] vec){close(0.475f, vec[0], "Spectrum length");}

	/** Checks the finite known spectrum and GC after the public scorer's conversion. */
	private static void checkVector(final float[] vec){
		checkLength(vec);
		close(0.5f, vec[1], "Spectrum GC");
		checkBins(vec, false);
	}

	private static void checkScore(final float score){close(0.25f, score, "Recorded score");}

	/** NaN or infinity must never make a literal fixture comparison pass. */
	private static void close(final float expected, final float actual, final String label){
		if(!(Math.abs(expected-actual)<=1e-6f)){
			throw new AssertionError(label+" expected "+expected+" but observed "+actual);
		}
	}

	private static byte[] ascii(final String s){return s.getBytes(StandardCharsets.US_ASCII);}

	/** Captures the vector supplied by actual scorer entry points, without a trained model. */
	private static final class RecordingNet extends CellNet{
		@Override
		public CellNet applyInput(float[] values){checkVector(values); inputs++; return this;}
		@Override
		public float feedForward(){
			if(inputs!=1){throw new AssertionError("Expected one encoded vector, observed "+inputs);}
			return 0.25f;
		}
		private int inputs;
	}

	private static final byte[] BASES=ascii("AAAAAAAAAACCCCCCCCCC");
	private SequenceSpectrumTest(){}
}
