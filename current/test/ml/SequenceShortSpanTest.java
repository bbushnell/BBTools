package test.ml;

import ml.CellNet;
import ml.NetFilter;
import ml.ScoreSequence;
import ml.SequenceToVector;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import structures.ByteBuilder;

/** Checks short nonempty vectors against literal base lanes and a direct entropy formula.
 * @author Shinobu */
public final class SequenceShortSpanTest{

	/** Runs one independent conversion or scorer case. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Expected one case name");}
		final String name=args[0];
		final int len=name.equals("raw1") ? 1 : name.equals("raw2") ? 2 :
			name.equals("raw15") ? 15 : name.endsWith("16") ? 16 : 7;
		final byte[] bases=new byte[len];
		for(int i=0; i<len; i++){bases[i]=(byte)"ACGT".charAt(i&3);}
		final float[] vec=new float[84];
		Arrays.fill(vec, 9);
		try{
			if(name.startsWith("raw") || name.equals("reuse7")){
				if(!(name.equals("raw1") || name.equals("raw2") || name.equals("raw7") ||
						name.equals("raw15") || name.equals("raw16") || name.equals("reuse7"))){
					throw new IllegalArgumentException("Unknown case: "+name);
				}
				if(name.equals("reuse7")){
					SequenceToVector.fillVector("AAAAAAAAAAAAAAAA".getBytes(StandardCharsets.US_ASCII), vec, 0);
				}
				SequenceToVector.fillVector(bases, vec, 0);
				checkRaw(bases, vec, 1e-6);
			}else if(name.equals("span7")){
				final byte[] padded=("GG"+new String(bases, StandardCharsets.US_ASCII)+"TT").getBytes(StandardCharsets.US_ASCII);
				SequenceToVector.fillVector(padded, vec, 0, 2, 8);
				checkRaw(bases, vec, 1e-6);
			}else if(name.equals("score7") || name.equals("netfilter7")){
				final RecordingNet net=new RecordingNet(bases);
				final float score=name.equals("score7") ? ScoreSequence.score(bases, vec, 0, net, false) :
					NetFilter.scoreSingle(bases, vec, 0, net);
				close(0.25, score, 0, "Scorer output");
				if(net.calls!=1){throw new AssertionError("Expected one inference, observed "+net.calls);}
			}else if(name.equals("spectrum7")){
				final float[] spectrum=new float[14];
				SequenceToVector.fillVector(bases, spectrum, 2);
				close(0.15, spectrum[0], 1e-6, "Six valid dimers encoded as length");
				close(entropy(bases), spectrum[2], 1e-6, "Spectrum entropy");
				// ACGTACG has AC/CG/GT/TA/AC/CG; reverse complements merge AC with GT.
				for(int i=0; i<10; i++){
					final int count=i==1 ? 3 : i==6 ? 2 : i==9 ? 1 : 0;
					close(count*2.5/6, spectrum[4+i], 1e-6, "Spectrum bin "+i);
				}
			}else if(name.equals("text7") || name.equals("text16")){
				final ByteBuilder bb=new ByteBuilder();
				final Method method=SequenceToVector.class.getDeclaredMethod("toVector", byte[].class,
					ByteBuilder.class, int.class, int.class, float.class, float[].class);
				method.setAccessible(true);
				method.invoke(null, bases, bb, 20, 0, 1f, null);
				final String[] fields=bb.toString().split("\t", -1);
				if(fields.length!=85){throw new AssertionError("Expected84features+label, observed "+fields.length);}
				for(int i=0; i<84; i++){vec[i]=Float.parseFloat(fields[i]);}
				checkRaw(bases, vec, 0.000051);
				close(1, Float.parseFloat(fields[84]), 0, "Target label");
			}else{
				throw new IllegalArgumentException("Unknown case: "+name);
			}
		}catch(NullPointerException e){
			throw new AssertionError("Short-span conversion failed: "+name, e);
		}catch(InvocationTargetException e){
			throw new AssertionError("Short-span conversion failed: "+name, e.getCause());
		}
		System.err.println("SequenceShortSpanTest PASS "+name);
	}

	/** Checks all raw lanes, zero padding, and four independently computed features. */
	private static void checkRaw(final byte[] bases, final float[] vec, final double tolerance){
		int gc=0;
		for(byte b : bases){if(b=='C' || b=='G'){gc++;}}
		close(bases.length/25.0, vec[0], tolerance, "Length");
		close(gc/(double)bases.length, vec[1], tolerance, "GC");
		close(entropy(bases), vec[2], tolerance, "Entropy");
		close(1.0/6, vec[3], tolerance, "Homopolymer");
		for(int i=0; i<20; i++){
			for(int j=0; j<4; j++){
				close(i<bases.length && j==(i&3) ? 1 : 0, vec[4+4*i+j], 0, "Base lane "+i+":"+j);
			}
		}
	}

	/** Direct sum of trinucleotide probabilities with the approved 14-kmer denominator. */
	private static double entropy(final byte[] bases){
		final int[] counts=new int[64];
		for(int i=0; i+2<bases.length; i++){
			counts[((i&3)<<4)|(((i+1)&3)<<2)|((i+2)&3)]++;
		}
		double sum=0;
		for(int count : counts){
			if(count>0){final double p=count/14.0; sum-=p*Math.log(p)/Math.log(14);}
		}
		return sum;
	}

	private static void close(final double expected, final double actual, final double tolerance, final String label){
		if(!(Math.abs(expected-actual)<=tolerance)){
			throw new AssertionError(label+" expected "+expected+" but observed "+actual);
		}
	}

	/** Receives the real caller's vector without loading or training a model. */
	private static final class RecordingNet extends CellNet{
		RecordingNet(byte[] bases_){bases=bases_;}
		@Override public CellNet applyInput(float[] vec){checkRaw(bases, vec, 1e-6); return this;}
		@Override public float feedForward(){calls++; return 0.25f;}
		final byte[] bases;
		int calls;
	}

	private SequenceShortSpanTest(){}
}
