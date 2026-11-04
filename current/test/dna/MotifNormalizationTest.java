package test.dna;

import java.nio.charset.StandardCharsets;

import dna.MotifMulti;
import dna.MotifProbsN;

/** Finite constructor/importance rejection cases and unchanged motif-score controls.
 * Each case runs in a fresh JVM; no matrix resource files or reference genome are needed.
 * @author Shinobu
 */
public final class MotifNormalizationTest{

	/** Runs the selected case, keeping oracle checks active with assertions disabled too. */
	public static void main(String[] args){
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Require case=controls|background|constant|invalid|underflow|overflow|importance");
		}
		final String mode=args[0].substring(5);
		int failures=0;
		if(mode.equals("controls")){controls();}
		else if(mode.equals("importance")){
			final MotifProbsN motif=ordinary();
			boolean rejected=false;
			try{
				final float[] weights=motif.positionImportance(new float[][]{BACKGROUND.clone()});
				System.out.println("Returned background weight="+weights[0]);
			}catch(IllegalArgumentException e){
				rejected=true;
				require(e.getMessage()!=null && e.getMessage().contains("ordinary"), "Diagnostic must identify the motif");
				System.out.println(e.getMessage());
			}
			require(rejected, "FAIL_REJECTION importance: background-only weights must not become NaN");
		}else if(mode.equals("invalid")){
			for(float value : new float[]{Float.NaN, Float.POSITIVE_INFINITY, Float.NEGATIVE_INFINITY, -1f, Float.MAX_VALUE}){
				final float[][] p={BACKGROUND.clone(), BACKGROUND.clone()};
				p[0][0]*=2;
				p[1][1]=value;
				failures+=(reject(p, "invalid-"+value) ? 0 : 1);
			}
		}else{
			final float[][] p;
			if(mode.equals("background")){p=new float[][]{BACKGROUND.clone(), BACKGROUND.clone()};}
			else if(mode.equals("constant")){
				p=new float[][]{BACKGROUND.clone()};
				for(int j=0; j<4; j++){p[0][j]*=2;}
			}else if(mode.equals("underflow") || mode.equals("overflow")){
				p=new float[200][4];
				for(int i=0; i<p.length; i++){
					for(int j=0; j<4; j++){
						p[i][j]=BACKGROUND[j]*(mode.equals("underflow") ? (j+1)*0.1f : (j+2));
					}
				}
			}else{throw new IllegalArgumentException(mode);}
			failures+=(reject(p, mode) ? 0 : 1);
			if(mode.equals("background")){failures+=(reject(new float[1][4], "all-zero") ? 0 : 1);}
			if(mode.equals("underflow")){
				final float[][] tiny={BACKGROUND.clone()};
				for(int j=0; j<4; j++){tiny[0][j]*=(j+1)*1e-22f;}
				failures+=(reject(tiny, "reciprocal-overflow") ? 0 : 1);
			}
		}
		require(failures==0, "FAIL_REJECTION "+mode+": "+failures+" unusable matrices returned");
		System.out.println("PASS "+mode);
	}

	/** Returns whether construction itself rejects; scoring exceptions cannot satisfy this oracle. */
	private static boolean reject(float[][] p, String name){
		final MotifProbsN motif;
		try{
			motif=new MotifProbsN(name, p, 0, 1);
		}catch(IllegalArgumentException e){
			require(e.getMessage()!=null && e.getMessage().contains(name), "Diagnostic must identify motif "+name);
			System.out.println(e.getMessage());
			return true;
		}
		System.out.println("Returned "+name+" min="+motif.minProb+" max="+motif.maxProb+
				" inverse="+motif.invProbDif+" normalizedMin="+motif.normalize(motif.minProb));
		return false;
	}

	/** Checks literal ratios, background-row tolerance, tuple widths, zeros and composite use. */
	private static void controls(){
		final MotifProbsN motif=ordinary();
		require(motif.importance[0]==1 && motif.importance[1]==0,
				"A background-only row is valid alongside an informative row");
		final float[] expected={3.4822023f, 1f, 1f, 0.28717458f};
		for(int i=0; i<4; i++){
			final byte[] sequence={(byte)"ACGT".charAt(i), 'A'};
			final float score=motif.matchStrength(sequence, 0);
			require(Math.abs(score-expected[i])<0.000001f, "Explicit likelihood-ratio control changed at "+i);
			System.out.println("ordinary "+i+" "+Float.toHexString(score)+" "+Float.toHexString(motif.normalize(score)));
		}
		require(motif.normalize(motif.minProb)==0, "Minimum strength must normalize to zero");
		require(Math.abs(motif.normalize(motif.maxProb)-1)<0.000001f, "Maximum strength must normalize to one");
		final byte[] best="AA".getBytes(StandardCharsets.US_ASCII);
		final MotifMulti multi=new MotifMulti("composite", motif);
		require(multi.matchStrength(best, 0)==motif.normalize(motif.matchStrength(best, 0)),
				"Composite must retain the child's normalized score");
		require(motif.matchStrength(new byte[]{'N', 'A'}, 0)==motif.minProb, "Undefined tuple fallback changed");
		require(motif.matchStrength(best, -1)==motif.minProb, "Out-of-range fallback changed");
		for(int n=1; n<=3; n++){
			final int width=1<<(2*n);
			final float[][] p=new float[1][width];
			for(int j=0; j<width; j++){p[0][j]=(j+1f)/(width*(width+1)/2);}
			final MotifProbsN tuple=new MotifProbsN("tuple"+n, p, 0, n);
			for(int code=0; code<width; code++){
				final byte[] sequence=new byte[n];
				for(int j=0; j<n; j++){sequence[j]=(byte)"ACGT".charAt((code>>(2*(n-j-1)))&3);}
				final float score=tuple.matchStrength(sequence, 0), normalized=tuple.normalize(score);
				require(Float.isFinite(score) && Float.isFinite(normalized), "Supported tuple must produce finite scores");
				System.out.println("tuple "+n+" "+code+" "+Float.toHexString(score)+" "+Float.toHexString(normalized));
			}
		}
		final float[][] zero={BACKGROUND.clone()};
		zero[0][0]=0;
		final MotifProbsN withZero=new MotifProbsN("valid-zero", zero, 0, 1);
		require(withZero.minProb==0 && withZero.normalize(0)==0, "An individual zero probability is valid");
		System.out.println("zero "+Float.toHexString(withZero.maxProb)+" "+Float.toHexString(withZero.normalize(0)));
	}

	/** Uses exact background constants with explicit ratios 2,1,1,0.5 in its first row. */
	private static MotifProbsN ordinary(){
		final float[][] p={BACKGROUND.clone(), BACKGROUND.clone()};
		p[0][0]*=2;
		p[0][3]*=0.5f;
		return new MotifProbsN("ordinary", p, 0, 1);
	}

	/** Throws independently of -ea/-da so both configurations check the same contract. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Literal monomer background from Motif; fixture inputs, not values obtained from the tested method. */
	private static final float[] BACKGROUND={0.256614f, 0.226617f, 0.238012f, 0.278756f};
}
