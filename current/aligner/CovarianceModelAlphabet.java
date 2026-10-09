package aligner;

import java.util.Arrays;

/** Infernal1.1.5 ambiguity semantics, precomputed outside the DP hot loops.
 * src/cm.c:3603-3695 uses Easel FAvgScore for singlets and src/alphabet.c
 * FastPairScore* for pairs: equal-weight averages of canonical BIT SCORES.
 * Do not replace this with log-summed probabilities or a best-base score.
 * @author Brian Bushnell, Raiden
 */
final class CovarianceModelAlphabet {
	/** Standalone CM tools validate with this alphabet, not Read's DNA alphabet
	 * (which excludes Easel's I=A alias). Never call from a scorer constructor:
	 * embedded callers retain their own input policy. No residues are rewritten. */
	static void configureStandaloneInput(){
		require(!stream.Read.IUPAC_TO_N && !stream.Read.DOT_DASH_X_TO_N && !stream.Read.LOWER_CASE_TO_N,
			"CM comparisons require original ambiguity symbols; input-to-N transformations are incompatible");
		require(stream.Read.JUNK_MODE==stream.Read.CRASH_JUNK || stream.Read.JUNK_MODE==stream.Read.IGNORE_JUNK,
			"CM tools validate every residue themselves; fixing or discarding input would change the scoring population");
		stream.Read.JUNK_MODE=stream.Read.IGNORE_JUNK;
	}
	/** Canonical codes stay0..3; pair tables use a fixed16 stride. */
	static int code(byte b){return CODES[b&255];}
	static byte normalized(byte b){
		final int c=code(b);if(c<0){throw invalid(b, -1);}return LETTERS[c];
	}
	static byte complement(byte b){
		final int c=code(b);if(c<0){throw invalid(b, -1);}return COMPLEMENT[c];
	}
	static byte[] encode(byte[] sequence){
		assert(sequence!=null):"CM encoding must preserve an existing input and report its actual bad position";
		final byte[] out=new byte[sequence.length+1];
		for(int i=0; i<sequence.length; i++){
			final int c=code(sequence[i]);if(c<0){throw invalid(sequence[i], i+1);}out[i+1]=(byte)c;
		}
		return out;
	}
	static float[][] expand(CovarianceModel model){
		assert(model!=null):"Expanded emissions belong to one fully constructed immutable CM";
		final float[][] out=new float[model.states()][];
		for(int v=0; v<out.length; v++){
			final float[] original=model.emissionScore[v];
			final int expected=CovarianceModelCyk.delta(model.type[v]);
			if(expected==0){require(original==null || original.length==0, "Nonemitting states cannot acquire emission scores, state="+v);out[v]=original;continue;}
			require(original!=null && original.length==(expected==2 ? 16 : 4),
				"Canonical emission dimensions must agree with state type before expansion, state="+v);
			final float[] expanded=new float[expected==2 ? 256 : 16];Arrays.fill(expanded, Float.NEGATIVE_INFINITY);
			for(int a=0; a<LETTERS.length; a++){
				if(expected==1){expanded[a]=singlet(original, a);}
				else{for(int b=0; b<LETTERS.length; b++){expanded[(a<<4)|b]=pair(original, a, b);}}
			}
			out[v]=expanded;
		}
		return out;
	}
	static float singlet(float[] scores, int code){
		assert(scores.length==4 && code>=0 && code<MASKS.length):"Singlet expansion uses canonical A,C,G,T scores and a valid residue mask";
		if(code<4){return scores[code];}
		float sum=0;final int mask=MASKS[code];
		for(int a=0; a<4; a++){if((mask&(1<<a))!=0){sum+=scores[a];}}
		return sum/Integer.bitCount(mask);
	}
	static float pair(float[] scores, int left, int right){
		assert(scores.length==16 && left>=0 && right>=0 && left<MASKS.length && right<MASKS.length):
			"Pair expansion uses canonical4x4 scores and valid left/right residue masks";
		if(left<4 && right<4){return scores[(left<<2)|right];}
		final int lm=MASKS[left], rm=MASKS[right];
		final float lw=1f/Integer.bitCount(lm), rw=1f/Integer.bitCount(rm);float sum=0;
		// Exclude zero-weight terms: Java uses -Infinity for impossible emissions,
		// whereas Infernal's finite IMPOSSIBLE sentinel permits multiplication by0.
		// The ascending loop and float multiplication order match FastPairScore*.
		for(int a=0; a<4; a++){if((lm&(1<<a))==0){continue;}
			for(int b=0; b<4; b++){if((rm&(1<<b))==0){continue;}
				float term=scores[(a<<2)|b];if(left>=4){term*=lw;}if(right>=4){term*=rw;}sum+=term;
			}
		}
		return sum;
	}
	private static CovarianceModelCyk.UnsupportedSequenceException invalid(byte b, int position){
		return new CovarianceModelCyk.UnsupportedSequenceException("Unsupported nucleotide"+(position>0 ? " at position "+position : "")+": "+(char)(b&255));
	}
	private static byte[] codes(){
		final byte[] table=new byte[256];Arrays.fill(table, (byte)-1);
		for(int i=0; i<LETTERS.length; i++){table[LETTERS[i]]=(byte)i;table[LETTERS[i]+32]=(byte)i;}
		// Bundled Easel esl_alphabet.c:create_rna, lines188-190.
		table['U']=table['u']=3;table['X']=table['x']=14;table['I']=table['i']=0;
		return table;
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final byte[] LETTERS={'A','C','G','T','R','Y','M','K','S','W','H','B','V','D','N'};
	private static final byte[] COMPLEMENT={'T','G','C','A','Y','R','K','M','S','W','D','V','B','H','N'};
	private static final int[] MASKS={1,2,4,8,5,10,3,12,6,9,11,14,7,13,15};
	private static final byte[] CODES=codes();
	private CovarianceModelAlphabet(){}
}
