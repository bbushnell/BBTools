package consensus;

import java.util.Arrays;

import dna.AminoAcid;
import shared.Tools;
import structures.ByteBuilder;

/** Mutable observations at one reference, insertion or deletion position.
 * A/C/G/T have separate counts and weights; other symbols contribute only to
 * the totals. DEL nodes have no allele arrays. BaseGraph adds its reference
 * pseudocount; this constructor does not. Updates and cached probabilities are
 * unsynchronized, so callers finish accumulation before scoring or traversal.
 * Integer accumulators do not saturate on overflow.
 * @author Brian Bushnell */
public class BaseNode extends BaseGraphPart implements Comparable<BaseNode>{

	/** Retained model serialization identifier. */
	private static final long serialVersionUID=7932097131372307182L;
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates an empty node, narrowing the reference symbol to a byte.
	 * @param refBase_ Reference symbol
	 * @param type_ REF, INS or DEL
	 * @param rpos_ Reference coordinate */
	public BaseNode(char refBase_, int type_, int rpos_){
		this((byte)refBase_, type_, rpos_);
	}
	
	/** Creates a node with no observations or cached probabilities.
	 * @param refBase_ Reference symbol, or an insertion placeholder
	 * @param type_ REF, INS or DEL; DEL omits allele arrays
	 * @param rpos_ Reference coordinate */
	public BaseNode(byte refBase_, int type_, int rpos_){
		super(type_);
		refBase=refBase_;
		rpos=rpos_;
		acgtWeight=(type==DEL ? null : new int[4]);
		acgtCount=(type==DEL ? null : new int[4]);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Adds one weighted observation; ambiguous symbols still contribute to totals.
	 * The caller may include quality, mapping or identity adjustments in the weight.
	 * Does not invalidate an already constructed probability cache.
	 * @param base Symbol indexing AminoAcid's base map
	 * @param quality Observation weight, not an ASCII-encoded quality character */
	public void add(byte base, int quality){
		int num=AminoAcid.baseToNumber[base];
		if(num>=0 && type!=DEL){
			acgtWeight[num]+=quality;
			acgtCount[num]++;
		}
		countSum++;
		weightSum+=quality;
	}
	
//	public byte consensus(){
//		assert(type!=DEL);
//		int maxPos=AminoAcid.baseToNumber0[refBase];
//		int max=acgtWeight[maxPos];
//		if(max*2>=weight && AminoAcid.isFullyDefined(refBase)){return refBase;}//Common case
//		if(maxPos<0){maxPos=0;}
//		
//		for(int i=0; i<acgtWeight.length; i++){
//			int weight=acgtWeight[i];
//			int depth=acgtCount[i];
//			if(weight>max && depth>=minDepth){
//				max=weight;
//				maxPos=i;
//			}
//		}
//		return AminoAcid.numberToBase[maxPos];
//	}
	
	/** Selects a REF/INS allele and writes its base and numeric confidence score.
	 * A reference with at least half the weight keeps the tie shortcut; otherwise
	 * weight wins, then observation count. Depth/allele-frequency policy can retain
	 * a REF symbol. Scores use weighted allele odds clamped to 2..41; fallback
	 * scores remain 20 for protected references, 2 for low support, and 0 for N.
	 * These are internal numeric Phred-scale scores, not FASTQ characters or a
	 * claim of independently calibrated error probabilities.
	 * @param r Caller-owned scratch with at least two bytes: base, numeric score
	 * @return Selected or retained base; also stored at r[0] */
	public byte consensus(byte[] r){
		return consensusInternal(r, true);
	}

	/** Selects an allele when traversal will calculate its own diversity quality.
	 * Exact BaseNode instances leave r[1] unchanged. Subclasses use their public
	 * consensus implementation so existing overrides and side effects are retained.
	 * @param r Caller-owned scratch with at least two bytes for subclass delegation
	 * @return Selected or retained base; also stored at r[0] */
	byte consensusBase(byte[] r){
		return getClass()==BaseNode.class ? consensusInternal(r, false) : consensus(r);
	}

	/** Shares allele selection between full-confidence and internal base-only calls.
	 * @param r Caller-owned base/quality scratch
	 * @param calculateQuality Write confidence, including fallbacks, when true
	 * @return Selected or retained base */
	private byte consensusInternal(byte[] r, boolean calculateQuality){
		assert(type!=DEL) : "DEL nodes have no allele arrays; consensus requires REF/INS at position "+rpos;
		if(onlyConvertNs && type==REF && refBase!='N'){
			r[0]=refBase;
			if(calculateQuality){r[1]=20;}
			return refBase;
		}
		int maxPos=AminoAcid.baseToNumber0[refBase];
		int maxWeight=acgtWeight[maxPos];
		int maxDepth=acgtCount[maxPos];
//		long sum=shared.Vector.sum(acgtWeight);
		
		if(acgtWeight[maxPos]*2<weightSum){//Uncommon case of ref being a minority
			for(int i=0; i<acgtWeight.length; i++){
				int x=acgtWeight[i];
				int y=acgtCount[i];
				if(x>maxWeight || (x==maxWeight && y>maxDepth)){
					maxWeight=x;
					maxDepth=y;
					maxPos=i;
				}
			}
		}
//		maxWeight=acgtWeight[maxPos];
//		maxDepth=acgtCount[maxPos];
		
		if(type==REF){
			final byte b=AminoAcid.numberToBase[maxPos];
			float af=maxDepth/(float)countSum;
			float maf=(refBase=='N' ? MAF_noref : MAF_sub);
			if(af<maf || maxDepth<minDepth){
				r[0]=refBase;
				if(calculateQuality){r[1]=(refBase=='N' ? (byte)0 : (byte)2);}
			}else{
				r[0]=b;
				if(calculateQuality){
					//CON-004: use selected-versus-competing weighted odds and retain numeric units.
					double quality=Tools.mid(2, 41, 10*Math.log10(maxWeight/Tools.max(0.01, (double)weightSum-maxWeight)));
					r[1]=(byte)quality;
				}
			}
			assert(r[0]!='N' || maxDepth<minDepth || af<maf) : "\n"+Arrays.toString(acgtCount)+", max="+maxDepth+", sum="+countSum+", maf="+maf+", af="+af+"\n"
					+ Arrays.toString(acgtWeight)+", max="+maxWeight+", weightSum="+weightSum+"\n";
		}else{
			assert(type==INS);
			r[0]=AminoAcid.numberToBase[maxPos];
			if(calculateQuality){
				double quality=Tools.mid(2, 41, 10*Math.log10(maxWeight/Tools.max(0.01, (double)weightSum-maxWeight)));
				r[1]=(byte)quality;
			}
		}
		return r[0];
	}
	
	/** Returns a cached concrete-base frequency, or 0.25 for an ambiguous query.
	 * @param b Query symbol
	 * @return Cached A/C/G/T frequency; an empty concrete-base cache contains zeros */
	public float baseProb(byte b){
		assert(acgtProb!=null) : "Scoring requires calcProbs after accumulation at reference position "+rpos;
		int x=AminoAcid.baseToNumber[b];
		return x>=0 ? acgtProb[x] : 0.25f;
	}
	
	//inflection at x=.25,y=.25
//	public float baseScore(byte b){
//		float prob=baseProb(b);
//		return prob>0.25f ? prob : 5*prob-1;
//	}

	/** Returns the current untransformed base probability; older transforms remain below.
	 * @param b Query symbol
	 * @return Same value as baseProb(b) */
	public float baseScore(byte b){
		float prob=baseProb(b);
//		return prob>0.25f ? (prob-0.25f)*slope : 4*prob-1;
		
//		byte maxBase=maxBase();
//		byte secondBase=secondBase();
//		float prob1=baseProb(maxBase);
//		float prob2=baseProb(secondBase);
		
		return prob;
	}
	
	/** Returns a minimum-count concrete base; ties retain the first array index. */
	byte minBase(){
		int idx=0, count=acgtCount[0];
		for(int i=1; i<4; i++){
			if(acgtCount[i]<count){
				idx=i;
				count=acgtCount[i];
			}
		}
		return AminoAcid.numberToBase[idx];
	}
	
	/** Returns a maximum-count concrete base; ties retain the first array index. */
	byte maxBase(){
		int idx=0, count=acgtCount[0];
		for(int i=1; i<4; i++){
			if(acgtCount[i]>count){
				idx=i;
				count=acgtCount[i];
			}
		}
		return AminoAcid.numberToBase[idx];
	}
	
	/** Returns a base with the second-largest count; ties follow array order. */
	byte secondBase(){
		int idx=0, count=acgtCount[0];
		int idx2=-1, count2=-1;
		for(int i=1; i<4; i++){
			int c=acgtCount[i];
			if(c>count){
				idx2=idx;
				count2=count;
				idx=i;
				count=c;
			}else if(c>count2){
				idx2=i;
				count2=c;
			}
		}
		return AminoAcid.numberToBase[idx2];
	}
	
	/** Adds the observation through the same counter path as add(base, quality).
	 * @param base Observed symbol
	 * @param quality Observation weight */
	void increment(byte base, int quality){add(base, quality);}
	
	/** Returns "Node" as the part type identifier */
	@Override
	public final String partString(){return "Node";}
	
	@Override
	public ByteBuilder toText(){
		return toTextCount();
	}
	
	/** Formats allele weights and edge weight totals without changing state.
	 * @return A newly allocated diagnostic builder */
	public ByteBuilder toTextWeight(){
		ByteBuilder bb=new ByteBuilder();
		bb.append('(');
		bb.append(partString()).comma().append(rpos).comma().append(typeString());
		if(acgtCount==null){bb.append("[]");}
		else{
			for(int i=0; i<4; i++){bb.comma().append(acgtWeight[i]);}
		}
		bb.space();
		if(refEdge!=null){bb.comma().append("REF:").append(refEdge.weightSum);}
		if(insEdge!=null){bb.comma().append("INS:").append(insEdge.weightSum);}
		if(delEdge!=null){bb.comma().append("DEL:").append(delEdge.weightSum);}
		bb.append(')');
		return bb;
	}
	
	/** Formats allele counts and edge count totals without changing state.
	 * @return A newly allocated diagnostic builder */
	public ByteBuilder toTextCount(){
		ByteBuilder bb=new ByteBuilder();
		bb.append('(');
		bb.append(partString()).comma().append(rpos).comma().append(typeString());
		if(acgtCount==null){bb.append("[]");}
		else{
			for(int i=0; i<4; i++){bb.comma().append(acgtCount[i]);}
		}
		bb.space();
		if(refEdge!=null){bb.comma().append("REF:").append(refEdge.countSum);}
		if(insEdge!=null){bb.comma().append("INS:").append(insEdge.countSum);}
		if(delEdge!=null){bb.comma().append("DEL:").append(delEdge.countSum);}
		bb.append(')');
		return bb;
	}
	
	/** Compares weight, then observation count, then type using integer differences.
	 * @param b Other node
	 * @return The first nonzero difference, or zero for tied values */
	@Override
	public int compareTo(BaseNode b){
		int dif=weightSum-b.weightSum;
		if(dif!=0){return dif;}
		dif=countSum-b.countSum;
		if(dif!=0){return dif;}
		return type-b.type;
	}
	
	/** Constructs the one-time count-frequency cache for a REF/INS node.
	 * Concrete-base counts form the denominator; ambiguous observations are excluded.
	 * Callers must finish updates first and clear the cache explicitly before rebuilding. */
	void calcProbs(){
		assert(acgtProb==null) : "Probability caches are initialized once after counts are loaded at position "+rpos;
		acgtProb=new float[4];
		float mult=1f/Tools.max(1, simd.Vector.sum(acgtCount));
		for(int i=0; i<acgtProb.length; i++){
			acgtProb[i]=acgtCount[i]*mult;
		}
	}

	/** Returns (largest count - half the second-largest count) / concrete-base total.
	 * @return Diversity-based support, or zero when no concrete bases were observed */
	public float alleleDif(){
		int max=0, second=0, sum=0;
		for(int x : acgtCount){
			if(x>max){
				second=max;
				max=x;
			}else if(x>second){
				second=x;
			}
			sum+=x;
		}
//		return (max-second)/(float)Tools.max(1, sum);
		return (max-0.5f*second)/(float)Tools.max(1, sum);
//		return max/(float)Tools.max(1, sum);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Total observations, including ambiguous symbols and any caller pseudocount. */
	public int countSum;
	/** Sum of observation weights, including ambiguous symbols. */
	public int weightSum;
	
	/** Reference coordinate associated with this node. */
	public final int rpos;
	
	/** Stored reference symbol or insertion placeholder. */
	public final byte refBase;
	/** A/C/G/T weight totals; null for DEL. */
	public final int[] acgtWeight;
	/** A/C/G/T observation totals; null for DEL. */
	public final int[] acgtCount;
	/** One-time concrete-base frequencies, null until calcProbs. */
	public float acgtProb[];
	
	//These fields are not really necessary
	/** Next reference node, when connected by the owning graph. */
	public BaseNode refEdge;
	/** First/next insertion node at this position. */
	public BaseNode insEdge;
	/** Associated deletion node, when connected by the owning graph. */
	public BaseNode delEdge;
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Retained coefficient for the commented historical score transform. */
	private static final float slope=(4f/3f);
	
}
