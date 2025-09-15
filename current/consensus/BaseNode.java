package consensus;

import java.util.Arrays;

import dna.AminoAcid;
import shared.Tools;
import stream.FASTQ;
import structures.ByteBuilder;

/**
 * A placeholder for a base.
 * Tracks counts of bases seen at that position.
 * Maintains edges to observed next bases.
 * 
 * @author Brian Bushnell
 * @date September 6, 2019
 *
 */
public class BaseNode extends BaseGraphPart implements Comparable<BaseNode> {

	/**
	 * 
	 */
	private static final long serialVersionUID = 7932097131372307182L;
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public BaseNode(char refBase_, int type_, int rpos_){
		this((byte)refBase_, type_, rpos_);
	}
	
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
	
	/** Add a traversal of the designated base and quality */
	public void add(byte base, int quality){
		int num=AminoAcid.baseToNumber[base];
		if(num>=0 && type!=DEL){
			acgtWeight[num]+=quality;
			acgtCount[num]++;
		}
		countSum++;
		weightSum+=quality;
	}
	
//	public byte consensus() {
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
	
	public byte consensus(byte[] r) {
		assert(type!=DEL);
		if(onlyConvertNs && type==REF && refBase!='N'){
			r[0]=refBase;
			r[1]=20;
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
				r[1]=(refBase=='N' ? (byte)0 : (byte)2);
			}else{
				r[0]=b;
				double quality=Tools.mid(2, 41, 10*Math.log10(maxWeight/Tools.max(0.01f, weightSum)));
				r[1]=(byte)(quality+FASTQ.ASCII_OFFSET);
			}
			assert(r[0]!='N' || maxDepth<minDepth || af<maf) : "\n"+Arrays.toString(acgtCount)+", max="+maxDepth+", sum="+countSum+", maf="+maf+", af="+af+"\n"
					+ Arrays.toString(acgtWeight)+", max="+maxWeight+", weightSum="+weightSum+"\n";
		}else{
			assert(type==INS);
			r[0]=AminoAcid.numberToBase[maxPos];
			double quality=Tools.mid(2, 41, 10*Math.log10(maxWeight/Tools.max(0.01f, weightSum)));
			r[1]=(byte)(quality+FASTQ.ASCII_OFFSET);
		}
		return r[0];
	}
	
	public float baseProb(byte b){
		assert(acgtProb!=null) : this;
		int x=AminoAcid.baseToNumber[b];
		return x>=0 ? acgtProb[x] : 0.25f;
	}
	
	//inflection at x=.25,y=.25
//	public float baseScore(byte b){
//		float prob=baseProb(b);
//		return prob>0.25f ? prob : 5*prob-1;
//	}

	//inflection at x=.25,y=0.  Range is 1 to -1.
	public float baseScore(byte b){
		float prob=baseProb(b);
//		return prob>0.25f ? (prob-0.25f)*slope : 4*prob-1;
		
//		byte maxBase=maxBase();
//		byte secondBase=secondBase();
//		float prob1=baseProb(maxBase);
//		float prob2=baseProb(secondBase);
		
		return prob;
	}
	
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
	
	/** Add a traversal of the designated quality */
	void increment(byte base, int quality){add(base, quality);}
	
	@Override
	public final String partString(){return "Node";}
	
	@Override
	public ByteBuilder toText(){
		return toTextCount();
	}
	
	public ByteBuilder toTextWeight(){
		ByteBuilder bb=new ByteBuilder();
		bb.append('(');
		bb.append(partString()).comma().append(rpos).comma().append(typeString());
		if(acgtCount==null){bb.append("[]");}
		else {
			for(int i=0; i<4; i++){bb.comma().append(acgtWeight[i]);}
		}
		bb.space();
		if(refEdge!=null){bb.comma().append("REF:").append(refEdge.weightSum);}
		if(insEdge!=null){bb.comma().append("INS:").append(insEdge.weightSum);}
		if(delEdge!=null){bb.comma().append("DEL:").append(delEdge.weightSum);}
		bb.append(')');
		return bb;
	}
	
	public ByteBuilder toTextCount(){
		ByteBuilder bb=new ByteBuilder();
		bb.append('(');
		bb.append(partString()).comma().append(rpos).comma().append(typeString());
		if(acgtCount==null){bb.append("[]");}
		else {
			for(int i=0; i<4; i++){bb.comma().append(acgtCount[i]);}
		}
		bb.space();
		if(refEdge!=null){bb.comma().append("REF:").append(refEdge.countSum);}
		if(insEdge!=null){bb.comma().append("INS:").append(insEdge.countSum);}
		if(delEdge!=null){bb.comma().append("DEL:").append(delEdge.countSum);}
		bb.append(')');
		return bb;
	}
	
	@Override
	public int compareTo(BaseNode b) {
		int dif=weightSum-b.weightSum;
		if(dif!=0){return dif;}
		dif=countSum-b.countSum;
		if(dif!=0){return dif;}
		return type-b.type;
	}
	
	void calcProbs(){
		assert(acgtProb==null);
		acgtProb=new float[4];
		float mult=1f/Tools.max(1, shared.Vector.sum(acgtCount));
		for(int i=0; i<acgtProb.length; i++){
			acgtProb[i]=acgtCount[i]*mult;
		}
	}

	//Difference between first and second most common bases, as a fraction of the total traversals
	public float alleleDif() {
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
	
	/** Number of times this node has been traversed */
	public int countSum;
	/** Sum of scores of traversals.  Generally, sum of quality scores of this or adjacent bases.
	 * Probably equal to sum of acgt. */
	public int weightSum;
	
	public final int rpos;
	
	public final byte refBase;
	public final int[] acgtWeight;
	public final int[] acgtCount;
	public float acgtProb[];
	
	//These fields are not really necessary
	public BaseNode refEdge;
	public BaseNode insEdge;
	public BaseNode delEdge;
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	private static final float slope=(4f/3f);
	
}
