package cardinality;

import shared.Tools;

/**
 * SetSketch with 16-bit registers (Ertl 2021).
 * Encodes register values as floor(-log_b(u)) where u is a uniform hash.
 * With b=1.001 and q=65534, uses the same memory as DDL (2 bytes/register).
 * Implemented for empirical collision rate and speed comparison with DDL.
 *
 * @author Noire
 * @date May 2026
 */
public final class SetSketch16 extends CardinalityTracker {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	SetSketch16(){
		this(2048, 31, -1, 0);
	}

	SetSketch16(int buckets_, int k_, long seed, float minProb_){
		super(buckets_, k_, seed, minProb_);
		maxArray=new char[buckets];
		valueBits=64-bucketBits;
		valueScale=1.0/(1L<<valueBits);
	}

	@Override
	public SetSketch16 copy(){return new SetSketch16(buckets, k, -1, minProb);}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final long cardinality(){
		if(lastCardinality>=0){return lastCardinality;}
		double sum=0;
		int filled=0;
		for(int i=0; i<buckets; i++){
			final int v=maxArray[i];
			if(v>0){sum+=v; filled++;}
		}
		if(filled==0){lastCardinality=0; return 0;}
		long card=(long)(buckets*Math.pow(BASE, sum/buckets));
		card=Math.min(clampToAdded ? added : Long.MAX_VALUE, card);
		lastCardinality=card;
		return card;
	}

	@Override
	public final void add(CardinalityTracker log){
		assert(log.getClass()==this.getClass());
		add((SetSketch16)log);
	}

	public void add(SetSketch16 log){
		added+=log.added;
		lastCardinality=-1;
		if(maxArray!=log.maxArray){
			for(int i=0; i<buckets; i++){
				if(log.maxArray[i]>maxArray[i]){
					maxArray[i]=log.maxArray[i];
				}
			}
		}
	}

	@Override
	public final void hashAndStore(final long number){
		final long rawKey=number^hashXor;
		final long key=Tools.hash64shift(rawKey);
		final int bucket=(int)(key&bucketMask);

		final long micro=(key>>>bucketBits)&0x3FL;
		microIndex|=(1L<<micro);

		long remaining=key>>>bucketBits;
		if(remaining<=0){remaining=1;}
		final double u=remaining*valueScale;
		int score=(int)(-Math.log(u)/LOG_BASE);
		if(score<1){score=1;}
		if(score>MAX_REGISTER){score=MAX_REGISTER;}

		final char old=maxArray[bucket];
		if(score<=old){return;}
		lastCardinality=-1;
		maxArray[bucket]=(char)score;
	}

	@Override
	public int filledBuckets(){
		int count=0;
		for(char v : maxArray){if(v>0){count++;}}
		return count;
	}

	@Override
	public double occupancy(){return (double)filledBuckets()/buckets;}

	public int[] bucketValues(){
		int[] out=new int[buckets];
		for(int i=0; i<buckets; i++){out[i]=maxArray[i];}
		return out;
	}

	@Override
	public final float[] compensationFactorLogBucketsArray(){return null;}

	@Override
	public double[] rawEstimates(){
		throw new UnsupportedOperationException("SetSketch16 does not support rawEstimates()");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final char[] maxArray;
	private final int valueBits;
	private final double valueScale;

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	static final double BASE=1.001;
	static final double LOG_BASE=Math.log(BASE);
	static final int MAX_REGISTER=65534;

}
