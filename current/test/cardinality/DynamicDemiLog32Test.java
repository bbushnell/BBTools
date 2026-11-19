package test.cardinality;

import java.util.Arrays;
import java.util.Random;

import cardinality.DynamicDemiLog;
import cardinality.DynamicDemiLog32;
import parse.LineParser1;
import shared.Shared;
import shared.Tools;
import simd.Vector;
import structures.ByteBuilder;

/** Independent bit-by-bit score oracle and unsigned register invariants.
 * @author Sayu
 */
public final class DynamicDemiLog32Test{

	public static void main(String[] args){
		assert(args!=null) : "Finite tests do not read external corpora";
		if(args.length!=0){throw new IllegalArgumentException("No arguments");}
		final Random random=new Random(32023);
		for(int exponent=1; exponent<=8; exponent++){
			for(int bit=0; bit<64; bit++){
				check(DynamicDemiLog32.scoreHash(1L<<bit,exponent)==oracle(1L<<bit,exponent),"Single-bit hash boundary");
			}
			check(DynamicDemiLog32.scoreHash(0,exponent)==oracle(0,exponent),"Zero hash must not wrap its rank or shift distance");
			for(int trial=0; trial<10000; trial++){
				final long hash=random.nextLong();
				check(DynamicDemiLog32.scoreHash(hash,exponent)==oracle(hash,exponent),"Independent score oracle");
			}
		}
		for(int width : new int[]{1,2,8,512,1024}){
			final DynamicDemiLog32 whole=DynamicDemiLog32.create(width,5,19,0,true,true,6);
			final DynamicDemiLog32 left=whole.copy(),right=whole.copy();
			final int[] expected=new int[width];
			for(int i=0; i<20000; i++){
				final long word=random.nextLong(), hash=Tools.hash64shift(word^19);
				final int bucket=(int)(hash&(width-1)),score=oracle(hash,6);
				if(Integer.compareUnsigned(score,expected[bucket])>0){expected[bucket]=score;}
				whole.add(word); (i%2==0 ? left : right).add(word);
			}
			check(Arrays.equals(expected,whole.maxArray32()),"Adaptive early exit equals nonadaptive maxima");
			left.add(right);
			check(Arrays.equals(expected,left.maxArray32()),"Partition merge equals serial sketch");
			int filled=0,floor=64;
			for(int value : expected){if(value!=0){filled++;}floor=Math.min(floor,value>>>26);}
			check(whole.filledBuckets()==filled && whole.getGlobalNLZ()==floor,"Filled and floor invariants");
			check(whole.cardinality()>=filled && whole.cardinality()<=20000,"Cardinality bounds with real add counters");
			if(width>1){
				final int[] folded=new int[width/2];
				for(int i=0; i<width; i++){int b=i&(folded.length-1);if(Integer.compareUnsigned(expected[i],folded[b])>0){folded[b]=expected[i];}}
				check(Arrays.equals(folded,whole.condenseSketch(width/2,false).maxArray32()),"Unsigned folding");
			}
		}
		final int[] boundaries={0,1,65535,65536,0x7fffffff,0x80000000,0xfffffffe,0xffffffff};
		final DynamicDemiLog32 wide=DynamicDemiLog32.fromArray(boundaries,5,0,6,null);
		check(wide.filledBuckets()==7 && wide.registerValue(1)>>>16==0,"Populated prefix0 is not empty");
		final ByteBuilder bytes=wide.toBytes(new ByteBuilder());
		final LineParser1 lp=new LineParser1('\t');lp.set(bytes.toBytes());
		check(lp.terms()==boundaries.length,"A48 width count");
		for(int i=0; i<boundaries.length; i++){check(lp.parseLongA48(i)==Integer.toUnsignedLong(boundaries[i]),"A48 uint32 lossless including high bit");}
		boolean refused=false;
		try{wide.maxArray();}catch(IllegalArgumentException expected){refused=true;}
		check(refused,"Wrong-width access throws, never truncates");
		refused=false;
		try{wide.as16();}catch(IllegalArgumentException expected){refused=true;}
		check(refused,"Legacy clade attachment must reject a wide sketch rather than truncate it");
		final DynamicDemiLog narrow=DynamicDemiLog.create(8,5,0,0,false);
		check(narrow.as16()==narrow,"Checked legacy conversion preserves the original sketch object");
		refused=false;
		try{wide.requireCompatible(DynamicDemiLog32.fromArray(boundaries,5,1,6,null));}catch(IllegalArgumentException expected){refused=true;}
		check(refused,"Seed mismatch rejected");
		refused=false;
		try{wide.requireCompatible(DynamicDemiLog.create(8,5,0,0,false));}catch(IllegalArgumentException expected){refused=true;}
		check(refused,"Mixed width rejected");
		final map.LongHashSet filter=new map.LongHashSet();filter.add(11);
		final DynamicDemiLog32 folding=DynamicDemiLog32.fromArray(new int[]{0x80000001,0,0x40000001,0},5,0,6,new long[]{11,0,12,0});
		final DynamicDemiLog32 preferred=folding.condenseSketch(2,filter);
		check(preferred.registerValue(0)==0x40000001 && preferred.kmerArray()[0]==12,"Explicit blacklist folding prefers a nonblacklisted witness without signed comparison");
		check(folding.registerValue(0)==0x80000001,"Folding must not mutate its source");
		final int[] fineA=new int[8],fineB=new int[8];
		Arrays.fill(fineA,1<<26);Arrays.fill(fineB,(1<<26)|65535);
		final double a=DynamicDemiLog32.fromArray(fineA,5,0,6,null).rawEstimates()[2];
		final double b=DynamicDemiLog32.fromArray(fineB,5,0,6,null).rawEstimates()[2];
		check(a!=b && Double.isFinite(a) && Double.isFinite(b),"Statistics retain lower16 mantissa precision");
		final boolean previous=Shared.SIMD;
		for(int n=0; n<70; n++){
			final int[] x=new int[n],y=new int[n];
			for(int i=0; i<n; i++){x[i]=boundaries[i%8];y[i]=(i%4==0 ? x[i] : boundaries[(i*3+1)%8]);}
			final int[] expected=new int[4];
			for(int i=0; i<n; i++){
				final long u=Integer.toUnsignedLong(x[i]),v=Integer.toUnsignedLong(y[i]);
				expected[u==0 && v==0 ? 3 : u<v ? 0 : u>v ? 2 : 1]++;
			}
			Shared.SIMD=false;check(Arrays.equals(expected,Vector.compareDDL(x,y)),"Unsigned scalar counts");
			Shared.SIMD=true;check(Arrays.equals(expected,Vector.compareDDL(x,y)),"Unsigned SIMD or capability-gated scalar counts");
			if(Vector.simd256){
				final int[] signed=new int[4];
				for(int i=0; i<n; i++){signed[x[i]==0 && y[i]==0 ? 3 : x[i]<y[i] ? 0 : x[i]>y[i] ? 2 : 1]++;}
				check(Arrays.equals(signed,simd.SIMDLogLog.compareDetailed(x,y)),"Existing signed public comparison remains signed");
			}
		}
		Shared.SIMD=previous;
		System.err.println("DYNAMIC_DEMILOG32_TEST_PASS checks="+checks+" simd="+Vector.simd256);
	}

	/** Build exponent and inverted mantissa one bit at a time; no normalization shifts. */
	private static int oracle(long hash, int exponent){
		assert(exponent>=1 && exponent<=8) : "Finite oracle uses the documented exponent domain";
		int zeros=0;
		while(zeros<64 && ((hash>>>(63-zeros))&1)==0){zeros++;}
		final int bits=32-exponent,limit=(1<<exponent)-1;
		long mantissa=0;
		for(int i=0; i<bits; i++){
			final int position=62-zeros-i;
			final long bit=position<0 ? 0 : (hash>>>position)&1;
			mantissa=(mantissa<<1)|(bit^1);
		}
		if(zeros>limit){mantissa=(1L<<bits)-1;}
		return (int)Math.max(1L,((long)Math.min(zeros,limit)<<bits)|mantissa);
	}
	private static void check(boolean condition,String message){
		checks++;if(!condition){throw new AssertionError(message);}
	}
	private static int checks;
}
