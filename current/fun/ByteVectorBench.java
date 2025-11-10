package fun;

import java.util.Random;

import dna.AminoAcid;
import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorSpecies;
import shared.Parse;
import shared.Timer;

/**
 * Benchmark quality offset application with different SIMD widths.
 * @author Isla&Brian
 */
public class ByteVectorBench {
	
	private static final VectorSpecies<Byte> BSPECIES64=ByteVector.SPECIES_64;
	private static final VectorSpecies<Byte> BSPECIES128=ByteVector.SPECIES_128;
	private static final VectorSpecies<Byte> BSPECIES256=ByteVector.SPECIES_256;
	
	private static final int BWIDTH64=BSPECIES64.length();
	private static final int BWIDTH128=BSPECIES128.length();
	private static final int BWIDTH256=BSPECIES256.length();
	
	public static void main(String[] args){
		int iterations=10000;
		int numArrays=1000;
		int minlen=100;
		int maxlen=160;
		int delta=33;
		if(args.length>0) {iterations=Parse.parseIntKMG(args[0]);}
		if(args.length>1) {numArrays=Parse.parseIntKMG(args[1]);}
		if(args.length>1) {minlen=Parse.parseIntKMG(args[2]);}
		if(args.length>1) {maxlen=Parse.parseIntKMG(args[3]);}
		
		
		System.out.println("Generating test data...");
		
		// Generate random-length arrays (100-160 bp)
		Random rand=new Random(42);
		byte[][] qualsArrays=new byte[numArrays][];
		byte[][] basesArrays=new byte[numArrays][];
		byte[] baseChoices={'A', 'C', 'G', 'T', 'N'};
		int range=maxlen-minlen+1;
		for(int i=0; i<numArrays; i++){
			int len=minlen+rand.nextInt(range);
			qualsArrays[i]=new byte[len];
			basesArrays[i]=new byte[len];
			
			for(int j=0; j<len; j++){
				qualsArrays[i][j]=(byte)(33+rand.nextInt(41)); // 33-73
				basesArrays[i][j]=baseChoices[rand.nextInt(5)];
			}
		}
		
		System.out.println("Generated "+numArrays+" arrays of length "+minlen+"-"+maxlen);
		System.out.println("Running "+iterations+" iterations...\n");
		
		// Make copies for each test
		byte[][][] qualsCopies=new byte[7][numArrays][];
		for(int test=0; test<7; test++){
			for(int i=0; i<numArrays; i++){
				qualsCopies[test][i]=qualsArrays[i].clone();
			}
		}
		
		// Warmup
		for(int iter=0; iter<iterations/8; iter++){
			for(int i=0; i<numArrays; i++){
				applyQualOffsetScalar(qualsCopies[0][i], basesArrays[i], delta);
				applyQualOffsetScalar(qualsCopies[0][i], basesArrays[i], -delta);
			}
		}
		
		// Warmup2
		for(int iter=0; iter<iterations/8; iter++){
			for(int i=0; i<numArrays; i++){
				applyQualOffset64(qualsCopies[1][i], basesArrays[i], delta);
				applyQualOffset64(qualsCopies[1][i], basesArrays[i], -delta);
				applyQualOffset128(qualsCopies[2][i], basesArrays[i], delta);
				applyQualOffset128(qualsCopies[2][i], basesArrays[i], -delta);
				applyQualOffset256(qualsCopies[3][i], basesArrays[i], delta);
				applyQualOffset256(qualsCopies[3][i], basesArrays[i], -delta);
				applyQualOffsetDual(qualsCopies[4][i], basesArrays[i], delta);
				applyQualOffsetDual(qualsCopies[4][i], basesArrays[i], -delta);
				applyQualOffsetDual2(qualsCopies[6][i], basesArrays[i], delta);
				applyQualOffsetDual2(qualsCopies[6][i], basesArrays[i], -delta);
				applyQualOffsetDispatch(qualsCopies[5][i], basesArrays[i], delta);
				applyQualOffsetDispatch(qualsCopies[5][i], basesArrays[i], -delta);
			}
		}
		
		// Test 1: Scalar
		Timer t1=new Timer();
		for(int iter=0; iter<iterations; iter++){
			for(int i=0; i<numArrays; i++){
				applyQualOffsetScalar(qualsCopies[0][i], basesArrays[i], delta);
				applyQualOffsetScalar(qualsCopies[0][i], basesArrays[i], -delta);
			}
		}
		t1.stop();
		
		// Test 2: SIMD 64-bit
		Timer t2=new Timer();
		for(int iter=0; iter<iterations; iter++){
			for(int i=0; i<numArrays; i++){
				applyQualOffset64(qualsCopies[1][i], basesArrays[i], delta);
				applyQualOffset64(qualsCopies[1][i], basesArrays[i], -delta);
			}
		}
		t2.stop();
		
		// Test 3: SIMD 128-bit
		Timer t3=new Timer();
		for(int iter=0; iter<iterations; iter++){
			for(int i=0; i<numArrays; i++){
				applyQualOffset128(qualsCopies[2][i], basesArrays[i], delta);
				applyQualOffset128(qualsCopies[2][i], basesArrays[i], -delta);
			}
		}
		t3.stop();
		
		// Test 4: SIMD 256-bit
		Timer t4=new Timer();
		for(int iter=0; iter<iterations; iter++){
			for(int i=0; i<numArrays; i++){
				applyQualOffset256(qualsCopies[3][i], basesArrays[i], delta);
				applyQualOffset256(qualsCopies[3][i], basesArrays[i], -delta);
			}
		}
		t4.stop();
		
		// Test 5: Dual 256+128
		Timer t5=new Timer();
		for(int iter=0; iter<iterations; iter++){
			for(int i=0; i<numArrays; i++){
				applyQualOffsetDual(qualsCopies[4][i], basesArrays[i], delta);
				applyQualOffsetDual(qualsCopies[4][i], basesArrays[i], -delta);
			}
		}
		t5.stop();
		
		// Test 6: Intelligent dispatch
		Timer t6=new Timer();
		for(int iter=0; iter<iterations; iter++){
			for(int i=0; i<numArrays; i++){
				applyQualOffsetDispatch(qualsCopies[5][i], basesArrays[i], delta);
				applyQualOffsetDispatch(qualsCopies[5][i], basesArrays[i], -delta);
			}
		}
		t6.stop();
		
		// Test 7: Dual 256+64
		Timer t7=new Timer();
		for(int iter=0; iter<iterations; iter++){
			for(int i=0; i<numArrays; i++){
				applyQualOffsetDual2(qualsCopies[6][i], basesArrays[i], delta);
				applyQualOffsetDual2(qualsCopies[6][i], basesArrays[i], -delta);
			}
		}
		t7.stop();
		
		// Results
		System.out.println("Scalar:             "+t1);
		System.out.println("SIMD 64-bit:        "+t2+String.format(" (%.2fx)", t1.elapsed/(double)t2.elapsed));
		System.out.println("SIMD 128-bit:       "+t3+String.format(" (%.2fx)", t1.elapsed/(double)t3.elapsed));
		System.out.println("SIMD 256-bit:       "+t4+String.format(" (%.2fx)", t1.elapsed/(double)t4.elapsed));
		System.out.println("Dual 256+128:       "+t5+String.format(" (%.2fx)", t1.elapsed/(double)t5.elapsed));
		System.out.println("Dual 256+64:        "+t7+String.format(" (%.2fx)", t1.elapsed/(double)t7.elapsed));
		System.out.println("Smart dispatch3:    "+t6+String.format(" (%.2fx)", t1.elapsed/(double)t6.elapsed));
		
		// Verify correctness
		System.out.println("\nVerifying correctness (first 3 arrays):");
		for(int i=0; i<3; i++){
			boolean match=true;
			byte[] reference=qualsCopies[0][i];
			for(int test=1; test<6; test++){
				if(!java.util.Arrays.equals(reference, qualsCopies[test][i])){
					match=false;
					System.out.println("Array "+i+" test "+test+" MISMATCH!");
				}
			}
			if(match){
				System.out.println("Array "+i+" (len="+reference.length+"): ✓");
			}
		}
	}
	
	// Scalar version
	static void applyQualOffsetScalar(final byte[] quals, final byte[] bases, final int delta){
		if(quals == null){return;}
		final int length=quals.length;
		
		for(int i=0; i<length; i++){
			byte b=bases[i];
			int q=quals[i]+delta;
			q=(AminoAcid.baseToNumber[b]<0 ? 0 : Math.max(2, q));
			quals[i]=(byte)q;
		}
	}
	
	// SIMD 64-bit
	static void applyQualOffset64(final byte[] quals, final byte[] bases, final int delta){
		if(quals == null){return;}
		final int length=quals.length;
		final int limit=BSPECIES64.loopBound(length);
		
		int i=0;
		ByteVector vdelta=ByteVector.broadcast(BSPECIES64, (byte)delta);
		ByteVector vn=ByteVector.broadcast(BSPECIES64, (byte)'N');
		ByteVector vzero=ByteVector.broadcast(BSPECIES64, (byte)0);
		ByteVector vcap2=ByteVector.broadcast(BSPECIES64, (byte)2);
		
		for(; i<limit; i += BWIDTH64){
			ByteVector vquals=ByteVector.fromArray(BSPECIES64, quals, i);
			ByteVector vbases=ByteVector.fromArray(BSPECIES64, bases, i);
			ByteVector vresult=vquals.add(vdelta);
			vresult=vresult.max(vcap2);
			VectorMask<Byte> maskN=vbases.eq(vn);
			vresult=vresult.blend(vzero, maskN);
			vresult.intoArray(quals, i);
		}
		
		for(; i<length; i++){
			byte b=bases[i];
			int q=quals[i]+delta;
			q=(AminoAcid.baseToNumber[b]<0 ? 0 : Math.max(2, q));
			quals[i]=(byte)q;
		}
	}
	
	// SIMD 128-bit
	static void applyQualOffset128(final byte[] quals, final byte[] bases, final int delta){
		if(quals == null){return;}
		final int length=quals.length;
		final int limit=BSPECIES128.loopBound(length);
		
		int i=0;
		ByteVector vdelta=ByteVector.broadcast(BSPECIES128, (byte)delta);
		ByteVector vn=ByteVector.broadcast(BSPECIES128, (byte)'N');
		ByteVector vzero=ByteVector.broadcast(BSPECIES128, (byte)0);
		ByteVector vcap2=ByteVector.broadcast(BSPECIES128, (byte)2);
		
		for(; i<limit; i += BWIDTH128){
			ByteVector vquals=ByteVector.fromArray(BSPECIES128, quals, i);
			ByteVector vbases=ByteVector.fromArray(BSPECIES128, bases, i);
			ByteVector vresult=vquals.add(vdelta);
			vresult=vresult.max(vcap2);
			VectorMask<Byte> maskN=vbases.eq(vn);
			vresult=vresult.blend(vzero, maskN);
			vresult.intoArray(quals, i);
		}
		
		for(; i<length; i++){
			byte b=bases[i];
			int q=quals[i]+delta;
			q=(AminoAcid.baseToNumber[b]<0 ? 0 : Math.max(2, q));
			quals[i]=(byte)q;
		}
	}
	
	// SIMD 256-bit
	static void applyQualOffset256(final byte[] quals, final byte[] bases, final int delta){
		if(quals == null){return;}
		final int length=quals.length;
		final int limit=BSPECIES256.loopBound(length);
		
		int i=0;
		ByteVector vdelta=ByteVector.broadcast(BSPECIES256, (byte)delta);
		ByteVector vn=ByteVector.broadcast(BSPECIES256, (byte)'N');
		ByteVector vzero=ByteVector.broadcast(BSPECIES256, (byte)0);
		ByteVector vcap2=ByteVector.broadcast(BSPECIES256, (byte)2);
		
		for(; i<limit; i += BWIDTH256){
			ByteVector vquals=ByteVector.fromArray(BSPECIES256, quals, i);
			ByteVector vbases=ByteVector.fromArray(BSPECIES256, bases, i);
			ByteVector vresult=vquals.add(vdelta);
			vresult=vresult.max(vcap2);
			VectorMask<Byte> maskN=vbases.eq(vn);
			vresult=vresult.blend(vzero, maskN);
			vresult.intoArray(quals, i);
		}
		
		for(; i<length; i++){
			byte b=bases[i];
			int q=quals[i]+delta;
			q=(AminoAcid.baseToNumber[b]<0 ? 0 : Math.max(2, q));
			quals[i]=(byte)q;
		}
	}
	
	// Dual 256+128
	static void applyQualOffsetDual(final byte[] quals, final byte[] bases, final int delta){
		if(quals == null){return;}
		final int length=quals.length;
		
		int i=0;
		
		// 256-bit pass
		ByteVector vdelta256=ByteVector.broadcast(BSPECIES256, (byte)delta);
		ByteVector vn256=ByteVector.broadcast(BSPECIES256, (byte)'N');
		ByteVector vzero256=ByteVector.broadcast(BSPECIES256, (byte)0);
		ByteVector vcap2_256=ByteVector.broadcast(BSPECIES256, (byte)2);
		
		for(; i<length - BWIDTH256+1; i += BWIDTH256){
			ByteVector vquals=ByteVector.fromArray(BSPECIES256, quals, i);
			ByteVector vbases=ByteVector.fromArray(BSPECIES256, bases, i);
			ByteVector vresult=vquals.add(vdelta256);
			vresult=vresult.max(vcap2_256);
			VectorMask<Byte> maskN=vbases.eq(vn256);
			vresult=vresult.blend(vzero256, maskN);
			vresult.intoArray(quals, i);
		}
		
		// 128-bit pass
		ByteVector vdelta128=ByteVector.broadcast(BSPECIES128, (byte)delta);
		ByteVector vn128=ByteVector.broadcast(BSPECIES128, (byte)'N');
		ByteVector vzero128=ByteVector.broadcast(BSPECIES128, (byte)0);
		ByteVector vcap2_128=ByteVector.broadcast(BSPECIES128, (byte)2);
		
		for(; i<length - BWIDTH128+1; i += BWIDTH128){
			ByteVector vquals=ByteVector.fromArray(BSPECIES128, quals, i);
			ByteVector vbases=ByteVector.fromArray(BSPECIES128, bases, i);
			ByteVector vresult=vquals.add(vdelta128);
			vresult=vresult.max(vcap2_128);
			VectorMask<Byte> maskN=vbases.eq(vn128);
			vresult=vresult.blend(vzero128, maskN);
			vresult.intoArray(quals, i);
		}
		
		// Scalar tail
		for(; i<length; i++){
			byte b=bases[i];
			int q=quals[i]+delta;
			q=(AminoAcid.baseToNumber[b]<0 ? 0 : Math.max(2, q));
			quals[i]=(byte)q;
		}
	}
	
	// Dual 256+64
	static void applyQualOffsetDual2(final byte[] quals, final byte[] bases, final int delta){
		if(quals == null){return;}
		final int length=quals.length;
		
		int i=0;
		
		// 256-bit pass
		ByteVector vdelta256=ByteVector.broadcast(BSPECIES256, (byte)delta);
		ByteVector vn256=ByteVector.broadcast(BSPECIES256, (byte)'N');
		ByteVector vzero256=ByteVector.broadcast(BSPECIES256, (byte)0);
		ByteVector vcap2_256=ByteVector.broadcast(BSPECIES256, (byte)2);
		
		for(; i<length - BWIDTH256+1; i += BWIDTH256){
			ByteVector vquals=ByteVector.fromArray(BSPECIES256, quals, i);
			ByteVector vbases=ByteVector.fromArray(BSPECIES256, bases, i);
			ByteVector vresult=vquals.add(vdelta256);
			vresult=vresult.max(vcap2_256);
			VectorMask<Byte> maskN=vbases.eq(vn256);
			vresult=vresult.blend(vzero256, maskN);
			vresult.intoArray(quals, i);
		}
		
		// 64-bit pass
		ByteVector vdelta64=ByteVector.broadcast(BSPECIES64, (byte)delta);
		ByteVector vn64=ByteVector.broadcast(BSPECIES64, (byte)'N');
		ByteVector vzero64=ByteVector.broadcast(BSPECIES64, (byte)0);
		ByteVector vcap2_64=ByteVector.broadcast(BSPECIES64, (byte)2);
		
		for(; i<length - BWIDTH64+1; i += BWIDTH64){
			ByteVector vquals=ByteVector.fromArray(BSPECIES64, quals, i);
			ByteVector vbases=ByteVector.fromArray(BSPECIES64, bases, i);
			ByteVector vresult=vquals.add(vdelta64);
			vresult=vresult.max(vcap2_64);
			VectorMask<Byte> maskN=vbases.eq(vn64);
			vresult=vresult.blend(vzero64, maskN);
			vresult.intoArray(quals, i);
		}
		
		// Scalar tail
		for(; i<length; i++){
			byte b=bases[i];
			int q=quals[i]+delta;
			q=(AminoAcid.baseToNumber[b]<0 ? 0 : Math.max(2, q));
			quals[i]=(byte)q;
		}
	}
	
	// Smart dispatch
	static void applyQualOffsetDispatch(final byte[] quals, final byte[] bases, final int delta){
		if(quals == null){return;}
		int simd=pickByteSIMD3(quals.length);
		if(simd == 2){
			applyQualOffsetDual2(quals, bases, delta);
		}else if(simd == 1){
			applyQualOffset64(quals, bases, delta);
		}else{
			applyQualOffsetScalar(quals, bases, delta);
		}
	}
	
	static int pickByteSIMD(int len){
		int cost256=(len>>4)+(len&31);  // 2*(len/32)+len%32
		int cost128=(len>>3)+(len&15);  // 2*(len/16)+len%16
		return cost256<cost128 ? 2 : cost128<len ? 1 : 0;
	}
	
	static int pickByteSIMD3(int len){
//		int cost256=(len>>4)+((len&31)>>2)+(len&7);
//		int cost64=(len>>2)+(len&7);
//		return cost256<cost64 ? 2 : cost64<len ? 1 : 0;
		return len<8 ? 0 : len<32 ? 1 : 2;
	}
}