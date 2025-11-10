package fun;

import java.util.Arrays;
import java.util.Random;

import dna.AminoAcid;
import shared.Timer;
import shared.Tools;
import shared.Vector;

/**
 * Benchmark and validate in-place reverse operations.
 * @author Isla & Brian
 */
public class BenchmarkReverse{
	
	public static void main(String[] args){
		int iterations=1000;
		int arraySize=1000;
		
		if(args.length>0){iterations=Integer.parseInt(args[0]);}
		if(args.length>1){arraySize=Integer.parseInt(args[1]);}
		
		System.out.println("Validating correctness across different array sizes...\n");
		
		//Test various sizes including edge cases
		int[] testSizes={0, 1, 7, 8, 9, 31, 32, 33, 63, 64, 65, 127, 128, 255, 256, 1000, 10000};
		boolean allPassed=true;
		
		for(int size : testSizes){
			boolean passed=validateReverse(size) && validateReverseComplement(size);
			System.out.println("Size "+size+": "+(passed ? "✓" : "✗ FAILED"));
			allPassed&=passed;
		}
		
		if(!allPassed){
			System.err.println("\nValidation FAILED! Not running benchmarks.");
			return;
		}
		
		System.out.println("\nAll validation passed!\n");
		System.out.println("Benchmarking with "+iterations+" iterations, array size "+arraySize+"...\n");
		
		//Generate test data
		Random rand=new Random(42);
		byte[][] testArraysReverse=new byte[arraySize][];
		byte[][] testArraysRevComp=new byte[arraySize][];
		
		for(int i=0; i<1000; i++){
			int len=1+Math.min(rand.nextInt(1000), rand.nextInt(1000));
			testArraysReverse[i]=randomBytes(rand, len);
			testArraysRevComp[i]=randomBases(rand, len);
		}
		
		//Warmup
		for(int i=0; i<10000; i++){
			byte[] temp=testArraysRevComp[i%1000];
			assert(Tools.min(temp)>0);
			reverseScalar(temp);
			assert(Tools.min(temp)>0);
			Vector.reverseInPlace(temp);
			assert(Tools.min(temp)>0);
			reverseComplementScalar(temp);
			assert(Tools.min(temp)>0);
			Vector.reverseComplementInPlace(temp);
			assert(Tools.min(temp)>0);
			Vector.reverseComplementInPlaceFast(temp);
			assert(Tools.min(temp)>0);
		}
		
		//Benchmark reverse
		Timer t1=new Timer();
		for(int i=0; i<iterations; i++){
			for(byte[] temp : testArraysReverse){
				reverseScalar(temp);
			}
		}
		t1.stop();
		
		Timer t2=new Timer();
		for(int i=0; i<iterations; i++){
			for(byte[] temp : testArraysReverse){
				Vector.reverseInPlace(temp);
			}
		}
		t2.stop();
		
		//Benchmark reverse-complement
		Timer t3=new Timer();
		for(int i=0; i<iterations; i++){
			for(byte[] temp : testArraysRevComp){
				reverseComplementScalar(temp);
			}
		}
		t3.stop();
		
		Timer t4=new Timer();
		for(int i=0; i<iterations; i++){
			for(byte[] temp : testArraysRevComp){
				Vector.reverseComplementInPlace(temp);
			}
		}
		t4.stop();
		
		Timer t5=new Timer();
		for(int i=0; i<iterations; i++){
			for(byte[] temp : testArraysRevComp){
				Vector.reverseComplementInPlaceFast(temp);
			}
		}
		t5.stop();
		
		//Results
		System.out.println("Reverse:");
		System.out.println("  Scalar:   "+t1+String.format(" (%.2fx)", t1.elapsed/(double)t1.elapsed));
		System.out.println("  SIMD:     "+t2+String.format(" (%.2fx)", t1.elapsed/(double)t2.elapsed));
		System.out.println();
		System.out.println("Reverse-complement:");
		System.out.println("  Scalar:   "+t3+String.format(" (%.2fx)", t3.elapsed/(double)t3.elapsed));
		System.out.println("  SIMD:     "+t4+String.format(" (%.2fx)", t3.elapsed/(double)t4.elapsed));
		System.out.println("  SIMDFast: "+t5+String.format(" (%.2fx)", t3.elapsed/(double)t5.elapsed));
	}
	
	static boolean validateReverse(int size){
		Random rand=new Random(size);
		byte[] original=randomBytes(rand, size);
		byte[] scalar=Arrays.copyOf(original, size);
		byte[] simd=Arrays.copyOf(original, size);
		
		reverseScalar(scalar);
		Vector.reverseInPlace(simd);
		
		return Arrays.equals(scalar, simd);
	}
	
	static boolean validateReverseComplement(int size){
		Random rand=new Random(size+1000);
		byte[] original=randomBases(rand, size);
		byte[] scalar=Arrays.copyOf(original, size);
		byte[] simd=Arrays.copyOf(original, size);
		
		reverseComplementScalar(scalar);
		Vector.reverseComplementInPlace(simd);
		
		boolean match=Arrays.equals(scalar, simd);
		if(!match && size<20){
			System.err.println("  Original: "+new String(original));
			System.err.println("  Scalar:   "+new String(scalar));
			System.err.println("  SIMD:     "+new String(simd));
		}
		return match;
	}
	
	static byte[] randomBytes(Random rand, int size){
		byte[] array=new byte[size];
		rand.nextBytes(array);
		return array;
	}
	
	static byte[] randomBases(Random rand, int size){
//		byte[] bases={'A','C','G','T','a','c','g','t','N'};
		byte[] bases={'A','C','G','T','N'};
		byte[] array=new byte[size];
		for(int i=0; i<size; i++){
			array[i]=bases[rand.nextInt(bases.length)];
		}
		return array;
	}
	
	static void reverseScalar(byte[] array){
		if(array==null || array.length<2){return;}
		int left=0;
		int right=array.length-1;
		while(left<right){
			byte temp=array[left];
			array[left]=array[right];
			array[right]=temp;
			left++;
			right--;
		}
	}
	
	static void reverseComplementScalar(byte[] array){
		if(array==null || array.length<1){return;}
		int left=0;
		int right=array.length-1;
		while(left<right){
			byte bLeft=array[left];
			byte bRight=array[right];
			array[left]=AminoAcid.baseToComplementExtended[bRight];
			array[right]=AminoAcid.baseToComplementExtended[bLeft];
			left++;
			right--;
		}
		if(left==right){
			array[left]=AminoAcid.baseToComplementExtended[array[left]];
		}
	}
}