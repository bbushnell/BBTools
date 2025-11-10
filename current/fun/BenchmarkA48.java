package fun;

import java.util.Random;

import shared.Timer;
import structures.ByteBuilder;

/**
 * Benchmark A48 encoding: reverse vs forward approaches.
 * @author Isla
 */
public class BenchmarkA48 {
	
	public static void main(String[] args){
		int iterations=1000000;
		int arraySize=1000;
		
		// Generate test data: random positive integers with bias toward smaller values
		Random rand=new Random(12345);
		long[] testValues=new long[arraySize];
		for(int i=0; i<arraySize; i++){
			double r=rand.nextDouble();
			testValues[i]=(long)(Integer.MAX_VALUE * r * r); // Squared for smaller bias
		}
		
		System.out.println("Benchmarking "+iterations+" iterations on "+arraySize+" values...\n");
		
		byte[] temp=new byte[11]; // Max length for 64-bit value in base-64
		
		// Warmup
		for(int i=0; i<20000; i++){
			ByteBuilder bb=new ByteBuilder();
			for(long v : testValues){
				bb.clear();
				appendA48Reverse(bb, v, temp);
			}
		}
		
		// Test 1: Reverse (original)
		Timer t1=new Timer();
		for(int i=0; i<iterations; i++){
			ByteBuilder bb=new ByteBuilder();
			for(long v : testValues){
				bb.clear();
				appendA48Reverse(bb, v, temp);
			}
		}
		t1.stop();
		
		// Test 2: Forward (no reverse)
		Timer t2=new Timer();
		for(int i=0; i<iterations; i++){
			ByteBuilder bb=new ByteBuilder();
			for(long v : testValues){
				bb.clear();
				appendA48Forward(bb, v);
			}
		}
		t2.stop();
		
		// Results
		System.out.println("Reverse (original): "+t1);
		System.out.println("Forward (no temp):  "+t2);
		System.out.println();
		
		double speedup=t1.elapsed/(double)t2.elapsed;
		if(speedup>1.0){
			System.out.printf("Forward is %.2fx faster\n", speedup);
		}else{
			System.out.printf("Reverse is %.2fx faster\n", 1.0/speedup);
		}
		
		// Verify correctness
		System.out.println("\nVerifying correctness on sample values:");
		long[] samples={0, 1, 63, 64, 4095, 4096, Integer.MAX_VALUE, (long)Integer.MAX_VALUE+1};
		for(long v : samples){
			ByteBuilder bb1=new ByteBuilder();
			ByteBuilder bb2=new ByteBuilder();
			appendA48Reverse(bb1, v, temp);
			appendA48Forward(bb2, v);
			String s1=bb1.toString();
			String s2=bb2.toString();
			boolean match=s1.equals(s2);
			System.out.println(v+": \""+s1+"\" "+(match?"✓":"✗ MISMATCH: \""+s2+"\""));
		}
	}
	
	// Original: build reversed, then flip
	static void appendA48Reverse(ByteBuilder bb, long value, byte[] temp){
		int i=0;
		while(value!=0){
			byte b=(byte)(value&0x3F);
			temp[i]=b;
			value=value>>6;
			i++;
		}
		if(i==0){
			bb.append((byte)'0');
		}else{
			bb.expand(i);
			final byte[] array=bb.array;
			for(int j=bb.length, k=i-1, lim=bb.length+i; j<lim; j++, k--){
				array[j]=((byte)(temp[k]+48));
			}
			bb.length+=i;
		}
	}
	
	// Forward: calculate length, write directly
	static void appendA48Forward(ByteBuilder bb, long value){
		if(value==0){
			bb.append((byte)'0');
			return;
		}
		
		int highBit=63-Long.numberOfLeadingZeros(value);
		int symbols=(highBit/6)+1;
		
		bb.expand(symbols);
		for(int shift=(symbols-1)*6; shift>=0; shift-=6){
			byte b=(byte)((value>>shift)&0x3F);
			bb.array[bb.length++]=(byte)(b+48);
		}
	}
}