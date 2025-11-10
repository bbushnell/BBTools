package fun;

import java.util.Random;

import shared.Timer;
import structures.ByteBuilder;

/**
 * Benchmark different integer-to-string append strategies.
 * @author Isla & Brian
 */
public class BenchIntAppend2{
	
	public static void main(String[] args){
		int iterations=100_000;
		int numValues=1000;
		
		if(args.length>0){iterations=Integer.parseInt(args[0]);}
		if(args.length>1){numValues=Integer.parseInt(args[1]);}
		
		System.out.println("Generating "+numValues+" random integers...");
		
		// Generate random test values with distribution: (rand^2)*MAX_VALUE
		Random rand=new Random(42);
		int[] testValues=new int[numValues];
		for(int i=0; i<numValues; i++){
			double r=rand.nextDouble();
			int val=(int)(r*r*Integer.MAX_VALUE);
			// 1/4 negative
			if(rand.nextInt(4)==0){val=-val;}
			testValues[i]=val;
		}
		
		System.out.println("Benchmarking "+iterations+" iterations...\n");
		
		// Warmup
		for(int i=0; i<10000; i++){
			ByteBuilder bb=new ByteBuilder();
			for(int v : testValues){
				bb.clear();
				appendOriginal(bb, v);
				appendViaToString(bb, v);
				appendForward(bb, v);
				appendBinarySearch(bb, v);
			}
		}
		
		// Test 1: Original (reverse with numbuffer)
		Timer t1=new Timer();
		for(int i=0; i<iterations; i++){
			ByteBuilder bb=new ByteBuilder();
			for(int v : testValues){
				bb.clear();
				appendOriginal(bb, v);
			}
		}
		t1.stop();
		
		// Test 2: Integer.toString approach
		Timer t2=new Timer();
		for(int i=0; i<iterations; i++){
			ByteBuilder bb=new ByteBuilder();
			for(int v : testValues){
				bb.clear();
				appendViaToString(bb, v);
			}
		}
		t2.stop();
		
		// Test 3: Forward write with digit counting
		Timer t3=new Timer();
		for(int i=0; i<iterations; i++){
			ByteBuilder bb=new ByteBuilder();
			for(int v : testValues){
				bb.clear();
				appendForward(bb, v);
			}
		}
		t3.stop();
		
		// Test 4: Binary search digit counting
		Timer t4=new Timer();
		for(int i=0; i<iterations; i++){
			ByteBuilder bb=new ByteBuilder();
			for(int v : testValues){
				bb.clear();
				appendBinarySearch(bb, v);
			}
		}
		t4.stop();
		
		// Results
		System.out.println("Original (reverse):     "+t1+String.format(" (%.2fx)", t1.elapsed/(double)t1.elapsed));
		System.out.println("Integer.toString():     "+t2+String.format(" (%.2fx)", t1.elapsed/(double)t2.elapsed));
		System.out.println("Forward (linear):       "+t3+String.format(" (%.2fx)", t1.elapsed/(double)t3.elapsed));
		System.out.println("Binary search:          "+t4+String.format(" (%.2fx)", t1.elapsed/(double)t4.elapsed));
		
		// Verify correctness (first 5 values)
		System.out.println("\nVerifying correctness (first 5 values):");
		for(int i=0; i<5; i++){
			int v=testValues[i];
			ByteBuilder bb1=new ByteBuilder();
			ByteBuilder bb2=new ByteBuilder();
			ByteBuilder bb3=new ByteBuilder();
			ByteBuilder bb4=new ByteBuilder();
			appendOriginal(bb1, v);
			appendViaToString(bb2, v);
			appendForward(bb3, v);
			appendBinarySearch(bb4, v);
			String s1=bb1.toString();
			String s2=bb2.toString();
			String s3=bb3.toString();
			String s4=bb4.toString();
			boolean match=(s1.equals(s2) && s2.equals(s3) && s3.equals(s4));
			System.out.println(v+": "+s1+" "+(match?"✓":"✗"));
		}
	}
	
	// Original method (current ByteBuilder implementation)
	private static byte[] numbuffer=new byte[20];
	private static byte[] ones100, tens100;
	static{
		ones100=new byte[100];
		tens100=new byte[100];
		for(int i=0; i<100; i++){
			ones100[i]=(byte)('0'+(i%10));
			tens100[i]=(byte)('0'+(i/10));
		}
	}
	
	static void appendOriginal(ByteBuilder bb, int x){
		bb.expand(11);
		if(x<0){
			if(x<=Integer.MIN_VALUE){
				bb.append((long)x);
				return;
			}else{
				bb.array[bb.length]='-';
				bb.length++;
				x=-x;
			}
		}else if(x==0){
			bb.array[bb.length]='0';
			bb.length++;
			return;
		}
		
		int pos=0;
		while(x>9){
			int y=x%100;
			x=x/100;
			numbuffer[pos]=ones100[y];
			pos++;
			numbuffer[pos]=tens100[y];
			pos++;
		}
		while(x>0){
			int y=x%10;
			x=x/10;
			numbuffer[pos]=ones100[y];
			pos++;
		}
		
		while(pos>0){
			pos--;
			bb.array[bb.length]=numbuffer[pos];
			bb.length++;
		}
	}
	
	// Method 2: Integer.toString
	static void appendViaToString(ByteBuilder bb, int x){
		String s=Integer.toString(x);
		int len=s.length();
		bb.expand(len);
		for(int i=0; i<len; i++){
			bb.array[bb.length++]=(byte)s.charAt(i);
		}
	}
	
	// Method 3: Forward write with linear digit counting
	static void appendForward(ByteBuilder bb, int x){
		if(x==0){
			bb.expand(1);
			bb.array[bb.length++]='0';
			return;
		}
		
		boolean negative=(x<0);
		if(negative){
			if(x==Integer.MIN_VALUE){
				appendViaToString(bb, x);
				return;
			}
			x=-x;
		}
		
		// Count digits - linear comparison
		int digits=(x<10) ? 1 : (x<100) ? 2 : (x<1000) ? 3 : 
		           (x<10000) ? 4 : (x<100000) ? 5 : (x<1000000) ? 6 :
		           (x<10000000) ? 7 : (x<100000000) ? 8 : 
		           (x<1000000000) ? 9 : 10;
		
		bb.expand(digits+(negative?1:0));
		
		if(negative){
			bb.array[bb.length++]='-';
		}
		
		int pos=bb.length+digits-1;
		while(x>=100){
			int q=x/100;
			int r=x-(q*100);
			bb.array[pos--]=ones100[r];
			bb.array[pos--]=tens100[r];
			x=q;
		}
		while(x>0){
			bb.array[pos--]=(byte)('0'+(x%10));
			x/=10;
		}
		
		bb.length+=digits;
	}
	
	// Method 4: Binary search digit counting
	static void appendBinarySearch(ByteBuilder bb, int x){
		if(x==0){
			bb.expand(1);
			bb.array[bb.length++]='0';
			return;
		}
		
		boolean negative=(x<0);
		if(negative){
			if(x==Integer.MIN_VALUE){
				appendViaToString(bb, x);
				return;
			}
			x=-x;
		}
		
		// Count digits - binary search
		int digits;
//		if(x<100000){
//			if(x<100){
//				digits=(x<10) ? 1 : 2;
//			}else{
//				if(x<1000){digits=3;}
//				else{digits=(x<10000) ? 4 : 5;}
//			}
//		}else{
//			if(x<10000000){
//				digits=(x<1000000) ? 6 : 7;
//			}else{
//				if(x<100000000){digits=8;}
//				else{digits=(x<1000000000) ? 9 : 10;}
//			}
//		}
		if(x<100000){
			digits=(x<10) ? 1 : (x<100) ? 2 : (x<1000) ? 3 : 
	           (x<10000) ? 4 : 5;
		}else{
			digits=(x<1000000) ? 6 :
	           (x<10000000) ? 7 : (x<100000000) ? 8 : 
	           (x<1000000000) ? 9 : 10;
		}
		
		bb.expand(digits+(negative?1:0));
		
		if(negative){
			bb.array[bb.length++]='-';
		}
		
		int pos=bb.length+digits-1;
		while(x>=100){
			int q=x/100;
			int r=x-(q*100);
			bb.array[pos--]=ones100[r];
			bb.array[pos--]=tens100[r];
			x=q;
		}
		while(x>0){
			bb.array[pos--]=(byte)('0'+(x%10));
			x/=10;
		}
		
		bb.length+=digits;
	}
}