package fun;

import java.util.Random;

import shared.Timer;
import structures.ByteBuilder;

/**
 * Benchmark different long-to-string append strategies.
 * @author Isla & Brian
 */
public class BenchLongAppend{
	
	public static void main(String[] args){
		int iterations=100_000;
		int numValues=1000;
		
		if(args.length>0){iterations=Integer.parseInt(args[0]);}
		if(args.length>1){numValues=Integer.parseInt(args[1]);}
		
		System.out.println("Generating "+numValues+" random longs...");
		
		// Generate random test values with distribution: (rand^2)*MAX_VALUE
		Random rand=new Random(42);
		long[] testValues=new long[numValues];
		for(int i=0; i<numValues; i++){
			double r=rand.nextDouble();
			long val=(long)(r*r*3000000/*Long.MAX_VALUE*/);
			// 1/4 negative
			if(rand.nextInt(4)==0){val=-val;}
			testValues[i]=val;
		}
		
		System.out.println("Benchmarking "+iterations+" iterations...\n");
		
		// Warmup
		for(int i=0; i<10000; i++){
			ByteBuilder bb=new ByteBuilder();
			for(long v : testValues){
				bb.clear();
				appendOriginal(bb, v);
				bb.clear();
				appendViaToString(bb, v);
				bb.clear();
				appendForward(bb, v);
				bb.clear();
				appendBinarySearch(bb, v);
			}
		}
		
		// Test 1: Original (reverse with numbuffer)
		Timer t1=new Timer();
		for(int i=0; i<iterations; i++){
			ByteBuilder bb=new ByteBuilder();
			for(long v : testValues){
				bb.clear();
				appendOriginal(bb, v);
			}
		}
		t1.stop();
		
		// Test 2: Long.toString approach
		Timer t2=new Timer();
		for(int i=0; i<iterations; i++){
			ByteBuilder bb=new ByteBuilder();
			for(long v : testValues){
				bb.clear();
				appendViaToString(bb, v);
			}
		}
		t2.stop();
		
		// Test 3: Forward write with digit counting
		Timer t3=new Timer();
		for(int i=0; i<iterations; i++){
			ByteBuilder bb=new ByteBuilder();
			for(long v : testValues){
				bb.clear();
				appendForward(bb, v);
			}
		}
		t3.stop();
		
		// Test 4: Binary search digit counting
		Timer t4=new Timer();
		for(int i=0; i<iterations; i++){
			ByteBuilder bb=new ByteBuilder();
			for(long v : testValues){
				bb.clear();
				appendBinarySearch(bb, v);
			}
		}
		t4.stop();
		
		// Results
		System.out.println("Original (reverse):     "+t1+String.format(" (%.2fx)", t1.elapsed/(double)t1.elapsed));
		System.out.println("Long.toString():        "+t2+String.format(" (%.2fx)", t1.elapsed/(double)t2.elapsed));
		System.out.println("Forward (linear):       "+t3+String.format(" (%.2fx)", t1.elapsed/(double)t3.elapsed));
		System.out.println("Binary search:          "+t4+String.format(" (%.2fx)", t1.elapsed/(double)t4.elapsed));
		
		// Verify correctness (first 5 values)
		System.out.println("\nVerifying correctness (first 5 values):");
		for(int i=0; i<5; i++){
			long v=testValues[i];
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
	
	static void appendOriginal(ByteBuilder bb, long x){
		bb.expand(20);
		if(x<0){
			if(x==Long.MIN_VALUE){
				bb.append("-9223372036854775808");
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
			long y=x%100;
			x=x/100;
			numbuffer[pos]=ones100[(int)y];
			pos++;
			numbuffer[pos]=tens100[(int)y];
			pos++;
		}
		while(x>0){
			long y=x%10;
			x=x/10;
			numbuffer[pos]=ones100[(int)y];
			pos++;
		}
		
		while(pos>0){
			pos--;
			bb.array[bb.length]=numbuffer[pos];
			bb.length++;
		}
	}
	
	// Method 2: Long.toString
	static void appendViaToString(ByteBuilder bb, long x){
		String s=Long.toString(x);
		int len=s.length();
		bb.expand(len);
		for(int i=0; i<len; i++){
			bb.array[bb.length++]=(byte)s.charAt(i);
		}
	}
	
	// Method 3: Forward write with linear digit counting
	static void appendForward(ByteBuilder bb, long x){
		if(x==0){
			bb.expand(1);
			bb.array[bb.length++]='0';
			return;
		}
		
		boolean negative=(x<0);
		if(negative){
			if(x==Long.MIN_VALUE){
				bb.append("-9223372036854775808");
				return;
			}
			x=-x;
		}
		
		// Count digits - linear comparison
		int digits=(x<10L) ? 1 : (x<100L) ? 2 : (x<1000L) ? 3 : 
		           (x<10000L) ? 4 : (x<100000L) ? 5 : (x<1000000L) ? 6 :
		           (x<10000000L) ? 7 : (x<100000000L) ? 8 : 
		           (x<1000000000L) ? 9 : (x<10000000000L) ? 10 :
		           (x<100000000000L) ? 11 : (x<1000000000000L) ? 12 :
		           (x<10000000000000L) ? 13 : (x<100000000000000L) ? 14 :
		           (x<1000000000000000L) ? 15 : (x<10000000000000000L) ? 16 :
		           (x<100000000000000000L) ? 17 : (x<1000000000000000000L) ? 18 : 19;
		
		bb.expand(digits+(negative?1:0));
		
		if(negative){
			bb.array[bb.length++]='-';
		}
		
		int pos=bb.length+digits-1;
		while(x>=100){
			long q=x/100;
			int r=(int)(x-(q*100));
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
	static void appendBinarySearch(ByteBuilder bb, long x){
		if(x==0){
			bb.expand(1);
			bb.array[bb.length++]='0';
			return;
		}
		
		boolean negative=(x<0);
		if(negative){
			if(x==Long.MIN_VALUE){
				bb.append("-9223372036854775808");
				return;
			}
			x=-x;
		}
		
		// Count digits - binary search
		int digits;
		if(x<100000000L){ // 1-8 digits
			if(x<10000L){
				if(x<100L){digits=(x<10L) ? 1 : 2;}
				else{digits=(x<1000L) ? 3 : 4;}
			}else{
				if(x<1000000L){digits=(x<100000L) ? 5 : 6;}
				else{digits=(x<10000000L) ? 7 : 8;}
			}
		}else{ // 9-19 digits
			if(x<10000000000000L){
				if(x<1000000000L){digits=9;}
				else{
					if(x<100000000000L){digits=(x<10000000000L) ? 10 : 11;}
					else{digits=(x<1000000000000L) ? 12 : 13;}
				}
			}else{
				if(x<1000000000000000L){
					digits=(x<100000000000000L) ? 14 : 15;
				}else{
					if(x<100000000000000000L){digits=(x<10000000000000000L) ? 16 : 17;}
					else{digits=(x<1000000000000000000L) ? 18 : 19;}
				}
			}
		}
		
//	// Version 1: Single if, two long ternaries
//		if(x<100000000L){ // 1-8 digits
//			digits=(x<100L) ? ((x<10L) ? 1 : 2) :
//			       (x<10000L) ? ((x<1000L) ? 3 : 4) :
//			       (x<1000000L) ? ((x<100000L) ? 5 : 6) :
//			       ((x<10000000L) ? 7 : 8);
//		}else{ // 9-19 digits
//			digits=(x<1000000000L) ? 9 :
//			       (x<100000000000L) ? ((x<10000000000L) ? 10 : 11) :
//			       (x<10000000000000L) ? ((x<1000000000000L) ? 12 : 13) :
//			       (x<1000000000000000L) ? ((x<100000000000000L) ? 14 : 15) :
//			       (x<100000000000000000L) ? ((x<10000000000000000L) ? 16 : 17) :
//			       ((x<1000000000000000000L) ? 18 : 19);
//		}

//		// Version 2: Two if levels before ternaries
//		if(x<100000000L){ // 1-8 digits
//			if(x<10000L){ // 1-4 digits
//				digits=(x<100L) ? ((x<10L) ? 1 : 2) : ((x<1000L) ? 3 : 4);
//			}else{ // 5-8 digits
//				digits=(x<1000000L) ? ((x<100000L) ? 5 : 6) : ((x<10000000L) ? 7 : 8);
//			}
//		}else{ // 9-19 digits
//			if(x<10000000000000L){ // 9-13 digits
//				digits=(x<1000000000L) ? 9 :
//				       (x<100000000000L) ? ((x<10000000000L) ? 10 : 11) :
//				       ((x<1000000000000L) ? 12 : 13);
//			}else{ // 14-19 digits
//				digits=(x<1000000000000000L) ? ((x<100000000000000L) ? 14 : 15) :
//				       (x<100000000000000000L) ? ((x<10000000000000000L) ? 16 : 17) :
//				       ((x<1000000000000000000L) ? 18 : 19);
//			}
//		}
		
		bb.expand(digits+(negative?1:0));
		
		if(negative){
			bb.array[bb.length++]='-';
		}
		
		int pos=bb.length+digits-1;
		while(x>=100){
			long q=x/100;
			int r=(int)(x-(q*100));
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