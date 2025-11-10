package fun;

import shared.Timer;
import structures.ByteBuilder;

/**
 * Benchmark different integer-to-string append strategies.
 * @author Isla
 */
public class BenchmarkIntAppend {
	
	public static void main(String[] args){
		int iterations=10_000_000;
		
		// Test values covering different ranges
		int[] testValues={0, 5, -5, 42, -42, 999, -999, 12345, -12345, 
		                  999999, -999999, Integer.MAX_VALUE, Integer.MIN_VALUE};
		
		System.out.println("Benchmarking "+iterations+" iterations...\n");
		
		// Warmup
		for(int i=0; i<100000; i++){
			ByteBuilder bb=new ByteBuilder();
			for(int v : testValues){
				bb.clear();
				appendOriginal(bb, v);
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
		
		// Results
		System.out.println("Original (reverse):     "+t1);
		System.out.println("Integer.toString():     "+t2);
		System.out.println("Forward (digit count):  "+t3);
		
		// Verify correctness
		System.out.println("\nVerifying correctness:");
		for(int v : testValues){
			ByteBuilder bb1=new ByteBuilder();
			ByteBuilder bb2=new ByteBuilder();
			ByteBuilder bb3=new ByteBuilder();
			appendOriginal(bb1, v);
			appendViaToString(bb2, v);
			appendForward(bb3, v);
			String s1=bb1.toString();
			String s2=bb2.toString();
			String s3=bb3.toString();
			boolean match=(s1.equals(s2) && s2.equals(s3));
			System.out.println(v+": "+s1+" "+(match?"✓":"✗ MISMATCH: "+s2+", "+s3));
		}
	}
	
	// Original method (current ByteBuilder implementation)
	private static byte[] numbuffer=new byte[20];
	private static byte[] ones100, tens100;
	static {
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
	
	// Method 3: Forward write
	static void appendForward(ByteBuilder bb, int x){
		if(x==0){
			bb.expand(1);
			bb.array[bb.length++]='0';
			return;
		}
		
		boolean negative=(x<0);
		if(negative){
			if(x==Integer.MIN_VALUE){
				bb.append((long)x);
				return;
			}
			x=-x;
		}
		
		// Count digits
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
}