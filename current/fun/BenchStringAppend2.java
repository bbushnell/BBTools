package fun;

import java.lang.invoke.MethodHandles;
import java.lang.invoke.VarHandle;
import java.nio.charset.StandardCharsets;
import java.util.Random;

import shared.Timer;
import shared.Vector;
import structures.ByteBuilder;

/**
 * Benchmark different String append strategies for ByteBuilder.
 * @author Isla & Brian
 */
public class BenchStringAppend2 {
	
	private static final VarHandle STRING_VALUE_HANDLE;
	private static final VarHandle STRING_CODER_HANDLE;
	private static final boolean USE_VARHANDLES;

	static {

		// Try VarHandles (Java 9+)
		VarHandle tempValueHandle = null;
		VarHandle tempCoderHandle = null;
		boolean successVarHandles = false;

		try {
			MethodHandles.Lookup lookup = MethodHandles.privateLookupIn(String.class, MethodHandles.lookup());
			tempValueHandle = lookup.findVarHandle(String.class, "value", byte[].class);
			tempCoderHandle = lookup.findVarHandle(String.class, "coder", byte.class);
			System.out.println("VarHandles available: true");
			successVarHandles = true;
		} catch (Exception e) {
			System.out.println("VarHandles available: false");
		}

		STRING_VALUE_HANDLE = tempValueHandle;
		STRING_CODER_HANDLE = tempCoderHandle;
		USE_VARHANDLES = successVarHandles;
	}

	public static void main(String[] args){
		int iterations = 100000;
		int numStrings = 1000;
		int minLen = 8;
		int maxLen = 50;

		if(args.length > 0) iterations = Integer.parseInt(args[0]);
		if(args.length > 1) numStrings = Integer.parseInt(args[1]);
		if(args.length > 2) minLen = Integer.parseInt(args[2]);
		if(args.length > 3) maxLen = Integer.parseInt(args[3]);

		System.out.println("Generating " + numStrings + " random strings (length " + minLen + "-" + maxLen + ")...");

		// Generate random strings
		Random rand = new Random(42);
		String[] strings = new String[numStrings];
		for(int i = 0; i < numStrings; i++){
			int len = minLen + rand.nextInt(maxLen - minLen + 1);
			StringBuilder sb = new StringBuilder(len);
			for(int j = 0; j < len; j++){
				char c = (char)('0' + rand.nextInt('z' - '0' + 1)); // '0' to 'z'
				sb.append(c);
			}
			strings[i] = sb.toString();
		}

		System.out.println("Running " + iterations + " iterations...\n");

		// Warmup
		{
			ByteBuilder bb = new ByteBuilder(1000);
			for(int iter = 0; iter < 10000; iter++){
				bb.clear();
				for(String s : strings){
					appendOriginal(bb, s);
				}
			}
		}

		// Test 1: Original (charAt loop)
		Timer t1 = new Timer();
		{
			ByteBuilder bb = new ByteBuilder(1000);
			for(int iter = 0; iter < iterations; iter++){
				bb.clear();
				for(String s : strings){
					appendOriginal(bb, s);
				}
			}
		}
		t1.stop();

		// Test 2: getBytes with arraycopy
		Timer t2 = new Timer();
		{
			ByteBuilder bb = new ByteBuilder(1000);
			for(int iter = 0; iter < iterations; iter++){
				bb.clear();
				for(String s : strings){
					appendGetBytes(bb, s);
				}
			}
		}
		t2.stop();

		// Test 4: VarHandles
		Timer t4 = new Timer();
		{
			ByteBuilder bb = new ByteBuilder(1000);
			for(int iter = 0; iter < iterations; iter++){
				bb.clear();
				for(String s : strings){
					appendVarHandles(bb, s);
				}
			}
		}
		t4.stop();
		


		// Test 1: Original (charAt loop)
		t1.start();
		{
			ByteBuilder bb = new ByteBuilder(1000);
			for(int iter = 0; iter < iterations; iter++){
				bb.clear();
				for(String s : strings){
					appendOriginal(bb, s);
				}
			}
		}
		t1.stop();

		// Test 2: getBytes with arraycopy
		t2.start();
		{
			ByteBuilder bb = new ByteBuilder(1000);
			for(int iter = 0; iter < iterations; iter++){
				bb.clear();
				for(String s : strings){
					appendGetBytes(bb, s);
				}
			}
		}
		t2.stop();

		// Test 4: VarHandles
		t4.start();
		{
			ByteBuilder bb = new ByteBuilder(1000);
			for(int iter = 0; iter < iterations; iter++){
				bb.clear();
				for(String s : strings){
					appendVarHandles(bb, s);
				}
			}
		}
		t4.stop();

		// Results
		System.out.println("Original (charAt):  " + t1 + String.format(" (%.2fx)", t1.elapsed / (double)t1.elapsed));
		System.out.println("getBytes():         " + t2 + String.format(" (%.2fx)", t1.elapsed / (double)t2.elapsed));
		System.out.println("VarHandles:         " + t4 + String.format(" (%.2fx)", t1.elapsed / (double)t4.elapsed));

		// Verify correctness
		System.out.println("\nVerifying correctness (first 3 strings):");
		for(int i = 0; i < 3; i++){
			ByteBuilder bb1 = new ByteBuilder();
			ByteBuilder bb2 = new ByteBuilder();
			ByteBuilder bb4 = new ByteBuilder();

			appendOriginal(bb1, strings[i]);
			appendGetBytes(bb2, strings[i]);
			appendVarHandles(bb4, strings[i]);

			String s1 = bb1.toString();
			String s2 = bb2.toString();
			String s4 = bb4.toString();

			boolean match = s1.equals(s2) && s2.equals(s4);
			System.out.println("\"" + strings[i] + "\": " + (match ? "✓" : "✗ MISMATCH: " + s1 + " | " + s2 + " | " + s4));
		}
	}

	// Original: charAt loop
	static void appendOriginal(ByteBuilder bb, String x){
		if(x == null) return;
		bb.expand(x.length());
		for(int i = 0; i < x.length(); i++){
			bb.array[bb.length] = (byte)x.charAt(i);
			bb.length++;
		}
	}

	// getBytes with arraycopy
	static void appendGetBytes(ByteBuilder bb, String s){
		if(s == null) return;
		byte[] x = s.getBytes(StandardCharsets.ISO_8859_1);
		bb.expand(x.length);
		System.arraycopy(x, 0, bb.array, bb.length, x.length);
		bb.length += x.length;
//		bb.append(s);
	}

	// VarHandles (Java 9+)
	static void appendVarHandles(ByteBuilder bb, String x){
		Vector.append(bb, x);
//		if(x == null) return;
//
//		if(!USE_VARHANDLES){
//			appendOriginal(bb, x);
//			return;
//		}
//
//		byte[] value = (byte[]) STRING_VALUE_HANDLE.get(x);
//		byte coder = (byte) STRING_CODER_HANDLE.get(x);
//
//		if(coder == 0){ // LATIN1
//			bb.expand(value.length);
//			System.arraycopy(value, 0, bb.array, bb.length, value.length);
//			bb.length += value.length;
//		} else { // UTF-16 fallback
//			bb.expand(x.length());
//			for(int i = 0; i < x.length(); i++){
//				bb.array[bb.length++] = (byte)x.charAt(i);
//			}
//		}
	}
}