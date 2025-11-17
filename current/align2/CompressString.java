package align2;

import dna.ChromosomeArray;
import dna.Data;
import structures.IntList;

/**
 * Utility class for compressing strings by identifying and removing repetitive patterns.
 * Provides various algorithms for detecting and compacting repeated sequences with
 * different period lengths and compression strategies. Used for reducing memory
 * footprint of genomic data with repetitive elements.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class CompressString {
	
	/**
	 * Program entry point for testing compression algorithms.
	 * Demonstrates various compression methods on command-line input and chromosome data.
	 * @param args Command-line arguments where args[0] is the string to compress
	 */
	public static void main(String[] args){
		
		String s;
		
		s=compressRepeats(args[0].getBytes(), 1);
		s=compress(args[0]);
		s=compressRepeatsUltra(args[0].getBytes(), 1, 3, null);
		System.out.println(args[0]+"\n"+s);
		
		System.exit(0);
		
		ChromosomeArray cha=Data.getChromosome(1);
		byte[] bytes=cha.array;
		
		int letters=0;
		for(int i=0; i<bytes.length; i++){if(bytes[i]!='N'){letters++;}}

		System.out.println("cha bytes length = "+bytes.length);
		System.out.println("cha letters length = "+letters);
		System.out.println("min="+cha.minIndex+", max="+cha.maxIndex+", length="+(cha.maxIndex-cha.minIndex+1));
		
		s=compressRepeatsUltra(bytes, 1, 1, null);
		System.out.println("compress(1) length: "+s.length());
		
		s=compressRepeatsUltra(bytes, 2, 2, null);
		System.out.println("compress(2) length: "+s.length());
		
		s=compressRepeatsUltra(bytes, 3, 3, null);
		System.out.println("compress(3) length: "+s.length());
		
		s=compressRepeatsUltra(bytes, 1, 2, null);
		System.out.println("compress(1,2) length: "+s.length());
		
		s=compressRepeatsUltra(bytes, 1, 3, null);
		System.out.println("compress(1,3) length: "+s.length());
		
//		s=compressRepeatsMultiperiod(bytes, 1, 1, null);
//		System.out.println("compress(1) length: "+s.length());
//
//		s=compressRepeatsMultiperiod(bytes, 2, 2, null);
//		System.out.println("compress(2) length: "+s.length());
//
//		s=compressRepeatsMultiperiod(bytes, 3, 3, null);
//		System.out.println("compress(3) length: "+s.length());
//
//		s=compressRepeatsMultiperiod(bytes, 1, 2, null);
//		System.out.println("compress(1,2) length: "+s.length());
//
//		s=compressRepeatsMultiperiod(bytes, 1, 3, null);
//		System.out.println("compress(1,3) length: "+s.length());
	}
	
	/**
	 * Compresses a string using multiple passes with increasing period lengths.
	 * Applies compression with periods 1, 2, and 3 sequentially for maximum compression.
	 * @param s The string to compress
	 * @return Compressed string with repetitive patterns removed
	 */
	public static String compress(String s){
		String s1=compressRepeats(s.getBytes(), 1);
		String s2=compressRepeats(s1.getBytes(), 2);
		String s3=compressRepeats(s2.getBytes(), 3);
		return s3;
	}
	
	/**
	 * Compresses repetitive patterns in byte array using single fixed period length.
	 * Identifies repeated sequences of specified period length and applies logarithmic
	 * compression based on repeat count. Skips compression for 0-2 repeats.
	 *
	 * @param array The byte array to compress
	 * @param period The period length for pattern detection
	 * @return Compressed string representation
	 */
	public static String compressRepeats(byte[] array, int period){
		
		StringBuilder sb=new StringBuilder(array.length);
		
		for(int base=0; base<array.length; base++){
			//Test for repeats of current pattern (array[i] to array[period-1])
			int repeats=countRepeats(array, base, period);
			int occurances=repeats+1;

//			System.out.println("base = "+base+"\t, repeats = "+repeats);
			
			if(repeats==0){
				//Advance pointer by 1
				sb.append((char)array[base]);
			}else if(repeats==1){
				//Still advance pointer by 1
				sb.append((char)array[base]);
			}else if(repeats==2){
				//Jump ahead
				base+=period-1;
			}else{
				//Compress and advance pointer by a factor of period
				
				int log=(32-Integer.numberOfLeadingZeros(repeats+1))-1;  //+1 is optional; gives lower compression
//				System.out.println("repeats="+repeats+
//						", Integer.highestOneBit("+repeats+")="+Integer.highestOneBit(repeats)+
//						", log="+log+", "+Integer.toBinaryString(repeats));
				assert(log>0 && log<=31);
				
				
				//Append
				for(int i=1; i<log; i++){
					for(int j=0; j<period; j++){
						sb.append((char)array[base+j]);
					}
				}
				
				base=base+(period*(repeats))-1;
				
			}
			
		}
		
		return sb.toString();
	}
	
	/**
	 * Compresses repetitive patterns using variable period lengths within specified range.
	 * Scans for repeats from minimum to maximum period length, compressing first match found.
	 * Uses logarithmic compression strategy for significant repeat counts.
	 *
	 * @param array The byte array to compress
	 * @param minPeriod Minimum period length to test
	 * @param maxPeriod Maximum period length to test
	 * @param list Optional IntList to track compressed positions (may be null)
	 * @return Compressed string representation
	 */
	public static String compressRepeatsMultiperiod(byte[] array, int minPeriod, int maxPeriod, IntList list){
		
		StringBuilder sb=new StringBuilder(array.length);
		
		for(int base=0; base<array.length; base++){
			//Test for repeats of current pattern (array[i] to array[period-1])
			
			
			int period=0;
			int repeats=0;
//			for(int x=maxPeriod; x>=minPeriod; x--){
//				int temp=countRepeats(array, base, x);
//				if(temp>1){
//					repeats=temp;
//					period=x;
//					break;
//				}
//			}
			for(int x=minPeriod; x<=maxPeriod; x++){
				int temp=countRepeats(array, base, x);
				if(temp>1){
					repeats=temp;
					period=x;
					break;
				}
			}
			int occurances=repeats+1;

//			System.out.println("base = "+base+"\t, repeats = "+repeats+"\t, period = "+period);
			
			if(repeats==0){
				//Advance pointer by 1
				sb.append((char)array[base]);
				if(list!=null){list.add(base);}
			}else if(repeats==1){
				//Still advance pointer by 1
				sb.append((char)array[base]);
				if(list!=null){list.add(base);}
			}
//			else if(repeats==2){
//				for(int j=0; j<period; j++){
//					sb.append((char)array[base+j]);
//				}
//				//Jump ahead
//				base+=2*period-1;
//			}
			else{
				//Compress and advance pointer by a factor of period
				
				int log=(32-Integer.numberOfLeadingZeros(repeats+1))-1;  //+1 is optional; gives lower compression
//				System.out.println("repeats="+repeats+
//						", Integer.highestOneBit("+repeats+")="+Integer.highestOneBit(repeats)+
//						", log="+log+", "+Integer.toBinaryString(repeats));
				assert(log>0 && log<=31);
				
				
				//Append
				for(int i=0; i<log; i++){
					for(int j=0; j<period; j++){
						sb.append((char)array[base+j]);
						if(list!=null){list.add(base+j);}
					}
				}
				
				base=base+(period*(repeats))-1;
				
			}
			
		}
		
		return sb.toString();
	}
	
	/**
	 * Optimized compression algorithm for repetitive patterns with minimal overhead.
	 * Uses simplified compression strategy that stores only one copy of repeated patterns
	 * regardless of repeat count. More aggressive compression than other methods.
	 *
	 * @param array The byte array to compress
	 * @param minPeriod Minimum period length to test
	 * @param maxPeriod Maximum period length to test
	 * @param list Optional IntList to track compressed positions (may be null)
	 * @return Compressed string representation with maximum compression
	 */
	public static String compressRepeatsUltra(byte[] array, int minPeriod, int maxPeriod, IntList list){
		
		StringBuilder sb=new StringBuilder(array.length);
		
		for(int base=0; base<array.length; base++){
			//Test for repeats of current pattern (array[i] to array[period-1])
			
			
			int period=0;
			int repeats=0;
			
			for(int x=minPeriod; x<=maxPeriod; x++){
				int temp=countRepeats(array, base, x);
//				System.out.println("*** temp="+temp+" for "+base+", "+x);
				if(temp>1){
					repeats=temp;
					period=x;
					break;
				}
			}
//			System.err.println(repeats);
			if(repeats==0){
				//Advance pointer by 1
				sb.append((char)array[base]);
				if(list!=null){list.add(base);}
			}else if(repeats==1){
				//Still advance pointer by 1
				sb.append((char)array[base]);
				if(list!=null){list.add(base);}

//				System.err.println(base);
				base=base+(period*(repeats))-1;
//				System.err.println(base);
			}
//			else if(repeats==2){
//				for(int j=0; j<period; j++){
//					sb.append((char)array[base+j]);
//				}
//				//Jump ahead
//				base+=2*period-1;
//			}
			else{
				//Compress and advance pointer
				
				//Append
				for(int j=0; j<period; j++){
					sb.append((char)array[base+j]);
					if(list!=null){list.add(base+j);}
				}
				
				base=base+(period*(repeats))-1;
				
			}
			
		}
		
		return sb.toString();
	}
	
	/**
	 * Counts consecutive repetitions of a pattern starting at specified base position.
	 * Compares pattern of given period length against subsequent occurrences in array
	 * until mismatch is found or end of array is reached.
	 *
	 * @param array The byte array to analyze
	 * @param base Starting position for pattern comparison
	 * @param period Length of pattern to match
	 * @return Number of complete repetitions found (0 if no repeats)
	 */
	public static int countRepeats(byte[] array, int base, int period){
		
		int max=array.length-period+1;
		
		int matches=0;
		boolean fail=false;
		for(int loc=base+period; loc<max && !fail; loc+=period){
			for(int i=0; i<period && !fail; i++){
				if(array[base+i]==array[loc+i]){
					matches++;
//					System.out.println("base = "+base+", loc = "+loc+", period = "+period+", and "+(base+i)+" == "+(loc+i));
				}else{
//					System.err.println("failed");
					fail=true;
				}
			}
		}
//		System.err.println("matches = "+matches);
		int repeats=matches/period;
//		System.err.println("repeats = "+repeats);
		return repeats;
	}
	
	
}
