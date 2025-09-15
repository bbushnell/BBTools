package fun;

import java.util.HashSet;
import java.util.Random;

import dna.AminoAcid;

public class ProbShared2 {

	public static void main(String args[]){
		int k=Integer.parseInt(args[0]);
		int len1=Integer.parseInt(args[1]);
		int len2=Integer.parseInt(args[2]);
		int rounds=Integer.parseInt(args[3]);
		
		System.out.println("Probability:   "+simulate(k, len1, len2, rounds));
	}
	
	static double simulate(int k, int len1, int len2, int rounds){
		int successes=0;
		final HashSet<Long> set=new HashSet<Long>();
		for(int i=0; i<rounds; i++){
			successes+=simulateOnePair(k, len1, len2, set);
		}
		return successes/(double)rounds;
	}
	
	static int simulateOnePair(int k, int len1, int len2, HashSet<Long> set){
		set.clear();
		
		final int shift=2*k;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		long kmer=0;
		int len=0;
		
		byte[] bases=randomSequence(len2);
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x<0){len=0;}else{len++;}
			if(len>=k){
				set.add(kmer);
			}
		}
		
		bases=randomSequence(len1);
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=baseToNumber[b];
			kmer=((kmer<<2)|x)&mask;
			if(x<0){len=0;}else{len++;}
			if(len>=k){
				if(set.contains(kmer)){return 1;}
			}
		}
		return 0;
	}
	
	static byte[] randomSequence(int len){
		byte[] array=new byte[len];
		for(int i=0; i<len; i++){
			int number=randy.nextInt(4);
			array[i]=AminoAcid.numberToBase[number];
		}
		return array;
	}
	
	static final Random randy=new Random();
	static final byte[] numberToBase=AminoAcid.numberToBase;
	static final byte[] baseToNumber=AminoAcid.baseToNumber;
	
}
