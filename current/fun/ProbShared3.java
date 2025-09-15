package fun;

import java.util.HashSet;
import java.util.Random;

public class ProbShared3 {

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
		fillRandomSet(k, len2, set);
		final long space=(long)Math.pow(4, k);
		final int kmers=len1-k+1;
		for(int i=0; i<kmers; i++){
			long kmer=(randy.nextLong()&Long.MAX_VALUE)%space;
			if(set.contains(kmer)){return 1;}
		}
		return 0;
	}
	
	static void fillRandomSet(int k, int len, HashSet<Long> set){
		set.clear();
		final long space=(long)Math.pow(4, k);
		final int kmers=len-k+1;
		for(int i=0; i<kmers; i++){
			set.add((randy.nextLong()&Long.MAX_VALUE)%space);
		}
	}
	
	static final Random randy=new Random();
	
}
