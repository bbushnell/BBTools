package fun;

import java.util.Random;

import shared.Parse;
import shared.Shared;
import shared.Tools;

public class Chance {
	
	//Probability of something with a chance of X happening at least Y times in Z chances
	public static void main(String[] args){
		
		int draws;
		int minSuccess;
		float prob;
		long rounds;
		try {
			draws = Parse.parseIntKMG(args[0]);
			minSuccess = Parse.parseIntKMG(args[1]);
			prob = Float.parseFloat(args[2]);
			rounds = Parse.parseKMG(args[3]);
		} catch (Exception e) {
			System.err.println("Chance (int)draws (int)minSuccess (float)prob (int)rounds");
			System.exit(1);
			throw new RuntimeException();
		}
		
		Random randy=Shared.threadLocalRandom();
		
		long passes=0;
		for(long i=0; i<rounds; i++){
			int pass=runOneRound(randy, draws, minSuccess, prob);
			passes+=pass;
		}
		
		double odds=passes*1.0/rounds;
		System.err.println("Probability: "+Tools.format("%.6f%%", 100*odds));
	}

	private static int runOneRound(Random randy, int draws, int minSuccess, float prob) {
		int success=0;
		for(int i=0; i<draws && success<minSuccess; i++){
			if(randy.nextFloat()<=prob){success++;}
		}
		return (success>=minSuccess ? 1 : 0);
	}
	
}
