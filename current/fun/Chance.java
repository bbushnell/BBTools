package fun;

import shared.Random;

import shared.Parse;
import shared.Shared;
import shared.Tools;

/**
 * Monte Carlo probability calculator for discrete success scenarios.
 * Simulates scenarios with multiple draws to determine the likelihood of achieving
 * a minimum number of successes given individual event probabilities.
 * @author Brian Bushnell
 */
public class Chance {
	
	//Probability of something with a chance of X happening at least Y times in Z chances
	/**
	 * Runs Monte Carlo simulation to calculate probability of success.
	 * Expects 4 arguments: draws, minSuccess, prob, rounds.
	 * @param args Command line arguments: [draws] [minSuccess] [prob] [rounds]
	 */
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

	/**
	 * Executes a single simulation round with early termination optimization.
	 * Stops drawing once minimum successes are reached to improve performance.
	 *
	 * @param randy Random number generator for simulation
	 * @param draws Maximum number of attempts allowed
	 * @param minSuccess Minimum successes needed to pass
	 * @param prob Individual success probability per draw
	 * @return 1 if minimum successes achieved, 0 otherwise
	 */
	private static int runOneRound(Random randy, int draws, int minSuccess, float prob) {
		int success=0;
		for(int i=0; i<draws && success<minSuccess; i++){
			if(randy.nextFloat()<=prob){success++;}
		}
		return (success>=minSuccess ? 1 : 0);
	}
	
}
