package fun;

import java.util.Arrays;
import shared.Random;

/**
 * Genetic algorithm implementation for solving optimization problems using
 * evolutionary computation techniques. Uses population-based evolution with
 * fitness-proportionate selection, crossover, and mutation to find optimal
 * solutions in binary-encoded search spaces.
 *
 * @author Brian Bushnell
 */
public class Genetic {
	
	/**
	 * Program entry point that runs the genetic algorithm optimization.
	 * Creates a Genetic instance, solves the optimization problem, and prints
	 * the best solution found in both binary and decimal representation.
	 *
	 * @param args Command-line arguments: [population_size] [bits] [iterations]
	 * [mutation_probability]
	 */
	public static void main(String[] args){
		Genetic g=new Genetic(args);
		long answer=g.solve();
		System.out.println(Long.toBinaryString(answer)+" \t-> "+f(answer));
	}
	
	/**
	 * Constructs a Genetic algorithm instance with parameters from command-line
	 * arguments. Sets default values if arguments are not provided.
	 *
	 * @param args Command-line arguments: [0]=population size (default 20),
	 * [1]=bits per individual (default 8), [2]=iterations (default 20),
	 * [3]=mutation probability (default 0.01)
	 */
	public Genetic(String[] args){
		pop=(args.length>0 ? Integer.parseInt(args[0]) : 20);
		bits=(args.length>1 ? Integer.parseInt(args[1]) : 8);
		iters=(args.length>2 ? Integer.parseInt(args[2]) : 20);
		mutProb=(args.length>3 ? Double.parseDouble(args[3]) : 0.01);
		mask=(bits>63 ? -1L : ~((-1L)<<bits));
	}
	
	/**
	 * Executes the genetic algorithm to find the optimal solution.
	 * Initializes random population, evolves through specified iterations using
	 * selection, crossover, and mutation, tracking the best individual across
	 * all generations.
	 *
	 * @return The DNA (bit string) of the best individual found
	 */
	public long solve(){
		final long mask=(bits>63 ? -1L : ~((-1L)<<bits));
		double[] prob=new double[pop];
		Bug[] current=new Bug[pop];
		for(int i=0; i<pop; i++){
			current[i]=new Bug(randy.nextLong()&mask);
		}
		Arrays.sort(current);
		Bug best=current[current.length-1];
		
		for(int i=0; i<iters; i++){
			
			if(true){
				System.out.println("Iteration "+i+": "+current[current.length-1]);
			}
			
			current=iterate(current, prob);
			Arrays.sort(current);
			if(best.compareTo(current[current.length-1])<0){best=current[current.length-1];}
		}
		
		return best.dna;
	}
	
	/**
	 * Performs one generation of evolution on the current population.
	 * Calculates fitness-proportionate selection probabilities and creates
	 * the next generation through breeding and mutation.
	 *
	 * @param current The current generation of individuals
	 * @param prob Probability array for selection (reused for efficiency)
	 * @return The next generation of individuals
	 */
	public Bug[] iterate(Bug[] current, double[] prob){
		Arrays.fill(prob,  0);
		double sum=0;
		for(int i=0; i<current.length; i++){
			Bug b=current[i];
			sum+=b.fitness;
			prob[i]=b.fitness;
		}
		double mult=1/sum;
		prob[0]*=mult;
		
		for(int i=1; i<prob.length; i++){
			prob[i]=prob[i-1]+prob[i]*mult;
		}
		
		Bug[] next=new Bug[current.length];
		for(int i=0; i<next.length; i++){
			long babyDna=breed(current, prob, mutProb);
			next[i]=new Bug(babyDna);
		}
		return next;
	}
	
	/**
	 * Creates offspring through crossover and optional mutation.
	 * Selects two parents using fitness-proportionate selection, performs
	 * uniform crossover with random mask, and applies bit-flip mutation.
	 *
	 * @param current Population to select parents from
	 * @param prob Cumulative selection probabilities
	 * @param mutProb Probability of mutation occurring
	 * @return DNA of the offspring individual
	 */
	public long breed(Bug[] current, double[] prob, double mutProb){
		double fa=randy.nextDouble();
		double fb=randy.nextDouble();
		Bug a=current[findIndex(fa, prob)];
		Bug b=current[findIndex(fb, prob)];
		long crossover=randy.nextLong();
		long baby=(a.dna&crossover)|(b.dna&~crossover);
		if(mutProb>0 && randy.nextDouble()<mutProb){
			long bit=(1L<<randy.nextInt(bits));
			baby^=bit;
		}
		return baby;
	}
	
	/**
	 * Finds the population index for fitness-proportionate selection.
	 * Uses binary search-like approach to locate individual based on
	 * random fitness threshold and cumulative probabilities.
	 *
	 * @param f Random value between 0.0 and 1.0 for selection
	 * @param prob Cumulative probability distribution
	 * @return Index of selected individual
	 */
	public int findIndex(double f, double[] prob){
		for(int i=pop-1; i>0; i--){
			if(prob[i-1]<f){return i;}
		}
		return 0;
	}
	
	/**
	 * Fitness function for evaluating individuals.
	 * Currently implements simple quadratic function f(x) = x^2.
	 * Higher values indicate better fitness.
	 *
	 * @param x Individual's DNA as long integer
	 * @return Fitness score (higher is better)
	 */
	public static double f(long x){
		return x*x;
	}
	
	/**
	 * Represents an individual in the genetic algorithm population.
	 * Encapsulates genetic material (DNA) and corresponding fitness value
	 * with comparison capability for sorting by fitness.
	 */
	private static class Bug implements Comparable<Bug> {
		
		/**
		 * Creates a new individual with specified genetic material.
		 * Automatically calculates fitness using the global fitness function.
		 * @param dna_ The genetic material as bit string (long integer)
		 */
		public Bug(long dna_){
			dna=dna_;
			fitness=f(dna);
		}
		
		@Override
		public int compareTo(Bug b){
			return (fitness<b.fitness ? -1 : fitness>b.fitness ? 1 : 0);
		}
		
		@Override
		public String toString(){
			return Long.toBinaryString(dna)+" \t-> "+fitness;
		}
		
		/** The genetic material represented as a bit string */
		final long dna;
		/** The fitness value calculated from the DNA */
		final double fitness;
	}
	
	/** Random number generator used throughout the genetic algorithm */
	public static final Random randy=shared.Shared.random();

	/** Population size for the genetic algorithm */
	final int pop;
	/** Number of bits in each individual's genetic representation */
	final int bits;
	/** Number of generations to evolve the population */
	final int iters;
	/** Probability of mutation occurring during breeding */
	final double mutProb;
	
	/** Bit mask to constrain DNA values to the specified bit length */
	final long mask;
	
}
