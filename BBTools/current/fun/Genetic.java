package fun;

import java.util.Arrays;
import java.util.Random;

public class Genetic {
	
	public static void main(String[] args){
		Genetic g=new Genetic(args);
		long answer=g.solve();
		System.out.println(Long.toBinaryString(answer)+" \t-> "+f(answer));
	}
	
	public Genetic(String[] args){
		pop=(args.length>0 ? Integer.parseInt(args[0]) : 20);
		bits=(args.length>1 ? Integer.parseInt(args[1]) : 8);
		iters=(args.length>2 ? Integer.parseInt(args[2]) : 20);
		mutProb=(args.length>3 ? Double.parseDouble(args[3]) : 0.01);
		mask=(bits>63 ? -1L : ~((-1L)<<bits));
	}
	
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
	
	public int findIndex(double f, double[] prob){
		for(int i=pop-1; i>0; i--){
			if(prob[i-1]<f){return i;}
		}
		return 0;
	}
	
	public static double f(long x){
		return x*x;
	}
	
	private static class Bug implements Comparable<Bug> {
		
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
		
		final long dna;
		final double fitness;
	}
	
	public static final Random randy=new Random();

	final int pop;
	final int bits;
	final int iters;
	final double mutProb;
	
	final long mask;
	
}
