package jgi;

import java.util.ArrayList;
import java.util.Random;

import barcode.CountBarcodes;
import dna.AminoAcid;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Jul 10, 2014
 *
 */
public class GreedyBarCodeFinder {
	
	/**
	 * Program entry point for barcode finding.
	 * Expects arguments: k-mer length, minimum hamming distance, optional rounds count.
	 * Reports the maximum number of valid barcodes found and execution time.
	 * @param args Command-line arguments [k, hamming_distance, rounds]
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		
		GreedyBarCodeFinder finder=new GreedyBarCodeFinder(args);
		int best=finder.find(finder.rounds);
		
		t.stop();
		System.err.println("There are at least "+best+" codes of length "+finder.k+" with mutual hamming distance at least "+finder.hdist);
		System.err.println("Time: \t"+t);
	}
	
	/**
	 * Constructs a GreedyBarCodeFinder with parameters from command-line arguments.
	 * Parses k-mer length, minimum hamming distance, and optionally number of rounds.
	 * @param args Array containing [k, hamming_distance, rounds] where rounds defaults to 20
	 */
	public GreedyBarCodeFinder(String[] args){
		k=Integer.parseInt(args[0]);
		hdist=Integer.parseInt(args[1]);
		rounds=(args.length>2 ? Integer.parseInt(args[2]) : 20);
	}
	
	/**
	 * Finds the maximum number of valid barcodes using multiple search strategies.
	 * Runs the deterministic algorithm once, then performs multiple randomized rounds
	 * to find the best result across all attempts.
	 *
	 * @param rounds Number of randomized search rounds to perform
	 * @return Maximum number of valid barcodes found across all rounds
	 */
	public int find(int rounds){
		ArrayList<String> list=new ArrayList<String>(1024);
		final int space=1<<(2*k);
		
		int[] set=new int[(int)space];
		if(set!=null){
			set=new int[(int)space];
			for(int i=0; i<set.length; i++){set[i]=i;}
		}
		
		int best=mainOld(k, hdist, list);
		for(int i=0; i<rounds; i++){
			best=Tools.max(best, test(k, hdist, set, list));
		}
		return best;
	}
	
	/**
	 * Deterministic greedy barcode selection algorithm.
	 * Iterates through all possible k-mers in sequential order, adding each
	 * k-mer to the barcode set only if it maintains the minimum hamming distance
	 * from all previously selected barcodes.
	 *
	 * @param k Length of k-mers to generate as barcodes
	 * @param hdist Minimum hamming distance required between any two barcodes
	 * @param list Collection to store selected barcode strings (cleared if not null)
	 * @return Number of valid barcodes found
	 */
	static int mainOld(int k, int hdist, ArrayList<String> list){

		final long space=1L<<(2*k);
		if(list==null){list=new ArrayList<String>(1024);}
		else{list.clear();}

		for(long kmer=0; kmer<space; kmer++){
			String s=AminoAcid.kmerToString(kmer, k);
			int dist=CountBarcodes.calcHdist(s, list);
			if(dist>=hdist){list.add(s);}
		}

		return list.size();

	}
	
	/**
	 * Randomized greedy barcode selection algorithm.
	 * Creates a shuffled array of all possible k-mers, then applies the same
	 * greedy selection process as mainOld but in random order. Different random
	 * orderings can yield different numbers of valid barcodes.
	 *
	 * @param k Length of k-mers to generate as barcodes
	 * @param hdist Minimum hamming distance required between any two barcodes
	 * @param set Pre-allocated array for shuffling k-mer indices
	 * @param list Collection to store selected barcode strings (cleared if not null)
	 * @return Number of valid barcodes found in this randomized round
	 */
	static int test(int k, int hdist, int[] set, ArrayList<String> list){
		
		final int space=1<<(2*k);
		if(set!=null){
			set=new int[(int)space];
			for(int i=0; i<set.length; i++){set[i]=i;}
		}
		Random randy=Shared.threadLocalRandom();
		for(int i=0; i<set.length; i++){
			int x=i+randy.nextInt(set.length-i);
			int temp=set[i];
			set[i]=set[x];
			set[x]=temp;
		}
		
		if(list==null){list=new ArrayList<String>(1024);}
		else{list.clear();}
		
		for(long kmer : set){
			String s=AminoAcid.kmerToString(kmer, k);
			int dist=CountBarcodes.calcHdist(s, list);
			if(dist>=hdist){list.add(s);}
		}
		
		return list.size();
	}
	
	/** Length of k-mers to use as barcodes */
	private final int k;
	/** Minimum hamming distance required between any two barcodes */
	private final int hdist;
	/** Number of randomized search rounds to perform */
	private int rounds;
	
	/** Maximum length of homopolymer sequences allowed (currently unused) */
	static int MAX_HOMOPOLYMER_LENGTH=99;
	
}
