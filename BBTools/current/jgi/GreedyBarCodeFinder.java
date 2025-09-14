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
	
	public static void main(String[] args){
		Timer t=new Timer();
		
		GreedyBarCodeFinder finder=new GreedyBarCodeFinder(args);
		int best=finder.find(finder.rounds);
		
		t.stop();
		System.err.println("There are at least "+best+" codes of length "+finder.k+" with mutual hamming distance at least "+finder.hdist);
		System.err.println("Time: \t"+t);
	}
	
	public GreedyBarCodeFinder(String[] args){
		k=Integer.parseInt(args[0]);
		hdist=Integer.parseInt(args[1]);
		rounds=(args.length>2 ? Integer.parseInt(args[2]) : 20);
	}
	
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
	
	private final int k;
	private final int hdist;
	private int rounds;
	
	static int MAX_HOMOPOLYMER_LENGTH=99;
	
}
