package jgi;

import java.util.Arrays;
import java.util.Comparator;

import dna.AminoAcid;
import fileIO.FileFormat;
import shared.Parse;
import shared.Timer;
import shared.Tools;
import stream.Read;
import template.BBTool_ST;

/**
 * @author Brian Bushnell
 * @date Feb 19, 2015
 *
 */
public class SmallKmerFrequency extends BBTool_ST {
	
	/**
	 * Code entrance from the command line.
	 * Must be overridden; the commented body is an example.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		FileFormat.PRINT_WARNING=false;
		SmallKmerFrequency bbt=new SmallKmerFrequency(args);
		bbt.process(t);
	}
	
	@Override
	protected void setDefaults(){
		k=2;
		display=3;
		addNumbers=false;
	}

	/**
	 * @param args
	 */
	public SmallKmerFrequency(String[] args) {
		super(args);
		reparse(args);
		
		kmerIndex=makeKmerIndex(k);
		maxKmer=Tools.max(kmerIndex);
		counts=new int[maxKmer+1];
		display=Tools.min(display, counts.length);
		if(out1!=null){
			ffout1=FileFormat.testOutput(out1, FileFormat.ATTACHMENT, ".info", true, overwrite, append, false);
		}
		kmers=new Kmer[counts.length];
		for(int i=0; i<kmerIndex.length; i++){
			int index=kmerIndex[i];
			if(kmers[index]==null){
				kmers[index]=new Kmer();
				kmers[index].s=AminoAcid.kmerToString(i, k);
				kmers[index].num=i;
			}
		}
//		System.err.println(Arrays.toString(kmers));
	}

	/* (non-Javadoc)
	 * @see jgi.BBTool_ST#parseArgument(java.lang.String, java.lang.String, java.lang.String)
	 */
	@Override
	public boolean parseArgument(String arg, String a, String b) {
		if(a.equals("k")){
			k=Integer.parseInt(b);
			return true;
		}else if(a.equals("display")){
			display=Integer.parseInt(b);
			return true;
		}else if(a.equals("addnumbers") || a.equals("number") || a.equals("count") || a.equals("numbers") || a.equals("counts")){
			addNumbers=Parse.parseBoolean(b);
			return true;
		}
		return false;
	}
	
	@Override
	protected boolean processReadPair(Read r1, Read r2) {
		if(r1!=null){
			makeKmerProfile(r1.bases, counts, true);
			sb.append(r1.id);
			Arrays.sort(kmers, numComparator);
			for(int i=0; i<counts.length; i++){
				kmers[i].count=counts[i];
			}
			Arrays.sort(kmers, countComparator);
			for(int i=0; i<display && kmers[i].count>0; i++){
				sb.append('\t');
				sb.append(kmers[i].s);
				if(addNumbers){sb.append('=').append(kmers[i].count);}
			}
//			sb.append('\n');
			r1.obj=sb.toString();
			sb.setLength(0);
		}
		if(r2!=null){
			makeKmerProfile(r2.bases, counts, true);
			sb.append(r2.id);
			Arrays.sort(kmers, numComparator);
			for(int i=0; i<counts.length; i++){
				kmers[i].count=counts[i];
			}
			Arrays.sort(kmers, countComparator);
			for(int i=0; i<display; i++){
				sb.append('\t');
				sb.append(kmers[i].s);
				if(addNumbers){sb.append('=').append(kmers[i].count);}
			}
//			sb.append('\n');
			r2.obj=sb.toString();
			sb.setLength(0);
		}
		return true;
	}
	
	/** Makes a kmer (e.g., tetramer) profile of a cluster */
	private final int[] makeKmerProfile(byte[] bases, int[] array_, boolean clear){
		final int nbits=2*k;
		final int[] array=(array_==null ? new int[maxKmer+1] : array_);
		final int mask=~((-1)<<(nbits));
		if(clear){Arrays.fill(array, 0);} //TODO: Can be cleared faster using an IntList.
		
		int keysCounted=0;
		
		int len=0;
		int kmer=0;
		for(byte b : bases){
			int x=AminoAcid.baseToNumber[b];
			if(x<0){
				len=0;
				kmer=0;
			}else{
				kmer=((kmer<<2)|x)&mask;
				len++;
				if(len>=k){
					int rkmer=AminoAcid.reverseComplementBinaryFast(kmer, k);
					keysCounted++;
					array[kmerIndex[Tools.min(kmer, rkmer)]]++;
				}
			}
		}
		return array;
	}
	
	@Override
	protected void startupSubclass() {}
	
	@Override
	protected void shutdownSubclass() {}
	
	@Override
	protected void showStatsSubclass(Timer t, long readsIn, long basesIn) {}
	
	private class Kmer{
		
		String s;
		int count=0;
		int num;
		
		@Override
		public String toString(){return "("+s+","+num+","+count+")";}
		
	}
	
	private static class NumComparator implements Comparator<Kmer>{
		
		@Override
		public int compare(Kmer a, Kmer b) {
			return a.num-b.num;
		}
		
	}
	
	private static class CountComparator implements Comparator<Kmer>{
		
		@Override
		public int compare(Kmer a, Kmer b) {
			return b.count-a.count;
		}
		
	}
	
	public static final int[] makeKmerIndex(final int n){
		final int max=(1<<(2*n))-1;
		int[] array=new int[max+1];
		
		int count=0;
		for(int i=0; i<=max; i++){
			final int a=i, b=AminoAcid.reverseComplementBinaryFast(i, n);
			int min=Tools.min(a, b);
			if(min==a){
				array[a]=array[b]=count;
				count++;
			}
		}
//		assert(false) : Arrays.toString(array);
		return array;
	}
	
	@Override
	protected final boolean useSharedHeader(){return true;}

	private static final NumComparator numComparator=new NumComparator();
	private static final CountComparator countComparator=new CountComparator();
	
	private int k;
	private int display;
	private boolean addNumbers;
	private final int maxKmer;
	private final int[] kmerIndex;
	private final int[] counts;
	private final StringBuilder sb=new StringBuilder();
	
	private final Kmer[] kmers;
	
}
