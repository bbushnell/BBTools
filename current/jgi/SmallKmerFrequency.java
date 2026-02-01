package jgi;

import java.util.Arrays;
import java.util.Comparator;

import dna.AminoAcid;
import fileIO.FileFormat;
import parse.Parse;
import shared.Timer;
import shared.Tools;
import stream.Read;
import template.BBTool_ST;

/**
 * Generates small-k k-mer frequency profiles for reads and writes the top hits.
 * Builds canonical k-mer indices, counts occurrences, and formats results per read.
 * @author Brian Bushnell
 * @date Feb 19, 2015
 */
public class SmallKmerFrequency extends BBTool_ST {
	
	/** Entry point that constructs the tool and runs processing with timing.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		FileFormat.PRINT_WARNING=false;
		SmallKmerFrequency bbt=new SmallKmerFrequency(args);
		bbt.process(t);
	}
	
	/** Sets default parameters: k=2, display=3, numeric counts disabled. */
	@Override
	protected void setDefaults(){
		k=2;
		display=3;
		addNumbers=false;
	}

	/** Parses arguments, builds k-mer index and count arrays, and prepares output.
	 * @param args Command-line arguments */
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
	/**
	 * Parses tool-specific arguments (k, display, addnumbers/number/counts).
	 *
	 * @param arg Full argument string
	 * @param a Argument key
	 * @param b Argument value
	 * @return true if the argument was handled
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
	
	/**
	 * Creates k-mer frequency profiles for each read in a pair and attaches results.
	 * Counts canonical k-mers, sorts by frequency, and emits top entries with optional counts.
	 *
	 * @param r1 First read (may be null)
	 * @param r2 Second read (may be null)
	 * @return true on success
	 */
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
			for(int i=0; i<display; i++){ //Possible bug: should check kmers[i].count>0 like r1
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
	
	/**
	 * Generates a canonical k-mer count array for a sequence.
	 * Resets counts on ambiguous bases and uses a rolling hash for efficiency.
	 *
	 * @param bases Input bases
	 * @param array_ Reusable count array (allocated if null)
	 * @param clear Whether to zero the array before counting
	 * @return The populated count array
	 */
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
	
	/** Subclass hook for startup (no-op). */
	@Override
	protected void startupSubclass() {}
	
	/** Subclass hook for shutdown (no-op). */
	@Override
	protected void shutdownSubclass() {}
	
	/**
	 * Subclass hook for stats reporting (no-op).
	 * @param t Timer tracking execution
	 * @param readsIn Reads processed
	 * @param basesIn Bases processed
	 */
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
	
	/**
	 * Builds canonical index mapping for all k-mers of length n, pairing forward and reverse complements.
	 * @param n K-mer length
	 * @return Mapping from k-mer value to canonical index
	 */
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
	
	/** Indicates this tool uses the shared header format.
	 * @return true */
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
