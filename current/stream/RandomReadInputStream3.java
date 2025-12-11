package stream;

import java.util.ArrayList;

import dna.Data;
import shared.Shared;
import shared.Tools;
import synth.RandomReads3;

/**
 * Generates random synthetic reads for testing using RandomReads3 with configurable errors, lengths, and pairing.
 * @author Brian Bushnell
 * @date Sep 10, 2014
 */
public class RandomReadInputStream3 extends ReadInputStream {
	
	/**
	 * Creates a random read stream with default error/quality settings.
	 * @param number_ Total reads to generate
	 * @param paired_ True to generate paired-end reads
	 */
	public RandomReadInputStream3(long number_, boolean paired_){
		Data.setGenome(Data.GENOME_BUILD);
		number=number_;
		paired=paired_;
		maxChrom=Data.numChroms;
		minQual=6;
		midQual=18;
		maxQual=30;
		restart();
	}
	
	/**
	 * Creates a fully parameterized random read generator.
	 * Allows control over read lengths, error counts/rates, chrom range, quality scores, and pairing.
	 *
	 * @param number_ Total reads to generate
	 * @param minreadlen_ Minimum read length
	 * @param maxreadlen_ Maximum read length
	 * @param maxSnps_ Max SNPs per read
	 * @param maxInss_ Max insertions per read
	 * @param maxDels_ Max deletions per read
	 * @param maxSubs_ Max substitutions per read
	 * @param snpRate_ SNP probability
	 * @param insRate_ Insertion probability
	 * @param delRate_ Deletion probability
	 * @param subRate_ Substitution probability
	 * @param maxInsertionLen_ Max insertion length
	 * @param maxDeletionLen_ Max deletion length
	 * @param maxSubLen_ Max substitution length
	 * @param minChrom_ Minimum chromosome to sample
	 * @param maxChrom_ Maximum chromosome to sample
	 * @param paired_ Generate paired-end reads
	 * @param minQual_ Minimum quality score
	 * @param midQual_ Midpoint quality score
	 * @param maxQual_ Maximum quality score
	 */
	public RandomReadInputStream3(long number_, int minreadlen_,  int maxreadlen_,
			int maxSnps_, int maxInss_, int maxDels_, int maxSubs_,
			float snpRate_, float insRate_, float delRate_, float subRate_,
			int maxInsertionLen_, int maxDeletionLen_,  int maxSubLen_,
			int minChrom_, int maxChrom_, boolean paired_,
			int minQual_, int midQual_, int maxQual_){
		Data.setGenome(Data.GENOME_BUILD);
		number=number_;
		minreadlen=minreadlen_;
		maxreadlen=maxreadlen_;

		maxInsertionLen=maxInsertionLen_;
		maxSubLen=maxSubLen_;
		maxDeletionLen=maxDeletionLen_;


		minInsertionLen=1;
		minSubLen=1;
		minDeletionLen=1;
		minNLen=1;
		
		minChrom=minChrom_;
		maxChrom=maxChrom_;
		
		maxSnps=maxSnps_;
		maxInss=maxInss_;
		maxDels=maxDels_;
		maxSubs=maxSubs_;

		snpRate=snpRate_;
		insRate=insRate_;
		delRate=delRate_;
		subRate=subRate_;
		
		paired=paired_;
		
		minQual=(byte) minQual_;
		midQual=(byte) midQual_;
		maxQual=(byte) maxQual_;
		
		restart();
	}
	
	/** Returns true while additional random reads remain to be generated. */
	@Override
	public boolean hasMore() {
		return number>consumed;
	}
	
	/** Generates and returns the next buffered batch of random reads.
	 * @return List of generated reads, or null if exhausted */
	@Override
	public synchronized ArrayList<Read> nextList() {
		if(next!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");}
		if(consumed>=number){return null;}
		if(buffer==null || next>=buffer.size()){fillBuffer();}
		ArrayList<Read> r=buffer;
		buffer=null;
		if(r!=null && r.size()==0){r=null;}
		consumed+=(r==null ? 0 : r.size());
//		assert(false) : r.size();
		return r;
	}
	
	private synchronized void fillBuffer(){
		buffer=null;
		next=0;
		
		long toMake=number-generated;
		if(toMake<1){return;}
		toMake=Tools.min(toMake, BUF_LEN);
		
		ArrayList<Read> reads=rr.makeRandomReadsX((int)toMake, minreadlen, maxreadlen, -1,
				maxSnps, maxInss, maxDels, maxSubs, maxNs,
				snpRate, insRate, delRate, subRate, NRate,
				minInsertionLen, minDeletionLen, minSubLen, minNLen,
				maxInsertionLen, maxDeletionLen, maxSubLen, maxNLen,
				minChrom, maxChrom,
				minQual, midQual, maxQual);
		
		generated+=reads.size();
		assert(generated<=number);
		buffer=reads;
//		assert(false) : reads.size()+", "+toMake;
	}
	
	/**
	 * Resets counters and regenerates the RandomReads3 generator for a fresh run.
	 */
	@Override
	public synchronized void restart(){
		next=0;
		buffer=null;
		consumed=0;
		generated=0;
		rr=new RandomReads3(1, paired);
	}

	/** No-op close for synthetic generation.
	 * @return false (nothing to close) */
	@Override
	public boolean close() {return false;}

	@Override
	public boolean paired() {
		return paired;
	}
	
	/** Returns an identifier for this synthetic stream ("random").
	 * @return Stream name */
	@Override
	public String fname(){return "random";}
	
	private ArrayList<Read> buffer=null;
	private int next=0;
	
	private final int BUF_LEN=Shared.bufferLen();;

	public long generated=0;
	public long consumed=0;
	
	public long number=100000;
	public int minreadlen=100;
	public int maxreadlen=100;

	public int maxInsertionLen=6;
	public int maxSubLen=6;
	public int maxDeletionLen=100;
	public int maxNLen=6;

	public int minInsertionLen=1;
	public int minSubLen=1;
	public int minDeletionLen=1;
	public int minNLen=1;
	
	public int minChrom=1;
	public int maxChrom=22;
	
	public int maxSnps=4;
	public int maxInss=2;
	public int maxDels=2;
	public int maxSubs=2;
	public int maxNs=2;

	public float snpRate=0.5f;
	public float insRate=0.25f;
	public float delRate=0.25f;
	public float subRate=0.10f;
	public float NRate=0.10f;
	
	public final boolean paired;

	public final byte minQual;
	public final byte midQual;
	public final byte maxQual;
	
	private RandomReads3 rr;

}
