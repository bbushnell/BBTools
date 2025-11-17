package stream;

import java.util.ArrayList;

import dna.Data;
import shared.Shared;
import shared.Tools;
import synth.RandomReads3;

/**
 * @author Brian Bushnell
 * @date Sep 10, 2014
 *
 */
public class RandomReadInputStream3 extends ReadInputStream {
	
	/**
	 * Creates a random read stream with default parameters.
	 * Sets genome build, read count, and pairing mode while using standard
	 * error rates and quality score ranges.
	 *
	 * @param number_ Total number of reads to generate
	 * @param paired_ Whether to generate paired-end reads
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
	 * Creates a random read stream with fully customizable parameters.
	 * Allows precise control over read lengths, error types and rates,
	 * chromosome range, quality scores, and pairing mode.
	 *
	 * @param number_ Total number of reads to generate
	 * @param minreadlen_ Minimum read length in bases
	 * @param maxreadlen_ Maximum read length in bases
	 * @param maxSnps_ Maximum SNPs per read
	 * @param maxInss_ Maximum insertions per read
	 * @param maxDels_ Maximum deletions per read
	 * @param maxSubs_ Maximum substitutions per read
	 * @param snpRate_ Probability of SNPs occurring
	 * @param insRate_ Probability of insertions occurring
	 * @param delRate_ Probability of deletions occurring
	 * @param subRate_ Probability of substitutions occurring
	 * @param maxInsertionLen_ Maximum length of insertion events
	 * @param maxDeletionLen_ Maximum length of deletion events
	 * @param maxSubLen_ Maximum length of substitution events
	 * @param minChrom_ Minimum chromosome number to sample from
	 * @param maxChrom_ Maximum chromosome number to sample from
	 * @param paired_ Whether to generate paired-end reads
	 * @param minQual_ Minimum quality score
	 * @param midQual_ Middle quality score for distribution
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
	
	@Override
	public boolean hasMore() {
		return number>consumed;
	}
	
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
	
	/**
	 * Fills the internal buffer with randomly generated reads.
	 * Calculates remaining reads needed and delegates generation to RandomReads3.
	 * Buffer size is limited to BUF_LEN to control memory usage.
	 */
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
	
	@Override
	public synchronized void restart(){
		next=0;
		buffer=null;
		consumed=0;
		generated=0;
		rr=new RandomReads3(1, paired);
	}

	@Override
	public boolean close() {return false;}

	@Override
	public boolean paired() {
		return paired;
	}
	
	@Override
	public String fname(){return "random";}
	
	/** Buffer for storing generated reads before consumption */
	private ArrayList<Read> buffer=null;
	/** Index of next read to return from buffer */
	private int next=0;
	
	/** Buffer size limit from shared configuration */
	private final int BUF_LEN=Shared.bufferLen();;

	/** Total number of reads generated so far */
	public long generated=0;
	/** Total number of reads consumed by caller */
	public long consumed=0;
	
	/** Target total number of reads to generate */
	public long number=100000;
	/** Minimum read length in bases */
	public int minreadlen=100;
	/** Maximum read length in bases */
	public int maxreadlen=100;

	/** Maximum length of insertion error events */
	public int maxInsertionLen=6;
	/** Maximum length of substitution error events */
	public int maxSubLen=6;
	/** Maximum length of deletion error events */
	public int maxDeletionLen=100;
	/** Maximum length of N-base (ambiguous) sequences */
	public int maxNLen=6;

	/** Minimum length of insertion error events */
	public int minInsertionLen=1;
	/** Minimum length of substitution error events */
	public int minSubLen=1;
	/** Minimum length of deletion error events */
	public int minDeletionLen=1;
	/** Minimum length of N-base (ambiguous) sequences */
	public int minNLen=1;
	
	/** Minimum chromosome number to sample reads from */
	public int minChrom=1;
	/** Maximum chromosome number to sample reads from */
	public int maxChrom=22;
	
	/** Maximum number of SNPs per read */
	public int maxSnps=4;
	/** Maximum number of insertions per read */
	public int maxInss=2;
	/** Maximum number of deletions per read */
	public int maxDels=2;
	/** Maximum number of substitutions per read */
	public int maxSubs=2;
	/** Maximum number of N-base regions per read */
	public int maxNs=2;

	/** Probability rate for SNP error generation */
	public float snpRate=0.5f;
	/** Probability rate for insertion error generation */
	public float insRate=0.25f;
	/** Probability rate for deletion error generation */
	public float delRate=0.25f;
	/** Probability rate for substitution error generation */
	public float subRate=0.10f;
	/** Probability rate for N-base (ambiguous) generation */
	public float NRate=0.10f;
	
	/** Whether to generate paired-end reads (true) or single-end (false) */
	public final boolean paired;

	/** Minimum quality score for generated bases */
	public final byte minQual;
	/** Middle quality score used in quality distribution */
	public final byte midQual;
	/** Maximum quality score for generated bases */
	public final byte maxQual;
	
	/** Generator instance for creating random reads */
	private RandomReads3 rr;

}
