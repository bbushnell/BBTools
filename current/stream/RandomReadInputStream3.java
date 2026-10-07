package stream;

import java.util.ArrayList;

import dna.Data;
import shared.Shared;
import shared.Tools;
import synth.RandomReads3;

/**
 * Batch adapter for RandomReads3 using the process-wide configured reference genome.
 * Targets and counters count fragments: each list entry is one read with an optional
 * attached mate. Nonpositive targets are exhausted; negative does not mean unlimited.
 * Construction still initializes genome/generator state for an empty target.
 * Batch generation and restart are synchronized, but public configuration, counters
 * and status queries require caller ownership. Configure between uses, not concurrently.
 * Generation also uses shared settings and updates FASTQ.TAG_CUSTOM through RandomReads3.
 * Returned lists and reads are handed to the caller without copying.
 * @author Brian Bushnell
 * @date Sep 10, 2014
 */
public class RandomReadInputStream3 extends ReadInputStream{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Initializes Data.GENOME_BUILD and creates a seed-1 generator with default settings.
	 * Requires the currently configured positive genome build. Uses Data.numChroms as
	 * the maximum chromosome and quality parameters 6, 18 and 30.
	 * @param number_ Finite fragment target; nonpositive values produce no batches
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
	 * Initializes the configured genome and a seed-1 generator with supplied parameters.
	 * N-event controls retain their field defaults. Quality integers are narrowed to bytes;
	 * generator preconditions still apply. This class does not validate every parameter.
	 *
	 * @param number_ Finite fragment target; nonpositive values produce no batches
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
	 * @param minQual_ Minimum quality parameter, narrowed to byte
	 * @param midQual_ Midpoint quality parameter, narrowed to byte
	 * @param maxQual_ Maximum quality parameter, narrowed to byte
	 */
	public RandomReadInputStream3(long number_, int minreadlen_, int maxreadlen_,
			int maxSnps_, int maxInss_, int maxDels_, int maxSubs_,
			float snpRate_, float insRate_, float delRate_, float subRate_,
			int maxInsertionLen_, int maxDeletionLen_, int maxSubLen_,
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

		minQual=(byte)minQual_;
		midQual=(byte)midQual_;
		maxQual=(byte)maxQual_;

		restart();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Compares the fragment target to consumed entries without generating or reading ahead. */
	@Override
	public boolean hasMore(){return number>consumed;}

	/** Generates and transfers the next batch using current generation parameters.
	 * List entries are fragments; attached mates are not counted separately. Releases the
	 * internal list reference and increments consumed by the returned entry count.
	 * @return Caller-owned list and Read references, or null at the target or for an empty batch
	 */
	@Override
	public synchronized ArrayList<Read> nextList(){
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

	/** Requests up to the captured batch size from the remaining fragment quota.
	 * Clears the buffer/index first and counts the generator's returned list entries.
	 * Configuration fields are forwarded on every fill; the generator may use shared state.
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
		buffer=reads;//buffer IS assigned the batch here (cf. SequentialReadInputStream#001, where the
		//analogous fillBuffer left buffer null and self-referenced buffer.size() -> NPE). Keep this assignment.
//		assert(false) : reads.size()+", "+toMake;
	}

	/** Clears buffered state and counters, then creates a new generator with seed 1.
	 * Retains this adapter's configuration and inherited errorState; does not reload the
	 * genome or restore generator-wide settings. Fixed seed is not a promise of replay
	 * across changed global settings, genome state or batching configuration.
	 */
	@Override
	public synchronized void restart(){
		next=0;
		buffer=null;
		consumed=0;
		generated=0;
		//TODO: Probable bug in RandomReads3.fillRandomChrom: total reference length below 8192
		//makes its total/8192 divisor zero. Small-genome reachability is unverified; generator review needed.
		rr=new RandomReads3(1, paired);
	}

	/** No-op that neither stops generation nor clears state or the inherited error flag.
	 * @return false, independently of the inherited error state
	 */
	@Override
	public boolean close(){return false;}

	/** Returns the fixed pairing mode passed to each newly created generator. */
	@Override
	public boolean paired(){return paired;}

	/** Returns an identifier for this synthetic stream ("random").
	 * @return Stream name */
	@Override
	public String fname(){return "random";}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Pending list, released to the caller by nextList without copying. */
	private ArrayList<Read> buffer=null;
	/** Legacy element index; currently only reset to zero and checked before batch access. */
	private int next=0;

	/** Maximum requested fragments per batch, captured during instance initialization. */
	private final int BUF_LEN=Shared.bufferLen();

	/** Fragments returned by the generator since restart. */
	public long generated=0;
	/** Fragments handed to callers since restart. */
	public long consumed=0;

	/** Finite target in fragments; nonpositive values are exhausted after construction. */
	public long number=100000;
	/** Minimum read length forwarded to the generator. */
	public int minreadlen=100;
	/** Maximum read length forwarded to the generator. */
	public int maxreadlen=100;

	/** Maximum insertion length forwarded for each fill. */
	public int maxInsertionLen=6;
	/** Maximum substitution-run length forwarded for each fill. */
	public int maxSubLen=6;
	/** Maximum deletion length forwarded for each fill. */
	public int maxDeletionLen=100;
	/** Maximum N-run length forwarded for each fill. */
	public int maxNLen=6;

	/** Minimum insertion length forwarded for each fill. */
	public int minInsertionLen=1;
	/** Minimum substitution-run length forwarded for each fill. */
	public int minSubLen=1;
	/** Minimum deletion length forwarded for each fill. */
	public int minDeletionLen=1;
	/** Minimum N-run length forwarded for each fill. */
	public int minNLen=1;

	/** Minimum chromosome number forwarded to the generator. */
	public int minChrom=1;
	/** Maximum chromosome number captured or supplied at construction; remains mutable. */
	public int maxChrom=22;

	/** Maximum SNP-event count forwarded to the generator. */
	public int maxSnps=4;
	/** Maximum insertion-event count forwarded to the generator. */
	public int maxInss=2;
	/** Maximum deletion-event count forwarded to the generator. */
	public int maxDels=2;
	/** Maximum substitution-run count forwarded to the generator. */
	public int maxSubs=2;
	/** Maximum N-run count forwarded to the generator. */
	public int maxNs=2;

	/** SNP probability parameter forwarded to the generator. */
	public float snpRate=0.5f;
	/** Insertion probability parameter forwarded to the generator. */
	public float insRate=0.25f;
	/** Deletion probability parameter forwarded to the generator. */
	public float delRate=0.25f;
	/** Substitution-run probability parameter forwarded to the generator. */
	public float subRate=0.10f;
	/** N-run probability parameter forwarded to the generator. */
	public float NRate=0.10f;

	/** Fixed pairing mode; mates are attached to returned primary entries. */
	public final boolean paired;

	/** Retained minimum quality parameter. */
	public final byte minQual;
	/** Retained midpoint quality parameter. */
	public final byte midQual;
	/** Retained maximum quality parameter. */
	public final byte maxQual;

	/** Generator replaced by restart with seed 1; uses shared genome/settings. */
	private RandomReads3 rr;

}
