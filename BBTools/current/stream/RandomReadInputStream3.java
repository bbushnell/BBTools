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
	public void start() {}
	
	
	@Override
	public boolean hasMore() {
		return number>consumed;
	}

	@Override
	public Read next() {
		if(consumed>=number){return null;}
		if(buffer==null || next>=buffer.size()){fillBuffer();}
		Read r=buffer.get(next);
		buffer.set(next, null);
		next++;
		consumed++;
		return r;
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
