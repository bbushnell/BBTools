package stream;

import java.util.ArrayList;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;

/**
 * Generates synthetic reads by sliding across reference chromosomes in order.
 * Produces controlled coverage for testing, with optional overlap and strand alternation.
 * @author Brian Bushnell
 * @date 2013
 */
public class SequentialReadInputStream extends ReadInputStream {
	
	/**
	 * Creates a sequential read generator with bounds, read length, overlap, and strand options.
	 *
	 * @param maxReads_ Maximum reads to generate (-1 for unlimited)
	 * @param readlen_ Target read length
	 * @param minreadlen_ Minimum length after trimming undefined bases
	 * @param overlap_ Overlap between consecutive reads
	 * @param alternateStrand_ Whether to reverse-complement odd-numbered reads
	 */
	public SequentialReadInputStream(long maxReads_, int readlen_, int minreadlen_, int overlap_, boolean alternateStrand_){
		
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);
		readlen=readlen_;
		minReadlen=minreadlen_;
		POSITION_INCREMENT=readlen;
		overlap=overlap_;
		alternateStrand=alternateStrand_;
		assert(overlap<POSITION_INCREMENT);
		
		maxPosition=Data.chromLengths[1];
		maxChrom=Data.numChroms;
		
		restart();
	}
	
	/**
	 * Resets positions, counters, and buffers to start generating reads from the beginning.
	 */
	@Override
	public void restart(){
		position=0;
		chrom=1;
		generated=0;
		consumed=0;
		next=0;
		buffer=null;
	}

	/** Indicates whether generated reads are paired; always false.
	 * @return false */
	@Override
	public boolean paired() {
		return false;
	}

	/** No-op close for synthetic generation.
	 * @return false (nothing to close) */
	@Override
	public boolean close() {return false;}
	
	/** Reports whether additional synthetic reads are available based on counters and chromosome bounds.
	 * @return true if more reads can be generated */
	@Override
	public boolean hasMore() {
		if(verbose){
			System.out.println("Called hasMore(): "+(id>=maxReads)+", "+(chrom<maxChrom)+", "+(position<=maxPosition)+", "+(buffer==null || next>=BUF_LEN));
			System.out.println(id+", "+maxReads+", "+chrom+", "+maxChrom+", "+position+", "+maxPosition+", "+buffer+", "+next+", "+(buffer==null ? -1 : BUF_LEN));
		}
//		if(buffer==null || next>=buffer.size()){
//			if(tf.isOpen()){
//				fillBuffer();
//			}else{
//				assert(generated>0) : "Was the file empty?";
//			}
//		}
//		return (buffer!=null && next<buffer.size());
		if(id>=maxReads){return false;}
		if(chrom<maxChrom){return true;}
		if(position<=maxPosition){return true;}
		if(buffer==null || next>=buffer.size()){return false;}
		return true;
	}
	
	/** Returns the next buffered batch of synthetic reads, generating as needed.
	 * @return List of generated reads, or null if exhausted */
	@Override
	public synchronized ArrayList<Read> nextList() {
		if(next!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");} //Possible bug: Should use buffer.isEmpty() not buffer.size()
		if(!hasMore()){return null;}
		if(buffer==null || next>=buffer.size()){fillBuffer();}
		ArrayList<Read> r=buffer;
		buffer=null;
		if(r!=null && r.size()==0){r=null;}
		consumed+=(r==null ? 0 : r.size());
		return r;
	}
	
	private synchronized void fillBuffer(){
//		System.out.println("fill "+chrom+", "+position);
		buffer=null;
		if(chrom>maxChrom){return;}
		ChromosomeArray cha=Data.getChromosome(chrom);
		next=0;
		
		if(position==0){
			while(position<=maxPosition && !AminoAcid.isFullyDefined((char)cha.get(position))){position++;} //Skip initial undefined bases at chromosome start //Possible bug: Undefined bases not being properly skipped
		}
		
		ArrayList<Read> reads=new ArrayList<Read>(BUF_LEN);
		int index=0;
		
		while(position<=maxPosition && index<buffer.size() && id<maxReads){
			int start=position;
			int stop=Tools.min(position+readlen-1, cha.maxIndex);
			byte[] s=cha.getBytes(start, stop);
//			assert(s.length==readlen) : s.length+", "+readlen;
			
			if(s.length<1 || !AminoAcid.isFullyDefined(s)){
				int firstGood=-1, lastGood=-1;
				for(int i=0; i<s.length; i++){ //Find longest contiguous defined region //Possible bug: Generic syntax may cause type safety issues
					if(AminoAcid.isFullyDefined(s[i])){
						lastGood=i;
						if(firstGood==-1){firstGood=i;}
					}
				}
				if(lastGood-firstGood+1>=minReadlen){
					start=start+firstGood;
					stop=stop-(s.length-lastGood-1);
					s=KillSwitch.copyOfRange(s, firstGood, lastGood+1);
					assert(s.length==lastGood-firstGood+1);
				}else{
					s=null;
				}
			}
			
			if(s!=null){
				Read r=new Read(s, null, id, chrom, start, stop, Shared.PLUS);
				if(alternateStrand && (r.numericID&1)==1){r.reverseComplement();}
				r.setSynthetic(true);
//				System.out.println("Made read: "+r);
//				assert(id!=54406) : "\n"+r.toString()+"\nbases: "+s.length+"\nstart: "+start+"\nstop: "+stop+"\nminlen: "+minReadlen+"\n";
				
				reads.add(r);
				index++;
				position+=(POSITION_INCREMENT-overlap);
				id++;
			}else{
				//Move to the next defined position
				while(AminoAcid.isFullyDefined((char)cha.get(position))){position++;}
				while(position<=maxPosition && !AminoAcid.isFullyDefined((char)cha.get(position))){position++;}
			}
		}
//		System.out.println("got "+index+" from "+chrom+", "+position);
		
		if(index==0){
			if(UNLOAD && chrom>0){Data.unload(chrom, true);}
			chrom++;
			position=0;
			buffer=null;
			fillBuffer();
			return;
		}
		
		generated+=index;
		
		buffer=reads;
	}
	
	/** Returns an identifier for this synthetic stream ("sequential").
	 * @return Stream name */
	@Override
	public String fname(){return "sequential";}
	
	private long id=0;
	
	public int position=0;
	public int maxPosition;
	
	private int chrom;
	
	private ArrayList<Read> buffer=null;
	private int next=0;
	
	private final int BUF_LEN=Shared.bufferLen();;
	public static boolean UNLOAD=false;

	public long generated=0;
	public long consumed=0;
	
	public final long maxReads;
	public final int readlen;
	public final int POSITION_INCREMENT;
	public final int minReadlen;
	public final int maxChrom;
	public final int overlap;
	public final boolean alternateStrand;
	
	public static boolean verbose=false;
	
}
