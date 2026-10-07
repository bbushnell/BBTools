package stream;

import java.util.ArrayList;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;

/**
 * Generates unpaired synthetic reads by sliding across loaded reference chromosomes.
 * Reference metadata in Data must be initialized before construction; chromosomes are
 * obtained through Data.getChromosome in numeric order, starting at chromosome 1.
 * Reads have numeric names, reference coordinates, no qualities and the synthetic flag.
 * Optional reverse complementation changes bases and strand, not reference coordinates.
 * Callers must coordinate access to the reader and shared reference/configuration state;
 * synchronized generation does not synchronize restart, configuration or external reads
 * of public counters.
 * @author Brian Bushnell
 * @date 2013
 */
public class SequentialReadInputStream extends ReadInputStream{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a generator using the current Data reference metadata and Shared buffer size.
	 *
	 * @param maxReads_ Maximum individual reads per pass; any negative value is unlimited
	 * @param readlen_ Target window length and initial position increment
	 * @param minreadlen_ Minimum retained span when a candidate contains undefined bases;
	 * fully defined candidates are retained without this minimum-length check
	 * @param overlap_ Amount subtracted from the position increment after an emitted read;
	 * must be less than readlen_, with negative values leaving gaps between windows
	 * @param alternateStrand_ Whether to reverse-complement reads with odd numeric IDs
	 */
	public SequentialReadInputStream(long maxReads_, int readlen_, int minreadlen_, int overlap_, boolean alternateStrand_){
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);
		readlen=readlen_;
		minReadlen=minreadlen_;
		POSITION_INCREMENT=readlen;
		overlap=overlap_;
		alternateStrand=alternateStrand_;
		assert(overlap<POSITION_INCREMENT);
		
		//[stream/SequentialReadInputStream#004] Initialize the first bound; parsing refreshes it at each reference entry.
		maxPosition=Data.chromLengths[1];
		maxChrom=Data.numChroms;
		
		restart();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Resets IDs, strand alternation, positions, counters and buffers for a new pass.
	 * Starts at chromosome 1 using the existing metadata/configuration; the reference
	 * bound is refreshed on entry. Does not clear the inherited error flag or preload
	 * chromosomes. An unloaded chromosome is obtained again through Data on demand.
	 */
	@Override
	public void restart(){
		//[stream/SequentialReadInputStream#003] Reset identity so read limits and strand alternation start over with each pass.
		id=0;
		position=0;
		chrom=1;
		generated=0;
		consumed=0;
		next=0;
		buffer=null;
	}

	/** @return false; generated records have no mates */
	@Override
	public boolean paired(){return false;}

	/**
	 * Performs no cleanup and does not prevent later generation or unload references.
	 * @return false; this method does not report the inherited error flag
	 */
	@Override
	public boolean close(){return false;}
	
	/**
	 * Checks the read limit, reference position and existing buffer without generating.
	 * A true result may precede a null batch if remaining references yield no usable reads.
	 * @return Whether more records may be available; false after known exhaustion
	 */
	@Override
	public boolean hasMore(){
		if(verbose){
			System.out.println("Called hasMore(): "+(id>=maxReads)+", "+(chrom<maxChrom)+", "+(position<=maxPosition)+", "+(buffer==null || next>=BUF_LEN));
			System.out.println(id+", "+maxReads+", "+chrom+", "+maxChrom+", "+position+", "+maxPosition+", "+buffer+", "+next+", "+(buffer==null ? -1 : BUF_LEN));
		}
		if(id>=maxReads){return false;}
		if(chrom<maxChrom){return true;}
		//[stream/SequentialReadInputStream#008] Terminal advancement resets position to zero; require a valid chromosome before advertising more.
		if(chrom<=maxChrom && position<=maxPosition){return true;}
		if(buffer==null || next>=buffer.size()){return false;}
		return true;
	}
	
	/**
	 * Generates and hands off the next batch, updating the consumed-record count.
	 * Each returned list and its records remain owned by the caller; the reader does not
	 * reuse them. A batch contains at most the captured buffer size and belongs to one
	 * reference chromosome. Empty chromosomes may be skipped before producing a batch.
	 * @return Nonempty list of unpaired synthetic reads, or null when exhausted
	 */
	@Override
	public synchronized ArrayList<Read> nextList(){
		if(next!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");}
		if(!hasMore()){return null;}
		if(buffer==null || next>=buffer.size()){fillBuffer();}
		ArrayList<Read> r=buffer;
		buffer=null;
		if(r!=null && r.size()==0){r=null;}
		consumed+=(r==null ? 0 : r.size());
		return r;
	}
	
	/**
	 * Fills one reference-local batch, skipping initial undefined bases and unsuitable
	 * windows. Undefined candidate ends are trimmed to the first and last defined bases;
	 * undefined bases inside that span remain. Iterates across exhausted references until
	 * a batch is ready, optionally unloading cached chromosome arrays as it advances.
	 */
	private synchronized void fillBuffer(){
		buffer=null;
		while(chrom<=maxChrom){
			ChromosomeArray cha=Data.getChromosome(chrom);
			next=0;

			if(position==0){
				//[stream/SequentialReadInputStream#004] Unequal references need their own bound, including re-entry after restart.
				maxPosition=Data.chromLengths[chrom];
				while(position<=maxPosition && !AminoAcid.isFullyDefined((char)cha.get(position))){position++;}//Skip initial undefined bases
			}

			ArrayList<Read> reads=new ArrayList<Read>(BUF_LEN);
			int index=0;

			//[stream/SequentialReadInputStream#001] cap the batch at BUF_LEN. Was index<buffer.size() — but
			//buffer is set null at the top of fillBuffer() and never reassigned before here, so buffer.size()
			//was a GUARANTEED NPE on any non-trivial chromosome (the local 'reads' is the batch being built,
			//sized BUF_LEN). NEEDS end-to-end validation with a real reference (synthetic-read mode) before trusting.
			while(position<=maxPosition && index<BUF_LEN && id<maxReads){
				int start=position;
				int stop=Tools.min(position+readlen-1, cha.maxIndex);
				byte[] s=cha.getBytes(start, stop);

				if(s.length<1 || !AminoAcid.isFullyDefined(s)){
					int firstGood=-1, lastGood=-1;
					//[stream/SequentialReadInputStream#005] Find the first-to-last defined span; retain interior undefined bases.
					for(int i=0; i<s.length; i++){
						if(AminoAcid.isFullyDefined(s[i])){
							lastGood=i;
							if(firstGood==-1){firstGood=i;}
						}
					}
					//TODO: Probable bug [stream/SequentialReadInputStream#006] - no defined bases gives span 1 from -1/-1; minReadlen<=1 can copy from -1. Known production callers use at least 50.
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

					reads.add(r);
					index++;
					position+=(POSITION_INCREMENT-overlap);
					id++;
				}else{
					//Move to the next defined position
					//ChromosomeArray.get returns N beyond maxIndex, so the first scan stops without an out-of-range array access.
					while(AminoAcid.isFullyDefined((char)cha.get(position))){position++;}
					while(position<=maxPosition && !AminoAcid.isFullyDefined((char)cha.get(position))){position++;}
				}
			}

			if(index==0){
				//Resolved #007: iterate to the next reference without retaining a call frame for each skip.
				if(UNLOAD && chrom>0){Data.unload(chrom, true);}
				chrom++;
				position=0;
				buffer=null;
				continue;
			}

			generated+=index;

			buffer=reads;
			return;
		}
	}
	
	/** Returns an identifier for this synthetic stream ("sequential").
	 * @return Stream name */
	@Override
	public String fname(){return "sequential";}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Next numeric read ID; also counts records constructed during the pass toward its limit. */
	private long id=0;
	
	/** Current reference position; zero marks entry and triggers the reference-bound refresh. */
	public int position=0;
	/** Reference-length bound refreshed on reference entry; manual pre-entry assignments are overwritten. */
	public int maxPosition;
	
	/** Current 1-based chromosome number; greater than maxChrom after terminal advancement. */
	private int chrom;
	
	/** Current generated batch until handed off; discarded on restart. */
	private ArrayList<Read> buffer=null;
	/** Legacy buffer cursor; blockwise access requires zero and does not advance it. */
	private int next=0;
	
	/** Batch entry capacity captured from Shared when this instance is initialized. */
	private final int BUF_LEN=Shared.bufferLen();
	/** Whether advancing past a reference unloads its cached chromosome; read during generation. */
	public static boolean UNLOAD=false;

	/** Individual reads placed in completed batches during this pass. */
	public long generated=0;
	/** Individual reads handed off by nextList during this pass. */
	public long consumed=0;
	
	/** Per-pass individual-read limit; negative constructor values become Long.MAX_VALUE. */
	public final long maxReads;
	/** Target candidate window length before reference-end clipping or undefined-end trimming. */
	public final int readlen;
	/** Position increment before subtracting overlap; equal to readlen. */
	public final int POSITION_INCREMENT;
	/** Minimum trimmed span for candidates containing undefined bases, not a universal output minimum. */
	public final int minReadlen;
	/** Last chromosome number captured from Data.numChroms at construction. */
	public final int maxChrom;
	/** Amount subtracted from POSITION_INCREMENT after each emitted read. */
	public final int overlap;
	/** Whether reads with odd numeric IDs are reverse-complemented. */
	public final boolean alternateStrand;
	
	/** Enables diagnostic availability messages on standard output. */
	public static boolean verbose=false;
	
}
