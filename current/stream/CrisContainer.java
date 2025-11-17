package stream;

import java.util.ArrayList;
import java.util.Comparator;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import sort.ReadComparatorClump;
import sort.ReadComparatorTopological5Bit;
import structures.ListNum;

/**
 * Container that wraps a ConcurrentReadInputStream with comparison capabilities.
 * Provides buffered read access with configurable sorting and pre-processing
 * based on the specified comparator type. Supports k-mer generation and clumping.
 * @author Brian Bushnell
 */
public class CrisContainer implements Comparable<CrisContainer> {
	
	/**
	 * Creates a CrisContainer from a filename with specified comparison behavior.
	 * Opens the file as a read stream and configures pre-processing based on comparator type.
	 *
	 * @param fname Input filename to read from
	 * @param comparator_ Comparator for ordering reads
	 * @param allowSubprocess Whether to allow subprocess execution for file access
	 */
	public CrisContainer(String fname, Comparator<Read> comparator_, boolean allowSubprocess){
		comparator=comparator_;
		genKmer=(comparator==ReadComparatorTopological5Bit.comparator);
		clump=(comparator==ReadComparatorClump.comparator);
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, allowSubprocess, true);
		cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ff, null, null, null);
//		System.err.println(genKmer+", "+clump+", "+comparator.getClass());
		cris.start();
		fetch();
	}
	
	/**
	 * Creates a CrisContainer from an existing ConcurrentReadInputStream.
	 * Configures pre-processing based on comparator type without opening new files.
	 * @param cris_ Existing read input stream to wrap
	 * @param comparator_ Comparator for ordering reads
	 */
	public CrisContainer(ConcurrentReadInputStream cris_, Comparator<Read> comparator_){
		comparator=comparator_;
		genKmer=(comparator==ReadComparatorTopological5Bit.comparator);
		clump=(comparator==ReadComparatorClump.comparator);
		cris=cris_;
//		System.err.println(genKmer+", "+clump+", "+comparator.getClass());
		fetch();
	}
	
	/**
	 * Fetches the next batch of reads and returns the previous batch.
	 * Applies k-mer generation or clumping preprocessing as configured.
	 * @return Previous list of reads, or null if this is the first fetch
	 */
	public ArrayList<Read> fetch(){
		final ArrayList<Read> old=list;
		fetchInner();
		return old;
	}
	
	/** Internal method that performs the actual read fetching and preprocessing.
	 * Handles list management, k-mer generation, clumping, and stream coordination. */
	private void fetchInner(){
		ListNum<Read> ln=cris.nextList();
		list=(ln==null ? null : ln.list);
		if(list.size()<1){list=null;}
		else if(genKmer){
			for(Read r : list){ReadComparatorTopological5Bit.genKmer(r);}
		}else if(clump){
			for(Read r : list){ReadComparatorClump.set(r);}
		}
		read=(list==null ? null : list.get(0));
		if(lastNum>=0){cris.returnList(lastNum, list==null);}
		if(ln!=null){lastNum=ln.id;}
		assert((read==null)==(list==null || list.size()==0));
//		if(count>0 && list!=null){
//			for(Read r : list){
//				assert(remainingReads>=0) : remainingReads+", "+count+", "+r.numericID;
//				double remaining=(count-sum);
//				double mult=2*(remaining/remainingReads);
//				sum=sum+randy.nextDouble()*mult;
//				r.rand=sum;
////				System.err.println(r.rand);
//				remainingReads--;
//			}
//		}
	}
	
	/** Closes the underlying read input stream.
	 * @return true if the stream was successfully closed */
	public boolean close(){
		return ReadWrite.closeStream(cris);
	}
	
	/** Returns the current read without consuming it.
	 * @return Current read at the head of the container, or null if empty */
	public Read peek(){return read;}
	
	@Override
	public int compareTo(CrisContainer other) {
		assert(read!=null);
		assert(other.read!=null);
		return comparator.compare(read, other.read);
	}
	
	/**
	 * Compares this container's current read to a specific read.
	 * Uses the configured comparator to determine ordering.
	 * @param other Read to compare against
	 * @return Negative, zero, or positive value indicating relative order
	 */
	public int compareTo(Read other) {
		return comparator.compare(read, other);
	}
	
	/** Checks if there are more reads available in this container.
	 * @return true if there is a current read available, false if exhausted */
	public boolean hasMore(){
		return read!=null;
	}
	
	/** Returns the underlying ConcurrentReadInputStream.
	 * @return The wrapped read input stream */
	public ConcurrentReadInputStream cris(){return cris;}
	
	/** The underlying concurrent read input stream */
	final ConcurrentReadInputStream cris;
	/** Current read at the head of the container */
	private Read read;
	/** ID number of the last list returned to the stream */
	private long lastNum=-1;
	/** Current list of reads from the input stream */
	private ArrayList<Read> list;
	/** Comparator used for ordering reads */
	private final Comparator<Read> comparator;
	/** Whether to apply clumping preprocessing for reads */
	/** Whether to generate k-mers for reads using topological comparator */
	private final boolean genKmer, clump;
//	private double sum=0;
//	final int count;
//	private final Random randy;
//	private int remainingReads;
	
}
