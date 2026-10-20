package clump;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import sort.SortByName;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ListNum;

/**
 * Drains a read stream into borrowed, already-started output streams. Owns input
 * closure during processing; output closure remains the caller's responsibility.
 * Supports optional name sorting through a temporary FASTQ file and distribution
 * by pivot. Instances are single-use and require exclusive input consumption.
 *
 * @author Brian Bushnell
 * @date January 2025
 */
public class StreamToOutput{
	
	/**
	 * Creates a StreamToOutput from file format specifications.
	 * Creates, retains and starts the input stream immediately.
	 *
	 * @param ffin1 Primary input file format
	 * @param ffin2 Secondary input file format (may be null for single-end)
	 * @param rosa_ Array of output streams for read distribution
	 * @param old K-mer comparator for hashing (may be incremented)
	 * @param sortByName_ Whether to sort reads by name before processing
	 * @param incrementComparator Whether to increment the seed and decrement the border
	 */
	public StreamToOutput(FileFormat ffin1, FileFormat ffin2, ConcurrentReadOutputStream[] rosa_, KmerComparator old, boolean sortByName_, boolean incrementComparator){
		cris=ConcurrentReadInputStream.getReadInputStream(-1, false, ffin1, ffin2, null, null);
		cris.start();
		rosa=rosa_;
		kc=(incrementComparator ? new KmerComparator(old.k, old.seed+1, old.border-1, old.hashes, false, false) : old);
		sortByName=sortByName_;
	}
	
	/**
	 * Creates a StreamToOutput from an existing input stream.
	 * Transfers input consumption/closure to this bridge without starting it.
	 *
	 * @param cris_ Already-started input stream to process
	 * @param rosa_ Array of output streams for read distribution
	 * @param old K-mer comparator for hashing (may be incremented)
	 * @param sortByName_ Whether to sort reads by name before processing
	 * @param incrementComparator Whether to increment the seed and decrement the border
	 */
	public StreamToOutput(ConcurrentReadInputStream cris_, ConcurrentReadOutputStream[] rosa_, KmerComparator old, boolean sortByName_, boolean incrementComparator){
		cris=cris_;
		rosa=rosa_;
		kc=(incrementComparator ? new KmerComparator(old.k, old.seed+1, old.border-1, old.hashes, false, false) : old);
		sortByName=sortByName_;
	}
	
	/**
	 * Main processing method that streams reads from input to outputs.
	 * Optionally sorts reads by name using temporary files before distribution.
	 * Routes to single or multi-stream processing based on output array size.
	 * With no destinations, still consumes/counts reads and closes the input.
	 * @return true if errors occurred during processing, false otherwise
	 */
	public boolean process(){
		//The optional sort spills as needed and writes a temporary FASTQ. SortByName's
		//legacy processInner accepts our input stream but does not close it.
		File temp=null;
		try{
			if(sortByName){
				try{
					temp=File.createTempFile("temp_namesort_", ".fq.gz");
				}catch(IOException e){
					throw new RuntimeException("Could not create a temporary name-sort file", e);
				}
				//stdin.fq supplies required input metadata only: processInner consumes cris,
				//so no stdin reader is created. Preserve its pairing in the temporary file reader.
				SortByName sbn=new SortByName(new String[] {"in=stdin.fq", "out="+temp.getAbsolutePath(),
						"overwrite=t", "interleaved="+cris.paired()});
				sbn.processInner(cris);
				errorState|=sbn.errorState;
				errorState|=ReadWrite.closeStream(cris);
				cris=null;
				if(errorState){return true;}
				FileFormat ff=FileFormat.testInput(temp.getAbsolutePath(), null, false);
				cris=ConcurrentReadInputStream.getReadInputStream(-1, false, ff, null, null, null);
				cris.start();
			}
			if(rosa==null || rosa.length<=1){
				processSingle(cris);
			}else{
				processMulti(cris);
			}
		}finally{
			errorState|=ReadWrite.closeStream(cris);
			if(temp!=null){temp.delete();}
		}
		return errorState;
	}
	
	/**
	 * Drains reads to one output stream, or counts them when no output is supplied.
	 * Forwards whole batches without distribution when an output is present.
	 * Tracks read and base counts during processing.
	 * @param cris Input stream to read from
	 */
	public void processSingle(ConcurrentReadInputStream cris){
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			if(rosa!=null && rosa.length>0){rosa[0].add(reads, ln.id);}
			
			for(Read r : reads){
				readsIn+=r.pairCount();
				basesIn+=r.pairLength();
			}
			
			cris.returnList(ln);
			
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
	}
	
	/**
	 * Processes reads for multiple output stream scenario.
	 * Distributes reads across output streams based on k-mer hash values.
	 * Uses modulo operation to assign reads to groups deterministically.
	 * @param cris Input stream to read from
	 */
	public void processMulti(ConcurrentReadInputStream cris){
		final int groups=rosa.length;
		
		@SuppressWarnings("unchecked")
		final ArrayList<Read>[] out=new ArrayList[groups];
		for(int i=0; i<out.length; i++){
			out[i]=new ArrayList<Read>();
		}
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			for(Read r : reads){
				//Fixed [clump/StreamToOutput#001]: hash the selected pivot just as the normal
				//KmerSort hash-and-split path does. KmerSort3 supplies next-pass comparator
				//selection when repartitioning; raw-pivot modulo used different destinations.
				long kmer=kc.hash(r, null, 0, false);
				int group=(int)(kc.hash(kmer)%groups);
				out[group].add(r);
				
				readsIn+=r.pairCount();
				basesIn+=r.pairLength();
			}
			for(int group=0; group<groups; group++){
				rosa[group].add(out[group], ln.id);
				out[group]=new ArrayList<Read>();
			}
			
			cris.returnList(ln);
			
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
	}
	
	/** Total number of reads transferred, including linked mates. */
	long readsIn=0;
	/** Total number of bases transferred, including linked mates. */
	long basesIn=0;
	
//	final FileFormat ffin1;
//	final FileFormat ffin2;
	/** Input stream for reading sequence data */
	ConcurrentReadInputStream cris;
	/** Array of output streams for distributing processed reads */
	ConcurrentReadOutputStream[] rosa;
	/** K-mer comparator used for hashing reads during multi-stream distribution */
	final KmerComparator kc;
	/** Whether to sort reads by name before processing */
	final boolean sortByName;
	
	/** Tracks whether any errors occurred during processing */
	boolean errorState=false;
	
}
