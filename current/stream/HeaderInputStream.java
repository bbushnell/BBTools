package stream;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.FileFormat;
import shared.Shared;

/**
 * Reads header lines from files, treating each line as a read with only a name field.
 * Extends ReadInputStream to provide streaming access to header-only sequences.
 * Useful for processing files containing just sequence identifiers or metadata.
 *
 * @author Brian Bushnell
 * @date June 1, 2016
 */
public class HeaderInputStream extends ReadInputStream {
	
	public static void main(String[] args){
		
		HeaderInputStream his=new HeaderInputStream(args[0], true);
		
		Read r=his.nextList().get(0);
		System.out.println(r.toText(false));
		his.close();
		
	}
	
	public HeaderInputStream(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.FASTQ, null, allowSubprocess_, false));
	}

	
	public HeaderInputStream(FileFormat ff){
		if(verbose){System.err.println("FastqReadInputStream("+ff+")");}
		
		stdin=ff.stdio();
		
		tf=new ByteFile1(ff);
	}
	
	@Override
	public boolean hasMore() {
		if(buffer==null || next>=buffer.size()){
			if(tf.isOpen()){
				fillBuffer();
			}else{
				assert(generated>0) : "Was the file empty?";
			}
		}
		return (buffer!=null && next<buffer.size());
	}
	
	/**
	 * Returns the entire current buffer as a list of reads.
	 * Used for block-wise processing instead of individual read access.
	 * Clears the buffer after returning and updates consumed count.
	 *
	 * @return ArrayList of Read objects or null if no reads available
	 * @throws RuntimeException if called after using next() method
	 */
	@Override
	public synchronized ArrayList<Read> nextList() {
		if(next!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");}
		if(buffer==null || next>=buffer.size()){fillBuffer();}
		ArrayList<Read> list=buffer;
		buffer=null;
		if(list!=null && list.size()==0){list=null;}
		consumed+=(list==null ? 0 : list.size());
		return list;
	}
	
	private synchronized void fillBuffer(){
		
		assert(buffer==null || next>=buffer.size());
		
		buffer=null;
		next=0;
		
		buffer=toReadList(tf, BUF_LEN, nextReadID);
		int bsize=(buffer==null ? 0 : buffer.size());
		nextReadID+=bsize;
		if(bsize<BUF_LEN){tf.close();}
		
		generated+=bsize;
		if(buffer==null){
			if(!errorState){
				errorState=true;
				System.err.println("Null buffer in FastqReadInputStream.");
			}
		}
	}
	
	/**
	 * Closes the input stream and underlying file.
	 * Reports verbose output if enabled and updates error state.
	 * @return true if any errors occurred during closing, false otherwise
	 */
	@Override
	public boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		errorState|=tf.close();
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		return errorState;
	}

	/** Resets the stream to initial state for re-reading.
	 * Clears counters, buffer, and resets the underlying file reader. */
	@Override
	public synchronized void restart() {
		generated=0;
		consumed=0;
		next=0;
		nextReadID=0;
		buffer=null;
		tf.reset();
	}
	
	public static ArrayList<Read> toReadList(ByteFile tf, int maxReadsToReturn, long numericID){
		byte[] line=null;
		ArrayList<Read> list=new ArrayList<Read>(Data.min(8192, maxReadsToReturn));
		int added=0;
		
//		Read prev=null;
		
		for(line=tf.nextLine(); line!=null && added<maxReadsToReturn; line=tf.nextLine()){
			
			Read r=new Read(null, null, new String(line, StandardCharsets.US_ASCII), numericID);

//			if(interleaved){
//				if(prev==null){prev=r;}
//				else{
//					prev.mate=r;
//					r.mate=prev;
//					r.setPairnum(1);
//					list.add(prev);
//					added++;
//					numericID++;
//					prev=null;
//				}
//			}else
			{
				list.add(r);
				added++;
				numericID++;
			}

			if(added>=maxReadsToReturn){break;}
		}
		assert(list.size()<=maxReadsToReturn);
		return list;
	}
	
	/** Gets the filename of the input source.
	 * @return Filename from the underlying ByteFile */
	@Override
	public String fname(){return tf.name();}

	@Override
	public boolean paired() {return false;}
	
	/** Reports whether the stream has encountered an error.
	 * @return true if an error has been detected, false otherwise */
	@Override
	public boolean errorState(){return errorState;}

	private ArrayList<Read> buffer=null;
	private int next=0;
	
	private final ByteFile tf;
//	private final boolean interleaved;

	private final int BUF_LEN=Shared.bufferLen();;
	
	public long generated=0;
	public long consumed=0;
	private long nextReadID=0;
	
	public final boolean stdin;
	public static boolean verbose=false;

}
