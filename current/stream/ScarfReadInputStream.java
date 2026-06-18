package stream;

import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;

public class ScarfReadInputStream extends ReadInputStream {
	
	public static void main(String[] args){
		
		ScarfReadInputStream fris=new ScarfReadInputStream(args[0], true);
		
		Read r=fris.nextList().get(0);
		System.out.println(r.toText(false));
		
	}
	
	public ScarfReadInputStream(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.SCARF, null, allowSubprocess_, false));
	}
	
	public ScarfReadInputStream(FileFormat ff){
		if(verbose){System.err.println("ScarfReadInputStream("+ff.name()+")");}
		
		stdin=ff.stdio();
		if(!ff.scarf()){
			System.err.println("Warning: Did not find expected scarf file extension for filename "+ff.name());
		}
		
		tf=ByteFile.makeByteFile(ff);
		
		//SCARF pairing is driven solely by the global FORCE_INTERLEAVED flag; per-file auto-detect (the commented alternative) is disabled.
		interleaved=FASTQ.FORCE_INTERLEAVED;//((tf.is()==System.in || stdin) ? FASTQ.FORCE_INTERLEAVED : FASTQ.isInterleaved(tf.name));
//		assert(false) : interleaved;
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
	
	@Override
	public synchronized ArrayList<Read> nextList() {
		if(next!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");}
		if(buffer==null || next>=buffer.size()){fillBuffer();}
		ArrayList<Read> list=buffer;
		buffer=null;
		if(list!=null && list.size()==0){list=null;}
		consumed+=(list==null ? 0 : list.size());
//		System.err.println(hashCode()+" produced "+r[0].numericID);
		return list;
	}
	
	private synchronized void fillBuffer(){
		
		assert(buffer==null || next>=buffer.size());
		
		buffer=null;
		next=0;
		
		buffer=FASTQ.toScarfReadList(tf, BUF_LEN, nextReadID, interleaved);
		int bsize=(buffer==null ? 0 : buffer.size());
		nextReadID+=bsize;
		if(bsize<BUF_LEN){tf.close();}
		
		generated+=bsize;
		//defensive/unreachable: FASTQ.toScarfReadList always returns a (possibly empty) list, never null. Clean EOF -> empty buffer (bsize=0<BUF_LEN closes tf), and nextList() maps empty->null as the EOF signal.
		if(buffer==null){
			if(!errorState){
				errorState=true;
				System.err.println("Null buffer in ScarfReadInputStream.");
			}
		}
	}
	
	/**
	 * Closes the input stream and releases associated resources.
	 * Updates error state based on the success of the close operation.
	 * @return true if an error occurred during closing, false otherwise
	 */
	@Override
	public boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		errorState|=tf.close();
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		return errorState;
	}

	/** Resets the stream to the beginning for re-reading the file.
	 * Clears all counters, buffers, and resets the underlying file reader. */
	@Override
	public synchronized void restart() {
		generated=0;
		consumed=0;
		next=0;
		nextReadID=0;
		buffer=null;
		tf.reset();
	}

	/** Indicates whether this stream contains interleaved paired-end reads.
	 * @return true if reads are interleaved pairs, false for single-end reads */
	@Override
	public boolean paired() {return interleaved;}
	
	/**
	 * Checks if the stream is in an error state.
	 * Combines local error state with FASTQ parser error state.
	 * @return true if any errors have been detected, false otherwise
	 */
	@Override
	public boolean errorState(){return errorState || FASTQ.errorState();}//anti-swallow: ORs the static (process-global) FASTQ parser error state so a parse failure isn't dropped -> contrast SamReadInputStream#001. Over-reports rather than under (safe direction).
	
	@Override
	public String fname(){return tf.name();}

	private ArrayList<Read> buffer=null;
	private int next=0;//vestigial: never incremented (blockwise-only reader via nextList) -> always 0, so the next-based single-read guards above are constant.
	
	private final ByteFile tf;
	private final boolean interleaved;

	private final int BUF_LEN=Shared.bufferLen();;
	private final long MAX_DATA=Shared.bufferData(); //TODO - lot of work for unlikely case of super-long scarf reads.  Must be disabled for paired-ends.

	public long generated=0;
	public long consumed=0;
	private long nextReadID=0;
	
	public final boolean stdin;
	public static boolean verbose=false;

}
