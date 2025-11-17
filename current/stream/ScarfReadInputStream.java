package stream;

import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Shared;

/**
 * Input stream for reading sequence data from SCARF format files.
 * SCARF is a tab-delimited sequence format containing read name, sequence,
 * quality scores, and metadata in columnar format.
 * Provides buffered reading with support for interleaved paired-end data.
 *
 * @author Brian Bushnell
 * @date July 2, 2025
 */
public class ScarfReadInputStream extends ReadInputStream {
	
	/**
	 * Test method for demonstrating ScarfReadInputStream functionality.
	 * Reads the first record from a SCARF file and prints it as text.
	 * @param args Command-line arguments where args[0] is the SCARF file path
	 */
	public static void main(String[] args){
		
		ScarfReadInputStream fris=new ScarfReadInputStream(args[0], true);
		
		Read r=fris.nextList().get(0);
		System.out.println(r.toText(false));
		
	}
	
	/**
	 * Creates a ScarfReadInputStream for the specified file.
	 * @param fname Path to the SCARF format input file
	 * @param allowSubprocess_ Whether to allow subprocess decompression for compressed files
	 */
	public ScarfReadInputStream(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.SCARF, null, allowSubprocess_, false));
	}
	
	/**
	 * Creates a ScarfReadInputStream from a FileFormat object.
	 * Validates the file format and initializes the underlying ByteFile reader.
	 * Sets interleaving mode based on FASTQ settings.
	 * @param ff FileFormat object containing file metadata and format information
	 */
	public ScarfReadInputStream(FileFormat ff){
		if(verbose){System.err.println("ScarfReadInputStream("+ff.name()+")");}
		
		stdin=ff.stdio();
		if(!ff.scarf()){
			System.err.println("Warning: Did not find expected scarf file extension for filename "+ff.name());
		}
		
		tf=ByteFile.makeByteFile(ff);
		
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
	
	/**
	 * Fills the internal read buffer by parsing SCARF data from the file.
	 * Reads up to BUF_LEN reads at a time and converts them using FASTQ parser.
	 * Closes the file automatically when fewer reads than buffer size are returned.
	 * Updates read ID counter and sets error state if buffer creation fails.
	 */
	private synchronized void fillBuffer(){
		
		assert(buffer==null || next>=buffer.size());
		
		buffer=null;
		next=0;
		
		buffer=FASTQ.toScarfReadList(tf, BUF_LEN, nextReadID, interleaved);
		int bsize=(buffer==null ? 0 : buffer.size());
		nextReadID+=bsize;
		if(bsize<BUF_LEN){tf.close();}
		
		generated+=bsize;
		if(buffer==null){
			if(!errorState){
				errorState=true;
				System.err.println("Null buffer in ScarfReadInputStream.");
			}
		}
	}
	
	@Override
	public boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		errorState|=tf.close();
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		return errorState;
	}

	@Override
	public synchronized void restart() {
		generated=0;
		consumed=0;
		next=0;
		nextReadID=0;
		buffer=null;
		tf.reset();
	}

	@Override
	public boolean paired() {return interleaved;}
	
	/** Return true if this stream has detected an error */
	@Override
	public boolean errorState(){return errorState || FASTQ.errorState();}
	
	@Override
	public String fname(){return tf.name();}

	/** Internal buffer holding the current batch of reads */
	private ArrayList<Read> buffer=null;
	/** Index of the next read to return from the current buffer */
	private int next=0;
	
	/** Underlying file reader for accessing the SCARF file data */
	private final ByteFile tf;
	/** Whether the input contains interleaved paired-end reads */
	private final boolean interleaved;

	/** Maximum number of reads to buffer at once */
	private final int BUF_LEN=Shared.bufferLen();;
	/** Maximum data size for buffering (currently unused for SCARF format) */
	private final long MAX_DATA=Shared.bufferData(); //TODO - lot of work for unlikely case of super-long scarf reads.  Must be disabled for paired-ends.

	/** Total number of reads generated from the input file */
	public long generated=0;
	/** Total number of reads consumed by the client */
	public long consumed=0;
	/** Numeric ID to assign to the next read parsed from the file */
	private long nextReadID=0;
	
	/** Whether the input is being read from standard input */
	public final boolean stdin;
	/** Controls verbose output for debugging stream operations */
	public static boolean verbose=false;

}
