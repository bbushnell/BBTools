package stream;

import java.util.ArrayList;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.LineParser1;
import shared.Shared;

/**
 * Reads GFA formatted sequence data and converts it to Read objects.
 * @author Brian Bushnell
 * @date November 21, 2025
 */
public class GfaReadInputStream extends ReadInputStream {
	
	/** Test method that reads and displays first read from a GFA file.
	 * @param args Command-line arguments; expects filename as first argument */
	public static void main(String[] args){
		
		GfaReadInputStream fris=new GfaReadInputStream(args[0], true);
		
		Read r=fris.nextList().get(0);
		System.out.println(r.toText(false));
		
	}
	
	/**
	 * Creates a GFA reader for the specified file.
	 * @param fname Input GFA filename
	 * @param allowSubprocess_ Whether to allow subprocess execution for compressed files
	 */
	public GfaReadInputStream(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.GFA, null, allowSubprocess_, false));
	}
	
	/**
	 * Creates a GFA reader from a FileFormat specification.
	 * Sets the amino-acid read flag (when Shared.AMINO_IN), opens the underlying ByteFile,
	 * and warns if the filename lacks a .gfa extension. GFA input is always single-ended.
	 * @param ff FileFormat object specifying input source and options
	 */
	public GfaReadInputStream(FileFormat ff){
		if(verbose){System.err.println("GfaReadInputStream("+ff+")");}
		flag=(Shared.AMINO_IN ? Read.AAMASK : 0);
		stdin=ff.stdio();
		if(!ff.gfa()){
			System.err.println("Warning: Did not find expected gfa file extension for filename "+ff.name());
		}
		tf=ByteFile.makeByteFile(ff);
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
		return list;
	}
	
	/**
	 * Refills the read buffer from the underlying GFA file.
	 * Reads up to BUF_LEN reads and applies header shrinking if enabled.
	 * Closes file when fewer reads than buffer size are returned.
	 */
	private synchronized void fillBuffer(){
		
		assert(buffer==null || next>=buffer.size());
		
		buffer=null;
		next=0;
		buffer=toReadList(tf, BUF_LEN, nextReadID, flag);
		int bsize=(buffer==null ? 0 : buffer.size());
		nextReadID+=bsize;
		if(bsize<BUF_LEN){tf.close();}
		
		generated+=bsize;
		//defensive/unreachable: toReadList always returns a (possibly empty) list, never null. An empty buffer -> nextList() maps empty->null as the EOF signal.
		if(buffer==null){
			if(!errorState){
				errorState=true;
				System.err.println("Null buffer in GfaReadInputStream.");
			}
		}
	}
	
	private ArrayList<Read> toReadList(final ByteFile bf, final int maxReadsToReturn,
			long numericID, final int flag){
		ArrayList<Read> list=new ArrayList<Read>(Data.min(400, maxReadsToReturn));
		//Only GFA 'S' (segment) lines become Reads; all other record types (H/L/P/C/W) are skipped. numericID advances only per added read.
		for(byte[] line=bf.nextLine(); line!=null && list.size()<maxReadsToReturn; line=bf.nextLine()) {
			if(line.length>0 && line[0]=='S') {
				lp.set(line);
				String id=lp.parseString(1);
				byte[] bases=lp.parseByteArray(2);
				//#002-refuted [stream/GfaReadInputStream#002]: a GFA '*' (no-sequence) segment does NOT yield a spurious 1bp read -> Read.validate (always -ea) rejects the junk base and crashes LOUD here (verified: reformat exit=1, no hang). This is the backstop that makes the unvalidated inline parse safe.
				Read r=new Read(bases, null, id, numericID++, flag);
				list.add(r);
			}
		}
		return list;
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
	public boolean paired() {return false;}
	
	@Override
	public String fname(){return tf.name();}
	
	/** Return true if this stream has detected an error */
	@Override
	//local-only is correct here: GFA has no delegated sub-parser carrying a separate error flag to OR-in (contrast Scarf's static FASTQ.errorState / Sam's streamer.errorState swallow #001). Inline parse errors surface as exceptions (e.g. Read.validate rejecting junk bases -> verified loud crash), not a swallowed flag.
	public boolean errorState(){return errorState;}

	/** Current buffer of reads loaded from file */
	private ArrayList<Read> buffer=null;
	/** Index of next read to return from buffer */
	private int next=0;//vestigial: never incremented (blockwise-only reader via nextList) -> always 0.
	
	/** Underlying file reader for GFA data */
	private final ByteFile tf;
	/** Read flags for amino acid mode or other special processing */
	private final int flag;

	/** Buffer size in number of reads to load at once */
	private final int BUF_LEN=Shared.bufferLen();
	/** Maximum data size for buffer (currently unused for super-long reads) */
	private final long MAX_DATA=Shared.bufferData(); //TODO - lot of work for unlikely case of super-long gfa reads.  Must be disabled for paired-ends.

	private final LineParser1 lp=new LineParser1('\t');
	
	/** Total number of reads loaded from file */
	public long generated=0;
	/** Total number of reads returned to caller */
	public long consumed=0;
	/** ID number to assign to next read loaded from file */
	private long nextReadID=0;
	
	/** Whether input is from standard input stream */
	public final boolean stdin;
	/** Whether to print verbose debugging information */
	public static boolean verbose=false;

}
