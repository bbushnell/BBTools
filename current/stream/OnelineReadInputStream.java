package stream;

import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;

/**
 * Read input stream for tab-delimited oneline format files.
 * Each line contains a sequence identifier, tab character, then sequence bases.
 * Supports both single and paired-end (interleaved) reads.
 * @author Brian Bushnell
 */
public class OnelineReadInputStream extends ReadInputStream {
	
	/** Program entry point for testing oneline format reading.
	 * @param args Command-line arguments containing filename */
	public static void main(String[] args){
		
		OnelineReadInputStream fris=new OnelineReadInputStream(args[0], true);
		
		Read r=fris.nextList().get(0);
		System.out.println(r.toText(false));
		
	}
	
	/**
	 * Constructs a OnelineReadInputStream from a filename.
	 * @param fname Input filename
	 * @param allowSubprocess_ Whether to allow subprocess execution for compressed files
	 */
	public OnelineReadInputStream(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.ONELINE, null, allowSubprocess_, false));
	}

	
	/**
	 * Constructs a OnelineReadInputStream from a FileFormat object.
	 * Validates oneline format extension and initializes the underlying ByteFile.
	 * @param ff FileFormat specifying the input file characteristics
	 */
	public OnelineReadInputStream(FileFormat ff){
		if(verbose){System.err.println("FastqReadInputStream("+ff+")");}
		
		stdin=ff.stdio();
		if(!ff.oneline()){
			System.err.println("Warning: Did not find expected oneline file extension for filename "+ff.name());
		}
		
//		interleaved=false;
//		assert(false) : "TODO: Detect interleaved.";
		interleaved=FASTQ.FORCE_INTERLEAVED; //(ff.stdio()) ? FASTQ_X.FORCE_INTERLEAVED : FASTQ_X.isInterleaved(ff.name(), false);
		
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
	
	/** Fills the internal buffer with reads by parsing oneline format data.
	 * Closes the file when fewer reads than buffer size are returned. */
	private synchronized void fillBuffer(){
		
		assert(buffer==null || next>=buffer.size());
		
		buffer=null;
		next=0;
		
		buffer=toReadList();
		int bsize=(buffer==null ? 0 : buffer.size());
//		nextReadID+=bsize;
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
	 * Parses oneline format data into Read objects.
	 * Each line is split at the last tab character to separate ID from sequence.
	 * Handles both single-end and paired-end (interleaved) reads.
	 * @return ArrayList of parsed Read objects
	 */
	private ArrayList<Read> toReadList(){
		ArrayList<Read> list=new ArrayList<Read>(BUF_LEN);
		Read r1=null, r2=null;
		long sum=0;
		for(byte[] line=tf.nextLine(); line!=null; line=tf.nextLine()){
			int index=Tools.lastIndexOf(line, (byte)'\t');
			String id=new String(line, 0, index);
			byte[] bases=KillSwitch.copyOfRange(line, index+1, line.length);
			sum+=bases.length;
			Read r=new Read(bases, null, id, nextReadID);
			if(r1==null){
				r1=r;
			}else{
				r2=r;
				r1.mate=r2;
				r2.mate=r1;
			}
			if(interleaved==(r2!=null)){
				list.add(r1);
				r1=r2=null;
				nextReadID++;
				if(list.size()>=BUF_LEN || sum>=MAX_DATA){break;}
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
	public boolean paired() {return interleaved;}
	
	/** Return true if this stream has detected an error */
	@Override
	public boolean errorState(){return errorState;}
	
	@Override
	public String fname(){return tf.name();}

	/** Buffer for storing parsed reads */
	private ArrayList<Read> buffer=null;
	/** Index of next read to return from buffer */
	private int next=0;
	
	/** Underlying file reader for oneline format data */
	private final ByteFile tf;
	/** Whether reads are paired and interleaved */
	private final boolean interleaved;

	/** Maximum number of reads to buffer at once */
	private final int BUF_LEN=Shared.bufferLen();;
	/** Maximum amount of sequence data to buffer in bytes */
	private final long MAX_DATA=Shared.bufferData(); //TODO - lot of work for unlikely case of super-long fastq reads.  Must be disabled for paired-ends.

	/** Total number of reads generated from this stream */
	public long generated=0;
	/** Total number of reads consumed by caller */
	public long consumed=0;
	/** ID number for the next read to be created */
	private long nextReadID=0;
	
	/** Whether input is from standard input */
	public final boolean stdin;
	/** Whether to print verbose debugging information */
	public static boolean verbose=false;

}
