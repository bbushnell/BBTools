package stream;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.FileFormat;
import shared.Shared;

/**
 * @author Brian Bushnell
 * @date June 1, 2016
 *
 */
public class HeaderInputStream extends ReadInputStream {
	
	/**
	 * Program entry point for testing HeaderInputStream functionality.
	 * Reads first header from specified file and prints it.
	 * @param args Command-line arguments where args[0] is the input filename
	 */
	public static void main(String[] args){
		
		HeaderInputStream his=new HeaderInputStream(args[0], true);
		
		Read r=his.nextList().get(0);
		System.out.println(r.toText(false));
		his.close();
		
	}
	
	/**
	 * Creates HeaderInputStream from filename with subprocess option.
	 * @param fname Input filename to read headers from
	 * @param allowSubprocess_ Whether to allow subprocess for compressed files
	 */
	public HeaderInputStream(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.FASTQ, null, allowSubprocess_, false));
	}

	
	/**
	 * Creates HeaderInputStream from FileFormat specification.
	 * Initializes ByteFile reader and sets stdin flag based on format.
	 * @param ff FileFormat containing input source and format information
	 */
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
	 * Fills the internal buffer with reads from the input file.
	 * Reads up to BUF_LEN header lines and converts them to Read objects.
	 * Closes file when fewer than BUF_LEN reads are obtained.
	 */
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
	
	/**
	 * Converts lines from ByteFile into a list of Read objects.
	 * Each line becomes a Read with only the name field populated from the line content.
	 * Reads up to maxReadsToReturn lines or until end of file.
	 *
	 * @param tf ByteFile to read lines from
	 * @param maxReadsToReturn Maximum number of reads to create
	 * @param numericID Starting numeric ID for the first read
	 * @return ArrayList of Read objects created from input lines
	 */
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
	
	@Override
	public String fname(){return tf.name();}

	@Override
	public boolean paired() {return false;}
	
	/** Return true if this stream has detected an error */
	@Override
	public boolean errorState(){return errorState;}

	/** Buffer holding current batch of Read objects */
	private ArrayList<Read> buffer=null;
	/** Index of next Read to return from buffer */
	private int next=0;
	
	/** ByteFile reader for the input source */
	private final ByteFile tf;
//	private final boolean interleaved;

	/** Buffer size for batch reading from shared configuration */
	private final int BUF_LEN=Shared.bufferLen();;
	
	/** Total number of reads generated from input */
	public long generated=0;
	/** Total number of reads consumed by caller */
	public long consumed=0;
	/** Numeric ID to assign to the next read created */
	private long nextReadID=0;
	
	/** Whether input is coming from standard input */
	public final boolean stdin;
	/** Controls verbose output during stream operations */
	public static boolean verbose=false;

}
