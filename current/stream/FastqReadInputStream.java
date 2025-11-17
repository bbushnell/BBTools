package stream;

import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import hiseq.IlluminaHeaderParser2;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Reads FASTQ formatted sequence data and converts it to Read objects.
 * Provides buffered access to FASTQ files with support for interleaved paired-end reads.
 * Handles header shrinking for Illumina format optimization and custom parsing.
 * @author Brian Bushnell
 */
public class FastqReadInputStream extends ReadInputStream {
	
	/** Test method that reads and displays first read from a FASTQ file.
	 * @param args Command-line arguments; expects filename as first argument */
	public static void main(String[] args){
		
		FastqReadInputStream fris=new FastqReadInputStream(args[0], true);
		
		Read r=fris.nextList().get(0);
		System.out.println(r.toText(false));
		
	}
	
	/**
	 * Creates a FASTQ reader for the specified file.
	 * @param fname Input FASTQ filename
	 * @param allowSubprocess_ Whether to allow subprocess execution for compressed files
	 */
	public FastqReadInputStream(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.FASTQ, null, allowSubprocess_, false));
	}
	
	/**
	 * Creates a FASTQ reader from a FileFormat specification.
	 * Initializes buffering, interleaving detection, and header processing options.
	 * Handles custom parsing for synthetic data if FASTQ.PARSE_CUSTOM is enabled.
	 * @param ff FileFormat object specifying input source and options
	 */
	public FastqReadInputStream(FileFormat ff){
		if(verbose){System.err.println("FastqReadInputStream("+ff+")");}
		flag=(Shared.AMINO_IN ? Read.AAMASK : 0);
		stdin=ff.stdio();
		shrinkHeaders=FASTQ.SHRINK_HEADERS;
		if(!ff.fastq()){
			System.err.println("Warning: Did not find expected fastq file extension for filename "+ff.name());
		}
		
		if(FASTQ.PARSE_CUSTOM){
			try {
				String s[]=ff.name().split("_");
//				maxSnps=toNumber(s[3]);
//				maxInss=toNumber(s[4]);
//				maxDels=toNumber(s[5]);
//				maxSubs=toNumber(s[6]);
				
//				s=s[8].split("\\.");
//
//				s=s[0].split("-");
				
//				if(s.length!=8 && s.length!=9){
//					if(Shared.WINDOWS){System.err.println("Note: Filename indicates non-synthetic data, but FASTQ.PARSE_CUSTOM="+FASTQ.PARSE_CUSTOM);}
//				}
				
//				minChrom=Gene.toChromosome(s[0]);
//				maxChrom=Gene.toChromosome(s[1]);

			} catch (Exception e) {
				// TODO Auto-generated catch block
				//			e.printStackTrace();
//				if(Shared.WINDOWS){System.err.println("Note: Filename indicates non-synthetic data, but FASTQ.PARSE_CUSTOM="+FASTQ.PARSE_CUSTOM);}
			}
		}
//		interleaved=false;
		interleaved=ff.interleaved();//(ff.stdio()) ? FASTQ.FORCE_INTERLEAVED : FASTQ.isInterleaved(ff.name(), false);
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
	 * Refills the read buffer from the underlying FASTQ file.
	 * Reads up to BUF_LEN reads and applies header shrinking if enabled.
	 * Closes file when fewer reads than buffer size are returned.
	 */
	private synchronized void fillBuffer(){
		
		assert(buffer==null || next>=buffer.size());
		
		buffer=null;
		next=0;
		buffer=FASTQ.toReadList(tf, BUF_LEN, nextReadID, interleaved, flag);
		int bsize=(buffer==null ? 0 : buffer.size());
		nextReadID+=bsize;
		if(bsize<BUF_LEN){tf.close();}
		
		generated+=bsize;
		if(buffer==null){
			if(!errorState){
				errorState=true;
				System.err.println("Null buffer in FastqReadInputStream.");
			}
		}else if(shrinkHeaders) {//Kind of slow to regenerate strings but oh well
			for(Read r : buffer) {shrinkHeader(r, r.mate);}
		}
	}
	
	/**
	 * Shrinks Illumina FASTQ headers to compact coordinate format.
	 * Parses original header and replaces with shortened version if possible.
	 * Handles both single and paired reads with appropriate read numbers.
	 *
	 * @param r1 First read (or single read)
	 * @param r2 Second read in pair, or null for single reads
	 */
	private void shrinkHeader(Read r1, Read r2) {
		ihp.parse(r1.id);
		if(!ihp.canShrink()) {return;}
		bbh.clear().colon().colon().colon();
		ihp.appendCoordinates(bbh).space().append(1).colon();
		r1.id=bbh.toString();
		if(r2!=null) {
			bbh.set(bbh.length-2, (byte)'2');
			r2.id=bbh.toString();
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
	
	@Override
	public String fname(){return tf.name();}
	
	/** Return true if this stream has detected an error */
	@Override
	public boolean errorState(){return errorState || FASTQ.errorState();}

	/** Current buffer of reads loaded from file */
	private ArrayList<Read> buffer=null;
	/** Index of next read to return from buffer */
	private int next=0;
	
	/** Underlying file reader for FASTQ data */
	private final ByteFile tf;
	/** Whether input contains interleaved paired-end reads */
	private final boolean interleaved;
	/** Read flags for amino acid mode or other special processing */
	private final int flag;
	/** Whether to apply header shrinking to Illumina format reads */
	private boolean shrinkHeaders;
	
	/** ByteBuilder for constructing shrunken headers */
	private final ByteBuilder bbh=new ByteBuilder(128);
	/** Parser for Illumina FASTQ header format */
	private final IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();

	/** Buffer size in number of reads to load at once */
	private final int BUF_LEN=Shared.bufferLen();
	/** Maximum data size for buffer (currently unused for super-long reads) */
	private final long MAX_DATA=Shared.bufferData(); //TODO - lot of work for unlikely case of super-long fastq reads.  Must be disabled for paired-ends.

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
