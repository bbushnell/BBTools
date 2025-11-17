package stream;

import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import dna.Data;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Shared;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Abstract base class for writing read data to various file formats in a separate thread.
 * Supports output to FASTQ, FASTA, SAM/BAM, and other formats with buffered writing.
 * Uses a producer-consumer pattern with a blocking queue for thread-safe operation.
 * @author Brian Bushnell
 */
public abstract class ReadStreamWriter extends Thread {
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	protected ReadStreamWriter(FileFormat ff, String qfname_, boolean read1_, int bufferSize, 
			CharSequence header, boolean buffered, boolean useSharedHeader){
//		assert(false) : useSharedHeader+", "+header;
		assert(ff!=null);
		assert(ff.write()) : "FileFormat is not in write mode for "+ff.name();
		
		assert(!ff.text() && !ff.unknownFormat()) : "Unknown format for "+ff;
		OUTPUT_FASTQ=ff.fastq();
		OUTPUT_FASTA=ff.fasta();
		OUTPUT_FASTR=ff.fastr();
//		boolean bread=(ext==TestFormat.txt);
		OUTPUT_SAM=ff.samOrBam();
		OUTPUT_BAM=ff.bam();
		OUTPUT_ATTACHMENT=ff.attachment();
		OUTPUT_HEADER=ff.header();
		OUTPUT_ONELINE=ff.oneline();
		SITES_ONLY=ff.sites();
		OUTPUT_STANDARD_OUT=ff.stdio();
		FASTA_WRAP=Shared.FASTA_WRAP;
		assert(((OUTPUT_SAM ? 1 : 0)+(OUTPUT_FASTQ ? 1 : 0)+(OUTPUT_FASTA ? 1 : 0)+(OUTPUT_ATTACHMENT ? 1 : 0)+
				(OUTPUT_HEADER ? 1 : 0)+(OUTPUT_ONELINE ? 1 : 0)+(SITES_ONLY ? 1 : 0))<=1) :
			OUTPUT_SAM+", "+SITES_ONLY+", "+OUTPUT_FASTQ+", "+OUTPUT_FASTA+", "+OUTPUT_ATTACHMENT+", "+OUTPUT_HEADER+", "+OUTPUT_ONELINE;
		
		fname=ff.name();
		qfname=qfname_;
		read1=read1_;
		allowSubprocess=ff.allowSubprocess();
		boolean append=ff.append();
//		assert(fname==null || (fname.contains(".sam") || fname.contains(".bam"))==OUTPUT_SAM) : "Outfile name and sam output mode flag disagree: "+fname;
		assert(read1 || !OUTPUT_SAM) : "Attempting to output paired reads to different sam files.";
		
		if(qfname==null){
			myQOutstream=null;
		}else{
			myQOutstream=ReadWrite.getOutputStream(qfname, (ff==null ? false : ff.append()), buffered, allowSubprocess);
		}

		final boolean supressHeader=(NO_HEADER || (ff.append() && ff.exists()));
		final boolean supressHeaderSequences=(NO_HEADER_SEQUENCES || supressHeader);
		final boolean RSSamWriter=ff.samOrBam() && ReadWrite.USE_READ_STREAM_SAM_WRITER;
		
		if(fname==null && !OUTPUT_STANDARD_OUT){
			myOutstream=null;
		}else if(RSSamWriter) {
			myOutstream=null;
		}else{
			if(OUTPUT_STANDARD_OUT){myOutstream=System.out;}
			
			else if(!ff.bam() || !Data.BAM_SUPPORT_OUT()){
				myOutstream=ReadWrite.getOutputStream(fname, append, buffered, allowSubprocess);
			}else{
				myOutstream=ReadWrite.getBamOutputStream(fname, append);
			}
			
			if(header!=null && !supressHeader){
					byte[] temp=new byte[header.length()];
					for(int i=0; i<temp.length; i++){temp[i]=(byte)header.charAt(i);}
					try {
						myOutstream.write(temp);
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
			}else if(OUTPUT_SAM && !supressHeader){
				if(useSharedHeader){
//					assert(false);
					ArrayList<byte[]> list=SamReadInputStream.getSharedHeader(true);
					if(list==null){
						System.err.println("Header was null.");
					}else{
						try {
							if(supressHeaderSequences){
								for(byte[] line : list){
									boolean sq=(line!=null && line.length>2 && line[0]=='@' && line[1]=='S' && line[2]=='Q' && line[3]=='\t');
									if(!sq){
										myOutstream.write(line);
										myOutstream.write('\n');
									}
								}
							}else{
								for(byte[] line : list){
									myOutstream.write(line);
									myOutstream.write('\n');
								}
							}
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}else{
						ByteBuilder bb=new ByteBuilder(4096);
						SamHeader.header0B(bb);
						bb.nl();
						int a=(MINCHROM==-1 ? 1 : MINCHROM);
						int b=(MAXCHROM==-1 ? Data.numChroms : MAXCHROM);
						if(!supressHeaderSequences){
							for(int chrom=a; chrom<=b; chrom++){
								SamHeader.printHeader1B(chrom, chrom, bb, myOutstream);
							}
						}
						SamHeader.header2B(bb);
						bb.nl();


						try {
							if(bb.length>0){myOutstream.write(bb.array, 0, bb.length);}
						} catch (IOException e) {
							KillSwitch.exceptionKill(e);
						}
				}
			}else if(ff.bread() && !supressHeader){
					try {
						myOutstream.write(("#"+Read.header()).getBytes());
					} catch (IOException e) {
						KillSwitch.exceptionKill(e);
					}
			}
		}
		
		assert(bufferSize>=1);
		queue=new ArrayBlockingQueue<Job>(bufferSize);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public abstract void run();

	/** Sends a poison pill to shutdown the writer thread gracefully.
	 * Creates a job that signals the thread to terminate and close streams. */
	public final synchronized void poison(){
		addJob(new Job(null, false, true, nextID++));
	}

	public final synchronized void addList(ListNum<Read> ln){
		assert(ln.id==nextID) : ln.id+", "+nextID;
		Job j=new Job(ln.list, ln.last(), ln.poison(), ln.id);
		nextID=ln.id+1;
		addJob(j);
	}

	/** Adds a list of reads to the write queue using default output streams.
	 * @param list List of reads to write */
	public final synchronized void addList(ArrayList<Read> list){
		Job j=new Job(list, false, false, nextID++);
		addJob(j);
	}
	
	/**
	 * Adds a job to the write queue, blocking until space is available.
	 * Handles InterruptedException by retrying the operation.
	 * @param j Job containing read list and output configuration
	 */
	private final synchronized void addJob(Job j){
//		System.err.println("Got job "+(j.list==null ? "null" : j.list.size()));
		boolean success=false;
		while(!success){
			try {
				queue.put(j);
				success=true;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				assert(!queue.contains(j)); //Hopefully it was not added.
			}
		}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	protected static final ByteBuilder toQualityB(final byte[] quals, final int len, 
			final int wrap, final ByteBuilder bb){
		if(quals==null){return fakeQualityB(30, len, wrap, bb);}
		assert(quals.length==len);
		bb.ensureExtra(NUMERIC_QUAL ? len*3+1 : len+1);
		if(NUMERIC_QUAL){
			if(len>0){bb.append((int)quals[0]);}
			for(int i=1, w=1; i<len; i++, w++){
				if(w>=wrap){
					bb.nl();
					w=0;
				}else{
					bb.append(' ');
				}
				bb.append((int)quals[i]);
			}
		}else{
			final byte b=FASTQ.ASCII_OFFSET_OUT;
			for(int i=0; i<len; i++){
				bb.append(b+quals[i]);
			}
		}
		return bb;
	}
	
	protected static final ByteBuilder fakeQualityB(final int q, final int len, 
			final int wrap, final ByteBuilder bb){
		bb.ensureExtra(NUMERIC_QUAL ? len*3+1 : len+1);
		if(NUMERIC_QUAL){
			if(len>0){bb.append(q);}
			for(int i=1, w=1; i<len; i++, w++){
				if(w>=wrap){
					bb.nl();
					w=0;
				}else{
					bb.append(' ');
				}
				bb.append(q);
			}
		}else{
			byte c=(byte)(q+FASTQ.ASCII_OFFSET_OUT);
			for(int i=0; i<len; i++){bb.append(c);}
		}
		return bb;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/** Returns the output file name */
	public String fname(){return fname;}
	/** Returns the number of reads written to output */
	public long readsWritten(){return readsWritten;}
	/** Returns the number of bases written to output */
	public long basesWritten(){return basesWritten;}

	/** Return true if this stream has detected an error */
	public final boolean errorState(){return errorState;}
	/** Return true if this stream has finished */
	public final boolean finishedSuccessfully(){return finishedSuccessfully;}
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** TODO */
	protected boolean errorState=false;
	/** Flag indicating if writing finished without errors */
	protected boolean finishedSuccessfully=false;
	
	/** True if output format is SAM */
	public final boolean OUTPUT_SAM;
	/** True if output format is BAM */
	public final boolean OUTPUT_BAM;
	/** True if output format is FASTQ */
	public final boolean OUTPUT_FASTQ;
	/** True if output format is FASTA */
	public final boolean OUTPUT_FASTA;
	/** True if output format is FASTR */
	public final boolean OUTPUT_FASTR;
	/** True if output format includes header information */
	public final boolean OUTPUT_HEADER;
	/** True if output format supports attachments */
	public final boolean OUTPUT_ATTACHMENT;
	/** True if output format uses single-line format */
	public final boolean OUTPUT_ONELINE;
	/** True if writing to standard output */
	public final boolean OUTPUT_STANDARD_OUT;
	/** True if outputting sites information only */
	public final boolean SITES_ONLY;
	/** True if output format is interleaved */
	public boolean OUTPUT_INTERLEAVED=false;
	
	/** Line wrap length for FASTA format output */
	protected final int FASTA_WRAP;
	
	/** Whether to allow spawning external processes like samtools */
	protected final boolean allowSubprocess;
	
	/** True if processing first reads in pairs */
	protected final boolean read1;
	/** Output file name */
	protected final String fname;
	/** Quality output file name */
	protected final String qfname;
	/** Primary output stream */
	protected final OutputStream myOutstream;
	/** Quality output stream */
	protected final OutputStream myQOutstream;
	/** Thread-safe queue for write jobs */
	protected final ArrayBlockingQueue<Job> queue;
	
	/** Count of reads written to output */
	protected long readsWritten=0;
	/** Count of bases written to output */
	protected long basesWritten=0;
	protected long nextID=0;
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Minimum chromosome number for SAM header generation */
	public static int MINCHROM=-1; //For generating sam header
	/** Maximum chromosome number for SAM header generation */
	public static int MAXCHROM=-1; //For generating sam header
	/** True to output quality scores in numeric format instead of ASCII */
	public static boolean NUMERIC_QUAL=true;
	/** True to output secondary alignments in SAM format */
	public static boolean OUTPUT_SAM_SECONDARY_ALIGNMENTS=false;
	
	/** True to ignore assertions about read pairing */
	public static boolean ignorePairAssertions=false;
	/** True to enable CIGAR string validation assertions */
	public static boolean ASSERT_CIGAR=false;
	/** True to suppress header output */
	public static boolean NO_HEADER=false;
	/** True to suppress sequence information in headers */
	public static boolean NO_HEADER_SEQUENCES=false;
	/** True to use attached SAM line information */
	public static boolean USE_ATTACHED_SAMLINE=false;
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	//TODO: Should be replaced with ListNum
	protected static class Job implements HasID{
		
		public Job(ArrayList<Read> list_, boolean closeWhenDone_,
				boolean poisonThread_, long id_){
			list=list_;
			close=closeWhenDone_;
			poison=poisonThread_;
			id=id_;
		}
		
		/*--------------------------------------------------------------*/
		
		/** Checks if this job has no reads to process.
		 * @return true if read list is null or empty */
		public boolean isEmpty(){return list==null || list.isEmpty();}
		/** List of reads to write */
		public final ArrayList<Read> list;
		/** Whether to close streams after processing this job */
		public final boolean close;
		/** Whether this job should shutdown the writer thread */
		public final boolean poison;
		public final long id;
		@Override
		public long id(){return id;}
		@Override
		public boolean poison(){return poison;}
		@Override
		public boolean last(){return close;}
		@Override
		public Job makePoison(long id_){return new Job(null, false, true, id_);}
		@Override
		public Job makeLast(long id_){return new Job(null, true, false, id_);}
		
	}
	
}
