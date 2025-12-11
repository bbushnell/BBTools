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
 * Base class for threaded read writers.
 * Wraps a blocking queue of write jobs and supports multiple output formats
 * (FASTQ/FASTA/SAM/BAM) with optional header emission and subsampling hooks.
 */
public abstract class ReadStreamWriter extends Thread {
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Builds a writer for the requested format and initializes output streams, headers,
	 * and the job queue.
	 *
	 * @param ff Output file format
	 * @param qfname_ Optional quality filename
	 * @param read1_ True if this writer handles read1
	 * @param bufferSize Queue capacity for pending write jobs
	 * @param header Optional header text to emit
	 * @param buffered True to buffer output streams
	 * @param useSharedHeader True to reuse a shared SAM header when present
	 */
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
	
	/** Entry point for the writer thread; implemented by concrete subclasses. */
	@Override
	public abstract void run();

	/** Enqueues a poison job to signal shutdown and close streams. */
	public final synchronized void poison(){
		addJob(new Job(null, false, true, nextID++));
	}

	/**
	 * Adds a list of reads from a ListNum wrapper, preserving ordering metadata.
	 */
	public final synchronized void addList(ListNum<Read> ln){
		assert(ln.id==nextID) : ln.id+", "+nextID;
		Job j=new Job(ln.list, ln.last(), ln.poison(), ln.id);
		nextID=ln.id+1;
		addJob(j);
	}

	/** Adds a list of reads to the queue using default output streams. */
	public final synchronized void addList(ArrayList<Read> list){
		Job j=new Job(list, false, false, nextID++);
		addJob(j);
	}
	
	/** Adds a job to the blocking queue, retrying on interruption. */
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
	
	/** Renders quality scores into a ByteBuilder in numeric or ASCII form. */
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
	
	
	/** @return Output filename (may be null for stdout). */
	public String fname(){return fname;}
	/** @return Number of reads written so far. */
	public long readsWritten(){return readsWritten;}
	/** @return Number of bases written so far. */
	public long basesWritten(){return basesWritten;}

	/** @return True if the writer encountered an error. */
	public final boolean errorState(){return errorState;}
	/** @return True when the writer finished without errors. */
	public final boolean finishedSuccessfully(){return finishedSuccessfully;}
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** True if an error was encountered while writing. */
	protected boolean errorState=false;
	/** True if writing completed without errors. */
	protected boolean finishedSuccessfully=false;
	
	/** True if the writer is emitting SAM format. */
	public final boolean OUTPUT_SAM;
	/** True if the writer is emitting BAM format. */
	public final boolean OUTPUT_BAM;
	/** True if the writer is emitting FASTQ format. */
	public final boolean OUTPUT_FASTQ;
	/** True if the writer is emitting FASTA format. */
	public final boolean OUTPUT_FASTA;
	/** True if the writer is emitting FASTR format. */
	public final boolean OUTPUT_FASTR;
	/** True if the output includes header records. */
	public final boolean OUTPUT_HEADER;
	/** True if the output includes attachment data. */
	public final boolean OUTPUT_ATTACHMENT;
	/** True if the output uses single-line formatting. */
	public final boolean OUTPUT_ONELINE;
	/** True if writing to standard output instead of a file. */
	public final boolean OUTPUT_STANDARD_OUT;
	/** True if only site information is being written. */
	public final boolean SITES_ONLY;
	/** True if interleaving read pairs in a single output stream. */
	public boolean OUTPUT_INTERLEAVED=false;
	
	/** Line wrap length for FASTA output. */
	protected final int FASTA_WRAP;
	
	/** True if subprocess-based compression is permitted. */
	protected final boolean allowSubprocess;
	
	/** True if this writer handles first-of-pair reads. */
	protected final boolean read1;
	/** Output filename for the primary stream. */
	protected final String fname;
	/** Output filename for qualities (when applicable). */
	protected final String qfname;
	/** Primary output stream for reads. */
	protected final OutputStream myOutstream;
	/** Output stream for quality data when separate files are used. */
	protected final OutputStream myQOutstream;
	/** Thread-safe queue of write jobs. */
	protected final ArrayBlockingQueue<Job> queue;
	
	/** Count of reads written to output so far. */
	protected long readsWritten=0;
	/** Count of bases written to output so far. */
	protected long basesWritten=0;
	/** Identifier for the next queued job to preserve ordering. */
	protected long nextID=0;
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Minimum chromosome number used for SAM header generation. */
	public static int MINCHROM=-1; //For generating sam header
	/** Maximum chromosome number used for SAM header generation. */
	public static int MAXCHROM=-1; //For generating sam header
	/** Output qualities in numeric form when true; otherwise ASCII. */
	public static boolean NUMERIC_QUAL=true;
	/** Emit secondary alignments in SAM output when true. */
	public static boolean OUTPUT_SAM_SECONDARY_ALIGNMENTS=false;
	
	/** Relax pairing assertions when writing reads. */
	public static boolean ignorePairAssertions=false;
	/** Enable assertions that validate CIGAR strings. */
	public static boolean ASSERT_CIGAR=false;
	/** Suppress header output when true. */
	public static boolean NO_HEADER=false;
	/** Suppress sequence lines in headers when true. */
	public static boolean NO_HEADER_SEQUENCES=false;
	/** Use attached SamLine data when emitting SAM records. */
	public static boolean USE_ATTACHED_SAMLINE=false;
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	//TODO: Should be replaced with ListNum
	/**
	 * Write job holding a read list plus control flags for shutdown/order preservation.
	 */
	protected static class Job implements HasID{
		
		public Job(ArrayList<Read> list_, boolean closeWhenDone_,
				boolean poisonThread_, long id_){
			list=list_;
			close=closeWhenDone_;
			poison=poisonThread_;
			id=id_;
		}
		
		/*--------------------------------------------------------------*/
		
		public boolean isEmpty(){return list==null || list.isEmpty();}
		public final ArrayList<Read> list;
		public final boolean close;
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
