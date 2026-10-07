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
 * Base class for threaded read writers with a bounded FIFO job queue.
 * Captures output-format settings and initializes legacy output streams and headers.
 * Concrete subclasses consume jobs and own record rendering, counters and completion.
 * When the SAM/BAM adapter flag is enabled, primary output initialization is delegated.
 */
public abstract class ReadStreamWriter extends Thread{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Captures format options, initializes outputs and creates the job queue; does not start this thread.
	 * A separate quality output is opened when requested. Primary output and header setup
	 * are skipped for the SAM/BAM adapter route selected by ReadWrite.USE_READ_STREAM_SAM_WRITER.
	 * Legacy headers are suppressed by NO_HEADER or append mode on an existing file.
	 * Otherwise explicit header text takes precedence, with each character narrowed to one byte;
	 * SAM/BAM output without explicit text uses shared or generated header records.
	 * NO_HEADER_SEQUENCES suppresses sequence records in those shared/generated headers.
	 * Default native-text preambles belong to ReadStreamByteWriter, after interleaving
	 * is selected, so its pairing marker can be the first line.
	 * @param ff Nonnull output descriptor in write mode
	 * @param qfname_ Optional separate quality filename, or null for no quality stream
	 * @param read1_ True for first-of-pair output; SAM/BAM requires true
	 * @param bufferSize Positive capacity for pending write jobs
	 * @param header Optional explicit legacy header text; no newline is added
	 * @param buffered Whether to buffer streams opened through ReadWrite.getOutputStream
	 * @param useSharedHeader Whether legacy SAM/BAM header setup requests the shared header
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
		}else if(RSSamWriter){
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
				try{
					myOutstream.write(temp);
				}catch(IOException e){
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
						try{
							if(supressHeaderSequences){
								for(byte[] line : list){
									boolean sq=(line!=null && line.length>3 && line[0]=='@' && line[1]=='S' && line[2]=='Q' && line[3]=='\t');
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
						}catch(IOException e){
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

					try{
						if(bb.length>0){myOutstream.write(bb.array, 0, bb.length);}
					}catch(IOException e){
						KillSwitch.exceptionKill(e);
					}
				}
			}
			//Native preambles are written by ReadStreamByteWriter after OUTPUT_INTERLEAVED is set.
			//A constructor prefix would hide the first-line marker required by RTextInputStream.isInterleaved.
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

	/** Enqueues a poison marker with the next ID; shutdown and stream cleanup belong to the consumer. */
	public final synchronized void poison(){
		addJob(new Job(null, false, true, nextID++));
	}

	/**
	 * Marks failure, discards queued jobs and offers a poison marker after acquiring this monitor.
	 * The queue offer does not wait for capacity; this method does not join the consumer.
	 * Unlike poison(), it bypasses the retrying addJob path. The consumer queue is FIFO,
	 * so the marker need not retain the ID of a discarded normal job.
	 */
	public final synchronized void abortNow(){
		errorState=true;
		finishedSuccessfully=false;
		queue.clear();
		queue.offer(new Job(null, false, true, nextID++));
	}

	/**
	 * Enqueues the wrapper's list, ID and terminal flags without copying the list.
	 * Requires the wrapper ID to equal the next expected ID, then advances that ID.
	 * @param ln Nonnull wrapper whose ID follows the previously submitted job
	 */
	public final synchronized void addList(ListNum<Read> ln){
		assert(ln.id==nextID) : ln.id+", "+nextID;
		Job j=new Job(ln.list, ln.last(), ln.poison(), ln.id);
		nextID=ln.id+1;
		addJob(j);
	}

	/** Enqueues the supplied list without copying it, assigning the next ID and no terminal flags. */
	public final synchronized void addList(ArrayList<Read> list){
		Job j=new Job(list, false, false, nextID++);
		addJob(j);
	}

	/** Adds a job to the blocking queue, retrying on interruption. */
	private final synchronized void addJob(Job j){
//		System.err.println("Got job "+(j.list==null ? "null" : j.list.size()));
		boolean success=false;
		while(!success){
			try{
				queue.put(j);
				success=true;
			}catch(InterruptedException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
				assert(!queue.contains(j)); //Hopefully it was not added.
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Appends quality scores to the supplied builder without changing the input array.
	 * Null qualities use fake score 30 for each requested position. Numeric mode separates
	 * scores with spaces and replaces separators with newlines at wrap boundaries;
	 * ASCII mode appends one offset-encoded byte per score and ignores wrap.
	 * Neither mode appends a trailing newline.
	 * @param quals Quality scores, or null to synthesize them
	 * @param len Number of scores; must equal quals.length when quals is nonnull
	 * @param wrap Numeric scores per line; values at most 1 place each score on its own line
	 * @param bb Nonnull destination, whose existing contents are retained
	 * @return The supplied builder
	 */
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
				//Select byte append; append(int) would print decimal digits, not a quality character.
				bb.append((byte)(b+quals[i]));
			}
		}
		return bb;
	}

	/**
	 * Appends len copies of a synthetic quality score using the current quality-output mode.
	 * Numeric mode uses spaces and wrap-boundary newlines; ASCII mode narrows the sum
	 * of q and FASTQ.ASCII_OFFSET_OUT to a byte. No trailing newline is appended.
	 * @param q Synthetic score to repeat
	 * @param len Number of scores to append
	 * @param wrap Numeric scores per line; ignored in ASCII mode
	 * @param bb Nonnull destination, whose existing contents are retained
	 * @return The supplied builder
	 */
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

	/** @return The name captured from the output descriptor, possibly null. */
	public String fname(){return fname;}
	/** @return The current read counter maintained by the concrete writer. */
	public long readsWritten(){return readsWritten;}
	/** @return The current base counter maintained by the concrete writer. */
	public long basesWritten(){return basesWritten;}

	/** @return The current error flag; this accessor does not wait for completion. */
	public final boolean errorState(){return errorState;}
	/** @return The completion flag maintained by the subclass; does not wait for completion. */
	public final boolean finishedSuccessfully(){return finishedSuccessfully;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** True if an error was encountered while writing. */
	protected boolean errorState=false;
	/** True if writing completed without errors. */
	protected boolean finishedSuccessfully=false;

	/** True for either SAM or BAM output. */
	public final boolean OUTPUT_SAM;
	/** True if the writer is emitting BAM format. */
	public final boolean OUTPUT_BAM;
	/** True if the writer is emitting FASTQ format. */
	public final boolean OUTPUT_FASTQ;
	/** True if the writer is emitting FASTA format. */
	public final boolean OUTPUT_FASTA;
	/** True if the writer is emitting FASTR format. */
	public final boolean OUTPUT_FASTR;
	/** True for the header-only output format selected by FileFormat.header(). */
	public final boolean OUTPUT_HEADER;
	/** True for attachment output selected by FileFormat.attachment(). */
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
	/** Legacy primary output, or null when absent or delegated to the SAM/BAM adapter. */
	protected final OutputStream myOutstream;
	/** Output stream for quality data when separate files are used. */
	protected final OutputStream myQOutstream;
	/** Thread-safe queue of write jobs. */
	protected final ArrayBlockingQueue<Job> queue;

	/** Count of reads written to output so far. */
	protected long readsWritten=0;
	/** Count of bases written to output so far. */
	protected long basesWritten=0;
	/** Next submission ID; checked for ListNum input and assigned for unnumbered jobs. */
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
	/** Suppress sequence records in legacy shared/generated SAM headers. */
	public static boolean NO_HEADER_SEQUENCES=false;
	/** Use attached SamLine data when emitting SAM records. */
	public static boolean USE_ATTACHED_SAMLINE=false;

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	//TODO: Should be replaced with ListNum
	/**
	 * Write job holding a list reference, ordering ID and independent close/poison flags.
	 * Marker factories allocate empty jobs; the consumer determines how flags affect output.
	 */
	protected static class Job implements HasID{

		/** Captures the list reference, close flag, poison flag and ID without copying records. */
		public Job(ArrayList<Read> list_, boolean closeWhenDone_,
				boolean poisonThread_, long id_){
			list=list_;
			close=closeWhenDone_;
			poison=poisonThread_;
			id=id_;
		}

		/*--------------------------------------------------------------*/

		/** Returns whether the list is null or empty, independently of terminal flags. */
		public boolean isEmpty(){return list==null || list.isEmpty();}
		/** Submitted list reference, or null for an empty control marker. */
		public final ArrayList<Read> list;
		/** Close/last flag exposed through last(). */
		public final boolean close;
		/** Poison flag exposed through poison(). */
		public final boolean poison;
		/** Ordering ID supplied at construction. */
		public final long id;
		@Override
		public long id(){return id;}
		@Override
		public boolean poison(){return poison;}
		@Override
		public boolean last(){return close;}
		/** Creates a new empty poison marker with the supplied ID and close=false. */
		@Override
		public Job makePoison(long id_){return new Job(null, false, true, id_);}
		/** Creates a new empty close marker with the supplied ID and poison=false. */
		@Override
		public Job makeLast(long id_){return new Job(null, true, false, id_);}

	}

}
