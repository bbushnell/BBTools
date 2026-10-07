package stream;

import java.util.ArrayList;
import java.util.Arrays;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import structures.ListNum;

/**
 * Common legacy base for numbered, recyclable Read batches and format-dispatch factories.
 * Concrete subclasses supply processing, pairing, sampling and resource/completion
 * behavior. This base supplies default thread startup, ListNum return adaptation,
 * and a convenience drain into memory. It does not impose a universal shutdown,
 * restart or final-counter guarantee. Newer I/O uses the separate Streamer API.
 *
 * Local factory construction can open format-specific inputs but does not start the
 * returned CRIS. MPI branches are unfinished scaffolding, not working distributed
 * transport: the true-MPI branch retains its deliberate assertion fence, and the
 * current receiving-rank wrapper cannot accept the null local source passed here.
 *
 * @author Brian Bushnell
 * @date Nov 26, 2014
 */
public abstract class ConcurrentReadInputStream implements ConcurrentReadStreamInterface{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Retains the input identifier; inherited buffer settings capture current Shared values.
	 * @param fname_ Input path or implementation-specific identifier */
	protected ConcurrentReadInputStream(String fname_){fname=fname_;}

	/** Resolves positional path arguments, normalizing case-insensitive "null" in place.
	 * Only the first four entries select files; extra entries are otherwise ignored.
	 * Uses current shared MPI defaults and requires a nonnull primary path.
	 * @param maxReads Limit forwarded unchanged; concrete readers determine exact units
	 * @param keepSamHeader Request header retention on the SAM/BAM path
	 * @param allowSubprocess Permit subprocesses during descriptor resolution
	 * @param args Mutable paths: primary, optional mate, optional QUAL1 and QUAL2
	 * @return Selected reader; see the descriptor overload for MPI limitations */
	protected static ConcurrentReadInputStream getReadInputStream(long maxReads, boolean keepSamHeader, boolean allowSubprocess, String...args){
		assert(args.length>0) : Arrays.toString(args);
		for(int i=0; i<args.length; i++){
			if("null".equalsIgnoreCase(args[i])){args[i]=null;}
		}
		assert(args[0]!=null) : Arrays.toString(args);

		assert(args.length<2 || !args[0].equalsIgnoreCase(args[1]));
		String in1=args[0], in2=null, qf1=null, qf2=null;
		if(args.length>1){in2=args[1];}
		if(args.length>2){qf1=args[2];}
		if(args.length>3){qf2=args[3];}

		final FileFormat ff1=FileFormat.testInput(in1, null, allowSubprocess);
		final FileFormat ff2=FileFormat.testInput(in2, null, allowSubprocess);

		return getReadInputStream(maxReads, keepSamHeader, ff1, ff2, qf1, qf2);
	}

	/**
	 * Selects a reader using shared MPI defaults and no separate quality files.
	 *
	 * @param maxReads Limit forwarded unchanged; normally fragments (reads or pairs)
	 * @param keepSamHeader Request shared header retention for SAM/BAM input
	 * @param ff1 Primary input file format (required)
	 * @param ff2 Secondary input file format (optional, for paired reads)
	 * @return Appropriate ConcurrentReadInputStream implementation
	 */
	public static ConcurrentReadInputStream getReadInputStream(long maxReads, boolean keepSamHeader, FileFormat ff1, FileFormat ff2){
		return getReadInputStream(maxReads, keepSamHeader, ff1, ff2, (String)null, (String)null, Shared.USE_MPI, Shared.MPI_KEEP_ALL);
	}

	/**
	 * Selects a reader with explicit MPI settings and no separate quality files.
	 * MPI scaffolding limitations are described on the fully configured overload.
	 *
	 * @param maxReads Limit forwarded unchanged; normally fragments (reads or pairs)
	 * @param keepSamHeader Request shared header retention for SAM/BAM input
	 * @param ff1 Primary input file format (required)
	 * @param ff2 Secondary input file format (optional, for paired reads)
	 * @param mpi Request the unfinished distributed-wrapper path
	 * @param keepAll Delivery policy forwarded to the distributed wrapper
	 * @return Appropriate ConcurrentReadInputStream implementation
	 */
	public static ConcurrentReadInputStream getReadInputStream(long maxReads, boolean keepSamHeader, FileFormat ff1, FileFormat ff2,
			final boolean mpi, final boolean keepAll){
		return getReadInputStream(maxReads, keepSamHeader, ff1, ff2, (String)null, (String)null, mpi, keepAll);
	}

//	/** @See primary method */
//	public static ConcurrentReadInputStream getReadInputStream(long maxReads, boolean keepSamHeader, FileFormat ff1, String qf1){
//		return getReadInputStream(maxReads, keepSamHeader, ff1, null, qf1, null, Shared.USE_MPI, Shared.MPI_KEEP_ALL);
//	}

	/**
	 * Selects a reader with optional quality files and shared MPI defaults.
	 * Quality paths are used only on the non-shredding FASTA path.
	 *
	 * @param maxReads Limit forwarded unchanged; normally fragments (reads or pairs)
	 * @param keepSamHeader Request shared header retention for SAM/BAM input
	 * @param ff1 Primary input file format (required)
	 * @param ff2 Secondary input file format (optional, for paired reads)
	 * @param qf1 Primary quality file path (optional)
	 * @param qf2 Secondary quality file path (optional)
	 * @return Appropriate ConcurrentReadInputStream implementation
	 */
	public static ConcurrentReadInputStream getReadInputStream(long maxReads, boolean keepSamHeader,
			FileFormat ff1, FileFormat ff2, String qf1, String qf2){
		return getReadInputStream(maxReads, keepSamHeader, ff1, ff2, qf1, qf2, Shared.USE_MPI, Shared.MPI_KEEP_ALL);
	}

	/**
	 * Selects legacy inputs from the primary descriptor and wraps them for concurrent delivery.
	 * FASTQ, ONELINE, FASTA, SCARF and HEADER construct matching reader types for an
	 * optional second file. Non-shredding FASTA uses qf1/qf2 when supplied; shredding
	 * ignores them. SAM/BAM, GFA, EMBL and GenBank construct only the primary input.
	 * Shared name assertions still inspect supplied descriptors/quality paths before dispatch.
	 * BREAD passes both descriptors to RTextInputStream and uses the legacy wrapper;
	 * synthetic sequential/random paths also bypass ordinary file-pair construction.
	 * CSFASTA is rejected. Format routing is based on ff1, not independently on ff2.
	 *
	 * With mpi=false, returns an unstarted wrapper; input constructors may open files.
	 * With mpi=true, rank zero first constructs and starts a local reader. USE_CRISMPI
	 * then hits the deliberate assertion fence (returns null with assertions disabled).
	 * The alternate distributed wrapper is unfinished and currently fails for a
	 * receiving rank's null source. Neither branch provides a complete distributed reader.
	 *
	 * @param maxReads Limit forwarded unchanged; normally fragments, with reader-specific counting
	 * @param keepSamHeader Request shared header retention on the SAM/BAM path only
	 * @param ff1 Primary input file format (required)
	 * @param ff2 Secondary input file format (optional, for paired reads)
	 * @param qf1 Primary quality file path (optional, for FASTA+QUAL)
	 * @param qf2 Secondary quality file path (optional, for FASTA+QUAL)
	 * @param mpi Request the unfinished distributed-wrapper path
	 * @param keepAll Delivery policy forwarded to the distributed wrapper
	 * @return Appropriate ConcurrentReadInputStream implementation
	 * @throws RuntimeException If the primary format is unsupported
	 */
	public static ConcurrentReadInputStream getReadInputStream(long maxReads, boolean keepSamHeader,
			FileFormat ff1, FileFormat ff2, String qf1, String qf2, final boolean mpi, final boolean keepAll){
		//MPI scaffolding: rank0 starts a local reader; other ranks pass a null source.
		//ConcurrentReadInputStreamD transport is unfinished; non-MPI uses the dispatch below.
		if(mpi){
			final int rank=Shared.MPI_RANK;
			final ConcurrentReadInputStream cris0;
			if(rank==0){
				cris0=getReadInputStream(maxReads, keepSamHeader, ff1, ff2, qf1, qf2, false, true);
				cris0.start();
			}else{
				cris0=null;
			}
			final ConcurrentReadInputStream crisD;
			if(Shared.USE_CRISMPI){
				//Deliberate fence (NOT a bug): the true-MPI cris is unfinished/disabled. The message IS the acceptance test - validate before uncommenting the line below; do not delete it as "obsolete".
				assert(false) : "To support MPI, uncomment this.";
//				crisD=new ConcurrentReadInputStreamMPI(cris0, rank==0, keepAll);
				crisD=null;
			}else{
				crisD=new ConcurrentReadInputStreamD(cris0, rank==0, keepAll);
			}
			return crisD;
		}

		assert(ff1!=null);
		assert(ff2==null || ff1.name()==null || !ff1.name().equalsIgnoreCase(ff2.name())) : ff1.name()+", "+ff2.name();
		//TODO: Probable guard typo - this qf1-gated check compares the primary name to qf2,
		//leaving a primary/qf1 collision unchecked here. Review the duplicate-file matrix separately.
		assert(qf1==null || ff1.name()==null || !ff1.name().equalsIgnoreCase(qf2));
		assert(qf1==null || qf2==null || !qf1.equalsIgnoreCase(qf2)); //#001-fix [stream/ConcurrentReadInputStream#001]: was missing '!' -> asserted qf1==qf2, crashing (-ea) on valid distinct paired qual files and passing the same-file error. Now asserts the two qual files DIFFER, matching the sibling asserts above; the same-file case is already blocked upstream by testForDuplicateFiles. Verified via reformat.sh paired FASTA+QUAL before/after, 2026-06-18.

		//Format dispatch: pick the format-specific ReadInputStream for ff1 (and ff2 if paired), then wrap it in the workhorse ConcurrentGenericReadInputStream (or ConcurrentLegacyReadInputStream for the older bread/sequential paths). Branches are mutually exclusive by FileFormat type; an unrecognized format crashes loud (final else).
		final ConcurrentReadInputStream cris;

		if(ff1.fastq()){

			ReadInputStream ris1, ris2;

			ris1=new FastqReadInputStream(ff1);
			try{
				ris2=(ff2==null ? null : new FastqReadInputStream(ff2));
			}catch(AssertionError e){//Handles problems with quality score autodetection
				ris1.close();
				throw e;
			}
			cris=new ConcurrentGenericReadInputStream(ris1, ris2, maxReads);

		}else if(ff1.oneline()){

			ReadInputStream ris1=new OnelineReadInputStream(ff1);
			ReadInputStream ris2=(ff2==null ? null : new OnelineReadInputStream(ff2));
			cris=new ConcurrentGenericReadInputStream(ris1, ris2, maxReads);

		}else if(ff1.fasta()){
			boolean amino=(ff1.amino() || Shared.AMINO_IN);

			ReadInputStream ris1;
			ReadInputStream ris2;
			if(ff1.preferShreds()){
				ris1=new FastaShredInputStream(ff1, amino, ff2==null ? Shared.bufferData() : -1);
				ris2=(ff2==null ? null : new FastaShredInputStream(ff2, amino, -1));
			}else{
				ris1=(qf1==null ? new FastaReadInputStream(ff1, (FASTQ.FORCE_INTERLEAVED && ff2==null), amino, ff2==null ? Shared.bufferData() : -1)
						: new FastaQualReadInputStream(ff1, qf1));
				ris2=(ff2==null ? null : qf2==null ? new FastaReadInputStream(ff2, false, amino, -1) : new FastaQualReadInputStream(ff2, qf2));
			}
			cris=new ConcurrentGenericReadInputStream(ris1, ris2, maxReads);

//			cris.start();
//			ListNum<Read> ln=cris.nextList();
//			System.out.println(ln);
//
//			assert(false) : ff1+", "+ff2;
		}else if(ff1.scarf()){

			ReadInputStream ris1=new ScarfReadInputStream(ff1);
			ReadInputStream ris2=(ff2==null ? null : new ScarfReadInputStream(ff2));
			cris=new ConcurrentGenericReadInputStream(ris1, ris2, maxReads);

		}else if(ff1.samOrBam()){
			int threads=Tools.mid(1, SamStreamer.DEFAULT_THREADS, Shared.threads());
			ReadInputStream ris1=new SamReadInputStream(ff1, keepSamHeader, threads, maxReads);
			cris=new ConcurrentGenericReadInputStream(ris1, null, maxReads);
			assert(!cris.paired()) : "\nff1="+ff1+"\nff2="+ff2+
				"\nris1="+ris1+"\nris2"+null+"\nris1paired"+ris1.paired()+
				"\np1="+cris.producers()[0]+"\np2=null"; //#002-fix [stream/ConcurrentReadInputStream#002]: was cris.producers()[2] - producers() is length-1 here (single-producer SAM cris), so the old [2] threw AIOOBE on assert-failure instead of printing the diagnostic; producer2 is null on this path.
		}else if(ff1.bread()){
//			assert(false) : ff1;
			RTextInputStream rtis=new RTextInputStream(ff1, ff2, maxReads);
			cris=new ConcurrentLegacyReadInputStream(rtis, maxReads);// TODO: Change to generic

		}else if(ff1.header()){

			HeaderInputStream ris1=new HeaderInputStream(ff1);
			HeaderInputStream ris2=(ff2==null ? null : new HeaderInputStream(ff2));
			cris=new ConcurrentGenericReadInputStream(ris1, ris2, maxReads);

		}else if(ff1.gfa()){

			GfaReadInputStream ris1=new GfaReadInputStream(ff1);
			cris=new ConcurrentGenericReadInputStream(ris1, null, maxReads);

		}else if(ff1.sequential()){

			SequentialReadInputStream ris=new SequentialReadInputStream(maxReads, 200, 50, 0, false);
			cris=new ConcurrentLegacyReadInputStream(ris, maxReads);

		}else if(ff1.csfasta()){

			throw new RuntimeException("csfasta is no longer supported.");

		}else if(ff1.random()){

			RandomReadInputStream3 ris=new RandomReadInputStream3(maxReads, FASTQ.FORCE_INTERLEAVED);
			cris=new ConcurrentGenericReadInputStream(ris, null, maxReads);

		}else if(ff1.embl()){

			EmblReadInputStream ris=new EmblReadInputStream(ff1);
			cris=new ConcurrentGenericReadInputStream(ris, null, maxReads);

		}else if(ff1.gbk()){

			GbkReadInputStream ris=new GbkReadInputStream(ff1);
			cris=new ConcurrentGenericReadInputStream(ris, null, maxReads);

		}else{
			cris=null;
			throw new RuntimeException(""+ff1);
		}

		return cris;
	}


	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructs, starts and drains a selected stream into one in-memory list.
	 * Uses shared MPI defaults; all factory limitations and concrete ownership rules apply.
	 * @param maxReads Limit forwarded to the selected reader
	 * @param keepSamHeader Request SAM/BAM shared header retention
	 * @param ff1 Required primary descriptor
	 * @param ff2 Optional secondary descriptor, used only by supported branches
	 * @param qf1 Optional primary QUAL path for non-shredding FASTA
	 * @param qf2 Optional secondary QUAL path for non-shredding FASTA
	 * @return New list containing retained Read references, without cloning the reads */
	public static ArrayList<Read> getReads(long maxReads, boolean keepSamHeader,
			FileFormat ff1, FileFormat ff2, String qf1, String qf2){
		ConcurrentReadInputStream cris=getReadInputStream(maxReads, keepSamHeader, ff1, ff2, qf1, qf2, Shared.USE_MPI, Shared.MPI_KEEP_ALL);
		cris.start();
		return cris.getReads();
	}

	/** Drains an already started stream, returning batches after copying their references.
	 * Stops on null wrapper, null payload or an empty payload; a supplied terminal
	 * wrapper is returned to the stream before closure through ReadWrite, not added
	 * to the result. An error reported by that helper invokes
	 * KillSwitch.kill rather than returning the accumulated data. No finally-close path
	 * is provided for an exception during consumption.
	 * @return New list of Read references; mate links are not flattened or cloned */
	public ArrayList<Read> getReads(){
		//Convenience drain: pull every list via nextList(), accumulate into one ArrayList, returning each (poison-returning the terminator). For callers that want all reads in memory rather than streaming. The trailing returnList after the loop handles the final size-0/null poison list that ended it.
		ListNum<Read> ln=nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		ArrayList<Read> out=new ArrayList<Read>();

		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			out.addAll(reads);
			returnList(ln.id, ln.list.isEmpty());
			ln=nextList();
			reads=(ln!=null ? ln.list : null);
		}
		if(ln!=null){
			returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
		boolean error=ReadWrite.closeStream(this);
		//[stream/ConcurrentReadInputStream getReads-error FIXED 2026-06-21] crash LOUD on a read error instead of returning PARTIAL reads
		//with only a stderr warning: a truncated/corrupt input silently yields fewer reads, and callers (e.g. aligner/AllToAll) load this as
		//REFERENCE data -> wrong results. KillSwitch.kill aborts non-zero (BBTools contract: crash, never silently wrong). Resolves the
		//AllToAll:223 "does not return the error state" TODO. Sibling of StreamerFactory.getReads#001.
		if(error){
			KillSwitch.kill("Error: a read error (corrupt or truncated input) occurred reading "+fname
				+" during getReads(); aborting rather than returning partial/incorrect read data.");
		}
		return out;
	}

	/**
	 * Default startup launches this Runnable on a new Thread, then sets started=true.
	 * Does not retain/join the thread, wait for input readiness or guard repeated calls.
	 * Concrete subclasses may override this policy; the historical collection-reader
	 * rationale below explains retaining the separate thread here.
	 */
	@Override
	public void start(){
//		System.err.println("Starting "+this);
		new Thread(this).start();// Prevents a strange deadlock in ConcurrentCollectionReadInputStream
		started=true;
	}

	/** Returns the observed startup flag, not a worker-readiness or completion barrier. */
	public final boolean started(){return started;}


	/*--------------------------------------------------------------*/
	/*----------------       Abstract Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Obtains a batch under the concrete stream's blocking/terminal/ownership contract.
	 * @return Batch or implementation-specific terminal result */
	@Override
	public abstract ListNum<Read> nextList();

	/** Forwards a nonnull wrapper's ID and isEmpty() result to the subclass return hook.
	 * Null wrappers do nothing. Does not validate IDs, clear lists or inspect marker type.
	 * @param ln Spent batch, or null to do nothing */
	@Override
	public final void returnList(ListNum<Read> ln){
		if(ln!=null){returnList(ln.id, ln.isEmpty());}
	}

	/**
	 * Returns a list container to the stream for reuse.
	 * Implementation depends on concrete subclass.
	 * @param listNum List identifier number
	 * @param poison True if this is a poison/termination signal
	 */
	@Override
	public abstract void returnList(long listNum, boolean poison);

	/** Executes the concrete stream's processing; ordinary callers use start(). */
	@Override
	public abstract void run();

	/** Requests shutdown through the concrete subclass; release and waiting semantics vary. */
	@Override
	public abstract void shutdown();

	/** Requests reuse/reset through the subclass; supported sources and startup sequencing vary. */
	@Override
	public abstract void restart();

	/** Closes through the subclass lifecycle; resource release and completion waiting vary.
	 * ReadWrite closure helpers query errorState afterward. */
	@Override
	public abstract void close();

	/**
	 * Returns true if this stream processes paired-end reads.
	 * Implementation depends on concrete subclass.
	 * @return true if stream handles paired reads
	 */
	@Override
	public abstract boolean paired();

	/** Returns the filename or identifier for this stream */
	@Override
	public String fname(){return fname;}

	/**
	 * Returns array of producer objects for this stream.
	 * Implementation depends on concrete subclass.
	 * @return Array of producer objects
	 */
	@Override
	public abstract Object[] producers();

	/**
	 * Returns true if the stream is in an error state.
	 * Implementation depends on concrete subclass.
	 * @return true if errors occurred during processing
	 */
	@Override
	public abstract boolean errorState();

	/**
	 * Sets sampling rate for subsampling reads during processing.
	 * Implementation depends on concrete subclass.
	 * @param rate Sampling rate between 0.0 and 1.0
	 * @param seed Seed interpreted by the concrete reader, including any special values
	 */
	@Override
	public abstract void setSampleRate(float rate, long seed);

	/**
	 * Reports the concrete stream's input-base counter.
	 * Counting stage, synchronization, snapshot visibility and finality depend on the subclass.
	 * @return Reported base count
	 */
	@Override
	public abstract long basesIn();

	/**
	 * Reports the concrete stream's input-read counter.
	 * Mate/fragment units, synchronization, counting stage and finality depend on the subclass.
	 * @return Reported read count in the implementation's units
	 */
	@Override
	public abstract long readsIn();

	/**
	 * Returns true if verbose output is enabled.
	 * Implementation depends on concrete subclass.
	 * @return true if verbose mode is active
	 */
	@Override
	public abstract boolean verbose();

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Per-instance snapshot of the shared batch-length setting. */
	final int BUF_LEN=Shared.bufferLen();
	/** Per-instance snapshot of the shared buffer-count setting. */
	final int NUM_BUFFS=Shared.numBuffers();
	/** Per-instance snapshot of the shared batch-data setting. */
	final long MAX_DATA=Shared.bufferData();
	/** Input name retained by the constructor. */
	public final String fname;
	/** Policy flag consulted by subclasses that support unequal paired inputs. */
	public boolean ALLOW_UNEQUAL_LENGTHS=false;
	/** Startup observation only; not volatile and not a readiness handshake. */
	boolean started=false;

	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/

	/** Shared progress-display request, interpreted by concrete readers. */
	public static boolean SHOW_PROGRESS=false;
	/** Requests elapsed seconds between progress dots in supporting readers. */
	public static boolean SHOW_PROGRESS2=false;
	/** Shared record interval for supporting progress displays. */
	public static long PROGRESS_INCR=1000000;
	/** Shared request to remove discarded reads where the concrete reader supports it. */
	public static boolean REMOVE_DISCARDED_READS=false;

}
