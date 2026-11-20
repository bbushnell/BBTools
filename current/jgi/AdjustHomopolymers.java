package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Expands or contracts exact byte-equal homopolymer runs in sequencing reads.
 * Advertised sequence inputs are single-end or paired-end FASTA/FASTQ. Editing
 * groups the bytes supplied by input/validation, including case in equality.
 * Defined bases are A/C/G/T/U in either case; undefined runs are retained.
 * The final run is adjusted like any other run.
 * <p>The required rate is finite, nonzero and at least -1. Rate {@code -1}
 * removes each defined run exactly; other rates use float multiplication and
 * int truncation toward zero. Expansion copies the run's final original quality;
 * contraction drops qualities with their bases. Null qualities remain null.
 * <p>One caller edits reads using reused instance scratch; stream wrappers manage
 * background I/O. Instances are not safe for concurrent editing. Construction
 * changes shared BBTools settings and does not restore the caller's configuration.
 * <p>Reads need nonnull sequence bases. Missing SAM sequence is not repaired, and
 * editing does not recalculate alignment metadata or guarantee mapped SAM/BAM
 * consistency. Arithmetic bounds do not guarantee allocation success for huge reads.
 * Captured streams receive cleanup attempts on processing failure; factory-internal
 * acquisition, partial startup and blocked underlying I/O retain their own limits.
 * @author Brian Bushnell
 * @date June 20, 2014
 */
public class AdjustHomopolymers{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructs the driver, processes reads and finally closes redirected diagnostics.
	 * The diagnostic finally block applies after construction succeeds; standard
	 * output and error remain open through {@link Shared#closeStream(PrintStream)}.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		//Include construction in the reported elapsed time.
		Timer t=new Timer();

		AdjustHomopolymers x=new AdjustHomopolymers(args);

		Throwable failure=null;
		try{x.process(t);}
		catch(RuntimeException | Error e){failure=e; throw e;}
		finally{
			try{Shared.closeStream(x.outstream);}
			catch(RuntimeException | Error e){
				if(failure==null){throw e;}
				retainFailure(failure, e);
			}
		}
	}

	/**
	 * Parses arguments, configures shared settings and prepares input/output formats.
	 * Validates the final rate after parsing, applies paired-path/interleaving rules
	 * and performs file preflight checks. Processing streams are created by
	 * {@link #process(Timer)}. Shared settings are not rolled back on return or failure.
	 * @param args Command-line configuration
	 * @throws IllegalArgumentException If a rate value is missing or malformed, or
	 * the final rate is nonfinite, zero or below -1
	 */
	public AdjustHomopolymers(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		Shared.capBuffers(4); //Only for singlethreaded programs

		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();

			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;

			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;
			extin=parser.extin;

			out1=parser.out1;
			out2=parser.out2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			extout=parser.extout;
		}

		if(!Float.isFinite(rate) || rate==0 || rate<-1){
			throw new IllegalArgumentException("rate is required and must be finite, nonzero, and at least -1; observed "+rate);
		}

		doPoundReplacement(); //Replace # with 1 and 2
		adjustInterleaving(); //Make sure interleaving agrees with number of input and output files
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses standard BBTools options plus this instance's rate and verbosity.
	 * Splits only the first equals sign. Final rate-domain and required-path checks
	 * belong to the constructor. Unknown options are printed and use the existing
	 * assertion policy; disabling assertions does not turn them into exceptions.
	 * @param args Command-line options
	 * @return Parser holding the standard options
	 * @throws IllegalArgumentException If a rate value is missing or malformed
	 */
	private Parser parse(String[] args){

		Parser parser=new Parser();

		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Keep any later equals signs in the value.
			String[] split=Parse.splitOnFirst(arg, '=');
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("rate")){
				if(b==null){throw new IllegalArgumentException("Missing value for rate; specify a finite nonzero value at least -1.");}
				rate=Float.parseFloat(b);
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		return parser;
	}

	/** Expands {@code #} paired paths to {@code 1}/{@code 2} when no mate path is set.
	 * Input expansion also requires that the literal primary path not exist.
	 * Requires a primary input and forbids a secondary output without a primary. */
	private void doPoundReplacement(){
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}

		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}

		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}

		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
	}

	/**
	 * Tries existing filename alternatives with {@code .gz}/{@code .bz2} added or
	 * removed for input and quality-input paths. Delegates to {@link Tools#fixExtension(String)},
	 * which respects the shared enable flag and leaves existing paths/stdin alone.
	 * This is filename lookup, not inspection of compression contents.
	 */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
		qfin1=Tools.fixExtension(qfin1);
		qfin2=Tools.fixExtension(qfin2);
	}

	/** Preflights readable inputs, writable outputs and duplicate specifications.
	 * Uses the existing Tools checks and overwrite/append flags; this is not an
	 * atomic guarantee about later opens or every filesystem alias. */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}

	/** Configures shared FASTQ interleaving from file layout and explicit options.
	 * Two inputs force interleaving off. Otherwise, two outputs imply interleaved
	 * input only when the interleaving option was not explicitly set. */
	private void adjustInterleaving(){
		//Adjust interleaved detection based on the number of input files
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}

		//Adjust interleaved settings based on number of output files
		if(!setInterleaved){
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}
	}

	/** Selects ByteFile mode 2 when neither mode is forced and threads exceed two,
	 * then asserts that the shared FASTA settings are consistent. */
	private static void checkStatics(){
		//Choose a default input implementation without overriding an explicit mode.
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}

		assert(FastaReadInputStream.settingsOK());
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Runs serial read editing with concurrent stream wrappers.
	 * Once the input factory returns, captures wrappers before their outer start
	 * calls and attempts independent cleanup in finally. Unchecked processing
	 * failures abort queued output and remain primary if cleanup also fails.
	 * Factory-internal acquisition, partial startup and blocked I/O are not repaired
	 * here. After normal processing and cleanup, prints counters/timing, then throws
	 * if aggregated statistics or stream error flags were set.
	 * @param t Running timer to stop and report after normal processing/cleanup
	 */
	void process(Timer t){

		final ConcurrentReadInputStream cris=makeCris();
		ConcurrentReadOutputStream ros=null;
		boolean outputStarted=false;
		Throwable failure=null;
		try{
			cris.start();
			if(verbose){outstream.println("Started cris");}
			final boolean paired=cris.paired();
			if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
			ros=makeCros(paired);
			if(ros!=null){ros.start(); outputStarted=true;}
			readsProcessed=readsOut=0;
			basesProcessed=basesOut=0;
			processInner(cris, ros);
			if(verbose){outstream.println("Finished; closing streams.");}
			errorState|=ReadStats.writeAll();
		}catch(RuntimeException | Error e){
			failure=e;
			throw e;
		}finally{
			finishProcessingStreams(cris, ros, outputStarted, failure);
		}

		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));

		//Reject reported statistics or stream failures after printing their totals.
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates the input wrapper; process owns its outer start and close.
	 * @return Unstarted concurrent input wrapper; factories may acquire sources internally */
	private ConcurrentReadInputStream makeCris(){
		return ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
	}

	/** Creates optional output, captured by process before its outer start call.
	 * Descriptors select the writers; factories may acquire/start internal resources.
	 * @param pairedInput Whether input is paired; used only to decide the interleaved-output status message
	 * @return Output wrapper, or null when no primary output is configured */
	private ConcurrentReadOutputStream makeCros(boolean pairedInput){
		if(ffout1==null){return null;}

		//Use four output buffers.
		final int buff=4;

		//Notify user of output mode
		if(pairedInput && out2==null && (in1!=null && !ffin1.samOrBam() && !ffout1.samOrBam())){
			outstream.println("Writing interleaved.");
		}

		return ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, false);
	}

	/**
	 * Borrows batches, checks pairing and delegates editing/submission.
	 * Owns one return attempt per borrowed batch, including the final empty batch or
	 * a batch left by failure. Clears the pending reference before attempting return,
	 * so a throwing return hook is not retried. A return failure is suppressed onto
	 * an existing processing failure. Does not close either borrowed stream.
	 * @param cris Input stream supplying batches
	 * @param ros Output stream for processed reads, or null to omit submission
	 */
	void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		ListNum<Read> ln=null;
		Throwable failure=null;
		try{
			ln=cris.nextList();
			if(ln!=null && !ln.isEmpty()){
				Read r=ln.get(0);
				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()) :
					"Non-SAM input mate links must agree with the stream's paired mode before processing.";
			}
			while(ln!=null && ln.size()>0){
				processList(ln, ros);
				final ListNum<Read> spent=ln;
				ln=null;
				cris.returnList(spent);
				ln=cris.nextList();
			}
		}catch(RuntimeException | Error e){
			failure=e;
			throw e;
		}finally{
			final ListNum<Read> pending=ln;
			ln=null;
			try{cris.returnList(pending);}
			catch(RuntimeException | Error e){
				if(failure==null){throw e;}
				retainFailure(failure, e);
			}
		}
	}

	/**
	 * Validates and edits the roots/mates in a borrowed batch.
	 * Updates logical pre-edit and retained post-edit counters before optional output
	 * submission; these counters do not prove successful writes. processInner owns
	 * returning the batch on success or failure. Does not close the streams.
	 * @param ln Borrowed nonempty read batch and its output ID
	 * @param ros Output stream, or null to count/edit without submission
	 */
	void processList(ListNum<Read> ln, final ConcurrentReadOutputStream ros){

		//Grab the actual read list from the ListNum
		final ArrayList<Read> reads=ln.list;

		//Loop through each read in the list
		for(int idx=0; idx<reads.size(); idx++){
			final Read r1=reads.get(idx);
			final Read r2=r1.mate;

			//Validate on the editing caller before using the reads.
			if(!r1.validated()){r1.validate(true);}
			if(r2!=null && !r2.validated()){r2.validate(true);}

			//Track the initial length for statistics
			final int initialLength1=r1.length();
			final int initialLength2=r1.mateLength();

			//Increment counters
			readsProcessed+=r1.pairCount();
			basesProcessed+=initialLength1+initialLength2;

			{
				boolean keep=processReadPair(r1, r2);

				if(!keep){reads.set(idx, null);}
				else{
					readsOut+=r1.pairCount();
					basesOut+=r1.pairLength();
				}
			}
		}

		//Output reads to the output stream
		if(ros!=null){ros.add(reads, ln.id);}

	}

	/** Attempts independent input/output shutdown and retains the first unchecked failure.
	 * Aborts output before input closure after processing failure; joins only after
	 * its outer start returned successfully. Reported stream errors enter errorState.
	 * Factories, partially started NEW backends and blocked I/O retain their own limits.
	 * @param cris Captured input wrapper to close
	 * @param ros Captured output wrapper, or null
	 * @param outputStarted Whether the output's outer start returned successfully
	 * @param primary Existing unchecked failure, or null for normal processing
	 * @throws RuntimeException If cleanup fails without an existing primary failure
	 * @throws Error If cleanup raises an Error without an existing primary failure
	 */
	private void finishProcessingStreams(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros,
			final boolean outputStarted, final Throwable primary){
		Throwable failure=primary;
		if(ros!=null && failure!=null){
			try{ros.abort();}
			catch(RuntimeException | Error e){failure=retainFailure(failure, e);}
		}
		try{errorState|=ReadWrite.closeStream(cris);}
		catch(RuntimeException | Error e){failure=retainFailure(failure, e);}
		if(ros!=null){
			try{
				if(failure==null && outputStarted){ros.close();}
				else{ros.abort();}
			}catch(RuntimeException | Error e){
				failure=retainFailure(failure, e);
				try{ros.abort();}
				catch(RuntimeException | Error secondary){failure=retainFailure(failure, secondary);}
			}
			if(outputStarted){
				try{
					ros.join();
					errorState|=ros.errorState() || !ros.finishedSuccessfully();
				}catch(RuntimeException | Error e){failure=retainFailure(failure, e);}
			}else{errorState=true;}
		}
		if(primary==null && failure!=null){
			if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
			throw (Error)failure;
		}
	}

	/** Retains the first failure, suppressing a different secondary object.
	 * @param primary Earlier failure, or null
	 * @param secondary Nonnull cleanup failure
	 * @return Primary when present, otherwise secondary; never self-suppresses */
	private static Throwable retainFailure(final Throwable primary, final Throwable secondary){
		if(primary==null){return secondary;}
		if(primary!=secondary){primary.addSuppressed(secondary);}
		return primary;
	}

	/**
	 * Edits each nonnull read and retains the pair on normal return.
	 * @param r1 First read, or null to skip
	 * @param r2 Mate, or null for single-end input
	 * @return true; this implementation never filters a normally processed pair
	 */
	boolean processReadPair(final Read r1, final Read r2){
		processRead(r1);
		processRead(r2);
		return true;
	}

	/**
	 * Adjusts each completed run, including the final run, using exact byte equality.
	 * Fully defined runs are removed exactly at rate=-1; other rates retain float
	 * multiplication and truncation toward zero. Undefined runs remain unchanged.
	 * Run ownership and builder capacity are checked before applying each adjustment.
	 * Added copies inherit the run's final original quality, and contraction removes
	 * bases/qualities from the run's end.
	 * Replaces base/quality arrays only when length changes. Null qualities stay
	 * null; identifiers, mate links and alignment metadata are not recalculated.
	 * Uses shared instance scratch and is not safe for concurrent calls. Rate-domain
	 * validation belongs to construction; a nonnull read must contain sequence bases
	 * and any nonnull quality array must cover those bases.
	 * @param r Caller-owned read to edit; null returns immediately
	 * @throws IllegalArgumentException If an adjustment erases a prior run or exceeds builder capacity
	 * @throws NullPointerException If a nonnull read has null bases
	 */
	void processRead(final Read r){
		if(r==null){return;}
		bbBases.clear();
		bbQuals.clear();
		final byte[] bases=r.bases;
		final byte[] quals=(r.quality==null ? fakeQuality(bases.length) : r.quality);

		int streak=0;
		for(int i=0; i<bases.length; i++){
			final byte b=bases[i];
			final byte q=quals[i];
			bbBases.append(b);
			bbQuals.append(q);
			streak++;
			if(i+1==bases.length || bases[i+1]!=b){
				//Validate the completed run's adjustment before changing either builder (JGI-052).
				int adjustment=AminoAcid.isFullyDefined(b) ? adjustmentForRun(rate, streak, bbBases.length) : 0;
				if(adjustment<0){
					bbBases.length+=adjustment;
					bbQuals.length+=adjustment;
				}else{
					for(; adjustment>0; adjustment--){
						bbBases.append(b);
						bbQuals.append(q);
					}
				}
				streak=0;
			}
		}
		assert(bbBases.length()==bbQuals.length()) : "Run editing must retain one scratch quality per base before replacing Read arrays: "+
			bbBases.length()+" bases, "+bbQuals.length()+" qualities";
		if(bbBases.length()!=r.length()){
			r.bases=bbBases.toBytes();
			r.quality=(r.quality==null ? null : bbQuals.toBytes());
		}
	}

	/**
	 * Computes a bounded adjustment without editing a read or allocating scratch storage.
	 * Rate=-1 deletes the run exactly; other rates retain float multiplication and int
	 * truncation. The adjustment may not erase prior output or exceed ByteBuilder.MAXLEN.
	 * Passing these arithmetic checks does not guarantee later allocation success.
	 * @param rate Configured finite nonzero rate, at least -1
	 * @param streak Positive length of the completed defined-base run
	 * @param currentLength Scratch length including prior output and the original run,
	 * at least streak
	 * @return Checked adjustment to apply equally to base and quality builders
	 * @throws IllegalArgumentException If the adjustment exceeds its run or builder capacity
	 */
	private static int adjustmentForRun(final float rate, final int streak, final int currentLength){
		assert(streak>0 && currentLength>=streak) :
			"processRead appends the completed run before adjustment, so scratch must contain it: streak="+streak+", currentLength="+currentLength;
		final int adjustment=(rate==-1f ? -streak : (int)(rate*streak));
		final long projected=(long)currentLength+adjustment;
		if(adjustment<-streak || projected>ByteBuilder.MAXLEN){
			throw new IllegalArgumentException("Homopolymer adjustment exceeds its run or builder capacity: rate="+rate+
				", streak="+streak+", currentLength="+currentLength+", adjustment="+adjustment+", projectedLength="+projected);
		}
		return adjustment;
	}

	/** Supplies cached synthetic qualities for null-quality read processing.
	 * Grows with ten spare entries and fills with Shared.FAKE_QUAL on allocation.
	 * Only the requested prefix is read; this does not give the Read a quality array.
	 * @param minlen Required number of usable entries
	 * @return Borrowed instance buffer of at least minlen entries; do not mutate */
	private byte[] fakeQuality(int minlen){
		if(fakeQuality.length<minlen){
			fakeQuality=KillSwitch.allocByte1D(minlen+10);
			Arrays.fill(fakeQuality, Shared.FAKE_QUAL);
		}
		return fakeQuality;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input path after paired-path/extension handling. */
	private String in1=null;
	/** Secondary input file path for paired-end data. */
	private String in2=null;

	/** Optional separate primary quality-input path. */
	private String qfin1=null;
	/** Optional separate mate quality-input path. */
	private String qfin2=null;

	/** Optional primary output path; null disables read submission. */
	private String out1=null;
	/** Secondary output file path for paired-end data. */
	private String out2=null;

	/** Optional separate primary quality-output path. */
	private String qfout1=null;
	/** Optional separate mate quality-output path. */
	private String qfout2=null;

	/** Optional input format extension override passed to FileFormat. */
	private String extin=null;
	/** Optional output format extension override passed to FileFormat. */
	private String extout=null;

	/** Whether interleaved processing was explicitly configured. */
	private boolean setInterleaved=false;

	/** Reused base scratch for the current read; owned by the single editing caller. */
	private ByteBuilder bbBases=new ByteBuilder();
	/** Parallel quality scratch, including synthetic entries when input quality is null. */
	private ByteBuilder bbQuals=new ByteBuilder();

	/** Required adjustment rate; construction validates the final parsed value. */
	private float rate=0.0f;

	/** Reused synthetic-quality buffer; does not replace a Read's null quality. */
	private byte[] fakeQuality=new byte[0];

	/*--------------------------------------------------------------*/

	/** Reads encountered before editing, including mates. */
	protected long readsProcessed=0;
	/** Bases encountered before editing, including mate bases. */
	protected long basesProcessed=0;

	/** Logical retained reads including mates, counted even without output; not write confirmation. */
	protected long readsOut=0;
	/** Logical retained post-edit bases including mates, counted before optional submission. */
	protected long basesOut=0;

	/** Limit forwarded to CRIS; normally fragments, with reader-specific units; -1 is unlimited. */
	private long maxReads=-1;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input descriptor. */
	private final FileFormat ffin1;
	/** Optional secondary input descriptor. */
	private final FileFormat ffin2;

	/** Primary output descriptor, or null when output is omitted. */
	private final FileFormat ffout1;
	/** Optional secondary output descriptor. */
	private final FileFormat ffout2;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Diagnostic stream selected by PreParser; main closes redirects after construction. */
	private PrintStream outstream=System.err;
	/** Process-wide verbose logging switch. */
	public static boolean verbose=false;
	/** Aggregated statistics/stream error flags; thrown failures are a separate error channel. */
	public boolean errorState=false;
	/** Whether existing outputs may be overwritten. */
	private boolean overwrite=true;
	/** Whether output is appended instead of overwritten. */
	private boolean append=false;
	/** Ordering-table request passed to output descriptors; fixed false in this driver. */
	private final boolean ordered=false;

}
