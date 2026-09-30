package stream;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Exercises factory-selected Streamer and Writer implementations with one consumer.
 * Reads and optionally writes sequence data, delegating format support, pairing,
 * ordering and worker selection to the factories. SAM/BAM input uses SamLines when
 * output is SAM/BAM or absent; other routes consume Reads. SAM-to-sequence output
 * requests Read conversion, while SAM-to-SAM output requests shared headers.
 * <p>
 * Limits and sampling are delegated to the reader. Skipping removes returned list
 * entries after sampling, so one skipped entry can represent a read pair. Reported
 * totals are counted here from remaining records and include mates in Read mode;
 * they are not the reader's aggregate counters. Even empty batches reach the writer
 * with their original IDs. Statistics precede the final reader error check.
 * <p>
 * This command changes shared parsing, compression and interleaving settings without
 * restoring them; it is not an isolated reusable configuration object. Its record
 * processing helpers are private, not subclass extension points.
 * <pre>
 * stream.sh in=reads.fq.gz out=sampled.fq.gz samplerate=0.1
 * stream.sh in=mapped.bam out=reads.fq.gz
 * stream.sh in1=r1.fq in2=r2.fq out=interleaved.fq
 * </pre>
 *
 * @author Brian Bushnell, Isla
 * @contributor Shinobu (documentation and formatting)
 * @date November 4, 2025
 */
public class StreamerWrapper{

	/*--------------------------------------------------------------*/
	/*----------------            Main              ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses arguments, runs the transfer, and closes a redirected status stream on normal return.
	 * @param args Command line arguments
	 */
	public static void main(final String[] args){
		final Timer t=new Timer();
		final StreamerWrapper x=new StreamerWrapper(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses settings, resolves file names and formats, and configures shared I/O flags.
	 * Does not create or start a Streamer or Writer. A primary input is required, and
	 * secondary output requires primary output. Standard parsing also changes globals.
	 * For SAM/BAM input without SAM/BAM output, selected SAM field parsing is disabled
	 * unless forceparse is set; forceparse does not restore previously disabled flags.
	 * @param args Command line arguments
	 */
	public StreamerWrapper(String[] args){
		{//Preparse block for help, config files, and outstream
			final PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();

			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;
			threadsIn=parser.threadsIn;
			threadsOut=parser.threadsOut;

			in1=parser.in1;
			in2=parser.in2;
			out1=parser.out1;
			out2=parser.out2;

			qfin1=parser.qfin1;
			qfin2=parser.qfin2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
		}

		doPoundReplacement();
		fixExtensions();
		adjustInterleaving();
		checkFileExistence();
		checkStatics();

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, true);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, append, true);

		final boolean samIn=(ffin1!=null && ffin1.samOrBam());
		final boolean samOut=(ffout1!=null && ffout1.samOrBam());
		SamLine.SET_FROM_OK=samIn;
		ReadStreamByteWriter.USE_ATTACHED_SAMLINE=samIn && samOut;
		//Determine if we need to parse SAM fields or can skip for performance
		if(!forceParse && samIn && !samOut){
			SamLine.PARSE_2=false;
			SamLine.PARSE_5=false;
			SamLine.PARSE_6=false;
			SamLine.PARSE_7=false;
			SamLine.PARSE_8=false;
			SamLine.PARSE_OPTIONAL=false;
		}
	}

	/**
	 * Parses wrapper options, standard Parser flags and the first two positional paths.
	 * Unknown options print a diagnostic and fail an assertion when assertions are enabled.
	 * @param args Command line arguments
	 * @return Parser object with standard flags processed
	 */
	private Parser parse(final String[] args){
		final Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			final String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			final String[] split=arg.split("=");
			final String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("samplerate") || a.equals("sample")){
				samplerate=Float.parseFloat(b);
			}else if(a.equals("sampleseed") || a.equals("seed")){
				sampleseed=Long.parseLong(b);
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("forceparse")){
				forceParse=Parse.parseBoolean(b);
			}else if(a.equals("skipreads")){
				skipreads=Parse.parseKMG(b);
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(i==0 && parser.in1==null && Tools.looksLikeInputSequenceStream(arg)){
				parser.in1=arg;
			}else if(i==1 && parser.in1!=null && parser.out1==null && Tools.looksLikeOutputSequenceStream(arg)){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		return parser;
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Expands a primary '#' path into mate paths when no secondary path was supplied.
	 * Input expands only if the literal name does not exist; output always expands.
	 * Also rejects missing primary input or secondary output without primary output.
	 */
	private void doPoundReplacement(){
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}

		//Do output file # replacement.
		//Asymmetry vs input above is intentional: input only expands # when no literal file by that
		//name exists (so a real file named with # still works); output always expands (it's being
		//created, nothing to collide with) - do not add a File.exists() guard here.
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}

		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}

		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
	}

	/** Resolves existing compressed/uncompressed alternatives for the two input paths. */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
	}

	/** Checks primary/mate input access, output permissions and duplicate sequence paths. */
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

	/**
	 * Adjust interleaved mode based on number of input and output files.
	 * Two input files forces non-interleaved mode.
	 * Two outputs with one input force interleaving only if no explicit mode was parsed.
	 * The shared FASTQ interleaving flags are not restored after processing.
	 */
	private void adjustInterleaving(){
		//Adjust interleaved detection based on the number of input files
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}

		//Adjust interleaved settings based on number of output files
		if(!setInterleaved){
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){//If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{//There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}
	}

	/** Reserved initialization hook; currently makes no changes. */
	private static void checkStatics(){
		//Empty
	}

	/*--------------------------------------------------------------*/
	/*----------------       Primary Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates the reader and optional writer, configures sampling, then consumes data.
	 * Limits are forwarded in the selected reader's units. Shared SAM headers are
	 * requested only for SAM/BAM input to SAM/BAM output; Read conversion is requested
	 * for sequence output. A missing output leaves the writer null.
	 * @param t Timer for tracking elapsed time
	 */
	private void process(final Timer t){
		final boolean inputReads=(ffin1!=null && !ffin1.samOrBam());
		final boolean inputSam=(ffin1!=null && ffin1.samOrBam());
		final boolean outputReads=(ffout1!=null && !ffout1.samOrBam());
		final boolean outputSam=(ffout1!=null && ffout1.samOrBam());
		final boolean saveHeader=inputSam && outputSam;

		final Streamer st=StreamerFactory.makeStreamer(ffin1, ffin2, qfin1, qfin2, ordered, maxReads,
			saveHeader, outputReads, threadsIn);
		st.setSampleRate(samplerate, sampleseed);
		final Writer fw=WriterFactory.makeWriter(ffout1, ffout2, qfout1, qfout2, threadsOut, null, saveHeader);

		process(st, fw, t, inputReads || outputReads);
	}

	/**
	 * Starts the reader/writer and consumes batches through null, preserving empty batches.
	 * Skips entries after sampling, counts remaining data, and forwards original batch IDs.
	 * Consumption failures request reader close and writer abort before being rethrown;
	 * startup and post-loop finalization are outside that catch. Normal completion drains
	 * the writer, prints totals, closes the reader, then checks observed error flags.
	 * Reader close has implementation-specific completion semantics, not an added join.
	 * For input thread hints 0 or 1, disables global constructor validation without
	 * restoration; Read mode still validates each retained read in processReadPair.
	 * @param st Unstarted input Streamer, configured before this call
	 * @param fw Unstarted output Writer, or null
	 * @param t Timer for tracking elapsed time
	 * @param readMode True for Read objects, false for SamLine objects
	 */
	private void process(final Streamer st, final Writer fw, final Timer t, final boolean readMode){
		if(threadsIn==0 || threadsIn==1){Read.VALIDATE_IN_CONSTRUCTOR=false;}
		st.start();
		if(fw!=null){fw.start();}
		try{
			if(readMode){
				for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()){
					if(skipreads>0){skipReads(ln.list);}
					for(final Read r : ln){
						processReadPair(r, r.mate);
					}
					if(fw!=null){fw.addReads(ln);}
				}
			}else{
				for(ListNum<SamLine> ln=st.nextLines(); ln!=null; ln=st.nextLines()){
					final ArrayList<SamLine> list=ln.list;
					if(skipreads>0){skipReads(list);}
					for(int i=0, len=list.size(); i<len; i++){
						final SamLine sl=list.get(i);
						final boolean keep=processSamLine(sl);
						//Current routing is SAM/BAM to SAM/BAM or null, so the private helper keeps all
						//records here. Retain its result handling; this is not a subclass override hook.
						if(!keep){list.set(i, null);}
					}
					if(fw!=null){fw.addLines(ln);}
				}
			}
		}catch(final Throwable x){
			//A consumer-side crash (e.g. an assertion on bad pairing) must not strand the pipeline:
			//without this, the non-daemon streamer/writer threads stayed blocked on full queues and the
			//JVM lived forever after main died (jstack-proven, 2026-09-05). Abort both sides, then
			//rethrow so the process still exits loud and nonzero.
			st.close();
			if(fw!=null){fw.finishError();}
			throw new RuntimeException("StreamerWrapper failed mid-stream; output is incomplete.", x);
		}
		boolean errorState=false;
		if(fw!=null){
			errorState=fw.poisonAndWait();
			assert(!readMode || readsIn==fw.readsWritten()) : readsIn+", "+fw.readsWritten()+", "+fw.getClass();
			assert(!readMode || basesIn==fw.basesWritten()) : basesIn+", "+fw.basesWritten()+", "+fw.getClass();
			readsOut=fw.readsWritten();
			basesOut=fw.basesWritten();
		}
		t.stop();
		System.err.println(Tools.timeReadsBasesProcessed(t, readsIn, basesIn, 8));
		st.close();//Prevents a BF4 hang with limited reads
		errorState|=st.errorState();

		if(errorState){throw new RuntimeException("Stream terminated in an error state; the output may be corrupt.");}
	}

	/**
	 * Removes leading list entries up to the remaining skip budget, modifying the list.
	 * Entries are already sampled and may each represent a pair; batch IDs are unchanged.
	 * @param list Mutable batch data
	 * @return Number of entries remaining
	 */
	private int skipReads(final ArrayList<?> list){
		if(skipreads>0){
			if(skipreads>=list.size()){skipreads-=list.size(); list.clear();}else{
				for(int i=0; i<skipreads; i++){list.set(i, null);}
				Tools.condenseStrict(list);
				skipreads=0;
			}
		}
		return list.size();
	}

	/**
	 * Counts a retained read and its mate, then validates any unvalidated objects.
	 * Counters include both reads of a pair and are incremented before validation.
	 * @param r1 Non-null primary read
	 * @param r2 Its mate, or null for unpaired data
	 */
	private void processReadPair(final Read r1, final Read r2){
		readsIn+=r1.pairCount();
		basesIn+=r1.pairLength();
		if(!r1.validated()){r1.validate(true);}
		if(r2!=null && !r2.validated()){r2.validate(true);}
	}

	/**
	 * Counts one SAM record and its sequence length, then evaluates the output condition.
	 * Current raw-line routing always uses SAM/BAM or no output, so this returns true.
	 * @param sl SamLine to process
	 * @return True for nonempty sequence, absent output, or SAM/BAM output
	 */
	private boolean processSamLine(final SamLine sl){
		final int len=sl.lengthOrZero();
		readsIn++;
		basesIn+=len;

		final boolean keep=(len>0 || ffout1==null || ffout1.samOrBam());
		return keep;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;
	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path */
	private String out2=null;

	/** Qual1 input file path */
	private String qfin1=null;
	/** Qual2 input file path */
	private String qfin2=null;
	/** Qual1 output file path */
	private String qfout1=null;
	/** Qual2 output file path */
	private String qfout2=null;

	/** Primary input file format */
	private FileFormat ffin1;
	/** Secondary input file format */
	private FileFormat ffin2;
	/** Primary output file format */
	private FileFormat ffout1;
	/** Secondary output file format */
	private FileFormat ffout2;

	/** Reader-selection thread hint; negative selects the format default. */
	private int threadsIn=-1;
	/** Writer-selection thread hint; negative selects the format default. */
	private int threadsOut=-1;
	/** Limit forwarded to each input reader in that reader's units; negative is unlimited. */
	private long maxReads=-1;
	/** Remaining returned entries to skip after sampling; paired entries count once. */
	private long skipreads=-1;
	/** Sampling rate forwarded to the reader; 1 keeps all otherwise eligible records. */
	private float samplerate=1f;
	/** Sampling seed forwarded to the reader, which defines the sampling algorithm. */
	private long sampleseed=17;
	/** Prevents local SAM parse-flag suppression; does not restore previously disabled flags. */
	private boolean forceParse=false;

	/** Individual reads counted after sampling/skipping, including mates in Read mode. */
	private long readsIn=0;
	/** Bases counted after sampling/skipping, including mate bases in Read mode. */
	private long basesIn=0;
	/** Writer-reported reads after normal writer completion; zero without output. */
	private long readsOut=0;
	/** Writer-reported bases after normal writer completion; zero without output. */
	private long basesOut=0;

	/** Overwrite existing output files */
	private boolean overwrite=false;
	/** Append to existing output files */
	private boolean append=false;
	/** Reader ordering request; the two-input factory path forces ordering. */
	private boolean ordered=true;
	/** Whether parsing explicitly set interleaving, suppressing output-based inference. */
	private boolean setInterleaved=false;
	/** Argument/setup status stream; final timing is printed directly to System.err. */
	private PrintStream outstream=System.err;

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/** Parsed verbosity flag; this class does not test it when printing messages. */
	public static boolean verbose=false;

}
