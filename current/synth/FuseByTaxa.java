package synth;

import java.io.File;
import java.io.PrintStream;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Fuses adjacent taxon-ID groups into N-padded FASTA chunks, retaining chunks
 * that meet a minimum length. Originally used to turn combined RefSeq clade
 * files into assembly-attempt sequences for contig-classification benchmarks,
 * dropping taxa represented only by short markers.
 * <p>
 * Input must keep each taxon's contigs adjacent. This class compares the literal
 * tokens extracted by {@link #parseTid(String)}; it neither sorts input nor
 * merges groups separated by another token. Use unpaired contigs: processing
 * visits primary entries returned by the reader, not their linked mates.
 * <p>
 * A nonempty chunk is flushed before adding a contig that would exceed maxlen,
 * including its padding. Individual contigs are not split, so one contig can
 * exceed maxlen; see SYN016. Each chunk is filtered independently by floor.
 * A null output still computes eligible-chunk and adjacent-group counters.
 * Instances are intended for one processing pass on one caller thread; argument
 * parsing also configures shared reader/compression settings.
 * <p>
 * Example: {@code fusebytaxa.sh in=sorted.fna.gz out=genomes.fna.gz floor=1m pad=10 maxlen=50m}
 *
 * @author Noire
 * @date July 16, 2026
 */
public class FuseByTaxa {

	/** Runs one fusion pass and closes redirected status output after normal completion.
	 * @param args Command-line options, including input/output and fusion limits
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		FuseByTaxa x=new FuseByTaxa(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	/** Parses options and constructs format descriptors without opening the processing streams.
	 * Format detection may inspect the input. Assertions require an input name and
	 * maxlen at least floor; other numerical ranges are not comprehensively checked.
	 * @param args Command-line options; PreParser also expands configuration files
	 */
	public FuseByTaxa(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("floor") || a.equals("minlen")){
				floor=Parse.parseKMG(b);
			}else if(a.equals("pad") || a.equals("npad")){
				pad=Integer.parseInt(b);
			}else if(a.equals("maxlen")){
				maxlen=Parse.parseKMG(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//standard flags
			}else if(i==0 && Tools.looksLikeInputSequenceStream(arg)){
				parser.in1=arg;
			}else if(i==1 && parser.in1!=null){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		Parser.processQuality();
		overwrite=parser.overwrite;
		in1=parser.in1;
		out1=parser.out1;
		assert(in1!=null) : "No input file specified.";
		assert(maxlen>=floor) : "maxlen ("+maxlen+") must be >= floor ("+floor+").";

		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, null, true, true);
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, null, true, overwrite, false, false);
	}

	/** Consumes adjacent groups, filters generated chunks, and prints final counters.
	 * Normal terminal consumption is followed by reader close and writer completion;
	 * reported errors cause failure before success statistics. Earlier thrown failures
	 * are not protected by a general finally-based cleanup path.
	 * @param t Caller-started timer, stopped after normal stream cleanup
	 * @throws RuntimeException If stream completion reports an input or output error
	 */
	void process(Timer t){
		checkOutputWritable();
		//Ordered batches and the single-reader hint preserve adjacency; the factory selects the concrete reader.
		Streamer st=StreamerFactory.makeStreamer(ffin1, null, true, -1L, false, true, 1);
		st.start();
		final ByteStreamWriter bsw=(ffout1==null ? null : new ByteStreamWriter(ffout1));
		if(bsw!=null){bsw.start();}

		final ByteBuilder bb=new ByteBuilder();   //current fused chunk (bases only)
		String curTid=null;
		int chunk=0;
		boolean emittedThisTaxon=false;

		for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()){
			for(Read r : ln){
				readsIn++;
				basesIn+=r.length();
				final String tid=parseTid(r.id);
				if(curTid==null){curTid=tid;}
				if(!tid.equals(curTid)){
					//New taxon: flush the current chunk, then account for whether the taxon was kept.
					if(emit(bsw, curTid, chunk, bb)){emittedThisTaxon=true;}
					if(emittedThisTaxon){taxaKept++;}else{taxaDropped++;}
					bb.clear(); curTid=tid; chunk=0; emittedThisTaxon=false;
				}else if(bb.length()>0 && (long)bb.length()+pad+r.length()>maxlen){
					//Same taxon but the chunk is full: flush and continue a new chunk.
					if(emit(bsw, curTid, chunk, bb)){emittedThisTaxon=true;}
					bb.clear(); chunk++;
				}
				if(bb.length()>0){for(int k=0; k<pad; k++){bb.append((byte)'N');}}
				//TODO: Probable bug - SYN016: a contig longer than maxlen is appended whole,
				//despite the documented chunk bound. Oversized-contig policy needs clarification.
				bb.append(r.bases);
			}
		}
		//Flush the final taxon.
		if(curTid!=null){
			if(emit(bsw, curTid, chunk, bb)){emittedThisTaxon=true;}
			if(emittedThisTaxon){taxaKept++;}else{taxaDropped++;}
		}

		st.close();
		boolean errorState=st.errorState();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		t.stop();
		if(errorState){
			throw new RuntimeException("FuseByTaxa encountered an input or output error; output may be incomplete.");
		}
		outstream.println(Tools.timeReadsBasesProcessed(t, readsIn, basesIn, 8));
		outstream.println("Taxa kept:    \t"+taxaKept);
		outstream.println("Taxa dropped: \t"+taxaDropped+" (below floor="+floor+")");
		outstream.println("Sequences out:\t"+seqsOut+"  ("+basesOut+" bp)");
	}

	/** Fails unwritable output before starting an input streamer that could block on a full queue. */
	private void checkOutputWritable(){
		if(ffout1==null || ffout1.canWrite()){return;}
		final File f=new File(ffout1.name());
		throw new RuntimeException("File "+ffout1.name()+" exists "+(f.canWrite() ?
				("and overwrite="+overwrite+".\nPlease add the flag ow to overwrite the file.") :
				"and is read-only."));
	}

	/** Counts a nonempty chunk meeting floor and queues its FASTA record when a writer exists.
	 * Counters and the return value describe eligibility, not confirmed output writes.
	 * @param bsw Started output writer, or null to count without writing
	 * @param tid Literal group token placed in the generated header
	 * @param chunk Zero-based chunk index within this adjacent group
	 * @param bb Current bases; not cleared or modified by this method
	 * @return true if the chunk is nonempty and meets floor, including with a null writer
	 */
	private boolean emit(ByteStreamWriter bsw, String tid, int chunk, ByteBuilder bb){
		if(bb.length()<floor || bb.length()==0){return false;}
		seqsOut++;
		basesOut+=bb.length();
		if(bsw!=null){
			header.clear();
			header.append('>').append("tid|").append(tid).append("|fused").append(chunk)
				.append(" len=").append(bb.length()).nl();
			bsw.print(header);
			bsw.print(bb);
			bsw.print(NEWLINE);
		}
		return true;
	}

	/** Extracts the token after the first pipe, ending at the next pipe or end of the ID.
	 * The prefix and token are not validated; the token may be empty or nonnumeric.
	 * @param id Read identifier, possibly null
	 * @return Extracted token, or "NA" when id is null or contains no pipe
	 */
	private static String parseTid(String id){
		if(id==null){return "NA";}
		int a=id.indexOf('|');
		if(a<0){return "NA";}
		int b=id.indexOf('|', a+1);
		return b<0 ? id.substring(a+1) : id.substring(a+1, b);
	}

	/*--------------------------------------------------------------*/

	/** Input name; required by the constructor assertion. */
	private String in1=null;
	/** Output name; null keeps counting but suppresses FASTA writing. */
	private String out1=null;
	/** Detected input format used by StreamerFactory. */
	private FileFormat ffin1;
	/** Output descriptor, or null when output is disabled. */
	private FileFormat ffout1;

	/** Minimum eligible nonempty chunk length in bases, inclusive. */
	private long floor=1000000;
	/** Requested Ns between contigs retained in the same chunk. */
	private int pad=10;
	/** Pre-append chunk limit in bases; does not split an oversized individual contig. */
	private long maxlen=50000000;

	/** Primary reader entries and their bases consumed by process. */
	private long readsIn=0, basesIn=0;
	/** Eligible generated chunks and bases, including when output is disabled. */
	private long seqsOut=0, basesOut=0;
	/** Adjacent groups with at least one eligible chunk, and groups with none. */
	private long taxaKept=0, taxaDropped=0;

	/** Reused header scratch for emitted FASTA records. */
	private final ByteBuilder header=new ByteBuilder(64);
	/** FASTA sequence-line terminator. */
	private static final byte[] NEWLINE={(byte)'\n'};

	/** Whether an existing output may be overwritten. */
	private boolean overwrite=false;
	/** Status destination selected by PreParser. */
	private PrintStream outstream=System.err;
	/** Parsed verbosity option; this class does not currently gate diagnostics on it. */
	public static boolean verbose=false;
}
