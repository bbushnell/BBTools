package jgi;

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
import stream.SamLine;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.ListNum;

/**
 * Merges multiple SAM or BAM files into a single output file.
 * Unlike the original {@link MergeSam}, which concatenated files line-by-line
 * via ByteFile and therefore only worked on text SAM, this version reads each
 * input through a Streamer (native BGZF/BAM decoding) and writes through a
 * Writer, so BAM inputs and outputs are handled correctly.
 *
 * The output header is taken from the first input file (all inputs are assumed
 * to share the same reference dictionary).  Records are streamed as raw SamLines
 * and passed through unchanged, so the merge never reconstructs reads from the
 * CIGAR/MD tags and therefore needs no reference.  A single drain thread assigns
 * globally-ascending batch ids to the writer; the heavy work (BGZF inflate on
 * input, deflate on output) runs in the Streamer's and Writer's own thread pools.
 *
 * @author UMP45
 * @date June 7, 2026
 */
public class MergeSam2 {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Program entry point.
	 * @param args Command line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		MergeSam2 x=new MergeSam2(args);
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructs a MergeSam2 instance from command-line arguments.
	 * Collects input files (positional or in=), the output file, and standard flags.
	 * @param args Command line arguments
	 */
	public MergeSam2(String[] args){

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

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("in")){
				//Allow comma-delimited lists as well as repeated in= flags
				for(String s : b.split(",")){in.add(s);}
			}else if(a.equals("out")){
				out=b;
			}else if(b==null && new File(arg).exists()){
				in.add(arg);//Positional input file (matches mergesam.sh usage)
			}else if(parser.parse(arg, a, b)){//Parse standard flags (ow, zl, t, etc.)
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		{//Process parser fields
			overwrite=parser.overwrite;
			append=parser.append;
		}

		if(out!=null && out.equalsIgnoreCase("null")){out=null;}

		if(in.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}

		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
		}
		if(!Tools.testInputFiles(false, true, in.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");
		}

		//SAM as the default format; the actual .sam/.bam extension on each file overrides it.
		ffout=FileFormat.testOutput(out, FileFormat.SAM, null, true, overwrite, append, false);
		ffin=FileFormat.testInputList(in, FileFormat.SAM, null, true, true);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Reads every input file and writes all alignment records to the output.
	 * The output header is loaded once from the first input file.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){

		//Take the header from the first input file, forced to SO:unsorted (this is a
		//concatenation merge, so the output is never sorted even if an input claimed to be).
		final ArrayList<byte[]> header=forceUnsorted(StreamerFactory.loadSharedHeader(ffin[0]));

		//Single shared output writer for all inputs
		final Writer fw=(ffout==null ? null : WriterFactory.makeWriter(ffout, null, -1, header, false));
		if(fw!=null){fw.start();}

		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;

		//Globally-ascending batch id across all input files (required by the ordered Writer)
		long id=0;

		//Drain each input file in turn into the shared writer
		for(FileFormat ff : ffin){
			final Streamer ss=StreamerFactory.makeStreamer(ff, 0, false, maxReads, false, false, -1);
			ss.start();
			if(verbose){outstream.println("Started streamer for "+ff.name());}

			//makeReads=false: stream raw SamLines and pass them through unchanged
			for(ListNum<SamLine> ln=ss.nextLines(); ln!=null && ln.size()>0; ln=ss.nextLines()){
				for(SamLine sl : ln){
					readsProcessed++;
					basesProcessed+=sl.lengthOrZero();
				}
				if(fw!=null){fw.addLines(new ListNum<SamLine>(ln.list, id));}
				id++;
			}

			ReadWrite.closeStream(ss);
		}

		if(fw!=null){
			errorState|=fw.poisonAndWait();
			readsOut=fw.readsWritten();
			basesOut=fw.basesWritten();
		}

		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/**
	 * Forces the @HD line's sort-order (SO:) field to "unsorted" in a copy of the header.
	 * A concatenation merge does not preserve any input's sort order, so the output must
	 * not advertise one.  Rewrites an existing SO: token, appends SO:unsorted if the @HD
	 * line lacks one, and synthesizes a minimal @HD line if the header has none.
	 * @param header Header lines from the first input file (not modified)
	 * @return A new header list with SO:unsorted guaranteed
	 */
	private static ArrayList<byte[]> forceUnsorted(ArrayList<byte[]> header){
		if(header==null){return header;}
		final ArrayList<byte[]> out=new ArrayList<byte[]>(header.size()+1);
		boolean foundHD=false;
		for(byte[] lineBytes : header){
			String line=new String(lineBytes);
			if(line.startsWith("@HD")){
				foundHD=true;
				String[] fields=line.split("\t");
				boolean foundSO=false;
				StringBuilder sb=new StringBuilder();
				for(int j=0; j<fields.length; j++){
					if(j>0){sb.append('\t');}
					if(fields[j].startsWith("SO:")){
						sb.append("SO:unsorted");
						foundSO=true;
					}else{
						sb.append(fields[j]);
					}
				}
				if(!foundSO){sb.append("\tSO:unsorted");}
				out.add(sb.toString().getBytes());
			}else{
				out.add(lineBytes);
			}
		}
		if(!foundHD){out.add(0, "@HD\tVN:1.6\tSO:unsorted".getBytes());}
		return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final ArrayList<String> in=new ArrayList<String>();
	private String out="stdout.sam";

	private long maxReads=-1;

	protected long readsProcessed=0;
	protected long basesProcessed=0;
	protected long readsOut=0;
	protected long basesOut=0;

	private final FileFormat[] ffin;
	private final FileFormat ffout;

	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;

}
