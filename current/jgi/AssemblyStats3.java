package jgi;

import java.io.File;
import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;

/**
 * Sequential driver for a nine-column, tab-delimited assembly statistics table.
 * Each input is loaded by Assembly up to its per-input nonempty-record limit before
 * its row is written; the default is unlimited. Lengths
 * count sequence bytes under Assembly's permissive FASTA rules; the four cutoff
 * columns sum lengths of contigs meeting inclusive 5k, 10k, 25k and 50k cutoffs.
 * GC is a fraction of A/C/G/T/U counts, excluding N and other ambiguity counts.
 * The header is emitted once per processing call, including when appending to an
 * existing output. Input names are written as supplied, without tab/newline escaping.
 * Construction parses options and prepares output metadata; processing acquires
 * and finishes the writer. Failure can leave a partial output table. A null output
 * description suppresses writing but still processes every input. Instances are mutable
 * and are not intended for concurrent processing.
 *
 * @author Brian Bushnell
 * @date January 21, 2025
 */
public class AssemblyStats3{

	/** Constructs and processes the driver, then releases owned diagnostic output.
	 * Constructor failures occur before the processing finally block is entered.
	 * @param args Input names and Parser-compatible flag=value arguments */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		AssemblyStats3 x=new AssemblyStats3(args);
		
		try{
			x.process(t);
		}finally{
			//Close owned diagnostics even when processing fails; shared streams stay open.
			Shared.closeStream(x.outstream);
		}
	}
	
	/** Parses inputs and output options without opening the primary output writer.
	 * An in/input/in1/input1 option replaces the accumulated input list; existing
	 * whole filenames take precedence over comma splitting. Bare existing filenames
	 * append to the list. Generic Parser options are accepted, but not all are used
	 * by this driver. The final extout value must be null or empty because the table
	 * format is fixed. reads/maxreads limits
	 * nonempty sequence records separately for each input, before sorting by length.
	 * @param args Arguments processed by PreParser before local and shared parsing
	 * @throws IllegalArgumentException If the final extout value is nonempty */
	public AssemblyStats3(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Parser parser=new Parser();
		parser.out1=out1;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			//JGI-059: preserve '=' characters within input and output values.
			String[] split=Parse.splitOnFirst(arg, '=');
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				in.clear();
				String[] b2=(b==null) ? null : (new File(b).exists() ? new String[] {b} : b.split(","));
				for(String b3 : b2){in.add(b3);}
			}else if(b==null && new File(arg).exists()){
				in.add(arg);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
		}
		
		if(parser.extout!=null && !parser.extout.isEmpty()){
			throw new IllegalArgumentException("extout is not supported by stats3: "+parser.extout+"; output is fixed nine-column TSV.");
		}

		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			out1=parser.out1;
		}
		
		//JGI-064: honor explicit overwrite/append settings with the existing Parser defaults.
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, parser.overwrite, parser.append, false);
	}
	
	/** Loads each input and writes its statistics, finishing any acquired writer on failure.
	 * Normal completion stops the supplied timer without printing a timing report.
	 * A writer acquisition failure occurs before this method's cleanup block.
	 * @param t Nonnull timer started by the caller; stopped only after successful cleanup */
	void process(Timer t){

		final ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout1);
		Throwable failure=null;
		try{
			if(bsw!=null){bsw.println(makeHeader());}
			for(String s : in){processInner(s, bsw);}
		}catch(RuntimeException|Error e){
			failure=e;
			throw e;
		}finally{
			finishWriter(bsw, failure);
		}
		if(verbose){outstream.println("Finished.");}
		
		t.stop();
	}

	/** Finishes an acquired writer and rejects a reported write error.
	 * Preserves a primary RuntimeException/Error, suppressing only a distinct cleanup
	 * failure on it. A null writer requires no cleanup.
	 * @param bsw Nullable writer acquired by process
	 * @param primary Original processing failure, or null for normal completion */
	private void finishWriter(final ByteStreamWriter bsw, final Throwable primary){
		if(bsw==null){return;}
		try{
			final boolean error=bsw.poisonAndWait();
			assert(!bsw.isAlive()) : "ByteStreamWriter.poisonAndWait must join before output is released: "+out1;
			if(error){
				throw new RuntimeException("Error writing assembly statistics to "+out1);
			}
		}catch(RuntimeException|Error e){
			if(primary==null){throw e;}
			if(e!=primary){primary.addSuppressed(e);}
		}
	}
	
	/** Builds the fixed column names without a trailing newline.
	 * @return Header in the same order as processInner's row values */
	public static String makeHeader(){
		ByteBuilder bb=new ByteBuilder();
		bb.append("fname");
		bb.tab().append("size");
		bb.tab().append("contigs");
		bb.tab().append("gc");
		bb.tab().append("maxContig");
		bb.tab().append("5kplus");
		bb.tab().append("10kplus");
		bb.tab().append("25kplus");
		bb.tab().append("50kplus");
		return bb.toString();
	}
	
	/** Loads the requested prefix of one assembly, then writes one row when output is enabled.
	 * Loading still occurs for a null writer. Empty assemblies have maximum contig
	 * length zero; GC uses Assembly's float fraction and the writer's numeric formatting.
	 * @param fname Input name passed to Assembly and written unchanged in the first column
	 * @param bsw Nullable destination; owned and finished by process */
	void processInner(String fname, ByteStreamWriter bsw){
		Assembly a=new Assembly(fname, maxReads);
		if(bsw==null){return;}
		int max=a.contigs.size>0 ? a.contigs.get(0) : 0;
		bsw.print(fname).tab().print(a.length).tab().print(a.contigs.size);
		bsw.tab().print(a.gc(), 3).tab().print(max);
		bsw.tab().print(a.lengthAtLeast(5000));
		bsw.tab().print(a.lengthAtLeast(10000));
		bsw.tab().print(a.lengthAtLeast(25000));
		bsw.tab().print(a.lengthAtLeast(50000));
		bsw.nl();
	}
	
	/*--------------------------------------------------------------*/
	
	/** Ordered input names; each explicit input option replaces earlier entries. */
	private ArrayList<String> in=new ArrayList<String>();
	/** Primary output name, defaulting to standard output through FileFormat. */
	private String out1="stdout.txt";
	
	/** Prepared output metadata; null disables writing, not input loading. */
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/

	/** Per-input nonempty-record limit; negative is unlimited, zero produces empty statistics. */
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	/** Diagnostic stream selected by PreParser; main releases owned streams after processing. */
	private java.io.PrintStream outstream=System.err;
	/** Enables the successful-completion message; not a throughput-report switch. */
	public static boolean verbose=false;
	
}
