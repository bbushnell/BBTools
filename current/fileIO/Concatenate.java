package fileIO;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ListNum;

/**
 * Concatenates byte lines from input files in list order, appending LF to each output line.
 * This class does not parse or transform sequence records; the selected I/O backends
 * handle compression and stream formats. Output defaults to stdout; {@code out=null}
 * consumes and counts input without creating an output writer.
 * <p>An {@code in=} argument replaces the input list. An existing literal path takes
 * precedence over comma splitting; existing bare paths append to the list.
 * The parsed {@code reads=} limit is currently unused.
 * <p>Instances have unsynchronized, cumulative counters and are intended for one
 * processing caller. Argument parsing also updates shared parser settings.
 *
 * @author Brian Bushnell
 * @date January 21, 2025
 */
public class Concatenate{

	/** Runs concatenation and closes owned diagnostics, leaving System.out and System.err open. */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		Concatenate x=new Concatenate(args);
		
		//Run the object
		try{
			x.process(t);
		}finally{
			//Close the print stream if it was redirected
			Shared.closeStream(x.outstream);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses input paths and output settings without opening the output writer.
	 * After successful preparsing, a constructor failure closes owned diagnostics;
	 * a distinct closing failure is suppressed on the original cause.
	 * On successful construction, main owns diagnostic-stream cleanup.
	 * @param args Command-line arguments, including shared parser options
	 */
	public Concatenate(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		try{
			Parser parser=new Parser();
			parser.out1=out1;
			for(int i=0; i<args.length; i++){
				String arg=args[i];
				String[] split=arg.split("=");
				String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1] : null;
				if(b!=null && b.equalsIgnoreCase("null")){b=null;}

				if(a.equals("parse_flag_goes_here")){
					//Unused template option; no setting is changed.
				}else if(a.equals("verbose")){
					verbose=Parse.parseBoolean(b);
				}else if(a.equals("in")){
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

			{//Process parser fields
				Parser.processQuality();

				maxReads=parser.maxReads;
				out1=parser.out1;
			}

			ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, parser.overwrite, parser.append, false);
		}catch(RuntimeException | Error cause){
			try{
				Shared.closeStream(outstream);
			}catch(RuntimeException | Error closing){
				if(cause!=closing){cause.addSuppressed(closing);}
			}
			throw cause;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Opens the optional writer, drains each input in order, then finishes output.
	 * Rejects matching raw names and existing ordinary-file aliases before opening output.
	 * A processing failure remains primary if writer finishing also fails; reported input
	 * or output errors prevent successful return. Diagnostic cleanup belongs to main.
	 * @param t Timer stopped on successful completion
	 */
	void process(Timer t){

		final ByteStreamWriter bsw;
		if(out1!=null){
			final int buff=4;//Unused template setting.

			final File output=new File(ffout1.name());
			for(String s : in){
				assert(!out1.equalsIgnoreCase(s)) : "Input file and output file have same name.";
				if(ffout1.file() && output.exists()){
					final File input=new File(s.trim().replace('\\', '/'));
					if(input.exists()){
						try{
							if(Files.isSameFile(input.toPath(), output.toPath())){
								throw new IllegalArgumentException("Input file and output file refer to the same file: "+s+" and "+out1);
							}
						}catch(IOException e){
							throw new RuntimeException("Could not compare input and output before opening output: "+s+" and "+out1, e);
						}
					}
				}
			}
			
			bsw=ByteStreamWriter.makeBSW(ffout1);
		}else{bsw=null;}
		
		boolean inputError=false;
		Throwable failure=null;
		try{
			for(String s : in){
				inputError|=processInner(s, bsw);
			}
		}catch(RuntimeException | Error cause){
			failure=cause;
			throw cause;
		}finally{
			try{
				if(bsw!=null && bsw.poisonAndWait()){
					throw new RuntimeException("An error occurred while finishing output; the output may be incomplete: "+out1);
				}
			}catch(RuntimeException | Error closing){
				// Keep the caller's failure primary when output finishing also fails.
				if(failure==null){throw closing;}
				if(failure!=closing){failure.addSuppressed(closing);}
			}
		}
		if(inputError){
			throw new RuntimeException("An error occurred while reading input; the output may be incomplete.");
		}
		if(verbose){outstream.println("Finished.");}
		
		t.stop();
		if(verbose){
			outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 12));
			outstream.println(bytesProcessed);
		}
	}
	
	/**
	 * Drains one newly created reader and closes it after processing, including failure.
	 * Counts each returned line and its length plus one logical LF byte, even without output.
	 * Processing failures propagate with distinct reader-closing failures suppressed.
	 * @param fname Input path or stream name interpreted by FileFormat
	 * @param bsw Borrowed writer, or null to discard lines; this method does not finish it
	 * @return Reader-reported error state after otherwise normal processing and close
	 */
	boolean processInner(String fname, ByteStreamWriter bsw){
		FileFormat ffin=FileFormat.testInput(fname, FileFormat.TXT, null, true, true);
		
		final ByteFile bf=ByteFile.makeByteFile(ffin);
		
		boolean inputError=false;
		Throwable failure=null;
		try{
			ListNum<byte[]> ln;
			while((ln=bf.nextList())!=null){
				for(byte[] line : ln.list){
					linesProcessed++;
					bytesProcessed+=(line.length+1);
					if(bsw!=null){bsw.println(line);}
				}
			}
		}catch(RuntimeException | Error cause){
			failure=cause;
			throw cause;
		}finally{
			try{
				inputError=bf.close();
				if(inputError && failure!=null){
					throw new RuntimeException("An error occurred while closing input after processing failed: "+fname);
				}
			}catch(RuntimeException | Error closing){
				if(failure==null){throw closing;}
				if(failure!=closing){failure.addSuppressed(closing);}
			}
		}
		return inputError;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Input names in processing order; each in= argument replaces earlier entries. */
	private ArrayList<String> in=new ArrayList<String>();
	/** Output name; null disables writing while stdout.txt selects standard output. */
	private String out1="stdout.txt";
	
	/** Output descriptor retaining parsed overwrite/append settings; null when output is disabled. */
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/

	/** Parsed reads= limit, currently unused by the line-processing loop. */
	private long maxReads=-1;
	/** Cumulative returned-line count and logical bytes (line length plus LF), not physical file size. */
	private long linesProcessed=0, bytesProcessed=0;
	
	/*--------------------------------------------------------------*/
	
	/** Diagnostic stream; redirected streams are owned, standard streams are borrowed. */
	private java.io.PrintStream outstream=System.err;
	/** Shared switch enabling completion and throughput diagnostics. */
	public static boolean verbose=false;
	
}
