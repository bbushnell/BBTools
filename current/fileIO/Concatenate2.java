package fileIO;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.util.ArrayList;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;

/**
 * Copies bytes from input streams in list order to one optional output stream.
 * ReadWrite selects compression and stdio backends; the copy loop preserves byte
 * order and line endings. With out=null, inputs are consumed and counted without
 * writing. Counts represent processed bytes and LF symbols.
 * Instances have single-caller state; construction also changes shared parser settings.
 * @author Brian Bushnell, Shinobu
 * @date January 21, 2025
 */
public class Concatenate2{

	/**
	 * Runs concatenation and closes owned diagnostics after processing.
	 * Borrowed standard diagnostic streams remain open.
	 * @param args Input/output arguments and shared parser options
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();

		//Create an instance of this class
		Concatenate2 x=new Concatenate2(args);

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
	 * Parses input names and output options without opening concatenation streams.
	 * Each in= replaces the input list: an existing literal path is kept whole,
	 * otherwise comma-separated names are used. Bare existing paths append to the list.
	 * After successful preparsing, constructor rejection closes owned diagnostics
	 * and preserves the primary failure, suppressing a distinct closing failure.
	 * On successful construction, main owns diagnostic cleanup.
	 * @param args Command-line arguments, including shared parser options
	 */
	public Concatenate2(String[] args){

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
				overwrite=parser.overwrite;
				append=parser.append;
			}
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
	 * Checks input/output names, opens optional output, and copies inputs in order.
	 * Returned input errors accumulate while later inputs are processed; output is
	 * finished before a reported completion error is thrown. Counters accumulate
	 * in this instance, including when output is disabled.
	 * @param t Timer stopped on successful completion
	 * @throws RuntimeException If an existing-file check or reported completion fails
	 * @throws AssertionError If raw input/output names match ignoring case
	 */
	void process(Timer t){

		final FileFormat ffout=FileFormat.testOutput(out1, FileFormat.TEXT, null, true, overwrite, append, true);
		if(out1!=null){
			// Check the raw backend path before output can truncate any existing input.
			final String name=out1.replace('\\', '/');
			final boolean stdio=name.equals("stdout") || name.startsWith("stdout.") ||
					name.equals("stderr") || name.startsWith("stderr.");
			final File output=new File(name.replace('|', '_'));
			for(String s : in){
				assert(!out1.equalsIgnoreCase(s)) : "Input file and output file have same name.";
				if(!stdio && !name.equals("/dev/null") && output.exists()){
					// Raw input normalizes backslashes but preserves pipes and whitespace.
					final File input=new File(s.replace('\\', '/'));
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
		}
		if(out1!=null && !overwrite && !append){
			// Match ReadWrite's raw output name and stream aliases; FileFormat trims names.
			final String name=out1.replace('\\', '/');
			final boolean stdio=name.equals("stdout") || name.startsWith("stdout.") ||
					name.equals("stderr") || name.startsWith("stderr.");
			if(!stdio && !name.equals("/dev/null") && new File(name.replace('|', '_')).exists()){
				throw new RuntimeException("Output file exists and overwrite=false: "+out1);
			}
		}
		final OutputStream os=(ffout!=null ? ReadWrite.getOutputStream(out1, ffout.append(), false, true) : null);

		boolean error=false;
		for(String s : in){
			error|=processInner(s, os);
		}

		if(out1!=null){error|=ReadWrite.finishWriting(null, os, out1, ffout.subprocess);}
		if(error){throw new RuntimeException("Concatenate2 encountered an input or output error; output may be corrupt.");}
		if(verbose){outstream.println("Finished.");}

		t.stop();
		if(verbose){
			outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 12));
			outstream.println(bytesProcessed);
		}
	}

	/**
	 * Copies one input through a 256 KiB buffer, counting bytes and LF symbols
	 * after each optional write. Borrows output and finishes only the input here.
	 * Exceptions in the copy loop invoke KillSwitch.exceptionKill; setup and
	 * unchecked completion failures propagate.
	 * @param fname Input path or stream name supplied to FileFormat and ReadWrite
	 * @param os Borrowed output, or null to consume and count without writing
	 * @return true if input closure or its subprocess reports an error
	 */
	boolean processInner(String fname, OutputStream os){
		final FileFormat ffin=FileFormat.testInput(fname, FileFormat.TXT, null, true, true);
		final byte[] buffer=new byte[1024*256];
		final InputStream is=ReadWrite.getInputStream(fname, false, true);

		try{
			for(int r=is.read(buffer); r>0; r=is.read(buffer)){
				if(os!=null){os.write(buffer, 0, r);}
				bytesProcessed+=r;
				linesProcessed+=Vector.countSymbols(buffer, 0, r, (byte)'\n');
			}
		}catch(Exception e){
			KillSwitch.exceptionKill(e);
		}
		return ReadWrite.finishReading(is, fname, ffin.subprocess);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Input paths or stream names in traversal order. */
	private ArrayList<String> in=new ArrayList<String>();
	/** Output name; stdout.txt selects stdout, and null disables writing. */
	private String out1="stdout.txt";
	/** Parsed permission to replace an existing output file. */
	private final boolean overwrite;
	/** Parsed request to append to output; takes precedence over overwrite. */
	private final boolean append;

	/*--------------------------------------------------------------*/

	/** Parsed common reads= setting; the byte-copy loop does not apply it. */
	private long maxReads=-1;
	/** Accumulated LF-symbol count and byte count across processed input chunks. */
	private long linesProcessed=0, bytesProcessed=0;

	/*--------------------------------------------------------------*/

	/** Owned redirected diagnostics or a borrowed standard stream. */
	private java.io.PrintStream outstream=System.err;
	/** Controls successful-completion diagnostics, including counter reporting. */
	public static boolean verbose=false;

}
