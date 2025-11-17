package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;

/**
 * @author Brian Bushnell
 * @date March 8, 2017
 *
 */
public class MergeSam {
	
	/** Program entry point for merging SAM files.
	 * @param args Command-line arguments including input/output paths and options */
	public static void main(String[] args){
		Timer t=new Timer();
		MergeSam x=new MergeSam(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor that parses command-line arguments and initializes file formats.
	 * Sets up input/output streams, validates file paths, and configures processing options.
	 * @param args Command-line arguments array
	 * @throws RuntimeException if no input files specified or output files cannot be written
	 */
	public MergeSam(String[] args){
		
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

			if(a.equals("invalid")){
				outInvalid=b;
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("in")){
				in.add(b);
			}else if(a.equals("out")){
				out=b;
			}else if(b==null && new File(arg).exists()){
				in.add(arg);
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			overwrite=parser.overwrite;
			append=parser.append;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in==null){throw new RuntimeException("Error - at least one input file is required.");}

		if(out!=null && out.equalsIgnoreCase("null")){out=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out+"\n");
		}

		ffout=FileFormat.testOutput(out, FileFormat.SAM, null, true, overwrite, append, false);
		ffoutInvalid=FileFormat.testOutput(outInvalid, FileFormat.SAM, null, true, overwrite, append, false);
		ffin=FileFormat.testInputList(in, FileFormat.SAM, null, true, true);
	}
	
	/**
	 * Main processing method that merges SAM files and filters invalid entries.
	 * Processes files sequentially, maintaining header order by only allowing
	 * header lines (@-prefixed) while in header mode, then switching to alignment mode.
	 *
	 * @param t Timer for tracking execution time
	 * @throws RuntimeException if processing encounters errors
	 */
	void process(Timer t){
		
		ByteStreamWriter bsw=null;
		if(ffout!=null){
			bsw=new ByteStreamWriter(ffout);
			bsw.start();
		}
		
		ByteStreamWriter bswInvalid=null;
		if(ffoutInvalid!=null){
			bswInvalid=new ByteStreamWriter(ffoutInvalid);
			bswInvalid.start();
		}

		boolean headerMode=true;
		for(int fnum=0; fnum<ffin.length; fnum++){

			ByteFile bf=ByteFile.makeByteFile(ffin[fnum]);
			
			byte[] line=bf.nextLine();
			
			while(line!=null){
				if(line.length>0){
					if(maxLines>0 && linesProcessed>=maxLines){break;}
					linesProcessed++;
					bytesProcessed+=line.length;
					
					boolean valid=true;
					if(line[0]=='@'){valid=headerMode;}
					else{headerMode=false;}

					if(valid){
						linesValid++;
						if(bsw!=null){
							bsw.println(line);
						}
					}else{
						if(bswInvalid!=null){
							bswInvalid.println(line);
						}
					}
				}
				line=bf.nextLine();
			}

			errorState|=bf.close();
		}
		
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		if(bswInvalid!=null){errorState|=bswInvalid.poisonAndWait();}
		
		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
		outstream.println();
		outstream.println("Valid Lines:       \t"+linesValid);
		outstream.println("Invalid Lines:     \t"+(linesProcessed-linesValid));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/** List of input SAM file paths to be merged */
	private ArrayList<String> in=new ArrayList<String>();
	/** Output file path for merged SAM data */
	private String out="stdout.sam";
	/** Output file path for invalid SAM entries (optional) */
	private String outInvalid=null;
	
	/*--------------------------------------------------------------*/
	
	/** Total number of lines read from all input files */
	private long linesProcessed=0;
	/** Number of lines determined to be valid SAM entries */
	private long linesValid=0;
	/** Total bytes read from all input files */
	private long bytesProcessed=0;
	
	/** Maximum number of lines to process before stopping */
	private long maxLines=Long.MAX_VALUE;
	
	/*--------------------------------------------------------------*/
	
	/** Array of FileFormat objects for input SAM files */
	private final FileFormat[] ffin;
	/** FileFormat object for the main output SAM file */
	private final FileFormat ffout;
	/** FileFormat object for invalid entries output file */
	private final FileFormat ffoutInvalid;
	
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages and logging */
	private PrintStream outstream=System.err;
	/** Global flag controlling verbose output across multiple classes */
	public static boolean verbose=false;
	/** Flag indicating whether processing encountered errors */
	public boolean errorState=false;
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files instead of overwriting */
	private boolean append=false;
	
}
