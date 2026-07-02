package driver;

import java.io.PrintStream;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * Renames NCBI sequence headers to TID format in FASTA files.
 * Processes text files line by line, converting ">ncbi" headers to ">tid" format
 * and inserting a pipe character between the identifier and description.
 *
 * @author Brian Bushnell
 * @date Oct 17, 2014
 */
public class RenameNcbiToTid {
	
	public static void main(String[] args){
		Timer t=new Timer();
		RenameNcbiToTid x=new RenameNcbiToTid(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public RenameNcbiToTid(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ReadWrite.verbose=verbose;
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;

			out1=parser.out1;
		}
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TEXT, null, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.TEXT, null, true, true);
	}
	
	void process(Timer t){
		
		final TextFile tf;
		{
			tf=new TextFile(ffin1);
			if(verbose){outstream.println("Started tf");}
		}
		
		final TextStreamWriter tsw;
		{
			//NOTE [driver/RenameNcbiToTid#002] LOW/dev: out= is optional-looking (only in1 required, ctor L81) but NOT
			//optional — omitting out= gives ffout1=testOutput(null)=null → this new TextStreamWriter(null) NPEs at ctor
			//(TextStreamWriter L41 ff.fastq()). Dead `tsw!=null` guard L119 = optional-output half-wired. Pattern (f),
			//same as MergeBigelow#001/FilterAssemblySummary#001. Dev one-off; a real user passes out=. errorState below
			//IS correctly propagated (L123-124 capture + L129 throw) — studied praise.
			tsw=new TextStreamWriter(ffout1);
			tsw.start();
			if(verbose){outstream.println("Started tsw");}
		}
		
		long linesProcessed=0;
		long charsProcessed=0;
		
		{
			String line;
			while((maxReads<0 || linesProcessed<maxReads) && (line=tf.nextLine())!=null){
				linesProcessed++;
				charsProcessed+=line.length();
				String result=processLine(line);
				if(tsw!=null && result!=null){tsw.println(result);}
			}
		}
		
		errorState|=tsw.poisonAndWait();
		errorState|=tf.close();
		
		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, charsProcessed, 8));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	private static String processLine(String line){
		if(line.startsWith(">ncbi")){
			line=line.replaceFirst(">ncbi", ">tid"); //n ">ncbi" has no regex metachars → replaceFirst is safe/literal here.
			int firstSpace=line.indexOf(' ');
			//TODO: Possible bug [driver/RenameNcbiToTid#001] LOW/dev: a ">ncbi..." header with NO space → firstSpace=-1 →
			//substring(0,-1) throws StringIndexOutOfBoundsException (crash, not a clean skip). Real for an accession-only
			//header like ">ncbi12345" (no description). Unguarded structured-field index (pattern e). Recommended trivial
			//guard: `if(firstSpace>=0){ line=line.substring(0,firstSpace)+"|"+line.substring(firstSpace+1); }` — i.e. a
			//spaceless header becomes ">tid12345" with no pipe (no description to separate). Left un-applied: dev one-off
			//(no .sh, no callers) and Brian's real NCBI headers always have a space; the else-behavior is a mild judgment call.
			line=line.substring(0, firstSpace)+"|"+line.substring(firstSpace+1);
		}
		return line;
	}
	
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
