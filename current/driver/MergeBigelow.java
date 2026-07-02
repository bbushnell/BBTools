package driver;

import java.io.PrintStream;
import java.util.HashMap;

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
 * Merges two tab-delimited text files based on matching first columns.
 * Processes text by removing specific patterns like SCGC identifiers,
 * converts to lowercase, and replaces commas with underscores.
 *
 * @author Brian Bushnell
 * @date Oct 17, 2014
 */
public class MergeBigelow {
	
	public static void main(String[] args){
		Timer t=new Timer();
		MergeBigelow x=new MergeBigelow(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public MergeBigelow(String[] args){
		
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
			in2=parser.in2;

			out1=parser.out1;
		}
		
		if(in1==null || in2==null){throw new RuntimeException("Error - two input files are required.");}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}

		assert(Tools.testInputFiles(false, true, in1, in2));
		assert(Tools.testForDuplicateFiles(true, in1, in2, out1));
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TEXT, null, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.TEXT, null, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.TEXT, null, true, true);
	}
	
	void process(Timer t){
		
		table=hash(ffin2);
		
		final TextFile tf;
		{
			tf=new TextFile(ffin1);
			if(verbose){outstream.println("Started tf");}
		}
		
		final TextStreamWriter tsw;
		{
			//NOTE [driver/MergeBigelow#001] LOW/dev: out= is optional-LOOKING (out1 defaults null; only in1&in2 are required
			//at ctor L83) but is NOT actually optional here — if out= is omitted, ffout1=FileFormat.testOutput(null,...)=null,
			//and this unconditional new TextStreamWriter(null) NPEs at its ctor (TextStreamWriter L41: ff.fastq() on null ff).
			//The dead `if(tsw!=null...)` guard at L128 shows optional-output WAS intended but was only half-wired. Same shape
			//as FilterAssemblySummary#001 / pattern (f). Dev one-off (no .sh, no callers) so LOW; a real user always passes out=.
			tsw=new TextStreamWriter(ffout1);
			tsw.start();
			if(verbose){outstream.println("Started tsw");}
		}
		
		long linesProcessed=0;
		long charsProcessed=0;
		
		{
			String line;
			while((line=tf.nextLine())!=null){
//				System.err.println("Processing "+line);
				linesProcessed++;
				charsProcessed+=line.length();
				CharSequence result=processLine(line);
				if(tsw!=null && result!=null){tsw.println(result);}
				if(maxReads>0 && linesProcessed>=maxReads){break;}
			}
		}
		
		//n GOOD (studied praise): errorState is fully propagated on the MAIN path — both the writer drain AND the in1
		//n reader close are captured, and L139 throws if set. This is the correct plumbing that SummarizeCrossblock#001
		//n and SummarizeSealCrosstalk#001 both LACK. (The gap is only #002: the in2 reader in hash() is never closed.)
		errorState|=tsw.poisonAndWait();
		errorState|=tf.close();

		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, charsProcessed, 8));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	private CharSequence processLine(String line){
		String[] split=line.split(delimiter);
		String[] split2=table.get(split[0]);
		if(split2==null){return line;} //Header
		StringBuilder sb=new StringBuilder();
		String tab="";
//		assert(false) : split.length+", "+split2.length;
//		System.err.println(split[1]);
		if(split.length>1){
			if(split[1].contains(" SCGC")){
				split[1]=split[1].substring(0, split[1].indexOf(" SCGC"));
//				System.err.println(split[1]);
			}
			if(split[1].contains(" "+split[0])){
				split[1]=split[1].substring(0, split[1].indexOf(" "+split[0]));
//				System.err.println(split[1]);
			}
			split[1]=split[1].toLowerCase();
//			System.err.println(split[1]);
		}
		for(int i=0; i<split.length; i++){
			sb.append(tab);
			sb.append(split[i].replace(',','_'));
			tab="\t";
		}
		for(int i=1; i<split2.length; i++){
			sb.append(tab);
			sb.append(split2[i].replace(',','_'));
			tab="\t";
		}
		return sb;
	}
	
	private HashMap<String, String[]> hash(FileFormat ff){
		final HashMap<String, String[]> table=new HashMap<String, String[]>();
		final TextFile tf;
		{
			tf=new TextFile(ff);
			if(verbose){outstream.println("Started tf");}
		}
		{
			String line;
			while((line=tf.nextLine())!=null){
				String[] split=line.split(delimiter);
				table.put(split[0], split); //n split[0] is always safe: String.split on any input yields length>=1.
			}
		}
		//NOTE [driver/MergeBigelow#002] LOW/resource/dev: this in2 TextFile is read to EOF then NEVER closed. TextFile does
		//NOT auto-close at EOF (readLine returns null without close(), verified TextFile L374), so the in2 InputStream leaks
		//until process exit, and TextFile.close()'s errorState (errorState|=ReadWrite.finishReading, L311) is discarded for in2.
		//Asymmetric with the main-loop tf which IS closed (process() L134). Should be `tf.close();` before return (dropping its
		//return is fine — the in1 close already feeds the propagated errorState). Dev one-off so LOW; matters only if reused.
		return table;
	}
	
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/

	private String in1=null;
	private String in2=null;
	private String out1=null;
	
	private String delimiter="\t";
	private HashMap<String, String[]> table;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/

	private final FileFormat ffin1;
	private final FileFormat ffin2;
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
