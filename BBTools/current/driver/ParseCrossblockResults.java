package driver;

import java.io.File;
import java.io.PrintStream;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date June 9, 2016
 *
 */
public class ParseCrossblockResults {
	
	public static void main(String[] args){
		Timer t=new Timer();
		ParseCrossblockResults x=new ParseCrossblockResults(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public ParseCrossblockResults(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		
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

		ffin1=FileFormat.testInput(in1, FileFormat.TEXT, null, false, false);
	}
	
	void process(Timer t){
		
		final TextFile tf;
		{
			tf=new TextFile(ffin1);
			if(verbose){outstream.println("Started tf");}
		}
		
		long linesProcessed=0;
		long charsProcessed=0;
		
		{
			String line;
			while((maxReads<0 || linesProcessed<maxReads) && (line=tf.nextLine())!=null){
				linesProcessed++;
				charsProcessed+=line.length();
				if(!line.startsWith("#")){
					processLine(line);
				}
			}
		}
		errorState|=tf.close();

		if(ffout1!=null){
			final TextStreamWriter tsw;
			{
				tsw=new TextStreamWriter(ffout1);
				tsw.start();
				if(verbose){outstream.println("Started tsw");}
				errorState|=tsw.poisonAndWait();
			}
		}
		
		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, charsProcessed, 8));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	private void processLine(String line){
		ResultsLine rl=new ResultsLine(line);
		if(rl.removed){
			basesDiscarded+=rl.length;
			contigsDiscarded++;
		}else{
			basesKept+=rl.length;
			contigsKept++;
		}
	}
	
	
	
	/*--------------------------------------------------------------*/
	
	private static class ResultsLine{
		
		public ResultsLine(String s){
			String[] split=s.split("\t");
			length=Integer.parseInt(split[3]);
			removed=Integer.parseInt(split[2])==1;
		}
		
		final int length;
		final boolean removed;
		
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;

	private long basesKept=0;
	private long basesDiscarded=0;
	private long contigsKept=0;
	private long contigsDiscarded=0;

	public long basesKept(){return basesKept;}
	public long basesDiscarded(){return basesDiscarded;}
	public long contigsKept(){return contigsKept;}
	public long contigsDiscarded(){return contigsDiscarded;}
	
	public long contigs(){return contigsKept+contigsDiscarded;}
	public long bases(){return basesKept+basesDiscarded;}
	
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
