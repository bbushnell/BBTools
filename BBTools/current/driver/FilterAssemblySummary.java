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
import tax.TaxFilter;

/**
 * @author Brian Bushnell
 * @date May 1, 2016
 *
 */
public class FilterAssemblySummary {
	
	public static void main(String[] args){
		Timer t=new Timer();
		FilterAssemblySummary x=new FilterAssemblySummary(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public FilterAssemblySummary(String[] args){
		
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
			}else if(TaxFilter.validArgument(a)){
				//do nothing
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
		
		//Make the actual filter
		filter=TaxFilter.makeFilter(args);
	}
	
	void process(Timer t){
		
		final TextFile tf;
		{
			tf=new TextFile(ffin1);
			if(verbose){outstream.println("Started tf");}
		}
		
		final TextStreamWriter tsw;
		{
			tsw=new TextStreamWriter(ffout1);
			tsw.start();
			if(verbose){outstream.println("Started tsw");}
		}

		long linesProcessed=0;
		long linesRetained=0;
		long charsProcessed=0;
		
		{
			String line;
			while((maxReads<0 || linesProcessed<maxReads) && (line=tf.nextLine())!=null){
				linesProcessed++;
				charsProcessed+=line.length();
				String result=processLine(line);
				if(result!=null){
					linesRetained++;
					if(tsw!=null){tsw.println(result);}
				}
			}
		}
		
		errorState|=tsw.poisonAndWait();
		errorState|=tf.close();
		
		t.stop();
		
		String kpstring=(linesRetained<100000 ? ""+linesRetained : linesRetained<100000000 ? (linesRetained/1000)+"k" : (linesRetained/1000000)+"m");
		while(kpstring.length()<8){kpstring=" "+kpstring;}
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, charsProcessed, 8));
		outstream.println("Lines Retained:     "+kpstring);
		
		
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	private String processLine(String line){
//		System.out.println("Processing line "+line);
		if(line.startsWith("#")){return null;}
		String[] split=line.split("\t");
		assert(split.length>6) : split.length+"\n"+"'"+line+"'";
		String id=split[6];
		int number=Integer.parseInt(id);
//		System.out.println("Found number "+number+" from "+id+"; node: "+filter.tree().getNode(number));
		boolean b=filter.passesFilter((int)number);
//		System.out.println("passesFilter? "+b);
		return b ? line : null;
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
	
	/** The actual filter */
	private final TaxFilter filter;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
