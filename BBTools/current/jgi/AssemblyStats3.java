package jgi;

import java.io.File;
import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Generates some stats from multiple files.
 * Uses the new Assembly class.
 * @author Brian Bushnell
 * @date January 21, 2025
 *
 */
public class AssemblyStats3 {

	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		AssemblyStats3 x=new AssemblyStats3(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
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
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
				//Set a variable here
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
		
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, true, false, false);
	}
	
	void process(Timer t){

		final ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout1);
		if(bsw!=null) {bsw.println(makeHeader());}
		
		for(String s : in) {
			processInner(s, bsw);
		}
		
		if(bsw!=null) {bsw.poisonAndWait();}
		if(verbose){outstream.println("Finished.");}
		
		t.stop();
		Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 12);
	}
	
	public static String makeHeader() {
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
	
	void processInner(String fname, ByteStreamWriter bsw) {
		Assembly a=new Assembly(fname);
		if(bsw==null) {return;}
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
	
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> in=new ArrayList<String>();
	private String out1="stdout.txt";
	
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	private long linesProcessed=0, bytesProcessed=0;
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
