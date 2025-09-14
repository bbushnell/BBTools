package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Tools;

/**
 * @author Brian Bushnell
 * @date Nov 21, 2014
 *
 */
public class CorrelateIdentity {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Create a new CorrelateIdentity instance
		CorrelateIdentity x=new CorrelateIdentity(args);
		
		///And run it
		x.process();
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CorrelateIdentity(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		/* Set global defaults */
		ReadWrite.ZIPLEVEL=6;
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_PIGZ=true;
		
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("out") || a.equals("out1")){
				out=b;
			}else if(a.equals("samplerate")){
				samplerate=Float.parseFloat(b);
				assert(samplerate<=1f && samplerate>=0f) : "samplerate="+samplerate+"; should be between 0 and 1";
			}else if(a.equals("sampleseed")){
				sampleseed=Long.parseLong(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
		}
		

		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			throw new RuntimeException("\nCan't write to some output files; overwrite="+overwrite+"\n");
		}
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		if(!Tools.testForDuplicateFiles(true, in1, in2, out)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}

		assert(in1==null || in1.toLowerCase().startsWith("stdin") || in1.toLowerCase().startsWith("standardin") || new File(in1).exists()) : "Can't find "+in1;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void process(){
		final String[][] matrix1, matrix2;
		
		{
			TextFile tf=new TextFile(in1);
			String[] s=tf.toStringLines();
			tf.close();
			matrix1=tf.doublesplitWhitespace(s, true);
		}
		
		{
			TextFile tf=new TextFile(in2);
			String[] s=tf.toStringLines();
			tf.close();
			matrix2=tf.doublesplitWhitespace(s, true);
		}
		
		ArrayList<String[]> list=new ArrayList<String[]>();
		for(int i=0; i<matrix1.length; i++){
			for(int j=1; j<=i; j++){
				list.add(new String[] {matrix1[i][j], matrix2[i][j]});
			}
		}
		
		Collections.shuffle(list);
		
		TextStreamWriter tsw=new TextStreamWriter(out, overwrite, append, true);
		tsw.start();
		for(String[] pair : list){
			tsw.print(pair[0]+"\t"+pair[1]+"\n");
		}
		tsw.poisonAndWait();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Input files */
	public String in1, in2;
	/** Output file */
	public String out;
	
//	private Random randy=new Random();

	private float samplerate=1;
	private float sampleseed=-1;
	private int columnLength=Integer.MAX_VALUE;
	private boolean overwrite=true;
	private boolean append=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Verbose messages */
	public static final boolean verbose=false; //123
	
	/** Print messages to this stream */
	private static PrintStream outstream=System.err;
	
}
