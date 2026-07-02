package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;

/**
 * Correlates identity values from two tab-separated matrix files.
 * Reads corresponding entries from two identity matrices and outputs paired values
 * for statistical analysis or correlation studies.
 *
 * @author Brian Bushnell
 * @date Nov 21, 2014
 */
public class CorrelateIdentity {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Program entry point.
	 * Creates a CorrelateIdentity instance and executes the correlation process.
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
	 * Constructor that parses command-line arguments and validates input/output files.
	 * Sets up compression defaults and validates file accessibility.
	 * @param args Command line arguments including input files (in1, in2) and output file (out)
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
		//Extracts the strictly-lower triangle + diagonal (j:1..i), skipping col 0 (row-label column) and the redundant
		//upper triangle of a symmetric identity matrix. Comprehension: pairs matrix1[i][j] with matrix2[i][j] as (x,y).
		for(int i=0; i<matrix1.length; i++){
			//TODO: Possible bug [driver/CorrelateIdentity#002] LOW/precondition: matrix2 is indexed by matrix1's dimensions
			//(loop bound matrix1.length, and matrix2[i][j]) with NO shape check that matrix2.length>=matrix1.length or that
			//matrix2[i] has >=i+1 cols. If in2 is a smaller/ragged matrix than in1, matrix2[i][j] throws AIOOBE. Correlating
			//two identity matrices assumes same taxa set / same shape, so this is a loud precondition crash — LOW, not fixed.
			for(int j=1; j<=i; j++){
				list.add(new String[] {matrix1[i][j], matrix2[i][j]});
			}
		}

		Collections.shuffle(list); //NOTE: uses Collections' own Random, NOT sampleseed => sampleseed is doubly-dead (see #003). Order is non-reproducible.

		//TODO: Possible bug [driver/CorrelateIdentity#001] LOW: out= omitted => new TextStreamWriter(null,...) => testOutput(null)
		//returns null => FileFormat-ctor NPE (TextStreamWriter.java:41). No out-required guard (in files ARE guarded L100-105).
		//Same shape as FilterAssemblySummary#001, and here there is not even the vestigial null-guard. Escalated, not fixed.
		TextStreamWriter tsw=new TextStreamWriter(out, overwrite, append, true);
		tsw.start();
		for(String[] pair : list){
			tsw.print(pair[0]+"\t"+pair[1]+"\n");
		}
		tsw.poisonAndWait(); //NOTE [#001 cont'd]: poisonAndWait() return DROPPED and this class has NO errorState field/terminal throw => a write failure is silently swallowed and the tool still exits 0 (pattern b). LOW.
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public String in1, in2;
	/** Output file for correlated value pairs */
	public String out;
	
//	private Random randy=shared.Shared.random();

	//TODO: Possible bug [driver/CorrelateIdentity#003] LOW/dead-code: samplerate, sampleseed, columnLength are ALL
	//parse-only/never-read in process() — subsampling was never implemented (`//private Random randy` above is commented
	//out; Collections.shuffle uses its own Random). So `samplerate=0.1` is silently accepted (even assert-validated L85)
	//yet has ZERO effect — the user gets every pair, not a sample. Not documented in matrixtocolumns.sh usage, so it's a
	//silent-accept-and-ignore of an undocumented arg rather than a broken doc promise. Also sampleseed is a `float` field
	//assigned Long.parseLong(b) (L87) => long→float widening (cosmetic, unused anyway). columnLength never referenced.
	private float samplerate=1;
	private float sampleseed=-1;
	private int columnLength=Integer.MAX_VALUE;
	private boolean overwrite=true;
	private boolean append=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Controls verbose output messages */
	public static final boolean verbose=false; //123
	
	private static PrintStream outstream=System.err;
	
}
