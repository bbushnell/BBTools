package covid;

import java.io.File;
import java.io.PrintStream;
import java.util.LinkedHashSet;

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
import structures.ByteBuilder;

/**
 * @author Brian Bushnell
 * @date April 5, 2020
 *
 */
public class SummarizeCoverage {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		SummarizeCoverage x=new SummarizeCoverage(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public SummarizeCoverage(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, /*getClass()*/null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		{//Parse the arguments
			final Parser parser=parse(args);
			overwrite=parser.overwrite;
			append=parser.append;
			
			if(in.isEmpty() && parser.in1!=null){
				in.add(parser.in1);
			}

			out1=parser.out1;
		}
		
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program

		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		Parser parser=new Parser();
		parser.out1=out1;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("refbases") || a.equals("refsize") || a.equals("reflen")){
				refLen=Long.parseLong(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(b==null && new File(arg).exists()){
				in.add(arg);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
//		in1=Tools.fixExtension(in1); //TODO
		if(in.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
//		if(!Tools.testForDuplicateFiles(true, out1, in.toArray(new String[0]))){//TODO
//			throw new RuntimeException("\nSome file names were specified multiple times.\n");
//		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
//		if(!ByteFile.FORCE_MODE_BF2){
//			ByteFile.FORCE_MODE_BF2=false;
//			ByteFile.FORCE_MODE_BF1=true;
//		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	void process(Timer t){
		
		ByteStreamWriter bsw=makeBSW(ffout1);
		
		if(bsw!=null){
//			bsw.println("#Sample\tAvgCov\t%<1x\t%<2x\t%<3x\t%<4x\t%<5x\t%<10x\t%<20x");
			String header="#Sample\tAvgCov\t%>=1x\t%>=2x\t%>=3x\t%>=4x\t%>=5x\t%>=10x\t%>=20x";
			bsw.println(header);
			linesOut++;
			bytesOut+=header.length();
		}
		
		for(String fname : in){
			FileFormat ffin=FileFormat.testInput(fname, FileFormat.TXT, null, true, true);
			ByteFile bf=ByteFile.makeByteFile(ffin);
			String name=ReadWrite.stripToCore(fname);
			if(name.toLowerCase().endsWith("_basecov")){
				name=name.substring(0, name.lastIndexOf('_'));
			}
			processInner(bf, bsw, name);
			errorState|=bf.close();
		}
		
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		
		t.stop();
		
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		
		outstream.println();
		outstream.println("Valid Lines:       \t"+linesOut);
		outstream.println("Invalid Lines:     \t"+(linesProcessed-linesOut));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void processInner(ByteFile bf, ByteStreamWriter bsw, String name){
		byte[] line=bf.nextLine();
		final int max=20;
		long[] hist=new long[max+1];
		long covSum=0;
		while(line!=null){
			if(line.length>0){
				if(maxLines>0 && linesProcessed>=maxLines){break;}
				linesProcessed++;
				bytesProcessed+=(line.length+1);
				
				final boolean valid=(line[0]!='#');
				
				if(valid){
					int tab=Tools.lastIndexOf(line, (byte)'\t');
					int cov=Parse.parseInt(line, tab+1);
					covSum+=cov;
					hist[Tools.min(max, cov)]++;
				}else{
					//Parse it if desired
				}
			}
			line=bf.nextLine();
		}
		long len=shared.Vector.sum(hist);
		if(refLen>0){
			assert(refLen>=len);
			len=refLen;
		}
		final double mult=1.0/Tools.max(1, len);
		
		final ByteBuilder bb=new ByteBuilder();
		//Change to cumulative
		for(int i=hist.length-1; i>0; i--){
			hist[i-1]+=hist[i];
		}
		bb.append(name).tab();
		bb.append(covSum*mult, 2).tab();
		bb.append(hist[1]*mult*100, 2).tab();
		bb.append(hist[2]*mult*100, 2).tab();
		bb.append(hist[3]*mult*100, 2).tab();
		bb.append(hist[4]*mult*100, 2).tab();
		bb.append(hist[5]*mult*100, 2).tab();
		bb.append(hist[10]*mult*100, 2).tab();
		bb.append(hist[20]*mult*100, 2).tab();
		bb.nl();
		linesOut++;
		bytesOut+=bb.length;
		if(bsw!=null){bsw.print(bb);}
	}
	
	private static ByteStreamWriter makeBSW(FileFormat ff){
		if(ff==null){return null;}
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		return bsw;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private LinkedHashSet<String> in=new LinkedHashSet<String>();
	private String out1="stdout.txt";
	
	/*--------------------------------------------------------------*/
	
	private long linesProcessed=0;
	private long linesOut=0;
	private long bytesProcessed=0;
	private long bytesOut=0;
	
	private long maxLines=Long.MAX_VALUE;
	
	private long refLen=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
//	private final FileFormat ffin1;
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
