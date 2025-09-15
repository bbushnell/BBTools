package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashSet;

import fileIO.ByteFile;
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
import tracker.ReadStats;

/**
 * Generates smaller sketches from input sketches.
 * 
 * @author Brian Bushnell
 * @date July 23, 2018
 *
 */
public class SubSketch extends SketchObject {
	
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
		
		final boolean oldUnpigz=ReadWrite.USE_UNPIGZ;
		final int oldBufLen=Shared.bufferLen();
		
		//Create an instance of this class
		SubSketch x=new SubSketch(args);
		
		//Run the object
		x.process(t);
		
		ReadWrite.USE_UNPIGZ=oldUnpigz;
		Shared.setBufferLen(oldBufLen);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
		
		assert(!x.errorState) : "This program ended in an error state.";
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public SubSketch(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		ReadWrite.USE_UNPIGZ=true;
		KILL_OK=true;
		
		//Create a parser object
		Parser parser=new Parser();
		
		defaultParams.printRefFileName=true;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("in")){
				addFiles(b, in);
			}else if(a.equals("files")){
				files=Integer.parseInt(b);
			}else if(parseSketchFlags(arg, a, b)){
				//Do nothing
			}else if(defaultParams.parse(arg, a, b)){
				//Do nothing
			}
//			else if(a.equals("size")){
//				size=Parse.parseIntKMG(b);
//			}
			
			else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}
			
			else if(a.equals("out") || a.equals("outsketch") || a.equals("outs") || a.equals("sketchout") || a.equals("sketch")){
				outSketch=b;
			}
			
			else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}
			
			else if(b==null && new File(arg).exists()){
				in.add(arg);
			}
			
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		assert(targetSketchSize>0) : "Must set size.";
		
		{//Expand # symbol
			LinkedHashSet<String> expanded=new LinkedHashSet<String>();
			for(String s : in){SketchSearcher.addFiles(s, expanded);}
			in.clear();
			in.addAll(expanded);
		}
		
		postParse();
		
		{//Process parser fields
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
		}
		
		//Ensure there is an input file
		if(in.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(!Tools.testOutputFiles(overwrite, append, false, outSketch)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+outSketch+"\n");
		}
//		assert(false) : ffout;
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in.toArray(new String[0]))){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		tool=new SketchTool(targetSketchSize, defaultParams);
		
//		assert(false) : defaultParams.toString()+"\n"+k+", "+amino+", "+HASH_VERSION;
		if(verbose || true){
			if(useWhitelist){outstream.println("Using a whitelist.");}
			if(blacklist!=null){outstream.println("Using a blacklist.");}
		}
		
		defaultParams.postParse(false, false);
		allowMultithreadedFastq=(in.size()==1 && Shared.threads()>2);
		if(!allowMultithreadedFastq){Shared.capBufferLen(40);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void process(Timer t){
		Timer ttotal=new Timer();
		
		t.start();
		inSketches=tool.loadSketches_MT(defaultParams, in);
		final int numLoaded=(inSketches.size());
		long sum=0;
		for(Sketch sk : inSketches){
			sum+=sk.length();
		}
		t.stop();
		outstream.println("Loaded "+numLoaded+" sketch"+(numLoaded==1 ? "" : "es")+" of total size "+sum+" in "+t);
		t.start();
		if(verbose && numLoaded>0){
			System.err.println("First sketch:\n"+inSketches.get(0));
		}
//		outstream.println(inSketches.get(0));
		
		int sizeOut=Sketch.targetSketchSize;
		{
			if(Sketch.SET_TARGET_SIZE){Sketch.AUTOSIZE=false;}
			Sketch.targetSketchSize=sizeOut;
			Sketch.maxGenomeFraction=1;
		}
		
		if(outSketch!=null && outSketch.indexOf('#')>=1 && files>1){
			ByteStreamWriter[] bswArray=new ByteStreamWriter[files];
			for(int i=0; i<files; i++){
				FileFormat ffout=FileFormat.testOutput(outSketch.replace("#", ""+i), FileFormat.SKETCH, null, false, overwrite, append, false);
				ByteStreamWriter bsw=new ByteStreamWriter(ffout);
				bsw.start();
				bswArray[i]=bsw;
			}

			processInner(inSketches, bswArray);

			for(ByteStreamWriter bsw : bswArray){
				bsw.poisonAndWait();
				errorState|=bsw.errorState;
			}
		}else{
			FileFormat ffout=FileFormat.testOutput(outSketch, FileFormat.SKETCH, null, false, overwrite, append, false);
			ByteStreamWriter bsw=null;
			if(ffout!=null){
				bsw=new ByteStreamWriter(ffout);
				bsw.start();
			}

			processInner(inSketches, bsw);

			if(bsw!=null){
				bsw.poisonAndWait();
				errorState|=bsw.errorState;
			}
		}
		
		t.stop();
		if(blacklist!=null){outstream.println("Evicted "+blackKeys+" blacklisted keys.");}
		outstream.println("Wrote "+sketchesOut+" sketches of total size "+keysOut+" in "+t);

		t.stop();
		ttotal.stop();
		outstream.println("Total Time: \t"+ttotal);
	}
	
	void processInner(ArrayList<Sketch> sketches, ByteStreamWriter bsw){
		ByteBuilder bb=new ByteBuilder();
		for(Sketch sk : sketches){
			final int target=Sketch.AUTOSIZE ? toSketchSize(sk.genomeSizeBases, sk.genomeSizeKmers, sk.genomeSizeEstimate(), targetSketchSize) : targetSketchSize;
//			if(!defaultParams.trackCounts()){sk.keyCounts=null;}
			if(blacklist!=null){blackKeys+=sk.applyBlacklist();}
			if(sk.length()>target){
				sk.resize(target);
				if(verbose){System.err.println("Resized to:\n"+sk);}
			}
			if(sk.length()>=minSketchSize){
				keysOut+=sk.length();
				sketchesOut++;
				sk.toBytes(bb);
				if(verbose){System.err.println("toBytes:\n"+bb);}
				if(bsw!=null){bsw.print(bb);}
				bb.clear();
			}
		}
	}
	
	void processInner(ArrayList<Sketch> sketches, ByteStreamWriter bswa[]){
		ByteBuilder bb=new ByteBuilder();
		for(Sketch sk : sketches){
			//final int target=Sketch.AUTOSIZE ? toSketchSize(sk.genomeSizeBases, sk.genomeSizeKmers, sk.genomeSizeEstimate(), targetSketchSize) : targetSketchSize;
//			if(!defaultParams.trackCounts()){sk.keyCounts=null;}
			if(blacklist!=null){blackKeys+=sk.applyBlacklist();}
			
			//Calculating target after applying blacklist gives better consistency with actual usage
			final int target=Sketch.AUTOSIZE ? toSketchSize(sk.genomeSizeBases, sk.genomeSizeKmers, sk.genomeSizeEstimate(), targetSketchSize) : targetSketchSize;
			
			if(sk.length()>target){
				sk.resize(target);
				if(verbose){System.err.println("Resized to:\n"+sk);}
			}
			if(sk.length()>=minSketchSize){
				keysOut+=sk.length();
				sketchesOut++;
				
				if(bswa!=null){
					ByteStreamWriter bsw=bswa[sk.sketchID%files];
					if(sk.fname()!=null && sk.fname().endsWith(".sketch")){sk.setFname(bsw.fname);}
					sk.toBytes(bb);//This is the time-limiting factor; could be multithreaded.
					if(verbose){System.err.println("toBytes:\n"+bb);}
					bsw.print(bb);
				}
				bb.clear();
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static boolean addFiles(String a, Collection<String> list){
		int initial=list.size();
		if(a==null){return false;}
		File f=null;
		if(a.indexOf(',')>=0){f=new File(a);}
		if(f==null || f.exists()){
			list.add(a);
		}else{
			for(String s : a.split(",")){
				list.add(s);
			}
		}
		return list.size()>initial;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private LinkedHashSet<String> in=new LinkedHashSet<String>();
	
	private String outSketch=null;
	
	private final SketchTool tool;
	
	private ArrayList<Sketch> inSketches;

	private long keysOut=0;
	private long sketchesOut=0;
	private long blackKeys=0;
	
	private int files=31;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Don't print caught exceptions */
	public static boolean suppressErrors=false;
	
}
