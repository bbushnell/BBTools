package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;

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
 * Combines multiple sketches into a single sketch.
 * 
 * @author Brian Bushnell
 * @date July 23, 2018
 *
 */
public class MergeSketch extends SketchObject {
	
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
		MergeSketch x=new MergeSketch(args);
		
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
	public MergeSketch(String[] args){
		
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
		parser.out1="stdout.txt";
		
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
			
			else if(a.equals("name") || a.equals("taxname")){
				outTaxName=b;
			}else if(a.equals("name0")){
				outName0=b;
			}else if(a.equals("fname")){
				outFname=b;
			}else if(a.equals("taxid") || a.equals("tid")){
				outTaxID=Integer.parseInt(b);
			}else if(a.equals("spid")){
				outSpid=Integer.parseInt(b);
			}else if(a.equals("imgid")){
				outImgID=Integer.parseInt(b);
			}else if((a.startsWith("meta_") || a.startsWith("mt_")) && b!=null){
				if(outMeta==null){outMeta=new ArrayList<String>();}
				int underscore=a.indexOf('_', 0);
				outMeta.add(a.substring(underscore+1)+":"+b);
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
		outMeta=SketchObject.fixMeta(outMeta);
		
		blacklist=null;
		
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
		
		ffout=FileFormat.testOutput(outSketch, FileFormat.SKETCH, null, false, overwrite, append, false);
		if(ffout!=null && !ffout.stdio() && !defaultParams.setColors){defaultParams.printColors=false;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, outSketch)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+outSketch+"\n");
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in.toArray(new String[0]))){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		tool=new SketchTool(targetSketchSize, defaultParams);
		
//		assert(false) : defaultParams.toString()+"\n"+k+", "+amino+", "+HASH_VERSION;
		if(verbose){
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
//		outstream.println(inSketches.get(0));
		
		ByteBuilder bb=new ByteBuilder();
		
		int sizeOut=(int)(Sketch.AUTOSIZE ? sum : Tools.min(Sketch.targetSketchSize, sum));
		{
			Sketch.AUTOSIZE=false;
			Sketch.targetSketchSize=sizeOut;
			Sketch.maxGenomeFraction=1;
		}
		SketchHeap heap=new SketchHeap(sizeOut, 0, tool.trackCounts);
		for(Sketch sk : inSketches){
			heap.add(sk);
		}
		heap.genomeSizeKmers=Tools.max(heap.genomeSizeKmers, sizeOut);
		ArrayList<String> meta=inSketches.get(0).meta;
		if(meta==null){meta=outMeta;}
		else if(outMeta!=null){meta.addAll(outMeta);}
		Sketch union=new Sketch(heap, false, tool.trackCounts, outMeta);

		if(outTaxName!=null){union.setTaxName(outTaxName);}
		if(outFname!=null){union.setFname(outFname);}
		if(outName0!=null){union.setName0(outName0);}

		if(outTaxID>=0){union.taxID=(outTaxID);}
		if(outSpid>=0){union.spid=(outSpid);}
		if(outImgID>=0){union.imgID=(outImgID);}
		
		if(outSketch!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(outSketch, overwrite, append, true, FileFormat.SKETCH);
			bsw.start();
			union.toBytes(bb);
			bsw.print(bb);
			bb.clear();
			bsw.poisonAndWait();
			errorState|=bsw.errorState;
			t.stop();
			outstream.println("Wrote "+1+" sketch of total size "+union.length()+" in \t"+t);
		}
		
		t.stop();
//		outstream.println("\nRan "+(inSketches.size()*refSketches.size())+" comparisons in \t"+t);
		ttotal.stop();
		outstream.println("Total Time: \t"+ttotal);
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
	
	private ArrayList<String> in=new ArrayList<String>();
	
	private String outSketch=null;
	
	private final SketchTool tool;
	
	private ArrayList<Sketch> inSketches;
	
	/*Override metadata */
	private String outTaxName=null;
	private String outFname=null;
	private String outName0=null;
	private int outTaxID=-1;
	private long outSpid=-1;
	private long outImgID=-1;
	private ArrayList<String> outMeta=null;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary output file */
	private final FileFormat ffout;
	
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
