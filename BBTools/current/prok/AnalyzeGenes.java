package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

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
import structures.IntList;

/**
 * This class is designed to analyze paired prokaryotic fna and gff files
 * to calculate the patterns in coding and noncoding frames, start and stop sites.
 * It outputs a pgm file.
 * @author Brian Bushnell
 * @date Sep 27, 2018
 *
 */
public class AnalyzeGenes {
	
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
		AnalyzeGenes x=new AnalyzeGenes(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public AnalyzeGenes(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null/*getClass()*/, false);
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

			out=parser.out1;
		}
		
		if(alignRibo){
			//Load sequences
			ProkObject.loadConsensusSequenceFromFile(false, false);
		}
		
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program
		
		//Determine how many threads may be used
		threads=Tools.min(fnaList.size(), Shared.threads(), Tools.max(32, Shared.CALC_LOGICAL_PROCESSORS()/2));
		
		ffout=FileFormat.testOutput(out, FileFormat.PGM, null, true, overwrite, append, false);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		
		Parser parser=new Parser();
		parser.overwrite=overwrite;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

//			outstream.println(arg+", "+a+", "+b);
			if(PGMTools.parseStatic(arg, a, b)){
				//do nothing
			}else if(a.equals("in") || a.equals("infna") || a.equals("fnain") || a.equals("fna") || a.equals("ref")){
				assert(b!=null);
				Tools.addFiles(b, fnaList);
			}else if(a.equals("gff") || a.equals("ingff") || a.equals("gffin")){
				assert(b!=null);
				Tools.addFiles(b, gffList);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ReadWrite.verbose=verbose;
			}else if(a.equals("alignribo") || a.equals("align")){
				alignRibo=Parse.parseBoolean(b);
			}else if(a.equals("adjustendpoints")){
				adjustEndpoints=Parse.parseBoolean(b);
			}
			
			else if(ProkObject.parse(arg, a, b)){}
			
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(arg.indexOf('=')<0 && new File(arg).exists() && FileFormat.isFastaFile(arg)){
				fnaList.add(arg);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		if(gffList.isEmpty()){
			for(String s : fnaList){
				String prefix=ReadWrite.stripExtension(s);
				String gff=prefix+".gff";
				File f=new File(gff);
				if(!f.exists()){
					String gz=gff+".gz";
					f=new File(gz);
					assert(f.exists() && f.canRead()) : "Can't read file "+gff;
					gff=gz;
				}
				gffList.add(gff);
			}
		}
		assert(gffList.size()==fnaList.size()) : "Number of fna and gff files do not match: "+fnaList.size()+", "+gffList.size();
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		fnaList=Tools.fixExtension(fnaList);
		gffList=Tools.fixExtension(gffList);
		if(fnaList.isEmpty()){throw new RuntimeException("Error - at least one input file is required.");}
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
		}
		
		//Ensure input files can be read
		ArrayList<String> foo=new ArrayList<String>();
		foo.addAll(fnaList);
		foo.addAll(gffList);
		if(!Tools.testInputFiles(false, true, foo.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		foo.add(out);
		if(!Tools.testForDuplicateFiles(true, foo.toArray(new String[0]))){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	void process(Timer t){
		
		final GeneModel pgm;
		if(Shared.threads()<2 || fnaList.size()<2){
			pgm=makeModelST();
		}else{
			pgm=spawnThreads();
		}
		
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout);
		
		ByteBuilder bb=new ByteBuilder();
		pgm.appendTo(bb);
		bytesOut+=bb.length;
		
		if(bsw!=null){
			bsw.addJob(bb);
			errorState|=bsw.poisonAndWait();
		}
		
		t.stop();
		
		outstream.println(timeReadsBasesGenesProcessed(t, pgm.readsProcessed, pgm.basesProcessed, pgm.genesProcessed, pgm.filesProcessed, 8));
		
		outstream.println();
		outstream.println(typesProcessed(pgm, 12));
		
		//outstream.println("Bytes Out:         \t"+bytesOut);
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private static String timeReadsBasesGenesProcessed(Timer t, long readsProcessed, long basesProcessed, long genesProcessed, long filesProcessed, int pad){
		return ("Time:                         \t"+t+"\n"+readsBasesGenesProcessed(t.elapsed, readsProcessed, basesProcessed, genesProcessed, filesProcessed, pad));
	}
	
	private static String readsBasesGenesProcessed(long elapsed, long reads, long bases, long genes, long files, int pad){
		double rpnano=reads/(double)elapsed;
		double bpnano=bases/(double)elapsed;
		double gpnano=genes/(double)elapsed;
		double fpnano=files/(double)elapsed;

		String rstring=Tools.padKMB(reads, pad);
		String bstring=Tools.padKMB(bases, pad);
		String gstring=Tools.padKMB(genes, pad);
		String fstring=Tools.padKMB(files, pad);
		ByteBuilder sb=new ByteBuilder();
		sb.append("Files Processed:    ").append(fstring).append(Tools.format(" \t%.2f  files/sec", fpnano*1000000000)).append('\n');
		sb.append("Sequences Processed:").append(rstring).append(Tools.format(" \t%.2fk seqs/sec", rpnano*1000000)).append('\n');
		sb.append("Genes Processed:    ").append(gstring).append(Tools.format(" \t%.2fk genes/sec", gpnano*1000000)).append('\n');
		sb.append("Bases Processed:    ").append(bstring).append(Tools.format(" \t%.2fm bases/sec", bpnano*1000));
		return sb.toString();
	}
	
	private static String typesProcessed(GeneModel pgm, int pad){
		
		ByteBuilder sb=new ByteBuilder();
		sb.append("CDS:   "+Tools.padLeft(pgm.statsCDS.lengthCount, pad)).nl();
		sb.append("tRNA:  "+Tools.padLeft(pgm.statstRNA.lengthCount, pad)).nl();
		sb.append("16S:   "+Tools.padLeft(pgm.stats16S.lengthCount, pad)).nl();
		sb.append("23S:   "+Tools.padLeft(pgm.stats23S.lengthCount, pad)).nl();
		sb.append("5S:    "+Tools.padLeft(pgm.stats5S.lengthCount, pad)).nl();
		sb.append("18S:   "+Tools.padLeft(pgm.stats18S.lengthCount, pad));
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	//TODO: Process each file in a thread.
	private GeneModel makeModelST(){
		GeneModel pgmSum=new GeneModel(true);
		
		for(int i=0; i<fnaList.size(); i++){
			String fna=fnaList.get(i);
			String gff=gffList.get(i);
			pgmSum.process(fna, gff);
		}
		return pgmSum;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private GeneModel spawnThreads(){
		
		//Do anything necessary prior to processing
		
		final AtomicInteger aint=new AtomicInteger(0);
		
		//Fill a list with FileThreads
		ArrayList<FileThread> alpt=new ArrayList<FileThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new FileThread(aint));
		}
		
		//Start the threads
		for(FileThread pt : alpt){
			pt.start();
		}
		
		//Wait for threads to finish
		GeneModel pgm=waitForThreads(alpt);
		
		//Do anything necessary after processing
		return pgm;
	}
	
	private GeneModel waitForThreads(ArrayList<FileThread> alpt){
		
		GeneModel pgm=new GeneModel(false);
		
		//Wait for completion of all threads
		boolean success=true;
		for(FileThread pt : alpt){
			
			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			
			//Accumulate per-thread statistics
			pgm.add(pt.pgm);
			
			success&=pt.success;
			errorState|=pt.errorStateT;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		return pgm;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class FileThread extends Thread {
		
		FileThread(AtomicInteger fnum_){
			fnum=fnum_;
			pgm=new GeneModel(true);
		}
		
		@Override
		public void run(){
			for(int i=fnum.getAndIncrement(); i<fnaList.size(); i=fnum.getAndIncrement()){
				String fna=fnaList.get(i);
				String gff=gffList.get(i);
				errorStateT=pgm.process(fna, gff)|errorState;
//				System.err.println("Processed "+fna+" in "+this.toString());
			}
			success=true;
		}
		
		private final AtomicInteger fnum;
		private final GeneModel pgm;
		boolean errorStateT=false;
		boolean success=false;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private ArrayList<String> fnaList=new ArrayList<String>();
	private ArrayList<String> gffList=new ArrayList<String>();
	private IntList taxList=new IntList();
	private String out=null;
	
	/*--------------------------------------------------------------*/
	
	private long bytesOut=0;
	static boolean alignRibo=true;
	static boolean adjustEndpoints=true;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffout;
	private final int threads;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}

