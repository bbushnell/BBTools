package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import kmer.AbstractKmerTableSet;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import tax.TaxFilter;
import tax.TaxTree;
import tracker.ReadStats;

/**
 * Compares one or more input sketches to a set of reference sketches.
 * 
 * @author Brian Bushnell
 * @date July 29, 2016
 *
 */
public class CompareSketch extends SketchObject {
	

	
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
		
		final int oldBufLen=Shared.bufferLen();
		
		//Create an instance of this class
		CompareSketch x=new CompareSketch(args);
		
		//Run the object
		x.process(t);

		Shared.setBufferLen(oldBufLen);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
		
		alignerPool.poison();
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CompareSketch(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null, false);
//			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
			silent=PreParser.silent;
			if(silent){AbstractKmerTableSet.DISPLAY_PROGRESS=false;}
		}
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		KILL_OK=true;
		TaxFilter.REQUIRE_PRESENT=false;
		defaultParams.mode=PER_FILE;
		
		//Create a parser object
		Parser parser=new Parser();
		parser.out1="stdout.txt";
		
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
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}else if(a.equals("alltoall") || a.equals("ata")){
				allToAll=Parse.parseBoolean(b);
			}else if(a.equals("skipcompare") || a.equals("sketchonly")){
				skipCompare=Parse.parseBoolean(b);
			}else if(a.equals("compareself") || a.equals("includeself")){
				compareSelf=Parse.parseBoolean(b);
			}else if(a.equals("printmemory")){
				printMemory=Parse.parseBoolean(b);
			}else if(a.equals("parsesubunit")){
				SketchMaker.parseSubunit=Parse.parseBoolean(b);
			}
			
			else if(a.equals("taxtree") || a.equals("tree")){
				taxTreeFile=b;
			}
			
			else if(a.equals("name") || a.equals("taxname")){
				outTaxName=b;
			}else if(a.equals("name0")){
				outName0=b;
			}else if(a.equals("fname")){
				outFname=b;
			}else if(a.equals("outsketch") || a.equals("sketchout") || a.equals("outs") || a.equals("sketch")){
				outSketch=b;
			}else if(a.equals("files")){
				sketchFiles=Integer.parseInt(b);
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
			
			else if(searcher.parse(arg, a, b, false)){
//				System.err.println("*"+arg);
				parser.parse(arg, a, b); //Catches shared flags like "threads"
				Blacklist.parseBlacklist(arg, a, b); //Catches flags like "nt" or "refseq"
			}
			
			else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}
			
			else if(searcher.parse(arg, a, b, true)){
//				System.err.println("**"+arg);
				//do nothing
			}
			
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		
		outMeta=SketchObject.fixMeta(outMeta);
		SketchObject.postParse();
		
		if(skipCompare){
			allToAll=false;
			searcher.autoIndex=false;
			makeIndex=false;
			in.addAll(searcher.refFiles);
			searcher.refFiles.clear();
		}else if(in.isEmpty() && args.length>0 && !allToAll){ //Allows first argument to be used as the input file without in= flag
			String x=args[0];
			if(x.indexOf('=')<0 && new File(x).exists() && searcher.refFiles.contains(x)){
				searcher.refFiles.remove(x);
				in.add(x);
			}
		}
		
		{//Process parser fields
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			out=parser.out1;
		}
		
//		assert(false) : in+"\n"+searcher.refFiles;
		
		if(allToAll){
			 LinkedHashSet<String> set=new LinkedHashSet<String>();
			 set.addAll(in);
			 set.addAll(searcher.refFiles);
			 in.clear();
			 searcher.refFiles.clear();
			 in.addAll(set);
			 searcher.refFiles.addAll(set);
		}
		
		//Ensure there is an input file
		if(in.isEmpty() && !skipCompare){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Ensure there is an ref file
		if(searcher.refFiles.isEmpty() && !skipCompare){
			if(outSketch==null){throw new RuntimeException("Error - at least one reference file is required.");}
		}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		ffout=FileFormat.testOutput(out, FileFormat.TEXT, null, false, overwrite, append, ordered);
		if(!ffout.stdio() && !defaultParams.setColors){defaultParams.printColors=false;}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, taxTreeFile)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		if(!Tools.testInputFiles(true, false, in.toArray(new String[0]))){
			if(in.size()==1){
				String s=in.get(0);
				String s1=s.replaceFirst("#", "1"), s2=s.replaceFirst("#", "2");
				Tools.testInputFiles(true, false, s1, s2);
			}else{
				throw new RuntimeException("\nCan't read some input files.\n");  
			}
		}
		
//		assert(makeIndex || defaultParams.printContam2==false) : "Contam2 requires the flag index=t";
		
		SSUMap.load(outstream);
		if(taxTreeFile!=null){setTaxtree(taxTreeFile, silent ? null : outstream);}
		defaultParams.postParse(true, true);
		if(!defaultParams.printSSU){processSSU=false;}
		allowMultithreadedFastq=in.size()<2 && !allToAll;
		if(!allowMultithreadedFastq){Shared.capBufferLen(40);}
//		assert(defaultParams.checkValid());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Main processing method that loads sketches and performs comparisons.
	 * Handles different comparison modes (all-to-all, whitelist, standard).
	 * Uses multithreading for large sketch sets and outputs results.
	 * @param t Timer for tracking execution time
	 */
	public void process(Timer t){
		Timer ttotal=new Timer();
		
		t.start();
		
		if(!silent){outstream.println("Loading sketches.");}
		searcher.makeTool(1, false, defaultParams.mergePairs);
		SketchTool tool=new SketchTool(targetSketchSize, defaultParams);
		
		final int mode2=(defaultParams.mode==PER_FILE ? PER_FILE : PER_TAXA);
		if(skipCompare){
			makeIndex=false;
			inSketches=tool.loadSketches_MT(defaultParams, in);
		}else if(!useWhitelist || allToAll){
			if(allToAll){
				makeIndex=searcher.refFileCount()>0 && (makeIndex || defaultParams.needIndex() || searcher.autoIndex);
				searcher.loadReferences(mode2, defaultParams);
				inSketches=(ArrayList<Sketch>) searcher.refSketches.clone();
			}else{
				inSketches=tool.loadSketches_MT(defaultParams, in);
				
				for(Sketch sk : inSketches){
					if(sk.taxID<1 || sk.taxID>=minFakeID || outTaxID>0){sk.taxID=outTaxID;}
					if(outSpid>0){sk.spid=outSpid;}
					if(outImgID>0){sk.imgID=outImgID;}
					if(outTaxName!=null){sk.setTaxName(outTaxName);}
					if(outFname!=null){sk.setFname(outFname);}
					if(outName0!=null){sk.setName0(outName0);}
					if(SketchMaker.parseSubunit && sk.name0()!=null){
						if(outMeta!=null){
							sk.meta=(ArrayList<String>)sk.meta.clone();
						}else if(sk.meta==null){
							if(sk.name0().contains("SSU_")){
								sk.addMeta("subunit:ssu");
							}else if(sk.name0().contains("LSU_")){
								sk.addMeta("subunit:lsu");
							}
						}
					}
					sk.setMeta(outMeta);
					if(defaultParams.printSSU()){sk.loadSSU();}//since taxID was just set
				}
				
				if(outTaxID>0){
					for(Sketch sk : inSketches){
						if(sk.taxID<1 || sk.taxID>=minFakeID){sk.taxID=outTaxID;}
					}
				}
				makeIndex=searcher.refFileCount()>0 && ((searcher.autoIndex && inSketches.size()>8) || defaultParams.needIndex() || (makeIndex && !searcher.autoIndex));
				searcher.loadReferences(mode2, defaultParams);
				if(mode2==PER_FILE){
					int max=inSketches.size();
					for(int i=0; i<searcher.refSketches.size(); i++){
						searcher.refSketches.get(i).sketchID=max+i+1;
					}
				}
			}
		}else{
			//assert(searcher.makeIndex && !searcher.autoIndex) : "whitelist=t requires index=t";
			makeIndex=true; //(searcher.refFileCount()>0); //Index is required in whitelist mode.
			searcher.loadReferences(mode2, defaultParams);
			inSketches=tool.loadSketches_MT(defaultParams, in);
		}
		
		if(outSketch!=null){
			writeSketches(outSketch, sketchFiles);
		}
		
		final int numLoaded=(inSketches.size()+searcher.refSketchCount())/(allToAll ? 2 : 1);
		t.stop();
		if(!silent){outstream.println("Loaded "+numLoaded+" sketch"+(numLoaded==1 ? "" : "es")+" in "+t.toString());}
		if(printMemory){
			System.gc();
			Shared.printMemory();
		}
		
		if(skipCompare) {
			ttotal.stop("Total Time: \t");
			return;
		}
		
		t.start();

		
		ByteStreamWriter tsw=(ffout==null ? null : new ByteStreamWriter(ffout));
		if(tsw!=null){
			tsw.start();
			if(defaultParams.format==DisplayParams.FORMAT_QUERY_REF_ANI || defaultParams.format==DisplayParams.FORMAT_CONSTELLATION){
				String s=defaultParams.header()+"\n";
				tsw.forcePrint(s.getBytes());
			}
		}

		boolean success=true;
		final int inSize=inSketches.size();
		if(inSize==1 || Shared.threads()<2 || inSize<4){
			ByteBuilder sb=new ByteBuilder();
			success=searcher.compare(inSketches, sb, defaultParams, Shared.threads());
			success&=(!searcher.errorState);
			if(tsw!=null){
				sb.append('\n');
				if(ordered){
					tsw.addJob(sb);
				}else{
					tsw.println(sb);
				}
			}
		}else{//More sketches than threads, and more than one thread
			final int threads=Tools.min(Shared.threads(), inSize);
			
			ArrayList<CompareThread> alct=new ArrayList<CompareThread>(threads);
			AtomicInteger next=new AtomicInteger(0);
			for(int i=0; i<threads; i++){
				alct.add(new CompareThread(i, next, tsw));
			}
			for(CompareThread ct : alct){ct.start();}
			for(CompareThread ct : alct){

				//Wait until this thread has terminated
				while(ct.getState()!=Thread.State.TERMINATED){
					try {
						//Attempt a join operation
						ct.join();
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}

				synchronized(ct){
					success&=ct.success;
				}
			}
			alct=null;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		if(tsw!=null){errorState|=tsw.poisonAndWait();}
		
		t.stop();
//		long comparisons=(makeIndex ? searcher.comparisons.get() : 
//			allToAll ? (inSketches.size()*(long)(inSketches.size()-(compareSelf ? 0 : 1)))
//					: inSketches.size()*(long)searcher.refSketchCount());
		long comparisons=searcher.comparisons.get();
		if(!skipCompare && !silent) {outstream.println("\nRan "+comparisons+" comparison"+(comparisons==1 ? "" : "s")+" in "+t);}
		ttotal.stop();
		if(!silent){outstream.println("Total Time: \t"+ttotal);}
	}
	
	/**
	 * Writes loaded sketches to output files.
	 * Delegates to single file or multiple file writing based on filename pattern.
	 * @param fname Output filename pattern (may contain # for multiple files)
	 * @param files Number of files to split output across
	 */
	void writeSketches(String fname, int files){
		if(fname==null){return;}
		if(files==1 || fname.indexOf('#')<0){
			writeOneSketchFile(fname);
		}else{
			writeManySketchFiles(fname, files);
		}
	}
	
	/** Writes all sketches to a single output file.
	 * @param fname Output filename */
	void writeOneSketchFile(String fname){
		if(fname==null){return;}
		ByteBuilder bb=new ByteBuilder();
		ByteStreamWriter bsw=new ByteStreamWriter(outSketch, overwrite, append, true, FileFormat.SKETCH);
		bsw.start();
		for(Sketch sk : inSketches){
			sk.toBytes(bb);
			bsw.print(bb);
			bb.clear();
		}
		bsw.poisonAndWait();
		errorState|=bsw.errorState;
	}
	
	/**
	 * Writes sketches to multiple output files based on sketch ID.
	 * Uses modulo distribution to assign sketches to files.
	 * @param fname Filename pattern with # placeholder
	 * @param files Number of output files to create
	 */
	void writeManySketchFiles(String fname, int files){
		if(fname==null){return;}
		assert(fname.indexOf('#')>=0) : fname;
		assert(files>0) : files;
		
		ByteStreamWriter[] bswa=new ByteStreamWriter[files];
		for(int i=0; i<files; i++){
			ByteStreamWriter bsw=new ByteStreamWriter(outSketch.replaceFirst("#", ""+i), overwrite, append, true, FileFormat.SKETCH);
			bsw.start();
			bswa[i]=bsw;
		}
		for(Sketch sk : inSketches){
			ByteBuilder bb=new ByteBuilder(4096);
			sk.toBytes(bb);
			bswa[sk.sketchID%files].addJob(bb);
		}
		for(ByteStreamWriter bsw : bswa){
			bsw.poisonAndWait();
			errorState|=bsw.errorState;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Adds file paths to a collection, handling comma-separated lists.
	 * If the input contains commas and doesn't represent an existing file,
	 * splits on commas and adds each part separately.
	 *
	 * @param a File path or comma-separated file paths
	 * @param list Collection to add file paths to
	 */
	private static void addFiles(String a, Collection<String> list){
		if(a==null){return;}
		File f=null;
		if(a.indexOf(',')>=0){f=new File(a);}
		if(f==null || f.exists()){
			list.add(a);
		}else{
			for(String s : a.split(",")){list.add(s);}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Worker thread for parallel sketch comparison.
	 * Each thread processes a subset of input sketches against the reference set. */
	private class CompareThread extends Thread {
		
		/**
		 * Creates a comparison thread.
		 * @param tid_ Thread ID for this worker
		 * @param nextSketch_ Shared counter for sketch assignment
		 * @param tsw_ Writer for output results
		 */
		CompareThread(final int tid_, final AtomicInteger nextSketch_, ByteStreamWriter tsw_){
			tid=tid_;
			nextSketch=nextSketch_;
			tsw=tsw_;
		}
		
		@Override
		public void run(){
			success=false;
			final int inLim=inSketches.size();
			final boolean json=defaultParams.json();
			
			for(int inNum=nextSketch.getAndIncrement(); inNum<inLim; inNum=nextSketch.getAndIncrement()){
				Sketch a=inSketches.get(inNum);
				assert(buffer.cbs==null); //Because this sketch will only be used by one thread at a time, so per-buffer bitsets are not needed.
				SketchResults sr=searcher.processSketch(a, buffer, fakeID, map, defaultParams, 1);
				a.clearRefHitCounts();
				
				if(tsw!=null){
					ByteBuilder sb=sr.toText(defaultParams);
					synchronized(tsw){
						if(ordered){
							if(json){
								if(inNum==0){
									sb.insert(0, (byte)'[');//Rare, slow case
								}
								if(inNum<inLim-1){
									sb.append(',');
								}else{
									sb.append(']');
								}
							}
							tsw.add(sb, inNum);
						}else{
							if(json){
								if(resultsPrinted==0){
									tsw.print('[');
								}else{
									sb.insert(0, (byte)',');
								}
							}
							tsw.print(sb);
						}
						resultsPrinted++;
					}
				}
			}
			synchronized(this){success=true;}
		}
		
		/** Thread ID for this comparison worker */
		private final int tid;
		/** Reusable buffer for sketch comparisons */
		private final CompareBuffer buffer=new CompareBuffer(false);

		/** Shared counter for thread-safe sketch assignment */
		private final AtomicInteger nextSketch;
		/** Counter for generating fake taxonomic IDs */
		private final AtomicInteger fakeID=new AtomicInteger(minFakeID);
		/** Thread-safe map for storing comparison results */
		private ConcurrentHashMap<Integer, Comparison> map=new ConcurrentHashMap<Integer, Comparison>(101);
		/** Writer for outputting comparison results */
		final ByteStreamWriter tsw;
		
		/** Whether this thread completed successfully */
		boolean success=false;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Input sketch file paths */
	private ArrayList<String> in=new ArrayList<String>();
	
	/** Primary output file path */
	private String out="stdout.txt";
	
	/** Taxonomic tree file path for taxonomic analysis */
	private String taxTreeFile=null;
	
	/** Loaded input sketches for comparison */
	private ArrayList<Sketch> inSketches;
	
	/** Handles sketch searching and reference loading */
	public final SketchSearcher searcher=new SketchSearcher();
	
	/** Whether to print memory usage statistics */
	private boolean printMemory=false;
	/** Whether to suppress status messages */
	private boolean silent=false;
	
	/*Override metadata */
	/** Override taxonomic name for output sketches */
	private String outTaxName=null;
	/** Override filename for output sketches */
	private String outFname=null;
	/** Override primary name for output sketches */
	private String outName0=null;
	/** Output file path for writing sketches */
	private String outSketch=null;
	/** Number of files to split sketch output across */
	private int sketchFiles=1;
	/** Override taxonomic ID for output sketches */
	private int outTaxID=-1;
	/** Override species ID for output sketches */
	private long outSpid=-1;
	/** Override IMG ID for output sketches */
	private long outImgID=-1;
	/** Additional metadata to add to output sketches */
	private ArrayList<String> outMeta=null;
	/** Counter for number of results printed (for JSON formatting) */
	private long resultsPrinted=0;
	
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
	/** Whether to maintain output order for multithreaded processing */
	private boolean ordered=true;
	
}
