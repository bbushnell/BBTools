package sketch;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.concurrent.atomic.AtomicInteger;

import bloom.KmerCountAbstract;
import bloom.ReadCounter;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import kmer.HashArrayHybridFast;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FastaReadInputStream;
import structures.IntList;
import structures.IntListCompressor;
import structures.LongList;
import tax.TaxNode;
import tax.TaxTree;

/**
 * Makes blacklists from existing sketches.
 * 
 * @author Brian Bushnell
 * @date November 12, 2019
 *
 */
public class BlacklistMaker2 extends SketchObject {
	
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
		
		ReadCounter.maxShortKmerLength=32;
		
		//Create an instance of this class
		BlacklistMaker2 x=new BlacklistMaker2(args);
		
		//Run the object
		x.process(t);
		
		ReadCounter.maxShortKmerLength=31;
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	@SuppressWarnings("unchecked")
	public BlacklistMaker2(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		if(Shared.threads()>=32){Shared.capThreads(Shared.threads()/2);}//Deals with hyperthreading, which causes excessive memory usage
		Shared.capBuffers(Tools.max(4, Shared.threads()+1));
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		KmerCountAbstract.SKETCH_MODE=true;
		KmerCountAbstract.STORE_HASHED=true;
		KmerCountAbstract.KEEP_DUPLICATE_KMERS=true;
		setKeyFraction(1); //Prevents an assertion error by ignoring min hash value
		
		//Create a parser object
		Parser parser=new Parser();
		
		int mode_=PER_TAXA;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parseMode(arg, a, b)>-1){
				mode_=parseMode(arg, a, b);
			}else if(a.equals("tossjunk")){
				tossJunk=Parse.parseBoolean(b);
			}
			
			else if(a.equals("taxtree") || a.equals("tree")){
				taxTreeFile=b;
			}else if(a.equals("imgfile") || a.equals("imgdump")){
				imgFile=b;
			}
			
			
			else if(a.equals("mincount") || a.equals("mintaxcount")){
				minTaxCount=Parse.parseIntKMG(b);
			}else if(a.equals("maxkeys") || a.equals("keys") || a.equals("length") || a.equals("size")){
				maxKeys=Parse.parseIntKMG(b);
			}else if(a.equals("name")){
				outName=b;
			}else if(a.equalsIgnoreCase("name0") || a.equalsIgnoreCase("nm0")){
				sketchName=b;
			}else if(a.equals("hist")){
				outHist=b;
			}
			
			else if(a.equals("taxlevel") || a.equals("tl") || a.equals("level") || a.equals("lv")){
				if(b==null){taxLevel=-1;}
				else if(Tools.isDigit(b.charAt(0))){
					taxLevel=Integer.parseInt(b);
				}else{
					taxLevel=TaxTree.parseLevel(b);
				}
			}
			

			
			else if(searcher.parse(arg, a, b, false)){
//				System.err.println("*"+arg);
				parser.parse(arg, a, b); //Catches shared flags like "threads"
			}else if(parseSketchFlags(arg, a, b)){
				//do nothing
			}else if(defaultParams.parse(arg, a, b)){
				//do nothing
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(searcher.parse(arg, a, b, true)){
//				System.err.println("**"+arg);
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		if("auto".equalsIgnoreCase(imgFile)){imgFile=TaxTree.defaultImgFile();}
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		
		mode=mode_;
		assert((mode!=PER_TAXA && mode!=PER_IMG) || taxTreeFile!=null) : "Please specify a TaxTree.";
		assert(mode!=PER_IMG || imgFile!=null);
		assert(mode==PER_TAXA || mode==PER_SEQUENCE || mode==PER_IMG);
		
		{//Process parser fields
			Parser.processQuality();
			
			overwrite=parser.overwrite;
			append=parser.append;

			outSketch=parser.out1;
		}
		
		postParse();
		
		assert(FastaReadInputStream.settingsOK());
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, outSketch)){
			outstream.println((outSketch==null)+", "+outSketch);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+outSketch+"\n");
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, outSketch, outHist)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}

		if(sketchName==null && outSketch!=null){
			sketchName=outSketch;
			sketchName=ReadWrite.stripToCore(sketchName);
		}
		
		//Create output FileFormat objects
		ffsketch=FileFormat.testOutput(outSketch, FileFormat.SKETCH, null, true, overwrite, append, false);
		ffhist=FileFormat.testOutput(outHist, FileFormat.TXT, null, true, overwrite, append, false);
		
		if(taxTreeFile!=null){setTaxtree(taxTreeFile, outstream);}
		
		if(imgFile!=null){
			TaxTree.loadIMG(imgFile, false, outstream);
		}
		
		maps=new HashMap[ways];
		for(int i=0; i<ways; i++){
			maps[i]=new HashMap<Long, IntListCompressor>();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		makeIndex=true; //(searcher.refFileCount()>0); //Index is required in whitelist mode.
		searcher.loadReferences(mode, defaultParams);
		sketchesProcessed=searcher.refSketchCount();
		
		//Reset counters
		keysProcessed=0;
		
		//Process the reads in separate threads
		spawnThreads();
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Report timing and results
		t.stop();
		outstream.println("Blacklist size: \t"+resultingSize+"\n");
		outstream.println(Tools.timeSketchesKeysProcessed(t, sketchesProcessed, keysProcessed, 8));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Spawn process threads */
	private void spawnThreads(){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(i, threads));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		//Wait for completion of all threads
		boolean success=true;
		for(ProcessThread pt : alpt){
			
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
			keysProcessed+=pt.keysProcessedT;
			success&=pt.success;
		}
		
		shrinkListsAndWriteHist();
		writeSketch(true);
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		//Do anything necessary after processing
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void writeSketch(boolean destroy){
		Sketch sk=toSketch(destroy);
		if(ffsketch!=null){errorState|=SketchTool.write(sk, ffsketch);}
	}
	
	private void shrinkListsAndWriteHist(){
		int max=1000000;
		long[] counts=new long[max+1];
		for(int i=0; i<ways; i++){
			for(Entry<Long, IntListCompressor> entry : maps[i].entrySet()){
				IntListCompressor value=entry.getValue();
				value.sortAndShrink();
				IntList list=value.list;
				int index=Tools.min(max, list.size);
				counts[index]++;
			}
		}
		if(ffhist!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(ffhist);
			bsw.start();
			bsw.print("#count\tkmers\n".getBytes());
			for(int i=0; i<counts.length; i++){
				long count=counts[i];
				if(count>0){
					bsw.print(i);
					bsw.print('\t');
					bsw.print(count);
					bsw.print('\n');
				}
			}
			errorState|=bsw.poisonAndWait();
		}
	}

	private Sketch toSketch(boolean destroy){
		long[] array=toArray(destroy);
		hashArrayToSketchArray(array);
		ArrayList<String> meta=new ArrayList<String>();
		meta.add("minTaxCount:"+minTaxCount);
		meta.add("taxLevel:"+taxLevel);
		Sketch sk=new Sketch(array, null, null, null, null, outTaxid, -1, -1, -1, -1, -1, outName, sketchName, ffsketch.simpleName(), meta);
		return sk;
	}
	
	private long[] toArray(boolean destroy){
		ArrayList<KeyPair> pairs=new ArrayList<KeyPair>();
		long entries=0;
		
		for(int i=0; i<ways; i++){
			for(Entry<Long, IntListCompressor> entry : maps[i].entrySet()){
				Long key=entry.getKey();
				IntList value=entry.getValue().list;
				entries++;
				if(value.size()>=minTaxCount){
					KeyPair p=new KeyPair(key.longValue(), value.size());
					pairs.add(p);
				}
			}
			if(destroy){maps[i]=null;}
		}
		System.err.println("Pruned "+entries+" high-count keys down to "+pairs.size()+" after tax level promotion.");
		if(pairs.isEmpty()){
			System.err.println("*** Warning - no keys retained! ***");
			assert(false);
		}
		LongList list=new LongList(Tools.mid(1, pairs.size(), maxKeys));
		Collections.sort(pairs);
		for(KeyPair p : pairs){
			list.add(p.key);
			if(list.size()>=maxKeys){break;}
		}
		list.sort();
		resultingSize=list.size();
		return list.toArray();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final int tid_, final int threads_){
			threadID=tid_;
			threads=threads_;
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Iterate through the reads */
		void processInner(){
			SketchIndex index=searcher.index;
			for(int i=threadID; i<index.tableArray.length; i+=threads) {
				HashArrayHybridFast table=(HashArrayHybridFast)index.tableArray[i];
				long[] array=table.array();
				for(long key : array){
					if(key>=0){
						processKey(key, table);
					}
				}
			}
		}
		
		void processKey(final long key0, HashArrayHybridFast table){
			keysProcessedT++;
			final int[] sketchIds=table.getValues(key0, singleton); //searcher.index.getSketchIdsMap(key, singleton);
			if(sketchIds==null){//Keys should be 63 bits so this should not happen
				assert(key0<0) : Long.toHexString(key0);
				return;
			}
			if(sketchIds.length<minTaxCount){return;}
			
			long key=Long.MAX_VALUE-key0;
//			System.err.println(key+", "+minHashValue);
//			System.err.print(key<minHashValue ? "L" : "G");
			for(int sid : sketchIds){
				if(sid<0){break;}
				final int trueID=sid-1;//Minimum id is 1, indicating sketch 0.
				if(mode==PER_SEQUENCE){
//					System.err.print("A");
					addToMap(key, trueID);
				}else{
//					System.err.print("B");
					Sketch sk=searcher.refSketches.get(trueID);
					int taxID=sk.taxID;
					taxID=taxtree.promote(taxID, taxLevel);
					
					boolean ok=true;
					if(tossJunk){
						if(taxID<=1){ok=false;}
						else{
							TaxNode tn=taxtree.getNode(taxID);
							if(tn==null || !tn.isRanked() || tn.pid==TaxTree.LIFE_ID){ok=false;}
						}
					}
					if(ok){
						if(taxID<0){
							taxID=nextUnknown.getAndIncrement();
						}
//						System.err.print("C");
						addToMap(key, taxID);
					}
				}
			}
		}
		
		void addToMap(long key0, int value){
			keysAddedT++;
			Long key=Long.valueOf(key0);
//			Long key=new Long(key0);
			HashMap<Long, IntListCompressor> map=maps[(int)(key%ways)];
			IntListCompressor lh=map.get(key);
			if(lh==null){
				synchronized(map){
					lh=map.get(key);
					if(lh==null){
						lh=new IntListCompressor();
						map.put(key, lh);
					}
				}
			}
			synchronized(lh){
				lh.add(value);
			}
		}
		
		/** Number of bases processed by this thread */
		protected long keysProcessedT=0;
		
		protected long keysAddedT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Thread ID */
		final int threadID;
		final int threads;
		
		private final int[] singleton=new int[1];
	}
	
	/*--------------------------------------------------------------*/
	
	private static class KeyPair implements Comparable<KeyPair> {

		KeyPair(long key_, int count_){
			key=key_;
			count=count_;
		}
		
		@Override
		public int compareTo(KeyPair o) {
			if(count!=o.count){return o.count-count;}
			return key>o.key ? -1 : key<o.key ? 1 : 0;
		}
		
		long key;
		int count;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final SketchSearcher searcher=new SketchSearcher();
	
	private final int mode;
	
	private String taxTreeFile=null;
	private String imgFile=null;

	private String outName="blacklist";
	private String sketchName=null;
	private int outTaxid=-1;
	
	private int taxLevel=1;
	private boolean tossJunk=true;
	private int minTaxCount=20;
	private int maxKeys=300000;
	
	private HashMap<Long, IntListCompressor>[] maps;
	
	final int ways=63;
	
	int resultingSize=-1;
	
	private final AtomicInteger nextUnknown=new AtomicInteger(SketchObject.minFakeID);
	
	/*--------------------------------------------------------------*/

	/** Primary output file path */
	private String outSketch=null;

	/** Histogram output file path */
	private String outHist=null;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long sketchesProcessed=0;
	/** Number of bases processed */
	protected long keysProcessed=0;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary output file */
	private final FileFormat ffsketch;
	/** Histogram output file */
	private final FileFormat ffhist;
	
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
	
}
