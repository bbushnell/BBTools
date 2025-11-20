package bbdukGemini;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLongArray;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import kmer.AbstractKmerTable;
import kmer.ScheduleMaker;
import shared.KillSwitch;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * Main driver for BBDuk. Orchestrates initialization, k-mer loading, 
 * worker processing, and statistics aggregation, using the multi-threaded 
 * Streamer/Writer/Accumulator pattern.
 * @author Brian Bushnell
 * @contributor Gemini
 * @date November 19, 2025
 */
public class BBDuk2 implements Accumulator<BBDuk2.ProcessThread>{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void main(String[] args){
		BBDuk2 x=new BBDuk2(args);
		x.process();
		Shared.closeStream(outstream);
	}
	
	public BBDuk2(String[] args){
		PreParser pp=new PreParser(args, getClass(), true);
		args=pp.args;
		outstream=pp.outstream;
		
		BBDukParser bbdp=new BBDukParser();
		Parser parser=new Parser();
		
		// Set global defaults
		ReadWrite.ZIPLEVEL=2;
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_PIGZ=true;
		if(!fileIO.ByteFile.FORCE_MODE_BF1 && !fileIO.ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			fileIO.ByteFile.FORCE_MODE_BF2=true;
		}
		
		// Parse arguments
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseZip(arg, a, b)){}
			else if(Parser.parseHist(arg, a, b)){}
			else if(Parser.parseCommonStatic(arg, a, b)){}
			else if(Parser.parseQualityAdjust(arg, a, b)){}
			else if(Parser.parseQuality(arg, a, b)){}
			else if(Parser.parseFasta(arg, a, b)){}
			else if(Parser.parseSam(arg, a, b)){}
			else if(parser.parseInterleaved(arg, a, b)){}
			else if(parser.parseTrim(arg, a, b)){}
			else if(parser.parseCommon(arg, a, b)){}
			else if(parser.parseCardinality(arg, a, b)){}
			else if(parser.parseFiles(arg, a, b)){} 
			else if(bbdp.parse(arg, a, b, i)){}      
			else throw new RuntimeException("Unknown parameter "+args[i]);
		}
		
		bbdp.postParse(Shared.AMINO_IN);
		parser.processQuality();
		
		// Handle Pound replacement for matched/bad output
		if(bbdp.outb1!=null && bbdp.outb2==null && bbdp.outb1.contains("#")){
			bbdp.outb2=bbdp.outb1.replace("#", "2");
			bbdp.outb1=bbdp.outb1.replace("#", "1");
		}
		
		// Capture BBDuk fields needed for harness logic
		altref=bbdp.altref;
		ktrimLeft=bbdp.ktrimLeft;
		ktrimRight=bbdp.ktrimRight;
		ktrimN=bbdp.ktrimN;
		ksplit=bbdp.ksplit;
		
		// Capture Parser fields needed for harness logic
		ordered=bbdp.ordered;
		maxReads=parser.maxReads;
		sampleseed=parser.sampleseed;
		samplerate=parser.samplerate;
		
		// Stat configuration
		makeReadStats=ReadStats.collectingStats();
		hitCounts=null; 

		// Initialize global statistics arrays
		int arraySize=initialSizeDefault; 
		
		refScafCounts=new int[bbdp.ref==null ? 1 : bbdp.ref.length];
		scaffoldReadCounts=new AtomicLongArray(arraySize);
		scaffoldBaseCounts=new AtomicLongArray(arraySize);
		
		THREADS=Shared.threads();
		
		// Initialize Loader
		loader=new BBDukIndexLoader(bbdp, parser, refScafCounts, REPLICATE_AMBIGUOUS);

		// Set processor defaults
		processor=new BBDukProcessor(bbdp, parser, null, scaffoldReadCounts, scaffoldBaseCounts);
		
		// I/O FileFormats
		ffin1=FileFormat.testInput(parser.in1, FileFormat.FASTQ, null, true, true);
		ffin2=FileFormat.testInput(parser.in2, FileFormat.FASTQ, null, true, true);
		
		// Output for "Clean" (Unmatched) reads
		ffout1=FileFormat.testOutput(parser.out1, FileFormat.FASTQ, null, true, parser.overwrite, parser.append, ordered);
		ffout2=FileFormat.testOutput(parser.out2, FileFormat.FASTQ, null, true, parser.overwrite, parser.append, ordered);

		// Output for "Dirty" (Matched) reads
		ffoutb1=FileFormat.testOutput(bbdp.outb1, FileFormat.FASTQ, null, true, parser.overwrite, parser.append, ordered);
		ffoutb2=FileFormat.testOutput(bbdp.outb2, FileFormat.FASTQ, null, true, parser.overwrite, parser.append, ordered);
		
		// Output for Singletons
		ffouts=FileFormat.testOutput(parser.outsingle, FileFormat.FASTQ, null, true, parser.overwrite, parser.append, ordered);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	public void process(){
		Timer t=new Timer();
		
		// 1. Initialize Kmer Index Tables
		ScheduleMaker scheduleMaker=new ScheduleMaker(WAYS, 12, prealloc_, preallocFraction);
		int[] schedule=scheduleMaker.makeSchedule();
		outstream.print("Allocating kmer table:\t");
		AbstractKmerTable[] keySets=AbstractKmerTable.preallocate(WAYS, AbstractKmerTable.ARRAY1D, schedule, -1L);
		
		// 2. Load References
		loader.setKeySets(keySets); 
		long storedKmers=loader.spawnLoadThreads(false);
		BBDukIndex index=loader.getIndex(); 
		
		if(storedKmers<1 && altref!=null && altref.length>0){
			outstream.println("Ref had no kmers; using alt ref.");
			storedKmers=loader.spawnLoadThreads(true);
		}
		
		if(storedKmers<1 && (ktrimLeft || ktrimRight || ktrimN || ksplit)){
			KillSwitch.kill("Kmer operation chosen but no kmers loaded.");
		}
		
		// 3. Update Processor with loaded Index
		BBDukProcessor sharedProcessor=processor.copyWithIndex(index);
		
		// 4. Setup Streamer and Writers
		final boolean saveHeader=false; 
		
		Streamer streamer=StreamerFactory.makeStreamer(ffin1, ffin2, ordered, maxReads, saveHeader, true);
		streamer.setSampleRate(samplerate, sampleseed);
		
		Writer writerUnmatched=WriterFactory.makeWriter(ffout1, ffout2, -1, null, saveHeader);
		Writer writerMatched=WriterFactory.makeWriter(ffoutb1, ffoutb2, -1, null, saveHeader);
		Writer writerSingles=WriterFactory.makeWriter(ffouts, null, -1, null, saveHeader);
		
		// Start streams
		streamer.start();
		if(writerUnmatched!=null){writerUnmatched.start();}
		if(writerMatched!=null){writerMatched.start();}
		if(writerSingles!=null){writerSingles.start();}
		
		// 5. Spawn Worker Threads
		spawnProcessThreads(streamer, writerUnmatched, writerMatched, writerSingles, sharedProcessor);
		
		// 6. Cleanup and Stats
		if(writerUnmatched!=null){writerUnmatched.poisonAndWait();}
		if(writerMatched!=null){writerMatched.poisonAndWait();}
		if(writerSingles!=null){writerSingles.poisonAndWait();}
		
		errorState|=ReadStats.writeAll();
		errorState|=ReadWrite.closeStreams(streamer, writerUnmatched, writerMatched, writerSingles);
		
		t.stop();
		
		// Aggregate final stats
		readsIn=sharedProcessor.readsInT;
		basesIn=sharedProcessor.basesInT;
		
		outstream.println(Tools.timeReadsBasesProcessed(t, readsIn, basesIn, 8));
	}
	
	private void spawnProcessThreads(Streamer streamer, Writer writerUnmatched, Writer writerMatched, Writer writerSingles, BBDukProcessor sharedProcessor){
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(THREADS);
		for(int i=0; i<THREADS; i++){
			alpt.add(new ProcessThread(streamer, writerUnmatched, writerMatched, writerSingles, sharedProcessor.clone(), i));
		}
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState|=!success;
	}
	
	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt){
			errorState|=!pt.success;
			processor.add(pt.processorT); 
		}
	}
	
	@Override
	public final ReadWriteLock rwlock(){return rwlock;}

	@Override
	public final boolean success(){return !errorState;}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	class ProcessThread extends Thread{

		ProcessThread(Streamer streamer_, Writer writerUnmatched_, Writer writerMatched_, Writer writerSingles_, BBDukProcessor processorT_, int tid_){
			streamer=streamer_;
			writerUnmatched=writerUnmatched_;
			writerMatched=writerMatched_;
			writerSingles=writerSingles_;
			processorT=processorT_;
			tid=tid_;
		}

		@Override
		public void run(){
			processInner();
			success=true;
		}

		void processInner(){
			for(ListNum<Read> ln=streamer.nextList(); ln!=null; ln=streamer.nextList()){
				ArrayList<Read> reads=ln.list;
				ArrayList<Read> listUnmatched=(reads!=null ? new ArrayList<Read>(reads.size()) : null); 
				ArrayList<Read> listMatched=(writerMatched==null || reads==null ? null : new ArrayList<Read>(reads.size()));
				ArrayList<Read> listSingles=(writerSingles==null || reads==null ? null : new ArrayList<Read>(reads.size()));
				
				if(reads!=null && !reads.isEmpty()){
					for(int i=0; i<reads.size(); i++){
						Read r1=reads.get(i);
						Read r2=r1.mate;

						int keep=processorT.processReadPair(r1, r2);

						if(keep==3){
							listUnmatched.add(r1);
						}else if(keep==1){
							if(r2!=null){r1.mate=null; r2.mate=null;}
							if(listSingles!=null){listSingles.add(r1);}
							else{listUnmatched.add(r1);}

							if(listMatched!=null && r2!=null){listMatched.add(r2);}
						}else if(keep==2){
							if(r1!=null){r1.mate=null; r2.mate=null;}
							if(listSingles!=null){listSingles.add(r2);}
							else{listUnmatched.add(r2);}

							if(listMatched!=null && r1!=null){listMatched.add(r1);}
						}else{
							if(listMatched!=null){
								if(r1!=null){listMatched.add(r1);}
							}
						}
					}
				}
				
				// Always send list to writers to maintain order, even if empty
				if(writerUnmatched!=null){
					writerUnmatched.addReads(new ListNum<Read>(listUnmatched, ln.id));
				}
				if(writerMatched!=null){
					writerMatched.addReads(new ListNum<Read>(listMatched, ln.id));
				}
				if(writerSingles!=null){
					writerSingles.addReads(new ListNum<Read>(listSingles, ln.id));
				}
			}
		}

		boolean success=false;
		final int tid;

		private final Streamer streamer;
		private final Writer writerUnmatched;
		private final Writer writerMatched;
		private final Writer writerSingles;
		final BBDukProcessor processorT; 
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/

	private static PrintStream outstream=System.err;
	public static final int WAYS=7;
	private static final int initialSizeDefault=128000;
	public static int THREADS=Shared.threads();
	public static boolean REPLICATE_AMBIGUOUS=false;
	public static boolean prealloc_=false;
	public static double preallocFraction=0.9;
	public static AtomicLongArray scaffoldReadCounts;
	public static AtomicLongArray scaffoldBaseCounts;
	
	// Static config for processor
	public static boolean makeReadStats;
	public static long[] hitCounts;
	public static final int HITCOUNT_LEN=1000;
	
	/*--------------------------------------------------------------*/
	/*----------------         Local Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final BBDukProcessor processor;
	private final BBDukIndexLoader loader; 
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

	private final int[] refScafCounts;
	private final FileFormat ffin1, ffin2;
	private final FileFormat ffout1, ffout2;
	private final FileFormat ffoutb1, ffoutb2;
	private final FileFormat ffouts;
	
	// Fields captured from Parsers
	private final String[] altref;
	private final boolean ktrimLeft, ktrimRight, ktrimN, ksplit;
	private final boolean ordered;
	private final long maxReads;
	private final long sampleseed;
	private final float samplerate;
	
	long readsIn=0;
	long basesIn=0;
	boolean errorState=false;
}