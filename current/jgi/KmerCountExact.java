package jgi;

import java.io.PrintStream;
import java.util.ArrayList;

import assemble.Shaver;
import assemble.Tadpole;
import bloom.KmerCountAbstract;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import kmer.AbstractKmerTable;
import kmer.AbstractKmerTableSet;
import kmer.DumpThread;
import kmer.HashArray1D;
import kmer.HashForest;
import kmer.KmerNode;
import kmer.KmerNode1D;
import kmer.KmerTableSet;
import kmer.Walker;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sketch.Sketch;
import sketch.SketchObject;
import sketch.SketchTool;
import stream.FastaReadInputStream;
import tracker.ReadStats;
import ukmer.AbstractKmerTableU;
import ukmer.HashArrayU1D;
import ukmer.HashForestU;
import ukmer.Kmer;
import ukmer.KmerNodeU;
import ukmer.KmerNodeU1D;
import ukmer.KmerTableSetU;
import ukmer.WalkerU;

/**
 * @author Brian Bushnell
 * @date Nov 22, 2013
 *
 */
public class KmerCountExact {
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer(), t2=new Timer();
		t.start();
		t2.start();
		
		//Create a new CountKmersExact instance
		KmerCountExact x=new KmerCountExact(args);
		t2.stop();
//		outstream.println("Initialization Time:      \t"+t2);
		
		///And run it
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public KmerCountExact(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		/* Set global defaults */
		ReadWrite.ZIPLEVEL=2;
		ReadWrite.USE_UNPIGZ=true;
		
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		/* Initialize local variables with defaults */
		boolean useForest_=false, useTable_=false, useArray_=true;
		Parser parser=new Parser();
		
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(Parser.parseQuality(arg, a, b)){
				//do nothing
			}else if(Parser.parseFasta(arg, a, b)){
				//do nothing
			}else if(parser.parseInterleaved(arg, a, b)){
				//do nothing
			}else if(parser.parseTrim(arg, a, b)){
				//do nothing
			}else if(a.equals("out") || a.equals("out1") || a.equals("outkmers") || a.equals("outk") || a.equals("dump")){
				outKmers=b;
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("intersection") || a.equals("refout") || a.equals("outref")){
				intersectionFile=b;
			}else if(a.equals("mincounttodump") || a.equals("mindump") || a.equals("mincount")){
				minToDump=Integer.parseInt(b);
			}else if(a.equals("maxcounttodump") || a.equals("maxdump") || a.equals("maxcount")){
				maxToDump=Integer.parseInt(b);
			}else if(a.equals("dumpthreads")){
				DumpThread.NUM_THREADS=Integer.parseInt(b);
			}else if(a.equals("hist") || a.equals("khist")){
				outHist=b;
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("forest")){
				useForest_=Parse.parseBoolean(b);
				if(useForest_){useTable_=useArray_=false;}
				assert(false) : a+" is deprecated.";
			}else if(a.equals("table")){
				useTable_=Parse.parseBoolean(b);
				if(useTable_){useForest_=useArray_=false;}
				assert(false) : a+" is deprecated.";
			}else if(a.equals("array")){
				useArray_=Parse.parseBoolean(b);
				if(useArray_){useTable_=useForest_=false;}
				assert(false) : a+" is deprecated.";
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=(b==null || b.equalsIgnoreCase("auto") ? Shared.threads() : Integer.parseInt(b));
			}else if(a.equals("verbose")){
				assert(false) : "Verbose flag is currently static final; must be recompiled to change.";
//				verbose=Parse.parseBoolean(b);
			}else if(a.equals("shave")){
				shave=Parse.parseBoolean(b);
			}else if(a.equals("rinse")){
				rinse=Parse.parseBoolean(b);
			}else if(a.equals("shavedepth")){
				shaveDepth=Integer.parseInt(b);
			}else if(a.equals("histcolumns")){
				histColumns=Integer.parseInt(b);
			}else if(a.equals("histmax") || a.equals("histlen") || a.equals("khistlen") || a.equals("histsize") || a.equals("khistsize")){
				histMax=Parse.parseIntKMG(b);
			}else if(a.equals("refmax")){
				refMax=Parse.parseIntKMG(b);
			}else if(a.equals("histheader")){
				histHeader=Parse.parseBoolean(b);
			}else if(a.equals("nzo") || a.equals("nonzeroonly")){
				histZeros=!Parse.parseBoolean(b);
			}else if(a.equals("gchist")){
				gcHist=Parse.parseBoolean(b);
			}else if(a.equals("logscale")){
				doLogScale=Parse.parseBoolean(b);
			}else if(a.equals("logwidth")){
				logWidth=Double.parseDouble(b);
			}else if(a.equals("logpasses")){
				logPasses=Integer.parseInt(b);
			}
			
			else if(a.equals("minheight")){
				minHeight=Long.parseLong(b);
			}else if(a.equals("minvolume")){
				minVolume=Long.parseLong(b);
			}else if(a.equals("minwidth")){
				minWidth=Integer.parseInt(b);
			}else if(a.equals("minpeak")){
				minPeak=Integer.parseInt(b);
			}else if(a.equals("maxpeak")){
				maxPeak=Integer.parseInt(b);
			}else if(a.equals("maxpeakcount") || a.equals("maxpc") || a.equals("maxpeaks")){
				maxPeakCount=Integer.parseInt(b);
			}else if(a.equals("ploidy")){
				ploidy=Integer.parseInt(b);
			}else if(a.equals("peaks") || a.equals("peaksout")){
				outPeaks=b;
			}else if(a.equals("smooth") || a.equals("smoothe")){
				smoothKhist=smoothPeaks=Parse.parseBoolean(b);
			}else if(a.equals("smoothkhist") || a.equals("smoothhist")){
				smoothKhist=Parse.parseBoolean(b);
			}else if(a.equals("smoothpeaks")){
				smoothPeaks=Parse.parseBoolean(b);
			}else if(a.equals("smoothradius") || a.equals("smootheradius")){
				smoothRadius=Integer.parseInt(b);
			}else if(a.equals("maxradius")){
				CallPeaks.maxRadius=Integer.parseInt(b);
			}else if(a.equals("progressivemult")){
				CallPeaks.progressiveMult=Float.parseFloat(b);
			}
			
			else if(KmerTableSet.isValidArgument(a)){
				//Do nothing
			}else if(a.equals("decimals")){
				decimals=Integer.parseInt(b);
			}
			
			else if(a.equals("sketchmode")){
				KmerCountAbstract.SKETCH_MODE=Parse.parseBoolean(b);
			}else if(a.equals("sketch")){
				sketchPath=b;
			}else if(a.equals("sketchlen") || a.equals("sketchlength") || a.equals("sketchsize")){
				sketchLength=Parse.parseIntKMG(b);
			}else if(a.equals("sketchname")){
				sketchName=b;
			}else if(a.equals("sketchid")){
				sketchID=Integer.parseInt(b);
			}else if(SketchObject.parseSketchFlags(arg, a, b)){
				//Do nothing
			}
			
			else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		if(sketchPath!=null){
			SketchObject.postParse();
		}
		
		{//Process parser fields
			Parser.processQuality();
		}
		
		/* Adjust I/O settings and filenames */
		
		assert(FastaReadInputStream.settingsOK());

		if(outKmers!=null && !Tools.canWrite(outKmers, overwrite)){throw new RuntimeException("Output file "+outKmers+" already exists, and overwrite="+overwrite);}
		if(sketchPath!=null && !Tools.canWrite(sketchPath, overwrite)){throw new RuntimeException("Output file "+sketchPath+" already exists, and overwrite="+overwrite);}
		
		assert(THREADS>0);
		
		if(DISPLAY_PROGRESS){
			outstream.println("Initial:");
			Shared.printMemory();
			outstream.println();
		}
		
//		final int tableType=(useForest ? AbstractKmerTable.FOREST1D : useTable ? AbstractKmerTable.TABLE : useArray ? AbstractKmerTable.ARRAY1D : 0);
		k=Tadpole.preparseK(args);
		
		if(k<=31){//TODO: 123 add "false" to the clause to force KmerTableSetU usage.
			tables=new KmerTableSet(args, 12);
			if(ref!=null){tables2=new KmerTableSet(new String[] {"k="+k, "in="+ref, "rcomp="+tables.rcomp()}, 12);}
		}else{
			tables=new KmerTableSetU(args, 0);
			if(ref!=null){tables2=new KmerTableSetU(new String[] {"k="+k, "in="+ref, "rcomp="+tables.rcomp()}, 0);}
		}
		if(tables.prefilter){tables.minProbMain=false;}
		
		ffSketch=FileFormat.testOutput(sketchPath, FileFormat.TXT, null, true, overwrite, append, false);
		

//		shift=bitsPerBase*k;
//		shift2=shift-bitsPerBase;
//		mask=(shift>63 ? -1L : ~((-1L)<<shift));
//		kmask=lengthMasks[k];
	}

	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Main processing method that executes the complete k-mer counting workflow.
	 * Counts k-mers, generates histograms, calls peaks, and creates sketches as configured.
	 * @param t Timer for tracking execution time across all phases
	 */
	public void process(Timer t){
		
		/* Check for output file collisions */
		Tools.testOutputFiles(overwrite, append, false, outKmers, outHist, outPeaks, sketchPath);
		
		/* Count kmers */
		process2();
		
		makeKhistAndPeaks();
		
		if(ffSketch!=null){
			makeSketch();
		}
		
		/* Stop timer and calculate speed statistics */
		t.stop();
		
		/* Throw an exception if errors were detected */
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	/**
	 * Core k-mer counting phase that fills hash tables with k-mer data.
	 * Loads input sequences, extracts k-mers, applies optional shaving/rinsing,
	 * and computes statistics on k-mer distribution.
	 */
	public void process2(){
		
		/* Start phase timer */
		Timer t=new Timer();
		
		AbstractKmerTableSet.DISPLAY_STATS=false;
		
		/* Fill tables with kmers */
		tables.process(t);
		
		if(DISPLAY_PROGRESS){
			outstream.println("After loading:");
			Shared.printMemory();
			outstream.println();
		}
		
		errorState|=tables.errorState;
		
		t.stop();
		outstream.println("Input:                      \t"+tables.readsIn+" reads \t\t"+tables.basesIn+" bases.");
		
		if(ref!=null){
			t.start();
			processRef();
			t.stop();
//			assert(false) : "Display stuff";
		}
		
		if(tables.qtrimLeft() || tables.qtrimRight()){
			outstream.println("QTrimmed:               \t"+tables.readsTrimmed+" reads ("+Tools.format("%.2f",tables.readsTrimmed*100.0/tables.readsIn)+"%) \t"+
					tables.basesTrimmed+" bases ("+Tools.format("%.2f",tables.basesTrimmed*100.0/tables.basesIn)+"%)");
		}
		if(tables.minAvgQuality()>0){
			outstream.println("Low quality discards:   \t"+tables.lowqReads+" reads ("+Tools.format("%.2f",tables.lowqReads*100.0/tables.readsIn)+"%) \t"+
					tables.lowqBases+" bases ("+Tools.format("%.2f",tables.lowqBases*100.0/tables.basesIn)+"%)");
		}
		
		if(shave || rinse){
			kmersRemoved=shave(shave, rinse, shaveDepth);
		}
		
		outstream.println("\nFor K="+tables.kbig());
		outstream.println("Unique Kmers:               \t"+tables.kmersLoaded);
		if(shave || rinse){
			outstream.println("After Shaving:              \t"+(tables.kmersLoaded-kmersRemoved));
		}

		averageCount=tables.kmersIn*1.0/tables.kmersLoaded;
		double actualDepth=Tools.observedToActualCoverage(averageCount);
		double readDepth=(actualDepth*tables.basesIn)/(tables.kmersIn);
		
		outstream.println("Average Kmer Count:         \t"+Tools.format("%."+decimals+"f", averageCount));
		outstream.println("Estimated Kmer Depth:       \t"+Tools.format("%."+decimals+"f", actualDepth));
		outstream.println("Estimated Read Depth:       \t"+Tools.format("%."+decimals+"f", readDepth));
		outstream.println();
		
		outstream.println("Load Time:                  \t"+t);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Ref Processing        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Processes reference sequences for intersection analysis.
	 * Creates a 2D matrix comparing k-mer counts between input reads and reference,
	 * useful for contamination detection and coverage analysis.
	 */
	void processRef(){
		long[][] counts=intersectionST(tables, tables2, histMax, refMax, bidirectional);
		ByteStreamWriter bsw=new ByteStreamWriter(intersectionFile, overwrite, append, false);
		bsw.start();
		bsw.print("#Depth");
		for(int i=0; i<counts[0].length; i++){
			bsw.print("\t"+i+"-Copy");
		}
		bsw.nl();
		
		int maxDepth=counts.length-1;
		while(maxDepth>0 && shared.Vector.sum(counts[maxDepth])<1){
			maxDepth--;
		}
		
		for(int depth=0; depth<=maxDepth; depth++){
			if(histZeros || shared.Vector.sum(counts[depth])>0){
				bsw.print(depth);
				for(int copy=0; copy<counts[depth].length; copy++){
					bsw.tab();
					bsw.print(counts[depth][copy]);
				}
				bsw.nl();
			}
		}
		bsw.poisonAndWait();
	}
	
	//qualities asdkjasdkasdladeudns,u384Gnah&dhhsd
	/**
	 * Creates intersection matrix between two k-mer tables using array iteration.
	 * Compares k-mer counts from reads vs reference for contamination analysis.
	 * This is the first implementation approach using direct array access.
	 *
	 * @param tables Primary k-mer table set (typically from reads)
	 * @param tables2 Reference k-mer table set
	 * @param xMax Maximum count value for primary table dimension
	 * @param yMax Maximum count value for reference table dimension
	 * @param bidirectional Whether to count k-mers unique to reference
	 * @return 2D matrix of count intersections [primary_count][reference_count]
	 */
	public static long[][] intersectionST_first(final AbstractKmerTableSet tables, 
			final AbstractKmerTableSet tables2, final int xMax, final int yMax, final boolean bidirectional){

		long[][] counts=new long[xMax+1][yMax+1];
		
		Timer t=new Timer();
		tables2.process(t);
		if(tables.getClass()==KmerTableSet.class){
			
			KmerTableSet ktsRead=(KmerTableSet)tables;
			KmerTableSet ktsRef=(KmerTableSet)tables2;
			for(AbstractKmerTable akt : ktsRead.tables()){
				assert(akt.getClass()==HashArray1D.class);
				HashArray1D ha=(HashArray1D)akt;
				long[] kmers=ha.array();
				int[] values=ha.values();
				for(int i=0; i<kmers.length; i++){
					final int readCount=Tools.mid(0, values[i], xMax);
					if(readCount>0){
						final long kmer=kmers[i];
						final int refCount=Tools.mid(0, ktsRef.getCount(kmer), yMax);
						counts[readCount][refCount]++;
					}
				}
				HashForest victims=ha.victims();
				for(KmerNode node : victims){
					KmerNode1D n=(KmerNode1D)node;
					final int readCount=Tools.mid(0, n.value(), xMax);
					if(readCount>0){
						final long kmer=n.pivot();
						final int refCount=Tools.mid(0, ktsRef.getCount(kmer), yMax);
						counts[readCount][refCount]++;
					}
				}
			}
			
			if(bidirectional){
				for(AbstractKmerTable akt : ktsRef.tables()){
					assert(akt.getClass()==HashArray1D.class);
					HashArray1D ha=(HashArray1D)akt;
					long[] kmers=ha.array();
					int[] values=ha.values();
					for(int i=0; i<kmers.length; i++){
						final int refCount=Tools.mid(0, values[i], yMax);
						if(refCount>0){
							final long kmer=kmers[i];
							final int readCount=Tools.mid(0, ktsRead.getCount(kmer), xMax);
							if(readCount<1){counts[readCount][refCount]++;}
						}
					}
					HashForest victims=ha.victims();
					for(KmerNode node : victims){
						KmerNode1D n=(KmerNode1D)node;
						final int refCount=Tools.mid(0, n.value(), yMax);
						if(refCount>0){
							final long kmer=n.pivot();
							final int readCount=Tools.mid(0, ktsRead.getCount(kmer), xMax);
							if(readCount<1){counts[readCount][refCount]++;}
						}
					}
				}
			}
			
		}else if(tables.getClass()==KmerTableSetU.class){
			
			Kmer kmer=new Kmer(tables.kbig());
			KmerTableSetU ktsRead=(KmerTableSetU)tables;
			KmerTableSetU ktsRef=(KmerTableSetU)tables2;
			for(AbstractKmerTableU akt : ktsRead.tables()){
				assert(akt.getClass()==HashArrayU1D.class);
				HashArrayU1D ha=(HashArrayU1D)akt;
				int[] values=ha.values();
				for(int i=0; i<values.length; i++){
					final int readCount=Tools.mid(0, values[i], xMax);
					if(readCount>0){
						ha.fillKmer(i, kmer);
						final int refCount=Tools.mid(0, ktsRef.getCount(kmer), yMax);
						counts[readCount][refCount]++;
					}
				}
				HashForestU victims=ha.victims();
				for(KmerNodeU node : victims){
					KmerNodeU1D n=(KmerNodeU1D)node;
					final int readCount=Tools.mid(0, n.value(), xMax);
					if(readCount>0){
						n.fillKmer(kmer);
						final int refCount=Tools.mid(0, ktsRef.getCount(kmer), yMax);
						counts[readCount][refCount]++;
					}
				}
			}
			
			if(bidirectional){
				for(AbstractKmerTableU akt : ktsRef.tables()){
					assert(akt.getClass()==HashArrayU1D.class);
					HashArrayU1D ha=(HashArrayU1D)akt;
					int[] values=ha.values();
					for(int i=0; i<values.length; i++){
						final int refCount=Tools.mid(0, values[i], yMax);
						if(refCount>0){
							ha.fillKmer(i, kmer);
							final int readCount=Tools.mid(0, ktsRead.getCount(kmer), xMax);
							if(readCount<1){counts[refCount][readCount]++;}
						}
					}
					HashForestU victims=ha.victims();
					for(KmerNodeU node : victims){
						KmerNodeU1D n=(KmerNodeU1D)node;
						final int refCount=Tools.mid(0, n.value(), yMax);
						if(refCount>0){
							n.fillKmer(kmer);
							final int readCount=Tools.mid(0, ktsRead.getCount(kmer), xMax);
							if(readCount<1){counts[refCount][readCount]++;}
						}
					}
				}
			}
			
		}else{
			assert(false) : tables.getClass();
		}
		return counts;
	}
	
	/**
	 * Creates intersection matrix between two k-mer tables using walker iteration.
	 * Alternative implementation using table walkers for more efficient traversal.
	 * This is the second implementation approach using iterator pattern.
	 *
	 * @param tables Primary k-mer table set (typically from reads)
	 * @param tables2 Reference k-mer table set
	 * @param xMax Maximum count value for primary table dimension
	 * @param yMax Maximum count value for reference table dimension
	 * @param bidirectional Whether to count k-mers unique to reference
	 * @return 2D matrix of count intersections [primary_count][reference_count]
	 */
	public static long[][] intersectionST_second(final AbstractKmerTableSet tables, 
			final AbstractKmerTableSet tables2, final int xMax, final int yMax, final boolean bidirectional){

		long[][] counts=new long[xMax+1][yMax+1];
		
		Timer t=new Timer();
		tables2.process(t);
		if(tables.getClass()==KmerTableSet.class){
			
			KmerTableSet ktsRead=(KmerTableSet)tables;
			KmerTableSet ktsRef=(KmerTableSet)tables2;
			for(AbstractKmerTable akt : ktsRead.tables()){
				assert(akt.getClass()==HashArray1D.class);//If Walker is made abstract, these lines can be omitted
				HashArray1D ha=(HashArray1D)akt;
				
				Walker w=ha.walk();
				while(w.next()){
					final int readCount=Tools.mid(0, w.value(), xMax);
					final int refCount=Tools.mid(0, ktsRef.getCount(w.kmer()), yMax);
					counts[readCount][refCount]++;
				}
			}
			
			if(bidirectional){
				for(AbstractKmerTable akt : ktsRef.tables()){
					assert(akt.getClass()==HashArray1D.class);//If Walker is made abstract, these lines can be omitted
					HashArray1D ha=(HashArray1D)akt;
					
					Walker w=ha.walk();
					while(w.next()){
						final int refCount=Tools.mid(0, w.value(), yMax);
						final int readCount=Tools.mid(0, ktsRead.getCount(w.kmer()), xMax);
						if(readCount<1){counts[readCount][refCount]++;}
					}
				}
			}
			
		}else if(tables.getClass()==KmerTableSetU.class){
			
			KmerTableSetU ktsRead=(KmerTableSetU)tables;
			KmerTableSetU ktsRef=(KmerTableSetU)tables2;
			for(AbstractKmerTableU akt : ktsRead.tables()){
				assert(akt.getClass()==HashArrayU1D.class);
				HashArrayU1D ha=(HashArrayU1D)akt;
				
				WalkerU w=ha.walk();
				while(w.next()){
					final int readCount=Tools.mid(0, w.value(), xMax);
					final int refCount=Tools.mid(0, ktsRef.getCount(w.kmer()), yMax);
					counts[readCount][refCount]++;
				}
			}
			
			if(bidirectional){
				for(AbstractKmerTableU akt : ktsRef.tables()){
					assert(akt.getClass()==HashArrayU1D.class);
					HashArrayU1D ha=(HashArrayU1D)akt;
					
					WalkerU w=ha.walk();
					while(w.next()){
						final int refCount=Tools.mid(0, w.value(), yMax);
						final int readCount=Tools.mid(0, ktsRead.getCount(w.kmer()), xMax);
						if(readCount<1){counts[readCount][refCount]++;}
					}
				}
			}
			
		}else{
			assert(false) : tables.getClass();
		}
		return counts;
	}
	
	/**
	 * Creates intersection matrix between two k-mer tables using optimized walker approach.
	 * This is the current implementation that provides the best performance
	 * by using table-level walkers instead of individual table iteration.
	 *
	 * @param tables Primary k-mer table set (typically from reads)
	 * @param tables2 Reference k-mer table set
	 * @param xMax Maximum count value for primary table dimension
	 * @param yMax Maximum count value for reference table dimension
	 * @param bidirectional Whether to count k-mers unique to reference
	 * @return 2D matrix of count intersections [primary_count][reference_count]
	 */
	public static long[][] intersectionST(final AbstractKmerTableSet tables, 
			final AbstractKmerTableSet tables2, final int xMax, final int yMax, final boolean bidirectional){

		long[][] counts=new long[xMax+1][yMax+1];
		
		Timer t=new Timer();
		tables2.process(t);
		if(tables.getClass()==KmerTableSet.class){
			
			KmerTableSet ktsRead=(KmerTableSet)tables;
			KmerTableSet ktsRef=(KmerTableSet)tables2;
			{
				Walker w=ktsRead.walk();
				while(w.next()){
					final int readCount=Tools.mid(0, w.value(), xMax);
					final int refCount=Tools.mid(0, ktsRef.getCount(w.kmer()), yMax);
					counts[readCount][refCount]++;
				}
			}
			if(bidirectional){
				Walker w=ktsRef.walk();
				while(w.next()){
					final int refCount=Tools.mid(0, w.value(), yMax);
					final int readCount=Tools.mid(0, ktsRead.getCount(w.kmer()), xMax);
					if(readCount<1){counts[readCount][refCount]++;}
				}
			}
			
		}else if(tables.getClass()==KmerTableSetU.class){
			
			KmerTableSetU ktsRead=(KmerTableSetU)tables;
			KmerTableSetU ktsRef=(KmerTableSetU)tables2;
			{
				WalkerU w=ktsRead.walk();
				while(w.next()){
					final int readCount=Tools.mid(0, w.value(), xMax);
					final int refCount=Tools.mid(0, ktsRef.getCount(w.kmer()), yMax);
					counts[readCount][refCount]++;
				}
			}
			if(bidirectional){
				WalkerU w=ktsRef.walk();
				while(w.next()){
					final int refCount=Tools.mid(0, w.value(), yMax);
					final int readCount=Tools.mid(0, ktsRead.getCount(w.kmer()), xMax);
					if(readCount<1){counts[readCount][refCount]++;}
				}
			}
			
		}else{
			assert(false) : tables.getClass();
		}
		return counts;
	}
	
//	void processRef(){
//		FileFormat ffref=FileFormat.testInput(ref, FileFormat.FA, null, true, false);
//		Streamer cris=makeCris(ffref);
//		//Do anything necessary prior to processing
//
//		{
//			//Grab the first ListNum of reads
//			ListNum<Read> ln=cris.nextList();
//
//			//Check to ensure pairing is as expected
//			if(ln!=null && !ln.isEmpty()){
//				Read r=ln.get(0);
//				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired());
//			}
//
//			//As long as there is a nonempty read list...
//			while(ln!=null && ln.size()>0){
//				
//				processList(ln, cris);
//
//				//Fetch a new list
//				ln=cris.nextList();
//			}
//
//			//Notify the input stream that the final list was used
//			if(ln!=null){
//				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
//			}
//		}
//
//		//Do anything necessary after processing
//		ReadWrite.closeStream(cris);
//	}
//	
//	/**
//	 * Process a list of Reads.
//	 * @param ln The list.
//	 * @param cris Read Input Stream
//	 * @param ros Read Output Stream for reads that will be retained
//	 */
//	void processList(ListNum<Read> ln, final Streamer cris){
//
//		//Grab the actual read list from the ListNum
//		final ArrayList<Read> reads=ln.list;
//		
//		//Loop through each read in the list
//		for(int idx=0; idx<reads.size(); idx++){
//			final Read r1=reads.get(idx);
//			final Read r2=r1.mate;
//			
//			//Validate reads in worker threads
//			if(!r1.validated()){r1.validate(true);}
//			if(r2!=null && !r2.validated()){r2.validate(true);}
//
//			//Track the initial length for statistics
//			final int initialLength1=r1.length();
//			final int initialLength2=r1.mateLength();
//
//			//Increment counters
//			readsProcessed+=r1.pairCount();
//			basesProcessed+=initialLength1+initialLength2;
//			
//			{
//				//Reads are processed in this block.
//				boolean keep=processReadPair(r1, r2);
//				
//				if(!keep){reads.set(idx, null);}
//				else{
//					readsOut+=r1.pairCount();
//					basesOut+=r1.pairLength();
//				}
//			}
//		}
//
//		//Output reads to the output stream
//		if(ros!=null){ros.add(reads, ln.id);}
//
//		//Notify the input stream that the list was used
//		cris.returnList(ln);
////		if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access
//	}
//	
//	private Streamer makeCris(FileFormat ff){
//		Streamer cris=Streamer.getReadInputStream(maxReads, true, ff, null);
//		cris.start(); //Start the stream
//		if(verbose){outstream.println("Started ref cris");}
//		boolean paired=cris.paired();
//		assert(!paired);
//		return cris;
//	}
//	
//	
//	/**
//	 * Process a single contig.
//	 * @param r The contig
//	 */
//	void processContig(final Read r){
//		
//		final byte[] bases=r.bases;
//		long kmer=0;
//		long rkmer=0;
//		long added=0;
//		int len=0;
//		
//		for(int i=0; i<bases.length; i++){
//			final byte b=bases[i];
//			final long x=AminoAcid.baseToNumber[b];
//			final long x2=AminoAcid.baseToComplementNumber[b];
////			assert(x!=x2) : x+", "+x2+", "+Character.toString((char)b)+"\n"+Arrays.toString(symbolToNumber0)+"\n"+Arrays.toString(symbolToComplementNumber);
//			kmer=((kmer<<bitsPerBase)|x)&mask;
//			//10000, 1111111111, 16, 16, 2, 10, 8
//			rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
//			if(x>=0){len++;}else{len=0; rkmer=0;}
//			if(verbose){
//				if(verbose){
//					String fwd=new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k));
//					String rev=AminoAcid.reverseComplementBases(fwd);
//					String fwd2=kmerToString(kmer, Tools.min(len, k));
//					outstream.println("fwd="+fwd+", fwd2="+fwd2+", rev="+rev+", kmer="+kmer+", rkmer="+rkmer);
//					outstream.println("b="+(char)b+", x="+x+", x2="+x2+", bitsPerBase="+bitsPerBase+", shift2="+shift2);
//					if(!amino){
//						assert(AminoAcid.stringToKmer(fwd)==kmer) : fwd+", "+AminoAcid.stringToKmer(fwd)+", "+kmer+", "+len;
//						if(len>=k){
//							assert(AminoAcid.reverseComplementBinaryFast(kmer, Tools.min(len, k))==rkmer);
//							assert(AminoAcid.reverseComplementBinaryFast(rkmer, Tools.min(len, k))==kmer);
//							assert(AminoAcid.kmerToString(kmer, Tools.min(len, k)).equals(fwd));
//							assert(AminoAcid.kmerToString(rkmer, Tools.min(len, k)).equals(rev)) : AminoAcid.kmerToString(rkmer, Tools.min(len, k))+" != "+rev+" (rkmer)";
//						}
//						assert(fwd.equalsIgnoreCase(fwd2)) : fwd+", "+fwd2; //may be unsafe
//					}
//					outstream.println("Scanning6 i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+", bases="+fwd+", rbases="+rev);
//				}
//			}
//			if(len>=k){
////				assert(kmer==AminoAcid.reverseComplementBinaryFast(rkmer, k)) : Long.toBinaryString(kmer)+", "+Long.toBinaryString(rkmer)+", "+Long.toBinaryString(mask)+", x="+x+", x2="+x2+", bits="+bitsPerBase+", s="+shift+", s2="+shift2+", b="+Character.toString((char)b);
//				refKmersT++;
//				final long extraBase=(i>=bases.length-1 ? -1 : symbolToNumber[bases[i+1]]);
//				final long atm=addToMap(kmer, rkmer, k, extraBase, id, kmask, hammingDistance, editDistance);
//				added+=atm;
////				assert(false) : atm+", "+map.contains(toValue(kmer, rkmer, kmask));
//				if(useShortKmers){
//					if(i==k2){added+=addToMapRightShift(kmer, rkmer, id);}
//					if(i==bases.length-1){added+=addToMapLeftShift(kmer, rkmer, extraBase, id);}
//				}
//			}
//		}
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Removes low-abundance k-mers that are likely errors.
	 * Shaving removes tip structures while rinsing removes bubbles in the k-mer graph.
	 * This iterative process improves data quality by removing sequencing errors.
	 *
	 * @param shave Whether to perform tip removal
	 * @param rinse Whether to perform bubble removal
	 * @param maxShaveDepth Maximum number of shaving iterations to perform
	 * @return Total number of k-mers removed across all iterations
	 */
	long shave(boolean shave, boolean rinse, int maxShaveDepth){
		long sum=0;

		for(int i=0; i<maxShaveDepth; i++){
			int a=i+1, b=maxShaveDepth, c=i+1;
			//				if(i>3){Shaver2.verbose2=true;}
			outstream.println("\nShave("+a+", "+b+", "+c+")");
			final Shaver shaver=Shaver.makeShaver(tables, THREADS, a, b, c, 1, 3, 100, 100, shave, rinse);
			long removed=shaver.shave(a, b);
			sum+=removed;
		}

		System.err.println();
		return sum;
	}
	
	/**
	 * Generates k-mer count histogram and optionally calls peaks.
	 * Creates histogram showing distribution of k-mer frequencies,
	 * applies optional smoothing, and identifies coverage peaks for genome size estimation.
	 *
	 * @param fname Histogram output file path
	 * @param peaks Peak analysis output file path
	 * @param cols Number of histogram columns to generate
	 * @param max Maximum histogram length
	 * @param printHeader Whether to include column headers
	 * @param printZeros Whether to include zero-count bins
	 * @param printTime Whether to include timing information
	 * @param smoothKhist Whether to smooth the histogram data
	 * @param smoothPeaks Whether to smooth data before peak calling
	 * @return Average k-mer count from the histogram
	 */
	private double makeKhist(String fname, String peaks, int cols, int max, boolean printHeader, boolean printZeros, boolean printTime, boolean smoothKhist, boolean smoothPeaks){
		if(fname==null && peaks==null){return -1;}
		
		final long[][] arrays=tables.makeKhist(fname, cols, max, printHeader, printZeros, printTime, smoothKhist, gcHist, doLogScale, logWidth, logPasses, smoothRadius);
		final long[] array=arrays[0];
		final long[] gcArray=arrays[1];
		
		double avg=Tools.averageHistogram(array);
		
		if(peaks!=null){
			CallPeaks.printClass=false;
			ArrayList<String> args=new ArrayList<String>();
			if((smoothPeaks && !smoothKhist) && smoothRadius>0){//!smoothKhist because if smoothKhist is true the array will already be smoothed
				args.add("smoothradius="+smoothRadius);
				args.add("smoothprogressive=t");
			}
			CallPeaks.printPeaks(array, gcArray, peaks, overwrite, minHeight, minVolume, minWidth, 
					Tools.max(tables.filterMax()+2, minPeak), maxPeak, maxPeakCount, k, ploidy, doLogScale, logWidth, args);
		}
		return avg;
	}
	
	/**
	 * Creates a MinHash sketch from the k-mer count data.
	 * Sketches provide compact representations of k-mer sets for rapid comparison
	 * and are useful for taxonomic classification and sequence similarity estimation.
	 */
	private void makeSketch(){
		Timer ts=new Timer();
		outstream.println("Generating sketch.");
		SketchObject.maxGenomeFraction=1;
		SketchObject.k=k;
		SketchTool sketcher=new SketchTool(sketchLength, minToDump, false, false, SketchObject.rcomp);
		Sketch sketch=sketcher.toSketch((KmerTableSet)tables, true);
		if(sketch==null){
			errorState=true;
			System.err.println("WARNING: No sketch was produced, presumably because no kmers passed the filter criteria.");
			assert(false);
			return;
		}
		sketch.setName0(ReadWrite.stripToCore(ffSketch.name()));
		SketchTool.write(sketch, ffSketch);
		ts.stop();
		outstream.println("Sketch Time:                \t"+ts);
	}
	
	/**
	 * Coordinates histogram generation and k-mer dumping operations.
	 * Uses multithreading when possible to parallelize histogram creation
	 * and k-mer output for improved performance on large datasets.
	 */
	private void makeKhistAndPeaks(){
		if(THREADS>1 && (outHist!=null || outPeaks!=null) && outKmers!=null){
			Timer tout=new Timer();
			tout.start();
			Thread a=new DumpKmersThread();
			Thread b=new MakeKhistThread();
			a.start();
			b.start();
			while(a.getState()!=Thread.State.TERMINATED){
				try {
					a.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			while(b.getState()!=Thread.State.TERMINATED){
				try {
					b.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			tout.stop();
			outstream.println("Write Time:                 \t"+tout);
		}else{
			if(outHist!=null || outPeaks!=null){
				averageCount=makeKhist(outHist, outPeaks, histColumns, histMax, histHeader, histZeros, true, smoothKhist, smoothPeaks);
			}
			if(outKmers!=null){
				//			tables.dumpKmersAsText(outKmers, minToDump, true);
				tables.dumpKmersAsBytes_MT(outKmers, minToDump, maxToDump, true, null);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Worker thread for dumping k-mers to output files.
	 * Allows k-mer output to proceed in parallel with histogram generation
	 * for improved overall performance.
	 */
	private class DumpKmersThread extends Thread {
		
		/** Default constructor for k-mer dumping thread */
		DumpKmersThread(){}
		
		@Override
		public void run(){
			tables.dumpKmersAsBytes_MT(outKmers, minToDump, maxToDump, false, null);
		}
		
	}
	
	/**
	 * Worker thread for histogram generation.
	 * Allows histogram creation to proceed in parallel with k-mer output
	 * for improved overall performance.
	 */
	private class MakeKhistThread extends Thread {
		
		/** Default constructor for histogram generation thread */
		MakeKhistThread(){}
		
		@Override
		public void run(){
			makeKhist(outHist, outPeaks, histColumns, histMax, histHeader, histZeros, false, smoothKhist, smoothPeaks);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Hold kmers. */
	private final AbstractKmerTableSet tables;
	/**
	 * Secondary k-mer table set for reference sequences used in intersection analysis
	 */
	private AbstractKmerTableSet tables2;//ref kmers
	
	/** Whether to perform tip shaving to remove likely sequencing errors */
	private boolean shave=false;
	/** Whether to perform bubble rinsing to remove likely sequencing errors */
	private boolean rinse=false;
	/** Maximum depth for iterative shaving operations */
	private int shaveDepth=1;
	
	/** Average count per unique k-mer in the dataset */
	private double averageCount=-1;
	/** Total number of bases processed from input sequences */
	private long basesIn=-1;
	/** Total number of reads processed from input files */
	private long readsIn=-1;
	/** Number of decimal places for displaying statistics */
	private int decimals=3;
	
	/** Number of k-mers removed during shaving and rinsing operations */
	private long kmersRemoved=0;
	
	/** Kmer count output file */
	private String outKmers=null;
	/** Histogram output file */
	private String outHist=null;
	/** Histogram peak output file */
	private String outPeaks=null;
	
	/** Radius for smoothing operations on histograms and peaks */
	private int smoothRadius=1;
	/** Whether to apply smoothing to the k-mer histogram */
	private boolean smoothKhist=false;
	/** Whether to apply smoothing before peak calling */
	private boolean smoothPeaks=false;
	
	/** Flag indicating whether errors occurred during processing */
	private boolean errorState=false;
	
	/** Histogram columns */
	private int histColumns=2;
	/** Histogram rows */
	private int histMax=100000;
	/** Print a histogram header */
	private boolean histHeader=true;
	/** Histogram show rows with 0 count */
	private boolean histZeros=false;
	/** Add gc information to kmer histogram */
	protected boolean gcHist=false;
	
	/** Whether to use logarithmic scaling for histogram bins */
	boolean doLogScale=true;
	/** Width parameter for logarithmic histogram binning */
	double logWidth=0.1;
	/** Number of passes for logarithmic histogram processing */
	int logPasses=1;
	
	/** Minimum height requirement for peak calling */
	private long minHeight=2;
	/** Minimum volume (area under curve) requirement for peak calling */
	private long minVolume=5;
	/** Minimum width requirement for peak calling */
	private int minWidth=3;
	/** Minimum count value for peak detection */
	private int minPeak=2;
	/** Maximum count value for peak detection */
	private int maxPeak=Integer.MAX_VALUE;
	/** Maximum number of peaks to identify */
	private int maxPeakCount=12;
	
	/** Expected ploidy level for peak calling analysis */
	private int ploidy=-1;
	
	/** Output file path for MinHash sketch */
	private String sketchPath=null;
	/** Target length for MinHash sketch */
	private int sketchLength=10000;
	/** Name identifier for the generated sketch */
	private String sketchName;
	/** Numeric identifier for the generated sketch */
	private int sketchID;
	/** File format handler for sketch output */
	private final FileFormat ffSketch;
	
	/*--------------------------------------------------------------*/
	/*----------------         Ref Counting         ----------------*/
	/*--------------------------------------------------------------*/

	/** Path to reference file for intersection analysis */
	private String ref=null;
	/** Output file path for intersection analysis results */
	private String intersectionFile=null;
	/** Whether to perform bidirectional intersection analysis */
	private boolean bidirectional=true;
	

	/** Maximum reference count value for intersection matrix */
	private int refMax=6;
	
	/*--------------------------------------------------------------*/
	/*----------------       Final Primitives       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** min kmer count to dump to text */
	private int minToDump=1;
	/** Maximum k-mer count allowed for inclusion in output */
	private int maxToDump=Integer.MAX_VALUE;

	/** K-mer length used for counting operations */
	final int k;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print messages to this stream */
	private static PrintStream outstream=System.err;
	/** Permission to overwrite existing files */
	public static boolean overwrite=true;
	/** Permission to append to existing files */
	public static boolean append=false;
	/** Display progress messages such as memory usage */
	public static boolean DISPLAY_PROGRESS=true;
	/** Verbose messages */
	public static final boolean verbose=false;
	/** Number of ProcessThreads */
	public static int THREADS=Shared.threads();

	
}
