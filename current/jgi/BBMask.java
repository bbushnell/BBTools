package jgi;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ConcurrentHashMap;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import structures.CoverageArray;
import structures.CoverageArray2;
import structures.CoverageArray3;
import structures.IntList;
import structures.ListNum;
import tracker.EntropyTracker;
import tracker.ReadStats;

/**
 * Loads reference sequences and marks optional repeats, low entropy/diversity,
 * and regions selected by SAM alignments or aggregate coverage thresholds.
 * Entropy masking is enabled by default; repeat masking is optional. Marking
 * accumulates in per-reference BitSets before bases are converted to N/lowercase
 * or masked spans are removed by splitting. Counts exclude marks already present
 * when each stage starts, including preexisting N/lowercase marks seeded on input.
 * 
 * @author Brian Bushnell
 * @date Feb 18, 2014
 *
 */
public class BBMask{

	/** Program entry point. Creates BBMask instance and processes input files.
	 * @param args Command-line arguments specifying input files and parameters */
	public static void main(String[] args){
		Timer t=new Timer();
		BBMask x=new BBMask(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructor that parses command-line arguments and initializes masking parameters.
	 * Sets up file formats, k-mer ranges, coverage thresholds, and processing modes.
	 * Also changes shared parser, compression, interleaving and SAM parsing settings;
	 * construction does not load the reference, start masking workers or write output.
	 * @param args Command-line arguments for configuration
	 */
	public BBMask(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(16);
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		
		boolean setEntropyMode=false;
		
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			//Only the first '=' separates the option key; later '=' characters belong to its value.
			String[] split=Parse.splitOnFirst(arg, '=');
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
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("t") || a.equals("threads")){
				Shared.setThreads(b);
			}else if(a.equals("sampad") || a.equals("sampadding") || a.equals("sp")){
				samPad=Integer.parseInt(b);
			}else if(a.equals("entropymode")){
				entropyMode=Parse.parseBoolean(b);
				setEntropyMode=true;
			}else if(a.equals("maskrepeats") || a.equals("mr")){
				processRepeats=Parse.parseBoolean(b);
			}else if(a.equals("masklowentropy") || a.equals("masklowcomplexity") || a.equals("mlc") || a.equals("mle") || a.equals("me")){
				processEntropy=Parse.parseBoolean(b);
			}else if(a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1") || a.equals("ref")){
				inRef=b;
			}else if(a.equals("insam") || a.equals("samin") || a.equals("sam")){
				if(b!=null){
					for(String s : b.split(",")){
						inSam.add(s);
					}
				}
			}else if(a.equals("out") || a.equals("output") || a.equals("out1") || a.equals("output1") || a.equals("output1")){
				outRef=b;
			}else if(a.equals("qfin") || a.equals("qfin1")){
				qfinRef=b;
			}else if(a.equals("qfout") || a.equals("qfout1")){
				qfoutRef=b;
			}else if(a.equals("extin")){
				extinRef=b;
			}else if(a.equals("extout")){
				extoutRef=b;
			}else if(a.equals("split")){
				splitMode=Parse.parseBoolean(b);
			}else if(a.equals("mink") || a.equals("kmin")){
				mink=mink2=Integer.parseInt(b);
			}else if(a.equals("maxk") || a.equals("kmax")){
				maxk=maxk2=Integer.parseInt(b);
			}else if(a.equals("k")){
				mink=maxk=mink2=maxk2=Integer.parseInt(b);
			}else if(a.equals("minkr") || a.equals("krmin")){
				mink=Integer.parseInt(b);
			}else if(a.equals("maxkr") || a.equals("krmax")){
				maxk=Integer.parseInt(b);
			}else if(a.equals("kr")){
				mink=maxk=Integer.parseInt(b);
			}else if(a.equals("mink2") || a.equals("kmin2") || a.equals("minke") || a.equals("kemin")){
				mink2=Integer.parseInt(b);
			}else if(a.equals("maxk2") || a.equals("kmax2") || a.equals("maxke") || a.equals("kemax")){
				maxk2=Integer.parseInt(b);
			}else if(a.equals("k2") || a.equals("ke")){
				mink2=maxk2=Integer.parseInt(b);
			}else if(a.equals("mincov")){
				mincov=Integer.parseInt(b);
			}else if(a.equals("maxcov")){
				maxcov=Integer.parseInt(b);
			}else if(a.equals("delcov") || a.equals("delcoverage")){
				includeDeletionCoverage=Parse.parseBoolean(b);
			}else if(a.equals("window") || a.equals("w")){
				window=Integer.parseInt(b);
			}else if(a.equals("ratio")){
				ratio=Float.parseFloat(b);
				if(!setEntropyMode){entropyMode=false;}
			}else if(a.equals("entropy") || a.equals("e")){
				entropyCutoff=Float.parseFloat(b);
				if(!setEntropyMode){entropyMode=true;}
			}else if(a.equals("lowercase") || a.equals("lc")){
				MaskByLowercase=Parse.parseBoolean(b);
			}else if(a.equals("minlen")){
				minlen=Integer.parseInt(b);
			}else if(a.equals("mincount")){
				mincount=Integer.parseInt(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.startsWith("minscaf") || a.startsWith("mincontig")){
				stream.FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
			}else if(!arg.contains("=") && FileFormat.hasSamOrBamExtension(arg)){
				inSam.add(arg);
			}else if(inRef==null && i==0 && Tools.looksLikeInputStream(arg)){
				inRef=arg;
			}else{
				System.err.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		bits32=(mincov>=Character.MAX_VALUE || maxcov>=Character.MAX_VALUE);
		
		{//Process parser fields
			Parser.processQuality();
		}

		assert(FastaReadInputStream.settingsOK());

		if(inRef==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}

		if(outRef!=null && outRef.equalsIgnoreCase("null")){outRef=null;}

		if(!Tools.testOutputFiles(overwrite, append, false, outRef)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+outRef+"\n");
		}
		
		ffoutRef=FileFormat.testOutput(outRef, FileFormat.FASTA, extoutRef, true, overwrite, append, false);

		ffinRef=FileFormat.testInput(inRef, FileFormat.FASTA, extinRef, true, true);
		
		if(inSam!=null && inSam.size()>0){
			ffinSam=new FileFormat[inSam.size()];
			for(int i=0; i<inSam.size(); i++){
				ffinSam[i]=FileFormat.testInput(inSam.get(i), FileFormat.SAM, null, true, false);
			}
		}else{
			ffinSam=null;
		}
		
		SamLine.CONVERT_CIGAR_TO_MATCH=false;
		
		EntropyTracker.verify=verify;
	}
	
	/*--------------------------------------------------------------*/


	/**
	 * Loads references, then adds repeat, entropy/diversity and SAM marks in that order.
	 * Each stage counts newly set bits, so overlapping marks count only once. Input
	 * N/lowercase marks are excluded from these totals. Conversion always runs, even
	 * with no new marks: it changes bases or replaces the map with split fragments.
	 * Unmarked non-ACGTN conversion can change bases without increasing the total.
	 * Writes the resulting map only when output is configured. This is a mutating
	 * processing invocation, not a reset/reuse or general failure-cleanup contract.
	 * @param t0 Started timer stopped when final totals are reported
	 */
	public void process(Timer t0){
		
		Timer t=new Timer();
		{
			t.start();
			outstream.println("Loading input");
			
			boolean oldTI=FASTQ.TEST_INTERLEAVED;
			boolean oldFI=FASTQ.FORCE_INTERLEAVED;
			FASTQ.TEST_INTERLEAVED=false;
			FASTQ.FORCE_INTERLEAVED=false;
			map=hashRef();
			FASTQ.TEST_INTERLEAVED=oldTI;
			FASTQ.FORCE_INTERLEAVED=oldFI;
			t.stop();
			
			outstream.println("Loading Time:                 \t"+t);
		}
		
		long repeats=0, lowcomplexity=0;
		long mapping=0;

		if(processRepeats && maxk>0){
			t.start();
			outstream.println("\nMasking repeats (to disable, set 'mr=f')");
			repeats=maskRepeats();
			t.stop();
			
			double rpnano=refReads/(double)(t.elapsed);
			double bpnano=refBases/(double)(t.elapsed);
			
			String rpstring=""+refReads;
			String bpstring=""+refBases;
			String bmstring=""+repeats;
	
			while(rpstring.length()<12){rpstring=" "+rpstring;}
			while(bpstring.length()<12){bpstring=" "+bpstring;}
			while(bmstring.length()<12){bmstring=" "+bmstring;}
	
			outstream.println("Repeat Masking Time:          \t"+t);
			outstream.println("Ref Bases:              "+bpstring+" \t"+Tools.format("%.2fm bases/sec", bpnano*1000));
			outstream.println("Repeat Bases Masked:    "+bmstring);
		}

		if(processEntropy && maxk2>0){
			t.start();
			if(entropyMode){
				outstream.println("\nMasking low-entropy (to disable, set 'mle=f')");
				lowcomplexity=maskLowEntropy();
			}else{
				outstream.println("\nMasking low-complexity (to disable, set 'mlc=f')");
				lowcomplexity=maskLowComplexity(null);
			}
			t.stop();
			
			double rpnano=refReads/(double)(t.elapsed);
			double bpnano=refBases/(double)(t.elapsed);
			
			String rpstring=""+refReads;
			String bpstring=""+refBases;
			String bmstring=""+lowcomplexity;
	
			while(rpstring.length()<12){rpstring=" "+rpstring;}
			while(bpstring.length()<12){bpstring=" "+bpstring;}
			while(bmstring.length()<12){bmstring=" "+bmstring;}
	
			outstream.println("Low Complexity Masking Time:  \t"+t);
			outstream.println("Ref Bases:              "+bpstring+" \t"+Tools.format("%.2fm bases/sec", bpnano*1000));
			outstream.println("Low Complexity Bases:   "+bmstring);
		}
		
		if(ffinSam!=null){
			t.start();
			outstream.println("\nMasking from sam");
			mapping=maskSam();
			t.stop();
			
			double rpnano=samReads/(double)(t.elapsed);
			double bpnano=samBases/(double)(t.elapsed);
			
			String rpstring=""+samReads;
			String bpstring=""+samBases;
			String bmstring=""+mapping;
	
			while(rpstring.length()<12){rpstring=" "+rpstring;}
			while(bpstring.length()<12){bpstring=" "+bpstring;}
			while(bmstring.length()<12){bmstring=" "+bmstring;}
			
			outstream.println("Sam Masking Time:             \t"+t);
			outstream.println("Sam Reads Processed:    "+rpstring+" \t"+Tools.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Sam Bases Processed:    "+bpstring+" \t"+Tools.format("%.2fm bases/sec", bpnano*1000));
			outstream.println("Sam Bases Masked:       "+bmstring);
		}
		long total=repeats+mapping+lowcomplexity, masked=0;
		
		if(total>0 || true){
			t.start();
			if(splitMode){
				masked=splitFromBitsets();
			}else{
				masked=maskFromBitsets(MaskByLowercase);
			}
			t.stop();
			outstream.println("Conversion Time:              \t"+t);
		}
		
		assert(total==masked) : repeats+", "+mapping+", "+lowcomplexity+", "+total+", "+masked;
		
		if(outRef!=null){
			t.start();
			outstream.println("\nWriting output");
			writeOutput();
			t.stop();
			outstream.println("Writing Time:                 \t"+t);
		}
		{
			t0.stop();
			String tstring=""+total;
			while(tstring.length()<12){tstring=" "+tstring;}
			outstream.println("\nTotal Bases Masked:     "+tstring+"/"+refBases+Tools.format("\t%.3f%%", total*100.0/refBases));
			outstream.println("Total Time:                   \t"+t0);
		}
		
		
		
		if(errorState){
			throw new RuntimeException("\nBBMask terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Sets bits in BitSet for positions exceeding maximum coverage threshold.
	 *
	 * @param bs BitSet to mark high coverage positions
	 * @param ca Coverage array containing depth values
	 * @param maxAllowedCoverage Maximum allowed coverage value
	 * @param maxLen Length of sequence to process
	 * @return Number of threshold-matching positions, including bits already set
	 */
	private static int setHighCoverage(BitSet bs, CoverageArray ca, int maxAllowedCoverage, int maxLen){
		int numSet=0;
		for(int i=0; i<maxLen; i++){
			if(ca.get(i)>maxAllowedCoverage){
				bs.set(i);
				numSet++;
			}
		}
		return numSet;
	}
	
	/**
	 * Sets bits in BitSet for positions below minimum coverage threshold.
	 *
	 * @param bs BitSet to mark low coverage positions
	 * @param ca Coverage array containing depth values
	 * @param minAllowedCoverage Minimum required coverage value
	 * @param maxLen Length of sequence to process
	 * @return Number of threshold-matching positions, including bits already set
	 */
	private static int setLowCoverage(BitSet bs, CoverageArray ca, int minAllowedCoverage, int maxLen){
		int numSet=0;
		for(int i=0; i<maxLen; i++){
			if(ca.get(i)<minAllowedCoverage){
				bs.set(i);
				numSet++;
			}
		}
		return numSet;
	}
	
	/**
	 * Mutates bases in the current reference map using each Read's attached BitSet.
	 * The BitSets are unchanged; counts follow maskBases, including its exclusions
	 * for preexisting masks and unmarked non-ACGTN conversion.
	 * @param lowercase If true, mask with lowercase letters; if false, use 'N'
	 * @return Sum of counted marked-base conversions across all sequences
	 */
	private long maskFromBitsets(final boolean lowercase){
		System.err.println("\nConverting masked bases to "+(lowercase ? "lower case" : "N")); //123
		long sum=0;
		for(Read r : map.values()){
			sum+=maskRead(r, (BitSet)r.obj, lowercase);
		}
		System.err.println("Done Masking");
		return sum;
	}
	
	/**
	 * Applies maskBases to a Read's bases without replacing its arrays or metadata.
	 *
	 * @param r Nonnull Read with a nonnull bases array to mutate
	 * @param bs Nonnull mask, read without modification
	 * @param lowercase If true, convert to lowercase; if false, use 'N'
	 * @return Count from maskBases, excluding its preexisting masks and off-bit conversions
	 */
	public static int maskRead(final Read r, final BitSet bs, final boolean lowercase){
		return maskBases(r.bases, bs, lowercase);
	}
	
	/** Mutates marked bases in place, leaving the BitSet unchanged.
	 * N mode counts each marked byte other than uppercase N, then writes uppercase N.
	 * Lowercase mode counts marked bytes that are neither lowercase nor uppercase N,
	 * then lowercases every marked byte, including N. If CONVERT_NON_ACGTN is enabled,
	 * unmarked non-ACGTN bytes become uppercase N without contributing to the count.
	 * Only positions within bases are visited; higher mask bits are ignored.
	 * @param bases Nonnull sequence array to mutate
	 * @param bs Nonnull mask, read without modification
	 * @param lowercase Whether to lowercase marked bytes instead of writing N
	 * @return Counted marked positions, excluding preexisting masks as described above
	 */
	public static int maskBases(final byte[] bases, final BitSet bs, final boolean lowercase){
		int sum=0;
		if(!lowercase){
			for(int i=0; i<bases.length; i++){
				if(bs.get(i)){
					if(bases[i]!='N'){sum++;}
					bases[i]='N';
				}else if(CONVERT_NON_ACGTN && !AminoAcid.isACGTN(bases[i])){
					bases[i]='N';
				}
			}
		}else{
			for(int i=0; i<bases.length; i++){
				if(bs.get(i)){
					if(!Tools.isLowerCase(bases[i]) && bases[i]!='N'){sum++;}
					bases[i]=(byte)Tools.toLowerCase(bases[i]);
				}else if(CONVERT_NON_ACGTN && !AminoAcid.isACGTN(bases[i])){
					bases[i]='N';
				}
			}
		}
		return sum;
	}
	
	/**
	 * Splits reads by removing masked regions, creating multiple shorter reads
	 * from each original read where masked regions act as break points.
	 * Generated suffixes skip all original IDs and previously emitted IDs.
	 * @return Newly masked bases, excluding preexisting N and, in lowercase mode, lowercase bases
	 */
	private long splitFromBitsets(){
		System.err.println("\nSplitting reads by removing masked bases"); //123
		long sum=0;
		
		LinkedHashMap<String, Read> map2=new LinkedHashMap<String, Read>();

		for(String key : map.keySet()){
			Read r=map.get(key);
			BitSet bs=((BitSet)r.obj);
			long rnum=0;
			if(bs.isEmpty()){
				Read old=map2.put(key, r);
				assert(old==null) : "Original ID must remain reserved while splitting: "+key;
			}else{
				byte[] bases=r.bases;
				byte[] quals=r.quality;

				int lastGood, lastBad;
				if(bs.get(0)){
					if(bases[0]!='N' && (!MaskByLowercase || !Tools.isLowerCase(bases[0]))){sum++;}
					lastGood=-1;
					lastBad=0;
				}else{
					lastGood=0;
					lastBad=-1;
				}
				
				int i=1;
				for(; i<bases.length; i++){
					if(bs.get(i)){
						if(bases[i]!='N' && (!MaskByLowercase || !Tools.isLowerCase(bases[i]))){sum++;}
						if(lastGood==i-1){
							int len=lastGood-lastBad;
							if(len>0){
								byte[] bases2=KillSwitch.copyOfRange(bases, lastBad+1, i);
								byte[] quals2=(quals==null ? null : KillSwitch.copyOfRange(quals, lastBad+1, i));
								String id=r.id+"_"+rnum++;
								while(map.containsKey(id) || map2.containsKey(id)){id=r.id+"_"+rnum++;}
								Read r2=new Read(bases2, quals2, id, r.numericID, r.flags);
								Read old=map2.put(r2.id, r2);
								assert(old==null) : "Split ID selection must prevent replacing an emitted fragment: "+r2.id;
							}
						}
						lastBad=i;
					}else{
						lastGood=i;
					}
				}
				if(lastGood==i-1){
					int len=lastGood-lastBad;
					if(len>0){
						byte[] bases2=KillSwitch.copyOfRange(bases, lastBad+1, i);
						byte[] quals2=(quals==null ? null : KillSwitch.copyOfRange(quals, lastBad+1, i));
						String id=r.id+"_"+rnum++;
						while(map.containsKey(id) || map2.containsKey(id)){id=r.id+"_"+rnum++;}
						Read r2=new Read(bases2, quals2, id, r.numericID, r.flags);
						Read old=map2.put(r2.id, r2);
						assert(old==null) : "Split ID selection must prevent replacing an emitted fragment: "+r2.id;
					}
				}
			}
		}
		
		map.clear();
		map.putAll(map2);
		map2.clear();
		
		System.err.println("Done Splitting");
		return sum;
	}
	
	/** Writes processed sequences to output file using concurrent stream writer. */
	private void writeOutput(){
		
		ConcurrentReadOutputStream ros=null;
		if(ffoutRef!=null){
			final int buff=16;
			ros=ConcurrentReadOutputStream.getStream(ffoutRef, null, qfoutRef, null, buff, null, false);
			ros.start();
		}
		
		long i=0;
		for(String name : map.keySet()){
			Read r=map.get(name);
			ArrayList<Read> list=new ArrayList<Read>(1);
			list.add(r);
			if(ros!=null){ros.add(list, i);}
			i++;
		}
		errorState|=ReadWrite.closeStream(ros);
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Processes all SAM files to mark mapped regions for masking.
	 * Accumulates coverage across all files before applying coverage thresholds.
	 * @return Number of additional bases marked for masking
	 */
	private long maskSam(){
		long before=0, after=0;
		for(Read r : map.values()){
			before+=((BitSet)r.obj).cardinality();
		}
		for(FileFormat ff : ffinSam){
			maskSam_MT(ff);
		}
		if(covmap!=null){//Move aggregate coverage information to BitSets after all inputs close.
			for(String rs : map.keySet()){
				Read r=map.get(rs);
				BitSet bs=(BitSet)r.obj;
				CoverageArray ca=covmap.remove(rs);
				if(maxcov>-1){
					setHighCoverage(bs, ca, maxcov, r.length());
				}
				if(mincov>-1){
					setLowCoverage(bs, ca, mincov, r.length());
				}
			}
		}
		for(Read r : map.values()){
			after+=((BitSet)r.obj).cardinality();
		}
		return after-before;
	}
	
	/** Multi-threaded SAM processing for a single file format.
	 * Observed owner interruption is reported after joining workers and closing input.
	 * @param ff FileFormat specification for SAM input file */
	private void maskSam_MT(FileFormat ff){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff, null, null, null);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		
		MaskSamThread[] threads=new MaskSamThread[Shared.threads()];
		for(int i=0; i<threads.length; i++){threads[i]=new MaskSamThread(cris);}
		for(int i=0; i<threads.length; i++){threads[i].start();}
		final InterruptedException interrupted=awaitWorkers(threads);
		try{
			errorState|=ReadWrite.closeStreams(cris);
			if(errorState){
				throw new RuntimeException("BBMask terminated in an error state; the output may be corrupt.");
			}
		}catch(RuntimeException | Error failure){
			if(interrupted!=null){failure.addSuppressed(interrupted);}
			throw failure;
		}finally{
			//Input cleanup may itself wait; restore only after that cleanup has finished.
			if(interrupted!=null){Thread.currentThread().interrupt();}
		}
		if(interrupted!=null){
			throw new RuntimeException("BBMask interrupted while waiting for SAM workers.", interrupted);
		}
	}
	
	/**
	 * Joins this invocation's already-started workers before caller cleanup.
	 * Does not restore interrupt status; callers restore it after closing their resources.
	 * @param threads Nonnull workers started by the caller, excluding the calling thread
	 * @return First interruption observed during joins, with later ones suppressed, or null
	 */
	private static InterruptedException awaitWorkers(final Thread[] threads){
		InterruptedException interrupted=null;
		for(final Thread thread : threads){
			assert(thread!=null && thread!=Thread.currentThread() && thread.getState()!=Thread.State.NEW) :
				"BBMask starts every owned worker before awaitWorkers; a null, self, or NEW entry cannot be joined: "+thread;
			while(thread.getState()!=Thread.State.TERMINATED){
				try{
					thread.join();
				}catch(InterruptedException e){
					if(interrupted==null){interrupted=e;}else if(interrupted!=e){interrupted.addSuppressed(e);}
				}
			}
		}
		return interrupted;
	}

	/** Worker thread for multi-threaded SAM file processing.
	 * Each thread processes SAM reads from shared input stream. */
	private class MaskSamThread extends Thread{
		
		/** Constructor for SAM processing thread.
		 * @param cris_ Concurrent input stream for reading SAM data */
		MaskSamThread(ConcurrentReadInputStream cris_){
			cris=cris_;
		}
		
		@Override
		public void run(){
			try{
				maskSam(cris);
			}catch(Throwable t){
				//Abandoned consumers can leave raw producers blocked, making ordinary stream close unsafe.
				t.printStackTrace(System.err);
				KillSwitch.kill("BBMask SAM worker failed; aborting because input cleanup cannot safely recover from worker failure.");
			}
		}
		
		/** Input stream for reading SAM format data */
		final ConcurrentReadInputStream cris;
		
		
	}
	
	/** Single-threaded SAM processing for a single file format.
	 * @param ff FileFormat specification for SAM input file */
	private void maskSam_ST(FileFormat ff){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff, null, null, null);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		
		maskSam(cris);

		errorState|=ReadWrite.closeStreams(cris);

		if(errorState){
			throw new RuntimeException("BBMask terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * Core SAM processing logic that reads mapped reads and marks reference
	 * positions for masking based on alignment coordinates.
	 * Attempts to return each borrowed batch once, including on processing failure. If returning
	 * the batch also fails, that cleanup error is suppressed on the processing error.
	 * @param cris Input stream for reading SAM format data
	 */
	private void maskSam(ConcurrentReadInputStream cris){
		
		long samReads=0;
		long samBases=0;
		IntList ranges=new IntList(16);
		ListNum<Read> ln=null;
		Throwable failure=null;
		try{
			ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

				for(int idx=0; idx<reads.size(); idx++){
					final Read r=reads.get(idx);
					assert(r.mate==null);

					final int initialLength1=r.length();

					{
						samReads++;
						samBases+=initialLength1;
						
						if(r.mapped()){
							SamLine sl=r.samline;
							assert(sl!=null) : "No sam line for read "+r;
							byte[] rname=sl.rname();
							assert(rname!=null) : "No rname for sam line "+sl;
							final String rs=new String(rname);
							Read ref=map.get(rs);
							if(ref==null){
								handleNoRef(rs);
							}else{
								final int reflen=ref.length();
								assert(ref!=null) : "Could not find reference scaffold '"+rs+"' for samline \n"+sl+"\n in set \n"+map.keySet()+"\n";
								if(covmap==null){
									BitSet bs=(BitSet)ref.obj;
									mask(bs, sl, reflen);
								}else{
									CoverageArray ca=covmap.get(rs);
									increment(ca, sl, r.match, reflen, ranges, includeDeletionCoverage, samPad);
								}
							}
						}
					}
				}

				final ListNum<Read> returned=ln;
				ln=null;
				cris.returnList(returned);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
		}catch(Throwable t){
			failure=t;
			throw t;
		}finally{
			if(ln!=null){
				try{
					cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
				}catch(Throwable t){
					if(failure==null){throw t;}
					if(failure!=t){failure.addSuppressed(t);}
				}
			}
		}
		
		synchronized(this){
			this.samBases+=samBases;
			this.samReads+=samReads;
		}
	}
	
	/**
	 * Marks reference positions covered by a SAM alignment for masking.
	 * Reconstructs the full soft-clipped span before applying padding and reference bounds.
	 * @param bs BitSet to mark positions for masking
	 * @param sl SAM line containing alignment information
	 * @param reflen Length of reference sequence
	 */
	private void mask(BitSet bs, SamLine sl, int reflen){
		final int rawStart=sl.start(true, false);
		final int start=Tools.max(0, rawStart-samPad);
		final int stop=Tools.min(sl.stop(rawStart, true, false)+1+samPad, reflen);
		if(stop>start){
			synchronized(bs){//Potential bottleneck
				bs.set(start, stop);
			}
		}
	}
	
	/**
	 * Increments coverage values for aligned regions in a coverage array.
	 * Uses the unclipped origin while decoding matches, then bounds emitted ranges.
	 *
	 * @param ca Coverage array to increment
	 * @param sl SAM line with alignment information
	 * @param match Match string for alignment details
	 * @param reflen Reference sequence length
	 * @param ranges IntList for storing alignment ranges
	 * @param includeDels Whether to include deletion coverage
	 * @param samPad Padding around aligned regions
	 */
	public static void increment(CoverageArray ca, SamLine sl, byte[] match, int reflen, IntList ranges, boolean includeDels, int samPad){
		final int start=sl.start(true, false);
		final int stop=sl.stop(start, true, false)+1;
		if(stop>start && stop>0 && start<reflen){
			ranges.clear();
			boolean numeric=false;
			if(match==null){
				assert(sl.cigar!=null);
				match=sl.toShortMatch(true);
				numeric=true;
			}else{
				for(byte b : match){
					if(Tools.isDigit(b)){
						numeric=true;
						break;
					}
				}
			}
			if(numeric){match=Read.toLongMatchString(match);}
			fillRanges(match, start, stop, ranges, includeDels);
			if(ranges.size>0){
				if(samPad!=0 || start<0 || stop>reflen){//Pad and bound the ranges, but don't let them overlap
					ranges.set(0, Tools.mid(0, reflen, ranges.get(0)-samPad));
					ranges.set(1, Tools.mid(0, reflen, ranges.get(1)+samPad));
					for(int i=2; i<ranges.size; i+=2){
						ranges.set(i, Tools.mid(ranges.get(i-1), reflen, ranges.get(i)-samPad));
						ranges.set(i+1, Tools.mid(0, reflen, ranges.get(i+1)+samPad));
					}
				}
				synchronized(ca){//Potential bottleneck
					ca.incrementRanges(ranges, 1);
				}
			}
		}
	}
	
	/**
	 * Parses alignment match string to identify covered reference ranges.
	 *
	 * @param longmatch Long-format match string from alignment
	 * @param start Zero-based inclusive starting position in the reference
	 * @param stop Exclusive ending position in the reference
	 * @param ranges Initially empty list receiving half-open pairs; adjacent ranges remain separate
	 * @param includeDels Whether to include deletion ranges
	 */
	public static void fillRanges(byte[] longmatch, int start, int stop, IntList ranges, boolean includeDels){
		assert(ranges.size==0);
		byte mode='?', lastMode='?';
		int rpos=start;
		int lastRpos=start;
		int rstart=start;
		for(int mpos=0; mpos<longmatch.length; mpos++){
			byte m=longmatch[mpos];
			if(m=='m' || m=='s' || m=='S' || m=='N' || m=='B'){//Little 's' is for a match classified as a sub to improve the affine score.
				mode='m';
				rpos++;
			}else if(m=='I' || m=='X' || m=='Y'){
				mode='I';
			}else if(m=='D'){
				mode='D';
				rpos++;
			}else if(m=='C'){
				mode='C';
				rpos++;
			}else{
				throw new RuntimeException("Invalid match string character '"+(char)m+"' = "+m+" (ascii).  " +
						"Match string should be in long format here.");
			}
			if(mode!=lastMode){
				if(mpos>0){
					if(lastMode=='m'){
						ranges.add(rstart);
						ranges.add(lastRpos);
					}else if(lastMode=='D' && includeDels){
						ranges.add(rstart);
						ranges.add(lastRpos);
					}
				}
				rstart=lastRpos;
			}
			lastMode=mode;
			lastRpos=rpos;
		}
		//Final cycle
		if(lastMode=='m'){
			ranges.add(rstart);
			ranges.add(lastRpos);
		}else if(mode=='D' && includeDels){
			ranges.add(rstart);
			ranges.add(lastRpos);
		}
		assert(rpos==stop) : start+", "+stop+", "+rpos+", "+rstart+", "+lastRpos+"\n"+new String(longmatch)+"\n"+ranges;
		assert((ranges.size&1)==0);
	}
	
	/** Handles cases where SAM reference name is not found in assembly.
	 * @param rname Reference name from SAM file */
	private void handleNoRef(String rname){
		assert(rname!=null);
		String ret=norefSet.putIfAbsent(rname, rname);
		if(ret==null){
			System.err.println("Warning! Scaffold not found in assembly: "+rname);
		}
	}
	
	/*--------------------------------------------------------------*/

	/**
	 * Loads references in input order, replacing each Read's obj with an owned BitSet.
	 * Seeds uppercase N positions, and all lowercase positions in lowercase mode,
	 * so later stage counts exclude those existing masks. Duplicate full IDs assert.
	 * Allocates coverage arrays only for SAM input with an enabled coverage threshold;
	 * their counts accumulate across SAM files before thresholds are applied.
	 * Updates reference counters before any splitting or base conversion. Owns and
	 * closes its input on the normal path; this does not add failure-cleanup guarantees.
	 * @return Insertion-ordered map of original IDs to loaded Reads with attached masks
	 */
	private LinkedHashMap<String, Read> hashRef(){


		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffinRef, null, qfinRef, null);
			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		
		final LinkedHashMap<String, Read> hmr=new LinkedHashMap<String, Read>();
		if(ffinSam!=null && (mincov>=0 || maxcov>=0)){
			covmap=new HashMap<String, CoverageArray>();
		}
		
		{

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

				for(int idx=0; idx<reads.size(); idx++){
					final Read r=reads.get(idx);
					final int len=r.length();
					final byte[] bases=r.bases;
					
					final BitSet bs=new BitSet(len);
					r.obj=bs;
					
					if(covmap!=null){
						if(bits32){
							covmap.put(r.id, new CoverageArray3(covmap.size(), len));
						}else{
							covmap.put(r.id, new CoverageArray2(covmap.size(), len));
						}
					}
					
					if(MaskByLowercase){
						for(int i=0; i<len; i++){
							if(bases[i]=='N' || Tools.isLowerCase(bases[i])){bs.set(i);}
						}
					}else{
						for(int i=0; i<len; i++){
							if(bases[i]=='N'){bs.set(i);}
						}
					}
					
					refReads++;
					refBases+=len;
					Read old=hmr.put(r.id, r);
					assert(old==null) : "Duplicate reference scaffold name "+r.id;
				}

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}

		errorState|=ReadWrite.closeStreams(cris);

		if(errorState){
			throw new RuntimeException("BBMask terminated in an error state; the output may be corrupt.");
		}
		
		return hmr;
	}
	
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Adds low-diversity marks across the reference map without changing bases.
	 * Reuses one matrix across reads; missing per-k arrays are allocated as needed.
	 * @param matrix Mutable per-k scratch arrays, or null to allocate a matrix
	 * @return Sum of newly set bits, excluding marks already present
	 */
	private long maskLowComplexity(short[][] matrix){
		long sum=0;
		if(matrix==null){matrix=new short[16][];}
		for(Read r : map.values()){
			sum+=maskLowComplexity(r, mink2, maxk2, window, ratio, matrix);
		}
		return sum;
	}
	
	/**
	 * Adds low-diversity marks to the Read's attached BitSet without changing bases.
	 * Examines the inclusive k range, reusing per-k counts and allocating missing arrays.
	 *
	 * @param r Read with bases and an attached BitSet to update
	 * @param mink Minimum k-mer size
	 * @param maxk Maximum k-mer size
	 * @param window Sliding window size
	 * @param ratio Minimum diversity ratio required
	 * @param matrix Nonnull mutable scratch matrix indexed by k; existing counts start at zero
	 * @return BitSet cardinality increase across all requested k values
	 */
	private static int maskLowComplexity(Read r, int mink, int maxk, int window, float ratio, short[][] matrix){
		
		final byte[] bases=r.bases;
		final BitSet bs=(BitSet)r.obj;
		
		int before=bs.cardinality();
		
		for(int k=mink; k<=maxk; k++){
			if(matrix[k]==null){matrix[k]=new short[(1<<(2*k))];}
		}
		
		for(int k=mink; k<=maxk; k++){
			final short[] counts=matrix[k];
			final int kmerspace=(1<<(2*k));
			final int mincount=(int)Math.ceil(ratio*Tools.min(window, kmerspace));
			maskLowComplexity(bases, bs, k, window, mincount, counts);
		}

		int after=bs.cardinality();
		
		return after-before;
	}
	
	
	/**
	 * Adds marks for full windows with insufficient k-mer diversity and no undefined bases.
	 * Reads bases without modification; counts is mutable rolling scratch expected to
	 * start at zero and drained on normal completion. It must not be shared by concurrent calls.
	 *
	 * @param bases Sequence bases to analyze
	 * @param bs BitSet to mark positions for masking
	 * @param k K-mer size for analysis
	 * @param window Size of sliding window
	 * @param mincount Minimum unique k-mer count required
	 * @param counts Zeroed frequency scratch indexed by encoded k-mer
	 */
	private static void maskLowComplexity(final byte[] bases, final BitSet bs, final int k, final int window, final int mincount, final short[] counts){
		assert(k>0) : "k must be greater than 0";
		
		if(verify){
			for(int c : counts){assert(c==0);}
		}
		
		final int mask=(k>15 ? -1 : ~((-1)<<(2*k)));
		int current=0, ns=0;
		int kmer=0, kmer2=0;
		
		for(int i=0, i2=-window; i2<bases.length; i++, i2++){
			
			if(i<bases.length){
				final byte b=bases[i];
				final int n=Dedupe.baseToNumber[b];
				kmer=((kmer<<2)|n)&mask;
				
				if(!AminoAcid.isFullyDefined(b)){ns++;}
				if(counts[kmer]<1){
					assert(counts[kmer]==0);
					current++;
				}
				counts[kmer]++;
				if(verify){assert(current==Tools.cardinality(counts)) : current+", "+Tools.cardinality(counts)+"\n"+Arrays.toString(counts);}
			}
			
			if(i2>=0){
				final byte b2=bases[i2];
				final int n2=Dedupe.baseToNumber[b2];
				kmer2=((kmer2<<2)|n2)&mask;

				if(!AminoAcid.isFullyDefined(b2)){
					ns--;
					assert(ns>=0);
				}
				counts[kmer2]--;
				if(counts[kmer2]<1){
					assert(counts[kmer2]==0) : Arrays.toString(counts);
					current--;
				}
				if(verify){assert(current==Tools.cardinality(counts)) : current+", "+Tools.cardinality(counts)+"\n"+Arrays.toString(counts);}
			}
			
			if(verify && i2>-1 && i<bases.length){
				assert(simd.Vector.sum(counts)==window);
			}
			
			if(current<mincount && ns<1 && i2>=-1 && i<bases.length){
				bs.set(i2+1, i+1);
			}
		}
		
	}
	
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	
	/** Adds entropy marks using a finite queue containing each reference once.
	 * Each worker owns its trackers and each polled Read's mask; bases are unchanged.
	 * Empty maps return zero. Counts are merged after all workers terminate.
	 * Reports observed owner interruption only after every worker has finished.
	 * @return Sum of newly set bits from the workers */
	private long maskLowEntropy(){
		if(map.isEmpty()){return 0;}
		ArrayBlockingQueue<Read> queue=new ArrayBlockingQueue<Read>(map.size());
		for(Read r : map.values()){queue.add(r);}
		int numThreads=Tools.min(Shared.threads(), queue.size());
		MaskLowEntropyThread[] threads=new MaskLowEntropyThread[numThreads];
		long sum=0;
		for(int i=0; i<threads.length; i++){threads[i]=new MaskLowEntropyThread(queue, mink2, maxk2, window, entropyCutoff);}
		for(int i=0; i<threads.length; i++){threads[i].start();}
		final InterruptedException interrupted=awaitWorkers(threads);
		for(int i=0; i<threads.length; i++){
			sum+=threads[i].masked;
		}
		if(interrupted!=null){
			Thread.currentThread().interrupt();
			throw new RuntimeException("BBMask interrupted while waiting for entropy workers.", interrupted);
		}
		return sum;
	}
	
	/** Worker thread for multi-threaded low-entropy masking.
	 * Each thread processes reads from shared queue using entropy analysis. */
	private class MaskLowEntropyThread extends Thread{
		
		/**
		 * Constructor for entropy masking thread.
		 *
		 * @param queue_ Queue containing reads to process
		 * @param mink_ Minimum k-mer size for entropy calculation
		 * @param maxk_ Maximum k-mer size for entropy calculation
		 * @param window_ Window size for entropy analysis
		 * @param cutoff_ Entropy threshold for masking
		 */
		MaskLowEntropyThread(ArrayBlockingQueue<Read> queue_, int mink_, int maxk_, int window_, float cutoff_){
			queue=queue_;
			minkT=mink_;
			maxkT=maxk_;
			//Preserve the requested window for EntropyTracker's short-counter bounds check (JGI-078).
			//Narrowing here would silently turn unsupported windows such as 65540 into 4.
			windowT=window_;
			cutoff=cutoff_;
			trackers=new EntropyTracker[maxkT+1];
			for(int i=minkT; i<=maxkT; i++){
				trackers[i]=new EntropyTracker(i, windowT, false, cutoff, true);
			}
		}
		
		@Override
		public void run(){
			for(Read r=queue.poll(); r!=null; r=queue.poll()){
				masked+=maskLowEntropy(r, minkT, maxkT, trackers);
			}
		}
		
		/** Queue containing reads to process */
		final ArrayBlockingQueue<Read> queue;
		/** Minimum k-mer size for this thread */
		final int minkT;
		/** Maximum k-mer size for this thread */
		final int maxkT;
		/** Entropy cutoff threshold for this thread */
		final float cutoff;
		/** Window size for entropy analysis */
		final int windowT;
		/** Array of entropy trackers for different k-mer sizes */
		final EntropyTracker[] trackers;
		/** Count of bases masked by this thread */
		long masked=0;
		
	}
	
	/**
	 * Adds low-entropy marks to the attached BitSet without changing Read bases.
	 * Reuses mutable trackers across the inclusive k range; they are not concurrent scratch.
	 *
	 * @param r Read with bases and an attached BitSet to update
	 * @param mink Minimum k-mer size for entropy calculation
	 * @param maxk Maximum k-mer size for entropy calculation
	 * @param trackers Nonnull mutable tracker array with entries for each requested k
	 * @return BitSet cardinality increase across all requested k values
	 */
	private int maskLowEntropy(Read r, int mink, int maxk, EntropyTracker[] trackers){
		final byte[] bases=r.bases;
		final BitSet bs=(BitSet)r.obj;
		
		int before=bs.cardinality();
		
		for(int k=mink; k<=maxk; k++){
			maskLowEntropy(bases, bs, trackers[k]);
		}

		int after=bs.cardinality();
		
		return after-before;
	}
	
	
	/**
	 * Adds marks for complete, undefined-base-free windows that fail the tracker cutoff.
	 * Leaves bases unchanged. A sequence shorter than the window returns without clearing
	 * the tracker; otherwise clears and reuses it, leaving its final window state on return.
	 * @param bases Sequence bases to analyze
	 * @param bs BitSet to mark positions for masking
	 * @param et Mutable tracker configured for k/window/cutoff, exclusively borrowed by this call
	 */
	public static void maskLowEntropy(final byte[] bases, final BitSet bs, final EntropyTracker et){
		final int window=et.windowBases();
		if(bases.length<window){return;}
		
		et.clear();
		for(int i=0, min=window-1; i<bases.length; i++){
			et.add(bases[i]);
			if(i>=min && et.ns()<1 && !et.passes()){bs.set(et.leftPos(), et.rightPos()+1);}
		}
	}
	
	
	/*--------------------------------------------------------------*/
	


	/** Adds repeat marks serially across the map without changing bases.
	 * @return Sum of newly set bits, excluding marks already present */
	private long maskRepeats_ST(){
		long sum=0;
		for(Read r : map.values()){
			sum+=maskRepeats(r, mink, maxk, mincount, minlen);
		}
		return sum;
	}
	
	/** Adds repeat marks using a finite queue containing each reference once.
	 * Each polled Read's mask is updated by one worker; bases are unchanged.
	 * Empty maps return zero. Counts are merged after all workers terminate.
	 * Reports observed owner interruption only after every worker has finished.
	 * @return Sum of newly set bits from the workers */
	private long maskRepeats(){
		if(map.isEmpty()){return 0;}
		ArrayBlockingQueue<Read> queue=new ArrayBlockingQueue<Read>(map.size());
		for(Read r : map.values()){queue.add(r);}
		int numThreads=Tools.min(Shared.threads(), queue.size());
		MaskRepeatThread[] threads=new MaskRepeatThread[numThreads];
		long sum=0;
		for(int i=0; i<threads.length; i++){threads[i]=new MaskRepeatThread(queue, mink, maxk, mincount, minlen);}
		for(int i=0; i<threads.length; i++){threads[i].start();}
		final InterruptedException interrupted=awaitWorkers(threads);
		for(int i=0; i<threads.length; i++){
			sum+=threads[i].masked;
		}
		if(interrupted!=null){
			Thread.currentThread().interrupt();
			throw new RuntimeException("BBMask interrupted while waiting for repeat workers.", interrupted);
		}
		return sum;
	}
	
	/** Worker thread for multi-threaded repeat masking.
	 * Each thread processes reads from shared queue using k-mer repeat detection. */
	private class MaskRepeatThread extends Thread{
		
		/**
		 * Constructor for repeat masking thread.
		 *
		 * @param queue_ Queue containing reads to process
		 * @param mink_ Minimum k-mer size for repeat detection
		 * @param maxk_ Maximum k-mer size for repeat detection
		 * @param mincount_ Minimum count threshold for repeats
		 * @param minlen_ Minimum length of repeats to mask
		 */
		MaskRepeatThread(ArrayBlockingQueue<Read> queue_, int mink_, int maxk_, int mincount_, int minlen_){
			queue=queue_;
			minkT=mink_;
			maxkT=maxk_;
			mincountT=mincount_;
			minlenT=minlen_;
		}
		
		@Override
		public void run(){
			for(Read r=queue.poll(); r!=null; r=queue.poll()){
				masked+=maskRepeats(r, minkT, maxkT, mincountT, minlenT);
			}
		}
		
		/** Queue containing reads to process */
		final ArrayBlockingQueue<Read> queue;
		/** Minimum k-mer size for this thread */
		final int minkT;
		/** Maximum k-mer size for this thread */
		final int maxkT;
		/** Minimum count threshold for this thread */
		final int mincountT;
		/** Minimum repeat length for this thread */
		final int minlenT;
		/** Count of bases masked by this thread */
		long masked=0;
		
	}
	
	/**
	 * Adds repeat marks to the attached BitSet without changing Read bases.
	 * For each k in the inclusive range, requires at least max(minlen, k*mincount) bases.
	 *
	 * @param r Read with bases and an attached BitSet to update
	 * @param mink Minimum k-mer size for repeat detection
	 * @param maxk Maximum k-mer size for repeat detection
	 * @param mincount Minimum repeat count threshold
	 * @param minlen Minimum length of repeats to mask
	 * @return BitSet cardinality increase across all requested k values
	 */
	private static int maskRepeats(Read r, int mink, int maxk, int mincount, int minlen){
		final byte[] bases=r.bases;
		final BitSet bs=(BitSet)r.obj;
		
		int before=bs.cardinality();
		
		for(int k=mink; k<=maxk; k++){
			maskRepeats(bases, bs, k, Tools.max(minlen, k*mincount));
		}

		int after=bs.cardinality();
		
		return after-before;
	}
	
	
	/**
	 * Adds qualifying tandem-repeat spans to the BitSet without modifying bases.
	 * A candidate at loc uses the preceding seed [loc-k,loc); marked spans include it.
	 *
	 * @param bases Sequence bases to analyze
	 * @param bs BitSet to mark positions for masking
	 * @param k K-mer size for repeat detection
	 * @param minlen Minimum length of repeats to mask
	 */
	public static void maskRepeats(final byte[] bases, final BitSet bs, final int k, final int minlen){
		final int lim=bases.length-k;
		final int mask=(k>15 ? -1 : ~((-1)<<(2*k)));
		for(int loc=0; loc<lim; loc++){
			int len=repeatLength(bases, k, mask, loc);
			if(len>=minlen){
				int a=loc-k, b=loc-k+len;
				bs.set(a, b);
				loc=Tools.max(loc, b-minlen);
			}else{
			}
		}
		
	}
	
	
	/**
	 * Extends the seed [loc-k,loc) while its encoded k-mer recurs without a k-base gap.
	 * The returned length includes the seed and ends after the last matching k-mer.
	 *
	 * @param bases Sequence bases to analyze
	 * @param k K-mer size for pattern matching
	 * @param mask Bit mask for k-mer extraction
	 * @param loc Exclusive end of the seed and first extension position
	 * @return Seed-inclusive span length, or zero if no matching extension is found
	 */
	private static int repeatLength(final byte[] bases, final int k, final int mask, final int loc){
		
		final int lim=bases.length;
		final int key=getInitialKey(bases, loc, k);
		if(key<0){return 0;}
		int kmer=key;
		int gap=0, last=-1;
		for(int i=loc; i<lim && gap<k; i++){
			final byte b=bases[i];
			final int n=Dedupe.baseToNumber[b];
			kmer=((kmer<<2)|n)&mask;
			if(kmer==key){
				last=i;
				gap=0;
			}else{
				gap++;
			}
		}
		
		return (last<0 ? 0 : last-loc+k+1);
	}
	
	/**
	 * Encodes the preceding k bases [loc-k,loc) as a repeat seed.
	 *
	 * @param bases Sequence bases
	 * @param loc Exclusive end of the seed in bases
	 * @param k K-mer size
	 * @return Encoded seed, or -1 when loc is too small to contain k preceding bases
	 */
	private static int getInitialKey(byte[] bases, int loc, int k){
		assert(k<16);
		int start=loc-k;
		int key=0;
		if(start<0){return -1;}
		for(int i=start; i<loc; i++){
			final byte b=bases[i];
			final int n=Dedupe.baseToNumber[b];
			key=(key<<2)|n;
		}
		assert(key>=0);
		return key;
	}

	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/

	/** Map storing reference sequences with associated masking BitSets */
	private LinkedHashMap<String, Read> map=null;
	/** Set tracking reference names not found in assembly */
	private ConcurrentHashMap<String, String> norefSet=new ConcurrentHashMap<String, String>(256, .75f, 16);
	/** Coverage arrays for SAM-based masking analysis */
	private HashMap<String, CoverageArray> covmap=null;
	
	/** Loaded reference record count, before splitting changes output record count. */
	private long refReads=0;
	/** Loaded reference bases, before masking or splitting. */
	private long refBases=0;
	
	/** Count of SAM reads processed */
	private long samReads=0;
	/** Total bases in SAM reads processed */
	private long samBases=0;
	
	/** Flag indicating whether processing encountered errors */
	public boolean errorState=false;

	/** Input reference file path */
	private String inRef=null;
	/** List of input SAM file paths */
	private ArrayList<String> inSam=new ArrayList<String>();

	/** Quality file for reference input */
	private String qfinRef=null;

	/** Output reference file path */
	private String outRef=null;

	/** Quality file for reference output */
	private String qfoutRef=null;

	/** File extension override for reference input */
	private String extinRef=null;
	/** File extension override for reference output */
	private String extoutRef=null;
	
	/** Whether to overwrite existing output files */
	private boolean overwrite=true;
	/** Whether to append to existing output files */
	private boolean append=false;

	/** Input limit passed separately to the reference stream and each SAM stream. */
	private long maxReads=-1;

	/** Whether to perform repeat masking */
	private boolean processRepeats=false;
	/** Minimum k-mer size for repeat detection */
	private int mink=5;
	/** Maximum k-mer size for repeat detection */
	private int maxk=5;
	/** Minimum length of repeats to mask */
	private int minlen=40;
	/** Minimum count threshold for repeat detection */
	private int mincount=4;

	/** Whether to perform entropy/complexity masking */
	private boolean processEntropy=true;
	/** Whether to use entropy-based (true) or complexity-based (false) masking */
	private boolean entropyMode=true;
	/** Whether to split reads at masked regions instead of masking */
	private boolean splitMode=false;
	/** Minimum k-mer size for entropy/complexity analysis */
	private int mink2=5;
	/** Maximum k-mer size for entropy/complexity analysis */
	private int maxk2=5;
	/** Window size for sliding window entropy/complexity analysis */
	private int window=80;
	/** Minimum k-mer diversity ratio for complexity masking */
	private float ratio=0.35f; //For complexity, if not in entropyMode
	/** Entropy threshold below which regions are masked */
	private float entropyCutoff=0.70f;
	
	/** Use 32-bit coverage arrays */
	private boolean bits32=true;
	/** Include deletions when calculating coverage */
	private boolean includeDeletionCoverage=true;

	/** If nonnegative, mask bases with coverage strictly below this minimum. */
	private int mincov=-1;
	/** If nonnegative, mask bases with coverage strictly above this maximum. */
	private int maxcov=-1;
	
	/** Padding to add around SAM alignment regions */
	private int samPad=0;

	/** File format specification for reference input */
	private final FileFormat ffinRef;
	/** Array of file format specifications for SAM inputs */
	private final FileFormat[] ffinSam;

	/** File format specification for reference output */
	private final FileFormat ffoutRef;

	/** Output stream for logging and status messages */
	private PrintStream outstream=System.err;

	/*--------------------------------------------------------------*/
	
	/** Global flag for verbose output */
	public static boolean verbose=false;
	/** Whether to convert non-ACGTN characters to N during masking */
	public static boolean CONVERT_NON_ACGTN=true;
	/** Debug flag for verification checks */
	private static boolean verify=false;
	/** Whether to mask by converting to lowercase instead of N */
	private static boolean MaskByLowercase=false;

}
