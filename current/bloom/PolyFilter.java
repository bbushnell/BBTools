package bloom;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;

import cardinality.CardinalityTracker;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBDuk;
import jgi.BBMerge;
import shared.KillSwitch;
import shared.MetadataWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import structures.ListNum;
import structures.LongHashSet;
import structures.Quantizer;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.EntropyTracker;
import tracker.ReadStats;

/**
 * Filters reads with artificial homopolymers.
 * 
 * @author Brian Bushnell
 * @date August 21, 2024
 *
 */
public class PolyFilter implements Accumulator<PolyFilter.ProcessThread> {
	
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
		PolyFilter x=new PolyFilter(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public PolyFilter(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean setInterleaved=false; //Whether interleaved was explicitly set.
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Tools.max(Shared.threads()>1 ? 2 : 1, Shared.threads()>20 ? Shared.threads()/2 : Shared.threads()));
		BBMerge.strict=true;
		KCountArray.LOCKED_INCREMENT=true;
		
		//Create a parser object
		Parser parser=new Parser();
		parser.loglog=parser.loglogOut=true;
		
		boolean setBits=false;
		int k_=31;
		int ksmall_=-1;
		int hashes_=3;
		int bits_=2;
		boolean rcomp_=true;
		boolean merge_=true;
		boolean tossjunk_=false;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("ordered")){
				ordered=Parse.parseBoolean(b);
			}
			
			else if(a.equalsIgnoreCase("k") || a.equalsIgnoreCase("bloomK") || a.equalsIgnoreCase("kbloom") || a.equalsIgnoreCase("bloomFilterK") || a.equalsIgnoreCase("kbig")){
				k_=Integer.parseInt(b);
				assert(k_>0 && k_<=31) : "K must be between 1 and 31, inclusive.";
			}else if(a.equalsIgnoreCase("ksmall") || a.equalsIgnoreCase("bloomKsmall") || a.equalsIgnoreCase("bloomFilterKsmall")){
				ksmall_=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("hashes") || a.equalsIgnoreCase("bloomHashes") || a.equalsIgnoreCase("bloomFilterHashes")){
				hashes_=Integer.parseInt(b);
			}else if(a.equals("rcomp")){
				rcomp_=Parse.parseBoolean(b);
			}else if(a.equals("bits")){
				setBits=true;
				bits_=Integer.parseInt(b);
				assert(bits_>0);
			}
			
			else if(a.equals("kpoly") || a.equals("kset") || a.equals("polyk")) {
				kpoly=Integer.parseInt(b);
			}else if(a.equals("hdist")) {
				hdist=Integer.parseInt(b);
			}else if(a.equals("mm") || a.equals("maskmiddle")) {
				maskMiddle=Parse.parseBoolean(b);
			}
			
			else if(a.equals("mincount") || a.equals("mindepth")){
				minCount=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("minPolymer") || a.equalsIgnoreCase("minPoly")){
				minPolymer=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("minPolymer2") || a.equalsIgnoreCase("minPoly2")){
				minPolymer2=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("lowCountFraction") || a.equals("lcf") || 
					a.equalsIgnoreCase("lowDepthFraction") || a.equals("ldf")){
				lowCountFraction=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("lowCountFraction2") || a.equals("lcf2") || 
					a.equalsIgnoreCase("lowDepthFraction2") || a.equals("ldf2")){
				lowCountFraction2=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("entropyCutoff") || a.equals("entropy") || a.equals("minentropy")){
				entropyCutoff=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("entropyCutoff2") || a.equals("entropy2") || a.equals("minentropy2")){
				entropyCutoff2=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("qualityCutoff") || a.equals("maq") || 
					a.equals("quality") || a.equals("minquality")){
				qualityCutoff=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("qualityCutoff2") || a.equals("maq2") ||
					a.equals("quality2") || a.equals("minquality2")){
				qualityCutoff2=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("purity")){
				nonPolyFraction=1-Float.parseFloat(b);
				assert(nonPolyFraction>=0 && nonPolyFraction<1) : b;
			}else if(a.equals("mers") || a.equals("polymers") || a.equals("polymer") || a.equals("bases")){
				polymers=b.toUpperCase().getBytes();
			}else if(a.equals("trimpolymers") || a.equals("trimpoly")){
				trimPolymers=(b==null ? "" : b.toUpperCase()).getBytes();
			}
			
			else if(a.equals("minprob")){
				ReadCounter.minProb=Float.parseFloat(b);
			}else if(a.equals("maxload")){
				maxLoad=Float.parseFloat(b);
			}else if(a.equals("merge")){
				merge_=Parse.parseBoolean(b);
			}else if(a.equals("tossjunk")){
				tossjunk_=Parse.parseBoolean(b);
			}else if(a.equals("memfraction") || a.equals("memmult") || a.equals("memratio")){
				memFraction=Float.parseFloat(b);
			}
			
			else if(a.equals("trim")){
				trimLeft=trimRight=parseIntOrBool(b, 6, 0);
			}else if(a.equals("trimleft")){
				trimLeft=parseIntOrBool(b, 6, 0);
			}else if(a.equals("trimright")){
				trimLeft=parseIntOrBool(b, 6, 0);
			}else if(a.equals("minlen") || a.equals("minlength")){
				minLength=Integer.parseInt(b);
			}else if(a.equals("maxnonpoly")){
				maxNonPoly=Integer.parseInt(b);
			}else if(a.equals("quantize") || a.equals("quantizesticky")){
				quantizeQuality=Quantizer.parse(arg, a, b);
			}
			
			else if(a.equals("cells")){
				BloomFilter.OVERRIDE_CELLS=Parse.parseKMG(b);
			}else if(a.equals("seed")){
				KCountArray7MTA.setSeed(Parse.parseKMG(b));
			}else if(a.equals("smoothwidth") || a.equals("junkwidth") || a.equals("width")){
				junkWidth=Integer.parseInt(b);
			}
			
			else if(a.equals("ref")){
				addFiles(b, ref);
			}else if(a.equals("extra")){
				addFiles(b, extra);
			}else if(a.equals("outm") || a.equals("outm1") || a.equals("out") || a.equals("out1")){
				out1=b;
			}else if(a.equals("outm2") || a.equals("out2")){
				out2=b;
			}else if(a.equals("outb") || a.equals("outb1") || a.equals("outbad") || a.equals("outbad1") || a.equals("outlow") || a.equals("outlow1")){
				outbad1=b;
			}else if(a.equals("outb2") || a.equals("outbad2") || a.equals("outlow2")){
				outbad2=b;
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		while(minCount>0 && ((1L<<bits_)-1)<minCount){bits_*=2;}
		
		if(ksmall_<=0){ksmall_=k_;}
		assert(ksmall_<=k_) : k_+", "+ksmall_;
		
		kbloom=k_;
		ksmall=Tools.min(kbloom, ksmall_);
		bits=bits_;
		hashes=hashes_;
		rcomp=rcomp_;
		tossjunk=tossjunk_;
		merge=merge_;
		KmerCountAbstract.CANONICAL=rcomp; //rcomp, or true, perhaps...  hmmm.
		polymerFlagged=new long[polymers.length];
		if(trimPolymers==null){trimPolymers=polymers;}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			
			extin=parser.extin;
			extout=parser.extout;
			
//			parser.loglogk=k;
			parser.loglogMinprob=ReadCounter.minProb;
			loglogOut=((parser.loglog&parser.loglogOut) ? CardinalityTracker.makeTracker(parser) : null);
		}
		
		SamLine.SET_FROM_OK=true;
		
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		
		//Do output file # replacement
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
		ref=Tools.fixExtension(ref);
		extra=Tools.fixExtension(extra);
		
		//Do output file # replacement
		if(outbad1!=null && outbad2==null && outbad1.indexOf('#')>-1){
			outbad2=outbad1.replace("#", "2");
			outbad1=outbad1.replace("#", "1");
		}
		
		//Adjust interleaved detection based on the number of input files
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		//Ensure out2 is not set without out1
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
		
		//Adjust interleaved settings based on number of output files
		if(!setInterleaved){
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outbad1, outbad2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2, outbad1, outbad2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		
		//Create output FileFormat objects
		ffoutm1=FileFormat.testOutput(outbad1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffoutm2=FileFormat.testOutput(outbad2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		if(kpoly>0) {
			kpoly=Tools.min(kpoly, 32);
			midMask=makeMidMask(maskMiddle, kpoly);
			kmerMask=(kpoly==32 ? 0 : ~((-1L)<<(2*kpoly)));
			set=makeSet(polymers, kpoly, hdist, midMask);
			System.err.println("Added "+set.size()+" kmer mutants to filter set.");
		}
		
		if(lowCountFraction>1.0f){
			filter=null;
		}else{
			outstream.println("Using "+bits+" bits per cell.");
			{
				Timer t=new Timer(outstream, true);
				if(ref.isEmpty()){
					filter=new BloomFilter(in1, in2, extra, ksmall, kbloom, bits, hashes, 1,
							rcomp, false, false, memFraction);
				}else{
					ref.addAll(extra);
					filter=new BloomFilter(null, null, ref, ksmall, kbloom, bits, hashes, 1,
							rcomp, false, false, memFraction);
				}
				t.stop("Filter creation: \t\t");
				outstream.println(filter.filter.toShortString());
			}

			{
				double a=filter.filter.estimateUniqueKmers(hashes);
				outstream.println("Estimated kmers of depth 1+: \t"+(long)a);
				if(bits>1){
					double usedFraction2=filter.filter.usedFraction(2);
					double b=filter.filter.estimateUniqueKmersFromUsedFraction(hashes, usedFraction2);
					outstream.println("Estimated kmers of depth 2+: \t"+(long)b);
					outstream.println("Used fraction for depth 2+:  \t"+Tools.format("%.3f%%", usedFraction2*100));
					if(usedFraction2>maxLoad){
						KillSwitch.kill("Max load exceeded, quitting: "+usedFraction2+" > "+maxLoad);
					}
				}else if(maxLoad<1.0f){
					double usedFraction=filter.filter.usedFraction(2);
					if(usedFraction>maxLoad){
						KillSwitch.kill("Max load exceeded, quitting: "+usedFraction+" > "+maxLoad);
					}
				}
			}
		}
	}
	
	private static void addFiles(String b, ArrayList<String> list){
		if(b==null){list.clear();}
		else{
			if(new File(b).exists()){list.add(b);}
			else{
				for(String s : b.split(",")){list.add(s);}
			}
		}
	}
	
	private static final int parseIntOrBool(String b, int defaultTrue, int defaultFalse) {
		if(b==null || Tools.startsWithLetter(b)) {
			boolean x=Parse.parseBoolean(b);
			return (x ? defaultTrue : defaultFalse);
		}
		return Integer.parseInt(b);
	}
	
	public static long makeMidMask(boolean maskMiddle, int k) {
		if(!maskMiddle) {return -1L;}
		int bitsPerBase=2;
		int midMaskLen=2-(k&1);
		assert(k>midMaskLen+1);
		int bits=midMaskLen*bitsPerBase;
		int shift=((k-midMaskLen)/2)*bitsPerBase;
		long middleMask=~((~((-1L)<<bits))<<shift);
		return middleMask;
	}
	
	public static LongHashSet makeSet(byte[] polymers, int k, int hdist, long midMask) {
		if(polymers==null || polymers.length<1) {return null;}
		LongHashSet set=new LongHashSet(1024);
		for(byte b : polymers) {
			long x=AminoAcid.baseToNumber[b];
			long kmer=0;
			for(int i=0; i<k; i++) {kmer=((kmer<<2)|x);}
			addToSet(kmer, k, midMask, hdist, set);
		}
		return set;
	}
	
	public static void addToSet(final long kmer, final long k, final long midMask, final int hdist,
			final LongHashSet set){
		set.add(kmer&midMask);
		if(hdist<1) {return;}
		
		for(int i=0; i<k; i++) {
			final int shift=2*i;
			final long masked=kmer&~(3L<<shift);
			for(long x=0; x<4; x++) {
				final long mutant=(masked|(x<<shift))&midMask;
				if(mutant!=kmer) {addToSet(mutant, k, midMask, hdist-1, set);}
			}
		}
	}
	
	public static boolean kmerScan(Read r, int k, long kmerMask, long midMask, LongHashSet set) {
		if(r==null || r.length()<k) {return false;}
		final byte[] bases=r.bases;
		long kmer=0;
		int len=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			kmer=((kmer<<2)|x)&kmerMask;
			if(x<0){len=0;}else{len++;}
			if(len>=k && set.contains(kmer&midMask)){return true;}
		}
		return false;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	public void process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		//Optionally create read output streams
		final ConcurrentReadOutputStream ros, rosb;
		if(ffout1!=null){
			//Select output buffer size based on whether it needs to be ordered
			final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);
			
			//Notify user of output mode
			if(cris.paired() && out2==null && (in1!=null && !ffin1.samOrBam() && !ffout1.samOrBam())){
				outstream.println("Writing interleaved.");
			}
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, null, null, buff, null, false);
			ros.start(); //Start the stream
		}else{ros=null;}
		
		if(ffoutm1!=null){
			//Select output buffer size based on whether it needs to be ordered
			final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);
			
			rosb=ConcurrentReadOutputStream.getStream(ffoutm1, ffoutm2, null, null, buff, null, false);
			rosb.start(); //Start the stream
		}else{rosb=null;}
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		Timer t2=new Timer(outstream, true);
		
		//Process the reads in separate threads
		spawnThreads(cris, ros, rosb);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, ros, rosb);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		t2.stop("\nFiltering Time:  \t\t");
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, true));
		if(merge) {
			outstream.println(Tools.numberPercent("Merged Reads:", readsMerged, readsMerged*100.0/readsProcessed, 4, 8));
		}
		if(polymers!=null && polymers.length>0) {
			for(int i=0; i<polymers.length; i++) {
				final String b=Character.toString((char)polymers[i]);
				final long x=polymerFlagged[i];
				outstream.println(Tools.numberPercent("Poly-"+b+" Reads:", x, x*100.0/readsProcessed, 4, 8));
			}
		}
		if(entropyCutoff2>0) {
			outstream.println(Tools.numberPercent("Low Entropy Reads:", entropyFlagged, entropyFlagged*100.0/readsProcessed, 4, 8));
		}
		if(qualityCutoff2>0) {
			outstream.println(Tools.numberPercent("Low Quality Reads:", qualityFlagged, qualityFlagged*100.0/readsProcessed, 4, 8));
		}
		if(filter!=null && lowCountFraction2<=1) {
			outstream.println(Tools.numberPercent("Low Depth Reads:", depthFlagged, depthFlagged*100.0/readsProcessed, 4, 8));
		}
		if(trimLeft>0 || trimRight>0) {
			outstream.println(Tools.numberPercent("Trimmed Reads:", readsTrimmed, readsTrimmed*100.0/readsProcessed, 4, 8));
		}
		{
			lastReadsOut=readsOut;
			lastReadsRemoved=readsProcessed-readsOut;
			lastBasesOut=basesOut;
			lastBasesRemoved=basesProcessed-basesOut;
			
			long x=readsProcessed-readsOut;
			outstream.println(Tools.numberPercent("Total Discards:", x, x*100.0/readsProcessed, 4, 8));
		}
		if(loglogOut!=null){
			outstream.println("Unique "+loglogOut.k+"-mers out:     \t"+loglogOut.cardinality());
		}
		
		MetadataWriter.write(null, readsProcessed, basesProcessed, readsOut, basesOut, false);
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros, final ConcurrentReadOutputStream rosb){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, ros, rosb, i));
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		
//		//Do anything necessary after processing
		
	}
	
	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt) {
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			readsMerged+=pt.readsMergedT;
			readsTrimmed+=pt.readsTrimmedT;
			basesTrimmed+=pt.basesTrimmedT;
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			errorState|=(!pt.success);

			Tools.add(polymerFlagged, pt.polymerFlaggedT);
			entropyFlagged+=pt.entropyFlaggedT;
			qualityFlagged+=pt.qualityFlaggedT;
			depthFlagged+=pt.depthFlaggedT;
		}
	}
	
	@Override
	public final boolean success(){return !errorState;}

	@Override
	public ReadWriteLock rwlock() {return null;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Calculate the length of the longest homopolymer region in this sequence,
	 * allowing a specified error rate, using the Phred algorithm.
	 * There may be longer regions, but they will have flanking subregions with
	 * greater than the specified error rate.
	 * @param bases Sequence to examine.
	 * @param mer Subunit of the homopolymer to count (e.g. 'G').
	 * @param nonFraction Maximal allowed fraction of nonmatching bases in the region.
	 * @return Length of the longest homopolymer.
	 */
	public static int polymerLen_old(final byte[] bases, final byte mer, final float nonFraction) {
		if(nonFraction<=0) {return polymerLen(bases, mer);}
		
		float score=0;
		float maxScore=0;
		int maxLen=0;
		int lastZero=-1;
		
		for(int i=0; i<bases.length; i++) {
			final byte b=bases[i];
			if(b!=mer) {//common case
				score=Tools.max(0, score-1);
				if(score<=0) {
					lastZero=i;
					maxScore=0;
				}
			}else {
				score+=nonFraction;
				if(score>=maxScore) {
					maxScore=score;
					maxLen=Tools.max(maxLen, i-lastZero);
				}
			}
		}
		
		return maxLen;
	}

	//This version does not look for a peak, just a maximal distance above 0.
	//It is not bidirectionally symmetric though, until the second loop.
	public static int polymerLen(final byte[] bases, final byte mer, final float nonFraction) {
		if(nonFraction<=0) {return polymerLen(bases, mer);}
		
		float score=0;
		int maxLen=0;
		int lastZero=-1;
		
		for(int i=0; i<bases.length; i++) {
			final byte b=bases[i];
			if(b!=mer) {//common case
				score=Tools.max(0, score-1);
				if(score<=0) {
					lastZero=i;
				}
			}else {
				score+=nonFraction;
				maxLen=Tools.max(maxLen, i-lastZero);
			}
		}
		
		score=0;
		lastZero=-1;
		for(int i=0; i<bases.length; i++) {
			final byte b=bases[bases.length-i-1];
			if(b!=mer) {//common case
				score=Tools.max(0, score-1);
				if(score<=0) {
					lastZero=i;
				}
			}else {
				score+=nonFraction;
				maxLen=Tools.max(maxLen, i-lastZero);
			}
		}
		
		return maxLen;
	}
	
	/**
	 * Calculate the length of the longest homopolymer in this sequence.
	 * @param bases Sequence to examine.
	 * @param mer Subunit of the homopolymer to count (e.g. 'G').
	 * @return Length of the longest homopolymer.
	 */
	public static int polymerLen(final byte[] bases, final byte mer) {
		int max=0, current=0;
		for(byte b : bases) {
			if(b==mer) {current++;}
			else {
				max=Tools.max(max, current);
				current=0;
			}
		}
		return Tools.max(current, max);
	}
	
	public int trimRead(Read r) {
		if(r==null) {return 0;}
		int trimmed=0;
		for(int i=0; i<trimPolymers.length; i++) {
			final byte b=trimPolymers[i];
			int x=BBDuk.trimPoly(r, trimLeft, trimRight, maxNonPoly, b);
			trimmed+=x;
		}
		return trimmed;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final ConcurrentReadOutputStream ros_, final ConcurrentReadOutputStream rosb_, final int tid_){
			cris=cris_;
			ros=ros_;
			rosb=rosb_;
			tid=tid_;
			
			eTracker=new EntropyTracker(false, Tools.max(0, entropyCutoff), entropyHighpass);
		}
		
		//Called by start()
		@Override
		public void run(){
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Iterate through the reads */
		void processInner(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}

			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
				ArrayList<Read> keepList=new ArrayList<Read>(reads.size());
				ArrayList<Read> tossList=new ArrayList<Read>(reads.size());

				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					Read r1=reads.get(idx);
					Read r2=r1.mate;
					
					//Validate reads in worker threads
					if(!r1.validated()){r1.validate(true);}
					if(r2!=null && !r2.validated()){r2.validate(true);}

					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=r1.mateLength();

					//Increment counters
					readsProcessedT+=r1.pairCount();
					basesProcessedT+=initialLength1+initialLength2;

					final int insert=(r2==null || !merge ? -1 : 
						BBMerge.findOverlapVStrict(r1, r2, false));
					final boolean junk1=isJunk(r1, insert);
					final boolean junk2=isJunk(r2, insert);
					boolean keep=!junk1 && !junk2;
					if(keep && (trimLeft>0 || trimRight>0)) {
						int trimmed1=trimRead(r1);
						int trimmed2=trimRead(r2);
						keep=r1.length()>minLength && (r2==null || r2.length()>=minLength);
						readsTrimmedT+=(trimmed1>0 ? 1 : 0)+(trimmed2>0 ? 1 : 0);
						basesTrimmedT+=trimmed1+trimmed2;
					}
					
					
					readsMergedT+=(insert>0 ? 2 : 0);
					assert(r1.mate==r2);
					
					if(keep){
						if(quantizeQuality) {
							Quantizer.quantize(r1);
							Quantizer.quantize(r2);
						}
						if(loglogOut!=null){loglogOut.hash(r1);}
						readsOutT+=r1.pairCount();
						basesOutT+=r1.pairLength();
						keepList.add(r1);
					}else{
						tossList.add(r1);
					}
				}

				//Output reads to the output streams
				if(ros!=null){ros.add(keepList, ln.id);}
				if(rosb!=null){rosb.add(tossList, ln.id);}

				//Notify the input stream that the list was used
				cris.returnList(ln);

				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		private boolean isJunk(Read r, int insert) {
			if(r==null) {return false;}
			
			//Reads merged with a short insert; probably OK regardless of poly-G
			if(insert>0 && insert<r.length()+10) {
				return false;
			}
			
			final float entropy=(eTracker==null ? 1 : eTracker.averageEntropy(r.bases, true));
			final float quality=(float)r.avgQualityByProbabilityDouble(false, r.length());
			final float lcf=(filter==null ? 0 : filter.lowCountFraction(r, minCount, true));
			
			int worstBase=0, worstLen=0;
			for(int i=0; i<polymers.length; i++) {
				final byte b=polymers[i];
				int len=polymerLen(r.bases, b, nonPolyFraction);
				if(len>worstLen) {
					worstLen=len;
					worstBase=i;
				}
			}
			
			if(worstLen>minPolymer2) {
				polymerFlaggedT[worstBase]++;
				return true;
			}
			
			if(set!=null && kmerScan(r, kpoly, kmerMask, midMask, set)) {
				polymerFlaggedT[worstBase]++;
				return true;
			}
			
			//Reads merged with a long insert; increase stringency
			if(insert>0 && insert>=r.length()+10) {
				return false;
			}
			
			if(worstLen>minPolymer && 
					(entropy<entropyCutoff || lcf>lowCountFraction || quality<qualityCutoff)) {
				polymerFlaggedT[worstBase]++;
				return true;
			}
			if(entropy<entropyCutoff2) {
				entropyFlaggedT++;
				return true;
			}
			if(quality<qualityCutoff2) {
				qualityFlaggedT++;
				return true;
			}
			if(lcf>=lowCountFraction2) {
				depthFlaggedT++;
				return true;
			}
			return false;
		}
		
		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		protected long readsTrimmedT=0;
		protected long basesTrimmedT=0;
		protected long readsMergedT=0;
		
		/** Number of reads retained by this thread */
		protected long readsOutT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Shared output stream */
		private final ConcurrentReadOutputStream ros;
		/** Matched output stream */
		private final ConcurrentReadOutputStream rosb;
		/** Thread ID */
		final int tid;
		
		final EntropyTracker eTracker;
		
		final long[] polymerFlaggedT=new long[polymers.length];
		long entropyFlaggedT=0;
		long qualityFlaggedT=0;
		long depthFlaggedT=0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private ArrayList<String> ref=new ArrayList<String>();
	private ArrayList<String> extra=new ArrayList<String>();
	
	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;

	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path */
	private String out2=null;

	/** Output file path for bad reads */
	private String outbad1=null;
	/** Secondary output file path for bad reads */
	private String outbad2=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	/** For calculating kmer cardinality in output */
	final CardinalityTracker loglogOut;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	public long readsProcessed=0;
	/** Number of bases processed */
	public long basesProcessed=0;
	/** Number of reads trimmed */
	public long readsTrimmed=0;
	/** Number of bases trimmed */
	public long basesTrimmed=0;
	protected long readsMerged=0;

	/** Number of reads retained */
	public long readsOut=0;
	/** Number of bases retained */
	public long basesOut=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;

	public static long lastReadsOut=-1;
	public static long lastReadsRemoved=-1;
	public static long lastBasesOut=-1;
	public static long lastBasesRemoved=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	final FileFormat ffin1;
	/** Secondary input file */
	final FileFormat ffin2;
	
	/** Primary output file */
	private final FileFormat ffout1;
	/** Secondary output file */
	private final FileFormat ffout2;
	
	/** Primary output file for matching reads */
	private final FileFormat ffoutm1;
	/** Secondary output file for matching reads */
	private final FileFormat ffoutm2;
	
	final BloomFilter filter;

	final int kbloom;
	final int ksmall;
	final int hashes;
	final int bits;
	final boolean rcomp;
	final boolean tossjunk;
	final boolean merge;
	int junkWidth=1;
	float memFraction=1.0f;
	float maxLoad=1.0f; //Crash if this is exceeded

	int trimLeft=6;
	int trimRight=6;
	int minLength=50;
	int maxNonPoly=2;
	
	int kpoly=29;
	int hdist=2;
	long midMask=-1;
	long kmerMask=-1;
	boolean maskMiddle=true;
	LongHashSet set=null;
	
	int minCount=2;
	float lowCountFraction=0.24f;
	float lowCountFraction2=1.1f;
	int minPolymer=20;
	int minPolymer2=29;
	float entropyCutoff=0.67f;
	float entropyCutoff2=0.20f;
	float qualityCutoff=12.5f;
	float qualityCutoff2=7.5f;
	float nonPolyFraction=1-0.85f;//0.15f;//These numbers are slightly different.
	byte[] polymers="GC".getBytes();
	byte[] trimPolymers=null;
	
	long[] polymerFlagged;
	long entropyFlagged;
	long qualityFlagged;
	long depthFlagged;
	
	boolean entropyHighpass=true;

	/** Quantize quality scores to reduce file size */
	boolean quantizeQuality=false;
	
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
	/** Reads are output in input order */
	private boolean ordered=false;
	
}
