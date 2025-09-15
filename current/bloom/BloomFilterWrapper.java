package bloom;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.BBMerge;
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
import structures.LongList;
import tracker.ReadStats;

/**
 * Wraps a BloomFilter to filter reads.
 * 
 * @author Brian Bushnell
 * @date April 23, 2018
 *
 */
public class BloomFilterWrapper {
	
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
		BloomFilterWrapper x=new BloomFilterWrapper(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public BloomFilterWrapper(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean setInterleaved=false; //Whether interleaved was explicitly set.
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Tools.max(Shared.threads()>1 ? 2 : 1, Shared.threads()>20 ? Shared.threads()/2 : Shared.threads()));
		SamLine.SET_FROM_OK=true;
		
		//Create a parser object
		Parser parser=new Parser();
		
		boolean setBits=false;
		int k_=31;
		int hashes_=2;
		int minConsecutiveMatches_=3;
		int bits_=1;
		int minCount_=1;
		boolean rcomp_=true;
		boolean requireBothToMatch_=false;
		boolean ecco_=false;
		boolean merge_=false;
		
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
			
			else if(a.equalsIgnoreCase("k") || a.equalsIgnoreCase("bloomK") || a.equalsIgnoreCase("bloomFilterK")){
				k_=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("hashes") || a.equalsIgnoreCase("bloomHashes") || a.equalsIgnoreCase("bloomFilterHashes")){
				hashes_=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("minhits") || a.equalsIgnoreCase("bloomMinHits") || a.equalsIgnoreCase("bloomFilterMinHits")){
				minConsecutiveMatches_=Integer.parseInt(b);
			}else if(a.equals("rcomp")){
				rcomp_=Parse.parseBoolean(b);
			}else if(a.equals("bits")){
				setBits=true;
				bits_=Integer.parseInt(b);
				assert(bits_>0);
			}else if(a.equals("mincount")){
				minCount_=Integer.parseInt(b);
//				assert(minCount_>=1) : "mincount must be at least 1.";
//				minCount_=Tools.max(1, minCount_);
			}else if(a.equals("minprob")){
				ReadCounter.minProb=Float.parseFloat(b);
			}else if(a.equals("requireboth")){
				requireBothToMatch_=Parse.parseBoolean(b);
			}else if(a.equals("ecco")){
				ecco_=Parse.parseBoolean(b);
			}else if(a.equals("merge")){
				merge_=Parse.parseBoolean(b);
			}else if(a.equals("memfraction") || a.equals("memmult") || a.equals("memratio")){
				memFraction=Float.parseFloat(b);
			}else if(a.equals("cells")){
				BloomFilter.OVERRIDE_CELLS=Parse.parseKMG(b);
			}else if(a.equals("seed")){
				KCountArray7MTA.setSeed(Parse.parseKMG(b));
			}
			
			else if(a.equals("serialin")){
				serialIn=b;
			}else if(a.equals("serialout")){
				serialOut=b;
			}
			
			else if(a.equals("ref")){
				boolean result=Tools.addFiles(b, ref);
				assert(result) : "Could not add file(s): "+b;
			}else if(a.equals("outm") || a.equals("outm1") || a.equals("outbad") || a.equals("outbad1") || a.equals("outmatch") || a.equals("outmatch1")){
				outm1=b;
			}else if(a.equals("outm2") || a.equals("outbad2") || a.equals("outmatch2")){
				outm2=b;
			}else if(a.equals("outu") || a.equals("outu1") || a.equals("out") || a.equals("out1")){
				out1=b;
			}else if(a.equals("outu2") || a.equals("out2")){
				out2=b;
			}else if(a.equals("outc") || a.equals("counts")){
				outc=b;
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
		
		while(minCount_>0 && (1L<<bits_)-1<minCount_){bits_*=2;}
		
		k=k_;
		bits=bits_;
		hashes=hashes_;
		minConsecutiveMatches=minConsecutiveMatches_;
		minCount=minCount_;
		rcomp=rcomp_;
		requireBothToMatch=requireBothToMatch_;
		ecco=ecco_;
		merge=merge_;
		KmerCountAbstract.CANONICAL=rcomp; //rcomp, or true, perhaps...  hmmm.
		
		shift=bitsPerBase*k;
		shift2=shift-bitsPerBase;
		mask=(shift>63 ? -1L : ~((-1L)<<shift));
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;
			
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
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
		
		//Do output file # replacement
		if(outm1!=null && outm2==null && outm1.indexOf('#')>-1){
			outm2=outm1.replace("#", "2");
			outm1=outm1.replace("#", "1");
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
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outm1, outm2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1, out2, outm1, outm2)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		
		//Create output FileFormat objects
		ffoutm1=FileFormat.testOutput(outm1, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		ffoutm2=FileFormat.testOutput(outm2, FileFormat.FASTQ, extout, true, overwrite, append, ordered);
		
		ffoutc=FileFormat.testOutput(outc, FileFormat.FASTA, extout, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTA, extin, true, true);
		
		{
			Timer t=new Timer(outstream, true);
			if(serialIn==null){
				filter=new BloomFilter(null, null, ref, k, k, bits, hashes, minConsecutiveMatches,
						rcomp, ecco, merge, memFraction);
				if(serialOut!=null){
					ReadWrite.writeObjectInThread(filter, serialOut, true);
				}
			}else{
				filter=ReadWrite.read(BloomFilter.class, serialIn, true);
			}
			t.stop("Filter creation: \t\t");
			outstream.println(filter.filter.toShortString());
			if(verbose){
				outstream.println("Keys Counted: "+KmerCountAbstract.keysCounted);
				outstream.println("Increments:   "+KmerCountAbstract.increments);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		//Optionally create read output streams
		final ConcurrentReadOutputStream rosu, rosm, rosc;
		if(ffout1!=null){
			//Select output buffer size based on whether it needs to be ordered
			final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);
			
			//Notify user of output mode
			if(cris.paired() && out2==null && (in1!=null && !ffin1.samOrBam() && !ffout1.samOrBam())){
				outstream.println("Writing interleaved.");
			}
			
			rosu=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, false);
			rosu.start(); //Start the stream
		}else{rosu=null;}
		
		if(ffoutm1!=null){
			//Select output buffer size based on whether it needs to be ordered
			final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);
			
			rosm=ConcurrentReadOutputStream.getStream(ffoutm1, ffoutm2, qfoutm1, qfoutm2, buff, null, false);
			rosm.start(); //Start the stream
		}else{rosm=null;}
		
		if(ffoutc!=null){
			//Select output buffer size based on whether it needs to be ordered
			final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);
			
			rosc=ConcurrentReadOutputStream.getStream(ffoutc, null, buff, null, false);
			rosc.start(); //Start the stream
		}else{rosc=null;}
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		Timer t2=new Timer(outstream, true);
		
		//Process the reads in separate threads
		spawnThreads(cris, rosu, rosm, rosc);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, rosu, rosm, rosc);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		t.stop();
		t2.stop("\nFiltering Time:  \t\t");
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		{
			long readsMatched=(readsProcessed-readsOut);
			long basesMatched=(basesProcessed-basesOut);
			double rpct=readsMatched*100.0/readsProcessed;
			double bpct=basesMatched*100.0/basesProcessed;
			String rstring=Tools.padKMB(readsMatched, 8);
			String bstring=Tools.padKMB(basesMatched, 8);
			StringBuilder sb=new StringBuilder();
			sb.append("Reads Matched:      ").append(rstring).append(Tools.format(" \t%.2f%%", rpct)).append('\n');
			sb.append("Bases Matched:      ").append(bstring).append(Tools.format(" \t%.2f%%", bpct));
			outstream.println(sb);
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, 
			final ConcurrentReadOutputStream rosu, final ConcurrentReadOutputStream rosm,
			final ConcurrentReadOutputStream rosc){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, rosu, rosm, rosc, i));
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
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		//Do anything necessary after processing
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, 
				final ConcurrentReadOutputStream rosu_, final ConcurrentReadOutputStream rosm_, 
				final ConcurrentReadOutputStream rosc_, final int tid_){
			cris=cris_;
			rosu=rosu_;
			rosm=rosm_;
			rosc=rosc_;
			tid=tid_;
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
				ArrayList<Read> matched=new ArrayList<Read>(reads.size());
				ArrayList<Read> counts=(rosc==null ? null : new ArrayList<Read>(reads.size()));

				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					Read r1=reads.get(idx);
					Read r2=r1.mate;
					final String r2id=r1.mateId();
					
					//Validate reads in worker threads
					if(!r1.validated()){r1.validate(true);}
					if(r2!=null && !r2.validated()){r2.validate(true);}

					//Track the initial length for statistics
					final int initialLength1=r1.length();
					final int initialLength2=r1.mateLength();

					//Increment counters
					readsProcessedT+=r1.pairCount();
					basesProcessedT+=initialLength1+initialLength2;
					
					
					//Reads are processed in this block.
					{
						if((ecco || merge) && r1!=null && r2!=null){
							if(merge){
								final int insert=BBMerge.findOverlapStrict(r1, r2, false);
								if(insert>0){
									r2.reverseComplement();
									r1=r1.joinRead(insert);
									r2=null;
								}
							}else if(ecco){
								BBMerge.findOverlapStrict(r1, r2, true);
							}
						}

						if(rosc!=null){
							addCounts(r1, counts);
							addCounts(r2, counts);
						}
						final boolean match=matchesFilter(r1, r2);
						if(match){
							reads.set(idx, null);
							matched.add(r1);
						}else{
							readsOutT+=r1.pairCount();
							basesOutT+=r1.pairLength();
						}
					}
				}

				//Output reads to the output streams
				if(rosu!=null){rosu.add(reads, ln.id);}
				if(rosm!=null){rosm.add(matched, ln.id);}
				if(rosc!=null){rosc.add(counts, ln.id);}

				//Notify the input stream that the list was used
				cris.returnList(ln);
//				if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access

				//Fetch a new list
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}

		/**
		 * Breaks a read into kmers and adds them, with their counts, to the new list.
		 * @param r Read to process
		 * @param counts Output list of kmers (bases) and counts (read name)
		 */
		private void addCounts(Read r, ArrayList<Read> counts){
			if(r==null){return;}
			final byte[] bases=r.bases;

			long kmer=0;
			long rkmer=0;
			int len=0;

			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				long x=AminoAcid.baseToNumber[b];
				long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<bitsPerBase)|x)&mask;
				rkmer=((rkmer>>>bitsPerBase)|(x2<<shift2))&mask;
				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{
					len++;
					if(len>=k){
						long key=rcomp ? Tools.max(kmer, rkmer) : kmer;
						int count=filter.getCount(key);
						if(count>=minCount){
							String s=AminoAcid.kmerToString(key, k);
							Read cread=new Read(s.getBytes(), null, ""+count, r.numericID);
							counts.add(cread);
						}
					}
				}
			}
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r1 Read 1
		 * @param r2 Read 2 (may be null)
		 * @return True if the reads should be kept (meaning they don't match the filter), false if they should be discarded.
		 */
		boolean matchesFilter(final Read r1, final Read r2){
			if(r2!=null && requireBothToMatch){
//				assert(false) : minCount+", "+filter.matches(r1, buffer, minCount)+", "+filter.matches(r2, buffer, minCount)+" "+filter.matchesEither(r1, r2, buffer, minCount);
				return filter.matches(r1, buffer, minCount) && filter.matches(r2, buffer, minCount);
			}else{
//				assert(false) : minCount+", "+filter.matches(r1, buffer, minCount)+", "+filter.matches(r2, buffer, minCount)+" "+filter.matchesEither(r1, r2, buffer, minCount)+"\n"+r1+"\n"+r2;
				return filter.matchesEither(r1, r2, buffer, minCount);
			}
		}
		
		private final LongList buffer=minConsecutiveMatches>1 ? new LongList(150) : null;
		
		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		/** Number of reads retained by this thread */
		protected long readsOutT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Shared output stream */
		private final ConcurrentReadOutputStream rosu;
		/** Matched output stream */
		private final ConcurrentReadOutputStream rosm;
		/** Count output stream */
		private final ConcurrentReadOutputStream rosc;
		/** Thread ID */
		final int tid;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private ArrayList<String> ref=new ArrayList<String>();
	
	/** Primary input file path */
	private String in1=null;
	/** Secondary input file path */
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	/** Primary output file path */
	private String out1=null;
	/** Secondary output file path */
	private String out2=null;
	
	/** Output for kmer counts */
	private String outc=null;

	private String qfout1=null;
	private String qfout2=null;

	/** Output file path for matching reads */
	private String outm1=null;
	/** Secondary output file path for matching reads */
	private String outm2=null;

	private String qfoutm1=null;
	private String qfoutm2=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;

	private String serialIn=null;
	private String serialOut=null;
	
	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of reads retained */
	protected long readsOut=0;
	/** Number of bases retained */
	protected long basesOut=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
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
	
	/** Output for kmers and counts */
	private final FileFormat ffoutc;
	
	final BloomFilter filter;
	
	final int k;
	final int hashes;
	final int bits;
	final int minConsecutiveMatches;
	final boolean rcomp;
	final boolean requireBothToMatch;
	final boolean ecco;
	final boolean merge;
	final int minCount;
	float memFraction=1.0f;
	
	final int bitsPerBase=BloomFilter.bitsPerBase;
	final int shift;
	final int shift2;
	final long mask;
	
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
