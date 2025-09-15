package var2;

import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import stream.SamReadStreamer;
import stream.SamStreamer;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Removes lines with unsupported variations from a sam file.
 * 
 * @author Brian Bushnell
 * @date January 25, 2018
 *
 */
public class FilterSam {
	
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
		FilterSam x=new FilterSam(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public FilterSam(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		SamLine.SET_FROM_OK=true;
		Var.CALL_INS=false;
		Var.CALL_DEL=false;
		
		//Create a parser object
		Parser parser=new Parser();
		
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
			}else if(a.equals("ss") || a.equals("streamer")){
				useStreamer=Parse.parseBoolean(b);
			}else if(a.equals("ploidy")){
				ploidy=Integer.parseInt(b);
			}else if(a.equals("prefilter")){
				prefilter=Parse.parseBoolean(b);
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("outbad") || a.equals("outb")){
				outBad=b;
			}else if(a.equals("vars") || a.equals("variants") || a.equals("varfile") || a.equals("inv")){
				varFile=b;
			}else if(a.equals("vcf") || a.equals("vcffile")){
				vcfFile=b;
			}else if(a.equals("maxsubs") || a.equals("subfilter")){
				maxSubs=Integer.parseInt(b);
			}else if(a.equals("maxvars") || a.equals("varfilter")){
				maxVars=Integer.parseInt(b);
			}else if(a.equals("maxbadsubs") || a.equals("maxbadvars") || a.equals("mbv")){
				maxBadVars=Integer.parseInt(b);
			}else if(a.equals("maxbadsubdepth") || a.equals("maxbadvardepth") || a.equals("maxbadsuballeledepth")
					|| a.equals("maxbadvaralleledepth") || a.equals("mbsad") || a.equals("mbvad") || a.equals("mbad")){
				maxBadAlleleDepth=Integer.parseInt(b);
			}else if(a.equals("maxbadallelefraction") || a.equals("mbaf")){
				maxBadAlleleFraction=Float.parseFloat(b);
			}else if(a.equals("minbadsubreaddepth") || a.equals("minbadreaddepth") || a.equals("mbsrd") || a.equals("mbrd")){
				minBadReadDepth=Integer.parseInt(b);
			}else if(a.equals("sub") || a.equals("subs")){
				Var.CALL_SUB=Parse.parseBoolean(b);
			}else if(a.equals("ins") || a.equals("inss")){
				Var.CALL_INS=Parse.parseBoolean(b);
			}else if(a.equals("del") || a.equals("dels")){
				Var.CALL_DEL=Parse.parseBoolean(b);
			}else if(a.equals("indel") || a.equals("indels")){
				Var.CALL_INS=Var.CALL_DEL=Parse.parseBoolean(b);
			}else if(a.equals("minedist") || a.equals("minenddist") || a.equals("border")){
				minEDist=Integer.parseInt(b);
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
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;

			outGood=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		//Ensure there is an input file
//		if(in1==null){throw new RuntimeException("Error - an input file and a VCF file are required.");}
		if(in1==null){throw new RuntimeException("Error - an input file is required.");}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, outGood, outBad)){
			outstream.println((outGood==null)+", "+(outBad==null));
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+outGood+", "+outBad+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, vcfFile, varFile)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, vcfFile, varFile, outBad, outGood)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		//Create output FileFormat objects
		ffoutGood=FileFormat.testOutput(outGood, FileFormat.SAM, null, true, overwrite, append, ordered);
		ffoutBad=FileFormat.testOutput(outBad, FileFormat.SAM, null, true, overwrite, append, ordered);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.SAM, null, true, true);
		
		Timer t=new Timer(outstream, true);
		{
			t.start();
			outstream.print("Loading scaffolds:  ");
			if(ref==null){
				scafMap=ScafMap.loadSamHeader(in1);
			}else{
				scafMap=ScafMap.loadReference(ref, true);
			}
			assert(scafMap!=null && scafMap.size()>0) : "No scaffold names were loaded.";
			t.stop(pad(scafMap.size(), 12)+" \t");
		}
		t.start();
		if(varFile!=null){
			outstream.print("Loading vars:       ");
			varMap=VcfLoader.loadVarFile(varFile, scafMap);
			t.stop(pad(varMap.size(), 12)+" \t");
		}else if(vcfFile!=null){t.start();
			outstream.print("Loading vcf:        ");
			varMap=VcfLoader.loadVcfFile(vcfFile, scafMap, true, false);
			t.stop(pad(varMap.size(), 12)+" \t");
		}else{
			outstream.print("Calling variants:\n");
			String inString="in="+in1;
			
//			private int maxBadAlleleDepth=2;
//			/** Maximum allele fraction for variants considered unsupported */
//			private float maxBadAlleleFraction=0.01f;
			
			String[] cvargs=new String[] {inString, "ref="+ref, "clearfilters", "minreads="+(maxBadAlleleDepth), 
					"minallelefraction="+maxBadAlleleFraction, "printexecuting=f"};

			boolean trd=Shared.TRIM_READ_COMMENTS;
			boolean trn=Shared.TRIM_RNAME;
			CallVariants cv=new CallVariants(cvargs);
			Shared.TRIM_READ_COMMENTS=trd;
			Shared.TRIM_RNAME=trn;
			
			cv.prefilter=prefilter;
			cv.ploidy=ploidy;

			varMap=cv.process(new Timer());
			scafMap=cv.scafMap;
			t.stop(/*pad(varMap.size(), 12)+" \t"*/);
		}
		
		assert(Var.CALL_INS || Var.CALL_DEL || Var.CALL_SUB) : "Must enable at least one of subs, insertions, or deletions.";
		subsOnly=!Var.CALL_INS && !Var.CALL_DEL;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		
		final SamReadStreamer ss;
		//Create a read input stream
		final ConcurrentReadInputStream cris;
		if(useStreamer){
			cris=null;
			ss=new SamReadStreamer(ffin1, streamerThreads, true, maxReads);
			ss.start();
			if(verbose){outstream.println("Started streamer");}
		}else{
			ss=null;
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			cris.start(); //Start the stream
			if(verbose){outstream.println("Started cris");}
		}
		
		//Optionally create a read output stream
		final ConcurrentReadOutputStream ros, rosb;
		if(ffoutGood!=null){
			//Select output buffer size based on whether it needs to be ordered
			final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);
			
			ros=ConcurrentReadOutputStream.getStream(ffoutGood, null, null, null, buff, null, true);
			ros.start(); //Start the stream
		}else{ros=null;}
		if(ffoutBad!=null){
			//Select output buffer size based on whether it needs to be ordered
			final int buff=(ordered ? Tools.mid(16, 128, (Shared.threads()*2)/3) : 8);
			
			rosb=ConcurrentReadOutputStream.getStream(ffoutBad, null, null, null, buff, null, false);
			rosb.start(); //Start the stream
		}else{rosb=null;}
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the reads in separate threads
		spawnThreads(cris, ss, ros, rosb);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris, ros, rosb);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		//Report timing and results
		{
			t.stop();
			
			//Calculate units per nanosecond
			double rpnano=readsProcessed/(double)(t.elapsed);
			double bpnano=basesProcessed/(double)(t.elapsed);

			final long rg=readsOut, rb=readsProcessed-readsOut;
			final long bg=basesOut, bb=basesProcessed-basesOut;
			
			//Add "k" and "m" for large numbers
//			String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
//			String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");
//			String rgstring=(rg<100000 ? ""+rg : rg<100000000 ? (rg/1000)+"k" : (rg/1000000)+"m");
//			String bgstring=(bg<100000 ? ""+bg : bg<100000000 ? (bg/1000)+"k" : (bg/1000000)+"m");
//			String rbstring=(rb<100000 ? ""+rb : rb<100000000 ? (rb/1000)+"k" : (rb/1000000)+"m");
//			String bbstring=(bb<100000 ? ""+bb : bb<100000000 ? (bb/1000)+"k" : (bb/1000000)+"m");
			
			final int len=12;
			
			final String rpstring=pad(readsProcessed, len);
			final String bpstring=pad(basesProcessed, len);
			final String rgstring=pad(rg, len);
			final String bgstring=pad(bg, len);
			final String rbstring=pad(rb, len);
			final String bbstring=pad(bb, len);
			
			final double mappedReadsDiscarded=mappedReadsProcessed-mappedReadsRetained;
			double avgQGood=(bqSumGood/(double)mappedReadsRetained);
			double avgQBad=(bqSumBad/(double)mappedReadsDiscarded);
			double avgMapQGood=(mapqSumGood/(double)mappedReadsRetained);
			double avgMapQBad=(mapqSumBad/(double)mappedReadsDiscarded);
			double avgVarsGood=(varSumGood/(double)mappedReadsRetained);
			double avgVarsBad=(varSumBad/(double)mappedReadsDiscarded);
			final String avgQGoodS=pad(Tools.format("%.2f", avgQGood), len);
			final String avgQBadS=pad(Tools.format("%.2f", avgQBad), len);
			final String avgMapQGoodS=pad(Tools.format("%.2f", avgMapQGood), len);
			final String avgMapQBadS=pad(Tools.format("%.2f", avgMapQBad), len);
			final String avgVarsGoodS=pad(Tools.format("%.2f", avgVarsGood), len);
			final String avgVarsBadS=pad(Tools.format("%.2f", avgVarsBad), len);
			
			outstream.println("Time:                         \t\t"+t);
			outstream.println("Reads Processed:    "+rpstring+" \t"+Tools.format("%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Processed:    "+bpstring+" \t"+Tools.format("%.2fm bases/sec", bpnano*1000));
			outstream.println();
			outstream.println("Reads Retained:     "+rgstring+" \t"+Tools.format("%.2f%%", (rg*100.0/readsProcessed)));
			outstream.println("Bases Retained:     "+bgstring+" \t"+Tools.format("%.2f%%", (bg*100.0/basesProcessed)));
			outstream.println("Avg. Qual Retained: "+avgQGoodS);
			outstream.println("Avg. MapQ Retained: "+avgMapQGoodS);
			outstream.println("Avg. Vars Retained: "+avgVarsGoodS);
			outstream.println();
			outstream.println("Reads Discarded:    "+rbstring+" \t"+Tools.format("%.2f%%", (rb*100.0/readsProcessed)));
			outstream.println("Bases Discarded:    "+bbstring+" \t"+Tools.format("%.2f%%", (bb*100.0/basesProcessed)));
			outstream.println("Avg. Qual Discarded:"+avgQBadS);
			outstream.println("Avg. MapQ Discarded:"+avgMapQBadS);
			outstream.println("Avg. Vars Discarded:"+avgVarsBadS);
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private static String pad(long s, int len){
		return pad(""+s, len);
	}
	
	private static String pad(String s, int len){
		while(s.length()<len){s=" "+s;}
		return s;
	}
	
	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, final SamStreamer ss, final ConcurrentReadOutputStream ros, final ConcurrentReadOutputStream rosb){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, ss, ros, rosb, i));
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
			mappedReadsProcessed+=pt.mappedReadsProcessedT;
			mappedBasesProcessed+=pt.mappedBasesProcessedT;
			mappedReadsRetained+=pt.mappedReadsRetainedT;
			mappedBasesRetained+=pt.mappedBasesRetainedT;
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			bqSumGood+=pt.qSumGoodT;
			bqSumBad+=pt.qSumBadT;
//			adSumGood+=pt.adSumGoodT;
//			adSumBad+=pt.adSumBadT;
//			rdSumGood+=pt.rdSumGoodT;
//			rdSumBad+=pt.rdSumBadT;
			varSumGood+=pt.varSumGoodT;
			varSumBad+=pt.varSumBadT;
			mapqSumGood+=pt.mapqSumGoodT;
			mapqSumBad+=pt.mapqSumBadT;
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
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final SamStreamer ss_, final ConcurrentReadOutputStream ros_, final ConcurrentReadOutputStream rosb_, final int tid_){
			cris=cris_;
			ss=ss_;
			ros=ros_;
			rosb=rosb_;
			tid=tid_;
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the reads
			if(useStreamer){
				processSS();
			}else{
				processCris();
			}
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Iterate through the reads */
		void processCris(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
//				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}

			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access

				ArrayList<Read> good=new ArrayList<Read>(reads.size());
				ArrayList<Read> bad=new ArrayList<Read>(Tools.max(1, reads.size()/4));
				
				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					//Validate reads in worker threads
					if(!r1.validated()){r1.validate(true);}

					//Track the initial length for statistics
					final int initialLength1=r1.length();

					//Increment counters
					readsProcessedT++;
					basesProcessedT+=initialLength1;
					
					{
						//Reads are processed in this block.
						boolean keep=processRead(r1);
						if(keep){
							//Increment counters
							readsOutT++;
							basesOutT+=initialLength1;
							good.add(r1);
						}else{
							bad.add(r1);
						}
					}
				}

				//Output reads to the output stream
				if(ros!=null){ros.add(good, ln.id);}
				if(rosb!=null){rosb.add(bad, ln.id);}

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
		
		/** Iterate through the reads */
		void processSS(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=ss.nextList();
			//Grab the actual read list from the ListNum
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//Check to ensure pairing is as expected
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
//				assert(ffin1.samOrBam() || (r.mate!=null)==cris.paired()); //Disabled due to non-static access
			}

			//As long as there is a nonempty read list...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access

				ArrayList<Read> good=new ArrayList<Read>(reads.size());
				ArrayList<Read> bad=new ArrayList<Read>(Tools.max(1, reads.size()/4));
				
				//Loop through each read in the list
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					//Validate reads in worker threads
					if(!r1.validated()){r1.validate(true);}

					//Track the initial length for statistics
					final int initialLength1=r1.length();

					//Increment counters
					readsProcessedT++;
					basesProcessedT+=initialLength1;
					
					{
						//Reads are processed in this block.
						boolean keep=processRead(r1);
						if(keep){
							//Increment counters
							readsOutT++;
							basesOutT+=initialLength1;
							good.add(r1);
						}else{
							bad.add(r1);
						}
					}
				}

				//Output reads to the output stream
				if(ros!=null){ros.add(good, ln.id);}
				if(rosb!=null){rosb.add(bad, ln.id);}

				//Fetch a new list
				ln=ss.nextList();
				reads=(ln!=null ? ln.list : null);
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		/**
		 * Process a read.
		 * @param r1 Read 1
		 * @return True if the read should be kept, false if it should be discarded.
		 */
		boolean processRead(final Read r1){
			return passesVariantFilter(r1);
		}
		
		private final boolean passesVariantFilter(Read r){
			if(!r.mapped() || r.bases==null || r.samline==null || r.match==null){return true;}
			final int vars=(subsOnly ? Read.countSubs(r.match) : Read.countVars(r.match, Var.CALL_SUB, Var.CALL_INS, Var.CALL_DEL));
			final int len=r.length();
			final double q=r.avgQualityByProbabilityDouble(false, r.length());
			final SamLine sl=r.samline;
			mappedReadsProcessedT++;
			mappedBasesProcessedT+=len;
			
			if(vars==0){
				varSumGoodT+=vars;
				qSumGoodT+=q;
				mapqSumGoodT+=sl.mapq;
				mappedReadsRetainedT++;
				mappedBasesRetainedT+=len;
				return true;
			}
			
			if((maxVars>=0 && vars>maxVars) || 
					(maxSubs>=0 && vars>maxSubs && ((subsOnly ? vars : Read.countSubs(r.match))>maxSubs))){
				varSumBadT+=vars;
				qSumBadT+=q;
				mapqSumBadT+=sl.mapq;
				return false;
			}
			
			if(vars<=maxBadVars){
				varSumGoodT+=vars;
				qSumGoodT+=q;
				mapqSumGoodT+=sl.mapq;
				mappedReadsRetainedT++;
				mappedBasesRetainedT+=len;
				return true;
			}
			
			final ArrayList<Var> list;
			if(subsOnly){
				list=AnalyzeVars.findUniqueSubs(r, sl, varMap, scafMap, maxBadAlleleDepth, maxBadAlleleFraction, minBadReadDepth, minEDist);
			}else{
				list=AnalyzeVars.findUniqueVars(r, sl, varMap, scafMap, maxBadAlleleDepth, maxBadAlleleFraction, minBadReadDepth, minEDist);
			}
			if(list==null || list.size()<=maxBadVars){
				varSumGoodT+=vars;
				qSumGoodT+=q;
				mapqSumGoodT+=sl.mapq;
				mappedReadsRetainedT++;
				mappedBasesRetainedT+=len;
				return true;
			}else{
				assert(list.size()>maxBadVars) : list.size()+", "+maxBadVars;
				for(Var v : list){
					assert(v!=null) : sl.cigar+"\n"+vars+"\n"+v+"\n"+list;
					assert((v.type()==Var.SUB && Var.CALL_SUB) || (v.type()==Var.INS && Var.CALL_INS) || (v.type()==Var.DEL && Var.CALL_DEL)) : 
						list.size()+", "+maxBadVars+"\n"+sl.cigar+"\n"+v+"\n"+list;
				}
				varSumBadT+=vars;
				qSumBadT+=q;
				mapqSumBadT+=sl.mapq;
				return false;
			}
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of reads processed by this thread */
		protected long basesProcessedT=0;
		/** Number of mapped reads processed by this thread */
		protected long mappedReadsProcessedT=0;
		/** Number of mapped bases processed by this thread */
		protected long mappedBasesProcessedT=0;
		/** Number of mapped reads retained by this thread */
		protected long mappedReadsRetainedT=0;
		/** Number of mapped bases retained by this thread */
		protected long mappedBasesRetainedT=0;
		/** Number of good reads processed by this thread */
		protected long readsOutT=0;
		/** Number of good bases processed by this thread */
		protected long basesOutT=0;
		protected double qSumGoodT=0;
		protected double qSumBadT=0;
		protected long varSumGoodT=0;
		protected long varSumBadT=0;
		protected long mapqSumGoodT=0;
		protected long mapqSumBadT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;

		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Shared input stream */
		private final SamStreamer ss;
		/** Good output stream */
		private final ConcurrentReadOutputStream ros;
		/** Bad output stream */
		private final ConcurrentReadOutputStream rosb;
		/** Thread ID */
		final int tid;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;
	
	/** Optional reference */
	private String ref=null;

	/** Good output file path */
	private String outGood=null;
	/** Bad output file path */
	private String outBad=null;

	/** VCF file path */
	private String varFile=null;
	/** Variant file path */
	private String vcfFile=null;
	private VarMap varMap=null;
	private ScafMap scafMap=null;
	
	/*--------------------------------------------------------------*/
	
	/** Maximum allowed substitutions in a read */
	private int maxSubs=-1;
	/** Maximum allowed variations in a read */
	private int maxVars=-1;
	
	/** Maximum allowed unsupported substitutions in a read */
	private int maxBadVars=1;
	/** Maximum variant depth for a variant to be considered unsupported */
	private int maxBadAlleleDepth=2;
	/** Maximum allele fraction for variants considered unsupported */
	private float maxBadAlleleFraction=0.01f;
	/** Minimum read depth for a variant to be considered unsupported */
	private int minBadReadDepth=2;
	/** Ignore vars within this distance of the ends */
	private int minEDist=5;
	private int ploidy=1;
	private boolean prefilter=false;

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;
	/** Number of mapped reads processed by this thread */
	protected long mappedReadsProcessed=0;
	/** Number of mapped bases processed by this thread */
	protected long mappedBasesProcessed=0;
	/** Number of mapped reads retained by this thread */
	protected long mappedReadsRetained=0;
	/** Number of mapped bases retained by this thread */
	protected long mappedBasesRetained=0;
	/** Number of good reads processed */
	protected long readsOut=0;
	/** Number of good bases processed */
	protected long basesOut=0;

	protected double bqSumGood=0;
	protected double bqSumBad=0;
//	protected double adSumGood=0;
//	protected double adSumBad=0;
//	protected double rdSumGood=0;
//	protected double rdSumBad=0;
	
	protected long varSumGood=0;
	protected long varSumBad=0;
	protected long mapqSumGood=0;
	protected long mapqSumBad=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	static boolean useStreamer=true;
	static int streamerThreads=3;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	
	/** Primary output file */
	private final FileFormat ffoutGood;
	/** Secondary output file */
	private final FileFormat ffoutBad;
	
	private final boolean subsOnly;
	
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
	private boolean ordered=true;
	
}
