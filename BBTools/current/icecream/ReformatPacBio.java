package icecream;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;

import aligner.Aligner;
import consensus.BaseGraph;
import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import json.JsonObject;
import prok.GeneCaller;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import structures.ByteBuilder;
import structures.IntHashSet;
import structures.IntList;
import tracker.EntropyTracker;
import tracker.ReadStats;

/**
 * Version of Reformat designed for PacBio data.
 * Supports some of Reformat's capability, like subsampling,
 * in a ZMW-aware tool.
 * 
 * @author Brian Bushnell
 * @date June 5, 2019
 *
 */
public final class ReformatPacBio {
	
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
		ReformatPacBio x=new ReformatPacBio(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public ReformatPacBio(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_BGZIP=ReadWrite.USE_UNBGZIP=ReadWrite.PREFER_BGZIP=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		SamLine.SET_FROM_OK=true;
		Shared.setBufferData(1000000);
//		Shared.FASTA_WRAP=511;
		Data.USE_SAMBAMBA=false;//Sambamba changes PacBio headers.
		Read.CHANGE_QUALITY=false;
		EntropyTracker.defaultK=3;
		
		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in1=parser.in1;
			extin=parser.extin;

			if(outg==null){outg=parser.out1;}
			extout=parser.extout;
		}
		
		sampleReadsExact=sampleReadsTarget>-1;
		sampleBasesExact=sampleBasesTarget>-1;
		sampleZMWsExact=sampleZMWsTarget>-1;
		sampleExact=(sampleReadsExact || sampleBasesExact || sampleZMWsExact);
		
		//Determine how many threads may be used
		threads=(sampleExact || (samplerate>0 && seed>=0)) ? 1 : Shared.threads();
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffoutg=FileFormat.testOutput(outg, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffoutb=FileFormat.testOutput(outb, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffstats=FileFormat.testOutput(outstats, FileFormat.TXT, null, false, overwrite, append, false);
		ffschist=FileFormat.testOutput(schist, FileFormat.TXT, null, false, overwrite, append, false);
		
		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		if(verbose){System.err.println("Finished constructor for "+getClass().getName());}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Parse arguments from the command line */
	private Parser parse(String[] args){
		//Create a parser object
		Parser parser=new Parser();
		
		//Set any necessary Parser defaults here
		//parser.foo=bar;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("format")){
				if(b==null){
					assert(false) : arg;
				}else if(Tools.isDigit(b.charAt(i))){
					format=Integer.parseInt(b);
				}else if(b.equalsIgnoreCase("json")){
					format=FORMAT_JSON;
				}else if(b.equalsIgnoreCase("text")){
					format=FORMAT_TEXT;
				}else{
					assert(false) : arg;
				}
				assert(format>=1 && format<=2) : arg;
			}else if(a.equals("json")){
				boolean x=Parse.parseBoolean(b);
				format=(x ? FORMAT_JSON : FORMAT_TEXT);
			}else if(a.equals("ss") || a.equals("samstreamer") || a.equals("streamer")){
				if(b!=null && Tools.isDigit(b.charAt(0))){
					ZMWStreamer.useStreamer=true;
					assert(Integer.parseInt(b)==1) : "ZMWStreamer threads currently capped at 1.";
//					ZMWStreamer.streamerThreads=Tools.max(1, Integer.parseInt(b));
				}else{
					ZMWStreamer.useStreamer=Parse.parseBoolean(b);
				}
			}
			
			else if(a.equals("keepshortreads") || a.equals("ksr")){
				keepShortReads=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("keepzmwstogether") || a.equals("kzt") || a.equals("keepreadstogether") || a.equals("krt")){
				keepZMWsTogether=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("subsampleFromEnds")){
				subsampleFromEnds=Parse.parseBoolean(b);
			}
			
			else if(a.equals("ccsinput") || a.equals("ccsin")){
				CCSInput=Parse.parseBoolean(b);
			}else if(a.equals("minlength") || a.equals("minlen")){
				minLengthAfterTrimming=Integer.parseInt(b);
			}else if(a.equals("flaglongreads")){
				flagLongReads=Parse.parseBoolean(b);
			}else if(a.equals("longreadmult")){
				longReadMult=Float.parseFloat(b);
			}else if(a.equals("trimreads") || a.equals("trim")){
				trimReads=Parse.parseBoolean(b);
			}else if(a.equals("parsecustom")){
				parseCustom=Parse.parseBoolean(b);
			}else if(a.equals("outg") || a.equals("outgood")){
				outg=b;
			}else if(a.equals("outb") || a.equals("outbad")){
				outb=b;
			}else if(a.equals("outs") || a.equals("outstats") || a.equals("stats")){
				outstats=b;
			}else if(a.equals("schist") || a.equals("shist")){
				schist=b;
			}
			
			else if(a.equals("shredlen")){
				shredLength=Integer.parseInt(b);
			}else if(a.equals("overlap")){
				overlap=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("minShredIdentity") || a.equalsIgnoreCase("minShredId") ||
					 a.equalsIgnoreCase("shredId")){
				minShredIdentity=Tools.max(0.01f, Float.parseFloat(b));
				if(minShredIdentity>1){minShredIdentity*=0.01f;}
				assert(minShredIdentity>0 && minShredIdentity<=1) : 
					"minShredIdentity ranges from 0 (exclusive) to 1 (inclusive).";
			}else if(a.equals("ccs") || a.equals("makeccs") || a.equals("consensus")){
				makeCCS=Parse.parseBoolean(b);
			}else if(a.equals("findorientation") || a.equals("orient") || a.equals("reorient")){
				findOrientation=Parse.parseBoolean(b);
			}
			
			else if(a.equals("minpasses")){
				minPasses=Float.parseFloat(b);
			}else if(a.equals("minsubreads")){
				minSubreads=Integer.parseInt(b);
			}
			
			else if(a.equalsIgnoreCase("samplerate")){
				samplerate=Float.parseFloat(b);
			}else if(a.equalsIgnoreCase("sampleReadsTarget") || a.equals("srt")){
				sampleReadsTarget=Parse.parseKMG(b);
			}else if(a.equalsIgnoreCase("sampleBasesTarget") || a.equals("sbt")){
				sampleBasesTarget=Parse.parseKMG(b);
			}else if(a.equalsIgnoreCase("sampleZMWsTarget") || a.equals("szt")){
				sampleZMWsTarget=Parse.parseKMG(b);
			}else if(a.equals("bestpass") || a.equals("keepbestpass") || 
					a.equals("keepbest") || a.equals("bestpassonly")){
				keepBestPass=Parse.parseBoolean(b);
			}else if(a.equals("longestpass") || a.equals("keeplongestpass") || 
					a.equals("keeplongest") || a.equals("longestpassonly")){
				keepLongestPass=Parse.parseBoolean(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}
			
			else if(a.equalsIgnoreCase("zmws") || a.equalsIgnoreCase("maxzmws")){
				maxZMWs=Parse.parseKMG(b);
			}
			
			else if(a.equals("whitelist")){
				whitelist=b;
			}else if(a.equals("whitelist")){
				blacklist=b;
			}
			
			else if(a.equalsIgnoreCase("trimpolya")){
				trimPolyA=Parse.parseBoolean(b);
			}else if(PolymerTrimmer.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("minentropy") || a.equals("entropy") || a.equals("entropyfilter") || a.equals("efilter")){
				if(b==null || Character.isLetter(b.charAt(0))){
					if(Parse.parseBoolean(b)){
						entropyCutoff=0.55f;
					}else{
						entropyCutoff=-1;
					}
				}else{
					entropyCutoff=Float.parseFloat(b);
				}
			}else if(a.equals("entropyblock") || a.equals("entropylength") || a.equals("entropylen") || a.equals("entlen")){
				entropyLength=Parse.parseIntKMG(b);
			}else if(a.equals("entropyfraction") || a.equals("entfraction")){
				entropyFraction=Float.parseFloat(b);
			}else if(a.equals("monomerfraction") || a.equals("maxmonomerfraction") || a.equals("mmf")){
				maxMonomerFraction=Float.parseFloat(b);
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

		if(verbose){System.err.println("Finished parser for "+getClass().getName());}
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, outg, outb, outstats, schist)){
			outstream.println((outg==null)+", "+(outb==null)+", "+outg+", "+outb+", "+outstats);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+
					outg+", "+outb+", "+outstats+", "+schist+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, whitelist, blacklist)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, outg, outb, outstats, schist, whitelist, blacklist)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Create read streams and process all data */
	void process(Timer t){
		if(verbose){System.err.println("Entered process()");}

		whiteSet=Tools.loadIntSet(whitelist);
		blackSet=Tools.loadIntSet(blacklist);

		//Count reads in the original file
		if(sampleExact){countInitial();}
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Create a read input stream
		ZMWStreamer zstream=new ZMWStreamer(ffin1, Shared.threads(), maxReads, maxZMWs);
		
		//Optionally create read output streams
		final ConcurrentReadOutputStream rosg=makeCros(ffoutg);
		final ConcurrentReadOutputStream rosb=makeCros(ffoutb);
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		
		//Process the reads in separate threads
		spawnThreads(zstream, rosg, rosb);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(null, rosg, rosb);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		writeHistogram(ffschist, subreadCounts);
		
		//Report timing and results
		t.stop();
		
		String stats=null;
		if(format==FORMAT_TEXT){
			ByteBuilder bb=toText(t);
			stats=bb.nl().toString();
		}else if(format==FORMAT_JSON){
			JsonObject jo=toJson(t);
			stats=jo.toStringln();
		}else{
			assert(false) : "Bad format: "+format;
		}
		if(ffstats==null){
			outstream.print(stats);
		}else{
			ReadWrite.writeString(stats, outstats);
		}
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	private ByteBuilder toText(Timer t){
		ByteBuilder bb=new ByteBuilder();
		bb.appendln(Tools.timeZMWsReadsBasesProcessed(t, ZMWs, readsProcessed, basesProcessed, 8));
		bb.appendln(Tools.ZMWsReadsBasesOut(ZMWs, readsProcessed, basesProcessed, ZMWsOut, readsOut, basesOut, 8, true));
		
		long readsFiltered=readsProcessed-readsOut;
		bb.appendln(Tools.numberPercent("Reads Filtered:", readsFiltered, readsFiltered*100.0/(readsProcessed), 3, 8));
		if(trimReads || trimPolyA){
			bb.appendln(Tools.numberPercent("Reads Trimmed:", readsTrimmed, readsTrimmed*100.0/(readsProcessed), 3, 8));
			bb.appendln(Tools.numberPercent("Bases Trimmed:", basesTrimmed, basesTrimmed*100.0/(basesProcessed), 3, 8));
		}
//		bb.appendln(Tools.number("Total ZMWs:", ZMWs, 8));
		bb.appendln(Tools.numberPercent("Partial ZMWs:", partiallyDiscardedZMWs, partiallyDiscardedZMWs*100.0/(ZMWs), 3, 8));
		bb.appendln(Tools.numberPercent("Discarded ZMWs:", fullyDiscardedZMWs, fullyDiscardedZMWs*100.0/(ZMWs), 3, 8));
//		bb.appendln(Tools.numberPercent("Low Entropy:", lowEntropyReads, lowEntropyReads*100.0/(readsProcessed), 3, 8));
		if(entropyCutoff>=0){
			bb.appendln(Tools.numberPercent("Low Entropy:", lowEntropyZMWs, lowEntropyZMWs*100.0/(ZMWs), 3, 8));
		}
		
		if(parseCustom){}
		return bb;
	}
	
	private JsonObject toJson(Timer t){
		JsonObject jo=new JsonObject();
		long readsFiltered=readsProcessed-readsOut;
		
//		asdf
		jo.add("Time", t.timeInSeconds());
		jo.add("ZMWs_Processed", ZMWs);
		jo.add("Reads_Processed", readsProcessed);
		jo.add("Bases_Processed", basesProcessed);
		jo.add("Reads_Out", readsOut);
		jo.add("Bases_Out", basesOut);
		jo.add("ZMWs_Out", ZMWsOut);
		jo.add("Reads_Filtered", readsFiltered);
		jo.add("Reads_Filtered_Pct", readsFiltered*100.0/(readsProcessed));
		if(trimReads){
			jo.add("Reads_Trimmed", readsTrimmed);
			jo.add("Reads_Trimmed_Pct", readsTrimmed*100.0/(readsProcessed));
			jo.add("Bases_Trimmed", basesTrimmed);
			jo.add("Bases_Trimmed_Pct", basesTrimmed*100.0/(basesProcessed));
		}
//		jo.add("Total_ZMWs", ZMWs);fullyDiscardedZMWs
		jo.add("Partial_ZMWs", partiallyDiscardedZMWs);
		jo.add("Partial_ZMWs_Pct", partiallyDiscardedZMWs*100.0/(ZMWs));
		jo.add("Discarded_ZMWs", fullyDiscardedZMWs);
		jo.add("Discarded_ZMWs_Pct", fullyDiscardedZMWs*100.0/(ZMWs));
//		jo.add("Low_Entropy", lowEntropyReads);
//		jo.add("Low_Entropy_Pct", lowEntropyReads*100.0/(readsProcessed));
		if(entropyCutoff>=0){
			jo.add("Low_Entropy", lowEntropyZMWs);
			jo.add("Low_Entropy_Pct", lowEntropyZMWs*100.0/(ZMWs));
		}
		if(parseCustom){
			//Special stats if desired
		}
		return jo;
	}
	
	private static void writeHistogram(FileFormat ff, long[] hist){
		if(ff==null){return;}
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();

		bsw.print("#Counted\t").println(shared.Vector.sum(hist));
		bsw.print("#Mean\t").println(Tools.averageHistogram(hist), 3);
		bsw.print("#Median\t").println(Tools.medianHistogram(hist));
		bsw.print("#Mode\t").println(Tools.calcModeHistogram(hist));
		bsw.print("#STDev\t").println(Tools.standardDeviationHistogram(hist), 3);
		bsw.print("#Value\tOccurances\n");
		
		for(int i=0; i<hist.length; i++){
			bsw.print(i).tab().println(hist[i]);
		}
		bsw.poisonAndWait();
	}
	
	private ConcurrentReadOutputStream makeCros(FileFormat ff){
		if(ff==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=16;

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(
				ff, null, buff, null, ff.samOrBam() && ffin1.samOrBam());
		ros.start(); //Start the stream
		return ros;
	}
	
	/** Spawn process threads */
	private void spawnThreads(final ZMWStreamer zstream,
			final ConcurrentReadOutputStream rosg,
			final ConcurrentReadOutputStream rosb){
		
		//Do anything necessary prior to processing
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(zstream, rosg, rosb, i));
		}
		
		//Start the threads
		for(ProcessThread pt : alpt){
			pt.start();
		}
		
		zstream.runStreamer(false);
		
		//Wait for threads to finish
		waitForThreads(alpt);
		
		//Do anything necessary after processing
		ZMWs=zstream.ZMWs;
	}
	
	private void waitForThreads(ArrayList<ProcessThread> alpt){
		
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
			ZMWs+=pt.ZMWsT;
			
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			ZMWsOut+=pt.ZMWsOutT;

			partiallyDiscardedZMWs+=pt.partiallyDiscardedZMWsT;
			fullyDiscardedZMWs+=pt.fullyDiscardedZMWsT;
			lowEntropyZMWs+=pt.lowEntropyZMWsT;
			lowEntropyReads+=pt.lowEntropyReadsT;
			
			basesTrimmed+=pt.basesTrimmedT;
			readsTrimmed+=pt.readsTrimmedT;
			
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          CCS Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	public Read makeConsensus(ZMW zmw){
		//Apply correct strand
		for(int i=0; i<zmw.size(); i++){
			zmw.get(i).setStrand(i&1);
		}
		final Read ref=zmw.medianRead(false);
		if(zmw.size()<3 || zmw.estimatePasses()<2.1){return ref;}
		
		BaseGraph bg=new BaseGraph(ref.id, ref.bases, ref.quality, ref.numericID, 0);
		float idSum=0;
		int added=0;
		for(Read r : zmw){
			if(r!=ref){
				final boolean rcomp=(r.strand()!=ref.strand());
				if(rcomp){r.reverseComplement();}
				float id=shredAndAdd(r, bg, null);
				idSum+=id;
				added++;
				if(rcomp){r.reverseComplement();}
			}
		}
		float avgId=idSum/(Tools.max(1, added));
		
		//TODO: Also interesting to get a traversal quality here...
		Read c=bg.traverse();
		c.obj=avgId;
		
		//If not enough subreads aligned, discard the result
		if(minPasses>1 && bg.baseTotal/(float)ref.length()<minPasses-1){c.setDiscarded(true);}
		return c;
	}
	
	private float shredAndAdd(Read subread, BaseGraph bg, Aligner ssa){
		int added=0;
		float idSum=0;
		if(ssa==null){ssa=GeneCaller.getSSA();}
		if(subread.length()<=shredLength){
			float id=bg.alignAndGenerateMatch(subread, ssa, findOrientation, minShredIdentity);
			added++;
			idSum+=id;
			if(id>=minShredIdentity){
				bg.add(subread);
			}
		}else{
			ArrayList<Read> shreds=shred(subread, shredLength, overlap);
			for(int i=0; i<shreds.size(); i++){
				Read shred=shreds.get(i);
				float id=bg.alignAndGenerateMatch(shred, ssa, findOrientation, minShredIdentity);
				added++;
				idSum+=id;
				if(id>=minShredIdentity){
					bg.add(shred, (i>0 ? overlap-1 : 0), 0);
				}
			}
		}
		float avgId=idSum/(Tools.max(1, added));
		return avgId;
	}
	
	private static ArrayList<Read> shred(final Read r, final int shredLength, final int overlap){
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		final String name=r.id;
		final int increment=shredLength-overlap;
		final double incMult=1.0f/increment;
		final int chunks=bases.length<=shredLength ? 1 : (int)Math.ceil((bases.length-overlap)*incMult);
		assert(chunks>0);
		final double inc2=bases.length/(double)chunks;
		
		final ArrayList<Read> list=new ArrayList<Read>(chunks);
		if(chunks==1){
			list.add(r);
			return list;
		}
		
		for(int chunk=0; chunk<chunks; chunk++){
			int a=(int)Math.floor(inc2*chunk);
			int b=(chunk==chunks-1 ? bases.length : overlap+(int)Math.floor(inc2*(chunk+1)));
			b=Tools.min(b, a+shredLength);
			final int length=b-a;
//			if(length<minLength){return;}
			final byte[] bases2=KillSwitch.copyOfRange(bases, a, b);
			final byte[] quals2=(quals==null ? null : KillSwitch.copyOfRange(quals, a, b));
			Read shred=new Read(bases2, quals2, name+"_"+a+"-"+(b-1), r.numericID);
//			readsOut++;
//			basesOut+=shred.length();
			list.add(shred);
		}
		return list;
	}
	
//	/** Generates match string and returns identity */
//	public float alignAndGenerateMatch(Read r, Aligner ssa){
//		if(ssa==null){ssa=GeneCaller.getSSA();}
//		byte[] query=r.bases;
//		int a=0, b=ref.length-1;
//		int[] max=ssa.fillUnlimited(query, original, a, b, -9999);
//		if(max==null){return 0;}
//		
//		final int rows=max[0];
//		final int maxCol=max[1];
//		final int maxState=max[2];
//		final byte[] match=ssa.traceback(query, original, a, b, rows, maxCol, maxState);
//		int[] score=ssa.score(query, original, a, b, rows, maxCol, 0);
//		r.match=match;
//		r.start=score[1];
//		r.stop=score[2];
//		r.setMapped(true);
//		final float identity=Read.identity(match);
//		return identity;
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	void countInitial(){
		if(verbose){System.err.println("Entered countInitial()");}
		assert(!ffin1.stdio()) : "Target subsampling can't be used with stdin.";
		initialReads=0;
		initialZMWs=0;
		initialBases=0;
		ZMWStreamer zstream=new ZMWStreamer(ffin1, Shared.threads(), maxReads, maxZMWs);
		zstream.runStreamer(true);
		for(ZMW zmw=zstream.nextZMW(); zmw!=null; zmw=zstream.nextZMW()){
			initialZMWs++;
			initialReads+=zmw.size();
			initialBases+=zmw.countBases();
		}
		readsRemaining=initialReads;
		ZMWsRemaining=initialZMWs;
		basesRemaining=initialBases;
		if(verbose){System.err.println("Finished countInitial()");}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Processes reads. */
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ZMWStreamer zstream_, 
				final ConcurrentReadOutputStream ros_,  
				final ConcurrentReadOutputStream rosb_, final int tid_){
			zstream=zstream_;
			rosg=ros_;
			rosb=rosb_;
			tid=tid_;
			
			if(entropyCutoff>=0){
				eTracker=new EntropyTracker(false, entropyCutoff, true);
			}else{
				eTracker=null;
			}
		}
		
		//Called by start()
		@Override
		public void run(){
			//Do anything necessary prior to processing
			randy=Shared.threadLocalRandom(seed<0 ? seed : seed+tid);
			//Randy can't be initialized in the constructor because it is threadlocal;
			//that leads to nonrandom results.
			
			//Process the reads
			processInner();
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/** Iterate through the reads */
		void processInner(){

			//As long as there is a nonempty read list...
			for(ZMW reads=zstream.nextZMW(); reads!=null; reads=zstream.nextZMW()){
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
				
				processList(reads);
			}
		}
		
		private int medianLength(ZMW list){
			if(list.size()<3){return -1;}
			IntList lengths=new IntList(list.size()-2);
			
			for(int i=1; i<list.size()-1; i++){
				lengths.add(list.get(i).length());
			}
			lengths.sort();
			int median=lengths.get((lengths.size-1)/2);
			return median;
		}
		
		int flagLowEntropyReads(final ZMW reads, final float minEnt, 
				final int minLen0, final float minFract){
			int found=0;
			for(Read r : reads){
				if(!r.discarded()){
					int minLen=Tools.min((int)(r.length()*minFract), minLen0);
					int maxBlock=eTracker.longestLowEntropyBlock(r.bases, false, maxMonomerFraction);
					if(maxBlock>=minLen){
						r.setDiscarded(true);
						r.setJunk(true);
						found++;
//						System.err.println(r.toFasta());
					}
				}
			}
			return found;
		}
		
		int flagLongReads(final ZMW reads, int median){
			int found=0;
			for(Read r : reads){
				if(r.length()>longReadMult*median){
					r.setDiscarded(true);
					r.setHasAdapter(true);
					found++;
				}
			}
			return found;
		}
		
		void discardEndReads(ZMW zmw, int toDiscard){
			if(toDiscard<1){return;}
			if(toDiscard==1){
				if(zmw.first().length()<zmw.last().length()){
					zmw.first().setDiscarded(true);
				}else{
					zmw.last().setDiscarded(true);
				}
				return;
			}else{
				zmw.first().setDiscarded(true);
				zmw.last().setDiscarded(true);
				for(int i=2; i<toDiscard; i++){
					zmw.get(zmw.size()-i).setDiscarded(true);
				}
			}
		}
		
		/** Each list is presumed to be all reads from a ZMW, in order */
		void processList(final ZMW reads){
			ZMWsT++;
			{
				int idx=Tools.min(subreadCounts.length-1, reads.size());
				synchronized(subreadCounts){
					//TODO: Slow
					//But other statistics could be gathered here as well
					subreadCounts[idx]++;
				}
			}
			long numBases=0;
			
			//Loop through each read in the list for statistics
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);
				
				//Validate reads in worker threads
				if(!r1.validated()){r1.validate(true);}

				//Track the initial length for statistics
				final int initialLength1=r1.length();

				//Increment counters
				readsProcessedT++;
				basesProcessedT+=initialLength1;
				numBases+=initialLength1;
			}
			final boolean fullPass=CCSInput || reads.size()>=3;

			if(whiteSet!=null || blackSet!=null){
				int zid=reads.zid();
				assert(zid>=0) : "Can't parse ZMW ID from header "+reads.get(0).id;
				if(whiteSet!=null && !whiteSet.contains(zid)){reads.setDiscarded(true);}
				if(blackSet!=null && blackSet.contains(zid)){reads.setDiscarded(true);}
			}

			final int subreads=reads.size();
			final float passes=reads.estimatePasses();
			if((minPasses>0 && passes<minPasses) || (minSubreads>0 && reads.size()<minSubreads)){
				reads.setDiscarded(true);
			}
			
			if(samplerate<1){
				if(keepZMWsTogether){
					if(randy.nextFloat()>samplerate){
						for(Read r : reads){r.setDiscarded(true);}
					}
				}else if(subsampleFromEnds){
					//Roll a fraction of the reads
					int remove=0;
					for(Read r : reads){
						if(randy.nextFloat()>samplerate){
							remove++;
						}
					}
					discardEndReads(reads, remove);
				}else{
					for(Read r : reads){
						if(randy.nextFloat()>samplerate){
							r.setDiscarded(true);
						}
					}
				}
			}
			
			if(keepBestPass || keepLongestPass){
				Read best=keepBestPass ? reads.medianRead(false) : reads.longestRead(false);
//				assert(reads.size()<5) : keepBestPass+", "+keepLongestPass+
//					", "+best.length()+", "+Arrays.toString(reads.lengths());
				for(Read r : reads){
					if(r!=best){
						r.setDiscarded(true);
					}
				}
			}else if(makeCCS && !reads.discarded()){
				Read r=makeConsensus(reads);
				reads.clear();
				reads.add(r);
				//TODO: Fix read header so that it looks like a CCS read.
				PBHeader pbh=new PBHeader(r.id);
				float avgId=(r.obj==null ? 0 : (Float)r.obj);
				r.id=pbh.runID+"/"+pbh.zmwID+"/ccs"
						+ "\tsubreads="+subreads+"\tpasses="+String.format("%.2f", passes)
						+"\tavgID="+String.format("%.4f", avgId);
				
//				ZMW.fixReadHeader(r, 0, trimmed);
			}
			
			if(sampleExact){
				if(sampleReadsExact){
					if(keepZMWsTogether){
						assert(readsRemaining>0) : readsRemaining;
						double prob=sampleReadsTarget/(double)(readsRemaining);
						if(randy.nextDouble()<prob){sampleReadsTarget-=reads.size();}
						else{reads.setDiscarded(true);}
						readsRemaining-=reads.size();
					}else{
						for(Read r : reads){
							if(r!=null && !r.discarded()){
								assert(readsRemaining>0) : readsRemaining;
								double prob=sampleReadsTarget/(double)(readsRemaining);
								if(randy.nextDouble()<prob){sampleReadsTarget--;}
								else{r.setDiscarded(true);}
							}
							readsRemaining--;
						}
					}
				}else if(sampleBasesExact){
					if(keepZMWsTogether){
						assert(basesRemaining>0) : basesRemaining;
						final long bases=reads.countBases();
						double prob=sampleBasesTarget/(double)(basesRemaining);
						if(randy.nextDouble()<prob){sampleBasesTarget-=bases;}
						else{reads.setDiscarded(true);}
						basesRemaining-=bases;
					}else{
						for(Read r : reads){
							if(r!=null && !r.discarded()){
								assert(basesRemaining>0) : basesRemaining;
								final int bases=r.length();
								double prob=sampleBasesTarget/(double)(basesRemaining);
								if(randy.nextDouble()<prob){sampleBasesTarget-=bases;}
								else{r.setDiscarded(true);}
								basesRemaining-=bases;
							}
						}
					}
				}else if(sampleZMWsExact){
					assert(ZMWsRemaining>0) : ZMWsRemaining;
					double prob=sampleZMWsTarget/(double)(ZMWsRemaining);
					if(randy.nextDouble()<prob){sampleZMWsTarget--;}
					else{reads.setDiscarded(true);}
					ZMWsRemaining--;
				}else{assert(false) : "No sampling mode.";}
			}
			
			if(trimReads || trimPolyA){
				int removed=0;
				for(int i=0; i<reads.size(); i++){
					Read r=reads.get(i);
					if(!r.discarded()){
						byte a=r.bases[0], b=r.bases[r.length()-1];
						int trimmed=0;
						if(!AminoAcid.isFullyDefined(a) || !AminoAcid.isFullyDefined(b)){
							trimmed+=trimRead(r, 80);
						}
						if(trimPolyA){
							int x=trimPolyA(r);
							trimmed+=x;
							if(x>0){r.setTrimmed(true);}
						}

						if(trimmed>0){
							basesTrimmedT+=trimmed;
							readsTrimmedT++;
							//TODO: technically this should track the lefty and right trim amounts
							//but it doesn't really matter as long as the sum is correct
							ZMW.fixReadHeader(r, 0, trimmed);
						}

						//TODO: Note again, removing intermediate reads messes up forward-rev ordering
						if(r.length()<minLengthAfterTrimming){//Discard short trash
							//						reads.set(i, null);
							r.setDiscarded(true);
							removed++;
						}
					}
				}
			}
			
			if(entropyCutoff>0){
				int bad=flagLowEntropyReads(reads, entropyCutoff, entropyLength, entropyFraction);
				if(bad>0){
					lowEntropyZMWsT++;
					lowEntropyReadsT+=bad;
					if(bad>=reads.size()){
						if(!reads.isEmpty()){outputReads(reads);}
						return; //No point in continuing...
					}
				}
			}
			
			if(reads.isEmpty()){return;}

			Read sample=null;//Read that will be searched for inverted repeat, typically median length
			Read shortest=null;//shortest read in the middle, or overall if no middle reads
			final int medianLength=medianLength(reads);
			int longReads=0;
			int sampleNum=0;
			
			if(reads.size()>=3){
				if(flagLongReads){longReads=flagLongReads(reads, medianLength);}
				for(int i=1; i<reads.size()-1; i++){
					Read r=reads.get(i);
					if(sample==null && r.length()==medianLength){
						sample=r;
						sampleNum=i;
					}
					if(shortest==null || r.length()<shortest.length()){shortest=r;}
				}
			}else{
				for(int i=0; i<reads.size(); i++){
					Read r=reads.get(i);
					if(sample==null || sample.length()<r.length()){
						sample=r;
						sampleNum=i;
					}
					if(shortest==null || r.length()<shortest.length()){shortest=r;}
				}
			}
			assert(sample!=null);
			
			if(keepShortReads){
				
			}else if(sample.discarded() || (longReads>0)){
				for(Read r : reads){
					r.setDiscarded(true);
				}
			}
			
			if(minLengthAfterTrimming>0){removeShortTrash(reads);}
			
			if(!reads.isEmpty()){outputReads(reads);}
		}
		
		private int removeShortTrash(ZMW reads) {
			int removed=0;
			for(int i=0; i<reads.size(); i++){
				Read r=reads.get(i);
				if(r.length()<minLengthAfterTrimming){//Discard short trash
					if(!r.discarded()){
						removed++;
						r.setDiscarded(true);
					}
				}
			}
			return removed;
		}
		
		private void fixReadHeader(Read r, int leftTrim, int rightTrim){
			leftTrim=Tools.max(0, leftTrim);
			rightTrim=Tools.max(0, rightTrim);
			if(leftTrim<1 && rightTrim<1){return;}
			final int idx=r.id.lastIndexOf('/');
			if(idx>0 && idx<r.id.length()-3){
				String prefix=r.id.substring(0, idx+1);
				String suffix=r.id.substring(idx+1);
				if(suffix.indexOf('_')>0){
					String coords=suffix, comment="";
					int tab=suffix.indexOf('\t');
					if(tab<0){tab=suffix.indexOf(' ');}
					if(tab>0){
						coords=coords.substring(0, tab);
						comment=coords.substring(tab);
					}
					String[] split=Tools.underscorePattern.split(coords);
					int left=Integer.parseInt(split[0]);
					int right=Integer.parseInt(split[1]);
					left+=leftTrim;
					right-=rightTrim;
					if(left>right){left=right;}
					
					if(right-left!=r.length()){right=left+r.length();}
//					System.err.println(r.length()+", "+(right-left));
					
					r.id=prefix+left+"_"+right+comment;
					final SamLine sl=r.samline;
					if(sl!=null){
						sl.qname=r.id;
						if(sl.optional!=null){
							for(int i=0; i<sl.optional.size(); i++){
								String s=sl.optional.get(i);
								if(s.startsWith("qe:i:")){
									s="qe:i:"+right;
									sl.optional.set(i, s);
								}else if(s.startsWith("qs:i:")){
									s="qs:i:"+left;
									sl.optional.set(i, s);
								}
							}
						}
					}
				}
			}
		}
		
		int trimRead(Read r, int lookahead){
			final byte[] bases=r.bases;
			
			int left=calcLeftTrim(bases, lookahead);
			int right=calcRightTrim(bases, lookahead);
			int trimmed=0;
			
			if(left+right>0){
//				System.err.println(r.length()+", "+left+", "+right+", "+r.id);
				trimmed=TrimRead.trimByAmount(r, left, right, 1, false);
//				System.err.println(r.length()+", "+left+", "+right+", "+r.id);
				fixReadHeader(r, left, right);
//				System.err.println(r.length()+", "+left+", "+right+", "+r.id);
			}
			
			if(r.samline!=null){
				r.samline.seq=r.bases;
				r.samline.qual=r.quality;
			}
			return trimmed;
		}
		
		int trimPolyA(Read r){
			final byte[] bases=r.bases;

			int left=Tools.max(PolymerTrimmer.testLeft(bases, 'A'), PolymerTrimmer.testLeft(bases, 'T'));
			int right=Tools.max(PolymerTrimmer.testRight(bases, 'A'), PolymerTrimmer.testRight(bases, 'T'));
			int trimmed=0;
			
			if(left+right>0){
//				System.err.println(r.length()+", "+left+", "+right+", "+r.id);
				trimmed=TrimRead.trimByAmount(r, left, right, 1, false);
//				System.err.println(r.length()+", "+left+", "+right+", "+r.id);
				fixReadHeader(r, left, right);
//				System.err.println(r.length()+", "+left+", "+right+", "+r.id);
			}
			
			if(r.samline!=null){
				r.samline.seq=r.bases;
				r.samline.qual=r.quality;
			}
			return trimmed;
		}
		
		final int calcLeftTrim(final byte[] bases, int lookahead){
			final int len=bases.length;
			int lastUndef=-1;
			for(int i=0, defined=0; i<len && defined<lookahead; i++){
				if(AminoAcid.isFullyDefined(bases[i])){
					defined++;
				}else{
					lastUndef=i;
					defined=0;
				}
			}
			return lastUndef+1;
		}
		
		final int calcRightTrim(final byte[] bases, int lookahead){
			final int len=bases.length;
			int lastUndef=len;
			for(int i=len-1, defined=0; i>=0 && defined<lookahead; i--){
				if(AminoAcid.isFullyDefined(bases[i])){
					defined++;
				}else{
					lastUndef=i;
					defined=0;
				}
			}
			return len-lastUndef-1;
		}
		
		private void outputReads(ZMW reads){
			final int size=reads.size();
			final ArrayList<Read> good=new ArrayList<Read>(size);
			final ArrayList<Read> bad=new ArrayList<Read>(size);
//			final ArrayList<Read> bad=(rosb==null ? null : new ArrayList<Read>(size));
			
			int discardedReads=0;
			int trimmedReads=0;
			boolean sendAllToDiscarded=false;
			
			//Check to see if any reads are discarded or ambiguous
			for(Read r : reads){
				if(r.discarded()){
					discardedReads++;
				}else if(r.trimmed()){
					trimmedReads++;
				}
			}
			
			//Unify flags on all reads
			if(keepZMWsTogether){
				if(discardedReads>0){sendAllToDiscarded=true;}
			}
			if(discardedReads>0){
				if(discardedReads==size){
					fullyDiscardedZMWsT++;
				}else{
					partiallyDiscardedZMWsT++;
				}
			}
			
			for(Read r : reads){
				if(r.discarded() || sendAllToDiscarded){
//					assert(false);
					if(bad!=null){bad.add(r);}
//					System.err.println("d:t="+r.tested()+",ad="+r.hasAdapter()+",ir="+r.invertedRepeat()+","+r.id+", reads="+reads.size());
				}else{
					if(good!=null){
						good.add(r);
					}
					readsOutT++;
					basesOutT+=r.length();
				}
			}

			if(good!=null && !good.isEmpty()){ZMWsOutT++;}
			if(rosg!=null && good!=null && !good.isEmpty()){rosg.add(good, 0);}
			if(rosb!=null && bad!=null && !bad.isEmpty()){rosb.add(bad, 0);}
		}
		
		/*--------------------------------------------------------------*/
		/*----------------      Inner Class Fields      ----------------*/
		/*--------------------------------------------------------------*/
		
		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		/** Number of ZMWs processed by this thread */
		protected long ZMWsT=0;
		
		protected long basesTrimmedT=0;
		protected long readsTrimmedT=0;
		
		/** Number of reads retained by this thread */
		protected long readsOutT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;
		/** Number of ZMWs retained by this thread */
		protected long ZMWsOutT=0;

		protected long partiallyDiscardedZMWsT=0;
		protected long fullyDiscardedZMWsT=0;
		
		protected long lowEntropyZMWsT=0;
		protected long lowEntropyReadsT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared output stream */
		private final ConcurrentReadOutputStream rosg;
		/** Shared output stream for bad reads */
		private final ConcurrentReadOutputStream rosb;
		/** Thread ID */
		final int tid;
		
		final ZMWStreamer zstream;
		final EntropyTracker eTracker;
		Random randy;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary input file path */
	private String in1=null;

	/** Primary output file path */
	private String outg=null;
	/** Bad output file path */
	private String outb=null;
	
	/** Stats (screen) output file path */
	private String outstats=null;
	
	/** Subread count histogram */
	private String schist=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;
	
	private float longReadMult=1.5f;
	
	/** For grading synthetic data */
	private boolean parseCustom;
	
	/** Input reads are CCS (full pass) */
	private boolean CCSInput;
	
	//Note: These flags are very similar... they need to be better-defined or merged.
	private boolean keepZMWsTogether=false;
	private boolean keepShortReads=true;
	
	private boolean subsampleFromEnds=false;
	
	private int format=FORMAT_TEXT;
	private static final int FORMAT_TEXT=1, FORMAT_JSON=2;
	
	/*--------------------------------------------------------------*/
	
	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of reads retained */
	protected long readsOut=0;
	/** Number of bases retained */
	protected long basesOut=0;
	/** Number of ZMWs retained */
	protected long ZMWsOut=0;

	protected long partiallyDiscardedZMWs=0;
	protected long fullyDiscardedZMWs=0;
	
	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	protected long basesTrimmed=0;
	protected long readsTrimmed=0;
	
	protected long lowEntropyZMWs=0;
	protected long lowEntropyReads=0;
	
	/** Total ZMWs observed */
	protected long ZMWs=0;
	
	/** Histogram */
	protected final long[] subreadCounts=new long[101];
	
	private boolean flagLongReads=false;
	private boolean trimReads=false;
	private int minLengthAfterTrimming=0;
	
	boolean trimPolyA=false;
	
	/*--------------------------------------------------------------*/
	/*----------------      Subsampling Fields      ----------------*/
	/*--------------------------------------------------------------*/
	
	//If a target is chosen, these will be initialized with the original counts
	long initialReads=0;
	long initialZMWs=0;
	long initialBases=0;
	
	long readsRemaining=0;
	long ZMWsRemaining=0;
	long basesRemaining=0;
	
	float samplerate=1.0f;
	long sampleReadsTarget=-1;
	long sampleBasesTarget=-1;
	long sampleZMWsTarget=-1;

	final boolean sampleReadsExact;
	final boolean sampleBasesExact;
	final boolean sampleZMWsExact;
	final boolean sampleExact;
	
	boolean keepBestPass=false;
	boolean keepLongestPass=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Subset Fields         ----------------*/
	/*--------------------------------------------------------------*/

	String whitelist;
	String blacklist;
	
	IntHashSet whiteSet;
	IntHashSet blackSet;
	
	/*--------------------------------------------------------------*/
	/*----------------        Entropy Fields        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Minimum entropy to be considered "complex", on a scale of 0-1 */
	float entropyCutoff=-1; //suggested 0.55
	/** Minimum length of a low-entropy block to fail a read */
	int entropyLength=350;
	/** Minimum length of a low-entropy block as a fraction of read length;
	 * the minimum of the two will be used */
	float entropyFraction=0.5f;
	
	float maxMonomerFraction=0.74f; //Suggested...  0.74
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;
	
	/** Primary output file */
	private final FileFormat ffoutg;
	
	/** Bad output file */
	private final FileFormat ffoutb;

	/** Stats output file */
	private final FileFormat ffstats;

	/** Subread count histogram */
	private final FileFormat ffschist;
	
	private final int threads;
	
	private long seed=-1;
	private long maxZMWs=-1;
	
	private int shredLength=500;
	private int overlap=10;//Helps make the shreds concur at their borders.
	private float minShredIdentity=0.6f;
	private boolean makeCCS=false;
	private boolean findOrientation=false;
	
	private float minPasses=0;
	private int minSubreads=0;
	
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
