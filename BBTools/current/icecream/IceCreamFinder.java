package icecream;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import aligner.AlignmentResult;
import aligner.FlatAligner;
import aligner.MultiStateAligner9PacBioAdapter2;
import aligner.SingleStateAlignerPacBioAdapter;
import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import json.JsonObject;
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
import tracker.EntropyTracker;
import tracker.ReadStats;

/**
 * Detects inverted repeats in PacBio reads.
 * Attempts to determine whether reads are chimeric
 * due to a missing adapter.
 * 
 * @author Brian Bushnell
 * @date June 5, 2019
 *
 */
public final class IceCreamFinder {
	
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
		IceCreamFinder x=new IceCreamFinder(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public IceCreamFinder(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Parser.setQuality(33);
		
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
		
		//Determine how many threads may be used
		threads=Shared.threads();
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffoutg=FileFormat.testOutput(outg, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffouta=FileFormat.testOutput(outa, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffoutb=FileFormat.testOutput(outb, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffoutj=FileFormat.testOutput(outj, FileFormat.FASTA, extout, true, overwrite, append, false);
		ffstats=FileFormat.testOutput(outstats, FileFormat.TXT, null, false, overwrite, append, false);
		ffasrhist=FileFormat.testOutput(asrhist, FileFormat.TXT, null, false, overwrite, append, false);
		ffirsrhist=FileFormat.testOutput(irsrhist, FileFormat.TXT, null, false, overwrite, append, false);
		
		if(!setAmbig && ffouta!=null){
			sendAmbigToGood=sendAmbigToBad=false;
		}
		
		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		
		final int alen=(adapter==null ? 20 : adapter.length);
		stride=(int)(alen*1.9f);
		window=(int)(alen*3.8f+10);
		ALIGN_ROWS=alen+1;
		ALIGN_COLUMNS=window+2;
		
		maxSwScore=aligner.MultiStateAligner9PacBioAdapter.maxQuality(alen);
		minSwScore=(int)(maxSwScore*adapterRatio2);
		minSwScoreSuspect=(int)(maxSwScore*Tools.min(adapterRatio2*suspectRatio, adapterRatio2-((1-suspectRatio)*.2f)));
		maxImperfectSwScore=aligner.MultiStateAligner9PacBioAdapter.maxImperfectScore(alen);
		suspectMidpoint=(minSwScoreSuspect+minSwScore)/2;
		if(adapter==null){alignAdapter=false;}
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
			}else if(a.equals("icecreamonly") || a.equals("ico")){
				filterIceCreamOnly=Parse.parseBoolean(b);
			}else if(a.equals("keepshortreads") || a.equals("ksr")){
				keepShortReads=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("keepzmwstogether") || a.equals("kzt") || a.equals("keepreadstogether") || a.equals("krt")){
				keepZMWsTogether=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("samplerate")){
				float f=Float.parseFloat(b);
				assert(false) : "TODO"; //TODO
			}else if(a.equals("modifyheader") || a.equals("modifyheaders") || a.equals("changeheader") || a.equals("changeheaders")){
				modifyHeader=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("ccs")){
				CCS=Parse.parseBoolean(b);
			}else if(a.equals("npad")){
				npad=Integer.parseInt(b);
			}else if(a.equals("minlength") || a.equals("minlen")){
				minLengthAfterTrimming=Integer.parseInt(b);
			}else if(a.equals("realign")){
				realign=Parse.parseBoolean(b);
			}else if(a.equals("aligninverse") || a.equals("alignrc") || a.equals("findicecream")){
				alignRC=Parse.parseBoolean(b);
			}else if(a.equals("alignadapter")){
				alignAdapter=Parse.parseBoolean(b);
			}else if(a.equals("timeless")){
				timeless=Parse.parseBoolean(b);
			}else if(a.equals("flaglongreads")){
				flagLongReads=Parse.parseBoolean(b);
			}else if(a.equals("trimreads") || a.equals("trim")){
				trimReads=Parse.parseBoolean(b);
			}else if(a.equals("adapter")){
				adapter=b==null ? null : b.getBytes();
			}else if(a.equals("targetqlen") || a.equals("qlen")){
				targetQlen=Integer.parseInt(b);
			}else if(a.equals("maxqlenfraction") || a.equals("maxfraction") || a.equals("qlenfraction")){
				maxQlenFraction=Float.parseFloat(b);
			}else if(a.equals("junctionfraction") || a.equals("shortfraction")){
				minJunctionFraction=Float.parseFloat(b);
			}else if(a.equals("minratio1") || a.equals("ratio1") || a.equals("id1")){
				minRatio1=Float.parseFloat(b);
			}else if(a.equals("minratio2") || a.equals("ratio2") || a.equals("id2")){
				minRatio2=Float.parseFloat(b);
			}else if(a.equals("minratio") || a.equals("ratio") || a.equals("id")){
				minRatio1=minRatio2=Float.parseFloat(b);
			}else if(a.equals("adapterratio") || a.equals("adapterratio1") || a.equals("ratior") || a.equals("ratior1")){
				adapterRatio=Float.parseFloat(b);
			}else if(a.equals("adapterratio2") || a.equals("ratior2")){
				adapterRatio2=Float.parseFloat(b);
			}else if(a.equals("suspectratio")){
				suspectRatio=Float.parseFloat(b);
			}else if(a.equals("minqlen")){
				minQlen=Integer.parseInt(b);
			}else if(a.equals("minscore")){
				minScore=Integer.parseInt(b);
			}else if(a.equals("parsecustom")){
				parseCustom=Parse.parseBoolean(b);
			}else if(a.equals("printtiming") || a.equals("extended")){
				printTiming=Parse.parseBoolean(b);
			}else if(a.equals("queuelen") || a.equals("qlen")){
				queuelen=Integer.parseInt(b);
			}else if(a.equals("outg") || a.equals("outgood")){
				outg=b;
			}else if(a.equals("outa") || a.equals("outambig")){
				outa=b;
			}else if(a.equals("outb") || a.equals("outbad")){
				outb=b;
			}else if(a.equals("outj") || a.equals("outjunctions") || a.equals("junctions")){
				outj=b;
			}else if(a.equals("outs") || a.equals("outstats") || a.equals("stats")){
				outstats=b;
			}else if(a.equals("asrhist") || a.equals("ahist")){
				asrhist=b;
			}else if(a.equals("irsrhist") || a.equals("irhist")){
				irsrhist=b;
			}else if(a.equals("ambig")){
				sendAmbigToGood=sendAmbigToBad=false;
				if(b!=null){
					String[] split2=b.split(",");
					for(String s2 : split2){
						if(s2.equalsIgnoreCase("good")){sendAmbigToGood=true;}
						else if(s2.equalsIgnoreCase("bad") || s2.equalsIgnoreCase("toss")){sendAmbigToBad=true;}
						else if(s2.equalsIgnoreCase("ambig")){}
						else{assert(false) : "Bad argument: '"+s2+"' in '"+arg+"'; should be good or bad";}
					}
				}
				setAmbig=true;
			}else if(a.equalsIgnoreCase("trimpolya")){//Parse standard flags in the parser
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
		
		return parser;
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, outg, outa, outb, outj, outstats, asrhist, irsrhist)){
			outstream.println((outg==null)+", "+(outb==null)+", "+outg+", "+outa+", "+outb+", "+outj+", "+outstats);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+outg+", "+outa+", "+outb+", "+outj+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, outg, outa, outb, outj, outstats, asrhist, irsrhist)){
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
		
		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;
		
		//Create a read input stream
		ZMWStreamer zstream=new ZMWStreamer(ffin1, Shared.threads(), maxReads, -1);
		
		//Optionally create read output streams
		final ConcurrentReadOutputStream rosg=makeCros(ffoutg);
		final ConcurrentReadOutputStream rosa=makeCros(ffouta);
		final ConcurrentReadOutputStream rosb=makeCros(ffoutb);
		final ConcurrentReadOutputStream rosj=makeCros(ffoutj);
		
		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;
		junctionsOut=0;
		
		//Process the reads in separate threads
		spawnThreads(zstream, rosg, rosa, rosb, rosj);
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(null, rosg, rosa, rosb, rosj);
		
		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;
		
		writeScoreRatioHistogram(ffasrhist, adapterScores);
		writeScoreRatioHistogram(ffirsrhist, repeatScores);
		
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
		
		bb.appendln(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		bb.appendln(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		
		long readsFiltered=readsProcessed-readsOut;
		bb.appendln(Tools.numberPercent("Reads Filtered:", readsFiltered, readsFiltered*100.0/(readsProcessed), 3, 8));
		if(trimReads || trimPolyA){
			bb.appendln(Tools.numberPercent("Reads Trimmed:", readsTrimmed, readsTrimmed*100.0/(readsProcessed), 3, 8));
			bb.appendln(Tools.numberPercent("Bases Trimmed:", basesTrimmed, basesTrimmed*100.0/(basesProcessed), 3, 8));
		}
		bb.appendln(Tools.number("Total ZMWs:", ZMWs, 8));
		bb.appendln(Tools.numberPercent("Bad ZMWs:", iceCreamZMWs, iceCreamZMWs*100.0/(ZMWs), 3, 8));
		bb.appendln(Tools.numberPercent("Absent Adapter:", missingAdapterZMWs, missingAdapterZMWs*100.0/(ZMWs), 3, 8));
		bb.appendln(Tools.numberPercent("Hidden Adapter:", hiddenAdapterZMWs, hiddenAdapterZMWs*100.0/(ZMWs), 3, 8));
		bb.appendln(Tools.numberPercent("Ambiguous IR:", ambiguousZMWs, ambiguousZMWs*100.0/(ZMWs), 3, 8));
//		bb.appendln(Tools.numberPercent("Low Entropy:", lowEntropyReads, lowEntropyReads*100.0/(readsProcessed), 3, 8));
		bb.appendln(Tools.numberPercent("Low Entropy:", lowEntropyZMWs, lowEntropyZMWs*100.0/(ZMWs), 3, 8));
		
		bb.appendln(Tools.number("Avg Score Ratio:", iceCreamRatio/ratiosCounted, 3, 8));
		bb.appendln(Tools.number("Score Cutoff:", minRatio2, 3, 8));
		
		if(printTiming){
			bb.appendln("Iterations:         "+alignmentIters/1000000+"m");
			bb.appendln("Short Iterations:   "+alignmentItersShort/1000000+"m");
			bb.appendln("Elapsed:            "+elapsed/1000000+"ms");
			bb.appendln("Elapsed Short:      "+elapsedShort/1000000+"ms");
		}
		
		if(parseCustom){
			{
				bb.appendln("\nReads:");
				bb.appendln(Tools.numberPercent("True Positive:", truePositiveReads, truePositiveReads*100.0/(readsProcessed), 2, 8));
				bb.appendln(Tools.numberPercent("True Negative:", trueNegativeReads, trueNegativeReads*100.0/(readsProcessed), 2, 8));
				bb.appendln(Tools.numberPercent("False Positive:", falsePositiveReads, falsePositiveReads*100.0/(readsProcessed), 2, 8));
				bb.appendln(Tools.numberPercent("False Negative:", falseNegativeReads, falseNegativeReads*100.0/(readsProcessed), 2, 8));
				bb.appendln(Tools.numberPercent("Ambiguous:", ambiguousReads, ambiguousReads*100.0/(readsProcessed), 2, 8));

				double snr=(truePositiveReads+trueNegativeReads+ambiguousReads)/Tools.max(1, falsePositiveReads+falseNegativeReads+ambiguousReads);
				snr=10*Math.log10(snr);
				bb.appendln(Tools.number("SNR:", snr, 2, 8));
			}
			{
				bb.appendln("\nZMWs:");
				bb.appendln(Tools.numberPercent("True Positive:", truePositiveZMWs, truePositiveZMWs*100.0/(ZMWs), 2, 8));
				bb.appendln(Tools.numberPercent("True Negative:", trueNegativeZMWs, trueNegativeZMWs*100.0/(ZMWs), 2, 8));
				bb.appendln(Tools.numberPercent("False Positive:", falsePositiveZMWs, falsePositiveZMWs*100.0/(ZMWs), 2, 8));
				bb.appendln(Tools.numberPercent("False Negative:", falseNegativeZMWs, falseNegativeZMWs*100.0/(ZMWs), 2, 8));
				bb.appendln(Tools.numberPercent("Ambiguous:", ambiguousZMWs, ambiguousZMWs*100.0/(ZMWs), 2, 8));

				double snr=(truePositiveZMWs+trueNegativeZMWs+ambiguousReads)/Tools.max(1, falsePositiveZMWs+falseNegativeZMWs+ambiguousZMWs);
				snr=10*Math.log10(snr);
				bb.appendln(Tools.number("SNR:", snr, 2, 8));
			}
		}
		return bb;
	}
	
	private JsonObject toJson(Timer t){
		JsonObject jo=new JsonObject();
		long readsFiltered=readsProcessed-readsOut;
		
		jo.add("Time", t.timeInSeconds());
		jo.add("Reads_Processed", readsProcessed);
		jo.add("Bases_Processed", basesProcessed);
		jo.add("Reads_Out", readsOut);
		jo.add("Bases_Out", basesOut);
		jo.add("Reads_Filtered", readsFiltered);
		jo.add("Reads_Filtered_Pct", readsFiltered*100.0/(readsProcessed));
		if(trimReads){
			jo.add("Reads_Trimmed", readsTrimmed);
			jo.add("Reads_Trimmed_Pct", readsTrimmed*100.0/(readsProcessed));
			jo.add("Bases_Trimmed", basesTrimmed);
			jo.add("Bases_Trimmed_Pct", basesTrimmed*100.0/(basesProcessed));
		}
		jo.add("Total_ZMWs", ZMWs);
		jo.add("Bad_ZMWs", iceCreamZMWs);
		jo.add("Bad_ZMWs_Pct", iceCreamZMWs*100.0/(ZMWs));
		jo.add("Absent_Adapter", missingAdapterZMWs);
		jo.add("Absent_Adapter_Pct", missingAdapterZMWs*100.0/(ZMWs));
		jo.add("Hidden_Adapter", hiddenAdapterZMWs);
		jo.add("Hidden_Adapter_Pct", hiddenAdapterZMWs*100.0/(ZMWs));
//		jo.add("Low_Entropy", lowEntropyReads);
//		jo.add("Low_Entropy_Pct", lowEntropyReads*100.0/(readsProcessed));
		jo.add("Low_Entropy", lowEntropyZMWs);
		jo.add("Low_Entropy_Pct", lowEntropyZMWs*100.0/(ZMWs));
		jo.add("Ambiguous_Inverted_Repeat", ambiguousZMWs);
		jo.add("Ambiguous_Inverted_Repeat_Pct", ambiguousZMWs*100.0/(ZMWs));
		jo.add("Avg_Score_Ratio_IR", iceCreamRatio/ratiosCounted);
		jo.add("Score_Cutoff_IR", minRatio2);
		
		jo.add("Alignment_Iterations_IR", alignmentIters);
		jo.add("Short_Alignment_Iterations_IR", alignmentItersShort);
		
		if(parseCustom){
			{
				double snr=(truePositiveReads+trueNegativeReads+ambiguousReads)/Tools.max(1, falsePositiveReads+falseNegativeReads+ambiguousReads);
				snr=10*Math.log10(snr);
				jo.add("TP_Reads", truePositiveReads);
				jo.add("TN_Reads", trueNegativeReads);
				jo.add("FP_Reads", falsePositiveReads);
				jo.add("FN_Reads", falseNegativeReads);
				jo.add("AM_Reads", ambiguousReads);

				jo.add("TP_Reads_Pct", truePositiveReads*100.0/(readsProcessed));
				jo.add("TN_Reads_Pct", trueNegativeReads*100.0/(readsProcessed));
				jo.add("FP_Reads_Pct", falsePositiveReads*100.0/(readsProcessed));
				jo.add("FN_Reads_Pct", falseNegativeReads*100.0/(readsProcessed));
				jo.add("AM_Reads_Pct", ambiguousReads*100.0/(readsProcessed));
				
				jo.add("SNR_Reads", snr);
			}
			{
				double snr=(truePositiveZMWs+trueNegativeZMWs+ambiguousZMWs)/Tools.max(1, falsePositiveZMWs+falseNegativeZMWs+ambiguousZMWs);
				snr=10*Math.log10(snr);
				jo.add("TP_ZMWs", truePositiveZMWs);
				jo.add("TN_ZMWs", trueNegativeZMWs);
				jo.add("FP_ZMWs", falsePositiveZMWs);
				jo.add("FN_ZMWs", falseNegativeZMWs);
				jo.add("AM_ZMWs", ambiguousZMWs);

				jo.add("TP_ZMWs_Pct", truePositiveZMWs*100.0/(ZMWs));
				jo.add("TN_ZMWs_Pct", trueNegativeZMWs*100.0/(ZMWs));
				jo.add("FP_ZMWs_Pct", falsePositiveZMWs*100.0/(ZMWs));
				jo.add("FN_ZMWs_Pct", falseNegativeZMWs*100.0/(ZMWs));
				jo.add("AM_ZMWs_Pct", ambiguousZMWs*100.0/(ZMWs));
				
				jo.add("SNR_ZMWs", snr);
			}
		}
		return jo;
	}
	
	private static void writeScoreRatioHistogram(FileFormat ff, long[] hist){
		if(ff==null){return;}
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		final float mult=1.0f/(hist.length-1);

		bsw.print("#Counted\t").println(shared.Vector.sum(hist));
		bsw.print("#Mean\t").println(Tools.averageHistogram(hist)*mult, 3);
		bsw.print("#Median\t").println(Tools.medianHistogram(hist)*mult, 3);
		bsw.print("#Mode\t").println(Tools.calcModeHistogram(hist)*mult, 3);
		bsw.print("#STDev\t").println(Tools.standardDeviationHistogram(hist)*mult, 3);
		bsw.print("#Value\tOccurances\n");
		
		for(int i=0; i<hist.length; i++){
			bsw.print(i*mult, 3).tab().println(hist[i]);
		}
		bsw.poisonAndWait();
	}
	
	private ConcurrentReadOutputStream makeCros(FileFormat ff){
		if(ff==null){return null;}

		//Select output buffer size based on whether it needs to be ordered
		final int buff=16;

		final ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ff, null, buff, null, ff.samOrBam() && ffin1.samOrBam());
		ros.start(); //Start the stream
		return ros;
	}
	
	/** Spawn process threads */
	private void spawnThreads(final ZMWStreamer zstream,
			final ConcurrentReadOutputStream rosg, final ConcurrentReadOutputStream rosa,
			final ConcurrentReadOutputStream rosb, final ConcurrentReadOutputStream rosj){
		
		//Do anything necessary prior to processing
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(zstream, rosg, rosa, rosb, rosj, i));
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

			truePositiveReads+=pt.truePositiveReadsT;
			trueNegativeReads+=pt.trueNegativeReadsT;
			falsePositiveReads+=pt.falsePositiveReadsT;
			falseNegativeReads+=pt.falseNegativeReadsT;
			ambiguousReads+=pt.ambiguousReadsT;

			truePositiveZMWs+=pt.truePositiveZMWsT;
			trueNegativeZMWs+=pt.trueNegativeZMWsT;
			falsePositiveZMWs+=pt.falsePositiveZMWsT;
			falseNegativeZMWs+=pt.falseNegativeZMWsT;
			ambiguousZMWs+=pt.ambiguousZMWsT;
			
			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			junctionsOut+=pt.junctionsOutT;
			
			alignmentIters+=pt.ica.iters();
			alignmentItersShort+=pt.ica.itersShort();
			elapsed+=pt.elapsedT;
			elapsedShort+=pt.elapsedShortT;
			
			iceCreamReads+=pt.iceCreamReadsT;
			iceCreamBases+=pt.iceCreamBasesT;
			iceCreamZMWs+=pt.iceCreamZMWsT;
			iceCreamRatio+=pt.iceCreamRatioT;
			ratiosCounted+=pt.ratiosCountedT;
			missingAdapterZMWs+=pt.missingAdapterZMWsT;
			hiddenAdapterZMWs+=pt.hiddenAdapterZMWsT;
			lowEntropyZMWs+=pt.lowEntropyZMWsT;
			lowEntropyReads+=pt.lowEntropyReadsT;
			
			basesTrimmed+=pt.basesTrimmedT;
			readsTrimmed+=pt.readsTrimmedT;
			
			for(int i=0; i<adapterScores.length; i++){
				adapterScores[i]+=pt.adapterScoresT[i];
				repeatScores[i]+=pt.repeatScoresT[i];
			}
			
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Processes reads. */
	private class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ZMWStreamer zstream_, 
				final ConcurrentReadOutputStream ros_, final ConcurrentReadOutputStream rosa_, 
				final ConcurrentReadOutputStream rosb_, final ConcurrentReadOutputStream rosj_, final int tid_){
			zstream=zstream_;
			rosg=ros_;
			rosa=rosa_;
			rosb=rosb_;
			rosj=rosj_;
			tid=tid_;
			
			Arrays.fill(tipBufferLeft, (byte)'N');
			Arrays.fill(tipBufferRight, (byte)'N');
			
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
		
		/** Each list is presumed to be all reads from a ZMW, in order */
		void processList(final ZMW reads){
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
			final boolean fullPass=CCS || reads.size()>=3;
			
			if(trimReads || trimPolyA){
				int removed=0;
				for(int i=0; i<reads.size(); i++){
					Read r=reads.get(i);
					byte a=r.bases[0], b=r.bases[r.length()-1];
					int trimmed=0;
					if(!AminoAcid.isFullyDefined(a) || !AminoAcid.isFullyDefined(b)){
						trimmed+=trimRead(r, 80);
					}
					
					if(trimReads && adapter!=null){
						int leftAdapterBases=alignLeftTipAdapter(r);
						int rightAdapterBases=alignRightTipAdapter(r);
						if(leftAdapterBases+rightAdapterBases>0){
							trimmed+=trimRead(r, adapterTipLen);
							r.setTrimmed(true);
						}
					}
					if(trimPolyA){
						int x=trimPolyA(r);
						trimmed+=x;
						if(x>0){r.setTrimmed(true);}
					}
					
					if(trimmed>0){
						basesTrimmedT+=trimmed;
						readsTrimmedT++;
					}
					
					//TODO: Note again, removing intermediate reads messes up forward-rev ordering
					if(r.length()<minLengthAfterTrimming){//Discard short trash
						reads.set(i, null);
						removed++;
					}
				}
				if(removed>0){
					Tools.condenseStrict(reads);
				}
			}
			
			if(entropyCutoff>0){
				int bad=flagLowEntropyReads(reads, entropyCutoff, entropyLength, entropyFraction);
				if(bad>0){
					lowEntropyZMWsT++;
					lowEntropyReadsT+=bad;
					if(bad>=reads.size()){
						if(!reads.isEmpty()){outputReads(reads, null);}
						return; //No point in continuing...
					}
				}
			}
			
			if(reads.isEmpty()){return;}

			Read sample=null;//Read that will be searched for inverted repeat, typically median length
			Read shortest=null;//shortest read in the middle, or overall if no middle reads
			final int medianLength=reads.medianLength(true);
			boolean foundInvertedRepeat=false;
			int longReads=0;
			int adapterReads=0;
			int maxAdapters=0;
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
			final AlignmentResult ar=align(sample, fullPass, reads.size(), sampleNum);
			
			if(ar!=null){
				foundInvertedRepeat=true;
				sample.setInvertedRepeat(true);
				if(ar.icecream || !filterIceCreamOnly){
					sample.setDiscarded(true);
				}else if(ar.ambiguous){
					sample.setAmbiguous(true);
				}
			}
			
			if(alignAdapter){
				double mult=foundInvertedRepeat ? 0.9 : 1.0;
				if(needsAdapterTest(sample)){
					int x=lookForAdapter(sample, mult);
					adapterReads+=(x>0 ? 1 : 0);
					maxAdapters=Tools.max(x, maxAdapters);
				}
				
				if(reads.size()>2){
					Read a=reads.get(0), b=reads.get(reads.size()-1);
					Read r=a.length()>b.length() ? a : b;
					if(needsAdapterTest(r)){
						int x=lookForAdapter(r, mult);
						adapterReads+=(x>0 ? 1 : 0);
						maxAdapters=Tools.max(x, maxAdapters);
					}
				}
				
				for(Read r : reads){
					if((r.length()>=shortReadMult*shortest.length() || adapterReads>0 || longReads>0 || foundInvertedRepeat) 
							&& needsAdapterTest(r)){
						int x=lookForAdapter(r, mult);
						adapterReads+=(x>0 ? 1 : 0);
						maxAdapters=Tools.max(x, maxAdapters);
					}
				}
			}
			
			if(ar!=null && ar.icecream){
				iceCreamRatioT+=ar.ratio;
				ratiosCountedT++;
				int idx=(int)(ar.ratio*200);
				repeatScoresT[idx]++;
				if(longReads+adapterReads==0){missingAdapterZMWsT++;}
			}
			if(longReads+adapterReads>0){
				hiddenAdapterZMWsT++;
			}
			
			if(keepShortReads && maxAdapters<2){
				if(foundInvertedRepeat && !sample.hasAdapter()){
					if(reads.size()>2){
						for(int i=1; i<reads.size()-1; i++){
							reads.get(i).setDiscarded(true);
						}
						Read r=reads.get(0);
						if(r.length()>=veryShortMult*medianLength){r.setDiscarded(true);}
						r=reads.get(reads.size()-1);
						if(r.length()>=veryShortMult*medianLength){r.setDiscarded(true);}
					}else if(reads.size()==2){
						for(Read r : reads){
							if(r.length()>=veryShortMult*sample.length()){
								if(ar.icecream){
									r.setDiscarded(true);
								}else if(ar.ambiguous){
									r.setAmbiguous(true);
								}
							}
						}
					}
				}
			}else if(sample.discarded() || (longReads+adapterReads>0)){
				for(Read r : reads){
					r.setDiscarded(true);
				}
			}
			
			ArrayList<Read> junctions=null;
			if(ar!=null){
				if(rosj!=null && !sample.hasAdapter()){
					Read r=ar.alignedRead;
					int width=Tools.min(200, ar.junctionLoc, r.length()-ar.junctionLoc);
					int a=ar.junctionLoc-width, b=ar.junctionLoc+width;
					byte[] bases=Arrays.copyOfRange(r.bases, a, b);
					Read junction=new Read(bases, null, r.id+"\tjunction:"+a+"-"+b, r.numericID);
					junctions=new ArrayList<Read>(1);
					junctions.add(junction);
					junctionsOutT++;
				}
				if(modifyHeader){
					sample.id=sample.id+"\tratio="+ar.ratio+"\tjunction="+ar.junctionLoc+
							"\tIR="+sample.invertedRepeat()+"\tAD="+sample.hasAdapter()+"\tFP="+fullPass+"\tsubreads="+reads.size();
				}
			}
			
			removeShortTrash(reads);
			
			if(!reads.isEmpty()){outputReads(reads, junctions);}
		}
		
		//TODO: Now that I think about it.  The order of the reads is important.
		//Since they go forward-rev-forward-rev it's imprudent to discard inner reads.
		private void removeShortTrash(ZMW reads) {
			int removed=0;
			for(int i=0; i<reads.size(); i++){
				Read r=reads.get(i);
				if(r.length()<minLengthAfterTrimming){//Discard short trash
					reads.set(i, null);
					removed++;
				}
			}
			if(removed>0){Tools.condenseStrict(reads);}
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
				ZMW.fixReadHeader(r, left, right);
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
				ZMW.fixReadHeader(r, left, right);
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
		
		final int alignLeftTipAdapter(Read r){
			assert(adapter.length<adapterTipLen); //Time to increase adapterTipLen
			if(r.length()<adapterTipLen){return 0;}
			final byte[] array=tipBufferLeft;
			
			for(int i=adapterTipPad, j=0; i<array.length; i++, j++){array[i]=r.bases[j];}
			int[] rvec=ssa.fillAndScoreLimited(adapter, array, 0, array.length, minSwScore);

			if(rvec==null || rvec[0]<minSwScore){return 0;}
			final int score=rvec[0];
			final int start=Tools.max(0, rvec[1]-adapterTipPad);
			final int stop=rvec[2]-adapterTipPad;
			for(int i=start; i<=stop; i++){r.bases[i]='X';}
			return stop-start+1;
		}
		
		final int alignRightTipAdapter(Read r){
			final byte[] bases=r.bases;
			assert(adapter.length<adapterTipLen); //Time to increase adapterTipLen
			if(r.length()<adapterTipLen){return 0;}
			final byte[] array=tipBufferRight;
			
			for(int i=0, j=bases.length-adapterTipLen; i<adapterTipLen; i++, j++){array[i]=bases[j];}
			int[] rvec=ssa.fillAndScoreLimited(adapter, array, 0, array.length, minSwScore);

			if(rvec==null || rvec[0]<minSwScore){return 0;}
			final int score=rvec[0];
			final int start=Tools.max(0, rvec[1]-adapterTipPad);
			final int stop=rvec[2]-adapterTipPad;
			for(int i=start; i<=stop; i++){r.bases[i]='X';}
			return stop-start+1;
		}
		
		boolean needsAdapterTest(Read r){
			if(r.tested() || r.hasAdapter()){return false;}
			if(adapterRatio<=0 || r.discarded()){return true;}
			AlignmentResult aa=fla.alignForwardShort(adapter, r.bases, 0, r.bases.length-1, adapterRatio);
			return aa!=null;
		}
		
		private void outputReads(ZMW reads, ArrayList<Read> junctions){
			final int size=reads.size();
			final ArrayList<Read> good=(rosg==null ? null : new ArrayList<Read>(size));
			final ArrayList<Read> ambig=(rosa==null ? null : new ArrayList<Read>(size));
			final ArrayList<Read> bad=(rosb==null ? null : new ArrayList<Read>(size));
			
			int discardedReads=0;
			int ambigReads=0;
			int trimmedReads=0;
			boolean sendAllToDiscarded=false;
			boolean sendAllToAmbiguous=false;
			
			//Check to see if any reads are discarded or ambiguous
			for(Read r : reads){
				if(r.discarded()){
					discardedReads++;
				}else if(r.ambiguous()){
					ambigReads++;
				}else if(r.trimmed()){
					trimmedReads++;
				}
			}
			
			//Unify flags on all reads
			if(keepZMWsTogether){
				if(discardedReads>0){sendAllToDiscarded=true;}
				else if(ambigReads>0){sendAllToAmbiguous=true;}
			}
//			if(discardedReads>0 || ambigReads>0){System.err.println("\nd="+discardedReads+", a="+ambigReads);}
			if(discardedReads>0){iceCreamZMWsT++;}
			
			int trueIceCreamReads=0;
			for(Read r : reads){
				boolean trueIceCream=(parseCustom ? ReadBuilder.isIceCream(r.id) : false);
				trueIceCreamReads+=(trueIceCream ? 1 : 0);
				if(r.discarded() || sendAllToDiscarded){
//					assert(false);
					if(bad!=null){bad.add(r);}
					if(trueIceCream){truePositiveReadsT++;}
					else{falsePositiveReadsT++;}
//					System.err.println("d:t="+r.tested()+",ad="+r.hasAdapter()+",ir="+r.invertedRepeat()+","+r.id+", reads="+reads.size());
				}else if(r.ambiguous() || sendAllToAmbiguous){
					if(ambig!=null){ambig.add(r);}
					if(sendAmbigToGood){
						readsOutT++;
						basesOutT+=r.length();
						if(good!=null) {good.add(r);}
					}
					if(sendAmbigToBad && bad!=null){bad.add(r);}
					ambiguousReadsT++;
//					System.err.println("a*:t="+r.tested()+",ad="+r.hasAdapter()+",ir="+r.invertedRepeat()+","+r.id+", reads="+reads.size());
				}else{
					if(good!=null){
						good.add(r);
					}
					readsOutT++;
					basesOutT+=r.length();
					if(trueIceCream && !r.trimmed()){falseNegativeReadsT++;}
					else{trueNegativeReadsT++;}
//					if(discardedReads>0 || ambigReads>0){System.err.println("g*:t="+r.tested()+",ad="+r.hasAdapter()+",ir="+r.invertedRepeat()+","+r.id+", reads="+reads.size());}
				}
			}
			
			if(trueIceCreamReads>0){
				if(discardedReads>0 || trimmedReads>0){
					truePositiveZMWsT++;
				}else if(ambigReads>0){
					ambiguousZMWs++;
				}else{
					falseNegativeZMWsT++;
//					StringBuilder sb=new StringBuilder();
//					for(Read r : reads) {sb.append("\n").append(r.id);}
//					System.err.println(sb);
				}
			}else{
				if(discardedReads>0){
					falsePositiveZMWsT++;
				}else if(ambigReads>0){
					ambiguousZMWs++;
				}else{
					trueNegativeZMWsT++;
				}
			}

			if(rosg!=null && good!=null && !good.isEmpty()){rosg.add(good, 0);}
			if(rosa!=null && ambig!=null && !ambig.isEmpty()){rosa.add(ambig, 0);}
			if(rosb!=null && bad!=null && !bad.isEmpty()){rosb.add(bad, 0);}
			if(rosj!=null && junctions!=null && !junctions.isEmpty()){rosj.add(junctions, 0);}
		}
		
		/**
		 * Align part of a read to itself to look for inverted repeats.
		 */
		AlignmentResult align(final Read r, boolean fullPass, int passes, int readNum){
			if(!alignRC){return null;}
			final byte[] bases=r.bases;
			
			int qlen=(int)Tools.max(minQlen, Tools.min(targetQlen, bases.length*maxQlenFraction));
			if(qlen>0.45f*bases.length){return null;}//Ignore short stuff
			
			//Perform an initial scan using the tips of the reads to look for an inverted repeat
			boolean tipOnly=filterIceCreamOnly && fullPass;
//			System.err.println(filterIceCreamOnly+", "+fullPass+", "+tipOnly);
			AlignmentResult a=alignLeft(bases, qlen, minRatio1, true, tipOnly);
			AlignmentResult b=(fullPass && false ? null : alignRight(bases, qlen, minRatio1, true, tipOnly));//A second alignment is not needed for a full pass.
			AlignmentResult ar=(a==null ? b : b==null ? a : a.maxScore>=b.maxScore ? a : b);
			
			//If nothing was detected, return
			if(ar==null){return null;}
			ar.alignedRead=r;
			
			//At this point, the read contains an inverted repeat of length qlen.
			final int expectedJunction=ar.rLen/2;
			
			if(ar.left){
				ar.junctionLoc=ar.maxRpos/2;
			}else{
				int innerLeft=ar.maxRpos;
				int innerRight=ar.rLen-ar.qLen;
				ar.junctionLoc=(innerLeft+innerRight)/2;
			}
			
//			if(fullPass){//This code doesn't seem to have any effect and it's not clear why it is present
//				int dif=Tools.absdif(expectedJunction, ar.junctionLoc);
//				if(dif>expectedJunction*0.1) {
//					if(filterIceCreamOnly){return ar;}
//				}
//			}
			
			if(realign){
				if(ar.junctionLoc<expectedJunction){
					int qlen2=(int)(ar.junctionLoc*0.9);
					if(qlen2>=qlen){
						ar=alignLeft(bases, qlen2, minRatio2, false, false);
					}
				}else{
					int qlen2=(int)((ar.rLen-ar.junctionLoc)*0.9);
					if(qlen2>=qlen){
						ar=alignRight(bases, qlen2, minRatio2, false, false);
					}
				}
				if(ar==null){return null;}
				ar.alignedRead=r;
			}
			
			//At this point, the read contains an inverted repeat mirrored across a junction.
			final float junctionFraction;
			if(ar.left){
				ar.junctionLoc=ar.maxRpos/2;
				junctionFraction=ar.junctionLoc/(float)ar.rLen;
			}else{
				int innerLeft=ar.maxRpos;
				int innerRight=ar.rLen-ar.qLen;
				ar.junctionLoc=(innerLeft+innerRight)/2;
				junctionFraction=(ar.rLen-ar.junctionLoc)/(float)ar.rLen;
			}
			
			final int dif=Tools.absdif(expectedJunction, ar.junctionLoc);
			if(fullPass){
				if(junctionFraction<minJunctionFraction){
					ar.icecream=false;
				}else{
					ar.icecream=true;
				}
			}else if(passes==2){
				if(junctionFraction<minJunctionFraction){
					if(readNum==0){
						//First read, the junction should be closer to the beginning for real ice cream
						if(ar.junctionLoc>expectedJunction){
							ar.icecream=false;
						}else{
							ar.ambiguous=true;
						}
					}else{
						//Last read, the junction should be closer to the end for real ice cream
						assert(readNum==1) : readNum;
						if(ar.junctionLoc<expectedJunction){
							ar.icecream=false;
						}else{
							ar.ambiguous=true;
						}
					}
					return ar;
				}else{
					ar.icecream=true;
				}
			}else{
				if(junctionFraction<minJunctionFraction){
					ar.icecream=true;
				}else{
					ar.ambiguous=true;
				}
			}
			return ar;
		}
		
		/** Align the left qlen bases to the rest of the read. */
		private AlignmentResult alignLeft(final byte[] bases, final int qlen, final float minRatio, boolean v2, boolean tipOnly){
			final byte[] query=Arrays.copyOfRange(bases, 0, qlen);
			AminoAcid.reverseComplementBasesInPlace(query);
			final AlignmentResult ar;
			final int rstop=bases.length-1;
			final int rstart=(tipOnly ? Tools.max(qlen, rstop-(int)(tipRatio*qlen)) : qlen);
			if(v2){
//				elapsedShortT-=System.nanoTime();
				ar=ica.alignForwardShort(query, bases, rstart, rstop, minScore, minRatio);
//				elapsedShortT+=System.nanoTime();
			}else{
//				elapsedT-=System.nanoTime();
				ar=ica.alignForward(query, bases, rstart, rstop, minScore, minRatio);
//				elapsedT+=System.nanoTime();
			}
			if(ar!=null){ar.left=true;}
			return ar;
		}
		
		/** Align the right qlen bases to the rest of the read. */
		private AlignmentResult alignRight(final byte[] bases, final int qlen, final float minRatio, boolean v2, boolean tipOnly){
			final byte[] query=Arrays.copyOfRange(bases, bases.length-qlen, bases.length);
			AminoAcid.reverseComplementBasesInPlace(query);
			final AlignmentResult ar;
			final int rstop=(tipOnly ? Tools.min((int)(tipRatio*qlen), bases.length-qlen-1) : bases.length-qlen-1);
			final int rstart=0;
			if(v2){
//				elapsedShortT-=System.nanoTime();
				ar=ica.alignForwardShort(query, bases, rstart, rstop, minScore, minRatio);
//				elapsedShortT+=System.nanoTime();
			}else{
//				elapsedT-=System.nanoTime();
				ar=ica.alignForward(query, bases, rstart, rstop, minScore, minRatio);
//				elapsedT+=System.nanoTime();
			}
			if(ar!=null){ar.left=false;}
			return ar;
		}

		/** Align the adapter sequence to the read, piecewise, to count matches. */
		private int lookForAdapter(Read r, double mult) {
			assert(!r.hasAdapter() && !r.tested());
			int begin=0;
			while(begin<r.length() && r.bases[begin]=='N'){begin++;} //Skip reads made of 'N'
			if(begin>=r.length()){return 0;}

			int suspectThresh=(int)(minSwScoreSuspect*mult);
			int scoreThresh=(int)(minSwScore*mult);
			
//			final byte[] array=npad(r.bases, pad ? npad : 0);
			final byte[] array=npad(r.bases, npad);
			
//			assert(array==r.bases) : npad;
			
			int lim=array.length-npad-stride;
			
			int found=0;

			int lastSuspect=-1;
			int lastConfirmed=-1;
			int maxScore=0;
			
//			final int revisedStride=(r.length()>1800 ? stride : (int)(stride*0.6f));
			for(int i=begin; i<lim; i+=stride){
				int j=Tools.min(i+window, array.length-1);
				if(j-i<window){
					lim=0; //Last loop cycle
//					i=Tools.max(0, array.length-2*query1.length);
				}
				
				{
					
					int[] rvec;
//					if(timeless){//A speed-optimized version
					rvec=ssa.fillAndScoreLimited(adapter, array, i, j, suspectThresh);
//					}else{rvec=msa.fillAndScoreLimited(adapter, array, i, j, suspectThresh);}
					if(rvec!=null && rvec[0]>=suspectThresh){
						final int score=rvec[0];
						final int start=rvec[1];
						final int stop=rvec[2];
						assert(score>=suspectThresh);
						if((i==0 || start>i) && (j==array.length-1 || stop<j)){
							boolean kill=(score>=scoreThresh ||
									(score>=suspectMidpoint && lastSuspect>0 && start>=lastSuspect && start-lastSuspect<suspectDistance) ||
									(lastConfirmed>0 && start>=lastConfirmed && start-lastConfirmed<suspectDistance));
							
							if(!kill && useLocality && array.length-stop>window){//Look ahead
								rvec=ssa.fillAndScoreLimited(adapter, array, stop, stop+window, suspectThresh);
								if(rvec!=null){
									if(score>=suspectMidpoint && rvec[0]>=suspectThresh && rvec[1]-stop<suspectDistance){kill=true;}
									else if(score>=suspectThresh && rvec[0]>=scoreThresh && rvec[1]-stop<suspectDistance){kill=true;}
								}
							}
							
							if(!kill && useAltMsa){//Try alternate msa
								rvec=msa2.fillAndScoreLimited(adapter, array, Tools.max(0, start-4), Tools.min(stop+4, array.length-1), suspectThresh);
								if(rvec!=null && rvec[0]>=(scoreThresh)){kill=true;}
							}
							
							if(kill){
								maxScore=Tools.max(maxScore, score);
//								if(print) {System.err.println("Found adapter at distance "+Tools.min(start, array.length-stop));}
//								System.out.println("-:"+score+", "+scoreThresh+", "+suspectThresh+"\t"+lastSuspect+", "+start+", "+stop);
								found++;
								for(int x=Tools.max(0, start); x<=stop; x++){array[x]='X';}
//								assert(Shared.threads()!=1) : new String(array)+"\n"+start+", "+stop;
								if(useLocality && score>=scoreThresh){lastConfirmed=Tools.max(lastConfirmed, stop);}
							}
						}
//						System.out.println("Set lastSuspect="+stop+" on score "+score);
						if(useLocality){lastSuspect=Tools.max(lastSuspect, stop);}
					}
				}
			}

			r.setTested(true);
			if(found>0){
				r.setHasAdapter(true);
				r.setDiscarded(true);
				r.setAmbiguous(false);
				
				int idx=(int)((maxScore*200.0)/maxSwScore);
				adapterScoresT[idx]++;
			}else{
				r.setHasAdapter(false);
			}
			
			return found;
		}

		/** Align the adapter sequence to the read ends, and trim if needed. */
		private int trimTerminalAdapters(Read r, double mult) {
			assert(!r.hasAdapter() && !r.tested());
			int begin=0;
			while(begin<r.length() && r.bases[begin]=='N'){begin++;} //Skip reads made of 'N'
			if(begin>=r.length()){return 0;}

			int suspectThresh=(int)(minSwScoreSuspect*mult);
			int scoreThresh=(int)(minSwScore*mult);
			
//			final byte[] array=npad(r.bases, pad ? npad : 0);
			final byte[] array=npad(r.bases, npad);
			
//			assert(array==r.bases) : npad;
			
			int lim=array.length-npad-stride;
			
			int found=0;

			int lastSuspect=-1;
			int lastConfirmed=-1;
			int maxScore=0;
			
//			final int revisedStride=(r.length()>1800 ? stride : (int)(stride*0.6f));
			for(int i=begin; i<lim; i+=stride){
				int j=Tools.min(i+window, array.length-1);
				if(j-i<window){
					lim=0; //Last loop cycle
//					i=Tools.max(0, array.length-2*query1.length);
				}
				
				{
					
					int[] rvec;
//					if(timeless){//A speed-optimized version
					rvec=ssa.fillAndScoreLimited(adapter, array, i, j, suspectThresh);
//					}else{rvec=msa.fillAndScoreLimited(adapter, array, i, j, suspectThresh);}
					if(rvec!=null && rvec[0]>=suspectThresh){
						final int score=rvec[0];
						final int start=rvec[1];
						final int stop=rvec[2];
						assert(score>=suspectThresh);
						if((i==0 || start>i) && (j==array.length-1 || stop<j)){
							boolean kill=(score>=scoreThresh ||
									(score>=suspectMidpoint && lastSuspect>0 && start>=lastSuspect && start-lastSuspect<suspectDistance) ||
									(lastConfirmed>0 && start>=lastConfirmed && start-lastConfirmed<suspectDistance));
							
							if(!kill && useLocality && array.length-stop>window){//Look ahead
								rvec=ssa.fillAndScoreLimited(adapter, array, stop, stop+window, suspectThresh);
								if(rvec!=null){
									if(score>=suspectMidpoint && rvec[0]>=suspectThresh && rvec[1]-stop<suspectDistance){kill=true;}
									else if(score>=suspectThresh && rvec[0]>=scoreThresh && rvec[1]-stop<suspectDistance){kill=true;}
								}
							}
							
							if(!kill && useAltMsa){//Try alternate msa
								rvec=msa2.fillAndScoreLimited(adapter, array, Tools.max(0, start-4), Tools.min(stop+4, array.length-1), suspectThresh);
								if(rvec!=null && rvec[0]>=(scoreThresh)){kill=true;}
							}
							
							if(kill){
								maxScore=Tools.max(maxScore, score);
//								if(print) {System.err.println("Found adapter at distance "+Tools.min(start, array.length-stop));}
//								System.out.println("-:"+score+", "+scoreThresh+", "+suspectThresh+"\t"+lastSuspect+", "+start+", "+stop);
								found++;
								for(int x=Tools.max(0, start); x<=stop; x++){array[x]='X';}
//								assert(Shared.threads()!=1) : new String(array)+"\n"+start+", "+stop;
								if(useLocality && score>=scoreThresh){lastConfirmed=Tools.max(lastConfirmed, stop);}
							}
						}
//						System.out.println("Set lastSuspect="+stop+" on score "+score);
						if(useLocality){lastSuspect=Tools.max(lastSuspect, stop);}
					}
				}
			}

			r.setTested(true);
			if(found>0){
				r.setHasAdapter(true);
				r.setDiscarded(true);
				r.setAmbiguous(false);
				
				int idx=(int)((maxScore*200.0)/maxSwScore);
				adapterScoresT[idx]++;
			}else{
				r.setHasAdapter(false);
			}
			
			return found;
		}
		
		private byte[] npad(final byte[] array, final int pad){
			if(pad<=0){return array;}
			final int len=array.length+2*pad;
			byte[] r=new byte[len];
			for(int i=0; i<pad; i++){r[i]='N';}
			for(int i=pad, j=0; j<array.length; i++, j++){r[i]=array[j];}
			for(int i=array.length+pad, limit=Tools.min(r.length, len+50); i<limit; i++){r[i]='N';}
			return r;
		}
		
		/*--------------------------------------------------------------*/
		/*----------------      Inner Class Fields      ----------------*/
		/*--------------------------------------------------------------*/

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		

		protected long truePositiveReadsT=0;
		protected long falsePositiveReadsT=0;
		protected long trueNegativeReadsT=0;
		protected long falseNegativeReadsT=0;
		protected long ambiguousReadsT=0;
		
		protected long truePositiveZMWsT=0;
		protected long falsePositiveZMWsT=0;
		protected long trueNegativeZMWsT=0;
		protected long falseNegativeZMWsT=0;
		protected long ambiguousZMWsT=0;
		
		
		protected long elapsedT=0;
		protected long elapsedShortT=0;
		
		//Unused
		protected long iceCreamReadsT=0;
		protected long iceCreamBasesT=0;
		
		protected long iceCreamZMWsT=0;
		protected double iceCreamRatioT=0;
		protected long ratiosCountedT=0;
		protected long missingAdapterZMWsT=0;
		protected long hiddenAdapterZMWsT=0;
		
		protected long basesTrimmedT=0;
		protected long readsTrimmedT=0;
		
		/** Number of reads retained by this thread */
		protected long readsOutT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;
		/** Number of junctions detected in IR reads by this thread */
		protected long junctionsOutT=0;
		
		protected long lowEntropyZMWsT=0;
		protected long lowEntropyReadsT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Shared output stream */
		private final ConcurrentReadOutputStream rosg;
		/** Shared output stream for ambiguous reads */
		private final ConcurrentReadOutputStream rosa;
		/** Shared output stream for bad reads */
		private final ConcurrentReadOutputStream rosb;
		/** Shared output stream for junctions */
		private final ConcurrentReadOutputStream rosj;
		/** Thread ID */
		final int tid;
		
		/* Aligners for this thread */
		
		/** For inverted repeat alignment */
		final IceCreamAligner ica=IceCreamAligner.makeAligner(32);
		/** For quickly aligning adapter sequence to whole read */
		final FlatAligner fla=new FlatAligner();
//		final MultiStateAligner9PacBioAdapter msa=alignAdapter ? new MultiStateAligner9PacBioAdapter(ALIGN_ROWS, ALIGN_COLUMNS) : null;
		/** For slowly aligning adapter sequence sectionwise */
		final SingleStateAlignerPacBioAdapter ssa=(alignAdapter || trimReads) ? new SingleStateAlignerPacBioAdapter(ALIGN_ROWS, ALIGN_COLUMNS, adapter.length) : null;
		/** Alternate scoring for slow adapter alignment */
		final MultiStateAligner9PacBioAdapter2 msa2=(alignAdapter || trimReads) ? new MultiStateAligner9PacBioAdapter2() : null;

		final ZMWStreamer zstream;
		final EntropyTracker eTracker;
		
		final long[] adapterScoresT=new long[201];
		final long[] repeatScoresT=new long[201];
		final byte[] tipBufferLeft=new byte[adapterTipLen+adapterTipPad];
		final byte[] tipBufferRight=new byte[adapterTipLen+adapterTipPad];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1=null;

	/** Primary output file path */
	private String outg=null;
	/** Ambiguous output file path */
	private String outa=null;
	/** Bad output file path */
	private String outb=null;
	/** Junction output file path */
	private String outj=null;
	/** Stats (screen) output file path */
	private String outstats=null;

	/** Adapter score ratio histogram */
	private String asrhist=null;

	/** Inverted repeat score ratio histogram */
	private String irsrhist=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;

	private int targetQlen=352;
	private int minQlen=100;
	
	/** Make a query at most this fraction of read length */
	private float maxQlenFraction=0.15f;
	
	/** Exit alignment early if score drops below this.
	 * An aggressive optimization that may miss some low-quality inverted repeats.
	 * -700 seems safe. */
	private int minScore=-800;
	
	/** Fraction of maximum alignment score to consider as matching for initial alignment */ 
	private float minRatio1=0.59f;
	
	/** Fraction of maximum alignment score to consider as matching for realignment */ 
	private float minRatio2=0.64f;
	
	private float adapterRatio=0.18f; //.18 for fa, or 0.57/0.6 for ica
	private float adapterRatio2=0.325f; //0.31f normal, 0.325 timeless
	private float suspectRatio=0.85f;
	private boolean useLocality=true;
	private boolean useAltMsa=true;
	
	
	private float tipRatio=1.5f;
	
	private float longReadMult=1.5f;
	
	private float shortReadMult=1.5f;
	
	private float veryShortMult=0.35f;
	
	/** Short half of inverted repeat must be at least this fraction of read length to be classed as a triangle */
	private float minJunctionFraction=0.4f;
	
	/** Only filter triangle reads, not all inverted repeats */
	private boolean filterIceCreamOnly=true;
	
	/** Align again once the junction position is provisionally identified */
	private boolean realign=true;
	
	/** For internal read array transfers */ 
	private int queuelen=80;
	
	/** For grading synthetic data */
	private boolean parseCustom;
	
	/** Input reads are CCS (full pass) */
	private boolean CCS;
	
	private boolean modifyHeader=false;

	private boolean sendAmbigToGood=true;
	private boolean sendAmbigToBad=false;
	private boolean setAmbig=false;
	
	//Note: These flags are very similar... they need to be better-defined or merged.
	private boolean keepZMWsTogether=false;
	private boolean keepShortReads=true;
	
	private int format=FORMAT_TEXT;
	private static final int FORMAT_TEXT=1, FORMAT_JSON=2;
	
	/*--------------------------------------------------------------*/

	/** Alignment iterations for inverted repeat calculation with ref columns and query rows */
	protected long alignmentIters=0;

	/** Alignment iterations for inverted repeat calculation with query columns and ref rows */
	protected long alignmentItersShort=0;
	
	/** Time spent in long iterations */
	protected long elapsed=0;
	
	/** Time spent in short iterations */
	protected long elapsedShort=0;
	
	/** Print iteration time statistics */
	protected boolean printTiming=false;
	
	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Number of reads retained */
	protected long readsOut=0;
	/** Number of bases retained */
	protected long basesOut=0;
	/** Number of junctions detected in IR reads */
	protected long junctionsOut=0;
	
	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;

	//Unused
	protected long iceCreamReads=0;
	protected long iceCreamBases=0;
	
	/** ZMWs with discarded reads */
	protected long iceCreamZMWs=0;
	
	/** Sum of IR alignment ratios for IR reads not containing adapters */
	protected double iceCreamRatio=0;
	
	/** Number of ratios in iceCreamRatio */
	protected long ratiosCounted=0;
	
	/** Histogram */
	protected final long[] adapterScores=new long[201];
	
	/** Histogram */
	protected final long[] repeatScores=new long[201];
	
	/** ZMWs with a triangle read but no hidden adapters */
	protected long missingAdapterZMWs=0;
	
	/** ZMWs with hidden adapters */
	protected long hiddenAdapterZMWs=0;
	
	protected long basesTrimmed=0;
	protected long readsTrimmed=0;
	
	protected long lowEntropyZMWs=0;
	protected long lowEntropyReads=0;
	
	/** Total ZMWs observed */
	protected long ZMWs=0;
	
	protected long truePositiveReads=0;
	protected long falsePositiveReads=0;
	protected long trueNegativeReads=0;
	protected long falseNegativeReads=0;
	protected long ambiguousReads=0;
	
	protected long truePositiveZMWs=0;
	protected long falsePositiveZMWs=0;
	protected long trueNegativeZMWs=0;
	protected long falseNegativeZMWs=0;
	protected long ambiguousZMWs=0;
	
	protected final int stride;
	protected final int window;
	protected final int ALIGN_ROWS;
	protected final int ALIGN_COLUMNS;

	private boolean timeless=true;
	protected final int maxSwScore;
	protected final int minSwScore;
	protected final int minSwScoreSuspect;
	protected final int maxImperfectSwScore;
	protected final int suspectMidpoint;
	protected final int suspectDistance=100;
	protected int npad=0;  //This is for catching adapters on the tips, which are not very relevant to ice cream.  Set to 35 for tip apdapters.
	
	private byte[] adapter="ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT".getBytes(); //This one seems to be correct
//	private byte[] adapter="ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT".getBytes();
	private boolean alignAdapter=true;
	private boolean alignRC=true;
	private boolean flagLongReads=true;
	private boolean trimReads=true;
	private int minLengthAfterTrimming=40;
	private int adapterTipLen=100;
	private int adapterTipPad=35;
	
	boolean trimPolyA=false;
	
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
	
	/** Ambiguous output file */
	private final FileFormat ffouta;
	
	/** Bad output file */
	private final FileFormat ffoutb;
	
	/** Junction output file */
	private final FileFormat ffoutj;

	/** Stats output file */
	private final FileFormat ffstats;

	/** Adapter score ratio histogram */
	private final FileFormat ffasrhist;

	/** Inverted repeat score ratio histogram */
	private final FileFormat ffirsrhist;
	
	private final int threads;
	
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
