package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

import cardinality.CardinalityTracker;
import dna.AminoAcid;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import hiseq.IlluminaHeaderParser2;
import hiseq.ReadHeaderParser;
import shared.KillSwitch;
import shared.MetadataWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import structures.ListNum;
import structures.LongList;
import structures.Quantizer;
import structures.SuperLongList;
import tracker.ReadStats;
import var2.ScafMap;

/**
 * @author Brian Bushnell
 * @date Sep 11, 2012
 *
 */
public class ReformatReads {

	public static void main(String[] args){
		Timer t=new Timer();
		ReformatReads x=new ReformatReads(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public ReformatReads(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean setInterleaved=false; //Whether it was explicitly set.

		
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		
		SamLine.SET_FROM_OK=true;
//		SamLine.CONVERT_CIGAR_TO_MATCH=true;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=Parse.splitOnFirst(arg, '=');
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("passes")){
				assert(false) : "'passes' is disabled.";
//				passes=Integer.parseInt(b);
			}else if(a.equals("path")){
				Data.setPath(b);
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("sample") || a.equals("samplereads") || a.equals("samplereadstarget") || a.equals("srt")){
				sampleReadsTarget=Parse.parseKMG(b);
				sampleReadsExact=(sampleReadsTarget>0);
			}else if(a.equals("samplebases") || a.equals("samplebasestarget") || a.equals("sbt")){
				sampleBasesTarget=Parse.parseKMG(b);
				sampleBasesExact=(sampleBasesTarget>0);
			}else if(a.equals("prioritizelength")){
				prioritizeLength=Parse.parseBoolean(b);
			}else if(a.equals("upsample")){
				allowUpsample=Parse.parseBoolean(b);
			}else if(a.equals("addslash")){
				addslash=Parse.parseBoolean(b);
			}else if(a.equals("addcolon")){
				addcolon=Parse.parseBoolean(b);
			}else if(a.equals("slashspace") || a.equals("spaceslash")){
				boolean x=Parse.parseBoolean(b);
				if(x){
					slash1=" /1";
					slash2=" /2";
				}else{
					slash1="/1";
					slash2="/2";
				}
			}else if(a.equals("addunderscore") || a.equals("underscore")){
				addunderscore=Parse.parseBoolean(b);
			}else if(a.equals("uniquenames")){
				uniqueNames=Parse.parseBoolean(b);
			}else if(a.equals("verifyinterleaved") || a.equals("verifyinterleaving") || a.equals("vint")){
				verifyinterleaving=Parse.parseBoolean(b);
			}else if(a.equals("verifypaired") || a.equals("verifypairing") || a.equals("vpair")){
				verifypairing=Parse.parseBoolean(b);
			}else if(a.equals("allowidenticalnames") || a.equals("ain")){
				allowIdenticalPairNames=Parse.parseBoolean(b);
			}else if(a.equals("rcompmate") || a.equals("rcm")){
				reverseComplementMate=Parse.parseBoolean(b);
				outstream.println("Set RCOMPMATE to "+reverseComplementMate);
			}else if(a.equals("rcomp") || a.equals("rc") || a.equals("reversecomplement")){
				reverseComplement=Parse.parseBoolean(b);
				outstream.println("Set RCOMP to "+reverseComplement);
			}else if(a.equals("comp") || a.equals("complement")){
				complement=Parse.parseBoolean(b);
				outstream.println("Warning: Setting non-standard parameter 'complement'.");
			}
			
			//Unneeded due to pad flag
//			else if(a.equals("npad") || a.equals("padns")){
//				assert(b!=null) : "npad requires an int value, e.g. npad=10";
//				assert(Character.isDigit(b.charAt(0)) || b.charAt(0)=='f') : "npad requires an int value, e.g. npad=10";
//				if(Character.isDigit(b.charAt(0)){npad=Integer.parseInt(b);}
//				else {npad=0;}
//			}
			
			else if(a.equals("deleteempty") || a.equals("deletempty") || a.equals("delempty") || a.equals("def")){
				deleteEmptyFiles=Parse.parseBoolean(b);
			}else if(a.equals("mappedonly")){
				mappedOnly=Parse.parseBoolean(b);
			}else if(a.equals("pairedonly")){
				pairedOnly=Parse.parseBoolean(b);
			}else if(a.equals("unpairedonly")){
				unpairedOnly=Parse.parseBoolean(b);
			}else if(a.equals("unmappedonly")){
				unmappedOnly=Parse.parseBoolean(b);
			}else if(a.equals("requiredbits") || a.equals("rbits")){
				requiredBits=Parse.parseIntHexDecOctBin(b);
			}else if(a.equals("filterbits") || a.equals("fbits")){
				filterBits=Parse.parseIntHexDecOctBin(b);
			}else if(a.equals("primaryonly")){
				primaryOnly=Parse.parseBoolean(b);
			}else if(a.equals("remap1")){
				remap1=Parse.parseRemap(b);
			}else if(a.equals("remap2")){
				remap2=Parse.parseRemap(b);
			}else if(a.equals("remap")){
				remap1=remap2=Parse.parseRemap(b);
			}else if(a.equals("minmapq")){
				minMapq=Integer.parseInt(b);
			}else if(a.equals("maxmapq")){
				maxMapq=Integer.parseInt(b);
			}else if(a.equals("skipreads")){
				skipreads=Parse.parseKMG(b);
			}else if(a.equals("undefinedton") || a.equals("iupacton") || a.equals("itn")){
				iupacToN=Parse.parseBoolean(b);
			}else if(a.equals("bottom")){
				bottom=Parse.parseBoolean(b);
			}else if(a.equals("top")){
				top=Parse.parseBoolean(b);
			}
			
			else if(a.equals("quantize") || a.equals("quantizesticky")){
				quantizeQuality=Quantizer.parse(arg, a, b);
			}else if(a.equals("samversion") || a.equals("samv") || a.equals("sam")){
				parser.parse(arg, a, b);
				fixCigar=true;
			}else if(a.equals("deleteinput")){
				deleteInput=Parse.parseBoolean(b);
			}else if(a.equals("invert") || a.equals("invertfilters")){
				invertFilters=Parse.parseBoolean(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(a.equals("top")){
				top=Parse.parseBoolean(b);
			}else if(a.equals("bottom")){
				bottom=Parse.parseBoolean(b);
			}
			
			else if(a.equalsIgnoreCase("lastUnderscoreToSpace") || a.equalsIgnoreCase("underscoreToSpace")){
				lastUnderscoreToSpace=Parse.parseBoolean(b);
			}
			
			else if(a.equals("fixheader") || a.equals("fixheaders") || a.equals("fixnames")){
				fixHeaders=Parse.parseBoolean(b);
			}
			
			else if(a.equals("padleft")){
				padLeft=Integer.parseInt(b);
			}else if(a.equals("padright")){
				padRight=Integer.parseInt(b);
			}else if(a.equals("pad")){
				assert(b!=null) : "No value for pad.";
				if(Character.isLetter(b.charAt(0))){padSymbol=(byte)b.charAt(0);}
				else{padLeft=padRight=Integer.parseInt(b);}
			}else if(a.equals("padsymbol")){
				padSymbol=(byte)b.charAt(0);
			}else if(a.equals("padq")){
				padQ=(byte)Integer.parseInt(b);
			}
			
			else if(a.equals("minsubs")){
				minsubs=Integer.parseInt(b);
			}
			
			else if(a.equals("tag")){
				tag=b;
			}else if(a.equals("delimiter")){
//				delimiterS=null;
				delimiter=Parse.parseSymbolToCharacter(b);
//				if(delimiter==0) {delimiterS=Parse.parseSymbol(b);}
				assert(delimiter>0) : delimiter+", '"+b+"'";
			}else if(a.equals("minvalue")){
				minValue=Integer.parseInt(b);
			}else if(a.equals("maxvalue")){
				maxValue=Integer.parseInt(b);
			}else if(a.equals("value")){
				requiredValue=b;
			}
			
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		if(k>0){parser.loglogk=k;}
		if(parser.loglog && k<1){k=parser.loglogk;}
		pad=padLeft>0 || padRight>0;
		if(AminoAcid.isFullyDefined(padSymbol)){padQ=Tools.max(padQ, (byte)2);}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			samplerate=parser.samplerate;
			sampleseed=parser.sampleseed;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			testsize=parser.testsize;
			trimBadSequence=parser.trimBadSequence;
			breakLength=parser.breakLength;
			stoptag=SamLine.MAKE_STOP_TAG;
			
			forceTrimModulo=parser.forceTrimModulo;
			forceTrimLeft=parser.forceTrimLeft;
			forceTrimRight=parser.forceTrimRight;
			forceTrimRight2=parser.forceTrimRight2;
			qtrimLeft=parser.qtrimLeft;
			qtrimRight=parser.qtrimRight;
			trimq=parser.trimq;
			trimE=parser.trimE();
			minAvgQuality=parser.minAvgQuality;
			minAvgQualityBases=parser.minAvgQualityBases;
			chastityFilter=parser.chastityFilter;
			failBadBarcodes=parser.failBadBarcodes;
			removeBadBarcodes=parser.removeBadBarcodes;
			failIfNoBarcode=parser.failIfNoBarcode;
			barcodes=parser.barcodes;
			maxNs=parser.maxNs;
			minConsecutiveBases=parser.minConsecutiveBases;
			minReadLength=parser.minReadLength;
			maxReadLength=parser.maxReadLength;
			minLenFraction=parser.minLenFraction;
			requireBothBad=parser.requireBothBad;
			minGC=parser.minGC;
			maxGC=parser.maxGC;
			filterGC=(minGC>0 || maxGC<1);
			usePairGC=parser.usePairGC;
			tossJunk=parser.tossJunk;
			recalibrateQuality=parser.recalibrateQuality;

			minIdFilter=parser.minIdFilter;
			maxIdFilter=parser.maxIdFilter;
			subfilter=parser.subfilter;
			clipfilter=parser.clipfilter;
			delfilter=parser.delfilter;
			insfilter=parser.insfilter;
			indelfilter=parser.indelfilter;
			dellenfilter=parser.dellenfilter;
			inslenfilter=parser.inslenfilter;
			editfilter=parser.editfilter;
			
			USE_EDIT_FILTER=(subfilter>-1 || minsubs>-1 || delfilter>-1 || insfilter>-1 || indelfilter>-1 || dellenfilter>-1 || inslenfilter>-1 || editfilter>-1 || clipfilter>-1);

			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;
			
			out1=parser.out1;
			out2=parser.out2;
			outsingle=parser.outsingle;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
			
			loglog=(parser.loglog ? CardinalityTracker.makeTracker(parser) : null);
			silent=Parser.silent;
		}
		
		if(recalibrateQuality){CalcTrueQuality.initializeMatrices();}
		
		if(SamLine.setxs && !SamLine.setintron){SamLine.INTRON_LIMIT=10;}
		qtrim=qtrimLeft||qtrimRight;
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		
		if(out1!=null && in1!=null && out1.indexOf('%')>-1){
			out1=out1.replace("%", ReadWrite.stripExtension(in1));
		}
		if(out2!=null && out2.indexOf('%')>-1){
			if(in2!=null){
				out2=out2.replace("%", ReadWrite.stripExtension(in2));
			}else if(in1!=null){
				out2=out2.replace("%", ReadWrite.stripExtension(in1));
			}
		}
		
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		if(verifyinterleaving || (verifypairing && in2==null)){
			verifypairing=true;
			setInterleaved=true;
//			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=true;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
		qfin1=Tools.fixExtension(qfin1);
		qfin2=Tools.fixExtension(qfin2);
		ref=Tools.fixExtension(ref);
		if(!Tools.testInputFiles(false, true, in1, in2, qfin1, qfin2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1==null){
			if(out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
			
			if(!parser.setOut){
				outstream.println("No output stream specified.  To write to stdout, please specify 'out=stdout.fq' or similar.");
//				out1="stdout";
			}
		}
		
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
		
		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
		if(outsingle!=null && outsingle.equalsIgnoreCase("null")){outsingle=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, qfout1, qfout2, outsingle)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		if(!Tools.testForDuplicateFiles(true, in1, in2, qfin1, qfin2, out1, out2, qfout1, qfout2, outsingle) || !ReadStats.testFiles(false)){
			throw new RuntimeException("Duplicate filenames are not allowed.");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffoutsingle=FileFormat.testOutput(outsingle, FileFormat.FASTQ, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		parser.validateStdio(ffin1, ffout1, ffoutsingle);

		assert(ReadStats.testFiles(true)) : "Existing output files specified, but overwrite==false";
		assert(ReadStats.testFiles(false)) : "Duplicate or output files specified";

//		outstream.println("\n"+ReadWrite.USE_PIGZ+", "+ReadWrite.USE_UNPIGZ+", "+Data.PIGZ()+", "+Data.UNPIGZ()+", "+ffin1+"\n");
//		assert(false) : ReadWrite.USE_PIGZ+", "+ReadWrite.USE_UNPIGZ+", "+Data.PIGZ()+", "+Data.UNPIGZ()+", "+ffin1;

		nameMap1=(uniqueNames ? new HashMap<String, Integer>() : null);
		nameMap2=(uniqueNames ? new HashMap<String, Integer>() : null);
		
		if(ref!=null){ScafMap.loadReference(ref, true);}
		if(tag!=null) {
			assert(delimiter>0/* || delimiterS!=null*/) :
				"When using a tag, a delimiter must be set; e.g. "
				+ "delimiter=X, delimiter=' ' or delimiter=space\n"
				+ "Most problematic symbols can be spelled out, like tab, "
				+ "pipe, asterisk, greaterthan, etc.";
//			if(delimiterS!=null) {
//				splitter=Pattern.compile(delimiterS);
//			}
//			lp=new LineParser2S(delimiter);
		}
//		lp=null;//lp is not currently used
	}

	
	
	void process(Timer t){
		
		long readsRemaining=0;
		long basesRemaining=0;
		
		if(sampleReadsExact || sampleBasesExact){
			if(prioritizeLength){
				SuperLongList sll=makeLengthHist(maxReads);
				LongList list=sll.list();
				long[] array=sll.array();
				if(sampleReadsExact){
					long sum=0;
					for(int i=list.size()-1; i>=0 && sum<sampleReadsTarget; i--){
						long num=list.get(i);
						sum++;
						if(sum>=sampleReadsTarget){
							minReadLength=Tools.max(minReadLength, (int)num);
						}
					}
					for(int i=array.length-1; i>=0 && sum<sampleReadsTarget; i--){
						long count=array[i];
						sum+=count;
						if(sum>=sampleReadsTarget){
							minReadLength=Tools.max(minReadLength, i);
						}
					}
				}else{
					long sum=0;
					for(int i=list.size()-1; i>=0 && sum<sampleBasesTarget; i--){
						long num=list.get(i);
						sum+=num;
						if(sum>=sampleBasesTarget){
							minReadLength=Tools.max(minReadLength, (int)num);
						}
					}
					for(int i=array.length-1; i>=0 && sum<sampleBasesTarget; i--){
						long count=array[i];
						sum+=(count*i);
						if(sum>=sampleBasesTarget){
							minReadLength=Tools.max(minReadLength, i);
						}
					}
//					System.err.println(sum+", "+minReadLength+", "+sll.median()+", "+sll.mean());
				}
//				System.err.println(sll.toString());
				sampleReadsExact=sampleBasesExact=false;
			}else{
				long[] counts=countReads(maxReads);
				readsRemaining=counts[0];
				basesRemaining=counts[2];
				setSampleSeed(sampleseed);
				double prob=(sampleReadsExact ? sampleReadsTarget/(double)(readsRemaining)
						: sampleBasesTarget/(double)(basesRemaining));
				if(!allowUpsample) {prob=Tools.min(prob, 1.0);}
				System.err.println("Initial samplerate set to "+String.format("%.6f", prob));
			}
		}
		
		
		final ConcurrentReadInputStream cris;
		{
			useSharedHeader=(ffin1.samOrBam() && ffout1!=null && ffout1.samOrBam());
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, useSharedHeader, ffin1, ffin2, qfin1, qfin2);
			cris.setSampleRate(samplerate, sampleseed);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
		if(!silent){
			outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));
		}
		
		assert(!paired || breakLength<1) : "Paired input cannot be broken with 'breaklength'";

		final ConcurrentReadOutputStream ros;
		if(ffout1!=null){
			final int buff=4;
			
			if(cris.paired() && ffout2==null && ffout1!=null && !ffout1.samOrBam()){
				outstream.println("Writing interleaved.");
			}
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, useSharedHeader);
			ros.start();
		}else{ros=null;}
		
		final ConcurrentReadOutputStream rosb;
		if(ffoutsingle!=null){
			final int buff=4;
			
			rosb=ConcurrentReadOutputStream.getStream(ffoutsingle, null, buff, null, useSharedHeader);
			rosb.start();
		}else{rosb=null;}
		final boolean discardTogether=(!paired || (ffoutsingle==null && !requireBothBad));
		
		long basesFTrimmedT=0;
		long readsFTrimmedT=0;
		
		long basesQTrimmedT=0;
		long readsQTrimmedT=0;
		
		long lowqBasesT=0;
		long lowqReadsT=0;
		
		long badGcBasesT=0;
		long badGcReadsT=0;
		
		long readShortDiscardsT=0;
		long baseShortDiscardsT=0;
		
		long unmappedReadsT=0;
		long unmappedBasesT=0;
		
		long idfilteredReadsT=0;
		long idfilteredBasesT=0;
		
		long basesSwappedT=0;
		long readsSwappedT=0;

		final boolean MAKE_IHIST=ReadStats.COLLECT_INSERT_STATS;
		
		final boolean sam=ffin1.samOrBam();
		
		final ReadStats readstats=ReadStats.collectingStats() ? new ReadStats() : null;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			Read prevRead=null;
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				if(skipreads>0){
					int removed=0;
					for(int i=0; i<reads.size(); i++){
						Read r=reads.get(i);
						if(r.numericID<skipreads){
							reads.set(i, null);
							removed++;
						}else{
							skipreads=-1;
							break;
						}
					}
					if(removed>0){
						Tools.condenseStrict(reads);
					}
				}
				
				ArrayList<Read> singles=(rosb==null ? null : new ArrayList<Read>(32));
				
				if(breakLength>0){
					Tools.breakReads(reads, breakLength, minReadLength, verbose? outstream : null);
				}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());

					final int minlen1=(int)Tools.max(initialLength1*minLenFraction, minReadLength);
					final int minlen2=(int)Tools.max(initialLength2*minLenFraction, minReadLength);
					
					final SamLine sl1;
					if(sam){
						sl1=r1.samline;
						assert(r2==null);
//						sl2=(r2==null ? null : r2.samline);
//						if(sl1!=null && Shared.TRIM_RNAME){
//							sl1.setRname(Tools.trimToWhitespace(sl1.rname()));
//							sl1.setRnext(Tools.trimToWhitespace(sl1.rnext()));
//						}
					}else{
						sl1=null;//=sl2=null;
					}
					
					if(fixHeaders){
						fixHeader(r1);
						fixHeader(r2);
					}
					
					if(lastUnderscoreToSpace){
						underscoreToSpace(r1);
						underscoreToSpace(r2);
					}
					
					if(readstats!=null){
						readstats.addToHistograms(r1);
						
						if(MAKE_IHIST && sl1!=null && !r1.secondary() && sl1.pairnum()==0){
							readstats.addToInsertHistogram(sl1);
						}
						
//						if(MAKE_IHIST && !r1.secondary()){
//							if(sl1!=null && sl2!=null && sl1.qname.equals(sl2.qname) && sl1.pairnum()!=sl2.pairnum()){
//								readstats.addToInsertHistogram(sl1, sl2);
//							}
//						}
					}
					if(!r1.secondary()){prevRead=r1;}
					
					if(loglog!=null){loglog.hash(r1);}
					
					if(k>0){
						{
							final int kmers=Tools.countKmers(r1.bases, k);
							kmersProcessed+=kmers;
							correctKmers+=(r1.quality==null ? kmers : Tools.countCorrectKmers(r1.quality, k));
						}
						if(r2!=null){
							final int kmers=Tools.countKmers(r2.bases, k);
							kmersProcessed+=kmers;
							correctKmers+=(r2.quality==null ? kmers : Tools.countCorrectKmers(r2.quality, k));
						}
					}
					
					{
						readsProcessed++;
						basesProcessed+=initialLength1;
						if(reverseComplement){r1.reverseComplement();}
						if(complement){r1.complement();}
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=initialLength2;
						if(reverseComplement || reverseComplementMate){r2.reverseComplement();}
						if(complement){r2.complement();}
					}
					
					if(!bottom || !top) {
						rhp.parse(r1);
						if(!bottom && rhp.surface()==1) {
							r1.setDiscarded(true);
							r2.setDiscarded(true);
						}
						if(!top && rhp.surface()==2) {
							r1.setDiscarded(true);
							r2.setDiscarded(true);
						}
					}
					
					if(verifypairing){
						String s1=r1==null ? null : r1.id;
						String s2=r2==null ? null : r2.id;
						boolean b=FASTQ.testPairNames(s1, s2, allowIdenticalPairNames);
						if(!b){
							outstream.println("Names do not appear to be correctly paired.\n"+s1+"\n"+s2+"\n");
							ReadWrite.closeStreams(cris, ros);
							System.exit(1);
						}
					}
					
					if(tossJunk){
						if(r1!=null && r1.junk()){
							lowqBasesT+=r1.length();
							lowqReadsT++;
							r1.setDiscarded(true);
						}
						if(r2!=null && r2.junk()){
							lowqBasesT+=r2.length();
							lowqReadsT++;
							r2.setDiscarded(true);
						}
					}
					
					if(iupacToN){
						if(r1!=null){r1.convertUndefinedTo((byte)'N');}
						if(r2!=null){r2.convertUndefinedTo((byte)'N');}
					}
					
					if(remap1!=null && r1!=null){
						int swaps=r1.remapAndCount(remap1);
						if(swaps>0){
							basesSwappedT+=swaps;
							readsSwappedT++;
						}
					}
					if(remap2!=null && r2!=null){
						int swaps=r2.remapAndCount(remap2);
						if(swaps>0){
							basesSwappedT+=swaps;
							readsSwappedT++;
						}
					}
					
					if(trimBadSequence){//Experimental
						if(r1!=null){
							int x=TrimRead.trimBadSequence(r1);
							basesQTrimmedT+=x;
							readsQTrimmedT+=(x>0 ? 1 : 0);
						}
						if(r2!=null){
							int x=TrimRead.trimBadSequence(r2);
							basesQTrimmedT+=x;
							readsQTrimmedT+=(x>0 ? 1 : 0);
						}
					}
					
					if(chastityFilter){
						if(r1!=null && r1.failsChastity()){
							lowqBasesT+=r1.pairLength();
							lowqReadsT+=r1.pairCount();
							r1.setPairDiscarded(true);
						}
					}
					
					if(removeBadBarcodes){
						if(r1!=null && !r1.discarded() && r1.failsBarcode(barcodes, failIfNoBarcode)){
							if(failBadBarcodes){KillSwitch.kill("Invalid barcode detected: "+r1.id+"\nThis can be disabled with the flag barcodefilter=f");}
							lowqBasesT+=r1.pairLength();
							lowqReadsT+=r1.pairCount();
							r1.setPairDiscarded(true);
						}
					}
					
					if(tag!=null && r1!=null && !r1.discarded()) {
						boolean pass=passesTagFilter(r1.id, tag, delimiter);
						if(!pass){r1.setPairDiscarded(true);}
					}
					
					if(filterBits!=0 || requiredBits!=0){
						if(r1!=null && !r1.discarded()){
							assert(sl1!=null) : "filterbits and requiredbits only work on sam/bam input.";
							if(((sl1.flag&filterBits)!=0) || ((sl1.flag&requiredBits)!=requiredBits)){
								r1.setDiscarded(true);
								unmappedBasesT+=initialLength1;
								unmappedReadsT++;
							}
						}
					}
					
					if(minMapq>=0 || maxMapq>=0){
						if(r1!=null && !r1.discarded()){
							assert(sl1!=null) : "mapq filters only work on sam/bam input.";
							final int mapq=sl1.mapped() ? sl1.mapq : 0;
							if((minMapq>=0 && mapq<minMapq) || (maxMapq>=0 && mapq>maxMapq)){
								r1.setDiscarded(true);
								unmappedBasesT+=initialLength1;
								unmappedReadsT++;
							}
						}
					}
					
					if(minIdFilter>=0 || maxIdFilter<1 || USE_EDIT_FILTER){
						if(r1!=null && !r1.discarded()){
							assert(r1.match!=null || r1.samline!=null) : "idfilter requires sam/bam input.";
							boolean pass=passesIDFilter(r1, minIdFilter, maxIdFilter, false);
							if(USE_EDIT_FILTER){
								pass=pass&&passesEditFilter(r1, false);
							}
							if(!pass){
								r1.setDiscarded(true);
								idfilteredBasesT+=initialLength1;
								idfilteredReadsT++;
							}
						}
						if(r2!=null && !r2.discarded()){
							assert(r2.match!=null || r2.samline!=null) : "idfilter requires sam/bam input.";
							boolean pass=passesIDFilter(r2, minIdFilter, maxIdFilter, false);
							if(USE_EDIT_FILTER){
								pass=pass&&passesEditFilter(r2, false);
							}
							if(!pass){
								r2.setDiscarded(true);
								idfilteredBasesT+=initialLength2;
								idfilteredReadsT++;
							}
						}
					}
					
					if(fixCigar){
						if(SamLine.VERSION==1.3f){
							if(r1!=null && !r1.discarded()){
								assert(sl1!=null) : "Cigar string adjustment only works on sam/bam input.";
								sl1.cigar=SamLine.toCigar13(sl1.cigar);
							}
						}else{
							if(r1!=null && !r1.discarded()){
								assert(sl1!=null) : "Cigar string adjustment only works on sam/bam input.";
								if(r1.match!=null){
									r1.toLongMatchString(false);
									int start=sl1.pos-1;
									int stop=start+Read.calcMatchLength(r1.match)-1;
									sl1.cigar=SamLine.toCigar14(r1.match, start, stop, Integer.MAX_VALUE, r1.bases);
								}
							}
						}
					}
					
					if(stoptag){
						if(r1!=null && !r1.discarded()){
							assert(sl1!=null) : "stoptag only works on sam/bam input.";
							if(sl1.mapped() && sl1.cigar!=null){
								if(sl1.optional==null){sl1.optional=new ArrayList<String>(2);}
								sl1.optional.add(SamLine.makeStopTag(sl1.pos, sl1.calcCigarLength(false, false), sl1.cigar, r1.perfect()));
							}
						}
					}
					
					if(pairedOnly || unpairedOnly){
						assert(sl1!=null) : "pairedonly requires sam/bam input.";
						
						if(r1!=null && !r1.discarded() && sl1.properPair()!=pairedOnly){
							r1.setDiscarded(true);
							unmappedBasesT+=initialLength1;
							unmappedReadsT++;
						}
					}
					
					if(mappedOnly || unmappedOnly){
						if(r1!=null && !r1.discarded() && (r1.mapped()!=mappedOnly || r1.bases==null || r1.secondary())){
							r1.setDiscarded(true);
							unmappedBasesT+=initialLength1;
							unmappedReadsT++;
						}
						if(r2!=null && !r2.discarded() && (r2.mapped()!=mappedOnly || r2.bases==null || r2.secondary())){
							r2.setDiscarded(true);
							unmappedBasesT+=initialLength2;
							unmappedReadsT++;
						}
					}
					
					if(primaryOnly){
						if(r1!=null && (r1.bases==null || r1.secondary())){
							r1.setDiscarded(true);
							unmappedBasesT+=initialLength1;
							unmappedReadsT++;
						}
						if(r2!=null && (r2.bases==null || r2.secondary())){
							r2.setDiscarded(true);
							unmappedBasesT+=initialLength2;
							unmappedReadsT++;
						}
					}
					
					if(filterGC && (initialLength1>0 || initialLength2>0)){
						float gc1=(initialLength1>0 ? r1.gc() : -1);
						float gc2=(initialLength2>0 ? r2.gc() : gc1);
						if(gc1==-1){gc1=gc2;}
						if(usePairGC){
							final float gc;
							if(r2==null){
								gc=gc1;
							}else{
								gc=(gc1*initialLength1+gc2*initialLength2)/(initialLength1+initialLength2);
							}
							gc1=gc2=gc;
						}
						if(r1!=null && !r1.discarded() && (gc1<minGC || gc1>maxGC)){
							r1.setDiscarded(true);
							badGcBasesT+=initialLength1;
							badGcReadsT++;
						}
						if(r2!=null && !r2.discarded() && (gc2<minGC || gc2>maxGC)){
							r2.setDiscarded(true);
							badGcBasesT+=initialLength2;
							badGcReadsT++;
						}
					}
					
					if(recalibrateQuality){
						if(r1!=null && !r1.discarded()){
							CalcTrueQuality.recalibrate(r1);
						}
						if(r2!=null && !r2.discarded()){
							CalcTrueQuality.recalibrate(r2);
						}
					}
					
					if(quantizeQuality){
						final byte[] quals1=r1.quality, quals2=(r2==null ? null : r2.quality);
						Quantizer.quantize(quals1);
						Quantizer.quantize(quals2);
					}
					
					if(forceTrimLeft>0 || forceTrimRight>=0 || forceTrimModulo>0 || forceTrimRight2>0){
						if(r1!=null && !r1.discarded()){
							final int len=r1.length();
							final int a=forceTrimLeft>0 ? forceTrimLeft : 0;
							final int b0=forceTrimModulo>0 ? len-1-len%forceTrimModulo : len;
							final int b1=forceTrimRight>=0 ? forceTrimRight : len;
							final int b2=forceTrimRight2>0 ? len-1-forceTrimRight2 : len;
							final int b=Tools.min(b0, b1, b2);
							final int x=TrimRead.trimToPosition(r1, a, b, 1);
							basesFTrimmedT+=x;
							readsFTrimmedT+=(x>0 ? 1 : 0);
							if(r1.length()<minlen1){r1.setDiscarded(true);}
						}
						if(r2!=null && !r2.discarded()){
							final int len=r2.length();
							final int a=forceTrimLeft>0 ? forceTrimLeft : 0;
							final int b0=forceTrimModulo>0 ? len-1-len%forceTrimModulo : len;
							final int b1=forceTrimRight>0 ? forceTrimRight : len;
							final int b2=forceTrimRight2>0 ? len-1-forceTrimRight2 : len;
							final int b=Tools.min(b0, b1, b2);
							final int x=TrimRead.trimToPosition(r2, a, b, 1);
							basesFTrimmedT+=x;
							readsFTrimmedT+=(x>0 ? 1 : 0);
							if(r2.length()<minlen2){r2.setDiscarded(true);}
						}
					}
					
					if(qtrim){
						if(r1!=null && !r1.discarded()){
							int x=TrimRead.trimFast(r1, qtrimLeft, qtrimRight, trimq, trimE, 1);
							basesQTrimmedT+=x;
							readsQTrimmedT+=(x>0 ? 1 : 0);
						}
						if(r2!=null && !r2.discarded()){
							int x=TrimRead.trimFast(r2, qtrimLeft, qtrimRight, trimq, trimE, 1);
							basesQTrimmedT+=x;
							readsQTrimmedT+=(x>0 ? 1 : 0);
						}
					}
					
					if(minAvgQuality>0){
						if(r1!=null && !r1.discarded() && r1.avgQuality(false, minAvgQualityBases)<minAvgQuality){
							lowqBasesT+=r1.length();
							lowqReadsT++;
							r1.setDiscarded(true);
						}
						if(r2!=null && !r2.discarded() && r2.avgQuality(false, minAvgQualityBases)<minAvgQuality){
							lowqBasesT+=r2.length();
							lowqReadsT++;
							r2.setDiscarded(true);
						}
					}
					
					if(maxNs>=0){
						if(r1!=null && !r1.discarded() && r1.countUndefined()>maxNs){
							lowqBasesT+=r1.length();
							lowqReadsT++;
							r1.setDiscarded(true);
						}
						if(r2!=null && !r2.discarded() && r2.countUndefined()>maxNs){
							lowqBasesT+=r2.length();
							lowqReadsT++;
							r2.setDiscarded(true);
						}
					}
					
					if(minConsecutiveBases>0){
						if(r1!=null && !r1.discarded() && !r1.hasMinConsecutiveBases(minConsecutiveBases)){
							lowqBasesT+=r1.length();
							lowqReadsT++;
							r1.setDiscarded(true);
						}
						if(r2!=null && !r2.discarded() && !r2.hasMinConsecutiveBases(minConsecutiveBases)){
							lowqBasesT+=r2.length();
							lowqReadsT++;
							r2.setDiscarded(true);
						}
					}
					
					if(minlen1>0 || minlen2>0 || maxReadLength>0){
//						assert(false) : minlen1+", "+minlen2+", "+maxReadLength+", "+r1.length();
						if(r1!=null && !r1.discarded()){
							int rlen=r1.length();
							if(rlen<minlen1 || (maxReadLength>0 && rlen>maxReadLength)){
								r1.setDiscarded(true);
								readShortDiscardsT++;
								baseShortDiscardsT+=rlen;
							}
						}
						if(r2!=null && !r2.discarded()){
							int rlen=r2.length();
							if(rlen<minlen1 || (maxReadLength>0 && rlen>maxReadLength)){
								r2.setDiscarded(true);
								readShortDiscardsT++;
								baseShortDiscardsT+=rlen;
							}
						}
					}
					
					boolean remove=false;
					if(r2==null){
						remove=r1.discarded();
					}else{
						remove=requireBothBad ? (r1.discarded() && r2.discarded()) : (r1.discarded() || r2.discarded());
					}
					if(invertFilters){
						remove=!remove;
						r1.setDiscarded(!r1.discarded());
						if(r2!=null){r2.setDiscarded(!r2.discarded());}
					}
					
					if(remove){reads.set(idx, null);}
					else{
						
						if(pad){
							pad(r1);
							pad(r2);
						}
						
						if(uniqueNames || addunderscore || addslash || addcolon){

							if(r1.id==null){r1.id=""+r1.numericID;}
							if(r2!=null && r2.id==null){r2.id=r1.id;}

							if(uniqueNames){

								{
									Integer v=nameMap1.get(r1.id);
									if(v==null){
										nameMap1.put(r1.id, 1);
									}else{
										v++;
										nameMap1.put(r1.id, v);
										r1.id=r1.id+"_"+v;
									}
								}
								if(r2!=null){
									Integer v=nameMap2.get(r2.id);
									if(v==null){
										nameMap2.put(r2.id, 1);
									}else{
										v++;
										nameMap2.put(r2.id, v);
										r2.id=r2.id+"_"+v;
									}
								}
							}
							if(addunderscore){
								r1.id=Tools.whitespace.matcher(r1.id).replaceAll("_");
								if(r2!=null){r2.id=Tools.whitespace.matcher(r2.id).replaceAll("_");}
							}
							if(addcolon){
								if(!r1.id.contains(colon1)){r1.id+=colon1;}
								if(r2!=null){
									if(!r2.id.contains(colon2)){r2.id+=colon2;}
								}
							}else if(addslash){
								if(!r1.id.contains(slash1)){r1.id+=slash1;}
								if(r2!=null){
									if(!r2.id.contains(slash2)){r2.id+=slash2;}
								}
							}
						}
					}
					
					if(singles!=null){
						if(r1.discarded() || (r2!=null && r2.discarded())){
							if(!r1.discarded()){
								Read r=r1.clone();
								r.mate=null;
								r.setPairnum(0);
								singles.add(r);
							}else if(r2!=null && !r2.discarded()){
								Read r=r2.clone();
								r.mate=null;
								r.setPairnum(0);
								singles.add(r);
							}
						}
					}
				}
				
				final ArrayList<Read> listOut;
				
//				assert(false) : sampleReadsExact+", "+sampleBasesExact;
				if(sampleReadsExact || sampleBasesExact){
					listOut=new ArrayList<Read>();
					if(sampleReadsExact){
						for(Read r : reads){
							if(r!=null){
								assert(readsRemaining>0) : readsRemaining;
								double prob=sampleReadsTarget/(double)(readsRemaining);
//								outstream.println("sampleReadsTarget="+sampleReadsTarget+", readsRemaining="+readsRemaining+", prob="+prob);
								while(allowUpsample && prob>1){
									listOut.add(r);
									sampleReadsTarget--;
									prob--;
								}
								if(randy.nextDouble()<prob){
									listOut.add(r);
									sampleReadsTarget--;
								}
							}
							readsRemaining--;
						}
					}else if(sampleBasesExact){
						for(Read r : reads){
							if(r!=null){
								assert(basesRemaining>0) : basesRemaining;
								final int bases=r.pairLength();
								double prob=sampleBasesTarget/(double)(basesRemaining);
								while(allowUpsample && prob>1){
									listOut.add(r);
									sampleBasesTarget-=bases;
									prob--;
								}
								if(randy.nextDouble()<prob){
									listOut.add(r);
									sampleBasesTarget-=bases;
								}
								basesRemaining-=bases;
							}
						}
					}
				}else{
					listOut=reads;
				}
//				if(deleteEmptyFiles){
					for(Read r : listOut){
						if(r!=null){
							readsOut1++;
							basesOut1+=r.length();
							if(r.mate!=null){
								readsOut2++;
								basesOut2+=r.mateLength();
							}
						}
					}
					if(singles!=null){
						for(Read r : singles){
							if(r!=null){
								readsOutSingle++;
								basesOutSingle+=r.length();
							}
						}
					}
//				}
				if(ros!=null){ros.add(listOut, ln.id);}
				if(rosb!=null){rosb.add(singles, ln.id);}

				cris.returnList(ln.id, false);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
//				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
				assert(ln.list.isEmpty());
				cris.returnList(ln.id, true);
			}
		}
		
		errorState|=ReadStats.writeAll();
		
		{
			//Prevent a spurious error message in the event of a race condition when maxReads is set.
			boolean b=ReadWrite.closeStream(cris);
			if(maxReads<1 || maxReads==Long.MAX_VALUE || (maxReads!=readsProcessed && maxReads*2!=readsProcessed && samplerate<1)){errorState|=b;}
		}
		errorState|=ReadWrite.closeOutputStreams(ros, rosb);
		
		if(deleteEmptyFiles){
			deleteEmpty(readsOut1, readsOut2, readsOutSingle);
		}
		
//		outstream.println(cris.errorState()+", "+(ros==null ? "null" : (ros.errorState()+", "+ros.finishedSuccessfully())));
//		if(ros!=null){
//			ReadStreamWriter rs1=ros.getRS1();
//			ReadStreamWriter rs2=ros.getRS2();
//			outstream.println(rs1==null ? "null" : rs1.finishedSuccessfully());
//			outstream.println(rs2==null ? "null" : rs2.finishedSuccessfully());
//		}
		
		t.stop();
		
		final long rawReadsIn=cris.readsIn(), rawBasesIn=cris.basesIn();
		final double rmult=100.0/rawReadsIn, bmult=100.0/rawBasesIn;
		final double rpmult=100.0/readsProcessed, bpmult=100.0/basesProcessed;
		
		if(silent){
			if(errorState){
				throw new RuntimeException("ReformatReads terminated in an error state; the output may be corrupt.");
			}
			return;
		}
		
		outstream.println("Input:                  \t"+cris.readsIn()+" reads          \t"+
				cris.basesIn()+" bases");
		if(samplerate!=1f){
			outstream.println("Processed:              \t"+readsProcessed+" reads          \t"+
					basesProcessed+" bases");
		}
		
		if(remap1!=null || remap2!=null){
			outstream.println("Base Transforms:        \t"+readsSwappedT+" reads ("+Tools.format("%.2f",readsSwappedT*rpmult)+"%) \t"+
					basesSwappedT+" bases ("+Tools.format("%.2f",basesSwappedT*bpmult)+"%)");
		}
		if(qtrim || trimBadSequence){
			outstream.println("QTrimmed:               \t"+readsQTrimmedT+" reads ("+Tools.format("%.2f",readsQTrimmedT*rpmult)+"%) \t"+
					basesQTrimmedT+" bases ("+Tools.format("%.2f",basesQTrimmedT*bpmult)+"%)");
		}
		if(forceTrimLeft>0 || forceTrimRight>0 || forceTrimRight2>0 || forceTrimModulo>0){
			outstream.println("FTrimmed:               \t"+readsFTrimmedT+" reads ("+Tools.format("%.2f",readsFTrimmedT*rpmult)+"%) \t"+
					basesFTrimmedT+" bases ("+Tools.format("%.2f",basesFTrimmedT*bpmult)+"%)");
		}
		if(minReadLength>0 || maxReadLength>0){
			outstream.println("Short Read Discards:    \t"+readShortDiscardsT+" reads ("+Tools.format("%.2f",readShortDiscardsT*rpmult)+"%) \t"+
					baseShortDiscardsT+" bases ("+Tools.format("%.2f",baseShortDiscardsT*bpmult)+"%)");
		}
		if(minAvgQuality>0 || maxNs>=0 || chastityFilter || tossJunk || removeBadBarcodes){
			outstream.println("Low quality discards:   \t"+lowqReadsT+" reads ("+Tools.format("%.2f",lowqReadsT*rpmult)+"%) \t"+
					lowqBasesT+" bases ("+Tools.format("%.2f",lowqBasesT*bpmult)+"%)");
		}
		if(minIdFilter>=0 || maxIdFilter<1 || USE_EDIT_FILTER){
			outstream.println("Identity/edit discards: \t"+idfilteredReadsT+" reads ("+Tools.format("%.2f",idfilteredReadsT*rpmult)+"%) \t"+
					idfilteredBasesT+" bases ("+Tools.format("%.2f",idfilteredBasesT*bpmult)+"%)");
		}
		
		if(filterGC){
			outstream.println("GC content discards:    \t"+badGcReadsT+" reads ("+Tools.format("%.2f",badGcReadsT*rpmult)+"%) \t"+
					badGcBasesT+" bases ("+Tools.format("%.2f",badGcBasesT*bpmult)+"%)");
		}
		final long ro=readsOut1+readsOut2+readsOutSingle, bo=basesOut1+basesOut2+basesOutSingle;
		outstream.println("Output:                 \t"+ro+" reads ("+Tools.format("%.2f",ro*rmult)+"%) \t"+
				bo+" bases ("+Tools.format("%.2f",bo*bmult)+"%)");
		
		if(k>0){
			outstream.println(k+"-mers processed:        \t"+kmersProcessed);
			outstream.println("Correct "+k+"-mers:          \t"+Tools.format("%.2f%%", correctKmers*100.0/kmersProcessed));
		}
		if(loglog!=null){
			outstream.println("Unique "+loglog.k+"-mers:           \t"+loglog.cardinality());
		}
		
		if(deleteInput && !errorState && out1!=null && in1!=null){
			try {
				new File(in1).delete();
				if(in2!=null){new File(in2).delete();}
			} catch (Exception e) {
				outstream.println("WARNING: Failed to delete input files.");
			}
		}
		
		outstream.println();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		MetadataWriter.write(null, readsProcessed, basesProcessed, readsOut1+readsOut2, basesOut1+basesOut2, false);
		if(testsize){
			long bytesProcessed=(new File(in1).length()+(in2==null ? 0 : new File(in2).length())+
					(qfin1==null ? 0 : new File(qfin1).length())+(qfin2==null ? 0 : new File(qfin2).length()));//*passes
			double xpnano=bytesProcessed/(double)(t.elapsed);
			String xpstring=(bytesProcessed<100000 ? ""+bytesProcessed : bytesProcessed<100000000 ? (bytesProcessed/1000)+"k" : (bytesProcessed/1000000)+"m");
			while(xpstring.length()<8){xpstring=" "+xpstring;}
			outstream.println("Bytes Processed:    "+xpstring+" \t"+Tools.format("%.2fm bytes/sec", xpnano*1000));
		}
		
		if(verifypairing){
			outstream.println("Names appear to be correctly paired.");
		}
		
		if(outstream!=System.err && outstream!=System.out){outstream.close();}
		
		if(errorState){
			throw new RuntimeException("ReformatReads terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/

	private void pad(Read r){
		pad(r, padLeft, padRight, padSymbol, padQ);
	}
	
	public static final void pad(Read r, int padLeft, int padRight, byte padSymbol, byte padQ){
		if(r==null || r.length()==0 || (padLeft<1 && padRight<1)){return;}
		padLeft=Tools.max(0, padLeft);
		padRight=Tools.max(0, padRight);
		r.bases=pad(r.bases, padLeft, padRight, padSymbol);
		r.quality=pad(r.quality, padLeft, padRight, padQ);
	}
	
	private static final byte[] pad(byte[] old, int padLeft, int padRight, byte padSymbol){
		if(old==null){return null;}
		final int innerLimit=old.length+padLeft;
		byte[] array=new byte[innerLimit+padRight];
		Arrays.fill(array, 0, padLeft, padSymbol);
		Arrays.fill(array, innerLimit, innerLimit+padRight, padSymbol);
		for(int i=0; i<old.length; i++){array[i+padLeft]=old[i];}
		return array;
	}
	
	public static final boolean passesIDFilter(Read r, float minId, float maxId, boolean requireMapped){
		if(!passesMinIDFilter(r, minId, requireMapped)){return false;}
		return passesMaxIDFilter(r, maxId);
	}
	
	public static final boolean passesMinIDFilter(Read r, float minId, boolean requireMapped){
		if(minId<=0 || r.perfect()){return true;}
		if(r.match==null && r.samline!=null){
			r.match=r.samline.toShortMatch(false);
		}
		if(r.match==null){return !requireMapped;}
		return Read.identityFlat(r.match, true)>=minId;
	}
	
	public static final boolean passesMaxIDFilter(Read r, float maxId){
		if(maxId>=1){return true;}
		if(r.match==null && r.samline!=null){
			r.match=r.samline.toShortMatch(false);
		}
		if(r.match==null){return true;}
		return Read.identityFlat(r.match, true)<=maxId;
	}
	
	public final boolean passesEditFilter(Read r, boolean requireMapped){
		if(r.perfect()){return true;}
		if(r.match==null && r.samline!=null){
			r.match=r.samline.toShortMatch(false);
		}
		if(r.match==null){return !requireMapped;}
		r.toLongMatchString(false);
		
		final int sub=Read.countSubs(r.match);
		final int ins=Read.countInsertions(r.match);
		final int del=Read.countDeletions(r.match);
		final int inscount=Read.countInsertionEvents(r.match);
		final int delcount=Read.countDeletionEvents(r.match);
		final int clip=SamLine.countLeadingClip(r.match)+SamLine.countTrailingClip(r.match);
		
		boolean bad=false;
		bad=bad||(subfilter>=0 && sub>subfilter);
		bad=bad||(minsubs>=0 && sub<minsubs);
		bad=bad||(clipfilter>=0 && clip>clipfilter);
		bad=bad||(insfilter>=0 && inscount>insfilter);
		bad=bad||(delfilter>=0 && delcount>delfilter);
		bad=bad||(inslenfilter>=0 && r.hasLongInsertion(inslenfilter));
		bad=bad||(dellenfilter>=0 && r.hasLongDeletion(dellenfilter));
		bad=bad||(indelfilter>=0 && inscount+delcount>indelfilter);
		bad=bad||(editfilter>=0 && sub+ins+del>editfilter);
		
		return !bad;
	}
	
	public boolean passesTagFilter(String s, String tag, char delimiter) {
		if(requiredValue!=null) {
			String value=Parse.parseString(s, tag, delimiter);
			return (value==null ? requiredValue==null : value.equals(requiredValue));
		}
		double value=Parse.parseDouble(s, tag, delimiter);
		return value>=minValue && value<=maxValue;
	}
	
//	public boolean passesTagFilter(String s, Pattern splitter) {
//		String[] split=s.split(s);
//		if(requiredValue!=null) {
//			String value=Parse.parseString(s, tag, delimiter);
//			return (value==null ? requiredValue==null : value.equals(requiredValue));
//		}
//		double value=Parse.parseDouble(s, tag, delimiter);
//		return value>=minValue && value<=maxValue;
//	}
	
	/*--------------------------------------------------------------*/
	
	private void deleteEmpty(long readsOut1, long readsOut2, long readsOutSingle){
		deleteEmpty(readsOut1, ffout1, qfout1);
		deleteEmpty(readsOut2, ffout2, qfout2);
		deleteEmpty(readsOutSingle, ffoutsingle, null);
	}
	
	private static void deleteEmpty(long count, FileFormat ff, String qf){
		try {
			if(ff!=null && count<1){
				String s=ff.name();
				if(s!=null && !ff.stdio() && !ff.devnull()){
					File f=new File(ff.name());
					if(f.exists()){
						f.delete();
					}
				}
				if(qf!=null){
					File f=new File(qf);
					if(f.exists()){
						f.delete();
					}
				}
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private long[] countReads(long maxReads){
		if(ffin1.stdio()){
			throw new RuntimeException("Can't precount reads from standard in, only from a file.");
		}
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			if(verbose){outstream.println("Counting Reads");}
			cris.start(); //4567
		}
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		long count=0, count2=0, bases=0;
		
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			count+=reads.size();
			for(Read r : reads){
				bases+=r.length();
				count2++;
				if(r.mate!=null){
					bases+=r.mateLength();
					count2++;
				}
			}
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln);
		errorState|=ReadWrite.closeStream(cris);
		return new long[] {count, count2, bases};
	}
	
	private SuperLongList makeLengthHist(long maxReads){
		if(ffin1.stdio()){
			throw new RuntimeException("Can't precount reads from standard in, only from a file.");
		}
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			if(verbose){outstream.println("Counting Reads");}
			cris.start(); //4567
		}
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		SuperLongList sll=new SuperLongList(200000);
		
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			for(Read r : reads){
				sll.add(r.length());
				if(r.mate!=null){
					sll.add(r.mateLength());
				}
			}
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln);
		errorState|=ReadWrite.closeStream(cris);
		sll.sort();
		return sll;
	}
	
	public void setSampleSeed(long seed){
		randy=Shared.threadLocalRandom(seed);
	}

	
	/*--------------------------------------------------------------*/
	
	public static final void underscoreToSpace(Read r){
		if(r!=null){
			r.id=underscoreToSpace(r.id);
		} 
	}
	
	public static String underscoreToSpace(String header){
		int x=header.lastIndexOf('_');
		if(x<0){return header;}
		byte[] bytes=header.getBytes();
		bytes[x]=' ';
		return new String(bytes);
	}
	
	public static final void fixHeader(Read r){
		if(r!=null){
			r.id=fixHeader(r.id);
			if(r.samline!=null){
				r.samline.qname=r.id;
			}
		} 
	}
	
	public static final String fixHeader(String header){
		if(header==null || header.length()<1){return header;}
		byte[] array=new byte[header.length()];
		boolean changed=false;
		for(int i=0; i<header.length(); i++){
			char c=header.charAt(i);
			byte b=headerSymbols[c];
			array[i]=b;
			changed|=(b!=c);
		}
		if(changed){header=new String(array);}
		return header;
	}
	
	/** For converting headers to filesystem-valid Strings */
	private static final byte[] headerSymbols=new byte[128];
	
	static{
		Arrays.fill(headerSymbols, (byte)'_');
		for(int i=0; i<128; i++){
			if(Character.isLetterOrDigit(i)){
				headerSymbols[i]=(byte)i;
			}
		}
		char[] acceptable=new char[] {'_', '.', '#', '-', '(', ')', '~'};
		for(char c : acceptable){
			headerSymbols[c]=(byte)c;
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	private String out1=null;
	private String out2=null;
	private String outsingle=null;

	private String qfout1=null;
	private String qfout2=null;
	
	private String extin=null;
	private String extout=null;
	
	private String ref=null;
	
	/*--------------------------------------------------------------*/
	
	public long readsProcessed=0;
	public long basesProcessed=0;
	public long kmersProcessed=0;
	public double correctKmers=0;
	
	//Only used with deleteEmptyFiles flag
	public long readsOut1=0;
	public long readsOut2=0;
	public long readsOutSingle=0;
	
	public long basesOut1=0;
	public long basesOut2=0;
	public long basesOutSingle=0;
	
	/*--------------------------------------------------------------*/
	
	/** For calculating kmer cardinality */
	private final CardinalityTracker loglog;
	
	private int k=0;
	
	/** Tracks names to ensure no duplicate names. */
	private final HashMap<String,Integer> nameMap1, nameMap2;
	private boolean uniqueNames=false;

	private boolean reverseComplementMate=false;
	private boolean reverseComplement=false;
	private boolean complement=false;
	private boolean verifyinterleaving=false;
	private boolean verifypairing=false;
	private boolean allowIdenticalPairNames=true;
	private boolean trimBadSequence=false;
	private boolean chastityFilter=false;
	/** Crash if a barcode is encountered that contains Ns or is not in the table */
	private final boolean failBadBarcodes;
	/** Remove reads with Ns in barcodes or that are not in the table */
	private final boolean removeBadBarcodes;
	/** Fail reads missing a barcode */
	private final boolean failIfNoBarcode;
	/** A set of valid barcodes; null if unused */
	private final HashSet<String> barcodes;
	private boolean deleteEmptyFiles=false;
	private boolean mappedOnly=false;
	private boolean pairedOnly=false;
	private boolean unpairedOnly=false;
	private boolean unmappedOnly=false;
	private boolean primaryOnly=false;
	/** For sam file filtering: These bits must be set. */
	private int requiredBits=0;
	/** For sam file filtering: These bits must be unset */
	private int filterBits=0;
	/** Add /1 and /2 to read names */
	private boolean addslash=false;
	/** Add 1: and 2: to read names */
	private boolean addcolon=false;
	/** Change read name whitespace to underscores */
	private boolean addunderscore=false;
	private boolean stoptag=false;
	private boolean iupacToN=false;
	
	boolean bottom=true;
	boolean top=true;

	private float minIdFilter=-1;
	private float maxIdFilter=-1;
	private int subfilter=-1;
	private int minsubs=-1;
	private int clipfilter=-1;
	private int delfilter=-1;
	private int insfilter=-1;
	private int indelfilter=-1;
	private int dellenfilter=-1;
	private int inslenfilter=-1;
	private int editfilter=-1;
	private int minMapq=-1;
	private int maxMapq=-1;
	
	boolean USE_EDIT_FILTER=false;

	private long maxReads=-1;
	private long skipreads=-1;
	private float samplerate=1f;
	private long sampleseed=-1;
	private boolean sampleReadsExact=false;
	private boolean sampleBasesExact=false;
	private boolean allowUpsample=false;
	private boolean prioritizeLength=false;
	private long sampleReadsTarget=0;
	private long sampleBasesTarget=0;
	
	/** Recalibrate quality scores using matrices */
	private boolean recalibrateQuality=false;
	private boolean qtrimRight=false;
	private boolean qtrimLeft=false;
	private final int forceTrimLeft;
	private final int forceTrimRight;
	private final int forceTrimRight2;
	/** Trim right bases of the read modulo this value.
	 * e.g. forceTrimModulo=50 would trim the last 3bp from a 153bp read. */
	private final int forceTrimModulo;
	private float trimq=6;
	/** Error rate for trimming (derived from trimq) */
	private final float trimE;
	private float minAvgQuality=0;
	private int minAvgQualityBases=0;
	private int maxNs=-1;
	private int minConsecutiveBases=0;
	private int breakLength=0;
	private int maxReadLength=0;
	private int minReadLength=0;
	private float minLenFraction=0;
	private float minGC=0;
	private float maxGC=1;
	private boolean filterGC=false;
	/** Average GC for paired reads. */
	private boolean usePairGC;
	private boolean tossJunk=false;
	/** Toss pair only if both reads are shorter than limit */
	private boolean requireBothBad=false;
	/** Invert filters. */
	private boolean invertFilters;
	
	/** Replace problematic symbols in read headers. */
	private boolean fixHeaders;
	
	/** Change the last underscore to a space, to fix UCSF headers */
	private boolean lastUnderscoreToSpace=false;
	
	private boolean useSharedHeader;
	
	private byte[] remap1=null, remap2=null;

	private boolean quantizeQuality=false;
	private boolean fixCigar=false;
	
	private final boolean pad;
	private byte padSymbol='N';
	private byte padQ=0;
	private int padLeft=0;
	private int padRight=0;
	
	//For scanning for specific terms in headers, like "count=53".
	/** Parse headers for this tag, like "count=" */
	private String tag;
	/** Use this delimiter when parsing the header */
	private char delimiter=0;
	private String delimiterS=null;
	/** Min value for numeric terms */
	private float minValue=Float.MIN_VALUE;
	/** Max value for numeric terms */
	private float maxValue=Float.MAX_VALUE;
	/** A specific String value for terms */
	private String requiredValue=null;
//	private final LineParser2S lp;
//	private Pattern splitter;
	
	ReadHeaderParser rhp=new IlluminaHeaderParser2();
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;

	private final FileFormat ffout1;
	private final FileFormat ffout2;
	private final FileFormat ffoutsingle;
	
	private final boolean qtrim;
	
	
	/*--------------------------------------------------------------*/
	
	private static String slash1=" /1";
	private static String slash2=" /2";
	private static final String colon1=" 1:";
	private static final String colon2=" 2:";
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean silent=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	private boolean testsize=false;
	private boolean deleteInput=false;
	
	private Random randy;
	
}
