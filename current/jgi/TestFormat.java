package jgi;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import align2.BandedAlignerConcrete;
import align2.QualityTools;
import barcode.Barcode;
import barcode.BarcodeStats;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import json.JsonObject;
import json.JsonParser;
import server.ServerTools;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import sketch.DisplayParams;
import sketch.SendSketch;
import sketch.Sketch;
import sketch.SketchMakerMini;
import sketch.SketchObject;
import sketch.SketchTool;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.Read;
import structures.ByteBuilder;
import structures.IntHashMapBinary;
import structures.ListNum;
import structures.LongPair;
import structures.Range;
import structures.SuperLongList;
import tracker.AdapterTracker;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Jan 6, 2018
 *
 */
public class TestFormat {

	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		TestFormat x=new TestFormat(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public TestFormat(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Read.JUNK_MODE=Read.FLAG_JUNK;
		Read.CHANGE_QUALITY=false;
		Read.NULLIFY_BROKEN_QUALITY=true;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		SketchObject.defaultParams.minProb=0;
		FileFormat.PRINT_WARNING=false;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("full")){
				full=Parse.parseBoolean(b);
			}else if(a.equals("in") || a.equals("in1")){
				in.add(b);
			}else if(a.equals("sketchsize")){
				sketchSize=Parse.parseIntKMG(b);
			}
			
			else if(a.equalsIgnoreCase("parseZmw") || a.equals("zmw") || a.equals("zmws")){
				parseZmw=Parse.parseBoolean(b);
			}
			
			else if(a.equals("barcodes") || a.equals("barcodefile") || a.equals("indexes") || a.equals("indices")){
				barcodeFile=b;
			}
			
			else if(a.equals("expectedbarcodes") || a.equals("expectedindexes") || a.equals("expectedindices") || a.equals("barcodelist") || a.equals("indexlist")){
				expectedCodeFile=b;
			}
			
			else if(a.equals("hist") || a.equals("hists")){
				if(!Parse.parseBoolean(b)){
					qhistFile=ihistFile=khistFile=bhistFile=lhistFile=gchistFile=zmwhistFile=null;
				}
			}else if(a.equals("qhist") || a.equals("qhistfile")){
				qhistFile=b;
			}else if(a.equals("ihist") || a.equals("ihistfile")){
				ihistFile=b;
			}else if(a.equals("khist") || a.equals("khistfile")){
				khistFile=b;
			}else if(a.equals("bhist") || a.equals("bhistFile")){
				if(b==null || b.equalsIgnoreCase("f") || b.equalsIgnoreCase("false")){
					bhistFile=null;
//					forceMakeBhist=false;
				}else if(b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true")){
//					forceMakeBhist=true;
				}else{
					bhistFile=b;
//					forceMakeBhist=true;
				}
			}else if(a.equals("maxbhistlen") || a.equals("bhistlen")){
				maxBhistLen=Parse.parseIntKMG(b);
			}else if(a.equals("lhist") || a.equals("lhistfile")){
				lhistFile=b;
			}else if(a.equals("zmwhist") || a.equals("zmwhistfile")){
				zmwhistFile=b;
			}else if(a.equals("gchist") || a.equals("gchistfile")){
				gchistFile=b;
			}else if(a.equals("junk") || a.equals("junkfile")){
				junkFile=b;
			}else if(a.equalsIgnoreCase("printBarcodes")){
				printBarcodes=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("printQhist")){
				printQhist=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("printIhist")){
				printIhist=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("printJunk") || a.equals("junk")){
				printJunk=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("printSpeed") || a.equals("speed")){
				printSpeed=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("fast")){
				fast=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("slow")){
				fast=!Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("edist") || a.equalsIgnoreCase("calcedist") || a.equalsIgnoreCase("edits")){
				BarcodeStats.calcEdist=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("sketch") || a.equalsIgnoreCase("card")){
				makeSketch=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("merge")){
				doMerge=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("trim")){
				doTrim=Parse.parseBoolean(b);
			}
			
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				in.add(arg);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();

//			assert(false) : FASTQ.DETECT_QUALITY+", "+FASTQ.IGNORE_BAD_QUALITY+", "+FASTQ.ASCII_OFFSET;
//			System.err.println(FASTQ.DETECT_QUALITY+", "+FASTQ.IGNORE_BAD_QUALITY+", "+FASTQ.ASCII_OFFSET);
			
			maxReads=parser.maxReads;
			makeSketch=(makeSketch || parser.loglog);
			SketchObject.targetSketchSize=sketchSize;
//			loglog=(parser.loglog ? CardinalityTracker.makeTracker(parser) : null);
		}
//		assert(false) : sketch+", "+parser.loglog;
		if(makeSketch){
			SketchObject.AUTOSIZE=false;
//			SketchObject.defaultParams.minKeyOccuranceCount=2;
			SketchObject.defaultParams.parse("trackcounts", "trackcounts", null);
//			SketchObject.defaultParams.minProb=0;
			SketchObject.postParse();
			tool=new SketchTool(sketchSize, SketchObject.defaultParams);
			smm=new SketchMakerMini(tool, SketchObject.ONE_SKETCH, SketchObject.defaultParams);
		}else{
			tool=null;
			smm=null;
		}
		
		makeBhist=bhistFile!=null;
		makeLhist=lhistFile!=null;//TODO
		makeGChist=gchistFile!=null;
		ReadStats.COLLECT_BASE_STATS=makeBhist;
		ReadStats.COLLECT_GC_STATS=makeGChist;
		ReadStats.BASE_HIST_FILE=bhistFile;
		ReadStats.GC_HIST_FILE=gchistFile;
		ReadStats.GC_BINS=4000;
		ReadStats.GC_BINS_AUTO=true;
		
		initialQin=FASTQ.ASCII_OFFSET;
		initialDetectQuality=FASTQ.DETECT_QUALITY;
		
		if(barcodeFile!=null){
			if(barcodeFile.indexOf('#')>0){
//				barcodeFile1=barcodeFile.replace('#', '1');
//				barcodeFile2=barcodeFile.replace('#', '2');
				barcodeFile=barcodeFile.replace('#', 's');
			}
		}
		
		//TODO: in is not considered here.
		if(!Tools.testForDuplicateFiles(true, qhistFile, ihistFile, khistFile, bhistFile, lhistFile, gchistFile, zmwhistFile, 
				barcodeFile, /*barcodeFile1, barcodeFile2,*/ expectedCodeFile, junkFile)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	
	
	void process(Timer t){
		boolean sequence=false, variant=false;
		for(String fname : in){
			final FileFormat ff=test(fname);
			if(full){
				if(ff.isSequence()){
					sequence=true;
					processReads(ff);
				}else if(ff.var()){
					variant=true;
					loadVars(ff);
				}else if(ff.vcf()){
					variant=true;
					loadVcf(ff);
				}else{
					System.err.println("Does not seem to be a sequence or variant format: "+ff.rawExtension());
				}
			}
//			System.err.println(ff);
		}
		
		
		if(sequence){
			printSequenceResults();
		}else if(variant){
			printVariantResults();
		}
		
		t.stop();
		if(printSpeed){
			outstream.println("Time:                         \t"+t);
			if(sequence){
				outstream.println("Reads Processed:    "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
			}else if(variant){
				outstream.println("Vars Processed:     "+variantsProcessed+" \t"+Tools.format("%.2fk vars/sec", (variantsProcessed/(double)(t.elapsed))*1000000));
			}
		}
	}
	
	void printVariantResults(){
		println("Format\t\t"+FileFormat.FORMAT_ARRAY[format]);
		println("Compression\t"+FileFormat.COMPRESSION_ARRAY[compression]);
		println("HeaderLines\t"+headerLinesProcessed);
		println("Variants\t"+variantsProcessed);
		if(ploidy>0){println("Ploidy\t\t"+ploidy);}
		if(pairingRate>0){println("PairingRate\t"+Tools.format("%.4f", pairingRate));}
		if(mapqAvg>0){println("MapqAvg\t\t"+Tools.format("%.2f", mapqAvg));}
		if(totalQualityAvg>0){println("QualityAvg\t\t"+Tools.format("%.2f", totalQualityAvg));}
		if(readLengthAvg>0){println("ReadLengthAvg\t"+Tools.format("%.2f", readLengthAvg));}
	}
	
	void printSequenceResults(){

		println("Format\t\t"+FileFormat.FORMAT_ARRAY[format]);
		println("Compression\t"+FileFormat.COMPRESSION_ARRAY[compression]);
		println("Interleaved\t"+interleaved);
		println("MaxLen\t\t"+maxLen);
		println("MinLen\t\t"+(minLen<Integer.MAX_VALUE ? minLen : 0));
		println("AvgLen\t\t"+Tools.format("%.2f",basesProcessed/Tools.max(1.0, readsProcessed)));
		sll.sort();
		println("StdevLen\t"+Tools.format("%.2f",sll.stdev()));
		println("ModeLen\t\t"+sll.mode());
		
		if(format!=FileFormat.FASTA && format!=FileFormat.UNKNOWN){
			println("QualOffset\t"+offset);
			long negatives=shared.Vector.sum(qhist, 0, 127);
			println("NegativeQuals\t"+negatives);
		}
		if(differs){}
		
		if(!full){return;}
		
		amino=acidsNotBasesProcessed*8>basesProcessed;
//		System.out.println(acidsProcessed+","+basesProcessed);
		println("");
		
		ReadStats.overwrite=true;
		errorState|=ReadStats.writeAll();
		if(amino){
			printAminoTop();
		}else{
			printNucleotideTop();
		}
		
		if(shared.Vector.sum(qhist)>0){
			printQhist();
		}
		if(doMerge && pairsProcessed>0 && !amino && shared.Vector.sum(ihist)>0){
			printIhist();
		}
		if(!amino && barcodeStats.size()>0){
			printBarcodes();
		}
		if(!amino && aTrack.makeSequence()){
			aTrack.writeAdapterConsensus(adapterFile);
		}
		if(!amino && junkProcessed>0){
			printJunk();
		}
	}
	
	void printAminoTop(){
		println("Content\t\tAminoAcids");
		println("FirstHeader\t"+firstHeader);
		println("Sequences\t"+readsProcessed);
		println("Residues\t"+basesProcessed);
		println("-Lowercase\t"+lowerUpperSymbol[0]);
		println("-Uppercase\t"+lowerUpperSymbol[1]);
		println("-Non-Letter\t"+lowerUpperSymbol[2]);
		println("-FullyDefined\t"+AXEGO[0]);
		println("-Stop\t\t"+AXEGO[2]);
		println("-No-call\t"+AXEGO[1]);
		println("-Gap\t\t"+AXEGO[3]);
		println("-Invalid\t"+AXEGO[4]);
	}
	
	void printNucleotideTop(){

		long GC=ACGTUNIGO[1]+ACGTUNIGO[2];
		long ATU=ACGTUNIGO[0]+ACGTUNIGO[3]+ACGTUNIGO[4];
		long T=ACGTUNIGO[3];
		long U=ACGTUNIGO[4];
		long N=ACGTUNIGO[5];
		long I=ACGTUNIGO[6];
		long G=ACGTUNIGO[7];
		long O=ACGTUNIGO[8];
		
		println("Content\t\tNucleotides");
		println("Type\t\t"+(U==0 ? "DNA" : T==0 ? "RNA" : "Mixed"));
		println("FirstHeader\t"+firstHeader);
		println("Reads\t\t"+readsProcessed);
		println("-JunkReads\t"+junkProcessed);
		println("-ChastityFail\t"+chastityFail);
		println("-BadPairNames\t"+badPairs);
		println("");
		println("Bases\t\t"+basesProcessed);
		println("-Lowercase\t"+lowerUpperSymbol[0]);
		println("-Uppercase\t"+lowerUpperSymbol[1]);
		println("-Non-Letter\t"+lowerUpperSymbol[2]);
		println("-FullyDefined\t"+(GC+ATU));
		println("-No-call\t"+(N));
		println("-Degenerate\t"+(I));
		println("-Gap\t\t"+(G));
		println("-Invalid\t"+(O));
		println("");
		println("GC\t\t"+Tools.format("%.3f", GC*1.0/(GC+ATU)));
		if(makeGChist){
			println("-GCMedian\t"+Tools.format("%.3f", ReadStats.GCMedian));
			println("-GCMode\t\t"+Tools.format("%.3f", ReadStats.GCMode));
			println("-GCSTDev\t"+Tools.format("%.3f", ReadStats.GCSTDev));
			println("");
		}
		
//		if(loglog!=null){println("Cardinality\t"+loglog.cardinality());}
		if(smm!=null){
			sketch=smm.toSketch(smm.pacBioDetected ? 2 : 1);
//			System.err.println(smm.pacBioDetected);
			println("Cardinality\t"+(sketch==null ? 0 : sketch.genomeSizeEstimate()));
			println("Cardinality_3+\t"+(sketch==null ? 0 : sketch.genomeSizeEstimate(3)));
			if(khistFile!=null){
				if(sketch!=null){
					ArrayList<LongPair> list=sketch.toKhist();
					ByteStreamWriter bsw=new ByteStreamWriter(khistFile, true, false, false);
					bsw.start();
					bsw.println("#Depth\tCount");
					for(LongPair lp : list){
						bsw.print(lp.a).print('\t').print(lp.b).nl();
					}
					bsw.poisonAndWait();
				}else{
					try {
						if(!FileFormat.isStdio(khistFile) && new File(khistFile).exists()){
							new File(khistFile).delete();
						}
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			if(zmwhistFile!=null){
				if(parseZmw && zmwMap.size()>0){
					IntHashMapBinary counts=zmwMap.toCountHistogram();
					int[] keys=counts.toKeyArray();
					Arrays.sort(keys);
					ByteStreamWriter bsw=new ByteStreamWriter(zmwhistFile, true, false, false);
					bsw.start();
					bsw.println("#Passes\tCount");
					for(int i=0, prev=0; i<keys.length; i++){
						int key=keys[i];
						int value=counts.get(key);
						while(prev<key-1){
							prev++;
							bsw.print(prev).print('\t').print(0).nl();
						}
						bsw.print(key).print('\t').print(value).nl();
						prev=key;
					}
					bsw.poisonAndWait();
					
				}else{
					ReadWrite.delete(zmwhistFile, false);
				}
			}
			if(sketch!=null){
				ServerTools.suppressErrors=true;
				String results=SendSketch.sendSketch(sketch, "refseq", DisplayParams.FORMAT_JSON, 0);
//				assert(results!=null) : results+", "+sketch.toString();
				if(results!=null){
					JsonObject all=JsonParser.parseJsonObjectStatic(results);
//					System.err.println(all);
					if(all!=null && all.jmapSize()>0){
						JsonObject topHit=null;
						for(String key : all.jmap.keySet()){
							JsonObject hit=all.jmap.get(key);
							topHit=hit;
							break;
						}
						println("Organism\t"+topHit.getString("taxName"));
						println("TaxID   \t"+topHit.getLong("TaxID"));
					}
				}
			}
		}
		println("Barcodes\t"+barcodeStats.size());
		println("ZMWs    \t"+zmwMap.cardinality());
		if(doMerge && pairsProcessed>0){
			final long numMerged=shared.Vector.sum(ihist);
			final double insertAvg=Tools.averageHistogram(ihist);
			final int insertMode=Tools.maxIndex(ihist);
			final double mergeFraction=(numMerged/(1.0*Tools.max(mergeAttempts, 1)));
			final double adapterBaseFraction=(adapterBases/(1.0*Tools.max(basesProcessed, 1)));
			final double adapterReadFraction=(adapterReads/(1.0*Tools.max(readsProcessed, 1)));
			println("\nMergable\t"+Tools.format("%.2f%%", 100*mergeFraction));
			if(mergeFraction>0.02){
				println("-InsertMean\t"+Tools.format("%.2f", insertAvg));
				println("-InsertMode\t"+insertMode);
				println("-AdapterReads\t"+Tools.format("%.3f%%", 100*adapterReadFraction));
				println("-AdapterBases\t"+Tools.format("%.3f%%", 100*adapterBaseFraction));
				if(aTrack.makeSequence()){
					println("-Adapter1\t"+aTrack.seq1);
					println("-Adapter2\t"+aTrack.seq2);
				}
			}
		}
	}
	
	void printQhist(){
		long qSum=0;
		double errorSum=0;
		long qCalled=0;
		for(int q=1, qo=1+qOffset; qo<qhist.length; q++, qo++){
			long count=qhist[qo];
			qCalled+=count;
			qSum+=q*count;
			errorSum+=(count*QualityTools.PROB_ERROR[q]);
		}
		qCalled=Tools.max(1, qCalled);
		double avg=qSum/qCalled;
		double errorAvg=errorSum/qCalled;
		double logAvg=QualityTools.probErrorToPhredDouble(errorAvg);
		double trimMult=100.0/(Tools.max(basesProcessed, 1));
		println("\nQErrorRate\t"+Tools.format("%.3f%%", 100*errorAvg));
		println("-QAvgLog\t"+Tools.format("%.2f", logAvg));
		println("-QAvgLinear\t"+Tools.format("%.2f", avg));
		println("-qMinUncalled\t"+qMinUncalled);
		println("-qMaxUncalled\t"+qMaxUncalled);
		println("-qMinCalled\t"+qMinCalled);
		println("-qMaxCalled\t"+qMaxCalled);
		println("-qScoreRange\t"+qhistToRangeString(qhist));
		
		if(doTrim){
			println("-TrimmedAtQ5\t"+Tools.format("%.2f%%", trimhist[5]*trimMult));
			println("-TrimmedAtQ10\t"+Tools.format("%.2f%%", trimhist[10]*trimMult));
			println("-TrimmedAtQ15\t"+Tools.format("%.2f%%", trimhist[15]*trimMult));
			println("-TrimmedAtQ20\t"+Tools.format("%.2f%%", trimhist[20]*trimMult));
//			println("-TrimmedAtQ25\t"+Tools.format("%.2f%%", trimhist[25]*trimMult));
//			println("-TrimmedAtQ30\t"+Tools.format("%.2f%%", trimhist[30]*trimMult));
		}

		if(printQhist){
			println("\nQhist:");
			for(int i=0; i<qhist.length; i++){
				long q=qhist[i];
				if(q>0){println((i-qOffset)+"\t\t"+q);}
			}
		}
		if(qhistFile!=null && shared.Vector.sum(qhist)>0){
			try {
				StringBuilder sb=new StringBuilder();
				sb.append("#QErrorRate\t"+Tools.format("%.3f%%\n", 100*errorAvg));
				sb.append("#QAvgLog\t"+Tools.format("%.2f\n", logAvg));
				sb.append("#QAvgLinear\t"+Tools.format("%.2f", avg));
				printToFileOffset(qhist, true, sb.toString(), qhistFile, qOffset);
			} catch (Throwable e) {
				System.err.println("ERROR - Could not write qhist: "+e.toString());
				errorState=true;
			}
		}
	}
	
	static String qhistToRangeString(long[] qhist) {
		ArrayList<Range> list=qhistToRanges(qhist);
		if(list.isEmpty()) {return "";}
		ByteBuilder bb=new ByteBuilder();
		for(Range r : list){
			if(!bb.isEmpty()){bb.comma();}
			r.appendTo(bb);
		}
		return bb.toString();
	}
	
	static ArrayList<Range> qhistToRanges(long[] qhist){
		ArrayList<Range> list=new ArrayList<Range>();
		Range current=null;
		for(int i=0; i<qhist.length; i++){
			final long count=qhist[i];
			final int q=i-128;
			if(count>0){
				if(current==null){
					current=new Range(q, q);
					list.add(current);
				}else if(current.touches(q)){
					current.absorb(q);
				}else{
					current=new Range(q, q);
					list.add(current);
				}
			}
		}
		return list;
	}
	
	void printIhist(){

		final long numMerged=shared.Vector.sum(ihist);
		final double insertAvg=Tools.averageHistogram(ihist);
		final int insertMode=Tools.maxIndex(ihist);
		final double mergeFraction=(numMerged/(1.0*Tools.max(mergeAttempts, 1)));
		final double adapterBaseFraction=(adapterBases/(1.0*Tools.max(basesProcessed, 1)));
		final double adapterReadFraction=(adapterReads/(1.0*Tools.max(readsProcessed, 1)));
		if(printIhist){
			println("\nIhist:");
			for(int i=0; i<ihist.length; i++){
				long q=ihist[i];
				if(q>0){println(i+"\t\t"+q);}
			}
		}

		if(ihistFile!=null){
			try {
				StringBuilder sb=new StringBuilder();
				sb.append("#InsertMean\t"+Tools.format("%.2f\n", insertAvg));
				sb.append("#InsertMode\t"+insertMode+"\n");
				sb.append("#AdapterReads\t"+Tools.format("%.2f%%\n", 100*adapterReadFraction));
				sb.append("#AdapterBases\t"+Tools.format("%.2f%%\n", 100*adapterBaseFraction));
				printToFile(ihist, true, sb.toString(), ihistFile);
			} catch (Throwable e) {
				System.err.println("ERROR - Could not write ihist: "+e.toString());
				errorState=true;
			}
		}
	}
	
	void printBarcodes(){
		final BarcodeStats bs=barcodeStats;
		
		if(barcodeFile!=null){
			bs.printToFile(barcodeFile, true);
		}
		
		println("");
		if(barcodesPerRead>1) {
			bs.leftStats=bs.makeLeft();
			bs.rightStats=bs.makeRight();
			bs.calcStats();
			println(bs.toStats("Barcodes"));
			bs.leftStats.calcStats();
			println(bs.leftStats.toStats("Barcode1"));
			bs.rightStats.calcStats();
			println(bs.rightStats.toStats("Barcode2"));
		}else{
			bs.calcStats();
			println(bs.toStats(""));
		}
		
		if(printBarcodes){
			ArrayList<Barcode> list=bs.toList();
			println("\nBarcodeList:");
			for(Barcode bc : list){println(bc);}
		}
	}
	
	void printJunk(){
		if(printJunk){
			println("\nJunkList:");
			for(String s : invalidHeaders){
				println(s);
			}
		}
		
		if(junkFile!=null){
			try {
				TextStreamWriter tsw=new TextStreamWriter(junkFile, true, false, false);
				tsw.start();
				for(String s : invalidHeaders){
					tsw.println(s);
				}
				errorState|=tsw.poisonAndWait();
			} catch (Throwable e) {
				System.err.println("ERROR - Could not write junk file: "+e.toString());
				errorState=true;
			}
		}
	}
	
	void println(Object o){System.out.println(o);}
	
	private FileFormat test(String fname){
		FileFormat ffName=FileFormat.testInput(fname, FileFormat.FASTQ, null, false, false, false);
		FileFormat ffContent=FileFormat.testInput(fname, ffName.format(), null, true, true, true);
		FileFormat ff=ffContent;
		
//		offset=33;
//		maxLen=0;
//		minLen=Integer.MAX_VALUE;
//		interleaved=false;
		
//		assert(false) : ffName+"\n"+ffContent;
		if(ff==null){
			System.out.println("null");
		}else{
			format=ff.format();
			compression=ff.compression();
			if(ff.fastq()){
				byte qold=stream.FASTQ.ASCII_OFFSET;
				stream.FASTQ.ASCII_OFFSET=33;
				int[] qi=FileFormat.testInterleavedAndQuality(fname, false);
				offset=qi[0];
				interleaved=(qi[1]==FileFormat.INTERLEAVED);
				maxLen=Tools.max(maxLen, qi[2]);
				minLen=Tools.min(minLen, qi[2]);
				barcodesPerRead=qi[4];
				barcodeDelimiter=(byte)qi[5];
				stream.FASTQ.ASCII_OFFSET=qold;
			}else if(ff.fasta()){
				interleaved=stream.FASTQ.testInterleavedFasta(fname, false);
			}
		}
		
		return ff;
	}
	
	void processReads(FileFormat ff){
		
		if(initialQin>=0){
			FASTQ.ASCII_OFFSET=initialQin;
			FASTQ.DETECT_QUALITY=initialDetectQuality;
		}
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff, null);
			cris.start();
		}
		final boolean paired=cris.paired();
		barcodeStats=new BarcodeStats((char)barcodeDelimiter, barcodesPerRead, "Full");
		barcodeStats.length1=ff.barcodeLength(1);
		barcodeStats.length2=ff.barcodeLength(2);
		if(expectedCodeFile!=null){
			barcodeStats.loadBarcodeList(expectedCodeFile, barcodeDelimiter, false, false);
		}
		
		spawnThreads(ff, cris);
		
		ReadWrite.closeStream(cris);
		if(verbose){outstream.println("Finished.");}
	}
	
	void loadVars(FileFormat ff){
		final ByteFile bf=ByteFile.makeByteFile(ff);
		final byte delimiter='\t';
		byte[] line=bf.nextLine();
		while(line!=null && line.length>0){
			if(line[0]!='#'){
				variantsProcessed++;
//				Var v=new Var(line, delimiter);
			}else{
				headerLinesProcessed++;
				String[] split=new String(line).split("\t");
				String a=split[0], b=(split.length>1 ? split[1] : null);
				assert(split.length>1) : new String(line);
				if(a.equalsIgnoreCase("#ploidy")){
					ploidy=Integer.parseInt(b);
				}else if(a.equalsIgnoreCase("#pairingRate")){
					pairingRate=Double.parseDouble(b);
				}else if(a.equalsIgnoreCase("#totalQualityAvg")){
					totalQualityAvg=Double.parseDouble(b);
				}else if(a.equalsIgnoreCase("#mapqAvg")){
					mapqAvg=Double.parseDouble(b);
				}else if(a.equalsIgnoreCase("#readLengthAvg")){
					readLengthAvg=Double.parseDouble(b);
				}
			}
			line=bf.nextLine();
		}
		bf.close();
	}
	
	void loadVcf(FileFormat ff){
		ByteFile bf=ByteFile.makeByteFile(ff);
		byte[] line=bf.nextLine();
		while(line!=null && line.length>0){
			headerLinesProcessed++;
			if(line[0]!='#'){
				variantsProcessed++;
//				Var v;
//				try {
//					v = Var.fromVCF(line, null);
//				} catch (Exception e) {
//					System.err.println("Unable to parse VCF line: '"+new String(line)+"'");
//				}
			}else{
				String[] split=new String(line).split("=");
				if(split.length==2){
					String a=split[0], b=split[1];
					if(a.equalsIgnoreCase("##ploidy")){
						ploidy=Integer.parseInt(b);
					}else if(a.equalsIgnoreCase("##properPairRate")){
						pairingRate= Double.parseDouble(b);
					}else if(a.equalsIgnoreCase("##totalQualityAvg")){
						totalQualityAvg= Double.parseDouble(b);
					}else if(a.equalsIgnoreCase("##mapqAvg")){
						mapqAvg= Double.parseDouble(b);
					}else if(a.equalsIgnoreCase("##readLengthAvg")){
						readLengthAvg= Double.parseDouble(b);
					}
				}
			}
			line=bf.nextLine();
		}
		bf.close();
	}
	
	/*--------------------------------------------------------------*/
	
	private static final byte[] makeToNum(){
		byte[] array=new byte[128];
		Arrays.fill(array, (byte)8);
		array['a']=array['A']=0;
		array['c']=array['C']=1;
		array['g']=array['G']=2;
		array['t']=array['T']=3;
		array['u']=array['U']=4;
		array['n']=array['N']=5;
		array['-']=7;
		for(byte b : AminoAcid.degenerateBases){
			if(Character.isLetter(b)){
				array[b]=array[Tools.toLowerCase(b)]=6;
			}
		}
		return array;
	}
	
	private static final byte[] makeToAmino(){
		byte[] array=new byte[128];
		Arrays.fill(array, (byte)4);
		for(AminoAcid aa : AminoAcid.AlphabeticalAAs){
			array[aa.letter]=0;
			array[Tools.toLowerCase(aa.letter)]=0;
		}
		array['X']=array['x']=array['.']=1;
		array['*']=2;
		array['-']=3;
		return array;
	}
	
	private static final byte[] makeAminoOnly(){
		byte[] array=new byte[128];
		Arrays.fill(array, (byte)0);
		for(int i=0; i<128; i++){
			if(Character.isLetter(i) && AminoAcid.acidToNumberExtended[i]>=0 && AminoAcid.baseToNumberExtended[i]<0){array[i]=1;}
		}
		return array;
	}
	
	private static final byte[] makeToLUS(){
		byte[] array=new byte[128];
		for(int i=0; i<128; i++){
			if(Tools.isLowerCase(i)){
				array[i]=0;
			}else if(Tools.isUpperCase(i)){
				array[i]=1;
			}else{
				array[i]=2;
			}
		}
		return array;
	}
	
	private static byte toNum(byte b){
		return b<0 ? 7 : toNum[b];
	}
	
	private static byte toLUS(byte b){
		return b<0 ? 2 : toLUS[b];
	}
	
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private void spawnThreads(final FileFormat ff, final ConcurrentReadInputStream cris){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		int threads=Shared.threads();
		if(threads>=48){
			threads=Tools.min(40, threads/2);
		}
		
		//Fill a list with ProcessThreads
		ArrayList<TestThread> alpt=new ArrayList<TestThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new TestThread(ff, cris));
		}
		
		boolean success=true;
		if(threads<2){
			alpt.get(0).run();
		}else{

			//Start the threads
			for(TestThread pt : alpt){
				pt.start();
			}

			//Wait for completion of all threads
			for(TestThread pt : alpt){

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
			}
		}
		
		//Accumulate per-thread statistics
		for(TestThread pt : alpt){
			//Accumulate per-thread statistics
			readsProcessed+=pt.readsProcessed_T;
			pairsProcessed+=pt.pairsProcessed_T;
			basesProcessed+=pt.basesProcessed_T;
			mergeAttempts+=pt.mergeAttempts_T;
			success&=pt.success_T;
			
			acidsNotBasesProcessed+=pt.acidsNotBasesProcessed_T;
			junkProcessed+=pt.junkProcessed_T;
			chastityFail+=pt.chastityFail_T;
			badPairs+=pt.badPairs_T;
			adapterBases+=pt.adapterBases_T;
			adapterReads+=pt.adapterReads_T;

			minLen=Tools.min(minLen, pt.minLen_T);
			maxLen=Tools.max(maxLen, pt.maxLen_T);
			sll.add(pt.sllT);
			if(parseZmw && pt.parseZmwT){
				zmwMap.incrementAll(pt.zmwMapT);
			}else{
				parseZmw=false;
			}

			qMinUncalled=Tools.min(qMinUncalled, pt.qMinUncalledT);
			qMaxUncalled=Tools.max(qMaxUncalled, pt.qMaxUncalledT);
			qMinCalled=Tools.min(qMinCalled, pt.qMinCalledT);
			qMaxCalled=Tools.max(qMaxCalled, pt.qMaxCalledT);
			
			add(ACGTUNIGO, pt.ACGTUNIGO_T);
			add(AXEGO, pt.AXEGO_T);
			add(lowerUpperSymbol, pt.lowerUpperSymbol_T);
			add(qhist, pt.qhist_T);
			add(ihist, pt.ihist_T);
			add(trimhist, pt.trimhist_T);
			
			barcodeStats.merge(pt.barcodeStats_T);
			aTrack.merge(pt.aTrack_T);
			invalidHeaders.addAll(pt.invalidHeaders_T);
			
			if(makeSketch){
				smm.add(pt.smm_T);
			}
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}

		//Do anything necessary after processing
	}
	
	private void printToFile(long[] hist, boolean nzo, String header, String fname){
		printToFileOffset(hist, nzo, header, fname, 0);
	}
	
	private void printToFileOffset(long[] hist, boolean nzo, String header, String fname, int offset){
		ByteStreamWriter bsw=new ByteStreamWriter(fname, true, false, false);
		bsw.start();
		bsw.println(header.getBytes());
		for(int i=0; i<hist.length; i++){
			long x=hist[i];
			if(!nzo || x>0){
				bsw.print(i-offset);
				bsw.print('\t');
				bsw.print(x);
				bsw.print('\n');
			}
		}
		errorState|=bsw.poisonAndWait();
	}
	
	private static void add(long[] dest, long[] source){
		for(int i=0; i<source.length; i++){
			dest[i]+=source[i];
		}
	}
	
	private final class TestThread extends Thread {

		TestThread(FileFormat ff_, ConcurrentReadInputStream cris_){
			ff=ff_;
			cris=cris_;
			if(makeSketch){
				smm_T=new SketchMakerMini(tool, 0, SketchObject.defaultParams);
			}else{
				smm_T=null;
			}
		}
		
		@Override
		public void run(){
			processInThread();
			success_T=true;
		}
		
		void processInThread(){
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ff==null || ff.samOrBam()) || (r.mate!=null)==cris.paired());
				if(parseZmw) {
					try {
						int zmw=Parse.parseZmw(r.id);
						parseZmwT=(zmw>=0 && zmw<=1000000000);
					} catch (Exception e) {}
				}
				if(r.numericID==0){
					synchronized(cris) {
						assert(firstHeader==null);
						firstHeader=r.id;
					}
				}
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					processPair(r1, r2);
				}
				
				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		void processPair(Read r1, Read r2){
			assert(r1!=null);
			processRead(r1);
			processRead(r2);
			
			if(r2!=null){
				pairsProcessed_T++;
				String s1=r1==null ? null : r1.id;
				String s2=r2==null ? null : r2.id;
				boolean b=FASTQ.testPairNames(s1, s2, true);
				if(!b){badPairs_T++;}
				
				if(doMerge && r1.length()>10 && r1.mateLength()>10){
					mergeAttempts_T++;
					final int insert=BBMerge.findOverlapLoose(r1, r2, false);
					if(insert>0){
						ihist_T[Tools.min(insert, ihist_T.length-1)]++;
						int max=Tools.max(r1.length(), r2.length());
						if(insert<max){
							int trim1=(Tools.max(r1.length()-insert, 0));
							int trim2=(Tools.max(r2.length()-insert, 0));
							if(trim1>0){
								adapterBases_T+=trim1;
								adapterReads_T++;
							}
							if(trim2>0){
								adapterBases_T+=trim2;
								adapterReads_T++;
							}
							aTrack_T.storeAdapterSequence(r1, insert);
							aTrack_T.storeAdapterSequence(r2, insert);
						}
					}
				}
			}
		}
		
		void processRead(Read r){
			if(r==null){return;}
			final byte[] bases=r.bases;
			final byte[] quals=r.quality;
			final int len=r.length();
			readsProcessed_T++;
			basesProcessed_T+=len;
			maxLen_T=Tools.max(len, maxLen_T);
			minLen_T=Tools.min(len, minLen_T);
			sllT.add(len);
			
			if(parseZmwT && r.pairnum()==0){
				int zmw=-1;
				try {
					zmw=Parse.parseZmw(r.id);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(zmw>=0){
					zmwMapT.increment(zmw);
				}else{
					parseZmwT=false;
				}
			}
			
			if(makeBhist) {
				if(maxBhistLen<1 || r.length()<=maxBhistLen){
					readstatsT.addToBaseHistogram(r);
				}
			}
//			if(makeLhist){readstatsT.addToLengthHistogram(r);}
			if(makeGChist){readstatsT.addToGCHistogram(r);}
			
			boolean cf=r.failsChastity(false);
			if(r.junk() || cf){
				if(r.junk()){junkProcessed_T++;}
				if(cf){chastityFail_T++;}
				if(printJunk){invalidHeaders_T.add(r.id);}
			}
			
			if(r.pairnum()==0){addBarcode(r);}
			
			if(bases!=null){
				
				if(doTrim){
					testTrim(bases, quals);
				}
				
//				if(loglog!=null){loglog.hash(r);}
				if(smm_T!=null){smm_T.processRead(r);}
//				if(fast){
//					for(byte b : bases){
//						ACGTUNIGO[toNum[b]]++;
//						lowerUpperSymbol[toLUS[b]]++;
//						AXEGO[toAmino[b]]++;
//						acidsNotBasesProcessed+=aminoOnly[b];
//					}
//				}else{
					for(byte b : bases){
						if(b>=0){
							ACGTUNIGO_T[toNum[b]]++;
							lowerUpperSymbol_T[toLUS[b]]++;
							AXEGO_T[toAmino[b]]++;
							acidsNotBasesProcessed_T+=aminoOnly[b];
						}else{
							ACGTUNIGO_T[8]++;
							lowerUpperSymbol_T[2]++;
							AXEGO_T[4]++;
						}
					}
//				}
			}
			if(quals!=null){
				for(int i=0; i<quals.length; i++){
					byte q=quals[i];
					byte b=bases[i];
					qhist_T[q+128]++;
					
					if(AminoAcid.isFullyDefined(b)){
						qMinCalledT=Tools.min(q, qMinCalledT);
						qMaxCalledT=Tools.max(q, qMaxCalledT);
					}else{
						qMinUncalledT=Tools.min(q, qMinUncalledT);
						qMaxUncalledT=Tools.max(q, qMaxUncalledT);
					}
				}
			}
		}
		
		void testTrim(byte[] bases, byte[] quals){
			if(quals!=null){
				for(int trimq=5; trimq<=20; trimq++){
					long packed=TrimRead.testOptimal(bases, quals, QualityTools.PROB_ERROR[trimq]);
					int a0=(int)((packed>>32)&0xFFFFFFFFL);
					int b0=(int)((packed)&0xFFFFFFFFL);
					trimhist_T[trimq]+=a0+b0;
				}
			}else{
				long packed=TrimRead.testOptimal(bases, quals, QualityTools.PROB_ERROR[5]);
				int a0=(int)((packed>>32)&0xFFFFFFFFL);
				int b0=(int)((packed)&0xFFFFFFFFL);
				trimhist_T[10]+=a0+b0;
				trimhist_T[15]+=a0+b0;
				trimhist_T[20]+=a0+b0;
				trimhist_T[25]+=a0+b0;
				trimhist_T[30]+=a0+b0;
			}
		}
		
		void addBarcode(Read r){
			String code=r.barcode(false);
			if(code==null){return;}
			barcodeStats_T.increment(code, 1);
		}
		
		boolean parseZmwT=false;
		private boolean success_T=false;
		
		private final FileFormat ff;
		private final ConcurrentReadInputStream cris;

		private long readsProcessed_T=0;
		private long pairsProcessed_T=0;
		private long basesProcessed_T=0;
		private long mergeAttempts_T=0;
		private long acidsNotBasesProcessed_T=0;
		private long junkProcessed_T=0;
		private long chastityFail_T=0;
		private long badPairs_T=0;
		private long adapterBases_T=0;
		private long adapterReads_T=0;
		
		private int qMinUncalledT=999;
		private int qMaxUncalledT=-999;
		private int qMinCalledT=999;
		private int qMaxCalledT=-999;

		private long[] ACGTUNIGO_T=new long[9];
		private long[] AXEGO_T=new long[5];
		private long[] lowerUpperSymbol_T=new long[3];
		private long[] qhist_T=new long[256];
		private long[] ihist_T=new long[1000];
		private long[] trimhist_T=new long[51];
		private int minLen_T=Integer.MAX_VALUE;
		private int maxLen_T=0;
		
//		private HashMap<String, Barcode> barcodes_T=new HashMap<String, Barcode>();
		private BarcodeStats barcodeStats_T=barcodeStats.copySpecial();
		private AdapterTracker aTrack_T=new AdapterTracker();
		private ArrayList<String> invalidHeaders_T=new ArrayList<String>();
		private final SketchMakerMini smm_T;

		SuperLongList sllT=new SuperLongList(lengthLimit);
		private final ReadStats readstatsT=new ReadStats();

//		private final BitSet zmwsT=new BitSet(8192);
		private IntHashMapBinary zmwMapT=new IntHashMapBinary(2048);
		
		private final BandedAlignerConcrete bandy=(BarcodeStats.calcEdist ? new BandedAlignerConcrete(35) : null);
	}
	
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> in=new ArrayList<String>();
	
	/*--------------------------------------------------------------*/

	private long variantsProcessed=0;
	private long headerLinesProcessed=0;
	private long readsProcessed=0;
	private long pairsProcessed=0;
	private long basesProcessed=0;
	private long mergeAttempts=0;
	private long acidsNotBasesProcessed=0;
	private long junkProcessed=0;
	private long chastityFail=0;
	private long badPairs=0;
	private long adapterBases=0;
	private long adapterReads=0;

	private long[] ACGTUNIGO=new long[9];
	private long[] AXEGO=new long[5];
	private long[] lowerUpperSymbol=new long[3];
	private long[] qhist=new long[256];
	private long[] ihist=new long[1000];
	private long[] trimhist=new long[51];
	private int minLen=Integer.MAX_VALUE;
	private int maxLen=0;
	
	private int qMinUncalled=999;
	private int qMaxUncalled=-999;
	private int qMinCalled=999;
	private int qMaxCalled=-999;
	
	private final int lengthLimit=100000;
	private SuperLongList sll=new SuperLongList(lengthLimit);
//	private BitSet zmws=new BitSet(8192);
	private IntHashMapBinary zmwMap=new IntHashMapBinary(2048);
	private boolean parseZmw=true;
	
	private ArrayList<String> invalidHeaders=new ArrayList<String>();
	private final SketchTool tool;
	private final SketchMakerMini smm;
	private Sketch sketch=null;

	private final byte initialQin;
	private final boolean initialDetectQuality;

//	private HashMap<String, Barcode> barcodes=new HashMap<String, Barcode>();
	private int barcodesPerRead;
	private byte barcodeDelimiter;
	private BarcodeStats barcodeStats;
	
	private AdapterTracker aTrack=new AdapterTracker();
	
	/*--------------------------------------------------------------*/

	int ploidy=-1;
	double pairingRate=-1;
	double mapqAvg=-1;
	double totalQualityAvg=-1;
	double readLengthAvg=-1;
	
	/*--------------------------------------------------------------*/

	private int format=FileFormat.UNKNOWN;
	private int compression=FileFormat.UNKNOWN;
	private boolean amino=false;
	private boolean differs=false;
	private boolean interleaved=false;
	private int offset=33;
	private boolean makeSketch=true;
	private boolean doMerge=true;
	private boolean doTrim=true;
	private int sketchSize=40000;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
	private boolean full=true;
	private boolean fast=true;
	
	private boolean printSpeed=false;
	private boolean errorState=false;
	
	private boolean printBarcodes=false;
	private boolean printQhist=false;
	private boolean printIhist=false;
	private boolean printJunk=false;
	private boolean makeBhist;
//	private boolean forceMakeBhist=false;;
	private int maxBhistLen=10000;
	private final boolean makeLhist;
	private final boolean makeGChist;
	
	private String qhistFile="qhist.txt";
	private String ihistFile="ihist.txt";
	private String khistFile="khist.txt";
	private String bhistFile="bhist.txt";
	private String lhistFile="lhist.txt";
	private String gchistFile="gchist.txt";
	private String zmwhistFile="zmwhist.txt";
	
	private String adapterFile="adapters.fa";

	private String firstHeader=null;
	
	private String barcodeFile="barcode#.txt";
	private String junkFile="junk.txt";
	private String expectedCodeFile=null;//List of expected barcodes
	
	/*--------------------------------------------------------------*/

	private final static int qOffset=128;
	private static final byte[] toNum=makeToNum();
	private static final byte[] toLUS=makeToLUS();
	private static final byte[] toAmino=makeToAmino();
	private static final byte[] aminoOnly=makeAminoOnly();
	
	/*--------------------------------------------------------------*/
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
