package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;

import align2.QualityTools;
import dna.ChromosomeArray;
import dna.Data;
import dna.Scaffold;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.KillSwitch;
import shared.LineParser1;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.Read;
import stream.SamLine;
import stream.SamLineStreamer;
import stream.ScaffoldCoordinates;
import stream.SiteScore;
import structures.ByteBuilder;
import structures.CoverageArray;
import structures.CoverageArray2;
import structures.CoverageArray3;
import structures.ListNum;
import structures.LongList;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date Jan 4, 2013
 *
 */
public class CoveragePileup {
	
	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void main(String[] args){
		Timer t=new Timer();
		
		CoveragePileup x=new CoveragePileup(args);
		
		x.process();
		
		t.stop();
		if(!KEY_VALUE){
			x.outstream.println();
			x.outstream.println("Time: \t"+t);
		}
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public CoveragePileup(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, printCommand ? getClass() : null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		int vectorMode=-1;
		boolean outset=false;
		ReadWrite.USE_UNPIGZ=true;
//		SamLine.RNAME_AS_BYTES=false;
		
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("ref") || a.equals("reference") || a.equals("fasta")){
				reference=b;
			}else if(a.equals("addfromref")){
				ADD_FROM_REF=Parse.parseBoolean(b);
			}else if(a.equals("addfromreads")){
				ADD_FROM_READS=Parse.parseBoolean(b);
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("covin") || a.equals("incov")){
				covIn.clear();
				if(b!=null) {
					for(String cf : Tools.commaPattern.split(b)) {covIn.add(cf);}
				}
			}else if(a.equals("out") || a.equals("coveragestats") || a.equals("covstats") || a.equals("stats")){
				covstats=b;
				outset=true;
			}else if(a.equals("minscaf") || a.equals("covminscaf")){
				minscaf=Integer.parseInt(b);
			}else if(a.equals("mindepth") || a.equals("mincov")){
				minDepthToBeCovered=Integer.parseInt(b);
			}else if(a.equals("border")){
				border=Integer.parseInt(b);
			}else if(a.equals("qtrim")/* || a.equals("trim")*/){
				if(b==null || b.length()==0){qtrimRight=qtrimLeft=true;}
				else if(b.equalsIgnoreCase("left") || b.equalsIgnoreCase("l")){qtrimLeft=true;qtrimRight=false;}
				else if(b.equalsIgnoreCase("right") || b.equalsIgnoreCase("r")){qtrimLeft=false;qtrimRight=true;}
				else if(b.equalsIgnoreCase("both") || b.equalsIgnoreCase("rl") || b.equalsIgnoreCase("lr")){qtrimLeft=qtrimRight=true;}
				else if(Tools.isDigit(b.charAt(0))){
					trimq=Float.parseFloat(b);
					qtrimRight=trimq>0;
				}else{qtrimRight=qtrimLeft=Parse.parseBoolean(b);}
			}else if(a.equals("trimq") || a.equals("trimquality")){
				trimq=Float.parseFloat(b);
			}else if(a.equals("ss") || a.equals("samstreamer")){
				if(b!=null && Tools.isDigit(b.charAt(0))){
					useStreamer=true;
					streamerThreads=Tools.max(1, Integer.parseInt(b));
				}else{
					useStreamer=Parse.parseBoolean(b);
				}
			}else if(a.equals("minq") || a.equals("minmapq")){
				minMapq=Integer.parseInt(b);
			}else if(a.equals("outsam")){
				outsam=b;
			}else if(a.equals("rpkm") || a.equals("fpkm") || a.equals("outrpkm")){
				outrpkm=b;
			}else if(a.equals("outorf")){
				outorf=b;
			}else if(a.equals("orffasta") || a.equals("fastaorf")){
				orffasta=b;
			}else if(a.equals("basecov") || a.equals("outcov")){
				basecov=b;
			}else if(a.equals("ranges") || a.equals("rangecov")){
				rangecov=b;
//			}else if(a.equals("mincov") || a.equals("mindepth")){
//				mincov=Integer.parseInt(b);
			}else if(a.equals("bincov") || a.equals("outbinned")){
				bincov=b;
			}else if(a.equals("normcov") || a.equals("outnormalized")){
				normcov=b;
			}else if(a.equals("normcovo") || a.equals("outnormalizedoverall")){
				normcovOverall=b;
			}else if(a.equals("delta")){
				DELTA_ONLY=Parse.parseBoolean(b);
			}else if(a.equals("physical") || a.equals("physicalcoverage") || a.equals("physcov")){
				PHYSICAL_COVERAGE=Parse.parseBoolean(b);
			}else if(a.equals("tlen")){
				USE_TLEN=Parse.parseBoolean(b);
			}else if(a.equals("hist") || a.equals("histogram") || a.equals("covhist")){
				histogram=b;
			}else if(a.equals("reads")){
				maxReads=Parse.parseKMG(b);
			}else if(a.equals("scafs") || a.equals("scaffolds")){
				initialScaffolds=Tools.max(128, (int)(Tools.min(Long.parseLong(b),2000000000)));
			}else if(a.equals("binsize")){
				binsize=Integer.parseInt(b);
			}else if(a.equals("32bit")){
				bits32=Parse.parseBoolean(b);
			}else if(a.equals("bitset") || a.equals("usebitset") || a.equals("bitsets") || a.equals("usebitsets")){
//				if(Parse.parseBoolean(b)){arrayMode=BITSET_MODE;}
				vectorMode=Parse.parseBoolean(b) ? BITSET_MODE : NOTHING_MODE;
			}else if(a.equals("array") || a.equals("arrays") || a.equals("usearrays")){
				vectorMode=Parse.parseBoolean(b) ? ARRAY_MODE : NOTHING_MODE;
			}else if(a.equals("median") || a.equals("calcmedian")){
				if(Parse.parseBoolean(b)){
					vectorMode=ARRAY_MODE;
				}
			}else if(a.startsWith("nonzero") || a.equals("nzo")){
				NONZERO_ONLY=Parse.parseBoolean(b);
				if(verbose){outstream.println("Set NONZERO_ONLY to "+NONZERO_ONLY);}
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
				if(verbose){outstream.println("Set overwrite to "+overwrite);}
			}else if(a.equalsIgnoreCase("twocolumn")){
				TWOCOLUMN=Parse.parseBoolean(b);
				if(verbose){outstream.println("Set TWOCOLUMN to "+TWOCOLUMN);}
			}else if(a.equalsIgnoreCase("keyvalue") || a.equalsIgnoreCase("machineout")){
				KEY_VALUE=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("countgc")){
				COUNT_GC=Parse.parseBoolean(b);
				if(verbose){outstream.println("Set COUNT_GC to "+COUNT_GC);}
			}else if(a.equals("secondary") || a.equals("usesecondary")){
				USE_SECONDARY=Parse.parseBoolean(b);
				if(verbose){outstream.println("Set USE_SECONDARY_ALIGNMENTS to "+USE_SECONDARY);}
			}else if(a.equals("softclip") || a.equals("includesoftclip")){
				INCLUDE_SOFT_CLIP=Parse.parseBoolean(b);
				if(verbose){outstream.println("Set INCLUDE_SOFT_CLIP to "+INCLUDE_SOFT_CLIP);}
			}else if(a.equals("keepshortbins") || a.equals("ksb")){
				KEEP_SHORT_BINS=Parse.parseBoolean(b);
				if(verbose){outstream.println("Set KEEP_SHORT_BINS to "+KEEP_SHORT_BINS);}
			}else if(a.equals("strandedcoverage") || a.equals("strandedcov") || a.equals("covstranded") || a.equals("stranded")){
				STRANDED=Parse.parseBoolean(b);
			}else if(a.equals("startcov") || a.equals("covstart") || a.equals("startonly")){
				START_ONLY=Parse.parseBoolean(b);
			}else if(a.equals("stopcov") || a.equals("covstop") || a.equals("stoponly")){
				STOP_ONLY=Parse.parseBoolean(b);
			}else if(a.equals("concise")){
				CONCISE=Parse.parseBoolean(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("normc") || a.equals("normalizecoverage")){
				NORMALIZE_COVERAGE=Parse.parseBoolean(b);
			}else if(a.equals("header") || a.equals("hdr")){
				printHeader=Parse.parseBoolean(b);
			}else if(a.equals("headerpound") || a.equals("#")){
				headerPound=Parse.parseBoolean(b);
			}else if(a.equals("stdev")){
				calcCovStdev=Parse.parseBoolean(b);
			}else if(a.equals("delcov") || a.equals("dels") || a.equals("includedels") || a.equals("includedeletions") || a.equals("delcoverage")){
				INCLUDE_DELETIONS=Parse.parseBoolean(b);
			}else if(a.equals("dupecoverage") || a.equals("dupecov") || a.equals("dupes") || a.equals("duplicates") || a.equals("includeduplicates")){
				INCLUDE_DUPLICATES=Parse.parseBoolean(b);
			}else if(a.equals("ignoredupes") || a.equals("ignoreduplicates")){
				INCLUDE_DUPLICATES=!Parse.parseBoolean(b);
			}else if(a.equals("normb") || a.equals("normalizebins")){
				try {
					NORMALIZE_LENGTH_BINS=Integer.parseInt(b);
				} catch (NumberFormatException e) {
					boolean x=Parse.parseBoolean(b);
					NORMALIZE_LENGTH_BINS=x ? 100 : -1;
				}
			}else if(a.equals("covwindow")){
				if(b==null || b.length()<1 || Character.isLetter(b.charAt(0))){
					USE_WINDOW=Parse.parseBoolean(b);
				}else{
					LOW_COV_WINDOW=Integer.parseInt(b);
					USE_WINDOW=(LOW_COV_WINDOW>0);
				}
			}else if(a.equals("covwindowavg")){
				LOW_COV_DEPTH=Double.parseDouble(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(Parser.parseCommonStatic(arg, a, b)){
				//do nothing
			}else if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(in1==null && arg.indexOf('=')<0 && new File(arg).exists()){
				in1=arg;
			}else{
				throw new RuntimeException("Unknown parameter: "+args[i]);
			}
			
		}
		trimE=(float)QualityTools.phredToProbError(trimq);
//		assert(false) : qtrimLeft+", "+qtrimRight+", "+trimq+", "+trimE;
		
		ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_UNMAPPED;
		if(!USE_SECONDARY){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_SECONDARY;}
		if(!INCLUDE_DUPLICATES){ReadWrite.SAMTOOLS_IGNORE_FLAG|=ReadWrite.SAM_DUPLICATE;}
		
		if(outsam==null){
			SamLine.PARSE_0=false;
			SamLine.PARSE_6=false;
			SamLine.PARSE_7=false;
			SamLine.PARSE_8=false;
			if(k<1 && trimq<1){SamLine.PARSE_10=false;}
			SamLine.PARSE_OPTIONAL=false;
		}
		
		if(vectorMode>-1){
			USE_BITSETS=(vectorMode==BITSET_MODE);
			USE_COVERAGE_ARRAYS=(vectorMode==ARRAY_MODE);
		}else{
			if(histogram==null && basecov==null && bincov==null && normcov==null && rangecov==null &&
					normcovOverall==null && outorf==null && !calcCovStdev){//No need for coverage array!
				USE_COVERAGE_ARRAYS=false;
				if(TWOCOLUMN){//No need for bitset, either!
					USE_BITSETS=false;
				}else{
					USE_BITSETS=true;
				}
			}
		}
		
		if(verbose){
			outstream.println("Set USE_COVERAGE_ARRAYS to "+USE_COVERAGE_ARRAYS);
			outstream.println("Set USE_BITSETS to "+USE_BITSETS);
		}
		
		if(maxReads<0){maxReads=Long.MAX_VALUE;}
		{
			final String a=(args.length>0 ? args[0] : null);
//			final String b=(args.length>1 ? args[1] : null);
			if(in1==null && a!=null && a.indexOf('=')<0 && (a.startsWith("stdin") || new File(a).exists())){in1=a;}
//			if(covstats==null && b!=null && b.indexOf('=')<0){covstats=b;}
			if(in1==null){in1="stdin";}
//			if(covstats==null && !outset){
////				out="stdout";
////				outstream.println("Warning: output destination not set; producing no output.  To print to standard out, set 'out=stdout'");
//				outstream=System.err;
//			}
		}
		assert(in1!=null);
//		assert(out!=null || outset) : "Output file was not set.";
		
		if(STRANDED){
			assert(basecov==null || basecov.indexOf('#')>=0) : "Output filenames must contain '#' symbol for strand-specific output.";
			assert(bincov==null || bincov.indexOf('#')>=0) : "Output filenames must contain '#' symbol for strand-specific output.";
			assert(normcov==null || normcov.indexOf('#')>=0) : "Output filenames must contain '#' symbol for strand-specific output.";
			assert(normcovOverall==null || normcovOverall.indexOf('#')>=0) : "Output filenames must contain '#' symbol for strand-specific output.";
			assert(histogram==null || histogram.indexOf('#')>=0) : "Output filenames must contain '#' symbol for strand-specific output.";
			assert(covstats==null || covstats.indexOf('#')>=0) : "Output filenames must contain '#' symbol for strand-specific output.";
		}
		
		if(!Tools.testInputFiles(false, true, in1, in2, reference)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		if(!Tools.testOutputFiles(overwrite, append, false, basecov, rangecov, bincov, normcov, normcovOverall, histogram, covstats, outrpkm)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+
					basecov+", "+rangecov+", "+bincov+", "+normcov+", "+normcovOverall+", "+histogram+", "+covstats+", "+outrpkm+"\n");
		}
	}

	
	/*--------------------------------------------------------------*/
	/*----------------       Data Structures        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/** The goal if this is to garbage-collect unnecessary objects, not really for reusing the object */
	public void clear(){
		list=null;
		table=null;
		pairTable=null;
		
		program=null;
		version=null;
		
		in1=null;
		covstats=null;
		outsam=null;
		outorf=null;
		outrpkm=null;
		reference=null;
		histogram=null;
		basecov=null;
		bincov=null;
		normcov=null;
		normcovOverall=null;
		orffasta=null;
		
		error=false;

		refBases=0;
		mappedBases=0;
		mappedNonClippedBases=0;
		mappedBasesWithDels=0;
		mappedReads=0;
		properPairs=0;
		readsProcessed=0;
		basesProcessed=0;
		kmersProcessed=0;
		mappedKmers=0;
		totalCoveredBases1=0;
		totalCoveredBases2=0;
		scaffoldsWithCoverage1=0;
		scaffoldsWithCoverage2=0;
		totalScaffolds=0;
	}
	
	public void createDataStructures(){
		refBases=0;
		mappedBases=0;
		mappedNonClippedBases=0;
		mappedBasesWithDels=0;
		mappedReads=0;
		properPairs=0;
		readsProcessed=0;
		basesProcessed=0;
		kmersProcessed=0;
		mappedKmers=0;
		totalCoveredBases1=0;
		totalCoveredBases2=0;
		scaffoldsWithCoverage1=0;
		scaffoldsWithCoverage2=0;
		totalScaffolds=0;
		error=false;
		list=new ArrayList<Scaffold>(initialScaffolds);
		table=new HashMap<String, Scaffold>(initialScaffolds);
		
		if(PHYSICAL_COVERAGE){
			pairTable=new HashMap<String, SamLine>();
			if(COUNT_GC){
				COUNT_GC=false;
				outstream.println("COUNT_GC disabled for physical coverage mode.");
			}
			if(USE_SECONDARY){
				USE_SECONDARY=false;
				outstream.println("USE_SECONDARY disabled for physical coverage mode.");
			}
			
			SamLine.PARSE_0=true;
			SamLine.PARSE_6=true;
			SamLine.PARSE_7=true;
			SamLine.PARSE_8=true;
			SamLine.PARSE_10=false;
			SamLine.PARSE_OPTIONAL=false;
		}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Read and process all input data. */
	public void process(){
		createDataStructures();

//		System.err.println("A");
		if(in1!=null && FileFormat.isBamFile(in1)){
			if(Data.SAMTOOLS() && useStreamer){ReadWrite.USE_SAMBAMBA=false;} //Disable because it takes forever to read the header
		}
		ByteFile tf=ByteFile.makeByteFile(in1, false);

//		System.err.println("B");
		
		final ByteStreamWriter tsw=(outsam==null ? null : new ByteStreamWriter(outsam, overwrite, false, true));
		if(tsw!=null){tsw.start();}
		ReadWrite.USE_SAMBAMBA=true;

//		System.err.println("C");
		processHeader(tf, tsw);

//		System.err.println("D");
		processReference();
//		System.err.println("E");
		if(maxReads<0){maxReads=Long.MAX_VALUE;}
		
		if(covIn!=null && !covIn.isEmpty()) {
			for(String s : covIn) {
				loadCoverageFile(s);
			}
		}else if(useStreamer && !FileFormat.isStdio(in1)){
			tf.close();
//			System.err.println("F");
			processViaStreamer(tsw);
//			System.err.println("G");
		}else{
//			processViaByteFile(tf, line, tsw);
			processViaByteFile(tf, tsw);
//			System.err.println("H");
		}
//		System.err.println("I");
		
		printOutput();
		
		if(orffasta!=null){
			processOrfsFasta(orffasta, outorf, table);
		}
		
		if(tsw!=null){tsw.waitForFinish();}
	}
	
	private void processViaStreamer(ByteStreamWriter tsw){
		SamLineStreamer ss=new SamLineStreamer(in1, streamerThreads, false, maxReads);
		ss.start();
		ByteBuilder bb=new ByteBuilder(33000);
		for(ListNum<SamLine> ln=ss.nextLines(); ln!=null && ln.size()>0; ln=ss.nextLines()){
			ArrayList<SamLine> list=(ln==null ? null : ln.list);
			for(SamLine sl : list){
				if(tsw!=null){
					sl.toBytes(bb);
					bb.nl();
					if(bb.length>=16384){
						tsw.print(bb);
						bb.clear();
					}
				}
				processSamLine(sl);
			}
		}
		if(tsw!=null){
			if(bb.length>0){tsw.print(bb);}
			tsw.poison();
		}
	}
	
	private void processViaByteFile(ByteFile tf, ByteStreamWriter tsw){
		for(byte[] line=tf.nextLine(); line!=null && readsProcessed<maxReads; line=tf.nextLine()){
			if(tsw!=null){tsw.println(line);}
			processSamLine(line);
		}
		
		tf.close();
		if(tsw!=null){tsw.poison();}
	}
	
	//Untested
	//See also CoverageArray.loadDepth; that implementation is better (also untested)
	private void loadCoverageFile(String fname) {
		ByteFile bf=ByteFile.makeByteFile(fname, true);
		
		byte[] line=bf.nextLine();
		for(; line!=null && Tools.startsWith(line, '#'); line=bf.nextLine()) {} //Skip header
		if(line==null) {return;}
		lp.set(line);
		String name=lp.parseString(0);
		Scaffold scaf=table.get(name);
		
		for(; line!=null; line=bf.nextLine()) {
			assert(!Tools.startsWith(line, '#'));
			lp.set(line);
			int len=lp.length(0);
			if(len>0 && !lp.termEquals(name, 0)) {
				name=lp.parseString(0);
				scaf=table.get(name);
			}
			int pos=lp.parseInt(1);
			int cov=lp.parseInt(2);
			if(cov>0) {
				if(USE_COVERAGE_ARRAYS){
					if(scaf.obj0==null){
						scaf.obj0=(bits32 ? new CoverageArray3(table.size(), scaf.length) : new CoverageArray2(table.size(), scaf.length));
					}
					CoverageArray ca=(CoverageArray)scaf.obj0;
					ca.increment(pos, cov);
				}else if(USE_BITSETS){
					if(scaf.obj0==null){scaf.obj0=new BitSet(scaf.length+1);}
					BitSet bs=(BitSet)scaf.obj0;
					bs.set(pos);
				}
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Setup            ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/** Process all sam header lines from the tf.
	 * Once a non-header line is encountered, return it.
	 * If non-null, print all lines to the tsw. */
	public void processHeader(ByteFile tf, ByteStreamWriter tsw){
		byte[] line=null;
		for(line=tf.nextLine(); line!=null && (line.length==0 || line[0]=='@'); line=tf.nextLine()){
//			System.err.println(new String(line));
			if(tsw!=null){tsw.println(line);}

			if(line.length>2){
				final byte a=line[1], b=line[2];

				if(a=='S' && b=='Q'){
					lp.set(line);
					Scaffold scaf=new Scaffold(lp);
					if(COUNT_GC){scaf.basecount=KillSwitch.allocLong1D(8);}
					assert(!table.containsKey(scaf.name)) : "\nDuplicate scaffold name!\n"+scaf+"\n\n"+table.get(scaf.name);
					table.put(scaf.name, scaf);
					list.add(scaf);
					refBases+=scaf.length;
//					sc.obj=new CoverageArray2(table.size(), sc.length+1);
//					outstream.println("Made scaffold "+sc.name+" of length "+sc.length);
				}else if(a=='P' && b=='G'){
					lp.set(line);
					for(int i=1, terms=lp.terms(); i<terms; i++) {
						if(lp.termStartsWith("PN:", i)){
							if(program==null){
								lp.incrementA(3);
								program=lp.parseStringFromCurrentField();
							}
						}else if(lp.termStartsWith("VN:", i)){
							if(version==null){
								lp.incrementA(3);
								version=lp.parseStringFromCurrentField();
							}
						}
					}
				}else if(a=='R' && b=='G'){
					//Do nothing
				}else if(a=='H' && b=='D'){
					//Do nothing
				}else if(a=='C' && b=='O'){
					//Do nothing
				}else{
					//				assert(false) : line;
				}
			}
		}
		if(line!=null){tf.pushBack(line);}
//		return line;
	}
	
	
	public void loadScaffoldsFromIndex(int minChrom, int maxChrom){
		
		final int[][] lengths=Data.scaffoldLengths;
		final int[][] locs=Data.scaffoldLocs;
		final byte[][][] names=Data.scaffoldNames;
		final int[] counts=new int[8];
		for(int chrom=minChrom; chrom<=maxChrom; chrom++){
			final ChromosomeArray ca=Data.getChromosome(chrom);
//			assert(false) : lengths[chrom]+", "+lengths.length+", "+names[chrom]+", "+locs[chrom]+", "+(ca==null);
			if(lengths[chrom]!=null){
				final int[] clengths=lengths[chrom];
				final int[] clocs=locs[chrom];
				final byte[][] cnames=names[chrom];
				for(int idx=0; idx<clengths.length; idx++){
					final int length=clengths[idx];
					final int loc=clocs[idx];
					final String name=new String(cnames[idx]);
					final Scaffold scaf=new Scaffold(name, length);
					if(ca!=null){
						scaf.gc=ca.calcGC(loc, length, counts);
					}
					if(COUNT_GC){scaf.basecount=KillSwitch.allocLong1D(8);}
					assert(!table.containsKey(scaf.name)) : "\nDuplicate scaffold name!\n"+scaf+"\n\n"+table.get(scaf.name);
					table.put(scaf.name, scaf);
					list.add(scaf);
					refBases+=scaf.length;
				}
			}
		}
	}
	
	
	public void processReference(){
		if(reference==null){return;}

		ByteFile bf=ByteFile.makeByteFile(reference, false);
		Scaffold scaf=null;
		int len=0;
		final long[] acgtn=KillSwitch.allocLong1D(8);
		boolean addLen=false;
		for(byte[] s=bf.nextLine(); s!=null; s=bf.nextLine()){
			if(s.length>0 && s[0]=='>'){
				if(scaf!=null){
					if(scaf.length>0 && scaf.length!=len){
						outstream.println("ERROR: Scaffold "+scaf.name+" has contradictory lengths of "+scaf.length+" and "+len+"\n"
								+ "This probably indicates a corrupt or incorrect reference.");
						errorState=true;
						if(Shared.EA()){KillSwitch.kill();}
					}
//					else{outstream.println(scaf.name+", "+scaf.length+", "+len);}
					scaf.length=len;
					if(addLen){
						refBases+=scaf.length;
						addLen=false;
					}
					scaf.gc=(float)((acgtn[1]+acgtn[2])*1d/Data.max(1, acgtn[0]+acgtn[1]+acgtn[2]+acgtn[3]));
					scaf=null;
					len=0;
					Arrays.fill(acgtn, 0);
				}
				
				String name=new String(s, 1, s.length-1);
				scaf=table.get(name);
				if(ADD_FROM_REF && scaf==null){
					scaf=new Scaffold(name, 0);
					if(!warned){
						outstream.println("Warning - SAM header did not include "+name+"\nAbsent scaffolds will be added; further warnings will be suppressed.");
						warned=true;
					}
					if(COUNT_GC){scaf.basecount=KillSwitch.allocLong1D(8);}
					table.put(name, scaf);
					list.add(scaf);
					addLen=true;
				}
			}else{
				len+=s.length;
				for(int i=0; i<s.length; i++){
					acgtn[charToNum[s[i]]]++;
				}
			}
		}
		if(scaf!=null){
			if(scaf.length>0 && scaf.length!=len){
				outstream.println("ERROR: Scaffold "+scaf.name+" has contradictory lengths of "+scaf.length+" and "+len+"\n"
						+ "This probably indicates a corrupt or incorrect reference.");
				errorState=true;
				if(Shared.EA()){KillSwitch.kill();}
			}
//			else{outstream.println(scaf.name+", "+scaf.length+", "+len);}
			scaf.length=len;
			if(addLen){
				refBases+=scaf.length;
				addLen=false;
			}
			scaf.gc=(float)((acgtn[1]+acgtn[2])*1d/Data.max(1, acgtn[0]+acgtn[1]+acgtn[2]+acgtn[3]));
			scaf=null;
			len=0;
			Arrays.fill(acgtn, 0);
		}
	}
	
	
	public void processOrfsFasta(String fname_in, String fname_out, HashMap<String, Scaffold> map){
		TextFile tf=new TextFile(fname_in, false);
		assert(!fname_in.equalsIgnoreCase(fname_out));
		TextStreamWriter tsw=new TextStreamWriter(fname_out, overwrite, false, true);
		tsw.start();
		
		if(printHeader){
			String pound=(headerPound ? "#" : "");
			tsw.print(pound+"mappedBases="+mappedBases+"\n");
			tsw.print(pound+"mappedNonClippedBases="+mappedNonClippedBases+"\n");
			tsw.print(pound+"mappedBasesWithDels="+mappedBasesWithDels+"\n");
			tsw.print(pound+"mappedReads="+mappedReads+"\n");
			tsw.print(pound+"name\tlength\tdepthSum\tavgDepth\tavgDepth/mappedBases\tminDepth\tmaxDepth\tmedianDepth\tstdDevDepth\tfractionCovered\n");
		}
		
		String line;
		final StringBuilder sb=new StringBuilder(500);
//		Formatter formatter=new Formatter(sb);
		
		while((line=tf.nextLine())!=null){
			if(line.length()>1 && line.charAt(0)=='>'){
				
				String[] split=line.split(" # "); //' # ' used as delimiters
				
				String orfname=split[0].substring(1).trim(); //In case there are spaces around the ' # ' delimiters
				String scafname=orfname;
				if(scafname.contains("_")){//PRODIGAL pads _1 to the name of the first orf of a scaffold, and etc
					int last=scafname.lastIndexOf('_');
					boolean numeric=false;
					for(int i=last+1; i<scafname.length(); i++){
						if(Tools.isDigit(scafname.charAt(i))){numeric=true;}
						else{numeric=false; break;}
					}
					if(numeric){scafname=scafname.substring(0, last);}
				}
				
				int start=Integer.parseInt(split[1].trim());
				int stop=Integer.parseInt(split[2].trim());
				int strand=Integer.parseInt(split[3].trim());
				if(strand==1){strand=Shared.PLUS;}else{strand=Shared.MINUS;}
				Orf orf=new Orf(orfname, start, stop, (byte)strand);
				
				Scaffold scaf=map.get(scafname);
//				if(scaf==null){scaf=map.get(orfname);}
				
//				assert(scaf!=null) : "\nCan't find scaffold for ("+orf+")\nfrom line\n"+line+"\n";
//				assert(orf.start>=0 && orf.stop<scaf.length) : "\norf goes out of scaffold bounds.\n"+orf+"\n"+scaf+"\n";
				
				if(scaf==null){
					outstream.println("Can't find scaffold for ("+orf+")\nfrom line\n"+line+"\nscafname='"+scafname+"'\norfname='"+orfname+"'");
					if(ABORT_ON_ERROR){
						tsw.poison();
						throw new RuntimeException("Aborting.");
					}
				}else{
					if(orf.start<0 && orf.stop>=scaf.length){
						outstream.println("orf goes out of scaffold bounds.\n"+orf+"\n"+scaf);
						if(ABORT_ON_ERROR){
							tsw.poison();
							throw new RuntimeException("Aborting.");
						}
					}
					
					CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0); //TODO:  Strand logic here depends on stranding protocol.
					orf.readCoverageArray(ca);
				}
				
//				{
//					tsw.print(Tools.format("%s\t%d\t", args));
//				}
				
				sb.append(orf.name).append('\t');
				sb.append(orf.length()).append('\t');
				sb.append(orf.baseDepth).append('\t');
				sb.append(Tools.format("%.4f", orf.avgCoverage())).append('\t');
				sb.append(orf.avgCoverage()/mappedNonClippedBases);

				sb.append('\t');
				sb.append(orf.minDepth).append('\t');
				sb.append(orf.maxDepth).append('\t');
				sb.append(orf.medianDepth).append('\t');
				sb.append(Tools.format("%.4f",orf.stdevDepth)).append('\t');
				sb.append(Tools.format("%.4f",orf.fractionCovered()));

				sb.append('\n');
				tsw.print(sb.toString());
				sb.setLength(0);
			}
		}
		
		tsw.poisonAndWait();
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean addCoverage(final String scafName, final byte[] seq, byte[] match, final int start0, final int stop0, final int readlen, 
			final int nonClippedBases, final int strand, int incrementFrags, boolean properPair, SamLine sl){//sl is optional
		Scaffold scaf=table.get(scafName);
		if(scaf==null){
			if(ADD_FROM_READS){
				if(!warned){
					outstream.println("Warning: A read was mapped to unknown reference sequence "+scafName+
							"\nFurther warnings will be suppressed.  These scaffolds will be included, but \n"
							+ "some coverage statistics will be inaccurate as scaffold lengths are not known.");
					warned=true;
				}
				scaf=new Scaffold(scafName, 0);
				if(COUNT_GC){scaf.basecount=KillSwitch.allocLong1D(8);}
				table.put(scafName, scaf);
				list.add(scaf);
				return addCoverage(scaf, seq, match, start0, stop0, readlen, nonClippedBases, strand, incrementFrags, properPair, sl);
			}else if(EA){
				KillSwitch.kill("ERROR: A read was mapped to unknown reference sequence "+scafName);
			}else if(!warned){
				outstream.println("Warning: Can't find "+scafName+"\nThis scaffold will not be included.  Further warnings will be suppressed.");
				warned=true;
				error=true;
			}
			return false;
		}
		return addCoverage(scaf, seq, match, start0, stop0, readlen, nonClippedBases, strand, incrementFrags, properPair, sl);
	}
	
	public boolean addCoverage(final Scaffold scaf, final byte[] seq, byte match[], final int start0, final int stop0, final int readlen, final int nonClippedBases, 
			final int strand, int incrementFrags, boolean properPair, SamLine sl){//sl is optional
		if(scaf==null){
			assert(false) : "Adding coverage to a null Scaffold.";
			return false;
		}
		
		if(ADD_FROM_READS){
			int x=scaf.length;
			scaf.length=Tools.max(scaf.length, stop0+1);
			if(scaf.length!=x){
				refBases=refBases+scaf.length-x;
			}
		}
		
		final int start=Tools.max(start0, 0);
		final int stop=Tools.min(stop0, scaf.length-1);
		
		assert(start>=0 && stop>=0) : "\nAn error was encountered when processing a read. Output will not be valid.\n"+
			"\nscafName="+scaf.name+"\nseq="+new String(seq)+"\nstart="+start+
			"\nstop="+stop+"\nreadlen="+readlen+"\nstrand="+strand+"\nscaf.length="+scaf.length+"\nscaf="+scaf;
		
		mappedBases+=readlen;
		mappedNonClippedBases+=nonClippedBases;
		mappedBasesWithDels+=(stop-start+1);
		mappedReads++;
		if(properPair){properPairs++;}

		scaf.readhits++;
		scaf.fraghits+=incrementFrags;
		if(strand==1){scaf.readhitsMinus++;}
		
		if(seq!=null && scaf.basecount!=null){
			final long[] counts=scaf.basecount;
			for(int i=0; i<seq.length; i++){
				counts[charToNum[seq[i]]]++;
			}
		}
		
		if(!INCLUDE_DELETIONS && !START_ONLY && !STOP_ONLY){
			assert(match!=null) : "Coverage excluding deletions cannot be calculated without a match string.";
			return addCoverageIgnoringDeletions(scaf, seq, match, start, stop, readlen, strand, incrementFrags);
		}
		
		final int basehits=stop-start+1;
		scaf.basehits+=basehits;
		
		if(USE_COVERAGE_ARRAYS){
			if(scaf.obj0==null){
				scaf.obj0=(bits32 ? new CoverageArray3(table.size(), scaf.length) : new CoverageArray2(table.size(), scaf.length));
				if(STRANDED){
					scaf.obj1=(bits32 ? new CoverageArray3(table.size(), scaf.length) : new CoverageArray2(table.size(), scaf.length));
				}
			}
			CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
			if(START_ONLY){
				ca.increment(start);
			}else if(STOP_ONLY){
				ca.increment(stop);
			}else{
				ca.incrementRange(start, stop);
//				mappedBases+=(stop-start+1);
//				mappedBases+=(readlen);
//				assert(readlen==(stop-start+1)) : "readlen="+readlen+", (stop-start+1)="+(stop-start+1)+", start="+start+", stop="+stop+", start0="+start0+", stop0="+stop0+
//				"\n"+sl+"\n";
				//123
			}
		}else if(USE_BITSETS){
			if(scaf.obj0==null){
				scaf.obj0=new BitSet(scaf.length+1);
				if(STRANDED){
					scaf.obj1=new BitSet(scaf.length+1);
				}
			}
			BitSet bs=(BitSet)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
			if(START_ONLY){
				bs.set(start);
			}else if(STOP_ONLY){
				bs.set(stop);
			}else{
				bs.set(start, stop+1);
			}
		}
		
		return true;
	}
	
	private boolean addCoverageIgnoringDeletions(final Scaffold scaf, final byte[] seq, byte match[], final int start, final int stop, final int readlen, final int strand, int incrementFrags){
		assert(!INCLUDE_DELETIONS && !START_ONLY && !STOP_ONLY);
		assert(match!=null) : "Coverage excluding deletions cannot be calculated without a match string.";
		
		if(Read.isShortMatchString(match)){
			match=Read.toLongMatchString(match);
		}
		
		int basehits=0;
		if(USE_COVERAGE_ARRAYS){
			if(scaf.obj0==null){
				scaf.obj0=(bits32 ? new CoverageArray3(table.size(), scaf.length) : new CoverageArray2(table.size(), scaf.length));
				if(STRANDED){
					scaf.obj1=(bits32 ? new CoverageArray3(table.size(), scaf.length) : new CoverageArray2(table.size(), scaf.length));
				}
			}
			CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
			for(int rpos=start, mpos=0; mpos<match.length && rpos<=stop; mpos++){
				byte m=match[mpos];
				if(m=='m' || m=='S' || m=='N'){
					ca.increment(rpos, 1);
					basehits++;
					rpos++;
				}else if(m=='X' || m=='Y' || m=='C' || m=='I'){
					//do nothing
				}else if(m=='D'){
					rpos++;
				}else{
					assert(false) : "Unhandled symbol "+m;
				}
			}
		}else if(USE_BITSETS){
			if(scaf.obj0==null){
				scaf.obj0=new BitSet(scaf.length+1);
				if(STRANDED){
					scaf.obj1=new BitSet(scaf.length+1);
				}
			}
			BitSet bs=(BitSet)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
			for(int rpos=start, mpos=0; mpos<match.length && rpos<=stop; mpos++){
				byte m=match[mpos];
				if(m=='m' || m=='S' || m=='N'){
					bs.set(rpos);
					basehits++;
					rpos++;
				}else if(m=='X' || m=='Y' || m=='C' || m=='I'){
					//do nothing
				}else if(m=='D'){
					rpos++;
				}else{
					assert(false) : "Unhandled symbol "+m;
				}
			}
		}
		scaf.basehits+=basehits;
		
		return true;
	}
	
	
	public boolean processSamLine(byte[] line){
		if(line==null || line.length<3){
			return false;
		}else if(line[0]=='@'){
			if(!error){
				outstream.println("Unexpected header line: "+line);
				outstream.println("This should not cause problems, and is probably due to concatenated sam files.\n" +
						"Supressing future unexpected header warnings.");
				error=true;
			}
			
			if(line[1]=='S' && line[2]=='Q'){
				lp.set(line);
				String name=Scaffold.name(lp);
				if(!table.containsKey(name)){
					Scaffold scaf=new Scaffold(lp);
					if(COUNT_GC){scaf.basecount=KillSwitch.allocLong1D(8);}
					table.put(scaf.name, scaf);
					list.add(scaf);
					refBases+=scaf.length;
				}
			}
		}else{
			SamLine sl=new SamLine(line);
			return processSamLine(sl);
		}
		return false;
	}
	
	private int trim(SamLine sl){
		assert(border>0 || (trimq>=0 && (qtrimLeft || qtrimRight)));
		Read r=null;
		Scaffold sc=null;
		int leftTrimAmount=border, rightTrimAmount=border;
		if(border>0){
			r=sl.toRead(false);
			sc=table.get(sl.rnameS());
			assert(sc!=null) : sl+"\n\n"+sl.rnameS()+"\n\n"+table;
			int skipTrimRange=Tools.max(10, border+5);
			if(r.start<skipTrimRange){
				if(r.strand()==Shared.PLUS){leftTrimAmount=0;}
				else{rightTrimAmount=0;}
			}
			if(r.stop>sc.length-skipTrimRange){
				if(r.strand()==Shared.PLUS){rightTrimAmount=0;}
				else{leftTrimAmount=0;}
			}
		}
		final int len0=sl.length();
		if(qtrimLeft || qtrimRight){
			long packed=TrimRead.testOptimal(sl.seq, sl.qual, trimE);
			if(qtrimLeft){leftTrimAmount=Tools.max(leftTrimAmount, (int)((packed>>32)&0xFFFFFFFFL));}
			if(qtrimRight){rightTrimAmount=Tools.max(rightTrimAmount, (int)((packed)&0xFFFFFFFFL));}
//			assert(false) : qtrimLeft+", "+qtrimRight+", "+trimq+", "+trimE+", "+rightTrimAmount+"\n"+sl+"\n";
		}
		final int trimmed;
		if(leftTrimAmount<1 && rightTrimAmount<1){trimmed=0;}
		else{
			if(r==null){r=sl.toRead(false);}
			if(sc==null){sc=table.get(sl.rnameS());}
			int scaflen=(sc==null ? 1999999999 : sc.length);
			trimmed=TrimRead.trimReadWithMatch(r, sl, leftTrimAmount, rightTrimAmount, 0, scaflen, false);
		}
//		assert(trimmed==len0-sl.length()) : trimmed+", "+len0+", "+sl.length();
//		assert(rightTrimAmount>0) : qtrimLeft+", "+qtrimRight+", "+trimq+", "+trimE+", "+rightTrimAmount+"\n"+sl+"\n";
		return trimmed;
	}
	
	public boolean processSamLine(SamLine sl){
		readsProcessed++;
		basesProcessed+=sl.length();
		if(sl.duplicate() && !INCLUDE_DUPLICATES){return false;}
		
		if(border>0 || (trimq>=0 && (qtrimLeft || qtrimRight))){
			if(sl.mapped() && sl.seq!=null){
				int trimmed=trim(sl);
				if(sl.length()<1 || trimmed<0){return false;}//trimmed<0 implies everything was trimmed
			}
		}
		
		final int kmers=Tools.countKmers(sl.seq, k);
		final double cKmers=sl.qual==null ? kmers : Tools.countCorrectKmers(sl.qual, k);
		kmersProcessed+=kmers;
		correctKmers+=cKmers;
		final boolean properPair=(sl.hasMate() && sl.mapped() && sl.primary() && sl.properPair());
		if(PHYSICAL_COVERAGE && properPair){
			SamLine mate=pairTable.remove(sl.qname);
			if(mate==null){pairTable.put(sl.qname, sl);}
			else{
				final int start1=sl.start(INCLUDE_SOFT_CLIP, false);
				final int stop1=sl.stop(start1, INCLUDE_SOFT_CLIP, false);
				final int start2=mate.start(INCLUDE_SOFT_CLIP, false);
				final int stop2=mate.stop(start2, INCLUDE_SOFT_CLIP, false);
				final int strand=(sl.pairnum()==0 ? sl.strand() : mate.strand());
				final int length=USE_TLEN ? sl.tlen : Tools.max(stop1, stop2)-Tools.min(start1, start2)+1;
				mappedKmers+=kmers;
				addCoverage(sl.rnameS(), null, null, Tools.min(start1, start2), Tools.max(stop1, stop2), length, sl.mappedNonClippedBases(), strand, 2, sl.properPair(), sl);
			}
		}else if(sl.mapped() && (USE_SECONDARY || sl.primary()) && sl.mapq>=minMapq){
			assert(sl.seq!=null || sl.cigar!=null) : "This program requires bases or a cigar string for every sam line.  Problem line:\n"+sl+"\n";
//			assert(sl.seq!=null) : sl.toString();
			final int length=sl.length();
			final int start=sl.start(INCLUDE_SOFT_CLIP, false);
			final int stop=sl.stop(start, INCLUDE_SOFT_CLIP, false);
//			assert(false && length==stop-start+1) : length+", "+start+", "+stop+", "+(stop-start+1);
//			assert(false) : "'"+new String(sl.rname())+"', '"+sl.rnameS()+"'";
//			assert(false) : "'"+sl.rnameS()+"'";
			final byte[] match=(INCLUDE_DELETIONS ? null : sl.toShortMatch(true));
			mappedKmers+=kmers;
			return addCoverage(sl.rnameS(), sl.seq, match, start, stop, length, sl.mappedNonClippedBases(), sl.strand(), sl.hasMate() ? 1 : 2, sl.properPair(), sl);
		}
		return false;
	}
	
	
	public boolean processRead(Read r){
		readsProcessed++;
		final int kmers=Tools.countKmers(r.bases, k);
		final double cKmers=r.quality==null ? kmers : Tools.countCorrectKmers(r.quality, k);
		kmersProcessed+=kmers;
		correctKmers+=cKmers;
		if(r.mapped() && r.bases!=null){
			if(USE_SECONDARY && r.sites!=null && r.sites.size()>0){
				boolean b=false;
				for(SiteScore ss : r.sites){
					b=processRead(r, ss) || b;
				}
				return b;
			}else{
				final Read mate=r.mate;
				final boolean set1=coords.set(r);
				final boolean set2=(PHYSICAL_COVERAGE && r.paired() && r.pairnum()==0 && coords2.set(mate));
				if(set1 && set2 && Tools.equals(coords.name, coords2.name)){
					final int start1=coords.start;
					final int stop1=coords.stop;
					final int start2=coords2.start;
					final int stop2=coords2.stop;
					final int strand=r.strand();
					final int length=Tools.max(stop1, stop2)-Tools.min(start1, start2)+1;
					mappedKmers+=kmers;
					addCoverage(new String(coords.name), null, null, Tools.min(start1, start2), Tools.max(stop1, stop2), length, 
							r.mappedNonClippedBases(), strand, 2-r.mateCount(), r.paired(), null);
				}else{
					if(set1){
						mappedKmers+=kmers;
						return addCoverage(new String(coords.name), r.bases, r.match, coords.start, coords.stop, r.length(), 
								r.mappedNonClippedBases(), coords.strand, 2-r.mateCount(), r.paired(), null);
					}
				}
			}
		}
		return false;
	}
	
	
	public boolean processRead(Read r, SiteScore ss){
		if(ss!=null && r.bases!=null){
			if(coords.set(ss)){
				return addCoverage(new String(coords.name), r.bases, ss.match, coords.start, coords.stop, r.length(), 
						r.mappedNonClippedBases(), coords.strand, 1, r.paired(), null);
			}
		}
		return false;
	}

	
	/*--------------------------------------------------------------*/
	/*----------------        Output Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public void printOutput(){
		
		totalScaffolds=list.size();
		
		long refKmers=0;
		long refBases2=0;
		for(Scaffold sc : list){
			refBases2+=sc.length;
			if(sc.length>=k){refKmers+=sc.length-k+1;}
		}

		String basecov1=(basecov==null ? null : (STRANDED ? basecov.replaceFirst("#", "1") : basecov));
		String rangecov1=(rangecov==null ? null : (STRANDED ? rangecov.replaceFirst("#", "1") : rangecov));
		String bincov1=(bincov==null ? null : (STRANDED ? bincov.replaceFirst("#", "1") : bincov));
		String normcov1=(normcov==null ? null : (STRANDED ? normcov.replaceFirst("#", "1") : normcov));
		String normcovOverall1=(normcovOverall==null ? null : (STRANDED ? normcovOverall.replaceFirst("#", "1") : normcovOverall));
		String histogram1=(histogram==null ? null : (STRANDED ? histogram.replaceFirst("#", "1") : histogram));
		String stats1=(covstats==null ? null : (STRANDED ? covstats.replaceFirst("#", "1") : covstats));

		String basecov2=(basecov==null || !STRANDED ? null : basecov.replaceFirst("#", "2"));
		String rangecov2=(rangecov==null || !STRANDED ? null : rangecov.replaceFirst("#", "2"));
		String bincov2=(bincov==null || !STRANDED ? null : bincov.replaceFirst("#", "2"));
		String normcov2=(normcov==null || !STRANDED ? null : normcov.replaceFirst("#", "2"));
		String normcovOverall2=(normcovOverall==null ? null : (STRANDED ? normcovOverall.replaceFirst("#", "2") : normcovOverall));
		String histogram2=(histogram==null || !STRANDED ? null : histogram.replaceFirst("#", "2"));
		String stats2=(covstats==null || !STRANDED ? null : covstats.replaceFirst("#", "2"));
		
		if(CONCISE){
			writeCoveragePerBaseConcise(basecov1, list, 0, minscaf);
			writeCoveragePerBaseConcise(basecov2, list, 1, minscaf);
		}else{
			writeCoveragePerBase(basecov1, list, DELTA_ONLY, 0, minscaf);
			writeCoveragePerBase(basecov2, list, DELTA_ONLY, 1, minscaf);
		}
		if(KEEP_SHORT_BINS){
			writeCoveragePerBaseBinned2(bincov1, list, binsize, 0, minscaf);
			writeCoveragePerBaseBinned2(bincov2, list, binsize, 1, minscaf);
		}else{
			writeCoveragePerBaseBinned(bincov1, list, binsize, 0, minscaf);
			writeCoveragePerBaseBinned(bincov2, list, binsize, 1, minscaf);
		}
		if(rangecov!=null) {
			writeCoverageRanges(rangecov1, list, 0, minscaf);
			writeCoverageRanges(rangecov2, list, 1, minscaf);
		}
		if(normcov!=null){
			writeCoveragePerBaseNormalized(normcov1, list, binsize, 0, minscaf);
			writeCoveragePerBaseNormalized(normcov2, list, binsize, 1, minscaf);
		}
		if(normcovOverall!=null){
			writeCoveragePerBaseNormalizedOverall(normcovOverall1, list, binsize, 0, minscaf);
			writeCoveragePerBaseNormalizedOverall(normcovOverall2, list, binsize, 1, minscaf);
		}
		if(outrpkm!=null){
			writeRPKM(outrpkm, in1, null, readsProcessed, NONZERO_ONLY,list);
		}
		
		{
			long[] hist=writeStats(stats1, 0);
			if(hist!=null){writeHist(histogram1, hist);}
			
			if(STRANDED){
				hist=writeStats(stats2, 1);
				if(hist!=null){writeHist(histogram2, hist);}
			}
		}
		
		final double mult=1.0/refBases;
		double depthCovered=mappedBases*mult;
		double depthCovered2=mappedNonClippedBases*mult;
		double depthCovered3=mappedBasesWithDels*mult;
		double pctScaffoldsWithCoverage=scaffoldsWithCoverage1*100.0/totalScaffolds;
		double pctCovered=totalCoveredBases1*100*mult;
		
		if(!KEY_VALUE){
			outstream.println("Reads:                               \t"+readsProcessed);
			outstream.println("Mapped reads:                        \t"+mappedReads);
			//		outstream.println("Mapped bases:                        \t"+mappedBases);
			//		outstream.println("Mapped non-clipped bases:            \t"+mappedNonClippedBases);
			outstream.println("Mapped bases:                        \t"+mappedNonClippedBases);
			outstream.println("Ref scaffolds:                       \t"+totalScaffolds);
			outstream.println("Ref bases:                           \t"+refBases);

			if(k>0){
				String kcovS=k+"-mer coverage:";
				String kcorrectS="Percent correct "+k+"-mers:";
				while(kcovS.length()<26){kcovS=kcovS+" ";}
				while(kcovS.length()<26){kcorrectS=kcorrectS+" ";}
				outstream.println(Tools.format("\n"+kcovS+"           \t%.3f", mappedKmers*1.0/refKmers));
				outstream.println(Tools.format(kcorrectS+"           \t%.3f", 100*correctKmers/kmersProcessed));
				//			outstream.println(kmersProcessed+", "+correctKmers);
			}

			outstream.println(Tools.format("\nPercent mapped:                      \t%.3f", mappedReads*100f/readsProcessed));
			outstream.println(Tools.format("Percent proper pairs:                \t%.3f", properPairs*100f/readsProcessed));
//			outstream.println(depthCovered);
			outstream.println(Tools.format("Average coverage:                    \t%.3f", depthCovered2));
			outstream.println(Tools.format("Average coverage with deletions:     \t%.3f", depthCovered3));
			if(USE_COVERAGE_ARRAYS && calcCovStdev){
				double[] stdev=standardDeviation(list, 0, minscaf);
				outstream.println(Tools.format("Standard deviation:                    \t%.3f", stdev[1]));
			}
			outstream.println(Tools.format("Percent scaffolds with any coverage: \t%.2f", pctScaffoldsWithCoverage));
			if(USE_COVERAGE_ARRAYS || USE_BITSETS){
				outstream.println(Tools.format("Percent of reference bases covered:  \t%.2f", pctCovered));
			}
		}else{
			outstream.println("reads="+readsProcessed);
			outstream.println("mappedReads="+mappedReads);
			outstream.println("mappedBases="+mappedNonClippedBases);
			outstream.println("mappedBasesWithDels="+mappedBasesWithDels);
			outstream.println("refScaffolds="+totalScaffolds);
			outstream.println("refBases="+refBases);
			outstream.println(Tools.format("percentMapped=%.3f", mappedReads*100f/readsProcessed));
			outstream.println(Tools.format("percentPaired=%.3f", properPairs*100f/readsProcessed));
			outstream.println(Tools.format("averageCoverage=%.3f", depthCovered2));
			if(USE_COVERAGE_ARRAYS && calcCovStdev){
				double[] stdev=standardDeviation(list, 0, minscaf);
				outstream.println(Tools.format("standardDeviation=%.3f", stdev[1]));
			}
			outstream.println(Tools.format("percentCoveredScaffolds=%.2f", pctScaffoldsWithCoverage));
			if(USE_COVERAGE_ARRAYS || USE_BITSETS){
				outstream.println(Tools.format("percentCoveredBases=%.2f", pctCovered));
			}
		}
	}
	
	public int basesUnderAverageCoverage(final int[] array, final double avg, final int window){
		if(array.length<window){return 0;}
		final long limit=(long)Math.ceil(window*avg);
		long covSum=0;
		int baseCount=0;
		for(int i=0; i<window; i++){
			covSum+=array[i];
		}
		
		boolean below=false;
		int lastStop=-1, lastStart=0;
		for(int a=0, b=window; b<array.length; a++, b++){
			if(covSum>=limit){
				if(below){//end range
					baseCount=b-Tools.max(lastStop+1, lastStart);
					lastStop=b-1;
					below=false;
				}
			}else{
				if(!below){//start range
					lastStart=a;
					below=true;
				}
			}
			covSum-=array[a];
			assert(covSum>=0);
			covSum+=array[b];
		}
		
		if(below){//end range
			baseCount=array.length-Tools.max(lastStop, lastStart);
		}
		
		assert(baseCount>=0);
		return baseCount;
	}
	
	public int basesUnderAverageCoverage(final char[] array, final double avg, final int window){
		if(array.length<window){return 0;}
		final long limit=(long)Math.ceil(window*avg);
		long covSum=0;
		int baseCount=0;
		for(int i=0; i<window; i++){
			covSum+=array[i];
		}
		
//		outstream.println("limit: "+limit);
		
		boolean below=false;
		int lastStop=-1, lastStart=0;
		for(int a=0, b=window; b<array.length; a++, b++){
			if(covSum>=limit){
				if(below){//end range
					baseCount+=b-Tools.max(lastStop+1, lastStart);
					
//					outstream.println("\nprev: "+lastStop+", "+lastStart);
//					outstream.println("end range at "+a+", "+b);
//					outstream.println("baseCount: "+baseCount+", covSum="+covSum);
					
					lastStop=b-1;
					below=false;
				}
			}else{
				if(!below){//start range
					
//					outstream.println("\nprev: "+lastStop+", "+lastStart);
//					outstream.println("start range at "+a+", "+b);
//					outstream.println("baseCount: "+baseCount+", covSum="+covSum);
					
					lastStart=a;
					below=true;
				}
			}
			covSum-=array[a];
			assert(covSum>=0);
			covSum+=array[b];
		}
		
		if(below){//end range
			baseCount+=array.length-Tools.max(lastStop+1, lastStart);
			
//			outstream.println("\nprev: "+lastStop+", "+lastStart);
//			outstream.println("end range at "+array.length);
//			outstream.println("baseCount: "+baseCount+", covSum="+covSum);
		}
		
		assert(baseCount>=0);
		return baseCount;
	}
	
	public long[] writeStats(String fname, int strand){
//		outstream.println("Writing stats for "+fname+", "+strand);
		final TextStreamWriter tsw=(fname==null ? null : new TextStreamWriter(fname, overwrite, false, true));
		
		if(tsw!=null){
			tsw.start();
			if(printHeader){
				String pound=(headerPound ? "#" : "");
				if(TWOCOLUMN){
					tsw.println(pound+"ID\tAvg_fold");
				}else{
					tsw.println(pound+"ID\tAvg_fold\tLength\tRef_GC\tCovered_percent\tCovered_bases\tPlus_reads\tMinus_reads"+
							(COUNT_GC ? "\tRead_GC" : "")+
							(USE_COVERAGE_ARRAYS ? ("\tMedian_fold\tStd_Dev") : "")+
							(USE_WINDOW ? "\tUnder_"+Tools.format("%.0f",LOW_COV_DEPTH)+"/"+LOW_COV_WINDOW : ""));
				}
			}
			//Maximally:
			//"ID\tAvg_fold\tLength\tRef_GC\tCovered_percent\tCovered_bases\tPlus_reads\tMinus_reads\tRead_GC\tMedian_fold\tStd_Dev\tUnder_X"
		}
		
		final int histmax=(bits32 ? 1000000 : Character.MAX_VALUE);
		final LongList hist=new LongList(Character.MAX_VALUE);
		
		long coveredScafTemp=0;
		long coveredBaseTemp=0;
		for(Scaffold scaf : list){
			final long sum=scaf.basehits;
			int covered=0;
			int median=0;
			int underWindowAverage=0;
			final double stdev;
			if(USE_COVERAGE_ARRAYS){
				CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
				if(ca!=null){
					for(int i=0; i<scaf.length; i++){
						int x=ca.get(i);
						hist.increment(Tools.min(x, histmax));
//						sum+=x;
						if(x>=minDepthToBeCovered){covered++;}
					}
					if(bits32){
						int[] array=((CoverageArray3)ca).array;
						stdev=Tools.standardDeviation(array);
						underWindowAverage=basesUnderAverageCoverage(array, LOW_COV_DEPTH, LOW_COV_WINDOW);
						Shared.sort(array);
						Tools.reverseInPlace(array);
						median=ca.get(scaf.length/2);
					}else{
						char[] array=((CoverageArray2)ca).array;
						stdev=Tools.standardDeviation(array);
						underWindowAverage=basesUnderAverageCoverage(array, LOW_COV_DEPTH, LOW_COV_WINDOW);
						Arrays.sort(array);
						Tools.reverseInPlace(array);
						median=ca.get(scaf.length/2);
					}
				}else{
					hist.increment(0, scaf.length);
					stdev=0;
				}
			}else if(USE_BITSETS){
//				sum+=scaf.basehits;
				BitSet bs=(BitSet)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
				covered=(bs==null ? 0 : bs.cardinality());
				stdev=-1;
			}else{
				stdev=-1;
			}
			
			if(sum>0){
				coveredScafTemp++;
			}
			//			pw.print(scaf.name);
			if(tsw!=null && (sum>0 || !NONZERO_ONLY) && scaf.length>=minscaf){
				if(TWOCOLUMN){
					tsw.print(Tools.format("%s\t%.4f\n", scaf.name, sum/(double)scaf.length));
				}else if(COUNT_GC){
					long[] bc=scaf.basecount;
					double gc=(bc[1]+bc[2])*1d/Data.max(1, bc[0]+bc[1]+bc[2]+bc[3]);
					if(USE_COVERAGE_ARRAYS){
						if(USE_WINDOW){
							
							//Variable portion:
							//"\tRead_GC\tMedian_fold\tStd_Dev\tUnder_X\n"
							//\t%.4f\t%d\t%.2f\t%d\n
							
							tsw.print(Tools.format("%s\t%.4f\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%.4f\t%d\t%.2f\t%d\n", scaf.name, sum/(double)scaf.length, scaf.length,
									scaf.gc, covered*100d/scaf.length, covered, (scaf.readhits-scaf.readhitsMinus), scaf.readhitsMinus,
									gc, median, stdev, underWindowAverage));
						}else{
							tsw.print(Tools.format("%s\t%.4f\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%.4f\t%d\t%.2f\n", scaf.name, sum/(double)scaf.length, scaf.length,
									scaf.gc, covered*100d/scaf.length, covered, (scaf.readhits-scaf.readhitsMinus), scaf.readhitsMinus,
									gc, median, stdev));
						}
					}else{
						tsw.print(Tools.format("%s\t%.4f\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%.4f\n", scaf.name, sum/(double)scaf.length, scaf.length,
								scaf.gc, covered*100d/scaf.length, covered, (scaf.readhits-scaf.readhitsMinus), scaf.readhitsMinus,
								gc));
					}
				}else{
					if(USE_COVERAGE_ARRAYS){
						if(USE_WINDOW){
							tsw.print(Tools.format("%s\t%.4f\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%.2f\t%d\n", scaf.name, sum/(double)scaf.length, scaf.length,
									scaf.gc, covered*100d/scaf.length, covered, (scaf.readhits-scaf.readhitsMinus), scaf.readhitsMinus,
									median, stdev, underWindowAverage));
						}else{
							tsw.print(Tools.format("%s\t%.4f\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%.2f\n", scaf.name, sum/(double)scaf.length, scaf.length,
									scaf.gc, covered*100d/scaf.length, covered, (scaf.readhits-scaf.readhitsMinus), scaf.readhitsMinus,
									median, stdev));
						}
					}else{
						tsw.print(Tools.format("%s\t%.4f\t%d\t%.4f\t%.4f\t%d\t%d\t%d\n", scaf.name, sum/(double)scaf.length, scaf.length,
								scaf.gc, covered*100d/scaf.length, covered, (scaf.readhits-scaf.readhitsMinus), scaf.readhitsMinus));
					}
				}
			}
			coveredBaseTemp+=covered;
		}
		
		if(strand==0){
			scaffoldsWithCoverage1+=coveredScafTemp;
			totalCoveredBases1+=coveredBaseTemp;
		}else{
			scaffoldsWithCoverage2+=coveredScafTemp;
			totalCoveredBases2+=coveredBaseTemp;
		}
		
		if(tsw!=null){tsw.poisonAndWait();}
		return hist==null ? null : hist.array;
	}
	
	/**
	 * Write a histogram of number of bases covered to each depth
	 * @param fname Output filename
	 * @param counts counts[X] stores the number of bases with coverage X
	 */
	public static void writeHist(String fname, long[] counts){
		if(fname==null){return;}
		assert(counts!=null) : "Can't write a histogram with null counts.";
		ByteStreamWriter tsw=new ByteStreamWriter(fname, overwrite, false, false);
		tsw.start();
		if(printHeader){
			if(headerPound){tsw.print('#');}
			tsw.println("Coverage\tnumBases");
		}
		int max=0;
		for(max=counts.length-1; max>0 && counts[max]==0; max--){}
		for(int i=0; i<=max; i++){
			long x=counts[i];
			tsw.print(i);
			tsw.print('\t');
			tsw.println(x);
		}
		
		tsw.poisonAndWait();
	}
	
	/**
	 * Prints coverage in this format:
	 * scafname TAB position TAB coverage
	 * scafname TAB position TAB coverage
	 * @param fname Output filename
	 * @param list List of reference scaffolds
	 * @param deltaOnly Only write lines when coverage changes
	 * @param strand Only use coverage from reads mapped to this strand (0=plus, 1=minus)
	 */
	public void writeCoveragePerBase(String fname, ArrayList<Scaffold> list, boolean deltaOnly, int strand, int minscaf){
		if(fname==null || (!STRANDED && strand>0)){return;}
		
		if(verbose){outstream.println("Starting tsw "+fname);}
		ByteStreamWriter tsw=new ByteStreamWriter(fname, overwrite, false, true);
		if(verbose){outstream.println("Created tsw "+fname);}
		tsw.start();
//		if(verbose){outstream.println("Started tsw "+fname);}
		if(printHeader){
			if(headerPound){tsw.print('#');}
			tsw.println("RefName\tPos\tCoverage");
		}
		
		for(Scaffold scaf : list){
			int last=-1;
			CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
			if(scaf.length>=minscaf){
				for(int i=0, len=scaf.length; i<len; i++){
					int x=(ca==null ? 0 : ca.get(i));
					if(!deltaOnly || x!=last){
						if(x>0 || !NONZERO_ONLY){
							tsw.print(scaf.name);
							tsw.print('\t');
							tsw.print(i);
							tsw.print('\t');
							tsw.println(x);
						}
						last=x;
					}
				}
			}
		}

		if(verbose){outstream.println("Closing tsw "+fname);}
		tsw.poisonAndWait();
		if(verbose){outstream.println("Closed tsw "+fname);}
	}
	
	/**
	 * Prints coverage in this format, skipping zero-coverage positions:
	 * #scafname
	 * position TAB coverage
	 * position TAB coverage
	 * @param fname Output filename
	 * @param list List of reference scaffolds
	 * @param strand Only use coverage from reads mapped to this strand (0=plus, 1=minus)
	 */
	public void writeCoveragePerBaseConcise(String fname, ArrayList<Scaffold> list, int strand, int minscaf){
		if(fname==null || (!STRANDED && strand>0)){return;}
		
		if(verbose){outstream.println("Starting tsw "+fname);}
		ByteStreamWriter tsw=new ByteStreamWriter(fname, overwrite, false, true);
		tsw.start();
		if(verbose){outstream.println("Started tsw "+fname);}
//		tsw.print(pound+"RefName\tPos\tCoverage\n");
		
		for(Scaffold scaf : list){
			tsw.print('#');
			tsw.println(scaf.name);
			CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
			if(scaf.length>=minscaf){
				for(int i=0; i<scaf.length; i++){
					int x=(ca==null ? 0 : ca.get(i));
					if(x>0){
						tsw.print(i);
						tsw.print('\t');
						tsw.println(x);
					}
				}
			}
		}

		if(verbose){outstream.println("Closing tsw "+fname);}
		tsw.poisonAndWait();
		if(verbose){outstream.println("Closed tsw "+fname);}
//		assert(false);
	}
	
	/**
	 * Prints coverage in this format, skipping zero-coverage positions:
	 * #scafname
	 * position TAB coverage
	 * position TAB coverage
	 * @param fname Output filename
	 * @param list List of reference scaffolds
	 * @param strand Only use coverage from reads mapped to this strand (0=plus, 1=minus)
	 */
	public void writeCoverageRanges(String fname, ArrayList<Scaffold> list, int strand, int minscaf){
		if(fname==null || (!STRANDED && strand>0)){return;}
		
		if(verbose){outstream.println("Starting tsw "+fname);}
		ByteStreamWriter tsw=new ByteStreamWriter(fname, overwrite, false, true);
		tsw.start();
		if(verbose){outstream.println("Started tsw "+fname);}
//		tsw.print(pound+"RefName\tPos\tCoverage\n");
		
		for(Scaffold scaf : list){
			tsw.print('#');
			tsw.println(scaf.name);
			CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
			if(scaf.length>=minscaf){
				int start=-1;
				int prev=-1;
				long sum=0;
				
				for(int i=0; i<scaf.length; i++){
					int x=(ca==null ? 0 : ca.get(i));
					if(x>=minDepthToBeCovered){
						if(prev<minDepthToBeCovered) {//open range
							start=i;
						}
						sum+=x;
					}else if(prev>=minDepthToBeCovered){//close range
						float len=i-start+1;
						tsw.print(start).print('-').print(i-1).tab().println(sum/len, 2);
						sum=0;
					}
					prev=x;
				}
				if(prev>=minDepthToBeCovered){//close final range
					float len=scaf.length-start;
					tsw.print(start).print('-').print(scaf.length-1).tab().println(sum/len, 2);
					sum=0;
				}
			}
		}

		if(verbose){outstream.println("Closing tsw "+fname);}
		tsw.poisonAndWait();
		if(verbose){outstream.println("Closed tsw "+fname);}
//		assert(false);
	}
	
	/**
	 * Note.  As written, this will truncate all trailing bases of each scaffold's length modulo binsize.
	 * For example, with binsize 1000, the last 500 bases of a 1500 base scaffold will be ignored.
	 * @param fname Output filename
	 * @param list List of reference scaffolds
	 * @param binsize Width of coverage bins in bp
	 * @param strand Only use coverage from reads mapped to this strand (0=plus, 1=minus)
	 */
	public static void writeCoveragePerBaseBinned(String fname, ArrayList<Scaffold> list, int binsize, int strand, int minscaf){
		if(fname==null || (!STRANDED && strand>0)){return;}
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		if(printHeader){
			String pound=(headerPound ? "#" : "");
			if(calcCovStdev){
				double[] stdev=standardDeviation(list, strand, minscaf);
				if(stdev!=null){
					tsw.print(pound+"Mean\t"+Tools.format("%.3f", stdev[0])+"\n");
					tsw.print(pound+"STDev\t"+Tools.format("%.3f", stdev[1])+"\n");
				}
			}
			tsw.print(pound+"RefName\tCov\tPos\tRunningPos\n");
		}
		
		long running=0;
		final float invbin=1f/binsize;
		for(Scaffold scaf : list){
			if(scaf.length>=binsize && scaf.length>=minscaf){
				CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
				int lastPos=-1, nextPos=binsize-1;
				long sum=0;
				for(int i=0; i<scaf.length; i++){
					int x=(ca==null ? 0 : ca.get(i));
					sum+=x;
					if(i>=nextPos){
						if(sum>0 || !NONZERO_ONLY){
							tsw.print(Tools.format("%s\t%.2f\t%d\t%d\n", scaf.name, sum*invbin, (i+1), running));
						}
						lastPos=i;
						running+=binsize;
						nextPos+=binsize;
						sum=0;
					}
				}
			}
		}
		
		tsw.poisonAndWait();
	}
	
	/**
	 * This version will NOT truncate all trailing bases of each scaffold's length modulo binsize.
	 * @param fname Output filename
	 * @param list List of reference scaffolds
	 * @param binsize Width of coverage bins in bp
	 * @param strand Only use coverage from reads mapped to this strand (0=plus, 1=minus)
	 */
	public static void writeCoveragePerBaseBinned2(String fname, ArrayList<Scaffold> list, int binsize, int strand, int minscaf){
		if(fname==null || (!STRANDED && strand>0)){return;}
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		if(printHeader){
			String pound=(headerPound ? "#" : "");
			if(calcCovStdev){
				double[] stdev=standardDeviationBinned(list, binsize, strand, minscaf);
				if(stdev!=null){
					tsw.print(pound+"Mean\t"+Tools.format("%.3f", stdev[0])+"\n");
					tsw.print(pound+"STDev\t"+Tools.format("%.3f", stdev[1])+"\n");
				}
			}
			tsw.print(pound+"RefName\tCov\tPos\tRunningPos\n");
		}
		
		long running=0;
		for(Scaffold scaf : list){
			CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
			int lastPos=-1, nextPos=binsize-1;
			long sum=0;
			final int lim=scaf.length-1;
			for(int i=0; i<scaf.length; i++){
				int x=(ca==null ? 0 : ca.get(i));
				sum+=x;
				if(i>=nextPos || i==lim){
					int bin=(i-lastPos);
					if(scaf.length>=minscaf){
						if(sum>0 || !NONZERO_ONLY){
							tsw.print(Tools.format("%s\t%.2f\t%d\t%d\n", scaf.name, sum/(float)bin, (i+1), running));
						}
					}
					running+=bin;
					nextPos+=binsize;
					lastPos=i;
					sum=0;
				}
			}
		}
		
		tsw.poisonAndWait();
	}

	//Unsafe because it will fail if there are over 2 billion bins
	public static double[] standardDeviationBinnedUnsafe(String fname, ArrayList<Scaffold> scaffolds, int binsize, int strand, int minscaf){
		
		LongList list=new LongList();
		for(Scaffold scaf : scaffolds){
			CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
			int lastPos=-1, nextPos=binsize-1;
			long sum=0;
			final int lim=scaf.length-1;
			for(int i=0; i<scaf.length; i++){
				int x=(ca==null ? 0 : ca.get(i));
				sum+=x;
				if(i>=nextPos || i==lim){
					int bin=(i-lastPos);
					if(scaf.length>=minscaf){
						list.add((int)(10*(sum/(double)bin)));
					}
					nextPos+=binsize;
					lastPos=i;
					sum=0;
				}
			}
		}
		list.sort();
		
		double mean=0.1*list.mean();
		double median=0.1*list.median();
		double mode=0.1*list.mode();
		double stdev=0.1*list.stdev();
		return new double[] {mean, median, mode, stdev};
	}

	public static double[] standardDeviationBinned(ArrayList<Scaffold> scaffolds, int binsize, int strand, int minscaf){
		double totalSum=0;
		long bins=0;
		
		for(Scaffold scaf : scaffolds){
			if(scaf.length>=minscaf){
				CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
				int lastPos=-1, nextPos=binsize-1;
				long tempSum=0;
				final int lim=scaf.length-1;
				for(int i=0; i<scaf.length; i++){
					int x=(ca==null ? 0 : ca.get(i));
					tempSum+=x;
					if(i>=nextPos || i==lim){
						int bin=(i-lastPos);
						double depth=(tempSum/(double)bin);
						totalSum+=depth;
						bins++;
						nextPos+=binsize;
						lastPos=i;
						tempSum=0;
					}
				}
			}
		}
		
		if(bins<1){return new double[] {0, 0};}
		final double mean=totalSum/(double)bins;
		double sumdev2=0;
		
		for(Scaffold scaf : scaffolds){
			if(scaf.length>=minscaf){
				CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
				int lastPos=-1, nextPos=binsize-1;
				long tempSum=0;
				final int lim=scaf.length-1;
				for(int i=0; i<scaf.length; i++){
					int x=(ca==null ? 0 : ca.get(i));
					tempSum+=x;
					if(i>=nextPos || i==lim){
						int bin=(i-lastPos);
						double depth=(tempSum/(double)bin);
						double dev=mean-depth;
						sumdev2+=(dev*dev);
						nextPos+=binsize;
						lastPos=i;
						tempSum=0;
					}
				}
			}
		}
		
		final double stdev=Math.sqrt(sumdev2/bins);
		return new double[] {mean, stdev};
	}

	public static double[] standardDeviation(ArrayList<Scaffold> scaffolds, int strand, int minscaf){
		long totalSum=0, bins=0, scafsCounted=0;
		for(Scaffold scaf : scaffolds){
//			System.err.println(scaf.name+", "+scaf.length+"/"+minscaf);
			if(scaf.length>=minscaf){
				final CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
				bins+=scaf.length;
				scafsCounted++;
				if(ca!=null) {
					for(int i=0; i<scaf.length; i++){
						int depth=(ca==null ? 0 : ca.get(i));
						totalSum+=depth;
					}
				}
			}
		}
		
		if(bins<1){return new double[] {0, 0};}
		final double mean=totalSum/(double)bins;
		double sumdev2=0;
//		System.err.println("sum="+totalSum+", bins="+bins+", mean="+mean+", scaffolds="+scafsCounted);
		for(Scaffold scaf : scaffolds){
			if(scaf.length>=minscaf){
				final CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
				double sumTemp=0;
				if(ca!=null) {
					for(int i=0; i<scaf.length; i++){
						int depth=(ca==null ? 0 : ca.get(i));
						double dev=mean-depth;
						sumTemp+=(dev*dev);
					}
				}else {
					sumTemp+=(mean*mean*scaf.length);
				}
				sumdev2+=sumTemp;
			}
		}

		final double stdev=Math.sqrt(sumdev2/bins);
//		System.err.println("sumdev2="+sumdev2+", stdev="+stdev);
		return new double[] {mean, stdev};
	}
	
	
	/**
	 * @param fname Output filename
	 * @param list List of reference scaffolds
	 * @param binsize Width of coverage bins in bp
	 * @param strand Only use coverage from reads mapped to this strand (0=plus, 1=minus)
	 */
	public static void writeCoveragePerBaseNormalized(String fname, ArrayList<Scaffold> list, double binsize, int strand, int minscaf){
		if(fname==null || (!STRANDED && strand>0)){return;}
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		if(printHeader){
			String pound=(headerPound ? "#" : "");
			tsw.print(pound+"RefName\tBin\tCov\tPos\tRunningPos\n");
		}
		
		double running=0;
		double invbin=1.0/binsize;
		final double invbincount=1.0/NORMALIZE_LENGTH_BINS;
		for(Scaffold scaf : list){
			if(NORMALIZE_LENGTH_BINS>0){
				binsize=scaf.length*invbincount;
				invbin=1.0/binsize;
			}
			
			if(scaf.length>=binsize && scaf.length>=minscaf){
				long max=-1;
				
				final CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
				double lastPos=-1, nextPos=binsize-1;
				long sum=0;

				if(NORMALIZE_COVERAGE){
					for(int i=0; i<scaf.length; i++){
						int x=(ca==null ? 0 : ca.get(i));
						sum+=x;
						if(i>=nextPos){
							max=Tools.max(sum, max);
							running+=binsize;
							nextPos+=binsize;
							sum=0;
						}
					}
					lastPos=-1;
					nextPos=binsize-1;
					sum=0;
					assert(max>-1) : max;
				}
				max=Tools.max(max, 1);
				final double binmult=(NORMALIZE_COVERAGE ? 1d/max : invbin);
				
//				assert(false) : NORMALIZE_COVERAGE+", "+binmult+", "+invbin+", "+max+", "+binsize;
				
				final String formatString=NORMALIZE_COVERAGE ? "%s\t%d\t%.5f\t%d\t%d\n" : "%s\t%d\t%.2f\t%d\t%d\n";
				int bin=1;
				for(int i=0; i<scaf.length; i++){
					int x=(ca==null ? 0 : ca.get(i));
					sum+=x;
					if(i>=nextPos){
//						outstream.println(x+", "+i+", "+nextPos+", "+sum+", "+(sum*binmult));
						tsw.print(String.format(formatString, scaf.name, bin, sum*binmult, (i+1), (long)running));
						bin++;
						lastPos=i;
						running+=binsize;
						nextPos+=binsize;
						sum=0;
					}
				}
			}
		}
		
		tsw.poisonAndWait();
	}
	

	
	/**
	 * @param fname Output filename
	 * @param list List of reference scaffolds
	 * @param binsize Width of coverage bins in bp
	 * @param strand Only use coverage from reads mapped to this strand (0=plus, 1=minus)
	 */
	public static void writeCoveragePerBaseNormalizedOverall(String fname, ArrayList<Scaffold> list, double binsize, int strand, int minscaf){
		if(fname==null || (!STRANDED && strand>0)){return;}
		
		assert(NORMALIZE_LENGTH_BINS>0) : "Must set 'normalizebins' flag to a positive integer.";
		double running=0;
		double invbin=1.0/binsize;
		long usedScafs=0;
		final double invbincount=1.0/NORMALIZE_LENGTH_BINS;

		double[] normalized=new double[NORMALIZE_LENGTH_BINS+1];
		double[] absolute=new double[NORMALIZE_LENGTH_BINS+1];
		
		for(Scaffold scaf : list){
			if(NORMALIZE_LENGTH_BINS>0){
				binsize=scaf.length*invbincount;
				invbin=1.0/binsize;
			}
			
			if(scaf.length>=binsize && scaf.length>=minscaf){
				usedScafs++;

				if(scaf.readhits>0){
					long max=-1;
					final CoverageArray ca=(CoverageArray)(STRANDED && strand==1 ? scaf.obj1 : scaf.obj0);
					double lastPos=-1, nextPos=binsize-1;
					long sum=0;

					{
						for(int i=0; i<scaf.length; i++){
							int x=(ca==null ? 0 : ca.get(i));
							sum+=x;
							if(i>=nextPos){
								max=Tools.max(sum, max);
								running+=binsize;
								nextPos+=binsize;
								sum=0;
							}
						}
						lastPos=-1;
						nextPos=binsize-1;
						sum=0;
						assert(max>-1) : max;
					}
					max=Tools.max(max, 1);
					final double binmult=1d/max;
					
					//				assert(false) : NORMALIZE_COVERAGE+", "+binmult+", "+invbin+", "+max+", "+binsize;

					int bin=1;
					for(int i=0; i<scaf.length; i++){
						int x=(ca==null ? 0 : ca.get(i));
						sum+=x;
						if(i>=nextPos){
							normalized[bin]+=(sum*binmult);
							absolute[bin]+=(sum*invbin);
							bin++;
							lastPos=i;
							running+=binsize;
							nextPos+=binsize;
							sum=0;
						}
					}
				}
			}
		}
		
		TextStreamWriter tsw=new TextStreamWriter(fname, overwrite, false, false);
		tsw.start();
		if(printHeader){
			String pound=(headerPound ? "#" : "");
			tsw.print(pound+"RefName\tBin\tAbsCov\tNormCov\n");
		}
		double invScafs=1d/Tools.max(1, usedScafs);
		
		final double maxNorm=Tools.max(normalized);
		final double normMult=1/maxNorm;
		
		for(int bin=1; bin<normalized.length; bin++){
//			assert((absolute[bin]*invScafs)!=Double.NaN && (normalized[bin]*invScafs)!=Double.NaN) : invScafs+", "+absolute[bin]+", "+normalized[bin];
//			assert(false) : invScafs+", "+absolute[bin]+", "+normalized[bin]+", "+(absolute[bin]*invScafs)+", "+(normalized[bin]*invScafs);
			tsw.print(Tools.format("%s\t%d\t%.5f\t%.5f\n", "all", bin, absolute[bin]*invScafs, normalized[bin]*normMult));
		}
		
		tsw.poisonAndWait();
	}
	

	
	/**
	 * Write RPKM statistics.
	 */
	public static void writeRPKM(String out, String in1, String in2, long readsIn, boolean printNonZeroOnly, ArrayList<Scaffold> list){
		if(out==null){return;}
		final TextStreamWriter tsw=new TextStreamWriter(out, overwrite, false, false);
		tsw.start();

		/* Count mapped reads */
		long mappedReads=0;
		long mappedFrags=0;
		for(Scaffold scaf : list){
			mappedReads+=scaf.readhits;
			mappedFrags+=scaf.fraghits;
		}
		mappedFrags/=2;
		
		/* Print header */
		tsw.print("#File\t"+(in1==null ? "" : in1)+(in2==null ? "" : "\t"+in2)+"\n");
		tsw.print(Tools.format("#Reads\t%d\n",readsIn));
		tsw.print(Tools.format("#Mapped\t%d\n",mappedReads));
		tsw.print(Tools.format("#RefSequences\t%d\n",list.size()));
		tsw.print("#Name\tLength\tBases\tCoverage\tReads\tRPKM\tFrags\tFPKM\n");
		
		final float readMult=1000000000f/Tools.max(1, mappedReads);
		final float fragMult=1000000000f/Tools.max(1, mappedFrags);
		
		/* Print data */
		for(final Scaffold scaf : list){
			final long reads=scaf.readhits;
			final long frags=scaf.fraghits/2;
			final long bases=scaf.basehits;
			final String s=scaf.name;
			final int len=scaf.length;
			final double invlen=1.0/Tools.max(1, len);
			final double readMult2=readMult*invlen;
			final double fragMult2=fragMult*invlen;
			if(reads>0 || !printNonZeroOnly){
				tsw.print(Tools.format("%s\t%d\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f\n",s,len,bases,bases*invlen,reads,reads*readMult2,frags,frags*fragMult2));
			}
		}
		tsw.poisonAndWait();
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** The list of all scaffolds */
	private ArrayList<Scaffold> list;
	/** Maps names to scaffolds */
	private HashMap<String, Scaffold> table;
	/** Converts BBMap index coordinates to scaffold coordinates */
	private final ScaffoldCoordinates coords=new ScaffoldCoordinates(), coords2=new ScaffoldCoordinates();
	
	/** Mapping program name */
	private String program=null;
	/** Mapping program version */
	private String version=null;

	//Inputs
	/** Primary input file (typically sam) */
	public String in1=null;
	/** Secondary input file (typically for coverage directly from BBMap) */
	public String in2=null;
	
	public ArrayList<String> covIn=new ArrayList<String>();
	/** Optional, for calculating GC */
	public String reference=null;
	public String orffasta=null;

	//Outputs
	/** Stream unaltered sam input to this output */
	public String outsam=null;
	/** Coverage statistics, one line per scaffold */
	public String covstats=null;
	public String outorf=null;
	/** Coverage histogram, one line per depth and one point per base */
	public String histogram=null;
	/** Coverage with one line per base */
	public String basecov=null;
	/** Coverage with one file per scaffold */
	public String basecov_ps=null;
	/** Coverage with one line per bin */
	public String bincov=null;
	/** Coverage with one line per bin, normalized by length and/or height */
	public String normcov=null;
	/** Coverage with one line per bin, normalized by length and/or height, for combined reference */
	public String normcovOverall=null;
	/** rpkm/fpkm output, similar to Seal */
	public String outrpkm=null;
	
	public String rangecov=null;
	
	/** Typically indicates that a header line was encountered in an unexpected place, e.g. with concatenated sam files. */
	private boolean error=false;
	
	private boolean warned=false;
	private final boolean EA=Shared.EA();
	
	/** Total length of reference */
	public long refBases=0;
	public long mappedBases=0;
	public long mappedNonClippedBases=0;
	public long mappedBasesWithDels=0;
	public long mappedReads=0;
	public long properPairs=0;
	public long readsProcessed=0;
	public long basesProcessed=0;
	public long kmersProcessed=0;
	public long mappedKmers=0;
	public double correctKmers=0;
	public long totalCoveredBases1=0;
	public long totalCoveredBases2=0;
	public long scaffoldsWithCoverage1=0;
	public long scaffoldsWithCoverage2=0;
	public long totalScaffolds=0;

	public int k=0;
	
	//Don't reset these variables when clearing.
	public long maxReads=-1;
	public int initialScaffolds=4096;
	public int binsize=1000;
	public boolean bits32=false;
	public int minMapq=0;
	public boolean useStreamer=true;
	public int streamerThreads=2;
	public int minDepthToBeCovered=1;
	
	private boolean qtrimLeft=false;
	private boolean qtrimRight=false;
	private float trimq=-1;
	private final float trimE;
	private int border=0;
	
	/** Don't print coverage info for scaffolds shorter than this */
	public int minscaf=0;
	
	public HashMap<String, SamLine> pairTable=new HashMap<String, SamLine>();
	
	public PrintStream outstream=System.err;
	
	private boolean errorState=false;
	
	private final LineParser1 lp=new LineParser1('\t');
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	

	/** Permission to overwrite existing files */
	public static boolean overwrite=true;
	/** Permission to append to existing files */
	public static boolean append=false;
	/** Print verbose log messages */
	public static boolean verbose=false;
	/** Print the arguments to main */
	public static boolean printCommand=true;
	
	/** Print headers in output files */
	public static boolean printHeader=true;
	/** Prepend '#' symbol to header lines */
	public static boolean headerPound=true;
	/** Calculate standard deviation of coverage */
	public static boolean calcCovStdev=true;
	
	/** Window size to use when calculating average coverage,
	 * for detecting contiguous low-coverage areas */
	public static int LOW_COV_WINDOW=500;
	/** Min average coverage to not be classified as low-depth */
	public static double LOW_COV_DEPTH=5;
	/** Print number of bases below a certain average coverage in a window */
	public static boolean USE_WINDOW=false;
	
	/** Track base composition of reads covering each scaffold */
	public static boolean COUNT_GC=true;
	/** Output in 2-column format ("#ID\tAvg_fold\n") */
	public static boolean TWOCOLUMN=false;
	/** Track coverage for strands independently */
	public static boolean STRANDED=false;
	/** Add scaffold information from the reference (in addition to sam header) */
	public static boolean ADD_FROM_REF=true;
	/** Add scaffold information from reads mapped to unknown scaffolds (in addition to sam header) */
	public static boolean ADD_FROM_READS=false;
	/** Only print scaffolds with nonzero coverage */
	public static boolean NONZERO_ONLY=false;
	/** Store coverage info in numeric arrays */
	public static boolean USE_COVERAGE_ARRAYS=true;
	/** Store coverage info in bitsets */
	public static boolean USE_BITSETS=false;
	/** Only print lines when coverage changes (for compatibility with Jigsaw) */
	public static boolean DELTA_ONLY=false;
	/** Process secondary alignments */
	public static boolean USE_SECONDARY=true;
	/** Include coverage of unsequenced middle portion of pairs */
	public static boolean PHYSICAL_COVERAGE=false;
	/** Use 'tlen' field when calculating physical coverage */
	public static boolean USE_TLEN=true;
	/** Abort on error; otherwise, errors may be ignored */
	public static boolean ABORT_ON_ERROR=true;
	/** Print coverage for the last bin of a scaffold, even if it is shorter than binsize */
	public static boolean KEEP_SHORT_BINS=true;
	/** Only track coverage for start location */
	public static boolean START_ONLY=false;
	/** Only track coverage for stop location */
	public static boolean STOP_ONLY=false;
	/** This appears to be the same as nonzeroonly... */
	public static boolean CONCISE=false;
	/** Normalize coverage by expression contig coverage as a fraction of its max coverage */
	public static boolean NORMALIZE_COVERAGE=false;
	/** Normalize contig length by binning into this many bins per contig */
	public static int NORMALIZE_LENGTH_BINS=100;
	/** Include soft-clipped bases in coverage */
	public static boolean INCLUDE_SOFT_CLIP=false;
	/** Include deletions/introns in coverage */
	public static boolean INCLUDE_DELETIONS=true;
	/** Include reads flagged as duplicates in coverage */
	public static boolean INCLUDE_DUPLICATES=true;

	public static boolean KEY_VALUE;
	
	/** Translation array for tracking base counts */
	private static final byte[] charToNum=AssemblyStats2.makeCharToNum();
	
	private static final int NOTHING_MODE=0, BITSET_MODE=1, ARRAY_MODE=2;
	
}
