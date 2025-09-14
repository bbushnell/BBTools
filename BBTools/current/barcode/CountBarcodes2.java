package barcode;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map.Entry;

import barcode.stub.PCRMatrixProb;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import hiseq.IlluminaHeaderParser2;
import shared.LineParser;
import shared.LineParser1;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import structures.LongPair;
import tracker.ReadStats;

/**
 * Counts barcodes in a fastq file and produces a summary.
 * 
 * Sample command:
 * countbarcodes2.sh pcrmatrix type=hdist in=sub005.fq.gz ow expected=expected.txt -Xmx31g t=32 
 * countbarcodes2.sh pcrmatrix type=hdist in=sub005.fq.gz ow expected=expected.txt -Xmx31g t=32 
 *    outcontam=contam2.txt quantset=quantset.txt
 * 
 * @author Brian Bushnell
 * @date June 20, 2014
 *
 */
public class CountBarcodes2 {
	
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
		CountBarcodes2 x=new CountBarcodes2(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CountBarcodes2(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		Shared.capBuffers(4); //Only for singlethreaded programs
		
		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			setInterleaved=parser.setInterleaved;

			if(parser.in1!=null) {
				for(String s : parser.in1.split(",")) {
					in1.add(s);
				}
			}
			if(parser.in2!=null) {
				for(String s : parser.in2.split(",")) {
					in2.add(s);
				}
			}
			extin=parser.extin;

			out1=parser.out1;
			extout=parser.extout;
		}

		PCRMatrix.postParseStatic();
		doPoundReplacement(); //Replace # with 1 and 2
		adjustInterleaving(); //Make sure interleaving agrees with number of input and output files
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, extout, true, overwrite, append, false);

		//Create input FileFormat objects
		for(String s : in1) {
			ffin1.add(FileFormat.testInput(s, FileFormat.FASTQ, extin, true, true));
		}
		for(String s : in2) {
			ffin2.add(FileFormat.testInput(s, FileFormat.FASTQ, extin, true, true));
		}
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
			}else if(a.equals("expected") || a.equals("valid") || a.equals("names") || 
					a.equals("barcodes") || a.equals("expectedbarcodes")){
				expectedBarcodesFile=b;
			}else if(a.equals("delimiter")){
				delimiter=(b==null ? -1 : (byte)b.charAt(0));
				barcodesPerRead=(delimiter<0 ? 1 : 2);
			}else if(a.equals("out") || a.equals("counts")){
				parser.out1=b;
			}else if(a.equals("countsin")){
				countsIn=b;
			}else if(a.equals("mincount") || a.equals("mincounta") || a.equals("mincountf")){
				minCountA=Long.parseLong(b);
			}else if(a.equals("mincountr")){
				minCountR=Long.parseLong(b);
			}else if(a.equals("mincount0")){
				minCount0=Long.parseLong(b);
			}else if(a.equals("mincountpercent") || a.equals("mincountpercentile")){
				minCountPercentile=Float.parseFloat(b);
			}else if(a.equals("pcrmatrix") || a.equals("usematrix") || a.equals("matrix")){
				useMatrix=Parse.parseBoolean(b);
			}else if(a.equals("client") || a.equals("server") || a.equals("useserver")){
				useServer=Parse.parseBoolean(b);
			}
			
			else if(a.equals("probsout") || a.equals("outprobs")){
				probsOut=b;
			}else if(a.equals("accuracyout") || a.equals("outaccuracy")){
				accuracyOut=b;
			}else if(a.equals("mapout") || a.equals("outmap")){
				mapOut=b;
			}else if(a.equals("barcodesout") || a.equals("outbarcodes")){
				barcodesOut=b;
			}else if(a.equals("contamout") || a.equals("outcontam")){
				contamOut=b;
			}else if(a.equals("transitions") || a.equals("errors") || a.equals("substitutions")){
				transitionsOut=b;
			}else if(a.equals("matrixout") || a.equals("countmatrix") || a.equals("outmatrix")){
				matrixOut=b;
			}else if(a.equalsIgnoreCase("printStats")) {
				printStats=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("excludeUnknownInAccuracy")) {
				excludeUnknownInAccuracy=Parse.parseBoolean(b);
			}
			
			else if(a.equals("quantifysource")){
				quantifySource=Parse.parseBoolean(b);
			}else if(a.equals("quantifysink")){
				quantifySource=!Parse.parseBoolean(b);
			}else if(a.equals("quantset")){
				quantSetFile=b;
			}
			
			else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(PCRMatrix.parseStatic(arg, a, b)){
				//Flag was captured by PCRMatrix; do nothing
			}else if(parser.parse(arg, a, b)){
				//Flag was captured by the parser; do nothing
			}else if(b==null && new File(arg).exists()){
				in1.add(arg);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		return parser;
	}
	
	/** Replace # with 1 and 2 in headers */
	private void doPoundReplacement(){
		//Do input file # replacement
		if(in2.isEmpty()) {
			for(int i=0; i<in1.size(); i++) {
				String s=in1.get(i);
				if(s.indexOf('#')>-1 && !new File(s).exists()){
					String a=s.replace("#", "1");
					String b=s.replace("#", "2");
					in1.set(i, a);
					in2.add(b);
				}
			}
		}
		
		//Ensure there is an input file
		if(in1.isEmpty() && countsIn==null){
			throw new RuntimeException("Error - at least one input file is required.");
		}
	}
	
	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
	}
	
	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		ArrayList<String> temp=new ArrayList<String>();
		temp.addAll(in1);
		temp.addAll(in2);
		if(expectedBarcodesFile!=null) {temp.add(expectedBarcodesFile);}
		
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, temp.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		if(out1!=null) {temp.add(out1);}
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, temp.toArray(new String[0]))){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}
	
	/** Make sure interleaving agrees with number of input and output files */
	private void adjustInterleaving(){
		//Adjust interleaved detection based on the number of input files
		if(!in2.isEmpty()){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=false; //FASTQ.TEST_INTERLEAVED=false; interferes with barcode detection
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
		
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read stream
		if(countsIn==null) {
			for(int i=0; i<ffin1.size(); i++){
				processInner(ffin1.get(i), (ffin2.size()<=i ? null : ffin2.get(i)));
			}
		}else {
			processCounts(countsIn);
		}
		
		if(calcStats) {
			if(bs.barcodesPerRead>1) {
				bs.leftStats=bs.makeLeft();
				bs.rightStats=bs.makeRight();
				bs.leftStats.calcStats();
				bs.rightStats.calcStats();
			}
			bs.calcStats();
			
			if(printStats) {
				String s=bs.toStats("Barcodes");
				System.err.println(s);
				if(bs.leftStats!=null){
					System.err.println(bs.leftStats.toStats("\nBarcode1"));
				}
				if(bs.rightStats!=null){
					System.err.println(bs.rightStats.toStats("\nBarcode2"));
				}
			}
		}
		
		if(ffout1!=null){
			bs.printToFile(ffout1, minCount0);
		}
		
		if(transitionsOut!=null) {
			ByteBuilder bb=bs.printTransitions();
			ReadWrite.writeString(bb, transitionsOut);
		}
		
		if(matrixOut!=null || probsOut!=null || accuracyOut!=null || 
				mapOut!=null || barcodesOut!=null || useMatrix) {
			assert(!bs.expectedCodeList.isEmpty()) : 
				"For matrix output expected barcodes are required.";
			Timer t2=new Timer();

			final Collection<Barcode> counts=bs.codeMap.values();
			if(minCountPercentile>0) {
				long thresh=BarcodeCounter.barcodeCountPercentile(counts, minCountPercentile);
			}
			
			final HashMap<String, String> assignmentMap;
			{//New code block for using server

				if(PCRMatrix.matrixType0==PCRMatrix.PROB_TYPE && PCRMatrixProb.clientside()
						&& !setUseServer && !useServer) {
					useServer=true;
				}
				if(useServer) {
					System.err.println("Using client-server mode for barcode analysis.");
					DemuxData dd=new DemuxData(bs.length1, bs.length2, bs.delimiter);
					dd.codeCounts=counts;
					dd.expectedList=new LinkedHashSet<String>(bs.expectedCodeMap.keySet());
					DemuxClient client=new DemuxClient();
					assignmentMap=client.getMap(dd, verboseClient);
					t2.stop("Assignment Time:\t\t");
					if(mapOut!=null) {
						PCRMatrix.printAssignmentMapStatic(assignmentMap, 
								mapOut, counts, overwrite, append);
					}
				}else {
					PCRMatrix pcrMatrix=PCRMatrix.create(bs.length1, bs.length2, bs.delimiter);
					assert(bs.expectedCodeList!=null && !bs.expectedCodeList.isEmpty());
					pcrMatrix.populateExpected(bs.expectedCodeMap.keySet());
					pcrMatrix.populateSplitCodes();
					pcrMatrix.initializeData();
					pcrMatrix.refine(counts, minCountR);
					//			pcrMatrix.verbose=true;
					assignmentMap=pcrMatrix.makeAssignmentMap(counts, minCountA);
					t2.stop("Assignment Time:\t\t");
					if(mapOut!=null) {
						pcrMatrix.printAssignmentMap(assignmentMap, 
								mapOut, counts, overwrite, append);
					}
					
					if(contamOut==null) {//Otherwise yield will be printed later
//						long processed=0;
//						for(Barcode b : counts) {processed+=b.count();}
//						System.err.println("totalBarcodes="+totalBarcodes+
//								", pcrMatrix.totalAssignedToExpected="+pcrMatrix.totalAssignedToExpected);
						double yield=pcrMatrix.totalAssignedToExpected/(double)totalBarcodes;
						outstream.println(String.format("Total Yield:    \t\t%.4f", yield));
					}
					
					if(matrixOut!=null) {
						ByteBuilder bb=pcrMatrix.toBytes(null);
						ReadWrite.writeString(bb, matrixOut);
					}
					
					if(probsOut!=null) {
						ByteBuilder bb=pcrMatrix.toBytesProb(null);
						ReadWrite.writeString(bb, probsOut);
					}
					
					if(accuracyOut!=null) {
						ByteBuilder bb=pcrMatrix.toAccuracy(null, excludeUnknownInAccuracy);
						ReadWrite.writeString(bb, accuracyOut);
					}
					
//					pcrMatrix=null;
				}
			}
			
			if(barcodesOut!=null) {
				ByteStreamWriter bsw=new ByteStreamWriter(barcodesOut, overwrite, append, true);
				bsw.start();
				ArrayList<Barcode> list=Barcode.summateAssignments(assignmentMap, bs.expectedCodeList, bs.codeMap);
				for(Barcode b : list) {
					bsw.print(b.name).tab().print(b.frequency, 5).tab().print(b.count()).nl();
				}
				bsw.poisonAndWait();
			}
			
			//This block requires headers are tagged with correct barcode,
			//using rename.sh suffix="correctbarcode" for each Seal-binned library.
			if(contamOut!=null) {
				contamMap=new LinkedHashMap<String, LongPair>(2*bs.expectedCodeList.size());
				for(Barcode bc : bs.expectedCodeList) {
					if(bc.expected==1 && !bc.isHomopolymer('G')) {
						contamMap.put(bc.name, new LongPair());
					}
				}
				HashSet<String> quantSet=null;
				if(quantSetFile!=null) {
					quantSet=new HashSet<String>();
					for(String s : TextFile.toStringLines(quantSetFile)) {
						quantSet.add(s);
					}
				}
				
				IlluminaHeaderParser2.PARSE_COMMENT=true;
				long[] sum=new long[5];
				for(int i=0; i<ffin1.size(); i++){
					long[] ret=processContam(ffin1.get(i), assignmentMap, contamMap, quantSet, PCRMatrix.byTile);
					Tools.add(sum, ret);
				}
				long processed=sum[0];
				long numAssignedToBarcode=sum[1];
				long numAssignedToValid=sum[2];
				long correct=sum[3];
				long incorrect=sum[4];
				
				long good=0, bad=0;
				double logsum=0;
				double logsumD=0;
				double logsumN=0;
				for(Entry<String, LongPair> e : contamMap.entrySet()) {
					if(quantSet==null || quantSet.contains(e.getKey())) {
						LongPair lp=e.getValue();
						good+=lp.a;
						bad+=lp.b;
						logsumD+=Math.log(Tools.max(1.0, (lp.a+lp.b)));
						if(lp.b>0) {
							logsumN+=Math.log(lp.b);
							double ppm=lp.b*1000000.0/Tools.max(1.0, (lp.a+lp.b));
							logsum+=Math.log(ppm);
//							System.err.println("log="+Math.log(ppm)+", sum="+logsum);
						}
					}
				}
				assert(good==correct);
				final int denominator=(quantSet==null ? contamMap.size() : quantSet.size());
				double avgppm=bad*1000000.0/(double)(good+bad);
				double geoavg=Math.exp(logsum/denominator);
				double geoavgD=Math.exp(logsumD/denominator);
				double geoavgN=Math.exp(logsumN/denominator);
				double geoavg2=1000000*geoavgN/geoavgD;
				double yield=numAssignedToValid/(double)processed;
				double cyield=good/(double)processed;
				double iyield=bad/(double)processed;
				outstream.println(String.format("Contam Avg PPM: \t%.2f", avgppm));
				outstream.println(String.format("Contam Geo Avg: \t%.2f", geoavg));
				outstream.println(String.format("Contam Geo Avg2:\t%.2f", geoavg2));
				
				outstream.println(String.format("Total Yield:    \t%.4f", yield));
				outstream.println(String.format("Correct Yield:  \t%.4f", cyield));
				outstream.println(String.format("Incorrect Yield:\t%.7f", iyield));
//				System.err.println("geoavg="+geoavg+", logsum="+logsum+", contamMap.size()="+contamMap.size());
				ByteStreamWriter bsw=new ByteStreamWriter(contamOut, overwrite, append, true);
				bsw.start();
				bsw.print("#Total\t").println(good+bad);
				bsw.print("#Good\t").println(good);
				bsw.print("#Bad\t").println(bad);
				bsw.print("#AvgPPM\t").println(avgppm, 2);
				bsw.print("#GeoPPM\t").println(geoavg, 2);
				bsw.print("#TotalYield\t").println(yield, 5);
				bsw.print("#CorrectYield\t").println(cyield, 5);
				bsw.println("#Code\tGood\tBad\tPPM");
				for(Entry<String, LongPair> e : contamMap.entrySet()) {
					String key=e.getKey();
					LongPair lp=e.getValue();
					double ppm=lp.b*1000000.0/Tools.max(1.0, (lp.a+lp.b));
					bsw.print(key).tab().print(lp.a).tab().print(lp.b).tab().println(ppm, 2);
				}
				bsw.poisonAndWait();
			}
		}
		
		if(verbose){outstream.println("Finished; closing streams.");}
		
		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		
		//Report timing and results
		t.stop();
		outstream.println();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private ConcurrentReadInputStream makeCris(FileFormat ff1, FileFormat ff2, boolean printPaired){
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, ff2, 
				null, null);
		cris.start(); //Start the stream
		if(verbose){outstream.println("Started cris");}
		boolean paired=cris.paired();
		if(printPaired && !ff1.samOrBam()){
			outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));
		}
		return cris;
	}
	
	/** Iterate through the reads */
	void processCounts(final String countsFile){
		
		ByteFile bf=ByteFile.makeByteFile(countsFile, true);
		LineParser lp=new LineParser1('\t');
		
		byte[] line=bf.nextLine();
		while(line!=null && Tools.startsWith(line, '#')) {line=bf.nextLine();}
		assert(line!=null) : "Empty counts file.";
		{
			lp.set(line);
			byte[] code=lp.parseByteArray(0);
			delimiter=Barcode.delimiter(code);
			barcodesPerRead=delimiter>0 ? 2 : 1;
			
			//TODO: This does not take tile into account.  EDIT:  Should be fine now
			final int bcLen1, bcLen2;
			if(barcodesPerRead==1) {
				bcLen1=code.length-Tools.trailingDigits(code);
				bcLen2=0;
			}else {
				bcLen1=Tools.indexOf(code, delimiter);
				bcLen2=code.length-bcLen1-1-Tools.trailingDigits(code);
			}
//			assert(false) : bcLen1+", "+bcLen2+", "+PCRMatrix.byTile;
			assert(bs==null); //Otherwise why call this, unless you have multiple counts files
			if(bs==null) {
				bs=new BarcodeStats(delimiter, barcodesPerRead, extin);
				bs.length1=bcLen1;
				bs.length2=bcLen2;
				if(expectedBarcodesFile!=null) {
					bs.loadBarcodeList(expectedBarcodesFile, delimiter, false, false);
				}
			}
		}
		
		for(; line!=null; line=bf.nextLine()) {
			lp.set(line);
			final String prefix=lp.parseString(0);
			final String code;
			final int count=lp.parseInt(1);
			final int tile;
			if(PCRMatrix.byTile) {
				int digits=Tools.trailingDigits(prefix);
				assert(digits>0);
				int idx=prefix.length()-digits;
				code=new String(line, 0, idx);
				tile=Parse.parseInt(prefix, idx);
				assert(tile>0);
//				assert(false) : digits+", "+idx+", "+line.length+", '"+new String(code)+"'";
			}else {
				code=prefix;
				tile=0;
			}
//			assert(false) : new String(code)+", "+tile+", "+count+", "+new String(line);
			totalBarcodes+=count;
			if(count>=minCount0) {
				bs.increment(code, count, tile);
			}
		}
		
		
		//Close the read streams
		errorState|=bf.close();
		
		//Do anything necessary after processing
		
	}
	
	/** Iterate through the reads */
	void processInner(final FileFormat ff1, final FileFormat ff2){
		
		//Do anything necessary prior to processing
		
		//Create a read input stream
		final ConcurrentReadInputStream cris=makeCris(ff1, ff2, true);
		
		if(bs==null) {
			if(delimiter<0) {
				delimiter=(byte)ff1.barcodeDelimiter();
				barcodesPerRead=ff1.barcodesPerRead();
			}
			bs=new BarcodeStats(delimiter, barcodesPerRead, extin);
			if(expectedBarcodesFile!=null) {
				bs.loadBarcodeList(expectedBarcodesFile, delimiter, false, false);
			}
		}
		if(bs.length1<1 && bs.length2<1) {
			bs.length1=ff1.barcodeLength(1);
			bs.length2=ff1.barcodeLength(2);
		}
		readFile(cris);
		
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
		
		//Do anything necessary after processing
		
	}
	
	/** Iterate through the reads */
	long[] processContam(FileFormat ff, HashMap<String, String> assignmentMap, 
			HashMap<String, LongPair> counts, HashSet<String> quantSet, boolean addTile){
		
		//Do anything necessary prior to processing
		
		//Create a read input stream
		final ConcurrentReadInputStream cris=makeCris(ff, null, false);
		
		//Grab the first ListNum of reads
		ListNum<Read> ln=cris.nextList();
		
		long processed=0, numAssignedToBarcode=0, numAssignedToValid=0, correct=0, incorrect=0;

		//As long as there is a nonempty read list...
		while(ln!=null && ln.size()>0){
//			if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
			
			for(Read r : ln) {
				processed++;
				ihp.parse(r);
				final String barcode=ihp.barcode();
				assert(barcode!=null) : r.id;
//				totalBarcodes++; //This should already be handled when reading counts or the reads
				final int tile=ihp.tile();
				String key=(addTile ? barcode+tile : barcode);
				String assigned=assignmentMap.get(key);
				final String actual=ihp.extra();
				if(assigned!=null) {
					numAssignedToBarcode++;
					if(counts.containsKey(assigned)) {numAssignedToValid++;}
				}
				if(actual!=null) {
					if(quantifySource) {
						LongPair lp=counts.get(actual);
						if(lp==null) {
							assert(false) : "Missing key "+actual+" in:\t"+counts.keySet();
							lp=new LongPair();
							counts.put(actual, lp);
						}
						if(lp!=null) {
							if(assigned==null || (quantSet!=null && !quantSet.contains(actual))) {}//Do nothing
							else if(assigned.equals(actual)) {lp.a++;correct++;}
							else {lp.b++;incorrect++;}
						}
					}else {//Quantify sink
						LongPair lp=(assigned==null ? null : counts.get(assigned));
						if(assigned==null) {}//do nothing - it was not assigned to anything
						else if(actual==null) {} //do nothing - we don't know where it was supposed to go
						else if((quantSet!=null && !quantSet.contains(assigned))) {} //do nothing 
						else if(lp==null) {//do nothing - it was assigned to a dummy barcode like GGGGGGGG
							assert(Barcode.isHomopolymer(assigned)) : 
								r.id+" -> "+assigned+"\n"+counts.keySet()+"\n";
						}else if(quantSet==null || (quantSet.contains(assigned))){
							if(assigned.equals(actual)) {lp.a++;correct++;}
							else {lp.b++;incorrect++;}
						}
					}
				}
			}
			
			cris.returnList(ln);
			
			//Fetch a new list
			ln=cris.nextList();
		}
		
		//Notify the input stream that the final list was used
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
		
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);
		
		//Do anything necessary after processing
		return new long[] {processed, numAssignedToBarcode, numAssignedToValid, correct, incorrect};
	}
	
	void readFile(final ConcurrentReadInputStream cris) {
		{
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();

			//Check to ensure pairing is as expected
			if(ln!=null && !ln.isEmpty()){
				Read r=ln.get(0);
				assert(r.samline!=null || (r.mate!=null)==cris.paired());
			}

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
				
				processList(ln, cris);

				//Fetch a new list
				ln=cris.nextList();
			}
			
			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
	}
	
	/**
	 * Process a list of Reads.
	 * @param ln The list.
	 * @param cris Read Input Stream
	 * @param ros Read Output Stream for reads that will be retained
	 */
	void processList(ListNum<Read> ln, final ConcurrentReadInputStream cris){

		//Grab the actual read list from the ListNum
		final ArrayList<Read> reads=ln.list;
		
		//Loop through each read in the list
		for(int idx=0; idx<reads.size(); idx++){
			final Read r1=reads.get(idx);
			final Read r2=r1.mate;
			
			//Validate reads in worker threads
			if(!r1.validated()){r1.validate(true);}
			if(r2!=null && !r2.validated()){r2.validate(true);}

			//Track the initial length for statistics
			final int initialLength1=r1.length();
			final int initialLength2=r1.mateLength();

			//Increment counters
			readsProcessed+=r1.pairCount();
			basesProcessed+=initialLength1+initialLength2;
			
			{
				//Reads are processed in this block.
				processReadPair(r1, r2);
			}
		}

		//Notify the input stream that the list was used
		cris.returnList(ln);
//		if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access
	}
	
	
	/**
	 * Process a single read pair.
	 * @param r1 Read 1
	 * @param r2 Read 2 (may be null)
	 * @return True if the reads should be kept, false if they should be discarded.
	 */
	void processReadPair(final Read r1, final Read r2){
		ihp.parse(r1);
		String code=ihp.barcode();
		if(code==null){return;}
		totalBarcodes++;
		int tile=0;
		if(PCRMatrix.byTile) {
			tile=ihp.tile();
			assert(tile>1000 && tile<3000);//For Illumina
		}
		bs.increment(code, 1, tile);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private ArrayList<String> in1=new ArrayList<String>();
	/** Secondary input file path */
	private ArrayList<String> in2=new ArrayList<String>();
	
	/** Expected Barcodes */
	private String expectedBarcodesFile=null;
	
	/** Pre-counted barcodes */
	private String countsIn=null;
	private long minCount0=0;
	private long minCountR=4;
	private long minCountA=4;
	private float minCountPercentile=0.0f;

	/** Primary output file path */
	private String out1=null;
	
	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;

	private String transitionsOut=null;
	private String matrixOut=null;
	private String probsOut=null;
	private String accuracyOut=null;
	private String mapOut=null;
	private String barcodesOut=null;
	private String contamOut=null;
	private String quantSetFile=null;
	private boolean quantifySource=false;
	private boolean excludeUnknownInAccuracy=true;
	
	private boolean useServer=false;
	private boolean setUseServer=false;
	
	/** Whether interleaved was explicitly set. */
	private boolean setInterleaved=false;

	private HashMap<String, LongPair> contamMap;
	
	private BarcodeStats bs;
	private boolean calcStats=true;
	private boolean printStats=true;
	
	byte delimiter=-1;
	int barcodesPerRead=-1;
	
//	private PCRMatrix pcrMatrix;
	private boolean useMatrix=false;
	
//	private LineParserS2 lp=new LineParserS2(':');
	private IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();
	
	/*--------------------------------------------------------------*/
	
	long totalBarcodes=0;

	/** Number of reads processed */
	protected long readsProcessed=0;
	/** Number of bases processed */
	protected long basesProcessed=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final ArrayList<FileFormat> ffin1=new ArrayList<FileFormat>();
	/** Secondary input file */
	private final ArrayList<FileFormat> ffin2=new ArrayList<FileFormat>();
	
	/** Primary output file */
	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	public static boolean verboseClient=true;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	
}
