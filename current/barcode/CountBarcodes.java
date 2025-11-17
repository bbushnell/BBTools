package barcode;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import align2.BandedAlignerConcrete;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import structures.StringNum;
import tracker.ReadStats;

/**
 * @author Brian Bushnell
 * @date June 20, 2014
 *
 */
public class CountBarcodes {

	/** Program entry point.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		CountBarcodes x=new CountBarcodes(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructor that parses command-line arguments and initializes file formats.
	 * Sets up input/output streams, barcode maps, and processing parameters.
	 * @param args Command-line arguments for configuration
	 */
	public CountBarcodes(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean setInterleaved=false; //Whether it was explicitly set.
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		

		expectedCodeMap=new HashMap<String,String>(200);
		validCodeMap=new HashMap<String,String>(200);
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("countundefined")){
				countUndefined=Parse.parseBoolean(b);
			}else if(a.equals("printheader") || a.equals("header")){
				printHeader=Parse.parseBoolean(b);
			}else if(a.equals("printrows") || a.equals("rows") || a.equals("maxrows")){
				maxRows=Integer.parseInt(b);
			}else if(a.equals("expected")){
				if(b!=null){
					for(String code : b.split(",")){
						expectedCodeMap.put(code, code);
						validCodeMap.put(code, code);
					}
				}
			}else if(a.equals("valid")){
				if(b!=null){
					for(String code : b.split(",")){
						validCodeMap.put(code, code);
					}
				}
			}else if(a.equals("counts") || a.equals("outc")){
				outCounts=b;
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
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

			out1=parser.out1;
			out2=parser.out2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
		
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
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+"\n");
		}
		
		validCodes=new ArrayList<String>(validCodeMap.size());
		validCodes.addAll(validCodeMap.values());
		expectedCodes=new ArrayList<String>(expectedCodeMap.size());
		expectedCodes.addAll(expectedCodeMap.values());
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffout2=FileFormat.testOutput(out2, FileFormat.FASTQ, extout, true, overwrite, append, false);
		ffCounts=FileFormat.testOutput(outCounts, FileFormat.TEXT, null, false, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}
	
	/**
	 * Main processing method that reads sequences, counts barcodes, and generates output.
	 * Extracts barcodes from read headers, maintains counts in a hash map, calculates
	 * distances to expected codes, and writes results to output files.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, qfin1, qfin2);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
//		if(verbose){
			if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
//		}

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=4;
			
			if(cris.paired() && out2==null && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) : "out1 and out2 have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, ffout2, qfout1, qfout2, buff, null, false);
			ros.start();
		}else{ros=null;}
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		HashMap<String, StringNum> map=new HashMap<String, StringNum>();
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					
					final int initialLength1=r1.length();
					final int initialLength2=(r1.mateLength());
					
					{
						readsProcessed++;
						basesProcessed+=initialLength1;
					}
					if(r2!=null){
						readsProcessed++;
						basesProcessed+=initialLength2;
					}
					
//					String id=r1.id;
//					int colon=id.lastIndexOf(':');
//					String barcode=id.substring(colon+1);
					final String barcode=r1.barcode(true);
					
					if(countUndefined || AminoAcid.isFullyDefined(barcode)){
						StringNum value=map.get(barcode);
						if(value==null){
							value=new StringNum(barcode, 0);
							map.put(barcode, value);
						}
						value.increment();
					}
				}
				
				final ArrayList<Read> listOut=reads;
				
				if(ros!=null){ros.add(listOut, ln.id);}

				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		ArrayList<StringNum> list=new ArrayList<StringNum>(map.size());
		list.addAll(map.values());
		Shared.sort(list);
		Collections.reverse(list);
		
		TextStreamWriter tsw=new TextStreamWriter(ffCounts);
		tsw.start();
		if(printHeader){
			tsw.print("#code\tcount\tHamming_dist\tedit_dist\tvalid\n");
		}
		for(StringNum sn : list){
			if(maxRows==0){break;}
			maxRows--;
			int hdist=calcHdist(sn.s, expectedCodes);
			int edist=hdist;
			if(hdist>1){
				try {
					edist=calcEdist(sn.s, expectedCodes);
				} catch (Exception e) {
					edist=hdist;
				}
			}
			boolean valid=validCodes.contains(sn.s);
			tsw.print(sn+"\t"+hdist+"\t"+edist+"\t"+(valid ? "valid" : "")+"\n");
		}
		tsw.poisonAndWait();
		errorState|=tsw.errorState;
		
		errorState|=ReadStats.writeAll();
		
		errorState|=ReadWrite.closeStreams(cris, ros);
		
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/**
	 * @param s
	 * @param expectedCodes
	 * @return
	 */
	public static int calcHdist(String s, ArrayList<String> expectedCodes) {
		int min=s.length();
		for(String code : expectedCodes){
			min=Tools.min(min, hdist(s, code));
			if(min<1){break;}
		}
		return min;
	}
	
	/**
	 * Calculates Hamming distance between two strings by counting substitutions.
	 * Only counts mismatches up to the length of the shorter string.
	 *
	 * @param s First string to compare
	 * @param code Second string to compare
	 * @return Number of character substitutions needed
	 */
	static int hdist(String s, String code) {
		final int min=Tools.min(s.length(), code.length());
		int subs=0;
		for(int i=0; i<min; i++){
			if(s.charAt(i)!=code.charAt(i)){subs++;}
		}
		return subs;
	}
	
	/**
	 * @param s
	 * @param expectedCodes
	 * @return
	 */
	int calcEdist(String s, ArrayList<String> expectedCodes) {
		int min=s.length();
		for(String code : expectedCodes){
			min=Tools.min(min, edist(s, code));
			if(min<1){break;}
		}
		return min;
	}
	
	/**
	 * Calculates edit distance between two strings using banded alignment.
	 * Allows insertions, deletions, and substitutions in the alignment.
	 *
	 * @param s First string to align
	 * @param code Second string to align
	 * @return Edit distance between the strings
	 */
	int edist(String s, String code) {
		int x=bandy.alignForward(s.getBytes(), code.getBytes(), 0, 0, s.length(), true);
		return x;
	}
	
	/*--------------------------------------------------------------*/
	
	/** First input file path */
	private String in1=null;
	/** Second input file path for paired reads */
	private String in2=null;
	
	/** Quality file for first input */
	private String qfin1=null;
	/** Quality file for second input */
	private String qfin2=null;

	/** First output file path */
	private String out1=null;
	/** Second output file path for paired reads */
	private String out2=null;
	/** Output file path for barcode count results */
	private String outCounts=null;

	/** Quality output file for first stream */
	private String qfout1=null;
	/** Quality output file for second stream */
	private String qfout2=null;
	
	/** Input file extension override */
	private String extin=null;
	/** Output file extension override */
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	/** Whether to reverse complement the mate read */
	private boolean reverseComplementMate=false;
	/** Whether to reverse complement reads */
	private boolean reverseComplement=false;
	/** Whether to count barcodes containing undefined bases (default true) */
	private boolean countUndefined=true;
	/** Whether to print column headers in output (default true) */
	private boolean printHeader=true;
	/** Maximum number of barcode rows to output (-1 for unlimited) */
	private int maxRows=-1;

	/** Maximum number of reads to process (-1 for unlimited) */
	private long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	
	/** File format handler for first input file */
	private final FileFormat ffin1;
	/** File format handler for second input file */
	private final FileFormat ffin2;

	/** File format handler for first output file */
	private final FileFormat ffout1;
	/** File format handler for second output file */
	private final FileFormat ffout2;
	/** File format handler for barcode count output */
	private final FileFormat ffCounts;

	/** Map of expected barcode codes for distance calculations */
	private final HashMap<String,String> expectedCodeMap;
	/** Map of valid barcode codes for validation marking */
	private final HashMap<String,String> validCodeMap;
	/** List of expected barcode sequences for distance calculations */
	private final ArrayList<String> expectedCodes;
	/** List of valid barcode sequences for validation */
	private final ArrayList<String> validCodes;
	
	/** Banded aligner for calculating edit distances with k=21 */
	private final BandedAlignerConcrete bandy=new BandedAlignerConcrete(21);
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for progress and error messages (default System.err) */
	private PrintStream outstream=System.err;
	/** Enable verbose output for debugging */
	public static boolean verbose=false;
	/** Flag indicating if an error occurred during processing */
	public boolean errorState=false;
	/** Whether to overwrite existing output files (default true) */
	private boolean overwrite=true;
	/** Whether to append to existing output files (default false) */
	private boolean append=false;
	
}
