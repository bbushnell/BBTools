package barcode;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashSet;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import parse.LineParserS4;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.TrimRead;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Processes demultiplexed sequencing files by tagging reads with the barcode parsed
 * from each input filename.
 * Files are processed in argument order and merged into {@code out}; the filename
 * extension selects any supported output format, with FASTQ as the fallback. For
 * example, reads from "52866.4.475040.GAGGCCGCCA-TTATCTAGCT.filter-DNA.fastq.gz"
 * receive the suffix "(tab)GAGGCCGCCA+TTATCTAGCT".
 *
 * @author Brian Bushnell
 * @date May 10, 2024
 */
public class TagAndMerge {

/**
	 * Program entry point for barcode tagging and merging operations.
	 * Creates a TagAndMerge instance, processes all input files, and handles timing.
	 * @param args Command-line arguments specifying input files and options
	 */
	public static void main(String[] args){
//Start a timer immediately upon code entrance.
		Timer t=new Timer();

//Create an instance of this class
		TagAndMerge x=new TagAndMerge(args);

//Run the object
		x.process(t);

//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

/**
	 * Constructs a TagAndMerge instance and parses command-line arguments.
	 * Initializes the input file list, output settings, trimming options, and
	 * paired character remapping for parsed barcodes. {@code in=} clears any
	 * earlier inputs and either accepts one existing literal path or splits a
	 * comma-delimited list; bare existing paths append to the current list.
	 * {@code reads}/{@code maxreads} applies separately to each input stream.
	 * @param args Command-line arguments including input files, output path, and processing options
	 */
	public TagAndMerge(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
//Set a variable here
			}else if(a.equals("trim")){
				trimLen=(b==null ? 1 : Integer.parseInt(b));
			}else if(a.equals("r1only") || a.equals("read1only") || a.equals("dropr2")){
				dropR2=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("shrinkHeader") || a.equals("shrink")){
				shrinkHeader=Parse.parseBoolean(b);
			}else if(a.equals("barcodes") || a.equals("barcodesout")){
				barcodesOutFile=b;
			}else if(a.equals("remap")){
				remap=(b==null ? null : b.toCharArray());
			}else if(a.equals("in")){
				in.clear();
				String[] b2=(b==null) ? null : (new File(b).exists() ? new String[] {b} : b.split(","));
				for(String b3 : b2){in.add(b3);}
			}else if(b==null && new File(arg).exists()){
				in.add(arg);
			}else if(parser.parse(arg, a, b)){
//do nothing
			}else{
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
		}

		{//Process parser fields
			Parser.processQuality();

			maxReads=parser.maxReads;
			out1=parser.out1;
			overwrite=parser.overwrite;
			append=parser.append;
		}

//Generic output construction does not reject existing ordinary files.
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			throw new RuntimeException("\n\noverwrite="+overwrite+
"; Can't write to output file "+out1+"\n");
		}

		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, null, true, overwrite, append, false);
	}

/**
	 * Creates the optional merged output stream, processes each input file with barcode
	 * tagging, writes distinct observed barcodes to {@code barcodesOutFile} if requested,
	 * reports processing statistics, and throws if a read or write stream reported an
	 * error during close.
	 * @param t Timer for tracking total execution time
	 */
	void process(Timer t){

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=4;

			for(String s : in){
				assert(!out1.equalsIgnoreCase(s)) : "Input file and output file have same name.";
			}

			ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, false);
			ros.start();
		}else{ros=null;}

		for(String s : in){
			processInner(s, ros);
		}

		errorState|=ReadWrite.closeStream(ros);
		outstream.println("Finished processing "+in.size()+" files.");

		if(barcodesOutFile!=null){
			TextStreamWriter tsw=new TextStreamWriter(barcodesOutFile, overwrite, append, false);
			tsw.start();
			for(String tag : barcodes){tsw.println(tag);}
			errorState|=tsw.poisonAndWait();
		}

		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

/**
	 * Processes a single input file by extracting its barcode and tagging all reads.
	 * Parses the barcode from {@code fname}, applies any character remapping,
	 * consumes reads up to the per-input {@code reads}/{@code maxreads} limit,
	 * tags each header, and writes the tagged lists when an output stream was
	 * requested.
	 *
	 * @param fname Input filename containing barcode information
	 * @param ros Output stream for writing tagged reads, or null to consume input only
	 */
	void processInner(String fname, ConcurrentReadOutputStream ros){
		FileFormat ffin=FileFormat.testInput(fname, FileFormat.FASTQ, null, true, true);
		String tag=Barcode.parseBarcodeFromFname(fname);
		assert(tag!=null) : "Can't find barcode in filename "+fname;//Possible bug: assertion disabled in production
		if(remap!=null){
			for(int i=0; i<remap.length; i+=2){
				tag=tag.replace(remap[i], remap[i+1]);
			}
		}
		if(tag!=null && Barcode.isBarcode(tag)){
			barcodes.add(tag);
		}

		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin, null);
			cris.start();
		}
		boolean paired=cris.paired();

		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin==null || ffin.samOrBam()) || (r.mate!=null)==cris.paired());
			}

//ln!=null prevents a compiler potential null access warning
			while(ln!=null && reads!=null && reads.size()>0){
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}

				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					processPair(r1, tag);
				}

				if(ros!=null){ros.add(reads, ln.id);}

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		errorState|=ReadWrite.closeStream(cris);
	}

/**
	 * Accounts one read pair, optionally discards R2, tags the remaining reads,
	 * and updates output counters after the optional mate removal.
	 *
	 * @param r Non-null forward read (R1)
	 * @param tag Barcode string to append to read headers
	 */
	private void processPair(Read r, String tag){
		readsProcessed+=r.pairCount();
		basesProcessed+=r.pairLength();
		if(dropR2){r.mate=null;}
		processRead(r, tag);
		processRead(r.mate, tag);
		readsOut+=r.pairCount();
		basesOut+=r.pairLength();
	}

/**
	 * Optionally trims a read, then replaces its header with the full original
	 * header or compact Illumina fields followed by a tab and barcode.
	 *
	 * @param r Read to process, or null for an absent/discarded mate
	 * @param tag Barcode string to append to the read header
	 */
	private void processRead(Read r, String tag){
		if(r==null){return;}
		if(trimLen>=0 && trimLen<r.length()){
			final int x=TrimRead.trimToPosition(r, 0, trimLen-1, 1);
//COMPREHENSION (HK416 2026-07-02): the guard above (trimLen<r.length()) already EXCLUDES trimLen>=length,
//so the old note's scenario can't occur. The real boundary is trimLen==0 (settable via trim=0): then
//trimToPosition(r,0,-1,1) trims to an empty read and x is likely <=0 -> this assert fires. trim=0 is a
//nonsensical flag value (trim to zero length), so this is a crash-loud on bad input, not a live bug.
			assert(x>0) : x+", "+trimLen+", "+r.length();
		}
		bb.clear();
		if(!shrinkHeader){
			bb.append(r.id).tab().append(tag);
		}else{
//LH00223:36:22JWNJLT3:8:1101:31031:1064 1:N:0:AAAGACCAGG+GAGTTTGAGA
//			ihp.parse(r.id);
//			final LineParser lp=ihp.lp();
			lp.set(r.id);
			bb.colon().colon().colon().appendTerm(lp, 3).colon().appendTerm(lp, 4).colon();
			bb.space().appendTerm(lp, 7);
//			bb.colon().appendTerm(lp,  8).colon().appendTerm(lp, 9).colon().appendTerm(lp,  10);
			bb.tab().append(tag);
		}
		r.id=bb.toString();
	}

/*--------------------------------------------------------------*/
/*----------------            Fields            ----------------*/
/*--------------------------------------------------------------*/

/** Input files in processing order. */
	private ArrayList<String> in=new ArrayList<String>();
/** Merged read output path, or null to consume inputs without writing reads. */
	private String out1=null;
/** Optional table of distinct observed barcode strings. */
	private String barcodesOutFile=null;
/** Distinct valid barcodes parsed from input filenames, in first-seen order. */
	private LinkedHashSet<String> barcodes=new LinkedHashSet<String>();

/** If nonnegative, trim each retained read to at most this many bases. */
	private int trimLen=-1;
/** If true, discard R2 after accounting the input pair. */
	private boolean dropR2=false;
/** If true, shrink Illumina headers before appending a tab and barcode. */
	private boolean shrinkHeader=false;
/** Alternating old/new characters applied to parsed barcodes before tagging. */
	private char[] remap=new char[] {'-','+'};

/** Merged read output format and overwrite/append policy. */
	private final FileFormat ffout1;
/** Whether output creation may replace existing files. */
	private boolean overwrite=true;
/** Whether output creation appends to existing raw/gzip files. */
	private boolean append=false;
/** Latched when any read/output close reports a corrupt or incomplete stream. */
	private boolean errorState=false;

/*--------------------------------------------------------------*/
/*----------------           Counters           ----------------*/
/*--------------------------------------------------------------*/

/** Per-input fragment limit passed to ConcurrentReadInputStream, or -1 for no limit. */
	private long maxReads=-1;
/** Input reads and bases seen before optional mate removal or trimming. */
	private long readsProcessed=0, basesProcessed=0;
/** Reads and bases retained after optional mate removal and trimming. */
	private long readsOut=0, basesOut=0;
/** Reusable builder for per-read header replacement. */
	private final ByteBuilder bb=new ByteBuilder();
//	private final IlluminaHeaderParser2 ihp=new IlluminaHeaderParser2();
/** Reusable parser for compact Illumina header fields. */
	private final LineParserS4 lp=new LineParserS4(":::::: ");

/*--------------------------------------------------------------*/
/*----------------             State            ----------------*/
/*--------------------------------------------------------------*/

/** Destination for progress and summary messages. */
	private java.io.PrintStream outstream=System.err;
/** Whether to print per-list fetch counts while processing. */
	public static boolean verbose=false;

}
