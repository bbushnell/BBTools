package jgi;

import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.MetadataWriter;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sort.ReadLengthComparator;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Reformats a fungal assembly for release. Also creates contig and agp files.
 * 
 * @author Brian Bushnell
 * @date December 9, 2015
 *
 */
public class FungalRelease {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Code entrance from the command line.
	 * 
	 * @param args
	 *            Command line arguments
	 */
	public static void main(String[] args) {
		Timer t = new Timer();
		FungalRelease x = new FungalRelease(args);
		x.process(t);

		// Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructor.
	 * 
	 * @param args
	 *            Command line arguments
	 */
	public FungalRelease(String[] args) {

		{// Preparse block for help, config files, and outstream
			PreParser pp = new PreParser(args, getClass(), false);
			args = pp.args;
			outstream = pp.outstream;
		}

		FASTQ.FORCE_INTERLEAVED = FASTQ.TEST_INTERLEAVED = false;
		Shared.FASTA_WRAP = 60;

		// Set shared static variables
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ = ReadWrite.USE_UNPIGZ = true;
		ReadWrite.setZipThreads(Shared.threads());

		Read.TO_UPPER_CASE = true;

		// Create a parser object
		Parser parser = new Parser();

		// Parse each argument
		for (int i = 0; i < args.length; i++) {
			String arg = args[i];

			// Break arguments into their constituent parts, in the form of
			// "a=b"
			String[] split = arg.split("=");
			String a = split[0].toLowerCase();
			String b = split.length > 1 ? split[1] : null;

			if (parser.parse(arg, a, b)) {// Parse standard flags in the parser
				// do nothing
			} else if (a.equals("verbose")) {
				verbose = Parse.parseBoolean(b);
			} else if (a.equals("mingapin")) {
				minGapIn = (int) Parse.parseKMG(b);
			} else if (a.equals("mingap") || a.equals("mingapout")) {
				minGapOut = (int) Parse.parseKMG(b);
			} else if (a.equals("minlen") || a.equals("minlength") || a.equals("minscaf")) {
				minScaf = (int) Parse.parseKMG(b);
			} else if (a.equals("mincontig")) {
				minContig = (int) Parse.parseKMG(b);
			} else if (a.equals("outc") || a.equals("contigs")) {
				outC = b;
			} else if (a.equals("qfoutc")) {
				qfoutC = b;
			} else if (a.equals("sortcontigs")) {
				sortContigs = Parse.parseBoolean(b);
			} else if (a.equals("sortcscaffolds")) {
				sortScaffolds = Parse.parseBoolean(b);
			} else if (a.equals("baniupac")) {
				banIupac = Parse.parseBoolean(b);
			} else if (a.equals("agp")) {
				agpFile = b;
			} else if (a.equals("legend")) {
				legendFile = b;
			} else if (a.equals("scafnum")) {
				scafNum = Parse.parseKMG(b);
			} else if (a.equals("renamescaffolds") || a.equals("rename")) {
				renameScaffolds = Parse.parseBoolean(b);
			} else if (a.equals("scafnum")) {
				contigNum = Parse.parseKMG(b);
			} else if (a.equals("renamecontigs")) {
				renameContigs = Parse.parseBoolean(b);
			} else if (a.equals("parse_flag_goes_here")) {
				// Set a variable here
			} else {
				outstream.println("Unknown parameter " + args[i]);
				assert (false) : "Unknown parameter " + args[i];
				// throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}

		{// Process parser fields
			Parser.processQuality();

			maxReads = parser.maxReads;

			overwrite = ReadStats.overwrite = parser.overwrite;
			append = ReadStats.append = parser.append;

			in1 = parser.in1;
			qfin1 = parser.qfin1;

			out1 = parser.out1;
			qfout1 = parser.qfout1;

			extin = parser.extin;
			extout = parser.extout;
		}

		assert (FastaReadInputStream.settingsOK());

		// Ensure there is an input file
		if (in1 == null) {
			throw new RuntimeException("Error - at least one input file is required.");
		}

		// Adjust the number of threads for input file reading
		if (!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads() > 2) {
			ByteFile.FORCE_MODE_BF2 = true;
		}

		// Ensure output files can be written
		if (!Tools.testOutputFiles(overwrite, append, false, out1, outC)) {
			outstream.println((out1 == null) + ", " + out1);
			throw new RuntimeException("\n\noverwrite=" + overwrite + "; Can't write to output files " + out1 + "\n");
		}

		// Ensure input files can be read
		if (!Tools.testInputFiles(false, true, in1)) {
			throw new RuntimeException("\nCan't read some input files.\n");  
		}

		// Ensure that no file was specified multiple times
		if (!Tools.testForDuplicateFiles(true, in1, out1, outC)) {
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}

		// Create output FileFormat objects
		ffout1 = FileFormat.testOutput(out1, FileFormat.FASTA, extout, true, overwrite, append, false);

		// Create output FileFormat objects
		ffoutC = FileFormat.testOutput(outC, FileFormat.FASTA, extout, true, overwrite, append, false);

		// Create input FileFormat objects
		ffin1 = FileFormat.testInput(in1, FileFormat.FASTA, extin, true, true);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t) {

		// Create a read input stream
		final ConcurrentReadInputStream cris;
		{
			cris = ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null, qfin1, null);
			cris.start(); // Start the stream
			if (verbose) {
				outstream.println("Started cris");
			}
		}

		// Optionally create a read output stream
		final ConcurrentReadOutputStream ros, rosc;
		final int buff = 4;
		if (ffout1 != null) {
			ros = ConcurrentReadOutputStream.getStream(ffout1, null, qfout1, null, buff, null, false);
			ros.start(); // Start the stream
		} else {
			ros = null;
		}
		if (ffoutC != null) {
			rosc = ConcurrentReadOutputStream.getStream(ffoutC, null, qfoutC, null, 4, null, false);
			rosc.start(); // Start the stream
		} else {
			rosc = null;
		}

		// Reset counters
		readsProcessed = 0;
		basesProcessed = 0;
		readsOut = 0;
		basesOut = 0;

		// Process the read stream
		processInner(cris, ros, rosc);

		if (verbose) {
			outstream.println("Finished; closing streams.");
		}

		// Write anything that was accumulated by ReadStats
		errorState |= ReadStats.writeAll();
		// Close the read streams
		errorState |= ReadWrite.closeStreams(cris, ros, rosc);

		// Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		MetadataWriter.write(null, readsProcessed, basesProcessed, readsOut, basesOut, false);
		
		// Throw an exception of there was an error in a thread
		if (errorState) {
			throw new RuntimeException(
					getClass().getName() + " terminated in an error state; the output may be corrupt.");
		}
	}

	/** Iterate through the reads */
	void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros,
			final ConcurrentReadOutputStream rosc) {

		ArrayList<Read> scaffolds = getReads(cris);

		final boolean makeLegend = (legendFile != null);
		ByteStreamWriter bswLegend = (makeLegend ? new ByteStreamWriter(legendFile, overwrite, append, false) : null);
		if (bswLegend != null) {
			bswLegend.start();
		}
		
		if (ros != null) {
			if (sortScaffolds) {
				Shared.sort(scaffolds, ReadLengthComparator.comparator);
			}
			if (renameScaffolds) {
				for (Read r : scaffolds) {
					String old = r.id;
					r.id = "scaffold_" + scafNum;
					if (bswLegend != null) {
						bswLegend.println(old + "\t" + r.id);
					}
					scafNum++;
				}
			}
			ros.add(scaffolds, 0);
		}
		if (bswLegend != null) {
			bswLegend.poisonAndWait();
		}

		final boolean makeAgp = (agpFile != null);
		ByteStreamWriter bswAgp = (makeAgp ? new ByteStreamWriter(agpFile, overwrite, append, false) : null);
		if (bswAgp != null) {
			bswAgp.start();
		}

		if (rosc != null || makeAgp) {// Process contigs
			ArrayList<Read> contigs = new ArrayList<Read>();
			for (Read r : scaffolds) {
				ArrayList<Read> temp = r.breakAtGaps(makeAgp, minContig);
				if (bswAgp != null) {
					bswAgp.print((byte[]) r.obj);
					r.obj=null;
				}
				contigs.addAll(temp);
			}
			if (sortContigs) {
				Shared.sort(contigs, ReadLengthComparator.comparator);
			}
			if (renameContigs) {
				for (Read r : contigs) {
					r.id = "contig_" + contigNum;
					contigNum++;
				}
			}
			if (rosc != null) {
				rosc.add(contigs, 0);
			}
		}

		if (bswAgp != null) {
			bswAgp.poisonAndWait();
		}

	}

	/** Iterate through the reads */
	private ArrayList<Read> getReads(final ConcurrentReadInputStream cris) {

		ArrayList<Read> all = new ArrayList<Read>(10000);

		{
			// Grab the first ListNum of reads
			ListNum<Read> ln = cris.nextList();
			// Grab the actual read list from the ListNum
			ArrayList<Read> reads = (ln != null ? ln.list : null);

			// Check to ensure pairing is as expected
			if (reads != null && !reads.isEmpty()) {
				Read r = reads.get(0);
				assert ((ffin1 == null || ffin1.samOrBam()) || (r.mate != null) == cris.paired());
			}

			// As long as there is a nonempty read list...
			while (ln != null && reads != null && reads.size() > 0) {// ln!=null
																		// prevents
																		// a
																		// compiler
																		// potential
																		// null
																		// access
																		// warning
				if (verbose) {
					outstream.println("Fetched " + reads.size() + " reads.");
				}

				// Loop through each read in the list
				for (int idx = 0; idx < reads.size(); idx++) {
					final Read r1 = reads.get(idx);
					assert (r1.mate == null);

					// Track the initial length for statistics
					final int initialLength1 = r1.length();

					// Increment counters
					readsProcessed+=1;
					basesProcessed+=initialLength1;

					boolean keep = processRead(r1);
					if (keep) {
						all.add(r1);
						readsOut+=1;
						basesOut+=initialLength1;
					}
				}

				// Notify the input stream that the list was used
				cris.returnList(ln);
				if (verbose) {
					outstream.println("Returned a list.");
				}

				// Fetch a new list
				ln = cris.nextList();
				reads = (ln != null ? ln.list : null);
			}

			// Notify the input stream that the final list was used
			if (ln != null) {
				cris.returnList(ln.id, ln.list == null || ln.list.isEmpty());
			}
		}

		return all;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Process a single read pair.
	 * 
	 * @param r1
	 *            Read 1
	 * @return True if the reads should be kept, false if they should be
	 *         discarded.
	 */
	boolean processRead(final Read r1) {
		//assert (!banIupac || !r1.containsNonACGTN()) : "Non-ACGTN base found in scaffold " + r1.id;
		if(banIupac){
			if(r1.containsNonACGTN()){
				KillSwitch.exceptionKill(new RuntimeException("Non-ACGTN base found in scaffold " + r1.id));
			}
		}
		r1.inflateGaps(minGapIn, minGapOut);
		return r1.length() >= minScaf;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private int minGapIn = 1;
	private int minGapOut = 10;
	private int minScaf = 1;
	private int minContig = 1;
	private long scafNum = 1;
	private long contigNum = 1;

	private boolean sortScaffolds = true;
	private boolean sortContigs = false;
	private boolean banIupac = true;
	private boolean renameScaffolds = true;
	private boolean renameContigs = false;

	/*--------------------------------------------------------------*/
	/*----------------          I/O Fields          ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	private String in1 = null;

	private String qfin1 = null;

	/** Primary output file path */
	private String out1 = null;
	private String outC = null;

	private String qfout1 = null;
	private String qfoutC = null;

	private String agpFile = null;
	private String legendFile = null;

	/** Override input file extension */
	private String extin = null;
	/** Override output file extension */
	private String extout = null;

	/*--------------------------------------------------------------*/

	/** Number of reads processed */
	protected long readsProcessed = 0;
	/** Number of bases processed */
	protected long basesProcessed = 0;

	/** Number of reads processed */
	protected long readsOut = 0;
	/** Number of bases processed */
	protected long basesOut = 0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads = -1;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file */
	private final FileFormat ffin1;

	/** Primary output file */
	private final FileFormat ffout1;
	private final FileFormat ffoutC;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Print status messages to this output stream */
	private PrintStream outstream = System.err;
	/** Print verbose messages */
	public static boolean verbose = false;
	/** True if an error was encountered */
	public boolean errorState = false;
	/** Overwrite existing output files */
	private boolean overwrite = false;
	/** Append to existing output files */
	private boolean append = false;
	/** This flag has no effect on singlethreaded programs */
	private final boolean ordered = false;

}
