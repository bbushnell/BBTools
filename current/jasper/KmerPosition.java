package jasper;

import java.util.ArrayList;
import java.util.HashSet;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;
import structures.LongList;

/**
 * Analyzes positional k-mer distribution in high-throughput sequencing reads.
 * Processes DNA sequencing data to identify k-mer occurrences at specific
 * positions within reads, comparing against a reference sequence k-mer set.
 * Supports paired-end read analysis with detailed positional statistics.
 *
 * @author Brian Bushnell
 * @author Jasper Toscani Field
 * @date Jun 4, 2020
 */
public class KmerPosition{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs {@link KmerPosition} from the command line.
	 * @param args command-line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();

		//Create an instance of this class
		KmerPosition x=new KmerPosition(args);

		//Run the object
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Parses command-line arguments and validates input, reference and output files.
	 * @param args BBTools-style {@code flag=value} arguments; values may contain '='
	 */
	public KmerPosition(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
				// Recognized legacy no-op placeholder.
			}else if(a.equals("ref")){
				ref=b;

			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				throw new IllegalArgumentException("Unknown parameter "+arg);
			}
		}

		{//Process parser fields
			Parser.processQuality();

			maxReads=parser.maxReads;
			in1=parser.in1;
			in2=parser.in2;
			out1=parser.out1;
		}

		if(in1==null){
			throw new IllegalArgumentException("Please specify an input file with in=<reads>.");
		}
		if(ref==null){
			throw new IllegalArgumentException("Please specify a reference file with ref=<fasta>.");
		}
		if(k<1){throw new IllegalArgumentException("k must be at least 1; got k="+k);}

		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
		ffref=FileFormat.testInput(ref, FileFormat.FASTA, null, true, true);
	}

	/**
	 * Loads reference k-mers, scans all reads, writes one row per read start position,
	 * and prints runtime statistics.
	 * @param t active timer to stop after output has finished
	 */
	void process(Timer t){
		HashSet<String> kr=kmerReturn();

		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2);
			cris.start();
		}
		boolean paired=cris.paired();

		long readsProcessed=0, basesProcessed=0;
		{

			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);

				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			// ln!=null prevents a compiler potential null access warning.
			while(ln!=null && reads!=null && reads.size()>0){
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}


				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;
					readsProcessed+=r1.pairCount();
					basesProcessed+=r1.pairLength();

					//  *********  Process reads here  *********
					processRead(r1, kr, counts1, totalEncounter1);
					if(r1.mate!=null){
						processRead(r2, kr, counts2, totalEncounter2);
					}
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
		errorState=ReadWrite.closeStreams(cris) | errorState;
		if(verbose){outstream.println("Finished reading data.");}

		outputResults(counts1, totalEncounter1, counts2, totalEncounter2);

		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+
				Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
		if(errorState){throw new RuntimeException("An error was encountered while closing streams.");}
	}

	/**
	 * Writes positional match counts and percentages for first and second mates.
	 * @param posCounts1 reference-matching k-mer counts by first-mate start position
	 * @param readCounts1 all String k-mer windows by first-mate start position
	 * @param posCounts2 reference-matching k-mer counts by second-mate start position
	 * @param readCounts2 all String k-mer windows by second-mate start position
	 */
	private void outputResults(LongList posCounts1, LongList readCounts1, LongList posCounts2,
			LongList readCounts2){
		if(ffout1==null){return;}
		ByteStreamWriter bsw=new ByteStreamWriter(ffout1);
		bsw.start();

		long[] readArray1=readCounts1.toArray();
		long[] countArray1=posCounts1.toArray();
		long[] readArray2=readCounts2.toArray();
		long[] countArray2=posCounts2.toArray();

		bsw.println("#pos\tread1_count\tread1_perc\tread2_count\tread2_perc");

		int maxLen=Tools.max(readArray1.length, readArray2.length);

		for(int i=0; i<maxLen; i++){

			bsw.print(i);
			bsw.print('\t');
			bsw.print(countArray1.length>i ? countArray1[i] : 0);
			bsw.print('\t');
			bsw.print(percent(countArray1, readArray1, i), 3);
			bsw.print('\t');
			bsw.print(countArray2.length>i ? countArray2[i] : 0);
			bsw.print('\t');
			bsw.print(percent(countArray2, readArray2, i), 3);
			bsw.println();
		}
		errorState=bsw.poisonAndWait() | errorState;
	}

	/** Returns a matched-window percentage, or zero when no windows reached {@code pos}.
	 * @param countArray match counts by position
	 * @param readArray total window counts by position
	 * @param pos position to report
	 * @return percent of windows that matched the reference */
	private static float percent(long[] countArray, long[] readArray, int pos){
		if(countArray.length<=pos || readArray.length<=pos || readArray[pos]<1){return 0;}
		return (countArray[pos]/(float)readArray[pos])*100;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads the reference FASTA and returns every forward String k-mer of length k.
	 * @return set of reference k-mers in their literal read orientation
	 */
	private HashSet<String> kmerReturn(){
		HashSet<String> hs=new HashSet<String>();
		ArrayList<Read> readArray=
				ConcurrentReadInputStream.getReads(maxReads, false, ffref, null, null, null);
		for(Read r : readArray){
			addToSet(hs, r);
		}

		return hs;
	}

	/**
	 * Adds all String k-mers from a reference sequence to a set.
	 * @param hs destination set for literal forward k-mers
	 * @param r reference sequence to scan
	 * @return number of k-mer windows visited, including duplicates already in the set
	 */
	private int addToSet(HashSet<String> hs, Read r){
		int countRead=0;
		for(int i=0, j=k; j<=r.length(); i++, j++){
			//String(byte[] bytes, int offset, int length)
			String s=new String(r.bases, i, k);
			hs.add(s);
			countRead++;
		}
		return countRead;
	}

	/**
	 * Increments total and matching k-mer counters by start position for one read.
	 * @param r read sequence to scan in forward orientation
	 * @param hs set of forward reference k-mers
	 * @param count matching k-mer count list to increment by start position
	 * @param readCount total k-mer opportunity list to increment by start position
	 * @return {@code count}
	 */
	private LongList processRead(Read r, HashSet<String> hs, LongList count, LongList readCount){
		for(int i=0, j=k; j<=r.length(); i++, j++){
			//String(byte[] bytes, int offset, int length)

			String s=new String(r.bases, i, k);
			readCount.increment(i);
			if(hs.contains(s)){
				count.increment(i);
			}
		}
		return count;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary read input path. */
	private String in1=null;

	/** Optional second-mate read input path. */
	private String in2=null;

	/** Positional count output path. */
	private String out1=null;

	/** Reference FASTA path. */
	private String ref=null;

	/** Primary read input format. */
	private final FileFormat ffin1;

	/** Optional second-mate input format. */
	private final FileFormat ffin2;

	/** Tabular count output format. */
	private final FileFormat ffout1;

	/** Reference FASTA input format. */
	private final FileFormat ffref;

	/*--------------------------------------------------------------*/

	/** Maximum reference fragments and input reads to process, or -1 for all. */
	private long maxReads=-1;

	/** True if any stream reported an error. */
	private boolean errorState=false;

	/** K-mer length in bases. */
	private int k=6;

	/** Positional counts of matching first-mate k-mers. */
	private LongList counts1=new LongList();

	/** Positional counts of all first-mate String k-mer windows, including ambiguous bases. */
	private LongList totalEncounter1=new LongList();

	/** Positional counts of matching second-mate k-mers. */
	private LongList counts2=new LongList();

	/** Positional counts of all second-mate String k-mer windows, including ambiguous bases. */
	private LongList totalEncounter2=new LongList();

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Destination for progress and timing messages. */
	private java.io.PrintStream outstream=System.err;

	/** Print extra stream activity when true. */
	public static boolean verbose=false;

}
