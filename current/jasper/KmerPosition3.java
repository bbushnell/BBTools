package jasper;

import java.util.ArrayList;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.LongHashSet;
import parse.Parse;
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
 * Identifies k-mer sequence matches between reads and reference sequences.
 * Processes high-throughput read files to find and report k-mer positions that match
 * a reference sequence, tracking occurrence frequencies at each position.
 * Uses 2-bit nucleotide encoding for memory-efficient k-mer comparison.
 * Supports paired-end read analysis with position-specific statistics.
 * Example: with read ACGTA, reference ATGTACC, and k=3 the shared k-mer GTA
 * begins at read position 2 and the output reports counts and percentages for
 * that position.
 * @author Jasper Toscani Field
 * @author Brian Bushnell
 * @date Jun 4, 2020
 */
public class KmerPosition3{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs {@link KmerPosition3} from the command line.
	 * @param args command-line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();

		//Create an instance of this class
		KmerPosition3 x=new KmerPosition3(args);

		//Run the object
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Parses command-line arguments and validates input, reference and output files.
	 * @param args BBTools-style {@code flag=value} arguments; values may contain '='
	 */
	public KmerPosition3(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		//Primary parsing of standard BBTools arguments.
		Parser parser=new Parser();

		//Loop through arguments up to the maximum number of arguments input.
		//process all remaining arguments.
		for(int i=0; i<args.length; i++){

			//Grab argument string at index.
			String arg=args[i];

			//Split argument string on "=".
			String[] split=arg.split("=", 2);

			//Convert the left side to lowercase.
			String a=split[0].toLowerCase();

			//Use the right side of the split if present.
			String b=split.length>1 ? split[1] : null;

			//If b isn't null but a string "null" was input, convert b to null.
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			// Recognized legacy no-op placeholder.
			if(a.equals("parse_flag_goes_here")){

			//Handle reference variable assignment.
			}else if(a.equals("ref")){
				ref=b;

			//Handle kmer variable assignment.
			}else if(a.equals("k")){
				k=Integer.parseInt(b);

			//Parse the reverse-complemented reference k-mer toggle.
			}else if(a.equals("rcomp")){
				rcomp=Parse.parseBoolean(b);

			//Parse common input/output flags and shared parser options.
			}else if(parser.parse(arg, a, b)){

			//If not one of the known parameters, let the user know they made a mistake.
			}else{
				throw new IllegalArgumentException("Unknown parameter "+arg);
			}
		}

		{//Handle quality scoring by identifying quality scoring method.
			Parser.processQuality();

			//appropriate argument passing
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
		//[jasper/KmerPosition3#001] 2-bit packing fits only k in [1,32] (64-bit long).
		if(k<1 || k>32){
			throw new IllegalArgumentException(
					"k must be 1-32 for 2-bit (64-bit long) encoding; got k="+k);
		}

		//File format handling for each file.
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, true, false, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
		ffref=FileFormat.testInput(ref, FileFormat.FASTA, null, true, true);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Primary processing function for k-mer position analysis.
	 * Loads reference k-mers, creates read input stream, processes all reads
	 * to find matching k-mers at each position, validates paired vs single-end state,
	 * returns processed read lists, writes tab-delimited statistics, and reports timing.
	 * @param t Timer object for tracking program execution time
	 */
	void process(Timer t){

		//Creates a empty LongHashSet.
		//This hash set will eventually hold kmers found in the reference of length k,
		//advancing 1 position each step.
		LongHashSet refKmerSet=loadReference();

		//Instantiate the input stream
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2);
			cris.start();
		}
		//Paired indicates whether input stream is processing the data as paired.
		boolean paired=cris.paired();


		long readsProcessed=0, basesProcessed=0;
		{
			//Returns listNum of 200 reads,
			//minimize function calling to reduce inter-thread communication.
			ListNum<Read> ln=cris.nextList();

			//Prevents null pointer exception if ListNum is null
			//meaning file was empty.
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			//Pulls first read in list, checks if paired input was selected and
			//determines if this is correct.
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);

				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			//Loop while more reads available.
			// ln!=null prevents a compiler potential null access warning.
			while(ln!=null && reads!=null && reads.size()>0){
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}

				//Loop through every read in list.
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					final Read r2=r1.mate;

					//Handles incrementing the number of reads processed if paired or unpaired.
					readsProcessed+=r1.pairCount();
					basesProcessed+=r1.pairLength();

					//  *********  Process reads here  *********
					//Pass read 1 (and 2 if paired) and the reference kmer set, with the positional count list.
					processRead(r1, refKmerSet, matchCounts1, totalCounts1);
					if(r1.mate!=null){
						processRead(r2, refKmerSet, matchCounts2, totalCounts2);
					}
				}

				//When done processing list, return to input stream to
				//notify we're ready for new list of reads.
				cris.returnList(ln);

				if(verbose){outstream.println("Returned a list.");}

				//Grab next list.
				ln=cris.nextList();

				//Prevent null pointer exception if list of reads is empty.
				reads=(ln!=null ? ln.list : null);
			}

			//If no list returned on line 176, return list now.
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}

		//Close stream after processing file, collect error if encountered.
		errorState=ReadWrite.closeStreams(cris) | errorState;
		if(verbose){outstream.println("Finished reading data.");}

		//Writes output to file indicated by user.
		outputResults(matchCounts1, totalCounts1, matchCounts2, totalCounts2);

		//Stop timer after all processes have competed and before printing runtime.
		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+
				Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
		outstream.println("Bases Processed:    "+basesProcessed+" \t"+
				Tools.format("%.2fk bases/sec", (basesProcessed/(double)(t.elapsed))*1000000));

		if(errorState){throw new RuntimeException("An error was encountered while closing streams.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Converts k-mer count lists to arrays and writes analysis statistics to output file.
	 * Outputs position-specific counts and percentages for both read sets in tabular format.
	 * Handles cases where read sets have different lengths by using the maximum length.
	 * @param posCounts1 List tracking k-mer counts at each position for read set 1
	 * @param readCounts1 List tracking total reads with nucleotides at each position for read set 1
	 * @param posCounts2 List tracking k-mer counts at each position for read set 2
	 * @param readCounts2 List tracking total reads with nucleotides at each position for read set 2
	 */
	private void outputResults(LongList posCounts1, LongList readCounts1, LongList posCounts2,
			LongList readCounts2){
		//Makes sure a valid output file name exists.
		if(ffout1==null){return;}

		//Begins output writer method.
		ByteStreamWriter bsw=new ByteStreamWriter(ffout1);
		bsw.start();

		//Converts count lists to arrays for indexing.
		long[] readArray1=readCounts1.toArray();
		long[] countArray1=posCounts1.toArray();
		long[] readArray2=readCounts2.toArray();
		long[] countArray2=posCounts2.toArray();

		//Writes the header line to the output file.
		bsw.println("#pos\tread1_count\tread1_perc\tread2_count\tread2_perc");

		// A valid k-mer can appear after leading ambiguous bases; emit zero
		// explicitly for positions with no valid denominator.
		//Finds the maximum of both read set lengths.
		//This handles if one read set is longer than the other.
		//Its important to use the longer length to avoid iterating out of bounds.
		int maxLen=Tools.max(readArray1.length, readArray2.length);

		//Iterate to the length of the longest read (if the value exists, otherwise report 0) and
		//write counts to output file.
		for(int i=0; i<maxLen; i++){

			//Write to file the position in reads
			bsw.print(i);
			bsw.print('\t');

			//Write to file the number of kmers found at position i.
			bsw.print(countArray1.length>i ? countArray1[i] : 0);
			bsw.print('\t');

			//Write to file the percentage of reads with a kmer at this particular position.
			bsw.print(percent(countArray1, readArray1, i), 3);
			bsw.print('\t');

			//Same statistics as above for the read 2 counts.
			bsw.print(countArray2.length>i ? countArray2[i] : 0);
			bsw.print('\t');

			//Same percentage statistic as above for read 2 kmers at this position.
			bsw.print(percent(countArray2, readArray2, i), 3);
			bsw.println();
		}

		//Stop the output stream by telling separate thread no more data is incoming
		//but to finish writing already passed data.
		errorState=bsw.poisonAndWait() | errorState;
	}

	/** Returns a matched-window percentage, or zero when no valid k-mers reached pos.
	 * @param countArray match counts by position
	 * @param readArray total valid k-mer counts by position
	 * @param pos position to report
	 * @return percent of valid k-mers that matched the reference */
	private static float percent(long[] countArray, long[] readArray, int pos){
		if(countArray.length<=pos || readArray.length<=pos || readArray[pos]<1){return 0;}
		return (countArray[pos]/(float)readArray[pos])*100;
	}

	/**
	 * Creates a set of binary reference k-mers of length k.
	 * Reads reference sequences and converts k-mers to 2-bit binary representation
	 * for memory-efficient storage and fast comparison with read k-mers.
	 * When {@code rcomp} is true, reverse complements are stored with the
	 * forward k-mers.
	 * @return set of reference k-mers in binary format
	 */
	private LongHashSet loadReference(){
		//Initialize empty LongHashSet to accept reference kmers.
		LongHashSet hs=new LongHashSet();

		//Get the reference sequences from the ref file.
		ArrayList<Read> readArray=
				ConcurrentReadInputStream.getReads(maxReads, false, ffref, null, null, null);

		//iterate over sequences pulled from reference file,
		//and pass the sequence to the byte conversion method
		for(Read r : readArray){
			addToSet(hs, r);
		}

		return hs;
	}

	/**
	 * Converts reference sequence k-mers to 2-bit binary representation and adds to hashset.
	 * Uses sliding window approach with bit-shifting to generate consecutive k-mers.
	 * Skips k-mers containing degenerate bases and optionally includes reverse complements.
	 * @param hs HashSet to store the converted k-mer values
	 * @param r Reference sequence read to process
	 * @return number of k-mer windows visited, including duplicates already in the set
	 */
	private int addToSet(LongHashSet hs, Read r){
		int proccessedKmers=0;
		//This is the old string-based code. Useful for comparisons to the new code below.
		/*for(int i=0, j=k; j<=r.length(); i++, j++) {
			//String(byte[] bytes, int offset, int length)
			String s=new String(r.bases, i, k);
			hs.add(s);
			countRead++;
		}*/

		//Convert kmer sequences to 2-bit notation, shifting each kmer by one nucleotide
		final int shift=2*k;

		//Make mask of 1's for the last k*2 bits
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		if(verbose){System.out.println(Long.toBinaryString(mask));}

		//number of consecutive, valid (non-degenerate) bases
		int len=0;

		//access read objects bases
		byte[] bases=r.bases;

		//binary representation of current kmer
		long kmer=0;

		//iterate over bases in read
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];

			//Fast nucleotide conversion from letter to 2-bit encoding
			long x=AminoAcid.baseToNumber[b];

			//Shift kmer to left, dropping oldest base, OR in new base on right
			kmer=((kmer<<2)|x)&mask;


			//ACGTTNTGCGC

			//This section handles creating kmers without degenerate bases at any position
			//Check new base is not degenerate (-1) and increment length of kmer
			if(x>=0){
				len++;
			}else{
				//If base is degenerate (-1), restart kmer construction by setting length to 0
				len=0;
				kmer=0;
			}

			//Once kmer reaches length k, add kmer to hashset and increment number of processed kmers
			if(len>=k){
				hs.add(kmer);
				if(rcomp){hs.add(AminoAcid.reverseComplementBinaryFast(kmer, k));}
				proccessedKmers++;
			}
		}

		return proccessedKmers;
	}

	/**
	 * Processes individual read using 2-bit encoding for fast k-mer comparison.
	 * Generates k-mers from read sequence, compares against reference k-mer set,
	 * and increments position-specific counters for matches and total coverage.
	 * Uses bit-shifting for efficient k-mer generation without string operations.
	 * @param r Read sequence to analyze
	 * @param hs Set of reference k-mers in binary format
	 * @param matchCounts List to increment when k-mer matches are found at positions
	 * @param totalCounts List to increment for fully defined k-mers at each position
	 */
	private void processRead(Read r, LongHashSet hs, LongList matchCounts, LongList totalCounts){

		//Convert kmer sequences to 2-bit notation, shifting each kmer by one nucleotide
		final int shift=2*k;

		//Make mask of 1's for the last k*2 bits
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		if(verbose){System.out.println(Long.toBinaryString(mask));}

		//number of consecutive, valid (non-degenerate) bases
		int len=0;

		//access read objects bases
		byte[] bases=r.bases;

		//binary representation of current kmer
		long kmer=0;

		//iterate over bases in read up to the length of the read
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];

			//Fast nucleotide conversion from letter to 2-bit encoding
			long x=AminoAcid.baseToNumber[b];

			//Shift kmer to left, dropping oldest base, OR in new base on right
			kmer=((kmer<<2)|x)&mask;

			//This section handles creating kmers without degenerate bases at any position
			//Check new base is not degenerate (-1) and increment length of kmer
			if(x>=0){
				len++;
			}else{

				//If base is degenerate (-1), restart kmer construction by setting length to 0
				len=0;
				kmer=0;
			}
			if(len>=k){

				//increment list of counts for reads containing nucleotides at each position
				//i - k + 1 is the first base of kmer, i is last base in kmer
				//we want start positions
				totalCounts.increment(i-k+1);

				//if the read kmer is in the hashset of kmers from the reference,
				//increment count of positions corresponding to the start of a kmer at that position
				if(hs.contains(kmer)){
					matchCounts.increment(i-k+1);
				}
			}
		}

	}



	//A -> 0 -> 00
	//C -> 1 -> 01
	//G -> 2 -> 10
	//T -> 3 -> 11
	//N -> -1 -> 11111111111111111111111111111111111111111

	//00110110 -> ATCG


	// kmer=00000000
	//Add T, 11                           G A A A        G A A A A
	//left shift:  kmer=kmer<<2  -> kmer=10000000<<2 -> 1000000000
	//kmer=00000000                                     0011111111
	//Or it with the new code                           0000000000
	//kmer=kmer|x  ->  kmer=00000000 | 11 -> 0000000011
	//mask it with the mask:
	//kmer=kmer&mask -> kmer=00000000 & 11111111 -> 00000011

	//kmer=00000011
	//Add C, 01
	//left shift : -> 00001100
	//Or it with the new code:                A A T C
	//kmer=kmer|x  ->  kmer=00001100 | 01 -> 00001101
	//mask (does nothing)

	// kmer= TATC = 11001101
	//Add G, 10
	// left-shift:  1100110100 (TATCA)
	//Or with the new code: kmer=1100110100 | 10 = 1100110110 = TATCG
	//Mask with 11111111 (TTTT): 1100110110
	//                             11111111
	//yields                       00110110 = ATCG

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file for high-throughput read sequences */
	private String in1=null;

	/** Paired-end read input file containing mates of in1 sequences */
	private String in2=null;

	/** Output file for k-mer positioning statistics and counts */
	private String out1=null;

	/** Reference sequence file in FASTA format for k-mer matching */
	private String ref=null;

	/** File format handler for primary input file parsing and validation */
	private final FileFormat ffin1;

	/** File format handler for paired-end input file parsing and validation */
	private final FileFormat ffin2;

	/** File format handler for output file structure and writing methods */
	private final FileFormat ffout1;

	/** File format handler for reference file assuming FASTA format */
	private final FileFormat ffref;

	/*--------------------------------------------------------------*/

	/** Maximum reference fragments and input reads to analyze; -1 indicates all. */
	private long maxReads=-1;

	/** Error flag that causes program to exit with status 1 when true */
	private boolean errorState=false;

	/** Flag to include reverse-complemented k-mers in the reference hashset */
	private boolean rcomp=true;

	/**
	 * K-mer length for sequence analysis; configurable via k=# command-line flag
	 */
	private int k=19;

	/** Position-specific counts of matching k-mers found in read set 1 */
	private LongList matchCounts1=new LongList();

	/** Position-specific counts of fully defined k-mers for read set 1 */
	private LongList totalCounts1=new LongList();

	/** Position-specific counts of matching k-mers found in read set 2 */
	private LongList matchCounts2=new LongList();

	/** Position-specific counts of fully defined k-mers for read set 2 */
	private LongList totalCounts2=new LongList();

	/*--------------------------------------------------------------*/

	/** Output stream for piping program statistics to output destination */
	private java.io.PrintStream outstream=System.err;

	/** Verbose mode flag for printing additional program execution information */
	public static boolean verbose=false;

}
