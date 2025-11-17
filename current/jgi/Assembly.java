package jgi;

import java.util.Arrays;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Tools;
import structures.IntList;

/**
 * Parses and analyzes FASTA genome assembly files to extract sequence statistics.
 * Loads FASTA files, computes nucleotide counts, tracks contig lengths, and
 * calculates genome composition metrics including GC content and length distributions.
 * @author Brian Bushnell
 */
public class Assembly {
	
	/** Creates a new Assembly by loading and parsing the specified FASTA file.
	 * @param fname_ Path to the FASTA file to analyze */
	public Assembly(String fname_) {
		fname=fname_;
		load();
	}
	
	/**
	 * Loads and parses the FASTA file, extracting contig statistics and base composition.
	 * Processes each sequence header and body, accumulating contig lengths and
	 * nucleotide counts. Sorts contigs by length in descending order after loading.
	 */
	void load() {
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		assert(ff.fasta());
		ByteFile bf=ByteFile.makeByteFile(ff);

		clear();
		acgtnio=new long[7];
		int contigLen=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()) {
			if(line[0]=='>') {
				headerLength+=(line.length-1);
				if(firstHeader==null && line.length>1) {
					firstHeader=new String(line, 1, line.length-1);
				}
				if(contigLen>0) {
					contigs.add(contigLen);
				}
				contigLen=0;
			}else {
				addToACGTNIO(line);
				contigLen+=line.length;
			}
		}
		bf.close();
		contigs.sort();
		contigs.reverse();
		length=contigs.sumLong();
	}
	
	/**
	 * Counts nucleotide bases from a sequence line and updates composition statistics.
	 * Processes each byte in the line using the baseToACGTNIO mapping table.
	 * @param line Sequence data as byte array
	 */
	void addToACGTNIO(byte[] line) {
		for(byte b : line) {
			byte x=baseToACGTNIO[b];
			acgtnio[x]++;
		}
	}
	
	/** Resets all assembly statistics to initial state.
	 * Clears contig list, resets lengths and counters to zero. */
	void clear() {
		contigs.clear();
		length=0;
		headerLength=0;
		firstHeader=null;
		acgtnio=null;
	}
	
	/** Calculates GC content as fraction of total non-N bases.
	 * @return GC content ratio (G+C)/(A+T+U+G+C) */
	float gc() {
		float AT=acgtnio[A]+acgtnio[T]+acgtnio[U];
		float GC=acgtnio[G]+acgtnio[C];
		return GC/Tools.max(1, AT+GC);
	}
	
	/**
	 * Calculates total length of contigs that are at least the specified minimum length.
	 * Utilizes sorted contig list for efficient early termination.
	 * @param minimum Minimum contig length threshold
	 * @return Sum of lengths for contigs >= minimum length
	 */
	long lengthAtLeast(int minimum) {
		long sum=0;
		for(int i=0; i<contigs.size; i++) {
			int len=contigs.get(i);
			if(len<minimum) {break;}
			sum+=len;
		}
		return sum;
	}
	
	/** Path to the FASTA file being analyzed */
	final String fname;
	/** List of contig lengths, sorted in descending order after loading */
	IntList contigs=new IntList();
	/** Total assembly length (sum of all contig lengths) */
	long length=0;
	/** Total length of all sequence headers in the file */
	long headerLength=0;
	/** First sequence header found in the file (without '>' character) */
	String firstHeader=null;
	/**
	 * Array counting occurrences of each base type: A, C, G, T, U, N, IUPAC, OTHER
	 */
	long[] acgtnio;
	
	/**
	 * Lookup table mapping ASCII values to base type indices (A=0, C=1, G=2, T=3, U=4, N=5, IUPAC=6, OTHER=7)
	 */
	public static final byte[] baseToACGTNIO=makeBaseToACGTUNIO();
	/** Index constant for non-standard bases in acgtnio array */
	/** Index constant for IUPAC ambiguous bases in acgtnio array */
	/** Index constant for ambiguous N bases in acgtnio array */
	/** Index constant for uracil bases in acgtnio array */
	/** Index constant for thymine bases in acgtnio array */
	/** Index constant for guanine bases in acgtnio array */
	/** Index constant for cytosine bases in acgtnio array */
	/** Index constant for adenine bases in acgtnio array */
	private static final byte A=0, C=1, G=2, T=3, U=4, N=5, IUPAC=6, OTHER=7;
	
	/**
	 * Creates the base-to-index mapping table for nucleotide classification.
	 * Maps ASCII values to appropriate indices: standard bases (ACGTU), ambiguous (N),
	 * IUPAC codes, and others. Handles both uppercase and lowercase letters.
	 * @return Lookup table mapping ASCII values to base type indices
	 */
	private static final byte[] makeBaseToACGTUNIO() {
		final byte[] array=new byte[128];
		Arrays.fill(array, OTHER);
		array['a']=array['A']=A;
		array['c']=array['C']=C;
		array['g']=array['G']=G;
		array['t']=array['T']=T;
		array['u']=array['U']=U;
		array['n']=array['N']=N;
		for(int i=0; i<array.length; i++) {
			if(AminoAcid.baseToNumberExtended[i]>=0 && array[i]==OTHER) {
				array[i]=IUPAC;
			}
		}
		return array;
	}
	
}
