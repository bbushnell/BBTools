package prok;

import java.util.ArrayList;
import java.util.Arrays;

import gff.GffLine;
import shared.Vector;
import stream.Read;
import structures.IntList;

/**
 * Tracks information about a scaffold for AnalyzeGenes.
 * @author Brian Bushnell
 * @date Sep 24, 2018
 *
 */
class ScafData {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates ScafData from a Read object.
	 * Initializes with read ID, bases, and empty frame array.
	 * @param r Read object containing sequence data
	 */
	ScafData(Read r){
		this(r.id, r.bases, new byte[r.length()]);
	}
	
	/**
	 * Creates ScafData with specified name, bases, and frame annotations.
	 * Initializes empty ArrayLists for CDS and RNA annotations on both strands.
	 *
	 * @param name_ Scaffold identifier
	 * @param bases_ Nucleotide sequence as byte array
	 * @param frames_ Frame annotation array parallel to bases
	 */
	ScafData(String name_, byte[] bases_, byte[] frames_){
		name=name_;
		bases=bases_;
		frames=frames_;
		cdsLines[0]=new ArrayList<GffLine>();
		cdsLines[1]=new ArrayList<GffLine>();
		rnaLines[0]=new ArrayList<GffLine>();
		rnaLines[1]=new ArrayList<GffLine>();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Clears frame annotations and start/stop position lists.
	 * Resets frames array to zeros and empties coordinate tracking lists. */
	void clear(){
		Arrays.fill(frames, (byte)0);
		starts.clear();
		stops.clear();
	}
	
	/** Reverse complements the scaffold sequence and flips strand orientation.
	 * Modifies bases array in-place and toggles strand between 0 and 1. */
	void reverseComplement(){
		Vector.reverseComplementInPlaceFast(bases);
		strand=1^strand;
	}
	
	/** Adds a coding sequence annotation to the appropriate strand collection.
	 * @param gline GFF line representing a CDS feature with valid strand */
	void addCDS(GffLine gline){
		assert(gline.strand>=0) : gline+"\n"+gline.strand;
		cdsLines[gline.strand].add(gline);
	}
	
	/** Adds an RNA annotation to the appropriate strand collection.
	 * @param gline GFF line representing an RNA feature with valid strand */
	void addRNA(GffLine gline){
		assert(gline.strand>=0) : gline+"\n"+gline.strand;
		rnaLines[gline.strand].add(gline);
	}
	
	/**
	 * Extracts a subsequence from the scaffold bases.
	 * Returns inclusive range from start to stop positions.
	 *
	 * @param start Starting position (inclusive)
	 * @param stop Ending position (inclusive)
	 * @return Subsequence as byte array
	 */
	byte[] fetch(int start, int stop){
		assert(start>=0 && stop<bases.length);
		assert(start<stop);
		return Arrays.copyOfRange(bases, start, stop+1);
	}
	
	/** Returns current strand orientation (0 or 1) */
	int strand(){return strand;}

	/** Returns scaffold length, or 0 if bases is null */
	public int length() {return bases==null ? 0 : bases.length;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Scaffold identifier */
	final String name;
	/** Nucleotide sequence as byte array */
	final byte[] bases;
	/** Frame annotation array parallel to bases */
	final byte[] frames;
	/** List of start codon positions for gene detection */
	final IntList starts=new IntList(8);
	/** List of stop codon positions for gene detection */
	final IntList stops=new IntList(8);
	/** Current strand orientation: 0 for forward, 1 for reverse */
	private int strand=0;
	
	/** gLines[strand] holds the GffLines for that strand */
	@SuppressWarnings("unchecked")
	ArrayList<GffLine>[] cdsLines=new ArrayList[2];
	/** RNA annotations indexed by strand [0=forward, 1=reverse] */
	@SuppressWarnings("unchecked")
	ArrayList<GffLine>[] rnaLines=new ArrayList[2];
}
