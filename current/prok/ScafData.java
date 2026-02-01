package prok;

import java.util.ArrayList;
import java.util.Arrays;

import gff.GffLine;
import simd.Vector;
import stream.Read;
import structures.IntList;

/**
 * Tracks information about a scaffold for AnalyzeGenes.
 * Manages scaffold sequence data, frame annotations, and gene feature collections
 * for prokaryotic genome analysis. Supports strand-specific annotation storage
 * and sequence manipulation operations including reverse complementing.
 *
 * @author Brian Bushnell
 * @date Sep 24, 2018
 */
class ScafData {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Constructs ScafData from a Read, copying name and bases and allocating frame array.
	 * @param r Read containing scaffold sequence */
	ScafData(Read r){
		this(r.id, r.bases, new byte[r.length()]);
	}
	
	/**
	 * Constructs ScafData with explicit name, bases, and frame annotations; initializes CDS/RNA lists per strand.
	 * @param name_ Scaffold name
	 * @param bases_ Sequence bases
	 * @param frames_ Frame annotation array
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
	
	/** Clears frame annotations and resets start/stop lists. */
	void clear(){
		Arrays.fill(frames, (byte)0);
		starts.clear();
		stops.clear();
	}
	
	/** Reverse-complements bases in place and flips strand. */
	void reverseComplement(){
		Vector.reverseComplementInPlaceFast(bases);
		strand=1^strand;
	}
	
	/** Adds a CDS annotation to the strand-appropriate list.
	 * @param gline CDS GFF line */
	void addCDS(GffLine gline){
		assert(gline.strand>=0) : gline+"\n"+gline.strand;
		cdsLines[gline.strand].add(gline);
	}
	
	/** Adds an RNA annotation to the strand-appropriate list.
	 * @param gline RNA GFF line */
	void addRNA(GffLine gline){
		assert(gline.strand>=0) : gline+"\n"+gline.strand;
		rnaLines[gline.strand].add(gline);
	}
	
	/**
	 * Returns subsequence from start..stop (inclusive) from the scaffold bases.
	 * @param start Start index (inclusive)
	 * #param stop End index (inclusive)
	 * #return Subsequence bytes
	 */
	byte[] fetch(int start, int stop){
		assert(start>=0 && stop<bases.length);
		assert(start<stop);
		return Arrays.copyOfRange(bases, start, stop+1);
	}
	
	/** Returns current strand (0 forward, 1 reverse). */
	int strand(){return strand;}

	/** Returns scaffold length (0 if bases is null).
	 * @return Length in bases */
	public int length() {return bases==null ? 0 : bases.length;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	final String name;
	final byte[] bases;
	final byte[] frames;
	final IntList starts=new IntList(8);
	final IntList stops=new IntList(8);
	private int strand=0;
	
	@SuppressWarnings("unchecked")
	ArrayList<GffLine>[] cdsLines=new ArrayList[2];
	@SuppressWarnings("unchecked")
	ArrayList<GffLine>[] rnaLines=new ArrayList[2];
}
