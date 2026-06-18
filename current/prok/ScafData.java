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
		//both strands' lists allocated up front so addCDS/addRNA can index [strand] without null checks
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
		//TODO: Possible bug [prok/ScafData#003] - QUESTION: resets frames/starts/stops but NOT cdsLines/rnaLines (annotations persist across clear()). Intentional (per-strand re-scan keeps annotations) or a missing reset? Trace callers to confirm.
		Arrays.fill(frames, (byte)0);
		starts.clear();
		stops.clear();
	}
	
	/** Reverse-complements bases in place and flips strand. */
	void reverseComplement(){
		Vector.reverseComplementInPlaceFast(bases);
		strand=1^strand;//strand provably in {0,1}: init 0, ctors never set it, this is the only mutation
	}
	
	/** Adds a CDS annotation to the strand-appropriate list.
	 * @param gline CDS GFF line */
	void addCDS(GffLine gline){
		//TODO: Possible bug [prok/ScafData#001] - assert guards strand>=0 (catches -1) but NOT the upper bound. GffLine.strand domain is {-1,0,1,2,3} (PLUS,MINUS,QMARK '?',DOT '.'; GffLine.java:697); strand 2 or 3 passes the assert then indexes cdsLines[2/3] -> AIOOBE (array size 2). Reachable (MEDIUM) iff a CDS with '?'/'.' strand reaches here via AnalyzeGenes; latent LOW if only internal 0/1 ORFs do. Trace AnalyzeGenes. Same issue in addRNA.
		assert(gline.strand>=0) : gline+"\n"+gline.strand;
		cdsLines[gline.strand].add(gline);//requires strand in {0,1}; see ScafData#001
	}
	
	/** Adds an RNA annotation to the strand-appropriate list.
	 * @param gline RNA GFF line */
	void addRNA(GffLine gline){
		assert(gline.strand>=0) : gline+"\n"+gline.strand;
		rnaLines[gline.strand].add(gline);//same as ScafData#001: strand 2/3 ('?'/'.') would AIOOBE here too
	}
	
	/**
	 * Returns subsequence from start..stop (inclusive) from the scaffold bases.
	 * @param start Start index (inclusive)
	 * @param stop End index (inclusive)
	 * @return Subsequence bytes
	 */
	byte[] fetch(int start, int stop){
		assert(start>=0 && stop<bases.length);
		//TODO: Possible bug [prok/ScafData#002] - assert(start<stop) is STRICTER than the inclusive contract: start==stop is a valid 1-base fetch (copyOfRange handles it) but asserts-out. Latent; reachable iff a caller requests a 1-base subsequence. Trace GeneCaller/callers.
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
