package prok;

import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import gff.GffLine;
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
	
	ScafData(Read r){
		this(r.id, r.bases, new byte[r.length()]);
	}
	
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
	
	void clear(){
		Arrays.fill(frames, (byte)0);
		starts.clear();
		stops.clear();
	}
	
	void reverseComplement(){
		AminoAcid.reverseComplementBasesInPlace(bases);
		strand=1^strand;
	}
	
	void addCDS(GffLine gline){
		assert(gline.strand>=0) : gline+"\n"+gline.strand;
		cdsLines[gline.strand].add(gline);
	}
	
	void addRNA(GffLine gline){
		assert(gline.strand>=0) : gline+"\n"+gline.strand;
		rnaLines[gline.strand].add(gline);
	}
	
	byte[] fetch(int start, int stop){
		assert(start>=0 && stop<bases.length);
		assert(start<stop);
		return Arrays.copyOfRange(bases, start, stop+1);
	}
	
	int strand(){return strand;}

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
	
	/** gLines[strand] holds the GffLines for that strand */
	@SuppressWarnings("unchecked")
	ArrayList<GffLine>[] cdsLines=new ArrayList[2];
	@SuppressWarnings("unchecked")
	ArrayList<GffLine>[] rnaLines=new ArrayList[2];
}
