package prok;

import java.util.ArrayList;
import java.util.HashMap;

import dna.AminoAcid;
import gff.GffLine;
import ml.CellNet;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * A scored open reading frame (gene candidate) on a scaffold: coordinates, reading frame, packed
 * start/stop codons, the component scores (start/stop/kmer) and the combined orfScore, plus the
 * per-strand path-finding fields (pathScore/prev/pathLength) used to chain ORFs into a gene model.
 * Extends PFeature for coordinate/strand handling.
 *
 * @author Brian
 */
public class Orf extends PFeature {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Constructs an ORF (Open Reading Frame) with genomic coordinates and sequence data.
	 * Extracts start and stop codons from the provided bases and optionally flips
	 * coordinates for minus-strand features.
	 *
	 * @param scafName_ Scaffold/contig name containing this ORF
	 * @param start_ 0-based start position (first base of start codon)
	 * @param stop_ 0-based stop position (last base of stop codon)
	 * @param strand_ Strand orientation (plus or minus)
	 * @param frame_ Reading frame (0, 1, or 2)
	 * @param bases Nucleotide sequence containing the ORF
	 * @param flip Whether to flip coordinates if on minus strand
	 * @param type_ Feature type (CDS, tRNA, rRNA, etc.)
	 */
	public Orf(String scafName_, int start_, int stop_, int strand_, int frame_, byte[] bases, boolean flip, int type_) {
		this(scafName_, start_, stop_, strand_, frame_, bases, flip, type_, true, true);
	}

	/** Carries known biological CDS boundaries independently of coordinate flipping.
	 * The legacy constructor leaves both boundaries unknown/truncated; RNA ignores these flags. */
	public Orf(String scafName_, int start_, int stop_, int strand_, int frame_, byte[] bases,
			boolean flip, int type_, boolean startTruncated_, boolean stopTruncated_){
		super(scafName_, start_, stop_, strand_, bases.length);
		frame=frame_;
		final boolean rna=type_!=CDS;
		if(rna && (start<0 || stop<start || stop>=bases.length)){
			throw new IllegalArgumentException("RNA coordinates must remain inside the contig for path/GFF output: "+start+".."+stop+" length="+bases.length);
		}
		// A valid1/2-base RNA edge candidate need not contain a CDS codon.
		// Preserve all in-bounds legacy values and the unchanged CDS read path.
		startCodon=rna && start>bases.length-3 ? -1 : getCodon(start, bases);
		stopCodon=rna && stop<2 ? -1 : getCodon(stop-2, bases);
		startTruncated=startTruncated_;
		stopTruncated=stopTruncated_;
		type=type_;
		
		if(flip && strand==Shared.MINUS){flip();}
	}

	/** Restores a retained path node without sequence or codon recomputation.
	 * Coordinates are already in the saved current frame; codons and truncation
	 * flags retain their original sense-frame meaning, independently of flip parity.
	 * Component scores and optional reporting attributes are restored separately.
	 */
	public Orf(final String scafName_, final int start_, final int stop_, final int strand_,
			final int frame_, final int scaflen_, final int flipped_, final int type_,
			final int startCodon_, final int stopCodon_, final boolean startTruncated_,
			final boolean stopTruncated_){
		super(scafName_, start_, stop_, strand_, scaflen_, flipped_);
		if(scafName_==null || start_<0 || stop_<start_ || stop_>=scaflen_){
			throw new IllegalArgumentException("Saved path nodes need a named physical contig span: "+
				scafName_+":"+start_+".."+stop_+" length="+scaflen_);
		}
		if(strand_<0 || strand_>1 || frame_<0 || frame_>2 || type_<CDS || type_>RNA){
			throw new IllegalArgumentException("Saved path-list and GFF identity is invalid: strand="+
				strand_+" frame="+frame_+" type="+type_);
		}
		frame=frame_;
		type=type_;
		startCodon=startCodon_;
		stopCodon=stopCodon_;
		startTruncated=startTruncated_;
		stopTruncated=stopTruncated_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Init Helpers         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Extracts a 3-base codon from the sequence starting at the specified position.
	 * Encodes the codon as a packed integer using base-to-number conversion.
	 *
	 * @param from Starting position in the sequence
	 * @param bases Nucleotide sequence array
	 * @return Packed integer representation of the codon
	 */
	private static int getCodon(int from, byte[] bases){
		int codon=0;
		for(int i=0; i<3; i++){
//			assert(i+from<bases.length) : i+", "+from+", "+bases.length;
			byte b=bases[i+from];//Unavailable RNA codons are excluded by the constructor; CDS keeps its ordinary bounds checks.
			int x=AminoAcid.baseToNumber[b];
			codon=(codon<<2)|x;
		}
		return codon;
	}

	public float calcOrfScore(){
		return calcOrfScore(0);
	}

	/** Combines start/stop/inner-kmer sub-scores and ORF length into the heuristic orfScore.
	 * @param overlap Bases overlapping the previous ORF (subtracted as a length penalty).
	 * @return The combined ORF score. */
	public float calcOrfScore(int overlap){
		//Heuristic: geometric-ish blend of start-codon (a), stop-codon (b) and inner-kmer (c) sub-scores,
		//each floored by f1/f2/f3, then scaled by length^2.5 minus an overlap penalty.
		double a=Math.sqrt(Tools.max(f1, e1+startScore));
//		double b=Math.sqrt(f2/*Tools.max(f2, e2+stopScore)*/);//This is better, ignoring stopscore completely
		double b=Math.sqrt(Tools.max(f2, e2+0.35f*stopScore));
		double c=Tools.max(f3, e3+averageKmerScore());
		assert(Double.isFinite(a)) : shared.KillSwitch.assertDie("Nonfinite start subscore: "+a);
		assert(Double.isFinite(b)) : shared.KillSwitch.assertDie("Nonfinite stop subscore: "+b);
		assert(Double.isFinite(c)) : shared.KillSwitch.assertDie("Nonfinite inner subscore: "+c);
		c=4*Math.pow(c, 2.2);
		double d=(0.1*a*b*c*(Math.pow(length()-overlap, 2.5)
			-(overlap<1 ? 0 : Math.pow(overlap+50, 2))));//TODO: Adjust these constants
		if(d>0){d=Math.sqrt(d);}
		assert(Double.isFinite(d)) : shared.KillSwitch.assertDie("Nonfinite ORF score: "+d);
		return (float)d;
	}
	
	/** Adjusts the heuristic orfScore using a neural network prediction.
	 * The NN outputs 0-1 (confidence this is a real gene). Scores above
	 * NN_THRESHOLD boost the heuristic; below penalize it. */
	public float calcOrfScoreHybrid(CellNet net, float gc){
		float heuristic=calcOrfScore();
		float[] in=new float[15];
		in[0]=startScore;
		in[1]=stopScore;
		in[2]=kmerScore;
		in[3]=averageKmerScore();
		in[4]=(float)Math.log10(Tools.max(1, length()));
		in[5]=gc;
		in[6]=(frame==0 ? 1 : 0);
		in[7]=(frame==1 ? 1 : 0);
		in[8]=(frame==2 ? 1 : 0);
		String sc=AminoAcid.codonToString(startCodon);
		in[9]=("ATG".equals(sc) ? 1 : 0);
		in[10]=("GTG".equals(sc) ? 1 : 0);
		in[11]=("TTG".equals(sc) ? 1 : 0);
		String ec=AminoAcid.codonToString(stopCodon);
		in[12]=("TAG".equals(ec) ? 1 : 0);
		in[13]=("TAA".equals(ec) ? 1 : 0);
		in[14]=("TGA".equals(ec) ? 1 : 0);
		net.applyInput(in);
		float nnOut=net.feedForward();
		float adjustment=1.0f+NN_STRENGTH*(nnOut-NN_THRESHOLD);
		return heuristic*Tools.max(0.01f, adjustment);
	}

	/** Adjusts heuristic orfScore using a StartNet (one-hot bases around start site).
	 * Needs access to scaffold bases for one-hot encoding. */
	public float calcOrfScoreHybridStart(CellNet net, float gc, byte[] bases){
		float heuristic=calcOrfScore();
		int numInputs=net.numInputs();
		float[] in=new float[numInputs];
		int idx=0;
		if(strand==Shared.PLUS){
			for(int offset=-START_UPSTREAM; offset<START_DOWNSTREAM; offset++){
				int pos=start+offset;
				int x=(pos>=0 && pos<bases.length) ? AminoAcid.baseToNumber[bases[pos]] : -1;
				in[idx++]=(x==0 ? 1 : 0);
				in[idx++]=(x==1 ? 1 : 0);
				in[idx++]=(x==2 ? 1 : 0);
				in[idx++]=(x==3 ? 1 : 0);
			}
		}else{
			for(int offset=START_UPSTREAM-1; offset>=-START_DOWNSTREAM; offset--){
				int pos=stop+offset;
				int x=(pos>=0 && pos<bases.length) ? AminoAcid.baseToNumber[bases[pos]] : -1;
				if(x>=0){x=3-x;}
				in[idx++]=(x==0 ? 1 : 0);
				in[idx++]=(x==1 ? 1 : 0);
				in[idx++]=(x==2 ? 1 : 0);
				in[idx++]=(x==3 ? 1 : 0);
			}
		}
		if(idx<numInputs){
			in[idx++]=startScore;
			in[idx++]=gc;
			in[idx++]=(float)(0.1*Math.log(Tools.max(1, length())));
		}
		net.applyInput(in);
		float nnOut=net.feedForward();
		float threshold=(NN_THRESHOLD>=0 ? NN_THRESHOLD : net.cutoff);
		if(threshold<0){threshold=0;}
		float adjustment=1.0f+NN_STRENGTH*(nnOut-threshold);
		return heuristic*Tools.max(0.01f, adjustment);
	}

	/** Adjusts heuristic orfScore using a StopNet (one-hot bases around stop site). */
	public float calcOrfScoreHybridStop(CellNet net, float gc, byte[] bases){
		float heuristic=calcOrfScore();
		int numInputs=net.numInputs();
		float[] in=new float[numInputs];
		int idx=0;
		int stopStart=stop-2;
		if(strand==Shared.PLUS){
			for(int offset=-STOP_UPSTREAM; offset<STOP_DOWNSTREAM; offset++){
				int pos=stopStart+offset;
				int x=(pos>=0 && pos<bases.length) ? AminoAcid.baseToNumber[bases[pos]] : -1;
				in[idx++]=(x==0?1:0); in[idx++]=(x==1?1:0); in[idx++]=(x==2?1:0); in[idx++]=(x==3?1:0);
			}
		}else{
			int stopAnchor=start+2;
			for(int offset=STOP_UPSTREAM-1; offset>=-STOP_DOWNSTREAM; offset--){
				int pos=stopAnchor+offset;
				int x=(pos>=0 && pos<bases.length) ? AminoAcid.baseToNumber[bases[pos]] : -1;
				if(x>=0){x=3-x;}
				in[idx++]=(x==0?1:0); in[idx++]=(x==1?1:0); in[idx++]=(x==2?1:0); in[idx++]=(x==3?1:0);
			}
		}
		if(idx<numInputs){
			in[idx++]=stopScore;
			in[idx++]=gc;
			in[idx++]=(float)(0.1*Math.log(Tools.max(1, length())));
		}
		net.applyInput(in);
		float nnOut=net.feedForward();
		float threshold=(NN_STOP_THRESHOLD>=0 ? NN_STOP_THRESHOLD : net.cutoff);
		if(threshold<0){threshold=0;}
		float adjustment=1.0f+NN_STOP_STRENGTH*(nnOut-threshold);
		return heuristic*Tools.max(0.01f, adjustment);
	}

	/** Adjusts using both StartNet and StopNet. */
	public float calcOrfScoreHybridBoth(CellNet startNet, CellNet stopNet, float gc, byte[] bases){
		float startAdj=calcOrfScoreHybridStart(startNet, gc, bases);
		float heuristic=calcOrfScore();
		float ratio=startAdj/Tools.max(0.001f, heuristic);

		int numInputs=stopNet.numInputs();
		float[] in=new float[numInputs];
		int idx=0;
		int stopStart=stop-2;
		if(strand==Shared.PLUS){
			for(int offset=-STOP_UPSTREAM; offset<STOP_DOWNSTREAM; offset++){
				int pos=stopStart+offset;
				int x=(pos>=0 && pos<bases.length) ? AminoAcid.baseToNumber[bases[pos]] : -1;
				in[idx++]=(x==0?1:0); in[idx++]=(x==1?1:0); in[idx++]=(x==2?1:0); in[idx++]=(x==3?1:0);
			}
		}else{
			int stopAnchor=start+2;
			for(int offset=STOP_UPSTREAM-1; offset>=-STOP_DOWNSTREAM; offset--){
				int pos=stopAnchor+offset;
				int x=(pos>=0 && pos<bases.length) ? AminoAcid.baseToNumber[bases[pos]] : -1;
				if(x>=0){x=3-x;}
				in[idx++]=(x==0?1:0); in[idx++]=(x==1?1:0); in[idx++]=(x==2?1:0); in[idx++]=(x==3?1:0);
			}
		}
		if(idx<numInputs){
			in[idx++]=stopScore;
			in[idx++]=gc;
			in[idx++]=(float)(0.1*Math.log(Tools.max(1, length())));
		}
		stopNet.applyInput(in);
		float nnOut=stopNet.feedForward();
		float threshold=(NN_STOP_THRESHOLD>=0 ? NN_STOP_THRESHOLD : stopNet.cutoff);
		if(threshold<0){threshold=0;}
		float stopAdj=1.0f+NN_STOP_STRENGTH*(nnOut-threshold);
		return startAdj*Tools.max(0.01f, stopAdj);
	}

	/** Adjusts heuristic orfScore using a MetaNet that combines StartNet+StopNet outputs
	 * with heuristic scores, codon one-hots, length, GC, and leader distance. */
	public float calcOrfScoreHybridMeta(CellNet startNet, CellNet stopNet, CellNet metaNet,
			float gc, byte[] bases, int maxLenForStop){
		float heuristic=calcOrfScore();

		float startNetOut=runSubNet(startNet, true, gc, bases);
		float stopNetOut=runSubNet(stopNet, false, gc, bases);

		int numInputs=metaNet.numInputs();
		float[] in=new float[numInputs];
		int idx=0;
		in[idx++]=startNetOut;
		in[idx++]=stopNetOut;
		in[idx++]=startScore;
		in[idx++]=stopScore;
		in[idx++]=kmerScore;
		in[idx++]=averageKmerScore();
		in[idx++]=(float)Math.log10(Tools.max(1, length()));
		in[idx++]=gc;

		if(numInputs>10){
			for(int pos=2; pos>=0; pos--){
				int base=(startCodon>>(pos*2))&3;
				in[idx++]=(base==0?1:0); in[idx++]=(base==1?1:0);
				in[idx++]=(base==2?1:0); in[idx++]=(base==3?1:0);
			}
			for(int pos=2; pos>=0; pos--){
				int base=(stopCodon>>(pos*2))&3;
				in[idx++]=(base==0?1:0); in[idx++]=(base==1?1:0);
				in[idx++]=(base==2?1:0); in[idx++]=(base==3?1:0);
			}
		}

		float leaderDist=(maxLenForStop>0 ? maxLenForStop : length());
		in[idx++]=(float)Math.log10(Tools.max(1, leaderDist));
		in[idx++]=length()/(float)Tools.max(1, leaderDist);

		metaNet.applyInput(in);
		float nnOut=metaNet.feedForward();
		float threshold=(NN_META_THRESHOLD>=0 ? NN_META_THRESHOLD : metaNet.cutoff);
		if(threshold<0){threshold=0;}
		float adjustment=1.0f+NN_META_STRENGTH*(nnOut-threshold);
		return heuristic*Tools.max(0.01f, adjustment);
	}

	private float runSubNet(CellNet net, boolean isStart, float gc, byte[] bases){
		int numInputs=net.numInputs();
		float[] in=new float[numInputs];
		int idx=0;
		if(isStart){
			if(strand==Shared.PLUS){
				for(int offset=-START_UPSTREAM; offset<START_DOWNSTREAM; offset++){
					int pos=start+offset;
					int x=(pos>=0 && pos<bases.length) ? AminoAcid.baseToNumber[bases[pos]] : -1;
					in[idx++]=(x==0?1:0); in[idx++]=(x==1?1:0); in[idx++]=(x==2?1:0); in[idx++]=(x==3?1:0);
				}
			}else{
				for(int offset=START_UPSTREAM-1; offset>=-START_DOWNSTREAM; offset--){
					int pos=stop+offset;
					int x=(pos>=0 && pos<bases.length) ? AminoAcid.baseToNumber[bases[pos]] : -1;
					if(x>=0){x=3-x;}
					in[idx++]=(x==0?1:0); in[idx++]=(x==1?1:0); in[idx++]=(x==2?1:0); in[idx++]=(x==3?1:0);
				}
			}
			if(idx<numInputs){in[idx++]=startScore; in[idx++]=gc; in[idx++]=(float)(0.1*Math.log(Tools.max(1, length())));}
		}else{
			int stopStart=stop-2;
			if(strand==Shared.PLUS){
				for(int offset=-STOP_UPSTREAM; offset<STOP_DOWNSTREAM; offset++){
					int pos=stopStart+offset;
					int x=(pos>=0 && pos<bases.length) ? AminoAcid.baseToNumber[bases[pos]] : -1;
					in[idx++]=(x==0?1:0); in[idx++]=(x==1?1:0); in[idx++]=(x==2?1:0); in[idx++]=(x==3?1:0);
				}
			}else{
				int stopAnchor=start+2;
				for(int offset=STOP_UPSTREAM-1; offset>=-STOP_DOWNSTREAM; offset--){
					int pos=stopAnchor+offset;
					int x=(pos>=0 && pos<bases.length) ? AminoAcid.baseToNumber[bases[pos]] : -1;
					if(x>=0){x=3-x;}
					in[idx++]=(x==0?1:0); in[idx++]=(x==1?1:0); in[idx++]=(x==2?1:0); in[idx++]=(x==3?1:0);
				}
			}
			if(idx<numInputs){in[idx++]=stopScore; in[idx++]=gc; in[idx++]=(float)(0.1*Math.log(Tools.max(1, length())));}
		}
		net.applyInput(in);
		return net.feedForward();
	}

	static int START_UPSTREAM=21;
	static int START_DOWNSTREAM=9;
	static int STOP_UPSTREAM=9;
	static int STOP_DOWNSTREAM=21;
	static float NN_STRENGTH=2.0f;
	/** Negative = use net.cutoff (trained optimal crossover). Override with nnthreshold=X. */
	static float NN_THRESHOLD=-1f;
	static float NN_STOP_STRENGTH=0.2f;
	static float NN_STOP_THRESHOLD=-1f;
	static float NN_META_STRENGTH=0.05f;
	static float NN_META_THRESHOLD=-1f;

	/** Mean inner-kmer score: total kmerScore divided by the number of inner kmers spanned by the CDS.
	 * @return kmerScore normalized by (length() - kInnerCDS - 2). */
	public float averageKmerScore(){
		final int innerKmers=length()-GeneModel.kInnerCDS-2;
		if(innerKmers<=0){throw new IllegalStateException("Cannot score inner kmers in a "+length()+" bp ORF");}
		return kmerScore/innerKmers; //This slightly affects score if kInnerCDS is changed
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Whether prev is a valid earlier ORF that this one may chain after, during path-finding.
	 * @param prev Candidate predecessor ORF. @param maxOverlap Max bases prev may overlap into this ORF.
	 * @return True iff prev legally precedes this ORF. */
	public boolean isValidPrev(Orf prev, int maxOverlap){
		assert(prev!=null) : "GeneCaller path lists must supply a real predecessor for an edge";
		//valid only if prev strictly precedes: ends before this ends, overlaps by <=maxOverlap, and starts earlier
		if(prev.stop>=stop || prev.stop>=start+maxOverlap || prev.start>=start){return false;}
		//Reading-frame exclusion applies only to coding genes. NcrnaScavenger
		//assigns RNA a placeholder frame0; applying this rule to RNA dropped two
		//verified Tsukamurella SRPs after endpoint refinement introduced2bp CDS overlap.
		if(type==CDS && prev.type==CDS && prev.frame==frame && prev.strand==strand && prev.stop>=start){return false;}
		return true;
	}

	//path-finding accessors below: the prevStrand arg uses 0=plus, anything else=minus (matches Shared.PLUS=0)
	public float pathScore() {return Tools.max(pathScorePlus, pathScoreMinus);}
	public float pathScore(int prevStrand) {return prevStrand==0 ? pathScorePlus : pathScoreMinus;}

	public Orf prev(){return pathScorePlus>=pathScoreMinus ? prevPlus : prevMinus;}
	public Orf prev(int prevStrand){return prevStrand==0 ? prevPlus : prevMinus;}

	public int pathLength(int prevStrand){return prevStrand==0 ? pathLengthPlus : pathLengthMinus;}
	public int pathLength(){return pathScorePlus>=pathScoreMinus ? pathLengthPlus : pathLengthMinus;}
	
	/*--------------------------------------------------------------*/
	/*----------------           ToString           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Converts ORFs to GFF format lines grouped by feature type.
	 * Creates separate lists for each requested type (CDS, rRNA, tRNA, etc.).
	 *
	 * @param orfs List of ORFs to convert
	 * @param types Comma-separated string of feature types to include
	 * @return Array of GFF line lists, one per requested type
	 */
	public static ArrayList<GffLine>[] toGffLinesByType(ArrayList<Orf> orfs, String types){
		String[] typeArray=types.split(",");
		@SuppressWarnings("unchecked")
		ArrayList<GffLine>[] lists=new ArrayList[typeArray.length];
		HashMap<String, ArrayList<GffLine>> map=new HashMap<String, ArrayList<GffLine>>(3*typeArray.length);
		for(int i=0; i<typeArray.length; i++){
			String type=typeArray[i];
			lists[i]=new ArrayList<GffLine>();
			//TODO: Possible bug [prok/Orf#003] - a duplicate token in `types` makes this map.put overwrite the
			//earlier mapping, so matching ORFs route only to the LAST list for that type and the earlier
			//duplicate's list is returned empty (silent data loss, no warning). Latent LOW: callers pass fixed,
			//non-duplicate type sets.
			map.put(type,  lists[i]);
		}
		for(Orf orf : orfs){
			String type=ProkObject.typeStrings2[orf.type];
			ArrayList<GffLine> glist=map.get(type);
			if(glist!=null) {//null when this ORF's type wasn't among the requested `types` -> silently skip it
				glist.add(new GffLine(orf));
			}
		}
		return lists;
	}
	
	public static ArrayList<GffLine> toGffLines(ArrayList<Orf> orfs){
		if(orfs==null) {return null;}
		ArrayList<GffLine> list=new ArrayList<GffLine>(orfs.size());
		for(Orf orf : orfs){list.add(new GffLine(orf));}
		return list;
	}
	
	public String toStringFlipped(){
		if(strand==flipped()){
			return toString();
		}
		flip();
		String s=toString();
		flip();
		return s;
	}
	
	@Override
	public String toString(){
		return toGff();
	}
	
	public String toGff(){
		ByteBuilder bb=new ByteBuilder();
		appendGff(bb);
		return bb.toString();
	}
	
	//TODO: Probable bug - ByteBuilder.append(char/String) truncates non-ASCII
	//contig/model text instead of UTF-8 encoding it (structures/ByteBuilder#002).
	//Reproduced by NcrnaPolicyIterationTest with an Omega contig: the saved-GFF
	//reader rejects invalid UTF-8. Live query projection rejects these names too.
	//A shared UTF-8 writer change is separate from shortlist iteration.
	public ByteBuilder appendGff(ByteBuilder bb){
		final float innerScore=innerScoreForGff();
		if(scafName==null){
			bb.append('.').tab();
		}else{
			for(int i=0, max=scafName.length(); i<max; i++){
				char c=scafName.charAt(i);
				if(c==' ' || c=='\t'){break;}
				bb.append(c);
			}
			bb.tab();
		}
		bb.append("BBTools").append('\t');
		bb.append(typeStrings2[type]).append('\t');
		bb.append(start+1).append('\t');
		bb.append(stop+1).append('\t');
		
		if(orfScore<0){bb.append('.').append('\t');}
		else{bb.append(orfScore, 2).append('\t');}
		
		bb.append(strand<0 ? '.' : Shared.strandCodes2[strand]).append('\t');
		
		bb.append('0').append('\t');

		//bb.append('.');
		bb.append(typeStrings[type]).append(',');
		if(type==0){
			bb.append("fr").append(frame).append(',');
		}
//		bb.append(startCodon).append(',');
//		bb.append(stopCodon).append(',');
		bb.append("startScr:").append(startScore, 3).append(',');
		bb.append("stopScr:").append(stopScore, 3).append(',');
		bb.append("innerScr:").append(innerScore, 3).append(',');
		bb.append("len:").append(length());
		if(type==0){
			bb.append(',');
			bb.append("start:").append(AminoAcid.codonToString(startCodon)).append(',');
			bb.append("stop:").append(AminoAcid.codonToString(stopCodon));
		}
		if(trnaModel!=null){
			bb.append(",model:").append(trnaModel);
		}
		if(trnaAnticodon!=null){
			bb.append(",anticodon:").append(trnaAnticodon);
		}
		if(ncrnaFamily!=null){
			// Keep the internal registry key stable while displaying the biological name.
			bb.append(",family:").append(NcrnaFamilyRegistry.displayName(ncrnaFamily));
		}
		final int partial=partialMask();
		if(partial!=0){
			bb.append(";partial=").append((partial&1)!=0 ? '1' : '0').append((partial&2)!=0 ? '1' : '0');
		}
		if(proteinEvidence!=null){
			proteinEvidence.append(bb);
			bb.append(";bbgap_candidate_origin=").append(proteinOnlyCandidate ? "rescued_start" : "native");
		}
		if(stopTailEvidence!=null){
			stopTailEvidence.append(bb);
			bb.append(";bbgap_candidate_origin=").append(proteinOnlyCandidate ? "rescued_start" : "native");
		}
		if(tailCascadeAssignment!=null){
			bb.append(";bbgap_mode=tailcascade;bbgap_status=")
				.append(tailCascadeAssignment.assigned() ? "stop_supported" : "no_hit");
			bb.append(";bbgap_tail_stage="); AnnotationGff.escaped(bb, tailCascadeAssignment.terminal);
			if(tailCascadePin!=null){bb.append(";bbgap_reference=").append(tailCascadePin);}
			bb.append(";bbgap_native_score=").append(tailCascadeNativeScore, 6);
			bb.append(";bbgap_call_bonus=").append(proteinBonus, 6);
			bb.append(";bbgap_pruning_protected=").append(tailCascadeProtected ? 't' : 'f');
			if(tailCascadeAssignment.assigned()){
				bb.append(";bbgap_tail_family="); AnnotationGff.escaped(bb, tailCascadeAssignment.familyId);
			}
		}
		return bb;
	}

	private float innerScoreForGff(){
		if(!Float.isFinite(orfScore) || !Float.isFinite(startScore) || !Float.isFinite(stopScore)){
			throw new IllegalStateException("GFF ORF scores must be finite: orf="+orfScore+
				" start="+startScore+" stop="+stopScore+" type="+type+" scaf="+scafName);
		}
		final float innerScore=averageKmerScore();
		if(!Float.isFinite(innerScore)){
			throw new IllegalStateException("GFF inner score must be finite: inner="+innerScore+
				" kmer="+kmerScore+" len="+length()+" type="+type+" scaf="+scafName);
		}
		return innerScore;
	}

	/** Prodigal left/right bits in the same forward-genomic coordinates emitted by appendGff. */
	public int partialMask(){
		assert(ncrnaPartial>=0 && ncrnaPartial<=3) : "Scavenger model-end provenance uses only sense 5-prime/3-prime bits";
		int sense=0;
		if(type==CDS){
			final int first=(strand==Shared.PLUS ? start : scaflen-1-stop);
			final int last=(strand==Shared.PLUS ? stop : scaflen-1-start);
			// GeneCaller.makeOrfsForFrame starts a synthetic edge ORF at exactly
			// startFrame. Its terminal branch ends at the last complete codon;
			// the preserved legacy branch removes one extra base. These exact
			// generator endpoints plus truncation provenance exclude internal gaps.
			final int terminal=scaflen-1-(scaflen-frame)%3;
			if(startTruncated && first==frame){sense|=1;}
			if(stopTruncated && (last==terminal || last==terminal-1)){sense|=2;}
		}else{
			sense=ncrnaPartial;
		}
		final int genomic=(strand==Shared.MINUS ? ((sense&1)<<1)|((sense&2)>>>1) : sense);
		// RNA model clipping can be internal; only an emitted physical-end call
		// gets partial. CDS has the exact frame-bound terminal evidence above.
		return type==CDS ? genomic : genomic&((start==0 ? 1 : 0)|(stop==scaflen-1 ? 2 : 0));
	}
	
	public boolean isSSU(){return type==r16S || type==r18S;}
	public boolean is5S(){return type==r5S;}
	public boolean is16S(){return type==r16S;}
	public boolean is18S(){return type==r18S;}
	public boolean is23S(){return type==r23S;}
	public boolean isCDS(){return type==CDS;}
	public boolean isRRNA(){return type==r18S || type==r16S || type==r5S || type==r23S;}
	public boolean isTRNA(){return type==tRNA;}
	
	/*--------------------------------------------------------------*/
	/*----------------          Overrides           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public float score() {
		return orfScore;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final int frame;

	//These are not needed but nice for printing
	public final int startCodon;
	public final int stopCodon;
	/** True for a synthetic/unknown biological 5' boundary; never swapped by flip(). */
	public final boolean startTruncated;
	/** True for a synthetic/unknown biological 3' boundary, including gaps and contig ends. */
	public final boolean stopTruncated;
	
	public float startScore;
	public float stopScore;
	public float kmerScore;
	
	public float orfScore;
	/** Sum of independent homology stop/ORF contributions, applied once after native scoring.
	 * Evidence retains both components; raw stopScore and NN inputs stay unchanged. */
	public float proteinBonus;
	/** Post-pruning opposite-strand preference; never included in native or support scores. */
	public float tailCascadeOppositePenalty;
	/** True only for starts admitted beyond the native start-score/breakLimit gates. */
	boolean proteinOnlyCandidate;
	ProteinSupport.Evidence proteinEvidence;
	StopTailSupport.Evidence stopTailEvidence;
	/** One immutable classifier decision shared by every native candidate at this stop. */
	prot.FamilyTailCascade.Assignment tailCascadeAssignment;
	String tailCascadePin;
	/** Native score before the flat bonus; meaningful only with tailCascadeAssignment. */
	float tailCascadeNativeScore;
	/** Protects the best native candidate from soft pruning, without inventing a start. */
	boolean tailCascadeProtected;

	//Path scores are for pathfinding phase
	
	public float pathScorePlus;
	public int pathLengthPlus=1;
	public Orf prevPlus;
	
	public float pathScoreMinus;
	public int pathLengthMinus=1;
	public Orf prevMinus;
	
	public final int type;
	public String trnaModel;
	/** Model ends absent beyond physical contig boundaries: bit1=5-prime, bit2=3-prime.
	 * Sense-oriented provenance remains unchanged when genomic coordinates flip. */
	public int ncrnaPartial;
	/** Anticodon triplet extracted structurally from the called sequence, or null */
	public String trnaAnticodon;
	/** Explicit generic-ncRNA family identity (e.g. "s18", "r58", "lsu"), set by
	 * NcrnaScavenger from NcrnaFamily.name; null for non-ncRNA Orfs and for ncRNA Orfs
	 * created before this field existed. NOT a substitute for trnaModel (a model label
	 * is not a safe family identity -- Ganyu/Qiqi design note, 2026-09-17). Purely
	 * additive: unused by appendGff() so far, does not change any existing GFF byte. */
	public String ncrnaFamily;
	/** Selected16S/23S/prok5S verifier identity before endpoint refinement,
	 * or -1 when unavailable; diagnostic provenance, never reconstructed or used by DP. */
	public float acceptedRrnaIdentity=-1f;
	/** Immutable accepted consensus provenance for SSU competition, independent of GFF annotation. */
	public byte[] ssuFullModel;
	public String ssuFullModelName;
	/** Only the euk18S mitochondrial model accepted through modelClipRescue sets this. */
	public boolean ssuClippedMitoRescue;
	
	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/

	/* for kinnercds=6 */ 
//	static float e1=0.1f; 
//	static float e2=-0.04f; 
//	static float e3=0.01f;//Decreasing this decreases TP but increases SNR
//	
//	static float f1=0.08f; 
//	static float f2=0.06f; 
//	static float f3=0.09f;

	/* for kinnercds=7 */ 
	static float e1=0.35f; 
	static float e2=-0.1f; 
	static float e3=-0.01f;//Decreasing this decreases TP but increases SNR
	
	static float f1=0.08f; 
	static float f2=0.02f; 
	static float f3=0.09f;
	
}
