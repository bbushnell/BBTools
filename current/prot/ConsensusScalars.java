package prot;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Computes per-family physicochemical scalars from consensus sequences.
 * Input: combined consensus FASTA (one sequence per family, unique IDs).
 * Output: TSV with family_id, hydrophobicity, net_charge, aromaticity,
 * fraction_small, pI_tail, valid_residues.
 *
 * Hydrophobicity: mean Kyte-Doolittle score per valid residue.
 * Net charge: (K+R-D-E) / valid_residues.
 * Aromaticity: (F+W+Y) / valid_residues.
 * Fraction small: (G+A+S) / valid_residues.
 * pI_tail: bisection on Henderson-Hasselbalch (60 iterations). This is a
 *   C-terminal fragment estimate, NOT full-protein pI — the free N-terminus
 *   is an artifact of the tail boundary, not a native terminus. Standard pKa:
 *   D=3.65, E=4.25, C=8.18, Y=10.07, H=6.00, K=10.53, R=12.48,
 *   Nterm=9.69, Cterm=2.34.
 * Valid residues: standard 20 amino acids only; X, B, Z, U, O, J, stop, and gaps skipped.
 * Sequences with zero valid residues emit NA for all scalars and valid_residues=0.
 *
 * @author Chloe
 */
public class ConsensusScalars{

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, ConsensusScalars.class, false);
		try{run(pp.args);}finally{Shared.closeStream(pp.outstream);}
	}

	public static void run(String[] args){
		String inPath=null, outPath=null;
		final Parser parser=new Parser();
		for(String arg : args){
			final String[] pair=arg.split("=", 2);
			if(pair.length!=2){throw new IllegalArgumentException("Require flag=value: "+arg);}
			final String a=pair[0].toLowerCase(), b=pair[1];
			if(a.equals("in")){inPath=b;}
			else if(a.equals("out")){outPath=b;}
			else if(parser.parse(arg, a, b)){/* handled */}
			else{throw new IllegalArgumentException("Unknown flag: "+a);}
		}
		if(inPath==null){inPath=parser.in1;}
		Parser.processQuality();
		if(inPath==null || outPath==null){
			System.err.println("Usage: consensusscalars.sh in=consensus.faa.gz out=scalars.tsv");
			System.exit(1);
		}
		try{
			if(new java.io.File(inPath).getCanonicalPath().equals(new java.io.File(outPath).getCanonicalPath())){
				throw new IllegalArgumentException("Input and output paths must differ");
			}
		}catch(java.io.IOException e){throw new RuntimeException("Cannot resolve canonical paths", e);}
		Shared.AMINO_IN=true;
		Shared.setThreads(1);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;

		final FileFormat ffIn=FileFormat.testInput(inPath, FileFormat.FASTA, null, false, true);
		final FileFormat ffOut=FileFormat.testOutput(outPath, FileFormat.TEXT, null, false, parser.overwrite, parser.append, false);
		final Streamer input=StreamerFactory.makeStreamer(ffIn, null, true, -1, false, true, 1);
		int families=0, naRows=0;
		boolean inputStarted=false;
		try{
			input.start();
			inputStarted=true;
			final ByteStreamWriter writer=new ByteStreamWriter(ffOut);
			boolean writerStarted=false;
			try{
				writer.start();
				writerStarted=true;
				writer.print(HEADER);

				final ByteBuilder bb=new ByteBuilder();
				final double[] scratch=new double[SCALAR_COUNT];
				final java.util.HashSet<String> seen=new java.util.HashSet<String>();
				for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
					for(Read read : batch.list){
						if(read==null){continue;}
						if(read.bases==null){
							throw new IllegalArgumentException("Null bases for read: "+(read.id==null ? "(null ID)" : read.id));
						}
						final String id=firstToken(read.id);
						if(id==null || id.isEmpty()){
							throw new IllegalArgumentException("Empty or missing FASTA ID");
						}
						if(!seen.add(id)){
							throw new IllegalArgumentException("Duplicate family ID: "+id);
						}
						bb.clear();
						if(computeAndAppend(bb, id, read.bases, scratch)){naRows++;}
						writer.print(bb);
						families++;
					}
				}
			}finally{
				if(writerStarted && writer.poisonAndWait()){
					throw new IllegalStateException("ConsensusScalars output write failed");
				}
			}
		}finally{
			if(inputStarted && ReadWrite.closeStream(input)){
				throw new IllegalStateException("ConsensusScalars input read failed");
			}
		}
		System.err.println("CONSENSUS_SCALARS_PASS families="+families+" na_rows="+naRows);
	}

	/** Computes five physicochemical scalars from a protein sequence.
	 * Fills out[offset..offset+4] with hydro, charge, arom, small, pI in that order.
	 * Undefined values are finite zero; use the mask overload when availability matters.
	 * @return valid residue count (0 if no standard amino acids). */
	public static int computeScalars(byte[] seq, int seqLen, double[] out, int offset){
		return (int)(computeScalarState(seq,seqLen,out,offset)>>>SCALAR_COUNT);
	}

	/** Zero-allocation availability API. Bit i in validityMask[0] marks scalar i defined.
	 * A cleared bit always accompanies finite zero, never a NaN comparison or implied gate.
	 * @return valid residue count, preserving the original computeScalars contract. */
	public static int computeScalars(byte[] seq,int seqLen,double[] out,int offset,int[] validityMask){
		if(validityMask==null || validityMask.length<1){throw new IllegalArgumentException("Caller-owned validity mask required");}
		final long state=computeScalarState(seq,seqLen,out,offset);
		validityMask[0]=(int)(state&ALL_DEFINED);
		return (int)(state>>>SCALAR_COUNT);
	}
	/** Mask-only convenience for callers that do not need the residue census; no allocation. */
	public static int computeScalarMask(byte[] seq,int seqLen,double[] out,int offset){
		return (int)(computeScalarState(seq,seqLen,out,offset)&ALL_DEFINED);
	}

	/** Internal packed return: high bits valid residue count, low five bits availability. */
	private static long computeScalarState(byte[] seq,int seqLen,double[] out,int offset){
		if(seq==null || out==null){throw new IllegalArgumentException("Null buffer");}
		if(seqLen<0 || seqLen>seq.length){throw new IllegalArgumentException("seqLen out of range: "+seqLen+"/"+seq.length);}
		if(offset<0 || offset>out.length-SCALAR_COUNT){throw new IllegalArgumentException("Insufficient output span: offset="+offset+" out.length="+out.length);}
		double sumHydro=0;
		int nK=0, nR=0, nD=0, nE=0, nF=0, nW=0, nY=0, nG=0, nA=0, nS=0, nH=0, nC=0;
		int valid=0;

		for(int i=0; i<seqLen; i++){
			final int aa=aaIndex(seq[i]);
			if(aa<0){continue;}
			valid++;
			sumHydro+=KD[aa];
			switch(seq[i]&0xDF){
				case 'K': nK++; break; case 'R': nR++; break;
				case 'D': nD++; break; case 'E': nE++; break;
				case 'F': nF++; break; case 'W': nW++; break;
				case 'Y': nY++; break; case 'G': nG++; break;
				case 'A': nA++; break; case 'S': nS++; break;
				case 'H': nH++; break; case 'C': nC++; break;
			}
		}

		if(valid<1){
			for(int i=0; i<SCALAR_COUNT; i++){out[offset+i]=0;}
			return 0;
		}
		out[offset+HYDRO]=sumHydro/valid;
		out[offset+CHARGE]=(double)(nK+nR-nD-nE)/valid;
		out[offset+AROM]=(double)(nF+nW+nY)/valid;
		out[offset+SMALL]=(double)(nG+nA+nS)/valid;
		final double pi=computePI(nD, nE, nC, nY, nH, nK, nR);
		final boolean defined=pi!=UNDEFINED_PI;
		out[offset+PI]=defined ? pi : 0;
		return ((long)valid<<SCALAR_COUNT)|(defined ? ALL_DEFINED : ALL_DEFINED^(1<<PI));
	}

	/** Convenience overload; allocates a scratch buffer per call. */
	public static boolean computeAndAppend(ByteBuilder bb, String id, byte[] seq){
		return computeAndAppend(bb, id, seq, new double[SCALAR_COUNT]);
	}

	/** Returns true if the row was NA (zero valid residues). Caller-owned scratch avoids per-record allocation. */
	public static boolean computeAndAppend(ByteBuilder bb, String id, byte[] seq, double[] scratch){
		final long state=computeScalarState(seq,seq.length,scratch,0);
		final int valid=(int)(state>>>SCALAR_COUNT),mask=(int)(state&ALL_DEFINED);

		bb.append(id).tab();
		if(valid<1){
			bb.append("NA\tNA\tNA\tNA\tNA\t0").nl();
			return true;
		}
		bb.append(scratch[HYDRO], 6).tab();
		bb.append(scratch[CHARGE], 6).tab();
		bb.append(scratch[AROM], 6).tab();
		bb.append(scratch[SMALL], 6).tab();
		if((mask&(1<<PI))==0){bb.append("NA");}else{bb.append(scratch[PI], 4);}
		bb.tab();
		bb.append(valid).nl();
		return false;
	}

	public static final int SCALAR_COUNT=5;
	public static final int HYDRO=0, CHARGE=1, AROM=2, SMALL=3, PI=4;
	public static final int ALL_DEFINED=(1<<SCALAR_COUNT)-1;
	/** Standalone pI API finite unavailable state, outside its [0,20] root bracket. */
	public static final double UNDEFINED_PI=-1;

	/** Isoelectric point by bisection on Henderson-Hasselbalch net charge.
	 * Bracket [0,20] covers all standard amino acids including polyArg.
	 * Returns UNDEFINED_PI if no root exists within the bracket. */
	public static double computePI(int nD, int nE, int nC, int nY, int nH, int nK, int nR){
		if(nD<0 || nE<0 || nC<0 || nY<0 || nH<0 || nK<0 || nR<0){throw new IllegalArgumentException("Negative residue count");}
		double lo=0, hi=20;
		final double qLo=netCharge(lo, nD, nE, nC, nY, nH, nK, nR);
		final double qHi=netCharge(hi, nD, nE, nC, nY, nH, nK, nR);
		if(qLo<=0 || qHi>=0){return UNDEFINED_PI;}
		for(int i=0; i<60; i++){
			final double pH=(lo+hi)*0.5;
			final double q=netCharge(pH, nD, nE, nC, nY, nH, nK, nR);
			if(q>0){lo=pH;}else{hi=pH;}
		}
		return (lo+hi)*0.5;
	}

	public static double netCharge(double pH, int nD, int nE, int nC, int nY, int nH, int nK, int nR){
		return chargePos(pH, PK_NTERM, 1)+chargePos(pH, PK_K, nK)+chargePos(pH, PK_R, nR)+chargePos(pH, PK_H, nH)
			+chargeNeg(pH, PK_CTERM, 1)+chargeNeg(pH, PK_D, nD)+chargeNeg(pH, PK_E, nE)+chargeNeg(pH, PK_C, nC)+chargeNeg(pH, PK_Y, nY);
	}

	/** Positive charge contribution: n residues with given pKa at given pH (Henderson-Hasselbalch). */
	static double chargePos(double pH, double pKa, int n){
		assert(n>=0) : "Residue count must be nonneg: n="+n;
		return n/(1.0+Math.pow(10, pH-pKa));
	}
	/** Negative charge contribution: n residues with given pKa at given pH (Henderson-Hasselbalch). */
	static double chargeNeg(double pH, double pKa, int n){
		assert(n>=0) : "Residue count must be nonneg: n="+n;
		return -n/(1.0+Math.pow(10, pKa-pH));
	}

	/** Maps standard 20 amino acids to index 0-19; returns -1 for anything else. */
	static int aaIndex(byte b){
		switch(b&0xDF){
			case 'A': return 0;  case 'C': return 1;  case 'D': return 2;  case 'E': return 3;
			case 'F': return 4;  case 'G': return 5;  case 'H': return 6;  case 'I': return 7;
			case 'K': return 8;  case 'L': return 9;  case 'M': return 10; case 'N': return 11;
			case 'P': return 12; case 'Q': return 13; case 'R': return 14; case 'S': return 15;
			case 'T': return 16; case 'V': return 17; case 'W': return 18; case 'Y': return 19;
			default: return -1;
		}
	}

	static String firstToken(String id){
		if(id==null || id.isEmpty()){return null;}
		int end=0;
		while(end<id.length() && !Character.isWhitespace(id.charAt(end))){end++;}
		return end==0 ? null : id.substring(0, end);
	}

	/** Kyte-Doolittle hydrophobicity, indexed by aaIndex (A=0..Y=19). */
	static final double[] KD={
		 1.8,   2.5,  -3.5,  -3.5,   2.8,  // A C D E F
		-0.4,  -3.2,   4.5,  -3.9,   3.8,  // G H I K L
		 1.9,  -3.5,  -1.6,  -3.5,  -4.5,  // M N P Q R
		-0.8,  -0.7,   4.2,  -0.9,  -1.3,  // S T V W Y
	};

	/** Standard pKa values for Henderson-Hasselbalch pI estimation. */
	static final double PK_NTERM=9.69, PK_CTERM=2.34;
	static final double PK_D=3.65, PK_E=4.25, PK_C=8.18;
	static final double PK_Y=10.07, PK_H=6.00, PK_K=10.53, PK_R=12.48;

	static final byte[] HEADER="family_id\thydrophobicity\tnet_charge\taromaticity\tfraction_small\tpI_tail\tvalid_residues\n".getBytes();
}
