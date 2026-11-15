package illumina;

import java.util.Arrays;

import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

/**
 * Represents one cluster's data across all cycles.
 * Contains bases, quality scores, and physical position.
 *
 * @author Chloe
 * @date October 15, 2025
 */
public class ClusterData{

	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/

	public ClusterData(int lane_, int tile_, int clusterIndex_, float x_, float y_){
		lane=lane_;
		tile=tile_;
		clusterIndex=clusterIndex_;
		//Transform to Illumina FASTQ coordinates: round(10*raw + 1000)
		xIllumina=Math.round(10*x_ + 1000);
		yIllumina=Math.round(10*y_ + 1000);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Set Methods          ----------------*/
	/*--------------------------------------------------------------*/

//	/**
//	 * Set read data for this cluster.
//	 * @param readNum Read number (1=R1, 2=I1, 3=I2, 4=R2)
//	 * @param bases Base calls
//	 * @param quals Quality scores (binned 0-3)
//	 */
//	public void setRead(int readNum, byte[] bases, byte[] quals) {
//		switch(readNum) {
//			case 1: read1=bases; qual1=quals; break;
//			case 2: index1=bases; qualI1=quals; break;
//			case 3: index2=bases; qualI2=quals; break;
//			case 4: read2=bases; qual2=quals; break;
//			default: throw new IllegalArgumentException("Invalid read number: " + readNum);
//		}
//	}

	/**
	 * Set read data using the class-wide default read, index, and UMI term indexes.
	 * @param bases Concatenated base calls for all read-structure terms
	 * @param quals Concatenated quality scores for all read-structure terms
	 * @param lengths Per-term lengths, or null to treat the full arrays as one term
	 */
	public void setData(byte[] bases, byte[] quals, int[] lengths){
		setData(bases, quals, lengths, r1term, i1term, i2term, r2term, u1term, u2term);
	}

	/**
	 * Set read data from selected read-structure terms.
	 * <p>Each term index selects one slot from {@code lengths}. Negative indexes and
	 * indexes outside the split term count are ignored, leaving previous values in
	 * their target fields. When {@code lengths} is null or names one full-length term,
	 * the selected base and quality arrays are borrowed;
	 * otherwise each split term is copied.</p>
	 * @param bases Concatenated base calls for all terms
	 * @param quals Concatenated quality scores for all terms
	 * @param lengths Per-term lengths that sum to the base and quality array lengths
	 * @param r1term_ Term index for read 1, or an absent index
	 * @param i1term_ Term index for index 1, or an absent index
	 * @param i2term_ Term index for index 2, or an absent index
	 * @param r2term_ Term index for read 2, or an absent index
	 * @param u1term_ Term index for UMI 1, or an absent index
	 * @param u2term_ Term index for UMI 2, or an absent index
	 */
	public void setData(byte[] bases, byte[] quals, int[] lengths,
			int r1term_, int i1term_, int i2term_, int r2term_, int u1term_, int u2term_){
		byte[][] splitBases=split(bases, lengths);
		byte[][] splitQuals=split(quals, lengths);
		final int terms=splitBases.length;
		if(r1term_>=0 && r1term_<terms){
			basesR1=splitBases[r1term_];
			qualsR1=splitQuals[r1term_];
		}
		if(i1term_>=0 && i1term_<terms){
			basesI1=splitBases[i1term_];
			qualsI1=splitQuals[i1term_];
		}
		if(i2term_>=0 && i2term_<terms){
			basesI2=splitBases[i2term_];
			qualsI2=splitQuals[i2term_];
		}
		if(r2term_>=0 && r2term_<terms){
			basesR2=splitBases[r2term_];
			qualsR2=splitQuals[r2term_];
		}
		if(u1term_>=0 && u1term_<terms){
			basesU1=splitBases[u1term_];
			qualsU1=splitQuals[u1term_];
		}
		if(u2term_>=0 && u2term_<terms){
			basesU2=splitBases[u2term_];
			qualsU2=splitQuals[u2term_];
		}
	}

	public static byte[][] split(byte[] in, int[] lengths){
		if(lengths==null || (lengths.length==1 && lengths[0]==in.length)){
			return new byte[][] {in};
		}
		assert(in.length==Tools.sum(lengths)) : in.length+"!="+Tools.sum(lengths);
		final int terms=lengths.length;
		byte[][] out=new byte[terms][];
		int prev=0, next=0;
		for(int i=0; i<terms; i++){
			prev=next;
			next+=lengths[i];
			out[i]=Arrays.copyOfRange(in, prev, next);
		}
		return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Output Methods         ----------------*/
	/*--------------------------------------------------------------*/

	public ByteBuilder header(ByteBuilder bb, int rnum){
		if(bb==null) {bb=new ByteBuilder();}
		bb.append(machineID).colon();
		bb.append(runID).colon();
		bb.append(flowcellID).colon();
		bb.append(lane).colon();
		bb.append(tile).colon();
		bb.append(xIllumina).colon();
		bb.append(yIllumina).colon();
		bb.space();
		bb.append(rnum).colon();
		bb.append(passFilter ? 'N' : 'Y').colon();
		bb.append(controlBits).colon();
		appendSubreads(bb, basesI1, basesI2, basesU1, basesU2, 0);
		if(appendBarcodeQuality) {appendSubreads(bb.space(), qualsI1, qualsI2, qualsU1, qualsU2, 33);}
		return bb;
	}

	private static ByteBuilder appendSubreads(ByteBuilder bb, byte[] I1, byte[] I2, byte[] U1, byte[] U2, int offset){
		if(I1!=null) {append(bb, I1, offset);}
		if(I2!=null && I1!=null) {bb.append(barcodeDelimiter);}
		if(I2!=null) {append(bb, I2, offset);}
		if(U1!=null) {append(bb.colon(), U1, offset);}
		if(U2!=null) {append(bb.colon(), U2, offset);}
		return bb;
	}

	private static ByteBuilder append(ByteBuilder bb, byte[] array, int offset){
		for(byte b : array) {bb.append((byte)(b+offset));}
		return bb;
	}

	public Read toRead(){
		Read r1=new Read(basesR1, qualsR1, header(null, 1).toString(), clusterIndex);
		r1.mate=(basesR2==null ? null : new Read(basesR2, qualsR2, header(null, 2).toString(), clusterIndex));
		return r1;
	}

	/**
	 * Convert to FASTQ records.
	 * @param bb Optional destination buffer
	 * @return A FASTQ-format representation of this cluster's called reads
	 */
	public ByteBuilder toBytes(ByteBuilder bb){
		if(bb==null){bb=new ByteBuilder();}

		header(bb.append('@'), 1).nl();
		bb.append(basesR1).nl();
		bb.append('+').nl();
		bb.append(qualsR1).nl();
		if(basesR2==null) {return bb;}

		header(bb.append('@'), 2).nl();
		bb.append(basesR2).nl();
		bb.append('+').nl();
		bb.append(qualsR2).nl();
		return bb;
	}

	/**
	 * Convert to tab-delimited text line.
	 * @param bb Optional destination buffer
	 * @param splitReads If true, output reads comma-delimited (R1,I1,I2,R2); if false, concatenate all
	 * @return A tab-delimited representation of this cluster
	 */
	public ByteBuilder toBytesOld(ByteBuilder bb, boolean splitReads){
		if(bb==null){bb=new ByteBuilder();}

		bb.append(tile).tab();
		bb.append(xIllumina).tab();
		bb.append(yIllumina).tab();
		bb.append(passFilter ? "1" : "0").tab();

		if(splitReads){
			//Comma-delimited reads: R1,I1,I2,R2
			if(basesR1!=null){bb.append(basesR1);} bb.append(',');
			if(basesI1!=null){bb.append(basesI1);} bb.append(',');
			if(basesI2!=null){bb.append(basesI2);} bb.append(',');
			if(basesR2!=null){bb.append(basesR2);}
			bb.tab();

			//Comma-delimited quals
			if(qualsR1!=null){bb.append(qualsR1);} bb.append(',');
			if(qualsI1!=null){bb.append(qualsI1);} bb.append(',');
			if(qualsI2!=null){bb.append(qualsI2);} bb.append(',');
			if(qualsR2!=null){bb.append(qualsR2);}
		}else{
			//Concatenated (default)
			if(basesR1!=null){bb.append(basesR1);}
			if(basesI1!=null){bb.append(basesI1);}
			if(basesI2!=null){bb.append(basesI2);}
			if(basesR2!=null){bb.append(basesR2);}
			bb.tab();

			//Concatenated quals
			if(qualsR1!=null){bb.append(qualsR1);}
			if(qualsI1!=null){bb.append(qualsI1);}
			if(qualsI2!=null){bb.append(qualsI2);}
			if(qualsR2!=null){bb.append(qualsR2);}
		}

		bb.nl();
		return bb;
	}

	@Override
	public String toString(){
		return toBytes(null).toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	public final int tile;
	public final int clusterIndex;
	public final int xIllumina, yIllumina; //Transformed coordinates
	public final int lane;
	public int controlBits=0;
	public boolean passFilter=true;

	public byte[] basesR1, basesI1, basesI2, basesR2, basesU1, basesU2;
	public byte[] qualsR1, qualsI1, qualsI2, qualsR2, qualsU1, qualsU2;

	public static int r1term=0;
	public static int i1term=1;
	public static int i2term=2;
	public static int r2term=3;
	public static int u1term=-1;
	public static int u2term=-1;
	public static String machineID="machine";
	public static String runID="run";
	public static String flowcellID="flowcell";
	public static char barcodeDelimiter='+';
	public static boolean appendBarcodeQuality=false;

}
