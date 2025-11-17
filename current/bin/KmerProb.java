package bin;

import java.util.ArrayList;

import dna.Data;
import fileIO.ByteFile;
import shared.LineParser1;
import shared.Tools;
import structures.FloatList;

/**
 * Calculates probability that two sequences come from the same genome based on
 * k-mer frequency cosine differences. Uses precomputed matrices derived from
 * 5000 bacterial genomes with <95% ANI across ~700M sequence pairs.
 * @author Brian Bushnell
 */
public class KmerProb {
	
	/** Loads the k-mer probability matrix from shred4merFractions.tsv file.
	 * Thread-safe lazy initialization - only loads once. */
	public static synchronized void load() {
		if(matrix!=null) {return;}
		String fname=Data.findPath("?shred4merFractions.tsv");
		matrix=load(fname);
	}
	
	/**
	 * Loads probability matrix from specified TSV file.
	 * Parses tab-delimited data, skipping comment lines starting with '#'.
	 * Each row becomes a probability vector for a sequence length bin.
	 * @param fname Path to the TSV file containing probability data
	 * @return 2D array where [length_bin][difference_index] = probability
	 */
	private static float[][] load(String fname){
		ArrayList<byte[]> lines=ByteFile.toLines(fname);
		ArrayList<FloatList> vectors=new ArrayList<FloatList>();
		LineParser1 lp=new LineParser1('\t');
		for(byte[] line : lines) {
			lp.set(line);
			if(!lp.startsWith('#')) {
//				System.err.println("Making a vector: "+lp.terms());
				FloatList list=new FloatList(lp.terms()-2);
				for(int i=2; i<lp.terms(); i++) {
					float f=lp.parseFloat(i);
					assert(f>=0 && f<=1);
					list.add(f);
				}
				vectors.add(list);
			}
		}
		float[][] data=new float[vectors.size()][];
		for(int i=0; i<vectors.size(); i++) {
			FloatList list=vectors.get(i);
			assert(list.size()==list.array.length);
			list.shrink();
			data[i]=list.array;
		}
//		System.err.println("Read "+data.length+" vectors.");
		return data;
	}
	
	/** Probability that two sequences with this 
	 * kmer frequency cosine difference come from the same genome.
	 * The length is the length of the shorter sequence.
	 * Genomes used for this were 5000 bacteria, <95% ANI,
	 * around 700m pairs.
	 * @param length
	 * @param dif
	 * @return
	 */
	public static float prob(long length, float dif) {
		int idx=quantizeLength(length);
		float[] array=matrix[idx];
		int idx2=Math.min(array.length-1, (int)(dif*1024));
		return array[idx2];
	}
	
	//Bins contain everything UP TO that size.
	//For example, bin 0 is 256, which contains 1-256.
	//Bin 1 is 256-362, then 363-512, etc.
	/**
	 * Maps sequence length to matrix bin index using log2 scaling.
	 * Clamps size to 200-200000 range then applies logarithmic binning.
	 * @param size Sequence length to quantize
	 * @return Matrix row index for this length bin
	 */
	static int quantizeLength(long size) {
		size=Tools.mid(size, 200, 200000);
		int idx=(int)Math.ceil(2*Tools.log2(size));
		return idx-idxOffset;
	}
	
	/**
	 * Converts matrix bin index back to representative sequence length.
	 * Inverse operation of quantizeLength().
	 * @param idx Matrix row index
	 * @return Representative length for this bin
	 */
	static int dequantizeLength(int idx) {
		idx+=idxOffset;
		int size=(int)Math.pow(2, idx/2f);
		return size;
	}
	
	/** Offset for length quantization, precomputed from minimum length (200) */
	private static final int idxOffset=(int)Math.ceil(2*Tools.log2(200));
	
	/** Probability matrix where [length_bin][difference_index] = probability */
	public static float[][] matrix;
	
}
