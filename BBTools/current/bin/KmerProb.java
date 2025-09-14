package bin;

import java.util.ArrayList;

import dna.Data;
import fileIO.ByteFile;
import shared.LineParser1;
import shared.Tools;
import structures.FloatList;

public class KmerProb {
	
	public static synchronized void load() {
		if(matrix!=null) {return;}
		String fname=Data.findPath("?shred4merFractions.tsv");
		matrix=load(fname);
	}
	
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
	static int quantizeLength(long size) {
		size=Tools.mid(size, 200, 200000);
		int idx=(int)Math.ceil(2*Tools.log2(size));
		return idx-idxOffset;
	}
	
	static int dequantizeLength(int idx) {
		idx+=idxOffset;
		int size=(int)Math.pow(2, idx/2f);
		return size;
	}
	
	private static final int idxOffset=(int)Math.ceil(2*Tools.log2(200));
	
	public static float[][] matrix;
	
}
