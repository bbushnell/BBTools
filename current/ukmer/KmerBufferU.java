package ukmer;

import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;

/**
 * @author Brian Bushnell
 * @date Jul 9, 2015
 *
 */
public class KmerBufferU {
	
	public KmerBufferU(int buflen, int kbig, boolean initValues){
		k=Kmer.getK(kbig);
		mult=Kmer.getMult(kbig);
		kmers=new LongList(buflen*mult);
		values=(initValues ? new IntList(buflen) : null);
	}
	
	public int add(Kmer kmer){
//		System.err.println("Adding "+kmer+"; this="+this+"; kmers.size="+kmers.size);
		add(kmer.key());
		return kmers.size;
//		System.err.println("Added "+kmer+"; this="+this+"; kmers.size="+kmers.size);
	}
	
	public void add(Kmer kmer, int value){
		add(kmer.key(), value);
	}
	
	public void add(long[] kmer){
		assert(values==null);
		assert(kmer.length==mult) : kmer.length+", "+mult+", "+k;
		kmers.append(kmer);
	}
	
	public void add(long[] kmer, int value){
		assert(kmer.length==mult);
		kmers.append(kmer);
		values.add(value);
		assert(values.size*mult==kmers.size);
	}
	
	public void clear(){
		kmers.clear();
		if(values!=null){values.clear();}
	}
	
	//Returns raw size of kmers array, rather than actual number of kmers
	final int size(){return kmers.size;}
	
	@Override
	public String toString(){
		ByteBuilder bb=new ByteBuilder();
		for(int i=0; i<kmers.size; i++){
			if(i>0){bb.append(',');}
			bb.appendKmer(kmers.get(i), k);
		}
		return bb.toString();
	}
	
	private final int mult;
	private final int k;
	final LongList kmers;
	final IntList values;
	
}
