package stream;

import java.util.ArrayList;
import java.util.Comparator;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import sort.ReadComparatorClump;
import sort.ReadComparatorTopological5Bit;
import structures.ListNum;

public class CrisContainer implements Comparable<CrisContainer> {
	
	public CrisContainer(String fname, Comparator<Read> comparator_, boolean allowSubprocess){
		comparator=comparator_;
		genKmer=(comparator==ReadComparatorTopological5Bit.comparator);
		clump=(comparator==ReadComparatorClump.comparator);
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, allowSubprocess, true);
		cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ff, null, null, null);
//		System.err.println(genKmer+", "+clump+", "+comparator.getClass());
		cris.start();
		fetch();
	}
	
	public CrisContainer(ConcurrentReadInputStream cris_, Comparator<Read> comparator_){
		comparator=comparator_;
		genKmer=(comparator==ReadComparatorTopological5Bit.comparator);
		clump=(comparator==ReadComparatorClump.comparator);
		cris=cris_;
//		System.err.println(genKmer+", "+clump+", "+comparator.getClass());
		fetch();
	}
	
	public ArrayList<Read> fetch(){
		final ArrayList<Read> old=list;
		fetchInner();
		return old;
	}
	
	private void fetchInner(){
		ListNum<Read> ln=cris.nextList();
		list=(ln==null ? null : ln.list);
		if(list.size()<1){list=null;}
		else if(genKmer){
			for(Read r : list){ReadComparatorTopological5Bit.genKmer(r);}
		}else if(clump){
			for(Read r : list){ReadComparatorClump.set(r);}
		}
		read=(list==null ? null : list.get(0));
		if(lastNum>=0){cris.returnList(lastNum, list==null);}
		if(ln!=null){lastNum=ln.id;}
		assert((read==null)==(list==null || list.size()==0));
//		if(count>0 && list!=null){
//			for(Read r : list){
//				assert(remainingReads>=0) : remainingReads+", "+count+", "+r.numericID;
//				double remaining=(count-sum);
//				double mult=2*(remaining/remainingReads);
//				sum=sum+randy.nextDouble()*mult;
//				r.rand=sum;
////				System.err.println(r.rand);
//				remainingReads--;
//			}
//		}
	}
	
	public boolean close(){
		return ReadWrite.closeStream(cris);
	}
	
	public Read peek(){return read;}
	
	@Override
	public int compareTo(CrisContainer other) {
		assert(read!=null);
		assert(other.read!=null);
		return comparator.compare(read, other.read);
	}
	
	public int compareTo(Read other) {
		return comparator.compare(read, other);
	}
	
	public boolean hasMore(){
		return read!=null;
	}
	
	public ConcurrentReadInputStream cris(){return cris;}
	
	final ConcurrentReadInputStream cris;
	private Read read;
	private long lastNum=-1;
	private ArrayList<Read> list;
	private final Comparator<Read> comparator;
	private final boolean genKmer, clump;
//	private double sum=0;
//	final int count;
//	private final Random randy;
//	private int remainingReads;
	
}
