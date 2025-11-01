package stream;

import structures.ListNum;
import java.util.ArrayList;

/**
 * Pairs reads from two separate Streamers (typically R1 and R2 files).
 * Ensures mate references are set correctly.
 * 
 * @author Isla
 * @date October 31, 2025
 */
public class PairStreamer implements Streamer {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public PairStreamer(Streamer s1_, Streamer s2_){
		s1=s1_;
		s2=s2_;
		assert(s1.pairnum()==0) : "First stream must be R1 (pairnum 0)";
		assert(s2.pairnum()==1) : "Second stream must be R2 (pairnum 1)";
		assert(!s1.paired()) : "First stream should not be interleaved";
		assert(!s2.paired()) : "Second stream should not be interleaved";
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void start(){
		s1.start();
		s2.start();
	}
	
	@Override
	public void close(){
		s1.close();
		s2.close();
	}
	
	@Override
	public boolean hasMore(){
		return s1.hasMore();
	}
	
	@Override
	public boolean paired(){return true;}
	
	@Override
	public int pairnum(){return 0;} // Paired data returns as R1
	
	@Override
	public long readsProcessed(){
		return s1.readsProcessed()+s2.readsProcessed();
	}
	
	@Override
	public long basesProcessed(){
		return s1.basesProcessed()+s2.basesProcessed();
	}
	
	@Override
	public ListNum<Read> nextList(){
		ListNum<Read> list1=s1.nextList();
		ListNum<Read> list2=s2.nextList();
		
		if(list1==null && list2==null){return null;}
		
		// Handle mismatched list sizes
		assert(list1!=null && list2!=null) : "Paired files have different read counts!";
		assert(list1.size()==list2.size()) : 
			"List size mismatch: "+list1.size()+" vs "+list2.size();
		
		// Mate the reads
		ArrayList<Read> reads1=list1.list;
		ArrayList<Read> reads2=list2.list;
		for(int i=0; i<reads1.size(); i++){
			Read r1=reads1.get(i);
			Read r2=reads2.get(i);
			assert(r1.numericID==r2.numericID);
			r1.mate=r2;
			r2.mate=r1;
		}
		
		return list1; // Return R1 list (now with mates set)
	}
	
	@Override
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("PairStreamer does not support SamLine");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private final Streamer s1; // R1
	private final Streamer s2; // R2
}