package stream;

import java.util.ArrayList;

import structures.ListNum;

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
	public String fname() {return s1.fname()+","+s2.fname();}
	
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
	public void setSampleRate(float rate, long seed){
		if(seed<0) {seed=(long)(Long.MAX_VALUE*Math.random());}
		s1.setSampleRate(rate, seed);//Much faster to pass through than process in this stream
		s2.setSampleRate(rate, seed);
	}
	
	@Override
	public ListNum<Read> nextList(){
		ListNum<Read> ln1=s1.nextList();
		ListNum<Read> ln2=s2.nextList();
		
		if(ln1==null && ln2==null){return null;}
		
		// Handle mismatched list sizes
		assert(ln1!=null && ln2!=null) : "Paired files have different read counts!";
		assert(ln1.size()==ln2.size()) : 
			"List size mismatch: "+ln1.size()+" vs "+ln2.size();
		
		// Mate the reads
		ArrayList<Read> reads1=ln1.list;
		ArrayList<Read> reads2=ln2.list;
		for(int i=0; i<reads1.size(); i++){
			Read r1=reads1.get(i);
			Read r2=reads2.get(i);
			assert(r1.numericID==r2.numericID);
			r1.mate=r2;
			r2.mate=r1;
		}
		
//		// Apply subsampling if needed
//		if(samplerate<1f && randy!=null){
//			int nulled=0;
//			for(int i=0; i<reads1.size(); i++){
//				if(randy.nextFloat()>=samplerate){
//					reads1.set(i, null);
//					nulled++;
//				}
//			}
//			if(nulled>0) {Tools.condenseStrict(reads1);}
//		}
		
		return ln1; // Return R1 list (now with mates set)
	}
	
	@Override
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("PairStreamer does not support SamLine");
	}
	
	@Override
	public boolean errorState(){return s1.errorState() || s2.errorState();}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final Streamer s1; // R1
	private final Streamer s2; // R2
//	private float samplerate=1f;
//	private java.util.Random randy=null;
	
}