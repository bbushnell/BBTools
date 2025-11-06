package stream;

import java.util.ArrayList;

import structures.ListNum;

/**
 * Pairs writes to two separate Writers (typically R1 and R2 files).
 * Delegates the same read list to both writers, which filter by pairnum.
 * 
 * @author Isla
 * @date October 31, 2025
 */
public class PairedWriter implements Writer {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Creates a paired writer wrapping two output streams.
	 * w1 should be configured with writeR1=true, writeR2=false.
	 * w2 should be configured with writeR1=false, writeR2=true.
	 * 
	 * @param w1_ Writer for R1 reads
	 * @param w2_ Writer for R2 reads
	 */
	public PairedWriter(Writer w1_, Writer w2_){
		w1=w1_;
		w2=w2_;
		assert(w1!=null && w2!=null);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void start(){
		w1.start();
		w2.start();
	}
	
	@Override
	public void close(){
		w1.close();
		w2.close();
	}
	
	@Override
	public long readsWritten(){
		return w1.readsWritten()+w2.readsWritten();
	}
	
	@Override
	public long basesWritten(){
		return w1.basesWritten()+w2.basesWritten();
	}
	
	public final void add(ArrayList<Read> list, long id) {addReads(new ListNum<Read>(list, id));}
	
	@Override
	public void addReads(ListNum<Read> reads){
		// Both writers get the same list; they filter internally by pairnum
		w1.addReads(reads);
		w2.addReads(reads);
	}
	
	@Override
	public void addLines(ListNum<SamLine> lines){
		throw new UnsupportedOperationException("PairedWriter does not support SamLine");
	}
	
	@Override
	public void poison(){
		w1.poison();
		w2.poison();
	}
	
	@Override
	public boolean waitForFinish(){
		boolean error1=w1.waitForFinish();
		boolean error2=w2.waitForFinish();
		return error1 || error2;
	}
	
	@Override
	public boolean poisonAndWait(){
		poison();
		return waitForFinish();
	}
	
	@Override
	public boolean errorState(){return w1.errorState() || w2.errorState();}
	
	@Override
	public boolean finishedSuccessfully() {return w1.finishedSuccessfully() && w2.finishedSuccessfully();}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Writer for R1 reads (writeR1=true, writeR2=false) */
	private final Writer w1;
	/** Writer for R2 reads (writeR1=false, writeR2=true) */
	private final Writer w2;
}