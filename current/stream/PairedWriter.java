package stream;

import java.util.ArrayList;

import structures.ListNum;

/**
 * Delegates read batches to two Writers, typically for separate R1 and R2 files.
 * Both children receive the same batch and read references, without a defensive
 * copy. The caller configures each child's mate selection; this wrapper neither
 * partitions reads nor changes child configuration. Child calls occur in w1/w2
 * order. Construction does not start either child.
 *
 * @author Isla
 * @date October 31, 2025
 */
public class PairedWriter implements Writer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Stores two nonnull, already configured child writers.
	 * For separate mate files, configure w1 with writeR1=true/writeR2=false
	 * and w2 with writeR1=false/writeR2=true before wrapping them.
	 * @param w1_ First child, normally configured for R1
	 * @param w2_ Second child, normally configured for R2
	 */
	public PairedWriter(Writer w1_, Writer w2_){
		w1=w1_;
		w2=w2_;
		assert(w1!=null && w2!=null);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts the first child, then the second. */
	@Override
	public void start(){
		w1.start();
		w2.start();
	}

	/** Returns the sum of currently reported child read counts without waiting. */
	@Override
	public long readsWritten(){return w1.readsWritten()+w2.readsWritten();}

	/** Returns the sum of currently reported child base counts without waiting. */
	@Override
	public long basesWritten(){return w1.basesWritten()+w2.basesWritten();}

	/**
	 * Wraps the existing list in one ListNum and delegates through addReads.
	 * @param list Read list shared with both children; not copied
	 * @param id Batch number obeying both children's ordering requirements
	 */
	public final void add(ArrayList<Read> list, long id){addReads(new ListNum<Read>(list, id));}

	/**
	 * Passes the identical batch to the first child, then the second.
	 * Child submission and ownership requirements still apply; this wrapper
	 * does not create independent copies for the two writers.
	 * @param reads Shared batch with an ID valid for both children
	 */
	@Override
	public void addReads(ListNum<Read> reads){
		// Children select their configured mates from the same batch.
		w1.addReads(reads);
		w2.addReads(reads);
	}

	/**
	 * Rejects SamLine batches; this wrapper only supports Read batches.
	 * @param lines Unsupported batch
	 * @throws UnsupportedOperationException Always
	 */
	@Override
	public void addLines(ListNum<SamLine> lines){
		throw new UnsupportedOperationException("PairedWriter does not support SamLine");
	}

	/** Signals normal end of input to the first child, then the second. */
	@Override
	public void poison(){
		w1.poison();
		w2.poison();
	}

	/**
	 * Waits for each child in order, then combines their returned error flags.
	 * Both waits occur even if the first returns true.
	 * @return true if either child reports an error
	 */
	@Override
	public boolean waitForFinish(){
		boolean error1=w1.waitForFinish();
		boolean error2=w2.waitForFinish();
		return error1 || error2;
	}

	/** Signals both children, then waits for both; returns their combined error flag. */
	@Override
	public boolean poisonAndWait(){
		poison();
		return waitForFinish();
	}

	/**
	 * Delegates to both children in order. Each child's finishError is required
	 * to be nonblocking by the Writer contract, so these two sequential calls
	 * must be bounded by the calls themselves rather than either child's backlog.
	 * This wrapper relies on that child contract; it does not wait for completion.
	 */
	@Override
	public void finishError(){
		w1.finishError();
		w2.finishError();
	}

	/** Returns whether either child reports an error; skips w2 when w1 returns true. */
	@Override
	public boolean errorState(){return w1.errorState() || w2.errorState();}

	/** Returns whether both children report success; skips w2 when w1 returns false. */
	@Override
	public boolean finishedSuccessfully(){return w1.finishedSuccessfully() && w2.finishedSuccessfully();}

	/** Returns the two child filenames formatted as (name1,name2). */
	@Override
	public final String fname(){return "("+w1.fname()+","+w2.fname()+")";}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** First child; the caller normally configures it to select R1 only. */
	private final Writer w1;
	/** Second child; the caller normally configures it to select R2 only. */
	private final Writer w2;
}
