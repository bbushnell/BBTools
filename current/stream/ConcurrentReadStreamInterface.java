package stream;

import structures.ListNum;

public interface ConcurrentReadStreamInterface extends Runnable{
	
	public void start();
	
	public ListNum<Read> nextList();
	
	public void returnList(ListNum<Read> ln);
	
	public void returnList(long listNumber, boolean poison);
	
	@Override
	public void run();
	
	public void shutdown();

	/** Reset state to allow production of reads from the beginning of the input files.
	 * Does not work with stdin (may cause strange behavior). */
	public void restart();
	
	public void close();
	
	public boolean paired();
	
	public String fname();
	
	public Object[] producers();
	
	public boolean errorState();
	
	public void setSampleRate(float rate, long seed);
	
	public long basesIn();
	
	public long readsIn();
	
	public boolean verbose();
	
//	public String classname(); //e.g. getClass().getName()
	
}
