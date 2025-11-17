package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Tools;
import structures.ListNum;

/**
 * @author Brian Bushnell
 * @date Apr 3, 2015
 *
 */
public class DualCris extends ConcurrentReadInputStream {
	
	/**
	 * Test driver that demonstrates dual-stream read processing with basic validation.
	 * Creates DualCris from command-line file arguments and processes read lists to
	 * verify pairing functionality.
	 * @param args Input file paths (first required, second optional)
	 */
	public static void main(String[] args){
		String a=args[0];
		String b=args.length>1 ? args[1] : null;
		FileFormat ff1=FileFormat.testInput(a, null, false);
		FileFormat ff2=(b==null ? null : FileFormat.testInput(b, null, false));
		DualCris cris=getReadInputStream(-1, false, ff1, ff2, null, null);
		cris.start();
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=ln.list;
		
		boolean foundR1=false, foundR2=false;
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			for(Read r1 : reads){
				Read r2=r1.mate;
				if(r1.pairnum()==0){foundR1=true;}
				else{foundR2=true;}
				if(r2!=null){
					if(r2.pairnum()==0){foundR1=true;}
					else{foundR2=true;}
				}
			}
			
			System.err.print(ln.id);
			
			cris.returnList(ln.id, foundR1, foundR2);
			foundR1=foundR2=false;
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
			System.err.print(",");
		}
		System.err.print("Finished.");
		cris.returnList(ln.id, foundR1, foundR2);
		ReadWrite.closeStreams(cris);
	}

	/**
	 * Factory method for creating DualCris instances from file formats and quality files.
	 * Constructs underlying ConcurrentReadInputStream objects and wraps them in DualCris
	 * for synchronized paired processing.
	 * @param maxReads Maximum number of reads to process (-1 for unlimited)
	 * @param keepSamHeader Whether to preserve SAM headers
	 * @param ff1 Primary input file format
	 * @param ff2 Secondary input file format (may be null for single-ended)
	 * @param qf1 Quality file for primary input
	 * @param qf2 Quality file for secondary input
	 * @return Configured DualCris instance ready for processing
	 */
	public static DualCris getReadInputStream(long maxReads, boolean keepSamHeader,
			FileFormat ff1, FileFormat ff2, String qf1, String qf2){
		ConcurrentReadInputStream cris1=(ff1==null ? null : ConcurrentReadInputStream.getReadInputStream(maxReads, keepSamHeader, ff1, null, qf1, null));
		ConcurrentReadInputStream cris2=(ff2==null ? null : ConcurrentReadInputStream.getReadInputStream(maxReads, keepSamHeader, ff2, null, qf2, null));
		return new DualCris(cris1, cris2);
	}
	
	/**
	 * Constructs DualCris from two existing ConcurrentReadInputStream objects.
	 * Sets up dual-stream processing with automatic filename concatenation for
	 * identification purposes.
	 * @param cris1_ Primary input stream (may be null)
	 * @param cris2_ Secondary input stream (may be null)
	 */
	public DualCris(ConcurrentReadInputStream cris1_, ConcurrentReadInputStream cris2_){
		super((cris1_==null ? "null" : cris1_.fname())+(cris2_==null ? "null" : ","+cris2_.fname()));
		cris1=cris1_;
		cris2=cris2_;
	}

	/** Primary input stream for R1 reads */
	private final ConcurrentReadInputStream cris1;
	/** Secondary input stream for R2 reads */
	private final ConcurrentReadInputStream cris2;
	/** Whether the secondary input stream is still producing data */
	/** Whether the primary input stream is still producing data */
	private boolean cris1Active, cris2Active;
	/** Tracks error state from either input stream */
	private boolean errorState=false;
	/** Controls verbose output for debugging */
	private boolean verbose=false;
	
	@Override
	public ListNum<Read> nextList() {
		
		ListNum<Read> ln1=null, ln2=null;
		if(cris1Active && cris1!=null){
			ln1=cris1.nextList();
			if(ln1==null){
				synchronized(this){
					cris1Active=false;
					System.err.println("\nSet cris1Active="+cris1Active);
				}
			}
		}
		if(cris2Active && cris2!=null){
			ln2=cris2.nextList();
			if(ln2!=null){
				for(Read r : ln2.list){r.setPairnum(1);}
			}else{
				synchronized(this){
					cris2Active=false;
					System.err.println("\nSet cris2Active="+cris2Active);
				}
			}
		}
		
		if(ln1!=null && ln2!=null){
			final int size1=ln1.size(), size2=ln2.size();
			final int min=Tools.min(size1, size2);
			for(int i=0; i<min; i++){
				Read r1=ln1.get(i);
				Read r2=ln2.get(i);
				r1.mate=r2;
				r2.mate=r1;
			}
			if(size2>size1){
				for(int i=size1; i<size2; i++){
					ln1.add(ln2.get(i));
				}
			}
		}else if(ln2!=null){
			ln1=ln2;
		}
		
		return ln1;
	}
	
	@Override
	public void returnList(long listNum, boolean poison) {
		throw new RuntimeException("Unsupported.");
	}
	
	/**
	 * Returns a processed list to both underlying streams with read-specific feedback.
	 * Deactivates streams when no reads of their type are found, enabling efficient
	 * cleanup when one stream is exhausted before the other.
	 * @param listNum List identifier to return
	 * @param foundR1 Whether R1 reads were found in this list
	 * @param foundR2 Whether R2 reads were found in this list
	 */
	public void returnList(long listNum, boolean foundR1, boolean foundR2) {
		if(cris1!=null && cris1Active){
			cris1.returnList(listNum, !foundR1);
			if(!foundR1){cris1Active=false;}
		}
		if(cris2!=null && cris2Active){
			cris2.returnList(listNum, !foundR2);
			if(!foundR2){cris2Active=false;}
		}
	}
	
	@Override
	public void start() {
		started=true;
		if(cris1!=null){
			cris1.start();
			cris1Active=true;
		}
		if(cris2!=null){
			cris2.start();
			cris2Active=true;
		}
	}
	
	@Override
	public void run() {assert(false);}
	
	@Override
	public void shutdown() {
		if(cris1!=null){cris1.shutdown();}
		if(cris2!=null){cris2.shutdown();}
		cris1Active=cris2Active=false;
	}
	
	@Override
	public void restart() {
		if(cris1!=null){
			cris1.restart();
			cris1Active=true;
		}
		if(cris2!=null){
			cris2.restart();
			cris2Active=true;
		}
	}
	
	@Override
	public void close() {
		if(cris1!=null){cris1.close();}
		if(cris2!=null){cris2.close();}
		cris1Active=cris2Active=false;
	}
	
	@Override
	public boolean paired() {
		assert(cris1!=null || cris2!=null);
		if(cris2!=null){return true;}
		if(cris1!=null){return cris1.paired();}
		return false;
	}
	
	@Override
	public Object[] producers() {
		ArrayList<Object> list=new ArrayList<Object>();
		if(cris1!=null){
			for(Object o : cris1.producers()){list.add(o);}
		}
		if(cris2!=null){
			for(Object o : cris2.producers()){list.add(o);}
		}
		return list.toArray();
	}
	
	@Override
	public boolean errorState() {
		if(cris1!=null){errorState|=cris1.errorState();}
		if(cris2!=null){errorState|=cris2.errorState();}
		return errorState;
	}
	
	@Override
	public void setSampleRate(float rate, long seed) {
		throw new RuntimeException("Invalid.");
	}
	
	@Override
	public long basesIn() {
		return (cris1==null ? 0 : cris1.basesIn())+(cris2==null ? 0 : cris2.basesIn());
	}
	
	@Override
	public long readsIn() {
		return (cris1==null ? 0 : cris1.readsIn())+(cris2==null ? 0 : cris2.readsIn());
	}
	
	@Override
	public boolean verbose() {
		return verbose;
	}

}
