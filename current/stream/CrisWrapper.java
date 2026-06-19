package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import structures.ListNum;

/**
 * Wraps a ConcurrentReadInputStream to provide single-read iteration with backtrack capability.
 * Buffers read lists from the underlying stream and manages state for sequential access.
 * Allows going back one read position for reprocessing scenarios.
 *
 * @author Brian Bushnell
 * @date Jul 18, 2014
 */
public class CrisWrapper {
	
	public CrisWrapper(long maxReads, boolean keepSamHeader, FileFormat ff1, FileFormat ff2){
		this(maxReads, keepSamHeader, ff1, ff2, (String)null, (String)null);
	}
	
	public CrisWrapper(long maxReads, boolean keepSamHeader, FileFormat ff1, FileFormat ff2, String qf1, String qf2){
		this(ConcurrentReadInputStream.getReadInputStream(maxReads, ff1.samOrBam(), ff1, ff2, qf1, qf2), true);
	}
		

	public CrisWrapper(ConcurrentReadInputStream cris_, boolean start){
		initialize(cris_, start);
	}
	
	public void initialize(ConcurrentReadInputStream cris_, boolean start){
		cris=cris_;
		if(start){cris.start();}
		ln=cris.nextList();
		reads=(ln==null ? null : ln.list);
		if(reads==null || reads.size()==0){
			reads=null;
			//System.err.println("Empty.");
			//#001-fix [stream/CrisWrapper#001]: guard ln!=null. cris.nextList() can return null (shutdown), giving reads=null here; ln.id would then NPE. When ln==null there is no list to return - just close.
			if(ln!=null){cris.returnList(ln.id, true);}
			errorState|=ReadWrite.closeStream(cris);
		}
		index=0;
		//System.err.println("Initialized.");
	}
	
	public Read next(){
		//System.err.println("*******1");
		Read r=null;
		if(reads==null || index>=reads.size()){
			//System.err.println("*******2");
			if(reads==null){return null;}
			index=0;
			if(reads.size()==0){
				reads=null;
				cris.returnList(ln.id, true);
				errorState|=ReadWrite.closeStream(cris);
				return null;
			}
			cris.returnList(ln.id, false);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
			if(reads==null){
				//System.err.println("*******3");
				//#001-fix [stream/CrisWrapper#001]: ln can be null here (nextList returned null on shutdown); guard before ln.id to avoid NPE. No list to return when ln==null - just close.
				if(ln!=null){cris.returnList(ln.id, true);}
				errorState|=ReadWrite.closeStream(cris);
				//System.err.println("Returning null (2)");
				return null;
			}
		}
		//System.err.println("*******4");
		if(index<reads.size()){
			//System.err.println("*******5");
			r=reads.get(index);
			index++;
		}else{//Near-dead defensive recursion: after a refill the new list is non-empty (the size==0 case already returned null above) with index reset to 0, so index<reads.size() holds and this branch isn't normally reached; bounded to one level if it ever is.
			//System.err.println("*******6");
			//System.err.println("Recalling");
			return next();
		}
		//System.err.println("*******7");
		//System.err.println("Returning "+(r==null ? "null" : r.id));
		return r;
	}
	
	public void goBack(){
		assert(index>0);
		index--;
	}
	
	private ListNum<Read> ln;
	private ArrayList<Read> reads;
	private int index;
	public ConcurrentReadInputStream cris;
	public boolean errorState=false;
	
}
