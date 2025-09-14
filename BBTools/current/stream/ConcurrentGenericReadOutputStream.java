package stream;

import java.io.File;
import java.lang.Thread.State;
import java.util.ArrayList;
import java.util.HashMap;

import fileIO.FileFormat;

/**
 * @author Brian Bushnell
 * @date Jan 26, 2015
 *
 */
public final class ConcurrentGenericReadOutputStream extends ConcurrentReadOutputStream {
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	ConcurrentGenericReadOutputStream(FileFormat ff1_, FileFormat ff2_, String qf1, String qf2, int rswBuffers, CharSequence header, boolean useSharedHeader){
		super(ff1_, ff2_);
		
		if(verbose){
			System.err.println("ConcurrentGenericReadOutputStream("+ff1+", "+ff2+", "+qf1+", "+qf2+", "+rswBuffers+", "+useSharedHeader+")");
		}
		
		assert(ff1!=null);
		assert(!ff1.text() && !ff1.unknownFormat()) : "Unknown format for "+ff1;
		
		if(ff1.hasName() && ff1.devnull()){
			File f=new File(ff1.name());
			assert(ff1.overwrite() || !f.exists() || ff1.name().equals("/dev/null")) : f.getAbsolutePath()+" already exists; please delete it.";
			if(ff2!=null){assert(!ff1.name().equals(ff2.name())) : ff1.name()+"=="+ff2.name();}
		}
		
		readstream1=new ReadStreamByteWriter(ff1, qf1, true, rswBuffers, header, useSharedHeader);
		readstream2=ff1.stdio() || ff2==null ? null : new ReadStreamByteWriter(ff2, qf2, false, rswBuffers, header, useSharedHeader);
		
		if(readstream2==null && readstream1!=null){
//			System.out.println("ConcurrentReadOutputStream detected interleaved output.");
			readstream1.OUTPUT_INTERLEAVED=true;
		}
		
		table=(ordered ? new HashMap<Long, ArrayList<Read>>(MAX_CAPACITY) : null);
		
		assert(readstream1==null || readstream1.read1==true);
		assert(readstream2==null || (readstream2.read1==false));
	}
	
	@Override
	public synchronized void start(){
		if(started){
			System.err.println("Resetting output stream.");
			nextListID=0;
			throw new RuntimeException();
		}else{
			started=true;
			if(readstream1!=null){readstream1.start();}
			if(readstream2!=null){readstream2.start();}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public synchronized void add(ArrayList<Read> list, long listnum){
		
		if(ordered){
			int size=table.size();
//			System.err.print(size+", ");
			final boolean flag=(size>=HALF_LIMIT);
			if(listnum>nextListID && size>=ADD_LIMIT){
				if(printBufferNotification){
					System.err.println("Output buffer became full; key "+listnum+" waiting on "+nextListID+".");
					printBufferNotification=false;
				}
				while(listnum>nextListID && size>=HALF_LIMIT){
					try {
						this.wait(20000);
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
					size=table.size();
				}
				if(printBufferNotification){
					System.err.println("Output buffer became clear for key "+listnum+"; next="+nextListID+", size="+size);
				}
			}
			addOrdered(list, listnum);
			assert(listnum!=nextListID);
			if(flag && listnum<nextListID){this.notifyAll();}
		}else{
			addDisordered(list, listnum);
		}
	}
	
	@Override
	public synchronized void close(){
		
		if(table!=null && !table.isEmpty()){
			errorState=true;
			System.err.println("Error: An unfinished ReadOutputStream was closed.");
		}
		//assert(table==null || table.isEmpty()); //Seems like a race condition.  Probably, I should wait at this point until the condition is true before proceeding.
		
//		readstream1.addList(null);
//		if(readstream2!=null){readstream2.addList(null);}
		readstream1.poison();
		if(readstream2!=null){readstream2.poison();}
	}
	
	@Override
	public void join(){
		while(readstream1!=null && readstream1.getState()!=Thread.State.TERMINATED){
			try {
				readstream1.join();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		while(readstream2!=null && readstream2.getState()!=Thread.State.TERMINATED){
			try {
				if(readstream2!=null){readstream2.join();}
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		assert(table==null || table.isEmpty());
		finishedSuccessfully=true;
	}
	
	@Override
	public synchronized void resetNextListID(){
		for(int i=0; i<2000 && !table.isEmpty(); i++){
			try {this.wait(2000);}
			catch (InterruptedException e) {e.printStackTrace();}
		}
		if(!table.isEmpty()){
			System.err.println("WARNING! resetNextListID() waited a long time and the table never cleared.  Process may have stalled.");
		}
		while(!table.isEmpty()){
			try {this.wait(2000);}
			catch (InterruptedException e) {e.printStackTrace();}
		}
		nextListID=0;
	}
	
	@Override
	public final String fname(){
//		if(STANDARD_OUT){return "stdout";}
		return readstream1.fname();
	}
	
	@Override
	public boolean errorState(){
		return errorState || (readstream1!=null && readstream1.errorState()) || (readstream2!=null && readstream2.errorState());
	}
	
	@Override
	public boolean finishedSuccessfully(){
		return finishedSuccessfully && (readstream1==null || readstream1.finishedSuccessfully()) && (readstream2==null || readstream2.finishedSuccessfully());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	
	private synchronized void addOrdered(ArrayList<Read> list, long listnum){
//		System.err.println("RTOS got "+listnum+" of size "+(list==null ? "null" : list.size())+
//				" with first read id "+(list==null || list.isEmpty() || list.get(0)==null ? "null" : ""+list.get(0).numericID));
		assert(list!=null) : listnum;
		assert(listnum>=nextListID) : listnum+", "+nextListID;
//		assert(list.isEmpty() || list.get(0)==null || list.get(0).numericID>=nextReadID) : list.get(0).numericID+", "+nextReadID;
		assert(!table.containsKey(listnum));
		
		table.put(listnum, new ArrayList<Read>(list));
		
		while(table.containsKey(nextListID)){
//			System.err.println("Writing list "+first.get(0).numericID);
			ArrayList<Read> value=table.remove(nextListID);
			write(value);
			nextListID++;
		}
		if(table.isEmpty()){notifyAll();}
	}
	
	private synchronized void addDisordered(ArrayList<Read> list, long listnum){
		assert(list!=null);
		assert(table==null);
		write(new ArrayList<Read>(list));
	}
	
	private synchronized void write(ArrayList<Read> list){
		if(readstream1!=null){
			if(readstream1.getState()==State.TERMINATED){throw new RuntimeException("Writing to a terminated thread.");}
			readstream1.addList(list);
		}
		if(readstream2!=null){
			if(readstream1.getState()==State.TERMINATED){throw new RuntimeException("Writing to a terminated thread.");}
			readstream2.addList(list);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final ReadStreamWriter getRS1(){return readstream1;}
	@Override
	public final ReadStreamWriter getRS2(){return readstream2;}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	private final ReadStreamWriter readstream1;
	private final ReadStreamWriter readstream2;
	private long nextListID=0;
	
	/** Number of lists held before the stream blocks */
	private final int MAX_CAPACITY=256;
	private final int ADD_LIMIT=MAX_CAPACITY-2;
	private final int HALF_LIMIT=ADD_LIMIT/2;
	
	/** For ordered output */
	private final HashMap<Long, ArrayList<Read>> table;
	
	{if(HALF_LIMIT<1){throw new RuntimeException("Capacity too low.");}}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private boolean printBufferNotification=true;
	
}
