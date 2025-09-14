package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import structures.ListNum;

public abstract class ReadInputStream {
	
	public static final ArrayList<Read> toReads(String fname, int defaultFormat, long maxReads){
		if(fname==null){return null;}
		FileFormat ff=FileFormat.testInput(fname, defaultFormat, null, false, true);
		return toReads(ff, maxReads);
	}
	
	public static final Read[] toReadArray(FileFormat ff, long maxReads){
		ArrayList<Read> list=toReads(ff, maxReads);
		return list==null ? null : list.toArray(new Read[0]);
	}
	
	public static final ArrayList<Read> toReads(FileFormat ff, long maxReads){
		ArrayList<Read> list=new ArrayList<Read>();

		/* Start an input stream */
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ff, null);
		cris.start(); //4567
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		/* Iterate through read lists from the input stream */
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			list.addAll(reads);
			
			/* Dispose of the old list and fetch a new one */
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		/* Cleanup */
		cris.returnList(ln);
		ReadWrite.closeStream(cris);
		return list;
	}
	

	public abstract Read next();
	
//	public final ArrayList<Read> fetchAll(){
//		ArrayList<Read> out=new ArrayList<Read>();
//		for(ArrayList<Read> list=nextList(); list!=null && list.size()>0; list=nextList()){
//			out.addAll(list);
//		}
//		close();
//		return out;
//	}
	
	public abstract ArrayList<Read> nextList();
	
	public abstract boolean hasMore();

	public abstract void restart();
	
	/** Returns true if there was an error, false otherwise */
	public abstract boolean close();

	public abstract boolean paired();

	protected static final ArrayList<Read> toList(Read[] array){
		if(array==null || array.length==0){return null;}
		ArrayList<Read> list=new ArrayList<Read>(array.length);
		for(int i=0; i<array.length; i++){list.add(array[i]);}
		return list;
	}
	
	/** Return true if this stream has detected an error */
	public boolean errorState(){return errorState;}
	/** TODO */
	protected boolean errorState=false;
	
	public final boolean preferLists(){return true;}

	public abstract void start();
	
	public abstract String fname();
	
}
