package stream;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Set;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Tools;
import structures.ListNum;

/**
 * Allows output of reads to multiple different output streams.
 * @author Brian Bushnell
 * @date Apr 12, 2015
 *
 */
public class MultiCros {
	
	public static void main(String[] args){
		String in=args[0];
		String pattern=args[1];
		ArrayList<String> names=new ArrayList<String>();
		for(int i=2; i<args.length; i++){
			names.add(args[i]);
		}
		final int buff=Tools.max(16, 2*Shared.threads());
		MultiCros mcros=new MultiCros(pattern, null, false, false, false, false, false, FileFormat.FASTQ, buff);
		
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, false, in);
		cris.start();
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		ArrayListSet als=new ArrayListSet(false);
		
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

			for(Read r1 : reads){
				als.add(r1, names);
			}
			cris.returnList(ln);
			if(mcros!=null){mcros.add(als, ln.id);}
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln);
		if(mcros!=null){mcros.add(als, ln.id);}
		ReadWrite.closeStreams(cris);
		ReadWrite.closeStreams(mcros);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public MultiCros(String pattern1_, String pattern2_,
			boolean ordered_, boolean overwrite_, boolean append_, boolean allowSubprocess_, boolean useSharedHeader_, int defaultFormat_, int maxSize_){
		assert(pattern1_!=null && pattern1_.indexOf('%')>=0);
		assert(pattern2_==null || pattern1_.indexOf('%')>=0);
		if(pattern2_==null && pattern1_.indexOf('#')>=0){
			pattern1=pattern1_.replaceFirst("#", "1");
			pattern2=pattern1_.replaceFirst("#", "2");
		}else{
			pattern1=pattern1_;
			pattern2=pattern2_;
		}
		
		ordered=ordered_;
		overwrite= overwrite_;
		append=append_;
		allowSubprocess=allowSubprocess_;
		useSharedHeader=useSharedHeader_;
		
		defaultFormat=defaultFormat_;
		maxSize=maxSize_;

		streamList=new ArrayList<ConcurrentReadOutputStream>();
		streamMap=new LinkedHashMap<String, ConcurrentReadOutputStream>();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public void add(ArrayListSet set, long listnum){
		for(String s : set.getNames()){
			ArrayList<Read> list=set.getAndClear(s);
			if(list!=null){
				add(list, listnum, s);
			}
		}
	}
		
	public void add(ArrayList<Read> list, long listnum, String name){
		ConcurrentReadOutputStream ros=getStream(name);
		ros.add(list, listnum);
	}
	
	public void close(){
		for(ConcurrentReadOutputStream cros : streamList){cros.close();}
	}
	
	public void join(){
		for(ConcurrentReadOutputStream cros : streamList){cros.join();}
	}
	
	public void resetNextListID(){
		for(ConcurrentReadOutputStream cros : streamList){cros.resetNextListID();}
	}
	
	public String fname(){return pattern1;}
	
	/** Return true if this stream has detected an error */
	public boolean errorState(){
		boolean b=errorState;
		for(ConcurrentReadOutputStream cros : streamList){
			b=b&&cros.errorState();
		}
		return b;
	}

	public boolean finishedSuccessfully(){
		boolean b=true;
		for(ConcurrentReadOutputStream cros : streamList){
			b=b&&cros.finishedSuccessfully();
		}
		return b;
	}
	
	public Set<String> getKeys(){return streamMap.keySet();}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	private ConcurrentReadOutputStream makeStream(String name){
		String s1=pattern1.replaceFirst("%", name);
		String s2=pattern2==null ? null : pattern2.replaceFirst("%", name);
		final FileFormat ff1=FileFormat.testOutput(s1, defaultFormat, null, allowSubprocess, overwrite, append, ordered);
		final FileFormat ff2=FileFormat.testOutput(s2, defaultFormat, null, allowSubprocess, overwrite, append, ordered);
		ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ff1, ff2, maxSize, null, useSharedHeader);
		return ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	public ConcurrentReadOutputStream getStream(String name){
		ConcurrentReadOutputStream ros=streamMap.get(name);
		if(ros==null){
			synchronized(streamMap){
				ros=streamMap.get(name);
				if(ros==null){
					ros=makeStream(name);
					ros.start();
					streamList.add(ros);
					streamMap.put(name, ros);
				}
			}
		}
		return ros;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	public final String pattern1, pattern2;
	public final ArrayList<ConcurrentReadOutputStream> streamList;
	public final LinkedHashMap<String, ConcurrentReadOutputStream> streamMap;
	public final boolean ordered;
	
	boolean errorState=false;
	boolean started=false;
	final boolean overwrite;
	final boolean append;
	final boolean allowSubprocess;
	final int defaultFormat;
	final int maxSize;
	final boolean useSharedHeader;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean verbose=false;

}
