package clump;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import sort.SortByName;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ListNum;

public class StreamToOutput {
	
	public StreamToOutput(FileFormat ffin1, FileFormat ffin2, ConcurrentReadOutputStream[] rosa_, KmerComparator old, boolean sortByName_, boolean incrementComparator){
		final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, false, ffin1, ffin2, null, null);
		cris.start();
		rosa=rosa_;
		kc=(incrementComparator ? new KmerComparator(old.k, old.seed+1, old.border-1, old.hashes, false, false) : old);
		sortByName=sortByName_;
	}
	
	public StreamToOutput(ConcurrentReadInputStream cris_, ConcurrentReadOutputStream[] rosa_, KmerComparator old, boolean sortByName_, boolean incrementComparator){
		cris=cris_;
		rosa=rosa_;
		kc=(incrementComparator ? new KmerComparator(old.k, old.seed+1, old.border-1, old.hashes, false, false) : old);
		sortByName=sortByName_;
	}
	
	public boolean process(){
		if(rosa==null || rosa.length==0){return errorState;}
		
		File temp=null;
		if(sortByName){
			try {
				temp=File.createTempFile("temp_namesort_", ".fq.gz");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			SortByName sbn=new SortByName(new String[] {"out="+temp.getAbsolutePath()});
			sbn.processInner(cris);
			FileFormat ff=FileFormat.testInput(temp.getAbsolutePath(), null, false);
			cris=ConcurrentReadInputStream.getReadInputStream(-1, false, ff, null, null, null);
		}
		
		if(rosa.length==1){
			processSingle(cris);
		}else{
			processMulti(cris);
		}
		
		errorState|=ReadWrite.closeStream(cris);
		if(temp!=null){
			temp.delete();
		}
		return errorState;
	}
	
	public void processSingle(ConcurrentReadInputStream cris){
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			if(rosa!=null){rosa[0].add(reads, ln.id);}
			
			for(Read r : reads){
				readsIn+=r.pairCount();
				basesIn+=r.pairLength();
			}
			
			cris.returnList(ln);
			
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
	}
	
	public void processMulti(ConcurrentReadInputStream cris){
		final int groups=rosa.length;
		
		@SuppressWarnings("unchecked")
		final ArrayList<Read>[] out=new ArrayList[groups];
		for(int i=0; i<out.length; i++){
			out[i]=new ArrayList<Read>();
		}
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			for(Read r : reads){
				long kmer=kc.hash(r, null, 0, false);
				int group=(int)(kmer%groups);
				out[group].add(r);
				
				readsIn+=r.pairCount();
				basesIn+=r.pairLength();
			}
			for(int group=0; group<groups; group++){
				rosa[group].add(out[group], ln.id);
				out[group]=new ArrayList<Read>();
			}
			
			cris.returnList(ln);
			
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
	}
	
	long readsIn=0;
	long basesIn=0;
	
//	final FileFormat ffin1;
//	final FileFormat ffin2;
	ConcurrentReadInputStream cris;
	ConcurrentReadOutputStream[] rosa;
	final KmerComparator kc;
	final boolean sortByName;
	
	boolean errorState=false;
	
}
