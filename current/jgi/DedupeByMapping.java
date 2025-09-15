package jgi;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.PriorityQueue;

import shared.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.ReadStreamWriter;
import stream.SamLine;
import structures.ListNum;
import template.BBTool_ST;

/**
 * @author Brian Bushnell
 * @date Jan 30, 2015
 *
 */
public class DedupeByMapping extends BBTool_ST{
	
	public static void main(String[] args){
		Timer t=new Timer();
		DedupeByMapping bbt=new DedupeByMapping(args);
		bbt.process(t);
	}

	/**
	 * @param args
	 */
	public DedupeByMapping(String[] args) {
		super(args);
		reparse(args);
		SamLine.SET_FROM_OK=true;
		ReadStreamWriter.USE_ATTACHED_SAMLINE=true;
		if(sorted){queue=new PriorityQueue<Quad>(initialSize);}
	}

	/* (non-Javadoc)
	 * @see jgi.BBTool_ST#setDefaults()
	 */
	@Override
	protected void setDefaults() {
		keepUnmapped=true;
		keepSingletons=true;
		sorted=false;
		usePairOrder=true;
	}
	
	@Override
	public boolean parseArgument(String arg, String a, String b) {
		if(a.equals("keepunmapped") | a.equals("ku")){
			keepUnmapped=Parse.parseBoolean(b);
			return true;
		}else if(a.equals("keepsingletons") | a.equals("ks")){
			keepSingletons=Parse.parseBoolean(b);
			return true;
		}else if(a.equals("ignorepairorder") | a.equals("ipo")){
			usePairOrder=!Parse.parseBoolean(b);
			return true;
		}else if(a.equals("sorted")){
			sorted=Parse.parseBoolean(b);
			return true;
		}
		return false;
	}
	
	@Override
	protected final boolean useSharedHeader(){return true;}
	
	@Override
	protected boolean processReadPair(Read r1, Read r2) {
		assert(r2==null);
		return (sorted ? processReadPair_sorted(r1) : processReadPair_unsorted(r1));
	}
	
	boolean processReadPair_unsorted(Read r1) {
		SamLine sl=r1.samline;
		if(!sl.primary()){return false;}
		if(sl.mapped()){
			String rname=new String(sl.rname());
			Integer x=contigToNumber.get(rname);
			if(x==null){
				x=contigToNumber.size();
				contigToNumber.put(rname, x);
			}
			r1.chrom=x;
			r1.start=sl.start(true, false);
			r1.stop=sl.stop(r1.start, true, false);
			r1.setStrand(sl.strand());
		}else{
			r1.chrom=-1;
			r1.start=-1;
		}
		
		Read old=nameToRead.get(r1.id);
		if(old==null){
			nameToRead.put(r1.id, r1);
		}else{
			//assert(old.mate==null);//This triggers...  maybe, after filtering and leaving some unpaired reads
			old.mate=r1;
			r1.mate=old;
			SamLine sl2=old.samline;
			if(sl2.pairnum()==1){
				nameToRead.put(r1.id, r1);
			}
		}
		return true;
	}
	
	
	boolean processReadPair_sorted(Read r1) {
		assert(false) : "TODO";
		SamLine sl=r1.samline;
		if(!sl.primary()){return false;}
		if(sl.mapped()){
			String rname=new String(sl.rname());
			Integer x=contigToNumber.get(rname);
			if(x==null){
				x=contigToNumber.size();
				contigToNumber.put(rname, x);
			}
			r1.chrom=x;
			r1.start=sl.start(true, false);
			r1.stop=sl.stop(r1.start, true, false);
			r1.setStrand(sl.strand());
		}else{
			r1.chrom=-1;
			r1.start=-1;
		}
		
		Read old=nameToRead.get(r1.id);
		if(old==null){
			nameToRead.put(r1.id, r1);
		}else{
			assert(old.mate==null);
			old.mate=r1;
			r1.mate=old;
			SamLine sl2=old.samline;
			if(sl2.pairnum()==1){
				nameToRead.put(r1.id, r1);
			}
		}
		return true;
	}
	
	@Override
	protected void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		if(sorted){processInner_sorted(cris, ros);}
		else{processInner_unsorted(cris, ros);}
	}
	
	
	void processInner_unsorted(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		
		readsProcessed=0;
		basesProcessed=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					assert(r1.mate==null);
					assert(r1.samline!=null);
					
					final int initialLength1=r1.length();
					
					{
						readsProcessed++;
						basesProcessed+=initialLength1;
					}
					
					processReadPair(r1, null);
				}

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		{
			contigToNumber=null;
			ArrayList<Read> list=new ArrayList<Read>(nameToRead.size());
			for(String key : nameToRead.keySet()){
				list.add(nameToRead.get(key));
			}
			nameToRead=null;
			for(int i=0; i<list.size(); i++){
				Read r1=list.set(i, null);
				
				Read r2=r1.mate;
				if(!r1.mapped() && !r1.mateMapped()){
					unmappedReads+=r1.pairCount();
					unmappedBases+=r1.pairLength();
					if(keepUnmapped){
						retainedReads+=r1.pairCount();
						retainedBases+=r1.pairLength();
						unmapped.add(r1);
					}
				}else if(keepSingletons && r2!=null && (r1.mapped()!=r1.mateMapped())){
					retainedReads+=r1.pairCount();
					retainedBases+=r1.pairLength();
					unmapped.add(r1);
				}else{
//					System.err.println(r1.strandChar()+", "+r1.start+", "+r1.stop+", "+r2.strandChar()+", "+r2.start+", "+r2.stop);
					Quad q=toQuad(r1, r2);
					Read old1=quadToRead.get(q);
					if(old1==null){quadToRead.put(q, r1);}
					else{
						Read old2=old1.mate;
						float a=(r1.expectedErrors(true, 0)+(r2==null ? 0 : r2.expectedErrors(true, 0)))/r1.pairLength();
						float b=old1.expectedErrors(true, 0)+(old2==null ? 0 : old2.expectedErrors(true, 0))/(old1.length()+old1.mateLength());
						if(a<b){
							quadToRead.put(q, r1);
							duplicateReads+=1+old1.mateCount();
							duplicateBases+=old1.length()+old1.mateLength();
						}else{
							duplicateReads+=r1.pairCount();
							duplicateBases+=r1.pairLength();
						}
					}
				}
			}
			list=null;
			nameToRead=null;
		}
		
		{
			ArrayList<Read> list=new ArrayList<Read>(Shared.bufferLen());
			int num=0;
			for(Quad q : quadToRead.keySet()){
				Read r=quadToRead.get(q);
				if(keepUnmapped || r.mapped() || (r.mate!=null && r.mate.mapped())){
					retainedReads+=1+r.mateCount();
					retainedBases+=r.length()+r.mateLength();
					list.add(r);
					if(list.size()>=Shared.bufferLen()){
						if(ros!=null){
							ros.add(list, num);
							num++;
						}
						list=new ArrayList<Read>(Shared.bufferLen());
					}
				}
			}
			if(list.size()>0){
				if(ros!=null){
					ros.add(list, num);
					num++;
				}
				list=null;
			}
			if(ros!=null && unmapped.size()>0){
				ros.add(unmapped, num);
				num++;
			}
		}
		outstream.println("Duplicate reads:     "+duplicateReads+" \t("+duplicateBases+" bases)");
		outstream.println("Unconsidered reads:  "+unmappedReads+" \t("+unmappedBases+" bases)");
		outstream.println("Retained reads:      "+retainedReads+" \t("+retainedBases+" bases)");
	}
	
	
	
	void processInner_sorted(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream ros){
		assert(false) : "TODO";
		readsProcessed=0;
		basesProcessed=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					assert(r1.mate==null);
					assert(r1.samline!=null);
					
					final int initialLength1=r1.length();
					
					{
						readsProcessed++;
						basesProcessed+=initialLength1;
					}
					
					processReadPair(r1, null);
				}

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		{
			contigToNumber=null;
			ArrayList<Read> list=new ArrayList<Read>(nameToRead.size());
			for(String key : nameToRead.keySet()){
				list.add(nameToRead.get(key));
			}
			nameToRead=null;
			for(int i=0; i<list.size(); i++){
				Read r1=list.set(i, null);
				
				Read r2=r1.mate;
				if(!r1.mapped() && !r1.mateMapped()){
					unmappedReads+=r1.pairCount();
					unmappedBases+=r1.pairLength();
					if(keepUnmapped){unmapped.add(r1);}
				}else{
					Quad q=toQuad(r1, r2);
					Read old1=quadToRead.get(q);
					if(old1==null){quadToRead.put(q, r1);}
					else{
						Read old2=old1.mate;
						float a=(r1.expectedErrors(true, 0)+(r2==null ? 0 : r2.expectedErrors(true, 0)))/r1.pairLength();
						float b=old1.expectedErrors(true, 0)+(old2==null ? 0 : old2.expectedErrors(true, 0))/(old1.length()+old1.mateLength());
						if(a<b){
							quadToRead.put(q, r1);
							duplicateReads+=1+old1.mateCount();
							duplicateBases+=old1.length()+old1.mateLength();
						}else{
							duplicateReads+=r1.pairCount();
							duplicateBases+=r1.pairLength();
						}
					}
				}
			}
			list=null;
			nameToRead=null;
		}
		
		{
			ArrayList<Read> list=new ArrayList<Read>(Shared.bufferLen());
			int num=0;
			for(Quad q : quadToRead.keySet()){
				Read r=quadToRead.get(q);
				if(keepUnmapped || r.mapped() || (r.mate!=null && r.mate.mapped())){
					list.add(r);
					if(list.size()>=Shared.bufferLen()){
						if(ros!=null){
							ros.add(list, num);
							num++;
						}
						list=new ArrayList<Read>(Shared.bufferLen());
					}
				}
			}
			if(list.size()>0){
				if(ros!=null){
					ros.add(list, num);
					num++;
				}
				list=null;
			}
			if(ros!=null && unmapped.size()>0){
				ros.add(unmapped, num);
				num++;
			}
		}
		outstream.println("Duplicate reads:    "+duplicateReads+" \t("+duplicateBases+" bases)");
		outstream.println("Unmapped reads:     "+unmappedReads+" \t("+unmappedBases+" bases)");
	}
	
	@Override
	protected void startupSubclass() {}
	
	@Override
	protected void shutdownSubclass() {}
	
	@Override
	protected void showStatsSubclass(Timer t, long readsIn, long basesIn) {}
	
	private Quad toQuad(Read r1, Read r2){

//		if(usePairOrder){
//			start1=start1_;
//			start2=start2_;
//			chr1=chr1_;
//			chr2=chr2_;
//		}else{
//			start1=Tools.max(start1_,start2_);
//			start2=Tools.min(start1_,start2_);
//			chr1=Tools.max(chr1_,chr2_);
//			chr2=Tools.min(chr1_,chr2_);
//		}
//
//		int pos1, pos2, chrom1, chrom2;
//
//		if()

		final int s1=r1.strand(), a1=r1.start, b1=r1.stop, c1=r1.chrom;
		final int s2=(r2==null ? 0 : r2.strand()), a2=(r2==null ? 0 : r2.start), b2=(r2==null ? 0 : r2.stop), c2=(r2==null ? 0 : r2.chrom);
		final Quad q;
		if(usePairOrder){
			q=new Quad(s1==0 ? a1 : b1, c1, s2==0 ? a2 : b2, c2);
		}else{
			if(s1==0){
				q=new Quad(s1==0 ? a1 : b1, c1, s2==0 ? a2 : b2, c2);
			}else{
				q=new Quad(s2==0 ? a2 : b2, c2, s1==0 ? a1 : b1, c1);
			}
		}
		
		//q=new Quad((r1.strand()==0 ? r1.start : r1.stop), r1.chrom, r2==null ? -2 : (r2.strand()==0 ? r2.start : r2.stop), r2==null ? -2 : r2.chrom);
		return q;
	}
	
	private static class Quad implements Comparable<Quad>{
		
		Quad(int start1_, int start2_, int chr1_, int chr2_){
			start1=start1_;
			start2=start2_;
			chr1=chr1_;
			chr2=chr2_;
			
//			System.err.println(usePairOrder+", "+this);
		}
		
		@Override
		public String toString(){
			return "("+start1+","+start2+","+chr1+","+chr2+")";
		}
		
		@Override
		public int hashCode(){
			return start1^(Integer.rotateLeft(start2, 8))^(Integer.rotateLeft(chr1, 16))^(Integer.rotateLeft(chr2, 24));
		}
		
		@Override
		public boolean equals(Object o){
			return equals((Quad)o);
		}
		
		public boolean equals(Quad o){
			return start1==o.start1 && start2==o.start2 && chr1==o.chr1 && chr2==o.chr2;
		}
		
		@Override
		public int compareTo(Quad b) {
			int x;
			x=chr1-b.chr1;
			if(x!=0){return x;}
			x=start1-b.start1;
			if(x!=0){return x;}
			x=chr2-b.chr2;
			if(x!=0){return x;}
			x=start2-b.start2;
			return x;
		}
		
		final int start1, start2, chr1, chr2;
	}
	
	private boolean keepUnmapped;
	private boolean keepSingletons;
	private boolean sorted;
	private static boolean usePairOrder;

	private long duplicateReads=0;
	private long duplicateBases=0;
	private long unmappedReads=0;
	private long unmappedBases=0;
	private long retainedReads=0;
	private long retainedBases=0;
	
	private int initialSize=(int)Tools.min(2000000, Tools.max(80000, Shared.memAvailable(1)/4000));
	
	private HashMap<String, Integer> contigToNumber=new HashMap<String, Integer>(initialSize);
	private LinkedHashMap<Quad, Read> quadToRead=new LinkedHashMap<Quad, Read>(initialSize);
	private LinkedHashMap<String, Read> nameToRead=new LinkedHashMap<String, Read>(initialSize);
	private ArrayList<Read> unmapped=new ArrayList<Read>(initialSize/2);
	
	private PriorityQueue<Quad> queue;
	
}
