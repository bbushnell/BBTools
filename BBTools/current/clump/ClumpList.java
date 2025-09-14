package clump;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.Map.Entry;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
import shared.Shared;
import shared.Tools;
import stream.Read;
import structures.LongM;

/**
 * A list of clumps, meaning a list of lists of reads.
 * Allows adding reads by streaming and generating new clumps as needed.
 * The input reads must be correctly ordered.
 * @author Brian Bushnell
 * @date Nov 9, 2015
 *
 */
public class ClumpList extends ArrayList<Clump> {

	public ClumpList(int k_){
		this(null, k_, false);
	}
	
	//TODO: Add flag to remove duplicates, and remove optical duplicates
	public ClumpList(ArrayList<Read> list, int k_, boolean makeSimpleConsensus_){
		k=k_;
		makeSimpleConsensus=makeSimpleConsensus_;
		if(list!=null){addReads(list);}
//		Shared.sort(this);//Does nothing since reads are already sorted by kmer.
	}
	
	public void addReads(ArrayList<Read> list){
		assert(list.getClass()!=Clump.class) : list.getClass();
		if(list.size()<100000 || Shared.threads()<2){addReadsST(list);}
		else{addReadsMT(list);}
		
//		System.err.println(addedR+", "+addedC+", "+list.size());
	}
	
	public void reorderPaired(){//More efficient but requires unpair and repair, and paired reads.
		ArrayList<Clump> temp=new ArrayList<Clump>(size());
		for(Clump c : this){
			if(!c.added){
				temp.add(c);
				c.added=true;
				addFriends(c, temp, 0);
			}
		}
		assert(temp.size()==size());
		super.clear();
		this.addAll(temp);
	}
	
	private void addFriends(Clump c, ArrayList<Clump> temp, int depth){
		LinkedList<Clump> q=(depth<100 ? new LinkedList<Clump>() : null);
		for(Read r : c){
			Read r2=r.mate;
			if(r2!=null){
				ReadKey rk=(ReadKey)r2.obj;
				Clump c2=rk.clump;
				if(!c2.added){
					temp.add(c2);
					c2.added=true;
					if(q!=null){
						q.addFirst(c2);
//						addFriends(c2, temp, depth+1);//Not as good.
					}
				}
			}
		}
		while(q!=null && !q.isEmpty()){
			addFriends(q.pollLast(), temp, depth+1);
		}
	}
	
	public void reorder(){//Less efficient but works with unpaired reads
		ArrayList<Clump> temp=new ArrayList<Clump>(size());
//		temp.addAll(this);
		LinkedHashMap<LongM, ArrayList<Clump>> map=map();
		LinkedList<Clump> q=new LinkedList<Clump>();
		
		final LongM key=new LongM();
		for(Entry<LongM, ArrayList<Clump>> e : map.entrySet()){
			ArrayList<Clump> list=e.getValue();
			if(list.size()>1){
				Clump c0=list.get(0);
				if(!c0.added){
					temp.add(c0);
					c0.added=true;
				}
				for(int i=1; i<list.size(); i++){
					Clump c=list.get(i);
					if(!c.added){
						q.addFirst(c);
						c.added=true;
					}
				}
				while(!q.isEmpty()){
					Clump c=q.pollLast();
					temp.add(c);
					key.set(c.kmer);
					ArrayList<Clump> list2=map.get(key);
					if(list2.size()>1){
						for(int i=1; i<list2.size(); i++){
							Clump c2=list2.get(i);
							if(!c2.added){
								c2.added=true;
								q.addFirst(c2);
							}
						}
					}
				}
			}
		}
		for(Clump c : this){
			if(!c.added){
				c.added=true;
				temp.add(c);
			}
		}
		assert(temp.size()==size()) : temp.size()+", "+this.size();
		this.clear();
		this.addAll(temp);
	}
	
	public void addReadsMT(ArrayList<Read> list){
		int threads=Tools.mid(2, Shared.threads()/2, 16);
//		assert(false) : threads;
		final ArrayList<ClumpThread> alct=new ArrayList<ClumpThread>(threads);
		int incr=((list.size()+threads-1)/threads);
		if(incr*threads<list.size()-1){incr++;}
		assert(threads*incr>=list.size()) : threads+", "+incr+", "+list.size();
		for(int i=0; i<threads; i++){
			alct.add(new ClumpThread(list, i*incr, (i+1)*incr, k, makeSimpleConsensus));
		}
		
		if(verbose){outstream.println("Starting clump threads.");}
		for(ClumpThread ct : alct){ct.start();}
		
		if(verbose){outstream.println("Waiting for threads.");}
		/* Wait for threads to die */
		for(ClumpThread ct : alct){
			
			/* Wait for a thread to die */
			while(ct.getState()!=Thread.State.TERMINATED){
				try {
					ct.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			addAll(ct.storage);
			cAdded+=ct.cAddedT;
			rAdded+=ct.rAddedT;
		}
		assert(list.size()==rAdded) : list.size()+"!="+rAdded;
//		for(int i=0; i<list.size(); i++){assert(list.get(i).doscarded()) : i;}
	}
	
//	long addedR=0, addedC=0;
//
//	@Override
//	public boolean add(Clump c){
//		assert(c!=null && c.size()>0);
//		synchronized(this){
//			addedR+=c.size();
//			addedC++;
//		}
//		return super.add(c);
//	}
//
//	@Override
//	public boolean addAll(Collection<? extends Clump> cc){
//		for(Clump c : cc){add(c);}
//		return true;
//	}
	
	public void addReadsST(ArrayList<Read> list){
		Clump currentClump=null;
		long currentKmer=-1;
		for(final Read r : list){
			final ReadKey key=(ReadKey)r.obj;
			if(key.kmer!=currentKmer){
				if(currentClump!=null){
					if(makeSimpleConsensus){
						currentClump.consensusRead();
						currentClump.clearCounts();
					}
					add(currentClump);
				}
				currentKmer=key.kmer;
				currentClump=Clump.makeClump(key.kmer);
			}
			currentClump.add(r);
//			assert(!r.doscarded());
//			r.setDoscarded(true);
//			assert(r.doscarded());
		}
		if(currentClump!=null && !currentClump.isEmpty()){
			if(makeSimpleConsensus){
				currentClump.consensusRead();
				currentClump.clearCounts();
			}
			add(currentClump);
		}
//		for(int i=0; i<list.size(); i++){assert(list.get(i).doscarded()) : i;}
	}
	
	public ArrayList<Read> process(final int threads, final int mode, final long[] rvector){
		final ArrayList<ProcessThread> alct=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){alct.add(new ProcessThread(mode));}
		
		if(verbose){outstream.println("Starting condense/correction threads.");}
		for(ProcessThread ct : alct){ct.start();}
		
		if(verbose){outstream.println("Waiting for threads.");}
		long readsThisPass=0;
		long corrections=0;
		long duplicates=0;
		/* Wait for threads to die */
		for(ProcessThread ct : alct){
			
			/* Wait for a thread to die */
			while(ct.getState()!=Thread.State.TERMINATED){
				try {
					ct.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			readsThisPass+=ct.storage.size();
			corrections+=ct.corrections;
			duplicates+=ct.duplicates;
		}
		
		if(verbose){outstream.println("Gathering reads.");}
		ArrayList<Read> list=new ArrayList<Read>((int)readsThisPass);
		for(int i=0; i<threads; i++){
			ProcessThread ct=alct.set(i, null);
			list.addAll(ct.storage);
		}

		rvector[0]+=corrections;
		rvector[1]+=duplicates;
//		assert(false) : duplicates+", "+Arrays.toString(rvector);
		assert(list.size()==readsThisPass);
		return list;
	}
	
	@Override
	public void clear(){
		super.clear();
		ptr.set(0);
		map=null;
	}
	
	public LinkedHashMap<LongM, ArrayList<Clump>> map(){
		LinkedHashMap<LongM, ArrayList<Clump>> temp=map;
		if(temp==null){
			synchronized(this){
				if(map==null){
					temp=makeMap();
				}
				assert(map==null);
				map=temp;
			}
		}
		return map;
	}
	
	private synchronized LinkedHashMap<LongM, ArrayList<Clump>> makeMap(){
		assert(this.map==null);
		LinkedHashMap<LongM, ArrayList<Clump>> map=new LinkedHashMap<LongM, ArrayList<Clump>>(this.size());
		for(Clump c : this){
			ArrayList<Clump> list=new ArrayList<Clump>(4);
			list.add(c);
			map.put(new LongM(c.kmer), list);
		}
		fillMap(map);
		assert(this.map==null);
		return map;
	}
	
	private void fillMap(LinkedHashMap<LongM, ArrayList<Clump>> map){
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=(shift>63 ? -1L : ~((-1L)<<shift));
		final LongM key=new LongM();
		
		for(Clump c : this){
			Read r=c.consensusRead();
			byte[] bases=r.bases;
			long kmer=0, rkmer=0;
			int len=0;
			for(int i=0; i<bases.length; i++){
				byte b=bases[i];
				long x=AminoAcid.baseToNumber[b];
				long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;

				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{len++;}

				if(len>=k){
					final long kmax=Tools.max(kmer, rkmer);
					key.set(kmax);
					if(kmax!=c.kmer){
						//This assertion should be fine, but it fails on nt.
//						assert(kmer!=c.kmer && rkmer!=c.kmer) : kmax+", "+kmer+", "+rkmer+", "+c.kmer+", "+r.length()+", "+(r.length()<10000 ? r.toString() : "");
					}
					ArrayList<Clump> list=map.get(key);
					if(list!=null && !list.contains(c)){list.add(c);}
				}
			}
		}
	}
	
	private class ProcessThread extends Thread{
		
		public ProcessThread(int mode_){
			mode=mode_;
		}
		
		@Override
		public void run(){
			final int size=size();
			for(int i=ptr.getAndIncrement(); i<size; i=ptr.getAndIncrement()){
				Clump c=get(i);
				if(mode==CONDENSE){
					ArrayList<Read> list=c.makeConsensus();
					storage.addAll(list);
				}else if(mode==CORRECT){
					corrections+=c.splitAndErrorCorrect();
					if(UNRCOMP){
						for(Read r : c){
							if(r.swapped()){
								ReadKey key=(ReadKey)r.obj;
								key.flip(r, k);
//								r.reverseComplement();
//								r.setSwapped(false);
//								key.position=r.length()-key.position+k-2;
							}
						}
					}
					storage.addAll(c);
				}else if(mode==DEDUPE){
					duplicates+=c.removeDuplicates();
					storage.addAll(c);
				}else{
					throw new RuntimeException("Unknown mode "+mode);
				}
				c.clear();
				set(i, null);
			}
		}
		
		public long corrections=0;
		public long duplicates=0;
		ArrayList<Read> storage=new ArrayList<Read>();
		private final int mode;
	}
	
	private static class ClumpThread extends Thread {
		
		public ClumpThread(ArrayList<Read> input_, int startIndex_, int stopIndex_, int k, boolean makeSimpleConsensus_){
			input=input_;
			startIndex=startIndex_;
			stopIndex=Tools.min(input.size(), stopIndex_);
//			outstream.println("Processing "+startIndex_+"-"+stopIndex_);
			storage=new ClumpList(k);
			makeSimpleConsensusT=makeSimpleConsensus_;
		}
		
		@Override
		public void run(){
			if(startIndex>=input.size()){return;}
			int i=startIndex;

			long currentKmer=-1;
			if(startIndex>0){
				Read r=input.get(i-1);
				final ReadKey key=(ReadKey)r.obj;
				currentKmer=key.kmer;
			}
			while(i<stopIndex){
				Read r=input.get(i);
				final ReadKey key=(ReadKey)r.obj;
				if(key.kmer!=currentKmer){
					if(currentClump!=null){
						storage.add(currentClump);
						cAddedT++;
						rAddedT+=currentClump.size();
					}
					currentKmer=key.kmer;
					currentClump=Clump.makeClump(key.kmer);
				}
				if(currentClump!=null){
					currentClump.add(r);
//					assert(!r.doscarded()) : +i+", "+startIndex;
//					r.setDoscarded(true);
				}
//				if(i==0){assert(r.doscarded()) : currentClump;}
				i++;
			}
			while(currentClump!=null && i<input.size()){
				Read r=input.get(i);
				final ReadKey key=(ReadKey)r.obj;
				if(key.kmer!=currentKmer){
					storage.add(currentClump);
					cAddedT++;
					rAddedT+=currentClump.size();
					currentClump=null;
					break;
				}else{
					currentClump.add(r);
//					assert(!r.doscarded()) : +i+", "+startIndex;
//					r.setDoscarded(true);
				}
				i++;
			}
			if(currentClump!=null && currentClump.size()>0){
				cAddedT++;
				rAddedT+=currentClump.size();
				storage.add(currentClump);
			}
			if(makeSimpleConsensusT){
				for(Clump c : storage){
					c.consensusRead();
					c.clearCounts();
				}
			}
			currentClump=null;
		}
		
		final boolean makeSimpleConsensusT;
		final ArrayList<Read> input;
		final int startIndex;
		final int stopIndex;
		ClumpList storage;
		Clump currentClump;
		long cAddedT=0, rAddedT=0;
	}
	
	static final int CONDENSE=1, CORRECT=2, DEDUPE=3;
	
	private final boolean makeSimpleConsensus;
	long cAdded=0, rAdded=0;
	
	final int k;
	final AtomicInteger ptr=new AtomicInteger(0);
	private LinkedHashMap<LongM, ArrayList<Clump>> map;

	public static boolean UNRCOMP=true;
	private static boolean verbose=false;
	
	private static final long serialVersionUID = 1L;
	private static final PrintStream outstream=System.err;
	
}
