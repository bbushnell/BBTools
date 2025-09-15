package bin;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.locks.ReadWriteLock;

import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;
import stream.SamLine;
import stream.SamLineStreamer;
import structures.IntHashMap;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.EntropyTracker;

public class SamLoader implements Accumulator<SamLoader.LoadThread> {
	
	public SamLoader(PrintStream outstream_) {
		outstream=outstream_;
	}
	
	@Deprecated
	public void load(ArrayList<String> fnames, HashMap<String, Contig> contigMap, IntHashMap[] graph) {
		//Contig list should already be sorted and numbered.
		ArrayList<Contig> list=new ArrayList<Contig>(contigMap.values());
		Collections.sort(list);
		for(int i=0; i<list.size(); i++) {list.get(i).setID(i);}
		load(fnames, contigMap, list, graph);
	}
	
	/** Spawn process threads */
	public void load(ArrayList<String> fnames, HashMap<String, Contig> contigMap, 
			ArrayList<Contig> contigs, IntHashMap[] graph){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=fnames.size();
		
		//Fill a list with LoadThreads
		ArrayList<LoadThread> alpt=new ArrayList<LoadThread>(threads);
		for(int i=0; i<threads; i++){
			String fname=fnames.get(i);
			LoadThread lt=new LoadThread(fname, i, contigMap, contigs, graph);
			alpt.add(lt);
		}
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		
	}
	
	@Override
	public synchronized void accumulate(LoadThread t) {
		synchronized(t) {
			readsIn+=t.readsInT;
			readsUsed+=t.readsUsedT;
			basesIn+=t.basesInT;
			bytesIn+=t.bytesInT;
			errorState|=(t.success);
		}
	}

	@Override
	public ReadWriteLock rwlock() {return null;}

	@Override
	public synchronized boolean success() {
		return errorState;
	}
	
	class LoadThread extends Thread {
		
		LoadThread(final String fname_, final int sample_, HashMap<String, 
				Contig> contigMap_, ArrayList<Contig> contigs_, IntHashMap[] graph_) {
			fname=fname_;
			sample=sample_;
			contigMap=contigMap_;
			contigs=contigs_;
			graph=graph_;
		}
		
		@Override
		public void run() {
			synchronized(this) {runInner();}
		}
		
		private void runInner() {
			long[] depthArray=new long[contigs.size()];
			SamLineStreamer ss=null;
			outstream.println("Loading "+fname);
			final int streamerThreads=Tools.min(3, Shared.threads());
			FileFormat ff=FileFormat.testInput(fname, FileFormat.SAM, null, true, false);
			ss=new SamLineStreamer(ff, streamerThreads, false, -1);
			ss.start();
			processSam_Thread(ss, depthArray);
			
			postprocess(depthArray);
			depthArray=null;
			success=true;
		}
		
		private void postprocess(long[] depthArray) {
			for(int cnum=0; cnum<depthArray.length; cnum++) {
				Contig c=contigs.get(cnum);
				float depth=depthArray[cnum]*1f/Tools.max(1, c.size());
				synchronized(c) {c.setDepth(depth, sample);}
			}
		}
		
		void processSam_Thread(SamLineStreamer ss, long[] depthArray) {
			ListNum<SamLine> ln=ss.nextLines();
			ArrayList<SamLine> reads=(ln==null ? null : ln.list);

			while(ln!=null && reads!=null && reads.size()>0){

				for(int idx=0; idx<reads.size(); idx++){
					SamLine sl=reads.get(idx);
					if(sl.mapped()) {
						boolean used=addSamLine(sl, depthArray);
						readsUsedT+=(used ? 1 : 0);
						readsInT++;
						basesInT+=(sl.seq==null ? 0 : sl.length());
						bytesInT+=(sl.countBytes());
					}
				}
				ln=ss.nextLines();
				reads=(ln==null ? null : ln.list);
			}
		}
		
		private int calcAlignedBases(SamLine sl, int contigLen) {
			int aligned=sl.mappedNonClippedBases();
			if(contigLen<1.5f*tipLimit) {return aligned;}
			int limit=Tools.min(tipLimit, contigLen/4);
			final int lineStart=sl.start(false, false);
			final int lineStop=sl.stop(lineStart, false, false);
			final int contigStart=limit;
			final int contigStop=contigLen-limit;
			if(lineStart>=contigStart && lineStop<=contigStop) {return aligned;}
			return Tools.overlapLength(lineStart, lineStop, contigStart, contigStop);
		}
		
		private boolean addSamLine(SamLine sl, long[] depthArray) {
			if(!sl.mapped()) {return false;}
			if(maxSubs<999 && sl.countSubs()>maxSubs) {return false;}
			if(minID>0 && sl.calcIdentity()<minID) {return false;}
			final String rname=ContigRenamer.toShortName(sl.rname());
			final Contig c1=contigMap.get(rname);
			if(c1==null) {return false;}//Contig not found; possibly too short
			assert(c1!=null) : "Can't find contig for rname "+rname;
			final int cid=c1.id();
			final int aligned=calcAlignedBases(sl, (int)c1.size());
			depthArray[cid]+=aligned;
			
			if(graph==null || sl.ambiguous() || !sl.hasMate() || !sl.nextMapped() 
					|| sl.pairedOnSameChrom() || sl.mapq<minMapq || aligned<minAlignedBases) {return true;}
			if(minMateq>0) {
				int mateq=sl.mateq();
				if(mateq>=0 && mateq<minMateq) {
//					System.err.println("mateq too low: "+mateq);
					return true;
				}
			}
			if(minMateID>0){
				float mateid=sl.mateID();
				if(mateid>0 && mateid<100*minMateID) {
//					System.err.println("mateid too low: "+mateid);
					return true;
				}
			}
			final String rnext=ContigRenamer.toShortName(sl.rnext());
			assert(rnext!=null && !"*".equals(rnext) && !"=".equals(rnext));
			
			final Contig c2=contigMap.get(rnext);
			if(c2==null) {
				System.err.println("Can't find "+rnext);
				return true;
			}//Contig not found
//			System.err.println("Adding edge: "+rname+" - "+rnext);
			if(minEntropy>0 && sl.seq!=null && !et.passes(sl.seq, true)) {return true;}
//			if(minEntropy>0 && sl.seq!=null && EntropyTracker.calcEntropy(sl.seq, kmerCounts, 4)<minEntropy) {return true;}
			assert(c2!=null) : "Can't find contig for rnext "+rnext;
			
			IntHashMap destMap=graph[cid];
			if(destMap==null) {
				synchronized(graph) {
					if(graph[cid]==null) {graph[cid]=new IntHashMap(5);}
					destMap=graph[cid];
				}
			}
			synchronized(destMap) {
				destMap.increment(c2.id());
			}
			return true;
		}
		
		final String fname;
		final int sample;
		final HashMap<String, Contig> contigMap;
		final ArrayList<Contig> contigs;
		final IntHashMap[] graph;
		final EntropyTracker et=new EntropyTracker(5, 80, false, minEntropy, true);
		final int[] kmerCounts=new int[256];
		long readsInT=0;
		long readsUsedT=0;
		long basesInT=0;
		long bytesInT=0;
		boolean success=false;
	}
	
	public PrintStream outstream=System.err;
	public long readsIn=0;
	public long readsUsed=0;
	public long basesIn=0;
	public long bytesIn=0;
	public int minMapq=4;
	public int minMateq=4;
	public float minID=0f;
	public float minMateID=0f;
	public int maxSubs=999;
	public int tipLimit=100;
	public float minEntropy=0;
	public int minAlignedBases=0;
	
	public boolean errorState=false;
	
}
