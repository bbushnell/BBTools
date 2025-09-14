package bin;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;

import prok.CallGenes;
import prok.GeneCaller;
import prok.Orf;
import shared.KillSwitch;
import shared.LineParserS1;
import shared.LineParserS4;
import shared.Shared;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.IntLongHashMap;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.EntropyTracker;

public class SpectraCounter extends BinObject implements Accumulator<SpectraCounter.LoadThread> {
	
	public SpectraCounter(PrintStream outstream_, boolean parseDepth_, 
			boolean parseTID_, IntLongHashMap sizeMap_) {
		outstream=outstream_;
		parseDepth=parseDepth_;
		parseTID=parseTID_;
		sizeMap=sizeMap_;
		if(calcEntropy) {
			if(AdjustEntropy.kLoaded!=4 || AdjustEntropy.wLoaded!=150) {
				AdjustEntropy.load(4, 150);
			}
			assert(AdjustEntropy.kLoaded==4 && AdjustEntropy.wLoaded==150) : 
				AdjustEntropy.kLoaded+", "+calcEntropy;
		}
	}
	
	/** Spawn process threads */
	public void makeSpectra(ArrayList<Contig> contigs, ConcurrentReadInputStream cris, int minlen){
		
		//Do anything necessary prior to processing
//		sizeMap=(parseTax ? new IntLongHashMap(1021) : null);
		
		GeneTools.setMode(call16S, call18S, false, false, false, false);
		//Determine how many threads may be used
		int threads=Tools.mid(1, cris==null ? contigs.size()/4 : 128, Shared.threads());
		if(loadThreadsOverride>0) {threads=loadThreadsOverride;}
		//Fill a list with LoadThreads
		ArrayList<LoadThread> alpt=new ArrayList<LoadThread>(threads);
		for(int i=0; i<threads; i++){
			LoadThread lt=new LoadThread(contigs, cris, minlen, i, threads);
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
			errorState|=(t.success);
			contigsLoaded+=t.contigsLoadedT;
			basesLoaded+=t.basesLoadedT;
			contigsRetained+=t.contigsRetainedT;
			basesRetained+=t.basesRetainedT;
		}
	}

	@Override
	public ReadWriteLock rwlock() {return null;}

	@Override
	public synchronized boolean success() {return errorState;}
	
	class LoadThread extends Thread {
		
		LoadThread(ArrayList<Contig> contigs_, ConcurrentReadInputStream cris_, 
				int minlen_, int tid_, int threads_) {
			contigs=contigs_;
			cris=cris_;
			minlen=minlen_;
			tid=tid_;
			threads=threads_;
		}
		
		@Override
		public void run() {
			et=new EntropyTracker(entropyK, entropyWindow, false);
			if(call16S || call18S) {
				caller=GeneTools.makeGeneCaller();
			}
			synchronized(this) {
				runInner();
			}
		}
		
		private void runInner() {
			if(cris==null) {//Calculate data on existing contigs
				runOnContigs();
			}else {//Load contigs concurrently
				runOnCris();
			}
			success=true;
		}
		
		void runOnContigs() {
			for(int i=tid; i<contigs.size(); i+=threads) {
				Contig c=contigs.get(i);
				processContig(c);
			}
		}
		
		void runOnCris() {
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();
			
			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){
				ArrayList<Contig> localContigs=new ArrayList<Contig>(ln.size());
				for(Read r : ln) {
					Contig c=loadContig(r);
					if(c!=null) {
						processContig(c);
						localContigs.add(c);
					}
				}
				synchronized(contigs) {contigs.addAll(localContigs);}
				
				//Notify the input stream that the list was used
				cris.returnList(ln);
//				if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access
				
				//Fetch a new list
				ln=cris.nextList();
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		Contig processRead(Read r) {
			Contig c=loadContig(r);
			if(c==null) {return null;}
			processContig(c);
			return c;
		}
		
		Contig loadContig(Read r) {
			contigsLoadedT++;
			basesLoadedT+=r.length();
			int tid=-1;
			if(parseTID) {
				tid=parseTaxID(r.name());
				if(tid>0) {
					synchronized(sizeMap) {
						sizeMap.increment(tid, r.length());
					}
				}
			}
			if(r.length()<minlen) {return null;}
			contigsRetainedT++;
			basesRetainedT+=r.length();
			Contig c=new Contig(r.name(), r.bases, (int)r.numericID);
			byte[][] ssu=callSSU(r);
			synchronized(c) {
				c.labelTaxid=tid;
				if(ssu!=null) {
					c.r16S=ssu[0];
					c.r18S=ssu[1];
				}
			}
			return c;
		}
		
		byte[][] callSSU(Read r) {
			if(caller==null || r.length()<900) {return null;}
			assert(call16S || call18S);
			ArrayList<Orf> genes=caller.callGenes(r);
			if(genes==null || genes.isEmpty()) {return null;}
			byte[] r16s, r18s;
			for(Orf orf : genes) {
				if(orf.is16S()) {
					r16s=CallGenes.fetch(orf, r).bases;
					return new byte[][] {r16s, null};
				}else if(orf.is18S()) {
					r18s=CallGenes.fetch(orf, r).bases;
					return new byte[][] {null, r18s};
				}
			}
			return null;
		}
		
		void processContig(Contig c) {
			synchronized(c) {
//				System.err.println("Thread "+tid+" got lock on "+c.name+", "+c.id()+", "+c.size());
				contigsProcessedT++;
				basesProcessedT+=c.size();
				c.loadCounts();
				if(c.numDepths()>1) {c.fillNormDepth();}
				if(calcEntropy) {
					if(!calcEntropyFast) {
						c.entropy=et.averageEntropy(c.bases, false);
					}else {
						//This would need regeneration of nns but uses much less CPU time loading.
						c.entropy=EntropyTracker.calcEntropyFromCounts(c.trimers);
					}
					c.entropy=AdjustEntropy.compensate(c.gc(), c.entropy);
				}
				if(calcStrandedness) {
					c.dimers=new int[16];
					c.strandedness=EntropyTracker.strandedness(c.bases, c.dimers, 2);
				}
				if(parseDepth) {
					boolean b=DataLoader.parseAndSetDepth(c, lps, lpt);
					if(!b) {
						KillSwitch.kill("Could not parse depth from header "+c.name+
								"\nThis program needs a sam file, a cov file, or labeled contigs.");
					}
					assert(b) : "Could not parse depth from "+c.name;
				}
				
				assert(c.tetramers!=null && c.numTetramers>0);
			}
		}
		
		final int tid;
		final int threads;
		final int minlen;
		final ArrayList<Contig> contigs;
		final ConcurrentReadInputStream cris;
		private EntropyTracker et;
		private GeneCaller caller;
//		final int[] counts=(calcEntropy ? new int[1<<(entropyK*2)] : null);
		boolean success=false;
		int contigsProcessedT=0;
		long basesProcessedT=0;
		LineParserS1 lps=new LineParserS1('_');
		LineParserS4 lpt=new LineParserS4(",,=,");
		
		int contigsLoadedT=0;
		long basesLoadedT=0;
		int contigsRetainedT=0;
		long basesRetainedT=0;
	}
	
	public PrintStream outstream=System.err;
	
	public final boolean parseDepth;
	public final boolean parseTID;
	public final IntLongHashMap sizeMap;

	public int contigsLoaded=0;
	public long basesLoaded=0;
	public int contigsRetained=0;
	public long basesRetained=0;
	
	public boolean errorState=false;
	public static boolean calcEntropy=true;
	public static boolean calcEntropyFast=false;
	public static boolean calcStrandedness=true;
	public static boolean call16S=false;
	public static boolean call18S=false;
	public static int loadThreadsOverride=-1;
	
}
