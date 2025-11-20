package bbdukGemini;

import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import kmer.AbstractKmerTable;
import shared.KillSwitch;
import shared.Parser;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * Coordinates multi-threaded loading of reference sequences into BBDukIndex.
 * @author Brian Bushnell
 * @contributor Gemini
 * @date November 19, 2025
 */
public class BBDukIndexLoader extends BBDukObject{
	
	private static final ArrayList<Read> POISON=new ArrayList<Read>(0);
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public BBDukIndexLoader(BBDukParser bbdp_, Parser parser_, int[] refScafCounts_, boolean replicateAmbiguous_){
		super(bbdp_, parser_);
		bbdp=bbdp_;
		parser=parser_; // Keep reference to create Index later
		ref=bbdp.ref;
		literal=bbdp.literal;
		replicateAmbiguous=replicateAmbiguous_;
		refScafCounts=refScafCounts_;
		scaffoldNames.add(""); 
		scaffoldLengths.add(0);
	}
	
	public void setKeySets(AbstractKmerTable[] keySets_){
		keySets=keySets_;
		WAYS=keySets.length;
	}
	
	public BBDukIndex getIndex(){
		if(index==null){
			index=new BBDukIndex(keySets, scaffoldNames, bbdp, parser);
		}
		return index;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Fills tables with kmers from references. */
	public long spawnLoadThreads(boolean altRefUsed){
		if((ref==null || ref.length<1) && (literal==null || literal.length<1)){return 0;}
		
		if(index==null){getIndex();}
		
		LoadThread[] loaders=new LoadThread[WAYS];
		for(int i=0; i<loaders.length; i++){
			loaders[i]=new LoadThread(i);
			loaders[i].start();
		}
		
		processReferences(loaders, altRefUsed);
		
		boolean success=true;
		long added=0;
		long refKmersT=0, refBasesT=0, refReadsT=0;
		
		for(LoadThread lt : loaders){
			while(lt.getState()!=Thread.State.TERMINATED){
				try{
					lt.join();
				}catch(InterruptedException e){
					e.printStackTrace();
					Thread.currentThread().interrupt();
				}
			}
			added+=lt.addedT;
			refKmersT+=lt.refKmersT;
			refBasesT+=lt.refBasesT;
			refReadsT+=lt.refReadsT;
			success&=lt.success;
		}
		if(!success){KillSwitch.kill("Failed loading ref kmers; aborting.");}
		
		// Correct statistics for number of threads
		BBDukIndexConstants.refKmers=refKmersT/WAYS;
		BBDukIndexConstants.refBases=refBasesT/WAYS;
		BBDukIndexConstants.refReads=refReadsT/WAYS;
		
		return added;
	}
	
	/** Iterates through reference paths and literal sequences, sending data to LoadThreads. */
	private long processReferences(LoadThread[] loaders, boolean altRefUsed){
		
		String[] currentRef=altRefUsed ? bbdp.altref : ref;
		long totalAdded=0;
		int refNum=0;
		
		if(currentRef!=null){
			for(String refname : currentRef){
				FileFormat ff=FileFormat.testInput(refname, FileFormat.FASTA, null, false, true);
				ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1L, false, ff, null, null, null, shared.Shared.USE_MPI, true);
				cris.start();
				ListNum<Read> ln=cris.nextList();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);
				
				while(ln!=null && reads!=null && reads.size()>0){
					ArrayList<Read> reads2=assignScaffoldIDs(reads, refNum);
					
					if(replicateAmbiguous){
						reads2=Tools.replicateAmbiguous(reads2, Tools.min(k, mink));
					}

					for(LoadThread lt : loaders){
						try{
							lt.queue.put(reads2);
						}catch(InterruptedException e){
							throw new RuntimeException(e);
						}
					}
					cris.returnList(ln);
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				cris.returnList(ln);
				ReadWrite.closeStream(cris);
				refNum++;
			}
		}

		if(literal!=null){
			ArrayList<Read> list=assignScaffoldIDs(bbdp.literal, refNum); 
			
			if(replicateAmbiguous){
				list=Tools.replicateAmbiguous(list, Tools.min(k, mink));
			}

			for(LoadThread lt : loaders){
				try{
					lt.queue.put(list);
				}catch(InterruptedException e){
					throw new RuntimeException(e);
				}
			}
		}
		
		for(LoadThread lt : loaders){
			try{
				lt.queue.put(POISON);
			}catch(InterruptedException e){
				throw new RuntimeException(e);
			}
		}
		
		return totalAdded;
	}
	
	/** Assigns unique ID to each scaffold/read in the list. */
	private ArrayList<Read> assignScaffoldIDs(ArrayList<Read> reads, int refNum){
		ArrayList<Read> reads2=new ArrayList<Read>(reads);
		for(Read r1 : reads2){
			final Read r2=r1.mate;
			final Integer id=scaffoldNames.size();
			
			refScafCounts[refNum]++; 
			
			scaffoldNames.add(r1.id==null ? id.toString() : r1.id);
			int len=r1.length();
			r1.obj=id;
			if(r2!=null){
				r2.obj=id;
				len+=r2.length();
			}
			scaffoldLengths.add(len);
		}
		return reads2;
	}
	
	/** Assigns unique ID to each literal sequence in the array. */
	private ArrayList<Read> assignScaffoldIDs(String[] literals, int refNum){
		ArrayList<Read> list=new ArrayList<Read>(literals.length);
		for(int i=0; i<literals.length; i++){
			final Integer id=scaffoldNames.size();
			final Read r=new Read(literals[i].getBytes(), null, id);
			
			refScafCounts[refNum]++;
			
			scaffoldNames.add(id.toString());
			scaffoldLengths.add(r.length());
			r.obj=id;
			list.add(r);
		}
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads kmers into a table partition. */
	private class LoadThread extends Thread{
		
		public LoadThread(final int tnum_){
			tnum=tnum_;
			index.tnum=tnum; 
		}
		
		private ArrayList<Read> fetch(){
			ArrayList<Read> list=null;
			while(list==null){
				try{
					list=queue.take();
				}catch(InterruptedException e){
					Thread.currentThread().interrupt();
				}
			}
			return list;
		}
		
		@Override
		public void run(){
			ArrayList<Read> reads=fetch();
			while(reads!=POISON){
				for(Read r1 : reads){
					assert(r1.pairnum()==0);
					final Read r2=r1.mate;
					
					final int rblen=(r1==null ? 0 : r1.length());
					final int rblen2=r1.mateLength();
					
					addedT+=index.addReadKmers(r1, rblen>20000000 ? k : rblen>5000000 ? 11 : rblen>500000 ? 2 : 0, (Integer)r1.obj, forbidNs, bbdp.useShortKmers, bbdp.minSkip, bbdp.maxSkip);
					
					if(r2!=null){
						addedT+=index.addReadKmers(r2, rblen2>20000000 ? k : rblen2>5000000 ? 11 : rblen2>500000 ? 2 : 0, (Integer)r2.obj, forbidNs, bbdp.useShortKmers, bbdp.minSkip, bbdp.maxSkip);
					}
					
					refReadsT+=r1.pairCount();
					refBasesT+=r1.pairLength();
				}
				reads=fetch();
			}
			
			AbstractKmerTable map=index.getKeySets()[tnum];
			if(map.canRebalance() && map.size()>2L*map.arrayLength()){
				map.rebalance();
			}
			success=true;
		}

		public long addedT=0;
		public long refKmersT=0, refReadsT=0, refBasesT=0;
		public final int tnum;
		public final ArrayBlockingQueue<ArrayList<Read>> queue=new ArrayBlockingQueue<ArrayList<Read>>(32);
		boolean success=false;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private BBDukIndex index;
	private final BBDukParser bbdp;
	private final Parser parser;
	private AbstractKmerTable[] keySets;
	private int WAYS;
	
	public final ArrayList<String> scaffoldNames=new ArrayList<String>();
	public final structures.IntList scaffoldLengths=new structures.IntList();
	
	private final String[] ref;
	private final String[] literal;
	private final boolean replicateAmbiguous;
	private final int[] refScafCounts;
}