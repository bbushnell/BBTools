package assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
import jgi.BBMerge;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;
import structures.LongList;
import ukmer.AbstractKmerTableU;
import ukmer.HashArrayU1D;
import ukmer.HashForestU;
import ukmer.Kmer;
import ukmer.KmerNodeU;
import ukmer.KmerTableSetU;


/**
 * Long-kmer assembler based on KmerCountExact.
 * @author Brian Bushnell
 * @date May 15, 2015
 *
 */
public class Tadpole2 extends Tadpole {
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer(), t2=new Timer();
		t.start();
		t2.start();
		
		//Create a new CountKmersExact instance
		Tadpole2 wog=new Tadpole2(args, true);
		t2.stop();
		outstream.println("Initialization Time:      \t"+t2);
		
		///And run it
		wog.process(t);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Tadpole2(String[] args, boolean setDefaults){
		super(args, setDefaults);
		
		final int extraBytesPerKmer;
		{
			int x=0;
			if(useOwnership){x+=4;}
			if(processingMode==correctMode || processingMode==discardMode){}
			else if(processingMode==contigMode || processingMode==extendMode){x+=1;}
			extraBytesPerKmer=x;
		}
		
		tables=new KmerTableSetU(args, extraBytesPerKmer);
		assert(kbig==tables.kbig) : kbig+", "+tables.kbig;
//		kbig=tables.kbig;
		ksmall=tables.k;
//		k2=tables.k2;
//		ways=tables.ways;
	}

	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	void initializeOwnership(){
		tables.initializeOwnership();
	}
	
	@Override
	long shave(boolean shave, boolean rinse){
		long sum=0;

		for(int i=0; i<maxShaveDepth; i++){
			int a=1, b=maxShaveDepth, c=i+1;
			//				if(i>3){Shaver2.verbose2=true;}
			outstream.println("\nShave("+a+", "+b+", "+c+")");
			final Shaver shaver=Shaver.makeShaver(tables, THREADS, a, b, c, minCountExtend, branchMult2, Tools.max(minContigLen, shaveDiscardLen), shaveExploreDist, shave, rinse);
			long removed=shaver.shave(a, b);
			
			sum+=removed;
			if(removed<100 || i>2){break;}
		}

		outstream.println();
		return sum;
	}
	
	@Override
	public long loadKmers(Timer t){
		tables.process(t);
		return tables.kmersLoaded;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Recall Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Gets the count/depth for the specified k-mer */
	public final int getCount(Kmer kmer){return tables.getCount(kmer);}
	@Override
	final int bridgeCount(Kmer kmer){return getCount(kmer);}
	@Override
	final int bridgeRightCounts(Kmer kmer, int[] counts){return fillRightCounts(kmer, counts);}
	/** Attempts to claim ownership of a k-mer for the specified thread ID */
	final boolean claim(Kmer kmer, int id){return tables.claim(kmer, id);}
	/**
	 * Attempts to claim ownership of all k-mers in the sequence for thread safety
	 */
	final boolean doubleClaim(ByteBuilder bb, int id/*, long rid*/, Kmer kmer){return tables.doubleClaim(bb, id/*, rid*/, kmer);}
	/** Attempts to claim ownership of k-mers in sequence with early exit option */
	final boolean claim(ByteBuilder bb, int id, /*long rid, */boolean earlyExit, Kmer kmer){return tables.claim(bb, id/*, rid*/, earlyExit, kmer);}
	/**
	 * Attempts to claim ownership of k-mers in byte array with early exit option
	 */
	final boolean claim(byte[] array, int len, int id, /*long rid, */boolean earlyExit, Kmer kmer){return tables.claim(array, len, id/*, rid*/, earlyExit, kmer);}
	/** Returns the thread ID that owns the specified k-mer, or -1 if unowned */
	final int findOwner(Kmer kmer){return tables.findOwner(kmer);}
	/** Finds owner of k-mers in sequence using working k-mer for efficiency */
	final int findOwner(ByteBuilder bb, int id, Kmer kmer){return tables.findOwner(bb, id, kmer);}
	/** Finds owner of k-mers in byte array using working k-mer for efficiency */
	final int findOwner(byte[] array, int len, int id, Kmer kmer){return tables.findOwner(array, len, id, kmer);}
	/** Releases ownership of a k-mer from the specified thread ID */
	final void release(Kmer kmer, int id){tables.release(kmer, id);}
	/** Releases ownership of all k-mers in the sequence from thread */
	final void release(ByteBuilder bb, int id, Kmer kmer){tables.release(bb, id, kmer);}
	/** Releases ownership of all k-mers in byte array from thread */
	final void release(byte[] array, int len, int id, Kmer kmer){tables.release(array, len, id, kmer);}
	/** Fills array with counts for right-extension bases from the k-mer */
	final int fillRightCounts(Kmer kmer, int[] counts){return tables.fillRightCounts(kmer, counts);}
	/** Fills array with counts for left-extension bases from the k-mer */
	final int fillLeftCounts(Kmer kmer, int[] counts){return tables.fillLeftCounts(kmer, counts);}
	/** Converts k-mer to human-readable DNA sequence string */
	final static StringBuilder toText(Kmer kmer){return AbstractKmerTableU.toText(kmer);}
	/** Converts k-mer key array to human-readable DNA sequence string */
	final static StringBuilder toText(long[] key, int k){return AbstractKmerTableU.toText(key, k);}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------          BuildThread         ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	BuildThread makeBuildThread(int id, int mode, ConcurrentReadInputStream[] crisa){
		return new BuildThread(id, mode, crisa);
	}
	
	/**
	 * Builds contigs.
	 */
	private class BuildThread extends AbstractBuildThread{
		
		/**
		 * Constructs BuildThread with specified parameters.
		 * @param id_ Thread identifier
		 * @param mode_ Processing mode
		 * @param crisa_ Array of input streams for read processing
		 */
		public BuildThread(int id_, int mode_, ConcurrentReadInputStream[] crisa_){
			super(id_, mode_, crisa_);
		}
		
		@Override
		public void run(){
			if(crisa==null || crisa.length==0){
				//Build from kmers
				
				if(id==0){outstream.print("Seeding with min count = ");}
				String comma="";
				for(int i=contigPasses-1; i>0; i--){
					minCountSeedCurrent=(int)Tools.min(Integer.MAX_VALUE, Tools.max(minCountSeed+i, (long)Math.floor((minCountSeed)*Math.pow(contigPassMult, i)*0.92-0.25) ));
					if(id==0){
						outstream.print(comma+minCountSeedCurrent);
						comma=", ";
					}
					while(processNextTable(nextTable[i])){}
					while(processNextVictims(nextVictims[i])){}
				}
				//Final pass
				minCountSeedCurrent=minCountSeed;
				if(id==0){outstream.println(comma+minCountSeedCurrent);}
				while(processNextTable(nextTable[0])){}
				while(processNextVictims(nextVictims[0])){}
			}else{
				//Extend reads
				for(ConcurrentReadInputStream cris : crisa){
					synchronized(crisa){
						if(!cris.started()){
							cris.start();
						}
					}
					run(cris);
				}
			}
		}
		
		/**
		 * Processes the next available k-mer table for contig building.
		 * @param aint Atomic counter for table assignment
		 * @return true if a table was processed, false if no more tables
		 */
		private boolean processNextTable(AtomicInteger aint){
			final int tnum=aint.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
			final HashArrayU1D table=tables.getTable(tnum);
			final int max=table.arrayLength();
			if(verbose && id==0){outstream.println("Processing table "+tnum+", size "+table.size()+", length "+max);}
			for(int cell=0; cell<max; cell++){
				if(verbose && id==0){outstream.println("Processing cell "+cell);}
				int x=processCell(table, cell, myKmer);
			}
			return true;
		}
		
		/**
		 * Processes overflow k-mers from forest data structure.
		 * @param aint Atomic counter for forest assignment
		 * @return true if a forest was processed, false if no more forests
		 */
		private boolean processNextVictims(AtomicInteger aint){
			final int tnum=aint.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
			final HashArrayU1D table=tables.getTable(tnum);
			final HashForestU forest=table.victims();
			if(verbose && id==0){outstream.println("Processing forest "+tnum+", size "+forest.size());}
			final int max=forest.arrayLength();
			for(int cell=0; cell<max; cell++){
				KmerNodeU kn=forest.getNode(cell);
				int x=traverseKmerNodeU(kn);
			}
			return true;
		}
		
		/**
		 * Processes a single hash table cell for contig building.
		 * Checks count threshold and ownership before attempting assembly.
		 *
		 * @param table Hash table containing the cell
		 * @param cell Cell index to process
		 * @param kmer Working k-mer object
		 * @return Length of contig created, or 0 if none
		 */
		private int processCell(HashArrayU1D table, int cell, Kmer kmer){
			int count=table.readCellValue(cell);
			if(count<minCountSeedCurrent){
				if(verbose){outstream.println("For cell "+cell+", count="+count);}
				return 0;
			}
			
			kmer=table.fillKmer(cell, kmer);
//			assert(kmer.verify(false));
//			assert(kmer.verify(true));

			if(verbose){outstream.println("id="+id+" processing cell "+cell+"; \tkmer="+kmer);}
			if(useOwnership){
				int owner=table.getCellOwner(cell);
				if(verbose){outstream.println("Owner is initially "+owner);}
				if(owner>-1){return 0;}
				owner=table.setOwner(kmer, id, cell);
				if(verbose){outstream.println("Owner is now "+owner);}
				if(owner!=id){return 0;}
			}
			return processKmer(kmer);
		}
		
		/**
		 * Recursively traverses k-mer tree nodes for processing.
		 * @param kn Root node of subtree to traverse
		 * @return Total length of contigs created from this subtree
		 */
		private int traverseKmerNodeU(KmerNodeU kn){
			int sum=0;
			if(kn!=null){
				sum+=processKmerNodeU(kn);
				if(kn.left()!=null){
					sum+=traverseKmerNodeU(kn.left());
				}
				if(kn.right()!=null){
					sum+=traverseKmerNodeU(kn.right());
				}
			}
			return sum;
		}
		
		/**
		 * Processes a single k-mer node for contig building.
		 * Checks count and ownership before attempting to build contig.
		 * @param kn K-mer node to process
		 * @return Length of contig created, or 0 if none
		 */
		private int processKmerNodeU(KmerNodeU kn){
			final long[] key=kn.pivot();
			final int count=kn.getValue(key);
			if(count<minCountSeedCurrent){return 0;}

			if(verbose){outstream.println("id="+id+" processing KmerNodeU; \tkmer="+Arrays.toString(key)+"\t"+toText(key, kbig));}
			if(useOwnership){
				int owner=kn.getOwner(key);
				if(verbose){outstream.println("Owner is initially "+owner);}
				if(owner>-1){return 0;}
				owner=kn.setOwner(key, id);
				if(verbose){outstream.println("Owner is now "+owner);}
				if(owner!=id){return 0;}
			}
			
			myKmer.setFrom(key);
			return processKmer(myKmer);
		}
		
		/**
		 * Creates a contig from a seed k-mer with coverage filtering.
		 * @param kmer Seed k-mer for contig building
		 * @return Length of contig created, or 0 if filtered out
		 */
		private int processKmer(Kmer kmer){
			Contig contig=makeContig(builderT, kmer, true);
			if(contig!=null){
				float coverage=tables.calcCoverage(contig, kmer);
				if(coverage<minCoverage || coverage>maxCoverage){return 0;}
				if(verbose){outstream.println("Added "+contig.length());}
				contig.id=(int)contigNum.incrementAndGet();
				contigs.add(contig);
				return contig.length();
			}else{
				if(verbose){outstream.println("Created null contig.");}
			}
			return 0;
		}
		
		/** Processes reads from input stream for extension.
		 * @param cris Concurrent read input stream */
		private void run(ConcurrentReadInputStream cris){
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			//While there are more reads lists...
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				//For each read (or pair) in the list...
				for(int i=0; i<reads.size(); i++){
					final Read r1=reads.get(i);
					final Read r2=r1.mate;
					
					processReadPair(r1, r2);
				}
				
				//Fetch a new read list
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			cris.returnList(ln);
		}

		//TODO: This appears to do read extension but is very confusing.
		/**
		 * Processes a read pair for extension or insert size estimation.
		 * Handles insert mode, error correction, and contig creation from reads.
		 * @param r1 First read of the pair
		 * @param r2 Second read of the pair (may be null)
		 */
		private void processReadPair(Read r1, Read r2){
			if(verbose){outstream.println("Considering read "+r1.id+" "+new String(r1.bases));}
			
			readsInT++;
			basesInT+=r1.length();
			if(r2!=null){
				readsInT++;
				basesInT+=r2.length();
			}
			
			if(mode==insertMode){
				int x=BBMerge.findOverlapStrict(r1, r2, false);
				if(x<1){
					x=findInsertSize(r1, r2, rightCounts, myKmer, myKmer2);
				}
				insertSizes.increment(Tools.max(x, 0));
				return;
			}
			
			if(ecco && r1!=null && r2!=null && !r1.discarded() && !r2.discarded()){BBMerge.findOverlapStrict(r1, r2, true);}
			if(r1!=null){
				if(r1.discarded()){
					lowqBasesT+=r1.length();
					lowqReadsT++;
				}else{
					byte[] bases=makeContig(r1.bases, builderT, r1.numericID, myKmer);
					if(bases!=null){
						if(verbose){outstream.println("Added "+bases.length);}
						final long num=contigNum.incrementAndGet();
						Contig temp=new Contig(bases, "contig_"+num+"_length_"+bases.length, (int)num);
						contigs.add(temp);
					}
				}
			}
			if(r2!=null){
				if(r2.discarded()){
					lowqBasesT+=r2.length();
					lowqReadsT++;
				}else{
					byte[] bases=makeContig(r2.bases, builderT, r2.numericID, myKmer);
					if(bases!=null){
						if(verbose){outstream.println("Added "+bases.length);}
						final long num=contigNum.incrementAndGet();
						Contig temp=new Contig(bases, "contig_"+num+"_length_"+bases.length, (int)num);
						contigs.add(temp);
					}
				}
			}
		}
		
		/** From kmers */
		//TODO: Possible bug [assemble/Tadpole2#002] - this makeContig never checks IGNORE_BAD_OWNER on a BAD_OWNER extendToRight result (always release+return null), whereas Tadpole1.makeContig honors IGNORE_BAD_OWNER (trim last base / bb.length-- and keep the contig). So the `ibo=t` / ignorebadowner flag silently no-ops for k>31. Non-default flag => LOW; feature-parity gap vs the k<=31 path.
		private Contig makeContig(final ByteBuilder bb, Kmer kmer, boolean alreadyClaimed){
			bb.setLength(0);
			bb.appendKmer(kmer);
			if(verbose){outstream.println("Filled bb: "+bb);}
			
			final int initialLength=bb.length();
			assert(initialLength==kbig);
			if(initialLength<kbig){return null;}
			
			boolean success=(alreadyClaimed || !useOwnership ? true : claim(kmer, id));
			if(verbose){outstream.println("Thread "+id+" checking owner after setting: "+findOwner(bb, id, kmer));}
			if(!success){
				assert(bb.length()==kbig);
//				release(bb, id); //no need to release
				return null;
			}
			if(verbose  /*|| true*/){outstream.println("Thread "+id+" building contig; initial length "+bb.length());}
			if(verbose){outstream.println("Extending to right.");}
			final int rightStatus, leftStatus;
			float leftRatio=0, rightRatio=0;
			{
				final int status=extendToRight(bb, leftCounts, rightCounts, id, kmer);
				
				if(status==DEAD_END){
					//do nothing
				}else if(status==LOOP){//TODO
					//special case - handle specially, for a loop with no obvious junction, e.g. long tandem repeat.
					//Perhaps, the last kmer should be reclassified as a junction and removed.
				}else if(status==BAD_SEED){
					assert(bb.length()==kbig);
					release(kmer, id);
					return null;
				}else{
					if(bb.length()==kbig){
						if(status==BAD_OWNER){
							release(kmer, id);
							return null;
						}else if(isBranchCode(status)){
							release(kmer, id);
							return null;
						}else{
							throw new RuntimeException("Bad return value: "+status);
						}
					}else{
						if(status==BAD_OWNER){
							release(bb, id, kmer);
							return null;
						}else if(status==F_BRANCH || status==D_BRANCH){
							rightRatio=calcRatio(rightCounts);
						}else if(status==B_BRANCH){
							rightRatio=calcRatio(leftCounts);
						}else{
							throw new RuntimeException("Bad return value: "+status);
						}
					}
				}
				rightStatus=status;
			}
			
//			success=extendToRight(bb, leftCounts, rightCounts, id, kmer);
//			if(!success){
//				release(bb, id, kmer);
//				return null;
//			}
			bb.reverseComplementInPlace();
			if(verbose  /*|| true*/){outstream.println("Extending rcomp to right; current length "+bb.length());}
			
			{
				final int status=extendToRight(bb, leftCounts, rightCounts, id, kmer);
				
				if(status==DEAD_END){
					//do nothing
				}else if(status==LOOP){//TODO
					//special case - handle specially, for a loop with no obvious junction, e.g. long tandem repeat.
					//Perhaps, the last kmer should be reclassified as a junction and removed.
				}else if(status==BAD_SEED){
					assert(false) : bb;//This should never happen.
					assert(bb.length()==kbig);
					release(kmer, id);
					return null;
				}else{
					if(status==BAD_OWNER){
						release(bb, id, kmer);
						return null;
					}else if(status==F_BRANCH || status==D_BRANCH){
						leftRatio=calcRatio(rightCounts);
					}else if(status==B_BRANCH){
						leftRatio=calcRatio(leftCounts);
					}else{
						throw new RuntimeException("Bad return value: "+status);
					}
				}
				leftStatus=status;
			}
//			success=extendToRight(bb, leftCounts, rightCounts, id, kmer);
//			if(!success){
//				release(bb, id, kmer);
//				return null;
//			}
			if(verbose  /*|| true*/){outstream.println("Final length for thread "+id+": "+bb.length());}
			//				if(useOwnership && THREADS==1){assert(claim(bases, bases.length, id, rid));}
			success=(useOwnership ? doubleClaim(bb, id, kmer) : true);
			if(verbose  /*|| true*/){outstream.println("Success for thread "+id+": "+success);}
			
			if(trimEnds>0){bb.trimByAmount(trimEnds, trimEnds);}
			else if(trimCircular && leftStatus==LOOP && rightStatus==LOOP){bb.trimByAmount(0, kbig-1);}
			if(bb.length()>=initialLength+minExtension && (bb.length()>=minContigLen || popBubbles)){
				if(success){
					bb.reverseComplementInPlace();
					byte[] bases=bb.toBytes();
					Contig c=new Contig(bases);
					c.leftCode=leftStatus;
					c.rightCode=rightStatus;
					c.rightRatio=rightRatio;
					c.leftRatio=leftRatio;
					c.tid=taxID;
					if(!c.canonical()){c.rcomp();}
					return c;
				}else{
					//					assert(false) : bb.length()+", "+id;
					release(bb, id, kmer);
					return null;
				}
			}
			if(verbose  /*|| true*/){outstream.println("Contig was too short for "+id+": "+bb.length());}
			return null;
		}
		
		/** From a seed */
		private byte[] makeContig(final byte[] bases, final ByteBuilder bb, long rid, final Kmer kmer){
			if(bases==null || bases.length<kbig){return null;}
//			if(verbose  /*|| true*/){outstream.println("Thread "+id+" checking owner: "+findOwner(bases, bases.length, id));}
			int owner=useOwnership ? findOwner(bases, bases.length, id, kmer) : -1;
			if(owner>=id){return null;}
			boolean success=(useOwnership ? claim(bases, bases.length, id, true, kmer) : true);
			if(verbose  /*|| true*/){outstream.println("Thread "+id+" checking owner after setting: "+findOwner(bases, bases.length, id, kmer));}
			if(!success){
				release(bases, bases.length, id, kmer);
				return null;
			}
			if(verbose  /*|| true*/){outstream.println("Thread "+id+" building contig; initial length "+bases.length);}
			bb.setLength(0);
			bb.append(bases);
			if(verbose){outstream.println("Extending to right.");}
			{
				final int status=extendToRight(bb, leftCounts, rightCounts, id, kmer);
				
				if(status==DEAD_END){
					//do nothing
				}else if(status==LOOP){//TODO
					//special case - handle specially, for a loop with no obvious junction, e.g. long tandem repeat.
					//Perhaps, the last kmer should be reclassified as a junction and removed.
				}else if(status==BAD_SEED){
					//do nothing
				}else{
					if(status==BAD_OWNER){
						release(bb.array, bb.length(), id, kmer);
						return null;
					}else if(isBranchCode(status)){
						//do nothing
					}else{
						throw new RuntimeException("Bad return value: "+status);
					}
				}
			}
//			success=extendToRight(bb, leftCounts, rightCounts, id, kmer);
//			if(!success){
//				release(bb.array, bb.length(), id, kmer);
//				return null;
//			}
			bb.reverseComplementInPlace();
			if(verbose  /*|| true*/){outstream.println("Extending rcomp to right; current length "+bb.length());}
			{
				final int status=extendToRight(bb, leftCounts, rightCounts, id, kmer);
				
				if(status==DEAD_END){
					//do nothing
				}else if(status==LOOP){//TODO
					//special case - handle specially, for a loop with no obvious junction, e.g. long tandem repeat.
					//Perhaps, the last kmer should be reclassified as a junction and removed.
				}else if(status==BAD_SEED){
					//do nothing
				}else{
					if(status==BAD_OWNER){
						release(bb.array, bb.length(), id, kmer);
						return null;
					}else if(isBranchCode(status)){
						//do nothing
					}else{
						throw new RuntimeException("Bad return value: "+status);
					}
				}
			}
//			success=extendToRight(bb, leftCounts, rightCounts, id, kmer);
//			if(!success){
//				release(bb.array, bb.length(), id, kmer);
//				return null;
//			}
			if(verbose  /*|| true*/){outstream.println("Final length for thread "+id+": "+bb.length());}
			//				if(useOwnership && THREADS==1){assert(claim(bases, bases.length, id, rid));}
			success=(useOwnership ? doubleClaim(bb, id, kmer) : true);
			if(verbose  /*|| true*/){outstream.println("Success for thread "+id+": "+success);}
			if(bb.length()>=bases.length+minExtension && (bb.length()>=minContigLen || popBubbles)){
				if(success){
					bb.reverseComplementInPlace();
					return bb.toBytes();
				}else{
					//					assert(false) : bb.length()+", "+id;
					release(bb.array, bb.length(), id, kmer);
					return null;
				}
			}
			if(verbose  /*|| true*/){outstream.println("Contig was too short for "+id+": "+bb.length());}
			return null;
		}
		
		/*--------------------------------------------------------------*/
		
		/** Thread-local k-mer object for assembly operations */
		private final Kmer myKmer=new Kmer(kbig);
		/** Second thread-local k-mer object for paired operations */
		private final Kmer myKmer2=new Kmer(kbig);
		
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------       Contig Processing      ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	ProcessContigThread makeProcessContigThread(ArrayList<Contig> contigs, AtomicInteger next){
		return new ProcessContigThread(contigs, next);
	}
	
	@Override
	public void initializeContigs(ArrayList<Contig> contigs){
		tables.clearOwnership();
		tables.initializeOwnership();
		final Kmer kmer=new Kmer(kbig);
		if(crossKGraph){
			initializeCrossKContigs(contigs, kmer);
			return;
		}
		final int invalidOwner=contigs.size();
		for(int i=0; i<contigs.size(); i++){contigs.get(i).id=i;}
		for(Contig c : contigs){
			c.leftKmer(kmer);
			claimGraphEnd(kmer, c.id, invalidOwner);
			c.rightKmer(kmer);
			claimGraphEnd(kmer, c.id, invalidOwner);
		}
	}

	/** Claims one full-graph endpoint, invalidating kmers shared by different contigs. */
	private void claimGraphEnd(final Kmer kmer, final int owner, final int invalidOwner){
		final int old=tables.findOwner(kmer);
		if(old<0){tables.claim(kmer, owner);}
		else if(old!=owner){tables.claim(kmer, invalidOwner);}
	}

	@Override
	void fillOwnerPath(final byte[] bases, final IntList path, final LongList seen,
			final ReadThreadedXResolver resolver){
		path.clear();
		final Kmer kmer=getLocalKmer();
		for(byte b : bases){
			kmer.addRight(b);
			if(kmer.len<1){
				if(path.size>2){resolver.observePath(path, seen);}
				path.clear();
			}else if(kmer.len>=kbig){
				final int owner=tables.findOwner(kmer);
				if(owner>=0 && owner<resolver.contigCount()
						&& (path.size==0 || path.array[path.size-1]!=owner)){path.add(owner);}
			}
		}
		if(path.size>2){resolver.observePath(path, seen);}
		path.clear();
	}

	@Override
	boolean crossKGraph(){return crossKGraph;}

	/** Marks assembled short kmers invalid, then claims only unique eligible contig tips. */
	private void initializeCrossKContigs(ArrayList<Contig> contigs, Kmer kmer){
		final int invalidOwner=contigs.size();
		int eligible=0;
		for(int i=0; i<contigs.size(); i++){contigs.get(i).id=i;}
		for(Contig c : contigs){invalidateCrossKInternal(c.bases, invalidOwner, kmer);}
		for(Contig c : contigs){
			if(c.leftBridgeEndpoint){eligible++;}
			if(c.rightBridgeEndpoint){eligible++;}
		}
		int internalBlocked=0;
		for(Contig c : contigs){
			c.leftKmer(kmer);
			if(c.leftBridgeEndpoint && tables.findOwner(kmer)>=0){internalBlocked++;}
			c.rightKmer(kmer);
			if(c.rightBridgeEndpoint && tables.findOwner(kmer)>=0){internalBlocked++;}
		}
		for(Contig c : contigs){
			if(!c.leftBridgeEndpoint){c.leftKmer(kmer); claimCrossKEnd(kmer, false, c.id, invalidOwner);}
			if(!c.rightBridgeEndpoint){c.rightKmer(kmer); claimCrossKEnd(kmer, false, c.id, invalidOwner);}
		}
		int tipBlocked=0;
		for(Contig c : contigs){
			c.leftKmer(kmer);
			if(c.leftBridgeEndpoint && tables.findOwner(kmer)>=0){tipBlocked++;}
			c.rightKmer(kmer);
			if(c.rightBridgeEndpoint && tables.findOwner(kmer)>=0){tipBlocked++;}
		}
		tipBlocked-=internalBlocked;
		for(Contig c : contigs){
			if(c.leftBridgeEndpoint){c.leftKmer(kmer); claimCrossKEnd(kmer, true, c.id, invalidOwner);}
			if(c.rightBridgeEndpoint){c.rightKmer(kmer); claimCrossKEnd(kmer, true, c.id, invalidOwner);}
		}
		int unique=0;
		for(Contig c : contigs){
			c.leftKmer(kmer);
			c.leftBridgeEndpoint&=(tables.findOwner(kmer)==c.id);
			c.rightKmer(kmer);
			c.rightBridgeEndpoint&=(tables.findOwner(kmer)==c.id);
			if(c.leftBridgeEndpoint){unique++;}
			if(c.rightBridgeEndpoint){unique++;}
		}
		if(verbose){
			outstream.println("Cross-k endpoints: eligible="+eligible+", unique="+unique+
					", internal="+internalBlocked+", ineligibleTip="+tipBlocked+
					", duplicateEligible="+(eligible-unique-internalBlocked-tipBlocked)+".");
		}
	}

	/** Marks every short kmer except the two terminal positions as assembled and untraversable. */
	private void invalidateCrossKInternal(final byte[] bases, final int invalidOwner, final Kmer kmer){
		if(bases.length<=kbig+1){return;}
		kmer.clearFast();
		for(int i=0; i<bases.length; i++){
			kmer.addRight(bases[i]);
			final int start=i-kbig+1;
			if(kmer.len>=kbig && start>0 && start<bases.length-kbig){tables.claim(kmer, invalidOwner);}
		}
	}

	/** Claims one endpoint, invalidating it if ineligible or already seen anywhere. */
	private void claimCrossKEnd(final Kmer kmer, final boolean eligible, final int owner, final int invalidOwner){
		final int old=tables.findOwner(kmer);
		tables.claim(kmer, eligible && old<0 ? owner : invalidOwner);
	}

	/**
	 * Worker thread for processing contigs to find connecting edges.
	 * Explores potential connections between contigs by following k-mer paths
	 * and identifying junction points and branch structures.
	 */
	class ProcessContigThread extends AbstractProcessContigThread {

		/**
		 * Constructs ProcessContigThread with contig list and work counter.
		 * @param contigs_ List of contigs to process
		 * @param next_ Atomic counter for work distribution
		 */
		ProcessContigThread(ArrayList<Contig> contigs_, AtomicInteger next_){
			super(contigs_, next_);
			kmerA=new Kmer(kbig);
			kmerB=new Kmer(kbig);
			kmerC=new Kmer(kbig);
			visited=KillSwitch.allocLong1D(Math.multiplyExact(crossKMaxLen, kmerA.array1().length));
			lastExitCondition=BAD_SEED;
		}

		@Override
		public void processContigLeft(Contig c, int[] leftCounts, int[] rightCounts, int[] extraCounts, ByteBuilder bb){
//			System.err.println("processContigLeft: "+c.id+", "+c.leftCode+", "+c.leftEdges);
			if(crossKGraph ? !c.leftBridgeEndpoint : (!refreshGraphEndpoints && c.leftCode==DEAD_END)){return;}
			
			final Kmer kmer0=c.leftKmer(kmerA);
			final Kmer kmer=kmerB;
			final int owner=tables.findOwner(kmer0);
			if(owner<0){
				endpointSeedsMissingT++;
				if(refreshGraphEndpoints){
					if(c.leftCode!=DEAD_END){graphEndCodesChangedT++;}
					c.leftCode=DEAD_END;
					c.leftRatio=0;
					graphEndsRefreshedT++;
					graphEndCodeCountsT[c.leftCode]++;
					missingEndCodeCountsT[c.leftCode]++;
				}
				return;
			}
			int leftMaxPos=fillLeftCounts(kmer0, leftCounts);
			int leftMax=leftCounts[leftMaxPos];
			int leftSecondPos=Tools.secondHighestPosition(leftCounts);
			int leftSecond=leftCounts[leftSecondPos];
			if(refreshGraphEndpoints){
				fillRightCounts(kmer0, rightCounts);
				final int old=c.leftCode;
				c.leftCode=classifyGraphEnd(rightCounts, leftCounts);
				c.leftRatio=graphEndRatio(c.leftCode, rightCounts, leftCounts);
				graphEndsRefreshedT++;
				if(c.leftCode!=old){graphEndCodesChangedT++;}
				graphEndCodeCountsT[c.leftCode]++;
			}
			//[assemble/Tadpole2#004 FIXED 2026-08-28] Post-bridge contigs can share a terminal
			//graph-k kmer.  Such an endpoint is topologically ambiguous, so skip it instead of
			//assigning the shared kmer to whichever contig happened to claim it last.
			if(owner!=c.id){
				endpointSeedsAmbiguousT++;
				if(refreshGraphEndpoints){ambiguousEndCodeCountsT[c.leftCode]++;}
				return;
			}
			endpointSeedsUniqueT++;
			if(refreshGraphEndpoints){uniqueEndCodeCountsT[c.leftCode]++;}
			
			assert(tables.getCount(kmer0)>0);
			assert(tables.findOwner(kmer0)==c.id) : tables.findOwner(kmer0)+", "+c.id;
			
			for(int x=0; x<leftCounts.length; x++){
				bb.clear();
				final int count=leftCounts[x];
				int target=-1;
				if(count>0 && isJunction(leftMax, count)){
					if(crossKGraph || refreshGraphEndpoints){traversalAttemptsT++;}
					kmer.setFrom(kmer0);
					kmer.addLeftNumeric(x);
					//[assemble/Tadpole2#003 FIXED 2026-08-27] A cross-k left extension was passed
					//directly to exploreRight, so its first step returned to the source tip. Reverse-
					//complement the seed to traverse outward, matching Tadpole1. Keep the established
					//dense-graph representation unchanged; changing it requires separate validation.
					if(crossKGraph){kmer.rcomp();}
					assert(tables.getCount(kmer)==count) : count+", "+tables.getCount(kmer);
					bb.append(AminoAcid.numberToBase[crossKGraph ? 3-x : x]);
					target=exploreRight(kmer, extraCounts, rightCounts, bb, c.id);
					if(crossKGraph || refreshGraphEndpoints){exitCountsT[lastExitCondition]++;}
					if(verbose){
						outstream.println(c.id+"L_F: x="+x+", cnt="+count+", dest="+target
								+", "+codeStrings[lastExitCondition]+", len="+lastLength+", orient="+lastOrientation);
					}
				}
				if(target>=0){
					if(crossKGraph){bb.reverseComplementInPlace();}
					Edge se=new Edge(c.id, target, lastLength, lastOrientation, count, bb.toBytes());
//					System.err.println("Adding "+se+"; x="+x+"; bb="+bb);
					c.addLeftEdge(se);
					edgesMadeT++;
				}
			}
		}

		@Override
		public void processContigRight(Contig c, int[] leftCounts, int[] rightCounts, int[] extraCounts, ByteBuilder bb){
//			System.err.println("processContigRight: "+c.id+", "+c.rightCode);
			if(crossKGraph ? !c.rightBridgeEndpoint : (!refreshGraphEndpoints && c.rightCode==DEAD_END)){return;}

			final Kmer kmer0=c.rightKmer(kmerA);
			final Kmer kmer=kmerB;
			final int owner=tables.findOwner(kmer0);
			if(owner<0){
				endpointSeedsMissingT++;
				if(refreshGraphEndpoints){
					if(c.rightCode!=DEAD_END){graphEndCodesChangedT++;}
					c.rightCode=DEAD_END;
					c.rightRatio=0;
					graphEndsRefreshedT++;
					graphEndCodeCountsT[c.rightCode]++;
					missingEndCodeCountsT[c.rightCode]++;
				}
				return;
			}
			int rightMaxPos=fillRightCounts(kmer0, rightCounts);
			int rightMax=rightCounts[rightMaxPos];
			int rightSecondPos=Tools.secondHighestPosition(rightCounts);
			int rightSecond=rightCounts[rightSecondPos];
			if(refreshGraphEndpoints){
				fillLeftCounts(kmer0, leftCounts);
				final int old=c.rightCode;
				c.rightCode=classifyGraphEnd(leftCounts, rightCounts);
				c.rightRatio=graphEndRatio(c.rightCode, leftCounts, rightCounts);
				graphEndsRefreshedT++;
				if(c.rightCode!=old){graphEndCodesChangedT++;}
				graphEndCodeCountsT[c.rightCode]++;
			}
			if(owner!=c.id){
				endpointSeedsAmbiguousT++;
				if(refreshGraphEndpoints){ambiguousEndCodeCountsT[c.rightCode]++;}
				return;
			}
			endpointSeedsUniqueT++;
			if(refreshGraphEndpoints){uniqueEndCodeCountsT[c.rightCode]++;}

			for(int x=0; x<rightCounts.length; x++){
				bb.clear();
				final int count=rightCounts[x];
				int target=-1;
				if(count>0 && isJunction(rightMax, count)){
					if(crossKGraph || refreshGraphEndpoints){traversalAttemptsT++;}
					kmer.setFrom(kmer0);
					kmer.addRightNumeric(x);
					assert(tables.getCount(kmer)==count) : count+", "+tables.getCount(kmer);
					bb.append(AminoAcid.numberToBase[x]);
					target=exploreRight(kmer, leftCounts, extraCounts, bb, c.id);
					if(crossKGraph || refreshGraphEndpoints){exitCountsT[lastExitCondition]++;}
					if(verbose){
						outstream.println(c.id+"R_F: x="+x+", cnt="+count+", dest="+target+", "+codeStrings[lastExitCondition]+", len="+lastLength+", orient="+lastOrientation);
					}
				}
				if(target>=0){
					lastOrientation|=1;
					Edge se=new Edge(c.id, target, lastLength, lastOrientation, count, bb.toBytes());
					c.addRightEdge(se);
					edgesMadeT++;
				}
			}
		}

		/**
		 * Explores rightward extension from k-mer until reaching owned k-mer or junction.
		 * Builds path sequence and determines termination conditions and target contig.
		 *
		 * @param kmer Starting k-mer for exploration
		 * @param leftCounts Buffer for left extension counts
		 * @param rightCounts Buffer for right extension counts
		 * @param bb ByteBuilder for path sequence
		 * @return Owner ID of target contig, or -1 if no valid connection
		 */
		private int exploreRight(Kmer kmer, int[] leftCounts, int[] rightCounts,
				ByteBuilder bb, final int source){
			final Kmer temp=kmerC;
			int length=1;
			int owner=-1;
			int visitedSize=0;
			lastTarget=-1;
			for(; length<crossKMaxLen; length++){
				if((crossKGraph || refreshGraphEndpoints) && !addVisited(kmer, visitedSize++)){
					lastExitCondition=LOOP;
					lastLength=length;
					return -1;
				}
				owner=tables.findOwner(kmer);
				if(crossKGraph && owner==source){
					lastExitCondition=LOOP;
					lastLength=length;
					return -1;
				}else if(owner>=contigs.size()){
					lastExitCondition=BAD_OWNER;
					lastLength=length;
					return -1;
				}
				if(owner>=0 && !crossKGraph){break;}

				final int leftMaxPos=fillLeftCounts(kmer, leftCounts);
				final int leftMax=leftCounts[leftMaxPos];
				final int leftSecondPos=Tools.secondHighestPosition(leftCounts);
				final int leftSecond=leftCounts[leftSecondPos];
				if(isJunction(leftMax, leftSecond)){
					lastExitCondition=B_BRANCH;
					lastLength=length;
					return -1;
				}
				if(owner>=0){break;}
				
				final int rightMaxPos=fillRightCounts(kmer, rightCounts);
				final int rightMax=rightCounts[rightMaxPos];
				final int rightSecondPos=Tools.secondHighestPosition(rightCounts);
				final int rightSecond=rightCounts[rightSecondPos];

//				outstream.println("* "+Arrays.toString(leftCounts)+", "+Arrays.toString(rightCounts)+", "+rightMaxPos);
				
				if(rightMax<minCountExtend){
//					assert(false) : Arrays.toString(rightCounts);
					lastExitCondition=DEAD_END;
					lastLength=length;
					return -1;
				}else if(isJunction(rightMax, rightSecond)){
					lastExitCondition=F_BRANCH;
					lastLength=length;
					return -1;
				}
				bb.append(AminoAcid.numberToBase[rightMaxPos]);
				long x=rightMaxPos;
				kmer.addRightNumeric(x);
			}
			lastLength=length;
			lastTarget=owner;
			if(owner>=0){
				lastExitCondition=SUCCESS;
				Contig dest=contigs.get(owner);
				dest.leftKmer(temp);
				
//				if(temp.equals(kmer)){
//					lastOrientation=temp.sameOrientation(kmer) ? 0 : 1;
//				}else{
//					dest.rightKmer(temp);
//					if(temp.equals(kmer)){
//						lastOrientation=temp.sameOrientation(kmer) ? 2 : 3;
//					}
//				}
				
				if(temp.equals(kmer)){
					lastOrientation=0;
				}else{
					dest.rightKmer(temp);
					if(temp.equals(kmer)){
						lastOrientation=2;
					}else{
						assert(false);
					}
				}
			}else{
				lastExitCondition=TOO_LONG;
			}
			return owner;
		}

		/** Adds one exact oriented kmer state; returns false if it was already visited. */
		private boolean addVisited(final Kmer kmer, final int size){
			final long[] array=kmer.array1();
			final int mult=array.length;
			for(int old=0; old<size; old++){
				final int offset=old*mult;
				int word=0;
				while(word<mult && visited[offset+word]==array[word]){word++;}
				if(word==mult){return false;}
			}
			final int offset=size*mult;
			for(int word=0; word<mult; word++){visited[offset+word]=array[word];}
			return true;
		}

		final Kmer kmerA, kmerB, kmerC;
		final long[] visited;

	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------       Extension Methods      ----------------*/
	/*--------------------------------------------------------------*/


	/**
	 * Estimates insert size between paired reads by measuring k-mer distance.
	 *
	 * @param r1 First read of pair
	 * @param r2 Second read of pair
	 * @param rightCounts Buffer for extension counts
	 * @param kmer1 Working k-mer for first read
	 * @param kmer2 Working k-mer for second read
	 * @return Estimated insert size or -1 if measurement failed
	 */
	public int findInsertSize(Read r1, Read r2, int[] rightCounts, Kmer kmer1, Kmer kmer2){
		kmer1=tables.rightmostKmer(r1.bases, r1.length(), kmer1);
		kmer2=tables.rightmostKmer(r2.bases, r2.length(), kmer2);
		if(kmer1==null || kmer2==null){return -1;}
		final int x=measureInsert(kmer1, kmer2, 24000, rightCounts);
		if(x<0){return -1;}
		return r1.length()+r2.length()+x-kbig;//TODO: May be off by 1.
	}

	/* (non-Javadoc)
	 * @see assemble.Tadpole#extendRead(stream.Read, stream.ByteBuilder, int[], int[], int)
	 */
	@Override
	public int extendRead(Read r, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int distance) {
		return extendRead(r, bb, leftCounts, rightCounts, distance, getLocalKmer());
	}

	@Override
	public int extendRead(Read r, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int distance, final Kmer kmer){
		final int initialLen=r.length();
		if(initialLen<kbig){return 0;}
		bb.setLength(0);
		bb.append(r.bases);
		Kmer temp=tables.rightmostKmer(bb, kmer);
		if(temp==null){return 0;}
		final int extension=extendToRight2_inner(bb, leftCounts, rightCounts, distance, true, kmer);
		if(extension>0){
			r.bases=bb.toBytes();
			if(r.quality!=null){
				final byte q=Shared.FAKE_QUAL;
				r.quality=KillSwitch.copyOf(r.quality, r.bases.length);
				for(int i=initialLen; i<r.quality.length; i++){
					r.quality[i]=q;
				}
			}
		}
		assert(extension==r.length()-initialLen);
		return extension;
	}
	
	/** Returns distance between the two kmers, or -1 */
	public int measureInsert(final Kmer kmer1, final Kmer kmer2, final int maxlen, final int[] rightCounts){
		int len=0;
		
		{
			int count=tables.getCount(kmer2);
			if(count<minCountSeed){return -1;}
		}
		
		int count=tables.getCount(kmer1);
		if(count<minCountSeed){return -1;}
		if(count<minCountSeed){
			if(verbose){outstream.println("Returning because count was too low: "+count);}
			return -1;
		}
		
		int rightMaxPos=fillRightCounts(kmer1, rightCounts);
		int rightMax=rightCounts[rightMaxPos];
//		int rightSecondPos=Tools.secondHighestPosition(rightCounts);
//		int rightSecond=rightCounts[rightSecondPos];
		
		if(rightMax<minCountExtend){return -1;}
//		if(isJunction(rightMax, rightSecond)){return -1;}
		
		while(!kmer1.equals(kmer2) && len<maxlen){
			
			//Generate the new kmer
//			final byte b=AminoAcid.numberToBase[rightMaxPos];
			final long x=rightMaxPos;
			kmer1.addRightNumeric(x);
			
			assert(tables.getCount(kmer1)==rightMax);
			count=rightMax;
			
			assert(count>=minCountExtend) : count;
			
			rightMaxPos=fillRightCounts(kmer1, rightCounts);
			rightMax=rightCounts[rightMaxPos];
//			rightSecondPos=Tools.secondHighestPosition(rightCounts);
//			rightSecond=rightCounts[rightSecondPos];
			
			if(verbose){
				outstream.println("kmer: "+kmer1);
				outstream.println("Counts: "+count+", "+Arrays.toString(rightCounts));
				outstream.println("rightMaxPos="+rightMaxPos);
				outstream.println("rightMax="+rightMax);
//				outstream.println("rightSecondPos="+rightSecondPos);
//				outstream.println("rightSecond="+rightSecond);
			}
			
			if(rightMax<minCountExtend){
				if(verbose){outstream.println("Breaking because highest right was too low:"+rightMax);}
				break;
			}

//			if(isJunction(rightMax, rightSecond)){return -1;}
			
			len++;
		}
		return len>=maxlen ? -1 : len;
	}
	

	
	/**
	 * Extend these bases into a contig.
	 * Stops at both left and right junctions.
	 * Claims ownership.
	 */
	public int extendToRight(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int id, Kmer kmer){
		if(bb.length()<kbig){return BAD_SEED;}
		kmer.clear();
		
		kmer=tables.rightmostKmer(bb, kmer);
		if(kmer==null || kmer.len<kbig){return BAD_SEED;}
		assert(kmer.len==kbig);
		
		/* Now the trailing kmer has been initialized. */
		
		if(verbose){
			outstream.println("extendToRight kmer="+kmer+", bb="+bb);
		}
		
		HashArrayU1D table=tables.getTable(kmer);
		int count=table.getValue(kmer);
		if(count<minCountSeed){
			if(verbose){outstream.println("Returning because count was too low: "+count);}
			return BAD_SEED;
		}
		
		int owner=(useOwnership ? table.getOwner(kmer) : id);
		if(verbose){outstream.println("Owner: "+owner);}
		if(owner>id){return BAD_OWNER;}
		
		int leftMaxPos=0;
		int leftMax=minCountExtend;
		int leftSecondPos=1;
		int leftSecond=0;
		
		if(leftCounts!=null){
			leftMaxPos=fillLeftCounts(kmer, leftCounts);
			leftMax=leftCounts[leftMaxPos];
			leftSecondPos=Tools.secondHighestPosition(leftCounts);
			leftSecond=leftCounts[leftSecondPos];
		}
		
		int rightMaxPos=fillRightCounts(kmer, rightCounts);
		int rightMax=rightCounts[rightMaxPos];
		int rightSecondPos=Tools.secondHighestPosition(rightCounts);
		int rightSecond=rightCounts[rightSecondPos];
		
		if(verbose){
			outstream.println("kmer: "+toText(kmer));
			outstream.println("Counts: "+count+", "+(leftCounts==null ? "null" : Arrays.toString(leftCounts))+", "+Arrays.toString(rightCounts));
			outstream.println("leftMaxPos="+leftMaxPos);
			outstream.println("leftMax="+leftMax);
			outstream.println("leftSecondPos="+leftSecondPos);
			outstream.println("leftSecond="+leftSecond);
			outstream.println("rightMaxPos="+rightMaxPos);
			outstream.println("rightMax="+rightMax);
			outstream.println("rightSecondPos="+rightSecondPos);
			outstream.println("rightSecond="+rightSecond);
		}
		
		if(rightMax<minCountExtend){return DEAD_END;}
		if(isJunction(rightMax, rightSecond)){//Returning here is fine because nothing can be added
			if(verbose){outstream.println("B: Breaking because isJunction("+rightMax+", "+rightSecond+", "+leftMax+", "+leftSecond+")");}
			return isJunction(leftMax, leftSecond) ? D_BRANCH : F_BRANCH;
		}
		if(isJunction(leftMax, leftSecond)){//Returning here is necessary, but this should mean the the length is exactly K
			assert(bb.length()==kbig) : bb.length()+", "+kbig+", "+leftMax+", "+leftSecond;
			if(verbose){outstream.println("B: Breaking because isJunction("+rightMax+", "+rightSecond+", "+leftMax+", "+leftSecond+")");}
			return B_BRANCH;
		}
		
		if(useOwnership){
			owner=table.setOwner(kmer, id);
			if(verbose){outstream.println("A. Owner is now "+id+" for kmer "+kmer);}
			if(owner!=id){
				if(verbose){outstream.println("Returning early because owner was "+owner+" for thread "+id+".");}
				return BAD_OWNER;
			}
		}
		
		final int maxLen=Tools.min((extendRight<0 ? maxContigLen : bb.length()+extendRight), maxContigLen);
		
		while(owner==id && bb.length()<maxLen){
			
			//Generate the new kmer
			final byte b=AminoAcid.numberToBase[rightMaxPos];
			
			//Now consider the next kmer
			final long evicted=kmer.addRightNumeric(rightMaxPos);
			
			table=tables.getTable(kmer);
			
			assert(table.getValue(kmer)==rightMax || rightMax==0);
			count=rightMax;
			
			assert(count>=minCountExtend) : count;

			if(leftCounts!=null){
				leftMaxPos=fillLeftCounts(kmer, leftCounts);
				leftMax=leftCounts[leftMaxPos];
				leftSecondPos=Tools.secondHighestPosition(leftCounts);
				leftSecond=leftCounts[leftSecondPos];
			}
			
			rightMaxPos=fillRightCounts(kmer, rightCounts);
			rightMax=rightCounts[rightMaxPos];
			rightSecondPos=Tools.secondHighestPosition(rightCounts);
			rightSecond=rightCounts[rightSecondPos];
			
			if(verbose){
				outstream.println("kmer: "+toText(kmer));
				outstream.println("Counts: "+count+", "+(leftCounts==null ? "null" : Arrays.toString(leftCounts))+", "+Arrays.toString(rightCounts));
				outstream.println("leftMaxPos="+leftMaxPos);
				outstream.println("leftMax="+leftMax);
				outstream.println("leftSecondPos="+leftSecondPos);
				outstream.println("leftSecond="+leftSecond);
				outstream.println("rightMaxPos="+rightMaxPos);
				outstream.println("rightMax="+rightMax);
				outstream.println("rightSecondPos="+rightSecondPos);
				outstream.println("rightSecond="+rightSecond);
			}
			
			final boolean fbranch=isJunction(rightMax, rightSecond);
			final boolean bbranch=isJunction(leftMax, leftSecond);
			final boolean hbranch=(leftCounts!=null && leftMaxPos!=evicted && branchMult1>0);
			if(bbranch){
				if(verbose){outstream.println("B: Breaking - isJunction("+rightMax+", "+rightSecond+", "+leftMax+", "+leftSecond+"); "
						+ "("+fbranch+", "+bbranch+", "+hbranch+")");}
				return fbranch ? D_BRANCH : B_BRANCH;
			}else if(hbranch){
				if(verbose){outstream.println("B: Breaking - isJunction("+rightMax+", "+rightSecond+", "
						+ ""+leftMax+", "+leftSecond+"); ("+fbranch+", "+bbranch+", "+hbranch+")");}
				if(verbose){outstream.println("Hidden branch: leftMaxPos!=evicted ("+leftMaxPos+"!="+evicted+")" +
						"\nleftMaxPos="+leftMaxPos+", leftMax="+leftMax+", leftSecondPos="+leftSecondPos+", leftSecond="+leftSecond);}
				return fbranch ? D_BRANCH : B_BRANCH;
			}
			
			bb.append(b);
			if(verbose){outstream.println("Added base "+(char)b);}
			
			if(useOwnership){
				owner=table.getOwner(kmer);
				if(verbose){outstream.println("Owner is initially "+id+" for key "+kmer);}
				//[assemble/Tadpole2#001 FIXED 2026-07-02] LONE DIVERGENCE vs Tadpole1: this loop-detection returned LOOP unconditionally, but Tadpole1.extendToRight (canonical, k<=31) returns `fbranch ? F_BRANCH : LOOP`. A forward-branching loop-closure was mis-coded LOOP instead of F_BRANCH, which (a) can trigger trimCircular's kbig-1 trim when both ends read LOOP and (b) suppresses forward-branch handling in default-on popBubbles - either can shift the k>31 contig COUNT by +/-1. owner==id is multithreaded-timing-gated => matched Brian's rare ±1-only-multithreaded-k>31 nondeterminism. Fixed to match Tadpole1 (see return below).
				if(owner==id){//loop detection
					if(verbose  /*|| true*/){
//						outstream.println(new String(bb.array, bb.length()-31, 31));
						outstream.println(bb);
						outstream.println(toText(kmer));
						outstream.println("Breaking because owner was "+owner+" for thread "+id+".");
					}
					return fbranch ? F_BRANCH : LOOP;//was `return LOOP` - a forward-branching loop-closure must be coded F_BRANCH, matching Tadpole1.extendToRight [assemble/Tadpole2#001 FIXED]
				}
				owner=table.setOwner(kmer, id);
				if(verbose){outstream.println("B. Owner is now "+id+" for kmer "+kmer);}
			}
			
			if(fbranch){
				if(verbose){outstream.println("B: Breaking - isJunction("+rightMax+", "+rightSecond+", "+leftMax+", "+leftSecond+"); "
						+ "("+fbranch+", "+bbranch+", "+hbranch+")");}
				return F_BRANCH;
			}else if(rightMax<minCountExtend){
				if(verbose){outstream.println("B: Breaking because highest right was too low:"+rightMax);}
				return DEAD_END;
			}
		}
		assert(owner!=id);
		if(verbose  /*|| true*/){
			outstream.println("Current contig: "+bb+"\nReturning because owner was "+owner+" for thread "+id+".");
		}
		return BAD_OWNER;
	}
	
	@Override
	public int extendToRight2(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int distance, boolean includeJunctionBase){
		initializeThreadLocals();
		return extendToRight2(bb, leftCounts, rightCounts, distance, includeJunctionBase, getLocalKmer());
	}
	
	@Override
	public int extendToRight2(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int distance, boolean includeJunctionBase, Kmer kmer){
		if(verbose || verbose2){outstream.println("Entering extendToRight2 (no kmers).");}
		final int initialLength=bb.length();
		if(initialLength<kbig){return 0;}
		kmer.clear();
		
		kmer=tables.rightmostKmer(bb, kmer);
		if(kmer==null || kmer.len<kbig){return 0;}
		assert(kmer.len==kbig);
		
		return extendToRight2_inner(bb, leftCounts, rightCounts, distance, includeJunctionBase, kmer);
	}
	
	/**
	 * Extend these bases to the right by at most 'distance'.
	 * Stops at right junctions only.
	 * Does not claim ownership.
	 */
	private int extendToRight2_inner(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int distance, boolean includeJunctionBase, Kmer kmer){
		if(verbose || verbose2){outstream.println("Entering extendToRight2_inner (with kmers).");}
		final int initialLength=bb.length();
		assert(kmer.len==kbig) : kmer.len+", "+kbig+", "+bb.length();
		
		HashArrayU1D table=tables.getTable(kmer);
		int count=table.getValue(kmer);
		if(count<minCountSeed){
			if(verbose || verbose2){outstream.println("Returning because count was too low: "+count+"<"+minCountSeed);}
			return 0;
		}
		
		int leftMaxPos=0;
		int leftMax=minCountExtend;
		int leftSecondPos=1;
		int leftSecond=0;
		
		if(leftCounts!=null){
			leftMaxPos=fillLeftCounts(kmer, leftCounts);
			leftMax=leftCounts[leftMaxPos];
			leftSecondPos=Tools.secondHighestPosition(leftCounts);
			leftSecond=leftCounts[leftSecondPos];
		}
		
		int rightMaxPos=fillRightCounts(kmer, rightCounts);
		int rightMax=rightCounts[rightMaxPos];
		int rightSecondPos=Tools.secondHighestPosition(rightCounts);
		int rightSecond=rightCounts[rightSecondPos];
		
		if(verbose){
			outstream.println("kmer: "+toText(kmer));
			outstream.println("Counts: "+count+", "+Arrays.toString(rightCounts));
			outstream.println("rightMaxPos="+rightMaxPos);
			outstream.println("rightMax="+rightMax);
			outstream.println("rightSecondPos="+rightSecondPos);
			outstream.println("rightSecond="+rightSecond);
		}
		
		if(rightMax<minCountExtend){
			if(verbose || verbose2){outstream.println("Returning because rightMax was too low: "+rightMax+"<"+minCountExtend+"\n"+count+", "+Arrays.toString(rightCounts));}
			return 0;
		}
		if(isJunction(rightMax, rightSecond, leftMax, leftSecond)){
			if(verbose || verbose2){outstream.println("Returning because isJunction: "+rightMax+", "+rightSecond+"; "+leftMax+", "+leftSecond);}
			return 0;
		}
		
		final int maxLen=Tools.min(bb.length()+distance, maxContigLen);
		
		while(bb.length()<maxLen){
			
			//Generate the new kmer
			final byte b=AminoAcid.numberToBase[rightMaxPos];
			
			//Now consider the next kmer
			final long evicted=kmer.addRightNumeric(rightMaxPos);
			
			table=tables.getTable(kmer);
			
			assert(table.getValue(kmer)==rightMax || rightMax==0);
			count=rightMax;
			
			assert(count>=minCountExtend) : count;
			
			if(leftCounts!=null){
				leftMaxPos=fillLeftCounts(kmer, leftCounts);
				leftMax=leftCounts[leftMaxPos];
				leftSecondPos=Tools.secondHighestPosition(leftCounts);
				leftSecond=leftCounts[leftSecondPos];
			}
			
			rightMaxPos=fillRightCounts(kmer, rightCounts);
			rightMax=rightCounts[rightMaxPos];
			rightSecondPos=Tools.secondHighestPosition(rightCounts);
			rightSecond=rightCounts[rightSecondPos];
			
			if(verbose){
				outstream.println("kmer: "+toText(kmer));
				outstream.println("Counts: "+count+", "+Arrays.toString(rightCounts));
				outstream.println("rightMaxPos="+rightMaxPos);
				outstream.println("rightMax="+rightMax);
				outstream.println("rightSecondPos="+rightSecondPos);
				outstream.println("rightSecond="+rightSecond);
			}

			if(isJunction(rightMax, rightSecond, leftMax, leftSecond)){
				if(includeJunctionBase && kmer.key()==kmer.array2()){//TODO: Does not work on palindromes.
					bb.append(b);
					if(verbose){outstream.println("Added base "+(char)b);}
				}
				break;
			}
			
			if(leftCounts!=null && leftMaxPos!=evicted){
				if(verbose){outstream.println("B: Breaking because of hidden branch: leftMaxPos!=evicted ("+leftMaxPos+"!="+evicted+")" +
						"\nleftMaxPos="+leftMaxPos+", leftMax="+leftMax+", leftSecondPos="+leftSecondPos+", leftSecond="+leftSecond);}
				if(includeJunctionBase && kmer.key()==kmer.array2()){//TODO: Does not work on palindromes.
					bb.append(b);
					if(verbose){outstream.println("Added base "+(char)b);}
				}
				break;
			}
			
			bb.append(b);
			if(verbose){outstream.println("Added base "+(char)b);}
			
			if(rightMax<minCountExtend){
				if(verbose || verbose2){outstream.println("C: Breaking because highest right was too low: "+rightMax+"<"+minCountExtend);}
				break;
			}
		}
		if(verbose || verbose2){outstream.println("Extended by "+(bb.length()-initialLength));}
		return bb.length()-initialLength;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Junk Detection        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean isJunk(Read r){
		boolean junk=isJunk(r, localRightCounts.get(), getLocalKmer());
		return junk;
	}
	
	@Override
	public boolean isJunk(Read r, final int[] counts, Kmer kmer){
		final int blen=r.length();
		if(blen<kbig){return true;}
		final byte[] bases=r.bases;
		kmer.clearFast();
		assert(kmer.len==0);

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts, to get the leftmost kmer */
		for(int i=0; i<kbig; i++){
			kmer.addRight(bases[i]);
		}
		
		if(kmer.len>=kbig){
			int maxPos=fillLeftCounts(kmer, counts);
			if(counts[maxPos]>0){return false;}
		}
		
		final boolean paired=(r.mateLength()>=kbig);
		int maxDepth=0;
		{
			for(int i=kbig; i<blen; i++){
				kmer.addRight(bases[i]);
				if(kmer.len>=kbig){
					int depth=getCount(kmer);
					if(depth>maxDepth){
						maxDepth=depth;
						if(maxDepth>1 && (!paired || maxDepth>2)){return false;}
					}
				}
			}
		}

		if(kmer.len>=kbig && !paired){
			int maxPos=fillRightCounts(kmer, counts);
			if(counts[maxPos]>0){return false;}
		}
		return true;
	}
	
	@Override
	public boolean hasKmersAtOrBelow(Read r, int tooLow, final float fraction){
		return hasKmersAtOrBelow(r, tooLow, fraction, getLocalKmer());
	}
	
	@Override
	public boolean hasKmersAtOrBelow(Read r, final int tooLow, final float fraction, Kmer kmer){
		final int blen=r.length();
		if(blen<kbig){return true;}
		final byte[] bases=r.bases;
		kmer.clearFast();
		assert(kmer.len==0) : kmer.len+", "+kmer;
		
//		outstream.println("\n"+new String(r.bases)+":");
		
		final int limit=Tools.max(1, Math.round((bases.length-kbig+1)*fraction));
		int valid=0, invalid=0;
		{
			for(int i=0; i<blen; i++){
				kmer.addRight(bases[i]);
				if(kmer.len>=kbig){
					int depth=getCount(kmer);
//					outstream.println("depth="+depth+", kmer="+kmer);
					if(depth>tooLow){valid++;}
					else{
						invalid++;
						if(invalid>=limit){return true;}
					}
				}
			}
		}
		
		//Compensate for nocalls changing the expected kmer count
		final int limit2=Tools.max(1, Math.round((valid+invalid)*fraction));
		return valid<1 || invalid>=limit2;
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------       Error Correction       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public int errorCorrect(Read r){
		initializeThreadLocals();
		int corrected=errorCorrect(r, localLeftCounts.get(), localRightCounts.get(), localIntList.get(), localIntList2.get(),
				localByteBuilder.get(), localByteBuilder2.get(), localTracker.get(), localBitSet.get(), getLocalKmer(), getLocalKmer2());
		return corrected;
	}
	
	@Override
	public int errorCorrect(Read r, final int[] leftCounts, final int[] rightCounts, LongList kmers, IntList counts, IntList counts2,
			final ByteBuilder bb, final ByteBuilder bb2, final ErrorTracker tracker, final BitSet bs, Kmer kmer, Kmer kmer2){
		return errorCorrect(r, leftCounts, rightCounts, counts, counts2, bb, bb2, tracker, bs, kmer, kmer2);
	}
	
	/**
	 * Quick error detection by sampling k-mers at regular intervals.
	 * Looks for coverage drops indicating potential sequencing errors.
	 *
	 * @param bases Sequence to check
	 * @param kmer Working k-mer object
	 * @return true if errors are likely present
	 */
	boolean hasErrorsFast(byte[] bases, Kmer kmer){
		if(bases.length<=kbig){return false;}
		int prev=-99;
		
		kmer.clearFast();
		final int incr=Tools.mid(1, kbig/2, 9), mcc=minCountCorrect();
		for(int i=0, next=kbig-1; i<bases.length; i++){
			kmer.addRight(bases[i]);
			if(i+1>=kbig){
				if(kmer.len()<kbig){
					assert(new String(bases).indexOf('N')>=0);
					return true;
				}
				if(i==next){
					assert(kmer.len()>=kbig);
					int count=getCount(kmer);
					final int min=Tools.min(count, prev), max=Tools.max(count, prev);
					if(prev>-99 && (count<mcc || (i>0 && (isError(max+1, min-1))))){
						return true;
					}
					prev=count;
					next=Tools.min(bases.length-1, next+incr);
				}
			}else{
				assert(kmer.len()<kbig);
			}
		}
		return false;
	}
	
	/**
	 * Main error correction implementation using multiple strategies.
	 * Applies pincer correction, tail correction, and reassembly as needed.
	 * Includes rollback functionality to prevent over-correction.
	 *
	 * @param r Read to correct
	 * @param leftCounts Buffer for left extension counts
	 * @param rightCounts Buffer for right extension counts
	 * @param counts K-mer count list for the read
	 * @param counts2 Additional count buffer
	 * @param bb ByteBuilder for sequence work
	 * @param bb2 Additional ByteBuilder
	 * @param tracker Error tracking object
	 * @param bs BitSet for base marking
	 * @param kmer Working k-mer object
	 * @param regenKmer K-mer for count regeneration
	 * @return Total number of bases corrected
	 */
	public int errorCorrect(Read r, final int[] leftCounts, final int[] rightCounts, IntList counts, IntList counts2,
			final ByteBuilder bb, final ByteBuilder bb2, final ErrorTracker tracker, final BitSet bs, final Kmer kmer, final Kmer regenKmer){
		
//		verbose=r.numericID==0;
		
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		tracker.clear();
		if(!r.containsUndefined() && !hasErrorsFast(bases, kmer)){return 0;}
		
		int valid=tables.fillCounts(bases, counts, kmer);
		if(valid<2){return 0;}
		int possibleErrors=countErrors(counts, quals);
		if(possibleErrors<0){return 0;}
		final float expectedErrors=r.expectedErrors(true, r.length());
		final Rollback roll=ECC_ROLLBACK ? new Rollback(r, counts) : null;
		
		int correctedPincer=0;
		int correctedTail=0;
		int correctedBrute=0;
		int correctedReassemble=0;
		int correctedSubstitute=0;

		if(ECC_PINCER){
			correctedPincer+=errorCorrectPincer(bases, quals, leftCounts, rightCounts, counts, bb, tracker, errorExtensionPincer, kmer);
		}
		
		if(ECC_TAIL || ECC_ALL){
			int start=(ECC_ALL ? 0 : counts.size-kbig-1);
//			if(ECC_PINCER && tracker!=null && tracker.detected>correctedPincer){start=start-kbig;}
			correctedTail+=errorCorrectTail(bases, quals, leftCounts, rightCounts, counts, bb, tracker, start, errorExtensionTail, kmer);
			r.reverseComplement();
			
			counts.reverse();
			correctedTail+=errorCorrectTail(bases, quals, leftCounts, rightCounts, counts, bb, tracker, start, errorExtensionTail, kmer);
			r.reverseComplement();
			counts.reverse();
		}
		
		if(ECC_REASSEMBLE){
			if((correctedPincer<1 && correctedTail<1) || countErrors(counts, quals)>0){
				correctedReassemble=reassemble(bases, quals, rightCounts, counts, counts2, tracker, errorExtensionReassemble, bb, bb2, kmer, regenKmer, bs);
			}
		}

		//Item 3c (Tadpole2 port): substitute runs strictly after reassemble, only as cleanup for
		//errors reassemble left behind -- mirrors Tadpole1.errorCorrect exactly.
		if(ECC_SUBSTITUTE && hasLowCount(counts)){
			correctedSubstitute=errorCorrectSubstitute(bases, quals, counts, counts2, bb, tracker, kmer, regenKmer, bs);
		}

		assert(correctedPincer+correctedTail+correctedReassemble+correctedBrute+correctedSubstitute==tracker.corrected())
			: correctedPincer+", "+correctedTail+", "+correctedReassemble+", "+correctedBrute+", "+correctedSubstitute+", "+tracker;

		if(ECC_ROLLBACK && (tracker.corrected()>0 || tracker.rollback)){
			
			if(!tracker.rollback && quals!=null && tracker.corrected()>3){
				float mult=Tools.max(1, 0.5f*(0.5f+0.01f*r.length()));//1 for a 150bp read.
				if(countErrors(counts, quals)>0 && tracker.corrected()>mult+expectedErrors){tracker.rollback=true;}
				else if(tracker.corrected()>2.5f*mult+expectedErrors){tracker.rollback=true;}
				//Item 3c (diagnostic, Amber directing): Trigger 1 = over-correction heuristic.
				if(tracker.rollback){tracker.rollbackTrigger=1;}
			}

//			boolean printed=false;
			IntList counts0=roll.counts0;
			for(int i=0; !tracker.rollback && i<counts.size; i++){
				int a=Tools.max(0, counts0.get(i));
				int b=Tools.max(0, counts.get(i));
//				assert(b+1>=a) : "Z: RID="+r.numericID+"; "+a+"->"+b+"\n"+counts0+"\n"+counts;
				if(b<a-1 && !isSimilar(a, b)){
//					assert(false) : "Y: RID="+r.numericID+"; "+a+"->"+b+"\n"+counts0+"\n"+counts;
					if(verbose){outstream.println("Y: RID="+r.numericID+"; "+a+"->"+b+"\n"+counts0+"\n"+counts);}
					tracker.rollback=true;
					//Item 3c (diagnostic): Trigger 2 = count-contradiction check.
					tracker.rollbackTrigger=2;
				}
//				else if(b<a-1 && !printed){
//					assert(false);
//					if(verbose){outstream.println("X: RID="+r.numericID+"; "+a+"->"+b+"\n"+counts0+"\n"+counts);}
//					printed=true;
//				}
			}

			if(tracker.rollback){
				//Item 3c: snapshot before clearCorrected() wipes it -- feeds the rollback
				//corrected-count histogram and the substitute/non-substitute rollback split.
				tracker.rollbackCorrectedSnapshot=tracker.corrected();
				tracker.rollbackSubstitute=(correctedSubstitute>0);
				roll.rollback(r, counts);
				tracker.clearCorrected();
				return 0;
			}
		}
		
		if(MARK_BAD_BASES>0 && (!MARK_ERROR_READS_ONLY || countErrors(counts, quals)>0 ||
				r.expectedErrors(false, r.length())>3)){
			int marked=markBadBases(bases, quals, counts, bs, MARK_BAD_BASES, MARK_DELTA_ONLY, MARK_QUALITY);
			tracker.marked=marked;
		}
		
		return tracker.corrected();
	}

	/** Item 3c (Tadpole2 port): true iff every kmer spanning base position p has a low count
	 * (below minCountCorrect). Mirrors Tadpole1's version, using kbig instead of k.
	 * @see assemble.Tadpole1#isErrorCovered */
	private boolean isErrorCovered(final int p, final IntList counts){
		final int lo=Tools.max(0, p-kbig+1);
		final int hi=Tools.min(p, counts.size-1);
		for(int i=lo; i<=hi; i++){
			if(counts.get(i)>=minCountCorrect()){return false;}
		}
		return true;
	}

	/** Item 3c (Fix #3, Noelle's perf review, Amber directing): builds the kmer spanning
	 * bases[windowStart..windowStart+kbig-1] ONCE, with bases[p] forced to a placeholder
	 * base (0/A) rather than its true original value -- so validity (the returned boolean)
	 * reflects only OTHER positions in the window, never p's own identity, exactly matching
	 * the old getCountWithSubstitution's behavior of NEVER looking up a kmer containing p's
	 * true original base (it was always overwritten with a real alternateBase before lookup).
	 * The caller then uses {@link Kmer#substituteBase} to swap position p (window-relative
	 * offset p-windowStart) between the placeholder and each real alternate for O(1) lookups,
	 * instead of rebuilding the whole O(kbig) window per alternate (the actual perf win --
	 * this method now runs once per span per position, not once per (span,alternate) pair). */
	private boolean buildWindow(final byte[] bases, final int windowStart, final int p, final Kmer scratch){
		scratch.clear();
		for(int i=windowStart; i<windowStart+kbig; i++){
			final byte b=(i==p ? AminoAcid.numberToBase[0] : bases[i]);
			scratch.addRight(b);
		}
		return scratch.len>=kbig;
	}

	/** Item 3c (Fix #4): tests `testBase` at one already-built span window. Returns a packed
	 * int: bit0 set if the count clears minCountCorrect() (a "confirmation"), bit1 set if the
	 * count exceeds the position's original depth (Brian's detectedSubstitute rule). If
	 * `valid` is false (buildWindow found an N elsewhere in the window), contributes 0 without
	 * touching `scratch` at all -- matches the original getCountWithSubstitution's behavior of
	 * never producing a confirmation for an N-containing window. Shared by both the per-
	 * position scan (called once per alternate base) and the post-batch validation (called
	 * once for the already-chosen base) -- kept as one method specifically so the two phases
	 * can never silently diverge in how a span is scored. */
	private int confirmSpan(final Kmer scratch, final int posInWindow, final int testBase, final boolean valid, final int windowStart, final IntList counts){
		if(!valid){return 0;}
		final long placeholder=scratch.substituteBase(posInWindow, testBase);
		assert(placeholder==0) : "buildWindow's placeholder invariant violated: "+placeholder;
		final int count=Tools.max(0, getCount(scratch));
		final long restored=scratch.substituteBase(posInWindow, placeholder);
		assert(restored==testBase) : "restore mismatch: expected "+testBase+", got "+restored;
		int r=0;
		if(count>=minCountCorrect()){r|=1;}
		if(count>Tools.max(0, counts.get(windowStart))){r|=2;}
		return r;
	}

	/** Item 3c (Tadpole2 port, 2026-08-23, Brian green-lit, Amber directing, Nepgear
	 * implementing): substitute ECC for k&gt;31. Faithful mechanism port of
	 * {@link Tadpole1#errorCorrectSubstitute} -- NO new heuristics or threshold tuning (Brian's
	 * explicit scope boundary; ratio-threshold work is his to do). Same algorithm: for each
	 * error-covered position, test all 3 single-base substitutions against up to 4
	 * reflection-closed spanning kmers, accept a correction only if exactly one alternate base
	 * is independently confirmed (count &gt;= minCountCorrect) by at least CONFIRM_MIN spans
	 * (1 here -- exact table, no collision risk), 0 or &gt;1 clearing the bar means skip as
	 * ambiguous. detectedSubstitute counts positions where some alternate base's confirmed depth exceeds
	 * the original's at the same window (Brian's ruling, ported verbatim from Tadpole1). Bounded
	 * re-pass (SUBSTITUTE_MAX_PASSES) resolves the "trapped" two-close-errors case.
	 * @see Tadpole1#errorCorrectSubstitute
	 *
	 * Fix #2 (Noelle's perf review, Amber directing): span selection is inlined as scalar
	 * locals instead of allocating a fresh IntList per position per pass -- the old
	 * selectSpanPositions helper method is gone. Zero allocation.
	 *
	 * Fix #3 (Noelle's perf review, Amber directing): each span's window is built ONCE per
	 * position (via {@link #buildWindow}), not once per (span,alternate) pair.
	 *
	 * Fix #4 (Noelle's finding: forward vs reverse-complement of the same read produced
	 * different corrections; Amber directing; design reviewed by both before implementation):
	 * two independent RC-dependence sources, both fixed here.
	 * (a) Reflection-closed span selection: the old single `mid=(lo+hi)/2` floors when
	 *     `hi-lo` is odd, so the RC-mirrored window lands on the ceiling, not the floor -- an
	 *     asymmetric choice under reversal. Now: `d=hi-lo` (still gated on `d>=4`); `d` even
	 *     keeps the single center `lo+d/2` (self-symmetric, unchanged); `d` odd uses BOTH
	 *     `lo+d/2` (floor) and `lo+d/2+1` (ceil) -- exact reflections of each other, so the
	 *     span SET is closed under reversal even though neither member is. `kmerCeilMid`
	 *     (another lazily-initialized per-thread scratch, {@link #getSubstituteKmerCeilMid()})
	 *     holds this fourth span when present.
	 * (b) Snapshot pass, batched proposal application: the old code committed each accepted
	 *     correction immediately (`bases[p]=...`, then `regenerateCounts`) before continuing
	 *     to p+1 in the SAME pass, so an early correction could change the evidence available
	 *     to later positions -- and "later" means something different depending on scan
	 *     direction. Now: the scan is non-mutating (bases/quals only read, against the
	 *     pass-start `counts` snapshot); accepted positions are staged in `bs` (position) and
	 *     `counts2` (chosen base, repurposed here as a per-BASE-position buffer, distinct from
	 *     its normal per-count-window role in reassemble -- verified free, reassemble is done
	 *     with both `counts2` and `bs` by the time substitute runs). If `bs` ends up empty
	 *     after the scan, break immediately (no proposals, no progress, matches the existing
	 *     zero-progress exit). Otherwise: snapshot `bases` into the already-per-thread
	 *     `ByteBuilder bb` (`bb.clear(); bb.append(bases);` -- reassemble is finished with it
	 *     too, so no new ThreadLocal state is needed; quals need no snapshot at all since
	 *     they're never written before validation succeeds), apply every staged proposal,
	 *     regenerate counts ONCE via the batch BitSet form, then re-validate ONLY the
	 *     already-chosen base per proposed position (via {@link #confirmSpan}, not a fresh
	 *     alternate search) against the new counts. If every proposal still clears
	 *     CONFIRM_MIN: apply qualities now, keep the batch. If ANY proposal fails: restore
	 *     every touched base from the `bb` snapshot (byte-exact -- preserves N/case with no
	 *     special-casing), regenerate counts AGAIN against the reverted bases, and stop this
	 *     read's substitute pass entirely (whole-pass revert, Option A -- trivially
	 *     RC-invariant by construction since revert-all is symmetric; a connected-component
	 *     alternative is held in reserve only if the aggregate value gate shows this costs
	 *     real corrections in practice). */
	private int errorCorrectSubstitute(final byte[] bases, final byte[] quals, final IntList counts, final IntList counts2, final ByteBuilder bb, final ErrorTracker tracker, final Kmer kmer, final Kmer regenKmer, final BitSet bs){
		int totalCorrected=0;
		final int readLen=bases.length;
		final Kmer kmerHi=getSubstituteKmerHi();
		final Kmer kmerMid=getSubstituteKmerMid();
		final Kmer kmerCeilMid=getSubstituteKmerCeilMid();
		bs.clear();
		counts2.clear();

		for(int pass=0; pass<SUBSTITUTE_MAX_PASSES; pass++){

			for(int p=0; p<readLen; p++){
				if(!isErrorCovered(p, counts)){continue;}

				final int lo=Tools.max(0, p-kbig+1);
				final int hi=Tools.min(p, readLen-kbig);
				if(hi<=lo){continue;}//Can't confirm with fewer than 2 spanning kmers
				final int d=hi-lo;
				final int floorMid=(d>=4) ? lo+d/2 : -1;
				final int ceilMid=(d>=4 && (d&1)!=0) ? floorMid+1 : -1;

				final boolean loValid=buildWindow(bases, lo, p, kmer);
				final boolean hiValid=buildWindow(bases, hi, p, kmerHi);
				final boolean floorMidValid=(floorMid>=0) && buildWindow(bases, floorMid, p, kmerMid);
				final boolean ceilMidValid=(ceilMid>=0) && buildWindow(bases, ceilMid, p, kmerCeilMid);
				final int posLo=p-lo, posHi=p-hi;
				final int posFloorMid=(floorMid>=0 ? p-floorMid : -1);
				final int posCeilMid=(ceilMid>=0 ? p-ceilMid : -1);

				final int origBase=AminoAcid.baseToNumber[bases[p]];
				int bestBase=-1, bestConfirmations=0, clearedCount=0;
				boolean anyIncrease=false;
				for(int alternateBase=0; alternateBase<4; alternateBase++){
					if(alternateBase==origBase){continue;}
					int confirmations=0, r;

					r=confirmSpan(kmer, posLo, alternateBase, loValid, lo, counts);
					confirmations+=(r&1); if((r&2)!=0){anyIncrease=true;}

					r=confirmSpan(kmerHi, posHi, alternateBase, hiValid, hi, counts);
					confirmations+=(r&1); if((r&2)!=0){anyIncrease=true;}

					if(floorMid>=0){
						r=confirmSpan(kmerMid, posFloorMid, alternateBase, floorMidValid, floorMid, counts);
						confirmations+=(r&1); if((r&2)!=0){anyIncrease=true;}
					}
					if(ceilMid>=0){
						r=confirmSpan(kmerCeilMid, posCeilMid, alternateBase, ceilMidValid, ceilMid, counts);
						confirmations+=(r&1); if((r&2)!=0){anyIncrease=true;}
					}

					if(confirmations>=CONFIRM_MIN){
						clearedCount++;
						if(confirmations>bestConfirmations){
							bestBase=alternateBase;
							bestConfirmations=confirmations;
						}
					}
				}
				if(anyIncrease){tracker.detectedSubstitute++;}

				if(clearedCount==1){//Exactly one alternate base confirmed; 0 or >1 is ambiguous, skip
					bs.set(p);
					counts2.set(p, bestBase);
				}
			}

			if(bs.isEmpty()){break;}//No proposals this pass -- later passes can't help either

			bb.clear();
			bb.append(bases);
			for(int p=bs.nextSetBit(0); p>=0; p=bs.nextSetBit(p+1)){
				bases[p]=AminoAcid.numberToBase[counts2.get(p)];
			}
			tables.regenerateCounts(bases, counts, regenKmer, bs);

			boolean allValid=true;
			for(int p=bs.nextSetBit(0); allValid && p>=0; p=bs.nextSetBit(p+1)){
				final int lo=Tools.max(0, p-kbig+1);
				final int hi=Tools.min(p, readLen-kbig);
				final int d=hi-lo;
				final int floorMid=(d>=4) ? lo+d/2 : -1;
				final int ceilMid=(d>=4 && (d&1)!=0) ? floorMid+1 : -1;

				final boolean loValid=buildWindow(bases, lo, p, kmer);
				final boolean hiValid=buildWindow(bases, hi, p, kmerHi);
				final boolean floorMidValid=(floorMid>=0) && buildWindow(bases, floorMid, p, kmerMid);
				final boolean ceilMidValid=(ceilMid>=0) && buildWindow(bases, ceilMid, p, kmerCeilMid);

				final int chosenBase=counts2.get(p);
				int confirmations=0, r;
				r=confirmSpan(kmer, p-lo, chosenBase, loValid, lo, counts); confirmations+=(r&1);
				r=confirmSpan(kmerHi, p-hi, chosenBase, hiValid, hi, counts); confirmations+=(r&1);
				if(floorMid>=0){ r=confirmSpan(kmerMid, p-floorMid, chosenBase, floorMidValid, floorMid, counts); confirmations+=(r&1); }
				if(ceilMid>=0){ r=confirmSpan(kmerCeilMid, p-ceilMid, chosenBase, ceilMidValid, ceilMid, counts); confirmations+=(r&1); }

				if(confirmations<CONFIRM_MIN){allValid=false;}
			}

			if(allValid){
				int passCorrected=0;
				for(int p=bs.nextSetBit(0); p>=0; p=bs.nextSetBit(p+1)){
					byte q=(quals==null ? 30 : quals[p]);
					q=(byte)Tools.mid(q+qIncreasePincer, qMinPincer, qMaxPincer);
					if(quals!=null){quals[p]=q;}
					passCorrected++;
				}
				totalCorrected+=passCorrected;
				bs.clear();
			}else{
				for(int p=bs.nextSetBit(0); p>=0; p=bs.nextSetBit(p+1)){
					bases[p]=bb.array[p];
				}
				tables.regenerateCounts(bases, counts, regenKmer, bs);
				bs.clear();
				break;
			}
		}

		tracker.correctedSubstitute+=totalCorrected;
		return totalCorrected;
	}

	/**
	 * Corrects errors using pincer approach with flanking k-mer validation.
	 * Identifies substitutions by comparing left and right k-mer counts around gaps.
	 *
	 * @param bases Sequence bases to correct
	 * @param quals Quality scores for bases
	 * @param leftBuffer Buffer for left extension counts
	 * @param rightBuffer Buffer for right extension counts
	 * @param counts K-mer counts for sequence
	 * @param bb ByteBuilder for extension work
	 * @param tracker Error tracking object
	 * @param errorExtension Extension distance for validation
	 * @param kmer Working k-mer object
	 * @return Number of bases corrected
	 */
	public int errorCorrectPincer(final byte[] bases, final byte[] quals, final int[] leftBuffer, final int[] rightBuffer,
			final IntList counts, final ByteBuilder bb, final ErrorTracker tracker, final int errorExtension, final Kmer kmer){
		
		int detected=0;
		int corrected=0;
		
		//a is the index of the left kmer
		//b is a+1
		//c is d-1
		//d is the index of the right kmer
		//the base between the kmers is at a+k
		for(int a=0, d=kbig+1; d<counts.size; a++, d++){
			final int aCount=counts.get(a);
			final int bCount=counts.get(a+1);
			final int cCount=counts.get(d-1);
			final int dCount=counts.get(d);
			final byte qb=(quals==null ? 20 : quals[a+kbig]);
			if(isError(aCount, bCount, qb) && isError(dCount, cCount, qb) && isSimilar(aCount, dCount)){
				if(verbose){
					outstream.println("Found error: "+aCount+", "+bCount+", "+cCount+", "+dCount);
				}
				//Looks like a 1bp substitution; attempt to correct.
				detected++;
				int ret=correctSingleBasePincer(a, d, bases, quals, leftBuffer, rightBuffer, counts, bb, errorExtension, kmer);
				corrected+=ret;
				if(verbose){
					outstream.println("Corrected error.");
				}
			}else{
				if(verbose){
					outstream.println("Not an error: "+aCount+", "+bCount+", "+cCount+", "+dCount+
							";  "+isError(aCount, bCount, qb)+", "+isError(dCount, cCount, qb)+", "+isSimilar(aCount, dCount));
				}
			}
		}
		
//		if(detected==0 && counts.get(0)>2 && counts.get(counts.size-1)>2){
//			assert(!verbose);
//			verbose=true;
//			outstream.println("\n"+counts);
//			errorCorrectPincer(bases, quals, leftBuffer, rightBuffer, kmers, counts, bb, tracker);
//			assert(false);
//		}
		
		{
			tracker.detectedPincer+=detected;
			tracker.correctedPincer+=corrected;
		}
		
		return corrected;
	}

	/**
	 * Corrects errors from sequence ends using tail-based approach.
	 * Validates corrections by extending and comparing with original sequence.
	 *
	 * @param bases Sequence bases to correct
	 * @param quals Quality scores for bases
	 * @param leftBuffer Buffer for left extension counts
	 * @param rightBuffer Buffer for right extension counts
	 * @param counts K-mer counts for sequence
	 * @param bb ByteBuilder for extension work
	 * @param tracker Error tracking object
	 * @param startPos Starting position for correction scanning
	 * @param errorExtension Extension distance for validation
	 * @param kmer Working k-mer object
	 * @return Number of bases corrected
	 */
	public int errorCorrectTail(final byte[] bases, final byte[] quals, final int[] leftBuffer, final int[] rightBuffer,
			final IntList counts, final ByteBuilder bb, final ErrorTracker tracker, final int startPos, final int errorExtension, final Kmer kmer){
		if(bases.length<kbig+2+errorExtension+deadZone){return 0;}
		int detected=0;
		int corrected=0;
		
		//a is the index of the left kmer
		//b is a+1
		//the base between the kmers is at a+k
		for(int a=Tools.max(startPos, errorExtension), lim=counts.size-deadZone-1; a<lim; a++){//errorExtension-1
			final int aCount=counts.get(a);
			final int bCount=counts.get(a+1);
			final byte qb=(quals==null ? 20 : quals[a+kbig]);
			if(isError(aCount, bCount, qb) && isSimilar(aCount, a-errorExtension, a-1, counts) && isError(aCount, a+2, a+kbig, counts)){
				if(verbose){
					outstream.println("Found error: "+aCount+", "+bCount);
				}
				//Assume like a 1bp substitution; attempt to correct.
				detected++;
				int ret=correctSingleBaseRight(a, bases, quals, leftBuffer, rightBuffer, counts, bb, errorExtension, kmer);
				corrected+=ret;
				if(verbose){
					outstream.println("Corrected error.");
				}
			}else{
				if(verbose){
					outstream.println("Not an error: "+aCount+", "+bCount+
							";  "+isError(aCount, bCount, qb)+", "+isSimilar(aCount, a-errorExtension, a-1, counts)+", "+isError(aCount, a+2, a+kbig, counts));
				}
			}
		}
		
//		if(detected==0 && counts.get(0)>2 && counts.get(counts.size-1)>2){
//			assert(!verbose);
//			verbose=true;
//			outstream.println("\n"+counts);
//			errorCorrectPincer(bases, quals, leftBuffer, rightBuffer, kmers, counts, bb, tracker);
//			assert(false);
//		}
		
		{
			tracker.detectedTail+=detected;
			tracker.correctedTail+=corrected;
		}
		
		return corrected;
	}
	
	@Override
	public int reassemble_inner(final ByteBuilder bb, final byte[] quals, final int[] rightCounts, final IntList counts,
			final int errorExtension, final Kmer kmer, final Kmer regenKmer){
		final int length=bb.length();
		if(length<kbig+1+deadZone){return 0;}
		final byte[] bases=bb.array;
		int detected=0;
		int corrected=0;
		
		//Initialize the kmer
		kmer.clear();
		
		//a is the index of the first base of the left kmer
		//b=a+1 is the index of the next base
		//ca=a-kbig is the index of the next base
		//cb=a-kbig is the index of the next base
		//the base between the kmers is at a+k
		for(int a=0, lim=length-deadZone-1; a<lim; a++){
			//Update the kmer
			kmer.addRight(bases[a]);
			
			if(verbose){
				outstream.println("kmer.len(): "+kmer.len()+" vs "+kbig+"; a="+a);
			}
			
			if(kmer.len()>=kbig){
				
				final int b=a+1;
				final int ca=a-kbig+1;
				final int cb=ca+1;
				
				final int aCount=counts.get(ca);
				final int bCount=counts.get(cb);
				final byte qb=(quals==null ? 20 : quals[b]);

				if(verbose){
					outstream.println("ca="+ca+", cb="+cb+"; aCount="+aCount+", bCount="+bCount);
					outstream.println(isError(aCount, bCount, qb)+", "+isSimilar(aCount, ca-errorExtension, ca-1, counts)+
							", "+isError(aCount, ca+2, ca+kbig, counts));
				}
				
//				if(isError(aCount, bCount) && isSimilar(aCount, ca-errorExtension, ca-1, counts) && isError(aCount, ca+2, ca+kbig, counts)){
				if(isSubstitution(ca, errorExtension, qb, counts)){
					if(verbose){
						outstream.println("***Found error: "+aCount+", "+bCount);
					}
					//Assume like a 1bp substitution; attempt to correct.

					int rightMaxPos=fillRightCounts(kmer, rightCounts);
					int rightMax=rightCounts[rightMaxPos];
					int rightSecondPos=Tools.secondHighestPosition(rightCounts);
					int rightSecond=rightCounts[rightSecondPos];

					byte base=bases[b];
					byte num=AminoAcid.baseToNumber[base];
					
					if(rightMax>=minCountExtend){
						detected++;
						if(num==rightMax){
							detected--;
//							bases2[b]=base;
						}else if((isError(rightMax, rightSecond, qb) || !isJunction(rightMax, rightSecond)) && isSimilar(aCount, rightMax)){
							bases[b]=AminoAcid.numberToBase[rightMaxPos];
							corrected++;
							tables.regenerateCounts(bases, counts, ca, regenKmer);
							if(verbose){outstream.println("Corrected error: "+num+"->"+rightMaxPos+". New counts:\n"+counts);}
						}
						
//						else if(rightSecond>=minCountExtend && isJunction(rightMax, rightSecond) && isSimilar(aCount, rightSecond)
//								&& !isSimilar(aCount, rightMax)){//This branch may not be very safe.
//							bases2[b]=AminoAcid.numberToBase[rightSecondPos];
//							corrected++;
//							if(verbose){outstream.println("Corrected error.");}
//						}
					}
					
				}else{
					if(verbose){
						outstream.println("Not an error: "+aCount+", "+bCount+
								";  "+isError(aCount, bCount, qb)+", "+isSimilar(aCount, a-errorExtension, a-1, counts)+", "+isError(aCount, a+2, a+kbig, counts));
					}
				}
			}
		}
		
		return corrected;
	}
	
	/**
	 * Corrects single base substitution using bidirectional validation.
	 * Extends from both flanking k-mers to confirm correction consistency.
	 *
	 * @param a Left k-mer position
	 * @param d Right k-mer position
	 * @param bases Sequence to correct
	 * @param quals Quality scores
	 * @param leftBuffer Buffer for left extensions
	 * @param rightBuffer Buffer for right extensions
	 * @param counts K-mer count list
	 * @param bb ByteBuilder for extensions
	 * @param errorExtension Validation extension distance
	 * @param kmer0 Working k-mer object
	 * @return 1 if corrected, 0 if correction failed
	 */
	private int correctSingleBasePincer(final int a, final int d, final byte[] bases, final byte[] quals, final int[] leftBuffer, final int[] rightBuffer,
			final IntList counts, final ByteBuilder bb, final int errorExtension, final Kmer kmer0){
		final byte leftReplacement, rightReplacement;
		final int loc=a+kbig;
		{
			bb.clear();
			Kmer kmer=getKmer(bases, a, kmer0);
			if(kmer==null){return 0;}
			int extension=extendToRight2_inner(bb, null, rightBuffer, errorExtension, true, kmer);
			if(extension<errorExtension){return 0;}
			for(int i=1; i<extension; i++){
				if(bb.get(i)!=bases[loc+i]){return 0;}
			}
			leftReplacement=bb.get(0);
		}
		{
			bb.clear();
			Kmer kmer=getKmer(bases, d, kmer0);
			if(kmer==null){return 0;}
			kmer.rcomp();
			int extension=extendToRight2_inner(bb, null, rightBuffer, errorExtension, true, kmer);
			if(extension<errorExtension){return 0;}
			bb.reverseComplementInPlace();
			for(int i=0; i<extension-1; i++){
				if(bb.get(i)!=bases[loc+i+1-extension]){return 0;}
			}
			rightReplacement=bb.get(extension-1);
		}
		if(leftReplacement!=rightReplacement){return 0;}
		if(bases[loc]==leftReplacement){return 0;}
		if(!isSimilar(bases, a, leftReplacement, counts, kmer0)){return 0;}
		
		bases[loc]=leftReplacement;
		assert(d==a+kbig+1);
		tables.regenerateCounts(bases, counts, a, kmer0);
		return 1;
	}
	
	/**
	 * Corrects single base substitution using rightward extension validation.
	 * Extends from left k-mer and validates against original sequence.
	 *
	 * @param a Left k-mer position
	 * @param bases Sequence to correct
	 * @param quals Quality scores
	 * @param leftBuffer Buffer for left extensions
	 * @param rightBuffer Buffer for right extensions
	 * @param counts K-mer count list
	 * @param bb ByteBuilder for extensions
	 * @param errorExtension0 Maximum extension distance
	 * @param kmer0 Working k-mer object
	 * @return 1 if corrected, 0 if correction failed
	 */
	private int correctSingleBaseRight(final int a, final byte[] bases, final byte[] quals, final int[] leftBuffer, final int[] rightBuffer,
			final IntList counts, final ByteBuilder bb, final int errorExtension0, final Kmer kmer0){
		final byte leftReplacement;
		final int loc=a+kbig;
		final int errorExtension=Tools.min(errorExtension0, bases.length-loc);
		{
			bb.clear();
			Kmer kmer=getKmer(bases, a, kmer0);
			if(kmer==null){return 0;}
			int extension=extendToRight2_inner(bb, null, rightBuffer, errorExtension, true, kmer);
			if(extension<errorExtension){return 0;}
			for(int i=1; i<extension; i++){
				if(bb.get(i)!=bases[loc+i]){
					return 0;
				}
			}
			leftReplacement=bb.get(0);
		}
		
		if(bases[loc]==leftReplacement){return 0;}
		if(!isSimilar(bases, a, leftReplacement, counts, kmer0)){return 0;}
		
		bases[loc]=leftReplacement;
		tables.regenerateCounts(bases, counts, a, kmer0);
		return 1;
	}
	
	/**
	 * Checks if replacing a base results in similar k-mer coverage.
	 * Used to validate error corrections by coverage consistency.
	 *
	 * @param bases Sequence containing k-mer
	 * @param a Position of k-mer start
	 * @param newBase Replacement base to test
	 * @param counts K-mer count list
	 * @param kmer0 Working k-mer object
	 * @return true if coverage is similar with new base
	 */
	private final boolean isSimilar(byte[] bases, int a, byte newBase, IntList counts, final Kmer kmer0){
		Kmer kmer=getKmer(bases, a, kmer0);
		if(kmer==null){
			assert(false); //Should never happen
			return false;
		}
		kmer.addRight(newBase);
		int count=getCount(kmer);
		int aCount=counts.get(a);
		boolean similar=isSimilar(aCount, count);
		return similar;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------  Inherited Abstract Methods  ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	final void makeKhist(){
		tables.makeKhist(outHist, histColumns, histMax, histHeader, histZeros, true, smoothHist, gcHist, false, 0.01, 1, 1);
	}
	@Override
	final void dumpKmersAsText(){
		tables.dumpKmersAsBytes_MT(outKmers, minToDump, maxToDump, true, null);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final KmerTableSetU tables(){return tables;}
	/**
	 * Unified k-mer table set for storing and accessing k-mer counts and ownership
	 */
	public final KmerTableSetU tables;
	
	/** Normal kmer length */
	final int ksmall;

	/** Experimental: short-k ownership represents longer-k assembled contigs. */
	boolean crossKGraph=false;

	/** Item 3c (Fix #3): lazily-initialized per-thread scratch Kmer for errorCorrectSubstitute's
	 * "hi" span -- see {@link #getSubstituteKmerHi()}. Tadpole2-private and independent of the
	 * shared initializeThreadLocals()/getLocalKmer() mechanism (which Tadpole.ExtendThread, the
	 * real hot-path caller, never invokes), so it works correctly regardless of which entry
	 * path reaches errorCorrectSubstitute, with no change to Tadpole1 or the shared abstract
	 * errorCorrect signature. */
	private final ThreadLocal<Kmer> substituteKmerHi=new ThreadLocal<Kmer>();
	/** Item 3c (Fix #3): same as {@link #substituteKmerHi} but for the "mid" span. */
	private final ThreadLocal<Kmer> substituteKmerMid=new ThreadLocal<Kmer>();
	/** Item 3c (Fix #4): same as {@link #substituteKmerHi} but for the "ceil-mid" span --
	 * reflection-closed span selection needs a second, distinct central window whenever
	 * hi-lo is odd (see errorCorrectSubstitute's span-geometry comment). */
	private final ThreadLocal<Kmer> substituteKmerCeilMid=new ThreadLocal<Kmer>();

	/** Returns this thread's scratch Kmer for errorCorrectSubstitute's "hi" span, creating it
	 * (one `new Kmer(kbig)`) on this thread's first call only. */
	private Kmer getSubstituteKmerHi(){
		Kmer k=substituteKmerHi.get();
		if(k==null){k=new Kmer(kbig); substituteKmerHi.set(k);}
		return k;
	}

	/** Returns this thread's scratch Kmer for errorCorrectSubstitute's "mid" span, creating it
	 * (one `new Kmer(kbig)`) on this thread's first call only. */
	private Kmer getSubstituteKmerMid(){
		Kmer k=substituteKmerMid.get();
		if(k==null){k=new Kmer(kbig); substituteKmerMid.set(k);}
		return k;
	}

	/** Returns this thread's scratch Kmer for errorCorrectSubstitute's "ceil-mid" span,
	 * creating it (one `new Kmer(kbig)`) on this thread's first call only. */
	private Kmer getSubstituteKmerCeilMid(){
		Kmer k=substituteKmerCeilMid.get();
		if(k==null){k=new Kmer(kbig); substituteKmerCeilMid.set(k);}
		return k;
	}

}
