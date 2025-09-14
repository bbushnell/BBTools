package assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.concurrent.atomic.AtomicInteger;

import dna.AminoAcid;
import jgi.BBMerge;
import kmer.AbstractKmerTable;
import kmer.HashArray1D;
import kmer.HashForest;
import kmer.KmerNode;
import kmer.KmerTableSet;
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
import ukmer.Kmer;


/**
 * Short-kmer assembler based on KmerCountExact.
 * @author Brian Bushnell
 * @date May 15, 2015
 *
 */
public class Tadpole1 extends Tadpole {
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer(), t2=new Timer();
		t.start();
		t2.start();
		
		//Create a new CountKmersExact instance
		Tadpole1 wog=new Tadpole1(args, true);
		t2.stop();
		outstream.println("Initialization Time:      \t"+t2);
		
		///And run it
		wog.process(t);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public Tadpole1(String[] args, boolean setDefaults){
		super(args, setDefaults);
		
		final int bytesPerKmer;
		{
			int mult=12;
			if(useOwnership){mult+=4;}
			if(processingMode==correctMode || processingMode==discardMode){}
			else if(processingMode==contigMode || processingMode==extendMode){mult+=1;}
			bytesPerKmer=mult;
		}
		
		tables=new KmerTableSet(args, bytesPerKmer);
		k=tables.k;
		k2=tables.k2;
		
		shift=2*k;
		shift2=shift-2;
		mask=(shift>63 ? -1L : ~((-1L)<<shift));
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
			//				if(i>3){Shaver.verbose2=true;}
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
	
	final long rcomp(long kmer){return AminoAcid.reverseComplementBinaryFast(kmer, k);}
	final long toValue(long kmer, long rkmer){return tables.toValue(kmer, rkmer);}
	public final int getCount(long kmer, long rkmer){return tables.getCount(kmer, rkmer);}
	final boolean claim(long kmer, int id){return claim(kmer, rcomp(kmer), id);}
	final boolean claim(long kmer, long rkmer, int id){return tables.claim(kmer, rkmer, id);}
	final boolean doubleClaim(ByteBuilder bb, int id/*, long rid*/){return tables.doubleClaim(bb, id/*, rid*/);}
	final boolean claim(ByteBuilder bb, int id, /*long rid, */boolean earlyExit){return tables.claim(bb, id/*, rid*/, earlyExit);}
	final boolean claim(byte[] array, int len, int id, /*long rid, */boolean earlyExit){return tables.claim(array, len, id/*, rid*/, earlyExit);}
	final int findOwner(long kmer){return tables.findOwner(kmer);}
	final int findOwner(ByteBuilder bb, int id){return tables.findOwner(bb, id);}
	final int findOwner(byte[] array, int len, int id){return tables.findOwner(array, len, id);}
	final void release(long key, int id){tables.release(key, id);}
	final void release(ByteBuilder bb, int id){tables.release(bb, id);}
	final void release(byte[] array, int len, int id){tables.release(array, len, id);}
	final int fillRightCounts(long kmer, long rkmer, int[] counts){return tables.fillRightCounts(kmer, rkmer, counts, mask, shift2);}
	final int fillLeftCounts(long kmer, long rkmer, int[] counts){return tables.fillLeftCounts(kmer, rkmer, counts, mask, shift2);}
	final StringBuilder toText(long kmer){return AbstractKmerTable.toText(kmer, k);}
	
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
		
		private boolean processNextTable(AtomicInteger aint){
			final int tnum=aint.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
			final HashArray1D table=tables.getTable(tnum);
			if(verbose2 && id==0){outstream.println("Processing table "+tnum+", size "+table.size());}
			final int max=table.arrayLength();
			for(int cell=0; cell<max; cell++){
				int x=processCell(table, cell);
			}
			return true;
		}
		
		private boolean processNextVictims(AtomicInteger aint){
			final int tnum=aint.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
			final HashArray1D table=tables.getTable(tnum);
			final HashForest forest=table.victims();
			if(verbose2 && id==0){outstream.println("Processing forest "+tnum+", size "+forest.size());}
			final int max=forest.arrayLength();
			for(int cell=0; cell<max; cell++){
				KmerNode kn=forest.getNode(cell);
				int x=traverseKmerNode(kn);
			}
			return true;
		}
		
		private int processCell(HashArray1D table, int cell){
			int count=table.readCellValue(cell);
			if(count<minCountSeedCurrent){return 0;}
			
			long key=table.getKmer(cell);

			if(verbose){outstream.println("id="+id+" processing cell "+cell+"; \tkmer="+key+"\t"+toText(key));}
			if(useOwnership){
				int owner=table.getCellOwner(cell);
				if(verbose){outstream.println("Owner is initially "+owner);}
				if(owner>-1){return 0;}
				owner=table.setOwner(key, id, cell);
				if(verbose){outstream.println("Owner is now "+owner);}
				if(owner!=id){return 0;}
			}
			return processKmer(key);
		}
		
		private int traverseKmerNode(KmerNode kn){
			int sum=0;
			if(kn!=null){
				sum+=processKmerNode(kn);
				if(kn.left()!=null){
					sum+=traverseKmerNode(kn.left());
				}
				if(kn.right()!=null){
					sum+=traverseKmerNode(kn.right());
				}
			}
			return sum;
		}
		
		private int processKmerNode(KmerNode kn){
			final long key=kn.pivot();
			final int count=kn.getValue(key);
			if(count<minCountSeedCurrent){return 0;}

			if(verbose){outstream.println("id="+id+" processing KmerNode; \tkmer="+key+"\t"+toText(key));}
			if(useOwnership){
				int owner=kn.getOwner(key);
				if(verbose){outstream.println("Owner is initially "+owner);}
				if(owner>-1){return 0;}
				owner=kn.setOwner(key, id);
				if(verbose){outstream.println("Owner is now "+owner);}
				if(owner!=id){return 0;}
			}
			return processKmer(key);
		}
		
		private int processKmer(long key){
			Contig contig=makeContig(key, builderT, true);
			if(contig!=null){
				float coverage=tables.calcCoverage(contig);
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
					x=findInsertSize(r1, r2, rightCounts);
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
					byte[] bases=makeContig(r1.bases, builderT, r1.numericID);
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
					byte[] bases=makeContig(r2.bases, builderT, r2.numericID);
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
		private Contig makeContig(final long key, final ByteBuilder bb, boolean alreadyClaimed){
			builderT.setLength(0);
			builderT.appendKmer(key, k);
			if(verbose){outstream.println("Filled builder: "+builderT);}
			
			final int initialLength=bb.length();
			assert(initialLength==k);
			if(initialLength<k){return null;}
//			outstream.print("A");
			
			boolean success=(alreadyClaimed || !useOwnership ? true : claim(key, id));
			if(verbose){outstream.println("Thread "+id+" checking owner after setting: "+findOwner(bb, id));}
			if(!success){
				assert(bb.length()==k);
//				release(bb, id); //no need to release
				return null;
			}
//			outstream.print("B");
			if(verbose  /*|| true*/){outstream.println("Thread "+id+" building contig; initial length "+bb.length());}
			if(verbose){outstream.println("Extending to right.");}
			final int rightStatus, leftStatus;
			float leftRatio=0, rightRatio=0;
			{
				final int status=extendToRight(bb, leftCounts, rightCounts, id);
				
				if(status==DEAD_END){
					//do nothing
				}else if(status==LOOP){//TODO
					//special case - handle specially, for a loop with no obvious junction, e.g. long tandem repeat.
					//Perhaps, the last kmer should be reclassified as a junction and removed.
				}else if(status==BAD_SEED){
					assert(bb.length()==k);
					release(key, id);
//					outstream.print("B1");
					return null;
				}else{
					if(bb.length()==k){
						if(status==BAD_OWNER){
							if(IGNORE_BAD_OWNER){
								//do nothing
							}else{
								release(key, id);
//								outstream.print("B2");
								return null;
							}
						}else if(isBranchCode(status)){
							release(key, id);
//							outstream.print("B3");
							return null;
						}else{
							throw new RuntimeException("Bad return value: "+status);
						}
					}else{
						if(status==BAD_OWNER){
							if(IGNORE_BAD_OWNER){
								bb.length--;
							}else{
								release(bb, id);
//								outstream.print("B4");
								return null;
							}
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
//			outstream.print("C");
			bb.reverseComplementInPlace();
			if(verbose  /*|| true*/){outstream.println("Extending rcomp to right; current length "+bb.length());}
			{
				final int status=extendToRight(bb, leftCounts, rightCounts, id);
				
				if(status==DEAD_END){
					//do nothing
				}else if(status==LOOP){//TODO
					//special case - handle specially, for a loop with no obvious junction, e.g. long tandem repeat.
					//Perhaps, the last kmer should be reclassified as a junction and removed.
				}else if(status==BAD_SEED){
					assert(false) : bb;//This should never happen.
					assert(bb.length()==k);
					release(key, id);
					return null;
				}else{
					if(status==BAD_OWNER){
						if(IGNORE_BAD_OWNER){
							if(bb.length()>k){bb.length--;}
						}else{
							release(bb, id);
//							outstream.print("C1");
							return null;
						}
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
//			outstream.print("D");

			if(verbose  /*|| true*/){outstream.println("A: Final length for thread "+id+": "+bb.length());}
			
			//				if(useOwnership && THREADS==1){assert(claim(bases, bases.length, id, rid));}
			success=(useOwnership ? doubleClaim(bb, id) : true);
			if(verbose  /*|| true*/){outstream.println("Success for thread "+id+": "+success);}
			
			if(trimEnds>0){bb.trimByAmount(trimEnds, trimEnds);}
			else if(trimCircular && leftStatus==LOOP && rightStatus==LOOP){bb.trimByAmount(0, k-1);}
			if(joinContigs || (bb.length()>=initialLength+minExtension && (bb.length()>=minContigLen || popBubbles))){
				if(success){
					bb.reverseComplementInPlace();
					byte[] bases=bb.toBytes();
					Contig c=new Contig(bases);
					c.leftCode=leftStatus;
					c.rightCode=rightStatus;
					c.rightRatio=rightRatio;
					c.leftRatio=leftRatio;
					if(!c.canonical()){c.rcomp();}
					return c;
				}else{
//					outstream.print("E");
					//					assert(false) : bb.length()+", "+id;
					release(bb, id);
					return null;
				}
			}
			if(verbose  /*|| true*/){outstream.println("A: Contig was too short for "+id+": "+bb.length());}
//			assert(false) : bb.length()+", "+initialLength+", "+minExtension+", "+minContigLen;
//			outstream.print("F");
			return null;
		}
		
		/** From a seed */
		private byte[] makeContig(final byte[] bases, final ByteBuilder bb, long rid){
			if(bases==null || bases.length<k){return null;}
//			if(verbose  /*|| true*/){outstream.println("Thread "+id+" checking owner: "+findOwner(bases, bases.length, id));}
			int owner=useOwnership ? findOwner(bases, bases.length, id) : -1;
			if(owner>=id){return null;}
			boolean success=(useOwnership ? claim(bases, bases.length, id, true) : true);
			if(verbose  /*|| true*/){outstream.println("Thread "+id+" checking owner after setting: "+findOwner(bases, bases.length, id));}
			if(!success){
				release(bases, bases.length, id);
				return null;
			}
			if(verbose  /*|| true*/){outstream.println("Thread "+id+" building contig; initial length "+bases.length);}
			bb.setLength(0);
			bb.append(bases);
			if(verbose){outstream.println("Extending to right.");}
			{
				final int status=extendToRight(bb, leftCounts, rightCounts, id);
				
				if(status==DEAD_END){
					//do nothing
				}else if(status==LOOP){//TODO
					//special case - handle specially, for a loop with no obvious junction, e.g. long tandem repeat.
					//Perhaps, the last kmer should be reclassified as a junction and removed.
				}else if(status==BAD_SEED){
					//do nothing
				}else{
					if(status==BAD_OWNER){
						release(bb, id);
						return null;
					}else if(isBranchCode(status)){
						//do nothing
					}else{
						throw new RuntimeException("Bad return value: "+status);
					}
				}
			}
			bb.reverseComplementInPlace();
			if(verbose  /*|| true*/){outstream.println("Extending rcomp to right; current length "+bb.length());}
			{
				final int status=extendToRight(bb, leftCounts, rightCounts, id);
				
				if(status==DEAD_END){
					//do nothing
				}else if(status==LOOP){//TODO
					//special case - handle specially, for a loop with no obvious junction, e.g. long tandem repeat.
					//Perhaps, the last kmer should be reclassified as a junction and removed.
				}else if(status==BAD_SEED){
					//do nothing
				}else{
					if(status==BAD_OWNER){
						release(bb, id);
						return null;
					}else if(isBranchCode(status)){
						//do nothing
					}else{
						throw new RuntimeException("Bad return value: "+status);
					}
				}
			}
			if(verbose  /*|| true*/){outstream.println("B: Final length for thread "+id+": "+bb.length());}
			
			//				if(useOwnership && THREADS==1){assert(claim(bases, bases.length, id, rid));}
			success=(useOwnership ? doubleClaim(bb, id) : true);
			if(verbose  /*|| true*/){outstream.println("Success for thread "+id+": "+success);}
			if(bb.length()>=bases.length+minExtension && (bb.length()>=minContigLen || popBubbles)){
				if(success){
					bb.reverseComplementInPlace();
					return bb.toBytes();
				}else{
					//					assert(false) : bb.length()+", "+id;
					release(bb.array, bb.length(), id);
					return null;
				}
			}
			
			if(verbose  /*|| true*/){outstream.println("B: Contig was too short for "+id+": "+bb.length());}
			return null;
		}
		
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
		{
			int cnum=0;
			for(Contig c : contigs){
				c.id=cnum;
//				if(c.leftBranch()){
				if(true){
					long kmer=c.leftKmer(k);
					tables.claim(kmer, rcomp(kmer), cnum);
				}
//				if(c.rightBranch()){
				if(true){
					long kmer=c.rightKmer(k);
					tables.claim(kmer, rcomp(kmer), cnum);
				}
				cnum++;
			}
		}
	}
	
	class ProcessContigThread extends AbstractProcessContigThread {
		
		ProcessContigThread(ArrayList<Contig> contigs_, AtomicInteger next_){
			super(contigs_, next_);
			lastExitCondition=BAD_SEED;
		}
		
		@Override
		public void processContigLeft(Contig c, int[] leftCounts, int[] rightCounts, int[] extraCounts, ByteBuilder bb){
			if(c.leftCode==DEAD_END){return;}
			
			final long kmer0=c.leftKmer(k);
			final long rkmer0=rcomp(kmer0);
			
			assert(tables.getCount(kmer0, rkmer0)>0);
			assert(tables.findOwner(kmer0)==c.id) : tables.findOwner(kmer0)+", "+c.id;
//			System.err.println(tables.findOwner(kmer0));
			
			int leftMaxPos=fillLeftCounts(kmer0, rkmer0, leftCounts);
			int leftMax=leftCounts[leftMaxPos];
			int leftSecondPos=Tools.secondHighestPosition(leftCounts);
			int leftSecond=leftCounts[leftSecondPos];
			
			for(int x=0; x<leftCounts.length; x++){
				bb.clear();
				final int count=leftCounts[x];
				int target=-1;
				if(count>0 && isJunction(leftMax, count)){
					long x2=3-x;
					long rkmer=((rkmer0<<2)|(long)x2)&mask;
					long kmer=(kmer0>>>2)|(((long)x)<<shift2);
					assert(kmer==rcomp(rkmer));
					assert(tables.getCount(kmer, rkmer)==count) : count+", "+tables.getCount(kmer, rkmer);
					bb.append(AminoAcid.numberToBase[x]);
					target=exploreRight(rkmer, kmer, extraCounts, rightCounts, bb);
					if(verbose){
						outstream.println(c.id+"L_F: x="+x+", cnt="+count+", dest="+target
								+", "+codeStrings[lastExitCondition]+", len="+lastLength+", orient="+lastOrientation);
					}
				}
				if(target>=0){
					Edge se=new Edge(c.id, target, lastLength, lastOrientation, count, bb.toBytes());
					c.addLeftEdge(se);
					edgesMadeT++;
				}
			}
		}

		@Override
		public void processContigRight(Contig c, int[] leftCounts, int[] rightCounts, int[] extraCounts, ByteBuilder bb){
			if(c.rightCode==DEAD_END){return;}

			final long kmer0=c.rightKmer(k);
			final long rkmer0=rcomp(kmer0);

			int rightMaxPos=fillRightCounts(kmer0, rkmer0, rightCounts);
			int rightMax=rightCounts[rightMaxPos];
			int rightSecondPos=Tools.secondHighestPosition(rightCounts);
			int rightSecond=rightCounts[rightSecondPos];

			for(int x=0; x<rightCounts.length; x++){
				bb.clear();
				final int count=rightCounts[x];
				int target=-1;
				if(count>0 && isJunction(rightMax, count)){
					long x2=3-x;
					long kmer=((kmer0<<2)|(long)x)&mask;
					long rkmer=(rkmer0>>>2)|(((long)x2)<<shift2);
					assert(kmer==rcomp(rkmer));
					assert(tables.getCount(kmer, rkmer)==count) : count+", "+tables.getCount(kmer, rkmer);
					bb.append(AminoAcid.numberToBase[x]);
					target=exploreRight(kmer, rkmer, leftCounts, extraCounts, bb);
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

		private int exploreRight(long kmer, long rkmer, int[] leftCounts, int[] rightCounts, ByteBuilder bb){
			int length=1;
			int owner=-1;
			lastTarget=-1;
			for(; length<500; length++){
				owner=tables.findOwner(kmer, rkmer);
				if(owner>=0){break;}

				final int leftMaxPos=fillLeftCounts(kmer, rkmer, leftCounts);
				final int leftMax=leftCounts[leftMaxPos];
				final int leftSecondPos=Tools.secondHighestPosition(leftCounts);
				final int leftSecond=leftCounts[leftSecondPos];
				if(isJunction(leftMax, leftSecond)){
					lastExitCondition=B_BRANCH;
					lastLength=length;
					return -1;
				}
				
				final int rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts);
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
				long x2=3-x;
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			}
			lastLength=length;
			lastTarget=owner;
			if(owner>=0){
				lastExitCondition=SUCCESS;
				Contig dest=contigs.get(owner);
				long left=dest.leftKmer(k);
//				if(left==kmer){
//					lastOrientation=0;
//				}else if(left==rkmer){
//					lastOrientation=1;
//				}else{
//					long right=dest.rightKmer(k);
//					if(right==kmer){
//						lastOrientation=2;
//					}else if(right==rkmer){
//						lastOrientation=3;
//					}
//				}
				if(left==kmer || left==rkmer){
					lastOrientation=0;
				}else{
					long right=dest.rightKmer(k);
					if(right==kmer || right==rkmer){
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

	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Extension Methods      ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public int findInsertSize(Read r1, Read r2, int[] rightCounts){
		final long kmer1=tables.rightmostKmer(r1.bases, r1.length());
		final long kmer2=tables.rightmostKmer(r2.bases, r2.length());
		if(kmer1<0 || kmer2<0){return -1;}
		final long rkmer1=rcomp(kmer1);
		final long rkmer2=rcomp(kmer2);
		final int x=measureInsert(kmer1, rkmer1, kmer2, rkmer2, 24000, rightCounts);
		if(x<0){return -1;}
		return r1.length()+r2.length()+x-k;//TODO: May be off by 1.
	}
	
	@Override
	public int extendRead(Read r, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int distance, final Kmer kmer){
		return extendRead(r, bb, leftCounts, rightCounts, distance);
	}

	@Override
	public int extendRead(Read r, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int distance){
		final int initialLen=r.length();
		if(initialLen<k){return 0;}
		bb.setLength(0);
		bb.append(r.bases);
		final int extension=extendToRight2(bb, leftCounts, rightCounts, distance, true);
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
	public int measureInsert(final long kmer1, final long rkmer1, final long kmer2, final long rkmer2, final int maxlen, final int[] rightCounts){
		final long key2=toValue(kmer2, rkmer2);
		long kmer=kmer1;
		long rkmer=rkmer1;
		int len=0;
		
		{
			int count=tables.getCount(key2);
			if(count<minCountSeed){return -1;}
		}
		
		long key=toValue(kmer, rkmer);
		int count=tables.getCount(key);
		if(count<minCountSeed){return -1;}
		if(count<minCountSeed){
			if(verbose){outstream.println("Returning because count was too low: "+count);}
			return -1;
		}
		
		int rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts);
		int rightMax=rightCounts[rightMaxPos];
//		int rightSecondPos=Tools.secondHighestPosition(rightCounts);
//		int rightSecond=rightCounts[rightSecondPos];
		
		if(rightMax<minCountExtend){return -1;}
//		if(isJunction(rightMax, rightSecond)){return -1;}
		
		while(key!=key2 && len<maxlen){
			
			//Generate the new kmer
//			final byte b=AminoAcid.numberToBase[rightMaxPos];
			final long x=rightMaxPos;
			final long x2=AminoAcid.numberToComplement[(int)x];
			
			//Now consider the next kmer
			kmer=((kmer<<2)|(long)x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			
			assert(tables.getCount(kmer, rkmer)==rightMax);
			count=rightMax;
			
			assert(count>=minCountExtend) : count;
			
			rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts);
			rightMax=rightCounts[rightMaxPos];
//			rightSecondPos=Tools.secondHighestPosition(rightCounts);
//			rightSecond=rightCounts[rightSecondPos];
			
			if(verbose){
				outstream.println("kmer: "+toText(kmer)+", "+toText(rkmer));
				outstream.println("Counts: "+count+", "+Arrays.toString(rightCounts));
				outstream.println("rightMaxPos="+rightMaxPos);
				outstream.println("rightMax="+rightMax);
//				outstream.println("rightSecondPos="+rightSecondPos);
//				outstream.println("rightSecond="+rightSecond);
			}
			
			if(rightMax<minCountExtend){
				if(verbose){outstream.println("A: Breaking because highest right was too low:"+rightMax);}
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
	public int extendToRight(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int id){
		if(bb.length()<k){return BAD_SEED;}
		long kmer=0;
		long rkmer=0;
		int len=0;
		
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts, to get the rightmost kmer */
		{
			final int bblen=bb.length();
			final byte[] bases=bb.array;
			for(int i=bblen-k; i<bblen; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{len++;}
				if(verbose){outstream.println("A: Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			}
		}
		
		if(len<k){return BAD_SEED;}
		else{assert(len==k);}
		
		/* Now the trailing kmer has been initialized. */
		
		long key=toValue(kmer, rkmer);
		HashArray1D table=tables.getTableForKey(key);
		int count=table.getValue(key);
		if(count<minCountSeed){
			if(verbose){outstream.println("Returning because count was too low: "+count);}
			return BAD_SEED;
		}
		
		int owner=(useOwnership ? table.getOwner(key) : id);
		if(verbose){outstream.println("Owner: "+owner);}
		if(owner>id){return BAD_OWNER;}
		
		int leftMaxPos=0;
		int leftMax=minCountExtend;
		int leftSecondPos=1;
		int leftSecond=0;
		
		if(leftCounts!=null){
			leftMaxPos=fillLeftCounts(kmer, rkmer, leftCounts);
			leftMax=leftCounts[leftMaxPos];
			leftSecondPos=Tools.secondHighestPosition(leftCounts);
			leftSecond=leftCounts[leftSecondPos];
		}
		
		int rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts);
		int rightMax=rightCounts[rightMaxPos];
		int rightSecondPos=Tools.secondHighestPosition(rightCounts);
		int rightSecond=rightCounts[rightSecondPos];
		
		if(verbose){
			outstream.println("kmer: "+toText(kmer)+", "+toText(rkmer));
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
			assert(bb.length()==k) : bb.length()+", "+k+", "+leftMax+", "+leftSecond;
			if(verbose){outstream.println("B: Breaking because isJunction("+rightMax+", "+rightSecond+", "+leftMax+", "+leftSecond+")");}
			return B_BRANCH;
		}
		
		if(useOwnership){
			owner=table.setOwner(key, id);
			if(verbose){outstream.println("A. Owner is now "+id+" for key "+key);}
			if(owner!=id){
				if(verbose){outstream.println("Returning early because owner was "+owner+" for thread "+id+".");}
				return BAD_OWNER;
			}
		}
		
		final int maxLen=Tools.min((extendRight<0 ? maxContigLen : bb.length()+extendRight), maxContigLen);
		
		while(owner==id && bb.length()<maxLen){
			
			//Generate the new kmer
			final byte b=AminoAcid.numberToBase[rightMaxPos];
			final long x=rightMaxPos;
			final long x2=AminoAcid.numberToComplement[(int)x];
			
			final long evicted=(kmer>>>shift2); //Binary value that falls off the end.
			
			//Now consider the next kmer
			kmer=((kmer<<2)|(long)x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			
			key=toValue(kmer, rkmer);
			table=tables.getTableForKey(key);
			
			assert(table.getValue(key)==rightMax || rightMax==0) : key+", "+table.getValue(key)+", "+rightMax+", "+Arrays.toString(rightCounts);
			count=rightMax;
			
			assert(count>=minCountExtend) : count;

			if(leftCounts!=null){
				leftMaxPos=fillLeftCounts(kmer, rkmer, leftCounts);
				leftMax=leftCounts[leftMaxPos];
				leftSecondPos=Tools.secondHighestPosition(leftCounts);
				leftSecond=leftCounts[leftSecondPos];
			}
			
			rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts);
			rightMax=rightCounts[rightMaxPos];
			rightSecondPos=Tools.secondHighestPosition(rightCounts);
			rightSecond=rightCounts[rightSecondPos];
			
			if(verbose){
				outstream.println("kmer: "+toText(kmer)+", "+toText(rkmer));
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
				owner=table.getOwner(key);
				if(verbose){outstream.println("Owner is initially "+id+" for key "+key);}
				if(owner==id){//loop detection
					if(verbose  /*|| true*/){
//						outstream.println(new String(bb.array, bb.length()-31, 31));
						outstream.println(bb);
						outstream.println(toText(kmer));
						outstream.println(toText(rkmer));
						outstream.println("Breaking because owner was "+owner+" for thread "+id+".");
					}
					return fbranch ? F_BRANCH : LOOP;
				}
				owner=table.setOwner(key, id);
				if(verbose){outstream.println("B. Owner is now "+id+" for key "+key);}
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
		if(verbose || verbose2){outstream.println("Entering extendToRight2 (no kmers).");}
		final int initialLength=bb.length();
		if(initialLength<k){return 0;}
		long kmer=0;
		long rkmer=0;
		
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts, to get the rightmost kmer */
		{
			int len=0;
			final byte[] bases=bb.array;
			for(int i=initialLength-k; i<initialLength; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{len++;}
				if(verbose){outstream.println("B: Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			}
			if(len<k){
				if(verbose || verbose2){outstream.println("Returning because len<k: "+len+"<"+k);}
				return 0;
			}
			else{assert(len==k);}
		}
		return extendToRight2(bb, leftCounts, rightCounts, distance, includeJunctionBase, kmer, rkmer);
	}
	
	/**
	 * Extend these bases to the right by at most 'distance'.
	 * Stops at right junctions only.
	 * Does not claim ownership.
	 */
	public int extendToRight2(final ByteBuilder bb, final int[] leftCounts, final int[] rightCounts, final int distance, boolean includeJunctionBase,
			long kmer, long rkmer){
		if(verbose || verbose2){outstream.println("Entering extendToRight2 (with kmers).");}
		final int initialLength=bb.length();
		
		/* Now the trailing kmer has been initialized. */
		
		long key=toValue(kmer, rkmer);
		HashArray1D table=tables.getTableForKey(key);
		int count=table.getValue(key);
		if(count<minCountSeed){
			if(verbose || verbose2){outstream.println("Returning because count was too low: "+count+"<"+minCountSeed);}
			return 0;
		}
		
		int leftMaxPos=0;
		int leftMax=minCountExtend;
		int leftSecondPos=1;
		int leftSecond=0;
		
		if(leftCounts!=null){
			leftMaxPos=fillLeftCounts(kmer, rkmer, leftCounts);
			leftMax=leftCounts[leftMaxPos];
			leftSecondPos=Tools.secondHighestPosition(leftCounts);
			leftSecond=leftCounts[leftSecondPos];
		}
		
		int rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts);
		int rightMax=rightCounts[rightMaxPos];
		int rightSecondPos=Tools.secondHighestPosition(rightCounts);
		int rightSecond=rightCounts[rightSecondPos];
		
		if(verbose){
			outstream.println("kmer: "+toText(kmer)+", "+toText(rkmer));
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
			final long x=rightMaxPos;
			final long x2=AminoAcid.numberToComplement[(int)x];
			
			final long evicted=(kmer>>>shift2); //Binary value that falls off the end.
			
			//Now consider the next kmer
			kmer=((kmer<<2)|(long)x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			
			key=toValue(kmer, rkmer);
			table=tables.getTableForKey(key);
			
			assert(table.getValue(key)==rightMax || rightMax==0);
			count=rightMax;
			
			assert(count>=minCountExtend) : count;
			
			if(leftCounts!=null){
				leftMaxPos=fillLeftCounts(kmer, rkmer, leftCounts);
				leftMax=leftCounts[leftMaxPos];
				leftSecondPos=Tools.secondHighestPosition(leftCounts);
				leftSecond=leftCounts[leftSecondPos];
			}
			
			rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts);
			rightMax=rightCounts[rightMaxPos];
			rightSecondPos=Tools.secondHighestPosition(rightCounts);
			rightSecond=rightCounts[rightSecondPos];
			
			if(verbose){
				outstream.println("kmer: "+toText(kmer)+", "+toText(rkmer));
				outstream.println("Counts: "+count+", "+Arrays.toString(rightCounts));
				outstream.println("rightMaxPos="+rightMaxPos);
				outstream.println("rightMax="+rightMax);
				outstream.println("rightSecondPos="+rightSecondPos);
				outstream.println("rightSecond="+rightSecond);
			}

			if(isJunction(rightMax, rightSecond, leftMax, leftSecond)){
				if(verbose){outstream.println("B: Breaking because isJunction("+rightMax+", "+rightSecond+", "+leftMax+", "+leftSecond+")");}
				if(includeJunctionBase && kmer>rkmer){
					bb.append(b);
					if(verbose){outstream.println("Added base "+(char)b);}
				}
				break;
			}
			
			if(leftCounts!=null && leftMaxPos!=evicted){
				if(verbose){outstream.println("B: Breaking because of hidden branch: leftMaxPos!=evicted ("+leftMaxPos+"!="+evicted+")" +
						"\nleftMaxPos="+leftMaxPos+", leftMax="+leftMax+", leftSecondPos="+leftSecondPos+", leftSecond="+leftSecond);}
				if(includeJunctionBase && kmer>rkmer){
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
		boolean junk=isJunk(r, localRightCounts.get());
		return junk;
	}
	
	@Override
	public boolean isJunk(Read r, final int[] counts, Kmer kmer){
		return isJunk(r, counts);
	}
	
	public boolean isJunk(Read r, final int[] counts){
		final int blen=r.length();
		if(blen<k){return true;}
		final byte[] bases=r.bases;
		
		long kmer=0;
		long rkmer=0;
		int len=0;

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts, to get the leftmost kmer */
		{
			for(int i=0; i<k; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{len++;}
				if(verbose){outstream.println("C: Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			}
		}
		
		if(len>=k){
			int maxPos=fillLeftCounts(kmer, rkmer, counts);
			if(counts[maxPos]>0){
				//outstream.println("Not junk, left counts[maxPos]="+counts[maxPos]+"; right="+fillRightCounts(kmer, rkmer, counts, mask, shift2));
				return false;
			}
		}
		
		final boolean paired=(r.mateLength()>=k);
		int maxDepth=0;
		{
			for(int i=k; i<blen; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{
					len++;
					if(len>=k){
						int depth=getCount(kmer, rkmer);
						if(depth>maxDepth){
							maxDepth=depth;
							if(maxDepth>1 && (!paired || maxDepth>2)){
								//outstream.println("Not junk, maxDepth="+maxDepth);
								return false;
							}
						}
					}
				}
				if(verbose){outstream.println("D: Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
			}
		}
		
		if(len>=k && !paired){
			int maxPos=fillRightCounts(kmer, rkmer, counts);
			if(counts[maxPos]>0){
				//outstream.println("Not junk, right counts[maxPos]="+counts[maxPos]);
				return false;
			}
		}
		return true;
	}
	
	@Override
	public boolean hasKmersAtOrBelow(Read r, int tooLow, final float fraction, Kmer kmer){
		return hasKmersAtOrBelow(r, tooLow, fraction);
	}
	
	@Override
	public boolean hasKmersAtOrBelow(Read r, final int tooLow, final float fraction){
		final int blen=r.length();
		if(blen<k){return true;}
		final byte[] bases=r.bases;
		
//		outstream.println("\n"+new String(r.bases)+":");
		
		long kmer=0;
		long rkmer=0;
		int len=0;
		
		final int limit=Tools.max(1, Math.round((bases.length-kbig+1)*fraction));
		int valid=0, invalid=0;

		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts, to get the leftmost kmer */
		{
			for(int i=0; i<blen; i++){
				final byte b=bases[i];
				final long x=AminoAcid.baseToNumber[b];
				final long x2=AminoAcid.baseToComplementNumber[b];
				kmer=((kmer<<2)|x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				if(x<0){
					len=0;
					kmer=rkmer=0;
				}else{
					len++;
					if(len>=k){
						int depth=getCount(kmer, rkmer);
//						outstream.println("depth="+depth+", kmer="+toText(kmer));
						if(depth>tooLow){valid++;}
						else{
							invalid++;
							if(invalid>=limit){return true;}
						}
					}
				}
				if(verbose){outstream.println("E: Scanning i="+i+", len="+len+", kmer="+kmer+", rkmer="+rkmer+"\t"+new String(bases, Tools.max(0, i-k2), Tools.min(i+1, k)));}
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
		int corrected=errorCorrect(r, localLeftCounts.get(), localRightCounts.get(), localLongList.get(),
				localIntList.get(), localIntList2.get(), localByteBuilder.get(), localByteBuilder2.get(), localTracker.get(), localBitSet.get());
		return corrected;
	}
	
	@Override
	public int errorCorrect(Read r, final int[] leftCounts, final int[] rightCounts, LongList kmers, IntList counts, IntList counts2,
			final ByteBuilder bb, final ByteBuilder bb2, final ErrorTracker tracker, final BitSet bs, Kmer kmer, Kmer kmer2){
		return errorCorrect(r, leftCounts, rightCounts, kmers, counts, counts2, bb, bb2, tracker, bs);
	}
	
	boolean hasErrorsFast(LongList kmers){
		if(kmers.size<1){return false;}
		int prev=-1;
		
		final int incr=Tools.mid(1, k/2, 9), mcc=minCountCorrect();
		for(int i=0; i<kmers.size; i+=incr){
			long kmer=kmers.get(i);
			if(kmer<0){
				return true;
			}
			long rkmer=rcomp(kmer);
			int count=getCount(kmer, rkmer);
			final int min=Tools.min(count, prev), max=Tools.max(count, prev);
			if(count<mcc || (i>0 && (isError(max+1, min-1)))){return true;}
			prev=count;
		}
		
		long kmer=kmers.get(kmers.size()-1);
		if(kmer<0){return true;}
		long rkmer=rcomp(kmer);
		int count=getCount(kmer, rkmer);
		final int min=Tools.min(count, prev), max=Tools.max(count, prev);
		return count<mcc || isError(max+1, min-1);
	}
	
	public int errorCorrect(Read r, final int[] leftCounts, final int[] rightCounts, LongList kmers, IntList counts, IntList counts2,
			final ByteBuilder bb, final ByteBuilder bb2, final ErrorTracker tracker, final BitSet bs){
		
		final byte[] bases=r.bases;
		final byte[] quals=r.quality;
		tracker.clear();
		int valid=tables.fillKmers(bases, kmers);
		if(valid<2){return 0;}
		if(!r.containsUndefined() && !hasErrorsFast(kmers)){return 0;}
		
		tables.fillCounts(kmers, counts);
		final int possibleErrors=tracker.suspected=countErrors(counts, quals);
		if(possibleErrors<0){return 0;}
		final float expectedErrors=r.expectedErrors(true, r.length());
		final Rollback roll=ECC_ROLLBACK ? new Rollback(r, counts) : null;
		
		assert(counts.size>0);
		
		int correctedPincer=0;
		int correctedTail=0;
		int correctedBrute=0;
		int correctedReassemble=0;
		
		if(ECC_PINCER){
			correctedPincer+=errorCorrectPincer(bases, quals, leftCounts, rightCounts, kmers, counts, bb, tracker, errorExtensionPincer);
		}
		
		if(ECC_TAIL || ECC_ALL){
			int start=(ECC_ALL ? 0 : counts.size-k-1);
//			if(ECC_PINCER && tracker!=null && tracker.detected>correctedPincer){start=start-k;}
			correctedTail+=errorCorrectTail(bases, quals, leftCounts, rightCounts, kmers, counts, bb, tracker, start, errorExtensionTail);
			r.reverseComplement();
			valid=tables.fillKmers(bases, kmers);
			counts.reverse();
			correctedTail+=errorCorrectTail(bases, quals, leftCounts, rightCounts, kmers, counts, bb, tracker, start, errorExtensionTail);
			r.reverseComplement();
			counts.reverse();
		}
		
		if(ECC_REASSEMBLE){
			if(verbose){outstream.println("Correcting "+possibleErrors+" errors.  Counts:\n"+counts);}
			if((correctedPincer<1 && correctedTail<1) || countErrors(counts, quals)>0){
				correctedReassemble=reassemble(bases, quals, rightCounts, counts, counts2, tracker, errorExtensionReassemble, bb, bb2, null, null, bs);
			}
//			assert(tracker.detectedReassemble>0) : counts;
			if(verbose){outstream.println("Corrected  "+correctedReassemble+" errors.  Counts:\n"+counts);}
		}
		assert(counts.size>0);
		
		//123 For testing.
		if(false && tracker.detected()>tracker.corrected()){
			correctedBrute+=errorCorrectBruteForce(bases, quals, leftCounts, rightCounts, kmers, counts, bb, tracker, errorExtensionPincer);
		}
		
		assert(correctedPincer+correctedTail+correctedReassemble+correctedBrute==tracker.corrected())
			: correctedPincer+", "+correctedTail+", "+correctedReassemble+", "+correctedBrute+", "+tracker;

		if(ECC_ROLLBACK && (tracker.corrected()>0 || tracker.rollback)){
			
			if(!tracker.rollback && quals!=null && tracker.corrected()>3){
				float mult=Tools.max(1, 0.5f*(0.5f+0.01f*r.length()));//1 for a 150bp read.
				if(countErrors(counts, quals)>0 && tracker.corrected()>mult+expectedErrors){tracker.rollback=true;}
				else if(tracker.corrected()>2.5f*mult+expectedErrors){tracker.rollback=true;}
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
				}
//				else if(b<a-1 && !printed){
//					assert(false);
//					if(verbose){outstream.println("X: RID="+r.numericID+"; "+a+"->"+b+"\n"+counts0+"\n"+counts);}
//					printed=true;
//				}
			}
			
			if(tracker.rollback){
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
	
	private int findBestMutant(final byte[] bases, final int a, final LongList kmers, final long[] bestMutant){
		Arrays.fill(bestMutant, -1);
		final long kmer0=kmers.get(a);
		if(kmer0<0){return -1;}
		long mask=3;
		int max=0, second=0;
		for(int i=0; i<k; i++){
			for(long num=0; num<4; num++){
				final long kmer=(kmer0&~mask)|(num<<(2*i));
				if(kmer!=kmer0){
					final long rkmer=rcomp(kmer);
					int count=getCount(kmer, rkmer);
					if(count>max){
						second=max;
						max=count;
						bestMutant[0]=kmer;
						bestMutant[1]=rkmer;
						bestMutant[2]=i;
					}else if(count>second){
						second=count;
					}
				}
			}
			mask<<=2;
		}
		bestMutant[3]=max;
		bestMutant[4]=second;
		return max>second ? max : -1;
	}
	
	private boolean fixBestMutant(final byte[] bases, final int a, final LongList kmers, final IntList counts, final int pos, final long kmer){
		int shift=2*pos;
		long mask=(3L)<<shift;
		byte num=(byte)((kmer>>shift)&mask);
		byte base=AminoAcid.numberToBase[num];
		byte old=bases[a+pos];
		bases[a+pos]=base;
		
//		long[] copy=counts.toArray();
		
//		tables.regenerateKmers(bases, kmers, counts, a);
		tables.fillKmers(bases, kmers);
		tables.fillCounts(kmers, counts);
		
//		for(int i=0; i<copy.length; i++){
//			if(counts.get(i)<copy[i]){
//				bases[a+pos]=old;
//				tables.fillKmers(bases, kmers);
//				tables.fillCounts(kmers, counts);
//				return false;
//			}
//		}
		
		return true;
	}
	
	public int errorCorrectBruteForce(final byte[] bases, final byte[] quals, final int[] leftBuffer, final int[] rightBuffer,
			final LongList kmers, final IntList counts, final ByteBuilder bb, final ErrorTracker tracker, final int errorExtension){
		
		int detected=0;
		int corrected=0;
		
		final long[] bestMutant=new long[5];
		final int[] countCopy=counts.toArray();
		final byte[] basesCopy=bases.clone();
		
		final int maxCount=Tools.max(counts.array, 0, counts.size-1);
		if(maxCount<minCountCorrect()){
			return 0;
		}
		int thresh=Tools.max(3, maxCount/4);
		
		for(int a=0, d=k+1; a<counts.size; a++, d++){
			final int aCount=counts.get(a);
			if(aCount<thresh){
				detected++;
				boolean fixed=false;
				int mCount=findBestMutant(bases, a, kmers, bestMutant);
				if(mCount>thresh && mCount>=minCountCorrect() && (isSimilar(mCount, maxCount) || (mCount<maxCount && mCount>(1+aCount)*3))){
					long kmer=bestMutant[0];
					long rkmer=bestMutant[1];
					int pos=(int)bestMutant[2];
					assert(mCount==bestMutant[3]);
					int secondBestCount=(int)bestMutant[4];
					if(mCount>(1+secondBestCount*3) && pos>=errorExtension && pos<=(k-errorExtension)){
						fixBestMutant(bases, a, kmers, counts, pos, kmer);
						fixed=true;
						corrected++;
					}
				}
			}
		}

		for(int i=0; i<countCopy.length; i++){
			if(counts.get(i)<countCopy[i]){
				for(int j=0; j<bases.length; j++){bases[j]=basesCopy[j];}
				counts.clear();
				for(int count : countCopy){counts.add(count);}
				tables.fillKmers(bases, kmers);
				corrected=0;
			}
		}
		
		{
			tracker.detectedBrute+=detected;
			tracker.correctedBrute+=corrected;
		}
		
		return corrected;
	}

	public int errorCorrectPincer(final byte[] bases, final byte[] quals, final int[] leftBuffer, final int[] rightBuffer,
			final LongList kmers, final IntList counts, final ByteBuilder bb, final ErrorTracker tracker, final int errorExtension){
		
		int detected=0;
		int corrected=0;
		
		//a is the index of the left kmer
		//b is a+1 (right-extension of left kmer)
		//c is d-1 (left-extension of right kmer)
		//d is the index of the right kmer
		//the base between the kmers is at a+k
		for(int a=0, d=k+1; d<counts.size; a++, d++){
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
				int ret=correctSingleBasePincer(a, d, bases, quals, leftBuffer, rightBuffer, kmers, counts, bb, errorExtension);
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

	public int errorCorrectTail(final byte[] bases, final byte[] quals, final int[] leftBuffer, final int[] rightBuffer,
			final LongList kmers, final IntList counts, final ByteBuilder bb, final ErrorTracker tracker, final int startPos, final int errorExtension){
		if(bases.length<k+2+errorExtension+deadZone){return 0;}
		int detected=0;
		int corrected=0;
		
		//a is the index of the left kmer
		//b is a+1
		//the base between the kmers is at a+k
		for(int a=Tools.max(startPos, errorExtension), lim=counts.size-deadZone-1; a<lim; a++){//errorExtension-1
			final int aCount=counts.get(a);
			final int bCount=counts.get(a+1);
			final byte qb=(quals==null ? 20 : quals[a+kbig]);
			if(isError(aCount, bCount, qb) && isSimilar(aCount, a-errorExtension, a-1, counts) && isError(aCount, a+2, a+k, counts)){
				if(verbose){
					outstream.println("Found error: "+aCount+", "+bCount);
				}
				//Assume like a 1bp substitution; attempt to correct.
				detected++;
				int ret=correctSingleBaseRight(a, bases, quals, leftBuffer, rightBuffer, kmers, counts, bb, errorExtension);
				corrected+=ret;
				if(verbose){
					outstream.println("Corrected error.");
				}
			}else{
				if(verbose){
					outstream.println("Not an error: "+aCount+", "+bCount+
							";  "+isError(aCount, bCount, qb)+", "+isSimilar(aCount, a-errorExtension, a-1, counts)+", "+isError(aCount, a+2, a+k, counts));
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
	
//	public int reassemble_inner(final ByteBuilder bases, final byte[] quals, final int[] rightCounts, final IntList counts,
//			final int errorExtension, final Kmer kmer, final Kmer kmer2){
//		throw new RuntimeException("TODO");
//	}
	
	@Override
	public int reassemble_inner(final ByteBuilder bb, final byte[] quals, final int[] rightCounts, final IntList counts,
			final int errorExtension, final Kmer kmer, final Kmer regenKmer){
		return reassemble_inner(bb, quals, rightCounts, counts, errorExtension);
	}
	
	public int reassemble_inner(final ByteBuilder bb, final byte[] quals, final int[] rightCounts, final IntList counts,
			final int errorExtension){
		final int length=bb.length();
		if(length<k+1+deadZone){return 0;}
		final byte[] bases=bb.array;
		
		long kmer=0, rkmer=0;

		int detected=0;
		int corrected=0;
		int len=0;
		
		//a is the index of the first base of the left kmer
		//b=a+1 is the index of the next base
		//ca=a-k is the index of the next base
		//cb=a-k is the index of the next base
		//the base between the kmers is at a+k
		for(int a=0, lim=length-deadZone-1; a<lim; a++){
			
			//Generate the new kmer
			final byte aBase=bases[a];
			final long x=AminoAcid.baseToNumber[aBase];
			final long x2=AminoAcid.baseToComplementNumber[aBase];
			
			if(x<0){
				len=0;
				kmer=rkmer=0;
			}else{


				//Now consider the next kmer
				kmer=((kmer<<2)|(long)x)&mask;
				rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
				len++;

				if(verbose){
					outstream.println("len: "+len+" vs "+k+"; a="+a);
				}

				if(len>=k){

					final int b=a+1;
					final int ca=a-k+1;
					final int cb=ca+1;

					final int aCount=counts.get(ca);
					final int bCount=counts.get(cb);
					final byte qb=(quals==null ? 20 : quals[b]);

					if(verbose){
						outstream.println("ca="+ca+", cb="+cb+"; aCount="+aCount+", bCount="+bCount);
						outstream.println(isError(aCount, bCount, qb)+", "+isSimilar(aCount, ca-errorExtension, ca-1, counts)+
								", "+isError(aCount, ca+2, ca+k, counts));
					}

//					if(isError(aCount, bCount) && isSimilar(aCount, ca-errorExtension, ca-1, counts) && isError(aCount, ca+2, ca+k, counts)){
					if(isSubstitution(ca, errorExtension, qb, counts)){
						if(verbose){
							outstream.println("***Found error: "+aCount+", "+bCount);
						}
						//Assume a 1bp substitution; attempt to correct.

						int rightMaxPos=fillRightCounts(kmer, rkmer, rightCounts);
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
								tables.regenerateCounts(bases, counts, ca);
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
									";  "+isError(aCount, bCount, qb)+", "+isSimilar(aCount, a-errorExtension, a-1, counts)+", "+isError(aCount, a+2, a+k, counts));
						}
					}
				}
			}
		}
		
		return corrected;
	}
	
	private int correctSingleBasePincer(final int a, final int d, final byte[] bases, final byte[] quals, final int[] leftBuffer, final int[] rightBuffer,
			final LongList kmers, final IntList counts, final ByteBuilder bb, final int errorExtension){
		final byte leftReplacement, rightReplacement;
		final int loc=a+k;
		{
			bb.clear();
			final long kmer=kmers.get(a);
			final long rkmer=rcomp(kmer);
			int extension=extendToRight2(bb, null, rightBuffer, errorExtension, true, kmer, rkmer);
			if(extension<errorExtension){return 0;}
			for(int i=1; i<extension; i++){
				if(bb.get(i)!=bases[loc+i]){
					return 0;
				}
			}
			leftReplacement=bb.get(0);
		}
		{
			bb.clear();
			final long rkmer=kmers.get(d);
			final long kmer=rcomp(rkmer);
			int extension=extendToRight2(bb, null, rightBuffer, errorExtension, true, kmer, rkmer);
			if(extension<errorExtension){return 0;}
			bb.reverseComplementInPlace();
			for(int i=0; i<extension-1; i++){
				if(bb.get(i)!=bases[loc+i+1-extension]){
					return 0;
				}
			}
			rightReplacement=bb.get(extension-1);
		}
		if(leftReplacement!=rightReplacement){return 0;}
		if(bases[loc]==leftReplacement){return 0;}
		if(!isSimilar(a, leftReplacement, kmers, counts)){return 0;}
		
		bases[loc]=leftReplacement;
		assert(d==a+k+1);
		tables.regenerateKmers(bases, kmers, counts, a);
		return 1;
	}
	
	private int correctSingleBaseRight(final int a, final byte[] bases, final byte[] quals, final int[] leftBuffer, final int[] rightBuffer,
			final LongList kmers, final IntList counts, final ByteBuilder bb, final int errorExtension0){
		final byte leftReplacement;
		final int loc=a+k;
		final int errorExtension=Tools.min(errorExtension0, bases.length-loc);
		{
			bb.clear();
			final long kmer=kmers.get(a);
			final long rkmer=rcomp(kmer);
			int extension=extendToRight2(bb, null, rightBuffer, errorExtension, true, kmer, rkmer);
			if(extension<errorExtension){return 0;}
			for(int i=1; i<extension; i++){
				if(bb.get(i)!=bases[loc+i]){
					return 0;
				}
			}
			leftReplacement=bb.get(0);
		}
		
		if(bases[loc]==leftReplacement){return 0;}
		if(!isSimilar(a, leftReplacement, kmers, counts)){return 0;}
		
		bases[loc]=leftReplacement;
		tables.regenerateKmers(bases, kmers, counts, a);
		return 1;
	}
	
	private boolean isSimilar(int a, byte newBase, LongList kmers, IntList counts){
		long kmer=kmers.get(a);
				
		final long x=AminoAcid.baseToNumber[newBase];
		kmer=((kmer<<2)|x)&mask;
		long rkmer=rcomp(kmer);
		int count=getCount(kmer, rkmer);
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
	public final KmerTableSet tables(){return tables;}
	public final KmerTableSet tables;
	
	/** Normal kmer length */
	final int k;
	/** k-1; used in some expressions */
	final int k2;
	
	final int shift;
	final int shift2;
	final long mask;
	
}
