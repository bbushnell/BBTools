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
	
	public final int getCount(Kmer kmer){return tables.getCount(kmer);}
	final boolean claim(Kmer kmer, int id){return tables.claim(kmer, id);}
	final boolean doubleClaim(ByteBuilder bb, int id/*, long rid*/, Kmer kmer){return tables.doubleClaim(bb, id/*, rid*/, kmer);}
	final boolean claim(ByteBuilder bb, int id, /*long rid, */boolean earlyExit, Kmer kmer){return tables.claim(bb, id/*, rid*/, earlyExit, kmer);}
	final boolean claim(byte[] array, int len, int id, /*long rid, */boolean earlyExit, Kmer kmer){return tables.claim(array, len, id/*, rid*/, earlyExit, kmer);}
	final int findOwner(Kmer kmer){return tables.findOwner(kmer);}
	final int findOwner(ByteBuilder bb, int id, Kmer kmer){return tables.findOwner(bb, id, kmer);}
	final int findOwner(byte[] array, int len, int id, Kmer kmer){return tables.findOwner(array, len, id, kmer);}
	final void release(Kmer kmer, int id){tables.release(kmer, id);}
	final void release(ByteBuilder bb, int id, Kmer kmer){tables.release(bb, id, kmer);}
	final void release(byte[] array, int len, int id, Kmer kmer){tables.release(array, len, id, kmer);}
	final int fillRightCounts(Kmer kmer, int[] counts){return tables.fillRightCounts(kmer, counts);}
	final int fillLeftCounts(Kmer kmer, int[] counts){return tables.fillLeftCounts(kmer, counts);}
	final static StringBuilder toText(Kmer kmer){return AbstractKmerTableU.toText(kmer);}
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
			final HashArrayU1D table=tables.getTable(tnum);
			final int max=table.arrayLength();
			if(verbose && id==0){outstream.println("Processing table "+tnum+", size "+table.size()+", length "+max);}
			for(int cell=0; cell<max; cell++){
				if(verbose && id==0){outstream.println("Processing cell "+cell);}
				int x=processCell(table, cell, myKmer);
			}
			return true;
		}
		
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
		
		private int processKmerNodeU(KmerNodeU kn){
			final long[] key=kn.pivot();
			final int count=kn.getValue(key);
			if(count<minCountSeedCurrent){return 0;}

			if(verbose){outstream.println("id="+id+" processing KmerNodeU; \tkmer="+Arrays.toString(key)+"\t"+toText(key, ksmall));}
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
		
		private final Kmer myKmer=new Kmer(kbig);
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
		{
			int cnum=0;
			for(Contig c : contigs){
				c.id=cnum;
//				if(c.leftBranch()){
				if(true){
					c.leftKmer(kmer);
					tables.claim(kmer, cnum);
				}
//				if(c.rightBranch()){
				if(true){
					c.rightKmer(kmer);
					tables.claim(kmer, cnum);
				}
				cnum++;
			}
		}
	}

	class ProcessContigThread extends AbstractProcessContigThread {

		ProcessContigThread(ArrayList<Contig> contigs_, AtomicInteger next_){
			super(contigs_, next_);
			kmerA=new Kmer(kbig);
			kmerB=new Kmer(kbig);
			kmerC=new Kmer(kbig);
			lastExitCondition=BAD_SEED;
		}

		@Override
		public void processContigLeft(Contig c, int[] leftCounts, int[] rightCounts, int[] extraCounts, ByteBuilder bb){
//			System.err.println("processContigLeft: "+c.id+", "+c.leftCode+", "+c.leftEdges);
			if(c.leftCode==DEAD_END){return;}
			
			final Kmer kmer0=c.leftKmer(kmerA);
			final Kmer kmer=kmerB;
			
			assert(tables.getCount(kmer0)>0);
			assert(tables.findOwner(kmer0)==c.id) : tables.findOwner(kmer0)+", "+c.id;

			int leftMaxPos=fillLeftCounts(kmer0, leftCounts);
			int leftMax=leftCounts[leftMaxPos];
			int leftSecondPos=Tools.secondHighestPosition(leftCounts);
			int leftSecond=leftCounts[leftSecondPos];
			
			for(int x=0; x<leftCounts.length; x++){
				bb.clear();
				final int count=leftCounts[x];
				int target=-1;
				if(count>0 && isJunction(leftMax, count)){
					kmer.setFrom(kmer0);
					kmer.addLeftNumeric(x);
					assert(tables.getCount(kmer)==count) : count+", "+tables.getCount(kmer);
					bb.append(AminoAcid.numberToBase[x]);
					target=exploreRight(kmer, extraCounts, rightCounts, bb);
					if(verbose){
						outstream.println(c.id+"L_F: x="+x+", cnt="+count+", dest="+target
								+", "+codeStrings[lastExitCondition]+", len="+lastLength+", orient="+lastOrientation);
					}
				}
				if(target>=0){
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
			if(c.rightCode==DEAD_END){return;}

			final Kmer kmer0=c.rightKmer(kmerA);
			final Kmer kmer=kmerB;

			int rightMaxPos=fillRightCounts(kmer0, rightCounts);
			int rightMax=rightCounts[rightMaxPos];
			int rightSecondPos=Tools.secondHighestPosition(rightCounts);
			int rightSecond=rightCounts[rightSecondPos];

			for(int x=0; x<rightCounts.length; x++){
				bb.clear();
				final int count=rightCounts[x];
				int target=-1;
				if(count>0 && isJunction(rightMax, count)){
					kmer.setFrom(kmer0);
					kmer.addRightNumeric(x);
					assert(tables.getCount(kmer)==count) : count+", "+tables.getCount(kmer);
					bb.append(AminoAcid.numberToBase[x]);
					target=exploreRight(kmer, leftCounts, extraCounts, bb);
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

		private int exploreRight(Kmer kmer, int[] leftCounts, int[] rightCounts, ByteBuilder bb){
			final Kmer temp=kmerC;
			int length=1;
			int owner=-1;
			lastTarget=-1;
			for(; length<500; length++){
				owner=tables.findOwner(kmer);
				if(owner>=0){break;}

				final int leftMaxPos=fillLeftCounts(kmer, leftCounts);
				final int leftMax=leftCounts[leftMaxPos];
				final int leftSecondPos=Tools.secondHighestPosition(leftCounts);
				final int leftSecond=leftCounts[leftSecondPos];
				if(isJunction(leftMax, leftSecond)){
					lastExitCondition=B_BRANCH;
					lastLength=length;
					return -1;
				}
				
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

		final Kmer kmerA, kmerB, kmerC;

	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------       Extension Methods      ----------------*/
	/*--------------------------------------------------------------*/


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
				if(owner==id){//loop detection
					if(verbose  /*|| true*/){
//						outstream.println(new String(bb.array, bb.length()-31, 31));
						outstream.println(bb);
						outstream.println(toText(kmer));
						outstream.println("Breaking because owner was "+owner+" for thread "+id+".");
					}
					return LOOP;
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
	public final KmerTableSetU tables;
	
	/** Normal kmer length */
	final int ksmall;
	
}
