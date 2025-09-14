package clump;

import java.util.ArrayList;

import dna.AminoAcid;
import hiseq.FlowcellCoordinate;
import shared.KillSwitch;
import shared.Parse;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

/**
 * A list of reads sharing a kmer.
 * @author Brian Bushnell
 * @date Nov 7, 2015
 *
 */
public class Clump extends ArrayList<Read> implements Comparable<Clump> {
	
	public static Clump makeClump(long kmer){
		try {
			return new Clump(kmer);
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
			throw new RuntimeException();
		}
	}
	
	private Clump(long kmer_){
		this(kmer_, 4);
	}

	private Clump(long kmer_, int size){
		super(size);
		kmer=kmer_;
	}

	@Override
	public boolean add(Read r){
		ReadKey key=(ReadKey) r.obj;
		assert(key.kmer==kmer);
		key.clump=this;
		return super.add(r);
	}
	
	private void setMaxima(){
		maxLeft=-1;
		maxRight=-1;
		width=-1;
		for(Read r : this){
			ReadKey key=(ReadKey) r.obj;
			final int pos=key.position;
			maxLeft=Tools.max(maxLeft, pos);
			maxRight=Tools.max(maxRight, r.length()-pos);
		}
		width=maxLeft+maxRight;
	}
	
	/** This will create counts of bases, or sums of qualities, at each position in the cluster. */
	private int[][] count(final boolean qualityScores){
		if(width<0){setMaxima();}
		
//		System.err.println("\n\n");
		final int[][] counts=new int[4][width];
		for(Read r : this){
			ReadKey key=(ReadKey) r.obj;
			final int pos=key.position;
			byte[] bases=r.bases, quals=r.quality;
			if(quals==null){useQuality=false;}
			
//			System.err.println("pos="+pos+", maxLeft="+maxLeft);
			for(int cloc=0, rloc=maxLeft-pos; cloc<bases.length; cloc++, rloc++){
//				System.err.println("cloc="+cloc+"/"+bases.length+", rloc="+rloc+"/"+width);
				int x=AminoAcid.baseToNumber[bases[cloc]];
				if(x>-1){
					final int q;
					if(qualityScores){q=(quals==null ? 20 : quals[cloc]);}
					else{q=1;}
					counts[x][rloc]+=q;
				}
			}
		}
		return counts;
	}
	
	/*--------------------------------------------------------------*/
	
	public Read makeSimpleConsensus(){
//		assert(Splitter.findBestPivot(this)<0) : Splitter.findBestPivot(this); //TODO: Slow
		if(size()==1){
			Read r=get(0);
			if(renameConsensus){r.id=r.numericID+" size=1";}
			return r;
		}
		final int[][] bcounts=baseCounts();
		final int[][] qcounts=qualityCounts();
		
		final byte[] bases=new byte[width], quals=new byte[width];
		for(int i=0; i<width; i++){
			int x=getConsensusAtPosition(qcounts, i);
			int y=getSecondHighest(qcounts, i);
			if(x>=0 && qcounts[x][i]==qcounts[y][i]){//Tie-breaker
				x=getConsensusAtPosition(bcounts, i);
				y=getSecondHighest(bcounts, i);
			}
			
			
			if(x<0){
//				System.err.println("q="+0+", x="+x+"; A="+counts[0][i]+", C="+counts[1][i]+", G="+counts[2][i]+", T="+counts[3][i]);
				bases[i]='N';
				quals[i]=0;
				assert(getSumAtPosition(qcounts, i)==0) : "\n"+bcounts[0][i]+", "+bcounts[1][i]+", "+bcounts[2][i]+", "+bcounts[3][i]+
														  "\n"+qcounts[0][i]+", "+qcounts[1][i]+", "+qcounts[2][i]+", "+qcounts[3][i]+
														  "\nwidth="+width+", i="+i+", size="+size()+"\n"+new String(bases, 0, i)+"\n"+get(0)+"\n"+get(1)+"\n";
//				assert(getSumAtPosition(bcounts, i)==0) : "\n"+bcounts[0][i]+", "+bcounts[1][i]+", "+bcounts[2][i]+", "+bcounts[3][i]+
//														  "\n"+qcounts[0][i]+", "+qcounts[1][i]+", "+qcounts[2][i]+", "+qcounts[3][i]+
//														  "\nwidth="+width+", i="+i+", size="+size()+"\n"+new String(bases, 0, i)+"\n"+get(0)+"\n"+get(1)+"\n";
			}else{
				final long bsum=getSumAtPosition(bcounts, i);
				final long bmajor=bcounts[x][i];
				final long bminor=bsum-bcounts[x][i];
				final long bsecond=bcounts[y][i];

				final long qsum=getSumAtPosition(qcounts, i);
				final long qmajor=qcounts[x][i];
				final long qminor=qsum-qcounts[x][i];
				final long qsecond=qcounts[y][i];
				
				float bratio=bminor/(float)(bmajor+bminor);
				double phred=(bminor==0 ? Read.MAX_CALLED_QUALITY() : -10*Math.log10(bratio));
				phred=Tools.min(phred, qmajor-qminor);
//				if(phred<0 || phred>127){
//					assert(false) :  i+","+x+","+bsum+","+qsum+","+bmajor+","+bminor+","+bratio+","+phred+"\n"+
//							bcounts[0][i]+","+bcounts[1][i]+","+bcounts[2][i]+","+bcounts[3][i]+"\n"+
//							qcounts[0][i]+","+qcounts[1][i]+","+qcounts[2][i]+","+qcounts[3][i]+"\n"+
//							this.toStringStaggered()+"\n";
//				}
//				assert(phred>=0 && phred<=127) : phred+","+x+","+i+","+bratio+","+bcounts[x][i]+","+bcounts[0][i]+
//					","+bcounts[1][i]+","+bcounts[2][i]+","+bcounts[3][i];
//				assert(phred>=0 && phred<=127) : bmajor+", "+bminor+", "+phred+", "+Read.MAX_CALLED_QUALITY;
				byte q=Read.capQuality((long)Math.round(phred));
				bases[i]=AminoAcid.numberToBase[x];
				quals[i]=q;
				assert(q>0);
				assert(x>-1);
				assert(bases[i]!='N');
			}
		}
		Read leftmost=this.get(0);//TODO:  Actually rightmost!
		Read r=new Read(bases, quals, leftmost.numericID+" size="+size(), 0);
		//TODO: Attach the long pair, and make sure the kmer location is correct.
//		assert(false) : "\n"+r.toFastq()+"\nCheck kmer location.";
//		assert(size()==1) : "\n"+r.toFastq()+"\n"+get(0).toFastq()+"\n"+get(size()-1).toFastq()+"\n";
		return r;
	}
	
	/*--------------------------------------------------------------*/
	
	public int removeDuplicates(){
		assert(KmerComparator.compareSequence);
		if(size()<2){return 0;}
		
		int removedTotal=0, removed=0;
		
		final boolean sortXY=(forceSortXY || sortYEarly() || (opticalOnly && (sortX || sortY) && size()>=sortXYSize));
		
		final int maxDiscarded;
		final int scan;
		
		if(opticalOnly && sortXY){
			scan=Tools.max(scanLimit, 200);
			maxDiscarded=scan+10;
		}else if(!sortXY && ((maxSubstitutions<1 && dupeSubRate<=0) || scanLimit<0)){
			scan=0;
			maxDiscarded=0;
		}else{
			scan=scanLimit;
			maxDiscarded=scan+10;
		}
		
		if(sortXY){
			assert(sortX || sortY);

			if(sortY){
				if(!sortYEarly()){
					this.sort(KmerComparatorY.comparator);
				}else{
//					for(int i=1; i<this.size(); i++){
//						Read a=get(i-1);
//						Read b=get(i);
//						assert(KmerComparatorY.comparator.compare(a, b)<=0) : a.obj+" vs "+b.obj;
//					}
//					//Otherwise, it was already Y-sorted.
				}
//				assert(false) : sortY();
				removed=removeDuplicates_inner(maxSubstitutions, dupeSubRate, scan, maxDiscarded, opticalOnly, true, markOnly, markAll, renameByCount, maxOpticalDistance);
				removedTotal+=removed;
//				System.err.println("RemovedY: "+removed);
				while((maxSubstitutions>0 || dupeSubRate>0) && scanLimit>0 && removed>maxDiscarded){
					removed=removeDuplicates_inner(maxSubstitutions, dupeSubRate, scan+10, maxDiscarded*2+20, opticalOnly, true, markOnly, markAll, renameByCount, maxOpticalDistance);
					removedTotal+=removed;
//					System.err.println("RemovedY: "+removed);
				}
			}
			if(sortX && (ReadKey.spanTilesX || !sortY)){
				this.sort(KmerComparatorX.comparator);
				removed=removeDuplicates_inner(maxSubstitutions, dupeSubRate, scan, maxDiscarded, opticalOnly, true, markOnly, markAll, renameByCount, maxOpticalDistance);
				removedTotal+=removed;
//				System.err.println("RemovedX: "+removed);
				while((maxSubstitutions>0 || dupeSubRate>0) && scanLimit>0 && removed>maxDiscarded){
					removed=removeDuplicates_inner(maxSubstitutions, dupeSubRate, scan+10, maxDiscarded*2+20, opticalOnly, true, markOnly, markAll, renameByCount, maxOpticalDistance);
					removedTotal+=removed;
//					System.err.println("RemovedX: "+removed);
				}
			}
		}else{
			removed=removeDuplicates_inner(maxSubstitutions, dupeSubRate, scan, maxDiscarded, opticalOnly, false, markOnly, markAll, renameByCount, maxOpticalDistance);
			removedTotal+=removed;
			while((maxSubstitutions>0 || dupeSubRate>0) && scanLimit>0 && removed>maxDiscarded){
				removed=removeDuplicates_inner(maxSubstitutions, dupeSubRate, scan+10, maxDiscarded*2+20, opticalOnly, false, markOnly, markAll, renameByCount, maxOpticalDistance);
				removedTotal+=removed;
			}
		}
		
		return removedTotal;
	}
	
	private int removeDuplicates_inner(final int maxSubs, final float subRate, final int scanLimit, final int maxDiscarded,
			final boolean optical, final boolean xySorted, final boolean mark, final boolean markAll, final boolean rename, final float dist){
		final int size=size();
		if(size<2){return 0;}
		int dupePairs=0;
		int dupeReads=0;
		
//		final boolean breakOnTile=(optical && !FlowcellCoordinate.spanTiles);
		
//		final long start=System.nanoTime();
		
		for(int i=0, lim=size-1; i<lim; i++){
			final Read a=get(i);
			if(!a.discarded()){
				final ReadKey keyA=(ReadKey) a.obj;
				final int aLen=a.length();
				int unequals=0;
				int discarded=0;
				for(int j=i+1; j<size && unequals<=scanLimit && discarded<=maxDiscarded && (!a.discarded() || markAll); j++){
					final Read b=get(j);
					if(b.discarded()){
						discarded++;
					}else{
						final int bLen=b.length();
						final ReadKey keyB=(ReadKey) b.obj;
						if(!containment && !keyA.equals(keyB)){break;}
//						if(containment && affix && !keyA.physicalAffix(keyB, aLen, bLen)){break;}
//						if(optical && keyA.lane!=keyB.lane){break;} //Already in equals method
//						if(breakOnTile && keyA.tile!=keyB.tile){break;} //Already in equals method
						if(optical && xySorted && !keyA.nearXY(keyB, dist)){break;}
						if(compareUMI && !keyA.umiMatches(keyB, umiSubs)) {break;}
//						if(System.nanoTime()-start>200000000000L){
//							TextStreamWriter tsw=new TextStreamWriter("foo.fq", true, false, false);
//							tsw.start();
//							for(Read x : this){
//								tsw.println(x.toFastq());
//							}
//							tsw.poisonAndWait();
//							assert(false) : "kmer="+kmer+", size="+size();
//						}
						if(equals(a, b, maxSubs, subRate)){
							if(!optical || keyA.near(keyB, dist)){
								if(printDuplicates){
									System.out.println(a.toFasta());
									System.out.println(b.toFasta());
								}
								float errA=a.expectedErrorsIncludingMate(true);
								float errB=b.expectedErrorsIncludingMate(true);
								if(markAll){
									b.setDiscarded(true);
									assert(!a.discarded() || markAll);
									dupePairs++;
									dupeReads+=1+b.mateCount();
									unequals=0;
									if(!a.discarded()){
										a.setDiscarded(true);
										dupePairs++;
										dupeReads+=1+a.mateCount();
									}
								}else if(containment || errB>=errA){
									b.setDiscarded(true);
									assert(!a.discarded() || markAll);
									a.copies+=b.copies+parseExtraCopies(b);
									dupePairs++;
									dupeReads+=1+b.mateCount();
									unequals=0;
								}else{
									a.setDiscarded(true);
									assert(!b.discarded() || markAll);
									b.copies+=a.copies+parseExtraCopies(a);
									dupePairs++;
									dupeReads+=1+a.mateCount();
								}
							}
						}else{
							unequals++;
						}
					}
				}
			}
		}
		if(dupeReads>0){
			for(int i=0; i<size; i++){
				Read a=get(i);
				if(a.discarded()){
					if(mark){
						if(!a.id.endsWith(" duplicate")){
							a.id+=" duplicate";
							if(a.mate!=null){a.mate.id+=" duplicate";}
						}
					}else{
						set(i, null);
					}
				}else if(rename && a.copies>1){
					renameFromCount(a);
				}
				a.copies=1;
			}
			if(!mark){
				int x=Tools.condenseStrict(this);
				assert(x==dupePairs) : size()+", "+size+", "+dupePairs+", "+x;
				assert((size()>0 || markAll) && size()==size-dupePairs) : size()+", "+size+", "+dupePairs;
			}
		}
		
		if(containment){
			dupeReads+=removeDuplicates_backwards(maxSubs, subRate, scanLimit, maxDiscarded, optical, xySorted, mark, markAll, rename, dist);
		}
		
		return dupeReads;
	}
	
	/** Only for containments */
	private int removeDuplicates_backwards(final int maxSubs, final float subRate, final int scanLimit, final int maxDiscarded,
			final boolean optical, final boolean xySorted, final boolean mark, final boolean markAll, final boolean rename, final float dist){
		final int size=size();
		if(size<2){return 0;}
		int dupePairs=0;
		int dupeReads=0;
		
//		final boolean breakOnTile=(optical && !FlowcellCoordinate.spanTiles);
		
//		final long start=System.nanoTime();
		
		for(int i=size-1; i>0; i--){
			final Read a=get(i);
			if(!a.discarded()){
				final ReadKey keyA=(ReadKey) a.obj;
				final int aLen=a.length();
				int unequals=0;
				int discarded=0;
				for(int j=i-1; j>=0 && unequals<=scanLimit && discarded<=maxDiscarded && (!a.discarded() || markAll); j--){
					final Read b=get(j);
					if(b.discarded()){
						discarded++;
					}else{
						final int bLen=b.length();
						final ReadKey keyB=(ReadKey) b.obj;
						if(!containment && !keyA.equals(keyB)){break;}
//						if(containment && affix && !keyA.physicalAffix(keyB, aLen, bLen)){break;}
//						if(optical && keyA.lane!=keyB.lane){break;} //Already in equals method
//						if(breakOnTile && keyA.tile!=keyB.tile){break;} //Already in equals method
						if(optical && xySorted && !keyA.nearXY(keyB, dist)){break;}
//						if(System.nanoTime()-start>200000000000L){
//							TextStreamWriter tsw=new TextStreamWriter("foo.fq", true, false, false);
//							tsw.start();
//							for(Read x : this){
//								tsw.println(x.toFastq());
//							}
//							tsw.poisonAndWait();
//							assert(false) : "kmer="+kmer+", size="+size();
//						}
						if(equals(a, b, maxSubs, subRate)){
							if(!optical || keyA.near(keyB, dist)){
								if(printDuplicates){
									System.out.println(a.toFasta());
									System.out.println(b.toFasta());
								}
								float errA=a.expectedErrorsIncludingMate(true);
								float errB=b.expectedErrorsIncludingMate(true);
								if(markAll){
									b.setDiscarded(true);
									assert(!a.discarded() || markAll);
									dupePairs++;
									dupeReads+=1+b.mateCount();
									unequals=0;
									if(!a.discarded()){
										a.setDiscarded(true);
										dupePairs++;
										dupeReads+=1+a.mateCount();
									}
								}else if(containment || errB>=errA){
									b.setDiscarded(true);
									assert(!a.discarded() || markAll);
									a.copies+=b.copies+parseExtraCopies(b);
									dupePairs++;
									dupeReads+=1+b.mateCount();
									unequals=0;
								}else{
									a.setDiscarded(true);
									assert(!b.discarded() || markAll);
									b.copies+=a.copies+parseExtraCopies(a);
									dupePairs++;
									dupeReads+=1+a.mateCount();
								}
							}
						}else{
							unequals++;
						}
					}
				}
			}
		}
		if(dupeReads>0){
			for(int i=0; i<size; i++){
				Read a=get(i);
				if(a.discarded()){
					if(mark){
						if(!a.id.endsWith(" duplicate")){
							a.id+=" duplicate";
							if(a.mate!=null){a.mate.id+=" duplicate";}
						}
					}else{
						set(i, null);
					}
				}else if(rename && a.copies>1){
					renameFromCount(a);
				}
				a.copies=1;
			}
			if(!mark){
				int x=Tools.condenseStrict(this);
				assert(x==dupePairs) : size()+", "+size+", "+dupePairs+", "+x;
				assert((size()>0 || markAll) && size()==size-dupePairs) : size()+", "+size+", "+dupePairs;
			}
		}
		return dupeReads;
	}
	
	public int parseExtraCopies(final Read a){
		final String id=a.id;
		final int space=id.lastIndexOf(' ');
		int extraCopies=0;
		if(space<0){return 0;}
		if(space>=0 && Tools.contains(id, "copies=", space+1)){
			extraCopies=Integer.parseInt(id.substring(space+8))-1;
		}
		return extraCopies;
	}
	
	private void renameFromCount(final Read a){
		final int newExtraCopies=a.copies-1;
		assert(newExtraCopies>0) : newExtraCopies;
		final int oldExtraCopies=parseExtraCopies(a);
		final int total=1+newExtraCopies+oldExtraCopies;
		renameToTotal(a, total);
		if(a.pairnum()==0 && a.mate!=null){
			assert(a.mate.pairnum()==1);
			renameToTotal(a.mate, total);
		}
	}
	
	private static void renameToTotal(final Read a, final int total){
		final String id=a.id;
		final int space=id.lastIndexOf(' ');
		if(space>=0 && Tools.contains(id, "copies=", space+1)){
			a.id=a.id.substring(0, space);
		}
		a.id=a.id+" copies="+total;
	}
	
//	public boolean nearby(Read a, Read b, FlowcellCoordinate fca, FlowcellCoordinate fcb, float maxDist){
////		fca.setFrom(a.id);
//		fcb.setFrom(b.id);
//		float dist=fca.distance(fcb);
//		return dist<=maxDist;
//	}
	
//	public boolean nearby(ReadKey ka, ReadKey kb, float maxDist){
//		float dist=ka.distance(kb);
//		return dist<=maxDist;
//	}
	
	public static boolean equals(Read a, Read b, int maxSubs, float dupeSubRate){
		if(a.numericID==b.numericID){return false;}
		if(dupeSubRate>0){maxSubs=Tools.max(maxSubs, (int)(dupeSubRate*Tools.min(a.length(), b.length())));}
		if(containment){
			return contains(a, b, maxSubs);
		}
		if(!equals(a.bases, b.bases, maxSubs)){return false;}
		if(a.mate!=null && !equals(a.mate.bases, b.mate.bases, maxSubs)){return false;}
		return true;
	}
	
	public static boolean equals(byte[] a, byte[] b, int maxSubs){
		if(a==b){return false;}//Nothing should subsume itself
		if(a==null || b==null){return false;}
		if(a.length!=b.length){return false;}
		int subs=0;
		if(allowNs){
			for(int i=0; i<a.length; i++){
				if(a[i]!=b[i] && (a[i]!='N' && b[i]!='N')){
					subs++;
					if(subs>maxSubs){return false;}
				}
			}
		}else{
			for(int i=0; i<a.length; i++){
				if(a[i]!=b[i]){
					subs++;
					if(subs>maxSubs){return false;}
				}
			}
		}
		return true;
	}
	
	public static boolean contains(Read a, Read b, int maxSubs){
		if(a.numericID==b.numericID){return false;}
		boolean ok=contains_inner(a, b, maxSubs);
		if(!ok || a.mate==null){return ok;}
		ok=contains_inner(a.mate, b.mate, maxSubs);
		if(!ok){return false;}
		ReadKey rka1=(ReadKey)a.obj;
		ReadKey rkb1=(ReadKey)b.obj;
		ReadKey rka2=(ReadKey)a.mate.obj; //TODO: In containment mode, mates need to always get keys.
		ReadKey rkb2=(ReadKey)b.mate.obj;
		return ((rka1.kmerMinusStrand==rkb1.kmerMinusStrand) && (rka2.kmerMinusStrand==rkb2.kmerMinusStrand)); //Ensures that both reads have the same directionality.
	}
	
	public static boolean contains_inner(Read a, Read b, int maxSubs){
//		if(a.length()==b.length()){return equals(a.bases, b.bases, maxSubs);}
		ReadKey rka=(ReadKey)a.obj;
		ReadKey rkb=(ReadKey)b.obj;
		if(affix ? !rka.physicalAffix(rkb, a.length(), b.length()) : !rka.physicallyContains(rkb, a.length(), b.length())){return false;}
		boolean flipped=false;
//		if(a.mate!=null){//In paired mode, need synchronization if the reads are in difference clumps.  But...  just don't do that.
//			Read max, min;
//			if(a.numericID>b.numericID){max=a; min=b;}//Protects against deadlocks.
//			else{max=b; min=a;}
//			synchronized(min){
//				synchronized(max){
//					if(rka.kmerMinusStrand!=rkb.kmerMinusStrand){
//						rkb.flip(b, k);
//						flipped=true;
//					}
//					boolean ok=contains(a.bases, b.bases, rka.position, rkb.position, maxSubs);
//					if(flipped){rkb.flip(b, k);}
//					return ok;
//				}
//			}
//		}
		if(rka.kmerMinusStrand!=rkb.kmerMinusStrand){
			rkb.flip(b, k);
			flipped=true;
		}
		boolean ok=contains(a.bases, b.bases, rka.position, rkb.position, maxSubs);
		if(flipped){rkb.flip(b, k);}
		return ok;
	}
	
	public static boolean contains(byte[] a, byte[] b, int posA, int posB, int maxSubs){
		if(a==null || b==null){return false;}
		assert(a.length>=b.length);
		if(a==b){return false;} //Nothing should contain itself
		int subs=0;
		
		int ai, bi;
		final int dif=posA-posB;
		if(dif>0){
			ai=0;
			bi=-dif;
		}else{
			ai=dif;
			bi=0;
		}
		
		if(allowNs){
			for(; ai<a.length && bi<b.length; ai++, bi++){
				if(ai>=0 && bi>=0){
					final byte na=a[ai], nb=b[bi];
					if(na!=nb && na!='N' && nb!='N'){
						subs++;
						if(subs>maxSubs){return false;}
					}
				}
			}
		}else{
			for(; ai<a.length && bi<b.length; ai++, bi++){
				if(ai>=0 && bi>=0 && a[ai]!=b[bi]){
					subs++;
					if(subs>maxSubs){return false;}
				}
			}
		}
		return true;
	}
	
	/*--------------------------------------------------------------*/
	
	public long splitAndErrorCorrect(){
		if(size()<Splitter.minSizeSplit){
			return errorCorrect();
		}
		long sum=0;
		ArrayList<Clump> list=Splitter.splitOnPivot(this);
		for(Clump c : list){
			sum+=c.errorCorrect();
		}
		return sum;
	}
	
	public long errorCorrect(){
		if(size()<=minCountCorrect){return 0;}
//		assert(Splitter.findBestPivot(this)<0); //TODO: Slow
		Read consen=makeSimpleConsensus();
		long sum=0;
		final int[] rvector=KillSwitch.allocInt1D(2);
		for(Read r : this){
			sum+=errorCorrect(r, consen, rvector);
		}
		return sum;
	}
	
	private int errorCorrect(Read call, Read ref, int[] rvector){

//		assert(call.validated());
//		assert(call.checkQuality());
//		assert(ref.checkQuality());
		
		final float identity=identity(call, ref.bases, rvector);
		if((identity<minIdentity && (rvector[1]>0 || rvector[0]<50)) || (identity==1 && !call.containsNocalls()/* && !ref.containsNocalls()*/)){
			//TODO: Strange, the identity==1 clause causes different behavior, though it does give a speedup.
			return 0;
		}
		final byte[] cbases=call.bases, cquals=call.quality;
		final byte[] rbases=ref.bases, rquals=ref.quality;
		
		ReadKey key=(ReadKey) call.obj;
		final int pos=key.position;
		
		final int[][] bcounts=baseCounts();
		final int[][] qcounts=qualityCounts();
		final float[][] qAvgs=qualityAverages();
		
		final int minDepth=(int)Tools.max(minCountCorrect, minSizeFractionCorrect*size());
		
		int corrections=0;
		
		final int cStart=0, rStart=maxLeft-pos, max=cbases.length;
		for(int cloc=cStart, rloc=rStart; cloc<max; cloc++, rloc++){
			//Called base, ref base
			final byte cb=cbases[cloc], rb=rbases[rloc];
			//Called quality, ref quality
			final byte cq=(cquals==null ? 20 : cquals[cloc]), rq=rquals[rloc];
			//Called number
			final byte cx=AminoAcid.baseToNumber[cb];
			//Ref number
			final byte rx=AminoAcid.baseToNumber[rb];
			
//			assert((cb=='N') == (cquals[cloc]==0));
			
			final byte b, q;
			if(cx<0){
				b=rb;
				q=(byte)Tools.min(20, rq);
			}else if(cb==rb){
				b=cb;
				q=(byte)Tools.mid(cq, cq+1, rq);//Adjust upward
				assert(q>=cq && (q<=rq || q<=cq));
			}else{
				final int cCount=bcounts[cx][rloc];
				final int rCount=bcounts[rx][rloc];
				final int altQSum=qcounts[cx][rloc];
				final int rQSum=qcounts[rx][rloc];
				final float rQAvg=qAvgs[rx][rloc];
				if(cCount<=maxCIncorrect && cq<=maxQIncorrect && cq*minQRatio<rQSum && cq*minAQRatio<8+rQAvg){
					final byte pminor=getSecondHighest(bcounts, rloc);

					assert(rx>=0 && rx<bcounts.length) : rx+", "+rloc+", "+bcounts.length+"\n"+call.toFastq()+"\n"+ref.toFastq();
					assert(rloc>=0 && rloc<bcounts[rx].length) : rx+", "+rloc+", "+bcounts[rloc].length+"\n"+call.toFastq()+"\n"+ref.toFastq();
					final int minorCount=bcounts[pminor][rloc];

					final long sum=getSumAtPosition(bcounts, rloc);
					//				final long altCount=sum-rCount;

					boolean change=false;
					if(rCount>=minCountCorrect && sum>=minDepth){
						final float ratioNeeded=Tools.min(minRatio, minRatioMult*minorCount+minRatioOffset+minRatioQMult*cq);
//						final float qratioNeeded=Tools.min(minRatio, minRatioMult*altQSum+minRatioOffset+minRatioQMult*cq); //altQSum is slightly different than minorQCount
						if(minorCount*ratioNeeded<=rCount && altQSum*ratioNeeded<=rQSum){
							change=true;
						}
						
//						else if(minorCount*ratioNeeded<=rCount){
//							assert(false) : "\n"+altQSum+", "+rQSum+", "+qratioNeeded+"\n"+cCount+", "+rCount+", "+sum+", "+ratioNeeded+"\n"+(altQSum*qratioNeeded);
//						}
					}
					if(change){
						corrections++;
						b=rb;
						q=(byte)Tools.min(cq+1, 25, rq);
						//					assert(q==25 || (q<=rq && q>=cq)) : q+", "+cq+", "+rq;
					}else{
						b=cb;
						q=(byte)Tools.mid(cq, cq-1, 6);//Adjust downward
						assert(q<=cq || q>=6) : q+","+cq;
					}
				}else{
					b=cb;
					q=cq;
				}
			}
			cbases[cloc]=b;
			if(cquals!=null){
				byte q2=(byte)Tools.mid(q, cq+maxQAdjust, cq-maxQAdjust);
				if(q2==0 && AminoAcid.isFullyDefined(b)){
					assert(!AminoAcid.isFullyDefined(cb));
					q2=(byte)Tools.mid(2, 25, (rq+25)/2);
				}else if(!AminoAcid.isFullyDefined(b)){
					q2=0;
				}
				cquals[cloc]=q2;
				assert((b=='N') == (cquals[cloc]==0)) : "b="+(char)b+", cb="+(char)cb+", rb="+(char)rb+", cx="+cx+", "
						+ "new="+cquals[cloc]+", q="+q+", old="+cq+", rq="+rq+", loc="+rloc+"\n"+call.toFastq()+"\n"+ref.toFastq();
			}
		}
		return corrections;
	}
	
	/*--------------------------------------------------------------*/
	
	//Only used by condense mode.
	public ArrayList<Read> makeConsensus(){
		if(size()==1){
			Read r=get(0);
			r.id=r.numericID+" size=1";
			return this;
		}
		
		ArrayList<Clump> clumps=Splitter.splitOnPivot(this);//TODO: Really, this should return null if there is no pivot

		ArrayList<Read> list=new ArrayList<Read>(clumps.size());
		for(Clump c : clumps){
			Read temp=c.makeSimpleConsensus();
			list.add(temp);
		}
		return list;
	}
	
	/*--------------------------------------------------------------*/
	
	private float identity(Read call, byte[] rbases, int[] rvector){
		ReadKey key=(ReadKey) call.obj;
		final int pos=key.position;
		byte[] cbases=call.bases, quals=call.quality;
		int good=0, bad=0;
		
		final int cStart=0, rStart=maxLeft-pos, max=cbases.length;
		for(int cloc=cStart, rloc=rStart; cloc<max; cloc++, rloc++){
			final byte cb=cbases[cloc], rb=rbases[rloc];
			if(AminoAcid.isFullyDefined(cb) && AminoAcid.isFullyDefined(rb)){
				if(cb==rb){good++;}
				else{bad++;}
			}
		}
		rvector[0]=good;
		rvector[1]=bad;
		return good==0 ? 0 : good/(float)(good+bad);
	}
	
	/*--------------------------------------------------------------*/
	
	long getSumAtPosition(int[][] counts, int pos){
		long sum=0;
		for(int x=0; x<4; x++){
			sum+=counts[x][pos];
		}
		return sum;
	}
	
	byte getConsensusAtPosition(int[][] counts, int pos){
		byte xMax=0;
		for(byte x=1; x<4; x++){
//			System.err.println("x="+x+", max="+max+", Checking "+counts[x][pos]+" vs "+counts[x][max]);
			if(counts[x][pos]>counts[xMax][pos]){xMax=x;}
		}
//		assert(counts[max][pos]>=counts[0][pos]);
//		assert(counts[max][pos]>=counts[1][pos]);
//		assert(counts[max][pos]>=counts[2][pos]) : max+", "+counts[max][pos]+", ["+counts[0][pos]+", "+counts[1][pos]+", "+counts[2][pos]+", "+counts[3][pos]+"]";
//		assert(counts[max][pos]>=counts[3][pos]);
		return (counts[xMax][pos]>0 ? xMax : -1);
	}
	
	byte getSecondHighest(int[][] counts, int pos){
		byte first=0;
		byte second=1;
		if(counts[first][pos]<counts[second][pos]){
			first=1; second=0;
		}
		for(byte x=2; x<4; x++){
//			System.err.println("x="+x+", max="+max+", Checking "+counts[x][pos]+" vs "+counts[x][max]);
			if(counts[x][pos]>counts[first][pos]){
				second=first;
				first=x;
			}else if(counts[x][pos]>counts[second][pos]){
				second=x;
			}
		}
//		assert(counts[max][pos]>=counts[0][pos]);
//		assert(counts[max][pos]>=counts[1][pos]);
//		assert(counts[max][pos]>=counts[2][pos]) : max+", "+counts[max][pos]+", ["+counts[0][pos]+", "+counts[1][pos]+", "+counts[2][pos]+", "+counts[3][pos]+"]";
//		assert(counts[max][pos]>=counts[3][pos]);
		
		return second; //may be actually 0.
		//return (counts[second][pos]>0 ? second : -1);
	}
	
	/*--------------------------------------------------------------*/
	
	public String toStringStaggered(){
		ByteBuilder sb=new ByteBuilder();
		for(Read r : this){
			ReadKey key=(ReadKey) r.obj;
			final int pos=key.position;
			byte[] bases=r.bases, quals=r.quality;
			int rloc=0, cloc=rloc+pos-maxLeft;
			while(cloc<0){sb.append(' '); cloc++;}
			sb.append(bases);
			sb.append('\n');
		}
		return sb.toString();
	}
	
	public Read consensusRead(){
		if(consensusRead==null){
			consensusRead=makeSimpleConsensus();
		}
		return consensusRead;
	}
	
	public int width(){
		assert(width>=0) : width;
		return width;
	}
	
	public int maxLeft(){
		assert(maxLeft>=0);
		return maxLeft;
	}
	
	public int maxRight(){
		assert(maxRight>=0);
		return maxRight;
	}
	
	/*--------------------------------------------------------------*/
	
	int[][] baseCounts(){
		if(baseCounts==null){
			baseCounts=count(false);
			int len=baseCounts[0].length;
			assert(width==-1 || width==len);
		}
		return baseCounts;
	}
	
	int[][] qualityCounts(){
		if(qualityCounts==null){
			qualityCounts=count(true);
			int len=qualityCounts[0].length;
			assert(width==-1 || width==len);
		}
		return qualityCounts;
	}
	
	float[][] qualityAverages(){
		if(qualityAverages==null){
			qualityAverages=new float[4][width];
			for(int i=0; i<4; i++){
				for(int j=0; j<width; j++){
					int b=baseCounts[i][j];
					int q=qualityCounts[i][j];
					qualityAverages[i][j]=(b==0 ? 0 : q/(float)b);
				}
			}
		}
		return qualityAverages;
	}

	void clearCounts(){
		baseCounts=qualityCounts=null;
		qualityAverages=null;
	}
	
	private void clearConsensus(){
		consensusRead=null;
	}
	
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean equals(Object b){
		return this==b;
	}
	
	@Override
	public int hashCode(){
		return Long.hashCode(kmer);
	}

	@Override
	public int compareTo(Clump o) {
		int x=Long.compare(kmer, o.kmer);
		return x!=0 ? x : o.size()-size();
	}
	
	/*--------------------------------------------------------------*/
	
	public final long kmer;
	
	private int width=-1;
	private int maxLeft=-1;
	private int maxRight=-1;
	
	private int[][] baseCounts;
	private int[][] qualityCounts;
	private float[][] qualityAverages;
	
	private Read consensusRead;
	
	boolean useQuality(){return useQuality;}
	private boolean useQuality=true;
	
	boolean added=false;
	
	public static final int overhead=overhead();
	private static int overhead(){
		return 16+ //self
				16+ //Backing array
				4+ //Backing array size
				4*(8)+ //Backing array initial capacity
				1*(8)+ //kmer
				3*(4)+ //int fields
				4*(8)+ //pointers
				2*(4); //booleans
	}
	
	/*--------------------------------------------------------------*/
	
	public static boolean parseStatic(String arg, String a, String b){
		if(a.equals("mincountcorrect") || a.equals("mincc")){
			minCountCorrect=Integer.parseInt(b);
		}else if(a.equals("minsizesplit") || a.equals("minss")){
			Splitter.minSizeSplit=Integer.parseInt(b);
		}else if(a.equals("minsizefractionsplit") || a.equals("minsfs")){
			Splitter.minSizeFractionSplit=Float.parseFloat(b);
		}else if(a.equals("minsizefractioncorrect") || a.equals("minsfc")){
			minSizeFractionCorrect=Float.parseFloat(b);
		}else if(a.equals("minratio") || a.equals("minr")){
			minRatio=Float.parseFloat(b);
		}else if(a.equals("minqratio") || a.equals("minqr")){
			minQRatio=Float.parseFloat(b);
		}else if(a.equals("minaqratio") || a.equals("minaqr")){
			minAQRatio=Float.parseFloat(b);
		}else if(a.equals("minratiooffset") || a.equals("minro")){
			minRatioOffset=Float.parseFloat(b);
		}else if(a.equals("minratiomult") || a.equals("minrm")){
			minRatioMult=Float.parseFloat(b);
		}else if(a.equals("minratioqmult") || a.equals("minrqm")){
			minRatioQMult=Float.parseFloat(b);
		}else if(a.equals("minidentity") || a.equals("minid")){
			minIdentity=Float.parseFloat(b);
		}else if(a.equals("maxqadjust")){
			maxQAdjust=(byte)Integer.parseInt(b);
		}else if(a.equals("maxqi") || a.equals("maxqincorrect") || a.equals("maxqualityincorrect")){
			maxQIncorrect=Integer.parseInt(b);
			if(maxCIncorrect<0){maxQIncorrect=Integer.MAX_VALUE;}
		}else if(a.equals("maxci") || a.equals("maxcincorrect") || a.equals("maxcountincorrect")){
			maxCIncorrect=Integer.parseInt(b);
			if(maxCIncorrect<0){maxCIncorrect=Integer.MAX_VALUE;}
		}else if(a.equals("border")){
			KmerComparator.defaultBorder=Integer.parseInt(b);
		}else if(a.equals("conservative")){
			conservativeFlag=Parse.parseBoolean(b);
		}else if(a.equals("aggressive")){
			aggressiveFlag=Parse.parseBoolean(b);
		}else if(a.equals("forceprocess")){
			forceProcess=Parse.parseBoolean(b);
		}else if(a.equals("mergefirst")){
			KmerComparator.mergeFirst=Parse.parseBoolean(b);
		}else if(a.equals("findcorrelations")){
			Splitter.FIND_CORRELATIONS=Parse.parseBoolean(b);
		}else if(a.equals("maxcorrelations")){
			Splitter.MAX_CORRELATIONS=Integer.parseInt(b);
		}
		
		else if(a.equals("sortx")){
			sortX=Parse.parseBoolean(b);
		}else if(a.equals("sorty")){
			sortY=Parse.parseBoolean(b);
		}else if(a.equals("sortxy")){
			sortX=sortY=Parse.parseBoolean(b);
		}else if(a.equals("forcesortxy") || a.equals("forcexy")){
			forceSortXY=Parse.parseBoolean(b);
		}else if(a.equals("sortxysize") || a.equals("xysize")){
			sortXYSize=Integer.parseInt(b);
		}
		
		else if(a.equals("removeallduplicates") || a.equals("allduplicates")){
			markAll=Parse.parseBoolean(b);
		}else if(a.equals("addcount") || a.equals("renamebycount")){
			renameByCount=Parse.parseBoolean(b);
		}else if(a.equals("optical") || a.equals("opticalonly")){
			opticalOnly=Parse.parseBoolean(b);
		}else if(a.equals("dupesubs") || a.equals("duplicatesubs") || a.equals("dsubs") || a.equals("subs") || a.equals("s")){
			maxSubstitutions=Integer.parseInt(b);
		}else if(a.equals("dupedist") || a.equals("duplicatedistance") || a.equals("ddist") || a.equals("dist") || a.equals("opticaldist") || a.equals("distance")){
			try {
				maxOpticalDistance=Float.parseFloat(b);
			} catch (NumberFormatException e) {
				maxOpticalDistance=Parse.parseKMG(b);
			}
			opticalOnly=maxOpticalDistance>=0;
		}else if(a.equals("scanlimit") || a.equals("scan")){
			scanLimit=Integer.parseInt(b);
		}else if(a.equals("allowns")){
			allowNs=Parse.parseBoolean(b);
		}else if(a.equals("containment") || a.equals("absorbcontainment") || a.equals("ac") || a.equals("contains")){
			containment=Parse.parseBoolean(b);
		}else if(a.equals("compareumi") || a.equals("umi")){
			compareUMI=FlowcellCoordinate.parseUMI=Parse.parseBoolean(b);
		}else if(a.equals("umisubs")){
			umiSubs=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("prefixOrSuffix") || a.equalsIgnoreCase("suffixOrPrefix") || a.equals("affix") || a.equals("pos")){
			affix=Parse.parseBoolean(b);
		}else if(a.equals("printduplicates")){
			printDuplicates=Parse.parseBoolean(b);
		}else if(a.equals("dupeidentity")){
			float dupeIdentity=Float.parseFloat(b);
			if(dupeIdentity>1){dupeIdentity=dupeIdentity/100;}
			assert(dupeIdentity<=1);
			dupeSubRate=1-dupeIdentity;
		}else if(a.equals("dupesubrate") || a.equals("dsr") || a.equals("subrate")){
			dupeSubRate=Float.parseFloat(b);
			if(dupeSubRate>1){dupeSubRate=dupeSubRate/100;}
			assert(dupeSubRate<=1);
		}
		
		else if(a.equals("minsizesplit")){
			Splitter.minSizeSplit=Integer.parseInt(b);
		}else if(a.equals("minsizefractionsplit")){
			Splitter.minSizeFractionSplit=Float.parseFloat(b);
		}else{
			return false;
		}
		
		return true;
	}
	
	static void setConservative(boolean newState){
		
		if(aggressiveFlag){return;}
		if(newState==conservativeMode){return;}
		conservativeMode=newState;
		
		Splitter.conservative=conservativeMode;
		
		if(conservativeMode){
			minCountCorrect++;
			minSizeFractionCorrect*=1.5f;
			minRatio*=1.25f;
			minQRatio*=1.5f;
			minAQRatio*=1.4f;
			minRatioOffset*=1.25f;
			minRatioQMult*=1.50f;
			minRatioMult*=1.4f;
			minIdentity=1-((1-minIdentity)*0.7f);
			if(maxQIncorrect==Integer.MAX_VALUE){maxQIncorrect=35;}
		}else{
			minCountCorrect--;
			minSizeFractionCorrect/=1.5f;
			minRatio/=1.25f;
			minQRatio/=1.5f;
			minAQRatio/=1.4f;
			minRatioOffset/=1.25f;
			minRatioQMult/=1.50f;
			minRatioMult/=1.4f;
			minIdentity=1-((1-minIdentity)/0.7f);
			if(maxQIncorrect==35){maxQIncorrect=Integer.MAX_VALUE;}
		}
	}
	
	/*--------------------------------------------------------------*/

	public static void setXY() {
		if(ReadKey.spanTilesX){sortX=true;}
		if(ReadKey.spanTilesY){sortY=true;}
	}

	static boolean allowNs=true;
	static boolean markAll=false;
	static boolean markOnly=false;
	static boolean opticalOnly=false;
	static boolean containment=false;
	static boolean compareUMI=false;
	static int umiSubs=0;
	static boolean affix=false;
	static boolean printDuplicates=false; //For debugging
	
	private static boolean renameByCount=false;
	private static int scanLimit=5;
	private static int maxSubstitutions=2;
	private static float dupeSubRate=0;
	private static float maxOpticalDistance=40;
	
	static boolean forceProcess=false;
	static boolean conservativeFlag=false;
	static boolean aggressiveFlag=false;
	static boolean conservativeMode=false;
	static boolean renameConsensus=false;
	static int minCountCorrect=4; //mcc=4 was slightly better than 3
	static float minSizeFractionCorrect=0.20f; //0.11 is very slightly better than 0.14
	static float minRatio=30.0f;
	static float minQRatio=2.8f; //Does nothing?
	static float minAQRatio=0.7f;
	static float minRatioOffset=1.9f; //3 is far worse than 2; 1 is better
	static float minRatioQMult=0.08f;
	static float minRatioMult=1.80f; //2.5 is far worse than 2; 1.5 is better
	static float minIdentity=0.97f; //0.98 is slightly better than 0.96; 0.94 is substantially worse
	static byte maxQAdjust=0;
	static int maxQIncorrect=Integer.MAX_VALUE;
	static int maxCIncorrect=Integer.MAX_VALUE;
	
	static boolean sortX=false; //Not needed for NextSeq
	static boolean sortY=true;
	static boolean forceSortXY=false; //Mainly for testing
	static int sortXYSize=6;
	
	/** May slightly increase speed for optical dedupe.  Can be safely disabled. */
	static boolean sortYEarly(){return sortY && (forceSortXY || opticalOnly);}
	
//	private static final boolean countQuality=false;
	public static int k=31;
	private static final long serialVersionUID = 1L;
	
}
