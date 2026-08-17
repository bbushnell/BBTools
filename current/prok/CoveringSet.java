package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.LongHashSet;
import map.LongIntMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.Read;
import structures.LongList;
import structures.ListNum;

/**
 * Greedy set-cover kmer selection for tRNA (or any short-gene) covering sets.
 * Replaces the bash kmercountexact+bbduk descent loop with a single in-memory
 * pass: count all kmers across a pool of sequences, then iteratively select the
 * top-N by count, evict covered sequences, recount, and repeat until a coverage
 * target or kmer budget is reached.
 *
 * Two-tier ranking (Brian's design): each round takes the top 2*step kmers by
 * CURRENT count on the remaining pool, re-sorts them by ORIGINAL prevalence
 * (recorded before any eviction), and keeps the top step. This concentrates on
 * kmers shared between rare and common tRNAs rather than random sequence that
 * happens to be common among the rare stragglers.
 *
 * @author Neptune, Brian Bushnell
 */
public class CoveringSet {

	public static void main(String[] args){
		Timer t=new Timer();
		CoveringSet x=new CoveringSet(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public CoveringSet(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(a.equals("kdesign")){
				kDesign=Integer.parseInt(b);
			}else if(a.equals("step")){
				step=(int)Parse.parseKMG(b);
			}else if(a.equals("maxkmers") || a.equals("budget")){
				maxKmers=(int)Parse.parseKMG(b);
			}else if(a.equals("mincov") || a.equals("mincoverage") || a.equals("target")){
				minCovFraction=Float.parseFloat(b);
			}else if(a.equals("extra") || a.equals("extra1")){
				extra=b;
			}else if(a.equals("copies")){
				copies=Integer.parseInt(b);
			}else if(a.equals("rcomp")){
				rcomp=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//handled
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		in=parser.in1;
		out=parser.out1;
		if(in==null){throw new RuntimeException("Error - an input file is required.");}
		if(out==null){throw new RuntimeException("Error - an output file is required.");}

		if(kDesign<0){kDesign=k;}
		assert(kDesign>=k) : "kdesign must be >= k";
		assert(k>0 && k<=31) : "k must be 1..31";
		assert(kDesign<=31) : "kdesign must be <= 31";

		overwrite=parser.overwrite;
		ffin=FileFormat.testInput(in, FileFormat.FASTA, null, true, true);
		ffout=FileFormat.testOutput(out, FileFormat.FASTA, null, true, overwrite, false, false);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	void process(Timer t){
		ArrayList<byte[]> pool=loadPool();
		final int totalSeqs=pool.size();
		outstream.println("Pool: "+totalSeqs+" sequences, k="+k+
				(kDesign!=k ? " (design k="+kDesign+")" : "")+
				", step="+step+", rcomp="+rcomp);

		final long designMask=~((-1L)<<(2*kDesign));
		final LongIntMap originalCounts=countKmers(pool, kDesign, designMask);
		outstream.println("Distinct "+kDesign+"-mers in pool: "+originalCounts.size());

		LongList selectedKmers=new LongList();
		LongHashSet selectedSet=new LongHashSet(maxKmers>0 ? maxKmers : 16384);
		boolean[] alive=new boolean[pool.size()];
		Arrays.fill(alive, true);
		int aliveCount=totalSeqs;
		int round=0;

		while(aliveCount>0){
			if(maxKmers>0 && selectedKmers.size>=maxKmers){break;}
			float covFrac=1f-(aliveCount/(float)totalSeqs);
			if(covFrac>=minCovFraction){break;}

			LongIntMap currentCounts=countKmersAlive(pool, alive, kDesign, designMask);
			if(currentCounts.isEmpty()){break;}

			long[] candidates=topNByCount(currentCounts, step*2);
			long[] ranked=rankByOriginal(candidates, originalCounts, step);

			int added=0;
			for(long kmer : ranked){
				if(selectedSet.contains(kmer)){continue;}
				selectedSet.add(kmer);
				selectedKmers.add(kmer);
				added++;
			}
			if(added==0){break;}

			int evicted=evictCovered(pool, alive, selectedSet, kDesign, designMask);
			aliveCount-=evicted;
			round++;

			float cov=1f-(aliveCount/(float)totalSeqs);
			if(verbose || round%10==0){
				outstream.println("  Round "+round+": +"+added+" kmers (total "+
						selectedKmers.size+"), evicted "+evicted+
						", alive="+aliveCount+", coverage="+String.format("%.4f", cov));
			}
		}

		float finalCov=1f-(aliveCount/(float)totalSeqs);
		outstream.println("Selected "+selectedKmers.size+" "+kDesign+"-mers in "+round+
				" rounds, coverage="+String.format("%.6f", finalCov)+
				" ("+aliveCount+" uncovered of "+totalSeqs+")");

		ArrayList<Read> output;
		if(kDesign==k){
			output=kmersToReads(selectedKmers, k);
		}else{
			output=designToUseKmers(selectedKmers, kDesign, k);
		}
		outstream.println("Output: "+output.size()+" "+k+"-mer sequences");

		if(ffout!=null){
			ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ffout, null, 4, null, false);
			ros.start();
			ros.add(output, 0);
			errorState|=ReadWrite.closeStream(ros);
		}

		t.stop();
		outstream.println("Time:   "+t);
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private ArrayList<byte[]> loadPool(){
		ArrayList<byte[]> pool=new ArrayList<>();
		loadFasta(in, pool);
		if(extra!=null){
			int base=pool.size();
			ArrayList<byte[]> extraSeqs=new ArrayList<>();
			loadFasta(extra, extraSeqs);
			for(int c=0; c<copies; c++){
				for(byte[] seq : extraSeqs){pool.add(seq);}
			}
			outstream.println("Added "+extraSeqs.size()+" extra sequences x"+copies+
					" copies (pool "+base+" -> "+pool.size()+")");
		}
		return pool;
	}

	private void loadFasta(String fname, ArrayList<byte[]> dest){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ff, null);
		cris.start();
		for(ListNum<Read> ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()){
			for(Read r : ln){
				if(r.bases!=null && r.length()>0){
					Tools.toUpperCase(r.bases);
					dest.add(r.bases);
					if(rcomp){
						byte[] rc=AminoAcid.reverseComplementBases(r.bases);
						dest.add(rc);
					}
				}
			}
			cris.returnList(ln);
		}
		ReadWrite.closeStream(cris);
	}

	private LongIntMap countKmers(ArrayList<byte[]> pool, int kLen, long mask){
		LongIntMap map=new LongIntMap((int)Tools.min(1<<20, pool.size()*50L));
		for(byte[] bases : pool){
			long kmer=0;
			int len=0;
			for(byte b : bases){
				int num=AminoAcid.baseToNumber[b];
				if(num>=0){
					kmer=((kmer<<2)|num)&mask;
					len++;
					if(len>=kLen){map.increment(kmer);}
				}else{len=0; kmer=0;}
			}
		}
		return map;
	}

	private LongIntMap countKmersAlive(ArrayList<byte[]> pool, boolean[] alive,
			int kLen, long mask){
		LongIntMap map=new LongIntMap((int)Tools.min(1<<20, pool.size()*10L));
		for(int i=0; i<pool.size(); i++){
			if(!alive[i]){continue;}
			byte[] bases=pool.get(i);
			long kmer=0;
			int len=0;
			for(byte b : bases){
				int num=AminoAcid.baseToNumber[b];
				if(num>=0){
					kmer=((kmer<<2)|num)&mask;
					len++;
					if(len>=kLen){map.increment(kmer);}
				}else{len=0; kmer=0;}
			}
		}
		return map;
	}

	/** Returns the top N kmers by count from the map, sorted descending. */
	private long[] topNByCount(LongIntMap map, int n){
		long[] keys=map.keys();
		int[] counts=new int[keys.length];
		for(int i=0; i<keys.length; i++){
			counts[i]=map.get(keys[i]);
		}
		//Sort by count descending: co-sort keys by counts
		int valid=keys.length;
		if(valid<=n){
			coSortDescending(keys, counts, valid);
			return keys;
		}
		coSortDescending(keys, counts, valid);
		return Arrays.copyOf(keys, n);
	}

	/** From the candidate list, re-rank by original prevalence, keep top step. */
	private long[] rankByOriginal(long[] candidates, LongIntMap originalCounts, int keep){
		int[] origCounts=new int[candidates.length];
		for(int i=0; i<candidates.length; i++){
			origCounts[i]=originalCounts.get(candidates[i]);
		}
		coSortDescending(candidates, origCounts, candidates.length);
		return candidates.length<=keep ? candidates : Arrays.copyOf(candidates, keep);
	}

	/** Marks sequences as not-alive if they contain any selected kmer.
	 * Returns number of newly evicted sequences. */
	private int evictCovered(ArrayList<byte[]> pool, boolean[] alive,
			LongHashSet selected, int kLen, long mask){
		int evicted=0;
		for(int i=0; i<pool.size(); i++){
			if(!alive[i]){continue;}
			byte[] bases=pool.get(i);
			long kmer=0;
			int len=0;
			for(byte b : bases){
				int num=AminoAcid.baseToNumber[b];
				if(num>=0){
					kmer=((kmer<<2)|num)&mask;
					len++;
					if(len>=kLen && selected.contains(kmer)){
						alive[i]=false;
						evicted++;
						break;
					}
				}else{len=0; kmer=0;}
			}
		}
		return evicted;
	}

	/** Converts selected kmer longs into single-kmer FASTA reads. */
	private ArrayList<Read> kmersToReads(LongList kmers, int kLen){
		ArrayList<Read> list=new ArrayList<>(kmers.size);
		for(int i=0; i<kmers.size; i++){
			byte[] bases=longToBytes(kmers.get(i), kLen);
			list.add(new Read(bases, null, "kmer_"+i, i));
		}
		return list;
	}

	/** For k+1 design: each selected (k+1)-mer yields two overlapping k-mers.
	 * Deduplicate and output only the unique k-mers, preserving selection order. */
	private ArrayList<Read> designToUseKmers(LongList designKmers, int kD, int kU){
		final long useMask=~((-1L)<<(2*kU));
		LongHashSet seen=new LongHashSet(designKmers.size*2);
		LongList useKmers=new LongList();
		for(int i=0; i<designKmers.size; i++){
			long dkmer=designKmers.get(i);
			// A (kD)-mer contains (kD-kU+1) overlapping (kU)-mers
			for(int offset=0; offset<=kD-kU; offset++){
				long ukmer=(dkmer>>(2*offset))&useMask;
				if(!seen.contains(ukmer)){
					seen.add(ukmer);
					useKmers.add(ukmer);
				}
			}
		}
		outstream.println("Design "+kD+"-mers: "+designKmers.size+" -> use "+kU+"-mers: "+useKmers.size);
		ArrayList<Read> list=new ArrayList<>(useKmers.size);
		for(int i=0; i<useKmers.size; i++){
			byte[] bases=longToBytes(useKmers.get(i), kU);
			list.add(new Read(bases, null, "kmer_"+i, i));
		}
		return list;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Utility Methods       ----------------*/
	/*--------------------------------------------------------------*/

	private static byte[] longToBytes(long kmer, int kLen){
		byte[] bases=new byte[kLen];
		for(int i=kLen-1; i>=0; i--){
			bases[i]=AminoAcid.numberToBase[(int)(kmer&3)];
			kmer>>>=2;
		}
		return bases;
	}

	/** Co-sort keys by counts, descending, using a simple insertion sort for
	 * small arrays or Arrays.sort on an index array for large ones. */
	private static void coSortDescending(long[] keys, int[] counts, int len){
		if(len<=64){
			for(int i=1; i<len; i++){
				long kk=keys[i]; int cc=counts[i]; int j=i;
				while(j>0 && counts[j-1]<cc){
					keys[j]=keys[j-1]; counts[j]=counts[j-1]; j--;
				}
				keys[j]=kk; counts[j]=cc;
			}
		}else{
			Integer[] idx=new Integer[len];
			for(int i=0; i<len; i++){idx[i]=i;}
			Arrays.sort(idx, (a, b)->counts[b]-counts[a]);
			long[] kk=new long[len]; int[] cc=new int[len];
			for(int i=0; i<len; i++){kk[i]=keys[idx[i]]; cc[i]=counts[idx[i]];}
			System.arraycopy(kk, 0, keys, 0, len);
			System.arraycopy(cc, 0, counts, 0, len);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in;
	private String out;
	private String extra;
	private final FileFormat ffin;
	private final FileFormat ffout;
	private boolean overwrite=true;
	private boolean errorState=false;

	private int k=17;
	private int kDesign=-1;
	private int step=500;
	private int maxKmers=0;
	private float minCovFraction=0.999f;
	private int copies=10;
	private boolean rcomp=false;

	private PrintStream outstream=System.err;
	private static boolean verbose=false;
}
