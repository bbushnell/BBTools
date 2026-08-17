package prot;

import java.util.ArrayList;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.IntHashMap;
import map.IntHashSet;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;
import tracker.KmerTracker;

/**
 * Computes per-organism dinucleotide-composition features (HH and CAGA from
 * {@link tracker.KmerTracker}) from a contig FASTA whose headers embed the taxon id
 * ({@code ..._tid_<int>}), then MIN-MAX scales each feature so that 1 = the maximum
 * and 0 = the minimum observed across the TRAINING organisms only (val organisms are
 * scaled by the same train-derived range, clamped to [0,1]). The train/val split is
 * reproduced exactly from {@link MagQCVectorMaker}'s procedure: usable = tids present
 * in the per-contig cache with at least one contig row of length &gt;= minlen (mirrors
 * {@code loadCache}'s per-row minlen filter) AND present in the sizemap file (mirrors
 * {@code genomeSize}, loaded from a file separate from the cache) -- sorted, shuffled
 * by {@code new Random(seed)}, first {@code round(valfrac*n)} held out -- so the
 * scaling never leaks validation organisms.
 *
 * <p>{@code sizemap=} is REQUIRED (not optional): {@link MagQCVectorMaker} itself
 * requires it, and a tid present in the cache but absent from the sizemap is exactly
 * the divergence that let validation organisms bleed into the training scaling range
 * (fixed 2026-08-16 -- the original cache-only usable set silently diverged from
 * MagQCVectorMaker's real one whenever the two files' tid coverage differed, and
 * {@code Collections.shuffle} is size/order-sensitive, so ANY membership difference
 * produces a completely different train/val permutation, not just a locally-different
 * one).</p>
 *
 * <p>Output: {@code tid<TAB>HH<TAB>CAGA}, one usable organism per line, consumed by
 * {@code MagQCVectorMaker subnet=ncrna snhhcaga=t kmerfile=<this>}.</p>
 *
 * <p>Usage: {@code java prot.GenomeDimerFeatures in=renamed.fa cache=percontig_cache.tsv
 * sizemap=sizemap.tsv out=kmerfeat.tsv [seed=1] [valfrac=0.10] [minlen=0]}</p>
 *
 * @author Eru
 */
public class GenomeDimerFeatures {

	public static void main(String[] args){
		String in=null, cache=null, out=null, sizemap=null;
		long seed=1;
		double valfrac=0.10;
		int minlen=0;
		for(String a : args){
			int e=a.indexOf('=');
			if(e<0){continue;}
			String k=a.substring(0, e).toLowerCase(), v=a.substring(e+1);
			if(k.equals("in")){in=v;}
			else if(k.equals("cache")){cache=v;}
			else if(k.equals("sizemap")){sizemap=v;}
			else if(k.equals("out")){out=v;}
			else if(k.equals("seed")){seed=Long.parseLong(v);}
			else if(k.equals("valfrac")){valfrac=Double.parseDouble(v);}
			else if(k.equals("minlen")){minlen=Integer.parseInt(v);}
		}
		if(in==null || cache==null || out==null || sizemap==null){
			throw new RuntimeException("Required: in=<contig.fa> cache=<percontig_cache.tsv> "
				+"sizemap=<tid bp, same file MagQCVectorMaker uses> out=<tid HH CAGA> [minlen=0]");
		}

		//1) usable tids = distinct cache tids (col f1) with >=1 contig row of length (col f3)
		//   >= minlen, reproducing MagQCVectorMaker's loadCache() per-row minlen filter exactly.
		final IntHashSet cacheUsable=new IntHashSet();
		final ByteFile bf=ByteFile.makeByteFile(cache, true);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			final int t1=indexOf(line, (byte)'\t', 0);
			if(t1<0){continue;}
			final int t2=indexOf(line, (byte)'\t', t1+1);
			if(t2<0){continue;}
			final int t3=indexOf(line, (byte)'\t', t2+1);
			if(t3<0){continue;}
			final int t4=indexOf(line, (byte)'\t', t3+1);
			if(t4<0){continue;}
			final int length=Integer.parseInt(new String(line, t3+1, t4-t3-1));
			if(length<minlen){continue;}
			cacheUsable.add(Integer.parseInt(new String(line, t1+1, t2-t1-1)));
		}
		bf.close();

		//2) sizemap tids (col 0), reproducing MagQCVectorMaker's genomeSize key set exactly.
		final IntHashSet sizemapTids=new IntHashSet();
		final ByteFile sbf=ByteFile.makeByteFile(sizemap, true);
		for(byte[] line=sbf.nextLine(); line!=null; line=sbf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			final int t1=indexOf(line, (byte)'\t', 0);
			if(t1<0){continue;}
			sizemapTids.add(Integer.parseInt(new String(line, 0, t1)));
		}
		sbf.close();

		//3) usable = cache-usable AND sizemap-present -- MagQCVectorMaker's real usable set.
		final IntHashSet usableSet=new IntHashSet();
		{
			final int[] cacheArray=cacheUsable.toArray();
			for(int tid : cacheArray){
				if(sizemapTids.contains(tid)){usableSet.add(tid);}
			}
		}

		//2) reproduce the exact train/val split (sorted, then Fisher-Yates shuffled by Random(seed)
		//   using the identical algorithm java.util.Collections.shuffle uses on a RandomAccess list --
		//   for(i=size; i>1; i--) swap(i-1, rnd.nextInt(i)) -- so the resulting permutation is
		//   byte-identical to the original boxed-list version given the same seed and starting order.
		final int[] usable=usableSet.toArray();
		java.util.Arrays.sort(usable);
		{
			final Random rnd=new Random(seed);
			for(int i=usable.length; i>1; i--){
				final int j=rnd.nextInt(i);
				final int tmp=usable[i-1]; usable[i-1]=usable[j]; usable[j]=tmp;
			}
		}
		final int nVal=(int)Math.round(usable.length*valfrac);
		System.err.println("usable orgs="+usable.length+", train="+(usable.length-nVal)+", val="+nVal);

		//3) accumulate dimer counts per tid from the contig FASTA (only for usable tids).
		//   tidToIndex maps tid -> position in the parallel `trackers` list (IntHashMap+parallel
		//   array instead of a boxed HashMap<Integer,KmerTracker>).
		final IntHashMap tidToIndex=new IntHashMap();
		final ArrayList<KmerTracker> trackers=new ArrayList<KmerTracker>();
		final IntList presentTids=new IntList();
		final FileFormat ff=FileFormat.testInput(in, FileFormat.FASTA, null, true, true);
		final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, false, ff, null);
		cris.start();
		long contigs=0;
		for(ListNum<Read> ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()){
			for(Read r : ln){
				if(r.bases==null || r.bases.length<2){continue;}
				final int tid=parseTid(r.id);
				if(tid<0 || !usableSet.contains(tid)){continue;}
				int idx=tidToIndex.get(tid);
				if(idx<0){
					idx=trackers.size();
					tidToIndex.put(tid, idx);
					trackers.add(new KmerTracker(2));
					presentTids.add(tid);
				}
				trackers.get(idx).add(r.bases);
				contigs++;
			}
			cris.returnList(ln);
		}
		cris.close();
		System.err.println("accumulated "+contigs+" contigs across "+trackers.size()+" orgs");

		//4) per-org HH, CAGA, indexed identically to `trackers`/`tidToIndex`.
		final float[] hhArr=new float[trackers.size()], caArr=new float[trackers.size()];
		for(int i=0; i<trackers.size(); i++){
			final long[] c=trackers.get(i).counts;
			hhArr[i]=KmerTracker.HH(c);
			caArr[i]=KmerTracker.CAGA(c);
		}

		//5) min/max over TRAIN orgs only (usable[nVal..end)).
		float hhMin=Float.MAX_VALUE, hhMax=-Float.MAX_VALUE, caMin=Float.MAX_VALUE, caMax=-Float.MAX_VALUE;
		for(int ui=nVal; ui<usable.length; ui++){
			final int idx=tidToIndex.get(usable[ui]);
			if(idx<0){continue;}
			final float hh=hhArr[idx], ca=caArr[idx];
			if(Float.isNaN(hh) || Float.isNaN(ca)){continue;}
			hhMin=Math.min(hhMin, hh); hhMax=Math.max(hhMax, hh);
			caMin=Math.min(caMin, ca); caMax=Math.max(caMax, ca);
		}
		final float hhRange=Math.max(1e-9f, hhMax-hhMin), caRange=Math.max(1e-9f, caMax-caMin);
		System.err.println("train HH ["+hhMin+","+hhMax+"]  CAGA ["+caMin+","+caMax+"]");

		//6) scale every usable org by the train range, clamp [0,1], write.
		presentTids.sort();
		final FileFormat outFF=FileFormat.testOutput(out, FileFormat.TEXT, null, false, true, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(outFF);
		bsw.start();
		final ByteBuilder bb=new ByteBuilder();
		int written=0;
		for(int pi=0; pi<presentTids.size; pi++){
			final int tid=presentTids.get(pi);
			final int idx=tidToIndex.get(tid);
			final float hh0=hhArr[idx], ca0=caArr[idx];
			if(Float.isNaN(hh0) || Float.isNaN(ca0)){continue;}
			final float hh=clamp01((hh0-hhMin)/hhRange);
			final float ca=clamp01((ca0-caMin)/caRange);
			//Float.toString() (via append(String), not a fixed-decimals numeric append) preserves
			//the original "tid+"\t"+hh+"\t"+ca+"\n"" string-concatenation formatting byte-for-byte
			//-- ByteBuilder has no shortest-round-trip float append, only fixed-decimals ones.
			bb.clear();
			bb.append(tid).tab().append(Float.toString(hh)).tab().append(Float.toString(ca)).nl();
			bsw.print(bb.toBytes());
			written++;
		}
		bsw.poisonAndWait();
		System.err.println("wrote "+written+" orgs to "+out);
	}

	/** Parses the integer following "_tid_" in a header, or -1 if absent. */
	static int parseTid(String header){
		if(header==null){return -1;}
		final int i=header.indexOf("_tid_");
		if(i<0){return -1;}
		int j=i+5, k=j;
		while(k<header.length() && Character.isDigit(header.charAt(k))){k++;}
		return (k>j ? Integer.parseInt(header.substring(j, k)) : -1);
	}

	private static float clamp01(float v){return v<0 ? 0 : (v>1 ? 1 : v);}

	private static int indexOf(byte[] a, byte b, int from){
		for(int i=from; i<a.length; i++){if(a[i]==b){return i;}}
		return -1;
	}
}
