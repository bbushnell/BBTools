package prot;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.FileFormat;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;
import tracker.KmerTracker;

/**
 * Computes per-organism dinucleotide-composition features (HH and CAGA from
 * {@link tracker.KmerTracker}) from a contig FASTA whose headers embed the taxon id
 * ({@code ..._tid_<int>}), then MIN-MAX scales each feature so that 1 = the maximum
 * and 0 = the minimum observed across the TRAINING organisms only (val organisms are
 * scaled by the same train-derived range, clamped to [0,1]). The train/val split is
 * reproduced exactly from {@link MagQCVectorMaker}'s procedure (usable = the distinct
 * tids in the per-contig cache, sorted, shuffled by {@code new Random(seed)}, first
 * {@code round(valfrac*n)} held out) so the scaling never leaks validation organisms.
 *
 * <p>Output: {@code tid<TAB>HH<TAB>CAGA}, one usable organism per line, consumed by
 * {@code MagQCVectorMaker subnet=ncrna snhhcaga=t kmerfile=<this>}.</p>
 *
 * <p>Usage: {@code java prot.GenomeDimerFeatures in=renamed.fa cache=percontig_cache.tsv
 * out=kmerfeat.tsv [seed=1] [valfrac=0.10]}</p>
 *
 * @author Eru
 */
public class GenomeDimerFeatures {

	public static void main(String[] args){
		String in=null, cache=null, out=null;
		long seed=1;
		double valfrac=0.10;
		for(String a : args){
			int e=a.indexOf('=');
			if(e<0){continue;}
			String k=a.substring(0, e).toLowerCase(), v=a.substring(e+1);
			if(k.equals("in")){in=v;}
			else if(k.equals("cache")){cache=v;}
			else if(k.equals("out")){out=v;}
			else if(k.equals("seed")){seed=Long.parseLong(v);}
			else if(k.equals("valfrac")){valfrac=Double.parseDouble(v);}
		}
		if(in==null || cache==null || out==null){
			throw new RuntimeException("Required: in=<contig.fa> cache=<percontig_cache.tsv> out=<tid HH CAGA>");
		}

		//1) usable tids = distinct tids in the cache (col f1), reproducing MagQCVectorMaker's usable set.
		final HashSet<Integer> usableSet=new HashSet<Integer>();
		final ByteFile bf=ByteFile.makeByteFile(cache, true);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			int t1=indexOf(line, (byte)'\t', 0);
			if(t1<0){continue;}
			int t2=indexOf(line, (byte)'\t', t1+1);
			if(t2<0){continue;}
			usableSet.add(Integer.parseInt(new String(line, t1+1, t2-t1-1)));
		}
		bf.close();

		//2) reproduce the exact train/val split (sorted, shuffled by Random(seed), first round(valfrac*n)=val).
		final ArrayList<Integer> usable=new ArrayList<Integer>(usableSet);
		Collections.sort(usable);
		Collections.shuffle(usable, new Random(seed));
		final int nVal=(int)Math.round(usable.size()*valfrac);
		final HashSet<Integer> trainTids=new HashSet<Integer>(usable.subList(nVal, usable.size()));
		System.err.println("usable orgs="+usable.size()+", train="+trainTids.size()+", val="+nVal);

		//3) accumulate dimer counts per tid from the contig FASTA (only for usable tids).
		final HashMap<Integer,KmerTracker> byTid=new HashMap<Integer,KmerTracker>();
		final FileFormat ff=FileFormat.testInput(in, FileFormat.FASTA, null, true, true);
		final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, false, ff, null);
		cris.start();
		long contigs=0;
		for(ListNum<Read> ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()){
			for(Read r : ln){
				if(r.bases==null || r.bases.length<2){continue;}
				final int tid=parseTid(r.id);
				if(tid<0 || !usableSet.contains(tid)){continue;}
				KmerTracker kt=byTid.get(tid);
				if(kt==null){byTid.put(tid, kt=new KmerTracker(2));}
				kt.add(r.bases);
				contigs++;
			}
			cris.returnList(ln);
		}
		cris.close();
		System.err.println("accumulated "+contigs+" contigs across "+byTid.size()+" orgs");

		//4) per-org HH, CAGA.
		final HashMap<Integer,float[]> feat=new HashMap<Integer,float[]>();
		for(Integer tid : byTid.keySet()){
			final long[] c=byTid.get(tid).counts;
			feat.put(tid, new float[]{KmerTracker.HH(c), KmerTracker.CAGA(c)});
		}

		//5) min/max over TRAIN orgs only.
		float hhMin=Float.MAX_VALUE, hhMax=-Float.MAX_VALUE, caMin=Float.MAX_VALUE, caMax=-Float.MAX_VALUE;
		for(Integer tid : trainTids){
			final float[] f=feat.get(tid);
			if(f==null || Float.isNaN(f[0]) || Float.isNaN(f[1])){continue;}
			hhMin=Math.min(hhMin, f[0]); hhMax=Math.max(hhMax, f[0]);
			caMin=Math.min(caMin, f[1]); caMax=Math.max(caMax, f[1]);
		}
		final float hhRange=Math.max(1e-9f, hhMax-hhMin), caRange=Math.max(1e-9f, caMax-caMin);
		System.err.println("train HH ["+hhMin+","+hhMax+"]  CAGA ["+caMin+","+caMax+"]");

		//6) scale every usable org by the train range, clamp [0,1], write.
		try{
			final BufferedWriter bw=new BufferedWriter(new FileWriter(out));
			final ArrayList<Integer> sorted=new ArrayList<Integer>(feat.keySet());
			Collections.sort(sorted);
			int written=0;
			for(Integer tid : sorted){
				final float[] f=feat.get(tid);
				if(f==null || Float.isNaN(f[0]) || Float.isNaN(f[1])){continue;}
				final float hh=clamp01((f[0]-hhMin)/hhRange);
				final float ca=clamp01((f[1]-caMin)/caRange);
				bw.write(tid+"\t"+hh+"\t"+ca+"\n");
				written++;
			}
			bw.close();
			System.err.println("wrote "+written+" orgs to "+out);
		}catch(Exception e){throw new RuntimeException(e);}
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
