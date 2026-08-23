package gff;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.util.zip.GZIPInputStream;

/**
 * ONE-OFF diagnostic tool, not a shipped feature. Histograms (called_boundary -
 * true_boundary) for tRNA start/stop calls on a benchmark, matching called vs
 * reference GFF loci by overlap, to determine whether the boundary-precision
 * NN's +-2bp candidate search range is sufficient. Delete once this question
 * is answered and absorbed into the tRNA boundary-NN plan.
 *
 * Usage: java gff.TrnaBoundaryOffsetHist <calledGff> <refDir1> [refDir2 ...]
 * calledGff: plain-text GFF3 (e.g. callgenes.sh batch output) with tRNA lines.
 * refDir*: directories containing tid_<N>_*.gff.gz reference annotation files
 *          (searched non-recursively), matched to called contigs by the
 *          tid|<N>|... prefix shared by both called and reference GFFs.
 */
public class TrnaBoundaryOffsetHist {

	public static void main(String[] args) throws IOException {
		String calledPath=args[0];
		Map<String,String> refPaths=new HashMap<String,String>();
		for(int i=1; i<args.length; i++){
			collectRefPaths(new File(args[i]), refPaths);
		}
		System.err.println("Reference genomes found: "+refPaths.size());

		Map<String,List<int[]>> calledByContig=parseGff(new BufferedReader(new InputStreamReader(new FileInputStream(calledPath))));

		Set<String> tids=new TreeSet<String>();
		for(String contig : calledByContig.keySet()){
			String tid=extractTid(contig);
			if(tid!=null){tids.add(tid);}
		}
		tids.addAll(refPaths.keySet());

		TreeMap<Integer,Integer> startHist=new TreeMap<Integer,Integer>();
		TreeMap<Integer,Integer> stopHist=new TreeMap<Integer,Integer>();
		int matched=0, unmatchedRef=0, genomesProcessed=0;

		for(String tid : tids){
			String refPath=refPaths.get(tid);
			if(refPath==null){continue;}
			genomesProcessed++;
			Map<String,List<int[]>> refByContig=parseGff(new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(refPath)))));
			for(Map.Entry<String,List<int[]>> e : refByContig.entrySet()){
				String contig=e.getKey();
				List<int[]> refLoci=e.getValue();
				List<int[]> calledLoci=calledByContig.get(contig);
				if(calledLoci==null){calledLoci=Collections.emptyList();}
				boolean[] claimed=new boolean[calledLoci.size()];
				for(int[] ref : refLoci){
					int bestIdx=-1, bestOverlap=-1;
					for(int i=0; i<calledLoci.size(); i++){
						if(claimed[i]){continue;}
						int[] c=calledLoci.get(i);
						if(c[2]!=ref[2]){continue;}
						int ov=overlap(c[0], c[1], ref[0], ref[1]);
						if(ov>0 && ov>bestOverlap){bestOverlap=ov; bestIdx=i;}
					}
					if(bestIdx>=0){
						claimed[bestIdx]=true;
						int[] c=calledLoci.get(bestIdx);
						int startOff, stopOff;
						if(ref[2]>=0){ // '+' strand: 5'=gStart, 3'=gEnd
							startOff=ref[0]-c[0];
							stopOff=c[1]-ref[1];
						}else{ // '-' strand: 5'=gEnd, 3'=gStart
							startOff=c[1]-ref[1];
							stopOff=ref[0]-c[0];
						}
						bump(startHist, startOff);
						bump(stopHist, stopOff);
						matched++;
					}else{
						unmatchedRef++;
					}
				}
			}
		}

		System.out.println("Genomes with a matched reference: "+genomesProcessed+" / "+tids.size());
		System.out.println("Matched loci (overlap-matched, offset computed): "+matched);
		System.out.println("Unmatched reference loci (no overlapping same-strand call, i.e. false negatives): "+unmatchedRef);
		printHist("START (5')", startHist, matched);
		printHist("STOP (3')", stopHist, matched);
	}

	private static void collectRefPaths(File dir, Map<String,String> out){
		File[] files=dir.listFiles();
		if(files==null){return;}
		Pattern p=Pattern.compile("^tid_(\\d+)_.*\\.gff\\.gz$");
		for(File f : files){
			Matcher m=p.matcher(f.getName());
			if(m.matches()){
				out.put(m.group(1), f.getAbsolutePath());
			}
		}
	}

	private static String extractTid(String contig){
		// format: tid|100177|NZ_CP124050.1
		String[] split=contig.split("\\|");
		return split.length>=2 ? split[1] : null;
	}

	/** Parses tRNA lines (column 3) from a GFF3 stream, grouped by contig (column 1).
	 * Each locus is stored as int[]{gStart, gEnd, strand} (strand: +1 for '+', -1 for '-'). */
	private static Map<String,List<int[]>> parseGff(BufferedReader br) throws IOException {
		Map<String,List<int[]>> map=new HashMap<String,List<int[]>>();
		String line;
		while((line=br.readLine())!=null){
			if(line.isEmpty() || line.charAt(0)=='#'){continue;}
			String[] f=line.split("\t");
			if(f.length<8 || !f[2].equals("tRNA")){continue;}
			int strand=f[6].equals("-") ? -1 : 1;
			int gStart=Integer.parseInt(f[3]);
			int gEnd=Integer.parseInt(f[4]);
			List<int[]> list=map.get(f[0]);
			if(list==null){list=new ArrayList<int[]>(); map.put(f[0], list);}
			list.add(new int[]{gStart, gEnd, strand});
		}
		br.close();
		return map;
	}

	private static int overlap(int a1, int a2, int b1, int b2){
		return Math.max(0, Math.min(a2, b2)-Math.max(a1, b1)+1);
	}

	private static void bump(TreeMap<Integer,Integer> hist, int offset){
		Integer v=hist.get(offset);
		hist.put(offset, v==null ? 1 : v+1);
	}

	private static void printHist(String label, TreeMap<Integer,Integer> hist, int total){
		System.out.println();
		System.out.println("=== "+label+" offset histogram (called - true, signed, read-orientation) ===");
		int within1=0, within2=0, sumAbs=0, minOff=Integer.MAX_VALUE, maxOff=Integer.MIN_VALUE;
		for(Map.Entry<Integer,Integer> e : hist.entrySet()){
			int off=e.getKey(), c=e.getValue();
			if(Math.abs(off)<=1){within1+=c;}
			if(Math.abs(off)<=2){within2+=c;}
			sumAbs+=Math.abs(off)*c;
			minOff=Math.min(minOff, off);
			maxOff=Math.max(maxOff, off);
		}
		for(int off=Math.max(minOff,-15); off<=Math.min(maxOff,15); off++){
			Integer c=hist.get(off);
			int count=(c==null) ? 0 : c;
			if(count>0 || (off>=-10 && off<=10)){
				System.out.println((off>=0?"+":"")+off+"\t"+count);
			}
		}
		int tailBelow=0, tailAbove=0;
		for(Map.Entry<Integer,Integer> e : hist.entrySet()){
			if(e.getKey()<-15){tailBelow+=e.getValue();}
			if(e.getKey()>15){tailAbove+=e.getValue();}
		}
		if(tailBelow>0){System.out.println("< -15\t"+tailBelow);}
		if(tailAbove>0){System.out.println("> +15\t"+tailAbove);}
		System.out.println("range: ["+minOff+", "+maxOff+"]");
		System.out.printf("within +-1bp: %d / %d = %.4f%n", within1, total, within1/(double)total);
		System.out.printf("within +-2bp: %d / %d = %.4f%n", within2, total, within2/(double)total);
		System.out.printf("mean |offset|: %.4f%n", sumAbs/(double)total);
	}
}
