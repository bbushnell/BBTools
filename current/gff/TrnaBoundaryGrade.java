package gff;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.util.zip.GZIPInputStream;

/**
 * ONE-OFF grading tool, not a shipped feature. Reports both-exact/5'-exact/
 * 3'-exact tRNA boundary accuracy on a benchmark -- the standard metric names
 * this project already uses (see the trna skill's eval_anticodon3.sh figures,
 * e.g. "81.0% both-exact / 91.9% 5'-exact"), computed the fast way (batch
 * callgenes + overlap matching, no per-genome JVM). Built for Part 3's
 * with-vs-without NN-refiner grading decision. Shares its overlap-matching
 * method with gff.TrnaBoundaryOffsetHist (same house rule: written once, not
 * duplicated, so both tools' notion of "matched locus" cannot drift).
 *
 * Usage: java gff.TrnaBoundaryGrade <calledGff> <refDir1> [refDir2 ...]
 */
public class TrnaBoundaryGrade {

	public static void main(String[] args) throws IOException {
		String calledPath=args[0];
		Map<String,String> refPaths=new HashMap<String,String>();
		for(int i=1; i<args.length; i++){
			collectRefPaths(new File(args[i]), refPaths);
		}

		Map<String,List<int[]>> calledByContig=parseGff(new BufferedReader(new InputStreamReader(new FileInputStream(calledPath))));

		Set<String> tids=new TreeSet<String>();
		for(String contig : calledByContig.keySet()){
			String tid=extractTid(contig);
			if(tid!=null){tids.add(tid);}
		}
		tids.addAll(refPaths.keySet());

		int matched=0, unmatchedRef=0, bothExact=0, startExact=0, stopExact=0, genomesProcessed=0;

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
						final boolean sExact, eExact;
						if(ref[2]>=0){sExact=(c[0]==ref[0]); eExact=(c[1]==ref[1]);}
						else{sExact=(c[1]==ref[1]); eExact=(c[0]==ref[0]);}
						if(sExact){startExact++;}
						if(eExact){stopExact++;}
						if(sExact && eExact){bothExact++;}
						matched++;
					}else{
						unmatchedRef++;
					}
				}
			}
		}

		final int totalRef=matched+unmatchedRef;
		System.out.println("Genomes with a matched reference: "+genomesProcessed+" / "+tids.size());
		System.out.println("Total reference tRNA loci: "+totalRef);
		System.out.println("Matched (overlap-matched, same-strand) loci: "+matched+" ("
			+String.format("%.4f", matched/(double)totalRef)+" recall)");
		System.out.println("Unmatched reference loci (false negatives): "+unmatchedRef);
		System.out.println();
		System.out.printf("both-exact: %d / %d = %.4f%n", bothExact, matched, bothExact/(double)matched);
		System.out.printf("5'-exact (start): %d / %d = %.4f%n", startExact, matched, startExact/(double)matched);
		System.out.printf("3'-exact (stop): %d / %d = %.4f%n", stopExact, matched, stopExact/(double)matched);
		System.out.println();
		System.out.printf("both-exact (of total ref, i.e. accounts for FNs too): %.4f%n", bothExact/(double)totalRef);
	}

	private static void collectRefPaths(File dir, Map<String,String> out){
		File[] files=dir.listFiles();
		if(files==null){return;}
		Pattern p=Pattern.compile("^tid_(\\d+)_.*\\.gff\\.gz$");
		for(File f : files){
			Matcher m=p.matcher(f.getName());
			if(m.matches()){out.put(m.group(1), f.getAbsolutePath());}
		}
	}

	private static String extractTid(String contig){
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
}
