package gff;

import java.io.*;
import java.util.*;

/**
 * ONE-OFF diagnostic, not a shipped feature. Histograms (called - true) tRNA
 * start/stop offsets using the full ~2M-record leak-fixed flanked training
 * corpus as ground truth (each record's true tRNA span is lflank+1 ..
 * length-rflank, from CutGff's own lflank=/rflank= header fields), instead of
 * the smaller 203-genome benchmark's overlap-matched external references.
 * Much larger N for a statistically robust boundary-offset distribution.
 * Delete once this question is fully absorbed into the boundary-NN plan.
 *
 * Usage: java gff.TrnaBoundaryOffsetHistFlanked <flankedCorpusFasta> <calledGff>
 * calledGff must be callgenes.sh output run directly on flankedCorpusFasta
 * (one record = one "contig"; contig IDs in the GFF are the first
 * whitespace-delimited token of each FASTA header, same truncation CallGenes
 * itself applies -- verified directly against real output before writing this).
 */
public class TrnaBoundaryOffsetHistFlanked {

	public static void main(String[] args) throws IOException {
		String fastaPath=args[0];
		String calledPath=args[1];

		Map<String,int[]> truth=parseTruth(fastaPath); // id -> {lflank, rflank, totalLen}
		System.err.println("Truth records parsed: "+truth.size());

		Map<String,List<int[]>> calledByContig=parseCalledGff(calledPath); // id -> list of {start,end}
		System.err.println("Contigs with >=1 called tRNA: "+calledByContig.size());

		TreeMap<Integer,Integer> startHist=new TreeMap<Integer,Integer>();
		TreeMap<Integer,Integer> stopHist=new TreeMap<Integer,Integer>();
		int matched=0, missed=0;

		for(Map.Entry<String,int[]> e : truth.entrySet()){
			String id=e.getKey();
			int[] t=e.getValue();
			int trueStart=t[0]+1; // lflank+1, 1-based
			int trueEnd=t[2]-t[1]; // totalLen - rflank
			List<int[]> calls=calledByContig.get(id);
			if(calls==null || calls.isEmpty()){missed++; continue;}
			int bestIdx=-1, bestOverlap=-1;
			for(int i=0; i<calls.size(); i++){
				int[] c=calls.get(i);
				int ov=overlap(c[0], c[1], trueStart, trueEnd);
				if(ov>bestOverlap){bestOverlap=ov; bestIdx=i;}
			}
			int[] c=calls.get(bestIdx);
			int startOff=trueStart-c[0];
			int stopOff=c[1]-trueEnd;
			bump(startHist, startOff);
			bump(stopHist, stopOff);
			matched++;
		}

		System.out.println("Matched: "+matched+"  Missed (no call on this record): "+missed);
		printHist("START (5')", startHist, matched);
		printHist("STOP (3')", stopHist, matched);
	}

	private static Map<String,int[]> parseTruth(String path) throws IOException {
		Map<String,int[]> map=new HashMap<String,int[]>();
		BufferedReader br=new BufferedReader(new FileReader(path));
		String line;
		String id=null; int lflank=-1, rflank=-1; StringBuilder seq=new StringBuilder();
		while((line=br.readLine())!=null){
			if(line.isEmpty()){continue;}
			if(line.charAt(0)=='>'){
				if(id!=null){
					map.put(id, new int[]{lflank, rflank, seq.length()});
				}
				String header=line.substring(1);
				id=header.split("\\s+")[0];
				lflank=parseIntField(header, "lflank=");
				rflank=parseIntField(header, "rflank=");
				seq.setLength(0);
			}else{
				seq.append(line.trim());
			}
		}
		if(id!=null){
			map.put(id, new int[]{lflank, rflank, seq.length()});
		}
		br.close();
		return map;
	}

	private static int parseIntField(String header, String key){
		int idx=header.indexOf(key);
		assert(idx>=0) : "Missing "+key+" in header (CutGff always appends lflank=/rflank= to this corpus, TrnaBoundaryOffsetHistFlanked assumes it): "+header;
		int start=idx+key.length();
		int end=start;
		while(end<header.length() && Character.isDigit(header.charAt(end))){end++;}
		return Integer.parseInt(header.substring(start, end));
	}

	private static Map<String,List<int[]>> parseCalledGff(String path) throws IOException {
		Map<String,List<int[]>> map=new HashMap<String,List<int[]>>();
		BufferedReader br=new BufferedReader(new FileReader(path));
		String line;
		while((line=br.readLine())!=null){
			if(line.isEmpty() || line.charAt(0)=='#'){continue;}
			String[] f=line.split("\t");
			if(f.length<8 || !f[2].equals("tRNA")){continue;}
			int gStart=Integer.parseInt(f[3]);
			int gEnd=Integer.parseInt(f[4]);
			List<int[]> list=map.get(f[0]);
			if(list==null){list=new ArrayList<int[]>(); map.put(f[0], list);}
			list.add(new int[]{gStart, gEnd});
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
		System.out.println("=== "+label+" offset histogram (called - true) ===");
		int within1=0, within2=0, within4=0; long sumAbs=0;
		int minOff=Integer.MAX_VALUE, maxOff=Integer.MIN_VALUE;
		for(Map.Entry<Integer,Integer> e : hist.entrySet()){
			int off=e.getKey(), c=e.getValue();
			if(Math.abs(off)<=1){within1+=c;}
			if(Math.abs(off)<=2){within2+=c;}
			if(Math.abs(off)<=4){within4+=c;}
			sumAbs+=(long)Math.abs(off)*c;
			minOff=Math.min(minOff, off);
			maxOff=Math.max(maxOff, off);
		}
		for(int off=Math.max(minOff,-10); off<=Math.min(maxOff,10); off++){
			Integer c=hist.get(off);
			System.out.println((off>=0?"+":"")+off+"\t"+(c==null?0:c));
		}
		long tailBelow=0, tailAbove=0;
		for(Map.Entry<Integer,Integer> e : hist.entrySet()){
			if(e.getKey()<-10){tailBelow+=e.getValue();}
			if(e.getKey()>10){tailAbove+=e.getValue();}
		}
		if(tailBelow>0){System.out.println("< -10\t"+tailBelow);}
		if(tailAbove>0){System.out.println("> +10\t"+tailAbove);}
		System.out.println("range: ["+minOff+", "+maxOff+"]");
		System.out.printf("within +-1bp: %d / %d = %.4f%n", within1, total, within1/(double)total);
		System.out.printf("within +-2bp: %d / %d = %.4f%n", within2, total, within2/(double)total);
		System.out.printf("within +-4bp: %d / %d = %.4f%n", within4, total, within4/(double)total);
		System.out.printf("mean |offset|: %.4f%n", sumAbs/(double)total);
	}
}
