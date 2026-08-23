package gff;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.util.zip.GZIPInputStream;

/**
 * ONE-OFF diagnostic, not shipped. For Brian's failure-mode analysis (Aug 22 2026):
 * joins the WITHOUT-refiner and WITH-refiner called GFFs against the reference,
 * per locus, and classifies each matched locus as FIXED / REGRESSED / STILL_WRONG /
 * UNCHANGED_CORRECT / UNCHANGED_WRONG, independently for start and stop. Reuses
 * the same overlap-matching method as gff.TrnaBoundaryGrade/TrnaBoundaryOffsetHist
 * (can't drift -- one matching definition for the whole project).
 *
 * Usage: java gff.TrnaBoundaryDiff <beforeGff> <afterGff> <refDir1> [refDir2 ...]
 *
 * Output (stdout, TSV): contig, trueStart, trueStop, strand, beforeStart, beforeStop,
 * afterStart, afterStop, startClass, stopClass -- one row per reference locus matched
 * in BOTH called gffs (loci missed entirely by either run are reported separately,
 * not classified).
 */
public class TrnaBoundaryDiff {

	public static void main(String[] args) throws IOException {
		String beforePath=args[0], afterPath=args[1];
		Map<String,String> refPaths=new HashMap<String,String>();
		for(int i=2; i<args.length; i++){collectRefPaths(new File(args[i]), refPaths);}

		Map<String,List<int[]>> beforeByContig=parseGff(new BufferedReader(new InputStreamReader(new FileInputStream(beforePath))));
		Map<String,List<int[]>> afterByContig=parseGff(new BufferedReader(new InputStreamReader(new FileInputStream(afterPath))));

		Set<String> tids=new TreeSet<String>();
		for(String contig : beforeByContig.keySet()){final String t=extractTid(contig); if(t!=null){tids.add(t);}}
		for(String contig : afterByContig.keySet()){final String t=extractTid(contig); if(t!=null){tids.add(t);}}
		tids.addAll(refPaths.keySet());

		int bothMatched=0, missedEither=0;
		final Map<String,Integer> counts=new TreeMap<String,Integer>();

		System.out.println("contig\ttrueStart\ttrueStop\tstrand\tbeforeStart\tbeforeStop\tafterStart\tafterStop"
			+"\tbeforeStartOff\tbeforeStopOff\tafterStartOff\tafterStopOff\tstartClass\tstopClass");

		for(String tid : tids){
			String refPath=refPaths.get(tid);
			if(refPath==null){continue;}
			Map<String,List<int[]>> refByContig=parseGff(new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(refPath)))));
			for(Map.Entry<String,List<int[]>> e : refByContig.entrySet()){
				final String contig=e.getKey();
				final List<int[]> refLoci=e.getValue();
				final List<int[]> beforeLoci=beforeByContig.getOrDefault(contig, Collections.emptyList());
				final List<int[]> afterLoci=afterByContig.getOrDefault(contig, Collections.emptyList());
				final boolean[] beforeClaimed=new boolean[beforeLoci.size()];
				final boolean[] afterClaimed=new boolean[afterLoci.size()];
				for(int[] ref : refLoci){
					final int[] b=bestMatch(beforeLoci, beforeClaimed, ref);
					final int[] a=bestMatch(afterLoci, afterClaimed, ref);
					if(b==null || a==null){missedEither++; continue;}
					bothMatched++;
					final boolean plus=(ref[2]>=0);
					final int trueStart=(plus ? ref[0] : ref[1]), trueStop=(plus ? ref[1] : ref[0]);
					final int bStart=(plus ? b[0] : b[1]), bStop=(plus ? b[1] : b[0]);
					final int aStart=(plus ? a[0] : a[1]), aStop=(plus ? a[1] : a[0]);
					final boolean bExactStart=(bStart==trueStart), aExactStart=(aStart==trueStart);
					final boolean bExactStop=(bStop==trueStop), aExactStop=(aStop==trueStop);
					final String startClass=classify(bExactStart, aExactStart);
					final String stopClass=classify(bExactStop, aExactStop);
					counts.merge("start:"+startClass, 1, Integer::sum);
					counts.merge("stop:"+stopClass, 1, Integer::sum);
					//STRAND-AWARE signed offset (positive = called extends OUTWARD past true, same
					//formula as gff.TrnaBoundaryOffsetHist -- computed from the RAW genomic ref/b/a
					//arrays, not the reoriented trueStart/bStart above: those are read-frame
					//coordinates whose NUMERIC DIRECTION already flips between strands (for minus
					//strand, "stop" = the smaller genomic coordinate, so a raw bStop-trueStop
					//subtraction on the reoriented values gets the sign BACKWARDS for minus-strand
					//loci -- caught before it corrupted the failure-mode analysis, see plan doc).
					final int bStartOff, bStopOff, aStartOff, aStopOff;
					if(plus){
						bStartOff=ref[0]-b[0]; bStopOff=b[1]-ref[1];
						aStartOff=ref[0]-a[0]; aStopOff=a[1]-ref[1];
					}else{
						bStartOff=b[1]-ref[1]; bStopOff=ref[0]-b[0];
						aStartOff=a[1]-ref[1]; aStopOff=ref[0]-a[0];
					}
					System.out.println(contig+"\t"+trueStart+"\t"+trueStop+"\t"+(plus?"+":"-")+"\t"
						+bStart+"\t"+bStop+"\t"+aStart+"\t"+aStop+"\t"
						+bStartOff+"\t"+bStopOff+"\t"+aStartOff+"\t"+aStopOff+"\t"+startClass+"\t"+stopClass);
				}
			}
		}

		System.err.println("Matched in both runs: "+bothMatched+"  (missed by at least one run: "+missedEither+")");
		for(Map.Entry<String,Integer> ce : counts.entrySet()){
			System.err.println(ce.getKey()+"\t"+ce.getValue());
		}
	}

	private static String classify(boolean before, boolean after){
		if(before && after){return "UNCHANGED_CORRECT";}
		if(!before && after){return "FIXED";}
		if(before && !after){return "REGRESSED";}
		return "STILL_WRONG";
	}

	private static int[] bestMatch(List<int[]> loci, boolean[] claimed, int[] ref){
		int bestIdx=-1, bestOverlap=-1;
		for(int i=0; i<loci.size(); i++){
			if(claimed[i]){continue;}
			final int[] c=loci.get(i);
			if(c[2]!=ref[2]){continue;}
			final int ov=overlap(c[0], c[1], ref[0], ref[1]);
			if(ov>0 && ov>bestOverlap){bestOverlap=ov; bestIdx=i;}
		}
		if(bestIdx<0){return null;}
		claimed[bestIdx]=true;
		return loci.get(bestIdx);
	}

	private static void collectRefPaths(File dir, Map<String,String> out){
		final File[] files=dir.listFiles();
		if(files==null){return;}
		final Pattern p=Pattern.compile("^tid_(\\d+)_.*\\.gff\\.gz$");
		for(File f : files){
			final Matcher m=p.matcher(f.getName());
			if(m.matches()){out.put(m.group(1), f.getAbsolutePath());}
		}
	}

	private static String extractTid(String contig){
		final String[] split=contig.split("\\|");
		return split.length>=2 ? split[1] : null;
	}

	private static Map<String,List<int[]>> parseGff(BufferedReader br) throws IOException {
		final Map<String,List<int[]>> map=new HashMap<String,List<int[]>>();
		String line;
		while((line=br.readLine())!=null){
			if(line.isEmpty() || line.charAt(0)=='#'){continue;}
			final String[] f=line.split("\t");
			if(f.length<8 || !f[2].equals("tRNA")){continue;}
			final int strand=f[6].equals("-") ? -1 : 1;
			final int gStart=Integer.parseInt(f[3]);
			final int gEnd=Integer.parseInt(f[4]);
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
