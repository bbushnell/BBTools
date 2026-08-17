package gff;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import fileIO.FileFormat;

/**
 * Extracts and classifies the FALSE-NEGATIVE tRNA set from a benchmark eval, for the
 * recall-gap investigation (Neptune's FN characterization, Step 1-2). A reference tRNA is
 * an FN when NO called tRNA overlaps it by reciprocal &gt;50% (2*overlap &gt; max(refLen,qryLen))
 * on the SAME seqid — the IDENTICAL predicate + loader (GffLine.loadGffFile) used by
 * gff.CompareGff's overlap metrics, applied ref-centrically so we recover the distinct
 * unmatched references (not the query-centric refCount-tp count). FN count per genome
 * cross-checks against CompareGff's reported tRNA 'fn' (they agree when matching is ~1:1).
 *
 * Emits outPrefix_fn.tsv (one row per FN ref tRNA) and prints a classification summary
 * (by amino acid, by length bin, by domain, and count of &gt;100nt likely-intron-bearing).
 *
 * Usage: TrnaFN outPrefix  domain benchList calledDir  [domain benchList calledDir ...]
 *   benchList  = list of genome .fna.gz paths (ref GFF derived: .fna.gz -> .gff.gz alongside)
 *   calledDir  = dir of callgenes output GFFs, one per genome as <name>.gff
 *
 * @author Noire
 * @date August 16, 2026
 */
public class TrnaFN {

	static final Pattern REF_NOTE=Pattern.compile("Note=tRNA-(\\w+)\\((\\w+)\\)");
	static final Pattern REF_PRODUCT=Pattern.compile("product=tRNA-(\\w+)");

	public static void main(String[] args){
		if(args.length<4 || (args.length-1)%3!=0){
			System.err.println("Usage: TrnaFN outPrefix  domain benchList calledDir  [domain benchList calledDir ...]");
			return;
		}
		final String outPrefix=args[0];

		final StringBuilder tsv=new StringBuilder();
		tsv.append("domain\tgenome\tseqid\tstart\tstop\tstrand\tamino\tanticodon\tlength\tintron_likely\n");

		// aggregation
		final TreeMap<String,int[]> byAmino=new TreeMap<>();      // amino -> {fn}
		final TreeMap<String,int[]> byDomain=new TreeMap<>();     // domain -> {refTotal, fn, fnIntron, callTotal, fp}
		final int[] lenBins=new int[6];                            // <70,70-79,80-89,90-99,100-119,>=120
		int totalRef=0, totalFN=0, totalIntronFN=0, missingGenomes=0, processedGenomes=0;

		for(int t=1; t+2<args.length; t+=3){
			final String domain=args[t], benchList=args[t+1], calledDir=args[t+2];
			byDomain.putIfAbsent(domain, new int[5]);
			//COMBINED mode: if calledDir is a single GFF file (one-JVM batched callgenes output over ALL
			//genomes at once), load it once keyed by RAW seqid; else the original per-genome <name>.gff dir.
			//seqids are globally unique, so a genome's calls are exactly the combined calls on its ref
			//seqids -- robust, no tid parsing.  (Recall is exact; a spurious call on a ref-LESS contig is
			//not attributed to a genome here, a negligible precision effect.)
			final boolean combined=new File(calledDir).isFile();
			final HashMap<String,ArrayList<GffLine>> combinedBySeqid=(combined ? loadCombinedBySeqid(calledDir) : null);
			final ArrayList<String> genomes=readList(benchList);
			for(String fna : genomes){
				final String name=stripFna(new File(fna).getName());
				final String refGff=fna.replace(".fna.gz", ".gff.gz");
				if(!new File(refGff).exists()){missingGenomes++; continue;}
				final ArrayList<GffLine> refs=loadTrna(refGff);
				final ArrayList<GffLine> calls;
				if(combined){
					final java.util.HashSet<String> mySeqids=new java.util.HashSet<>();
					for(GffLine r : refs){mySeqids.add(r.seqid);}
					calls=new ArrayList<>();
					for(String s : mySeqids){
						final ArrayList<GffLine> c=combinedBySeqid.get(s);
						if(c!=null){calls.addAll(c);}
					}
				}else{
					final String calledGff=calledDir+"/"+name+".gff";
					if(!new File(calledGff).exists()){missingGenomes++; continue;}
					calls=loadTrna(calledGff);
				}
				processedGenomes++;
				// index calls by seqid, sorted by start
				final HashMap<String,ArrayList<GffLine>> callsBySeqid=new HashMap<>();
				for(GffLine q : calls){
					callsBySeqid.computeIfAbsent(q.seqid, k->new ArrayList<>()).add(q);
				}
				for(ArrayList<GffLine> list : callsBySeqid.values()){Collections.sort(list, START_COMP);}

				// query-centric precision pass: a CALL is FP if it overlaps no ref (same predicate)
				final HashMap<String,ArrayList<GffLine>> refBySeqid=new HashMap<>();
				for(GffLine r : refs){refBySeqid.computeIfAbsent(r.seqid, k->new ArrayList<>()).add(r);}
				for(ArrayList<GffLine> list : refBySeqid.values()){Collections.sort(list, START_COMP);}
				byDomain.get(domain)[3]+=calls.size();
				for(GffLine q : calls){
					if(!overlapsAny(q, refBySeqid.get(q.seqid))){byDomain.get(domain)[4]++;}
				}

				for(GffLine r : refs){
					totalRef++;
					byDomain.get(domain)[0]++;
					if(overlapsAny(r, callsBySeqid.get(r.seqid))){continue;}
					// FALSE NEGATIVE
					final String[] aa=refAmino(r.attributes);
					final String amino=aa[0], anticodon=(aa[1]==null ? "NA" : aa[1]);
					final int length=r.stop-r.start+1;
					final boolean intron=length>100;
					final char strand=(r.strand==GffLine.MINUS ? '-' : '+');
					tsv.append(domain).append('\t').append(name).append('\t').append(r.seqid).append('\t')
						.append(r.start).append('\t').append(r.stop).append('\t').append(strand).append('\t')
						.append(amino).append('\t').append(anticodon).append('\t').append(length).append('\t')
						.append(intron?"1":"0").append('\n');
					totalFN++;
					byDomain.get(domain)[1]++;
					byAmino.computeIfAbsent(amino, k->new int[1])[0]++;
					lenBins[lenBin(length)]++;
					if(intron){totalIntronFN++; byDomain.get(domain)[2]++;}
				}
			}
		}

		fileIO.ReadWrite.writeString(tsv, outPrefix+"_fn.tsv");

		System.err.println("=== tRNA FALSE-NEGATIVE characterization ===");
		System.err.println("Genomes processed: "+processedGenomes+" (missing/skipped: "+missingGenomes+")");
		System.err.printf("Total reference tRNAs: %d | FN (unmatched): %d (%.1f%% miss = recall %.1f%%)%n",
			totalRef, totalFN, 100.0*totalFN/Math.max(1,totalRef), 100.0*(totalRef-totalFN)/Math.max(1,totalRef));
		System.err.println("FN written: "+outPrefix+"_fn.tsv");

		int totalCall=0, totalFP=0;
		System.err.println("\n-- by domain (ref / FN / recall | calls / FP / precision | FN-intron>100nt) --");
		for(Map.Entry<String,int[]> e : byDomain.entrySet()){
			final int[] v=e.getValue();
			final int tp=v[3]-v[4];
			totalCall+=v[3]; totalFP+=v[4];
			System.err.printf("  %-10s ref=%-5d FN=%-4d recall=%.1f%% | calls=%-5d FP=%-4d prec=%.1f%% | FN_intron=%d%n",
				e.getKey(), v[0], v[1], 100.0*(v[0]-v[1])/Math.max(1,v[0]),
				v[3], v[4], 100.0*tp/Math.max(1,v[3]), v[2]);
		}
		System.err.printf("  %-10s ref=%-5d FN=%-4d recall=%.1f%% | calls=%-5d FP=%-4d prec=%.1f%%%n",
			"TOTAL", totalRef, totalFN, 100.0*(totalRef-totalFN)/Math.max(1,totalRef),
			totalCall, totalFP, 100.0*(totalCall-totalFP)/Math.max(1,totalCall));

		System.err.println("\n-- FN by amino acid (which families are missed) --");
		// sort by count desc
		final ArrayList<Map.Entry<String,int[]>> aminos=new ArrayList<>(byAmino.entrySet());
		aminos.sort((x,y)->y.getValue()[0]-x.getValue()[0]);
		for(Map.Entry<String,int[]> e : aminos){
			System.err.printf("  %-6s %4d (%.1f%% of FN)%n", e.getKey(), e.getValue()[0], 100.0*e.getValue()[0]/Math.max(1,totalFN));
		}

		System.err.println("\n-- FN by length bin --");
		final String[] binLabel={"<70","70-79","80-89","90-99","100-119",">=120"};
		for(int i=0; i<6; i++){
			System.err.printf("  %-8s %4d (%.1f%%)%n", binLabel[i], lenBins[i], 100.0*lenBins[i]/Math.max(1,totalFN));
		}
		System.err.printf("%n>100nt (likely intron-bearing): %d / %d FN (%.1f%%)%n",
			totalIntronFN, totalFN, 100.0*totalIntronFN/Math.max(1,totalFN));
	}

	/** Ref is matched iff any call on the same seqid overlaps by reciprocal >50%
	 * (2*overlap > max(refLen, qryLen)) — identical predicate to CompareGff.overlapMetrics. */
	private static boolean overlapsAny(GffLine r, ArrayList<GffLine> calls){
		if(calls==null){return false;}
		final int rs=r.start, re=r.stop, refLen=re-rs;//no-+1 convention (matches CompareGff)
		for(GffLine q : calls){
			if(q.start>re){break;}//sorted by start
			final int ov=Math.min(re, q.stop)-Math.max(rs, q.start);
			final int qryLen=q.stop-q.start;
			if(2*ov>Math.max(refLen, qryLen)){return true;}
		}
		return false;
	}

	private static String[] refAmino(String attr){
		if(attr!=null){
			Matcher m=REF_NOTE.matcher(attr);
			if(m.find()){return new String[]{m.group(1), m.group(2)};}
			m=REF_PRODUCT.matcher(attr);
			if(m.find()){return new String[]{m.group(1), null};}
		}
		return new String[]{"UNK", null};
	}

	private static int lenBin(int len){
		if(len<70){return 0;}
		if(len<80){return 1;}
		if(len<90){return 2;}
		if(len<100){return 3;}
		if(len<120){return 4;}
		return 5;
	}

	private static ArrayList<GffLine> loadTrna(String path){
		final FileFormat ff=FileFormat.testInput(path, FileFormat.GFF, null, true, true);
		return GffLine.loadGffFile(ff, "tRNA", false);
	}

	/** Combined-GFF cache: path -> (raw seqid -> its tRNA calls). Loaded once per file. */
	private static final HashMap<String,HashMap<String,ArrayList<GffLine>>> COMBINED=new HashMap<>();
	private static HashMap<String,ArrayList<GffLine>> loadCombinedBySeqid(String path){
		HashMap<String,ArrayList<GffLine>> m=COMBINED.get(path);
		if(m!=null){return m;}
		m=new HashMap<>();
		for(GffLine q : loadTrna(path)){
			m.computeIfAbsent(q.seqid, k->new ArrayList<GffLine>()).add(q);
		}
		COMBINED.put(path, m);
		return m;
	}

	private static String stripFna(String name){
		if(name.endsWith(".fna.gz")){return name.substring(0, name.length()-7);}
		if(name.endsWith(".fna")){return name.substring(0, name.length()-4);}
		return name;
	}

	private static ArrayList<String> readList(String path){
		final ArrayList<String> out=new ArrayList<>();
		try(java.io.BufferedReader br=new java.io.BufferedReader(new java.io.FileReader(path))){
			String line;
			while((line=br.readLine())!=null){
				line=line.trim();
				if(line.length()>0 && !line.startsWith("#")){out.add(line);}
			}
		}catch(Exception e){e.printStackTrace();}
		return out;
	}

	private static final java.util.Comparator<GffLine> START_COMP=
		(x,y)->Integer.compare(x.start, y.start);
}
