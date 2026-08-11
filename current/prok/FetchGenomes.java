package prok;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;

import tax.TaxNode;
import tax.TaxTree;

/**
 * Modern replacement for FetchProks: selects a phylogenetically diverse set of RefSeq
 * genomes from NCBI's assembly_summary.txt (no FTP directory crawling) and writes a
 * hardened download script.
 *
 * <p>Selection: rows are filtered (latest version, full genome_rep, valid ftp_path,
 * optionally not excluded_from_refseq), ranked within each species (refseq_category
 * &gt; assembly_level &gt; fewer contigs &gt; newer), then sampled with taxonomy-aware
 * quotas (maxperspecies/maxpergenus/maxperfamily) using the TaxTree, so diversity is
 * measured against real lineage rather than name tokens. Columns are resolved from the
 * header by NAME, so format drift does not break parsing.
 *
 * <p>The emitted script fixes FetchProks' silent-blank-file failure mode: every file is
 * fetched to a temp name with curl --fail + retries, integrity-checked with gzip -t,
 * renamed only on success; fna+gff are kept as a pair or not at all; failures append to
 * failed.txt. Sequence headers gain a _tid_&lt;taxid&gt; suffix (no server dependency —
 * the taxid comes from the summary itself).
 *
 * <p>Usage: java prok.FetchGenomes summary=bacteria.txt,archaea.txt tree=auto
 *   out=fetch.sh [maxperspecies=1] [maxpergenus=2] [maxperfamily=0] [skip=tids.txt]
 *
 * <p>Get the summaries with:
 *   curl -O https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
 *
 * @author Eru
 * @date 2026-08-10
 */
public class FetchGenomes {

	public static void main(String[] args){
		String summaryFiles=null, treeFile=null, out=null, skipFile=null;
		int maxPerSpecies=1, maxPerGenus=2, maxPerFamily=0;
		boolean allowExcluded=false, rename=true;
		long minSize=0;
		for(String arg : args){
			int eq=arg.indexOf('=');
			if(eq<0){continue;}
			String a=arg.substring(0, eq).toLowerCase(), b=arg.substring(eq+1);
			if(a.equals("summary") || a.equals("in")){summaryFiles=b;}
			else if(a.equals("tree")){treeFile=b;}
			else if(a.equals("out")){out=b;}
			else if(a.equals("skip")){skipFile=b;}
			else if(a.equals("maxperspecies")){maxPerSpecies=Integer.parseInt(b);}
			else if(a.equals("maxpergenus")){maxPerGenus=Integer.parseInt(b);}
			else if(a.equals("maxperfamily")){maxPerFamily=Integer.parseInt(b);}
			else if(a.equals("allowexcluded")){allowExcluded=parseBool(b);}
			else if(a.equals("rename")){rename=parseBool(b);}
			else if(a.equals("minsize")){minSize=Long.parseLong(b);}
			else{System.err.println("Warning: unknown arg "+arg);}
		}
		if(summaryFiles==null || out==null){
			throw new RuntimeException("Required: summary=<assembly_summary.txt[,file2]> out=<script.sh> "
				+"[tree=auto] [maxperspecies=1] [maxpergenus=2] [maxperfamily=0] [skip=tids.txt]");
		}

		TaxTree tree=null;
		if(treeFile!=null){
			if(treeFile.equalsIgnoreCase("auto")){treeFile=TaxTree.defaultTreeFile();}
			tree=TaxTree.loadTaxTree(treeFile, System.err, false, false);
		}

		final HashSet<Integer> skip=new HashSet<Integer>();
		if(skipFile!=null){
			try{
				BufferedReader br=new BufferedReader(new FileReader(skipFile));
				String line;
				while((line=br.readLine())!=null){
					line=line.trim();
					if(line.length()>0){skip.add(Integer.parseInt(line));}
				}
				br.close();
				System.err.println("skip list: "+skip.size()+" taxids");
			}catch(Exception e){throw new RuntimeException(e);}
		}

		//1) parse + filter all summaries
		ArrayList<Row> rows=new ArrayList<Row>();
		long parsed=0, filtered=0;
		for(String f : summaryFiles.split(",")){
			try{
				BufferedReader br=new BufferedReader(new FileReader(f), 1<<20);
				HashMap<String,Integer> col=null;
				String line;
				while((line=br.readLine())!=null){
					if(line.length()==0){continue;}
					if(line.charAt(0)=='#'){
						if(line.startsWith("#assembly_accession") || line.startsWith("# assembly_accession")){
							col=headerMap(line);
						}
						continue;
					}
					if(col==null){throw new RuntimeException("No header line before data in "+f);}
					parsed++;
					Row r=parseRow(line.split("\t", -1), col);
					if(r==null){filtered++; continue;}
					if(minSize>0 && r.genomeSize>0 && r.genomeSize<minSize){filtered++; continue;}
					if(!allowExcluded && r.excluded){filtered++; continue;}
					if(skip.contains(r.taxid)){filtered++; continue;}
					rows.add(r);
				}
				br.close();
			}catch(RuntimeException e){throw e;
			}catch(Exception e){throw new RuntimeException(e);}
		}
		System.err.println("parsed="+parsed+" rows, filtered="+filtered+", candidates="+rows.size());

		//2) best assembly per species (species_taxid), up to maxPerSpecies
		HashMap<Integer,ArrayList<Row>> bySpecies=new HashMap<Integer,ArrayList<Row>>();
		for(Row r : rows){
			ArrayList<Row> l=bySpecies.get(r.speciesTaxid);
			if(l==null){bySpecies.put(r.speciesTaxid, l=new ArrayList<Row>());}
			l.add(r);
		}
		final Comparator<Row> byQuality=new Comparator<Row>(){
			@Override
			public int compare(Row a, Row b){
				if(a.categoryScore!=b.categoryScore){return b.categoryScore-a.categoryScore;}
				if(a.levelScore!=b.levelScore){return b.levelScore-a.levelScore;}
				if(a.contigCount!=b.contigCount){return a.contigCount-b.contigCount;}
				return b.date.compareTo(a.date);
			}
		};
		ArrayList<Row> speciesBest=new ArrayList<Row>();
		for(ArrayList<Row> l : bySpecies.values()){
			Collections.sort(l, byQuality);
			for(int i=0; i<l.size() && i<Math.max(1, maxPerSpecies); i++){speciesBest.add(l.get(i));}
		}
		System.err.println("species="+bySpecies.size()+", after per-species cap="+speciesBest.size());

		//3) taxonomy-aware quotas: walk in quality order, cap per genus/family via the TaxTree.
		//A species whose lineage is unresolvable falls back to the name's genus token.
		Collections.sort(speciesBest, byQuality);
		HashMap<String,Integer> genusCount=new HashMap<String,Integer>();
		HashMap<Integer,Integer> familyCount=new HashMap<Integer,Integer>();
		TreeMap<String,Integer> phylumKept=new TreeMap<String,Integer>();
		ArrayList<Row> kept=new ArrayList<Row>();
		int newTaxids=0;
		for(Row r : speciesBest){
			String genusKey=null;
			int familyId=-1;
			String phylum="unresolved";
			if(tree!=null){
				//A fresh summary can hold taxids newer than the tree dump; probe with
				//skipAssertion and fall back taxid -> species_taxid -> name token.
				final int lookup=(tree.getNode(r.taxid, true)!=null ? r.taxid
					: tree.getNode(r.speciesTaxid, true)!=null ? r.speciesTaxid : -1);
				if(lookup<0){newTaxids++;}
				else{
					TaxNode gn=tree.getNodeAtLevel(lookup, TaxTree.GENUS);
					TaxNode fn=tree.getNodeAtLevel(lookup, TaxTree.FAMILY);
					TaxNode pn=tree.getNodeAtLevel(lookup, TaxTree.PHYLUM);
					if(gn!=null){genusKey="g"+gn.id;}
					if(fn!=null){familyId=fn.id;}
					if(pn!=null && pn.name!=null){phylum=pn.name;}
				}
			}
			if(genusKey==null){genusKey="n"+nameGenus(r.orgName);}
			int gc=genusCount.getOrDefault(genusKey, 0);
			if(maxPerGenus>0 && gc>=maxPerGenus){continue;}
			if(maxPerFamily>0 && familyId>=0){
				int fc=familyCount.getOrDefault(familyId, 0);
				if(fc>=maxPerFamily){continue;}
				familyCount.put(familyId, fc+1);
			}
			genusCount.put(genusKey, gc+1);
			phylumKept.put(phylum, phylumKept.getOrDefault(phylum, 0)+1);
			kept.add(r);
		}
		System.err.println("kept="+kept.size()+" genomes across "+genusCount.size()+" genera");
		if(newTaxids>0){
			System.err.println(newTaxids+" taxids newer than the tree dump (name-token genus fallback used)");
		}
		System.err.println("\nKept per phylum:");
		for(String p : phylumKept.keySet()){
			System.err.println(String.format("%-32s\t%d", p, phylumKept.get(p)));
		}

		//4) emit the hardened download script
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out), 1<<20);
			bw.write("#!/bin/bash\n");
			bw.write("#Generated by prok.FetchGenomes: "+kept.size()+" fna+gff pairs.\n");
			bw.write("#Safe to re-run: existing non-empty pairs are skipped. Failures append to failed.txt.\n");
			bw.write("set -u\n");
			bw.write("RENAME="+(rename ? 1 : 0)+"\n\n");
			bw.write("fetchone(){ #url out\n");
			bw.write("\tlocal url=\"$1\" out=\"$2\" i\n");
			bw.write("\tfor i in 1 2 3; do\n");
			bw.write("\t\tif curl -sS --fail --connect-timeout 30 --retry 2 -o \"$out.tmp\" \"$url\" "
				+"&& gzip -t \"$out.tmp\" 2>/dev/null; then\n");
			bw.write("\t\t\tmv \"$out.tmp\" \"$out\"; return 0\n");
			bw.write("\t\tfi\n");
			bw.write("\t\trm -f \"$out.tmp\"; sleep $((i*5))\n");
			bw.write("\tdone\n");
			bw.write("\trm -f \"$out.tmp\"; return 1\n");
			bw.write("}\n\n");
			bw.write("fetchpair(){ #fnaUrl gffUrl prefix tid\n");
			bw.write("\tlocal fna=\"$1\" gff=\"$2\" prefix=\"$3\" tid=\"$4\"\n");
			bw.write("\tif [ -s \"$prefix.fna.gz\" ] && [ -s \"$prefix.gff.gz\" ]; then return 0; fi\n");
			bw.write("\tif fetchone \"$fna\" \"$prefix.fna.gz\" && fetchone \"$gff\" \"$prefix.gff.gz\"; then\n");
			bw.write("\t\tif [ \"$RENAME\" = 1 ]; then\n");
			bw.write("\t\t\tzcat \"$prefix.fna.gz\" | sed -E \"s/^>([^ \\t]+).*/>\\\\1_tid_${tid}/\" "
				+"| gzip > \"$prefix.fna.gz.tmp\" && mv \"$prefix.fna.gz.tmp\" \"$prefix.fna.gz\"\n");
			bw.write("\t\tfi\n");
			bw.write("\t\techo \"ok: $prefix\"\n");
			bw.write("\telse\n");
			bw.write("\t\trm -f \"$prefix.fna.gz\" \"$prefix.gff.gz\"\n");
			bw.write("\t\techo \"$prefix\" >> failed.txt\n");
			bw.write("\t\techo \"FAILED: $prefix\" >&2\n");
			bw.write("\tfi\n");
			bw.write("}\n\n");
			for(Row r : kept){
				String p=r.ftpPath;
				while(p.endsWith("/")){p=p.substring(0, p.length()-1);}
				String base=p.substring(p.lastIndexOf('/')+1);
				String root=p+"/"+base;
				if(root.startsWith("ftp://")){root="https://"+root.substring(6);}
				String prefix="tid_"+r.taxid+"_"+sanitize(r.orgName);
				bw.write("fetchpair "+root+"_genomic.fna.gz "+root+"_genomic.gff.gz "+prefix+" "+r.taxid+"\n");
			}
			bw.write("\necho \"done; failures (if any) are listed in failed.txt\"\n");
			bw.close();
			System.err.println("\nwrote "+out+" ("+kept.size()+" pairs)");
		}catch(Exception e){throw new RuntimeException(e);}
	}

	/** Maps header column names to indices; the leading '#' is stripped from the first name. */
	private static HashMap<String,Integer> headerMap(String line){
		String[] h=line.substring(1).trim().split("\t");
		HashMap<String,Integer> m=new HashMap<String,Integer>();
		for(int i=0; i<h.length; i++){m.put(h[i].trim(), i);}
		for(String req : new String[]{"assembly_accession","refseq_category","taxid","species_taxid",
			"organism_name","version_status","assembly_level","genome_rep","seq_rel_date","ftp_path"}){
			if(!m.containsKey(req)){throw new RuntimeException("Summary header missing column: "+req);}
		}
		return m;
	}

	/** Parses one data row; returns null if it fails the basic filters. */
	private static Row parseRow(String[] f, HashMap<String,Integer> col){
		if(f.length<col.size()-4){return null;}//tolerate a few trailing-column absences
		String status=get(f, col, "version_status");
		if(!"latest".equals(status)){return null;}
		String rep=get(f, col, "genome_rep");
		if(!"Full".equals(rep)){return null;}
		String ftp=get(f, col, "ftp_path");
		if(ftp==null || ftp.length()<8 || ftp.equals("na")){return null;}
		Row r=new Row();
		try{
			r.taxid=Integer.parseInt(get(f, col, "taxid"));
			r.speciesTaxid=Integer.parseInt(get(f, col, "species_taxid"));
		}catch(NumberFormatException e){return null;}
		r.orgName=get(f, col, "organism_name");
		r.ftpPath=ftp;
		r.date=nvl(get(f, col, "seq_rel_date"));
		String cat=nvl(get(f, col, "refseq_category"));
		r.categoryScore=(cat.contains("reference") ? 2 : cat.contains("representative") ? 1 : 0);
		String level=nvl(get(f, col, "assembly_level"));
		r.levelScore=("Complete Genome".equals(level) ? 3 : "Chromosome".equals(level) ? 2
			: "Scaffold".equals(level) ? 1 : 0);
		String ex=get(f, col, "excluded_from_refseq");
		r.excluded=(ex!=null && ex.length()>0 && !ex.equals("na"));
		r.contigCount=parseIntSafe(get(f, col, "contig_count"), Integer.MAX_VALUE);
		r.genomeSize=parseLongSafe(get(f, col, "genome_size"), 0);
		return r;
	}

	private static String get(String[] f, HashMap<String,Integer> col, String name){
		Integer i=col.get(name);
		return (i==null || i>=f.length) ? null : f[i];
	}
	private static String nvl(String s){return s==null ? "" : s;}
	private static int parseIntSafe(String s, int dflt){
		if(s==null || s.length()==0 || s.equals("na")){return dflt;}
		try{return Integer.parseInt(s);}catch(NumberFormatException e){return dflt;}
	}
	private static long parseLongSafe(String s, long dflt){
		if(s==null || s.length()==0 || s.equals("na")){return dflt;}
		try{return Long.parseLong(s);}catch(NumberFormatException e){return dflt;}
	}
	private static boolean parseBool(String s){
		return s==null || s.equals("t") || s.equals("true") || s.equals("1") || s.equals("yes");
	}

	/** First name token, Candidatus-aware, for lineage fallback when the tree lacks the taxid. */
	private static String nameGenus(String name){
		if(name==null){return "unknown";}
		String s=name.trim();
		if(s.startsWith("Candidatus ")){s=s.substring("Candidatus ".length());}
		int sp=s.indexOf(' ');
		return sp>0 ? s.substring(0, sp) : s;
	}

	/** Filesystem-safe organism name: alnum preserved, runs of anything else collapse to '_'. */
	private static String sanitize(String name){
		StringBuilder sb=new StringBuilder();
		boolean us=false;
		for(int i=0; i<name.length() && sb.length()<80; i++){
			char c=name.charAt(i);
			if(Character.isLetterOrDigit(c)){sb.append(c); us=false;}
			else if(!us && sb.length()>0){sb.append('_'); us=true;}
		}
		while(sb.length()>0 && sb.charAt(sb.length()-1)=='_'){sb.setLength(sb.length()-1);}
		return sb.length()>0 ? sb.toString() : "unnamed";
	}

	/** One filtered assembly_summary row's selection-relevant fields. */
	private static final class Row {
		int taxid, speciesTaxid, categoryScore, levelScore, contigCount;
		long genomeSize;
		boolean excluded;
		String orgName, ftpPath, date;
	}
}
