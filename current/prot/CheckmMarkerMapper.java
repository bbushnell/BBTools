package prot;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import parse.PreParser;
import shared.KillSwitch;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Maps CheckM v1 domain-level single-copy marker genes onto our family set, producing
 * the Tier-B subnet subset files checkm1_bac.txt / checkm1_arc.txt (one family rank per
 * line, the exact contract MagQCVectorMaker.loadRanks consumes).
 *
 * <p>CheckM's markers are Pfam/TIGRFAM HMM accessions; our families are mmseqs clusters
 * with no Pfam label, so the correspondence is by homology: hmmsearch (run upstream, with
 * --cut_ga = CheckM's own per-model gathering thresholds) of our family representative
 * sequences against CheckM's HMM library, whose per-hit query accession is a marker. A
 * family is an input to the domain subnet iff its rep hits at least one marker in that
 * domain's set. Nothing CheckM-derived enters the trained product — only which of OUR
 * family ranks become inputs.
 *
 * <p>Inputs:
 * <ul>
 * <li>familylist= our familylist (#rank\trep_id\tocc_total; header lines start '#').
 * <li>taxonsets= CheckM's taxon_marker_sets.tsv (7 tab cols; no header). The two rows
 *     with col1=="domain" and col2 in {Bacteria,Archaea} carry col5=numMarkerGenes and
 *     col7=the marker sets as Python "[set([...]),...]" of PF#####(.v)/TIGR##### accessions.
 * <li>hits= hmmsearch --tblout output (whitespace-delimited; '#' comment lines). col1=target
 *     sequence name (our rep_id); col4=query accession (the marker, e.g. PF01000.21).
 * </ul>
 * Outputs: outbac=/outarc= the two rank-list subset files; a coverage report to stderr
 * (and report= if given): markers-covered / declared-count per domain, families per subset,
 * and the uncovered markers by accession.
 *
 * <p>Accession matching is version-stripped (PF01000.21 -> PF01000) so a Pfam version tick
 * between the marker-set table and the HMM library can never silently zero the intersection.
 *
 * @author UMP45
 */
public class CheckmMarkerMapper {

	public static void main(String[] args){
		final long t0=System.nanoTime();
		{
			PreParser pp=new PreParser(args, CheckmMarkerMapper.class, false);
			args=pp.args;
		}
		String familylist=null, taxonsets=null, hits=null, outbac=null, outarc=null, report=null;
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final int eq=arg.indexOf('=');
			final String a=(eq<0 ? arg : arg.substring(0, eq)).toLowerCase();
			final String b=(eq<0 ? null : arg.substring(eq+1));
			if(a.equals("familylist")){familylist=b;}
			else if(a.equals("taxonsets") || a.equals("taxon_marker_sets") || a.equals("markersets")){taxonsets=b;}
			else if(a.equals("hits") || a.equals("tblout")){hits=b;}
			else if(a.equals("outbac")){outbac=b;}
			else if(a.equals("outarc")){outarc=b;}
			else if(a.equals("report")){report=b;}
			else{throw new RuntimeException("Unknown parameter: "+arg);}
		}
		if(familylist==null || taxonsets==null || hits==null || outbac==null || outarc==null){
			throw new RuntimeException("Required: familylist= taxonsets= hits= outbac= outarc= [report=]");
		}

		CheckmMarkerMapper mapper=new CheckmMarkerMapper();
		mapper.repToRank=mapper.loadFamilyList(familylist);
		final int numFam=mapper.numFam;

		DomainSet bac=mapper.loadDomain(taxonsets, "Bacteria");
		DomainSet arc=mapper.loadDomain(taxonsets, "Archaea");

		final boolean[] bacMask=new boolean[numFam];
		final boolean[] arcMask=new boolean[numFam];
		mapper.mapHits(hits, bac, arc, bacMask, arcMask);

		final int[] bacRanks=writeSubset(outbac, bacMask);
		final int[] arcRanks=writeSubset(outarc, arcMask);

		final StringBuilder rep=new StringBuilder();
		rep.append("CheckmMarkerMapper coverage:\n");
		rep.append(domainReport("Bacteria", bac, bacRanks.length));
		rep.append(domainReport("Archaea",  arc, arcRanks.length));
		final String reportStr=rep.toString();
		System.err.print(reportStr);
		if(report!=null){
			ByteStreamWriter bsw=new ByteStreamWriter(FileFormat.testOutput(report, FileFormat.TXT, null, true, true, false, false));
			bsw.start(); bsw.print(new ByteBuilder().append(reportStr)); bsw.poisonAndWait();
		}
		System.err.println("Wrote "+bacRanks.length+" bac ranks -> "+outbac+", "+arcRanks.length
				+" arc ranks -> "+outarc+" in "+String.format("%.2fs", (System.nanoTime()-t0)/1e9));
	}

	/*--------------------------------------------------------------*/
	/*----------------          Loaders             ----------------*/
	/*--------------------------------------------------------------*/

	/** rep_id -> rank; also sets numFam (== number of data rows). */
	private HashMap<String,Integer> loadFamilyList(String file){
		final HashMap<String,Integer> map=new HashMap<String,Integer>();
		int maxRank=-1;
		final ByteFile bf=ByteFile.makeByteFile(file, true);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			final int tab=Tools.indexOf(line, (byte)'\t');
			assert(tab>0) : "familylist row lacks a tab: "+new String(line);
			final int rank=Parse.parseInt(line, 0, tab);
			int tab2=tab+1;
			while(tab2<line.length && line[tab2]!='\t'){tab2++;}
			final String rep=new String(line, tab+1, tab2-(tab+1));
			final Integer prev=map.put(rep, Integer.valueOf(rank));
			assert(prev==null) : "duplicate rep_id in familylist: "+rep+" (ranks "+prev+" and "+rank+")";
			maxRank=Tools.max(maxRank, rank);
		}
		bf.close();
		numFam=maxRank+1;
		assert(numFam==map.size()) : "familylist ranks not contiguous 0.."+(numFam-1)
				+": maxRank="+maxRank+" but "+map.size()+" distinct rep_ids (familylist "+file+")";
		System.err.println("familylist: "+map.size()+" families (ranks 0.."+maxRank+")");
		return map;
	}

	/** Extracts one domain's marker accessions (version-stripped) from taxon_marker_sets.tsv. */
	private DomainSet loadDomain(String file, String taxon){
		final DomainSet ds=new DomainSet(taxon);
		final ByteFile bf=ByteFile.makeByteFile(file, true);
		int rows=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			final String[] f=new String(line).split("\t");
			if(f.length<7 || !f[0].equals("domain") || !f[1].equals(taxon)){continue;}
			rows++;
			ds.declared=Integer.parseInt(f[3+1]); // col5 = numMarkerGenes (0-based index 4)
			final Matcher m=ACC.matcher(f[6]);   // col7 = the "[set([...]),...]" marker sets
			while(m.find()){ds.acc.add(stripVersion(m.group()));}
		}
		bf.close();
		assert(rows==1) : "expected exactly 1 domain/"+taxon+" row in "+file+", found "+rows;
		// The collocated sets partition the markers (disjoint), so |accessions| == numMarkerGenes.
		// A mismatch means the col7 Python-literal parse is wrong -> crash before it silently
		// corrupts the subset. Source of the invariant: taxon_marker_sets.tsv col5.
		if(ds.acc.size()!=ds.declared){
			KillSwitch.kill("domain/"+taxon+": extracted "+ds.acc.size()+" marker accessions but col5 "
					+"(numMarkerGenes) declares "+ds.declared+" in "+file+"; col7 set() parse is wrong.");
		}
		System.err.println("CheckM domain/"+taxon+": "+ds.acc.size()+" marker accessions (matches declared count)");
		return ds;
	}

	/** Streams hmmsearch tblout; records, per domain, which family ranks and which markers were hit. */
	private void mapHits(String file, DomainSet bac, DomainSet arc, boolean[] bacMask, boolean[] arcMask){
		final ByteFile bf=ByteFile.makeByteFile(file, true);
		long lines=0, accShaped=0, resolved=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			lines++;
			final String repId=token(line, 0);
			final String rawAcc=token(line, 3); // col4 = query (profile) accession
			final String acc=stripVersion(rawAcc);
			if(SHAPE.matcher(rawAcc).matches()){accShaped++;}
			final Integer rankI=repToRank.get(repId);
			if(rankI==null){continue;}
			resolved++;
			final int rank=rankI.intValue();
			assert(rank>=0 && rank<bacMask.length) : "rank "+rank+" out of [0,"+bacMask.length+") for rep "+repId;
			if(bac.acc.contains(acc)){bacMask[rank]=true; bac.hit.add(acc);}
			if(arc.acc.contains(acc)){arcMask[rank]=true; arc.hit.add(acc);}
		}
		bf.close();
		// Guards against the silent-garbage class: a wrong column mapping or wrong file gives
		// col4-not-an-accession and/or col1-not-a-rep_id. Crash loud rather than emit an empty
		// or nonsensical subset.
		if(lines>0 && accShaped<0.9*lines){
			KillSwitch.kill("hmmsearch tblout col4 doesn't look like marker accessions ("+accShaped+"/"+lines
					+" match PF/TIGR) in "+file+"; column mapping or input file is wrong.");
		}
		if(lines>0 && resolved<0.9*lines){
			KillSwitch.kill("hmmsearch tblout col1 rarely resolves to a family rank ("+resolved+"/"+lines
					+") in "+file+"; the seqdb was not consensus_reps_v4 or column mapping is wrong.");
		}
		System.err.println("hits: "+lines+" hit lines ("+accShaped+" accession-shaped, "+resolved+" resolved to a family)");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Output             ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes the set-bits of mask as ascending ranks, one per line; returns the ranks. */
	private static int[] writeSubset(String file, boolean[] mask){
		final structures.IntList l=new structures.IntList();
		for(int r=0; r<mask.length; r++){if(mask[r]){l.add(r);}}
		final ByteStreamWriter bsw=new ByteStreamWriter(FileFormat.testOutput(file, FileFormat.TXT, null, true, true, false, false));
		bsw.start();
		final ByteBuilder bb=new ByteBuilder();
		for(int i=0; i<l.size; i++){bb.append(l.get(i)).nl();}
		bsw.print(bb);
		bsw.poisonAndWait();
		return l.toArray();
	}

	private static String domainReport(String taxon, DomainSet ds, int families){
		final ArrayList<String> missing=new ArrayList<String>();
		for(String a : ds.acc){if(!ds.hit.contains(a)){missing.add(a);}}
		Collections.sort(missing);
		final StringBuilder sb=new StringBuilder();
		sb.append("  ").append(taxon).append(": ").append(ds.hit.size()).append('/').append(ds.acc.size())
			.append(" markers covered by >=1 family; ").append(families).append(" families in subset; ")
			.append(missing.size()).append(" markers uncovered.\n");
		if(!missing.isEmpty()){sb.append("    uncovered: ").append(missing).append('\n');}
		return sb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	/** The idx-th whitespace-delimited token of a byte-line (0-based), as a String. */
	private static String token(byte[] line, int idx){
		int a=0, n=line.length, field=0;
		while(a<n && line[a]<=' '){a++;}
		while(field<idx && a<n){
			while(a<n && line[a]>' '){a++;}   // skip this token
			while(a<n && line[a]<=' '){a++;}   // skip whitespace
			field++;
		}
		if(a>=n){return "";}
		int b=a;
		while(b<n && line[b]>' '){b++;}
		return new String(line, a, b-a);
	}

	/** PF01000.21 -> PF01000; TIGR01080 -> TIGR01080. */
	private static String stripVersion(String acc){
		final int dot=acc.indexOf('.');
		return dot<0 ? acc : acc.substring(0, dot);
	}

	/*--------------------------------------------------------------*/

	private int numFam;
	private HashMap<String,Integer> repToRank;

	/** Any PF/TIGR accession token, optionally Pfam-versioned; used to scan the set() literal. */
	private static final Pattern ACC=Pattern.compile("PF\\d{5}(?:\\.\\d+)?|TIGR\\d{5}");
	/** Whole-token shape test for a tblout accession field. */
	private static final Pattern SHAPE=Pattern.compile("(?:PF\\d{5}(?:\\.\\d+)?|TIGR\\d{5})");

	private static final class DomainSet {
		DomainSet(String taxon){this.taxon=taxon;}
		final String taxon;
		int declared=-1;
		final HashSet<String> acc=new HashSet<String>();
		final HashSet<String> hit=new HashSet<String>();
	}
}
