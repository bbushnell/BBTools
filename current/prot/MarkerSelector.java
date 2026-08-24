package prot;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;

/**
 * Per-phylum universal single-copy marker selection (Tier B subnet precompute) -- a Java port of
 * the historical {@code marker_select.py} (verified line-for-line against that source before
 * porting; see mag-qc {@code plans/TIER_B_DESIGN.md}). Reads the SPARSE per-org cache
 * ({@link PerOrgRollup}'s {@code out2=} output) and {@code taxpgm.tsv} (tid-&gt;phylum), and for
 * every phylum group with &gt;= minorgs organisms (plus a cross-phylum {@code ALL_<domain>} row),
 * selects family ranks that are MARKERS: present in &gt;= minprev fraction of the group's
 * organisms, AND single-copy in &gt;= minsc fraction of the organisms where present.
 *
 * <p>Unlike the original (which ran once per domain via a {@code bac}/{@code arc} CLI arg), this
 * processes BOTH domains in one pass over the per-org cache -- phylum names never collide across
 * domains in NCBI taxonomy, so a single {@code phylum -> Group} map is safe and simpler than two
 * invocations. The {@code ALL_<domain>} cross-phylum row is still emitted per domain, matching the
 * original semantics exactly.
 *
 * <p>Family ranks in the sparse per-org cache are ALREADY familylist_v4 ranks (0-based, assigned
 * by {@code CacheBuilder.loadFamilyList} from the same familylist file) -- no rep-name remap is
 * needed here, unlike the old pipeline's rep-id intermediate step.
 *
 * <p>Two separate runs make the two subset families used by Tier B:
 * {@code minprev=0.95 minsc=0.90} (marker) and {@code minprev=0.90 minsc=0.90} (m90) -- same tool,
 * different floor, per TIER_B_DESIGN.md.
 *
 * <p>Usage: markerselector.sh perorg=perorg_sparse_v4.tsv taxpgm=taxpgm.tsv out=markersets_v4.tsv
 *          outdir=subsets/ tag=marker [minprev=0.95] [minsc=0.90] [minorgs=3]
 *
 * @author Eru
 */
public class MarkerSelector {

	public static void main(String[] args){
		Timer t=new Timer();
		MarkerSelector x=new MarkerSelector(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public MarkerSelector(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("perorg") || a.equals("in")){perorgFile=b;}
			else if(a.equals("taxpgm")){taxpgmFile=b;}
			else if(a.equals("out")){outFile=b;}
			else if(a.equals("outdir")){outDir=b;}
			else if(a.equals("tag")){tag=b;}
			else if(a.equals("minprev")){minPrev=Double.parseDouble(b);}
			else if(a.equals("minsc")){minSc=Double.parseDouble(b);}
			else if(a.equals("minorgs")){minOrgs=Integer.parseInt(b);}
			else if(a.equals("ow") || a.equals("overwrite")){overwrite=Parse.parseBoolean(b);}
			else{outstream.println("Unknown parameter "+arg); assert(false) : "Unknown parameter "+arg;}
		}
		assert(perorgFile!=null) : "perorg= is required.";
		assert(taxpgmFile!=null) : "taxpgm= is required.";
		assert(outFile!=null) : "out= is required.";
		assert(outDir!=null) : "outdir= is required (per-phylum rank files).";
		assert(tag!=null) : "tag= is required (e.g. marker or m90 -- names the per-phylum rank files).";
	}

	/** Per-phylum-group accumulator: n orgs, present[rank]->#orgs-with-it,
	 *  single[rank]->#orgs-with-it-as-EXACTLY-1-copy. HashMap-keyed (not a fixed F-sized array) so
	 *  the tool needs no a-priori family-count parameter -- ranks simply appear as they're seen. */
	static final class Group {
		int nOrgs;
		final HashMap<Integer, Integer> present=new HashMap<Integer, Integer>();
		final HashMap<Integer, Integer> single=new HashMap<Integer, Integer>();

		void add(List<int[]> presRankCount){
			nOrgs++;
			for(int[] rc : presRankCount){
				final int rank=rc[0], count=rc[1];
				present.merge(rank, 1, Integer::sum);
				if(count==1){single.merge(rank, 1, Integer::sum);}
			}
		}
	}

	/** tid<TAB>phylum, matching FamilyListBuilder.loadTaxpgm / MagQCVectorMaker.loadAux exactly. */
	HashMap<Integer, String> loadTaxpgm(){
		HashMap<Integer, String> map=new HashMap<Integer, String>(1<<16);
		final ByteFile bf=ByteFile.makeByteFile(taxpgmFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			lp.set(line);
			if(lp.terms()<2){continue;}
			map.put(lp.parseInt(0), lp.parseString(1));
		}
		bf.close();
		return map;
	}

	void process(Timer t){
		final HashMap<Integer, String> tid2phylum=loadTaxpgm();
		outstream.println("taxpgm tids: "+tid2phylum.size());

		final HashMap<String, Group> groups=new HashMap<String, Group>();
		final ArrayList<int[]> presRankCount=new ArrayList<int[]>();

		final ByteFile bf=ByteFile.makeByteFile(perorgFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		long rowsRead=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			lp.set(line);
			assert(lp.terms()>=16) : "Malformed sparse per-org row: "+lp.terms()
				+" fields (need >=16): "+new String(line);

			final int tid=lp.parseInt(0);
			final String domain=lp.parseString(1);

			presRankCount.clear();
			final int flen=lp.length(15);
			if(flen>0){parseFamCounts(lp.line(), lp.a(), lp.b(), presRankCount);}

			final String allGroup="ALL_"+domain;
			groupFor(groups, allGroup).add(presRankCount);

			final String phy=tid2phylum.get(tid);
			if(phy!=null){groupFor(groups, phy).add(presRankCount);}

			rowsRead++;
		}
		bf.close();
		outstream.println("Rows read: "+rowsRead+", groups: "+groups.size());

		writeMarkerSets(groups);

		t.stop();
		outstream.println("Time: \t"+t);
	}

	private static Group groupFor(HashMap<String, Group> groups, String name){
		Group g=groups.get(name);
		if(g==null){g=new Group(); groups.put(name, g);}
		return g;
	}

	/** Selects markers per qualifying group, writes the summary TSV, and emits one per-phylum rank
	 *  file per group (skips ALL_<domain> rows for the per-phylum files -- those are summary-only,
	 *  matching the original pipeline where ALL_<domain> fed its own explicit subset name, not a
	 *  per-phylum split). */
	void writeMarkerSets(HashMap<String, Group> groups){
		final ArrayList<String> ordered=new ArrayList<String>(groups.keySet());
		ordered.sort(Comparator.<String>comparingInt(g -> -groups.get(g).nOrgs).thenComparing(g -> g));

		final FileFormat ff=FileFormat.testOutput(outFile, FileFormat.TXT, null, true, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		int nGroups=0;
		try{
			final ByteBuilder bb=new ByteBuilder(1<<16);
			bb.append("#phylum\tn_orgs\tn_markers\tranks").nl();
			bsw.print(bb);
			bb.clear();

			for(String name : ordered){
				final Group g=groups.get(name);
				final int n=g.nOrgs;
				if(n<minOrgs){continue;}

				final ArrayList<Integer> marks=new ArrayList<Integer>();
				for(Map.Entry<Integer, Integer> e : g.present.entrySet()){
					final int rank=e.getKey(), p=e.getValue();
					if(p<=0){continue;}
					if(p<minPrev*n){continue;}
					final int s=g.single.getOrDefault(rank, 0);
					if(s<minSc*p){continue;}
					marks.add(rank);
				}
				marks.sort(null);

				bb.append(name).append('\t').append(n).append('\t').append(marks.size()).append('\t');
				for(int i=0; i<marks.size(); i++){
					if(i>0){bb.append(',');}
					bb.append(marks.get(i).intValue());
				}
				bb.nl();
				bsw.print(bb);
				bb.clear();
				nGroups++;

				if(!name.startsWith("ALL_")){writeRankFile(name, marks);}
				outstream.println(name+": n="+n+" markers="+marks.size());
			}
		}finally{
			bsw.poisonAndWait();
		}
		outstream.println("Wrote "+nGroups+" qualifying groups (n_orgs>="+minOrgs+") to "+outFile);
	}

	/** One rank per line, e.g. subsets/marker_Pseudomonadota.txt or subsets/m90_Pseudomonadota.txt. */
	private void writeRankFile(String phylum, List<Integer> marks){
		final String path=outDir+(outDir.endsWith("/") ? "" : "/")+tag+"_"+phylum+".txt";
		final FileFormat ff=FileFormat.testOutput(path, FileFormat.TXT, null, true, overwrite, false, false);
		final ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		final ByteBuilder bb=new ByteBuilder(1<<12);
		for(int r : marks){bb.append(r).nl();}
		bsw.print(bb);
		bsw.poisonAndWait();
	}

	/** Parses "rank:count;rank:count;..." within [a,b) of line into (rank,count) pairs, appended
	 *  to the (caller-cleared) output list. Mirrors CacheBuilder/PerOrgRollup's identical parser. */
	private static void parseFamCounts(byte[] line, int a, int b, List<int[]> out){
		int start=a;
		for(int i=a; i<=b; i++){
			if(i==b || line[i]==';'){
				if(i>start){
					int colon=-1;
					for(int j=start; j<i; j++){if(line[j]==':'){colon=j; break;}}
					assert(colon>start) : "Malformed famcounts field (missing ':' in a rank:count pair).";
					final int rank=Parse.parseInt(line, start, colon);
					final int count=Parse.parseInt(line, colon+1, i);
					out.add(new int[]{rank, count});
				}
				start=i+1;
			}
		}
	}

	private String perorgFile=null;
	private String taxpgmFile=null;
	private String outFile=null;
	private String outDir=null;
	private String tag=null;
	private double minPrev=0.95;
	private double minSc=0.90;
	private int minOrgs=3;
	private boolean overwrite=true;
	private java.io.PrintStream outstream=System.err;
}
