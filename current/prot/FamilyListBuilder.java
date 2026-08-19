package prot;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import tax.TaxTree;

/**
 * Builds a per-phylum-selected family list from mmseqs cluster output (STEP 1
 * of the v3 net-derivation rebuild, replacing the old flat top-N-by-global-
 * prevalence selection).
 *
 * <p>Selection method (Brian, 2026-08-17; spec locked with UMP45):
 * <ul>
 * <li><b>occ</b> (prevalence) = count of DISTINCT tids with &gt;=1 member in
 *     the family -- genome count, NOT copy count. A genome contributing 3
 *     members to one family still counts once.
 * <li>Every EXCLUDED tid (corpus tid-collision casualties this session found
 *     in the foundation cache -- see records/EXCLUDED_TIDS.tsv /
 *     records/TIDCOLLISION_RECONCILIATION.md) is dropped from ALL occ
 *     counting: it never contributes a member, never counts toward a
 *     phylum's genome total.
 * <li>For each phylum with &gt;=MINPHYLUMGENOMES (default 3) non-excluded
 *     genomes: rank its candidate families (any family with &gt;=1 member of
 *     that phylum) by (occ_in_phylum desc, occ_total desc, rep_id asc), take
 *     the top PERPHYLUMTOP (default 100).
 * <li>Union every qualifying phylum's top-N set (a family can qualify via
 *     multiple phyla; dedup by rep_id).
 * <li>FLOOR, not a cap: if the union already has &gt;=FLOOR (default 6000)
 *     families, keep it whole -- no trimming. Otherwise pad with the
 *     globally-most-prevalent (occ_total desc) families not already in the
 *     union until reaching FLOOR exactly.
 * </ul>
 *
 * <p>Output familylist.tsv now contains ONLY the selected set -- a change
 * from the old all-clusters-ranked contract, where familylist.tsv held every
 * cluster and a separate topn= cutoff (applied by both this tool's repsout
 * and by CacheBuilder) picked the top slice. Rank 0..K-1, ordered by
 * occ_total descending (rep_id ascending tiebreak -- this final ordering is
 * otherwise arbitrary for training, it only needs to be deterministic). K may
 * exceed FLOOR if the per-phylum union alone already did.
 * <b>CacheBuilder's topn= must be set to the actual K this run reports, not
 * assumed to equal FLOOR.</b>
 *
 * <p>Usage: familylistbuilder.sh cluster=families_cluster.tsv reps=families_rep_seq.fasta
 *        taxpgm=taxpgm.tsv excluded=EXCLUDED_TIDS.tsv out=familylist.tsv
 *        repsout=consensus_reps.fasta [perphylumtop=100] [floor=6000]
 *        [minphylumgenomes=3] [tidsout=tids.txt]
 *
 * @author Brian Bushnell, Eru
 */
public class FamilyListBuilder {

	public static void main(String[] args){
		Timer t=new Timer();
		FamilyListBuilder x=new FamilyListBuilder(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public FamilyListBuilder(String[] args){
		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("cluster") || a.equals("clustertsv")){clusterFile=b;}
			else if(a.equals("reps") || a.equals("repfasta")){repFastaFile=b;}
			else if(a.equals("taxpgm")){taxpgmFile=b;}
			else if(a.equals("excluded") || a.equals("excludedtids")){excludedFile=b;}
			else if(a.equals("out") || a.equals("familylist")){outFile=b;}
			else if(a.equals("repsout") || a.equals("topreps")){repsOutFile=b;}
			else if(a.equals("perphylumtop") || a.equals("topx")){perPhylumTop=Integer.parseInt(b);}
			else if(a.equals("floor") || a.equals("targety")){floor=Integer.parseInt(b);}
			else if(a.equals("minphylumgenomes") || a.equals("minphylum")){minPhylumGenomes=Integer.parseInt(b);}
			else if(a.equals("tidsout") || a.equals("tidsfile")){tidsOutFile=b;}
			else if(a.equals("ow") || a.equals("overwrite")){overwrite=Parse.parseBoolean(b);}
			else{outstream.println("Unknown parameter "+arg); assert(false) : "Unknown parameter "+arg;}
		}
		assert(clusterFile!=null) : "cluster= is required.";
		assert(repFastaFile!=null) : "reps= is required.";
		assert(taxpgmFile!=null) : "taxpgm= is required.";
		assert(excludedFile!=null) : "excluded= is required.";
		assert(outFile!=null) : "out= is required.";
	}

	void process(Timer t){
		HashSet<Integer> excluded=loadExcluded();
		outstream.println("Excluded tids (never counted): "+excluded.size());

		HashMap<Integer, String> tid2phylum=loadTaxpgm();
		outstream.println("taxpgm tids: "+tid2phylum.size());

		HashMap<String, Integer> phylumGenomeCount=countPhylumGenomes(tid2phylum, excluded);
		int qualifyingPhyla=0;
		for(int c : phylumGenomeCount.values()){if(c>=minPhylumGenomes){qualifyingPhyla++;}}
		outstream.println("Distinct phyla: "+phylumGenomeCount.size()+"; qualifying (>="+minPhylumGenomes+" genomes): "+qualifyingPhyla);

		HashMap<String, HashSet<Integer>> familyTids=loadClusterMembership(excluded, tid2phylum);

		// occ_total per family, and a global rank-by-prevalence array (also used for padding).
		ArrayList<String> allReps=new ArrayList<String>(familyTids.size());
		ArrayList<Integer> allOccTotal=new ArrayList<Integer>(familyTids.size());
		for(HashMap.Entry<String, HashSet<Integer>> e : familyTids.entrySet()){
			allReps.add(e.getKey());
			allOccTotal.add(e.getValue().size());
		}
		outstream.println("Total families (non-excluded members only): "+allReps.size());

		// Per-phylum candidate lists: for every family, break its tid set down by
		// phylum once, and register it as a candidate for each phylum it touches.
		HashMap<String, ArrayList<int[]>> candidatesByPhylum=new HashMap<String, ArrayList<int[]>>(); // phylum -> [repIndex, occInPhylum]
		for(int i=0; i<allReps.size(); i++){
			HashSet<Integer> tids=familyTids.get(allReps.get(i));
			HashMap<String, Integer> perPhylum=new HashMap<String, Integer>(4);
			for(int tid : tids){
				String phy=tid2phylum.get(tid);
				assert(phy!=null) : "tid "+tid+" has no taxpgm phylum entry and is not excluded -- corpus/taxpgm mismatch";
				Integer cur=perPhylum.get(phy);
				perPhylum.put(phy, cur==null ? 1 : cur+1);
			}
			for(HashMap.Entry<String, Integer> e : perPhylum.entrySet()){
				String phy=e.getKey();
				if(phylumGenomeCount.get(phy)<minPhylumGenomes){continue;}
				ArrayList<int[]> list=candidatesByPhylum.get(phy);
				if(list==null){list=new ArrayList<int[]>(); candidatesByPhylum.put(phy, list);}
				list.add(new int[]{i, e.getValue()});
			}
		}

		// Rank each qualifying phylum's candidates by (occ_in_phylum desc, occ_total desc, rep asc); take top X.
		HashSet<String> selected=new HashSet<String>();
		for(HashMap.Entry<String, ArrayList<int[]>> e : candidatesByPhylum.entrySet()){
			ArrayList<int[]> list=e.getValue();
			list.sort((a, b) -> {
				int c=Integer.compare(b[1], a[1]); // occ_in_phylum desc
				if(c!=0){return c;}
				c=Integer.compare(allOccTotal.get(b[0]), allOccTotal.get(a[0])); // occ_total desc
				if(c!=0){return c;}
				return allReps.get(a[0]).compareTo(allReps.get(b[0])); // rep asc
			});
			int take=Tools.min(perPhylumTop, list.size());
			for(int i=0; i<take; i++){selected.add(allReps.get(list.get(i)[0]));}
		}
		outstream.println("Union of per-phylum top-"+perPhylumTop+" sets: "+selected.size()+" families ("+qualifyingPhyla+" phyla)");

		// Global rank order (by occ_total desc, rep asc) -- used both for padding
		// (walk it in order, skip already-selected) and for the final output order.
		Integer[] globalOrder=new Integer[allReps.size()];
		for(int i=0; i<globalOrder.length; i++){globalOrder[i]=i;}
		Arrays.sort(globalOrder, (a, b) -> {
			int c=Integer.compare(allOccTotal.get(b), allOccTotal.get(a));
			if(c!=0){return c;}
			return allReps.get(a).compareTo(allReps.get(b));
		});

		int paddedIn=0;
		if(selected.size()<floor){
			for(int idx : globalOrder){
				if(selected.size()>=floor){break;}
				String rep=allReps.get(idx);
				if(selected.add(rep)){paddedIn++;}
			}
		}
		outstream.println("Padded in "+paddedIn+" additional globally-prevalent families (floor="+floor+")");
		outstream.println("FINAL selected family count: "+selected.size());

		// Final output order: walk globalOrder (already occ_total desc / rep asc),
		// emit only the ones actually selected -- guarantees the union members
		// (which may not be globally top-ranked) still sort correctly relative to
		// each other and to the padding.
		ArrayList<String> finalReps=new ArrayList<String>(selected.size());
		ArrayList<Integer> finalOcc=new ArrayList<Integer>(selected.size());
		for(int idx : globalOrder){
			String rep=allReps.get(idx);
			if(selected.contains(rep)){
				finalReps.add(rep);
				finalOcc.add(allOccTotal.get(idx));
			}
		}
		assert(finalReps.size()==selected.size());

		writeFamilyList(finalReps, finalOcc);
		if(repsOutFile!=null){writeRepsOut(new HashSet<String>(finalReps));}
		if(tidsOutFile!=null){writeTidsOut(familyTids, excluded);}

		outstream.println("Prevalence distribution (selected set only):");
		int[] thresholds={1, 2, 5, 10, 50, 100, 500, 1000, 5000};
		for(int thresh : thresholds){
			int count=0;
			for(int occ : finalOcc){if(occ>=thresh){count++;}}
			outstream.println("  occ_total>="+thresh+": "+count);
		}

		t.stop();
		outstream.println("Time: \t"+t);
	}

	/** excluded_tid<TAB>class<TAB>reason, header starts with '#'. Only column 0 matters. */
	HashSet<Integer> loadExcluded(){
		HashSet<Integer> set=new HashSet<Integer>();
		final ByteFile bf=ByteFile.makeByteFile(excludedFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0 || line[0]=='#'){continue;}
			lp.set(line);
			set.add(lp.parseInt(0));
		}
		bf.close();
		return set;
	}

	/** tid<TAB>phylum, two-field contract matching MagQCVectorMaker.loadAux. */
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

	/** phylum -> count of distinct NON-EXCLUDED tids with that phylum. */
	static HashMap<String, Integer> countPhylumGenomes(HashMap<Integer, String> tid2phylum, HashSet<Integer> excluded){
		HashMap<String, Integer> counts=new HashMap<String, Integer>();
		for(HashMap.Entry<Integer, String> e : tid2phylum.entrySet()){
			if(excluded.contains(e.getKey())){continue;}
			String phy=e.getValue();
			Integer cur=counts.get(phy);
			counts.put(phy, cur==null ? 1 : cur+1);
		}
		return counts;
	}

	/** mmseqs cluster.tsv: rep_id<TAB>member_id, NOT grouped by rep. Excluded-tid
	 *  members are dropped entirely; a non-excluded member whose tid has no
	 *  taxpgm entry is a corpus-integrity failure (assert, see loop body). */
	HashMap<String, HashSet<Integer>> loadClusterMembership(HashSet<Integer> excluded, HashMap<Integer, String> tid2phylum){
		outstream.println("Reading cluster membership from "+clusterFile+"...");
		HashMap<String, HashSet<Integer>> familyTids=new HashMap<String, HashSet<Integer>>(1<<22);
		long totalMembers=0, excludedMembers=0;
		final ByteFile bf=ByteFile.makeByteFile(clusterFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			lp.set(line);
			if(lp.terms()<2){continue;}
			String rep=lp.parseString(0);
			String member=lp.parseString(1);
			int tid=parseTid(member);
			totalMembers++;
			if(excluded.contains(tid)){excludedMembers++; continue;}
			assert(tid2phylum.containsKey(tid)) : "Member tid "+tid+" ("+member
				+") is not excluded but has no taxpgm.tsv entry -- corpus/taxpgm coverage mismatch";
			HashSet<Integer> tids=familyTids.get(rep);
			if(tids==null){tids=new HashSet<Integer>(); familyTids.put(rep, tids);}
			tids.add(tid);
			if(totalMembers%5000000==0){
				outstream.println("  "+totalMembers/1000000+"M members, "+familyTids.size()+" families");
			}
		}
		bf.close();
		outstream.println("Total members: "+totalMembers+" ("+excludedMembers+" from excluded tids, skipped)");
		outstream.println("Total families with >=1 non-excluded member: "+familyTids.size());
		return familyTids;
	}

	void writeFamilyList(ArrayList<String> reps, ArrayList<Integer> occ){
		outstream.println("Writing familylist ("+reps.size()+" SELECTED families) to "+outFile+"...");
		FileFormat ff=FileFormat.testOutput(outFile, FileFormat.TXT, null, true, overwrite, false, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		try{
			ByteBuilder bb=new ByteBuilder(1<<16);
			bb.append("#rank\trep_id\tocc_total").nl();
			bsw.print(bb); bb.clear();
			for(int i=0; i<reps.size(); i++){
				bb.append(i).append('\t').append(reps.get(i)).append('\t').append(occ.get(i)).nl();
				bsw.print(bb); bb.clear();
			}
		}finally{
			bsw.poisonAndWait();
		}
	}

	void writeRepsOut(HashSet<String> selected){
		outstream.println("Extracting "+selected.size()+" selected rep sequences from "+repFastaFile+"...");
		FileFormat ff=FileFormat.testOutput(repsOutFile, FileFormat.FA, null, true, overwrite, false, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		int written=0;
		try{
			final ByteFile bf=ByteFile.makeByteFile(repFastaFile, true);
			boolean keep=false;
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length>0 && line[0]=='>'){
					int end=line.length;
					for(int i=1; i<line.length; i++){if(line[i]==' ' || line[i]=='\t'){end=i; break;}}
					String name=new String(line, 1, end-1);
					keep=selected.contains(name);
					if(keep){written++;}
				}
				if(keep){bsw.println(line);}
			}
			bf.close();
		}finally{
			bsw.poisonAndWait();
		}
		outstream.println("Wrote "+written+" rep sequences to "+repsOutFile);
	}

	void writeTidsOut(HashMap<String, HashSet<Integer>> familyTids, HashSet<Integer> excluded){
		HashSet<Integer> allTids=new HashSet<Integer>();
		for(HashSet<Integer> tids : familyTids.values()){allTids.addAll(tids);}
		Integer[] tidArr=allTids.toArray(new Integer[allTids.size()]);
		Arrays.sort(tidArr);
		FileFormat ff=FileFormat.testOutput(tidsOutFile, FileFormat.TXT, null, true, overwrite, false, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		ByteBuilder bb=new ByteBuilder();
		for(Integer tid : tidArr){bb.append(tid.intValue()).nl();}
		bsw.print(bb);
		bsw.poisonAndWait();
		outstream.println("Wrote "+tidArr.length+" distinct (non-excluded) tids to "+tidsOutFile);
	}

	/** Parse the taxid from a cluster member id. A non-positive result means a
	 *  corrupt or unfixed header -- a corpus-integrity failure that aborts the
	 *  run. Unconditional throw (NOT an -ea assert): this is a data contract on
	 *  external input, and the one context that runs without -ea is exactly
	 *  where a silent -1 would poison the foundation. Costs one compare on the
	 *  happy path. */
	static int parseTid(String name){
		int tid=TaxTree.parseTaxID(name);
		if(tid<=0){throw new RuntimeException("Non-positive tid parsed from member header — corpus "
			+"should carry only valid tids after the tid|-1| source fix: "+name);}
		return tid;
	}

	private String clusterFile=null;
	private String repFastaFile=null;
	private String taxpgmFile=null;
	private String excludedFile=null;
	private String outFile=null;
	private String repsOutFile=null;
	private String tidsOutFile=null;
	private int perPhylumTop=100;
	private int floor=6000;
	private int minPhylumGenomes=3;
	private boolean overwrite=true;
	private java.io.PrintStream outstream=System.err;
}
