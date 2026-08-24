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
 * Builds a per-phylum-selected family list from mmseqs cluster output (v4 of
 * the net-derivation family selection, replacing v3's top-100-plus-global-
 * padding scheme).
 *
 * <p>Selection method (Brian, 2026-08-18; spec locked with UMP45):
 * <ul>
 * <li><b>occ</b> (prevalence) = count of DISTINCT tids with &gt;=1 member in
 *     the family -- genome count, NOT copy count.
 * <li>Every EXCLUDED tid is dropped from ALL occ counting (see
 *     records/EXCLUDED_TIDS.tsv / records/TIDCOLLISION_RECONCILIATION.md).
 * <li>For each phylum P with &gt;=MINPHYLUMGENOMES (default 3) non-excluded
 *     genomes: candidate pool = families with &gt;=1 member in P AND
 *     occ_in_phylum(F,P) &gt;= MINOCCFRAC (default 0.50) * genomes_in_P --
 *     this floor keeps the "phylum-specific" half from degenerating into
 *     1-2-genome noise.
 * <li>TOP_N (default 100, the "universal backbone" half) = the pool's HIGHEST
 *     occ_single members -- <b>occ_single(F) = count of genomes where F
 *     appears as EXACTLY ONE copy</b> (per-(family,genome) member count==1),
 *     NOT occ_total. Tiebreak (occ_total desc, occ_in_phylum desc, rep asc).
 *     Rationale (Brian, 2026-08-18): a family with copy&gt;1 in many genomes
 *     only signals contamination risk if it is NATIVELY single-copy, so the
 *     backbone marker set must itself be built from clean single-copy
 *     families, not merely globally-prevalent ones. BOT_N (default 100, the
 *     "phylum-specific" half, UNCHANGED from the original design) = the
 *     pool's LOWEST occ_total members (tiebreak occ_in_phylum desc, rep asc)
 *     -- multi-copy is fine here, this half is not a copy-number reference.
 *     If the pool is smaller than TOP_N+BOT_N, the two halves overlap and P
 *     simply contributes its whole (smaller) pool.
 * <li>P's contribution = TOP_N union BOT_N (dedup within P). Union every
 *     qualifying phylum's contribution across phyla (dedup by rep_id) --
 *     THIS UNION IS THE FINAL SET, no padding to any target count.
 * </ul>
 *
 * <p>Output familylist.tsv contains ONLY the selected set, rank 0..K-1,
 * ordered by occ_total descending (rep_id ascending tiebreak for
 * determinism). K is whatever the union naturally comes out to (expected
 * ~4000, NOT a fixed target). <b>Downstream topn= (CacheBuilder etc.) must be
 * set to the actual K this run reports.</b>
 *
 * <p>sweep= mode: given a comma-separated list of ntop:nbot:minoccfrac
 * triples, computes and prints the resulting union size + per-phylum
 * contribution distribution for EACH combo from a single pass over
 * cluster.tsv (no output files written) -- lets Brian/UMP45 pick a
 * (ntop,nbot,minoccfrac) that lands near a target size without re-scanning
 * the 10.5GB cluster file once per candidate. <b>Sweep mode ranks BOTH halves
 * by occ_total</b> (the original, cheaper metric) -- it is a coarse union-
 * SIZE estimator only, used to pick the (ntop,nbot,minoccfrac) triple before
 * the occ_single refinement was added; a real build (this class run without
 * sweep=) always uses the occ_single-ranked backbone half described above.
 *
 * <p>Usage: familylistbuilder.sh cluster=families_cluster.tsv reps=families_rep_seq.fasta
 *        taxpgm=taxpgm.tsv excluded=EXCLUDED_TIDS.tsv out=familylist.tsv
 *        repsout=consensus_reps.fasta [ntop=100] [nbot=100] [minoccfrac=0.50]
 *        [minphylumgenomes=3] [tidsout=tids.txt]
 *    OR: sweep=100:100:0.50,100:100:0.70,150:150:0.50,150:150:0.70 (report-only, no output files)
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
			else if(a.equals("ntop") || a.equals("topn")){nTop=Integer.parseInt(b);}
			else if(a.equals("nbot") || a.equals("bottomn")){nBot=Integer.parseInt(b);}
			else if(a.equals("minoccfrac") || a.equals("occfloor")){minOccFrac=Double.parseDouble(b);}
			else if(a.equals("minphylumgenomes") || a.equals("minphylum")){minPhylumGenomes=Integer.parseInt(b);}
			else if(a.equals("tidsout") || a.equals("tidsfile")){tidsOutFile=b;}
			else if(a.equals("sweep")){sweepSpec=b;}
			else if(a.equals("ow") || a.equals("overwrite")){overwrite=Parse.parseBoolean(b);}
			else{outstream.println("Unknown parameter "+arg); assert(false) : "Unknown parameter "+arg;}
		}
		assert(clusterFile!=null) : "cluster= is required.";
		assert(repFastaFile!=null) : "reps= is required.";
		assert(taxpgmFile!=null) : "taxpgm= is required.";
		assert(excludedFile!=null) : "excluded= is required.";
		if(sweepSpec==null){assert(outFile!=null) : "out= is required (unless sweep= is used).";}
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

		ArrayList<String> allReps=new ArrayList<String>(familyTids.size());
		ArrayList<Integer> allOccTotal=new ArrayList<Integer>(familyTids.size());
		for(HashMap.Entry<String, HashSet<Integer>> e : familyTids.entrySet()){
			allReps.add(e.getKey());
			allOccTotal.add(e.getValue().size());
		}
		outstream.println("Total families (non-excluded members only): "+allReps.size());

		// Per-phylum candidate lists: family -> occ_in_phylum, for every qualifying
		// phylum it touches. Independent of ntop/nbot/minoccfrac -- computed once,
		// reused by every sweep combo.
		HashMap<String, ArrayList<int[]>> candidatesByPhylum=
			buildCandidatesByPhylum(allReps, familyTids, tid2phylum, phylumGenomeCount);

		if(sweepSpec!=null){
			runSweep(candidatesByPhylum, allReps, allOccTotal, phylumGenomeCount);
			t.stop();
			outstream.println("Time: \t"+t);
			return;
		}

		// BUILD mode: the backbone (TOP_N) half is ranked by occ_single, which
		// needs per-(family,genome) copy counts -- a heavier structure than the
		// presence-only familyTids above. Bound the cost by computing it ONLY
		// for families that actually pass some qualifying phylum's gate (the
		// candidate universe), via a second pass over cluster.tsv.
		Combo combo=new Combo(nTop, nBot, minOccFrac);
		HashMap<String, ArrayList<int[]>> pools=gatePools(combo, candidatesByPhylum, phylumGenomeCount);
		HashSet<String> candidateReps=collectCandidateReps(pools, allReps);
		outstream.println("Candidate universe (pre top/bot truncation, needs occ_single): "
			+candidateReps.size()+" families across "+pools.size()+" qualifying phyla");

		HashMap<String, HashMap<Integer, Integer>> memberCounts=loadCandidateMemberCounts(candidateReps, excluded);
		HashMap<String, Integer> occSingle=deriveOccSingle(memberCounts);
		checkOccSingleCoverage(pools, allReps, occSingle);

		SelectionResult result=selectWithOccSingle(combo, pools, allReps, allOccTotal, occSingle);
		outstream.println("Qualifying phyla contributing: "+result.phylaContributing
			+"; hit full "+(nTop+nBot)+": "+result.phylaFull);
		outstream.println("FINAL selected family count: "+result.selected.size());
		printPerPhylumDistribution(result);

		ArrayList<String> finalReps=new ArrayList<String>(result.selected.size());
		ArrayList<Integer> finalOcc=new ArrayList<Integer>(result.selected.size());
		orderFinal(result.selected, allReps, allOccTotal, finalReps, finalOcc);

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

		outstream.println("occ_single sanity (top of final list, occ_total-desc order, for independent spot-check):");
		int showTop=Tools.min(5, finalReps.size());
		for(int i=0; i<showTop; i++){
			String rep=finalReps.get(i);
			outstream.println("  rank "+i+"\t"+rep+"\tocc_total="+finalOcc.get(i)
				+"\tocc_single="+occSingle.get(rep));
		}

		t.stop();
		outstream.println("Time: \t"+t);
	}

	/** phylum -> list of [repIndex, occ_in_phylum] for every family with >=1
	 *  member in that phylum, restricted to phyla meeting minPhylumGenomes. */
	HashMap<String, ArrayList<int[]>> buildCandidatesByPhylum(ArrayList<String> allReps,
			HashMap<String, HashSet<Integer>> familyTids, HashMap<Integer, String> tid2phylum,
			HashMap<String, Integer> phylumGenomeCount){
		HashMap<String, ArrayList<int[]>> candidatesByPhylum=new HashMap<String, ArrayList<int[]>>();
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
				Integer genomeCount=phylumGenomeCount.get(phy);
				if(genomeCount==null || genomeCount<minPhylumGenomes){continue;}
				ArrayList<int[]> list=candidatesByPhylum.get(phy);
				if(list==null){list=new ArrayList<int[]>(); candidatesByPhylum.put(phy, list);}
				list.add(new int[]{i, e.getValue()});
			}
		}
		return candidatesByPhylum;
	}

	/** Applies one (ntop, nbot, minoccfrac) combo to the precomputed per-phylum
	 *  candidate lists. No file IO -- pure in-memory selection, safe to call
	 *  repeatedly per sweep combo. */
	SelectionResult select(Combo combo, HashMap<String, ArrayList<int[]>> candidatesByPhylum,
			ArrayList<String> allReps, ArrayList<Integer> allOccTotal, HashMap<String, Integer> phylumGenomeCount){
		SelectionResult result=new SelectionResult();
		for(HashMap.Entry<String, ArrayList<int[]>> e : candidatesByPhylum.entrySet()){
			String phy=e.getKey();
			int genomesInP=phylumGenomeCount.get(phy);
			double gate=combo.minOccFrac*genomesInP;

			ArrayList<int[]> pool=new ArrayList<int[]>();
			for(int[] cand : e.getValue()){
				if(cand[1]>=gate){pool.add(cand);}
			}
			if(pool.isEmpty()){continue;}
			result.phylaContributing++;

			ArrayList<int[]> topSorted=new ArrayList<int[]>(pool);
			topSorted.sort((a, b) -> {
				int c=Integer.compare(allOccTotal.get(b[0]), allOccTotal.get(a[0])); // occ_total desc
				if(c!=0){return c;}
				c=Integer.compare(b[1], a[1]); // occ_in_phylum desc
				if(c!=0){return c;}
				return allReps.get(a[0]).compareTo(allReps.get(b[0])); // rep asc
			});
			ArrayList<int[]> botSorted=new ArrayList<int[]>(pool);
			botSorted.sort((a, b) -> {
				int c=Integer.compare(allOccTotal.get(a[0]), allOccTotal.get(b[0])); // occ_total asc
				if(c!=0){return c;}
				c=Integer.compare(b[1], a[1]); // occ_in_phylum desc
				if(c!=0){return c;}
				return allReps.get(a[0]).compareTo(allReps.get(b[0])); // rep asc
			});

			HashSet<String> contribution=new HashSet<String>();
			int takeTop=Tools.min(combo.nTop, topSorted.size());
			for(int i=0; i<takeTop; i++){contribution.add(allReps.get(topSorted.get(i)[0]));}
			int takeBot=Tools.min(combo.nBot, botSorted.size());
			for(int i=0; i<takeBot; i++){contribution.add(allReps.get(botSorted.get(i)[0]));}

			if(contribution.size()>=combo.nTop+combo.nBot){result.phylaFull++;}
			result.selected.addAll(contribution);
			result.perPhylumCount.put(phy, contribution.size());
		}
		return result;
	}

	void runSweep(HashMap<String, ArrayList<int[]>> candidatesByPhylum, ArrayList<String> allReps,
			ArrayList<Integer> allOccTotal, HashMap<String, Integer> phylumGenomeCount){
		outstream.println("SWEEP MODE -- report-only, no output files written.");
		outstream.println("ntop\tnbot\tminoccfrac\tunion_size\tphyla_contributing\tphyla_full");
		for(String token : sweepSpec.split(",")){
			String[] parts=token.split(":");
			assert(parts.length==3) : "sweep combo must be ntop:nbot:minoccfrac, got "+token;
			int nt=Integer.parseInt(parts[0]);
			int nb=Integer.parseInt(parts[1]);
			double frac=Double.parseDouble(parts[2]);
			Combo combo=new Combo(nt, nb, frac);
			SelectionResult result=select(combo, candidatesByPhylum, allReps, allOccTotal, phylumGenomeCount);
			outstream.println(nt+"\t"+nb+"\t"+frac+"\t"+result.selected.size()+"\t"
				+result.phylaContributing+"\t"+result.phylaFull);
		}
	}

	/** Filters each qualifying phylum's candidate list down to the gate
	 *  (occ_in_phylum &gt;= minoccfrac*genomes_in_P), WITHOUT taking top/bot yet
	 *  -- exposes the raw pool so its membership can be unioned into the
	 *  occ_single candidate universe before ranking. */
	HashMap<String, ArrayList<int[]>> gatePools(Combo combo, HashMap<String, ArrayList<int[]>> candidatesByPhylum,
			HashMap<String, Integer> phylumGenomeCount){
		HashMap<String, ArrayList<int[]>> pools=new HashMap<String, ArrayList<int[]>>();
		for(HashMap.Entry<String, ArrayList<int[]>> e : candidatesByPhylum.entrySet()){
			String phy=e.getKey();
			int genomesInP=phylumGenomeCount.get(phy);
			double gate=combo.minOccFrac*genomesInP;
			ArrayList<int[]> pool=new ArrayList<int[]>();
			for(int[] cand : e.getValue()){
				if(cand[1]>=gate){pool.add(cand);}
			}
			if(!pool.isEmpty()){pools.put(phy, pool);}
		}
		return pools;
	}

	static HashSet<String> collectCandidateReps(HashMap<String, ArrayList<int[]>> pools, ArrayList<String> allReps){
		HashSet<String> reps=new HashSet<String>();
		for(ArrayList<int[]> pool : pools.values()){
			for(int[] cand : pool){reps.add(allReps.get(cand[0]));}
		}
		return reps;
	}

	/** Second pass over cluster.tsv, restricted to candidateReps (bounds memory
	 *  -- see class javadoc): rep -> {tid -> member count for that genome}.
	 *  Needed because familyTids/occ_total only recorded PRESENCE, not the
	 *  per-genome copy count occ_single requires. */
	HashMap<String, HashMap<Integer, Integer>> loadCandidateMemberCounts(HashSet<String> candidateReps, HashSet<Integer> excluded){
		outstream.println("Second pass: counting per-(family,genome) copies for "
			+candidateReps.size()+" candidate families...");
		HashMap<String, HashMap<Integer, Integer>> counts=new HashMap<String, HashMap<Integer, Integer>>(candidateReps.size()*2);
		long totalLines=0, matched=0;
		final ByteFile bf=ByteFile.makeByteFile(clusterFile, true);
		final LineParser1 lp=new LineParser1((byte)'\t');
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			lp.set(line);
			if(lp.terms()<2){continue;}
			totalLines++;
			String rep=lp.parseString(0);
			if(candidateReps.contains(rep)){
				String member=lp.parseString(1);
				int tid=parseTid(member);
				if(!excluded.contains(tid)){
					matched++;
					HashMap<Integer, Integer> tidCounts=counts.get(rep);
					if(tidCounts==null){tidCounts=new HashMap<Integer, Integer>(); counts.put(rep, tidCounts);}
					Integer cur=tidCounts.get(tid);
					tidCounts.put(tid, cur==null ? 1 : cur+1);
				}
			}
			if(totalLines%20000000==0){
				outstream.println("  "+totalLines/1000000+"M lines scanned, "+matched+" matched candidate families");
			}
		}
		bf.close();
		outstream.println("Second pass done: "+matched+" matching members, "+counts.size()+" candidate families with data");
		return counts;
	}

	static HashMap<String, Integer> deriveOccSingle(HashMap<String, HashMap<Integer, Integer>> memberCounts){
		HashMap<String, Integer> occSingle=new HashMap<String, Integer>(memberCounts.size()*2);
		for(HashMap.Entry<String, HashMap<Integer, Integer>> e : memberCounts.entrySet()){
			int single=0;
			for(int c : e.getValue().values()){if(c==1){single++;}}
			occSingle.put(e.getKey(), single);
		}
		return occSingle;
	}

	/** Every family in every gated pool MUST have an occ_single entry -- pass 1
	 *  and pass 2 read the identical file, so a miss here means the two passes
	 *  disagreed (a real bug), not a data-coverage gap. */
	static void checkOccSingleCoverage(HashMap<String, ArrayList<int[]>> pools, ArrayList<String> allReps,
			HashMap<String, Integer> occSingle){
		for(ArrayList<int[]> pool : pools.values()){
			for(int[] cand : pool){
				String rep=allReps.get(cand[0]);
				assert(occSingle.containsKey(rep)) : "Candidate family "+rep+" missing from pass-2 "
					+"occ_single map -- pass 1/2 file read mismatch (should be impossible, same file)";
			}
		}
	}

	/** BUILD-mode selection: TOP_N ranked by occ_single desc (tiebreak occ_total
	 *  desc, occ_in_phylum desc, rep asc) -- the clean single-copy backbone.
	 *  BOT_N ranked by occ_total asc (tiebreak occ_in_phylum desc, rep asc),
	 *  UNCHANGED from the original design. Pools are pre-gated (see gatePools). */
	SelectionResult selectWithOccSingle(Combo combo, HashMap<String, ArrayList<int[]>> pools,
			ArrayList<String> allReps, ArrayList<Integer> allOccTotal, HashMap<String, Integer> occSingle){
		SelectionResult result=new SelectionResult();
		for(HashMap.Entry<String, ArrayList<int[]>> e : pools.entrySet()){
			String phy=e.getKey();
			ArrayList<int[]> pool=e.getValue();
			result.phylaContributing++;

			ArrayList<int[]> topSorted=new ArrayList<int[]>(pool);
			topSorted.sort((a, b) -> {
				String repA=allReps.get(a[0]), repB=allReps.get(b[0]);
				int c=Integer.compare(occSingle.get(repB), occSingle.get(repA)); // occ_single desc
				if(c!=0){return c;}
				c=Integer.compare(allOccTotal.get(b[0]), allOccTotal.get(a[0])); // occ_total desc
				if(c!=0){return c;}
				c=Integer.compare(b[1], a[1]); // occ_in_phylum desc
				if(c!=0){return c;}
				return repA.compareTo(repB); // rep asc
			});
			ArrayList<int[]> botSorted=new ArrayList<int[]>(pool);
			botSorted.sort((a, b) -> {
				int c=Integer.compare(allOccTotal.get(a[0]), allOccTotal.get(b[0])); // occ_total asc
				if(c!=0){return c;}
				c=Integer.compare(b[1], a[1]); // occ_in_phylum desc
				if(c!=0){return c;}
				return allReps.get(a[0]).compareTo(allReps.get(b[0])); // rep asc
			});

			HashSet<String> contribution=new HashSet<String>();
			int takeTop=Tools.min(combo.nTop, topSorted.size());
			for(int i=0; i<takeTop; i++){contribution.add(allReps.get(topSorted.get(i)[0]));}
			int takeBot=Tools.min(combo.nBot, botSorted.size());
			for(int i=0; i<takeBot; i++){contribution.add(allReps.get(botSorted.get(i)[0]));}

			if(contribution.size()>=combo.nTop+combo.nBot){result.phylaFull++;}
			result.selected.addAll(contribution);
			result.perPhylumCount.put(phy, contribution.size());
		}
		return result;
	}

	/** min/median/max markers-per-phylum + the most-truncated phyla -- the
	 *  starvation check UMP45 gates a high floor on (a phylum with very few
	 *  contributed markers gives a coarse completeness estimate). */
	void printPerPhylumDistribution(SelectionResult result){
		ArrayList<HashMap.Entry<String, Integer>> entries=
			new ArrayList<HashMap.Entry<String, Integer>>(result.perPhylumCount.entrySet());
		entries.sort((a, b) -> Integer.compare(a.getValue(), b.getValue()));
		int n=entries.size();
		outstream.println("Per-phylum contribution distribution ("+n+" phyla):");
		if(n==0){outstream.println("  (no contributing phyla)"); return;}
		int min=entries.get(0).getValue();
		int max=entries.get(n-1).getValue();
		double median=(n%2==1) ? entries.get(n/2).getValue()
			: (entries.get(n/2-1).getValue()+entries.get(n/2).getValue())/2.0;
		outstream.println("  min="+min+" median="+median+" max="+max);
		int showN=Tools.min(10, n);
		outstream.println("  "+showN+" most-truncated phyla:");
		for(int i=0; i<showN; i++){
			outstream.println("    "+entries.get(i).getKey()+"\t"+entries.get(i).getValue());
		}
	}

	void orderFinal(HashSet<String> selected, ArrayList<String> allReps, ArrayList<Integer> allOccTotal,
			ArrayList<String> finalReps, ArrayList<Integer> finalOcc){
		Integer[] order=new Integer[allReps.size()];
		for(int i=0; i<order.length; i++){order[i]=i;}
		Arrays.sort(order, (a, b) -> {
			int c=Integer.compare(allOccTotal.get(b), allOccTotal.get(a)); // occ_total desc
			if(c!=0){return c;}
			return allReps.get(a).compareTo(allReps.get(b)); // rep asc
		});
		for(int idx : order){
			String rep=allReps.get(idx);
			if(selected.contains(rep)){
				finalReps.add(rep);
				finalOcc.add(allOccTotal.get(idx));
			}
		}
		assert(finalReps.size()==selected.size());
	}

	static class Combo {
		Combo(int nTop_, int nBot_, double minOccFrac_){nTop=nTop_; nBot=nBot_; minOccFrac=minOccFrac_;}
		final int nTop, nBot;
		final double minOccFrac;
	}

	static class SelectionResult {
		HashSet<String> selected=new HashSet<String>();
		int phylaContributing=0;
		int phylaFull=0;
		/** phylum -> count of families it contributed to the union (its TOP_N/BOT_N
		 *  overlap-deduped size) -- for the min/median starvation-check UMP45 gates on. */
		HashMap<String, Integer> perPhylumCount=new HashMap<String, Integer>();
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
	private String sweepSpec=null;
	private int nTop=100;
	private int nBot=100;
	private double minOccFrac=0.50;
	private int minPhylumGenomes=3;
	private boolean overwrite=true;
	private java.io.PrintStream outstream=System.err;
}
