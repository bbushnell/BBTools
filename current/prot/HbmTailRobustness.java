package prot;

import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Direct known-family tail diagnostics with reproducible artificial endings.
 * This measures alignment robustness, not biological precision or retrieval recall.
 * Inputs are encoded by the ordinary protein reader; no caller scores are changed.
 * Replacement names specify maximum spans: replace40 changes the entire query
 * when it is shorter than forty residues. Length remains explicit in every row.
 * @author Keqing
 */
public final class HbmTailRobustness{

	public static void main(String[] args){
		try{run(args);}catch(Throwable failure){failure.printStackTrace(System.err); System.exit(1);}
	}

	private static void run(String[] args) throws Exception{
		final HashMap<String,String> opt=HmmComparisonData.options(args);
		final HashSet<String> allowed=new HashSet<String>(Arrays.asList("reference", "in", "families", "out", "count"));
		for(String key : opt.keySet()){require(allowed.contains(key), "Unknown diagnostic option: "+key);}
		require(opt.containsKey("in")==opt.containsKey("families"), "Real queries require both in=FASTA and families=query/family_id TSV");
		final ProteinFamilyAssigner pack=ProteinFamilyAssigner.load(Paths.get(HmmComparisonData.required(opt, "reference")));
		final HbmCompetitiveAssign.Library library=pack.library();
		final ShortlistProbe shortlist=new ShortlistProbe(library.sidecar);
		final int[][] tallies=new int[VARIANTS.length][3];
		final ByteStreamWriter out=HmmComparisonData.writer(HmmComparisonData.required(opt, "out"));
		out.println("query\tquery_source\tfamily_id\tvariant\tlength\tstatus\tscore64\tidentity\tpaired_fraction\tcore_coverage\tpaired\tpaired_core\tterminal_paired\tterminal_identical\tquery_terminal_paired\tquery_terminal_identical\tquery_suffix\tmodel_suffix\tquery_start\tmodel_start\tmodel_length\tpath\tshortlist_rank\tterminal20_score64\tterminal40_score64\tpaired_identity\tcounterfactual_pass");
		int count=0;
		if(opt.containsKey("in")){
			final HashMap<String,Integer> families=new HashMap<String,Integer>();
			final List<String[]> rows=HmmComparisonData.rows(opt.get("families"));
			require(!rows.isEmpty() && Arrays.equals(rows.get(0), new String[]{"query", "family_id"}), "Expected query/family_id mapping header");
			for(int i=1; i<rows.size(); i++){
				final String[] row=rows.get(i);
				require(row.length==2 && families.put(row[0], Integer.parseInt(row[1]))==null, "Malformed or duplicate query mapping at row "+i);
			}
			for(ProteinSequence p : ProteinSearch.readFasta(opt.get("in"))){
				final Integer family=families.remove(p.id);
				require(family!=null, "Missing or repeated mapping for query "+p.id);
				int rank=-1;
				for(int i=0; i<library.familyIds.length; i++){if(library.familyIds[i]==family){rank=i; break;}}
				require(rank>=0, "Family absent from the frozen reference: "+family);
				emitVariants(out, library, shortlist, tallies, p.id, "input", rank, p.enc); count++;
			}
			require(families.isEmpty(), "Mapping contains queries missing from FASTA");
		}else{
			final int n=Integer.parseInt(opt.getOrDefault("count", "200"));
			require(n>0 && n<=pack.familyCount(), "Consensus diagnostic count must fit the reference");
			for(int i=0; i<n; i++){
				final int rank=(int)((long)i*pack.familyCount()/n);
				emitVariants(out, library, shortlist, tallies, library.ids[rank], "consensus_self", rank, library.refs[rank]); count++;
			}
		}
		HmmComparisonData.close(out);
		require(count>0, "No diagnostic queries were evaluated");
		for(int v=0; v<VARIANTS.length; v++){
			require(tallies[v][0]==count, "Every query must contribute each fixed perturbation exactly once");
			System.err.println("TAIL_POLICY_DIAGNOSTIC variant="+VARIANTS[v]+" queries="+count
				+" current_pass="+tallies[v][1]+" counterfactual_pass="+tallies[v][2]);
		}
		System.err.println("HBM_TAIL_ROBUSTNESS_PASS queries="+count+" variants="+VARIANTS.length+" direct_known_family=true biological_accuracy=false");
	}

	/** Fixed seed is reset per query so execution order cannot alter perturbations. */
	private static void emitVariants(ByteStreamWriter out, HbmCompetitiveAssign.Library library, ShortlistProbe shortlist, int[][] tallies,
			String id, String source, int rank, byte[] original){
		final byte[] base=Arrays.copyOfRange(original, Math.max(0, original.length-POLICY.maxLength), original.length);
		require(base.length>20, "Diagnostic needs over twenty residues so trim20 remains nonempty: "+id);
		for(int v=0; v<VARIANTS.length; v++){
			final byte[] query=variant(base, v);
			final HbmPositionModel.TailResult a=library.gate.alignTail(query, rank, POLICY.maxLength, true);
			final HbmTailGate.Metrics m=new HbmTailGate.Metrics();
			final int status=library.gate.evaluateTail(query, rank, a, POLICY, m);
			int jointIdentical=0, queryPaired=0, queryIdentical=0, identities=0;
			if(a!=null){
				int q=a.queryStart, r=a.start;
				for(byte op : a.path){
					if(op=='m'){
						final boolean identical=query[q]<20 && query[q]==library.refs[rank][r];
						if(identical){identities++;}
						if(q>=query.length-POLICY.terminalWindow){
							queryPaired++; if(identical){queryIdentical++;}
							if(r>=library.refs[rank].length-POLICY.terminalWindow && identical){jointIdentical++;}
						}
						q++; r++;
					}else if(op=='I'){q++;}else if(op=='D'){r++;}
					else{throw new IllegalArgumentException("Invalid tail trace operation");}
				}
				require(q==query.length && r==library.refs[rank].length, "Terminal diagnostic must consume both declared endpoints");
			}
			final int score20=a==null ? Integer.MIN_VALUE : library.gate.tailWindowScore(query, rank, a, 20);
			final int score40=a==null ? Integer.MIN_VALUE : library.gate.tailWindowScore(query, rank, a, 40);
			final float pairedIdentity=m.paired==0 ? Float.NaN : (float)(100.0*identities/m.paired);
			// Counterfactual only: same alignment, coverage and suffix bounds; replace
			// whole-path identity/joint-terminal pairing with paired identity and query-end signal.
			final boolean candidate=a!=null && m.paired>=POLICY.minPaired && m.pairedCore>=POLICY.minPaired
				&& pairedIdentity>=POLICY.minIdentity && m.pairedFraction>=POLICY.minCoverage && m.coreCoverage>=POLICY.minCoverage
				&& m.querySuffix<=POLICY.maxTerminalGap && m.referenceSuffix<=POLICY.maxTerminalGap
				&& queryPaired>=POLICY.minTerminalPaired && score20>0;
			tallies[v][0]++; if(status==HbmTailGate.SUPPORTED){tallies[v][1]++;} if(candidate){tallies[v][2]++;}
			out.println(new ByteBuilder().append(id).tab().append(source).tab().append(library.familyIds[rank])
				.tab().append(VARIANTS[v]).tab().append(query.length).tab().append(HbmTailGate.statusName(status))
				.tab().append(a==null ? 0 : a.score).tab().append(Float.toString(m.identity))
				.tab().append(Float.toString(m.pairedFraction)).tab().append(Float.toString(m.coreCoverage))
				.tab().append(m.paired).tab().append(m.pairedCore).tab().append(m.terminalPaired).tab().append(jointIdentical)
				.tab().append(queryPaired).tab().append(queryIdentical).tab().append(m.querySuffix).tab().append(m.referenceSuffix)
				.tab().append(m.queryStart).tab().append(m.referenceStart).tab().append(library.refs[rank].length)
				.tab().append(a==null ? new byte[0] : a.path).tab().append(shortlist.rank(id, query, rank))
				.tab().append(a==null ? "NA" : Integer.toString(score20)).tab().append(a==null ? "NA" : Integer.toString(score40))
				.tab().append(Float.toString(pairedIdentity)).tab().append(candidate));
		}
	}

	/** Replacements change up to the named number of residues; shuffling preserves composition. */
	public static byte[] variant(byte[] base, int v){
		assert(base.length>20 && v>=0 && v<VARIANTS.length) : "Fixed trimming variants require a nonempty remaining query";
		if(v==0){return base.clone();}
		if(v==5 || v==6){return Arrays.copyOf(base, base.length-(v==5 ? 10 : 20));}
		final Random random=new Random(20261018L+v);
		if(v==7 || v==8){
			final int n=v==7 ? 10 : 20;
			final byte[] expanded=Arrays.copyOf(base, base.length+n);
			for(int i=base.length; i<expanded.length; i++){expanded[i]=(byte)random.nextInt(20);}
			return Arrays.copyOfRange(expanded, Math.max(0, expanded.length-POLICY.maxLength), expanded.length);
		}
		final byte[] out=base.clone();
		final int n=Math.min(base.length, v==1 ? 10 : v==3 ? 40 : 20);
		if(v==4){
			for(int i=out.length-1; i>out.length-n; i--){
				final int j=out.length-n+random.nextInt(i-(out.length-n)+1);
				final byte b=out[i]; out[i]=out[j]; out[j]=b;
			}
		}else{
			for(int i=out.length-n; i<out.length; i++){out[i]=(byte)((out[i]+1+random.nextInt(19))%20);}
		}
		return out;
	}
	/** Single-pack full rank with exactly the caller's global tail normalization and tie order. */
	private static final class ShortlistProbe{
		ShortlistProbe(FamilyShortlistSidecar sidecar_){
			sidecar=sidecar_;
			scratch=new ProteinSearcher.ShortlistScratch(sidecar.dims, sidecar.nFamilies, 50);
			hits=new int[][]{scratch.sk}; composition=new float[][]{scratch.compScores}; scores=new float[][]{scratch.f4};
			medians=new double[][]{sidecar.medLen}; spreads=new double[][]{sidecar.sdLnLen};
		}
		int rank(String id, byte[] encoded, int target){
			final ProteinSequence query=new ProteinSequence(id, HmmComparisonData.ascii(encoded));
			ProteinSearcher.scoreFamiliesF4(query, sidecar, scratch);
			FamilyShortlistScores.scoreTail(hits, composition, medians, spreads, query.length(), scores);
			int rank=1;
			for(int i=0; i<scratch.f4.length; i++){
				if(scratch.f4[i]>scratch.f4[target] || (scratch.f4[i]==scratch.f4[target] && i<target)){rank++;}
			}
			assert(rank>=1 && rank<=sidecar.nFamilies) : "Exactly one rank per bound family, including canonical ties";
			return rank;
		}
		final FamilyShortlistSidecar sidecar;
		final ProteinSearcher.ShortlistScratch scratch;
		final int[][] hits;
		final float[][] composition, scores;
		final double[][] medians, spreads;
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final HbmTailGate.Policy POLICY=HbmTailGate.EXPERIMENTAL;
	private static final String[] VARIANTS={"original", "replace10", "replace20", "replace40", "shuffle20", "trim10", "trim20", "append10", "append20"};
	private HbmTailRobustness(){}
}
