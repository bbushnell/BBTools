package prot;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Consensus-reassignment protein family clustering: Neptune's 4-step tRNA
 * algorithm ported to amino-acid space with a real {@link AAGraph} consensus.
 *
 * <ol>
 * <li><b>Greedy seed</b> ({@link ProteinClusterer}) -- a k-mer-seeded CD-HIT pass
 * gives an order-sensitive initial partition.</li>
 * <li><b>Consensus centroid</b> ({@link AAGraph#traverse}) -- each cluster's
 * members are piled into a graph and traversed, yielding a consensus that
 * includes insertions and (optionally) is pulled toward divergent members.</li>
 * <li><b>Global reassign</b> ({@link ProteinSearcher}) -- every input sequence is
 * re-placed against the consensus of every cluster. This is what makes the result
 * order-independent: the greedy order stops mattering once each sequence picks its
 * globally-best consensus. Sequences with no consensus above threshold orphan.</li>
 * <li><b>Recruit orphans</b> -- consensuses are rebuilt from the reassigned
 * groups and orphans get a second placement pass; whatever still fails to place
 * becomes its own singleton family.</li>
 * </ol>
 *
 * <p>Thresholds default to the amino-acid family band (identity ~30%, shorter
 * reduced-alphabet seeds), not the near-duplicate regime of the greedy clusterer.
 * Correctness-first MVP: single-threaded, consensuses rebuilt rather than
 * incrementally updated. A speed pass is deferred.</p>
 *
 * @author Eru
 */
public final class ProteinFamilyClusterer {

	/*--------------------------------------------------------------*/
	/*----------------            Config            ----------------*/
	/*--------------------------------------------------------------*/

	/** Greedy Step-1 join identity, percent [0,100]. */
	public double clusterIdentity=30.0;
	/** Aligned-fraction floor for both greedy joins and reassignment [0,1]. */
	public double minCoverage=0.7;
	/** Reassign/recruit identity floor, percent [0,100]. */
	public double reassignIdentity=30.0;
	/** X-padding on each consensus pivot end (allows the consensus to grow). */
	public int pad=20;
	/** Consensus refinement passes (rebuild using the prior consensus as pivot). */
	public int consensusPasses=1;
	/** Pass the identity-inverted weighting through to each consensus graph. */
	public boolean weightByIdentity=false;
	/** Identity ceiling for identity-inverted weighting (percent). */
	public float identityCeiling=40f;
	/** Consensus substitution / deletion / insertion allele-fraction gates. */
	public float MAF_sub=0.25f, MAF_del=0.5f, MAF_ins=0.5f;
	/** Trim consensus ends below this fraction of max depth (strips smeared low-depth padding tails). */
	public float trimDepthFraction=0.1f;
	/** K-mer length for seeding (shorter than nucleotide; AA space is less specific). */
	public int k=4;
	/** Use the amino8 reduced alphabet for seeding (more sensitive at low identity). */
	public boolean reducedSeed=true;
	/** Emit per-step progress to {@link #outstream}. */
	public boolean verbose=false;
	/** Progress stream. */
	public PrintStream outstream=System.err;

	/*--------------------------------------------------------------*/
	/*----------------           Pipeline           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Clusters the input sequences into families via the 4-step algorithm.
	 * @param seqs In-memory protein sequences (unique ids).
	 * @return Families, each with a consensus representative and its members.
	 */
	public List<ProteinFamily> cluster(final List<ProteinSequence> seqs){
		if(seqs==null){throw new RuntimeException("Null sequence list.");}

		//Step 1: greedy seed.
		final ProteinClusterer greedy=new ProteinClusterer();
		greedy.minIdentity=clusterIdentity;
		greedy.minCoverage=minCoverage;
		greedy.k=k;
		greedy.reducedSeed=reducedSeed;
		final List<ProteinCluster> clusters=greedy.cluster(seqs);
		if(verbose){log("Step 1 greedy: "+clusters.size()+" clusters from "+seqs.size()+" seqs.");}

		//Step 2: a consensus centroid per greedy cluster.
		final ArrayList<ProteinSequence> consensus=new ArrayList<ProteinSequence>();
		for(int i=0; i<clusters.size(); i++){
			consensus.add(buildConsensus(membersOf(clusters.get(i)), i));
		}

		//Step 3: global reassignment against those consensuses.
		final HashMap<Integer, ArrayList<ProteinSequence>> groups=
				new HashMap<Integer, ArrayList<ProteinSequence>>();
		ArrayList<ProteinSequence> orphans=new ArrayList<ProteinSequence>();
		assignBySearch(seqs, consensus, groups, orphans);
		if(verbose){log("Step 3 reassign: "+groups.size()+" groups, "+orphans.size()+" orphans.");}

		//Rebuild consensus from the reassigned groups (stable positions in memberLists).
		final ArrayList<ProteinSequence> rebuilt=new ArrayList<ProteinSequence>();
		final ArrayList<ArrayList<ProteinSequence>> memberLists=
				new ArrayList<ArrayList<ProteinSequence>>();
		for(Map.Entry<Integer, ArrayList<ProteinSequence>> e : groups.entrySet()){
			final ArrayList<ProteinSequence> mem=e.getValue();
			rebuilt.add(buildConsensus(mem, memberLists.size()));
			memberLists.add(mem);
		}

		//Step 4: recruit orphans against the rebuilt consensuses.
		if(!orphans.isEmpty() && !rebuilt.isEmpty()){
			final HashMap<Integer, ArrayList<ProteinSequence>> recruitGroups=
					new HashMap<Integer, ArrayList<ProteinSequence>>();
			final ArrayList<ProteinSequence> stillOrphan=new ArrayList<ProteinSequence>();
			assignBySearch(orphans, rebuilt, recruitGroups, stillOrphan);
			int recruited=0;
			for(Map.Entry<Integer, ArrayList<ProteinSequence>> e : recruitGroups.entrySet()){
				memberLists.get(e.getKey().intValue()).addAll(e.getValue());
				recruited+=e.getValue().size();
			}
			orphans=stillOrphan;
			if(verbose){log("Step 4 recruit: "+recruited+" recruited, "+orphans.size()+" remain.");}
		}

		//Assemble: a final consensus per group, plus each leftover orphan as a singleton.
		final ArrayList<ProteinFamily> families=new ArrayList<ProteinFamily>();
		for(ArrayList<ProteinSequence> mem : memberLists){
			final int id=families.size();
			families.add(new ProteinFamily(id, buildConsensus(mem, id), mem));
		}
		for(ProteinSequence o : orphans){
			final int id=families.size();
			final ArrayList<ProteinSequence> mem=new ArrayList<ProteinSequence>();
			mem.add(o);
			families.add(new ProteinFamily(id, buildConsensus(mem, id), mem));
		}
		if(verbose){log("Result: "+families.size()+" families.");}
		return families;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Consensus           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Builds a consensus sequence for a set of members via an {@link AAGraph}.
	 * The longest member is the initial pivot; optional extra passes re-pivot on the
	 * prior consensus. Falls back to the pivot if traversal trims to nothing.
	 *
	 * @param members Members to summarize (at least one).
	 * @param famId Id used for the consensus sequence's identifier.
	 * @return The consensus as a ProteinSequence ("fam"+famId).
	 */
	private ProteinSequence buildConsensus(final ArrayList<ProteinSequence> members, final int famId){
		ProteinSequence pivot=members.get(0);
		for(ProteinSequence s : members){if(s.length()>pivot.length()){pivot=s;}}
		byte[] cons=graphConsensus(pivot.enc, members);
		for(int p=1; p<consensusPasses && cons.length>0; p++){cons=graphConsensus(cons, members);}
		if(cons.length==0){cons=pivot.enc;}//never emit an empty consensus
		return new ProteinSequence("fam"+famId, AAGraph.decode(cons));
	}

	/** Piles every member onto an X-padded pivot graph and returns the trace. */
	private byte[] graphConsensus(final byte[] pivotEnc, final ArrayList<ProteinSequence> members){
		final AAGraph g=new AAGraph(pivotEnc, pad);
		g.weightByIdentity=weightByIdentity;
		g.identityCeiling=identityCeiling;
		g.MAF_sub=MAF_sub; g.MAF_del=MAF_del; g.MAF_ins=MAF_ins;
		g.trimDepthFraction=trimDepthFraction;
		for(ProteinSequence s : members){
			final AAAlignment aln=AAAligner.alignGlocal(s.enc, g.pivot, true);
			if(aln!=null){g.add(s.enc, aln);}
		}
		return g.traverse();
	}

	/*--------------------------------------------------------------*/
	/*----------------         Assignment           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Places each query on its best consensus target (by bitscore, subject to the
	 * identity and coverage floors), filling {@code groups} keyed by target index
	 * and collecting unplaceable queries into {@code orphans}.
	 *
	 * @param queries Sequences to place.
	 * @param targets Consensus centroids (unique ids); index = group key.
	 * @param groups Output map: target index -> assigned queries.
	 * @param orphans Output list: queries with no target above threshold.
	 */
	private void assignBySearch(final List<ProteinSequence> queries,
			final List<ProteinSequence> targets,
			final HashMap<Integer, ArrayList<ProteinSequence>> groups,
			final ArrayList<ProteinSequence> orphans){
		if(targets.isEmpty() || queries.isEmpty()){orphans.addAll(queries); return;}

		final HashMap<String, Integer> tIndex=new HashMap<String, Integer>();
		for(int i=0; i<targets.size(); i++){tIndex.put(targets.get(i).id, Integer.valueOf(i));}

		final ProteinSearcher ps=new ProteinSearcher();
		ps.k=k;
		ps.reducedSeed=reducedSeed;
		ps.minPident=reassignIdentity;
		ps.evalueCutoff=Double.MAX_VALUE;//rely on identity + coverage, not significance
		final List<ProteinHit> hits=ps.search(queries, targets);

		//Best hit per query that also clears the coverage floor.
		final HashMap<String, ProteinSequence> qById=byId(queries);
		final HashMap<String, ProteinHit> best=new HashMap<String, ProteinHit>();
		for(ProteinHit h : hits){
			final ProteinSequence q=qById.get(h.query);
			final double cov=(h.qend-h.qstart+1)/(double)q.length();
			if(cov<minCoverage){continue;}
			final ProteinHit cur=best.get(h.query);
			if(cur==null || h.bitscore>cur.bitscore){best.put(h.query, h);}
		}

		for(ProteinSequence q : queries){
			final ProteinHit h=best.get(q.id);
			if(h==null){orphans.add(q);}
			else{
				final Integer idx=tIndex.get(h.target);
				ArrayList<ProteinSequence> list=groups.get(idx);
				if(list==null){list=new ArrayList<ProteinSequence>(); groups.put(idx, list);}
				list.add(q);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	/** Extracts the member sequences from a greedy cluster. */
	private static ArrayList<ProteinSequence> membersOf(final ProteinCluster c){
		final ArrayList<ProteinSequence> list=new ArrayList<ProteinSequence>(c.members.size());
		for(ClusterMember m : c.members){list.add(m.seq);}
		return list;
	}

	/** Indexes sequences by id for coverage lookups. */
	private static HashMap<String, ProteinSequence> byId(final List<ProteinSequence> seqs){
		final HashMap<String, ProteinSequence> map=new HashMap<String, ProteinSequence>();
		for(ProteinSequence s : seqs){map.put(s.id, s);}
		return map;
	}

	private void log(final String s){outstream.println(s);}
}
