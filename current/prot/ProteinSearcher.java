package prot;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import dna.AminoAcid;

/**
 * In-memory, blastp-style protein similarity search (Phase-1 MVP).
 *
 * <p><b>This is the load-bearing API.</b> {@link #search(List, List)} takes
 * query and database protein sequences as in-memory objects and returns a list
 * of {@link ProteinHit} objects — no disk round-trip. It is callable from other
 * BBTools code (e.g. so a binner could later assess an in-memory bin against a
 * marker set).</p>
 *
 * <p>Algorithm (basic version, correctness-first):</p>
 * <ul>
 * <li><b>K-mer seed:</b> exact (or amino8-reduced) k-mer overlap selects
 * candidate targets for each query — a target is a candidate when it shares at
 * least {@link #minSeedHits} distinct k-mers with the query.
 * <li><b>Extend/score:</b> each candidate is aligned with BLOSUM62 affine-gap
 * local Smith-Waterman ({@link AAAligner}), yielding the single best HSP per
 * query-target pair.
 * <li><b>Statistics:</b> rigorous gapped-BLOSUM62 bitscore, plus an approximate
 * Karlin-Altschul E-value (no edge-length correction; flagged approximate).
 * <li><b>Filter + sort:</b> E-value cutoff and optional min raw score / min
 * pident; {@code maxTargetSeqs} caps distinct targets by best-HSP score; output
 * is sorted by the frozen total order (query, E-value, bitscore, target, tstart,
 * qstart) so results are deterministic.
 * </ul>
 *
 * <p>Deferred (reported, not built): rigorous edge-corrected E-values, indexed /
 * multithreaded fast search, multi-HSP per pair, and blastx / six-frame
 * translated search.</p>
 *
 * @author Eru
 */
public final class ProteinSearcher {

	/** K-mer length for seeding. */
	public int k=5;
	/** Minimum distinct shared k-mers for a target to be a candidate. */
	public int minSeedHits=1;
	/** E-value significance cutoff (frozen default 10). */
	public double evalueCutoff=10.0;
	/** Minimum raw score to report (0 = disabled). */
	public int minRawScore=0;
	/** Minimum percent identity to report (0 = disabled). */
	public double minPident=0;
	/** Cap on distinct targets per query, ranked by best-HSP score. */
	public int maxTargetSeqs=Integer.MAX_VALUE;
	/** Use the amino8 reduced alphabet for seeding (more sensitive). */
	public boolean reducedSeed=false;
	/**
	 * True: the emitted E-values omit the KA edge-length correction and are
	 * therefore approximate (bitscores remain rigorous). Fixed true for the MVP.
	 */
	public static final boolean EVALUE_APPROXIMATE=true;

	/** Maps extended-alphabet codes (0-19) to amino8 reduced groups (0-6). */
	private static final int[] EXT_TO_AMINO8=buildReducedMap();

	private static int[] buildReducedMap(){
		final int[] map=new int[Blosum62.X_CODE+1];
		java.util.Arrays.fill(map, -1);
		for(char c='A'; c<='Z'; c++){
			final int ext=AminoAcid.acidToNumber[c];
			if(ext>=0 && ext<=19){map[ext]=AminoAcid.acidToNumber8[c];}
		}
		return map;
	}

	/**
	 * Searches every query against the database and returns all passing hits,
	 * sorted in the frozen deterministic total order.
	 *
	 * @param queries In-memory query protein sequences.
	 * @param targets In-memory database protein sequences.
	 * @return Hits passing the significance filters, deterministically sorted.
	 */
	public List<ProteinHit> search(final List<ProteinSequence> queries,
			final List<ProteinSequence> targets){
		if(queries==null || targets==null){
			throw new RuntimeException("Null query or target list.");
		}
		checkDuplicateIds(queries, "query");
		checkDuplicateIds(targets, "database");

		long totalDbResidues=0;
		for(ProteinSequence t : targets){totalDbResidues+=t.length();}

		//Build the target k-mer index: kmer -> distinct target indices.
		final HashMap<Long, ArrayList<Integer>> index=buildIndex(targets);

		final ArrayList<ProteinHit> all=new ArrayList<ProteinHit>();
		final int[] hitCount=new int[targets.size()];
		final ArrayList<Integer> touched=new ArrayList<Integer>();

		for(ProteinSequence q : queries){
			final double searchSpace=(double)q.length()*(double)totalDbResidues;

			//Seed: tally shared k-mers per target.
			for(Integer idx : touched){hitCount[idx]=0;}
			touched.clear();
			final HashSet<Long> qkmers=kmerSet(q.enc);
			if(qkmers.isEmpty()){
				//Query too short for a k-mer: fall back to aligning against all targets.
				for(int i=0; i<targets.size(); i++){hitCount[i]=minSeedHits;}
			}else{
				for(Long km : qkmers){
					final ArrayList<Integer> list=index.get(km);
					if(list!=null){
						for(Integer idx : list){
							if(hitCount[idx]==0){touched.add(idx);}
							hitCount[idx]++;
						}
					}
				}
			}

			//Extend + score every candidate target.
			final ArrayList<ProteinHit> qhits=new ArrayList<ProteinHit>();
			final int nTargets=targets.size();
			for(int i=0; i<nTargets; i++){
				if(qkmers.isEmpty() ? true : hitCount[i]>=minSeedHits){
					final ProteinSequence t=targets.get(i);
					final AAAlignment aln=AAAligner.align(q.enc, t.enc);
					if(aln==null){continue;}
					if(aln.rawScore<minRawScore){continue;}
					final double pid=aln.pident();
					if(pid<minPident){continue;}
					final double e=aln.evalue(searchSpace);
					if(e>evalueCutoff){continue;}
					qhits.add(new ProteinHit(q.id, t.id, aln, e, EVALUE_APPROXIMATE));
				}
			}

			//maxTargetSeqs: keep top-N distinct targets by best-HSP score (bitscore).
			if(maxTargetSeqs<qhits.size()){
				Collections.sort(qhits, BY_SCORE_DESC);
				while(qhits.size()>maxTargetSeqs){qhits.remove(qhits.size()-1);}
			}
			all.addAll(qhits);
		}

		Collections.sort(all, TOTAL_ORDER);
		return all;
	}

	/**
	 * Convenience: single-pair search returning the best HSP or null.
	 * @param q Query sequence.
	 * @param t Target sequence.
	 * @param searchSpace Effective search space for the E-value.
	 * @return Best hit, or null if none scores positively.
	 */
	public ProteinHit searchPair(final ProteinSequence q, final ProteinSequence t,
			final double searchSpace){
		final AAAlignment aln=AAAligner.align(q.enc, t.enc);
		if(aln==null){return null;}
		return new ProteinHit(q.id, t.id, aln, aln.evalue(searchSpace), EVALUE_APPROXIMATE);
	}

	/** Builds the target k-mer index (kmer packed as a long -> target indices). */
	private HashMap<Long, ArrayList<Integer>> buildIndex(final List<ProteinSequence> targets){
		final HashMap<Long, ArrayList<Integer>> index=new HashMap<Long, ArrayList<Integer>>();
		for(int i=0; i<targets.size(); i++){
			final Integer idx=Integer.valueOf(i);
			for(Long km : kmerSet(targets.get(i).enc)){
				ArrayList<Integer> list=index.get(km);
				if(list==null){list=new ArrayList<Integer>(); index.put(km, list);}
				list.add(idx);//distinct per target, so no duplicate index per kmer
			}
		}
		return index;
	}

	/**
	 * Extracts the set of distinct k-mers from an encoded sequence. K-mers
	 * containing X (ambiguous) are skipped. Uses the amino8 reduced alphabet
	 * when {@link #reducedSeed} is set.
	 * @param enc Encoded residues.
	 * @return Set of packed k-mers.
	 */
	private HashSet<Long> kmerSet(final byte[] enc){
		final HashSet<Long> set=new HashSet<Long>();
		final int bits=reducedSeed ? 3 : 5;
		assert(bits*k<=62) : "k too large for packed k-mer: k="+k+", bits="+bits;
		if(enc.length<k){return set;}
		final long mask=(1L<<(bits*k))-1;
		long kmer=0;
		int valid=0;
		for(int i=0; i<enc.length; i++){
			final byte e=enc[i];
			int code;
			if(reducedSeed){
				code=(e>=0 && e<EXT_TO_AMINO8.length) ? EXT_TO_AMINO8[e] : -1;
			}else{
				code=Blosum62.isStandard(e) ? e : -1;
			}
			if(code<0){kmer=0; valid=0; continue;}//reset on X/ambiguous residue
			kmer=((kmer<<bits)|code)&mask;
			valid++;
			if(valid>=k){set.add(Long.valueOf(kmer));}
		}
		return set;
	}

	/** Fails loudly on duplicate identifiers within one input (frozen contract). */
	private static void checkDuplicateIds(final List<ProteinSequence> seqs, final String which){
		final HashSet<String> seen=new HashSet<String>();
		for(ProteinSequence s : seqs){
			if(!seen.add(s.id)){
				throw new RuntimeException("Duplicate "+which+" identifier: '"+s.id+"'.");
			}
		}
	}

	/** Orders hits by best-HSP score descending (for maxTargetSeqs culling). */
	private static final Comparator<ProteinHit> BY_SCORE_DESC=new Comparator<ProteinHit>(){
		@Override
		public int compare(ProteinHit a, ProteinHit b){
			if(a.bitscore!=b.bitscore){return a.bitscore>b.bitscore ? -1 : 1;}
			return a.target.compareTo(b.target);
		}
	};

	/**
	 * The frozen total output order: query asc, E-value asc, bitscore desc,
	 * target asc, tstart asc, qstart asc. Deterministic (a total order).
	 */
	private static final Comparator<ProteinHit> TOTAL_ORDER=new Comparator<ProteinHit>(){
		@Override
		public int compare(ProteinHit a, ProteinHit b){
			int c=a.query.compareTo(b.query);
			if(c!=0){return c;}
			if(a.evalue!=b.evalue){return a.evalue<b.evalue ? -1 : 1;}
			if(a.bitscore!=b.bitscore){return a.bitscore>b.bitscore ? -1 : 1;}
			c=a.target.compareTo(b.target);
			if(c!=0){return c;}
			if(a.tstart!=b.tstart){return a.tstart<b.tstart ? -1 : 1;}
			return Integer.compare(a.qstart, b.qstart);
		}
	};
}
