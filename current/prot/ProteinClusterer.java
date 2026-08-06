package prot;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import dna.AminoAcid;

/**
 * In-memory greedy identity-threshold protein clustering (CD-HIT / linclust
 * set-cover shape).
 *
 * <p><b>This is the load-bearing API.</b> {@link #cluster(List)} takes a set of
 * protein sequences as in-memory {@link ProteinSequence} objects and returns
 * {@link ProteinCluster} objects — each a representative plus its members with
 * per-member identity and coverage. No disk round-trip; callable from other
 * BBTools code (e.g. so a marker-gene factory or QuickBin can cluster proteins
 * in memory).</p>
 *
 * <p>Algorithm (basic version, correctness-first):</p>
 * <ul>
 * <li><b>Order:</b> sequences are processed longest-first (ties broken by id),
 * the CD-HIT convention that makes the longest sequence of a family its
 * representative.
 * <li><b>K-mer seed:</b> a growing index maps each k-mer to the representatives
 * that contain it. A new sequence is only aligned against representatives with
 * which it shares at least {@link #minSeedHits} distinct k-mers, avoiding
 * all-vs-all where the seeds allow.
 * <li><b>Extend/score:</b> each candidate representative is aligned to the
 * sequence with BLOSUM62 affine-gap local alignment ({@link AAAligner}).
 * <li><b>Assign:</b> the sequence joins the best candidate (highest identity)
 * that meets both {@link #minIdentity} and {@link #minCoverage}; if none
 * qualifies it becomes a new representative and its k-mers enter the index.
 * </ul>
 *
 * <p>Identity is the alignment percent identity (identical columns / aligned
 * columns). Coverage of a member is the fraction of its residues spanned by the
 * alignment; both the member's and the representative's coverage must meet the
 * threshold to guard against a high-identity hit over only a tiny region.</p>
 *
 * <p>Deferred (reported, not built): cross-run stable ids and incremental
 * update; a substitution-aware similar-k-mer prefilter (only exact or amino8
 * seeds here); orthology/paralog distinction; and any indexing / multithreading
 * for speed.</p>
 *
 * @author Eru
 */
public final class ProteinClusterer {

	/** K-mer length for seeding. */
	public int k=5;
	/** Minimum distinct shared k-mers with a representative to align against it. */
	public int minSeedHits=1;
	/** Minimum percent identity to the representative to join a cluster [0,100]. */
	public double minIdentity=90.0;
	/** Minimum aligned fraction of both member and representative [0,1]. */
	public double minCoverage=0.8;
	/** Use the amino8 reduced alphabet for seeding (more sensitive). */
	public boolean reducedSeed=false;

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
	 * Clusters the input sequences and returns the resulting families.
	 *
	 * @param seqs In-memory protein sequences to cluster.
	 * @return Clusters in representative-creation order (id ascending).
	 */
	public List<ProteinCluster> cluster(final List<ProteinSequence> seqs){
		if(seqs==null){throw new RuntimeException("Null sequence list.");}
		checkDuplicateIds(seqs);

		//Longest-first, ties broken by id, for a deterministic representative choice.
		final ArrayList<ProteinSequence> order=new ArrayList<ProteinSequence>(seqs);
		Collections.sort(order, LONGEST_FIRST);

		final ArrayList<ProteinCluster> clusters=new ArrayList<ProteinCluster>();
		//Growing seed index: k-mer -> distinct cluster indices whose rep holds it.
		final HashMap<Long, ArrayList<Integer>> index=new HashMap<Long, ArrayList<Integer>>();
		//Reusable per-sequence candidate tally over existing representatives.
		final ArrayList<Integer> hitCounts=new ArrayList<Integer>();
		final ArrayList<Integer> touched=new ArrayList<Integer>();

		for(final ProteinSequence s : order){
			final HashSet<Long> kmers=kmerSet(s.enc);

			//Tally distinct shared k-mers per existing representative.
			for(final Integer idx : touched){hitCounts.set(idx.intValue(), Integer.valueOf(0));}
			touched.clear();
			if(!kmers.isEmpty()){
				for(final Long km : kmers){
					final ArrayList<Integer> list=index.get(km);
					if(list!=null){
						for(final Integer idx : list){
							final int i=idx.intValue();
							final int c=hitCounts.get(i).intValue();
							if(c==0){touched.add(idx);}
							hitCounts.set(i, Integer.valueOf(c+1));
						}
					}
				}
			}

			//Find the best qualifying representative among the candidates.
			ProteinCluster best=null;
			double bestIdentity=-1, bestCoverage=0;
			final int nClusters=clusters.size();
			for(int i=0; i<nClusters; i++){
				//No seeds anywhere means k-mer seeding cannot help; fall back to all reps.
				final boolean candidate=kmers.isEmpty() || hitCounts.get(i).intValue()>=minSeedHits;
				if(!candidate){continue;}
				final ProteinCluster c=clusters.get(i);
				//Align member (query) to representative (target).
				final AAAlignment aln=AAAligner.align(s.enc, c.representative.enc);
				if(aln==null){continue;}
				final double identity=aln.pident();
				if(identity<minIdentity){continue;}
				final double memberCov=(aln.qStop-aln.qStart+1)/(double)s.length();
				final double repCov=(aln.tStop-aln.tStart+1)/(double)c.representative.length();
				if(memberCov<minCoverage || repCov<minCoverage){continue;}
				if(identity>bestIdentity){best=c; bestIdentity=identity; bestCoverage=memberCov;}
			}

			if(best!=null){
				best.add(s, bestIdentity, bestCoverage);
			}else{
				//New representative: create a cluster and register its k-mers.
				final int id=clusters.size();
				clusters.add(new ProteinCluster(id, s));
				hitCounts.add(Integer.valueOf(0));
				final Integer idxObj=Integer.valueOf(id);
				for(final Long km : kmers){
					ArrayList<Integer> list=index.get(km);
					if(list==null){list=new ArrayList<Integer>(); index.put(km, list);}
					list.add(idxObj);//distinct per rep, so no duplicate index per k-mer
				}
			}
		}

		return clusters;
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

	/** Fails loudly on duplicate identifiers within the input. */
	private static void checkDuplicateIds(final List<ProteinSequence> seqs){
		final HashSet<String> seen=new HashSet<String>();
		for(final ProteinSequence s : seqs){
			if(!seen.add(s.id)){
				throw new RuntimeException("Duplicate sequence identifier: '"+s.id+"'.");
			}
		}
	}

	/** Orders sequences longest-first, ties broken by id ascending (deterministic). */
	private static final Comparator<ProteinSequence> LONGEST_FIRST=new Comparator<ProteinSequence>(){
		@Override
		public int compare(final ProteinSequence a, final ProteinSequence b){
			if(a.length()!=b.length()){return a.length()>b.length() ? -1 : 1;}
			return a.id.compareTo(b.id);
		}
	};
}
