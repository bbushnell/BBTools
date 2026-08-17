package prot;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;

import dna.AminoAcid;
import map.LongArrayListHashMap;
import map.LongHashSet;
import structures.IntList;
import structures.LongList;

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
 * <p>The seed index (k-mer -&gt; representative indices), per-sequence k-mer
 * membership test, and per-representative hit tally reuse instance-level
 * primitive scratch across calls to {@link #cluster(List)} instead of boxing
 * every k-mer/count in {@code java.util} collections -- zero heap allocation
 * per sequence beyond growth-triggered scratch resizing. This is a storage-
 * strategy change only: candidate order, seeding, and assignment logic are
 * unmodified, so output is identical to the boxed-collection implementation.</p>
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

	/** Reused per-sequence scratch: membership test for {@link #kmerSet}'s dedup pass. */
	private final LongHashSet kmerSeen=new LongHashSet(64);
	/** Reused per-sequence scratch: this sequence's distinct k-mers, first-seen order. */
	private final LongList kmerScratch=new LongList(64);

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
		final LongArrayListHashMap<Integer> index=new LongArrayListHashMap<Integer>();
		//Reusable per-sequence candidate tally over existing representatives.
		final IntList hitCounts=new IntList();
		final IntList touched=new IntList();

		for(final ProteinSequence s : order){
			kmerSeen.clear();
			kmerScratch.clear();
			kmerSet(s.enc, kmerSeen, kmerScratch);

			//Tally distinct shared k-mers per existing representative.
			for(int ti=0; ti<touched.size; ti++){hitCounts.set(touched.get(ti), 0);}
			touched.clear();
			if(kmerScratch.size>0){
				for(int ki=0; ki<kmerScratch.size; ki++){
					final long km=kmerScratch.get(ki);
					final ArrayList<Integer> list=index.get(km);
					if(list!=null){
						for(final Integer idxObj : list){
							final int i=idxObj.intValue();
							final int c=hitCounts.get(i);
							if(c==0){touched.add(i);}
							hitCounts.set(i, c+1);
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
				final boolean candidate=(kmerScratch.size==0) || hitCounts.get(i)>=minSeedHits;
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
				hitCounts.add(0);
				for(int ki=0; ki<kmerScratch.size; ki++){
					index.put(kmerScratch.get(ki), id);
				}
			}
		}

		return clusters;
	}

	/**
	 * Populates {@code out} with the distinct k-mers of an encoded sequence, in
	 * first-seen order (iteration order does not affect clustering output --
	 * only set membership does). K-mers containing X (ambiguous) are skipped.
	 * Uses the amino8 reduced alphabet when {@link #reducedSeed} is set.
	 * @param enc Encoded residues.
	 * @param seen Reused dedup set; caller must clear before this call.
	 * @param out Reused output list; caller must clear before this call.
	 */
	private void kmerSet(final byte[] enc, final LongHashSet seen, final LongList out){
		final int bits=reducedSeed ? 3 : 5;
		assert(bits*k<=62) : "k too large for packed k-mer: k="+k+", bits="+bits;
		if(enc.length<k){return;}
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
			if(valid>=k && seen.add(kmer)){out.add(kmer);}
		}
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
