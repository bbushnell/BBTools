package assemble;

import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import map.LongHashMap;
import ukmer.Kmer;

/**
 * Enumerates read-supported exact K-1 endpoint overlaps without selecting or merging them.
 * Tied/competing destinations remain in the census. Counts are borrowed from the
 * original assembly table; no reference sequence or second table is consulted.
 * The borrowed contig list, sequences and count evidence must stay stable during
 * collection. Reused count/key scratch makes an instance non-reentrant; callbacks
 * must consume or copy scratch counts before returning.
 * @author Fischl
 */
public final class SameKJoinCollector{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Receives completed, eligible unitigs; no mincontig filter is applied here. */
	public SameKJoinCollector(final ArrayList<Contig> contigs_, final TadpoleGraph.Counts counts_, final int k_){
		if(contigs_==null || counts_==null || k_<2){throw new IllegalArgumentException("Same-K census requires contigs, counts and K>=2.");}
		contigs=contigs_; counts=counts_; k=k_; overlap=k-1;
		key=new Kmer(k); neighbor=new Kmer(k);
		if(key.kbig!=k){throw new IllegalArgumentException("Collector K must match the actual table K.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Indexes incoming endpoint words with collision chains and verifies every base.
	 * Physical endpoint ordering removes reverse-complement duplicate pairs, not
	 * alternate biological connections. Same-contig pairs are retained and tagged
	 * by the caller; a same-end hairpin is excluded as an unsupported self pairing.
	 * Each invocation numbers emitted events from one. Physical endpoint IDs are
	 * twice the contig list index plus zero for its left end or one for its right end;
	 * they are not Contig.id values. No edge or replacement sequence is created.
	 * Both complete endpoint K-mers must have positive table counts. Undefined endpoint
	 * words are skipped; interior bases outside these words are not validated here.
	 * Sources are visited in physical-end order and matching destinations in index-chain
	 * order. Event IDs describe this traversal, not a sequence-sorted or biological ranking.
	 * @param sink Synchronous consumer; must not mutate borrowed inputs or recursively collect
	 * @return Number of candidates delivered during this invocation
	 */
	public long collect(final Sink sink){
		assert(sink!=null) : "A census must consume each candidate before scratch counts are reused.";
		final int ends=Math.toIntExact(2L*contigs.size());
		final int[] heads=new int[ends], next=new int[ends];
		Arrays.fill(heads, -1); Arrays.fill(next, -1);
		final LongHashMap index=new LongHashMap(Math.max(16, ends));
		for(int end=0; end<ends; end++){
			final Contig c=contigs.get(end/2);
			final boolean reverse=(end&1)==1;
			if(!eligible(c) || !defined(c, reverse, 0, k)){continue;}
			final long hash=hash(c, reverse, 0);
			int group=index.get(hash);
			if(group<0){group=end; index.put(hash, group);}
			next[end]=heads[group]; heads[group]=end;
		}
		long events=0;
		for(int sourceEnd=0; sourceEnd<ends; sourceEnd++){
			final Contig a=contigs.get(sourceEnd/2);
			final boolean ar=(sourceEnd&1)==0;
			if(!eligible(a) || !defined(a, ar, a.length()-k, k)){continue;}
			final int group=index.get(hash(a, ar, a.length()-overlap));
			for(int destEnd=group<0 ? -1 : heads[group]; destEnd>=0; destEnd=next[destEnd]){
				if(sourceEnd>=destEnd){continue;}
				final Contig b=contigs.get(destEnd/2);
				final boolean br=(destEnd&1)==1;
				if(!exact(a, ar, b, br)){continue;}
				neighborhood(a, ar, a.length()-k, sourceCounts);
				neighborhood(b, br, 0, destCounts);
				// Two actual endpoint Kmers form the claimed one-step DBG adjacency.
				// A K-1 string match alone, with missing endpoint evidence, is not enough.
				if(sourceCounts[0]<1 || destCounts[0]<1){continue;}
				sink.accept(++events, sourceEnd, destEnd, a, ar, b, br, sourceCounts, destCounts);
			}
		}
		return events;
	}

	/** Requires non-null bases, unused/non-associate state and length at least K; no coverage filter. */
	private boolean eligible(final Contig c){return c!=null && c.bases!=null && !c.used() && !c.associate() && c.length()>=k;}

	/** Rejects undefined endpoint words without changing their source contig. */
	private static boolean defined(final Contig c, final boolean reverse, final int start, final int length){
		assert(start>=0 && start+length<=c.length()) : "Endpoint validation stays within an eligible unitig.";
		for(int i=0; i<length; i++){if(AminoAcid.baseToNumber[baseAt(c, reverse, start+i)]<0){return false;}}
		return true;
	}

	/** The hash only proposes candidates; exact() handles all hash collisions. */
	private long hash(final Contig c, final boolean reverse, final int start){
		assert(start>=0 && start+overlap<=c.length()) : "K-1 anchor must fit within the endpoint.";
		long hash=0;
		for(int i=0; i<overlap; i++){hash=hash*0x9E3779B185EBCA87L+baseAt(c, reverse, start+i)+1;}
		return hash;
	}

	/** Verifies a suffix/prefix overlap in the actual proposed join orientations. */
	private boolean exact(final Contig a, final boolean ar, final Contig b, final boolean br){
		assert(a.length()>=k && b.length()>=k) : "Both arms must retain at least one non-overlap base.";
		for(int i=0; i<overlap; i++){
			if(baseAt(a, ar, a.length()-overlap+i)!=baseAt(b, br, i)){return false;}
		}
		return true;
	}

	/** Queries current/right4/left4 counts in the requested orientation; absent keys become zero.
	 * Mutates only this collector's key scratch and the supplied nine-count output. */
	private void neighborhood(final Contig c, final boolean reverse, final int start, final int[] out){
		SameKJoinFeatures.neighborhood(counts, c, reverse, start, k, key, neighbor, out);
	}

	/** Copies the oriented half-open interval [start,stop) into a new independent array
	 * for reference labeling, not a replacement for the full contig. Coordinates refer
	 * to the requested orientation; bases use the same normalization as endpoint matching. */
	public static byte[] context(final Contig c, final boolean reverse, final int start, final int stop){
		assert(start>=0 && stop>start && stop<=c.length()) : "Truth context must consist of retained contig bases.";
		final byte[] result=new byte[stop-start];
		for(int p=start; p<stop; p++){result[p-start]=baseAt(c, reverse, p);}
		return result;
	}

	/**
	 * Normalizes case and orientation through the native DNA lookup without changing source bases.
	 * U/u is treated as T; undefined and negative byte values become N.
	 * pos is zero-based in the requested orientation; the result is uppercase A/C/G/T or N.
	 */
	static byte baseAt(final Contig c, final boolean reverse, final int pos){
		assert(pos>=0 && pos<c.length()) : "Candidate geometry may not read beyond its contig.";
		final byte raw=c.bases[reverse ? c.length()-1-pos : pos];
		if(raw<0){return 'N';}
		final int code=AminoAcid.baseToNumber[raw];
		return code<0 ? (byte)'N' : AminoAcid.numberToBase[reverse ? 3-code : code];
	}

	/** Receives read-only contigs and count scratch valid only until accept returns. */
	public interface Sink{
		/** Receives one event, physical endpoint IDs, borrowed contigs/orientations and
		 * current/right4/left4 count arrays. Both count arrays are reused by the collector. */
		void accept(long event, int sourceEnd, int destEnd, Contig source, boolean sourceReverse,
				Contig dest, boolean destReverse, int[] sourceCounts, int[] destCounts);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed unitig list; indices determine stable physical endpoint identities. */
	private final ArrayList<Contig> contigs;
	/** Borrowed original evidence, queried without adding candidate sequence counts. */
	private final TadpoleGraph.Counts counts;
	/** Exact endpoint word length and required K-1 suffix/prefix overlap, in bases. */
	private final int k, overlap;
	/** Reused complete-word and neighbor-query scratch, owned by this collector. */
	private final Kmer key, neighbor;
	/** Nine raw nucleotide-indexed depths for each callback, overwritten by later candidates. */
	private final int[] sourceCounts=new int[9], destCounts=new int[9];
}
