package assemble;

import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import map.LongHashMap;
import ukmer.Kmer;

/**
 * Enumerates all exact K-1 endpoint overlaps without selecting or merging them.
 * Tied/competing destinations remain in the census. Counts are borrowed from the
 * original assembly table; no reference sequence or second table is consulted.
 * @author Fischl
 */
final class SameKJoinCollector {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Receives completed, eligible unitigs; no mincontig filter is applied here. */
	SameKJoinCollector(final ArrayList<Contig> contigs_, final TadpoleGraph.Counts counts_, final int k_){
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
	 */
	long collect(final Sink sink){
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

	/** Native formation guards already ran; retain every non-retired complete endpoint. */
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

	/** Count order is current, right A/C/G/T, left A/C/G/T; identities are not sorted away. */
	private void neighborhood(final Contig c, final boolean reverse, final int start, final int[] out){
		SameKJoinFeatures.neighborhood(counts, c, reverse, start, k, key, neighbor, out);
	}

	/** Bounded persistent context for reference labeling, not a replacement for the full contig. */
	static byte[] context(final Contig c, final boolean reverse, final int start, final int stop){
		assert(start>=0 && stop>start && stop<=c.length()) : "Truth context must consist of retained contig bases.";
		final byte[] result=new byte[stop-start];
		for(int p=start; p<stop; p++){result[p-start]=baseAt(c, reverse, p);}
		return result;
	}

	/** Normalizes ASCII case and orientation without mutating stored sequences. */
	static byte baseAt(final Contig c, final boolean reverse, final int pos){
		assert(pos>=0 && pos<c.length()) : "Candidate geometry may not read beyond its contig.";
		final byte raw=c.bases[reverse ? c.length()-1-pos : pos];
		if(raw<0){return 'N';}
		final int code=AminoAcid.baseToNumber[raw];
		return code<0 ? (byte)'N' : AminoAcid.numberToBase[reverse ? 3-code : code];
	}

	/** Receives read-only contigs and count scratch valid only until accept returns. */
	interface Sink {
		void accept(long event, int sourceEnd, int destEnd, Contig source, boolean sourceReverse,
				Contig dest, boolean destReverse, int[] sourceCounts, int[] destCounts);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final ArrayList<Contig> contigs;
	private final TadpoleGraph.Counts counts;
	private final int k, overlap;
	private final Kmer key, neighbor;
	private final int[] sourceCounts=new int[9], destCounts=new int[9];
}
