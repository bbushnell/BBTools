package assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import dna.AminoAcid;
import map.LongHashSet;
import ukmer.Kmer;

/**
 * Samples internal accepted steps of formed unitigs, independently of truth.
 * Virtual arms are training controls, never new assembly contigs or merge inputs.
 * Memory is proportional to the requested sample, not total sequence length.
 * @author Fischl
 */
final class SameKAcceptedControls {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrows immutable raw unitigs and their original count table. */
	SameKAcceptedControls(final ArrayList<Contig> contigs_, final TadpoleGraph.Counts counts_, final int k_){
		if(contigs_==null || counts_==null || k_<2){throw new IllegalArgumentException("Controls require unitigs, counts and K>=2.");}
		contigs=contigs_; counts=counts_; k=k_;
		key=new Kmer(k); neighbor=new Kmer(k);
		if(key.kbig!=k){throw new IllegalArgumentException("Control K must equal the actual table K.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Uniform sampling without replacement over eligible cuts, not over contigs.
	 * Report sampled/population before truth filtering; label failures do not redraw.
	 */
	Result collect(final int requested, final long seed, final Sink sink){
		if(requested<0 || sink==null){throw new IllegalArgumentException("A nonnegative sample count and sink are required.");}
		long population=0;
		for(Contig c : contigs){population=Math.addExact(population, available(c));}
		final int sampled=(int)Math.min(population, requested);
		final long[] selected=sample(population, sampled, seed);
		long offset=0;
		int next=0;
		for(int id=0; id<contigs.size() && next<sampled; id++){
			final Contig c=contigs.get(id);
			final long available=available(c), limit=offset+available;
			final int first=next;
			while(next<sampled && selected[next]<limit){next++;}
			if(next>first){emit(c, id, selected, first, next, offset, sink);}
			offset=limit;
		}
		assert(next==sampled) : "Every sampled global cut must map to exactly one eligible contig.";
		return new Result(population, sampled);
	}

	/** Both virtual arms must have at least K+2 bases; keep the native near-K guard. */
	private long available(final Contig c){
		return c==null || c.bases==null || c.used() || c.associate() ? 0 : Math.max(0L, (long)c.length()-k-4);
	}

	/** One rolling count traversal supplies all selected arm means in this contig. */
	private void emit(final Contig c, final int id, final long[] selected, final int first,
			final int end, final long offset, final Sink sink){
		assert(first<end && available(c)>0) : "Only contigs with sampled eligible internal cuts need a count scan.";
		final long[] prefix=new long[end-first];
		long total=0;
		int next=first;
		key.clear();
		for(int pos=0; pos<c.length(); pos++){
			final byte base=SameKJoinCollector.baseAt(c, false, pos);
			if(AminoAcid.baseToNumber[base]<0){throw new IllegalArgumentException("Undefined base in a formed control unitig: "+id);}
			key.addRight(base);
			if(pos<k-1){continue;}
			final int count=counts.count(key);
			if(count<1){throw new IllegalStateException("Missing internal unitig Kmer at unitig="+id+" end="+pos);}
			total+=count;
			if(next<end && pos+1==cut(selected[next], offset)){
				prefix[next-first]=total; next++;
			}
		}
		assert(next==end) : "Sorted cuts must all be reached by the single rolling count traversal.";
		for(int i=first; i<end; i++){
			final int cut=cut(selected[i], offset), sourcePositions=cut-k+1, destPositions=c.length()-cut;
			assert(sourcePositions>=3 && destPositions>=3) : "K+2 virtual arms each retain three Kmer positions.";
			final double sourceMean=prefix[i-first]/(double)sourcePositions;
			final double destMean=(total-prefix[i-first])/(double)destPositions;
			SameKJoinFeatures.neighborhood(counts, c, false, cut-k, k, key, neighbor, sourceCounts);
			SameKJoinFeatures.neighborhood(counts, c, false, cut-k+1, k, key, neighbor, destCounts);
			final int sourceInward=AminoAcid.baseToNumber[SameKJoinCollector.baseAt(c, false, cut-k-1)];
			final int destInward=AminoAcid.baseToNumber[SameKJoinCollector.baseAt(c, false, cut+1)];
			sink.accept(id, cut, c, sourceMean, destMean, sourceCounts, destCounts, sourceInward, destInward);
		}
	}

	/** Source is [0,cut); destination is [cut-K+1,L), sharing exactly K-1 bases. */
	private int cut(final long index, final long offset){return Math.toIntExact(index-offset+k+2);}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Floyd's uniform subset algorithm, with primitive storage and sorted output. */
	static long[] sample(final long population, final int count, final long seed){
		if(population<0 || count<0 || count>population){throw new IllegalArgumentException("Sample must fit its population.");}
		if(count==0){return new long[0];}
		final LongHashSet selected=new LongHashSet(Math.max(16, count));
		final Random random=new Random(seed);
		for(long j=population-count; j<population; j++){
			final long candidate=bounded(random, j+1);
			if(!selected.add(candidate)){
				final boolean added=selected.add(j);
				assert(added) : "Floyd sampling's new upper endpoint cannot already be selected.";
			}
		}
		final long[] result=selected.toArray();
		assert(result.length==count) : "Sampling without replacement must preserve the requested cardinality.";
		Arrays.sort(result);
		return result;
	}

	/** Java-8-compatible unbiased bounded-long draw using rejection of the incomplete tail. */
	private static long bounded(final Random random, final long bound){
		assert(bound>0) : "Floyd sampling draws in [0,j] for a nonnegative population index.";
		long bits, value;
		do{bits=random.nextLong()>>>1; value=bits%bound;}while(bits-value+(bound-1)<0);
		return value;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Scratch neighborhoods remain valid only during the callback; no labels are implied. */
	interface Sink {
		void accept(int contigId, int cut, Contig contig, double sourceMean, double destMean,
				int[] sourceCounts, int[] destCounts, int sourceInwardBase, int destInwardBase);
	}

	/** Pre-label sampling metadata for selection weighting and reproducibility. */
	static final class Result {
		Result(final long population_, final int sampled_){population=population_; sampled=sampled_;}
		double probability(){return population==0 ? 0 : sampled/(double)population;}
		final long population;
		final int sampled;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final ArrayList<Contig> contigs;
	private final TadpoleGraph.Counts counts;
	private final int k;
	private final Kmer key, neighbor;
	private final int[] sourceCounts=new int[9], destCounts=new int[9];
}
