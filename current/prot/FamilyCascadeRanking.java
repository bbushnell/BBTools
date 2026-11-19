package prot;

import java.util.Arrays;
import java.util.HashMap;

import ddl.DDLFamilySearch;

/** Shared DDL-plus-cosine ranking for calibrated family cascades.
 * @author Yelan
 */
public final class FamilyCascadeRanking{

	private FamilyCascadeRanking(){}

	/*--------------------------------------------------------------*/
	/*----------------             Methods          ----------------*/
	/*--------------------------------------------------------------*/

	/** Ranks families by raw DDL matches + dimer cosine descending, then family rank ascending. */
	public static int rank(final DDLFamilySearch.Workspace work,
			final ProteinSpectrum.Spectrum querySpectrum, final double queryNorm,
			final ProteinSpectrum.Spectrum[] spectra, final double[] norms,
			final int limit, final int[] out, final Scratch scratch){
		return rank(work, querySpectrum, queryNorm, spectra, norms, null, limit, out, scratch);
	}

	/** Ranks by raw DDL + cosine; optional canonical ranks break ties independent of profile storage order. */
	public static int rank(final DDLFamilySearch.Workspace work,
			final ProteinSpectrum.Spectrum querySpectrum, final double queryNorm,
			final ProteinSpectrum.Spectrum[] spectra, final double[] norms, final int[] canonicalRanks,
			final int limit, final int[] out, final Scratch scratch){
		if(work==null || querySpectrum==null || spectra==null || norms==null || out==null || scratch==null ||
				spectra.length!=norms.length || work.scores.length!=spectra.length || work.touched.length!=spectra.length ||
				(canonicalRanks!=null && canonicalRanks.length!=spectra.length) || limit<1 || out.length<limit){
			throw new IllegalArgumentException("Invalid DDL/cosine rank dimensions");
		}
		if(scratch.capacity()<limit){throw new IllegalArgumentException("Scratch capacity "+scratch.capacity()+" < limit "+limit);}
		scratch.clear();
		for(int family=0; family<spectra.length; family++){
			final int rank=canonicalRanks==null ? family : canonicalRanks[family];
			if(rank<0){throw new IllegalArgumentException("Negative canonical rank for family index "+family);}
			final ProteinSpectrum.Spectrum spectrum=spectra[family];
			final double norm=norms[family];
			// Retain cosine's fail-loud preconditions even for candidates whose bound loses.
			if(spectrum==null){throw new IllegalArgumentException("Missing calibrated spectrum for family "+family);}
			if(querySpectrum.k!=spectrum.k){throw new IllegalArgumentException("Cosine needs the same alphabet and k");}
			if(!Double.isFinite(queryNorm) || !Double.isFinite(norm) || queryNorm<0 || norm<0){
				throw new IllegalArgumentException("Invalid cosine norms");
			}
			scratch.familiesConsidered++;
			final int ddl=work.scores[family];
			if(queryNorm==0 || norm==0){
				scratch.zeroNormShortcuts++;
				scratch.offer(family, rank, ddl+0.0, limit);
			}else if(scratch.size==limit && !better((double)ddl+1.0, rank, scratch.scores[0], scratch.ranks[0])){
				// cosine<=1: this offer cannot change the heap, including canonical-rank ties.
				scratch.boundSkips++;
			}else{
				scratch.cosineCalls++;
				scratch.offer(family, rank, score(querySpectrum, queryNorm, ddl, spectrum, norm), limit);
			}
		}
		return scratch.drain(out);
	}

	/** Converts one storage-order family ID vector to canonical rank tie-breaks. */
	public static int[] canonicalRanks(final String[] storageIds, final String[] canonicalIds){
		if(storageIds==null || canonicalIds==null || storageIds.length!=canonicalIds.length || storageIds.length<1){
			throw new IllegalArgumentException("Family ID vectors must have identical nonzero lengths");
		}
		final HashMap<String,Integer> ranks=new HashMap<String,Integer>(canonicalIds.length*2);
		for(int i=0; i<canonicalIds.length; i++){
			if(canonicalIds[i]==null || ranks.put(canonicalIds[i], i)!=null){
				throw new IllegalArgumentException("Duplicate or null canonical family ID: "+canonicalIds[i]);
			}
		}
		final int[] out=new int[storageIds.length];
		for(int i=0; i<storageIds.length; i++){
			final Integer rank=ranks.remove(storageIds[i]);
			if(rank==null){throw new IllegalArgumentException("Storage family missing from canonical roster: "+storageIds[i]);}
			out[i]=rank.intValue();
		}
		if(!ranks.isEmpty()){throw new IllegalArgumentException("Canonical families absent from storage roster: "+ranks.keySet());}
		return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Fixed-capacity worst-root heap for the best candidate families. */
	public static final class Scratch{
		public Scratch(final int limit){
			if(limit<1){throw new IllegalArgumentException("Invalid shortlist limit: "+limit);}
			families=new int[limit]; ranks=new int[limit]; scores=new double[limit];
		}
		public int capacity(){return families.length;}
		public void clear(){size=0; familiesConsidered=cosineCalls=zeroNormShortcuts=boundSkips=0;}
		void offer(final int family, final int rank, final double score, final int limit){
			if(!Double.isFinite(score)){throw new IllegalArgumentException("Invalid cascade score for family "+family+": "+score);}
			if(size<limit){
				families[size]=family; ranks[size]=rank; scores[size]=score; siftUp(size++);
			}else if(better(score, rank, scores[0], ranks[0])){
				families[0]=family; ranks[0]=rank; scores[0]=score; siftDown(0);
			}
		}
		int drain(final int[] out){
			final int n=size;
			for(int i=n-1; i>=0; i--){
				out[i]=families[0];
				size--;
				if(size>0){
					families[0]=families[size]; scores[0]=scores[size];
					ranks[0]=ranks[size];
					siftDown(0);
				}
			}
			Arrays.fill(out, n, out.length, -1);
			return n;
		}
		private void siftUp(int child){
			while(child>0){
				final int parent=(child-1)>>>1;
				if(!worse(child, parent)){break;}
				swap(child, parent); child=parent;
			}
		}
		private void siftDown(int parent){
			while(true){
				final int left=parent*2+1;
				if(left>=size){break;}
				final int right=left+1;
				int worst=left;
				if(right<size && worse(right, left)){worst=right;}
				if(!worse(worst, parent)){break;}
				swap(parent, worst); parent=worst;
			}
		}
		private boolean worse(final int a, final int b){
			return scores[a]!=scores[b] ? scores[a]<scores[b] : ranks[a]>ranks[b];
		}
		private void swap(final int a, final int b){
			final int tf=families[a]; final double ts=scores[a];
			final int tr=ranks[a];
			families[a]=families[b]; scores[a]=scores[b];
			ranks[a]=ranks[b];
			families[b]=tf; scores[b]=ts;
			ranks[b]=tr;
		}
		/** Per-rank counters; on success considered equals calls + zero shortcuts + bound skips. */
		public long familiesConsidered(){return familiesConsidered;}
		/** Actual dot-product evaluations in this rank call, excluding later acceptance gates. */
		public long cosineCalls(){return cosineCalls;}
		public long zeroNormShortcuts(){return zeroNormShortcuts;}
		public long boundSkips(){return boundSkips;}
		private final int[] families, ranks;
		private final double[] scores;
		private int size;
		private long familiesConsidered, cosineCalls, zeroNormShortcuts, boundSkips;
	}

	/** Fixed-capacity worst-root heap for one family's hardest negative query ordinals. */
	public static final class QueryScratch{
		public QueryScratch(final int limit){
			if(limit<1){throw new IllegalArgumentException("Invalid query limit: "+limit);}
			ordinals=new int[limit]; scores=new double[limit];
		}
		public int capacity(){return ordinals.length;}
		public void clear(){size=0;}
		public void offer(final int ordinal, final double score, final int limit){
			if(limit<1 || limit>capacity()){throw new IllegalArgumentException("Invalid query limit "+limit+" for capacity "+capacity());}
			if(ordinal<0 || !Double.isFinite(score)){throw new IllegalArgumentException("Invalid negative query "+ordinal+": "+score);}
			if(size<limit){
				ordinals[size]=ordinal; scores[size]=score; siftUp(size++);
			}else if(better(score, ordinal, scores[0], ordinals[0])){
				ordinals[0]=ordinal; scores[0]=score; siftDown(0);
			}
		}
		public int drain(final int[] out){
			if(out.length<size){throw new IllegalArgumentException("Output length "+out.length+" < ranked queries "+size);}
			final int n=size;
			for(int i=n-1; i>=0; i--){
				out[i]=ordinals[0];
				size--;
				if(size>0){
					ordinals[0]=ordinals[size]; scores[0]=scores[size];
					siftDown(0);
				}
			}
			Arrays.fill(out, n, out.length, -1);
			return n;
		}
		private void siftUp(int child){
			while(child>0){
				final int parent=(child-1)>>>1;
				if(!worse(child, parent)){break;}
				swap(child, parent); child=parent;
			}
		}
		private void siftDown(int parent){
			while(true){
				final int left=parent*2+1;
				if(left>=size){break;}
				final int right=left+1;
				int worst=left;
				if(right<size && worse(right, left)){worst=right;}
				if(!worse(worst, parent)){break;}
				swap(parent, worst); parent=worst;
			}
		}
		private boolean worse(final int a, final int b){
			return scores[a]!=scores[b] ? scores[a]<scores[b] : ordinals[a]>ordinals[b];
		}
		private void swap(final int a, final int b){
			final int to=ordinals[a]; final double ts=scores[a];
			ordinals[a]=ordinals[b]; scores[a]=scores[b];
			ordinals[b]=to; scores[b]=ts;
		}
		private final int[] ordinals;
		private final double[] scores;
		private int size;
	}

	/** Flat fixed-capacity worst-root heaps for each family's hardest negative query ordinals. */
	public static final class QueryTopK{
		public QueryTopK(final int families_, final int limit_){
			if(families_<1 || limit_<1 || (long)families_*limit_>Integer.MAX_VALUE){
				throw new IllegalArgumentException("Invalid per-family query top-k dimensions: "+families_+" x "+limit_);
			}
			families=families_; limit=limit_;
			final int cells=families*limit;
			sizes=new int[families]; ordinals=new int[cells]; scores=new double[cells];
			Arrays.fill(ordinals, -1);
		}
		public int families(){return families;}
		public int limit(){return limit;}
		public int size(final int family){
			checkFamily(family);
			return sizes[family];
		}
		public void clear(){Arrays.fill(sizes, 0);}
		public void offer(final int family, final int ordinal, final double score){
			checkFamily(family);
			if(ordinal<0 || !Double.isFinite(score)){
				throw new IllegalArgumentException("Invalid negative query "+ordinal+" for family "+family+": "+score);
			}
			final int offset=family*limit;
			int size=sizes[family];
			if(size<limit){
				final int cell=offset+size;
				ordinals[cell]=ordinal; scores[cell]=score;
				sizes[family]=++size; siftUp(offset, size-1);
			}else if(better(score, ordinal, scores[offset], ordinals[offset])){
				ordinals[offset]=ordinal; scores[offset]=score;
				siftDown(offset, size, 0);
			}
		}
		public int drain(final int family, final int[] out){
			checkFamily(family);
			final int offset=family*limit, n=sizes[family];
			if(out.length<n){throw new IllegalArgumentException("Output length "+out.length+" < ranked queries "+n);}
			int size=n;
			for(int i=n-1; i>=0; i--){
				out[i]=ordinals[offset];
				size--;
				if(size>0){
					final int last=offset+size;
					ordinals[offset]=ordinals[last]; scores[offset]=scores[last];
					siftDown(offset, size, 0);
				}
			}
			sizes[family]=0;
			Arrays.fill(out, n, out.length, -1);
			return n;
		}
		private void siftUp(final int offset, int child){
			while(child>0){
				final int parent=(child-1)>>>1;
				if(!worse(offset+child, offset+parent)){break;}
				swap(offset+child, offset+parent); child=parent;
			}
		}
		private void siftDown(final int offset, final int size, int parent){
			while(true){
				final int left=parent*2+1;
				if(left>=size){break;}
				final int right=left+1;
				int worst=left;
				if(right<size && worse(offset+right, offset+left)){worst=right;}
				if(!worse(offset+worst, offset+parent)){break;}
				swap(offset+parent, offset+worst); parent=worst;
			}
		}
		private boolean worse(final int a, final int b){
			return scores[a]!=scores[b] ? scores[a]<scores[b] : ordinals[a]>ordinals[b];
		}
		private void swap(final int a, final int b){
			final int to=ordinals[a]; final double ts=scores[a];
			ordinals[a]=ordinals[b]; scores[a]=scores[b];
			ordinals[b]=to; scores[b]=ts;
		}
		private void checkFamily(final int family){
			if(family<0 || family>=families){throw new IllegalArgumentException("Invalid family index: "+family);}
		}
		private final int families, limit;
		private final int[] sizes, ordinals;
		private final double[] scores;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Exact scalar for both query-family ranking and per-family hard-negative ranking. */
	public static double score(final ProteinSpectrum.Spectrum query, final double queryNorm, final int ddl,
			final ProteinSpectrum.Spectrum familySpectrum, final double familyNorm){
		return ddl+ProteinSpectrum.cosine(query, familySpectrum, queryNorm, familyNorm);
	}

	private static boolean better(final double scoreA, final int familyA, final double scoreB, final int familyB){
		return scoreA!=scoreB ? scoreA>scoreB : familyA<familyB;
	}
}
