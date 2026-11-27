package prok;

import java.util.Arrays;

/** Worker-owned native-first proposal interleaving with duplicate removal.
 * Native order contains positive shared-key counts only; sketch order contains
 * every original model. Each round proposes native[i], then sketch[i]. A repeated
 * model consumes that proposal slot but never another attempted-model position.
 * This is the stable zip of the two orders followed by duplicate removal, not
 * alternation of successful emissions. All arrays are reused after construction.
 * @author Brian, Raiden
 */
public final class NcrnaDdlInterleave{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaDdlInterleave(final int models){
		require(models>0, "Interleaving requires the complete positive model count");
		nativeKeys=new long[models];
		sketchKeys=new long[models];
		nativeScores=new float[models];
		sketchScores=new float[models];
		ranks=new int[models];
		order=new int[models];
		seen=new boolean[models];
		sketchFirst=new boolean[models];
		filtered=new boolean[models];
	}

	/*--------------------------------------------------------------*/
	/*----------------          Query Order         ----------------*/
	/*--------------------------------------------------------------*/

	/** Replaces both score orders and the hybrid permutation. Counts are unsigned
	 * 16-bit, matching TrnaKmerIndex.MAX_COUNT=65535; conversion to float is exact.
	 * Both inputs are borrowed read-only and neither array is retained. */
	public void build(final char[] counts, final float[] scores){
		ready=false;
		size=0;
		require(counts!=null && scores!=null && counts.length==ranks.length && scores.length==ranks.length,
			"Native counts and sketch scores must retain the same complete original model order");
		Arrays.fill(seen, false);
		Arrays.fill(sketchFirst, false);
		int positives=0;
		for(int model=0; model<ranks.length; model++){
			final float score=scores[model];
			require(Float.isFinite(score) && score>=0, "Negative or nonfinite sketch scores cannot define a scientific order");
			nativeScores[model]=counts[model];
			sketchScores[model]=score;
			if(counts[model]>0){nativeKeys[positives++]=pack(counts[model], model);}
			// Positive finite float bits sort numerically; normalize negative zero
			// because NcrnaDdlRanks.best treats both zero representations as tied.
			sketchKeys[model]=pack(Float.floatToIntBits(score==0 ? 0 : score), model);
		}
		sort(nativeKeys, positives);
		sort(sketchKeys, sketchKeys.length);
		int next=1;
		for(int position=0; position<ranks.length; position++){
			if(position<positives){next=emit(nativeKeys[positives-1-position], next, false);}
			next=emit(sketchKeys[ranks.length-1-position], next, true);
		}
		assert(next==ranks.length+1) : "The complete sketch order must place every model exactly once after duplicate removal";
		size=ranks.length;
		ready=true;
	}

	/** Experimental caller policy, separate from the complete saved diagnostic.
	 * Native proposals pass the absolute floor, with the same positive top1 fallback
	 * as TrnaKmerIndex.shortlist. Sketch proposals require score>0, regardless of
	 * native count. Filter before indexed native-first zip, then deduplicate and cap.
	 * Ties use original model order. Returns the actual emitted length, possibly0.
	 * Neither cap nor a duplicate causes refill from later proposals after the cap.
	 */
	public int buildCaller(final char[] counts, final float[] scores, final int floor, final int cap){
		ready=false;
		size=0;
		require(counts!=null && scores!=null && counts.length==ranks.length && scores.length==ranks.length
			&& floor>=0 && floor<=65535 && cap>0 && cap<=ranks.length,
			"Caller ordering needs the complete model counts/scores, an unsigned16bit floor and a positive model-position cap");
		Arrays.fill(seen, false);
		Arrays.fill(sketchFirst, false);
		Arrays.fill(ranks, 0);
		int natives=0, sketches=0, best=-1;
		for(int model=0; model<ranks.length; model++){
			final int count=counts[model];
			final float score=scores[model];
			require(Float.isFinite(score) && score>=0, "Caller sketch proposals require finite nonnegative scores");
			nativeScores[model]=count;
			sketchScores[model]=score;
			if(count>0){
				if(best<0 || count>counts[best]){best=model;}
				if(count>=floor){nativeKeys[natives++]=pack(count, model);}
			}
			if(score>0){sketchKeys[sketches++]=pack(Float.floatToIntBits(score), model);}
		}
		if(natives==0 && best>=0){nativeKeys[natives++]=pack(counts[best], best);}
		sort(nativeKeys, natives);
		sort(sketchKeys, sketches);
		int next=1;
		for(int position=0, rounds=Math.max(natives, sketches); position<rounds && next<=cap; position++){
			if(position<natives){next=emit(nativeKeys[natives-1-position], next, false);}
			if(position<sketches && next<=cap){next=emit(sketchKeys[sketches-1-position], next, true);}
		}
		size=next-1;
		assert(size<=cap && size<=natives+sketches)
			: "The caller emits only distinct eligible proposals within its attempted-model cap";
		ready=true;
		return size;
	}

	/** Orders only the supplied native-eligible prefix, then applies the same
	 * caller cap. Zero DDL scores remain eligible; ties use original model ID.
	 * Native cutoff, keep-best fallback and counts belong to TrnaKmerIndex. */
	public int buildEligibleCaller(final float[] scores, final int[] eligible, final int count, final int cap){
		return buildEligibleCaller(scores, eligible, count, cap, false);
	}
	/** A reviewed secondary sorter may supply the complete native-eligible order.
	 * Raw scores and caller cap still retain their original meanings. */
	public int buildEligibleCaller(final float[] scores, final int[] eligible, final int count, final int cap,
			final boolean presorted){
		require(scores!=null && scores.length==ranks.length && eligible!=null && count>=0
			&& count<=eligible.length && count<=ranks.length && cap>0 && cap<=ranks.length,
			"Ordering-only DDL ranking needs the complete score vector and actual native-eligible prefix");
		ready=false;
		Arrays.fill(seen, false);
		Arrays.fill(ranks, 0);
		Arrays.fill(sketchFirst, false);
		for(int position=0; position<count; position++){
			final int model=eligible[position];
			require(model>=0 && model<ranks.length && !seen[model],
				"Native eligibility must contain distinct original model IDs");
			seen[model]=true;
			final float score=scores[model];
			require(Float.isFinite(score) && score>=0, "Equal-bin match counts must be finite and nonnegative");
			sketchScores[model]=score;
			sketchKeys[position]=pack(Float.floatToIntBits(score==0 ? 0 : score), model);
		}
		Arrays.fill(seen, false);
		if(!presorted){sort(sketchKeys, count);}
		size=Math.min(cap, count);
		for(int position=0; position<size; position++){emit(sketchKeys[presorted ? position : count-1-position], position+1, true);}
		ready=true;
		return size;
	}

	/** Complete sketch-only order, including zero-score ties by original ordinal.
	 * No native counts, floor, positive-only filter or interleave are consulted. */
	public int buildSketchCaller(final float[] scores, final int cap){
		ready=false;
		size=0;
		require(scores!=null && scores.length==ranks.length && cap>0 && cap<=ranks.length,
			"Sketch-only ordering needs all original model scores and a positive position cap");
		Arrays.fill(seen, false);
		Arrays.fill(sketchFirst, false);
		Arrays.fill(ranks, 0);
		Arrays.fill(nativeScores, Float.NaN);
		for(int model=0; model<scores.length; model++){
			final float score=scores[model];
			require(Float.isFinite(score) && score>=0, "Sketch-only order requires finite nonnegative scores, including zero");
			sketchScores[model]=score;
			sketchKeys[model]=pack(Float.floatToIntBits(score==0 ? 0 : score), model);
		}
		sort(sketchKeys, sketchKeys.length);
		for(int position=0; position<cap; position++){emit(sketchKeys[sketchKeys.length-1-position], position+1, true);}
		size=cap;
		ready=true;
		return size;
	}

	/** Explicit GC-adjusted order: finite negative scores remain eligible. Hard
	 * exclusions never consume an attempted position. The old nonnegative routes
	 * keep their original validation/packing. Original ordinal breaks every tie. */
	public int buildFilteredSketchCaller(final float[] scores, final boolean[] eligible, final int cap){
		ready=false;
		size=0;
		require(scores!=null && eligible!=null && scores.length==ranks.length && eligible.length==ranks.length
			&& cap>0 && cap<=ranks.length, "Filtered sketch ordering needs complete scores/eligibility and a positive cap");
		Arrays.fill(seen, false);
		Arrays.fill(sketchFirst, false);
		Arrays.fill(ranks, 0);
		Arrays.fill(nativeScores, Float.NaN);
		int count=0;
		for(int model=0; model<scores.length; model++){
			final float score=scores[model];
			require(Float.isFinite(score), "GC-adjusted proposal scores must stay finite; hard exclusions use their mask");
			sketchScores[model]=score;
			if(!eligible[model]){continue;}
			// IEEE negative bits run backwards. This transform preserves signed
			// numerical order, including subnormals, without allocating a sorter.
			final int bits=Float.floatToIntBits(score==0 ? 0 : score);
			final int sortable=bits^((bits>>31)&0x7fffffff);
			sketchKeys[count++]=((long)sortable<<32)|(0xffffffffL-model);
		}
		sort(sketchKeys, count);
		final int limit=Math.min(cap, count);
		for(int position=0; position<limit; position++){emit(sketchKeys[count-1-position], position+1, true);}
		size=limit;
		ready=true;
		return size;
	}

	/** Actual live prefix; zero is a valid query with neither hybrid stream. */
	public int size(){
		require(ready, "Build a query order before reading its live length");
		return size;
	}

	/** Original model at a zero-based emitted position; tail scratch is never live. */
	public int modelAt(final int position){
		require(ready && position>=0 && position<size, "A caller position must lie inside the current emitted prefix");
		return order[position];
	}

	/** Which stream first emitted this model, not whether it also had native support. */
	public boolean sketchFirst(final int model){
		require(ready && model>=0 && model<ranks.length && ranks[model]>0,
			"First-emission provenance exists only for models in the current emitted prefix");
		return sketchFirst[model];
	}

	/** Original positive-native ranking with complete score ties, for cache parity. */
	public void nativeBest(final boolean[] acceptable, final int[] tuple){
		check(acceptable);
		for(int model=0; model<filtered.length; model++){filtered[model]=acceptable[model] && nativeScores[model]>0;}
		NcrnaDdlRanks.best(nativeScores, filtered, tuple);
	}

	/** Original sketch ranking with complete score ties, before interleaving. */
	public void sketchBest(final boolean[] acceptable, final int[] tuple){
		check(acceptable);
		NcrnaDdlRanks.best(sketchScores, acceptable, tuple);
	}

	/** Hybrid positions form a total order; its first acceptable model has a
	 * one-position tie interval. Empty acceptability writes four -1 fields. */
	public void hybridBest(final boolean[] acceptable, final int[] tuple){
		check(acceptable);
		require(tuple!=null && tuple.length>=4, "Hybrid results use rank/tie-first/tie-last/model fields");
		int selected=-1;
		for(int model=0; model<ranks.length; model++){
			if(acceptable[model] && ranks[model]>0 && (selected<0 || ranks[model]<ranks[selected])){selected=model;}
		}
		if(selected<0){Arrays.fill(tuple, 0, 4, -1);}
		else{
			tuple[0]=ranks[selected];
			tuple[1]=ranks[selected];
			tuple[2]=ranks[selected];
			tuple[3]=selected;
		}
	}

	/** One-based deduplicated position;0 means absent from a filtered/capped order. */
	public int rank(final int model){
		require(ready && model>=0 && model<ranks.length, "Build the complete order before reading an original model position");
		return ranks[model];
	}

	/** Copies the first distinct model positions into caller-owned reusable scratch.
	 * The destination length is the position cap, not a successful-alignment count. */
	public void copyFirst(final int[] destination){
		require(ready && destination!=null && destination.length>0 && destination.length<=size,
			"A hybrid shortlist needs a positive prefix of the actual emitted model order");
		System.arraycopy(order, 0, destination, 0, destination.length);
	}

	private void check(final boolean[] acceptable){
		require(ready && acceptable!=null && acceptable.length==ranks.length,
			"Build a complete query order before applying its locus-specific oracle flags");
	}
	private int emit(final long key, final int next, final boolean fromSketch){
		final int model=~(int)key;
		assert(model>=0 && model<ranks.length) : "Packed model "+model+" must be in [0,"+ranks.length+") before indexing hybrid scratch";
		if(seen[model]){return next;}
		seen[model]=true;
		sketchFirst[model]=fromSketch;
		ranks[model]=next;
		order[next-1]=model;
		return next+1;
	}
	private static long pack(final int score, final int model){
		assert(score>=0 && model>=0) : "Nonnegative score/model fields preserve signed packed-key order; score="+score+", model="+model;
		return ((long)score<<32)|(0xffffffffL-model);
	}
	/** In-place heap sorting keeps query ordering free of temporary sort arrays. */
	private static void sort(final long[] keys, final int length){
		assert(length>=0 && length<=keys.length) : "Live prefix "+length+" must fit key capacity "+keys.length+" before heap indexing";
		for(int root=(length>>>1)-1; root>=0; root--){sift(keys, root, length);}
		for(int end=length-1; end>0; end--){
			final long largest=keys[0];
			keys[0]=keys[end];
			keys[end]=largest;
			sift(keys, 0, end);
		}
	}
	private static void sift(final long[] keys, int root, final int end){
		assert(root>=0 && root<end && end<=keys.length) : "Heap root "+root+" must be inside [0,"+end+") with capacity "+keys.length+" before child indexing";
		final long value=keys[root];
		final int half=end>>>1;
		while(root<half){
			int child=(root<<1)+1;
			if(child+1<end && keys[child]<keys[child+1]){child++;}
			if(value>=keys[child]){break;}
			keys[root]=keys[child];
			root=child;
		}
		keys[root]=value;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final long[] nativeKeys, sketchKeys;
	private final float[] nativeScores, sketchScores;
	private final int[] ranks, order;
	private final boolean[] seen, filtered, sketchFirst;
	private int size;
	private boolean ready;
}
