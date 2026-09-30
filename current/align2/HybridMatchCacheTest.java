package align2;

import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.atomic.AtomicLong;
import stream.SiteScore;

/** Ownership/miss tests and reused native-MSA scratch oracle.
 * With no arguments, runs bounded in-memory tests. Mapper-audit runs require
 * expectcache=t or expectcache=f as the first argument, followed by BBMapS
 * arguments. A positive expectation requires a fixture that actually hits the
 * cache; merely enabling a cache does not guarantee reuse on arbitrary reads.
 * Runs in a dedicated JVM because the mapper and audit hook use static state.
 * @author Collei */
public final class HybridMatchCacheTest{

	public static void main(final String[] args){
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Run with assertions enabled");
		if(args.length>0){
			final boolean expectCache;
			if(args[0].equalsIgnoreCase("expectcache=t")){expectCache=true;}
			else if(args[0].equalsIgnoreCase("expectcache=f")){expectCache=false;}
			else{throw new IllegalArgumentException("Mapper audit requires expectcache=t or expectcache=f before BBMapS arguments");}
			require(args.length>1, "Mapper audit requires BBMapS arguments after expectcache");
			require(HybridMatchCache.auditHits==null, "A mapper audit must own the global cache-hit hook exclusively");
			HybridMatchCache.auditHits=new AtomicLong();
			try{
				BBMapS.main(Arrays.copyOfRange(args, 1, args.length));
				final long hits=HybridMatchCache.auditHits.get();
				System.err.println("CACHE_AUDIT_HITS\t"+hits);
				// Public hybridmaxindel enables caching independently of the old research
				// flag (BBMapThread constructor). The fixture states its expected result.
				require(expectCache ? hits>0 : hits==0,
						"Cache audit expected "+(expectCache ? "positive hits" : "zero hits")+", observed "+hits);
			}finally{HybridMatchCache.auditHits=null;}
			return;
		}
		ownershipAndMisses();
		reusedWorker();
		System.out.println("HYBRID_MATCH_CACHE_TEST_PASS\ttrue");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Checks owned array snapshots, one-use hits, and rejection of changed inputs. */
	private static void ownershipAndMisses(){
		final byte[] bases={'A', 'C', 'G', 'T'}, minus={'A', 'C', 'G', 'T'};
		SiteScore site=site(4), result=site.clone(); result.match=new byte[]{'m', 'm', 'm', 'm'};
		final HybridMatchCache cache=new HybridMatchCache();
		cache.remember(0, 7, site, result, bases, minus, 200, 300);
		result.match[0]='S';// The result buffer is not cache-owned.
		require(!cache.restore(7, site, bases, minus, 200, 300, false), "Observation phase must not replay");
		cache.replaying=true;
		require(cache.restore(7, site, bases, minus, 200, 300, false), "Expected unchanged-site hit");
		require(site.match[0]=='m', "Captured result must own its bytes");
		site.match=null;
		require(!cache.restore(7, site, bases, minus, 200, 300, false), "Cache hit must be consumed");
		for(int strand=0; strand<2; strand++){
			cache.clear(); site=site(4); result=site.clone(); result.match=new byte[]{'m', 'm', 'm', 'm'};
			cache.remember(0, 7, site, result, bases, minus, 200, 300); cache.replaying=true;
			final byte[] changed=strand==0 ? bases : minus;
			final byte old=changed[0]; changed[0]='T';
			require(!cache.restore(7, site, bases, minus, 200, 300, false), "Captured bases must be owned snapshots for strand "+strand);
			changed[0]=old;
		}
		for(int failure=0; failure<8; failure++){
			cache.clear(); site=site(4); result=site.clone(); result.match=new byte[]{'m', 'm', 'm', 'm'};
			cache.remember(0, 7, site, result, bases, minus, 200, 300); cache.replaying=true;
			SiteScore target=site; byte[] query=bases; int maximum=300; long id=7; boolean secondary=false;
			switch(failure){
			case 0: site.start++; break;
			case 1: site.slowScore++; break;
			case 2: query=bases.clone(); query[0]='T'; break;
			case 3: maximum++; break;
			case 4: id++; break;
			case 5: secondary=true; break;
			case 6: target=site.clone(); break;
			case 7: cache.clear(); cache.replaying=true; break;
			default: throw new AssertionError(failure);
			}
			require(!cache.restore(id, target, query, minus, 200, maximum, secondary), "Unsafe cache hit for mutation "+failure);
		}
		cache.clear(); site=site(4); result=site.clone(); result.match=new byte[]{'m', 'm', 'm', 'm'};
		result.stop++; cache.remember(0, 7, site, result, bases, minus, 200, 300); cache.replaying=true;
		require(!cache.restore(7, site, bases, minus, 200, 300, false), "Changed observation coordinates must not be admitted");
		System.out.println("ownership_and_misses\tPASS");
	}
	/** Alternates alignments and widths to check scratch state after skipping a cached fill. */
	private static void reusedWorker(){
		final Random random=new Random(441);
		final byte[] reference=new byte[340];
		for(int i=0; i<reference.length; i++){reference[i]=(byte)"ACGT".charAt(random.nextInt(4));}
		final byte[] a=Arrays.copyOfRange(reference, 80, 200), b=Arrays.copyOfRange(reference, 20, 140),
				c=Arrays.copyOfRange(reference, 30, 150);
		a[50]=(byte)(a[50]=='A' ? 'C' : 'A'); b[30]=(byte)(b[30]=='G' ? 'T' : 'G');
		final MultiStateAligner11ts cached=new MultiStateAligner11ts(160, 420), ordinary=new MultiStateAligner11ts(160, 420);
		for(int round=0; round<12; round++){
			final int width=round%2==0 ? 150 : 330;
			final Result first=align(cached, a, reference, 70, 220, 6000), control=align(ordinary, a, reference, 70, 220, 6000);
			require(first.same(control), "Initial A mismatch");
			align(cached, b, reference, 0, width, round%3==0 ? 0 : 6000);
			align(ordinary, b, reference, 0, width, round%3==0 ? 0 : 6000);
			final SiteScore input=site(a.length), generated=input.clone(); generated.match=first.match;
			final HybridMatchCache memo=new HybridMatchCache();
			memo.remember(0, round, input, generated, a, null, 10000, 12000); memo.replaying=true;
			require(memo.restore(round, input, a, null, 10000, 12000, false), "Required A cache hit");
			final Result recomputed=align(ordinary, a, reference, 70, 220, 6000);
			require(Arrays.equals(input.match, recomputed.match) && first.same(recomputed), "Recomputed A differs");
			final Result afterHit=align(cached, c, reference, 20, 175, 5000), afterRecompute=align(ordinary, c, reference, 20, 175, 5000);
			require(afterHit.same(afterRecompute), "Subsequent C depends on skipped A at round "+round);
		}
		System.out.println("reused_worker_12_sequences\tPASS");
	}
	/** Returns score/endpoints and expanded native match for a nonempty oracle alignment. */
	private static Result align(final MultiStateAligner11ts msa, final byte[] query, final byte[] ref,
			final int start, final int end, final int minimum){
		final int[] max=msa.fillLimited(query, ref, start, end, minimum, null);
		require(max!=null, "Oracle fixture must yield an alignment");
		final int[] score=msa.score(query, ref, start, end, max[0], max[1], max[2], false);
		final byte[] match=msa.traceback(query, ref, start, end, max[0], max[1], max[2], false);
		return new Result(score, match);
	}
	/** Makes an isolated candidate whose scalar state is stable during observation. */
	private static SiteScore site(final int length){
		final SiteScore site=new SiteScore(1, (byte)0, 10, 10+length-1, 3, 100);
		site.slowScore=100; site.pairedScore=150; return site;
	}
	/** Keeps test failures loud even if an invocation accidentally disables assertions. */
	private static void require(final boolean pass, final String message){if(!pass){throw new AssertionError(message);}}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns the two arrays returned by an oracle alignment; compares their contents. */
	private static final class Result{
		Result(final int[] score_, final byte[] match_){score=score_; match=match_;}
		boolean same(final Result r){return Arrays.equals(score, r.score) && Arrays.equals(match, r.match);}
		final int[] score;
		final byte[] match;
	}
}
