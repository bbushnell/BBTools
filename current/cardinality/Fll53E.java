package cardinality;

import dna.Data;

/**
 * Fll53E: production facade for the FLL53 estimator stack — the exact lane
 * that won the equal-total-bytes crowns at 512B and 1KB (see
 * fll53_chron_report.md).  Registered in CardinalityTracker.makeTracker as
 * loglogType=fll53.  Wraps:
 *   Fll53 sketch  +  Fll53MLE likelihood  +  window complexity net (blend)
 *   +  Fll53TrajC chronicle (37B trajectory state)  +  fusion net
 *   +  v5a envelope guard (chron is NEVER used unguarded: the fusion net
 *      exploded +60% beyond its training envelope on a real 1.3B-add
 *      stream; out-of-envelope falls back to the snapshot blend).
 *
 * Resource files (per Fll52's Data.findPath convention, falling back to
 * plain working-directory names; tables may be .tsv or .tsv.gz):
 *   fll53mle_<orgs>o_deep.tsv[.gz]  likelihood tables (192/384/768/1536)
 *   netwin53.txt                    window complexity net
 *   netchron53d_b33.txt             chronicle fusion net
 *   netchron53d_b33.txt.env         its training envelope (guard; REQUIRED
 *                                   for the chron lane to activate)
 *
 * Degradation ladder (each step warns once on stderr):
 *   full stack -> no fusion net or no .env: snapshot blend
 *              -> no window net: pure MLE (duplicate-blind!)
 *              -> no MLE table: construction fails (nothing to stand on).
 *
 * Memory: orgs*8/3 sketch bytes + Fll53TrajC.STATE_BYTES (37) + O(1).
 * No merge support (chronicle state is order-sensitive): add(tracker)
 * throws, like Fll52 — restrict to one thread.
 *
 * @author Amber (Brian Bushnell's design review)
 * @date July 2026
 */
public final class Fll53E extends CardinalityTracker {

	/** Supported organism counts (one likelihood table each). */
	public static final int[] SIZES={192, 384, 768, 1536};

	private final Fll53 fll;
	private final int orgs, fllBuckets;
	private final Fll53MLE mle;
	private final ComplexityFarm.BigNet win;    // null -> pure MLE
	private final ComplexityFarm.BigNet fuse;   // null -> blend only
	private final double[][] env;               // non-null iff fuse!=null
	private final Fll53TrajC tr;                // null iff win==null
	private long adds=0;

	Fll53E(){this(1920, 31, -1, 0);}   // 384 organisms = 1KB default

	Fll53E(parse.Parser p){
		this(p.loglogbuckets, p.loglogk, p.loglogseed, p.loglogMinprob);
	}

	/** sizeBytes rounded to the nearest supported class:
	 * 512B/1KB/2KB/4KB (192/384/768/1536 organisms). */
	public Fll53E(int sizeBytes){
		this(sizeBytes*3/8*Fll53.BUCKETS_PER_WORD, 31, -1, 0);
	}

	/** buckets_ (organisms*5, the factory convention) is rounded to the
	 * nearest supported organism class. */
	Fll53E(int buckets_, int k_, long seed, float minProb_){
		super(64, k_, seed, minProb_);   // own state unused; the wrapped sketch holds everything
		int best=SIZES[0];
		final int wantOrgs=buckets_/Fll53.BUCKETS_PER_WORD;
		for(int o : SIZES){
			if(Math.abs(o-wantOrgs)<Math.abs(best-wantOrgs)){best=o;}
		}
		orgs=best;
		fllBuckets=orgs*Fll53.BUCKETS_PER_WORD;
		fll=new Fll53(fllBuckets, k_, seed, minProb_);

		String tablePath=find("fll53mle_"+orgs+"o_deep.tsv");
		if(tablePath==null){tablePath=find("fll53mle_"+orgs+"o_deep.tsv.gz");}
		if(tablePath==null){
			throw new RuntimeException("Fll53E: likelihood table fll53mle_"
				+orgs+"o_deep.tsv[.gz] not found; cannot construct an estimator");
		}
		Fll53MLE mle0;
		try{mle0=new Fll53MLE(tablePath);}
		catch(Exception e){throw new RuntimeException("Fll53E: bad table "+tablePath, e);}
		mle=mle0;

		win=loadNet(find("netwin53.txt"), "window net (falling back to pure "
			+"MLE: duplicate-heavy streams WILL overcount)");
		final String fusePath=(win==null) ? null : find("netchron53d_b33.txt");
		final ComplexityFarm.BigNet fuse0=loadNet(fusePath,
			"chronicle net (falling back to snapshot blend)");
		final double[][] env0=(fuse0==null) ? null
			: Fll53Chron.loadEnv(fusePath+".env");
		if(fuse0!=null && env0==null){
			System.err.println("WARNING: Fll53E: "+fusePath+".env not found; "
				+"chron lane disabled (never run unguarded), using blend");
		}
		fuse=(env0==null) ? null : fuse0;
		env=env0;
		tr=(win==null) ? null : new Fll53TrajC(win, mle, 4L*fllBuckets);
	}

	private static String find(String name){
		final String p=Data.findPath("?"+name, false);
		if(p!=null && new java.io.File(p).isFile()){return p;}
		return new java.io.File(name).isFile() ? name : null;
	}

	private static ComplexityFarm.BigNet loadNet(String path, String warn){
		if(path==null){
			System.err.println("WARNING: Fll53E: missing "+warn);
			return null;
		}
		try{return ComplexityFarm.BigNet.load(path);}
		catch(Exception e){
			System.err.println("WARNING: Fll53E: unreadable "+warn+" ["+e+"]");
			return null;
		}
	}

	/** Per-element hook: all external adds arrive here via
	 * CardinalityTracker.add(number); the chronicle ticks here. */
	@Override
	public void hashAndStore(final long number){
		fll.hashAndStore(number);
		adds++;
		if(tr!=null){tr.tick(fll, adds);}
	}

	@Override
	public Fll53E copy(){return new Fll53E(fllBuckets, k, -1, minProb);}

	/** No merge support (chronicle state is order-sensitive); restrict to
	 * one thread, like Fll52. */
	@Override
	public void add(CardinalityTracker log){throw new UnsupportedOperationException(
		"FLL53 cannot merge; use threads=1");}

	@Override
	public float[] compensationFactorLogBucketsArray(){return null;}

	/** Current distinct-count estimate: the guarded chronicle lane when the
	 * full stack is loaded, degrading per the class javadoc. */
	@Override
	public long cardinality(){
		if(adds<1){return 0;}
		final double eM=mle.estimate(fll);
		// Cold sketch: the likelihood lane covers the low range on its own,
		// and the window nets have no window to read.
		if(win==null || eM<32.0*fllBuckets){return finish(eM);}
		double ySum=0;
		int k2=0;
		for(int w0=0; w0+Fll53FarmW.WIN<=orgs; w0+=Fll53FarmW.WIN){
			ySum+=win.predict(Fll53FarmW.features(fll, adds, w0));
			k2++;
		}
		if(k2==0){return finish(eM);}
		final double cHat=Math.pow(2.0, ySum/k2);
		final double w=Math.max(0, Math.min(1, (cHat-0.85)/0.10));
		final double blend=w*eM+(1-w)*adds*cHat;
		if(fuse==null){return finish(blend);}
		final double[] ff=Fll53TrajC.fuseC(tr, fll, adds, fllBuckets);
		if(!Fll53Chron.inEnvelope(env, ff, 0.05)){return finish(blend);}
		final double cC=Math.pow(2.0, fuse.predict(ff));
		final double wc=Math.max(0, Math.min(1, (cC-0.85)/0.10));
		return finish(wc*eM+(1-wc)*adds*cC);
	}

	/** House convention: estimates never exceed the add count. */
	private long finish(double est){
		long c=Math.round(est);
		if(clampToAdded){c=Math.min(c, adds);}
		return Math.max(0, c);
	}

	/** Stream-complexity diagnostic ĉ = distinct/adds in (0,1]; NaN before
	 * the sketch is warm enough to read. */
	public double complexity(){
		if(win==null || adds<4L*fllBuckets){return Double.NaN;}
		double ySum=0;
		int k2=0;
		for(int w0=0; w0+Fll53FarmW.WIN<=orgs; w0+=Fll53FarmW.WIN){
			ySum+=win.predict(Fll53FarmW.features(fll, adds, w0));
			k2++;
		}
		return (k2==0) ? Double.NaN : Math.pow(2.0, ySum/k2);
	}

	public long adds(){return adds;}
	public int sizeBytes(){return orgs*8/3+Fll53TrajC.STATE_BYTES;}
}
