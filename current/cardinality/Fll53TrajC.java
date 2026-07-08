package cardinality;

/**
 * Compact chronicle recorder (v5, Brian's priority 2026-07-05): the honest
 * production form of Fll53Chron.Traj.  Traj keeps full trajectory VECTORS
 * (~2KB of doubles per sketch) but the fusion features only ever read the
 * last four adds-clock lags, the first checkpoint, running aggregates, the
 * last two promotions, and the last MLE witness — a fixed set of rolling
 * registers.  Point values are stored as 16-bit floats (Half); running sums
 * stay float32 (fp16 accumulation drifts); adds-clock checkpoint TIMES are
 * not stored at all, because ticks fire at exactly firstCk*2^i (adds
 * increments by one, so the crossing is exact) and la[i] is derivable from
 * the count.
 *
 * Persistent state: 11 shorts + 3 floats + 3 bytes = 37 bytes/sketch
 * (STATE_BYTES=37), vs ~2144 for Traj — memory that buys organisms instead.
 *
 * fuseC() reproduces Fll53Chron.fuse()'s 36 features in the same order;
 * TrajCAblate measures the fp16 feature/estimate drift side by side.
 *
 * @author Amber (design pressure: Brian — "it's a toy until it's efficient")
 * @date July 2026
 */
public class Fll53TrajC {

	static final double LOG2=Math.log(2.0);
	/** Honest persistent bytes per sketch (11*2 + 3*4 + 3).  yM1 is float32,
	 * not fp16: mSlope reconstructs dn=nHat-2^yM1*a1, a difference of large
	 * numbers where fp16's 2^-11 relative error can flip a small dn (measured:
	 * feature-27 drift up to 4.15 in log2 space; float32 kills it). */
	public static final int STATE_BYTES=37;
	static final byte G_UNSET=-128;

	final ComplexityFarm.BigNet base;
	final Fll53MLE mle;            // independent complexity witness (may be null)
	final long firstCk;            // derivable from B in production (4*B)

	// Persistent registers (the 37 bytes).
	short y1, y2, y3, y4;          // last four adds-clock estimates (newest first)
	short yF;                      // first checkpoint estimate
	short mn, mx;                  // trajectory envelope
	float ySum;                    // running sum of recorded estimates
	short yE1, yE2;                // estimates at the last two promotions
	short laE1, laE2;              // log2-adds at the last two promotions
	float yESum;                   // running sum over promotions
	float yM1;                     // MLE witness at the last checkpoint (float32: see STATE_BYTES)
	byte n=0, nE=0;                // clock counts (capped 64/40 like Traj)
	byte lastG=G_UNSET;

	// Transient (derivable: next = firstCk<<n while n<64).
	long next;

	public Fll53TrajC(ComplexityFarm.BigNet base, Fll53MLE mle, long firstCk){
		assert(firstCk>0) : "firstCk="+firstCk;
		this.base=base; this.mle=mle; this.firstCk=firstCk; next=firstCk;
	}

	/** Called after every add; fires either clock when due.  Mirrors
	 * Fll53Chron.Traj.tick exactly, modulo fp16 storage. */
	public void tick(Fll53 f, long adds){
		if(adds>=next){
			if(n<64){
				final float y=(float)Fll53Chron.yNow(base, f, adds);
				final float yM=(float)Fll53Chron.yMle(mle, f, adds);
				y4=y3; y3=y2; y2=y1; y1=Half.f2h(y);
				if(n==0){yF=y1; mn=y1; mx=y1;}
				else{
					if(y<Half.h2f(mn)){mn=Half.f2h(y);}
					if(y>Half.h2f(mx)){mx=Half.f2h(y);}
				}
				ySum+=Half.h2f(y1);   // sum what was STORED, so mean matches reads
				yM1=yM;
				n++;
			}
			next*=2;
			if(lastG==G_UNSET){lastG=(byte)f.getGlobalExp();}
		}
		final int g=f.getGlobalExp();
		if(lastG!=G_UNSET && g!=lastG){
			lastG=(byte)g;
			if(nE<40){
				final float la=(float)(Math.log(adds)/LOG2);
				final float y=(float)Fll53Chron.yNow(base, f, adds);
				yE2=yE1; yE1=Half.f2h(y);
				laE2=laE1; laE1=Half.f2h(la);
				yESum+=Half.h2f(yE1);
				nE++;
			}
		}
	}

	/** log2 of the i-th adds-clock checkpoint time — exact, no storage. */
	double laAt(int i){return Math.log(firstCk)/LOG2+i;}

	/** 36 fusion features, same order as Fll53Chron.fuse(). */
	public static double[] fuseC(Fll53TrajC tr, Fll53 f, long adds, int B){
		final double now=Fll53Chron.yNow(tr.base, f, adds);
		final double la=Math.log(adds)/LOG2;
		final int n=tr.n;
		final double y1=(n>0) ? Half.h2f(tr.y1) : now;
		final double y2=(n>1) ? Half.h2f(tr.y2) : y1;
		final double y3=(n>2) ? Half.h2f(tr.y3) : y2;
		final double y4=(n>3) ? Half.h2f(tr.y4) : y3;
		final double yF=(n>0) ? Half.h2f(tr.yF) : now;
		final double laF=(n>0) ? tr.laAt(0) : la;
		final double la1=(n>0) ? tr.laAt(n-1) : la;
		final double mn=(n>0) ? Half.h2f(tr.mn) : now;
		final double mx=(n>0) ? Half.h2f(tr.mx) : now;
		final double mean=(n>0) ? tr.ySum/n : now;
		final int nE=tr.nE;
		final double yE1=(nE>0) ? Half.h2f(tr.yE1) : now;
		final double yE2=(nE>1) ? Half.h2f(tr.yE2) : yE1;
		final double laE1=(nE>0) ? Half.h2f(tr.laE1) : la;
		final double laE2=(nE>1) ? Half.h2f(tr.laE2) : laE1;
		final double meanE=(nE>0) ? tr.yESum/nE : now;
		final double yMnow=Fll53Chron.yMle(tr.mle, f, adds);
		final double yM1=(n>0) ? tr.yM1 : yMnow;
		double mSlope=yMnow;
		if(tr.mle!=null && n>0){
			final double a1=Math.pow(2, la1);
			if(adds>a1+1){
				final double dn=tr.mle.estimate(f)-Math.pow(2, yM1)*a1;
				mSlope=Math.max(-8, Math.min(0.5,
					Math.log(Math.max(1e-3, dn/(adds-a1)))/LOG2));
			}
		}
		final double[] raw=new double[16];
		{
			int k=0;
			for(int w0=0; w0+Fll53Chron.WIN<=f.getNumOrgs(); w0+=Fll53Chron.WIN){
				final double[] wf=Fll53FarmW.features(f, adds, w0);
				for(int i=0; i<16; i++){raw[i]+=wf[i];}
				k++;
			}
			for(int i=0; i<16; i++){raw[i]/=Math.max(1, k);}
		}
		return new double[]{
			la-Math.log(B)/LOG2,
			now,
			y1, y2, y3, y4,
			now-y1,
			y1-y3,
			now-2*y1+y2,
			la-la1,
			mn, mx, mean,
			n/32.0,
			yF,
			laF-Math.log(B)/LOG2,
			yE1, yE2,
			la-laE1,
			laE1-laE2,
			nE/16.0,
			meanE,
			now-yE1,
			f.getGlobalExp()/32.0,
			yMnow,
			yM1,
			yMnow-yM1,
			mSlope,
			raw[1], raw[2], raw[3],
			raw[5], raw[6], raw[7],
			raw[4], raw[15]};
	}
}
