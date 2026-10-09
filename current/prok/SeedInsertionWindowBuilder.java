package prok;

import java.util.ArrayList;
import java.util.Arrays;
import map.LongHashSet;
import structures.IntList;
import structures.LongList;

/** Proposal-only joining across one declared CM boundary. Reads no truth and
 * performs no alignment. Coordinates are on one sense-oriented contig. Diagonal
 * radius is an explicit geometric tolerance, not a learned confidence interval.
 * @author Brian Bushnell, Raiden */
final class SeedInsertionWindowBuilder {

	Result build(int[] centers, long[] keys, int seqLen, SeedModelPositionTable table,
			int sideSupport, int diagonalRadius, int endpointSlack, int minGap, int maxGap){
		require(centers!=null && keys!=null && centers.length==keys.length && seqLen>0 && table!=null,
			"Joined windows require a bound position table and the actual shared-scan key/center stream");
		require(sideSupport>0 && diagonalRadius>=0 && endpointSlack>=0 && minGap>0 && maxGap>=minGap,
			"Support and gap envelope must be explicit; tolerance cannot silently reverse coordinate intervals");
		final Result result=new Result();
		for(int i=0; i<centers.length; i++){
			require(centers[i]>=8 && (long)centers[i]+8<seqLen && (i==0 || centers[i]>centers[i-1]),
				"One oriented K17 scan supplies ordered unique centers with17 real genomic bases");
			final SeedModelPositionTable.Entry e=table.get(keys[i]);require(e!=null, "Every scanned key must belong to the full position resource");
			if(e.size()==0){result.unmappedHits++;}else{result.mappedHits++;}
		}
		for(int model=0; model<table.modelCount(); model++){
			for(int site=0; site<table.siteCount(); site++){
				final ArrayList<Group> left=groups(centers, keys, table, model, site, -1, diagonalRadius);
				final ArrayList<Group> right=groups(centers, keys, table, model, site, 1, diagonalRadius);
				final SiteResult sr=new SiteResult(model, site, left, right);result.sites.add(sr);
				for(Group l:left){
					if(l.support()<sideSupport){continue;}
					for(Group r:right){
						if(r.support()<sideSupport){continue;}
						final long gapLow=r.low-l.high, gapHigh=r.high-l.low;
						if(gapHigh<minGap || gapLow>maxGap){sr.outOfEnvelopePairs++;continue;}
						if((long)l.maxX+SeedModelPositionTable.K>r.minX || l.maxP>=r.minP){sr.unorderedPairs++;continue;}
						if(l.support()==1 && r.support()==1 && l.keys[0]==r.keys[0]){sr.sameKeyPairs++;continue;}
						final boolean tandem=completePartner(l, right, sideSupport) && completePartner(r, left, sideSupport);
						final boolean reset=modelReset(l, r, left, right, sideSupport);
						final long rawStart=Math.min(l.low, l.minX)-endpointSlack;
						final long rawStop=Math.max(r.high+table.modelLength(model)-1L, r.maxX+16L)+endpointSlack;
						final long start=Math.max(0L, rawStart), stop=Math.min(seqLen-1L, rawStop);
						require(start<=stop && start<=l.minX && stop>=r.maxX+16L,
							"A joined outer envelope must retain both complete anchor groups after physical clipping");
						result.proposals.add(new Proposal(model, site, l, r, gapLow, gapHigh, rawStart, rawStop,
							(int)start, (int)stop, tandem, reset));sr.pairs++;
					}
				}
			}
		}
		return result;
	}
	/** Group by an intersection shared by EVERY member, never single-linkage.
	 * Exact alternative placements survive in separate groups or as an explicit
	 * repeated-position flag. Distinct seed keys, not occurrences, supply support. */
	private ArrayList<Group> groups(int[] centers, long[] keys, SeedModelPositionTable table,
			int model, int site, int side, int radius){
		assert(side==-1 || side==1):"Crossing seeds cannot anchor either side of a boundary";
		x.clear();p.clear();key.clear();ambiguous.clear();
		for(int i=0; i<centers.length; i++){
			final SeedModelPositionTable.Entry e=table.get(keys[i]);
			for(int j=0; j<e.size(); j++){
				final SeedModelPositionTable.Position pos=e.get(j);
				if(pos.model!=model || table.side(pos, site)!=side){continue;}
				final boolean repeat=j>0 && e.get(j-1).model==model || j+1<e.size() && e.get(j+1).model==model;
				x.add(centers[i]-8);p.add(pos.start);key.add(keys[i]);
				ambiguous.add(repeat || pos.rfStartLow!=pos.rfStartHigh || pos.rfEndLow!=pos.rfEndHigh ? 1 : 0);
			}
		}
		final int n=x.size;if(order.length<n){order=new long[Math.max(n, 16)];}
		for(int i=0; i<n; i++){order[i]=((long)(x.get(i)-p.get(i))<<32)|(i&0xffffffffL);}
		Arrays.sort(order, 0, n);final ArrayList<Group> out=new ArrayList<Group>();
		if(n==0){return out;}
		int first=0;long lo=(order[0]>>32)-radius, hi=(order[0]>>32)+radius;
		for(int at=1; at<n; at++){
			final long d=order[at]>>32, nextLo=Math.max(lo, d-radius), nextHi=Math.min(hi, d+radius);
			if(nextLo<=nextHi){lo=nextLo;hi=nextHi;}
			else{out.add(group(first, at, lo, hi));first=at;lo=d-radius;hi=d+radius;}
		}
		out.add(group(first, n, lo, hi));return out;
	}
	private Group group(int from, int to, long low, long high){
		assert(from>=0 && from<to && to<=x.size && low<=high):"A seed group retains a nonempty membership and common diagonal interval";
		final LongHashSet distinct=new LongHashSet(Math.max(4, (to-from)*2));
		final int[] xs=new int[to-from], ps=new int[to-from];final long[] ks=new long[to-from];
		int minX=Integer.MAX_VALUE, maxX=-1, minP=Integer.MAX_VALUE, maxP=-1;boolean uncertain=false;
		for(int at=from; at<to; at++){
			final int i=(int)order[at], j=at-from;xs[j]=x.get(i);ps[j]=p.get(i);ks[j]=key.get(i);distinct.add(ks[j]);
			minX=Math.min(minX, xs[j]);maxX=Math.max(maxX, xs[j]);minP=Math.min(minP, ps[j]);maxP=Math.max(maxP, ps[j]);uncertain|=ambiguous.get(i)!=0;
		}
		final long[] unique=distinct.toArray();Arrays.sort(unique);
		return new Group(low, high, minX, maxX, minP, maxP, unique, xs, ps, ks, uncertain);
	}
	private static boolean completePartner(Group g, ArrayList<Group> opposite, int support){
		assert(support>0):"Tandem safeguards use the same declared side-support threshold as proposals";
		for(Group other:opposite){
			if(other.low>g.high){break;}
			if(other.support()<support || other.high<g.low){continue;}
			final Group left=g.minP<other.minP ? g : other, right=left==g ? other : g;
			if(left.maxP<right.minP && (long)left.maxX+17<=right.minX
				&& (left.support()>1 || right.support()>1 || left.keys[0]!=right.keys[0])){return true;}
		}
		return false;
	}
	/** A supported right-side cluster followed by a supported left-side cluster
	 * BETWEEN the proposed flanks exposes a model-position reset. Keep the
	 * hypothesis flagged for diagnostic counting, rather than emit it as a gene. */
	private static boolean modelReset(Group a, Group z, ArrayList<Group> left, ArrayList<Group> right, int support){
		assert(a.maxX<z.minX):"Reset inspection is limited to ordered proposal flanks";
		for(Group r:right){
			if(r.support()<support || r.minX<=a.maxX || r.maxX>=z.minX){continue;}
			for(Group l:left){
				if(l.support()>=support && l.minX>r.maxX && l.maxX<z.minX){return true;}
			}
		}
		return false;
	}
	static final class Group {
		Group(long lo, long hi, int x0, int x1, int p0, int p1, long[] k, int[] xs, int[] ps, long[] ks, boolean a){
			low=lo;high=hi;minX=x0;maxX=x1;minP=p0;maxP=p1;keys=k;genomicStarts=xs;modelStarts=ps;memberKeys=ks;ambiguous=a;
		}
		int support(){return keys.length;}
		final long low, high;final int minX, maxX, minP, maxP;final long[] keys, memberKeys;
		final int[] genomicStarts, modelStarts;final boolean ambiguous;
	}
	static final class Proposal {
		Proposal(int m, int s, Group a, Group z, long gl, long gh, long lo, long hi, int x, int y, boolean t, boolean r){
			model=m;site=s;left=a;right=z;gapLow=gl;gapHigh=gh;rawStart=lo;rawStop=hi;start=x;stop=y;tandem=t;reset=r;
		}
		boolean eligible(){return !tandem && !reset;}
		boolean ambiguous(){return left.ambiguous || right.ambiguous;}
		final int model, site, start, stop;final long gapLow, gapHigh, rawStart, rawStop;
		final Group left, right;final boolean tandem, reset;
	}
	static final class SiteResult {
		SiteResult(int m, int s, ArrayList<Group> a, ArrayList<Group> z){model=m;site=s;left=a;right=z;}
		final int model, site;final ArrayList<Group> left, right;
		long outOfEnvelopePairs, unorderedPairs, sameKeyPairs, pairs;
	}
	static final class Result {
		final ArrayList<Proposal> proposals=new ArrayList<Proposal>();
		final ArrayList<SiteResult> sites=new ArrayList<SiteResult>();int mappedHits, unmappedHits;
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private final IntList x=new IntList(), p=new IntList(), ambiguous=new IntList();
	private final LongList key=new LongList();private long[] order=new long[0];
}
