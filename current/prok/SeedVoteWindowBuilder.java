package prok;

import java.util.ArrayList;
import java.util.Arrays;
import structures.IntList;

/** Groups votes with compatible endpoint ranges, then votes and pads.
 * Scratch is reused per scavenger. Only returned windows allocate per result.
 * @author Raiden
 */
final class SeedVoteWindowBuilder {
	ArrayList<int[]> build(int[] centers, long[] keys, int seqLen, SeedOffsetTable table,
			int k, int minLen, int fallbackPad, int slack, float fraction, int fallbackCoreLength){
		return build(centers, keys, seqLen, table, k, minLen, fallbackPad, slack, fraction, fallbackCoreLength, null);
	}
	/** The optional observer receives only executed groups, before claimed-region subtraction. */
	ArrayList<int[]> build(int[] centers, long[] keys, int seqLen, SeedOffsetTable table,
			int k, int minLen, int fallbackPad, int slack, float fraction, int fallbackCoreLength, NcrnaVoteDiagnostics trace){
		require(centers!=null && keys!=null && centers.length==keys.length && table!=null && table.k()==k,
			"Vote-first windows require the bound seed table and one key per position");
		require(k>0 && (k&1)==1 && seqLen>0 && minLen>0 && fallbackPad>=0 && slack>=0
			&& Float.isFinite(fraction) && fraction>0 && fraction<=1 && fallbackCoreLength>0,
			"Window geometry requires odd K, positive lengths and a core-overlap fraction in (0,1]");
		ensure(centers.length);fallback.clear();int n=0;voted=0;fallbackWindows=0;
		for(int i=0; i<centers.length; i++){
			final int center=centers[i];require(center>=k/2 && (long)center+k/2<seqLen && (i==0 || center>centers[i-1]),
				"Seed centers must be unique, ordered and contain their full K bases");
			final SeedOffsetTable.KmerInfo info=table.get(keys[i]);require(info!=null, "Every scanned seed must belong to its bound offset table");
			if(!info.trained()){fallback.add(center);continue;}
			left[n]=center-k/2-info.leftOffset;right[n]=center+k/2+info.rightOffset;
			leftWeight[n]=NcrnaScavenger.voteWeight(info.leftSD, info.leftOffset);
			rightWeight[n]=NcrnaScavenger.voteWeight(info.rightSD, info.rightOffset);original[n]=i;
			require(Double.isFinite(left[n]) && Double.isFinite(right[n]) && left[n]>=Integer.MIN_VALUE && left[n]<=Integer.MAX_VALUE
				&& right[n]>=left[n] && leftWeight[n]>0 && rightWeight[n]>0,
				"Finite ordered predictions and positive weights are required before coordinate sorting");
			// W2-D25 weights already downweight uncertain, distant endpoint predictions.
			// Their inverse square roots provide a dispersion scale, NOT a confidence interval.
			// Split the existing geometric tolerance between both sides of each endpoint.
			final double tolerance=(1.0-fraction)*(right[n]-left[n]+1)*0.5;
			leftRadius[n]=tolerance+Math.sqrt(1.0/leftWeight[n]);
			rightRadius[n]=tolerance+Math.sqrt(1.0/rightWeight[n]);
			// Rounded-left ties retain original seed order; negative predictions sort before positive ones.
			order[n]=(Math.round(left[n])<<32)|(n&0xffffffffL);n++;
		}
		Arrays.sort(order, 0, n);final ArrayList<int[]> out=new ArrayList<int[]>();
		if(n>0){int first=0, index=(int)order[0];
			double leftLo=left[index]-leftRadius[index], leftHi=left[index]+leftRadius[index];
			double rightLo=right[index]-rightRadius[index], rightHi=right[index]+rightRadius[index];
			for(int j=1; j<n; j++){
				index=(int)order[j];
				final double nextLeftLo=Math.max(leftLo, left[index]-leftRadius[index]);
				final double nextLeftHi=Math.min(leftHi, left[index]+leftRadius[index]);
				final double nextRightLo=Math.max(rightLo, right[index]-rightRadius[index]);
				final double nextRightHi=Math.min(rightHi, right[index]+rightRadius[index]);
				// Every member must admit a common ordered endpoint pair. Intersections only
				// shrink: a noisy bridging seed can never reconnect incompatible precise loci.
				if(nextLeftLo<=nextLeftHi && nextRightLo<=nextRightHi && nextLeftLo<=nextRightHi){
					leftLo=nextLeftLo;leftHi=nextLeftHi;rightLo=nextRightLo;rightHi=nextRightHi;
				}else{
					emit(first, j, centers, seqLen, minLen, slack, out, trace);first=j;
					leftLo=left[index]-leftRadius[index];leftHi=left[index]+leftRadius[index];
					rightLo=right[index]-rightRadius[index];rightHi=right[index]+rightRadius[index];
				}
			}
			emit(first, n, centers, seqLen, minLen, slack, out, trace);
		}
		fallback.sort();int start=-1, stop=-1, firstFallback=0;
		final double fallbackOverlap=fraction*Math.min(fallbackCoreLength, seqLen);
		for(int i=0; i<fallback.size; i++){
			final int center=fallback.get(i), lo=(int)Math.max(0L, (long)center-fallbackPad);
			final int hi=(int)Math.min(seqLen-1L, (long)center+fallbackPad+k);
			if(start<0){start=lo;stop=hi;continue;}
			final int a=Math.max(start, lo), z=Math.min(stop, hi);
			if((long)z-a+1>=fallbackOverlap){start=a;stop=z;}
			else{emitFallback(start, stop, minLen, out, trace, firstFallback, i);firstFallback=i;start=lo;stop=hi;}
		}
		if(start>=0){emitFallback(start, stop, minLen, out, trace, firstFallback, fallback.size);}
		out.sort((a,b)->a[0]!=b[0]?Integer.compare(a[0],b[0]):a[1]!=b[1]?Integer.compare(a[1],b[1]):Integer.compare(a[2],b[2]));
		assert(out.size()==voted+fallbackWindows):"Window provenance counters must conserve returned candidate windows";
		return out;
	}
	private void emit(int from, int to, int[] centers, int seqLen, int minLen, int slack, ArrayList<int[]> out, NcrnaVoteDiagnostics trace){
		assert(from>=0 && to>from && to<=centers.length):"Each emitted vote group owns a nonempty sorted seed range";
		double l=0, r=0, wl=0, wr=0;
		for(int j=from; j<to; j++){final int i=(int)order[j];l+=left[i]*leftWeight[i];wl+=leftWeight[i];r+=right[i]*rightWeight[i];wr+=rightWeight[i];}
		require(wl>0 && wr>0 && Double.isFinite(l) && Double.isFinite(r), "Weighted coordinates must be finite before clipping");
		final long start=Math.max(0L, Math.round(l/wl)-slack), stop=Math.min(seqLen-1L, Math.round(r/wr)+slack);
		if(trace!=null){trace.core(from, to, order, original, left, right, l/wl, r/wr, wl, wr, start, stop, stop-start+1>=minLen);}
		if(stop-start+1>=minLen){out.add(new int[]{(int)start,(int)stop,NcrnaScavenger.WINDOW_VOTED});voted++;}
		else{for(int j=from; j<to; j++){fallback.add(centers[original[(int)order[j]]]);}}
	}
	private void emitFallback(int start, int stop, int minLen, ArrayList<int[]> out, NcrnaVoteDiagnostics trace, int from, int to){
		assert(start>=0 && stop>=start):"Fallback intersection must retain valid inclusive source bounds";
		if(trace!=null){trace.fallback(fallback, from, to, start, stop, stop-start+1>=minLen);}
		if(stop-start+1>=minLen){out.add(new int[]{start,stop,NcrnaScavenger.WINDOW_FALLBACK});fallbackWindows++;}
	}
	private void ensure(int n){
		assert(n>=0):"Scratch capacity follows an actual seed-array length";
		if(n<=order.length){return;}final int size=Math.max(n, Math.max(16, order.length*2));
		order=new long[size];original=new int[size];left=new double[size];right=new double[size];leftWeight=new double[size];rightWeight=new double[size];
		leftRadius=new double[size];rightRadius=new double[size];
	}
	static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	int voted,fallbackWindows;
	private final IntList fallback=new IntList();
	private long[] order=new long[0];private int[] original=new int[0];
	private double[] left=new double[0],right=new double[0],leftWeight=new double[0],rightWeight=new double[0];
	private double[] leftRadius=new double[0],rightRadius=new double[0];
}
