package assemble;

import java.util.Arrays;
import structures.IntList;

/** Exhaustive interval geometry against independently enumerated window coverage.
 * @author Fischl */
public final class LocalEditTroughLocatorTest {
	public static void main(final String[] args){
		for(int k:new int[]{31,62,63}){
			for(int width=1;width<=k;width++){geometry(k,width);}
			boundaries(k);reuse(k);
		}
		System.out.println("LOCAL_EDIT_LOCATOR_TEST_OK checks="+checks+"; all trough widths1..K at K31/62/63, base/gap intersections, explicit RC coordinates, supported flanks, undefined/edge/wide/reset; geometry only.");
	}
	private static void geometry(final int k,final int width){
		final byte[] bases=sequence(3*k+17);final IntList counts=profile(bases.length,k);final int a=4,end=a+width;
		for(int j=a;j<end;j++){counts.set(j,j%3==0 ? -1 : j%3);}
		final LocalEditTroughLocator locator=new LocalEditTroughLocator(k,3);locator.reset(bases,counts);
		check(locator.next() && locator.depthStart==a && locator.depthEnd==end,"Expected exactly the bounded low-depth region.");
		for(int p=0;p<bases.length;p++){
			boolean inEvery=true;for(int j=a;j<end;j++){if(!(j<=p && p<j+k)){inEvery=false;}}
			check(inEvery==(p>=locator.baseStart && p<locator.baseEnd),"Base locus must equal direct intersection of every low kmer.");
		}
		for(int g=0;g<=bases.length;g++){
			boolean crossesEvery=true;for(int j=a;j<end;j++){if(!(j<g && g<j+k)){crossesEvery=false;}}
			check(crossesEvery==(g>=locator.gapStart && g<locator.gapEnd),"Gap locus must be strictly internal to every low kmer.");
		}
		final int bs=locator.baseStart,be=locator.baseEnd,gs=locator.gapStart,ge=locator.gapEnd;
		check(!locator.next() && locator.baseStart==-1,"Exhaustion must not expose the previous locus as a new result.");
		final IntList reverse=new IntList(counts.size);for(int j=counts.size-1;j>=0;j--){reverse.add(counts.get(j));}
		// Original all-A read has all-T RC; count-start indices reverse exactly.
		final byte[] rc=bases.clone();Arrays.fill(rc,(byte)'T');locator.reset(rc,reverse);
		check(locator.next(),"Mirrored profile must retain its bounded trough.");
		check(locator.baseStart==bases.length-be && locator.baseEnd==bases.length-bs,"Base coordinates reflect as n-1-p.");
		check(locator.gapStart==bases.length-ge+1 && locator.gapEnd==bases.length-gs+1,"Integer insertion boundaries reflect as n-g, not n-1-g.");
		if(width==k){check(be-bs==1 && ge==gs,"K-long trough localizes one base and no common internal gap.");}
		if(width==k-1){check(be-bs==2 && ge-gs==1,"K-1 trough localizes one insertion boundary, not one substitution base.");}
	}
	private static void boundaries(final int k){
		final byte[] bases=sequence(4*k+20);final IntList counts=profile(bases.length,k);
		final LocalEditTroughLocator locator=new LocalEditTroughLocator(k,3);
		for(int i=0;i<k;i++){counts.set(i,0);}locator.reset(bases,counts);
		check(!locator.next() && locator.skippedEdge==1,"Leading low region lacks a supported left anchor.");
		counts.clear();for(int i=0;i<bases.length-k+1;i++){counts.add(40);}
		for(int i=3;i<3+k+1;i++){counts.set(i,0);}locator.reset(bases,counts);
		check(!locator.next() && locator.skippedWide==1,"More than K low windows cannot localize a single base by intersection.");
		counts.set(3+k,3); // Exact threshold is a supported right flank.
		for(int npos:new int[]{2,3,3+k,3+2*k-1}){
			final byte[] invalid=bases.clone();invalid[npos]='N';locator.reset(invalid,counts);
			check(!locator.next() && locator.skippedUndefined==1,"Undefined bases in a trough or supporting window must defer localization.");
		}
		final byte[] outside=bases.clone();outside[0]='N';locator.reset(outside,counts);
		check(locator.next(),"An undefined base outside all evidence windows must not invalidate this locus.");
		counts.set(1,-2);locator.reset(bases,counts);boolean failed=false;
		try{locator.next();}catch(IllegalArgumentException e){failed=true;}check(failed,"Invalid depth must fail loudly, not become error evidence.");
	}
	private static void reuse(final int k){
		final byte[] bases=sequence(8*k);final IntList counts=profile(bases.length,k);
		for(int j=3;j<6;j++){counts.set(j,0);}for(int j=4*k;j<5*k;j++){counts.set(j,2);}
		bases[k]='N';final LocalEditTroughLocator locator=new LocalEditTroughLocator(k,3);locator.reset(bases,counts);
		check(locator.next() && locator.depthStart==4*k && locator.skippedUndefined==1,"Later trough must ignore an obsolete undefined-base witness.");
		check(!locator.next(),"Two-region profile has no third region.");
		locator.reset(new byte[0],new IntList());check(!locator.next(),"Empty read/profile must terminate cleanly after reuse.");
		boolean failed=false;try{locator.reset(new byte[k],new IntList());}catch(IllegalArgumentException e){failed=true;}
		check(failed,"Profile shape must match original read coordinates.");
		failed=false;try{locator.next();}catch(IllegalStateException e){failed=true;}
		check(failed,"Rejected reset must invalidate any previous profile.");
	}
	private static byte[] sequence(final int n){assert(n>=0) : "Test sequence length must be nonnegative.";final byte[] b=new byte[n];Arrays.fill(b,(byte)'A');return b;}
	private static IntList profile(final int n,final int k){assert(n>=k) : "This fixture requires at least one full window.";final IntList list=new IntList(n-k+1);for(int j=0;j<n-k+1;j++){list.add(40);}return list;}
	private static void check(final boolean ok,final String why){checks++;if(!ok){throw new AssertionError(why);}}
	private static long checks;
}
