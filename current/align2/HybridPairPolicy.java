package align2;

import java.util.ArrayList;
import dna.Data;
import stream.Read;
import stream.SamLine;
import stream.SiteScore;

/** Computational policy extracted from the qualified paired retry experiment.
 * Does not choose or commit an attempt, change index limits, or emit diagnostics.
 * @author Collei */
public final class HybridPairPolicy {

	private HybridPairPolicy(){}
	public static final int NO_SITE=1, HALF_ERRORS=2, UNKNOWN=4;

	/** Generate disposable top candidates so observation preserves owned sites. */
	public static Observation observe(BBMapThread owner,Read first,byte[] minus1,byte[] minus2,
			int imperfect1,int max1,int imperfect2,int max2,boolean geometryOnly){
		assert(first!=null && first.mate!=null) : "Paired retry observation requires both reads";
		final Read second=first.mate;
		final SiteScore a=generate(owner,first,minus1,imperfect1,max1,0);
		final SiteScore b=generate(owner,second,minus2,imperfect2,max2,1);
		return new Observation(geometryOnly ? 0 : flags(a,first.length())|flags(b,second.length()),
				geometryOnly ? geometry(a,first.length(),b,second.length()) : "NOT_RUN");
	}

	static SiteScore generate(BBMapThread owner,Read read,byte[] minus,int imperfect,int max){
		return generate(owner,read,minus,imperfect,max,-1);
	}
	private static SiteScore generate(BBMapThread owner,Read read,byte[] minus,int imperfect,int max,int slot){
		if(slot>=0 && owner.hybridMatchCache!=null){owner.hybridMatchCache.clear(slot);}
		final SiteScore original=read.topSite();
		if(original==null){return null;}
		final SiteScore copy=original.clone();
		copy.gaps=original.gaps==null ? null : original.gaps.clone();
		copy.match=original.match==null ? null : original.match.clone();
		owner.genMatchStringForSite(read.numericID,copy,read.bases,minus,imperfect,max,read.mate,false);
		if(slot>=0 && owner.hybridMatchCache!=null){
			owner.hybridMatchCache.remember(slot,read.numericID,original,copy,read.bases,minus,imperfect,max);
		}
		return copy;
	}

	static int flags(SiteScore site,int length){
		if(site==null){return NO_SITE;}
		if(site.match==null){return UNKNOWN;}
		final int[] counts=new int[5];
		count(site.match,length,!site.plus(),counts);
		if(counts[4]>0){return UNKNOWN;}
		return signal(counts[2],counts[3],length) ? HALF_ERRORS : 0;
	}

	/** Counters: substitutions by half, all discrepancies by half, unknown bases. */
	static void count(byte[] match,int length,boolean reverse,int[] out){
		assert(length>=2 && out.length==5) : "NativeHalfSignalProbe contract: two nonempty halves and five counters";
		int q=0;
		for(byte op:match){
			assert("msSVIDXYCNB".indexOf((char)op)>=0) : "Unsupported native operation in paired retry observation: "+op;
			if(op=='D'){continue;}
			assert(q<length) : "Native match cannot consume more than the query length";
			final int position=reverse ? length-1-q : q,half=position<length/2 ? 0 : 1;
			if(op=='S' || op=='V'){out[half]++;}
			if("SVIXYC".indexOf((char)op)>=0){out[2+half]++;}
			if(op=='N' || op=='B'){out[4]++;}
			q++;
		}
		assert(q==length) : "Native match must consume exactly the query for half-read classification";
	}

	static boolean signal(int a,int b,int length){
		assert(length>=2) : "Half-error ratios require two nonempty query halves";
		return a+b>=12 && Math.max(a/(double)(length/2),b/(double)(length-length/2))>=.40
				&& Math.min(a/(double)(length/2),b/(double)(length-length/2))<=.10;
	}

	static String geometry(SiteScore a,int length1,SiteScore b,int length2){
		final Shape one=shape(a,length1),two=shape(b,length2);
		if(!one.state.equals("OK")){return one.state;}
		if(!two.state.equals("OK")){return two.state;}
		if(one.chrom!=two.chrom || one.scaffold!=two.scaffold){return "UNEVALUABLE_REFERENCE";}
		return overlaps(one.deletions,two.aligned)+overlaps(two.deletions,one.aligned)>0 ? "CONFLICT" : "NO_CONFLICT";
	}

	private static Shape shape(SiteScore site,int length){
		if(site==null){return new Shape("UNEVALUABLE_UNMAPPED");}
		if(site.match==null){return new Shape("UNEVALUABLE_NATIVE_MATCH");}
		if(Data.scaffoldLocs==null || Data.scaffoldLengths==null){return new Shape("UNEVALUABLE_NATIVE_METADATA");}
		final int idx=Data.scaffoldIndex(site.chrom,site.start);
		final int left=Data.scaffoldLocs[site.chrom][idx],size=Data.scaffoldLengths[site.chrom][idx];
		if(site.start<left || site.stop<site.start || (long)site.stop>=(long)left+size){return new Shape("UNEVALUABLE_NATIVE_BOUNDS");}
		final Shape out=new Shape("OK");out.chrom=site.chrom;out.scaffold=idx;
		int r=site.start,q=0;
		for(byte op:site.match){if("CXY".indexOf((char)op)>=0){return new Shape("UNEVALUABLE_NATIVE_CLIP");}}
		for(int i=0;i<site.match.length;){
			final byte op=site.match[i];
			assert("msSVNBID".indexOf((char)op)>=0) : "Unsupported native geometry operation: "+op;
			int n=1;while(i+n<site.match.length && site.match[i+n]==op){n++;}i+=n;
			if(op=='D'){
				if(n>SamLine.INTRON_LIMIT){return new Shape("UNEVALUABLE_NATIVE_INTRON");}
				if(n>=101){out.deletions.add(new int[]{r,Math.addExact(r,n)});}
				r=Math.addExact(r,n);
			}else if(op=='I'){q=Math.addExact(q,n);}
			else{out.aligned.add(new int[]{r,Math.addExact(r,n)});r=Math.addExact(r,n);q=Math.addExact(q,n);}
		}
		assert(q==length && r==(long)site.stop+1) : "Native geometry consumption must match query length and candidate limits";
		return out;
	}

	private static int overlaps(ArrayList<int[]> gaps,ArrayList<int[]> blocks){
		int n=0;
		for(int[] gap:gaps){for(int[] block:blocks){if(gap[0]<block[1] && block[0]<gap[1]){n++;break;}}}
		return n;
	}

	public static final class Observation {
		Observation(int f,String g){flags=f;geometry=g;}
		public final int flags;
		public final String geometry;
		public boolean route(){return (flags&(NO_SITE|HALF_ERRORS))!=0;}
	}
	private static final class Shape {
		Shape(String s){state=s;}
		final String state;
		int chrom,scaffold;
		final ArrayList<int[]> aligned=new ArrayList<int[]>(),deletions=new ArrayList<int[]>();
	}
}
