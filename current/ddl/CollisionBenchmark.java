package ddl;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import stream.Read;

/**
 * Measures per-register false-match rates between random genomes.
 * Generates N random genomes internally, builds DDL and LL6 sketches,
 * and performs all-pairs comparison in a single JVM invocation.
 *
 * Usage: CollisionBenchmark [n=100] [len=4000000] [k=31] [buckets=2048]
 *
 * @author Noire
 * @date May 12, 2026
 */
public class CollisionBenchmark {

	public static void main(String[] args){
		int n=100, len=4000000, k=31, buckets=2048;
		long seed=42L;

		for(String arg : args){
			String[] ab=arg.split("=");
			String a=ab[0].toLowerCase(), b=ab.length>1 ? ab[1] : null;
			if(b==null) continue;
			if(a.equals("n")){n=Integer.parseInt(b);}
			else if(a.equals("len") || a.equals("length")){len=Integer.parseInt(b);}
			else if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("buckets")){buckets=Integer.parseInt(b);}
			else if(a.equals("seed")){seed=Long.parseLong(b);}
		}

		System.err.println("CollisionBenchmark: n="+n+" len="+len+" k="+k+" buckets="+buckets);

		// Build sketches from random genomes
		int[][] ddlVals=new int[n][];
		int[][] ll6Vals=new int[n][];
		long rng=seed;
		for(int i=0; i<n; i++){
			rng=mix(rng+i*99991L);
			byte[] bases=randomGenome(len, rng);
			Read r=new Read(bases, null, "genome_"+i, 0);

			CardinalityTracker ddl=CardinalityTracker.makeTracker("ddl", buckets, k, 12345L, 0f);
			ddl.hash(r);
			ddlVals[i]=ddl.bucketValues();

			CardinalityTracker ll6=CardinalityTracker.makeTracker("ll6", buckets, k, 12345L, 0f);
			ll6.hash(r);
			ll6Vals[i]=ll6.bucketValues();

			if((i+1)%10==0){System.err.println("  Built "+i+1+" sketches...");}
		}
		System.err.println("Built "+n+" sketches. Comparing all pairs...");

		// All-pairs comparison
		long pairs=0, ddlTotal=0, ll6Total=0;
		int ddlMax=0, ll6Max=0;

		System.out.println("IdxA\tIdxB\tDDL_equal\tLL6_equal");
		for(int a=0; a<n; a++){
			for(int b=a+1; b<n; b++){
				int ddlEq=countEqual(ddlVals[a], ddlVals[b]);
				int ll6Eq=countEqual(ll6Vals[a], ll6Vals[b]);
				pairs++;
				ddlTotal+=ddlEq;
				ll6Total+=ll6Eq;
				if(ddlEq>ddlMax){ddlMax=ddlEq;}
				if(ll6Eq>ll6Max){ll6Max=ll6Eq;}
				if(ddlEq>0 || ll6Eq>0){
					System.out.println(a+"\t"+b+"\t"+ddlEq+"\t"+ll6Eq);
				}
			}
		}

		System.err.println("Pairs: "+pairs);
		System.err.println("DDL avg false matches/pair: "+String.format("%.4f", (double)ddlTotal/pairs));
		System.err.println("DDL max false matches: "+ddlMax);
		System.err.println("DDL per-register rate: "+String.format("%.6f", (double)ddlTotal/(pairs*buckets)));
		System.err.println("LL6 avg false matches/pair: "+String.format("%.4f", (double)ll6Total/pairs));
		System.err.println("LL6 max false matches: "+ll6Max);
		System.err.println("LL6 per-register rate: "+String.format("%.6f", (double)ll6Total/(pairs*buckets)));
	}

	private static int countEqual(int[] a, int[] b){
		int eq=0;
		for(int i=0; i<a.length; i++){
			if(a[i]!=0 && a[i]==b[i]){eq++;}
		}
		return eq;
	}

	private static byte[] randomGenome(int len, long seed){
		byte[] bases=new byte[len];
		long state=seed;
		for(int i=0; i<len; i++){
			state=mix(state+i);
			bases[i]="ACGT".getBytes()[(int)(((state>>>32)&0x7FFFFFFFL)%4)];
		}
		return bases;
	}

	private static long mix(long x){
		x=((x>>>16)^x)*0x45d9f3bL;
		x=((x>>>16)^x)*0x45d9f3bL;
		x=(x>>>16)^x;
		return x;
	}
}
