package ddl;

import java.util.ArrayList;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;

/**
 * Compares two FASTA files using both DDL and LL6 sketches.
 * Outputs WKID and ANI estimates for each sketch type.
 * Used for benchmarking DDL vs LL6 on real nucleotide sequences
 * at known mutation rates.
 *
 * Usage: ANIBenchmark <ref.fa> <query.fa> [k=31] [buckets=2048]
 *
 * @author Noire
 * @date May 12, 2026
 */
public class ANIBenchmark {

	public static void main(String[] args){
		if(args.length<2){
			System.err.println("Usage: ANIBenchmark <ref.fa> <query.fa> [k=31] [buckets=2048]");
			System.exit(1);
		}

		String file1=null, file2=null;
		int k=31, buckets=2048;
		for(String arg : args){
			String[] ab=arg.split("=");
			String a=ab[0].toLowerCase(), b=ab.length>1 ? ab[1] : null;
			if(a.equals("k") && b!=null){k=Integer.parseInt(b);}
			else if(a.equals("buckets") && b!=null){buckets=Integer.parseInt(b);}
			else if(file1==null){file1=arg;}
			else if(file2==null){file2=arg;}
		}

		long seed=12345L;

		// Build DDL sketches
		CardinalityTracker ddlRef=CardinalityTracker.makeTracker("ddl", buckets, k, seed, 0f);
		CardinalityTracker ddlQry=CardinalityTracker.makeTracker("ddl", buckets, k, seed, 0f);
		hashFile(file1, ddlRef);
		hashFile(file2, ddlQry);

		// Build LL6 sketches
		CardinalityTracker ll6Ref=CardinalityTracker.makeTracker("ll6", buckets, k, seed, 0f);
		CardinalityTracker ll6Qry=CardinalityTracker.makeTracker("ll6", buckets, k, seed, 0f);
		hashFile(file1, ll6Ref);
		hashFile(file2, ll6Qry);

		// Compare DDL
		int[] ddlValsRef=ddlRef.bucketValues();
		int[] ddlValsQry=ddlQry.bucketValues();
		int[] ddlCmp=compareBuckets(ddlValsRef, ddlValsQry);
		float ddlWkid=DynamicDemiLog.wkid(ddlCmp[0], ddlCmp[1], ddlCmp[2]);
		float ddlAni=DynamicDemiLog.ani(ddlCmp[0], ddlCmp[1], ddlCmp[2], k);

		// Compare LL6
		int[] ll6ValsRef=ll6Ref.bucketValues();
		int[] ll6ValsQry=ll6Qry.bucketValues();
		int[] ll6Cmp=compareBuckets(ll6ValsRef, ll6ValsQry);
		float ll6Wkid=DynamicDemiLog.wkid(ll6Cmp[0], ll6Cmp[1], ll6Cmp[2]);
		float ll6Ani=DynamicDemiLog.ani(ll6Cmp[0], ll6Cmp[1], ll6Cmp[2], k);

		// Output: tab-separated for easy parsing
		// DDL_WKID DDL_ANI DDL_equal LL6_WKID LL6_ANI LL6_equal
		System.out.println(String.format("%.6f\t%.6f\t%d\t%.6f\t%.6f\t%d",
			ddlWkid, ddlAni, ddlCmp[1],
			ll6Wkid, ll6Ani, ll6Cmp[1]));
	}

	private static int[] compareBuckets(int[] a, int[] b){
		int lower=0, equal=0, higher=0, bothEmpty=0;
		for(int i=0; i<a.length; i++){
			if(a[i]==0 && b[i]==0){bothEmpty++; continue;}
			if(a[i]<b[i]){lower++;}
			else if(a[i]>b[i]){higher++;}
			else{equal++;}
		}
		return new int[]{lower, equal, higher, bothEmpty};
	}

	private static void hashFile(String path, CardinalityTracker ct){
		FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		Streamer cris=StreamerFactory.getReadInputStream(-1, false, ff, null, -1);
		cris.start();
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		while(ln!=null && reads!=null && reads.size()>0){
			for(Read r : reads){
				ct.hash(r);
			}
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln);
		ReadWrite.closeStreams(cris);
	}
}
