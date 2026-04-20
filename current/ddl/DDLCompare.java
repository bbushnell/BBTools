package ddl;

import java.util.ArrayList;

import cardinality.DynamicDemiLog;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Timer;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;

/**
 * Pairwise genome comparison using DynamicDemiLog bucket matching.
 * Creates a DDL for each input genome, compares them, and reports
 * WKID, ANI, cardinality, and bucket-level statistics.
 *
 * @author Brian Bushnell, Ady
 * @date April 17, 2026
 */
public class DDLCompare {

	public static void main(String[] args){
		if(args.length<1){
			System.err.println("Usage: DDLCompare <genome1.fa> <genome2.fa> [k=31] [buckets=2048]");
			System.err.println("   or: DDLCompare <query.fa> ref=<ddls.tsv> [k=31]");
			System.exit(1);
		}

		String file1=null, file2=null, refFile=null;
		int k=31, buckets=2048, maxRecords=Integer.MAX_VALUE, minHits=3;
		boolean collisionTest=false;
		for(int i=0; i<args.length; i++){
			String[] split=args[i].split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("buckets")){buckets=Integer.parseInt(b);}
			else if(a.equals("ref")){refFile=b;}
			else if(a.equals("records") || a.equals("maxrecords")){maxRecords=Integer.parseInt(b);}
			else if(a.equals("minhits")){minHits=Integer.parseInt(b);}
			else if(a.equals("collisiontest")){collisionTest=true;}
			else if(a.equals("in") || a.equals("in1") || a.equals("query")){file1=b;}
			else if(a.equals("in2")){file2=b;}
			else if(file1==null){file1=args[i];}
			else if(file2==null){file2=args[i];}
		}

		if(collisionTest && refFile!=null){
			collisionTest(refFile, k);
			return;
		}
		if(refFile!=null){
			compareToRefs(file1, refFile, k, buckets, maxRecords, minHits);
			return;
		}

		Timer t=new Timer();

		DynamicDemiLog ddlA=DynamicDemiLog.create(buckets, k, 12345L, 0f);
		DynamicDemiLog ddlB=DynamicDemiLog.create(buckets, k, 12345L, 0f);

		long basesA=hashFile(file1, ddlA, k);
		long basesB=hashFile(file2, ddlB, k);

		long cardA=ddlA.cardinality();
		long cardB=ddlB.cardinality();

		int[] cmp=ddlA.compareToDetailed(ddlB);
		int lower=cmp[0], equal=cmp[1], higher=cmp[2], bothEmpty=cmp[3];

		float c=DynamicDemiLog.wkid(lower, equal, higher);
		float cAB=DynamicDemiLog.containmentAB(lower, equal, higher);
		float cBA=DynamicDemiLog.containmentBA(lower, equal, higher);
		float ani=DynamicDemiLog.ani(lower, equal, higher, k);
		float compAB=DynamicDemiLog.completeness(lower, equal, higher);
		float compBA=DynamicDemiLog.completenessBA(lower, equal, higher);

		t.stop();

		System.out.println("File1:\t"+file1);
		System.out.println("File2:\t"+file2);
		System.out.println("Bases1:\t"+basesA);
		System.out.println("Bases2:\t"+basesB);
		System.out.println("Cardinality1:\t"+cardA);
		System.out.println("Cardinality2:\t"+cardB);
		System.out.println("K:\t"+k);
		System.out.println("Buckets:\t"+buckets);
		System.out.println("Lower:\t"+lower);
		System.out.println("Equal:\t"+equal);
		System.out.println("Higher:\t"+higher);
		System.out.println("BothEmpty:\t"+bothEmpty);
		System.out.println("WKID:\t"+String.format("%.6f", c));
		System.out.println("WKID(1in2):\t"+String.format("%.6f", cAB));
		System.out.println("WKID(2in1):\t"+String.format("%.6f", cBA));
		System.out.println("ANI:\t"+String.format("%.6f", ani));
		System.out.println("Completeness(1->2):\t"+String.format("%.6f", compAB));
		System.out.println("Completeness(2->1):\t"+String.format("%.6f", compBA));
		System.out.println("Time:\t"+t);

		System.err.println(String.format("ANI: %.2f%%  WKID: %.4f  Comp(1->2): %.4f  Comp(2->1): %.4f  (%d/%d buckets match)",
			ani*100, c, compAB, compBA, equal, lower+equal+higher));
	}

	/** Compares a query FASTA against pre-built reference DDLs. */
	private static void compareToRefs(String queryPath, String refPath, int k, int buckets, int maxRecords, int minHits){
		Timer t=new Timer();

		System.err.println("Loading references from "+refPath+"...");
		ArrayList<DDLRecord> refs=DDLLoader.loadFile(refPath, k);
		System.err.println("Loaded "+refs.size()+" reference DDLs.");

		DynamicDemiLog query=DynamicDemiLog.create(buckets, k, 12345L, 0f);
		long bases=hashFile(queryPath, query, k);
		long card=query.cardinality();
		System.err.println("Query: "+queryPath+"  bases="+bases+"  cardinality="+card);

		//Score all references
		final int n=refs.size();
		final float[] anis=new float[n];
		final int[] matches=new int[n];
		final int[] indices=new int[n];
		for(int i=0; i<n; i++){
			DDLRecord ref=refs.get(i);
			int[] cmp=query.compareToDetailed(ref.ddl);
			matches[i]=cmp[1];
			anis[i]=DynamicDemiLog.ani(cmp[0], cmp[1], cmp[2], k);
			indices[i]=i;
		}

		//Sort by ANI descending (insertion sort is fine for small output)
		for(int i=1; i<n; i++){
			for(int j=i; j>0 && anis[indices[j]]>anis[indices[j-1]]; j--){
				int tmp=indices[j]; indices[j]=indices[j-1]; indices[j-1]=tmp;
			}
		}

		//Print top records (skip those below minHits)
		System.out.println("ANI\tWKID\tComplt\tMatches\tBases\tTID\tName");
		int shown=0;
		final int toShow=Math.min(maxRecords, n);
		for(int r=0; r<n && shown<toShow; r++){
			final int i=indices[r];
			DDLRecord ref=refs.get(i);
			int[] cmp=query.compareToDetailed(ref.ddl);
			int lower=cmp[0], equal=cmp[1], higher=cmp[2];
			if(equal<minHits){continue;}
			float ani=DynamicDemiLog.ani(lower, equal, higher, k);
			float c=DynamicDemiLog.wkid(lower, equal, higher);
			float comp=DynamicDemiLog.completeness(lower, equal, higher);
			System.out.println(String.format("%.4f\t%.4f\t%.4f\t%d\t%d\t%d\t%s",
				ani, c, comp, equal, ref.bases, ref.taxID, ref.name));
			shown++;
		}

		t.stop();
		System.err.println("Time: \t"+t);
	}

	/** Measures collision rates across all pairs in a DDL file. */
	private static void collisionTest(String refPath, int k){
		System.err.println("Loading references from "+refPath+"...");
		ArrayList<DDLRecord> refs=DDLLoader.loadFile(refPath, k);
		final int n=refs.size();
		System.err.println("Loaded "+n+" DDLs. Testing "+n*(n-1)/2+" pairs...");

		long totalPairs=0, totalMatches=0;
		int maxMatches=0;
		System.out.println("RefA\tRefB\tMatches\tNameA\tNameB");
		for(int a=0; a<n; a++){
			for(int b=a+1; b<n; b++){
				int[] cmp=refs.get(a).ddl.compareToDetailed(refs.get(b).ddl);
				int equal=cmp[1];
				totalPairs++;
				totalMatches+=equal;
				if(equal>maxMatches){maxMatches=equal;}
				if(equal>0){
					System.out.println(a+"\t"+b+"\t"+equal+"\t"+refs.get(a).name+"\t"+refs.get(b).name);
				}
			}
		}
		System.err.println("Pairs: "+totalPairs);
		System.err.println("Total matches: "+totalMatches);
		System.err.println("Avg matches/pair: "+String.format("%.4f", (double)totalMatches/totalPairs));
		System.err.println("Max matches in any pair: "+maxMatches);
		System.err.println("Avg collision rate/bucket: "+String.format("%.6f", (double)totalMatches/(totalPairs*2048)));
	}

	private static long hashFile(String path, DynamicDemiLog ddl, int k){
		FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		Streamer cris=StreamerFactory.getReadInputStream(-1, false, ff, null, -1);
		cris.start();

		long bases=0;
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		while(ln!=null && reads!=null && reads.size()>0){
			for(Read r : reads){
				ddl.hash(r);
				bases+=r.length();
				if(r.mate!=null){bases+=r.mate.length();}
			}
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln);
		ReadWrite.closeStreams(cris);
		return bases;
	}
}
