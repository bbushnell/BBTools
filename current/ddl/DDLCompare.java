package ddl;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;

import cardinality.DynamicDemiLog;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Timer;
import simd.Vector;
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

		String file1=null, file2=null, refFile=null, queryFile=null;
		int k=31, buckets=2048, maxRecords=Integer.MAX_VALUE, minHits=3, threads=1;
		boolean collisionTest=false, useIndex=false;
		for(int i=0; i<args.length; i++){
			String[] split=args[i].split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("exponent") || a.equals("ebits")){DynamicDemiLog.setExponent(Integer.parseInt(b));}
			else if(a.equals("buckets")){buckets=Integer.parseInt(b);}
			else if(a.equals("ref")){refFile=b;}
			else if(a.equals("queryfile") || a.equals("qf")){queryFile=b;}
			else if(a.equals("records") || a.equals("maxrecords")){maxRecords=Integer.parseInt(b);}
			else if(a.equals("minhits")){minHits=Integer.parseInt(b);}
			else if(a.equals("collisiontest")){collisionTest=true;}
			else if(a.equals("index")){useIndex=b==null || b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true");}
			else if(a.equals("t") || a.equals("threads")){threads=Integer.parseInt(b);}
			else if(a.equals("in") || a.equals("in1") || a.equals("query")){file1=b;}
			else if(a.equals("in2")){file2=b;}
			else if(file1==null){file1=args[i];}
			else if(file2==null){file2=args[i];}
		}

		if(collisionTest && refFile!=null){
			collisionTest(refFile, k);
			return;
		}
		if(queryFile!=null && refFile!=null){
			compareMultiQuery(queryFile, refFile, k, maxRecords, minHits, useIndex, threads);
			return;
		}
		if(refFile!=null){
			compareToRefs(file1, refFile, k, buckets, maxRecords, minHits, useIndex, threads);
			return;
		}

		Timer t=new Timer();

		DynamicDemiLog ddlA, ddlB;
		long basesA, basesB;
		if(threads>1){
			ddlA=(DynamicDemiLog)MultithreadedSketchLoader.loadTrackerFromSequence(
				file1, "DDL", buckets, k, 12345L, 0f, false, threads);
			ddlB=(DynamicDemiLog)MultithreadedSketchLoader.loadTrackerFromSequence(
				file2, "DDL", buckets, k, 12345L, 0f, false, threads);
			basesA=-1; basesB=-1;
		}else{
			ddlA=DynamicDemiLog.create(buckets, k, 12345L, 0f, true);
			ddlB=DynamicDemiLog.create(buckets, k, 12345L, 0f, true);
			basesA=hashFile(file1, ddlA, k);
			basesB=hashFile(file2, ddlB, k);
		}

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

		int[] ll6cmp=compareExponentOnly(ddlA.maxArray(), ddlB.maxArray());
		int ll6lower=ll6cmp[0], ll6equal=ll6cmp[1], ll6higher=ll6cmp[2];
		float ll6wkid=DynamicDemiLog.wkid(ll6lower, ll6equal, ll6higher);
		float ll6ani=DynamicDemiLog.ani(ll6lower, ll6equal, ll6higher, k);

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
		System.out.println("LL6_Equal:\t"+ll6equal);
		System.out.println("LL6_WKID:\t"+String.format("%.6f", ll6wkid));
		System.out.println("LL6_ANI:\t"+String.format("%.6f", ll6ani));
		System.out.println("Time:\t"+t);

		System.err.println(String.format("ANI: %.2f%%  WKID: %.4f  Comp(1->2): %.4f  Comp(2->1): %.4f  (%d/%d buckets match)",
			ani*100, c, compAB, compBA, equal, lower+equal+higher));
	}

	/** Compares a query FASTA against pre-built reference DDLs. */
	private static void compareToRefs(String queryPath, String refPath, int k, int buckets,
			int maxRecords, int minHits, boolean useIndex, int threads){
		Timer t=new Timer();
		System.err.println("Mode: "+(useIndex ? "indexed" : "pairwise")+"  threads="+threads);

		System.err.println("Loading references from "+refPath+"...");
		long t0=System.nanoTime();
		ArrayList<DDLRecord> refs=DDLLoaderMT.loadFile(refPath, k, threads);
		long t1=System.nanoTime();
		System.err.println("Loaded "+refs.size()+" reference DDLs in "+String.format("%.3f", (t1-t0)*1e-9)+" seconds.");

		DDLIndex index=null;
		if(useIndex){
			long ti0=System.nanoTime();
			index=new DDLIndex();
			index.addAll(refs, threads);
			long ti1=System.nanoTime();
			System.err.println("Built inverted index in "+String.format("%.3f", (ti1-ti0)*1e-9)+" seconds.");
		}

		DynamicDemiLog query;
		long bases;
		if(threads>1){
			query=(DynamicDemiLog)MultithreadedSketchLoader.loadTrackerFromSequence(
				queryPath, "DDL", buckets, k, 12345L, 0f, false, threads);
			bases=-1;
		}else{
			query=DynamicDemiLog.create(buckets, k, 12345L, 0f, true);
			bases=hashFile(queryPath, query, k);
		}
		long t2=System.nanoTime();
		long card=query.cardinality();
		System.err.println("Query: "+queryPath+(bases>=0 ? "  bases="+bases : "")+"  cardinality="+card+"  sketch time: "+String.format("%.3f", (t2-t1)*1e-9)+"s  threads="+threads);

		final int n=refs.size();

		if(useIndex){
			long t3=System.nanoTime();
			int[][] topHits=index.topHits(query, Math.min(maxRecords*2, n));
			long t4=System.nanoTime();
			System.err.println("Index query returned "+topHits.length+" hits in "+String.format("%.6f", (t4-t3)*1e-9)+" seconds.");

			System.out.println("ANI\tWKID\tComplt\tMatches\tBases\tTID\tName");
			int shown=0;
			for(int[] hit : topHits){
				if(shown>=maxRecords){break;}
				int cladeID=hit[0], matchCount=hit[1];
				if(matchCount<minHits){continue;}
				DDLRecord ref=refs.get(cladeID);
				float ani=index.ani(matchCount, query.filledBuckets(), cladeID, k);
				System.out.println(String.format("%.4f\t-\t-\t%d\t%d\t%d\t%s",
					ani, matchCount, ref.bases, ref.taxID, ref.name));
				shown++;
			}
		}else{
			long t3=System.nanoTime();
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
			long t4=System.nanoTime();
			System.err.println("Compared against "+n+" references in "+String.format("%.3f", (t4-t3)*1e-9)+" seconds.");

			for(int i=1; i<n; i++){
				for(int j=i; j>0 && anis[indices[j]]>anis[indices[j-1]]; j--){
					int tmp=indices[j]; indices[j]=indices[j-1]; indices[j-1]=tmp;
				}
			}

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
		}

		t.stop();
		System.err.println("Total time: \t"+t);
	}

	/** Compares multiple pre-built DDL queries against pre-built DDL references. */
	private static void compareMultiQuery(String queryPath, String refPath, int k,
			int maxRecords, int minHits, boolean useIndex, int threads){
		Timer t=new Timer();
		System.err.println("Mode: multi-query "+(useIndex ? "indexed" : "pairwise")+"  threads="+threads);

		System.err.println("Loading references from "+refPath+"...");
		long t0=System.nanoTime();
		ArrayList<DDLRecord> refs=DDLLoaderMT.loadFile(refPath, k, threads);
		long t1=System.nanoTime();
		System.err.println("Loaded "+refs.size()+" reference DDLs in "+String.format("%.3f", (t1-t0)*1e-9)+" seconds.");

		System.err.println("Loading queries from "+queryPath+"...");
		long tq0=System.nanoTime();
		ArrayList<DDLRecord> queries=DDLLoaderMT.loadFile(queryPath, k, threads);
		long tq1=System.nanoTime();
		System.err.println("Loaded "+queries.size()+" query DDLs in "+String.format("%.3f", (tq1-tq0)*1e-9)+" seconds.");

		final int nRefs=refs.size();
		final int nQueries=queries.size();
		final int[] bestMatch=new int[nQueries];
		final float[] bestAni=new float[nQueries];
		final int[] bestHits=new int[nQueries];

		final int buckets=refs.get(0).ddl.maxArray().length;

		DDLIndex index=null;
		if(useIndex){
			long ti0=System.nanoTime();
			index=new DDLIndex();
			index.addAll(refs, threads);
			long ti1=System.nanoTime();
			System.err.println("Built inverted index in "+String.format("%.3f", (ti1-ti0)*1e-9)+" seconds.");
		}

		final AtomicLong nextQuery=new AtomicLong(0);
		final DDLIndex idx=index;
		long tc0=System.nanoTime();

		CompareThread[] workers=new CompareThread[threads];
		for(int wi=0; wi<threads; wi++){
			workers[wi]=new CompareThread(queries, refs, nextQuery,
				bestMatch, bestAni, bestHits, nQueries, nRefs, buckets, k, useIndex, idx);
			workers[wi].start();
		}
		for(CompareThread w : workers){try{w.join();}catch(InterruptedException e){}}

		long tc1=System.nanoTime();
		long totalComparisons=(long)nQueries*nRefs;
		System.err.println("Compared "+nQueries+" queries × "+nRefs+" refs ("+totalComparisons+" comparisons) in "+String.format("%.3f", (tc1-tc0)*1e-9)+" seconds.");

		System.out.println("ANI\tWKID\tMatches\tBases\tTID\tQueryName\tRefName");
		int shown=0;
		for(int qi=0; qi<nQueries && shown<maxRecords; qi++){
			if(bestHits[qi]<minHits) continue;
			DDLRecord ref=refs.get(bestMatch[qi]);
			DDLRecord q=queries.get(qi);
			System.out.println(String.format("%.4f\t-\t%d\t%d\t%d\t%s\t%s",
				bestAni[qi], bestHits[qi], ref.bases, ref.taxID, q.name, ref.name));
			shown++;
		}

		t.stop();
		System.err.println("Total time: \t"+t);
	}

	static class CompareThread extends Thread {
		final ArrayList<DDLRecord> queries, refs;
		final AtomicLong nextQuery;
		final int[] bestMatch, bestHits;
		final float[] bestAni;
		final int nQueries, nRefs, buckets, k;
		final boolean useIndex;
		final DDLIndex idx;

		CompareThread(ArrayList<DDLRecord> queries, ArrayList<DDLRecord> refs,
				AtomicLong nextQuery, int[] bestMatch, float[] bestAni, int[] bestHits,
				int nQueries, int nRefs, int buckets, int k, boolean useIndex, DDLIndex idx){
			this.queries=queries; this.refs=refs;
			this.nextQuery=nextQuery;
			this.bestMatch=bestMatch; this.bestAni=bestAni; this.bestHits=bestHits;
			this.nQueries=nQueries; this.nRefs=nRefs; this.buckets=buckets;
			this.k=k; this.useIndex=useIndex; this.idx=idx;
		}

		@Override
		public void run(){
			while(true){
				int qi=(int)nextQuery.getAndIncrement();
				if(qi>=nQueries) break;
				if(useIndex){
					DynamicDemiLog q=queries.get(qi).ddl;
					int[][] hits=idx.topHits(q, 1);
					if(hits.length>0){
						bestMatch[qi]=hits[0][0];
						bestHits[qi]=hits[0][1];
						bestAni[qi]=idx.ani(hits[0][1], q.filledBuckets(), hits[0][0], k);
					}
				}else{
					char[] qa=queries.get(qi).ddl.maxArray();
					float topAni=0;
					int topIdx=0, topHits=0;
					for(int ri=0; ri<nRefs; ri++){
						int[] cmp=Vector.compareDDL(qa, refs.get(ri).ddl.maxArray());
						float ani=DynamicDemiLog.ani(cmp[0], cmp[1], cmp[2], k);
						if(ani>topAni){topAni=ani; topIdx=ri; topHits=cmp[1];}
					}
					bestMatch[qi]=topIdx;
					bestAni[qi]=topAni;
					bestHits[qi]=topHits;
				}
			}
		}
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

	/** Compares two DDL register arrays using only the exponent bits,
	 *  simulating LL-equivalent comparison from DDL sketches. */
	private static int[] compareExponentOnly(char[] a, char[] b){
		final int mbits=16-DynamicDemiLog.exponentBits();
		int lower=0, equal=0, higher=0;
		for(int i=0; i<a.length; i++){
			int ea=a[i]>>mbits, eb=b[i]>>mbits;
			if(ea==0 && eb==0){/* both empty */}
			else if(ea<eb){lower++;}
			else if(ea==eb){equal++;}
			else{higher++;}
		}
		return new int[]{lower, equal, higher};
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
