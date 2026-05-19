package ddl;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;

import cardinality.DynamicDemiLog;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import idaligner.QuantumAligner;
import shared.Timer;
import simd.Vector;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Pairwise genome comparison using DynamicDemiLog bucket matching.
 * Creates a DDL for each input genome, compares them, and reports
 * WKID, ANI, cardinality, and bucket-level statistics.
 *
 * @author Brian Bushnell, Ady, Noire
 * @date April 17, 2026
 */
public class DDLCompare {

	public static void main(String[] args){
		if(args.length<1){
			System.err.println("Usage: DDLCompare <genome1.fa> <genome2.fa> [k=31] [buckets=2048]");
			System.err.println("   or: DDLCompare <query.fa> ref=<ddls.tsv> [records=20]");
			System.err.println("   or: DDLCompare qf=queries.tsv ref=ddls.tsv t=32 [format=json]");
			System.exit(1);
		}

		String file1=null, file2=null, refFile=null, queryFile=null;
		int k=31, buckets=2048, maxRecords=20, minHits=5, threads=1;
		boolean collisionTest=false, useIndex=false, useSSU=false;
		DDLFormatter formatter=new DDLFormatter();

		for(int i=0; i<args.length; i++){
			String[] split=args[i].split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(formatter.parse(args[i], a, b)){/* handled by formatter */}
			else if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("exponent") || a.equals("ebits")){DynamicDemiLog.setExponent(Integer.parseInt(b));}
			else if(a.equals("buckets")){buckets=Integer.parseInt(b);}
			else if(a.equals("ref")){refFile=b;}
			else if(a.equals("queryfile") || a.equals("qf")){queryFile=b;}
			else if(a.equals("records") || a.equals("maxrecords")){maxRecords=Integer.parseInt(b);}
			else if(a.equals("minhits")){minHits=Integer.parseInt(b);}
			else if(a.equals("collisiontest")){collisionTest=true;}
			else if(a.equals("index")){useIndex=b==null || b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true");}
			else if(a.equals("ssu")){useSSU=b==null || b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true");}
			else if(a.equals("t") || a.equals("threads")){threads=Integer.parseInt(b);}
			else if(a.equals("in") || a.equals("in1") || a.equals("query")){file1=b;}
			else if(a.equals("in2")){file2=b;}
			else if(file1==null){file1=args[i];}
			else if(file2==null){file2=args[i];}
		}

		if(queryFile!=null){
			formatter.printQueryName=true;
		}
		if(useSSU){
			formatter.printSSU=true;
		}

		if(collisionTest && refFile!=null){
			collisionTest(refFile, k);
			return;
		}
		if(queryFile!=null && refFile!=null){
			compareMultiQuery(queryFile, refFile, k, maxRecords, minHits, useIndex, useSSU, threads, formatter);
			return;
		}
		if(refFile!=null){
			compareToRefs(file1, refFile, k, buckets, maxRecords, minHits, useIndex, useSSU, threads, formatter);
			return;
		}

		comparePairwise(file1, file2, k, buckets, threads);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Pairwise Mode          ----------------*/
	/*--------------------------------------------------------------*/

	/** Detailed pairwise comparison of two genome files. */
	private static void comparePairwise(String file1, String file2, int k, int buckets, int threads){
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

	/*--------------------------------------------------------------*/
	/*----------------     Single-Query Mode        ----------------*/
	/*--------------------------------------------------------------*/

	/** Compares a query FASTA against pre-built reference DDLs. */
	private static void compareToRefs(String queryPath, String refPath, int k, int buckets,
			int maxRecords, int minHits, boolean useIndex, boolean useSSU, int threads, DDLFormatter formatter){
		Timer t=new Timer();
		System.err.println("Mode: "+(useIndex ? "indexed" : "pairwise")+"  threads="+threads+(useSSU ? "  ssu=t" : ""));

		System.err.println("Loading references from "+refPath+"...");
		long t0=System.nanoTime();
		ArrayList<DDLRecord> refs=DDLLoaderMT.loadFile(refPath, k, threads);
		long t1=System.nanoTime();
		System.err.println("Loaded "+refs.size()+" reference DDLs in "+String.format("%.3f", (t1-t0)*1e-9)+" seconds.");

		if(useSSU){DDLSSULoader.loadAndAttachDefaults(refs);}

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

		DDLRecord qRec=new DDLRecord(query, -1, -1, queryPath);
		qRec.bases=bases>=0 ? bases : 0;
		qRec.cardinality=card;

		final int n=refs.size();
		DDLComparisonHeap heap=new DDLComparisonHeap(maxRecords);
		DDLComparison working=new DDLComparison();

		long t3=System.nanoTime();
		if(useIndex){
			int[] counts=index.query(query);
			for(int ri=0; ri<n; ri++){
				if(counts[ri]<minHits){continue;}
				working.compare(qRec, refs.get(ri), k);
				heap.offer(working);
			}
		}else{
			for(int ri=0; ri<n; ri++){
				working.compare(qRec, refs.get(ri), k);
				heap.offer(working);
			}
		}
		long t4=System.nanoTime();
		System.err.println("Compared against "+n+" references in "+String.format("%.3f", (t4-t3)*1e-9)+" seconds.");

		ArrayList<DDLComparison> results=heap.toList();
		if(useSSU){alignSSU(results);}

		printResults(results, minHits, formatter);

		t.stop();
		System.err.println("Total time: \t"+t);
	}

	/*--------------------------------------------------------------*/
	/*----------------     Multi-Query Mode         ----------------*/
	/*--------------------------------------------------------------*/

	/** Compares multiple pre-built DDL queries against pre-built DDL references. */
	private static void compareMultiQuery(String queryPath, String refPath, int k,
			int maxRecords, int minHits, boolean useIndex, boolean useSSU, int threads, DDLFormatter formatter){
		Timer t=new Timer();
		System.err.println("Mode: multi-query "+(useIndex ? "indexed" : "pairwise")+"  threads="+threads+(useSSU ? "  ssu=t" : ""));

		System.err.println("Loading references from "+refPath+"...");
		long t0=System.nanoTime();
		ArrayList<DDLRecord> refs=DDLLoaderMT.loadFile(refPath, k, threads);
		long t1=System.nanoTime();
		System.err.println("Loaded "+refs.size()+" reference DDLs in "+String.format("%.3f", (t1-t0)*1e-9)+" seconds.");

		if(useSSU){DDLSSULoader.loadAndAttachDefaults(refs);}

		System.err.println("Loading queries from "+queryPath+"...");
		long tq0=System.nanoTime();
		ArrayList<DDLRecord> queries=DDLLoaderMT.loadFile(queryPath, k, threads);
		long tq1=System.nanoTime();
		System.err.println("Loaded "+queries.size()+" query DDLs in "+String.format("%.3f", (tq1-tq0)*1e-9)+" seconds.");

		if(useSSU){DDLSSULoader.loadAndAttachDefaults(queries);}

		final int nRefs=refs.size();
		final int nQueries=queries.size();

		DDLIndex index=null;
		if(useIndex){
			long ti0=System.nanoTime();
			index=new DDLIndex();
			index.addAll(refs, threads);
			long ti1=System.nanoTime();
			System.err.println("Built inverted index in "+String.format("%.3f", (ti1-ti0)*1e-9)+" seconds.");
		}

		final AtomicLong nextQuery=new AtomicLong(0);
		final AtomicLong totalComparisonsPerformed=new AtomicLong(0);
		final DDLComparisonHeap[] heaps=new DDLComparisonHeap[nQueries];

		long tc0=System.nanoTime();
		CompareThread[] workers=new CompareThread[threads];
		for(int wi=0; wi<threads; wi++){
			workers[wi]=new CompareThread(queries, refs, nextQuery, heaps,
				nQueries, nRefs, k, maxRecords, minHits, useIndex, index, totalComparisonsPerformed);
			workers[wi].start();
		}
		for(CompareThread w : workers){try{w.join();}catch(InterruptedException e){}}

		long tc1=System.nanoTime();
		long bruteForce=(long)nQueries*nRefs;
		long performed=totalComparisonsPerformed.get();
		System.err.println("Compared "+nQueries+" queries x "+nRefs+" refs in "+String.format("%.3f", (tc1-tc0)*1e-9)+" seconds.");
		System.err.println("Comparisons performed: "+performed+" / "+bruteForce+" brute-force"+
			(useIndex ? " (index efficiency: "+String.format("%.4f%%", (1.0-(double)performed/bruteForce)*100)+")" : ""));

		ByteBuilder bb=new ByteBuilder();
		boolean json=formatter.format==DDLFormatter.FORMAT_JSON;
		if(json){formatter.jsonStart(bb);}else{formatter.header(bb);}
		for(int qi=0; qi<nQueries; qi++){
			if(heaps[qi]==null){continue;}
			ArrayList<DDLComparison> results=heaps[qi].toList();
			if(useSSU){alignSSU(results);}
			for(DDLComparison c : results){
				if(c.equal<minHits){continue;}
				formatter.format(c, bb);
			}
		}
		if(json){formatter.jsonEnd(bb);}
		System.out.print(bb);

		t.stop();
		System.err.println("Total time: \t"+t);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Compare Thread         ----------------*/
	/*--------------------------------------------------------------*/

	static class CompareThread extends Thread {
		final ArrayList<DDLRecord> queries, refs;
		final AtomicLong nextQuery;
		final DDLComparisonHeap[] heaps;
		final int nQueries, nRefs, k, maxRecords, minHits;
		final boolean useIndex;
		final DDLIndex idx;
		final AtomicLong totalComparisonsPerformed;

		CompareThread(ArrayList<DDLRecord> queries, ArrayList<DDLRecord> refs,
				AtomicLong nextQuery, DDLComparisonHeap[] heaps,
				int nQueries, int nRefs, int k, int maxRecords, int minHits,
				boolean useIndex, DDLIndex idx, AtomicLong totalComparisonsPerformed){
			this.queries=queries; this.refs=refs;
			this.nextQuery=nextQuery; this.heaps=heaps;
			this.nQueries=nQueries; this.nRefs=nRefs;
			this.k=k; this.maxRecords=maxRecords; this.minHits=minHits;
			this.useIndex=useIndex; this.idx=idx;
			this.totalComparisonsPerformed=totalComparisonsPerformed;
		}

		@Override
		public void run(){
			long localComparisons=0;
			DDLComparison working=new DDLComparison();
			while(true){
				int qi=(int)nextQuery.getAndIncrement();
				if(qi>=nQueries){break;}
				DDLComparisonHeap heap=new DDLComparisonHeap(maxRecords);
				DDLRecord qRec=queries.get(qi);
				if(useIndex){
					int[] counts=idx.query(qRec.ddl);
					for(int ri=0; ri<nRefs; ri++){
						if(counts[ri]<minHits){continue;}
						working.compare(qRec, refs.get(ri), k);
						heap.offer(working);
						localComparisons++;
					}
				}else{
					for(int ri=0; ri<nRefs; ri++){
						working.compare(qRec, refs.get(ri), k);
						heap.offer(working);
						localComparisons++;
					}
				}
				heaps[qi]=heap;
			}
			totalComparisonsPerformed.addAndGet(localComparisons);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        SSU Alignment         ----------------*/
	/*--------------------------------------------------------------*/

	/** Aligns SSU sequences for finalized comparison results.
	 * Only aligns when both query and ref have matching SSU type (16S or 18S). */
	private static void alignSSU(ArrayList<DDLComparison> results){
		QuantumAligner aligner=new QuantumAligner();
		for(DDLComparison c : results){
			if(c.queryRecord==null || c.refRecord==null){continue;}
			byte[] q=null, r=null;
			if(c.queryRecord.r16S!=null && c.refRecord.r16S!=null){
				q=c.queryRecord.r16S; r=c.refRecord.r16S;
			}else if(c.queryRecord.r18S!=null && c.refRecord.r18S!=null){
				q=c.queryRecord.r18S; r=c.refRecord.r18S;
			}
			if(q!=null){c.ssuIdentity=aligner.align(q, r);}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Collision Test          ----------------*/
	/*--------------------------------------------------------------*/

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

	/*--------------------------------------------------------------*/
	/*----------------         Output Helper        ----------------*/
	/*--------------------------------------------------------------*/

	private static void printResults(ArrayList<DDLComparison> results, int minHits, DDLFormatter formatter){
		ByteBuilder bb=new ByteBuilder();
		if(formatter.format==DDLFormatter.FORMAT_JSON){
			formatter.jsonStart(bb);
			for(DDLComparison c : results){
				if(c.equal<minHits){continue;}
				formatter.format(c, bb);
			}
			formatter.jsonEnd(bb);
		}else{
			formatter.header(bb);
			for(DDLComparison c : results){
				if(c.equal<minHits){continue;}
				formatter.format(c, bb);
			}
		}
		System.out.print(bb);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Helpers             ----------------*/
	/*--------------------------------------------------------------*/

	/** Compares two DDL register arrays using only the exponent bits. */
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
