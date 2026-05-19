package ddl;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;

import bin.GeneTools;
import cardinality.DynamicDemiLog;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import idaligner.QuantumAligner;
import map.IntObjectMap;
import prok.CallGenes;
import prok.GeneCaller;
import prok.Orf;
import prok.ProkObject;
import shared.Timer;
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
		boolean perContig=false, parallelLoad=true;
		int minSketchLen=400;
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
			else if(a.equals("percontig") || a.equals("persequence")){perContig=b==null || b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true");}
			else if(a.equals("perfile")){perContig=!(b==null || b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true"));}
			else if(a.equals("minsketch") || a.equals("minsketchlen")){minSketchLen=Integer.parseInt(b);}
			else if(a.equals("parallelload") || a.equals("pload")){parallelLoad=b==null || b.equalsIgnoreCase("t") || b.equalsIgnoreCase("true");}
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
			compareToRefs(file1, refFile, k, buckets, maxRecords, minHits, useIndex,
				useSSU, perContig, minSketchLen, parallelLoad, threads, formatter);
			return;
		}

		comparePairwise(file1, file2, k, buckets, threads);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Pairwise Mode          ----------------*/
	/*--------------------------------------------------------------*/

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
			basesA=hashFile(file1, ddlA, k, null, null);
			basesB=hashFile(file2, ddlB, k, null, null);
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

	private static void compareToRefs(String queryPath, String refPath, int k, int buckets,
			int maxRecords, int minHits, boolean useIndex, boolean useSSU,
			boolean perContig, int minSketchLen, boolean parallelLoad, int threads, DDLFormatter formatter){
		Timer t=new Timer();
		System.err.println("Mode: "+(perContig ? "percontig " : "perfile ")
			+(useIndex ? "indexed" : "pairwise")+"  threads="+threads+(useSSU ? "  ssu=t" : "")
			+(parallelLoad && useSSU ? "  parallelLoad" : ""));

		/* --- Phase 1: Load refs, SSU maps, PGM --- */

		if(useSSU){
			ProkObject.callCDS=false;
			ProkObject.calltRNA=false;
			ProkObject.call23S=false;
			ProkObject.call5S=false;
		}

		SSUPGMLoadThread ssuThread=null;
		if(useSSU && parallelLoad){
			ssuThread=new SSUPGMLoadThread();
			ssuThread.start();
		}

		long ts=System.nanoTime();
		System.err.println("Loading references from "+refPath+"...");
		ArrayList<DDLRecord> refs=DDLLoaderMT.loadFile(refPath, k, threads);
		long tRefLoad=System.nanoTime()-ts;
		System.err.println("Loaded "+refs.size()+" reference DDLs in "+fmt(tRefLoad)+" seconds.");

		long tSSULoad=0;
		if(useSSU){
			ts=System.nanoTime();
			IntObjectMap<byte[]> map16, map18;
			if(ssuThread!=null){
				while(ssuThread.getState()!=Thread.State.TERMINATED){
					try{ssuThread.join();}catch(InterruptedException e){e.printStackTrace();}
				}
				map16=ssuThread.map16S;
				map18=ssuThread.map18S;
			}else{
				IntObjectMap<byte[]>[] maps=DDLSSULoader.loadSSUMapsDefaults();
				map16=maps[0]; map18=maps[1];
				GeneTools.loadPGM();
			}
			DDLSSULoader.attachSSU(refs, map16, map18);
			tSSULoad=System.nanoTime()-ts;
			System.err.println("SSU load"+(ssuThread!=null ? " (parallel)" : "")+": "+fmt(tSSULoad)+"s"
				+(ssuThread!=null ? " (attach only; maps loaded in parallel with refs)" : ""));
		}

		/* --- Phase 2: Build index --- */

		long tIndex=0;
		DDLIndex index=null;
		if(useIndex){
			ts=System.nanoTime();
			index=new DDLIndex();
			index.addAll(refs, threads);
			tIndex=System.nanoTime()-ts;
			System.err.println("Built inverted index in "+fmt(tIndex)+" seconds.");
		}

		/* --- Phase 3: Load query --- */

		ts=System.nanoTime();
		final ArrayList<DDLRecord> queries;
		if(perContig){
			queries=DDLQueryLoaderSF.loadPerContig(queryPath, buckets, k, useSSU, threads, minSketchLen);
			System.err.println("Loaded "+queries.size()+" contigs from "+queryPath
				+" (minsketch="+minSketchLen+") in "+fmt(System.nanoTime()-ts)+"s  threads="+threads);
		}else{
			DDLRecord qRec=DDLQueryLoaderSF.loadPerFile(queryPath, buckets, k, useSSU, threads);
			queries=new ArrayList<>(1);
			queries.add(qRec);
			System.err.println("Query: "+queryPath+"  bases="+qRec.bases+"  cardinality="+qRec.cardinality
				+"  sketch: "+fmt(System.nanoTime()-ts)+"s  threads="+threads
				+(useSSU ? "  ssu="+(qRec.hasSSU() ? (qRec.r16S!=null ? "16S("+qRec.r16S.length+"bp)" : "18S("+qRec.r18S.length+"bp)") : "none") : ""));
		}
		long tQueryLoad=System.nanoTime()-ts;

		/* --- Phase 4: Compare --- */

		final int nRefs=refs.size();
		final int nQueries=queries.size();

		if(nQueries==1 && !perContig){
			DDLRecord qRec=queries.get(0);
			DDLComparisonHeap heap=new DDLComparisonHeap(maxRecords);
			DDLComparison working=new DDLComparison();

			ts=System.nanoTime();
			if(useIndex){
				int[] counts=index.query(qRec.ddl);
				for(int ri=0; ri<nRefs; ri++){
					if(counts[ri]<minHits){continue;}
					working.compare(qRec, refs.get(ri), k);
					heap.offer(working);
				}
			}else{
				for(int ri=0; ri<nRefs; ri++){
					working.compare(qRec, refs.get(ri), k);
					heap.offer(working);
				}
			}
			long tCompare=System.nanoTime()-ts;
			System.err.println("Compared against "+nRefs+" references in "+fmt(tCompare)+" seconds.");

			ArrayList<DDLComparison> results=heap.toList();
			long tAlign=0;
			int alignCount=0;
			if(useSSU){
				ts=System.nanoTime();
				alignCount=alignSSU(results);
				tAlign=System.nanoTime()-ts;
				System.err.println("SSU alignments: "+alignCount+" in "+fmt(tAlign)+" seconds.");
			}

			printResults(results, minHits, formatter);

			t.stop();
			System.err.println("\nSubphase timing:");
			System.err.println("  Ref load:    \t"+fmt(tRefLoad)+"s");
			if(useSSU){System.err.println("  SSU load:    \t"+fmt(tSSULoad)+"s");}
			if(useIndex){System.err.println("  Index build: \t"+fmt(tIndex)+"s");}
			System.err.println("  Query load:  \t"+fmt(tQueryLoad)+"s");
			System.err.println("  Compare:     \t"+fmt(tCompare)+"s");
			if(useSSU){System.err.println("  SSU align:   \t"+fmt(tAlign)+"s  ("+alignCount+" alignments)");}
			System.err.println("Total time: \t"+t);
		}else{
			/* Multi-query compare (perContig or single query dispatched here) */
			if(perContig){formatter.printQueryName=true;}
			@SuppressWarnings("unchecked")
			final ArrayList<DDLComparison>[] allResults=new ArrayList[nQueries];
			final AtomicLong nextQuery=new AtomicLong(0);
			final AtomicLong totalComparisonsPerformed=new AtomicLong(0);

			ts=System.nanoTime();
			CompareThread[] workers=new CompareThread[threads];
			for(int wi=0; wi<threads; wi++){
				workers[wi]=new CompareThread(queries, refs, nextQuery, allResults,
					nQueries, nRefs, k, maxRecords, minHits, useIndex, index,
					totalComparisonsPerformed, useSSU);
				workers[wi].start();
			}
			for(CompareThread w : workers){try{w.join();}catch(InterruptedException e){}}
			long tCompare=System.nanoTime()-ts;

			long bruteForce=(long)nQueries*nRefs;
			long performed=totalComparisonsPerformed.get();
			long totalAlignments=0;
			for(CompareThread w : workers){totalAlignments+=w.alignCountT;}
			System.err.println("Compared "+nQueries+" queries x "+nRefs+" refs in "+fmt(tCompare)+" seconds."
				+(useSSU ? "  ("+totalAlignments+" SSU alignments inline)" : ""));
			System.err.println("Comparisons performed: "+performed+" / "+bruteForce+" brute-force"+
				(useIndex ? " (index efficiency: "+String.format("%.4f%%", (1.0-(double)performed/bruteForce)*100)+")" : ""));

			ts=System.nanoTime();
			ByteBuilder bb=new ByteBuilder();
			boolean json=formatter.format==DDLFormatter.FORMAT_JSON;
			if(json){formatter.jsonStart(bb);}else{formatter.header(bb);}
			for(int qi=0; qi<nQueries; qi++){
				if(allResults[qi]==null){continue;}
				for(DDLComparison c : allResults[qi]){
					if(c.equal<minHits){continue;}
					formatter.format(c, bb);
				}
			}
			if(json){formatter.jsonEnd(bb);}
			System.out.print(bb);
			long tFormat=System.nanoTime()-ts;

			t.stop();
			System.err.println("\nSubphase timing:");
			System.err.println("  Ref load:      \t"+fmt(tRefLoad)+"s");
			if(useSSU){System.err.println("  SSU load:      \t"+fmt(tSSULoad)+"s");}
			if(useIndex){System.err.println("  Index build:   \t"+fmt(tIndex)+"s");}
			System.err.println("  Query load:    \t"+fmt(tQueryLoad)+"s  ("+nQueries+" queries)");
			System.err.println("  Compare+align: \t"+fmt(tCompare)+"s"+(useSSU ? "  ("+totalAlignments+" alignments)" : ""));
			System.err.println("  Format:        \t"+fmt(tFormat)+"s");
			System.err.println("Total time: \t"+t);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------     Multi-Query Mode         ----------------*/
	/*--------------------------------------------------------------*/

	private static void compareMultiQuery(String queryPath, String refPath, int k,
			int maxRecords, int minHits, boolean useIndex, boolean useSSU, int threads, DDLFormatter formatter){
		Timer t=new Timer();
		System.err.println("Mode: multi-query "+(useIndex ? "indexed" : "pairwise")+"  threads="+threads+(useSSU ? "  ssu=t" : ""));

		long ts=System.nanoTime();
		System.err.println("Loading references from "+refPath+"...");
		ArrayList<DDLRecord> refs=DDLLoaderMT.loadFile(refPath, k, threads);
		long tRefLoad=System.nanoTime()-ts;
		System.err.println("Loaded "+refs.size()+" reference DDLs in "+fmt(tRefLoad)+" seconds.");

		long tSSULoad=0;
		if(useSSU){
			ts=System.nanoTime();
			DDLSSULoader.loadAndAttachDefaults(refs);
			tSSULoad=System.nanoTime()-ts;
			System.err.println("SSU ref load+attach: "+fmt(tSSULoad)+" seconds.");
		}

		ts=System.nanoTime();
		System.err.println("Loading queries from "+queryPath+"...");
		ArrayList<DDLRecord> queries=DDLLoaderMT.loadFile(queryPath, k, threads);
		long tQueryLoad=System.nanoTime()-ts;
		System.err.println("Loaded "+queries.size()+" query DDLs in "+fmt(tQueryLoad)+" seconds.");

		long tSSUQueryLoad=0;
		if(useSSU){
			ts=System.nanoTime();
			DDLSSULoader.loadAndAttachDefaults(queries);
			tSSUQueryLoad=System.nanoTime()-ts;
			System.err.println("SSU query load+attach: "+fmt(tSSUQueryLoad)+" seconds.");
		}

		final int nRefs=refs.size();
		final int nQueries=queries.size();

		long tIndex=0;
		DDLIndex index=null;
		if(useIndex){
			ts=System.nanoTime();
			index=new DDLIndex();
			index.addAll(refs, threads);
			tIndex=System.nanoTime()-ts;
			System.err.println("Built inverted index in "+fmt(tIndex)+" seconds.");
		}

		@SuppressWarnings("unchecked")
		final ArrayList<DDLComparison>[] allResults=new ArrayList[nQueries];
		final AtomicLong nextQuery=new AtomicLong(0);
		final AtomicLong totalComparisonsPerformed=new AtomicLong(0);

		ts=System.nanoTime();
		CompareThread[] workers=new CompareThread[threads];
		for(int wi=0; wi<threads; wi++){
			workers[wi]=new CompareThread(queries, refs, nextQuery, allResults,
				nQueries, nRefs, k, maxRecords, minHits, useIndex, index,
				totalComparisonsPerformed, useSSU);
			workers[wi].start();
		}
		for(CompareThread w : workers){try{w.join();}catch(InterruptedException e){}}
		long tCompare=System.nanoTime()-ts;

		long bruteForce=(long)nQueries*nRefs;
		long performed=totalComparisonsPerformed.get();
		long totalAlignments=0;
		for(CompareThread w : workers){totalAlignments+=w.alignCountT;}
		System.err.println("Compared "+nQueries+" queries x "+nRefs+" refs in "+fmt(tCompare)+" seconds."
			+(useSSU ? "  ("+totalAlignments+" SSU alignments inline)" : ""));
		System.err.println("Comparisons performed: "+performed+" / "+bruteForce+" brute-force"+
			(useIndex ? " (index efficiency: "+String.format("%.4f%%", (1.0-(double)performed/bruteForce)*100)+")" : ""));

		ts=System.nanoTime();
		ByteBuilder bb=new ByteBuilder();
		boolean json=formatter.format==DDLFormatter.FORMAT_JSON;
		if(json){formatter.jsonStart(bb);}else{formatter.header(bb);}
		for(int qi=0; qi<nQueries; qi++){
			if(allResults[qi]==null){continue;}
			for(DDLComparison c : allResults[qi]){
				if(c.equal<minHits){continue;}
				formatter.format(c, bb);
			}
		}
		if(json){formatter.jsonEnd(bb);}
		System.out.print(bb);
		long tFormat=System.nanoTime()-ts;

		t.stop();
		System.err.println("\nSubphase timing:");
		System.err.println("  Ref load:      \t"+fmt(tRefLoad)+"s");
		if(useSSU){System.err.println("  SSU ref load:  \t"+fmt(tSSULoad)+"s");}
		System.err.println("  Query load:    \t"+fmt(tQueryLoad)+"s");
		if(useSSU){System.err.println("  SSU query load:\t"+fmt(tSSUQueryLoad)+"s");}
		if(useIndex){System.err.println("  Index build:   \t"+fmt(tIndex)+"s");}
		System.err.println("  Compare+align: \t"+fmt(tCompare)+"s"+(useSSU ? "  ("+totalAlignments+" alignments)" : ""));
		System.err.println("  Format:        \t"+fmt(tFormat)+"s");
		System.err.println("Total time: \t"+t);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Compare Thread         ----------------*/
	/*--------------------------------------------------------------*/

	static class CompareThread extends Thread {
		final ArrayList<DDLRecord> queries, refs;
		final AtomicLong nextQuery;
		final ArrayList<DDLComparison>[] results;
		final int nQueries, nRefs, k, maxRecords, minHits;
		final boolean useIndex, alignSSU;
		final DDLIndex idx;
		final AtomicLong totalComparisonsPerformed;
		long alignCountT;

		CompareThread(ArrayList<DDLRecord> queries, ArrayList<DDLRecord> refs,
				AtomicLong nextQuery, ArrayList<DDLComparison>[] results,
				int nQueries, int nRefs, int k, int maxRecords, int minHits,
				boolean useIndex, DDLIndex idx, AtomicLong totalComparisonsPerformed,
				boolean alignSSU){
			this.queries=queries; this.refs=refs;
			this.nextQuery=nextQuery; this.results=results;
			this.nQueries=nQueries; this.nRefs=nRefs;
			this.k=k; this.maxRecords=maxRecords; this.minHits=minHits;
			this.useIndex=useIndex; this.idx=idx;
			this.totalComparisonsPerformed=totalComparisonsPerformed;
			this.alignSSU=alignSSU;
		}

		@Override
		public void run(){
			long localComparisons=0;
			long localAligns=0;
			DDLComparison working=new DDLComparison();
			QuantumAligner aligner=(alignSSU ? new QuantumAligner() : null);
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
				ArrayList<DDLComparison> list=heap.toList();
				if(aligner!=null){
					for(DDLComparison c : list){
						if(c.queryRecord==null || c.refRecord==null){continue;}
						byte[] q=null, r=null;
						if(c.queryRecord.r16S!=null && c.refRecord.r16S!=null){
							q=c.queryRecord.r16S; r=c.refRecord.r16S;
						}else if(c.queryRecord.r18S!=null && c.refRecord.r18S!=null){
							q=c.queryRecord.r18S; r=c.refRecord.r18S;
						}
						if(q!=null){c.ssuIdentity=aligner.align(q, r); localAligns++;}
					}
				}
				results[qi]=list;
			}
			totalComparisonsPerformed.addAndGet(localComparisons);
			alignCountT=localAligns;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------     SSU/PGM Load Thread      ----------------*/
	/*--------------------------------------------------------------*/

	static class SSUPGMLoadThread extends Thread {
		@Override
		public void run(){
			IntObjectMap<byte[]>[] maps=DDLSSULoader.loadSSUMapsDefaults();
			map16S=maps[0];
			map18S=maps[1];
			GeneTools.loadPGM();
			caller=GeneTools.makeGeneCaller();
		}
		IntObjectMap<byte[]> map16S;
		IntObjectMap<byte[]> map18S;
		GeneCaller caller;
	}

	/*--------------------------------------------------------------*/
	/*----------------        SSU Alignment         ----------------*/
	/*--------------------------------------------------------------*/

	/** Single-threaded SSU alignment for small result sets (single-query path). */
	private static int alignSSU(ArrayList<DDLComparison> results){
		QuantumAligner aligner=new QuantumAligner();
		int count=0;
		for(DDLComparison c : results){
			if(c.queryRecord==null || c.refRecord==null){continue;}
			byte[] q=null, r=null;
			if(c.queryRecord.r16S!=null && c.refRecord.r16S!=null){
				q=c.queryRecord.r16S; r=c.refRecord.r16S;
			}else if(c.queryRecord.r18S!=null && c.refRecord.r18S!=null){
				q=c.queryRecord.r18S; r=c.refRecord.r18S;
			}
			if(q!=null){c.ssuIdentity=aligner.align(q, r); count++;}
		}
		return count;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Collision Test          ----------------*/
	/*--------------------------------------------------------------*/

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

	private static String fmt(long nanos){return String.format("%.3f", nanos*1e-9);}

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

	/** Hashes a FASTA file into a DDL, optionally gene-calling for SSU inline.
	 * @param caller GeneCaller for SSU detection (null to skip gene-calling)
	 * @param rec DDLRecord to attach found SSU to (null to skip) */
	private static long hashFile(String path, DynamicDemiLog ddl, int k,
			GeneCaller caller, DDLRecord rec){
		FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		Streamer cris=StreamerFactory.getReadInputStream(-1, false, ff, null, -1);
		cris.start();

		boolean needSSU=(caller!=null && rec!=null);
		long bases=0;
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		while(ln!=null && reads!=null && reads.size()>0){
			for(Read r : reads){
				ddl.hash(r);
				bases+=r.length();
				if(r.mate!=null){bases+=r.mate.length();}
				if(needSSU && r.length()>=MIN_GENE_CALL_LENGTH){
					ArrayList<Orf> genes=caller.callGenes(r);
					if(genes!=null){
						for(Orf orf : genes){
							if(orf.is16S()){rec.r16S=CallGenes.fetch(orf, r).bases; needSSU=false; break;}
							else if(orf.is18S()){rec.r18S=CallGenes.fetch(orf, r).bases; needSSU=false; break;}
						}
					}
				}
			}
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln);
		ReadWrite.closeStreams(cris);
		return bases;
	}

	private static final int MIN_GENE_CALL_LENGTH=800;
}
