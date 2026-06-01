package ddl;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;

import bin.GeneTools;
import cardinality.DynamicDemiLog;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import idaligner.QuantumAligner;
import map.IntObjectMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.CallGenes;
import prok.GeneCaller;
import prok.Orf;
import prok.ProkObject;
import shared.Shared;
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

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		Timer t=new Timer();
		DDLCompare x=new DDLCompare(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public DDLCompare(String[] args){

		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;

		{
			final Parser parser=parse(args);
			Parser.processQuality();
			overwrite=parser.overwrite;
			if(parser.in1!=null){file1=parser.in1;}
			if(parser.in2!=null){file2=parser.in2;}
		}

		if("refseq".equals(refFile)){refFile=shared.Resources.find("?refseqSketchDDL_k25e5b4096.tsv.gz");}

		if(blacklistSet){
			if(blacklistFile!=null){DynamicDemiLog.loadBlacklist(blacklistFile);}
		}else if(!DynamicDemiLog.blacklistExists()){
			String blPath=dna.Data.findPath("?genomeDDLBlacklist_k25e5b4096.fa.gz", false);
			if(blPath!=null){DynamicDemiLog.loadBlacklist(blPath);}
		}

		if(queryFile!=null){formatter.printQueryName=true;}
		if(useSSU){formatter.printSSU=true;}
	}

	private Parser parse(String[] args){
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(formatter.parse(arg, a, b)){/* handled by formatter */}
			else if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("exponent") || a.equals("ebits")){DynamicDemiLog.setExponent(Integer.parseInt(b));}
			else if(a.equals("buckets")){buckets=Integer.parseInt(b);}
			else if(a.equals("ref")){refFile=b;}
			else if(a.equals("refseq")){refFile="refseq";}
			else if(a.equals("query")){file1=b;}
			else if(a.equals("queryfile") || a.equals("qf")){queryFile=b;}
			else if(a.equals("records") || a.equals("maxrecords")){maxRecords=Integer.parseInt(b);}
			else if(a.equals("minhits")){minHits=Integer.parseInt(b);}
			else if(a.equals("collisiontest")){collisionTest=true;}
			else if(a.equals("index")){useIndex=Parse.parseBoolean(b);}
			else if(a.equals("ssu")){useSSU=Parse.parseBoolean(b);}
			else if(a.equals("percontig") || a.equals("persequence")){perContig=Parse.parseBoolean(b);}
			else if(a.equals("perfile")){perContig=!Parse.parseBoolean(b);}
			else if(a.equals("minsketch") || a.equals("minsketchlen")){minSketchLen=Integer.parseInt(b);}
			else if(a.equals("parallelload") || a.equals("pload")){parallelLoad=Parse.parseBoolean(b);}
			else if(a.equals("blacklist")){blacklistFile=b; blacklistSet=true;}
			else if(parser.parse(arg, a, b)){/* handled by parser */}
			else if(file1==null){file1=arg;}
			else if(file2==null){file2=arg;}
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		return parser;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	void process(Timer t){
		final int threads=Shared.threads();

		if(collisionTest && refFile!=null){
			collisionTest(refFile, threads);
		}else if(queryFile!=null && refFile!=null){
			compareMultiQuery(threads);
		}else if(refFile!=null){
			compareToRefs(threads);
		}else{
			comparePairwise(threads);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Pairwise Mode          ----------------*/
	/*--------------------------------------------------------------*/

	private void comparePairwise(int threads){
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

		outstream.println(String.format("ANI: %.2f%%  WKID: %.4f  Comp(1->2): %.4f  Comp(2->1): %.4f  (%d/%d buckets match)",
			ani*100, c, compAB, compBA, equal, lower+equal+higher));
	}

	/*--------------------------------------------------------------*/
	/*----------------     Single-Query Mode        ----------------*/
	/*--------------------------------------------------------------*/

	private void compareToRefs(int threads){
		Timer t=new Timer();
		outstream.println("Mode: "+(perContig ? "percontig " : "perfile ")
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
		outstream.println("Loading references from "+refFile+"...");
		ArrayList<DDLRecord> refs=DDLLoaderMT.loadFile(refFile, k, threads);
		long tRefLoad=System.nanoTime()-ts;
		outstream.println("Loaded "+refs.size()+" reference DDLs in "+fmt(tRefLoad)+" seconds.");

		if(!refs.isEmpty()){
			final int refBuckets=refs.get(0).ddl.buckets;
			if(refBuckets!=buckets){
				outstream.println("Auto-setting query buckets="+refBuckets+" to match reference database.");
				buckets=refBuckets;
			}
		}

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
			outstream.println("SSU load"+(ssuThread!=null ? " (parallel)" : "")+": "+fmt(tSSULoad)+"s"
				+(ssuThread!=null ? " (attach only; maps loaded in parallel with refs)" : ""));
		}

		/* --- Phase 2: Build index --- */

		long tIndex=0;
		DDLIndex index=null;
		if(useIndex){
			ts=System.nanoTime();
			index=new DDLIndex(refs.get(0).ddl.buckets);
			index.addAll(refs, threads);
			tIndex=System.nanoTime()-ts;
			outstream.println("Built inverted index in "+fmt(tIndex)+" seconds.");
		}

		/* --- Phase 3: Load query --- */

		ts=System.nanoTime();
		final ArrayList<DDLRecord> queries;
		if(perContig){
			queries=DDLQueryLoaderSF.loadPerContig(file1, buckets, k, useSSU, threads, minSketchLen);
			outstream.println("Loaded "+queries.size()+" contigs from "+file1
				+" (minsketch="+minSketchLen+") in "+fmt(System.nanoTime()-ts)+"s  threads="+threads);
		}else{
			DDLRecord qRec=DDLQueryLoaderSF.loadPerFile(file1, buckets, k, useSSU, threads);
			queries=new ArrayList<>(1);
			queries.add(qRec);
			outstream.println("Query: "+file1+"  bases="+qRec.bases+"  cardinality="+qRec.cardinality
				+"  sketch: "+fmt(System.nanoTime()-ts)+"s  threads="+threads
				+(useSSU ? "  ssu="+(qRec.hasRibo() ? DDLRecord.riboName(qRec.riboType())+"("+qRec.riboBytes().length+"bp)" : "none") : ""));
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
			outstream.println("Compared against "+nRefs+" references in "+fmt(tCompare)+" seconds.");

			ArrayList<DDLComparison> results=heap.toList();
			long tAlign=0;
			int alignCount=0;
			if(useSSU){
				ts=System.nanoTime();
				alignCount=alignSSU(results);
				tAlign=System.nanoTime()-ts;
				outstream.println("SSU alignments: "+alignCount+" in "+fmt(tAlign)+" seconds.");
			}

			printResults(results);

			t.stop();
			outstream.println("\nSubphase timing:");
			outstream.println("  Ref load:    \t"+fmt(tRefLoad)+"s");
			if(useSSU){outstream.println("  SSU load:    \t"+fmt(tSSULoad)+"s");}
			if(useIndex){outstream.println("  Index build: \t"+fmt(tIndex)+"s");}
			outstream.println("  Query load:  \t"+fmt(tQueryLoad)+"s");
			outstream.println("  Compare:     \t"+fmt(tCompare)+"s");
			if(useSSU){outstream.println("  SSU align:   \t"+fmt(tAlign)+"s  ("+alignCount+" alignments)");}
			outstream.println("Total time: \t"+t);
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
			outstream.println("Compared "+nQueries+" queries x "+nRefs+" refs in "+fmt(tCompare)+" seconds."
				+(useSSU ? "  ("+totalAlignments+" SSU alignments inline)" : ""));
			outstream.println("Comparisons performed: "+performed+" / "+bruteForce+" brute-force"+
				(useIndex ? " (index efficiency: "+String.format("%.4f%%", (1.0-(double)performed/bruteForce)*100)+")" : ""));

			ts=System.nanoTime();
			ByteBuilder bb=new ByteBuilder();
			boolean json=formatter.format==DDLFormatter.FORMAT_JSON;
			if(json){formatter.jsonStart(bb);}else{formatter.header(bb);}
			for(int qi=0; qi<nQueries; qi++){
				if(allResults[qi]==null){continue;}
				for(DDLComparison dc : allResults[qi]){
					if(dc.equal<minHits){continue;}
					formatter.format(dc, bb);
				}
			}
			if(json){formatter.jsonEnd(bb);}
			System.out.print(bb);
			long tFormat=System.nanoTime()-ts;

			t.stop();
			outstream.println("\nSubphase timing:");
			outstream.println("  Ref load:      \t"+fmt(tRefLoad)+"s");
			if(useSSU){outstream.println("  SSU load:      \t"+fmt(tSSULoad)+"s");}
			if(useIndex){outstream.println("  Index build:   \t"+fmt(tIndex)+"s");}
			outstream.println("  Query load:    \t"+fmt(tQueryLoad)+"s  ("+nQueries+" queries)");
			outstream.println("  Compare+align: \t"+fmt(tCompare)+"s"+(useSSU ? "  ("+totalAlignments+" alignments)" : ""));
			outstream.println("  Format:        \t"+fmt(tFormat)+"s");
			outstream.println("Total time: \t"+t);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------     Multi-Query Mode         ----------------*/
	/*--------------------------------------------------------------*/

	private void compareMultiQuery(int threads){
		Timer t=new Timer();
		outstream.println("Mode: multi-query "+(useIndex ? "indexed" : "pairwise")+"  threads="+threads+(useSSU ? "  ssu=t" : ""));

		long ts=System.nanoTime();
		outstream.println("Loading references from "+refFile+"...");
		ArrayList<DDLRecord> refs=DDLLoaderMT.loadFile(refFile, k, threads);
		long tRefLoad=System.nanoTime()-ts;
		outstream.println("Loaded "+refs.size()+" reference DDLs in "+fmt(tRefLoad)+" seconds.");

		long tSSULoad=0;
		if(useSSU){
			ts=System.nanoTime();
			DDLSSULoader.loadAndAttachDefaults(refs);
			tSSULoad=System.nanoTime()-ts;
			outstream.println("SSU ref load+attach: "+fmt(tSSULoad)+" seconds.");
		}

		ts=System.nanoTime();
		outstream.println("Loading queries from "+queryFile+"...");
		ArrayList<DDLRecord> queries=DDLLoaderMT.loadFile(queryFile, k, threads);
		long tQueryLoad=System.nanoTime()-ts;
		outstream.println("Loaded "+queries.size()+" query DDLs in "+fmt(tQueryLoad)+" seconds.");

		long tSSUQueryLoad=0;
		if(useSSU){
			ts=System.nanoTime();
			DDLSSULoader.loadAndAttachDefaults(queries);
			tSSUQueryLoad=System.nanoTime()-ts;
			outstream.println("SSU query load+attach: "+fmt(tSSUQueryLoad)+" seconds.");
		}

		final int nRefs=refs.size();
		final int nQueries=queries.size();

		long tIndex=0;
		DDLIndex index=null;
		if(useIndex){
			ts=System.nanoTime();
			index=new DDLIndex(refs.get(0).ddl.buckets);
			index.addAll(refs, threads);
			tIndex=System.nanoTime()-ts;
			outstream.println("Built inverted index in "+fmt(tIndex)+" seconds.");
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
		outstream.println("Compared "+nQueries+" queries x "+nRefs+" refs in "+fmt(tCompare)+" seconds."
			+(useSSU ? "  ("+totalAlignments+" SSU alignments inline)" : ""));
		outstream.println("Comparisons performed: "+performed+" / "+bruteForce+" brute-force"+
			(useIndex ? " (index efficiency: "+String.format("%.4f%%", (1.0-(double)performed/bruteForce)*100)+")" : ""));

		ts=System.nanoTime();
		ByteBuilder bb=new ByteBuilder();
		boolean json=formatter.format==DDLFormatter.FORMAT_JSON;
		if(json){formatter.jsonStart(bb);}else{formatter.header(bb);}
		for(int qi=0; qi<nQueries; qi++){
			if(allResults[qi]==null){continue;}
			for(DDLComparison dc : allResults[qi]){
				if(dc.equal<minHits){continue;}
				formatter.format(dc, bb);
			}
		}
		if(json){formatter.jsonEnd(bb);}
		System.out.print(bb);
		long tFormat=System.nanoTime()-ts;

		t.stop();
		outstream.println("\nSubphase timing:");
		outstream.println("  Ref load:      \t"+fmt(tRefLoad)+"s");
		if(useSSU){outstream.println("  SSU ref load:  \t"+fmt(tSSULoad)+"s");}
		outstream.println("  Query load:    \t"+fmt(tQueryLoad)+"s");
		if(useSSU){outstream.println("  SSU query load:\t"+fmt(tSSUQueryLoad)+"s");}
		if(useIndex){outstream.println("  Index build:   \t"+fmt(tIndex)+"s");}
		outstream.println("  Compare+align: \t"+fmt(tCompare)+"s"+(useSSU ? "  ("+totalAlignments+" alignments)" : ""));
		outstream.println("  Format:        \t"+fmt(tFormat)+"s");
		outstream.println("Total time: \t"+t);
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

	private static int alignSSU(ArrayList<DDLComparison> results){
		QuantumAligner aligner=new QuantumAligner();
		int count=0;
		for(DDLComparison dc : results){
			if(dc.queryRecord==null || dc.refRecord==null){continue;}
			int qt=dc.queryRecord.riboType(), rt=dc.refRecord.riboType();
			if(qt!=DDLRecord.RIBO_NONE && qt==rt){
				dc.ssuIdentity=aligner.align(dc.queryRecord.riboBytes(), dc.refRecord.riboBytes());
				count++;
			}
		}
		return count;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Collision Test          ----------------*/
	/*--------------------------------------------------------------*/

	private void collisionTest(String refPath, int threads){
		outstream.println("Loading references from "+refPath+"...");
		ArrayList<DDLRecord> refs=DDLLoader.loadFile(refPath, k);
		final int n=refs.size();
		outstream.println("Loaded "+n+" DDLs. Testing "+n*(n-1)/2+" pairs...");

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
		outstream.println("Pairs: "+totalPairs);
		outstream.println("Total matches: "+totalMatches);
		outstream.println("Avg matches/pair: "+String.format("%.4f", (double)totalMatches/totalPairs));
		outstream.println("Max matches in any pair: "+maxMatches);
		outstream.println("Avg collision rate/bucket: "+String.format("%.6f", (double)totalMatches/(totalPairs*2048)));
	}

	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private void printResults(ArrayList<DDLComparison> results){
		ByteBuilder bb=new ByteBuilder();
		if(formatter.format==DDLFormatter.FORMAT_JSON){
			formatter.jsonStart(bb);
			for(DDLComparison dc : results){
				if(dc.equal<minHits){continue;}
				formatter.format(dc, bb);
			}
			formatter.jsonEnd(bb);
		}else{
			formatter.header(bb);
			for(DDLComparison dc : results){
				if(dc.equal<minHits){continue;}
				formatter.format(dc, bb);
			}
		}
		System.out.print(bb);
	}

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

	/** Hashes a FASTA file into a DDL, optionally gene-calling for SSU inline. */
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

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String file1=null;
	private String file2=null;
	private String refFile=null;
	private String queryFile=null;

	private int k=25;
	private int buckets=2048;
	private int maxRecords=20;
	private int minHits=5;

	private boolean collisionTest=false;
	private boolean useIndex=false;
	private boolean useSSU=false;
	private boolean perContig=false;
	private boolean parallelLoad=true;
	private int minSketchLen=400;
	private boolean overwrite=false;
	private String blacklistFile=null;
	private boolean blacklistSet=false;

	private final DDLFormatter formatter=new DDLFormatter();
	private PrintStream outstream=System.err;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	private static final int MIN_GENE_CALL_LENGTH=800;
}
