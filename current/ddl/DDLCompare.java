package ddl;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;

import bin.GeneTools;
import cardinality.DynamicDemiLog;
import idaligner.QuantumAligner;
import map.IntObjectMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.ProkObject;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;

/**
 * Genome comparison against pre-built DDL reference databases.
 * Supports single-query, per-contig, and multi-query batch modes.
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

		{
			final Parser parser=parse(args);
			Parser.processQuality();
			overwrite=parser.overwrite;
			if(parser.in1!=null){file1=parser.in1;}
			if(parser.in2!=null){file2=parser.in2;}
		}

		if("refseq".equals(refFile)){refFile=shared.Resources.find(
			"?refseqSketchDDL_k25e5b4096.tsv.gz,?refseqSketchDDL_k25e5b4096_merged.tsv.gz,"
			+"?refseqSketchDDL_k25e5b2048.tsv.gz,?refseqSketchDDL_k25e5b2048_merged.tsv.gz");}

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
			collisionTest(refFile);
		}else if(queryFile!=null && refFile!=null){
			multiQueryMode(threads);
		}else if(refFile!=null){
			queryMode(threads);
		}else{
			DDLBenchmark.pairwiseTest(file1, file2, k, buckets, Shared.threads(), outstream);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Reference Loading       ----------------*/
	/*--------------------------------------------------------------*/

	private ArrayList<DDLRecord> loadRefs(int threads){
		outstream.println("Loading references from "+refFile+"...");
		ArrayList<DDLRecord> refs=DDLLoaderMT.loadFile(refFile, k, threads);
		if(!refs.isEmpty()){
			final int refBuckets=refs.get(0).ddl.buckets;
			if(refBuckets!=buckets){
				outstream.println("Auto-setting query buckets="+refBuckets+" to match reference database.");
				buckets=refBuckets;
			}
		}
		return refs;
	}

	private long loadSSURefs(ArrayList<DDLRecord> refs, SSUPGMLoadThread ssuThread){
		if(!useSSU){return 0;}
		long ts=System.nanoTime();
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
		long elapsed=System.nanoTime()-ts;
		outstream.println("SSU load"+(ssuThread!=null ? " (parallel)" : "")+": "+fmt(elapsed)+"s"
			+(ssuThread!=null ? " (attach only; maps loaded in parallel with refs)" : ""));
		return elapsed;
	}

	private DDLIndex buildIndex(ArrayList<DDLRecord> refs, int threads){
		if(!useIndex){return null;}
		DDLIndex index=new DDLIndex(refs.get(0).ddl.buckets);
		index.addAll(refs, threads);
		return index;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Query Mode            ----------------*/
	/*--------------------------------------------------------------*/

	private void queryMode(int threads){
		Timer t=new Timer();
		outstream.println("Mode: "+(perContig ? "percontig " : "perfile ")
			+(useIndex ? "indexed" : "pairwise")+"  threads="+threads+(useSSU ? "  ssu=t" : "")
			+(parallelLoad && useSSU ? "  parallelLoad" : ""));

		SSUPGMLoadThread ssuThread=null;
		if(useSSU){
			ProkObject.callCDS=false; ProkObject.calltRNA=false;
			ProkObject.call23S=false; ProkObject.call5S=false;
			if(parallelLoad){ssuThread=new SSUPGMLoadThread(); ssuThread.start();}
		}

		long ts=System.nanoTime();
		ArrayList<DDLRecord> refs=loadRefs(threads);
		long tRefLoad=System.nanoTime()-ts;
		outstream.println("Loaded "+refs.size()+" reference DDLs in "+fmt(tRefLoad)+" seconds.");

		long tSSULoad=loadSSURefs(refs, ssuThread);

		ts=System.nanoTime();
		DDLIndex index=buildIndex(refs, threads);
		long tIndex=System.nanoTime()-ts;
		if(useIndex){outstream.println("Built inverted index in "+fmt(tIndex)+" seconds.");}

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

		if(queries.size()==1 && !perContig){
			long[] times=compareSingleQuery(queries.get(0), refs, index);
			t.stop();
			printTimingSingle(t, tRefLoad, tSSULoad, tIndex, tQueryLoad, times[0], times[1], (int)times[2]);
		}else{
			if(perContig){formatter.printQueryName=true;}
			ts=System.nanoTime();
			CompareResult cr=gatherResults(queries, refs, index, threads);
			long tCompare=System.nanoTime()-ts;
			formatAndPrint(cr.results, queries.size());
			t.stop();
			printTimingMulti(t, tRefLoad, tSSULoad, tIndex, tQueryLoad, tCompare,
				cr.comparisons, cr.alignments, queries.size(), refs.size());
		}
	}

	/** Single-query fast path: direct comparison loop without CompareThread.
	 * @return long[]{tCompare, tAlign, alignCount} */
	private long[] compareSingleQuery(DDLRecord qRec, ArrayList<DDLRecord> refs, DDLIndex index){
		final int nRefs=refs.size();
		DDLComparisonHeap heap=new DDLComparisonHeap(maxRecords);
		DDLComparison working=new DDLComparison();

		long ts=System.nanoTime();
		if(index!=null){
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
		long tAlign=0; int alignCount=0;
		if(useSSU){
			ts=System.nanoTime();
			alignCount=alignSSU(results);
			tAlign=System.nanoTime()-ts;
			outstream.println("SSU alignments: "+alignCount+" in "+fmt(tAlign)+" seconds.");
		}

		printResults(results);
		return new long[]{tCompare, tAlign, alignCount};
	}

	/*--------------------------------------------------------------*/
	/*----------------      Multi-Query Mode        ----------------*/
	/*--------------------------------------------------------------*/

	private void multiQueryMode(int threads){
		Timer t=new Timer();
		outstream.println("Mode: multi-query "+(useIndex ? "indexed" : "pairwise")
			+"  threads="+threads+(useSSU ? "  ssu=t" : ""));

		long ts=System.nanoTime();
		ArrayList<DDLRecord> refs=loadRefs(threads);
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

		ts=System.nanoTime();
		DDLIndex index=buildIndex(refs, threads);
		long tIndex=System.nanoTime()-ts;
		if(useIndex){outstream.println("Built inverted index in "+fmt(tIndex)+" seconds.");}

		ts=System.nanoTime();
		CompareResult cr=gatherResults(queries, refs, index, threads);
		long tCompare=System.nanoTime()-ts;

		outstream.println("Compared "+queries.size()+" queries x "+refs.size()+" refs in "+fmt(tCompare)+" seconds."
			+(useSSU ? "  ("+cr.alignments+" SSU alignments inline)" : ""));
		long bruteForce=(long)queries.size()*refs.size();
		outstream.println("Comparisons performed: "+cr.comparisons+" / "+bruteForce+" brute-force"
			+(useIndex ? " (index efficiency: "+String.format("%.4f%%", (1.0-(double)cr.comparisons/bruteForce)*100)+")" : ""));

		ts=System.nanoTime();
		formatAndPrint(cr.results, queries.size());
		long tFormat=System.nanoTime()-ts;

		t.stop();
		outstream.println("\nSubphase timing:");
		outstream.println("  Ref load:      \t"+fmt(tRefLoad)+"s");
		if(useSSU){outstream.println("  SSU ref load:  \t"+fmt(tSSULoad)+"s");}
		outstream.println("  Query load:    \t"+fmt(tQueryLoad)+"s");
		if(useSSU){outstream.println("  SSU query load:\t"+fmt(tSSUQueryLoad)+"s");}
		if(useIndex){outstream.println("  Index build:   \t"+fmt(tIndex)+"s");}
		outstream.println("  Compare+align: \t"+fmt(tCompare)+"s"+(useSSU ? "  ("+cr.alignments+" alignments)" : ""));
		outstream.println("  Format:        \t"+fmt(tFormat)+"s");
		outstream.println("Total time: \t"+t);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Comparison            ----------------*/
	/*--------------------------------------------------------------*/

	private CompareResult gatherResults(ArrayList<DDLRecord> queries,
			ArrayList<DDLRecord> refs, DDLIndex index, int threads){
		final int nQueries=queries.size(), nRefs=refs.size();
		@SuppressWarnings("unchecked")
		final ArrayList<DDLComparison>[] allResults=new ArrayList[nQueries];
		final AtomicLong nextQuery=new AtomicLong(0);
		final AtomicLong totalComparisons=new AtomicLong(0);

		CompareThread[] workers=new CompareThread[threads];
		for(int wi=0; wi<threads; wi++){
			workers[wi]=new CompareThread(queries, refs, nextQuery, allResults,
				nQueries, nRefs, k, maxRecords, minHits, useIndex, index,
				totalComparisons, useSSU);
			workers[wi].start();
		}
		for(CompareThread w : workers){try{w.join();}catch(InterruptedException e){}}

		long alignments=0;
		for(CompareThread w : workers){alignments+=w.alignCountT;}
		return new CompareResult(allResults, totalComparisons.get(), alignments);
	}

	private static class CompareResult {
		final ArrayList<DDLComparison>[] results;
		final long comparisons, alignments;
		CompareResult(ArrayList<DDLComparison>[] r, long c, long a){
			results=r; comparisons=c; alignments=a;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Output              ----------------*/
	/*--------------------------------------------------------------*/

	private void formatAndPrint(ArrayList<DDLComparison>[] allResults, int nQueries){
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
	}

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

	private void printTimingSingle(Timer t, long tRefLoad, long tSSULoad, long tIndex,
			long tQueryLoad, long tCompare, long tAlign, int alignCount){
		outstream.println("\nSubphase timing:");
		outstream.println("  Ref load:    \t"+fmt(tRefLoad)+"s");
		if(useSSU){outstream.println("  SSU load:    \t"+fmt(tSSULoad)+"s");}
		if(useIndex){outstream.println("  Index build: \t"+fmt(tIndex)+"s");}
		outstream.println("  Query load:  \t"+fmt(tQueryLoad)+"s");
		outstream.println("  Compare:     \t"+fmt(tCompare)+"s");
		if(useSSU){outstream.println("  SSU align:   \t"+fmt(tAlign)+"s  ("+alignCount+" alignments)");}
		outstream.println("Total time: \t"+t);
	}

	private void printTimingMulti(Timer t, long tRefLoad, long tSSULoad, long tIndex,
			long tQueryLoad, long tCompare, long comparisons, long alignments,
			int nQueries, int nRefs){
		long bruteForce=(long)nQueries*nRefs;
		outstream.println("Compared "+nQueries+" queries x "+nRefs+" refs in "+fmt(tCompare)+" seconds."
			+(useSSU ? "  ("+alignments+" SSU alignments inline)" : ""));
		outstream.println("Comparisons performed: "+comparisons+" / "+bruteForce+" brute-force"
			+(useIndex ? " (index efficiency: "+String.format("%.4f%%", (1.0-(double)comparisons/bruteForce)*100)+")" : ""));

		outstream.println("\nSubphase timing:");
		outstream.println("  Ref load:      \t"+fmt(tRefLoad)+"s");
		if(useSSU){outstream.println("  SSU load:      \t"+fmt(tSSULoad)+"s");}
		if(useIndex){outstream.println("  Index build:   \t"+fmt(tIndex)+"s");}
		outstream.println("  Query load:    \t"+fmt(tQueryLoad)+"s  ("+nQueries+" queries)");
		outstream.println("  Compare+align: \t"+fmt(tCompare)+"s"+(useSSU ? "  ("+alignments+" alignments)" : ""));
		outstream.println("Total time: \t"+t);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Collision Test          ----------------*/
	/*--------------------------------------------------------------*/

	private void collisionTest(String refPath){
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
	/*----------------     SSU/PGM Load Thread      ----------------*/
	/*--------------------------------------------------------------*/

	static class SSUPGMLoadThread extends Thread {
		@Override
		public void run(){
			IntObjectMap<byte[]>[] maps=DDLSSULoader.loadSSUMapsDefaults();
			map16S=maps[0];
			map18S=maps[1];
			GeneTools.loadPGM();
		}
		IntObjectMap<byte[]> map16S;
		IntObjectMap<byte[]> map18S;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/

	private static String fmt(long nanos){return String.format("%.3f", nanos*1e-9);}

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
