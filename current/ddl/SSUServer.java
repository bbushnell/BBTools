package ddl;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.BindException;
import java.net.InetSocketAddress;
import java.util.ArrayList;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicLong;

import com.sun.net.httpserver.Headers;
import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;

import bin.GeneTools;
import cardinality.DynamicDemiLog;
import fileIO.FileFormat;
import idaligner.QuantumAligner;
import prok.ProkObject;
import server.ServerTools;
import shared.Resources;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.StreamerFactory;
import structures.ByteBuilder;

/**
 * Persistent HTTP server for SSU (16S/18S) classification.
 * Preloads 276k DDL reference sketches and inverted index at startup,
 * then serves classification queries over HTTP.
 *
 * Modeled after CladeServer. Body-prefix routing:
 *   //JSON\n  — request JSON output (strip, then check next prefix)
 *   //Call\n  — gene-calling mode (strip, remaining body is genome FASTA)
 *   //Status  — health check
 *   (default) — raw FASTA, SSU compare mode
 *
 * @author Chloe, Brian Bushnell
 * @date May 20, 2026
 */
public class SSUServer {

	public static void main(String[] args){
		Timer t=new Timer();
		SSUServer server=new SSUServer(args);
		server.loadReferences();
		t.stop();
		System.err.println("References loaded in "+t);
		server.startServer();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Constructor           ----------------*/
	/*--------------------------------------------------------------*/

	public SSUServer(String[] args){
		for(int i=0; i<args.length; i++){
			String[] split=args[i].split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(a.equals("port")){port=Integer.parseInt(b);}
			else if(a.equals("kill") || a.equals("killcode")){killCode=b;}
			else if(a.equals("prefix") || a.equals("addressprefix")){addressPrefix=b;}
			else if(a.equals("domain")){domain=b;}
			else if(a.equals("ref")){refFile=b;}
			else if(a.equals("ref16s") || a.equals("ref16")){ref16sFile=b;}
			else if(a.equals("ref18s") || a.equals("ref18")){ref18sFile=b;}
			else if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("exponent") || a.equals("ebits")){DynamicDemiLog.setExponent(Integer.parseInt(b));}
			else if(a.equals("buckets")){buckets=Integer.parseInt(b);}
			else if(a.equals("records") || a.equals("maxrecords")){maxRecords=Integer.parseInt(b);}
			else if(a.equals("minhits")){minHits=Integer.parseInt(b);}
			else if(a.equals("buffer")){buffer=Integer.parseInt(b);}
			else if(a.equals("maxsize")){maxSize=Long.parseLong(b);}
			else if(a.equals("t") || a.equals("threads")){threads=Integer.parseInt(b);}
			else if(a.equals("verbose")){verbose=true;}
			else if(a.equals("verbose2")){verbose2=true; verbose=true;}
		}
		DynamicDemiLog.setExponent(exponent);
	}

	/*--------------------------------------------------------------*/
	/*----------------      Reference Loading       ----------------*/
	/*--------------------------------------------------------------*/

	private void loadReferences(){
		System.err.println("SSUServer starting on port "+port+"  k="+k
			+"  exponent="+DynamicDemiLog.exponentBits()+"  buckets="+buckets);

		if(refFile==null && ref16sFile==null && ref18sFile==null){
			refFile=Resources.find("?ssuSketchDDL.tsv.gz");
		}
		if(refFile==null && ref16sFile==null && ref18sFile==null){
			System.err.println("ERROR: Required resource not found: ssuSketchDDL.tsv.gz");
			System.err.println("Download it from:");
			System.err.println("  https://sourceforge.net/projects/bbmap/files/Resources/");
			System.err.println("Place it in BBTools/resources/ and try again.");
			System.exit(1);
		}

		ProkObject.callCDS=false;
		ProkObject.calltRNA=false;
		ProkObject.call23S=false;
		ProkObject.call5S=false;

		maps=DDLSSULoader.loadSSUMapsDefaults();
		map.IntObjectMap<byte[]> map16=maps[0], map18=maps[1];

		Thread pgmThread=new Thread(()->{GeneTools.loadPGM();});
		pgmThread.setDaemon(true);
		pgmThread.start();

		if(ref16sFile!=null || ref18sFile!=null){
			refs=new ArrayList<>();
			if(ref16sFile!=null){
				ArrayList<DDLRecord> r16=DDLLoaderMT.loadFile(ref16sFile, k, threads);
				for(DDLRecord rec : r16){
					if(rec.taxID>0 && map16!=null){
						byte[] seq=map16.get(rec.taxID);
						if(seq!=null){rec.r16S=seq;}
					}
				}
				refs.addAll(r16);
			}
			if(ref18sFile!=null){
				ArrayList<DDLRecord> r18=DDLLoaderMT.loadFile(ref18sFile, k, threads);
				for(DDLRecord rec : r18){
					if(rec.taxID>0 && map18!=null){
						byte[] seq=map18.get(rec.taxID);
						if(seq!=null){rec.r18S=seq;}
					}
				}
				refs.addAll(r18);
			}
		}else{
			System.err.println("Loading SSU references from "+refFile+"...");
			refs=DDLLoaderMT.loadFile(refFile, k, threads);
			int a16=0, a18=0;
			for(DDLRecord rec : refs){
				if(rec.taxID<=0){continue;}
				boolean is16S=(rec.filename!=null && rec.filename.contains("16S"));
				boolean is18S=(rec.filename!=null && rec.filename.contains("18S"));
				if(is16S && map16!=null){
					byte[] seq=map16.get(rec.taxID);
					if(seq!=null){rec.r16S=seq; a16++;}
				}else if(is18S && map18!=null){
					byte[] seq=map18.get(rec.taxID);
					if(seq!=null){rec.r18S=seq; a18++;}
				}else if(map16!=null){
					byte[] seq=map16.get(rec.taxID);
					if(seq!=null){rec.r16S=seq; a16++;}
					else if(map18!=null){
						seq=map18.get(rec.taxID);
						if(seq!=null){rec.r18S=seq; a18++;}
					}
				}
			}
			System.err.println("Attached "+a16+" 16S and "+a18+" 18S sequences to "+refs.size()+" refs.");
		}

		System.err.println("Total: "+refs.size()+" reference SSU DDLs.");

		index=new DDLIndex(refs.get(0).ddl.buckets);
		index.addAll(refs, threads);
		System.err.println("Built inverted index.");

		consensus16S=loadFirstSequence("?16S_consensus_sequence.fa");
		consensus18S=loadFirstSequence("?18S_consensus_sequence.fa");

		try{pgmThread.join();}catch(InterruptedException e){}
		System.err.println("PGM model loaded.");

		nRefs=refs.size();
		startTime=System.currentTimeMillis();
	}

	/*--------------------------------------------------------------*/
	/*----------------       Server Startup         ----------------*/
	/*--------------------------------------------------------------*/

	private void startServer(){
		for(int attempt=0; attempt<100; attempt++){
			try{
				InetSocketAddress isa=new InetSocketAddress(port);
				httpServer=HttpServer.create(isa, 0);
				break;
			}catch(BindException e){
				System.err.println("Port "+port+" in use, retrying in 2s... (attempt "+(attempt+1)+")");
				try{Thread.sleep(2000);}catch(InterruptedException ie){}
			}catch(IOException e){
				throw new RuntimeException(e);
			}
		}
		if(httpServer==null){
			System.err.println("ERROR: Could not bind to port "+port+" after 100 attempts.");
			System.exit(1);
		}

		int handlerThreads=Tools.max(2, threads);
		httpServer.setExecutor(Executors.newFixedThreadPool(handlerThreads));
		httpServer.createContext("/", new UniversalHandler());
		httpServer.createContext("/kill", new KillHandler());
		httpServer.start();
		System.err.println("SSUServer ready on port "+port+" with "+nRefs+" references, "
			+handlerThreads+" handler threads.");
	}

	/*--------------------------------------------------------------*/
	/*----------------      Universal Handler       ----------------*/
	/*--------------------------------------------------------------*/

	private class UniversalHandler implements HttpHandler {
		@Override
		public void handle(HttpExchange t) throws IOException {
			if("OPTIONS".equalsIgnoreCase(t.getRequestMethod())){
				addCorsHeaders(t);
				t.sendResponseHeaders(204, -1);
				t.close();
				return;
			}
			if("GET".equalsIgnoreCase(t.getRequestMethod())){
				addCorsHeaders(t);
				byte[] usage=usageString().getBytes();
				t.sendResponseHeaders(200, usage.length);
				OutputStream os=t.getResponseBody();
				os.write(usage);
				os.close();
				t.close();
				return;
			}
			if(!hasPermission(t)){return;}

			ByteBuilder request=getRequest(t);
			if(request==null){
				reply(t, "# ERROR: Failed to read request body\n".getBytes(), 400);
				return;
			}
			if(maxSize>0 && request.length()>maxSize){
				String msg="{\"error\":\"Input exceeds maximum size of "+maxSize+" bytes\"}";
				reply(t, msg.getBytes(), 413);
				return;
			}

			queryCount.incrementAndGet();
			long tStart=System.nanoTime();

			try{
				byte[] response=getResponse(request);
				long elapsed=System.nanoTime()-tStart;
				if(verbose){
					System.err.println("Query from "+t.getRemoteAddress()
						+" size="+request.length()+" elapsed="+String.format("%.3f", elapsed*1e-9)+"s");
				}
				reply(t, response, 200);
			}catch(Exception e){
				e.printStackTrace();
				String msg="{\"error\":\"Internal server error: "+e.getMessage()+"\"}";
				reply(t, msg.getBytes(), 500);
			}
		}

		private byte[] getResponse(ByteBuilder request){
			boolean jsonOutput=false;
			boolean callMode=false;
			int reqRecords=maxRecords;
			int reqMinHits=minHits;
			int reqBuffer=buffer;

			String body=request.toString();
			if(body.startsWith("//JSON\n")){
				jsonOutput=true;
				body=body.substring(7);
			}
			if(body.startsWith("//Call\n")){
				callMode=true;
				body=body.substring(7);
			}
			if(body.startsWith("//Status")){
				return statusResponse().getBytes();
			}

			DDLFormatter formatter=makeFormatter(jsonOutput);

			while(body.startsWith("//")){
				int nl=body.indexOf('\n');
				if(nl<0){break;}
				String line=body.substring(2, nl);
				body=body.substring(nl+1);
				String[] split=line.split("=");
				if(split.length!=2){continue;}
				String a=split[0].toLowerCase();
				String b=split[1];
				if(formatter.parse(line, a, b)){/* handled */}
				else if(a.equals("records") || a.equals("maxrecords")){reqRecords=Integer.parseInt(b);}
				else if(a.equals("minhits")){reqMinHits=Integer.parseInt(b);}
				else if(a.equals("buffer")){reqBuffer=Integer.parseInt(b);}
			}

			ArrayList<DDLRecord> queries;
			if(callMode){
				queries=processCallMode(body);
			}else{
				queries=processSSUMode(body);
			}

			if(queries.isEmpty()){
				if(jsonOutput){
					return "{\"query_count\":0,\"ref_count\":0,\"elapsed_ms\":0,\"results\":[]}".getBytes();
				}
				return "# No valid SSU sequences found in input\n".getBytes();
			}

			return compareAndFormat(queries, formatter, jsonOutput, reqRecords, reqMinHits, reqBuffer);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Kill Handler          ----------------*/
	/*--------------------------------------------------------------*/

	private class KillHandler implements HttpHandler {
		@Override
		public void handle(HttpExchange t) throws IOException {
			String path=t.getRequestURI().getPath();
			String code=path.substring(path.lastIndexOf('/')+1);
			if(killCode!=null && killCode.equals(code)){
				reply(t, "SSUServer shutting down.\n".getBytes(), 200);
				System.err.println("Kill code received. Shutting down.");
				httpServer.stop(1);
				System.exit(0);
			}else{
				reply(t, "Invalid kill code.\n".getBytes(), 403);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Query Processing       ----------------*/
	/*--------------------------------------------------------------*/

	private ArrayList<DDLRecord> processSSUMode(String fastaBody){
		QuantumAligner aligner=new QuantumAligner();
		ArrayList<DDLRecord> queries=new ArrayList<>();

		String[] lines=fastaBody.split("\n");
		StringBuilder seqBuilder=new StringBuilder();
		String currentName=null;

		for(String line : lines){
			line=line.trim();
			if(line.startsWith(">")){
				if(currentName!=null && seqBuilder.length()>0){
					DDLRecord rec=makeSSURecord(currentName, seqBuilder.toString(), aligner, "http_query");
					if(rec!=null){queries.add(rec);}
				}
				currentName=line.substring(1).trim();
				seqBuilder.setLength(0);
			}else if(!line.isEmpty()){
				seqBuilder.append(line);
			}
		}
		if(currentName!=null && seqBuilder.length()>0){
			DDLRecord rec=makeSSURecord(currentName, seqBuilder.toString(), aligner, "http_query");
			if(rec!=null){queries.add(rec);}
		}
		return queries;
	}

	private DDLRecord makeSSURecord(String name, String seqStr, QuantumAligner aligner, String filename){
		byte[] bases=seqStr.getBytes();
		if(bases.length<50){return null;}
		float id16=(consensus16S!=null ? aligner.align(bases, consensus16S) : 0);
		float id18=(consensus18S!=null ? aligner.align(bases, consensus18S) : 0);
		boolean is16S=(id16>=id18);

		DynamicDemiLog ddl=DynamicDemiLog.create(buckets, k, 12345L, 0f, true);
		Read r=new Read(bases, null, name, 0);
		ddl.hash(r);
		DDLRecord rec=new DDLRecord(ddl, -1, -1, name);
		rec.bases=bases.length;
		rec.contigs=1;
		rec.cardinality=ddl.cardinality();
		rec.filename=filename;
		rec.contigName=name;
		rec.ssuStart=0;
		rec.ssuStrand=(byte)'+';
		if(is16S){rec.r16S=bases;}else{rec.r18S=bases;}
		return rec;
	}

	private ArrayList<DDLRecord> processCallMode(String fastaBody){
		java.io.File tmpFile=null;
		try{
			tmpFile=java.io.File.createTempFile("ssuserver_call_", ".fa");
			java.io.FileWriter fw=new java.io.FileWriter(tmpFile);
			fw.write(fastaBody);
			fw.close();

			ArrayList<SSURecord> ssus=DDLQueryLoaderSF.findAllSSU(tmpFile.getAbsolutePath(), threads, MIN_GENE_CALL_LENGTH);
			ArrayList<DDLRecord> queries=new ArrayList<>();
			for(SSURecord ssu : ssus){
				String cname=ssu.contigName;
				if(cname!=null && cname.indexOf(' ')>0){cname=cname.substring(0, cname.indexOf(' '));}
				Read r=new Read(ssu.bases, null, cname+":"+ssu.start+(char)ssu.strand, 0);
				DynamicDemiLog ddl=DynamicDemiLog.create(buckets, k, 12345L, 0f, true);
				ddl.hash(r);
				DDLRecord rec=new DDLRecord(ddl, -1, -1, r.id);
				rec.bases=ssu.bases.length;
				rec.contigs=1;
				rec.cardinality=ddl.cardinality();
				rec.filename=ssu.fileName!=null ? new java.io.File(ssu.fileName).getName() : "http_query";
				rec.contigName=cname;
				rec.ssuStart=ssu.start;
				rec.ssuStrand=ssu.strand;
				if(ssu.is16S()){rec.r16S=ssu.bases;}else{rec.r18S=ssu.bases;}
				queries.add(rec);
			}
			return queries;
		}catch(IOException e){
			e.printStackTrace();
			return new ArrayList<>();
		}finally{
			if(tmpFile!=null){tmpFile.delete();}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------    Compare and Format        ----------------*/
	/*--------------------------------------------------------------*/

	private byte[] compareAndFormat(ArrayList<DDLRecord> queries, DDLFormatter formatter, boolean jsonOutput,
			int reqRecords, int reqMinHits, int reqBuffer){
		final int nQueries=queries.size();
		@SuppressWarnings("unchecked")
		final ArrayList<DDLComparison>[] allResults=new ArrayList[nQueries];
		final AtomicLong nextQuery=new AtomicLong(0);
		final AtomicLong totalComparisons=new AtomicLong(0);

		long ts=System.nanoTime();
		int compareThreads=Tools.min(threads, nQueries);
		CompareThread[] workers=new CompareThread[compareThreads];
		for(int wi=0; wi<compareThreads; wi++){
			workers[wi]=new CompareThread(queries, refs, nextQuery, allResults,
				nQueries, nRefs, k, reqRecords, reqBuffer, reqMinHits, true, index,
				totalComparisons, true, false);
			workers[wi].start();
		}
		for(CompareThread w : workers){try{w.join();}catch(InterruptedException e){}}
		long tCompare=System.nanoTime()-ts;

		ByteBuilder bb=new ByteBuilder();
		if(jsonOutput){
			bb.append("{\"query_count\":").append(nQueries);
			bb.append(",\"ref_count\":").append(nRefs);
			bb.append(",\"elapsed_ms\":").append((int)(tCompare/1000000));
			bb.append(",\"results\":");
			formatter.jsonStart(bb);
		}else{
			formatter.header(bb);
		}
		for(int qi=0; qi<nQueries; qi++){
			if(allResults[qi]==null){continue;}
			int rank=0;
			for(DDLComparison c : allResults[qi]){
				if(c.equal<reqMinHits){continue;}
				c.rank=++rank;
				formatter.format(c, bb);
			}
		}
		if(jsonOutput){
			formatter.jsonEnd(bb);
			bb.append('}');
		}

		return bb.toBytes();
	}

	/*--------------------------------------------------------------*/
	/*----------------        HTTP Helpers          ----------------*/
	/*--------------------------------------------------------------*/

	private ByteBuilder getRequest(HttpExchange t){
		try{
			InputStream is=t.getRequestBody();
			ByteBuilder bb=new ByteBuilder(8192);
			byte[] buf=new byte[8192];
			int len;
			while((len=is.read(buf))>0){
				bb.append(buf, 0, len);
			}
			is.close();
			return bb;
		}catch(IOException e){
			return null;
		}
	}

	private void reply(HttpExchange t, byte[] response, int code){
		try{
			addCorsHeaders(t);
			t.sendResponseHeaders(code, response.length);
			OutputStream os=t.getResponseBody();
			os.write(response);
			os.close();
			t.close();
		}catch(IOException e){
			if(verbose){e.printStackTrace();}
		}
	}

	private void addCorsHeaders(HttpExchange t){
		Headers h=t.getResponseHeaders();
		h.add("Access-Control-Allow-Origin", domain!=null ? domain : "*");
		h.add("Access-Control-Allow-Methods", "GET, POST, OPTIONS");
		h.add("Access-Control-Allow-Headers", "Content-Type");
	}

	private boolean hasPermission(HttpExchange t) throws IOException {
		if(addressPrefix==null){return true;}
		String address=t.getRemoteAddress().toString();
		if(address.startsWith(addressPrefix) || address.startsWith("/127.0.0.1")){return true;}
		reply(t, "Access denied.\n".getBytes(), 403);
		return false;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Formatter Setup       ----------------*/
	/*--------------------------------------------------------------*/

	private DDLFormatter makeFormatter(boolean jsonOutput){
		DDLFormatter f=new DDLFormatter();
		if(jsonOutput){f.format=DDLFormatter.FORMAT_JSON;}
		f.useAlignmentANI=true;
		f.printSSU=false;
		f.printBases=false;
		f.printCompleteness=false;
		f.printQueryName=false;
		f.printRank=true;
		f.printType=true;
		f.printQLen=true;
		f.printRLen=true;
		f.printFile=true;
		f.printContig=true;
		f.printStart=true;
		f.printStrand=true;
		return f;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Status              ----------------*/
	/*--------------------------------------------------------------*/

	private String statusResponse(){
		long uptime=(System.currentTimeMillis()-startTime)/1000;
		return "SSUServer status: OK\n"
			+"Port: "+port+"\n"
			+"References: "+nRefs+"\n"
			+"Queries served: "+queryCount.get()+"\n"
			+"Uptime: "+uptime+"s\n"
			+"Memory: "+(Runtime.getRuntime().totalMemory()/(1024*1024))+"MB\n";
	}

	private String usageString(){
		return "SSUServer — SSU (16S/18S) classification server\n"
			+"POST raw FASTA to classify SSU sequences.\n"
			+"Prefix body with //JSON\\n for JSON output.\n"
			+"Prefix body with //Call\\n for gene-calling mode.\n"
			+"Per-request parameters as //key=value\\n prefixes:\n"
			+"  //records=10\\n  //minhits=4\\n  //lineage=t\\n  //rank=t\\n\n"
			+"  All DDLFormatter flags supported (printANI, printWKID, etc.)\n"
			+"Port: "+port+"  References: "+nRefs+"\n";
	}

	/*--------------------------------------------------------------*/
	/*----------------          Helpers             ----------------*/
	/*--------------------------------------------------------------*/

	private static byte[] loadFirstSequence(String path){
		String resolved=dna.Data.findPath(path);
		if(resolved==null){System.err.println("Warning: consensus file not found: "+path); return null;}
		FileFormat ff=FileFormat.testInput(resolved, FileFormat.FASTA, null, true, true);
		ArrayList<Read> reads=StreamerFactory.getReads(-1, false, ff, null, null, null);
		return (reads!=null && !reads.isEmpty()) ? reads.get(0).bases : null;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	/* Configuration */
	private int port=3070;
	private int k=19;
	private int exponent=4;
	private int buckets=128;
	private int maxRecords=5;
	private int minHits=8;
	private int buffer=0;
	private int threads=Shared.threads();
	private long maxSize=100*1024*1024;
	private String killCode=null;
	private String addressPrefix=null;
	private String domain=null;
	private boolean verbose=false;
	private boolean verbose2=false;

	/* Reference file paths */
	private String refFile=null;
	private String ref16sFile=null;
	private String ref18sFile=null;

	/* Preloaded state (immutable after startup) */
	private ArrayList<DDLRecord> refs;
	private DDLIndex index;
	private map.IntObjectMap<byte[]>[] maps;
	private byte[] consensus16S;
	private byte[] consensus18S;
	private int nRefs;

	/* Server state */
	private HttpServer httpServer;
	private long startTime;
	private final AtomicLong queryCount=new AtomicLong(0);

	private static final int MIN_GENE_CALL_LENGTH=800;
}
