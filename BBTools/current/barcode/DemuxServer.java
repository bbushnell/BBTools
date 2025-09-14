package barcode;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.net.InetSocketAddress;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicLongArray;

import com.sun.net.httpserver.Headers;
import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;
import com.sun.net.httpserver.HttpsServer;

import barcode.stub.PCRMatrixProb;
import dna.Data;
import fileIO.ReadWrite;
import server.ServerTools;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sketch.Sketch;
import structures.ByteBuilder;
import structures.StringNum;

/**
 * Server for barcode demux queries.
 * @author Brian Bushnell
 * @date June 24, 2024
 *
 */
public class DemuxServer {
	
	/*--------------------------------------------------------------*/
	/*----------------            Startup           ----------------*/
	/*--------------------------------------------------------------*/

	/** Command line entrance */
	public static void main(String[] args) throws Exception {
		Timer t=new Timer();
		@SuppressWarnings("unused")
		DemuxServer ds=new DemuxServer(args);
		
		t.stop("Time: ");
		
		System.err.println("Ready!");
		
		//ts.begin();
	}
	
	/** Constructor */
	private DemuxServer(String[] args) throws Exception {
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_UNPIGZ=true;
		
		int port_=4096; //Demux server
		String killCode_=null;
		boolean allowLocalHost_=false;
		String addressPrefix_="128."; //LBL
		boolean https=false;
		
		PCRMatrix.devMode=true;
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("verbose2")){
				verbose2=Parse.parseBoolean(b);
			}else if(a.equals("html")){
				useHtml=Parse.parseBoolean(b);
			}else if(a.equals("https")){
				https=Parse.parseBoolean(b);
			}else if(a.equals("http")){
				https=!Parse.parseBoolean(b);
			}else if(a.equals("domain")){
				domain=b;
				while(domain!=null && domain.endsWith("/")){
					domain=domain.substring(0, domain.length()-1);
				}
			}else if(a.equals("port")){
				port_=Integer.parseInt(b);
			}else if(a.equals("kill") || a.equals("killcode")){
				killCode_=b;
			}else if(a.equals("printip")){
				printIP=Parse.parseBoolean(b);
			}else if(a.equals("printheaders")){
				printHeaders=Parse.parseBoolean(b);
			}else if(a.equals("countqueries")){
				countQueries=Parse.parseBoolean(b);
			}else if(a.equals("allowlocalhost")){
				allowLocalHost_=Parse.parseBoolean(b);
			}else if(a.equals("addressprefix")){
				addressPrefix_=b;
			}else if(a.equals("matrixthreads")){
				PCRMatrix.matrixThreads=Integer.parseInt(b);
			}else if(a.equals("devmode")){
				PCRMatrix.devMode=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				throw new RuntimeException(arg);
			}
		}
		
		port=port_;
		allowLocalHost=allowLocalHost_;
		addressPrefix=addressPrefix_;
		
		ReadWrite.USE_UNPIGZ=false;
//		ReadWrite.USE_UNBGZIP=false;
		
		//Wait for server initialization
		httpServer=initializeServer(1000, 8, https);
		assert(httpServer!=null);
		
		//Initialize handlers
		httpServer.createContext("/", new BaseHandler());
		httpServer.createContext("/demux", new DemuxHandler());
		if(killCode!=null){
			httpServer.createContext("/kill", new KillHandler());
		}

		httpServer.createContext("/stats", new StatsHandler());
		httpServer.createContext("/favicon.ico", new IconHandler());
		
		handlerThreads=handlerThreads>0 ? handlerThreads : Tools.max(2, Shared.threads());
		httpServer.setExecutor(java.util.concurrent.Executors.newFixedThreadPool(handlerThreads)); // Creates a multithreaded executor
//		httpServer.setExecutor(java.util.concurrent.Executors.newCachedThreadPool()); // Creates a multithreaded executor
//		httpServer.setExecutor(null); // Creates a singlethreaded executor
		
		//Start the server
		httpServer.start();
	}
	

	
	/** Iterative wait for server initialization */
	private HttpServer initializeServer(int millis0, int iterations, boolean https){
		HttpServer server=null;
		InetSocketAddress isa=new InetSocketAddress(port);
		Exception ee=null;
		for(int i=0, millis=millis0; i<iterations && server==null; i++){
			try {
				if(https){
					server=HttpsServer.create(isa, 0);
				}else{
					server=HttpServer.create(isa, 0);
				}
			} catch (java.net.BindException e) {//Expected
				System.err.println(e);
				System.err.println("\nWaiting "+millis+" ms");
				ee=e;
				ServerTools.pause(millis);
				millis=millis*2;
			} catch (IOException e) {//Not sure when this would occur...  it would be unexpected
				System.err.println(e);
				System.err.println("\nWaiting "+millis+" ms");
				ee=e;
				ServerTools.pause(millis);
				millis=millis*2;
			}
		}
		if(server==null){throw new RuntimeException(ee);}
		return server;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Handlers           ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private class BaseHandler implements HttpHandler {
		
		@SuppressWarnings("restriction")
		@Override
		public void handle(HttpExchange t) throws IOException {
			if(verbose2){
				System.err.println("Base handler");
				printQuery(t, addressPrefix, allowLocalHost, printIP, printHeaders);
				InetSocketAddress remote=t.getRemoteAddress();
				System.err.println("Ping from "+remote);
			}

			baseQueries.incrementAndGet();
			ServerTools.reply("Shush", "text/plain", t, verbose2, 444, true);
		}
		
	}
	
	/*--------------------------------------------------------------*/
	
	/** Handles queries for favicon.ico */
	private class IconHandler implements HttpHandler {
		
		@Override
		public void handle(HttpExchange t) throws IOException {
			if(verbose2){
				System.err.println("Icon handler");
				printQuery(t, addressPrefix, allowLocalHost, printIP, printHeaders);
			}
			iconQueries.incrementAndGet();
			ServerTools.reply(favIcon, "image/x-icon", t, verbose2, 200, true);
		}
		
	}
	
	/*--------------------------------------------------------------*/
	
	/** Handles queries that fall through other handlers */
	private class StatsHandler implements HttpHandler {
		
		@Override
		public void handle(HttpExchange t) throws IOException {
			if(verbose2){
				System.err.println("Stats handler");
				printQuery(t, addressPrefix, allowLocalHost, printIP, printHeaders);
			}
			final long startTime=System.nanoTime();
			returnStats(startTime, t);
		}
		
	}
	
	/*--------------------------------------------------------------*/
	
	/** Handles requests to kill the server */
	private class KillHandler implements HttpHandler {
		
		@Override
		public void handle(HttpExchange t) throws IOException {
			if(verbose2){
				System.err.println("Kill handler");
				printQuery(t, addressPrefix, allowLocalHost, printIP, printHeaders);
			}
			
			//Parse the query from the URL
			String rparam=getRParam(t);
			InetSocketAddress remote=t.getRemoteAddress();
			
			if(testCode(t, rparam)){
				ServerTools.reply("Success.", "text/plain", t, verbose2, 200, true);
				System.err.println("Killed by remote address "+remote);
				//TODO: Perhaps try to close open resources such as the server
				KillSwitch.killSilent();
			}
			
			if(verbose){System.err.println("Bad kill from address "+remote);}
			ServerTools.reply("Bad Code", "text/plain", t, verbose2, 403, true);
		}
		
		/** Determines whether kill code was correct */
		private boolean testCode(HttpExchange t, String rparam){
			String[] params = rparam.split("/");
			if(verbose2){System.err.println(Arrays.toString(params));}
			
			if(killCode!=null){
				if(params.length>1){//URL mode
					return (params[1].equals(killCode));
				}else{//Body mode
					try {
						String code=ServerTools.receive(t);
						return (code!=null && code.equals(killCode));
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			return false;
		}
	}
	
	/** Parse the query from the URL */
	private static String getRParam(HttpExchange t){
		if(verbose2){System.err.println("getRParam");}
		String rparam = t.getRequestURI().toString();
		
		//Trim leading slashes
		while(rparam.startsWith("/")){
			rparam = rparam.substring(1);
		}
		
		//Trim trailing slashes
		while(rparam.endsWith("/")){
			rparam = rparam.substring(0, rparam.length()-1);
		}
		
		if(verbose){System.err.println(rparam==null || rparam.trim().length()<1 ? "usage" : rparam+"\t"+System.currentTimeMillis());}
		return rparam;
	}
	
	/*--------------------------------------------------------------*/

	/** Handles demux queries */
	private class DemuxHandler implements HttpHandler {
		
		DemuxHandler(){}
		
		@Override
		public void handle(HttpExchange t) throws IOException {
			final long startTime=System.nanoTime();
			if(verbose2 || true){
				System.err.println("Demux handler");
				printQuery(t, addressPrefix, allowLocalHost, printIP, printHeaders);
			}
			queries.incrementAndGet();
			
			//Parse the query from the URL
			ArrayList<byte[]> body=getBody(t);
			final long fetchTime=System.nanoTime();
			System.err.println("Fetched "+body.size()+" chunks in "+
					(fetchTime-startTime)/1000000+" ms.");
			DemuxData dd=new DemuxData(body);
			final long parseTime=System.nanoTime();
			System.err.println("Parsed "+body.size()+" chunks in "+
					(parseTime-fetchTime)/1000000+" ms.");
			HashMap<String, String> map=toAssignmentMap(dd);
			ArrayList<Pair> list=toPairList(map);
			map=null;
			final long assignTime=System.nanoTime();
			System.err.println("Assigned "+list.size()+" barcodes in "+
					(assignTime-parseTime)/1000000+" ms.");
			
			reply(list, dd.coding, dd.length1, dd.length2, "text/plain", t, verbose2, 200, true);
			final long replyTime=System.nanoTime();
			System.err.println("Replied with "+list.size()+" barcodes in "+
					(replyTime-assignTime)/1000000+" ms.");
			
			final long stopTime=System.nanoTime();
			final long elapsed=stopTime-startTime;
//			if(response.startsWith("Welcome to ")){
//				timeMeasurementsUsage.incrementAndGet();
//				elapsedTimeUsage.addAndGet(elapsed);
//				lastTimeUsage.set(elapsed);
//			}
//			else{
				timeMeasurementsRemote.incrementAndGet();
				elapsedTimeRemote.addAndGet(elapsed);
				lastTimeRemote.set(elapsed);
//			}
		}
	}
	
//	/** Write to the body of an incoming HTTP session */
//	@SuppressWarnings("restriction")
//	private static boolean reply(HashMap<String, String> map, String type, HttpExchange t, 
//			boolean verbose, int code, boolean close){
//		if(verbose){System.err.println("Sending: "+map.size());}
//		
//		{
//			Headers h = t.getResponseHeaders();
////			String type="text/plain";
//			h.add("Content-Type", type);
//		}
//		try {
//			t.sendResponseHeaders(code, 0);
//			OutputStream os = t.getResponseBody();
//			ByteBuilder bb=new ByteBuilder();
//			for(Entry<String, String> e : map.entrySet()) {
//				bb.append(e.getKey()).tab().append(e.getValue()).nl();
//				if(bb.length>8000000) {
//					os.write(bb.array, 0, bb.length);
//					bb.clear();
//				}
//			}
//			if(bb.length>0) {
//				os.write(bb.array, 0, bb.length);
//				bb.clear();
//			}
//			os.close();
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//			if(close){t.close();}
//			return false;
//		}
//		if(close){t.close();}
//		return true;
//	}
	
	/** Write to the body of an incoming HTTP session */
	@SuppressWarnings("restriction")
	private static boolean reply(ArrayList<Pair> map, int encoding, int len1, int len2,
			String type, HttpExchange t, boolean verbose, int code, boolean close){
		if(verbose){System.err.println("Sending: "+map.size());}
		
		{
			Headers h = t.getResponseHeaders();
//			String type="text/plain";
			h.add("Content-Type", type);
		}
		try {
			t.sendResponseHeaders(code, 0);
			OutputStream os=t.getResponseBody();
			
			writeToStream(map, os, encoding, len1, len2);
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			if(close){t.close();}
			return false;
		}
		if(close){t.close();}
		return true;
	}
	
	public static boolean writeToStream(ArrayList<Pair> map, OutputStream os, final int encoding, int len1, int len2) {
		final byte[] buffer=new byte[64];
		ByteBuilder bb=new ByteBuilder();
		Pair prev=new Pair("", "");
		try {
			for(Pair p : map) {
				if(encoding==Sketch.RAW) {
					bb.append(p.a).tab().append(p.b).nl();
				}else {//Allows map compression
					assert(encoding==Sketch.A48);
					appendA48(bb, p, prev, len1, len2, buffer);
					prev=p;
				}
				if(bb.length>8000000) {
					os.write(bb.array, 0, bb.length);
					bb.clear();
				}
			}
			if(bb.length>0) {
				os.write(bb.array, 0, bb.length);
				bb.clear();
			}
			if(os!=System.out && os!=System.err) {os.close();}//For testing
			return true;
		} catch (IOException e) {
			e.printStackTrace();
		}
		try {
			if(os!=System.out && os!=System.err) {os.close();}
			return false;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return false;
		}
	}
	
	private static void appendA48(ByteBuilder bb, Pair p, Pair prev, int length1, int length2, byte[] buffer) {
		final String a1=p.a, b1=p.b, a0=prev.a, b0=prev.b;
		
		//Key
		if(!Tools.equalsSubstring(a1, a0, 0, length1)) {
			final long x=DemuxData.encodeACGTN(a1, 0, length1);
			DemuxData.appendA48(x, bb, buffer);
		}
		if(length2>0) {
			bb.tab();
			if(!Tools.equalsSubstring(a1, a0, length1+1, length1+1+length2)) {
				final long x=DemuxData.encodeACGTN(a1, length1+1, length1+1+length2);
				DemuxData.appendA48(x, bb, buffer);
			}
		}
		bb.tab();
		
		//Value
		if(!Tools.equalsSubstring(b1, b0, 0, length1)) {
			final long x=DemuxData.encodeACGTN(b1, 0, length1);
			DemuxData.appendA48(x, bb, buffer);
		}
		if(length2>0) {
			bb.tab();
			if(!Tools.equalsSubstring(b1, b0, length1+1, length1+1+length2)) {
				final long x=DemuxData.encodeACGTN(b1, length1+1, length1+1+length2);
				DemuxData.appendA48(x, bb, buffer);
			}
		}
		
		//TODO
//		TODO: I need to decode this on the other end.
		bb.nl();
	}
	
	private HashMap<String, String> toAssignmentMap(DemuxData dd){
		System.err.println("dd: type="+dd.type+", len1="+dd.length1+", len2="+dd.length2+", delimiter="+dd.barcodeDelimiter+
				", hdistsum="+dd.hdistSum+", expected="+dd.expectedList.size()+", counts="+dd.codeCounts.size());
		PCRMatrixProb matrix=(PCRMatrixProb)PCRMatrix.create(dd.type, dd.length1, dd.length2, dd.barcodeDelimiter, dd.hdistSum);
		matrix.populateExpected(dd.expectedList);
		matrix.populateSplitCodes(dd);
		matrix.initializeData();
		matrix.refine(dd.codeCounts, 4);
		HashMap<String, String> map=matrix.makeAssignmentMap(dd.codeCounts, 4);
		return map;
	}
	
	static final ArrayList<Pair> toPairList(HashMap<String, String> map){
		ArrayList<Pair> list=new ArrayList<Pair>(map.size());
		for(Entry<String, String> e : map.entrySet()) {
			Pair p=new Pair(e.getKey(), e.getValue());
			list.add(p);
		}
		Collections.sort(list);
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Helpers           ----------------*/
	/*--------------------------------------------------------------*/
	
	static ArrayList<byte[]> getBody(HttpExchange t){
		if(verbose2){System.err.println("getBody");}
		InputStream is=t.getRequestBody();
		
		ArrayList<byte[]> list=ServerTools.readStreamToList(is);
		return list;
	}
	
	private void returnStats(long startTime, HttpExchange t){
		if(logUsage){System.err.println("stats");}
		String stats=basicStats();
		bytesOut.addAndGet(stats.length());
		ServerTools.reply(stats, "text/plain", t, verbose2, 200, true);
		final long stopTime=System.nanoTime();
		final long elapsed=stopTime-startTime;
		timeMeasurementsUsage.incrementAndGet();
		elapsedTimeUsage.addAndGet(elapsed);
		lastTimeUsage.set(elapsed);
	}
	
	private String basicStats(){
		if(!countQueries){return "";}
		StringBuilder sb=new StringBuilder(500);
		
		final long mq=malformedQueries.get();
		final long iq=internalQueries.get();
		final long q=queries.get();
		final double avgTimeDL=.000001*(elapsedTimeLocal.get()/(Tools.max(1.0, timeMeasurementsLocal.get())));//in milliseconds
		final double lastTimeDL=.000001*lastTimeLocal.get();
		final double avgTimeDR=.000001*(elapsedTimeRemote.get()/(Tools.max(1.0, timeMeasurementsRemote.get())));//in milliseconds
		final double lastTimeDR=.000001*lastTimeRemote.get();
		final double avgTimeDU=.000001*(elapsedTimeUsage.get()/(Tools.max(1.0, timeMeasurementsUsage.get())));//in milliseconds
		final double lastTimeDU=.000001*lastTimeUsage.get();
		final long exq=q-iq;

		sb.append('\n').append("Queries:   ").append(q);
		sb.append('\n').append("Avg time:  ").append(Tools.format("%.3f ms", avgTimeDR));
		sb.append('\n').append("Last time: ").append(Tools.format("%.3f ms", lastTimeDR));
		sb.append('\n').append("Avg time:  ").append(Tools.format("%.3f ms (usage queries)", avgTimeDU));
		sb.append('\n').append("Last time: ").append(Tools.format("%.3f ms (usage queries)", lastTimeDU));

		sb.append('\n');
		sb.append('\n').append("Internal:  ").append(iq);
		sb.append('\n').append("External:  ").append(exq);
		sb.append('\n');
		sb.append('\n').append("BytesIn:   ").append(bytesIn.get());
		sb.append('\n').append("BytesOut:  ").append(bytesOut.get());
		sb.append('\n');
		
		return sb.toString();
	}
	
	private static boolean printQuery(HttpExchange t, 
			String prefix, boolean allowLocalHost, boolean printIP, boolean printHeaders){
		
		InetSocketAddress client=t.getRemoteAddress();
		InetSocketAddress server=t.getLocalAddress();
		
		Headers requestHeaders=t.getRequestHeaders();
		final String xff=requestHeaders.getFirst("X-forwarded-for");
		
		if(printIP){
			final String country=requestHeaders.getFirst("Cf-ipcountry");
			final String agent=requestHeaders.getFirst("User-agent");
			System.err.println("IP: "+client+" "+server);
			System.err.println("Xff="+xff);
			System.err.println("Ipc="+country);
			System.err.println("Agent="+agent);
		}
		
		//This is for IPv4, class A.  Probably extends outside of Berkeley.
		String clientAddress=client.toString();
		String serverAddress=server.toString();
		
		if(clientAddress.contains("127.0.0.1") || 
				clientAddress.contains("/0:0:0:0:0:0:0:1:")){//TODO: contains versus startsWith?
			
			if(printHeaders){
				Headers responseHeaders=t.getResponseHeaders();
				System.err.println("\nRequest: ");
				for(Entry<String, List<String>> entry : requestHeaders.entrySet()){
					System.err.println(entry.getKey()+" -> "+entry.getValue());
				}
				System.err.println("\nResponse: ");
				for(Entry<String, List<String>> entry : responseHeaders.entrySet()){
					System.err.println(entry.getKey()+" -> "+entry.getValue());
				}
			}
			
//			final String xff=requestHeaders.getFirst("X-forwarded-for");
			if(xff!=null){
				if(xff.startsWith(prefix)){
					System.err.println("@Local");
					return true;
				}
				clientAddress=xff;
			}else{
				System.err.println((allowLocalHost ? "@Local" : "@Unknown"));
				return allowLocalHost;
			}
		}else{
			if(clientAddress.startsWith(prefix)){
				System.err.println("@Local");
				return true;
			}
		}
		
		//Makes sure they match up to the first delimiter
		//TODO: This needs to be reviewed
		for(int i=0, max=Tools.max(clientAddress.length(), serverAddress.length()); i<max; i++){
			char cc=clientAddress.charAt(i), sc=serverAddress.charAt(i);
			if(cc!=sc){break;}
			if(cc=='.'){//IPv4
				System.err.println("@Local");
				return true; 
			}else if(cc==':'){//IPv6; probably depends on how long the mask is
				System.err.println("@Local");
				return true;
			}
		}

		System.err.println("@Remote");
		return false;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	static class Pair implements Comparable<Pair> {
		
		Pair(String a_, String b_){
			a=a_;
			b=b_;
		}
		
		@Override
		public int compareTo(Pair o) {
			int x=b.compareTo(o.b);
			return x!=0 ? x : a.compareTo(o.a);
		}
		
		public final String toString() {return a+","+b;}
		
		final String a, b;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Counters           ----------------*/
	/*--------------------------------------------------------------*/
	
	private HashMap<String, StringNum> versionMap=new HashMap<String, StringNum>();
	private AtomicLongArray timesByCount=new AtomicLongArray(10000);
	private AtomicLongArray queryCounts=new AtomicLongArray(10000);
	
	private AtomicLong queries=new AtomicLong(0);
	/** Same IP address mask */
	private AtomicLong internalQueries=new AtomicLong(0);
	private AtomicLong iconQueries=new AtomicLong(0);
	private AtomicLong baseQueries=new AtomicLong(0);
	
	private AtomicLong bytesIn=new AtomicLong(0);
	private AtomicLong bytesOut=new AtomicLong(0);
	
	private AtomicLong elapsedTimeUsage=new AtomicLong(0);
	private AtomicLong timeMeasurementsUsage=new AtomicLong(0);
	private AtomicLong lastTimeUsage=new AtomicLong(0);
	
	private AtomicLong elapsedTimeRemote=new AtomicLong(0);
	private AtomicLong timeMeasurementsRemote=new AtomicLong(0);
	private AtomicLong lastTimeRemote=new AtomicLong(0);
	
	private AtomicLong elapsedTimeLocal=new AtomicLong(0);
	private AtomicLong timeMeasurementsLocal=new AtomicLong(0);
	private AtomicLong lastTimeLocal=new AtomicLong(0);

	private AtomicLong malformedQueries=new AtomicLong(0);
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	private boolean printIP=true;
	private boolean printHeaders=false;
	private boolean countQueries=true;
	
	final boolean allowLocalHost;
	final String addressPrefix;
	/** Address of current server instance (optional) */
	private String domain=null;

	
	
	private final String startTime=new Date().toString();
	
	/** Listen on this port */
	private final int port;
	int handlerThreads=0;
	
	private final HttpServer httpServer;
	private boolean useHtml=false;
	private String killCode=null;
	
	private final String favIconPath=Data.findPath("?favicon.ico");
	private final byte[] favIcon=ReadWrite.readRaw(favIconPath);
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	private static boolean verbose=false, verbose2=false, logUsage=false;
	/** True if an error was encountered */
	private boolean errorState=false;
	
}
