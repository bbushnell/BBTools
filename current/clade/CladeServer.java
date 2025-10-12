package clade;

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.net.InetSocketAddress;
import java.util.ArrayList;
import java.util.Date;
import java.util.concurrent.atomic.AtomicLong;

import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;

import fileIO.ReadWrite;
import server.ServerTools;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * HTTP server for taxonomic classification using QuickClade.
 * Loads reference database once and handles multiple client requests.
 *
 * @author Chloe
 * @date September 14, 2025
 */
public class CladeServer {

	/*--------------------------------------------------------------*/
	/*----------------            Startup           ----------------*/
	/*--------------------------------------------------------------*/

	/** Command line entrance */
	public static void main(String[] args) throws Exception {
		Timer t=new Timer();
		CladeServer cs=new CladeServer(args);

		t.stop("Time: ");

		System.err.println("Ready!");

		//Server runs until killed
	}

	/** Constructor */
	public CladeServer(String[] args) throws Exception {

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		//Set shared static variables
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		//Default values
		int port_=3069;
		String killCode_=null;
		boolean allowRemoteFileAccess_=false;
		boolean allowLocalHost_=true;
		String addressPrefix_=null;
		ArrayList<String> ref_=new ArrayList<String>();

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
			}else if(a.equals("port")){
				port_=Integer.parseInt(b);
			}else if(a.equals("kill") || a.equals("killcode")){
				killCode_=b;
			}else if(a.equals("localhost") || a.equals("allowlocalhost")){
				allowLocalHost_=Parse.parseBoolean(b);
			}else if(a.equals("prefix") || a.equals("addressprefix")){
				addressPrefix_=b;
			}else if(a.equals("remotefileaccess") || a.equals("allowremotefileaccess")){
				allowRemoteFileAccess_=Parse.parseBoolean(b);
			}else if(a.equals("ref") || a.equals("reference")){
				Tools.getFileOrFiles(b, ref_, true, false, false, false);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		{//Process parser fields
			in=parser.in1;
		}

		//Adjust final fields
		port=port_;
		killCode=killCode_;
		allowRemoteFileAccess=allowRemoteFileAccess_;
		allowLocalHost=allowLocalHost_;
		addressPrefix=addressPrefix_;

		//Load reference
		if(ref_.isEmpty() && in!=null){ref_.add(in);}
		if(ref_.isEmpty()){
			throw new RuntimeException("No reference specified. Use ref=<file>");
		}

		outstream.println("Loading reference database...");
		Timer refTimer=new Timer();

		//Initialize CladeIndex with reference files
		if(verbose){System.err.println("[" + new Date() + "] Loading reference database from: " + ref_);}
		index=CladeIndex.loadIndex(ref_);
		if(verbose){System.err.println("[" + new Date() + "] Database loaded successfully with " + index.size() + " clades");}

		refTimer.stop();
		outstream.println("Loaded "+index.size()+" reference clades in "+refTimer);

		//Initialize the server
		initializeServer();
		if(verbose){
			System.err.println("[" + new Date() + "] CladeServer initialized on port " + port);
			System.err.println("[" + new Date() + "] Verbose mode: " + verbose + ", Verbose2: " + verbose2);
			System.err.println("[" + new Date() + "] Kill code: " + (killCode != null ? "enabled" : "disabled"));
		}
		outstream.println("Clade server started on port "+port);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Server Setup         ----------------*/
	/*--------------------------------------------------------------*/

	/** Initialize and start the HTTP server */
	private void initializeServer() throws Exception {

		//Try to bind the server to the port; repeat until successful
		for(int i=0; i<1000; i++){
			Exception ee=tryInitialize(2000);
			if(ee==null){
				int handlerThreads=Tools.max(2, Shared.threads());
				server.setExecutor(java.util.concurrent.Executors.newFixedThreadPool(handlerThreads));
				server.start();
				serverStartTime=System.currentTimeMillis();
				return;
			}else if(i>6){
				throw ee;
			}
		}
		KillSwitch.kill("Failed to bind to port "+port);
	}

	/** Try to bind to the port */
	private Exception tryInitialize(int millis){
		InetSocketAddress isa=new InetSocketAddress(port);
		Exception ee=null;
		try {
			server=HttpServer.create(isa, 0);
		} catch (java.net.BindException e) {//Expected
			System.err.println(e);
			System.err.println("\nWaiting "+millis+" ms");
			ee=e;
			try {
				Thread.sleep(millis);
			} catch (InterruptedException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			return ee;
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		//Add handlers
		server.createContext("/", new UniversalHandler());
		server.createContext("/clade", new CladeHandler());
		server.createContext("/kill", new KillHandler());
		server.createContext("/stats", new StatsHandler());
		server.createContext("/favicon.ico", new IconHandler());
		return null;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Handlers             ----------------*/
	/*--------------------------------------------------------------*/

	/** Handles requests for favicon.ico */
	class IconHandler implements HttpHandler {

		@Override
		public void handle(HttpExchange t) throws IOException {
			if(verbose2){System.err.println("Icon handler");}
			iconQueries.incrementAndGet();
			ServerTools.reply("", "text/plain", t, verbose2, 404, true);
		}
	}

	/** Handles both GET (help) and POST (clade processing) at root */
	class UniversalHandler implements HttpHandler {

		@Override
		public void handle(HttpExchange t) throws IOException {
			String method = t.getRequestMethod();
			if(verbose2){System.err.println("Universal handler - method: " + method);}

			if("POST".equals(method)) {
				// Route POST requests to CladeHandler
				new CladeHandler().handle(t);
			} else {
				// Route GET requests to help
				if(verbose2){System.err.println("Help handler");}
				final long startTime=System.nanoTime();
				String response = usage();
				ServerTools.reply(response, "text/plain", t, verbose2, 200, true);
				logQuery(t, null, System.nanoTime()-startTime, null);
			}
		}
	}

	/** Handles queries that fall through other handlers */
	class HelpHandler implements HttpHandler {

		@Override
		public void handle(HttpExchange t) throws IOException {
			if(verbose2){System.err.println("Help handler");}
			final long startTime=System.nanoTime();
			String response = usage();
			ServerTools.reply(response, "text/plain", t, verbose2, 200, true);
			logQuery(t, null, System.nanoTime()-startTime, null);
		}
	}

	/** Handles requests for server statistics */
	class StatsHandler implements HttpHandler {

		@Override
		public void handle(HttpExchange t) throws IOException {
			if(verbose2){System.err.println("Stats handler");}
			final long startTime=System.nanoTime();
			ByteBuilder bb=new ByteBuilder();
			long uptimeMillis=System.currentTimeMillis()-serverStartTime;
			long uptimeSeconds=uptimeMillis/1000;
			bb.append("Server uptime: ").append(uptimeSeconds).append(" seconds").append('\n');
			bb.append("Total queries: ").append(queryCount.get()).append('\n');
			bb.append("Clade queries: ").append(cladeQueries.get()).append('\n');
			bb.append("Icon queries: ").append(iconQueries.get()).append('\n');
			bb.append("Reference clades: ").append(index.size()).append('\n');
			String response = bb.toString();
			ServerTools.reply(response, "text/plain", t, verbose2, 200, true);
			logQuery(t, null, System.nanoTime()-startTime, null);
		}
	}

	/** Handles requests to kill the server */
	class KillHandler implements HttpHandler {

		@Override
		public void handle(HttpExchange t) throws IOException {
			if(verbose2){System.err.println("Kill handler");}

			if(killCode==null){
				String response = "Kill code not enabled.";
				ServerTools.reply(response, "text/plain", t, verbose2, 403, true);
				return;
			}

			String query=t.getRequestURI().toString();
			if(verbose2){System.err.println("query="+query);}

			String[] parts=query.split("/");
			String code=(parts.length>2 ? parts[2] : null);

			if(code!=null && code.equals(killCode)){
				String response = "Shutting down server.";
				ServerTools.reply(response, "text/plain", t, verbose2, 200, true);
				System.exit(0);
			}else{
				String response = "Invalid kill code.";
				ServerTools.reply(response, "text/plain", t, verbose2, 403, true);
			}
		}
	}

	/** Handles clade classification requests */
	class CladeHandler implements HttpHandler {

		@Override
		public void handle(HttpExchange t) throws IOException {
			final long handlerStart = System.nanoTime();
			if(verbose2){
				System.err.println("[" + new Date() + "] CladeHandler.handle() - New request from " + t.getRemoteAddress());
				System.err.println("[" + new Date() + "]   Method: " + t.getRequestMethod());
				System.err.println("[" + new Date() + "]   URI: " + t.getRequestURI());
				System.err.println("[" + new Date() + "]   Headers: " + t.getRequestHeaders().entrySet());
			}
			if(verbose){outstream.println("CladeHandler.handle() called");}
			if(verbose2){outstream.println("Got a clade request.");}
			CladeInstance ci=new CladeInstance(t);
			ci.respond();
			if(verbose2){System.err.println("[" + new Date() + "] CladeHandler completed in " + ((System.nanoTime() - handlerStart) / 1000000) + "ms");}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Helper Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/** Generate usage information */
	private String usage(){
		StringBuilder sb=new StringBuilder();
		sb.append("Usage:\n");
		sb.append("  POST to /clade with clade data in request body\n");
		sb.append("  Parameters in URL: /format=oneline/hits=5/\n\n");
		sb.append("Available parameters:\n");
		sb.append("  format={human|oneline} - Output format\n");
		sb.append("  hits=<int> - Number of results per query\n");
		sb.append("  printqtid={t|f} - Print query TaxID\n");
		sb.append("  banself={t|f} - Ban self-matches\n");
		sb.append("  bandupes={t|f} - Ban duplicate matches\n");
		sb.append("  heap=<int> - Heap size for comparisons\n");
		return sb.toString();
	}

	/** Log query information */
	void logQuery(HttpExchange t, String type, long nanos, String response){
		queryCount.incrementAndGet();

		if(verbose){
			ByteBuilder bb=new ByteBuilder();
			bb.append("[").append(new Date().toString()).append("] ");
			bb.append(t.getRemoteAddress().toString()).append(" ");
			if(type!=null){bb.append(type).append(" ");}
			bb.append(String.format("%.3f", nanos/1000000.0)).append("ms");
			if(response!=null && response.length()<100){
				bb.append(" response: ").append(response.trim());
			}
			outstream.println(bb);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	/** HTTP server instance */
	private HttpServer server;

	/** Reference clade index */
	private final CladeIndex index;

	/** Input file (alternative to ref) */
	private String in=null;

	/*--------------------------------------------------------------*/
	/*----------------        Final Fields          ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary output stream */
	private final PrintStream outstream;

	/** Server port */
	private final int port;

	/** Code required to shut down server */
	private final String killCode;

	/** Allow connections from localhost */
	private final boolean allowLocalHost;

	/** Allow remote file access */
	private final boolean allowRemoteFileAccess;

	/** Required address prefix */
	private final String addressPrefix;

	/** Display verbose output */
	private boolean verbose=false;

	/** Display extra verbose output */
	private boolean verbose2=false;

	/** Server start time */
	private long serverStartTime;

	/** Query counters */
	private final AtomicLong queryCount=new AtomicLong();
	private final AtomicLong cladeQueries=new AtomicLong();
	private final AtomicLong iconQueries=new AtomicLong();

	/** Inner class that handles individual clade requests */
	class CladeInstance {

		CladeInstance(HttpExchange t_){
			t=t_;
			startTime=System.nanoTime();
			// Create connection-specific context
			context=new CladeContext();
		}

		void respond() throws IOException {
			final long respondStart = System.nanoTime();
			if(verbose2){System.err.println("[" + new Date() + "] CladeInstance.respond() ENTRY");}
			if(verbose2){System.err.println("DEBUG: respond() called - new HTTP request " + System.currentTimeMillis());}
			try {
				if(verbose){outstream.println("CladeInstance.respond() called");}

				//Security checks
				String address=t.getRemoteAddress().toString();
				if(addressPrefix!=null && !address.startsWith(addressPrefix)){
					if(!allowLocalHost || !address.startsWith("/127.0.0.1")){
						String response="Access denied from "+address;
						ServerTools.reply(response, "text/plain", t, verbose2, 403, true);
						return;
					}
				}

				//Parse URL parameters
				String query=t.getRequestURI().toString();
				if(verbose){outstream.println("Query: "+query);}
				parseParams(query, context);

				//Read request body
				if(verbose2){System.err.println("[" + new Date() + "] Reading request body...");}
				if(verbose){outstream.println("Reading request body...");}
				InputStream is=t.getRequestBody();
				ByteBuilder bb=new ByteBuilder();
				byte[] buffer=new byte[8192];
				int len;
				int totalBytes = 0;
				while((len=is.read(buffer))!=-1){
					bb.append(buffer, 0, len);
					totalBytes += len;
				}
				is.close();
				if(verbose2){System.err.println("[" + new Date() + "] Read " + totalBytes + " bytes from request body");}
				if(verbose){outstream.println("Read "+bb.length()+" bytes from body");}

			//Parse parameters and clades from request body
			String requestBody = bb.toString();
			if(verbose2){System.err.println("[" + new Date() + "] Request body (first 200 chars): " + requestBody.substring(0, Math.min(200, requestBody.length())));}
			parseRequestBody(requestBody, context);
			if(verbose2){System.err.println("[" + new Date() + "] Context after parsing: format=" + context.format + ", hits=" + context.hits + ", heap=" + context.heap);}

			// Parse clades from request body (standard Clade format from SendClade)
			ArrayList<Clade> clades = parseClades(requestBody);

			if(verbose2){System.err.println("DEBUG: Parsing returned " + (clades != null ? clades.size() : "null") + " clades");}
			if(clades==null || clades.isEmpty()){
				String response="No valid clades found in request";
				if(verbose2){System.err.println("DEBUG: Sending 'no clades' response");}
				ServerTools.reply(response, "text/plain", t, verbose2, 400, true);
				return;
			}
			if(verbose2){System.err.println("DEBUG: Proceeding with " + clades.size() + " clades - starting analysis");}

			//Process clades and generate response
			ByteBuilder response=new ByteBuilder();

			if(verbose2){System.err.println("DEBUG: Starting clade processing loop with " + clades.size() + " clades");}
			int queryNumber = 1;
			for(Clade clade : clades){//TODO - multithread this, at least 4 threads.
				if(verbose2){System.err.println("DEBUG: Processing clade " + queryNumber + ": " + (clade != null ? clade.toString() : "null"));}
				// Use thread-safe findBest method with context-specific hits parameter
				if(verbose2){
					System.err.println("[" + new Date() + "] Searching for clade: " + clade.name + " (taxID=" + clade.taxID + ", bases=" + clade.bases + ")");
					System.err.println("DEBUG: Calling index.findBest() with hits=" + context.hits);
				}
				final long searchStart = System.nanoTime();
				ArrayList<Comparison> results=index.findBest(clade, context.hits);
				final long searchTime = (System.nanoTime() - searchStart) / 1000000;
				if(verbose2){
					System.err.println("[" + new Date() + "] Search completed in " + searchTime + "ms, found " + (results != null ? results.size() : 0) + " results");
					System.err.println("DEBUG: findBest() returned " + (results != null ? results.size() : "null") + " results");
				}
				if(verbose2){System.err.println("DEBUG: Calling formatResults() for query " + queryNumber);}
				formatResults(clade, results, response, context, queryNumber);
				if(verbose2){System.err.println("DEBUG: formatResults() completed for query " + queryNumber);}
				queryNumber++;
			}
			if(verbose2){System.err.println("DEBUG: Clade processing loop completed");}

			//Send response
			String responseStr=response.toString();
			if(verbose2){
				System.err.println("[" + new Date() + "] Sending response: " + responseStr.length() + " bytes");
				if(responseStr.length() < 500) {
					System.err.println("[" + new Date() + "] Response content: " + responseStr);
				}
			}
			ServerTools.reply(responseStr, "text/plain", t, verbose2, 200, true);
			if(verbose){System.err.println("[" + new Date() + "] Response sent successfully, total time: " + ((System.nanoTime() - respondStart) / 1000000) + "ms");}

			cladeQueries.incrementAndGet();
			logQuery(t, "clade", System.nanoTime()-startTime, null);

			} catch (Exception e) {
				if(verbose || verbose2){
					System.err.println("[" + new Date() + "] ERROR in CladeInstance.respond(): " + e.getMessage());
					System.err.println("[" + new Date() + "] Stack trace:");
					e.printStackTrace(System.err);
				}
				outstream.println("ERROR in CladeInstance.respond(): "+e.getMessage());
				e.printStackTrace(outstream);
				String errorResponse = "Internal server error: "+e.getMessage();
				try {
					ServerTools.reply(errorResponse, "text/plain", t, verbose2, 500, true);
				} catch (Exception ioe) {
					outstream.println("ERROR sending error response: "+ioe.getMessage());
				}
			}
		}

		/** Parse parameters from URL into context */
		void parseParams(String query, CladeContext ctx){
			if(verbose2){System.err.println("[" + new Date() + "] Parsing URL parameters: " + query);}
			if(query==null || query.length()<2){return;}

			//Remove leading /clade
			if(query.startsWith("/clade")){
				query=query.substring(6);
			}

			//Parse parameters
			String[] parts=query.split("/");
			for(String part : parts){
				if(part.isEmpty()){continue;}
				String[] kv=part.split("=");
				if(kv.length!=2){continue;}

				String key=kv[0].toLowerCase();
				String value=kv[1];

				if(key.equals("format")){
					if(value.equals("oneline") || value.equals("machine")){
						ctx.format=CladeSearcher.MACHINE;
					}else{
						ctx.format=CladeSearcher.HUMAN;
					}
				}else if(key.equals("hits")){
					ctx.hits=Integer.parseInt(value);
				}else if(key.equals("heap")){
					ctx.heap=Integer.parseInt(value);
				}else if(key.equals("printqtid")){
					ctx.printQTID=Parse.parseBoolean(value);
				}else if(key.equals("banself")){
					ctx.banSelf=Parse.parseBoolean(value);
				}else if(key.equals("bandupes")){
					ctx.banDupes=Parse.parseBoolean(value);
				}
			}
		}

		/** Parse parameters from request body (SendClade format) */
		void parseRequestBody(String requestBody, CladeContext ctx){
			if(verbose2){System.err.println("[" + new Date() + "] Parsing request body parameters");}
			if(requestBody==null || requestBody.isEmpty()){return;}

			// SendClade sends parameters in first line: format=oneline/hits=10/heap=1/
			String[] lines = requestBody.split("\n", 2);
			if(lines.length < 1){return;}

			String paramLine = lines[0];

			//Parse parameters from the parameter line
			String[] parts=paramLine.split("/");
			for(String part : parts){
				if(part.isEmpty()){continue;}
				String[] kv=part.split("=");
				if(kv.length!=2){continue;}

				String key=kv[0].toLowerCase();
				String value=kv[1];

				if(key.equals("format")){
					if(value.equals("oneline") || value.equals("machine")){
						ctx.format=CladeSearcher.MACHINE;
					}else{
						ctx.format=CladeSearcher.HUMAN;
					}
				}else if(key.equals("hits")){
					ctx.hits=Integer.parseInt(value);
				}else if(key.equals("heap")){
					ctx.heap=Integer.parseInt(value);
				}else if(key.equals("printqtid")){
					ctx.printQTID=Parse.parseBoolean(value);
				}else if(key.equals("banself")){
					ctx.banSelf=Parse.parseBoolean(value);
				}else if(key.equals("bandupes")){
					ctx.banDupes=Parse.parseBoolean(value);
				}
			}
		}


		/** Parse clades from request body - supports both standard Clade and PreClade formats */
		ArrayList<Clade> parseClades(String data){
			if(verbose2){
				System.err.println("[" + new Date() + "] parseClades() ENTRY - parsing " + data.length() + " bytes");
				if(verbose2){System.err.println("DEBUG: parseClades() called - entry point " + System.currentTimeMillis());}
			}
			ArrayList<Clade> list=new ArrayList<Clade>();

			// Format detection: check first and second lines to distinguish PreClade from standard formats
			String[] allLines = data.split("\n");
			if(allLines.length == 0) {
				if(verbose2){System.err.println("DEBUG: No lines found in data");}
				return list;
			}

			String firstLine = allLines[0].trim();
			if(firstLine.startsWith("//PreClade")){
				if(verbose2){System.err.println("DEBUG: Detected PreClade format (first line)");}
				return parsePreCladeFormat(allLines);
			} else {
				// Check if PreClade header is on second line (after parameters)
				if(allLines.length > 1 && allLines[1].trim().startsWith("//PreClade")){
					if(verbose2){System.err.println("DEBUG: Detected PreClade format (second line)");}
					// Skip the parameter line and pass the rest
					String[] precladeLines = new String[allLines.length - 1];
					System.arraycopy(allLines, 1, precladeLines, 0, allLines.length - 1);
					return parsePreCladeFormat(precladeLines);
				}

				// Standard Clade format or SendClade format - skip parameter line (first line) and only process clade data
				String[] lines=data.split("\n", 2);
				if(lines.length < 2){
					if(verbose2){System.err.println("DEBUG: No clade data found after parameter line");}
					return list;
				}

				String cladeData = lines[1]; // Skip first line (parameters)
				if(verbose2){System.err.println("DEBUG: parsing clade data: " + cladeData.substring(0, Math.min(100, cladeData.length())) + "...");}
				assert(cladeData != null && cladeData.length() > 0) : "CladeData is null or empty - parsing should not continue";

				String[] cladeLines=cladeData.split("\n");
				if(cladeLines.length == 0) {
					if(verbose2){System.err.println("DEBUG: No clade lines found");}
					return list;
				}

				if(verbose2){System.err.println("DEBUG: Detected standard Clade format");}
				return parseStandardCladeFormat(cladeLines);
			}
		}

		/** Parse standard Clade format (existing logic) */
		private ArrayList<Clade> parseStandardCladeFormat(String[] cladeLines){
			ArrayList<Clade> list=new ArrayList<Clade>();

			//Convert clade data to byte array list
			ArrayList<byte[]> byteLines=new ArrayList<byte[]>();
			for(String line : cladeLines){
				byteLines.add(line.getBytes());
			}

			//Create LineParser for parsing (tab-delimited)
			shared.LineParser1 lp=new shared.LineParser1('\t');

			//Use CladeLoader parsing approach: collect lines from one '#' to the next
			ArrayList<byte[]> currentClade = new ArrayList<byte[]>(20);
			for(int pos=0; pos<byteLines.size(); pos++){
				byte[] line=byteLines.get(pos);
				String lineStr=new String(line).trim();

				if(lineStr.startsWith("#") && !currentClade.isEmpty()){
					//Found new clade header and we have collected lines - process previous clade
					Clade c=Clade.parseClade(currentClade, lp);
					if(c!=null){
						c.finish();
						list.add(c);
					}
					currentClade.clear();
				}
				//Always add current line to collection
				currentClade.add(line);
			}

			//Process final clade if any lines remain
			if(currentClade.size() > 1){
				Clade c=Clade.parseClade(currentClade, lp);
				if(c!=null){
					c.finish();
					list.add(c);
				}
			}

			if(verbose2){System.err.println("DEBUG: parseStandardCladeFormat() parsed " + list.size() + " total clades");}
			return list;
		}

		/** Parse PreClade v2.0 format - Brian's idiot-proof 7-line format */
		private ArrayList<Clade> parsePreCladeFormat(String[] lines){
			ArrayList<Clade> list=new ArrayList<Clade>();

			try {
				// PreClade v2.0 format:
				// Line 0: //PreClade Format 2.0 (only once at beginning)
				// Then for each sequence:
				// Line 1: #
				// Line 2: name
				// Line 3: 1-mers (5 values: A,C,G,T,N)
				// Line 4: 2-mers (16 values)
				// Line 5: 3-mers (64 values)
				// Line 6: 4-mers (256 values)
				// Line 7: 5-mers (1024 values)

				int i = 0;

				// Skip the header line if present
				if(i < lines.length && lines[i].trim().startsWith("//PreClade")){
					i++; // Skip the header
				}

				// Now parse entries that start with #
				while(i < lines.length){
					String line = lines[i].trim();

					// Skip empty lines
					if(line.isEmpty()){
						i++;
						continue;
					}

					// Look for entry separator
					if(line.equals("#")){
						// Start of a PreClade entry - must have exactly 6 more lines
						if(i + 6 >= lines.length){
							if(verbose2){
								System.err.println("ERROR: Incomplete PreClade v2.0 entry - need exactly 6 lines after #");
							}
							break;
						}

						// Parse the 6 lines after #
						String name = lines[i+1].trim();       // Sequence name
						String monomers = lines[i+2].trim();   // 1-mers (5 values)
						String dimers = lines[i+3].trim();     // 2-mers (16 values)
						String trimers = lines[i+4].trim();    // 3-mers (64 values)
						String tetramers = lines[i+5].trim();  // 4-mers (256 values)
						String pentamers = lines[i+6].trim();  // 5-mers (1024 values)

						// Parse this PreClade entry
						Clade c = parsePreCladeV2Entry(name, monomers, dimers, trimers, tetramers, pentamers);
						if(c != null){
							list.add(c);
						}

						// Move to next potential record
						i += 7; // Skip the # and 6 data lines
					} else {
						// Skip non-separator lines (shouldn't happen in valid format)
						i++;
					}
				}

			} catch (Exception e) {
				if(verbose2){
					System.err.println("ERROR parsing PreClade v2.0 format: " + e.getMessage());
					e.printStackTrace();
				}
			}

			if(verbose2){System.err.println("DEBUG: parsePreCladeFormat() parsed " + list.size() + " PreClade v2.0 entries");}
			return list;
		}

		/** Parse a single PreClade v2.0 entry with idiot-proof format */
		private Clade parsePreCladeV2Entry(String name, String monomersLine, String dimersLine,
				String trimersLine, String tetramersLine, String pentamersLine) {
			try {
				// Use new PreClade class for proper canonical mapping and GC compensation
				PreClade preClade = new PreClade(name, monomersLine, dimersLine,
					trimersLine, tetramersLine, pentamersLine);
				return preClade.toClade();

			} catch (Exception e) {
				if(verbose2){
					System.err.println("ERROR parsing PreClade v2.0 entry: " + e.getMessage());
					e.printStackTrace();
				}
				return null;
			}
		}


		/** Format results for output */
		void formatResults(Clade query, ArrayList<Comparison> results, ByteBuilder bb, CladeContext ctx, int queryNumber){
			// Version check to verify we're running the fixed code
			if(verbose2){System.err.println("DEBUG: formatResults() ENTRY - Sept 16 2025 version with null pointer fix active");}
			if(verbose2){System.err.println("DEBUG: formatResults() parameters - query=" + (query != null ? query.toString() : "null") +
				" results=" + (results != null ? results.size() : "null") + " bb=" + (bb != null ? "valid" : "null") +
				" ctx=" + (ctx != null ? ctx.toString() : "null"));}
			//Filter out null comparisons before calculating maxHits
			if(verbose2){System.err.println("DEBUG: formatResults() filtering null comparisons from results.size=" + results.size());}
			ArrayList<Comparison> validResults = new ArrayList<>();
			for(Comparison comp : results) {
				if(comp != null && comp.ref != null) {
					validResults.add(comp);
				}
			}
			if(verbose2){System.err.println("DEBUG: formatResults() after filtering: validResults.size=" + validResults.size());}

			//Limit results to requested hits
			if(verbose2){System.err.println("DEBUG: formatResults() calculating maxHits - ctx.hits=" + ctx.hits + " validResults.size=" + validResults.size());}
			int maxHits=Math.min(ctx.hits, validResults.size());
			if(verbose2){System.err.println("DEBUG: formatResults() maxHits=" + maxHits);}

			if(verbose2){System.err.println("DEBUG: formatResults() checking format - ctx.format=" + ctx.format + " MACHINE=" + CladeSearcher.MACHINE);}
			if(ctx.format==CladeSearcher.MACHINE){
				if(verbose2){System.err.println("DEBUG: formatResults() using MACHINE format - entering loop with maxHits=" + maxHits);}
				//One-line format with #Query markers
				bb.append("#Query").append(queryNumber).append('\n');
				for(int i=0; i<maxHits; i++){
					Comparison comp=validResults.get(i);
					// No need for null checks since validResults only contains valid comparisons
					Clade ref=comp.ref;

					bb.append(query.name).tab();
					bb.append(String.format("%.3f", query.gc)).tab();
					bb.append(query.bases).tab();
					bb.append(query.contigs).tab();
					bb.append(ref.name != null ? ref.name : "Unknown_TaxID_" + ref.taxID).tab();
					bb.append(ref.taxID).tab();
					bb.append(String.format("%.3f", ref.gc)).tab();
					bb.append(ref.bases).tab();
					bb.append(ref.contigs).tab();
					bb.append(ref.level).tab();
					bb.append(String.format("%.3f", comp.gcdif)).tab();
					bb.append(String.format("%.3f", comp.strdif)).tab();
					bb.append(String.format("%.3f", comp.k3dif)).tab();
					bb.append(String.format("%.3f", comp.k4dif)).tab();
					bb.append(String.format("%.3f", comp.k5dif)).tab();
					bb.append(ref.lineage()).nl();
				}
			}else{
				if(verbose2){System.err.println("DEBUG: formatResults() using HUMAN format - entering human-readable section");}
				//Human-readable format
				bb.append("Query: ").append(query.name).nl();
				bb.append("GC: ").append(String.format("%.3f", query.gc)).nl();
				bb.append("Bases: ").append(query.bases).nl();
				bb.append("Contigs: ").append(query.contigs).nl();
				bb.nl();

				if(verbose2){System.err.println("DEBUG: formatResults() human format starting loop with maxHits=" + maxHits);}
				for(int i=0; i<maxHits; i++){
					if(verbose2){System.err.println("DEBUG: formatResults() human format processing hit " + i + "/" + maxHits);}
					Comparison comp=validResults.get(i);
					// No need for null checks since validResults only contains valid comparisons
					Clade ref=comp.ref;
					if(ref.name == null) {
						System.err.println("NOTICE: ref.name is null for index " + i + " taxID=" + ref.taxID + " - using fallback");
					}
					if(verbose2){System.err.println("DEBUG: formatResults() human format hit " + i + " validated - formatting output");}

					bb.append("Hit ").append(i+1).append(":\n");
					bb.append("  Name: ").append(ref.name != null ? ref.name : "Unknown_TaxID_" + ref.taxID).nl();
					bb.append("  TaxID: ").append(ref.taxID).nl();
					bb.append("  Level: ").append(ref.level).nl();
					bb.append("  k5dif: ").append(String.format("%.3f", comp.k5dif)).nl();
					bb.append("  Lineage: ").append(ref.lineage()).nl();
					bb.nl();
				}
			}
			if(verbose2){System.err.println("DEBUG: formatResults() COMPLETED SUCCESSFULLY");}
		}

		/** HTTP exchange */
		private final HttpExchange t;

		/** Request start time */
		private final long startTime;

		/** Connection-specific context containing all parameters */
		private final CladeContext context;
	}

}