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
import shared.LineParser1;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * HTTP server for taxonomic classification using QuickClade.
 * Refactored with clean handler architecture supporting multiple request types.
 *
 * Supports:
 * - Standard Clade format classification
 * - PreClade format classification (privacy-preserving)
 * - FetchClade by taxID/organism name
 * - FetchSSU (16S/18S) by taxID/organism name
 * - CompareSSU alignment
 *
 * @author Chloe
 * @date October 14, 2025
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
			}else if(a.equals("prefix") || a.equals("addressprefix")){
				addressPrefix=b;
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

	/** Try to bind to the port and create handler instances */
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

		//Create permanent handler instances
		cladeClassificationHandler=new CladeClassificationHandler();
		fetchCladeHandler=new FetchCladeHandler();
		fetchSSUHandler=new FetchSSUHandler();
		compareSSUHandler=new CompareSSUHandler();

		//Add handlers to server
		server.createContext("/", new UniversalHandler());
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
			reply(t, "".getBytes(), 404);
		}
	}

	/** Routes requests based on body content */
	class UniversalHandler implements HttpHandler {

		@Override
		public void handle(HttpExchange t) throws IOException {
			String method=t.getRequestMethod();
			if(verbose2){System.err.println("Universal handler - method: " + method);}

			//GET returns usage
			if("GET".equals(method)) {
				final long startTime=System.nanoTime();
				reply(t, usage().getBytes(), 200);
				logQuery(t, "help", System.nanoTime()-startTime, null);
				return;
			}

			//POST - route by request type
			final long startTime=System.nanoTime();
			String address=t.getRemoteAddress().toString();
			if(!hasPermission(t, address)){return;}

			ByteBuilder request=getRequest(t);
			if(request==null){
				reply(t, "Failed to read request".getBytes(), 400);
				return;
			}

			int type=queryType(request);
			byte[] response=null;

			switch(type) {
				case CLADE:
				case PRECLADE:
					response=cladeClassificationHandler.getResponse(request);
					cladeQueries.incrementAndGet();
					break;
				case FETCH_CLADE:
					response=fetchCladeHandler.getResponse(request);
					break;
				case FETCH_SSU:
					response=fetchSSUHandler.getResponse(request);
					break;
				case COMPARE_SSU:
					response=compareSSUHandler.getResponse(request);
					break;
				default:
					response="Invalid request type".getBytes();
			}

			if(response!=null){
				reply(t, response, 200);
			}else{
				reply(t, "Handler returned null response".getBytes(), 500);
			}
			logQuery(t, "request", System.nanoTime()-startTime, null);
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
			reply(t, bb.toBytes(), 200);
			logQuery(t, "stats", System.nanoTime()-startTime, null);
		}
	}

	/** Handles requests to kill the server */
	class KillHandler implements HttpHandler {

		@Override
		public void handle(HttpExchange t) throws IOException {
			if(verbose2){System.err.println("Kill handler");}

			if(killCode==null){
				reply(t, "Kill code not enabled.".getBytes(), 403);
				return;
			}

			String query=t.getRequestURI().toString();
			if(verbose2){System.err.println("query="+query);}

			String[] parts=query.split("/");
			String code=(parts.length>2 ? parts[2] : null);

			if(code!=null && code.equals(killCode)){
				reply(t, "Shutting down server.".getBytes(), 200);
				System.exit(0);
			}else{
				reply(t, "Invalid kill code.".getBytes(), 403);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------    Request Type Handlers     ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Handles Clade and PreClade classification requests.
	 * Parses incoming k-mer signatures and classifies against reference database.
	 */
	class CladeClassificationHandler implements HttpHandler {

		/**
		 * Standard HTTP handler interface - checks permission, gets request, processes, replies.
		 */
		@Override
		public void handle(HttpExchange t) throws IOException {
			String address=t.getRemoteAddress().toString();
			if(!hasPermission(t, address)){return;}

			ByteBuilder request=getRequest(t);
			if(request==null){
				reply(t, "Failed to read request".getBytes(), 400);
				return;
			}

			byte[] response=getResponse(request);
			if(response!=null){
				reply(t, response, 200);
			}else{
				reply(t, "Classification failed".getBytes(), 500);
			}
		}

		/**
		 * Processes Clade or PreClade request and returns classification results.
		 *
		 * @param request ByteBuilder containing request body
		 * @return byte[] classification results, or null on error
		 */
		byte[] getResponse(ByteBuilder request) {
			try {
				//Parse URL parameters into context
				CladeContext context=new CladeContext();
				//TODO: Parse URL parameters if needed

				//Parse request body parameters
				context.parseRequestBody(request.toBytes());

				//Parse clades from request body - byte[] parsing for efficiency
				ArrayList<Clade> clades=parseClades(request.toBytes());

				if(clades==null || clades.isEmpty()){
					return "No valid clades found in request".getBytes();
				}

				//Process clades and generate response
				ByteBuilder response=new ByteBuilder();

				int queryNumber=1;
				for(Clade clade : clades){
					//Use thread-safe findBest method with context-specific hits parameter
					ArrayList<Comparison> results=index.findBest(clade, context.hits);
					formatResults(clade, results, response, context, queryNumber);
					queryNumber++;
				}

				return response.toBytes();

			} catch (Exception e) {
				if(verbose || verbose2){
					System.err.println("[" + new Date() + "] ERROR in CladeClassificationHandler: " + e.getMessage());
					e.printStackTrace(System.err);
				}
				return ("Internal server error: "+e.getMessage()).getBytes();
			}
		}

		/**
		 * Parse clades from request body - supports both standard Clade and PreClade formats.
		 * Uses LineParser1 for efficient byte[] parsing without String conversion.
		 */
		private ArrayList<Clade> parseClades(byte[] data){
			if(verbose2){
				System.err.println("[" + new Date() + "] parseClades() ENTRY - parsing " + data.length + " bytes");
			}
			ArrayList<Clade> list=new ArrayList<Clade>();

			//Use LineParser1 to split by newlines - NO String conversion
			LineParser1 newlineParser=new LineParser1('\n');
			newlineParser.set(data);

			if(newlineParser.terms()==0) {
				if(verbose2){System.err.println("DEBUG: No lines found in data");}
				return list;
			}

			//Format detection: check first line for PreClade header
			if(newlineParser.termStartsWith("//PreClade", 0)){
				if(verbose2){System.err.println("DEBUG: Detected PreClade format (first line)");}
				return parsePreCladeFormat(newlineParser);
			} else {
				//Check if PreClade header is on second line (after parameters)
				if(newlineParser.terms()>1 && newlineParser.termStartsWith("//PreClade", 1)){
					if(verbose2){System.err.println("DEBUG: Detected PreClade format (second line)");}
					return parsePreCladeFormat(newlineParser);
				}

				//Standard Clade format - skip first line (parameters), process rest as Clade data
				if(verbose2){System.err.println("DEBUG: Detected standard Clade format");}
				return parseStandardCladeFormat(newlineParser);
			}
		}

		/**
		 * Parse standard Clade format - optimized with byte[] parsing.
		 * Collects lines between '#' markers and passes to Clade.parseClade().
		 */
		private ArrayList<Clade> parseStandardCladeFormat(LineParser1 newlineParser){
			ArrayList<Clade> list=new ArrayList<Clade>();

			//Create LineParser for parsing individual clade lines (tab-delimited)
			LineParser1 lp=new LineParser1('\t');

			//Use CladeLoader parsing approach: collect lines from one '#' to the next
			//Skip first term (parameter line)
			ArrayList<byte[]> currentClade=new ArrayList<byte[]>(20);
			for(int term=1; term<newlineParser.terms(); term++){
				byte[] line=newlineParser.parseByteArray(term);

				if(newlineParser.termStartsWith("#", term) && !currentClade.isEmpty()){
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
			if(currentClade.size()>1){
				Clade c=Clade.parseClade(currentClade, lp);
				if(c!=null){
					c.finish();
					list.add(c);
				}
			}

			if(verbose2){System.err.println("DEBUG: parseStandardCladeFormat() parsed " + list.size() + " total clades");}
			return list;
		}

		/**
		 * Parse PreClade v2.0 format - optimized with byte[] and LineParser1.
		 * PreClade format contains 7 lines per sequence (header, name, 5 k-mer count lines).
		 */
		private ArrayList<Clade> parsePreCladeFormat(LineParser1 newlineParser){
			ArrayList<Clade> list=new ArrayList<Clade>();

			try {
				int term=0;

				//Skip the header line if present
				if(term<newlineParser.terms() && newlineParser.termStartsWith("//PreClade", term)){
					term++; //Skip the header
				}

				//Parse entries that start with #
				while(term<newlineParser.terms()){
					//Skip empty lines
					if(newlineParser.length(term)==0){
						term++;
						continue;
					}

					//Look for entry separator
					if(newlineParser.termEquals("#", term)){
						//Start of a PreClade entry - must have exactly 6 more lines
						if(term+6>=newlineParser.terms()){
							if(verbose2){
								System.err.println("ERROR: Incomplete PreClade v2.0 entry - need exactly 6 lines after #");
							}
							break;
						}

						//Extract the 6 byte[] lines after #
						byte[] nameLine=newlineParser.parseByteArray(term+1);     //Sequence name
						byte[] line1=newlineParser.parseByteArray(term+2);        //1-mers (5 values)
						byte[] line2=newlineParser.parseByteArray(term+3);        //2-mers (16 values)
						byte[] line3=newlineParser.parseByteArray(term+4);        //3-mers (64 values)
						byte[] line4=newlineParser.parseByteArray(term+5);        //4-mers (256 values)
						byte[] line5=newlineParser.parseByteArray(term+6);        //5-mers (1024 values)

						//Parse this PreClade entry using byte[] parsing
						Clade c=parsePreCladeV2Entry(nameLine, line1, line2, line3, line4, line5);
						if(c!=null){
							list.add(c);
						}

						//Move to next potential record
						term+=7; //Skip the # and 6 data lines
					} else {
						//Skip non-separator lines (shouldn't happen in valid format)
						term++;
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

		/**
		 * Parse a single PreClade v2.0 entry - optimized with byte[] and LineParser1.
		 * Uses LineParser1.parseLongArray() which automatically detects comma delimiter.
		 */
		private Clade parsePreCladeV2Entry(byte[] nameLine, byte[] line1, byte[] line2,
				byte[] line3, byte[] line4, byte[] line5) {
			try {
				//Parse name from byte[] - only String we need
				String name=new String(nameLine).trim();

				//Use LineParser1 to parse comma-delimited count lines
				LineParser1 lp=new LineParser1(',');

				//Build k-mer counts array [1-5]
				long[][] counts=new long[6][]; //Index 0 unused, 1-5 for k=1 through k=5

				//Parse each k-mer line using LineParser1.parseLongArray()
				lp.set(line1);
				counts[1]=lp.parseLongArray(0);  //1-mers (5 values)

				lp.set(line2);
				counts[2]=lp.parseLongArray(0);  //2-mers (16 values)

				lp.set(line3);
				counts[3]=lp.parseLongArray(0);  //3-mers (64 values)

				lp.set(line4);
				counts[4]=lp.parseLongArray(0);  //4-mers (256 values)

				lp.set(line5);
				counts[5]=lp.parseLongArray(0);  //5-mers (1024 values)

				//Use byte[]-optimized PreClade constructor
				PreClade preClade=new PreClade(name, counts);
				return preClade.toClade();

			} catch (Exception e) {
				if(verbose2){
					System.err.println("ERROR parsing PreClade v2.0 entry: " + e.getMessage());
					e.printStackTrace();
				}
				return null;
			}
		}

		/**
		 * Format classification results for output.
		 * Supports both HUMAN and MACHINE formats.
		 */
		void formatResults(Clade query, ArrayList<Comparison> results, ByteBuilder bb, CladeContext ctx, int queryNumber){
			//Filter out null comparisons before calculating maxHits
			ArrayList<Comparison> validResults=new ArrayList<>();
			for(Comparison comp : results) {
				if(comp!=null && comp.ref!=null) {
					validResults.add(comp);
				}
			}

			//Limit results to requested hits
			int maxHits=Math.min(ctx.hits, validResults.size());

			if(ctx.format==CladeSearcher.MACHINE){
				//One-line format with #Query markers
				bb.append("#Query").append(queryNumber).append('\n');
				for(int i=0; i<maxHits; i++){
					Comparison comp=validResults.get(i);
					Clade ref=comp.ref;

					bb.append(query.name).tab();
					bb.append(String.format("%.3f", query.gc)).tab();
					bb.append(query.bases).tab();
					bb.append(query.contigs).tab();
					bb.append(ref.name!=null ? ref.name : "Unknown_TaxID_" + ref.taxID).tab();
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
				//Human-readable format
				bb.append("Query: ").append(query.name).nl();
				bb.append("GC: ").append(String.format("%.3f", query.gc)).nl();
				bb.append("Bases: ").append(query.bases).nl();
				bb.append("Contigs: ").append(query.contigs).nl();
				bb.nl();

				for(int i=0; i<maxHits; i++){
					Comparison comp=validResults.get(i);
					Clade ref=comp.ref;

					bb.append("Hit ").append(i+1).append(":\n");
					bb.append("  Name: ").append(ref.name!=null ? ref.name : "Unknown_TaxID_" + ref.taxID).nl();
					bb.append("  TaxID: ").append(ref.taxID).nl();
					bb.append("  Level: ").append(ref.level).nl();
					bb.append("  k5dif: ").append(String.format("%.3f", comp.k5dif)).nl();
					bb.append("  Lineage: ").append(ref.lineage()).nl();
					bb.nl();
				}
			}
		}
	}

	/**
	 * Handles FetchClade requests - fetches Clade file by taxID or organism name.
	 * Format: //FetchClade\nE.coli\nH.sapiens\n562
	 * Returns Clade file data for requested organisms.
	 */
	class FetchCladeHandler implements HttpHandler {

		@Override
		public void handle(HttpExchange t) throws IOException {
			String address=t.getRemoteAddress().toString();
			if(!hasPermission(t, address)){return;}

			ByteBuilder request=getRequest(t);
			byte[] response=getResponse(request);
			reply(t, response, 200);
		}

		/**
		 * Process FetchClade request and return Clade file data.
		 * TODO: Implementation pending - returns null stub for now.
		 */
		byte[] getResponse(ByteBuilder request) {
			//TODO: Parse taxID/organism names from request
			//TODO: Lookup Clade objects from index
			//TODO: Format as Clade file
			return null; //Stub
		}
	}

	/**
	 * Handles FetchSSU requests - fetches 16S/18S sequences by taxID or organism name.
	 * Format: //FetchSSU\nE.coli\nH.sapiens\n562
	 * Returns FASTA format SSU sequences.
	 */
	class FetchSSUHandler implements HttpHandler {

		@Override
		public void handle(HttpExchange t) throws IOException {
			String address=t.getRemoteAddress().toString();
			if(!hasPermission(t, address)){return;}

			ByteBuilder request=getRequest(t);
			byte[] response=getResponse(request);
			reply(t, response, 200);
		}

		/**
		 * Process FetchSSU request and return FASTA sequences.
		 * TODO: Implementation pending - returns null stub for now.
		 */
		byte[] getResponse(ByteBuilder request) {
			//TODO: Parse taxID/organism names from request
			//TODO: Lookup Clade objects from index
			//TODO: Extract r16S or r18S sequences
			//TODO: Format as FASTA
			return null; //Stub
		}
	}

	/**
	 * Handles CompareSSU requests - aligns query SSU to reference and finds best match.
	 * Format: //CompareSSU\n>query\nACGT...
	 * Returns alignment results.
	 */
	class CompareSSUHandler implements HttpHandler {

		@Override
		public void handle(HttpExchange t) throws IOException {
			String address=t.getRemoteAddress().toString();
			if(!hasPermission(t, address)){return;}

			ByteBuilder request=getRequest(t);
			byte[] response=getResponse(request);
			reply(t, response, 200);
		}

		/**
		 * Process CompareSSU request and return alignment results.
		 * TODO: Implementation pending - returns null stub for now.
		 */
		byte[] getResponse(ByteBuilder request) {
			//TODO: Parse query SSU sequence
			//TODO: Align against reference SSUs
			//TODO: Find best matches
			//TODO: Format results
			return null; //Stub
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Helper Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Check if requester has permission to access server.
	 * Sends 403 response if permission denied.
	 */
	private static boolean hasPermission(HttpExchange t, String address) {
		if(addressPrefix!=null && !address.startsWith(addressPrefix)){
			if(!address.startsWith("/127.0.0.1")){
				String response="Access denied from "+address;
				reply(t, response.getBytes(), 403);
				return false;
			}
		}
		return true;
	}

	/**
	 * Read request body from HttpExchange into ByteBuilder.
	 */
	private static final ByteBuilder getRequest(HttpExchange t) {
		//Read request body
		if(verbose2){System.err.println("[" + new Date() + "] Reading request body...");}
		if(verbose){outstream.println("Reading request body...");}
		InputStream is=t.getRequestBody();
		ByteBuilder bb=new ByteBuilder();
		byte[] buffer=new byte[8192];
		int len;
		int totalBytes=0;
		try{
			while((len=is.read(buffer))!=-1){
				bb.append(buffer, 0, len);
				totalBytes+=len;
			}
			is.close();
		}catch(IOException e){
			e.printStackTrace();
			return null;
		}
		if(verbose2){System.err.println("[" + new Date() + "] Read " + totalBytes + " bytes from request body");}
		if(verbose){outstream.println("Read "+bb.length()+" bytes from body");}
		return bb;
	}

	/**
	 * Send reply to client.
	 */
	private static final void reply(HttpExchange t, byte[] response, int code) {
		ServerTools.reply(response, "text/plain", t, verbose2, code, true);
	}

	/**
	 * Detect request type from message body format.
	 * Checks header line for format markers.
	 */
	private static final int queryType(ByteBuilder data) {
		if(data==null || data.length()<1) {return INVALID;}
		if(!data.startsWith("//") || data.startsWith("//Clade")) {return CLADE;}
		if(data.startsWith("//PreClade")) {return PRECLADE;}
		if(data.startsWith("//FetchClade")) {return FETCH_CLADE;}
		if(data.startsWith("//FetchSSU")) {return FETCH_SSU;}
		if(data.startsWith("//CompareSSU")) {return COMPARE_SSU;}
		return INVALID;
	}

	/**
	 * Generate usage/help information.
	 */
	private String usage(){
		StringBuilder sb=new StringBuilder();
		sb.append("CladeServer - Taxonomic Classification Server\n\n");
		sb.append("Supported request types:\n");
		sb.append("  Standard Clade format - Taxonomic classification\n");
		sb.append("  //PreClade - Privacy-preserving classification\n");
		sb.append("  //FetchClade - Fetch Clade file by taxID/organism\n");
		sb.append("  //FetchSSU - Fetch 16S/18S sequence by taxID/organism\n");
		sb.append("  //CompareSSU - Align query SSU to references\n\n");
		sb.append("Parameters in URL: /format=oneline/hits=5/\n");
		sb.append("  format={human|oneline} - Output format\n");
		sb.append("  hits=<int> - Number of results per query\n");
		sb.append("  printqtid={t|f} - Print query TaxID\n");
		sb.append("  banself={t|f} - Ban self-matches\n");
		sb.append("  bandupes={t|f} - Ban duplicate matches\n");
		sb.append("  heap=<int> - Heap size for comparisons\n");
		return sb.toString();
	}

	/**
	 * Log query information for monitoring.
	 */
	static void logQuery(HttpExchange t, String type, long nanos, String response){
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

	/** Input file (alternative to ref) */
	private String in=null;

	/** Server start time */
	private long serverStartTime;

	/** Permanent handler instances */
	private CladeClassificationHandler cladeClassificationHandler;
	private FetchCladeHandler fetchCladeHandler;
	private FetchSSUHandler fetchSSUHandler;
	private CompareSSUHandler compareSSUHandler;

	/*--------------------------------------------------------------*/
	/*----------------        Final Fields          ----------------*/
	/*--------------------------------------------------------------*/

	/** Server port */
	private final int port;

	/** Code required to shut down server */
	private final String killCode;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary output stream */
	private static PrintStream outstream;

	/** Display verbose output */
	private static boolean verbose=false;
	/** Display extra verbose output */
	private static boolean verbose2=false;

	/** Required address prefix */
	private static String addressPrefix;

	/** Reference clade index */
	private static CladeIndex index;

	/** Query counters */
	private static final AtomicLong queryCount=new AtomicLong();
	private static final AtomicLong cladeQueries=new AtomicLong();
	private static final AtomicLong iconQueries=new AtomicLong();

	/** Request type constants */
	private static final int INVALID=0, CLADE=1, PRECLADE=2, FETCH_CLADE=3, COMPARE_SSU=4, FETCH_SSU=5;
}
