package tax;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.net.InetSocketAddress;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicLongArray;

import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;
import com.sun.net.httpserver.HttpsServer;

import dna.Data;
import fileIO.ReadWrite;
import json.JsonObject;
import server.PercentEncoding;
import server.ServerTools;
import shared.KillSwitch;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sketch.CompareBuffer;
import sketch.Comparison;
import sketch.DisplayParams;
import sketch.Sketch;
import sketch.SketchMakerMini;
import sketch.SketchObject;
import sketch.SketchResults;
import sketch.SketchSearcher;
import sketch.SketchTool;
import sketch.Whitelist;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;
import structures.StringNum;

/**
 * Server for taxonomy or Sketch queries.
 * @author Shijie Yao, Brian Bushnell
 * @date Dec 13, 2016
 *
 */
public class TaxServer {
	
	/*--------------------------------------------------------------*/
	/*----------------            Startup           ----------------*/
	/*--------------------------------------------------------------*/

	/** Command line entrance */
	public static void main(String[] args) throws Exception {
		Timer t=new Timer();
		@SuppressWarnings("unused")
		TaxServer ts=new TaxServer(args);
		
		t.stop("Time: ");
		
		System.err.println("Ready!");
		
		//ts.begin();
	}
	
	/** Constructor */
	public TaxServer(String[] args) throws Exception {
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_UNPIGZ=true;
		TaxFilter.printNodesAdded=false;
		TaxFilter.REQUIRE_PRESENT=false; //Due to missing entries in TaxDump.
		Read.JUNK_MODE=Read.FIX_JUNK;
		SketchObject.compareSelf=true;
		
		int port_=3068; //Taxonomy server
		String killCode_=null;
		boolean allowRemoteFileAccess_=false;
		boolean allowLocalHost_=false;
		String addressPrefix_="128."; //LBL
		long defaultSketchReads=200000;
		boolean https=false;
		
		int serverNum_=0;
		int serverCount_=1;
		
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
				verbose2=SketchObject.verbose2=Parse.parseBoolean(b);
			}else if(a.equals("html")){
				useHtml=Parse.parseBoolean(b);
			}else if(a.equals("https")){
				https=Parse.parseBoolean(b);
			}else if(a.equals("http")){
				https=!Parse.parseBoolean(b);
			}else if(a.equals("servers") || a.equals("numservers") || a.equals("servercount")){
				serverCount_=Integer.parseInt(b);
				assert(serverCount_>0) : arg;
			}else if(a.equals("servernum")){
				serverNum_=Integer.parseInt(b);
				assert(serverNum_>=0) : arg;
			}else if(a.startsWith("slave") && Tools.isDigit(a.charAt(a.length()-1))){
				int num=Integer.parseInt(a.substring(5));
				if(slaveAddress==null){slaveAddress=new ArrayList<String>(serverCount_);}
				while(slaveAddress.size()<=num){slaveAddress.add(null);}
				slaveAddress.set(num, b);
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				giTableFile=b;
			}else if(a.equals("tree") || a.equals("taxtree")){
				taxTreeFile=b;
			}else if(a.equals("accession")){
				accessionFile=b;
			}else if(a.equals("pattern")){
				patternFile=b;
			}else if(a.equals("size") || a.equals("sizefile")){
				sizeFile=b;
			}else if(a.equalsIgnoreCase("img")){
				imgFile=b;
			}else if(a.equals("domain")){
				domain=b;
				while(domain!=null && domain.endsWith("/")){domain=domain.substring(0, domain.length()-1);}
			}else if(a.equals("port")){
				port_=Integer.parseInt(b);
			}else if(a.equals("kill") || a.equals("killcode")){
				killCode_=b;
			}else if(a.equals("oldcode")){
				oldKillCode=b;
			}else if(a.equals("oldaddress")){
				oldAddress=b;
			}else if(a.equals("sketchonly")){
				sketchOnly=Parse.parseBoolean(b);
			}else if(a.equals("sketchreads")){
				defaultSketchReads=Parse.parseKMG(b);
			}else if(a.equals("handlerthreads")){
				handlerThreads=Integer.parseInt(b);
			}else if(a.equals("sketchthreads") || a.equals("sketchcomparethreads")){
				maxConcurrentSketchCompareThreads=Integer.parseInt(b);
			}else if(a.equals("sketchloadthreads")){
				maxConcurrentSketchLoadThreads=Integer.parseInt(b);
			}else if(a.equals("hashnames")){
				hashNames=Parse.parseBoolean(b);
			}else if(a.equals("hashdotformat")){
				hashDotFormat=Parse.parseBoolean(b);
			}else if(a.equals("printip")){
				printIP=Parse.parseBoolean(b);
			}else if(a.equals("printheaders")){
				printHeaders=Parse.parseBoolean(b);
			}else if(a.equals("countqueries")){
				countQueries=Parse.parseBoolean(b);
			}else if(a.equals("clear") || a.equals("clearmem")){
				clearMem=Parse.parseBoolean(b);
			}else if(a.equals("dbname")){
				SketchObject.defaultParams.dbName=b;
			}else if(a.equals("allowremotefileaccess")){
				allowRemoteFileAccess_=Parse.parseBoolean(b);
			}else if(a.equals("allowlocalhost")){
				allowLocalHost_=Parse.parseBoolean(b);
			}else if(a.equals("addressprefix")){
				addressPrefix_=b;
			}else if(a.equals("maxpigzprocesses")){
				AccessionToTaxid.maxPigzProcesses=Integer.parseInt(b);
			}else if(a.equals("path") || a.equals("treepath") || a.equals("basepath")){
				basePath=b;
			}else if(a.equalsIgnoreCase("prealloc")){
				if(b==null || Character.isLetter(b.charAt(0))){
					if(Parse.parseBoolean(b)){
						//The hashtable for strings is getting huge.
						//Setting this number lower gives more space for hashing long strings,
						//and less space for hashing numbers
						prealloc=0.75f;
					}else{
						prealloc=0;
					}
				}else{
					prealloc=Float.parseFloat(b);
				}
				SketchObject.prealloc=prealloc;
			}else if(searcher.parse(arg, a, b, true)){
				//do nothing
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				throw new RuntimeException(arg);
			}
		}
		if("auto".equalsIgnoreCase(imgFile)){imgFile=TaxTree.defaultImgFile();}
		if("auto".equalsIgnoreCase(taxTreeFile)){taxTreeFile=TaxTree.defaultTreeFile();}
		if("auto".equalsIgnoreCase(giTableFile)){giTableFile=TaxTree.defaultTableFile();}
		if("auto".equalsIgnoreCase(accessionFile)){accessionFile=TaxTree.defaultAccessionFile();}
		if("auto".equalsIgnoreCase(patternFile)){patternFile=TaxTree.defaultPatternFile();}
		if("auto".equalsIgnoreCase(sizeFile)){sizeFile=TaxTree.defaultSizeFile();}
		
		serverNum=AccessionToTaxid.serverNum=serverNum_;
		serverCount=AccessionToTaxid.serverCount=serverCount_;
		distributed=AccessionToTaxid.distributed=serverCount>1;
		assert(serverNum<serverCount && serverNum>=0);
		if(distributed && serverNum==0){
			assert(slaveAddress!=null);
			assert(slaveAddress.size()==serverCount);
			for(int i=1; i<slaveAddress.size(); i++){
				assert(slaveAddress.get(i)!=null);
			}
		}
		
		maxConcurrentSketchCompareThreads=Tools.mid(1, maxConcurrentSketchCompareThreads, Shared.threads());
		maxConcurrentSketchLoadThreads=Tools.mid(1, maxConcurrentSketchLoadThreads, Shared.threads());
		assert(maxConcurrentSketchCompareThreads>=1);
		assert(maxConcurrentSketchLoadThreads>=1);
		
		if(basePath==null || basePath.trim().length()==0){basePath="";}
		else{
			basePath=basePath.trim().replace('\\', '/').replaceAll("/+", "/");
			if(!basePath.endsWith("/")){basePath=basePath+"/";}
		}
		
		//Adjust SketchSearch rcomp and amino flags
		SketchObject.postParse();
		
		if(sketchOnly){
//			hashNames=false;
			giTableFile=null;
			accessionFile=null;
			imgFile=null;
			patternFile=null;
		}
		
		port=port_;
		killCode=killCode_;
		allowRemoteFileAccess=allowRemoteFileAccess_;
		allowLocalHost=allowLocalHost_;
		addressPrefix=addressPrefix_;
		
		//Fill some data objects
		USAGE=makeUsagePrefix();
		rawHtml=(useHtml ? loadRawHtml() : null);
		typeMap=makeTypeMap();
		commonMap=makeCommonMap();
		
		//Load the GI table
		if(giTableFile!=null){
			outstream.println("Loading gi table.");
			GiToTaxid.initialize(giTableFile);
		}
		
		//Load the taxTree
		if(taxTreeFile!=null){
			tree=TaxTree.loadTaxTree(taxTreeFile, outstream, hashNames, hashDotFormat);
			if(hashNames){tree.hashChildren();}
			assert(tree.nameMap!=null || sketchOnly);
		}else{//The tree is required
			tree=null;
			throw new RuntimeException("No tree specified.");
		}
		//Set a default taxtree for sketch-related usage
		SketchObject.taxtree=tree;
		
		if(sizeFile!=null){
			Timer t=new Timer();
			outstream.println("Loading size file.");
			tree.loadSizeFile(sizeFile);
			t.stopAndPrint();
		}
		
		if(imgFile!=null){
			TaxTree.loadIMG(imgFile, false, outstream);
		}
		
		if(patternFile!=null){
			Timer t=new Timer();
			AnalyzeAccession.loadCodeMap(patternFile);
			outstream.println("Loading pattern table.");
			t.stopAndPrint();
		}
		
		//Load accession files
		if(accessionFile!=null){
			Timer t=new Timer();
			AccessionToTaxid.tree=tree;
			AccessionToTaxid.prealloc=prealloc;
			outstream.println("Loading accession table.");
			AccessionToTaxid.load(accessionFile);
			t.stopAndPrint();
//			if(searcher.refFiles.isEmpty()){System.gc();}
		}
		
//		assert(false) : searcher.refFileCount();
		
		//Load reference sketches
		hasSketches=searcher.refFileCount()>0;
		if(hasSketches){
			outstream.println("Loading sketches.");
			Timer t=new Timer();
			searcher.loadReferences(SketchObject.PER_TAXA, SketchObject.defaultParams);
			t.stopAndPrint();
//			System.gc();
		}
		
		SketchObject.allowMultithreadedFastq=(maxConcurrentSketchLoadThreads>1);
		SketchObject.defaultParams.maxReads=defaultSketchReads;
		ReadWrite.USE_UNPIGZ=false;
//		ReadWrite.USE_UNBGZIP=false;
		
		if(clearMem){
			System.err.println("Clearing memory.");
			System.gc();
			Shared.printMemory();
		}
		
		//If there is a kill code, kill the old instance
		if(oldKillCode!=null && oldAddress!=null){
			outstream.println("Killing old instance.");
			killOldInstance();
		}
		
		//Wait for server initialization
		httpServer=initializeServer(1000, 8, https);
		assert(httpServer!=null);
		
		//Initialize handlers
		if(!sketchOnly){
			httpServer.createContext("/", new TaxHandler(false));
			httpServer.createContext("/tax", new TaxHandler(false));
			httpServer.createContext("/stax", new TaxHandler(true));
			httpServer.createContext("/simpletax", new TaxHandler(true));
		}else{
			httpServer.createContext("/", new SketchHandler());
		}
		httpServer.createContext("/sketch", new SketchHandler());
		if(killCode!=null){
			httpServer.createContext("/kill", new KillHandler());
		}

		httpServer.createContext("/help", new HelpHandler());
		httpServer.createContext("/usage", new HelpHandler());
		httpServer.createContext("/stats", new StatsHandler());
		httpServer.createContext("/favicon.ico", new IconHandler());
		
		handlerThreads=handlerThreads>0 ? handlerThreads : Tools.max(2, Shared.threads());
		httpServer.setExecutor(java.util.concurrent.Executors.newFixedThreadPool(handlerThreads)); // Creates a multithreaded executor
//		httpServer.setExecutor(java.util.concurrent.Executors.newCachedThreadPool()); // Creates a multithreaded executor
//		httpServer.setExecutor(null); // Creates a singlethreaded executor
		
		//Start the server
		httpServer.start();
	}
	
	/** Kill a prior server instance */
	private void killOldInstance(){
		StringNum result=null;
		try {
			result=ServerTools.sendAndReceive(oldKillCode.getBytes(), oldAddress);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			if(e!=null){e.printStackTrace();}
			System.err.println("\nException suppressed; continuing.\n");
			return;
		}
		if(result==null || result.s==null || !"Success.".equals(result.s)){
//			KillSwitch.kill("Bad kill result: "+result+"\nQuitting.\n");
			System.err.println("Bad kill result: "+result+"\nContinuing.\n");
		}
		ServerTools.pause(1000);
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
	
	public void returnUsage(long startTime, HttpExchange t){
		if(useHtml){
			returnUsageHtml(startTime, t);
			return;
		}
		if(logUsage){System.err.println("usage");}
//		String usage=USAGE(USAGE);
		bytesOut.addAndGet(USAGE.length());
		ServerTools.reply(USAGE, "text/plain", t, verbose2, 200, true);
		final long stopTime=System.nanoTime();
		final long elapsed=stopTime-startTime;
		timeMeasurementsUsage.incrementAndGet();
		elapsedTimeUsage.addAndGet(elapsed);
		lastTimeUsage.set(elapsed);
	}
	
	public void returnUsageHtml(long startTime, HttpExchange t){
		if(logUsage){System.err.println("usage");}
		String s=makeUsageHtml();
		bytesOut.addAndGet(s.length());
		ServerTools.reply(s, "html", t, verbose2, 200, true);
		final long stopTime=System.nanoTime();
		final long elapsed=stopTime-startTime;
		timeMeasurementsUsage.incrementAndGet();
		elapsedTimeUsage.addAndGet(elapsed);
		lastTimeUsage.set(elapsed);
	}
	
	public void returnStats(long startTime, HttpExchange t){
		if(logUsage){System.err.println("stats");}
		String stats=makeStats();
		bytesOut.addAndGet(stats.length());
		ServerTools.reply(stats, "text/plain", t, verbose2, 200, true);
		final long stopTime=System.nanoTime();
		final long elapsed=stopTime-startTime;
		timeMeasurementsUsage.incrementAndGet();
		elapsedTimeUsage.addAndGet(elapsed);
		lastTimeUsage.set(elapsed);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Handlers           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Handles queries for favicon.ico */
	class IconHandler implements HttpHandler {
		
		@Override
		public void handle(HttpExchange t) throws IOException {
			if(verbose2){System.err.println("Icon handler");}
			iconQueries.incrementAndGet();
			ServerTools.reply(favIcon, "image/x-icon", t, verbose2, 200, true);
		}
		
	}
	
	/*--------------------------------------------------------------*/
	
	/** Handles queries that fall through other handlers */
	class HelpHandler implements HttpHandler {
		
		@Override
		public void handle(HttpExchange t) throws IOException {
			if(verbose2){System.err.println("Help handler");}
			final long startTime=System.nanoTime();
			returnUsage(startTime, t);
		}
		
	}
	
	/*--------------------------------------------------------------*/
	
	/** Handles queries that fall through other handlers */
	class StatsHandler implements HttpHandler {
		
		@Override
		public void handle(HttpExchange t) throws IOException {
			if(verbose2){System.err.println("Http handler");}
			final long startTime=System.nanoTime();
			returnStats(startTime, t);
		}
		
	}
	
	/*--------------------------------------------------------------*/
	
	/** Handles requests to kill the server */
	class KillHandler implements HttpHandler {
		
		@Override
		public void handle(HttpExchange t) throws IOException {
			if(verbose2){System.err.println("Kill handler");}
			
			//Parse the query from the URL
			String rparam=getRParam(t, false);
			InetSocketAddress remote=t.getRemoteAddress();
			
			if(testCode(t, rparam)){
				ServerTools.reply("Success.", "text/plain", t, verbose2, 200, true);
				System.err.println("Killed by remote address "+remote);
				//TODO: Perhaps try to close open resources such as the server
				KillSwitch.killSilent();
			}
			
			if(verbose){System.err.println("Bad kill from address "+remote);}
			ServerTools.reply(BAD_CODE, "text/plain", t, verbose2, 403, true);
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
	
	/*--------------------------------------------------------------*/

	/** Listens for sketch comparison requests */
	class SketchHandler implements HttpHandler {
		
		@Override
		public void handle(HttpExchange t) throws IOException {
			if(verbose2){outstream.println("Got a request.");}
			SketchInstance si=new SketchInstance(t);
			if(verbose2){outstream.println("Made si.");}
			si.handleInner();
			if(verbose2){outstream.println("Done.");}
		}
		
		private ArrayList<Sketch> loadSketchesFromBody(String body){
			//List of query sketches
			ArrayList<Sketch> sketches=null;
			
			if(body!=null && body.length()>0){
				sketches=searcher.loadSketchesFromString(body);
				if(Whitelist.exists()){
					for(Sketch sk : sketches){
						Whitelist.apply(sk);
					}
				}
			}
			return sketches;
		}
		
		private ArrayList<Sketch> loadSketchesFromFile(String fname, DisplayParams params){
			//List of query sketches
			ArrayList<Sketch> sketches=null;
			
			SketchTool tool=searcher.tool;
			if(tool.minKeyOccuranceCount!=params.minKeyOccuranceCount || params.trackCounts()){
				tool=new SketchTool(SketchObject.targetSketchSize, params);
			}
			
			if(verbose2){System.err.println("Loading sketches from file "+fname);}
			sketches=tool.loadSketchesFromFile(fname, (SketchMakerMini)null, maxConcurrentSketchLoadThreads, params.maxReads, params.mode, params, true);
			if(verbose2){System.err.println("Loaded "+(sketches==null ? "null" : sketches.size())+" sketches from file "+fname);}
			return sketches;
		}
		
		//Created in handle()
		private class SketchInstance {
			
			SketchInstance(HttpExchange t_){
				t=t_;
				instanceStartTime=System.nanoTime();
			}
			
			void handleInner(){
				
				if(!hasSketches){
					if(verbose2){System.err.println("No sketches.");}
					ServerTools.reply("\nERROR: This server has no sketches loaded.\n"
							+ "Please download the latest BBTools version to use SendSketch.\n", "text/plain", t, verbose2, 400, true);
					return;
				}
				
				String rparam=parseRparamSketch(t);
				if(verbose2){System.err.println("Parsed rparam.");}
				if(rparam==null){
					returnUsage(instanceStartTime, t);
					return;
				}
				final boolean internal=incrementQueries(t, fileMode, refMode, false, false, false, false, false, false, false, false, -1);
				if(verbose2){System.err.println("Incremented queries rparam.");}
				
				if(verbose2){System.err.println(rparam);}
				if(verbose2){System.err.println("fileMode="+fileMode+", refMode="+refMode);}
				
				if(fileMode && (!internal || !allowRemoteFileAccess)){
					if(verbose){System.err.println("Illegal file query from "+ServerTools.getClientAddress(t));}
					malformedQueries.incrementAndGet();
//					if(verbose){System.err.println("test1");}
//					String body=getBody(t);//123
					ServerTools.reply("\nERROR: This server does not allow remote file access. "
							+ "You may only use the 'local' flag from within the local intranet, for files on the server's filesystem.\n", "text/plain", t, verbose2, 400, true);
//					if(verbose){System.err.println("test2");}
					return;
				}
				
				String body=getBody(t);
				if(verbose2){System.err.println("Got body.");}
				
				if(body!=null){bytesIn.addAndGet(body.length());}
				
				if(verbose2){System.err.println("Found body: "+body);}
				if(body!=null && body.length()>0){
					if((fileMode || refMode) && !body.startsWith("##")){
						body="##"+body;
					}
					try {
						params=params.parseDoubleHeader(body);
						if(verbose2){System.err.println("Passed parse params.");}
						
					} catch (Throwable e) {
						String s=Tools.toString(e);
						ServerTools.reply("\nERROR: \n"+ s,
								"text/plain", t, verbose2, 400, true);
						return;
					}
					if(!params.compatible()){
						ServerTools.reply("\nERROR: The sketch is not compatible with this server.\n"
								+ "Server settings: k="+SketchObject.k+(SketchObject.k2>0 ? ","+SketchObject.k2 : "")
								+" amino="+SketchObject.amino+" hash_version="+SketchObject.HASH_VERSION+"\n"
								+ "You may need to download a newer version of BBTools; this server is running version "+Shared.BBTOOLS_VERSION_STRING,
								"text/plain", t, verbose2, 400, true);
						return;
					}
				}
				if(params.trackCounts()){
					depthQueries.incrementAndGet();
				}
				
				if(verbose2){System.err.println("Parsed params: "+params.toString());}
				
				//List of query sketches
				ArrayList<Sketch> sketches;
				
				if(fileMode){
					assert(internal && allowRemoteFileAccess);
					if(!allowRemoteFileAccess) {throw new RuntimeException("Bad access to "+rparam);}
					File f=new File(rparam);
					if(!f.exists() && !rparam.startsWith("/")){
						String temp="/"+rparam;
						f=new File(temp);
						if(f.exists()){rparam=temp;}
					}
					if(!f.exists() || !f.canRead()) {
						ServerTools.reply("\nERROR: This file cannot be read.\n", "text/plain", t, verbose2, 400, true);
						return;
					}
//					if(f.exists()){
//						if(f.length()>100000000L){
//							if(params.reads<0){
//								params.reads=200000;//Cap default number of reads at 200000
//							}
//						}
//					}
					sketches=loadSketchesFromFile(rparam, params);
				}else if(refMode){
					String[] split=rparam.split(",");
					sketches=new ArrayList<Sketch>(split.length);
					for(String s : split){
						Sketch sk=findRefSketch(s);
						if(sk!=null){sketches.add(sk);}
					}
				}else{
					sketches=loadSketchesFromBody(body);
				}
				if(verbose2){System.err.println("Loaded "+sketches.size()+" sketches.");}
				
				final int numSketches=sketches==null ? 0 : sketches.size();
				if(params.chunkNum<0){
					if(numSketches<2){
						unknownChunkSingle.incrementAndGet();
					}else{
						unknownChunkMulti.incrementAndGet();
					}
				}else if(params.chunkNum==0){
					if(numSketches<2){
						firstChunkSingle.incrementAndGet();
					}else{
						firstChunkMulti.incrementAndGet();
					}
				}else{
					if(numSketches<2){
						nthChunkSingle.incrementAndGet();
					}else{
						nthChunkMulti.incrementAndGet();
					}
				}
				bulkCount.addAndGet(numSketches);
				
				if(params.inputVersion==null){params.inputVersion="unknown";}
				synchronized(versionMap){
					StringNum sn=versionMap.get(params.inputVersion);
					if(sn==null){versionMap.put(params.inputVersion, new StringNum(params.inputVersion, 1));}
					else{sn.increment();}
				}
				
				String response=null;
				if(sketches==null || sketches.isEmpty()){
					malformedQueries.incrementAndGet();
					response="Error.";
					if(verbose){
						StringBuilder sb=new StringBuilder();
						sb.append("Malformed query from ").append(ServerTools.getClientAddress(t)).append(". body:");
						if(body==null){
							sb.append(" null");
						}else{
							String[] split = body.split("\n");
							sb.append(" ").append(split.length).append(" lines total, displaying ").append(Tools.min(3, split.length)).append('.');
							for(int i=0; i<3 && i<split.length; i++){
								String s=split[i];
								int len=s.length();
								if(s.length()>1000){s=s.substring(0, 1000)+" [truncated, "+len+" total]";}
								sb.append('\n');
								sb.append(s);
							}
						}
						System.err.println(sb);
					}
				}else{
					if(verbose2){
						System.err.println("Received "+sketches.get(0).name()+", size "+sketches.get(0).keys.length);
						System.err.println("params: "+params);
						System.err.println("postparsed: "+params.postParsed());
						System.err.println("taxwhitelist: "+params.taxFilterWhite);
					}
					response=compare(sketches, params);
//					searcher.compare(sketches, response, params, maxConcurrentSketchCompareThreads); //This is where it gets stuck if comparing takes too long
					if(verbose2){System.err.println("Result: '"+response+"'");}
				}
				
				bytesOut.addAndGet(response.length());
				ServerTools.reply(response, "text/plain", t, verbose2, 200, true);
				
				final long stopTime=System.nanoTime();
				final long elapsed=stopTime-instanceStartTime;
				if(fileMode){
					timeMeasurementsLocal.incrementAndGet();
					elapsedTimeLocal.addAndGet(elapsed);
					lastTimeLocal.set(elapsed);
				}else if(refMode){
					timeMeasurementsReference.incrementAndGet();
					elapsedTimeReference.addAndGet(elapsed);
					lastTimeReference.set(elapsed);
				}else{
					timeMeasurementsRemote.incrementAndGet();
					elapsedTimeRemote.addAndGet(elapsed);
					lastTimeRemote.set(elapsed);

					queryCounts.incrementAndGet(Tools.min(numSketches, queryCounts.length()-1));
					timesByCount.addAndGet(Tools.min(numSketches, queryCounts.length()-1), elapsed);
				}
			}
			
			private String parseRparamSketch(HttpExchange t){
				//Parse the query from the URL
				String rparam=getRParam(t, false);
				if(rparam!=null){bytesIn.addAndGet(rparam.length());}

				if(rparam.length()<1 || rparam.equalsIgnoreCase("help") || rparam.equalsIgnoreCase("usage") || rparam.equalsIgnoreCase("help/") || rparam.equalsIgnoreCase("usage/")){
					return null;
				}
				
				if(rparam.startsWith("sketch/")){rparam=rparam.substring(7);}
				else if(rparam.equals("sketch")){rparam="";}
				while(rparam.startsWith("/")){rparam=rparam.substring(1);}
				
				//Toggle between local files and sketch transmission
				
				if(rparam.length()<2){
					params=SketchObject.defaultParams;
				}else{
					params=SketchObject.defaultParams.clone();
					String[] args=rparam.split("/");
					int trimmed=0;
					for(int i=0; i<args.length; i++){//parse rparam
						String arg=args[i];
						if(arg.length()>0){
							String[] split=arg.split("=");
							String a=split[0].toLowerCase();
							String b=split.length>1 ? split[1] : null;

							if(a.equals("file")){
								fileMode=true;
								trimmed+=5;
								break;
							}else if(a.equals("ref") || a.equals("taxid") || a.equals("tid")){
								refMode=true;
								trimmed+=4;
								break;
							}else if(params.parse(arg, a, b)){
								trimmed+=arg.length()+1;
							}else{
								assert(false) : "Bad argument:'"+arg+"'"+"\n"+Arrays.toString(args)+"\n"+rparam;
							}
						}
					}
					params.postParse(true, true);
//					System.err.println("Trimmed="+trimmed+", rparam="+rparam);
					if(trimmed>0){
						rparam=rparam.substring(Tools.min(trimmed, rparam.length()));
					}
//					System.err.println("rparam="+rparam);
				}
				
				if(verbose2){
					System.err.println(rparam);
					System.err.println("rparam.startsWith(\"file/\"):"+rparam.startsWith("file/"));
				}
				
				return rparam;
			}
			
			private final HttpExchange t;
			private DisplayParams params;
			private final long instanceStartTime;
			private boolean fileMode=false;
			private boolean refMode=false;
		}
		
	}
	
	private Sketch findRefSketch(String s){
		assert(s!=null);
		if(s==null || s.length()<1){return null;}
		int tid=-1;
		if(Tools.isDigit(s.charAt(0))){tid=Integer.parseInt(s);}
		else{
			TaxNode tn=getTaxNodeByName(s);
			tid=tn==null ? -1 : tn.id;
		}
		Sketch sk=tid<0 ? null : searcher.findReferenceSketch(tid);
		if(sk!=null){sk=(Sketch)sk.clone();}
		return sk;
	}
	
	/*--------------------------------------------------------------*/

	/** Handles taxonomy lookups */
	class TaxHandler implements HttpHandler {
		
		public TaxHandler(boolean skipNonCanonical_){
			skipNonCanonical=skipNonCanonical_;
		}
		
		@Override
		public void handle(HttpExchange t) throws IOException {
			if(verbose2){System.err.println("Tax handler");}
			final long startTime=System.nanoTime();
			
			if(sketchOnly){
				ServerTools.reply("\nERROR: This server is tunning in sketch mode and should not be used for taxonomic lookups.\n"
						+ "The taxonomy server is at "+Shared.taxServer()+"\n", "text/plain", t, verbose2, 400, true);
				return;
			}
			
			//Parse the query from the URL
			String rparam=getRParam(t, true);
			

			boolean simple=skipNonCanonical;
			
			{//Legacy support for old style of invoking simple
				if(rparam.startsWith("simpletax/")){rparam=rparam.substring(7); simple=true;}
				else if(rparam.startsWith("stax/")){rparam=rparam.substring(5); simple=true;}
				else if(rparam.startsWith("tax/")){rparam=rparam.substring(4);}
				else if(rparam.equals("simpletax") || rparam.equals("stax")){rparam=""; simple=true;}
				else if(rparam.equals("tax")){rparam="";}
			}
			while(rparam.startsWith("/")){rparam=rparam.substring(1);}
			if(rparam.length()<1 || rparam.equalsIgnoreCase("help") || rparam.equalsIgnoreCase("usage")){
				returnUsage(startTime, t);
				return;
			}

			String[] params = rparam.split("/");
			if(verbose2){System.err.println(Arrays.toString(params));}
			
			final String response=toResponse(simple, params, t);
			final String type=response.startsWith("{") ? "application/json" : "text/plain";
			
			ServerTools.reply(response, type, t, verbose2, 200, true);
			
			final long stopTime=System.nanoTime();
			final long elapsed=stopTime-startTime;
			if(response.startsWith("Welcome to ")){
				timeMeasurementsUsage.incrementAndGet();
				elapsedTimeUsage.addAndGet(elapsed);
				lastTimeUsage.set(elapsed);
			}else{
				timeMeasurementsRemote.incrementAndGet();
				elapsedTimeRemote.addAndGet(elapsed);
				lastTimeRemote.set(elapsed);
			}
		}
		
		//TODO: Integrate something like this to improve parsing
//		String parse(String rparam){
//
//			if(rparam.length()<2){return rparam;}
//
//			String[] args=rparam.split("/");
//			int trimmed=0;
//			for(int i=0; i<args.length; i++){//parse rparam
//				String arg=args[i];
//				if(arg.length()>0){
//					String[] split=arg.split("=");
//					String a=split[0].toLowerCase();
//					String b=split.length>1 ? split[1] : null;
//
//					if(a.equals("file")){
//						fileMode=true;
//						trimmed+=5;
//						break;
//					}else if(params.parse(arg, a, b)){
//						trimmed+=arg.length()+1;
//					}else{
//						assert(false) : "Bad argument:'"+arg+"'"+"\n"+Arrays.toString(args)+"\n"+rparam;
//					}
//				}
//			}
//			params.postParse(true);
//			//			System.err.println("Trimmed="+trimmed+", rparam="+rparam);
//			if(trimmed>0){
//				rparam=rparam.substring(Tools.min(trimmed, rparam.length()));
//			}
//		}
		
		/** Only print nodes at canonical tax levels */
		public final boolean skipNonCanonical;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Helpers           ----------------*/
	/*--------------------------------------------------------------*/
	
	static String getBody(HttpExchange t){
		if(verbose2){System.err.println("getBody");}
		InputStream is=t.getRequestBody();
		String s=ServerTools.readStream(is);
		return s;
	}
	
	/** Parse the query from the URL */
	static String getRParam(HttpExchange t, boolean allowPost){
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
		
		if(allowPost && ("$POST".equalsIgnoreCase(rparam) || "POST".equalsIgnoreCase(rparam))){
			String body=getBody(t);
			rparam=body;
		}
		
		if(verbose){System.err.println(rparam==null || rparam.trim().length()<1 ? "usage" : rparam+"\t"+System.currentTimeMillis());}
		return rparam;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Taxonomy Formatting     ----------------*/
	/*--------------------------------------------------------------*/
	
	/** All tax queries enter here from the handler */
	String toResponse(boolean simple, String[] params, HttpExchange t){
		if(verbose2){System.err.println("toResponse");}

		boolean printNumChildren=false;
		boolean printChildren=false;
		boolean printPath=false;
		boolean printSize=false;
		boolean printRange=false;
		boolean silvaHeader=false;
		boolean plaintext=false, semicolon=false, path=false;
		boolean mononomial=false;
		boolean ancestor=false;
		int source=SOURCE_REFSEQ;
		
		ArrayList<String> newParams=new ArrayList<String>(params.length);
		for(int i=0; i<params.length-1; i++){
			String s=params[i];
			if(s.length()==0 || s.equals("tax")){
				//Do nothing
			}else if(s.equals("printnumchildren") || s.equals("numchildren")){
				printNumChildren=true;
			}else if(s.equals("printchildren") || s.equals("children")){
				printChildren=true;
			}else if(s.equals("mononomial") || s.equals("cn") 
					|| s.equals("fixname") || s.equals("fn") || s.equals("mono") || s.equals("mononomial")){
				mononomial=true;
			}else if(s.equals("printpath") || s.equals("pp")){
				printPath=true;
			}else if(s.equals("printsize") || s.equals("size") || s.equals("ps")){
				printSize=true;
			}else if(s.equals("printrange") || s.equals("range")){
				printRange=true;
			}else if(s.equals("plaintext") || s.equals("pt")){
				plaintext=true; semicolon=false; path=false;
			}else if(s.equals("semicolon") || s.equals("sc")){
				semicolon=true; plaintext=false; path=false;
			}else if(s.equals("path") || s.equals("pa")){
				path=true; semicolon=false; plaintext=false;
			}else if(s.equals("silva")){
				silvaHeader=true;
				source=SOURCE_SILVA;
			}else if(s.equals("refseq")){
				source=SOURCE_REFSEQ;
			}else if(s.equals("ancestor")){
				ancestor=true;
			}else if(s.equals("simple")){
				simple=true;
			}else{
				newParams.add(s);
			}
		}
		newParams.add(params[params.length-1]);
		params=newParams.toArray(new String[newParams.size()]);
		
//		System.err.println("mononomial="+mononomial);
		
		if(params.length<2){
			if(params.length==1 && "advice".equalsIgnoreCase(params[0])){return TAX_ADVICE;}
			if(logUsage){System.err.println("usage");}
			return USAGE(USAGE);
		}
		if(params.length>3){
			if(logUsage){System.err.println("usage");}
			return USAGE(USAGE);
		}
		
		final String query=params[params.length-1];
		final String[] names=query.split(",");
		if(names==null){return USAGE(USAGE);}
		
		if(names==null || names.length<2){
			firstChunkSingle.incrementAndGet();
		}else{
			firstChunkMulti.incrementAndGet();
		}
		bulkCount.addAndGet(names.length);
//		System.err.println(params[2]+", "+ancestor);
		
		//Raw query type code
		final int type;
		//Type code excluding formatting
		final int type2;
		{
			String typeS=params[0];
			Integer value=typeMap.get(typeS);
			if(value==null){
				if(typeS.equalsIgnoreCase("advice")){
					return TAX_ADVICE;
				}else{
					return "{\"error\": \"Bad type ("+typeS+"); should be gi, taxid, or name.\"}";
				}
			}
			int x=value.intValue();
			if((x&15)==HEADER && silvaHeader){x=SILVAHEADER;}
			type=x;
			type2=x&15;
			if(type2==IMG){source=SOURCE_IMG;}
		}
		
		plaintext=(type>=PT_BIT || plaintext);
		semicolon=(type>=SC_BIT || semicolon);
		path=(type>=PA_BIT || path);
		if(semicolon || path){plaintext=false;}
		if(path){semicolon=false;}
		
		final boolean internal=incrementQueries(t, false, false, simple, ancestor, 
				plaintext, semicolon, path, printChildren, printPath, printSize, type); //Ignores usage information.
		
//		if(type2==GI){//123
//			return "{\"error\": \"GI number support is temporarily suspended due to conflicts in NCBI databases.  "
//					+ "It may come back split into nucleotide and protein GI numbers, which currently are not exclusive.\"}";
//		}
		
		if(!internal && !allowRemoteFileAccess){
			path=printPath=false;
		}
		
		if(verbose2){System.err.println("Type: "+type);}
		if(type2==NAME || type2==HEADER || type2==SILVAHEADER){
			for(int i=0; i<names.length; i++){
				names[i]=PercentEncoding.codeToSymbol(names[i]);
				if(type2==HEADER || type2==SILVAHEADER){
					if(names[i].startsWith("@") || names[i].startsWith(">")){names[i]=names[i].substring(1);}
				}
			}
			if(verbose2){System.err.println("Revised: "+Arrays.toString(names));}
		}
		
		if(ancestor){
			if(verbose2){System.err.println("toAncestor: "+Arrays.toString(names));}
			return toAncestor(type, names, plaintext, semicolon, path, query, simple, !simple, printNumChildren, printChildren, printPath, printSize, printRange, mononomial, source);
		}
		
		if(semicolon){
			return toSemicolon(type, names, simple, mononomial);
		}else if(plaintext){
			return toText(type, names);
		}else if(path){
			return toPath(type, names, source);
		}

		JsonObject j=new JsonObject();
		for(String name : names){
			j.add(name, toJson(type, name, simple, !simple, printNumChildren, printChildren, 
					printPath, printSize, printRange, mononomial, source));
		}
		return j.toString();
	}
	
	/** Look up common ancestor of terms */
	String toAncestor(final int type, final String[] names, boolean plaintext, boolean semicolon, boolean path,
			String query, final boolean skipNonCanonical, boolean originalLevel, boolean printNumChildren,
			boolean printChildren, boolean printPath, boolean printSize, boolean printRange, boolean mononomial,
			int source){
		IntList ilist=toIntList(type, names);
		int id=FindAncestor.findAncestor(tree, ilist);
		TaxNode tn=(id>-1 ? tree.getNode(id) : null);
		if(tn==null){
			return new JsonObject("error","Not found.").toString(query);
		}
		if(semicolon){
			return tree.toSemicolon(tn, skipNonCanonical, mononomial);
		}else if(plaintext){
			return ""+id;
		}else if(path){
			return toPath(tn, source);
		}
		
		JsonObject j=new JsonObject();
//		j.add("name", mononomial ? tree.mononomial(tn) : tn.name);
		j.add("name", tn.name);
		if(mononomial || true){
			String mono=tree.mononomial(tn);
			if(tn.name!=mono){j.add("mononomial", mono);}
		}
		j.add("tax_id", tn.id);
		if(printNumChildren){j.add("num_children", tn.numChildren);}
		if(printPath){j.add("path", toPath(tn, source));}
		if(printSize){
			j.add("size", tree.toSize(tn));
			j.add("cumulative_size", tree.toSizeC(tn));
			j.add("seqs", tree.toSeqs(tn));
			j.add("cumulative_seqs", tree.toSeqsC(tn));
			j.add("cumulative_nodes", tree.toNodes(tn));
		}
		j.add("level", tn.levelStringExtended(originalLevel));
		if(tn.levelExtended<1 && printRange){
			j.add("maxDescendent", TaxTree.levelToStringExtended(tn.maxChildLevelExtended));
			j.add("minAncestor", TaxTree.levelToStringExtended(tn.minParentLevelExtended));
		}
//		if(printChildren){j.add(getChildren(id, originalLevel, printRange));}
		while(tn!=null && tn.levelExtended!=TaxTree.LIFE_E && tn.id!=TaxTree.CELLULAR_ORGANISMS_ID){
			if(!skipNonCanonical || tn.isSimple()){
				j.addAndRename(tn.levelStringExtended(originalLevel), toJson(tn, originalLevel, printNumChildren, printChildren, printPath, printSize, printRange, mononomial, source, -1));
			}
			if(tn.pid==tn.id){break;}
			tn=tree.getNode(tn.pid);
		}
		return j.toString();
	}
	
	JsonObject getChildren(final int id, boolean originalLevel, boolean printRange, boolean mononomial){
		TaxNode x=tree.getNode(id);
		if(x==null || x.numChildren==0){return null;}
		ArrayList<TaxNode> list=tree.getChildren(x);
		return makeChildrenObject(list, originalLevel, printRange, mononomial);
	}
	
	JsonObject makeChildrenObject(ArrayList<TaxNode> list, boolean originalLevel, boolean printRange, boolean mononomial){
		if(list==null || list.isEmpty()){return null;}
		JsonObject j=new JsonObject();
		for(TaxNode tn : list){
			JsonObject child=new JsonObject();
//			child.add("name", mononomial ? tree.mononomial(tn) : tn.name);
			child.add("name", tn.name);
			if(mononomial || true){
				String mono=tree.mononomial(tn);
				if(tn.name!=mono){child.add("mononomial", mono);}
			}
			child.add("tax_id", tn.id);
			child.add("num_children", tn.numChildren);
			child.add("level", tn.levelStringExtended(originalLevel));
			if(tn.levelExtended<1 && printRange){
				child.add("maxDescendent", TaxTree.levelToStringExtended(tn.maxChildLevelExtended));
				child.add("minAncestor", TaxTree.levelToStringExtended(tn.minParentLevelExtended));
			}
			j.add(tn.id+"", child);
		}
		return j;
	}
	
	/** Format a reply as plaintext, comma-delimited, TaxID only */
	String toText(final int type, final String[] names){
		
		StringBuilder sb=new StringBuilder();
		String comma="";

		int type2=type&15;
		if(type2==GI){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeGi(Long.parseLong(name));
				if(tn==null){sb.append("-1");}
				else{sb.append(tn.id);}
				comma=",";
			}
		}else if(type2==NAME){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeByName(name);
				if(tn==null){sb.append("-1");}
				else{sb.append(tn.id);}
				comma=",";
			}
		}else if(type2==TAXID){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeTaxid(Integer.parseInt(name));
				if(tn==null){sb.append("-1");}
				else{sb.append(tn.id);}
				comma=",";
			}
		}else if(type2==ACCESSION){
			for(String name : names){
				sb.append(comma);
				int ncbi=accessionToTaxid(name);
				sb.append(ncbi);
				comma=",";
			}
		}else if(type2==HEADER || type2==SILVAHEADER){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeHeader(name, type2==SILVAHEADER);
				if(tn==null){sb.append("-1");}
				else{sb.append(tn.id);}
				comma=",";
			}
		}else if(type2==IMG){
			for(String name : names){
				sb.append(comma);
				int ncbi=imgToTaxid(Long.parseLong(name));
				sb.append(ncbi);
				comma=",";
			}
		}else{
			return "Bad type; should be pt_gi or pt_name; e.g. /pt_gi/1234";
		}
		
		return sb.toString();
	}
	
	private TaxNode toNode(final int type, final String name){
		int type2=type&15;
		final TaxNode tn;
		if(type2==GI){
			tn=getTaxNodeGi(Long.parseLong(name));
		}else if(type2==NAME){
			tn=getTaxNodeByName(name);
		}else if(type2==TAXID){
			tn=getTaxNodeTaxid(Integer.parseInt(name));
		}else if(type2==ACCESSION){
			int ncbi=accessionToTaxid(name);
			tn=(ncbi<0 ? null : tree.getNode(ncbi));
		}else if(type2==HEADER || type2==SILVAHEADER){
			tn=getTaxNodeHeader(name, type2==SILVAHEADER);
		}else if(type2==IMG){
			int ncbi=imgToTaxid(Long.parseLong(name));
			tn=(ncbi<0 ? null : tree.getNode(ncbi));
		}else{
			tn=null;
		}
		return tn;
	}
	
	/** Format a reply as paths, comma-delimited*/
	String toPath(final int type, final String[] names, final int source){
		
		StringBuilder sb=new StringBuilder();
		String comma="";

		int type2=type&15;
		
		for(String name : names){
			sb.append(comma);
			if(type2==IMG){
				sb.append(toPathIMG(Long.parseLong(name)));
			}else{
				TaxNode tn=toNode(type2, name);
				sb.append(toPath(tn, source));
			}
			comma=",";
		}
		
		return sb.toString();
	}
	
	/** Format a reply as plaintext, semicolon-delimited, full lineage */
	String toSemicolon(final int type, final String[] names, boolean skipNonCanonical, boolean mononomial){
		
		StringBuilder sb=new StringBuilder();
		String comma="";
		
		int type2=type&15;
		if(type2==GI){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeGi(Long.parseLong(name));
				if(tn==null){sb.append("Not found");}
				else{sb.append(tree.toSemicolon(tn, skipNonCanonical, mononomial));}
				comma=",";
			}
		}else if(type2==NAME){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeByName(name);
				if(tn==null){sb.append("Not found");}
				else{sb.append(tree.toSemicolon(tn, skipNonCanonical, mononomial));}
				comma=",";
			}
		}else if(type2==TAXID){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeTaxid(Integer.parseInt(name));
//				if(verbose2){outstream.println("name="+name+", tn="+tn);}
				if(tn==null){sb.append("Not found");}
				else{sb.append(tree.toSemicolon(tn, skipNonCanonical, mononomial));}
				comma=",";
			}
		}else if(type2==ACCESSION){
			for(String name : names){
				sb.append(comma);
				final int tid=accessionToTaxid(name);
				TaxNode tn=tree.getNode(tid, true);
				if(tn==null){sb.append("Not found");}
				else{sb.append(tree.toSemicolon(tn, skipNonCanonical, mononomial));}
				comma=",";
			}
		}else if(type2==HEADER || type2==SILVAHEADER){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeHeader(name, type2==SILVAHEADER);
				if(tn==null){sb.append("Not found");}
				else{sb.append(tree.toSemicolon(tn, skipNonCanonical, mononomial));}
				comma=",";
			}
		}else if(type2==IMG){
			for(String name : names){
				sb.append(comma);
				final int tid=imgToTaxid(Long.parseLong(name));
				TaxNode tn=tree.getNode(tid, true);
				if(tn==null){sb.append("Not found");}
				else{sb.append(tree.toSemicolon(tn, skipNonCanonical, mononomial));}
				comma=",";
			}
		}else{
			return "Bad type; should be sc_gi or sc_name; e.g. /sc_gi/1234";
		}
		
//		if(verbose2){outstream.println("In toSemicolon; type="+type+", type2="+type2+", made "+sb);}
		
		return sb.toString();
	}
	
	/** Create a JsonObject from a String, including full lineage */
	JsonObject toJson(final int type, final String name, boolean skipNonCanonical, boolean originalLevel, 
			boolean printNumChildren, boolean printChildren, boolean printPath, boolean printSize, 
			boolean printRange, boolean mononomial, int source){
		final TaxNode tn0;
		TaxNode tn;
		
		long img=-1;
		if(type==GI){
			tn0=getTaxNodeGi(Long.parseLong(name));
		}else if(type==NAME){
			tn0=getTaxNodeByName(name);
		}else if(type==TAXID){
			tn0=getTaxNodeTaxid(Integer.parseInt(name));
		}else if(type==ACCESSION){
			int ncbi=accessionToTaxid(name);
			tn0=(ncbi>=0 ? tree.getNode(ncbi) : null);
		}else if(type==HEADER || type==SILVAHEADER){
			tn0=getTaxNodeHeader(name, type==SILVAHEADER);
		}else if(type==IMG){
			img=Long.parseLong(name);
			final int tid=imgToTaxid(img);
			tn0=tree.getNode(tid, true);
		}else{
			JsonObject j=new JsonObject("error","Bad type; should be gi, taxid, or name; e.g. /name/homo_sapiens");
			j.add("name", name);
			j.add("type", type);
			return j;
		}
		tn=tn0;
		if(verbose2){System.err.println("Got node: "+tn);}
		
		if(tn!=null){
			JsonObject j=new JsonObject();
//			j.add("name", mononomial ? tree.mononomial(tn) : tn.name);
			j.add("name", tn.name);
			if(mononomial || true){
				String mono=tree.mononomial(tn);
				if(tn.name!=mono){j.add("mononomial", mono);}
			}
			j.add("tax_id", tn.id);
			if(printNumChildren){j.add("num_children", tn.numChildren);}
			if(printPath){j.add("path", type==IMG ? toPathIMG(img) : toPath(tn, source));}
			if(printSize){
				j.add("size", tree.toSize(tn));
				j.add("cumulative_size", tree.toSizeC(tn));
				j.add("seqs", tree.toSeqs(tn));
				j.add("cumulative_seqs", tree.toSeqsC(tn));
				j.add("cumulative_nodes", tree.toNodes(tn));
			}
			j.add("level", tn.levelStringExtended(originalLevel));
			if(tn.levelExtended<1 && printRange){
				j.add("maxDescendent", TaxTree.levelToStringExtended(tn.maxChildLevelExtended));
				j.add("minAncestor", TaxTree.levelToStringExtended(tn.minParentLevelExtended));
			}
			if(printChildren && (tn.id==tn.pid || tn.id==TaxTree.CELLULAR_ORGANISMS_ID)){j.add("children", getChildren(tn.id, originalLevel, printRange, mononomial));}
			while(tn!=null && tn.levelExtended!=TaxTree.LIFE_E && tn.id!=TaxTree.CELLULAR_ORGANISMS_ID){
//				System.err.println(tn+", "+(!skipNonCanonical)+", "+tn.isSimple());
				if(!skipNonCanonical || tn.isSimple()){
					j.addAndRename(tn.levelStringExtended(originalLevel), toJson(tn, originalLevel, printNumChildren, printChildren, printPath && tn==tn0, printSize, printRange, mononomial, source, img));
//					System.err.println(j);
				}
				if(tn.pid==tn.id){break;}
				tn=tree.getNode(tn.pid);
			}
			return j;
		}
		{
			JsonObject j=new JsonObject("error","Not found.");
			j.add("name", name);
			j.add("type", type);
			return j;
		}
	}
	
	/** Create a JsonObject from a TaxNode, at that level only */
	JsonObject toJson(TaxNode tn, boolean originalLevel, boolean printNumChildren, 
			boolean printChildren, boolean printPath, boolean printSize, boolean printRange, 
			boolean mononomial, int source, long img){
		JsonObject j=new JsonObject();
//		j.add("name", mononomial ? tree.mononomial(tn) : tn.name);
		j.add("name", tn.name);
		if(mononomial || true){
			String mono=tree.mononomial(tn);
			if(tn.name!=mono){j.add("mononomial", mono);}
		}
		j.add("tax_id", tn.id);
		if(printNumChildren){j.add("num_children", tn.numChildren);}
		if(printPath){j.add("path", source==SOURCE_IMG ? toPathIMG(img) : toPath(tn, source));}
		if(printSize){
			j.add("size", tree.toSize(tn));
			j.add("cumulative_size", tree.toSizeC(tn));
			j.add("seqs", tree.toSeqs(tn));
			j.add("cumulative_seqs", tree.toSeqsC(tn));
			j.add("cumulative_nodes", tree.toNodes(tn));
		}
		if(tn.levelExtended<1 && printRange){
			j.add("maxDescendent", TaxTree.levelToStringExtended(tn.maxChildLevelExtended));
			j.add("minAncestor", TaxTree.levelToStringExtended(tn.minParentLevelExtended));
		}
		if(printChildren){
			JsonObject children=getChildren(tn.id, originalLevel, printRange, mononomial);
			if(children!=null){j.add("children", children);}
		}
		return j;
	}
	
	String toPath(TaxNode tn, int source){
		if(tn==null){return "null";}
		String path;
		if(source==SOURCE_REFSEQ){
			path=tree.toDir(tn, basePath)+"refseq_"+tn.id+".fa.gz";
		}else if(source==SOURCE_SILVA){
			path=tree.toDir(tn, basePath)+"silva_"+tn.id+".fa.gz";
		}else if(source==SOURCE_IMG){
			assert(false);
			path="null";
		}else{
			assert(false);
			path="null";
		}
		if(!path.equals("null") && !new File(path).exists()){path="null";}
		return path;
	}
	
	String toPathIMG(long imgID){
		String path="/global/dna/projectdirs/microbial/img_web_data/taxon.fna/"+imgID+".fna";
		if(!new File(path).exists()){path="null";}
		return path;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Taxonomy Lookup       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Convert a list of terms to a list of TaxIDs */
	IntList toIntList(final int type, final String[] names){
		IntList list=new IntList(names.length);
		int type2=type&15;
		if(type2==GI){
			for(String name : names){
				TaxNode tn=getTaxNodeGi(Long.parseLong(name));
				if(tn!=null){list.add(tn.id);}
				else{notFound.incrementAndGet();}
			}
		}else if(type2==NAME){
			for(String name : names){
				TaxNode tn=getTaxNodeByName(name);
				if(tn!=null){list.add(tn.id);}
				else{notFound.incrementAndGet();}
			}
		}else if(type2==TAXID){
			for(String name : names){
				TaxNode tn=getTaxNodeTaxid(Integer.parseInt(name));
				if(tn!=null){list.add(tn.id);}
				else{notFound.incrementAndGet();}
			}
		}else if(type2==ACCESSION){
			for(String name : names){
				int ncbi=accessionToTaxid(name);
				if(ncbi>=0){list.add(ncbi);}
				else{notFound.incrementAndGet();}
			}
		}else if(type2==IMG){
			for(String name : names){
				final int tid=imgToTaxid(Long.parseLong(name));
				if(tid>=0){list.add(tid);}
				else{notFound.incrementAndGet();}
			}
		}else{
			throw new RuntimeException("{\"error\": \"Bad type\"}");
		}
		return list;
	}
	
	public static final String stripAccession(String s){
		if(s==null){return null;}
		s=s.toUpperCase();
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c=='.' || c==':'){return s.substring(0, i);}
		}
		return s;
	}

	private final int imgToTaxid(String imgID) {
		return imgToTaxid(Long.parseLong(imgID));
	}
	private final int imgToTaxid(long imgID) {
		int ncbi=-1;
		if(imgFile!=null) {
			ncbi=TaxTree.imgToTaxid(imgID);
			if(ncbi>=0) {return ncbi;}
		}
		return AccessionToTaxid.getImg(imgID);
	}
	
	private int accessionToTaxid(String accession){
		if(accession==null){return -1;}
		int tid=AccessionToTaxid.get(accession);
		if(tid<0 && distributed && serverNum==0){
			accession=stripAccession(accession);
			int slaveNum=accession.hashCode()%serverCount;
			if(slaveNum!=serverNum){
				String path=slaveAddress.get(slaveNum);
				tid=TaxClient.accessionToTaxidSpecificServer(path, accession);
			}
		}
		return tid;
	}
	
	/** Look up a TaxNode by parsing the organism name */
	TaxNode getTaxNodeByName(String name){
		if(verbose2){System.err.println("Fetching node for "+name);}
		List<TaxNode> list=tree.getNodesByNameExtended(name);
		if(verbose2){System.err.println("Fetched "+list);}
		if(list==null){
			if(verbose2){System.err.println("Fetched in common map "+name);}
			String name2=commonMap.get(name);
			if(verbose2){System.err.println("Fetched "+name2);}
			if(name2!=null){list=tree.getNodesByName(name2);}
		}
		if(list==null){notFound.incrementAndGet();}
		return list==null ? null : list.get(0);
	}
	
	/** Look up a TaxNode from the gi number */
	TaxNode getTaxNodeGi(long gi){
		int ncbi=-1;
		try {
			ncbi=GiToTaxid.getID(gi);
		} catch (Throwable e) {
			if(verbose){e.printStackTrace();}
		}
		if(ncbi<0){notFound.incrementAndGet();}
		return ncbi<0 ? null : getTaxNodeTaxid(ncbi);
	}
	
	/** Look up a TaxNode by parsing the full header */
	TaxNode getTaxNodeHeader(String header, boolean silvaMode){
		TaxNode tn=silvaMode ? tree.getNodeSilva(header, true) : tree.parseNodeFromHeader(header, true);
		if(tn==null){notFound.incrementAndGet();}
		return tn;
	}
	
	/** Look up a TaxNode from the ncbi TaxID */
	TaxNode getTaxNodeTaxid(int ncbi){
		TaxNode tn=null;
		try {
			tn=tree.getNode(ncbi);
		} catch (Throwable e) {
			if(verbose){e.printStackTrace();}
		}
		if(tn==null){notFound.incrementAndGet();}
		return tn;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------     Data Initialization      ----------------*/
	/*--------------------------------------------------------------*/

	private static HashMap<String, Integer> makeTypeMap() {
		HashMap<String, Integer> map=new HashMap<String, Integer>(63);
		map.put("gi", GI);
//		map.put("ngi", NGI);
//		map.put("pgi", PGI);
		map.put("name", NAME);
		map.put("tax_id", TAXID);
		map.put("ncbi", TAXID);
		map.put("taxid", TAXID);
		map.put("id", TAXID);
		map.put("tid", TAXID);
		map.put("header", HEADER);
		map.put("accession", ACCESSION);
		map.put("img", IMG);
		map.put("silvaheader", SILVAHEADER);
		
		map.put("pt_gi", PT_GI);
//		map.put("pt_ngi", PT_NGI);
//		map.put("pt_pgi", PT_PGI);
		map.put("pt_name", PT_NAME);
		map.put("pt_tax_id", PT_TAXID);
		map.put("pt_id", PT_TAXID);
		map.put("pt_tid", PT_TAXID);
		map.put("pt_ncbi", PT_TAXID);
		map.put("pt_taxid", PT_TAXID);
		map.put("pt_header", PT_HEADER);
		map.put("pt_header", PT_HEADER);
		map.put("pt_accession", PT_ACCESSION);
		map.put("pt_img", PT_IMG);
		map.put("pt_silvaheader", PT_SILVAHEADER);

		map.put("sc_gi", SC_GI);
//		map.put("sc_ngi", SC_NGI);
//		map.put("sc_pgi", SC_PGI);
		map.put("sc_name", SC_NAME);
		map.put("sc_tax_id", SC_TAXID);
		map.put("sc_id", SC_TAXID);
		map.put("sc_tid", SC_TAXID);
		map.put("sc_ncbi", SC_TAXID);
		map.put("sc_taxid", SC_TAXID);
		map.put("sc_header", SC_HEADER);
		map.put("sc_header", SC_HEADER);
		map.put("sc_accession", SC_ACCESSION);
		map.put("sc_silvaheader", SC_SILVAHEADER);
	
		return map;
	}
	
	public static HashMap<String, String> makeCommonMap(){
		HashMap<String, String> map=new HashMap<String, String>();
		map.put("human", "homo sapiens");
		map.put("cat", "felis catus");
		map.put("dog", "canis lupus familiaris");
		map.put("mouse", "mus musculus");
		map.put("cow", "bos taurus");
		map.put("bull", "bos taurus");
		map.put("horse", "Equus ferus");
		map.put("pig", "Sus scrofa domesticus");
		map.put("sheep", "Ovis aries");
		map.put("goat", "Capra aegagrus");
		map.put("turkey", "Meleagris gallopavo");
		map.put("fox", "Vulpes vulpes");
		map.put("chicken", "Gallus gallus domesticus");
		map.put("wolf", "canis lupus");
		map.put("fruitfly", "drosophila melanogaster");
		map.put("zebrafish", "Danio rerio");
		map.put("catfish", "Ictalurus punctatus");
		map.put("trout", "Oncorhynchus mykiss");
		map.put("salmon", "Salmo salar");
		map.put("tilapia", "Oreochromis niloticus");
		map.put("e coli", "Escherichia coli");
		map.put("e.coli", "Escherichia coli");

		map.put("lion", "Panthera leo");
		map.put("tiger", "Panthera tigris");
		map.put("bear", "Ursus arctos");
		map.put("deer", "Odocoileus virginianus");
		map.put("coyote", "Canis latrans");

		map.put("corn", "Zea mays subsp. mays");
		map.put("maize", "Zea mays subsp. mays");
		map.put("oat", "Avena sativa");
		map.put("wheat", "Triticum aestivum");
		map.put("rice", "Oryza sativa");
		map.put("potato", "Solanum tuberosum");
		map.put("barley", "Hordeum vulgare");
		map.put("poplar", "Populus alba");
		map.put("lettuce", "Lactuca sativa");
		map.put("beet", "Beta vulgaris");
		map.put("strawberry", "Fragaria x ananassa");
		map.put("orange", "Citrus sinensis");
		map.put("lemon", "Citrus limon");
		map.put("soy", "Glycine max");
		map.put("soybean", "Glycine max");
		map.put("grape", "Vitis vinifera");
		map.put("olive", "Olea europaea");
		map.put("cotton", "Gossypium hirsutum");
		map.put("apple", "Malus pumila");
		map.put("bannana", "Musa acuminata");
		map.put("tomato", "Solanum lycopersicum");
		map.put("sugarcane", "Saccharum officinarum");
		map.put("bean", "Phaseolus vulgaris");
		map.put("onion", "Allium cepa");
		map.put("garlic", "Allium sativum");
		
		map.put("pichu", "mus musculus");
		map.put("pikachu", "mus musculus");
		map.put("vulpix", "Vulpes vulpes");
		map.put("ninetails", "Vulpes vulpes");
		map.put("mareep", "Ovis aries");
		
		return map;
	}
	
	//Customize usage message to include domain
	private String makeUsagePrefix(){
		if(!sketchOnly){
			return "Welcome to the JGI taxonomy server!\n"
					+ "This service provides taxonomy information from NCBI taxID numbers, gi numbers, organism names, and accessions.\n"
					+ "The output is formatted as a Json object.\n\n"
					+ "Usage:\n\n"
					+ "All addresses below are assumed to be prefixed by "+domain+", e.g. /name/homo_sapiens implies a full URL of:\n"
					+ domain+"/name/homo_sapiens\n"
					+ "\n"
					+ "/name/homo_sapiens will give taxonomy information for an organism name.\n"
					+ "Names are case-insensitive and underscores are equivalent to spaces.\n"
					+ "/id/9606 will give taxonomy information for an NCBI taxID.\n"
					+ "/gi/1234 will give taxonomy information from an NCBI gi number.\n"
					
//					+ "\n****NOTICE**** gi number support is temporarily suspended due to conflicts in NCBI data.\n"
//					+ "Support may be restored, altered, or discontinued pending a response from NCBI.\n"
//					+ "Currently, it is not possible to ensure correct results when looking up a GI number, because some map to multiple organisms.\n\n"
					
					+ "/accession/NZ_AAAA01000057.1 will give taxonomy information from an accession.\n"
					+ "/header/ will accept an NCBI sequence header such as gi|7|emb|X51700.1| Bos taurus\n"
					+ "/silvaheader/ will accept a Silva sequence header such as KC415233.1.1497 Bacteria;Spirochaetae;Spirochaetes\n"
					+ "/img/ will accept an IMG id such as 2724679250\n"
					+ "Vertical bars (|) may cause problems on the command line and can be replaced by tilde (~).\n"
					+ "\nComma-delimited lists are accepted for bulk queries, such as tax/gi/1234,7000,42\n"
					+ "For plaintext (non-Json) results, add the term /pt/ or /sc/.\n"
					+ "pt will give just the taxID, while sc will give the whole lineage, semicolon-delimited. For example:\n"
					+ "/pt/name/homo_sapiens\n"
					+ "/sc/gi/1234\n\n"
					+ "Additional supported display options are children, numchildren, range, simple, path, size, and ancestor.\n"
					+ "The order is not important but they need to come before the query term.  For example:\n"
					+ "/children/numchildren/range/gi/1234\n"
					+ "\nTo find the common ancestor of multiple organisms, add /ancestor/. For example:\n"
					+ "/id/ancestor/1234,5678,42\n"
					+ "/name/ancestor/homo_sapiens,canis_lupus,bos_taurus\n"
					+ "\nFor a simplified taxonomic tree, add simple.\n"
					+ "This will ignore unranked or uncommon levels like tribe and parvorder, and only display the following levels:\n"
					+ "SUBSPECIES, SPECIES, GENUS, FAMILY, ORDER, CLASS, PHYLUM, KINGDOM, SUPERKINGDOM, DOMAIN\n"
					+ "For example:\n"
					+ "/simple/id/1234\n"
					+ "\nTo print taxonomy from the command line in Linux, use curl:\n"
					+ "curl https://taxonomy.jgi.doe.gov/id/9606\n"
					+ "\nQueries longer than around 8kB can be sent via POST: curl https://taxonomy..doe.gov/POST"
					+ "\n...where the data sent is, for example: name/e.coli,h.sapiens,c.lupus\n"
					+ "\nLast restarted "+startTime+"\n"
					+ "Running BBTools version "+Shared.BBTOOLS_VERSION_STRING+"\n";
		}else{
			StringBuilder sb=new StringBuilder();
			sb.append("Welcome to the JGI"+(SketchObject.defaultParams.dbName==null ? "" : " "+SketchObject.defaultParams.dbName)+" sketch server!\n");
//			if(dbName!=null){
//				sb.append("This server has the "+dbName+ " database loaded.\n");
//			}
			sb.append("\nUsage:\n\n");
			sb.append("sendsketch.sh in=file.fasta"+(SketchObject.defaultParams.dbName==null ? "" : " "+SketchObject.defaultParams.dbName.toLowerCase())+"\n\n");
			sb.append("SendSketch creates a sketch from a local sequence file, and sends the sketch to this server.\n");
			sb.append("The server receives the sketch, compares it to all sketches in memory, and returns the results.\n");
			sb.append("For files on the same system as the server, the 'local' flag may be used to offload sketch creation to the server.\n");
			sb.append("For more details and parameters please run sendsketch.sh with no arguments.\n");
			sb.append("\n");
			if(SketchObject.useWhitelist()){
				sb.append("This server is running in whitelist mode; for best results, use local queries.\n");
				sb.append("Remote queries should specify a larger-than-normal sketch size.\n\n");
			}else if(SketchObject.blacklist()!=null){
				sb.append("This server is running in blacklist mode, using "+new File(SketchObject.blacklist()).getName()+".\n\n");
			}
			sb.append("Last restarted "+startTime+"\n");
			sb.append("Running BBTools version "+Shared.BBTOOLS_VERSION_STRING+"\n");
			sb.append("Settings:\tk="+SketchObject.k+(SketchObject.k2>0 ? ","+SketchObject.k2 : ""));
			if(SketchObject.amino){sb.append(" amino");}
			if(SketchObject.makeIndex){sb.append(" index");}
			if(SketchObject.useWhitelist()){sb.append(" whitelist");}
			if(SketchObject.blacklist()!=null){sb.append(" blacklist="+new File(SketchObject.blacklist()).getName());}
			sb.append('\n');
			return sb.toString();
		}
	}
	
	private String makeUsageHtml(){
		String html=rawHtml;
		html=html.replace("STATISTICSSTRING", makeStats());
//		html=html.replace("TIMESTAMPSTRING", startTime);
//		html=html.replace("VERSIONSTRING", "Running BBTools version "+Shared.BBMAP_VERSION_STRING);
		return html;
	}
	
	private String loadRawHtml(){
		String path=Data.findPath("?tax_server.html");
		String html=ReadWrite.readString(path);
		return html;
	}
	
	private String makeStats(){
		ByteBuilder sb=new ByteBuilder();
		
		if(!sketchOnly){
			sb.append("JGI taxonomy server stats:\n"
					+ "\nLast restarted "+startTime+"\n"
					+ "Running BBTools version "+Shared.BBTOOLS_VERSION_STRING+"\n");
		}else{
			sb.append("JGI"+(SketchObject.defaultParams.dbName==null ? "" : " "+SketchObject.defaultParams.dbName)+" sketch server stats:\n\n");
			
			if(domain!=null) {sb.append("Domain: "+domain+"\n");}
			if(SketchObject.useWhitelist()){
				sb.append("This server is running in whitelist mode; for best results, use local queries.\n");
				sb.append("Remote queries should specify a larger-than-normal sketch size.\n\n");
			}else if(SketchObject.blacklist()!=null){
				sb.append("This server is running in blacklist mode, using "+new File(SketchObject.blacklist()).getName()+".\n\n");
			}
			sb.append("Last restarted "+startTime+"\n");
			sb.append("Running BBTools version "+Shared.BBTOOLS_VERSION_STRING+"\n");
			sb.append("Settings: k="+SketchObject.k+(SketchObject.k2>0 ? ","+SketchObject.k2 : ""));
			if(SketchObject.amino){sb.append(" amino");}
			if(SketchObject.makeIndex){sb.append(" index");}
			if(SketchObject.useWhitelist()){sb.append(" whitelist");}
			if(SketchObject.blacklist()!=null){sb.append(" blacklist="+new File(SketchObject.blacklist()).getName());}
		}
		sb.nl().nl();
		sb.append(basicStats());
		if(sketchOnly){sb.append(makeExtendedStats());}
		
		return sb.toString();
	}
	
	public String makeExtendedStats(){
		ByteBuilder sb=new ByteBuilder();
		sb.append('\n');

		{
			sb.append("\nVersion\tCount\n");
			ArrayList<String> list=new ArrayList<String>();
			for(Entry<String, StringNum> e : versionMap.entrySet()){
				list.add(e.getValue().toString());
			}
			Collections.sort(list);
			for(String s : list){
				sb.append(s).append('\n');
			}
		}

		{
			sb.append("\nSketchs\tCount\tAvgTime\n");
			for(int i=0; i<timesByCount.length(); i++){
				double a=timesByCount.get(i)/1000000.0;
				long b=queryCounts.get(i);
				if(b>0){
					sb.append(i).append('\t').append(b).append('\t').append(a/b, 3).append('\n');
				}
			}
			sb.append('\n');
		}
		return sb.toString();
	}
	
	public String USAGE(String prefix){
		if(!countQueries){return prefix;}
		String basicStats=basicStats();
		return (prefix==null ? basicStats : prefix+"\n"+basicStats);
	}
	
	public String basicStats(){
		if(!countQueries){return "";}
		StringBuilder sb=new StringBuilder(500);
		
		final long uq=usageQueries.getAndIncrement();
		final long mq=malformedQueries.get();
		final long pt=plaintextQueries.get(), sc=semicolonQueries.get(), pa=pathQueries.get(), pp=printPathQueries.get(), ps=printSizeQueries.get();
		final long iq=internalQueries.get();
		final long lq=localQueries.get();
		final long rfq=refQueries.get();
		final long q=queries.get();
		final long nf=notFound.get();
		final double avgTimeDL=.000001*(elapsedTimeLocal.get()/(Tools.max(1.0, timeMeasurementsLocal.get())));//in milliseconds
		final double lastTimeDL=.000001*lastTimeLocal.get();
		final double avgTimeDR=.000001*(elapsedTimeRemote.get()/(Tools.max(1.0, timeMeasurementsRemote.get())));//in milliseconds
		final double lastTimeDR=.000001*lastTimeRemote.get();
		final double avgTimeDRF=.000001*(elapsedTimeReference.get()/(Tools.max(1.0, timeMeasurementsReference.get())));//in milliseconds
		final double lastTimeDRF=.000001*lastTimeReference.get();
		final double avgTimeDU=.000001*(elapsedTimeUsage.get()/(Tools.max(1.0, timeMeasurementsUsage.get())));//in milliseconds
		final double lastTimeDU=.000001*lastTimeUsage.get();
		final long exq=q-iq;
		final long rmq=q-lq;

		sb.append('\n').append("Queries:   ").append(q);
		sb.append('\n').append("Usage:     ").append(uq);
		if(sketchOnly){
			sb.append('\n').append("Invalid:   ").append(mq);
			sb.append('\n').append("Avg time:  ").append(Tools.format("%.3f ms (local queries)", avgTimeDL));
			sb.append('\n').append("Last time: ").append(Tools.format("%.3f ms (local queries)", lastTimeDL));
			sb.append('\n').append("Avg time:  ").append(Tools.format("%.3f ms (remote queries)", avgTimeDR));
			sb.append('\n').append("Last time: ").append(Tools.format("%.3f ms (remote queries)", lastTimeDR));
			sb.append('\n').append("Avg time:  ").append(Tools.format("%.3f ms (ref queries)", avgTimeDRF));
			sb.append('\n').append("Last time: ").append(Tools.format("%.3f ms (ref queries)", lastTimeDRF));
		}else{
			sb.append('\n').append("Avg time:  ").append(Tools.format("%.3f ms", avgTimeDR));
			sb.append('\n').append("Last time: ").append(Tools.format("%.3f ms", lastTimeDR));
			sb.append('\n').append("Avg time:  ").append(Tools.format("%.3f ms (usage queries)", avgTimeDU));
			sb.append('\n').append("Last time: ").append(Tools.format("%.3f ms (usage queries)", lastTimeDU));
		}
		sb.append('\n');
		sb.append('\n').append("Internal:  ").append(iq);
		sb.append('\n').append("External:  ").append(exq);
		if(!sketchOnly){sb.append('\n').append("NotFound:  ").append(nf);}
		sb.append('\n');
		
		if(sketchOnly){
			sb.append('\n').append("Local:     ").append(lq);
			sb.append('\n').append("Remote:    ").append(rmq);
			sb.append('\n').append("Reference: ").append(rfq);
			sb.append('\n');
			sb.append('\n').append("Depth:     ").append(depthQueries.get());
			sb.append('\n');
			sb.append('\n').append("Sketches:  ").append(querySketches.get());
			sb.append('\n').append("BytesIn:   ").append(bytesIn.get());
			sb.append('\n').append("BytesOut:  ").append(bytesOut.get());
			sb.append('\n');
			sb.append('\n').append("Single:    ").append(firstChunkSingle.get());
			sb.append('\n').append("Bulk:      ").append(firstChunkMulti.get());
			sb.append('\n').append("UnknownS:  ").append(unknownChunkSingle.get());
			sb.append('\n').append("UnknownB:  ").append(unknownChunkMulti.get());
			sb.append('\n').append("Total:     ").append(bulkCount.get());
		}else{
			sb.append('\n').append("gi:        ").append(giQueries.get());
			sb.append('\n').append("Name:      ").append(nameQueries.get());
			sb.append('\n').append("TaxID:     ").append(taxidQueries.get());
			sb.append('\n').append("Header:    ").append(headerQueries.get());
			sb.append('\n').append("Accession: ").append(accessionQueries.get());
			sb.append('\n').append("IMG:       ").append(imgQueries.get());
			sb.append('\n').append("Silva:     ").append(silvaHeaderQueries.get());
			sb.append('\n');
			sb.append('\n').append("Simple:    ").append(simpleQueries.get());
			sb.append('\n').append("Ancestor:  ").append(ancestorQueries.get());
			sb.append('\n').append("Children:  ").append(childrenQueries.get());
			sb.append('\n');
			sb.append('\n').append("Json:      ").append(q-pt-sc-pa);
			sb.append('\n').append("Plaintext: ").append(pt);
			sb.append('\n').append("Semicolon: ").append(sc);
			sb.append('\n').append("Path:      ").append(pa+pp);
			sb.append('\n').append("Size:      ").append(ps);
			sb.append('\n').append("Single:    ").append(firstChunkSingle.get());
			sb.append('\n').append("Bulk:      ").append(firstChunkMulti.get());
			sb.append('\n').append("Total:     ").append(bulkCount.get());
		}
		sb.append('\n');
		return sb.toString();
	}
	
	public boolean incrementQueries(HttpExchange t, boolean local, boolean refMode, boolean simple, boolean ancestor, 
			boolean plaintext, boolean semicolon, boolean path, boolean printChildren, boolean printPath, boolean printSize, int type){
		final boolean internal=ServerTools.isInternalQuery(t, addressPrefix, allowLocalHost, printIP, printHeaders);
		
		if(!countQueries){return internal;}
		queries.incrementAndGet();
		if(local){localQueries.incrementAndGet();}
		else if(refMode){localQueries.incrementAndGet();}
		
		if(type>=0){
			int type2=type&15;
			if(type2==GI){
				giQueries.incrementAndGet();
			}else if(type2==NAME){
				nameQueries.incrementAndGet();
			}else if(type2==TAXID){
				taxidQueries.incrementAndGet();
			}else if(type2==ACCESSION){
				accessionQueries.incrementAndGet();
			}else if(type2==IMG){
				imgQueries.incrementAndGet();
			}else if(type2==HEADER){
				headerQueries.incrementAndGet();
			}else if(type2==UNKNOWN){
				unknownQueries.incrementAndGet();
			}else if(type2==SILVAHEADER){
				silvaHeaderQueries.incrementAndGet();
			}
			
			if(simple){simpleQueries.incrementAndGet();}
			if(ancestor){ancestorQueries.incrementAndGet();}
			
			if(plaintext){plaintextQueries.incrementAndGet();}
			else if(semicolon){semicolonQueries.incrementAndGet();}
			else if(path){pathQueries.incrementAndGet();}
			
			if(printChildren){childrenQueries.incrementAndGet();}
			if(printPath){printPathQueries.incrementAndGet();}
			if(printSize){printSizeQueries.incrementAndGet();}
		}
		
		if(internal){internalQueries.incrementAndGet();}
		
		return internal;
	}
	
	/*--------------------------------------------------------------*/
	
	String compare(ArrayList<Sketch> inSketches, DisplayParams params){
		boolean success=true;
		final int inSize=inSketches.size();
		querySketches.addAndGet(inSize);
		if(Shared.threads()<2 || maxConcurrentSketchCompareThreads<2 || inSize<4){
			ByteBuilder sb=new ByteBuilder();
			success=searcher.compare(inSketches, sb, params, maxConcurrentSketchCompareThreads);
			return sb.toString();
		}else{//More sketches than threads, and more than one thread
			final int threads=Tools.min(maxConcurrentSketchCompareThreads, (inSize+4)/4);
			
			ByteBuilder[] out=new ByteBuilder[inSize];
			ArrayList<CompareThread> alct=new ArrayList<CompareThread>(threads);
			AtomicInteger next=new AtomicInteger(0);
			for(int i=0; i<threads; i++){
				alct.add(new CompareThread(inSketches, i, next, out, params));
			}
			for(CompareThread ct : alct){ct.start();}
			for(CompareThread ct : alct){

				//Wait until this thread has terminated
				while(ct.getState()!=Thread.State.TERMINATED){
					try {
						//Attempt a join operation
						ct.join();
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}

				synchronized(ct){
					success&=ct.success;
				}
			}
			alct=null;
			
			int len=0;
			for(ByteBuilder bb : out){len=len+bb.length;}
			ByteBuilder bb2=new ByteBuilder(len);
			for(int i=0; i<out.length; i++){
				ByteBuilder bb=out[i];
				bb2.append(bb);
				out[i]=null;
			}
			return bb2.toString();
		}
	}
	
	private class CompareThread extends Thread {
		
		CompareThread(final ArrayList<Sketch> inSketches_, final int tid_, final AtomicInteger nextSketch_, ByteBuilder[] out_, DisplayParams params_){
			inSketches=inSketches_;
			tid=tid_;
			nextSketch=nextSketch_;
			out=out_;
			params=params_;
		}
		
		@Override
		public void run(){
			success=false;
			final int inLim=inSketches.size();
			final boolean json=params.json();
			
			for(int inNum=nextSketch.getAndIncrement(); inNum<inLim; inNum=nextSketch.getAndIncrement()){
				Sketch a=inSketches.get(inNum);
				assert(buffer.cbs==null); //Because this sketch will only be used by one thread at a time, so per-buffer bitsets are not needed.
				SketchResults sr=searcher.processSketch(a, buffer, fakeID, map, params, 1);
				a.clearRefHitCounts();

				ByteBuilder bb=sr.toText(params);
				if(out!=null){
					if(json && inLim>1){
						if(inNum==0){
							bb.insert(0, (byte)'[');
						}
						if(inNum<inLim-1){
							bb.append(',');
						}else{
							bb.append(']');
						}
					}
					synchronized(out){
						out[inNum]=bb;
					}
				}
			}
			synchronized(this){success=true;}
		}
		
		private final ArrayList<Sketch> inSketches;
		private final int tid;
		private final CompareBuffer buffer=new CompareBuffer(false);
		private final DisplayParams params;
		private final ByteBuilder[] out;

		private final AtomicInteger nextSketch;
		private final AtomicInteger fakeID=new AtomicInteger(SketchObject.minFakeID);
		private ConcurrentHashMap<Integer, Comparison> map=new ConcurrentHashMap<Integer, Comparison>(101);
		
		boolean success=false;
		
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean sketchOnly=false;
	
	/*--------------------------------------------------------------*/
	/*----------------           Counters           ----------------*/
	/*--------------------------------------------------------------*/
	
	private HashMap<String, StringNum> versionMap=new HashMap<String, StringNum>();
	private AtomicLongArray timesByCount=new AtomicLongArray(10000);
	private AtomicLongArray queryCounts=new AtomicLongArray(10000);

	private AtomicLong notFound=new AtomicLong(0);
	private AtomicLong queries=new AtomicLong(0);
	/** Same IP address mask */
	private AtomicLong internalQueries=new AtomicLong(0);
	/** Local filesystem sketch */
	private AtomicLong localQueries=new AtomicLong(0);
	private AtomicLong refQueries=new AtomicLong(0);

	private AtomicLong depthQueries=new AtomicLong(0);
	
	private AtomicLong iconQueries=new AtomicLong(0);

	private AtomicLong querySketches=new AtomicLong(0);

	private AtomicLong unknownChunkSingle=new AtomicLong(0);
	private AtomicLong unknownChunkMulti=new AtomicLong(0);
	private AtomicLong firstChunkSingle=new AtomicLong(0);
	private AtomicLong firstChunkMulti=new AtomicLong(0);
	private AtomicLong nthChunkSingle=new AtomicLong(0);
	private AtomicLong nthChunkMulti=new AtomicLong(0);
	
	private AtomicLong singleQueries=new AtomicLong(0);
	private AtomicLong bulkQueries=new AtomicLong(0);
	private AtomicLong bulkCount=new AtomicLong(0);
	
	private AtomicLong giQueries=new AtomicLong(0);
	private AtomicLong nameQueries=new AtomicLong(0);
	private AtomicLong taxidQueries=new AtomicLong(0);
	private AtomicLong headerQueries=new AtomicLong(0);
	private AtomicLong accessionQueries=new AtomicLong(0);
	private AtomicLong imgQueries=new AtomicLong(0);
	private AtomicLong unknownQueries=new AtomicLong(0);
	private AtomicLong silvaHeaderQueries=new AtomicLong(0);
	
	private AtomicLong plaintextQueries=new AtomicLong(0);
	private AtomicLong semicolonQueries=new AtomicLong(0);
	private AtomicLong pathQueries=new AtomicLong(0);
	private AtomicLong printPathQueries=new AtomicLong(0);
	private AtomicLong printSizeQueries=new AtomicLong(0);
	private AtomicLong childrenQueries=new AtomicLong(0);
	
	private AtomicLong simpleQueries=new AtomicLong(0);
	private AtomicLong ancestorQueries=new AtomicLong(0);

	private AtomicLong usageQueries=new AtomicLong(0);
	private AtomicLong bytesIn=new AtomicLong(0);
	private AtomicLong bytesOut=new AtomicLong(0);

//	private AtomicLong elapsedTime=new AtomicLong(0);
//	private AtomicLong timeMeasurements=new AtomicLong(0);
//	private AtomicLong lastTime=new AtomicLong(0);
	
	private AtomicLong elapsedTimeUsage=new AtomicLong(0);
	private AtomicLong timeMeasurementsUsage=new AtomicLong(0);
	private AtomicLong lastTimeUsage=new AtomicLong(0);
	
	private AtomicLong elapsedTimeRemote=new AtomicLong(0);
	private AtomicLong timeMeasurementsRemote=new AtomicLong(0);
	private AtomicLong lastTimeRemote=new AtomicLong(0);
	
	private AtomicLong elapsedTimeLocal=new AtomicLong(0);
	private AtomicLong timeMeasurementsLocal=new AtomicLong(0);
	private AtomicLong lastTimeLocal=new AtomicLong(0);
	
	private AtomicLong elapsedTimeReference=new AtomicLong(0);
	private AtomicLong timeMeasurementsReference=new AtomicLong(0);
	private AtomicLong lastTimeReference=new AtomicLong(0);

	private AtomicLong malformedQueries=new AtomicLong(0);
	
	/*--------------------------------------------------------------*/
	/*----------------            Params            ----------------*/
	/*--------------------------------------------------------------*/

	public boolean printIP=false;
	public boolean printHeaders=false;
	public boolean countQueries=true;
	public float prealloc=0;
	public boolean useHtml=false;

	/** Location of GiTable file */
	private String giTableFile=null;
	/** Location of TaxTree file */
	private String taxTreeFile="auto";
	/** Comma-delimited locations of Accession files */
	private String accessionFile=null;
	/** Location of IMG dump file */
	private String imgFile=null;
	/** Location of accession pattern file */
	private String patternFile=null;
	
	private String sizeFile=null;

	/** Location of sequence directory tree */
	private String basePath="/global/cfs/cdirs/bbtools/tree/";
	
	/** Used for taxonomic tree traversal */
	private final TaxTree tree;
	
	/** Maps URL Strings to numeric query types */
	private final HashMap<String, Integer> typeMap;
	/** Maps common organism names to scientific names */
	private final HashMap<String, String> commonMap;
	
	/** Hash taxonomic names for lookup */
	private boolean hashNames=true;
	private boolean hashDotFormat=true;

	/** Kill code of prior server instance (optional) */
	private String oldKillCode=null;
	/** Address of prior server instance (optional) */
	private String oldAddress=null;

	/** Address of current server instance (optional) */
	public String domain=null;

	public int maxConcurrentSketchCompareThreads=8;//TODO: This might be too high when lots of concurrent sessions are active
	public int maxConcurrentSketchLoadThreads=4;//TODO: This might be too high when lots of concurrent sessions are active
	public int handlerThreads=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private final boolean distributed;
	private final int serverNum;
	private final int serverCount;
	private ArrayList<String> slaveAddress;
	
	public final String favIconPath=Data.findPath("?favicon.ico");
	public final byte[] favIcon=ReadWrite.readRaw(favIconPath);
	
	private final String startTime=new Date().toString();
	
	/** Listen on this port */
	public final int port;
	/** Code to validate kill requests */
	public final String killCode;
	
	public final HttpServer httpServer;

	/** Bit to set for plaintext query types */
	public static final int PT_BIT=16;
	/** Bit to set for semicolon-delimited query types */
	public static final int SC_BIT=32;
	/** Bit to set for path query types */
	public static final int PA_BIT=64;
	/** Request query types */
	public static final int UNKNOWN=0, GI=1, NAME=2, TAXID=3, HEADER=4, ACCESSION=5, IMG=6, SILVAHEADER=7;
	/** Plaintext-response query types */
	public static final int PT_GI=GI+PT_BIT, PT_NAME=NAME+PT_BIT, PT_TAXID=TAXID+PT_BIT,
			PT_HEADER=HEADER+PT_BIT, PT_ACCESSION=ACCESSION+PT_BIT, PT_IMG=IMG+PT_BIT, PT_SILVAHEADER=SILVAHEADER+PT_BIT;
	/** Semicolon-response query types */
	public static final int SC_GI=GI+SC_BIT, SC_NAME=NAME+SC_BIT, SC_TAXID=TAXID+SC_BIT,
			SC_HEADER=HEADER+SC_BIT, SC_ACCESSION=ACCESSION+SC_BIT, SC_IMG=IMG+SC_BIT, SC_SILVAHEADER=SILVAHEADER+PT_BIT;
	
	public static final int SOURCE_REFSEQ=1, SOURCE_SILVA=2, SOURCE_IMG=3;
	
	/** Generic response when asking for tax advice */
	public static final String TAX_ADVICE="This site does not give tax advice.";
	/** Generic response for incorrect kill code */
	public static final String BAD_CODE="Incorrect code.";
	/** Generic response for badly-formatted queries */
	public final String USAGE;
	/** HTML version */
//	public final String USAGE_HTML;
	public final String rawHtml;
	
	/** Tool for comparing query sketches to reference sketches */
	public final SketchSearcher searcher=new SketchSearcher();

	public final boolean hasSketches;

	final boolean allowRemoteFileAccess;
	final boolean allowLocalHost;
	final String addressPrefix;
	private boolean clearMem=true;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false, verbose2=false, logUsage=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	
}
