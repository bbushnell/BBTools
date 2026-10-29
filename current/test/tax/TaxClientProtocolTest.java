package test.tax;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.io.PrintStream;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLDecoder;
import java.net.URLStreamHandler;
import java.net.URLStreamHandlerFactory;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import shared.Shared;
import tax.Query;
import tax.TaxClient;

/** Tests taxonomy client routing and CLI payload preservation without opening sockets.
 * A JVM-local protocol handler supplies deterministic responses and records URLs.
 * @author Brian Bushnell, Shinobu */
public final class TaxClientProtocolTest{

	/** Runs one isolated baseline/fixed client contract suite.
	 * @param args Optional expect=baseline/fixed and batch=baseline/fixed
	 * @throws Exception If the isolated fixture cannot capture client output */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		Boolean batchBaseline=null;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else if(arg.equals("batch=baseline")){batchBaseline=Boolean.TRUE;}
			else if(arg.equals("batch=fixed")){batchBaseline=Boolean.FALSE;}
			else{throw new IllegalArgumentException(arg);}
		}
		URL.setURLStreamHandlerFactory(new URLStreamHandlerFactory(){
			@Override
			public URLStreamHandler createURLStreamHandler(String protocol){
				if(!protocol.equals("taxfixture")){return null;}
				return new URLStreamHandler(){
					@Override
					protected URLConnection openConnection(URL url){return new Reply(url);}
				};
			}
		});
		Shared.setTaxServer("taxfixture://local/");
		int header=TaxClient.headerToTaxid("tid|4 read0");
		check(header==(baseline ? 2 : 4), "Single header must use the header route after repair");
		int[] array=TaxClient.headerToTaxidArray("tid|4 read0");
		check(array.length==1 && array[0]==4, "Existing array-header route is the unaffected control");
		for(boolean slow : new boolean[]{true, false}){
			int before=opened;
			String output=runMain("header=tid|4 read=extra", "slow="+slow);
			if(baseline){check(opened==before, "Baseline CLI must expose its ignored header list");}
			else{
				check(opened==before+1, "CLI must submit each header in both modes");
				check(output.contains("tid|4 read=extra\t4"), "CLI must retain the complete header value: "+output);
			}
		}
		runMain("name=Alpha=Beta", "slow=t");
		check(lastPath.equals(baseline ? "/pt/name/Alpha" : "/pt/name/Alpha=Beta"),
				"CLI name payload must preserve equals signs: "+lastPath);
		int emptyFailures=0;
		for(int kind=0; kind<4; kind++){
			int before=opened;
			try{
				int[] result;
				if(kind==0){result=TaxClient.accessionToTaxidArray(new ArrayList<String>());}
				else if(kind==1){result=TaxClient.giToTaxidArray(new ArrayList<Integer>());}
				else if(kind==2){result=TaxClient.nameToTaxidArray(new ArrayList<String>());}
				else{result=TaxClient.headerToTaxidArray(new ArrayList<String>());}
				check(!baseline && result==null, "Empty query list must return the no-result sentinel");
			}catch(NullPointerException e){
				if(!baseline){throw e;}
				emptyFailures++;
			}
			check(opened==before, "Empty lists must not make requests");
		}
		check(emptyFailures==(baseline ? 4 : 0), "Unexpected empty-list exception count");
		int[] gi=Query.getGi("7");
		check(gi.length==1 && gi[0]==2, "Legacy response parser control failed");
		// URL parsers may retain or discard the embedded tab; neither creates a path delimiter.
		check(baseline ? lastPath.equals("/pt_gi\t7") || lastPath.equals("/pt_gi7") : lastPath.equals("/pt_gi/7"),
				"Legacy query path separator: "+lastPath);
		if(batchBaseline!=null){checkShortBatches(batchBaseline.booleanValue());}
		check(opened==closed, "All acquired response streams must be closed");
		System.out.println("PASS client_protocol baseline="+baseline+" requests="+opened+" empty_failures="+emptyFailures);
	}

	/** Prevents construction of this standalone fixture. */
	private TaxClientProtocolTest(){}

	/** Checks complete two-item replies and truncated replies for every fast CLI route.
	 * @param baseline Whether a truncated reply should expose the old silent omission
	 * @throws Exception If CLI output cannot be captured */
	private static void checkShortBatches(final boolean baseline) throws Exception{
		final String[][] cases={
			{"name=Alpha,Beta", "Alpha", "Beta", "2", "3"},
			{"gi=7,8", "7", "8", "2", "3"},
			{"accession=AAA,BBB", "AAA", "BBB", "2", "3"},
			{"header=tid|4 one,tid|5 two", "tid|4 one", "tid|5 two", "4", "5"}
		};
		try{
			for(String[] c : cases){
				final String kind=c[0].substring(0, c[0].indexOf('='));
				forcedResponse=c[3]+","+c[4];
				String output=runMain(c[0], "slow=f");
				check(countLine(output, c[1]+"\t"+c[3])==1 &&
						countLine(output, c[2]+"\t"+c[4])==1,
						"Complete "+kind+" reply must emit both requested mappings: "+output);
				forcedResponse=c[3];
				Throwable failure=null;
				try{output=runMain(c[0], "slow=f");}
				catch(RuntimeException|Error e){failure=e; output=lastCliOutput;}
				if(baseline){
					check(failure==null, "Baseline short "+kind+" reply should return normally");
					check(countLine(output, c[1]+"\t"+c[3])==1 &&
							countLine(output, c[2]+"\t"+c[4])==0,
							"Baseline must expose the omitted second "+kind+" result: "+output);
				}else{
					check(failure instanceof RuntimeException,
							"Short "+kind+" reply must fail clearly: "+failure);
					check(failure.getMessage()!=null &&
							failure.getMessage().startsWith("Expected 2 taxonomy IDs") &&
							failure.getMessage().endsWith("received 1."),
							"Failure must identify expected and actual counts: "+failure);
					check(countLine(output, c[1]+"\t"+c[3])==0,
							"Reject a short "+kind+" batch before emitting partial mappings: "+output);
				}
				System.out.println("CLIENT_BATCH_RECEIPT kind="+kind+
						" expected=2 normal=2 short_received=1 rejected="+!baseline);
			}
		}finally{
			forcedResponse=null;
		}
	}

	/** Counts exact output lines without treating diagnostic text as result records.
	 * @param output Captured CLI output
	 * @param expected Complete result line to count
	 * @return Number of exact matches */
	private static int countLine(final String output, final String expected){
		int count=0;
		for(String line : output.split("\\r?\\n")){
			if(line.equals(expected)){count++;}
		}
		return count;
	}

	/** Captures the CLI's diagnostic/result streams in this isolated fixture JVM.
	 * @param args Client arguments
	 * @return Captured UTF-8 output
	 * @throws Exception If stream construction fails */
	private static String runMain(String... args) throws Exception{
		PrintStream originalErr=System.err, originalOut=System.out, originalClient=TaxClient.outstream;
		ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		PrintStream capture=new PrintStream(bytes, true, "UTF-8");
		try{
			System.setErr(capture);
			System.setOut(capture);
			TaxClient.main(args);
			capture.flush();
			return new String(bytes.toByteArray(), StandardCharsets.UTF_8);
		}finally{
			capture.flush();
			lastCliOutput=new String(bytes.toByteArray(), StandardCharsets.UTF_8);
			System.setErr(originalErr);
			System.setOut(originalOut);
			TaxClient.outstream=originalClient;
			capture.close();
		}
	}

	/** Checks stay active under -da.
	 * @param condition Required result
	 * @param message Failure explanation */
	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** In-memory URL connection; no DNS, listener, HTTP client or network is involved. */
	private static final class Reply extends URLConnection{
		/** Binds this response to the requested URL.
		 * @param url Requested fixture URL */
		Reply(URL url){super(url);}
		@Override
		public void connect(){}
		@Override
		public InputStream getInputStream(){
			try{lastPath=URLDecoder.decode(url.getPath(), "UTF-8");}
			catch(java.io.UnsupportedEncodingException e){throw new AssertionError(e);}
			opened++;
			String response=forcedResponse!=null ? forcedResponse : lastPath.startsWith("/pt/header/") ? "4" : "2";
			byte[] body=response.getBytes(StandardCharsets.UTF_8);
			return new ByteArrayInputStream(body){
				@Override
				public void close(){closed++;}
			};
		}
	}

	/** Recorded decoded path of the most recent in-memory request. */
	private static String lastPath;
	/** Optional deterministic response for the batch-length cases. */
	private static String forcedResponse;
	/** Most recent captured CLI output, retained even when main throws. */
	private static String lastCliOutput;
	/** Counts response streams opened and closed by production client code. */
	private static int opened, closed;
}
