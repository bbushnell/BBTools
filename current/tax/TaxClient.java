package tax;

import java.io.PrintStream;
import java.util.ArrayList;

import parse.Parse;
import parse.PreParser;
import server.PercentEncoding;
import server.ServerTools;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.StringNum;

/**
 * Blocking client for resolving identifiers through the configured taxonomy server.
 * Single and comma-separated batch queries return integer taxonomy IDs; negative
 * entries are server-provided missing-result sentinels. Short requests use encoded
 * GET paths; requests of at least 2000 message characters use POST.
 * Null or empty transport responses trigger up to twelve attempts with backoff.
 * Methods allocate request/response storage and consult global server configuration;
 * the CLI also mutates the shared diagnostic stream and verbosity flag.
 * @author Brian Bushnell
 */
public class TaxClient{

	/** Constructs an otherwise stateless client; lookup methods are static. */
	public TaxClient(){}

	/**
	 * Program entry point for taxonomic ID lookup queries.
	 * Processes names, GI numbers, accessions and headers in separate groups, writing each
	 * input followed by its taxonomy ID to the preparser-selected diagnostic stream.
	 * slow=t sends individual requests; slow=f sends one batch per identifier type.
	 * Fast mode validates the response count before printing that identifier group.
	 * @param args Command-line arguments specifying query parameters and options
	 * @throws IllegalStateException If a fast-query group receives no result or the wrong count
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null, false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ArrayList<Integer> gi=new ArrayList<Integer>();
		ArrayList<String> acc=new ArrayList<String>();
		ArrayList<String> name=new ArrayList<String>();
		ArrayList<String> header=new ArrayList<String>();

		boolean slow=true;

		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase();
			String b=split.length>1 && !split[1].isEmpty() ? split[1] : null;

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("name") || a.equals("names")){
				if(b!=null){
					for(String s : b.split(",")){
						name.add(s);
					}
				}
			}else if(a.equals("gi")){
				if(b!=null){
					for(String s : b.split(",")){
						gi.add(Integer.parseInt(s));
					}
				}
			}else if(a.equals("accession")){
				if(b!=null){
					for(String s : b.split(",")){
						acc.add(s);
					}
				}
			}else if(a.equals("header")){
				if(b!=null){
					for(String s : b.split(",")){
						header.add(s);
					}
				}
			}else if(a.equals("slow")){
				slow=Parse.parseBoolean(b);
			}else if(a.equals("fast")){
				slow=!Parse.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Legacy template option: the parsed value has no effect.
				long ignoredValue=Parse.parseKMG(b);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		if(slow){
			for(String s : name){
				outstream.println(s+"\t"+nameToTaxid(s));
			}
			for(Integer s : gi){
				outstream.println(s+"\t"+giToTaxid(s));
			}
			for(String s : acc){
				outstream.println(s+"\t"+accessionToTaxid(s));
			}
			for(String s : header){
				outstream.println(s+"\t"+headerToTaxid(s));
			}
		}else{
			if(name.size()>0){
				int[] ids=nameToTaxidArray(name);
				requireBatchCount(ids, name.size(), "name");
				for(int i=0; i<ids.length; i++){
					outstream.println(name.get(i)+"\t"+ids[i]);
				}
			}
			if(gi.size()>0){
				int[] ids=giToTaxidArray(gi);
				requireBatchCount(ids, gi.size(), "gi");
				for(int i=0; i<ids.length; i++){
					outstream.println(gi.get(i)+"\t"+ids[i]);
				}
			}
			if(acc.size()>0){
				int[] ids=accessionToTaxidArray(acc);
				requireBatchCount(ids, acc.size(), "accession");
				for(int i=0; i<ids.length; i++){
					outstream.println(acc.get(i)+"\t"+ids[i]);
				}
			}
			if(header.size()>0){
				int[] ids=headerToTaxidArray(header);
				requireBatchCount(ids, header.size(), "header");
				for(int i=0; i<ids.length; i++){
					outstream.println(header.get(i)+"\t"+ids[i]);
				}
			}
		}

		t.stopAndPrint();
	}

	/** Rejects incomplete or oversized replies before associating IDs with queries (TAX-072).
	 * @param ids Parsed server result, or null for an unusable response
	 * @param expected Positive number of queries in this CLI group
	 * @param kind Identifier type for diagnostics
	 * @throws IllegalStateException If the response does not contain exactly one ID per query */
	private static void requireBatchCount(final int[] ids, final int expected, final String kind){
		assert(expected>0) : "Only nonempty CLI groups reach this check (main size guards); expected="+expected;
		if(ids==null || ids.length!=expected){
			final String received=ids==null ? "no response" : Integer.toString(ids.length);
			throw new IllegalStateException("Expected "+expected+" taxonomy IDs for "+kind+
					" queries, received "+received+".");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Single Queries       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Converts a single accession ID to taxonomic ID using default server.
	 * @param accession The accession identifier to look up
	 * @return Parsed taxID, or -1 if the response is empty, null or lacks an initial digit/sign
	 * @throws NumberFormatException If a numeric-leading response is not a valid int
	 */
	public static int accessionToTaxid(String accession){
		String s=sendAndReceive("pt/accession/", accession);
		if(s==null || s.length()<1 || !Tools.isDigitOrSign(s.charAt(0))){return -1;}
		return Integer.parseInt(s);
	}

	/**
	 * Converts a single accession ID to taxonomic ID using specified server.
	 * @param path Base server URL, including any required trailing slash
	 * @param accession The accession identifier to look up
	 * @return Parsed taxID, or -1 if the response is empty, null or lacks an initial digit/sign
	 * @throws NumberFormatException If a numeric-leading response is not a valid int
	 */
	public static int accessionToTaxidSpecificServer(String path, String accession){
		String s=sendAndReceiveSpecificServer(path, "pt/accession/", accession);
		if(s==null || s.length()<1 || !Tools.isDigitOrSign(s.charAt(0))){return -1;}
		return Integer.parseInt(s);
	}

	/**
	 * Converts a GI number to taxonomic ID.
	 * @param gi The GI number to look up
	 * @return Parsed taxID, or -1 if the response is empty, null or lacks an initial digit/sign
	 * @throws NumberFormatException If a numeric-leading response is not a valid int
	 */
	public static int giToTaxid(int gi){
		String s=sendAndReceive("pt/gi/", Integer.toString(gi));
		if(s==null || s.length()<1 || !Tools.isDigitOrSign(s.charAt(0))){return -1;}
		return Integer.parseInt(s);
	}

	/**
	 * Converts an organism name to taxonomic ID.
	 * Spaces in the name are replaced with underscores before lookup.
	 * @param name The organism name to look up
	 * @return Parsed taxID, or -1 if the response is empty, null or lacks an initial digit/sign
	 * @throws NumberFormatException If a numeric-leading response is not a valid int
	 */
	public static int nameToTaxid(String name){
		String s=sendAndReceive("pt/name/", name.replace(' ', '_'));
		if(s==null || s.length()<1 || !Tools.isDigitOrSign(s.charAt(0))){return -1;}
		return Integer.parseInt(s);
	}

	/**
	 * Converts a sequence header to taxonomic ID.
	 * @param header The sequence header to look up
	 * @return Parsed taxID, or -1 if the response is empty, null or lacks an initial digit/sign
	 * @throws NumberFormatException If a numeric-leading response is not a valid int
	 */
	public static int headerToTaxid(String header){
		String s=sendAndReceive("pt/header/", header);
		if(s==null || s.length()<1 || !Tools.isDigitOrSign(s.charAt(0))){return -1;}
		return Integer.parseInt(s);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Batch Queries       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Converts comma-separated accession IDs to array of taxonomic IDs.
	 * @param accession Comma-separated accession identifiers
	 * @return Parsed IDs, or null if the response is empty, null or lacks an initial digit/sign
	 * @throws NumberFormatException If a response field is malformed or out of range
	 */
	public static int[] accessionToTaxidArray(String accession){
		String s=sendAndReceive("pt/accession/", accession);
		return splitOutput(s);
	}

	/**
	 * Converts comma-separated GI numbers to array of taxonomic IDs.
	 * @param gi Comma-separated GI numbers as string
	 * @return Parsed IDs, or null if the response is empty, null or lacks an initial digit/sign
	 * @throws NumberFormatException If a response field is malformed or out of range
	 */
	public static int[] giToTaxidArray(String gi){
		String s=sendAndReceive("pt/gi/", gi);
		return splitOutput(s);
	}

	/**
	 * Converts comma-separated organism names to array of taxonomic IDs.
	 * Spaces in names are replaced with underscores before lookup.
	 * @param name Comma-separated organism names
	 * @return Parsed IDs, or null if the response is empty, null or lacks an initial digit/sign
	 * @throws NumberFormatException If a response field is malformed or out of range
	 */
	public static int[] nameToTaxidArray(String name){
		String s=sendAndReceive("pt/name/", name.replace(' ', '_'));
		return splitOutput(s);
	}

	/**
	 * Converts comma-separated sequence headers to array of taxonomic IDs.
	 * @param header Comma-separated sequence headers
	 * @return Parsed IDs, or null if the response is empty, null or lacks an initial digit/sign
	 * @throws NumberFormatException If a response field is malformed or out of range
	 */
	public static int[] headerToTaxidArray(String header){
		String s=sendAndReceive("pt/header/", header);
		return splitOutput(s);
	}


	/**
	 * Converts list of accession IDs to array of taxonomic IDs.
	 * Null or empty lists return null without making a request.
	 * @param accession List of accession identifiers, or null
	 * @return Parsed IDs, or null if the response is empty, null or lacks an initial digit/sign
	 * @throws NumberFormatException If a response field is malformed or out of range
	 */
	public static int[] accessionToTaxidArray(ArrayList<String> accession){
		if(accession==null || accession.isEmpty()){return null;}
		String s=sendAndReceive("pt/accession/", fuse(accession));
		return splitOutput(s);
	}

	/**
	 * Converts list of GI numbers to array of taxonomic IDs.
	 * Null or empty lists return null without making a request.
	 * @param gi List of GI numbers, or null
	 * @return Parsed IDs, or null if the response is empty, null or lacks an initial digit/sign
	 * @throws NumberFormatException If a response field is malformed or out of range
	 */
	public static int[] giToTaxidArray(ArrayList<Integer> gi){
		if(gi==null || gi.isEmpty()){return null;}
		String s=sendAndReceive("pt/gi/", fuse(gi));
		return splitOutput(s);
	}

	/**
	 * Converts list of organism names to array of taxonomic IDs.
	 * Spaces in names are replaced with underscores before lookup.
	 * Null or empty lists return null without making a request.
	 * @param name List of organism names, or null
	 * @return Parsed IDs, or null if the response is empty, null or lacks an initial digit/sign
	 * @throws NumberFormatException If a response field is malformed or out of range
	 */
	public static int[] nameToTaxidArray(ArrayList<String> name){
		if(name==null || name.isEmpty()){return null;}
		String s=sendAndReceive("pt/name/", fuse(name).replace(' ', '_'));
		return splitOutput(s);
	}

	/**
	 * Converts list of sequence headers to array of taxonomic IDs.
	 * Null or empty lists return null without making a request.
	 * @param header List of sequence headers, or null
	 * @return Parsed IDs, or null if the response is empty, null or lacks an initial digit/sign
	 * @throws NumberFormatException If a response field is malformed or out of range
	 */
	public static int[] headerToTaxidArray(ArrayList<String> header){
		if(header==null || header.isEmpty()){return null;}
		String s=sendAndReceive("pt/header/", fuse(header));
		return splitOutput(s);
	}

	/**
	 * Parses comma-separated taxonomic IDs, dropping trailing empty fields.
	 * Only the first character is screened before integer parsing; malformed later
	 * fields and integer overflow throw rather than returning a missing-result array.
	 * @param s Server response containing comma-separated taxonomic IDs
	 * @return Parsed IDs, or null if the response is empty, null or lacks an initial digit/sign
	 * @throws NumberFormatException If any numeric field is malformed or out of range
	 */
	private static final int[] splitOutput(String s){
		if(s==null || s.length()<1 || !Tools.isDigitOrSign(s.charAt(0))){return null;}
		String[] split=s.split(",");
		int[] ret=new int[split.length];
		for(int i=0; i<split.length; i++){
			ret[i]=Integer.parseInt(split[i]);
		}
		return ret;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Transport          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Sends query to default taxonomy server and receives response.
	 * @param prefix URL prefix for the query type
	 * @param message Query message to send
	 * @return Server response string, or null if communication failed
	 */
	private static final String sendAndReceive(String prefix, String message){
		String path=Shared.taxServer();
		return sendAndReceiveSpecificServer(path, prefix, message);
	}

	/**
	 * Sends query to specified taxonomy server with retry logic.
	 * Makes at most 12 attempts while the response is null or empty. There are at
	 * most 11 sleeps, starting at 200ms and doubling after each retry. Interrupted
	 * sleeps currently print a stack trace and continue the retry loop.
	 *
	 * @param path Server URL path
	 * @param prefix URL prefix for the query type
	 * @param message Query message to send
	 * @return Last response, possibly null or empty when retries are exhausted
	 */
	private static final String sendAndReceiveSpecificServer(String path, String prefix, String message){
		String response=null;
		for(int i=0, millis=200; i<12 && (response==null || response.length()==0); i++){
			if(i>0){
				System.err.println("Null response, retrying after "+millis+" ms.");
				try{
					Thread.sleep(millis);
				}catch(InterruptedException e){
					e.printStackTrace();
				}
				millis*=2;
			}
			response=sendAndReceiveOnceSpecificServer(path, prefix, message);
		}
		return response;
	}

	/**
	 * Sends single query to default taxonomy server without retries.
	 * @param prefix URL prefix for the query type
	 * @param message Query message to send
	 * @return Server response string, or null if communication failed
	 */
	private static String sendAndReceiveOnce(String prefix, String message){
		String path=Shared.taxServer();
		return sendAndReceiveOnceSpecificServer(path, prefix, message);
	}

	/**
	 * Sends single query to specified server using GET or POST based on message length.
	 * Uses GET below 2000 message characters before percent encoding, and POST
	 * otherwise. POST sends prefix+message bytes to the base URL followed by $POST.
	 *
	 * @param path Server URL path
	 * @param prefix URL prefix for the query type
	 * @param message Non-null query text
	 * @return Server response string, or null if communication failed
	 */
	private static String sendAndReceiveOnceSpecificServer(String path, String prefix, String message){
		final String response;
		if(message.length()<2000){//NERSC Apache limit is around 8kb
			message=prefix+PercentEncoding.symbolToCode(message);
			ByteBuilder bb=ServerTools.readPage(path+message, false);
			response=(bb==null ? null : bb.toString());
		}else{
			StringNum sn=ServerTools.sendAndReceive((prefix+message).getBytes(), path+"$POST");
			response=sn.s;
		}
		return response;
	}

	/**
	 * Concatenates ArrayList elements into comma-separated string.
	 * @param list ArrayList of objects to concatenate
	 * @return Comma-separated string representation, or null if list is null/empty
	 */
	private static String fuse(ArrayList<?> list){
		//FIXED [tax/TaxClient#001]: size()<0 never rejected an empty list,
		//so removing the final comma called setLength(-1). Return null instead;
		//public list-query methods also skip empty requests before transport.
		if(list==null || list.size()<1){return null;}
		StringBuilder sb=new StringBuilder();
		for(Object s : list){
			sb.append(s).append(',');
		}
		sb.setLength(sb.length()-1);
		return sb.toString();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Shared CLI result/diagnostic stream, initially System.err. */
	public static PrintStream outstream=System.err;
	/** Parsed CLI flag; currently unused by this class's request logic. */
	public static boolean verbose=false;

}
