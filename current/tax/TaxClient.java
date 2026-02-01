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
 * Client for querying taxonomic information from remote taxonomy server.
 * Converts GI numbers, accession IDs, and organism names to taxonomic IDs.
 * Supports both single queries and batch processing for efficiency.
 * @author Brian Bushnell
 */
public class TaxClient {
	
	/**
	 * Program entry point for taxonomic ID lookup queries.
	 * Parses command-line arguments and processes taxonomic queries in batch or individual mode.
	 * @param args Command-line arguments specifying query parameters and options
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
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
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
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
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
		}else{
			if(name.size()>0){
				int[] ids=nameToTaxidArray(name);
				for(int i=0; i<ids.length; i++){
					outstream.println(name.get(i)+"\t"+ids[i]);
				}
			}
			if(gi.size()>0){
				int[] ids=giToTaxidArray(gi);
				for(int i=0; i<ids.length; i++){
					outstream.println(gi.get(i)+"\t"+ids[i]);
				}
			}
			if(acc.size()>0){
				int[] ids=accessionToTaxidArray(acc);
				for(int i=0; i<ids.length; i++){
					outstream.println(acc.get(i)+"\t"+ids[i]);
				}
			}
		}
		
		t.stopAndPrint();
	}
	
	/**
	 * Converts a single accession ID to taxonomic ID using default server.
	 * @param accession The accession identifier to look up
	 * @return Taxonomic ID, or -1 if not found or invalid response
	 */
	public static int accessionToTaxid(String accession){
		String s=sendAndReceive("pt/accession/",accession);
		if(s==null || s.length()<1 || !Tools.isDigitOrSign(s.charAt(0))){return -1;}
		return Integer.parseInt(s);
	}
	
	/**
	 * Converts a single accession ID to taxonomic ID using specified server.
	 * @param path Server path/URL to query
	 * @param accession The accession identifier to look up
	 * @return Taxonomic ID, or -1 if not found or invalid response
	 */
	public static int accessionToTaxidSpecificServer(String path, String accession){
		String s=sendAndReceiveSpecificServer(path,"pt/accession/",accession);
		if(s==null || s.length()<1 || !Tools.isDigitOrSign(s.charAt(0))){return -1;}
		return Integer.parseInt(s);
	}
	
	/**
	 * Converts a GI number to taxonomic ID.
	 * @param gi The GI number to look up
	 * @return Taxonomic ID, or -1 if not found or invalid response
	 */
	public static int giToTaxid(int gi){
		String s=sendAndReceive("pt/gi/",Integer.toString(gi));
		if(s==null || s.length()<1 || !Tools.isDigitOrSign(s.charAt(0))){return -1;}
		return Integer.parseInt(s);
	}
	
	/**
	 * Converts an organism name to taxonomic ID.
	 * Spaces in the name are replaced with underscores before lookup.
	 * @param name The organism name to look up
	 * @return Taxonomic ID, or -1 if not found or invalid response
	 */
	public static int nameToTaxid(String name){
		String s=sendAndReceive("pt/name/",name.replace(' ', '_'));
		if(s==null || s.length()<1 || !Tools.isDigitOrSign(s.charAt(0))){return -1;}
		return Integer.parseInt(s);
	}
	
	/**
	 * Converts a sequence header to taxonomic ID.
	 * @param header The sequence header to look up
	 * @return Taxonomic ID, or -1 if not found or invalid response
	 */
	public static int headerToTaxid(String header){
		String s=sendAndReceive("pt/name/",header);
		if(s==null || s.length()<1 || !Tools.isDigitOrSign(s.charAt(0))){return -1;}
		return Integer.parseInt(s);
	}
	
	
	/**
	 * Converts comma-separated accession IDs to array of taxonomic IDs.
	 * @param accession Comma-separated accession identifiers
	 * @return Array of taxonomic IDs, or null if invalid response
	 */
	public static int[] accessionToTaxidArray(String accession){
		String s=sendAndReceive("pt/accession/",accession);
//		System.err.println("Sent "+accession+"\nReceived "+s);
		return splitOutput(s);
	}
	
	/**
	 * Converts comma-separated GI numbers to array of taxonomic IDs.
	 * @param gi Comma-separated GI numbers as string
	 * @return Array of taxonomic IDs, or null if invalid response
	 */
	public static int[] giToTaxidArray(String gi){
		String s=sendAndReceive("pt/gi/",gi);
		return splitOutput(s);
	}
	
	/**
	 * Converts comma-separated organism names to array of taxonomic IDs.
	 * Spaces in names are replaced with underscores before lookup.
	 * @param name Comma-separated organism names
	 * @return Array of taxonomic IDs, or null if invalid response
	 */
	public static int[] nameToTaxidArray(String name){
		String s=sendAndReceive("pt/name/",name.replace(' ', '_'));
		return splitOutput(s);
	}
	
	/**
	 * Converts comma-separated sequence headers to array of taxonomic IDs.
	 * @param header Comma-separated sequence headers
	 * @return Array of taxonomic IDs, or null if invalid response
	 */
	public static int[] headerToTaxidArray(String header){
		String s=sendAndReceive("pt/header/",header);
		return splitOutput(s);
	}
	
	
	/**
	 * Converts list of accession IDs to array of taxonomic IDs.
	 * @param accession List of accession identifiers
	 * @return Array of taxonomic IDs, or null if invalid response
	 */
	public static int[] accessionToTaxidArray(ArrayList<String> accession){
		String s=sendAndReceive("pt/accession/",fuse(accession));
		return splitOutput(s);
	}
	
	/**
	 * Converts list of GI numbers to array of taxonomic IDs.
	 * @param gi List of GI numbers
	 * @return Array of taxonomic IDs, or null if invalid response
	 */
	public static int[] giToTaxidArray(ArrayList<Integer> gi){
		String s=sendAndReceive("pt/gi/",fuse(gi));
		return splitOutput(s);
	}
	
	/**
	 * Converts list of organism names to array of taxonomic IDs.
	 * Spaces in names are replaced with underscores before lookup.
	 * @param name List of organism names
	 * @return Array of taxonomic IDs, or null if invalid response
	 */
	public static int[] nameToTaxidArray(ArrayList<String> name){
		String s=sendAndReceive("pt/name/",fuse(name).replace(' ', '_'));
		return splitOutput(s);
	}
	
	/**
	 * Converts list of sequence headers to array of taxonomic IDs.
	 * @param header List of sequence headers
	 * @return Array of taxonomic IDs, or null if invalid response
	 */
	public static int[] headerToTaxidArray(ArrayList<String> header){
		String s=sendAndReceive("pt/header/",fuse(header));
		return splitOutput(s);
	}
	
	/**
	 * Parses comma-separated taxonomic IDs from server response.
	 * @param s Server response containing comma-separated taxonomic IDs
	 * @return Array of parsed taxonomic IDs, or null if invalid format
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
	 * Retries up to 12 times with exponential backoff starting at 200ms.
	 *
	 * @param path Server URL path
	 * @param prefix URL prefix for the query type
	 * @param message Query message to send
	 * @return Server response string, or null if all retries failed
	 */
	private static final String sendAndReceiveSpecificServer(String path, String prefix, String message){
		String response=null;
		for(int i=0, millis=200; i<12 && (response==null || response.length()==0); i++) {
			if(i>0){
				System.err.println("Null response, retrying after "+millis+" ms.");
				try {
					Thread.sleep(millis);
				} catch (InterruptedException e) {
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
	 * Uses GET for messages under 2000 characters, POST for longer messages.
	 *
	 * @param path Server URL path
	 * @param prefix URL prefix for the query type
	 * @param message Query message to send
	 * @return Server response string, or null if communication failed
	 */
	private static String sendAndReceiveOnceSpecificServer(String path, String prefix, String message){
		final String response;
		if(message.length()<2000){//NERSC Apache limit is around 8kb
			message=prefix+PercentEncoding.symbolToCode(message);
			ByteBuilder bb=ServerTools.readPage(path+message, false);
			response=(bb==null ? null : bb.toString());
			//			System.err.println("S&R1 send:\n"+path+message);
			//			System.err.println("S&R1 receive:\n"+response);
		}else{
			StringNum sn=ServerTools.sendAndReceive((prefix+message).getBytes(), path+"$POST");
			response=sn.s;
			//			System.err.println("S&R2: "+message+"\n"+path+"$POST");
		}
		return response;
	}
	
	/**
	 * Concatenates ArrayList elements into comma-separated string.
	 * @param list ArrayList of objects to concatenate
	 * @return Comma-separated string representation, or null if list is null/empty
	 */
	private static String fuse(ArrayList<?> list){
		if(list==null || list.size()<0){return null;} //Possible bug: size() is never negative
		StringBuilder sb=new StringBuilder();
		for(Object s : list){
			sb.append(s).append(',');
		}
		sb.setLength(sb.length()-1);
		return sb.toString();
	}
	
	/** Output stream for results, defaults to System.err */
	public static PrintStream outstream=System.err;
	/** Enable verbose output for debugging */
	public static boolean verbose=false;
	
}
