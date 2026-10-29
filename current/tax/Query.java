package tax;

import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Arrays;

import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * Legacy blocking URL client for GI/accession-to-taxID lookups and throughput probes.
 * Unlike TaxClient, this class does not encode request payloads, retry transport
 * failures, or validate response shape before parsing. Server configuration is
 * global; individual calls allocate their own request and response storage.
 * @author Brian Bushnell
 */
public class Query{

	/** Constructs the stateless legacy client; all operations are static. */
	public Query(){}

	/** Times repeated reads of one complete URL, discarding the returned content.
	 * Reports overall throughput and a second interval excluding the first request.
	 * Transport failures are printed but do not stop the timing loop.
	 * @param args Complete URL followed by an optional positive request count */
	public static void main(String[] args){
		String x=args[0];
		int requests=(args.length>1 ? Integer.parseInt(args[1]) : 1);
		Timer t=new Timer();
		Timer t2=new Timer();//Retained legacy timer; not used by the reports.
		Timer t3=new Timer();
		for(int i=0; i<requests; i++){
			String s=request(x);
			if(i==0){t3.start();}
		}
		t3.stop();
		t.stop("Time:  \t");
		double qps=requests*1000000000.0/(t.elapsed);
		System.err.println(Tools.format("Qps:   \t%.2f", qps));
		System.err.println("Time2: \t"+t3);
		double qps2=(requests-1)*1000000000.0/(t3.elapsed);
		System.err.println(Tools.format("Qps2:  \t%.2f", qps2));
	}

	/*--------------------------------------------------------------*/
	/*----------------         TaxID Queries        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolves GI numbers through the configured taxonomy server.
	 * @param args GI numbers, represented as strings
	 * @return Parsed server taxonomy IDs in response order */
	public static int[] getGi(String... args){return get("pt_gi", args);}

	/** Resolves accession identifiers through the configured taxonomy server.
	 * @param args Accession strings
	 * @return Parsed server taxonomy IDs in response order */
	public static int[] getAccession(String... args){return get("pt_accession", args);}

	/** Builds a legacy query and parses its comma-separated integer response.
	 * Payloads are appended literally; callers must supply URL-compatible values.
	 * No response-count check or null transport-failure guard is performed here.
	 * @param type Legacy route such as pt_gi or pt_accession
	 * @param args Identifier strings to join with commas
	 * @return Newly allocated array of server taxonomy IDs
	 * @throws NumberFormatException If any response field is not a valid int
	 * @throws NullPointerException If the transport returns null */
	public static int[] get(String type, String... args){
		String host=Shared.taxServer();
		StringBuilder sb=new StringBuilder(host.length()+32);
		sb.append(host);
		sb.append(type);
		sb.append('/');
		String comma="";
		for(String arg : args){
			sb.append(comma);
			sb.append(arg);
			comma=",";
		}
		String result=request(sb.toString());
		String[] split=result.split(",");
		int[] ret=new int[split.length];
		for(int i=0; i<split.length; i++){
			ret[i]=Integer.parseInt(split[i]);
		}
		return ret;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Transport          ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads an entire URL response, closing the stream on normal and handled I/O exits.
	 * The byte buffer grows geometrically; the result uses the platform default
	 * charset. On an IOException, prints the failure, attempts closure and returns
	 * null. Other unchecked failures can bypass closure. No timeout or size limit
	 * is configured by this method.
	 * @param x Complete URL
	 * @return Response text, an empty string for EOF, or null after a handled failure */
	public static String request(String x){
		InputStream is=toStream(x);
		if(is==null){return null;}
		try{
			byte[] buffer=new byte[256];
			int count=is.read(buffer);
			int next=0;
			while(count>-1){
				next+=count;
				if(next>=buffer.length){
					buffer=Arrays.copyOf(buffer, buffer.length*2);
				}
				count=is.read(buffer, next, buffer.length-next);
			}
			is.close();
			return new String(buffer, 0, next);
		}catch(IOException e){
			e.printStackTrace();
		}
		try{
			is.close();
		}catch(IOException e){
			e.printStackTrace();
		}
		return null;
	}

	/** Opens a URL stream, printing handled URL/connection failures.
	 * The caller owns and must close a returned stream.
	 * @param x Complete URL
	 * @return Open stream, or null on MalformedURLException or IOException */
	public static InputStream toStream(String x){
		URL url=null;
		try{
			url=new URL(x);
		}catch(MalformedURLException e1){
			e1.printStackTrace();
			return null;
		}
		try{
			InputStream is=url.openStream();
			return is;
		}catch(IOException e){
			e.printStackTrace();
		}
		return null;
	}
}
