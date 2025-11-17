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
 * Provides HTTP-based query interface for taxonomic database lookups and identifier resolution.
 * Handles communication with taxonomic servers to retrieve GI numbers and accession mappings.
 * Supports performance testing and throughput measurement for database queries.
 * @author Brian Bushnell
 */
public class Query {
	
	/** For testing throughput */
	public static void main(String[] args){
		
//		sun.net.http.errorstream.enableBuffering=true;

//		System.setProperty("sun.net.http.errorstream.enableBuffering", "true");
//		System.setProperty("http.errorstream.enableBuffering", "true");
//////		System.setProperty("", "");
//		System.setProperty("http.keepAlive", "true");
//		System.setProperty("http.maxConnections", "50");
////		System.setProperty("", "");
////		System.setProperty("", "");
//
//		System.out.println("sun.net.http.errorstream.enableBuffering="+System.getProperty("sun.net.http.errorstream.enableBuffering"));
//		System.out.println("sun.net.http.errorstream.timeout="+System.getProperty("sun.net.http.errorstream.timeout"));
//		System.out.println("sun.net.http.errorstream.bufferSize="+System.getProperty("sun.net.http.errorstream.bufferSize"));
//		System.out.println("http.keepAlive="+System.getProperty("http.keepAlive"));
//		System.out.println("http.maxConnections="+System.getProperty("http.maxConnections"));
//		System.out.println(System.getProperties());
//		System.out.println("="+System.getProperty(""));
//		System.out.println("="+System.getProperty(""));
//		System.out.println("="+System.getProperty(""));
		
		String x=args[0];
		int requests=(args.length>1 ? Integer.parseInt(args[1]) : 1);

		Timer t=new Timer();
		Timer t2=new Timer();
		Timer t3=new Timer();
		
//		URL url=toUrl(x);
		for(int i=0; i<requests; i++){
			String s=request(x);
//			System.out.println(s);
			if(i==0){t3.start();}
			//t2.stop("");
			//t2.start();
		}
		
		t3.stop();
		
		t.stop("Time:  \t");
		double qps=requests*1000000000.0/(t.elapsed);
		System.err.println(Tools.format("Qps:   \t%.2f", qps));
		
		System.err.println("Time2: \t"+t3);
		double qps2=(requests-1)*1000000000.0/(t3.elapsed);
		System.err.println(Tools.format("Qps2:  \t%.2f", qps2));
		
//		System.out.println(s);
		
//		BufferedReader reader;
//		try {
//			reader = new BufferedReader(new InputStreamReader(url.openStream(), "UTF-8"));
//
//			for (String line; (line = reader.readLine()) != null;) {
//				System.out.println(line);
//			}
//
//		} catch (UnsupportedEncodingException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
	}
	
	/**
	 * Retrieves GI numbers for given taxonomic identifiers.
	 * @param args Variable number of taxonomic identifiers to query
	 * @return Array of GI numbers corresponding to the input identifiers
	 */
	public static int[] getGi(String...args){
		return get("pt_gi", args);
	}
	
	/**
	 * Retrieves accession numbers for given taxonomic identifiers.
	 * @param args Variable number of taxonomic identifiers to query
	 * @return Array of accession numbers corresponding to the input identifiers
	 */
	public static int[] getAccession(String...args){
		return get("pt_accession", args);
	}
	
	/**
	 * Generic method for retrieving taxonomic data of specified type.
	 * Constructs query URL and parses comma-separated response into integer array.
	 *
	 * @param type Query type identifier (e.g., "pt_gi", "pt_accession")
	 * @param args Variable number of taxonomic identifiers to query
	 * @return Array of integers parsed from server response
	 */
	public static int[] get(String type, String...args){
		String host=Shared.taxServer();
		StringBuilder sb=new StringBuilder(host.length()+32);
		sb.append(host);
		sb.append(type);
		sb.append('\t');
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
	
	/**
	 * Executes HTTP request and returns response as string.
	 * Reads input stream into dynamically expanding buffer and handles IOException.
	 * @param x URL string for the HTTP request
	 * @return Response content as string, or null if request fails
	 */
	public static String request(String x){
		InputStream is=toStream(x);
		if(is==null){return null;}
		try {
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
			is.close(); //TODO: Why am I not closing this stream?
			return new String(buffer, 0, next);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
			is.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	
	/**
	 * Converts URL string to InputStream for reading HTTP response.
	 * Handles MalformedURLException and IOException with error printing.
	 * @param x URL string to convert
	 * @return InputStream for reading response, or null if URL is invalid or connection fails
	 */
	public static InputStream toStream(String x){
		URL url=null;
		
		try {
			url = new URL(x);
		} catch (MalformedURLException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
			return null;
		}
		
		try {
			InputStream is=url.openStream();
			return is;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	
//	public static URL toUrl(String x){
//		URL url=null;
//		try {
//			url = new URL(x);
//		} catch (MalformedURLException e1) {
//			// TODO Auto-generated catch block
//			e1.printStackTrace();
//		}
//
//		return url;
//	}
//
//	public static URLConnection toUrlConnection(URL url, String x){
//		URLConnection urlc=null;
//
//		try {
//			urlc = url.openConnection();
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		return urlc;
//	}
	
}
