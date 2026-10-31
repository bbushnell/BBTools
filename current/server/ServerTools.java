package server;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.InetSocketAddress;
import java.net.MalformedURLException;
import java.net.URI;
import java.net.URL;
import java.net.URLConnection;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map.Entry;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;
import java.util.zip.GZIPOutputStream;

import com.sun.net.httpserver.Headers;
import com.sun.net.httpserver.HttpExchange;

import fileIO.ReadWrite;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import structures.StringNum;

/**
 * Utilities for small BBTools HTTP clients and HttpExchange response helpers.
 *
 * @author Brian Bushnell
 * @date June 3, 2025
 */
public class ServerTools{

	/** Constructs a ServerTools helper; current helpers are static. */
	public ServerTools(){}

	/**
	 * Diagnostic client that sends the same optional message several times.
	 * Prints each request duration in nanoseconds, the final response body,
	 * the mean duration, throughput for timed sends excluding reporting,
	 * and the upper median duration.
	 * @param args {@code <address> [rounds] [message]}
	 */
	public static void main(String[] args){
		if(args.length<1 || args.length>3){
			throw new IllegalArgumentException("Usage: ServerTools <address> [rounds] [message]");
		}
		int rounds=1;
		final String address=args[0];
		if(args.length>1){rounds=Integer.parseInt(args[1]);}
		if(rounds<1){throw new IllegalArgumentException("rounds must be positive: "+rounds);}
		final String message=args.length>2 ? args[2] : "";
		final byte[] messageBytes=message.getBytes();

		final long[] times=new long[rounds];
		StringNum response=null;
		for(int i=0; i<rounds; i++){
			final long start=System.nanoTime();
			response=sendAndReceive(messageBytes, address);
			times[i]=System.nanoTime()-start;
			System.out.println(times[i]);
		}

		System.out.println(response.s);

		Arrays.sort(times);
		long sum=0;
		for(long time : times){sum+=time;}
		System.out.println("Avg:    \t"+(sum/(rounds*1000000.0))+" ms");
		System.out.println("QPS:    \t"+(rounds*1000000000.0/sum)+" qps");
		System.out.println("Median: \t"+(times[rounds/2]/1000000.0)+" ms");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Client Helpers        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Reads one URL into a ByteBuilder, opening and closing each retry's
	 * stream. Failed attempts clear partial bytes before the next retry.
	 *
	 * @param address URL to read; optionally percent-encoded before opening
	 * @param convert true to convert common symbols to percent escapes first
	 * @return page bytes, or an empty builder if all ten attempts fail
	 */
	public static ByteBuilder readPage(String address, boolean convert){
		if(convert){address=PercentEncoding.commonSymbolToCode(address);}
		final ByteBuilder bb=new ByteBuilder(256);
		boolean success=false;
		for(int i=0; i<10 && !success; i++){
			//FIXED [server/ServerTools#005]: read failures now close the URL stream before retrying.
			try{
				try(InputStream is=new URL(address).openStream()){
					byte[] buffer=new byte[4096];
					for(int len=is.read(buffer); len>0; len=is.read(buffer)){
						bb.append(buffer, 0, len);
					}
				}
				success=true;
			}catch(MalformedURLException e){
				e.printStackTrace();
				bb.clear();
				Tools.pause(1000);
				System.err.println("Retrying; attempt "+(i+1)+", URL "+address);
			}catch(IOException e){
				e.printStackTrace();
				bb.clear();
				Tools.pause(1000);
				System.err.println("Retrying; attempt "+(i+1)+", URL "+address);
			}
		}
		return bb;
	}

	/**
	 * Sends one message with the Java 11 client when available,
	 * otherwise via the legacy Java 8 URLConnection path.
	 * @param message request body, or null for no body
	 * @param address target URL
	 * @return response body and status code */
	public static StringNum sendAndReceive(byte[] message, String address){
		if(Shared.javaVersion>=11){
			return sendAndReceive_httpClient(message, address);
		}
		return sendAndReceive_httpUrlConnection(message, address);
	}


	/**
	 * Legacy Java 8 HTTP POST path; a null message sends no body and
	 * status 1 means no HTTP response was read.
	 * @param message request body, or null for no body
	 * @param address target URL
	 * @return response body and status code */
	public static StringNum sendAndReceive_httpUrlConnection(byte[] message, String address){
		address=PercentEncoding.commonSymbolToCode(address);
		URL url=null;
		InputStream is=null;
		HttpURLConnection connection=null;
		OutputStream os=null;
		try{
			url=new URL(address);
			connection=(HttpURLConnection) url.openConnection();
			connection.setDoOutput(true);
			connection.setConnectTimeout(40000); //For testing
			os=connection.getOutputStream();
		}catch(IOException e1){
			if(!suppressErrors){e1.printStackTrace();}
		}

		if(os!=null){
			try{
				//TODO: It may be useful to set a timeout here.
				if(message!=null){os.write(message);}
			}catch(IOException e){
				if(!suppressErrors){e.printStackTrace();}
			}
			try{
				os.close();
			}catch(IOException e){
				e.printStackTrace();
			}
		}

		String result=null;
		int responseCode=1;
		if(connection!=null){
			IOException noInputStream=null;
			try{
				is=connection.getInputStream();
			}catch(IOException e){
				noInputStream=e;
			}
			if(is==null){is=connection.getErrorStream();}

			try{
				responseCode=connection.getResponseCode();
				if(!suppressErrors && (responseCode<200 || responseCode>299)){
					System.err.println("Error: Server returned response code "+responseCode);
					new Exception().printStackTrace();
				}
			}catch(IOException e){
				if(!suppressErrors){e.printStackTrace();}
			}

			if(is!=null){
				result=readStream(is);
				try{
					is.close();
				}catch(IOException e){
					e.printStackTrace();
				}
			}else if(noInputStream!=null && !suppressErrors){
				noInputStream.printStackTrace();
			}
		}
		if(result!=null && result.length()<1){result=null;}
		return new StringNum(result, responseCode);
	}

	/**
	 * Modern Java 11 HTTP POST path; retries I/O and non-2xx
	 * responses but preserves interrupt cancellation.
	 * @param message request body, or null for no body
	 * @param address target URL
	 * @return response body and status code */
	public static StringNum sendAndReceive_httpClient(byte[] message, String address){
		address=PercentEncoding.commonSymbolToCode(address);

		//FIXED [server/ServerTools#006]: null means no request body, as documented and as
		//the Java 8 path already allowed.
		final HttpRequest.BodyPublisher body=message==null ?
				HttpRequest.BodyPublishers.noBody() : HttpRequest.BodyPublishers.ofByteArray(message);
		final HttpRequest request=HttpRequest.newBuilder()
				  .uri(URI.create(address))
				  .POST(body)
				  .build();
		HttpResponse<String> response=null;

		int responseCode=-1;
		String result=null;

		for(int i=0; i<14 && (result==null || responseCode<200 || responseCode>299); i++){
			try{
				if(i>0){Tools.sleep(40*i*i);}
				response=client.send(request, HttpResponse.BodyHandlers.ofString());
				result=response.body();
				responseCode=response.statusCode();
			}catch(InterruptedException e2){
				//FIXED [server/ServerTools#007]: preserve cancellation and do not retry an
				//interrupted send as ordinary I/O.
				Thread.currentThread().interrupt();
				e2.printStackTrace();
				break;
			}catch(IOException e2){
				e2.printStackTrace();
			}
			if(!suppressErrors && (responseCode<200 || responseCode>299) && i>0){
				System.err.println("Warning: Server returned response code "+responseCode);
				new Exception().printStackTrace();
			}
		}
		if(result!=null && result.length()<1){result=null;}
		return new StringNum(result, responseCode);
	}

	/** Compresses all input bytes and returns a complete GZIP member.
	 * @param data uncompressed bytes
	 * @return GZIP-compressed bytes, or the bytes written before an unexpected compression error */
	public static byte[] gzipCompress(byte[] data){
		//FIXED [server/ServerTools#001]: GZIPOutputStream is now closed before
		//toByteArray(), flushing the deflater's buffered body and gzip trailer.
		final ByteArrayOutputStream bos=new ByteArrayOutputStream();
		try(GZIPOutputStream gzip=new GZIPOutputStream(bos)){
			gzip.write(data);
		}catch(IOException e){
			e.printStackTrace();
		}
		return bos.toByteArray();
	}

	/**
	 * Returns the first Content-Encoding header value from a Java 11
	 * response, or empty string if absent.
	 * @param httpResponse response to inspect
	 * @return first Content-Encoding header, or an empty string if absent */
	public static String determineContentEncoding(HttpResponse<?> httpResponse){
		return httpResponse.headers().firstValue("Content-Encoding").orElse("");
	}


	/**
	 * Legacy URLConnection POST path for chunked request bodies; consumes
	 * and closes the server response.
	 * @param messages request body chunks; null chunks are skipped
	 * @param address target URL
	 * @param verbose true to print send and receive timings
	 * @return response body and status code */
	public static StringNum sendAndReceive(ArrayList<byte[]> messages, String address,
			boolean verbose){
		address=PercentEncoding.commonSymbolToCode(address);
		URL url=null;
		InputStream is=null;
		HttpURLConnection connection=null;
		OutputStream os=null;
		try{
			url=new URL(address);
			connection=(HttpURLConnection) url.openConnection();
			connection.setDoOutput(true);
			connection.setConnectTimeout(1600000); //For testing
			os=connection.getOutputStream();
		}catch(IOException e1){
			if(!suppressErrors){e1.printStackTrace();}
		}

		long sent=0, received=0;
		final Timer t=new Timer();
		if(os!=null){
			try{
				//TODO: It may be useful to set a timeout here.
				for(byte[] message : messages){
					if(message!=null){
						os.write(message);
						sent+=message.length;
					}
				}
			}catch(IOException e){
				if(!suppressErrors){e.printStackTrace();}
			}
			try{
				os.close();
			}catch(IOException e){
				e.printStackTrace();
			}
		}
		if(verbose){t.stopAndStart("Sent "+sent+" bytes in ");}

		String result=null;
		int responseCode=1;
		if(connection!=null){
			IOException noInputStream=null;
			try{
				is=connection.getInputStream();
			}catch(IOException e){
				noInputStream=e;
			}
			if(is==null){is=connection.getErrorStream();}

			try{
				responseCode=connection.getResponseCode();
				if(!suppressErrors && (responseCode<200 || responseCode>299)){
					System.err.println("Error: Server returned response code "+responseCode);
					new Exception().printStackTrace();
				}
			}catch(IOException e){
				if(!suppressErrors){e.printStackTrace();}
			}

			if(is!=null){
				result=readStream(is);
				received+=(result==null ? 0 : result.length());
				try{
					is.close();
				}catch(IOException e){
					e.printStackTrace();
				}
			}else if(noInputStream!=null && !suppressErrors){
				noInputStream.printStackTrace();
			}
		}
		if(verbose){t.stop("Recieved "+received+" bytes in ");}

		return new StringNum(result, responseCode);
	}


	/**
	 * Legacy FTP URLConnection send/read helper retained for historical callers;
	 * it always returns status 400.
	 * @param message request body, or null for no body
	 * @param address target URL
	 * @return response body with a historical status code of 400 */
	public static StringNum sendAndReceiveFTP(byte[] message, String address){
		address=PercentEncoding.commonSymbolToCode(address);
		URL url=null;
		InputStream is=null;
		URLConnection connection=null;
		OutputStream os=null;
		try{
			url=new URL(address);
			connection=url.openConnection();
			connection.setDoOutput(true);
			connection.setConnectTimeout(40000); //For testing
			os=connection.getOutputStream();
		}catch(IOException e1){
			if(!suppressErrors){e1.printStackTrace();}
		}

		if(os!=null){
			try{
				//TODO: It may be useful to set a timeout here.
				if(message!=null){os.write(message);}
			}catch(IOException e){
				if(!suppressErrors){e.printStackTrace();}
			}
			try{
				os.close();
			}catch(IOException e){
				e.printStackTrace();
			}
		}

		String result=null;
		if(connection!=null){
			IOException noInputStream=null;
			try{
				is=connection.getInputStream();
			}catch(IOException e){
				noInputStream=e;
			}

			if(is!=null){
				result=readStream(is);
				try{
					is.close();
				}catch(IOException e){
					e.printStackTrace();
				}
			}else if(noInputStream!=null && !suppressErrors){
				noInputStream.printStackTrace();
			}
		}

		return new StringNum(result, 400);
	}

	/*--------------------------------------------------------------*/
	/*----------------  Stream and Response Helpers ----------------*/
	/*--------------------------------------------------------------*/

	/** Consumes and closes the body of an incoming HTTP session.
	 * @param t HTTP exchange to read from
	 * @return request body text, or null on read or close failure */
	public static String receive(HttpExchange t){
		return readStream(t.getRequestBody());
	}

	/**
	 * Consumes and closes an InputStream using the default charset;
	 * returns null on read or close failure.
	 * @param is stream to consume
	 * @return stream contents, or null on read or close failure */
	public static String readStream(InputStream is){
		if(is==null){return null;}
		//FIXED [server/ServerTools#008]: close even when a read throws.
		try(InputStream in=is){
			byte[] buffer=new byte[256];
			int count=in.read(buffer);
			int next=0;
			while(count>-1){
				next+=count;
				if(next>=buffer.length){
					buffer=Arrays.copyOf(buffer, buffer.length*2);
				}
				count=in.read(buffer, next, buffer.length-next);
			}
			return new String(buffer, 0, next);
		}catch(IOException e){
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * Consumes and closes an InputStream, returning byte chunks of at most
	 * about 8 MB or null on read or close failure.
	 * @param is stream to consume
	 * @return non-empty byte chunks, or null on read or close failure */
	public static ArrayList<byte[]> readStreamToList(InputStream is){
		if(is==null){return null;}
		ArrayList<byte[]> list=new ArrayList<byte[]>();
		ByteBuilder bb=new ByteBuilder();
		long next=0;
		final byte[] buffer=new byte[4096];
		//FIXED [server/ServerTools#009]: close even when a chunk read throws.
		try(InputStream in=is){
			int count=in.read(buffer);
			while(count>-1){
				next+=count;
				bb.append(buffer, count);
				if(bb.length()>8000000){
					list.add(bb.toBytes());
					bb.clear();
				}
				count=in.read(buffer);
			}
			if(bb.length()>0){
				list.add(bb.toBytes());
				bb.clear();
			}
			System.err.println("Read "+next+" bytes into "+list.size()+" chunks.");
			return list;
		}catch(IOException e){
			e.printStackTrace();
		}
		return null;
	}

	/** Writes a String response using the platform-default charset.
	 * @param response response text
	 * @param type Content-Type header value
	 * @param t HTTP exchange to write
	 * @param verbose true to print the response text
	 * @param code HTTP status code
	 * @param close true to close the exchange after writing or on failure
	 * @return true if the response bytes were written successfully */
	public static boolean reply(String response, String type, HttpExchange t, boolean verbose,
			int code, boolean close){
		if(verbose){System.err.println("Sending: "+response);}

		return reply(response.getBytes(), type, t, false, code, close);
	}

	/** Writes a file response, substituting an empty 400 response if the file read fails.
	 * @param path file to read into the response body
	 * @param type Content-Type header value
	 * @param t HTTP exchange to write
	 * @param verbose true to print the source path
	 * @param code HTTP status code for a successful file read
	 * @param close true to close the exchange after writing or on failure
	 * @return true if the response bytes were written successfully */
	public static boolean replyWithFile(String path, String type, HttpExchange t, boolean verbose,
			int code, boolean close){
		if(verbose){System.err.println("Sending: "+path);}

		byte[] response=null;
		try{
			response=ReadWrite.readRaw(path);
		}catch(IOException e1){
			e1.printStackTrace();
			response=new byte[0];
			code=400; //Kept as the historical read-failure status.
		}

		return reply(response, type, t, false, code, close);
	}

	/**
	 * Writes exactly response.length bytes, closes the body on success
	 * and closes the exchange if requested.
	 * @param response response bytes
	 * @param type Content-Type header value
	 * @param t HTTP exchange to write
	 * @param verbose true to print the response array identity
	 * @param code HTTP status code
	 * @param close true to close the exchange after writing or on failure
	 * @return true if the response bytes were written successfully */
	public static boolean reply(byte[] response, String type, HttpExchange t, boolean verbose,
			int code, boolean close){
		if(verbose){System.err.println("Sending: "+response);}

		{
			final Headers h=t.getResponseHeaders();
			h.add("Content-Type", type);
		}
		try{
			t.sendResponseHeaders(code, response.length);
			final OutputStream os=t.getResponseBody();
			os.write(response);
			os.close();
		}catch(IOException e){
			e.printStackTrace();
			if(close){t.close();}
			return false;
		}
		if(close){t.close();}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Address Helpers        ----------------*/
	/*--------------------------------------------------------------*/


	/**
	 * Waits until the requested interval elapses, printing interruptions
	 * but preserving the full wait.
	 * @param millis milliseconds to wait */
	public static void pause(long millis){
		String lock=new String("1");
		synchronized(lock){
			final long time=System.currentTimeMillis()+millis;
			while(System.currentTimeMillis()<time){
				try{
					lock.wait(Tools.max(100, time-System.currentTimeMillis()));
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
		}
	}

	/** Returns X-Forwarded-For for loopback proxies, otherwise the remote socket address string.
	 * @param t HTTP exchange to inspect
	 * @return forwarded client or remote socket address string */
	public static String getClientAddress(HttpExchange t){

		final InetSocketAddress client=t.getRemoteAddress();

		//This is for IPv4, class A.  Probably extends outside of Berkeley.
		String clientAddress=client.toString();

		if(clientAddress.contains("127.0.0.1")){
			final Headers clientRequestHeaders=t.getRequestHeaders();

			final String xff=clientRequestHeaders.getFirst("X-forwarded-for");
			if(xff!=null){clientAddress=xff;}
		}
		return clientAddress;
	}

	/**
	 * Applies the current internal-address policy. Loopback requests may be
	 * allowed directly, or may inherit X-Forwarded-For. Non-loopback requests
	 * are accepted on a literal prefix match. Finally, the client and server
	 * socket strings may match through a '.' or ':' delimiter.
	 *
	 * @param t HTTP exchange to examine
	 * @param prefix IP prefix for internal networks
	 * @param allowLocalHost whether to accept loopback with no X-Forwarded-For
	 * @param printIP whether to print client and server addresses
	 * @param printHeaders whether to print all request/response headers
	 * @return true if request is from internal network, false otherwise
	 */
	public static boolean isInternalQuery(HttpExchange t,
			String prefix, boolean allowLocalHost, boolean printIP, boolean printHeaders){

		final InetSocketAddress client=t.getRemoteAddress();
		final InetSocketAddress server=t.getLocalAddress();

		if(printIP){System.err.println(client+"\t"+server);}

		//This is for IPv4, class A.  Probably extends outside of Berkeley.
		String clientAddress=client.toString();
		String serverAddress=server.toString();

		if(clientAddress.contains("127.0.0.1")){//TODO: contains versus startsWith?
			final Headers requestHeaders=t.getRequestHeaders();

			if(printHeaders){
				final Headers responseHeaders=t.getResponseHeaders();
				System.err.println("\nRequest: ");
				for(Entry<String, List<String>> entry : requestHeaders.entrySet()){
					System.err.println(entry.getKey()+" -> "+entry.getValue());
				}
				System.err.println("\nResponse: ");
				for(Entry<String, List<String>> entry : responseHeaders.entrySet()){
					System.err.println(entry.getKey()+" -> "+entry.getValue());
				}
			}

			final String xff=requestHeaders.getFirst("X-forwarded-for");
			if(xff!=null){
				if(xff.startsWith(prefix)){return true;}
				clientAddress=xff;
			}else{
				return allowLocalHost;
			}
		}else{
			if(clientAddress.startsWith(prefix)){return true;}
		}

		//Makes sure they match up to the first delimiter
		//FIXED [server/ServerTools#004]: was Tools.max; now compare only
		//the overlapping region of the client and server address strings.
		for(int i=0, min=Tools.min(clientAddress.length(), serverAddress.length()); i<min; i++){
			final char cc=clientAddress.charAt(i), sc=serverAddress.charAt(i);
			if(cc!=sc){break;}
			if(cc=='.'){//IPv4
				return true;
			}else if(cc==':'){//IPv6; probably depends on how long the mask is
				return true;
			}
		}

		return false;
	}

	/*--------------------------------------------------------------*/
	/*----------------         FTP Helpers          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Lists an FTP directory via ";type=d", retrying and returning only a fully
	 * read listing. Each raw line is prefixed with the parent URL of baseAddress;
	 * FetchProks consumes those full paths recursively. This assumes JDK FTP
	 * directory listings return names relative to that parent path.
	 * @param baseAddress FTP directory URL
	 * @param retries number of retries after the first attempt
	 * @return full child paths from the listing, or an empty list if all attempts fail
	 */
	public static ArrayList<String> listDirectory(String baseAddress, final int retries){
		while(baseAddress.endsWith("/")){baseAddress=baseAddress.substring(0, baseAddress.length()-1);}
		final String baseAddress2=baseAddress.substring(0, baseAddress.lastIndexOf('/')+1);
		final String address=baseAddress+";type=d";
		final ArrayList<String> list=new ArrayList<String>();
		boolean success=false;
		for(int i=Tools.min(0, retries); !success && i<=retries; i++){
			try{
				final URL url=new URL(address);
				final URLConnection connection=url.openConnection();
				//FIXED [server/ServerTools#010]: close the FTP directory reader on read failure.
				try(BufferedReader br=new BufferedReader(new InputStreamReader(connection.getInputStream()))){
					for(String line=br.readLine(); line!=null; line=br.readLine()){
						list.add(baseAddress2+line);
					}
				}

				success=true;
			}catch(IOException e){
				e.printStackTrace();
				//FIXED [server/ServerTools#011]: a final failed attempt must not return partial entries.
				list.clear();
				if(i<retries){
					try{
						final int x=Tools.mid(30000, i*1000, 1000);
						System.err.println("Sleeping for "+x+"ms...");
						Thread.sleep(x);
					}catch(InterruptedException e1){
						e1.printStackTrace();
					}
					System.err.println("** Retrying ** "+address);
				}else{
					System.err.println("*** Gave up ***");
				}
			}
		}
		return list;
	}

	/** Reads a percent-encoded FTP URL into one line per list element.
	 * @param address FTP URL to read
	 * @return file lines
	 * @throws Exception if URL construction, opening, reading, or closing fails */
	public static ArrayList<String> readFTPFile(String address) throws Exception{
		//FIXED [server/ServerTools#002]: reader now closed via try-with-resources.
		address=PercentEncoding.commonSymbolToCode(address);
		final URL url=new URL(address);
		final URLConnection conn=url.openConnection();

		final ArrayList<String> list=new ArrayList<String>();
		try(BufferedReader reader=new BufferedReader(new InputStreamReader(conn.getInputStream()))){
			for(String s=reader.readLine(); s!=null; s=reader.readLine()){
				list.add(s);
			}
		}
		return list;
	}

	/*--------------------------------------------------------------*/
	/*----------------    Server-Side Utilities     ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses FASTA text into Reads, skipping sequence before the first
	 * header and using default-charset bases.
	 * @param fastaBody FASTA text
	 * @return parsed reads with generated numeric ids */
	public static ArrayList<Read> parseFastaBody(String fastaBody){
		final ArrayList<Read> reads=new ArrayList<Read>();
		final String[] lines=fastaBody.split("\n");
		final StringBuilder seqBuilder=new StringBuilder();
		String currentName=null;
		int id=0;
		for(String line : lines){
			line=line.trim();
			if(line.startsWith(">")){
				if(currentName!=null && seqBuilder.length()>0){
					reads.add(new Read(seqBuilder.toString().getBytes(), null, currentName, id++));
				}
				currentName=line.substring(1).trim();
				seqBuilder.setLength(0);
			}else if(!line.isEmpty()){
				seqBuilder.append(line);
			}
		}
		if(currentName!=null && seqBuilder.length()>0){
			reads.add(new Read(seqBuilder.toString().getBytes(), null, currentName, id++));
		}
		return reads;
	}

	/** Appends a JSON-string-safe ASCII encoding of a UTF-16 Java String.
	 * @param bb destination builder
	 * @param s string to escape */
	//FIXED [server/ServerTools#003]: now escapes every control char below 0x20.
	//FIXED [server/ServerTools#012]: non-ASCII UTF-16 code units are now written
	//as ASCII backslash-u escapes rather than being truncated by ByteBuilder.append(char).
	public static void escapeJson(ByteBuilder bb, String s){
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c=='"'){bb.append('\\').append('"');}
			else if(c=='\\'){bb.append('\\').append('\\');}
			else if(c=='\n'){bb.append('\\').append('n');}
			else if(c=='\t'){bb.append('\\').append('t');}
			else if(c=='\r'){bb.append('\\').append('r');}
			else if(c=='\b'){bb.append('\\').append('b');}
			else if(c=='\f'){bb.append('\\').append('f');}
			else if(c<0x20){//Remaining C0 control chars -> \\u00XX
				appendJsonUnicode(bb, c);
			}
			else if(c>0x7F){appendJsonUnicode(bb, c);}
			else{bb.append(c);}
		}
	}

	/** Append one UTF-16 code unit as an ASCII JSON Unicode escape.
	 * @param bb destination builder
	 * @param c code unit to escape */
	private static void appendJsonUnicode(ByteBuilder bb, char c){
		bb.append('\\').append('u');
		for(int shift=12; shift>=0; shift-=4){
			final int x=(c>>shift)&0xF;
			bb.append((char)(x<10 ? '0'+x : 'A'+x-10));
		}
	}

	/** Adds permissive CORS headers, using '*' when domain is null.
	 * @param t HTTP exchange to mutate
	 * @param domain Access-Control-Allow-Origin value, or null for any origin */
	public static void addCorsHeaders(HttpExchange t, String domain){
		Headers h=t.getResponseHeaders();
		h.add("Access-Control-Allow-Origin", domain!=null ? domain : "*");
		h.add("Access-Control-Allow-Methods", "GET, POST, OPTIONS");
		h.add("Access-Control-Allow-Headers", "Content-Type");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Suppresses some legacy caught-exception diagnostics. */
	public static boolean suppressErrors=false;

	/** Shared Java 11+ HTTP client whose executor uses daemon threads so diagnostics can exit. */
	private static HttpClient client;
	static{
		try{

			//Without daemon threads the process won't terminate after the shared client starts.
			final ThreadFactory daemonThreadFactory=new ThreadFactory(){
				public Thread newThread(Runnable r){
					final Thread t=new Thread(r);
					t.setDaemon(true);
					return t;
				}
			};

			final ExecutorService executor=Executors.newFixedThreadPool(4, daemonThreadFactory);
			client=HttpClient.newBuilder()
				.version(HttpClient.Version.HTTP_2)
				.executor(executor)
				.build();
		}catch(Throwable e){
			System.err.println("You may need a more recent version of Java, specifically, "+
					"11 or higher.  Your version: "+Shared.javaVersion);
			e.printStackTrace();
		}
	}

}
