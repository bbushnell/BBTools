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
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;
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
import structures.ByteBuilder;
import structures.StringNum;

public class ServerTools {
	
	public static void main(String[] args){
		String address=args[0];
		int rounds=1;
		String message="";
		if(args.length>1){rounds=Integer.parseInt(args[1]);}
		if(args.length>2){message=args[2];}
		byte[] messageBytes=message.getBytes();
		
		long[] times=new long[rounds];
		StringNum response=null;
		long prevTime=System.nanoTime();
		for(int i=0; i<rounds; i++){
			response=sendAndReceive(messageBytes, address);
			long currentTime=System.nanoTime();
			times[i]=currentTime-prevTime;
			prevTime=currentTime;
			System.out.println(times[i]);
		}
		
		System.out.println(response.s);
		
		Arrays.sort(times);
		long sum=shared.Vector.sum(times);
		System.out.println("Avg:    \t"+sum/1000000.0+" ms");
		System.out.println("QPS:    \t"+(rounds*1000000000/sum)+" ms");
		System.out.println("Median: \t"+(times[rounds/2]/1000000.0)+" ms");
		
	}
	
	public static ByteBuilder readPage(String address, boolean convert){
    	if(convert){address=PercentEncoding.commonSymbolToCode(address);}
//    	assert(false) : address;
		ByteBuilder bb=new ByteBuilder(256);
		boolean success=false;
		for(int i=0; i<10 && !success; i++){
			try {
				URL url=new URL(address);
				InputStream is=url.openStream();

				byte[] buffer=new byte[4096];
				for(int len=is.read(buffer); len>0; len=is.read(buffer)){
					bb.append(buffer, 0, len);
				}
				is.close();
				success=true;
			} catch (MalformedURLException e) {
				e.printStackTrace();
				bb.clear();
				Tools.pause(1000);
				System.err.println("Retrying; attempt "+(i+1)+", URL "+address);
			} catch (IOException e) {
				e.printStackTrace();
				bb.clear();
				Tools.pause(1000);
				System.err.println("Retrying; attempt "+(i+1)+", URL "+address);
			}
		}
        return bb;
    }
	
	/** Send a message to a remote URL, and return the response.
	 * Set message to null if there is no message. */
	public static StringNum sendAndReceive(byte[] message, String address){
		//This may or may not work since the client is statically initialized.
		if(Shared.javaVersion>=11) {
			return sendAndReceive_httpClient(message, address);
		}
		return sendAndReceive_httpUrlConnection(message, address);
	}
	
	
	/** Send a message to a remote URL, and return the response.
	 * Set message to null if there is no message. */
	public static StringNum sendAndReceive_httpUrlConnection(byte[] message, String address){
    	address=PercentEncoding.commonSymbolToCode(address);
		URL url=null;
		InputStream is=null;
		HttpURLConnection connection=null;
		OutputStream os=null;
		try {
			url=new URL(address);
			connection=(HttpURLConnection) url.openConnection();
			connection.setDoOutput(true);
			connection.setConnectTimeout(40000); //For testing
			os=connection.getOutputStream();
		} catch (IOException e1) {
//			System.err.println("A:\t"+address+" -> "+url+" -> "+connection+" -> "+os);
			// TODO Auto-generated catch block
			if(!suppressErrors){e1.printStackTrace();}
		}
		
		if(os!=null){
			try {
				//TODO: It may be useful to set a timeout here.
				if(message!=null){os.write(message);}
			} catch (IOException e) {
//				System.err.println("B:\t"+connection+" -> "+os);
				// TODO Auto-generated catch block
				if(!suppressErrors){e.printStackTrace();}
			}
			try {
				os.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		String result=null;
		int responseCode=1;
		if(connection!=null){
			IOException noInputStream=null;
			try {
				is=connection.getInputStream();
			} catch (IOException e) {
				noInputStream=e;
			}
			if(is==null){is=connection.getErrorStream();}
			
			try {
				responseCode=connection.getResponseCode();
				if(!suppressErrors && (responseCode<200 || responseCode>299)){
					System.err.println("Error: Server returned response code "+responseCode);
//					System.err.println(connection.getRequestProperties());//Can't do this - throws exception while connected
//					System.err.println(connection.getResponseMessage());
					new Exception().printStackTrace();
				}
			} catch (IOException e) {
				if(!suppressErrors) {
					e.printStackTrace();
				}
			}
			
			if(is!=null){
				result=readStream(is);
				try {
//					System.err.println("C:\t"+connection+" -> "+os+" -> "+is+" -> "+(result!=null));
					is.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}else if(noInputStream!=null && !suppressErrors){
				noInputStream.printStackTrace();
			}
		}
		if(result!=null && result.length()<1) {result=null;}
		return new StringNum(result, responseCode);
	}
	
	/** Send a message to a remote URL, and return the response.
	 * Set message to null if there is no message. */
	public static StringNum sendAndReceive_httpClient(byte[] message, String address){
    	address=PercentEncoding.commonSymbolToCode(address);
    	
		HttpRequest request = HttpRequest.newBuilder()
				  .uri(URI.create(address))
//				  .header("Accept-Encoding", "gzip")
				  .POST(HttpRequest.BodyPublishers.ofByteArray(message))
				  .build();
		HttpResponse<String> response=null;
		
		int responseCode=-1;
		String result=null;
		
		for(int i=0; i<14 && (result==null || responseCode<200 || responseCode>299); i++) {
			try {
				if(i>0) {Tools.sleep(40*i*i);}
				response = client.send(request, HttpResponse.BodyHandlers.ofString());
				result=response.body();
				responseCode=response.statusCode();
			} catch (IOException | InterruptedException e2) {
				e2.printStackTrace();
//				return null;
			}
			if(!suppressErrors && (responseCode<200 || responseCode>299) && i>0){
				System.err.println("Warning: Server returned response code "+responseCode);
				new Exception().printStackTrace();
			}
		}
		if(result!=null && result.length()<1) {result=null;}
		return new StringNum(result, responseCode);
	}
	
	public static byte[] gzipCompress(byte[] data) {
		ByteArrayOutputStream bos = new ByteArrayOutputStream();
		try {
			 GZIPOutputStream gzip = new GZIPOutputStream(bos);
			gzip.write(data);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        return bos.toByteArray();
    }
	
	public static String determineContentEncoding(
	        HttpResponse<?> httpResponse) {
	    return httpResponse.headers().firstValue("Content-Encoding").orElse("");
	}
	
	
	/** Send a message to a remote URL, and return the response.
	 * Set message to null if there is no message. */
	public static StringNum sendAndReceive(ArrayList<byte[]> messages, String address, 
			boolean verbose){
    	address=PercentEncoding.commonSymbolToCode(address);
		URL url=null;
		InputStream is=null;
		HttpURLConnection connection=null;
		OutputStream os=null;
		try {
			url=new URL(address);
			connection=(HttpURLConnection) url.openConnection();
			connection.setDoOutput(true);
			connection.setConnectTimeout(1600000); //For testing
			os=connection.getOutputStream();
		} catch (IOException e1) {
//			System.err.println("A:\t"+address+" -> "+url+" -> "+connection+" -> "+os);
			// TODO Auto-generated catch block
			if(!suppressErrors){e1.printStackTrace();}
		}
		
		long sent=0, received=0;
		Timer t=new Timer();
		if(os!=null){
			try {
				//TODO: It may be useful to set a timeout here.
				for(byte[] message : messages) {
					if(message!=null){
						os.write(message);
						sent+=message.length;
					}
				}
			} catch (IOException e) {
//				System.err.println("B:\t"+connection+" -> "+os);
				// TODO Auto-generated catch block
				if(!suppressErrors){e.printStackTrace();}
			}
			try {
				os.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		if(verbose) {t.stopAndStart("Sent "+sent+" bytes in ");}

		String result=null;
		int responseCode=1;
		if(connection!=null){
			IOException noInputStream=null;
			try {
				is=connection.getInputStream();
			} catch (IOException e) {
				noInputStream=e;
			}
			if(is==null){is=connection.getErrorStream();}
			
			try {
				responseCode=connection.getResponseCode();
				if(!suppressErrors && (responseCode<200 || responseCode>299)){
					System.err.println("Error: Server returned response code "+responseCode);
//					System.err.println(connection.getRequestProperties());//Can't do this - throws exception while connected
//					System.err.println(connection.getResponseMessage());
					new Exception().printStackTrace();
				}
			} catch (IOException e) {
				if(!suppressErrors) {
					e.printStackTrace();
				}
			}
			
			if(is!=null){
				result=readStream(is);
				received+=(result==null ? 0 : result.length());
				try {
//					System.err.println("C:\t"+connection+" -> "+os+" -> "+is+" -> "+(result!=null));
					is.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}else if(noInputStream!=null && !suppressErrors){
				noInputStream.printStackTrace();
			}
		}
		if(verbose) {t.stop("Recieved "+received+" bytes in ");}
		
		return new StringNum(result, responseCode);
	}
	
	
	/** Send a message to a remote URL, and return the response.
	 * Set message to null if there is no message. */
	public static StringNum sendAndReceiveFTP(byte[] message, String address){
    	address=PercentEncoding.commonSymbolToCode(address);
		URL url=null;
		InputStream is=null;
		URLConnection connection=null;
		OutputStream os=null;
		try {
			url=new URL(address);
			connection=url.openConnection();
			connection.setDoOutput(true);
			connection.setConnectTimeout(40000); //For testing
			os=connection.getOutputStream();
		} catch (IOException e1) {
//			System.err.println("A:\t"+address+" -> "+url+" -> "+connection+" -> "+os);
			// TODO Auto-generated catch block
			if(!suppressErrors){e1.printStackTrace();}
		}
		
		if(os!=null){
			try {
				//TODO: It may be useful to set a timeout here.
				if(message!=null){os.write(message);}
			} catch (IOException e) {
//				System.err.println("B:\t"+connection+" -> "+os);
				// TODO Auto-generated catch block
				if(!suppressErrors){e.printStackTrace();}
			}
			try {
				os.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		String result=null;
//		int responseCode=1;
		if(connection!=null){
			IOException noInputStream=null;
			try {
				is=connection.getInputStream();
			} catch (IOException e) {
				noInputStream=e;
			}
			
//			try {
				//responseCode=connection.getResponseCode();
//				if(!suppressErrors && (responseCode<200 || responseCode>299)){
//					System.err.println("Error: Server returned response code "+responseCode);
//				}
//			} 
//			catch (IOException e) {
//				if(!suppressErrors) {
//					e.printStackTrace();
//				}
//			}
			
			if(is!=null){
				result=readStream(is);
				try {
//					System.err.println("C:\t"+connection+" -> "+os+" -> "+is+" -> "+(result!=null));
					is.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}else if(noInputStream!=null && !suppressErrors){
				noInputStream.printStackTrace();
			}
		}
		
		return new StringNum(result, 400);
	}

	/** Read the body of an incoming HTTP session */
	public static String receive(HttpExchange t){
		InputStream is=t.getRequestBody();
		String s=readStream(is);
		return s;
	}
	
	/** Completely read an InputStream into a String */
	public static String readStream(InputStream is){
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
			is.close();
			return new String(buffer, 0, next);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	
	/** Completely read an InputStream into a String */
	public static ArrayList<byte[]> readStreamToList(InputStream is){
		if(is==null){return null;}
		ArrayList<byte[]> list=new ArrayList<byte[]>();
		ByteBuilder bb=new ByteBuilder();
		long next=0;
		final byte[] buffer=new byte[4096];
		try {
			int count=is.read(buffer);
			while(count>-1){
				next+=count;
				bb.append(buffer, count);
				if(bb.length()>8000000) {
					list.add(bb.toBytes());
					bb.clear();
				}
				count=is.read(buffer);
			}
			if(bb.length()>0) {
				list.add(bb.toBytes());
				bb.clear();
			}
			System.err.println("Read "+next+" bytes into "+list.size()+" chunks.");
			is.close();
			return list;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	
	/** Write to the body of an incoming HTTP session */
	public static boolean reply(String response, String type, HttpExchange t, boolean verbose, int code, boolean close){
		if(verbose){System.err.println("Sending: "+response);}

		return reply(response.getBytes(), type, t, false, code, close);
	}
	
	/** Write to the body of an incoming HTTP session */
	public static boolean replyWithFile(String path, String type, HttpExchange t, boolean verbose, int code, boolean close){
		if(verbose){System.err.println("Sending: "+path);}
		
		byte[] response=null;
		try {
			response=ReadWrite.readRaw(path);
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
			response=new byte[0];
			code=400; //Node sure about this
		}
		
		return reply(response, type, t, false, code, close);
	}
	
	/** Write to the body of an incoming HTTP session */
	public static boolean reply(byte[] response, String type, HttpExchange t, boolean verbose, int code, boolean close){
		if(verbose){System.err.println("Sending: "+response);}
		
		{
			Headers h = t.getResponseHeaders();
//			String type="text/plain";
			h.add("Content-Type", type);
		}
		try {
			t.sendResponseHeaders(code, response.length);
			OutputStream os = t.getResponseBody();
			os.write(response);
			os.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			if(close){t.close();}
			return false;
		}
		if(close){t.close();}
		return true;
	}
	
	
	/**
	 * Wait for a set amount of time
	 * @param millis Time to wait
	 */
	public static void pause(long millis){
		String lock=new String("1");
		synchronized(lock){
			final long time=System.currentTimeMillis()+millis;
			while(System.currentTimeMillis()<time){
				try {
					lock.wait(Tools.max(100, time-System.currentTimeMillis()));
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	
	public static String getClientAddress(HttpExchange t) {
		
		InetSocketAddress client=t.getRemoteAddress();
//		InetSocketAddress server=t.getLocalAddress();
		
		//This is for IPv4, class A.  Probably extends outside of Berkeley.
		String clientAddress=client.toString();
//		String ls=server.toString();
		
		if(clientAddress.contains("127.0.0.1")){
			Headers clientRequestHeaders=t.getRequestHeaders();
//			Headers resh=t.getResponseHeaders();
	
			String xff=clientRequestHeaders.getFirst("X-forwarded-for");
			if(xff!=null){clientAddress=xff;}
			
//			System.err.println("\nRequest: ");
//			for(Entry<String, List<String>> entry : reqh.entrySet()){
//				System.err.println(entry.getKey()+" -> "+entry.getValue());
//			}
		}
		return clientAddress;
	}
	
	public static boolean isInternalQuery(HttpExchange t, 
			String prefix, boolean allowLocalHost, boolean printIP, boolean printHeaders){
		
		InetSocketAddress client=t.getRemoteAddress();
		InetSocketAddress server=t.getLocalAddress();
		
		if(printIP){System.err.println(client+"\t"+server);}
		
		//This is for IPv4, class A.  Probably extends outside of Berkeley.
		String clientAddress=client.toString();
		String serverAddress=server.toString();
		
		if(clientAddress.contains("127.0.0.1")){//TODO: contains versus startsWith?
			Headers requestHeaders=t.getRequestHeaders();
			
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
		//TODO: This needs to be reviewed
		for(int i=0, max=Tools.max(clientAddress.length(), serverAddress.length()); i<max; i++){
			char cc=clientAddress.charAt(i), sc=serverAddress.charAt(i);
			if(cc!=sc){break;}
			if(cc=='.'){//IPv4
				return true; 
			}else if(cc==':'){//IPv6; probably depends on how long the mask is
				return true;
			}
		}
	
		return false;
	}
	
	public static ArrayList<String> listDirectory(String baseAddress, final int retries){
//		System.err.println("listDirectory '"+baseAddress+"'");
		while(baseAddress.endsWith("/")){baseAddress=baseAddress.substring(0, baseAddress.length()-1);}
		String baseAddress2=baseAddress.substring(0, baseAddress.lastIndexOf('/')+1);
		String address=baseAddress+";type=d";
		ArrayList<String> list=new ArrayList<String>();
		boolean success=false;
		for(int i=Tools.min(0, retries); !success && i<=retries; i++) {
			try {
				//			System.err.println("URL '"+address+"'");
				URL url=new URL(address);
				URLConnection connection=url.openConnection();
				InputStream is=connection.getInputStream();
				BufferedReader br=new BufferedReader(new InputStreamReader(is));

				String line=null;
				while((line=br.readLine())!=null) {
					list.add(baseAddress2+line);
				}

				is.close();
				success=true;
			} catch (IOException e) {
				e.printStackTrace();
				if(i<retries){
					try {
						int x=Tools.mid(30000, i*1000, 1000);
						System.err.println("Sleeping for "+x+"ms...");
						Thread.sleep(x);
					} catch (InterruptedException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					}
					System.err.println("** Retrying ** "+address);
					list.clear();
				}else{
					System.err.println("*** Gave up ***");
				}
			}
		}
		return list;
	}
	
	public static ArrayList<String> readFTPFile(String address) throws Exception{
    	address=PercentEncoding.commonSymbolToCode(address);
		URL url = new URL(address);
		URLConnection conn = url.openConnection();

		BufferedReader reader = new BufferedReader(new InputStreamReader(
		                conn.getInputStream()));
		ArrayList<String> list=new ArrayList<String>();
		for(String s=reader.readLine(); s!=null; s=reader.readLine()){
			list.add(s);
		}
		return list;
	}
	
	/** Don't print caught exceptions */
	public static boolean suppressErrors=false;
	
	private static HttpClient client;
	static{
		try {

			//Without make the threads Daemons the process won't terminate
			ThreadFactory daemonThreadFactory = new ThreadFactory() {
				public Thread newThread(Runnable r) {
					Thread t = new Thread(r);
					t.setDaemon(true);
					return t;
				}
			};

			ExecutorService executor = Executors.newFixedThreadPool(4, daemonThreadFactory);
//			ExecutorService executor = Executors.newFixedThreadPool(8);
	        client = HttpClient.newBuilder()
	                .version(HttpClient.Version.HTTP_2)
	                .executor(executor)
	                .build();
			
//			client=HttpClient.newHttpClient();
		} catch (Throwable e) {
			System.err.println("You may need a more recent version of Java, specifically, 11 or higher.  Your version: "+Shared.javaVersion);
			e.printStackTrace();
		}
	}
	
}
