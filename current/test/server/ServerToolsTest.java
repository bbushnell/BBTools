package test.server;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.net.Authenticator;
import java.net.CookieHandler;
import java.net.InetSocketAddress;
import java.net.ProxySelector;
import java.net.URI;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLStreamHandler;
import java.net.URLStreamHandlerFactory;
import java.net.http.HttpClient;
import java.net.http.HttpHeaders;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.nio.charset.StandardCharsets;
import java.time.Duration;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;

import javax.net.ssl.SSLContext;
import javax.net.ssl.SSLParameters;
import javax.net.ssl.SSLSession;

import com.sun.net.httpserver.Headers;
import com.sun.net.httpserver.HttpContext;
import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpPrincipal;

import server.ServerTools;
import structures.ByteBuilder;
import structures.StringNum;

/** Finite ServerTools checks using fake streams, exchanges and clients; no listener by default. */
public final class ServerToolsTest {

	/** Runs against either the saved baseline or repaired ServerTools. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1 || !(args[0].equals("baseline") || args[0].equals("fixed"))){
			throw new IllegalArgumentException("Expected baseline or fixed");
		}
		final boolean fixed=args[0].equals("fixed");
		checkReadStreamClose("readstream-close-on-failure", fixed);
		checkReadStreamToListClose("readstreamtolist-close-on-failure", fixed);
		checkReadPageRetryClose("readpage-retry-close-on-failure", fixed);
		checkReadPageCloseFailureRetry("readpage-close-failure-retry-control");
		checkListDirectoryFinalFailure("listdirectory-final-failure-empty-close", fixed);
		checkNullHttpClient("java11-null-body", fixed);
		checkInterruptedHttpClient("httpclient-interrupted-retry-restores-status", fixed);
		checkEscapeJsonUnicode("escapejson-unicode-ascii", fixed);
		checkMainOutput("main-mean-qps-validation", fixed);
		checkReceive("receive-delegates-body");
		checkReplyFailure("reply-close-on-write-failure");
		checkAddressControls();
	}

	/** A read failure must not bypass InputStream.close(). */
	private static void checkReadStreamClose(String name, boolean expected){
		final ThrowingInputStream stream=new ThrowingInputStream("abc".getBytes());
		final String result=ServerTools.readStream(stream);
		expect(name, result==null && stream.closed, expected);
	}

	/** The chunked request-body reader has the same ownership requirement as readStream. */
	private static void checkReadStreamToListClose(String name, boolean expected){
		final ThrowingInputStream stream=new ThrowingInputStream("abcdef".getBytes());
		final Object result=ServerTools.readStreamToList(stream);
		expect(name, result==null && stream.closed, expected);
	}

	/** readPage should close a failed URL stream before clearing and retrying. */
	private static void checkReadPageRetryClose(String name, boolean expected) throws Exception{
		final ThrowingInputStream first=new ThrowingInputStream("bad".getBytes(StandardCharsets.UTF_8));
		final Scenario scenario=new Scenario(first, new ByteArrayInputStream("ok".getBytes(StandardCharsets.UTF_8)));
		final ByteBuilder result=ServerTools.readPage(scenarioUrl(scenario), false);
		final boolean good="ok".equals(result.toString()) && scenario.opens==2 && first.closed;
		expect(name, good, expected);
	}

	/** A close-time failure after a complete read must still be retried. */
	private static void checkReadPageCloseFailureRetry(String name) throws Exception{
		final CloseThrowingInputStream first=
				new CloseThrowingInputStream("bad".getBytes(StandardCharsets.UTF_8));
		final Scenario scenario=new Scenario(first, new ByteArrayInputStream("ok".getBytes(StandardCharsets.UTF_8)));
		final ByteBuilder result=ServerTools.readPage(scenarioUrl(scenario), false);
		final boolean good="ok".equals(result.toString()) && scenario.opens==2 && first.closed;
		expect(name, good, true);
	}

	/** listDirectory must close and discard partial entries after its final failed read. */
	private static void checkListDirectoryFinalFailure(String name, boolean expected) throws Exception{
		final ThrowingInputStream stream=new ThrowingInputStream("first\nsecond".getBytes(StandardCharsets.UTF_8));
		final Scenario scenario=new Scenario(stream);
		final ArrayList<String> result=ServerTools.listDirectory(scenarioUrl(scenario), 0);
		final boolean good=result.isEmpty() && scenario.opens==1 && stream.closed;
		expect(name, good, expected);
	}

	/** sendAndReceive documents null as a no-body request and the Java-8 path honored that. */
	private static void checkNullHttpClient(String name, boolean expected) throws Exception{
		final ClientState state=new ClientState(new FakeHttpClient(false, "ok"));
		try{
			boolean good=false;
			try{
				final StringNum result=ServerTools.sendAndReceive_httpClient(null, "http://127.0.0.1/null");
				good=result!=null && result.n==200 && "ok".equals(result.s) &&
						state.client.sends==1 && state.client.lastBodyLength==0;
			}catch(NullPointerException e){good=false;}
			expect(name, good, expected);
		}finally{state.restore();}
	}

	/** Cancellation should be preserved instead of being consumed and retried as a normal failure. */
	private static void checkInterruptedHttpClient(String name, boolean expected) throws Exception{
		final ClientState state=new ClientState(new FakeHttpClient(true, "ok"));
		try{
			Thread.interrupted();
			final StringNum result=ServerTools.sendAndReceive_httpClient(new byte[] {1}, "http://127.0.0.1/interrupted");
			final boolean interrupted=Thread.currentThread().isInterrupted();
			final boolean good=result!=null && result.n<0 && state.client.sends==1 && interrupted;
			Thread.interrupted();
			expect(name, good, expected);
		}finally{
			Thread.interrupted();
			state.restore();
		}
	}

	/** JSON escapes must stay ASCII when writing through ByteBuilder bytes. */
	private static void checkEscapeJsonUnicode(String name, boolean expected){
		final char escape='\\';
		final String input="\"\\\n\t\r\b\f"+((char)1)+((char)0x00E9)+((char)0x6F22)+
				new String(Character.toChars(0x1F600));
		final String oracle="\\\"\\\\\\n\\t\\r\\b\\f"+
				escape+"u0001"+escape+"u00E9"+escape+"u6F22"+escape+"uD83D"+escape+"uDE00";
		final ByteBuilder bb=new ByteBuilder();
		ServerTools.escapeJson(bb, input);
		final byte[] bytes=bb.toBytes();
		boolean ascii=true;
		for(byte b : bytes){
			if((b&0x80)!=0){ascii=false;}
		}
		expect(name, ascii && oracle.equals(bb.toString()), expected);
	}

	/** Diagnostic main should print a true mean, label QPS as throughput and reject bad args. */
	private static void checkMainOutput(String name, boolean expected) throws Exception{
		final ClientState state=new ClientState(new FakeHttpClient(false, "ok"));
		final PrintStream old=System.out;
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		final PrintStream capture=new PrintStream(bytes, true, "UTF-8");
		boolean good=false;
		try{
			System.setOut(capture);
			boolean badCount=false, badRounds=false, goodOutput=false;
			try{
				ServerTools.main(new String[] {});
			}catch(IllegalArgumentException e){badCount=e.getMessage().contains("Usage:");}catch(RuntimeException e){badCount=false;}
			try{
				ServerTools.main(new String[] {"http://127.0.0.1/main", "0"});
			}catch(IllegalArgumentException e){badRounds=e.getMessage().contains("rounds");}catch(RuntimeException e){badRounds=false;}
			bytes.reset();
			ServerTools.main(new String[] {"http://127.0.0.1/main", "3", "msg"});
			final String[] lines=bytes.toString("UTF-8").replace("\r", "").split("\n");
			if(lines.length==7){
				final long a=Long.parseLong(lines[0]), b=Long.parseLong(lines[1]);
				final long c=Long.parseLong(lines[2]), sum=a+b+c;
				final String expectedMean="Avg:    \t"+(sum/3000000.0)+" ms";
				goodOutput="ok".equals(lines[3]) && expectedMean.equals(lines[4]) &&
						lines[5].startsWith("QPS:    \t") && lines[5].endsWith(" qps") &&
						state.client.sends==3 && state.client.lastBodyLength==3;
			}
			good=badCount && badRounds && goodOutput;
		}finally{
			System.setOut(old);
			capture.close();
			state.restore();
		}
		expect(name, good, expected);
	}

	private static void checkReceive(String name){
		final Exchange exchange=new Exchange("body".getBytes(), false, "127.0.0.1", "10.1.2.3");
		expect(name, "body".equals(ServerTools.receive(exchange)), true);
	}

	private static void checkReplyFailure(String name){
		final Exchange exchange=new Exchange(new byte[0], true, "127.0.0.1", "10.1.2.3");
		final boolean ok=ServerTools.reply("abc".getBytes(), "text/plain", exchange, false, 200, true);
		expect(name, !ok && exchange.closed, true);
	}

	/** Address cases are controls for the current prefix/header policy, not a policy rewrite. */
	private static void checkAddressControls(){
		final Exchange loopback=new Exchange(new byte[0], false, "127.0.0.1", "10.1.2.3");
		loopback.request.add("X-forwarded-for", "10.4.5.6");
		expect("getclientaddress-xff", ServerTools.getClientAddress(loopback).equals("10.4.5.6"), true);
		expect("internal-xff-prefix", ServerTools.isInternalQuery(loopback, "10.", false, false, false), true);

		final Exchange allowedLocal=new Exchange(new byte[0], false, "127.0.0.1", "10.1.2.3");
		expect("internal-loopback-allowed", ServerTools.isInternalQuery(allowedLocal, "10.", true, false, false), true);
		expect("internal-loopback-denied", ServerTools.isInternalQuery(allowedLocal, "10.", false, false, false), false);

		final Exchange external=new Exchange(new byte[0], false, "8.8.8.8", "10.1.2.3");
		expect("external-prefix-miss", ServerTools.isInternalQuery(external, "10.", false, false, false), false);
	}

	private static void expect(String name, boolean actual, boolean expected){
		System.out.println(name+"\t"+actual+"\t"+expected);
		if(actual!=expected){throw new AssertionError(name+": expected "+expected+", observed "+actual);}
	}

	private static String scenarioUrl(Scenario scenario){
		installUrlHandler();
		synchronized(SCENARIOS){
			final String host="case"+(nextScenario++);
			SCENARIOS.put(host, scenario);
			return TEST_PROTOCOL+"://"+host+"/entry";
		}
	}

	private static synchronized void installUrlHandler(){
		if(handlerInstalled){return;}
		URL.setURLStreamHandlerFactory(new URLStreamHandlerFactory() {
			@Override
			public URLStreamHandler createURLStreamHandler(String protocol){
				if(!TEST_PROTOCOL.equals(protocol)){return null;}
				return new URLStreamHandler() {
					@Override
					protected URLConnection openConnection(URL url) throws IOException{
						return new ScenarioConnection(url);
					}
				};
			}
		});
		handlerInstalled=true;
	}

	private static synchronized Scenario scenario(URL url) throws IOException{
		final Scenario scenario=SCENARIOS.get(url.getHost());
		if(scenario==null){throw new IOException("No scenario for "+url);}
		return scenario;
	}

	/** Saves and replaces the private static HttpClient, then restores it. */
	private static final class ClientState {
		ClientState(FakeHttpClient client_) throws Exception{
			client=client_;
			field=ServerTools.class.getDeclaredField("client");
			field.setAccessible(true);
			original=(HttpClient)field.get(null);
			field.set(null, client);
		}
		void restore() throws Exception{field.set(null, original);}
		final Field field;
		final HttpClient original;
		final FakeHttpClient client;
	}

	/** Returns one successful response unless configured to interrupt once first. */
	private static final class FakeHttpClient extends HttpClient {
		FakeHttpClient(boolean interrupt_, String body_){
			interrupt=interrupt_;
			body=body_;
		}

		@Override
		public <T> HttpResponse<T> send(HttpRequest request, HttpResponse.BodyHandler<T> handler)
				throws IOException, InterruptedException{
			sends++;
			lastBodyLength=request.bodyPublisher().isPresent() ?
					request.bodyPublisher().get().contentLength() : -1;
			if(interrupt && sends==1){throw new InterruptedException("injected");}
			@SuppressWarnings("unchecked")
			final T cast=(T)body;
			return new Response<T>(request, cast);
		}

		@Override
		public <T> CompletableFuture<HttpResponse<T>> sendAsync(
				HttpRequest request, HttpResponse.BodyHandler<T> handler){
			throw new UnsupportedOperationException();
		}

		@Override
		public <T> CompletableFuture<HttpResponse<T>> sendAsync(HttpRequest request,
				HttpResponse.BodyHandler<T> handler, HttpResponse.PushPromiseHandler<T> pushPromiseHandler){
			throw new UnsupportedOperationException();
		}

		@Override public Optional<CookieHandler> cookieHandler(){return Optional.empty();}
		@Override public Optional<Duration> connectTimeout(){return Optional.empty();}
		@Override public Redirect followRedirects(){return Redirect.NEVER;}
		@Override public Optional<ProxySelector> proxy(){return Optional.empty();}
		@Override public SSLContext sslContext(){return null;}
		@Override public SSLParameters sslParameters(){return null;}
		@Override public Optional<Authenticator> authenticator(){return Optional.empty();}
		@Override public Version version(){return Version.HTTP_1_1;}
		@Override public Optional<Executor> executor(){return Optional.empty();}

		final boolean interrupt;
		final String body;
		int sends;
		long lastBodyLength;
	}

	private static final class Response<T> implements HttpResponse<T> {
		Response(HttpRequest request_, T body_){
			request=request_;
			body=body_;
		}

		@Override public int statusCode(){return 200;}
		@Override public HttpRequest request(){return request;}
		@Override public Optional<HttpResponse<T>> previousResponse(){return Optional.empty();}
		@Override public HttpHeaders headers(){return HttpHeaders.of(java.util.Collections.<String, java.util.List<String>>emptyMap(), (a, b)->true);}
		@Override public T body(){return body;}
		@Override public Optional<SSLSession> sslSession(){return Optional.empty();}
		@Override public URI uri(){return request.uri();}
		@Override public HttpClient.Version version(){return HttpClient.Version.HTTP_1_1;}

		final HttpRequest request;
		final T body;
	}

	private static final class ScenarioConnection extends URLConnection {
		ScenarioConnection(URL url){super(url);}

		@Override
		public void connect(){}

		@Override
		public InputStream getInputStream() throws IOException{
			return scenario(url).open();
		}
	}

	private static final class Scenario {
		Scenario(InputStream... streams_){streams=streams_;}

		synchronized InputStream open() throws IOException{
			if(opens>=streams.length){throw new IOException("No stream "+opens);}
			return streams[opens++];
		}

		final InputStream[] streams;
		int opens;
	}

	/** Gives callers data once, then throws while still expecting close(). */
	private static final class ThrowingInputStream extends InputStream {
		ThrowingInputStream(byte[] data_){data=data_;}

		@Override
		public int read() throws IOException{
			if(pos<data.length){return data[pos++];}
			throw new IOException("injected");
		}

		@Override
		public int read(byte[] buffer, int off, int len) throws IOException{
			if(pos>0){throw new IOException("injected");}
			final int count=Math.min(len, data.length);
			System.arraycopy(data, 0, buffer, off, count);
			pos+=count;
			return count;
		}

		@Override
		public void close(){closed=true;}

		final byte[] data;
		int pos;
		boolean closed;
	}

	/** Gives callers all data and then throws only on close(). */
	private static final class CloseThrowingInputStream extends ByteArrayInputStream {
		CloseThrowingInputStream(byte[] data){super(data);}

		@Override
		public void close() throws IOException{
			closed=true;
			throw new IOException("injected close");
		}

		boolean closed;
	}

	/** In-memory exchange with optional output failure and controllable addresses. */
	private static final class Exchange extends HttpExchange {
		Exchange(byte[] body, final boolean fail, String remote, String local){
			input=new ByteArrayInputStream(body);
			remoteAddress=new InetSocketAddress(remote, 1);
			localAddress=new InetSocketAddress(local, 2);
			output=new OutputStream() {
				@Override
				public void write(int b) throws IOException{
					if(fail){throw new IOException("injected");}
					bytes.write(b);
				}
				@Override
				public void close(){outputClosed=true;}
			};
		}

		@Override public Headers getRequestHeaders(){return request;}
		@Override public Headers getResponseHeaders(){return response;}
		@Override public URI getRequestURI(){return URI.create("/");}
		@Override public String getRequestMethod(){return "POST";}
		@Override public HttpContext getHttpContext(){return null;}
		@Override public void close(){closed=true; outputClosed=true;}
		@Override public InputStream getRequestBody(){return input;}
		@Override public OutputStream getResponseBody(){return output;}
		@Override public void sendResponseHeaders(int code, long length){status=code; declared=length;}
		@Override public InetSocketAddress getRemoteAddress(){return remoteAddress;}
		@Override public InetSocketAddress getLocalAddress(){return localAddress;}
		@Override public String getProtocol(){return "HTTP/1.1";}
		@Override public Object getAttribute(String name){return null;}
		@Override public void setAttribute(String name, Object value){}
		@Override public void setStreams(InputStream in, OutputStream out){throw new UnsupportedOperationException();}
		@Override public HttpPrincipal getPrincipal(){return null;}
		@Override public int getResponseCode(){return status;}

		final Headers request=new Headers();
		final Headers response=new Headers();
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		final InputStream input;
		final OutputStream output;
		final InetSocketAddress remoteAddress;
		final InetSocketAddress localAddress;
		int status;
		long declared;
		boolean closed, outputClosed;
	}

	private static final String TEST_PROTOCOL="servertest";
	private static final Map<String, Scenario> SCENARIOS=new HashMap<String, Scenario>();
	private static boolean handlerInstalled=false;
	private static int nextScenario=0;
}
