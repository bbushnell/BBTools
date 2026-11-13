package test.barcode;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.net.InetSocketAddress;
import java.net.URI;
import java.util.ArrayList;
import java.util.Locale;

import com.sun.net.httpserver.Headers;
import com.sun.net.httpserver.HttpContext;
import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpPrincipal;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import sketch.Sketch;

/** Public-only in-memory fixture for DemuxServer assignment-reply write failures. */
public final class DemuxServerReplyErrorTest{

	/** Requires mode=baseline|fixed and case=... . */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, DemuxServerReplyErrorTest.class, false);
		try{
			String mode=null, testCase=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(key.equals("case")){testCase=value;}
				else if(!Parser.parseStatic(arg, key, value)){
					throw new IllegalArgumentException("Unknown argument: "+arg);
				}
			}
			check("baseline".equals(mode) || "fixed".equals(mode),
					"Set mode=baseline or mode=fixed.");
			check("raw-success".equals(testCase) || "body-write-io".equals(testCase) ||
					"body-close-io".equals(testCase) || "headers-io".equals(testCase),
					"Set case=raw-success, body-write-io, body-close-io, or headers-io.");
			final CaseResult result=invoke(testCase);
			checkCase(mode, testCase, result);
			System.out.println("DEMUX_SERVER_REPLY_ERROR mode="+mode+" case="+testCase+
					" returned="+result.returned+" status="+result.status+
					" bodyCloseCalls="+result.bodyCloseCalls+
					" exchangeCloseCalls="+result.exchangeCloseCalls+
					" failed=0");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Invokes the private DemuxServer reply({@code ArrayList<Pair>}, ...) method. */
	private static CaseResult invoke(final String testCase) throws Exception{
		final RecordingExchange exchange=new RecordingExchange(
				"headers-io".equals(testCase),
				"body-write-io".equals(testCase),
				"body-close-io".equals(testCase));
		final ArrayList<Object> pairs=new ArrayList<Object>(1);
		pairs.add(pair("ACGT", "TTTT"));
		final Method reply=Class.forName("barcode.DemuxServer").getDeclaredMethod(
				"reply", ArrayList.class, int.class, int.class, int.class, String.class,
				HttpExchange.class, boolean.class, int.class, boolean.class);
		reply.setAccessible(true);
		final boolean returned=((Boolean)reply.invoke(null, pairs, Sketch.RAW, 4, 0,
				"text/plain", exchange, false, 200, true)).booleanValue();
		return new CaseResult(returned, exchange.status, exchange.headers.getFirst("Content-Type"),
				exchange.body.toString("UTF-8"), exchange.body.closeCalls,
				exchange.closeCalls);
	}

	/** Constructs DemuxServer.Pair without broadening its production visibility. */
	private static Object pair(final String key, final String value) throws Exception{
		final Constructor<?> constructor=Class.forName("barcode.DemuxServer$Pair").
				getDeclaredConstructor(String.class, String.class);
		constructor.setAccessible(true);
		return constructor.newInstance(key, value);
	}

	/** Checks exact status, body, return value, and owned-close observations. */
	private static void checkCase(final String mode, final String testCase,
			final CaseResult result){
		if("raw-success".equals(testCase)){
			check(result.returned, "RAW success should return true: "+result);
			check(result.status==200, "RAW success should send status 200: "+result);
			check("text/plain".equals(result.contentType),
					"RAW success should set text/plain: "+result);
			check("ACGT\tTTTT\n".equals(result.body),
					"RAW success body mismatch: "+result);
			check(result.bodyCloseCalls==1 && result.exchangeCloseCalls==1,
					"RAW success should close body and exchange once: "+result);
		}else if("body-write-io".equals(testCase)){
			check(result.status==200, "Body-write failure should send headers first: "+result);
			check(result.body.length()==0, "Body-write failure must not emit a partial body: "+result);
			check(result.bodyCloseCalls==1 && result.exchangeCloseCalls==1,
					"Body-write failure should close body and exchange once: "+result);
			check("baseline".equals(mode)==result.returned,
					"Body-write return mismatch for "+mode+": "+result);
		}else if("body-close-io".equals(testCase)){
			check(result.status==200, "Body-close failure should send headers first: "+result);
			check("ACGT\tTTTT\n".equals(result.body),
					"Body-close failure should write the complete RAW body: "+result);
			check(result.bodyCloseCalls==2 && result.exchangeCloseCalls==1,
					"Body-close failure should attempt body close twice and exchange close once: "+
					result);
			check("baseline".equals(mode)==result.returned,
					"Body-close return mismatch for "+mode+": "+result);
		}else{
			check(!result.returned, "Header failure should already return false: "+result);
			check(result.status==200, "Header failure should record attempted status 200: "+result);
			check(result.body.length()==0, "Header failure must not write a body: "+result);
			check(result.bodyCloseCalls==0 && result.exchangeCloseCalls==1,
					"Header failure should close only the exchange: "+result);
		}
	}

	/** Always enforce fixture oracles independent of assertion configuration. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private DemuxServerReplyErrorTest(){}

	/** Captured observations from one private reply invocation. */
	private static final class CaseResult{

		CaseResult(final boolean returned_, final int status_, final String contentType_,
				final String body_, final int bodyCloseCalls_, final int exchangeCloseCalls_){
			returned=returned_;
			status=status_;
			contentType=contentType_;
			body=body_;
			bodyCloseCalls=bodyCloseCalls_;
			exchangeCloseCalls=exchangeCloseCalls_;
		}

		@Override
		public String toString(){
			return "returned="+returned+", status="+status+", contentType="+
					contentType+", body='"+body.replace("\n", "\\n")+"', bodyCloseCalls="+
					bodyCloseCalls+", exchangeCloseCalls="+exchangeCloseCalls;
		}

		final boolean returned;
		final int status;
		final String contentType;
		final String body;
		final int bodyCloseCalls;
		final int exchangeCloseCalls;
	}

	/** HttpExchange with controllable header and body failure points. */
	private static final class RecordingExchange extends HttpExchange{

		RecordingExchange(final boolean throwHeaders_, final boolean throwWrite_,
				final boolean throwClose_){
			throwHeaders=throwHeaders_;
			body=new RecordingOutputStream(throwWrite_, throwClose_);
		}

		@Override
		public Headers getRequestHeaders(){return requestHeaders;}

		@Override
		public Headers getResponseHeaders(){return headers;}

		@Override
		public URI getRequestURI(){return URI.create("/demux");}

		@Override
		public String getRequestMethod(){return "POST";}

		@Override
		public HttpContext getHttpContext(){return null;}

		@Override
		public void close(){closeCalls++;}

		@Override
		public InputStream getRequestBody(){return null;}

		@Override
		public OutputStream getResponseBody(){return body;}

		@Override
		public void sendResponseHeaders(final int code, final long length) throws IOException{
			status=code;
			responseLength=length;
			if(throwHeaders){throw new IOException("header failure");}
		}

		@Override
		public InetSocketAddress getRemoteAddress(){return null;}

		@Override
		public int getResponseCode(){return status;}

		@Override
		public InetSocketAddress getLocalAddress(){return null;}

		@Override
		public String getProtocol(){return "HTTP/1.1";}

		@Override
		public Object getAttribute(final String name){return null;}

		@Override
		public void setAttribute(final String name, final Object value){}

		@Override
		public void setStreams(final InputStream input, final OutputStream output){}

		@Override
		public HttpPrincipal getPrincipal(){return null;}

		final Headers headers=new Headers();
		final Headers requestHeaders=new Headers();
		final RecordingOutputStream body;
		final boolean throwHeaders;
		int status=-1;
		long responseLength=-1;
		int closeCalls=0;
	}

	/** OutputStream that records bytes and optionally fails writes or closes. */
	private static final class RecordingOutputStream extends OutputStream{

		RecordingOutputStream(final boolean throwWrite_, final boolean throwClose_){
			throwWrite=throwWrite_;
			throwClose=throwClose_;
		}

		@Override
		public void write(final int value) throws IOException{
			if(throwWrite){throw new IOException("body write failure");}
			bytes.write(value);
		}

		@Override
		public void write(final byte[] buffer, final int offset, final int length) throws IOException{
			if(throwWrite){throw new IOException("body write failure");}
			bytes.write(buffer, offset, length);
		}

		@Override
		public void close() throws IOException{
			closeCalls++;
			if(throwClose){throw new IOException("body close failure");}
		}

		String toString(final String charsetName) throws IOException{
			return bytes.toString(charsetName);
		}

		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		final boolean throwWrite;
		final boolean throwClose;
		int closeCalls=0;
	}
}
