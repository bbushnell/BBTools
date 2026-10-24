package test.server;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.lang.reflect.Constructor;
import java.net.InetSocketAddress;
import java.net.URI;
import java.nio.charset.StandardCharsets;
import java.util.LinkedHashMap;
import java.util.Map;

import com.sun.net.httpserver.Headers;
import com.sun.net.httpserver.HttpContext;
import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpPrincipal;

/** Finite handler checks using fake exchanges; never opens a listening socket. */
public class SimpleHttpServerTest {

	/** Runs against either the saved baseline or repaired handler. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1 || !(args[0].equals("baseline") || args[0].equals("fixed"))){
			throw new IllegalArgumentException("Expected baseline or fixed");
		}
		final boolean fixed=args[0].equals("fixed");
		final Constructor<?> ctor=Class.forName("server.SimpleHttpServer$GetHandler").getDeclaredConstructor();
		ctor.setAccessible(true);
		final HttpHandler handler=(HttpHandler)ctor.newInstance();
		checkResponse(handler, "ordinary", "/gi/123456", "gi", "123456", fixed);
		checkResponse(handler, "trailing-slash", "/gi/123456/", "gi", "123456", fixed);
		checkResponse(handler, "unicode", "/gi/\u00e9\u6f22", "gi", "\u00e9\u6f22", fixed);
		checkResponse(handler, "raw-percent", "/gi/a%2Fb", "gi", "a%2Fb", fixed);
		checkResponse(handler, "query-preserved", "/gi/a?x=y", "gi", "a?x=y", fixed);
		checkResponse(handler, "reserved-key", "/tax_id/1", "tax_id", "654321", fixed);
		checkResponse(handler, "bad-path-control", "/onlykey", "error",
				"need restful-style param like gi/123456", true);
		final Exchange failed=new Exchange("/gi/123456", true);
		boolean threw=false;
		try{handler.handle(failed);}catch(IOException expected){threw=true;}
		if(!threw){throw new AssertionError("Injected output failure was swallowed");}
		expect("failure-closes-exchange", failed.closed, fixed);
	}

	/** Independent flat-JSON oracle also verifies exact HTTP byte length. */
	private static void checkResponse(HttpHandler handler, String name, String uri,
			String key, String value, boolean expected) throws IOException{
		final Exchange exchange=new Exchange(uri, false);
		handler.handle(exchange);
		boolean good=exchange.code==200 && exchange.declared==exchange.bytes.size() && exchange.outputClosed;
		try{
			final Map<String,String> fields=parse(new String(exchange.bytes.toByteArray(), StandardCharsets.UTF_8));
			good&=value.equals(fields.get(key));
			if(!key.equals("error")){
				good&="654321".equals(fields.get("tax_id")) && "e. coli".equals(fields.get("organism"));
				good&=fields.size()==(key.equals("tax_id") || key.equals("organism") ? 2 : 3);
			}else{good&=fields.size()==1;}
		}catch(IllegalArgumentException e){good=false;}
		expect(name, good, expected);
	}

	private static void expect(String name, boolean actual, boolean expected){
		System.out.println(name+"\t"+actual+"\t"+expected);
		if(actual!=expected){throw new AssertionError(name+": expected valid="+expected+", observed="+actual);}
	}

	/** Parses only the fixture's JSON object of strings, rejecting missing separators. */
	private static Map<String,String> parse(String text){
		final Cursor c=new Cursor(text);
		final Map<String,String> map=new LinkedHashMap<String,String>();
		c.take('{');
		if(c.peek()!='}'){
			while(true){
				String key=c.string(); c.take(':'); String value=c.string();
				if(map.put(key, value)!=null){throw new IllegalArgumentException("Duplicate key");}
				if(c.peek()!=','){break;}
				c.take(',');
			}
		}
		c.take('}');
		if(c.peek()!=-1){throw new IllegalArgumentException("Trailing input");}
		return map;
	}

	/** Small independent syntax reader; no production JSON implementation is reused. */
	private static class Cursor {
		Cursor(String text_){text=text_;}
		int peek(){
			while(pos<text.length() && Character.isWhitespace(text.charAt(pos))){pos++;}
			return pos==text.length() ? -1 : text.charAt(pos);
		}
		void take(char ch){
			if(peek()!=ch){throw new IllegalArgumentException("Expected "+ch+" at "+pos);}
			pos++;
		}
		String string(){
			take('"');
			StringBuilder sb=new StringBuilder();
			while(pos<text.length()){
				char ch=text.charAt(pos++);
				if(ch=='"'){return sb.toString();}
				if(ch<32){throw new IllegalArgumentException("Unescaped control");}
				if(ch=='\\'){
					if(pos==text.length()){throw new IllegalArgumentException("Truncated escape");}
					ch=text.charAt(pos++);
					switch(ch){
						case '"': case '\\': case '/': break;
						case 'b': ch='\b'; break;
						case 'f': ch='\f'; break;
						case 'n': ch='\n'; break;
						case 'r': ch='\r'; break;
						case 't': ch='\t'; break;
						case 'u':
							if(pos+4>text.length()){throw new IllegalArgumentException("Truncated unicode escape");}
							ch=(char)Integer.parseInt(text.substring(pos, pos+4),16); pos+=4; break;
						default: throw new IllegalArgumentException("Bad escape");
					}
				}
				sb.append(ch);
			}
			throw new IllegalArgumentException("Unterminated string");
		}
		final String text;
		int pos;
	}

	/** In-memory exchange with a deterministic output-failure option. */
	private static class Exchange extends HttpExchange {
		Exchange(String path, final boolean fail){
			uri=URI.create(path);
			output=new OutputStream(){
				@Override public void write(int b) throws IOException{
					if(fail){throw new IOException("Injected output failure");}
					bytes.write(b);
				}
				@Override public void close(){outputClosed=true;}
			};
		}
		@Override public Headers getRequestHeaders(){return request;}
		@Override public Headers getResponseHeaders(){return response;}
		@Override public URI getRequestURI(){return uri;}
		@Override public String getRequestMethod(){return "GET";}
		@Override public HttpContext getHttpContext(){return null;}
		@Override public void close(){closed=true; outputClosed=true;}
		@Override public InputStream getRequestBody(){return input;}
		@Override public OutputStream getResponseBody(){return output;}
		@Override public void sendResponseHeaders(int status, long length){code=status; declared=length;}
		@Override public InetSocketAddress getRemoteAddress(){return new InetSocketAddress("127.0.0.1",1);}
		@Override public InetSocketAddress getLocalAddress(){return new InetSocketAddress("127.0.0.1",2);}
		@Override public String getProtocol(){return "HTTP/1.1";}
		@Override public Object getAttribute(String name){return null;}
		@Override public void setAttribute(String name, Object value){}
		@Override public void setStreams(InputStream in, OutputStream out){throw new UnsupportedOperationException();}
		@Override public HttpPrincipal getPrincipal(){return null;}
		@Override public int getResponseCode(){return code;}
		final URI uri;
		final Headers request=new Headers(), response=new Headers();
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		final InputStream input=new ByteArrayInputStream(new byte[0]);
		final OutputStream output;
		int code;
		long declared;
		boolean closed, outputClosed;
	}
}
