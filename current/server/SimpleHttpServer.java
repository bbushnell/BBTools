package server;

import java.io.IOException;
import java.io.OutputStream;
import java.net.InetSocketAddress;
import java.nio.charset.StandardCharsets;
import java.util.LinkedHashMap;
import java.util.Map;

import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;

/**
 * Demonstration HTTP server returning one requested key/value plus fixed mock data.
 * This is not a taxonomy lookup service. It preserves raw URI text, including
 * percent escapes and query text, and returns status 200 with an error object
 * when the URI does not supply exactly two slash-separated components.
 *
 * @author Shijie Yao
 * @date Dec 13, 2016
 */
public class SimpleHttpServer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructs the stateless demonstration entry point; does not start a listener. */
	public SimpleHttpServer(){}

	/** Starts the demonstration listener on {@link #port}; arguments are unused.
	 * @param args Unused command-line arguments
	 * @throws IOException If binding the listener fails */
	public static void main(String[] args) throws IOException{
		final HttpServer server=HttpServer.create(new InetSocketAddress(port), 0);
		server.createContext("/", new GetHandler());
		server.setExecutor(null);
		server.start();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Handler            ----------------*/
	/*--------------------------------------------------------------*/

	/** Stateless handler; all response data belongs to the current request. */
	static class GetHandler implements HttpHandler{

		/** Creates a handler with no request-specific instance state. */
		GetHandler(){}

		/**
		 * Writes a UTF-8 JSON object and closes the exchange even if writing fails.
		 * The historical handler accepts any HTTP method and does not read its body.
		 * @param exchange Incoming exchange owned and closed by this handler
		 * @throws IOException If sending headers or the response body fails
		 */
		@Override
		public void handle(HttpExchange exchange) throws IOException{
			try{
				String paramsText=exchange.getRequestURI().toString();
				if(paramsText.startsWith("/")){paramsText=paramsText.substring(1);}
				// Fixed [server/SimpleHttpServer#001]: remove the final slash, not the prefix.
				if(paramsText.endsWith("/")){paramsText=paramsText.substring(0, paramsText.length()-1);}
				final String[] params=paramsText.split("/");
				final Map<String, String> fields=new LinkedHashMap<String, String>();
				if(params.length==2){
					fields.put(params[0], params[1]);
					// Deliberate mock values; preserve replacement of these requested keys.
					fields.put("tax_id", "654321");
					fields.put("organism", "e. coli");
				}else{
					fields.put("error", "need restful-style param like gi/123456");
				}

				final StringBuilder json=new StringBuilder(128);
				json.append('{');
				boolean first=true;
				for(Map.Entry<String, String> entry : fields.entrySet()){
					// Fixed [server/SimpleHttpServer#001]: separate every JSON member.
					if(!first){json.append(',');}
					appendQuoted(json, entry.getKey());
					json.append(':');
					appendQuoted(json, entry.getValue());
					first=false;
				}
				json.append('}');
				final byte[] response=json.toString().getBytes(StandardCharsets.UTF_8);
				exchange.getResponseHeaders().set("Content-Type", "application/json; charset=UTF-8");
				exchange.sendResponseHeaders(200, response.length);
				try(OutputStream output=exchange.getResponseBody()){
					output.write(response);
				}
			}finally{
				exchange.close();
			}
		}
	}

	/** Quotes JSON string text without the ASCII-only conversion of ByteBuilder.
	 * @param out Destination builder, retaining its existing contents
	 * @param value Non-null field text; control characters are escaped */
	private static void appendQuoted(StringBuilder out, String value){
		assert(value!=null) : "Demo response fields must have values before JSON serialization";
		out.append('"');
		for(int i=0; i<value.length(); i++){
			final char c=value.charAt(i);
			if(c=='"' || c=='\\'){
				out.append('\\').append(c);
			}else if(c<0x20){
				out.append("\\u00").append(HEX[c>>>4]).append(HEX[c&15]);
			}else{out.append(c);}
		}
		out.append('"');
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Demonstration listener port, set before calling main. */
	public static int port=8321;
	/** Hexadecimal digits for JSON control-character escapes. */
	private static final char[] HEX="0123456789abcdef".toCharArray();
}
