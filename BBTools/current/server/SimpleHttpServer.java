package server;

import java.io.IOException;
import java.io.OutputStream;
import java.net.InetSocketAddress;
import java.util.HashMap;
import java.util.Map;

import com.sun.net.httpserver.Headers;
import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;

import shared.Tools;

/**
 * @author Shijie Yao
 * @date Dec 13, 2016
 *
 */
public class SimpleHttpServer {
	public static int port = 8321;

	public static void main(String[] args) throws Exception {
		HttpServer server = HttpServer.create(new InetSocketAddress(SimpleHttpServer.port), 0);
		server.createContext("/", new GetHandler());
		server.setExecutor(null); // creates a default executor
		server.start();
	}

	static class GetHandler implements HttpHandler {
		@Override
		public void handle(HttpExchange t) throws IOException {
			Headers h = t.getResponseHeaders();
			h.add("Content-Type", "application/json");
			//String query = t.getRequestURI().getQuery(); //the KEY=VAL&KEY=VAL params in URL
			String rparam = t.getRequestURI().toString();   //restful style params, KEY/VAL in URL
			if (rparam.startsWith("/")){
				rparam = rparam.substring(1);
			}
			if (rparam.endsWith("/")){
				rparam = rparam.substring(rparam.length()-1);
			}
			System.out.println(rparam);



			String[] params = rparam.split("/");
			Map<String,String> map = new HashMap<String,String>();  //should be <String, Object>

			if (params.length == 2){
				map.put(params[0], params[1]);
				// fill in other data here

				map.put("tax_id", "654321");
				map.put("organism", "e. coli");


			} else {
				map.put("error", "need restful-style param like gi/123456");
			}

			String response = "{";

			for (Map.Entry<String, String> entry : map.entrySet()){
				response+=Tools.format("\"%s\":\"%s\"", entry.getKey(), entry.getValue());
			}
			response+="}";

			//String response = "{\"name\": \"a json string\"}";
			t.sendResponseHeaders(200, response.length());
			OutputStream os = t.getResponseBody();
			os.write(response.getBytes());
			os.close();
		}
	}
}
