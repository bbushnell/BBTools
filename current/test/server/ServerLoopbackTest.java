package test.server;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
import java.net.InetAddress;
import java.net.InetSocketAddress;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;

import server.ServerTools;
import structures.StringNum;

/** Explicitly-started finite loopback checks; transport mode needs no modern stress-client classes. */
public final class ServerLoopbackTest {

	/** Runs a bounded loopback listener only with explicit loopback arguments. */
	public static void main(String[] args) throws Exception{
		if(args.length!=3 || !(args[0].equals("loopback") || args[0].equals("transport")) ||
				!(args[1].equals("baseline") || args[1].equals("fixed"))){
			throw new IllegalArgumentException("Usage: java -Xmx128m test.server.ServerLoopbackTest <loopback|transport> <baseline|fixed> <info.tsv>");
		}
		new ServerLoopbackTest(args[1].equals("fixed"), Paths.get(args[2]), args[0].equals("transport")).run();
	}

	private ServerLoopbackTest(boolean fixed_, Path infoPath_, boolean transportOnly_){
		fixed=fixed_; infoPath=infoPath_; transportOnly=transportOnly_;
	}

	private void run() throws Exception{
		try{
			start();
			checkNullPost();
			checkEmptyPost();
			if(transportOnly){
				expect("ok-request-count", okHandler.requests.get()==(fixed ? 2 : 1), true);
				expect("error-request-count", errorHandler.requests.get()==0, true);
				System.out.println("SERVER_TRANSPORT_LOOPBACK_PASS");
				return;
			}
			checkStressClient("StressTest", true);
			checkStressClient("StressTest2", true);
			checkStressClient("StressTest3", true);
			checkStressClient("StressTest-500", errorAddress(), false, "completed\t6",
					"non2xx\t6", "exceptions\t0");
			checkStressClient("StressTest2-500", errorAddress(), false, "completed\t6",
					"non2xx\t6", "exceptions\t0");
			checkStressClient("StressTest3-500", errorAddress(), false, "completed\t6",
					"non2xx\t6", "exceptions\t0");
			checkStressClient("StressTest-drop", dropAddress(), false, "completed\t0",
					"non2xx\t0", "exceptions\t6");
			checkStressClient("StressTest2-drop", dropAddress(), false, "completed\t0",
					"non2xx\t0", "exceptions\t6");
			checkStressClient("StressTest3-drop", dropAddress(), false, "completed\t0",
					"non2xx\t0", "exceptions\t6");
			expect("ok-request-count", okHandler.requests.get()==(fixed ? 20 : 19), true);
			expect("error-request-count", errorHandler.requests.get()==18, true);
		}finally{
			stop();
		}
	}

	private void start() throws IOException{
		okHandler=new ReplyHandler(200, "ok");
		errorHandler=new ReplyHandler(500, "fail");
		server=HttpServer.create(new InetSocketAddress(InetAddress.getByName("127.0.0.1"), 0), 0);
		server.createContext("/ok", okHandler);
		server.createContext("/fail", errorHandler);
		server.createContext("/drop", new DropHandler());
		executor=Executors.newFixedThreadPool(2);
		server.setExecutor(executor);
		server.start();
		startWatchdog();
		writeInfo(infoPath, server.getAddress().getPort());
	}

	private void stop(){
		if(watchdog!=null){watchdog.interrupt();}
		if(server!=null){server.stop(0);}
		if(executor!=null){executor.shutdownNow();}
	}

	private void startWatchdog(){
		final HttpServer serverRef=server;
		final ExecutorService executorRef=executor;
		watchdog=new Thread(new Runnable() {
			@Override
			public void run(){
				try {
					Thread.sleep(60000);
				} catch (InterruptedException e) {
					return;
				}
				serverRef.stop(0);
				executorRef.shutdownNow();
				System.err.println("ServerLoopbackTest watchdog stopped loopback server after 60s.");
				Runtime.getRuntime().halt(124);
			}
		}, "ServerLoopbackTest-watchdog");
		watchdog.setDaemon(true);
		watchdog.start();
	}

	private void checkNullPost(){
		boolean good=false;
		try{
			final StringNum result=ServerTools.sendAndReceive(null, okAddress());
			good=result!=null && result.n==200 && "ok".equals(result.s);
		}catch(NullPointerException e){good=false;}
		expect("real-httpclient-null-post", good, fixed);
	}

	private void checkEmptyPost(){
		final StringNum result=ServerTools.sendAndReceive(new byte[0], okAddress());
		expect("real-httpclient-empty-post", result!=null && result.n==200 && "ok".equals(result.s), true);
	}

	private void checkStressClient(String name, boolean success) throws Exception{
		final String endpoint=success ? okAddress() : errorAddress();
		final String non2xx=success ? "non2xx\t0" : "non2xx\t6";
		checkStressClient(name, endpoint, success, "completed\t6", non2xx, "exceptions\t0");
	}

	private void checkStressClient(String name, String endpoint, boolean success, String completed,
			String non2xx, String exceptions) throws Exception{
		final String[] args=name.startsWith("StressTest3") ?
				new String[] {endpoint, "6"} : new String[] {endpoint, "6", "2"};
		// Optional modern-only validation dependencies are resolved solely in full mode.
		final String className=name.startsWith("StressTest2") ? "test.server.StressTest2"
				: name.startsWith("StressTest3") ? "test.server.StressTest3" : "test.server.StressTest";
		final String output=runClient(new StressRunner(){
			@Override
			public void run(String[] args){
				try{Class.forName(className).getMethod("main", String[].class).invoke(null, (Object)args);}
				catch(InvocationTargetException e){
					final Throwable cause=e.getCause();
					if(cause instanceof RuntimeException){throw (RuntimeException)cause;}
					if(cause instanceof Error){throw (Error)cause;}
					throw new IllegalStateException("Unexpected checked failure from "+className, cause);
				}catch(ReflectiveOperationException e){
					throw new IllegalStateException("Missing full-scope stress fixture: "+className, e);
				}
			}
		}, args, success);
		expect(name, output.contains(completed) && output.contains(non2xx) &&
				output.contains(exceptions), true);
	}

	private String runClient(StressRunner runner, String[] args, boolean success) throws Exception{
		final PrintStream oldOut=System.out;
		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final PrintStream capture=new PrintStream(buffer, true, "UTF-8");
		Throwable thrown=null;
		String output;
		try{
			System.setOut(capture);
			try{
				runner.run(args);
			}catch(Throwable t){thrown=t;}
		}finally{
			System.setOut(oldOut);
			capture.close();
			output=buffer.toString("UTF-8");
		}
		if(success && thrown!=null){throw new AssertionError("Unexpected stress failure", thrown);}
		if(!success && !(thrown instanceof IllegalStateException)){
			throw new AssertionError("Expected stress IllegalStateException", thrown);
		}
		return output;
	}

	private String okAddress(){return "http://127.0.0.1:"+server.getAddress().getPort()+"/ok";}
	private String errorAddress(){return "http://127.0.0.1:"+server.getAddress().getPort()+"/fail";}
	private String dropAddress(){return "http://127.0.0.1:"+server.getAddress().getPort()+"/drop";}

	private static void expect(String name, boolean actual, boolean expected){
		System.out.println(name+"\t"+actual+"\t"+expected);
		if(actual!=expected){throw new AssertionError(name+": expected "+expected+", observed "+actual);}
	}

	private static void writeInfo(Path path, int port) throws IOException{
		final Path parent=path.toAbsolutePath().getParent();
		if(parent!=null){Files.createDirectories(parent);}
		Files.write(path, ("host\t127.0.0.1\nport\t"+port+"\n").getBytes(StandardCharsets.UTF_8),
				StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.WRITE);
	}

	private static void drain(InputStream input) throws IOException{
		final byte[] buffer=new byte[256];
		while(input.read(buffer)>=0){}
		input.close();
	}

	private interface StressRunner {
		void run(String[] args);
	}

	private static final class DropHandler implements HttpHandler {
		@Override
		public void handle(HttpExchange exchange) throws IOException{
			try{
				drain(exchange.getRequestBody());
			}finally{
				exchange.close();
			}
		}
	}

	private static final class ReplyHandler implements HttpHandler {
		ReplyHandler(int code_, String body_){
			code=code_;
			body=body_.getBytes(StandardCharsets.UTF_8);
		}

		@Override
		public void handle(HttpExchange exchange) throws IOException{
			try{
				drain(exchange.getRequestBody());
				requests.incrementAndGet();
				exchange.sendResponseHeaders(code, body.length);
				OutputStream out=exchange.getResponseBody();
				try{
					out.write(body);
				}finally{
					out.close();
				}
			}finally{
				exchange.close();
			}
		}

		final int code;
		final byte[] body;
		final AtomicInteger requests=new AtomicInteger();
	}

	private final boolean fixed;
	private final boolean transportOnly;
	private final Path infoPath;
	private HttpServer server;
	private ExecutorService executor;
	private Thread watchdog;
	private ReplyHandler okHandler;
	private ReplyHandler errorHandler;
}
