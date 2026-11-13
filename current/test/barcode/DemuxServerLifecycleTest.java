package test.barcode;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.Locale;

import com.sun.net.httpserver.HttpServer;

import barcode.Barcode;
import barcode.DemuxData;
import barcode.PCRMatrix;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import sketch.Sketch;
import shared.Shared;

/** Exercises disposable localhost DemuxServer kill and request-logging controls. */
public final class DemuxServerLifecycleTest{

	/** Requires mode=baseline|fixed, case=..., and a new dir=path. */
	public static void main(final String[] args) throws Exception{
		if(childMode(args)){
			runChild(stripChildMode(args));
			return;
		}
		final PreParser pp=new PreParser(args, DemuxServerLifecycleTest.class, false);
		try{
			String mode=null, testCase=null, directory=null, installation=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(key.equals("case")){testCase=value;}
				else if(key.equals("dir")){directory=value;}
				else if(key.equals("installation")){installation=value;}
				else if(!Parser.parseStatic(arg, key, value)){
					throw new IllegalArgumentException("Unknown argument: "+arg);
				}
			}
			check("baseline".equals(mode) || "fixed".equals(mode),
					"Set mode=baseline or mode=fixed.");
			check("kill-omitted".equals(testCase) || "kill-negative".equals(testCase) ||
					"kill-positive".equals(testCase) || "metadata-off".equals(testCase) ||
					"metadata-on".equals(testCase),
					"Set case=kill-omitted, kill-negative, kill-positive, metadata-off, or metadata-on.");
			check(directory!=null, "Set dir to a nonexistent output directory.");
			check(installation!=null, "Set installation to the BBTools root.");
			final Path dir=Files.createDirectory(Paths.get(directory));

			final boolean needsKill=testCase.startsWith("kill-") && !"kill-omitted".equals(testCase);
			final boolean verbose2="metadata-on".equals(testCase);
			final Path launcher=Paths.get(installation, "dev", "demuxserverlifecycletest.sh");
			final ChildServer child=ChildServer.start(dir, launcher, needsKill, verbose2);
			int caseExit=Integer.MIN_VALUE;
			try{
				if(testCase.startsWith("kill-")){
					checkKill(mode, testCase, child);
				}else{
					checkMetadata(mode, testCase, child);
				}
				caseExit=child.exitCode();
			}finally{
				child.cleanup();
			}
			System.out.println("DEMUX_SERVER_LIFECYCLE mode="+mode+" case="+testCase+
					" case_exit="+formatExit(caseExit)+" cleanup_exit="+
					formatExit(child.exitCode())+" failed=0");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void runChild(final String[] args) throws Exception{
		final Constructor<?> constructor=Class.forName("barcode.DemuxServer").
				getDeclaredConstructor(String[].class);
		constructor.setAccessible(true);
		final Object demuxServer=constructor.newInstance((Object)args);
		final Field httpServerField=demuxServer.getClass().getDeclaredField("httpServer");
		httpServerField.setAccessible(true);
		final HttpServer httpServer=(HttpServer)httpServerField.get(demuxServer);
		System.out.println("PORT\t"+httpServer.getAddress().getPort());
		System.out.flush();
	}

	private static void checkKill(final String mode, final String testCase,
			final ChildServer child) throws Exception{
		final String suffix="kill-negative".equals(testCase) ? "BAD" : SECRET;
		final Response response=child.get("/kill/"+suffix);
		if("fixed".equals(mode) && "kill-negative".equals(testCase)){
			check(response.code==403 && response.body.contains("Bad Code"),
					"Bad fixed kill code should return 403/Bad Code, observed "+response);
			check(child.alive(), "Bad fixed kill code must not terminate the server child.");
		}else if("fixed".equals(mode) && "kill-positive".equals(testCase)){
			check(response.code==200 && response.body.contains("Success."),
					"Good fixed kill code should return 200/Success, observed "+response);
			check(child.waitForExit(EXPECTED_KILL_EXIT, 5000),
					"Good fixed kill code must halt the child with exit "+EXPECTED_KILL_EXIT+
					"; observed "+child.exitCode());
		}else{
			check(response.code==444 && response.body.contains("Shush"),
					"Baseline and omitted kill cases should fall through to the base handler, observed "+
					response);
			check(child.alive(), "Fall-through kill case must leave the server child alive.");
		}
	}

	private static void checkMetadata(final String mode, final String testCase,
			final ChildServer child) throws Exception{
		final Response response=child.post("/demux", hdistBody());
		check(response.code==200 && "ACGT\tACGT\n".equals(response.body),
				"Tiny Hamming demux request should return one literal assignment, observed "+response);
		child.awaitStderr("Replied ", 5000);
		final String stderr=child.stderr();
		final boolean demuxBanner=stderr.contains("Demux handler");
		final boolean localMarker=stderr.contains("@Local");
		if("metadata-off".equals(testCase) && "fixed".equals(mode)){
			check(!demuxBanner, "verbose2=f fixed server should not print Demux handler; stderr:\n"+
					stderr);
			check(!localMarker, "verbose2=f fixed server should not dump @Local; stderr:\n"+
					stderr);
		}else{
			check(demuxBanner && localMarker,
					"baseline or verbose2=t server should dump demux request metadata; stderr:\n"+
					stderr);
		}
		check(stderr.contains("Fetched ") && stderr.contains("Parsed ") &&
				stderr.contains("Assigned ") && stderr.contains("Replied "),
				"The BAR-006 fixture preserves existing phase timing logs; stderr:\n"+stderr);
		check(child.alive(), "Demux metadata cases must leave the server child alive for cleanup.");
	}

	private static byte[] hdistBody(){
		final DemuxData data=new DemuxData(4, 0, 0);
		data.type=PCRMatrix.HDIST_TYPE;
		data.coding=Sketch.RAW;
		data.addPolyA=data.addPolyC=data.addPolyG=data.addPolyT=false;
		data.expectedList=new LinkedHashSet<String>();
		data.expectedList.add("ACGT");
		final ArrayList<Barcode> counts=new ArrayList<Barcode>();
		counts.add(new Barcode("ACGT", 2));
		data.codeCounts=counts;
		final ByteArrayOutputStream out=new ByteArrayOutputStream();
		try{
			for(final byte[] chunk : data.encode()){out.write(chunk);}
		}catch(final IOException e){throw new AssertionError(e);}
		return out.toByteArray();
	}

	/** Always enforce fixture oracles independent of assertion configuration. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private DemuxServerLifecycleTest(){}

	private static final String SECRET="SECRET";
	private static final int EXPECTED_KILL_EXIT=1;

	/** One disposable localhost server process owned by this fixture JVM. */
	private static final class ChildServer{

		static ChildServer start(final Path dir, final Path launcher, final boolean kill,
				final boolean verbose2) throws Exception{
			final Path stdout=dir.resolve("server.stdout");
			final Path stderr=dir.resolve("server.stderr");
			final ArrayList<String> command=new ArrayList<String>();
			command.add(launcher.toString());
			command.add("--child");
			command.add("port=0");
			command.add("allowlocalhost=t");
			command.add("addressprefix=127.");
			command.add("verbose=f");
			command.add("verbose2="+verbose2);
			command.add("printip=f");
			command.add("printheaders=f");
			command.add("matrixthreads=1");
			command.add("simd=f");
			if(kill){command.add("kill="+SECRET);}
			Files.write(dir.resolve("server.command"), join(command).getBytes(StandardCharsets.UTF_8),
					StandardOpenOption.CREATE_NEW);
			final ProcessBuilder pb=new ProcessBuilder(command);
			pb.redirectOutput(stdout.toFile());
			pb.redirectError(stderr.toFile());
			final ChildServer child=new ChildServer(pb.start(), stdout, stderr);
			try{
				child.awaitReady();
				return child;
			}catch(final Throwable t){
				child.cleanup();
				if(t instanceof Exception){throw (Exception)t;}
				if(t instanceof Error){throw (Error)t;}
				throw new RuntimeException(t);
			}
		}

		private ChildServer(final Process process_, final Path stdout_, final Path stderr_){
			process=process_;
			stdout=stdout_;
			stderr=stderr_;
			port=-1;
		}

		Response get(final String path) throws IOException{return request("GET", path, null);}

		Response post(final String path, final byte[] body) throws IOException{
			return request("POST", path, body);
		}

		private Response request(final String method, final String path, final byte[] body)
				throws IOException{
			HttpURLConnection connection=null;
			try{
				connection=(HttpURLConnection)new URL("http://127.0.0.1:"+port+path).
						openConnection();
				connection.setRequestMethod(method);
				connection.setConnectTimeout(2000);
				connection.setReadTimeout(5000);
				if(body!=null){
					connection.setDoOutput(true);
					connection.setFixedLengthStreamingMode(body.length);
					try(OutputStream out=connection.getOutputStream()){
						out.write(body);
					}
				}
				final int code=connection.getResponseCode();
				InputStream in=(code<400 ? connection.getInputStream() : connection.getErrorStream());
				if(in==null){in=connection.getInputStream();}
				try(InputStream stream=in){
					return new Response(code, read(stream));
				}
			}finally{
				if(connection!=null){connection.disconnect();}
			}
		}

		void awaitReady() throws Exception{
			final long deadline=System.currentTimeMillis()+10000;
			while(System.currentTimeMillis()<deadline){
				if(!alive()){throw new AssertionError("DemuxServer child exited early with "+exitCode()+
						"; stdout:\n"+stdout()+"\nstderr:\n"+stderr());}
				port=parsePort(stdout());
				if(port>0){return;}
				Thread.sleep(50);
			}
			throw new AssertionError("Timed out waiting for DemuxServer bound port; stdout:\n"+
					stdout()+"\nstderr:\n"+stderr());
		}

		void awaitStderr(final String needle, final long millis) throws Exception{
			final long deadline=System.currentTimeMillis()+millis;
			while(System.currentTimeMillis()<deadline){
				if(stderr().contains(needle)){return;}
				check(alive(), "Server child exited while waiting for stderr marker "+needle+
						"; stdout:\n"+stdout()+"\nstderr:\n"+stderr());
				Thread.sleep(50);
			}
			throw new AssertionError("Timed out waiting for stderr marker "+needle+"; stderr:\n"+
					stderr());
		}

		boolean waitForExit(final int expected, final long millis) throws InterruptedException{
			final long deadline=System.currentTimeMillis()+millis;
			while(alive() && System.currentTimeMillis()<deadline){Thread.sleep(50);}
			return !alive() && exitCode()==expected;
		}

		boolean waitForDead(final long millis) throws InterruptedException{
			final long deadline=System.currentTimeMillis()+millis;
			while(alive() && System.currentTimeMillis()<deadline){Thread.sleep(50);}
			return !alive();
		}

		boolean alive(){
			try{
				process.exitValue();
				return false;
			}catch(final IllegalThreadStateException itse){
				return true;
			}
		}

		int exitCode(){
			try{return process.exitValue();}
			catch(final IllegalThreadStateException itse){return Integer.MIN_VALUE;}
		}

		String stdout() throws IOException{
			return readPath(stdout);
		}

		String stderr() throws IOException{
			return readPath(stderr);
		}

		void cleanup() throws InterruptedException{
			if(alive()){
				process.destroy();
				if(!waitForDead(1000) && alive()){
					process.destroyForcibly();
					if(!waitForDead(2000)){
						throw new AssertionError("Failed to stop owned DemuxServer child.");
					}
				}
			}
		}

		int port;
		final Process process;
		final Path stdout;
		final Path stderr;
	}

	/** Captured HTTP status and body text. */
	private static final class Response{
		Response(final int code_, final String body_){
			code=code_;
			body=body_;
		}
		@Override
		public String toString(){return "code="+code+", body="+body;}
		final int code;
		final String body;
	}

	private static int parsePort(final String stdout){
		final String[] lines=stdout.split("\\R");
		for(final String line : lines){
			if(line.startsWith("PORT\t")){return Integer.parseInt(line.substring(5));}
		}
		return -1;
	}

	private static String readPath(final Path path) throws IOException{
		if(!Files.exists(path)){return "";}
		return new String(Files.readAllBytes(path), StandardCharsets.UTF_8);
	}

	private static String read(final InputStream in) throws IOException{
		final ByteArrayOutputStream out=new ByteArrayOutputStream();
		final byte[] buffer=new byte[8192];
		for(int len=in.read(buffer); len>=0; len=in.read(buffer)){
			out.write(buffer, 0, len);
		}
		return out.toString("UTF-8");
	}

	private static String formatExit(final int exitCode){
		return exitCode==Integer.MIN_VALUE ? "running" : Integer.toString(exitCode);
	}

	private static boolean childMode(final String[] args){
		for(final String arg : args){
			if("child=t".equalsIgnoreCase(arg)){return true;}
		}
		return false;
	}

	private static String[] stripChildMode(final String[] args){
		final ArrayList<String> list=new ArrayList<String>();
		for(final String arg : args){
			if(!"child=t".equalsIgnoreCase(arg)){list.add(arg);}
		}
		return list.toArray(new String[list.size()]);
	}

	private static String join(final ArrayList<String> command){
		final StringBuilder sb=new StringBuilder();
		for(int i=0; i<command.size(); i++){
			if(i>0){sb.append('\n');}
			sb.append(command.get(i));
		}
		sb.append('\n');
		return sb.toString();
	}
}
