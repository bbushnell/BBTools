package test.barcode;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.lang.reflect.Field;
import java.net.InetAddress;
import java.net.InetSocketAddress;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;

import barcode.Barcode;
import barcode.CountBarcodes2;
import barcode.DemuxData;
import barcode.NovaDemux;
import barcode.PCRMatrix;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import sketch.Sketch;

/**
 * Records one real CLI request without constructing a matrix or contacting a public server.
 * Run each CLI/type/sum or alias-order combination in a fresh JVM. Allowed PROB
 * requests must reach the recorder with the requested sum mode, while explicit
 * non-PROB requests must fail before any network request or map write.
 * @author Brian Bushnell, Sucrose
 */
public final class DemuxCallerSettingsTest{

	/** Requires mode, cli, type, sum and a new dir; Count2 input defaults to counts, or select input=fastq. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, DemuxCallerSettingsTest.class, false);
		try{
			String mode=null, cli=null, type=null, sum=null, directory=null, input=null,
					coding="a48", expect=null, expectedType=null, expectedSum=null;
			boolean defaultMatrixArg=true;
			final ArrayList<String> matrixArgs=new ArrayList<String>();
			final ArrayList<String> errorMust=new ArrayList<String>();
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(key.equals("cli")){cli=value;}
				else if(key.equals("type")){type=value;}
				else if(key.equals("sum")){sum=value;}
				else if(key.equals("dir")){directory=value;}
				else if(key.equals("input")){input=value;}
				else if(key.equals("coding")){coding=value;}
				else if(key.equals("expect")){expect=value;}
				else if(key.equals("expectedtype")){expectedType=value;}
				else if(key.equals("expectedsum")){expectedSum=value;}
				else if(key.equals("matrixarg")){matrixArgs.add(value);}
				else if(key.equals("errormust")){errorMust.add(value);}
				else if(key.equals("defaultmatrix")){defaultMatrixArg=Parse.parseBoolean(value);}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			check("baseline".equals(mode) || "fixed".equals(mode), "Set mode=baseline or fixed.");
			check("count2".equals(cli) || "nova".equals(cli), "Set cli=count2 or nova.");
			if(input==null){input="count2".equals(cli) ? "counts" : "fastq";}
			check("counts".equals(input) || "fastq".equals(input), "Set input=counts or fastq.");
			check(!"nova".equals(cli) || "fastq".equals(input), "Nova uses FASTQ input only.");
			check("hdist".equals(type) || "prob".equals(type) || "tile".equals(type), "Set type=hdist, prob or tile.");
			check("t".equals(sum) || "f".equals(sum), "Set sum=t or f.");
			check("a48".equals(coding), "This twelve-case fixture uses coding=a48.");
			if(expect==null){expect=("fixed".equals(mode) && !"prob".equals(type) ? "reject" : "request");}
			check("request".equals(expect) || "reject".equals(expect), "Set expect=request or reject.");
			if(expectedType==null){expectedType=("baseline".equals(mode) ? "prob" : type);}
			if(expectedSum==null){expectedSum=("fixed".equals(mode) && "t".equals(sum) ? "t" : "f");}
			check("reject".equals(expect) || "prob".equals(expectedType),
					"Allowed remote requests should be expected as PROB mode.");
			check("t".equals(expectedSum) || "f".equals(expectedSum), "Set expectedsum=t or f.");
			check(directory!=null, "Set dir to a nonexistent directory.");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Caller-settings fixture requires -ea for production invariants.");
			check(!probLoaded(), "Use the public-only classpath; probability implementation must be absent.");
			Shared.setThreads(1);
			final Path dir=Files.createDirectory(Paths.get(directory));
			if(input.equals("counts")){
				//Use CountBarcodes2UseServerTest's counts route; FASTQ metadata failure is a separate finding.
				write(dir.resolve("counts.tsv"), "tile".equals(type) ? "ACGT1101\t1\n" : "ACGT\t1\n");
			}else{
				write(dir.resolve("input.fq"), "@INST:1:FC:1:1101:100:100 1:N:0:ACGT\nACGT\n+\nIIII\n");
			}
			write(dir.resolve("expected.txt"), "ACGT\n");
			final String[] command=command(dir, cli, type, sum, coding, input, matrixArgs, defaultMatrixArg);
			write(dir.resolve("cli.args"), String.join("\n", command)+"\n");
			final Recorder recorder=new Recorder(dir);
			final Throwable cliFailure=runCli(dir, cli, command, recorder);
			if("reject".equals(expect)){
				check(cliFailure instanceof IllegalArgumentException,
						"Expected a non-PROB option rejection before contacting the recorder, observed "+cliFailure);
				if(errorMust.isEmpty()){errorMust.add("server"); errorMust.add("probability");}
				try(java.io.PrintStream failureOut=new java.io.PrintStream(
						Files.newOutputStream(dir.resolve("failure.txt"), StandardOpenOption.CREATE_NEW), false, "UTF-8")){
					cliFailure.printStackTrace(failureOut);
				}
				verifyFailure(cliFailure, errorMust);
				recorder.verifyRejected();
				check(!Files.exists(dir.resolve("map.tsv")), "Rejected request must not write map.tsv.");
			}else{
				if(cliFailure!=null){throw new AssertionError("Expected allowed PROB request to succeed.", cliFailure);}
				recorder.verify(mode, type, sum, expectedType, expectedSum);
				check(Arrays.equals("ACGT\tACGT\t1\t0\n".getBytes(StandardCharsets.US_ASCII),
						Files.readAllBytes(dir.resolve("map.tsv"))),
						"CLI map must contain exactly ACGT->ACGT, count1, distance0.");
			}
			check(!probLoaded(), "Server-mode CLI must not load an optional probability implementation.");
			System.out.println("DEMUX_CALLER_SETTINGS mode="+mode+" cli="+cli+" input="+input+" requested_type="+type+
					" requested_sum="+sum+" expect="+expect+
					" observed_type="+(recorder.decoded==null ? -1 : recorder.decoded.type)+
					" observed_sum="+(recorder.decoded==null ? false : recorder.decoded.hdistSum)+
					" coding="+coding+" requests="+recorder.requests);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Uses real parser flags and one ACGT observation; tile mode carries tile1101 in either input. */
	private static String[] command(final Path dir, final String cli, final String type,
			final String sum, final String coding, final String input, final ArrayList<String> matrixArgs,
			final boolean defaultMatrixArg){
		final ArrayList<String> args=new ArrayList<String>();
		args.add(input.equals("counts") ? "countsin="+dir.resolve("counts.tsv") : "in="+dir.resolve("input.fq"));
		args.add("expected="+dir.resolve("expected.txt"));
		args.add("mapout="+dir.resolve("map.tsv"));
		args.add("useserver=t");
		if(matrixArgs.isEmpty()){
			if(defaultMatrixArg){args.add("matrixtype="+type);}
		}else{args.addAll(matrixArgs);}
		args.add("hdistsum="+sum);
		args.add("coding="+coding);
		args.add("t=1");
		args.add("matrixthreads=1");
		args.add("interleaved=f");
		args.add("qin=33");
		args.add("ow=f");
		args.add("pigz=f");
		args.add("unpigz=f");
		args.add("loglog=f");
		args.add("trimreadcomments=f");
		if(cli.equals("nova")){args.add("nosplit=t");}
		return args.toArray(new String[args.size()]);
	}

	/** Owns the loopback listener, executor, deadline and temporary Shared URL override. */
	private static Throwable runCli(final Path dir, final String cli, final String[] command,
			final Recorder recorder) throws Exception{
		final Field field=Shared.class.getDeclaredField("demuxServer");
		field.setAccessible(true);
		final Object previous=field.get(null);
		final HttpServer server=HttpServer.create(new InetSocketAddress(InetAddress.getByName("127.0.0.1"), 0), 0);
		final ExecutorService executor=Executors.newSingleThreadExecutor();
		final Thread watchdog=new Thread(new Runnable(){
			@Override
			public void run(){
				try{Thread.sleep(60000);}
				catch(final InterruptedException e){return;}
				System.err.println("DemuxCallerSettingsTest exceeded its 60-second finite-case deadline.");
				Runtime.getRuntime().halt(124);
			}
		}, "DemuxCallerSettingsTest-watchdog");
		watchdog.setDaemon(true);
		try{
			server.setExecutor(executor);
			server.createContext("/", recorder);
			server.start();
			watchdog.start();
			final String address="http://127.0.0.1:"+server.getAddress().getPort()+"/";
			write(dir.resolve("address.txt"), address+"\n");
			try{
				field.set(null, address);
				try{
					if(cli.equals("count2")){CountBarcodes2.main(command);}
					else{NovaDemux.main(command);}
					return null;
				}catch(final Throwable t){return t;}
			}finally{field.set(null, previous);}
		}finally{
			try{
				server.stop(0);
				executor.shutdownNow();
				check(executor.awaitTermination(5, TimeUnit.SECONDS), "Recording handler did not terminate after server shutdown.");
			}finally{
				//Keep the finite-case deadline armed if an uninterruptible handler survived shutdown.
				if(executor.isTerminated()){
					watchdog.interrupt();
					watchdog.join(1000);
				}
			}
		}
	}

	/** Single executor owns writes; volatile fields publish handler results to the CLI thread. */
	private static final class Recorder implements HttpHandler{

		Recorder(final Path dir_){
			dir=dir_;
			//ACGT in base5 A/C/G/N/T is 0*125+1*25+2*5+4=39; A48 digit39 is 'W'.
			response="W\tW\n".getBytes(StandardCharsets.US_ASCII);
		}

		/** Always attempts a complete literal response, even when recording/decoding fails. */
		@Override
		public void handle(final HttpExchange exchange) throws IOException{
			try{
				try{
					requests++;
					check("POST".equals(exchange.getRequestMethod()), "Expected CLI POST, observed "+exchange.getRequestMethod());
					check("/demux".equals(exchange.getRequestURI().toString()), "Unexpected CLI request target: "+exchange.getRequestURI());
					raw=readRequest(exchange.getRequestBody());
					Files.write(dir.resolve("request.bin"), raw, StandardOpenOption.CREATE_NEW);
					final ArrayList<byte[]> chunks=new ArrayList<byte[]>(1);
					chunks.add(raw);
					decoded=new DemuxData(chunks);
				}catch(final Throwable t){failure=t;}
				exchange.getResponseHeaders().set("Content-Type", "text/plain; charset=US-ASCII");
				exchange.sendResponseHeaders(200, response.length);
				try(OutputStream out=exchange.getResponseBody()){out.write(response);}
			}catch(final IOException e){failure=e; throw e;}
			finally{exchange.close();}
		}

		/** Baseline requires the exact lost-setting pair, not merely a round-trip mismatch. */
		void verify(final String mode, final String type, final String sum,
				final String expectedTypeName, final String expectedSumName) throws Exception{
			if(failure!=null){throw new AssertionError("Recording server failed while capturing the actual CLI request.", failure);}
			check(requests==1 && decoded!=null && raw!=null, "Expected exactly one completely decoded CLI request; observed "+requests);
			final int expectedType=parseType(expectedTypeName);
			final boolean expectedSum="t".equals(expectedSumName);
			check(decoded.type==expectedType && decoded.hdistSum==expectedSum,
					"Requested "+type+"/"+sum+", expected wire pair "+expectedType+"/"+expectedSum+
					", observed "+decoded.type+"/"+decoded.hdistSum);
			check(decoded.length1==4 && decoded.length2==0, "Expected single4-base barcode metadata, observed "+decoded.length1+"/"+decoded.length2);
			check(decoded.coding==Sketch.A48, "CLI did not request the literal response's A48 encoding.");
			check(decoded.expectedList.size()==1 && decoded.expectedList.contains("ACGT"), "Expected set must contain only ACGT.");
			check(decoded.codeCounts.size()==1, "Request must contain exactly one observed barcode.");
			final Barcode barcode=decoded.codeCounts.iterator().next();
			check("ACGT".equals(barcode.name) && barcode.count()==1, "Request must carry literal ACGT count1.");
			final String wire=new String(raw, StandardCharsets.US_ASCII);
			check(wire.contains("\ntype=")==(expectedType!=PCRMatrix.PROB_TYPE), "Default type flag must be omitted for old servers.");
			check(wire.contains("\nhdistsum=")==expectedSum, "False sum flag must be omitted for old servers.");
			write(dir.resolve("observed.txt"), "type="+decoded.type+"\nhdistsum="+decoded.hdistSum+
					"\nlength1="+decoded.length1+"\nlength2="+decoded.length2+"\ncount=ACGT:1\nexpected=ACGT\n");
		}

		/** Non-PROB remote requests must reject before any HTTP request or map write. */
		void verifyRejected() throws Exception{
			if(failure!=null){throw new AssertionError("Recording server failed despite expected rejection.", failure);}
			check(requests==0 && decoded==null && raw==null, "Rejected request should not reach recorder; observed "+requests);
			check(!Files.exists(dir.resolve("request.bin")), "Rejected request must not write request.bin.");
		}

		private final Path dir;
		private final byte[] response;
		private volatile int requests;
		private volatile byte[] raw;
		private volatile DemuxData decoded;
		private volatile Throwable failure;
	}

	/** Bounds this tiny fixture's captured HTTP body; construction and decoding are never timed. */
	private static byte[] readRequest(final InputStream in) throws IOException{
		final ByteArrayOutputStream out=new ByteArrayOutputStream();
		final byte[] buffer=new byte[1024];
		try(InputStream input=in){
			for(int n; (n=input.read(buffer))>=0; ){
				check(out.size()+n<=65536, "One-barcode request exceeded the 64KiB fixture bound.");
				out.write(buffer, 0, n);
			}
		}
		return out.toByteArray();
	}

	/** Mirrors the existing public-only caller fixture without invoking optional constructors. */
	private static boolean probLoaded() throws Exception{
		final Field field=PCRMatrix.class.getDeclaredField("probLoaded");
		field.setAccessible(true);
		return field.getBoolean(null);
	}

	private static void write(final Path path, final String text) throws IOException{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}
	private static int parseType(final String type){
		if("hdist".equals(type)){return PCRMatrix.HDIST_TYPE;}
		if("prob".equals(type)){return PCRMatrix.PROB_TYPE;}
		if("tile".equals(type)){return PCRMatrix.TILE_TYPE;}
		throw new AssertionError("Unknown matrix type: "+type);
	}
	private static void verifyFailure(final Throwable failure, final ArrayList<String> errorMust){
		final StringBuilder sb=new StringBuilder();
		for(Throwable t=failure; t!=null; t=t.getCause()){sb.append(t.toString()).append('\n');}
		final String observed=sb.toString().toLowerCase(Locale.ROOT);
		for(final String term : errorMust){
			check(observed.indexOf(term.toLowerCase(Locale.ROOT))>=0,
					"Expected rejection to mention '"+term+"', observed "+sb);
		}
	}
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	private DemuxCallerSettingsTest(){}
}
