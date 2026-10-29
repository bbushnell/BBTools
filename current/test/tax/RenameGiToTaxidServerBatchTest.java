package test.tax;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLDecoder;
import java.net.URLStreamHandler;
import java.net.URLStreamHandlerFactory;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

import shared.Shared;
import shared.Timer;
import tax.RenameGiToTaxid;

/** Reproduces server-mode RenameGiToTaxid GFF prefix and cross-batch read-limit bugs without opening sockets.
 * A JVM-local protocol handler maps accessions to deterministic taxIDs.
 * @author Brian Bushnell, Sucrose */
public final class RenameGiToTaxidServerBatchTest{

	/** Runs finite server-mode RenameGiToTaxid baselines or fixed expectations.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If fixture I/O or reflective tuning fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		String testCase="all";
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else if(arg.startsWith("case=")){testCase=arg.substring("case=".length());}
			else{throw new IllegalArgumentException(arg);}
		}
		final boolean all="all".equals(testCase);
		final boolean gffPrefix="gff-prefix".equals(testCase);
		final boolean readsLimit="reads-limit".equals(testCase);
		if(!all && !gffPrefix && !readsLimit){
			throw new IllegalArgumentException("Unknown case="+testCase);
		}

		URL.setURLStreamHandlerFactory(new URLStreamHandlerFactory(){
			@Override
			public URLStreamHandler createURLStreamHandler(String protocol){
				if(!protocol.equals("taxfixture")){return null;}
				return new URLStreamHandler(){
					@Override
					protected URLConnection openConnection(URL url){return new Reply(url);}
				};
			}
		});

		Path dir=Files.createTempDirectory("renamegi-server-batch");
		try{
			if(all || gffPrefix){
				checkGffPrefixTrueHealthy(dir);
				checkGffPrefixFalse(dir, baseline);
			}
			if(all || readsLimit){
				checkCrossBatchReadLimit(dir, baseline);
			}
		}finally{
			deleteTree(dir);
		}
		System.out.println("PASS renamegi_server_batch baseline="+baseline+
				" case="+testCase+" requests="+opened);
	}

	/** Prevents construction of this standalone fixture. */
	private RenameGiToTaxidServerBatchTest(){}

	/** Confirms prefix=f currently drops every GFF feature column after the seqid.
	 * Fixed mode rejects GFF prefix=f before any server request or output creation.
	 * @param dir temporary fixture directory
	 * @param baseline whether to run this unresolved baseline-only check
	 * @throws Exception If the tool invocation fails */
	private static void checkGffPrefixFalse(Path dir, boolean baseline) throws Exception{
		final Path in=dir.resolve("prefix.gff"), out=dir.resolve("prefix.out.gff");
		Files.write(in, lines(
				"##gff-version 3",
				"##sequence-region ACC1 1 5",
				"ACC1\tProdigal\tCDS\t1\t5\t.\t+\t0\tID=cds1"),
				StandardCharsets.UTF_8);
		if(!baseline){
			opened=0;
			expectGffPrefixFalseRejected(in.toString(), out);

			final Path fa=dir.resolve("prefix.fa"), mixedOut=dir.resolve("prefix-mixed.out.gff");
			Files.write(fa, lines(">ACC1", "A"), StandardCharsets.UTF_8);
			opened=0;
			expectGffPrefixFalseRejected(fa+","+in, mixedOut);
			return;
		}
		opened=0;
		run(in, out, -1, "prefix=f");
		List<String> actual=Files.readAllLines(out, StandardCharsets.UTF_8);
		List<String> expected=lines(
				"##gff-version 3",
				"##sequence-region tid|42 1 5",
				"tid|42|1");
		check(actual.equals(expected), "Unexpected prefix=f GFF output: "+actual);
		check(opened==1, "One GFF batch should make one server request, not "+opened);
	}

	/** Requires GFF prefix=f to fail during construction without touching output or the server.
	 * @param in input path or path list
	 * @param out primary output path that must not be created
	 * @throws Exception If rejection fails */
	private static void expectGffPrefixFalseRejected(String in, Path out) throws Exception{
		try{
			run(in, out, -1, "prefix=f");
		}catch(RuntimeException e){
			check(e.getMessage().equals("GFF input requires prefix=t; prefix=f would drop feature columns."),
					"Unexpected prefix=f rejection: "+e.getMessage());
			check(opened==0, "Rejected GFF prefix=f should not query the server, but opened "+opened);
			check(!Files.exists(out), "Rejected GFF prefix=f created output: "+out);
			return;
		}
		throw new AssertionError("GFF prefix=f should be rejected during construction.");
	}

	/** Confirms prefix=t keeps GFF feature and sequence-region seqids aligned.
	 * @param dir temporary fixture directory
	 * @throws Exception If the tool invocation fails */
	private static void checkGffPrefixTrueHealthy(Path dir) throws Exception{
		final Path in=dir.resolve("prefix-t.gff"), out=dir.resolve("prefix-t.out.gff");
		Files.write(in, lines(
				"##gff-version 3",
				"##sequence-region ACC1 1 5",
				"ACC1\tProdigal\tCDS\t1\t5\t.\t+\t0\tID=cds1"),
				StandardCharsets.UTF_8);
		opened=0;
		run(in, out, -1, "prefix=t");
		List<String> actual=Files.readAllLines(out, StandardCharsets.UTF_8);
		List<String> expected=lines(
				"##gff-version 3",
				"##sequence-region tid|42|ACC1 1 5",
				"tid|42|ACC1\tProdigal\tCDS\t1\t5\t.\t+\t0\tID=cds1");
		check(actual.equals(expected), "Unexpected prefix=t GFF output: "+actual);
		check(opened==1, "One GFF prefix=t batch should make one server request, not "+opened);
	}

	/** Confirms server-mode reads= currently caps only the current batch, not later input.
	 * @param dir temporary fixture directory
	 * @param baseline whether to expect baseline or fixed output
	 * @throws Exception If the tool invocation fails */
	private static void checkCrossBatchReadLimit(Path dir, boolean baseline) throws Exception{
		final Path in=dir.resolve("reads.fa"), out=dir.resolve("reads.out.fa");
		Files.write(in, lines(
				">ACC1", "A",
				">ACC2", "C",
				">ACC3", "G"),
				StandardCharsets.UTF_8);
		opened=0;
		run(in, out, 1, "prefix=t", "reads=1");
		List<String> actual=Files.readAllLines(out, StandardCharsets.UTF_8);
		List<String> expected=baseline ? lines(">tid|42|ACC1", "A", "C", "G") : lines(">tid|42|ACC1", "A");
		check(actual.equals(expected), "Unexpected reads=1 cross-batch output: "+actual);
		check(opened==(baseline ? 3 : 1), "Unexpected reads=1 request count: "+opened);
	}

	/** Invokes RenameGiToTaxid with a fixture URL and optional reflected maxStoredBytes override.
	 * @param in input path
	 * @param out primary output path
	 * @param maxStoredBytes reflected maxStoredBytes value, or negative to keep the default
	 * @param extra extra tool arguments
	 * @throws Exception If construction or invocation fails */
	private static void run(Path in, Path out, long maxStoredBytes, String... extra) throws Exception{
		run(in.toString(), out, maxStoredBytes, extra);
	}

	/** Invokes RenameGiToTaxid with a fixture URL and optional reflected maxStoredBytes override.
	 * @param in input path or path list
	 * @param out primary output path
	 * @param maxStoredBytes reflected maxStoredBytes value, or negative to keep the default
	 * @param extra extra tool arguments
	 * @throws Exception If construction or invocation fails */
	private static void run(String in, Path out, long maxStoredBytes, String... extra) throws Exception{
		Shared.setTaxServer("taxfixture://local/");
		ArrayList<String> args=new ArrayList<String>();
		args.add("in="+in);
		args.add("out="+out);
		args.add("ow=t");
		args.add("server=t");
		args.add("mode=accession");
		for(String s : extra){args.add(s);}
		RenameGiToTaxid tool=new RenameGiToTaxid(args.toArray(new String[args.size()]));
		if(maxStoredBytes>=0){
			Field f=RenameGiToTaxid.class.getDeclaredField("maxStoredBytes");
			f.setAccessible(true);
			f.setLong(tool, maxStoredBytes);
		}
		Method process=RenameGiToTaxid.class.getDeclaredMethod("process", Timer.class);
		process.setAccessible(true);
		process.invoke(tool, new Timer());
	}

	/** Builds a mutable line list for Files helpers and exact expected output.
	 * @param values desired lines
	 * @return lines as an ArrayList */
	private static ArrayList<String> lines(String... values){
		ArrayList<String> list=new ArrayList<String>(values.length);
		for(String s : values){list.add(s);}
		return list;
	}

	/** Checks stay active under -da.
	 * @param condition Required result
	 * @param message Failure explanation */
	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Best-effort recursive cleanup of finite local fixtures.
	 * @param root root to delete */
	private static void deleteTree(Path root){
		if(root==null || !Files.exists(root)){return;}
		try{
			ArrayList<Path> paths=new ArrayList<Path>();
			try(Stream<Path> stream=Files.walk(root)){
				stream.forEach(paths::add);
			}
			for(int i=paths.size()-1; i>=0; i--){Files.deleteIfExists(paths.get(i));}
		}catch(IOException e){
			throw new RuntimeException(e);
		}
	}

	/** In-memory URL connection; no DNS, listener, HTTP client or network is involved. */
	private static final class Reply extends URLConnection{
		/** Binds this response to the requested URL.
		 * @param url Requested fixture URL */
		Reply(URL url){super(url);}
		@Override
		public void connect(){}
		@Override
		public InputStream getInputStream(){
			try{lastPath=URLDecoder.decode(url.getPath(), "UTF-8");}
			catch(java.io.UnsupportedEncodingException e){throw new AssertionError(e);}
			opened++;
			String query=lastPath.substring(lastPath.lastIndexOf('/')+1);
			String[] split=query.split(",");
			StringBuilder sb=new StringBuilder();
			for(String s : split){
				if(sb.length()>0){sb.append(',');}
				if("ACC1".equals(s)){sb.append(42);}
				else if("ACC2".equals(s)){sb.append(43);}
				else if("ACC3".equals(s)){sb.append(44);}
				else{sb.append(-1);}
			}
			final byte[] body=sb.toString().getBytes(StandardCharsets.UTF_8);
			return new ByteArrayInputStream(body);
		}
	}

	/** Decoded path of the most recent fixture request. */
	private static String lastPath;
	/** Count of fixture responses opened by production client code. */
	private static int opened;
}
