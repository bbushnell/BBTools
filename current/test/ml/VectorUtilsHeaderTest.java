package test.ml;

import ml.VectorUtils;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.stream.Stream;

import shared.Timer;

/** Literal multi-file header cases for {@link VectorUtils}; run with assertions enabled.
 * Checks matching headers, unweighted aliases, conflicting layouts and explicit overrides.
 * @author Brian Bushnell, Shinobu */
public class VectorUtilsHeaderTest{

	/** Runs all four cases or one named case, using a separate temporary directory for each. */
	public static void main(String[] args) throws Exception{
		if(args.length==1){
			run(args[0]);
			return;
		}
		assert(args.length==0) : "This fixture accepts one optional case name";
		run("matching_headers");
		run("unweighted_alias");
		run("conflicting_layouts");
		run("explicit_override");
	}

	/** Runs one isolated input/output example with literal expected bytes. */
	private static void run(String name) throws Exception{
		final Path dir=Files.createTempDirectory("vectorutils-headers-");
		try{
			final String firstHeader="#dims\t2\t1\n";
			final String secondHeader;
			final boolean override=name.equals("explicit_override");
			final boolean rejects=name.equals("conflicting_layouts");
			if(name.equals("matching_headers")){
				secondHeader=firstHeader;
			}else if(name.equals("unweighted_alias")){
				secondHeader="#dims\t2\t1\t0\n";
			}else if(rejects || override){
				secondHeader="#dims\t1\t1\t1\n";
			}else{
				throw new IllegalArgumentException("Unknown fixture case: "+name);
			}
			final Path first=dir.resolve("first.tsv"), second=dir.resolve("second.tsv");
			final Path out=dir.resolve("combined.tsv");
			Files.write(first, (firstHeader+FIRST_ROWS).getBytes(StandardCharsets.US_ASCII));
			Files.write(second, (secondHeader+SECOND_ROWS).getBytes(StandardCharsets.US_ASCII));
			final String[] options=override ?
				new String[]{"in="+first+","+second, "out="+out, "ow=t", "dims=1,1,1"} :
				new String[]{"in="+first+","+second, "out="+out, "ow=t"};
			boolean rejected=false;
			try{
				new VectorUtils(options).process(new Timer());
			}catch(AssertionError e){
				if(!rejects || e.getMessage()==null || !e.getMessage().contains("Conflicting #dims")){
					throw e;
				}
				rejected=true;
			}
			if(rejects){
				assert(rejected) : "Conflicting #dims headers were accepted without an explicit override";
				assert(!Files.exists(out)) : "Layout rejection must precede output creation: "+out;
			}else{
				final String expected=(override ? "#dims\t1\t1\t1\n" : firstHeader)+FIRST_ROWS+SECOND_ROWS;
				final String actual=new String(Files.readAllBytes(out), StandardCharsets.US_ASCII);
				assert(expected.equals(actual)) : "Unexpected output for "+name+":\n"+actual;
			}
			System.err.println("VectorUtilsHeaderTest PASS "+name);
		}finally{
			try(Stream<Path> paths=Files.walk(dir)){
				paths.sorted(Comparator.reverseOrder()).forEach(path -> {
					try{Files.deleteIfExists(path);}
					catch(Exception e){throw new RuntimeException(e);}
				});
			}
		}
	}

	/** Both layouts have three fields, so column-count validation cannot distinguish them. */
	private static final String FIRST_ROWS="10\t20\t0\n11\t21\t1\n";
	private static final String SECOND_ROWS="30\t2\t0\n31\t3\t1\n";
}
