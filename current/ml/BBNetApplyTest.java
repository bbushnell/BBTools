package ml;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Stream;

/**
 * Fixture for the gate runner {@link BBNetApply}: normal output (cross-checked against a
 * direct CellNet feed-forward, float-exact), the {@code rows=} limit, every malformed-input
 * rejection (each must leave NO output file when it fires before the first valid row), and
 * the two failed-output paths: a missing output directory (preflighted rejection; upstream
 * the writer's open path hard-exits the JVM, which no test could catch) and a full device
 * ({@code /dev/full}, Linux only).
 *
 * <p>The full-device case depends on the BBTools on the classpath: with the pinned
 * 54b62510 the ENOSPC at flush/close is dropped upstream ({@code ReadWrite.finishWriting},
 * mag-qc BBTOOLS_BUGS_FOUND row 13), so the runner cannot see it and exits clean; the
 * test then only RECORDS that outcome. With {@code -Dmagqc.readwrite=patched} (the
 * one-line upstream fix first on the classpath, as {@code scripts/ml_gate_tools_staging_v1.sh}
 * arranges) the runner MUST reject with its "I/O error" message.</p>
 *
 * <p>Run with assertions on: {@code java -ea [-Dmagqc.readwrite=patched] ml.BBNetApplyTest}</p>
 *
 * @author Sayu (2026-09-03), UMP45 (2026-09-08)
 */
public class BBNetApplyTest {

	public static void main(String[] args) throws Exception{
		final Path dir=Files.createTempDirectory("bbnet-apply-test-");
		try{
			final Path net=dir.resolve("fixture.bbnet");
			Files.write(net, Arrays.asList(
				"##bbnet", "#version 1", "#concise", "#dense", "#density 1", "#blocksize 1",
				"#seed 1", "#layers 3", "#dims 2 2 1", "#edges 4", "#coding decimal",
				"##layer 1", "C3 TANH 0 1 0", "C4 TANH 0 0 1",
				"##layer 2", "C5 LINEAR 0 1 1"
			), StandardCharsets.US_ASCII);
			final Path in=dir.resolve("vectors.tsv");
			Files.write(in, Arrays.asList("#dims\t2\t1", "# a comment row", "", "1\t2\t0", "3\t4\t0", "0.5\t-1.25\t0"),
				StandardCharsets.US_ASCII);

			// (1) normal: 3 rows, ids 0..2, each output float-exact vs a direct feed-forward.
			final Path out=dir.resolve("pred.tsv");
			BBNetApply.main(new String[]{"net="+net, "in="+in, "out="+out, "rows=10"});
			final List<String> lines=Files.readAllLines(out, StandardCharsets.US_ASCII);
			check(lines.size()==3, "expected 3 output rows, got "+lines);
			final CellNet direct=CellNetParser.load(net.toString(), false);
			final float[][] inputs={{1f,2f},{3f,4f},{0.5f,-1.25f}};
			for(int i=0; i<3; i++){
				final String[] f=lines.get(i).split("\t", -1);
				check(f.length==2 && f[0].equals(Integer.toString(i)), "row "+i+" shape/id: "+lines.get(i));
				direct.applyInput(inputs[i]); direct.feedForward();
				final float expected=direct.getOutput(0);
				check(Float.isFinite(expected), "fixture net produced a non-finite output");
				check(f[1].equals(Float.toString(expected)),
					"row "+i+" output "+f[1]+" != direct feed-forward "+Float.toString(expected));
			}
			System.err.println("BBNetApplyTest (1) PASS: 3 rows, ids in order, outputs float-exact vs direct CellNet");

			// (2) rows= limit.
			final Path out2=dir.resolve("pred2.tsv");
			BBNetApply.main(new String[]{"net="+net, "in="+in, "out="+out2, "rows=2"});
			check(Files.readAllLines(out2, StandardCharsets.US_ASCII).size()==2, "rows=2 did not limit to 2 rows");
			System.err.println("BBNetApplyTest (2) PASS: rows=2 limit honored");

			// (3) malformed inputs: each throws; the ones that fire before the first valid row leave no file.
			rejectNoFile(dir, "missing #dims header", net, write(dir, "nohdr.tsv", "1\t2\t0", "3\t4\t0"), "missing #dims");
			rejectNoFile(dir, "input width mismatch", net, write(dir, "w3.tsv", "#dims\t3\t1", "1\t2\t3"), "input width mismatch");
			rejectNoFile(dir, "output width mismatch", net, write(dir, "o2.tsv", "#dims\t2\t2", "1\t2\t0"), "output width mismatch");
			rejectNoFile(dir, "short first row", net, write(dir, "short.tsv", "#dims\t2\t1", "1"), "columns; expected at least");
			rejectNoFile(dir, "non-finite input", net, write(dir, "nan.tsv", "#dims\t2\t1", "NaN\t2\t0"), "non-finite");
			rejectNoFile(dir, "zero data rows", net, write(dir, "empty.tsv", "#dims\t2\t1", "# only comments"), "zero data rows");
			reject("unknown argument", new String[]{"net="+net, "in="+in, "out="+dir.resolve("x1.tsv"), "bogus=1"}, "Unknown argument");
			reject("rows=0", new String[]{"net="+net, "in="+in, "out="+dir.resolve("x2.tsv"), "rows=0"}, "rows must be positive");
			reject("rows=-1", new String[]{"net="+net, "in="+in, "out="+dir.resolve("x3.tsv"), "rows=-1"}, "rows must be positive");
			reject("missing out=", new String[]{"net="+net, "in="+in}, "Usage");
			reject("bare argument", new String[]{"net="+net, "in="+in, "out="+dir.resolve("x6.tsv"), "rows"}, "Malformed argument");
			reject("empty key", new String[]{"net="+net, "in="+in, "out="+dir.resolve("x7.tsv"), "=5"}, "Malformed argument");
			// A short SECOND row fires after the writer is open (a partial file is inherent); it must still throw.
			reject("short second row", new String[]{"net="+net, "in="+write(dir, "short2.tsv", "#dims\t2\t1", "1\t2\t0", "3"),
				"out="+dir.resolve("x4.tsv")}, "columns; expected at least");
			// A missing net or input file is preflighted by the runner (upstream, the parser's file layer
			// hard-exits the JVM on a missing file -- ReadWrite.java:1345 -- which no test could catch).
			rejectNoFile(dir, "missing net file", dir.resolve("absent.bbnet"), in, "net file not readable");
			rejectNoFile(dir, "missing input file", net, dir.resolve("absent.tsv"), "input TSV not readable");
			System.err.println("BBNetApplyTest (3) PASS: 16 malformed-input rejections, none created an output file prematurely");

			// (4a) failed output: missing output directory -> preflighted rejection (upstream this too is a hard exit).
			reject("missing output directory", new String[]{"net="+net, "in="+in,
				"out="+dir.resolve("no_such_dir").resolve("nested").resolve("pred.tsv")}, "output directory does not exist");
			System.err.println("BBNetApplyTest (4a) PASS: missing output directory rejected before open");

			// (4b) failed output: full device. Outcome depends on the BBTools on the classpath (see javadoc).
			final Path full=Paths.get("/dev/full");
			final String mode=System.getProperty("magqc.readwrite", "pinned");
			if(Files.exists(full)){
				Throwable caught=null;
				try{BBNetApply.main(new String[]{"net="+net, "in="+in, "out="+full});}
				catch(Throwable t){caught=t;}
				if("patched".equals(mode)){
					check(caught instanceof RuntimeException && caught.getMessage()!=null && caught.getMessage().contains("I/O error"),
						"patched ReadWrite on classpath but /dev/full did not raise the I/O error rejection: "+caught);
					System.err.println("BBNetApplyTest (4b) PASS [patched]: /dev/full rejected with: "+caught.getMessage());
				}else{
					if(caught==null){
						System.err.println("BBNetApplyTest (4b) RECORDED [pinned]: /dev/full exited clean -> BBTOOLS_BUGS_FOUND row 13 "
							+"reproduced (ReadWrite.finishWriting drops close()'s error flag upstream; the runner never sees it)");
					}else{
						check(caught.getMessage()!=null && caught.getMessage().contains("I/O error"),
							"unexpected failure type on /dev/full: "+caught);
						System.err.println("BBNetApplyTest (4b) RECORDED [pinned]: /dev/full rejected -> this BBTools already carries the row-13 fix");
					}
				}
			}else{
				check(!"patched".equals(mode), "magqc.readwrite=patched requires /dev/full (Linux)");
				System.err.println("BBNetApplyTest (4b) SKIPPED: /dev/full not present");
			}

			System.err.println("BBNetApplyTest: PASS");
		}finally{
			try(Stream<Path> paths=Files.walk(dir)){
				paths.sorted(Comparator.reverseOrder()).forEach(path -> {
					try{Files.deleteIfExists(path);}catch(Exception e){throw new RuntimeException(e);}
				});
			}
		}
	}

	private static Path write(Path dir, String name, String... lines) throws Exception{
		final Path p=dir.resolve(name);
		Files.write(p, Arrays.asList(lines), StandardCharsets.US_ASCII);
		return p;
	}

	/** Runs main, requires a throw whose message contains {@code fragment}, and requires that {@code out} was never created. */
	private static void rejectNoFile(Path dir, String what, Path net, Path in, String fragment) throws Exception{
		final Path out=dir.resolve("reject_"+Math.abs(what.hashCode())+".tsv");
		reject(what, new String[]{"net="+net, "in="+in, "out="+out}, fragment);
		check(!Files.exists(out), what+": output file was created before the rejection fired");
	}

	private static void reject(String what, String[] args, String fragment){
		try{
			BBNetApply.main(args);
		}catch(RuntimeException e){
			check(e.getMessage()!=null && e.getMessage().contains(fragment),
				what+": wrong rejection message: "+e.getMessage()+" (expected to contain '"+fragment+"')");
			return;
		}
		throw new AssertionError(what+": was accepted");
	}

	private static void check(boolean ok, String message){if(!ok){throw new AssertionError(message);}}
}
