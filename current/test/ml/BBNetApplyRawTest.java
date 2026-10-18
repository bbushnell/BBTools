package test.ml;

import ml.BBNetApply;
import ml.BBNetApplyRaw;
import ml.CellNet;
import ml.CellNetParser;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Stream;

/**
 * Fixture for {@link BBNetApplyRaw}: header/row format ({@code #row<TAB>pred0}, then
 * {@code i<TAB>%.8f}), comment/blank lines skipped, values equal to a direct CellNet
 * feed-forward within the %.8f rounding, agreement with {@link BBNetApply} (the
 * {@code #dims} runner) on the same vectors, malformed-row modes, and the
 * core rejections for zero data rows, wrong argument count and missing files.
 *
 * <p>Run with assertions on: {@code java -ea test.ml.BBNetApplyRawTest}</p>
 *
 * @author UMP45 (2026-09-08)
 */
public class BBNetApplyRawTest{

	/** Runs the aggregate fixture or one retained baseline case selected by name. */
	public static void main(String[] args) throws Exception{
		if(args.length==1){runSelectable(args[0]); return;}
		final Path dir=Files.createTempDirectory("bbnet-apply-raw-test-");
		try{
			final Path net=dir.resolve("fixture.bbnet");
			Files.write(net, Arrays.asList(
				"##bbnet", "#version 1", "#concise", "#dense", "#density 1", "#blocksize 1",
				"#seed 1", "#layers 3", "#dims 2 2 1", "#edges 4", "#coding decimal",
				"##layer 1", "C3 TANH 0 1 0", "C4 TANH 0 0 1",
				"##layer 2", "C5 LINEAR 0 1 1"
			), StandardCharsets.US_ASCII);
			final Path raw=dir.resolve("raw.tsv");
			Files.write(raw, Arrays.asList("# comment", "1 2", "", "3\t4", "0.5 -1.25"), StandardCharsets.US_ASCII);

			// (1) format + values.
			final String captured=capture(new String[]{net.toString(), raw.toString()});
			final String[] lines=captured.split("\n", -1);
			check(lines.length==5 && lines[4].length()==0, "expected header + 3 rows + trailing newline, got: "+captured);
			check("#row\tpred0".equals(lines[0]), "bad header: "+lines[0]);
			final CellNet direct=CellNetParser.load(net.toString(), false);
			final float[][] inputs={{1f, 2f}, {3f, 4f}, {0.5f, -1.25f}};
			for(int i=0; i<3; i++){
				final String[] f=lines[i+1].split("\t", -1);
				check(f.length==2 && f[0].equals(Integer.toString(i)), "row "+i+" shape/id: "+lines[i+1]);
				check(f[1].matches("-?[0-9]+\\.[0-9]{8}"), "row "+i+" value not %.8f: "+f[1]);
				direct.applyInput(inputs[i]); direct.feedForward();
				final double expected=direct.getOutput(0);
				check(Math.abs(Double.parseDouble(f[1])-expected)<=5.0e-9+1e-12,
					"row "+i+" value "+f[1]+" != direct "+expected+" beyond %.8f rounding");
			}
			System.err.println("BBNetApplyRawTest (1) PASS: header, 3 rows, %.8f values equal direct feed-forward");

			// (2) agreement with the #dims runner on the same vectors.
			final Path dims=dir.resolve("dims.tsv");
			Files.write(dims, Arrays.asList("#dims\t2\t1", "1\t2\t0", "3\t4\t0", "0.5\t-1.25\t0"), StandardCharsets.US_ASCII);
			final Path pred=dir.resolve("pred.tsv");
			BBNetApply.main(new String[]{"net="+net, "in="+dims, "out="+pred, "rows=3"});
			final List<String> dimsOut=Files.readAllLines(pred, StandardCharsets.US_ASCII);
			check(dimsOut.size()==3, "BBNetApply rows: "+dimsOut);
			for(int i=0; i<3; i++){
				// The raw runner prints the float at %.8f; the #dims runner prints Float.toString (shortest
				// text that round-trips the SAME float). Parse the #dims text back to that float and widen,
				// then the two can only differ by the %.8f rounding.
				final double a=Double.parseDouble(lines[i+1].split("\t", -1)[1]);
				final double b=Float.parseFloat(dimsOut.get(i).split("\t", -1)[1]);
				check(Math.abs(a-b)<=5.0e-9+1e-12, "row "+i+": raw runner "+a+" vs #dims runner (as float) "+b);
			}
			System.err.println("BBNetApplyRawTest (2) PASS: raw runner agrees with ml.BBNetApply on 3 rows");

			// (3) rejections.
			final Path empty=dir.resolve("empty.tsv");
			Files.write(empty, Arrays.asList("# nothing", ""), StandardCharsets.US_ASCII);
			reject("zero data rows", new String[]{net.toString(), empty.toString()}, "zero data rows");
			reject("wrong argument count", new String[]{net.toString()}, "Usage");
			reject("missing net file", new String[]{dir.resolve("absent.bbnet").toString(), raw.toString()}, "net file not readable");
			reject("missing rows file", new String[]{net.toString(), dir.resolve("absent.tsv").toString()}, "rows file not readable");
			System.err.println("BBNetApplyRawTest (3) PASS: zero rows, bad arg count, missing net, missing rows file rejected");

			System.err.println("BBNetApplyRawTest: PASS");
		}finally{
			try(Stream<Path> paths=Files.walk(dir)){
				paths.sorted(Comparator.reverseOrder()).forEach(path -> {
					try{Files.deleteIfExists(path);}catch(Exception e){throw new RuntimeException(e);}
				});
			}
		}
	}

	/** Runs one retained malformed/whitespace case so baseline failures cannot mask later rows. */
	private static void runSelectable(String name) throws Exception{
		final Path dir=Files.createTempDirectory("bbnet-apply-raw-case-");
		try{
			final Path net=dir.resolve("fixture.bbnet");
			Files.write(net, Arrays.asList(
				"##bbnet", "#version 1", "#concise", "#dense", "#density 1", "#blocksize 1",
				"#seed 1", "#layers 3", "#dims 2 2 1", "#edges 4", "#coding decimal",
				"##layer 1", "C3 TANH 0 1 0", "C4 TANH 0 0 1",
				"##layer 2", "C5 LINEAR 0 1 1"
			), StandardCharsets.US_ASCII);
			final String clean=capture(new String[]{net.toString(),
				write(dir, "clean.tsv", "1 2").toString()});
			if(name.equals("spaced_blank")){
				final String out=capture(new String[]{net.toString(),
					write(dir, "spaced_blank.tsv", "   ", "1 2").toString()});
				expectExactSinglePrediction(name, out, clean);
			}else if(name.equals("spaced_comment")){
				final String out=capture(new String[]{net.toString(),
					write(dir, "spaced_comment.tsv", "   # comment", "1 2").toString()});
				expectExactSinglePrediction(name, out, clean);
			}else if(name.equals("nan_first")){
				rejectNoOutput(name, new String[]{net.toString(),
					write(dir, "nan_first.tsv", "NaN 2").toString()}, "non-finite");
			}else if(name.equals("infinity_first")){
				rejectNoOutput(name, new String[]{net.toString(),
					write(dir, "infinity_first.tsv", "Infinity 2").toString()}, "non-finite");
			}else if(name.equals("short_first")){
				rejectNoOutput(name, new String[]{net.toString(),
					write(dir, "short_first.tsv", "1").toString()}, "input width");
			}else if(name.equals("wide_first")){
				rejectNoOutput(name, new String[]{net.toString(),
					write(dir, "wide_first.tsv", "1 2 3").toString()}, "input width");
			}else if(name.equals("short_later")){
				final String out=rejectCapture(name, new String[]{net.toString(),
					write(dir, "short_later.tsv", "1 2", "3").toString()}, "input width");
				final String[] lines=out.split("\n", -1);
				check(lines.length==3 && lines[0].equals("#row\tpred0") &&
					lines[1].startsWith("0\t") && lines[2].length()==0,
					name+": rejection after one valid row should not emit a stale-tail second row: "+out);
			}else{
				throw new IllegalArgumentException("Unknown fixture mode: "+name);
			}
			System.err.println("BBNetApplyRawTest PASS "+name);
		}finally{
			try(Stream<Path> paths=Files.walk(dir)){
				paths.sorted(Comparator.reverseOrder()).forEach(path -> {
					try{Files.deleteIfExists(path);}catch(Exception e){throw new RuntimeException(e);}
				});
			}
		}
	}

	/** Temporarily redirects System.out to capture BBNetApplyRaw; requires no concurrent console writers. */
	static String capture(String[] args) throws Exception{
		final PrintStream old=System.out;
		final ByteArrayOutputStream bo=new ByteArrayOutputStream();
		final PrintStream ps=new PrintStream(bo, true, "UTF-8");
		System.setOut(ps);
		try{BBNetApplyRaw.main(args);}
		finally{System.setOut(old); ps.flush();}
		return new String(bo.toByteArray(), StandardCharsets.UTF_8);
	}

	/** Writes one US-ASCII whitespace-vector fixture file in the temp directory. */
	private static Path write(Path dir, String name, String... lines) throws Exception{
		final Path p=dir.resolve(name);
		Files.write(p, Arrays.asList(lines), StandardCharsets.US_ASCII);
		return p;
	}

	/** Requires a whitespace/comment case to exactly match a clean single-row control. */
	private static void expectExactSinglePrediction(String what, String captured, String expected){
		check(captured.equals(expected),
			what+": output differs from clean one-row control; expected: "+expected+" got: "+captured);
		final String[] lines=captured.split("\n", -1);
		check(lines.length==3 && lines[2].length()==0,
			"expected header + one row + trailing newline, got: "+captured);
		check("#row\tpred0".equals(lines[0]), "bad header: "+lines[0]);
		check(lines[1].startsWith("0\t"), "bad prediction row: "+lines[1]);
	}

	/** Runs main and requires a rejection that leaves stdout empty. */
	private static void reject(String what, String[] args, String fragment) throws Exception{
		final String out=rejectCapture(what, args, fragment);
		check(out.length()==0, what+": expected no stdout, got: "+out);
	}

	/** Temporarily redirects System.out to capture a rejection; requires no concurrent console writers. */
	private static String rejectCapture(String what, String[] args, String fragment) throws Exception{
		final PrintStream old=System.out;
		final ByteArrayOutputStream bo=new ByteArrayOutputStream();
		final PrintStream ps=new PrintStream(bo, true, "UTF-8");
		System.setOut(ps);
		try{
			BBNetApplyRaw.main(args);
		}catch(RuntimeException e){
			final String out=new String(bo.toByteArray(), StandardCharsets.UTF_8);
			check(e.getMessage()!=null && e.getMessage().contains(fragment),
				what+": wrong rejection message: "+e.getMessage()+" (expected '"+fragment+"')");
			return out;
		}finally{
			System.setOut(old);
			ps.flush();
		}
		throw new AssertionError(what+": was accepted; stdout="+
			new String(bo.toByteArray(), StandardCharsets.UTF_8));
	}

	/** Requires a malformed first row to reject before emitting the header. */
	private static void rejectNoOutput(String what, String[] args, String fragment) throws Exception{
		final String out=rejectCapture(what, args, fragment);
		check(out.length()==0, what+": malformed first row emitted stdout: "+out);
	}

	/** Fails one fixture predicate. */
	private static void check(boolean ok, String message){if(!ok){throw new AssertionError(message);}}
}
