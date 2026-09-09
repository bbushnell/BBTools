package ml;

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
 * {@code #dims} runner) on the same vectors, and the two rejections (zero data rows,
 * wrong argument count).
 *
 * <p>Run with assertions on: {@code java -ea ml.BBNetApplyRawTest}</p>
 *
 * @author UMP45 (2026-09-08)
 */
public class BBNetApplyRawTest {

	public static void main(String[] args) throws Exception{
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
			final float[][] inputs={{1f,2f},{3f,4f},{0.5f,-1.25f}};
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

	/** Runs BBNetApplyRaw.main with stdout captured; restores System.out even on failure. */
	static String capture(String[] args) throws Exception{
		final PrintStream old=System.out;
		final ByteArrayOutputStream bo=new ByteArrayOutputStream();
		final PrintStream ps=new PrintStream(bo, true, "UTF-8");
		System.setOut(ps);
		try{BBNetApplyRaw.main(args);}
		finally{System.setOut(old); ps.flush();}
		return new String(bo.toByteArray(), StandardCharsets.UTF_8);
	}

	private static void reject(String what, String[] args, String fragment) throws Exception{
		final PrintStream old=System.out;
		System.setOut(new PrintStream(new ByteArrayOutputStream(), true, "UTF-8"));
		try{
			BBNetApplyRaw.main(args);
		}catch(RuntimeException e){
			check(e.getMessage()!=null && e.getMessage().contains(fragment),
				what+": wrong rejection message: "+e.getMessage()+" (expected '"+fragment+"')");
			return;
		}finally{
			System.setOut(old);
		}
		throw new AssertionError(what+": was accepted");
	}

	private static void check(boolean ok, String message){if(!ok){throw new AssertionError(message);}}
}
