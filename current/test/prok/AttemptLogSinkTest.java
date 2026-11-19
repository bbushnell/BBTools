package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import fileIO.ByteStreamWriter;
import prok.AttemptLogSink;
import prok.RrnaResourceIO;

/** Attempt identity is optional, while real zero is a measured value. */
public final class AttemptLogSinkTest {
	public static void main(String[] args)throws Exception{
		final Path out=Files.createTempDirectory("attempt_log_").resolve("attempt.tsv");
		final ByteStreamWriter writer=RrnaResourceIO.writer(out.toString());
		final AttemptLogSink sink=new AttemptLogSink(writer);
		writer.println(AttemptLogSink.HEADER);
		sink.onAttempt(0, 3, 10, 20, 1.5f, false, 10, 20, 1, "ctg\tone",
			"model\none", 123, -1f, "NO_ALIGNMENT_FOUND", "cand0");
		sink.onAttempt(1, 3, 11, 21, 2.5f, true, 12, 22, 0, "ctg2",
			"model2", 124, 0f, "ACCEPTED", "cand1");
		boolean failed=false;
		try{
			sink.onAttempt(2, 3, 12, 22, 3.5f, false, 12, 22, 0, "ctg3",
				"model3", 125, Float.NaN, "BAD", "cand2");
		}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Measured attempt identity must reject nonfinite values");
		failed=false;
		try{
			sink.onAttempt(2, 3, 12, 22, 3.5f, true, 12, 22, 0, "ctg3",
				"model3", 125, -1f, "BAD", "cand3");
		}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Accepted attempts must carry the measured winning identity");
		RrnaResourceIO.finish(writer);
		final java.util.List<String> lines=Files.readAllLines(out, StandardCharsets.US_ASCII);
		check(lines.size()==4 && lines.get(0).equals(AttemptLogSink.SCHEMA)
				&& lines.get(1).equals(AttemptLogSink.COLUMNS),
			"Attempt log schema and column header must be owned by the writer");
		check(lines.get(2).equals("cand0\t0\t3\t10\t20\t10\t20\t1\tctg one\tmodel one\t123\t1.5\t\tNO_ALIGNMENT_FOUND\tfalse"),
			"Unmeasured attempt identity must be an empty optional cell: "+lines.get(2));
		check(lines.get(3).equals("cand1\t1\t3\t11\t21\t12\t22\t0\tctg2\tmodel2\t124\t2.5\t0.0\tACCEPTED\ttrue"),
			"Measured zero identity must remain numeric: "+lines.get(3));
		System.out.println("PASS AttemptLogSinkTest optional_identity zero_identity finite_guard");
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
