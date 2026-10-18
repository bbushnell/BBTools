package test.ml;

import ml.RegCalEval;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Locale;

import structures.ByteBuilder;

/** Checks empty validation data and bin-count rejection through the real evaluator.
 * Retains literal inputs, the linear model, ranking output and any rejection.
 * Ordinary controls preserve skipped malformed rows and positive bin counts.
 * Also checks rejection of nonfinite targets without changing finite target ranges.
 * @author Shinobu */
public final class RegCalEvalInputTest{

	/** Runs one case in a new artifact directory; assertions below are unconditional. */
	public static void main(String[] args) throws Exception{
		if(args.length!=2){throw new IllegalArgumentException("Expected case and new artifact directory");}
		final String name=args[0];
		final boolean empty=name.equals("empty"), malformed=name.equals("malformed_only");
		final boolean zero=name.equals("bins_zero"), negative=name.equals("bins_negative");
		final boolean nonfinite=name.equals("target_nan") || name.equals("target_posinf") || name.equals("target_neginf");
		final boolean outside=name.equals("target_outside_range");
		final boolean reject=empty || malformed || zero || negative || nonfinite;
		if(!reject && !name.equals("ordinary") && !name.equals("one_bin") && !name.equals("skip_then_valid") &&
				!name.equals("skip_bad_target") && !outside){
			throw new IllegalArgumentException("Unknown case: "+name);
		}
		final int bins=zero ? 0 : negative ? -1 : name.equals("one_bin") ? 1 : 2;
		final String valid="0\t0.75\n", bad="bad\t0.75\n";
		final String target=name.equals("target_nan") ? "NaN" : name.equals("target_posinf") ? "Infinity" : "-Infinity";
		final String input=nonfinite ? valid+"0\t"+target+"\n" : outside ? "0\t2\n" :
			empty ? "" : malformed ? bad : name.equals("skip_then_valid") ? bad+valid :
			name.equals("skip_bad_target") ? "0\tbad\n"+valid : valid;
		final Path root=Files.createDirectory(Paths.get(args[1]));
		final Path nets=Files.createDirectory(root.resolve("nets"));
		final Path validation=root.resolve("validation.tsv");
		final ByteBuilder model=new ByteBuilder();
		model.append("##bbnet\n#version 1\n#concise\n#dense\n"+
			"#density 1\n#blocksize 1\n#seed 1\n#layers 2\n#dims 1 1\n#coding A48\nC2 LINEAR ");
		model.appendFloatA48(0.25f).space().appendFloatA48(0f).nl();
		Files.write(nets.resolve("net_fixture.bbnet"), model.toBytes());
		Files.write(validation, input.getBytes(StandardCharsets.UTF_8));

		final PrintStream originalOut=System.out;
		final Locale originalLocale=Locale.getDefault();
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		RuntimeException rejection=null;
		try(PrintStream captured=new PrintStream(bytes, true, "UTF-8")){
			System.setOut(captured);
			Locale.setDefault(Locale.ROOT);
			RegCalEval.main(new String[]{"nets="+nets, "val="+validation, "bins="+bins});
		}catch(RuntimeException e){
			rejection=e;
		}finally{
			System.setOut(originalOut);
			Locale.setDefault(originalLocale);
		}
		final String output=new String(bytes.toByteArray(), StandardCharsets.UTF_8);
		Files.write(root.resolve("ranking.stdout"), bytes.toByteArray());
		final String diagnostic=rejection==null ? "NONE" : rejection.getClass().getName()+": "+rejection.getMessage();
		Files.write(root.resolve("rejection.txt"), (diagnostic+"\n").getBytes(StandardCharsets.UTF_8));
		originalOut.print(output);
		if(reject){
			final String expected=nonfinite ? "Non-finite target" : empty || malformed ? "No usable validation rows" : "bins must be positive";
			if(!(rejection instanceof IllegalArgumentException) || rejection.getMessage()==null ||
					!rejection.getMessage().contains(expected)){
				throw new AssertionError("Expected input rejection: "+name+"; observed "+diagnostic, rejection);
			}
			if(nonfinite && (!diagnostic.contains(validation.toString()) || !diagnostic.contains(target))){
				throw new AssertionError("Target rejection must identify input and value: "+diagnostic);
			}
			if(!output.isEmpty()){throw new AssertionError("Rejected input emitted a ranking: "+name);}
		}else{
			if(rejection!=null){throw new AssertionError("Ordinary input rejected: "+name, rejection);}
			checkMetrics(output, outside ? 1.75 : 0.5, outside ? 3.0625 : 0.25);
		}
		System.err.println("RegCalEvalInputTest PASS "+name);
	}

	/** Literal one-row oracles: prediction .25 against targets .75 or 2. */
	private static void checkMetrics(final String output, final double ece, final double mse){
		int rows=0;
		for(String line : output.split("\n")){
			final String[] fields=line.trim().split("\\s+");
			if(fields.length>0 && fields[0].equals("net_fixture")){
				if(fields.length!=4 || Double.parseDouble(fields[2])!=ece || Double.parseDouble(fields[3])!=mse){
					throw new AssertionError("Unexpected calibration metrics: "+line);
				}
				rows++;
			}
		}
		if(rows!=1){throw new AssertionError("Expected one model row, observed "+rows);}
	}

	private RegCalEvalInputTest(){}

}
