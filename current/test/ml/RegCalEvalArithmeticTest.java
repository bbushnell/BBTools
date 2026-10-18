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

/** Checks reported calibration MSE through the real evaluator with finite data.
 * Literal linear networks produce known constant predictions or a finite-operand
 * inference overflow that must be rejected; no training runs.
 * Input files and captured output remain available for inspection.
 * @author Shinobu */
public final class RegCalEvalArithmeticTest{

	/** Runs one case in a new directory beneath the caller's result root. */
	public static void main(String[] args) throws Exception{
		if(args.length!=2){throw new IllegalArgumentException("Expected case and new input directory");}
		final String name=args[0];
		final boolean rejectOverflow=name.equals("finite_inference_overflow");
		final String prediction, target;
		final double expected;
		if(name.equals("ordinary")){
			prediction="0.25"; target="0.75"; expected=0.25;
		}else if(name.equals("large_square")){
			prediction="0"; target="1e20"; expected=1e40;
		}else if(name.equals("large_difference")){
			prediction="3e38"; target="-3e38"; expected=3.6e77;
		}else if(rejectOverflow){
			prediction="0"; target="0"; expected=0;
		}else{
			throw new IllegalArgumentException("Unknown case: "+name);
		}
		final Path root=Files.createDirectory(Paths.get(args[1]));
		final Path nets=Files.createDirectory(root.resolve("nets"));
		final Path validation=root.resolve("validation.tsv");
		final ByteBuilder model=new ByteBuilder();
		model.append("##bbnet\n#version 1\n#concise\n#dense\n"+
			"#density 1\n#blocksize 1\n#seed 1\n#layers 2\n#dims 1 1\n#coding A48\nC2 LINEAR ");
		// Native A48 preserves finite float bits without the decimal parser's
		// exponent/long-mantissa limits.
		model.appendFloatA48(Float.parseFloat(prediction)).space()
			.appendFloatA48(rejectOverflow ? 3e38f : 0f).nl();
		Files.write(nets.resolve("net_fixture.bbnet"), model.toBytes());
		Files.write(validation, ((rejectOverflow ? "2" : "0")+"\t"+target+"\n").getBytes(StandardCharsets.UTF_8));

		final PrintStream originalOut=System.out;
		final Locale originalLocale=Locale.getDefault();
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		RuntimeException rejection=null;
		try(PrintStream captured=new PrintStream(bytes, true, "UTF-8")){
			System.setOut(captured);
			Locale.setDefault(Locale.ROOT);
			RegCalEval.main(new String[]{"nets="+nets, "val="+validation, "bins=2"});
		}catch(RuntimeException e){
			if(!rejectOverflow){throw e;}
			rejection=e;
		}finally{
			System.setOut(originalOut);
			Locale.setDefault(originalLocale);
		}
		final String output=new String(bytes.toByteArray(), StandardCharsets.UTF_8);
		originalOut.print(output);
		if(rejectOverflow){
			if(rejection==null){
				if(readMse(output)!=Double.POSITIVE_INFINITY){
					throw new AssertionError("Expected finite operands to overflow inference before calibration");
				}
				throw new AssertionError("Non-finite prediction was not rejected for "+name);
			}
			final String message=rejection.getMessage();
			if(message==null || !message.contains("Non-finite prediction") ||
					!message.contains("net_fixture.bbnet") || !message.contains("validation row 1") ||
					!message.contains("Infinity")){
				throw new AssertionError("Unexpected overflow rejection", rejection);
			}
			if(!output.isEmpty()){throw new AssertionError("Rejected prediction emitted a calibration ranking: "+output);}
			System.err.println("RegCalEvalArithmeticTest PASS "+name);
			return;
		}
		final double observed=readMse(output);
		if(!Double.isFinite(observed)){
			throw new AssertionError("MSE must be finite for "+name+" but observed "+observed);
		}
		// Float-encoded model/data values round slightly; the independent decimal
		// oracle allows that rounding, while retaining a meaningful relative bound.
		if(!(Math.abs(observed-expected)<=expected*1e-6)){
			throw new AssertionError("MSE expected "+expected+" but observed "+observed);
		}
		System.err.println("RegCalEvalArithmeticTest PASS "+name);
	}

	/** Reads the MSE column from exactly one model's displayed result row. */
	private static double readMse(final String output){
		double mse=Double.NaN;
		int rows=0;
		for(String line : output.split("\n")){
			final String[] fields=line.trim().split("\\s+");
			if(fields.length>0 && fields[0].equals("net_fixture")){
				if(fields.length!=4){throw new AssertionError("Unexpected result row: "+line);}
				mse=Double.parseDouble(fields[3]);
				rows++;
			}
		}
		if(rows!=1){throw new AssertionError("Expected one model result, observed "+rows);}
		return mse;
	}

	private RegCalEvalArithmeticTest(){}

}
