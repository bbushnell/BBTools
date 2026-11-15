package test.icecream;

import java.lang.reflect.Field;
import java.lang.reflect.Method;

import icecream.IceCreamFinder;
import icecream.ReformatPacBio;
import json.JsonObject;
import shared.Timer;

/** Literal empty/positive report checks without starting input or worker threads.
 * @author Brian Bushnell, Shinobu
 */
public final class IceCreamReportTest{

	/** Requires an existing tiny input for normal constructor validation.
	 * @param args Exactly in=path
	 * @throws Exception If an existing private report seam changes
	 */
	public static void main(String[] args) throws Exception{
		if(args.length!=1 || !args[0].startsWith("in=")){
			throw new IllegalArgumentException("Use in=existing-tiny-input");
		}
		final String[] options={args[0], "threads=1", "trim=t", "jni=f"};
		final Object[] programs={new IceCreamFinder(options), new ReformatPacBio(options)};
		int failures=0;
		System.out.println("program\tcase\texpected\tobserved\tresult");
		for(Object program : programs){
			final Report empty=new Report(program);
			failures+=percent(program, "empty-filtered", empty, "Reads_Filtered_Pct", "Reads Filtered:", 0.0, "0.000%");
			failures+=percent(program, "empty-trimmed", empty, "Reads_Trimmed_Pct", "Reads Trimmed:", 0.0, "0.000%");
			failures+=percent(program, "empty-bases", empty, "Bases_Trimmed_Pct", "Bases Trimmed:", 0.0, "0.000%");
			final boolean finder=program instanceof IceCreamFinder;
			failures+=percent(program, "empty-wells", empty, finder ? "Bad_ZMWs_Pct" : "Discarded_ZMWs_Pct",
					finder ? "Bad ZMWs:" : "Discarded ZMWs:", 0.0, "0.000%");
			if(finder){failures+=average(program, "empty-average", empty, null, "NA");}
			else{
				failures+=observe(program, "empty-wells-out-text", "0.00%", empty.lastToken("ZMWs Out:"));
				failures+=observe(program, "empty-reads-out-text", "0.00%", empty.lastToken("Reads Out:"));
				failures+=observe(program, "empty-bases-out-text", "0.00%", empty.lastToken("Bases Out:"));
				set(program, "ZMWsOut", 1L);
			}

			set(program, "readsProcessed", 4L); set(program, "readsOut", 1L);
			set(program, "basesProcessed", 16L); set(program, "basesOut", 4L);
			set(program, "ZMWs", 2L); set(program, "readsTrimmed", 2L); set(program, "basesTrimmed", 4L);
			final Report positive=new Report(program);
			failures+=percent(program, "positive-filtered", positive, "Reads_Filtered_Pct", "Reads Filtered:", 75.0, "75.000%");
			failures+=percent(program, "positive-trimmed", positive, "Reads_Trimmed_Pct", "Reads Trimmed:", 50.0, "50.000%");
			failures+=percent(program, "positive-bases", positive, "Bases_Trimmed_Pct", "Bases Trimmed:", 25.0, "25.000%");
			if(!finder){
				failures+=observe(program, "positive-wells-out-text", "50.00%", positive.lastToken("ZMWs Out:"));
				failures+=observe(program, "positive-reads-out-text", "25.00%", positive.lastToken("Reads Out:"));
				failures+=observe(program, "positive-bases-out-text", "25.00%", positive.lastToken("Bases Out:"));
			}
			if(finder){
				failures+=average(program, "no-comparison-average", positive, null, "NA");
				set(program, "iceCreamRatio", 1.0); set(program, "ratiosCounted", 2L);
				failures+=average(program, "positive-average", new Report(program), 0.5, "0.500");
			}

			set(program, "trimReads", false); set(program, "trimPolyA", true);
			failures+=percent(program, "poly-only-trimmed", new Report(program),
					"Reads_Trimmed_Pct", "Reads Trimmed:", 50.0, "50.000%");
		}
		if(failures!=0){throw new AssertionError(failures+" report observations differ from literal expectations");}
	}

	/** Checks both public encodings of a percentage against independently chosen counts. */
	private static int percent(Object program, String name, Report report, String key, String label,
			double expected, String expectedText){
		return observe(program, name+"-json", expected, report.json.getDouble(key))+
				observe(program, name+"-text", expectedText, report.lastToken(label));
	}

	/** Undefined averages retain their JSON field with null rather than inventing a zero observation. */
	private static int average(Object program, String name, Report report, Double expected, String text){
		final String key="Avg_Score_Ratio_IR";
		final Object value=report.json.toString().contains('"'+key+'"') ? report.json.getDouble(key) : "missing field";
		return observe(program, name+"-json", expected, value)+
				observe(program, name+"-text", text, report.lastToken("Avg Score Ratio:"));
	}

	/** Assigns existing counters/settings, without exercising the classifier or input pipeline. */
	private static void set(Object program, String name, Object value) throws Exception{
		final Field field=program.getClass().getDeclaredField(name);
		field.setAccessible(true);
		field.set(program, value);
	}

	/** Prints all mismatches so the baseline exposes each independent report defect. */
	private static int observe(Object program, String name, Object expected, Object actual){
		final boolean pass=expected==null ? actual==null : expected.equals(actual);
		System.out.println(program.getClass().getSimpleName()+"\t"+name+"\t"+expected+"\t"+actual+"\t"+(pass ? "PASS" : "FAIL"));
		return pass ? 0 : 1;
	}

	/** Captures both reports from the same stopped timer and counter state. */
	private static final class Report{
		Report(Object program) throws Exception{
			final Timer timer=new Timer();
			timer.stop();
			final Method toText=program.getClass().getDeclaredMethod("toText", Timer.class);
			final Method toJson=program.getClass().getDeclaredMethod("toJson", Timer.class);
			toText.setAccessible(true); toJson.setAccessible(true);
			text=toText.invoke(program, timer).toString();
			json=(JsonObject)toJson.invoke(program, timer);
		}
		String lastToken(String prefix){
			for(String line : text.split("\n")){
				if(line.startsWith(prefix)){
					final String[] words=line.trim().split("\\s+");
					return words[words.length-1];
				}
			}
			return "missing line";
		}
		final String text;
		final JsonObject json;
	}

	/** Static-only fixture. */
	private IceCreamReportTest(){}
}
