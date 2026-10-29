package test.tax;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import tax.AccessionToTaxid;
import tax.AnalyzeAccession;

/** Finite oracle for accession comma truncation in disk-table row digitization.
 * @author Jean */
public final class AnalyzeAccessionCommaTest{

	/** Runs baseline or fixed expectations for comma-delimited alternate accessions.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If fixture construction or static cleanup fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		final Path root=Files.createTempDirectory("analyze-accession-comma-");
		try{
			resetAnalyzeAccessionPatterns();
			final Path patterns=root.resolve("patterns.tsv");
			Files.write(patterns, ("#Pattern\tCount\tCombos\tBits\n"+
					"LLL\t1\t17576\t14.10\n").getBytes(StandardCharsets.UTF_8));
			AnalyzeAccession.loadCodeMap(patterns.toString());

			final long primary=AnalyzeAccession.digitize("ABC");
			final long dotted=AnalyzeAccession.digitize(row("ABC.1"));
			final long colon=AnalyzeAccession.digitize(row("ABC:1"));
			final long comma=AnalyzeAccession.digitize(row("ABC,ALT"));
			check(primary>0, "Primary accession should digitize: "+primary);
			check(AccessionToTaxid.parseLineToTaxid(row("ABC,ALT"), (byte)'\t')==123,
					"Comma-row TaxID parse should still read the third column");
			check(dotted==primary, "Dotted row changed digitization: "+dotted+
					" versus "+primary);
			check(colon==primary, "Colon row changed digitization: "+colon+
					" versus "+primary);

			if(baseline){
				check(comma<=0, "Baseline must expose comma row rejection: "+comma);
			}else{
				check(comma==primary, "Comma row should digitize to the primary key: "+
						comma+" versus "+primary);
			}
		}finally{
			resetAnalyzeAccessionPatterns();
			deleteTempTree(root);
		}
		System.out.println("ANALYZE_ACCESSION_COMMA_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private AnalyzeAccessionCommaTest(){}

	/** Constructs a minimal accession2taxid row around a primary accession field.
	 * @param accession Primary accession field
	 * @return Tab-delimited row bytes */
	private static byte[] row(final String accession){
		return (accession+"\t"+accession+".1\t123\t0").getBytes(StandardCharsets.US_ASCII);
	}

	/** Clears AnalyzeAccession's static pattern state for same-JVM fixture hygiene.
	 * @throws Exception If reflection fails */
	private static void resetAnalyzeAccessionPatterns() throws Exception{
		set("codeMap", null);
		set("patternCodeMap", null);
		set("patternList", null);
		set("codeBits0", -1);
		set("codeBits1", -1);
		set("patternShift", null);
		set("patternLowBits", null);
		set("longestPattern", -1);
	}

	/** Assigns one private or public static AnalyzeAccession field.
	 * @param name Field name
	 * @param value Value to assign
	 * @throws Exception If reflection fails */
	private static void set(final String name, final Object value) throws Exception{
		final Field field=AnalyzeAccession.class.getDeclaredField(name);
		field.setAccessible(true);
		field.set(null, value);
	}

	/** Deletes the owned temporary tree, visiting children before parents.
	 * @param root Temporary fixture root
	 * @throws Exception If walking or deletion fails */
	private static void deleteTempTree(final Path root) throws Exception{
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths.sorted(java.util.Comparator.reverseOrder())::iterator){
				Files.deleteIfExists(path);
			}
		}
	}

	/** Keeps fixture checks active under -da.
	 * @param condition Required condition
	 * @param message Failure message */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
