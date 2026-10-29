package test.tax;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import tax.GiToTaxid;

/** Finite oracle for GiToTaxid raw-table empty input and overflow handling.
 * @author Jean */
public final class GiToTaxidRawBaselineTest{

	/** Runs baseline or fixed expectations for legacy sparse GI loading.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If fixture construction or cleanup fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		final Path root=Files.createTempDirectory("gi-raw-");
		try{
			checkExplicitTaxidControls(baseline);
			checkExplicitTaxidOverflow(baseline);
			checkEmptyRawInput(root, baseline);
			checkRawLoaderOverflow(root, baseline);
		}finally{
			GiToTaxid.DISK_MODE=false;
			GiToTaxid.unload();
			deleteTempTree(root);
		}
		System.out.println("GI_TO_TAXID_RAW_BASELINE_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private GiToTaxidRawBaselineTest(){}

	/** Checks valid explicit IDs and assertion-dependent malformed-input handling.
	 * @param baseline True to skip tightened malformed-input controls */
	private static void checkExplicitTaxidControls(final boolean baseline){
		final String[] headers={"tid|0|", "tid|7|", "tid|999999999",
				"tid|999999999|", "tid|1000000000|",
				"tid|2147483647|", ">ncbi|0000000000000000000007|",
				"tid|00000000002147483647|"};
		final int[] expected={0, 7, 999999999, 999999999, 1000000000,
				Integer.MAX_VALUE, 7, Integer.MAX_VALUE};
		for(int i=0; i<headers.length; i++){
			final String header=headers[i];
			check(GiToTaxid.parseTaxidNumber(header, '|')==expected[i],
					"Explicit String parser control failed: "+header);
			check(GiToTaxid.getID(header)==expected[i],
					"Explicit String lookup control failed: "+header);
			check(GiToTaxid.getID(toBytes(header))==expected[i],
					"Explicit byte lookup control failed: "+header);
		}
		if(!baseline){
			final String[] malformed={"tid|99999999:|", "tid|123456789a|"};
			for(String header : malformed){
				checkMalformedExplicitTaxid(header);
			}
		}
		check(GiToTaxid.parseTaxidNumber("tid|1000000000 trailing", '|')==1000000000,
				"Explicit String parser failed the 10-digit space terminator");
		check(GiToTaxid.getID("tid|1000000000 trailing")==1000000000,
				"Explicit String lookup failed the 10-digit space terminator");
	}

	/** Verifies that malformed explicit TaxIDs reject or assert.
	 * @param header Header with a non-digit in the explicit decimal field */
	private static void checkMalformedExplicitTaxid(final String header){
		final boolean requireAssertion=GiToTaxidRawBaselineTest.class.desiredAssertionStatus();
		boolean rejected=false;
		int parsed=0;
		try{
			parsed=GiToTaxid.parseTaxidNumber(header, '|');
		}catch(AssertionError expected){
			rejected=true;
		}
		checkMalformedExplicitTaxid(header, "String parser", requireAssertion,
				rejected, parsed);
		rejected=false;
		parsed=0;
		try{
			parsed=GiToTaxid.getID(header);
		}catch(AssertionError expected){
			rejected=true;
		}
		checkMalformedExplicitTaxid(header, "String lookup", requireAssertion,
				rejected, parsed);
		rejected=false;
		parsed=0;
		try{
			parsed=GiToTaxid.getID(toBytes(header));
		}catch(AssertionError expected){
			rejected=true;
		}
		checkMalformedExplicitTaxid(header, "byte lookup", requireAssertion,
				rejected, parsed);
	}

	/** Verifies one malformed result in the active assertion mode.
	 * @param header Header with a non-digit in the explicit decimal field
	 * @param label Parser path label
	 * @param requireAssertion True when assertions must reject the header
	 * @param rejected True if parsing threw the expected assertion
	 * @param parsed Parsed value when assertions are disabled */
	private static void checkMalformedExplicitTaxid(final String header,
			final String label, final boolean requireAssertion,
			final boolean rejected, final int parsed){
		if(requireAssertion){
			check(rejected, "Malformed explicit "+label+
					" did not assert: "+header+", parsed="+parsed);
		}else{
			check(!rejected && parsed==-1,
					"Malformed explicit "+label+" accepted "+header+
					" as "+parsed);
		}
	}

	/** Verifies overflow behavior for explicit TaxID parsers.
	 * @param baseline True to expect current int-wrapped values */
	private static void checkExplicitTaxidOverflow(final boolean baseline){
		final int expected=baseline ? 705032704 : -1;
		final int oneAboveExpected=baseline ? Integer.MIN_VALUE : -1;
		final int oneAboveMaxString=GiToTaxid.parseTaxidNumber(
				"tid|2147483648|", '|');
		final int oneAboveMaxByte=GiToTaxid.getID(toBytes("tid|2147483648|"));
		final int stringOverflow=GiToTaxid.parseTaxidNumber(
				"tid|5000000000|", '|');
		final int byteOverflow=GiToTaxid.getID(toBytes("tid|5000000000|"));
		check(oneAboveMaxString==oneAboveExpected,
				"String one-past-max explicit TaxID parsed to "+
				oneAboveMaxString+" instead of "+oneAboveExpected);
		check(oneAboveMaxByte==oneAboveExpected,
				"byte[] one-past-max explicit TaxID parsed to "+
				oneAboveMaxByte+" instead of "+oneAboveExpected);
		check(stringOverflow==expected,
				"String overflowing explicit TaxID parsed to "+stringOverflow+
				" instead of "+expected);
		check(byteOverflow==expected,
				"byte[] overflowing explicit TaxID parsed to "+byteOverflow+
				" instead of "+expected);
	}

	/** Verifies raw empty-table behavior through initialize().
	 * @param root Temporary fixture root
	 * @param baseline True to expect the empty-input null-slice crash
	 * @throws Exception If file creation fails */
	private static void checkEmptyRawInput(final Path root,
			final boolean baseline) throws Exception{
		GiToTaxid.DISK_MODE=false;
		GiToTaxid.unload();
		final Path table=root.resolve("empty.tsv");
		Files.write(table, new byte[0]);
		try{
			GiToTaxid.initialize(table.toString());
			check(!baseline, "Baseline empty raw input should crash");
		}catch(NullPointerException expected){
			if(!baseline){throw expected;}
		}
		if(!baseline){
			check(GiToTaxid.isInitialized(),
					"Fixed empty raw input should initialize an empty table");
			check(GiToTaxid.getID(0, false)==-2,
					"Fixed empty raw table should have max GI -1");
		}
	}

	/** Verifies raw table GI insertion across the small and long GI table paths.
	 * @param root Temporary fixture root
	 * @param baseline True to expect int-wrapped GI insertion
	 * @throws Exception If file creation fails */
	private static void checkRawLoaderOverflow(final Path root,
			final boolean baseline) throws Exception{
		GiToTaxid.DISK_MODE=false;
		GiToTaxid.unload();
		final Path table=root.resolve("overflow.tsv");
		Files.write(table,
				("A\tA.1\t11\t5\n"+
				"B\tB.1\t7\t4294967299\n"+
				"C\tC.1\t8\t00000000004294967300\n").getBytes(StandardCharsets.US_ASCII));
		GiToTaxid.initialize(table.toString());
		check(GiToTaxid.getID(5)==11,
				"Raw loader should retain an ordinary small GI mapping");
		if(baseline){
			check(GiToTaxid.getID(3)==7,
					"Baseline raw loader should wrap GI 4294967299 to slot 3");
			check(GiToTaxid.getID(4)==8,
					"Baseline raw loader should wrap GI 4294967300 to slot 4");
		}else{
			check(GiToTaxid.getID(3, false)==0,
					"Fixed raw loader should leave low slot 3 empty");
			check(GiToTaxid.getID(4, false)==0,
					"Fixed raw loader should leave low slot 4 empty");
			check(GiToTaxid.getID(4294967299L)==7,
					"Fixed raw loader should retain long GI 4294967299");
			check(GiToTaxid.getID(4294967300L)==8,
					"Fixed raw loader should retain long GI 4294967300");
		}
	}

	/** Encodes a fixture header as ASCII bytes.
	 * @param s Header text
	 * @return ASCII bytes */
	private static byte[] toBytes(final String s){
		return s.getBytes(StandardCharsets.US_ASCII);
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
