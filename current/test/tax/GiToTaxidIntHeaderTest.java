package test.tax;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import tax.GiToTaxidInt;

/** Finite oracle for GiToTaxidInt String leading-greater-than GI headers,
 * explicit-TaxID overflow, raw-loader GI overflow, and failed same-path
 * initialize retries.
 * @author Jean, Shinobu */
public final class GiToTaxidIntHeaderTest{

	/** Runs baseline or fixed expectations for legacy dense GI parsing and loading.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If reflection setup or temporary file cleanup fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		final Path root=Files.createTempDirectory("gi-int-header-");
		try{
			final int[] table=new int[10];
			table[3]=7;
			table[7]=9;
			ARRAY_FIELD.set(null, table);

			check(GiToTaxidInt.getID(3)==7,
					"Direct GI 3 control failed");
			check(GiToTaxidInt.getID("gi|3|")==7,
					"String GI control should map gi 3 once");
			check(GiToTaxidInt.getID(toBytes(">gi|3|"))==7,
					"byte[] leading-greater-than GI should map gi 3 once");
			check(GiToTaxidInt.parseGiToTaxid(toBytes(">gi|3|"))==7,
					"byte[] parseGiToTaxid should map gi 3 once");

			final int stringHeader=GiToTaxidInt.getID(">gi|3|");
			final int stringParse=GiToTaxidInt.parseGiToTaxid(">gi|3|");
			final int expected=baseline ? 9 : 7;
			check(stringHeader==expected, "String GI header mapped to "+
					stringHeader+" instead of "+expected);
			check(stringParse==expected, "String parseGiToTaxid mapped to "+
					stringParse+" instead of "+expected);

			checkExplicitTaxidControls(baseline);
			checkExplicitTaxidOverflow(baseline);
			final int overflowExpected=baseline ? 705032704 : -1;
			final int stringOverflow=GiToTaxidInt.parseTaxidNumber(
					"tid|5000000000|", '|');
			final int byteOverflow=GiToTaxidInt.getID(toBytes("tid|5000000000|"));
			check(stringOverflow==overflowExpected,
					"String overflowing explicit TaxID parsed to "+stringOverflow+
					" instead of "+overflowExpected);
			check(byteOverflow==overflowExpected,
					"byte[] overflowing explicit TaxID parsed to "+byteOverflow+
					" instead of "+overflowExpected);

			checkFailedInitializeRetry(root, baseline);
			checkRawLoaderOverflow(root, baseline);
		}finally{
			GiToTaxidInt.unload();
			deleteTempTree(root);
		}
		System.out.println("GI_TO_TAXID_INT_HEADER_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private GiToTaxidIntHeaderTest(){}

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
			check(GiToTaxidInt.parseTaxidNumber(header, '|')==expected[i],
					"Explicit String parser control failed: "+header);
			check(GiToTaxidInt.getID(header)==expected[i],
					"Explicit String lookup control failed: "+header);
			check(GiToTaxidInt.getID(toBytes(header))==expected[i],
					"Explicit byte lookup control failed: "+header);
		}
		if(!baseline){
			final String[] malformed={"tid|99999999:|", "tid|123456789a|"};
			for(String header : malformed){
				checkMalformedExplicitTaxid(header);
			}
		}
		check(GiToTaxidInt.parseTaxidNumber("tid|1000000000 trailing", '|')==1000000000,
				"Explicit String parser failed the 10-digit space terminator");
		check(GiToTaxidInt.getID("tid|1000000000 trailing")==1000000000,
				"Explicit String lookup failed the 10-digit space terminator");
	}

	/** Verifies that malformed explicit TaxIDs reject or assert.
	 * @param header Header with a non-digit in the explicit decimal field */
	private static void checkMalformedExplicitTaxid(final String header){
		final boolean requireAssertion=GiToTaxidIntHeaderTest.class.desiredAssertionStatus();
		boolean rejected=false;
		int parsed=0;
		try{
			parsed=GiToTaxidInt.parseTaxidNumber(header, '|');
		}catch(AssertionError expected){
			rejected=true;
		}
		checkMalformedExplicitTaxid(header, "String parser", requireAssertion,
				rejected, parsed);
		rejected=false;
		parsed=0;
		try{
			parsed=GiToTaxidInt.getID(header);
		}catch(AssertionError expected){
			rejected=true;
		}
		checkMalformedExplicitTaxid(header, "String lookup", requireAssertion,
				rejected, parsed);
		rejected=false;
		parsed=0;
		try{
			parsed=GiToTaxidInt.getID(toBytes(header));
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

	/** Verifies the first overflowing value on both explicit TaxID parsers.
	 * @param baseline True to expect historical int wrapping */
	private static void checkExplicitTaxidOverflow(final boolean baseline){
		final int expected=baseline ? Integer.MIN_VALUE : -1;
		final int stringOverflow=GiToTaxidInt.parseTaxidNumber(
				"tid|2147483648|", '|');
		final int byteOverflow=GiToTaxidInt.getID(toBytes("tid|2147483648|"));
		check(stringOverflow==expected,
				"String one-past-max explicit TaxID parsed to "+
				stringOverflow+" instead of "+expected);
		check(byteOverflow==expected,
				"byte[] one-past-max explicit TaxID parsed to "+
				byteOverflow+" instead of "+expected);
	}

	/** Verifies same-path initialize behavior after a failed load.
	 * @param root Temporary fixture root
	 * @param baseline True to expect the poisoned same-path retry
	 * @throws Exception If file creation fails */
	private static void checkFailedInitializeRetry(final Path root,
			final boolean baseline) throws Exception{
		GiToTaxidInt.unload();
		final Path table=root.resolve("late.tsv");
		try{
			GiToTaxidInt.initialize(table.toString());
			throw new AssertionError("Missing table should fail to initialize");
		}catch(RuntimeException expected){
			final String message=expected.getMessage();
			final String expectedMessage="Invalid file: "+table;
			check(expectedMessage.equals(message),
					"Missing table failure was '"+message+
					"' instead of '"+expectedMessage+"'");
		}
		Files.write(table, "A\tA.1\t7\t3\n".getBytes(StandardCharsets.US_ASCII));
		GiToTaxidInt.initialize(table.toString());
		if(baseline){
			check(!GiToTaxidInt.isInitialized(),
					"Baseline must expose same-path initialize retry poisoning");
			check(ARRAY_FIELD.get(null)==null,
					"Baseline retry should leave the dense table unloaded");
		}else{
			check(GiToTaxidInt.isInitialized(),
					"Fixed same-path initialize retry should load the table");
			check(GiToTaxidInt.getID(3)==7,
					"Fixed same-path retry table did not map GI 3");
		}
	}

	/** Verifies raw table GI overflow behavior through initialize().
	 * @param root Temporary fixture root
	 * @param baseline True to expect int-wrapped GI insertion
	 * @throws Exception If file creation or array reflection fails */
	private static void checkRawLoaderOverflow(final Path root,
			final boolean baseline) throws Exception{
		checkRawLoaderRow(root, "gi3.tsv", "3", 7,
				"normal GI 3");
		checkRawLoaderRow(root, "leading-zero-gi3.tsv",
				"0000000000000000000000000000003", 7,
				"leading-zero GI 3");
		checkOverflowRawLoaderRow(root, "int-overflow.tsv",
				"4294967299", baseline);
		checkOverflowRawLoaderRow(root, "int-overflow-leading-zero.tsv",
				"00000000004294967299", baseline);
		checkOverflowRawLoaderRow(root, "long-overflow.tsv",
				"18446744073709551619", baseline);
	}

	/** Verifies one accepted raw GI row through initialize().
	 * @param root Temporary fixture root
	 * @param name Fixture file name
	 * @param gi Raw GI token to write
	 * @param tid Expected TaxID
	 * @param label Case label for diagnostics
	 * @throws Exception If file creation or array reflection fails */
	private static void checkRawLoaderRow(final Path root, final String name,
			final String gi, final int tid, final String label) throws Exception{
		GiToTaxidInt.unload();
		final Path table=root.resolve(name);
		Files.write(table,
				("A\tA.1\t"+tid+"\t"+gi+"\n").getBytes(StandardCharsets.US_ASCII));
		GiToTaxidInt.initialize(table.toString());
		final int[] array=(int[])ARRAY_FIELD.get(null);
		check(array.length>3 && array[3]==tid,
				"Raw loader should store "+label+" at slot 3");
	}

	/** Verifies one out-of-range raw GI row through initialize().
	 * @param root Temporary fixture root
	 * @param name Fixture file name
	 * @param gi Raw GI token to write
	 * @param baseline True to expect int-wrapped GI insertion
	 * @throws Exception If file creation or array reflection fails */
	private static void checkOverflowRawLoaderRow(final Path root,
			final String name, final String gi, final boolean baseline)
			throws Exception{
		GiToTaxidInt.unload();
		final Path table=root.resolve(name);
		Files.write(table,
				("A\tA.1\t7\t"+gi+"\n").getBytes(StandardCharsets.US_ASCII));
		GiToTaxidInt.initialize(table.toString());
		final int[] array=(int[])ARRAY_FIELD.get(null);
		if(baseline){
			check(array.length>3 && array[3]==7,
					"Baseline raw loader should wrap GI "+gi+" to slot 3");
		}else{
			check(array.length==0,
					"Fixed raw loader should reject GI "+gi+"; length="+
					array.length);
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

	/** Reflected GiToTaxidInt.array field. */
	private static final Field ARRAY_FIELD;

	static{
		try{
			ARRAY_FIELD=GiToTaxidInt.class.getDeclaredField("array");
			ARRAY_FIELD.setAccessible(true);
		}catch(Exception e){throw new RuntimeException(e);}
	}
}
