package test.tax;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import tax.GiToTaxid;

/** Finite oracle for String GI headers with a leading greater-than character.
 * @author Jean */
public final class GiToTaxidHeaderTest{

	/** Runs baseline or fixed expectations for String and byte[] GI header lookup.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If fixture construction or cleanup fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		final Path root=Files.createTempDirectory("gi-header-");
		try{
			final Path table=root.resolve("gi.tsv");
			Files.write(table, (
					"A\tA.1\t999\t3\n"+
					"B\tB.1\t111\t7\n"+
					"C\tC.1\t333\t111\n").getBytes(StandardCharsets.US_ASCII));

			GiToTaxid.DISK_MODE=false;
			GiToTaxid.unload();
			GiToTaxid.initialize(table.toString());

			check(GiToTaxid.getID("gi|3|")==999,
					"String GI control failed for gi 3");
			check(GiToTaxid.getID("gi|7|")==111,
					"String GI control failed for gi 7");
			check(GiToTaxid.getID(toBytes(">gi|7|"))==111,
					"byte[] GI header should map gi 7 once");
			check(GiToTaxid.parseGiToTaxid(toBytes(">gi|7|"))==111,
					"byte[] parseGiToTaxid should map gi 7 once");

			final int stringHeader=GiToTaxid.getID(">gi|7|");
			final int stringParse=GiToTaxid.parseGiToTaxid(">gi|7|");
			final int expected=baseline ? 333 : 111;
			check(stringHeader==expected, "String GI header mapped to "+
					stringHeader+" instead of "+expected);
			check(stringParse==expected, "String parseGiToTaxid mapped to "+
					stringParse+" instead of "+expected);
		}finally{
			GiToTaxid.DISK_MODE=false;
			GiToTaxid.unload();
			deleteTempTree(root);
		}
		System.out.println("GI_TO_TAXID_HEADER_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private GiToTaxidHeaderTest(){}

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
