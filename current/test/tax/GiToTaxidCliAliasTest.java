package test.tax;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import fileIO.ReadWrite;
import tax.GiToTaxid;

/** Finite oracle for GiToTaxid command-line input/output aliasing.
 * @author Jean, Shinobu */
public final class GiToTaxidCliAliasTest{

	/** Runs the positional second-input/output alias baseline or fixed oracle.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If fixture file setup or cleanup fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		final Path root=Files.createTempDirectory("gi-cli-alias-");
		final boolean pigz=ReadWrite.USE_PIGZ;
		final boolean unpigz=ReadWrite.USE_UNPIGZ;
		final int zipLevel=ReadWrite.ZIPLEVEL;
		final int pigzBlockSize=ReadWrite.PIGZ_BLOCKSIZE;
		try{
			checkHealthyAbsoluteOutput(root);
			if(!baseline){
				checkHealthyCommaList(root);
				checkHealthyLiteralComma(root);
				checkHealthyCompressionFlag(root);
				checkDirectAliasRejected(root);
				checkCanonicalAliasRejected(root);
				checkHardLinkAliasRejected(root);
				checkExtraOutputRejected(root);
			}
			checkSecondInputAlias(root, baseline);
		}finally{
			GiToTaxid.DISK_MODE=false;
			GiToTaxid.unload();
			ReadWrite.USE_PIGZ=pigz;
			ReadWrite.USE_UNPIGZ=unpigz;
			ReadWrite.ZIPLEVEL=zipLevel;
			ReadWrite.PIGZ_BLOCKSIZE=pigzBlockSize;
			deleteTempTree(root);
		}
		System.out.println("GI_TO_TAXID_CLI_ALIAS_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private GiToTaxidCliAliasTest(){}

	/** Verifies the ordinary input plus absolute-output CLI path.
	 * @param root Temporary fixture root
	 * @throws Exception If fixture file setup or table reading fails */
	private static void checkHealthyAbsoluteOutput(final Path root)
			throws Exception{
		final byte[] input=
				"A\tA.1\t7\t3\n".getBytes(StandardCharsets.US_ASCII);
		final Path table=root.resolve("healthy.accession2taxid");
		final Path out=root.resolve("healthy").resolve("gitable.int2d");
		Files.createDirectories(out.getParent());
		Files.write(table, input);

		run(table, out);

		check(matches(table, input),
				"GiToTaxid CLI should preserve ordinary input bytes");
		checkTable(out, 3, 7, "ordinary absolute-output table");
	}

	/** Verifies the documented comma-list input syntax.
	 * @param root Temporary fixture root
	 * @throws Exception If fixture file setup or table reading fails */
	private static void checkHealthyCommaList(final Path root) throws Exception{
		final byte[] first=
				"C\tC.1\t13\t7\n".getBytes(StandardCharsets.US_ASCII);
		final byte[] second=
				"D\tD.1\t17\t11\n".getBytes(StandardCharsets.US_ASCII);
		final Path input1=root.resolve("comma1.accession2taxid");
		final Path input2=root.resolve("comma2.accession2taxid");
		final Path out=root.resolve("comma").resolve("gitable.int2d");
		Files.createDirectories(out.getParent());
		Files.write(input1, first);
		Files.write(input2, second);

		run(new String[] {input1+","+input2, out.toString()});

		check(matches(input1, first),
				"GiToTaxid CLI should preserve comma-list input 1");
		check(matches(input2, second),
				"GiToTaxid CLI should preserve comma-list input 2");
		checkTable(out, 7, 13, "comma-list requested output");
		checkTable(out, 11, 17, "comma-list requested output");
	}

	/** Verifies an existing comma-containing name stays one input path.
	 * @param root Temporary fixture root
	 * @throws Exception If fixture file setup or table reading fails */
	private static void checkHealthyLiteralComma(final Path root)
			throws Exception{
		final byte[] input=
				"H\tH.1\t31\t23\n".getBytes(StandardCharsets.US_ASCII);
		final Path table=root.resolve("literal,comma.tsv");
		final Path out=root.resolve("literal-comma").resolve("gitable.int2d");
		Files.createDirectories(out.getParent());
		Files.write(table, input);

		run(table, out);

		check(matches(table, input),
				"GiToTaxid CLI should preserve literal comma input bytes");
		checkTable(out, 23, 31, "literal comma requested output");
	}

	/** Verifies compression flags do not consume the positional output.
	 * @param root Temporary fixture root
	 * @throws Exception If fixture file setup or table reading fails */
	private static void checkHealthyCompressionFlag(final Path root)
			throws Exception{
		final byte[] input=
				"E\tE.1\t19\t13\n".getBytes(StandardCharsets.US_ASCII);
		final Path table=root.resolve("compressed.accession2taxid");
		final Path out=root.resolve("compressed").resolve("gitable.int2d");
		Files.createDirectories(out.getParent());
		Files.write(table, input);

		run(new String[] {table.toString(), "zl=1", out.toString()});

		check(matches(table, input),
				"GiToTaxid CLI should preserve compression-flag input");
		checkTable(out, 13, 19, "compression-flag requested output");
	}

	/** Verifies direct input/output alias rejection happens before writing.
	 * @param root Temporary fixture root
	 * @throws Exception If fixture file setup fails */
	private static void checkDirectAliasRejected(final Path root)
			throws Exception{
		final byte[] input=
				"F\tF.1\t23\t17\n".getBytes(StandardCharsets.US_ASCII);
		final Path table=root.resolve("direct-alias.tsv");
		Files.write(table, input);

		boolean thrown=false;
		try{
			run(new String[] {table.toString(), table.toString()});
		}catch(IllegalArgumentException e){
			thrown=e.getMessage().contains("Input and output files alias");
		}

		check(thrown, "Fixed GiToTaxid CLI should reject a direct alias");
		check(matches(table, input),
				"Fixed GiToTaxid CLI should reject direct aliases before writing");
	}

	/** Verifies canonical alias rejection before writing.
	 * @param root Temporary fixture root
	 * @throws Exception If fixture file setup fails */
	private static void checkCanonicalAliasRejected(final Path root)
			throws Exception{
		final byte[] input=
				"I\tI.1\t37\t29\n".getBytes(StandardCharsets.US_ASCII);
		final Path dir=root.resolve("canonical");
		Files.createDirectories(dir);
		final Path table=dir.resolve("alias.tsv");
		Files.write(table, input);

		boolean thrown=false;
		try{
			run(new String[] {
					dir.resolve(".").resolve("alias.tsv").toString(),
					table.toString()
			});
		}catch(IllegalArgumentException e){
			thrown=e.getMessage().contains("Input and output files alias");
		}

		check(thrown, "Fixed GiToTaxid CLI should reject a canonical alias");
		check(matches(table, input),
				"Fixed GiToTaxid CLI should reject canonical aliases before writing");
	}

	/** Verifies existing hard-link alias rejection where the filesystem permits it.
	 * @param root Temporary fixture root
	 * @throws Exception If fixture file setup fails */
	private static void checkHardLinkAliasRejected(final Path root)
			throws Exception{
		final byte[] input=
				"J\tJ.1\t41\t31\n".getBytes(StandardCharsets.US_ASCII);
		final Path dir=root.resolve("hard-link");
		Files.createDirectories(dir);
		final Path table=dir.resolve("source.tsv");
		final Path link=dir.resolve("link.tsv");
		Files.write(table, input);
		try{
			Files.createLink(link, table);
		}catch(UnsupportedOperationException e){
			System.out.println("SKIP hard-link alias: "+e);
			return;
		}

		boolean thrown=false;
		try{
			run(new String[] {link.toString(), table.toString()});
		}catch(IllegalArgumentException e){
			thrown=e.getMessage().contains("Input and output files alias");
		}

		check(thrown, "Fixed GiToTaxid CLI should reject a hard-link alias");
		check(matches(table, input),
				"Fixed GiToTaxid CLI should reject hard-link aliases before writing");
		System.out.println("hard_link_alias_checked=true");
	}

	/** Verifies extra bare output tokens are rejected as ambiguous.
	 * @param root Temporary fixture root
	 * @throws Exception If fixture file setup fails */
	private static void checkExtraOutputRejected(final Path root)
			throws Exception{
		final byte[] input=
				"G\tG.1\t29\t19\n".getBytes(StandardCharsets.US_ASCII);
		final Path table=root.resolve("extra-output.accession2taxid");
		final Path out1=root.resolve("extra1").resolve("gitable.int2d");
		final Path out2=root.resolve("extra2").resolve("gitable.int2d");
		Files.write(table, input);

		boolean thrown=false;
		try{
			run(new String[] {table.toString(), out1.toString(), out2.toString()});
		}catch(IllegalArgumentException e){
			thrown=e.getMessage().contains("Multiple output paths");
		}

		check(thrown,
				"Fixed GiToTaxid CLI should reject extra bare output tokens");
		check(matches(table, input),
				"Fixed GiToTaxid CLI should reject extra outputs before writing");
		check(!Files.exists(out1),
				"Fixed GiToTaxid CLI should not create the first extra output");
		check(!Files.exists(out2),
				"Fixed GiToTaxid CLI should not create the second extra output");
	}

	/** Verifies that a second accession argument is not also the output path.
	 * @param root Temporary fixture root
	 * @param baseline True to expect the second input to be overwritten
	 * @throws Exception If fixture file setup or table reading fails */
	private static void checkSecondInputAlias(final Path root,
			final boolean baseline) throws Exception{
		final byte[] first=
				"A\tA.1\t7\t3\n".getBytes(StandardCharsets.US_ASCII);
		final byte[] second=
				"B\tB.1\t11\t5\n".getBytes(StandardCharsets.US_ASCII);
		final Path input1=root.resolve("first.accession2taxid");
		final Path input2=root.resolve("second.accession2taxid");
		final Path outDir=root.resolve("out");
		final Path out=outDir.resolve("gitable.int2d");
		Files.createDirectories(outDir);
		Files.write(input1, first);
		Files.write(input2, second);

		run(input1, input2, out);

		check(matches(input1, first),
				"GiToTaxid CLI should preserve the first input");
		final boolean secondUnchanged=matches(input2, second);
		final boolean outExists=Files.exists(out);
		if(baseline){
			check(!secondUnchanged,
					"Baseline should overwrite the second accession input");
			check(!outExists,
					"Baseline should leave the requested output unwritten");
			checkTable(input2, 3, 7, "baseline alias output");
			checkTable(input2, 5, 11, "baseline alias output");
		}else{
			check(secondUnchanged,
					"Fixed GiToTaxid CLI should preserve the second input");
			check(outExists,
					"Fixed GiToTaxid CLI should create the requested output");
			checkTable(out, 3, 7, "fixed requested output");
			checkTable(out, 5, 11, "fixed requested output");
		}
	}

	/** Runs GiToTaxid.main with one input and one output.
	 * @param input Raw GI mapping input
	 * @param out Serialized table output */
	private static void run(final Path input, final Path out){
		GiToTaxid.DISK_MODE=false;
		GiToTaxid.unload();
		GiToTaxid.main(new String[] {input.toString(), out.toString()});
	}

	/** Runs GiToTaxid.main with two inputs and one output.
	 * @param input1 First raw GI mapping input
	 * @param input2 Second raw GI mapping input
	 * @param out Serialized table output */
	private static void run(final Path input1, final Path input2, final Path out){
		run(new String[] {input1.toString(), input2.toString(), out.toString()});
	}

	/** Runs GiToTaxid.main with explicit arguments.
	 * @param args Command-line arguments */
	private static void run(final String[] args){
		GiToTaxid.DISK_MODE=false;
		GiToTaxid.unload();
		GiToTaxid.main(args);
	}

	/** Checks one stored GI mapping in a serialized int[][] table.
	 * @param path Serialized table path
	 * @param gi GI index to check
	 * @param tid Expected TaxID
	 * @param label Case label for diagnostics */
	private static void checkTable(final Path path, final int gi, final int tid,
			final String label){
		final int[][] table=ReadWrite.read(int[][].class, path.toString(), true);
		check(table.length>0 && table[0]!=null && table[0].length>gi,
				label+" missing GI "+gi);
		check(table[0][gi]==tid,
				label+" stored GI "+gi+" -> "+table[0][gi]+
				" instead of "+tid);
	}

	/** Compares a file to expected bytes.
	 * @param path File to read
	 * @param expected Expected bytes
	 * @return True if the file exists and exactly matches {@code expected}
	 * @throws Exception If reading fails */
	private static boolean matches(final Path path, final byte[] expected)
			throws Exception{
		return Files.exists(path) &&
				java.util.Arrays.equals(Files.readAllBytes(path), expected);
	}

	/** Deletes an owned temporary tree.
	 * @param root Temporary fixture root
	 * @throws Exception If walking or deletion fails */
	private static void deleteTempTree(final Path root) throws Exception{
		if(!Files.exists(root)){return;}
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
