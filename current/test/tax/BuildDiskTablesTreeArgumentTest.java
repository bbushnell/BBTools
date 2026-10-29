package test.tax;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import tax.BuildDiskTables;
import tax.TaxTree;

/** Constructor-only BuildDiskTables tree and nullable file argument parsing checks.
 * @author Jean */
public final class BuildDiskTablesTreeArgumentTest{

	/** Runs baseline or fixed expectations for raw and config-expanded tree args.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If fixture setup or reflection fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		final Path root=Files.createTempDirectory("build-disk-tables-tree-");
		try{
			for(String alias : new String[]{"tree", "taxtree"}){
				checkAlias(root, alias, baseline, false);
				checkAlias(root, alias, baseline, true);
			}
			checkNullableField(root, ACCESSION_FILE, baseline, "accession");
			checkNullableField(root, GI_TABLE_FILE, baseline, "table", "gi",
					"gitable");
			checkNullableField(root, PATTERN_FILE, baseline, "pattern");
			checkNullableField(root, ACC_OUT, baseline, "accout", "accessionout");
			checkNullableField(root, GI_OUT, baseline, "giout");
		}finally{
			deleteTempTree(root);
		}
		System.out.println("BUILD_DISK_TABLES_TREE_ARGUMENT_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private BuildDiskTablesTreeArgumentTest(){}

	/** Checks one tree-key spelling through direct or config-expanded args.
	 * @param root Temporary directory for config fixtures
	 * @param alias Tree argument key
	 * @param baseline True to expect the original empty-arg behavior
	 * @param config True to expand the argument from a config file
	 * @throws Exception If constructor reflection or config setup fails */
	private static void checkAlias(final Path root, final String alias,
			final boolean baseline, final boolean config) throws Exception{
		final String custom=root.resolve(alias+".tsv").toString();
		final Path literalEquals=root.resolve(alias+"=literal.tsv");
		final String literal=literalEquals.toString();
		final int equals=literal.indexOf('=');
		check(equals>=0, "Literal equals fixture path lacks equals: "+literal);

		checkValue(root, alias, null, TaxTree.defaultTreeFile(), baseline, config,
				"bare");
		final String actualEmpty=checkValue(root, alias, "",
				baseline ? TaxTree.defaultTreeFile() : null, baseline, config,
				"empty");
		final String actualNull=checkValue(root, alias, "null",
				baseline ? TaxTree.defaultTreeFile() : null, baseline, config,
				"null");
		checkValue(root, alias, "none", null, baseline, config, "none");
		checkValue(root, alias, "f", null, baseline, config, "false");
		checkValue(root, alias, "0", null, baseline, config, "zero");
		checkValue(root, alias, "t", TaxTree.defaultTreeFile(), baseline, config,
				"true");
		checkValue(root, alias, "1", TaxTree.defaultTreeFile(), baseline, config,
				"one");
		checkValue(root, alias, "auto", TaxTree.defaultTreeFile(), baseline, config,
				"auto");
		checkValue(root, alias, custom, custom, baseline, config, "custom");
		final String actualLiteral=checkValue(root, alias, literal,
				baseline ? literal.substring(0, equals) : literal,
				baseline, config, "literal-equals");
		receipt(alias, config, actualEmpty, actualNull, actualLiteral);
	}

	/** Checks nullable file options through direct and config-expanded args.
	 * @param root Temporary directory for config fixtures
	 * @param field Reflected BuildDiskTables field to check
	 * @param baseline True to expect the original literal equals behavior
	 * @param aliases Equivalent argument spellings for the reflected field
	 * @throws Exception If constructor reflection or config setup fails */
	private static void checkNullableField(final Path root, final Field field,
			final boolean baseline, final String... aliases) throws Exception{
		for(String alias : aliases){
			checkNullableOption(root, alias, field, baseline, false);
			checkNullableOption(root, alias, field, baseline, true);
		}
	}

	/** Checks one nullable non-tree file option spelling.
	 * @param root Temporary directory for config fixtures
	 * @param alias Argument key
	 * @param field Reflected BuildDiskTables field to check
	 * @param baseline True to expect the original literal equals behavior
	 * @param config True to expand the argument from a config file
	 * @throws Exception If constructor reflection or config setup fails */
	private static void checkNullableOption(final Path root, final String alias,
			final Field field, final boolean baseline, final boolean config)
			throws Exception{
		final String ordinary=root.resolve(alias+"-file.tsv").toString();
		final Path literalEquals=root.resolve(alias+"=file.tsv");
		final String literal=literalEquals.toString();
		final int equals=literal.indexOf('=');
		check(equals>=0, "Literal equals fixture path lacks equals: "+literal);

		final String actualEmpty=checkFieldValue(root, alias, "", field, null,
				config, "empty");
		final String actualNull=checkFieldValue(root, alias, "null", field, null,
				config, "null");
		final String actualOrdinary=checkFieldValue(root, alias, ordinary, field,
				ordinary, config, "ordinary");
		final String actualLiteral=checkFieldValue(root, alias, literal, field,
				baseline ? literal.substring(0, equals) : literal, config,
				"literal-equals");
		nullableReceipt(alias, config, actualEmpty, actualNull, actualOrdinary,
				actualLiteral);
	}

	/** Checks one nullable non-tree file option value.
	 * @param root Temporary directory for config fixtures
	 * @param alias Argument key
	 * @param value Exact value after '='
	 * @param field Reflected BuildDiskTables field to check
	 * @param expected Expected reflected field value
	 * @param config True to expand the argument from a config file
	 * @param label Case label for diagnostics
	 * @return Reflected constructor value after validation
	 * @throws Exception If constructor reflection or config setup fails */
	private static String checkFieldValue(final Path root, final String alias,
			final String value, final Field field, final String expected,
			final boolean config, final String label) throws Exception{
		final String arg=alias+"="+value;
		final String actual=reflectedValue(root, arg, config, field);
		check(equals(expected, actual), label+" "+alias+" field="+field.getName()+
				" config="+config+" expected "+expected+" found "+actual);
		return actual;
	}

	/** Checks one value for one tree-key spelling.
	 * @param root Temporary directory for config fixtures
	 * @param alias Tree argument key
	 * @param value Null for a bare key, otherwise the exact value after '='
	 * @param expected Expected reflected taxTreeFile value
	 * @param baseline True for the baseline expectation label
	 * @param config True to expand the argument from a config file
	 * @param label Case label for diagnostics
	 * @return Reflected constructor value after validation
	 * @throws Exception If constructor reflection or config setup fails */
	private static String checkValue(final Path root, final String alias,
			final String value, final String expected, final boolean baseline,
			final boolean config, final String label) throws Exception{
		final String arg=(value==null ? alias : alias+"="+value);
		final String actual=taxTreeFile(root, arg, config);
		check(equals(expected, actual), label+" "+alias+" config="+config+
				" baseline="+baseline+" expected "+expected+" found "+actual);
		return actual;
	}

	/** Prints the changed-case values for the saved baseline/fixed logs.
	 * @param alias Tree argument key
	 * @param config True when expanded from a config file
	 * @param empty Actual {@code alias=} value
	 * @param explicitNull Actual {@code alias=null} value
	 * @param literal Actual value for a filename containing '=' */
	private static void receipt(final String alias, final boolean config,
			final String empty, final String explicitNull, final String literal){
		System.out.println("BUILD_DISK_TABLES_TREE_ARGUMENT_VALUES alias="+alias+
				" config="+config+" empty="+format(empty)+" explicit_null="+
				format(explicitNull)+" literal_equals="+format(literal));
	}

	/** Formats reflected nulls visibly for receipt lines.
	 * @param value Reflected taxTreeFile value
	 * @return Printable receipt value */
	private static String format(final String value){
		return value==null ? "<null>" : value;
	}

	/** Prints nullable file-option changed-case values for saved logs.
	 * @param alias Argument key
	 * @param config True when expanded from a config file
	 * @param empty Actual {@code alias=} value
	 * @param explicitNull Actual {@code alias=null} value
	 * @param ordinary Actual value for an ordinary filename
	 * @param literal Actual value for a filename containing '=' */
	private static void nullableReceipt(final String alias, final boolean config,
			final String empty, final String explicitNull, final String ordinary,
			final String literal){
		System.out.println("BUILD_DISK_TABLES_NULLABLE_FILE_VALUES alias="+alias+
				" config="+config+" empty="+format(empty)+" explicit_null="+
				format(explicitNull)+" ordinary="+format(ordinary)+
				" literal_equals="+format(literal));
	}

	/** Reflects taxTreeFile after BuildDiskTables construction.
	 * Intentionally does not call process(), so no tree or server is opened.
	 * @param root Temporary directory for config fixtures
	 * @param arg Tree argument to pass directly or through config=
	 * @param config True to expand {@code arg} from a config file
	 * @return Reflected constructor output
	 * @throws Exception If fixture setup or reflection fails */
	private static String taxTreeFile(final Path root, final String arg,
			final boolean config) throws Exception{
		return reflectedValue(root, arg, config, TAX_TREE_FILE);
	}

	/** Reflects one field after BuildDiskTables construction.
	 * Intentionally does not call process(), so no tree or server is opened.
	 * @param root Temporary directory for config fixtures
	 * @param arg Argument to pass directly or through config=
	 * @param config True to expand {@code arg} from a config file
	 * @param field Reflected BuildDiskTables field to read
	 * @return Reflected constructor output
	 * @throws Exception If fixture setup or reflection fails */
	private static String reflectedValue(final Path root, final String arg,
			final boolean config, final Field field) throws Exception{
		String[] args={arg, "outstream=null"};
		if(config){
			final Path file=Files.createTempFile(root, "tree-arg-", ".config");
			Files.write(file, (arg+"\n").getBytes(StandardCharsets.US_ASCII));
			args=new String[]{"config="+file, "outstream=null"};
		}
		final BuildDiskTables tables=new BuildDiskTables(args);
		return (String)field.get(tables);
	}

	/** Deletes the owned temporary tree, visiting children before parents.
	 * @param root Temporary fixture root
	 * @throws Exception If walking or deletion fails */
	private static void deleteTempTree(final Path root) throws Exception{
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths.sorted(
					java.util.Comparator.reverseOrder())::iterator){
				Files.deleteIfExists(path);
			}
		}
	}

	/** Null-safe equality.
	 * @param a First value
	 * @param b Second value
	 * @return True if both are null or equal */
	private static boolean equals(final Object a, final Object b){
		return a==null ? b==null : a.equals(b);
	}

	/** Keeps fixture checks active under -da.
	 * @param condition Required condition
	 * @param message Failure message */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Reflected tree input path. */
	private static final Field TAX_TREE_FILE=field("taxTreeFile");
	/** Reflected accession input path. */
	private static final Field ACCESSION_FILE=field("accessionFile");
	/** Reflected GI table input path. */
	private static final Field GI_TABLE_FILE=field("giTableFile");
	/** Reflected pattern input path. */
	private static final Field PATTERN_FILE=field("patternFile");
	/** Reflected accession output path. */
	private static final Field ACC_OUT=field("accOut");
	/** Reflected GI output path. */
	private static final Field GI_OUT=field("giOut");

	/** Returns an accessible BuildDiskTables field.
	 * @param name Field name
	 * @return Accessible field */
	private static Field field(final String name){
		try{
			final Field field=BuildDiskTables.class.getDeclaredField(name);
			field.setAccessible(true);
			return field;
		}catch(NoSuchFieldException e){throw new ExceptionInInitializerError(e);}
	}
}
