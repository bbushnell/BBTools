package test.tax;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import parse.Parse;
import tax.AnalyzeAccession;

/** Repeated AnalyzeAccession digitization workload over prebuilt valid accessions.
 * @author Jean */
public final class AnalyzeAccessionBenchmark{

	/** Measures String and accession2taxid byte-row digitization.
	 * @param args Optional batch=65536, rounds=64, warmup=auto
	 * @throws Exception If setup or shared-state reset fails */
	public static void main(String[] args) throws Exception{
		int batch=65536, rounds=64, warmup=-1;
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			if(split.length!=2){throw new IllegalArgumentException("Expected key=value: "+arg);}
			final String key=split[0].toLowerCase(), value=split[1];
			if(key.equals("batch")){batch=Parse.parseIntKMG(value);}
			else if(key.equals("rounds")){rounds=Parse.parseIntKMG(value);}
			else if(key.equals("warmup")){
				warmup=value.equalsIgnoreCase("auto") ? -1 : Parse.parseIntKMG(value);
			}else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		if(batch<1 || rounds<1){throw new IllegalArgumentException("batch and rounds must be positive");}
		if(warmup<0){warmup=warmupRounds(batch);}
		if(warmup<1){throw new IllegalArgumentException("warmup must be positive");}

		final Path root=Files.createTempDirectory("analyze-accession-benchmark-");
		try{
			resetAnalyzeAccessionPatterns();
			final Path patterns=root.resolve("patterns.tsv");
			Files.write(patterns, ("#Pattern\tCount\tCombos\tBits\n"+
					"LLL\t1\t17576\t14.10\n").getBytes(StandardCharsets.US_ASCII));
			AnalyzeAccession.loadCodeMap(patterns.toString());
			checkControls();

			final String[] roots=new String[batch];
			final String[] strings=new String[batch];
			final byte[][] rows=new byte[batch][];
			buildCorpus(roots, strings, rows);
			checkCorpus(roots, strings, rows);

			final long stringExpected=sumStrings(strings);
			final long rowExpected=sumRows(rows);
			runStrings(strings, warmup, stringExpected);
			runRows(rows, warmup, rowExpected);

			final long stringStart=System.nanoTime();
			final long stringChecksum=runStrings(strings, rounds, stringExpected);
			final long stringElapsed=System.nanoTime()-stringStart;
			final long rowStart=System.nanoTime();
			final long rowChecksum=runRows(rows, rounds, rowExpected);
			final long rowElapsed=System.nanoTime()-rowStart;
			System.out.println("string_checksum="+stringChecksum+
					" row_checksum="+rowChecksum+" string_expected="+
					(stringExpected*rounds)+" row_expected="+(rowExpected*rounds)+
					" batch="+batch+" rounds="+rounds+" warmup="+warmup+
					" string_compute_ns="+stringElapsed+
					" row_compute_ns="+rowElapsed);
		}finally{
			resetAnalyzeAccessionPatterns();
			deleteTempTree(root);
		}
	}

	/** Prevents construction of this standalone benchmark. */
	private AnalyzeAccessionBenchmark(){}

	/** Verifies that plain, dotted and colon suffixed accessions share one key. */
	private static void checkControls(){
		final long plain=AnalyzeAccession.digitize("ABC");
		final long dotted=AnalyzeAccession.digitize("ABC.1");
		final long colon=AnalyzeAccession.digitize("ABC:ALT");
		check(plain>0, "Plain control did not digitize: "+plain);
		check(dotted==plain, "Dotted control changed digitization: "+dotted);
		check(colon==plain, "Colon control changed digitization: "+colon);
	}

	/** Builds fixed roots, String keys and tab-delimited raw rows before timing.
	 * @param roots Destination for unsuffixed accession roots
	 * @param strings Destination for String accessions
	 * @param rows Destination for accession2taxid rows */
	private static void buildCorpus(final String[] roots, final String[] strings,
			final byte[][] rows){
		for(int i=0; i<strings.length; i++){
			final String root=root(i);
			roots[i]=root;
			final String accession;
			switch(i%3){
				case 0: accession=root; break;
				case 1: accession=root+".1"; break;
				default: accession=root+":ALT"; break;
			}
			strings[i]=accession;
			rows[i]=(accession+"\t"+root+".1\t"+(1000+i)+"\t0").getBytes(StandardCharsets.US_ASCII);
		}
	}

	/** Checks every suffixed corpus entry against its unsuffixed root.
	 * @param roots Unsuffixed accession roots
	 * @param strings String accession corpus
	 * @param rows accession2taxid byte-row corpus */
	private static void checkCorpus(final String[] roots, final String[] strings,
			final byte[][] rows){
		for(int i=0; i<roots.length; i++){
			final long expected=AnalyzeAccession.digitize(roots[i]);
			final long stringValue=AnalyzeAccession.digitize(strings[i]);
			final long rowValue=AnalyzeAccession.digitize(rows[i]);
			if(expected<0){throw new AssertionError("Invalid root: "+roots[i]);}
			check(stringValue==expected, strings[i]+" digitized to "+
					stringValue+" instead of root value "+expected);
			check(rowValue==expected, new String(rows[i], StandardCharsets.US_ASCII)+
					" digitized to "+rowValue+" instead of root value "+expected);
		}
	}

	/** Creates a three-letter key for this corpus index.
	 * @param i Corpus index
	 * @return Uppercase accession root */
	private static String root(final int i){
		final char a=(char)('A'+((i/676)%26));
		final char b=(char)('A'+((i/26)%26));
		final char c=(char)('A'+(i%26));
		return new String(new char[] {a, b, c});
	}

	/** Runs repeated String digitization and checks each round.
	 * @param strings Prebuilt accession corpus
	 * @param rounds Number of measured passes
	 * @param expected Expected checksum per pass
	 * @return Sum of all digitized keys */
	private static long runStrings(final String[] strings, final int rounds,
			final long expected){
		long sum=0;
		for(int round=0; round<rounds; round++){
			final long roundSum=sumStrings(strings);
			if(roundSum!=expected){
				throw new AssertionError("String checksum mismatch: "+roundSum+
						" != "+expected);
			}
			sum+=roundSum;
		}
		return sum;
	}

	/** Runs repeated byte-row digitization and checks each round.
	 * @param rows Prebuilt accession2taxid rows
	 * @param rounds Number of measured passes
	 * @param expected Expected checksum per pass
	 * @return Sum of all digitized keys */
	private static long runRows(final byte[][] rows, final int rounds,
			final long expected){
		long sum=0;
		for(int round=0; round<rounds; round++){
			final long roundSum=sumRows(rows);
			if(roundSum!=expected){
				throw new AssertionError("byte[] row checksum mismatch: "+roundSum+
						" != "+expected);
			}
			sum+=roundSum;
		}
		return sum;
	}

	/** Sums one pass over String accession digitization.
	 * @param strings Prebuilt accession corpus
	 * @return Sum of all digitized keys */
	private static long sumStrings(final String[] strings){
		long sum=0;
		for(String s : strings){
			final long value=AnalyzeAccession.digitize(s);
			if(value<0){throw new AssertionError("Invalid String digitization: "+s);}
			sum+=value;
		}
		return sum;
	}

	/** Sums one pass over byte-row digitization.
	 * @param rows Prebuilt accession2taxid rows
	 * @return Sum of all digitized keys */
	private static long sumRows(final byte[][] rows){
		long sum=0;
		for(byte[] row : rows){
			final long value=AnalyzeAccession.digitize(row);
			if(value<0){throw new AssertionError("Invalid row digitization: "+
					new String(row, StandardCharsets.US_ASCII));}
			sum+=value;
		}
		return sum;
	}

	/** Selects enough untimed calls to compile both hot digitization paths.
	 * @param batch Number of accessions per round
	 * @return Untimed warmup round count */
	private static int warmupRounds(final int batch){
		return Math.max(MIN_WARMUP_ROUNDS,
				(MIN_WARMUP_LOOKUPS+batch-1)/batch);
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

	/** Minimum untimed rounds before measurement. */
	private static final int MIN_WARMUP_ROUNDS=8;
	/** Minimum untimed String and byte-row digitizations before measurement. */
	private static final int MIN_WARMUP_LOOKUPS=50000;
}
