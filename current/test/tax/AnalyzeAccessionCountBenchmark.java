package test.tax;

import java.io.BufferedOutputStream;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;

import parse.Parse;
import tax.DiskAccessionTable;

/** Generates two-pattern accession input and validates counting or disk-table output.
 * @author Sucrose */
public final class AnalyzeAccessionCountBenchmark{

	/** Generates bounded input or checks counting summaries and disk-table output.
	 * @param args mode=generate out=&lt;input.tsv&gt; rows=&lt;even count&gt;, or
	 *        mode=check in=&lt;output.tsv&gt; rows=&lt;rows per file&gt; copies=&lt;input repetitions&gt;,
	 *        or mode=checkdisk in=&lt;table&gt;
	 * @throws Exception If file creation, output reading, or validation fails */
	public static void main(String[] args) throws Exception{
		String mode=null, in=null, out=null;
		int rows=-1, copies=1;
		for(String arg : args){
			String[] split=arg.split("=", 2);
			if(split.length!=2){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String key=split[0].toLowerCase();
			final String value=split[1];
			if(key.equals("mode")){mode=value;}
			else if(key.equals("in")){in=value;}
			else if(key.equals("out")){out=value;}
			else if(key.equals("rows")){rows=Parse.parseIntKMG(value);}
			else if(key.equals("copies")){copies=Parse.parseIntKMG(value);}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		if("generate".equals(mode)){
			generate(requiredPath(out, "out"), checkedRows(rows));
		}else if("check".equals(mode)){
			checkOutput(requiredPath(in, "in"), checkedRows(rows), checkedCopies(copies));
		}else if("checkdisk".equals(mode)){
			checkDisk(requiredPath(in, "in"));
		}else{
			throw new IllegalArgumentException("mode must be generate, check, or checkdisk; mode="+mode);
		}
	}

	/** Prevents construction of this standalone fixture. */
	private AnalyzeAccessionCountBenchmark(){}

	/** Writes {@code rows} alternating ABC and NZ_ accession rows.
	 * @param out Destination TSV path
	 * @param rows Positive even row count
	 * @throws Exception If the output cannot be written */
	private static void generate(final Path out, final int rows) throws Exception{
		final long bytes=INPUT_HEADER.length+(rows/2L)*PAIR_BYTES;
		if(bytes>MAX_BYTES){
			throw new IllegalArgumentException("Requested "+bytes+
					" bytes, above the 128 MiB input cap; rows="+rows+
					" maxRows="+MAX_ROWS);
		}
		final int chunkPairs=Math.min(rows/2, CHUNK_PAIRS);
		final byte[] chunk=makeChunk(chunkPairs);
		int pairsLeft=rows/2;
		try(OutputStream os=new BufferedOutputStream(Files.newOutputStream(out))){
			os.write(INPUT_HEADER);
			while(pairsLeft>=chunkPairs){
				os.write(chunk);
				pairsLeft-=chunkPairs;
			}
			if(pairsLeft>0){os.write(chunk, 0, pairsLeft*PAIR_BYTES);}
		}
		System.out.println("ANALYZE_ACCESSION_COUNT_BENCHMARK_GENERATE_PASS rows="+
				rows+" bytes="+bytes);
	}

	/** Builds one reusable buffer of alternating rows.
	 * @param pairs Number of ABC/NZ_ row pairs
	 * @return Reusable ASCII chunk */
	private static byte[] makeChunk(final int pairs){
		final byte[] chunk=new byte[pairs*PAIR_BYTES];
		for(int offset=0, i=0; i<pairs; i++, offset+=PAIR_BYTES){
			System.arraycopy(ABC_ROW, 0, chunk, offset, ABC_ROW.length);
			System.arraycopy(NZ_ROW, 0, chunk, offset+ABC_ROW.length, NZ_ROW.length);
		}
		return chunk;
	}

	/** Validates the exact two-pattern AnalyzeAccession summary.
	 * @param in AnalyzeAccession output file
	 * @param rows Rows per generated input file
	 * @param copies Number of comma-separated repetitions of that input
	 * @throws Exception If the output is missing or wrong */
	private static void checkOutput(final Path in, final int rows, final int copies)
			throws Exception{
		final long expected=Math.multiplyExact(rows/2L, copies);
		final List<String> lines=Files.readAllLines(in, StandardCharsets.US_ASCII);
		check(lines.size()==3, "Expected header plus exactly two pattern rows, found "+
				lines.size()+": "+lines);
		check(HEADER.equals(lines.get(0)), "Unexpected header: "+lines.get(0));
		final HashMap<String, String[]> map=new HashMap<String, String[]>();
		for(int i=1; i<lines.size(); i++){
			final String[] split=lines.get(i).split("\t", -1);
			check(split.length==4, "Expected four tab-delimited columns: "+lines.get(i));
			check(map.put(split[0], split)==null, "Duplicate pattern row: "+split[0]);
		}
		checkPattern(map, LLL_PATTERN, expected, LLL_COMBOS, LLL_BITS);
		checkPattern(map, NZ_PATTERN, expected, NZ_COMBOS, NZ_BITS);
		check(map.isEmpty(), "Unexpected pattern rows: "+map.keySet());
		System.out.println("ANALYZE_ACCESSION_COUNT_BENCHMARK_CHECK_PASS rows="+
				rows+" copies="+copies+" count="+expected);
	}

	/** Validates one pattern row.
	 * @param map Pattern row map
	 * @param pattern Expected pattern
	 * @param count Expected count
	 * @param combos Expected combination count string
	 * @param bits Expected log2 combination string */
	private static void checkPattern(final HashMap<String, String[]> map,
			final String pattern, final long count, final String combos,
			final String bits){
		final String[] row=map.remove(pattern);
		check(row!=null, "Missing pattern "+pattern);
		check(Long.toString(count).equals(row[1]), pattern+" count "+row[1]+" != "+count);
		check(combos.equals(row[2]), pattern+" combos "+row[2]+" != "+combos);
		check(bits.equals(row[3]), pattern+" bits "+row[3]+" != "+bits);
	}

	/** Validates the BuildDiskTables disk accession table from the generated input.
	 * @param in DiskAccessionTable file
	 * @throws Exception If the table is missing, corrupt, or contains wrong values */
	private static void checkDisk(final Path in) throws Exception{
		DiskAccessionTable table=null;
		try{
			table=new DiskAccessionTable(in.toString());
			check(table.get(ABC_KEY)==ABC_TAXID,
					"ABC key "+ABC_KEY+" expected "+ABC_TAXID+" found "+table.get(ABC_KEY));
			check(table.get(NZ_KEY)==NZ_TAXID,
					"NZ key "+NZ_KEY+" expected "+NZ_TAXID+" found "+table.get(NZ_KEY));
		}finally{
			if(table!=null){table.close();}
		}
		final long entryCount=entryCount(in);
		check(entryCount==DISK_ENTRY_COUNT,
				"entryCount "+entryCount+" != "+DISK_ENTRY_COUNT);
		final long bytes=Files.size(in);
		check(bytes==DISK_BYTES, "size "+bytes+" != "+DISK_BYTES);
		System.out.println("ANALYZE_ACCESSION_COUNT_BENCHMARK_CHECKDISK_PASS bytes="+
				bytes+" entryCount="+entryCount);
	}

	/** Reads the table header entry count independently from DiskAccessionTable.
	 * @param in DiskAccessionTable file
	 * @return Second big-endian long from the table header
	 * @throws Exception If the table header cannot be read */
	private static long entryCount(final Path in) throws Exception{
		try(RandomAccessFile raf=new RandomAccessFile(in.toFile(), "r")){
			raf.readLong();
			return raf.readLong();
		}
	}

	/** Converts a required filename argument to a path.
	 * @param name Filename from the command line
	 * @param flag Flag name for diagnostics
	 * @return Path for {@code name} */
	private static Path requiredPath(final String name, final String flag){
		if(name==null || name.length()<1){
			throw new IllegalArgumentException("Missing required "+flag+"=<path>");
		}
		return Paths.get(name);
	}

	/** Checks the requested row count.
	 * @param rows Requested accession row count
	 * @return {@code rows} when positive and even */
	private static int checkedRows(final int rows){
		if(rows<2 || (rows&1)!=0){
			throw new IllegalArgumentException("rows must be a positive even count: "+rows);
		}
		return rows;
	}

	/** Checks the requested input repetition count.
	 * @param copies Requested comma-separated input repetitions
	 * @return {@code copies} when positive */
	private static int checkedCopies(final int copies){
		if(copies<1){throw new IllegalArgumentException("copies must be positive: "+copies);}
		return copies;
	}

	/** Keeps fixture checks active under -da.
	 * @param condition Required condition
	 * @param message Failure message */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** AnalyzeAccession output header. */
	private static final String HEADER="#Pattern\tCount\tCombos\tBits";
	/** Letter-digit pattern from ABC rows. */
	private static final String LLL_PATTERN="LLLDDDDDDDDD";
	/** Special-case NZ pattern from NZ_ rows. */
	private static final String NZ_PATTERN="NZ-DDDDDDDDD";
	/** Combination count for {@link #LLL_PATTERN}. */
	private static final String LLL_COMBOS="17576000000000";
	/** log2 combination count for {@link #LLL_PATTERN}. */
	private static final String LLL_BITS="44.00";
	/** Combination count for {@link #NZ_PATTERN}. */
	private static final String NZ_COMBOS="1000000000";
	/** log2 combination count for {@link #NZ_PATTERN}. */
	private static final String NZ_BITS="29.90";
	/** Literal key for ABC000000000: {@code ((28*1000000000L+1)<<2)|0}. */
	private static final long ABC_KEY=112000000004L;
	/** Literal key for NZ_000000000: {@code ((0+1)<<2)|2}. */
	private static final long NZ_KEY=6L;
	/** TaxID for generated ABC rows. */
	private static final int ABC_TAXID=7;
	/** TaxID for generated NZ rows. */
	private static final int NZ_TAXID=8;
	/** Expected distinct table entries. */
	private static final long DISK_ENTRY_COUNT=2;
	/** Expected minimum table size with two entries and 16-slot partitions. */
	private static final long DISK_BYTES=3407904;
	/** accession2taxid header expected by AnalyzeAccession. */
	private static final byte[] INPUT_HEADER="accession\taccession.version\ttaxid\tgi\n".getBytes(StandardCharsets.US_ASCII);
	/** One fixed LLL pattern accession2taxid row. */
	private static final byte[] ABC_ROW="ABC000000000\tABC000000000.1\t7\t0\n".getBytes(StandardCharsets.US_ASCII);
	/** One fixed NZ pattern accession2taxid row. */
	private static final byte[] NZ_ROW="NZ_000000000\tNZ_000000000.1\t8\t0\n".getBytes(StandardCharsets.US_ASCII);
	/** Bytes per generated pair of rows. */
	private static final int PAIR_BYTES=ABC_ROW.length+NZ_ROW.length;
	/** Maximum generated input size. */
	private static final int MAX_BYTES=128*1024*1024;
	/** Maximum generated rows after rounding down to a complete row pair. */
	private static final int MAX_ROWS=2*((MAX_BYTES-INPUT_HEADER.length)/PAIR_BYTES);
	/** Row pairs buffered per write chunk. */
	private static final int CHUNK_PAIRS=8192;
}
