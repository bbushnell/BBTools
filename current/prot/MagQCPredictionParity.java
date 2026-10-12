package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;

import structures.ByteBuilder;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import map.ObjectIntMap;

/** Compares native ProkCC rows to named CPU/GPU panel predictions for final capstone gates. */
public final class MagQCPredictionParity {

	public static void main(String[] args) throws Exception{
		final LinkedHashMap<String,String> o=new LinkedHashMap<String,String>();
		for(String arg:parse.Parser.parseConfig(args)){
			final int eq=arg.indexOf('=');
			require(eq>0, "Expected key=value: "+arg);
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT);
			require(KEYS.contains(key) && o.put(key, arg.substring(eq+1))==null, "Unknown or duplicate option: "+key);
		}
		if(parse.Parse.parseBoolean(o.getOrDefault("selftest", "f"))){selftest();return;}
		final int expected=Integer.parseInt(required(o, "expected"));
		final double valueTolerance=Double.parseDouble(required(o, "valuetolerance"));
		final double binaryTolerance=Double.parseDouble(required(o, "binarytolerance"));
		final boolean floatBits=parse.Parse.parseBoolean(o.getOrDefault("floatbits", "f"));
			require(expected>=0 && finite(valueTolerance) && finite(binaryTolerance) &&
			valueTolerance>=0 && binaryTolerance>=0, "Expected count and tolerances must be finite and nonnegative");
		require(!floatBits || (valueTolerance==0 && binaryTolerance==0), "Float-bit comparison requires zero tolerances");
		final LinkedHashMap<String,Double> heads=new LinkedHashMap<String,Double>();
		addHeads(heads, required(o, "valueheads"), valueTolerance);
		addHeads(heads, o.get("binaryheads"), binaryTolerance);
		require(!heads.isEmpty(), "At least one output head is required");
		final LinkedHashMap<String,double[]> nativeRows=load(required(o, "native"), heads);
		final LinkedHashMap<String,double[]> expectedRows=load(required(o, "expectedpredictions"), heads);
		require(nativeRows.size()==expected && expectedRows.size()==expected, "Prediction row count differs from expected");
		final double[] maximum=new double[heads.size()];
		String firstFailure=null;
		long mismatches=0;
		for(Map.Entry<String,double[]> entry:nativeRows.entrySet()){
			final double[] e=expectedRows.get(entry.getKey());
			require(e!=null, "Missing expected prediction for "+entry.getKey());
			final double[] n=entry.getValue();
			int index=0;
			for(double tolerance:heads.values()){
				final double diff=Math.abs(n[index]-e[index]);
				require(finite(diff), "Nonfinite difference at bin="+entry.getKey()+" head="+head(heads,index));
				if(!(floatBits ? sameFloatBits(n[index], e[index]) : diff<=tolerance)){
					mismatches++;
					if(firstFailure==null){
						firstFailure="Head differs beyond tolerance or float-bit identity: bin="+entry.getKey()+
							" head="+head(heads,index)+" native="+n[index]+" expected="+e[index]+" abs="+diff+" tolerance="+tolerance;
					}
				}
				maximum[index]=Math.max(maximum[index], diff); index++;
			}
		}
		require(expectedRows.keySet().equals(nativeRows.keySet()), "Prediction bin sets differ");
		final ByteBuilder out=new ByteBuilder();
		if(floatBits){out.append("#comparison\tfloat_raw_bits_after_decimal_roundtrip\n");}
		out.append("head\ttolerance\tmax_absolute_difference\n");
		int index=0;
		for(Map.Entry<String,Double> head:heads.entrySet()){
			out.append(head.getKey()).tab().appendSlow(head.getValue().doubleValue()).tab().appendSlow(maximum[index++]).nl();
		}
		if(o.containsKey("out")){writeText(o.get("out"), out);}
		else{System.out.print(out.toString());}
		// Preserve all per-head maxima on failure without relaxing the acceptance bounds.
		require(firstFailure==null, firstFailure+"; mismatched_values="+mismatches);
		System.out.println("MAGQC_PREDICTION_PARITY_PASS rows="+nativeRows.size()+" heads="+heads.size());
	}

	/** Loads named predictions and rejects collisions after filename normalization. */
	private static LinkedHashMap<String,double[]> load(String path, LinkedHashMap<String,Double> heads){
		final LinkedHashMap<String,double[]> result=new LinkedHashMap<String,double[]>();
		try(Rows rows=new Rows(path)){
			rows.col("bin_id");
			while(rows.next()){
				final String id=identity(rows.text("bin_id"));
				final double[] values=new double[heads.size()];
				int index=0;
				for(String head:heads.keySet()){values[index++]=number(rows, head);}
				require(result.put(id, values)==null, "Duplicate prediction: "+id);
			}
		}
		return result;
	}

	/** Resolves the supported report aliases and parses a finite value from its byte range. */
	private static double number(Rows rows, String head){
		int col=rows.columns.get(head);
		if(col<0 && head.endsWith("_error")){
			col=rows.columns.get("abs_residual_"+head.substring(0, head.length()-6));
		}
		if(col<0 && head.equals("contamination_ge_0_01")){
			col=rows.columns.get("contamination_ge_0.01");
		}
		if(col<0 && head.equals("completeness_gt_0_99")){
			col=rows.columns.get("completeness_gt_0.99");
		}
		if(col<0){rows.col(head);}
		rows.p.length(col);
		final double value=parse.Parse.parseDoubleSlow(rows.p.line(), rows.p.a(), rows.p.b());
		require(finite(value), "Nonfinite prediction field: "+head);
		return value;
	}

	/** Matches reports by basename, removing one gzip and one FASTA suffix. */
	private static String identity(final String path){
		final Path name=Paths.get(path).getFileName();
		require(name!=null, "Prediction path has no bin filename: "+path);
		String id=name.toString();
		if(id.endsWith(".gz")){id=id.substring(0, id.length()-3);}
		for(String suffix:FASTA_SUFFIXES){
			if(id.endsWith(suffix)){id=id.substring(0, id.length()-suffix.length()); break;}
		}
		require(!id.isEmpty(), "Empty bin identity: "+path);
		return id;
	}

	/** Writes a fresh summary and waits for I/O errors before reporting acceptance. */
	private static void writeText(final String path, final ByteBuilder contents){
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, true);
		writer.start();
		writer.print(contents);
		if(writer.poisonAndWait()){throw new IllegalStateException("Error writing "+path);}
	}

	/** Byte-oriented table reader supporting plain headers and native #columns headers. */
	private static final class Rows implements AutoCloseable {
		/** Reads exactly one header, skipping leading metadata and blank lines. */
		Rows(final String path){
			file=ByteFile.makeByteFile(path, false);
			int n=-1;
			try{
				for(byte[] line=file.nextLine(); line!=null; line=file.nextLine()){
					if(line.length==0){continue;}
					p.set(line);
					final boolean marked=p.termEquals("#columns", 0);
					if(line[0]=='#' && !marked){continue;}
					final int skip=marked ? 1 : 0;
					n=p.terms()-skip;
					for(int i=0; i<n; i++){
						final String name=p.parseString(i+skip);
						require(!name.isEmpty() && !columns.contains(name),
							"Empty or repeated table column in "+path+": "+name);
						columns.put(name, i);
					}
					break;
				}
				require(n>0, "Missing table header: "+path);
			}catch(RuntimeException|Error e){file.close(); throw e;}
			width=n;
		}

		/** Advances only over full-width data rows; malformed rows cannot be silently skipped. */
		boolean next(){
			final byte[] row=file.nextLine();
			if(row==null){return false;}
			p.set(row);
			require(p.terms()==width, "Ragged prediction row: expected "+width+" fields, got "+p.terms());
			return true;
		}

		/** Requires a named column so a missing field cannot be interpreted as another head. */
		int col(final String name){
			final int n=columns.get(name);
			require(n>=0, "Missing prediction column: "+name);
			return n;
		}

		/** Materializes only retained string fields such as bin identity. */
		String text(final String name){return p.parseString(col(name));}
		/** Propagates asynchronous reader failure on both successful and rejected input. */
		@Override
		public void close(){require(!file.close(), "Prediction table read failed");}

		final ByteFile file;
		final int width;
		final LineParser1 p=new LineParser1('\t');
		final ObjectIntMap<String> columns=new ObjectIntMap<String>(String.class);
	}

	/** Reconstructs the reported FP32 values and retains the distinction between signed zeros. */
	private static boolean sameFloatBits(final double left, final double right){
		final float a=(float)left, b=(float)right;
		require(Float.isFinite(a) && Float.isFinite(b), "Float-bit comparison requires finite representable FP32 values");
		return Float.floatToRawIntBits(a)==Float.floatToRawIntBits(b);
	}

	private static void addHeads(LinkedHashMap<String,Double> heads, String text, double tolerance){
		if(text==null || text.isEmpty()){return;}
		for(String head:text.split(",", -1)){
			require(!head.isEmpty() && heads.put(head, tolerance)==null, "Empty or duplicate head: "+head);
		}
	}

	private static String head(LinkedHashMap<String,Double> heads, int index){
		int i=0; for(String head:heads.keySet()){if(i++==index){return head;}}
		throw new AssertionError(index);
	}

	private static void selftest() throws Exception{
		final Path root=Files.createTempDirectory("magqc-prediction-parity-");
		final Path nativeFile=root.resolve("native.tsv"), gpu=root.resolve("gpu.tsv"), out=root.resolve("out.tsv");
		Files.write(nativeFile, ("#columns\tbin_id\tgene_completeness\tcontamination_ge_0.01\n"+
			"/tmp/a.fa\t0.5\t0.25001\n/tmp/b.fna.gz\t0.75\t0\n").getBytes(StandardCharsets.UTF_8));
		Files.write(gpu, ("bin_id\tgene_completeness\tcontamination_ge_0_01\n"+
			"a\t0.500001\t0.25002\nb\t0.749999\t0.00001\n").getBytes(StandardCharsets.UTF_8));
		main(new String[]{"native="+nativeFile, "expectedpredictions="+gpu, "expected=2",
			"valueheads=gene_completeness", "binaryheads=contamination_ge_0_01",
			"valuetolerance=0.000002", "binarytolerance=0.00002", "out="+out});
		require(new String(Files.readAllBytes(out), StandardCharsets.UTF_8).contains("gene_completeness"),
			"Selftest did not write a per-head summary");
		boolean rejected=false;
		try{main(new String[]{"native="+nativeFile, "expectedpredictions="+gpu, "expected=2",
			"valueheads=gene_completeness", "valuetolerance=0.0000001", "binarytolerance=0.00002"});}
		catch(IllegalArgumentException expected){rejected=true;}
		require(rejected, "Beyond-tolerance value head was accepted");
		final Path failureReport=root.resolve("failed.tsv");
		Files.write(nativeFile, "bin_id\tvalue\na\t0.2\nb\t0.9\n".getBytes(StandardCharsets.UTF_8));
		Files.write(gpu, "bin_id\tvalue\na\t0.1\nb\t0.1\n".getBytes(StandardCharsets.UTF_8));
		rejected=false;
		try{main(new String[]{"native="+nativeFile, "expectedpredictions="+gpu, "expected=2",
			"valueheads=value", "valuetolerance=0.05", "binarytolerance=0.05", "out="+failureReport});}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("mismatched_values=2");}
		require(rejected && new String(Files.readAllBytes(failureReport), StandardCharsets.UTF_8).contains("0.8"),
			"A rejected comparison must retain the later maximum, not stop after the first mismatch");
		require(sameFloatBits(0.1, 0.100000001), "Equivalent decimal forms must reconstruct the same reported FP32 value");
		require(!sameFloatBits(0.0, -0.0), "Float-bit comparison must reject opposite zero signs");
		require(!sameFloatBits(1.0f, Math.nextUp(1.0f)), "Float-bit comparison must reject adjacent FP32 values");
		Files.write(nativeFile, "bin_id\tvalue\na\t0.0\n".getBytes(StandardCharsets.UTF_8));
		Files.write(gpu, "bin_id\tvalue\na\t-0.0\n".getBytes(StandardCharsets.UTF_8));
		rejected=false;
		try{main(new String[]{"native="+nativeFile, "expectedpredictions="+gpu, "expected=1",
			"valueheads=value", "valuetolerance=0", "binarytolerance=0", "floatbits=t"});}
		catch(IllegalArgumentException expected){rejected=true;}
		require(rejected, "Reported opposite signed zeros passed the float-bit CLI");
		final LinkedHashMap<String,Double> aliases=new LinkedHashMap<String,Double>();
		addHeads(aliases, "gene_completeness_error,completeness_gt_0_99", 0);
		Files.write(nativeFile, ("#fixture\tparser\n\n#columns\tbin_id\tabs_residual_gene_completeness\tcompleteness_gt_0.99\n"+
			"/tmp/a.fasta.gz\t1e-3\t0.9\n").getBytes(StandardCharsets.UTF_8));
		final double[] aliasValues=load(nativeFile.toString(), aliases).get("a");
		require(aliasValues!=null && aliasValues[0]==0.001 && aliasValues[1]==0.9,
			"Metadata, FASTA identity, exponent parsing or supported head aliases changed");
		final LinkedHashMap<String,Double> oneHead=new LinkedHashMap<String,Double>();
		addHeads(oneHead, "value", 0);
		final String[] malformed={"#metadata\tonly\n", "bin_id\tvalue\tvalue\na\t1\t2\n",
			"bin_id\tvalue\na\n", "bin_id\tother\na\t1\n", "bin_id\tvalue\n.fa\t1\n",
			"bin_id\tvalue\n/tmp/a.fa\t1\na.fasta.gz\t2\n", "bin_id\tvalue\na\tNaN\n"};
		for(String table:malformed){
			Files.write(nativeFile, table.getBytes(StandardCharsets.UTF_8));
			rejected=false;
			try{load(nativeFile.toString(), oneHead);}
			catch(IllegalArgumentException expected){rejected=true;}
			require(rejected, "Malformed prediction table was accepted: "+table);
		}
		System.out.println("MAGQC_PREDICTION_PARITY_SELFTEST_PASS malformed="+malformed.length);
	}

	private static String required(LinkedHashMap<String,String> o, String key){
		final String value=o.get(key); require(value!=null && !value.isEmpty(), "Required: "+key); return value;
	}
	private static boolean finite(double value){return !Double.isNaN(value) && !Double.isInfinite(value);}
	private static void require(boolean good, String why){if(!good){throw new IllegalArgumentException(why);}}
	private static final String[] FASTA_SUFFIXES={".fasta", ".fna", ".fa"};
	private static final java.util.HashSet<String> KEYS=new java.util.HashSet<String>(Arrays.asList("native",
		"expectedpredictions", "expected", "valueheads", "binaryheads", "valuetolerance", "binarytolerance", "floatbits", "out", "selftest"));
}
