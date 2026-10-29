package test.tax;

import java.nio.file.Files;
import java.nio.file.Path;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import structures.ByteBuilder;
import tax.GiToTaxid;
import tax.GiToTaxidInt;

/** Measures raw GI table loading from a deterministic modest TSV corpus.
 * Input generation uses BBTools output code before warmup and measurement.
 * The reported compute_ns spans unload/load plus exact mapping validation.
 * @author Jean */
public final class GiRawLoaderBenchmark{

	/** Builds a small raw GI table and repeatedly unloads and reloads it.
	 * Reported compute_ns is timed after warmup and includes only unload/load
	 * plus exact mapping validation, not raw table generation.
	 * @param args mode=sparse|dense, keys=65536, copies=10, rounds=8, warmup=2
	 * @throws Exception If temporary fixture creation or deletion fails */
	public static void main(String[] args) throws Exception{
		String mode="sparse";
		int keys=65536, copies=10, rounds=8, warmup=2;
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			if(split.length!=2){throw new IllegalArgumentException("Expected key=value: "+arg);}
			final String key=split[0].toLowerCase(), value=split[1];
			if(key.equals("mode") || key.equals("kind")){mode=value.toLowerCase();}
			else if(key.equals("keys")){keys=Parse.parseIntKMG(value);}
			else if(key.equals("copies")){copies=Parse.parseIntKMG(value);}
			else if(key.equals("rounds")){rounds=Parse.parseIntKMG(value);}
			else if(key.equals("warmup")){warmup=Parse.parseIntKMG(value);}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		final int classMode=classMode(mode);
		validateSizes(keys, copies, rounds, warmup);

		final Path root=Files.createTempDirectory("gi-raw-loader-benchmark-");
		final Path raw=root.resolve("gi.tsv");
		try{
			writeRawTable(raw, keys, copies);
			final long expected=sumExpected(keys);
			run(raw, classMode, warmup, keys, expected);
			final long start=System.nanoTime();
			final long checksum=run(raw, classMode, rounds, keys, expected);
			final long elapsed=System.nanoTime()-start;
			final long expectedTotal=expected*rounds;
			if(checksum!=expectedTotal){
				throw new AssertionError("Timed checksum mismatch: "+checksum+
						" != "+expectedTotal);
			}
			System.out.println("mode="+mode+" checksum="+checksum+
					" expected="+expectedTotal+" keys="+keys+" copies="+
					copies+" rows="+((long)keys*copies)+" input_bytes="+
					Files.size(raw)+" rounds="+rounds+" warmup="+warmup+
					" compute_ns="+elapsed);
		}finally{
			GiToTaxid.unload();
			GiToTaxidInt.unload();
			Files.deleteIfExists(raw);
			Files.deleteIfExists(root);
		}
	}

	/** Prevents construction of this standalone benchmark. */
	private GiRawLoaderBenchmark(){}

	/** Converts the benchmark class selector into a dispatch code.
	 * @param mode Requested raw loader implementation
	 * @return Zero for sparse {@link GiToTaxid}, one for dense {@link GiToTaxidInt} */
	private static int classMode(final String mode){
		if(mode.equals("sparse")){return SPARSE;}
		if(mode.equals("dense")){return DENSE;}
		throw new IllegalArgumentException("Unknown mode: "+mode);
	}

	/** Rejects parameter values that would leave the intended bounded workload.
	 * @param keys Number of distinct small GIs
	 * @param copies Number of duplicate rows for each GI
	 * @param rounds Number of timed load repetitions
	 * @param warmup Number of untimed load repetitions */
	private static void validateSizes(final int keys, final int copies,
			final int rounds, final int warmup){
		if(keys<1 || keys>MAX_KEYS){
			throw new IllegalArgumentException("keys must be in [1,"+MAX_KEYS+"]: "+keys);
		}
		if(copies<1 || rounds<1 || warmup<1){
			throw new IllegalArgumentException("copies, rounds and warmup must be positive");
		}
		if(((long)keys)*copies>MAX_ROWS){
			throw new IllegalArgumentException("keys*copies must stay under "+MAX_ROWS+
					": "+(((long)keys)*copies));
		}
	}

	/** Writes the repeated raw GI rows outside the timed region.
	 * @param raw Destination TSV path under the owned temporary directory
	 * @param keys Number of deterministic GI keys
	 * @param copies Repetitions per key
	 * @throws Exception If the BBTools writer reports failure */
	private static void writeRawTable(final Path raw, final int keys,
			final int copies) throws Exception{
		final ByteStreamWriter out=new ByteStreamWriter(FileFormat.testOutput(
				raw.toString(), FileFormat.TEXT, null, true, true, false, false));
		out.start();
		Throwable primary=null;
		try{
			final ByteBuilder bb=new ByteBuilder(64);
			for(int copy=0; copy<copies; copy++){
				for(int gi=0; gi<keys; gi++){
					bb.clear().append("GI").append(gi).tab().append("GI").append(gi)
							.append(".1").tab().append(taxid(gi)).tab().append(gi).nl();
					out.print(bb);
				}
			}
		}catch(RuntimeException e){
			primary=e;
			throw e;
		}catch(Error e){
			primary=e;
			throw e;
		}finally{
			try{
				if(out.poisonAndWait()){
					final RuntimeException failure=
							new RuntimeException("Failed to write "+raw);
					if(primary==null){throw failure;}
					primary.addSuppressed(failure);
				}
			}catch(RuntimeException e){
				if(primary==null){throw e;}
				if(primary!=e){primary.addSuppressed(e);}
			}catch(Error e){
				if(primary==null){throw e;}
				if(primary!=e){primary.addSuppressed(e);}
			}
		}
	}

	/** Runs repeated raw loader reloads and validates every mapped key.
	 * @param raw Raw input table
	 * @param mode Sparse or dense dispatch code
	 * @param rounds Number of unload/load repetitions
	 * @param keys Number of exact GI keys to validate after each load
	 * @param expected Expected checksum for one complete mapping
	 * @return Sum of all validated TaxIDs */
	private static long run(final Path raw, final int mode, final int rounds,
			final int keys, final long expected){
		long sum=0;
		for(int round=0; round<rounds; round++){
			final long roundSum;
			if(mode==SPARSE){
				GiToTaxid.unload();
				GiToTaxid.initialize(raw.toString());
				roundSum=sumSparse(keys);
			}else{
				GiToTaxidInt.unload();
				GiToTaxidInt.initialize(raw.toString());
				roundSum=sumDense(keys);
			}
			if(roundSum!=expected){
				throw new AssertionError("Round "+round+" checksum mismatch: "+
						roundSum+" != "+expected);
			}
			sum+=roundSum;
		}
		return sum;
	}

	/** Checks the sparse two-dimensional table exactly.
	 * @param keys Number of GI keys
	 * @return Sum of all TaxIDs */
	private static long sumSparse(final int keys){
		long sum=0;
		for(int gi=0; gi<keys; gi++){
			final int expected=taxid(gi);
			final int actual=GiToTaxid.getID(gi);
			if(actual!=expected){
				throw new AssertionError("Sparse GI "+gi+" mapped to "+actual+
						" instead of "+expected);
			}
			sum+=actual;
		}
		return sum;
	}

	/** Checks the dense one-dimensional table exactly.
	 * @param keys Number of GI keys
	 * @return Sum of all TaxIDs */
	private static long sumDense(final int keys){
		long sum=0;
		for(int gi=0; gi<keys; gi++){
			final int expected=taxid(gi);
			final int actual=GiToTaxidInt.getID(gi);
			if(actual!=expected){
				throw new AssertionError("Dense GI "+gi+" mapped to "+actual+
						" instead of "+expected);
			}
			sum+=actual;
		}
		return sum;
	}

	/** Computes the expected checksum for a complete key range.
	 * @param keys Number of GI keys
	 * @return Sum of all deterministic TaxIDs */
	private static long sumExpected(final int keys){
		long sum=0;
		for(int gi=0; gi<keys; gi++){sum+=taxid(gi);}
		return sum;
	}

	/** Assigns a stable literal TaxID for one GI.
	 * @param gi Small GI key
	 * @return Deterministic TaxID */
	private static int taxid(final int gi){return 1000+((gi*31)%1000000);}

	/** Sparse {@link GiToTaxid} benchmark selector. */
	private static final int SPARSE=0;
	/** Dense {@link GiToTaxidInt} benchmark selector. */
	private static final int DENSE=1;
	/** Largest GI corpus accepted by this bounded benchmark. */
	private static final int MAX_KEYS=65536;
	/** Largest generated row count accepted by this bounded benchmark. */
	private static final long MAX_ROWS=1000000;
}
