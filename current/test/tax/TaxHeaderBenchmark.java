package test.tax;

import parse.Parse;
import tax.TaxTree;

/** Allocation-free timed header-ID parsing over a deterministic mixed-format corpus.
 * Strings and literal expected IDs are constructed before warmup and measurement.
 * @author Brian Bushnell, Shinobu */
public final class TaxHeaderBenchmark{

	/** Measures repeated parsing and checks the result against construction-time IDs.
	 * @param args Optional batch=65536 and rounds=64 */
	public static void main(String[] args){
		int batch=65536, rounds=64;
		for(String arg : args){
			String[] split=arg.split("=", 2);
			if(split.length!=2){throw new IllegalArgumentException("Expected key=value: "+arg);}
			if(split[0].equals("batch")){batch=Parse.parseIntKMG(split[1]);}
			else if(split[0].equals("rounds")){rounds=Parse.parseIntKMG(split[1]);}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		if(batch<1 || rounds<1){throw new IllegalArgumentException("batch and rounds must be positive");}
		final String[] headers=new String[batch];
		long expected=0;
		for(int i=0; i<headers.length; i++){
			final int id=1+(int)(((long)i*31)%2000000);
			switch(i&3){
				case 0: headers[i]="tid|"+id+"|sequence"; break;
				case 1: headers[i]="genome_ncbi_"+id+".fna"; break;
				case 2: headers[i]="plastid_name_tid_"+id; break;
				default: headers[i]="prefix ncbi|"+id+"|sequence"; break;
			}
			final int actual=TaxTree.parseTaxID(headers[i]);
			if(actual!=id){throw new AssertionError("Header oracle mismatch: "+headers[i]+" -> "+actual+", expected "+id);}
			expected+=id;
		}
		final long warm=run(headers, 8);
		if(warm!=expected*8){throw new AssertionError("Warmup checksum mismatch");}
		final long start=System.nanoTime();
		final long sum=run(headers, rounds);
		final long elapsed=System.nanoTime()-start;
		if(sum!=expected*rounds){throw new AssertionError("Timed checksum mismatch");}
		System.out.println("checksum="+sum+" expected="+(expected*rounds)+" batch="+batch+
				" rounds="+rounds+" compute_ns="+elapsed);
	}

	/** Prevents construction of this standalone benchmark. */
	private TaxHeaderBenchmark(){}

	/** Parses prebuilt strings without allocating benchmark objects in the timed loop.
	 * @param headers Deterministic corpus
	 * @param rounds Number of passes
	 * @return Sum of all parsed IDs */
	private static long run(String[] headers, int rounds){
		long sum=0;
		for(int round=0; round<rounds; round++){
			for(String header : headers){sum+=TaxTree.parseTaxID(header);}
		}
		return sum;
	}
}
