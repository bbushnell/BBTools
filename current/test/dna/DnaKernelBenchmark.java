package test.dna;

import java.nio.charset.StandardCharsets;

import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import shared.Shared;

/** Single-thread timing driver for unchanged, fully defined DNA/AA inputs.
 * Correctness edge cases live in the separate contract tests. A timing is useful
 * only after calibration reaches ten seconds of measured kernel work.
 * @author Shinobu
 */
public final class DnaKernelBenchmark{

	/** Requires case=amino|reverse|extract|scaffold and positive iterations=. */
	public static void main(String[] args){
		String mode=null;
		long iterations=0;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("iterations=")){iterations=Long.parseLong(arg.substring(11));}
			else{throw new IllegalArgumentException(arg);}
		}
		if(iterations<1 || !("amino".equals(mode) || "reverse".equals(mode) ||
				"extract".equals(mode) || "scaffold".equals(mode))){
			throw new IllegalArgumentException("Require a supported case= and positive iterations=");
		}
		Shared.setThreads(1);
		final int[][] oldLocs=Data.scaffoldLocs, oldLengths=Data.scaffoldLengths;
		final int oldPadding=Data.interScaffoldPadding;
		try{
			final DnaKernelBenchmark benchmark=new DnaKernelBenchmark();
			benchmark.run(mode, Math.min(500000L, Math.max(10000L, iterations/20)));
			final long start=System.nanoTime();
			final long checksum=benchmark.run(mode, iterations);
			final long nanos=System.nanoTime()-start;
			System.out.println("case\titerations\tkernel_nanos\tchecksum");
			System.out.println(mode+"\t"+iterations+"\t"+nanos+"\t"+checksum);
		}finally{
			Data.scaffoldLocs=oldLocs;
			Data.scaffoldLengths=oldLengths;
			Data.interScaffoldPadding=oldPadding;
		}
	}

	/** Builds fixed synthetic inputs without depending on the helpers being timed. */
	private DnaKernelBenchmark(){
		final byte[] residues="ACDEFGHIKLMNPQRSTVWY*".getBytes(StandardCharsets.US_ASCII);
		final byte[] bases={'A', 'C', 'G', 'T'}, complements={'T', 'G', 'C', 'A'};
		long state=19;
		for(int r=0; r<amino.length; r++){
			for(int i=0; i<amino[r].length; i++){
				amino[r][i]=residues[(i+r)%residues.length];
				state=state*6364136223846793005L+1442695040888963407L;
				final int code=(int)(state>>>62);
				dna[r][i]=bases[code];
				reverse[r][dna[r].length-i-1]=complements[code];
			}
		}
		final byte[] reference=new byte[4096];
		for(int i=0; i<reference.length; i++){reference[i]=bases[i&3];}
		chromosome.set(0, reference);
		Data.scaffoldLocs=new int[][]{new int[1024]};
		Data.scaffoldLengths=new int[][]{new int[1024]};
		Data.interScaffoldPadding=128;
		for(int i=0; i<1024; i++){
			Data.scaffoldLocs[0][i]=i*1152;
			Data.scaffoldLengths[0][i]=1024;
		}
	}

	/** Runs one kernel, retaining allocated results and a value-dependent checksum. */
	private long run(String mode, long iterations){
		long checksum=0;
		if(mode.equals("amino")){
			for(long i=0; i<iterations; i++){
				final byte[] result=AminoAcid.toNTs(amino[(int)i&63]);
				sink=result;
				checksum+=result[(int)(i%result.length)];
			}
		}else if(mode.equals("reverse")){
			for(long i=0; i<iterations; i++){
				final int index=(int)i&63;
				checksum+=AminoAcid.equalsReverseComp(dna[index], reverse[index]) ? 1 : 0;
			}
			if(checksum!=iterations){throw new AssertionError("Literal reverse complements must all match");}
		}else if(mode.equals("extract")){
			for(long i=0; i<iterations; i++){
				final int start=(int)i&2047;
				final byte[] result=chromosome.getBytes(start, start+255);
				sink=result;
				checksum+=result[(int)i&255];
			}
		}else{
			for(long i=0; i<iterations; i++){
				final int start=((int)i&1023)*1152+10;
				checksum+=Data.scaffoldOverlapLength(0, start, start+63);
			}
			if(checksum!=iterations*64){throw new AssertionError("Interior queries must retain all 64 bases");}
		}
		return checksum;
	}

	/** Forces each allocated output to escape; the same store is present in both builds. */
	private static volatile byte[] sink;
	private final byte[][] amino=new byte[64][1024], dna=new byte[64][1024], reverse=new byte[64][1024];
	private final ChromosomeArray chromosome=new ChromosomeArray(1, Shared.PLUS);
}
