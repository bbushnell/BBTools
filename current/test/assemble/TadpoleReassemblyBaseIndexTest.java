package test.assemble;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.Random;

import assemble.Tadpole;
import assemble.Tadpole1;
import assemble.Tadpole2;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;
import structures.IntList;
import ukmer.Kmer;

/** Finite reassembly contracts where a wrong base's numeric code equals support depth.
 * Counts come from actual repeated reference records; no table or correction method is mocked.
 * Each invocation uses a fresh JVM and new output directory. This is not a speed benchmark.
 * @author Brian Bushnell, Shinobu
 */
public final class TadpoleReassemblyBaseIndexTest{

	/** Records all six outcomes before reporting any correction-contract failure. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, TadpoleReassemblyBaseIndexTest.class, false);
		try{
			int k=0, depth=0;
			Path root=null;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(pair.length==2, "Expected flag=value: "+arg);
				if(pair[0].equalsIgnoreCase("k")){k=Integer.parseInt(pair[1]);}
				else if(pair[0].equalsIgnoreCase("depth")){depth=Integer.parseInt(pair[1]);}
				else if(pair[0].equalsIgnoreCase("outdir")){root=Paths.get(pair[1]).toAbsolutePath();}
				else{throw new IllegalArgumentException("Unknown fixture option: "+arg);}
			}
			check(root!=null && (k==31 || k==63) && (depth==2 || depth==3),
					"Specify new outdir=DIR, k=31|63 and depth=2|3.");
			Files.createDirectory(root);
			final byte[] truth=new byte[4*k+21];
			final Random random=new Random(8192741L);
			for(int i=0; i<truth.length; i++){truth[i]=BASES[random.nextInt(4)];}
			final int position=2*k+5;
			truth[position]='A';
			final String sequence=new String(truth, StandardCharsets.US_ASCII);
			final StringBuilder fasta=new StringBuilder();
			for(int copy=0; copy<depth; copy++){fasta.append('>').append(copy).append('\n').append(sequence).append('\n');}
			final Path reference=root.resolve("counts.fa");
			Files.write(reference, fasta.toString().getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
			final String[] options={"in="+reference, "out=null", "mode=correct", "k="+k, "packed=f",
					"t=1", "buildthreads=1", "ways=1", "initialsize=1024", "prealloc=f", "prefilter=0",
					"minprob=0", "minprobmain=f", "rcomp=t", "maskmiddle=f", "maskcore=f", "qin=33",
					"mce=2", "shave=f", "rinse=f", "interleaved=f", "pigz=f", "unpigz=f", "showstats=f",
					"aggressive=f", "conservative=f", "errorpath=1", "errormult1=16", "errormult2=2.6",
					"errormultqfactor=0.002", "errorlowerconst=4", "mincountcorrect=3", "deadzone=0",
					"ecc=f", "ecco=f", "merge=f", "fixindels=f", "hpindel=f", "markerrors=f",
					"hollow=f", "hollowrescue=f"};
			Files.write(root.resolve("table.args"), Arrays.asList(options), StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
			final Tadpole tad=Tadpole.makeTadpole(options, false);
			try{
				check(tad.tables().kbig()==k && !Kmer.PACKED && !Kmer.MASK_CORE &&
						(k==31 ? tad instanceof Tadpole1 : tad instanceof Tadpole2), "Unexpected count-table geometry.");
				tad.loadKmers(new Timer());
				verifyCounts(tad, truth, k, depth, -1);
				for(int base=1; base<BASES.length; base++){
					final byte[] mutation=truth.clone();
					mutation[position]=BASES[base];
					verifyCounts(tad, mutation, k, depth, position);
				}
				final StringBuilder report=new StringBuilder("case\tquality\texpected\tcorrected\tbases\tcounts\tqualities\tpassed\n");
				int failures=0;
				for(boolean quality : new boolean[]{false, true}){
					failures+=runCase(tad, truth, k, depth, position, (byte)'A', "clean", quality, report);
					failures+=runCase(tad, truth, k, depth, position, (byte)'C', "control", quality, report);
					failures+=runCase(tad, truth, k, depth, position, BASES[depth], "index-equals-depth", quality, report);
				}
				Files.write(root.resolve("checks.tsv"), report.toString().getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
				System.out.print(report);
				System.out.println("REASSEMBLY_BASE_INDEX k="+k+" depth="+depth+" cases=6 failed="+failures);
				if(failures>0){throw new AssertionError("Reassembly base-index contract failures: "+failures);}
			}finally{tad.clearData();}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Uses coordinate-defined truth to check every loaded original or mutated k-mer window. */
	private static IntList verifyCounts(final Tadpole tad, final byte[] bases, final int k,
			final int depth, final int mutation){
		final IntList counts=new IntList(bases.length);
		tad.tables().fillCounts(bases, counts, new Kmer(k));
		check(counts.size==bases.length-k+1, "Fixture lacks a count for every window.");
		for(int start=0; start<counts.size; start++){
			final int expected=start<=mutation && mutation<start+k ? 0 : depth;
			check(counts.get(start)==expected, "Nonunique or unexpected fixture kmer at "+start+
					": expected="+expected+", actual="+counts.get(start));
		}
		return counts;
	}

	/** Calls the real concrete inner method and checks full bases, counts and unchanged qualities. */
	private static int runCase(final Tadpole tad, final byte[] truth, final int k, final int depth,
			final int position, final byte replacement, final String name, final boolean quality,
			final StringBuilder report){
		assert(position>=k+5 && position+k<truth.length) : "Substitution needs complete unaffected left and right count flanks.";
		final byte[] input=truth.clone();
		input[position]=replacement;
		final int expected=replacement=='A' ? 0 : 1;
		final IntList counts=verifyCounts(tad, input, k, depth, expected==0 ? -1 : position);
		final ByteBuilder bases=new ByteBuilder(input.length);
		bases.append(input);
		final byte[] quals=quality ? new byte[input.length] : null;
		if(quality){Arrays.fill(quals, (byte)30);}
		final int corrected=tad.reassemble_inner(bases, quals, new int[4], counts, 5, new Kmer(k), new Kmer(k));
		boolean basesMatch=bases.length()==truth.length, countsMatch=counts.size==truth.length-k+1, qualsMatch=true;
		for(int i=0; i<truth.length; i++){
			basesMatch&=bases.array[i]==truth[i];
			if(quality){qualsMatch&=quals[i]==30;}
		}
		for(int i=0; i<counts.size; i++){countsMatch&=counts.get(i)==depth;}
		final boolean passed=corrected==expected && basesMatch && countsMatch && qualsMatch;
		report.append(name).append('\t').append(quality ? "Q30" : "null").append('\t').append(expected)
				.append('\t').append(corrected).append('\t').append(basesMatch).append('\t').append(countsMatch)
				.append('\t').append(qualsMatch).append('\t').append(passed).append('\n');
		return passed ? 0 : 1;
	}

	/** Setup failures must remain fatal even with Java assertions disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Numeric base codes used to construct the literal A/C/G/T mutation cases. */
	private static final byte[] BASES={'A', 'C', 'G', 'T'};
}
