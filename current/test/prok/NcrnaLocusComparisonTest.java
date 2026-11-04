package test.prok;

import java.lang.management.ManagementFactory;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;

import prok.NcrnaCallerProtection;
import prok.NcrnaLocusComparison;
import prok.NcrnaLocusComparison.Locus;
import prok.NcrnaTextIO.Output;

/** Literal checks for one-to-one matching, augmenting paths and overlap geometry.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaLocusComparisonTest{

	public static void main(final String[] args) throws Exception{
		for(final String arg:args){
			if(arg.startsWith("-Xss")){
				assert(ManagementFactory.getRuntimeMXBean().getInputArguments().contains(arg)) :
					"The shell launcher must forward the requested JVM stack setting";
			}
		}
		checkCallerPolicy();
		checkProtectedOrder();
		final ArrayList<Locus> truth=new ArrayList<Locus>();
		truth.add(new Locus("contig", 1, 100, 0, null));
		truth.add(new Locus("contig", 101, 200, 0, null));
		truth.add(new Locus("contig", 201, 300, 1, null));
		final ArrayList<Locus> calls=new ArrayList<Locus>();
		calls.add(new Locus("contig", 51, 150, 0, "flexible"));
		calls.add(new Locus("contig", 1, 100, 0, "first_only"));
		calls.add(new Locus("contig", 1, 100, 0, "duplicate"));
		calls.add(new Locus("contig", 201, 300, 0, "wrong_strand"));
		calls.add(new Locus("another_contig", 201, 300, 1, "wrong_contig"));
		final int[] matches=NcrnaLocusComparison.match(truth, calls);
		assert(matches[0]>=1 && matches[0]<=2 && matches[1]==0 && matches[2]==-1) :
			"Augmentation must free the first reference for its specific call; duplicates cannot create a third match";
		calls.add(new Locus("contig", 201, 300, 1, "correct_strand"));
		final int[] complete=NcrnaLocusComparison.match(truth, calls);
		assert(complete[2]==5) : "The new same-contig, same-strand call must match the third reference";
		truth.clear();
		calls.clear();
		truth.add(new Locus("contig", 1, 5, 0, null));
		calls.add(new Locus("contig", 4, 8, 0, "two_bases"));
		assert(NcrnaLocusComparison.match(truth, calls)[0]==-1) : "Two of five bases is less than half the shorter span";
		calls.add(new Locus("contig", 3, 7, 0, "three_bases"));
		assert(NcrnaLocusComparison.match(truth, calls)[0]==1) : "Three of five bases satisfies the overlap rule";
		checkLongChain();
		checkCoordinates();
		System.err.println("Locus comparison checks passed: augmentation, duplicates, contigs, strands and odd-length overlap.");
	}

	/** External text retains strict positive signed-int coordinates and model tokens. */
	private static void checkCoordinates() throws Exception{
		final Path root=Files.createTempDirectory("locus-coordinates-");
		final String[] starts={"1", "0", "2147483648", "00000000001", "+1", " 1", "", "1.0"};
		try(AutoCloseable cleanup=()->{
			for(int i=0; i<starts.length; i++){Files.deleteIfExists(root.resolve("truth_"+i+".tsv"));}
			Files.deleteIfExists(root.resolve("calls.gff"));
			Files.delete(root);
		}){
			for(int i=0; i<starts.length; i++){
				final Path path=root.resolve("truth_"+i+".tsv");
				try(Output output=new Output(path.toString())){
					output.writer.print("contig\ttest\t"+starts[i]+"\t2147483647\t+\n");
				}
				try{
					final ArrayList<Locus> loci=NcrnaLocusComparison.readTruth(path.toString(), "test");
					assert(i==0 && loci.size()==1 && loci.get(0).start==1 && loci.get(0).stop==Integer.MAX_VALUE)
						: "Only the declared positive signed-int coordinate is valid: start="+starts[i];
				}catch(IllegalArgumentException failure){
					assert(i>0) : "The valid extreme coordinate must parse: "+failure;
				}
			}
			final Path calls=root.resolve("calls.gff");
			try(Output output=new Output(calls.toString())){
				output.writer.print("contig\tfixture\tRNA\t1\t4\t.\t+\t.\tfamily:test;model:model0 description words\n");
			}
			final ArrayList<Locus> parsed=NcrnaLocusComparison.readCalls(calls.toString(), "test");
			assert(parsed.size()==1 && parsed.get(0).model.equals("model0"))
				: "Caller provenance must retain the model identifier without its description";
		}
	}

	/** Adjacent flexible calls leave the last reference free. The final call can
	 * take only the first reference, so every earlier match must move one slot.
	 * The expected permutation follows geometry, independently of the search code.
	 */
	private static void checkLongChain(){
		final int count=4096;
		final ArrayList<Locus> truth=new ArrayList<Locus>(count);
		final ArrayList<Locus> calls=new ArrayList<Locus>(count);
		for(int i=0; i<count; i++){
			truth.add(new Locus("chain", 3*i+1, 3*i+2, 0, null));
			if(i+1<count){calls.add(new Locus("chain", 3*i+1, 3*i+5, 0, "flexible"));}
		}
		calls.add(new Locus("chain", 1, 2, 0, "first_only"));
		final int[] matches=NcrnaLocusComparison.match(truth, calls);
		assert(matches[0]==count-1) : "The final specific call must consume the first reference: "+matches[0];
		for(int i=1; i<count; i++){
			assert(matches[i]==i-1) : "The only complete matching shifts each flexible call right: reference="+i+", call="+matches[i];
		}
		System.err.println("Long-chain matching passed: "+count+" references matched in the literal expected permutation.");
	}

	/** Two call orders describe the same biology and must preserve the same protected set. */
	private static void checkProtectedOrder(){
		final ArrayList<Locus> truth=new ArrayList<Locus>();
		truth.add(new Locus("contig", 100, 200, 0, null));
		truth.add(new Locus("contig", 180, 280, 0, null));
		truth.add(new Locus("contig", 20, 120, 0, null));
		final ArrayList<Locus> calls=new ArrayList<Locus>();
		final Locus a=new Locus("contig", 60, 160, 0, "A");
		final Locus b=new Locus("contig", 140, 240, 0, "B");
		calls.add(a);
		calls.add(b);
		final int[] baseline=NcrnaLocusComparison.match(truth, calls);
		assert(baseline[0]>=0 && baseline[2]>=0 && baseline[1]<0) : "Literal baseline protects references zero and two";
		final boolean[] fixedMask={true, false, true};
		final int[] sameOrder=NcrnaLocusComparison.match(truth, calls, fixedMask);
		assert(sameOrder[0]>=0 && sameOrder[2]>=0) : "The fixed baseline mask must be preserved in original call order";
		calls.set(0, b);
		calls.set(1, a);
		final int[] reordered=NcrnaLocusComparison.match(truth, calls, fixedMask);
		assert(reordered[0]>=0 && reordered[2]>=0 && reordered[1]<0) : "Reordering identical calls must not create a protected loss";
		final int[] extended=NcrnaLocusComparison.match(truth, calls, new boolean[]{true, false, false});
		int matched=0;
		for(final int call:extended){if(call>=0){matched++;}}
		assert(extended[0]>=0 && matched==2) : "The all-reference pass must augment through a protected match while retaining it";
		calls.remove(1);
		final int[] loss=NcrnaLocusComparison.match(truth, calls, fixedMask);
		assert(loss[0]>=0 && loss[2]<0) : "Removing A really loses the protected left locus; priority must not fabricate its recovery";
	}

	/** Actual CallGenes aliases must not add unrequested genomes or side outputs. */
	private static void checkCallerPolicy(){
		NcrnaCallerProtection.validatePolicy(new String[]{"sam1=t", "taxonomy=t", "ncrnaidpass=0.70"});
		for(final String alias:new String[]{"infna", "fnain", "fna", "aminoout", "aaout", "amino",
				"16sout", "18sout", "itsout", "json_out", "json", "hist", "lengthhist", "lhist", "ncrnadiag", "hbmdiag"}){
			boolean rejected=false;
			try{
				NcrnaCallerProtection.validatePolicy(new String[]{alias+"=unowned_file"});
			}catch(IllegalArgumentException expected){rejected=true;}
			assert(rejected) : "Caller policy must not accept a driver-owned I/O alias: "+alias;
		}
	}
}
