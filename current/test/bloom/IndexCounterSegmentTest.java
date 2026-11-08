package test.bloom;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import bloom.IndexCounter;
import bloom.KCountArray3;
import bloom.KmerCountAbstract;
import dna.AminoAcid;
import dna.ChromosomeArray;
import dna.Data;
import shared.Shared;

/**
 * Regression checks for IndexCounter fixed-segment boundary accounting.
 * Runs sequentially in a fresh JVM and installs trimmed in-memory reference
 * arrays rather than building an on-disk reference index.
 *
 * @author Brian Bushnell, Jean
 */
public final class IndexCounterSegmentTest{

	/** Runs the exact-table segment accounting cases and saves every raw outcome. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1 || !args[0].startsWith("outdir=")){
			throw new IllegalArgumentException("Usage: outdir=<new directory>");
		}
		final Path root=Paths.get(args[0].substring(7)).toAbsolutePath();
		check(!Files.exists(root), "Result directory must be new: "+root);
		Files.createDirectories(root);

		assertKey("ACG", ACG, 6);
		assertKey("ACT", ACT, 7);
		assertKey("CGT", CGT, 27);
		assertKey("CTT", CTT, 31);
		assertKey("GTA", GTA, 44);
		assertKey("TAC", TAC, 49);

		final CaseSpec[] cases={
				new CaseSpec("terminal-remainder", "NACGTACTT",
						expected(ACG, CGT, GTA, TAC, ACT, CTT), 6),
				new CaseSpec("divisible", "NACGTACT",
						expected(ACG, CGT, GTA, TAC, ACT), 5),
				new CaseSpec("one-window", "NACG",
						expected(ACG), 1)
		};

		int failures=0;
		for(CaseSpec spec : cases){failures+=checkCase(root, spec);}
		if(failures>0){
			throw new AssertionError("IndexCounter segment test failed "+failures+
					" of "+cases.length+" cases.");
		}
		System.out.println("IndexCounterSegmentTest PASS: "+checks+" checks");
	}

	/** Counts one literal chromosome, saves its raw outcome, and checks every k=3 raw key. */
	private static int checkCase(final Path root, final CaseSpec spec) throws Exception{
		final Observation actual=countLiteralChromosome(spec.bases);
		final ArrayList<String> failures=new ArrayList<String>();
		for(int key=0; key<spec.expected.length; key++){
			checks++;
			if(actual.counts[key]!=spec.expected[key]){
				failures.add("wrong count for key "+key+"; expected "+spec.expected[key]+
						", observed "+actual.counts[key]);
			}
		}
		checks++;
		if(actual.positiveKeys!=spec.expectedReads){
			failures.add("wrong distinct table key count; expected "+spec.expectedReads+
					", observed "+actual.positiveKeys);
		}
		checks++;
		if(actual.readsProcessed!=spec.expectedReads){
			failures.add("wrong accepted kmer window count; expected "+spec.expectedReads+
					", observed "+actual.readsProcessed);
		}
		checks++;
		if(!actual.inputPreserved){failures.add("literal input bytes were not preserved: "+actual);}
		checks++;
		if(!actual.restored){failures.add("shared globals were not restored: "+actual);}
		writeCase(root.resolve(spec.name+".txt"), spec, actual, failures);
		System.out.println(spec.name+" pass="+failures.isEmpty()+" reads="+
				actual.readsProcessed+" positiveKeys="+actual.positiveKeys+
				" restored="+actual.restored+" inputPreserved="+actual.inputPreserved);
		return failures.isEmpty() ? 0 : 1;
	}

	/** Runs the real IndexCounter worker path on a single in-memory chromosome. */
	private static Observation countLiteralChromosome(final String bases){
		final int oldThreads=Shared.threads();
		final int oldBuffers=Shared.READ_BUFFER_NUM_BUFFERS();
		final int oldNumChroms=Data.numChroms;
		final ChromosomeArray[] oldChromosomes=Data.chromosomePlusMatrix;
		final long oldReads=KmerCountAbstract.readsProcessed;
		final long oldKeys=KmerCountAbstract.keysCounted;
		final long oldIncrements=KmerCountAbstract.increments;
		final boolean oldCanonical=KmerCountAbstract.CANONICAL;
		Observation actual=null;
		try{
			Shared.setThreads(1);
			Data.numChroms=1;
			final ChromosomeArray ca=makeChromosome(bases);
			Data.chromosomePlusMatrix=new ChromosomeArray[]{null, ca};
			KmerCountAbstract.readsProcessed=0;
			KmerCountAbstract.keysCounted=0;
			KmerCountAbstract.increments=0;
			KmerCountAbstract.CANONICAL=true;

			final KCountArray3 counts=new KCountArray3(2048, 4);
			new IndexCounter(3, false).countFromIndex(counts);
			actual=new Observation(bases, ca, counts);
		}catch(Exception e){
			throw new RuntimeException("IndexCounter countFromIndex failed for "+bases, e);
		}finally{
			Data.numChroms=oldNumChroms;
			Data.chromosomePlusMatrix=oldChromosomes;
			KmerCountAbstract.readsProcessed=oldReads;
			KmerCountAbstract.keysCounted=oldKeys;
			KmerCountAbstract.increments=oldIncrements;
			KmerCountAbstract.CANONICAL=oldCanonical;
			Shared.setThreads(oldThreads);
			Shared.setBuffers(oldBuffers);
			if(actual!=null){
				actual.recordRestoration(oldThreads, oldBuffers, oldNumChroms, oldChromosomes,
						oldReads, oldKeys, oldIncrements, oldCanonical);
			}
		}
		return actual;
	}

	/** Installs the literal indexed-array shape required by Data.getChromosome. */
	private static ChromosomeArray makeChromosome(final String bases){
		final ChromosomeArray ca=new ChromosomeArray(1, Shared.PLUS);
		ca.set(0, bases);
		ca.resize(ca.maxIndex+1);
		check(ca.minIndex==0, "ChromosomeArray must start at zero for IndexCounter: minIndex="+
				ca.minIndex+", bases="+bases);
		check(ca.maxIndex==bases.length()-1,
				"ChromosomeArray maxIndex must match the final installed base for IndexCounter: "+
				"maxIndex="+ca.maxIndex+", bases="+bases);
		check(ca.array.length==bases.length(),
				"IndexCounter scans ChromosomeArray.array.length directly; expected exact converter "+
				"trim length "+bases.length()+", observed "+ca.array.length);
		check(Arrays.equals(ca.array, bases.getBytes(StandardCharsets.US_ASCII)),
				"Installed ChromosomeArray bytes changed before counting: expected "+bases+
				", observed "+new String(ca.array, StandardCharsets.US_ASCII));
		return ca;
	}

	/** Returns a complete expected table for all cells in the exact KCountArray3. */
	private static int[] expected(final int... presentKeys){
		final int[] expected=new int[CELLS];
		for(int key : presentKeys){expected[key]=1;}
		return expected;
	}

	/** Verifies the documented literal k=3 key values used in the table oracle. */
	private static void assertKey(final String bases, final int expected, final int literal){
		check(expected==literal, "Literal key for "+bases+" changed; expected "+
				literal+", observed "+expected);
	}

	/** Saves the full actual table before a fixed-oracle mismatch can fail the run. */
	private static void writeCase(final Path path, final CaseSpec spec, final Observation actual,
			final List<String> failures) throws Exception{
		final ArrayList<String> lines=new ArrayList<String>();
		lines.add("name="+spec.name);
		lines.add("input="+spec.bases);
		lines.add("pass="+failures.isEmpty());
		lines.add("expectedReads="+spec.expectedReads);
		lines.add("actualReads="+actual.readsProcessed);
		lines.add("actualPositiveKeys="+actual.positiveKeys);
		lines.add("actualKeysCounted="+actual.keysCounted);
		lines.add("actualIncrements="+actual.increments);
		lines.add("installedBytes="+actual.installedBytes);
		lines.add("inputPreserved="+actual.inputPreserved);
		lines.add("globalsRestored="+actual.restored);
		lines.add("expectedCounts="+Arrays.toString(spec.expected));
		lines.add("actualCounts="+Arrays.toString(actual.counts));
		lines.add("failures="+failures.size());
		lines.addAll(failures);
		Files.write(path, lines, StandardCharsets.UTF_8);
	}

	/** Unit failures remain fatal without relying on the launcher's assertion flag. */
	private static void check(final boolean condition, final String message){
		checks++;
		if(!condition){throw new AssertionError(message);}
	}

	/** Raw table counts and diagnostics saved before any case assertion can fail. */
	private static final class Observation{

		/** Captures exact k=3 counts and shared diagnostics immediately after counting. */
		Observation(final String input_, final ChromosomeArray ca, final KCountArray3 counts_){
			input=input_;
			installedBytes=new String(ca.array, StandardCharsets.US_ASCII);
			inputPreserved=input.equals(installedBytes);
			counts=new int[CELLS];
			int positive=0;
			for(int i=0; i<counts.length; i++){
				counts[i]=counts_.read(i);
				if(counts[i]>0){positive++;}
			}
			positiveKeys=positive;
			readsProcessed=KmerCountAbstract.readsProcessed;
			keysCounted=KmerCountAbstract.keysCounted;
			increments=KmerCountAbstract.increments;
		}

		/** Records whether shared globals were restored before the observation escaped. */
		void recordRestoration(final int oldThreads, final int oldBuffers, final int oldNumChroms,
				final ChromosomeArray[] oldChromosomes, final long oldReads, final long oldKeys,
				final long oldIncrements, final boolean oldCanonical){
			restored=Shared.threads()==oldThreads &&
					Shared.READ_BUFFER_NUM_BUFFERS()==oldBuffers &&
					Data.numChroms==oldNumChroms &&
					Data.chromosomePlusMatrix==oldChromosomes &&
					KmerCountAbstract.readsProcessed==oldReads &&
					KmerCountAbstract.keysCounted==oldKeys &&
					KmerCountAbstract.increments==oldIncrements &&
					KmerCountAbstract.CANONICAL==oldCanonical;
		}

		/** Renders the full saved observation for boundary-failure diagnosis. */
		@Override
		public String toString(){
			return "input="+input+
					", ACG="+counts[ACG]+
					", ACT="+counts[ACT]+
					", CGT="+counts[CGT]+
					", CTT="+counts[CTT]+
					", GTA="+counts[GTA]+
					", TAC="+counts[TAC]+
					", positiveKeys="+positiveKeys+
					", readsProcessed="+readsProcessed+
					", keysCounted="+keysCounted+
					", increments="+increments+
					", installedBytes="+installedBytes+
					", inputPreserved="+inputPreserved+
					", restored="+restored;
		}

		final String input;
		final String installedBytes;
		final int[] counts;
		final int positiveKeys;
		final long readsProcessed;
		final long keysCounted;
		final long increments;
		final boolean inputPreserved;
		boolean restored;
	}

	/** One literal input and its exact post-fix count oracle. */
	private static final class CaseSpec{

		/** Stores the case name, installed bases and fixed expected exact table. */
		CaseSpec(final String name_, final String bases_, final int[] expected_,
				final long expectedReads_){
			name=name_;
			bases=bases_;
			expected=expected_;
			expectedReads=expectedReads_;
		}

		final String name;
		final String bases;
		final int[] expected;
		final long expectedReads;
	}

	private static final int ACG=(int)AminoAcid.stringToKmer("ACG");
	private static final int ACT=(int)AminoAcid.stringToKmer("ACT");
	private static final int CGT=(int)AminoAcid.stringToKmer("CGT");
	private static final int CTT=(int)AminoAcid.stringToKmer("CTT");
	private static final int GTA=(int)AminoAcid.stringToKmer("GTA");
	private static final int TAC=(int)AminoAcid.stringToKmer("TAC");
	private static final int CELLS=2048;

	private static long checks=0;
}
