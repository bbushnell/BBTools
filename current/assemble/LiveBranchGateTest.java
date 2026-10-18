package assemble;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Random;
import java.util.concurrent.atomic.AtomicReference;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import ml.CellNet;
import ml.CellNetParser;
import shared.Shared;
import structures.ByteBuilder;

/** Guard, thread isolation and scalar inference checks for the opt-in live gate. @author Fischl */
public final class LiveBranchGateTest {

	/** Writes constant-model fixtures and optionally checks the frozen candidate model. */
	public static void main(final String[] args) throws InterruptedException{
		if(args.length==3 && args[0].equals("truth")){truth(args[1], args[2]); return;}
		if(args.length==3 && args[0].equals("reverse")){reverseReads(args[1], args[2]); return;}
		if(args.length<1 || args.length>2 || !new File(args[0]).mkdir()){
			throw new IllegalArgumentException("Supply NEW_fixture_directory [candidate.bbnet].");
		}
		final boolean oldFma=Shared.SIMD_FMA;
		try{
			Shared.SIMD_FMA=false;
			final SameKJoinFeatures.Moments moments=new SameKJoinFeatures.Moments();
			moments.add(10); moments.add(20);
			final CellNet yes=constant(args[0]+"/accept.bbnet", 1), no=constant(args[0]+"/reject.bbnet", 0);
			final LiveBranchGate gate=new LiveBranchGate(yes, 31, 2, 1, moments);
			final int[] left={10, 2, 0, 0}, right={12, 2, 0, 0};
			check(!gate.rescue(10, left, right, false, false, 0), "Native accepts must bypass rescue.");
			check(!gate.rescue(1, left, right, true, true, 0), "Node depth below mce cannot be rescued.");
			check(!gate.rescue(10, left, new int[]{0, 0, 0, 0}, true, false, 0), "Missing preferred continuation cannot be rescued.");
			check(!gate.rescue(10, new int[]{10, 10, 0, 0}, right, true, true, 0), "Equal maxima cannot be rescued.");
			check(!gate.rescue(10, new int[]{10, 5, 0, 0}, right, true, true, 0), "Half-depth competitor is a hard stop.");
			check(!gate.rescue(10, left, right, true, true, 1), "Wrong actual predecessor cannot enter the preferred route.");
			check(gate.totals()[0]==0, "Ineligible calls must never evaluate a network.");
			check(gate.rescue(10, left, right, true, true, 0), "An eligible route at the cutoff is accepted.");
			check(gate.rescue(10, right, left, true, true, 0), "Swapping route sides preserves the decision.");
			check(gate.rescue(10, new int[]{0, 0, 2, 12}, new int[]{0, 0, 2, 10}, true, true, 3), "Reverse-complement traversal preserves the decision.");
			check(gate.totals()[0]==3 && gate.totals()[1]==3, "Guarded score counters disagree.");
			final LiveBranchGate reject=new LiveBranchGate(no, 31, 2, 0.5f, moments);
			check(!reject.rescue(10, left, right, true, true, 0), "Low network scores must retain native stops.");
			check(reject.totals()[0]==1 && reject.totals()[1]==0, "Rejected scores must still count as evaluations.");
			if(args.length==2){parity(args[1], moments);}
			System.err.println("LIVE_BRANCH_GATE_TEST_PASS guarded_calls=10 candidate_parity="+(args.length==2));
		}finally{Shared.SIMD_FMA=oldFma;}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Every fixture contig must be an exact substring of either reference strand. */
	private static void truth(final String referencePath, final String assemblyPath){
		final ArrayList<Contig> references=TadpoleGraph.readContigs(referencePath);
		final ArrayList<Contig> contigs=TadpoleGraph.readContigs(assemblyPath);
		check(references.size()==1 && !contigs.isEmpty(), "The exact fixture requires one reference and nonempty assembly.");
		final String reference=new String(references.get(0).bases, StandardCharsets.US_ASCII);
		long bases=0;
		int longest=0;
		for(Contig contig : contigs){
			final String sequence=new String(contig.bases, StandardCharsets.US_ASCII);
			check(reference.contains(sequence) || reference.contains(reverseComplement(sequence)),
					"Fixture contig is not an exact reference substring: "+contig.id+" length="+sequence.length());
			bases+=sequence.length(); longest=Math.max(longest, sequence.length());
		}
		System.err.println("LIVE_BRANCH_TRUTH_PASS contigs="+contigs.size()+" bases="+bases+" longest="+longest);
	}

	/** Reverses input order and strand to perturb seeding without changing counts. */
	private static void reverseReads(final String input, final String output){
		final ArrayList<Contig> reads=TadpoleGraph.readContigs(input);
		check(!reads.isEmpty() && !new File(output).exists(), "Reverse fixture needs reads and a new output.");
		final ByteStreamWriter writer=new ByteStreamWriter(output, false, false, false);
		writer.start();
		final ByteBuilder text=new ByteBuilder();
		for(int i=reads.size()-1; i>=0; i--){
			text.clear().append('>').append(i).nl();
			text.append(reverseComplement(new String(reads.get(i).bases, StandardCharsets.US_ASCII))).nl();
			writer.print(text);
		}
		if(writer.poisonAndWait()){throw new IllegalStateException("Reverse fixture write failed.");}
	}

	/** Literal fixture-only reverse complement, independent of graph orientation code. */
	private static String reverseComplement(final String sequence){
		final StringBuilder reverse=new StringBuilder(sequence.length());
		for(int i=sequence.length()-1; i>=0; i--){
			final int base="ACGT".indexOf(sequence.charAt(i));
			check(base>=0, "Fixture sequence must be unambiguous DNA.");
			reverse.append("TGCA".charAt(base));
		}
		return reverse.toString();
	}

	/** Builds explicit dense12-to1 constant functions for assembly and eligibility fixtures. */
	private static CellNet constant(final String path, final int value){
		assert(value==0 || value==1) : "Constant fixtures exercise the two probability endpoints.";
		final ByteBuilder text=new ByteBuilder("##bbnet\n#version 1\n#concise\n#dense\n#density 1\n#blocksize 1\n#seed 1\n#layers 2\n#dims 12 1\n#edges 12\n#coding decimal\n##layer 1\nC13 LINEAR ");
		text.append(value);
		for(int i=0; i<12; i++){text.append(" 0");}
		text.nl();
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start(); writer.print(text);
		if(writer.poisonAndWait()){throw new IllegalStateException("Constant fixture write failed.");}
		return CellNetParser.loadInferenceFromLines(ByteFile.toLines(path));
	}

	/** Compares private inference with the evaluator's legacy dense path on10000 generated inputs. */
	private static void parity(final String path, final SameKJoinFeatures.Moments moments) throws InterruptedException{
		final CellNet reference=CellNetParser.load(path, false);
		reference.simdFF=false;
		final CellNet isolated=CellNetParser.loadInferenceFromLines(ByteFile.toLines(path));
		final LiveBranchGate gate=new LiveBranchGate(isolated, 31, 2, 0.52372265f, moments);
		final float[][] vectors=new float[10000][12];
		final int[] expected=new int[vectors.length];
		final Random random=new Random(3801);
		final int[] left=new int[4], right=new int[4];
		for(int i=0; i<vectors.length; i++){
			for(int j=0; j<4; j++){left[j]=random.nextInt(300); right[j]=random.nextInt(300);}
			LiveBranchPolicy.fill(vectors[i], 24+random.nextInt(104), random.nextInt(300), left, right,
					1+random.nextDouble()*150, 151+random.nextDouble()*150);
			expected[i]=Float.floatToRawIntBits(reference.applyInput(vectors[i]).feedForwardDense());
		}
		// A later ordinary sparse load must not change this gate's dispatch or copies.
		final ArrayList<byte[]> sparseLines=new ArrayList<byte[]>();
		for(String line:("##bbnet\n#version 1\n#concise\n#coding decimal\n#seed 23\n"+
				"#layers 2\n#blocksize 1\n#dims 3 1\n#sparse\nI4 0 2\nC4 LINEAR 0.25 0.5 -0.25\n").split("\n")){
			sparseLines.add(line.getBytes(StandardCharsets.US_ASCII));
		}
		final CellNet sparse=CellNetParser.loadFromLines(sparseLines);
		check(!sparse.dense && sparse.applyInput(new float[]{2, 1, 4}).feedForward()==0.25f,
				"The unrelated model must exercise a real sparse load with gapped inputs.");
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final ArrayList<Thread> workers=new ArrayList<Thread>();
		for(int worker=0; worker<4; worker++){
			final Thread thread=new Thread(new Runnable(){
				@Override public void run(){
					try{
						for(int i=0; i<vectors.length; i++){
							if(Float.floatToRawIntBits(gate.score(vectors[i]))!=expected[i]){
								throw new AssertionError("Private/scalar score mismatch at generated row "+i);
							}
						}
					}catch(Throwable problem){failure.compareAndSet(null, problem);}
				}
			});
			workers.add(thread); thread.start();
		}
		for(Thread thread : workers){thread.join();}
		if(failure.get()!=null){throw new AssertionError("Concurrent candidate inference failed.", failure.get());}
		System.err.println("LIVE_BRANCH_GATE_PARITY_PASS vectors=10000 workers=4 exact_float_bits=true unrelated_sparse_load=true");
	}

	/** Keeps test failures active even if the launcher accidentally loses assertions. */
	private static void check(final boolean pass, final String message){
		if(!pass){throw new AssertionError(message);}
	}
}
