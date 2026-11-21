package test.jgi;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.Objects;

import jgi.Assembly;
import parse.Parse;
import structures.IntList;

/**
 * Seven finite FASTA cases for Assembly's omitted final positive-length record.
 * Constructs the public production class once per case and observes package-private
 * state by reflection; no production access or reader settings are changed.
 * Both phases require exact literal metrics. Baseline expects five incomplete
 * length lists; fixed expects complete lists. An empty trailing header and empty
 * input are controls unchanged by flushing the final positive length.
 * Inputs contain no blank separator lines. Blank-line policy, error cleanup,
 * overflow, allocation limits and AssemblyStats3 output are outside this fixture.
 * @author Sucrose
 */
public final class AssemblyBoundaryTest{

	/** Requires phase=baseline|fixed and a new outdir whose parent already exists. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown Assembly fixture option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/assemblyboundarytest.sh with assertions enabled.");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		Files.createDirectory(root);
		final AssemblyBoundaryTest test=new AssemblyBoundaryTest(root, phase.equals("baseline"));
		for(final Row row : ROWS){test.run(row);}
		test.check(test.cases==7, "Frozen EOF matrix requires seven cases: "+test.cases);
		test.check(test.constructions==7, "Construct actual Assembly once per case: "+test.constructions);
		test.check(test.incompleteLists==(test.baseline ? 5 : 0),
			"Expected five baseline EOF omissions and zero fixed omissions: "+test.incompleteLists);
		test.check(test.checks==143, "Seven cases with twenty checks plus three summary checks: "+test.checks);
		final String summary="ASSEMBLY_BOUNDARY_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" constructions="+test.constructions+
			" incompleteLists="+test.incompleteLists+" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(root.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+root);}
		System.out.print(summary);
	}

	/** Binds observation-only reflection and creates one table retaining both oracles. */
	private AssemblyBoundaryTest(final Path root_, final boolean baseline_) throws Exception{
		root=root_; baseline=baseline_;
		length=field("length"); contigs=field("contigs"); composition=field("acgtnio");
		headerLength=field("headerLength"); firstHeader=field("firstHeader"); fname=field("fname");
		gc=Assembly.class.getDeclaredMethod("gc"); gc.setAccessible(true);
		atLeast=Assembly.class.getDeclaredMethod("lengthAtLeast", int.class); atLeast.setAccessible(true);
		Files.write(root.resolve("cases.tsv"), ("case\tinputBytes\tbaselineSize\tbaselineCount\tbaselineMax"+
			"\tbaselineLengths\tbaselineAtLeast\tfixedSize\tfixedCount\tfixedMax\tfixedLengths\tfixedAtLeast"+
			"\tactualSize\tactualCount\tactualMax\tactualLengths\tactualAtLeast"+
			"\texpectedACGTUNIO\tactualACGTUNIO\texpectedGC\tactualGC\texpectedHeaderBytes"+
			"\tactualHeaderBytes\texpectedFirstHeader\tactualFirstHeader\tinputUnchanged\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Writes a tiny retained literal input and observes the real constructor's results. */
	private void run(final Row row) throws Exception{
		final byte[] input=row.input.getBytes(StandardCharsets.US_ASCII);
		final Path file=root.resolve(row.name+".fa");
		Files.write(file, input, StandardOpenOption.CREATE_NEW);
		check(Arrays.equals(input, Files.readAllBytes(file)), row.name+": written input differs from literal bytes");
		cases++;
		constructions++;
		final Assembly assembly=new Assembly(file.toString());
		final long actualSize=length.getLong(assembly);
		final IntList list=(IntList)contigs.get(assembly);
		final int[] actualLengths=Arrays.copyOf(list.array, list.size);
		final int actualMax=list.size==0 ? 0 : list.get(0);
		final long[] actualAtLeast=new long[MINIMUMS.length];
		for(int i=0; i<MINIMUMS.length; i++){actualAtLeast[i]=(Long)atLeast.invoke(assembly, MINIMUMS[i]);}
		final long[] actualComposition=((long[])composition.get(assembly)).clone();
		final float actualGC=(Float)gc.invoke(assembly);
		final long actualHeaderBytes=headerLength.getLong(assembly);
		final String actualHeader=(String)firstHeader.get(assembly);
		final boolean unchanged=Arrays.equals(input, Files.readAllBytes(file));
		final Metrics expected=baseline ? row.before : row.complete;
		if(!Arrays.equals(actualLengths, row.complete.lengths)){incompleteLists++;}
		final String observed=row.name+"\t"+input.length+"\t"+row.before.tsv()+"\t"+row.complete.tsv()+
			"\t"+actualSize+"\t"+list.size+"\t"+actualMax+"\t"+Arrays.toString(actualLengths)+
			"\t"+Arrays.toString(actualAtLeast)+"\t"+Arrays.toString(row.bases)+"\t"+Arrays.toString(actualComposition)+
			"\t"+row.gc+"\t"+actualGC+"\t"+row.headerBytes+"\t"+actualHeaderBytes+
			"\t"+row.first+"\t"+actualHeader+"\t"+unchanged+"\n";
		Files.write(root.resolve("cases.tsv"), observed.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
		check(actualSize==expected.size, row.name+": size expected "+expected.size+", saw "+actualSize);
		check(list.size==expected.count, row.name+": count expected "+expected.count+", saw "+list.size);
		check(actualMax==expected.max, row.name+": maximum expected "+expected.max+", saw "+actualMax);
		check(Arrays.equals(actualLengths, expected.lengths), row.name+": descending lengths expected "+Arrays.toString(expected.lengths)+", saw "+Arrays.toString(actualLengths));
		for(int i=0; i<MINIMUMS.length; i++){
			check(actualAtLeast[i]==expected.thresholds[i], row.name+": lengthAtLeast("+MINIMUMS[i]+
				") expected "+expected.thresholds[i]+", saw "+actualAtLeast[i]);
		}
		check(Arrays.equals(actualComposition, row.bases), row.name+": A,C,G,T,U,N,IUPAC,OTHER expected "+Arrays.toString(row.bases)+", saw "+Arrays.toString(actualComposition));
		check(Float.floatToIntBits(actualGC)==Float.floatToIntBits(row.gc), row.name+": GC expected "+row.gc+", saw "+actualGC);
		check(actualHeaderBytes==row.headerBytes, row.name+": header bytes expected "+row.headerBytes+", saw "+actualHeaderBytes);
		check(Objects.equals(actualHeader, row.first), row.name+": first header expected "+row.first+", saw "+actualHeader);
		check(file.toString().equals(fname.get(assembly)), row.name+": constructor did not retain the supplied path");
		check(unchanged, row.name+": production construction changed input bytes");
	}

	/** Resolves a named production field without changing its contents or visibility. */
	private static Field field(final String name) throws Exception{
		final Field field=Assembly.class.getDeclaredField(name);
		field.setAccessible(true);
		return field;
	}

	/** Accumulates observable mismatches so each normal-returning case leaves evidence. */
	private void check(final boolean correct, final String message){
		checks++;
		if(!correct){failures++; System.err.println("FAIL: "+message);}
	}

	/** Rejects invalid invocation before interpreting production observations. */
	private static void require(final boolean correct, final String message){
		if(!correct){throw new IllegalArgumentException(message);}
	}

	/** Makes the sole 5,000-base literal-pattern input; never generates bulk data. */
	private static String longInput(){
		final char[] bases=new char[5000];
		Arrays.fill(bases, 'G');
		return ">a\nT\n>long\n"+new String(bases)+"\n";
	}

	/** Independent literal size, count, maximum, ordered lengths and threshold sums. */
	private static final class Metrics{
		Metrics(final long size_, final int count_, final int max_, final int[] lengths_, final long[] thresholds_){
			size=size_; count=count_; max=max_; lengths=lengths_; thresholds=thresholds_;
		}
		String tsv(){return size+"\t"+count+"\t"+max+"\t"+Arrays.toString(lengths)+"\t"+Arrays.toString(thresholds);}
		final long size;
		final int count, max;
		final int[] lengths;
		final long[] thresholds;
	}

	/** Input and metadata expectations; composition and headers are identical in both phases. */
	private static final class Row{
		Row(final String name_, final String input_, final Metrics before_, final Metrics complete_,
				final long[] bases_, final float gc_, final long headerBytes_, final String first_){
			name=name_; input=input_; before=before_; complete=complete_;
			bases=bases_; gc=gc_; headerBytes=headerBytes_; first=first_;
		}
		final String name, input, first;
		final Metrics before, complete;
		final long[] bases;
		final float gc;
		final long headerBytes;
	}

	private final Path root;
	private final boolean baseline;
	private final Field length, contigs, composition, headerLength, firstHeader, fname;
	private final Method gc, atLeast;
	private int cases, constructions, incompleteLists, checks, failures;
	/** Positions in every literal threshold array, in bases; includes each published Stats3 cutoff. */
	private static final int[] MINIMUMS={0, 1, 2, 3, 5000, 5001, 10000, 25000, 50000};
	private static final Metrics EMPTY=new Metrics(0, 0, 0, new int[]{}, new long[]{0, 0, 0, 0, 0, 0, 0, 0, 0});
	private static final Metrics SINGLE=new Metrics(8, 1, 8, new int[]{8}, new long[]{8, 8, 8, 8, 0, 0, 0, 0, 0});
	private static final Metrics FIRST=new Metrics(2, 1, 2, new int[]{2}, new long[]{2, 2, 2, 0, 0, 0, 0, 0, 0});
	private static final Metrics TWO=new Metrics(5, 2, 3, new int[]{3, 2}, new long[]{5, 5, 5, 3, 0, 0, 0, 0, 0});
	private static final Row[] ROWS={
		new Row("one-newline", ">one\nACGTUNRY\n", EMPTY, SINGLE,
			new long[]{1, 1, 1, 1, 1, 1, 2, 0}, 0.4f, 3, "one"),
		new Row("one-no-newline", ">one\nACGTUNRY", EMPTY, SINGLE,
			new long[]{1, 1, 1, 1, 1, 1, 2, 0}, 0.4f, 3, "one"),
		new Row("two", ">first\nAA\n>last\nCCG\n", FIRST, TWO,
			new long[]{2, 2, 1, 0, 0, 0, 0, 0}, 0.6f, 9, "first"),
		new Row("two-empty-header", ">first\nAA\n>last\nCCG\n>\n", TWO, TWO,
			new long[]{2, 2, 1, 0, 0, 0, 0, 0}, 0.6f, 9, "first"),
		new Row("multiline-final", ">multi\nAC\nGTU\nNRY\n", EMPTY, SINGLE,
			new long[]{1, 1, 1, 1, 1, 1, 2, 0}, 0.4f, 5, "multi"),
		new Row("final-5000", longInput(),
			new Metrics(1, 1, 1, new int[]{1}, new long[]{1, 1, 0, 0, 0, 0, 0, 0, 0}),
			new Metrics(5001, 2, 5000, new int[]{5000, 1}, new long[]{5001, 5001, 5000, 5000, 5000, 0, 0, 0, 0}),
			new long[]{0, 0, 5000, 1, 0, 0, 0, 0}, 5000f/5001f, 5, "a"),
		new Row("empty", "", EMPTY, EMPTY,
			new long[]{0, 0, 0, 0, 0, 0, 0, 0}, 0f, 0, null)
	};
}
