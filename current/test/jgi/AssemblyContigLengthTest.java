package test.jgi;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.InputStream;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import jgi.Assembly;
import parse.Parse;
import shared.Shared;
import structures.IntList;

/**
 * One clean public constructor supplies a tiny Assembly holder, then eight
 * reflective calls exercise its private addContigLength(int,int) arithmetic.
 * These are arithmetic-helper cases, not huge FASTA reads or failed constructors.
 * Baseline returns three exact negative wrapped lengths; fixed rejects those
 * additions with an input-named RuntimeException and no returned length.
 * Every added operand is at most three. No giant arrays, workers or OS failure
 * injection are used. Metric snapshots describe only the tiny successful holder
 * before and after the arithmetic calls, never a materialized huge assembly.
 * @author Sucrose
 */
public final class AssemblyContigLengthTest{

	/** Requires phase=baseline|fixed, a fresh outdir and the one-CPU launcher. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown contig-length fixture option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Use dev/assemblycontiglengthtest.sh with assertions enabled.");
		require(Shared.threads()==1, "The raw bootstrap requires one tool thread.");
		require(!new File(INPUT).exists(), "A real stdin.fa would override borrowed standard input.");
		require(ByteFile1.bufferlen>=PREFIX.length, "The bootstrap prefix must fit in the initial BF1 buffer.");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		Files.createDirectory(root);
		final AssemblyContigLengthTest test=new AssemblyContigLengthTest(root, phase.equals("baseline"));
		final Assembly holder=test.bootstrap();
		for(final Row row : ROWS){test.run(holder, row);}
		test.inspect(holder, "same-small-holder-after-arithmetic");
		test.check(test.cases==8, "Expected eight literal arithmetic cases: "+test.cases);
		test.check(test.reflectiveCalls==8, "Expected eight explicit reflective calls: "+test.reflectiveCalls);
		test.check(test.constructions==1, "Only one clean public constructor is exercised: "+test.constructions);
		test.check(test.negativeReturns==(test.baseline ? 3 : 0), "Unexpected negative returned lengths: "+test.negativeReturns);
		test.check(test.rejections==(test.baseline ? 0 : 3), "Unexpected arithmetic rejections: "+test.rejections);
		test.check(test.returns==(test.baseline ? 8 : 5), "Unexpected arithmetic returns: "+test.returns);
		test.check(test.metricSnapshots==2, "Expected two snapshots of the same tiny holder: "+test.metricSnapshots);
		final String summary="ASSEMBLY_CONTIG_LENGTH_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" publicBootstraps="+test.constructions+
			" reflectiveCalls="+test.reflectiveCalls+" negativeReturns="+test.negativeReturns+
			" rejections="+test.rejections+" returns="+test.returns+" metricSnapshots="+test.metricSnapshots+
			" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(root.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+root);}
		System.out.print(summary);
	}

	/** Binds private observations and the actual production arithmetic helper. */
	private AssemblyContigLengthTest(final Path root_, final boolean baseline_) throws Exception{
		root=root_; baseline=baseline_;
		add=Assembly.class.getDeclaredMethod("addContigLength", Integer.TYPE, Integer.TYPE);
		require(Modifier.isPrivate(add.getModifiers()) && !Modifier.isStatic(add.getModifiers()) &&
			add.getReturnType()==Integer.TYPE, "Expected private instance int addContigLength(int,int).");
		add.setAccessible(true);
		length=field("length"); contigs=field("contigs"); bases=field("acgtnio");
		headerLength=field("headerLength"); firstHeader=field("firstHeader"); fname=field("fname");
		gc=Assembly.class.getDeclaredMethod("gc"); gc.setAccessible(true);
		Files.write(root.resolve("prefix.fa"), PREFIX, StandardOpenOption.CREATE_NEW);
		Files.write(root.resolve("cases.tsv"), ("case\tcoverage\tphase\tcurrent\tadded\texpectedOutcome"+
			"\texpectedReturn\tactualReturn\texceptionClass\texceptionMessage\tsuppressedCount\tglobalsUnchanged\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		Files.write(root.resolve("metrics.tsv"), ("stage\tsize\tcount\tlengths\tmax\tbases\tgc\theaderBytes\tfirst\tfname\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Runs the sole real constructor with a copied clean prefix and restores all bootstrap globals. */
	private Assembly bootstrap() throws Exception{
		final InputStream originalIn=System.in;
		final PrintStream originalErr=System.err;
		final boolean force=ByteFile.FORCE_MODE_BF1, simd=Shared.SIMD, buffered=ByteFile1.BUFFERED;
		final BootstrapInput input=new BootstrapInput();
		final ByteArrayOutputStream diagnostics=new ByteArrayOutputStream();
		final PrintStream capture=new PrintStream(diagnostics, true, "UTF-8");
		final Assembly holder;
		try{
			ByteFile.FORCE_MODE_BF1=true;
			Shared.SIMD=false;
			ByteFile1.BUFFERED=false;
			System.setIn(input);
			System.setErr(capture);
			constructions++;
			holder=new Assembly(INPUT);
		}finally{
			System.setIn(originalIn);
			System.setErr(originalErr);
			ByteFile.FORCE_MODE_BF1=force;
			Shared.SIMD=simd;
			ByteFile1.BUFFERED=buffered;
			capture.close();
		}
		final boolean restored=globalsMatch(originalIn, originalErr, force, simd, buffered);
		check(input.closeCalls==0, "Bootstrap standard input is borrowed: "+input.closeCalls);
		check(input.available()==0, "Bootstrap prefix was not fully consumed: "+input.available());
		check(input.unchanged(), "Bootstrap prefix storage changed");
		check(restored, "Bootstrap did not restore streams and reader flags");
		check(diagnostics.size()==0, "Clean bootstrap emitted diagnostics: "+diagnostics.toString("UTF-8"));
		Files.write(root.resolve("bootstrap.stderr"), diagnostics.toByteArray(), StandardOpenOption.CREATE_NEW);
		Files.write(root.resolve("bootstrap.tsv"), ("coverage\tpublicConstructions\tborrowedCloseCalls\tglobalsRestored\tprefixUnchanged\n"+
			"clean-public-constructor\t"+constructions+"\t"+input.closeCalls+"\t"+restored+"\t"+input.unchanged()+"\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		inspect(holder, "clean-public-bootstrap");
		return holder;
	}

	/** Checks an exact literal sum or rejection without constructing an enormous sequence. */
	private void run(final Assembly holder, final Row row) throws Exception{
		final InputStream originalIn=System.in;
		final PrintStream originalErr=System.err;
		final boolean force=ByteFile.FORCE_MODE_BF1, simd=Shared.SIMD, buffered=ByteFile1.BUFFERED;
		Object result=null;
		Throwable error=null;
		cases++; reflectiveCalls++;
		try{result=add.invoke(holder, row.current, row.added);}
		catch(final InvocationTargetException caught){error=caught.getCause();}
		final boolean returned=error==null && result instanceof Integer;
		final boolean reject=!baseline && row.overflow;
		final boolean globals=globalsMatch(originalIn, originalErr, force, simd, buffered);
		if(returned){returns++; if(((Integer)result).intValue()<0){negativeReturns++;}}
		if(error!=null){rejections++;}
		check(reject ? error!=null && result==null : returned, row.name+": unexpected return/rejection: "+result+" / "+error);
		check(reject ? error!=null && error.getClass()==RuntimeException.class &&
			String.valueOf(error.getMessage()).contains(INPUT) && String.valueOf(error.getMessage()).contains("2147483647") :
			returned && ((Integer)result).intValue()==row.expectedBaseline,
			row.name+": expected "+(reject ? "input-named overflow rejection" : Integer.toString(row.expectedBaseline))+", saw "+result+" / "+error);
		check(error==null || error.getSuppressed().length==0, row.name+": arithmetic-only rejection acquired suppressed errors");
		check(globals, row.name+": arithmetic helper changed System streams or reader flags");
		final String observed=row.name+"\tprivate-addContigLength-reflection\t"+(baseline ? "baseline" : "fixed")+
			"\t"+row.current+"\t"+row.added+"\t"+(reject ? "reject" : "return")+
			"\t"+(reject ? "unavailable" : Integer.toString(row.expectedBaseline))+"\t"+(result==null ? "unavailable" : result.toString())+
			"\t"+(error==null ? "none" : error.getClass().getName())+"\t"+message(error)+
			"\t"+(error==null ? 0 : error.getSuppressed().length)+"\t"+globals+"\n";
		Files.write(root.resolve("cases.tsv"), observed.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
	}

	/** Checks nine literal metrics of the tiny holder, not metrics for the large arithmetic operands. */
	private void inspect(final Assembly holder, final String stage) throws Exception{
		metricSnapshots++;
		final long size=length.getLong(holder), headers=headerLength.getLong(holder);
		final IntList list=(IntList)contigs.get(holder);
		final int[] lengths=Arrays.copyOf(list.array, list.size);
		final int max=list.size==0 ? 0 : list.get(0);
		final long[] composition=((long[])bases.get(holder)).clone();
		final float fraction=(Float)gc.invoke(holder);
		final String first=(String)firstHeader.get(holder), path=(String)fname.get(holder);
		check(size==2, stage+": expected two bases: "+size);
		check(list.size==1, stage+": expected one contig: "+list.size);
		check(Arrays.equals(lengths, new int[]{2}), stage+": expected lengths [2]: "+Arrays.toString(lengths));
		check(max==2, stage+": expected maximum two: "+max);
		check(Arrays.equals(composition, new long[]{1, 1, 0, 0, 0, 0, 0, 0}), stage+": expected A1/C1: "+Arrays.toString(composition));
		check(Float.floatToIntBits(fraction)==Float.floatToIntBits(0.5f), stage+": expected GC0.5: "+fraction);
		check(headers==2, stage+": expected two header bytes: "+headers);
		check("ok".equals(first), stage+": expected first header ok: "+first);
		check(INPUT.equals(path), stage+": expected input stdin.fa: "+path);
		final String observed=stage+"\t"+size+"\t"+list.size+"\t"+Arrays.toString(lengths)+"\t"+max+
			"\t"+Arrays.toString(composition)+"\t"+fraction+"\t"+headers+"\t"+first+"\t"+path+"\n";
		Files.write(root.resolve("metrics.tsv"), observed.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
	}

	/** Resolves a private observation field without changing production visibility. */
	private static Field field(final String name) throws Exception{
		final Field value=Assembly.class.getDeclaredField(name);
		value.setAccessible(true);
		return value;
	}
	/** Observes the process globals used only by the public bootstrap. */
	private static boolean globalsMatch(final InputStream in, final PrintStream err, final boolean force, final boolean simd, final boolean buffered){
		return System.in==in && System.err==err && ByteFile.FORCE_MODE_BF1==force && Shared.SIMD==simd && ByteFile1.BUFFERED==buffered;
	}
	/** Keeps exception text in a single TSV cell. */
	private static String message(final Throwable error){
		return error==null ? "" : String.valueOf(error.getMessage()).replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');
	}
	/** Accumulates independent mismatches without skipping later arithmetic rows. */
	private void check(final boolean correct, final String message){
		checks++;
		if(!correct){failures++; System.err.println("FAIL: "+message);}
	}
	/** Rejects fixture misuse. */
	private static void require(final boolean correct, final String message){
		if(!correct){throw new IllegalArgumentException(message);}
	}

	/** A tiny copied clean input whose close count observes borrowed-stream ownership only. */
	private static final class BootstrapInput extends ByteArrayInputStream{
		BootstrapInput(){super(PREFIX.clone());}
		@Override
		public void close(){closeCalls++;}
		boolean unchanged(){return Arrays.equals(PREFIX, buf);}
		int closeCalls;
	}

	/** Literal operands and an independently stated baseline return value. */
	private static final class Row{
		Row(final String name_, final int current_, final int added_, final int expectedBaseline_, final boolean overflow_){
			require(current_>=0 && added_>=0 && added_<=3, "Fixture arithmetic uses a nonnegative current and at most three added bases: "+name_);
			name=name_; current=current_; added=added_; expectedBaseline=expectedBaseline_; overflow=overflow_;
		}
		final String name;
		final int current, added, expectedBaseline;
		final boolean overflow;
	}

	private final Path root;
	private final boolean baseline;
	private final Method add, gc;
	private final Field length, contigs, bases, headerLength, firstHeader, fname;
	private int cases, reflectiveCalls, constructions, negativeReturns, rejections, returns, metricSnapshots, checks, failures;
	private static final String INPUT="stdin.fa";
	private static final byte[] PREFIX=">ok\nAC\n".getBytes(StandardCharsets.US_ASCII);
	private static final Row[] ROWS={
		new Row("zero", 0, 0, 0, false),
		new Row("ordinary", 6, 2, 8, false),
		new Row("below-limit", 2147483645, 1, 2147483646, false),
		new Row("exact-limit", 2147483646, 1, 2147483647, false),
		new Row("limit-plus-zero", 2147483647, 0, 2147483647, false),
		new Row("limit-plus-one", 2147483647, 1, -2147483648, true),
		new Row("near-limit-plus-two", 2147483646, 2, -2147483648, true),
		new Row("near-limit-plus-three", 2147483646, 3, -2147483647, true)
	};
}
