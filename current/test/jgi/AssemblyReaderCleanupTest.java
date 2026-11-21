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
import fileIO.FileFormat;
import jgi.Assembly;
import parse.Parse;
import shared.Shared;
import structures.IntList;

/**
 * One clean public constructor bootstraps an Assembly, followed by six reflective
 * calls to its private load(ByteFile) helper. These are helper ownership checks,
 * not six failed public constructors. The tracked readers own no OS resource or
 * worker; close counts prove attempts, not real-stream closure or factory safety.
 * Baseline is the already checked-close implementation before exceptional cleanup:
 * three read-exception rows skip close. Fixed must close once, preserve the exact
 * primary exception and suppress a flagged or thrown secondary cleanup failure.
 * Normal close=false, normal close=true and normal throwing-close are controls.
 * Each row has fresh exceptions; failed calls have no reported metrics.
 * @author Sucrose
 */
public final class AssemblyReaderCleanupTest{

	/** Requires a fresh output directory and the one-CPU development launcher. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown reader-cleanup fixture option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Use dev/assemblyreadercleanuptest.sh with assertions enabled.");
		require(Shared.threads()==1, "The bootstrap requires one tool thread.");
		require(!new File(INPUT).exists(), "A real stdin.fa would override borrowed standard input.");
		require(ByteFile1.bufferlen>=PREFIX.length, "The bootstrap prefix must fit in the initial BF1 buffer.");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		Files.createDirectory(root);
		final AssemblyReaderCleanupTest test=new AssemblyReaderCleanupTest(root, phase.equals("baseline"));
		final Assembly holder=test.bootstrap();
		for(int row=0; row<6; row++){test.run(holder, row);}
		test.check(test.cases==6, "Expected six private-helper rows: "+test.cases);
		test.check(test.helperCalls==6, "Expected six private-helper invocations: "+test.helperCalls);
		test.check(test.constructions==1, "Only one clean public constructor is exercised: "+test.constructions);
		test.check(test.metricSnapshots==2, "Metrics belong to the bootstrap and successful helper only: "+test.metricSnapshots);
		final String summary="ASSEMBLY_READER_CLEANUP_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" publicBootstraps="+test.constructions+
			" privateHelperCalls="+test.helperCalls+" metricSnapshots="+test.metricSnapshots+
			" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(root.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+root);}
		System.out.print(summary);
	}

	/** Resolves observation fields and the private helper without widening production APIs. */
	private AssemblyReaderCleanupTest(final Path root_, final boolean baseline_) throws Exception{
		root=root_; baseline=baseline_;
		length=field("length"); contigs=field("contigs"); bases=field("acgtnio");
		headerLength=field("headerLength"); firstHeader=field("firstHeader"); fname=field("fname");
		gc=Assembly.class.getDeclaredMethod("gc"); gc.setAccessible(true);
		load=Assembly.class.getDeclaredMethod("load", ByteFile.class);
		require(Modifier.isPrivate(load.getModifiers()) && load.getReturnType()==Void.TYPE,
			"This fixture requires the private void load(ByteFile) helper.");
		load.setAccessible(true);
		Files.write(root.resolve("prefix.fa"), PREFIX, StandardOpenOption.CREATE_NEW);
		Files.write(root.resolve("cases.tsv"), ("case\tcoverage\tphase\texpectedOutcome\tactualOutcome"+
			"\texpectedCloseCalls\tcloseCalls\tcloseReturns\ttrueReturns\tcloseThrows\tnextCalls\tlinesReturned"+
			"\treadErrors\tsamePrimary\tsameSecondary\tsuppressedCount\tsuppressedSameSecondary"+
			"\texceptionClass\texceptionMessage\tsuppressedClass\tsuppressedMessage\tlinesUnchanged\tglobalsRestored\tmetrics\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Acquires the sole real holder, always restoring the borrowed streams and reader globals. */
	private Assembly bootstrap() throws Exception{
		final InputStream originalIn=System.in;
		final PrintStream originalErr=System.err;
		final boolean originalForce=ByteFile.FORCE_MODE_BF1, originalSIMD=Shared.SIMD, originalBuffered=ByteFile1.BUFFERED;
		final BootstrapInput input=new BootstrapInput();
		final ByteArrayOutputStream diagnostics=new ByteArrayOutputStream();
		final PrintStream capture=new PrintStream(diagnostics, true, "UTF-8");
		final Assembly assembly;
		try{
			ByteFile.FORCE_MODE_BF1=true;
			Shared.SIMD=false;
			ByteFile1.BUFFERED=false;
			System.setIn(input);
			System.setErr(capture);
			constructions++;
			assembly=new Assembly(INPUT);
		}finally{
			System.setIn(originalIn);
			System.setErr(originalErr);
			ByteFile.FORCE_MODE_BF1=originalForce;
			Shared.SIMD=originalSIMD;
			ByteFile1.BUFFERED=originalBuffered;
			capture.close();
		}
		final boolean restored=globalsMatch(originalIn, originalErr, originalForce, originalSIMD, originalBuffered);
		check(input.closeCalls==0, "Bootstrap standard input is borrowed: "+input.closeCalls);
		check(input.available()==0, "Bootstrap did not consume its complete literal prefix: "+input.available());
		check(input.unchanged(), "Bootstrap prefix storage changed");
		check(restored, "Bootstrap did not restore streams and reader globals");
		check(diagnostics.size()==0, "Clean bootstrap emitted diagnostics: "+diagnostics.toString("UTF-8"));
		final String metrics=inspect(assembly, "public-bootstrap");
		Files.write(root.resolve("bootstrap.stderr"), diagnostics.toByteArray(), StandardOpenOption.CREATE_NEW);
		Files.write(root.resolve("bootstrap.tsv"), ("coverage\tpublicConstructions\tborrowedCloseCalls\tglobalsRestored\tprefixUnchanged\tmetrics\n"+
			"clean-public-constructor\t"+constructions+"\t"+input.closeCalls+"\t"+restored+"\t"+input.unchanged()+"\t"+metrics+"\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		return assembly;
	}

	/** Invokes real parsing on finite tracked lines and compares independent phase-specific oracles. */
	private void run(final Assembly holder, final int row) throws Exception{
		final boolean readError=row>=2 && row<=4, closeFlag=row==1 || row==3, closeError=row==4 || row==5;
		final String name=NAMES[row];
		final TrackedReader reader=new TrackedReader(name, readError, closeFlag, closeError);
		final InputStream originalIn=System.in;
		final PrintStream originalErr=System.err;
		final boolean originalForce=ByteFile.FORCE_MODE_BF1, originalSIMD=Shared.SIMD, originalBuffered=ByteFile1.BUFFERED;
		Throwable error=null;
		cases++; helperCalls++;
		try{load.invoke(holder, reader);}
		catch(final InvocationTargetException caught){error=caught.getCause();}
		final int expectedClose=baseline && readError ? 0 : 1;
		final int expectedSuppressed=!baseline && readError && (closeFlag || closeError) ? 1 : 0;
		final Throwable[] suppressed=error==null ? new Throwable[0] : error.getSuppressed();
		final String expected=readError ? "same-primary" : closeError ? "same-secondary" : closeFlag ? "input-error" : "return";
		final boolean outcome=readError ? error==reader.primary : closeError ? error==reader.secondary :
			closeFlag ? inputError(error) : error==null;
		final boolean unchanged=Arrays.equals(reader.header, new byte[]{'>', 'o', 'k'}) && Arrays.equals(reader.sequence, new byte[]{'A', 'C'});
		final boolean restored=globalsMatch(originalIn, originalErr, originalForce, originalSIMD, originalBuffered);
		check(reader.nextCalls==3, name+": expected three nextLine calls: "+reader.nextCalls);
		check(reader.linesReturned==2, name+": expected header and sequence lines: "+reader.linesReturned);
		check(reader.readErrors==(readError ? 1 : 0), name+": unexpected read exception count: "+reader.readErrors);
		check(reader.closeCalls==expectedClose, name+": expected close calls "+expectedClose+", saw "+reader.closeCalls);
		check(reader.closeReturns==(closeError ? 0 : expectedClose), name+": unexpected close returns: "+reader.closeReturns);
		check(reader.trueReturns==(closeFlag ? expectedClose : 0), name+": unexpected true close flags: "+reader.trueReturns);
		check(reader.closeThrows==(closeError ? expectedClose : 0), name+": unexpected thrown cleanup failures: "+reader.closeThrows);
		check(unchanged, name+": tracked line storage changed");
		check(restored, name+": helper changed bootstrap globals");
		check(outcome, name+": expected "+expected+", caught "+error);
		check(suppressed.length==expectedSuppressed, name+": expected suppressed count "+expectedSuppressed+", saw "+suppressed.length);
		check(expectedSuppressed==0 || (suppressed.length==1 && (closeError ? suppressed[0]==reader.secondary :
			inputError(suppressed[0]) && suppressed[0]!=reader.primary && suppressed[0]!=reader.secondary)),
			name+": secondary cleanup error lost identity or input diagnostic");
		final String metrics=error==null ? inspect(holder, name) : "unavailable";
		final Throwable secondary=suppressed.length==1 ? suppressed[0] : null;
		final String actual=error==null ? "return" : error==reader.primary ? "same-primary" : error==reader.secondary ? "same-secondary" :
			inputError(error) ? "input-error" : "other-error";
		final String observed=name+"\tprivate-load-ByteFile\t"+(baseline ? "baseline" : "fixed")+"\t"+expected+"\t"+actual+
			"\t"+expectedClose+"\t"+reader.closeCalls+"\t"+reader.closeReturns+"\t"+reader.trueReturns+"\t"+reader.closeThrows+
			"\t"+reader.nextCalls+"\t"+reader.linesReturned+"\t"+reader.readErrors+"\t"+(error==reader.primary)+"\t"+(error==reader.secondary)+
			"\t"+suppressed.length+"\t"+(secondary==reader.secondary)+"\t"+type(error)+"\t"+message(error)+
			"\t"+type(secondary)+"\t"+message(secondary)+"\t"+unchanged+"\t"+restored+"\t"+metrics+"\n";
		Files.write(root.resolve("cases.tsv"), observed.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
	}

	/** Checks literal metrics only after a successful public constructor or helper return. */
	private String inspect(final Assembly assembly, final String name) throws Exception{
		metricSnapshots++;
		final long size=length.getLong(assembly), headers=headerLength.getLong(assembly);
		final IntList list=(IntList)contigs.get(assembly);
		final int[] lengths=Arrays.copyOf(list.array, list.size);
		final int max=list.size==0 ? 0 : list.get(0);
		final long[] composition=((long[])bases.get(assembly)).clone();
		final float fraction=(Float)gc.invoke(assembly);
		final String first=(String)firstHeader.get(assembly), path=(String)fname.get(assembly);
		check(size==2, name+": expected two bases, saw "+size);
		check(list.size==1, name+": expected one contig, saw "+list.size);
		check(Arrays.equals(lengths, new int[]{2}), name+": expected lengths [2], saw "+Arrays.toString(lengths));
		check(max==2, name+": expected maximum two bases, saw "+max);
		check(Arrays.equals(composition, new long[]{1, 1, 0, 0, 0, 0, 0, 0}), name+": expected A1/C1, saw "+Arrays.toString(composition));
		check(Float.floatToIntBits(fraction)==Float.floatToIntBits(0.5f), name+": expected GC0.5, saw "+fraction);
		check(headers==2, name+": expected two header bytes, saw "+headers);
		check("ok".equals(first), name+": expected first header ok, saw "+first);
		check(INPUT.equals(path), name+": expected input stdin.fa, saw "+path);
		return "size="+size+";count="+list.size+";lengths="+Arrays.toString(lengths)+";max="+max+
			";bases="+Arrays.toString(composition)+";gc="+fraction+";headerBytes="+headers+";first="+first+";fname="+path;
	}

	/** Resolves an observation field without changing production visibility. */
	private static Field field(final String name) throws Exception{
		final Field field=Assembly.class.getDeclaredField(name);
		field.setAccessible(true);
		return field;
	}
	/** Requires the already accepted checked-close diagnostic, not an arbitrary exception. */
	private static boolean inputError(final Throwable error){
		return error!=null && error.getClass()==RuntimeException.class && String.valueOf(error.getMessage()).contains(INPUT);
	}
	/** Observes globals; helper calls do not intentionally change any of them. */
	private static boolean globalsMatch(final InputStream in, final PrintStream err, final boolean force, final boolean simd, final boolean buffered){
		return System.in==in && System.err==err && ByteFile.FORCE_MODE_BF1==force && Shared.SIMD==simd && ByteFile1.BUFFERED==buffered;
	}
	/** Makes one exception description fit in a TSV cell. */
	private static String message(final Throwable error){
		return error==null ? "" : String.valueOf(error.getMessage()).replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');
	}
	/** Records absent exceptions explicitly. */
	private static String type(final Throwable error){return error==null ? "none" : error.getClass().getName();}
	/** Accumulates independent failures without skipping the remaining rows. */
	private void check(final boolean correct, final String message){
		checks++;
		if(!correct){failures++; System.err.println("FAIL: "+message);}
	}
	/** Rejects fixture misuse before the public bootstrap changes process state. */
	private static void require(final boolean correct, final String message){
		if(!correct){throw new IllegalArgumentException(message);}
	}

	/** Finite borrowed bootstrap bytes; closing is observed separately from ByteFile cleanup. */
	private static final class BootstrapInput extends ByteArrayInputStream{
		BootstrapInput(){super(PREFIX.clone());}
		@Override
		public void close(){closeCalls++;}
		boolean unchanged(){return Arrays.equals(PREFIX, buf);}
		int closeCalls;
	}

	/** Supplies two immutable-oracle lines, then EOF or a fresh retained primary exception. */
	private static final class TrackedReader extends ByteFile{
		TrackedReader(final String name, final boolean readError_, final boolean closeFlag_, final boolean closeError_){
			super(FileFormat.testInput(INPUT, FileFormat.FASTA, null, false, false));
			readError=readError_; closeFlag=closeFlag_; closeError=closeError_;
			primary=new RuntimeException("Injected primary read failure: "+name);
			secondary=new RuntimeException("Injected secondary close failure: "+name);
		}
		@Override
		public byte[] nextLine(){
			nextCalls++;
			if(nextCalls==1){linesReturned++; return header;}
			if(nextCalls==2){linesReturned++; return sequence;}
			if(nextCalls>3){throw new AssertionError("Assembly exceeded the finite three-call reader contract: "+nextCalls);}
			if(readError){readErrors++; throw primary;}
			return null;
		}
		@Override
		public boolean close(){
			closeCalls++;
			if(closeError){closeThrows++; throw secondary;}
			closeReturns++;
			if(closeFlag){trueReturns++;}
			return closeFlag;
		}
		@Override
		public long lineNum(){return linesReturned;}
		@Override
		public boolean isOpen(){return closeCalls==0;}
		@Override
		public InputStream is(){throw new UnsupportedOperationException("Tracked reader owns no stream");}
		@Override
		public void reset(){throw new UnsupportedOperationException("Tracked reader is single use");}
		@Override
		public void pushBack(final byte[] line){throw new UnsupportedOperationException("Assembly must not push back fixture lines");}
		final byte[] header={'>', 'o', 'k'}, sequence={'A', 'C'};
		final RuntimeException primary, secondary;
		final boolean readError, closeFlag, closeError;
		int nextCalls, linesReturned, readErrors, closeCalls, closeReturns, trueReturns, closeThrows;
	}

	private final Path root;
	private final boolean baseline;
	private final Field length, contigs, bases, headerLength, firstHeader, fname;
	private final Method gc, load;
	private int cases, constructions, helperCalls, metricSnapshots, checks, failures;
	private static final String INPUT="stdin.fa";
	private static final byte[] PREFIX=">ok\nAC\n".getBytes(StandardCharsets.US_ASCII);
	private static final String[] NAMES={"normal-close-false", "normal-close-true", "primary-close-false",
		"primary-close-true", "primary-close-throws", "normal-close-throws"};
}
