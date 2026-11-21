package test.jgi;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
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
 * Two real Assembly constructor calls over borrowed raw stdin: clean EOF and a
 * ByteFile1-latched read IOException after the same seven-byte FASTA prefix.
 * Baseline accepts the partial assembly; fixed must reject it with the input name.
 * The clean two-base control retains its exact metrics. Reader selection, scalar
 * parsing, raw buffering and System streams are restored after each call.
 * No private reader helper, decoder, malformed FASTA, owned file, worker or timeout
 * is used. Borrowed-stream close counts do not establish ByteFile.close invocation
 * or owned-reader cleanup. Original IOException identity is not an oracle because
 * ByteFile1 catches it before Assembly receives the close-error flag.
 * @author Sucrose
 */
public final class AssemblyReadErrorTest{

	/** Requires phase=baseline|fixed, a fresh outdir and a one-CPU fixture JVM. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown reader fixture option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a new outdir.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Use dev/assemblyreaderrortest.sh with assertions enabled.");
		require(Shared.threads()==1, "The raw scalar fixture requires one tool thread.");
		require(!new File(INPUT).exists(), "A real stdin.fa would override borrowed standard input.");
		require(ByteFile1.bufferlen>=PREFIX.length, "The seven-byte prefix must fit in the initial reader buffer.");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		Files.createDirectory(root);
		final AssemblyReadErrorTest test=new AssemblyReadErrorTest(root, phase.equals("baseline"));
		test.run(false);
		test.run(true);
		test.check(test.cases==2, "Expected the clean and read-error rows: "+test.cases);
		test.check(test.constructions==2, "Call the public constructor once per row: "+test.constructions);
		test.check(test.partialAccepted==(test.baseline ? 1 : 0), "Unexpected accepted partial assemblies: "+test.partialAccepted);
		test.check(test.rejected==(test.baseline ? 0 : 1), "Unexpected read-error rejections: "+test.rejected);
		test.check(test.metricSnapshots==(test.baseline ? 2 : 1), "Metrics exist only for returned assemblies: "+test.metricSnapshots);
		test.check(test.checks==(test.baseline ? 39 : 30), "Unexpected pre-count check total: "+test.checks);
		final String summary="ASSEMBLY_READ_ERROR_"+(test.failures==0 ? "PASS" : "FAIL")+
			" phase="+phase+" cases="+test.cases+" constructions="+test.constructions+
			" partialAccepted="+test.partialAccepted+" rejected="+test.rejected+
			" metricSnapshots="+test.metricSnapshots+" checks="+test.checks+" failures="+test.failures+"\n";
		Files.write(root.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		if(test.failures>0){throw new AssertionError(summary+"Inspect "+root);}
		System.out.print(summary);
	}

	/** Binds metric observations without changing production field or method visibility. */
	private AssemblyReadErrorTest(final Path root_, final boolean baseline_) throws Exception{
		root=root_; baseline=baseline_;
		length=field("length"); contigs=field("contigs"); bases=field("acgtnio");
		headerLength=field("headerLength"); firstHeader=field("firstHeader"); fname=field("fname");
		gc=Assembly.class.getDeclaredMethod("gc"); gc.setAccessible(true);
		Files.write(root.resolve("prefix.fa"), PREFIX, StandardOpenOption.CREATE_NEW);
		Files.write(root.resolve("cases.tsv"), ("case\tphase\texpectedOutcome\tactualOutcome\tbulkCalls\tsingleCalls"+
			"\tbytesDelivered\tinjectedErrors\tborrowedCloseCalls\tglobalsRestored\tprefixUnchanged"+
			"\texceptionClass\texceptionMessage\tmetrics\n").getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Uses the public acquisition path; only the borrowed stream's finite reads are injected. */
	private void run(final boolean inject) throws Exception{
		final String name=inject ? "read-error" : "clean-eof";
		final InputStream originalIn=System.in;
		final PrintStream originalErr=System.err;
		final boolean originalForce=ByteFile.FORCE_MODE_BF1, originalSIMD=Shared.SIMD, originalBuffered=ByteFile1.BUFFERED;
		final TrackedInput input=new TrackedInput(inject);
		final ByteArrayOutputStream diagnostics=new ByteArrayOutputStream();
		final PrintStream capture=new PrintStream(diagnostics, true, "UTF-8");
		Assembly assembly=null;
		Throwable error=null;
		cases++;
		try{
			ByteFile.FORCE_MODE_BF1=true;
			Shared.SIMD=false;
			ByteFile1.BUFFERED=false;
			System.setIn(input);
			System.setErr(capture);
			constructions++;
			try{assembly=new Assembly(INPUT);}
			catch(final Throwable caught){error=caught;}
		}finally{
			System.setIn(originalIn);
			System.setErr(originalErr);
			ByteFile.FORCE_MODE_BF1=originalForce;
			Shared.SIMD=originalSIMD;
			ByteFile1.BUFFERED=originalBuffered;
			capture.close();
		}
		final boolean restored=System.in==originalIn && System.err==originalErr &&
			ByteFile.FORCE_MODE_BF1==originalForce && Shared.SIMD==originalSIMD && ByteFile1.BUFFERED==originalBuffered;
		final boolean unchanged=Arrays.equals(PREFIX, input.bytes);
		final String diagnostic=new String(diagnostics.toByteArray(), StandardCharsets.UTF_8);
		Files.write(root.resolve(name+".stderr"), diagnostics.toByteArray(), StandardOpenOption.CREATE_NEW);
		final boolean expectFailure=inject && !baseline;
		if(inject && assembly!=null){partialAccepted++;}
		if(inject && error!=null){rejected++;}
		check(expectFailure ? assembly==null && error instanceof RuntimeException &&
			String.valueOf(error.getMessage()).contains(INPUT) : assembly!=null && error==null,
			name+": unexpected constructor outcome: "+error);
		check(input.bulkCalls==2 && input.singleCalls==0, name+": expected two bulk reads and no single reads; saw "+input.bulkCalls+"/"+input.singleCalls);
		check(input.delivered==7, name+": production must consume the complete seven-byte literal prefix: "+input.delivered);
		check(input.errors==(inject ? 1 : 0), name+": unexpected injected read-error count: "+input.errors);
		check(input.closeCalls==0, name+": exact System.in is borrowed and must stay open: "+input.closeCalls);
		check(unchanged, name+": input prefix storage changed");
		check(restored, name+": System streams or reader settings were not restored");
		check(diagnostic.equals(inject ? WARNING+System.lineSeparator() : ""), name+": unexpected reader diagnostics: "+diagnostic);
		final String metrics=assembly==null ? "unavailable" : inspect(assembly, name);
		final String observed=name+"\t"+(baseline ? "baseline" : "fixed")+"\t"+(expectFailure ? "reject" : "return")+
			"\t"+(assembly==null ? "reject" : "return")+"\t"+input.bulkCalls+"\t"+input.singleCalls+
			"\t"+input.delivered+"\t"+input.errors+"\t"+input.closeCalls+"\t"+restored+"\t"+unchanged+
			"\t"+(error==null ? "none" : error.getClass().getName())+"\t"+clean(error==null ? "" : String.valueOf(error.getMessage()))+
			"\t"+metrics+"\n";
		Files.write(root.resolve("cases.tsv"), observed.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
	}

	/** Checks nine independent literal metrics only when the public constructor returns. */
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
		check(Arrays.equals(lengths, new int[]{2}), name+": expected descending lengths [2], saw "+Arrays.toString(lengths));
		check(max==2, name+": expected maximum two bases, saw "+max);
		check(Arrays.equals(composition, new long[]{1, 1, 0, 0, 0, 0, 0, 0}), name+": expected A1/C1, saw "+Arrays.toString(composition));
		check(Float.floatToIntBits(fraction)==Float.floatToIntBits(0.5f), name+": expected GC0.5, saw "+fraction);
		check(headers==2, name+": expected two header bytes, saw "+headers);
		check("ok".equals(first), name+": expected first header ok, saw "+first);
		check(INPUT.equals(path), name+": expected the literal stdin.fa path, saw "+path);
		return "size="+size+";count="+list.size+";lengths="+Arrays.toString(lengths)+";max="+max+
			";bases="+Arrays.toString(composition)+";gc="+fraction+";headerBytes="+headers+";first="+first+";fname="+path;
	}

	/** Resolves an observation field, never a reader factory or private load helper. */
	private static Field field(final String name) throws Exception{
		final Field field=Assembly.class.getDeclaredField(name);
		field.setAccessible(true);
		return field;
	}
	/** Keeps one exception message in a TSV cell. */
	private static String clean(final String value){return value.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');}
	/** Records failures without hiding the other row's outcome. */
	private void check(final boolean correct, final String message){
		checks++;
		if(!correct){failures++; System.err.println("FAIL: "+message);}
	}
	/** Rejects fixture misuse before changing process-wide state. */
	private static void require(final boolean correct, final String message){
		if(!correct){throw new IllegalArgumentException(message);}
	}

	/** Supplies a finite prefix, then real EOF or an IOException from the actual bulk read. */
	private static final class TrackedInput extends InputStream{
		TrackedInput(final boolean inject_){inject=inject_;}
		@Override
		public int read() throws IOException{
			singleCalls++;
			throw new IOException("Unexpected single-byte access to the BF1 fixture");
		}
		@Override
		public int read(final byte[] target, final int offset, final int length) throws IOException{
			if(target==null){throw new NullPointerException();}
			if(offset<0 || length<0 || offset>target.length-length){throw new IndexOutOfBoundsException();}
			if(length==0){return 0;}
			bulkCalls++;
			if(delivered<bytes.length){
				final int count=Math.min(length, bytes.length-delivered);
				System.arraycopy(bytes, delivered, target, offset, count);
				delivered+=count;
				return count;
			}
			if(inject){errors++; throw new IOException("Injected Assembly BF1 read failure");}
			return -1;
		}
		/** Observes borrowed-input ownership; it is not evidence of ByteFile closure. */
		@Override
		public void close(){closeCalls++;}
		final byte[] bytes=PREFIX.clone();
		final boolean inject;
		int bulkCalls, singleCalls, delivered, errors, closeCalls;
	}

	private final Path root;
	private final boolean baseline;
	private final Field length, contigs, bases, headerLength, firstHeader, fname;
	private final Method gc;
	private int cases, constructions, partialAccepted, rejected, metricSnapshots, checks, failures;
	private static final String INPUT="stdin.fa";
	private static final byte[] PREFIX=">ok\nAC\n".getBytes(StandardCharsets.US_ASCII);
	private static final String WARNING="Error: truncated or corrupt input for 'stdin.fa'; data may be incomplete.";
}
