package test.jgi;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import fileIO.ReadWrite;
import jgi.BBDuk;

/** Finite retained-PrintStream-error observations, one scenario per fresh JVM.
 * Direct cases exercise close or finishWriting; program cases use ordinary BBDuk
 * construction and public process(), including its legacy output writer and join.
 * The only fault injection is a bounded in-memory OutputStream throwing IOException.
 * Borrowed stdout is restored in finally. Reflection reads only final input/output
 * counters; no private fields, flags, reference tables or error states are seeded.
 * Baseline expects the currently swallowed errors. Fixed expects their propagation.
 * This fixture does not exercise broken pipes, subprocesses or compression codecs.
 * @author Shinobu
 */
public final class BBDukStdoutErrorTest{

	/** Requires scenario=NAME phase=baseline|fixed outdir=NEW through the dev launcher. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final int split=arg.indexOf('=');
			require(split>0 && split<arg.length()-1, "Use nonempty flag=value arguments: "+arg);
			require(options.put(arg.substring(0, split), arg.substring(split+1))==null, "Duplicate argument: "+arg);
		}
		require(options.size()==3 && options.keySet().containsAll(Arrays.asList("scenario", "phase", "outdir")), "Use exactly scenario, phase and outdir.");
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Run the dev launcher with assertions enabled.");
		final String name=options.get("scenario"), phase=options.get("phase");
		require(Arrays.asList(NAMES).contains(name), "Unknown scenario: "+name);
		require(phase.equals("baseline") || phase.equals("fixed"), "Use phase=baseline or phase=fixed.");
		final Path out=Paths.get(options.get("outdir")).toAbsolutePath().normalize();
		require(out.toString().indexOf('\t')<0 && out.toString().indexOf('\n')<0, "Artifact paths cannot contain tabs or newlines.");
		new BBDukStdoutErrorTest(name, phase, out).run();
	}

	private BBDukStdoutErrorTest(final String name_, final String phase_, final Path out_) throws IOException{
		name=name_; phase=phase_; out=out_;
		program=name.startsWith("bbduk_"); file=name.equals("bbduk_file_ok");
		borrowed=!file && (program || name.startsWith("borrowed_"));
		failWrite=name.endsWith("write_error"); failFlush=name.endsWith("flush_error"); failClose=name.endsWith("close_error");
		fault=failWrite || failFlush || failClose;
		Files.createDirectory(out);
		write("contract.txt", "One fresh JVM; finite in-memory write/flush/close failure; preserve borrowed stdout ownership.\n"+
			"Direct owned uses ReadWrite.close; direct borrowed uses finishWriting(null,...); BBDuk uses public process().\n"+
			"Baseline: retained errors are observable but the finalizer/program does not report them.\n"+
			"Fixed: the finalizer returns true or BBDuk sets errorState and throws its normal terminal error.\n");
	}

	private void run() throws Exception{
		final PrintStream saved=System.out;
		final ProbeOutputStream backend=(file ? null : new ProbeOutputStream(failWrite, failFlush, failClose));
		final PrintStream sink=(file ? null : new PrintStream(backend, false, "US-ASCII"));
		boolean reported=false, retained=false, ownerError=false, ownerThrew=false, reusable=false;
		int writeCalls=0, productionFlushes=0, productionCloses=0;
		byte[] payload=null, after=null;
		try{
			if(borrowed){System.setOut(sink);}
			if(program){
				final Path input=out.resolve("input.fq"), output=out.resolve("processed.fq");
				Files.write(input, FASTQ, StandardOpenOption.CREATE_NEW);
				final String[] args={"in="+input, "out="+(file ? output.toString() : "stdout.fq"),
					"t=1", "simd=f", "qin=33", "qout=33", "minlen=1", "qtrim=f", "maq=0", "maxns=-1",
					"reads=1", "array=t", "prealloc=f", "ow=f", "app=f"};
				write("constructor-args.txt", String.join("\n", args)+"\n");
				final BBDuk owner=new BBDuk(args);
				check(!owner.errorState, "An ordinary fresh BBDuk owner must not begin in an error state.");
				try{owner.process();}
				catch(RuntimeException e){
					check(e.getClass()==RuntimeException.class && (BBDuk.class.getName()+" terminated in an error state; the output may be corrupt.").equals(e.getMessage()),
						"Only BBDuk's final errorState exception is an expected propagation result: "+e);
					ownerThrew=true;
					write("owner-exception.txt", e.getClass().getName()+"\n"+e.getMessage()+"\n");
				}
				ownerError=owner.errorState;
				check(ownerThrew==ownerError, "BBDuk.process must throw when its final errorState is set, and otherwise return normally.");
				reported=ownerThrew && ownerError;
				final long readsIn=counter(owner, "readsIn"), basesIn=counter(owner, "basesIn"), readsOut=counter(owner, "readsOut"), basesOut=counter(owner, "basesOut");
				check(readsIn==1 && basesIn==8 && readsOut==1 && basesOut==8, "The same unfiltered one-read path must reach output aggregation even when the sink loses its bytes.");
				write("owner-counters.tsv", "reads_in\tbases_in\treads_out\tbases_out\n"+readsIn+"\t"+basesIn+"\t"+readsOut+"\t"+basesOut+"\n");
				if(file){payload=Files.readAllBytes(output);}
			}else{
				sink.write(DIRECT);
				reported=(borrowed ? ReadWrite.finishWriting(null, sink, "stdout.fq", false) : ReadWrite.close(sink));
			}
			if(sink!=null){
				writeCalls=backend.writeCalls; productionFlushes=backend.flushCalls; productionCloses=backend.closeCalls;
				retained=sink.checkError();
				payload=backend.bytes.toByteArray();
				check(writeCalls>0, "The probe must observe a real write attempt, not a seeded PrintStream error flag.");
				check(productionFlushes>0, "Finalization must flush the sink before checkError observation can itself flush it.");
				check(productionCloses==(borrowed ? 0 : 1), "Finalization must attempt the owned close once, including a failing close, and never close borrowed stdout; attempts="+productionCloses);
				check(retained==fault, "PrintStream must retain the injected write/flush/close error and stay clear for controls.");
				if(borrowed && !fault){
					sink.write(AFTER); sink.flush();
					after=backend.bytes.toByteArray();
					final byte[] expected=Arrays.copyOf(payload, payload.length+AFTER.length);
					System.arraycopy(AFTER, 0, expected, payload.length, AFTER.length);
					reusable=Arrays.equals(after, expected) && !sink.checkError();
					check(reusable, "Healthy borrowed stdout must still accept bytes after finalization.");
				}
			}
			final byte[] expected=(failWrite ? new byte[0] : program ? FASTQ : DIRECT);
			check(Arrays.equals(payload, expected), "Literal data differs: write faults retain no bytes; all other cases preserve the complete literal output.");
			check(reported==(fault && phase.equals("fixed")), "Baseline/fixed finalization must expose the stated retained-error contract; scenario="+name+" reported="+reported);
		}finally{
			if(borrowed){System.setOut(saved);}
			// The fixture owns this in-memory sink after observing the production ownership contract.
			if(sink!=null){sink.close();}
		}
		Files.write(out.resolve("payload.bin"), payload, StandardOpenOption.CREATE_NEW);
		if(after!=null){Files.write(out.resolve("after-finalization.bin"), after, StandardOpenOption.CREATE_NEW);}
		final int discrepancy=(fault && !reported ? 1 : 0);
		write("observations.tsv", "scenario\tphase\tprogram\tfile\tborrowed\tfault\treported\tretained\treusable\towner_error\towner_threw\twrite_calls\tproduction_flushes\tproduction_close_attempts\tbytes\tdiscrepancy\n"+
			name+"\t"+phase+"\t"+program+"\t"+file+"\t"+borrowed+"\t"+fault+"\t"+reported+"\t"+retained+"\t"+reusable+"\t"+ownerError+"\t"+ownerThrew+"\t"+writeCalls+"\t"+productionFlushes+"\t"+productionCloses+"\t"+payload.length+"\t"+discrepancy+"\n");
		final String summary="BBDUK_STDOUT_ERROR_PASS scenario="+name+" phase="+phase+" fault="+fault+" reported="+reported+" discrepancies="+discrepancy+" checks="+checks+"\n";
		write("summary.txt", summary); System.out.print(summary);
	}

	/** Reads normal owner aggregation; never alters the owner's state. */
	private static long counter(final BBDuk owner, final String name) throws Exception{
		final Field field=BBDuk.class.getDeclaredField(name); field.setAccessible(true);
		return ((Number)field.get(owner)).longValue();
	}
	private void check(final boolean value, final String why){checks++; require(value, why);}
	private static void require(final boolean value, final String why){if(!value){throw new IllegalArgumentException(why);}}
	private void write(final String name, final String value) throws IOException{
		Files.write(out.resolve(name), value.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Deliberately finite sink; faults throw once per attempted operation, never block. */
	private static final class ProbeOutputStream extends OutputStream{
		ProbeOutputStream(final boolean write_, final boolean flush_, final boolean close_){
			assert((write_ ? 1 : 0)+(flush_ ? 1 : 0)+(close_ ? 1 : 0)<=1) : "Each finite scenario isolates one write, flush or close fault.";
			failWrite=write_; failFlush=flush_; failClose=close_;
		}
		@Override
		public void write(final int value) throws IOException{
			assert(closeCalls==0) : "PrintStream must not delegate data after attempting its owned backend close.";
			writeCalls++;
			if(failWrite){throw new IOException("EXPECTED_WRITE_FAILURE");}
			bytes.write(value);
		}
		@Override
		public void write(final byte[] data, final int offset, final int length) throws IOException{
			assert(closeCalls==0 && data!=null && offset>=0 && length>=0 && offset<=data.length-length) : "PrintStream must delegate valid ranges before attempting backend close.";
			writeCalls++;
			if(failWrite){throw new IOException("EXPECTED_WRITE_FAILURE");}
			bytes.write(data, offset, length);
		}
		@Override
		public void flush() throws IOException{
			assert(closeCalls==0) : "PrintStream must finish flushing before attempting its owned backend close.";
			flushCalls++;
			if(failFlush){throw new IOException("EXPECTED_FLUSH_FAILURE");}
		}
		@Override
		public void close() throws IOException{
			closeCalls++;//An attempted close is recorded even when it fails.
			if(failClose){throw new IOException("EXPECTED_CLOSE_FAILURE");}
		}
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		final boolean failWrite, failFlush, failClose;
		int writeCalls, flushCalls, closeCalls;
	}

	private final String name, phase;
	private final Path out;
	private final boolean program, file, borrowed, failWrite, failFlush, failClose, fault;
	private int checks;
	private static final byte[] DIRECT="ABC\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] AFTER="AFTER\n".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] FASTQ="@stdout_probe\nACGTACGT\n+\nIIIIIIII\n".getBytes(StandardCharsets.US_ASCII);
	private static final String[] NAMES={"owned_ok", "owned_close_error", "borrowed_ok", "borrowed_write_error", "borrowed_flush_error", "bbduk_stdout_ok", "bbduk_stdout_write_error", "bbduk_file_ok"};
}
