package test.clump;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

import clump.KmerComparator;
import clump.StreamToOutput;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;

/** Finite startup checks and ordinary output controls for the overflow bridge.
 * The namesort baseline rejects the incomplete SortByName configuration. If a
 * replacement is reached, the probe records startup before allowing consumption.
 * All retained input/output handles are closed even when a processing call fails.
 * @author Shinobu
 */
public class StreamToOutputLifecycleTest {

	/** Runs baseline|fixed checks in a fresh directory on the selected test host. */
	public static void main(String[] args) throws Exception{
		assert(args.length==2) : "Expected baseline|fixed and a fresh directory";
		final boolean fixed=args[0].equals("fixed");
		assert(fixed || args[0].equals("baseline")) : "Unknown expectation: "+args[0];
		Shared.setThreads(1);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		FASTQ.ASCII_OFFSET=FASTQ.ASCII_OFFSET_OUT=33;
		FASTQ.DETECT_QUALITY=FASTQ.DETECT_QUALITY_OUT=false;
		final Path root=Paths.get(args[1]);
		Files.createDirectory(root);
		System.out.println("case\tinput_field_set\treplacement_started\toutput_matches");
		checkExisting(root, false, fixed);
		checkExisting(root, true, fixed);
		checkFileConstructor(root, fixed);
		if(fixed){checkPairedSort(root);}
	}

	/** Verifies exact records in ordinary and namesorted output. */
	private static void checkExisting(Path root, boolean sort, boolean fixed) throws Exception{
		final String name=sort ? "namesort" : "passthrough";
		final Path input=root.resolve(name+".fq"), output=root.resolve(name+".out.fq");
		Files.write(input, UNSORTED.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(
				-1, false, inputFormat(input), null, null, null);
		final ConcurrentReadOutputStream ros=outputStream(output);
		cris.start();
		ros.start();
		final StartupProbe probe=new StartupProbe(cris, ros, sort);
		boolean error=false, setupRejected=false;
		try{
			error=probe.process();
		}catch(RuntimeException e){
			if(!sort || fixed || !"Error - at least one input file is required.".equals(e.getMessage())){throw e;}
			setupRejected=true;
		}finally{
			// Retain this handle so the baseline's replaced or rejected input is closed.
			error|=ReadWrite.closeStream(cris);
			error|=ReadWrite.closeStream(ros);
		}
		assert(!error) : "Ordinary fixture I/O reported an error: "+name;
		assert(setupRejected==(sort && !fixed)) : "Unexpected SortByName constructor contract: "+name;
		if(setupRejected){
			assert(!probe.entered && Files.size(output)==0) : "Rejected setup must not emit any reads";
			System.out.println(name+"\ttrue\tNA\tfalse");
			return;
		}
		assert(probe.entered) : "Expected single-output dispatch: "+name;
		assert(probe.started) : "Input must start before nextList";
		final String expected=sort ? SORTED : UNSORTED;
		final String actual=new String(Files.readAllBytes(output), StandardCharsets.US_ASCII);
		assert(actual.equals(expected)) : "Record/base/quality mismatch for "+name+":\n"+actual;
		System.out.println(name+"\ttrue\t"+probe.started+"\ttrue");
	}

	/** An empty baseline source reaches EOF without a consumer. Never invokes
	 * process() with the baseline's missing input field. The fixed path drains it. */
	private static void checkFileConstructor(Path root, boolean fixed) throws Exception{
		final Path input=root.resolve("empty.fq"), output=root.resolve("empty.out.fq");
		Files.createFile(input);
		final ConcurrentReadOutputStream ros=outputStream(output);
		ros.start();
		final StreamToOutput tool=new StreamToOutput(inputFormat(input), null,
				new ConcurrentReadOutputStream[] {ros}, new KmerComparator(5, false, false), false, false);
		final Field field=StreamToOutput.class.getDeclaredField("cris");
		field.setAccessible(true);
		final ConcurrentReadInputStream actual=(ConcurrentReadInputStream)field.get(tool);
		assert((actual!=null)==fixed) : "File constructor must retain the stream it starts";
		if(actual!=null){assert(!tool.process()) : "Empty input processing failed";}
		assert(!ReadWrite.closeStream(ros)) : "Empty output close failed";
		assert(Files.size(output)==0) : "Empty input must not emit records";
		System.out.println("file-constructor\t"+(actual!=null)+"\tNA\ttrue");
	}

	private static FileFormat inputFormat(Path input){
		return FileFormat.testInput(input.toString(), FileFormat.FASTQ, null, false, true);
	}

	/** Fixed-path control: preserve both mates and their qualities through the temporary file. */
	private static void checkPairedSort(Path root) throws Exception{
		final String a1="@alpha/1\nACGTACGT\n+\nIIIIIIII\n";
		final String a2="@alpha/2\nTGCATGCA\n+\nGGGGGGGG\n";
		final String z1="@zeta/1\nTTGCAACC\n+\nHHHHHHHH\n";
		final String z2="@zeta/2\nGGCCAATT\n+\nFFFFFFFF\n";
		final Path input1=root.resolve("paired1.fq"), input2=root.resolve("paired2.fq");
		final Path output=root.resolve("paired.out.fq");
		Files.write(input1, (z1+a1).getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		Files.write(input2, (z2+a2).getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(
				-1, false, inputFormat(input1), inputFormat(input2), null, null);
		final ConcurrentReadOutputStream ros=outputStream(output);
		cris.start();
		ros.start();
		final StartupProbe probe=new StartupProbe(cris, ros, true);
		boolean error=false;
		try{error=probe.process();}finally{
			error|=ReadWrite.closeStream(cris);
			error|=ReadWrite.closeStream(ros);
		}
		assert(!error && probe.entered && probe.started) : "Paired name-sort lifecycle failed";
		final String actual=new String(Files.readAllBytes(output), StandardCharsets.US_ASCII);
		assert(actual.equals(a1+a2+z1+z2)) : "Name sorting lost pairing, bases or qualities:\n"+actual;
		System.out.println("paired-namesort\ttrue\ttrue\ttrue");
	}

	private static ConcurrentReadOutputStream outputStream(Path output){
		final FileFormat ff=FileFormat.testOutput(output.toString(), FileFormat.FASTQ,
				null, false, false, false, false);
		return ConcurrentReadOutputStream.getStream(ff, null, 4, null, false);
	}

	/** Observes startup without allowing the baseline to enter a blocking nextList. */
	private static class StartupProbe extends StreamToOutput{
		StartupProbe(ConcurrentReadInputStream cris, ConcurrentReadOutputStream ros, boolean sort){
			super(cris, new ConcurrentReadOutputStream[] {ros},
					new KmerComparator(5, false, false), sort, false);
		}
		@Override
		public void processSingle(ConcurrentReadInputStream input){
			entered=true;
			started=input.started();
			if(!started){input.start();}
			super.processSingle(input);
		}
		boolean entered, started;
	}

	private static final String A="@alpha\nACGTACGT\n+\nIIIIIIII\n";
	private static final String Z="@zeta\nTTGCAACC\n+\nHHHHHHHH\n";
	private static final String UNSORTED=Z+A;
	private static final String SORTED=A+Z;
}
