package test.ml;

import ml.CellNet;
import ml.Matrix;
import ml.ProcessBBMergeHeaders;
import ml.Sample;
import ml.SampleSet;
import ml.Trainer;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import fileIO.FileFormat;
import shared.Timer;

/**
 * Preserves ordinary ProcessBBMergeHeaders and Trainer helper output while
 * completion-result handling is tightened.
 * @author Sucrose
 */
public class CompletionIoSuccessTest{

	/** Runs one ordinary completion case. */
	public static void main(String[] args) throws Exception{
		assert(args.length==2) : "Usage: CompletionIoSuccessTest <case> <artifact-dir>";
		final String name=args[0];
		final Path root=Paths.get(args[1]);
		assert(Files.isDirectory(root)) : "Artifact directory must already exist: "+root;
		final Path dir=root.resolve(name);
		if(!(name.equals("process-bbmerge") || name.equals("trainer-network") ||
				name.equals("trainer-data"))){
			throw new IllegalArgumentException("Unknown completion I/O case: "+name);
		}
		Files.createDirectory(dir);
		if(name.equals("process-bbmerge")){
			processBBMerge(dir);
		}else if(name.equals("trainer-network")){
			trainerNetwork(dir);
		}else{
			trainerData(dir);
		}
		System.err.println("CompletionIoSuccessTest PASS "+name);
	}

	/** Runs ProcessBBMergeHeaders on one accepted BBMerge read and one rejected non-BBMerge read. */
	private static void processBBMerge(final Path dir) throws Exception{
		final Path in=dir.resolve("bbmerge.fq");
		final Path out=dir.resolve("bbmerge-processed.header");
		final String accepted="insert=202 synthetic mo=14_r1ee=5.2728_r2ee=3.4856"
			+"_bi=202_bo=98_bb=5.3063_br=0.0598_bbi=6_sbi=270_sbo=30"
			+"_sbb=12.4775_sbr=0.4343_sbbi=14_be=6.5990_pr=0.0007";
		Files.write(in, Arrays.asList(
			"@"+accepted,
			"ACGT",
			"+",
			"IIII",
			"@rejected",
			"TGCA",
			"+",
			"JJJJ"
		), StandardCharsets.US_ASCII);
		final String report=withErr(new ThrowingRunnable(){
			@Override
			public void run() throws Exception{
				new ProcessBBMergeHeaders(new String[]{"in="+in, "out="+out}).process(new Timer());
			}
		});
		Files.write(dir.resolve("report.stderr"), report.getBytes(StandardCharsets.US_ASCII));
		Files.write(dir.resolve("report-normalized.stderr"),
			normalizeProcessReport(report).getBytes(StandardCharsets.US_ASCII));

		final List<String> lines=Files.readAllLines(out, StandardCharsets.US_ASCII);
		final List<String> expected=Arrays.asList(
			BBMERGE_HEADER,
			"1\t14\t98\t6\t30\t14\t8.7584\t6.5990\t0.0598\t5.3063\t0.4343\t12.4775\t0.00070000"
		);
		assert(lines.equals(expected)) : "ProcessBBMergeHeaders header output drifted:\n"+
			"Expected:\n"+expected+"\nActual:\n"+lines;
		assert(report.contains("Reads Processed:    2 \t")) :
			"ProcessBBMergeHeaders must count accepted and rejected input reads exactly:\n"+report;
	}

	/** Directly exercises Trainer.writeNetwork with a deterministic tiny CellNet. */
	private static void trainerNetwork(final Path dir) throws Exception{
		final Path data=writeTrainerInput(dir);
		final Path out=dir.resolve("network.bbnet");
		final Trainer trainer=new Trainer(new String[]{"in="+data, "dims=1,1", "epochs=1"});
		final ArrayList<String> commands=new ArrayList<String>();
		commands.add("#CL CompletionIoSuccessTest trainer-network");
		final CellNet net=Trainer.randomNetwork(new int[]{1, 1}, 1, 1, 0, -1, 1, commands);
		invoke(trainer, "writeNetwork", net, out);

		final String text=read(out);
		assert(text.startsWith("##bbnet\n")) : "Trainer.writeNetwork must emit a bbnet header:\n"+text;
		assert(text.contains("#dims 1 1\n")) :
			"Trainer.writeNetwork must preserve the one-input, one-output topology:\n"+text;
		assert(text.contains("#CL CompletionIoSuccessTest trainer-network\n")) :
			"Trainer.writeNetwork must preserve command metadata:\n"+text;
		assert(((Long)getField(trainer, "linesOut")).longValue()==net.lastLinesWritten) :
			"writeNetwork must keep the legacy network line counter in sync";
		assert(((Long)getField(trainer, "bytesOut")).longValue()>0) :
			"writeNetwork must count serialized network bytes";
	}

	/** Directly exercises Trainer.writeData with deterministic materialized samples. */
	private static void trainerData(final Path dir) throws Exception{
		final Path data=writeTrainerInput(dir);
		final Path out=dir.resolve("vectors.tsv");
		final Trainer trainer=new Trainer(new String[]{"in="+data, "dims=1,1", "epochs=1"});
		setField(trainer, "data", makeSampleSet());
		invoke(trainer, "writeData", out);

		final String expected=
			"dims\t1\t1\n"+
			"feature\ttarget\n"+
			"0\t0\n"+
			"1\t1\n";
		final String actual=read(out);
		assert(actual.equals(expected)) : "Trainer.writeData output drifted:\n"+
			"Expected:\n"+expected+"Actual:\n"+actual;
		assert(((Long)getField(trainer, "linesOut")).longValue()==2) :
			"writeData must count two materialized samples";
		assert(((Long)getField(trainer, "bytesOut")).longValue()==expected.length()) :
			"writeData must count the header and first sample together, then later samples";
	}

	/** Invokes one private Trainer writer method with an explicit FileFormat. */
	private static void invoke(final Trainer trainer, final String name, final Path out)
			throws Exception{
		invoke(trainer, name, FileFormat.testOutput(out.toString(),
			FileFormat.TXT, null, true, true, false, false));
	}

	/** Invokes one private Trainer writer method with a network and explicit FileFormat. */
	private static void invoke(final Trainer trainer, final String name, final CellNet net,
			final Path out) throws Exception{
		invoke(trainer, name, net, FileFormat.testOutput(out.toString(),
			FileFormat.TXT, null, true, true, false, false));
	}

	/** Invokes one private Trainer writer method by reflection. */
	private static void invoke(final Trainer trainer, final String name, final Object... args)
			throws Exception{
		final Class<?>[] types=new Class<?>[args.length];
		for(int i=0; i<args.length; i++){types[i]=args[i].getClass();}
		final Method method=Trainer.class.getDeclaredMethod(name, types);
		method.setAccessible(true);
		method.invoke(trainer, args);
	}

	/** Writes a tiny training matrix that satisfies Trainer construction. */
	private static Path writeTrainerInput(final Path dir) throws Exception{
		final Path path=dir.resolve("trainer-input.tsv");
		Files.write(path, Arrays.asList(
			"#dims\t1\t1",
			"0\t0",
			"1\t1"
		), StandardCharsets.US_ASCII);
		return path;
	}

	/** Builds the exact rows used by the direct writeData case without shuffling them. */
	private static SampleSet makeSampleSet(){
		final Matrix m=new Matrix();
		m.columns=new ArrayList<String>();
		m.columns.add("feature");
		m.columns.add("target");
		m.numInputs=1;
		m.numOutputs=1;
		m.inputs=new float[][]{{0}, {1}};
		m.outputs=new float[][]{{0}, {1}};
		m.weights=new float[][]{{1}, {1}};
		m.data=new float[][][]{m.inputs, m.outputs, m.weights};
		m.validLines=m.inputs.length;
		m.initializeRange();
		final SampleSet set=new SampleSet(m);
		set.samples=new Sample[]{
			new Sample(m.inputs[0], m.outputs[0], 1, 0),
			new Sample(m.inputs[1], m.outputs[1], 1, 1)
		};
		return set;
	}

	/** Removes only elapsed time and rate variability from the ProcessBBMergeHeaders report. */
	private static String normalizeProcessReport(final String report){
		final StringBuilder sb=new StringBuilder();
		for(String line : report.split("\n", -1)){
			if(line.startsWith("Time:")){
				line="Time:                         \t<elapsed>";
			}else if(line.startsWith("Reads Processed:")){
				final int tab=line.indexOf('\t');
				assert(tab>=0) : "Expected tab before read-rate field: "+line;
				line=line.substring(0, tab+1)+"<rate> reads/sec";
			}
			if(sb.length()>0){sb.append('\n');}
			sb.append(line);
		}
		return sb.toString();
	}

	/** Reads a US-ASCII text file. */
	private static String read(final Path path) throws Exception{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}

	/** Assigns one private field by reflection. */
	private static void setField(final Object owner, final String name, final Object value)
			throws Exception{
		final Field field=owner.getClass().getDeclaredField(name);
		field.setAccessible(true);
		field.set(owner, value);
	}

	/** Returns one private field by reflection. */
	private static Object getField(final Object owner, final String name) throws Exception{
		final Field field=owner.getClass().getDeclaredField(name);
		field.setAccessible(true);
		return field.get(owner);
	}

	/** Captures stderr while a tool writes its report. */
	private static String withErr(final ThrowingRunnable runnable) throws Exception{
		final PrintStream old=System.err;
		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final PrintStream ps=new PrintStream(buffer, true, "UTF-8");
		try{
			System.setErr(ps);
			runnable.run();
		}finally{
			System.setErr(old);
			ps.close();
		}
		return new String(buffer.toByteArray(), StandardCharsets.UTF_8);
	}

	/** Runnable that can throw checked exceptions. */
	private static interface ThrowingRunnable{
		void run() throws Exception;
	}

	private static final String BBMERGE_HEADER=
		"#Correct\tminOverlap\tbestOverlap\tbestBadInt\tsecondBestOverlap\tsecondBestBadInt\t"+
		"expectedErrors\tbestExpectedErrors\tbestRatio\tbestBad\tsecondBestRatio\tsecondBestBad\t"+
		"probability";
}
