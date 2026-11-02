package test.scalar;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import fileIO.FileFormat;
import scalar.ScalarData;
import scalar.ScalarInterval;
import scalar.ScalarIntervals;
import structures.IntList;

/** Finite oracle for ScalarData generated Start/End coordinate contracts. */
public class ScalarWindowCoordinateTest{

	/** Runs all coordinate contract cases, or one case with {@code case=}. */
	public static void main(String[] args) throws Exception{
		String mode="all";
		Path dir=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("work=")){dir=java.nio.file.Paths.get(arg.substring(5));}
			else if(dir==null){dir=java.nio.file.Paths.get(arg);}
			else{throw new IllegalArgumentException(arg);}
		}
		if(dir==null){throw new IllegalArgumentException("Missing work=<freshdir>");}
		Files.createDirectories(dir);
		if(mode.equals("all")){
			run("rolling", dir.resolve("rolling"));
			run("finite-edges", dir.resolve("finite-edges"));
			run("reset", dir.resolve("reset"));
			run("unknown", dir.resolve("unknown"));
			run("tsv", dir.resolve("tsv"));
			run("merge-roundtrip", dir.resolve("merge-roundtrip"));
		}else{run(mode, dir.resolve(mode));}
		System.out.println("PASS "+mode);
	}

	/** Runs one independent coordinate contract case in a reset fixture directory. */
	private static void run(String mode, Path dir) throws Exception{
		Files.createDirectories(dir);
		resetStatics();
		if(mode.equals("rolling")){
			final ScalarData d=intervals(dir.resolve("rolling.fa"), ">one\nACGTACGT\n", 4, 2, 1, true);
			writeActual(dir.resolve("rolling.actual.tsv"), d);
			requireRows(d, 3);
			requireRow(d, 0, 5, 0, 5);
			requireRow(d, 1, 2, 2, 7);
			requireRow(d, 2, 1, 3, 8);
		}else if(mode.equals("finite-edges")){
			ScalarData d=intervals(dir.resolve("partial.fa"), ">partial\nACG\n", 4, 100, 1, true);
			writeActual(dir.resolve("partial.actual.tsv"), d);
			requireRows(d, 1);
			requireRow(d, 0, 3, 0, 3);

			d=intervals(dir.resolve("with_ns.fa"), ">with_ns\nACGTNNACGTNN\n", 4, 100, 1, true);
			writeActual(dir.resolve("with_ns.actual.tsv"), d);
			requireRows(d, 1);
			requireRow(d, 0, 8, 2, 10);

			d=intervals(dir.resolve("cross_ambiguity.fa"),
					">cross_ambiguity\nACGTNNACGTNNACGT\n", 4, 100, 1, true);
			writeActual(dir.resolve("cross_ambiguity.actual.tsv"), d);
			requireRows(d, 1);
			requireRow(d, 0, 12, 8, 16);

			d=intervals(dir.resolve("isolated.fa"), ">isolated\nAN\n", 4, 100, 1, true);
			writeActual(dir.resolve("isolated.actual.tsv"), d);
			requireRows(d, 1);
			requireRow(d, 0, 1, -1, -1);

			d=intervals(dir.resolve("alln.fa"), ">alln\nNNNN\n", 4, 100, 1, true);
			writeActual(dir.resolve("alln.actual.tsv"), d);
			requireRows(d, 0);
		}else if(mode.equals("reset")){
			final ScalarData d=intervals(dir.resolve("reset.fa"), ">a\nACGT\n>b\nACGT\n", 4, 100, 1, true);
			writeActual(dir.resolve("reset.actual.tsv"), d);
			requireRows(d, 2);
			requireRow(d, 0, 4, 0, 4);
			requireRow(d, 1, 4, 0, 4);
		}else if(mode.equals("unknown")){
			ScalarData d=intervals(dir.resolve("break_false.fa"), ">a\nACGT\n>b\nACGT\n", 4, 2, 1, false);
			writeActual(dir.resolve("break_false.actual.tsv"), d);
			require(d.gc.size()>0, "break=false fixture did not emit rows");
			requireUnknown(d);

			d=intervals(dir.resolve("window_zero.fa"), ">one\nACGT\n", 0, 100, 1, true);
			writeActual(dir.resolve("window_zero.actual.tsv"), d);
			requireRows(d, 1);
			requireRow(d, 0, 4, 0, 4);
		}else if(mode.equals("tsv")){
			final Path legacy=dir.resolve("legacy_start.tsv");
			write(legacy, "#Name\tLength\tGC\tHH\tCAGA\tDepth\tStart\tTaxID\tTaxID2\n"+
					"ctg\t6\t0.5\t0.25\t0.75\t3\t12\t7\t8\n");
			ScalarData d=data(true).readTSV(input(legacy));
			writeActual(dir.resolve("legacy_start.actual.tsv"), d);
			requireRows(d, 1);
			requireRow(d, 0, 6, 12, -1);

			final Path explicit=dir.resolve("explicit_end.tsv");
			write(explicit, "#Name\tLength\tGC\tHH\tCAGA\tDepth\tStart\tTaxID\tTaxID2\tEnd\n"+
					"ctg\t4\t0.5\t0.25\t0.75\t2\t2\t11\t12\t9\n");
			d=data(true).readTSV(input(explicit));
			writeActual(dir.resolve("explicit_end.actual.tsv"), d);
			requireRows(d, 1);
			requireRow(d, 0, 4, 2, 9);
			final ScalarInterval row=d.getRow(0, null);
			require(row.taxID==11 && row.taxID2==12,
					"End column shifted taxonomy fields");
		}else if(mode.equals("merge-roundtrip")){
			final Path a=dir.resolve("a.tsv"), b=dir.resolve("b.tsv");
			write(a, "#Name\tLength\tGC\tHH\tCAGA\tDepth\tStart\tTaxID\tTaxID2\tEnd\n"+
					"a\t3\t0.25\t0.25\t0.75\t1\t4\t5\t6\t7\n");
			write(b, "#Name\tLength\tGC\tHH\tCAGA\tDepth\tStart\tTaxID\tTaxID2\tEnd\n"+
					"b\t5\t0.5\t0.5\t0.5\t2\t11\t12\t13\t16\n");
			final ScalarData merged=data(true).readTSV(input(a));
			merged.add(data(true).readTSV(input(b)));
			final Path out=dir.resolve("merged.tsv");
			merged.print(out.toString(), true, true, true, 1);
			final java.util.List<String> lines=Files.readAllLines(out, StandardCharsets.UTF_8);
			final String[] header=lines.get(0).split("\t", -1);
			require(header.length==10 && header[7].equals("TaxID") &&
					header[8].equals("TaxID2") && header[9].equals("End"),
					"End must be appended after existing taxonomy columns when printPos=t");
			final ScalarData reread=data(true).readTSV(input(out));
			final ScalarData fromRow=data(true);
			fromRow.addRow(row("manual", 9, 21, 33));
			writeActual(dir.resolve("merge_roundtrip.actual.tsv"), reread);
			writeActual(dir.resolve("add_row.actual.tsv"), fromRow);
			requireRows(reread, 2);
			requireRow(reread, 0, 3, 4, 7);
			requireRow(reread, 1, 5, 11, 16);
			requireRows(fromRow, 1);
			requireRow(fromRow, 0, 9, 21, 33);
		}else{throw new IllegalArgumentException(mode);}
	}

	/** Generates a scalar table from a finite FASTA input. */
	private static ScalarData intervals(Path path, String fasta,
			int window, int interval, int minlen, boolean breakOnContig) throws Exception{
		write(path, fasta);
		return ScalarIntervals.toIntervals(
				FileFormat.testInput(path.toString(), FileFormat.FASTA, null, true, true),
				window, interval, minlen, breakOnContig, -1);
	}

	/** Creates an instance without adding a production-visible constructor. */
	private static ScalarData data(boolean names) throws Exception{
		final Constructor<ScalarData> c=ScalarData.class.getDeclaredConstructor(
				boolean.class, long.class);
		c.setAccessible(true);
		return c.newInstance(names, -1L);
	}

	/** Returns a literal row, including future End when ScalarInterval exposes it. */
	private static ScalarInterval row(String name, float length, int start, int end)
			throws Exception{
		final ScalarInterval row=new ScalarInterval();
		row.name=name; row.length=length; row.start=start;
		row.gc=0.5f; row.hh=0.25f; row.caga=0.75f; row.depth=6;
		row.taxID=7; row.taxID2=8;
		setEnd(row, end);
		return row;
	}

	/** Resolves scalar TSV inputs. */
	private static FileFormat input(Path path){
		return FileFormat.testInput(path.toString(), FileFormat.TXT, null, false, true);
	}

	/** Resets scalar globals touched by these fresh-JVM cases. */
	private static void resetStatics() throws Exception{
		set(ScalarIntervals.class, "printName", true);
		set(ScalarIntervals.class, "mt", false);
		set(ScalarIntervals.class, "sendInThread", false);
		set(ScalarIntervals.class, "covFile", null);
		set(ScalarIntervals.class, "depthFile", null);
		set(ScalarIntervals.class, "readsProcessed", 0L);
		set(ScalarIntervals.class, "basesProcessed", 0L);
		set(ScalarIntervals.class, "errorState", false);
		set(ScalarData.class, "depthMap", null);
		ScalarData.makeSketch=false;
		ScalarData.makeClade=false;
	}

	/** Sets package-private static fields without widening production access. */
	private static void set(Class<?> clazz, String name, Object value) throws Exception{
		final Field f=clazz.getDeclaredField(name);
		f.setAccessible(true);
		f.set(null, value);
	}

	/** Returns End if the production list exists, or unknown for the pre-End baseline. */
	private static int endAt(ScalarData data, int index) throws Exception{
		try{
			final Field f=ScalarData.class.getDeclaredField("end");
			f.setAccessible(true);
			final IntList ends=(IntList)f.get(data);
			return ends!=null && ends.size>index ? ends.get(index) : -1;
		}catch(NoSuchFieldException nsfe){return -1;}
	}

	/** Returns row.end when present in the fixed source, or unknown on the baseline. */
	private static int intervalEnd(ScalarInterval row) throws Exception{
		try{
			final Field f=ScalarInterval.class.getDeclaredField("end");
			f.setAccessible(true);
			return f.getInt(row);
		}catch(NoSuchFieldException nsfe){return -1;}
	}

	/** Sets ScalarInterval.end when present in the fixed source. */
	private static void setEnd(ScalarInterval row, int end) throws Exception{
		try{
			final Field f=ScalarInterval.class.getDeclaredField("end");
			f.setAccessible(true);
			f.setInt(row, end);
		}catch(NoSuchFieldException nsfe){}
	}

	/** Requires generated coordinates to be unknown for every row. */
	private static void requireUnknown(ScalarData data) throws Exception{
		for(int i=0; i<data.gc.size(); i++){requireRow(data, i, data.length.get(i), -1, -1);}
	}

	/** Requires the row count of every row-aligned core list. */
	private static void requireRows(ScalarData data, int rows){
		require(data.gc.size()==rows && data.hh.size()==rows && data.caga.size()==rows &&
				data.depth.size()==rows && data.length.size()==rows && data.start.size()==rows,
				"row count mismatch");
	}

	/** Requires one row's represented length and half-open metric span. */
	private static void requireRow(ScalarData data, int row, float length, int start, int end)
			throws Exception{
		final ScalarInterval si=data.getRow(row, null);
		require(si.length==length, "row "+row+" length: expected "+length+" got "+si.length);
		require(si.start==start, "row "+row+" start: expected "+start+" got "+si.start);
		final int observedEnd=endAt(data, row);
		require(observedEnd==end, "row "+row+" end: expected "+end+" got "+observedEnd);
		final int rowEnd=intervalEnd(si);
		require(rowEnd==end, "row "+row+" interval end: expected "+end+" got "+rowEnd);
	}

	/** Writes observed rows for baseline diagnostics before assertions can fail. */
	private static void writeActual(Path path, ScalarData data) throws Exception{
		final StringBuilder sb=new StringBuilder("#Row\tName\tLength\tStart\tEnd\tGC\tHH\tCAGA\n");
		for(int i=0; i<data.gc.size(); i++){
			final ScalarInterval row=data.getRow(i, null);
			sb.append(i).append('\t').append(row.name).append('\t').append(row.length).
					append('\t').append(row.start).append('\t').append(endAt(data, i)).
					append('\t').append(row.gc).append('\t').append(row.hh).
					append('\t').append(row.caga).append('\n');
		}
		write(path, sb.toString());
		System.out.print(path.getFileName()+"\n"+sb);
	}

	/** Writes deterministic UTF-8 fixture input. */
	private static void write(Path path, String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.UTF_8));
	}

	/** These checks remain active under both -ea and -da. */
	private static void require(boolean valid, String message){
		if(!valid){throw new AssertionError(message);}
	}
}
