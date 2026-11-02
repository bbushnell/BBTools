package test.scalar;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import scalar.CloudPlot;
import scalar.ScalarData;
import scalar.ScalarInterval;
import scalar.ScalarIntervals;
import stream.ConcurrentReadInputStream;

/** Finite, network-free scalar table and worker contract regressions. */
public class ScalarDataContractTest{

	/** Runs one independently specified case in a fresh JVM. */
	public static void main(String[] args) throws Exception{
		final String mode=args[0];
		final Path dir=java.nio.file.Paths.get(args[1]);
		Files.createDirectories(dir);
		if(mode.startsWith("roundtrip")){
			final boolean names=mode.contains("name"), pos=mode.contains("pos");
			final ScalarData original=data(true);
			original.addRow(row("contig_A", 120, 37));
			final Path file=dir.resolve("table.tsv");
			final ByteStreamWriter writer=ByteStreamWriter.makeBSW(
					FileFormat.testOutput(file.toString(), FileFormat.TXT, null, false, true, false, false));
			original.print(writer, true, names, pos, 10);
			require(!writer.poisonAndWait(), "test table writer failed");
			final ScalarData loaded=data(true).readTSV(input(file));
			final ScalarInterval result=loaded.getRow(0, null);
			require(result.gc==0.5f && result.hh==0.25f && result.caga==0.75f, "metric columns shifted");
			require(result.length==120 && result.depth==6 && result.taxID==7 && result.taxID2==8, "optional columns shifted");
			require(result.start==(pos ? 37 : 0), "start column mismatch");
			require(names ? "contig_A".equals(result.name) : result.name==null, "name column mismatch");
		}else if(mode.equals("null-name")){
			final ScalarData d=data(true);
			d.addRow(row(null, 120, 0));
			d.addRow(row("second", 140, 10));
			require(d.names.size()==2 && d.name(0)==null && "second".equals(d.name(1)), "null name shifted parallel columns");
		}else if(mode.equals("names-disabled-print")){
			final ScalarData d=data(false);
			d.addRow(row("discarded", 120, 37));
			require(d.names==null, "Fixture must exercise disabled name storage");
			final Path file=dir.resolve("names-disabled.tsv");
			d.print(file.toString(), true, true, true, 10);
			final java.util.List<String> lines=Files.readAllLines(file, StandardCharsets.UTF_8);
			require(lines.size()==2, "Expected one header and one data row");
			final String[] header=lines.get(0).split("\t", -1), row=lines.get(1).split("\t", -1);
			require(header.length==10 && row.length==10 && row[0].isEmpty(),
					"Requested Name column must retain an empty field when names were not stored");
			require(Float.parseFloat(row[1])==120 && Float.parseFloat(row[2])==0.5f
					&& Integer.parseInt(row[6])==37 && Integer.parseInt(row[8])==8,
					"Missing name shifted numeric fields");
		}else if(mode.equals("legacy-four") || mode.equals("legacy-named")){
			final Path file=dir.resolve("legacy.tsv");
			final boolean named=mode.equals("legacy-named");
			write(file, named ? "#GC\tHH\tCAGA\tTaxID\tName\n0.5\t0.25\t0.75\t7\tcontig\n" : "0.5\t0.25\t0.75\t7\n");
			final ScalarData d=data(true).readTSV(input(file));
			final ScalarInterval result=d.getRow(0, null);
			require(result.taxID==7 && result.depth==0 && result.length==0, "legacy columns shifted");
			require(named ? "contig".equals(result.name) : result.name==null, "legacy missing-name handling");
		}else if(mode.equals("empty-coverage")){
			final Path file=dir.resolve("empty.cov");
			write(file, "");
			try{data(false).loadDepth(file.toString(), null);}
			catch(IllegalArgumentException expected){System.out.println("PASS "+mode); return;}
			throw new AssertionError("Empty coverage must be rejected with a clear argument error");
		}else if(mode.equals("writer-error")){
			final ScalarData d=data(true);
			d.addRow(row("contig", 10, 0));
			try{d.print("/dev/full", true, true, true, 1);}
			catch(RuntimeException expected){System.out.println("PASS "+mode); return;}
			throw new AssertionError("Owned writer failure was not propagated");
		}else if(mode.equals("null-output")){
			data(false).print((FileFormat)null, true, false, false, 0);
		}else if(mode.equals("empty-collapse")){
			final ScalarData d=ScalarIntervals.collapse(new ArrayList<ScalarData>(), true);
			require(d!=null && d.gc.size()==0, "empty input must yield an empty table");
		}else if(mode.equals("order-control")){
			final int[] order=CloudPlot.parseOrder("gc,hh,caga");
			require(order.length==4 && order[0]==0 && order[1]==1 && order[2]==2 && order[3]==-1, "legacy order changed");
		}else if(mode.startsWith("order-")){
			final String value=mode.substring(6).replace('_', ',');
			try{CloudPlot.parseOrder(value);}
			catch(IllegalArgumentException expected){System.out.println("PASS "+mode); return;}
			throw new AssertionError("Invalid dimension count/empty metric accepted: "+value);
		}else if(mode.equals("paired-serial") || mode.equals("worker-no-send")){
			set("printName", true);
			set("mt", mode.equals("worker-no-send"));
			set("sendInThread", false);
			set("errorState", false);
			ScalarData.makeSketch=false;
			ScalarData.makeClade=false;
			final Path a=dir.resolve("a.fa"), b=dir.resolve("b.fa");
			write(a, ">first\nACGTACGT\n"); write(b, ">second\nAAAACCCC\n");
			final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, input(a), input(b));
			cris.start();
			final ArrayList<ScalarData> tables;
			try{tables=ScalarIntervals.toIntervals_multi(cris, 0, 0, 1, -1, true, -1);}
			finally{require(!fileIO.ReadWrite.closeStreams(cris), "fixture input failed");}
			final ScalarData d=ScalarIntervals.collapse(tables, true);
			require(d!=null && d.gc.size()==2, "both mates must be represented");
			require(!((Boolean)get("errorState")), "successful workers reported failure without in-thread labeling");
		}else{throw new IllegalArgumentException(mode);}
		System.out.println("PASS "+mode);
	}

	/** Creates an instance without changing production access just for the baseline. */
	private static ScalarData data(boolean names) throws Exception{
		final Constructor<ScalarData> c=ScalarData.class.getDeclaredConstructor(boolean.class, long.class);
		c.setAccessible(true);
		return c.newInstance(names, -1L);
	}
	/** Returns a literal row with distinct optional columns. */
	private static ScalarInterval row(String name, float length, int start){
		final ScalarInterval r=new ScalarInterval();
		r.name=name; r.length=length; r.start=start;
		r.gc=0.5f; r.hh=0.25f; r.caga=0.75f; r.depth=6;
		r.taxID=7; r.taxID2=8;
		return r;
	}
	/** Resolves the fixture format, including FASTA for paired tests. */
	private static FileFormat input(Path path){
		return FileFormat.testInput(path.toString(), FileFormat.TXT, null, false, true);
	}
	/** Writes deterministic UTF-8 fixture input. */
	private static void write(Path path, String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.UTF_8));
	}
	/** Sets an existing package-private configuration field for isolated worker tests. */
	private static void set(String name, Object value) throws Exception{
		final Field f=ScalarIntervals.class.getDeclaredField(name);
		f.setAccessible(true); f.set(null, value);
	}
	/** Reads the worker error result without changing its production visibility. */
	private static Object get(String name) throws Exception{
		final Field f=ScalarIntervals.class.getDeclaredField(name);
		f.setAccessible(true); return f.get(null);
	}
	/** These checks remain active under both -ea and -da. */
	private static void require(boolean valid, String message){if(!valid){throw new AssertionError(message);}}
}
