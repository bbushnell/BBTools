package test.hiseq;

import java.io.File;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ReadWrite;
import hiseq.FlowCell;
import hiseq.MicroTile;
import hiseq.Tile;
import hiseq.TileDump;

/** Literal old-format inputs and exact scalar-counter checks for the optional v4 extension.
 * Reflection keeps the same fixture compilable against the unmodified v3 writer.
 * This does not claim persistence of cycle arrays or all FlowCell metadata.
 * @author Brian Bushnell, Sucrose */
public class HiseqDumpCompatibilityTest{

	/** Requires phase=baseline|fixed and a fresh absolute outdir with an existing parent. */
	public static void main(String[] args) throws Exception{
		String phase=null, directory=null;
		for(String arg : args){
			final String[] pair=arg.split("=", 2);
			require(pair.length==2, "Use name=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException(arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "phase=baseline|fixed required");
		final File out=new File(directory);
		require(out.isAbsolute() && out.mkdir(), "Use a fresh absolute outdir: "+out);
		Locale.setDefault(Locale.ROOT);
		Tile.xSize=Tile.ySize=500;
		MicroTile.TRACK_CYCLES=false;
		ByteFile.FORCE_MODE_BF1=true; ByteFile.FORCE_MODE_BF2=false;
		ReadWrite.USE_PIGZ=false; ReadWrite.USE_UNPIGZ=false;
		ReadWrite.USE_BGZIP=false; ReadWrite.USE_UNBGZIP=false;
		legacy(out);
		final FlowCell cell=cell();
		final File defaultFile=new File(out, "default-v3.tsv");
		cell.dump(defaultFile.getPath(), false);
		require(Files.readAllLines(defaultFile.toPath(), StandardCharsets.US_ASCII).get(0).equals("#Version\t3"),
			"Existing writer signature must default to v3");
		final MicroTile old=TileDump.loadDump(defaultFile.getPath()).toList().get(0);
		require(old.misses==0 && cell.toList().get(0).misses==1, "Retain v3's documented approximate reconstruction");
		System.out.println("v3-default\tPASS\tmisses=1->0");
		if(phase.equals("fixed")){exact(out, cell, defaultFile);}
	}

	/** Independent positional literals exercise each old reader, including optional merged fields. */
	private static void legacy(File out) throws Exception{
		final String v1="0\t0\t0\t499\t0\t499\t3\t2\t8\t1\t2\t1\t2\t1\t0"
			+"\t50\t20.5\t90\t2\t0\t23\t0\t0\t0\t0\t2\t3\t3\t3\t2\t1\t1\t5\t0\n";
		final MicroTile one=load(out, "literal-v1.tsv", 1, v1).toList().get(0);
		require(one.readCount==3 && one.alignedReadCount==2 && one.baseErrorCount==2 && one.discard==2,
			"Literal v1 integer positions");
		require(one.misses==1 && one.hits==2 && one.depthSum==6 && one.readQualityByProbSum==61.5,
			"Literal v1 reconstruction");
		final String common="0\t0\t0\t499\t0\t499\t3\t12\t2\t8\t1\t2\t1\t2\t1\t0"
			+"\t0\t0\t50\t20.5\t90\t0.125\t2\t0\t0\t0\t0\t0\t23\t0\t0\t0\t0\t2"
			+"\t3\t3\t3\t2\t1\t1\t5\t0\t3\t2\t1\t0\t0";
		for(int version=2; version<=3; version++){
			final String row=common+(version==3 ? "\t5\t1.25\t0.25" : "");
			for(int merged=0; merged<=1; merged++){
				final MicroTile mt=load(out, "literal-v"+version+"-"+merged+".tsv", version,
					row+(merged==1 ? "\t2\t50\t6\t1\t50\t0.5\t0.1" : "")+"\n").toList().get(0);
				require(mt.readCount==3 && mt.baseCount==12 && mt.alignedBaseCount==8 && mt.discard==2,
					"Literal v"+version+" integer positions");
				require(mt.misses==1 && mt.hits==2 && mt.depthSum==6 && mt.baseErrorProbSum==1.5,
					"Literal old reconstruction");
				require(mt.validKmerSum==(version==3 ? 5 : 0) && mt.goodKmerSum==(version==3 ? 1.25 : 0),
					"Version-dependent kmer fields");
				require(value(mt, "barcodes").longValue()==3 && value(mt, "mergedReads").longValue()==2*merged
					&& value(mt, "insertSum").longValue()==50*merged, "Barcode and optional merge positions");
			}
		}
		System.out.println("literal-readers\tPASS\tv1,v2,v3,optional-merge");
	}

	/** Tests exact counters, old-reader prefix interpretation, both overloads, and fail-before-open checks. */
	private static void exact(File out, FlowCell cell, File defaultFile) throws Exception{
		final File explicit=new File(out, "explicit-v3.tsv"), four=new File(out, "exact-v4.tsv");
		write(cell, explicit, 3, false);
		require(Arrays.equals(Files.readAllBytes(defaultFile.toPath()), Files.readAllBytes(explicit.toPath())),
			"Explicit legacy writing must be byte-identical to the default");
		write(cell, four, 4, true);
		final FlowCell loaded=TileDump.loadDump(four.getPath());
		// Independent v4 input: a literal legacy row plus a literal raw suffix, never writer-produced.
		final String raw="2\t1\t9007199254740993\t61.123456789012345\t270.1234567890123\t1.0E-200\t1234.1234567890124";
		final MicroTile literal=load(out, "literal-v4.tsv", 4,
			data(new File(out, "literal-v3-1.tsv")).get(0)+"\t"+raw+"\n").toList().get(0);
		require(literal.hits==2 && literal.misses==1 && literal.depthSum==9007199254740993L,
			"Literal v4 exact integer positions");
		require(literal.readQualityByProbSum==61.123456789012345 && literal.probErrorFreeSum==270.1234567890123
			&& literal.baseErrorProbSum==1e-200 && literal.goodKmerSum==1234.1234567890124,
			"Literal v4 exact double positions");
		final FlowCell existing=new FlowCell(31);
		require(TileDump.loadDump(four.getPath(), existing)==existing, "Existing target overload returns its target");
		for(int i=0; i<cell.toList().size(); i++){
			equalCounters(cell.toList().get(i), loaded.toList().get(i));
			equalCounters(cell.toList().get(i), existing.toList().get(i));
		}
		final List<String> threeRows=data(defaultFile), fourRows=data(four);
		for(int i=0; i<threeRows.size(); i++){
			final String[] old=threeRows.get(i).split("\t"), current=fourRows.get(i).split("\t");
			require(current.length==64 && (old.length==50 || old.length==57), "57 legacy + 7 raw fields");
			for(int j=0; j<old.length; j++){require(old[j].equals(current[j]), "Legacy column changed at "+j);}
			final String[] expectedRaw=raw.split("\t");
			for(int j=0; j<7; j++){
				require(j<3 ? Long.parseLong(expectedRaw[j])==Long.parseLong(current[j+57])
					: Double.doubleToLongBits(Double.parseDouble(expectedRaw[j]))
					==Double.doubleToLongBits(Double.parseDouble(current[j+57])), "Raw suffix position "+j);
			}
		}
		// Emulate the pre-extension v3 reader exactly at its public row boundary: it ignores the raw suffix.
		final File rows=new File(out, "v4-rows-for-old-reader.tsv");
		Files.write(rows.toPath(), (String.join("\n", fourRows)+"\n").getBytes(StandardCharsets.US_ASCII));
		final ByteFile bf=ByteFile.makeByteFile(rows.getPath(), false);
		final FlowCell oldTarget=new FlowCell(31);
		try{require(TileDump.loadTiles2(bf, oldTarget, 3)==3, "Old-reader path consumes all v4 rows");}
		finally{require(!bf.close(), "Old-reader input closes successfully");}
		require(oldTarget.toList().get(0).misses==0, "Old-reader path retains approximate legacy semantics");
		require(value(oldTarget.toList().get(1), "mergedReads").longValue()==2, "Merged prefix survives old reading");
		final File again=new File(out, "v4-again.tsv");
		// Fits are deliberately outside the scalar-counter persistence contract; supply them for reporting.
		loaded.uniqueToReadErrorRateFormula=cell.uniqueToReadErrorRateFormula.clone();
		loaded.uniqueToBaseErrorRateFormula=cell.uniqueToBaseErrorRateFormula.clone();
		write(loaded, again, 4, false);
		final FlowCell twice=TileDump.loadDump(again.getPath());
		for(int i=0; i<cell.toList().size(); i++){equalCounters(cell.toList().get(i), twice.toList().get(i));}
		for(int version : new int[] {1, 2, 5}){
			final File bad=new File(out, "unsupported-"+version+".tsv");
			rejectedWrite(cell, bad, version);
		}
		for(String field : DOUBLES){
			final MicroTile mt=cell.toList().get(0);
			final Field f=MicroTile.class.getField(field);
			final double saved=f.getDouble(mt);
			for(double invalid : new double[] {Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY}){
				f.setDouble(mt, invalid);
				rejectedWrite(cell, new File(out, "nonfinite-"+field+"-"+Double.toString(invalid)+".tsv"), 4);
			}
			f.setDouble(mt, saved);
		}
		final String valid=fourRows.get(0);
		int invalidIndex=0;
		for(String suffix : new String[] {"\tNaN", "\tInfinity", "\t"}){
			final File bad=new File(out, "bad-input-"+(invalidIndex++)+".tsv");
			final String truncated=valid.substring(0, valid.lastIndexOf('\t'))+suffix;
			boolean threw=false;
			try{load(out, bad.getName(), 4, truncated+"\n");}catch(IllegalArgumentException e){threw=true;}
			require(threw, "Reject malformed/nonfinite exact counter");
		}
		System.out.println("exact-v4\tPASS\t3rows,all-scalar-counters,bit-exact-doubles,legacy-prefix,bad-inputs");
	}

	/** Three distinct counter states include unmerged, merged, and zero-read alignment-only tiles. */
	private static FlowCell cell() throws Exception{
		final FlowCell fc=new FlowCell(31);
		fc.uniqueToReadErrorRateFormula=new double[] {0.01, 0.0001};
		fc.uniqueToBaseErrorRateFormula=new double[] {0.01, 0.0001};
		for(int i=0; i<3; i++){
			final MicroTile mt=fc.getMicroTile(0, 0, i*500, 0);
			mt.readCount=(i==2 ? 0 : 3); mt.baseCount=(i==2 ? 0 : 12);
			mt.hits=2; mt.misses=1; mt.depthSum=9007199254740993L;
			mt.alignedReadCount=2; mt.alignedBaseCount=8; mt.readErrorCount=1; mt.baseErrorCount=2;
			mt.kmerReadErrorCount=1; mt.kmerBaseErrorCount=2; mt.readInsCount=1; mt.readDelCount=1;
			mt.readQualityByProbSum=61.123456789012345; mt.probErrorFreeSum=270.1234567890123;
			mt.baseErrorProbSum=1e-200; mt.validKmerSum=2000; mt.goodKmerSum=1234.1234567890124;
			mt.acgtn=new long[] {3, 3, 3, 2, 1}; mt.homoPolyGCount=1; mt.homoPolyGSum=5; mt.discard=2;
			set(mt, "barcodes", 3); set(mt, "barcodeHDistSum", 2); set(mt, "barcodePolymers", 1);
			set(mt, "mergedReads", i==1 ? 2 : 0); set(mt, "insertSum", i==1 ? 50 : 0);
			set(mt, "overlapSum", i==1 ? 6 : 0); set(mt, "mergeErrorSum", i==1 ? 1 : 0);
		}
		return fc;
	}

	private static FlowCell load(File out, String name, int version, String row) throws Exception{
		final File file=new File(out, name);
		Files.write(file.toPath(), ("#Version\t"+version+"\n#xSize\t500\n#ySize\t500\n"+row).getBytes(StandardCharsets.US_ASCII));
		return TileDump.loadDump(file.getPath());
	}

	/** Calls new APIs reflectively so baseline compilation and old-reader tests need no production edits. */
	private static void write(FlowCell fc, File file, int version, boolean flowCell) throws Exception{
		try{
			if(flowCell){FlowCell.class.getMethod("dump", String.class, boolean.class, int.class)
				.invoke(fc, file.getPath(), false, version);}
			else{TileDump.class.getMethod("write", FlowCell.class, String.class, boolean.class, int.class)
				.invoke(null, fc, file.getPath(), false, version);}
		}catch(InvocationTargetException e){
			if(e.getCause() instanceof RuntimeException){throw (RuntimeException)e.getCause();}
			throw e;
		}
	}

	private static void rejectedWrite(FlowCell cell, File out, int version) throws Exception{
		boolean threw=false;
		try{write(cell, out, version, false);}catch(IllegalArgumentException e){threw=true;}
		require(threw && !out.exists(), "Reject before creating output or starting its writer: "+out);
	}

	private static List<String> data(File file) throws Exception{
		final List<String> rows=new ArrayList<String>();
		for(String line : Files.readAllLines(file.toPath(), StandardCharsets.US_ASCII)){
			if(!line.startsWith("#")){rows.add(line);}
		}
		return rows;
	}

	/** Every scalar raw counter has an explicit name; derived state and cycle arrays are excluded. */
	private static void equalCounters(MicroTile expected, MicroTile actual) throws Exception{
		for(String field : LONGS){require(value(expected, field).longValue()==value(actual, field).longValue(), field);}
		for(String field : DOUBLES){require(Double.doubleToLongBits(value(expected, field).doubleValue())
			==Double.doubleToLongBits(value(actual, field).doubleValue()), field+" must retain exact double bits");}
		require(expected.discard==actual.discard && Arrays.equals(expected.acgtn, actual.acgtn), "Discard/base counts");
		require(expected.lane==actual.lane && expected.tile==actual.tile && expected.x1==actual.x1
			&& expected.x2==actual.x2 && expected.y1==actual.y1 && expected.y2==actual.y2, "Coordinates");
	}

	private static Number value(MicroTile tile, String name) throws Exception{
		final Field field=MicroTile.class.getDeclaredField(name);
		field.setAccessible(true);
		return (Number)field.get(tile);
	}

	private static void set(MicroTile tile, String name, long value) throws Exception{
		final Field field=MicroTile.class.getDeclaredField(name);
		field.setAccessible(true);
		field.setLong(tile, value);
	}

	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final String[] LONGS={"hits", "misses", "depthSum", "readCount", "baseCount", "alignedReadCount",
		"alignedBaseCount", "readErrorCount", "baseErrorCount", "kmerReadErrorCount", "kmerBaseErrorCount",
		"readInsCount", "readDelCount", "validKmerSum", "homoPolyGCount", "homoPolyGSum", "barcodes",
		"barcodeHDistSum", "barcodePolymers", "mergedReads", "insertSum", "overlapSum", "mergeErrorSum"};
	private static final String[] DOUBLES={"readQualityByProbSum", "probErrorFreeSum", "baseErrorProbSum", "goodKmerSum"};
}
