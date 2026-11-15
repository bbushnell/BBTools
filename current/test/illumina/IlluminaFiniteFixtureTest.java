package test.illumina;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Locale;
import java.util.zip.GZIPOutputStream;

import illumina.Cbcl2Text;
import illumina.CbclDecoder;
import illumina.ClusterData;
import illumina.FilterReader;

/** Finite Illumina CBCL/filter/RunInfo checks.
 * @author Brian Bushnell, Jean */
public class IlluminaFiniteFixtureTest{

	/** Accepts {@code phase=baseline|fixed}, a finite {@code scenario}, and {@code outdir=fresh}. */
	public static void main(String[] args) throws Exception{
		String phase=null, scenario=null, directory=null;
		for(String arg : args){
			final String[] pair=arg.split("=", 2);
			require(pair.length==2, "Use name=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("scenario")){scenario=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException(arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "phase=baseline|fixed required");
		require(scenario!=null, "scenario is required");
		final boolean fixed=phase.equals("fixed");
		Locale.setDefault(Locale.ROOT);
		resetClusterTerms();
		if(scenario.equals("cbcl")){
			require(directory!=null, "cbcl requires fresh outdir");
			cbcl(fresh(directory), fixed);
		}else if(scenario.equals("cbcl-remap")){
			require(directory!=null, "cbcl-remap requires fresh outdir");
			cbclRemap(fresh(directory), fixed);
		}else if(scenario.equals("cbcl-truncate")){
			require(directory!=null, "cbcl-truncate requires fresh outdir");
			cbclTruncate(fresh(directory), fixed);
		}else if(scenario.equals("cbcl-compacted")){
			require(directory!=null, "cbcl-compacted requires fresh outdir");
			cbclCompacted(fresh(directory), fixed);
		}else if(scenario.equals("runinfo-r2")){
			require(directory!=null, "runinfo-r2 requires fresh outdir");
			runInfoR2(fresh(directory), fixed);
		}else if(scenario.equals("runinfo-i2")){
			require(directory!=null, "runinfo-i2 requires fresh outdir");
			runInfoI2(fresh(directory));
		}else if(scenario.equals("filter-truncate")){
			require(directory!=null, "filter-truncate requires fresh outdir");
			filterTruncate(fresh(directory), fixed);
		}else{throw new IllegalArgumentException("Unknown scenario: "+scenario);}
	}

	/** Two gzip members in one CBCL file must decode the requested tile, not tile zero. */
	private static void cbcl(File out, boolean fixed) throws Exception{
		final File cbcl=new File(out, "two-tiles.cbcl");
		writeCbcl(cbcl, new int[]{1101, 2202}, new int[]{1, 1},
				new byte[][]{{call(0, 1)}, {call(3, 3)}});
		boolean threw=false;
		byte[][] actual=null;
		try{actual=CbclDecoder.readTile(cbcl.getPath(), 2202);}
		catch(IOException e){threw=true;}
		if(fixed){
			require(!threw, "Tile 2202 should decode after the fix");
			final String observed=ascii(actual[0])+"/"+ascii(actual[1]);
			require("T/$".equals(observed), "Tile 2202 decoded "+observed);
		}else{
			require(threw, "Vendor-shaped two-tile CBCL should fail before the fix");
		}
		System.out.println("cbcl-two-tile\tPASS\tthrew="+threw);
	}

	/** Nonzero CBCL quality remaps must be parsed as interleaved boundary/remap pairs. */
	private static void cbclRemap(File out, boolean fixed) throws Exception{
		final File cbcl=new File(out, "remap.cbcl");
		writeCbcl(cbcl, new int[]{1101}, new int[]{1}, new byte[][]{{call(2, 2)}},
				new int[]{0, 12, 30, 37}, false);
		final byte[][] actual=CbclDecoder.readTile(cbcl.getPath(), 1101);
		final String observed=ascii(actual[0])+"/"+ascii(actual[1]);
		require((fixed ? "G/?" : "G/$").equals(observed), "Remapped CBCL decoded "+observed);
		System.out.println("cbcl-remap\tPASS\t"+observed);
	}

	/** Header and block sizes must match available vendor tile records and gzip data. */
	private static void cbclTruncate(File out, boolean fixed) throws Exception{
		final File shortHeader=new File(out, "short-header.cbcl");
		writeShortHeaderCbcl(shortHeader);
		final File shortBlock=new File(out, "short-block.cbcl");
		writeCbcl(shortBlock, new int[]{1101}, new int[]{1}, new byte[][]{{call(1, 1)}},
				null, false, 100, -1);
		final File wrongRawSize=new File(out, "wrong-uncompressed.cbcl");
		writeCbcl(wrongRawSize, new int[]{1101}, new int[]{1}, new byte[][]{{call(1, 1)}},
				null, false, -1, 2);
		final boolean headerThrew=throwsCbcl(shortHeader);
		final boolean blockThrew=throwsCbcl(shortBlock);
		final boolean rawSizeThrew=throwsCbcl(wrongRawSize);
		require(headerThrew, "Truncated vendor tile record should throw");
		require(blockThrew==fixed, "Truncated CBCL block result changed: "+blockThrew);
		require(rawSizeThrew==fixed, "Declared uncompressed size result changed: "+rawSizeThrew);
		System.out.println("cbcl-truncate\tPASS\theader="+headerThrew+
				"\tblock="+blockThrew+"\trawsize="+rawSizeThrew);
	}

	/** Raw compacted tile decode is supported, but positioned text needs a cluster index map. */
	private static void cbclCompacted(File out, boolean fixed) throws Exception{
		final File run=singleCycleRun(out, true);
		final File cbcl=new File(run, "Data/Intensities/BaseCalls/L001/C1.1/L001_1.cbcl");
		boolean rawThrew=false;
		byte[][] actual=null;
		try{actual=CbclDecoder.readTile(cbcl.getPath(), 1101);}
		catch(IOException e){rawThrew=true;}
		require(!rawThrew, "Raw compacted CBCL decode should not assign coordinates or throw");
		final String observed=ascii(actual[0])+"/"+ascii(actual[1]);
		require("A/\"".equals(observed), "Raw compacted CBCL decoded "+observed);

		final boolean textThrew=throwsCompactedPosition(run);
		require(textThrew==fixed, "Compacted coordinate-assignment guard result changed: "+textThrew);
		System.out.println("cbcl-compacted\tPASS\trawThrew="+rawThrew+"\ttextThrew="+textThrew);
	}

	/** Invoke coordinate assignment directly; process() converts IOExceptions to KillSwitch kills. */
	private static boolean throwsCompactedPosition(File run) throws Exception{
		final Cbcl2Text converter=new Cbcl2Text(run.getPath(),
				new File(run, "compacted.fq").getPath(), 1, "1101", "auto");
		final Method processTile=Cbcl2Text.class.getDeclaredMethod(
				"processTile", int.class, int.class, float[][].class);
		processTile.setAccessible(true);
		try{
			processTile.invoke(converter, 1, 1101, new float[][]{{1f, 2f}});
			return false;
		}catch(InvocationTargetException e){
			Throwable cause=e.getCause();
			if(cause instanceof IOException){return true;}
			if(cause instanceof Exception){throw (Exception)cause;}
			if(cause instanceof Error){throw (Error)cause;}
			throw new AssertionError(cause);
		}
	}

	/** Three read terms R1,I1,R2 must not demote R2 into a second index barcode. */
	private static void runInfoR2(File out, boolean fixed) throws Exception{
		final File run=threeTermRun(out, new boolean[]{false, true, false});
		final File fastq=new File(out, "runinfo-r2.fq");
		new Cbcl2Text(run.getPath(), fastq.getPath(), 1, "1101", "auto").process();
		final String observed=readAscii(fastq);
		final String expected=fixed ?
				"@machine:run:flowcell:1:1101:1010:1020: 1:N:0:C\n"+
				"A\n+\n\"\n"+
				"@machine:run:flowcell:1:1101:1010:1020: 2:N:0:C\n"+
				"G\n+\n$\n"
				:
				"@machine:run:flowcell:1:1101:1010:1020: 1:N:0:C+G\n"+
				"A\n+\n\"\n";
		require(expected.equals(observed), "Unexpected RunInfo R2 output:\n"+observed);
		System.out.println("runinfo-r2\tPASS\trecords="+records(observed));
	}

	/** Three read terms R1,I1,I2 must not invent a second FASTQ record. */
	private static void runInfoI2(File out) throws Exception{
		final File run=threeTermRun(out, new boolean[]{false, true, true});
		final File fastq=new File(out, "runinfo-i2.fq");
		new Cbcl2Text(run.getPath(), fastq.getPath(), 1, "1101", "auto").process();
		final String observed=readAscii(fastq);
		final String expected=
				"@machine:run:flowcell:1:1101:1010:1020: 1:N:0:C+G\n"+
				"A\n+\n\"\n";
		require(expected.equals(observed), "Unexpected RunInfo I2 output:\n"+observed);
		System.out.println("runinfo-i2\tPASS\trecords="+records(observed));
	}

	/** The declared filter count must match the number of following flag bytes. */
	private static void filterTruncate(File out, boolean fixed) throws Exception{
		final File good=new File(out, "complete.filter");
		writeFilter(good, 3, new byte[]{1, 0, 1});
		require(Arrays.equals(FilterReader.readFilters(good.getPath()),
				new boolean[]{true, false, true}), "Complete filter control");

		final File shortFile=new File(out, "truncated.filter");
		writeFilter(shortFile, 3, new byte[]{1});
		boolean threw=false;
		boolean[] actual=null;
		try{actual=FilterReader.readFilters(shortFile.getPath());}
		catch(Exception e){threw=true;}
		if(fixed){
			require(threw, "Truncated filter should fail loudly");
		}else{
			require(!threw && Arrays.equals(actual, new boolean[]{true, false, false}),
					"Historical truncation result changed: "+Arrays.toString(actual));
		}
		System.out.println("filter-truncate\tPASS\tthrew="+threw);
	}

	private static File threeTermRun(File out, boolean[] indexed) throws Exception{
		final File run=new File(out, "run");
		final File intensities=new File(run, "Data/Intensities");
		final File basecalls=new File(intensities, "BaseCalls/L001");
		require(new File(basecalls, "C1.1").mkdirs(), "Could not create run directories");
		writeLocs(new File(intensities, "s.locs"), new float[][]{{1f, 2f}});
		writeFilter(new File(basecalls, "s_1_1101.filter"), 1, new byte[]{1});
		writeRunInfo(run, new int[]{1, 1, 1}, indexed);
		writeCbcl(new File(basecalls, "C1.1/L001_1.cbcl"), new int[]{1101}, new int[]{1},
				new byte[][]{{call(0, 1)}});
		writeCbcl(new File(basecalls, "C2.1/L001_1.cbcl"), new int[]{1101}, new int[]{1},
				new byte[][]{{call(1, 2)}});
		writeCbcl(new File(basecalls, "C3.1/L001_1.cbcl"), new int[]{1101}, new int[]{1},
				new byte[][]{{call(2, 3)}});
		return run;
	}

	private static File singleCycleRun(File out, boolean nonPfExcluded) throws Exception{
		final File run=new File(out, "run");
		final File intensities=new File(run, "Data/Intensities");
		final File basecalls=new File(intensities, "BaseCalls/L001");
		require(new File(basecalls, "C1.1").mkdirs(), "Could not create run directories");
		writeLocs(new File(intensities, "s.locs"), new float[][]{{1f, 2f}});
		writeFilter(new File(basecalls, "s_1_1101.filter"), 1, new byte[]{1});
		writeRunInfo(run, new int[]{1}, new boolean[]{false});
		writeCbcl(new File(basecalls, "C1.1/L001_1.cbcl"), new int[]{1101}, new int[]{1},
				new byte[][]{{call(0, 1)}}, null, nonPfExcluded);
		return run;
	}

	private static void writeRunInfo(File run, int[] cycles, boolean[] indexed) throws Exception{
		final StringBuilder xml=new StringBuilder();
		xml.append("<RunInfo><Run><Reads>\n");
		for(int i=0; i<cycles.length; i++){
			xml.append("<Read Number=\"").append(i+1).append("\" NumCycles=\"").append(cycles[i])
					.append("\" IsIndexedRead=\"").append(indexed[i] ? "Y" : "N").append("\"/>\n");
		}
		xml.append("</Reads></Run></RunInfo>\n");
		Files.write(new File(run, "RunInfo.xml").toPath(),
				xml.toString().getBytes(StandardCharsets.US_ASCII));
	}

	private static void writeLocs(File file, float[][] positions) throws Exception{
		final ByteArrayOutputStream out=new ByteArrayOutputStream();
		writeInt(out, 1);
		writeInt(out, Float.floatToIntBits(1f));
		writeInt(out, positions.length);
		for(float[] position : positions){
			writeInt(out, Float.floatToIntBits(position[0]));
			writeInt(out, Float.floatToIntBits(position[1]));
		}
		Files.write(file.toPath(), out.toByteArray());
	}

	private static void writeFilter(File file, int clusters, byte[] flags) throws Exception{
		final ByteArrayOutputStream out=new ByteArrayOutputStream();
		writeInt(out, 0);
		writeInt(out, 3);
		writeInt(out, clusters);
		out.write(flags);
		Files.write(file.toPath(), out.toByteArray());
	}

	private static void writeCbcl(File file, int[] tileNums, int[] clusterCounts,
			byte[][] rawBlocks) throws Exception{
		writeCbcl(file, tileNums, clusterCounts, rawBlocks, null, false);
	}

	private static void writeCbcl(File file, int[] tileNums, int[] clusterCounts,
			byte[][] rawBlocks, int[] qscoreRemap, boolean nonPfExcluded) throws Exception{
		writeCbcl(file, tileNums, clusterCounts, rawBlocks, qscoreRemap, nonPfExcluded,
				-1, -1);
	}

	private static void writeCbcl(File file, int[] tileNums, int[] clusterCounts,
			byte[][] rawBlocks, int[] qscoreRemap, boolean nonPfExcluded,
			int overrideCompressedSize, int overrideRawSize) throws Exception{
		require(tileNums.length==clusterCounts.length && tileNums.length==rawBlocks.length,
				"Mismatched CBCL metadata");
		final int bins=(qscoreRemap==null ? 0 : qscoreRemap.length);
		require(file.getParentFile().isDirectory() || file.getParentFile().mkdirs(),
				"Could not create "+file.getParent());
		final byte[][] compressedBlocks=new byte[rawBlocks.length][];
		for(int i=0; i<rawBlocks.length; i++){compressedBlocks[i]=gzip(rawBlocks[i]);}
		final ByteArrayOutputStream out=new ByteArrayOutputStream();
		writeShort(out, 1);
		writeInt(out, 8+4+bins*8+4+tileNums.length*16+1);
		out.write(2);
		out.write(2);
		writeInt(out, bins);
		for(int i=0; i<bins; i++){
			writeInt(out, i);
			writeInt(out, qscoreRemap[i]);
		}
		writeInt(out, tileNums.length);
		for(int i=0; i<tileNums.length; i++){
			writeInt(out, tileNums[i]);
			writeInt(out, clusterCounts[i]);
			writeInt(out, overrideRawSize>=0 ? overrideRawSize : rawBlocks[i].length);
			writeInt(out, overrideCompressedSize>=0 ? overrideCompressedSize : compressedBlocks[i].length);
		}
		out.write(nonPfExcluded ? 1 : 0);
		for(byte[] compressedBlock : compressedBlocks){out.write(compressedBlock);}
		Files.write(file.toPath(), out.toByteArray());
	}

	private static void writeShortHeaderCbcl(File file) throws Exception{
		final ByteArrayOutputStream out=new ByteArrayOutputStream();
		writeShort(out, 1);
		writeInt(out, 8+4+8+1);
		out.write(2);
		out.write(2);
		writeInt(out, 0);
		writeInt(out, 1);
		writeInt(out, 1101);
		Files.write(file.toPath(), out.toByteArray());
	}

	private static boolean throwsCbcl(File file) {
		try{
			CbclDecoder.readTile(file.getPath(), 1101);
			return false;
		}catch(IOException e){
			return true;
		}
	}

	private static byte[] gzip(byte[] raw) throws Exception{
		final ByteArrayOutputStream baos=new ByteArrayOutputStream();
		final GZIPOutputStream gzip=new GZIPOutputStream(baos);
		gzip.write(raw);
		gzip.close();
		return baos.toByteArray();
	}

	private static void writeShort(ByteArrayOutputStream out, int value){
		out.write(value&0xFF);
		out.write((value>>>8)&0xFF);
	}

	private static void writeInt(ByteArrayOutputStream out, int value){
		out.write(value&0xFF);
		out.write((value>>>8)&0xFF);
		out.write((value>>>16)&0xFF);
		out.write((value>>>24)&0xFF);
	}

	private static byte call(int base, int quality){
		return (byte)(base|(quality<<2));
	}

	private static String ascii(byte[] array){
		return new String(array, StandardCharsets.US_ASCII);
	}

	private static String readAscii(File file) throws Exception{
		return new String(Files.readAllBytes(file.toPath()), StandardCharsets.US_ASCII);
	}

	private static int records(String fastq){
		int count=0;
		for(int i=0; i<fastq.length(); i++){
			if(fastq.charAt(i)=='@'){count++;}
		}
		return count;
	}

	private static File fresh(String directory){
		final File out=new File(directory);
		require(out.isAbsolute() && out.mkdir(), "Cannot create fresh outdir: "+out);
		return out;
	}

	private static void resetClusterTerms(){
		ClusterData.r1term=0;
		ClusterData.i1term=1;
		ClusterData.i2term=2;
		ClusterData.r2term=3;
		ClusterData.u1term=-1;
		ClusterData.u2term=-1;
		ClusterData.machineID="machine";
		ClusterData.runID="run";
		ClusterData.flowcellID="flowcell";
		ClusterData.barcodeDelimiter='+';
		ClusterData.appendBarcodeQuality=false;
	}

	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
