package test.hiseq;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Locale;
import java.util.zip.GZIPOutputStream;

import fileIO.ByteFile;
import fileIO.ReadWrite;
import hiseq.CycleTracker;
import hiseq.FlowCell;
import hiseq.IlluminaHeaderParser1;
import hiseq.IlluminaHeaderParser2;
import hiseq.MicroTile;
import hiseq.ReadHeaderParser;
import hiseq.Tile;
import hiseq.TileDump;
import stream.Read;

/** Finite model and dump checks; root runs the fatal legacy-parser case in its own JVM.
 * @author Brian Bushnell, Sucrose */
public class HiseqModelBoundaryTest{

	/** Accepts phase=baseline|fixed, scenario=models|parser-no-comment|dump|writer-full, outdir=fresh. */
	public static void main(String[] args) throws Exception{
		String phase=null, scenario="models", directory=null;
		for(String arg : args){
			final String[] pair=arg.split("=", 2);
			require(pair.length==2, "Use name=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("scenario")){scenario=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException(arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "phase=baseline|fixed required");
		final boolean fixed=phase.equals("fixed");
		Locale.setDefault(Locale.ROOT);
		ReadHeaderParser.PARSE_COORDINATES=true;
		ReadHeaderParser.PARSE_COMMENT=true;
		MicroTile.TRACK_CYCLES=false;
		if(scenario.equals("models")){models(fixed);}
		else if(scenario.equals("parser-no-comment")){
			ReadHeaderParser.PARSE_COMMENT=false;
			coordinates(new IlluminaHeaderParser1().parse(HEADER));
			System.out.println("parser-no-comment\tPASS\t1,1101,12,34");
		}else if(scenario.equals("dump")){
			require(directory!=null, "dump requires fresh outdir");
			final File out=new File(directory);
			require(out.isAbsolute() && out.mkdir(), "Cannot create fresh outdir: "+out);
			dump(out, fixed);
		}else if(scenario.equals("writer-full")){
			require(new File("/dev/full").exists(), "This failure case requires /dev/full on Dori");
			boolean threw=false;
			try{TileDump.write(cell(), "/dev/full", true);}
			catch(RuntimeException e){threw=true; System.out.println("writer-cause\t"+e.getClass().getName());}
			require(threw==fixed, "writer failure reporting: threw="+threw+", fixed="+fixed);
			System.out.println("writer-full\tPASS\tthrew="+threw);
		}else{throw new IllegalArgumentException("Unknown scenario: "+scenario);}
	}

	/** Literal coordinates and two reads give independent parser/composition controls. */
	private static void models(boolean fixed){
		coordinates(new IlluminaHeaderParser1().parse(HEADER+" 1:N:0:ACGT"));
		final IlluminaHeaderParser2 parser=new IlluminaHeaderParser2().parse(HEADER+" 1:N:0:ACGT");
		coordinates(parser);
		require(parser.canShrink() && !parser.chastityFail() && parser.pairnum()==0, "Ordinary header metadata");
		parser.parse(":::1:1101:12:34 1:");
		coordinates(parser);
		require(parser.looksShrunk() && !parser.canShrink(), "Actual shrinker prefix must be recognized");
		parser.parse(HEADER+":UMI 2:Y:0:AACC+GGTT extra");
		require(parser.chastityFail() && parser.pairnum()==1 && "UMI".equals(parser.index3()), "UMI/pair/chastity");
		require(parser.numBarcodes()==2 && parser.barcodeLetters()==8 && "extra".equals(parser.extra()), "Dual barcode");
		System.out.println("headers\tPASS\tordinary,shrunk,UMI");

		final CycleTracker tracker=new CycleTracker(), other=new CycleTracker();
		tracker.clear(); tracker.process();
		tracker.add(read("ACGT")); tracker.process();
		require(tracker.avg('G')==0.25f && tracker.max('G')==1f && tracker.averages[5]==20f, "Single-read literal cycles");
		other.add(read("GG")); tracker.add(other); tracker.process();
		require(tracker.avg('G')==0.5f && tracker.averages[5]==20f, "Merged counts must be three G among six bases");
		System.out.println("cycles-control\tPASS\tG=3/6,Q=20");
		tracker.clear(); tracker.process();
		final boolean finiteEmpty=finite(tracker);
		require(finiteEmpty==fixed, "Reset/process finite="+finiteEmpty+", fixed="+fixed);
		tracker.clear(); tracker.add(read("AA")); tracker.process();
		final boolean finiteReuse=finite(tracker);
		require(finiteReuse==fixed && tracker.avg('A')==1f, "Shorter-read reset reuse");
		System.out.println("cycles-reset\tPASS\tfinite-empty="+finiteEmpty+",finite-shorter="+finiteReuse);

		final MicroTile tile=new MicroTile();
		tile.goodKmerSum=0.5; tile.validKmerSum=1;
		final double actual=tile.goodKmerFraction(), expected=fixed ? 0.5 : 1d/3;
		require(Math.abs(actual-expected)<1e-12, "Good/valid fraction: "+actual);
		System.out.println("good-fraction\tPASS\t"+actual);
	}

	/** Captures existing rounded serialization; its preservation policy is deliberately not chosen here. */
	private static void dump(File out, boolean fixed) throws Exception{
		ByteFile.FORCE_MODE_BF1=true; ByteFile.FORCE_MODE_BF2=false;
		ReadWrite.USE_PIGZ=false; ReadWrite.USE_UNPIGZ=false;
		ReadWrite.USE_BGZIP=false; ReadWrite.USE_UNBGZIP=false;
		final FlowCell original=cell();
		final File full=new File(out, "three-reads.dump");
		TileDump.write(original, full.getPath(), false);
		final FlowCell loaded=TileDump.loadDump(full.getPath());
		final MicroTile tile=loaded.toList().get(0);
		require(tile.readCount==3 && tile.baseCount==12 && tile.alignedReadCount==0, "Roundtrip integer controls");
		final double before=original.toList().get(0).uniquePercent(), after=tile.uniquePercent();
		require(Math.abs(before-100d/3)<1e-10, "Independent initial uniqueness oracle");
		if(!fixed){require(after==0, "Historical rounded-counter loss changed: "+after);}
		System.out.println("roundtrip-capture\tOBSERVED\tbefore="+before+",after="+after);
		final ByteFile reader=ByteFile.makeByteFile(full.getPath(), false);
		final long rows;
		try{
			final FlowCell target=new FlowCell(31);
			final int version=TileDump.parseHeader(reader, target);
			rows=TileDump.loadTiles2(reader, target, version);
			require(target.toList().size()==1, "Exactly one parsed microtile");
		}finally{require(!reader.close(), "Valid dump reader status");}
		require(rows==(fixed ? 1 : 0), "Loader must return actual row count: "+rows);
		System.out.println("loader-v3-count\tPASS\t"+rows);

		final File legacy=new File(out, "one-row-v1.tsv");
		final String row="0\t0\t0\t499\t0\t499\t3\t0\t0\t0\t0\t0\t0\t0\t0\t33.3333\t20\t90\t2\t0\t23\t0\t0\t0\t0\t0\t3\t3\t3\t0\t0\t0\t0\t0\n";
		Files.write(legacy.toPath(), row.getBytes(StandardCharsets.US_ASCII));
		final ByteFile oldReader=ByteFile.makeByteFile(legacy.getPath(), false);
		final long oldRows;
		try{oldRows=TileDump.loadTiles1(oldReader, new FlowCell(31));}
		finally{require(!oldReader.close(), "Valid legacy reader status");}
		require(oldRows==(fixed ? 1 : 0), "Legacy loader row count: "+oldRows);
		System.out.println("loader-v1-count\tPASS\t"+oldRows);
		final File empty=new File(out, "empty-rows.tsv");
		Files.write(empty.toPath(), new byte[0]);
		for(int version=1; version<=3; version+=2){
			final ByteFile emptyReader=ByteFile.makeByteFile(empty.getPath(), false);
			try{
				final long count=version==1 ? TileDump.loadTiles1(emptyReader, new FlowCell(31))
						: TileDump.loadTiles2(emptyReader, new FlowCell(31), version);
				require(count==0, "Empty loader control");
			}finally{require(!emptyReader.close(), "Empty reader status");}
		}
		System.out.println("loader-empty\tPASS\t0,0");

		final File gzip=new File(out, "valid.dump.gz"), corrupt=new File(out, "bad-crc.dump.gz");
		try(GZIPOutputStream stream=new GZIPOutputStream(Files.newOutputStream(gzip.toPath()))){
			stream.write(Files.readAllBytes(full.toPath()));
		}
		require(TileDump.loadDump(gzip.getPath()).toList().get(0).readCount==3, "Valid gzip control");
		final byte[] damaged=Files.readAllBytes(gzip.toPath());
		damaged[damaged.length-8]^=1;// Alter the CRC trailer, preserving every compressed payload byte.
		Files.write(corrupt.toPath(), damaged);
		boolean threw=false;
		try{TileDump.loadDump(corrupt.getPath());}
		catch(RuntimeException e){threw=true; System.out.println("reader-cause\t"+e.getClass().getName());}
		require(threw==fixed, "CRC failure reporting: threw="+threw+", fixed="+fixed);
		System.out.println("reader-crc\tPASS\tthrew="+threw);
	}

	/** Three reads, twelve bases, and one missing kmer among three give 33.3333% in the dump. */
	private static FlowCell cell(){
		Tile.xSize=Tile.ySize=500;
		final FlowCell cell=new FlowCell(31);
		final MicroTile tile=cell.getMicroTile(0, 0, 0, 0);
		tile.readCount=3; tile.baseCount=12; tile.hits=2; tile.misses=1; tile.depthSum=6;
		tile.acgtn[0]=12;
		cell.calcStats();
		return cell;
	}

	private static Read read(String bases){
		final byte[] qualities=new byte[bases.length()];
		Arrays.fill(qualities, (byte)20);
		return new Read(bases.getBytes(StandardCharsets.US_ASCII), qualities, HEADER, 0);
	}

	private static void coordinates(ReadHeaderParser parser){
		require(parser.lane()==1 && parser.tile()==1101 && parser.xPos()==12 && parser.yPos()==34, "Literal coordinate oracle");
	}

	private static boolean finite(CycleTracker tracker){
		if(tracker.averages!=null){for(float value : tracker.averages){if(!Float.isFinite(value)){return false;}}}
		for(float[] cycle : tracker.cycleAverages){
			if(cycle!=null){for(float value : cycle){if(!Float.isFinite(value)){return false;}}}
		}
		return true;
	}

	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final String HEADER="M:1:F:1:1101:12:34";
}
