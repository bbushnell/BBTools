package test.prot;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import parse.LineParser1;
import prot.FamilyCompactCodec;
import structures.ByteBuilder;

/** Exact ten-line compact family resource fixtures.
 * @author Yelan
 */
public final class FamilyCompactCodecTest{

	/** Verifies compact family header and record wire-format edge cases. */
	public static void main(String[] args){
		if(args.length!=0){throw new IllegalArgumentException("No compact-codec-test arguments");}
		final FamilyCompactCodec.Header header=header();
		final byte[] headerBytes=FamilyCompactCodec.appendHeader(header, new ByteBuilder()).toBytes();
		final byte[][] headerLines=lines(headerBytes);
		check(headerLines.length==10,"Header occupies exactly ten lines");
		check(FamilyCompactCodec.parseHeader(headerLines, 0, parser()).get("geometry").equals("k5.1"),
			"Packed key/value header roundtrips");
		final byte[][] padded=headerLines.clone();
		padded[8]=ascii("later=true");
		reject(()->FamilyCompactCodec.parseHeader(padded, 0, parser()));
		padded[8]=ascii("");
		padded[1]=ascii("view=SHORT\tview=TAILS");
		reject(()->FamilyCompactCodec.parseHeader(padded, 0, parser()));
		final java.util.LinkedHashMap<String, String> many=new java.util.LinkedHashMap<String, String>();
		for(int i=0; i<12; i++){many.put("k"+i, "v"+i);}
		check(FamilyCompactCodec.parseHeader(lines(FamilyCompactCodec.appendHeader(
			new FamilyCompactCodec.Header(many), new ByteBuilder()).toBytes()), 0, parser()).fields.size()==12,
			"Header supports more key/value pairs than physical metadata rows");

		final FamilyCompactCodec.Record record=record();
		final byte[] bytes=FamilyCompactCodec.appendRecord(record, new ByteBuilder()).toBytes();
		final byte[][] rows=lines(bytes);
		check(rows.length==10 && rows[5].length==1540, "Family occupies ten lines with a 1540-byte DDL row");
		final FamilyCompactCodec.Record copy=FamilyCompactCodec.parseRecord(rows, 0, parser());
		check(copy.familyId.equals("fam0") && copy.rank==7 && copy.stableId==123 &&
			copy.members==30 && copy.measured==30,"Family metadata roundtrips");
		check(Arrays.equals(copy.consensusRaw, ascii("ACDEFGHIK")) &&
			Arrays.equals(copy.ddl, record.ddl) && Arrays.equals(copy.monomers, record.monomers),
			"Compact payload arrays roundtrip");
		check(sameRawDoubles(copy.scalarBounds, record.scalarBounds) &&
			sameRawDouble(copy.lengthMedian, record.lengthMedian) &&
			sameRawDouble(copy.kdMean, record.kdMean) &&
			sameRawDouble(copy.meanGc, record.meanGc) &&
			sameRawDouble(copy.sampleSd, record.sampleSd), "Compact f64 fields preserve raw bits");
		check(Float.floatToRawIntBits(copy.monomers[0])==Float.floatToRawIntBits(-0.0f),
			"Monomer storage preserves exact f32 bits");
		check(Arrays.equals(bytes, FamilyCompactCodec.appendRecord(copy, new ByteBuilder()).toBytes()),
			"Compact family byte roundtrip");
		final byte[] emptyBytes=FamilyCompactCodec.appendRecord(emptyRecord(), new ByteBuilder()).toBytes();
		final FamilyCompactCodec.Record empty=FamilyCompactCodec.parseRecord(lines(emptyBytes), 0, parser());
		check(empty.members==0 && empty.consensusRaw.length==0 && empty.lengthMedian==0 &&
			empty.ddl[0]==0 && empty.monomers[0]==0,
			"Empty compact records retain inert payloads without rank imputation");
		final byte[][] badEmptyDdl=lines(emptyBytes);
		badEmptyDdl[5]=ascii("ddl\t000");
		reject(()->FamilyCompactCodec.parseRecord(badEmptyDdl, 0, parser()));
		final byte[][] badEmptyMono=lines(emptyBytes);
		badEmptyMono[6]=ascii("mono\t0");
		reject(()->FamilyCompactCodec.parseRecord(badEmptyMono, 0, parser()));
		final byte[][] badSingletonGc=lines(FamilyCompactCodec.appendRecord(singletonRecord(), new ByteBuilder()).toBytes());
		badSingletonGc[3]=new ByteBuilder("gc\t1\t").appendDoubleA48(0.5).tab()
			.appendDoubleA48(0.125).toBytes();
		reject(()->FamilyCompactCodec.parseRecord(badSingletonGc, 0, parser()));
		final byte[][] badConsensus=rows.clone();
		badConsensus[4]=ascii("seq\tACD*");
		reject(()->FamilyCompactCodec.parseRecord(badConsensus, 0, parser()));
		final FamilyCompactCodec.Record badScalar=record();
		badScalar.scalarBounds[0]=Double.POSITIVE_INFINITY;
		reject(()->FamilyCompactCodec.appendRecord(badScalar, new ByteBuilder()));
		reject(()->FamilyCompactCodec.appendRecord(record(scalars(), Double.NaN, 0.25, 0.5, 0.125, monomers()),
			new ByteBuilder()));
		reject(()->FamilyCompactCodec.appendRecord(record(scalars(), 50.5, Double.POSITIVE_INFINITY,
			0.5, 0.125, monomers()), new ByteBuilder()));
		reject(()->FamilyCompactCodec.appendRecord(record(scalars(), 50.5, 0.25, Double.NEGATIVE_INFINITY,
			0.125, monomers()), new ByteBuilder()));
		reject(()->FamilyCompactCodec.appendRecord(record(scalars(), 50.5, 0.25, 0.5, Double.NaN,
			monomers()), new ByteBuilder()));
		final float[] badMonomers=monomers();
		badMonomers[1]=Float.NaN;
		reject(()->FamilyCompactCodec.appendRecord(record(scalars(), 50.5, 0.25, 0.5, 0.125, badMonomers),
			new ByteBuilder()));

		final char[] ddl=new char[FamilyCompactCodec.DDL_WIDTH];
		ddl[0]=0;
		ddl[1]=1;
		ddl[2]=4096;
		ddl[3]=65535;
		ddl[511]=64;
		final String ddlLine=FamilyCompactCodec.appendFixedDdl(ddl, new ByteBuilder()).toString();
		check(ddlLine.startsWith("ddl\t000001100?oo"), "Fixed DDL uses padded radix-64 triplets");
		final byte[][] ddlRows=rows.clone();
		ddlRows[5]=ascii(ddlLine);
		check(FamilyCompactCodec.parseRecord(ddlRows, 0, parser()).ddl[3]==65535,
			"Fixed DDL parses high uint16 values");
		for(int badPos : new int[]{4, 5, 6}){
			final byte[][] bad=rows.clone();
			final byte[] line=bad[5].clone();
			line[badPos]=(byte)(badPos==4 ? '@' : 'p');
			bad[5]=line;
			reject(()->FamilyCompactCodec.parseRecord(bad, 0, parser()));
		}
		final byte[][] shortDdl=rows.clone();
		shortDdl[5]=ascii(new String(rows[5], 0, rows[5].length-1, StandardCharsets.US_ASCII));
		reject(()->FamilyCompactCodec.parseRecord(shortDdl, 0, parser()));

		for(int line=0; line<10; line++){
			final byte[][] bad=rows.clone();
			bad[line]=ascii(line==0 ? "seq\tACD" : line==7 ? "pos\textra" : "bad");
			reject(()->FamilyCompactCodec.parseRecord(bad, 0, parser()));
		}
		final byte[][] wrongRank=rows.clone();
		wrongRank[2]=new ByteBuilder("rank\t").appendDoubleA48(39.5).tab().appendDoubleA48(0.25).toBytes();
		reject(()->FamilyCompactCodec.parseRecord(wrongRank, 0, parser()));
		final byte[][] wrongGc=rows.clone();
		wrongGc[3]=new ByteBuilder("gc\t29\t").appendDoubleA48(0.5).tab().appendDoubleA48(0.125).toBytes();
		reject(()->FamilyCompactCodec.parseRecord(wrongGc, 0, parser()));
		final byte[][] shortHeader=new byte[9][];
		reject(()->FamilyCompactCodec.parseHeader(shortHeader, 0, parser()));
		System.out.println("FAMILY_COMPACT_CODEC_TEST_PASS checks="+checks);
	}

	/** Returns a valid sparse compact header with spare blank rows. */
	private static FamilyCompactCodec.Header header(){
		final byte[][] lines=new byte[10][];
		lines[0]=ascii("#bbgap_family_compact\t1");
		lines[1]=ascii("view=SHORT\tcensus=1");
		lines[2]=ascii("k=6\tseed=0\twidth=512");
		lines[3]=ascii("geometry=k5.1");
		for(int i=4; i<lines.length; i++){lines[i]=ascii("");}
		return FamilyCompactCodec.parseHeader(lines, 0, parser());
	}

	/** Returns a live compact record with every payload class populated. */
	private static FamilyCompactCodec.Record record(){
		final char[] ddl=new char[FamilyCompactCodec.DDL_WIDTH];
		for(int i=0; i<ddl.length; i++){ddl[i]=(char)(i*127);}
		ddl[511]=65535;
		return new FamilyCompactCodec.Record("fam0", 7, 123, 30, 270, 266,
			0.2f, 30f, -0.25f, 0.72f, 50f, 0.14f, 31, FamilyCompactCodec.AlignmentCalibration.ENABLED,
			scalars(), 50.5, 0.1, 30, 0.3, 0.2, ascii("ACDEFGHIK"), ddl, monomers());
	}

	/** Returns a live compact record with caller-supplied sidecar payloads. */
	private static FamilyCompactCodec.Record record(double[] scalars, double lengthMedian, double kdMean,
			double meanGc, double sampleSd, float[] monomers){
		final char[] ddl=new char[FamilyCompactCodec.DDL_WIDTH];
		for(int i=0; i<ddl.length; i++){ddl[i]=(char)(i*127);}
		return new FamilyCompactCodec.Record("fam0", 7, 123, 30, 270, 266,
			0.2f, 30f, -0.25f, 0.72f, 50f, 0.14f, 31, FamilyCompactCodec.AlignmentCalibration.ENABLED,
			scalars, lengthMedian, kdMean, 30, meanGc, sampleSd, ascii("ACDEFGHIK"), ddl, monomers);
	}

	private static double[] scalars(){
		return new double[]{0.1, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1, 9.1};
	}

	/** Returns exact f32 monomer payloads, including negative zero. */
	private static float[] monomers(){
		final float[] monomers=new float[FamilyCompactCodec.MONOMER_WIDTH];
		monomers[0]=-0.0f;
		monomers[1]=0.5f;
		monomers[2]=Float.intBitsToFloat(1);
		return monomers;
	}

	/** Returns a zero-member compact record with inert search payloads. */
	private static FamilyCompactCodec.Record emptyRecord(){
		return new FamilyCompactCodec.Record("empty", 8, 124, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			FamilyCompactCodec.AlignmentCalibration.UNAVAILABLE, new double[10], 0, 0, 0, 0, 0, ascii(""),
			new char[FamilyCompactCodec.DDL_WIDTH], new float[FamilyCompactCodec.MONOMER_WIDTH]);
	}

	/** Returns a one-member record for the zero-SD GC boundary. */
	private static FamilyCompactCodec.Record singletonRecord(){
		return new FamilyCompactCodec.Record("singleton", 9, 125, 1, 80, 76, 0, 0, 0, 0, 0, 0, 0,
			FamilyCompactCodec.AlignmentCalibration.UNAVAILABLE, new double[10], 80, 0, 1, 0.5, 0, ascii("ACD"),
			new char[FamilyCompactCodec.DDL_WIDTH], new float[FamilyCompactCodec.MONOMER_WIDTH]);
	}

	private static byte[][] lines(byte[] text){
		final String[] lines=new String(text, StandardCharsets.US_ASCII).split("\n", -1);
		final byte[][] out=new byte[lines.length-1][];
		for(int i=0; i<out.length; i++){out[i]=ascii(lines[i]);}
		return out;
	}

	private static boolean sameRawDoubles(double[] a, double[] b){
		if(a.length!=b.length){return false;}
		for(int i=0; i<a.length; i++){
			if(!sameRawDouble(a[i], b[i])){return false;}
		}
		return true;
	}

	private static boolean sameRawDouble(double a, double b){
		return Double.doubleToRawLongBits(a)==Double.doubleToRawLongBits(b);
	}

	private static LineParser1 parser(){return new LineParser1('\t');}
	private static byte[] ascii(String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static void check(boolean condition, String message){
		checks++;
		if(!condition){throw new AssertionError(message);}
	}
	private static void reject(Runnable r){
		checks++;
		try{r.run();}
		catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Expected compact codec rejection");
	}
	private static int checks;
}
