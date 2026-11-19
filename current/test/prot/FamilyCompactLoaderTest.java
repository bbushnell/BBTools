package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import prot.FamilyCompactCodec;
import prot.FamilyCompactLoader;
import prot.FamilyCompactLoader.Loaded;
import prot.FamilyReference;
import structures.ByteBuilder;

/** Ten-line compact-family loader framing and census fixtures.
 * @author Yelan
 */
public final class FamilyCompactLoaderTest{

	/** Verifies serial/parallel compact loading and malformed compact inputs. */
	public static void main(String[] args)throws Exception{
		if(args.length!=1 || !args[0].startsWith("fixture=")){
			throw new IllegalArgumentException("Require fixture=NEW_DIR on Dori");
		}
		final Path root=Paths.get(args[0].substring(8)).toAbsolutePath().normalize();
		Files.createDirectory(root);
		for(int size : new int[]{0, 1, 18, 19, 20, 21, 39, 40, 41, 81}){
			final FamilyCompactCodec.Record[] expected=records(size);
			final Path path=root.resolve("compact_"+size+".txt");
			write(path, encode(expected));
			final Loaded serial=FamilyCompactLoader.load(path.toString(), 1, size);
			verify(serial, expected);
			for(int threads : new int[]{2, 4, 8}){
				final Loaded parallel=FamilyCompactLoader.load(path.toString(), threads, size);
				verify(parallel, expected);
				check(Arrays.equals(stableBytes(serial), stableBytes(parallel)),
					"Serial and parallel compact loads agree at "+size+" families");
			}
		}

		final FamilyCompactCodec.Record[] expected=records(41);
		final byte[] good=encode(expected);
		final Path gzip=root.resolve("compact.gz");
		try(GZIPOutputStream out=new GZIPOutputStream(Files.newOutputStream(gzip))){out.write(good);}
		verify(FamilyCompactLoader.load(gzip.toString(), 4, 41), expected);

		final List<String> lines=lines(good);
		mutate(root, lines, "duplicate_id", 20, 2, expected[0].familyId, 41);
		mutate(root, lines, "duplicate_rank", 20, 1, Integer.toString(expected[0].rank), 41);
		mutate(root, lines, "duplicate_stable", 20, 3, Integer.toString(expected[0].stableId), 41);
		mutate(root, lines, "bad_census", 1, 0, "families=42", 41);
		reject(root, "header_limit", lines, 40);
		final ArrayList<String> badReserved=new ArrayList<String>(lines);
		badReserved.set(17, "pos\tunexpected");
		reject(root, "bad_reserved", badReserved, 41);
		for(int missing=1; missing<=10; missing++){
			reject(root, "truncated_"+missing, lines.subList(0, lines.size()-missing), 41);
		}
		reject(root, "absent_header", lines.subList(10, lines.size()), 41);
		reject(root, "empty_file", new ArrayList<String>(), 41);
		System.out.println("FAMILY_COMPACT_LOADER_TEST_PASS checks="+checks);
	}

	/** Verifies loaded physical and live-family census state. */
	private static void verify(Loaded loaded, FamilyCompactCodec.Record[] expected){
		check(loaded.all.length==expected.length && loaded.allStableIds.length==expected.length,
			"Complete compact census retained");
		check(loaded.blocks==1+expected.length/FamilyCompactCodec.FAMILIES_PER_BATCH,
			"Compact 200-line block count");
		int live=0;
		long members=0, residues=0, words=0;
		for(int i=0; i<expected.length; i++){
			final FamilyCompactCodec.Record e=expected[i];
			check(loaded.all[i].id.equals(e.familyId) && loaded.all[i].rank==e.rank,
				"Physical family identity retained");
			check(loaded.allStableIds[i]==e.stableId,"Physical stable ID retained");
			verify(loaded.all[i], e);
			members+=e.members;
			residues+=e.residues;
			words+=e.validWords;
			if(e.members>0){
				check(loaded.live[live]==loaded.all[i],"Live compact view reuses physical objects");
				check(loaded.liveStableIds[live]==e.stableId &&
					raw(loaded.lengthMedians[live])==raw(e.lengthMedian) &&
					raw(loaded.kdMeans[live])==raw(e.kdMean) &&
					loaded.gcMeasured[live]==e.measured &&
					raw(loaded.gcMeans[live])==raw(e.meanGc) &&
					raw(loaded.gcSampleSds[live])==raw(e.sampleSd),
					"Compact sidecar metadata follows live order without changing f64 bits");
				live++;
			}else{
				check(loaded.all[i].ddl==null && loaded.all[i].spectrum==null,
					"Empty compact rows stay nonsearchable");
			}
		}
		check(loaded.live.length==live && loaded.members==members &&
			loaded.residues==residues && loaded.validWords==words,"Compact totals are conserved");
	}

	/** Verifies one loaded compact family against its wire-format record. */
	private static void verify(FamilyReference family, FamilyCompactCodec.Record expected){
		check(family.members==expected.members && family.residues==expected.residues &&
			family.validWords==expected.validWords,"Compact family counts retained");
		check(raw(family.cosineLower)==raw(expected.cosineLower) &&
			raw(family.identityLower)==raw(expected.identityLower) &&
			raw(family.blosumLower)==raw(expected.blosumLower) &&
			raw(family.cosineAccept)==raw(expected.cosineAccept) &&
			raw(family.identityUpper)==raw(expected.identityUpper) &&
			raw(family.blosumFinal)==raw(expected.blosumFinal),
			"Compact gate f32 payload retained without changing bits");
		if(expected.members==0){
			check(family.scalarMask==0 && rawEquals(family.scalarBounds, new double[10]) &&
				family.consensus==null && family.consensusRaw==null &&
				family.ddl==null && family.spectrum==null,"Empty compact payload stays inert");
		}else{
			check(family.scalarMask==expected.scalarMask &&
					rawEquals(family.scalarBounds, expected.scalarBounds),
				"Compact scalar f64 payload retained without changing bits");
			check(family.alignmentCalibration.label.equals(expected.alignmentCalibration.label),
				"Compact alignment calibration state retained");
			check(Arrays.equals(family.consensusRaw, expected.consensusRaw),
				"Compact consensus raw bytes retained");
			check(Arrays.equals(family.ddl, expected.ddl), "Compact DDL registers retained");
			check(rawEquals(family.spectrum, expected.monomers),
				"Compact monomer f32 payload retained without changing bits");
		}
	}

	/** Returns a compact record set with live, empty and unsorted physical rows. */
	private static FamilyCompactCodec.Record[] records(int size){
		final FamilyCompactCodec.Record[] records=new FamilyCompactCodec.Record[size];
		for(int i=0; i<size; i++){
			final String id="family_"+i;
			final int rank=size-1-i, stableId=1000+i;
			if(i%7==6){
				records[i]=new FamilyCompactCodec.Record(id, rank, stableId, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
					FamilyCompactCodec.AlignmentCalibration.UNAVAILABLE, new double[10], 0, 0, 0, 0, 0,
					new byte[0], new char[FamilyCompactCodec.DDL_WIDTH],
					new float[FamilyCompactCodec.MONOMER_WIDTH]);
			}else{
				final char[] ddl=new char[FamilyCompactCodec.DDL_WIDTH];
				ddl[i%FamilyCompactCodec.DDL_WIDTH]=(char)(i+1);
				final float[] monomers=new float[FamilyCompactCodec.MONOMER_WIDTH];
				monomers[i%FamilyCompactCodec.MONOMER_WIDTH]=1;
				final long members=i+1, residues=members*200, validWords=residues-4;
				final double[] scalarBounds=new double[10];
				final int scalarMask=members>=30 ? 0x09 : 0;
				if(scalarMask!=0){
					scalarBounds[0]=Double.longBitsToDouble(0x3fc0000000000001L);
					scalarBounds[1]=Double.longBitsToDouble(0x3fe0000000000001L);
					scalarBounds[6]=Double.longBitsToDouble(0x3fd0000000000001L);
					scalarBounds[7]=Double.longBitsToDouble(0x3fe8000000000001L);
				}
				records[i]=new FamilyCompactCodec.Record(id, rank, stableId, members, residues, validWords,
					0.1f, 30f, -0.2f, 0.75f, 90f, 0.144f, scalarMask,
					FamilyCompactCodec.AlignmentCalibration.ENABLED, scalarBounds,
					50+(i%10)*0.5, 0.25, members, 0.5, members==1 ? 0 : 0.125,
					"ACDEFGHIK".getBytes(StandardCharsets.US_ASCII), ddl, monomers);
			}
		}
		return records;
	}

	/** Encodes records behind a header with matching census fields. */
	private static byte[] encode(FamilyCompactCodec.Record[] records){
		int live=0;
		for(FamilyCompactCodec.Record record : records){if(record.members>0){live++;}}
		final LinkedHashMap<String, String> fields=new LinkedHashMap<String, String>();
		fields.put("families", Integer.toString(records.length));
		fields.put("live", Integer.toString(live));
		fields.put("empty", Integer.toString(records.length-live));
		final ByteBuilder out=FamilyCompactCodec.appendHeader(new FamilyCompactCodec.Header(fields),
			new ByteBuilder());
		for(FamilyCompactCodec.Record record : records){FamilyCompactCodec.appendRecord(record, out);}
		return out.toBytes();
	}

	/** Emits loader-order stable IDs for serial/parallel equality checks. */
	private static byte[] stableBytes(Loaded loaded){
		final ByteBuilder out=new ByteBuilder();
		for(int stableId : loaded.allStableIds){out.append(stableId).nl();}
		for(int stableId : loaded.liveStableIds){out.append(stableId).nl();}
		return out.toBytes();
	}

	/** Mutates one compact cell and verifies the loader rejects the result. */
	private static void mutate(Path root, List<String> original, String name, int row, int term,
			String value, int maxFamilies)throws Exception{
		final ArrayList<String> changed=new ArrayList<String>(original);
		cell(changed, row, term, value);
		reject(root, name, changed, maxFamilies);
	}

	/** Replaces one tab-delimited cell in a generated compact row. */
	private static void cell(List<String> lines, int row, int term, String value){
		final String line=lines.get(row);
		int start=0;
		for(int i=0; i<term; i++){
			start=line.indexOf('\t', start)+1;
			check(start>0, "Fixture field exists");
		}
		int end=line.indexOf('\t', start);
		if(end<0){end=line.length();}
		lines.set(row, line.substring(0, start)+value+line.substring(end));
	}

	/** Verifies malformed compact text fails under serial and parallel loaders. */
	private static void reject(Path root, String name, List<String> lines, int maxFamilies)throws Exception{
		final Path path=root.resolve(name+".txt");
		write(path, (lines.isEmpty() ? "" : String.join("\n", lines)+"\n").getBytes(StandardCharsets.US_ASCII));
		for(int threads : new int[]{1, 4}){
			boolean failed=false;
			try{FamilyCompactLoader.load(path.toString(), threads, maxFamilies);}
			catch(IllegalArgumentException expected){
				failed=true;
				check(expected.getMessage()!=null && !expected.getMessage().isEmpty(),
					"Compact loader failures are diagnostic");
			}
			check(failed,"Malformed compact reference accepted");
		}
	}

	/** Splits generated ASCII compact text into complete lines. */
	private static List<String> lines(byte[] bytes){
		final ArrayList<String> out=new ArrayList<String>();
		int start=0;
		for(int i=0; i<bytes.length; i++){
			if(bytes[i]=='\n'){
				out.add(new String(bytes, start, i-start, StandardCharsets.US_ASCII));
				start=i+1;
			}
		}
		check(start==bytes.length,"Generated compact fixture has complete lines");
		return out;
	}

	private static void write(Path path, byte[] bytes)throws Exception{Files.write(path, bytes);}
	private static int raw(float value){return Float.floatToRawIntBits(value);}
	private static long raw(double value){return Double.doubleToRawLongBits(value);}
	private static boolean rawEquals(float[] a, float[] b){
		if(a.length!=b.length){return false;}
		for(int i=0; i<a.length; i++){if(raw(a[i])!=raw(b[i])){return false;}}
		return true;
	}
	private static boolean rawEquals(double[] a, double[] b){
		if(a.length!=b.length){return false;}
		for(int i=0; i<a.length; i++){if(raw(a[i])!=raw(b[i])){return false;}}
		return true;
	}
	private static void check(boolean condition, String message){
		checks++;
		if(!condition){throw new AssertionError(message);}
	}
	private static int checks;
}
