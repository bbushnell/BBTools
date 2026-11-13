package test.barcode;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.Locale;

import barcode.Barcode;
import barcode.DemuxData;
import barcode.PCRMatrix;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/** Finite wire checks for DemuxData matrix type and Hamming-distance sum flags. */
public final class DemuxDataMatrixWireTest{

	/** Requires mode=baseline|fixed; no optional matrix construction is attempted. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, DemuxDataMatrixWireTest.class, false);
		try{
			String mode=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(!Parser.parseStatic(arg, key, value)){
					throw new IllegalArgumentException("Unknown argument: "+arg);
				}
			}
			check("baseline".equals(mode) || "fixed".equals(mode),
					"Set mode=baseline or mode=fixed.");

			int matched=0, lost=0;
			for(final int type : TYPES){
				for(final boolean hdistSum : SUMS){
					final ArrayList<byte[]> wire=wire(type, hdistSum);
					final DemuxData decoded=new DemuxData(wire);
					final boolean same=(decoded.type==type && decoded.hdistSum==hdistSum);
					System.out.println("case type="+type+" hdistSum="+hdistSum+
							" decodedType="+decoded.type+" decodedHdistSum="+
							decoded.hdistSum);
					if("baseline".equals(mode)){
						check(decoded.type==PCRMatrix.PROB_TYPE && !decoded.hdistSum,
								"Baseline should lose matrix flags to PROB/false: type="+type+", hdistSum="+
								hdistSum+", decodedType="+decoded.type+", decodedHdistSum="+
								decoded.hdistSum);
					}else{
						check(same, "Matrix flags did not round-trip: type="+type+
								", hdistSum="+hdistSum+", decodedType="+decoded.type+
								", decodedHdistSum="+decoded.hdistSum);
					}
					checkWire(mode, wire, type, hdistSum);
					if(same){matched++;}else{lost++;}
				}
			}
			checkOldWireDefault();
			System.out.println("DEMUXDATA_MATRIX_WIRE mode="+mode+
					" matched="+matched+" lost="+lost+" oldDefault=ok");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private DemuxDataMatrixWireTest(){}

	/** Encodes a minimal request with one observed barcode. */
	private static ArrayList<byte[]> wire(final int type, final boolean hdistSum){
		final DemuxData data=new DemuxData(4, 0, 0);
		data.type=type;
		data.hdistSum=hdistSum;
		data.expectedList=new LinkedHashSet<String>();
		data.expectedList.add("ACGT");
		final ArrayList<Barcode> counts=new ArrayList<Barcode>();
		counts.add(new Barcode("ACGT", 1));
		data.codeCounts=counts;
		return data.encode();
	}

	/** Requires absent baseline fields and nondefault fixed fields. */
	private static void checkWire(final String mode, final ArrayList<byte[]> chunks,
			final int type, final boolean hdistSum){
		final String wire=toString(chunks);
		final boolean hasType=wire.contains("\ntype=");
		final boolean hasHdistSum=wire.contains("\nhdistsum=");
		final boolean shouldHaveType=type!=PCRMatrix.PROB_TYPE;
		final boolean shouldHaveHdistSum=hdistSum;
		if("baseline".equals(mode)){
			check(!hasType, "Baseline unexpectedly emitted type:\n"+wire);
			check(!hasHdistSum, "Baseline unexpectedly emitted hdistSum:\n"+wire);
		}else{
			check(hasType==shouldHaveType, "Fixed wire type presence mismatch for "+
					type+":\n"+wire);
			check(hasHdistSum==shouldHaveHdistSum,
					"Fixed wire hdistSum presence mismatch for "+hdistSum+":\n"+wire);
			if(shouldHaveType){
				check(wire.contains("\ntype="+type+"\n"),
						"Fixed wire did not emit numeric matrix type "+type+":\n"+wire);
			}
			if(shouldHaveHdistSum){
				check(wire.contains("\nhdistsum=true\n"),
						"Fixed wire did not emit hdistSum=true:\n"+wire);
			}
		}
	}

	/** Missing fields from existing clients must keep the decode-side default. */
	private static void checkOldWireDefault(){
		final ArrayList<byte[]> chunks=oldWire(wire(PCRMatrix.PROB_TYPE, false));
		final DemuxData decoded=new DemuxData(chunks);
		check(decoded.type==PCRMatrix.PROB_TYPE,
				"Old wire matrix type default changed: "+decoded.type);
		check(!decoded.hdistSum, "Old wire hdistSum default changed.");
	}

	/** Simulates an existing client by deleting any fixed type/hdistSum lines. */
	private static ArrayList<byte[]> oldWire(final ArrayList<byte[]> chunks){
		final String[] lines=toString(chunks).split("\n");
		final StringBuilder sb=new StringBuilder();
		for(final String line : lines){
			final String key=line.split("=", 2)[0].toLowerCase(Locale.ROOT);
			if("type".equals(key) || "hdistsum".equals(key)){continue;}
			sb.append(line).append('\n');
		}
		final ArrayList<byte[]> old=new ArrayList<byte[]>(1);
		old.add(sb.toString().getBytes(StandardCharsets.UTF_8));
		return old;
	}

	/** Joins DemuxData chunks as UTF-8 text. */
	private static String toString(final ArrayList<byte[]> chunks){
		final StringBuilder sb=new StringBuilder();
		for(final byte[] chunk : chunks){
			sb.append(new String(chunk, StandardCharsets.UTF_8));
		}
		return sb.toString();
	}

	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final int[] TYPES={PCRMatrix.HDIST_TYPE, PCRMatrix.PROB_TYPE,
			PCRMatrix.TILE_TYPE};
	private static final boolean[] SUMS={false, true};
}
