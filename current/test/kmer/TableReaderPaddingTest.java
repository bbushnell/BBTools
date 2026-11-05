package test.kmer;

import java.nio.charset.StandardCharsets;
import java.util.BitSet;

import kmer.HashArray1D;
import kmer.TableReader;
import stream.Read;

/**
 * Literal padding checks for TableReader's marked positions and masked-base count.
 * Each fixture contains one explicitly encoded reference key in seven tiny shards.
 * Expected strings specify the entire result independently of the masking code.
 * @author Shinobu
 */
public final class TableReaderPaddingTest{

	/** Runs one positive-padding case or an unchanged zero-padding/no-hit control. */
	public static void main(String[] args){
		if(args.length!=1 || !args[0].startsWith("case=")){
			throw new IllegalArgumentException("Expected one case=name argument.");
		}
		final String mode=args[0].substring(5);
		switch(mode){
			case "right-edge": checkMask(mode, 3, 0, 64, 1, "AAA", "NNN"); break;
			case "interior": checkMask(mode, 3, 0, 64, 1, "CCCAAACCC", "CCNNNNNCC"); break;
			case "short-edge": checkMask(mode, 5, 4, 256, 1, "AAAA", "NNNN"); break;
			case "zero-control": checkMask(mode, 3, 0, 64, 0, "CCCAAACCC", "CCCNNNCCC"); break;
			case "no-hit-control": checkMask(mode, 3, 0, 64, 1, "CCC", "CCC"); break;
			default: throw new IllegalArgumentException("Unknown padding case: "+mode);
		}
		System.out.println("PASS padding case="+mode);
	}

	/** Checks exact marked positions, output bases, qualities and returned count. */
	private static void checkMask(final String mode, final int k, final int mink,
			final long key, final int pad, final String input, final String expected){
		final HashArray1D[] tables=new HashArray1D[7];
		for(int i=0; i<tables.length; i++){tables[i]=new HashArray1D(new int[] {17}, -1L);}
		check(tables[(int)(key%7)].set(key, 1)==1, mode, "Fixture must create one reference key");
		final TableReader reader=new TableReader(k, mink, 0, 0, 0, false, false);
		reader.trimPad=pad;
		final Read read=read(input);
		final BitSet bits=reader.markBits(read, tables);
		check(bits==null || bits.length()<=input.length(), mode, "Marked bits extend beyond the read");
		int marked=0;
		for(int i=0; i<expected.length(); i++){
			final boolean mask=(expected.charAt(i)=='N');
			check((bits!=null && bits.get(i))==mask, mode, "Unexpected marked position "+i);
			if(mask){marked++;}
		}
		check((bits==null ? 0 : bits.cardinality())==marked, mode, "Marked cardinality differs");
		final int actual=reader.kMask(read, tables);
		check(actual==marked, mode, "Returned count differs from masked-base count");
		check(expected.equals(new String(read.bases, StandardCharsets.US_ASCII)), mode, "Masked bases differ");
		for(int i=0; i<expected.length(); i++){
			check(read.quality[i]==(expected.charAt(i)=='N' ? 0 : 40), mode, "Unexpected quality at "+i);
		}
		check(tables[(int)(key%7)].getValue(key)==1, mode, "Query changed the reference value");
	}

	/** Makes a fresh short read whose original bases are all defined and qualities are 40. */
	private static Read read(final String bases){
		final byte[] quality=new byte[bases.length()];
		for(int i=0; i<quality.length; i++){quality[i]=40;}
		return new Read(bases.getBytes(StandardCharsets.US_ASCII), quality, bases, 0);
	}

	/** Keeps exact-oracle failures visible with assertions enabled or disabled. */
	private static void check(final boolean condition, final String mode, final String message){
		if(!condition){throw new AssertionError("KMR113 case="+mode+": "+message);}
	}
}
