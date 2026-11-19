package test.prok;

import java.util.Arrays;
import java.util.Locale;

import cardinality.DynamicDemiLog;
import dna.AminoAcid;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlQueryKeys;
import shared.Shared;

/** Raw-key codec parity against fresh native sketches across the complete grid.
 * Isolated keys exercise the codec independently of bucket competition. Whole
 * sequences prove that per-bucket raw-key maxima equal native sketch maxima.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlQueryKeysTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaDdlQueryKeysTest.class, false);
		final int previous=DynamicDemiLog.exponentBits();
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown test parameter: "+arg);
			}
			check(parser.in1==null && parser.in2==null && parser.out1==null && parser.out2==null
				&& parser.maxReads<0 && parser.samplerate==1, "Codec parity uses complete synthetic sequences only");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Raw-key parity requires native assertions enabled");
			check(!DynamicDemiLog.blacklistExists(), "A process blacklist would change the native codec control");
			DynamicDemiLog.setExponent(4);
			check(NcrnaDdlQueryKeys.encodeKey(0)==4095 && NcrnaDdlQueryKeys.encodeKey(-1)==1
				&& NcrnaDdlQueryKeys.encodeKey(1)==65535
				&& NcrnaDdlQueryKeys.encodeKey(Long.MAX_VALUE)==4096
				&& NcrnaDdlQueryKeys.encodeKey(Long.MIN_VALUE)==4095,
				"Literal mixed-hash edges preserve native exponent4 NLZ masking, minimum1 and sixteen-bit keys");
			final byte[] sequence=dna(384);
			for(final int k:new int[]{5, 7, 9, 11, 13, 15, 19}){
				for(int width=32; width<=2048; width<<=1){
					final NcrnaDdlQueryKeys query=new NcrnaDdlQueryKeys(width, k);
					for(int strand=0; strand<2; strand++){
						final byte[] oriented=sequence.clone();
						if(strand==1){AminoAcid.reverseComplementBasesInPlace(oriented);}
						checkSequence(query, oriented, width, k, false);
						for(final int offset:new int[]{0, 39, 107}){
							checkSequence(query, Arrays.copyOfRange(oriented, offset, offset+k), width, k, true);
						}
						final byte[] interrupted=oriented.clone();
						Arrays.fill(interrupted, 93, 117, (byte)'N');
						for(int i=160; i<190; i++){interrupted[i]=(byte)Character.toLowerCase(interrupted[i]);}
						checkSequence(query, interrupted, width, k, false);
					}
					query.build(new byte[0]);
					check(query.packedKeys().length==0 && query.equalBuckets(new char[width])==0,
						"Replacing a query must discard all preceding raw keys; empty buckets never match");
				}
			}
			for(int k=7; k<=11; k++){
				for(int width=256; width<=2048; width<<=1){
					checkForwardMasked(sequence, width, k);
				}
			}
			pp.outstream.println("NCRNA_DDL_QUERY_KEYS_TEST_PASS");
		}finally{
			DynamicDemiLog.setExponent(previous);
			Shared.closeStream(pp.outstream);
		}
	}

	private static void checkSequence(final NcrnaDdlQueryKeys query, final byte[] bases,
			final int width, final int k, final boolean isolated){
		query.build(bases);
		final int[] raw=query.packedKeys();
		final char[] maximum=new char[width];
		for(final int packed:raw){
			final int bucket=packed>>>16;
			final char key=(char)packed;
			check(bucket<width && key!=0, "Packed query keys must preserve their native bucket and nonempty sixteen-bit value");
			if(key>maximum[bucket]){maximum[bucket]=key;}
		}
		final DynamicDemiLog nativeSketch=DynamicDemiLog.create(width, k, 12345, 0, false);
		nativeSketch.hash(bases, null);
		check(Arrays.equals(maximum, nativeSketch.maxArray()),
			"Raw-key maxima must equal a fresh native sketch for k="+k+", buckets="+width);
		if(isolated){check(raw.length==1, "An isolated unambiguous k-mer has exactly one native encoded key");}
		if(raw.length>0){
			final char[] profile=new char[width];
			profile[raw[0]>>>16]=(char)raw[0];
			check(query.equalBuckets(profile)==1, "One observed profile bucket contributes once despite repeated query occurrences");
		}
		if(width==32 && !isolated){
			int occupied=0;
			for(final char key:maximum){if(key!=0){occupied++;}}
			check(raw.length>occupied, "The multi-kmer fixture must retain competing raw keys below native bucket maxima");
		}
	}

	/** Shared insertion independently encodes the identical forward masked words. */
	private static void checkForwardMasked(final byte[] bases, final int width, final int k){
		final NcrnaDdlQueryKeys raw=new NcrnaDdlQueryKeys(width, k);
		final DynamicDemiLog nativeSketch=DynamicDemiLog.create(width, k, 12345, 0, false, false);
		raw.build(bases, nativeSketch, false, true);
		final char[] maximum=new char[width];
		for(final int packed:raw.packedKeys()){
			final int bucket=packed>>>16;
			final char key=(char)packed;
			check(bucket<width && key>0, "Forward masked keys need an in-range bucket and a nonempty native value");
			if(key>maximum[bucket]){maximum[bucket]=key;}
		}
		check(Arrays.equals(maximum, nativeSketch.maxArray()),
			"Prok-local encoding must match shared native insertion for forward masked k="+k+", buckets="+width);
	}

	private static byte[] dna(final int length){
		check(length>=200, "The parity fixture needs distinct isolated positions and a multi-key collision population");
		final byte[] bases=new byte[length];
		long state=729731;
		for(int i=0; i<bases.length; i++){
			state^=state<<13;
			state^=state>>>7;
			state^=state<<17;
			bases[i]=AminoAcid.numberToBase[(int)(state&3)];
		}
		return bases;
	}

	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
