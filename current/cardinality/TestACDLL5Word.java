package cardinality;

/**
 * Simulate one word (6 buckets) on both CDLL5 and ACDLL5, hash by hash,
 * and report where they diverge.
 */
public class TestACDLL5Word{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		// Use 6 buckets so everything fits in one word
		final CompressedDynamicLogLog5 cdll5=new CompressedDynamicLogLog5(6, 31, 42, 0);
		final ArithmeticCompressedDynamicLogLog5 acdll5=new ArithmeticCompressedDynamicLogLog5(6, 31, 42, 0);

		int divergences=0;
		final int totalHashes=50000;

		for(long i=0; i<totalHashes; i++){
			cdll5.hashAndStore(i);
			acdll5.hashAndStore(i);

			// Compare all 6 buckets
			boolean mismatch=false;
			for(int b=0; b<6; b++){
				final int cReg=cdll5.readBucket(b);
				final int cExp=(cReg>>>2)&0x7;
				final int cHist=cReg&0x3;
				final int aExp=acdll5.readExponent(b);
				final int aHist=acdll5.readHistory(b);
				if(cExp!=aExp || cHist!=aHist){
					mismatch=true;
				}
			}

			if(mismatch){
				divergences++;
				if(divergences<=30){
					System.out.printf("DIVERGE after hash %d (minZeros: C=%d A=%d):%n",
						i, cdll5.getMinZeros(), acdll5.getMinZeros());
					for(int b=0; b<6; b++){
						final int cReg=cdll5.readBucket(b);
						final int cExp=(cReg>>>2)&0x7;
						final int cHist=cReg&0x3;
						final int aExp=acdll5.readExponent(b);
						final int aHist=acdll5.readHistory(b);
						final String flag=(cExp!=aExp || cHist!=aHist) ? " <<<" : "";
						System.out.printf("  b%d: CDLL5(exp=%d,hist=%d) ACDLL5(exp=%d,hist=%d)%s%n",
							b, cExp, cHist, aExp, aHist, flag);
					}
					System.out.println();
				}
			}
		}
		System.out.printf("Total divergences: %d / %d hashes%n", divergences, totalHashes);

		// Final state
		System.out.println("\nFinal state:");
		System.out.printf("minZeros: CDLL5=%d ACDLL5=%d%n", cdll5.getMinZeros(), acdll5.getMinZeros());
		for(int b=0; b<6; b++){
			final int cReg=cdll5.readBucket(b);
			final int cExp=(cReg>>>2)&0x7;
			final int cHist=cReg&0x3;
			final int aExp=acdll5.readExponent(b);
			final int aHist=acdll5.readHistory(b);
			System.out.printf("  b%d: CDLL5(exp=%d,hist=%d) ACDLL5(exp=%d,hist=%d)%n",
				b, cExp, cHist, aExp, aHist);
		}
	}

}
