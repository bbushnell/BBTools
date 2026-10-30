package simd;

import java.util.Arrays;

import dna.AminoAcid;
import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.LongVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorShuffle;
import jdk.incubator.vector.VectorSpecies;
import shared.Random;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntList;

/**
 * Holds SIMD methods using dual-width strategy (128-bit + 64-bit).
 * This approach minimizes scalar tail operations across all array sizes.
 * @author Brian Bushnell, Isla
 * @author Sucrose - 128-bit port
 * @date Nov 7, 2025
 */
final class SIMDByte128{

	/** Runs a reverse-complement microbenchmark.
	 * @param args Ignored
	 */
	public static void main(String[] args){
		final int iterations=1000000000;
		final int k=31;
		final long mask=(1L<<(2*k))-1; // Mask for 31 bases (62 bits)

		// Generate random k-mers
		final int len=8192;
		final int mod=len-1;
		long[] kmers=new long[len];
		shared.Random randy=shared.Shared.random(12345);
		for(int i=0; i<kmers.length; i++){
			kmers[i]=randy.nextLong()&mask;
		}

		// Warmup
		for(int i=0; i<10000; i++){
			AminoAcid.reverseComplementBinaryFast(kmers[i&mod], k);
			reverseComplementBinaryFastBitwise2(kmers[i&mod], k);
		}

		// Test correctness on all test k-mers
		for(int i=0; i<kmers.length; i++){
			long result1=AminoAcid.reverseComplementBinaryFast(kmers[i], k);
			long result2=reverseComplementBinaryFastBitwise2(kmers[i], k);
			if(result1!=result2){
				System.err.println("Results don't match for kmer "+i+"!");
				System.err.println("Input:    "+kmers[i]);
				System.err.println("Original: "+result1);
				System.err.println("Bitwise:  "+result2);
				return;
			}
		}
		System.err.println("Correctness check passed on "+kmers.length+" k-mers");

		// Benchmark original
		long sum1=0;
		long t1=System.nanoTime();
		for(int i=0; i<iterations; i++){
			sum1+=AminoAcid.reverseComplementBinaryFast(kmers[i&mod], k);
		}
		long t2=System.nanoTime();
		System.err.println("Original: "+(t2-t1)/1000000+"ms, checksum="+sum1);

		// Benchmark bitwise
		long sum2=0;
		long t3=System.nanoTime();
		for(int i=0; i<iterations; i++){
			sum2+=reverseComplementBinaryFastBitwise2(kmers[i&mod], k);
		}
		long t4=System.nanoTime();
		System.err.println("Bitwise:  "+(t4-t3)/1000000+"ms, checksum="+sum2);

		double speedup=(double)(t2-t1)/(t4-t3);
		System.err.println("Speedup: "+String.format("%.2fx", speedup));
	}

//	public static void main(String[] args){
//		final int iterations=1000000000;
//		final int k=15;
//		final int mask=(1<<(2*k))-1; // Mask for 15 bases (30 bits)
//
//		// Generate random k-mers
//		int[] kmers=new int[10000];
//		shared.Random randy=shared.Shared.random(12345);
//		for(int i=0; i<kmers.length; i++){
//			kmers[i]=randy.nextInt()&mask;
//		}
//
//		// Warmup
//		for(int i=0; i<10000; i++){
//			AminoAcid.reverseComplementBinaryFast(kmers[i%kmers.length], k);
//			reverseComplementBinaryFastBitwise(kmers[i%kmers.length], k);
//		}
//
//		// Test correctness on all test k-mers
//		for(int i=0; i<kmers.length; i++){
//			int result1=(int)AminoAcid.reverseComplementBinaryFast(kmers[i], k);
//			int result2=reverseComplementBinaryFastBitwise(kmers[i], k);
//			if(result1!=result2){
//				System.err.println("Results don't match for kmer "+i+"!");
//				System.err.println("Input:    "+kmers[i]);
//				System.err.println("Original: "+result1);
//				System.err.println("Bitwise:  "+result2);
//				return;
//			}
//		}
//		System.err.println("Correctness check passed on "+kmers.length+" k-mers");
//
//		// Benchmark original
//		long sum1=0;
//		long t1=System.nanoTime();
//		for(int i=0; i<iterations; i++){
//			sum1+=AminoAcid.reverseComplementBinaryFast(kmers[i%kmers.length], k);
//		}
//		long t2=System.nanoTime();
//		System.err.println("Original: "+(t2-t1)/1000000+"ms, checksum="+sum1);
//
//		// Benchmark bitwise
//		long sum2=0;
//		long t3=System.nanoTime();
//		for(int i=0; i<iterations; i++){
//			sum2+=reverseComplementBinaryFastBitwise(kmers[i%kmers.length], k);
//		}
//		long t4=System.nanoTime();
//		System.err.println("Bitwise:  "+(t4-t3)/1000000+"ms, checksum="+sum2);
//
//		double speedup=(double)(t2-t1)/(t4-t3);
//		System.err.println("Speedup: "+String.format("%.2fx", speedup));
//	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	/** Package-private helper; operations are static. */
	SIMDByte128(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Builds a 16-lane lookup table for ASCII nucleotide predicates.
	 * @param half selected half of the folded ASCII low 5-bit index space
	 * @return Byte vector holding accepted uppercase bases at their low-4-bit indexes */
	private static ByteVector makeNucleotideLookup(final int half){
		final byte[] lookup=new byte[BWIDTH128];
		//0xFF has bit 5 set, so it cannot equal any byte folded with &~32.
		Arrays.fill(lookup, (byte)0xFF);
		final byte[] valid={'A', 'B', 'C', 'D', 'G', 'H', 'K', 'M', 'N',
				'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y'};
		for(final byte b : valid){
			if((b&16)==half){lookup[b&15]=b;}
		}
		return ByteVector.fromArray(BSPECIES128, lookup, 0);
	}

	/** Returns number of matches.
	 * @param s1 First sequence
	 * @param s2 Second sequence
	 * @param a1 First index in s1
	 * @param b1 Last index in s1
	 * @param a2 First index in s2
	 * @param b2 Last index in s2
	 * @return Number of equal positions
	 */
	static final int countMatches(final byte[] s1, final byte[] s2, int a1, int b1, int a2, int b2){
		int i=a1, j=a2;
		int matches=0;

		{//128-bit loop
			for(; j<=b2-BWIDTH128+1; i+=BWIDTH128, j+=BWIDTH128){
				ByteVector v1=ByteVector.fromArray(BSPECIES128, s1, i);
				ByteVector v2=ByteVector.fromArray(BSPECIES128, s2, j);
				VectorMask<Byte> x=v1.eq(v2);
				matches+=x.trueCount();
			}
		}

		{//64-bit loop
			for(; j<=b2-BWIDTH64+1; i+=BWIDTH64, j+=BWIDTH64){
				ByteVector v1=ByteVector.fromArray(BSPECIES64, s1, i);
				ByteVector v2=ByteVector.fromArray(BSPECIES64, s2, j);
				VectorMask<Byte> x=v1.eq(v2);
				matches+=x.trueCount();
			}
		}

		//Scalar tail
		for(; j<=b2; i++, j++){
			final byte x=s1[i], y=s2[j];
			final int m=((x==y) ? 1 : 0);
			matches+=m;
		}
		return matches;
	}

	/** Returns index of symbol.
	 * @param a Array to search
	 * @param symbol Byte to find
	 * @param from First index, inclusive
	 * @param to Last index, exclusive
	 * @return First matching index, or {@code to}
	 */
	static final int find(final byte[] a, final byte symbol, final int from, final int to){
		int pos=from;

		{//128-bit loop
			for(; pos<=to-BWIDTH128; pos+=BWIDTH128){
				ByteVector v=ByteVector.fromArray(BSPECIES128, a, pos);
				VectorMask<Byte> x=v.eq(symbol);
				int t=x.firstTrue();
				if(t<BWIDTH128){ return pos+t; }
			}
		}

		{//64-bit loop
			for(; pos<=to-BWIDTH64; pos+=BWIDTH64){
				ByteVector v=ByteVector.fromArray(BSPECIES64, a, pos);
				VectorMask<Byte> x=v.eq(symbol);
				int t=x.firstTrue();
				if(t<BWIDTH64){ return pos+t; }
			}
		}

		//Scalar tail
		while(pos<to && a[pos]!=symbol){ pos++; }
		return pos;
	}

	/**
	 * Sums the array.
	 * Widens bytes to int before reducing to avoid signed byte overflow.
	 * ByteVector.reduceLanesToLong(ADD) uses byte-precision arithmetic which
	 * wraps at +/-127, producing incorrect results for quality sums.
	 * @param a A byte array.
	 * @param from First index, inclusive
	 * @param to Last index, inclusive
	 * @return The sum.
	 */
	static final long sum(final byte[] a, final int from, final int to){
		int i=from;
		long c=0;

		{//128-bit loop: widen 16 bytes to 4 IntVectors of 4 ints each.
			for(; i<=to-BWIDTH128+1; i+=BWIDTH128){
				ByteVector va=ByteVector.fromArray(BSPECIES128, a, i);
				c+=((IntVector)va.convertShape(VectorOperators.B2I, ISPECIES128, 0)).reduceLanesToLong(VectorOperators.ADD);
				c+=((IntVector)va.convertShape(VectorOperators.B2I, ISPECIES128, 1)).reduceLanesToLong(VectorOperators.ADD);
				c+=((IntVector)va.convertShape(VectorOperators.B2I, ISPECIES128, 2)).reduceLanesToLong(VectorOperators.ADD);
				c+=((IntVector)va.convertShape(VectorOperators.B2I, ISPECIES128, 3)).reduceLanesToLong(VectorOperators.ADD);
			}
		}

		{//64-bit loop: widen 8 bytes to 2 IntVectors of 4 ints each.
			for(; i<=to-BWIDTH64+1; i+=BWIDTH64){
				ByteVector va=ByteVector.fromArray(BSPECIES64, a, i);
				c+=((IntVector)va.convertShape(VectorOperators.B2I, ISPECIES128, 0)).reduceLanesToLong(VectorOperators.ADD);
				c+=((IntVector)va.convertShape(VectorOperators.B2I, ISPECIES128, 1)).reduceLanesToLong(VectorOperators.ADD);
			}
		}

		//Scalar tail
		for(; i<=to; i++){c+=a[i];}
		return c;
	}

	/**
	 * Find positions of the given symbol in a byte array using SIMD.
	 * @param buffer The byte array to search
	 * @param from Starting position (inclusive)
	 * @param to Ending position (exclusive)
	 * @param symbol Character to find
	 * @param positions IntList to store newline positions
	 * @return Number of symbols found, including pre-existing ones
	 */
	static final int findSymbols(final byte[] buffer, final int from,
		final int to, final byte symbol, final IntList positions){
		int i=from;

		{//128-bit loop
			final ByteVector newlineVec128=ByteVector.broadcast(BSPECIES128, symbol);
			for(; i<=to-BWIDTH128; i+=BWIDTH128){
				ByteVector vec=ByteVector.fromArray(BSPECIES128, buffer, i);
				VectorMask<Byte> mask=vec.eq(newlineVec128);
				long bits=mask.toLong();

				while(bits!=0){
					int lane=Long.numberOfTrailingZeros(bits);
					positions.add(i+lane);
					bits&=(bits-1); //Clear lowest set bit
				}
			}
		}

		{//64-bit loop
			final ByteVector newlineVec64=ByteVector.broadcast(BSPECIES64, symbol);
			for(; i<=to-BWIDTH64; i+=BWIDTH64){
				ByteVector vec=ByteVector.fromArray(BSPECIES64, buffer, i);
				VectorMask<Byte> mask=vec.eq(newlineVec64);
				long bits=mask.toLong();

				while(bits!=0){
					int lane=Long.numberOfTrailingZeros(bits);
					positions.add(i+lane);
					bits&=(bits-1);
				}
			}
		}

		//Scalar tail
		for(; i<to; i++){
			if(buffer[i]==symbol){positions.add(i);}
		}

		return positions.size();
	}

	/**
	 * Count positions of the given symbol in a byte array using SIMD.
	 * @param buffer The byte array to search
	 * @param from Starting position (inclusive)
	 * @param to Ending position (exclusive)
	 * @param symbol Character to find
	 * @return Number of symbols found, including pre-existing ones
	 */
	static final int countSymbols(final byte[] buffer, final int from,
		final int to, final byte symbol){
		int positions=0;
		int i=from;

		{//128-bit loop
			final ByteVector newlineVec128=ByteVector.broadcast(BSPECIES128, symbol);
			for(; i<=to-BWIDTH128; i+=BWIDTH128){
				ByteVector vec=ByteVector.fromArray(BSPECIES128, buffer, i);
				VectorMask<Byte> mask=vec.eq(newlineVec128);
				long bits=mask.toLong();
				positions+=Long.bitCount(bits);
			}
		}

		{//64-bit loop
			final ByteVector newlineVec64=ByteVector.broadcast(BSPECIES64, symbol);
			for(; i<=to-BWIDTH64; i+=BWIDTH64){
				ByteVector vec=ByteVector.fromArray(BSPECIES64, buffer, i);
				VectorMask<Byte> mask=vec.eq(newlineVec64);
				long bits=mask.toLong();
				positions+=Long.bitCount(bits);
			}
		}

		//Scalar tail
		for(; i<to; i++){
			if(buffer[i]==symbol){positions++;}
		}

		return positions;
	}

	/**
	 * Find the last symbol in buffer by scanning backwards using SIMD.
	 * @param buffer Buffer to scan
	 * @param limit Scan backwards from this position (exclusive)
	 * @param symbol Symbol to find
	 * @return Position of last symbol, or -1 if none found
	 */
	static int findLastSymbol(byte[] buffer, int limit, final byte symbol){
		int i=limit;

		{//128-bit loop - work backwards
			final ByteVector newlineVec128=ByteVector.broadcast(BSPECIES128, symbol);
			for(; i>=BWIDTH128; i-=BWIDTH128){
				ByteVector vec=ByteVector.fromArray(BSPECIES128, buffer, i-BWIDTH128);
				VectorMask<Byte> mask=vec.eq(newlineVec128);
				long bits=mask.toLong();

				if(bits!=0){
					int lane=63-Long.numberOfLeadingZeros(bits);
					return (i-BWIDTH128)+lane;
				}
			}
		}

		{//64-bit loop
			final ByteVector newlineVec64=ByteVector.broadcast(BSPECIES64, symbol);
			for(; i>=BWIDTH64; i-=BWIDTH64){
				ByteVector vec=ByteVector.fromArray(BSPECIES64, buffer, i-BWIDTH64);
				VectorMask<Byte> mask=vec.eq(newlineVec64);
				long bits=mask.toLong();

				if(bits!=0){
					int lane=63-Long.numberOfLeadingZeros(bits);
					return (i-BWIDTH64)+lane;
				}
			}
		}

		//Scalar tail
		for(i--; i>=0; i--){
			if(buffer[i]==symbol){return i;}
		}

		return -1;
	}

	/**
	 * Find FASTA record boundaries (\n>) using SIMD.
	 * @param buffer The byte array to search
	 * @param from Starting position (inclusive)
	 * @param to Ending position (exclusive)
	 * @param positions IntList to store \n positions before >
	 * @return Number of boundaries found
	 */
	static final int findFastaHeaders(final byte[] buffer, final int from,
		final int to, final IntList positions){
		int i=from+1; // Need to check i-1, so start at from+1

		{//128-bit loop
			final ByteVector carrotVec128=ByteVector.broadcast(BSPECIES128, (byte)'>');
			final ByteVector slashnVec128=ByteVector.broadcast(BSPECIES128, (byte)'\n');

			for(; i<=to-BWIDTH128; i+=BWIDTH128){
				ByteVector vec=ByteVector.fromArray(BSPECIES128, buffer, i);
				ByteVector vecPrev=ByteVector.fromArray(BSPECIES128, buffer, i-1);

				VectorMask<Byte> isCarrot=vec.eq(carrotVec128);
				VectorMask<Byte> prevIsSlashn=vecPrev.eq(slashnVec128);
				VectorMask<Byte> isBoundary=isCarrot.and(prevIsSlashn);

				long bits=isBoundary.toLong();
				while(bits!=0){
					int lane=Long.numberOfTrailingZeros(bits);
					positions.add(i+lane-1); // Store \n position
					bits&=(bits-1);
				}
			}
		}

		{//64-bit loop
			final ByteVector carrotVec64=ByteVector.broadcast(BSPECIES64, (byte)'>');
			final ByteVector slashnVec64=ByteVector.broadcast(BSPECIES64, (byte)'\n');

			for(; i<=to-BWIDTH64; i+=BWIDTH64){
				ByteVector vec=ByteVector.fromArray(BSPECIES64, buffer, i);
				ByteVector vecPrev=ByteVector.fromArray(BSPECIES64, buffer, i-1);

				VectorMask<Byte> isCarrot=vec.eq(carrotVec64);
				VectorMask<Byte> prevIsSlashn=vecPrev.eq(slashnVec64);
				VectorMask<Byte> isBoundary=isCarrot.and(prevIsSlashn);

				long bits=isBoundary.toLong();
				while(bits!=0){
					int lane=Long.numberOfTrailingZeros(bits);
					positions.add(i+lane-1);
					bits&=(bits-1);
				}
			}
		}

		//Scalar tail
		for(; i<to; i++){
			if(buffer[i]=='>' && buffer[i-1]=='\n'){
				positions.add(i-1);
			}
		}

		return positions.size();
	}

	/**
	 * Adds a constant delta to all bytes in an array using SIMD.
	 * @param array The byte array to modify in-place.
	 * @param delta The value to add to each element.
	 */
	static final void add(final byte[] array, final byte delta){
		if(array==null){return;}

		final int length=array.length;
		int i=0;

		{//128-bit loop
			final ByteVector vdelta128=ByteVector.broadcast(BSPECIES128, delta);
			for(; i<=length-BWIDTH128; i+=BWIDTH128){
				ByteVector va=ByteVector.fromArray(BSPECIES128, array, i);
				ByteVector vresult=va.add(vdelta128);
				vresult.intoArray(array, i);
			}
		}

		{//64-bit loop
			final ByteVector vdelta64=ByteVector.broadcast(BSPECIES64, delta);
			for(; i<=length-BWIDTH64; i+=BWIDTH64){
				ByteVector va=ByteVector.fromArray(BSPECIES64, array, i);
				ByteVector vresult=va.add(vdelta64);
				vresult.intoArray(array, i);
			}
		}

		//Scalar tail
		for(; i<length; i++){
			array[i]+=delta;
		}
	}

	/**
	 * Adds in int precision, applies a lower bound, then stores byte casts.
	 * @param array The byte array to modify in-place.
	 * @param delta The value to add to each element.
	 * @param cap Lower bound before narrowing each result to byte.
	 * @return Byte cast of min(127, smallest pre-cap sum); 0 for null, 127 for empty.
	 */
	static final byte addAndCapMin(final byte[] array, final byte delta, final int cap){
		if(array==null){return 0;}

		final int length=array.length;
		int i=0;
		int min=127;
		final long threshold=(long)cap-delta;
		final byte thresholdByte=(byte)Math.max(Byte.MIN_VALUE, Math.min(Byte.MAX_VALUE, threshold));
		final byte outputDelta=(threshold>Byte.MAX_VALUE ? (byte)(cap-Byte.MAX_VALUE) : delta);

		// max(x, cap-delta)+delta equals max(x+delta, cap) when the threshold is a byte.
		// Thresholds at or below byte range make max an identity; thresholds above range
		// clamp every lane to 127, then adjust outputDelta to produce (byte)cap.
		{//128-bit loop
			ByteVector vdelta128=ByteVector.broadcast(BSPECIES128, outputDelta);
			ByteVector vthreshold128=ByteVector.broadcast(BSPECIES128, thresholdByte);
			ByteVector vmin=ByteVector.broadcast(BSPECIES128, (byte)127);
			final int start=i;

			for(; i<=length-BWIDTH128; i+=BWIDTH128){
				ByteVector va=ByteVector.fromArray(BSPECIES128, array, i);
				ByteVector vresult=va.max(vthreshold128).add(vdelta128);
				vmin=vmin.min(va);
				vresult.intoArray(array, i);
			}

			if(i>start){min=Math.min(min, vmin.reduceLanes(VectorOperators.MIN)+delta);}
		}

		{//64-bit loop
			ByteVector vdelta64=ByteVector.broadcast(BSPECIES64, outputDelta);
			ByteVector vthreshold64=ByteVector.broadcast(BSPECIES64, thresholdByte);
			ByteVector vmin=ByteVector.broadcast(BSPECIES64, (byte)127);
			final int start=i;

			for(; i<=length-BWIDTH64; i+=BWIDTH64){
				ByteVector va=ByteVector.fromArray(BSPECIES64, array, i);
				ByteVector vresult=va.max(vthreshold64).add(vdelta64);
				vmin=vmin.min(va);
				vresult.intoArray(array, i);
			}

			if(i>start){min=Math.min(min, vmin.reduceLanes(VectorOperators.MIN)+delta);}
		}

		//Scalar tail
		for(; i<length; i++){
			int b=array[i]+delta;
			min=Math.min(min, b);
			array[i]=(byte)Math.max(cap, b);
		}

		return (byte)min;
	}

	/**
	 * Applies quality offset delta, zeros quality for N bases, caps others at 2.
	 * @param quals Quality array to modify in-place.
	 * @param bases Sequence bases array.
	 * @param delta The offset to add to quality scores.
	 */
	static final void applyQualOffset(final byte[] quals, final byte[] bases, final int delta){
		if(quals==null){return;}

		final int length=quals.length;
		int i=0;

		{//128-bit loop
			ByteVector vdelta128=ByteVector.broadcast(BSPECIES128, (byte)delta);
			ByteVector vn128=ByteVector.broadcast(BSPECIES128, (byte)'N');
			ByteVector vzero128=ByteVector.broadcast(BSPECIES128, (byte)0);
			ByteVector vcap2_128=ByteVector.broadcast(BSPECIES128, (byte)2);

			for(; i<=length-BWIDTH128; i+=BWIDTH128){
				ByteVector vquals=ByteVector.fromArray(BSPECIES128, quals, i);
				ByteVector vbases=ByteVector.fromArray(BSPECIES128, bases, i);
				ByteVector vresult=vquals.add(vdelta128);
				vresult=vresult.max(vcap2_128);
				//KNOWN/WONTFIX [simd/SIMDByte256#002], inherited: this SIMD path zeros quality only for uppercase 'N', while the scalar twin
				//(Vector.applyQualOffset) zeros it for EVERY undefined base (baseToNumber[b]<0: lowercase 'n', IUPAC, '.', '-'). So
				//SIMD-on vs SIMD-off give DIFFERENT quality for a base+qscore pair whose base is lowercase/IUPAC - but that never
				//legitimately occurs (reads are uppercase ACGTN; contigs carry no qscores). WONTFIX (Brian): matching scalar would
				//add ~10 vector ops per 32-byte block to the original 256-bit every-read quality path, not worth it for input that cannot arise.
				VectorMask<Byte> maskN=vbases.eq(vn128);
				vresult=vresult.blend(vzero128, maskN);
				vresult.intoArray(quals, i);
			}
		}

		{//64-bit loop
			ByteVector vdelta64=ByteVector.broadcast(BSPECIES64, (byte)delta);
			ByteVector vn64=ByteVector.broadcast(BSPECIES64, (byte)'N');
			ByteVector vzero64=ByteVector.broadcast(BSPECIES64, (byte)0);
			ByteVector vcap2_64=ByteVector.broadcast(BSPECIES64, (byte)2);

			for(; i<=length-BWIDTH64; i+=BWIDTH64){
				ByteVector vquals=ByteVector.fromArray(BSPECIES64, quals, i);
				ByteVector vbases=ByteVector.fromArray(BSPECIES64, bases, i);
				ByteVector vresult=vquals.add(vdelta64);
				vresult=vresult.max(vcap2_64);
				VectorMask<Byte> maskN=vbases.eq(vn64);
				vresult=vresult.blend(vzero64, maskN);
				vresult.intoArray(quals, i);
			}
		}

		//Scalar tail
		for(; i<length; i++){
			byte b=bases[i];
			int q=quals[i]+delta;
			q=(AminoAcid.baseToNumber[b]<0 ? 0 : Math.max(2, q));
			quals[i]=(byte)q;
		}
	}

	/**
	 * Zeros quality for N bases, caps others between 2 and max.
	 * @param quals Quality array to modify in-place.
	 * @param bases Sequence bases array.
	 * @param max Maximum quality
	 */
	static final void capQuality(final byte[] quals, final byte[] bases, final byte max){
		if(quals==null){return;}

		final int length=quals.length;
		int i=0;

		{//128-bit loop
			ByteVector vn128=ByteVector.broadcast(BSPECIES128, (byte)'N');
			ByteVector vzero128=ByteVector.broadcast(BSPECIES128, (byte)0);
			ByteVector vcap2_128=ByteVector.broadcast(BSPECIES128, (byte)2);

			for(; i<=length-BWIDTH128; i+=BWIDTH128){
				ByteVector vquals=ByteVector.fromArray(BSPECIES128, quals, i);
				ByteVector vbases=ByteVector.fromArray(BSPECIES128, bases, i);
				ByteVector vresult=vquals.max(vcap2_128).min(ByteVector.broadcast(BSPECIES128, max));
				//KNOWN/WONTFIX [simd/SIMDByte256#002], inherited: same 'N'-only vs baseToNumber[b]<0 divergence as applyQualOffset (twin). Scalar
				//Vector.capQuality zeros quality for ALL undefined bases; this SIMD path only for uppercase 'N'. Differs only on a
				//lowercase/IUPAC base carrying a qscore, which never legitimately occurs (uppercase-ACGTN reads; contigs have no quals).
				//WONTFIX (Brian): ~10 extra vector ops per 32-byte block on the original 256-bit hot quality path was not worth handling
				//input that cannot arise.
				VectorMask<Byte> maskN=vbases.eq(vn128);
				vresult=vresult.blend(vzero128, maskN);
				vresult.intoArray(quals, i);
			}
		}

		{//64-bit loop
			ByteVector vn64=ByteVector.broadcast(BSPECIES64, (byte)'N');
			ByteVector vzero64=ByteVector.broadcast(BSPECIES64, (byte)0);
			ByteVector vcap2_64=ByteVector.broadcast(BSPECIES64, (byte)2);

			for(; i<=length-BWIDTH64; i+=BWIDTH64){
				ByteVector vquals=ByteVector.fromArray(BSPECIES64, quals, i);
				ByteVector vbases=ByteVector.fromArray(BSPECIES64, bases, i);
				ByteVector vresult=vquals.max(vcap2_64).min(ByteVector.broadcast(BSPECIES64, max));
				VectorMask<Byte> maskN=vbases.eq(vn64);
				vresult=vresult.blend(vzero64, maskN);
				vresult.intoArray(quals, i);
			}
		}

		//Scalar tail
		for(; i<length; i++){
			byte b=bases[i];
			int q=quals[i];
			q=(AminoAcid.baseToNumber[b]<0 ? 0 : Math.min(Math.max(2, q), max));
			quals[i]=(byte)q;
		}
	}

	/**
	 * Converts U to T and u to t in a byte array using SIMD.
	 * @param bases The base array to modify in-place.
	 */
	static final void uToT(final byte[] bases){
		if(bases==null){return;}

		final int length=bases.length;
		int i=0;

		{//128-bit loop
			ByteVector vU128=ByteVector.broadcast(BSPECIES128, (byte)'U');
			ByteVector vu128=ByteVector.broadcast(BSPECIES128, (byte)'u');
			ByteVector vT128=ByteVector.broadcast(BSPECIES128, (byte)'T');
			ByteVector vt128=ByteVector.broadcast(BSPECIES128, (byte)'t');

			for(; i<=length-BWIDTH128; i+=BWIDTH128){
				ByteVector vbases=ByteVector.fromArray(BSPECIES128, bases, i);
				VectorMask<Byte> maskU=vbases.eq(vU128);
				VectorMask<Byte> masku=vbases.eq(vu128);
				vbases=vbases.blend(vT128, maskU);
				vbases=vbases.blend(vt128, masku);
				vbases.intoArray(bases, i);
			}
		}

		{//64-bit loop
			ByteVector vU64=ByteVector.broadcast(BSPECIES64, (byte)'U');
			ByteVector vu64=ByteVector.broadcast(BSPECIES64, (byte)'u');
			ByteVector vT64=ByteVector.broadcast(BSPECIES64, (byte)'T');
			ByteVector vt64=ByteVector.broadcast(BSPECIES64, (byte)'t');

			for(; i<=length-BWIDTH64; i+=BWIDTH64){
				ByteVector vbases=ByteVector.fromArray(BSPECIES64, bases, i);
				VectorMask<Byte> maskU=vbases.eq(vU64);
				VectorMask<Byte> masku=vbases.eq(vu64);
				vbases=vbases.blend(vT64, maskU);
				vbases=vbases.blend(vt64, masku);
				vbases.intoArray(bases, i);
			}
		}

		//Scalar tail
		for(; i<length; i++){
			bases[i]=AminoAcid.uToT[bases[i]];
		}
	}

	/** Converts lowercase ASCII letters to N in-place.
	 * @param array byte array to modify; null is accepted
	 * @return true */
	static final boolean lowerCaseToN(final byte[] array){
		if(array==null){return true;}

		final int length=array.length;
		int i=0;
		final byte a='a', z='z', N='N';

		{//128-bit loop
			ByteVector va128=ByteVector.broadcast(BSPECIES128, a);
			ByteVector vz128=ByteVector.broadcast(BSPECIES128, z);
			ByteVector vN128=ByteVector.broadcast(BSPECIES128, N);

			for(; i<=length-BWIDTH128; i+=BWIDTH128){
				ByteVector vb=ByteVector.fromArray(BSPECIES128, array, i);
				VectorMask<Byte> maskLower=vb.compare(VectorOperators.GE, va128)
					.and(vb.compare(VectorOperators.LE, vz128));
				ByteVector vresult=vb.blend(vN128, maskLower);
				vresult.intoArray(array, i);
			}
		}

		{//64-bit loop
			ByteVector va64=ByteVector.broadcast(BSPECIES64, a);
			ByteVector vz64=ByteVector.broadcast(BSPECIES64, z);
			ByteVector vN64=ByteVector.broadcast(BSPECIES64, N);

			for(; i<=length-BWIDTH64; i+=BWIDTH64){
				ByteVector vb=ByteVector.fromArray(BSPECIES64, array, i);
				VectorMask<Byte> maskLower=vb.compare(VectorOperators.GE, va64)
					.and(vb.compare(VectorOperators.LE, vz64));
				ByteVector vresult=vb.blend(vN64, maskLower);
				vresult.intoArray(array, i);
			}
		}

		//Scalar tail
		for(; i<length; i++){
			array[i]=AminoAcid.lowerCaseToNocall[array[i]];
		}

		return true;
	}

	/** Converts dot, dash, X/x, and lowercase n to N in-place.
	 * @param array byte array to modify; null is accepted
	 * @return true */
	static final boolean dotDashXToN(final byte[] array){
		if(array==null){return true;}

		final int length=array.length;
		int i=0;
		final byte dot='.', dash='-', X='X', x='x', n='n', N='N';

		{//128-bit loop
			ByteVector vdot128=ByteVector.broadcast(BSPECIES128, dot);
			ByteVector vdash128=ByteVector.broadcast(BSPECIES128, dash);
			ByteVector vX128=ByteVector.broadcast(BSPECIES128, X);
			ByteVector vx128=ByteVector.broadcast(BSPECIES128, x);
			ByteVector vn128=ByteVector.broadcast(BSPECIES128, n);
			ByteVector vN128=ByteVector.broadcast(BSPECIES128, N);

			for(; i<=length-BWIDTH128; i+=BWIDTH128){
				ByteVector vb=ByteVector.fromArray(BSPECIES128, array, i);
				VectorMask<Byte> maskDot=vb.eq(vdot128);
				VectorMask<Byte> maskDash=vb.eq(vdash128);
				VectorMask<Byte> maskX=vb.eq(vX128);
				VectorMask<Byte> maskx=vb.eq(vx128);
				VectorMask<Byte> maskn=vb.eq(vn128);
				VectorMask<Byte> maskAny=maskDot.or(maskDash).or(maskX).or(maskx).or(maskn);
				ByteVector vresult=vb.blend(vN128, maskAny);
				vresult.intoArray(array, i);
			}
		}

		{//64-bit loop
			ByteVector vdot64=ByteVector.broadcast(BSPECIES64, dot);
			ByteVector vdash64=ByteVector.broadcast(BSPECIES64, dash);
			ByteVector vX64=ByteVector.broadcast(BSPECIES64, X);
			ByteVector vx64=ByteVector.broadcast(BSPECIES64, x);
			ByteVector vn64=ByteVector.broadcast(BSPECIES64, n);
			ByteVector vN64=ByteVector.broadcast(BSPECIES64, N);

			for(; i<=length-BWIDTH64; i+=BWIDTH64){
				ByteVector vb=ByteVector.fromArray(BSPECIES64, array, i);
				VectorMask<Byte> maskDot=vb.eq(vdot64);
				VectorMask<Byte> maskDash=vb.eq(vdash64);
				VectorMask<Byte> maskX=vb.eq(vX64);
				VectorMask<Byte> maskx=vb.eq(vx64);
				VectorMask<Byte> maskn=vb.eq(vn64);
				VectorMask<Byte> maskAny=maskDot.or(maskDash).or(maskX).or(maskx).or(maskn);
				ByteVector vresult=vb.blend(vN64, maskAny);
				vresult.intoArray(array, i);
			}
		}

		//Scalar tail
		for(; i<length; i++){
			array[i]=AminoAcid.dotDashXToNocall[array[i]];
		}

		return true;
	}

	/** Checks for amino-acid symbols that are not IUPAC nucleotide bases.
	 * Detects E/F/I/J/L/O/P/Q/Z case-insensitively, plus '*', '.', and '-'.
	 * @param array byte array to check; null is accepted
	 * @return true if any inspected byte is amino-only */
	static final boolean isProtein(final byte[] array){
		if(array==null){return false;}

		final int length=array.length;
		int i=0;
		final byte mask=(byte)~32;
		final byte E='E', F='F', I='I', J='J', L='L', O='O', P='P', Q='Q', Z='Z';
		final byte star='*', dot='.', dash='-';

		boolean protein=false;

		{//128-bit loop
			VectorMask<Byte> foundProtein=BSPECIES128.maskAll(false);
			ByteVector vmask128=ByteVector.broadcast(BSPECIES128, mask);
			ByteVector vE128=ByteVector.broadcast(BSPECIES128, E);
			ByteVector vF128=ByteVector.broadcast(BSPECIES128, F);
			ByteVector vI128=ByteVector.broadcast(BSPECIES128, I);
			ByteVector vJ128=ByteVector.broadcast(BSPECIES128, J);
			ByteVector vL128=ByteVector.broadcast(BSPECIES128, L);
			ByteVector vO128=ByteVector.broadcast(BSPECIES128, O);
			ByteVector vP128=ByteVector.broadcast(BSPECIES128, P);
			ByteVector vQ128=ByteVector.broadcast(BSPECIES128, Q);
			ByteVector vZ128=ByteVector.broadcast(BSPECIES128, Z);
			ByteVector vstar128=ByteVector.broadcast(BSPECIES128, star);
			ByteVector vdot128=ByteVector.broadcast(BSPECIES128, dot);
			ByteVector vdash128=ByteVector.broadcast(BSPECIES128, dash);

			for(; i<=length-BWIDTH128; i+=BWIDTH128){
				ByteVector vb=ByteVector.fromArray(BSPECIES128, array, i);
				ByteVector vbUpper=vb.and(vmask128);
				VectorMask<Byte> isProtein=vbUpper.eq(vE128).or(vbUpper.eq(vF128))
						.or(vbUpper.eq(vI128)).or(vbUpper.eq(vJ128)).or(vbUpper.eq(vL128))
						.or(vbUpper.eq(vO128)).or(vbUpper.eq(vP128)).or(vbUpper.eq(vQ128))
						.or(vbUpper.eq(vZ128)).or(vb.eq(vstar128)).or(vb.eq(vdot128))
						.or(vb.eq(vdash128));
				foundProtein=isProtein.or(foundProtein);
			}

			protein=foundProtein.anyTrue();
		}

		{//64-bit loop
			VectorMask<Byte> foundProtein=BSPECIES64.maskAll(false);
			ByteVector vmask64=ByteVector.broadcast(BSPECIES64, mask);
			ByteVector vE64=ByteVector.broadcast(BSPECIES64, E);
			ByteVector vF64=ByteVector.broadcast(BSPECIES64, F);
			ByteVector vI64=ByteVector.broadcast(BSPECIES64, I);
			ByteVector vJ64=ByteVector.broadcast(BSPECIES64, J);
			ByteVector vL64=ByteVector.broadcast(BSPECIES64, L);
			ByteVector vO64=ByteVector.broadcast(BSPECIES64, O);
			ByteVector vP64=ByteVector.broadcast(BSPECIES64, P);
			ByteVector vQ64=ByteVector.broadcast(BSPECIES64, Q);
			ByteVector vZ64=ByteVector.broadcast(BSPECIES64, Z);
			ByteVector vstar64=ByteVector.broadcast(BSPECIES64, star);
			ByteVector vdot64=ByteVector.broadcast(BSPECIES64, dot);
			ByteVector vdash64=ByteVector.broadcast(BSPECIES64, dash);

			for(; i<=length-BWIDTH64; i+=BWIDTH64){
				ByteVector vb=ByteVector.fromArray(BSPECIES64, array, i);
				ByteVector vbUpper=vb.and(vmask64);
				VectorMask<Byte> isProtein=vbUpper.eq(vE64).or(vbUpper.eq(vF64))
						.or(vbUpper.eq(vI64)).or(vbUpper.eq(vJ64)).or(vbUpper.eq(vL64))
						.or(vbUpper.eq(vO64)).or(vbUpper.eq(vP64)).or(vbUpper.eq(vQ64))
						.or(vbUpper.eq(vZ64)).or(vb.eq(vstar64)).or(vb.eq(vdot64))
						.or(vb.eq(vdash64));
				foundProtein=isProtein.or(foundProtein);
			}

			protein|=foundProtein.anyTrue();
		}

		//Scalar tail
		for(; i<length; i++){
			byte b=array[i];
			boolean nuc=AminoAcid.baseToNumberExtended[b]>=0;
			boolean amino=AminoAcid.acidToNumberExtended[b]>=0;
			protein|=(amino && !nuc);
		}

		return protein;
	}

	/** Converts lowercase ASCII letters to uppercase in-place.
	 * @param array byte array to modify; null is accepted
	 * @return true */
	static final boolean toUpperCase(final byte[] array){
		if(array==null){return true;}

		final int length=array.length;
		int i=0;
		final byte A='A', Z='Z';
		final byte mask=(byte) ~32; // Cast needed as ~32 is int

		{//128-bit loop
			ByteVector vmask128=ByteVector.broadcast(BSPECIES128, mask);
			ByteVector vA128=ByteVector.broadcast(BSPECIES128, A);
			ByteVector vZ128=ByteVector.broadcast(BSPECIES128, Z);

			for(; i<=length-BWIDTH128; i+=BWIDTH128){
				ByteVector vb0=ByteVector.fromArray(BSPECIES128, array, i);
				ByteVector vbMasked=vb0.and(vmask128); // Force potential uppercase

				// Check if the masked version is a valid letter (A-Z)
				// This correctly identifies both 'a'->'A' and 'A'->'A' as valid.
				VectorMask<Byte> isLetter=vbMasked.compare(VectorOperators.GE, vA128)
					.and(vbMasked.compare(VectorOperators.LE, vZ128));

				// BLEND: If isLetter (Mask True) -> Take vbMasked.
				//        If !isLetter (Mask False) -> Keep vb0 (original).
				vb0.blend(vbMasked, isLetter).intoArray(array, i);
			}
		}

		{//64-bit loop
			ByteVector vmask64=ByteVector.broadcast(BSPECIES64, mask);
			ByteVector vA64=ByteVector.broadcast(BSPECIES64, A);
			ByteVector vZ64=ByteVector.broadcast(BSPECIES64, Z);

			for(; i<=length-BWIDTH64; i+=BWIDTH64){
				ByteVector vb0=ByteVector.fromArray(BSPECIES64, array, i);
				ByteVector vbMasked=vb0.and(vmask64);

				VectorMask<Byte> isLetter=vbMasked.compare(VectorOperators.GE, vA64)
					.and(vbMasked.compare(VectorOperators.LE, vZ64));

				vb0.blend(vbMasked, isLetter).intoArray(array, i);
			}
		}

		//Scalar tail
		for(; i<length; i++){
			array[i]=AminoAcid.toUpperCase[array[i]];
		}

		return true;
	}

	/** Checks if all bytes are ASCII letters.
	 * @param array byte array to check; null is accepted
	 * @return true if every byte is A-Z or a-z */
	static final boolean allLetters(final byte[] array){
		if(array==null){return true;}

		final int length=array.length;
		int i=0;
		final byte A='A', Z='Z';
		final byte mask=~32;

		boolean success=true;

		{//128-bit loop
			VectorMask<Byte> invalid=BSPECIES128.maskAll(false);
			ByteVector vmask128=ByteVector.broadcast(BSPECIES128, mask);
			ByteVector vA128=ByteVector.broadcast(BSPECIES128, A);
			ByteVector vZ128=ByteVector.broadcast(BSPECIES128, Z);

			for(; i<=length-BWIDTH128; i+=BWIDTH128){
				ByteVector vb0=ByteVector.fromArray(BSPECIES128, array, i);
				ByteVector vb=vb0.and(vmask128);
				VectorMask<Byte> validLow=vb.compare(VectorOperators.GE, vA128);
				VectorMask<Byte> validHigh=vb.compare(VectorOperators.LE, vZ128);
				VectorMask<Byte> valid=validLow.and(validHigh);
				invalid=valid.not().or(invalid);
			}

			success=!invalid.anyTrue();
		}

		{//64-bit loop
			VectorMask<Byte> invalid=BSPECIES64.maskAll(false);
			ByteVector vmask64=ByteVector.broadcast(BSPECIES64, mask);
			ByteVector vA64=ByteVector.broadcast(BSPECIES64, A);
			ByteVector vZ64=ByteVector.broadcast(BSPECIES64, Z);

			for(; i<=length-BWIDTH64; i+=BWIDTH64){
				ByteVector vb0=ByteVector.fromArray(BSPECIES64, array, i);
				ByteVector vb=vb0.and(vmask64);
				VectorMask<Byte> validLow=vb.compare(VectorOperators.GE, vA64);
				VectorMask<Byte> validHigh=vb.compare(VectorOperators.LE, vZ64);
				VectorMask<Byte> valid=validLow.and(validHigh);
				invalid=valid.not().or(invalid);
			}

			success&=!invalid.anyTrue();
		}

		//Scalar tail
		for(; i<length; i++){
			int b=(array[i]&mask);
			success&=(b>=A && b<=Z);
		}

		return success;
	}

	/** Masks IUPAC ambiguity letters in-place.
	 * B/D/H/K/M/N/R/S/V/W/X/Y in either case become N; A/C/G/T/U
	 * and unrelated bytes remain unchanged.
	 * @param array byte array to modify; null is accepted
	 * @return true */
	static final boolean iupacToN(final byte[] array){
		if(array==null){return true;}

		final int length=array.length;
		int i=0;
		final byte B='B', D='D', H='H', K='K', M='M', N='N';
		final byte R='R', S='S', V='V', W='W', X='X', Y='Y';
		final byte mask=~32;

		{//128-bit loop
			ByteVector vmask128=ByteVector.broadcast(BSPECIES128, mask);
			ByteVector vB128=ByteVector.broadcast(BSPECIES128, B);
			ByteVector vD128=ByteVector.broadcast(BSPECIES128, D);
			ByteVector vH128=ByteVector.broadcast(BSPECIES128, H);
			ByteVector vK128=ByteVector.broadcast(BSPECIES128, K);
			ByteVector vM128=ByteVector.broadcast(BSPECIES128, M);
			ByteVector vN128=ByteVector.broadcast(BSPECIES128, N);
			ByteVector vR128=ByteVector.broadcast(BSPECIES128, R);
			ByteVector vS128=ByteVector.broadcast(BSPECIES128, S);
			ByteVector vV128=ByteVector.broadcast(BSPECIES128, V);
			ByteVector vW128=ByteVector.broadcast(BSPECIES128, W);
			ByteVector vX128=ByteVector.broadcast(BSPECIES128, X);
			ByteVector vY128=ByteVector.broadcast(BSPECIES128, Y);

			for(; i<=length-BWIDTH128; i+=BWIDTH128){
				ByteVector vb0=ByteVector.fromArray(BSPECIES128, array, i);
				ByteVector vb=vb0.and(vmask128);
				VectorMask<Byte> isIupac=vb.eq(vB128).or(vb.eq(vD128))
						.or(vb.eq(vH128)).or(vb.eq(vK128)).or(vb.eq(vM128))
						.or(vb.eq(vN128)).or(vb.eq(vR128)).or(vb.eq(vS128))
						.or(vb.eq(vV128)).or(vb.eq(vW128)).or(vb.eq(vX128))
						.or(vb.eq(vY128));
				vb0.blend(vN128, isIupac).intoArray(array, i);
			}
		}

		{//64-bit loop
			ByteVector vmask64=ByteVector.broadcast(BSPECIES64, mask);
			ByteVector vB64=ByteVector.broadcast(BSPECIES64, B);
			ByteVector vD64=ByteVector.broadcast(BSPECIES64, D);
			ByteVector vH64=ByteVector.broadcast(BSPECIES64, H);
			ByteVector vK64=ByteVector.broadcast(BSPECIES64, K);
			ByteVector vM64=ByteVector.broadcast(BSPECIES64, M);
			ByteVector vN64=ByteVector.broadcast(BSPECIES64, N);
			ByteVector vR64=ByteVector.broadcast(BSPECIES64, R);
			ByteVector vS64=ByteVector.broadcast(BSPECIES64, S);
			ByteVector vV64=ByteVector.broadcast(BSPECIES64, V);
			ByteVector vW64=ByteVector.broadcast(BSPECIES64, W);
			ByteVector vX64=ByteVector.broadcast(BSPECIES64, X);
			ByteVector vY64=ByteVector.broadcast(BSPECIES64, Y);

			for(; i<=length-BWIDTH64; i+=BWIDTH64){
				ByteVector vb0=ByteVector.fromArray(BSPECIES64, array, i);
				ByteVector vb=vb0.and(vmask64);
				VectorMask<Byte> isIupac=vb.eq(vB64).or(vb.eq(vD64))
						.or(vb.eq(vH64)).or(vb.eq(vK64)).or(vb.eq(vM64))
						.or(vb.eq(vN64)).or(vb.eq(vR64)).or(vb.eq(vS64))
						.or(vb.eq(vV64)).or(vb.eq(vW64)).or(vb.eq(vX64))
						.or(vb.eq(vY64));
				vb0.blend(vN64, isIupac).intoArray(array, i);
			}
		}

		//Scalar tail
		for(; i<length; i++){
			array[i]=AminoAcid.iupacToNocall[array[i]];
		}

		return true;
	}

	/** Counts query/reference mismatches from a reference offset.
	 * A query N always counts as a mismatch, including when the reference is also N.
	 * Stops after the vector or scalar block that brings {@code subs} above
	 * {@code maxSubs}, so the returned count may exceed {@code maxSubs}.
	 * @param query query bases
	 * @param ref reference bases
	 * @param maxSubs mismatch budget for early termination
	 * @param rStart starting offset in {@code ref}
	 * @return counted mismatches */
	static int alignSIMD(byte[] query, byte[] ref, final int maxSubs, final int rStart){
		int subs=0;
		int i=0;
		int j=rStart;

		{//128-bit loop
			ByteVector vN128=ByteVector.broadcast(BSPECIES128, 'N');
			for(final int lim=query.length-BWIDTH128; i<=lim && subs<=maxSubs; i+=BWIDTH128, j+=BWIDTH128){
				ByteVector qVec=ByteVector.fromArray(BSPECIES128, query, i);
				ByteVector rVec=ByteVector.fromArray(BSPECIES128, ref, j);
				VectorMask<Byte> subMask=qVec.eq(rVec).not().or(qVec.eq(vN128));
				long subBits=subMask.toLong();
				subs+=Long.bitCount(subBits);
			}
		}

		{//64-bit loop
			ByteVector vN64=ByteVector.broadcast(BSPECIES64, 'N');
			for(final int lim=query.length-BWIDTH64; i<=lim && subs<=maxSubs; i+=BWIDTH64, j+=BWIDTH64){
				ByteVector qVec=ByteVector.fromArray(BSPECIES64, query, i);
				ByteVector rVec=ByteVector.fromArray(BSPECIES64, ref, j);
				VectorMask<Byte> subMask=qVec.eq(rVec).not().or(qVec.eq(vN64));
				long subBits=subMask.toLong();
				subs+=Long.bitCount(subBits);
			}
		}

		//Scalar tail
		for(; i<query.length && subs<=maxSubs; i++, j++){
			subs+=(query[i]!=ref[j] || query[i]=='N') ? 1 : 0;
		}

		return subs;
	}

	/** Counts mismatches against a fused reference from a reference offset.
	 * Reference N padding increments a separate N count; non-N differences
	 * increment the substitution count.
	 * @param query query bases
	 * @param ref fused reference bases
	 * @param maxSubs reject when substitutions exceed this budget
	 * @param maxNs reject when N padding exceeds maxNs+maxSubs
	 * @param rStart starting offset in {@code ref}
	 * @return substitutions, or {@code maxSubs+1} if subs>maxSubs or ns>maxNs+maxSubs */
	static int alignFusedSIMD(byte[] query, byte[] ref, final int maxSubs,
			final int maxNs, final int rStart){
		int subs=0;
		int ns=0;
		int i=0;
		int j=rStart;
		final int qlen=query.length;
		final long limitNs=(long)maxNs+maxSubs;

		// 128-bit loop
		ByteVector vN128=ByteVector.broadcast(BSPECIES128, (byte)'N');
		for(final int lim=qlen-BWIDTH128; i<=lim; i+=BWIDTH128, j+=BWIDTH128){
			ByteVector qVec=ByteVector.fromArray(BSPECIES128, query, i);
			ByteVector rVec=ByteVector.fromArray(BSPECIES128, ref, j);

			// Identify Ns (Padding) vs Real Mismatches
			VectorMask<Byte> maskN=rVec.eq(vN128);
			VectorMask<Byte> maskDiff=qVec.eq(rVec).not();
			VectorMask<Byte> maskSub=maskDiff.and(maskN.not());

			ns+=maskN.trueCount();
			subs+=maskSub.trueCount();

			if(subs>maxSubs || ns>limitNs){return maxSubs+1;}
		}

		// 64-bit loop
		ByteVector vN64=ByteVector.broadcast(BSPECIES64, (byte)'N');
		for(final int lim=qlen-BWIDTH64; i<=lim; i+=BWIDTH64, j+=BWIDTH64){
			ByteVector qVec=ByteVector.fromArray(BSPECIES64, query, i);
			ByteVector rVec=ByteVector.fromArray(BSPECIES64, ref, j);

			VectorMask<Byte> maskN=rVec.eq(vN64);
			VectorMask<Byte> maskDiff=qVec.eq(rVec).not();
			VectorMask<Byte> maskSub=maskDiff.and(maskN.not());

			ns+=maskN.trueCount();
			subs+=maskSub.trueCount();

			if(subs>maxSubs || ns>limitNs){return maxSubs+1;}
		}

		// Scalar tail
		for(; i<qlen; i++, j++){
			byte q=query[i];
			byte r=ref[j];
			if(r=='N'){ns++;}
			else if(q!=r){subs++;}
		}

		if(subs>maxSubs || ns>limitNs){return maxSubs+1;}
		return subs;
	}

	/** Checks if all bytes are ASCII nucleotide letters.
	 * Allows IUPAC nucleotide letters including X/x and rejects amino-only
	 * letters E/F/I/J/L/O/P/Q/Z, case-insensitively.
	 * @param array byte array to check; null is accepted
	 * @return true if every byte is an accepted nucleotide letter */
	static final boolean isNucleotide(final byte[] array){
		if(array==null){return true;}

		final int length=array.length;
		if(length<BWIDTH128){
			boolean success=true;
			for(int i=0; i<length; i++){
				success&=(AminoAcid.baseToNumberExtended[array[i]]>=0);
			}
			return success;
		}

		int i=0;
		VectorMask<Byte> invalid=BSPECIES128.maskAll(false);
		final ByteVector caseMask=ByteVector.broadcast(BSPECIES128, (byte)~32);
		final ByteVector indexMask=ByteVector.broadcast(BSPECIES128, (byte)15);

		for(; i<=length-BWIDTH128; i+=BWIDTH128){
			final ByteVector folded=ByteVector.fromArray(BSPECIES128, array, i).and(caseMask);
			final ByteVector indexes=folded.and(indexMask);
			final VectorMask<Byte> valid=indexes.selectFrom(NUCLEOTIDE_LOOKUP_LOW_128).eq(folded)
					.or(indexes.selectFrom(NUCLEOTIDE_LOOKUP_HIGH_128).eq(folded));
			invalid=valid.not().or(invalid);
		}

		if(i<length){
			//Overlapping the final 16 bytes is safe for a read-only Boolean predicate.
			final ByteVector folded=ByteVector.fromArray(BSPECIES128, array,
					length-BWIDTH128).and(caseMask);
			final ByteVector indexes=folded.and(indexMask);
			final VectorMask<Byte> valid=indexes.selectFrom(NUCLEOTIDE_LOOKUP_LOW_128).eq(folded)
					.or(indexes.selectFrom(NUCLEOTIDE_LOOKUP_HIGH_128).eq(folded));
			invalid=valid.not().or(invalid);
		}

		return !invalid.anyTrue();
	}

	/**
	 * Appends each quality byte plus {@code delta} to {@code bb}.
	 * Bytes are cast after addition, matching the scalar FASTQ quality path.
	 * @param quals source quality bytes
	 * @param bb destination builder; capacity grows and length advances by {@code quals.length}
	 * @param delta value to add to each source byte before the final byte cast
	 */
	static void addAndAppend(byte[] quals, ByteBuilder bb, int delta){
		final int qlen=quals.length;
		bb.ensureExtra(qlen);
		final byte[] array=bb.array;
		final int offset=bb.length;

		int i=0;

		{//128-bit loop
			final ByteVector vDelta128=ByteVector.broadcast(BSPECIES128, (byte)delta);
			for(; i<=qlen-BWIDTH128; i+=BWIDTH128){
				ByteVector vq=ByteVector.fromArray(BSPECIES128, quals, i);
				ByteVector vResult=vq.add(vDelta128);
				vResult.intoArray(array, offset+i);
			}
		}

		{//64-bit loop
			final ByteVector vDelta64=ByteVector.broadcast(BSPECIES64, (byte)delta);
			for(; i<=qlen-BWIDTH64; i+=BWIDTH64){
				ByteVector vq=ByteVector.fromArray(BSPECIES64, quals, i);
				ByteVector vResult=vq.add(vDelta64);
				vResult.intoArray(array, offset+i);
			}
		}

		//Scalar tail
		for(; i<qlen; i++){
			array[offset+i]=(byte)(quals[i]+delta);
		}

		bb.length+=qlen;
	}

	/**
	 * Appends {@code qFake} for A/C/G/T/U bases, case-insensitively, and {@code qUndef} for all other bytes.
	 * @param bases source bases to classify
	 * @param bb destination builder; capacity grows and length advances by {@code bases.length}
	 * @param qFake byte value to append for fully defined bases
	 * @param qUndef byte value to append for ambiguous or invalid bases
	 */
	static void appendFake(byte[] bases, ByteBuilder bb, int qFake, int qUndef){
		final int blen=bases.length;
		bb.ensureExtra(blen);
		final byte[] array=bb.array;
		final int offset=bb.length;

		final byte mask=~32; //Uppercase mask

		int i=0;

		{//128-bit loop
			final ByteVector vQFake128=ByteVector.broadcast(BSPECIES128, (byte)qFake);
			final ByteVector vQUndef128=ByteVector.broadcast(BSPECIES128, (byte)qUndef);
			final ByteVector vmask128=ByteVector.broadcast(BSPECIES128, mask);
			final ByteVector vA128=ByteVector.broadcast(BSPECIES128, (byte)'A');
			final ByteVector vC128=ByteVector.broadcast(BSPECIES128, (byte)'C');
			final ByteVector vG128=ByteVector.broadcast(BSPECIES128, (byte)'G');
			final ByteVector vT128=ByteVector.broadcast(BSPECIES128, (byte)'T');
			final ByteVector vU128=ByteVector.broadcast(BSPECIES128, (byte)'U');

			for(; i<=blen-BWIDTH128; i+=BWIDTH128){
				ByteVector vBases=ByteVector.fromArray(BSPECIES128, bases, i);
				ByteVector vb=vBases.and(vmask128);
				VectorMask<Byte> isA=vb.eq(vA128);
				VectorMask<Byte> isC=vb.eq(vC128);
				VectorMask<Byte> isG=vb.eq(vG128);
				VectorMask<Byte> isT=vb.eq(vT128);
				VectorMask<Byte> isU=vb.eq(vU128);
				VectorMask<Byte> isDefined=isA.or(isC).or(isG).or(isT).or(isU);
				ByteVector vResult=vQFake128.blend(vQUndef128, isDefined.not());
				vResult.intoArray(array, offset+i);
			}
		}

		{//64-bit loop
			final ByteVector vQFake64=ByteVector.broadcast(BSPECIES64, (byte)qFake);
			final ByteVector vQUndef64=ByteVector.broadcast(BSPECIES64, (byte)qUndef);
			final ByteVector vmask64=ByteVector.broadcast(BSPECIES64, mask);
			final ByteVector vA64=ByteVector.broadcast(BSPECIES64, (byte)'A');
			final ByteVector vC64=ByteVector.broadcast(BSPECIES64, (byte)'C');
			final ByteVector vG64=ByteVector.broadcast(BSPECIES64, (byte)'G');
			final ByteVector vT64=ByteVector.broadcast(BSPECIES64, (byte)'T');
			final ByteVector vU64=ByteVector.broadcast(BSPECIES64, (byte)'U');

			for(; i<=blen-BWIDTH64; i+=BWIDTH64){
				ByteVector vBases=ByteVector.fromArray(BSPECIES64, bases, i);
				ByteVector vb=vBases.and(vmask64);
				VectorMask<Byte> isA=vb.eq(vA64);
				VectorMask<Byte> isC=vb.eq(vC64);
				VectorMask<Byte> isG=vb.eq(vG64);
				VectorMask<Byte> isT=vb.eq(vT64);
				VectorMask<Byte> isU=vb.eq(vU64);
				VectorMask<Byte> isDefined=isA.or(isC).or(isG).or(isT).or(isU);
				ByteVector vResult=vQFake64.blend(vQUndef64, isDefined.not());
				vResult.intoArray(array, offset+i);
			}
		}

		//Scalar tail
		for(; i<blen; i++){
			array[offset+i]=(byte)(AminoAcid.isFullyDefined(bases[i]) ? qFake : qUndef);
		}

		bb.length+=blen;
	}

	/**
	 * Appends qualities in reverse order after adding {@code delta} to each byte.
	 * Bytes are cast after addition, matching {@link #addAndAppend(byte[], ByteBuilder, int)}.
	 * @param quals source quality bytes
	 * @param bb destination builder; capacity grows and length advances by {@code quals.length}
	 * @param delta value to add to each source byte before the final byte cast
	 */
	static void addAndAppendReversed(byte[] quals, ByteBuilder bb, int delta){
		final int qlen=quals.length;
		bb.ensureExtra(qlen);
		final byte[] array=bb.array;
		final int bufferStart=bb.length;
		final int bufferStop=bb.length+qlen;

		int i=bufferStart;
		int j=qlen; // Position after last unused quality score
		final int limit128=bufferStart+(qlen/BWIDTH128)*BWIDTH128;
		final int limit64=bufferStart+(qlen/BWIDTH64)*BWIDTH64;

		{//128-bit loop
			final ByteVector vDelta128=ByteVector.broadcast(BSPECIES128, (byte)delta);

			for(; i<limit128; i+=BWIDTH128, j-=BWIDTH128){
				ByteVector vq=ByteVector.fromArray(BSPECIES128, quals, j-BWIDTH128);
				ByteVector vAdded=vq.add(vDelta128);
				ByteVector vReversed=vAdded.rearrange(B_REVERSE_SHUFFLE_128);
				vReversed.intoArray(array, i);
			}
		}

		{//64-bit loop
			final ByteVector vDelta64=ByteVector.broadcast(BSPECIES64, (byte)delta);

			for(; i<limit64; i+=BWIDTH64, j-=BWIDTH64){
				ByteVector vq=ByteVector.fromArray(BSPECIES64, quals, j-BWIDTH64);
				ByteVector vAdded=vq.add(vDelta64);
				ByteVector vReversed=vAdded.rearrange(B_REVERSE_SHUFFLE_64);
				vReversed.intoArray(array, i);
			}
		}

		//Scalar tail
		for(j--; i<bufferStop; i++, j--){
			array[i]=(byte)(quals[j]+delta);
		}

		bb.length+=qlen;
	}

	/**
	 * Reverses the first {@code len} bytes in place.
	 * @param array array to mutate; null is ignored
	 * @param len number of leading bytes to reverse
	 */
	static void reverseInPlace(final byte[] array, final int len){
		if(array==null || len<2){return;}

		int left=0;
		int right=len;

		{//128-bit loop - work from both ends toward middle
			for(; left+BWIDTH128<=len/2; left+=BWIDTH128, right-=BWIDTH128){
				//Load from both ends
				ByteVector vLeft=ByteVector.fromArray(BSPECIES128, array, left);
				ByteVector vRight=ByteVector.fromArray(BSPECIES128, array, right-BWIDTH128);

				//Reverse each vector
				ByteVector vLeftRev=vLeft.rearrange(B_REVERSE_SHUFFLE_128);
				ByteVector vRightRev=vRight.rearrange(B_REVERSE_SHUFFLE_128);

				//Write to opposite ends
				vRightRev.intoArray(array, left);
				vLeftRev.intoArray(array, right-BWIDTH128);
			}
		}

		{//64-bit loop
			for(; left+BWIDTH64<=len/2; left+=BWIDTH64, right-=BWIDTH64){
				ByteVector vLeft=ByteVector.fromArray(BSPECIES64, array, left);
				ByteVector vRight=ByteVector.fromArray(BSPECIES64, array, right-BWIDTH64);

				ByteVector vLeftRev=vLeft.rearrange(B_REVERSE_SHUFFLE_64);
				ByteVector vRightRev=vRight.rearrange(B_REVERSE_SHUFFLE_64);

				vRightRev.intoArray(array, left);
				vLeftRev.intoArray(array, right-BWIDTH64);
			}
		}

		//Scalar tail - swap remaining elements
		right--;
		while(left<right){
			byte temp=array[left];
			array[left]=array[right];
			array[right]=temp;
			left++;
			right--;
		}
	}

	/**
	 * Reverse-complements the first {@code len} bytes in place.
	 * A/C/G/T/N blocks are vectorized; the first non-ACGTN block transfers the
	 * remaining inward region to the scalar path so IUPAC bases receive their
	 * proper {@link AminoAcid#baseToComplementExtended} complements.
	 * @param array array to mutate; null is ignored
	 * @param len number of leading bases to reverse-complement
	 */
	static void reverseComplementInPlace(final byte[] array, final int len){
		if(array==null || len<1){return;}
		int left=0;
		int right=len;
		final byte caseMask=32; // 0x20 - difference between upper and lowercase
		boolean clean=true;//FIXED [simd/SIMDByte256#001], inherited: per-block non-ACGTN detection; on the first impure block, bail to the proper-complement scalar tail instead of coercing IUPAC to N (which diverged from the scalar/SIMD-off path).

		{//128-bit loop
			final ByteVector vCaseMask128=ByteVector.broadcast(BSPECIES128, caseMask);
			final ByteVector vA128=ByteVector.broadcast(BSPECIES128, (byte)'A');
			final ByteVector vC128=ByteVector.broadcast(BSPECIES128, (byte)'C');
			final ByteVector vG128=ByteVector.broadcast(BSPECIES128, (byte)'G');
			final ByteVector vT128=ByteVector.broadcast(BSPECIES128, (byte)'T');
			final ByteVector vN128=ByteVector.broadcast(BSPECIES128, (byte)'N');

			for(; left+BWIDTH128<=len/2; left+=BWIDTH128, right-=BWIDTH128){
				ByteVector vLeft=ByteVector.fromArray(BSPECIES128, array, left);
				ByteVector vRight=ByteVector.fromArray(BSPECIES128, array, right-BWIDTH128);

				//Process left vector
				ByteVector casesLeft=vLeft.and(vCaseMask128);
				ByteVector upperLeft=vLeft.and(vCaseMask128.not());
				VectorMask<Byte> isA_L=upperLeft.eq(vA128);
				VectorMask<Byte> isC_L=upperLeft.eq(vC128);
				VectorMask<Byte> isG_L=upperLeft.eq(vG128);
				VectorMask<Byte> isT_L=upperLeft.eq(vT128);
				ByteVector compLeft=vN128.blend(vT128, isA_L).blend(vG128, isC_L).blend(vC128, isG_L).blend(vA128, isT_L);
				compLeft=compLeft.or(casesLeft);
				ByteVector vLeftRevComp=compLeft.rearrange(B_REVERSE_SHUFFLE_128);

				//Process right vector
				ByteVector casesRight=vRight.and(vCaseMask128);
				ByteVector upperRight=vRight.and(vCaseMask128.not());
				VectorMask<Byte> isA_R=upperRight.eq(vA128);
				VectorMask<Byte> isC_R=upperRight.eq(vC128);
				VectorMask<Byte> isG_R=upperRight.eq(vG128);
				VectorMask<Byte> isT_R=upperRight.eq(vT128);
				ByteVector compRight=vN128.blend(vT128, isA_R).blend(vG128, isC_R).blend(vC128, isG_R).blend(vA128, isT_R);
				compRight=compRight.or(casesRight);
				ByteVector vRightRevComp=compRight.rearrange(B_REVERSE_SHUFFLE_128);

				//FIXED [simd/SIMDByte256#001], inherited: bail to the scalar tail (proper IUPAC complement) if either block holds a non-ACGTN.
				//isA/isC/isG/isT are already computed for the blend above; only isN is extra, so detection is near-free for the common all-ACGTN case.
				VectorMask<Byte> validL=isA_L.or(isC_L).or(isG_L).or(isT_L).or(upperLeft.eq(vN128));
				VectorMask<Byte> validR=isA_R.or(isC_R).or(isG_R).or(isT_R).or(upperRight.eq(vN128));
				if(!validL.allTrue() || !validR.allTrue()){clean=false; break;}

				//Write to opposite ends
				vRightRevComp.intoArray(array, left);
				vLeftRevComp.intoArray(array, right-BWIDTH128);
			}
		}

		if(clean){//64-bit loop (skipped once a non-ACGTN block has sent us to the scalar tail)
			final ByteVector vCaseMask64=ByteVector.broadcast(BSPECIES64, caseMask);
			final ByteVector vA64=ByteVector.broadcast(BSPECIES64, (byte)'A');
			final ByteVector vC64=ByteVector.broadcast(BSPECIES64, (byte)'C');
			final ByteVector vG64=ByteVector.broadcast(BSPECIES64, (byte)'G');
			final ByteVector vT64=ByteVector.broadcast(BSPECIES64, (byte)'T');
			final ByteVector vN64=ByteVector.broadcast(BSPECIES64, (byte)'N');

			for(; left+BWIDTH64<=len/2; left+=BWIDTH64, right-=BWIDTH64){
				ByteVector vLeft=ByteVector.fromArray(BSPECIES64, array, left);
				ByteVector vRight=ByteVector.fromArray(BSPECIES64, array, right-BWIDTH64);

				//Process left
				ByteVector casesLeft=vLeft.and(vCaseMask64);
				ByteVector upperLeft=vLeft.and(vCaseMask64.not());
				VectorMask<Byte> isA_L=upperLeft.eq(vA64);
				VectorMask<Byte> isC_L=upperLeft.eq(vC64);
				VectorMask<Byte> isG_L=upperLeft.eq(vG64);
				VectorMask<Byte> isT_L=upperLeft.eq(vT64);
				ByteVector compLeft=vN64.blend(vT64, isA_L).blend(vG64, isC_L).blend(vC64, isG_L).blend(vA64, isT_L);
				compLeft=compLeft.or(casesLeft);
				ByteVector vLeftRevComp=compLeft.rearrange(B_REVERSE_SHUFFLE_64);

				//Process right
				ByteVector casesRight=vRight.and(vCaseMask64);
				ByteVector upperRight=vRight.and(vCaseMask64.not());
				VectorMask<Byte> isA_R=upperRight.eq(vA64);
				VectorMask<Byte> isC_R=upperRight.eq(vC64);
				VectorMask<Byte> isG_R=upperRight.eq(vG64);
				VectorMask<Byte> isT_R=upperRight.eq(vT64);
				ByteVector compRight=vN64.blend(vT64, isA_R).blend(vG64, isC_R).blend(vC64, isG_R).blend(vA64, isT_R);
				compRight=compRight.or(casesRight);
				ByteVector vRightRevComp=compRight.rearrange(B_REVERSE_SHUFFLE_64);

				//FIXED [simd/SIMDByte256#001], inherited: same per-block non-ACGTN bail as the 128-bit loop.
				VectorMask<Byte> validL=isA_L.or(isC_L).or(isG_L).or(isT_L).or(upperLeft.eq(vN64));
				VectorMask<Byte> validR=isA_R.or(isC_R).or(isG_R).or(isT_R).or(upperRight.eq(vN64));
				if(!validL.allTrue() || !validR.allTrue()){clean=false; break;}

				vRightRevComp.intoArray(array, left);
				vLeftRevComp.intoArray(array, right-BWIDTH64);
			}
		}

		//FIXED [simd/SIMDByte256#001], inherited: this scalar tail (proper IUPAC complement via baseToComplementExtended) now also serves as the fallback for any block the SIMD loops flagged as non-ACGTN (clean==false).  Previously the SIMD loops coerced IUPAC->N while this tail complemented them, so the result was position-dependent and diverged from the SIMD-off path; now the whole array is proper-complemented (SIMD for the ACGTN prefix/suffix, scalar from the first impure block inward).
		//Scalar tail
		right--;
		while(left<right){
			byte bLeft=array[left];
			byte bRight=array[right];
			array[left]=AminoAcid.baseToComplementExtended[bRight];
			array[right]=AminoAcid.baseToComplementExtended[bLeft];
			left++;
			right--;
		}
		//Handle middle element if odd length
		if(left==right){
			array[left]=AminoAcid.baseToComplementExtended[array[left]];
		}
	}

	/**
	 * Checks whether the first {@code length} bytes are ASCII A/C/G/T/N, case-insensitively.
	 * @param array byte array to check; null is accepted
	 * @param length number of leading bytes to inspect
	 * @return true if {@code array} is null, {@code length} is zero, or every inspected byte is A/C/G/T/N
	 */
	static boolean isACGTN(byte[] array, int length){
		if(array==null || length==0){return true;}

		int i=0;
		final byte UPPERMASK=(byte)(~32);

		{//128-bit loop
			VectorMask<Byte> valid=BSPECIES128.maskAll(true);
			final ByteVector vA128=ByteVector.broadcast(BSPECIES128, (byte)'A');
			final ByteVector vC128=ByteVector.broadcast(BSPECIES128, (byte)'C');
			final ByteVector vG128=ByteVector.broadcast(BSPECIES128, (byte)'G');
			final ByteVector vT128=ByteVector.broadcast(BSPECIES128, (byte)'T');
			final ByteVector vN128=ByteVector.broadcast(BSPECIES128, (byte)'N');
			final ByteVector vUPPER128=ByteVector.broadcast(BSPECIES128, UPPERMASK);

			for(; i<=length-BWIDTH128; i+=BWIDTH128){
				ByteVector v=ByteVector.fromArray(BSPECIES128, array, i).and(vUPPER128);
				VectorMask<Byte> isA=v.eq(vA128);
				VectorMask<Byte> isC=v.eq(vC128);
				VectorMask<Byte> isG=v.eq(vG128);
				VectorMask<Byte> isT=v.eq(vT128);
				VectorMask<Byte> isN=v.eq(vN128);
				VectorMask<Byte> isAny=isA.or(isC).or(isG).or(isT).or(isN);
				valid=valid.and(isAny);
			}

			if(!valid.allTrue()){return false;}
		}

		{//64-bit loop
			VectorMask<Byte> valid=BSPECIES64.maskAll(true);
			final ByteVector vA64=ByteVector.broadcast(BSPECIES64, (byte)'A');
			final ByteVector vC64=ByteVector.broadcast(BSPECIES64, (byte)'C');
			final ByteVector vG64=ByteVector.broadcast(BSPECIES64, (byte)'G');
			final ByteVector vT64=ByteVector.broadcast(BSPECIES64, (byte)'T');
			final ByteVector vN64=ByteVector.broadcast(BSPECIES64, (byte)'N');
			final ByteVector vUPPER64=ByteVector.broadcast(BSPECIES64, UPPERMASK);

			for(; i<=length-BWIDTH64; i+=BWIDTH64){
				ByteVector v=ByteVector.fromArray(BSPECIES64, array, i).and(vUPPER64);
				VectorMask<Byte> isA=v.eq(vA64);
				VectorMask<Byte> isC=v.eq(vC64);
				VectorMask<Byte> isG=v.eq(vG64);
				VectorMask<Byte> isT=v.eq(vT64);
				VectorMask<Byte> isN=v.eq(vN64);
				VectorMask<Byte> isAny=isA.or(isC).or(isG).or(isT).or(isN);
				valid=valid.and(isAny);
			}

			if(!valid.allTrue()){return false;}
		}

		//Scalar tail
		for(; i<length; i++){
			byte b=array[i];
			if(!ACGTN[b]){return false;}
		}

		return true;
	}

	/** Reverse-complements a binary kmer with SIMD byte reversal.
	 * @param kmer Binary kmer
	 * @param k Kmer length
	 * @return Reverse-complemented kmer
	 */
	public static final long reverseComplementBinaryFastSIMD(long kmer, int k){
		// Swap adjacent 2-bit pairs
		long x=((kmer&0x3333333333333333L)<<2)|((kmer&0xCCCCCCCCCCCCCCCCL)>>2);

		// Swap nibbles
		x=((x&0x0F0F0F0F0F0F0F0FL)<<4)|((x&0xF0F0F0F0F0F0F0F0L)>>4);

		// Complement
		x=~x;

		// Reverse byte order with SIMD
		LongVector lvec=LongVector.broadcast(LSPECIES64, x);
		ByteVector bvec=lvec.reinterpretAsBytes();
		bvec=bvec.rearrange(B_REVERSE_SHUFFLE_64);
		lvec=bvec.reinterpretAsLongs();
		x=lvec.lane(0);

		// Right-align for k<32
		if(k<32){x=x>>>(2*(32-k));}

		return x;
	}

	/** Reverse-complements a binary kmer with fixed bitwise swaps.
	 * @param kmer Binary kmer
	 * @param k Kmer length
	 * @return Reverse-complemented kmer
	 */
	public static final long reverseComplementBinaryFastBitwise(long kmer, int k){
		// Complement first
		long x=~kmer;

		// Swap adjacent 2-bit pairs
		x=((x&0x3333333333333333L)<<2)|((x&0xCCCCCCCCCCCCCCCCL)>>>2);

		// Swap nibbles
		x=((x&0x0F0F0F0F0F0F0F0FL)<<4)|((x&0xF0F0F0F0F0F0F0F0L)>>>4);

		// Swap bytes
		x=((x&0x00FF00FF00FF00FFL)<<8)|((x&0xFF00FF00FF00FF00L)>>>8);

		// Swap 16-bit words
		x=((x&0x0000FFFF0000FFFFL)<<16)|((x&0xFFFF0000FFFF0000L)>>>16);

		// Swap 32-bit dwords
		x=(x<<32)|(x>>>32);

		// Right-align for k<32
		x=x>>>(2*(32-k));

		return x;
	}

	/** Reverse-complements a binary kmer with looped bitwise swaps.
	 * @param kmer Binary kmer
	 * @param k Kmer length
	 * @return Reverse-complemented kmer
	 */
	public static final long reverseComplementBinaryFastBitwise2(long kmer, int k){
		long x=~kmer;
		int kmerBits=2*k;
		int swapWidth=2;

		for(int i=0; i<SWAP_MASKS.length && swapWidth<=kmerBits; i++){
			int shift=1<<i; // 1, 2, 4, 8, 16
			x=((x&SWAP_MASKS[i])<<shift)|((x&~SWAP_MASKS[i])>>>shift);
			swapWidth<<=1;
		}

		if(k<32){x=x>>>(2*(32-k));}
		return x;
	}

	/** Reverse-complements a 32-bit binary kmer.
	 * @param kmer Binary kmer
	 * @param k Kmer length
	 * @return Reverse-complemented kmer
	 */
	public static final int reverseComplementBinaryFastBitwise(int kmer, int k){
		// Complement first
		int x=~kmer;

		// Swap adjacent 2-bit pairs
		x=((x&0x33333333)<<2)|((x&0xCCCCCCCC)>>>2);

		// Swap nibbles
		x=((x&0x0F0F0F0F)<<4)|((x&0xF0F0F0F0)>>>4);

		// Swap bytes
		x=((x&0x00FF00FF)<<8)|((x&0xFF00FF00)>>>8);

		// Swap 16-bit words
		x=(x<<16)|(x>>>16);

		// Right-align for k<16
		x=x>>>(2*(16-k));

		return x;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Preferred byte-vector species. */
	private static final VectorSpecies<Byte> BSPECIES128=ByteVector.SPECIES_128;
	/** Number of byte lanes in BSPECIES128. */
	private static final int BWIDTH128=16;

	/** Tail byte-vector species. */
	private static final VectorSpecies<Byte> BSPECIES64=ByteVector.SPECIES_64;
	/** Number of byte lanes in BSPECIES64. */
	private static final int BWIDTH64=8;

	/** Integer-vector species for widening byte reductions. */
	private static final VectorSpecies<Integer> ISPECIES128=IntVector.SPECIES_128;

	/** Single-lane long species for reinterpret-based byte shuffles. */
	private static final VectorSpecies<Long> LSPECIES64=LongVector.SPECIES_64;
	/** Number of long lanes in LSPECIES64. */
	private static final int LWIDTH64=1;

	/** 128-bit byte reverse shuffle. */
	private static final VectorShuffle<Byte> B_REVERSE_SHUFFLE_128;
	/** 64-bit byte reverse shuffle. */
	private static final VectorShuffle<Byte> B_REVERSE_SHUFFLE_64;
	/** True for ASCII A/C/G/T/N in either case. */
	private static final boolean[] ACGTN;
	/** Low-half lookup table for accepted IUPAC nucleotide ASCII letters. */
	private static final ByteVector NUCLEOTIDE_LOOKUP_LOW_128;
	/** High-half lookup table for accepted IUPAC nucleotide ASCII letters. */
	private static final ByteVector NUCLEOTIDE_LOOKUP_HIGH_128;

	static{
		assert(BSPECIES128.length()==BWIDTH128);
		assert(BSPECIES64.length()==BWIDTH64);
		assert(LSPECIES64.length()==LWIDTH64);

		//128-bit reverse shuffle
		int vlen128=BSPECIES128.length();
		int[] indices128=new int[vlen128];
		for(int i=0; i<vlen128; i++){
			indices128[i]=vlen128-1-i;
		}
		B_REVERSE_SHUFFLE_128=VectorShuffle.fromArray(BSPECIES128, indices128, 0);

		//64-bit reverse shuffle
		int vlen64=BSPECIES64.length();
		int[] indices64=new int[vlen64];
		for(int i=0; i<vlen64; i++){
			indices64[i]=vlen64-1-i;
		}
		B_REVERSE_SHUFFLE_64=VectorShuffle.fromArray(BSPECIES64, indices64, 0);
		ACGTN=new boolean[128];
		Arrays.fill(ACGTN, false);
		ACGTN['A']=ACGTN['C']=ACGTN['G']=ACGTN['T']=ACGTN['N']=true;
		ACGTN['a']=ACGTN['c']=ACGTN['g']=ACGTN['t']=ACGTN['n']=true;

		NUCLEOTIDE_LOOKUP_LOW_128=makeNucleotideLookup(0);
		NUCLEOTIDE_LOOKUP_HIGH_128=makeNucleotideLookup(16);
	}

	/** Masks for progressively wider binary-kmer swaps. */
	private static final long[] SWAP_MASKS={
		0x3333333333333333L, // 2-bit pairs
		0x0F0F0F0F0F0F0F0FL, // nibbles
		0x00FF00FF00FF00FFL, // bytes
		0x0000FFFF0000FFFFL, // 16-bit words
		0x00000000FFFFFFFFL  // 32-bit dwords
	};

}
