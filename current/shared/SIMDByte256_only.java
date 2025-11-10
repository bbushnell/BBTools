package shared;

import dna.AminoAcid;
import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorShuffle;
import jdk.incubator.vector.VectorSpecies;
import structures.ByteBuilder;
import structures.IntList;

/**
 * Holds SIMD methods.
 * @author Brian Bushnell
 * @date Sep 12, 2023?
 *
 */
final class SIMDByte256_only{

	private static final VectorSpecies<Byte> BSPECIES=ByteVector.SPECIES_256;
	private static final int BWIDTH=BSPECIES.length();
	
	private static final VectorShuffle<Byte> B_REVERSE_SHUFFLE;
	static {
		int vlen=BSPECIES.length();
		int[] indices=new int[vlen];
		for(int i=0; i<vlen; i++){
			indices[i]=vlen-1-i;
		}
		B_REVERSE_SHUFFLE=VectorShuffle.fromArray(BSPECIES, indices, 0);
	}
	
	/** Returns number of matches */
	static final int countMatches(final byte[] s1, final byte[] s2, int a1, int b1, int a2, int b2){
		final int length=b2-a2+1;
		final int limit0=BSPECIES.loopBound(length);
		final int limit=a2+limit0;

		int i=a1, j=a2;
		int matches=0;
		for(; j<limit; i+=BWIDTH, j+=BWIDTH){// SIMD loop
			ByteVector v1=ByteVector.fromArray(BSPECIES, s1, i);
			ByteVector v2=ByteVector.fromArray(BSPECIES, s2, j);
			VectorMask<Byte> x=v1.eq(v2);
			matches+=x.trueCount();// This might be slow, or might not
		}
		for(; j<=b2; i++, j++){
			final byte x=s1[i], y=s2[j];
			final int m=((x==y) ? 1 : 0);
			matches+=m;
		}
		return matches;
	}

	/** Returns index of symbol */
	static final int find(final byte[] a, final byte symbol, final int from, final int to){// 15% Slower than scalar code, at least for ByteFile1
		final int length=to-from;// Intentionally exclusive
		final int limit0=BSPECIES.loopBound(length);
		final int limit=from+limit0;

		int pos=from;
		for(; pos<limit; pos+=BWIDTH){// SIMD loop
			ByteVector v=ByteVector.fromArray(BSPECIES, a, pos);
			VectorMask<Byte> x=v.eq(symbol);
			int t=x.firstTrue();
			if(t<BWIDTH){ return pos+t; }
			//			if(x.anyTrue()) {break;}
		}
		while(pos<to && a[pos]!=symbol){ pos++; }
		return pos;
	}

	
	/**
	 * Sums the array.
	 * @param a A vector.
	 * @return The sum.
	 */
	static final long sum(final byte[] a, final int from, final int to){// Tested as 4x scalar speed
		// TODO: Test speed.
		final int length=to-from+1;
		final int limit0=BSPECIES.loopBound(length);
		final int limit=from+limit0;

		int i=from;
		long c=0;
		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector va=ByteVector.fromArray(BSPECIES, a, i);
			c+=va.reduceLanesToLong(VectorOperators.ADD);
		}
		for(; i<=to; i++){ c+=a[i]; }// Residual scalar loop
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
		final int limit=BSPECIES.loopBound(to-from);
		final ByteVector newlineVec=ByteVector.broadcast(BSPECIES, symbol);

		int i=from;

		// SIMD loop
		for(; i<from+limit; i+=BWIDTH){
			ByteVector vec=ByteVector.fromArray(BSPECIES, buffer, i);
			VectorMask<Byte> mask=vec.eq(newlineVec);
			//			if(!mask.anyTrue()) {continue;}//Hopefully common case - Not faster, maybe 1% slower
			// Convert mask to long bitmask
			long bits=mask.toLong();

			// Extract set bit positions using bit manipulation
			//Brian version
			//			for(int lane=0; bits!=0; lane++){//Looks strange but lane needs to be incremented
			//				int zeros=Long.numberOfTrailingZeros(bits);
			//				lane+=zeros;
			//				positions.add(lane+i);
			//				bits>>>=(zeros+1);
			//			}
			while(bits!=0){//Isla version - 5% faster
				int lane=Long.numberOfTrailingZeros(bits);
				positions.add(i+lane);
				bits&=(bits-1); // Clear lowest set bit
			}
		}

		// Residual scalar loop
		for(; i<to; i++){
			if(buffer[i]==symbol){positions.add(i);}
		}

		return positions.size();
	}
	
	/**
	 * Find the last symbol in buffer by scanning backwards using SIMD.
	 * @param buffer Buffer to scan
	 * @param limit Scan backwards from this position (exclusive)
	 * @return Position of last newline, or -1 if none found
	 */
	static int findLastSymbol(byte[] buffer, int limit, final byte symbol){
		final ByteVector newlineVec=ByteVector.broadcast(BSPECIES, symbol);
		
		// Start from the last aligned chunk and work backwards
		int i=((limit-1)/BWIDTH)*BWIDTH; // Round down to vector boundary
		
		// SIMD loop - work backwards
		for(; i>=0; i-=BWIDTH){
			ByteVector vec=ByteVector.fromArray(BSPECIES, buffer, i);
			VectorMask<Byte> mask=vec.eq(newlineVec);
			long bits=mask.toLong();
			
			if(bits!=0){
				// Found newline(s) in this chunk - find the highest bit
				int lane=63-Long.numberOfLeadingZeros(bits); // Highest set bit
				return i+lane;
			}
		}
		
		// Residual scalar loop for beginning of buffer
		for(i+=BWIDTH-1; i>=0; i--){
			if(buffer[i]==symbol){return i;}
		}
		
		return -1;
	}

	
	/**
	 * Adds a constant delta to all bytes in an array using SIMD.
	 * @param array The byte array to modify in-place.
	 * @param delta The value to add to each element.
	 */
	static final void add(final byte[] array, final byte delta){
		if(array==null){return;}

		final int length=array.length;
		final int limit=BSPECIES.loopBound(length);

		int i=0;
		final ByteVector vdelta=ByteVector.broadcast(BSPECIES, delta);
		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector va=ByteVector.fromArray(BSPECIES, array, i);
			ByteVector vresult=va.add(vdelta);
			vresult.intoArray(array, i);
		}
		for(; i<length; i++){// Residual scalar loop
			array[i]+=delta;
		}
	}

	
	/**
	 * Adds delta to all bytes, caps at minimum value, and returns the minimum encountered.
	 * @param array The byte array to modify in-place.
	 * @param delta The value to add to each element.
	 * @param cap The minimum allowed value after addition.
	 * @return The minimum value encountered after addition (before capping).
	 */
	static final byte addAndCapMin(final byte[] array, final byte delta, final int cap){
		if(array==null){return 0;}

		final int length=array.length;
		final int limit=BSPECIES.loopBound(length);

		int i=0;
		ByteVector vdelta=ByteVector.broadcast(BSPECIES, delta);
		ByteVector vcap=ByteVector.broadcast(BSPECIES, (byte)cap);
		ByteVector vmin=ByteVector.broadcast(BSPECIES, (byte)127);

		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector va=ByteVector.fromArray(BSPECIES, array, i);
			ByteVector vresult=va.add(vdelta);
			vmin=vmin.min(vresult);
			ByteVector vcapped=vresult.max(vcap);
			vcapped.intoArray(array, i);
		}

		int min=vmin.reduceLanes(VectorOperators.MIN);

		for(; i<length; i++){// Residual scalar loop
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
		final int limit=BSPECIES.loopBound(length);

		int i=0;
		ByteVector vdelta=ByteVector.broadcast(BSPECIES, (byte)delta);
		ByteVector vn=ByteVector.broadcast(BSPECIES, (byte)'N');
		ByteVector vzero=ByteVector.broadcast(BSPECIES, (byte)0);
		ByteVector vcap2=ByteVector.broadcast(BSPECIES, (byte)2);

		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vquals=ByteVector.fromArray(BSPECIES, quals, i);
			ByteVector vbases=ByteVector.fromArray(BSPECIES, bases, i);

			// Add delta
			ByteVector vresult=vquals.add(vdelta);

			// Cap at 2 for non-N bases
			vresult=vresult.max(vcap2);

			// Create mask: where bases == 'N'
			VectorMask<Byte> maskN=vbases.eq(vn);

			// Blend: if N then 0, else capped result
			vresult=vresult.blend(vzero, maskN);

			vresult.intoArray(quals, i);
		}
		
		for(; i<length; i++){// Residual scalar loop
			byte b=bases[i];
			int q=quals[i]+delta;
			q=(AminoAcid.baseToNumber[b]<0 ? 0 : Math.max(2, q));
			quals[i]=(byte)q;
		}
	}

	
	/**
	 * Zeros quality for N bases, caps others at 2.
	 * @param quals Quality array to modify in-place.
	 * @param bases Sequence bases array.
	 */
	static final void capQuality(final byte[] quals, final byte[] bases){
		if(quals==null){return;}

		final int length=quals.length;
		final int limit=BSPECIES.loopBound(length);

		int i=0;
		ByteVector vn=ByteVector.broadcast(BSPECIES, (byte)'N');
		ByteVector vzero=ByteVector.broadcast(BSPECIES, (byte)0);
		ByteVector vcap2=ByteVector.broadcast(BSPECIES, (byte)2);

		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vquals=ByteVector.fromArray(BSPECIES, quals, i);
			ByteVector vbases=ByteVector.fromArray(BSPECIES, bases, i);

			// Cap at 2 for non-N bases
			ByteVector vresult=vquals.max(vcap2);

			// Create mask: where bases == 'N'
			VectorMask<Byte> maskN=vbases.eq(vn);

			// Blend: if N then 0, else capped result
			vresult=vresult.blend(vzero, maskN);

			vresult.intoArray(quals, i);
		}
		
		for(; i<length; i++){// Residual scalar loop
			byte b=bases[i];
			int q=quals[i];
			q=(AminoAcid.baseToNumber[b]<0 ? 0 : Math.max(2, q));
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
		final int limit=BSPECIES.loopBound(length);

		int i=0;
		ByteVector vU=ByteVector.broadcast(BSPECIES, (byte)'U');
		ByteVector vu=ByteVector.broadcast(BSPECIES, (byte)'u');
		ByteVector vT=ByteVector.broadcast(BSPECIES, (byte)'T');
		ByteVector vt=ByteVector.broadcast(BSPECIES, (byte)'t');

		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vbases=ByteVector.fromArray(BSPECIES, bases, i);

			// Create masks for U and u
			VectorMask<Byte> maskU=vbases.eq(vU);
			VectorMask<Byte> masku=vbases.eq(vu);

			// Replace U with T, u with t
			vbases=vbases.blend(vT, maskU);
			vbases=vbases.blend(vt, masku);

			vbases.intoArray(bases, i);
		}

		for(; i<length; i++){// Residual scalar loop
			bases[i]=AminoAcid.uToT[bases[i]];
		}
	}

	
	/**
	 * Converts lowercase letters to N.
	 * @param array The byte array to modify in-place.
	 * @return Always true.
	 */
	static final boolean lowerCaseToN(final byte[] array){
		if(array==null){return true;}

		final int length=array.length;
		final int limit=BSPECIES.loopBound(length);

		int i=0;
		final byte a='a', N='N';

		ByteVector va=ByteVector.broadcast(BSPECIES, a);
		ByteVector vN=ByteVector.broadcast(BSPECIES, N);

		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vb=ByteVector.fromArray(BSPECIES, array, i);

			// Create mask: where b > a (lowercase)
			VectorMask<Byte> maskLower=vb.compare(VectorOperators.GE, va);

			// Replace with N where lowercase
			ByteVector vresult=vb.blend(vN, maskLower);

			vresult.intoArray(array, i);
		}

		for(; i<length; i++){// Residual scalar loop
			array[i]=AminoAcid.lowerCaseToNocall[array[i]];
		}

		return true;
	}

	
	/**
	 * Converts dot, dash, and X to N.
	 * @param array The byte array to modify in-place.
	 * @return Always true.
	 */
	static final boolean dotDashXToN(final byte[] array){
		if(array==null){return true;}

		final int length=array.length;
		final int limit=BSPECIES.loopBound(length);

		int i=0;
		final byte dot='.', dash='-', X='X', N='N';

		ByteVector vdot=ByteVector.broadcast(BSPECIES, dot);
		ByteVector vdash=ByteVector.broadcast(BSPECIES, dash);
		ByteVector vX=ByteVector.broadcast(BSPECIES, X);
		ByteVector vN=ByteVector.broadcast(BSPECIES, N);

		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vb=ByteVector.fromArray(BSPECIES, array, i);

			// Create masks for each character
			VectorMask<Byte> maskDot=vb.eq(vdot);
			VectorMask<Byte> maskDash=vb.eq(vdash);
			VectorMask<Byte> maskX=vb.eq(vX);

			// Combine masks: any of the three
			VectorMask<Byte> maskAny=maskDot.or(maskDash).or(maskX);

			// Replace with N where mask is true
			ByteVector vresult=vb.blend(vN, maskAny);

			vresult.intoArray(array, i);
		}

		for(; i<length; i++){// Residual scalar loop
			array[i]=AminoAcid.dotDashXToNocall[array[i]];
		}

		return true;
	}
	
	
	/**
	 * Checks if array contains common amino acids (E or L).
	 * @param array The byte array to check.
	 * @return true if E or L found (likely protein).
	 */
	static final boolean isProtein(final byte[] array){
		if(array==null){return false;}
		
		final int length=array.length;
		final int limit=BSPECIES.loopBound(length);
		
		int i=0;
		final byte E='E', L='L';
		
		ByteVector vE=ByteVector.broadcast(BSPECIES, E);
		ByteVector vL=ByteVector.broadcast(BSPECIES, L);
		
		VectorMask<Byte> foundProtein=BSPECIES.maskAll(false);
		
		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vb=ByteVector.fromArray(BSPECIES, array, i);
			
			VectorMask<Byte> isE=vb.eq(vE);
			VectorMask<Byte> isL=vb.eq(vL);
			
			foundProtein=isE.or(isL).or(foundProtein);
		}
		
		boolean protein=foundProtein.anyTrue();
		
		for(; i<length; i++){// Residual scalar loop
			byte b=array[i];
			boolean nuc=AminoAcid.baseToNumberExtended[b]>=0;
			boolean amino=AminoAcid.acidToNumberExtended[b]>=0;
//			protein|=(b==E || b==L);
			protein|=(amino && !nuc);
		}
		
		return protein;
	}
	
	/* */
	
	
	/**
	 * Converts to uppercase using bitmask. Returns false if non-letter found.
	 * @param array The byte array to modify in-place.
	 * @return false if any byte is outside A-Z range after conversion.
	 */
	static final boolean toUpperCase(final byte[] array){
		if(array==null){return true;}
		
		final int length=array.length;
		final int limit=BSPECIES.loopBound(length);
		
		int i=0;
		final byte A='A', Z='Z';
		final byte mask=~32;
		
		ByteVector vmask=ByteVector.broadcast(BSPECIES, mask);
		ByteVector vA=ByteVector.broadcast(BSPECIES, A);
		ByteVector vZ=ByteVector.broadcast(BSPECIES, Z);
		
		VectorMask<Byte> invalid=BSPECIES.maskAll(false);
		
		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vb0=ByteVector.fromArray(BSPECIES, array, i);
			
			// Apply uppercase mask
			ByteVector vb=vb0.and(vmask);
			
			vb.intoArray(array, i);
			
			// Check if in [A, Z]
			VectorMask<Byte> validLow=vb.compare(VectorOperators.GE, vA);
			VectorMask<Byte> validHigh=vb.compare(VectorOperators.LE, vZ);
			VectorMask<Byte> valid=validLow.and(validHigh);
			
			invalid=valid.not().or(invalid);
		}
		
		boolean success=!invalid.anyTrue();
		
		for(; i<length; i++){// Residual scalar loop
			array[i]=AminoAcid.toUpperCase[array[i]];
		}
		
		return success;
	}

	
	/**
	 * Checks if all bytes are letters (case-insensitive check via mask).
	 * @param array The byte array to check.
	 * @return false if any non-letter found.
	 */
	static final boolean allLetters(final byte[] array){
		if(array==null){return true;}
		
		final int length=array.length;
		final int limit=BSPECIES.loopBound(length);
		
		int i=0;
		final byte A='A', Z='Z';
		final byte mask=~32;
		
		ByteVector vmask=ByteVector.broadcast(BSPECIES, mask);
		ByteVector vA=ByteVector.broadcast(BSPECIES, A);
		ByteVector vZ=ByteVector.broadcast(BSPECIES, Z);
		
		VectorMask<Byte> invalid=BSPECIES.maskAll(false);
		
		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vb0=ByteVector.fromArray(BSPECIES, array, i);
			
			// Apply uppercase mask
			ByteVector vb=vb0.and(vmask);
			
			// Check if in [A, Z]
			VectorMask<Byte> validLow=vb.compare(VectorOperators.GE, vA);
			VectorMask<Byte> validHigh=vb.compare(VectorOperators.LE, vZ);
			VectorMask<Byte> valid=validLow.and(validHigh);
			
			invalid=valid.not().or(invalid);
		}
		
		boolean success=!invalid.anyTrue();
		
		for(; i<length; i++){// TODO: Change to lookup table
			int b=(array[i] & mask);
			success&=(b>=A && b<=Z);
		}
		
		return success;
	}

	
	/**
	 * Converts IUPAC ambiguity codes to N, preserves A/C/G/T/U (case-insensitive).
	 * @param array The byte array to modify in-place.
	 * @return Always true.
	 */
	static final boolean iupacToN(final byte[] array){
		if(array==null){return true;}
		
		final int length=array.length;
		final int limit=BSPECIES.loopBound(length);
		
		int i=0;
		final byte A='A', C='C', G='G', T='T', U='U', N='N';
		final byte mask=~32;
		
		ByteVector vmask=ByteVector.broadcast(BSPECIES, mask);
		ByteVector vA=ByteVector.broadcast(BSPECIES, A);
		ByteVector vC=ByteVector.broadcast(BSPECIES, C);
		ByteVector vG=ByteVector.broadcast(BSPECIES, G);
		ByteVector vT=ByteVector.broadcast(BSPECIES, T);
		ByteVector vU=ByteVector.broadcast(BSPECIES, U);
		ByteVector vN=ByteVector.broadcast(BSPECIES, N);
		
		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vb0=ByteVector.fromArray(BSPECIES, array, i);
			
			// Apply uppercase mask for comparison
			ByteVector vb=vb0.and(vmask);
			
			// Check if A, C, G, T, or U
			VectorMask<Byte> isA=vb.eq(vA);
			VectorMask<Byte> isC=vb.eq(vC);
			VectorMask<Byte> isG=vb.eq(vG);
			VectorMask<Byte> isT=vb.eq(vT);
			VectorMask<Byte> isU=vb.eq(vU);
			VectorMask<Byte> isValid=isA.or(isC).or(isG).or(isT).or(isU);
			
			// Keep original if valid, replace with N if not
			ByteVector vresult=vb0.blend(vN, isValid.not());
			
			vresult.intoArray(array, i);
		}
		
		for(; i<length; i++){// Residual scalar loop
			array[i]=AminoAcid.baseToACGTN[array[i]];
		}
		
		return true;
	}

	
	/**
	 * Checks if all letters and no E or L (nucleotide validation).
	 * @param array The byte array to check.
	 * @return true if valid nucleotide sequence.
	 */
	static final boolean isNucleotide(final byte[] array){
		if(array==null){return true;}
		
		final int length=array.length;
		final int limit=BSPECIES.loopBound(length);
		
		int i=0;
		final byte E='E', L='L';
		final byte A='A', Z='Z';
		final byte mask=~32;
		
		ByteVector vmask=ByteVector.broadcast(BSPECIES, mask);
		ByteVector vA=ByteVector.broadcast(BSPECIES, A);
		ByteVector vZ=ByteVector.broadcast(BSPECIES, Z);
		ByteVector vE=ByteVector.broadcast(BSPECIES, E);
		ByteVector vL=ByteVector.broadcast(BSPECIES, L);
		
		VectorMask<Byte> invalid=BSPECIES.maskAll(false);
		
		for(; i<limit; i+=BWIDTH){// SIMD loop
			ByteVector vb0=ByteVector.fromArray(BSPECIES, array, i);
			
			// Apply uppercase mask
			ByteVector vb=vb0.and(vmask);
			
			// Check if in [A, Z]
			VectorMask<Byte> validLow=vb.compare(VectorOperators.GE, vA);
			VectorMask<Byte> validHigh=vb.compare(VectorOperators.LE, vZ);
			VectorMask<Byte> isLetter=validLow.and(validHigh);
			
			// Check if E or L
			VectorMask<Byte> isE=vb.eq(vE);
			VectorMask<Byte> isL=vb.eq(vL);
			VectorMask<Byte> isProtein=isE.or(isL);
			
			// Invalid if not letter OR is protein marker
			invalid=isLetter.not().or(isProtein).or(invalid);
		}
		
		boolean success=!invalid.anyTrue();
		
		for(; i<length; i++){// Residual scalar loop
//			int b=(array[i] & mask);
//			success&=((b>=A && b<=Z) && (b!=E && b!=L));
			success&=(AminoAcid.baseToNumberExtended[array[i]]>=0);
		}
		
		return success;
	}
	
	/** SIMD version: Add delta to each qual and append to ByteBuilder */
	static void addAndAppend(byte[] quals, ByteBuilder bb, int delta) {
		final int qlen=quals.length;
		bb.ensureExtra(qlen);
		final byte[] array=bb.array;
		final int offset=bb.length;
		
		final ByteVector vDelta=ByteVector.broadcast(BSPECIES, (byte)delta);
		
		int i=0;
		// Vector loop
		for(; i<qlen-BWIDTH+1; i+=BWIDTH){
			ByteVector vq=ByteVector.fromArray(BSPECIES, quals, i);
			ByteVector vResult=vq.add(vDelta);
			vResult.intoArray(array, offset+i);
		}
		
		// Scalar tail
		for(; i<qlen; i++){
			array[offset+i]=(byte)(quals[i]+delta);
		}
		
		bb.length+=qlen;
	}

	/** SIMD version: Generate fake quals based on whether bases are defined */
	static void appendFake(byte[] bases, ByteBuilder bb, int qFake, int qUndef) {
		final int blen=bases.length;
		bb.ensureExtra(blen);
		final byte[] array=bb.array;
		final int offset=bb.length;
		
		final int limit=BSPECIES.loopBound(blen);
		
		final ByteVector vQFake=ByteVector.broadcast(BSPECIES, (byte)qFake);
		final ByteVector vQUndef=ByteVector.broadcast(BSPECIES, (byte)qUndef);
		
		// Vectors for valid bases
		final byte mask=~32; // Uppercase mask
		final ByteVector vmask=ByteVector.broadcast(BSPECIES, mask);
		final ByteVector vA=ByteVector.broadcast(BSPECIES, (byte)'A');
		final ByteVector vC=ByteVector.broadcast(BSPECIES, (byte)'C');
		final ByteVector vG=ByteVector.broadcast(BSPECIES, (byte)'G');
		final ByteVector vT=ByteVector.broadcast(BSPECIES, (byte)'T');
		final ByteVector vU=ByteVector.broadcast(BSPECIES, (byte)'U');
		
		int i=0;
		// Vector loop
		for(; i<limit; i+=BWIDTH){
			ByteVector vBases=ByteVector.fromArray(BSPECIES, bases, i);
			
			// Apply uppercase mask for comparison
			ByteVector vb=vBases.and(vmask);
			
			// Check if A, C, G, T, or U (fully defined)
			VectorMask<Byte> isA=vb.eq(vA);
			VectorMask<Byte> isC=vb.eq(vC);
			VectorMask<Byte> isG=vb.eq(vG);
			VectorMask<Byte> isT=vb.eq(vT);
			VectorMask<Byte> isU=vb.eq(vU);
			VectorMask<Byte> isDefined=isA.or(isC).or(isG).or(isT).or(isU);
			
			// Blend: if defined use qFake, else use qUndef
			ByteVector vResult=vQFake.blend(vQUndef, isDefined.not());
			vResult.intoArray(array, offset+i);
		}
		
		// Scalar tail
		for(; i<blen; i++){
			array[offset+i]=(byte)(AminoAcid.isFullyDefined(bases[i]) ? qFake : qUndef);
		}
		
		bb.length+=blen;
	}
	
	/** SIMD version: Add delta, reverse, and append to ByteBuilder */
	static void addAndAppendReversed(byte[] quals, ByteBuilder bb, int delta) {
		final int qlen=quals.length;
		bb.ensureExtra(qlen);
		final byte[] array=bb.array;
		final int offset=bb.length;
		
		final int BWIDTH=BSPECIES.length();
		final ByteVector vDelta=ByteVector.broadcast(BSPECIES, (byte)delta);
		
		int i=0;
		int outPos=offset+qlen-BWIDTH; // Start writing from the END
		
		// Vector loop
		for(; i<qlen-BWIDTH+1; i+=BWIDTH, outPos-=BWIDTH){
			ByteVector vq=ByteVector.fromArray(BSPECIES, quals, i);
			ByteVector vAdded=vq.add(vDelta);
			ByteVector vReversed=vAdded.rearrange(B_REVERSE_SHUFFLE);
			vReversed.intoArray(array, outPos);
		}
		
		// Scalar tail - i is where vector loop stopped
		for(int k=qlen-1, j=offset+(qlen-1-i); k>=i; k--, j--){
			array[j]=(byte)(quals[k]+delta);
		}

		bb.length+=qlen;
	}

}
