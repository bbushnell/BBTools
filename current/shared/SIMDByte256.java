package shared;

import java.util.Arrays;

import dna.AminoAcid;
import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorShuffle;
import jdk.incubator.vector.VectorSpecies;
import structures.ByteBuilder;
import structures.IntList;

/**
 * Holds SIMD methods using dual-width strategy (256-bit + 64-bit).
 * This approach minimizes scalar tail operations across all array sizes.
 * @author Brian Bushnell, Isla
 * @date Nov 7, 2025
 */
final class SIMDByte256{

	private static final VectorSpecies<Byte> BSPECIES256=ByteVector.SPECIES_256;
	private static final int BWIDTH256=BSPECIES256.length(); //32
	
	private static final VectorSpecies<Byte> BSPECIES64=ByteVector.SPECIES_64;
	private static final int BWIDTH64=BSPECIES64.length(); //8
	
	private static final VectorShuffle<Byte> B_REVERSE_SHUFFLE_256;
	private static final VectorShuffle<Byte> B_REVERSE_SHUFFLE_64;
	private static final boolean[] ACGTN;
	
	static {
		//256-bit reverse shuffle
		int vlen256=BSPECIES256.length();
		int[] indices256=new int[vlen256];
		for(int i=0; i<vlen256; i++){
			indices256[i]=vlen256-1-i;
		}
		B_REVERSE_SHUFFLE_256=VectorShuffle.fromArray(BSPECIES256, indices256, 0);
		
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
	}
	
	/** Returns number of matches */
	static final int countMatches(final byte[] s1, final byte[] s2, int a1, int b1, int a2, int b2){
		int i=a1, j=a2;
		int matches=0;
		
		{//256-bit loop
			for(; j<=b2-BWIDTH256+1; i+=BWIDTH256, j+=BWIDTH256){
				ByteVector v1=ByteVector.fromArray(BSPECIES256, s1, i);
				ByteVector v2=ByteVector.fromArray(BSPECIES256, s2, j);
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

	/** Returns index of symbol */
	static final int find(final byte[] a, final byte symbol, final int from, final int to){
		int pos=from;
		
		{//256-bit loop
			for(; pos<=to-BWIDTH256; pos+=BWIDTH256){
				ByteVector v=ByteVector.fromArray(BSPECIES256, a, pos);
				VectorMask<Byte> x=v.eq(symbol);
				int t=x.firstTrue();
				if(t<BWIDTH256){ return pos+t; }
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
	 * @param a A vector.
	 * @return The sum.
	 */
	static final long sum(final byte[] a, final int from, final int to){
		int i=from;
		long c=0;
		
		{//256-bit loop
			for(; i<=to-BWIDTH256+1; i+=BWIDTH256){
				ByteVector va=ByteVector.fromArray(BSPECIES256, a, i);
				c+=va.reduceLanesToLong(VectorOperators.ADD);
			}
		}
		
		{//64-bit loop
			for(; i<=to-BWIDTH64+1; i+=BWIDTH64){
				ByteVector va=ByteVector.fromArray(BSPECIES64, a, i);
				c+=va.reduceLanesToLong(VectorOperators.ADD);
			}
		}
		
		//Scalar tail
		for(; i<=to; i++){ c+=a[i]; }
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

		{//256-bit loop
			final ByteVector newlineVec256=ByteVector.broadcast(BSPECIES256, symbol);
			for(; i<=to-BWIDTH256; i+=BWIDTH256){
				ByteVector vec=ByteVector.fromArray(BSPECIES256, buffer, i);
				VectorMask<Byte> mask=vec.eq(newlineVec256);
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
	 * Find the last symbol in buffer by scanning backwards using SIMD.
	 * @param buffer Buffer to scan
	 * @param limit Scan backwards from this position (exclusive)
	 * @return Position of last newline, or -1 if none found
	 */
	static int findLastSymbol(byte[] buffer, int limit, final byte symbol){
		//Start from the last aligned chunk and work backwards
		int i=((limit-1)/BWIDTH256)*BWIDTH256; //Round down to 256-bit boundary
		
		{//256-bit loop - work backwards
			final ByteVector newlineVec256=ByteVector.broadcast(BSPECIES256, symbol);
			for(; i>=0; i-=BWIDTH256){
				ByteVector vec=ByteVector.fromArray(BSPECIES256, buffer, i);
				VectorMask<Byte> mask=vec.eq(newlineVec256);
				long bits=mask.toLong();
				
				if(bits!=0){
					int lane=63-Long.numberOfLeadingZeros(bits); //Highest set bit
					return i+lane;
				}
			}
		}
		
		{//64-bit loop for beginning of buffer
			final ByteVector newlineVec64=ByteVector.broadcast(BSPECIES64, symbol);
			for(i+=BWIDTH256-BWIDTH64; i>=0; i-=BWIDTH64){
				ByteVector vec=ByteVector.fromArray(BSPECIES64, buffer, i);
				VectorMask<Byte> mask=vec.eq(newlineVec64);
				long bits=mask.toLong();
				
				if(bits!=0){
					int lane=63-Long.numberOfLeadingZeros(bits);
					return i+lane;
				}
			}
		}
		
		//Scalar tail for beginning of buffer
		for(i+=BWIDTH64-1; i>=0; i--){
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
		
		{//256-bit loop
			final ByteVector carrotVec256=ByteVector.broadcast(BSPECIES256, (byte)'>');
			final ByteVector slashnVec256=ByteVector.broadcast(BSPECIES256, (byte)'\n');
			
			for(; i<=to-BWIDTH256; i+=BWIDTH256){
				ByteVector vec=ByteVector.fromArray(BSPECIES256, buffer, i);
				ByteVector vecPrev=ByteVector.fromArray(BSPECIES256, buffer, i-1);
				
				VectorMask<Byte> isCarrot=vec.eq(carrotVec256);
				VectorMask<Byte> prevIsSlashn=vecPrev.eq(slashnVec256);
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
		
		{//256-bit loop
			final ByteVector vdelta256=ByteVector.broadcast(BSPECIES256, delta);
			for(; i<=length-BWIDTH256; i+=BWIDTH256){
				ByteVector va=ByteVector.fromArray(BSPECIES256, array, i);
				ByteVector vresult=va.add(vdelta256);
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
	 * Adds delta to all bytes, caps at minimum value, and returns the minimum encountered.
	 * @param array The byte array to modify in-place.
	 * @param delta The value to add to each element.
	 * @param cap The minimum allowed value after addition.
	 * @return The minimum value encountered after addition (before capping).
	 */
	static final byte addAndCapMin(final byte[] array, final byte delta, final int cap){
		if(array==null){return 0;}

		final int length=array.length;
		int i=0;
		int min=127;

		{//256-bit loop
			ByteVector vdelta256=ByteVector.broadcast(BSPECIES256, delta);
			ByteVector vcap256=ByteVector.broadcast(BSPECIES256, (byte)cap);
			ByteVector vmin=ByteVector.broadcast(BSPECIES256, (byte)127);

			for(; i<=length-BWIDTH256; i+=BWIDTH256){
				ByteVector va=ByteVector.fromArray(BSPECIES256, array, i);
				ByteVector vresult=va.add(vdelta256);
				vmin=vmin.min(vresult);
				ByteVector vcapped=vresult.max(vcap256);
				vcapped.intoArray(array, i);
			}
			
			min=Math.min(min, vmin.reduceLanes(VectorOperators.MIN));
		}
		
		{//64-bit loop
			ByteVector vdelta64=ByteVector.broadcast(BSPECIES64, delta);
			ByteVector vcap64=ByteVector.broadcast(BSPECIES64, (byte)cap);
			ByteVector vmin=ByteVector.broadcast(BSPECIES64, (byte)127);
			
			for(; i<=length-BWIDTH64; i+=BWIDTH64){
				ByteVector va=ByteVector.fromArray(BSPECIES64, array, i);
				ByteVector vresult=va.add(vdelta64);
				vmin=vmin.min(vresult);
				ByteVector vcapped=vresult.max(vcap64);
				vcapped.intoArray(array, i);
			}
			
			min=Math.min(min, vmin.reduceLanes(VectorOperators.MIN));
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
		
		{//256-bit loop
			ByteVector vdelta256=ByteVector.broadcast(BSPECIES256, (byte)delta);
			ByteVector vn256=ByteVector.broadcast(BSPECIES256, (byte)'N');
			ByteVector vzero256=ByteVector.broadcast(BSPECIES256, (byte)0);
			ByteVector vcap2_256=ByteVector.broadcast(BSPECIES256, (byte)2);

			for(; i<=length-BWIDTH256; i+=BWIDTH256){
				ByteVector vquals=ByteVector.fromArray(BSPECIES256, quals, i);
				ByteVector vbases=ByteVector.fromArray(BSPECIES256, bases, i);
				ByteVector vresult=vquals.add(vdelta256);
				vresult=vresult.max(vcap2_256);
				VectorMask<Byte> maskN=vbases.eq(vn256);
				vresult=vresult.blend(vzero256, maskN);
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
	 * Zeros quality for N bases, caps others at 2.
	 * @param quals Quality array to modify in-place.
	 * @param bases Sequence bases array.
	 */
	static final void capQuality(final byte[] quals, final byte[] bases){
		if(quals==null){return;}

		final int length=quals.length;
		int i=0;
		
		{//256-bit loop
			ByteVector vn256=ByteVector.broadcast(BSPECIES256, (byte)'N');
			ByteVector vzero256=ByteVector.broadcast(BSPECIES256, (byte)0);
			ByteVector vcap2_256=ByteVector.broadcast(BSPECIES256, (byte)2);

			for(; i<=length-BWIDTH256; i+=BWIDTH256){
				ByteVector vquals=ByteVector.fromArray(BSPECIES256, quals, i);
				ByteVector vbases=ByteVector.fromArray(BSPECIES256, bases, i);
				ByteVector vresult=vquals.max(vcap2_256);
				VectorMask<Byte> maskN=vbases.eq(vn256);
				vresult=vresult.blend(vzero256, maskN);
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
				ByteVector vresult=vquals.max(vcap2_64);
				VectorMask<Byte> maskN=vbases.eq(vn64);
				vresult=vresult.blend(vzero64, maskN);
				vresult.intoArray(quals, i);
			}
		}
		
		//Scalar tail
		for(; i<length; i++){
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
		int i=0;
		
		{//256-bit loop
			ByteVector vU256=ByteVector.broadcast(BSPECIES256, (byte)'U');
			ByteVector vu256=ByteVector.broadcast(BSPECIES256, (byte)'u');
			ByteVector vT256=ByteVector.broadcast(BSPECIES256, (byte)'T');
			ByteVector vt256=ByteVector.broadcast(BSPECIES256, (byte)'t');

			for(; i<=length-BWIDTH256; i+=BWIDTH256){
				ByteVector vbases=ByteVector.fromArray(BSPECIES256, bases, i);
				VectorMask<Byte> maskU=vbases.eq(vU256);
				VectorMask<Byte> masku=vbases.eq(vu256);
				vbases=vbases.blend(vT256, maskU);
				vbases=vbases.blend(vt256, masku);
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

	
	/**
	 * Converts lowercase letters to N.
	 * @param array The byte array to modify in-place.
	 * @return Always true.
	 */
	static final boolean lowerCaseToN(final byte[] array){
		if(array==null){return true;}

		final int length=array.length;
		int i=0;
		final byte a='a', N='N';

		{//256-bit loop
			ByteVector va256=ByteVector.broadcast(BSPECIES256, a);
			ByteVector vN256=ByteVector.broadcast(BSPECIES256, N);

			for(; i<=length-BWIDTH256; i+=BWIDTH256){
				ByteVector vb=ByteVector.fromArray(BSPECIES256, array, i);
				VectorMask<Byte> maskLower=vb.compare(VectorOperators.GE, va256);
				ByteVector vresult=vb.blend(vN256, maskLower);
				vresult.intoArray(array, i);
			}
		}
		
		{//64-bit loop
			ByteVector va64=ByteVector.broadcast(BSPECIES64, a);
			ByteVector vN64=ByteVector.broadcast(BSPECIES64, N);

			for(; i<=length-BWIDTH64; i+=BWIDTH64){
				ByteVector vb=ByteVector.fromArray(BSPECIES64, array, i);
				VectorMask<Byte> maskLower=vb.compare(VectorOperators.GE, va64);
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

	
	/**
	 * Converts dot, dash, and X to N.
	 * @param array The byte array to modify in-place.
	 * @return Always true.
	 */
	static final boolean dotDashXToN(final byte[] array){
		if(array==null){return true;}

		final int length=array.length;
		int i=0;
		final byte dot='.', dash='-', X='X', N='N';

		{//256-bit loop
			ByteVector vdot256=ByteVector.broadcast(BSPECIES256, dot);
			ByteVector vdash256=ByteVector.broadcast(BSPECIES256, dash);
			ByteVector vX256=ByteVector.broadcast(BSPECIES256, X);
			ByteVector vN256=ByteVector.broadcast(BSPECIES256, N);

			for(; i<=length-BWIDTH256; i+=BWIDTH256){
				ByteVector vb=ByteVector.fromArray(BSPECIES256, array, i);
				VectorMask<Byte> maskDot=vb.eq(vdot256);
				VectorMask<Byte> maskDash=vb.eq(vdash256);
				VectorMask<Byte> maskX=vb.eq(vX256);
				VectorMask<Byte> maskAny=maskDot.or(maskDash).or(maskX);
				ByteVector vresult=vb.blend(vN256, maskAny);
				vresult.intoArray(array, i);
			}
		}
		
		{//64-bit loop
			ByteVector vdot64=ByteVector.broadcast(BSPECIES64, dot);
			ByteVector vdash64=ByteVector.broadcast(BSPECIES64, dash);
			ByteVector vX64=ByteVector.broadcast(BSPECIES64, X);
			ByteVector vN64=ByteVector.broadcast(BSPECIES64, N);

			for(; i<=length-BWIDTH64; i+=BWIDTH64){
				ByteVector vb=ByteVector.fromArray(BSPECIES64, array, i);
				VectorMask<Byte> maskDot=vb.eq(vdot64);
				VectorMask<Byte> maskDash=vb.eq(vdash64);
				VectorMask<Byte> maskX=vb.eq(vX64);
				VectorMask<Byte> maskAny=maskDot.or(maskDash).or(maskX);
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
	
	
	/**
	 * Checks if array contains common amino acids (E or L).
	 * @param array The byte array to check.
	 * @return true if E or L found (likely protein).
	 */
	static final boolean isProtein(final byte[] array){
		if(array==null){return false;}
		
		final int length=array.length;
		int i=0;
		final byte E='E', L='L';
		
		boolean protein=false;
		
		{//256-bit loop
			VectorMask<Byte> foundProtein=BSPECIES256.maskAll(false);
			ByteVector vE256=ByteVector.broadcast(BSPECIES256, E);
			ByteVector vL256=ByteVector.broadcast(BSPECIES256, L);
			
			for(; i<=length-BWIDTH256; i+=BWIDTH256){
				ByteVector vb=ByteVector.fromArray(BSPECIES256, array, i);
				VectorMask<Byte> isE=vb.eq(vE256);
				VectorMask<Byte> isL=vb.eq(vL256);
				foundProtein=isE.or(isL).or(foundProtein);
			}
			
			protein=foundProtein.anyTrue();
		}
		
		{//64-bit loop
			VectorMask<Byte> foundProtein=BSPECIES64.maskAll(false);
			ByteVector vE64=ByteVector.broadcast(BSPECIES64, E);
			ByteVector vL64=ByteVector.broadcast(BSPECIES64, L);
			
			for(; i<=length-BWIDTH64; i+=BWIDTH64){
				ByteVector vb=ByteVector.fromArray(BSPECIES64, array, i);
				VectorMask<Byte> isE=vb.eq(vE64);
				VectorMask<Byte> isL=vb.eq(vL64);
				foundProtein=isE.or(isL).or(foundProtein);
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
	
	
	/**
	 * Converts to uppercase using bitmask. Returns false if non-letter found.
	 * @param array The byte array to modify in-place.
	 * @return false if any byte is outside A-Z range after conversion.
	 */
	static final boolean toUpperCase(final byte[] array){
		if(array==null){return true;}
		
		final int length=array.length;
		int i=0;
		final byte A='A', Z='Z';
		final byte mask=~32;
		
		boolean success=true;
		
		{//256-bit loop
			VectorMask<Byte> invalid=BSPECIES256.maskAll(false);
			ByteVector vmask256=ByteVector.broadcast(BSPECIES256, mask);
			ByteVector vA256=ByteVector.broadcast(BSPECIES256, A);
			ByteVector vZ256=ByteVector.broadcast(BSPECIES256, Z);
			
			for(; i<=length-BWIDTH256; i+=BWIDTH256){
				ByteVector vb0=ByteVector.fromArray(BSPECIES256, array, i);
				ByteVector vb=vb0.and(vmask256);
				vb.intoArray(array, i);
				VectorMask<Byte> validLow=vb.compare(VectorOperators.GE, vA256);
				VectorMask<Byte> validHigh=vb.compare(VectorOperators.LE, vZ256);
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
				vb.intoArray(array, i);
				VectorMask<Byte> validLow=vb.compare(VectorOperators.GE, vA64);
				VectorMask<Byte> validHigh=vb.compare(VectorOperators.LE, vZ64);
				VectorMask<Byte> valid=validLow.and(validHigh);
				invalid=valid.not().or(invalid);
			}
			
			success&=!invalid.anyTrue();
		}
		
		//Scalar tail
		for(; i<length; i++){
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
		int i=0;
		final byte A='A', Z='Z';
		final byte mask=~32;
		
		boolean success=true;
		
		{//256-bit loop
			VectorMask<Byte> invalid=BSPECIES256.maskAll(false);
			ByteVector vmask256=ByteVector.broadcast(BSPECIES256, mask);
			ByteVector vA256=ByteVector.broadcast(BSPECIES256, A);
			ByteVector vZ256=ByteVector.broadcast(BSPECIES256, Z);
			
			for(; i<=length-BWIDTH256; i+=BWIDTH256){
				ByteVector vb0=ByteVector.fromArray(BSPECIES256, array, i);
				ByteVector vb=vb0.and(vmask256);
				VectorMask<Byte> validLow=vb.compare(VectorOperators.GE, vA256);
				VectorMask<Byte> validHigh=vb.compare(VectorOperators.LE, vZ256);
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
		int i=0;
		final byte A='A', C='C', G='G', T='T', U='U', N='N';
		final byte mask=~32;
		
		{//256-bit loop
			ByteVector vmask256=ByteVector.broadcast(BSPECIES256, mask);
			ByteVector vA256=ByteVector.broadcast(BSPECIES256, A);
			ByteVector vC256=ByteVector.broadcast(BSPECIES256, C);
			ByteVector vG256=ByteVector.broadcast(BSPECIES256, G);
			ByteVector vT256=ByteVector.broadcast(BSPECIES256, T);
			ByteVector vU256=ByteVector.broadcast(BSPECIES256, U);
			ByteVector vN256=ByteVector.broadcast(BSPECIES256, N);
			
			for(; i<=length-BWIDTH256; i+=BWIDTH256){
				ByteVector vb0=ByteVector.fromArray(BSPECIES256, array, i);
				ByteVector vb=vb0.and(vmask256);
				VectorMask<Byte> isA=vb.eq(vA256);
				VectorMask<Byte> isC=vb.eq(vC256);
				VectorMask<Byte> isG=vb.eq(vG256);
				VectorMask<Byte> isT=vb.eq(vT256);
				VectorMask<Byte> isU=vb.eq(vU256);
				VectorMask<Byte> isValid=isA.or(isC).or(isG).or(isT).or(isU);
				ByteVector vresult=vb0.blend(vN256, isValid.not());
				vresult.intoArray(array, i);
			}
		}
		
		{//64-bit loop
			ByteVector vmask64=ByteVector.broadcast(BSPECIES64, mask);
			ByteVector vA64=ByteVector.broadcast(BSPECIES64, A);
			ByteVector vC64=ByteVector.broadcast(BSPECIES64, C);
			ByteVector vG64=ByteVector.broadcast(BSPECIES64, G);
			ByteVector vT64=ByteVector.broadcast(BSPECIES64, T);
			ByteVector vU64=ByteVector.broadcast(BSPECIES64, U);
			ByteVector vN64=ByteVector.broadcast(BSPECIES64, N);
			
			for(; i<=length-BWIDTH64; i+=BWIDTH64){
				ByteVector vb0=ByteVector.fromArray(BSPECIES64, array, i);
				ByteVector vb=vb0.and(vmask64);
				VectorMask<Byte> isA=vb.eq(vA64);
				VectorMask<Byte> isC=vb.eq(vC64);
				VectorMask<Byte> isG=vb.eq(vG64);
				VectorMask<Byte> isT=vb.eq(vT64);
				VectorMask<Byte> isU=vb.eq(vU64);
				VectorMask<Byte> isValid=isA.or(isC).or(isG).or(isT).or(isU);
				ByteVector vresult=vb0.blend(vN64, isValid.not());
				vresult.intoArray(array, i);
			}
		}
		
		//Scalar tail
		for(; i<length; i++){
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
		int i=0;
		final byte E='E', L='L';
		final byte A='A', Z='Z';
		final byte mask=~32;
		
		boolean success=true;
		
		{//256-bit loop
			VectorMask<Byte> invalid=BSPECIES256.maskAll(false);
			ByteVector vmask256=ByteVector.broadcast(BSPECIES256, mask);
			ByteVector vA256=ByteVector.broadcast(BSPECIES256, A);
			ByteVector vZ256=ByteVector.broadcast(BSPECIES256, Z);
			ByteVector vE256=ByteVector.broadcast(BSPECIES256, E);
			ByteVector vL256=ByteVector.broadcast(BSPECIES256, L);
			
			for(; i<=length-BWIDTH256; i+=BWIDTH256){
				ByteVector vb0=ByteVector.fromArray(BSPECIES256, array, i);
				ByteVector vb=vb0.and(vmask256);
				VectorMask<Byte> validLow=vb.compare(VectorOperators.GE, vA256);
				VectorMask<Byte> validHigh=vb.compare(VectorOperators.LE, vZ256);
				VectorMask<Byte> isLetter=validLow.and(validHigh);
				VectorMask<Byte> isE=vb.eq(vE256);
				VectorMask<Byte> isL=vb.eq(vL256);
				VectorMask<Byte> isProtein=isE.or(isL);
				invalid=isLetter.not().or(isProtein).or(invalid);
			}
			
			success=!invalid.anyTrue();
		}
		
		{//64-bit loop
			VectorMask<Byte> invalid=BSPECIES64.maskAll(false);
			ByteVector vmask64=ByteVector.broadcast(BSPECIES64, mask);
			ByteVector vA64=ByteVector.broadcast(BSPECIES64, A);
			ByteVector vZ64=ByteVector.broadcast(BSPECIES64, Z);
			ByteVector vE64=ByteVector.broadcast(BSPECIES64, E);
			ByteVector vL64=ByteVector.broadcast(BSPECIES64, L);
			
			for(; i<=length-BWIDTH64; i+=BWIDTH64){
				ByteVector vb0=ByteVector.fromArray(BSPECIES64, array, i);
				ByteVector vb=vb0.and(vmask64);
				VectorMask<Byte> validLow=vb.compare(VectorOperators.GE, vA64);
				VectorMask<Byte> validHigh=vb.compare(VectorOperators.LE, vZ64);
				VectorMask<Byte> isLetter=validLow.and(validHigh);
				VectorMask<Byte> isE=vb.eq(vE64);
				VectorMask<Byte> isL=vb.eq(vL64);
				VectorMask<Byte> isProtein=isE.or(isL);
				invalid=isLetter.not().or(isProtein).or(invalid);
			}
			
			success&=!invalid.anyTrue();
		}
		
		//Scalar tail
		for(; i<length; i++){
			success&=(AminoAcid.baseToNumberExtended[array[i]]>=0);
		}
		
		return success;
	}
	
	/** Dual SIMD version: Add delta to each qual and append to ByteBuilder */
	static void addAndAppend(byte[] quals, ByteBuilder bb, int delta) {
		final int qlen=quals.length;
		bb.ensureExtra(qlen);
		final byte[] array=bb.array;
		final int offset=bb.length;
		
		int i=0;
		
		{//256-bit loop
			final ByteVector vDelta256=ByteVector.broadcast(BSPECIES256, (byte)delta);
			for(; i<=qlen-BWIDTH256; i+=BWIDTH256){
				ByteVector vq=ByteVector.fromArray(BSPECIES256, quals, i);
				ByteVector vResult=vq.add(vDelta256);
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

	/** Dual SIMD version: Generate fake quals based on whether bases are defined */
	static void appendFake(byte[] bases, ByteBuilder bb, int qFake, int qUndef) {
		final int blen=bases.length;
		bb.ensureExtra(blen);
		final byte[] array=bb.array;
		final int offset=bb.length;
		
		final byte mask=~32; //Uppercase mask
		
		int i=0;
		
		{//256-bit loop
			final ByteVector vQFake256=ByteVector.broadcast(BSPECIES256, (byte)qFake);
			final ByteVector vQUndef256=ByteVector.broadcast(BSPECIES256, (byte)qUndef);
			final ByteVector vmask256=ByteVector.broadcast(BSPECIES256, mask);
			final ByteVector vA256=ByteVector.broadcast(BSPECIES256, (byte)'A');
			final ByteVector vC256=ByteVector.broadcast(BSPECIES256, (byte)'C');
			final ByteVector vG256=ByteVector.broadcast(BSPECIES256, (byte)'G');
			final ByteVector vT256=ByteVector.broadcast(BSPECIES256, (byte)'T');
			final ByteVector vU256=ByteVector.broadcast(BSPECIES256, (byte)'U');
			
			for(; i<=blen-BWIDTH256; i+=BWIDTH256){
				ByteVector vBases=ByteVector.fromArray(BSPECIES256, bases, i);
				ByteVector vb=vBases.and(vmask256);
				VectorMask<Byte> isA=vb.eq(vA256);
				VectorMask<Byte> isC=vb.eq(vC256);
				VectorMask<Byte> isG=vb.eq(vG256);
				VectorMask<Byte> isT=vb.eq(vT256);
				VectorMask<Byte> isU=vb.eq(vU256);
				VectorMask<Byte> isDefined=isA.or(isC).or(isG).or(isT).or(isU);
				ByteVector vResult=vQFake256.blend(vQUndef256, isDefined.not());
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
	
	/** Dual SIMD version: Add delta, reverse, and append to ByteBuilder */
	static void addAndAppendReversed(byte[] quals, ByteBuilder bb, int delta){
		final int qlen=quals.length;
		bb.ensureExtra(qlen);
		final byte[] array=bb.array;
		final int bufferStart=bb.length;
		final int bufferStop=bb.length+qlen;
		
		int i=bufferStart;
		int j=qlen; // Position after last unused quality score
		final int limit256=bufferStart+(qlen/BWIDTH256)*BWIDTH256;
		final int limit64=bufferStart+(qlen/BWIDTH64)*BWIDTH64;
		
		{//256-bit loop
			final ByteVector vDelta256=ByteVector.broadcast(BSPECIES256, (byte)delta);
			
			for(; i<limit256; i+=BWIDTH256, j-=BWIDTH256){
				ByteVector vq=ByteVector.fromArray(BSPECIES256, quals, j-BWIDTH256);
				ByteVector vAdded=vq.add(vDelta256);
				ByteVector vReversed=vAdded.rearrange(B_REVERSE_SHUFFLE_256);
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
	
	/** Dual SIMD version: Reverse array in-place */
	static void reverseInPlace(final byte[] array, final int len){
		if(array==null || len<2){return;}
		
		int left=0;
		int right=len;
		
		{//256-bit loop - work from both ends toward middle
			for(; left+BWIDTH256<=len/2; left+=BWIDTH256, right-=BWIDTH256){
				//Load from both ends
				ByteVector vLeft=ByteVector.fromArray(BSPECIES256, array, left);
				ByteVector vRight=ByteVector.fromArray(BSPECIES256, array, right-BWIDTH256);
				
				//Reverse each vector
				ByteVector vLeftRev=vLeft.rearrange(B_REVERSE_SHUFFLE_256);
				ByteVector vRightRev=vRight.rearrange(B_REVERSE_SHUFFLE_256);
				
				//Write to opposite ends
				vRightRev.intoArray(array, left);
				vLeftRev.intoArray(array, right-BWIDTH256);
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

	/** Dual SIMD version: Reverse-complement array in-place (ACGTacgt, others→N) */
	static void reverseComplementInPlace(final byte[] array, final int len){
		if(array==null || len<1){return;}
		int left=0;
		int right=len;
		final byte caseMask=32; // 0x20 - difference between upper and lowercase
		
		{//256-bit loop
			final ByteVector vCaseMask256=ByteVector.broadcast(BSPECIES256, caseMask);
			final ByteVector vA256=ByteVector.broadcast(BSPECIES256, (byte)'A');
			final ByteVector vC256=ByteVector.broadcast(BSPECIES256, (byte)'C');
			final ByteVector vG256=ByteVector.broadcast(BSPECIES256, (byte)'G');
			final ByteVector vT256=ByteVector.broadcast(BSPECIES256, (byte)'T');
			final ByteVector vN256=ByteVector.broadcast(BSPECIES256, (byte)'N');
			
			for(; left+BWIDTH256<=len/2; left+=BWIDTH256, right-=BWIDTH256){
				ByteVector vLeft=ByteVector.fromArray(BSPECIES256, array, left);
				ByteVector vRight=ByteVector.fromArray(BSPECIES256, array, right-BWIDTH256);
				
				//Process left vector
				ByteVector casesLeft=vLeft.and(vCaseMask256);
				ByteVector upperLeft=vLeft.and(vCaseMask256.not());
				VectorMask<Byte> isA_L=upperLeft.eq(vA256);
				VectorMask<Byte> isC_L=upperLeft.eq(vC256);
				VectorMask<Byte> isG_L=upperLeft.eq(vG256);
				VectorMask<Byte> isT_L=upperLeft.eq(vT256);
				ByteVector compLeft=vN256.blend(vT256, isA_L).blend(vG256, isC_L).blend(vC256, isG_L).blend(vA256, isT_L);
				compLeft=compLeft.or(casesLeft);
				ByteVector vLeftRevComp=compLeft.rearrange(B_REVERSE_SHUFFLE_256);
				
				//Process right vector
				ByteVector casesRight=vRight.and(vCaseMask256);
				ByteVector upperRight=vRight.and(vCaseMask256.not());
				VectorMask<Byte> isA_R=upperRight.eq(vA256);
				VectorMask<Byte> isC_R=upperRight.eq(vC256);
				VectorMask<Byte> isG_R=upperRight.eq(vG256);
				VectorMask<Byte> isT_R=upperRight.eq(vT256);
				ByteVector compRight=vN256.blend(vT256, isA_R).blend(vG256, isC_R).blend(vC256, isG_R).blend(vA256, isT_R);
				compRight=compRight.or(casesRight);
				ByteVector vRightRevComp=compRight.rearrange(B_REVERSE_SHUFFLE_256);
				
				//Write to opposite ends
				vRightRevComp.intoArray(array, left);
				vLeftRevComp.intoArray(array, right-BWIDTH256);
			}
		}
		
		{//64-bit loop
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
				
				vRightRevComp.intoArray(array, left);
				vLeftRevComp.intoArray(array, right-BWIDTH64);
			}
		}
		
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
	 * Check if array contains only uppercase A/C/G/T/N using SIMD.
	 * @param array Byte array to check
	 * @param length Number of elements to check
	 * @return true if all bases are A/C/G/T/N (uppercase only), false otherwise
	 */
	static boolean isACGTN(byte[] array, int length){
		if(array==null || length==0){return true;}
		
		int i=0;
		final byte UPPERMASK=(byte)(~32);
		
		{//256-bit loop
			VectorMask<Byte> valid=BSPECIES256.maskAll(true);
			final ByteVector vA256=ByteVector.broadcast(BSPECIES256, (byte)'A');
			final ByteVector vC256=ByteVector.broadcast(BSPECIES256, (byte)'C');
			final ByteVector vG256=ByteVector.broadcast(BSPECIES256, (byte)'G');
			final ByteVector vT256=ByteVector.broadcast(BSPECIES256, (byte)'T');
			final ByteVector vN256=ByteVector.broadcast(BSPECIES256, (byte)'N');
			final ByteVector vUPPER256=ByteVector.broadcast(BSPECIES256, UPPERMASK);
			
			for(; i<=length-BWIDTH256; i+=BWIDTH256){
				ByteVector v=ByteVector.fromArray(BSPECIES256, array, i).and(vUPPER256);
				VectorMask<Byte> isA=v.eq(vA256);
				VectorMask<Byte> isC=v.eq(vC256);
				VectorMask<Byte> isG=v.eq(vG256);
				VectorMask<Byte> isT=v.eq(vT256);
				VectorMask<Byte> isN=v.eq(vN256);
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
			if(!ACGTN[b]) {return false;}
		}
		
		return true;
	}

}