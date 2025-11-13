package shared;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import dna.AminoAcid;
import ml.Cell;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;

/** 
 * Protects normal classes from seeing SIMD in case it doesn't compile or is absent.
 * ...in theory.
 * @author Brian Bushnell
 * @date Sep 12, 2013
 *
 */
public final class Vector {

	public static void main(String[] args) {
		int[] array=new int[19999];
		for(int i=0; i<array.length; i++) {
			array[i]=(byte)(i&127);
		}
		Timer t=new Timer();

		int loops=80000;

		long sum=0, sum2=0, max=0;
		for(int outer=0; outer<4; outer++) {
			t.start();
			sum=0; sum2=0;
			Shared.SIMD=false;
			for(int i=0; i<loops; i++) {
				sum=sum(array);
				sum2+=sum;
			}
			System.err.println(sum2);
			t.stop("Scalar: ");

			t.start();
			sum=0; sum2=0;
			Shared.SIMD=true;
			for(int i=0; i<loops; i++) {
				sum=sum(array);
				sum2+=sum;
			}
			System.err.println(sum2);
			t.stop("Vector: ");
		}

		for(int outer=0; outer<4; outer++) {
			t.start();
			max=0; sum=0; sum2=0;
			Shared.SIMD=false;
			for(int i=0; i<loops; i++) {
				max=max(array);
				sum2+=max;
			}
			System.err.println(max+", "+sum2);
			t.stop("Scalar: ");

			t.start();
			max=0; sum=0; sum2=0;
			Shared.SIMD=true;
			for(int i=0; i<loops; i++) {
				max=max(array);
				sum2+=max;
			}
			System.err.println(max+", "+sum2);
			t.stop("Vector: ");
		}
	}

	/** 
	 * Returns "c+=a[i]*b[i]" where a and b are equal-length arrays.
	 * @param a A vector to multiply.
	 * @param b A vector to multiply.
	 * @return Sum of products of vector elements.
	 */
	public static final float fma(final float[] a, final float[] b){
		assert(a.length==b.length);
		if(Shared.SIMD && a.length>=MINLEN32) {return SIMD.fma(a, b);}
		float c=0;
		for(int i=0; i<a.length; i++) {c+=a[i]*b[i];}
		return c;
	}

	/** 
	 * Returns "c+=a[i]*b[bSet[i]]".
	 * @param a A vector to multiply.
	 * @param b A vector to multiply.
	 * @param bSet Subset of B's indices.
	 * @param blockSize bSet should be in sets of consecutive indices of this length,
	 * for example, blockSize=8 would allow AVX256 vector operations.
	 * @return Sum of products of vector elements.
	 */
	public static final float fma(final float[] a, final float[] b, final int[] bSet, 
		final int blockSize, boolean allowSimd){
		assert(a.length==bSet.length);
		if(Shared.SIMD && a.length>=MINLEN32 && a.length==b.length) {return SIMD.fma(a, b);}
		if(Shared.SIMD && a.length>=MINLEN32 && allowSimd && ((blockSize&7)==0)) {//This ensures length-8 blocks
			return SIMD.fmaSparse(a, b, bSet);
		}
		float c=0;
		for(int i=0; i<a.length; i++) {c+=a[i]*b[bSet[i]];}
		return c;
	}

	public static final void feedForward(Cell[] layer, float[] valuesIn){
		//		assert(layer.length==valuesIn.length);
		if(Shared.SIMD && valuesIn.length>=MINLEN32) {
			SIMD.feedForward(layer, valuesIn);
			return;
		}

		for(int cnum=0; cnum<layer.length; cnum++) {
			Cell c=layer[cnum];
			float[] weights=c.weights;
			float sum=c.bias;
			assert(valuesIn.length==weights.length) : valuesIn.length+", "+weights.length;
			sum+=Vector.fma(valuesIn, weights);
			c.sum=sum;
			final float v=(float)c.activation(sum);
			c.setValue(v);
		}
	}

	public static final void feedForwardDense(Cell[] layer, float[] valuesIn){
		//		assert(layer.length==valuesIn.length);
		if(false && Shared.SIMD && valuesIn.length>=MINLEN32) {//Discovered anomaly here; very different results
			SIMD.feedForward(layer, valuesIn);
			return;
		}

		for(int cnum=0; cnum<layer.length; cnum++) {
			Cell c=layer[cnum];
			float[] weights=c.weights;
			float sum=c.bias;
			assert(valuesIn.length==weights.length) : valuesIn.length+", "+weights.length;
			sum+=Vector.fma(valuesIn, weights);
			c.sum=sum;
			final float v=(float)c.activation(sum);
			c.setValue(v);
		}
	}

	public static void backPropFma(Cell[] layer, float[] eOverNetNext, float[][] weightsOutLnum) {
		if(Shared.SIMD && eOverNetNext.length>=MINLEN32) {
			SIMD.backPropFma(layer, eOverNetNext, weightsOutLnum);
			return;
		}
		for(int i=0; i<layer.length; i++){
			Cell cell=layer[i];
			cell.eTotalOverOut=Vector.fma(weightsOutLnum[i], eOverNetNext);
		}
	}

	/** 
	 * Performs "a+=incr" where a and incr are equal-length arrays.
	 * @param a A vector to increment.
	 * @param incr Increment amount.
	 */
	public static final void add(final float[] a, final float[] incr){
		assert(a.length==incr.length);
		if(Shared.SIMD && a.length>=MINLEN32) {SIMD.add(a, incr); return;}
		for(int i=0; i<a.length; i++) {a[i]+=incr[i];}
	}

	public static final float absDifFloat(float[] a, float[] b){
		if(Shared.SIMD && a.length>=MINLEN32) {return SIMD.absDifFloat(a, b);}
		assert(a.length==b.length);
		float sum=0;
		for(int i=0; i<a.length; i++){
			sum+=Math.abs(a[i]-b[i]);
		}
		return (float)sum;
	}

	public static float[] compensate(long[] a, int k, int[] gcmap) {
		final float[] aSum=new float[k+1];
		final float inv=1f/(k+1);

		for(int i=0; i<a.length; i++) {
			int gc=gcmap[i];
			aSum[gc]+=a[i];
		}

		for(int i=0; i<aSum.length; i++) {
			aSum[i]=inv/Math.max(aSum[i], 1);
		}

		float[] comp=new float[a.length];
		for(int i=0; i<a.length; i++) {
			int gc=gcmap[i];
			comp[i]=a[i]*aSum[gc];
		}
		//This just needs to add to approximately 1.
		//    	assert(Tools.sum(comp)==1) : "k="+k+", "+Tools.sum(comp)+"\n"+Arrays.toString(aSum)+"\n"+Arrays.toString(gcmap)+"\n"+Arrays.toString(a)+"\n";
		return comp;
	}

	public static float[] compensate(int[] a, int k, int[] gcmap) {
		final float[] aSum=new float[k+1];
		final float inv=1f/(k+1);

		for(int i=0; i<a.length; i++) {
			int gc=gcmap[i];
			aSum[gc]+=a[i];
		}

		for(int i=0; i<aSum.length; i++) {
			aSum[i]=inv/Math.max(aSum[i], 1);
		}

		float[] comp=new float[a.length];
		for(int i=0; i<a.length; i++) {
			int gc=gcmap[i];
			comp[i]=a[i]*aSum[gc];
		}
		//This just needs to add to approximately 1.
		//    	assert(Tools.sum(comp)==1) : "k="+k+", "+Tools.sum(comp)+"\n"+Arrays.toString(aSum)+"\n"+Arrays.toString(gcmap)+"\n"+Arrays.toString(a)+"\n";
		return comp;
	}



	public static float absDifComp(long[] a, long[] b, int k, int[] gcmap) {
		//    	if(Shared.SIMD && a.length>Vector.MINLEN32) {return SIMD.absDifComp(a, b, k, gcmap);}
		float[] af=compensate(a, k, gcmap);
		float[] bf=compensate(b, k, gcmap);
		float ret=Vector.absDifFloat(af, bf);
		return Tools.mid(0, 1, (Float.isFinite(ret) && ret>0 ? ret : 0));
	}

	public static float cosineDifference(int[] a, int[] b) {
		float inva=1f/Math.max(1, sum(a));
		float invb=1f/Math.max(1, sum(b));
		float ret=1-cosineSimilarity(a, b, inva, invb);
		return (Float.isFinite(ret) && ret>0 ? ret : 0);
	}

	public static float cosineSimilarity(int[] a, int[] b, float inva, float invb) {
		assert(a.length==b.length);
		if(Shared.SIMD && a.length>=MINLEN32) {return SIMD.cosineSimilarity(a, b, inva, invb);}
		float dotProduct=0f;
		float normVec1=0f;
		float normVec2=0f;

		for (int i=0; i<a.length; i++) {
			float ai=a[i]*inva, bi=b[i]*invb;
			dotProduct+=ai*bi;
			normVec1+=ai*ai;
			normVec2+=bi*bi;
		}

		normVec1=Math.max(normVec1, 1e-15f);
		normVec2=Math.max(normVec2, 1e-15f);
		return (float)(dotProduct/(Math.sqrt(normVec1)*Math.sqrt(normVec2)));
	}

	/** 
	 * Performs "a[i]+=b[i]*mult" where a and b are equal-length arrays.
	 * @param a A vector to increment.
	 * @param b Increment amount.
	 * @param mult Increment multiplier.
	 */
	public static final void addProduct(final float[] a, final float[] b, final float mult){
		assert(a.length==b.length);
		if(Shared.SIMD && a.length>=MINLEN32) {SIMD.addProduct(a, b, mult); return;}
		for(int i=0; i<a.length; i++) {a[i]+=b[i]*mult;}
	}

	/** 
	 * Performs "a[i]+=b[bSet[i]]*mult".
	 * @param a A vector to increment.
	 * @param b Increment amount.
	 * @param bSet Subset of B's indices.
	 * @param mult Increment multiplier.
	 */
	public static final void addProduct(final float[] a, final float[] b, int[] bSet, final float mult, int blockSize){
		assert(a.length==bSet.length);
		if(Shared.SIMD && a.length>=MINLEN32 && a.length==b.length) {SIMD.addProduct(a, b, mult); return;}
		if(Shared.SIMD && a.length>=MINLEN32 && SIMD_MULT_SPARSE && ((blockSize&7)==0)) {SIMD.addProductSparse(a, b, bSet, mult); return;}
		for(int i=0; i<a.length; i++) {a[i]+=b[bSet[i]]*mult;}
	}

	public static void copy(float[] dest, float[] source) {
		//		assert(a.length==b.length);
		if(SIMDCOPY && Shared.SIMD && dest.length>=MINLEN32) {SIMD.copy(dest, source); return;}
		for(int i=0, max=Tools.min(dest.length, source.length); i<max; i++) {dest[i]=source[i];}
	}

	public static void copy(int[] dest, int[] source) {
		//		assert(a.length==b.length);
		//		if(SIMDCOPY && Shared.SIMD && dest.length>=MINLEN32) {SIMD.copy(dest, source); return;}//TODO
		for(int i=0, max=Tools.min(dest.length, source.length); i<max; i++) {dest[i]=source[i];}
	}

	/** Returns number of matches */
	public static final int countMatches(final byte[] s1, final byte[] s2, int a1, int b1, int a2, int b2){
		assert(b1-a1==b2-a2) : a1+"-"+b1+", "+a2+"-"+b2+", len="+s1.length+", "+(b1-a1)+"!="+(b2-a2);
		if(Shared.SIMD && b1-a1+1>=MINLEN8) {return SIMDByte256.countMatches(s1, s2, a1, b1, a2, b2);}
		int matches=0;
		for(int i=a1, j=a2; j<=b2; i++, j++) {
			final byte x=s1[i], y=s2[j];
			final int m=((x==y) ? 1 : 0);//Does not take into account capitalization or undefined bases
			matches+=m;
		}
		assert(matches>=0 && matches<=b1-a1+1);
		return matches;
	}

	public static final int countMismatches(final byte[] s1, final byte[] s2, int a1, int b1, int a2, int b2){
		return (b1-a1+1)-countMatches(s1, s2, a1, b1, a2, b2);
	}

	public static final int find(final byte[] a, final byte symbol, final int from, final int to){
		//		if(Shared.SIMD && to-from>=MINLEN8) {return SIMDByte.find(a, symbol, from, to);}
		int len=from;
		while(len<to && a[len]!=symbol){len++;}
		return len;
	}


	public static double sum(float[] array){//
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN32) {return SIMD.sum(array, 0, array.length-1);}
		double x=0;
		for(float y : array){x+=y;}
		return x;
	}

	public static long sum(byte[] array){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN8) {return SIMDByte256.sum(array, 0, array.length-1);}
		long x=0;
		for(byte y : array){x+=y;}
		return x;
	}

	public static long sum(char[] array){
		if(array==null){return 0;}
		//		if(Shared.SIMD && array.length>=SMINLEN) {return SIMD.sum(array, 0, array.length-1);}
		long x=0;
		for(char y : array){x+=y;}
		return x;
	}

	public static long sum(short[] array){
		if(array==null){return 0;}
		//		if(Shared.SIMD && array.length>=SMINLEN) {return SIMD.sum(array, 0, array.length-1);}
		long x=0;
		for(short y : array){x+=y;}
		return x;
	}

	public static long sum(int[] array){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN32) {return SIMD.sum(array, 0, array.length-1);}
		long x=0;
		for(int y : array){x+=y;}
		return x;
	}

	public static double sum(double[] array){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN64) {return SIMD.sum(array, 0, array.length-1);}
		double x=0;
		for(double y : array){x+=y;}
		return x;
	}

	public static long sum(long[] array){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN64) {return SIMD.sum(array, 0, array.length-1);}
		long x=0;
		for(long y : array){x+=y;}
		return x;
	}

	public static long sum(int[] array, int from, int to){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN32) {return SIMD.sum(array, 0, array.length-1);}
		long x=0;
		for(int i=from; i<=to; i++){x+=array[i];}
		return x;
	}

	public static long sum(long[] array, int from, int to){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN64) {return SIMD.sum(array, 0, array.length-1);}
		long x=0;
		for(int i=from; i<=to; i++){x+=array[i];}
		return x;
	}

	public static final int max(int[] array){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN32) {return SIMD.max(array, 0, array.length-1);}
		int max=array[0];
		for(int i=1; i<array.length; i++){
			int x=array[i];
			max=(x>max ? x : max);
		}
		return max;
	}

	public static final float max(float[] array){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN32) {return SIMD.max(array, 0, array.length-1);}
		float max=array[0];
		for(int i=1; i<array.length; i++){
			float x=array[i];
			max=(x>max ? x : max);
		}
		return max;
	}

	public static final long max(long[] array){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN32) {return SIMD.max(array, 0, array.length-1);}
		long max=array[0];
		for(int i=1; i<array.length; i++){
			long x=array[i];
			max=(x>max ? x : max);
		}
		return max;
	}

	/**
	 * Find positions of the given symbol in a byte array.
	 * @param buffer The byte array to search
	 * @param from Starting position (inclusive)
	 * @param to Ending position (exclusive)
	 * @param symbol Character to find
	 * @param positions IntList to store newline positions
	 * @return Number of symbols found, including pre-existing ones
	 */
	public static final int findSymbols(final byte[] array, 
		final int from, final int to, final byte symbol, final IntList positions){
		if(array==null){return positions.size();}
		if(Shared.SIMD && array.length>=MINLEN8) {
			return SIMDByte256.findSymbols(array, from, to, symbol, positions);
		}
		for(int i=from; i<to; i++){
			if(array[i]==symbol){positions.add(i);}
		}
		return positions.size();
	}

	/**
	 * Finds the last newline in buffer by scanning backwards.
	 * @param buffer Buffer to scan
	 * @param limit Scan backwards from this position (exclusive)
	 * @param symbol Symbol to find
	 * @return Position of last newline, or -1 if none found
	 */
	public static int findLastSymbol(final byte[] buffer, final int limit, final byte symbol){
		if(buffer==null){return -1;}
		if(Shared.SIMD && limit>=MINLEN8) {
			return SIMDByte256.findLastSymbol(buffer, limit, symbol);
		}
		for(int i=limit-1; i>=0; i--){
			if(buffer[i]==symbol){
				return i;
			}
		}
		return -1;
	}

	/**
	 * Find positions of FASTA record boundaries (\n>) in a byte array.
	 * Stores the position of the \n character before each >.
	 * @param buffer The byte array to search
	 * @param from Starting position (inclusive)
	 * @param to Ending position (exclusive)
	 * @param positions IntList to store boundary positions (\n positions)
	 * @return Number of boundaries found, including pre-existing ones
	 */
	public static final int findFastaHeaders(final byte[] array, 
		final int from, final int to, final IntList positions){
		if(array==null){return positions.size();}
		if(Shared.SIMD && array.length>=MINLEN8) {
			return SIMDByte256.findFastaHeaders(array, from, to, positions);
		}
		for(int i=from+1; i<to; i++){
			if(array[i]==carrot && array[i-1]==slashn){
				positions.add(i-1); // Store position of \n
			}
		}
		return positions.size();
	}

	/**
	 * Convert a FASTA record (header + sequence) to a Read.
	 * @param record Byte array containing ">header\nsequence" (may span multiple lines)
	 * @param readID Numeric ID for the read
	 * @param pairnum Pair number (0 or 1)
	 * @return Read object
	 */
	public static Read fastaRecordToRead(byte[] record, long readID, int pairnum, 
			ByteBuilder bb, IntList newlines){
		if(record==null || record.length==0){return null;}

		// Find end of header (first newline)
		int headerEnd=0;
		while(headerEnd<record.length && record[headerEnd]!='\n'){headerEnd++;}

		// Extract header (skip leading >)
		byte[] header=Arrays.copyOfRange(record, 1, headerEnd);

		// Extract and concatenate sequence lines
		headerEnd++;
		final int remaining=record.length-headerEnd;
		bb.clear();
		bb.ensureExtra(remaining);
		if(Shared.SIMD && remaining>=MINLEN8){
			newlines.clear();
			SIMDByte256.findSymbols(record, headerEnd, record.length, slashn, newlines);

			// Copy segments between newlines
			int start=headerEnd;
			for(int i=0; i<newlines.size(); i++){
				int nlPos=newlines.get(i);
				int len=nlPos-start;
				if(len>0){
					System.arraycopy(record, start, bb.array, bb.length, len);
					bb.length+=len;
				}
				start=nlPos+1; // Skip the newline
			}

			// Copy final segment after last newline
			int len=record.length-start;
			if(len>0){
				System.arraycopy(record, start, bb.array, bb.length, len);
				bb.length+=len;
			}
		}else{
			for(int i=headerEnd; i<record.length; i++){
				if(record[i]!='\n'){
					bb.array[bb.length++]=record[i];
				}
			}	
		}
//		assert(false) : bb.length+", "+headerEnd+", "+record.length;
		Read r=new Read(bb.toBytes(), null, new String(header, StandardCharsets.US_ASCII), readID);
		r.setPairnum(pairnum);
		return r;
	}

	public static void add(byte[] array, byte delta){
		if(array==null){return;}
		if(Shared.SIMD && array.length>=MINLEN8) {
			SIMDByte256.add(array, delta);
			return;
		}
		for(int i=0; i<array.length; i++) {array[i]+=delta;}
	}

	public static byte addAndCapMin(final byte[] array, final byte delta, final int cap){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN8) {
			return SIMDByte256.addAndCapMin(array, delta, cap);
		}
		int min=127;
		for(int i=0; i<array.length; i++) {
			int b=array[i]+delta;
			min=Math.min(min, b);
			array[i]=(byte)Math.max(cap, b);
		}
		return (byte)min;
	}

	public static void applyQualOffset(final byte[] quals, final byte[] bases, final int delta){
		if(quals==null){return;}
		if(Shared.SIMD && quals.length>=MINLEN8) {
			SIMDByte256.applyQualOffset(quals, bases, delta);
			return;
		}
		for(int i=0; i<quals.length; i++) {
			byte b=bases[i];
			int q=quals[i]+delta;
			q=(AminoAcid.baseToNumber[b]<0 ? 0 : Math.max(2, q));
			quals[i]=(byte)q;
		}
	}

	public static void capQuality(final byte[] quals, final byte[] bases){
		if(quals==null){return;}
		if(Shared.SIMD && quals.length>=MINLEN8) {
			SIMDByte256.capQuality(quals, bases);
			return;
		}
		for(int i=0; i<quals.length; i++) {
			byte b=bases[i];
			int q=quals[i];
			q=(AminoAcid.baseToNumber[b]<0 ? 0 : Math.max(2, q));
			quals[i]=(byte)q;
		}
	}

	public static void uToT(byte[] bases){
		if(bases==null){return;}
		if(Shared.SIMD && bases.length>=MINLEN8) {
			SIMDByte256.uToT(bases);
			return;
		}
		for(int i=0; i<bases.length; i++){
			bases[i]=AminoAcid.uToT[bases[i]];
		}
	}

	/** Returns false if there was a problem */
	public static boolean dotDashXToN(byte[] array) {
		if(array==null){return true;}
		if(Shared.SIMD && array.length>=MINLEN8) {
			return SIMDByte256.dotDashXToN(array);
		}
		for(int i=0; i<array.length; i++){
			array[i]=AminoAcid.dotDashXToNocall[array[i]];
		}
		return true;
	}

	/** Looks for common amino acids that are not IUPAC bases */
	public static boolean isProtein(byte[] array) {
		if(array==null){return true;}
		if(Shared.SIMD && array.length>=MINLEN8) {
			return SIMDByte256.isProtein(array);
		}
		boolean protein=false;
		for(int i=0; i<array.length; i++){
			byte b=array[i];
			boolean nuc=AminoAcid.baseToNumberExtended[b]>=0;
			boolean amino=AminoAcid.acidToNumberExtended[b]>=0;
			protein|=(amino && !nuc);
		}
		return protein;
	}

	/** Returns false if there was a problem */
	public static boolean lowerCaseToN(byte[] array) {
		if(array==null){return true;}
		if(Shared.SIMD && array.length>=MINLEN8) {
			return SIMDByte256.lowerCaseToN(array);
		}
		for(int i=0; i<array.length; i++){
			array[i]=AminoAcid.lowerCaseToNocall[array[i]];
		}
		return true;
	}

	/** Returns false if there was a problem */
	public static boolean toUpperCase(byte[] array) {
		if(array==null){return true;}
		if(Shared.SIMD && array.length>=MINLEN8) {
			return SIMDByte256.toUpperCase(array);
		}
		boolean success=true;
		for(int i=0; i<array.length; i++){
			array[i]=AminoAcid.toUpperCase[array[i]];
		}
		return success;
	}

	/** Returns false if there are nonletters */
	public static boolean allLetters(byte[] array) {
		if(array==null){return true;}
		if(Shared.SIMD && array.length>=MINLEN8) {
			return SIMDByte256.allLetters(array);
		}
		final byte A='A', Z='Z';
		final byte mask=~32;
		boolean success=true;
		for(int i=0; i<array.length; i++){//Could do lookup instead
			int b=(array[i] & mask);
			success&=(b>=A && b<=Z);
		}
		return success;
	}

	/** Returns false if there was a problem */
	public static boolean iupacToN(byte[] array) {
		if(array==null){return true;}
		if(Shared.SIMD && array.length>=MINLEN8) {
			return SIMDByte256.iupacToN(array);
		}
		for(int i=0; i<array.length; i++){
			array[i]=AminoAcid.baseToACGTN[array[i]];
		}
		return true;
	}

	/** Looks for common amino acids that are not IUPAC bases */
	public static boolean isNucleotide(byte[] array) {
		if(array==null){return true;}
		if(Shared.SIMD && array.length>=MINLEN8) {
			return SIMDByte256.isNucleotide(array);
		}
		boolean success=true;
		for(int i=0; i<array.length; i++){
			success&=(AminoAcid.baseToNumberExtended[array[i]]>=0);;
		}
		return success;
	}

	/** Scalar version: Add delta to each qual and append to ByteBuilder */
	public static void addAndAppend(byte[] quals, ByteBuilder bb, int delta) {
		if(quals==null){return;}
		if(Shared.SIMD && quals.length>=MINLEN8) {
			SIMDByte256.addAndAppend(quals, bb, delta);
			return;
		}
		final int qlen=quals.length;
		bb.ensureExtra(qlen);
		final byte[] array=bb.array;
		for(int i=0, j=bb.length; i<qlen; i++, j++){
			array[j]=(byte)(quals[i]+delta);
		}
		bb.length+=qlen;
	}

	/** Scalar version: Add delta to each qual and append to ByteBuilder */
	public static void addAndAppendReversed(byte[] quals, ByteBuilder bb, int delta) {
		if(quals==null){return;}
		if(Shared.SIMD && quals.length>=MINLEN8) {
			SIMDByte256.addAndAppendReversed(quals, bb, delta);
			return;
		}
		final int qlen=quals.length;
		bb.ensureExtra(qlen);
		final byte[] array=bb.array;
		for(int i=qlen-1, j=bb.length; i>=0; i--, j++){
			array[j]=(byte)(quals[i]+delta);
		}
		bb.length+=qlen;
	}

	/** Scalar version: Generate fake quals based on whether bases are defined */
	public static void appendFake(byte[] bases, ByteBuilder bb, int qFake, int qUndef) {
		if(bases==null){return;}
		if(Shared.SIMD && bases.length>=MINLEN8) {
			SIMDByte256.appendFake(bases, bb, qFake, qUndef);
			return;
		}
		final int blen=bases.length;
		bb.ensureExtra(blen);
		final byte[] array=bb.array;
		for(int i=0, j=bb.length; i<blen; i++, j++){
			array[j]=(byte)(AminoAcid.isFullyDefined(bases[i]) ? qFake : qUndef);
		}
		bb.length+=blen;
	}

	public static ByteBuilder append(ByteBuilder bb, String s) {
		if(varHandles) {
			return VarHandler.appendString(bb, s);
		} else {
			final int len=s.length();
			bb.expand(len);
			for(int i=0, j=bb.length; i<len; i++, j++){
				bb.array[j]=(byte)s.charAt(i);
			}
			bb.length+=len;
			return bb;
		}
	}

	public static void reverseInPlace(final byte[] array) {
		if(array==null){return;}
		if(Shared.SIMD && array.length>=MINLEN8) {
			SIMDByte256.reverseInPlace(array, array.length);
			return;
		}
		final int max=array.length/2, last=array.length-1;
		for(int i=0; i<max; i++){
			byte temp=array[i];
			array[i]=array[last-i];
			array[last-i]=temp;
		}
	}

	public static void reverseInPlace(final byte[] array, final int length) {
		if(array==null){return;}
		if(Shared.SIMD && array.length>=MINLEN8) {
			SIMDByte256.reverseInPlace(array, length);
			return;
		}
		final int max=length/2, last=length-1;
		for(int i=0; i<max; i++){
			byte temp=array[i];
			array[i]=array[last-i];
			array[last-i]=temp;
		}
	}

	public static void reverseComplementInPlace(final byte[] bases) {
		reverseComplementInPlace(bases, bases.length);
	}

	public static void reverseComplementInPlace(final byte[] bases, final int length) {
		if(bases==null){return;}
		if(Shared.SIMD && bases.length>=MINLEN8 && SIMDByte256.isACGTN(bases, length)) {
			SIMDByte256.reverseComplementInPlace(bases, length);
			return;
		}
		final int last=length-1;
		final int max=length/2;
		for(int i=0; i<max; i++){
			byte a=bases[i];
			byte b=bases[last-i];
			bases[i]=AminoAcid.baseToComplementExtended[b];
			bases[last-i]=AminoAcid.baseToComplementExtended[a];
		}
		if((length&1)==1){//Odd length; process middle
			bases[max]=AminoAcid.baseToComplementExtended[bases[max]];
		}
	}

	/** Converts IUPAC to N */
	public static void reverseComplementInPlaceFast(final byte[] bases) {
		reverseComplementInPlaceFast(bases, bases.length);
	}

	/** Converts IUPAC to N */
	public static void reverseComplementInPlaceFast(final byte[] bases, final int length) {
		if(bases==null){return;}
		if(Shared.SIMD && bases.length>=MINLEN8) {
			SIMDByte256.reverseComplementInPlace(bases, length);
			return;
		}
		final int last=length-1;
		final int max=length/2;
		for(int i=0; i<max; i++){
			byte a=bases[i];
			byte b=bases[last-i];
			bases[i]=AminoAcid.baseToComplementExtended[b];
			bases[last-i]=AminoAcid.baseToComplementExtended[a];
		}
		if((length&1)==1){//Odd length; process middle
			bases[max]=AminoAcid.baseToComplementExtended[bases[max]];
		}
	}

	public static int findKey(final int[] keys, final int key, final int initial, final int invalid){
		assert(keys!=null);
		if(key==invalid) {return -1;}
		if(Shared.SIMD && keys.length>=MINLEN32) {
			return SIMD.findKey(keys, key, initial, invalid);
		}
		final int limit=keys.length;
		int cell=initial;
		for(; cell<limit && keys[cell]!=key && keys[cell]!=invalid; cell++){}
		if(cell<limit) {return keys[cell]==key ? cell : -1;}
		cell=0;
		for(; cell<initial && keys[cell]!=key && keys[cell]!=invalid; cell++){}
		return keys[cell]==key ? cell : -1;
	}

	public static int findKeyOrInvalid(final int[] keys, final int key, final int initial, final int invalid){
		assert(keys!=null);
		assert(key!=invalid) : "Collision - this should have been intercepted.";
		if(Shared.SIMD && keys.length>=MINLEN32) {
			return SIMD.findKeyOrInvalid(keys, key, initial, invalid);
		}
		final int limit=keys.length;
		int cell=initial;
		for(; cell<limit && keys[cell]!=key && keys[cell]!=invalid; cell++){}
		if(cell<limit) {return cell;}
		cell=0;
		for(; cell<initial && keys[cell]!=key && keys[cell]!=invalid; cell++){}
		assert(cell<initial) : "Overflow at size "+limit+": key="+key+", initial="+initial+", invalid="+invalid;
		return cell<initial ? cell : -1;
	}

	public static int findKeyScalar(final int[] keys, final int key, final int initial, final int invalid){
		final int limit=keys.length;
		int cell=initial;
		for(; cell<limit && keys[cell]!=key && keys[cell]!=invalid; cell++){}
		if(cell<limit) {return keys[cell]==key ? cell : -1;}
		cell=0;
		for(; cell<initial && keys[cell]!=key && keys[cell]!=invalid; cell++){}
		return keys[cell]==key ? cell : -1;
	}

	public static int findKeyOrInvalidScalar(final int[] keys, final int key, final int initial, final int invalid){
		final int limit=keys.length;
		int cell=initial;
		for(; cell<limit && keys[cell]!=key && keys[cell]!=invalid; cell++){}
		if(cell<limit) {return cell;}
		cell=0;
		for(; cell<initial && keys[cell]!=key && keys[cell]!=invalid; cell++){}
		assert(cell<initial) : "Overflow at size "+limit+": key="+key+", initial="+initial+", invalid="+invalid;
		return cell<initial ? cell : -1;
	}

	public static void changeAll(int[] keys, int oldKey, int newKey){
		assert(keys!=null);
		if(Shared.SIMD && keys.length>=MINLEN32) {
			SIMD.changeAll(keys, oldKey, newKey);
			return;
		}
		for(int i=0; i<keys.length; i++){
			if(keys[i]==oldKey){keys[i]=newKey;}
		}
	}

	private static synchronized boolean vectorLoaded() {
		try{Class.forName("jdk.incubator.vector.ByteVector");}
		catch(ClassNotFoundException e){return false;}
		return true;
	}

	private static synchronized boolean simd256Avaliable() {
		try{return SIMD.maxVectorLength()>=256;}
		catch(Throwable e){return false;}
	}

	public static final int MINLEN8=8;//Due to dual SIMD
	public static final int MINLEN16=16;
	public static final int MINLEN32=16;//16 or 32 are optimal; 0, 24, and 48 are worse.
	public static final int MINLEN64=8;
	public static boolean SIMDCOPY=false;//Does not seem to affect speed, but could increase power usage.
	public static boolean SIMD_MULT_SPARSE=true;//Grants a speedup, and same results (but currently broken at ebs=1)
	public static boolean SIMD_FMA_SPARSE=true;//Grants a speedup, slightly different results
	public static final boolean vectorLoaded=vectorLoaded();
	public static final boolean simd256=simd256Avaliable();
	public static final boolean varHandles=(Shared.javaVersion>=9 && VarHandler.AVAILABLE);
	private final static byte slashr='\r', slashn='\n', carrot='>', plus='+', at='@';//, tab='\t';

}
