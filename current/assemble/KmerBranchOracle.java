package assemble;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import dna.AminoAcid;
import map.LongHashMap;

/**
 * Exact source-genome census for collected, oriented junction words.
 * Every occurrence contributes both flanks, including occurrences after the
 * second hit. Counts describe physical reference positions; reverse-complement
 * hits contribute flanks in query orientation. This supplies local topology,
 * not a claim that a whole contig merge is safe.
 * @author Fischl
 */
public final class KmerBranchOracle {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Utility class; reference evidence must never be constructed from read depths. */
	private KmerBranchOracle(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Indexes query seeds, then scans each reference once. All queries have one K.
	 * Circularity is explicit per reference, never inferred from a header. A query
	 * may span a circular origin, but K must not exceed that reference's length.
	 * Boundary/N flanks remain unknown rather than evidence of a true dead end.
	 */
	public static Result[] scan(final List<byte[]> queries, final List<byte[]> references,
			final boolean[] circular, final int k){
		if(queries==null || references==null || circular==null || circular.length!=references.size() || k<2){
			throw new IllegalArgumentException("Oracle needs queries, references, one topology flag per reference and K>=2.");
		}
		final int seedLength=Math.min(k, 31);
		final long mask=(1L<<(2*seedLength))-1;
		final LongHashMap index=new LongHashMap(Math.max(16, queries.size()));
		final ArrayList<Pattern> patterns=new ArrayList<Pattern>();
		final int[] heads=new int[Math.toIntExact(2L*queries.size())];
		Arrays.fill(heads, -1);
		final Result[] results=new Result[queries.size()];
		for(int i=0; i<queries.size(); i++){
			final byte[] word=normalize(queries.get(i), k);
			final byte[] reverse=reverseComplement(word);
			final boolean palindrome=Arrays.equals(word, reverse);
			results[i]=new Result();
			addPattern(new Pattern(word, false, palindrome, results[i], seedLength), patterns, index, heads);
			if(!palindrome){addPattern(new Pattern(reverse, true, false, results[i], seedLength), patterns, index, heads);}
		}
		for(int r=0; r<references.size(); r++){
			final byte[] reference=references.get(r);
			if(reference==null || reference.length==0){throw new IllegalArgumentException("Empty reference at index "+r);}
			if(circular[r] && k>reference.length){
				throw new IllegalArgumentException("Circular reference shorter than K needs an explicit multi-wrap policy: "+r);
			}
			if(reference.length<k || patterns.isEmpty()){continue;}
			long seed=0;
			int valid=0;
			final long scanLength=reference.length+(circular[r] ? k-1L : 0L);
			for(long p=0; p<scanLength; p++){
				final int base=code(reference[(int)(p%reference.length)]);
				if(base<0){seed=0; valid=0; continue;}
				seed=((seed<<2)|base)&mask;
				if(valid<seedLength){valid++;}
				if(valid<seedLength){continue;}
				final int group=index.get(seed);
				for(int id=group<0 ? -1 : heads[group]; id>=0; id=patterns.get(id).next){
					final Pattern pattern=patterns.get(id);
					final long start=p-seedLength+1-pattern.offset;
					if(start<0 || start>=reference.length || (!circular[r] && start+k>reference.length)){continue;}
					if(matches(reference, start, pattern.word)){
						int left=flank(reference, start-1, circular[r]);
						int right=flank(reference, start+k, circular[r]);
						if(pattern.reverse){final int old=left; left=complement(right); right=complement(old);}
						pattern.result.add(left, right, pattern.palindrome);
					}
				}
			}
		}
		return results;
	}

	/** Links patterns by immutable map group ID; LongHashMap.put is insert-if-absent. */
	private static void addPattern(final Pattern pattern, final ArrayList<Pattern> patterns,
			final LongHashMap index, final int[] heads){
		assert(patterns.size()<heads.length) : "At most two orientations are indexed per query.";
		int group=index.get(pattern.seed);
		if(group<0){group=patterns.size(); index.put(pattern.seed, group);}
		pattern.next=heads[group];
		heads[group]=patterns.size();
		patterns.add(pattern);
	}

	/** Full equality, not seed/hash equality, establishes an occurrence. */
	private static boolean matches(final byte[] reference, final long start, final byte[] word){
		assert(start>=0 && start<reference.length) : "Each candidate is anchored to one physical reference start.";
		for(int i=0; i<word.length; i++){
			if(code(reference[(int)((start+i)%reference.length)])!=code(word[i])){return false;}
		}
		return true;
	}

	/** Returns unknown at a linear boundary or an undefined reference base. */
	private static int flank(final byte[] reference, final long position, final boolean circular){
		assert(reference.length>0) : "Modulo wrapping requires a nonempty, validated reference.";
		if(!circular && (position<0 || position>=reference.length)){return -1;}
		return code(reference[(int)Math.floorMod(position, (long)reference.length)]);
	}

	/** Copies a persistent query once and rejects ambiguous query bases. */
	private static byte[] normalize(final byte[] input, final int k){
		if(input==null || input.length!=k){throw new IllegalArgumentException("Oracle query length must equal K="+k);}
		final byte[] word=new byte[k];
		for(int i=0; i<k; i++){
			final int base=code(input[i]);
			if(base<0){throw new IllegalArgumentException("Oracle query contains a non-ACGT base.");}
			word[i]=AminoAcid.numberToBase[base];
		}
		return word;
	}

	/** Reverse-complements a saved query; no per-reference-position sequence allocation. */
	private static byte[] reverseComplement(final byte[] word){
		assert(word.length>0) : "Oracle patterns must contain a complete positive-length word.";
		final byte[] reverse=new byte[word.length];
		for(int i=0; i<word.length; i++){reverse[word.length-1-i]=AminoAcid.numberToBase[3-code(word[i])];}
		return reverse;
	}

	/** Unknown stays unknown under strand reversal. */
	private static int complement(final int base){return base<0 ? -1 : 3-base;}

	/** Byte-safe DNA conversion; non-ASCII/undefined reference bases break seeds. */
	private static int code(final byte base){return base<0 ? -1 : AminoAcid.baseToNumber[base];}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Complete aggregate of flanks, retaining paired context rather than only marginals. */
	public static final class Result {
		/** Adds one physical occurrence, with both interpretations of a palindromic word. */
		private void add(final int left, final int right, final boolean palindrome){
			assert(left>=-1 && left<4 && right>=-1 && right<4) : "Flank masks encode A/C/G/T or explicit unknown.";
			occurrences++;
			if(left<0 || right<0){unknownFlanks++;}
			addMasks(left, right);
			if(palindrome){addMasks(complement(right), complement(left));}
		}

		/** A 16-bit pair mask retains which left and right bases co-occur at one locus. */
		private void addMasks(final int left, final int right){
			assert(left>=-1 && right>=-1) : "Unknown flanks must not be shifted into a truth mask.";
			if(left>=0){leftMask|=1<<left;}
			if(right>=0){rightMask|=1<<right;}
			if(left>=0 && right>=0){pairMask|=1<<(4*left+right);}
		}

		/** Strict origin-blind traversal truth, not whole-contig or context-aware pair truth. */
		public boolean unambiguous(final int left, final int right){
			if(left<0 || left>3 || right<0 || right>3){throw new IllegalArgumentException("Expected flanks must be A/C/G/T codes.");}
			return occurrences>0 && unknownFlanks==0 && leftMask==(1<<left) && rightMask==(1<<right);
		}

		/** Physical hits in either strand; palindrome hits are counted once. */
		public long occurrences;
		/** Hits with at least one linear-boundary or undefined flank. */
		public long unknownFlanks;
		public int leftMask, rightMask, pairMask;
	}

	/** One exact orientation, with a diverse seed chosen once per retained query. */
	private static final class Pattern {
		/** Chooses the window with most adjacent-base changes; ties use the earliest. */
		Pattern(final byte[] word_, final boolean reverse_, final boolean palindrome_,
				final Result result_, final int seedLength){
			word=word_; reverse=reverse_; palindrome=palindrome_; result=result_;
			assert(word.length>=seedLength && seedLength>0 && seedLength<=31) : "Packed seed must fit 62 bits.";
			int changes=0;
			for(int i=1; i<seedLength; i++){if(word[i]!=word[i-1]){changes++;}}
			int best=changes, selected=0;
			for(int i=1; i+seedLength<=word.length; i++){
				if(word[i]!=word[i-1]){changes--;}
				if(word[i+seedLength-1]!=word[i+seedLength-2]){changes++;}
				if(changes>best){best=changes; selected=i;}
			}
			offset=selected;
			long packed=0;
			for(int i=offset; i<offset+seedLength; i++){packed=(packed<<2)|code(word[i]);}
			seed=packed;
		}

		final byte[] word;
		final boolean reverse, palindrome;
		final Result result;
		final int offset;
		final long seed;
		int next=-1;
	}
}
