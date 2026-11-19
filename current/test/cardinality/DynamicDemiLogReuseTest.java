package test.cardinality;

import java.util.Arrays;
import java.util.Locale;

import cardinality.CardinalityTracker;
import cardinality.DynamicDemiLog;
import ddl.DDLComparison;
import dna.AminoAcid;
import kmer.AbstractKmerTable;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import simd.Vector;

/** Checks reusable native sketches against independently constructed fresh
 * sketches after a dense sequence has raised the adaptive early-exit floor.
 * Also checks literal comparison counts and caller-owned scratch reuse.
 * This is a finite correctness fixture, not a speed or memory benchmark.
 * @author Brian Bushnell, Raiden
 */
public final class DynamicDemiLogReuseTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, DynamicDemiLogReuseTest.class, false);
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
				&& parser.maxReads<0 && parser.samplerate==1, "Reuse parity uses complete synthetic sequences only");
			run();
			pp.outstream.println("DYNAMIC_DEMILOG_REUSE_TEST_PASS");
		}finally{
			DynamicDemiLog.setExponent(previous);
			Shared.closeStream(pp.outstream);
		}
	}

	private static void run(){
		DynamicDemiLog.setExponent(4);
		checkKeyEncoding();
		checkUnsignedMax();
		checkMiddleMasks();
		checkWordIterator();
		checkPrehashedPartitions();
		check(!DynamicDemiLog.blacklistExists(), "The finite fixture requires an unmasked native sketch process");
		boolean assertions=false;
		assert(assertions=true);
		check(assertions, "State parity requires native assertions enabled");
		for(final int width:new int[]{32, 256, 2048}){
			for(final int k:new int[]{5, 19, 63}){
				for(final boolean counts:new boolean[]{false, true}){
					for(final boolean kmers:new boolean[]{false, true}){checkReuse(width, k, counts, kmers);}
				}
			}
		}
		checkComparison();
		checkVectorComparison();
		checkScore();
	}

	/** Native insertion is the independent zero-bit oracle,including optional state. */
	private static void checkPrehashedPartitions(){
		final byte[] bases=dna(65536, 20261114);
		for(final boolean counts:new boolean[]{false,true}){
			for(final boolean kmers:new boolean[]{false,true}){
				final DynamicDemiLog nativeSketch=DynamicDemiLog.create(512, 9, 12345, 0, counts, kmers);
				final DynamicDemiLog prehashed=DynamicDemiLog.create(512, 9, 12345, 0, counts, kmers);
				CardinalityTracker.forEachKmer(bases, 9, 9, false, true, (word, k)->{
					nativeSketch.hashAndStore(word);
					prehashed.storeHashed(word, Tools.hash64shift(word^12345L), 0, 0);
				});
				checkEqual(nativeSketch, prehashed);
				check(nativeSketch.getGlobalNLZ()>0, "A real dense sequence exercises the adaptive floor in zero-bit native parity");
			}
		}
		final DynamicDemiLog nativeOne=DynamicDemiLog.create(512, 9, 12345, 0, true, true);
		final DynamicDemiLog partitioned=DynamicDemiLog.create(512, 9, 12345, 0, true, true);
		for(final int[] invalid:new int[][]{{0,-1},{0,10},{-1,1},{2,1},{4,2},{1,0}}){
			boolean rejected=false;
			try{partitioned.storeHashed(0, Tools.hash64shift(12345L), invalid[0], invalid[1]);}
			catch(AssertionError expected){rejected=true;}
			check(rejected && partitioned.filledBuckets()==0, "Invalid selector widths and values must fail before shifting or modifying native bins");
		}
		for(int bits=1; bits<=2; bits++){
			final int width=512>>bits;
			for(int partition=0; partition<(1<<bits); partition++){
				for(long word=0; word<1024; word++){
					nativeOne.clear();
					partitioned.clear();
					final long hash=Tools.hash64shift(word^12345L);
					nativeOne.hashAndStore(word);
					partitioned.storeHashed(word, hash, partition, bits);
					final char expected=nativeOne.maxArray()[(int)(hash&511)];
					int occupied=0;
					for(int bin=0; bin<512; bin++){
						if(partitioned.maxArray()[bin]==0){continue;}
						occupied++;
						check(bin>=partition*width && bin<(partition+1)*width
							&& partitioned.maxArray()[bin]==expected && partitioned.counts16()[bin]==1
							&& partitioned.kmerArray()[bin]==word
							&& partitioned.gcArray()[bin]==(byte)Long.bitCount(word&0x5555555555555555L),
							"Selector replacement confines one unchanged native key and its literal/count/GC payload to the requested half or quarter");
					}
					check(occupied==1 && partitioned.filledBuckets()==1, "One observed word cannot populate another partition or corrupt native occupancy");
				}
			}
		}
	}

	/** Literal hash patterns pin zero reservation, exponent wrapping and the
	 * native shift contract without constructing sketches or changing settings. */
	private static void checkKeyEncoding(){
		check(DynamicDemiLog.encodeKey(-1L, 4, 12)==1,
			"An all-one hash must retain a nonzero key when exponent and inverted mantissa are zero");
		check(DynamicDemiLog.encodeKey(0L, 4, 12)==4095,
			"An all-zero hash wraps exponent64 to zero and retains twelve inverted fraction bits");
		check(DynamicDemiLog.encodeKey(1L, 4, 12)==65535
			&& DynamicDemiLog.encodeKey(3L, 4, 12)==61439,
			"Extreme low hashes must wrap their exponents into the same four-bit fields as native insertion");
		check(DynamicDemiLog.encodeKey(3L, 6, 10)==64511,
			"Extracting the encoder must preserve native Java shift behavior even beyond the representable fraction");
		check(DynamicDemiLog.exponentBits()==4,
			"An explicit-width raw-key encoding must not change the process-global sketch recipe");
	}

	/** Preserve the existing explicit-key insertion contract independently of
	 * the removed pre-binned projection API. */
	private static void checkUnsignedMax(){
		final char[] insertion={32768,0};
		check(!DynamicDemiLog.addKey(insertion, 0, (char)32767)
			&& !DynamicDemiLog.addKey(insertion, 0, (char)32768)
			&& !DynamicDemiLog.addKey(insertion, 1, (char)0)
			&& DynamicDemiLog.addKey(insertion, 0, (char)65535)
			&& DynamicDemiLog.addKey(insertion, 1, (char)1)
			&& Arrays.equals(insertion, new char[]{65535,1}),
			"Explicit-bin insertion must raise only unsigned maxima while treating zero as empty");
	}

	/** Check each physical bit against a symbolic left-to-right word position,
	 * independently of the production block-mask arithmetic. */
	private static void checkMiddleMasks(){
		for(final boolean amino:new boolean[]{false,true}){
			final int bits=amino ? 5 : 2, limit=amino ? 12 : 31;
			for(int k=1; k<=limit; k++){
				for(int width=0; width<=k; width++){
					final long observed=AbstractKmerTable.makeMiddleMask(k, amino, width);
					final int first=(k-width+1)/2;
					long expected=0;
					for(int bit=0; bit<64; bit++){
						final int position=k-1-bit/bits;
						if(position<first || position>=first+width){expected|=1L<<bit;}
					}
					check(observed==expected, "Explicit middle masks must clear only the declared central symbols, retaining outside bits");
					if(width==2-(k&1)){
						check(observed==AbstractKmerTable.makeMiddleMask(k, amino),
							"The existing odd-one/even-two mask API must preserve its exact symbolic placement");
					}
				}
			}
		}
	}

	/** Strict spans use a literal per-position oracle; native ambiguity mode is
	 * checked against the unchanged hashSmall path, not a copied rolling loop. */
	private static void checkWordIterator(){
		final byte[][] sequences={new byte[0], dna(512, 81273),
			"acgtACGTNNACGTACGTACGTACGTACGTACGTACGTACGTACGT".getBytes(java.nio.charset.StandardCharsets.US_ASCII),
			"ACGTRYSWKMBDHVNACGTACGTACGTACGTACGTACGTACGTACGTACGT".getBytes(java.nio.charset.StandardCharsets.US_ASCII)};
		final WordBuffer observed=new WordBuffer(), expected=new WordBuffer();
		for(final int[] lengths:new int[][]{{1,1},{1,7},{7,10},{10,31},{31,31}}){
			for(final boolean canonical:new boolean[]{false,true}){
				for(final byte[] bases:sequences){
					observed.clear();
					CardinalityTracker.forEachKmer(bases, lengths[0], lengths[1], canonical, true, observed);
					expected.clear();
					literalWords(bases, lengths[0], lengths[1], canonical, expected);
					check(observed.size()==expected.size(), "Shared strict iteration must retain every complete span exactly once per declared length");
					for(int i=0; i<observed.size(); i++){
						check(observed.words.array[i]==expected.words.array[i] && observed.lengths.array[i]==expected.lengths.array[i],
							"Shared strict words must retain exact native bit codes, orientation and short-then-long emission order");
					}
					if(!canonical){continue;}
					observed.clear();
					CardinalityTracker.forEachKmer(bases, lengths[0], lengths[1], true, false, observed);
					for(final int k:lengths){
						final DynamicDemiLog nativeSketch=DynamicDemiLog.create(256, k, 12345, 0, false, false);
						final DynamicDemiLog emitted=DynamicDemiLog.create(256, k, 12345, 0, false, false);
						nativeSketch.hash(bases, null);
						for(int i=0; i<observed.size(); i++){
							if(observed.lengths.array[i]==k){emitted.add(observed.words.array[i]);}
						}
						check(Arrays.equals(nativeSketch.maxArray(), emitted.maxArray()),
							"Native-mode emitted words must reproduce unchanged hashSmall maxima even around ambiguous bytes");
					}
				}
			}
		}
	}

	private static void literalWords(final byte[] bases, final int shortK, final int longK,
			final boolean canonical, final WordBuffer output){
		assert(bases!=null && output.size()==0 && shortK>0 && shortK<=longK && longK<=31)
			: "Literal iteration starts with an empty oracle and reads complete declared DNA spans";
		final int[] lengths=shortK==longK ? new int[]{shortK} : new int[]{shortK,longK};
		for(int end=0; end<bases.length; end++){
			for(final int k:lengths){
				if(end+1<k){continue;}
				long forward=0, reverse=0;
				boolean valid=true;
				for(int i=0; i<k; i++){
					final int code="ACTG".indexOf(Character.toUpperCase((char)bases[end+1-k+i]));
					if(code<0){valid=false; break;}
					forward=4*forward+code;
					reverse|=(long)(code^2)<<(2*i);
				}
				if(valid){output.accept(canonical ? Math.max(forward, reverse) : forward, k);}
			}
		}
	}

	private static final class WordBuffer implements CardinalityTracker.KmerConsumer{
		@Override
		public void accept(final long word, final int k){
			assert(k>=1 && k<=31 && word>=0 && word<(1L<<(2*k)))
				: "The iterator must emit only unsigned packed bases inside the declared word length";
			words.add(word);
			lengths.add(k);
		}
		void clear(){
			words.clear();
			lengths.clear();
		}
		int size(){return words.size;}
		final structures.LongList words=new structures.LongList();
		final structures.IntList lengths=new structures.IntList();
	}

	/** Runs hundreds of independent sequences through one sketch. Native counts
	 * and GC are coupled by the factory; both settings cross with literal keys.
	 * Long-k mode also exercises the superclass's retained thread-local Kmer. */
	private static void checkReuse(final int width, final int k, final boolean counts,
			final boolean kmers){
		final DynamicDemiLog reused=DynamicDemiLog.create(width, k, 12345, 0, counts, kmers);
		final char[] maxima=reused.maxArray(), frequencies=reused.counts16();
		final byte[] gc=reused.gcArray();
		final long[] keys=reused.kmerArray();
		reused.id=17;
		reused.name="retained";
		for(int sequence=0; sequence<320; sequence++){
			try{
				reused.clear();
				final DynamicDemiLog fresh=DynamicDemiLog.create(width, k, 12345, 0, counts, kmers);
				check(reused.maxArray()==maxima && reused.counts16()==frequencies && reused.gcArray()==gc && reused.kmerArray()==keys,
					"Clearing must preserve every borrowed array view without allocating replacements");
				check(reused.id==17 && reused.name.equals("retained"), "Clearing observations must retain sketch metadata");
				check(reused.getLastCardinality()==fresh.getLastCardinality(), "Clear must invalidate a previously cached estimate");
				checkEqual(reused, fresh);
				final byte[] bases=sequence(k, sequence);
				//Direct bytes keep the lowercase case independent of Read normalization.
				reused.hash(bases, null);
				fresh.hash(bases, null);
				checkEqual(reused, fresh);
				if(sequence%16==0){
					check(reused.filledBuckets()>0, "The long sequence must populate native buckets");
					if(width==32){check(reused.getGlobalNLZ()>0, "Long then short must exercise a raised adaptive floor");}
				}
			}catch(AssertionError failure){
				throw new AssertionError("Reuse parity: buckets="+width+", k="+k+", counts/GC="+counts
					+", kmers="+kmers+", sequence="+sequence, failure);
			}
		}
	}

	/** Exposed native state and optional payloads must agree, including empty entries. */
	private static void checkEqual(final DynamicDemiLog a, final DynamicDemiLog b){
		check(Arrays.equals(a.maxArray(), b.maxArray()) && Arrays.equals(a.counts16(), b.counts16())
			&& Arrays.equals(a.gcArray(), b.gcArray()) && Arrays.equals(a.kmerArray(), b.kmerArray()),
			"Independent native sketches must agree in maxima, multiplicity, GC and literal keys");
		check(a.filledBuckets()==b.filledBuckets() && a.getGlobalNLZ()==b.getGlobalNLZ(),
			"Reused adaptive-floor and occupancy state must agree with fresh construction");
		check(Arrays.equals(a.bucketValues(), b.bucketValues()) && a.occupancy()==b.occupancy(),
			"Every exposed bucket and occupancy must match, including cleared entries");
		check(a.cardinality()==b.cardinality() && a.getLastCardinality()==b.getLastCardinality()
			&& Arrays.equals(a.rawEstimates(), b.rawEstimates()),
			"All raw estimators and the cached final estimate must match after every sequence");
		check(Double.doubleToLongBits(a.branch1Rate())==Double.doubleToLongBits(b.branch1Rate())
			&& Double.doubleToLongBits(a.branch2Rate())==Double.doubleToLongBits(b.branch2Rate()),
			"Diagnostic branch rates must not retain prior-sequence observations");
	}

	/** Alternates long high-occupancy inputs with short, empty, repetitive and
	 * interrupted inputs. Each cycle changes the deterministic sequence seed. */
	private static byte[] sequence(final int k, final int index){
		final int phase=index%16;
		final int length=phase==0 ? 65536 : phase==1 ? k : phase==2 ? k-1
			: phase==3 ? 0 : 17+((index*701)%4096);
		final byte[] bases=dna(length, 12345L+index);
		if(phase==4){Arrays.fill(bases, (byte)'A');}
		if(phase==5){Arrays.fill(bases, (byte)'N');}
		if(phase==6){
			for(int i=0; i<bases.length; i+=k+3){bases[i]='N';}
		}
		if(phase==7){
			for(int i=0; i<bases.length; i++){bases[i]=(byte)Character.toLowerCase(bases[i]);}
		}
		if(phase==8){AminoAcid.reverseComplementBasesInPlace(bases);}
		return bases;
	}

	private static void checkComparison(){
		final char[] left={0, 5, 8, 7}, right={0, 7, 8, 3};
		final int[] scratch={-1, -1, -1, -1, 123};
		check(DynamicDemiLog.compareDetailed(left, right, scratch)==scratch,
			"The reusable comparison overload must return the caller's same scratch array");
		check(Arrays.equals(scratch, new int[]{1, 1, 1, 1, 123}),
			"Literal buckets contain one lower, one equal, one higher and one both-empty pair");
		check(Arrays.equals(DynamicDemiLog.compareDetailed(left, right), new int[]{1, 1, 1, 1}),
			"The existing allocating overload must retain its four-field result contract");
		DynamicDemiLog.compareDetailed(new char[4], new char[4], scratch);
		check(Arrays.equals(scratch, new int[]{0, 0, 0, 4, 123}), "Repeated comparisons must overwrite all four counts");
	}

	/** A literal unsigned-bin oracle covers vector blocks, scalar tails, empty
	 * bins and reused output. The hardware gate remains mandatory on Java8. */
	private static void checkVectorComparison(){
		final boolean previous=Shared.SIMD;
		try{
			for(int mode=0; mode<2; mode++){
				if(mode==1 && !Vector.simd256){continue;}
				Shared.SIMD=mode==1;
				vectorPair(new char[]{0, 65535, 32768, 1, 0, 9}, new char[]{0, 32768, 65535, 0, 1, 9});
				final java.util.Random random=new java.util.Random(20261114);
				for(final int width:new int[]{0,1,15,16,17,31,32,33,85,256,1024}){
					final char[] left=new char[width], right=new char[width];
					for(int sample=0; sample<16; sample++){
						for(int i=0; i<width; i++){
							left[i]=(char)random.nextInt(65536);
							right[i]=i%3==0 ? left[i] : (char)random.nextInt(65536);
							if(i%5==0){left[i]=0;}
							if(i%7==0){right[i]=0;}
						}
						vectorPair(left, right);
					}
				}
			}
		}finally{Shared.SIMD=previous;}
	}

	private static void vectorPair(final char[] left, final char[] right){
		assert(left.length==right.length) : "The independent comparison oracle reads corresponding unsigned bins";
		final char[] beforeLeft=left.clone(), beforeRight=right.clone();
		final int[] expected=new int[4], scratch={-1,-1,-1,-1,123};
		for(int i=0; i<left.length; i++){
			final int a=left[i], b=right[i];
			if(a==0 && b==0){expected[3]++;}
			else if(a==b){expected[1]++;}
			else if(a<b){expected[0]++;}
			else{expected[2]++;}
		}
		check(Vector.compareDDL(left, right, scratch)==scratch && scratch[4]==123
			&& Arrays.equals(Arrays.copyOf(scratch, 4), expected),
			"Native scalar/SIMD counts must match literal unsigned bins and retain borrowed scratch/tail");
		check(Arrays.equals(Vector.compareDDL(left, right), expected)
			&& Arrays.equals(left, beforeLeft) && Arrays.equals(right, beforeRight),
			"Allocating and scratch comparisons must preserve both read-only sketches and all four categories");
		Vector.compareDDL(new char[left.length], new char[left.length], scratch);
		check(Arrays.equals(scratch, new int[]{0,0,0,left.length,123}),
			"An all-empty next pair must replace every old counter without touching the scratch tail");
	}

	/** Reusing the native comparison record must retain its exact float weighting
	 * after positive and zero-match pairs overwrite the same fields. */
	private static void checkScore(){
		final DDLComparison comparison=new DDLComparison();
		comparison.set(new int[]{0,16,0,0}, null, null, 10);
		check(comparison.cachedScore==4f, "Sixteen identical nonempty bins have unit similarity and square-root evidence four");
		for(final int[] counts:new int[][]{{0,0,0,16},{7,3,5,2},{0,4,0,0},{2,0,6,0}}){
			comparison.set(counts, null, null, 10);
			final float expected=DynamicDemiLog.wkid(counts[0], counts[1], counts[2])*(float)Math.sqrt(counts[1]);
			check(Float.floatToIntBits(comparison.cachedScore)==Float.floatToIntBits(expected)
				&& comparison.lower==counts[0] && comparison.equal==counts[1]
				&& comparison.higher==counts[2] && comparison.bothEmpty==counts[3],
				"The reusable full record must replace its prior counts and retain the native float product and square-root rounding");
		}
	}

	/** Fixed xorshift DNA supplies broad bucket occupancy without external files. */
	private static byte[] dna(final int length, final long seed){
		check(length>=0 && seed!=0, "Synthetic hashing needs a nonnegative length and nonzero xorshift seed");
		final byte[] bases=new byte[length];
		long state=seed;
		for(int i=0; i<bases.length; i++){
			state^=state<<13;
			state^=state>>>7;
			state^=state<<17;
			bases[i]=AminoAcid.numberToBase[(int)(state&3)];
		}
		return bases;
	}

	private static void check(final boolean valid, final String why){
		if(!valid){throw new AssertionError(why);}
	}
}
