package test.ukmer;

import java.util.Arrays;
import java.util.Random;

import assemble.KmerWindowDepth;
import dna.AminoAcid;
import ukmer.AbstractKmerTableU;
import ukmer.HashArrayH1D;
import ukmer.HollowKmerLookup;
import ukmer.Kmer;

/** Compares endpoint-only lookups with existing packed and alternate-span queries.
 * @author Brian Bushnell, Fischl */
public final class HollowKmerLookupTest {

	/** Exercises both strands, word boundaries, missing keys and partition routing. */
	public static void main(final String[] args){
		Kmer.PACKED=true;
		final Random random=new Random(42121);
		for(int k : new int[]{65, 66, 67, 72, 96, 124, 128, 160, 192, 224, 310, 342}){
			for(int ways : new int[]{1, 7}){
				for(int trial=0; trial<25; trial++){
					exercise(k, ways, random, trial==0);
				}
			}
		}
		invalidSources();
		System.out.println("HOLLOW_LOOKUP_TEST_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fixtures           ----------------*/
	/*--------------------------------------------------------------*/

	/** Inserts through Kmer, then looks up through independently encoded endpoint words. */
	private static void exercise(final int k, final int ways, final Random random, final boolean tie){
		assert(k>64 && ways>0) : "Fixtures must satisfy the production Hollow geometry and partition contract.";
		final AbstractKmerTableU[] tables=new AbstractKmerTableU[ways];
		for(int i=0; i<ways; i++){tables[i]=new HashArrayH1D(new int[]{1009});}
		final HollowKmerLookup lookup=new HollowKmerLookup(tables, k, true);
		final byte[] bases=new byte[k+3];
		for(int i=0; i<bases.length; i++){bases[i]=(byte)AminoAcid.numberToBase[random.nextInt(4)];}
		final byte[] nominal=Arrays.copyOf(bases, k);
		if(tie){
			for(int i=0; i<32; i++){
				nominal[k-1-i]=(byte)AminoAcid.numberToBase[3-AminoAcid.baseToNumber[nominal[i]]];
			}
		}
		final Kmer word=make(nominal);
		check(lookup.count(arm(nominal, 0), arm(nominal, k-32))==0, "Unloaded keys must return zero.");
		put(tables, word, 17);
		compare(lookup, tables, word, nominal);
		final byte[] rc=nominal.clone();
		AminoAcid.reverseComplementBasesInPlace(rc);
		compare(lookup, tables, make(rc), rc);
		final byte[] interior=nominal.clone();
		interior[32]=(byte)AminoAcid.numberToBase[(AminoAcid.baseToNumber[interior[32]]+1)&3];
		compare(lookup, tables, make(interior), interior);

		// Canonical forward anchors also address a forward-only table. Reversing
		// them must NOT hit that same entry unless the endpoints are palindromic.
		final HollowKmerLookup forward=new HollowKmerLookup(tables, k, false);
		final byte[] canonical=word.toString().getBytes(java.nio.charset.StandardCharsets.US_ASCII);
		if(word.key()==word.array2()){AminoAcid.reverseComplementBasesInPlace(canonical);}
		check(forward.count(arm(canonical, 0), arm(canonical, k-32))==17, "Forward-only lookup must keep the supplied strand.");
		AminoAcid.reverseComplementBasesInPlace(canonical);
		check(forward.count(arm(canonical, 0), arm(canonical, k-32))==(tie ? 17 : 0), "Forward-only lookup must not silently canonicalize.");

		for(int base=0; base<4; base++){
			final byte[] right=Arrays.copyOfRange(nominal, 1, k+1);
			right[k-1]=(byte)AminoAcid.numberToBase[base];
			final Kmer neighbor=word.clone();
			neighbor.addRightNumeric(base);
			put(tables, neighbor, 21+base);
			compare(lookup, tables, neighbor, right);
			final byte[] left=new byte[k];
			System.arraycopy(nominal, 0, left, 1, k-1);
			left[0]=(byte)AminoAcid.numberToBase[base];
			neighbor.setFrom(word); neighbor.addLeftNumeric(base);
			put(tables, neighbor, 31+base);
			compare(lookup, tables, neighbor, left);
		}
		if(k>=67){
			for(int span : new int[]{k-1, k+1}){alternate(lookup, tables, bases, k, span);}
		}
		incomingChoices(lookup, tables, nominal);
		incomingChoices(lookup, tables, rc);
	}

	/** Rebuilds all16 complete successor sequences, then checks the endpoint-only matrix.
	 * Canonical ties may pool cells, so compare against final table values rather
	 * than assuming every insertion created a distinct key. Both orientations run. */
	private static void incomingChoices(final HollowKmerLookup lookup, final AbstractKmerTableU[] tables,
			final byte[] bases){
		assert(bases.length>64) : "The new left-anchor end and the rightmost appended base must occupy different positions.";
		final byte[] original=bases.clone();
		final byte[] next=Arrays.copyOfRange(bases, 1, bases.length+1);
		for(int left=0; left<4; left++){
			for(int right=0; right<4; right++){
				next[31]=AminoAcid.numberToBase[left];
				next[next.length-1]=AminoAcid.numberToBase[right];
				put(tables, make(next), 1+(left<<2)+right);
			}
		}
		final int[] counts=new int[16];
		Arrays.fill(counts, -1);
		lookup.fillRightCounts16(arm(bases, 0), arm(bases, bases.length-32), counts);
		for(int left=0; left<4; left++){
			for(int right=0; right<4; right++){
				next[31]=AminoAcid.numberToBase[left];
				next[next.length-1]=AminoAcid.numberToBase[right];
				final Kmer key=make(next);
				check(counts[(left<<2)|right]==tables[key.mod(tables.length)].getValue(key),
						"Joint anchor counts must match independently rebuilt full successors in sequence orientation.");
			}
		}
		check(Arrays.equals(bases, original), "Enumerating alternatives must not edit the observed context.");
		boolean rejected=false;
		try{lookup.fillRightCounts16(0, 0, new int[15]);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "An undersized successor matrix must fail before any count is written.");
	}

	/** The diagnostic's synthetic interior must not change context or one-step counts. */
	private static void alternate(final HollowKmerLookup lookup, final AbstractKmerTableU[] tables,
			final byte[] bases, final int k, final int span){
		assert(span>=66 && bases.length>=span) : "Diagnostic alternate spans retain both 33-base contexts.";
		final Kmer word=new Kmer(k, true), neighbor=new Kmer(k, true);
		KmerWindowDepth.setHollowQuery(word, bases, 0, span);
		put(tables, word, 41);
		check(lookup.count(arm(bases, 0), arm(bases, span-32))==tables[word.mod(tables.length)].getValue(word),
				"Alternate observed spacing must use the loaded K salt.");
		for(int base=0; base<4; base++){
			neighbor.setFrom(word); neighbor.addRightNumeric(base);
			put(tables, neighbor, 51+base);
			check(lookup.count(arm(bases, 1), (arm(bases, span-32)<<2)|base)==tables[neighbor.mod(tables.length)].getValue(neighbor),
					"Alternate-span right neighbor differs from diagnostic reconstruction.");
			neighbor.setFrom(word); neighbor.addLeftNumeric(base);
			put(tables, neighbor, 61+base);
			check(lookup.count(((long)base<<62)|(arm(bases, 0)>>>2), arm(bases, span-33))==tables[neighbor.mod(tables.length)].getValue(neighbor),
					"Alternate-span left neighbor differs from diagnostic reconstruction.");
		}
	}

	/** Invalid representation must fail rather than return plausible wrong counts. */
	private static void invalidSources(){
		final AbstractKmerTableU[] pair={new HashArrayH1D(new int[]{101})};
		expectReject(pair, 64);
		expectReject(null, 128);
		expectReject(new AbstractKmerTableU[0], 128);
		expectReject(new AbstractKmerTableU[1], 128);
		expectReject(new AbstractKmerTableU[]{new HashArrayH1D(new int[]{101}, false)}, 128);
	}

	/** Checks public constructor rejection without depending on assertions being enabled. */
	private static void expectReject(final AbstractKmerTableU[] tables, final int k){
		boolean rejected=false;
		try{new HollowKmerLookup(tables, k, true);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Invalid endpoint lookup source was accepted.");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	/** Stores only through the original Kmer hashing and routing implementation. */
	private static void put(final AbstractKmerTableU[] tables, final Kmer word, final int value){
		assert(word.hollow) : "All fixture entries must use the same endpoint-only identity as the lookup.";
		tables[word.mod(tables.length)].set(word, value);
	}

	/** Compares real table values, not merely two independently formatted hash strings. */
	private static void compare(final HollowKmerLookup lookup, final AbstractKmerTableU[] tables,
			final Kmer word, final byte[] bases){
		assert(bases.length==word.kbig) : "Nominal fixtures need the same observed and table span.";
		check(lookup.count(arm(bases, 0), arm(bases, bases.length-32))==tables[word.mod(tables.length)].getValue(word),
				"Endpoint count differs from packed Kmer count at K="+word.kbig);
	}

	/** Encodes one arm directly from bases, independently of packed Kmer storage. */
	private static long arm(final byte[] bases, final int start){
		assert(start>=0 && start+32<=bases.length) : "Every fixture arm must contain exactly 32 defined bases.";
		long value=0;
		for(int i=start; i<start+32; i++){value=(value<<2)|AminoAcid.baseToNumber[bases[i]];}
		return value;
	}

	/** Makes the packed implementation serve as the independent fixture producer. */
	private static Kmer make(final byte[] bases){
		assert(bases.length>64) : "Hollow fixtures require disjoint endpoint anchors.";
		final Kmer word=new Kmer(bases.length, true);
		for(byte base : bases){word.addRight(base);}
		return word;
	}

	/** A disabled -ea must not turn failed comparisons into a passing test report. */
	private static void check(final boolean pass, final String message){
		if(!pass){throw new AssertionError(message);}
		checks++;
	}

	private static int checks;
}
