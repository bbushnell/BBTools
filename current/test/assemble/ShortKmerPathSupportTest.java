package test.assemble;

import java.util.Arrays;
import java.util.Random;

import assemble.ShortKmerPathSupport;
import assemble.ShortKmerPathSupport.CountLookup;
import dna.AminoAcid;
import ukmer.Kmer;

/** Checks the rolling short-word ring against independently rebuilt windows.
 * No count table, reference truth or production extension is needed for these units.
 * @author Brian Bushnell, Fischl */
public final class ShortKmerPathSupportTest {

	/** Exercises geometry, both orientations, ring wrap, reset and support boundaries. */
	public static void main(final String[] args){
		Kmer.PACKED=true;
		final Random random=new Random(102342);
		final CountLookup lookup=new CountLookup(){
			@Override
			public int count(final Kmer key){return (int)(key.xor()&255);}
		};
		for(int supportK : new int[]{64, 96}){
			for(int span : new int[]{64, 65, 67, 96, 97, 124, 128, 192, 224}){
				if(span<supportK){continue;}
				final ShortKmerPathSupport support=new ShortKmerPathSupport(lookup, span, supportK);
				for(int trial=0; trial<3; trial++){
					final byte[] bases=new byte[span+150];
					for(int i=0; i<bases.length; i++){bases[i]=AminoAcid.numberToBase[random.nextInt(4)];}
					exercise(support, lookup, bases, span, supportK);
					AminoAcid.reverseComplementBasesInPlace(bases);
					exercise(support, lookup, bases, span, supportK);
					if(supportK==64){
						exercise(new ShortKmerPathSupport(lookup, span), lookup, bases, span, supportK);
					}
				}
			}
			constantSupport(0, true, supportK);
			constantSupport(1, true, supportK);
			constantSupport(2, false, supportK);
			constantSupport(Integer.MAX_VALUE, false, supportK);
			choiceOnly(supportK);
			invalidInputs(supportK);
		}
		System.out.println("SHORT_KMER_PATH_SUPPORT_TEST_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Test Methods          ----------------*/
	/*--------------------------------------------------------------*/

	/** Tests every next-base choice repeatedly before committing the saved base. */
	private static void exercise(final ShortKmerPathSupport support, final CountLookup lookup,
			final byte[] bases, final int span, final int supportK){
		assert(bases.length>span+9) : "Nonzero context offsets must fit before sequential extension starts.";
		final byte[] original=bases.clone();
		support.reset(bases, span+9);
		check(support.chosenCount()==-1, "Reset must clear previous candidate evidence.");
		for(int end=span+9; end<bases.length; end++){
			final int expectedLow=bruteLow(lookup, bases, end, span, supportK);
			check(support.lowWindowCount()==expectedLow, "Rolling weak-window count differs from rebuilding every contained word.");
			for(int candidate=3; candidate>=0; candidate--){
				final Kmer next=word(bases, end-supportK+1, end, supportK);
				next.addRightNumeric(candidate);
				final int count=lookup.count(next);
				final boolean expected=expectedLow>0 || count<=1;
				check(support.allowsCut(candidate)==expected, "Short-support decision differs from independent context and candidate evidence.");
				check(support.chosenCount()==count, "Reported candidate support differs from its exact word.");
				check(support.allowsCut(candidate)==expected, "Repeated candidate checks must not advance the context.");
				check(support.lowWindowCount()==expectedLow, "Candidate scratch must not modify the context ring.");
			}
			support.advance(AminoAcid.baseToNumber[bases[end]]);
			check(support.chosenCount()==-1, "Advancing must invalidate previous candidate evidence.");
		}
		check(Arrays.equals(bases, original), "Support checks must not edit the saved path.");
		//A discontinuous move is legal only through reset; use an earlier prefix.
		support.reset(bases, span);
		check(support.lowWindowCount()==bruteLow(lookup, bases, span, span, supportK), "Reset must replace the previous read/ring contents completely.");
	}

	/** Reconstructs every contained word independently, without a rolling count ring. */
	private static int bruteLow(final CountLookup lookup, final byte[] bases, final int end, final int span, final int supportK){
		assert(span>=supportK && end>=span && end<=bases.length) : "The brute-force context must contain complete short-word windows.";
		int low=0;
		for(int start=end-span; start<=end-supportK; start++){
			if(lookup.count(word(bases, start, start+supportK, supportK))<=1){low++;}
		}
		return low;
	}

	/** Builds a query directly from the requested half-open base range. */
	private static Kmer word(final byte[] bases, final int start, final int end, final int supportK){
		assert(start>=0 && end>=start && end<=bases.length) : "The independent word must use only the requested saved prefix.";
		final Kmer key=new Kmer(supportK, false);
		for(int i=start; i<end; i++){key.addRight(bases[i]);}
		return key;
	}

	/** Tests depth0/1 eligibility, depth2 protection, high counts and repeated words. */
	private static void constantSupport(final int count, final boolean expected, final int supportK){
		final ShortKmerPathSupport support=new ShortKmerPathSupport(new CountLookup(){
			@Override
			public int count(final Kmer key){return count;}
		}, 192, supportK);
		final byte[] bases=new byte[400];
		Arrays.fill(bases, (byte)'A');
		support.reset(bases, 192);
		for(int i=192; i<bases.length; i++){
			check(support.lowWindowCount()==(count<=1 ? 192-supportK+1 : 0), "Identical words still occupy distinct context positions.");
			check(support.allowsCut(i&3)==expected, "Fixed support thresholds or large counts changed.");
			support.advance(0);
		}
	}

	/** An absent candidate alone may qualify a cut even when every prior word is supported. */
	private static void choiceOnly(final int supportK){
		final byte[] bases=new byte[192];
		final Random random=new Random(234235);
		for(int i=0; i<bases.length; i++){bases[i]=AminoAcid.numberToBase[random.nextInt(4)];}
		final Kmer chosen=word(bases, 192-supportK+1, 192, supportK);
		chosen.addRightNumeric(2);
		final long first=chosen.xor(), second=chosen.xor2();
		final CountLookup lookup=new CountLookup(){
			@Override
			public int count(final Kmer key){return key.xor()==first && key.xor2()==second ? 0 : 20;}
		};
		final ShortKmerPathSupport support=new ShortKmerPathSupport(lookup, 192, supportK);
		support.reset(bases, 192);
		check(bruteLow(lookup, bases, 192, 192, supportK)==0 && support.lowWindowCount()==0,
				"The choice-only fixture must not contain an earlier weak context word.");
		check(support.allowsCut(2) && support.chosenCount()==0, "A missing proposed successor must qualify an existing contradiction.");
		check(!support.allowsCut(0) && support.chosenCount()==20, "A supported alternative candidate must not inherit previous weak evidence.");
	}

	/** Requires loud failures for invalid geometry, state, sequence and count contracts. */
	private static void invalidInputs(final int supportK){
		final CountLookup lookup=new CountLookup(){
			@Override
			public int count(final Kmer key){return 2;}
		};
		boolean rejected=false;
		try{new ShortKmerPathSupport(lookup, supportK-1, supportK);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Contexts shorter than the selected word cannot contain complete support.");
		for(int invalidK : new int[]{0, 63, 65, 95, 97, 128}){
			rejected=false;
			try{new ShortKmerPathSupport(lookup, 192, invalidK);}
			catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "Unsupported short geometry must fail rather than round silently.");
		}
		final ShortKmerPathSupport support=new ShortKmerPathSupport(lookup, 192, supportK);
		rejected=false;
		try{support.allowsCut(0);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Queries before reset must fail.");
		final byte[] bases=new byte[192];
		Arrays.fill(bases, (byte)'A');
		rejected=false;
		try{support.reset(bases, 191);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Truncated long contexts must fail.");
		bases[100]='N'; rejected=false;
		try{support.reset(bases, 192);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Ambiguous bases must not silently remove short windows.");
		bases[100]='A'; support.reset(bases, 192); rejected=false;
		try{support.advance(4);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Undefined appended bases must fail.");
		final ShortKmerPathSupport negative=new ShortKmerPathSupport(new CountLookup(){
			@Override
			public int count(final Kmer key){return -1;}
		}, 192, supportK);
		rejected=false;
		try{negative.reset(bases, 192);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Unnormalized lookup failures must not be treated as absent support.");
	}

	/** Unit failures stay fatal even if assertions are accidentally disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	private static int checks;
}
