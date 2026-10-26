package test.ukmer;

import java.io.File;
import java.util.Arrays;
import java.util.Random;

import bloom.KmerCountAbstract;
import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import shared.Shared;
import shared.Timer;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/** Tests two-pass singleton filtering against known Solid64 occurrence counts.
 * Fixtures include absent words, reverse complements, palindromes, repeated
 * words in one read, ambiguous separators, Q0 reads and crowded Bloom lanes.
 * @author Brian Bushnell, Fischl */
public final class SingletonPrefilterTest {

	/** Generates a fresh fixture and tests native loaders, never cached counts. */
	public static void main(final String[] args){
		if(args.length!=1){throw new IllegalArgumentException("Supply a new fixture FASTQ path.");}
		final String path=args[0];
		if(new File(path).exists()){throw new IllegalArgumentException("Refusing to overwrite fixture: "+path);}
		final boolean oldPacked=Kmer.PACKED, oldBuffered=KmerCountAbstract.BUFFERED;
		final boolean oldSketch=KmerCountAbstract.SKETCH_MODE;
		final int oldThreads=Shared.threads();
		try{
			Kmer.PACKED=true;
			KmerCountAbstract.SKETCH_MODE=false;
			final Fixture fixture=new Fixture();
			fixture.write(path);
			final KmerTableSetU exact=load(path, false, 1, false, 32768);
			checkKnownCounts(exact, fixture);
			for(int threads : new int[]{1, 4}){
				for(boolean buffered : new boolean[]{false, true}){
					for(long filterBytes : new long[]{32768, 4194304}){
						final KmerTableSetU filtered=load(path, true, threads, buffered, filterBytes);
						compare(exact, filtered, fixture, threads, buffered, filterBytes);
						filtered.clear();
					}
				}
			}
			exact.clear();
			System.out.println("SINGLETON_PREFILTER_TEST_PASS checks="+checks);
		}finally{
			Kmer.PACKED=oldPacked;
			KmerCountAbstract.BUFFERED=oldBuffered;
			KmerCountAbstract.SKETCH_MODE=oldSketch;
			Shared.setThreads(oldThreads);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Native Loading        ----------------*/
	/*--------------------------------------------------------------*/

	/** Uses the proposed short-table recipe with an explicitly sized count-min filter. */
	private static KmerTableSetU load(final String path, final boolean filter,
			final int threads, final boolean buffered, final long filterBytes){
		assert(filterBytes>=25010) : "The native prefilter disables itself below100000 cells; the fixture must exercise it.";
		KmerCountAbstract.BUFFERED=buffered;
		final String[] args={"in="+path, "k=64", "hollow=f", "hashonly=f",
				"packed=t", "rcomp=t", "prefilter="+(filter ? "1" : "f"),
				"filtermem="+filterBytes, "prehashes=2", "prepasses=1", "onepass=f",
				"prealloc=f", "initialsize=1024", "ways=7", "minprob=0",
				"minprobmain=f", "qtrim=f", "maq=0", "ecco=f", "merge=f",
				"reads=-1", "interleaved=f", "qin=33", "t="+threads};
		final KmerTableSetU table=new KmerTableSetU(args, 0);
		table.process(new Timer());
		check(table.kbig==64 && !table.hollow && table.prefilter==filter,
				"Native loading must retain the requested Solid64 and filter mode.");
		return table;
	}

	/** First validates the unfiltered loader against deliberately constructed multiplicities. */
	private static void checkKnownCounts(final KmerTableSetU exact, final Fixture fixture){
		final Kmer word=new Kmer(64, false);
		for(int i=0; i<fixture.words.length; i++){
			fill(word, fixture.words[i]);
			check(count(exact, word)==fixture.expected[i],
					"Known occurrence count differs at fixture word "+i+" expected="+fixture.expected[i]);
		}
	}

	/** Tests both orientations; retained singletons prove collisions do not alter the predicate. */
	private static void compare(final KmerTableSetU exact, final KmerTableSetU filtered,
			final Fixture fixture, final int threads, final boolean buffered, final long filterBytes){
		final Kmer word=new Kmer(64, false);
		final byte[] reverse=new byte[64];
		int discarded=0, retained=0;
		for(int i=0; i<fixture.words.length; i++){
			final byte[] bases=fixture.words[i];
			for(int orientation=0; orientation<2; orientation++){
				if(orientation==0){fill(word, bases);}
				else{
					System.arraycopy(bases, 0, reverse, 0, 64);
					AminoAcid.reverseComplementBasesInPlace(reverse);
					fill(word, reverse);
				}
				final int before=count(exact, word), after=count(filtered, word);
				check(before==fixture.expected[i], "Reverse complements must preserve exact count at word "+i);
				check((before<=1)==(after<=1), "Filtering changed the weak/strong predicate at word "+i);
				check(before>=2 ? after==before : after<=before,
						"Filtering must retain exact strong counts and never create support at word "+i);
				if(orientation==0 && before==1){
					if(after==0){discarded++;}else{retained++;}
				}
			}
		}
		check(discarded>0, "The fixture must actually remove singletons, not pass with a no-op filter.");
		if(filterBytes==32768){
			check(retained>0, "The crowded fixture must retain collision-inflated singleton keys.");
		}
		System.out.println("SINGLETON_PREFILTER_CASE threads="+threads+" buffered="+buffered+
				" filterbytes="+filterBytes+" discarded_singletons="+discarded+" retained_singletons="+retained);
	}

	/** Normalizes the table's absent-key sentinel, not a stored count, to zero. */
	private static int count(final KmerTableSetU table, final Kmer word){
		return Math.max(0, table.getCount(word));
	}

	/** Rebuilds an exact complete word for independent table lookup. */
	private static void fill(final Kmer word, final byte[] bases){
		assert(bases.length==64) : "Every query must be a complete Solid64 word.";
		word.clear();
		for(byte base : bases){word.addRight(base);}
	}

	/** Fails even if the launcher accidentally disables assertions. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Fixture Data          ----------------*/
	/*--------------------------------------------------------------*/

	/** Keeps only small deterministic query arrays; streams native FASTQ through ByteStreamWriter. */
	private static final class Fixture {

		/** Constructs random count0/1/2/3 queries and four deliberate repeated-word cases. */
		Fixture(){
			final Random random=new Random(102942);
			for(int i=0; i<RANDOM_WORDS; i++){
				for(int j=0; j<64; j++){words[i][j]=AminoAcid.numberToBase[random.nextInt(4)];}
				expected[i]=i%4;
			}
			expected[0]=11;//Q0 support above the two-bit filter's saturation count must remain exact.
			Arrays.fill(words[RANDOM_WORDS], (byte)'A');
			for(int j=0; j<64; j++){
				words[RANDOM_WORDS+1][j]=(byte)((j&1)==0 ? 'A' : 'T');
				words[RANDOM_WORDS+2][j]=(byte)((j&1)==0 ? 'T' : 'A');
				words[RANDOM_WORDS+3][j]=AminoAcid.numberToBase[random.nextInt(4)];
			}
			expected[RANDOM_WORDS]=2;
			expected[RANDOM_WORDS+1]=2;
			expected[RANDOM_WORDS+2]=1;
			expected[RANDOM_WORDS+3]=2;
		}

		/** Writes known occurrences; odd copies use reverse complement and qualities include Q0. */
		void write(final String path){
			final ByteStreamWriter out=new ByteStreamWriter(path, false, false, true);
			out.start();
			final byte[] reverse=new byte[64];
			for(int i=0; i<RANDOM_WORDS; i++){
				for(int copy=0; copy<expected[i]; copy++){
					byte[] bases=words[i];
					if((copy&1)==1){
						System.arraycopy(bases, 0, reverse, 0, 64);
						AminoAcid.reverseComplementBasesInPlace(reverse);
						bases=reverse;
					}
					writeRead(out, bases, i%4);
				}
			}
			final byte[] repeated=new byte[65];
			Arrays.fill(repeated, (byte)'A');
			writeRead(out, repeated, 0);
			final byte[] palindrome=new byte[66];
			for(int i=0; i<palindrome.length; i++){palindrome[i]=(byte)((i&1)==0 ? 'A' : 'T');}
			writeRead(out, palindrome, 0);
			final byte[] separated=new byte[129];
			System.arraycopy(words[RANDOM_WORDS+3], 0, separated, 0, 64);
			separated[64]='N';
			System.arraycopy(words[RANDOM_WORDS+3], 0, separated, 65, 64);
			writeRead(out, separated, 0);
			check(!out.poisonAndWait(), "Fixture FASTQ must be completely written before native loading.");
		}

		/** Writes one short FASTQ record without per-base objects or Strings. */
		void writeRead(final ByteStreamWriter out, final byte[] bases, final int qualityIndex){
			assert(qualityIndex>=0 && qualityIndex<4) : "Fixture qualities are Q0, Q1, Q20 and Q40.";
			final int quality=(qualityIndex==0 ? 0 : qualityIndex==1 ? 1 : qualityIndex==2 ? 20 : 40)+33;
			out.print('@').print(nextId++).nl();
			out.print(bases).nl().print('+').nl();
			for(int i=0; i<bases.length; i++){out.print((byte)quality);}
			out.nl();
		}

		final byte[][] words=new byte[RANDOM_WORDS+4][64];
		final int[] expected=new int[RANDOM_WORDS+4];
		private long nextId;
	}

	private static final int RANDOM_WORDS=12000;
	private static int checks;
}
