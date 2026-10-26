package test.assemble;

import java.util.Arrays;
import java.util.Random;

import assemble.KmerWindowDepth;
import assemble.Tadpole2;
import dna.AminoAcid;
import parse.PreParser;
import shared.Timer;
import stream.Read;
import structures.ByteBuilder;
import ukmer.AbstractKmerTableU;
import ukmer.Kmer;

/** Exercises the actual Tadpole2 read-extension entry points with controlled counts.
 * @author Brian Bushnell, Fischl */
public final class HollowLiveExtensionTest {

	/** Loads one small fixture per engine, then installs disjoint deterministic paths. */
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, HollowLiveExtensionTest.class, false);
		args=pp.args;
		String input=null;
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			if(split.length==2 && split[0].equalsIgnoreCase("in")){input=split[1];}
			else{throw new IllegalArgumentException("Expected in=<small fixture FASTA>, not "+arg);}
		}
		if(input==null){throw new IllegalArgumentException("A small input fixture is required.");}
		final Tadpole2 baseline=make(input, false), guarded=make(input, true);
		baseline.loadKmers(new Timer());
		guarded.loadKmers(new Timer());
		final Random random=new Random(42123);
		for(int span : new int[]{K+1, K-1}){
			for(int stop : new int[]{0, 4, 19}){
				for(int orientation=0; orientation<2; orientation++){
					final byte[] path=new byte[SEED+EXTENSION];
					for(int i=0; i<path.length; i++){path[i]=(byte)AminoAcid.numberToBase[random.nextInt(4)];}
					if(orientation==1){AminoAcid.reverseComplementBasesInPlace(path);}
					populate(baseline, path, span, stop);
					populate(guarded, path, span, stop);
					exercise(baseline, guarded, path, stop);
				}
			}
		}
		System.out.println("HOLLOW_LIVE_EXTENSION_TEST_PASS checks="+checks);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fixtures           ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructor parsing is exercised rather than changing protected fields in the test. */
	private static Tadpole2 make(final String input, final boolean guard){
		assert(input!=null) : "The ordinary table loader requires a real small fixture, not fabricated loader state.";
		return new Tadpole2(new String[]{"in="+input, "k="+K, "mode=extend", "hollow=t",
				"hollowbranchstop="+guard, "hollowbranchratio=4", "hollowbranchsupport=3",
				"t=1", "ways=7", "initialsize=1009", "minprob=0", "prefilter=f",
				"pop=f", "sweeplen=0", "ordered=t", "interleaved=f", "extendrollback=0"}, false);
	}

	/** Makes a strong unbranched nominal path and one ambiguous alternate-span choice. */
	private static void populate(final Tadpole2 engine, final byte[] path, final int span, final int stop){
		assert(stop>=0 && stop<EXTENSION) : "The veto must lie inside the otherwise fully extendable suffix.";
		final AbstractKmerTableU[] tables=engine.tables.tables();
		final Kmer word=new Kmer(K, true), neighbor=new Kmer(K, true);
		for(int end=SEED; end<=path.length; end++){
			word.clear();
			for(int i=end-K; i<end; i++){word.addRight(path[i]);}
			put(tables, word, 30);
		}
		final int end=SEED+stop;
		KmerWindowDepth.setHollowQuery(word, path, end-span, span);
		put(tables, word, 10);
		final int chosen=AminoAcid.baseToNumber[path[end]];
		for(int base=0; base<4; base++){
			neighbor.setFrom(word); neighbor.addRightNumeric(base);
			put(tables, neighbor, base==chosen || base==((chosen+1)&3) ? 5 : 1);
		}
	}

	/** Checks both junction conventions, repeated calls and Read/quality reconstruction. */
	private static void exercise(final Tadpole2 baseline, final Tadpole2 guarded, final byte[] path, final int stop){
		assert(path.length==SEED+EXTENSION) : "The expected suffix length is fixed independently of the extension implementation.";
		for(boolean junction : new boolean[]{false, true}){
			final ByteBuilder plain=prefix(path), checked=prefix(path);
			final int[] counts=new int[4];
			final int added=baseline.extendToRight2(plain, null, counts, EXTENSION, junction, new Kmer(K, true));
			check(added==EXTENSION && Arrays.equals(plain.toBytes(), path), "Guard-off nominal path changed.");
			final int limited=guarded.extendToRight2(checked, null, counts, EXTENSION, junction, new Kmer(K, true));
			check(limited==stop && Arrays.equals(checked.toBytes(), Arrays.copyOf(path, SEED+stop)),
					"Pre-append guard must stop at the diagnostic boundary, including junction-base mode.");
			final ByteBuilder repeated=prefix(path);
			int total=0;
			for(int i=0; i<EXTENSION; i++){
				final int extension=guarded.extendToRight2(repeated, null, counts, 1, junction, new Kmer(K, true));
				total+=extension;
				if(extension==0){break;}
			}
			check(total==stop && Arrays.equals(repeated.toBytes(), checked.toBytes()), "Repeated short calls bypassed the stop or leaked guard state.");
			check(guarded.extendToRight2(repeated, null, counts, 1, junction, new Kmer(K, true))==0,
					"A stopped branch must remain stopped when the caller tries again.");
		}
		final byte[] quality=new byte[SEED];
		Arrays.fill(quality, (byte)30);
		final Read read=new Read(Arrays.copyOf(path, SEED), quality, "guard_fixture", 7);
		final int added=guarded.extendRead(read, new ByteBuilder(), null, new int[4], EXTENSION, new Kmer(K, true));
		check(added==stop && Arrays.equals(read.bases, Arrays.copyOf(path, SEED+stop)), "Read entry point changed the retained prefix.");
		check(read.quality.length==read.length() && Arrays.equals(Arrays.copyOf(read.quality, SEED), quality),
				"Read extension must preserve original qualities and match the extended sequence length.");
	}

	/** Returns a builder holding only the immutable original seed. */
	private static ByteBuilder prefix(final byte[] path){
		assert(path.length>=SEED) : "Every fixture must provide the complete seed before extension.";
		final ByteBuilder builder=new ByteBuilder(SEED);
		builder.append(path, 0, SEED);
		return builder;
	}

	/** Uses ordinary packed Kmer hashing to install exact fixture counts. */
	private static void put(final AbstractKmerTableU[] tables, final Kmer word, final int count){
		assert(count>0) : "Zero is the table's empty sentinel, not a populated fixture count.";
		tables[word.mod(tables.length)].set(word, count);
	}

	/** All comparisons remain fatal even if the invocation disables assertions. */
	private static void check(final boolean pass, final String message){
		if(!pass){throw new AssertionError(message);}
		checks++;
	}

	private static int checks;
	private static final int K=128, SEED=129, EXTENSION=20;
}
