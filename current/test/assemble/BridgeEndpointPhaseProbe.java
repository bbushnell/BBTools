package test.assemble;

import java.lang.reflect.Field;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import assemble.Contig;
import assemble.Tadpole;
import fileIO.ByteStreamWriter;
import shared.Timer;
import structures.ByteBuilder;

/** Isolates bridge-end eligibility loss without graph traversal or merging.
 * @author Brian Bushnell, Fischl */
public final class BridgeEndpointPhaseProbe {

	/** Requires a fresh evidence directory and the expected outcome: lost or preserved. */
	public static void main(final String[] args) throws Exception{
		check(args.length==2, "Specify an unused evidence directory and lost|preserved.");
		check(args[1].equals("lost") || args[1].equals("preserved"), "Expected outcome must be explicit.");
		final Path root=Files.createDirectory(Paths.get(args[0]));
		final boolean expectLoss=args[1].equals("lost");
		probe(root, 300, 64, 150, expectLoss);
		probe(root, 64, 31, 40, expectLoss);
		probe(root, 31, 21, 25, expectLoss);
		System.out.println("BRIDGE_ENDPOINT_PHASE_PROBE_PASS expected="+args[1]+" cases=3");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Probe             ----------------*/
	/*--------------------------------------------------------------*/

	/** Both short tips exist, but no read contains even one long-K word. */
	private static void probe(final Path root, final int high, final int low,
			final int readLength, final boolean expectLoss) throws Exception{
		assert(low<=readLength && readLength<high) :
				"The fixture must supply short-K tip evidence but no high-K words.";
		final byte[] bases=randomBases(600, 923L+high);
		final Path reads=root.resolve("reads_"+high+"_"+low+".fa");
		final ByteBuilder bb=new ByteBuilder();
		bb.append(">left\n").append(Arrays.copyOfRange(bases, 0, readLength)).nl();
		bb.append(">right\n").append(Arrays.copyOfRange(bases, bases.length-readLength, bases.length)).nl();
		final ByteStreamWriter writer=new ByteStreamWriter(reads.toString(), false, false, false);
		writer.start();
		writer.print(bb);
		check(!writer.poisonAndWait(), "Fixture write failed: "+reads);

		final Contig direct=contig(bases);
		initialize(reads, low, direct);
		check(eligible(direct), "Short-K reference must recognize both unique supported endpoints.");

		final Contig sequential=contig(bases);
		initialize(reads, high, sequential);
		final boolean afterHigh=eligible(sequential);
		initialize(reads, low, sequential);
		final boolean afterLow=eligible(sequential);
		check(Arrays.equals(bases, sequential.bases), "Ownership initialization must not modify contig sequence.");
		System.out.println("ENDPOINT_PHASE high="+high+" low="+low+" direct="+eligible(direct)+
				" afterHigh="+afterHigh+" afterLow="+afterLow);
		check(afterLow!=expectLoss, "Observed endpoint lifetime differs from the explicitly requested expectation.");

		// Permanent topology exclusions must not be reset by a later supported table.
		final Contig protectedEnd=contig(bases);
		endpoint("leftBridgeEndpoint").setBoolean(protectedEnd, false);
		endpoint("rightBridgeEndpoint").setBoolean(protectedEnd, false);
		initialize(reads, high, protectedEnd);
		initialize(reads, low, protectedEnd);
		check(!endpoint("leftBridgeEndpoint").getBoolean(protectedEnd)
				&& !endpoint("rightBridgeEndpoint").getBoolean(protectedEnd),
				"Phase-local ownership must not revive permanently excluded loop endpoints.");
	}

	/** Loads exactly one table at a time and calls the real cross-K ownership initializer. */
	private static void initialize(final Path reads, final int k, final Contig contig) throws Exception{
		assert(contig.length()>k) : "Both terminal words must fit; this is not a short-contig exclusion test.";
		final Tadpole tad=Tadpole.makeTadpole(new String[]{"in="+reads, "k="+k, "t=1",
				"mce=1", "mcs=1", "prefilter=f", "prealloc=f", "initialsize=1000",
				"minprob=0", "minprobmain=f", "hashkmers=explicit"}, true);
		try{
			final Field cross=tad.getClass().getDeclaredField("crossKGraph");
			cross.setAccessible(true);
			cross.setBoolean(tad, true);
			tad.loadKmers(new Timer());
			final ArrayList<Contig> contigs=new ArrayList<Contig>();
			contigs.add(contig);
			tad.getClass().getMethod("initializeContigs", ArrayList.class).invoke(tad, contigs);
			check(!tad.errorState, "Table loading/ownership initialization reported an error at K="+k);
		}finally{tad.clearData();}
	}

	/** Fresh eligible outer ends are the input contract of a cross-K bridge pass. */
	private static Contig contig(final byte[] bases) throws Exception{
		assert(bases.length>300) : "All fixture contigs must contain both long-K tips.";
		final Contig result=new Contig(bases.clone(), "fixture", 0);
		endpoint("leftBridgeEndpoint").setBoolean(result, true);
		endpoint("rightBridgeEndpoint").setBoolean(result, true);
		return result;
	}

	/** Reads production state without changing its visibility or adding a product flag. */
	private static boolean eligible(final Contig contig) throws Exception{
		return endpoint("leftBridgeEndpoint").getBoolean(contig)
				&& endpoint("rightBridgeEndpoint").getBoolean(contig);
	}

	/** Reflection keeps a diagnostic fixture from changing production access modifiers. */
	private static Field endpoint(final String name) throws Exception{
		final Field field=Contig.class.getDeclaredField(name);
		field.setAccessible(true);
		return field;
	}

	/** Generates deterministic, nonhomopolymeric input without an external corpus. */
	private static byte[] randomBases(final int length, final long seed){
		assert(length>=600) : "The two read fragments must not overlap in this fixture.";
		final Random random=new Random(seed);
		final byte[] alphabet={'A', 'C', 'G', 'T'}, bases=new byte[length];
		for(int i=0; i<length; i++){bases[i]=alphabet[random.nextInt(4)];}
		return bases;
	}

	/** A diagnostic expectation mismatch must fail loudly even if assertions were disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
