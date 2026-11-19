package test.prok;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Locale;

import dna.GeneticCode;
import prok.CallGenes;
import prok.CallGenesIteration;
import prok.GeneCaller;
import prok.GeneModel;
import prok.ProkObject;
import prok.StopTailCascadeResources;
import prok.StopTailCascadeSupport;
import prok.StopTailSupport;
import prot.StopTailDetector;
import shared.Shared;
import stream.Read;

/** Real private-tail objects must not enter public ncRNA iteration scopes.
 * Reflection only invokes real legacy constructors over a literal tiny FASTA;
 * it avoids the hash-only legacy loader without changing product admission.
 * No uninitialized objects, mocks, payload hashes or biological panel are used.
 * Run in an isolated JVM: this fixture disables RNA flags and sets the shared thread
 * count without restoring them. It retains generated tiny inputs under a fresh root.
 * Documentation and formatting by Fischl; original authorship is unresolved.
 */
public final class NcrnaPrivateTailIterationTest{

	/** Checks four protein-support/iteration rejections and three ordinary seed scopes.
	 * @param args Fresh absolute fixture directory and absolute existing PGM file;
	 * the root parent must exist and the JVM must have no preloaded ncRNA families */
	public static void main(String[] args)throws Exception{
		if(args.length!=2){
			throw new IllegalArgumentException("Supply a fresh absolute fixture directory and explicit existing PGM");
		}
		final Path root=Paths.get(args[0]);
		final Path pgm=Paths.get(args[1]);
		check(pgm.isAbsolute() && Files.isRegularFile(pgm), "Explicit PGM must exist");
		check(root.isAbsolute() && !Files.exists(root), "Fresh absolute fixture root required");
		Files.createDirectory(root);
		Shared.setThreads(1);
		ProkObject.call16S=ProkObject.call23S=ProkObject.call5S=ProkObject.call18S=ProkObject.calltRNA=false;
		check(GeneCaller.ncrnaFamilies.isEmpty(), "The fresh guard fixture must not have unrelated RNA families");
		final StringBuilder peptide=new StringBuilder();
		for(int i=0; i<5; i++){peptide.append("ACDEFGHIKLMNPQRSTVWY");}
		final Path tails=root.resolve("legacy.faa");
		Files.write(tails, (">literal100\n"+peptide+"\n").getBytes(StandardCharsets.US_ASCII));
		final Path dna=root.resolve("driver.fa");
		Files.write(dna, ">intact\nATGAAATAA\n".getBytes(StandardCharsets.US_ASCII));
		final CallGenes driver=new CallGenes(new String[]{"in="+dna, "pgm="+pgm, "annotation=off", "taxonomy=f",
			"cdsonly=t", "minlen=60", "passes=1", "t=1", "overwrite=f", "append=f"});
		final CallGenesIteration factory=(input, caller, model)->{
			throw new AssertionError("Guard fixture must not execute workers");
		};
		// The real constructor has a package-private parameter type; do not widen product visibility for this fixture.
		final Class<?> pairedType=Class.forName("prot.HbmTailPairedShortlist", false,
			StopTailDetector.Reference.class.getClassLoader());
		final Constructor<StopTailDetector.Reference> detector=StopTailDetector.Reference.class.getDeclaredConstructor(
			String.class, int.class, int.class, pairedType, StopTailDetector.HbmInput.class);
		detector.setAccessible(true);
		final StopTailDetector.Reference reference=detector.newInstance(tails.toString(), 0, 0, null, null);
		final Constructor<StopTailSupport> legacy=StopTailSupport.class.getDeclaredConstructor(
			StopTailDetector.Reference.class, String.class, String.class, StopTailDetector.Method.class,
			double.class, float.class, float.class);
		legacy.setAccessible(true);
		try(StopTailSupport support=legacy.newInstance(reference, "unmeasured", "unmeasured",
			StopTailDetector.Method.AAI, 70d, 0f, 0f)){
			final GeneCaller caller=caller();
			caller.setStopTailSupport(support);
			rejectCandidateReuse(caller, "legacy_stop_tail");
			rejectDriver(driver, factory, "stopTailReference", support);
			caller.setStopTailSupport(null);
			ordinaryScope(caller);
		}
		final Path manifest=StopTailCascadeStructuralTest.writeInformationResources(root.resolve("cascade"));
		final StopTailCascadeResources resources=StopTailCascadeResources.readStructural(manifest.toString(), "unmeasured");
		try(StopTailCascadeSupport support=resources.loadSupport(1)){
			final GeneCaller caller=caller();
			caller.setTailCascadeSupportStructural(support, 0, false, "unmeasured", false);
			rejectCandidateReuse(caller, "tail_cascade");
			rejectDriver(driver, factory, "tailCascadeReference", support);
			caller.setTailCascadeSupportStructural(null, 0, false, "unmeasured", false);
			ordinaryScope(caller);
		}
		ordinaryScope(caller());
		driver.setIteration(factory);
		System.out.println("NCRNA_PRIVATE_TAIL_ITERATION_TEST_PASS rejected=4 ordinary_seed_scopes=3 real_objects=true");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Scope checks         ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructs a fixture caller and fresh gene model with translation table11 and literal thresholds. */
	private static GeneCaller caller(){
		final GeneticCode code=GeneticCode.forTable(11);
		return new GeneCaller(60, 60, 60, -999, -999, -999, -999, -999, new GeneModel(true, code), code);
	}

	/** Returns a new literal nine-base read, so orientation checks never reuse a mutated input. */
	private static Read read(){
		return new Read("ATGAAATAA".getBytes(StandardCharsets.US_ASCII), null, "intact", 0);
	}

	/** Requires candidate-reuse rejection and verifies that the borrowed read bases are unchanged. */
	private static void rejectCandidateReuse(final GeneCaller caller, final String mode){
		final Read read=read();
		final byte[] before=read.bases.clone();
		reject(()->caller.beginNcrnaIteration(read, new int[]{0}), mode+" candidate reuse");
		check(Arrays.equals(before, read.bases), "Rejected scope changed input orientation");
	}

	/** Temporarily injects real protein support into one initially-null driver field.
	 * Requires iteration-factory rejection and clears the field in finally; support is caller-owned. */
	private static void rejectDriver(final CallGenes driver, final CallGenesIteration factory,
			final String fieldName, final Object support)throws Exception{
		final Field field=CallGenes.class.getDeclaredField(fieldName);
		field.setAccessible(true);
		check(field.get(driver)==null, "Driver must start without private tail support");
		field.set(driver, support);
		try{
			reject(()->driver.setIteration(factory), fieldName+" driver iteration");
		}finally{
			field.set(driver, null);
		}
	}

	/** Requires IllegalArgumentException mentioning protein and reuse/iteration; other failures propagate. */
	private static void reject(final Runnable action, final String label){
		boolean rejected=false;
		try{action.run();}catch(IllegalArgumentException expected){
			final String message=String.valueOf(expected.getMessage()).toLowerCase(Locale.ROOT);
			check(message.contains("protein") && (message.contains("reuse") || message.contains("iteration")),
				"Wrong rejection for "+label+": "+message);
			rejected=true;
		}
		check(rejected, "Private mode entered "+label);
	}

	/** Opens an ordinary seed scope, checks input preservation and closes that scope in finally. */
	private static void ordinaryScope(final GeneCaller caller){
		final Read read=read();
		final byte[] before=read.bases.clone();
		caller.beginNcrnaIteration(read);
		try{
			check(Arrays.equals(before, read.bases), "Ordinary seed scope changed the input");
		}finally{
			caller.endNcrnaIteration();
		}
	}

	/** Throws AssertionError on a failed expectation, independent of -ea. */
	private static void check(boolean value, String message){
		if(!value){throw new AssertionError(message);}
	}
}
