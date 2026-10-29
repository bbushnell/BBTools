package test.tax;

import java.io.FileNotFoundException;
import java.nio.file.Files;
import java.nio.file.Path;

import fileIO.ByteFile;
import fileIO.ReadWrite;
import tax.AccessionToTaxid;

/** Failure-path fixture for AccessionToTaxid.load global flag restoration.
 * @author Jean */
public final class AccessionToTaxidLoadFailureTest{

	/** Runs baseline or fixed expectations for failed load restoration.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If fixture setup or cleanup fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		final FlagSnapshot original=FlagSnapshot.capture();
		final float oldPrealloc=AccessionToTaxid.prealloc;
		final Path root=Files.createTempDirectory("accession-load-failure-");
		try{
			final Path inputDirectory=root.resolve("input-is-directory.tsv");
			Files.createDirectory(inputDirectory);

			AccessionToTaxid.prealloc=0;
			setSentinelFlags();
			checkDirectoryOpenFailure(inputDirectory, baseline);
		}finally{
			AccessionToTaxid.prealloc=oldPrealloc;
			original.restore();
			deleteTempTree(root);
		}
		System.out.println("ACCESSION_LOAD_FAILURE_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private AccessionToTaxidLoadFailureTest(){}

	/** Sets the four globals that AccessionToTaxid.load snapshots.
	 * FORCE_MODE_BF1 is disabled so the loader's USE_UNPIGZ branch leaves a
	 * visible baseline mutation before pre-record-processing ByteFile construction fails. */
	private static void setSentinelFlags(){
		ByteFile.FORCE_MODE_BF2=false;
		ByteFile.FORCE_MODE_BF1=false;
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_GUNZIP=false;
	}

	/** Verifies that a failed load restores or leaks its saved globals.
	 * @param inputDirectory Directory path used as a deterministic invalid file
	 * @param baseline True to expect legacy failure-path leakage */
	private static void checkDirectoryOpenFailure(final Path inputDirectory,
			final boolean baseline){
		final FlagSnapshot before=FlagSnapshot.capture();
		RuntimeException failure=null;
		try{
			AccessionToTaxid.load(inputDirectory.toString());
		}catch(RuntimeException t){
			failure=t;
		}
		final FlagSnapshot after=FlagSnapshot.capture();

		check(failure!=null,
				"Directory input should fail before any accession records are processed");
		check(directoryOpenFailure(failure, inputDirectory),
				"Expected FileNotFoundException for "+inputDirectory+
				", observed "+failure);

		if(baseline){
			check(before.forceModeBf1==false && after.forceModeBf1==true,
					"Baseline should leave FORCE_MODE_BF1 enabled after failed load; before="+
					before+" after="+after);
		}else{
			after.checkMatches(before,
					"load(String) should restore all saved globals after failed load");
		}
		System.out.println("load_failure_restored="+after.matches(before)+
				" expect="+(baseline ? "baseline" : "fixed")+
				" failure="+failure.getClass().getSimpleName());
	}

	/** Verifies the expected directory-open failure from ByteFile construction.
	 * @param failure Runtime failure from {@link AccessionToTaxid#load(String)}
	 * @param inputDirectory Owned directory that should have failed to open as a file
	 * @return True if the cause chain names this exact directory in a FileNotFoundException */
	private static boolean directoryOpenFailure(final RuntimeException failure,
			final Path inputDirectory){
		for(Throwable t=failure; t!=null; t=t.getCause()){
			if(t instanceof FileNotFoundException){
				final String message=t.getMessage();
				return message!=null && message.contains(inputDirectory.toString());
			}
		}
		return false;
	}

	/** Deletes the owned temporary tree, visiting children before parents.
	 * @param root Temporary fixture root
	 * @throws Exception If walking or deletion fails */
	private static void deleteTempTree(final Path root) throws Exception{
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths.sorted(java.util.Comparator.reverseOrder())::iterator){
				Files.deleteIfExists(path);
			}
		}
	}

	/** Keeps fixture checks active under -da.
	 * @param condition Required condition
	 * @param message Failure message */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Independent snapshot of globals AccessionToTaxid.load must restore. */
	private static final class FlagSnapshot{

		/** Captures the current flag state.
		 * @return Snapshot of the four saved AccessionToTaxid.load globals */
		static FlagSnapshot capture(){
			return new FlagSnapshot(ByteFile.FORCE_MODE_BF2, ByteFile.FORCE_MODE_BF1,
					ReadWrite.USE_UNPIGZ, ReadWrite.USE_GUNZIP);
		}

		/** Stores all four saved flags independently.
		 * @param forceModeBf2_ ByteFile.FORCE_MODE_BF2 value
		 * @param forceModeBf1_ ByteFile.FORCE_MODE_BF1 value
		 * @param useUnpigz_ ReadWrite.USE_UNPIGZ value
		 * @param useGunzip_ ReadWrite.USE_GUNZIP value */
		FlagSnapshot(final boolean forceModeBf2_, final boolean forceModeBf1_,
				final boolean useUnpigz_, final boolean useGunzip_){
			forceModeBf2=forceModeBf2_;
			forceModeBf1=forceModeBf1_;
			useUnpigz=useUnpigz_;
			useGunzip=useGunzip_;
		}

		/** Restores this snapshot to the process globals. */
		void restore(){
			ByteFile.FORCE_MODE_BF2=forceModeBf2;
			ByteFile.FORCE_MODE_BF1=forceModeBf1;
			ReadWrite.USE_UNPIGZ=useUnpigz;
			ReadWrite.USE_GUNZIP=useGunzip;
		}

		/** Reports whether this snapshot equals another.
		 * @param expected Expected flag snapshot
		 * @return True if all flags match */
		boolean matches(final FlagSnapshot expected){
			return forceModeBf2==expected.forceModeBf2 &&
					forceModeBf1==expected.forceModeBf1 &&
					useUnpigz==expected.useUnpigz &&
					useGunzip==expected.useGunzip;
		}

		/** Checks that this snapshot equals an expected snapshot.
		 * @param expected Expected flag snapshot
		 * @param message Failure message prefix */
		void checkMatches(final FlagSnapshot expected, final String message){
			check(matches(expected), message+"; before="+expected+" after="+this);
		}

		/** Describes all saved flags for assertion diagnostics. */
		@Override
		public String toString(){
			return "BF2="+forceModeBf2+
					",BF1="+forceModeBf1+
					",UNPIGZ="+useUnpigz+
					",GUNZIP="+useGunzip;
		}

		/** ByteFile.FORCE_MODE_BF2 value */
		final boolean forceModeBf2;
		/** ByteFile.FORCE_MODE_BF1 value */
		final boolean forceModeBf1;
		/** ReadWrite.USE_UNPIGZ value */
		final boolean useUnpigz;
		/** ReadWrite.USE_GUNZIP value */
		final boolean useGunzip;
	}
}
