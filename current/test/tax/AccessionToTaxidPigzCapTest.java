package test.tax;

import java.lang.reflect.Field;
import java.nio.file.Files;
import java.nio.file.Path;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import tax.AccessionToTaxid;

/** Constructor-only fixture for AccessionToTaxid's maxpigzprocesses cap.
 * @author Jean */
public final class AccessionToTaxidPigzCapTest{

	/** Runs baseline or fixed expectations for gzip input subprocess flags.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If fixture setup, reflection or cleanup fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		final FlagSnapshot originalFlags=FlagSnapshot.capture();
		final int originalMaxPigz=maxPigzProcesses();
		final float oldPrealloc=AccessionToTaxid.prealloc;
		final Path root=Files.createTempDirectory("accession-pigz-cap-");
		try{
			final Path gz0=touch(root, "00.tsv.gz");
			final Path gz1=touch(root, "01.tsv.gz");
			final Path gz2=touch(root, "02.tsv.gz");
			final Path gz3=touch(root, "03.tsv.gz");
			final Path plain0=touch(root, "10.tsv");
			final Path plain1=touch(root, "11.tsv");

			AccessionToTaxid.prealloc=0;
			checkGzipCap(new Path[] {gz0, gz1, gz2, gz3}, baseline);
			checkUncompressedControl(new Path[] {plain0, plain1});
		}finally{
			AccessionToTaxid.prealloc=oldPrealloc;
			setMaxPigzProcesses(originalMaxPigz);
			originalFlags.restore();
			deleteTempTree(root);
		}
		System.out.println("ACCESSION_PIGZ_CAP_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private AccessionToTaxidPigzCapTest(){}

	/** Verifies maxpigzprocesses chooses which gzip inputs request subprocesses.
	 * The effective FileFormat policy currently ignores that requested permission,
	 * so a counter repair alone cannot enforce the cap.
	 * @param files Empty .gz files passed to the constructor only
	 * @param baseline True to expect legacy all-requested behavior
	 * @throws Exception If construction or reflection fails */
	private static void checkGzipCap(final Path[] files, final boolean baseline)
			throws Exception{
		final AccessionToTaxid loader=construct(files, 2);
		final boolean[] requested=requestedSubprocess(loader);
		final boolean[] effective=effectiveSubprocess(loader);
		final boolean[] expected=(baseline ?
				new boolean[] {true, true, true, true} :
				new boolean[] {true, true, false, false});
		checkArray(requested, expected,
				"maxpigzprocesses=2 should request subprocesses only for first two gz inputs");
		check(allEqual(effective),
				"Effective FileFormat subprocess policy should ignore requested per-file cap");
		System.out.println("pigz_cap_all_requested="+allTrue(requested)+
				" effective_all_equal="+allEqual(effective)+
				" expect="+(baseline ? "baseline" : "fixed"));
	}

	/** Verifies uncompressed inputs request but do not effectively use subprocesses.
	 * @param files Empty uncompressed paths
	 * @throws Exception If construction or reflection fails */
	private static void checkUncompressedControl(final Path[] files) throws Exception{
		final AccessionToTaxid loader=construct(files, 0);
		checkArray(requestedSubprocess(loader), new boolean[] {true, true},
				"Plain text inputs should still request ordinary subprocess permission");
		checkArray(effectiveSubprocess(loader), new boolean[] {false, false},
				"Plain text inputs should not effectively use pigz/unpigz subprocesses");
	}

	/** Constructs a loader without calling process().
	 * @param files Input paths for the in= argument
	 * @param maxPigzProcesses maxpigzprocesses argument
	 * @return Constructed loader with populated FileFormat array */
	private static AccessionToTaxid construct(final Path[] files, final int maxPigzProcesses){
		setSentinelFlags();
		return new AccessionToTaxid(new String[] {
				"in="+join(files),
				"maxpigzprocesses="+maxPigzProcesses,
				"prealloc=0"});
	}

	/** Sets constructor globals so gzip input requests can record subprocess permission. */
	private static void setSentinelFlags(){
		ByteFile.FORCE_MODE_BF2=false;
		ByteFile.FORCE_MODE_BF1=false;
		ReadWrite.USE_UNPIGZ=true;
	}

	/** Extracts requested subprocess flags from the loader's private FileFormat array.
	 * @param loader Constructor-only AccessionToTaxid instance
	 * @return One FileFormat allowSubprocess field value per input
	 * @throws Exception If reflection fails */
	private static boolean[] requestedSubprocess(final AccessionToTaxid loader)
			throws Exception{
		final FileFormat[] formats=ffin(loader);
		final boolean[] flags=new boolean[formats.length];
		final Field field=FileFormat.class.getDeclaredField("allowSubprocess");
		field.setAccessible(true);
		for(int i=0; i<formats.length; i++){
			flags[i]=field.getBoolean(formats[i]);
		}
		return flags;
	}

	/** Extracts effective subprocess flags from the loader's private FileFormat array.
	 * @param loader Constructor-only AccessionToTaxid instance
	 * @return One FileFormat.allowSubprocess() value per input
	 * @throws Exception If reflection fails */
	private static boolean[] effectiveSubprocess(final AccessionToTaxid loader)
			throws Exception{
		final FileFormat[] formats=ffin(loader);
		final boolean[] flags=new boolean[formats.length];
		for(int i=0; i<formats.length; i++){
			flags[i]=formats[i].allowSubprocess();
		}
		return flags;
	}

	/** Reflects the private final ffin array.
	 * @param loader Constructor-only AccessionToTaxid instance
	 * @return FileFormat array populated by the constructor
	 * @throws Exception If reflection fails */
	private static FileFormat[] ffin(final AccessionToTaxid loader) throws Exception{
		final Field field=AccessionToTaxid.class.getDeclaredField("ffin");
		field.setAccessible(true);
		return (FileFormat[])field.get(loader);
	}

	/** Reads AccessionToTaxid.maxPigzProcesses.
	 * @return Current static max pigz process cap
	 * @throws Exception If reflection fails */
	private static int maxPigzProcesses() throws Exception{
		final Field field=AccessionToTaxid.class.getDeclaredField("maxPigzProcesses");
		field.setAccessible(true);
		return field.getInt(null);
	}

	/** Writes AccessionToTaxid.maxPigzProcesses.
	 * @param value Restored static max pigz process cap
	 * @throws Exception If reflection fails */
	private static void setMaxPigzProcesses(final int value) throws Exception{
		final Field field=AccessionToTaxid.class.getDeclaredField("maxPigzProcesses");
		field.setAccessible(true);
		field.setInt(null, value);
	}

	/** Creates an empty owned input file.
	 * @param root Temporary directory
	 * @param name File name to create
	 * @return Created file
	 * @throws Exception If file creation fails */
	private static Path touch(final Path root, final String name) throws Exception{
		return Files.createFile(root.resolve(name));
	}

	/** Joins paths into an AccessionToTaxid in= argument.
	 * @param files Paths to join
	 * @return Comma-delimited path list */
	private static String join(final Path[] files){
		final StringBuilder sb=new StringBuilder();
		for(Path file : files){
			if(sb.length()>0){sb.append(',');}
			sb.append(file.toString());
		}
		return sb.toString();
	}

	/** Reports whether every array value is true.
	 * @param array Array to inspect
	 * @return True if every value is true */
	private static boolean allTrue(final boolean[] array){
		for(boolean b : array){
			if(!b){return false;}
		}
		return true;
	}

	/** Reports whether every array value equals the first value.
	 * @param array Array to inspect
	 * @return True if values are uniform */
	private static boolean allEqual(final boolean[] array){
		for(int i=1; i<array.length; i++){
			if(array[i]!=array[0]){return false;}
		}
		return true;
	}

	/** Requires identical boolean arrays.
	 * @param actual Observed values
	 * @param expected Expected values
	 * @param message Failure message */
	private static void checkArray(final boolean[] actual, final boolean[] expected,
			final String message){
		check(actual.length==expected.length,
				message+"; length "+actual.length+" != "+expected.length);
		for(int i=0; i<actual.length; i++){
			check(actual[i]==expected[i],
					message+"; index "+i+" expected "+expected[i]+" found "+actual[i]);
		}
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

	/** Independent snapshot of globals the constructor mutates. */
	private static final class FlagSnapshot{

		/** Captures the current flag state.
		 * @return Snapshot of constructor-touched globals */
		static FlagSnapshot capture(){
			return new FlagSnapshot(ByteFile.FORCE_MODE_BF2, ByteFile.FORCE_MODE_BF1,
					ReadWrite.USE_UNPIGZ);
		}

		/** Stores constructor-touched flags.
		 * @param forceModeBf2_ ByteFile.FORCE_MODE_BF2 value
		 * @param forceModeBf1_ ByteFile.FORCE_MODE_BF1 value
		 * @param useUnpigz_ ReadWrite.USE_UNPIGZ value */
		FlagSnapshot(final boolean forceModeBf2_, final boolean forceModeBf1_,
				final boolean useUnpigz_){
			forceModeBf2=forceModeBf2_;
			forceModeBf1=forceModeBf1_;
			useUnpigz=useUnpigz_;
		}

		/** Restores this snapshot to process globals. */
		void restore(){
			ByteFile.FORCE_MODE_BF2=forceModeBf2;
			ByteFile.FORCE_MODE_BF1=forceModeBf1;
			ReadWrite.USE_UNPIGZ=useUnpigz;
		}

		/** ByteFile.FORCE_MODE_BF2 value. */
		final boolean forceModeBf2;
		/** ByteFile.FORCE_MODE_BF1 value. */
		final boolean forceModeBf1;
		/** ReadWrite.USE_UNPIGZ value. */
		final boolean useUnpigz;
	}
}
