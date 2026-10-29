package test.tax;

import java.io.RandomAccessFile;
import java.lang.reflect.Field;
import java.nio.MappedByteBuffer;
import java.nio.file.Files;
import java.nio.file.Path;

import tax.DiskAccessionTable;

/** Injected DiskAccessionTable.Builder close-failure cleanup check.
 * @author Jean */
public final class DiskAccessionBuilderCloseTest{

	/** Runs baseline or fixed expectations for Builder.close cleanup.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If fixture setup, reflection or cleanup fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		final Path root=Files.createTempDirectory("disk-accession-builder-close-");
		try{
			checkHealthyClose(root.resolve("healthy.bin"));
			checkInjectedEarlyCloseFailure(root.resolve("injected.bin"), baseline);
		}finally{
			deleteTempTree(root);
		}
		System.out.println("DISK_ACCESSION_BUILDER_CLOSE_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private DiskAccessionBuilderCloseTest(){}

	/** Verifies normal Builder.close writes a readable table.
	 * @param path Destination table
	 * @throws Exception If build or lookup fails */
	private static void checkHealthyClose(final Path path) throws Exception{
		try(DiskAccessionTable.Builder builder=DiskAccessionTable.build(path.toString(), 1)){
			check(builder.add(KEY, TAX_ID), "Healthy builder did not add the fixture key");
		}
		checkEntryCount(path, 1);
		DiskAccessionTable table=null;
		try{
			table=new DiskAccessionTable(path.toString());
			check(table.get(KEY)==TAX_ID, "Healthy table lookup failed");
			check(table.get(MISSING_KEY)==-1, "Healthy table missing-key lookup failed");
		}finally{
			if(table!=null){table.close();}
		}
	}

	/**
	 * Injects a null mapped segment to exercise the real close force loop.
	 *
	 * This proves cleanup and retry behavior after an artificial unchecked
	 * exception in {@link DiskAccessionTable.Builder#close()}. It does not prove
	 * that a real {@link MappedByteBuffer#force()} or header-write
	 * {@link java.io.IOException} is reachable, and it does not prove that
	 * production code can null the private segment array.
	 *
	 * @param path Destination table
	 * @param baseline True to expect the current leaked file descriptor
	 * @throws Exception If build, reflection or cleanup fails
	 */
	private static void checkInjectedEarlyCloseFailure(final Path path,
			final boolean baseline) throws Exception{
		final DiskAccessionTable.Builder builder=
				DiskAccessionTable.build(path.toString(), 1);
		RandomAccessFile raf=null;
		try{
			raf=randomAccessFile(builder);
			check(builder.add(KEY, TAX_ID), "Injected builder did not add the fixture key");
			final MappedByteBuffer[][] parts=parts(builder);
			parts[0][0]=null;

			expectNullPartFailure(builder, "initial");
			check(closed(builder), "Initial close failure should set the closed latch");
			check(raf.getFD().valid()==baseline,
					"Initial close failure descriptor state; baseline="+baseline+
					" valid="+raf.getFD().valid());

			builder.close();
			check(raf.getFD().valid()==baseline,
					"Repeated close descriptor state; baseline="+baseline+
					" valid="+raf.getFD().valid());
		}finally{
			if(raf==null){builder.close();}
			else{raf.close();}
		}
	}

	/** Verifies the serialized entry count written by Builder.close.
	 * @param path Complete table path
	 * @param expected Expected serialized entry count
	 * @throws Exception If the header cannot be read */
	private static void checkEntryCount(final Path path, final long expected)
			throws Exception{
		try(RandomAccessFile raf=new RandomAccessFile(path.toFile(), "r")){
			raf.seek(8);
			check(raf.readLong()==expected, "Serialized entry count mismatch");
		}
	}

	/** Expects the injected null segment to be the close failure.
	 * @param builder Builder with an injected null mapped segment
	 * @param label Call label for diagnostics
	 * @throws Exception If close fails with the wrong exception */
	private static void expectNullPartFailure(final DiskAccessionTable.Builder builder,
			final String label) throws Exception{
		try{
			builder.close();
			throw new AssertionError(label+" close should fail on the injected segment");
		}catch(NullPointerException expected){
			check(expected.getStackTrace().length>0 &&
					expected.getStackTrace()[0].getClassName().equals(
							"tax.DiskAccessionTable$Builder") &&
					expected.getStackTrace()[0].getMethodName().equals("close"),
					label+" close failed outside Builder.close");
		}
	}

	/** Reflects the builder output file.
	 * @param builder Builder to inspect
	 * @return Reflected random-access file
	 * @throws Exception If reflection fails */
	private static RandomAccessFile randomAccessFile(
			final DiskAccessionTable.Builder builder) throws Exception{
		return (RandomAccessFile)field("raf").get(builder);
	}

	/** Reflects mapped segments from the builder.
	 * @param builder Builder to inspect and fault-inject
	 * @return Reflected mapped segment matrix
	 * @throws Exception If reflection fails */
	private static MappedByteBuffer[][] parts(
			final DiskAccessionTable.Builder builder) throws Exception{
		return (MappedByteBuffer[][])field("parts").get(builder);
	}

	/** Reflects the close-entry latch from the builder.
	 * @param builder Builder to inspect
	 * @return True after Builder.close has entered
	 * @throws Exception If reflection fails */
	private static boolean closed(final DiskAccessionTable.Builder builder)
			throws Exception{
		return ((Boolean)field("closed").get(builder)).booleanValue();
	}

	/** Finds a private Builder field and makes it readable.
	 * @param name Private Builder field name
	 * @return Accessible field
	 * @throws NoSuchFieldException If the field is absent */
	private static Field field(final String name) throws NoSuchFieldException{
		Field field=DiskAccessionTable.Builder.class.getDeclaredField(name);
		field.setAccessible(true);
		return field;
	}

	/** Deletes the owned temporary tree, visiting children before parents.
	 * @param root Temporary fixture root
	 * @throws Exception If walking or deletion fails */
	private static void deleteTempTree(final Path root) throws Exception{
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths.sorted(
					java.util.Comparator.reverseOrder())::iterator){
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

	/** Fixture accession key. */
	private static final long KEY=123456789L;
	/** Missing accession key. */
	private static final long MISSING_KEY=987654321L;
	/** Fixture TaxID. */
	private static final int TAX_ID=4567;

}
