package test.tax;

import java.io.DataOutputStream;
import java.io.IOException;
import java.lang.reflect.Field;
import java.nio.file.Files;
import java.nio.file.Path;

import fileIO.ReadWrite;
import tax.DiskGiTable;
import tax.GiToTaxid;

/** Same-JVM GI table lifecycle checks for disk and in-memory reloads.
 * @author Jean */
public final class GiToTaxidLifecycleTest{

	/** Runs baseline or fixed expectations for disk-to-memory reload state.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If fixture construction fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		final Path root=Files.createTempDirectory("gi-lifecycle-");
		try{
			final Path disk=root.resolve("gitable.disk");
			final Path array=root.resolve("gitable.int2d");
			writeDiskGiTable(disk);
			writeArrayGiTable(array);
			checkCloseFailure(disk, baseline);

			GiToTaxid.DISK_MODE=false;
			GiToTaxid.unload();
			check(!GiToTaxid.isInitialized(),
					"GI table should be uninitialized after unload before disk load");
			GiToTaxid.initializeDisk(disk.toString());
			check(GiToTaxid.isInitialized(),
					"Disk GI table should be initialized after initializeDisk()");
			check(GiToTaxid.getID(3)==333 && GiToTaxid.getID(7)==777,
					"Disk GI table control failed");

			GiToTaxid.unload();
			check(!GiToTaxid.isInitialized(),
					"GI table should be uninitialized after disk unload");
			GiToTaxid.initialize(array.toString());
			check(GiToTaxid.isInitialized(),
					"In-memory GI table should be initialized after initialize()");

			final int gi3=GiToTaxid.getID(3);
			final int gi7=GiToTaxid.getID(7);
			if(baseline){
				check(gi3==333 && gi7==777,
						"Baseline must expose stale disk routing after unload: "+
						gi3+", "+gi7);
			}else{
				check(gi3==444 && gi7==888,
						"initialize() must load the new in-memory table after disk unload: "+
						gi3+", "+gi7);
			}
		}finally{
			closeRemainingDiskTable(null);
			GiToTaxid.DISK_MODE=false;
			GiToTaxid.unload();
			deleteTempTree(root);
		}
		System.out.println("GI_TO_TAXID_LIFECYCLE_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private GiToTaxidLifecycleTest(){}

	/** Writes a small memory-mappable GI table with known IDs.
	 * @param path Output disk-table path
	 * @throws Exception If the table cannot be written */
	private static void writeDiskGiTable(final Path path) throws Exception{
		try(DataOutputStream out=new DataOutputStream(Files.newOutputStream(path))){
			out.writeLong(10);
			out.writeLong(0);
			for(int gi=0; gi<=10; gi++){
				final int tid=(gi==3 ? 333 : gi==7 ? 777 : 0);
				out.writeInt(tid);
			}
		}
	}

	/** Writes a serialized in-memory GI table with IDs that differ from disk.
	 * @param path Output serialized table path
	 * @throws Exception If the table cannot be written */
	private static void writeArrayGiTable(final Path path) throws Exception{
		final int[][] array=new int[][] {new int[8]};
		array[0][3]=444;
		array[0][7]=888;
		ReadWrite.write(array, path.toString(), true);
	}

	/** Injects a disk-table close failure and checks unload cleanup.
	 * @param path Existing disk table path for subclass construction
	 * @param baseline True to expect legacy unload to ignore the disk table
	 * @throws Exception If reflection or construction fails */
	private static void checkCloseFailure(final Path path, final boolean baseline)
			throws Exception{
		DiskGiTable oldTable=null;
		try{
			GiToTaxid.initializeDisk(path.toString());
			oldTable=diskTable();
			setDiskTable(new ThrowingDiskGiTable(path.toString()));
			closeOwnedTable(oldTable);
			oldTable=null;
			try{
				GiToTaxid.unload();
				if(!baseline){throw new AssertionError("unload() must report disk close failure");}
				check(GiToTaxid.DISK_MODE,
						"Baseline must expose that unload() left disk mode enabled");
			}catch(RuntimeException e){
				if(baseline){throw e;}
				check(e.getCause() instanceof IOException,
						"Close failure should retain IOException cause, not "+e.getCause());
				check(!GiToTaxid.isInitialized(),
						"unload() should clear initialized state after close failure");
				check(!GiToTaxid.DISK_MODE,
						"unload() should leave disk mode disabled after close failure");
				check(diskTable()==null,
						"unload() should clear the disk table reference before close failure");
			}
		}finally{
			closeRemainingDiskTable(oldTable);
			GiToTaxid.DISK_MODE=false;
			GiToTaxid.unload();
		}
	}

	/** Clears and closes any disk table still owned by GiToTaxid.
	 * @param extraTable Extra table opened by this fixture and not currently
	 * attached to GiToTaxid; null is accepted
	 * @throws Exception If a non-injected close fails */
	private static void closeRemainingDiskTable(final DiskGiTable extraTable)
			throws Exception{
		final DiskGiTable current=diskTable();
		setDiskTable(null);
		if(current!=null){closeOwnedTable(current);}
		if(extraTable!=null && extraTable!=current){closeOwnedTable(extraTable);}
	}

	/** Closes a table opened by this fixture.
	 * The only swallowed failure is the fixture's deliberate injected close error.
	 * @param table Table to close
	 * @throws IOException If a non-injected close fails */
	private static void closeOwnedTable(final DiskGiTable table) throws IOException{
		try{
			table.close();
		}catch(IOException e){
			if(!(table instanceof ThrowingDiskGiTable)){throw e;}
		}
	}

	/** Reads GiToTaxid's private disk table for lifecycle assertions.
	 * @return Active disk table, or null
	 * @throws Exception If reflection fails */
	private static DiskGiTable diskTable() throws Exception{
		return (DiskGiTable)DISK_TABLE_FIELD.get(null);
	}

	/** Injects the active private disk table.
	 * @param table Replacement table
	 * @throws Exception If reflection fails */
	private static void setDiskTable(final DiskGiTable table) throws Exception{
		DISK_TABLE_FIELD.set(null, table);
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

	/** Disk table whose close() always reports a finite cleanup failure. */
	private static final class ThrowingDiskGiTable extends DiskGiTable{

		/** Opens a valid disk table for injected close-failure checks.
		 * @param path Disk GI table path
		 * @throws IOException If the table cannot be opened */
		ThrowingDiskGiTable(final String path) throws IOException{super(path);}

		@Override
		public void close() throws IOException{
			super.close();
			throw new IOException("injected close failure");
		}
	}

	/** Keeps fixture checks active under -da.
	 * @param condition Required condition
	 * @param message Failure message */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Reflected GiToTaxid.diskTable field. */
	private static final Field DISK_TABLE_FIELD;

	static{
		try{
			DISK_TABLE_FIELD=GiToTaxid.class.getDeclaredField("diskTable");
			DISK_TABLE_FIELD.setAccessible(true);
		}catch(Exception e){throw new RuntimeException(e);}
	}
}
