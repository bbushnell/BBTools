package test.tax;

import java.io.DataOutputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.lang.reflect.Field;
import java.nio.channels.FileChannel;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import tax.DiskGiTable;

/** Finite lookup and constructor-failure checks for DiskGiTable.
 * @author Jean */
public final class DiskGiTableFailureTest{

	/** Runs baseline or fixed expectations for constructor resource cleanup.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If fixture construction fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		final Path root=Files.createTempDirectory("disk-gi-table-");
		try{
			final Path valid=root.resolve("valid.gitable");
			final Path truncated=root.resolve("truncated.gitable");
			writeDiskGiTable(valid);
			writeTruncatedDiskGiTable(truncated);
			checkNormalOpenLookupClose(valid);
			checkTruncatedConstructorFailure(truncated, baseline);
		}finally{
			deleteTempTree(root);
		}
		System.out.println("DISK_GI_TABLE_FAILURE_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private DiskGiTableFailureTest(){}

	/** Writes a small complete disk table with known direct lookups.
	 * @param path Output table path
	 * @throws IOException If the table cannot be written */
	private static void writeDiskGiTable(final Path path) throws IOException{
		try(DataOutputStream out=new DataOutputStream(Files.newOutputStream(path))){
			out.writeLong(10);
			out.writeLong(0);
			for(int gi=0; gi<=10; gi++){
				final int tid=(gi==3 ? 333 : gi==7 ? 777 : 0);
				out.writeInt(tid);
			}
		}
	}

	/** Writes only the header of a table whose body is declared but absent.
	 * @param path Output table path
	 * @throws IOException If the table cannot be written */
	private static void writeTruncatedDiskGiTable(final Path path) throws IOException{
		try(DataOutputStream out=new DataOutputStream(Files.newOutputStream(path))){
			out.writeLong(10);
			out.writeLong(0);
		}
	}

	/** Verifies direct lookup and explicit close on a complete table.
	 * @param path Complete disk table path
	 * @throws Exception If construction, lookup or close fails */
	private static void checkNormalOpenLookupClose(final Path path) throws Exception{
		final DiskGiTable table=new DiskGiTable(path.toString());
		try{
			checkOwnedHandlesOpen(table, true);
			check(table.get(-1)==-1, "Negative GI should return -1");
			check(table.get(3)==333, "GI 3 should map to 333");
			check(table.get(7)==777, "GI 7 should map to 777");
			check(table.get(11)==-2, "Out-of-range GI should return -2");
		}finally{
			table.close();
		}
		checkOwnedHandlesOpen(table, false);
		System.out.println("normal_owned_handles_closed=true");
	}

	/** Verifies failed construction on a truncated table.
	 * @param path Truncated disk table path
	 * @param baseline True to expect legacy constructor cleanup leakage
	 * @throws Exception If descriptor inspection fails */
	private static void checkTruncatedConstructorFailure(final Path path,
			final boolean baseline) throws Exception{
		final int before=countOpenDescriptors(path);
		final int attempts=3;
		int failures=0;
		for(int i=0; i<attempts; i++){
			try{
				final DiskGiTable table=new DiskGiTable(path.toString());
				table.close();
				throw new AssertionError(
						"Truncated DiskGiTable construction should fail on attempt "+i);
			}catch(IOException expected){
				failures++;
			}
		}
		check(failures==attempts, "Expected "+attempts+
				" truncated constructor failures, observed "+failures);
		final int after=countOpenDescriptors(path);
		if(baseline){
			check(after>before,
					"Baseline should leak this fixture after constructor failure; before="+
					before+" after="+after);
		}else{
			check(after==before,
					"Constructor failure should close opened resources; before="+
					before+" after="+after);
		}
		System.out.println("constructor_failure_delta="+(after-before)+
				" expect="+(baseline ? "baseline" : "fixed"));
	}

	/** Counts this JVM's open descriptors for exactly one fixture path.
	 * @param target Fixture file to count
	 * @return Number of /proc/self/fd entries pointing to target
	 * @throws IOException If descriptor inspection is unavailable */
	private static int countOpenDescriptors(final Path target) throws IOException{
		check(Files.isDirectory(FD_ROOT),
				"DiskGiTableFailureTest requires /proc/self/fd");
		final String canonical=target.toRealPath().toString();
		int count=0;
		try(DirectoryStream<Path> fds=Files.newDirectoryStream(FD_ROOT)){
			for(Path fd : fds){
				try{
					if(Files.readSymbolicLink(fd).toString().equals(canonical)){count++;}
				}catch(IOException ignored){
					// Descriptor disappeared between directory listing and readlink.
				}
			}
		}
		return count;
	}

	/** Checks whether the table's owned channel and file handle are open.
	 * @param table Complete table to inspect
	 * @param expected Expected open state for both handles
	 * @throws Exception If reflection or descriptor inspection fails */
	private static void checkOwnedHandlesOpen(final DiskGiTable table,
			final boolean expected) throws Exception{
		final FileChannel channel=(FileChannel)field("channel").get(table);
		final RandomAccessFile raf=(RandomAccessFile)field("raf").get(table);
		final boolean channelOpen=channel.isOpen();
		final boolean rafOpen=raf.getFD().valid();
		check(channelOpen==expected && rafOpen==expected,
				"Owned DiskGiTable handles should both be "+(expected ? "open" : "closed")+
				"; channelOpen="+channelOpen+" rafOpen="+rafOpen);
	}

	/** Returns a private DiskGiTable field for read-only fixture inspection.
	 * @param name Field name
	 * @return Accessible field
	 * @throws NoSuchFieldException If the field is absent */
	private static Field field(final String name) throws NoSuchFieldException{
		final Field field=DiskGiTable.class.getDeclaredField(name);
		field.setAccessible(true);
		return field;
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

	/** Linux descriptor directory for this JVM only. */
	private static final Path FD_ROOT=Paths.get("/proc/self/fd");

}
