package test.tax;

import java.io.EOFException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.lang.reflect.Field;
import java.nio.BufferUnderflowException;
import java.nio.channels.AsynchronousFileChannel;
import java.nio.channels.FileChannel;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.util.Arrays;

import tax.AsyncDiskAccessor;
import tax.DiskAccessionTable;
import tax.NIODiskAccessor;

/** Offline constructor-failure checks for disk accession handle cleanup.
 * @author Jean */
public final class DiskAccessorConstructorCloseTest{

	/** Runs baseline or fixed constructor-close expectations.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If fixture setup fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		final Path root=Files.createTempDirectory("disk-accessor-constructor-close-");
		try{
			final Path healthy=root.resolve("healthy-accession.bin");
			writeHealthyTable(healthy);
			checkHealthyControls(healthy);

			checkFailedTableConstructor(root.resolve("malformed-table.bin"), baseline);
			checkFailedNioConstructor(root.resolve("malformed-nio.bin"), baseline);
			checkFailedAsyncConstructor(root.resolve("malformed-async.bin"), baseline);
		}finally{
			deleteTempTree(root);
		}
		System.out.println("DISK_ACCESSOR_CONSTRUCTOR_CLOSE_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private DiskAccessorConstructorCloseTest(){}

	/** Writes a valid tiny table through the production builder.
	 * @param path Destination table
	 * @throws IOException If the table cannot be written */
	private static void writeHealthyTable(final Path path) throws IOException{
		try(DiskAccessionTable.Builder builder=
				DiskAccessionTable.build(path.toString(), HEALTHY_KEYS.length)){
			for(int i=0; i<HEALTHY_KEYS.length; i++){
				check(builder.add(HEALTHY_KEYS[i], HEALTHY_VALUES[i]),
						"Failed to add healthy literal key "+HEALTHY_KEYS[i]);
			}
		}
	}

	/** Verifies all three constructors close owned handles on a valid tiny table.
	 * @param table Valid tiny accession table
	 * @throws Exception If construction, lookup or close fails */
	private static void checkHealthyControls(final Path table) throws Exception{
		DiskAccessionTable dat=null;
		try{
			dat=new DiskAccessionTable(table.toString());
			check(dat.get(HEALTHY_KEYS[0])==HEALTHY_VALUES[0],
					"Memory-mapped first-key lookup control failed");
			check(dat.get(HEALTHY_KEYS[1])==HEALTHY_VALUES[1],
					"Memory-mapped second-key lookup control failed");
			check(dat.get(MISSING_KEY)==-1,
					"Memory-mapped missing-key control failed");
			final FileChannel channel=fileChannel(dat, "channel");
			final RandomAccessFile raf=randomAccessFile(dat, "raf");
			dat.close();
			check(!channel.isOpen(), "Healthy DiskAccessionTable channel stayed open");
			check(!raf.getFD().valid(),
					"Healthy DiskAccessionTable file handle stayed open");
			dat=null;
		}finally{
			if(dat!=null){dat.close();}
		}

		NIODiskAccessor nio=null;
		try{
			nio=new NIODiskAccessor(table.toString(), 1);
			check(Arrays.equals(nio.getBatch(new long[] {HEALTHY_KEYS[0],
					HEALTHY_KEYS[2], MISSING_KEY}), new int[] {HEALTHY_VALUES[0],
					HEALTHY_VALUES[2], -1}), "NIO lookup control failed");
			final FileChannel channel=fileChannel(nio, "channel");
			final RandomAccessFile raf=randomAccessFile(nio, "raf");
			nio.close();
			check(!channel.isOpen(), "Healthy NIODiskAccessor channel stayed open");
			check(!raf.getFD().valid(),
					"Healthy NIODiskAccessor file handle stayed open");
			nio=null;
		}finally{
			if(nio!=null){nio.close();}
		}

		AsyncDiskAccessor async=null;
		try{
			async=new AsyncDiskAccessor(table.toString());
			check(Arrays.equals(async.getBatch(new long[] {HEALTHY_KEYS[1],
					HEALTHY_KEYS[2], MISSING_KEY}), new int[] {HEALTHY_VALUES[1],
					HEALTHY_VALUES[2], -1}), "Async lookup control failed");
			final AsynchronousFileChannel afc=asyncFileChannel(async, "afc");
			async.close();
			check(!afc.isOpen(), "Healthy AsyncDiskAccessor channel stayed open");
			async=null;
		}finally{
			if(async!=null){async.close();}
		}
	}

	/** Checks a failed memory-mapped table constructor for leaked handles.
	 * @param malformed Tiny malformed table path
	 * @param baseline True to expect the current leak
	 * @throws Exception If fixture setup fails */
	private static void checkFailedTableConstructor(final Path malformed,
			final boolean baseline) throws Exception{
		writeMalformedHeader(malformed);
		final int before=openHandles(malformed);
		expectConstructorFailure("DiskAccessionTable", EOFException.class,
				new ThrowingRunnable(){
					@Override
					public void run() throws Exception{
						new DiskAccessionTable(malformed.toString());
					}
				});
		checkConstructorLeak("DiskAccessionTable", malformed, before, 1, baseline);
	}

	/** Checks a failed positioned-read accessor constructor for leaked handles.
	 * @param malformed Tiny malformed table path
	 * @param baseline True to expect the current leak
	 * @throws Exception If fixture setup fails */
	private static void checkFailedNioConstructor(final Path malformed,
			final boolean baseline) throws Exception{
		writeMalformedHeader(malformed);
		final int before=openHandles(malformed);
		expectConstructorFailure("NIODiskAccessor", BufferUnderflowException.class,
				new ThrowingRunnable(){
					@Override
					public void run() throws Exception{
						new NIODiskAccessor(malformed.toString(), 1);
					}
				});
		checkConstructorLeak("NIODiskAccessor", malformed, before, 1, baseline);
	}

	/** Checks a failed asynchronous accessor constructor for leaked handles.
	 * @param malformed Tiny malformed table path
	 * @param baseline True to expect the current leak
	 * @throws Exception If fixture setup fails */
	private static void checkFailedAsyncConstructor(final Path malformed,
			final boolean baseline) throws Exception{
		writeMalformedHeader(malformed);
		final int before=openHandles(malformed);
		expectConstructorFailure("AsyncDiskAccessor", EOFException.class,
				new ThrowingRunnable(){
					@Override
					public void run() throws Exception{
						new AsyncDiskAccessor(malformed.toString());
					}
				});
		checkConstructorLeak("AsyncDiskAccessor", malformed, before, 2, baseline);
	}

	/** Writes a tiny header that opens successfully, then fails metadata reads.
	 * @param path Destination malformed table path
	 * @throws IOException If the header cannot be written */
	private static void writeMalformedHeader(final Path path) throws IOException{
		Files.write(path, new byte[] {0, 0, 0, 4});
	}

	/** Runs a constructor expected to fail while reading table metadata.
	 * @param label Constructor label
	 * @param expectedType Expected metadata-read failure type
	 * @param action Constructor invocation
	 * @throws Exception If the constructor throws an unexpected checked failure */
	private static void expectConstructorFailure(final String label,
			final Class<? extends Throwable> expectedType,
			final ThrowingRunnable action) throws Exception{
		Throwable thrown=null;
		try{
			action.run();
		}catch(Throwable t){
			thrown=t;
		}
		check(thrown!=null, label+" malformed metadata must fail");
		if(!expectedType.isInstance(thrown)){
			if(thrown instanceof Exception){throw (Exception)thrown;}
			if(thrown instanceof Error){throw (Error)thrown;}
			throw new AssertionError(thrown);
		}
	}

	/** Applies baseline or fixed file-handle expectations after construction fails.
	 * @param label Constructor label
	 * @param path Malformed table path
	 * @param before Open handles before construction
	 * @param baselineLeaks Exact number of file handles expected to leak now
	 * @param baseline True to expect the current leak
	 * @throws IOException If handle counting fails */
	private static void checkConstructorLeak(final String label, final Path path,
			final int before, final int baselineLeaks, final boolean baseline)
			throws IOException{
		final int after=openHandles(path);
		final int leaked=after-before;
		final int expectedLeaks=(baseline ? baselineLeaks : 0);
		System.out.println(label+" malformed constructor leaked_handles="+
				leaked+" expect="+expectedLeaks);
		check(leaked==expectedLeaks, label+" expected "+expectedLeaks+
				" leaked malformed-file handles but saw "+leaked+
				"; before="+before+" after="+after+" baseline="+baseline);
	}

	/** Reflects an owned FileChannel from a disk accessor.
	 * @param owner Object that owns the channel field
	 * @param name Private field name
	 * @return Reflected file channel
	 * @throws Exception If reflection fails */
	private static FileChannel fileChannel(final Object owner, final String name)
			throws Exception{
		return (FileChannel)field(owner, name).get(owner);
	}

	/** Reflects an owned RandomAccessFile from a disk accessor.
	 * @param owner Object that owns the file field
	 * @param name Private field name
	 * @return Reflected random-access file
	 * @throws Exception If reflection fails */
	private static RandomAccessFile randomAccessFile(final Object owner,
			final String name) throws Exception{
		return (RandomAccessFile)field(owner, name).get(owner);
	}

	/** Reflects an owned AsynchronousFileChannel from a disk accessor.
	 * @param owner Object that owns the channel field
	 * @param name Private field name
	 * @return Reflected asynchronous file channel
	 * @throws Exception If reflection fails */
	private static AsynchronousFileChannel asyncFileChannel(final Object owner,
			final String name) throws Exception{
		return (AsynchronousFileChannel)field(owner, name).get(owner);
	}

	/** Finds a private field and makes it readable.
	 * @param owner Object declaring the field
	 * @param name Private field name
	 * @return Accessible field
	 * @throws NoSuchFieldException If the field is absent */
	private static Field field(final Object owner, final String name)
			throws NoSuchFieldException{
		Field field=owner.getClass().getDeclaredField(name);
		field.setAccessible(true);
		return field;
	}

	/** Counts file descriptors in this process that target a path.
	 * @param path Target file
	 * @return Current open file descriptor count for that exact file
	 * @throws IOException If /proc/self/fd cannot be scanned */
	private static int openHandles(final Path path) throws IOException{
		final Path target=path.toRealPath();
		int count=0;
		try(DirectoryStream<Path> stream=Files.newDirectoryStream(PROC_FD)){
			for(Path fd : stream){
				try{
					if(Files.readSymbolicLink(fd).equals(target)){count++;}
				}catch(NoSuchFileException ignored){}
			}
		}
		return count;
	}

	/** Deletes temporary test files.
	 * @param root Temporary directory root
	 * @throws IOException If cleanup fails */
	private static void deleteTempTree(final Path root) throws IOException{
		try(DirectoryStream<Path> stream=Files.newDirectoryStream(root)){
			for(Path path : stream){Files.deleteIfExists(path);}
		}
		Files.deleteIfExists(root);
	}

	/** Throws AssertionError when a test condition fails.
	 * @param condition Condition that must be true
	 * @param message Failure message */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Constructor invocation that can throw any metadata-read failure. */
	private interface ThrowingRunnable{
		/** Runs the constructor under test.
		 * @throws Exception If the constructor reports a checked failure */
		void run() throws Exception;
	}

	/** Linux file-descriptor directory for the running fixture process. */
	private static final Path PROC_FD=new java.io.File("/proc/self/fd").toPath();
	/** Distinct literal keys in the healthy builder-produced table. */
	private static final long[] HEALTHY_KEYS={1, 2, 3};
	/** Expected taxIDs for {@link #HEALTHY_KEYS}. */
	private static final int[] HEALTHY_VALUES={11, 22, 33};
	/** Key absent from the healthy table. */
	private static final long MISSING_KEY=4;

}
