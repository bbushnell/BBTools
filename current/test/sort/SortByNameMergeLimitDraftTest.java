package test.sort;

import java.io.File;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;

import fileIO.FileFormat;
import sort.SortByName;

/** Draft coverage for SortByName recursive merge limits below configured maxfiles.
 * @author Jean
 */
public final class SortByNameMergeLimitDraftTest{

	/** Takes baseline/fixed, an unused runtime directory, and a fresh fixture directory. */
	public static void main(String[] args) throws Exception{
		if(args.length!=3){throw new IllegalArgumentException("Expected phase, runtime, and fresh fixture directory");}
		final String phase=args[0];
		final boolean fixed=phase.equals("fixed") || phase.equals("devfull-fixed");
		if(!fixed && !phase.equals("baseline") && !phase.equals("devfull-baseline")){
			throw new IllegalArgumentException("Unknown phase: "+phase);
		}
		final Path root=Paths.get(args[2]).toAbsolutePath();
		Files.createDirectory(root);
		System.out.println("case\tresult_size\tmerge_calls\tmax_batch");
		if(phase.startsWith("devfull-")){
			runDeleteAfterCloseFailure(root, fixed);
			return;
		}
		runComputedCap(root, fixed);
		runConfiguredCapControl(root);
		runConfiguredEntryControl(root, fixed);
		runMaxMergeFilesCap(fixed);
	}

	/** Ensures a computed merge cap below maxfiles still drives recursive reduction. */
	private static void runComputedCap(Path root, boolean fixed) throws Exception{
		final String name="computed-cap";
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Probe probe=new Probe(argsFor(dir, 100));
		try{
			withForcedObservedSize(new ReflectiveCall(){
				@Override
				public void run() throws Exception{
					final ArrayList<String> result=mergeRecursive(probe, fileNames(5));
					if(fixed){
						assert(result.size()==2 && probe.calls==5 && probe.maxBatch==2) :
							"Computed cap was not honored; result="+result.size()+
							" calls="+probe.calls+" maxBatch="+probe.maxBatch;
					}else{
						assert(result.size()==5 && probe.calls==0) :
							"Expected current mergeRecursive to ignore computed cap below maxfiles; result="+
							result.size()+" calls="+probe.calls;
					}
					System.out.println(name+'\t'+result.size()+'\t'+probe.calls+'\t'+probe.maxBatch);
				}
			});
		}finally{probe.deleteTemps();}
	}

	/** Shows ordinary recursion still works when the user-configured maxfiles cap is exceeded. */
	private static void runConfiguredCapControl(Path root) throws Exception{
		final String name="configured-cap-control";
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Probe probe=new Probe(argsFor(dir, 3));
		try{
			final ArrayList<String> result=mergeRecursive(probe, fileNames(5));
			assert(result.size()==2 && probe.calls==2 && probe.maxBatch==3) :
				"Configured maxfiles cap did not reduce five inputs to two groups; result="+
				result.size()+" calls="+probe.calls+" maxBatch="+probe.maxBatch;
			System.out.println(name+'\t'+result.size()+'\t'+probe.calls+'\t'+probe.maxBatch);
		}finally{probe.deleteTemps();}
	}

	/** Verifies the final merge entry point honors the configured maxfiles limit. */
	private static void runConfiguredEntryControl(Path root, boolean fixed) throws Exception{
		final String name="configured-entry-control";
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Probe probe=new Probe(argsFor(dir, 3), true);
		try{
			final boolean errorState=mergeAndDump(probe, realFileNames(dir, 5));
			final int expectedCalls=(fixed ? 2 : 0);
			final int expectedMaxBatch=(fixed ? 3 : 0);
			assert(!errorState && probe.calls==expectedCalls && probe.maxBatch==expectedMaxBatch) :
				"Configured final-merge gate did not reduce five inputs to two groups; errorState="+
				errorState+" calls="+probe.calls+" maxBatch="+probe.maxBatch;
			System.out.println(name+"\t0\t"+probe.calls+'\t'+probe.maxBatch);
		}finally{probe.deleteTemps();}
	}

	/** Preserves inputs when buffered output failure is only reported by poisonAndWait. */
	private static void runDeleteAfterCloseFailure(Path root, boolean fixed) throws Exception{
		final String name="delete-after-close-failure";
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Path input=dir.resolve("input.fq");
		write(input, READ);
		final ArrayList<String> names=new ArrayList<String>();
		names.add(input.toString());
		final FileFormat ffout=FileFormat.testOutput("/dev/full", FileFormat.FASTQ, null,
			true, true, false, false);
		final boolean errorState=SortByName.mergeAndDump(names, ffout, null, true,
			false, false, System.err, 4);
		final boolean retained=Files.exists(input);
		if(fixed){
			if(!errorState || !retained){
				throw new AssertionError("Expected failed merge to preserve input; errorState="+
					errorState+" retained="+retained);
			}
		}else if(!errorState || retained){
			throw new AssertionError("Expected current merge to delete before output close; errorState="+
				errorState+" retained="+retained);
		}
		System.out.println(name+"\t"+(retained ? 1 : 0)+"\t0\t0");
	}

	/** Constructs arguments for SortByName validation without processing reads. */
	private static String[] argsFor(Path dir, int maxfiles) throws Exception{
		final Path input=dir.resolve("input.fq"), output=dir.resolve("out.fq");
		write(input, READ);
		return new String[]{"in="+input, "out="+output, "delete=f", "overwrite=f",
			"maxfiles="+maxfiles};
	}

	/** Invokes the private recursive reducer under a no-I/O subclass hook. */
	@SuppressWarnings("unchecked")
	private static ArrayList<String> mergeRecursive(Probe probe, ArrayList<String> names){
		try{
			final Method method=SortByName.class.getDeclaredMethod("mergeRecursive",
				ArrayList.class, boolean.class);
			method.setAccessible(true);
			return (ArrayList<String>)method.invoke(probe, names, false);
		}catch(InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof RuntimeException){throw (RuntimeException)cause;}
			if(cause instanceof Error){throw (Error)cause;}
			throw new RuntimeException(cause);
		}catch(Exception e){
			throw new RuntimeException(e);
		}
	}

	/** Invokes the private final-merge wrapper that decides whether recursion is needed. */
	private static boolean mergeAndDump(Probe probe, ArrayList<String> names){
		try{
			final Method method=SortByName.class.getDeclaredMethod("mergeAndDump",
				ArrayList.class, Boolean.TYPE, Boolean.TYPE, Boolean.TYPE);
			method.setAccessible(true);
			return (Boolean)method.invoke(probe, names, false, false, false);
		}catch(InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof RuntimeException){throw (RuntimeException)cause;}
			if(cause instanceof Error){throw (Error)cause;}
			throw new RuntimeException(cause);
		}catch(Exception e){
			throw new RuntimeException(e);
		}
	}

	/** Forces maxMergeFiles(maxSizeObserved()) down to two without large fixture files. */
	private static void withForcedObservedSize(ReflectiveCall call) throws Exception{
		final Field maxLength=SortByName.class.getDeclaredField("maxLengthObservedStatic");
		final Field quality=SortByName.class.getDeclaredField("qualityObservedStatic");
		maxLength.setAccessible(true);
		quality.setAccessible(true);
		final long oldMax=maxLength.getLong(null);
		final boolean oldQuality=quality.getBoolean(null);
		try{
			maxLength.setLong(null, Long.MAX_VALUE/4);
			quality.setBoolean(null, false);
			call.run();
		}finally{
			maxLength.setLong(null, oldMax);
			quality.setBoolean(null, oldQuality);
		}
	}

	/** Exercises maxMergeFiles cap arithmetic without requiring a huge heap. */
	private static void runMaxMergeFilesCap(boolean fixed) throws Exception{
		if(!fixed){return;}
		final Method method=SortByName.class.getDeclaredMethod("maxMergeFiles", Long.TYPE, Long.TYPE);
		method.setAccessible(true);
		final int capped=(Integer)method.invoke(null, 0L, Long.MAX_VALUE);
		final int floor=(Integer)method.invoke(null, Long.MAX_VALUE, 0L);
		assert(capped==10000 && floor==2) :
			"Expected maxMergeFiles to cap in long space; capped="+capped+" floor="+floor;
		System.out.println("max-merge-files-cap\t"+capped+"\t0\t0");
	}

	/** Returns synthetic file names; the Probe hook prevents opening them. */
	private static ArrayList<String> fileNames(int count){
		final ArrayList<String> names=new ArrayList<String>();
		for(int i=0; i<count; i++){names.add("in"+i+".fq");}
		return names;
	}

	/** Writes sorted single-read FASTQs for the final-merge entry test. */
	private static ArrayList<String> realFileNames(Path dir, int count) throws Exception{
		final ArrayList<String> names=new ArrayList<String>();
		for(int i=0; i<count; i++){
			final Path file=dir.resolve("real"+i+".fq");
			write(file, "@read"+i+"\nACGT\n+\nIIII\n");
			names.add(file.toString());
		}
		return names;
	}

	private static void write(Path file, String text) throws Exception{
		Files.write(file, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	private static final String READ="@read0\nACGT\n+\nIIII\n";

	/** Callback that can throw checked reflection exceptions. */
	private static interface ReflectiveCall{
		void run() throws Exception;
	}

	/** Overrides SortByName's per-group merge to record recursion without opening inputs. */
	private static final class Probe extends SortByName{

		Probe(String[] args){
			this(args, false);
		}

		Probe(String[] args, boolean writeTemps_){
			super(args);
			writeTemps=writeTemps_;
		}

		@Override
		public void mergeRecursiveInner(ArrayList<String> inList, FileFormat ffoutTemp1,
				FileFormat ffoutTemp2, boolean deleteAfterMerge, boolean allowSubprocess){
			calls++;
			maxBatch=Math.max(maxBatch, inList.size());
			if(writeTemps){
				final boolean errorState=SortByName.mergeAndDump(inList, ffoutTemp1, ffoutTemp2,
					deleteAfterMerge, false, allowSubprocess, System.err, 4);
				if(errorState){throw new RuntimeException("Error while writing recursive temp merge.");}
			}
			if(ffoutTemp1!=null){
				final String name=ffoutTemp1.name();
				temps.add(name);
				if(!writeTemps){new File(name).delete();}
			}
		}

		/** Deletes temp files created by SortByName.getTempFile before this hook ran. */
		void deleteTemps(){
			for(String temp : temps){new File(temp).delete();}
		}

		int calls=0;
		int maxBatch=0;
		final boolean writeTemps;
		final ArrayList<String> temps=new ArrayList<String>();
	}
}
