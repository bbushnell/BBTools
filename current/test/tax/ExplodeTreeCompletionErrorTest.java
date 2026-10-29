package test.tax;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.zip.GZIPInputStream;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Timer;
import tax.ExplodeTree;
import tax.TaxTree;

/** Finite injected-reader check for ExplodeTree input close-status propagation.
 * @author Sucrose */
public final class ExplodeTreeCompletionErrorTest{

	/** Constructs this stateless fixture; checks are run through {@link #main(String[])}. */
	public ExplodeTreeCompletionErrorTest(){}

	/** Runs an ordinary success control and an injected close-failure case.
	 * @param args Optional {@code expect=baseline} or {@code expect=fixed} mode selector
	 * @throws Exception if setup, output reading, or cleanup fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		run(false, baseline);
		run(true, baseline);
		System.out.println("EXPLODETREE_COMPLETION_ERROR_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Runs one real ExplodeTree process over an injected finite ByteFile.
	 * @param closeError Whether the injected reader should report close failure
	 * @param baseline Whether to assert the baseline or fixed expectation
	 * @throws Exception if setup, output reading, or cleanup fails */
	private static void run(final boolean closeError, final boolean baseline) throws Exception{
		final Path root=Files.createTempDirectory("explodetree-completion-");
		try{
			final Path input=Files.createTempFile(root, "input", ".fa");
			Files.write(input, new byte[0]);
			final ByteFile reader=new InjectedByteFile(FileFormat.testInput(input.toString(),
					FileFormat.FASTA, null, true, true), closeError);
			final InjectedExplodeTree tree=new InjectedExplodeTree(input, root,
					TaxTreeContractTest.tree(TaxTree.GENUS_E), reader);

			Throwable failure=null;
			try{tree.process(new Timer());}
			catch(Throwable t){failure=t;}
			check(((InjectedByteFile)reader).closes==1, "ExplodeTree must close its input exactly once");
			check(readOnlyOutput(root).equals(EXPECTED_FASTA),
					"ExplodeTree must emit the injected FASTA before any input close-status failure");

			if(!closeError){
				check(failure==null, "Ordinary successful output must not fail: "+failure);
				check(!tree.errorState, "Ordinary successful output must not set errorState");
			}else if(baseline){
				check(failure==null, "Baseline must expose the discarded close-status failure");
				check(!tree.errorState, "Baseline must leave errorState false after a close failure");
			}else{
				check(failure instanceof RuntimeException,
						"Fixed ExplodeTree must fail loudly on input close failure: "+failure);
				check(FIXED_CLOSE_MESSAGE.equals(failure.getMessage()),
						"Fixed ExplodeTree must reuse the existing errorState exception text: "+failure);
				check(tree.errorState, "Fixed ExplodeTree must retain close failure in errorState");
			}
		}finally{
			deleteTree(root);
		}
	}

	/** Reads the single emitted FASTA gzip payload from the directory tree.
	 * @param root Temporary output root
	 * @return Decompressed FASTA text emitted by ExplodeTree
	 * @throws Exception if the output cannot be found or read */
	private static String readOnlyOutput(final Path root) throws Exception{
		Path out=null;
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths::iterator){
				if(path.getFileName().toString().equals("4.fa.gz")){
					check(out==null, "Expected exactly one output file named 4.fa.gz");
					out=path;
				}
			}
		}
		check(out!=null, "ExplodeTree must write taxid 4 output under the output tree");
		try(GZIPInputStream in=new GZIPInputStream(Files.newInputStream(out));
				ByteArrayOutputStream buffer=new ByteArrayOutputStream()){
			byte[] temp=new byte[256];
			for(int len=in.read(temp); len>=0; len=in.read(temp)){buffer.write(temp, 0, len);}
			return new String(buffer.toByteArray(), StandardCharsets.UTF_8);
		}
	}

	/** Deletes the owned temporary tree, visiting children before parents.
	 * @param root Temporary directory tree to delete
	 * @throws Exception if tree traversal or deletion fails */
	private static void deleteTree(final Path root) throws Exception{
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths.sorted(Comparator.reverseOrder())::iterator){
				Files.deleteIfExists(path);
			}
		}
	}

	/** Keeps checks active with assertions disabled.
	 * @param condition Invariant that must be true
	 * @param message Failure detail for the thrown assertion */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** ExplodeTree variant that returns the injected finite reader instead of opening input. */
	private static final class InjectedExplodeTree extends ExplodeTree{
		/** Creates an ExplodeTree wired to a finite injected reader.
		 * @param input Placeholder input file path for FileFormat validation
		 * @param out Temporary output root
		 * @param tree_ Injected taxonomy tree
		 * @param reader_ Injected input reader */
		InjectedExplodeTree(final Path input, final Path out, final TaxTree tree_, final ByteFile reader_){
			super(FileFormat.testInput(input.toString(), FileFormat.FASTA, null, true, true),
					tree_, out.toString()+File.separator, "", true);
			reader=reader_;
		}
		@Override
		protected ByteFile makeByteFile(final FileFormat ff){return reader;}
		/** Reader returned by {@link #makeByteFile(FileFormat)}. */
		final ByteFile reader;
	}

	/** Minimal ByteFile that supplies one FASTA record and then reports a chosen close status. */
	private static final class InjectedByteFile extends ByteFile{
		/** Creates a finite reader over {@link ExplodeTreeCompletionErrorTest#LINES}.
		 * @param ff Placeholder FileFormat passed to the ByteFile base class
		 * @param closeError_ Status to return from close */
		InjectedByteFile(final FileFormat ff, final boolean closeError_){
			super(ff);
			closeError=closeError_;
		}
		@Override
		public void reset(){next=0;}
		@Override
		public InputStream is(){return null;}
		@Override
		public long lineNum(){return next;}
		@Override
		public boolean close(){closes++; open=false; return closeError;}
		@Override
		public void pushBack(final byte[] line){throw new AssertionError("ExplodeTree does not use pushBack");}
		@Override
		public boolean isOpen(){return open;}
		@Override
		public byte[] nextLine(){
			if(!open){return null;}
			return next<LINES.length ? LINES[next++].getBytes(StandardCharsets.UTF_8) : null;
		}
		/** Status to return when closed. */
		final boolean closeError;
		/** Next input line index and number of close attempts. */
		int next=0, closes=0;
		/** Whether the reader can still return lines. */
		boolean open=true;
	}

	/** Injected FASTA lines read by InjectedByteFile. */
	private static final String[] LINES={">tid|4 read0", "ACGT"};
	/** Expected decompressed output for the injected FASTA. */
	private static final String EXPECTED_FASTA=">tid|4 read0\nACGT\n";
	/** Expected exception message from a fixed close-status failure. */
	private static final String FIXED_CLOSE_MESSAGE=InjectedExplodeTree.class.getName()+
			" terminated in an error state; the output may be corrupt.";
}
