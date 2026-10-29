package test.tax;

import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.Timer;
import tax.TaxSize;
import tax.TaxTree;
import tax.TaxTreeText;

/** Finite injected-reader check for TaxSize input close-status propagation.
 * @author Brian Bushnell, Sucrose */
public final class TaxSizeCompletionErrorTest{

	/** Runs an ordinary success control and an injected close-failure case.
	 * @param args Optional expect=baseline or expect=fixed, plus optional
	 * reader=baseline or reader=fixed for the throwing-reader oracle
	 * @throws Exception If fixture setup or processing fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		Boolean readerBaseline=null;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else if(arg.equals("reader=baseline")){readerBaseline=Boolean.TRUE;}
			else if(arg.equals("reader=fixed")){readerBaseline=Boolean.FALSE;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		final boolean bf1=ByteFile.FORCE_MODE_BF1, bf2=ByteFile.FORCE_MODE_BF2;
		try{
			run(false, baseline);
			run(true, baseline);
			if(readerBaseline!=null){runThrowingReader(readerBaseline.booleanValue());}
		}finally{
			ByteFile.FORCE_MODE_BF1=bf1;
			ByteFile.FORCE_MODE_BF2=bf2;
		}
		System.out.println("TAXSIZE_COMPLETION_ERROR_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private TaxSizeCompletionErrorTest(){}

	/** Runs one real TaxSize process over an injected finite ByteFile.
	 * @param closeError Whether the injected input reports close failure
	 * @param baseline Whether to assert baseline or fixed expectations
	 * @throws Exception If the temporary tree, process, or output read fails */
	private static void run(final boolean closeError, final boolean baseline) throws Exception{
		final Path root=Files.createTempDirectory("taxsize-completion-");
		try{
			final Path treeFile=root.resolve("tree.tsv");
			final Path input=root.resolve("input.fa");
			final Path output=root.resolve("taxsize.tsv");
			TaxTreeText.write(TaxTreeContractTest.tree(TaxTree.GENUS_E), treeFile.toString(), true);
			Files.write(input, new byte[0]);

			final ByteFile reader=new InjectedByteFile(FileFormat.testInput(input.toString(),
					FileFormat.FASTA, null, true, true), closeError);
			final InjectedTaxSize taxSize=new InjectedTaxSize(input, output, treeFile, reader);
			Throwable failure=null;
			try{taxSize.process(new Timer());}
			catch(Throwable t){failure=t;}
			check(((InjectedByteFile)reader).closes==1, "TaxSize must close its input exactly once");
			check(new String(Files.readAllBytes(output), StandardCharsets.UTF_8).equals(EXPECTED_OUTPUT),
					"TaxSize must emit the exact healthy table before any input close-status failure");

			if(!closeError){
				check(failure==null, "Ordinary successful input close must not fail: "+failure);
				check(!taxSize.errorState, "Ordinary successful input close must not set errorState");
			}else if(baseline){
				check(failure==null, "Baseline must expose the discarded close-status failure");
				check(!taxSize.errorState, "Baseline must leave errorState false after a close failure");
			}else{
				check(failure instanceof RuntimeException,
						"Fixed TaxSize must fail loudly on input close failure: "+failure);
				check(FIXED_CLOSE_MESSAGE.equals(failure.getMessage()),
						"Fixed TaxSize must reuse the existing errorState exception text: "+failure);
				check(taxSize.errorState, "Fixed TaxSize must retain close failure in errorState");
			}
		}finally{
			deleteTree(root);
		}
	}

	/** Runs one TaxSize process whose input throws before EOF.
	 * @param baseline Whether to assert baseline or fixed expectations
	 * @throws Exception If the temporary tree or process setup fails */
	private static void runThrowingReader(final boolean baseline) throws Exception{
		final Path root=Files.createTempDirectory("taxsize-throwing-");
		try{
			final Path treeFile=root.resolve("tree.tsv");
			final Path input=root.resolve("input.fa");
			final Path output=root.resolve("taxsize.tsv");
			TaxTreeText.write(TaxTreeContractTest.tree(TaxTree.GENUS_E), treeFile.toString(), true);
			Files.write(input, new byte[0]);

			final ByteFile reader=new ThrowingByteFile(FileFormat.testInput(
					input.toString(), FileFormat.FASTA, null, true, true));
			final InjectedTaxSize taxSize=new InjectedTaxSize(input, output, treeFile, reader);
			Throwable failure=null;
			try{taxSize.process(new Timer());}
			catch(Throwable t){failure=t;}

			check(failure==NEXT_LINE_FAILURE,
					"TaxSize must retain the exact primary nextLine failure: "+failure);
			final int expectedCloses=baseline ? 0 : 1;
			final int closes=((InjectedByteFile)reader).closes;
			check(closes==expectedCloses,
					"Throwing input close count expected "+expectedCloses+
					", observed "+closes);
			System.out.println("TAXSIZE_THROWING_READER_RECEIPT expect="+
					(baseline ? "baseline" : "fixed")+" primary=exact closes="+
					closes);
		}finally{
			deleteTree(root);
		}
	}

	/** Deletes the owned temporary tree, visiting children before parents.
	 * @param root owned temporary directory root
	 * @throws Exception if the tree cannot be walked or deleted */
	private static void deleteTree(final Path root) throws Exception{
		try(java.util.stream.Stream<Path> paths=Files.walk(root)){
			for(Path path : (Iterable<Path>)paths.sorted(java.util.Comparator.reverseOrder())::iterator){
				Files.deleteIfExists(path);
			}
		}
	}

	/** Test checks stay active under -da.
	 * @param condition condition required for test success
	 * @param message diagnostic message for failure */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** TaxSize variant that returns the injected finite reader instead of opening input. */
	private static final class InjectedTaxSize extends TaxSize{
		/** Constructs a TaxSize around real descriptors and a later injected reader.
		 * @param input FASTA input placeholder used for descriptor validation
		 * @param output TaxSize TSV output path
		 * @param tree TaxTree file path
		 * @param reader_ injected finite input reader */
		InjectedTaxSize(final Path input, final Path output, final Path tree, final ByteFile reader_){
			super(new String[] {"in="+input, "out="+output, "tree="+tree});
			reader=reader_;
		}
		@Override
		protected ByteFile makeByteFile(final FileFormat ff){return reader;}
		/** Injected finite input reader. */
		private final ByteFile reader;
	}

	/** Minimal ByteFile that supplies one FASTA record and then reports a chosen close status. */
	private static class InjectedByteFile extends ByteFile{
		/** Creates a finite reader over one synthetic FASTA record.
		 * @param ff validated FASTA input descriptor
		 * @param closeError_ whether close reports an error */
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
		public void pushBack(final byte[] line){throw new AssertionError("TaxSize does not use pushBack");}
		@Override
		public boolean isOpen(){return open;}
		@Override
		public byte[] nextLine(){
			if(!open){return null;}
			return next<LINES.length ? LINES[next++].getBytes(StandardCharsets.UTF_8) : null;
		}
		/** Whether close reports an error. */
		private final boolean closeError;
		/** Next synthetic line index and number of close calls. */
		private int next=0, closes=0;
		/** Whether the reader remains open. */
		private boolean open=true;
	}

	/** ByteFile that throws the exact injected primary exception from nextLine. */
	private static final class ThrowingByteFile extends InjectedByteFile{
		/** Creates a throwing reader with successful close status.
		 * @param ff validated FASTA input descriptor */
		ThrowingByteFile(final FileFormat ff){super(ff, false);}
		@Override
		public byte[] nextLine(){throw NEXT_LINE_FAILURE;}
	}

	/** One assigned FASTA record with four counted nonheader bytes. */
	private static final String[] LINES={">tid|4 read0", "ACGT"};

	/** Primary unchecked input failure thrown by the throwing-reader case. */
	private static final RuntimeException NEXT_LINE_FAILURE=
			new RuntimeException("injected nextLine failure");

	/** Exact TaxSize table for one four-base record at taxID 4 in the five-node tree. */
	private static final String EXPECTED_OUTPUT="#taxID\tbases\tbasesC\tseqs\tseqsC\tnodesC\n"
			+"1\t0\t4\t0\t1\t5\n"
			+"2\t0\t4\t0\t1\t4\n"
			+"3\t0\t4\t0\t1\t3\n"
			+"4\t4\t4\t1\t1\t1\n"
			+"5\t0\t0\t0\t0\t1\n";

	/** Existing process-level errorState exception text. */
	private static final String FIXED_CLOSE_MESSAGE=InjectedTaxSize.class.getName()+
			" terminated in an error state; the output may be corrupt.";
}
