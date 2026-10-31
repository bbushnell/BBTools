package test.jasper;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import jasper.KmerPosition;
import jasper.KmerPosition3;

/** Exact finite baseline fixtures for jasper k-mer position tools.
 * @author Jean */
public final class KmerPositionTest{

	/** This class only contains static test fixtures. */
	private KmerPositionTest(){
		throw new AssertionError();
	}

	/** Runs exact KmerPosition and KmerPosition3 output fixtures.
	 * @param args No arguments are expected
	 * @throws Exception if a temporary file or jasper run fails */
	public static void main(String[] args) throws Exception{
		if(args.length>0){throw new IllegalArgumentException("No arguments expected");}

		final Path dir=Files.createTempDirectory("jasper-kmer-position-");
		try{
			write(dir.resolve("ref.fa"), ">ref\nACGTA\n");
			write(dir.resolve("r1.fq"), "@r1/1\nACGTA\n+\nIIIII\n@r2/1\nTTAC\n+\nIIII\n");
			write(dir.resolve("r2.fq"), "@r1/2\nTACCC\n+\nIIIII\n@r2/2\nNNACG\n+\nIIIII\n");
			write(dir.resolve("ref=equals.fa"), ">ref\nACGTA\n");
			write(dir.resolve("r=1.fq"), "@r1/1\nACGTA\n+\nIIIII\n@r2/1\nTTAC\n+\nIIII\n");
			write(dir.resolve("r=2.fq"), "@r1/2\nTACCC\n+\nIIIII\n@r2/2\nNNACG\n+\nIIIII\n");
			write(dir.resolve("leading-ref.fa"), ">ref\nACG\n");
			write(dir.resolve("leading.fq"), "@leading\nNACG\n+\nIIII\n");

			int failures=0;
			failures+=runCase("kmerposition-forward", new TestCase(){
				@Override
				public void run() throws Exception{checkKmerPosition(dir);}
			});
			failures+=runCase("kmerposition-equals-filenames", new TestCase(){
				@Override
				public void run() throws Exception{checkEqualsFilenames(dir);}
			});
			failures+=runCase("kmerposition3-rcomp-true", new TestCase(){
				@Override
				public void run() throws Exception{
					checkKmerPosition3Rcomp(dir, true, KMER_POSITION3_RCOMP_TRUE);
				}
			});
			failures+=runCase("kmerposition3-rcomp-false", new TestCase(){
				@Override
				public void run() throws Exception{
					checkKmerPosition3Rcomp(dir, false, KMER_POSITION3_RCOMP_FALSE);
				}
			});
			failures+=runCase("kmerposition3-leading-n", new TestCase(){
				@Override
				public void run() throws Exception{checkKmerPosition3LeadingN(dir);}
			});
			if(failures>0){throw new AssertionError("jasper fixture failures="+failures);}
		}finally{
			deleteFixtures(dir);
		}
	}

	/** Checks current string-based exact forward k-mer matching.
	 * @param dir Temporary directory holding input FASTA/FASTQ files
	 * @throws IOException if the output file cannot be read */
	private static void checkKmerPosition(Path dir) throws IOException{
		final Path out=dir.resolve("kmer-position.tsv");
		KmerPosition.main(new String[] {
				"in="+dir.resolve("r1.fq"),
				"in2="+dir.resolve("r2.fq"),
				"ref="+dir.resolve("ref.fa"),
				"out="+out,
				"k=3"});
		expect("kmerposition-forward", slurp(out), KMER_POSITION);
	}

	/** Checks that paths containing '=' remain intact during argument splitting.
	 * @param dir Temporary directory holding input FASTA/FASTQ files
	 * @throws IOException if the output files cannot be read */
	private static void checkEqualsFilenames(Path dir) throws IOException{
		final Path out=dir.resolve("kmer=position.tsv");
		KmerPosition.main(new String[] {
				"in="+dir.resolve("r=1.fq"),
				"in2="+dir.resolve("r=2.fq"),
				"ref="+dir.resolve("ref=equals.fa"),
				"out="+out,
				"k=3"});
		expect("kmerposition-equals-filenames", slurp(out), KMER_POSITION);
		final Path out3=dir.resolve("kmer=position3.tsv");
		KmerPosition3.main(new String[] {
				"in="+dir.resolve("r=1.fq"),
				"in2="+dir.resolve("r=2.fq"),
				"ref="+dir.resolve("ref=equals.fa"),
				"out="+out3,
				"k=3",
				"rcomp=false"});
		expect("kmerposition3-equals-filenames", slurp(out3), KMER_POSITION3_RCOMP_FALSE);
	}

	/** Checks 2-bit exact matching with and without reverse-complemented reference k-mers.
	 * @param dir Temporary directory holding input FASTA/FASTQ files
	 * @param rcomp Whether KmerPosition3 should add reference reverse complements
	 * @param expected Expected tab-delimited output
	 * @throws IOException if the output file cannot be read */
	private static void checkKmerPosition3Rcomp(Path dir, boolean rcomp, String expected)
			throws IOException{
		final Path out=dir.resolve("kmer-position3-"+rcomp+".tsv");
		KmerPosition3.main(new String[] {
				"in="+dir.resolve("r1.fq"),
				"in2="+dir.resolve("r2.fq"),
				"ref="+dir.resolve("ref.fa"),
				"out="+out,
				"k=3",
				"rcomp="+rcomp});
		expect("kmerposition3-rcomp-"+rcomp, slurp(out), expected);
	}

	/** Checks that leading ambiguous bases produce explicit 0 percentages.
	 * @param dir Temporary directory holding input FASTA/FASTQ files
	 * @throws IOException if the output file cannot be read */
	private static void checkKmerPosition3LeadingN(Path dir) throws IOException{
		final Path out=dir.resolve("kmer-position3-leading-n.tsv");
		KmerPosition3.main(new String[] {
				"in="+dir.resolve("leading.fq"),
				"ref="+dir.resolve("leading-ref.fa"),
				"out="+out,
				"k=3",
				"rcomp=false"});
		expect("kmerposition3-leading-n", slurp(out), KMER_POSITION3_LEADING_N);
	}

	/** Runs one fixture case and allows later cases to keep running after a failure.
	 * @param label Case label
	 * @param test Test case to run
	 * @return 1 on failure, 0 on success */
	private static int runCase(String label, TestCase test){
		try{
			test.run();
			return 0;
		}catch(Throwable t){
			System.out.println(label+"\tFAIL\t"+t);
			t.printStackTrace(System.out);
		}
		return 1;
	}

	/** One independently run fixture. */
	private interface TestCase{
		/** Runs this fixture.
		 * @throws Exception if the fixture fails */
		void run() throws Exception;
	}

	/** Writes ASCII fixture text.
	 * @param path Destination path
	 * @param text Text to write
	 * @throws IOException if writing fails */
	private static void write(Path path, String text) throws IOException{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII));
	}

	/** Reads ASCII output text.
	 * @param path Source path
	 * @return File contents
	 * @throws IOException if reading fails */
	private static String slurp(Path path) throws IOException{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII).replace("\r", "");
	}

	/** Compares exact output strings.
	 * @param label Case label
	 * @param actual Observed output
	 * @param expected Expected output */
	private static void expect(String label, String actual, String expected){
		if(expected.equals(actual)){
			System.out.println(label+"\tPASS");
			return;
		}
		throw new AssertionError(label+" expected:\n"+expected+"\nobserved:\n"+actual);
	}

	/** Deletes the finite test's generated files and directory.
	 * @param dir Temporary fixture directory
	 * @throws IOException if deletion fails */
	private static void deleteFixtures(Path dir) throws IOException{
		delete(dir.resolve("kmer-position3-true.tsv"));
		delete(dir.resolve("kmer-position3-false.tsv"));
		delete(dir.resolve("kmer=position3.tsv"));
		delete(dir.resolve("kmer-position3-leading-n.tsv"));
		delete(dir.resolve("kmer-position.tsv"));
		delete(dir.resolve("kmer=position.tsv"));
		delete(dir.resolve("leading.fq"));
		delete(dir.resolve("leading-ref.fa"));
		delete(dir.resolve("r=2.fq"));
		delete(dir.resolve("r2.fq"));
		delete(dir.resolve("r=1.fq"));
		delete(dir.resolve("r1.fq"));
		delete(dir.resolve("ref=equals.fa"));
		delete(dir.resolve("ref.fa"));
		delete(dir);
	}

	/** Deletes one generated path if it exists.
	 * @param path Path to delete
	 * @throws IOException if deletion fails */
	private static void delete(Path path) throws IOException{
		Files.deleteIfExists(path);
	}

	/** Shared output header for both jasper tools. */
	private static final String HEADER="#pos\tread1_count\tread1_perc\tread2_count\tread2_perc\n";

	/** Current forward-only exact output from KmerPosition for the finite pair fixture. */
	private static final String KMER_POSITION=HEADER+
			"0\t1\t50\t0\t0\n"+
			"1\t1\t50\t0\t0\n"+
			"2\t1\t100\t1\t50\n";

	/** Current KmerPosition3 output when reference reverse complements are included. */
	private static final String KMER_POSITION3_RCOMP_TRUE=HEADER+
			"0\t1\t50\t1\t100\n"+
			"1\t2\t100\t0\t0\n"+
			"2\t1\t100\t1\t50\n";

	/** Current KmerPosition3 output without reference reverse-complement k-mers. */
	private static final String KMER_POSITION3_RCOMP_FALSE=HEADER+
			"0\t1\t50\t0\t0\n"+
			"1\t1\t50\t0\t0\n"+
			"2\t1\t100\t1\t50\n";

	/** KmerPosition3 output for a leading-N read whose first valid k-mer starts at 1. */
	private static final String KMER_POSITION3_LEADING_N=HEADER+
			"0\t0\t0\t0\t0\n"+
			"1\t1\t100\t0\t0\n";

}
