package test.assemble;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;

import assemble.AbstractExploreThread;
import assemble.AbstractShaveThread;
import assemble.ShaveObject;
import assemble.Shaver;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/**
 * Synthetic Shaver diagnostic aggregation witness for event codes outside 0..7.
 * Injects worker matrices directly; this is not a production graph reproduction.
 * @author Jean
 */
public final class ShaverEventCountsTest{

	/** Static test entry point only. */
	private ShaverEventCountsTest(){}

	/** Checks that worker branch cells appear in diagnostic matrices. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, ShaverEventCountsTest.class, false);
		try{
			Path root=null;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				if(pair.length<2){throw new IllegalArgumentException("Expected flag=value for "+arg);}
				if(pair[0].equalsIgnoreCase("outdir")){root=Paths.get(pair[1]).toAbsolutePath();}
				else{throw new IllegalArgumentException("Unknown test option: "+arg);}
			}
			if(root==null){throw new IllegalArgumentException("Specify new outdir=DIR.");}
			Files.createDirectory(root);

			final Path input=root.resolve("reads.fa");
			Files.write(input, FASTA.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
			final KmerTableSetU table=tables(input);
			try{
				final StringBuilder checks=new StringBuilder("case\tpassed\tdetail\n");
				final CaptureResult result=captureShave(table);
				Files.write(root.resolve("event-output.txt"),
						result.output.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
				final long[][] eventMatrix=matrix(result.output, "Event counts:");
				final long[][] removalMatrix=matrix(result.output, "Removal counts:");
				int failed=0;
				failed+=check(checks, "removed-count", result.removed==0, "removed="+result.removed);
				failed+=check(checks, "worker-bubble-total",
						result.shaver.bubblesFound==1, "bubbles="+result.shaver.bubblesFound);
				failed+=check(checks, "event-matrix", matchesExpectedEventMatrix(eventMatrix),
						matrixDetail(eventMatrix));
				failed+=check(checks, "removal-matrix", matchesExpectedRemovalMatrix(removalMatrix),
						matrixDetail(removalMatrix));
				Files.write(root.resolve("checks.tsv"),
						checks.toString().getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
				System.out.println("SHAVER_EVENT_COUNTS_TEST checks=4 failed="+failed);
				if(failed>0){throw new AssertionError("Shaver event-count diagnostic contract failures; inspect "+root);}
			}finally{table.clear();}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Creates a tiny real table set; synthetic workers only need ownership setup. */
	private static KmerTableSetU tables(final Path input){
		final String[] options={"in="+input, "k=63", "packed=t", "rcomp=t", "ways=1",
				"initialsize=128", "prealloc=f", "prefilter=f", "minprob=0", "minprobmain=f",
				"interleaved=f", "t=1"};
		final KmerTableSetU table=new KmerTableSetU(options, 0);
		try{
			table.process(new Timer());
		}catch(RuntimeException t){
			table.clear();
			throw t;
		}catch(Error t){
			table.clear();
			throw t;
		}
		return table;
	}

	/** Captures printEventCounts while restoring every static diagnostic sink. */
	private static CaptureResult captureShave(final KmerTableSetU table) throws Exception{
		final boolean oldPrint=ShaveObject.printEventCounts;
		final PrintStream oldStream=getOutstream();
		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final PrintStream capture=new PrintStream(buffer, true, "US-ASCII");
		try{
			ShaveObject.printEventCounts=true;
			setOutstream(capture);
			final MatrixShaver shaver=new MatrixShaver(table);
			final long removed=shaver.shave();
			capture.flush();
			final String output=new String(buffer.toByteArray(), StandardCharsets.US_ASCII);
			return new CaptureResult(removed, shaver, output);
		}finally{
			setOutstream(oldStream);
			ShaveObject.printEventCounts=oldPrint;
			capture.close();
		}
	}

	/** Parses the exact matrix printed after a diagnostic heading. */
	private static long[][] matrix(final String output, final String heading){
		final int start=output.indexOf(heading);
		if(start<0){throw new AssertionError("Missing heading: "+heading+"\n"+output);}
		final String[] lines=output.substring(start+heading.length()).split("\\n");
		final ArrayList<long[]> rows=new ArrayList<long[]>();
		for(String line : lines){
			line=line.trim();
			if(line.length()<1){
				if(rows.size()>0){break;}
				continue;
			}
			final String[] split=line.split("\\s+");
			final long[] row=new long[split.length];
			for(int col=0; col<split.length; col++){
				try{
					row[col]=Long.parseLong(split[col]);
				}catch(NumberFormatException e){
					if(rows.size()>0){return rows.toArray(new long[rows.size()][]);}
					throw new AssertionError("Expected matrix row after "+heading+": "+line);
				}
			}
			rows.add(row);
		}
		return rows.toArray(new long[rows.size()][]);
	}

	/** Correct diagnostic output includes every injected synthetic worker event. */
	private static boolean matchesExpectedEventMatrix(final long[][] matrix){
		final long[][] expected=expectedMatrix();
		expected[ShaveObject.F_BRANCH][ShaveObject.D_BRANCH]=1;
		return matches(matrix, expected);
	}

	/** Correct diagnostic output includes every injected synthetic removal event. */
	private static boolean matchesExpectedRemovalMatrix(final long[][] matrix){
		return matches(matrix, expectedMatrix());
	}

	/** Returns a 20x20 matrix containing all removable synthetic worker cells. */
	private static long[][] expectedMatrix(){
		final long[][] expected=new long[ShaveObject.MAX_CODE][ShaveObject.MAX_CODE];
		expected[ShaveObject.DEAD_END][ShaveObject.DEAD_END]=1;
		expected[ShaveObject.DEAD_END][ShaveObject.B_BRANCH]=1;
		expected[ShaveObject.B_BRANCH][ShaveObject.B_BRANCH]=1;
		return expected;
	}

	/** Returns true only for identical matrix dimensions and contents. */
	private static boolean matches(final long[][] matrix, final long[][] expected){
		if(matrix.length!=expected.length){return false;}
		for(int row=0; row<matrix.length; row++){
			if(matrix[row].length!=expected[row].length){return false;}
			for(int col=0; col<matrix[row].length; col++){
				if(matrix[row][col]!=expected[row][col]){return false;}
			}
		}
		return true;
	}

	/** Keeps every expectation active without depending on JVM assertions. */
	private static int check(final StringBuilder checks, final String name,
			final boolean passed, final String detail){
		checks.append(name).append('\t').append(passed ? 1 : 0).append('\t').append(detail).append('\n');
		return passed ? 0 : 1;
	}

	/** Compact one-line detail that preserves every printed matrix row. */
	private static String matrixDetail(final long[][] matrix){
		final StringBuilder sb=new StringBuilder();
		for(int row=0; row<matrix.length; row++){
			if(row>0){sb.append(';');}
			for(int col=0; col<matrix[row].length; col++){
				if(col>0){sb.append(',');}
				sb.append(matrix[row][col]);
			}
		}
		return sb.toString();
	}

	/** Returns a named long matrix from an exploration worker. */
	private static long[][] workerMatrix(final AbstractExploreThread worker, final String name){
		try{
			final Field field=AbstractExploreThread.class.getDeclaredField(name);
			field.setAccessible(true);
			return (long[][])field.get(worker);
		}catch(ReflectiveOperationException e){throw new AssertionError(e);}
	}

	/** Reads ShaveObject's package-local diagnostic stream. */
	private static PrintStream getOutstream() throws Exception{
		final Field field=ShaveObject.class.getDeclaredField("outstream");
		field.setAccessible(true);
		return (PrintStream)field.get(null);
	}

	/** Sets ShaveObject's package-local diagnostic stream. */
	private static void setOutstream(final PrintStream outstream) throws Exception{
		final Field field=ShaveObject.class.getDeclaredField("outstream");
		field.setAccessible(true);
		field.set(null, outstream);
	}

	/** Test-only Shaver whose exploration worker injects synthetic matrix cells. */
	private static final class MatrixShaver extends Shaver{
		MatrixShaver(final KmerTableSetU tables_){
			super(tables_, 1, 1, 1, 1, 1, 3, 100, 100, true, true);
			tables=tables_;
		}

		@Override
		protected AbstractExploreThread makeExploreThread(final int id){
			return new AbstractExploreThread(id, tables.kbig()){
				@Override
				protected boolean processNextTable(final Kmer kmer, final Kmer temp){
					final long[][] countMatrix=workerMatrix(this, "countMatrixT");
					final long[][] removeMatrix=workerMatrix(this, "removeMatrixT");
					countMatrix[DEAD_END][DEAD_END]=1;
					removeMatrix[DEAD_END][DEAD_END]=1;
					countMatrix[DEAD_END][B_BRANCH]=1;
					removeMatrix[DEAD_END][B_BRANCH]=1;
					countMatrix[B_BRANCH][B_BRANCH]=1;
					removeMatrix[B_BRANCH][B_BRANCH]=1;
					countMatrix[F_BRANCH][D_BRANCH]=1;
					return false;
				}

				@Override
				protected boolean processNextVictims(final Kmer kmer, final Kmer temp){return false;}
			};
		}

		@Override
		protected AbstractShaveThread makeShaveThread(final int id){
			return new AbstractShaveThread(id){
				@Override
				protected boolean processNextTable(){return false;}
			};
		}

		@Override
		protected KmerTableSetU tables(){return tables;}

		private final KmerTableSetU tables;
	}

	/** Captured result of one synthetic shave pass. */
	private static final class CaptureResult{
		CaptureResult(final long removed_, final MatrixShaver shaver_, final String output_){
			removed=removed_;
			shaver=shaver_;
			output=output_;
		}

		final long removed;
		final MatrixShaver shaver;
		final String output;
	}

	private static final String FASTA=">a\n"
			+"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
}
