package test.dna;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import dna.ChromosomeArray;
import dna.ChromToFasta;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Shared;

/** Finite line-width checks using tiny inputs and literal FASTA expectations.
 * Run each mode in a fresh JVM. Borrowed-writer probes stop after eight delegated
 * print calls and always close their caller-owned writer. Owned-writer probes use
 * only empty input, so the faulty baseline never enters an unbounded loop.
 * This fixture does not test arbitrary writer failure or large-coordinate overflow.
 * @author Sucrose
 */
public final class ChromToFastaWidthTest{

	/** Requires case=controls|borrowed-invalid|owned-invalid and a new work= directory. */
	public static void main(String[] args) throws Exception{
		String mode=null;
		Path work=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("work=")){work=Paths.get(arg.substring(5));}
			else{throw new IllegalArgumentException(arg);}
		}
		if(work==null || !("controls".equals(mode) || "borrowed-invalid".equals(mode) || "owned-invalid".equals(mode))){
			throw new IllegalArgumentException("Require case=controls|borrowed-invalid|owned-invalid and work=");
		}
		Files.createDirectory(work);
		Shared.setThreads(1);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		if("controls".equals(mode)){controls(work);}
		else if("borrowed-invalid".equals(mode)){borrowedInvalid(work);}
		else{ownedInvalid(work);}
		System.out.println("PASS "+mode);
	}

	/** Literal controls cover wrapping, empty input, split numbering and both ownership forms. */
	private static void controls(Path work) throws Exception{
		boolean pass=true;
		pass&=borrowed(work, "contig", "ACGT", 2, "AC\nGT\n", -1);
		pass&=borrowed(work, "contig", "", 4, "", -1);
		pass&=borrowed(work, "chrom", "NNACNNTGNN", 4, ">7\nNNAC\nNNTG\nNN\n", -1);
		pass&=borrowed(work, "chrom", "", 4, ">7\n", -1);
		pass&=borrowed(work, "split", "NNACNNTGNN", 4, ">1\nAC\n>2\nTG\n", 3);
		pass&=borrowed(work, "split", "", 4, "", 1);
		pass&=borrowed(work, "split", "NNNN", 4, "", 1);
		pass&=owned(work, "filename", 4, ">7\n");
		pass&=owned(work, "cli-whole", 4, "");
		pass&=owned(work, "cli-split", 4, "");
		pass&=owned(work, "cli-options-whole", 4, "");
		pass&=owned(work, "cli-options-split", 4, "");
		pass&=owned(work, "cli-config-whole", 4, "");
		pass&=owned(work, "cli-config-split", 4, "");
		require(pass, "Positive-width FASTA or ownership control failed");
	}

	/** Every borrowed entry must reject before printing, even when no sequence would be emitted. */
	private static void borrowedInvalid(Path work) throws Exception{
		boolean pass=true;
		for(int width : new int[]{0, -1}){
			pass&=borrowed(work, "contig", "ACGT", width, null, -1);
			pass&=borrowed(work, "contig", "", width, null, -1);
			pass&=borrowed(work, "chrom", "ACGT", width, null, -1);
			pass&=borrowed(work, "chrom", "", width, null, -1);
			pass&=borrowed(work, "split", "NNACNNTGNN", width, null, -1);
			pass&=borrowed(work, "split", "", width, null, -1);
			pass&=borrowed(work, "split", "NNNN", width, null, -1);
		}
		require(pass, "Nonpositive borrowed width must throw a diagnostic IAE before output");
	}

	/** Empty input keeps all owned-writer baseline calls finite; sentinel bytes expose early truncation. */
	private static void ownedInvalid(Path work) throws Exception{
		boolean pass=true;
		for(int width : new int[]{0, -1}){
			pass&=owned(work, "filename", width, null);
			pass&=owned(work, "cli-whole", width, null);
			pass&=owned(work, "cli-split", width, null);
			pass&=owned(work, "cli-options-whole", width, null);
			pass&=owned(work, "cli-options-split", width, null);
			pass&=owned(work, "cli-config-whole", width, null);
			pass&=owned(work, "cli-config-split", width, null);
		}
		require(pass, "Nonpositive owned width must throw a diagnostic IAE and preserve existing output");
	}

	/** Calls a real formatter with a capped, started writer and closes it even after a probe exception.
	 * A null expected output selects the invalid-width oracle; otherwise output and split return are literal.
	 */
	private static boolean borrowed(Path work, String kind, String sequence, int width, String expected, int expectedNext) throws Exception{
		final String label=kind+"-length"+sequence.length()+"-width"+width;
		final Path output=work.resolve(label+".fa");
		final StringBuilder builder=new StringBuilder(sequence);
		final ChromosomeArray chromosome=chromosome(sequence);
		final CappedWriter writer=new CappedWriter(output);
		RuntimeException failure=null;
		int next=-1;
		boolean leftOpen=false;
		writer.start();
		try{
			try{
				if("contig".equals(kind)){ChromToFasta.writeContig(builder, writer, width);}
				else if("chrom".equals(kind)){ChromToFasta.writeChrom(chromosome, writer, width);}
				else{next=ChromToFasta.writeContigs(chromosome, 1, 2, width, writer);}
			}catch(RuntimeException problem){failure=problem;}
			leftOpen=writer.isAlive();
			writer.appendCallerMarker();
		}finally{
			final boolean writerError=writer.poisonAndWait();
			require(!writer.isAlive() && writer.getState()==Thread.State.TERMINATED,
					"Caller cleanup must terminate the writer: "+label);
			require(!writerError, "Unexpected raw-file writer failure: "+label);
		}
		final String written=read(output);
		require(written.endsWith(CALLER_MARKER), "Borrowing formatter must leave writer usable by its caller: "+label);
		final String actual=written.substring(0, written.length()-CALLER_MARKER.length());
		final boolean preserved=builder.toString().equals(sequence);
		final boolean pass=leftOpen && preserved && (expected==null ?
				diagnosticIAE(failure) && writer.calls==0 && actual.isEmpty() :
				failure==null && expected.equals(actual) && next==expectedNext);
		System.out.println(label+" exception="+failureName(failure)+" printCalls="+writer.calls+
				" leftOpen="+leftOpen+" builderUnchanged="+preserved+" terminated=true writerError=false"+
				" next="+next+" output="+visible(actual)+" pass="+pass);
		return pass;
	}

	/** Tests owned output only with empty input; these calls cannot enter either formatting loop. */
	private static boolean owned(Path work, String kind, int width, String expected) throws Exception{
		final String label=kind+"-width"+width;
		final Path output=work.resolve(label+".fa");
		Files.write(output, SENTINEL.getBytes(StandardCharsets.US_ASCII));
		RuntimeException failure=null;
		if("filename".equals(kind)){
			final ChromosomeArray empty=chromosome("");
			require(empty.maxIndex==-1, "Owned probe must not enter the sequence loop");
			try{ChromToFasta.writeChrom(empty, output.toString(), width);}
			catch(RuntimeException problem){failure=problem;}
		}else{
			final Path input=work.resolve(label+"-empty-input");
			Files.createDirectory(input);
			final String[] arguments=cliArguments(work, label, kind, input, output, width);
			try{ChromToFasta.main(arguments);}
			catch(RuntimeException problem){failure=problem;}
		}
		final String actual=read(output);
		final boolean pass=(expected==null ? diagnosticIAE(failure) && SENTINEL.equals(actual) :
				failure==null && expected.equals(actual));
		System.out.println(label+" exception="+failureName(failure)+" sentinelPreserved="+SENTINEL.equals(actual)+
				" output="+visible(actual)+" pass="+pass);
		return pass;
	}

	/** Keeps the requested signed width in positional, option-prefixed or config-expanded input.
	 * The heap flag is passed to main as a removable token; it does not change this JVM's heap.
	 */
	private static String[] cliArguments(Path work, String label, String kind, Path input, Path output, int width) throws Exception{
		final String trigger=kind.endsWith("split") ? "2" : "0";
		final String widthText=Integer.toString(width);
		if(kind.startsWith("cli-options-")){
			return new String[]{"--outstream=stderr", "-Xmx512m", input.toString(), output.toString(), widthText, trigger};
		}
		if(kind.startsWith("cli-config-")){
			final Path config=work.resolve(label+".config");
			final String contents=input+"\n"+output+"\n"+widthText+"\n"+trigger+"\n";
			Files.write(config, contents.getBytes(StandardCharsets.UTF_8));
			return new String[]{"config="+config};
		}
		return new String[]{input.toString(), output.toString(), widthText, trigger};
	}

	/** Builds only tiny arrays, preserving the standard empty sentinel when sequence is empty. */
	private static ChromosomeArray chromosome(String sequence){
		final ChromosomeArray result=new ChromosomeArray(7, Shared.PLUS);
		if(!sequence.isEmpty()){result.set(0, sequence);}
		return result;
	}

	/** Requires an explicit diagnostic rather than accepting unrelated bounds or probe failures. */
	private static boolean diagnosticIAE(RuntimeException failure){
		return failure instanceof IllegalArgumentException && failure.getMessage()!=null && !failure.getMessage().isEmpty();
	}

	private static String failureName(RuntimeException failure){return failure==null ? "none" : failure.getClass().getSimpleName();}
	private static String visible(String text){return text.replace("\n", "|");}
	private static String read(Path path) throws Exception{return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);}

	/** Keeps all fixture oracles active with assertions disabled. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Bounds formatter progress on the producer thread; inherited shutdown bypasses print. */
	private static final class CappedWriter extends TextStreamWriter{
		CappedWriter(Path output){super(output.toString(), true, false, false);}
		@Override
		public void print(CharSequence text){
			calls++;
			if(calls>8){throw new ProbeStop();}
			super.print(text);
		}
		/** One finite caller write after formatting; excluded from the formatter's eight-call budget. */
		void appendCallerMarker(){super.print(CALLER_MARKER);}
		/** Includes the ninth, rejected call; zero proves validation preceded all text output. */
		int calls;
	}

	/** Distinct from IllegalArgumentException, so the finite cutoff can never satisfy the repair oracle. */
	private static final class ProbeStop extends RuntimeException{
		private static final long serialVersionUID=1L;
		ProbeStop(){super("Formatter exceeded eight delegated print calls");}
	}

	/** Existing bytes which invalid owned calls must not truncate or replace. */
	private static final String SENTINEL="KEEP THIS OUTPUT\n";
	/** Appended only by the caller, proving the borrowed writer remains usable. */
	private static final String CALLER_MARKER="CALLER STILL OWNS WRITER\n";
}
