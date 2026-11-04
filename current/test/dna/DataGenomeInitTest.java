package test.dna;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.zip.GZIPOutputStream;

import dna.Data;
import dna.FastaToChromArrays2;
import fileIO.ReadWrite;
import shared.Shared;

/** Checks completion-marker semantics after catchable genome metadata failures.
 * Each mode requires a fresh JVM and work directory. Disk metadata describes one
 * four-base chromosome; no chromosome objects or mapping indexes are generated.
 * Memory retries replace all three consumed lists. This does not test rollback,
 * concurrent readers, automatic production retries or exceptional reader closure.
 * @author Sucrose
 */
public final class DataGenomeInitTest{

	/** Runs case=disk|lists|controls work=absolute-new-directory. */
	public static void main(String[] args) throws Exception{
		String mode=null;
		Path work=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("work=")){work=Paths.get(arg.substring(5));}
			else{throw new IllegalArgumentException(arg);}
		}
		require(("disk".equals(mode) || "lists".equals(mode) || "controls".equals(mode)) &&
				work!=null && work.isAbsolute(), "Require case=disk|lists|controls and absolute work path");
		Files.createDirectory(work);
		require(Data.GENOME_BUILD==-1 && FastaToChromArrays2.LISTBUILD==-1 &&
				FastaToChromArrays2.SUMMARY_LIST==null && FastaToChromArrays2.INFO_LIST==null &&
				FastaToChromArrays2.SCAF_LIST==null, "Fixture requires fresh JVM reference state");
		require(FastaToChromArrays2.currentVersion()==5, "Literal metadata fixture requires format version 5");
		Shared.setThreads(1);
		ReadWrite.RAWMODE=false;
		ReadWrite.USE_UNPIGZ=ReadWrite.USE_GUNZIP=ReadWrite.USE_UNBGZIP=false;
		ReadWrite.ALLOW_NATIVE_BGZF=false;
		Data.setPath(work.toString());
		Data.LOAD_SCAFFOLDS=true;
		Data.allowZeroSizedGenome=false;
		final boolean lists="lists".equals(mode);
		if("controls".equals(mode)){
			install(work, false, false);
			Data.setGenome(BUILD);
			require(complete(), "Ordinary metadata load differs from literal fixture");
			checkGuards();
		}else{
			install(work, lists, true);
			final RuntimeException first=attempt();
			// All lists are fresh again: failed parsing consumed summary and part of info.
			install(work, lists, true);
			final RuntimeException repeated=attempt();
			install(work, lists, false);
			final RuntimeException recovered=attempt();
			final boolean loaded=recovered==null && complete();
			System.out.println("case="+mode+" first="+describe(first)+" repeated="+
					describe(repeated)+" recovery="+describe(recovered)+" complete="+loaded);
			// Record both retry outcomes before failing the baseline for silent returns.
			require(parseFailure(first) && parseFailure(repeated) && loaded,
					"Failed initialization must reject bad metadata again and load corrected metadata on retry");
			checkGuards();
		}
		if(lists){require(!Files.exists(work.resolve("ref")), "Memory metadata test wrote a reference tree");}
		System.out.println("PASS "+mode);
	}

	/** Installs complete independent inputs; setup failures cannot count as loader rejection. */
	private static void install(Path work, boolean lists, boolean bad) throws Exception{
		final String info=bad ? BAD_INFO : GOOD_INFO;
		if(lists){
			FastaToChromArrays2.LISTBUILD=BUILD;
			FastaToChromArrays2.SUMMARY_LIST=lines(SUMMARY);
			FastaToChromArrays2.INFO_LIST=lines(info);
			FastaToChromArrays2.SCAF_LIST=lines(SCAFFOLDS);
		}else{
			final Path genome=work.resolve("ref/genome/"+BUILD);
			Files.createDirectories(genome);
			Files.write(genome.resolve("summary.txt"), SUMMARY.getBytes(StandardCharsets.US_ASCII));
			Files.write(genome.resolve("info.txt"), info.getBytes(StandardCharsets.US_ASCII));
			try(GZIPOutputStream out=new GZIPOutputStream(Files.newOutputStream(genome.resolve("scaffolds.txt.gz")))){
				out.write(SCAFFOLDS.getBytes(StandardCharsets.US_ASCII));
			}
		}
	}

	/** Creates new mutable list storage without trailing empty lines. */
	private static ArrayList<String> lines(String text){
		return new ArrayList<String>(Arrays.asList(text.split("\n")));
	}

	/** Catches only loader exceptions; assertion errors and fixture setup failures remain failures. */
	private static RuntimeException attempt(){
		try{Data.setGenome(BUILD);}
		catch(RuntimeException problem){return problem;}
		return null;
	}

	/** Identifies the deliberately malformed numeric field rather than accepting unrelated failures. */
	private static boolean parseFailure(RuntimeException problem){
		return problem instanceof NumberFormatException && problem.getMessage()!=null &&
				problem.getMessage().contains("bad-length");
	}

	/** Formats the observed exception without depending on stack-trace line numbers. */
	private static String describe(RuntimeException problem){
		return problem==null ? "none" : problem.getClass().getSimpleName()+":"+problem.getMessage();
	}

	/** Checks literal public metadata, guarding shapes before indexing partially loaded arrays. */
	private static boolean complete(){
		return Data.GENOME_BUILD==BUILD && Data.numChroms==1 && Data.numBases==4 &&
				Data.numDefinedBases==4 && Data.numContigs==1 && Data.numScaffolds==1 &&
				Data.interScaffoldPadding==0 && "dna032".equals(Data.name) &&
				"fixture".equals(Data.genomeSource) && !Data.scaffoldPrefixes &&
				value(Data.chromLengths, 4) && value(Data.chromDefinedBases, 4) &&
				value(Data.chromUndefinedBases, 0) && value(Data.chromContigs, 1) &&
				value(Data.chromScaffolds, 1) && value(Data.chromStartPad, 0) &&
				Data.chromosomePlusMatrix!=null && Data.chromosomePlusMatrix.length==2 &&
				Data.chromosomePlusMatrix[1]==null &&
				Data.scaffoldNames!=null && Data.scaffoldNames.length==2 &&
				Data.scaffoldNames[1]!=null && Data.scaffoldNames[1].length==1 &&
				Arrays.equals(Data.scaffoldNames[1][0], new byte[]{'o', 'n', 'e'}) &&
				Data.scaffoldLocs!=null && Data.scaffoldLocs.length==2 &&
				Arrays.equals(Data.scaffoldLocs[1], new int[]{0}) &&
				Data.scaffoldLengths!=null && Data.scaffoldLengths.length==2 &&
				Arrays.equals(Data.scaffoldLengths[1], new int[]{4});
	}

	/** Checks the chromosome-one slot and the reserved chromosome-zero array shape. */
	private static boolean value(int[] array, int expected){
		return array!=null && array.length==2 && array[1]==expected;
	}

	/** Successful repeats are no-ops; another build must fail before any metadata I/O. */
	private static void checkGuards(){
		final Object[] before=identities();
		Data.setGenome(BUILD);
		require(sameIdentities(before) && complete(), "Successful same-build initialization must be a no-op");
		RuntimeException changed=null;
		try{Data.setGenome(BUILD+1);}
		catch(RuntimeException problem){changed=problem;}
		require(changed!=null && "Changing genomes is not currently supported.".equals(changed.getMessage()),
				"Different-build request must retain the existing guard diagnostic");
		require(sameIdentities(before) && complete(), "Rejected build change altered initialized metadata");
	}

	/** Captures references, not contents, to detect unwanted reinitialization. */
	private static Object[] identities(){
		return new Object[]{Data.chromosomePlusMatrix, Data.chromLengths, Data.chromDefinedBases,
				Data.chromUndefinedBases, Data.chromContigs, Data.chromScaffolds, Data.chromStartPad,
				Data.scaffoldNames, Data.scaffoldLocs, Data.scaffoldLengths};
	}

	/** Compares the exact metadata arrays captured before a no-op or rejected build change. */
	private static boolean sameIdentities(Object[] before){
		final Object[] after=identities();
		for(int i=0; i<before.length; i++){
			if(before[i]!=after[i]){return false;}
		}
		return true;
	}

	/** Always-on test assertion so -da exercises the same independent oracle. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Positive build identifier scoped to the fixture's private work directory. */
	private static final int BUILD=17;
	/** Literal version-5 totals for one four-base chromosome and one scaffold. */
	private static final String SUMMARY="#Version\t5\nchroms\t1\nbases\t4\ndefined\t4\n"
			+"undefined\t0\ncontigs\t1\nscaffolds\t1\ninterpad\t0\n"
			+"name\tdna032\nsource\tfixture\nscafprefixes\tfalse\n";
	/** Chromosome-one counts, length and zero padding used by the independent state oracle. */
	private static final String GOOD_INFO="#Version\t5\n1\t1\t1\t4\t4\t0\t0\t0\n";
	/** Same record with a nonnumeric length to fail after metadata arrays have been allocated. */
	private static final String BAD_INFO="#Version\t5\n1\t1\t1\tbad-length\t4\t0\t0\t0\n";
	/** One named scaffold on chromosome one, starting at zero and spanning four bases. */
	private static final String SCAFFOLDS="#Version\t5\n1\t1\t0\t4\tone\n";
}
