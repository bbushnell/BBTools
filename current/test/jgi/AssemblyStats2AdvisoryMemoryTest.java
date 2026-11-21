package test.jgi;

import java.io.ByteArrayInputStream;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import dna.FastaToChromArrays2;
import fileIO.ReadWrite;
import jgi.AssemblyStats2;
import parse.Parse;
import shared.Shared;

/**
 * Exact finite advisory-byte observations for AssemblyStats2's private
 * bbmapMemoryBytes helper, kept separate from minimum-heap or formatted-MB claims.
 * @author Jean
 */
public final class AssemblyStats2AdvisoryMemoryTest{

	/** Requires phase=baseline|fixed and a fresh output directory. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown advisory fixture option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a fresh outdir.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/assemblystats2advisorymemorytest.sh with assertions enabled.");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		Files.createDirectory(root);
		new AssemblyStats2AdvisoryMemoryTest(root, phase).run();
	}

	/** Resolves only the private advisory helper. */
	private AssemblyStats2AdvisoryMemoryTest(final Path root_, final String phase_)
			throws Exception{
		root=root_;
		phase=phase_;
		memory=AssemblyStats2.class.getDeclaredMethod(
			"bbmapMemoryBytes", long[].class, Long.TYPE, Long.TYPE, Integer.TYPE);
		memory.setAccessible(true);
		require(Modifier.isPrivate(memory.getModifiers()) &&
			Modifier.isStatic(memory.getModifiers()) && memory.getReturnType()==Long.TYPE,
			"The advisory oracle calls private static long bbmapMemoryBytes directly.");
		Files.write(root.resolve("bbmap-memory.tsv"), (HEADER+"\n").getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
	}

	/** Literal returned-byte expectations were selected independently in the canonical advisory plan. */
	private void run() throws Exception{
		for(final Row row : ROWS){runRow(row);}
		final boolean match=rows==7 && parserRows==7 &&
			advisoryRows==("baseline".equals(phase) ? 2 : 7) &&
			unexpected==0 && restored;
		final String summary="ASSEMBLY_STATS2_ADVISORY_MEMORY_CAPTURE_"+
			(match ? "PASS" : "FAIL")+" phase="+phase+" rows="+rows+
			" parserCorrectRows="+parserRows+" advisoryCorrectRows="+advisoryRows+
			" unexpectedObservations="+unexpected+
			" dependencyStateRestored="+restored+"\n";
		Files.write(root.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		System.out.print(summary);
		if(!match){throw new AssertionError(summary+"Inspect literal byte oracles in "+root);}
	}

	/** Consumes the suffix as FASTA sequence, not as a header, line ending, separator or quality. */
	private void runRow(final Row row) throws Exception{
		final String inputName=row.name+".input";
		final byte[] input=(">one\nAAAAAAAAAA"+row.suffix+"\n>four\nAAAAAAAAAA\n")
			.getBytes(StandardCharsets.US_ASCII);
		Files.write(root.resolve(inputName), input, StandardOpenOption.CREATE_NEW);

		final int[] before=padding();
		final StaticState saved=new StaticState();
		long[] counts=null;
		long observedBytes=-1, observedHeader=-1;
		int observedScaffolds=-1;
		boolean countsPreserved=false, parserCorrect=false, paddingRestored=false;
		int[] after=before;
		Throwable failure=null;
		try{
			require(Arrays.equals(before, PADDING),
				"Advisory literals pin MID/START/END/MAX_LENGTH: "+Arrays.toString(before));
			final AssemblyStats2 stats=newStats();
			counts=stats.countFasta(new ByteArrayInputStream(input), null);
			observedHeader=field("HEADERLENSUM").getLong(null);
			observedScaffolds=((ArrayList<?>)field("tlist").get(stats)).size();
			parserCorrect=Arrays.equals(row.correctCounts, counts) &&
				observedHeader==7 && observedScaffolds==2;
			if(parserCorrect){
				final long[] snapshot=counts.clone();
				observedBytes=((Long)memory.invoke(null, counts, 2L, 7L, 5)).longValue();
				countsPreserved=Arrays.equals(snapshot, counts);
			}
		}catch(final InvocationTargetException caught){
			failure=caught.getCause();
		}catch(final Throwable caught){
			failure=caught;
		}finally{
			after=padding();
			try{
				saved.restore();
			}finally{
				FastaToChromArrays2.MID_PADDING=before[0];
				FastaToChromArrays2.START_PADDING=before[1];
				FastaToChromArrays2.END_PADDING=before[2];
				FastaToChromArrays2.MAX_LENGTH=before[3];
				paddingRestored=Arrays.equals(before, padding());
			}
		}

		final boolean correct=failure==null && observedBytes==row.intendedBytes;
		final long expected="baseline".equals(phase) ? row.baselineBytes : row.intendedBytes;
		final boolean match=failure==null && parserCorrect && countsPreserved &&
			Arrays.equals(before, after) && paddingRestored && observedBytes==expected;
		rows++;
		if(parserCorrect){parserRows++;}
		if(correct){advisoryRows++;}
		if(!match){unexpected++;}
		restored&=paddingRestored;

		final String text=row.name+"\t"+phase+"\t"+inputName+
			"\t"+Arrays.toString(row.correctCounts)+"\t"+Arrays.toString(counts)+
			"\t2\t"+observedScaffolds+"\t7\t"+observedHeader+
			"\t5\t"+row.baselineBytes+"\t"+row.intendedBytes+"\t"+observedBytes+
			"\t"+countsPreserved+"\t"+parserCorrect+"\t"+correct+"\t"+match+
			"\t"+Arrays.toString(before)+"\t"+Arrays.toString(after)+
			"\t"+paddingRestored+"\t"+type(failure)+
			"\t"+clean(failure==null ? "" : String.valueOf(failure.getMessage()))+"\n";
		Files.write(root.resolve("bbmap-memory.tsv"), text.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.APPEND);
	}

	/** Creates a scanner instance with the literal cutoff and N-split options. */
	private static AssemblyStats2 newStats(){
		return new AssemblyStats2(new String[]{"in=stdin", "cutoff=4", "n=2",
			"format=1", "gcformat=1", "gcfraction=f", "showbbmap=f",
			"showspeed=f", "n_=f"});
	}

	/** Snapshot order is explicit in the TSV and independently pinned. */
	private static int[] padding(){
		return new int[]{FastaToChromArrays2.MID_PADDING, FastaToChromArrays2.START_PADDING,
			FastaToChromArrays2.END_PADDING, FastaToChromArrays2.MAX_LENGTH};
	}

	/** Resolves one private field for state restoration or exact observations. */
	private static Field field(final String name) throws Exception{
		final Field field=AssemblyStats2.class.getDeclaredField(name);
		field.setAccessible(true);
		return field;
	}

	/** Makes one exception message fit in a TSV cell. */
	private static String clean(final String value){
		return value.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');
	}

	/** Records absent exceptions explicitly. */
	private static String type(final Throwable error){
		return error==null ? "none" : error.getClass().getName();
	}

	/** Rejects fixture misuse. */
	private static void require(final boolean correct, final String message){
		if(!correct){throw new IllegalArgumentException(message);}
	}

	/** Saves and restores every mutable static field used by direct scanner rows. */
	private static final class StaticState{
		StaticState() throws Exception{
			for(final StaticRef ref : STATIC_FIELDS){
				final Field field=ref.clazz.getDeclaredField(ref.name);
				field.setAccessible(true);
				values.put(field, field.get(null));
			}
		}

		void restore() throws Exception{
			for(final Map.Entry<Field, Object> entry : values.entrySet()){
				entry.getKey().set(null, entry.getValue());
			}
		}

		private final LinkedHashMap<Field, Object> values=new LinkedHashMap<Field, Object>();
	}

	/** A mutable static field whose value must be restored after direct constructor calls. */
	private static final class StaticRef{
		StaticRef(final Class<?> clazz_, final String name_){
			clazz=clazz_;
			name=name_;
		}

		final Class<?> clazz;
		final String name;
	}

	/** One advisory-memory row with baseline and intended private-helper byte counts. */
	private static final class Row{
		Row(final String name_, final String suffix_, final long[] correctCounts_,
				final long baselineBytes_, final long intendedBytes_){
			name=name_;
			suffix=suffix_;
			correctCounts=correctCounts_;
			baselineBytes=baselineBytes_;
			intendedBytes=intendedBytes_;
		}

		final String name, suffix;
		final long[] correctCounts;
		final long baselineBytes, intendedBytes;
	}

	private final Path root;
	private final String phase;
	private final Method memory;
	private int rows, parserRows, advisoryRows, unexpected;
	private boolean restored=true;

	private static final String HEADER="case\tphase\tinput\tcorrectCounts\tactualCounts"+
		"\tscaffolds\tscaffoldsObserved\theaderBytes\theaderBytesObserved\tk"+
		"\tbaselineBytes\tintendedBytes\tactualBytes\tcountsPreserved\tparserCorrect"+
		"\tadvisoryCorrect\texpectedObservationMatched\tpaddingBefore\tpaddingAfter"+
		"\tpaddingRestored\texceptionClass\texceptionMessage";

	private static final Row[] ROWS={
		new Row("defined-only", "", new long[]{20, 0, 0, 0, 0, 0, 0, 0}, 7386, 7386),
		new Row("iupac-only", "RRRRRRRRRR", new long[]{20, 0, 0, 0, 10, 0, 0, 0}, 7401, 7401),
		new Row("n-only", "NNNNNNNNNN", new long[]{20, 0, 0, 0, 0, 10, 0, 0}, 7386, 7401),
		new Row("iupac3-n7", "RRRNNNNNNN", new long[]{20, 0, 0, 0, 3, 7, 0, 0}, 7390, 7401),
		new Row("other-only", "??????????", new long[]{20, 0, 0, 0, 0, 0, 10, 0}, 7386, 7401),
		new Row("consumed-control-only", "++++++++++", new long[]{20, 0, 0, 0, 0, 0, 0, 10}, 7386, 7401),
		new Row("all-undefined", "RRRNNN??++", new long[]{20, 0, 0, 0, 3, 3, 2, 2}, 7390, 7401)
	};

	private static final int[] PADDING={300, 8000, 8000, 536670912};

	private static final StaticRef[] STATIC_FIELDS={
		new StaticRef(AssemblyStats2.class, "GCFORMAT"),
		new StaticRef(AssemblyStats2.class, "GC_FRACTION"),
		new StaticRef(AssemblyStats2.class, "FORMAT"),
		new StaticRef(AssemblyStats2.class, "cutoff"),
		new StaticRef(AssemblyStats2.class, "LIMSUM"),
		new StaticRef(AssemblyStats2.class, "HEADERLENSUM"),
		new StaticRef(AssemblyStats2.class, "bbmapkmer"),
		new StaticRef(AssemblyStats2.class, "overwrite"),
		new StaticRef(AssemblyStats2.class, "append"),
		new StaticRef(AssemblyStats2.class, "useheader"),
		new StaticRef(AssemblyStats2.class, "addfilename"),
		new StaticRef(AssemblyStats2.class, "printL90"),
		new StaticRef(AssemblyStats2.class, "swapNL"),
		new StaticRef(AssemblyStats2.class, "printExtended"),
		new StaticRef(AssemblyStats2.class, "logSumOffset"),
		new StaticRef(AssemblyStats2.class, "powSumPower"),
		new StaticRef(AssemblyStats2.class, "logSumBase"),
		new StaticRef(AssemblyStats2.class, "squareLog"),
		new StaticRef(AssemblyStats2.class, "logPower"),
		new StaticRef(AssemblyStats2.class, "showspeed"),
		new StaticRef(AssemblyStats2.class, "printheadersize"),
		new StaticRef(AssemblyStats2.class, "skipDuplicateLines"),
		new StaticRef(AssemblyStats2.class, "N_UNDERSCORE"),
		new StaticRef(AssemblyStats2.class, "assemblyScoreMinLen"),
		new StaticRef(AssemblyStats2.class, "assemblyScoreMaxLen"),
		new StaticRef(AssemblyStats2.class, "alignedFraction"),
		new StaticRef(AssemblyStats2.class, "printAssemblyScore"),
		new StaticRef(AssemblyStats2.class, "lastLogSum"),
		new StaticRef(AssemblyStats2.class, "lastL90"),
		new StaticRef(AssemblyStats2.class, "lastL50"),
		new StaticRef(AssemblyStats2.class, "lastSize"),
		new StaticRef(AssemblyStats2.class, "lastContigs"),
		new StaticRef(AssemblyStats2.class, "lastMaxContig"),
		new StaticRef(AssemblyStats2.class, "outstream"),
		new StaticRef(ReadWrite.class, "USE_PIGZ"),
		new StaticRef(ReadWrite.class, "USE_UNPIGZ"),
		new StaticRef(Shared.class, "mainClass"),
		new StaticRef(Shared.class, "COMMAND_LINE")
	};
}
