package test.jgi;

import java.io.InputStream;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import fileIO.ReadWrite;
import jgi.AssemblyStats2;
import parse.Parse;
import shared.Shared;
import structures.LongList;

/**
 * Twelve finite FASTA EOF rows for JGI114 trailing-header finalization.
 * The baseline oracle captures duplicate storage for unterminated trailing headers.
 * @author Jean
 */
public final class AssemblyStats2FastaEofTest{

	/** Requires phase=baseline|fixed and a fresh output directory. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown FASTA EOF fixture option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null && !directory.isEmpty(), "Set a fresh outdir.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/assemblystats2fastaeoftest.sh with assertions enabled.");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		Files.createDirectory(root);
		new AssemblyStats2FastaEofTest(root, phase).run();
	}

	/** Opens the evidence table before any row writes its literal input. */
	private AssemblyStats2FastaEofTest(final Path root_, final String phase_)
			throws Exception{
		root=root_;
		phase=phase_;
		tlist=field("tlist");
		Files.write(root.resolve("fasta-eof-cases.tsv"), (HEADER+"\n")
			.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Exercises only full-buffer and one-byte fragmented reads of active or trailing-header EOF. */
	private void run() throws Exception{
		for(final int cap : new int[]{Integer.MAX_VALUE, 1}){
			for(final int length : new int[]{3, 4}){
				for(final String tail : new String[]{"none", "terminated", "partial"}){
					final String prefix=length==3 ? ">a\nAAA\n" : ">a\nAAAA\n";
					final String text=prefix+("none".equals(tail) ? "" :
						"terminated".equals(tail) ? ">tail\n" : ">tail");
					final long[] counts=length==3 ? THREE_A : FOUR_A;
					final String correct=length==3 ? THREE_STORAGE : FOUR_STORAGE;
					final String observed="partial".equals(tail) ?
						(length==3 ? THREE_DUPLICATE : FOUR_DUPLICATE) : correct;
					row("fasta-"+length+"-"+tail+"-"+cap,
						text.getBytes(StandardCharsets.US_ASCII), cap, counts, correct,
						counts, observed);
				}
			}
		}
		final boolean match=rows==12 && correctRows==("baseline".equals(phase) ? 8 : 12) &&
			unexpected==0 && restored;
		final String summary="ASSEMBLY_STATS2_FASTA_EOF_CAPTURE_"+(match ? "PASS" : "FAIL")+
			" phase="+phase+" rows="+rows+" parserCorrectRows="+correctRows+
			" unexpectedObservations="+unexpected+" dependencyStateRestored="+restored+"\n";
		Files.write(root.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		System.out.print(summary);
		if(!match){throw new AssertionError(summary+"Inspect saved FASTA EOF rows in "+root);}
	}

	/** Records exact correct and source-predicted observations without deriving either from output. */
	private void row(final String name, final byte[] inputBytes, final int cap,
			final long[] correctCounts, final String correctStorage,
			final long[] baselineCounts, final String baselineStorage) throws Exception{
		assert(correctCounts.length==8 && baselineCounts.length==8) :
			"AssemblyStats2 scanner counts have eight makeCharToNum categories.";
		final CappedInputStream input=new CappedInputStream(inputBytes, cap);
		final StaticState saved=new StaticState();
		AssemblyStats2 stats=null;
		long[] actual=null;
		String storage="unavailable";
		Throwable failure=null;
		boolean restoredRow=false;
		try{
			stats=new AssemblyStats2(new String[]{"in=stdin", "cutoff=4", "minscaf=1", "n=2",
				"format=1", "gcformat=1", "gcfraction=f", "showbbmap=f", "showspeed=f", "n_=f"});
			try{actual=stats.countFasta(input, null);}
			catch(final Throwable caught){failure=caught;}
			storage=storage(stats);
		}finally{
			saved.restore();
			restoredRow=true;
		}
		final boolean correct=failure==null && Arrays.equals(correctCounts, actual) &&
			correctStorage.equals(storage);
		final boolean readAccounting=input.singleReads==0 && input.delivered==inputBytes.length;
		final boolean match=failure==null &&
			Arrays.equals("baseline".equals(phase) ? baselineCounts : correctCounts, actual) &&
			("baseline".equals(phase) ? baselineStorage : correctStorage).equals(storage) &&
			readAccounting && restoredRow;
		rows++;
		if(correct){correctRows++;}
		if(!match){unexpected++;}
		restored&=restoredRow;
		Files.write(root.resolve(name+".input"), inputBytes, StandardOpenOption.CREATE_NEW);
		final String text=name+'\t'+phase+"\tfasta\t"+cap+'\t'+inputBytes.length+
			'\t'+Arrays.toString(correctCounts)+'\t'+correctStorage+
			'\t'+Arrays.toString(baselineCounts)+'\t'+baselineStorage+
			'\t'+Arrays.toString(actual)+'\t'+storage+'\t'+type(failure)+
			'\t'+clean(failure==null ? "" : String.valueOf(failure.getMessage()))+
			'\t'+input.bulkReads+'\t'+input.singleReads+'\t'+input.delivered+
			'\t'+readAccounting+'\t'+correct+'\t'+match+'\t'+restoredRow+'\n';
		Files.write(root.resolve("fasta-eof-cases.tsv"),
			text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
	}

	/** Reads all nonzero bins and the entire logical large-list family. */
	private String storage(final AssemblyStats2 stats) throws Exception{
		return "clist="+bins(stats, "clist")+";slist="+bins(stats, "slist")+
			";sclist1="+bins(stats, "sclist1")+";sclist2="+bins(stats, "sclist2")+
			";llist="+Arrays.toString(array(stats, field("llist")))+";triples="+triples(stats);
	}

	/** Sparse literal-friendly representation includes every nonzero logical bin. */
	private static String bins(final AssemblyStats2 stats, final String name) throws Exception{
		final long[] values=array(stats, field(name));
		final StringBuilder result=new StringBuilder();
		for(int i=0; i<values.length; i++){
			if(values[i]!=0){
				if(result.length()>0){result.append(',');}
				result.append(i).append(':').append(values[i]);
			}
		}
		return result.toString();
	}

	/** Returns the logical contents of a private AssemblyStats2 LongList field. */
	private static long[] array(final AssemblyStats2 stats, final Field field) throws Exception{
		final LongList list=(LongList)field.get(stats);
		return Arrays.copyOf(list.array, list.size);
	}

	/** Returns compact Triple length/contig-count/contig-length observations. */
	private String triples(final AssemblyStats2 stats) throws Exception{
		final StringBuilder sb=new StringBuilder();
		final ArrayList<?> list=(ArrayList<?>)tlist.get(stats);
		for(int i=0; i<list.size(); i++){
			if(i>0){sb.append(',');}
			final Object triple=list.get(i);
			final Class<?> clazz=triple.getClass();
			final Field length=clazz.getDeclaredField("length");
			final Field contigs=clazz.getDeclaredField("contigs");
			final Field contiglen=clazz.getDeclaredField("contiglen");
			length.setAccessible(true);
			contigs.setAccessible(true);
			contiglen.setAccessible(true);
			sb.append(length.getLong(triple)).append('/').append(contigs.getLong(triple)).
				append('/').append(contiglen.getLong(triple));
		}
		return sb.toString();
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

	/** Same positive-bulk-read cap contract as the accepted fragmented FASTA rows. */
	private static final class CappedInputStream extends InputStream{
		CappedInputStream(final byte[] data_, final int cap_){
			require(data_!=null && cap_>0,
				"Bounded fixture stream requires data and a positive read cap.");
			data=data_;
			cap=cap_;
		}

		@Override
		public int read(){
			singleReads++;
			if(pos>=data.length){return -1;}
			delivered++;
			return data[pos++]&0xFF;
		}

		@Override
		public int read(final byte[] target, final int offset, final int length){
			if(target==null){throw new NullPointerException();}
			if(offset<0 || length<0 || offset>target.length-length){
				throw new IndexOutOfBoundsException();
			}
			bulkReads++;
			if(length==0){return 0;}
			if(pos>=data.length){return -1;}
			final int count=Math.min(Math.min(length, cap), data.length-pos);
			System.arraycopy(data, pos, target, offset, count);
			pos+=count;
			delivered+=count;
			return count;
		}

		private final byte[] data;
		private final int cap;
		private int pos, bulkReads, singleReads, delivered;
	}

	private final Path root;
	private final String phase;
	private final Field tlist;
	private int rows, correctRows, unexpected;
	private boolean restored=true;

	private static final String HEADER="case\tphase\tformat\tchunkCap\tinputBytes"+
		"\tcorrectCounts\tcorrectStorage\tbaselineCounts\tbaselineStorage"+
		"\tactualCounts\tactualStorage\texceptionClass\texceptionMessage"+
		"\tbulkReads\tsingleReads\tdelivered\treadAccounting"+
		"\tparserCorrect\texpectedObservationMatched\tdependencyStateRestored";

	private static final long[] THREE_A={3, 0, 0, 0, 0, 0, 0, 0};
	private static final long[] FOUR_A={4, 0, 0, 0, 0, 0, 0, 0};
	private static final String THREE_STORAGE="clist=3:1;slist=3:1;"+
		"sclist1=3:1;sclist2=3:3;llist=[];triples=";
	private static final String FOUR_STORAGE=
		"clist=;slist=;sclist1=;sclist2=;llist=[4];triples=4/1/4";
	private static final String THREE_DUPLICATE="clist=3:3;slist=3:2;"+
		"sclist1=3:3;sclist2=3:9;llist=[];triples=";
	private static final String FOUR_DUPLICATE=
		"clist=;slist=;sclist1=;sclist2=;llist=[4, 4, 4];triples=4/1/4,4/2/8";

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
