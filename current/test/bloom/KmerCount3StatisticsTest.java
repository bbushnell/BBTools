package test.bloom;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import bloom.KCountArray;
import bloom.KCountArray2;
import bloom.KCountArray3;
import bloom.KmerCount3;
import bloom.KmerCount4;
import bloom.KmerCount5;
import bloom.KmerCount6;
import bloom.KmerCountAbstract;

/** Checks the actual KmerCount3/4/5/6 formatters with real small tables and a long histogram.
 * The boundary adapter represents a possible k16 histogram without allocating its
 * 2^32-cell table. It tests reporting, not an end-to-end giant counting workload.
 * All cases retain fixed success expectations and raw reports on the broken baseline.
 * Run once under each assertion mode; no fixture assertion depends on -ea.
 * @author Brian Bushnell, Sucrose */
public final class KmerCount3StatisticsTest{

	/** Accepts one fresh output directory and runs all eleven cases despite report failures. */
	public static void main(final String[] args) throws Exception{
		if(args.length!=1 || !args[0].startsWith("outdir=")){
			throw new IllegalArgumentException("Usage: outdir=<new directory>");
		}
		final Path root=Paths.get(args[0].substring(7)).toAbsolutePath();
		check(!Files.exists(root), "Fresh output prevents stale reports from satisfying a new capture.");
		Files.createDirectories(root);
		int failures=0;
		final KCountArray2 real=new KCountArray2(1024, 4);
		real.increment(42, 8);
		real.increment(43, 9);
		check(real.read(42)==8 && real.read(43)==9 && real.cells==1024,
				"The real control must contain two distinct unsaturated cells at depths8 and9 in1024cells.");
		failures+=runCase(root, "small-real-table", real, 1022, 2, "99.805%", "0.195%");

		final long[] frequencies=new long[100000];
		frequencies[0]=2147483641L;
		frequencies[8]=2147483648L;
		frequencies[9]=7;
		check(frequencies[0]+frequencies[8]+frequencies[9]==4294967296L,
				"The literal boundary histogram represents exactly4^16 cells; do not allocate that giant table.");
		failures+=runCase(root, "long-histogram-boundary", new HistogramAdapter(frequencies),
				2147483641L, 2147483655L, "50.000%", "50.000%");
		for(int version=4; version<=6; version++){
			final KCountArray2 direct=version==4 ? new KCountArray2(1024, 4) : null;
			final KCountArray generic=version==4 ? null : new KCountArray3(1024, 4);
			if(direct!=null){
				direct.increment(42, 8);
				direct.increment(43, 9);
			}else{
				generic.increment(42, 8);
				generic.increment(43, 9);
			}
			check((direct!=null ? direct.read(42) : generic.read(42))==8 &&
					(direct!=null ? direct.read(43) : generic.read(43))==9,
					"Each real summary control needs two distinct unsaturated cells at depths8 and9.");
			failures+=runCase(root, "kmer"+version+"-small-real-table", version, direct, generic,
					1022, 2, "99.805%", "0.195%", 17, summaries(17, 2, "8.500", "100.000%"));
			failures+=runCase(root, "kmer"+version+"-long-histogram-boundary", version,
					version==4 ? new HistogramAdapter(frequencies.clone()) : null,
					version==4 ? null : new GenericHistogramAdapter(frequencies.clone()),
					2147483641L, 2147483655L, "50.000%", "50.000%", 17179869247L,
					summaries(17179869247L, 2147483655L, "8.000", "100.000%"));
			failures+=runCase(root, "kmer"+version+"-empty-real-table", version,
					version==4 ? new KCountArray2(1024, 4) : null,
					version==4 ? null : new KCountArray3(1024, 4),
					1024, 0, "100.000%", "0.000%", 11, summaries(11, 0, "0.000", "0.000%"));
		}
		if(failures>0){throw new AssertionError("KmerCount statistics failed "+failures+" of11 cases.");}
		System.out.println("KMER_COUNT3_STATISTICS_TEST_PASS cases=11");
	}

	/** Preserves the original KmerCount3 case names, histogram inputs and literal output oracles. */
	private static int runCase(final Path root, final String name, final KCountArray2 counts,
			final long zero, final long grouped, final String zeroPercent, final String groupedPercent) throws Exception{
		return runCase(root, name, 3, counts, null, zero, grouped, zeroPercent, groupedPercent,
				KmerCountAbstract.keysCounted, new String[0]);
	}

	/** Captures the selected real formatter; exactly one table argument matches its public API. */
	private static int runCase(final Path root, final String name, final int version,
			final KCountArray2 direct, final KCountArray generic, final long zero, final long grouped,
			final String zeroPercent, final String groupedPercent, final long keys, final String[] summaries) throws Exception{
		check(version>=3 && version<=6 && (version<=4 ? direct!=null && generic==null : direct==null && generic!=null),
				"KmerCount3/4 accept KCountArray2; KmerCount5/6 accept KCountArray.");
		final Path dir=Files.createDirectory(root.resolve(name));
		final long[] before=(direct!=null ? direct.transformToFrequency() : generic.transformToFrequency()).clone();
		final ArrayList<String> histogram=new ArrayList<String>();
		histogram.add("depth\tcells");
		for(int i=0; i<before.length; i++){if(before[i]!=0){histogram.add(i+"\t"+before[i]);}}
		Files.write(dir.resolve("input-histogram.tsv"), histogram, StandardCharsets.UTF_8);
		final ArrayList<String> oracle=new ArrayList<String>();
		oracle.add("row\tpercent\tcells");
		for(int i=0; i<LABELS.length; i++){
			oracle.add(LABELS[i]+"\t"+percent(i, zeroPercent, groupedPercent)+"\t"+expected(i, zero, grouped));
		}
		oracle.addAll(Arrays.asList(summaries));
		Files.write(dir.resolve("oracle.tsv"), oracle, StandardCharsets.UTF_8);

		final ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		final PrintStream prior=System.out;
		final Locale priorLocale=Locale.getDefault();
		final Locale priorFormat=Locale.getDefault(Locale.Category.FORMAT);
		final Locale priorDisplay=Locale.getDefault(Locale.Category.DISPLAY);
		final long priorKeys=KmerCountAbstract.keysCounted;
		Throwable failure=null;
		try(PrintStream capture=new PrintStream(buffer, true, "UTF-8")){
			Locale.setDefault(Locale.ROOT);
			System.setOut(capture);
			KmerCountAbstract.keysCounted=keys;
			switch(version){
				case 3: KmerCount3.printStatistics(direct); break;
				case 4: KmerCount4.printStatistics(direct); break;
				case 5: KmerCount5.printStatistics(generic); break;
				case 6: KmerCount6.printStatistics(generic); break;
				default: throw new AssertionError("Unreachable formatter: "+version);
			}
		}catch(Throwable caught){failure=caught;}
		finally{
			System.setOut(prior);
			KmerCountAbstract.keysCounted=priorKeys;
			Locale.setDefault(priorLocale);
			Locale.setDefault(Locale.Category.FORMAT, priorFormat);
			Locale.setDefault(Locale.Category.DISPLAY, priorDisplay);
		}
		final String report=buffer.toString("UTF-8");
		Files.write(dir.resolve("report.txt"), report.getBytes(StandardCharsets.UTF_8));
		if(failure!=null){
			try(PrintStream errors=new PrintStream(Files.newOutputStream(dir.resolve("failure.txt")), true, "UTF-8")){
				failure.printStackTrace(errors);
			}
		}
		final boolean restored=System.out==prior && KmerCountAbstract.keysCounted==priorKeys && Locale.getDefault().equals(priorLocale) &&
				Locale.getDefault(Locale.Category.FORMAT).equals(priorFormat) &&
				Locale.getDefault(Locale.Category.DISPLAY).equals(priorDisplay);
		final boolean unchanged=Arrays.equals(before, direct!=null ? direct.transformToFrequency() : generic.transformToFrequency());
		final boolean finite=!report.contains("NaN") && !report.contains("Infinity");
		final boolean rows=checkRows(name, report, zero, grouped, zeroPercent, groupedPercent, summaries);
		final boolean pass=failure==null && restored && unchanged && finite && rows;
		final String summary="case="+name+" assertions="+KmerCount3StatisticsTest.class.desiredAssertionStatus()+
				" histogramBins="+before.length+" restored="+restored+" histogramUnchanged="+unchanged+
				" finite="+finite+" rows="+rows+" failure="+(failure==null ? "none" : failure)+" pass="+pass;
		Files.write(dir.resolve("outcome.txt"), Arrays.asList(summary), StandardCharsets.UTF_8);
		System.out.println(summary);
		return pass ? 0 : 1;
	}

	/** Checks one header, all histogram rows and every supplied summary/separator exactly once in order. */
	private static boolean checkRows(final String name, final String report, final long zero,
			final long grouped, final String zeroPercent, final String groupedPercent, final String[] summaries){
		final String[] lines=report.split("\\r?\\n", -1);
		if(lines.length!=LABELS.length+2+summaries.length || !lines[0].equals("Kmer fraction:") || !lines[lines.length-1].isEmpty()){
			System.out.println("case="+name+" invalid header/row count/terminal newline; lines="+lines.length);
			return false;
		}
		boolean pass=true;
		for(int i=0; i<LABELS.length; i++){
			final String[] fields=lines[i+1].split("\t", -1);
			final String wanted=LABELS[i]+"\t"+percent(i, zeroPercent, groupedPercent)+"\t"+expected(i, zero, grouped);
			final String actual=fields.length==3 ? fields[0].trim()+"\t"+fields[1].trim()+"\t"+fields[2].trim() : lines[i+1];
			if(fields.length!=3 || !actual.equals(wanted)){
				pass=false;
				System.out.println("case="+name+" row="+LABELS[i]+" expected="+wanted.replace('\t', '|')+
						" actual="+actual.replace('\t', '|'));
			}
		}
		for(int i=0; i<summaries.length; i++){
			final String actual=normalize(lines[LABELS.length+1+i]);
			if(!actual.equals(summaries[i])){
				pass=false;
				System.out.println("case="+name+" summary="+i+" expected="+summaries[i].replace('\t', '|')+
						" actual="+actual.replace('\t', '|'));
			}
		}
		return pass;
	}

	/** Literal summary expectations: all occupied cells exceed singleton depth. */
	private static String[] summaries(final long keys, final long occupied, final String average, final String usefulPercent){
		return new String[]{"", "Keys Counted:\t\t"+keys, "Unique:\t\t"+occupied, "Avg Sites/Key:\t\t"+average,
				"", "Singleton:\t0.000%\t0", "Useful:\t"+usefulPercent+"\t"+occupied};
	}

	/** Removes field padding only; preserves tab count and blank summary separators. */
	private static String normalize(final String line){
		final String[] fields=line.split("\t", -1);
		for(int i=0; i<fields.length; i++){fields[i]=fields[i].trim();}
		return String.join("\t", fields);
	}

	/** Literal output oracle: only zero-depth cells and the8–15 group are nonzero. */
	private static long expected(final int row, final long zero, final long grouped){
		return row==0 ? zero : row==8 ? grouped : 0;
	}

	/** Uses independently stated decimal percentages; does not reproduce production aggregation. */
	private static String percent(final int row, final String zero, final String grouped){
		return row==0 ? zero : row==8 ? grouped : "0.000%";
	}

	/** Supplies a fixed long histogram from a small valid backing table; no giant table is constructed. */
	private static final class HistogramAdapter extends KCountArray2{
		HistogramAdapter(final long[] frequencies_){
			super(1024, 4);
			frequencies=frequencies_;
		}

		@Override
		public long[] transformToFrequency(){return frequencies;}

		private final long[] frequencies;
	}

	/** Same fixed histogram boundary for the KCountArray-based KmerCount5/6 APIs. */
	private static final class GenericHistogramAdapter extends KCountArray3{
		GenericHistogramAdapter(final long[] frequencies_){
			super(1024, 4);
			frequencies=frequencies_;
		}

		@Override
		public long[] transformToFrequency(){return frequencies;}

		private final long[] frequencies;
	}

	/** Keeps setup invariants enabled in both assertion modes. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final String[] LABELS={"0", "1", "2", "3", "4", "5", "6", "7", "8-15", "16-31",
			"32-63", "64-127", "128-255", "256-511", "512-1023", "1024-2047", "2048-4095", "4096-8191",
			"8192-16383", "16384-32767", "32768-65535", "65536+"};
}
