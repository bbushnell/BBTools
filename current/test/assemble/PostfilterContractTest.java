package test.assemble;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import assemble.Postfilter;
import parse.Parse;
import parse.PreParser;
import shared.Shared;

/** Constructor-only Postfilter path and threshold contracts on isolated literal inputs.
 * Never invokes mapping or filtering; acceptance here measures early validation,
 * not demonstrated data loss or complete pipeline behavior.
 * @author Brian Bushnell, Shinobu
 */
public final class PostfilterContractTest{

	/** Saves every result before failing; outdir must be new. Optional modes select thresholds or locale handling. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, PostfilterContractTest.class, false);
		try{
			System.err.println("POSTFILTER_ASSERTIONS enabled="+PostfilterContractTest.class.desiredAssertionStatus());
			Path root=null;
			boolean thresholds=false, locales=false;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				if(pair.length!=2){throw new IllegalArgumentException("Expected named fixture option: "+arg);}
				if(pair[0].equalsIgnoreCase("outdir")){root=Paths.get(pair[1]).toAbsolutePath();}
				else if(pair[0].equalsIgnoreCase("thresholds")){thresholds=Parse.parseBoolean(pair[1]);}
				else if(pair[0].equalsIgnoreCase("locales")){locales=Parse.parseBoolean(pair[1]);}
				else{throw new IllegalArgumentException("Expected outdir=NEW_DIRECTORY, thresholds=t or locales=t: "+arg);}
			}
			if(root==null){throw new IllegalArgumentException("Specify outdir=NEW_DIRECTORY.");}
			if(thresholds && locales){throw new IllegalArgumentException("Choose only one fixture mode.");}
			Files.createDirectory(root);
			if(locales){runLocales(root); return;}
			if(thresholds){runThresholds(root); return;}
			final StringBuilder rows=new StringBuilder("case\tpassed\tdetail\n");
			int failed=0;
			for(int mode=0; mode<NAMES.length; mode++){failed+=runCase(root, mode, rows);}
			Files.write(root.resolve("checks.tsv"), rows.toString().getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
			System.out.println("POSTFILTER_PREFLIGHT checks="+NAMES.length+" failed="+failed);
			if(failed!=0){throw new AssertionError("Postfilter preflight failures; inspect "+root);}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Checks case-independent threshold names without normalizing filenames.
	 * Changes the process locale only within this isolated, serial fixture and restores it
	 * even on failure. All four constructor outcomes are saved before the final assertion. */
	private static void runLocales(final Path root) throws Exception{
		assert(root!=null && Files.isDirectory(root)) : "main must create a fresh root before locale cases write their literal inputs.";
		final Locale original=Locale.getDefault();
		final Locale display=Locale.getDefault(Locale.Category.DISPLAY);
		final Locale format=Locale.getDefault(Locale.Category.FORMAT);
		final String[] names={"english-lower", "english-upper", "turkish-lower", "turkish-upper"};
		final StringBuilder rows=new StringBuilder("case\tvalues\tpreserved\tnoOutput\tpassed\n");
		int failed=0;
		try{
			for(int i=0; i<names.length; i++){
				Locale.setDefault(i<2 ? Locale.ENGLISH : new Locale("tr", "TR"));
				final Path dir=Files.createDirectory(root.resolve(names[i]));
				final Path reads=dir.resolve("ReadsI.fq"), reference=dir.resolve("ReferenceI.fa");
				final Path output=dir.resolve("CleanI.fa"), coverage=dir.resolve("CoverageI.tsv");
				Files.write(reads, READS, StandardOpenOption.CREATE_NEW);
				Files.write(reference, REFERENCE, StandardOpenOption.CREATE_NEW);
				final boolean upper=(i&1)!=0;
				final String[] args={"in="+reads, "ref="+reference, "out="+output,
						"covstats="+coverage, "showstats=f", "overwrite=f",
						(upper ? "MINCOV=7.25" : "mincov=7.25"),
						(upper ? "MINPERCENT=73" : "minpercent=73")};
				final Postfilter filter=new Postfilter(args, false);
				final boolean values=sameDouble(filter, "minCoverage", 7.25) &&
						sameDouble(filter, "minCoveredPercent", 73);
				final boolean preserved=Arrays.equals(READS, Files.readAllBytes(reads)) &&
						Arrays.equals(REFERENCE, Files.readAllBytes(reference));
				final boolean noOutput=!Files.exists(output) && !Files.exists(coverage);
				final boolean passed=values && preserved && noOutput;
				rows.append(names[i]).append('\t').append(values).append('\t').append(preserved)
						.append('\t').append(noOutput).append('\t').append(passed).append('\n');
				if(!passed){failed++;}
			}
		}finally{
			Locale.setDefault(original);
			Locale.setDefault(Locale.Category.DISPLAY, display);
			Locale.setDefault(Locale.Category.FORMAT, format);
		}
		Files.write(root.resolve("checks.tsv"), rows.toString().getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		System.out.print(rows);
		System.out.println("POSTFILTER_LOCALES checks="+names.length+" failed="+failed);
		if(failed!=0){throw new AssertionError("Postfilter locale failures; inspect "+root);}
	}

	/** Rejects effective NaNs while retaining non-NaN values and last-option-wins parsing.
	 * Uses constructors only; no mapper, filter, read worker or output writer is started.
	 * Each case has independent paths and all results are saved before the final failure. */
	private static void runThresholds(final Path root) throws Exception{
		assert(root!=null && Files.isDirectory(root)) : "main must create a fresh root before independent threshold cases write their inputs.";
		final String[] names={"nan-coverage", "nan-percent", "finite-control", "disabled-control",
				"infinity-control", "last-option-control"};
		final String[][] options={{"mincov=NaN"}, {"minpercent=NaN"}, {"mincov=2.5", "minpercent=95"},
				{"mincov=-1", "minpercent=0"}, {"mincov=Infinity", "minpercent=Infinity"},
				{"mincov=NaN", "mincov=2", "minpercent=NaN", "minpercent=95"}};
		final double[] coverages={Double.NaN, 2, 2.5, -1, Double.POSITIVE_INFINITY, 2};
		final double[] percentages={95, Double.NaN, 95, 0, Double.POSITIVE_INFINITY, 95};
		final StringBuilder rows=new StringBuilder("case\taccepted\texpected\tvalues\tpreserved\tnoOutput\tdiagnostic\tpassed\n");
		int failed=0;
		for(int i=0; i<names.length; i++){
			final Path dir=Files.createDirectory(root.resolve(names[i]));
			final Path reads=dir.resolve("reads.fq"), reference=dir.resolve("reference.fa");
			final Path output=dir.resolve("clean.fa"), coverage=dir.resolve("coverage.tsv");
			Files.write(reads, READS, StandardOpenOption.CREATE_NEW);
			Files.write(reference, REFERENCE, StandardOpenOption.CREATE_NEW);
			final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+reads, "ref="+reference,
					"out="+output, "covstats="+coverage, "showstats=f", "overwrite=f"));
			args.addAll(Arrays.asList(options[i]));
			Postfilter filter=null;
			IllegalArgumentException failure=null;
			try{filter=new Postfilter(args.toArray(new String[0]), false);}
			catch(IllegalArgumentException e){failure=e;}
			final boolean accepted=filter!=null, expected=i>=2;
			final boolean values=!accepted || sameDouble(filter, "minCoverage", coverages[i]) &&
					sameDouble(filter, "minCoveredPercent", percentages[i]);
			final boolean preserved=Arrays.equals(READS, Files.readAllBytes(reads)) &&
					Arrays.equals(REFERENCE, Files.readAllBytes(reference));
			final boolean noOutput=!Files.exists(output) && !Files.exists(coverage);
			final boolean diagnostic=failure==null || "mincov and minpercent must not be NaN.".equals(failure.getMessage());
			final boolean passed=accepted==expected && values && preserved && noOutput && diagnostic;
			rows.append(names[i]).append('\t').append(accepted).append('\t').append(expected)
					.append('\t').append(values).append('\t').append(preserved).append('\t').append(noOutput)
					.append('\t').append(diagnostic).append('\t').append(passed).append('\n');
			if(failure!=null){System.err.println(names[i]); failure.printStackTrace(System.err);}
			if(!passed){failed++;}
		}
		Files.write(root.resolve("checks.tsv"), rows.toString().getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		System.out.print(rows);
		System.out.println("POSTFILTER_THRESHOLDS checks="+names.length+" failed="+failed);
		if(failed!=0){throw new AssertionError("Postfilter threshold preflight failures; inspect "+root);}
	}

	/** Observes the private value that process() forwards, without invoking either child tool. */
	private static boolean sameDouble(final Postfilter filter, final String name, final double expected) throws Exception{
		assert(filter!=null && name!=null) : "Only a successfully constructed filter exposes the effective threshold forwarded by process().";
		final java.lang.reflect.Field field=Postfilter.class.getDeclaredField(name);
		field.setAccessible(true);
		return Double.doubleToLongBits(field.getDouble(filter))==Double.doubleToLongBits(expected);
	}

	/** Requires early rejection of destructive aliases/overwrite refusal, and preserves all inputs.
	 * Ordinary, distinct-dirty and null/empty-dirty cases are controls; equals-input probes parsing.
	 * @return Zero for the expected constructor outcome and unchanged files, otherwise one
	 */
	private static int runCase(final Path root, final int mode, final StringBuilder rows) throws Exception{
		assert(mode>=0 && mode<NAMES.length) : "Each preflight mode must map to exactly one named result row.";
		final Path dir=Files.createDirectory(root.resolve(NAMES[mode]));
		final Path reads=dir.resolve(mode==2 ? "reads=sample.fq" : "reads.fq");
		final Path reference=dir.resolve("reference.fa"), output=dir.resolve("clean.fa"), coverage=dir.resolve("coverage.tsv");
		final Path existing=dir.resolve("existing-dirty.fa");
		Files.write(reads, READS, StandardOpenOption.CREATE_NEW);
		Files.write(reference, REFERENCE, StandardOpenOption.CREATE_NEW);
		if(mode==7){Files.write(existing, SENTINEL, StandardOpenOption.CREATE_NEW);}
		Path dirty=null;
		if(mode==1){dirty=dir.resolve("dirty.fa");}
		else if(mode==3){dirty=reads;}
		else if(mode==4){dirty=reference;}
		else if(mode==5){dirty=coverage;}
		else if(mode==6){dirty=output;}
		else if(mode==7){dirty=existing;}
		final ArrayList<String> options=new ArrayList<String>();
		options.add("in="+reads); options.add("ref="+reference);
		options.add("out="+output); options.add("covstats="+coverage);
		options.add("showstats=f"); options.add("overwrite="+(mode==7 ? "f" : "t"));
		if(mode==8){options.add("outdirty=null");}
		else if(mode==9){options.add("outdirty=");}
		else if(dirty!=null){options.add("outdirty="+dirty);}
		RuntimeException failure=null;
		try{new Postfilter(options.toArray(new String[0]), false);}
		catch(RuntimeException caught){
			// File-preflight rejection is expected; null/index/setup failures are not evidence.
			if(caught.getClass()!=RuntimeException.class && !(caught instanceof IllegalArgumentException)){throw caught;}
			failure=caught;
		}
		final boolean accepted=failure==null, expected=mode<3 || mode==8 || mode==9;
		final boolean preserved=Arrays.equals(READS, Files.readAllBytes(reads)) &&
				Arrays.equals(REFERENCE, Files.readAllBytes(reference)) &&
				(mode!=7 || Arrays.equals(SENTINEL, Files.readAllBytes(existing)));
		final boolean noOutput=!Files.exists(output) && !Files.exists(coverage) &&
				(mode!=1 || !Files.exists(dirty));
		final boolean explained=failure==null || failure.getMessage()!=null && !failure.getMessage().isEmpty();
		final boolean passed=accepted==expected && preserved && noOutput && explained;
		rows.append(NAMES[mode]).append('\t').append(passed ? 1 : 0)
				.append("\taccepted=").append(accepted).append(" expected=").append(expected)
				.append(" preserved=").append(preserved).append(" noOutput=").append(noOutput)
				.append(" explained=").append(explained).append(" exception=")
				.append(failure==null ? "none" : failure.getClass().getSimpleName()).append('\n');
		if(failure!=null){System.err.println(NAMES[mode]); failure.printStackTrace(System.err);}
		return passed ? 0 : 1;
	}

	/** Fixed scenarios; only constructor validation runs, never the production pipeline. */
	private static final String[] NAMES={"ordinary-control", "dirty-control", "equals-input",
			"dirty-alias-reads", "dirty-alias-reference", "dirty-alias-coverage", "dirty-alias-clean", "dirty-overwrite-false",
			"dirty-null-control", "dirty-empty-control"};
	/** Valid tiny read/reference records; exact bytes must survive every case. */
	private static final byte[] READS="@r\nACGT\n+\nIIII\n".getBytes(StandardCharsets.UTF_8);
	private static final byte[] REFERENCE=">ref\nACGT\n".getBytes(StandardCharsets.UTF_8);
	/** Nonempty existing output used only by the overwrite-refusal case. */
	private static final byte[] SENTINEL=">keep\nTGCA\n".getBytes(StandardCharsets.UTF_8);
}
