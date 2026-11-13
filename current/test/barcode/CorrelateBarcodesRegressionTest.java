package test.barcode;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

import barcode.CorrelateBarcodes;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.Read;

/**
 * Finite output-filtering and quality-boundary checks, one scenario per child JVM.
 * Expected failures can precede production stream cleanup, so this wrapper always
 * exits its JVM after reporting success or failure. Do not invoke it inside another test JVM.
 * @author Brian Bushnell, Sucrose
 */
public final class CorrelateBarcodesRegressionTest{

	/** Requires mode=baseline|fixed, case=filter|quality, scenario=name and a new dir=path. */
	public static void main(final String[] args){
		int exit=1;
		try{
			run(args);
			exit=0;
		}catch(final Throwable failure){
			failure.printStackTrace(System.err);
		}
		//An expected production exception may leave I/O threads alive; normal return is not containment.
		System.exit(exit);
	}

	/** Parses the independent test selectors without passing them to the production CLI. */
	private static void run(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, CorrelateBarcodesRegressionTest.class, false);
		try{
			String mode=null, name=null, scenario=null, directory=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(key.equals("case")){name=value;}
				else if(key.equals("scenario")){scenario=value;}
				else if(key.equals("dir")){directory=value;}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			check("baseline".equals(mode) || "fixed".equals(mode), "Set mode=baseline or mode=fixed.");
			check(directory!=null, "Set dir to a new directory whose parent exists.");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "The production baseline requires -ea.");
			check(Read.MAX_CALLED_QUALITY()==50, "Each scenario starts at default maximum 50; high changes it through the production CLI.");
			final Path dir=Files.createDirectory(Paths.get(directory));
			final boolean baseline=mode.equals("baseline");
			if("filter".equals(name)){checkFilter(dir, scenario, baseline);}
			else if("quality".equals(name)){checkQuality(dir, scenario, baseline);}
			else{throw new IllegalArgumentException("Set case=filter or case=quality.");}
			System.out.println("CORRELATE_BARCODES case="+name+" scenario="+scenario+" mode="+mode+" failed=0");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Uses literal retained records, including the minimum-threshold equality control. */
	private static void checkFilter(final Path dir, final String scenario, final boolean baseline) throws IOException{
		final boolean single="single".equals(scenario);
		final int average, minimum;
		final String expected1, expected2;
		if("none".equals(scenario)){
			average=minimum=0;
			expected1=ALL1; expected2=ALL2;
		}else if("minimum".equals(scenario) || single){
			average=0; minimum=10;
			expected1=baseline ? ALL1 : KEEP1+MID1+BOUNDARY1;
			expected2=baseline ? ALL2 : KEEP2+MID2+BOUNDARY2;
		}else if("average".equals(scenario) || "both".equals(scenario)){
			average=30; minimum="both".equals(scenario) ? 10 : 0;
			expected1=KEEP1; expected2=KEEP2;
		}else{throw new IllegalArgumentException("Filter scenario must be none, minimum, average, both or single.");}
		final Path in1=dir.resolve("input1.fq"), in2=dir.resolve("input2.fq");
		final Path out1=dir.resolve("output1.fq"), out2=dir.resolve("output2.fq");
		write(in1, ALL1);
		if(!single){write(in2, ALL2);}
		invoke(in1, single ? null : in2, out1, single ? null : out2, null, average, minimum);
		checkText("First output", expected1, read(out1));
		if(!single){checkText("Second output", expected2, read(out2));}
		else{check(!Files.exists(out2), "Single-end case must not create a second output.");}
		System.out.println("FILTER paired="+!single+" maq="+average+" mmq="+minimum+" output_records_per_file="+
				("none".equals(scenario) || (average==0 && baseline) ? 4 : average>0 ? 1 : 3));
	}

	/** Isolates Q50 axes; high separately exercises the configured upper raw-quality bound. */
	private static void checkQuality(final Path dir, final String scenario, final boolean baseline) throws IOException{
		final boolean high="high".equals(scenario);
		check("control".equals(scenario) || "barcode".equals(scenario) || "read1".equals(scenario) || "read2".equals(scenario) || high,
				"Quality scenario must be control, barcode, read1, read2 or high.");
		//Uniform Q93 has error probability 10^-9.3; the conversion saturates at 60 below 10^-6.
		//These are literal expectations, not values computed by the method being tested.
		final int bar=high ? 60 : "barcode".equals(scenario) ? 50 : 49;
		final int first=high ? 60 : "read1".equals(scenario) ? 50 : 49;
		final int second=high ? 60 : "read2".equals(scenario) ? 50 : 49;
		final String prefix="@ACGT_"+(high ? "~~~~" : bar==50 ? "SSSS" : "RRRR")+"_quality/";
		final String input1=prefix+"1\nACGT\n+\n"+(high ? "~~~~" : first==50 ? "SSSS" : "RRRR")+"\n";
		final String input2=prefix+"2\nTGCA\n+\n"+(high ? "~~~~" : second==50 ? "SSSS" : "RRRR")+"\n";
		final Path in1=dir.resolve("input1.fq"), in2=dir.resolve("input2.fq");
		final Path out1=dir.resolve("output1.fq"), out2=dir.resolve("output2.fq"), correlation=dir.resolve("correlation.tsv");
		final Path averageHistogram=dir.resolve("average.tsv"), minimumHistogram=dir.resolve("minimum.tsv");
		write(in1, input1); write(in2, input2);
		final boolean expectFailure=baseline && !"control".equals(scenario);
		boolean observedFailure=false;
		try{
			if(high){
				invoke(in1, in2, out1, out2, correlation, 0, 0,
						"maxcalledquality=93", "baqhist="+averageHistogram, "bmqhist="+minimumHistogram);
			}else{invoke(in1, in2, out1, out2, correlation, 0, 0);}
		}catch(final ArrayIndexOutOfBoundsException failure){
			if(!expectFailure){throw failure;}
			boolean inProcess=false;
			for(final StackTraceElement frame : failure.getStackTrace()){
				inProcess|=frame.getClassName().equals("barcode.CorrelateBarcodes") && frame.getMethodName().equals("process");
			}
			check(failure.getClass()==ArrayIndexOutOfBoundsException.class && inProcess,
					"Expected the exact array-bounds exception from CorrelateBarcodes.process, observed "+failure);
			observedFailure=true;
			System.out.println("QUALITY exception="+failure.getClass().getName()+" in_process="+inProcess);
		}
		check(observedFailure==expectFailure, "Quality-bound exception expectation did not match mode/scenario.");
		if(high){check(Read.MAX_CALLED_QUALITY()==93, "Production must apply maxcalledquality=93 for the upper-domain case.");}
		if(!observedFailure){
			checkText("First quality output", input1, read(out1));
			checkText("Second quality output", input2, read(out2));
			checkCorrelation(correlation, first, second, bar);
			if(high){
				checkHistogram(averageHistogram, 60);
				checkHistogram(minimumHistogram, 93);
				System.out.println("QUALITY_HIGH maximum_called=93 raw_quality=93 average_bin=60 minimum_bin=93");
			}
			System.out.println("QUALITY read1="+first+" read2="+second+" barcode="+bar+" observations_per_mate=1");
		}
	}

	/** One fragment contributes to one literal histogram bin; every other row must be empty. */
	private static void checkHistogram(final Path path, final int expectedBin) throws IOException{
		final List<String> lines=Files.readAllLines(path, StandardCharsets.US_ASCII);
		check(!lines.isEmpty() && lines.get(0).equals("#Quality\tcount\tfraction"), "Unexpected quality histogram header.");
		long total=0;
		int previous=-1;
		for(int i=1; i<lines.size(); i++){
			final String[] fields=lines.get(i).split("\t", -1);
			check(fields.length==3, "Histogram row must contain three fields: "+lines.get(i));
			final int quality=Integer.parseInt(fields[0]);
			check(quality>previous, "Histogram bins must be nonnegative and strictly increasing.");
			previous=quality;
			final long expected=quality==expectedBin ? 1 : 0;
			final long count=Long.parseLong(fields[1]);
			check(count==expected && Double.parseDouble(fields[2])==expected,
					"Only literal histogram bin "+expectedBin+" may contain the single observation.");
			total+=count;
		}
		check(total==1, "Histogram must retain exactly one fragment at bin "+expectedBin+"; observed "+total);
	}

	/** Requires exactly one populated bin per mate without constraining extra zero rows. */
	private static void checkCorrelation(final Path path, final int first, final int second, final int bar) throws IOException{
		final List<String> lines=Files.readAllLines(path, StandardCharsets.US_ASCII);
		check(!lines.isEmpty() && lines.get(0).equals("#Read1_Q\tBar_Q\tstdev\tcount\tRead2_Q\tBar_Q\tstdev\tcount"),
				"Correlation header must identify both mates and their barcode-quality counts.");
		long total1=0, total2=0;
		int previous=-1;
		for(int i=1; i<lines.size(); i++){
			final String[] fields=lines.get(i).split("\t", -1);
			check(fields.length==8, "Correlation row must contain eight fields: "+lines.get(i));
			final int row=Integer.parseInt(fields[0]);
			check(row>previous && Integer.parseInt(fields[4])==row, "Correlation rows must increase with matching mate indices.");
			previous=row;
			final long count1=Long.parseLong(fields[3]), count2=Long.parseLong(fields[7]);
			check(count1>=0 && count1<=1 && count2>=0 && count2<=1, "One input pair permits only zero or one count per bin.");
			if(count1==1){
				check(row==first && Double.parseDouble(fields[1])==bar && Double.parseDouble(fields[2])==0,
						"First mate must occupy its literal read/barcode quality bin with zero deviation.");
			}
			if(count2==1){
				check(row==second && Double.parseDouble(fields[5])==bar && Double.parseDouble(fields[6])==0,
						"Second mate must occupy its literal read/barcode quality bin with zero deviation.");
			}
			total1+=count1; total2+=count2;
		}
		check(total1==1 && total2==1, "Correlation must retain exactly one observation per mate: "+total1+", "+total2);
	}

	/** Pins input/output encodings and pairing; fixture-only selectors never enter production arguments. */
	private static void invoke(final Path in1, final Path in2, final Path out1, final Path out2,
			final Path correlation, final int average, final int minimum, final String... extra){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in1="+in1, "out1="+out1,
				"maq="+average, "mmq="+minimum, "qin=33", "qout=33", "interleaved=f", "trimreadcomments=f",
				"ow=f", "t=1", "simd=f", "pigz=f", "unpigz=f"));
		if(in2!=null){args.add("in2="+in2);}
		if(out2!=null){args.add("out2="+out2);}
		if(correlation!=null){args.add("outcor="+correlation);}
		args.addAll(Arrays.asList(extra));
		CorrelateBarcodes.main(args.toArray(new String[args.size()]));
	}

	/** Writes only inside the newly created case directory. */
	private static void write(final Path path, final String text) throws IOException{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII));
	}
	private static String read(final Path path) throws IOException{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}
	/** Full literal comparison includes record order, mate identity, bases and body qualities. */
	private static void checkText(final String label, final String expected, final String actual){
		check(expected.equals(actual), label+" expected:\n"+expected+"observed:\n"+actual);
	}
	/** Test oracles remain enforced even if assertions are accidentally disabled elsewhere. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	private CorrelateBarcodesRegressionTest(){}

	private static final String KEEP1="@ACGT_IIII_keep/1\nACGT\n+\nIIII\n";
	private static final String LOW1="@ACGT_!III_lowmin/1\nACGT\n+\nIIII\n";
	private static final String MID1="@ACGT_5555_mid/1\nACGT\n+\nIIII\n";
	private static final String BOUNDARY1="@ACGT_++++_boundary/1\nACGT\n+\nIIII\n";
	private static final String KEEP2="@ACGT_IIII_keep/2\nTGCA\n+\nIIII\n";
	private static final String LOW2="@ACGT_!III_lowmin/2\nTGCA\n+\nIIII\n";
	private static final String MID2="@ACGT_5555_mid/2\nTGCA\n+\nIIII\n";
	private static final String BOUNDARY2="@ACGT_++++_boundary/2\nTGCA\n+\nIIII\n";
	private static final String ALL1=KEEP1+LOW1+MID1+BOUNDARY1;
	private static final String ALL2=KEEP2+LOW2+MID2+BOUNDARY2;
}
