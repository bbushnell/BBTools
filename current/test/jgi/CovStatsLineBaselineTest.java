package test.jgi;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLStreamHandler;
import java.net.URLStreamHandlerFactory;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import driver.MergeCoverageOTU;
import jgi.CovStatsLine;
import parse.Parse;

/** Captures coverage-header defects and their bounded repairs with literal controls.
 * Direct cases configure CovStatsLine's process-wide header serially in a fresh JVM.
 * Header uses four baseline/eight fixed children, overflow two merger children, nonfinite two filter children.
 * Input-close uses two fresh fixture JVMs for controlled original-main injection, not native-file CLI failure.
 * Arithmetic uses literal read-count/scalar controls and eight small real CLI children in both phases.
 * Input-cleanup uses six isolated original-main children to check primary failure, close state and suppression.
 * Fixed aggregation oracles use the approved unavailable median/deviation contract;
 * numeric singleton values and historical baseline expectations remain separate.
 * Baseline records current defects; fixed retains the six cases and adds separately
 * labeled post-repair contracts for optional columns, labels and invalid schemas.
 * @author Sucrose, Jean
 */
public final class CovStatsLineBaselineTest{

	/** Requires baseline/fixed and a new outdir; controlled child modes remain isolated. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, scope="header";
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException(arg);}
			if("phase".equals(pair[0])){phase=pair[1];}
			else if("outdir".equals(pair[0])){directory=pair[1];}
			else if("scope".equals(pair[0])){scope=pair[1];}
			else{throw new IllegalArgumentException(arg);}
		}
		final boolean overflow="overflow".equals(scope);
		final boolean nonfinite="nonfinite".equals(scope);
		final boolean outputclose="outputclose".equals(scope);
		final boolean inputclose="inputclose".equals(scope);
		final boolean arithmetic="arithmetic".equals(scope);
		final boolean inputcleanup="inputcleanup".equals(scope);
		final boolean median="median".equals(scope);
		final int inputcleanupChild=inputCleanupChild(scope);
		final boolean inputcloseChild="inputclose-normal-child".equals(scope) || "inputclose-failure-child".equals(scope);
		if(!overflow && !nonfinite && !outputclose && !inputclose && !arithmetic && !inputcleanup && !median
				&& inputcleanupChild<0 && !inputcloseChild && !"header".equals(scope)){
			throw new IllegalArgumentException("Unknown scope: "+scope);
		}
		check("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline or phase=fixed.");
		final boolean fixed="fixed".equals(phase);
		check(directory!=null && directory.indexOf('=')<0,
				"Set a new outdir without equals signs; the legacy merger splits filenames at equals signs.");
		final Path filter=Paths.get("filterbycoverage.sh").toAbsolutePath().normalize();
		final Path merger=Paths.get("mergeOTUs.sh").toAbsolutePath().normalize();
		check(Files.isRegularFile(filter), "Missing actual filter launcher: "+filter);
		check(Files.isRegularFile(merger), "Missing actual merger launcher: "+merger);
		final Path root=Files.createDirectory(Paths.get(directory).toAbsolutePath());
		if(inputcleanupChild>=0){
			controlledInputCleanup(root, inputcleanupChild/2, (inputcleanupChild&1)!=0);
			return;
		}
		if(inputcleanup){
			final Path workingDirectory=Paths.get(".").toAbsolutePath().normalize();
			final Path launcher=workingDirectory.resolve("dev/covstatslinebaselinetest.sh");
			for(int kind=0; kind<3; kind++){
				inputCleanupCase(root, launcher, workingDirectory, kind, false, fixed);
				inputCleanupCase(root, launcher, workingDirectory, kind, true, fixed);
			}
			System.out.println("COV_STATS_INPUTCLEANUP_"+(fixed ? "FIXED" : "BASELINE")+"_PASS checks="+checks);
			return;
		}
		if(arithmetic){
			directReadTotal(root, false, fixed);
			directReadTotal(root, true, fixed);
			final ArithmeticCase[] cases=arithmeticCases();
			for(final ArithmeticCase example : cases){directArithmetic(root, example, fixed);}
			for(final int index : new int[]{0, 2, 4, 6}){mergeArithmetic(root, merger, cases[index], fixed);}
			filterReadTotal(root, filter, false, false, fixed);
			filterReadTotal(root, filter, false, true, fixed);
			filterReadTotal(root, filter, true, false, fixed);
			filterReadTotal(root, filter, true, true, fixed);
			System.out.println("COV_STATS_ARITHMETIC_"+(fixed ? "FIXED" : "BASELINE")+"_PASS checks="+checks);
			return;
		}
		if(median){
			directMedianNumeric(root);
			directMedianNA(root, fixed);
			directMedianAdd(root, false, fixed);
			directMedianAdd(root, true, fixed);
			directMedianRejectedAdd(root);
			mergeMedian(root, merger, false, fixed);
			mergeMedian(root, merger, true, fixed);
			filterMedianNA(root, filter, fixed);
			System.out.println("COV_STATS_MEDIAN_"+(fixed ? "FIXED" : "BASELINE")+"_PASS checks="+checks);
			return;
		}
		if(inputcloseChild){
			controlledInputClose(root, "inputclose-failure-child".equals(scope));
			return;
		}
		if(inputclose){
			final Path workingDirectory=Paths.get(".").toAbsolutePath().normalize();
			final Path launcher=workingDirectory.resolve("dev/covstatslinebaselinetest.sh");
			inputCloseCase(root, launcher, workingDirectory, false, fixed);
			inputCloseCase(root, launcher, workingDirectory, true, fixed);
			System.out.println("COV_STATS_INPUTCLOSE_"+(fixed ? "FIXED" : "BASELINE")+"_PASS checks="+checks);
			return;
		}
		if(nonfinite){
			directNonfiniteCases(root, fixed);
			filterNonfinite(root, filter, false, fixed);
			filterNonfinite(root, filter, true, fixed);
			System.out.println("COV_STATS_NONFINITE_"+(fixed ? "FIXED" : "BASELINE")+"_PASS checks="+checks);
			return;
		}
		if(overflow){
			directOverflow(root, false, fixed);
			directOverflow(root, true, fixed);
			mergeOverflow(root, merger, false, fixed);
			mergeOverflow(root, merger, true, fixed);
			if(fixed){longSumContract(root);}
			System.out.println("COV_STATS_OVERFLOW_"+(fixed ? "FIXED" : "BASELINE")+"_PASS checks="+checks);
			return;
		}
		if(outputclose){
			mergeOutputClose(root, merger, fixed);
			System.out.println("COV_STATS_OUTPUTCLOSE_"+(fixed ? "FIXED" : "BASELINE")+"_PASS checks="+checks);
			return;
		}
		direct(root, false, fixed);
		direct(root, true, fixed);
		filter(root, filter, false, fixed);
		filter(root, filter, true, fixed);
		merge(root, merger, false, fixed);
		merge(root, merger, true, fixed);
		if(fixed){
			write(root.resolve("fixed-contracts.txt"), "The four fixed-contract cases are post-repair contracts, not captured baselines.\n");
			mergeContract(root, merger, false);
			mergeContract(root, merger, true);
			rejectHeader(root, merger, false);
			rejectHeader(root, merger, true);
		}
		System.out.println("COV_STATS_"+(fixed ? "FIXED" : "BASELINE")+"_PASS checks="+checks);
	}

	/** Recognizes only the six isolated child scopes; arithmetic and older modes retain their dispatch. */
	private static int inputCleanupChild(final String scope){
		for(int kind=0; kind<3; kind++){
			if(scope.equals(inputCleanupScope(kind, false))){return 2*kind;}
			if(scope.equals(inputCleanupScope(kind, true))){return 2*kind+1;}
		}
		return -1;
	}
	/** Stable names distinguish primary failure type from the independently configured close outcome. */
	private static String inputCleanupScope(final int kind, final boolean failClose){
		return "inputcleanup-"+inputCleanupName(kind)+(failClose ? "-failure-child" : "-normal-child");
	}
	private static String inputCleanupName(final int kind){
		return kind==0 ? "numeric" : kind==1 ? "header" : "aggregate";
	}

	/** Literal primary failures happen inside the original input loop, before its normal close call. */
	private static String inputCleanupData(final int kind){
		if(kind==0){return HEADER_MIN+"\n"+arithmeticRow("ctg", "bad-fold", 1, 1, 0, "0.5")+"\n";}
		if(kind==1){return HEADER_MIN.substring(1)+"\n";}
		return HEADER_MIN+"\n"+arithmeticRow("a otu", "1", 1, 4611686018427387904L, 0, "0.5")
				+"\n"+arithmeticRow("b otu", "1", 1, 4611686018427387904L, 0, "0.5")+"\n";
	}
	/** Independent class/message oracles include an Error as well as two RuntimeExceptions. */
	private static String inputCleanupPrimary(final int kind){
		if(kind==0){return "java.lang.NumberFormatException: For input string: \"bad-fold\"\n";}
		if(kind==1){return "java.lang.AssertionError: Expected a header line starting with #\n";}
		return "java.lang.ArithmeticException: long overflow\n";
	}

	/** Checks the captured leak or fixed close/suppression while retaining each original primary exactly. */
	private static void inputCleanupCase(final Path root, final Path launcher, final Path workingDirectory,
			final int kind, final boolean failClose, final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve("inputcleanup-"+inputCleanupName(kind)+(failClose ? "-failure" : "-normal")));
		final Path child=work.resolve("child");
		final String data=inputCleanupData(kind), primary=inputCleanupPrimary(kind);
		write(work.resolve("input.tsv"), data);
		write(work.resolve("expected-primary.txt"), primary);
		write(work.resolve(fixed ? "expected-fixed.txt" : "expected-baseline.txt"), "Original primary; openCount=1; "
				+(fixed ? "closeCalls>0; closeFailed="+failClose : "closeCalls=0; closeFailed=false")
				+"; mainReturned=false; output absent; suppressed count="+(fixed && failClose ? 1 : 0)
				+". Bytes/EOF are observations, not an EOF contract.\n");
		final int exit=run(launcher, Arrays.asList("phase="+(fixed ? "fixed" : "baseline"), "scope="+inputCleanupScope(kind, failClose),
				"outdir="+child), work, workingDirectory);
		check(exit==1, "Input-cleanup child must throw its original primary, not time out; exit="+exit+"; "+work);
		check(Files.size(work.resolve("stdout"))==0, "Input-cleanup child wrote stdout: "+work);
		final String observation=read(child.resolve("observation.tsv"));
		check(primary.equals(read(child.resolve("main-exception.txt"))) && observation.contains("mainReturned=false\n"),
				"Input-cleanup original primary differs: "+work);
		check(inputCloseCount(observation, "openCount")==1, "Input-cleanup stream was not opened exactly once: "+work);
		check((fixed ? inputCloseCount(observation, "closeCalls")>0 : inputCloseCount(observation, "closeCalls")==0)
				&& observation.contains("closeFailed="+(fixed && failClose)+"\n"),
				"Input-cleanup close state differs: "+work);
		check(!Files.exists(child.resolve("output.tsv")) && observation.contains("outputExists=false\n")
				&& observation.contains("outputBytes=-1\n"), "Input-cleanup primary created output: "+work);
		final String expectedSuppressed=fixed && failClose
				? "count=1\nsuppressed[0]=java.lang.RuntimeException: "+INPUT_CLOSE_REJECTION+"\n" : "count=0\n";
		check(expectedSuppressed.equals(read(child.resolve("suppressed.txt"))), "Input-cleanup suppressed failure differs: "+work);
		final String diagnostics=read(work.resolve("stderr"));
		check(diagnostics.contains(primary.trim()) && (fixed && failClose
				? diagnostics.contains("java.io.IOException: "+INPUT_CLOSE_FAILURE)
				&& diagnostics.contains("Suppressed: java.lang.RuntimeException: "+INPUT_CLOSE_REJECTION)
				: !diagnostics.contains(INPUT_CLOSE_FAILURE)),
				"Input-cleanup primary/close diagnostics differ: "+work);
		final long inputBytes=data.getBytes(StandardCharsets.US_ASCII).length;
		final long bytesRead=inputCloseCount(observation, "bytesRead");
		check(inputCloseCount(observation, "inputBytes")==inputBytes && bytesRead>0 && bytesRead<=inputBytes
				&& inputCloseCount(observation, "readCalls")>0, "Input-cleanup read counters are inconsistent: "+work);
	}

	/** Runs only in a fresh child: URL factory and probe are private to that one original-main failure. */
	private static void controlledInputCleanup(final Path root, final int kind, final boolean failClose) throws Exception{
		final String input=inputCleanupData(kind);
		final byte[] data=input.getBytes(StandardCharsets.US_ASCII);
		final InputCloseProbe probe=new InputCloseProbe(data, failClose);
		final int[] openCount={0};
		final Path output=root.resolve("output.tsv");
		write(root.resolve("input.tsv"), input);
		write(root.resolve("scope.txt"), "Controlled original-main owned-stream primary="+inputCleanupName(kind)
				+" failClose="+failClose+"; not a native file-close or read-failure reproduction.\n");
		URL.setURLStreamHandlerFactory(new URLStreamHandlerFactory(){
			@Override
			public URLStreamHandler createURLStreamHandler(final String protocol){
				if(!"jar".equals(protocol)){return null;}
				return new URLStreamHandler(){
					@Override
					protected URLConnection openConnection(final URL url) throws IOException{
						if(!INPUT_CLOSE_URL.equals(url.toExternalForm())){throw new IOException("Unexpected controlled URL: "+url);}
						openCount[0]++;
						return new URLConnection(url){
							@Override
							public void connect(){connected=true;}
							@Override
							public InputStream getInputStream(){connect(); return probe;}
						};
					}
				};
			}
		});
		Throwable primary=null;
		boolean mainReturned=false;
		try{
			MergeCoverageOTU.main(new String[]{"in="+INPUT_CLOSE_URL, "out="+output});
			mainReturned=true;
		}catch(RuntimeException | Error failure){primary=failure;}
		try{
			write(root.resolve("main-exception.txt"), primary==null ? "none\n"
					: primary.getClass().getName()+": "+primary.getMessage()+"\n");
			final Throwable[] suppressed=primary==null ? new Throwable[0] : primary.getSuppressed();
			final StringBuilder suppression=new StringBuilder("count=").append(suppressed.length).append('\n');
			for(int i=0; i<suppressed.length; i++){
				suppression.append("suppressed[").append(i).append("]=").append(suppressed[i].getClass().getName())
						.append(": ").append(suppressed[i].getMessage()).append('\n');
			}
			write(root.resolve("suppressed.txt"), suppression.toString());
			final boolean outputExists=Files.exists(output);
			write(root.resolve("observation.tsv"), "openCount="+openCount[0]+"\ninputBytes="+data.length
					+"\nbytesRead="+probe.bytesRead+"\nreadCalls="+probe.readCalls+"\neofReads="+probe.eofReads
					+"\ncloseCalls="+probe.closeCalls+"\ncloseFailed="+probe.closeFailed+"\nmainReturned="+mainReturned
					+"\noutputExists="+outputExists+"\noutputBytes="+(outputExists ? Files.size(output) : -1)+"\n");
		}catch(Exception evidenceFailure){
			if(primary==null){throw evidenceFailure;}
			if(primary!=evidenceFailure){primary.addSuppressed(evidenceFailure);}
		}
		if(primary instanceof RuntimeException){throw (RuntimeException)primary;}
		if(primary instanceof Error){throw (Error)primary;}
	}

	/** Independent signed-long boundaries; invoking reads must not mutate the parsed row. */
	private static void directReadTotal(final Path root, final boolean control, final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve(control ? "arithmetic-reads-control" : "arithmetic-reads-wrap"));
		final long strand=control ? 2305843009213693952L : 4611686018427387904L;
		final long expected=control ? 4611686018427387904L : Long.MIN_VALUE;
		final String input=arithmeticRow("ctg", "1", 1, strand, strand, "0.5");
		write(work.resolve("input.tsv"), HEADER_MIN+"\n"+input+"\n");
		write(work.resolve(fixed ? "expected-fixed.txt" : "expected-baseline.txt"),
				fixed && !control ? "ArithmeticException; receiverUnchanged=true\n" : "reads="+expected+"\nreceiverUnchanged=true\n");
		CovStatsLine.initializeHeader(HEADER_MIN);
		final CovStatsLine row=new CovStatsLine(input);
		final String before=exactAggregateState(row);
		long observed=0;
		RuntimeException failure=null;
		try{observed=row.reads();}catch(RuntimeException e){failure=e;}
		write(work.resolve("exception.txt"), failure==null ? "none\n" : failure.toString()+"\n");
		write(work.resolve("observation.tsv"), "reads="+observed+"\n");
		write(work.resolve("receiver-before.tsv"), before+"\n");
		write(work.resolve("receiver-after.tsv"), exactAggregateState(row)+"\n");
		check(fixed && !control ? failure instanceof ArithmeticException : failure==null && observed==expected,
				"Read-sum outcome differs: "+work);
		check(before.equals(exactAggregateState(row)), "reads changed receiver state: "+work);
	}

	/** Literal expectations include wrapping, denominator-only overflow and finite extreme-value controls. */
	private static ArithmeticCase[] arithmeticCases(){
		final long large=4611686018427387904L, small=2305843009213693952L;
		final String max="1.7976931348623157E308";
		return new ArithmeticCase[]{
				new ArithmeticCase("plus-wrap", arithmeticRow("a otu", "1", 1, large, 0, "0.5"),
						arithmeticRow("b otu", "1", 1, large, 0, "0.5"), Long.MIN_VALUE, 0, Long.MIN_VALUE, 1, 0x1p62),
				new ArithmeticCase("minus-wrap", arithmeticRow("a otu", "1", 1, 0, large, "0.5"),
						arithmeticRow("b otu", "1", 1, 0, large, "0.5"), 0, Long.MIN_VALUE, Long.MIN_VALUE, 1, 0x1p62),
				new ArithmeticCase("total-wrap", arithmeticRow("a otu", "1", 1, large, 0, "0.5"),
						arithmeticRow("b otu", "1", 1, 0, large, "0.5"), large, large, Long.MIN_VALUE, 1, 0x1p62),
				new ArithmeticCase("count-control", arithmeticRow("a otu", "1", 1, small, 0, "0.25"),
						arithmeticRow("b otu", "1", 1, 0, small, "0.75"), small, small, large, 1, 0.5),
				new ArithmeticCase("max-fold-overflow", arithmeticRow("a otu", max, 1, 1, 0, "0.5"),
						arithmeticRow("b otu", max, 1, 1, 0, "0.5"), 2, 0, 2, Double.POSITIVE_INFINITY, 0.5),
				new ArithmeticCase("max-fold-control", arithmeticRow("a otu", max, 1, 1, 0, "0.5"),
						arithmeticRow("b otu", "0", 1, 1, 0, "0.5"), 2, 0, 2, 0x1.fffffffffffffp1022, 0.5),
				new ArithmeticCase("ordinary-control", arithmeticRow("a otu", "1", 1, 1, 0, "0.5"),
						arithmeticRow("b otu", "3", 1, 1, 0, "0.5"), 2, 0, 2, 2, 0.5)};
	}

	/** Serializes independent test input; nonzero median/stdDev make unintended policy changes observable. */
	private static String arithmeticRow(final String id, final String fold, final long length,
			final long plus, final long minus, final String readGC){
		return id+"\t"+fold+"\t"+length+"\t0.5\t100\t"+length+"\t"+plus+"\t"+minus+"\t"+readGC+"\t7\t1.5\t0";
	}

	/** Saves raw states and checks captured arithmetic or fixed rejection/unchanged finite controls. */
	private static void directArithmetic(final Path root, final ArithmeticCase example, final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve("arithmetic-direct-"+example.name));
		writeArithmeticInput(work, example, fixed);
		CovStatsLine.initializeHeader(HEADER_MIN);
		final CovStatsLine receiver=new CovStatsLine(example.first), other=new CovStatsLine(example.second);
		final String receiverBefore=exactAggregateState(receiver), otherBefore=exactAggregateState(other);
		final String receiverTextBefore=receiver.toString(), otherTextBefore=other.toString();
		RuntimeException failure=null;
		try{receiver.add(other);}catch(RuntimeException e){failure=e;}
		write(work.resolve("receiver-before.tsv"), receiverBefore+"\n");
		write(work.resolve("receiver-after.tsv"), exactAggregateState(receiver)+"\n");
		write(work.resolve("other-before.tsv"), otherBefore+"\n");
		write(work.resolve("other-after.tsv"), exactAggregateState(other)+"\n");
		write(work.resolve("exception.txt"), failure==null ? "none\n" : failure.toString()+"\n");
		if(fixed && example.invalidAggregate){
			check(failure instanceof ArithmeticException, "Expected checked aggregate rejection: "+work);
			check(receiverBefore.equals(exactAggregateState(receiver)) && receiverTextBefore.equals(receiver.toString()),
					"Rejected add mutated receiver: "+work);
			check(otherBefore.equals(exactAggregateState(other)) && otherTextBefore.equals(other.toString()),
					"Rejected add mutated source: "+work);
			check(failure.getMessage()!=null && ("max-fold-overflow".equals(example.name)
					? failure.getMessage().contains("Avg_fold") && failure.getMessage().contains("a otu")
					&& failure.getMessage().contains("b otu") : failure.getMessage().contains("overflow")),
					"Aggregate rejection lacks expected context: "+work);
			return;
		}
		check(failure==null, "Arithmetic add baseline unexpectedly threw: "+work);
		check("a otu".equals(receiver.id) && receiver.length==2 && receiver.coveredBases==2 && receiver.underMin==0
				&& receiver.plusReads==example.plus && receiver.minusReads==example.minus && receiver.reads()==example.reads
				&& receiver.median==(fixed ? 0 : 14) && receiver.toString().contains(fixed ? "\tNA\t" : "\t14\t"),
				"Arithmetic add integral baseline differs: "+work);
		check(sameScalar(receiver.avgFold, example.fold) && sameScalar(receiver.refGC, 0.5)
				&& sameScalar(receiver.readGC, example.readGC)
				&& (fixed ? summariesUnavailable(receiver) : sameScalar(receiver.stdDev, 1.5)),
				"Arithmetic add scalar baseline differs: "+work);
		check(otherBefore.equals(exactAggregateState(other)) && otherTextBefore.equals(other.toString()),
				"Arithmetic add mutated its source: "+work);
	}

	/** Writes literal input and the independent expected outcome for the selected phase. */
	private static Path writeArithmeticInput(final Path work, final ArithmeticCase example, final boolean fixed) throws Exception{
		final Path input=work.resolve("input.tsv");
		write(input, HEADER_MIN+"\n"+example.first+"\n"+example.second+"\n");
		if(fixed && example.invalidAggregate){
			write(work.resolve("expected-rejection.txt"), "ArithmeticException; direct receiver/source unchanged; merger output absent.\n");
			return input;
		}
		write(work.resolve(fixed ? "expected-fixed.tsv" : "expected-baseline.tsv"), "length=2\tcoveredBases=2\tunderMin=0"
				+"\tmedian="+(fixed ? "NA" : "14")+"\tstdDev="+(fixed ? "NA" : "1.5")
				+"\tplusReads="+example.plus+"\tminusReads="+example.minus+"\treads="+example.reads
				+"\tavgFold="+example.fold+"\trefGC=0.5\treadGC="+example.readGC+"\n");
		return input;
	}

	/** Parses output scalar terms directly; the production constructor must not hide an emitted Infinity. */
	private static void mergeArithmetic(final Path root, final Path launcher, final ArithmeticCase example,
			final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve("arithmetic-merge-"+example.name));
		final Path input=writeArithmeticInput(work, example, fixed), output=work.resolve("output.tsv");
		final int exit=run(launcher, Arrays.asList("in="+input, "out="+output, "-Xmx128m", "-Xms32m"), work);
		if(fixed && example.invalidAggregate){
			check(exit==1, "Arithmetic merger must reject clearly; exit="+exit+"; "+work);
			check(Files.size(work.resolve("stdout"))==0, "Rejected arithmetic merger wrote stdout: "+work);
			check(!Files.exists(output), "Rejected arithmetic merger created output: "+work);
			final String diagnostics=read(work.resolve("stderr"));
			check(diagnostics.contains("ArithmeticException") && diagnostics.contains("max-fold-overflow".equals(example.name)
					? "Nonfinite aggregated coverage value" : "long overflow"), "Arithmetic merger rejection differs: "+work);
			return;
		}
		check(exit==0, "Arithmetic merger baseline exit="+exit+"; inspect "+work);
		check(Files.size(work.resolve("stdout"))==0, "Arithmetic merger wrote stdout: "+work);
		final List<String> lines=Files.readAllLines(output, StandardCharsets.US_ASCII);
		final String[] terms=lines.size()==2 ? lines.get(1).split("\t", -1) : new String[0];
		check(lines.size()==2 && HEADER_MIN.equals(lines.get(0)) && terms.length==12 && "otu".equals(terms[0])
				&& Long.parseLong(terms[2])==2 && Long.parseLong(terms[5])==2 && Long.parseLong(terms[6])==example.plus
				&& Long.parseLong(terms[7])==example.minus && (fixed ? "NA".equals(terms[9]) : Long.parseLong(terms[9])==14)
				&& Long.parseLong(terms[11])==0,
				"Arithmetic merger header/integral fields differ: "+work);
		check(sameScalar(Double.parseDouble(terms[1]), example.fold) && sameScalar(Double.parseDouble(terms[3]), 0.5)
				&& sameScalar(Double.parseDouble(terms[4]), 100) && sameScalar(Double.parseDouble(terms[8]), example.readGC)
				&& (fixed ? "NA".equals(terms[10]) : sameScalar(Double.parseDouble(terms[10]), 1.5)),
				"Arithmetic merger scalar fields differ: "+work);
	}

	/** Isolates current-row filtering from a prior-row sum used only by the optional log. */
	private static void filterReadTotal(final Path root, final Path launcher, final boolean prior,
			final boolean control, final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve("arithmetic-filter-"+(prior ? "prior-" : "current-")
				+(control ? "control" : "wrap")));
		final Path input=work.resolve("input.fa"), coverage=work.resolve("coverage.tsv"), previous=work.resolve("prior.tsv");
		final Path clean=work.resolve("clean.fa"), dirty=work.resolve("dirty.fa"), log=work.resolve("log.tsv");
		final String fasta=">ctg\nACGT\n";
		final long strand=control ? 2305843009213693952L : 4611686018427387904L;
		final String counted=arithmeticRow("ctg", "1", 4, strand, strand, "0.5");
		write(input, fasta);
		write(coverage, HEADER_MIN+"\n"+(prior ? arithmeticRow("ctg", "1", 4, 1, 1, "0.5") : counted)+"\n");
		if(prior){write(previous, HEADER_MIN+"\n"+counted+"\n");}
		final boolean rejected=!prior && !control;
		final String expectedCount=control ? "4611686018427387904" : "-9223372036854775808";
		final String logHeader="#assembly\tcontig\tcontam\tlength\tavgFold\treads\tpercentCovered"
				+(prior ? "\tavgFold0\treads0\tnormRatio" : "")+"\n";
		final String expectedLog=logHeader+"input\tctg\t"+(rejected ? "1" : "0")
				+"\t4\t1.00\t"+(prior ? "2" : expectedCount)+"\t100.00"
				+(prior ? "\t1.00\t"+expectedCount+"\t1.00" : "")+"\n";
		if(fixed && !control){
			write(work.resolve("expected-rejection.txt"), "ArithmeticException: long overflow; no success summary or data log row.\n"
					+"Output/log files may exist; this is a processing failure, not constructor preflight.\n");
		}else{
			write(work.resolve("expected-clean.fa"), rejected ? "" : fasta);
			write(work.resolve("expected-dirty.fa"), rejected ? fasta : "");
			write(work.resolve("expected-log.tsv"), expectedLog);
		}
		final ArrayList<String> arguments=new ArrayList<String>(Arrays.asList("-Xmx128m", "-Xms32m", "in="+input,
				"cov="+coverage, "out="+clean, "outdirty="+dirty, "log="+log, "logheader=t", "logappend=f",
				"basesundermin=0", "mincoverage=0", "minpercent=0", "minreads=0", "minlength=1", "trim=0",
				"fastawrap=1000", "t=1", "int=f", "ow=f", "simd=f"));
		if(prior){arguments.add("cov0="+previous);}
		final int exit=run(launcher, arguments, work);
		if(fixed && !control){
			check(exit==1, "Read-total filter must reject clearly; exit="+exit+"; "+work);
			check(Files.size(work.resolve("stdout"))==0, "Rejected read-total filter wrote stdout: "+work);
			final String diagnostics=read(work.resolve("stderr"));
			check(diagnostics.contains("ArithmeticException"), "Read-total filter lacks typed rejection: "+work);
			check(diagnostics.contains("long overflow"), "Read-total filter lacks overflow diagnostic: "+work);
			check(!diagnostics.contains("\nTime:"), "Rejected read-total filter reported success: "+work);
			final String actualLog=Files.exists(log) ? read(log) : "";
			check(actualLog.isEmpty() || logHeader.equals(actualLog), "Rejected read-total filter emitted a data log row: "+work);
			return;
		}
		check(exit==0, "Read-total filter baseline exit="+exit+"; inspect "+work);
		check(Files.size(work.resolve("stdout"))==0, "Read-total filter wrote stdout: "+work);
		check((rejected ? "" : fasta).equals(read(clean)), "Read-total filter clean output differs: "+work);
		check((rejected ? fasta : "").equals(read(dirty)), "Read-total filter dirty output differs: "+work);
		check(expectedLog.equals(read(log)), "Read-total filter literal log differs: "+work);
		final String diagnostics=read(work.resolve("stderr"));
		check(stat(diagnostics, "Reads In", 1) && stat(diagnostics, "Reads Out", rejected ? 0 : 1)
				&& stat(diagnostics, "Reads Filtered", rejected ? 1 : 0), "Read-total filter counters differ: "+work);
	}

	/** Immutable independent input/expected tuple, never initialized from production aggregation. */
	private static final class ArithmeticCase{
		ArithmeticCase(final String name, final String first, final String second, final long plus,
				final long minus, final long reads, final double fold, final double readGC){
			this.name=name; this.first=first; this.second=second;
			this.plus=plus; this.minus=minus; this.reads=reads; this.fold=fold; this.readGC=readGC;
			invalidAggregate=reads==Long.MIN_VALUE || Double.isInfinite(fold);
		}
		final String name, first, second;
		final long plus, minus, reads;
		final double fold, readGC;
		final boolean invalidAggregate;
	}

	/** Captures baseline/fixed original-main behavior with an owned URL stream in a separate JVM for each case. */
	private static void inputCloseCase(final Path root, final Path launcher, final Path workingDirectory,
			final boolean failClose, final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve(failClose ? "inputclose-failure" : "inputclose-normal"));
		final Path child=work.resolve("child");
		final boolean rejected=failClose && fixed;
		final String expectedOutput=outputCloseExpected(fixed);
		write(work.resolve("input.tsv"), inputCloseData());
		if(rejected){
			write(work.resolve("expected-rejection.txt"), INPUT_CLOSE_REJECTION+"\nOutput must remain absent.\n");
		}else{write(work.resolve("expected-output.tsv"), expectedOutput);}
		write(work.resolve("scope.txt"), "Controlled original-main owned-stream injection; not a native file-close CLI reproduction.\n");
		final int exit=run(launcher, Arrays.asList("phase="+(fixed ? "fixed" : "baseline"),
				"scope="+(failClose ? "inputclose-failure-child" : "inputclose-normal-child"), "outdir="+child),
				work, workingDirectory);
		check(exit==(rejected ? 1 : 0), "Controlled input-close child exit="+exit+"; inspect "+work);
		check(Files.size(work.resolve("stdout"))==0, "Controlled input-close child wrote stdout: "+work);
		final String observation=read(child.resolve("observation.tsv"));
		final String mainException=read(child.resolve("main-exception.txt"));
		check(rejected ? observation.contains("mainReturned=false\n")
				&& ("java.lang.RuntimeException: "+INPUT_CLOSE_REJECTION+"\n").equals(mainException)
				: observation.contains("mainReturned=true\n") && "none\n".equals(mainException),
				"Original main's controlled outcome differs: "+work);
		check(inputCloseCount(observation, "openCount")==1, "Synthetic URL was not opened exactly once: "+work);
		check(inputCloseCount(observation, "bytesRead")==inputCloseData().getBytes(StandardCharsets.US_ASCII).length
				&& inputCloseCount(observation, "readCalls")>0 && inputCloseCount(observation, "eofReads")>0,
				"Controlled input was not consumed through EOF: "+work);
		check(inputCloseCount(observation, "closeCalls")>0 && observation.contains("closeFailed="+failClose+"\n"),
				"Controlled close was not observed as configured: "+work);
		final Path output=child.resolve("output.tsv");
		check(rejected ? !Files.exists(output) && observation.contains("outputExists=false\n")
				&& observation.contains("outputBytes=-1\n")
				: Files.isRegularFile(output) && expectedOutput.equals(read(output)),
				"Controlled input-close output differs from rejection/literal oracle: "+work);
		final String diagnostics=read(work.resolve("stderr"));
		check(failClose ? diagnostics.contains("java.io.IOException: "+INPUT_CLOSE_FAILURE)
				&& (!rejected || diagnostics.contains(INPUT_CLOSE_REJECTION))
				: !diagnostics.contains(INPUT_CLOSE_FAILURE), "Controlled close diagnostic differs: "+work);
	}

	/** The URL factory is one-install global state, so this method is reachable only in an isolated child mode. */
	private static void controlledInputClose(final Path root, final boolean failClose) throws Exception{
		final byte[] data=inputCloseData().getBytes(StandardCharsets.US_ASCII);
		final InputCloseProbe probe=new InputCloseProbe(data, failClose);
		final int[] openCount={0};
		final Path output=root.resolve("output.tsv");
		write(root.resolve("input.tsv"), inputCloseData());
		write(root.resolve("scope.txt"), "Controlled original-main injection through owned "+INPUT_CLOSE_URL+"\n");
		URL.setURLStreamHandlerFactory(new URLStreamHandlerFactory(){
			@Override
			public URLStreamHandler createURLStreamHandler(final String protocol){
				if(!"jar".equals(protocol)){return null;}
				return new URLStreamHandler(){
					@Override
					protected URLConnection openConnection(final URL url) throws IOException{
						if(!INPUT_CLOSE_URL.equals(url.toExternalForm())){throw new IOException("Unexpected controlled URL: "+url);}
						openCount[0]++;
						return new URLConnection(url){
							@Override
							public void connect(){connected=true;}
							@Override
							public InputStream getInputStream(){connect(); return probe;}
						};
					}
				};
			}
		});
		Throwable mainFailure=null;
		boolean mainReturned=false;
		try{
			MergeCoverageOTU.main(new String[]{"in="+INPUT_CLOSE_URL, "out="+output});
			mainReturned=true;
		}catch(RuntimeException | Error failure){mainFailure=failure;}
		try{
			write(root.resolve("main-exception.txt"), mainFailure==null ? "none\n"
					: mainFailure.getClass().getName()+": "+mainFailure.getMessage()+"\n");
			final boolean outputExists=Files.exists(output);
			write(root.resolve("observation.tsv"), "openCount="+openCount[0]+"\ninputBytes="+data.length
					+"\nbytesRead="+probe.bytesRead+"\nreadCalls="+probe.readCalls+"\neofReads="+probe.eofReads
					+"\ncloseCalls="+probe.closeCalls+"\ncloseFailed="+probe.closeFailed+"\nmainReturned="+mainReturned
					+"\noutputExists="+outputExists+"\noutputBytes="+(outputExists ? Files.size(output) : -1)+"\n");
		}catch(Exception evidenceFailure){
			if(mainFailure==null){throw evidenceFailure;}
			if(mainFailure!=evidenceFailure){mainFailure.addSuppressed(evidenceFailure);}
		}
		if(mainFailure instanceof RuntimeException){throw (RuntimeException)mainFailure;}
		if(mainFailure instanceof Error){throw (Error)mainFailure;}
	}

	/** Independent two-row input shared with the accepted normal-output oracle. */
	private static String inputCloseData(){
		return HEADER_MIN+"\n"+OUTPUT_CLOSE_ROWS[0]+"\n"+OUTPUT_CLOSE_ROWS[1]+"\n";
	}
	/** Reads an observed nonnegative counter; missing fields fail rather than silently defaulting. */
	private static long inputCloseCount(final String observation, final String name){
		final Matcher matcher=Pattern.compile("(?m)^"+Pattern.quote(name)+"=([0-9]+)$").matcher(observation);
		if(!matcher.find()){throw new AssertionError("Missing input-close observation "+name+": "+observation);}
		return Long.parseLong(matcher.group(1));
	}

	/** Owned bytes with observable reads and a one-shot checked close failure; no borrowed stdin is involved. */
	private static final class InputCloseProbe extends ByteArrayInputStream{
		InputCloseProbe(final byte[] data, final boolean failClose){super(data); this.failClose=failClose;}
		@Override
		public synchronized int read(){
			final int value=super.read();
			readCalls++;
			if(value<0){eofReads++;}else{bytesRead++;}
			return value;
		}
		@Override
		public synchronized int read(final byte[] data, final int offset, final int length){
			final int count=super.read(data, offset, length);
			readCalls++;
			if(count<0){eofReads++;}else{bytesRead+=count;}
			return count;
		}
		@Override
		public synchronized void close() throws IOException{
			closeCalls++;
			if(!closed){
				closed=true;
				if(failClose){closeFailed=true; throw new IOException(INPUT_CLOSE_FAILURE);}
			}
		}
		private final boolean failClose;
		private boolean closed, closeFailed;
		private long bytesRead, readCalls, eofReads, closeCalls;
	}

	/** Thirteen nonfinite cases and two finite controls; every expectation is independent of production parsing. */
	private static void directNonfiniteCases(final Path root, final boolean fixed) throws Exception{
		final String[] ordinary=NONFINITE_ROW.split("\t", -1);
		final double[] finite={10, 0.4, 0.25, 1.5};
		final String[] special={"NaN", "+Infinity", "-Infinity"};
		final String[] names={"nan", "positive-infinity", "negative-infinity"};
		final double[] expectedSpecial={Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY};
		for(int field=0; field<SCALAR_COLUMNS.length; field++){
			for(int kind=0; kind<special.length; kind++){
				final String[] columns=ordinary.clone();
				final double[] expected=finite.clone();
				columns[SCALAR_COLUMNS[field]]=special[kind];
				expected[field]=expectedSpecial[kind];
				directNonfinite(root, SCALAR_CASE_NAMES[field]+"-"+names[kind], columns, expected, fixed);
			}
		}
		final String[] exponent=ordinary.clone();
		exponent[1]="1e309";
		directNonfinite(root, "avg-fold-exponent-overflow", exponent,
				new double[]{Double.POSITIVE_INFINITY, 0.4, 0.25, 1.5}, fixed);
		directNonfinite(root, "ordinary-finite", ordinary, finite, fixed);
		final String[] policy=ordinary.clone();
		policy[1]="-0.0"; policy[3]="1.25"; policy[8]="-0.5"; policy[10]="-1";
		//Preserves current finite acceptance, not an endorsement of biological ranges or sentinels for every field.
		directNonfinite(root, "finite-policy", policy, new double[]{-0.0, 1.25, -0.5, -1}, fixed);
	}

	/** Captures baseline acceptance or fixed rejection, preserving all finite control fields and raw bits. */
	private static void directNonfinite(final Path root, final String name, final String[] columns,
			final double[] expected, final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve("nonfinite-direct-"+name));
		final String input=String.join("\t", columns);
		write(work.resolve("input.tsv"), HEADER_WINDOW+"\n"+input+"\n");
		CovStatsLine.initializeHeader(HEADER_WINDOW);
		CovStatsLine row=null;
		RuntimeException failure=null;
		try{row=new CovStatsLine(input);}
		catch(RuntimeException e){failure=e;}
		write(work.resolve("exception.txt"), failure==null ? "none\n"
				: failure.getClass().getName()+": "+failure.getMessage()+"\n");
		if(row!=null){
			write(work.resolve("observation.tsv"), overflowDetails(row)+"\treadGC="+row.readGC
					+"\tmedian="+row.median+"\tstdDev="+row.stdDev+"\n");
			write(work.resolve("scalar-bits.tsv"), exactAggregateState(row)+"\n");
		}
		int nonfiniteField=-1;
		for(int i=0; i<expected.length; i++){
			if(Double.isNaN(expected[i]) || Double.isInfinite(expected[i])){nonfiniteField=i; break;}
		}
		if(fixed && nonfiniteField>=0){
			check(failure instanceof NumberFormatException, "Expected nonfinite constructor rejection: "+work);
			final String message=failure.getMessage();
			check(message!=null && message.contains("ctg") && message.contains(SCALAR_HEADERS[nonfiniteField])
					&& message.contains(columns[SCALAR_COLUMNS[nonfiniteField]]),
					"Nonfinite rejection lacks record/column/original token: "+work+"; "+message);
			return;
		}
		if(failure!=null){throw new AssertionError("Unexpected constructor rejection: "+work, failure);}
		check(sameScalar(row.avgFold, expected[0]) && sameScalar(row.refGC, expected[1])
				&& sameScalar(row.readGC, expected[2]) && sameScalar(row.stdDev, expected[3]),
				"Constructor changed scalar values/classes/signs: "+work);
		check("ctg".equals(row.id) && row.length==1000 && row.coveredBases==900 && row.coveredPercent()==90
				&& row.plusReads==7 && row.minusReads==3 && row.reads()==10 && row.median==8 && row.underMin==0,
				"Constructor changed other coverage fields: "+work);
	}

	/** NaN has no payload requirement here; all other values include raw signed-zero/infinity bits. */
	private static boolean sameScalar(final double actual, final double expected){
		return Double.isNaN(expected) ? Double.isNaN(actual)
				: Double.doubleToRawLongBits(actual)==Double.doubleToRawLongBits(expected);
	}

	/** Captures NaN acceptance or rejection before outputs open, with the same finite filter/log control. */
	private static void filterNonfinite(final Path root, final Path launcher, final boolean finite,
			final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve(finite ? "nonfinite-filter-control" : "nonfinite-filter-nan"));
		final Path input=work.resolve("input.fa"), coverage=work.resolve("coverage.tsv");
		final Path clean=work.resolve("clean.fa"), dirty=work.resolve("dirty.fa"), log=work.resolve("log.tsv");
		final StringBuilder sequence=new StringBuilder(1000);
		for(int i=0; i<250; i++){sequence.append("ACGT");}
		final String fasta=">ctg\n"+sequence+"\n";
		final String[] columns=NONFINITE_ROW.split("\t", -1);
		if(!finite){columns[1]="NaN";}
		write(input, fasta);
		write(coverage, HEADER_WINDOW+"\n"+String.join("\t", columns)+"\n");
		final String expectedLog="#assembly\tcontig\tcontam\tlength\tavgFold\treads\tpercentCovered\n"
				+"input\tctg\t0\t1000\t"+(finite ? "10.00" : "NaN")+"\t10\t90.00\n";
		if(fixed && !finite){
			write(work.resolve("expected-rejection.txt"), "NumberFormatException naming ctg, Avg_fold and NaN; clean/dirty/log absent.\n");
		}else{
			write(work.resolve("expected-clean.fa"), fasta);
			write(work.resolve("expected-dirty.fa"), "");
			write(work.resolve("expected-log.tsv"), expectedLog);
		}
		final int exit=run(launcher, Arrays.asList("-Xmx128m", "-Xms32m", "in="+input, "cov="+coverage,
				"out="+clean, "outdirty="+dirty, "log="+log, "logheader=t", "logappend=f", "basesundermin=-1",
				"mincoverage=5", "minpercent=40", "minreads=0", "minlength=1", "trim=0", "fastawrap=1000",
				"t=1", "int=f", "ow=f", "simd=f"), work);
		if(fixed && !finite){
			check(exit==1, "Nonfinite filter must reject clearly; exit="+exit+"; "+work);
			check(Files.size(work.resolve("stdout"))==0, "Rejected nonfinite filter wrote stdout: "+work);
			check(!Files.exists(clean), "Rejected nonfinite filter created clean output: "+work);
			check(!Files.exists(dirty), "Rejected nonfinite filter created dirty output: "+work);
			check(!Files.exists(log), "Rejected nonfinite filter created log output: "+work);
			final String diagnostics=read(work.resolve("stderr"));
			check(diagnostics.contains("NumberFormatException") && diagnostics.contains("ctg")
					&& diagnostics.contains("Avg_fold") && diagnostics.contains("NaN") && !diagnostics.contains("\nTime:"),
					"Nonfinite filter lacks early descriptive rejection: "+work);
			return;
		}
		check(exit==0, "Nonfinite filter baseline exit="+exit+"; inspect "+work);
		check(Files.size(work.resolve("stdout"))==0, "Nonfinite filter wrote stdout: "+work);
		check(fasta.equals(read(clean)), "Nonfinite baseline clean output differs: "+work);
		check("".equals(read(dirty)), "Nonfinite baseline dirty output differs: "+work);
		check(expectedLog.equals(read(log)), "Nonfinite baseline exact log differs: "+work);
		final String diagnostics=read(work.resolve("stderr"));
		check(stat(diagnostics, "Reads In", 1) && stat(diagnostics, "Reads Out", 1)
				&& stat(diagnostics, "Reads Filtered", 0), "Nonfinite baseline filter counters differ: "+work);
	}

	/** Saves only two tiny summary rows, including a separate correct-large oracle that is not a captured result. */
	private static Path overflowInput(final Path work, final boolean control, final boolean fixed) throws Exception{
		final String[] rows=control ? OVERFLOW_CONTROL_ROWS : OVERFLOW_LARGE_ROWS;
		final Path input=work.resolve("input.tsv");
		write(input, OVERFLOW_HEADER+"\n"+rows[0]+"\n"+rows[1]+"\n");
		write(work.resolve(fixed ? "expected-merged-fixed.tsv" : "expected-merged-baseline.tsv"), OVERFLOW_HEADER+"\n"
				+(control ? OVERFLOW_CONTROL_OUTPUT : fixed ? OVERFLOW_CORRECT_OUTPUT : OVERFLOW_BASELINE_OUTPUT)+"\n");
		if(!control){write(work.resolve("expected-correct-large.tsv"), OVERFLOW_HEADER+"\n"+OVERFLOW_CORRECT_OUTPUT+"\n");}
		return input;
	}

	/** Checks wrapped baseline or exact repaired totals, with an exact small power-of-two control. */
	private static void directOverflow(final Path root, final boolean control, final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve(control ? "overflow-direct-control" : "overflow-direct-large"));
		overflowInput(work, control, fixed);
		final String[] rows=control ? OVERFLOW_CONTROL_ROWS : OVERFLOW_LARGE_ROWS;
		CovStatsLine.initializeHeader(OVERFLOW_HEADER);
		final CovStatsLine first=new CovStatsLine(rows[0]), second=new CovStatsLine(rows[1]);
		first.add(second);
		write(work.resolve("observation.tsv"), overflowDetails(first)+"\n");
		check(overflowCountsMatch(first, control, fixed, "a otu"), "Direct overflow counts differ: "+work);
		check(overflowMetricsMatch(first, control, fixed), "Direct overflow weighted metrics differ: "+work);
	}

	/** Uses the real legacy merger for both two-row shapes; exact files and reparsed values are independent checks. */
	private static void mergeOverflow(final Path root, final Path launcher, final boolean control, final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve(control ? "overflow-merge-control" : "overflow-merge-large"));
		final Path input=overflowInput(work, control, fixed), output=work.resolve("output.tsv");
		final int exit=run(launcher, Arrays.asList("in="+input, "out="+output, "-Xmx128m", "-Xms32m"), work);
		check(exit==0, "Overflow merger exit="+exit+"; inspect "+work);
		check(Files.size(work.resolve("stdout"))==0, "Overflow merger wrote stdout: "+work);
		check(read(work.resolve(fixed ? "expected-merged-fixed.tsv" : "expected-merged-baseline.tsv")).equals(read(output)),
				"Overflow merger exact output differs: "+work);
		final List<String> lines=Files.readAllLines(output, StandardCharsets.US_ASCII);
		CovStatsLine.initializeHeader(lines.get(0));
		final CovStatsLine merged=new CovStatsLine(lines.get(1));
		write(work.resolve("reparse-observation.tsv"), overflowDetails(merged)+"\n");
		check(overflowCountsMatch(merged, control, fixed, "otu") && overflowMetricsMatch(merged, control, fixed),
				"Overflow merger reparse differs: "+work);
	}

	/** Literal baseline wraps and fixed positive totals are independent of production add's implementation. */
	private static boolean overflowCountsMatch(final CovStatsLine row, final boolean control, final boolean fixed, final String id){
		return id.equals(row.id) && row.length==(control ? 2048L : fixed ? 3000000000L : -1294967296L)
				&& row.coveredBases==(control ? 1920L : fixed ? 2800000000L : -1494967296L)
				&& row.underMin==(control ? 1536L : fixed ? 2400000000L : -1894967296L)
				&& row.plusReads==2 && row.minusReads==0 && row.reads()==2;
	}
	private static boolean overflowMetricsMatch(final CovStatsLine row, final boolean control, final boolean fixed){
		if(fixed && !control){
			//Inputs are below2^53; four ULPs cover reciprocal/multiply rounding without hiding count errors.
			//Integral totals and the independent formatted CLI oracle are still required to match exactly.
			return withinFourUlps(row.avgFold, 1d) && withinFourUlps(row.refGC, 0.4)
					&& withinFourUlps(row.coveredPercent(), 93.33333333333333);
		}
		return row.avgFold==(control ? 1d : 3000000000d) && row.refGC==(control ? 0.4 : 1200000000d)
				&& row.coveredPercent()==(control ? 93.75 : -149496729600d);
	}
	private static boolean withinFourUlps(final double actual, final double expected){
		return Math.abs(actual-expected)<=4*Math.ulp(expected);
	}

	/** New post-repair contract: an unrepresentable long sum fails before changing any receiver field. */
	private static void longSumContract(final Path root) throws Exception{
		final Path work=Files.createDirectory(root.resolve("fixed-long-sum-contract"));
		write(work.resolve("contract.txt"), "Post-repair contract, not a captured baseline: Long.MAX_VALUE+1 must preserve receiver.\n");
		final String first="limit otu\t1.25\t9223372036854775807\t0.4\t0\t1\t2\t3\t1";
		final String second="extra otu\t2\t1\t0.6\t100\t1\t4\t5\t1";
		write(work.resolve("input.tsv"), OVERFLOW_HEADER+"\n"+first+"\n"+second+"\n");
		CovStatsLine.initializeHeader(OVERFLOW_HEADER);
		final CovStatsLine receiver=new CovStatsLine(first), other=new CovStatsLine(second);
		//Nonzero legacy fields make premature mean/read-count/summary mutation observable too.
		receiver.readGC=0.25; receiver.median=7; receiver.stdDev=1.5;
		final String before=exactAggregateState(receiver);
		ArithmeticException failure=null;
		try{receiver.add(other);}
		catch(ArithmeticException e){failure=e;}
		final String after=exactAggregateState(receiver);
		write(work.resolve("receiver-before.tsv"), before+"\n");
		write(work.resolve("receiver-after.tsv"), after+"\n");
		write(work.resolve("exception.txt"), failure==null ? "none\n" : failure.getClass().getName()+": "+failure.getMessage()+"\n");
		check(failure!=null, "Long.MAX_VALUE+1 did not fail: "+work);
		check(before.equals(after), "Checked-count failure changed receiver state: "+work);
	}

	/** Raw double bits plus all integral/string values detect even a one-ULP premature receiver mutation. */
	private static String exactAggregateState(final CovStatsLine row){
		return "id="+row.id+"\tlength="+row.length+"\tcoveredBases="+row.coveredBases+"\tunderMin="+row.underMin
				+"\tplusReads="+row.plusReads+"\tminusReads="+row.minusReads+"\tmedian="+row.median
				+"\tmedianAvailable="+row.medianAvailable+"\tmedianText="+row.medianText()
				+"\tavgFoldBits="+Double.doubleToRawLongBits(row.avgFold)+"\trefGCBits="+Double.doubleToRawLongBits(row.refGC)
				+"\treadGCBits="+Double.doubleToRawLongBits(row.readGC)+"\tstdDevBits="+Double.doubleToRawLongBits(row.stdDev)
				+"\tstdDevAvailable="+row.stdDevAvailable+"\tstdDevText="+row.stdDevText();
	}
	private static String overflowDetails(final CovStatsLine row){
		return "id="+row.id+"\tlength="+row.length+"\tcoveredBases="+row.coveredBases+"\tunderMin="+row.underMin
				+"\tavgFold="+row.avgFold+"\trefGC="+row.refGC+"\tcoveredPercent="+row.coveredPercent()
				+"\tplusReads="+row.plusReads+"\tminusReads="+row.minusReads;
	}

	/** Reuses the real merger to capture success on a literal file and ignored /dev/full close failure. */
	private static void mergeOutputClose(final Path root, final Path launcher, final boolean fixed) throws Exception{
		final Path devFull=Paths.get("/dev/full");
		final int devFullMode=((Number)Files.getAttribute(devFull, "unix:mode")).intValue();
		check((devFullMode&0170000)==0020000,
				"scope=outputclose requires the existing /dev/full character device.");
		mergeOutputCloseNormal(Files.createDirectory(root.resolve("outputclose-normal")), launcher, fixed);
		mergeOutputCloseFailure(Files.createDirectory(root.resolve("outputclose-dev-full")), launcher, devFull, fixed);
	}

	/** The two-row control retains all fields except the approved unavailable merged summaries. */
	private static void mergeOutputCloseNormal(final Path work, final Path launcher, final boolean fixed) throws Exception{
		final Path input=outputCloseInput(work), output=work.resolve("output.tsv");
		final String expected=outputCloseExpected(fixed);
		write(work.resolve("expected-output.tsv"), expected);
		final int exit=run(launcher, Arrays.asList("in="+input, "out="+output, "-Xmx128m", "-Xms32m"), work);
		check(exit==0, "Output-close normal merger exit="+exit+"; inspect "+work);
		check(Files.size(work.resolve("stdout"))==0, "Output-close normal merger wrote stdout: "+work);
		check(expected.equals(read(output)), "Output-close normal merger exact output differs: "+work);
	}

	/** Captures the ignored write/close error baseline and the repaired error-state failure. */
	private static void mergeOutputCloseFailure(final Path work, final Path launcher, final Path devFull,
			final boolean fixed) throws Exception{
		final Path input=outputCloseInput(work);
		final String expectedDiagnostic="driver.MergeCoverageOTU terminated in an error state; the output may be corrupt.";
		write(work.resolve(fixed ? "expected-fixed.txt" : "expected-baseline.txt"), fixed ? expectedDiagnostic+"\n"
				: "The current merger ignores TextStreamWriter.poisonAndWait(), so /dev/full still exits0.\n");
		final int exit=run(launcher, Arrays.asList("in="+input, "out="+devFull, "-Xmx128m", "-Xms32m"), work);
		check(exit==(fixed ? 1 : 0), "/dev/full output-close exit="+exit+"; inspect "+work);
		check(Files.size(work.resolve("stdout"))==0, "/dev/full output-close wrote stdout: "+work);
		if(fixed){
			check(read(work.resolve("stderr")).contains(expectedDiagnostic),
					"/dev/full fixed output-close lacks stable error-state diagnostic: "+work);
		}
	}

	/** Two contigs in one OTU exercise aggregation before output closure. */
	private static Path outputCloseInput(final Path work) throws Exception{
		final Path input=work.resolve("input.tsv");
		write(input, HEADER_MIN+"\n"+OUTPUT_CLOSE_ROWS[0]+"\n"+OUTPUT_CLOSE_ROWS[1]+"\n");
		return input;
	}

	/** Numeric singleton medians remain directly representable. */
	private static void directMedianNumeric(final Path root) throws Exception{
		final Path work=Files.createDirectory(root.resolve("median-direct-numeric"));
		write(work.resolve("input.tsv"), HEADER_MIN+"\n"+MEDIAN_FIRST+"\n");
		CovStatsLine.initializeHeader(HEADER_MIN);
		final CovStatsLine row=new CovStatsLine(MEDIAN_FIRST);
		write(work.resolve("receiver.tsv"), row.toString()+"\n");
		check("ctgA otuA".equals(row.id) && row.median==8 && medianSingletonOutput(false).equals(row.toString()),
				"Numeric median singleton parsed incorrectly: "+work);
		check(row.length==1000 && row.plusReads==7 && row.minusReads==3 && row.underMin==600
				&& sameScalar(row.readGC, 0.25) && sameScalar(row.stdDev, 1.5)
				&& row.medianAvailable && "8".equals(row.medianText())
				&& row.stdDevAvailable && "1.5000".equals(row.stdDevText()),
				"Numeric median singleton changed nonmedian fields: "+work);
	}

	/** Missing pooled medians are literal NA in the repaired row format and rejected by the baseline. */
	private static void directMedianNA(final Path root, final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve("median-direct-na"));
		final String input=medianNA(MEDIAN_FIRST);
		write(work.resolve("input.tsv"), HEADER_MIN+"\n"+input+"\n");
		CovStatsLine.initializeHeader(HEADER_MIN);
		NumberFormatException failure=null;
		CovStatsLine row=null;
		try{row=new CovStatsLine(input);}
		catch(NumberFormatException e){failure=e;}
		write(work.resolve("exception.txt"), failure==null ? "none\n" : failure.toString()+"\n");
		if(failure==null){write(work.resolve("receiver.tsv"), row.toString()+"\n");}
		if(fixed){
			check(failure==null, "Fixed NA median was rejected: "+work);
			check(medianSingletonOutput(true).equals(row.toString()) && !row.medianAvailable
					&& row.median==0 && "NA".equals(row.medianText())
					&& row.stdDevAvailable && "1.5000".equals(row.stdDevText()),
					"Fixed NA median changed availability or singleton deviation: "+work);
		}else{
			check(failure!=null, "Baseline NA median was accepted unexpectedly: "+work);
			check(failure.getMessage()!=null && failure.getMessage().contains("NA"),
					"Baseline NA median rejection lacks the original token: "+work);
		}
	}

	/** Every successful add makes both summaries unavailable while leaving the source row unchanged. */
	private static void directMedianAdd(final Path root, final boolean zero, final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve(zero ? "median-add-zero" : "median-add-positive"));
		final String first=zero ? MEDIAN_ZERO_FIRST : MEDIAN_FIRST;
		final String second=zero ? MEDIAN_ZERO_SECOND : MEDIAN_SECOND;
		write(work.resolve("input.tsv"), HEADER_MIN+"\n"+first+"\n"+second+"\n");
		CovStatsLine.initializeHeader(HEADER_MIN);
		final CovStatsLine receiver=new CovStatsLine(first), other=new CovStatsLine(second);
		final String otherBefore=exactAggregateState(other);
		final String otherTextBefore=other.toString();
		RuntimeException failure=null;
		try{receiver.add(other);}
		catch(RuntimeException e){failure=e;}
		write(work.resolve("exception.txt"), failure==null ? "none\n" : failure.toString()+"\n");
		write(work.resolve("receiver.tsv"), receiver.toString()+"\n");
		check(failure==null, "Median add unexpectedly failed: "+work);
		check(otherBefore.equals(exactAggregateState(other)) && otherTextBefore.equals(other.toString()),
				"Median add mutated its source row: "+work);
		check(medianAddOutput(zero, fixed).equals(receiver.toString()) && (!fixed || summariesUnavailable(receiver)),
				"Median/deviation add output or availability differs: "+work);
	}

	/** A rejected add keeps the receiver's singleton median available through serialization. */
	private static void directMedianRejectedAdd(final Path root) throws Exception{
		final Path work=Files.createDirectory(root.resolve("median-rejected-add"));
		final String first="limit otu\t1.25\t9223372036854775807\t0.4\t0\t1\t2\t3\t0.25\t8\t1.5\t0";
		final String second="extra otu\t2\t1\t0.6\t100\t1\t4\t5\t0.75\t9\t2.5\t1";
		final String expected="limit otu\t1.2500\t9223372036854775807\t0.4000\t0.0000\t1\t2\t3\t8\t0\t0.2500\t1.5000";
		write(work.resolve("input.tsv"), HEADER_MIN+"\n"+first+"\n"+second+"\n");
		CovStatsLine.initializeHeader(HEADER_MIN);
		final CovStatsLine receiver=new CovStatsLine(first), other=new CovStatsLine(second);
		final String before=exactAggregateState(receiver);
		ArithmeticException failure=null;
		try{receiver.add(other);}
		catch(ArithmeticException e){failure=e;}
		write(work.resolve("exception.txt"), failure==null ? "none\n" : failure.toString()+"\n");
		write(work.resolve("receiver.tsv"), receiver.toString()+"\n");
		check(failure!=null, "Rejected median add did not fail: "+work);
		check(before.equals(exactAggregateState(receiver)), "Rejected median add mutated receiver: "+work);
		check(expected.equals(receiver.toString()) && receiver.medianAvailable && "8".equals(receiver.medianText())
				&& receiver.stdDevAvailable && "1.5000".equals(receiver.stdDevText()),
				"Rejected add lost numeric singleton summaries: "+work);
	}

	/** Real merger output sums current medians or emits NA once aggregation marks them unavailable. */
	private static void mergeMedian(final Path root, final Path launcher, final boolean singleton,
			final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve(singleton ? "median-merge-singleton" : "median-merge-positive"));
		final Path input=work.resolve("input.tsv"), output=work.resolve("output.tsv");
		final String rows=singleton ? MEDIAN_FIRST+"\n" : MEDIAN_FIRST+"\n"+MEDIAN_SECOND+"\n";
		final String expected=HEADER_MIN+"\n"+(singleton ? MEDIAN_SINGLETON_OUTPUT : medianMergeOutput(fixed))+"\n";
		write(input, HEADER_MIN+"\n"+rows);
		write(work.resolve(fixed ? "expected-fixed.tsv" : "expected-baseline.tsv"), expected);
		final int exit=run(launcher, Arrays.asList("in="+input, "out="+output, "-Xmx128m", "-Xms32m"), work);
		check(exit==0, "Median merger exit="+exit+"; inspect "+work);
		check(Files.size(work.resolve("stdout"))==0, "Median merger wrote stdout: "+work);
		check(expected.equals(read(output)), "Median merger exact output differs: "+work);
	}

	/** FilterByCoverage must accept an unavailable median while retaining other finite filters. */
	private static void filterMedianNA(final Path root, final Path launcher, final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve("median-filter-na"));
		final Path input=work.resolve("input.fa"), coverage=work.resolve("coverage.tsv");
		final Path clean=work.resolve("clean.fa"), dirty=work.resolve("dirty.fa");
		final StringBuilder sequence=new StringBuilder(1000);
		for(int i=0; i<250; i++){sequence.append("ACGT");}
		final String fasta=">ctgA\n"+sequence+"\n";
		write(input, fasta);
		write(coverage, HEADER_MIN+"\n"+medianNA(MEDIAN_FIRST.replace("ctgA otuA", "ctgA"))+"\n");
		final int exit=run(launcher, Arrays.asList("-Xmx128m", "-Xms32m", "in="+input, "cov="+coverage,
				"out="+clean, "outdirty="+dirty, "basesundermin=-1", "mincoverage=5", "minpercent=40",
				"minreads=0", "minlength=1", "trim=0", "fastawrap=1000", "t=1", "int=f", "ow=f", "simd=f"), work);
		if(fixed){
			check(exit==0, "Fixed NA median filter failed: exit="+exit+"; "+work);
			check(Files.size(work.resolve("stdout"))==0, "Fixed NA median filter wrote stdout: "+work);
			check(fasta.equals(read(clean)), "Fixed NA median filter clean output differs: "+work);
			check("".equals(read(dirty)), "Fixed NA median filter dirty output differs: "+work);
			final String diagnostics=read(work.resolve("stderr"));
			check(stat(diagnostics, "Reads In", 1) && stat(diagnostics, "Reads Out", 1)
					&& stat(diagnostics, "Reads Filtered", 0), "Fixed NA median filter counters differ: "+work);
		}else{
			check(exit==1, "Baseline NA median filter did not reject: exit="+exit+"; "+work);
			check(Files.size(work.resolve("stdout"))==0, "Baseline NA median filter wrote stdout: "+work);
			check(!Files.exists(clean), "Baseline NA median filter created clean output: "+work);
			check(!Files.exists(dirty), "Baseline NA median filter created dirty output: "+work);
			final String diagnostics=read(work.resolve("stderr"));
			check(diagnostics.contains("NumberFormatException") && diagnostics.contains("NA")
					&& !diagnostics.contains("\nTime:"), "Baseline NA median filter lacks early rejection: "+work);
		}
	}
	private static String medianNA(final String row){
		return row.replace("\t8\t", "\tNA\t");
	}
	private static String medianSingletonOutput(final boolean na){
		return "ctgA otuA\t10.0000\t1000\t0.4000\t90.0000\t900\t7\t3\t"+(na ? "NA" : "8")
				+"\t600\t0.2500\t1.5000";
	}
	private static String medianAddOutput(final boolean zero, final boolean fixed){
		return zero ? "zeroA otuA\t0.0000\t0\t0.0000\t0.0000\t0\t0\t0\t"+(fixed ? "NA" : "16")
				+"\t0\t0.0000\t"+(fixed ? "NA" : "1.5000") : "ctgA otuA\t10.0000\t2000\t0.4000\t90.0000\t1800\t14\t6\t"
				+(fixed ? "NA" : "16")+"\t1200\t0.2500\t"+(fixed ? "NA" : "1.5000");
	}
	private static String medianMergeOutput(final boolean fixed){
		return "otuA\t10.0000\t2000\t0.4000\t90.0000\t1800\t14\t6\t0.2500\t"
				+(fixed ? "NA" : "16")+"\t"+(fixed ? "NA" : "1.5000")+"\t1200";
	}

	/** Approved aggregation loses both summaries; zero storage alone cannot distinguish known zero. */
	private static boolean summariesUnavailable(final CovStatsLine row){
		return !row.medianAvailable && row.median==0 && "NA".equals(row.medianText())
				&& !row.stdDevAvailable && sameScalar(row.stdDev, 0d) && "NA".equals(row.stdDevText());
	}
	/** Keeps historical close baselines separate from the approved fixed pooled-summary policy. */
	private static String outputCloseExpected(final boolean fixed){
		return fixed ? HEADER_MIN+"\n"
				+"otuA\t17.5000\t40\t0.3125\t62.5000\t25\t6\t4\t0.5000\tNA\tNA\t3\n" : OUTPUT_CLOSE_EXPECTED;
	}

	/** Checks Under_* normalization against the already-recognized Under_min control. */
	private static void direct(final Path root, final boolean control, final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve(control ? "direct-control" : "direct-under"));
		final String header=control ? HEADER_MIN : HEADER_WINDOW;
		write(work.resolve("input.tsv"), header+"\n"+ROW+"\n");
		CovStatsLine.initializeHeader(header);
		final CovStatsLine row=new CovStatsLine(ROW);
		write(work.resolve("observation.tsv"), "id="+row.id+"\tunderMin="+row.underMin
				+"\tavgFold="+row.avgFold+"\tlength="+row.length+"\tcoveredBases="+row.coveredBases
				+"\tplusReads="+row.plusReads+"\tminusReads="+row.minusReads+"\trefGC="+row.refGC
				+"\treadGC="+row.readGC+"\tmedian="+row.median+"\tstdDev="+row.stdDev+"\n");
		check(row.underMin==(control || fixed ? 600 : 0), "Unexpected Under value; inspect "+work);
		check("ctg".equals(row.id) && metricsMatch(row), "Unchanged direct metrics differ; inspect "+work);
	}

	/** Captures the wrong acceptance and repaired rejection of the same one-thousand-base contig. */
	private static void filter(final Path root, final Path launcher, final boolean control, final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve(control ? "filter-control" : "filter-under"));
		final StringBuilder sequence=new StringBuilder(1000);
		for(int i=0; i<250; i++){sequence.append("ACGT");}
		final String fasta=">ctg\n"+sequence+"\n";
		final Path input=work.resolve("input.fa"), coverage=work.resolve("coverage.tsv");
		final Path clean=work.resolve("clean.fa"), dirty=work.resolve("dirty.fa");
		write(input, fasta);
		write(coverage, (control ? HEADER_MIN : HEADER_WINDOW)+"\n"+ROW+"\n");
		final boolean rejected=control || fixed;
		final String expectedClean=rejected ? "" : fasta, expectedDirty=rejected ? fasta : "";
		write(work.resolve("expected-clean.fa"), expectedClean);
		write(work.resolve("expected-dirty.fa"), expectedDirty);
		final int exit=run(launcher, Arrays.asList("-Xmx128m", "-Xms32m", "in="+input,
				"cov="+coverage, "out="+clean, "outdirty="+dirty, "basesundermin=500",
				"mincoverage=5", "minpercent=40", "minreads=0", "minlength=1", "trim=0",
				"fastawrap=1000", "t=1", "int=f", "ow=f", "simd=f"), work);
		check(exit==0, "Filter exit="+exit+"; inspect "+work);
		check(Files.size(work.resolve("stdout"))==0, "Filter wrote unexpected stdout: "+work);
		check(expectedClean.equals(read(clean)), "Filter clean output differs: "+work);
		check(expectedDirty.equals(read(dirty)), "Filter dirty output differs: "+work);
		final String diagnostics=read(work.resolve("stderr"));
		check(stat(diagnostics, "Reads In", 1) && stat(diagnostics, "Reads Out", rejected ? 0 : 1)
				&& stat(diagnostics, "Reads Filtered", rejected ? 1 : 0), "Filter counters differ: "+work);
	}

	/** A single row avoids aggregation policy; the recognized Under_min spelling isolates serialization. */
	private static void merge(final Path root, final Path launcher, final boolean control, final boolean fixed) throws Exception{
		final Path work=Files.createDirectory(root.resolve(control ? "merge-control" : "merge-order"));
		final String header=control ? HEADER_COMPATIBLE : HEADER_MIN;
		final String inputRow=control ? MERGE_COMPATIBLE : MERGE_INPUT;
		final Path input=work.resolve("input.tsv"), output=work.resolve("output.tsv");
		write(input, header+"\n"+inputRow+"\n");
		final String expected=header+"\n"+(control || !fixed ? SERIALIZED : HEADER_PRESERVING)+"\n";
		write(work.resolve(fixed ? "expected-fixed.tsv" : "expected-baseline.tsv"), expected);
		//Keep the independent meaning oracle alongside the literal baseline/fixed expectations.
		write(work.resolve("expected-header-preserving.tsv"), header+"\n"
				+(control ? SERIALIZED : HEADER_PRESERVING)+"\n");
		//The legacy main reads only args[0]/args[1]; JVM launcher options must follow these two.
		final int exit=run(launcher, Arrays.asList("in="+input, "out="+output, "-Xmx128m", "-Xms32m"), work);
		check(exit==0, "Merger exit="+exit+"; inspect "+work);
		check(Files.size(work.resolve("stdout"))==0, "Merger wrote unexpected stdout: "+work);
		final String actual=read(output);
		check(expected.equals(actual), "Merger exact header/row differs: "+work);
		final List<String> lines=Files.readAllLines(output, StandardCharsets.US_ASCII);
		CovStatsLine.initializeHeader(lines.get(0));
		boolean parsed=false;
		String failure="none";
		try{
			final CovStatsLine row=new CovStatsLine(lines.get(1));
			parsed="otuA".equals(row.id) && row.underMin==600 && metricsMatch(row);
		}catch(NumberFormatException e){
			failure=e.getClass().getName()+": "+e.getMessage();
		}
		write(work.resolve("reparse-observation.txt"), "metricsMatch="+parsed+"\nexception="+failure+"\n");
		check(control || fixed ? parsed : !parsed && failure.contains("NumberFormatException") && failure.contains("1.5000"),
				"Unexpected merger reparse result: "+work+"; "+failure);
	}

	/** Fixed-only contracts: preserve an actual Under threshold/window label and omit absent optional fields. */
	private static void mergeContract(final Path root, final Path launcher, final boolean coreOnly) throws Exception{
		final Path work=Files.createDirectory(root.resolve(coreOnly ? "fixed-core-contract" : "fixed-under-contract"));
		final String header=coreOnly ? CORE : HEADER_WINDOW;
		final String inputRow=coreOnly ? MERGE_CORE_INPUT : MERGE_INPUT;
		final String expected=header+"\n"+(coreOnly ? CORE_OUTPUT : HEADER_PRESERVING)+"\n";
		final Path input=work.resolve("input.tsv"), output=work.resolve("output.tsv");
		write(input, header+"\n"+inputRow+"\n");
		write(work.resolve("expected-fixed.tsv"), expected);
		final int exit=run(launcher, Arrays.asList("in="+input, "out="+output, "-Xmx128m", "-Xms32m"), work);
		check(exit==0, "Fixed merger contract exit="+exit+"; inspect "+work);
		check(Files.size(work.resolve("stdout"))==0, "Fixed merger contract wrote stdout: "+work);
		check(expected.equals(read(output)), "Fixed merger contract changed header/order/omissions: "+work);
		final List<String> lines=Files.readAllLines(output, StandardCharsets.US_ASCII);
		CovStatsLine.initializeHeader(lines.get(0));
		final CovStatsLine row=new CovStatsLine(lines.get(1));
		final boolean fieldsMatch=coreOnly ? coreMetricsMatch(row) && row.readGC==0 && row.median==0
				&& row.stdDev==0 && row.underMin==0 : metricsMatch(row) && row.underMin==600;
		write(work.resolve("reparse-observation.txt"), "fieldsMatch="+fieldsMatch+"\tunderMin="+row.underMin+"\n");
		check("otuA".equals(row.id) && fieldsMatch, "Fixed merger contract reparse differs: "+work);
	}

	/** These new guards are post-repair contracts, not observations captured on unmodified production. */
	private static void rejectHeader(final Path root, final Path launcher, final boolean duplicate) throws Exception{
		final Path work=Files.createDirectory(root.resolve(duplicate ? "fixed-duplicate-contract" : "fixed-unsupported-contract"));
		final String header=CORE+(duplicate ? "\tUnder_5/500\tUnder_min" : "\tMystery");
		final String row=MERGE_CORE_INPUT+(duplicate ? "\t600\t600" : "\t123");
		final Path input=work.resolve("input.tsv"), output=work.resolve("must-not-exist.tsv");
		write(input, header+"\n"+row+"\n");
		final int exit=run(launcher, Arrays.asList("in="+input, "out="+output, "-Xmx128m", "-Xms32m"), work);
		check(exit==1, "Invalid header must fail clearly, rather than time out; exit="+exit+"; "+work);
		check(!Files.exists(output), "Invalid header created output before rejection: "+work);
		check(Files.size(work.resolve("stdout"))==0, "Invalid header wrote stdout: "+work);
		final String diagnostics=read(work.resolve("stderr"));
		check(diagnostics.contains("IllegalArgumentException") && diagnostics.contains(duplicate
				? "Duplicate coverage column meaning: Under_min" : "Unsupported coverage column: Mystery"),
				"Invalid header lacks the expected diagnostic: "+work);
	}

	/** Checks independent literal input values shared by direct and compatible-order controls. */
	private static boolean metricsMatch(final CovStatsLine row){
		return coreMetricsMatch(row)
				&& row.readGC==0.25 && row.median==8 && row.stdDev==1.5;
	}
	private static boolean coreMetricsMatch(final CovStatsLine row){
		return row.avgFold==10 && row.length==1000 && row.refGC==0.4 && row.coveredBases==900
				&& row.coveredPercent()==90 && row.plusReads==7 && row.minusReads==3 && row.reads()==10;
	}

	private static boolean stat(final String diagnostics, final String label, final long value){
		return Pattern.compile("(?m)^"+Pattern.quote(label)+":\\s+"+value+"(?:\\s|$)").matcher(diagnostics).find();
	}
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII));
	}
	private static String read(final Path path) throws Exception{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	/** Supervises only this invocation's owned process group; saves the complete literal command and diagnostics. */
	private static int run(final Path launcher, final List<String> arguments, final Path work) throws Exception{
		return run(launcher, arguments, work, launcher.getParent());
	}
	/** Keeps the same supervision with an explicit cwd for fixture scripts located under dev/. */
	private static int run(final Path launcher, final List<String> arguments, final Path work,
			final Path workingDirectory) throws Exception{
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("/usr/bin/timeout", "--verbose",
				"--signal=TERM", "--kill-after=2s", "30s", "bash", launcher.toString()));
		command.addAll(arguments);
		Files.write(work.resolve("args.txt"), command, StandardCharsets.UTF_8);
		final Process process=new ProcessBuilder(command).directory(workingDirectory.toFile())
				.redirectOutput(work.resolve("stdout").toFile()).redirectError(work.resolve("stderr").toFile()).start();
		Throwable primary=null;
		try{
			process.getOutputStream().close();
			if(!process.waitFor(60, TimeUnit.SECONDS)){
				Files.write(work.resolve("timeout.txt"), Collections.singletonList("Owned launcher exceeded60 seconds."), StandardCharsets.US_ASCII);
				throw new AssertionError("Finite coverage baseline exceeded60 seconds: "+work);
			}
			final int exit=process.exitValue();
			write(work.resolve("exit.txt"), exit+"\n");
			return exit;
		}catch(Exception e){
			primary=e; throw e;
		}catch(Error e){
			primary=e; throw e;
		}finally{
			try{finishOwnedChild(process);}
			catch(RuntimeException e){if(primary==null){throw e;}else if(primary!=e){primary.addSuppressed(e);}}
			catch(Error e){if(primary==null){throw e;}else if(primary!=e){primary.addSuppressed(e);}}
		}
	}

	/** Allows GNU timeout's two-second group KILL grace before forcing/reaping its exact owned supervisor. */
	private static void finishOwnedChild(final Process process){
		if(!process.isAlive()){return;}
		process.destroy();
		if(waitOwnedChild(process, 5)){return;}
		process.destroyForcibly();
		if(!waitOwnedChild(process, 5)){throw new AssertionError("Owned coverage launcher did not terminate after forced cleanup.");}
	}
	/** Keeps repeated interrupts from removing the cleanup deadline; restores interruption on return. */
	private static boolean waitOwnedChild(final Process process, final int seconds){
		final long deadline=System.nanoTime()+TimeUnit.SECONDS.toNanos(seconds);
		boolean interrupted=false;
		try{
			while(process.isAlive()){
				final long remaining=deadline-System.nanoTime();
				if(remaining<=0){return false;}
				try{return process.waitFor(remaining, TimeUnit.NANOSECONDS);}
				catch(InterruptedException e){interrupted=true;}
			}
			return true;
		}finally{if(interrupted){Thread.currentThread().interrupt();}}
	}

	private static final String CORE="#ID\tAvg_fold\tLength\tRef_GC\tCovered_percent\tCovered_bases\tPlus_reads\tMinus_reads";
	private static final String HEADER_WINDOW=CORE+"\tRead_GC\tMedian_fold\tStd_Dev\tUnder_5/500";
	private static final String HEADER_MIN=CORE+"\tRead_GC\tMedian_fold\tStd_Dev\tUnder_min";
	private static final String HEADER_COMPATIBLE=CORE+"\tMedian_fold\tUnder_min\tRead_GC\tStd_Dev";
	private static final String ROW="ctg\t10\t1000\t0.4\t90\t900\t7\t3\t0.25\t8\t1.5\t600";
	private static final String MEDIAN_FIRST="ctgA otuA\t10\t1000\t0.4\t90\t900\t7\t3\t0.25\t8\t1.5\t600";
	private static final String MEDIAN_SECOND="ctgB otuA\t10\t1000\t0.4\t90\t900\t7\t3\t0.25\t8\t2.5\t600";
	private static final String MEDIAN_ZERO_FIRST="zeroA otuA\t3\t0\t0.4\t0\t0\t0\t0\t0.25\t8\t1.5\t0";
	private static final String MEDIAN_ZERO_SECOND="zeroB otuA\t9\t0\t0.6\t0\t0\t0\t0\t0.75\t8\t2.5\t0";
	private static final String MEDIAN_SINGLETON_OUTPUT="otuA\t10.0000\t1000\t0.4000\t90.0000\t900\t7\t3\t0.2500\t8\t1.5000\t600";
	private static final String MERGE_INPUT="contigA otuA\t10\t1000\t0.4\t90\t900\t7\t3\t0.25\t8\t1.5\t600";
	private static final String MERGE_CORE_INPUT="contigA otuA\t10\t1000\t0.4\t90\t900\t7\t3";
	private static final String MERGE_COMPATIBLE="contigA otuA\t10\t1000\t0.4\t90\t900\t7\t3\t8\t600\t0.25\t1.5";
	private static final String SERIALIZED="otuA\t10.0000\t1000\t0.4000\t90.0000\t900\t7\t3\t8\t600\t0.2500\t1.5000";
	private static final String HEADER_PRESERVING="otuA\t10.0000\t1000\t0.4000\t90.0000\t900\t7\t3\t0.2500\t8\t1.5000\t600";
	private static final String CORE_OUTPUT="otuA\t10.0000\t1000\t0.4000\t90.0000\t900\t7\t3";
	private static final String OVERFLOW_HEADER=CORE+"\tUnder_5/500";
	private static final String[] OVERFLOW_LARGE_ROWS={
			"a otu\t1\t1500000000\t0.4\t93.3333\t1400000000\t1\t0\t1200000000",
			"b otu\t1\t1500000000\t0.4\t93.3333\t1400000000\t1\t0\t1200000000"};
	private static final String[] OVERFLOW_CONTROL_ROWS={
			"a otu\t1\t1024\t0.4\t93.7500\t960\t1\t0\t768",
			"b otu\t1\t1024\t0.4\t93.7500\t960\t1\t0\t768"};
	private static final String OVERFLOW_BASELINE_OUTPUT="otu\t3000000000.0000\t-1294967296\t1200000000.0000\t-149496729600.0000\t-1494967296\t2\t0\t-1894967296";
	private static final String OVERFLOW_CONTROL_OUTPUT="otu\t1.0000\t2048\t0.4000\t93.7500\t1920\t2\t0\t1536";
	private static final String OVERFLOW_CORRECT_OUTPUT="otu\t1.0000\t3000000000\t0.4000\t93.3333\t2800000000\t2\t0\t2400000000";
	private static final String[] OUTPUT_CLOSE_ROWS={
			"ctgA otuA\t10\t20\t0.25\t50\t10\t2\t3\t0.4\t10\t2\t1",
			"ctgB otuA\t25\t20\t0.375\t75\t15\t4\t1\t0.6\t20\t4\t2"};
	private static final String OUTPUT_CLOSE_EXPECTED=HEADER_MIN+"\n"
			+"otuA\t17.5000\t40\t0.3125\t62.5000\t25\t6\t4\t0.5000\t30\t2.0000\t3\n";
	private static final String NONFINITE_ROW="ctg\t10\t1000\t0.4\t90\t900\t7\t3\t0.25\t8\t1.5\t0";
	private static final int[] SCALAR_COLUMNS={1, 3, 8, 10};
	private static final String[] SCALAR_CASE_NAMES={"avg-fold", "ref-gc", "read-gc", "std-dev"};
	private static final String[] SCALAR_HEADERS={"Avg_fold", "Ref_GC", "Read_GC", "Std_Dev"};
	private static final String INPUT_CLOSE_URL="jar:file:/controlled-inputclose.jar!/coverage.tsv";
	private static final String INPUT_CLOSE_FAILURE="Controlled owned input close failure";
	private static final String INPUT_CLOSE_REJECTION="driver.MergeCoverageOTU terminated in an error state while closing the input; output was not written.";
	private static int checks;
}
