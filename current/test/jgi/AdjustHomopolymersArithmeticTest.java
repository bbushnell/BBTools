package test.jgi;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

import jgi.AdjustHomopolymers;
import parse.Parse;
import structures.ByteBuilder;

/**
 * Scalar-only adjustment baseline and fixed arithmetic expectations.
 * Invokes the production helper without constructing a driver, Read or ByteBuilder.
 * Large lengths are integers only: no sequence, quality or scratch allocation,
 * processing loop, output writer or extreme-rate read editing is performed.
 * Baseline records two incorrect full-deletion endpoints and five unsafe accepts.
 * Fixed expectations delete defined runs exactly at -1 and reject excess capacity.
 * This checks arithmetic decisions, not allocation success or huge-read execution.
 * @author Sucrose
 */
public final class AdjustHomopolymersArithmeticTest{

	/** Requires phase=baseline|fixed and a new evidence directory. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown arithmetic fixture option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null, "Set a new outdir.");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		Files.createDirectory(root);
		final AdjustHomopolymersArithmeticTest test=new AdjustHomopolymersArithmeticTest(root, phase.equals("baseline"));
		for(final Row row : ROWS){test.run(row);}
		test.check(test.cases==14, "Frozen scalar matrix requires 14 cases: "+test.cases);
		test.check(test.helperCalls==14, "Call the actual production helper once per case: "+test.helperCalls);
		test.check(test.prefixMismatches==(test.baseline ? 2 : 0), "Exact -1 deletion must preserve the prefix; mismatches="+test.prefixMismatches);
		test.check(test.unsafeAccepted==(test.baseline ? 5 : 0), "No fixed adjustment may erase an earlier run or exceed MAXLEN; accepts="+test.unsafeAccepted);
		test.check(test.rejected==(test.baseline ? 0 : 4), "Reject the four fixed expansion overflows; rejections="+test.rejected);
		if(test.failures>0){throw new AssertionError("Arithmetic fixture failed "+test.failures+" of "+test.checks+" checks; inspect "+root);}
		System.out.println("ADJUST_HOMOPOLYMERS_ARITHMETIC_PASS phase="+phase+" cases="+test.cases+
			" helperCalls="+test.helperCalls+" prefixMismatches="+test.prefixMismatches+
			" unsafeAccepted="+test.unsafeAccepted+" rejected="+test.rejected+" checks="+test.checks);
	}

	/** Binds the static production seam without creating a driver or changing globals. */
	private AdjustHomopolymersArithmeticTest(final Path root_, final boolean baseline_) throws Exception{
		root=root_; baseline=baseline_;
		helper=AdjustHomopolymers.class.getDeclaredMethod("adjustmentForRun", float.class, int.class, int.class);
		helper.setAccessible(true);
		Files.write(root.resolve("cases.tsv"), ("case\trate\tstreak\tcurrentLength\tbaselineAdjustment\tbaselineProjected"+
			"\tfixedExpectation\tactualAdjustment\tactualProjected\tactualException\tactualMessage\n").getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE);
	}

	/** Records the actual helper outcome and compares independent literal expectations. */
	private void run(final Row row) throws Exception{
		require(row.streak>0 && row.length>=row.streak && row.length<=MAXLEN &&
			Float.isFinite(row.rate) && row.rate!=0 && row.rate>=-1,
			"Fixture scalar state must satisfy the production caller's rate/run contract: "+row.name);
		cases++;
		Integer adjustment=null;
		Throwable error=null;
		helperCalls++;
		try{adjustment=(Integer)helper.invoke(null, row.rate, row.streak, row.length);}
		catch(InvocationTargetException e){error=e.getCause();}
		final Long projected=(adjustment==null ? null : (long)row.length+adjustment);
		if(error!=null){rejected++;}
		if(adjustment!=null){
			if(adjustment<-row.streak || projected>MAXLEN){unsafeAccepted++;}
			if(row.rate==-1 && projected!=row.fixedProjected){prefixMismatches++;}
		}
		final boolean expectRejection=!baseline && row.rejectFixed;
		final int expectedAdjustment=baseline ? row.baselineAdjustment : row.fixedAdjustment;
		final long expectedProjected=baseline ? row.baselineProjected : row.fixedProjected;
		final String message=(error==null ? "" : String.valueOf(error.getMessage()));
		final String observed=row.name+"\t"+row.rate+"\t"+row.streak+"\t"+row.length+
			"\t"+row.baselineAdjustment+"\t"+row.baselineProjected+"\t"+
			(row.rejectFixed ? "IllegalArgumentException" : row.fixedAdjustment+":"+row.fixedProjected)+
			"\t"+adjustment+"\t"+projected+"\t"+(error==null ? "none" : error.getClass().getName())+
			"\t"+message.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ')+"\n";
		Files.write(root.resolve("cases.tsv"), observed.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
		check(expectRejection ? error!=null && error.getClass()==IllegalArgumentException.class : error==null && adjustment!=null,
			row.name+": expected "+(expectRejection ? "IllegalArgumentException" : "scalar return")+", saw "+error);
		check(expectRejection ? message.contains(Float.toString(row.rate)) &&
			message.contains(Integer.toString(row.streak)) && message.contains(Integer.toString(row.length)) :
			adjustment!=null && adjustment==expectedAdjustment && projected==expectedProjected,
			row.name+": expected "+(expectRejection ? "rate/streak/currentLength in diagnostic" : expectedAdjustment+":"+expectedProjected)+
			", saw "+adjustment+":"+projected+" "+message);
	}

	/** Accumulates failures so the complete scalar matrix remains available for review. */
	private void check(final boolean correct, final String message){
		checks++;
		if(!correct){failures++; System.err.println("FAIL: "+message);}
	}

	/** Rejects malformed fixture setup before interpreting results. */
	private static void require(final boolean correct, final String message){
		if(!correct){throw new IllegalArgumentException(message);}
	}

	/** One literal numeric oracle; rejected fixed rows have no return-value expectation. */
	private static final class Row{
		Row(final String name_, final float rate_, final int streak_, final int length_,
				final int baselineAdjustment_, final long baselineProjected_, final int fixedAdjustment_,
				final long fixedProjected_, final boolean rejectFixed_){
			name=name_; rate=rate_; streak=streak_; length=length_;
			baselineAdjustment=baselineAdjustment_; baselineProjected=baselineProjected_;
			fixedAdjustment=fixedAdjustment_; fixedProjected=fixedProjected_; rejectFixed=rejectFixed_;
		}
		final String name;
		final float rate;
		final int streak, length, baselineAdjustment, fixedAdjustment;
		final long baselineProjected, fixedProjected;
		final boolean rejectFixed;
	}

	private final Path root;
	private final boolean baseline;
	private final Method helper;
	private int cases, helperCalls, prefixMismatches, unsafeAccepted, rejected, checks, failures;
	/** Compile-time capacity constant; no ByteBuilder is instantiated. */
	private static final int MAXLEN=ByteBuilder.MAXLEN;
	private static final Row[] ROWS={
		new Row("plus-two", 2f, 3, 13, 6, 19, 6, 19, false),
		new Row("minus-half", -0.5f, 3, 13, -1, 12, -1, 12, false),
		new Row("seven-tenths", 0.7f, 10, 20, 7, 27, 7, 27, false),
		new Row("quarter-truncates", 0.25f, 3, 13, 0, 13, 0, 13, false),
		new Row("minus-one-small", -1f, 3, 13, -3, 10, -3, 10, false),
		new Row("minus-one-2p24", -1f, 16777216, 16777226, -16777216, 10, -16777216, 10, false),
		new Row("minus-one-2p24-plus-one", -1f, 16777217, 16777227, -16777216, 11, -16777217, 10, false),
		new Row("minus-one-2p24-plus-three", -1f, 16777219, 16777229, -16777220, 9, -16777219, 10, false),
		new Row("capacity-exact", 1f, 1, MAXLEN-1, 1, 2147483627L, 1, 2147483627L, false),
		new Row("capacity-zero", 0.25f, 1, MAXLEN, 0, 2147483627L, 0, 2147483627L, false),
		new Row("capacity-plus-one", 1f, 1, MAXLEN, 1, 2147483628L, 0, 0, true),
		new Row("capacity-int-overflow", 32f, 1, MAXLEN, 32, 2147483659L, 0, 0, true),
		new Row("saturate-finite", Float.MAX_VALUE, 1, 11, Integer.MAX_VALUE, 2147483658L, 0, 0, true),
		new Row("saturate-infinite", Float.MAX_VALUE, 2, 12, Integer.MAX_VALUE, 2147483659L, 0, 0, true)
	};
}
