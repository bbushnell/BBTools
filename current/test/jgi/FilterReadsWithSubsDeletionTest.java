package test.jgi;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

import jgi.FilterReadsWithSubs;
import parse.LineParser1;
import parse.Parse;
import stream.Read;
import stream.SamLine;

/**
 * Four finite SAM cases for terminal-deletion bounds capture in
 * FilterReadsWithSubs.  Each process handles one literal case so the SAM
 * parser, private filter object and mutable stream globals start fresh.
 * Baseline mode requires observed old behavior; fixed mode requires intended
 * corrected behavior.
 * @author Jean
 */
public final class FilterReadsWithSubsDeletionTest{

	/** Requires baseline or fixed phase, a named case and a fresh output directory. */
	public static void main(final String[] args) throws Exception{
		String phase=null, caseName=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("case")){caseName=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown deletion fixture option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase),
			"Set phase=baseline or phase=fixed.");
		require(caseName!=null && !caseName.isEmpty(), "Set case=<literal>.");
		require(directory!=null && !directory.isEmpty(), "Set a fresh outdir.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Run through dev/filtersubsdeletiontest.sh with assertions enabled.");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		Files.createDirectory(root);
		new FilterReadsWithSubsDeletionTest(root, phase).row(CaseLiteral.byName(caseName));
	}

	/** Resolves the private processRead method once. */
	private FilterReadsWithSubsDeletionTest(final Path root_, final String phase_)
			throws Exception{
		root=root_;
		phase=phase_;
		processRead=FilterReadsWithSubs.class.getDeclaredMethod("processRead", Read.class);
		processRead.setAccessible(true);
		Files.write(root.resolve("filtersubs-deletion-cases.tsv"),
			(HEADER+"\n").getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Runs one literal SAM row through normal SamLine parsing and private filtering. */
	private void row(final CaseLiteral literal) throws Exception{
		Throwable failure=null;
		String failureOrigin="none";
		String failureTopFrame="none";
		Boolean actualKeep=null;
		String shortMatch="unavailable";
		String expandedMatch="unavailable";
		String matchCounts="unavailable";
		int basesLength=-1, qualsLength=-1, matchLength=-1, expandedLength=-1;
		try{
			final FilterReadsWithSubs filter=new FilterReadsWithSubs(new String[]{
				"in=stdin", "out=stdout.sam", "keepperfect=t", "countindels=t",
				"minsubs=1", "minq=0", "maxq=99"});
			final SamLine sam=new SamLine(LINE_PARSER.set(literal.samLine.
				getBytes(StandardCharsets.US_ASCII)));
			final Read read=sam.toRead(false);
			if(!read.validated()){read.validate(true);}
			basesLength=read.bases==null ? -1 : read.bases.length;
			qualsLength=read.quality==null ? -1 : read.quality.length;
			if(read.match!=null){
				matchLength=read.match.length;
				shortMatch=new String(read.match, StandardCharsets.US_ASCII);
				final byte[] expanded=read.shortmatch() ?
					Read.toLongMatchString(read.match) : read.match;
				expandedLength=expanded.length;
				expandedMatch=new String(expanded, StandardCharsets.US_ASCII);
				matchCounts=counts(expanded);
			}
			actualKeep=(Boolean)processRead.invoke(filter, read);
		}catch(final InvocationTargetException ite){
			failure=ite.getCause();
			failureOrigin="processRead";
		}catch(final Throwable caught){
			failure=caught;
			failureOrigin="setup";
		}
		final String actualException=type(failure);
		failureTopFrame=topFrame(failure);
		final boolean baseline=literal.oldException.equals(actualException) &&
			literal.oldOrigin.equals(failureOrigin) &&
			literal.oldTopFrame.equals(failureTopFrame) &&
			(literal.oldKeep==null ? actualKeep==null : literal.oldKeep.equals(actualKeep));
		final boolean correct=literal.intendedException.equals(actualException) &&
			literal.intendedOrigin.equals(failureOrigin) &&
			literal.intendedTopFrame.equals(failureTopFrame) &&
			(literal.intendedKeep==null ? actualKeep==null :
				literal.intendedKeep.equals(actualKeep));
		Files.write(root.resolve(literal.name+".sam"),
			(literal.samLine+"\n").getBytes(StandardCharsets.US_ASCII),
			StandardOpenOption.CREATE_NEW);
		final String row=literal.name+'\t'+phase+'\t'+literal.cigar+'\t'+literal.md+
			"\tACG\tIII\t"+literal.oldKeep+'\t'+literal.oldException+
			'\t'+literal.oldOrigin+'\t'+literal.oldTopFrame+
			'\t'+literal.intendedKeep+'\t'+literal.intendedException+
			'\t'+literal.intendedOrigin+'\t'+literal.intendedTopFrame+
			'\t'+shortMatch+'\t'+expandedMatch+'\t'+matchCounts+
			'\t'+basesLength+'\t'+qualsLength+'\t'+matchLength+
			'\t'+expandedLength+'\t'+actualKeep+'\t'+actualException+
			'\t'+failureOrigin+'\t'+failureTopFrame+'\t'+clean(failure==null ?
				"" : String.valueOf(failure.getMessage()))+'\t'+correct+'\t'+baseline+'\n';
		Files.write(root.resolve("filtersubs-deletion-cases.tsv"),
			row.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);
		final boolean expected="baseline".equals(phase) ? baseline : correct;
		final String label="baseline".equals(phase) ? "BASELINE_CAPTURE" : "FIXED_CAPTURE";
		final String summary="FILTERSUBS_DELETION_"+label+'_'+
			(expected ? "PASS" : "FAIL")+" phase="+phase+" case="+literal.name+
			" baselinePredictionMatched="+baseline+" intendedResultMatched="+correct+
			" exceptionOrigin="+failureOrigin+" topFrame="+failureTopFrame+
			" actualException="+actualException+"\n";
		Files.write(root.resolve("summary.txt"), summary.getBytes(StandardCharsets.UTF_8),
			StandardOpenOption.CREATE_NEW);
		System.out.print(summary);
		if(!expected){throw new AssertionError(summary+
			"Inspect saved filtersubs deletion row in "+root);}
	}

	/** Counts the match alphabet used by FilterReadsWithSubs plus a final other bin. */
	private static String counts(final byte[] match){
		int m=0, s=0, c=0, n=0, iCount=0, d=0, other=0;
		for(final byte symbol : match){
			if(symbol=='m'){m++;}
			else if(symbol=='S'){s++;}
			else if(symbol=='C'){c++;}
			else if(symbol=='N'){n++;}
			else if(symbol=='I'){iCount++;}
			else if(symbol=='D'){d++;}
			else{other++;}
		}
		return "m="+m+",S="+s+",C="+c+",N="+n+",I="+iCount+
			",D="+d+",other="+other;
	}

	/** Makes one exception message fit in a TSV cell. */
	private static String clean(final String value){
		return value.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');
	}

	/** Records absent exceptions explicitly. */
	private static String type(final Throwable error){
		return error==null ? "none" : error.getClass().getName();
	}

	/** Returns the top application frame without pinning a line number or JDK message. */
	private static String topFrame(final Throwable error){
		if(error==null){return "none";}
		final StackTraceElement[] trace=error.getStackTrace();
		if(trace.length<1){return "no-stack";}
		final StackTraceElement top=trace[0];
		return top.getClassName()+"."+top.getMethodName();
	}

	/** Rejects fixture misuse. */
	private static void require(final boolean correct, final String message){
		if(!correct){throw new IllegalArgumentException(message);}
	}

	/** One literal SAM record plus its old and intended one-row oracles. */
	private static final class CaseLiteral{
		static CaseLiteral byName(final String name){
			for(final CaseLiteral literal : CASES){
				if(literal.name.equals(name)){return literal;}
			}
			throw new IllegalArgumentException("Unknown case "+name);
		}

		CaseLiteral(final String name_, final String cigar_, final String md_,
				final Boolean oldKeep_, final String oldException_,
				final String oldOrigin_, final String oldTopFrame_,
				final Boolean intendedKeep_, final String intendedException_,
				final String intendedOrigin_, final String intendedTopFrame_){
			name=name_;
			cigar=cigar_;
			md=md_;
			oldKeep=oldKeep_;
			oldException=oldException_;
			oldOrigin=oldOrigin_;
			oldTopFrame=oldTopFrame_;
			intendedKeep=intendedKeep_;
			intendedException=intendedException_;
			intendedOrigin=intendedOrigin_;
			intendedTopFrame=intendedTopFrame_;
			samLine=name+"\t0\tref\t1\t60\t"+cigar+
				"\t*\t0\t0\tACG\tIII\tMD:Z:"+md;
		}

		final String name;
		final String cigar;
		final String md;
		final Boolean oldKeep;
		final String oldException;
		final String oldOrigin;
		final String oldTopFrame;
		final Boolean intendedKeep;
		final String intendedException;
		final String intendedOrigin;
		final String intendedTopFrame;
		final String samLine;
	}

	private final Path root;
	private final String phase;
	private final Method processRead;

	private static final LineParser1 LINE_PARSER=new LineParser1('\t');
	private static final String AIOOBE="java.lang.ArrayIndexOutOfBoundsException";
	private static final String PROCESS_READ="jgi.FilterReadsWithSubs.processRead";
	private static final String HEADER="case\tphase\tcigar\tmd\tseq\tqual"+
		"\tpredictedOldKeep\tpredictedOldException\tpredictedOldExceptionOrigin"+
		"\tpredictedOldTopFrame"+
		"\tintendedKeep\tintendedException\tintendedExceptionOrigin"+
		"\tintendedTopFrame\tshortMatch\texpandedMatch"+
		"\tmatchCounts\tbasesLength\tqualsLength\tshortMatchLength"+
		"\texpandedMatchLength\tactualKeep\texceptionClass\texceptionOrigin"+
		"\ttopFrame\texceptionMessage\tintendedResultMatched\tbaselinePredictionMatched";

	private static final CaseLiteral[] CASES={
		new CaseLiteral("perfect-3eq", "3=", "3", Boolean.TRUE, "none",
			"none", "none", Boolean.TRUE, "none", "none", "none"),
		new CaseLiteral("internal-1eq-1d-2eq", "1=1D2=", "1^A2",
			Boolean.TRUE, "none", "none", "none", Boolean.TRUE, "none",
			"none", "none"),
		new CaseLiteral("terminal-3eq-1d", "3=1D", "3^A0",
			null, AIOOBE, "processRead", PROCESS_READ, Boolean.TRUE, "none",
			"none", "none"),
		new CaseLiteral("terminal-3m-1d", "3M1D", "3^A0",
			null, AIOOBE, "processRead", PROCESS_READ, Boolean.TRUE, "none",
			"none", "none")
	};
}
