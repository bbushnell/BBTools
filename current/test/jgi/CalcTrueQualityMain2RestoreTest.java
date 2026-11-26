package test.jgi;

import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ReadWrite;
import jgi.CalcTrueQuality;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.FASTQ;
import stream.SamLine;
import tracker.ReadStats;

/**
 * Same-JVM fixture for CalcTrueQuality.main2 restoration after main throws.
 * The controlled missing-input failure exercises early static mutations without
 * starting readers, workers, variant callers, or matrix writers.
 * @author Jean
 */
public final class CalcTrueQualityMain2RestoreTest{

	/** Runs two independent saved-state profiles in the requested phase. */
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CalcTrueQualityMain2RestoreTest.class, false);
		final SavedState initial=new SavedState();
		final boolean initialParse8=SamLine.PARSE_8;
		final boolean initialParseOptionalMdOnly=SamLine.PARSE_OPTIONAL_MD_ONLY;
		try{
			final Args parsed=new Args(pp.args);
			new CalcTrueQualityMain2RestoreTest(parsed.fixed).process();
		}finally{
			initial.restoreAll();
			SamLine.PARSE_8=initialParse8;
			SamLine.PARSE_OPTIONAL_MD_ONLY=initialParseOptionalMdOnly;
			Shared.closeStream(pp.outstream);
		}
	}

	/** Stores validated arguments. */
	private CalcTrueQualityMain2RestoreTest(final boolean fixed_){
		fixed=fixed_;
	}

	/** Runs distinct profiles to prove restoration is to the immediate saved state. */
	private void process(){
		runProfile("profile1", SavedState.profile1(), 44);
		runProfile("profile2", SavedState.profile2(), 66);
		System.out.println("CALCTRUEQUALITY_MAIN2_RESTORE_PASS fixed="+fixed);
	}

	/** Applies one saved-state profile, calls main2, and checks baseline/fixed state. */
	private void runProfile(final String name, final SavedState saved, final int dirtyBuffer){
		saved.restoreAll();
		saved.checkCurrent(name+" pre-call saved");
		saved.checkWitnesses(name, dirtyBuffer);

		final String[] args={"readbufferlen="+dirtyBuffer, "bf1=t"};
		SavedState observed=null;
		try{
			CalcTrueQuality.main2(args);
			throw fail(name+" main2 unexpectedly returned");
		}catch(RuntimeException e){
			check(e.getClass()==RuntimeException.class,
					name+" unexpected main2 exception class: "+e.getClass().getName());
			check(EXPECTED_MISSING_INPUT.equals(e.getMessage()),
					name+" unexpected main2 exception: "+e);
			observed=new SavedState();
		}

		final SavedState expected=(fixed ? saved : saved.leaked(dirtyBuffer));
		observed.checkMatches(expected, name+" "+(fixed ? "fixed" : "baseline"));
		check(!observed.bf3, name+" BF3 changed");
		check(!observed.bf4, name+" BF4 changed");
		System.out.println("MAIN2_RESTORE_PROFILE_PASS "+name+" fixed="+fixed);
	}

	/** Snapshot of the saved main2 statics plus BF3/BF4 fixture-hygiene guards. */
	private static final class SavedState{
		SavedState(){
			collectQualityStats=ReadStats.COLLECT_QUALITY_STATS;
			overwrite=ReadStats.overwrite;
			bufferLen=Shared.bufferLen();
			usePigz=ReadWrite.USE_PIGZ;
			useUnpigz=ReadWrite.USE_UNPIGZ;
			zipLevel=ReadWrite.ZIPLEVEL;
			bf1=ByteFile.FORCE_MODE_BF1;
			bf2=ByteFile.FORCE_MODE_BF2;
			bf3=ByteFile.FORCE_MODE_BF3;
			bf4=ByteFile.FORCE_MODE_BF4;
			testInterleaved=FASTQ.TEST_INTERLEAVED;
			forceInterleaved=FASTQ.FORCE_INTERLEAVED;
		}

		SavedState(final boolean collectQualityStats_, final boolean overwrite_,
				final int bufferLen_, final boolean usePigz_, final boolean useUnpigz_,
				final int zipLevel_, final boolean bf1_, final boolean bf2_,
				final boolean testInterleaved_, final boolean forceInterleaved_){
			collectQualityStats=collectQualityStats_;
			overwrite=overwrite_;
			bufferLen=bufferLen_;
			usePigz=usePigz_;
			useUnpigz=useUnpigz_;
			zipLevel=zipLevel_;
			bf1=bf1_;
			bf2=bf2_;
			bf3=false;
			bf4=false;
			testInterleaved=testInterleaved_;
			forceInterleaved=forceInterleaved_;
		}

		/** First pre-call profile; every mutable early field differs from CTQ's writes. */
		static SavedState profile1(){
			return new SavedState(false, false, 33, true, false, 3,
					false, true, true, true);
		}

		/** Second pre-call profile, using distinct unchanged control and numeric sentinels. */
		static SavedState profile2(){
			return new SavedState(false, true, 55, true, false, 5,
					false, true, true, true);
		}

		/** Baseline state after main2 throws before its current restoration block. */
		SavedState leaked(final int dirtyBuffer){
			return new SavedState(true, overwrite, dirtyBuffer, false, true, 8,
					true, false, false, false);
		}

		/** Checks that this profile witnesses every early main2 mutation. */
		void checkWitnesses(final String label, final int dirtyBuffer){
			check(!collectQualityStats, label+" collectQualityStats precondition");
			check(bufferLen!=dirtyBuffer, label+" bufferLen precondition");
			check(usePigz, label+" USE_PIGZ precondition");
			check(!useUnpigz, label+" USE_UNPIGZ precondition");
			check(zipLevel!=8, label+" ZIPLEVEL precondition");
			check(!bf1, label+" BF1 precondition");
			check(bf2, label+" BF2 precondition");
			check(!bf3, label+" BF3 guard precondition");
			check(!bf4, label+" BF4 guard precondition");
			check(testInterleaved, label+" TEST_INTERLEAVED precondition");
			check(forceInterleaved, label+" FORCE_INTERLEAVED precondition");
		}

		/** Applies this state to the actual process statics. */
		void restoreAll(){
			ReadStats.COLLECT_QUALITY_STATS=collectQualityStats;
			ReadStats.overwrite=overwrite;
			Shared.setBufferLen(bufferLen);
			ReadWrite.USE_PIGZ=usePigz;
			ReadWrite.USE_UNPIGZ=useUnpigz;
			ReadWrite.ZIPLEVEL=zipLevel;
			ByteFile.FORCE_MODE_BF1=bf1;
			ByteFile.FORCE_MODE_BF2=bf2;
			ByteFile.FORCE_MODE_BF3=bf3;
			ByteFile.FORCE_MODE_BF4=bf4;
			FASTQ.TEST_INTERLEAVED=testInterleaved;
			FASTQ.FORCE_INTERLEAVED=forceInterleaved;
		}

		/** Checks the current process statics against this expected state. */
		void checkCurrent(final String label){
			new SavedState().checkMatches(this, label);
		}

		/** Checks this snapshot against the expected state. */
		void checkMatches(final SavedState expected, final String label){
			check(collectQualityStats==expected.collectQualityStats,
					label+" collectQualityStats");
			check(overwrite==expected.overwrite, label+" ReadStats.overwrite");
			check(bufferLen==expected.bufferLen, label+" bufferLen expected="+
					expected.bufferLen+" observed="+bufferLen);
			check(usePigz==expected.usePigz, label+" USE_PIGZ");
			check(useUnpigz==expected.useUnpigz, label+" USE_UNPIGZ");
			check(zipLevel==expected.zipLevel, label+" ZIPLEVEL expected="+
					expected.zipLevel+" observed="+zipLevel);
			check(bf1==expected.bf1, label+" BF1");
			check(bf2==expected.bf2, label+" BF2");
			check(bf3==expected.bf3, label+" BF3");
			check(bf4==expected.bf4, label+" BF4");
			check(testInterleaved==expected.testInterleaved, label+" TEST_INTERLEAVED");
			check(forceInterleaved==expected.forceInterleaved, label+" FORCE_INTERLEAVED");
		}

		boolean collectQualityStats;
		boolean overwrite;
		int bufferLen;
		boolean usePigz;
		boolean useUnpigz;
		int zipLevel;
		boolean bf1;
		boolean bf2;
		boolean bf3;
		boolean bf4;
		boolean testInterleaved;
		boolean forceInterleaved;
	}

	/** Parses and stores CalcTrueQualityMain2RestoreTest arguments. */
	private static final class Args{
		Args(final String[] args0){
			final String[] args=(args0==null ? new String[0] : args0);
			final Parser parser=new Parser();
			for(final String arg : args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String key=split[0].toLowerCase(Locale.ROOT);
				final String value=(split.length>1 ? split[1] : null);
				if("fixed".equals(key)){
					fixed=Parse.parseBoolean(value);
				}else if("phase".equals(key)){
					fixed="fixed".equalsIgnoreCase(value);
					check(fixed || "baseline".equalsIgnoreCase(value),
							"phase must be baseline or fixed: "+value);
				}else if(parser.parse(arg, key, value)){
					// Consumed by BBTools common parsing.
				}else{
					throw fail("Unknown argument: "+key);
				}
			}
		}

		/** Whether main2 should restore state in the product code. */
		boolean fixed;
	}

	/** Throws when a condition is false. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw fail(message);}
	}

	/** Returns a RuntimeException with a diagnostic message. */
	private static RuntimeException fail(final String message){
		return new RuntimeException(message);
	}

	/** Whether to require main2 restoration after the expected throw. */
	private final boolean fixed;

	/** Exact constructor failure used before pipeline startup. */
	private static final String EXPECTED_MISSING_INPUT="Error - at least one input file is required.";
}
