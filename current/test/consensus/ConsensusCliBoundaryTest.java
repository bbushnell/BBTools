package test.consensus;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.concurrent.TimeUnit;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Finite CLI regression oracles for consensus package boundary cases.
 * Each case runs the production main in a short child JVM to keep static parser
 * and consensus policy state from leaking between tools.
 *
 * @author Brian Bushnell, Jean
 */
public final class ConsensusCliBoundaryTest{

	/** Parses outdir and runs all cases, reporting every mismatched literal record. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, ConsensusCliBoundaryTest.class, false);
		try{
			String outdir=null;
			Phase phase=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("outdir")){outdir=value;}
				else if(key.equals("phase")){phase=Phase.parse(value);}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			boolean assertions=false;
			assert(assertions=true) : "Consensus CLI regression fixture requires -ea to retain fixture invariants";
			check(assertions, "Run ConsensusCliBoundaryTest with -ea.");
			check(outdir!=null, "Set outdir to a fresh directory whose parent exists.");
			check(phase!=null, "Set phase=baseline or phase=fixed.");
			new ConsensusCliBoundaryTest(Files.createDirectory(Paths.get(outdir)).toAbsolutePath(), phase).run();
		}finally{Shared.closeStream(pp.outstream);}
	}

	private ConsensusCliBoundaryTest(final Path root_, final Phase phase_){
		root=root_;
		phase=phase_;
	}

	/** Runs every independent oracle and fails only after observations are saved. */
	private void run() throws Exception{
		final List<String> failures=new ArrayList<String>();
		final List<CaseResult> results=new ArrayList<CaseResult>();
		for(final CaseSpec spec : CASES){
			try{
				results.add(spec.run(this));
			}catch(final Throwable failure){
				results.add(new CaseResult(spec.name, spec.expected,
						"THREW "+failure.getClass().getName()+": "+failure.getMessage(), "-", "-"));
				failures.add(spec.name+": "+failure);
			}
		}
		writeObservations(results);
		if(!failures.isEmpty()){
			throw new AssertionError(failures.size()+" consensus CLI case(s) failed: "+failures);
		}
		System.out.println("CONSENSUS_CLI_BOUNDARY_PASS phase="+phase.argument+" cases="+results.size());
	}

	private CaseResult noNsUnchanged(final CaseSpec spec) throws Exception{
		final Path dir=caseDir(spec.name);
		writeFasta(dir.resolve("ref.fa"), "ref", "ACGTACGT");
		writeSam(dir.resolve("mapped.sam"), 8, sam("long0", 99, 1, 60, "4M", 5, 8, "ACGT"));
		final Child child=runMain(dir, "consensus.FixScaffoldGaps",
				"in="+dir.resolve("mapped.sam"), "ref="+dir.resolve("ref.fa"),
				"out="+dir.resolve("fixed.fa"), "ow=t", "append=f", "t=1");
		final String fasta=read(dir.resolve("fixed.fa"));
		final String expected=">ref\nACGTACGT\n";
		check(child.exit==0, spec.name+" expected exit0: "+child);
		check(expected.equals(fasta), spec.name+" should keep an N-free reference unchanged; expected "
				+literal(expected)+" observed "+literal(fasta));
		return new CaseResult(spec.name, spec.expected, fasta, child.out, child.err);
	}

	private CaseResult shortInternalNs(final CaseSpec spec) throws Exception{
		final Path dir=caseDir(spec.name);
		writeFasta(dir.resolve("ref.fa"), "ref", "AAANNNCCC");
		writeSam(dir.resolve("mapped.sam"), 9, sam("short0", 99, 1, 60, "3M", 4, 4, "AAA"));
		return runFixGapsExpecting(spec, dir, "AAACCC", "AAANNNCCC");
	}

	private CaseResult edgeProximalNs(final CaseSpec spec) throws Exception{
		final Path dir=caseDir(spec.name);
		final String reference=repeat('A', 301)+repeat('N', 10)+"C";
		writeFasta(dir.resolve("ref.fa"), "ref", reference);
		writeSam(dir.resolve("mapped.sam"), reference.length(), sam("edge0", 99, 1, 60, "10M", 11, 20, repeat('A', 10)));
		return runFixGapsExpecting(spec, dir, repeat('A', 301)+"C", reference);
	}

	private CaseResult terminalNs(final CaseSpec spec) throws Exception{
		final Path dir=caseDir(spec.name);
		writeFasta(dir.resolve("ref.fa"), "ref", "AAA"+repeat('N', 10));
		writeSam(dir.resolve("mapped.sam"), 13, sam("terminal0", 99, 1, 60, "3M", 4, 3, "AAA"));
		return runFixGapsExpecting(spec, dir, "AAA", "AAA"+repeat('N', 10));
	}

	private CaseResult minDepthUnsupportedGap(final CaseSpec spec) throws Exception{
		final Path dir=caseDir(spec.name);
		final String reference=repeat('A', 350)+repeat('N', 10)+repeat('C', 350);
		writeFasta(dir.resolve("ref.fa"), "ref", reference);
		writeSameElevenRecordGapSam(dir.resolve("mapped.sam"), reference.length(), 100, 20);
		final CaseResult result=runFixGapsExpecting(spec, dir, reference, reference, "mindepth=10");
		if(phase==Phase.BASELINE){
			check(result.stderr.contains("Gaps Narrowed              1"),
					spec.name+" baseline should expose a below-mindepth resize attempt: "+literal(result.stderr));
			check(result.stderr.contains("Ns Removed                20"),
					spec.name+" baseline should count a hidden below-mindepth resize: "+literal(result.stderr));
		}else{
			check(result.stderr.contains("Gaps Narrowed              0"),
					spec.name+" fixed mode should gate the below-mindepth resize: "+literal(result.stderr));
			check(result.stderr.contains("Ns Removed                 0"),
					spec.name+" fixed mode should leave low-support gaps unmodified: "+literal(result.stderr));
		}
		return result;
	}

	private CaseResult zeroAcceptedPairs(final CaseSpec spec) throws Exception{
		final Path dir=caseDir(spec.name);
		final String reference=repeat('A', 350)+repeat('N', 10)+repeat('C', 350);
		writeFasta(dir.resolve("ref.fa"), "ref", reference);
		writeSam(dir.resolve("mapped.sam"), reference.length());
		final CaseResult result=runFixGapsExpecting(spec, dir, reference, reference);
		if(phase==Phase.BASELINE){
			check(result.stderr.contains("NaN"),
					spec.name+" baseline should expose 0/0 Average Insert: "+literal(result.stderr));
		}else{
			check(!result.stderr.contains("NaN") && !result.stderr.contains("Infinity"),
					spec.name+" must keep the empty accepted-pair summary finite: "+literal(result.stderr));
		}
		return result;
	}

	private CaseResult minDepthSupportedGap(final CaseSpec spec) throws Exception{
		final Path dir=caseDir(spec.name);
		final String reference=repeat('A', 350)+repeat('N', 10)+repeat('C', 350);
		final String expanded=repeat('A', 350)+repeat('N', 90)+repeat('C', 350);
		writeFasta(dir.resolve("ref.fa"), "ref", reference);
		writeSameElevenRecordGapSam(dir.resolve("mapped.sam"), reference.length(), 100, 20);
		final CaseResult result=runFixGapsExpecting(spec, dir, reference, expanded, "mindepth=1");
		if(phase==Phase.BASELINE){
			check(result.stderr.contains("Gaps Narrowed              1"),
					spec.name+" baseline should expose the sparse-histogram shrink attempt: "+literal(result.stderr));
			check(result.stderr.contains("Ns Removed                20"),
					spec.name+" baseline should count the hidden shrink attempt: "+literal(result.stderr));
		}else{
			check(result.stderr.contains("Gaps Widened               1"),
					spec.name+" fixed mode should widen a supported short gap: "+literal(result.stderr));
			check(result.stderr.contains("Ns Added                  80"),
					spec.name+" fixed mode should report the emitted widening: "+literal(result.stderr));
		}
		return result;
	}

	private CaseResult minimumClampStats(final CaseSpec spec) throws Exception{
		final Path dir=caseDir(spec.name);
		final String reference=repeat('A', 350)+repeat('N', 10)+repeat('C', 350);
		writeFasta(dir.resolve("ref.fa"), "ref", reference);
		writeSameElevenRecordGapSam(dir.resolve("mapped.sam"), reference.length(), 20, 100);
		final CaseResult result=runFixGapsExpecting(spec, dir, reference, reference, "mindepth=1");
		if(phase==Phase.FIXED){
			check(result.stderr.contains("Gaps Unchanged             1"),
					spec.name+" fixed mode should report the clamped emitted difference: "+literal(result.stderr));
			check(result.stderr.contains("Ns Removed                 0"),
					spec.name+" fixed mode should not count raw shrinkage hidden by the minimum-gap clamp: "
					+literal(result.stderr));
		}
		return result;
	}

	private CaseResult zeroMinDepthRejected(final CaseSpec spec) throws Exception{
		final Path dir=caseDir(spec.name);
		final String reference=repeat('A', 350)+repeat('N', 10)+repeat('C', 350);
		writeFasta(dir.resolve("ref.fa"), "ref", reference);
		writeDeepInsertSam(dir.resolve("mapped.sam"), reference.length(), 10);
		final Child child=runMain(dir, "consensus.FixScaffoldGaps",
				"in="+dir.resolve("mapped.sam"), "ref="+dir.resolve("ref.fa"),
				"out="+dir.resolve("fixed.fa"), "ow=t", "append=f", "t=1", "mindepth=0");
		check(child.exit!=0, spec.name+" should reject mindepth=0 before output: "+child);
		check(child.err.contains("AssertionError"),
				spec.name+" should preserve the positive mindepth parse assertion: "+literal(child.err));
		check(child.err.contains("consensus.FixScaffoldGaps.parse"),
				spec.name+" should reject mindepth=0 in FixScaffoldGaps.parse: "+literal(child.err));
		check(!Files.exists(dir.resolve("fixed.fa")),
				spec.name+" should not create an output FASTA after parse rejection.");
		return new CaseResult(spec.name, spec.expected, child.err, child.out, child.err);
	}

	private CaseResult lilypadZeroAcceptedPairs(final CaseSpec spec) throws Exception{
		final Path dir=caseDir(spec.name);
		writeFasta(dir.resolve("ref.fa"), "ref", "ACGTACGT");
		writeSam(dir.resolve("empty.sam"), 8);
		final Child child=runMain(dir, "consensus.Lilypad",
				"in="+dir.resolve("empty.sam"), "ref="+dir.resolve("ref.fa"),
				"out="+dir.resolve("lilypad.fa"), "ow=t", "append=f", "t=1");
		check(child.exit==0, spec.name+" expected exit0: "+child);
		if(phase==Phase.BASELINE){
			check(child.err.contains("NaN"),
					spec.name+" baseline should expose NaN Average Insert: "+literal(child.err));
		}else{
			check(!child.err.contains("NaN") && !child.err.contains("Infinity"),
					spec.name+" fixed mode must keep Average Insert finite: "+literal(child.err));
		}
		return new CaseResult(spec.name, spec.expected, child.err, child.out, child.err);
	}

	private CaseResult lilypadNoMd(final CaseSpec spec) throws Exception{
		final Path dir=caseDir(spec.name);
		writeFasta(dir.resolve("ref.fa"), "ref", "ACGTACGT");
		writeSam(dir.resolve("mapped.sam"), 8, samNoMd("nomd0", 0, 1, 60, "4M", 0, 0, "ACGT"));
		final Child child=runMain(dir, "consensus.Lilypad",
				"in="+dir.resolve("mapped.sam"), "ref="+dir.resolve("ref.fa"),
				"out="+dir.resolve("lilypad.fa"), "ow=t", "append=f", "t=1");
		if(phase==Phase.BASELINE){
			check(child.exit!=0, spec.name+" baseline should reject M-only SAM without MD: "+child);
			check(child.err.contains("no MD tag and no ScafMap loaded"),
					spec.name+" baseline should expose missing ScafMap publication: "+literal(child.err));
			return new CaseResult(spec.name, spec.expected, child.err, child.out, child.err);
		}
		final String fasta=read(dir.resolve("lilypad.fa"));
		final String expected=">ref\nACGTACGT\n";
		check(child.exit==0, spec.name+" expected exit0: "+child);
		check(expected.equals(fasta), spec.name+" should resolve M-only SAM with ref=; expected "
				+literal(expected)+" observed "+literal(fasta));
		return new CaseResult(spec.name, spec.expected, fasta, child.out, child.err);
	}

	private CaseResult lilypadMultiworkerReciprocalJoin(final CaseSpec spec) throws Exception{
		final Path dir=caseDir(spec.name);
		write(dir.resolve("ref.fa"), ">left\n"+repeat('A', 40)+"\n>right\n"+repeat('C', 40)+"\n");
		writeLilypadMultiworkerSam(dir.resolve("mapped.sam"));
		final int gap=(phase==Phase.BASELINE ? 69 : 68);
		final String expected=">left\n"+repeat('A', 40)+repeat('N', gap)+repeat('C', 40)+"\n";
		final Child single=runMain(dir, "consensus.Lilypad",
				"in="+dir.resolve("mapped.sam"), "ref="+dir.resolve("ref.fa"),
				"out="+dir.resolve("lilypad1.fa"), "ow=t", "append=f", "t=1");
		final String singleFasta=read(dir.resolve("lilypad1.fa"));
		check(single.exit==0, spec.name+" single-worker control expected exit0: "+single);
		check(expected.equals(singleFasta), spec.name+" single-worker control expected "
				+literal(expected)+" observed "+literal(singleFasta));
		if(phase==Phase.BASELINE){
			return new CaseResult(spec.name, spec.expected, singleFasta, single.out, single.err);
		}
		final Child multi=runMain(dir, "consensus.Lilypad",
				"in="+dir.resolve("mapped.sam"), "ref="+dir.resolve("ref.fa"),
				"out="+dir.resolve("lilypad4.fa"), "ow=t", "append=f", "t=4");
		final String multiFasta=read(dir.resolve("lilypad4.fa"));
		check(multi.exit==0, spec.name+" multiworker control expected exit0: "+multi);
		check(expected.equals(multiFasta), spec.name+" four workers must retain the reciprocal mixed-pair join; expected "
				+literal(expected)+" observed "+literal(multiFasta));
		return new CaseResult(spec.name, spec.expected, multiFasta, multi.out, multi.err);
	}

	private CaseResult consensusEmpty(final CaseSpec spec) throws Exception{
		final Path dir=caseDir(spec.name);
		writeFasta(dir.resolve("ref.fa"), "ref", "ACGTACGT");
		writeSam(dir.resolve("empty.sam"), 8);
		final Child child=runMain(dir, "consensus.ConsensusMaker",
				"in="+dir.resolve("empty.sam"), "ref="+dir.resolve("ref.fa"),
				"out="+dir.resolve("consensus.fa"), "ow=t", "append=f", "t=1");
		check(child.exit==0, spec.name+" expected exit0: "+child);
		if(phase==Phase.BASELINE){
			check(child.err.contains("NaN"),
					"Baseline empty input should expose NaN Avg Identity; stderr="+literal(child.err));
		}else{
			check(!child.err.contains("NaN") && !child.err.contains("Infinity"),
					"Fixed empty input must not print non-finite Avg Identity; stderr="+literal(child.err));
		}
		return new CaseResult(spec.name, spec.expected, child.err, child.out, child.err);
	}

	private CaseResult consensusPerfect(final CaseSpec spec) throws Exception{
		final Path dir=caseDir(spec.name);
		writeFasta(dir.resolve("ref.fa"), "ref", "ACGTACGT");
		writeSam(dir.resolve("mapped.sam"), 8, sam("perfect0", 0, 1, 60, "8M", 0, 0, "ACGTACGT"));
		final Child child=runMain(dir, "consensus.ConsensusMaker",
				"in="+dir.resolve("mapped.sam"), "ref="+dir.resolve("ref.fa"),
				"out="+dir.resolve("consensus.fa"), "name=consensus", "ow=t", "append=f", "t=1");
		final String fasta=read(dir.resolve("consensus.fa"));
		final String expected=">consensus\nACGTACGT\n";
		check(child.exit==0, spec.name+" expected exit0: "+child);
		check(expected.equals(fasta), spec.name+" should keep a perfect alignment unchanged; expected "
				+literal(expected)+" observed "+literal(fasta));
		check(!child.err.contains("NaN") && !child.err.contains("Infinity"),
				"Perfect alignment should report finite stderr counters; stderr="+literal(child.err));
		return new CaseResult(spec.name, spec.expected, fasta, child.out, child.err);
	}

	private CaseResult consensusFailingModel(final CaseSpec spec) throws Exception{
		final Path dir=caseDir(spec.name);
		writeFasta(dir.resolve("model_ref.fa"), "model_ref", "NNNNNNNN");
		writeFasta(dir.resolve("ref.fa"), "ref", "ACGTACGT");
		writeSam(dir.resolve("empty.sam"), 8);
		final Path model=dir.resolve("undefined_range.alm");
		final Child modelChild=runMain(dir, "consensus.ConsensusMaker",
				"in="+dir.resolve("empty.sam"), "ref="+dir.resolve("model_ref.fa"), "outmodel="+model,
				"ow=t", "append=f", "t=1");
		check(modelChild.exit==0, spec.name+" expected model creation exit0: "+modelChild);
		writeRepeatedPerfectSam(dir.resolve("mapped.sam"), 8, 4096);
		final Child child=runMain(dir, "consensus.ConsensusMaker",
				"in="+dir.resolve("mapped.sam"), "ref="+dir.resolve("ref.fa"), "inmodel="+model,
				"out="+dir.resolve("consensus.fa"), "ow=t", "append=f", "t=1");
		if(phase==Phase.BASELINE){
			check(child.exit==0, spec.name+" baseline should not throw on the undefined range: "+child);
		}else{
			check(child.exit!=0, spec.name+" fixed mode should report a model-scoring worker failure: "+child);
			check(child.err.contains("ArithmeticException"),
					spec.name+" fixed mode should preserve the normalizedScore exception: "+literal(child.err));
		}
		return new CaseResult(spec.name, spec.expected, child.err, child.out, child.err);
	}

	private CaseResult fixGapsNoMd(final CaseSpec spec) throws Exception{
		final Path dir=caseDir(spec.name);
		writeFasta(dir.resolve("ref.fa"), "ref", "ACGTACGT");
		writeSam(dir.resolve("mapped.sam"), 8, samNoMd("nomd0", 99, 1, 60, "4M", 5, 8, "ACGT"));
		final Child child=runMain(dir, "consensus.FixScaffoldGaps",
				"in="+dir.resolve("mapped.sam"), "ref="+dir.resolve("ref.fa"),
				"out="+dir.resolve("fixed.fa"), "ow=t", "append=f", "t=1");
		if(phase==Phase.BASELINE){
			check(child.exit!=0, spec.name+" baseline should reject M-only SAM without MD: "+child);
			check(child.err.contains("no MD tag and no ScafMap loaded"),
					spec.name+" baseline should expose missing ScafMap publication: "+literal(child.err));
			return new CaseResult(spec.name, spec.expected, child.err, child.out, child.err);
		}
		final String fasta=read(dir.resolve("fixed.fa"));
		final String expected=">ref\nACGTACGT\n";
		check(child.exit==0, spec.name+" expected exit0: "+child);
		check(expected.equals(fasta), spec.name+" should resolve M-only SAM with ref=; expected "
				+literal(expected)+" observed "+literal(fasta));
		return new CaseResult(spec.name, spec.expected, fasta, child.out, child.err);
	}

	private CaseResult insertListAllScaffolds(final CaseSpec spec) throws Exception{
		final Path dir=caseDir(spec.name);
		write(dir.resolve("ref.fa"), ">first\nACGTACGT\n>second\nTTTTGGGG\n");
		write(dir.resolve("mapped.sam"), String.format("@HD\tVN:1.6\tSO:unsorted%n")+
				"@SQ\tSN:first\tLN:8\n@SQ\tSN:second\tLN:8\n"+
				sam("first0", 99, "first", 1, 60, "4M", 5, 8, "ACGT")+
				sam("second0", 99, "second", 1, 60, "4M", 5, 8, "TTTT"));
		final Path insertList=dir.resolve("inserts.tsv");
		final Child child=runMain(dir, "consensus.FixScaffoldGaps",
				"in="+dir.resolve("mapped.sam"), "ref="+dir.resolve("ref.fa"),
				"out="+dir.resolve("fixed.fa"), "insertlist="+insertList, "ow=t", "append=f", "t=1");
		final String actual=read(insertList);
		final String secondOnly=insertRows(8);
		final String expected=(phase==Phase.BASELINE ? secondOnly : insertRows(8)+secondOnly);
		check(child.exit==0, spec.name+" expected exit0: "+child);
		check(expected.equals(actual), spec.name+" should retain insert diagnostics for both scaffolds; expected "
				+literal(expected)+" observed "+literal(actual));
		return new CaseResult(spec.name, spec.expected, actual, child.out, child.err);
	}

	private CaseResult runFixGapsExpecting(final CaseSpec spec, final Path dir, final String baselineBases,
			final String fixedBases, final String... extra) throws Exception{
		final List<String> args=new ArrayList<String>();
		args.add("in="+dir.resolve("mapped.sam"));
		args.add("ref="+dir.resolve("ref.fa"));
		args.add("out="+dir.resolve("fixed.fa"));
		args.add("ow=t");
		args.add("append=f");
		args.add("t=1");
		args.addAll(Arrays.asList(extra));
		final Child child=runMain(dir, "consensus.FixScaffoldGaps", args.toArray(new String[args.size()]));
		final String fasta=read(dir.resolve("fixed.fa"));
		final String expected=">ref\n"+(phase==Phase.BASELINE ? baselineBases : fixedBases)+"\n";
		check(child.exit==0, spec.name+" expected exit0: "+child);
		check(expected.equals(fasta), spec.name+" expected "+literal(expected)+" observed "+literal(fasta));
		return new CaseResult(spec.name, spec.expected, fasta, child.out, child.err);
	}

	private Child runMain(final Path dir, final String mainClass, final String... args) throws Exception{
		final String childName="child"+(++childSerial)+"."+mainClass.replace('.', '_');
		final List<String> command=new ArrayList<String>();
		command.add(new File(System.getProperty("java.home"), "bin/java").getPath());
		command.add("-Xmx256m");
		command.add("-Xms32m");
		command.add("-XX:ActiveProcessorCount=1");
		command.add("-ea");
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(mainClass);
		command.addAll(Arrays.asList(args));
		command.add("fastawrap=0");
		command.add("simd=f");
		final Process process=new ProcessBuilder(command).directory(dir.toFile()).start();
		final ByteArrayOutputStream stdout=new ByteArrayOutputStream();
		final ByteArrayOutputStream stderr=new ByteArrayOutputStream();
		final Thread outThread=copyThread(process.getInputStream(), stdout, mainClass+" stdout");
		final Thread errThread=copyThread(process.getErrorStream(), stderr, mainClass+" stderr");
		final boolean finished=process.waitFor(30, TimeUnit.SECONDS);
		if(!finished){
			process.destroyForcibly();
			process.waitFor();
		}
		outThread.join(5000);
		errThread.join(5000);
		check(!outThread.isAlive() && !errThread.isAlive(),
				"Child stream pumps must finish after process exit for "+mainClass);
		final Child child=new Child(finished ? process.exitValue() : -1,
				new String(stdout.toByteArray(), StandardCharsets.UTF_8),
				new String(stderr.toByteArray(), StandardCharsets.UTF_8), command);
		write(dir.resolve(childName+".command.txt"), child.command.toString()+"\n");
		write(dir.resolve(childName+".stdout.txt"), child.out);
		write(dir.resolve(childName+".stderr.txt"), child.err);
		check(finished, "Child did not finish within 30 seconds for "+mainClass+": "+child);
		return child;
	}

	private static Thread copyThread(final InputStream in, final ByteArrayOutputStream out, final String name){
		final Thread thread=new Thread(new Runnable(){
			@Override
			public void run(){
				final byte[] buffer=new byte[8192];
				int len;
				try{
					while((len=in.read(buffer))>=0){out.write(buffer, 0, len);}
				}catch(final IOException ioe){throw new RuntimeException(ioe);}
			}
		}, name);
		thread.setDaemon(true);
		thread.start();
		return thread;
	}

	private Path caseDir(final String name) throws IOException{
		return Files.createDirectory(root.resolve(name));
	}

	private static String sam(final String name, final int flag, final int pos, final int mapq,
			final String cigar, final int matePos, final int tlen, final String bases){
		return sam(name, flag, "ref", pos, mapq, cigar, matePos, tlen, bases);
	}

	private static String sam(final String name, final int flag, final String rname, final int pos,
			final int mapq, final String cigar, final int matePos, final int tlen, final String bases){
		return sam(name, flag, rname, pos, mapq, cigar, matePos, tlen, bases, true);
	}

	private static String samNoMd(final String name, final int flag, final int pos, final int mapq,
			final String cigar, final int matePos, final int tlen, final String bases){
		return sam(name, flag, "ref", pos, mapq, cigar, matePos, tlen, bases, false);
	}

	private static String sam(final String name, final int flag, final String rname, final int pos,
			final int mapq, final String cigar, final int matePos, final int tlen, final String bases,
			final boolean mdTag){
		return name+'\t'+flag+'\t'+rname+'\t'+pos+'\t'+mapq+'\t'+cigar+'\t'+"="+'\t'+matePos+'\t'+tlen+
				'\t'+bases+'\t'+repeat('I', bases.length())+(mdTag ? "\tMD:Z:"+bases.length() : "")+'\n';
	}

	private static void writeDeepInsertSam(final Path path, final int refLength,
			final int records) throws IOException{
		final String[] lines=new String[records];
		for(int i=0; i<records; i++){
			lines[i]=sam("deep"+i, 99, 551, 60, "10M", 651, 100, repeat('C', 10));
		}
		writeSam(path, refLength, lines);
	}

	private static void writeSameElevenRecordGapSam(final Path path, final int refLength,
			final int deepInsert, final int gapInsert) throws IOException{
		final List<String> records=new ArrayList<String>();
		for(int i=0; i<10; i++){
			records.add(sam("deep"+i, 99, 551, 60, "10M", 551+deepInsert, deepInsert, repeat('C', 10)));
		}
		final int gapPos=(gapInsert>20 ? 321 : 346);
		records.add(sam("gap0", 99, gapPos, 60, "10M", gapPos+gapInsert, gapInsert, repeat('A', 10)));
		writeSam(path, refLength, records.toArray(new String[records.size()]));
	}

	private static void writeLilypadMultiworkerSam(final Path path) throws IOException{
		final StringBuilder sb=new StringBuilder("@HD\tVN:1.6\tSO:unsorted\n"+
				"@SQ\tSN:left\tLN:40\n@SQ\tSN:right\tLN:40\n");
		for(int i=0; i<1100; i++){
			sb.append(sam("insert"+i, 99, "left", 1, 60, "10M", 91, 100, repeat('A', 10)));
		}
		for(int i=0; i<1024; i++){
			final int pnext=(i&7)+1;
			sb.append("mixed_lr").append(i).append('\t').append(99).append('\t')
					.append("left").append('\t').append(20).append('\t').append(60).append('\t')
					.append("10M\tright\t").append(pnext).append("\t0\t")
					.append("AAAAAAAAAA\tIIIIIIIIII\tMD:Z:10\n");
			sb.append("mixed_rl").append(i).append('\t').append(83).append('\t')
					.append("right").append('\t').append(8).append('\t').append(60).append('\t')
					.append("10M\tleft\t").append(20+(i&7)).append("\t0\t")
					.append("CCCCCCCCCC\tIIIIIIIIII\tMD:Z:10\n");
		}
		write(path, sb.toString());
	}

	private static void writeRepeatedPerfectSam(final Path path, final int refLength,
			final int records) throws IOException{
		final String[] lines=new String[records];
		for(int i=0; i<records; i++){
			lines[i]=sam("perfect"+i, 0, 1, 60, "8M", 0, 0, "ACGTACGT");
		}
		writeSam(path, refLength, lines);
	}

	private static void writeSam(final Path path, final int refLength, final String... records) throws IOException{
		final StringBuilder sb=new StringBuilder("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ref\tLN:");
		sb.append(refLength).append('\n');
		for(final String record : records){sb.append(record);}
		write(path, sb.toString());
	}

	private static void writeFasta(final Path path, final String name, final String bases) throws IOException{
		write(path, '>'+name+'\n'+bases+'\n');
	}

	private void writeObservations(final List<CaseResult> results) throws IOException{
		final StringBuilder sb=new StringBuilder("case\texpected\tobserved\tstdout\tstderr\n");
		for(final CaseResult result : results){
			sb.append(tsv(result.name)).append('\t').append(tsv(result.expected)).append('\t')
					.append(tsv(result.observed)).append('\t').append(tsv(result.stdout)).append('\t')
					.append(tsv(result.stderr)).append('\n');
		}
		write(root.resolve("observations.tsv"), sb.toString());
	}

	private static String insertRows(final int length){
		final StringBuilder sb=new StringBuilder(length*6);
		for(int i=0; i<length; i++){
			final int depth=(i>0 && i<length-1 ? 1 : 0);
			final int insert=(depth>0 ? length : 0);
			sb.append(i).append('\t').append(depth).append('\t').append(insert).append('\n');
		}
		return sb.toString();
	}

	private static String read(final Path path) throws IOException{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}

	private static void write(final Path path, final String text) throws IOException{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	private static String repeat(final char c, final int length){
		final char[] array=new char[length];
		Arrays.fill(array, c);
		return new String(array);
	}

	private static String tsv(final String text){
		return text.replace("\\", "\\\\").replace("\t", "\\t").replace("\r", "\\r").replace("\n", "\\n");
	}

	private static String literal(final String text){return tsv(text);}

	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private abstract static class CaseSpec{
		CaseSpec(final String name_, final String expected_){
			name=name_;
			expected=expected_;
		}
		abstract CaseResult run(ConsensusCliBoundaryTest test) throws Exception;
		final String name;
		final String expected;
	}

	private enum Phase{
		BASELINE("baseline"),
		FIXED("fixed");

		static Phase parse(final String value){
			check(value!=null, "Set phase=baseline or phase=fixed.");
			final String lower=value.toLowerCase(Locale.ROOT);
			for(final Phase phase : values()){
				if(phase.argument.equals(lower)){return phase;}
			}
			throw new IllegalArgumentException("Unknown phase: "+value);
		}

		Phase(final String argument_){
			argument=argument_;
		}
		final String argument;
	}

	private static final class CaseResult{
		CaseResult(final String name_, final String expected_, final String observed_,
				final String stdout_, final String stderr_){
			name=name_;
			expected=expected_;
			observed=observed_;
			stdout=stdout_;
			stderr=stderr_;
		}
		final String name;
		final String expected;
		final String observed;
		final String stdout;
		final String stderr;
	}

	private static final class Child{
		Child(final int exit_, final String out_, final String err_, final List<String> command_){
			exit=exit_;
			out=out_;
			err=err_;
			command=command_;
		}
		@Override
		public String toString(){
			return "exit="+exit+" command="+command+" stdout="+literal(out)+" stderr="+literal(err);
		}
		final int exit;
		final String out;
		final String err;
		final List<String> command;
	}

	private final Path root;
	private final Phase phase;
	private int childSerial=0;

	private static final CaseSpec[] CASES={
		new CaseSpec("fixgaps_no_ns_control", ">ref\\nACGTACGT\\n"){
			@Override
			CaseResult run(final ConsensusCliBoundaryTest test) throws Exception{return test.noNsUnchanged(this);}
		},
		new CaseSpec("fixgaps_short_internal_ns", ">ref\\nAAANNNCCC\\n"){
			@Override
			CaseResult run(final ConsensusCliBoundaryTest test) throws Exception{return test.shortInternalNs(this);}
		},
		new CaseSpec("fixgaps_edge_proximal_ns", "preserve 301A+10N+C"){
			@Override
			CaseResult run(final ConsensusCliBoundaryTest test) throws Exception{return test.edgeProximalNs(this);}
		},
		new CaseSpec("fixgaps_terminal_ns", ">ref\\nAAA+10N\\n"){
			@Override
			CaseResult run(final ConsensusCliBoundaryTest test) throws Exception{return test.terminalNs(this);}
		},
		new CaseSpec("fixgaps_mindepth_unsupported", "gap unchanged below mindepth=10"){
			@Override
			CaseResult run(final ConsensusCliBoundaryTest test) throws Exception{return test.minDepthUnsupportedGap(this);}
		},
		new CaseSpec("fixgaps_zero_accepted_pairs", "finite summary with no spanning pairs"){
			@Override
			CaseResult run(final ConsensusCliBoundaryTest test) throws Exception{return test.zeroAcceptedPairs(this);}
		},
		new CaseSpec("fixgaps_mindepth_supported", "supported gap widens under mindepth=1"){
			@Override
			CaseResult run(final ConsensusCliBoundaryTest test) throws Exception{return test.minDepthSupportedGap(this);}
		},
		new CaseSpec("fixgaps_minimum_clamp_stats", "clamped raw shrink leaves a gap unchanged"){
			@Override
			CaseResult run(final ConsensusCliBoundaryTest test) throws Exception{return test.minimumClampStats(this);}
		},
		new CaseSpec("fixgaps_m_cigar_no_md", "M-only SAM records resolve through ref="){
			@Override
			CaseResult run(final ConsensusCliBoundaryTest test) throws Exception{return test.fixGapsNoMd(this);}
		},
		new CaseSpec("fixgaps_zero_mindepth_rejected", "mindepth=0 is rejected"){
			@Override
			CaseResult run(final ConsensusCliBoundaryTest test) throws Exception{return test.zeroMinDepthRejected(this);}
		},
		new CaseSpec("lilypad_zero_accepted_pairs", "finite summary with no spanning pairs"){
			@Override
			CaseResult run(final ConsensusCliBoundaryTest test) throws Exception{return test.lilypadZeroAcceptedPairs(this);}
		},
		new CaseSpec("lilypad_m_cigar_no_md", "M-only SAM records resolve through ref="){
			@Override
			CaseResult run(final ConsensusCliBoundaryTest test) throws Exception{return test.lilypadNoMd(this);}
		},
		new CaseSpec("lilypad_multiworker_reciprocal_join", "four workers retain reciprocal mixed-pair join"){
			@Override
			CaseResult run(final ConsensusCliBoundaryTest test) throws Exception{
				return test.lilypadMultiworkerReciprocalJoin(this);
			}
		},
		new CaseSpec("consensus_empty_input_finite", "no NaN/Infinity in stderr"){
			@Override
			CaseResult run(final ConsensusCliBoundaryTest test) throws Exception{return test.consensusEmpty(this);}
		},
		new CaseSpec("consensus_perfect_control", ">consensus\\nACGTACGT\\n"){
			@Override
			CaseResult run(final ConsensusCliBoundaryTest test) throws Exception{return test.consensusPerfect(this);}
		},
		new CaseSpec("consensus_failing_model_exits", "worker scoring exception exits clearly"){
			@Override
			CaseResult run(final ConsensusCliBoundaryTest test) throws Exception{return test.consensusFailingModel(this);}
		},
		new CaseSpec("fixgaps_insertlist_two_scaffolds", "both scaffolds' insert rows"){
			@Override
			CaseResult run(final ConsensusCliBoundaryTest test) throws Exception{return test.insertListAllScaffolds(this);}
		},
	};
}
