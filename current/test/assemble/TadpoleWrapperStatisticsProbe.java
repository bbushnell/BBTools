package test.assemble;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.Random;

import assemble.TadpoleWrapper;
import fileIO.ReadWrite;
import jgi.AssemblyStats2;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Observes Wrapper statistics freshness after a real earlier statistics invocation.
 * A single 16-base prior contig seeds the public last* fields; Wrapper then assembles
 * four copies of a known 512-base read in the same JVM. The emitted FASTA is checked
 * independently of those fields. Ordinary mode records freshness rather than
 * asserting it; optional state checks require quiet reporting without leaked settings.
 * This does not demonstrate selection of a wrong K.
 * Each invocation requires a fresh directory and an isolated JVM because both tools
 * retain global settings. Use the development launcher on the normal validation host.
 * @author Brian Bushnell, Shinobu
 */
public final class TadpoleWrapperStatisticsProbe{

	/** Static fixture entry point only. */
	private TadpoleWrapperStatisticsProbe(){}

	/** Accepts outdir=NEW_DIRECTORY, showstats=t/f, and optional statecheck=t/statsalias=t.
	 * State checks exercise quiet reporting and temporary global restoration; the alias
	 * option sends stats= instead of showstats= to Wrapper. Ordinary observations retain
	 * their original input/settings so the saved baseline remains comparable.
	 */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, TadpoleWrapperStatisticsProbe.class, false);
		try{
			Path root=null;
			boolean showStats=true;
			boolean stateCheck=false, statsAlias=false;
			for(String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String name=pair[0], value=pair.length>1 ? pair[1] : null;
				if(name.equalsIgnoreCase("outdir")){
					require(value!=null && !value.isEmpty(), "outdir requires a new fixture directory.");
					root=Paths.get(value).toAbsolutePath();
				}else if(name.equalsIgnoreCase("showstats")){showStats=Parse.parseBoolean(value);}
				else if(name.equalsIgnoreCase("statecheck")){stateCheck=Parse.parseBoolean(value);}
				else if(name.equalsIgnoreCase("statsalias")){statsAlias=Parse.parseBoolean(value);}
				else{throw new IllegalArgumentException("Unknown probe option: "+arg);}
			}
			require(root!=null, "Specify outdir=NEW_DIRECTORY; existing directories are refused.");
			require(!stateCheck || !showStats, "State checks require showstats=f to exercise quiet reporting.");
			Files.createDirectory(root);
			run(root, showStats, stateCheck, statsAlias);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Uses actual prior statistics, then records the fields consumed by Wrapper after its native invocation. */
	private static void run(final Path root, final boolean showStats, final boolean stateCheck,
			final boolean statsAlias) throws Exception{
		final Path prior=root.resolve("prior.fa"), input=root.resolve("reads.fa");
		final Path output=root.resolve("assembly-31.fa"), report=root.resolve("observations.tsv");
		Files.write(prior, ">prior\nACGTACGTACGTACGT\n".getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		final byte[] original=fixture();
		Files.write(input, original, StandardOpenOption.CREATE_NEW);
		AssemblyStats2.main(new String[]{"in="+prior, "out="+root.resolve("prior-statistics.txt"), "ow=f"});
		require(metricsMatch(16), "The prior statistics fixture must establish one 16-base contig.");
		final ByteArrayOutputStream diagnostics=stateCheck ? new ByteArrayOutputStream() : null;
		final PrintStream stateStream=stateCheck ? new PrintStream(diagnostics) : null;
		final PrintStream previousStream=AssemblyStats2.outstream;
		if(stateCheck){
			AssemblyStats2.FORMAT=3; AssemblyStats2.showspeed=true;
			AssemblyStats2.printExtended=false; AssemblyStats2.printL90=false;
			AssemblyStats2.outstream=stateStream;
		}
		RuntimeException failure=null;
		try{
			final String[] command=new String[]{"in="+input, "out="+root.resolve("assembly-%.fa"),
					"k=31", "t=1", "buildthreads=1", "prealloc=f", "prefilter=0", "minprob=0",
					"mcs=1", "mce=1", "mincontig=1", "minextension=1", "trimends=0", "shave=f",
					"rinse=f", "ow=f", "interleaved=f", "assemblydefaults=f", "branchrescue=f",
					(statsAlias ? "stats=" : "showstats=")+showStats};
			final String[] invocation=stateCheck ? Arrays.copyOf(command, command.length+2) : command;
			if(stateCheck){invocation[command.length]="pigz=f"; invocation[command.length+1]="unpigz=f";}
			TadpoleWrapper.main(invocation);
			if(stateCheck){
				require(AssemblyStats2.FORMAT==3 && AssemblyStats2.showspeed &&
						!AssemblyStats2.printExtended && !AssemblyStats2.printL90 && AssemblyStats2.outstream==stateStream,
						"Quiet statistics leaked temporary report settings or replaced the borrowed stream.");
				require(!ReadWrite.USE_PIGZ && !ReadWrite.USE_UNPIGZ, "Statistics overrode the trial's explicit compression policy.");
				require(diagnostics.size()==0, "Suppressed statistics wrote diagnostics to the borrowed stream.");
				stateStream.print("still open"); stateStream.flush();
				require(!stateStream.checkError() && diagnostics.size()==10, "Statistics closed its borrowed diagnostic stream.");
				System.out.println("WRAPPER_STATS_STATE_PASS alias="+statsAlias);
			}
		}catch(RuntimeException e){failure=e;}
		finally{if(stateCheck){AssemblyStats2.outstream=previousStream; stateStream.close();}}
		require(Arrays.equals(original, Files.readAllBytes(input)), "Wrapper changed its generated read input.");
		final long outputBases=failure==null ? checkedOutputBases(output) : -1;
		final boolean fresh=failure==null && metricsMatch(READ_LENGTH);
		final ByteBuilder bb=new ByteBuilder(512);
		bb.append("showstats\twrapper_completed\toutput_bases\tL50\tL90\ttotal_bases\tcontigs\tmax_contig\tfresh\texception\n");
		bb.append(Boolean.toString(showStats)).tab().append(Boolean.toString(failure==null)).tab();
		bb.append(outputBases).tab().append(AssemblyStats2.lastL50).tab().append(AssemblyStats2.lastL90).tab();
		bb.append(AssemblyStats2.lastSize).tab().append(AssemblyStats2.lastContigs).tab().append(AssemblyStats2.lastMaxContig).tab();
		bb.append(Boolean.toString(fresh)).tab().append(failure==null ? "-" : failure.getClass().getName()).nl();
		Files.write(report, bb.toBytes(), StandardOpenOption.CREATE_NEW);
		if(failure!=null){throw failure;}
		System.out.println("WRAPPER_STATS_OBSERVED fresh="+fresh+" output_bases="+outputBases);
	}

	/** Requires one complete nonempty FASTA record of the independently known fixture length. */
	private static long checkedOutputBases(final Path output) throws Exception{
		require(Files.isRegularFile(output), "Wrapper did not emit the expected K31 assembly.");
		int headers=0;
		long bases=0;
		for(String line:Files.readAllLines(output, StandardCharsets.US_ASCII)){
			if(line.startsWith(">")){headers++;}
			else if(!line.isEmpty()){
				require(headers>0, "Assembly sequence preceded its FASTA header.");
				bases+=line.length();
			}
		}
		require(headers==1 && bases==READ_LENGTH, "Expected one "+READ_LENGTH+"-base contig; headers="+headers+", bases="+bases);
		return bases;
	}

	/** Checks all five fields used by Wrapper against a single known contig, without invoking its comparison heuristic. */
	private static boolean metricsMatch(final long length){
		return AssemblyStats2.lastL50==length && AssemblyStats2.lastL90==length &&
				AssemblyStats2.lastSize==length && AssemblyStats2.lastContigs==1 && AssemblyStats2.lastMaxContig==length;
	}

	/** Builds the same deterministic read sequence used by the existing finite wrapper fixture. */
	private static byte[] fixture(){
		final Random random=new Random(714031);
		final char[] alphabet={'A', 'C', 'G', 'T'}, sequence=new char[READ_LENGTH];
		for(int i=0; i<sequence.length; i++){sequence[i]=alphabet[random.nextInt(4)];}
		final StringBuilder fasta=new StringBuilder();
		for(int i=0; i<4; i++){fasta.append(">sample").append(i).append('\n').append(sequence).append('\n');}
		return fasta.toString().getBytes(StandardCharsets.US_ASCII);
	}

	/** Rejects invalid fixture setup independently of JVM assertion settings. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Known length of each repeated input read and the expected assembled contig. */
	private static final int READ_LENGTH=512;
}
