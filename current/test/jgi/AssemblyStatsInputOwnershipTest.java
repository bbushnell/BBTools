package test.jgi;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;

import jgi.AssemblyStats2;
import parse.PreParser;
import shared.Shared;

/** Checks that assembly statistics borrow standard input, including filename aliases.
 * Uses real statistics over one literal contig and observes close calls directly.
 * Run in an isolated JVM because native parsers retain process-wide settings.
 * The original System.in is restored on success and failure; no real input is read.
 * @author Brian Bushnell, Shinobu */
public final class AssemblyStatsInputOwnershipTest{

	/** Runs exact stdin and the stdin.fa alias; accepts no program options. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, AssemblyStatsInputOwnershipTest.class, false);
		try{
			require(pp.args.length==0, "Standard-input ownership tests take no program options.");
			check("stdin");
			check("stdin.fa");
			System.out.println("ASSEMBLY_STATS_STDIN_OWNERSHIP_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Static fixture entry point only. */
	private AssemblyStatsInputOwnershipTest(){}

	/** Reads the literal input through the supplied spelling without closing its owner. */
	private static void check(final String name){
		require(name.equals("stdin") || !new File(name).exists(),
				"A real file would take precedence over the standard-input alias: "+name);
		final InputStream original=System.in;
		final TrackedInput input=new TrackedInput();
		try{
			System.setIn(input);
			AssemblyStats2.main(new String[]{"in="+name, "format=0", "showspeed=f", "showbbmap=f"});
			require(!input.closed, "Assembly statistics closed borrowed standard input: "+name);
			require(input.available()==0, "Statistics did not consume the supplied contig: "+name);
			require(AssemblyStats2.lastSize==16 && AssemblyStats2.lastContigs==1 &&
					AssemblyStats2.lastL50==16 && AssemblyStats2.lastL90==16 && AssemblyStats2.lastMaxContig==16,
					"Expected fresh statistics for one literal 16-base contig: "+name);
			System.out.println("ASSEMBLY_STATS_STDIN_OPEN\t"+name);
		}finally{System.setIn(original);}
	}

	/** Enforces expectations under either JVM assertion setting. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** In-memory input whose otherwise harmless close operation is observable. */
	private static final class TrackedInput extends ByteArrayInputStream{
		/** Supplies one known, fully defined contig without filesystem input. */
		TrackedInput(){super(">fixture\nACGTACGTACGTACGT\n".getBytes(StandardCharsets.US_ASCII));}
		/** Records ownership violations even though ByteArrayInputStream remains readable. */
		@Override
		public void close() throws IOException{closed=true; super.close();}
		/** Whether the statistics invocation attempted to close this borrowed stream. */
		boolean closed;
	}
}
