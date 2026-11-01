package test.parse;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import fileIO.ByteFile;
import parse.Parser;

/** Independent checks for option precedence and configuration-file expansion.
 * Each invocation must run in its own JVM because Parser changes process-wide state.
 * @author Brian Bushnell, Shinobu */
public final class ParserStateTest {

	/** Runs finite literal-oracle checks, reporting every failure before exiting.
	 * @param args No arguments are supported
	 * @throws Exception if temporary fixture creation or deletion fails */
	public static void main(String[] args) throws Exception{
		if(args.length!=0){throw new IllegalArgumentException("No arguments expected");}
		final ParserStateTest test=new ParserStateTest();
		final Path dir=Files.createTempDirectory("parser-state-");
		try{
			final Path plain=dir.resolve("plain.config"), equals=dir.resolve("a=b.config");
			Files.write(plain, Arrays.asList("# ignored", "  t=2  ", "", "in=a=b.fq"), StandardCharsets.UTF_8);
			Files.write(equals, Arrays.asList("t=3"), StandardCharsets.UTF_8);
			for(int before=1; before<=4; before++){
				for(int after=1; after<=4; after++){
					final int first=before, last=after;
					test.check("bf"+first+" then bf"+last, () -> {
						resetModes();
						option("bf"+first, "t");
						option("bf"+last, "t");
						final boolean[] actual={ByteFile.FORCE_MODE_BF1, ByteFile.FORCE_MODE_BF2,
							ByteFile.FORCE_MODE_BF3, ByteFile.FORCE_MODE_BF4};
						for(int i=0; i<4; i++){
							require(actual[i]==(i==last-1), "only last enabled mode may remain: "+Arrays.toString(actual));
						}
					});
				}
			}
			test.check("disabling unrelated reader mode", () -> {
				resetModes(); option("bf3", "t"); option("bf2", "f");
				require(ByteFile.FORCE_MODE_BF3 && !ByteFile.FORCE_MODE_BF2, "false must not select a different mode");
			});
			test.check("config preserves order and values", () -> {
				final String[] actual=Parser.parseConfig(new String[]{"t=1", "config="+plain, "t=4"});
				require(Arrays.equals(actual, new String[]{"t=1", "t=2", "in=a=b.fq", "t=4"}), Arrays.toString(actual));
			});
			test.check("equals in config filename", () -> {
				final String[] actual=Parser.parseConfig(new String[]{"config="+equals});
				require(Arrays.equals(actual, new String[]{"t=3"}), Arrays.toString(actual));
			});
		}finally{
			resetModes();
			Files.deleteIfExists(dir.resolve("plain.config"));
			Files.deleteIfExists(dir.resolve("a=b.config"));
			Files.delete(dir);
		}
		System.out.println("PARSER_STATE_RESULT cases="+test.cases+" failures="+test.failures);
		if(test.failures!=0){throw new AssertionError("Parser state failures="+test.failures);}
	}

	/** Invokes the public entry point and rejects accidentally unrecognized options. */
	private static void option(String name, String value){
		require(Parser.parseStatic(name+"="+value, name, value), "unrecognized option "+name);
	}

	/** Resets only the four reader selection flags exercised by this fixture. */
	private static void resetModes(){
		ByteFile.FORCE_MODE_BF1=ByteFile.FORCE_MODE_BF2=ByteFile.FORCE_MODE_BF3=ByteFile.FORCE_MODE_BF4=false;
	}

	/** Records failures independently of whether Java assertions are enabled. */
	private void check(String name, Runnable action){
		cases++;
		try{action.run();}
		catch(RuntimeException | AssertionError e){
			failures++;
			System.out.println("FAIL "+name+": "+e);
		}
	}

	/** Fails an expectation even when the JVM runs with -da. */
	private static void require(boolean condition, String reason){
		if(!condition){throw new AssertionError(reason);}
	}

	/** Number of attempted independent cases. */
	private int cases;
	/** Number of failed cases. */
	private int failures;
}
