package test.dna;

import java.lang.reflect.Method;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.regex.Pattern;

import dna.Data;
import structures.ByteBuilder;

/** Finite exit-status checks for Data's subprocess-availability probe.
 * The executor supplies a saved shell fixture that prints "stdout\n" to stdout,
 * then "stderr\n" to stderr, then exits with its first argument as the status.
 * Invoke that fixture with /bin/sh; it must not spawn persistent children or wait
 * for input. This test writes no files and does not change tool-availability caches.
 * Repeated fast exits expose a race when the probe skips waiting for an exited child;
 * a passing baseline does not prove that scheduling-dependent defect absent.
 * @author Sucrose
 */
public final class DataProcessProbeTest{

	/** Requires case=controls|reject script=path; accepts capture=t|f and repeats=1..64. */
	public static void main(String[] args) throws Exception{
		String mode=null;
		Path script=null;
		boolean capture=true;
		int repeats=32;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("script=")){script=Paths.get(arg.substring(7)).toAbsolutePath();}
			else if(arg.equals("capture=t")){capture=true;}
			else if(arg.equals("capture=f")){capture=false;}
			else if(arg.startsWith("repeats=")){repeats=Integer.parseInt(arg.substring(8));}
			else{throw new IllegalArgumentException(arg);}
		}
		require(("controls".equals(mode) || "reject".equals(mode)) && script!=null,
				"Require case=controls|reject script=<saved-finite-shell-fixture>");
		require(repeats>=1 && repeats<=64, "Bounded probe count must be in1..64: "+repeats);
		require(Files.isRegularFile(script) && Files.isReadable(script), "Unreadable fixture: "+script);
		require(!Pattern.compile("\\s").matcher(script.toString()).find(),
				"Data.testExecute splits on whitespace; fixture path must contain none: "+script);
		require(Files.isExecutable(Paths.get("/bin/sh")), "This fixture requires /bin/sh");

		final Method probe=Data.class.getDeclaredMethod("testExecute", String.class, ByteBuilder.class);
		probe.setAccessible(true);
		final int[] statuses="controls".equals(mode) ? new int[]{0, 1, 125} : new int[]{126, 127};
		final ByteBuilder output=capture ? new ByteBuilder() : null;
		int mismatches=0;
		for(int status : statuses){
			final int expected=status<126 ? 1 : -1;
			int accepted=0, rejected=0;
			for(int trial=0; trial<repeats; trial++){
				if(output!=null){output.clear(); output.append("stale buffer contents");}
				final int actual=(Integer)probe.invoke(null, "/bin/sh "+script+" "+status, output);
				require(actual==1 || actual==-1, "Unexpected probe result: "+actual);
				if(actual==1){accepted++;}else{rejected++;}
				if(actual!=expected){mismatches++;}
				if(output!=null){
					require(output.toString().equals("stdout\nstderr\n"),
							"Probe must clear the buffer and capture both streams; status="+status+" trial="+trial);
				}
			}
			System.out.println("exit="+status+" capture="+capture+" trials="+repeats+
					" accepted="+accepted+" rejected="+rejected+" expected="+expected);
		}
		require(mismatches==0, "Availability disagreed with the child exit status in "+mismatches+" trials");
		System.out.println("PASS "+mode+" capture="+capture);
	}

	/** Keeps the oracle active with assertions either enabled or disabled. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
