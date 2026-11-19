package test.benchmark;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.nio.file.attribute.FileTime;
import java.util.Arrays;

import benchmark.CallGenesBenchmarkSupport;
import benchmark.CallGenesBenchmarkSupport.InputVerification;
import gff.GffLine;
import prok.ProkObject;

/** Literal live-support contracts, without loading deprecated diagnostic classes or hashing inputs.
 * Execute in the Dori validation build with a fresh fixture directory.
 * @author Sayu
 */
public final class CallGenesBenchmarkSupportTest{

	/** Creates fresh filesystem cases and checks the shared schema, flags, partial sides and CPU clock. */
	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh absolute fixture directory");
		final Path root=Paths.get(args[0]);
		check(root.isAbsolute() && !Files.exists(root), "Fixture directory must be fresh and absolute");
		Files.createDirectory(root);
		schemas();
		cdsFlags();
		partials();
		verification(root);
		final long before=CallGenesBenchmarkSupport.processCpuNanos();
		final long after=CallGenesBenchmarkSupport.processCpuNanos();
		check(before>=0 && after>=before, "Whole-process CPU nanoseconds must be available and monotonic");
		System.out.println("CALLGENES_BENCHMARK_SUPPORT_TEST_PASS checks="+checks+" hashes=0");
	}

	/** Checks the literal input-column and grade-counter orders used by live consumers. */
	private static void schemas(){
		final String header="genome\tfasta\tfasta_sha80\treference\treference_sha80\tnative\tnative_sha80"
			+"\tbaseline\tbaseline_sha80\tannotations\tannotations_sha80\tgrades\tgrades_sha80";
		check(header.equals(CallGenesBenchmarkSupport.INPUT_HEADER), "Literal accepted 13-column input schema");
		final String[] grades={"tp_start_reference", "tp_stop_reference", "fn_start", "fn_stop",
			"tp_start_query", "tp_stop_query", "fp_start", "fp_stop"};
		check(Arrays.equals(grades, CallGenesBenchmarkSupport.GRADE_NAMES), "Literal native eight-counter order");
	}

	/** Checks idempotent CDS-only selection and restores all inherited global call flags. */
	private static void cdsFlags(){
		final boolean[] previous=flags();
		try{
			setFlags(new boolean[]{false, true, true, true, true, true});
			CallGenesBenchmarkSupport.cdsOnly();
			check(Arrays.equals(flags(), new boolean[]{true, false, false, false, false, false}),
				"Reset every inherited RNA flag, including a disabled CDS flag");
			CallGenesBenchmarkSupport.cdsOnly();
			check(Arrays.equals(flags(), new boolean[]{true, false, false, false, false, false}), "CDS-only is idempotent");
		}finally{setFlags(previous);}
	}

	/** Literal biological sides include conflicts, absent flags and strand reversal. */
	private static void partials()throws Exception{
		final String[][] cases={
			{".", "complete", "complete"},
			{"partial=false", "complete", "complete"},
			{"partial=00", "complete", "complete"},
			{"partial=10", "five_prime", "three_prime"},
			{"partial=01", "three_prime", "five_prime"},
			{"partial=11", "both", "both"},
			{"partial=true", "unknown", "unknown"},
			{"partial=true;start_range=.,10", "five_prime", "three_prime"},
			{"partial=true;end_range=99,.", "three_prime", "five_prime"},
			{"start_range=.,10", "five_prime", "three_prime"},
			{"end_range=99,.", "three_prime", "five_prime"},
			{"start_range=10,10;end_range=99,99", "complete", "complete"},
			{"partial=false;start_range=.,10", "unknown", "unknown"},
			{"partial=10;start_range=10,10", "unknown", "unknown"},
			{"start_range=11,11", "unknown", "unknown"},
			{"start_range=.,.", "unknown", "unknown"}
		};
		for(String[] item:cases){
			for(int strand=0; strand<2; strand++){
				check(CallGenesBenchmarkSupport.partial(feature(item[0], strand)).equals(item[strand+1]),
					"Literal partial side: "+item[0]+" strand="+strand);
			}
		}
		final GffLine absent=feature(".", 0);
		absent.attributes=null;
		check(CallGenesBenchmarkSupport.partial(absent).equals("complete"), "Null attributes remain complete");
		check(CallGenesBenchmarkSupport.partial(feature("partial=10", 2)).equals("unknown"), "Unknown strand has no biological side");
		for(String bad:new String[]{"partial=", "partial=12", "partial=true;partial=false",
			"start_range=10", "start_range=,10", "start_range=10,", "start_range=10,10,10",
			"start_range=0,10", "start_range=.,10;start_range=.,10", "end_range=99,-1"}){
			reject(()->CallGenesBenchmarkSupport.partial(feature(bad, 0)), null);
		}
	}

	/** Structural identity rejects metadata changes, replaced leaves and parent retargets. */
	private static void verification(Path root)throws Exception{
		final Path file=write(root.resolve("stable.txt"), "abc");
		check(!new InputVerification(null).structural && new InputVerification(null).mode().equals("sha80"),
			"Default verification remains sha80");
		check(!new InputVerification("sha80").structural, "Explicit sha80 remains accepted");
		reject(()->new InputVerification("none"), "verification must");
		reject(()->new InputVerification(""), "verification must");
		reject(()->new InputVerification(null).pin(file, "unmeasured"), "Input requires sha80");
		final InputVerification stable=bound(file, "unmeasured");
		check(stable.structural && stable.mode().equals("structural"), "Structural mode is explicit");
		stable.pin(file, "unmeasured");
		stable.verifyComplete();
		checks++;
		bound(file, "00000000000000000000").verifyComplete();
		checks++;
		reject(()->stable.pin(file, "00000000000000000000"), "Conflicting input identities");
		reject(()->bound(file, "missing"), "Input requires sha80");
		reject(()->bound(root, "unmeasured"), "Regular nonsymlink");
		final Path link=root.resolve("leaf-link");
		Files.createSymbolicLink(link, file);
		reject(()->bound(link, "unmeasured"), "Regular nonsymlink");
		write(file, "abcd");
		reject(stable::verifyComplete, "Structural input changed");
		reject(()->stable.pin(file, "unmeasured"), "Structural input changed");
		final InputVerification modified=bound(file, "unmeasured");
		final FileTime time=Files.getLastModifiedTime(file);
		Files.setLastModifiedTime(file, FileTime.fromMillis(time.toMillis()+5000));
		reject(modified::verifyComplete, "Structural input changed");
		final InputVerification replaced=bound(file, "unmeasured");
		final Path replacement=write(root.resolve("replacement.txt"), "abcd");
		Files.setLastModifiedTime(replacement, Files.getLastModifiedTime(file));
		Files.move(replacement, file, StandardCopyOption.REPLACE_EXISTING);
		reject(replaced::verifyComplete, "Structural input changed");
		final Path a=Files.createDirectory(root.resolve("a")), b=Files.createDirectory(root.resolve("b"));
		final Path parent=root.resolve("parent-link");
		write(a.resolve("x"), "x");
		Files.createLink(b.resolve("x"), a.resolve("x"));
		Files.createSymbolicLink(parent, a);
		final InputVerification retarget=bound(parent.resolve("x"), "unmeasured");
		Files.delete(parent);
		Files.createSymbolicLink(parent, b);
		reject(retarget::verifyComplete, "Structural input changed");
		final InputVerification symlink=bound(file, "unmeasured");
		Files.move(file, root.resolve("moved.txt"));
		Files.createSymbolicLink(file, root.resolve("moved.txt"));
		reject(symlink::verifyComplete, "Require regular non-symlink");
	}

	/** Creates a structural verifier and captures one declared input identity. */
	private static InputVerification bound(Path file, String label)throws Exception{
		final InputVerification result=new InputVerification("structural");
		result.pin(file, label);
		return result;
	}

	/** Parses a fixed 10–99 CDS with the supplied attributes and plus, minus or unknown strand. */
	private static GffLine feature(String attributes, int strand){
		return new GffLine(("c\tfixture\tCDS\t10\t99\t.\t"+(strand==0 ? "+" : strand==1 ? "-" : "?")
			+"\t0\t"+attributes).getBytes(StandardCharsets.US_ASCII));
	}

	/** Snapshots CDS, tRNA, 16S, 23S, 5S and 18S selection flags in restoration order. */
	private static boolean[] flags(){
		return new boolean[]{ProkObject.callCDS, ProkObject.calltRNA, ProkObject.call16S,
			ProkObject.call23S, ProkObject.call5S, ProkObject.call18S};
	}

	/** Restores or replaces the six global flags in the order returned by flags(). */
	private static void setFlags(boolean[] values){
		ProkObject.callCDS=values[0];
		ProkObject.calltRNA=values[1];
		ProkObject.call16S=values[2];
		ProkObject.call23S=values[3];
		ProkObject.call5S=values[4];
		ProkObject.call18S=values[5];
	}

	/** Writes or replaces an ASCII mutation fixture and returns its path. */
	private static Path write(Path file, String text)throws Exception{
		Files.write(file, text.getBytes(StandardCharsets.US_ASCII));
		return file;
	}

	/** Requires an argument/state exception, optionally checking its diagnostic fragment. */
	private static void reject(Checked action, String message)throws Exception{
		try{action.run();}
		catch(IllegalArgumentException | IllegalStateException expected){
			check(message==null || expected.getMessage().contains(message), "Wrong rejection: "+expected);
			return;
		}
		throw new AssertionError("Invalid support input was accepted");
	}

	/** Invalid-input operation whose checked exceptions may propagate outside the expected rejection types. */
	private interface Checked{
		/** Executes the fixture operation. */
		void run()throws Exception;
	}
	/** Fails immediately, counting only successful assertions. */
	private static void check(boolean valid, String message){
		if(!valid){throw new AssertionError(message);}
		checks++;
	}
	/** Successful assertions and explicit successful verification calls in this sequential process. */
	private static int checks;
}
