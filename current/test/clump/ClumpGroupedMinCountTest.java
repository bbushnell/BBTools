package test.clump;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** End-to-end grouped mincount checks through clumpify.sh.
 * Baseline mode records the V2/V3 grouped multi-pass lost-table assertion;
 * fixed mode requires every case to conserve all literal FASTQ records.
 * @author Jean
 */
public final class ClumpGroupedMinCountTest{

	/** Takes baseline/fixed, runtime directory, and a fresh fixture directory. */
	public static void main(String[] args) throws Exception{
		assert(args.length==3) : "Expected phase, runtime, and fresh fixture directory";
		final boolean fixed=args[0].equals("fixed");
		assert(fixed || args[0].equals("baseline")) : "Unknown phase: "+args[0];
		final Path runtime=Paths.get(args[1]), root=Paths.get(args[2]);
		Files.createDirectory(root);
		System.out.println("case\texit\tmissing_table_assertion\texact_payload\tinput_retained");

		run(runtime, root, "v2-one-pass", fixed, false, "v2=t", "passes=1");
		run(runtime, root, "v3-one-pass", fixed, false, "v3=t", "passes=1");
		run(runtime, root, "v1-multipass", fixed, false, "v1=t", "passes=2");
		run(runtime, root, "v2-multipass", fixed, true, "v2=t", "passes=2");
		run(runtime, root, "v3-multipass", fixed, true, "v3=t", "passes=2");
	}

	/** Saves one clumpify.sh invocation and checks the legacy assertion or full
	 * payload conservation. */
	private static void run(Path runtime, Path root, String name, boolean fixed,
			boolean legacyFails, String version, String passCount) throws Exception{
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final Path in=dir.resolve("input.fq"), out=dir.resolve("out.fq");
		final byte[] original=INPUT.getBytes(StandardCharsets.US_ASCII);
		Files.write(in, original, StandardOpenOption.CREATE_NEW);

		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(runtime.resolve("clumpify.sh").toString());
		command.addAll(Arrays.asList("-Xmx512m", "-Xms64m", "t=1", "k=7", "mincount=2",
				"seed=1", "hashes=4", "border=0", "qin=33", "qout=33", "rcomp=f",
				"ecco=f", "correct=f", "dedupe=f", "condense=f", "quantize=f",
				"reorder=f", "overwrite=f", "groups=3", "unpair=f", "repair=f",
				"namesort=f", "interleaved=f", "compresstemp=f", "deleteinput=f",
				version, passCount, "in="+in, "out="+out));
		Files.write(dir.resolve("arguments.txt"), command, StandardCharsets.UTF_8,
				StandardOpenOption.CREATE_NEW);

		final ProcessBuilder builder=new ProcessBuilder(command);
		builder.environment().put("PATH",
				System.getProperty("java.home")+"/bin:"+builder.environment().get("PATH"));
		builder.redirectOutput(dir.resolve("stdout.txt").toFile());
		builder.redirectError(dir.resolve("stderr.txt").toFile());
		System.err.println("Running "+name+"; logs in "+dir);
		final int status=builder.start().waitFor();

		final String error=new String(Files.readAllBytes(dir.resolve("stderr.txt")),
				StandardCharsets.UTF_8);
		final boolean missingTable=error.contains("AssertionError")
				&& error.contains("Grouped minCount requires a KmerSplit-produced table");
		final boolean exact=(status==0 && Files.exists(out)
				&& parse(Files.readAllLines(out, StandardCharsets.US_ASCII)).equals(EXPECTED));
		final boolean retained=Files.exists(in) && Arrays.equals(original, Files.readAllBytes(in));

		System.out.println(name+'\t'+status+'\t'+missingTable+'\t'+exact+'\t'+retained);
		if(legacyFails && !fixed){
			assert(status!=0 && missingTable)
					: "Legacy grouped multipass mincount should fail when pass 2 has no "
					+"KmerSplit-built ClumpTools table; "+name+" exit="+status+", stderr="+error;
		}else{
			assert(status==0 && exact) : "Grouped mincount payload changed: "+name;
		}
		assert(retained) : "Original input changed or was deleted: "+name;
	}

	/** Strict four-line fixture FASTQ reader that preserves literal bases and qualities. */
	private static Map<String, String> parse(List<String> lines){
		final Map<String, String> records=new LinkedHashMap<String, String>();
		assert(lines.size()%4==0) : "Incomplete FASTQ fixture; 4-line parsing requires all lines";
		for(int i=0; i<lines.size(); i+=4){
			assert(lines.get(i).startsWith("@") && lines.get(i+2).equals("+"))
					: "Malformed FASTQ fixture at record "+(i/4);
			final String bases=lines.get(i+1), qualities=lines.get(i+3);
			assert(bases.length()==qualities.length())
					: "FASTQ sequence/quality lengths must match for "+lines.get(i);
			final String prior=records.put(lines.get(i).substring(1), bases+'\n'+qualities);
			assert(prior==null) : "Duplicate output identifier would hide a lost record: "+lines.get(i);
		}
		return records;
	}

	private static final String INPUT="@alpha\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"+
			"@beta\nACGTACGTACGT\n+\nHHHHHHHHHHHH\n"+
			"@gamma\nTGCATGCATGCA\n+\nGGGGGGGGGGGG\n"+
			"@delta\nTGCATGCATGCA\n+\nFFFFFFFFFFFF\n";
	private static final Map<String, String> EXPECTED=parse(Arrays.asList(INPUT.split("\n")));
}
