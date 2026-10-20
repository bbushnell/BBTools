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

/** Finite legacy clumpify.sh conservation checks; launches one isolated JVM per case.
 * Run on a compute host because the normal tools create I/O and worker threads.
 * @author Shinobu, Sucrose
 */
public class ClumpCliTest{

	/** Takes the runtime root, a fresh fixture directory, and optional comma-separated case names. */
	public static void main(String[] args) throws Exception{
		assert(args.length==2 || args.length==3) : "Expected runtime, fresh directory, and optional case names";
		selected=args.length<3 ? null : new java.util.HashSet<String>(Arrays.asList(args[2].split(",")));
		final Path runtime=Paths.get(args[0]), root=Paths.get(args[1]);
		Files.createDirectory(root);
		final Path single=root.resolve("single.fq"), left=root.resolve("left.fq"), right=root.resolve("right.fq");
		final Path empty=root.resolve("empty.fq");
		write(single, SINGLE); write(left, LEFT); write(right, RIGHT); write(empty, "");
		System.out.println("case\treads\tpayloads\tinput_unchanged");
		run(runtime, root, "single-g1", single, null, SINGLE, false, "groups=1");
		run(runtime, root, "single-g3", single, null, SINGLE, false, "groups=3");
		run(runtime, root, "single-g3-p2", single, null, SINGLE, false, "groups=3", "passes=2");
		run(runtime, root, "single-v2-p2", single, null, SINGLE, false, "groups=3", "passes=2", "v2=t");
		run(runtime, root, "single-v1-p2", single, null, SINGLE, false, "groups=3", "passes=2", "v1=t");
		run(runtime, root, "paired-g1", left, right, LEFT+RIGHT, true, "groups=1");
		run(runtime, root, "paired-g3-p2", left, right, LEFT+RIGHT, true, "groups=3", "passes=2");
		run(runtime, root, "paired-namesort", left, right, LEFT+RIGHT, true, "groups=3", "namesort=t");
		run(runtime, root, "paired-namesort-v1", left, right, LEFT+RIGHT, true, "groups=3", "namesort=t", "v1=t");
		run(runtime, root, "unpair-g3", left, right, LEFT+RIGHT, false, "groups=3", "unpair=t");
		run(runtime, root, "empty-g1", empty, null, "", false, "groups=1");
		run(runtime, root, "empty-g3-p2", empty, null, "", false, "groups=3", "passes=2");
		run(runtime, root, "shortnames-g3", single, null, SHORT_NAMES, false, "groups=3", "shortname=t");
		assert(selected==null || selected.isEmpty()) : "Unknown CLI case names: "+selected;
	}

	/** Compares exact named base/quality records; order is deliberately unrestricted. */
	private static void run(Path runtime, Path root, String name, Path in1, Path in2,
			String expectedText, boolean pairedOutput, String... options) throws Exception{
		if(selected!=null && !selected.remove(name)){return;}
		final Path dir=root.resolve(name);
		Files.createDirectory(dir);
		final byte[] input1=Files.readAllBytes(in1), input2=in2==null ? null : Files.readAllBytes(in2);
		final String output=dir.resolve(pairedOutput ? "out#.fq" : "out.fq").toString();
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(runtime.resolve("clumpify.sh").toString());
		command.addAll(Arrays.asList("-Xmx512m", "-Xms64m", "t=1", "k=7", "seed=1", "hashes=4",
				"border=0", "qin=33", "qout=33", "rcomp=f", "ecco=f", "correct=f", "dedupe=f",
				"condense=f", "quantize=f", "reorder=f", "overwrite=f", "deleteinput=f", "in="+in1,
				"out="+output));
		if(in2!=null){command.add("in2="+in2);}else{command.add("interleaved=f");}
		command.addAll(Arrays.asList(options));
		Files.write(dir.resolve("arguments.txt"), command, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final ProcessBuilder builder=new ProcessBuilder(command);
		builder.environment().put("PATH", System.getProperty("java.home")+"/bin:"+builder.environment().get("PATH"));
		builder.redirectOutput(dir.resolve("stdout.txt").toFile());
		builder.redirectError(dir.resolve("stderr.txt").toFile());
		System.err.println("Running "+name+"; logs in "+dir);
		final int status=builder.start().waitFor();
		assert(status==0) : name+" failed with exit "+status+"; inspect "+dir.resolve("stderr.txt");
		final Map<String, String> expected=parse(Arrays.asList(expectedText.split("\n")));
		final Map<String, String> observed=new LinkedHashMap<String, String>();
		if(pairedOutput){
			final Map<String, String> one=parse(Files.readAllLines(dir.resolve("out1.fq"), StandardCharsets.US_ASCII));
			final Map<String, String> two=parse(Files.readAllLines(dir.resolve("out2.fq"), StandardCharsets.US_ASCII));
			assert(one.equals(parse(Files.readAllLines(in1, StandardCharsets.US_ASCII))))
					: "Read-one membership changed: "+name;
			assert(two.equals(parse(Files.readAllLines(in2, StandardCharsets.US_ASCII))))
					: "Read-two membership changed: "+name;
			final java.util.Iterator<String> mates=two.keySet().iterator();
			for(String id : one.keySet()){
				final String mate=mates.next();
				assert(id.substring(0, id.length()-2).equals(mate.substring(0, mate.length()-2)))
						: "Paired output files disagree on fragment order: "+name;
			}
			if(name.startsWith("paired-namesort")){
				assert(one.keySet().iterator().next().equals("alpha/1")) : "Requested name order was not restored";
			}
			observed.putAll(one); observed.putAll(two);
		}else{
			observed.putAll(parse(Files.readAllLines(dir.resolve("out.fq"), StandardCharsets.US_ASCII)));
		}
		assert(observed.equals(expected)) : "Named sequence/quality records changed: "+name+"; "+observed;
		assert(Arrays.equals(input1, Files.readAllBytes(in1))) : "Input one changed: "+name;
		assert(in2==null || Arrays.equals(input2, Files.readAllBytes(in2))) : "Input two changed: "+name;
		System.out.println(name+'\t'+observed.size()+"\tEXACT\ttrue");
	}

	/** Strict four-line fixture FASTQ reader that rejects repeated identifiers. */
	private static Map<String, String> parse(List<String> lines){
		final Map<String, String> records=new LinkedHashMap<String, String>();
		if(lines.size()==1 && lines.get(0).isEmpty()){return records;}
		assert(lines.size()%4==0) : "Incomplete FASTQ fixture: "+lines;
		for(int i=0; i<lines.size(); i+=4){
			assert(lines.get(i).startsWith("@") && lines.get(i+2).equals("+")) : "Malformed FASTQ fixture";
			final String bases=lines.get(i+1), qualities=lines.get(i+3);
			assert(bases.length()==qualities.length()) : "Sequence/quality length mismatch";
			final String prior=records.put(lines.get(i).substring(1), bases+'\n'+qualities);
			assert(prior==null) : "Duplicate output identifier: "+lines.get(i);
		}
		return records;
	}

	private static void write(Path path, String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	private static java.util.Set<String> selected;
	private static final String SINGLE="@r0\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"+
			"@r1\nTTGCAACCGGTT\n+\nHHHHHHHHHHHH\n"+
			"@r2\nACG\n+\nFFF\n"+
			"@r3\nNNNNNNNN\n+\nGGGGGGGG\n";
	private static final String SHORT_NAMES=SINGLE.replace("@r0", "@0000000000 1:")
			.replace("@r1", "@0000000001 1:").replace("@r2", "@0000000002 1:").replace("@r3", "@0000000003 1:");
	private static final String LEFT="@zeta/1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"+
			"@alpha/1\nTTGCAACCGGTT\n+\nHHHHHHHHHHHH\n";
	private static final String RIGHT="@zeta/2\nTTTTCCCCAAAA\n+\nGGGGGGGGGGGG\n"+
			"@alpha/2\nAAAACCCCTTTT\n+\nFFFFFFFFFFFF\n";
}
