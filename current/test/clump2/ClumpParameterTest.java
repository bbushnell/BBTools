package test.clump2;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

import clump2.Clumpify;
import clump2.PivotSet;
import shared.Shared;
import shared.Timer;

/** Finite constructor/pass-count regressions; creates only disposable fixture inputs.
 * Invalid pass counts are processed only when the constructor accepts them, so the
 * baseline demonstrates the skipped output and optional deletion without workers.
 * @author Shinobu
 */
public class ClumpParameterTest {

	/** Runs baseline or repaired expectations in a previously nonexistent directory. */
	public static void main(String[] args) throws Exception{
		assert(args.length==2) : "Expected baseline|fixed and a fresh fixture directory";
		final boolean fixed=args[0].equals("fixed");
		assert(fixed || args[0].equals("baseline")) : "Unknown expectation: "+args[0];
		Shared.setThreads(1);
		final Path root=Paths.get(args[1]);
		Files.createDirectory(root);
		System.out.println("case\taccepted\tinput_present\toutput_present");
		for(int passes : new int[] {0, -1}){
			for(boolean delete : new boolean[] {false, true}){
				checkPasses(root, passes, delete, fixed);
			}
		}
		for(int groups : new int[] {1, 2}){
			for(int passes : new int[] {1, 2}){
				final Path input=writeInput(root, "valid-g"+groups+"-p"+passes);
				new Clumpify(new String[] {"in="+input, "groups="+groups, "passes="+passes});
				assert(Files.exists(input)) : "Constructor-only control must retain its input";
				System.out.println("valid-g"+groups+"-p"+passes+"\ttrue\ttrue\tfalse");
			}
		}
		checkAlias(root, "mincount", false);
		checkAlias(root, "mincr", !fixed);
	}

	/** Checks skipped processing and deletion solely on the input created here. */
	private static void checkPasses(Path root, int passes, boolean delete, boolean fixed) throws Exception{
		final String name="passes"+passes+"-delete"+delete;
		final Path input=writeInput(root, name);
		final Path output=root.resolve(name+".out.fq");
		Clumpify tool=null;
		boolean rejected=false;
		try{
			tool=new Clumpify(new String[] {"in="+input, "out="+output,
					"groups=2", "passes="+passes, "deleteinput="+delete});
		}catch(IllegalArgumentException e){
			assert(e.getMessage().contains("passes")) : "Unexpected parameter rejection: "+e;
			rejected=true;
		}
		assert(rejected==fixed) : "Invalid passes must be rejected before optional input deletion: "+name;
		if(tool!=null){tool.process(new Timer());}
		final boolean present=Files.exists(input), produced=Files.exists(output);
		assert(present==(fixed || !delete)) : "Unexpected disposable input state for "+name;
		assert(!produced) : "A rejected or zero-iteration pass must not create output: "+output;
		if(present){
			assert(new String(Files.readAllBytes(input), StandardCharsets.US_ASCII).equals(FASTQ))
					: "Input changed without a processing pass: "+input;
		}
		System.out.println(name+"\t"+!rejected+"\t"+present+"\t"+produced);
	}

	/** The supported long option is a control for its missing short alias. */
	private static void checkAlias(Path root, String option, boolean expectUnknown) throws Exception{
		final Path input=writeInput(root, option);
		boolean unknown=false;
		try{
			final PivotSet tool=new PivotSet(new String[] {"in="+input, option+"=3", "k=5"});
			final Field field=PivotSet.class.getDeclaredField("minCount");
			field.setAccessible(true);
			assert(field.getInt(tool)==3) : option+" must select the same minimum count as mincount";
		}catch(AssertionError e){
			if(e.getMessage()==null || !e.getMessage().contains("Unknown parameter "+option+"=3")){throw e;}
			unknown=true;
		}
		assert(unknown==expectUnknown) : "Unexpected alias recognition: "+option;
		System.out.println(option+"\t"+!unknown+"\t"+Files.exists(input)+"\tfalse");
	}

	/** Uses CREATE_NEW to ensure the test never replaces an existing input. */
	private static Path writeInput(Path root, String name) throws Exception{
		final Path path=root.resolve(name+".fq");
		Files.write(path, FASTQ.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		return path;
	}

	private static final String FASTQ="@fixture\nACGTACGT\n+\nIIIIIIII\n";
}
