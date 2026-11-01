package test.shared;

import java.io.ByteArrayOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.zip.GZIPOutputStream;

import fileIO.ByteFile;
import fileIO.ReadWrite;
import shared.Tools;

/** Literal collection-count and input-error checks for Tools.addNames.
 * Generated fixtures remain in the supplied directory for inspection.
 * @author Shinobu
 */
public final class ToolsAddNamesTest{

	/** Runs all independent cases; accepts exactly one work=directory argument. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1 || !args[0].startsWith("work=")){
			throw new IllegalArgumentException("Expected work=directory");
		}
		final Path work=Paths.get(args[0].substring(5));
		Files.createDirectories(work);
		ByteFile.FORCE_MODE_BF1=true;
		ByteFile.FORCE_MODE_BF2=ByteFile.FORCE_MODE_BF3=ByteFile.FORCE_MODE_BF4=false;
		ReadWrite.USE_UNPIGZ=false;
		ReadWrite.ALLOW_NATIVE_BGZF=false;
		final ToolsAddNamesTest test=new ToolsAddNamesTest();
		test.check("literal insertion and duplicate", () -> {
			final String literal=work.resolve("absent-name").toString();
			require(!Files.exists(Paths.get(literal)), "Literal control must not name a real file");
			final HashSet<String> names=new HashSet<String>();
			require(Tools.addNames(literal, names, false)==1 && names.contains(literal), "First literal insertion");
			require(Tools.addNames(literal, names, false)==0 && names.size()==1, "Set duplicate adds no entry");
		});
		final String[] extensions={"txt", "fa", "fq", "sam"};
		final String[] contents={"alpha\nbeta\nalpha\n",
				">alpha\nA\n>beta\nC\n>alpha\nG\n",
				"@alpha\nA\n+\nI\n@beta\nC\n+\nI\n@alpha\nG\n+\nI\n",
				"@HD\tVN:1.6\nalpha\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\n"+
				"beta\t4\t*\t0\t0\t*\t*\t0\t0\tC\tI\nalpha\t4\t*\t0\t0\t*\t*\t0\t0\tG\tI\n"};
		for(int i=0; i<extensions.length; i++){
			final Path path=work.resolve("names."+extensions[i]);
			Files.write(path, contents[i].getBytes(StandardCharsets.US_ASCII));
			test.check(extensions[i]+" list contents and count", () -> {
				final ArrayList<String> list=new ArrayList<String>();
				final int listed=Tools.addNames(path.toString(), list, false);
				require(list.equals(Arrays.asList("alpha", "beta", "alpha")), "Names or file order changed: "+list);
				require(listed==3, "Expected three list insertions, actual="+listed);
			});
			test.check(extensions[i]+" set contents and count", () -> {
				final HashSet<String> set=new HashSet<String>();
				set.add("alpha");
				final int added=Tools.addNames(path.toString(), set, false);
				require(added==1 && set.equals(new HashSet<String>(Arrays.asList("alpha", "beta"))),
						"Expected only beta to be newly inserted, actual count="+added+", names="+set);
			});
		}
		final Path empty=work.resolve("empty.txt");
		Files.write(empty, new byte[0]);
		test.check("empty file control", () -> {
			final ArrayList<String> names=new ArrayList<String>();
			require(Tools.addNames(empty.toString(), names, false)==0 && names.isEmpty(), "Empty file must add nothing");
		});
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		try(GZIPOutputStream gzip=new GZIPOutputStream(bytes)){
			gzip.write("alpha\nbeta\n".getBytes(StandardCharsets.US_ASCII));
		}
		final byte[] corrupt=bytes.toByteArray();
		final Path valid=work.resolve("valid.txt.gz");
		Files.write(valid, corrupt);
		test.check("gzip input control", () -> {
			final ArrayList<String> names=new ArrayList<String>();
			final int added=Tools.addNames(valid.toString(), names, false);
			require(names.equals(Arrays.asList("alpha", "beta")), "Valid gzip names changed: "+names);
			require(added==2, "Expected two gzip insertions, actual="+added);
		});
		corrupt[corrupt.length-8]^=1;// Change only the CRC; preserve the compressed payload.
		final Path bad=work.resolve("corrupt.txt.gz");
		Files.write(bad, corrupt);
		test.check("corrupt input rejection", () -> {
			boolean rejected=false;
			try{Tools.addNames(bad.toString(), new ArrayList<String>(), false);}
			catch(RuntimeException expected){rejected=true;}
			require(rejected, "Corrupt names must not be accepted as a complete collection");
		});
		System.out.println("ToolsAddNamesTest cases="+test.cases+" failures="+test.failures);
		if(test.failures>0){throw new AssertionError("Name collection failures="+test.failures);}
	}

	/** Reports every finite case rather than stopping at the first baseline failure. */
	private void check(String name, Action action){
		cases++;
		try{action.run();System.out.println(name+" PASS");}
		catch(Exception | AssertionError failure){failures++;System.out.println(name+" FAIL: "+failure);}
	}

	/** Keeps expectations active with assertions disabled. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** A finite case that may read or write a tiny fixture. */
	private interface Action{void run() throws Exception;}
	/** Attempted and failed groups. */
	private int cases, failures;
}
