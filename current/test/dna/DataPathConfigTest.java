package test.dna;

import java.util.Arrays;
import java.util.Locale;

import dna.Data;

/** Checks path configuration strings without opening reference or resource files.
 * Run each mode in a fresh JVM. ASCII alias checks use Locale.ROOT; this fixture
 * does not claim locale-independent production normalization or filesystem access.
 * Existing literal-null index behavior and retained custom-prefix temporary paths
 * are controls, not redesigned semantics.
 * @author Sucrose
 */
public final class DataPathConfigTest{

	/** Runs case=controls|typed-null|windows; no work directory or file I/O is needed. */
	public static void main(String[] args){
		require(args.length==1 && args[0].startsWith("case="), "Require one case= argument");
		final String mode=args[0].substring(5);
		require("controls".equals(mode) || "typed-null".equals(mode) || "windows".equals(mode),
				"Unknown path test case: "+mode);
		Locale.setDefault(Locale.ROOT);
		final DataPathConfigTest test=new DataPathConfigTest();
		int failures=0;
		if("typed-null".equals(mode)){
			if(!test.check("literal-null-control", "null", nullPaths())){failures++;}
			if(!test.check("typed-null", null, nullPaths())){failures++;}
		}else if("windows".equals(mode)){
			for(String alias : new String[]{"?Shared.WINDOWS", "?shared.windows", "?SHARED.WINDOWS"}){
				if(!test.check(alias, alias, relativePaths("C:/workspace/tempdir/"))){failures++;}
			}
		}else{
			for(String alias : new String[]{"?local", "?LOCAL", ".", "/.", "./", ".\\", "\\."}){
				if(!test.check(alias, alias, relativePaths(""))){failures++;}
			}
			for(String alias : new String[]{"null", "NuLl"}){
				if(!test.check(alias, alias, nullPaths())){failures++;}
			}
			final String temp=System.getenv("TEMPDIR");
			for(String alias : new String[]{"?unix", "?UNIX"}){
				if(!test.check(alias, alias, relativePaths(temp==null ? "" : temp+"/"))){failures++;}
			}
			for(String prefix : new String[]{"fixture/reference", "fixture/reference/", "fixture\\reference"}){
				if(!test.check(prefix, prefix, customPaths())){failures++;}
			}
			// Direct empty-string input is a custom root; CLI null conversion is outside this API.
			if(!test.check("empty-prefix", "", new String[]{"/", "/ref/", "/ref/genome/",
					"/ref/index/", "/ref/genes/", "/ref/chain/", "/ref/qual/", TEMP, CWD})){failures++;}
		}
		require(failures==0, mode+": "+failures+" path configuration checks failed");
		System.out.println("PASS "+mode);
	}

	/** Captures fixed installation and fallback roots before testing mutable reference settings. */
	private DataPathConfigTest(){installation=installationPaths();}

	/** Seeds the same prior temporary path, then compares every mutable root and fixed root. */
	private boolean check(String label, String input, String[] expected){
		Data.setPath("?local");
		Data.ROOT_TEMPDIR=TEMP;
		RuntimeException failure=null;
		try{Data.setPath(input);}
		catch(RuntimeException problem){failure=problem;}
		final String[] actual=referencePaths();
		final boolean retained=Arrays.equals(installation, installationPaths());
		final boolean ok=failure==null && Arrays.equals(expected, actual) && retained;
		System.out.println("case="+label+" exception="+(failure==null ? "none" : failure.getClass().getSimpleName())+
				" roots="+Arrays.toString(actual)+" installationUnchanged="+retained+" pass="+ok);
		return ok;
	}

	/** Current roots in the fixed order used by the independent literal expectations. */
	private static String[] referencePaths(){
		return new String[]{Data.ROOT_BASE, Data.ROOT_REF, Data.ROOT_GENOME, Data.ROOT_INDEX,
				Data.ROOT_GENE, Data.ROOT_CHAIN, Data.ROOT_QUALITY, Data.ROOT_TEMPDIR, Data.ROOT_CURRENT};
	}

	/** Local, Windows and unix selectors share these reference roots but choose different temp paths. */
	private static String[] relativePaths(String temp){
		return new String[]{"", "ref/", "ref/genome/", "ref/index/", "ref/genes/", "ref/chain/",
				"ref/qual/", temp, CWD};
	}

	/** Existing literal-null behavior deliberately retains the relative index prefix and current directory. */
	private static String[] nullPaths(){
		return new String[]{null, null, null, "ref/index/", null, null, null, null, CWD};
	}

	/** Fixed expected custom-prefix spelling, including slash normalization and retained temporary path. */
	private static String[] customPaths(){
		return new String[]{"fixture/reference/", "fixture/reference/ref/", "fixture/reference/ref/genome/",
				"fixture/reference/ref/index/", "fixture/reference/ref/genes/", "fixture/reference/ref/chain/",
				"fixture/reference/ref/qual/", TEMP, CWD};
	}

	/** Reads strings only; no resource resolution or filesystem existence checks are performed. */
	private static String[] installationPaths(){
		return new String[]{Data.ROOT(), Data.RESOURCES(), Data.NETWORKS(), Data.DORI_RESOURCES,
				Data.DORI_NETWORKS, Data.NERSC_RESOURCES, Data.NERSC_NETWORKS};
	}

	/** Always-on oracle for both -ea and -da runs. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Initial fixed installation and shared-fallback strings, excluded from setPath's mutable roots. */
	private final String[] installation;
	/** Sentinel proving custom prefixes retain, and selectors replace, the prior temporary path. */
	private static final String TEMP="retained-temp/";
	/** Existing ROOT_CURRENT contract is the JVM working-directory property. */
	private static final String CWD=System.getProperty("user.dir");
}
