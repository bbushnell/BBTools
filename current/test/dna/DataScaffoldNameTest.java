package test.dna;

import java.nio.charset.StandardCharsets;
import java.util.HashMap;

import dna.Data;
import dna.ScafLoc;

/** Finite sequential scaffold-name lookup checks with literal names and coordinates.
 * Run each mode in a fresh JVM. Installs only tiny private arrays through public
 * fields and restores the exact previous state. No reference loading, broad reset,
 * concurrent rename or previously borrowed-map update contract is tested.
 * @author Sucrose
 */
public final class DataScaffoldNameTest{

	/** Runs case=controls|cached-rename with assertions enabled or disabled. */
	public static void main(String[] args){
		require(args.length==1 && args[0].startsWith("case="), "Require one case= argument");
		final String mode=args[0].substring(5);
		final int oldBuild=Data.GENOME_BUILD, oldCount=Data.numScaffolds;
		final byte[][][] oldNames=Data.scaffoldNames;
		final int[][] oldLocs=Data.scaffoldLocs;
		final HashMap<String, ScafLoc> oldTable=Data.scaffoldNameTable;
		try{
			if("controls".equals(mode)){controls();}
			else if("cached-rename".equals(mode)){cachedRename();}
			else{throw new IllegalArgumentException(mode);}
		}finally{
			Data.GENOME_BUILD=oldBuild;
			Data.numScaffolds=oldCount;
			Data.scaffoldNames=oldNames;
			Data.scaffoldLocs=oldLocs;
			Data.scaffoldNameTable=oldTable;
		}
		require(Data.GENOME_BUILD==oldBuild && Data.numScaffolds==oldCount &&
				Data.scaffoldNames==oldNames && Data.scaffoldLocs==oldLocs && Data.scaffoldNameTable==oldTable,
				"Fixture must restore the exact borrowed public state");
		System.out.println("PASS "+mode+" stateRestored=true");
	}

	/** Checks trimming before lookup, repeated trimming, unchanged cached names and null metadata. */
	private static void controls(){
		install("alpha description");
		final byte[] original=Data.scaffoldNames[1][0];
		Data.trimScaffoldNames();
		boolean pass=trimmedLookup("trim-before-lookup");
		require("alpha description".equals(new String(original, StandardCharsets.US_ASCII)),
				"Trimming must replace the name entry without changing borrowed original bytes");
		Data.trimScaffoldNames();
		pass&=trimmedLookup("repeated-trim");

		install("alpha");
		pass&=trimmedLookup("unchanged-before");
		Data.trimScaffoldNames();
		pass&=trimmedLookup("unchanged-after");

		Data.GENOME_BUILD=-1;
		Data.numScaffolds=0;
		Data.scaffoldNames=null;
		Data.scaffoldLocs=null;
		Data.scaffoldNameTable=null;
		Data.trimScaffoldNames();
		final boolean nullState=Data.scaffoldNames==null && Data.scaffoldLocs==null &&
				Data.scaffoldNameTable==null && Data.getScafLoc("alpha")==null && Data.getScafLoc(ascii("alpha"))==null;
		System.out.println("null-metadata pass="+nullState);
		pass&=nullState;

		Data.GENOME_BUILD=7;
		Data.scaffoldNames=new byte[][][]{null, {null}};
		Data.scaffoldLocs=new int[][]{null, {0}};
		Data.trimScaffoldNames();
		final boolean nullEntries=Data.scaffoldNames[0]==null && Data.scaffoldNames[1][0]==null &&
				Data.getScafLoc("alpha")==null && Data.getScafLoc(ascii("alpha"))==null && Data.scafNameTable().isEmpty();
		System.out.println("null-row-and-name pass="+nullEntries);
		pass&=nullEntries;
		require(pass, "Scaffold-name control failed");
	}

	/** Builds the old-name cache before trimming; public lookups must then describe current names. */
	private static void cachedRename(){
		install("alpha description");
		require(matches(Data.getScafLoc("alpha description"), "alpha description", 10) &&
				matches(Data.getScafLoc(ascii("alpha description")), "alpha description", 10),
				"Precondition: both overloads must find the original name at chromosome1/location10");
		require(matches(Data.getScafLoc("plain"), "plain", 30), "Precondition: unchanged scaffold is at location30");
		Data.trimScaffoldNames();
		boolean pass=trimmedLookup("lookup-before-trim");
		Data.trimScaffoldNames();
		pass&=trimmedLookup("lookup-before-repeated-trim");
		require(pass, "Name cache retained stale keys after trimming scaffold names");
	}

	/** Installs two distinct names on chromosome1 with fixed zero-based starts10 and30. */
	private static void install(String first){
		Data.GENOME_BUILD=7;
		Data.numScaffolds=2;
		Data.scaffoldNames=new byte[][][]{null, {ascii(first), ascii("plain")}};
		Data.scaffoldLocs=new int[][]{null, {10, 30}};
		Data.scaffoldNameTable=null;
	}

	/** Compares both lookup overloads to independent literal names/coordinates and absent old keys. */
	private static boolean trimmedLookup(String label){
		final ScafLoc alphaString=Data.getScafLoc("alpha"), alphaBytes=Data.getScafLoc(ascii("alpha"));
		final ScafLoc plainString=Data.getScafLoc("plain"), plainBytes=Data.getScafLoc(ascii("plain"));
		final ScafLoc oldString=Data.getScafLoc("alpha description"), oldBytes=Data.getScafLoc(ascii("alpha description"));
		final boolean names="alpha".equals(new String(Data.scaffoldNames[1][0], StandardCharsets.US_ASCII)) &&
				"plain".equals(new String(Data.scaffoldNames[1][1], StandardCharsets.US_ASCII));
		final boolean pass=names && matches(alphaString, "alpha", 10) && matches(alphaBytes, "alpha", 10) &&
				matches(plainString, "plain", 30) && matches(plainBytes, "plain", 30) && oldString==null && oldBytes==null;
		System.out.println(label+" namesTrimmed="+names+" alphaString="+describe(alphaString)+
				" alphaBytes="+describe(alphaBytes)+" plainString="+describe(plainString)+
				" plainBytes="+describe(plainBytes)+" oldString="+describe(oldString)+
				" oldBytes="+describe(oldBytes)+" pass="+pass);
		return pass;
	}

	/** No map identity is required: already borrowed maps are outside this sequential lookup oracle. */
	private static boolean matches(ScafLoc loc, String name, int position){
		return loc!=null && name.equals(loc.name) && loc.chrom==1 && loc.loc==position;
	}
	private static String describe(ScafLoc loc){return loc==null ? "null" : loc.name+":"+loc.chrom+":"+loc.loc;}
	private static byte[] ascii(String text){return text.getBytes(StandardCharsets.US_ASCII);}

	/** Keeps the oracle active with assertions disabled. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
