package test.parse;

import java.util.Arrays;

import parse.PreParser;

/** Focused fixture for JVM-option removal and legacy argument normalization. */
public final class PreParserJvmFlagTest{

	/** Checks JVM-option removal and preservation of a similarly named tool argument.
	 * @param args Unused */
	public static void main(String[] args){
		final boolean oldPrint=PreParser.printExecuting;
		PreParser.printExecuting=false;
		try{
			checkFiltered("-Xss64m");
			checkFiltered("Xss=64m");
			checkFiltered("-Xmx2g");
			checkRetained("path=/tmp/xssdata");
			for(int constructor : new int[]{3, 4, 6}){
				checkLegacySpelling(constructor, false);
				checkLegacySpelling(constructor, true);
			}
			checkRetainedKnownOption();
		}finally{PreParser.printExecuting=oldPrint;}
		System.out.println("PASS PreParserJvmFlagTest");
	}

	/** Requires the JVM flag to disappear while an unrelated option remains. */
	private static void checkFiltered(String flag){
		final PreParser pp=new PreParser(new String[]{flag, "sample=value"}, null, null, false, true, false);
		if(!pp.jflag){throw new AssertionError("JVM flag was not recognized: "+flag);}
		if(pp.args.length!=1 || !"sample=value".equals(pp.args[0])){throw new AssertionError("JVM flag leaked to tool args: "+flag);}
	}

	/** Requires a tool argument containing JVM-like text to remain unchanged. */
	private static void checkRetained(String arg){
		final PreParser pp=new PreParser(new String[]{arg}, null, null, false, true, false);
		if(pp.jflag){throw new AssertionError("Tool argument was mistaken for a JVM flag: "+arg);}
		if(pp.args.length!=1 || !arg.equals(pp.args[0])){throw new AssertionError("Tool argument was removed: "+arg);}
	}

	/** Existing constructors must keep their original spelling and array-ownership defaults. */
	private static void checkLegacySpelling(int constructor, boolean filtered){
		final String[] input=(filtered ? new String[]{"-Xss64m", "--outstream=stderr", "-1", "--label", "--keep=null"} :
				new String[]{"-1", "--label", "--keep=null"});
		final PreParser pp=(constructor==3 ? new PreParser(input, null, false) :
				constructor==4 ? new PreParser(input, System.err, null, false) :
				new PreParser(input, System.err, null, false, true, false));
		final String[] expected={"1", "label", "keep="};
		if(!Arrays.equals(pp.args, expected) || pp.original!=input || pp.jflag!=filtered ||
				pp.outstream!=System.err || pp.hyphens!=(filtered ? 7 : 5) || (pp.args==input)==filtered){
			throw new AssertionError("Legacy constructor "+constructor+" filtered="+filtered+
					" changed normalization or ownership: "+Arrays.toString(pp.args));
		}
		System.out.println("constructor="+constructor+" filtered="+filtered+" args="+Arrays.toString(pp.args)+
				" hyphens="+pp.hyphens+" originalRetained=true inputArrayReturned="+(pp.args==input));
	}

	/** The six-argument constructor may process an outstream while retaining its normalized token. */
	private static void checkRetainedKnownOption(){
		final String[] input={"--outstream=stderr", "-1"};
		final PreParser pp=new PreParser(input, System.err, null, false, false, false);
		if(pp.args!=input || pp.outstream!=System.err || pp.hyphens!=3 ||
				!Arrays.equals(pp.args, new String[]{"outstream=stderr", "1"})){
			throw new AssertionError("removeKnown=false changed legacy option retention: "+Arrays.toString(pp.args));
		}
		System.out.println("removeKnown=false args="+Arrays.toString(pp.args)+" hyphens="+pp.hyphens);
	}
}
