package parse;

/** Focused fixture for JVM flags removed from application arguments. */
public final class PreParserJvmFlagTest {

	public static void main(String[] args){
		final boolean oldPrint=PreParser.printExecuting;
		PreParser.printExecuting=false;
		try{
			checkFiltered("-Xss64m");
			checkFiltered("Xss=64m");
			checkFiltered("-Xmx2g");
			checkRetained("path=/tmp/xssdata");
		}finally{PreParser.printExecuting=oldPrint;}
		System.out.println("PASS PreParserJvmFlagTest");
	}

	private static void checkFiltered(String flag){
		final PreParser pp=new PreParser(new String[]{flag,"sample=value"},null,null,false,true,false);
		if(!pp.jflag){throw new AssertionError("JVM flag was not recognized: "+flag);}
		if(pp.args.length!=1 || !"sample=value".equals(pp.args[0])){throw new AssertionError("JVM flag leaked to tool args: "+flag);}
	}

	private static void checkRetained(String arg){
		final PreParser pp=new PreParser(new String[]{arg},null,null,false,true,false);
		if(pp.jflag){throw new AssertionError("Tool argument was mistaken for a JVM flag: "+arg);}
		if(pp.args.length!=1 || !arg.equals(pp.args[0])){throw new AssertionError("Tool argument was removed: "+arg);}
	}
}
