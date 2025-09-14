package shared;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile1;
import json.JsonObject;
import stream.NullOutputStream;

/**
 * Pre-parses the command line to handle:
 * Leading hyphens
 * Java flags
 * Output stream redirection
 * Help/version information
 * 
 * @author Brian Bushnell
 * @date November 28, 2017
 *
 */
public class PreParser {
	
	/**
	 * Auxiliary constructor.
	 * @param args Command line arguments
	 * @param c Class object of the caller
	 * @param printVersion True if the BBTools version should be printed to the screen
	 */
	public PreParser(final String[] args_, final Class<?> c, boolean printVersion){
		this(args_, System.err, c, printVersion, true, true);
	}
	
	/**
	 * Auxiliary constructor.
	 * @param args Command line arguments
	 * @param defaultPrintStream Print stream that will be used if not overridden
	 * @param c Class object of the caller
	 * @param printVersion True if the BBTools version should be printed to the screen
	 */
	public PreParser(final String[] args_, PrintStream defaultPrintStream, final Class<?> c, boolean printVersion){
		this(args_, defaultPrintStream, c, printVersion, true, true);
	}
	
	/**
	 * Primary constructor.
	 * @param args0 Command line arguments
	 * @param defaultPrintStream Print stream that will be used if not overridden
	 * @param c Class object of the caller
	 * @param printVersion True if the BBTools version should be printed to the screen
	 * @param removeKnown Remove parameters from args that are parsed here
	 * @param autoExit Exit if there is a version or help flag
	 */
	public PreParser(final String[] args0, PrintStream defaultPrintStream, final Class<?> c, boolean printVersion, boolean removeKnown, boolean autoExit){
		original=args0;
		if(Shared.COMMAND_LINE==null){Shared.COMMAND_LINE=(original==null ? null : original.clone());}
		if(Shared.mainClass==null){Shared.mainClass=c;}
		
		Parser.parseHelp(args0, autoExit);
		
		String[] args=Parser.parseConfig(args0);
		this.config=(args!=args0);
		
		final ArrayList<String> unknown=new ArrayList<String>(args.length);
		int removed=0, hyphens=0;

		PrintStream outstream=((defaultPrintStream == null) ? System.err : defaultPrintStream);
		boolean help=false, jflag=false, json=false;
		for(int i=0; i<args.length; i++){
			String s=args[i];
			boolean remove=false;
			if(Parser.isJavaFlag(s)){
				jflag=true;
				remove=true;
			}else if(s==null){
				remove=true;
			}else{
				
//				while(a.charAt(0)=='-' && (a.indexOf('.')<0 || i>1 || !new File(a).exists())){a=a.substring(1);} //Another possibility
				
				if(s.length()>0 && s.charAt(0)=='-'){
					int cnt=1;
					while(cnt<s.length() && s.charAt(cnt)=='-'){cnt++;}
					s=s.substring(cnt);
					hyphens+=cnt;
					args[i]=s;//Note that this mutates the original.  Moral:  Don't use hyphens.
				}
				
//				final String[] split=s.split("=");
				final String[] split=Parse.splitOnFirst(s, '=');
				assert(split.length<3) : "To many '=' signs: "+s; //This can't happen anymore
				assert(split!=null && split.length>0) : "Please do not put spaces around = symbols.\nSyntax is 'arg=value' with no spaces.";
				final String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1] : null;
				if(b==null){
					//do nothing
				}else if("null".equalsIgnoreCase(b)){//Replace "null" with nothing
					b=null;
					args[i]=s=(a+"=");
				}
				
				if(a.equals("outstream")){
					outstream=parseOutstream(b);
					remove=removeKnown;
				}else if(a.equals("json")){
					json=Parse.parseBoolean(b);
				}else if(a.equals("silent")){
					silent=Parse.parseBoolean(b);
				}else if(a.equals("printexecuting")){
					printExecuting=Parse.parseBoolean(b);
					remove=true;
				}else if(a.equals("metadatafile")){
					MetadataWriter.fnameStatic=b;
					remove=true;
				}else if(a.equals("bufferbf") || a.equals("bufferbf1")){//for testing
					ByteFile1.BUFFERED=Parse.parseBoolean(b);
					remove=true;
				}
			}
			
			if(remove){removed++;}
			else{unknown.add(s);}
		}
		this.args=(removed==0 ? args : unknown.toArray(new String[0]));
		this.outstream=outstream;
		this.help=help;
		this.jflag=jflag;
		this.hyphens=hyphens;
		this.json=json;
		
		if(json){
			jsonObject=new JsonObject();
			if(c!=null && printExecuting){jsonObject.add("commandLine", c.getName()+" "+Arrays.toString(original));}
			if(printVersion){jsonObject.add("version", Shared.BBTOOLS_VERSION_STRING);}
			jsonObject.add("memory", Shared.memTotal());
			jsonObject.add("assertions", Shared.EA());
			jsonObject.add("javaVersion", Shared.javaVersion);
			jsonObject.add("JVM_args", Shared.JVM_ARGS());
		}else if(!silent){
			if(c!=null && printExecuting){outstream.println("Executing "+c.getName()+" "+Arrays.toString(original));}
			if(printVersion){outstream.println("Version "+Shared.BBTOOLS_VERSION_STRING );}
			if(c!=null || printVersion){outstream.println();}
		}
	}
	
	public static boolean isAmino(String[] args){
		boolean amino=false;
		for(String arg : args){
			if(arg!=null && arg.startsWith("amino")){
				final String[] split=arg.split("=");
				assert(split.length<3) : "To many '=' signs: "+arg;
				final String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1] : null;
				if(a.equals("amino")){amino=Parse.parseBoolean(b);}
			}
		}
		return amino;
	}

	public static int stripHyphens(String[] args){
		int stripped=0;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			if(Parser.isJavaFlag(arg)){
				//Do nothing
			}else if(arg.length()>0 && arg.charAt(0)=='-'){
				int cnt=1;
				while(cnt<arg.length() && arg.charAt(cnt)=='-'){cnt++;}
				arg=arg.substring(cnt);
				stripped+=cnt;
			}
		}
		return stripped;
	}
	
	private static PrintStream parseOutstream(final String b) {
		
		try {
			if(b==null || b.equalsIgnoreCase("null")){return new PrintStream(new NullOutputStream());}
			else if(b.equalsIgnoreCase("stdout") || b.startsWith("stdout.") || b.equals("System.out") || b.equalsIgnoreCase("sysout")){return System.out;}
			else if(b.equalsIgnoreCase("stderr") || b.startsWith("stderr.") || b.equals("System.err")){return System.err;}
			else{return new PrintStream(b);}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		KillSwitch.kill("Unable to process argument outstream="+b);
		return null; //Unreachable
	}
	
	public final String[] original;
	public final String[] args;
	public final PrintStream outstream;
	public final boolean help;
	public final boolean config;
	public final boolean jflag;
	public final boolean json;
	public final int hyphens;
	public JsonObject jsonObject; 

	public static boolean printExecuting=true;
	public static boolean silent=false;
	
}
