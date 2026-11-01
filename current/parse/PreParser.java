package parse;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile1;
import json.JsonObject;
import shared.KillSwitch;
import shared.MetadataWriter;
import shared.Shared;
import stream.NullOutputStream;

/**
 * Performs the first parsing pass over tool command lines.
 * <p>
 * This removes Java flags and selected process-wide options, records whether
 * options such as config files, JSON mode, or Java flags were encountered, and
 * leaves the remaining arguments for tool-specific parsers.
 */
public class PreParser{
	
	/**
	 * Auxiliary constructor with default print stream and auto-exit enabled.
	 * @param args_ Command line arguments
	 * @param c Class object of the caller
	 * @param printVersion True if the BBTools version should be printed
	 */
	public PreParser(final String[] args_, final Class<?> c, boolean printVersion){
		this(args_, System.err, c, printVersion, true, true);
	}
	
	/**
	 * Auxiliary constructor with auto-exit and parameter removal enabled.
	 *
	 * @param args_ Command line arguments
	 * @param defaultPrintStream Print stream that will be used if not overridden
	 * @param c Class object of the caller
	 * @param printVersion True if the BBTools version should be printed
	 */
	public PreParser(final String[] args_, PrintStream defaultPrintStream,
			final Class<?> c, boolean printVersion){
		this(args_, defaultPrintStream, c, printVersion, true, true);
	}
	
	/**
	 * Primary constructor that pre-parses command line arguments.
	 * Handles leading hyphens, Java flags, output stream redirection, and help/version info.
	 * Sets up JSON output mode if requested and initializes the filtered argument array.
	 *
	 * @param args0 Original command line arguments
	 * @param defaultPrintStream Print stream that will be used if not overridden
	 * @param c Class object of the caller for identification
	 * @param printVersion True if the BBTools version should be printed
	 * @param removeKnown Remove parameters from args that are parsed here
	 * @param autoExit Exit if there is a version or help flag
	 */
	public PreParser(final String[] args0, PrintStream defaultPrintStream,
			final Class<?> c, boolean printVersion, boolean removeKnown,
			boolean autoExit){
		original=args0;
		if(Shared.COMMAND_LINE==null){Shared.COMMAND_LINE=(original==null ? null : original.clone());}
		if(Shared.mainClass==null){Shared.mainClass=c;}
		
		Parser.parseHelp(args0, autoExit);
		
		String[] args=Parser.parseConfig(args0);
		this.config=(args!=args0);
		
		final ArrayList<String> unknown=new ArrayList<String>(args.length);
		int removed=0, hyphens=0;

		PrintStream outstream=((defaultPrintStream==null) ? System.err : defaultPrintStream);
		boolean outstreamOwned=false;
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
				
//				while(a.charAt(0)=='-' && (a.indexOf('.')<0 || i>1 ||
//						!new File(a).exists())){a=a.substring(1);} //Another possibility
				
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
				assert(split!=null && split.length>0) :
					"Please do not put spaces around = symbols.\n"+
					"Syntax is 'arg=value' with no spaces.";
				final String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1] : null;
				if(b==null){
					//do nothing
				}else if("null".equalsIgnoreCase(b)){//Replace "null" with nothing
					b=null;
					args[i]=s=(a+"=");
				}
				
				if(a.equals("outstream")){
					try{
						final OutstreamReplacement replacement=parseOutstream(b);
						if(outstreamOwned){outstream.close();}
						outstream=replacement.outstream;
						outstreamOwned=replacement.owned;
					}catch(FileNotFoundException e){
						e.printStackTrace();
						if(outstreamOwned){outstream.close();}
						KillSwitch.kill("Unable to process argument outstream="+b);
					}
					remove=removeKnown;
				}else if(a.equals("json")){
					json=Parse.parseBoolean(b);
				}else if(a.equals("proxyhost")){
					proxyhost=b;
					remove=true;
				}else if(a.equals("proxyport")){
					proxyport=b;
					remove=true;
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
			if(c!=null && printExecuting){
				jsonObject.add("commandLine", c.getName()+" "+Arrays.toString(original));
			}
			if(printVersion){jsonObject.add("version", Shared.BBTOOLS_VERSION_STRING);}
			jsonObject.add("memory", Shared.memTotal());
			jsonObject.add("assertions", Shared.EA());
			jsonObject.add("javaVersion", Shared.javaVersion);
			jsonObject.add("JVM_args", Shared.JVM_ARGS());
		}else if(!silent){
			if(c!=null && printExecuting){
				outstream.println("Executing "+c.getName()+" "+Arrays.toString(original));
			}
			if(printVersion){outstream.println("Version "+Shared.BBTOOLS_VERSION_STRING);}
			if(c!=null || printVersion){outstream.println();}
		}
	}
	
	/**
	 * Reports whether lowercase {@code amino} is explicitly enabled in an argument list.
	 * Arguments with mixed-case names or names that merely start with {@code amino},
	 * such as {@code aminopath}, are ignored.
	 *
	 * @param args Arguments to scan
	 * @return The last exact lowercase {@code amino} flag value, or false if absent
	 */
	public static boolean isAmino(String[] args){
		boolean amino=false;
		for(String arg : args){
			if(arg!=null && arg.startsWith("amino")){
				final String[] split=Parse.splitOnFirst(arg, '=');
				assert(split.length<3) : "To many '=' signs: "+arg;
				final String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1] : null;
				if(a.equals("amino")){amino=Parse.parseBoolean(b);}
			}
		}
		return amino;
	}

	/**
	 * Counts leading hyphens on non-Java arguments.
	 * Java flags are skipped so that leading JVM option hyphens are not counted.
	 *
	 * @param args Arguments to inspect
	 * @return Number of leading hyphen characters counted
	 */
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
	
	/**
	 * Resolves an outstream option to a print stream and its ownership.
	 * Null-like values silence output with constructor-owned streams,
	 * {@code stdout} and {@code stderr} borrow their system streams,
	 * and other values open constructor-owned file streams.
	 *
	 * @param b Stream specifier from the outstream flag
	 * @return PrintStream and ownership for the requested destination
	 * @throws FileNotFoundException If the output path cannot be opened
	 */
	private static OutstreamReplacement parseOutstream(final String b)
			throws FileNotFoundException{
		if(b==null || b.equalsIgnoreCase("null")){
			return new OutstreamReplacement(new PrintStream(new NullOutputStream()), true);
		}else if(b.equalsIgnoreCase("stdout") || b.startsWith("stdout.") ||
				b.equals("System.out") || b.equalsIgnoreCase("sysout")){
			return new OutstreamReplacement(System.out, false);
		}else if(b.equalsIgnoreCase("stderr") || b.startsWith("stderr.") ||
				b.equals("System.err")){
			return new OutstreamReplacement(System.err, false);
		}
		return new OutstreamReplacement(new PrintStream(b), true);
	}

	/** Stream resolved from an outstream option and whether this constructor owns it. */
	private static final class OutstreamReplacement{
		OutstreamReplacement(PrintStream outstream_, boolean owned_){
			outstream=outstream_;
			owned=owned_;
		}

		/** Selected output stream. */
		final PrintStream outstream;
		/** True when the constructor opened and must close this stream on replacement. */
		final boolean owned;
	}
	
	/** Original argument array reference, after any in-place normalization. */
	public final String[] original;
	/** Arguments remaining after recognized pre-parser options are removed. */
	public final String[] args;
	/**
	 * Selected output stream for status and version messages. File and null streams
	 * are constructor-created; default and explicit System streams are borrowed.
	 * The caller must close the final constructor-created stream; borrowed streams
	 * retain their original ownership.
	 */
	public final PrintStream outstream;
	/** Always false after construction; recognized help may exit when autoExit is enabled. */
	public final boolean help;
	/** True when a config file expanded the original arguments. */
	public final boolean config;
	/** True when at least one Java flag was removed. */
	public final boolean jflag;
	/** True when JSON status output is enabled. */
	public final boolean json;
	/** Number of leading hyphens removed from non-Java arguments. */
	public final int hyphens;
	/** Parsed proxy port from {@code proxyport}, or null when unset. */
	public String proxyport=null;
	/** Parsed proxy host from {@code proxyhost}, or null when unset. */
	public String proxyhost=null;
	/** JSON status object; populated only when {@link #json} is true. */
	public JsonObject jsonObject; 

	/** Process-global control for command-line echoing. */
	public static boolean printExecuting=true;
	/** Process-global control for startup status output. */
	public static boolean silent=false;
	
}
