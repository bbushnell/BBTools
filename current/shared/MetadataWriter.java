package shared;

import java.net.InetAddress;
import java.net.UnknownHostException;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;

import fileIO.ReadWrite;
import json.JsonObject;

/**
 * Writes execution metadata for BBTools programs to files in TSV or JSON format.
 * Captures runtime information including timestamps, system details, version info,
 * command line parameters, and processing statistics for program execution tracking.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class MetadataWriter {
	
	/**
	 * Writes execution metadata to a file in the configured format.
	 * Uses the global jsonMode flag to determine output format (JSON or TSV).
	 *
	 * @param fname Output filename, or null to use fnameStatic
	 * @param readsIn Number of input reads processed
	 * @param basesIn Number of input bases processed
	 * @param readsOut Number of output reads produced
	 * @param basesOut Number of output bases produced
	 * @param append Whether to append to existing file or overwrite
	 */
	public static void write(String fname, long readsIn, long basesIn, long readsOut, long basesOut, boolean append){
		if(fname==null){fname=fnameStatic;}
		if(fname==null){return;}
		fnameStatic=null;
		final String s;
		if(jsonMode){
			s=toJson(readsIn, basesIn, readsOut, basesOut);
		}else{
			s=toTsv(readsIn, basesIn, readsOut, basesOut);
		}
		ReadWrite.writeStringInThread(s, fname, append);
	}
	
	/**
	 * Formats execution metadata as tab-separated values.
	 * Includes timestamp, hostname, BBTools version, Java version, command line,
	 * shell script equivalent, and read/base counts.
	 *
	 * @param readsIn Number of input reads processed
	 * @param basesIn Number of input bases processed
	 * @param readsOut Number of output reads produced
	 * @param basesOut Number of output bases produced
	 * @return TSV-formatted metadata string
	 */
	public static String toTsv(long readsIn, long basesIn, long readsOut, long basesOut){
		Map<String,String> env=System.getenv();
		StringBuilder sb=new StringBuilder();
		
		sb.append("Time\t").append(new Date()).append('\n');
		try {
			sb.append("Host\t").append(InetAddress.getLocalHost().getHostName()).append('\n');
		} catch (UnknownHostException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		sb.append("BBToolsVersion\t").append(Shared.BBTOOLS_VERSION_STRING).append('\n');
		sb.append("JavaVersion\t").append(Shared.javaVersion).append('\n');
		sb.append("Command\t").append(Shared.fullCommandline()).append('\n');
		String script=commandToShellscript(Shared.fullCommandline());
		if(script!=null){sb.append("Script\t").append(script).append('\n');}
		sb.append("ReadsIn\t").append(readsIn).append('\n');
		sb.append("BasesIn\t").append(basesIn).append('\n');
		sb.append("ReadsOut\t").append(readsOut).append('\n');
		sb.append("BasesOut\t").append(basesOut).append('\n');
		return sb.toString();
	}
	
	/**
	 * Formats execution metadata as a JSON object.
	 * Contains the same metadata as TSV format but structured as JSON.
	 *
	 * @param readsIn Number of input reads processed
	 * @param basesIn Number of input bases processed
	 * @param readsOut Number of output reads produced
	 * @param basesOut Number of output bases produced
	 * @return JSON-formatted metadata string
	 */
	public static String toJson(long readsIn, long basesIn, long readsOut, long basesOut){
		Map<String,String> env=System.getenv();
		JsonObject jo=new JsonObject();
		
		jo.add("Time", new Date().toString());
		try {
			jo.add("Host", InetAddress.getLocalHost().getHostName());
		} catch (UnknownHostException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
//		jo.add("Host", env.get("HOSTNAME"));
//		assert(false) : env;InetAddress.getLocalHost().getHostName()
		jo.add("BBToolsVersion", Shared.BBTOOLS_VERSION_STRING);
		jo.add("JavaVersion", Shared.javaVersion);
		jo.add("Command", Shared.fullCommandline());
		String script=commandToShellscript(Shared.fullCommandline());
		if(script!=null){jo.add("Script", script);}
		jo.add("ReadsIn", readsIn);
		jo.add("BasesIn", basesIn);
		jo.add("ReadsOut", readsOut);
		jo.add("BasesOut", basesOut);
		return jo.toString();
	}
	
	/**
	 * Translates Java class command line to equivalent shell script command.
	 * Uses the shellMap to convert Java class names to their corresponding
	 * shell script names (e.g., "jgi.BBDuk" becomes "bbduk.sh").
	 *
	 * @param command Full command line string starting with Java class name
	 * @return Shell script equivalent command, or null if no mapping found
	 */
	private static String commandToShellscript(String command) {
		if(command==null) {return null;}
		final HashMap<String, String> map=shellMap();
		String[] split=command.split(" ");
		int pos=0;
		for(pos=0; pos<split.length; pos++){
			if(map.containsKey(split[pos])){
				split[pos]=map.get(split[pos]);
				break;
			}
		}
		if(pos>=split.length){return null;}
		StringBuilder sb=new StringBuilder();
		for(; pos<split.length; pos++){
			sb.append(split[pos]);
			sb.append(' ');
		}
		sb.setLength(sb.length()-1);
		return sb.toString();
	}
	
	/** Anything in this map will have the command line translated
	 * to the equivalent shell script command for the "Script" key */
	private static HashMap<String, String> shellMap(){
		HashMap<String, String> map=new HashMap<String, String>();
		map.put("bloom.BloomFilterCorrectorWrapper", "bbcms.sh");
		map.put("jgi.ReformatReads", "reformat.sh");
		map.put("jgi.FungalRelease", "fungalrelease.sh");
		map.put("jgi.MakeLengthHistogram", "readlength.sh");
		map.put("jgi.AssemblyStats2", "stats.sh");
		map.put("jgi.BBDuk", "bbduk.sh");
		map.put("assemble.Tadpole", "tadpole.sh");
		map.put("sketch.SendSketch", "sendsketch.sh");
		map.put("clump.Clumpify", "clumpify.sh");
		return map;
	}
	
	/** Static filename to use when no filename is explicitly provided to write() */
	public static String fnameStatic;
	/** Whether to output metadata in JSON format (true) or TSV format (false) */
	public static boolean jsonMode=true;
	
}
