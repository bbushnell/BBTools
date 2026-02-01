package align2;

import java.io.PrintStream;
import java.util.ArrayList;

import dna.Data;
import fileIO.TextFile;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import stream.Read;

/**
 * Wrapper class for executing multiple BBMap variants with batch processing capabilities.
 * Parses input/output file lists and coordinates sequential execution of mapping jobs.
 * Supports BBMap, BBMapPacBio, BBMapPacBioSkimmer, BBMap5, BBMapAcc, and BBSplitter.
 *
 * @author Brian Bushnell
 * @date Mar 27, 2014
 */
public class BBWrap {
	
	public static void main(String[] args){
		BBWrap wrapper=new BBWrap();
		ArrayList<String> list=wrapper.parse(args);
		wrapper.execute(list);
		
		//Close the print stream if it was redirected
		Shared.closeStream(outstream);
	}

	/**
	 * Parses command-line arguments and populates file lists for batch processing.
	 * Extracts file paths, mapper selection, and other parameters while preserving
	 * unrecognized arguments for passing to the underlying mapper.
	 *
	 * @param args Command-line arguments to parse
	 * @return ArrayList of unprocessed arguments to pass to the mapper
	 */
	private final ArrayList<String> parse(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Read.TO_UPPER_CASE=true;
		
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			final String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("path") || a.equals("root")){
				Data.setPath(b);
				args[i]=null;
			}else if(a.equals("mapper")){
				mapper=b;
				args[i]=null;
			}else if(a.equals("ref") || a.equals("reference") || a.equals("fasta")){
				ref=b;
				args[i]=null;
			}else if(a.equals("in") || a.equals("in1")){
				add(b, in1List);
				args[i]=null;
			}else if(a.equals("in2")){
				add(b, in2List);
				args[i]=null;
			}else if(a.equals("inlist") || a.equals("in1list")){
				addFileContents(b, in1List);
				args[i]=null;
			}else if(a.equals("in2list")){
				addFileContents(b, in2List);
				args[i]=null;
			}else if(a.equals("out") || a.equals("out1")){
				add(b, out1List);
				args[i]=null;
			}else if(a.equals("out2")){
				add(b, out2List);
				args[i]=null;
			}else if(a.equals("outlist") || a.equals("out1list")){
				addFileContents(b, out1List);
				args[i]=null;
			}else if(a.equals("out2list")){
				addFileContents(b, out2List);
				args[i]=null;
			}else if(a.equals("outm") || a.equals("outm1") || a.equals("outmapped") || a.equals("outmapped1")){
				add(b, outm1List);
				args[i]=null;
			}else if(a.equals("outm2") || a.equals("outmapped2")){
				add(b, outm2List);
				args[i]=null;
			}else if(a.equals("outu") || a.equals("outu1") || a.equals("outunmapped") || a.equals("outunmapped1")){
				add(b, outu1List);
				args[i]=null;
			}else if(a.equals("outu2") || a.equals("outunmapped2")){
				add(b, outu2List);
				args[i]=null;
			}else if(a.equals("outmlist") || a.equals("outm1list") || a.equals("outmappedlist") || a.equals("outmapped1list")){
				addFileContents(b, outm1List);
				args[i]=null;
			}else if(a.equals("outm2list") || a.equals("outmapped2list")){
				addFileContents(b, outm2List);
				args[i]=null;
			}else if(a.equals("outulist") || a.equals("outu1list") || a.equals("outunmappedlist") || a.equals("outunmapped1list")){
				addFileContents(b, outu1List);
				args[i]=null;
			}else if(a.equals("outu2list") || a.equals("outunmapped2list")){
				addFileContents(b, outu2List);
				args[i]=null;
			}else if(a.equals("outb") || a.equals("outb1") || a.equals("outblack") || a.equals("outblack1") || a.equals("outblacklist") || a.equals("outblacklist1")){
				add(b, outb1List);
				args[i]=null;
			}else if(a.equals("outb2") || a.equals("outblack2") || a.equals("outblacklist2")){
				add(b, outb2List);
				args[i]=null;
			}else if(a.equals("qualityhistogram") || a.equals("qualityhist") || a.equals("qhist")){
				add(b, qhistList);
				args[i]=null;
			}else if(a.equals("matchhistogram") || a.equals("matchhist") || a.equals("mhist")){
				add(b, mhistList);
				args[i]=null;
			}else if(a.equals("inserthistogram") || a.equals("inserthist") || a.equals("ihist")){
				add(b, ihistList);
				args[i]=null;
			}else if(a.equals("bamscript") || a.equals("bs")){
				add(b, bsList);
				args[i]=null;
			}else if(a.equals("append") || a.equals("app")){
				append=Parse.parseBoolean(b);
			}
		}
		
		ArrayList<String> list=new ArrayList<String>();
		for(String s : args){
			if(s!=null){
				list.add(s);
			}
		}
//		return list.toArray(new String[list.size()]);
		return list;
		
	}
	
	/**
	 * Adds comma-separated file paths to the specified list.
	 * Handles null and "null" strings by ignoring them.
	 * @param s Comma-separated string of file paths
	 * @param list Target list to populate with individual paths
	 */
	private static void add(String s, ArrayList<String> list){
		if(s!=null && !"null".equals(s.toLowerCase())){
			String[] sa=s.split(",");
			for(String ss : sa){
				list.add(ss);
			}
		}
	}
	
	/**
	 * Reads file contents line-by-line and adds each line to the specified list.
	 * Used for processing file lists where each line contains a file path.
	 * @param s Path to file containing list of file paths
	 * @param list Target list to populate with file paths from the file
	 */
	private static void addFileContents(String s, ArrayList<String> list){
		if(s!=null && !"null".equals(s.toLowerCase())){
			String[] sa=TextFile.toStringLines(s);
			for(String ss : sa){
				list.add(ss);
			}
		}
	}
	
	/**
	 * Executes the selected mapper for each input file in the batch.
	 * Coordinates sequential processing of file lists, ensuring reference indexing
	 * optimization by loading the index once and reusing it for subsequent runs.
	 * @param base Base arguments to pass to each mapper execution
	 */
	private void execute(ArrayList<String> base){
		for(int i=0; i<in1List.size(); i++){
			ArrayList<String> list=(ArrayList<String>) base.clone();
			
			if(i==0 && ref!=null){list.add("ref="+ref);}
			else if(i>0){list.add("indexloaded=t");}
			
			addToList(list, bsList, "bs", i);
			addToList(list, qhistList, "qhist", i);
			addToList(list, mhistList, "mhist", i);
			addToList(list, ihistList, "ihist", i);
			addToList(list, in1List, "in", i);
			addToList(list, out1List, "out", i);
			addToList(list, outu1List, "outu", i);
			addToList(list, outm1List, "outm", i);
			addToList(list, outb1List, "outb", i);
			addToList(list, in2List, "in2", i);
			addToList(list, out2List, "out2", i);
			addToList(list, outu2List, "outu2", i);
			addToList(list, outm2List, "outm2", i);
			addToList(list, outb2List, "outb2", i);
			
			String[] args=list.toArray(new String[list.size()]);
			if(mapper==null || mapper.equalsIgnoreCase("bbmap")){
				BBMap.main(args);
			}else if(mapper.equalsIgnoreCase("bbmappacbio") || mapper.equalsIgnoreCase("pacbio")){
				BBMapPacBio.main(args);
			}else if(mapper.equalsIgnoreCase("bbmappacbioskimmer") || mapper.equalsIgnoreCase("pacbioskimmer") || mapper.equalsIgnoreCase("skimmer") || mapper.equalsIgnoreCase("bbmapskimmer")){
				BBMapPacBioSkimmer.main(args);
			}else if(mapper.equalsIgnoreCase("bbmap5") || mapper.equalsIgnoreCase("5")){
				BBMap5.main(args);
			}else if(mapper.equalsIgnoreCase("bbmapacc") || mapper.equalsIgnoreCase("acc")){
				BBMapAcc.main(args);
			}else if(mapper.equalsIgnoreCase("bbsplit") || mapper.equalsIgnoreCase("bbsplitter")){
				BBSplitter.main(args);
			}
		}
	}
	
	/**
	 * Adds a key-value argument to the argument list based on batch index.
	 * Handles both indexed file lists and append mode for single-file reuse.
	 *
	 * @param list Target argument list
	 * @param source Source file list containing paths
	 * @param key Parameter name for the argument
	 * @param i Current batch index
	 */
	private void addToList(ArrayList<String> list, ArrayList<String> source, String key, int i){
		if(source.size()>i){
			list.add(key+"="+source.get(i));
		}else if(append && source.size()==1){
			list.add(key+"="+source.get(0));
		}
	}

	private String ref;
	private String mapper="bbmap";

	private ArrayList<String> bsList=new ArrayList<String>();
	private ArrayList<String> qhistList=new ArrayList<String>();
	private ArrayList<String> mhistList=new ArrayList<String>();
	private ArrayList<String> ihistList=new ArrayList<String>();
	
	private ArrayList<String> in1List=new ArrayList<String>();
	private ArrayList<String> out1List=new ArrayList<String>();
	private ArrayList<String> outu1List=new ArrayList<String>();
	private ArrayList<String> outm1List=new ArrayList<String>();
	private ArrayList<String> outb1List=new ArrayList<String>();

	private ArrayList<String> in2List=new ArrayList<String>();
	private ArrayList<String> out2List=new ArrayList<String>();
	private ArrayList<String> outu2List=new ArrayList<String>();
	private ArrayList<String> outm2List=new ArrayList<String>();
	private ArrayList<String> outb2List=new ArrayList<String>();
	
	private boolean append=false;
	
	static PrintStream outstream=System.err;
	
}
