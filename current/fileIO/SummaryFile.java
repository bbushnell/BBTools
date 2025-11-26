package fileIO;

import java.io.File;

import dna.Data;
import shared.Parse;
import shared.PreParser;

/**
 * Tests to see if a summary file matches a reference fasta file, based on date, size, and name
 * @author Brian Bushnell
 * @date Mar 11, 2013
 *
 */
public class SummaryFile {
	
	/**
	 * Program entry point for command-line summary file validation.
	 * Parses arguments to extract summary file and reference FASTA paths.
	 * @param args Command-line arguments including summary and reference files
	 */
	public static void main(String[] args){
		if(args.length==0){
			System.out.println("Usage: SummaryFile <summary file> <reference fasta>");
			System.exit(0);
		}

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		String summary=null, ref=null;
		
		for(int i=0; i<args.length; i++){

			if(args[i].contains("=")){
				final String arg=args[i];
				final String[] split=arg.split("=");
				String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1] : null;
				
				if(a.equals("summary")){
					summary=b;
				}else if(a.equals("ref") || a.equals("reference")){
					ref=b;
				}else{
					throw new RuntimeException("Unknown parameter: "+args[i]);
				}

			}else{
				if(args[i].endsWith("summary.txt")){
					summary=args[i];
				}else{
					ref=args[i];
				}
			}
		}
		
		if(summary==null && args.length>0){
			summary=args[0];
		}
		
		if(summary==null){
			System.out.println("Usage: SummaryFile <summary file> <reference fasta>");
			System.exit(0);
		}
		
		if(ref==null){
			
		}
	}
	
	/**
	 * Validates summary file metadata against a reference FASTA file.
	 * Checks source path, file size, and last modified timestamp for exact matches.
	 * Returns false for stdin inputs or when any validation check fails.
	 *
	 * @param refName Path to reference FASTA file to validate against
	 * @return true if summary matches reference file exactly, false otherwise
	 */
	public boolean compare(final String refName){
		try {
			File ref=new File(refName);
			if(!ref.exists()){
				if(refName.startsWith("stdin")){return false;}
				else{
					assert(false) : "No such file: "+refName;
				}
			}
//			if(!refName.equals(source) && !Files.isSameFile(ref.toPath(), new File(source).toPath())){ //This is Java-7 specific.
////				assert(false) : refName+", "+source+": "+(Files.isSameFile(ref.toPath(), new File(source).toPath()))+
////						"\n"+ref.getCanonicalPath()+", "+new File(source).getCanonicalPath()+": "+(ref.getCanonicalPath().equals(new File(source).getCanonicalPath()));
//				return false;
//
//			}
			if(!refName.equals(source) && !ref.getCanonicalPath().equals(new File(source).getCanonicalPath())){
//				assert(false) : refName+", "+source+": "+(Files.isSameFile(ref.toPath(), new File(source).toPath()))+
//						"\n"+ref.getCanonicalPath()+", "+new File(source).getCanonicalPath()+": "+(ref.getCanonicalPath().equals(new File(source).getCanonicalPath()));
				return false;
				
			}
			if(bytes!=ref.length()){
//				assert(false) : bytes+", "+ref.length();
				return false;
			}
			if(modified!=ref.lastModified()){
//				assert(false) : modified+", "+ref.lastModified();
				return false;
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return false;
		}
		return true;
	}
	
	/**
	 * Static comparison method that creates a SummaryFile instance and validates
	 * it against a reference file. Returns false if summary file doesn't exist.
	 *
	 * @param summaryName Path to summary file to load and validate
	 * @param refName Path to reference FASTA file for validation
	 * @return true if summary matches reference file, false otherwise
	 */
	public static boolean compare(final String summaryName, final String refName){
		assert(refName!=null) : "Null reference file name.";
		if(!new File(summaryName).exists()){
//			assert(false);
			return false;
		}
		SummaryFile sf=new SummaryFile(summaryName);
		return sf.compare(refName);
	}
	
	/** Generates default summary file path using current genome build from Data.
	 * @return Standard summary file path for current genome build */
	public static String getName(){
		return getName(Data.GENOME_BUILD);
	}
	
	/**
	 * Generates summary file path for specified genome build number.
	 * Constructs path using Data.ROOT_GENOME + build + "/summary.txt" format.
	 * @param build Genome build number for path generation
	 * @return Summary file path for specified build
	 */
	public static String getName(int build){
		return Data.ROOT_GENOME+build+"/summary.txt";
	}
	
	/**
	 * Constructs SummaryFile by parsing tab-delimited summary text file.
	 * Reads file line-by-line, extracting metadata values including chromosome
	 * counts, base counts, version info, and file characteristics. Handles both
	 * comment lines (#Version) and data lines (key\tvalue format).
	 *
	 * @param path Path to summary file to parse and load
	 */
	public SummaryFile(String path){
		summaryFname=path;
		String s;
		TextFile tf=new TextFile(summaryFname, false);
		for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(s.charAt(0)=='#'){
				if(s.startsWith("#Version")){
					String[] split=s.split("\t");
					version=(split.length>1 ? Integer.parseInt(split[1]) : 0);
				}
			}else{
				String[] split=s.split("\t");
				String a=split[0];
				String b=split[1];
				if(a.equalsIgnoreCase("chroms")){chroms=(int)Long.parseLong(b);}
				else if(a.equalsIgnoreCase("bases")){bases=Long.parseLong(b);}
				else if(a.equalsIgnoreCase("version")){version=Integer.parseInt(b);}
				else if(a.equalsIgnoreCase("defined")){definedBases=Long.parseLong(b);}
				else if(a.equalsIgnoreCase("contigs")){contigs=Integer.parseInt(b);}
				else if(a.equalsIgnoreCase("scaffolds")){scaffolds=Integer.parseInt(b);}
				else if(a.equalsIgnoreCase("interpad")){interpad=Integer.parseInt(b);}
				else if(a.equalsIgnoreCase("undefined")){undefinedBases=Long.parseLong(b);}
				else if(a.equalsIgnoreCase("name")){name=b;}
				else if(a.equalsIgnoreCase("source")){source=b;}
				else if(a.equalsIgnoreCase("bytes")){bytes=Long.parseLong(b);}
				else if(a.equalsIgnoreCase("last modified")){modified=Long.parseLong(b);}
				else if(a.equalsIgnoreCase("scafprefixes")){scafprefixes=Parse.parseBoolean(b);}
				else{throw new RuntimeException("In file "+tf.name+": Unknown term "+s);}
			}
		}
		tf.close();
	}

	/** Path to the summary file that was parsed */
	public final String summaryFname;

	/** Number of chromosomes in the genome build */
	public int chroms;
	/** Number of contigs in the genome assembly */
	public long contigs;
	/** Number of scaffolds in the genome assembly */
	public long scaffolds;
	/** Inter-scaffold padding length used in assembly */
	public int interpad;
	/** Total number of bases in the genome including undefined bases */
	public long bases;
	/** Number of defined bases (A, T, G, C) excluding N's and gaps */
	public long definedBases;
	/** Number of undefined bases (N's, gaps, and ambiguous bases) */
	public long undefinedBases;
	/** Name identifier for the genome build */
	public String name;
	/** Source file path that the summary was generated from */
	public String source;
	/** Version number of the genome build or summary format */
	public int version;
	/** File size in bytes of the original reference file */
	public long bytes;
	/** Last modified timestamp of the original reference file */
	public long modified;
	/** Whether scaffold names include standard prefixes in the assembly */
	public boolean scafprefixes;
	
}
