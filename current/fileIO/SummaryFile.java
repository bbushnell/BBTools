package fileIO;

import java.io.File;

import dna.Data;
import parse.Parse;
import parse.PreParser;

/**
 * Parses and validates genome summary files, extracting critical metadata about
 * genome builds and file characteristics. Tests summary file compatibility with
 * reference FASTA files based on date, size, and source path matching.
 *
 * @author Brian Bushnell
 * @date Mar 11, 2013
 */
public class SummaryFile {
	
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
	
	/** @return true if the given reference FASTA matches this summary by canonical path, byte size, and last-modified time. */
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
	
	/** Loads the summary file and compares it against the reference. @return true if they match (false if the summary file is missing). */
	public static boolean compare(final String summaryName, final String refName){
		assert(refName!=null) : "Null reference file name.";
		if(!new File(summaryName).exists()){
//			assert(false);
			return false;
		}
		SummaryFile sf=new SummaryFile(summaryName);
		return sf.compare(refName);
	}
	
	/** @return Path to the current genome build's summary.txt. */
	public static String getName(){
		return getName(Data.GENOME_BUILD);
	}
	
	/** @return Path to the given genome build's summary.txt. */
	public static String getName(int build){
		return Data.ROOT_GENOME+build+"/summary.txt";
	}
	
	/** Parses the summary file at the given path, populating the metadata fields. */
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

	/** Path to the parsed summary.txt file. */
	public final String summaryFname;

	/** Number of chromosomes in the genome build. */
	public int chroms;
	/** Number of contigs. */
	public long contigs;
	/** Number of scaffolds. */
	public long scaffolds;
	/** Padding length (N runs) inserted between scaffolds. */
	public int interpad;
	/** Total base count (defined + undefined). */
	public long bases;
	/** Count of defined (non-N) bases. */
	public long definedBases;
	/** Count of undefined (N) bases. */
	public long undefinedBases;
	/** Genome name. */
	public String name;
	/** Path to the source reference FASTA this summary was built from. */
	public String source;
	/** Summary-file format version. */
	public int version;
	/** Size in bytes of the source reference file (used to detect staleness). */
	public long bytes;
	/** Last-modified timestamp of the source reference file (used to detect staleness). */
	public long modified;
	/** True if scaffold names carry build-specific prefixes. */
	public boolean scafprefixes;
	
}
