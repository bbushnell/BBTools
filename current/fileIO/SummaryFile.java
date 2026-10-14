package fileIO;

import java.io.File;

import dna.Data;
import parse.Parse;
import parse.PreParser;

/**
 * Parses legacy genome-summary metadata used by reference-cache reuse checks.
 * RefToIndex and RandomReads3 use compare to match the stored source name, on-disk
 * size and modification time. This is a metadata heuristic: no sequence contents,
 * digest, format version or conversion options are compared.
 *
 * Missing fields retain Java defaults, repeated fields overwrite earlier values,
 * and the public metadata remains mutable after construction. Unknown data keys
 * and malformed values throw; parsing does not validate all required metadata.
 * The diagnostic main is unfinished and does not perform the comparison.
 *
 * @author Brian Bushnell
 * @contributor Shinobu
 * @date Mar 11, 2013
 */
public class SummaryFile{
	
	/**
	 * Legacy argument-parser stub; it does not open a summary or compare a reference.
	 * No arguments print usage and exit successfully; other arguments go through
	 * PreParser and local positional/key-value parsing before returning.
	 * @param args Legacy summary/ref arguments
	 */
	public static void main(String[] args){
		//TODO: Probable bug - FIO-065: this entry point parses options but never validates
		//or compares anything. Keep dormant until a supported diagnostic use is assigned.
		if(args.length==0){
			System.out.println("Usage: SummaryFile <summary file> <reference fasta>");
			System.exit(0);
		}

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object(){}.getClass().getEnclosingClass(), false);
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
	 * Compares source-name/path, on-disk byte length and modification time only.
	 * Equal raw names bypass canonical-path comparison. A missing stdin-prefixed
	 * name returns false; another missing name asserts when assertions are enabled.
	 * Caught Exceptions print diagnostics and return false; AssertionError propagates.
	 * Does not compare sequence content, summary version or other metadata fields.
	 * @param refName Candidate reference filename
	 * @return Whether the three stored metadata values match
	 */
	public boolean compare(final String refName){
		try{
			File ref=new File(refName);
			if(!ref.exists()){
				if(refName.startsWith("stdin")){return false;}
				else{
					assert(false) : "No such file: "+refName;
				}
			}
			if(!refName.equals(source) && !ref.getCanonicalPath().equals(new File(source).getCanonicalPath())){
				return false;
				
			}
			if(bytes!=ref.length()){
				return false;
			}
			if(modified!=ref.lastModified()){
				return false;
			}
		}catch(Exception e){
			e.printStackTrace();
			return false;
		}
		return true;
	}
	
	/**
	 * Parses an existing summary and applies the instance metadata comparison.
	 * Missing summary returns false. Construction/parse failures propagate rather
	 * than being converted to a cache miss; the instance method catches its own
	 * ordinary comparison Exceptions. No reference contents are inspected.
	 * @param summaryName Summary filename, tested with File.exists before opening
	 * @param refName Candidate reference filename, asserted nonnull
	 * @return Metadata match, or false when the summary is absent
	 */
	public static boolean compare(final String summaryName, final String refName){
		assert(refName!=null) : "Null reference file name.";
		if(!new File(summaryName).exists()){
			return false;
		}
		SummaryFile sf=new SummaryFile(summaryName);
		return sf.compare(refName);
	}
	
	/** @return Configured root plus current build number and /summary.txt; performs no I/O. */
	public static String getName(){
		return getName(Data.GENOME_BUILD);
	}
	
	/** Concatenates Data.ROOT_GENOME and build without checking existence or validity.
	 * @param build Genome build number
	 * @return Configured summary path */
	public static String getName(int build){
		return Data.ROOT_GENOME+build+"/summary.txt";
	}
	
	/**
	 * Streams a summary through TextFile with subprocess decompression disabled.
	 * TextFile skips its blank-line predicate; retained lines and values are not
	 * trimmed. A case-sensitive #Version prefix sets the version (zero without a
	 * following field); other #-prefixed lines are ignored. Data keys are matched
	 * case-insensitively and use the first two tab-separated fields. String.split
	 * drops trailing empty fields; fields after the value are ignored.
	 * Repeated keys overwrite, absent keys retain defaults, and unknown keys throw.
	 * Numeric values have no added range validation; chroms narrows a parsed long.
	 * Closes the input after parsing or failure; reported close errors also throw.
	 * A parse/read failure stays primary, with distinct cleanup failures suppressed.
	 * @param path Summary input name passed to TextFile
	 */
	public SummaryFile(String path){
		summaryFname=path;
		String s;
		TextFile tf=new TextFile(summaryFname, false);
		Throwable failure=null;
		try{
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
					//FIO-064: preserve the long totals written by FastaToChromArrays2.
					else if(a.equalsIgnoreCase("contigs")){contigs=Long.parseLong(b);}
					else if(a.equalsIgnoreCase("scaffolds")){scaffolds=Long.parseLong(b);}
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
		}catch(RuntimeException | Error t){
			failure=t;
			throw t;
		}finally{
			//FIO-063: exceptional parsing still owns this reader; retain the original failure.
			try{
				if(tf.close()){throw new RuntimeException("Error closing summary file: "+summaryFname);}
			}catch(RuntimeException | Error t){
				if(failure==null){throw t;}
				if(failure!=t){failure.addSuppressed(t);}
			}
		}
	}

	/** Path to the parsed summary.txt file. */
	public final String summaryFname;

	/** Number of chromosomes in the genome build. */
	public int chroms;
	/** Contig count, parsed as a long to match the summary producer. */
	public long contigs;
	/** Scaffold count, parsed as a long to match the summary producer. */
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
	/** Source File.lastModified value in milliseconds since the epoch, used for matching. */
	public long modified;
	/** True if scaffold names carry build-specific prefixes. */
	public boolean scafprefixes;
	
}
