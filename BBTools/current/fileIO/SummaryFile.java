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
	
	public static boolean compare(final String summaryName, final String refName){
		assert(refName!=null) : "Null reference file name.";
		if(!new File(summaryName).exists()){
//			assert(false);
			return false;
		}
		SummaryFile sf=new SummaryFile(summaryName);
		return sf.compare(refName);
	}
	
	public static String getName(){
		return getName(Data.GENOME_BUILD);
	}
	
	public static String getName(int build){
		return Data.ROOT_GENOME+build+"/summary.txt";
	}
	
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

	public final String summaryFname;

	public int chroms;
	public long contigs;
	public long scaffolds;
	public int interpad;
	public long bases;
	public long definedBases;
	public long undefinedBases;
	public String name;
	public String source;
	public int version;
	public long bytes;
	public long modified;
	public boolean scafprefixes;
	
}
