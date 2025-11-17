package fileIO;

import java.io.File;
import java.util.ArrayList;


/**
 * Utility class for finding files in directory trees based on filename patterns.
 * Searches files recursively and filters by prefix, suffix, and middle string matching.
 * Provides both static convenience methods and instance-based searching.
 * @author Brian Bushnell
 */
public class FindFiles {
	
	
	/**
	 * Command-line entry point for file finding operations.
	 * Arguments: root_directory prefix suffix [middle]
	 * Prints found files in two formats for comparison.
	 * @param args Command-line arguments: [0]=root, [1]=prefix, [2]=suffix, [3]=middle (optional)
	 */
	public static void main(String[] args){
		
		String root=args[0];
//		if(root.equals(".")){root=null;}
		String prefix=args[1];
		String suffix=(args[2].equals("null") ? null : args[2]);
		String middle=null;
		
		if(args.length>3){
			middle=(args[3].equals("null") ? null : args[3]);
		}

		boolean NEWLINE=true;
		boolean BOTH=true;

		ArrayList<String> results=findFiles(root, prefix, suffix, middle);
		for(String s : results){
			if(NEWLINE){
				System.out.println(s);
			}else{
				System.out.print(s+" ");
			}
		}


		if(BOTH){
			System.out.println();
			NEWLINE=!NEWLINE;
			for(String s : results){
				if(NEWLINE){
					System.out.println(s);
				}else{
					System.out.print(s+" ");
				}
			}
		}
	}
	
	
	/**
	 * Constructs a FindFiles instance with search patterns.
	 * Converts wildcard patterns (*,#) to null and normalizes to lowercase.
	 * Validates that * wildcards are replaced with # for command-line safety.
	 *
	 * @param pre Prefix pattern to match (null, *, or # for no constraint)
	 * @param suf Suffix pattern to match (null, *, or # for no constraint)
	 * @param mid Middle string pattern to match (null, *, or # for no constraint)
	 */
	public FindFiles(String pre, String suf, String mid){
		assert(!"*".equals(pre)) : "Use # instead of *, which has problems from the command line";
		assert(!"*".equals(suf)) : "Use # instead of *, which has problems from the command line";
		prefix=((pre==null || pre.equals("*") || pre.equals("#")) ? null : pre.toLowerCase());
		suffix=((suf==null || suf.equals("*") || suf.equals("#")) ? null : suf.toLowerCase());
		middle=((mid==null || mid.equals("*") || mid.equals("#")) ? null : mid.toLowerCase());
	}
	
	/**
	 * Static convenience method to find files with prefix and suffix patterns.
	 * Creates temporary FindFiles instance and searches from root directory.
	 *
	 * @param root Root directory path to search from
	 * @param prefix Filename prefix to match
	 * @param suffix Filename suffix to match
	 * @return List of absolute paths of matching files
	 */
	public static ArrayList<String> findFiles(String root, String prefix, String suffix){
		return findFiles(root, prefix, suffix, null);
	}
	
	/**
	 * Static convenience method to find files with prefix, suffix, and middle patterns.
	 * Creates temporary FindFiles instance and searches from root directory.
	 *
	 * @param root Root directory path to search from
	 * @param prefix Filename prefix to match
	 * @param suffix Filename suffix to match
	 * @param mid Middle string pattern to match within filename
	 * @return List of absolute paths of matching files
	 */
	public static ArrayList<String> findFiles(String root, String prefix, String suffix, String mid){
		FindFiles ff=new FindFiles(prefix, suffix, mid);
		return ff.findFiles(root);
	}
	
	/**
	 * Searches for matching files starting from specified path.
	 * Uses current directory if path is null.
	 * @param path Starting directory path (uses "." if null)
	 * @return List of absolute paths of files matching this instance's patterns
	 */
	public ArrayList<String> findFiles(String path){
		findFiles(new File(path==null ? "." : path));
		return results;
	}
	
	/**
	 * Recursively searches for files matching the configured patterns.
	 * If path is directory, recursively searches all contents.
	 * If path is file, considers it for pattern matching.
	 *
	 * @param path File or directory to search
	 * @return List of absolute paths of matching files
	 */
	public ArrayList<String> findFiles(File path){
		
		if(path.isDirectory()){
			File[] array=path.listFiles();
			if(array==null){System.err.println("null contents for "+path.getAbsolutePath());}
			else{for(File f : array){findFiles(f);}}
		}else{
			consider(path);
		}
		return results;
	}
	
	/**
	 * Evaluates a single file against the configured search patterns.
	 * Checks filename (not full path) against prefix, suffix, and middle constraints.
	 * Adds matching files to results list as absolute paths.
	 * @param in File to evaluate for pattern matching
	 */
	public void consider(File in){
//		System.out.println("Considering "+in.getAbsolutePath()+" versus '"+prefix+"' '"+suffix+"'");
		if(!in.exists()){return;}
		assert(in.exists()) : in;
		assert(in.isFile());
		String abs=in.getAbsolutePath();
//		System.out.println("Considering "+abs);
		String abs2=abs.toLowerCase();
		int slashLoc=abs2.lastIndexOf(slash);
		if(slashLoc>-1){
			abs2=abs2.substring(slashLoc+1);
		}
//		System.out.println("a");
		if(prefix!=null && !abs2.startsWith(prefix)){return;}
//		System.out.println("b");
		if(suffix!=null && !abs2.endsWith(suffix)){return;}
//		System.out.println("c");
		
		if(middle!=null && !abs2.contains(middle)){return;}
		
		results.add(abs);
	}
	
	
	/** List storing absolute paths of files that match search patterns */
	public ArrayList<String> results=new ArrayList<String>();
	/** Filename prefix pattern to match (null means no constraint) */
	public String prefix;
	/** Filename suffix pattern to match (null means no constraint) */
	public String suffix;
	/**
	 * Middle string pattern that must be contained in filename (null means no constraint)
	 */
	public String middle;
	/**
	 * Platform-specific file separator character obtained from system properties
	 */
	public static final char slash=System.getProperty("file.separator").charAt(0);
	
}
