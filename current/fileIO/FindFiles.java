package fileIO;

import java.io.File;
import java.util.ArrayList;

/**
 * Collects absolute file paths using literal prefix, suffix and substring constraints.
 * Matching lowercases names and constructor arguments with the default locale;
 * these are not regexes or general glob patterns. Traversal follows the order
 * returned by File.listFiles(), without sorting or deduplication.
 * Instance searches append to a shared mutable result list; static convenience
 * calls use a fresh instance. Concurrent searches or reconfiguration are unsupported.
 * A null directory listing is logged and leaves partial results; directory-link
 * cycles are not detected.
 * @author Brian Bushnell, Shinobu
 */
public class FindFiles{

	/*--------------------------------------------------------------*/
	/*----------------        Static Entry          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Searches once and prints the same results twice: one path per line, then
	 * a blank line followed by space-separated paths. No usage/count validation
	 * precedes argument indexing. Extra arguments after the optional middle are ignored.
	 * @param args Root, prefix, suffix and optional middle. Literal "null" disables
	 * suffix or middle, but not prefix; "#" disables any constraint. Quote it in a shell.
	 */
	public static void main(String[] args){
		String root=args[0];
		//The historical dot-to-null conversion is inactive; both forms select the working directory.
		String prefix=args[1];
		String suffix=(args[2].equals("null") ? null : args[2]);
		String middle=null;
		if(args.length>3){middle=(args[3].equals("null") ? null : args[3]);}
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

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Initializes literal matching constraints, lowercased using the default locale.
	 * Null and "#" disable a constraint. "*" also disables one, but prefix and suffix
	 * "*" are assertion-rejected to discourage shell expansion; middle "*" is allowed.
	 *
	 * @param pre Literal prefix, null or "#"; "*" only with assertions disabled
	 * @param suf Literal suffix, null or "#"; "*" only with assertions disabled
	 * @param mid Literal substring, null, "#" or "*"
	 */
	public FindFiles(String pre, String suf, String mid){
		assert(!"*".equals(pre)) : "Use # instead of *, which has problems from the command line";
		assert(!"*".equals(suf)) : "Use # instead of *, which has problems from the command line";
		prefix=((pre==null || pre.equals("*") || pre.equals("#")) ? null : pre.toLowerCase());
		suffix=((suf==null || suf.equals("*") || suf.equals("#")) ? null : suf.toLowerCase());
		middle=((mid==null || mid.equals("*") || mid.equals("#")) ? null : mid.toLowerCase());
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Searches with a fresh instance and no substring constraint.
	 *
	 * @param root File or directory path; null uses the working directory
	 * @param prefix Constructor prefix constraint
	 * @param suffix Constructor suffix constraint
	 * @return Fresh mutable result list, possibly partial after a directory-list failure
	 */
	public static ArrayList<String> findFiles(String root, String prefix, String suffix){
		return findFiles(root, prefix, suffix, null);
	}

	/**
	 * Searches with a fresh instance using the constructor's constraint rules.
	 *
	 * @param root File or directory path; null uses the working directory
	 * @param prefix Constructor prefix constraint
	 * @param suffix Constructor suffix constraint
	 * @param mid Constructor substring constraint
	 * @return Fresh mutable result list, possibly partial after a directory-list failure
	 */
	public static ArrayList<String> findFiles(String root, String prefix, String suffix, String mid){
		FindFiles ff=new FindFiles(prefix, suffix, mid);
		return ff.findFiles(root);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Appends matches to this instance's existing results; never clears the list.
	 * @param path File or directory path, or null for the working directory
	 * @return The live {@link #results} list, including matches from earlier calls
	 */
	public ArrayList<String> findFiles(String path){
		findFiles(new File(path==null ? "." : path));
		return results;
	}

	/**
	 * Recursively appends matches without clearing, sorting or deduplicating results.
	 * A null directory listing is reported to stderr and that subtree is skipped.
	 * Non-directories are passed to {@link #consider(File)}, which ignores missing paths.
	 * Links are followed by File's directory tests; no cycle guard is applied.
	 *
	 * @param path Non-null file or directory to search
	 * @return The live accumulated result list, which may be partial after listing failure
	 */
	public ArrayList<String> findFiles(File path){
		//TODO: Probable bug [FIO-118] - directory links can revisit ancestors; traversal has no cycle guard. Policy deferred.
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
	 * Appends one matching absolute path, retaining its original spelling.
	 * Missing paths are ignored; existing paths are assertion-checked as regular files.
	 * The lowercased filename must satisfy every non-null literal constraint.
	 * This method does not lowercase caller replacements of the public constraint fields.
	 * @param in Non-null file to consider; directories are assertion-rejected
	 */
	public void consider(File in){
		if(!in.exists()){return;}
		assert(in.exists()) : in;
		assert(in.isFile());
		String abs=in.getAbsolutePath();
		String abs2=abs.toLowerCase();
		int slashLoc=abs2.lastIndexOf(slash);
		if(slashLoc>-1){abs2=abs2.substring(slashLoc+1);}
		if(prefix!=null && !abs2.startsWith(prefix)){return;}
		if(suffix!=null && !abs2.endsWith(suffix)){return;}
		if(middle!=null && !abs2.contains(middle)){return;}
		results.add(abs);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Live mutable output, appended by all instance searches; callers may clear or replace it. */
	public ArrayList<String> results=new ArrayList<String>();
	/** Literal prefix, initially lowercased by the constructor; null disables it. */
	public String prefix;
	/** Literal suffix, initially lowercased by the constructor; null disables it. */
	public String suffix;
	/**
	 * Literal substring, initially lowercased by the constructor; null disables it.
	 */
	public String middle;
	/**
	 * First character of file.separator, cached when this class is initialized.
	 */
	public static final char slash=System.getProperty("file.separator").charAt(0);
}
