package fileIO;

import java.io.File;

/**
 * Traverses files to print Windows 7-Zip commands or perform compression copies.
 * The default PRINT_7Z_BATCH mode prints commands only, always with a .gz suffix.
 * With that shared flag disabled programmatically, zip selects .zip or .gz output
 * and ReadWrite.copyFile performs the copy through its compression-aware factories.
 *
 * Paths are processed in argument order; format keywords affect subsequent paths.
 * Existing lowercase .gz/.zip/.bz2 suffixes and paths containing familytree after
 * lowercasing are skipped. Input files remain in place. Printed commands use the
 * hardcoded executable path and emit file operands without added quoting.
 * Shared mode flags should be configured before traversal; traversal is not synchronized.
 * @author Brian Bushnell, Shinobu
 */
public class CompressFiles {

	/*--------------------------------------------------------------*/
	/*----------------           Startup            ----------------*/
	/*--------------------------------------------------------------*/

	/** Processes bare zip/gzip/gz keywords and paths from left to right.
	 * Keywords are case-insensitive; other arguments are paths. There is no CLI
	 * switch for PRINT_7Z_BATCH, so the default invocation only prints commands.
	 * @param args Format keywords and file or directory paths
	 */
	public static void main(String[] args){
		for(String s : args){
			if(s.equalsIgnoreCase("zip")){
				zip=true;
			}else if(s.equalsIgnoreCase("gzip") || s.equalsIgnoreCase("gz")){
				zip=false;
			}else{
				compressFiles(s);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Compression         ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a File for the supplied path and delegates traversal.
	 * @param path Nonnull file or directory path
	 */
	public static void compressFiles(String path){
		File f=new File(path);
		compressFiles(f);
	}

	/** Recurses through directories in their returned listing order, or handles one file.
	 * @param path Nonnull file or directory
	 * @throws RuntimeException If a directory listing cannot be read
	 */
	public static void compressFiles(File path){
		if(path.isDirectory()){
			File[] array=path.listFiles();
			if(array==null){throw new RuntimeException("Failed to list directory contents (I/O error): "+path.getAbsolutePath());}
			for(File f : array){compressFiles(f);}
		}else{
			compress(path);
		}
	}

	/** Applies suffix/path filters, then prints a command or delegates a compression copy.
	 * Lowercase .gz/.zip/.bz2 suffixes return before reporting the absolute path on
	 * stderr; the familytree filter runs afterward. Command mode ignores zip and
	 * prints a .gz destination. Copy mode selects its suffix from zip and delegates
	 * destination validation, compression and stream handling to ReadWrite.copyFile.
	 * @param in Existing regular file, checked by assertions
	 */
	public static void compress(File in){
		assert(in.exists());
		assert(in.isFile());
		String abs=in.getAbsolutePath();
//		System.out.println("Considering "+abs);
		if(abs.endsWith(".gz") || abs.endsWith(".zip") || abs.endsWith(".bz2")){return;}

//		if(!abs.contains("custom_summary_") || !abs.endsWith("Gene_build36.txt")){return;} //TODO ***TEMPORARY***
		System.err.println(abs);
//		if(!abs.endsWith(".gvla")){return;} //TODO ***TEMPORARY***
//		if(!abs.endsWith(".gvla") ||
//				!(abs.contains("seqGene") || abs.contains("refGene") || abs.contains("unionGene"))){return;} //TODO ***TEMPORARY***
		if(abs.toLowerCase().contains("familytree")){return;} //TODO ***TEMPORARY***

		if(PRINT_7Z_BATCH){
			//-mx=4 is fast; -mx=5 or 6 is slow; 7+ is very slow.
//			System.out.println("C:"+Data.SLASH+"\"Program Files\""+Data.SLASH+"7-Zip"+Data.SLASH+"7z a -mx=4 "+abs+".zip "+abs);
			System.out.println("C:\\\"Program Files\"\\7-Zip\\7z a -mx=4 "+abs+".gz "+abs);
		}else{
			System.out.println("Compressing "+abs+" to "+(zip ? "zip" : "gz"));
			ReadWrite.copyFile(abs, abs+(zip ? ".zip" : ".gz"));
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Configuration         ----------------*/
	/*--------------------------------------------------------------*/

	/** Selects .zip rather than .gz for direct copies; command mode does not use it. */
	public static boolean zip=true;

	/** Prints Windows 7-Zip commands instead of copying; configured programmatically. */
	public static boolean PRINT_7Z_BATCH=true;

}
