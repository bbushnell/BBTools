package fileIO;

import java.io.File;


/**
 * Utility for batch file compression operations with support for multiple
 * compression formats and recursive directory processing.
 * Provides command-line interface for compressing files using zip or gzip formats.
 * @author Brian Bushnell
 */
public class CompressFiles {
	
	
	/**
	 * Program entry point for batch file compression.
	 * Processes command-line arguments to set compression format and compress files.
	 * Supports "zip", "gzip", and "gz" format specifiers.
	 * @param args Command-line arguments containing format options and file paths
	 */
	public static void main(String[] args){
		for(String s : args){
			if(s.equalsIgnoreCase("zip")){
				zip=true;
				gzip=false;
			}else if(s.equalsIgnoreCase("gzip") || s.equalsIgnoreCase("gz")){
				zip=false;
				gzip=true;
			}else{
				compressFiles(s);
			}
		}
	}
	
	
	/**
	 * Compresses files at the specified path string.
	 * Creates File object and delegates to File-based compression method.
	 * @param path Path to file or directory to compress
	 */
	public static void compressFiles(String path){
		File f=new File(path);
		compressFiles(f);
	}
	
	/**
	 * Recursively compresses files in directories or single files.
	 * For directories, processes all contained files recursively.
	 * For individual files, calls compress method directly.
	 * @param path File or directory to compress
	 */
	public static void compressFiles(File path){
		
		if(path.isDirectory()){
			File[] array=path.listFiles();
			for(File f : array){compressFiles(f);}
		}else{
			compress(path);
		}
		
	}
	
	/**
	 * Compresses a single file using the selected compression format.
	 * Skips already compressed files (.gz, .zip, .bz2 extensions).
	 * Includes filtering logic to skip certain file types like familytree files.
	 * Supports both actual compression and 7-Zip batch command generation.
	 *
	 * @param in Input file to compress
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
	
	
	/** Flag indicating whether to use ZIP compression format */
	public static boolean zip=true;
	/** Flag indicating whether to use GZIP compression format */
	public static boolean gzip=!zip;
	
	/**
	 * Flag controlling whether to print 7-Zip batch commands instead of compressing
	 */
	public static boolean PRINT_7Z_BATCH=true;
	
}
