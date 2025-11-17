package fileIO;

import java.io.File;


/**
 * Utility class for processing and renaming files within directories.
 * Recursively traverses directories and copies files matching specific criteria,
 * particularly targeting chromosome text files for format conversion.
 * @author Brian Bushnell
 */
public class CopyFiles {
	
	
	/** Program entry point that processes each command-line argument as a file path.
	 * @param args File or directory paths to process */
	public static void main(String[] args){
		for(String s : args){
			renameFiles(s);
		}
	}
	
	
	/**
	 * Processes files at the specified path for renaming operations.
	 * Creates a File object and delegates to the File-based method.
	 * @param path String path to file or directory to process
	 */
	public static void renameFiles(String path){
		File f=new File(path);
		renameFiles(f);
	}
	
	/**
	 * Recursively processes files and directories for renaming operations.
	 * If the path is a directory, recursively processes all contained files.
	 * If the path is a file, delegates to rename method for processing.
	 * @param path File or directory to process recursively
	 */
	public static void renameFiles(File path){
		
		if(path.isDirectory()){
			File[] array=path.listFiles();
			for(File f : array){renameFiles(f);}
		}else{
			rename(path);
		}
		
	}
	
	/**
	 * Renames individual files based on specific criteria.
	 * Copies files that start with "chr" and end with ".txt" to new files
	 * with ".flow" extension, preserving the original file.
	 * @param in Input file to process (must exist and be a regular file)
	 */
	public static void rename(File in){
		assert(in.exists());
		assert(in.isFile());
		String abs=in.getAbsolutePath();
		
		
		int dot=abs.lastIndexOf('.');
		int slash=abs.lastIndexOf('/');
		
//		String[] split=Person.parsePath(abs.substring(0, slash));
//		String name=split[0];
//		String out=abs.substring(0, dot)+"_"+name+".txt";
		
		
		
		String fname=abs.substring(slash+1);
		
//		System.out.println(fname);
		
		
		if(fname.startsWith("chr") && fname.endsWith(".txt")){
			
			String out=abs.replace(".txt", ".flow");
			assert(!out.equals(abs)) : out+", "+abs;
			
			System.out.println("Renaming "+abs+" to "+out);
			ReadWrite.copyFile(abs, out);
		}
	}
	
}
