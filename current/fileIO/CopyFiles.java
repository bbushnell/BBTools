package fileIO;

import java.io.File;

/**
 * Copies chromosome text files beside their originals with a .flow extension.
 * Each input is a file or a directory tree; only basenames starting with "chr"
 * and ending with ".txt" match, with case-sensitive comparisons. The terminal
 * extension changes, but the file contents are copied without format conversion.
 * Despite the method names and historical "Renaming" message, originals remain.
 * <p>
 * Traversal follows directory links without cycle detection (FIO-111). A thrown
 * exception or assertion failure stops traversal; prior copies are not rolled back. Destination
 * nonexistence is assertion-checked by {@link ReadWrite#copyFile(String, String)}.
 * @author Brian Bushnell, Shinobu
 */
public class CopyFiles{

	/**
	 * Processes each argument as a literal file or directory path, in argument order.
	 * No flags are parsed; no arguments means no work.
	 * @param args File or directory paths to process
	 */
	public static void main(String[] args){
		for(String s : args){renameFiles(s);}
	}

	/**
	 * Wraps a path and delegates to the recursive file-based overload.
	 * @param path Non-null file or directory path
	 */
	public static void renameFiles(String path){
		File f=new File(path);
		renameFiles(f);
	}

	/**
	 * Visits directory children recursively, or passes a nondirectory to {@link #rename(File)}.
	 * Child order is the order returned by File.listFiles; it is not sorted here.
	 * Directory links are followed, with the unresolved cycle risk noted below.
	 * @param path Non-null file or directory to process
	 * @throws RuntimeException If directory contents cannot be listed, or a copy fails
	 */
	public static void renameFiles(File path){
		//TODO: Probable bug [FIO-111] - directory links can revisit ancestors; traversal has no cycle guard.
		if(path.isDirectory()){
			File[] array=path.listFiles();
			if(array==null){throw new RuntimeException("Failed to list directory contents (I/O error): "+path.getAbsolutePath());}
			for(File f : array){renameFiles(f);}
		}else{
			rename(path);
		}
	}

	/**
	 * Copies a matching file to a sibling path ending in .flow; skips other basenames.
	 * Only the terminal .txt suffix is replaced, leaving parent directories and any
	 * earlier .txt substring intact. The source is preserved and its bytes are copied
	 * through ReadWrite.copyFile; the printed "Renaming" message is historical.
	 * Input existence and regular-file status are assertion-checked before filtering.
	 * The destination must not exist; no parent directories are created.
	 * @param in Non-null existing regular file
	 * @throws RuntimeException If the underlying copy reports an I/O failure
	 */
	public static void rename(File in){
		assert(in.exists());
		assert(in.isFile());
		String abs=in.getAbsolutePath();
		int dot=abs.lastIndexOf('.');
		// Historical alternate naming used a parent-path name and a .txt suffix; it is inactive.
		String fname=in.getName();
		if(fname.startsWith("chr") && fname.endsWith(".txt")){
			String out=abs.substring(0, dot)+".flow";
			assert(!out.equals(abs)) : out+", "+abs;
			System.out.println("Renaming "+abs+" to "+out);
			ReadWrite.copyFile(abs, out);
		}
	}
}
