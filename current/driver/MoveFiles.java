package driver;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import shared.Tools;


/**
 * Utility for organizing files by chromosome number into subdirectories.
 * Processes files in a root directory and moves those matching chromosome patterns
 * (chr1-chr22) into corresponding subdirectories.
 * @author Brian Bushnell
 */
public class MoveFiles {
	
	
	/**
	 * Main entry point for file organization by chromosome.
	 * Takes a root directory path as argument, scans for files matching chromosome
	 * patterns (chr1-chr22), creates subdirectories for each chromosome, and moves
	 * matching files into appropriate subdirectories.
	 * @param args Command line arguments where args[0] is the root directory path
	 */
	public static void main(String[] args){
		
		String root=args[0].replace('\\', '/');
		
		File dir=new File(root);
		assert(dir.exists());
		assert(dir.isDirectory());

		
		File[] files=dir.listFiles();
		//NOTE [driver/MoveFiles#001] LOW/dev (no .sh, no callers): (a) despite the class name "MoveFiles" and the javadoc
		//("moves matching files"), copyFile() only COPIES — the source is never deleted, so this is a copy, not a move.
		//(b) dir.listFiles() can return null (I/O error) even past the exists()/isDirectory() asserts → NPE in the for-each
		//below (unguarded; and the asserts vanish under -da). (c) name-matching (strip ext → strip trailing non-digits →
		//endsWith("chr"+chrom)) assumes the chrom number is the last digit-run — fragile but endsWith avoids chr1/chr11
		//collisions. args[0] unguarded. Dead one-off → LOW.
		for(int chrom=1; chrom<=22; chrom++){
			
			String key="chr"+chrom;
			
			File dest=new File(root+"/"+key);
			if(!dest.exists()){
				dest.mkdir();
			}
			
			for(File f : files){
				String name=f.getName();
				if(name.contains("/")){
					name=name.substring(name.lastIndexOf("/")+1);
				}
				final String name2=name;
				
				if(name.contains(".")){
					name=name.substring(0,name.lastIndexOf("."));
				}
				
				while(name.length()>1 && !Tools.isDigit(name.charAt(name.length()-1))){
					name=name.substring(0, name.length()-1);
				}
				name=name.toLowerCase();
				
				if(f.isFile() && name.endsWith("chr"+chrom)){
					copyFile(f.getAbsolutePath(), dest.getAbsolutePath()+"/"+name2);
				}
			}
		}
		
	}

	
	/**
	 * @param srFile
	 * @param dtFile
	 * {@link from http://www.roseindia.net/java/beginners/CopyFile.shtml}
	 */
	private static void copyFile(String src, String dst){
//		assert(false) : src+" -> "+dst;
		try{
			File f1 = new File(src);
			File f2 = new File(dst);
			InputStream in = new FileInputStream(f1);
			//For Append the file.
			//	      OutputStream out = new FileOutputStream(f2,true);

			//For Overwrite the file.
			OutputStream out = new FileOutputStream(f2);

			byte[] buf = new byte[16384];
			int len;
			while ((len = in.read(buf)) > 0){
				out.write(buf, 0, len);
			}
			in.close();
			out.close();
		}catch(FileNotFoundException e){
			throw new RuntimeException(e);
		}catch(IOException e){
			throw new RuntimeException(e);
		}
	}
	
	
}
