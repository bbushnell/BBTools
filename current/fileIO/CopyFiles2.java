package fileIO;

import java.io.File;

import shared.Timer;


/**
 * Utility for copying files from input directories to output directories with filtering.
 * Supports directory traversal, file validation, and automatic compression for large TSV files.
 * Filters files based on directory names, filename patterns, and exclusion lists.
 * @author Brian Bushnell
 */
public class CopyFiles2 {
	
	
	/**
	 * Program entry point for file copying utility.
	 * Accepts optional command-line arguments for source and destination paths.
	 * @param args Command-line arguments: [input_root] [output_root] (optional)
	 */
	public static void main(String[] args){
		
		Timer t=new Timer();
		
		if(args.length>0){
			assert(args.length==2);
			inRoots=new String[] {args[0]};
			outRoot=args[1];
		}
		
		for(String inRoot : inRoots){
			copyFiles(inRoot, outRoot);
		}
		
		t.stop();
		System.out.println("Time:\t"+t);
	}
	
	
	/**
	 * Copies files from input path to output path using string paths.
	 * Creates File objects and delegates to File-based copy method.
	 * @param in Input directory or file path
	 * @param out Output directory or file path
	 */
	public static void copyFiles(String in, String out){
		File fin=new File(in);
		File fout=new File(out);
		copyFiles(fin, fout);
	}
	
	/**
	 * Recursively copies files and directories with path filtering and ASM removal.
	 * Skips files matching patterns in badNames array.
	 * Removes "ASM" components from output paths.
	 * Creates directories as needed and processes all subdirectories.
	 *
	 * @param in Input File object (file or directory)
	 * @param out Output File object (destination path)
	 */
	public static void copyFiles(File in, File out){
		
		String abs=in.getAbsolutePath();
		for(String s : badNames){
			if(abs.matches(s)){
				return;
			}
		}
		
		{
			String temp=out.getAbsolutePath();
			if(temp.endsWith("\\ASM")){
				temp=temp.replace("\\ASM", "");
			}else if(temp.contains("\\ASM\\")){
				temp=temp.replace("\\ASM\\", "");
			}
			out=new File(temp);
		}
		
		if(in.isDirectory()){
//			System.out.println("PATH: "+in.getAbsolutePath());
			if(!out.exists()){
				out.mkdir();
			}
			
			File[] array=in.listFiles();
			for(File f : array){
//				String outname=f.getAbsolutePath().replace(inRoot, outRoot);
				
				String outname=out.getAbsolutePath()+"\\"+f.getName();
				
				File f2=new File(outname);
				copyFiles(f, f2);
			}
		}
		
		else{
			copyFile(in, out);
		}
		
	}
	
	/**
	 * Copies a single file if it passes validation filters.
	 * Validates against badNames, dirNames, and fileNames arrays.
	 * Automatically compresses TSV files by adding .zip extension.
	 * Skips existing destination files.
	 *
	 * @param in Input file to copy
	 * @param out Destination file path
	 */
	public static void copyFile(File in, File out){
		assert(in.exists());
		assert(in.isFile());
		
		if(out.exists()){
			System.out.println("Skipping existing file "+out.getAbsolutePath());
			return;
		}
		
		String abs=in.getAbsolutePath();
		String fname=in.getName();
		
		boolean valid=false;
		
		for(String s : badNames){
			if(fname.matches(s)){
				valid=false;
				return;
			}
		}
		
		for(String s : dirNames){
			if(abs.contains(s)){
				valid=true;
				break;
			}
		}
		
		for(String s : fileNames){
			if(valid){break;}
			if(fname.matches(s)){
				valid=true;
			}
		}
		
		if(!valid){return;}
		
		if(abs.endsWith(".tsv")/* && in.length()>4000000*/){
			out=new File(out.getAbsolutePath()+".zip");
		}
		
//		if(abs.endsWith(".bz2")){
//			out=new File(out.getAbsolutePath().replace(".bz2", ".zip"));
//		}
		
		System.out.println("Copying file to "+out.getAbsolutePath());
		ReadWrite.copyFile(in.getAbsolutePath(), out.getAbsolutePath());
		
	}
	
//	public static String[] inRoots={"F:\\UTSW_batch_1\\", "F:\\UTSW_batch_2\\"};
	/** Input root directories to process for file copying */
	public static String[] inRoots={"F:\\UTSW_second_set\\"};
	/** Output root directory where files will be copied */
	public static String outRoot="C:\\Data\\OCT_8\\";
	
	/** Directory name patterns that validate files for copying */
	public static final String[] dirNames={"\\CNV\\", "\\SV\\"};
	
	/** Absolute path patterns for files to include in copying process */
	public static final String[] fileNamesAbsolute={
		".*\\\\gene-GS.+-ASM.*\\.tsv.*",
		".*\\\\geneVarSummary-GS.+-ASM.*\\.tsv.*",
		".*\\\\summary-GS.+-ASM.*\\.tsv.*",
		".*\\\\var-GS.+-ASM.*\\.tsv.*",
		".*\\\\manifest\\.all",
		".*\\\\README\\..*",
		".*\\\\version",
	};
	
	/** Filename patterns for files to include in copying process */
	public static final String[] fileNames={
		"gene-GS.+-ASM.*\\.tsv.*",
		"geneVarSummary-GS.+-ASM.*\\.tsv.*",
		"summary-GS.+-ASM.*\\.tsv.*",
		"var-GS.+-ASM.*\\.tsv.*",
		"manifest\\.all",
		"README\\..*",
		"version",
	};
	
	/**
	 * Filename patterns to exclude from copying (system files and evidence directories)
	 */
	public static final String[] badNames={
		".*AppleDouble.*",
		".*DS_Store.*",
		".*EVIDENCE.*"
	};
	
		
}
