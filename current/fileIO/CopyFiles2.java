package fileIO;

import java.io.File;
import java.nio.file.Path;

import shared.Timer;

/**
 * Copies selected files using a legacy Windows directory profile.
 * Traversal removes exact uppercase ASM components from destinations, excludes
 * matching input paths, and selects files by directory substring or filename regex.
 * Every selected lowercase .tsv input gets a .zip destination suffix, regardless
 * of size; existing base or final ZIP destinations are skipped.
 * Child paths and directory filters use backslashes. Configuration is public
 * mutable state intended for sequential use, not concurrent reconfiguration.
 * @author Brian Bushnell, Shinobu
 */
public class CopyFiles2{

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Copies each configured input root and prints elapsed time to stdout.
	 * A supplied pair replaces the static roots for subsequent calls as well.
	 * @param args No arguments to use the configured profile, or two positional
	 * paths: input and destination. The two-argument count is assertion-checked.
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
	 * Wraps paths as files and delegates to {@link #copyFiles(File, File)}.
	 * @param in Non-null input directory or file path
	 * @param out Non-null corresponding destination path
	 */
	public static void copyFiles(String in, String out){
		File fin=new File(in);
		File fout=new File(out);
		copyFiles(fin, fout);
	}

	/**
	 * Traverses an input tree, flattening exact ASM components in destination paths.
	 * A full absolute-input-path match in {@link #badNames} skips that subtree.
	 * Destination roots and other name components are retained. Missing directories
	 * get a single mkdir attempt; missing ancestors are not recursively created.
	 * A failed directory listing throws, and eligible files go to {@link #copyFile(File, File)}.
	 * Directory links and destinations inside the input tree are not guarded here.
	 *
	 * @param in Non-null input directory or file
	 * @param out Non-null destination before ASM removal
	 */
	public static void copyFiles(File in, File out){
		String abs=in.getAbsolutePath();
		for(String s : badNames){
			if(abs.matches(s)){return;}
		}
		{
			final Path path=out.toPath().toAbsolutePath();
			Path cleaned=path.getRoot();
			//Keep the drive or UNC root and every name other than an exact ASM component.
			for(Path part : path){
				if(!part.toString().equals("ASM")){cleaned=cleaned.resolve(part);}
			}
			out=cleaned.toFile();
		}
		//TODO: Probable bug [FIO-115] - traversal can revisit directory links or a destination created inside the input tree; no cycle/containment guard.
		if(in.isDirectory()){
			if(!out.exists()){out.mkdir();}
			File[] array=in.listFiles();
			if(array==null){throw new RuntimeException("Failed to list directory contents (I/O error): "+in.getAbsolutePath());}
			for(File f : array){
				String outname=out.getAbsolutePath()+"\\"+f.getName();
				File f2=new File(outname);
				copyFiles(f, f2);
			}
		}else{
			copyFile(in, out);
		}
	}

	/**
	 * Copies a selected file without applying the traversal's ASM path rewrite.
	 * An existing destination is skipped before selection. A full filename match
	 * in {@link #badNames} rejects the file; otherwise any directory substring in
	 * {@link #dirNames} or full filename regex in {@link #fileNames} selects it.
	 * Matching uses case-sensitive strings/regexes. Lowercase .tsv input paths append
	 * .zip to the destination, which is checked again for existence. There is no size
	 * threshold. {@link ReadWrite#copyFile(String, String)} handles the actual I/O
	 * and compression without creating parent directories; failures propagate.
	 * Copy and skip messages go to stdout.
	 *
	 * @param in Existing regular input file, assertion-checked
	 * @param out Non-null destination before an optional .zip suffix
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
			if(fname.matches(s)){valid=true;}
		}
		if(!valid){return;}
		//The historical four-megabyte threshold is inactive; all selected TSV inputs qualify.
		if(abs.endsWith(".tsv")){
			out=new File(out.getAbsolutePath()+".zip");
			//Keep the early base-path skip above, and also check the actual ZIP destination.
			if(out.exists()){
				System.out.println("Skipping existing file "+out.getAbsolutePath());
				return;
			}
		}
		//A historical .bz2-to-.zip naming alternative remains inactive.
		System.out.println("Copying file to "+out.getAbsolutePath());
		ReadWrite.copyFile(in.getAbsolutePath(), out.getAbsolutePath());
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	//Earlier input profiles used F:\UTSW_batch_1 and F:\UTSW_batch_2.
	/** Mutable Windows input profile, replaced by main when paths are supplied. */
	public static String[] inRoots={"F:\\UTSW_second_set\\"};
	/** Mutable destination shared by all configured input roots. */
	public static String outRoot="C:\\Data\\OCT_8\\";

	/** Case-sensitive absolute-path substrings selecting files, not regexes. */
	public static final String[] dirNames={"\\CNV\\", "\\SV\\"};

	/** Historical absolute-path regex scaffold; the copying paths do not consult it. */
	public static final String[] fileNamesAbsolute={
		".*\\\\gene-GS.+-ASM.*\\.tsv.*",
		".*\\\\geneVarSummary-GS.+-ASM.*\\.tsv.*",
		".*\\\\summary-GS.+-ASM.*\\.tsv.*",
		".*\\\\var-GS.+-ASM.*\\.tsv.*",
		".*\\\\manifest\\.all",
		".*\\\\README\\..*",
		".*\\\\version",
	};

	/** Case-sensitive full-filename regexes selecting files outside selected directories. */
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
	 * Full-match exclusion regexes: absolute paths during traversal, filenames in copyFile.
	 * Their wildcard patterns exclude system files and EVIDENCE paths case-sensitively.
	 */
	public static final String[] badNames={
		".*AppleDouble.*",
		".*DS_Store.*",
		".*EVIDENCE.*"
	};
}
