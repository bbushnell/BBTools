package fileIO;

import java.io.File;

import shared.Timer;


public class CopyFiles2 {
	
	
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
	
	
	public static void copyFiles(String in, String out){
		File fin=new File(in);
		File fout=new File(out);
		copyFiles(fin, fout);
	}
	
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
	public static String[] inRoots={"F:\\UTSW_second_set\\"};
	public static String outRoot="C:\\Data\\OCT_8\\";
	
	public static final String[] dirNames={"\\CNV\\", "\\SV\\"};
	
	public static final String[] fileNamesAbsolute={
		".*\\\\gene-GS.+-ASM.*\\.tsv.*",
		".*\\\\geneVarSummary-GS.+-ASM.*\\.tsv.*",
		".*\\\\summary-GS.+-ASM.*\\.tsv.*",
		".*\\\\var-GS.+-ASM.*\\.tsv.*",
		".*\\\\manifest\\.all",
		".*\\\\README\\..*",
		".*\\\\version",
	};
	
	public static final String[] fileNames={
		"gene-GS.+-ASM.*\\.tsv.*",
		"geneVarSummary-GS.+-ASM.*\\.tsv.*",
		"summary-GS.+-ASM.*\\.tsv.*",
		"var-GS.+-ASM.*\\.tsv.*",
		"manifest\\.all",
		"README\\..*",
		"version",
	};
	
	public static final String[] badNames={
		".*AppleDouble.*",
		".*DS_Store.*",
		".*EVIDENCE.*"
	};
	
		
}
