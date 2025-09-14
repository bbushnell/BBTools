package fileIO;

import java.io.File;


public class CompressFiles {
	
	
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
	
	
	public static void compressFiles(String path){
		File f=new File(path);
		compressFiles(f);
	}
	
	public static void compressFiles(File path){
		
		if(path.isDirectory()){
			File[] array=path.listFiles();
			for(File f : array){compressFiles(f);}
		}else{
			compress(path);
		}
		
	}
	
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
	
	
	public static boolean zip=true;
	public static boolean gzip=!zip;
	
	public static boolean PRINT_7Z_BATCH=true;
	
}
