package fileIO;

import java.io.File;


public class CopyFiles {
	
	
	public static void main(String[] args){
		for(String s : args){
			renameFiles(s);
		}
	}
	
	
	public static void renameFiles(String path){
		File f=new File(path);
		renameFiles(f);
	}
	
	public static void renameFiles(File path){
		
		if(path.isDirectory()){
			File[] array=path.listFiles();
			for(File f : array){renameFiles(f);}
		}else{
			rename(path);
		}
		
	}
	
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
