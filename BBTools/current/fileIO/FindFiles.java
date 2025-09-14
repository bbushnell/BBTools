package fileIO;

import java.io.File;
import java.util.ArrayList;


public class FindFiles {
	
	
	public static void main(String[] args){
		
		String root=args[0];
//		if(root.equals(".")){root=null;}
		String prefix=args[1];
		String suffix=(args[2].equals("null") ? null : args[2]);
		String middle=null;
		
		if(args.length>3){
			middle=(args[3].equals("null") ? null : args[3]);
		}

		boolean NEWLINE=true;
		boolean BOTH=true;

		ArrayList<String> results=findFiles(root, prefix, suffix, middle);
		for(String s : results){
			if(NEWLINE){
				System.out.println(s);
			}else{
				System.out.print(s+" ");
			}
		}


		if(BOTH){
			System.out.println();
			NEWLINE=!NEWLINE;
			for(String s : results){
				if(NEWLINE){
					System.out.println(s);
				}else{
					System.out.print(s+" ");
				}
			}
		}
	}
	
	
	public FindFiles(String pre, String suf, String mid){
		assert(!"*".equals(pre)) : "Use # instead of *, which has problems from the command line";
		assert(!"*".equals(suf)) : "Use # instead of *, which has problems from the command line";
		prefix=((pre==null || pre.equals("*") || pre.equals("#")) ? null : pre.toLowerCase());
		suffix=((suf==null || suf.equals("*") || suf.equals("#")) ? null : suf.toLowerCase());
		middle=((mid==null || mid.equals("*") || mid.equals("#")) ? null : mid.toLowerCase());
	}
	
	public static ArrayList<String> findFiles(String root, String prefix, String suffix){
		return findFiles(root, prefix, suffix, null);
	}
	
	public static ArrayList<String> findFiles(String root, String prefix, String suffix, String mid){
		FindFiles ff=new FindFiles(prefix, suffix, mid);
		return ff.findFiles(root);
	}
	
	public ArrayList<String> findFiles(String path){
		findFiles(new File(path==null ? "." : path));
		return results;
	}
	
	public ArrayList<String> findFiles(File path){
		
		if(path.isDirectory()){
			File[] array=path.listFiles();
			if(array==null){System.err.println("null contents for "+path.getAbsolutePath());}
			else{for(File f : array){findFiles(f);}}
		}else{
			consider(path);
		}
		return results;
	}
	
	public void consider(File in){
//		System.out.println("Considering "+in.getAbsolutePath()+" versus '"+prefix+"' '"+suffix+"'");
		if(!in.exists()){return;}
		assert(in.exists()) : in;
		assert(in.isFile());
		String abs=in.getAbsolutePath();
//		System.out.println("Considering "+abs);
		String abs2=abs.toLowerCase();
		int slashLoc=abs2.lastIndexOf(slash);
		if(slashLoc>-1){
			abs2=abs2.substring(slashLoc+1);
		}
//		System.out.println("a");
		if(prefix!=null && !abs2.startsWith(prefix)){return;}
//		System.out.println("b");
		if(suffix!=null && !abs2.endsWith(suffix)){return;}
//		System.out.println("c");
		
		if(middle!=null && !abs2.contains(middle)){return;}
		
		results.add(abs);
	}
	
	
	public ArrayList<String> results=new ArrayList<String>();
	public String prefix;
	public String suffix;
	public String middle;
	public static final char slash=System.getProperty("file.separator").charAt(0);
	
}
