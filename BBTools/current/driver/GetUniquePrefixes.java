package driver;

import java.util.HashSet;

import fileIO.TextFile;

public class GetUniquePrefixes {
	
	public static void main(String[] args){
		
		
		String fname=args[0];
		int prefix=Integer.parseInt(args[1]);
		TextFile tf=new TextFile(fname);
		
		HashSet<String> set=new HashSet<String>();
		
		String header=null;
		StringBuilder sequence=new StringBuilder();
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			if(line.startsWith(">")){
				if(sequence.length()>0){
					if(sequence.length()>prefix){sequence.setLength(prefix);}
					String s=sequence.toString();
					if(set.contains(s)){}
					else{
						set.add(s);
						System.out.println(header+"\n"+s);
					}
				}
				sequence.setLength(0);
				header=line;
			}else{
				sequence.append(line);
			}
		}
		
	}
	
}
