package driver;

import fileIO.TextFile;

public class Grep {
	
	public static void main(String[] args){
		
		TextFile tf=new TextFile(args[0], true);
		
		String s=null;
		
		for(s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(s.contains(args[1])){System.out.println(s);}
		}
		tf.close();
		
	}
	
}
