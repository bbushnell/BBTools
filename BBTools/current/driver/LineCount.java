package driver;

import fileIO.TextFile;

public class LineCount {
	
	public static void main(String[] args){
		
		TextFile tf=new TextFile(args[0], false);
		long lines=tf.countLines();
		tf.close();
		System.out.println(args[0]+" has "+lines+" non-blank lines.");
		
	}
	
}
