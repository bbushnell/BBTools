package jgi;

import fileIO.TextFile;

/**
 * Simple text file search utility that prints lines containing specified strings.
 * Reads a text file line by line and outputs any lines that contain one or more
 * of the search strings provided as command-line arguments.
 *
 * @author Brian Bushnell
 * @date Jun 18, 2013
 */
public class FindString {
	
	public static void main(String[] args){
		String fname=args[0];
		TextFile tf=new TextFile(fname, true);
		for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
			boolean b=false;
			for(int i=1; i<args.length; i++){
				if(line.contains(args[i])){b=true;break;}
			}
			if(b){System.out.println(line);}
		}
		tf.close();
	}
	
}
