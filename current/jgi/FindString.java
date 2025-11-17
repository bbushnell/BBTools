package jgi;

import fileIO.TextFile;

/**
 * @author Brian Bushnell
 * @date Jun 18, 2013
 *
 */
public class FindString {
	
	/**
	 * Program entry point that performs text search operation.
	 * Takes a filename as the first argument followed by one or more search strings.
	 * Reads the file line by line and prints any lines containing at least one search string.
	 *
	 * @param args Command-line arguments where args[0] is the filename and
	 * args[1...n] are search strings
	 */
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
