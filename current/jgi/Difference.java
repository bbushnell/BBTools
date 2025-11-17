package jgi;

import fileIO.TextFile;

/**
 * @author Brian Bushnell
 * @date Oct 9, 2013
 *
 */
public class Difference {
	
	/**
	 * Compares two text files specified as command-line arguments.
	 * Reads files line by line, reports differences to stderr, and exits
	 * if more than 4 differences are found or files have different lengths.
	 *
	 * @param args Command-line arguments: args[0] = first file, args[1] = second file
	 * @throws AssertionError If more than 4 differences found or files have different lengths
	 */
	public static void main(String[] args){

		TextFile tf1=new TextFile(args[0], false);
		TextFile tf2=new TextFile(args[1], false);

		String s1=tf1.readLine(false);
		String s2=tf2.readLine(false);
		
		int difs=0;
		int i=1;
		while(s1!=null && s2!=null){
			if(!s1.equals(s2)){
				difs++;
				System.err.println("Line "+i+":\n"+s1+"\n"+s2+"\n");
				assert(difs<5);
			}
			i++;
			s1=tf1.readLine(false);
			s2=tf2.readLine(false);
		}
		
		assert(s1==null && s2==null) : "Line "+i+":\n"+s1+"\n"+s2+"\n";
		
		tf1.close();
		tf2.close();
	}
	
}
