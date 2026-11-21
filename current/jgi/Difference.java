package jgi;

import fileIO.TextFile;

/**
 * Compares two text files by physical line and prints mismatches to stderr.
 * With assertions enabled, prints the fifth mismatch before failing; one to four
 * mismatches return normally if the files have equal line counts. Unequal line
 * counts also fail by assertion. Extra arguments are ignored.
 * Closure is attempted for every opened reader on normal and exceptional exits.
 * @author Brian Bushnell, Shinobu
 * @date Oct 9, 2013
 */
public class Difference{

	/**
	 * Compares the first two input names without skipping blank lines.
	 * @param args Two input filenames; additional arguments are ignored
	 * @throws IllegalArgumentException If fewer than two filenames are supplied
	 * @throws RuntimeException If input, output or reader closure fails
	 */
	public static void main(final String[] args){
		if(args==null || args.length<2){
			throw new IllegalArgumentException("Usage: Difference <file1> <file2>");
		}
		//JGI-155: the first reader also needs cleanup if the second open fails.
		try(final Input input1=new Input(args[0]); final Input input2=new Input(args[1])){
			final TextFile tf1=input1.reader, tf2=input2.reader;
			String s1=tf1.readLine(false), s2=tf2.readLine(false);
			int difs=0, i=1;
			while(s1!=null && s2!=null){
				if(!s1.equals(s2)){
					difs++;
					System.err.println("Line "+i+":\n"+s1+"\n"+s2+"\n");
					assert(difs<5) : "Difference stops after printing five mismatches; line="+i;
				}
				i++;
				s1=tf1.readLine(false);
				s2=tf2.readLine(false);
			}
			assert(s1==null && s2==null) : "Unequal physical line counts. Line "+i+":\n"+s1+"\n"+s2+"\n";
			//JGI-157: PrintStream records write failures instead of throwing them.
			if(System.err.checkError()){throw new RuntimeException("Error writing Difference output");}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns one reader; try-with-resources preserves primary and suppressed failures. */
	private static final class Input implements AutoCloseable{
		Input(final String name_){
			name=name_;
			reader=new TextFile(name, false);
		}

		/** Turns TextFile's reported close error into an unsuccessful main invocation. */
		@Override
		public void close(){
			if(reader.close()){throw new RuntimeException("Error closing "+name);}
		}

		final TextFile reader;
		private final String name;
	}
}
