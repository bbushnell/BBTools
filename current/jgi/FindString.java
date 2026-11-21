package jgi;

import fileIO.TextFile;

/**
 * Prints each retained text line containing any supplied literal substring once.
 * Matching is case-sensitive. TextFile's logical-line reader skips blank lines
 * and preserves whitespace in retained lines. An empty search term matches every
 * retained line; no search terms consume the input without printing anything.
 * Reader closure is attempted on normal and exceptional exits; stdout remains borrowed.
 * @author Brian Bushnell, Shinobu
 * @date Jun 18, 2013
 */
public class FindString{

	/**
	 * Searches the input named by the first argument.
	 * @param args Input filename followed by zero or more literal search terms
	 * @throws IllegalArgumentException If no input filename is supplied
	 * @throws RuntimeException If input, output or reader closure fails
	 */
	public static void main(final String[] args){
		if(args==null || args.length<1){
			throw new IllegalArgumentException("Usage: FindString <file> [substring ...]");
		}
		//JGI-155: cleanup must also run when reading or matching fails.
		try(final Input input=new Input(args[0])){
			final TextFile tf=input.reader;
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				boolean found=false;
				for(int i=1; i<args.length; i++){
					if(line.contains(args[i])){found=true; break;}
				}
				if(found){System.out.println(line);}
			}
			//JGI-157: check the borrowed stream without closing it.
			if(System.out.checkError()){throw new RuntimeException("Error writing FindString output");}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns one reader; try-with-resources preserves primary and suppressed failures. */
	private static final class Input implements AutoCloseable{
		Input(final String name_){
			name=name_;
			reader=new TextFile(name, true);
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
